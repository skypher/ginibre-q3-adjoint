// Diagnostic only: exact search for the fused-leaf oscillation estimate
//
//   |C_F(p,q;a)| <= sum_{v in p star a} B_F(q,v)
//
// on actual plus-word backgrounds F.  All arithmetic is exact.  The finite
// mode uses SU(2)_k fusion; the ordinary mode uses the unrestricted SU(2)
// tensor rule.  OpenMP concurrency is selected from CPU availability,
// current load, and available RAM rather than from a fixed thread count.
#include <algorithm>
#include <atomic>
#include <boost/multiprecision/cpp_int.hpp>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <unistd.h>

using Big = boost::multiprecision::cpp_int;

struct Fusion {
  int level;  // -1 means ordinary SU(2).
  int endpoint_limit;
  std::vector<std::vector<std::vector<int>>> channels;

  Fusion(int k, int limit) : level(k), endpoint_limit(limit) {
    const auto size = static_cast<std::size_t>(limit) + 1U;
    channels.resize(size, std::vector<std::vector<int>>(size));
    for (int a = 0; a <= limit; ++a) {
      for (int b = 0; b <= limit; ++b) {
        const int hi =
            level < 0 ? a + b : std::min(a + b, 2 * level - a - b);
        for (int c = std::abs(a - b); c <= hi; c += 2)
          channels[static_cast<std::size_t>(a)]
                  [static_cast<std::size_t>(b)]
                      .push_back(c);
      }
    }
  }

  const std::vector<int>& fuse(int a, int b) const {
    return channels[static_cast<std::size_t>(a)]
                   [static_cast<std::size_t>(b)];
  }
};

struct Table {
  int dim;
  std::vector<Big> x;

  explicit Table(int n)
      : dim(n),
        x(static_cast<std::size_t>(n) * static_cast<std::size_t>(n)) {}

  Big& at(int r, int s) {
    return x[static_cast<std::size_t>(r) *
                 static_cast<std::size_t>(dim) +
             static_cast<std::size_t>(s)];
  }
  const Big& at(int r, int s) const {
    return x[static_cast<std::size_t>(r) *
                 static_cast<std::size_t>(dim) +
             static_cast<std::size_t>(s)];
  }
  Big get(int r, int s) const {
    if (r < 0 || s < 0 || r >= dim || s >= dim) return 0;
    return at(r, s);
  }
};

struct Evaluation {
  Big c_first;
  Big c_second;
  Big c;
  std::vector<std::pair<int, Big>> rhs_terms;
  Big rhs;
  Big lhs;
  Big gap;
};

struct Witness {
  bool present = false;
  std::size_t word_index = std::numeric_limits<std::size_t>::max();
  std::vector<int> word;
  int p = -1;
  int q = -1;
  int a = -1;
  Evaluation eval;
  Table table{1};
};

struct RatioWitness {
  bool present = false;
  Big num;
  Big den;
  std::vector<int> word;
  int p = -1;
  int q = -1;
  int a = -1;
  Big gap;
};

struct BatchStats {
  std::uint64_t words = 0;
  std::uint64_t triples = 0;
  std::uint64_t equalities = 0;
  std::uint64_t nonzero_equalities = 0;
  Witness failure;
  Witness sign_failure;
  Witness negative_boundary;
  RatioWitness tight_ratio;
  RatioWitness tight_disjoint_ratio;
  RatioWitness disjoint_equality;
};

static Big abs_big(const Big& z) { return z < 0 ? -z : z; }

static Table table_for(const std::vector<int>& word, const Fusion& fusion) {
  const int dim =
      fusion.level < 0
          ? std::accumulate(word.begin(), word.end(), 0) + 1
          : fusion.level + 1;
  Table f(dim);
  f.at(0, 0) = 1;
  for (int label : word) {
    Table g(dim);
    for (int r = 0; r < dim; ++r) {
      for (int s = 0; s < dim; ++s) {
        if (f.at(r, s) == 0) continue;
        for (int u : fusion.fuse(r, label))
          if (u < dim) g.at(u, s) += f.at(r, s);
        for (int v : fusion.fuse(s, label))
          if (v < dim) g.at(r, v) += f.at(r, s);
      }
    }
    f = std::move(g);
  }
  return f;
}

static Big boundary(const Table& f, const Fusion& fusion, int q, int a) {
  Big ans = 0;
  for (int t : fusion.fuse(q, a)) ans += f.get(0, t);
  ans -= f.get(q, a);
  return ans;
}

static Big updated_entry(const Table& f, const Fusion& fusion, int leaf,
                         int r, int s) {
  Big ans = 0;
  for (int u : fusion.fuse(leaf, r)) ans += f.get(u, s);
  for (int v : fusion.fuse(leaf, s)) ans += f.get(r, v);
  return ans;
}

static Big updated_boundary(const Table& f, const Fusion& fusion, int leaf,
                            int q, int a) {
  Big ans = 0;
  for (int t : fusion.fuse(q, a))
    ans += updated_entry(f, fusion, leaf, 0, t);
  ans -= updated_entry(f, fusion, leaf, q, a);
  return ans;
}

static Evaluation evaluate(const Table& f, const Fusion& fusion, int p,
                           int q, int a) {
  Evaluation e;
  for (int t : fusion.fuse(q, a)) e.c_first += f.get(p, t);
  for (int u : fusion.fuse(p, q)) e.c_second += f.get(u, a);
  e.c = e.c_first - e.c_second;
  e.lhs = abs_big(e.c);
  for (int v : fusion.fuse(p, a)) {
    Big b = boundary(f, fusion, q, v);
    e.rhs_terms.emplace_back(v, b);
    e.rhs += b;
  }
  e.gap = e.rhs - e.lhs;
  return e;
}

static void enumerate_words_rec(int pos, int n, int max_label, int lo,
                                std::vector<int>& word,
                                std::vector<std::vector<int>>& out) {
  if (pos == n) {
    out.push_back(word);
    return;
  }
  for (int p = lo; p <= max_label; ++p) {
    word[static_cast<std::size_t>(pos)] = p;
    enumerate_words_rec(pos + 1, n, max_label, p, word, out);
  }
}

static std::vector<std::vector<int>> enumerate_words(int n, int max_label) {
  if (n == 0) return {{}};
  std::vector<std::vector<int>> out;
  std::vector<int> word(static_cast<std::size_t>(n));
  enumerate_words_rec(0, n, max_label, 1, word, out);
  return out;
}

static double load_average() {
  std::ifstream in("/proc/loadavg");
  double load = 0.0;
  in >> load;
  return load;
}

static std::uint64_t available_ram_bytes() {
  std::ifstream in("/proc/meminfo");
  std::string key;
  std::uint64_t kib = 0;
  std::string unit;
  while (in >> key >> kib >> unit) {
    if (key == "MemAvailable:") return kib * 1024;
  }
  return 0;
}

static int choose_threads(int dim) {
#ifdef _OPENMP
  const int cpus = std::max(1, omp_get_num_procs());
#else
  const int cpus =
      std::max(1u, std::thread::hardware_concurrency());
#endif
  const double load = load_average();
  // Leave CPUs already occupied by the one-minute load and one additional
  // scheduling margin.  If load exceeds capacity, retain one worker.
  const int by_load =
      std::max(1, cpus - static_cast<int>(std::ceil(load)) - 1);
  const std::uint64_t ram = available_ram_bytes();
  // cpp_int storage is dynamic.  This deliberately overestimates each
  // worker as six dense tables with 128 bytes per entry and uses half of
  // currently available RAM.
  const auto dim_u = static_cast<std::uint64_t>(dim);
  const std::uint64_t bytes_per_worker =
      std::max<std::uint64_t>(1, 6ULL * dim_u * dim_u * 128ULL);
  const auto cpus_u = static_cast<std::uint64_t>(cpus);
  const int by_ram =
      ram == 0
          ? cpus
          : static_cast<int>(std::max<std::uint64_t>(
                1, std::min(cpus_u, ram / 2 / bytes_per_worker)));
  return std::max(1, std::min({cpus, by_load, by_ram}));
}

static bool witness_before(const Witness& candidate, const Witness& current) {
  if (!current.present) return true;
  return std::tie(candidate.word_index, candidate.p, candidate.q, candidate.a) <
         std::tie(current.word_index, current.p, current.q, current.a);
}

static bool ratio_better(const RatioWitness& candidate,
                         const RatioWitness& current) {
  if (!current.present) return true;
  return candidate.num * current.den > current.num * candidate.den;
}

static BatchStats inspect_batch(const std::vector<std::vector<int>>& words,
                                const Fusion& fusion, int label_limit,
                                int threads) {
  BatchStats total;
  total.words = words.size();
  std::mutex merge_mutex;

#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    BatchStats local;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (std::size_t wi = 0; wi < words.size(); ++wi) {
      const auto& word = words[wi];
      Table f = table_for(word, fusion);

      // Check the induction premise independently.
      for (int q = 1; q <= label_limit; ++q) {
        for (int a = 0; a <= label_limit; ++a) {
          const Big b = boundary(f, fusion, q, a);
          if (b < 0) {
            Witness w;
            w.present = true;
            w.word_index = wi;
            w.word = word;
            w.p = 0;
            w.q = q;
            w.a = a;
            w.eval.rhs = b;
            w.table = f;
            if (witness_before(w, local.negative_boundary))
              local.negative_boundary = std::move(w);
          }
        }
      }

      for (int p = 1; p <= label_limit; ++p) {
        for (int q = 1; q <= label_limit; ++q) {
          for (int a = 0; a <= label_limit; ++a) {
            ++local.triples;
            Evaluation e = evaluate(f, fusion, p, q, a);
            const Big sign_test = Big(p - a) * e.c;
            if (sign_test > 0 || (p == a && e.c != 0)) {
              const bool earlier =
                  !local.sign_failure.present ||
                  std::tie(wi, p, q, a) <
                      std::tie(local.sign_failure.word_index,
                               local.sign_failure.p, local.sign_failure.q,
                               local.sign_failure.a);
              if (earlier) {
                Witness w;
                w.present = true;
                w.word_index = wi;
                w.word = word;
                w.p = p;
                w.q = q;
                w.a = a;
                w.eval = e;
                w.table = f;
                local.sign_failure = std::move(w);
              }
            }
            if (e.gap == 0) {
              ++local.equalities;
              if (e.lhs != 0) ++local.nonzero_equalities;
            }
            // Exclude a=0 from the extremizer report: there the estimate is
            // the forced identity C_F(p,q;0)=-B_F(q,p).
            if (a > 0 && e.rhs > 0 && e.lhs > 0 && e.gap > 0) {
              RatioWitness rw;
              rw.present = true;
              rw.num = e.lhs;
              rw.den = e.rhs;
              rw.word = word;
              rw.p = p;
              rw.q = q;
              rw.a = a;
              rw.gap = e.gap;
              if (ratio_better(rw, local.tight_ratio))
                local.tight_ratio = std::move(rw);
            }
            const bool disjoint =
                a > 0 && p != q && p != a &&
                std::find(word.begin(), word.end(), q) == word.end() &&
                std::find(word.begin(), word.end(), a) == word.end();
            if (disjoint && e.rhs > 0 && e.lhs > 0 && e.gap > 0) {
              RatioWitness rw;
              rw.present = true;
              rw.num = e.lhs;
              rw.den = e.rhs;
              rw.word = word;
              rw.p = p;
              rw.q = q;
              rw.a = a;
              rw.gap = e.gap;
              if (ratio_better(rw, local.tight_disjoint_ratio))
                local.tight_disjoint_ratio = std::move(rw);
            }
            if (disjoint && e.rhs > 0 && e.lhs > 0 && e.gap == 0 &&
                !local.disjoint_equality.present) {
              local.disjoint_equality.present = true;
              local.disjoint_equality.num = e.lhs;
              local.disjoint_equality.den = e.rhs;
              local.disjoint_equality.word = word;
              local.disjoint_equality.p = p;
              local.disjoint_equality.q = q;
              local.disjoint_equality.a = a;
              local.disjoint_equality.gap = 0;
            }
            if (e.gap < 0) {
              Witness w;
              w.present = true;
              w.word_index = wi;
              w.word = word;
              w.p = p;
              w.q = q;
              w.a = a;
              w.eval = std::move(e);
              w.table = f;
              if (witness_before(w, local.failure))
                local.failure = std::move(w);
            }
          }
        }
      }
    }

    std::lock_guard<std::mutex> lock(merge_mutex);
    total.triples += local.triples;
    total.equalities += local.equalities;
    total.nonzero_equalities += local.nonzero_equalities;
    if (local.failure.present &&
        witness_before(local.failure, total.failure))
      total.failure = std::move(local.failure);
    if (local.sign_failure.present &&
        witness_before(local.sign_failure, total.sign_failure))
      total.sign_failure = std::move(local.sign_failure);
    if (local.negative_boundary.present &&
        witness_before(local.negative_boundary, total.negative_boundary))
      total.negative_boundary = std::move(local.negative_boundary);
    if (local.tight_ratio.present &&
        ratio_better(local.tight_ratio, total.tight_ratio))
      total.tight_ratio = std::move(local.tight_ratio);
    if (local.tight_disjoint_ratio.present &&
        ratio_better(local.tight_disjoint_ratio,
                     total.tight_disjoint_ratio))
      total.tight_disjoint_ratio = std::move(local.tight_disjoint_ratio);
    if (local.disjoint_equality.present &&
        !total.disjoint_equality.present)
      total.disjoint_equality = std::move(local.disjoint_equality);
  }
  return total;
}

static std::string word_string(const std::vector<int>& word) {
  std::ostringstream out;
  out << '[';
  for (std::size_t i = 0; i < word.size(); ++i) {
    if (i) out << ',';
    out << word[i];
  }
  out << ']';
  return out.str();
}

static void print_table(const Table& f) {
  std::cout << "F_NONZERO";
  for (int r = 0; r < f.dim; ++r)
    for (int s = 0; s < f.dim; ++s)
      if (f.at(r, s) != 0)
        std::cout << " (" << r << ',' << s << ")=" << f.at(r, s);
  std::cout << '\n';
}

static void print_evaluation(const Evaluation& e) {
  std::cout << "C_FIRST=" << e.c_first << " C_SECOND=" << e.c_second
            << " C=" << e.c << " ABS_C=" << e.lhs << '\n';
  std::cout << "RHS_TERMS";
  for (const auto& [v, b] : e.rhs_terms)
    std::cout << " B(q," << v << ")=" << b;
  std::cout << '\n';
  std::cout << "RHS=" << e.rhs << " GAP=" << e.gap << '\n';
}

static int replay(int argc, char** argv) {
  // --replay finite k p q a [word labels...]
  // --replay ordinary endpoint_limit p q a [word labels...]
  if (argc < 7) {
    std::cerr << "replay usage: --replay finite|ordinary K_OR_LIMIT p q a "
                 "[word...]\n";
    return 2;
  }
  const std::string kind = argv[2];
  const int level_or_limit = std::atoi(argv[3]);
  const int p = std::atoi(argv[4]);
  const int q = std::atoi(argv[5]);
  const int a = std::atoi(argv[6]);
  std::vector<int> word;
  for (int i = 7; i < argc; ++i) word.push_back(std::atoi(argv[i]));
  const int level = kind == "finite" ? level_or_limit : -1;
  int limit = level_or_limit;
  if (level < 0) {
    limit = std::max({limit, p, q, a,
                      std::accumulate(word.begin(), word.end(), 0)});
  }
  Fusion fusion(level, limit);
  Table f = table_for(word, fusion);
  Evaluation e = evaluate(f, fusion, p, q, a);
  std::cout << "REPLAY mode=" << kind << " parameter=" << level_or_limit
            << " word=" << word_string(word) << " p=" << p << " q=" << q
            << " a=" << a << '\n';
  print_table(f);
  print_evaluation(e);
  const Big after_p = updated_boundary(f, fusion, p, q, a);
  const Big after_a = updated_boundary(f, fusion, a, q, p);
  std::cout << "IDENTITY B_AFTER_P=" << after_p
            << " RHS_PLUS_C=" << e.rhs + e.c
            << " B_AFTER_A=" << after_a
            << " RHS_MINUS_C=" << e.rhs - e.c << '\n';
  if (after_p != e.rhs + e.c || after_a != e.rhs - e.c) {
    std::cerr << "INTERNAL_IDENTITY_FAILURE\n";
    return 4;
  }
  const bool sign_failure =
      Big(p - a) * e.c > 0 || (p == a && e.c != 0);
  std::cout << "SIGN_TEST (p-a)*C=" << Big(p - a) * e.c
            << " result=" << (sign_failure ? "FAIL" : "PASS") << '\n';
  if (e.gap < 0) return 1;
  return sign_failure ? 5 : 0;
}

int main(int argc, char** argv) {
  if (argc > 1 && std::string(argv[1]) == "--replay")
    return replay(argc, argv);

  // finite [max_k=10] [max_n=7] [min_k=1] [min_n=0]
  // ordinary [max_label=8] [max_n=7] [ignored=1] [min_n=0]
  const std::string mode = argc > 1 ? argv[1] : "finite";
  const int first_limit =
      argc > 2 ? std::atoi(argv[2]) : (mode == "finite" ? 10 : 8);
  const int max_n = argc > 3 ? std::atoi(argv[3]) : 7;
  const int min_k = argc > 4 ? std::atoi(argv[4]) : 1;
  const int min_n = argc > 5 ? std::atoi(argv[5]) : 0;
  if (mode != "finite" && mode != "ordinary") {
    std::cerr << "usage: finite [max_k] [max_n] [min_k] [min_n] | "
                 "ordinary [max_label] [max_n] [ignored] [min_n]\n";
    return 2;
  }

  std::uint64_t grand_words = 0;
  std::uint64_t grand_triples = 0;
  bool sign_counterexample_reported = false;
  const int start = mode == "finite" ? min_k : first_limit;
  const int stop = first_limit;
  for (int parameter = start; parameter <= stop; ++parameter) {
    const int level = mode == "finite" ? parameter : -1;
    const int label_limit = parameter;
    Fusion fusion(level,
                  level < 0
                      ? std::max(max_n * label_limit, label_limit)
                      : label_limit);
    const int dim = level < 0 ? max_n * label_limit + 1 : level + 1;
    const int threads = choose_threads(dim);
    std::cout << "CONFIG mode=" << mode << " parameter=" << parameter
              << " max_n=" << max_n << " threads=" << threads
              << " cpus=" << sysconf(_SC_NPROCESSORS_ONLN)
              << " load1=" << std::fixed << std::setprecision(2)
              << load_average() << " mem_available=" << available_ram_bytes()
              << '\n';
    for (int n = min_n; n <= max_n; ++n) {
      const auto words = enumerate_words(n, label_limit);
      BatchStats stats =
          inspect_batch(words, fusion, label_limit, threads);
      grand_words += stats.words;
      grand_triples += stats.triples;
      std::cout << "BATCH parameter=" << parameter << " n=" << n
                << " words=" << stats.words << " triples=" << stats.triples
                << " equalities=" << stats.equalities
                << " nonzero_equalities=" << stats.nonzero_equalities;
      if (stats.tight_ratio.present) {
        const auto& t = stats.tight_ratio;
        std::cout << " max_strict_ratio=" << t.num << '/' << t.den
                  << " ratio_word=" << word_string(t.word)
                  << " ratio_pqa=" << t.p << ',' << t.q << ',' << t.a
                  << " ratio_gap=" << t.gap;
      }
      if (stats.tight_disjoint_ratio.present) {
        const auto& t = stats.tight_disjoint_ratio;
        std::cout << " max_disjoint_strict_ratio=" << t.num << '/' << t.den
                  << " disjoint_word=" << word_string(t.word)
                  << " disjoint_pqa=" << t.p << ',' << t.q << ',' << t.a
                  << " disjoint_gap=" << t.gap;
      }
      if (stats.disjoint_equality.present) {
        const auto& t = stats.disjoint_equality;
        std::cout << " disjoint_equality=" << t.num << '/' << t.den
                  << " equality_word=" << word_string(t.word)
                  << " equality_pqa=" << t.p << ',' << t.q << ',' << t.a;
      }
      std::cout << '\n';

      if (stats.negative_boundary.present) {
        const auto& w = stats.negative_boundary;
        std::cout << "NEGATIVE_BOUNDARY mode=" << mode
                  << " parameter=" << parameter
                  << " word=" << word_string(w.word) << " q=" << w.q
                  << " a=" << w.a << " B=" << w.eval.rhs << '\n';
        print_table(w.table);
        return 3;
      }
      if (stats.sign_failure.present) {
        if (!sign_counterexample_reported) {
          const auto& w = stats.sign_failure;
          std::cout << "SIGN_COUNTEREXAMPLE mode=" << mode
                    << " parameter=" << parameter
                    << " word=" << word_string(w.word) << " p=" << w.p
                    << " q=" << w.q << " a=" << w.a
                    << " (p-a)*C=" << Big(w.p - w.a) * w.eval.c << '\n';
          print_table(w.table);
          print_evaluation(w.eval);
          std::cout << "REPLAY_COMMAND " << argv[0] << " --replay " << mode
                    << ' ' << parameter << ' ' << w.p << ' ' << w.q << ' '
                    << w.a;
          for (int z : w.word) std::cout << ' ' << z;
          std::cout << '\n';
          sign_counterexample_reported = true;
        }
      }
      if (stats.failure.present) {
        const auto& w = stats.failure;
        std::cout << "COUNTEREXAMPLE mode=" << mode
                  << " parameter=" << parameter
                  << " word=" << word_string(w.word) << " p=" << w.p
                  << " q=" << w.q << " a=" << w.a << '\n';
        print_table(w.table);
        print_evaluation(w.eval);
        const Big after_p =
            updated_boundary(w.table, fusion, w.p, w.q, w.a);
        const Big after_a =
            updated_boundary(w.table, fusion, w.a, w.q, w.p);
        std::cout << "IDENTITY B_AFTER_P=" << after_p
                  << " RHS_PLUS_C=" << w.eval.rhs + w.eval.c
                  << " B_AFTER_A=" << after_a
                  << " RHS_MINUS_C=" << w.eval.rhs - w.eval.c << '\n';
        std::cout << "REPLAY_COMMAND " << argv[0] << " --replay " << mode
                  << ' ' << parameter << ' ' << w.p << ' ' << w.q << ' '
                  << w.a;
        for (int z : w.word) std::cout << ' ' << z;
        std::cout << '\n';
        return 1;
      }
    }
    if (mode == "ordinary") break;
  }
  std::cout << "PASS_OSCILLATION mode=" << mode << " parameters=" << start
            << ".." << stop << " min_n=" << min_n << " max_n=" << max_n
            << " words=" << grand_words << " triples=" << grand_triples
            << " sign_conjecture="
            << (sign_counterexample_reported ? "FAIL" : "PASS_IN_RANGE")
            << '\n';
  return sign_counterexample_reported ? 5 : 0;
}
