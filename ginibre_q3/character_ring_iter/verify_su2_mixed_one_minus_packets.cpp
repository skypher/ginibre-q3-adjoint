#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;
using Pair = std::pair<int, int>;
using State = std::map<Pair, cpp_int>;

static std::vector<int> fuse(int x, int p, int k) {
  std::vector<int> out;
  const int lo = std::abs(x - p);
  const int hi = std::min(x + p, 2 * k - x - p);
  for (int y = lo; y <= hi; y += 2) out.push_back(y);
  return out;
}

static void add_wedge(State& out, int x, int y, const cpp_int& c) {
  if (x == y || c == 0) return;
  if (x > y) out[{x, y}] += c;
  else out[{y, x}] -= c;
}

static State apply_wedge(const State& in, int p, int k) {
  State out;
  for (const auto& [uv, c] : in) {
    const auto [u, v] = uv;
    for (int x : fuse(u, p, k)) add_wedge(out, x, v, c);
    for (int y : fuse(v, p, k)) add_wedge(out, u, y, c);
  }
  for (auto it = out.begin(); it != out.end();) {
    if (it->second == 0) it = out.erase(it);
    else ++it;
  }
  return out;
}

static State apply_symmetric(const State& in, int p, int k) {
  State out;
  for (const auto& [uv, c] : in) {
    const auto [u, v] = uv;
    for (int x : fuse(u, p, k)) out[{x, v}] += c;
    for (int y : fuse(v, p, k)) out[{u, y}] += c;
  }
  return out;
}

static cpp_int coefficient(const State& x, int u, int v) {
  const auto it = x.find({u, v});
  return it == x.end() ? cpp_int(0) : it->second;
}

static cpp_int h_at(const State& x, int u, int v) {
  if (u <= v) return 0;
  const cpp_int c = coefficient(x, u, v);
  return (u & 1) ? c : -c;
}

static int fusion_entry(int target, int source, int p, int k) {
  const auto outputs = fuse(source, p, k);
  return std::find(outputs.begin(), outputs.end(), target) != outputs.end();
}

static int n1_squared_entry(int target, int source, int k) {
  int value = 0;
  for (int middle : fuse(source, 1, k))
    value += fusion_entry(target, middle, 1, k);
  return value;
}

[[noreturn]] static void fail(const std::string& message) {
  std::cerr << "SU2_MIXED_ONE_MINUS_PACKETS FAIL " << message << "\n";
  std::exit(1);
}

struct Counters {
  std::uint64_t words = 0;
  std::uint64_t formula_entries = 0;
  std::uint64_t empirical_entries = 0;
  std::uint64_t q1_identities = 0;
  std::uint64_t dominant_words = 0;
  std::uint64_t simple_current_cases = 0;
};

struct Diagnostic {
  bool found = false;
  int k = 0;
  int q = 0;
  int a = 0;
  cpp_int value = 0;
  std::vector<int> word;
};

static Diagnostic remainder_obstruction;
static Diagnostic stable_remainder_obstruction;
static Diagnostic append_remainder_obstruction;

static void record_obstruction(Diagnostic& obstruction, int k, int q, int a,
                               const cpp_int& value,
                               const std::vector<int>& word) {
  if (value >= 0) return;
  const auto candidate =
      std::tuple(word.size(), k, q, word, a);
  const auto current =
      std::tuple(obstruction.word.size(), obstruction.k,
                 obstruction.q, obstruction.word, obstruction.a);
  if (!obstruction.found || candidate < current) {
    obstruction = {true, k, q, a, value, word};
  }
}

static void inspect_word(const State& x, const State& even_product, int k,
                         int q, int suffix_sum,
                         const std::vector<int>& word, Counters& counters) {
  ++counters.words;
  for (const auto& [uv, unused] : x) {
    (void)unused;
    if (h_at(x, uv.first, uv.second) < 0)
      fail("opposite-parity B cone");
  }

  const State y = apply_wedge(apply_wedge(x, 1, k), 1, k);
  for (int a = 1; a <= k; a += 2) {
    cpp_int packet = 0;
    for (int u = 1; u <= k; u += 2) {
      packet += h_at(x, u, 0)
                * (n1_squared_entry(a, u, k) + (a == u));
    }
    for (int u = 2; u <= k; u += 2)
      packet -= 2 * h_at(x, u, 1) * fusion_entry(a, u, 1, k);
    packet += h_at(x, a, 2);
    if (a == 1) packet += -2 * h_at(x, 1, 0) + h_at(x, 2, 1);
    if (packet != coefficient(y, a, 0))
      fail("corrected P24.8 mismatch");
    ++counters.formula_entries;
    if (coefficient(y, a, 0) < 0)
      fail("empirical mixed packet");
    ++counters.empirical_entries;
  }

  if (q == 1) {
    std::vector<cpp_int> defect(static_cast<std::size_t>(k + 1));
    for (int a = 0; a <= k; a += 2) {
      cpp_int value = -coefficient(even_product, a, 2);
      for (int u = 0; u <= k; u += 2)
        value += fusion_entry(a, u, 2, k)
                 * coefficient(even_product, u, 0);
      defect[static_cast<std::size_t>(a)] = value;
    }
    for (int a = 1; a <= k; a += 2) {
      cpp_int rhs = 0;
      for (int u : fuse(a, 1, k))
        rhs += defect[static_cast<std::size_t>(u)];
      if (rhs != coefficient(y, a, 0))
        fail("q=1 boundary-defect identity");
      ++counters.q1_identities;
    }
  }

  if (q >= 3) {
    for (int a = 1; a <= k; a += 2) {
      cpp_int remainder = 0;
      for (int u = 0; u <= k; u += 2) {
        const cpp_int boundary =
            coefficient(even_product, u, 0)
            + coefficient(even_product, u, 2);
        remainder += fusion_entry(a, u, q, k) * boundary;
        const cpp_int internal =
            coefficient(even_product, u, q - 1)
            + coefficient(even_product, u, q + 1);
        remainder -= fusion_entry(a, u, 1, k) * internal;
      }
      record_obstruction(remainder_obstruction, k, q, a, remainder, word);
      if (q + suffix_sum <= k)
        record_obstruction(stable_remainder_obstruction, k, q, a,
                           remainder, word);
    }
  }

  if (q > suffix_sum && q + suffix_sum <= k) {
    for (int u = 2; u <= k; u += 2)
      if (h_at(x, u, 1) != 0)
        fail("dominant chamber reached negative row");
    ++counters.dominant_words;
  }
  (void)word;
}

static void recurse(State x, State even_product, int k, int q, int max_len,
                    int next_label, int suffix_sum, std::vector<int>& word,
                    Counters& counters) {
  inspect_word(x, even_product, k, q, suffix_sum, word, counters);
  if (static_cast<int>(word.size()) == max_len) return;
  for (int p = next_label; p <= k; p += 2) {
    word.push_back(p);
    const State child_x = apply_wedge(x, p, k);
    const State old_terminal =
        apply_wedge(apply_wedge(x, 1, k), 1, k);
    const State child_terminal =
        apply_wedge(apply_wedge(child_x, 1, k), 1, k);
    for (int a = 1; a <= k; a += 2) {
      cpp_int append_remainder = coefficient(child_terminal, a, 0);
      for (int u = 1; u <= k; u += 2)
        append_remainder -= fusion_entry(a, u, p, k)
                            * coefficient(old_terminal, u, 0);
      record_obstruction(append_remainder_obstruction, k, q, a,
                         append_remainder, word);
    }
    recurse(child_x, apply_symmetric(even_product, p, k), k, q, max_len,
            p, suffix_sum + p, word, counters);
    word.pop_back();
  }
}

static void verify_simple_current(int max_k, int max_power,
                                  Counters& counters) {
  for (int k = 2; k <= max_k; k += 2) {
    for (int q = 1; q <= k; q += 2) {
      State x;
      x[{q, 0}] = 1;
      for (int n = 0; n <= max_power; ++n) {
        const State y = apply_wedge(apply_wedge(x, 1, k), 1, k);
        for (int a = 1; a <= k; a += 2)
          if (coefficient(y, a, 0) < 0)
            fail("simple-current ray");
        ++counters.simple_current_cases;
        x = apply_wedge(x, k, k);
      }
    }
  }
}

int main(int argc, char** argv) {
  const int max_k = argc > 1 ? std::stoi(argv[1]) : 14;
  const int max_len = argc > 2 ? std::stoi(argv[2]) : 7;
  const int max_power = argc > 3 ? std::stoi(argv[3]) : 40;
  Counters counters;
  for (int k = 1; k <= max_k; ++k) {
    for (int q = 1; q <= k; q += 2) {
      State x;
      x[{q, 0}] = 1;
      State even_product;
      even_product[{0, 0}] = 1;
      std::vector<int> word;
      recurse(x, even_product, k, q, max_len, 2, 0, word, counters);
    }
  }
  verify_simple_current(std::max(max_k, 100), max_power, counters);
  std::cout << "SU2_MIXED_ONE_MINUS_PACKETS PASS"
            << " max_k=" << max_k
            << " max_len=" << max_len
            << " max_simple_current_k=" << std::max(max_k, 100)
            << " max_simple_current_power=" << max_power
            << " words=" << counters.words
            << " formula_entries=" << counters.formula_entries
            << " empirical_entries=" << counters.empirical_entries
            << " q1_identities=" << counters.q1_identities
            << " dominant_words=" << counters.dominant_words
            << " simple_current_cases=" << counters.simple_current_cases
            << "\n";
  if (remainder_obstruction.found) {
    std::cout << "NATURAL_REMAINDER_OBSTRUCTION"
              << " k=" << remainder_obstruction.k
              << " q=" << remainder_obstruction.q
              << " a=" << remainder_obstruction.a
              << " value=" << remainder_obstruction.value
              << " word=[";
    for (std::size_t i = 0; i < remainder_obstruction.word.size(); ++i) {
      if (i) std::cout << ",";
      std::cout << remainder_obstruction.word[i];
    }
    std::cout << "]\n";
  } else {
    std::cout << "NATURAL_REMAINDER_OBSTRUCTION none\n";
  }
  if (stable_remainder_obstruction.found) {
    std::cout << "STABLE_REMAINDER_OBSTRUCTION"
              << " k=" << stable_remainder_obstruction.k
              << " q=" << stable_remainder_obstruction.q
              << " a=" << stable_remainder_obstruction.a
              << " value=" << stable_remainder_obstruction.value
              << " word=[";
    for (std::size_t i = 0; i < stable_remainder_obstruction.word.size();
         ++i) {
      if (i) std::cout << ",";
      std::cout << stable_remainder_obstruction.word[i];
    }
    std::cout << "]\n";
  } else {
    std::cout << "STABLE_REMAINDER_OBSTRUCTION none\n";
  }
  if (append_remainder_obstruction.found) {
    std::cout << "APPEND_REMAINDER_OBSTRUCTION"
              << " k=" << append_remainder_obstruction.k
              << " q=" << append_remainder_obstruction.q
              << " a=" << append_remainder_obstruction.a
              << " value=" << append_remainder_obstruction.value
              << " word=[";
    for (std::size_t i = 0; i < append_remainder_obstruction.word.size();
         ++i) {
      if (i) std::cout << ",";
      std::cout << append_remainder_obstruction.word[i];
    }
    std::cout << "]\n";
  } else {
    std::cout << "APPEND_REMAINDER_OBSTRUCTION none\n";
  }
}
