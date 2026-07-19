#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Cell = std::pair<int, int>;
using Part = std::vector<int>;

constexpr int kMaxM = 10;

struct Key {
  int part = -1;
  uint8_t len = 0;
  std::array<uint32_t, kMaxM> strips{};

  bool operator==(const Key &other) const {
    if (part != other.part || len != other.len) return false;
    for (int i = 0; i < len; ++i) {
      if (strips[i] != other.strips[i]) return false;
    }
    return true;
  }
};

struct KeyHash {
  std::size_t operator()(const Key &k) const {
    uint64_t h = 1469598103934665603ULL;
    auto mix = [&](uint64_t x) {
      h ^= x;
      h *= 1099511628211ULL;
    };
    mix(static_cast<uint64_t>(k.part));
    mix(static_cast<uint64_t>(k.len));
    for (int i = 0; i < k.len; ++i) mix(k.strips[i]);
    return static_cast<std::size_t>(h);
  }
};

struct Trans {
  int to = -1;
  std::array<Cell, 2> strip{};
  uint32_t code = 0;
};

std::string part_key(const Part &p) {
  std::ostringstream out;
  for (std::size_t i = 0; i < p.size(); ++i) {
    if (i) out << ',';
    out << p[i];
  }
  return out.str();
}

uint32_t cell_code(Cell c) {
  return (static_cast<uint32_t>(c.first) << 6) | static_cast<uint32_t>(c.second);
}

uint32_t strip_code(std::array<Cell, 2> strip) {
  uint32_t a = cell_code(strip[0]);
  uint32_t b = cell_code(strip[1]);
  if (b < a) std::swap(a, b);
  return (a << 12) | b;
}

std::string strip_string(uint32_t code) {
  uint32_t a = code >> 12;
  uint32_t b = code & ((1u << 12) - 1u);
  auto cell = [](uint32_t x) {
    std::ostringstream out;
    out << '(' << (x >> 6) << ',' << (x & 63u) << ')';
    return out.str();
  };
  return cell(a) + ";" + cell(b);
}

int size_of(const Part &p) {
  int n = 0;
  for (int x : p) n += x;
  return n;
}

bool has_cell(const Part &p, Cell c) {
  int r = c.first;
  int col = c.second;
  return r >= 1 && r <= static_cast<int>(p.size()) && col >= 1 && col <= p[r - 1];
}

std::vector<Cell> cells_of(const Part &p) {
  std::vector<Cell> cells;
  cells.reserve(size_of(p));
  for (int r = 1; r <= static_cast<int>(p.size()); ++r) {
    for (int c = 1; c <= p[r - 1]; ++c) cells.push_back({r, c});
  }
  return cells;
}

bool contains_part(const Part &q, const Part &p) {
  const std::size_t n = std::max(q.size(), p.size());
  for (std::size_t i = 0; i < n; ++i) {
    int qi = i < q.size() ? q[i] : 0;
    int pi = i < p.size() ? p[i] : 0;
    if (qi < pi) return false;
  }
  return true;
}

void gen_parts_rec(int n, int max_first, Part &cur, std::vector<Part> &out) {
  if (n == 0) {
    out.push_back(cur);
    return;
  }
  for (int x = std::min(n, max_first); x >= 1; --x) {
    cur.push_back(x);
    gen_parts_rec(n - x, x, cur, out);
    cur.pop_back();
  }
}

std::vector<Part> gen_parts(int n) {
  std::vector<Part> out;
  Part cur;
  gen_parts_rec(n, n, cur, out);
  return out;
}

std::string key_strips(const Key &k) {
  std::ostringstream out;
  for (int i = 0; i < k.len; ++i) {
    if (i) out << '|';
    out << strip_string(k.strips[i]);
  }
  return out.str();
}

double elapsed_seconds(std::chrono::steady_clock::time_point t0) {
  using D = std::chrono::duration<double>;
  return D(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

int main(int argc, char **argv) {
  int max_m = 8;
  if (argc >= 2) max_m = std::atoi(argv[1]);
  if (max_m < 1 || max_m > kMaxM) {
    std::cerr << "max_m must be in [1," << kMaxM << "]\n";
    return 2;
  }

  const int max_size = 4 * max_m;
  std::vector<Part> parts;
  std::vector<std::vector<int>> by_size(max_size + 1);
  std::unordered_map<std::string, int> part_index;

  for (int n = 0; n <= max_size; ++n) {
    for (const Part &p : gen_parts(n)) {
      int idx = static_cast<int>(parts.size());
      parts.push_back(p);
      by_size[n].push_back(idx);
      part_index.emplace(part_key(p), idx);
    }
  }

  std::vector<std::vector<Trans>> horiz(parts.size());
  for (int pi = 0; pi < static_cast<int>(parts.size()); ++pi) {
    const Part &p = parts[pi];
    int n = size_of(p);
    if (n + 2 > max_size) continue;
    for (int qi : by_size[n + 2]) {
      const Part &q = parts[qi];
      if (!contains_part(q, p)) continue;
      std::vector<Cell> diff;
      for (Cell c : cells_of(q)) {
        if (!has_cell(p, c)) diff.push_back(c);
      }
      if (diff.size() != 2) continue;
      if (diff[0].second == diff[1].second) continue;
      Trans tr;
      tr.to = qi;
      tr.strip = {diff[0], diff[1]};
      tr.code = strip_code(tr.strip);
      horiz[pi].push_back(tr);
    }
  }

  int threads = 1;
#ifdef _OPENMP
  threads = omp_get_max_threads();
#endif
  std::cout << "TIGHT_OMEGA_PACKED_CERT max_m=" << max_m
            << " partitions=" << parts.size()
            << " threads=" << threads << std::endl;
  std::cout.flush();

  int empty = part_index.at("");
  for (int m = 1; m <= max_m; ++m) {
    std::unordered_map<Key, uint64_t, KeyHash> states;
    Key root;
    root.part = empty;
    root.len = 0;
    states[root] = 1;

    for (int r = 1; r <= m; ++r) {
      auto level_start = std::chrono::steady_clock::now();
      std::vector<std::pair<Key, uint64_t>> current;
      current.reserve(states.size());
      for (const auto &kv : states) current.push_back(kv);

      std::vector<std::unordered_map<Key, uint64_t, KeyHash>> locals(threads);
      for (auto &mp : locals) {
        mp.reserve((current.size() / std::max(1, threads)) * 4 + 1024);
      }

      std::atomic<std::size_t> done{0};
      const std::size_t report_stride = std::max<std::size_t>(1, current.size() / 20);
      Cell forced{2 * r - 1, 1};

#pragma omp parallel for schedule(dynamic, 4096)
      for (std::int64_t idx = 0; idx < static_cast<std::int64_t>(current.size()); ++idx) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto &local = locals[tid];
        const Key &st = current[static_cast<std::size_t>(idx)].first;
        uint64_t multiplicity = current[static_cast<std::size_t>(idx)].second;
        if (!has_cell(parts[st.part], forced)) {
          for (const Trans &odd : horiz[st.part]) {
            bool has_forced = false;
            Cell other{-1, -1};
            for (Cell c : odd.strip) {
              if (c == forced) {
                has_forced = true;
              } else {
                other = c;
              }
            }
            if (!has_forced || other.second == 1) continue;
            for (const Trans &even : horiz[odd.to]) {
              Key ns = st;
              ns.part = even.to;
              ns.strips[ns.len] = even.code;
              ++ns.len;
              local[ns] += multiplicity;
            }
          }
        }

        std::size_t after = done.fetch_add(1, std::memory_order_relaxed) + 1;
        if (after % report_stride == 0 || after == current.size()) {
#pragma omp critical
          {
            std::cout << "progress m=" << m << " r=" << r
                      << " processed=" << after << "/" << current.size()
                      << " elapsed=" << std::fixed << std::setprecision(2)
                      << elapsed_seconds(level_start) << "s" << std::endl;
            std::cout.flush();
          }
        }
      }

      std::unordered_map<Key, uint64_t, KeyHash> next;
      std::size_t reserve_total = 0;
      for (const auto &mp : locals) reserve_total += mp.size();
      next.reserve(reserve_total);
      std::size_t merged = 0;
      for (auto &mp : locals) {
        for (const auto &kv : mp) next[kv.first] += kv.second;
        merged += mp.size();
        std::cout << "merge m=" << m << " r=" << r
                  << " merged_local_states=" << merged << "/" << reserve_total
                  << " elapsed=" << std::fixed << std::setprecision(2)
                  << elapsed_seconds(level_start) << "s" << std::endl;
        std::cout.flush();
      }

      uint64_t total = 0;
      uint64_t max_count = 0;
      for (const auto &kv : next) {
        total += kv.second;
        max_count = std::max(max_count, kv.second);
      }
      std::cout << "LEVEL m=" << m << " r=" << r
                << " states=" << next.size()
                << " total_paths=" << total
                << " max_fiber=" << max_count
                << " elapsed=" << std::fixed << std::setprecision(2)
                << elapsed_seconds(level_start) << "s" << std::endl;
      std::cout.flush();
      states.swap(next);
    }

    uint64_t total = 0;
    uint64_t max_count = 0;
    Key best;
    for (const auto &kv : states) {
      total += kv.second;
      if (kv.second > max_count) {
        max_count = kv.second;
        best = kv.first;
      }
    }
    uint64_t pow2 = (m < 63 ? (1ULL << m) : 0);
    std::cout << "RESULT m=" << m
              << " states=" << states.size()
              << " total_paths=" << total
              << " max_omega=" << max_count
              << " pow2=" << pow2
              << " passes_pow2=" << (max_count <= pow2 ? "yes" : "no")
              << std::endl;
    std::cout << "BEST m=" << m
              << " partition=" << part_key(parts[best.part])
              << " strips=" << key_strips(best) << std::endl;
    std::cout.flush();
  }

  std::cout << "__EXIT_STATUS=0" << std::endl;
  return 0;
}
