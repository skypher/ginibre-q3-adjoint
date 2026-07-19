#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using Cell = std::pair<int, int>;
using Part = std::vector<int>;

std::string part_key(const Part &p) {
  std::ostringstream out;
  for (std::size_t i = 0; i < p.size(); ++i) {
    if (i) out << ',';
    out << p[i];
  }
  return out.str();
}

std::string strip_key(std::vector<Cell> strip) {
  std::sort(strip.begin(), strip.end());
  std::ostringstream out;
  out << '(';
  for (std::size_t i = 0; i < strip.size(); ++i) {
    if (i) out << ';';
    out << strip[i].first << ',' << strip[i].second;
  }
  out << ')';
  return out.str();
}

std::string append_strip_key(const std::string &prefix, const std::string &s) {
  if (prefix.empty()) return s;
  return prefix + "|" + s;
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

struct Trans {
  int to = -1;
  std::vector<Cell> strip;
  std::string strip_s;
};

struct State {
  int part = -1;
  std::string strips;

  bool operator==(const State &other) const {
    return part == other.part && strips == other.strips;
  }
};

struct StateHash {
  std::size_t operator()(const State &s) const {
    return (static_cast<std::size_t>(s.part) * 1315423911u) ^
           std::hash<std::string>{}(s.strips);
  }
};

std::string count_to_string(unsigned long long x) {
  return std::to_string(x);
}

}  // namespace

int main(int argc, char **argv) {
  int max_m = 7;
  if (argc >= 2) max_m = std::atoi(argv[1]);
  if (max_m < 1) return 2;

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
    auto pcells = cells_of(p);
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
      tr.strip = diff;
      tr.strip_s = strip_key(diff);
      horiz[pi].push_back(tr);
    }
  }

  std::unordered_map<State, unsigned long long, StateHash> states;
  int empty = part_index.at("");
  states[{empty, ""}] = 1;

  std::cout << "TIGHT_OMEGA_SMALL_DIAGNOSTIC max_m=" << max_m
            << " partitions=" << parts.size() << std::endl;

  for (int m = 1; m <= max_m; ++m) {
    // Continue from previous endpoint set when m increments is not valid, so
    // rebuild from scratch for each m.
    states.clear();
    states[{empty, ""}] = 1;

    for (int r = 1; r <= m; ++r) {
      std::unordered_map<State, unsigned long long, StateHash> next;
      Cell forced{2 * r - 1, 1};
      for (const auto &kv : states) {
        const State &st = kv.first;
        unsigned long long multiplicity = kv.second;
        if (has_cell(parts[st.part], forced)) continue;
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
            State ns;
            ns.part = even.to;
            ns.strips = append_strip_key(st.strips, even.strip_s);
            next[ns] += multiplicity;
          }
        }
      }
      unsigned long long total = 0;
      unsigned long long max_count = 0;
      for (const auto &kv : next) {
        total += kv.second;
        max_count = std::max(max_count, kv.second);
      }
      std::cout << "progress m=" << m << " r=" << r
                << " states=" << next.size()
                << " total_paths=" << count_to_string(total)
                << " max_fiber=" << count_to_string(max_count) << std::endl;
      states.swap(next);
    }

    unsigned long long total = 0;
    unsigned long long max_count = 0;
    State best;
    for (const auto &kv : states) {
      total += kv.second;
      if (kv.second > max_count) {
        max_count = kv.second;
        best = kv.first;
      }
    }
    unsigned long long pow2 = 1ULL << m;
    std::cout << "RESULT m=" << m
              << " states=" << states.size()
              << " total_paths=" << count_to_string(total)
              << " max_omega=" << count_to_string(max_count)
              << " pow2=" << pow2
              << " passes_pow2=" << (max_count <= pow2 ? "yes" : "no")
              << std::endl;
    std::cout << "BEST m=" << m
              << " partition=" << part_key(parts[best.part])
              << " strips=" << best.strips << std::endl;
  }

  std::cout << "__EXIT_STATUS=0" << std::endl;
  return 0;
}
