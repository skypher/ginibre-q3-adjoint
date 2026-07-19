#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <unistd.h>

namespace {

using Cell = std::pair<int, int>;
using Part = std::vector<int>;

struct Prefix {
    Part part;
    std::vector<Cell> fixed;
    std::uint64_t value;
};

class CappedMemo {
   public:
    explicit CappedMemo(std::size_t shard_count)
        : maps_(shard_count), locks_(shard_count) {
        for (omp_lock_t& lock : locks_) omp_init_lock(&lock);
    }

    ~CappedMemo() {
        for (omp_lock_t& lock : locks_) omp_destroy_lock(&lock);
    }

    bool get(const std::string& key, std::uint64_t& value) {
        const std::size_t shard = shard_for(key);
        omp_set_lock(&locks_[shard]);
        const auto it = maps_[shard].find(key);
        if (it == maps_[shard].end()) {
            omp_unset_lock(&locks_[shard]);
            return false;
        }
        value = it->second;
        omp_unset_lock(&locks_[shard]);
        return true;
    }

    void put(const std::string& key, std::uint64_t value) {
        const std::size_t shard = shard_for(key);
        omp_set_lock(&locks_[shard]);
        maps_[shard][key] = value;
        omp_unset_lock(&locks_[shard]);
    }

    std::size_t size_approx() {
        std::size_t total = 0;
        for (std::size_t i = 0; i < maps_.size(); ++i) {
            omp_set_lock(&locks_[i]);
            total += maps_[i].size();
            omp_unset_lock(&locks_[i]);
        }
        return total;
    }

   private:
    std::size_t shard_for(const std::string& key) const {
        return std::hash<std::string>{}(key) % maps_.size();
    }

    std::vector<std::unordered_map<std::string, std::uint64_t>> maps_;
    std::vector<omp_lock_t> locks_;
};

std::string part_key(const Part& part) {
    std::ostringstream out;
    for (std::size_t i = 0; i < part.size(); ++i) {
        if (i) out << ',';
        out << part[i];
    }
    return out.str();
}

std::string part_display(const Part& part) {
    std::ostringstream out;
    out << '(';
    for (std::size_t i = 0; i < part.size(); ++i) {
        if (i) out << ',';
        out << part[i];
    }
    out << ')';
    return out.str();
}

std::string cell_display(Cell cell) {
    std::ostringstream out;
    out << '(' << cell.first << ',' << cell.second << ')';
    return out.str();
}

std::string fixed_display(const std::vector<Cell>& fixed) {
    std::ostringstream out;
    out << '[';
    for (std::size_t i = 0; i < fixed.size(); ++i) {
        if (i) out << ',';
        out << cell_display(fixed[i]);
    }
    out << ']';
    return out.str();
}

std::string state_key(const Part& part, const std::vector<Cell>& fixed) {
    std::string out;
    out.reserve(2 + part.size() + 2 * fixed.size());
    out.push_back(static_cast<char>(part.size()));
    for (int row_length : part) out.push_back(static_cast<char>(row_length));
    out.push_back(static_cast<char>(fixed.size()));
    for (Cell cell : fixed) {
        out.push_back(static_cast<char>(cell.first));
        out.push_back(static_cast<char>(cell.second));
    }
    return out;
}

void gen_parts_rec(int n, int max_part, Part& current, std::vector<Part>& out) {
    if (n == 0) {
        out.push_back(current);
        return;
    }
    for (int x = std::min(n, max_part); x >= 1; --x) {
        current.push_back(x);
        gen_parts_rec(n - x, x, current, out);
        current.pop_back();
    }
}

std::vector<Part> gen_parts(int n) {
    std::vector<Part> out;
    Part current;
    gen_parts_rec(n, n, current, out);
    return out;
}

bool has_cell(const Part& part, Cell cell) {
    const int r = cell.first;
    const int c = cell.second;
    return 1 <= r && r <= static_cast<int>(part.size()) &&
           1 <= c && c <= part[static_cast<std::size_t>(r - 1)];
}

std::vector<Cell> cells_of(const Part& part) {
    std::vector<Cell> out;
    for (int r = 1; r <= static_cast<int>(part.size()); ++r) {
        for (int c = 1; c <= part[static_cast<std::size_t>(r - 1)]; ++c) {
            out.push_back({r, c});
        }
    }
    return out;
}

bool is_removable(const Part& part, Cell cell) {
    if (!has_cell(part, cell)) return false;
    const int r = cell.first;
    const int c = cell.second;
    if (part[static_cast<std::size_t>(r - 1)] != c) return false;
    if (r < static_cast<int>(part.size()) &&
        part[static_cast<std::size_t>(r)] >= c) {
        return false;
    }
    return true;
}

std::vector<Cell> removable_corners(const Part& part) {
    std::vector<Cell> out;
    for (int r = 1; r <= static_cast<int>(part.size()); ++r) {
        Cell cell{r, part[static_cast<std::size_t>(r - 1)]};
        if (is_removable(part, cell)) out.push_back(cell);
    }
    return out;
}

Part remove_cell(Part part, Cell cell) {
    --part[static_cast<std::size_t>(cell.first - 1)];
    while (!part.empty() && part.back() == 0) part.pop_back();
    return part;
}

mpz_class pow2(int exponent) {
    mpz_class out;
    mpz_ui_pow_ui(out.get_mpz_t(), 2UL, static_cast<unsigned long>(exponent));
    return out;
}

int bit_length(const mpz_class& value) {
    if (sgn(value) <= 0) return 0;
    return static_cast<int>(mpz_sizeinbase(value.get_mpz_t(), 2));
}

bool contains_cell(const std::vector<Cell>& fixed, Cell cell) {
    return std::find(fixed.begin(), fixed.end(), cell) != fixed.end();
}

mpz_class envelope_value(const Part& part,
                         const std::unordered_map<std::string, mpz_class>& envelope) {
    const auto it = envelope.find(part_key(part));
    if (it == envelope.end()) {
        std::cerr << "missing envelope for " << part_display(part) << "\n";
        std::exit(3);
    }
    return it->second;
}

mpz_class fixed_top_value(const Part& part, const std::vector<Cell>& fixed,
                          const std::unordered_map<std::string, mpz_class>& envelope) {
    if (fixed.empty()) return envelope_value(part, envelope);

    const Cell e = fixed.front();
    if (!has_cell(part, e)) return 0;
    for (Cell cell : fixed) {
        if (!has_cell(part, cell)) return 0;
    }

    std::vector<Cell> rest(fixed.begin() + 1, fixed.end());
    mpz_class total = 0;
    for (Cell y : removable_corners(part)) {
        if (contains_cell(fixed, y)) continue;
        Part after_y = remove_cell(part, y);
        if (!is_removable(after_y, e)) continue;
        Part child = remove_cell(after_y, e);
        bool rest_survives = true;
        for (Cell cell : rest) {
            if (!has_cell(child, cell)) {
                rest_survives = false;
                break;
            }
        }
        if (!rest_survives) continue;
        total += fixed_top_value(child, rest, envelope);
    }
    return total;
}

std::uint64_t fixed_top_value_bounded(
    const Part& part, const std::vector<Cell>& fixed,
    const std::unordered_map<std::string, mpz_class>& envelope,
    const mpz_class& budget, std::uint64_t budget_u64, std::uint64_t cap,
    CappedMemo& memo) {
    if (fixed.empty()) {
        const mpz_class value = envelope_value(part, envelope);
        return value > budget ? cap : value.get_ui();
    }

    const std::string key = state_key(part, fixed);
    std::uint64_t cached = 0;
    if (memo.get(key, cached)) return cached;

    const Cell e = fixed.front();
    if (!has_cell(part, e)) return 0;
    for (Cell cell : fixed) {
        if (!has_cell(part, cell)) return 0;
    }

    std::vector<Cell> rest(fixed.begin() + 1, fixed.end());
    std::uint64_t total = 0;
    for (Cell y : removable_corners(part)) {
        if (contains_cell(fixed, y)) continue;
        Part after_y = remove_cell(part, y);
        if (!is_removable(after_y, e)) continue;
        Part child = remove_cell(after_y, e);
        bool rest_survives = true;
        for (Cell cell : rest) {
            if (!has_cell(child, cell)) {
                rest_survives = false;
                break;
            }
        }
        if (!rest_survives) continue;

        total += fixed_top_value_bounded(
            child, rest, envelope, budget, budget_u64, cap, memo);
        if (total > budget_u64) {
            memo.put(key, cap);
            return cap;
        }
    }
    memo.put(key, total);
    return total;
}

std::vector<Prefix> extend_frontier(const std::vector<Prefix>& frontier,
                                    const mpz_class& budget,
                                    std::uint64_t budget_u64,
                                    std::uint64_t cap,
                                    const std::unordered_map<std::string, mpz_class>& envelope,
                                    CappedMemo& memo) {
    const int threads = omp_get_max_threads();
    std::vector<std::vector<Prefix>> locals(static_cast<std::size_t>(threads));

#pragma omp parallel for schedule(dynamic, 32)
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(frontier.size()); ++i) {
        const int tid = omp_get_thread_num();
        const Prefix& prefix = frontier[static_cast<std::size_t>(i)];
        for (Cell cell : cells_of(prefix.part)) {
            if (contains_cell(prefix.fixed, cell)) continue;
            std::vector<Cell> next_fixed = prefix.fixed;
            next_fixed.push_back(cell);
            std::uint64_t value = fixed_top_value_bounded(
                prefix.part, next_fixed, envelope, budget, budget_u64, cap, memo);
            if (value > budget_u64) {
                locals[static_cast<std::size_t>(tid)].push_back(
                    {prefix.part, std::move(next_fixed), value});
            }
        }
    }

    std::vector<Prefix> out;
    for (auto& local : locals) {
        out.insert(out.end(),
                   std::make_move_iterator(local.begin()),
                   std::make_move_iterator(local.end()));
    }
    return out;
}

}  // namespace

int main() {
    constexpr int target_m = 31;
    constexpr int target_n = 2 * target_m - 1;
    constexpr int max_depth = 14;
    constexpr bool assert_expected_counts = true;
    constexpr std::size_t expected_over_budget_shapes = 67038;
    const std::vector<std::size_t> expected_depth_counts{
        331969, 1229238, 2687932, 2263824, 499606,
        49486, 6142, 610, 24, 0};
    const Part expected_best{12, 10, 9, 7, 6, 5, 4, 3, 2, 2, 1};
    const mpz_class expected_max_w("22734772942354066748");
    const mpz_class budget = pow2(2 * target_m);
    const std::uint64_t budget_u64 = 1ULL << (2 * target_m);
    const std::uint64_t cap = budget_u64 + 1;
    int failures = 0;

    char hostname[256] = {};
    if (gethostname(hostname, sizeof(hostname) - 1) != 0) hostname[0] = '\0';

    std::cout << "B/C fixed-even recursive-envelope m=31 frontier GMP verifier\n"
              << "host=" << (hostname[0] == '\0' ? "unknown" : hostname) << "\n"
              << "target_m=" << target_m
              << " target_n=" << target_n
              << " budget=2^" << 2 * target_m
              << " budget_value=" << budget.get_str()
              << " cap_value=" << cap
              << " max_depth=" << max_depth
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;

    std::unordered_map<std::string, mpz_class> envelope;
    CappedMemo capped_memo(4096);
    std::vector<Prefix> frontier;
    std::size_t m31_over_budget_shapes = 0;
    mpz_class m31_max_w = 0;
    Part m31_best_part;

    for (int m = 1; m <= target_m; ++m) {
        const int n = 2 * m - 1;
        std::vector<Part> parts = gen_parts(n);
        std::vector<mpz_class> values(parts.size());
        std::vector<std::vector<Prefix>> top_frontiers;
        if (m == target_m) top_frontiers.resize(parts.size());

#pragma omp parallel for schedule(dynamic, 128)
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(parts.size()); ++i) {
            const Part& part = parts[static_cast<std::size_t>(i)];
            if (n == 1) {
                values[static_cast<std::size_t>(i)] = 1;
                continue;
            }

            mpz_class best = 0;
            std::vector<Prefix> local_frontier;
            for (Cell e : cells_of(part)) {
                mpz_class value = fixed_top_value(part, std::vector<Cell>{e}, envelope);
                if (value > best) best = value;
                if (m == target_m && value > budget) {
                    local_frontier.push_back({part, std::vector<Cell>{e}, cap});
                }
            }
            values[static_cast<std::size_t>(i)] = best;
            if (m == target_m) {
                top_frontiers[static_cast<std::size_t>(i)] = std::move(local_frontier);
            }
        }

        mpz_class max_w = 0;
        Part best_part;
        for (std::size_t i = 0; i < parts.size(); ++i) {
            const Part& part = parts[i];
            const mpz_class& value = values[i];
            envelope[part_key(part)] = value;
            if (value > max_w) {
                max_w = value;
                best_part = part;
            }
            if (m == target_m && value > budget) {
                ++m31_over_budget_shapes;
            }
            if (m == target_m && !top_frontiers[i].empty()) {
                frontier.insert(frontier.end(),
                                std::make_move_iterator(top_frontiers[i].begin()),
                                std::make_move_iterator(top_frontiers[i].end()));
            }
        }

        if (m == target_m) {
            m31_max_w = max_w;
            m31_best_part = best_part;
        }

        std::cout << "RESULT m=" << m
                  << " partitions=" << parts.size()
                  << " max_W=" << max_w.get_str()
                  << " max_W_bits=" << bit_length(max_w)
                  << " best_shape=" << part_display(best_part)
                  << std::endl;
    }

    std::cout << "M31_OVER_BUDGET_SHAPES count=" << m31_over_budget_shapes << "\n"
              << "M31_OVER_BUDGET_SHAPES expected_count="
              << (assert_expected_counts ? std::to_string(expected_over_budget_shapes)
                                          : std::string("unasserted"))
              << "\n"
              << "M31_MAX W=" << m31_max_w.get_str()
              << " expected_W="
              << (assert_expected_counts ? expected_max_w.get_str()
                                          : std::string("unasserted"))
              << " bits=" << bit_length(m31_max_w)
              << " shape=" << part_display(m31_best_part) << std::endl;

    std::vector<std::size_t> depth_counts;
    depth_counts.push_back(frontier.size());
    std::cout << "M31_TOP_DEPTH depth=1 count=" << frontier.size()
              << " expected_count="
              << (assert_expected_counts ? std::to_string(expected_depth_counts[0])
                                          : std::string("unasserted"))
              << std::endl;
    for (int depth = 2; depth <= max_depth && !frontier.empty(); ++depth) {
        frontier = extend_frontier(
            frontier, budget, budget_u64, cap, envelope, capped_memo);
        depth_counts.push_back(frontier.size());
        std::cout << "M31_TOP_DEPTH depth=" << depth
                  << " count=" << frontier.size()
                  << " expected_count="
                  << (assert_expected_counts &&
                              depth <= static_cast<int>(expected_depth_counts.size())
                          ? std::to_string(expected_depth_counts[static_cast<std::size_t>(depth - 1)])
                          : std::string("unasserted"))
                  << " capped_memo_entries=" << capped_memo.size_approx()
                  << std::endl;
    }

    if (m31_max_w <= budget) ++failures;
    if (!frontier.empty()) ++failures;
    if (assert_expected_counts) {
        if (m31_over_budget_shapes != expected_over_budget_shapes) ++failures;
        if (m31_max_w != expected_max_w) ++failures;
        if (m31_best_part != expected_best) ++failures;
        if (depth_counts.size() != expected_depth_counts.size()) {
            ++failures;
        } else {
            for (std::size_t i = 0; i < expected_depth_counts.size(); ++i) {
                if (depth_counts[i] != expected_depth_counts[i]) ++failures;
            }
        }
    }

    const std::size_t display_limit = 80;
    for (std::size_t i = 0; i < frontier.size() && i < display_limit; ++i) {
        const Prefix& row = frontier[i];
        std::cout << "M31_REMAINING_PREFIX shape=" << part_display(row.part)
                  << " fixed=" << fixed_display(row.fixed)
                  << " bounded_W_top" << row.fixed.size() << "=" << row.value
                  << " exceeds_budget=yes\n";
    }
    std::cout << "SUMMARY remaining_count=" << frontier.size()
              << " failures=" << failures << std::endl;
    std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << std::endl;
    std::cout.flush();
    // The memo table can contain hundreds of millions of states; after the
    // final status line it is intentionally left to the operating system.
    std::_Exit(failures == 0 ? 0 : 1);
}
