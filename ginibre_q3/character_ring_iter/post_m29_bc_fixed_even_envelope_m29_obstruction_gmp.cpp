#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <unistd.h>

namespace {

using Cell = std::pair<int, int>;
using Part = std::vector<int>;

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

int part_size(const Part& part) {
    int out = 0;
    for (int x : part) out += x;
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

bool is_column_even(const Part& part) {
    if (part.empty()) return true;
    const int width = part.front();
    for (int c = 1; c <= width; ++c) {
        int height = 0;
        for (int row_length : part) {
            if (row_length >= c) ++height;
        }
        if (height % 2 != 0) return false;
    }
    return true;
}

mpz_class top_value(const Part& part, Cell e,
                    const std::unordered_map<std::string, mpz_class>& envelope) {
    if (!has_cell(part, e)) return 0;

    mpz_class total = 0;
    const std::vector<Cell> corners = removable_corners(part);
    for (Cell y : corners) {
        if (y == e) continue;
        Part after_y = remove_cell(part, y);
        if (!is_removable(after_y, e)) continue;
        Part child = remove_cell(after_y, e);
        const auto it = envelope.find(part_key(child));
        if (it == envelope.end()) {
            std::cerr << "missing child envelope for " << part_display(child) << "\n";
            std::exit(3);
        }
        total += it->second;
    }
    return total;
}

mpz_class top_two_value(const Part& part, Cell e, Cell f,
                        const std::unordered_map<std::string, mpz_class>& envelope) {
    if (e == f || !has_cell(part, e) || !has_cell(part, f)) return 0;

    mpz_class total = 0;
    const std::vector<Cell> corners = removable_corners(part);
    for (Cell y : corners) {
        if (y == e || y == f) continue;
        Part after_y = remove_cell(part, y);
        if (!is_removable(after_y, e)) continue;
        Part child = remove_cell(after_y, e);
        if (!has_cell(child, f)) continue;
        total += top_value(child, f, envelope);
    }
    return total;
}

mpz_class top_three_value(const Part& part, Cell e, Cell f, Cell g,
                          const std::unordered_map<std::string, mpz_class>& envelope) {
    if (e == f || e == g || f == g ||
        !has_cell(part, e) || !has_cell(part, f) || !has_cell(part, g)) {
        return 0;
    }

    mpz_class total = 0;
    const std::vector<Cell> corners = removable_corners(part);
    for (Cell y : corners) {
        if (y == e || y == f || y == g) continue;
        Part after_y = remove_cell(part, y);
        if (!is_removable(after_y, e)) continue;
        Part child = remove_cell(after_y, e);
        if (!has_cell(child, f) || !has_cell(child, g)) continue;
        total += top_two_value(child, f, g, envelope);
    }
    return total;
}

}  // namespace

int main() {
    constexpr int target_m = 29;
    constexpr int target_n = 2 * target_m - 1;
    constexpr std::size_t expected_over_budget_count = 144;
    constexpr std::size_t expected_top_pair_over_budget_count = 346;
    constexpr std::size_t expected_top_two_over_budget_count = 62;
    constexpr std::size_t expected_top_three_over_budget_count = 0;
    const Part target{11, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1};
    const mpz_class expected_target_w("367676108829086275");
    const mpz_class budget = pow2(2 * target_m);

    char hostname[256] = {};
    if (gethostname(hostname, sizeof(hostname) - 1) != 0) hostname[0] = '\0';

    std::cout << "B/C fixed-even recursive-envelope m=29 obstruction GMP verifier\n"
              << "VERIFIES the recursive envelope of "
              << "def:bc-fixed-even-recursive-envelope and the m=29 over-budget frontier\n"
              << "host=" << (hostname[0] == '\0' ? "unknown" : hostname) << "\n"
              << "target_m=" << target_m
              << " target_shape=" << part_display(target)
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;
    std::fflush(stdout);

    std::unordered_map<std::string, mpz_class> envelope;
    int failures = 0;
    bool target_seen = false;
    mpz_class target_w = -1;
    mpz_class m29_max_w = 0;
    Part m29_best_part;
    std::size_t m29_over_budget_count = 0;
    std::size_t m29_top_pair_over_budget_count = 0;
    std::size_t m29_top_two_over_budget_count = 0;
    std::size_t m29_top_three_over_budget_count = 0;
    std::vector<std::pair<Part, mpz_class>> m29_over_budget_shapes;
    std::vector<std::tuple<Part, Cell, mpz_class>> m29_top_pair_over_budget;
    std::vector<std::tuple<Part, Cell, Cell, mpz_class>> m29_top_two_over_budget;
    std::vector<std::tuple<Part, Cell, Cell, Cell, mpz_class>> m29_top_three_over_budget;

    for (int m = 1; m <= target_m; ++m) {
        const int n = 2 * m - 1;
        std::vector<Part> parts = gen_parts(n);
        std::vector<mpz_class> values(parts.size());
        std::vector<std::vector<std::pair<Cell, mpz_class>>> top_pair_values;
        if (m == target_m) top_pair_values.resize(parts.size());

#pragma omp parallel for schedule(dynamic, 128)
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(parts.size()); ++i) {
            const Part& part = parts[static_cast<std::size_t>(i)];
            if (n == 1) {
                values[static_cast<std::size_t>(i)] = 1;
                continue;
            }

            mpz_class best = 0;
            std::vector<std::pair<Cell, mpz_class>> local_top_pair_values;
            const std::vector<Cell> cells = cells_of(part);
            const std::vector<Cell> corners = removable_corners(part);
            for (Cell e : cells) {
                mpz_class total = 0;
                for (Cell y : corners) {
                    if (y == e) continue;
                    Part after_y = remove_cell(part, y);
                    if (!is_removable(after_y, e)) continue;
                    Part child = remove_cell(after_y, e);
                    const auto it = envelope.find(part_key(child));
                    if (it == envelope.end()) {
                        std::cerr << "missing child envelope for " << part_display(child) << "\n";
                        std::exit(3);
                    }
                    total += it->second;
                }
                if (total > best) best = total;
                if (m == target_m && total > budget) {
                    local_top_pair_values.push_back({e, total});
                }
            }
            if (m == target_m) {
                top_pair_values[static_cast<std::size_t>(i)] = std::move(local_top_pair_values);
            }
            values[static_cast<std::size_t>(i)] = best;
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
            if (m == target_m && part == target) {
                target_seen = true;
                target_w = value;
            }
            if (m == target_m && value > budget) {
                ++m29_over_budget_count;
                m29_over_budget_shapes.push_back({part, value});
            }
            if (m == target_m && !top_pair_values[i].empty()) {
                for (const auto& [cell, top_value] : top_pair_values[i]) {
                    ++m29_top_pair_over_budget_count;
                    m29_top_pair_over_budget.push_back({part, cell, top_value});
                }
            }
        }

        if (m == target_m) {
            m29_max_w = max_w;
            m29_best_part = best_part;
        }

        std::cout << "RESULT m=" << m
                  << " partitions=" << parts.size()
                  << " max_W=" << max_w.get_str()
                  << " max_W_bits=" << bit_length(max_w)
                  << " best_shape=" << part_display(best_part)
                  << std::endl;
        std::fflush(stdout);
    }

    if (part_size(target) != target_n) ++failures;
    if (is_column_even(target)) ++failures;
    if (!target_seen) ++failures;
    if (target_w != expected_target_w) ++failures;
    if (target_w <= budget) ++failures;
    if (m29_max_w != target_w) ++failures;
    if (m29_best_part != target) ++failures;
    if (m29_over_budget_count != expected_over_budget_count) ++failures;
    if (m29_over_budget_shapes.size() != m29_over_budget_count) ++failures;
    if (m29_top_pair_over_budget_count != expected_top_pair_over_budget_count) ++failures;
    if (m29_top_pair_over_budget.size() != m29_top_pair_over_budget_count) ++failures;
    for (const auto& [part, e, top_value_for_e] : m29_top_pair_over_budget) {
        (void)top_value_for_e;
        for (Cell f : cells_of(part)) {
            if (f == e) continue;
            mpz_class value = top_two_value(part, e, f, envelope);
            if (value > budget) {
                ++m29_top_two_over_budget_count;
                m29_top_two_over_budget.push_back({part, e, f, value});
            }
        }
    }
    if (m29_top_two_over_budget_count != expected_top_two_over_budget_count) ++failures;
    if (m29_top_two_over_budget.size() != m29_top_two_over_budget_count) ++failures;
    for (const auto& [part, e, f, top_two_value_for_ef] : m29_top_two_over_budget) {
        (void)top_two_value_for_ef;
        for (Cell g : cells_of(part)) {
            if (g == e || g == f) continue;
            mpz_class value = top_three_value(part, e, f, g, envelope);
            if (value > budget) {
                ++m29_top_three_over_budget_count;
                m29_top_three_over_budget.push_back({part, e, f, g, value});
            }
        }
    }
    if (m29_top_three_over_budget_count != expected_top_three_over_budget_count) ++failures;
    if (m29_top_three_over_budget.size() != m29_top_three_over_budget_count) ++failures;
    const bool target_in_over_budget =
        std::any_of(m29_over_budget_shapes.begin(), m29_over_budget_shapes.end(),
                    [&](const auto& row) {
                        return row.first == target && row.second == target_w;
                    });
    const bool target_has_top_pair_over_budget =
        std::any_of(m29_top_pair_over_budget.begin(), m29_top_pair_over_budget.end(),
                    [&](const auto& row) {
                        return std::get<0>(row) == target && std::get<2>(row) == target_w;
                    });
    if (!target_in_over_budget) {
        ++failures;
    }
    if (!target_has_top_pair_over_budget) {
        ++failures;
    }

    std::cout << "TARGET size=" << part_size(target)
              << " column_even=" << (is_column_even(target) ? "yes" : "no")
              << " W=" << target_w.get_str()
              << " expected_W=" << expected_target_w.get_str()
              << " budget=2^" << (2 * target_m)
              << " budget_value=" << budget.get_str()
              << " exceeds_budget=" << (target_w > budget ? "yes" : "no")
              << " target_is_m29_maximizer=" << (m29_max_w == target_w && m29_best_part == target ? "yes" : "no")
              << std::endl;
    std::cout << "M29_OVER_BUDGET count=" << m29_over_budget_count
              << " expected_count=" << expected_over_budget_count
              << " target_in_over_budget=" << (target_in_over_budget ? "yes" : "no")
              << "\n";
    std::cout << "M29_TOP_PAIR_OVER_BUDGET count=" << m29_top_pair_over_budget_count
              << " expected_count=" << expected_top_pair_over_budget_count
              << " target_has_top_pair_over_budget=" << (target_has_top_pair_over_budget ? "yes" : "no")
              << "\n";
    std::cout << "M29_TOP_TWO_OVER_BUDGET count=" << m29_top_two_over_budget_count
              << " expected_count=" << expected_top_two_over_budget_count
              << "\n";
    std::cout << "M29_TOP_THREE_OVER_BUDGET count=" << m29_top_three_over_budget_count
              << " expected_count=" << expected_top_three_over_budget_count
              << "\n";
    for (const auto& [part, value] : m29_over_budget_shapes) {
        std::cout << "M29_OVER_BUDGET_SHAPE shape=" << part_display(part)
                  << " W=" << value.get_str()
                  << " exceeds_budget=yes\n";
    }
    for (const auto& [part, cell, value] : m29_top_pair_over_budget) {
        std::cout << "M29_TOP_PAIR_OVER_BUDGET shape=" << part_display(part)
                  << " e=" << cell_display(cell)
                  << " W_top=" << value.get_str()
                  << " exceeds_budget=yes\n";
    }
    for (const auto& [part, e, f, value] : m29_top_two_over_budget) {
        std::cout << "M29_TOP_TWO_OVER_BUDGET shape=" << part_display(part)
                  << " e=" << cell_display(e)
                  << " f=" << cell_display(f)
                  << " W_top2=" << value.get_str()
                  << " exceeds_budget=yes\n";
    }
    for (const auto& [part, e, f, g, value] : m29_top_three_over_budget) {
        std::cout << "M29_TOP_THREE_OVER_BUDGET shape=" << part_display(part)
                  << " e=" << cell_display(e)
                  << " f=" << cell_display(f)
                  << " g=" << cell_display(g)
                  << " W_top3=" << value.get_str()
                  << " exceeds_budget=yes\n";
    }
    std::cout << "NOTE obstruction applies only to this recursive envelope, not to actual broad-Omega cardinalities\n";
    std::cout << "SUMMARY failures=" << failures << std::endl;
    std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << std::endl;
    return failures == 0 ? 0 : 1;
}
