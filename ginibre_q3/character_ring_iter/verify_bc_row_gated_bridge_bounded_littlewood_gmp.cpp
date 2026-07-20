// Fast exact replacement for the state-expansion B/C row-gated bridge replay.
// Bounded-Littlewood determinants are evaluated independently modulo 32-bit
// primes, then uniquely reconstructed before the original 870 Chain checks.

#include <array>

#define GINIBRE_Q3_FULL_Q3_BCD_REMAINING_DATA_HPP
namespace full_q3_bcd_remaining {

enum class TailMethod : unsigned char {
    polynomial,
    rational_cap,
    directed_interval
};

struct RowCutoff {
    char family;
    int rank;
    int tail_onset;
    int moment_through;
    TailMethod tail_method;
};

inline constexpr std::array<RowCutoff, 37> row_cutoffs{{
    {'B', 2, 63, 61, TailMethod::directed_interval},
    {'B', 3, 63, 61, TailMethod::directed_interval},
    {'B', 4, 63, 61, TailMethod::directed_interval},
    {'B', 5, 49, 47, TailMethod::directed_interval},
    {'B', 6, 47, 45, TailMethod::directed_interval},
    {'B', 7, 45, 43, TailMethod::directed_interval},
    {'B', 8, 45, 43, TailMethod::directed_interval},
    {'B', 9, 43, 41, TailMethod::directed_interval},
    {'B', 10, 43, 41, TailMethod::directed_interval},
    {'B', 11, 43, 51, TailMethod::directed_interval},
    {'B', 12, 43, 53, TailMethod::directed_interval},
    {'B', 13, 43, 57, TailMethod::directed_interval},
    {'B', 14, 43, 41, TailMethod::directed_interval},
    {'B', 15, 43, 41, TailMethod::directed_interval},
    {'B', 16, 43, 41, TailMethod::directed_interval},
    {'B', 17, 43, 41, TailMethod::directed_interval},
    {'B', 18, 43, 41, TailMethod::directed_interval},
    {'B', 19, 43, 41, TailMethod::directed_interval},
    {'C', 2, 63, 61, TailMethod::directed_interval},
    {'C', 3, 63, 61, TailMethod::directed_interval},
    {'C', 4, 63, 61, TailMethod::directed_interval},
    {'C', 5, 49, 47, TailMethod::directed_interval},
    {'C', 6, 47, 45, TailMethod::directed_interval},
    {'C', 7, 45, 43, TailMethod::directed_interval},
    {'C', 8, 45, 43, TailMethod::directed_interval},
    {'C', 9, 43, 41, TailMethod::directed_interval},
    {'C', 10, 43, 41, TailMethod::directed_interval},
    {'C', 11, 43, 41, TailMethod::directed_interval},
    {'C', 12, 43, 41, TailMethod::directed_interval},
    {'C', 13, 43, 43, TailMethod::directed_interval},
    {'C', 14, 43, 51, TailMethod::directed_interval},
    {'C', 15, 43, 41, TailMethod::directed_interval},
    {'C', 16, 43, 41, TailMethod::directed_interval},
    {'C', 17, 43, 45, TailMethod::directed_interval},
    {'C', 18, 43, 41, TailMethod::directed_interval},
    {'C', 19, 43, 41, TailMethod::directed_interval},
    // The shared determinant engine constructs all three classical matrices;
    // this degree-zero sentinel keeps its D matrix nonempty and is not a
    // bridge row.
    {'D', 4, 1, 0, TailMethod::directed_interval},
}};

inline constexpr int required_maximum_moment = 61;
inline constexpr std::size_t polynomial_rows = 0;
inline constexpr std::size_t rational_cap_rows = 0;
inline constexpr std::size_t directed_interval_rows = row_cutoffs.size();
inline constexpr std::size_t full_residual_pairs = 0;

inline const RowCutoff* find_row(char family, int rank) {
    for (const RowCutoff& row : row_cutoffs) {
        if (row.family == family && row.rank == rank) return &row;
    }
    return nullptr;
}

}  // namespace full_q3_bcd_remaining

#define main full_q3_bounded_littlewood_main_disabled
#include "verify_full_q3_bcd_bounded_littlewood_gmp.cpp"
#undef main

#include <fstream>
#include <map>
#include <optional>
#include <regex>
#include <set>
#include <string>
#include <tuple>

namespace fast_bc_bridge {

using Interval = std::pair<BigInt, BigInt>;
using ClaimKey = std::tuple<char, int, int>;

std::map<ClaimKey, BigInt> read_source_claims(
    const std::vector<std::string>& paths
) {
    const std::regex delta_pattern(
        R"(([BC])_(\d+).*?Delta_(\d+)\s*=\s*(-?\d+))"
    );
    std::map<ClaimKey, BigInt> claims;
    for (const std::string& path : paths) {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("cannot open B/C source log: " + path);
        std::string line;
        while (std::getline(input, line)) {
            for (std::sregex_iterator it(line.begin(), line.end(), delta_pattern), end;
                 it != end; ++it) {
                const ClaimKey key{
                    (*it)[1].str()[0],
                    std::stoi((*it)[2].str()),
                    std::stoi((*it)[3].str()),
                };
                const BigInt value((*it)[4].str());
                const auto [found, inserted] = claims.emplace(key, value);
                if (!inserted && found->second != value) {
                    throw std::runtime_error("conflicting B/C source-log claims");
                }
            }
        }
    }
    return claims;
}

void compare_stable_table(const std::string& path, const std::vector<BigInt>& stable) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("could not open stable table: " + path);
    std::vector<bool> seen(stable.size(), false);
    int degree = -1;
    std::string value;
    while (input >> degree >> value) {
        if (degree < 0 || degree >= static_cast<int>(stable.size())) continue;
        if (seen[static_cast<std::size_t>(degree)]) {
            throw std::runtime_error("duplicate stable-table degree");
        }
        seen[static_cast<std::size_t>(degree)] = true;
        if (BigInt(value) != stable[static_cast<std::size_t>(degree)]) {
            throw std::runtime_error("stable recurrence/table mismatch");
        }
    }
    for (bool present : seen) {
        if (!present) throw std::runtime_error("stable table omits a required degree");
    }
}

BigInt q3(const std::vector<BigInt>& moments, int n) {
    BigInt total = 0;
    for (int k = 0; k <= n; ++k) {
        total += binomial(n, k)
            * (moments[static_cast<std::size_t>(k + 2)]
                   * moments[static_cast<std::size_t>(n - k)]
               - moments[static_cast<std::size_t>(k + 1)]
                   * moments[static_cast<std::size_t>(n - k + 1)]);
    }
    return 2 * total;
}

BigInt chain_difference(const std::vector<BigInt>& moments, int m) {
    return q3(moments, 2 * m + 3) - 4 * q3(moments, 2 * m + 1);
}

std::vector<BigInt> linear_coefficients(
    const std::vector<BigInt>& moments,
    int m
) {
    const int maximum = 2 * m + 5;
    std::vector<BigInt> answer(static_cast<std::size_t>(maximum + 1));
    for (const auto& [n, scale] :
         std::array<std::pair<int, int>, 2>{{{2 * m + 3, 1}, {2 * m + 1, -4}}}) {
        for (int k = 0; k <= n; ++k) {
            const BigInt coefficient = 2 * scale * binomial(n, k);
            int first = k + 2;
            int second = n - k;
            answer[static_cast<std::size_t>(first)] +=
                coefficient * moments[static_cast<std::size_t>(second)];
            answer[static_cast<std::size_t>(second)] +=
                coefficient * moments[static_cast<std::size_t>(first)];
            first = k + 1;
            second = n - k + 1;
            answer[static_cast<std::size_t>(first)] -=
                coefficient * moments[static_cast<std::size_t>(second)];
            answer[static_cast<std::size_t>(second)] -=
                coefficient * moments[static_cast<std::size_t>(first)];
        }
    }
    return answer;
}

std::map<std::pair<int, int>, BigInt> quadratic_coefficients(int m) {
    std::map<std::pair<int, int>, BigInt> answer;
    for (const auto& [n, scale] :
         std::array<std::pair<int, int>, 2>{{{2 * m + 3, 1}, {2 * m + 1, -4}}}) {
        for (int k = 0; k <= n; ++k) {
            const BigInt coefficient = 2 * scale * binomial(n, k);
            for (const auto& [left_raw, right_raw, sign] :
                 std::array<std::tuple<int, int, int>, 2>{{
                     {k + 2, n - k, 1}, {k + 1, n - k + 1, -1}}}) {
                const int left = std::min(left_raw, right_raw);
                const int right = std::max(left_raw, right_raw);
                answer[{left, right}] += sign * coefficient;
            }
        }
    }
    return answer;
}

BigInt product_lower(const Interval& left, const Interval& right, const BigInt& coefficient) {
    const std::array<BigInt, 4> products{{
        left.first * right.first,
        left.first * right.second,
        left.second * right.first,
        left.second * right.second,
    }};
    return coefficient >= 0
        ? coefficient * *std::min_element(products.begin(), products.end())
        : coefficient * *std::max_element(products.begin(), products.end());
}

const MomentRow* find_row(const std::vector<MomentRow>& rows, char family, int rank) {
    for (const MomentRow& row : rows) {
        if (row.family == family && row.rank == rank) return &row;
    }
    return nullptr;
}

int bridge_moment_through(char family, int rank) {
    if (family == 'B') {
        if (rank <= 4) return 61;
        if (rank == 5) return 47;
        if (rank == 6) return 45;
        if (rank <= 8) return 43;
        return 41;
    }
    if (family == 'C') {
        if (rank <= 4) return 61;
        if (rank == 5) return 47;
        if (rank == 6) return 45;
        if (rank <= 8) return 43;
        return 41;
    }
    return 0;
}

std::optional<BigInt> exact_delta(
    const std::vector<MomentRow>& rows,
    const std::vector<BigInt>& stable,
    char family,
    int rank,
    int degree
) {
    if (degree < 2 * rank + 2) return BigInt(0);
    if (rank >= 20) return std::nullopt;
    const MomentRow* row = find_row(rows, family, rank);
    if (row == nullptr || degree > bridge_moment_through(family, rank)) {
        return std::nullopt;
    }
    return row->moments[static_cast<std::size_t>(degree)]
        - stable[static_cast<std::size_t>(degree)];
}

void verify_bridge(const std::vector<MomentRow>& rows, const std::vector<BigInt>& stable) {
    int row_count = 0;
    int steps = 0;
    int exact_steps = 0;
    int interval_steps = 0;
    int correction_entries = 0;
    BigInt correction_checksum = 0;
    BigInt ledger_checksum = 0;
    BigInt minimum;
    std::string minimum_label;
    bool have_minimum = false;

    for (const MomentRow& row : rows) {
        const int boundary = 2 * row.rank + 2;
        const int bridge_maximum = bridge_moment_through(row.family, row.rank);
        for (int degree = 0; degree <= bridge_maximum; ++degree) {
            const BigInt& value = row.moments[static_cast<std::size_t>(degree)];
            if (degree < boundary && value != stable[static_cast<std::size_t>(degree)]) {
                throw std::runtime_error("pre-boundary stabilization mismatch");
            }
            if (degree >= boundary) {
                const BigInt delta = value - stable[static_cast<std::size_t>(degree)];
                if (delta > 0 || delta < -stable[static_cast<std::size_t>(degree)]) {
                    throw std::runtime_error("finite correction left [-s_j,0]");
                }
                correction_checksum +=
                    BigInt((row.family == 'B' ? 1 : 2) * 1000000
                           + row.rank * 1000 + degree)
                    * delta;
                ++correction_entries;
            }
        }
    }
    if (correction_entries != 832) {
        throw std::runtime_error("generated correction-entry count mismatch");
    }

    for (char family : std::array<char, 2>{{'B', 'C'}}) {
        for (int rank = 2; rank <= 30; ++rank) {
            ++row_count;
            for (int m = rank - 1; m <= 29; ++m) {
                ++steps;
                const int maximum = 2 * m + 5;
                BigInt margin = chain_difference(stable, m);
                if (margin <= 0) throw std::runtime_error("stable Chain difference is not positive");
                const auto linear = linear_coefficients(stable, m);
                const auto quadratic = quadratic_coefficients(m);
                std::vector<Interval> intervals(static_cast<std::size_t>(maximum + 1));
                bool all_exact = true;
                for (int degree = 0; degree <= maximum; ++degree) {
                    const auto delta = exact_delta(rows, stable, family, rank, degree);
                    if (delta.has_value()) {
                        intervals[static_cast<std::size_t>(degree)] = {*delta, *delta};
                        ledger_checksum +=
                            BigInt((family == 'B' ? 1 : 2) * 1000000
                                   + rank * 1000 + degree)
                            * (*delta);
                    } else {
                        all_exact = false;
                        intervals[static_cast<std::size_t>(degree)] = {
                            -stable[static_cast<std::size_t>(degree)], 0};
                    }
                }
                for (int degree = 0; degree <= maximum; ++degree) {
                    const Interval& interval = intervals[static_cast<std::size_t>(degree)];
                    margin += linear[static_cast<std::size_t>(degree)]
                        * (linear[static_cast<std::size_t>(degree)] >= 0
                               ? interval.first : interval.second);
                }
                for (const auto& [indices, coefficient] : quadratic) {
                    margin += product_lower(
                        intervals[static_cast<std::size_t>(indices.first)],
                        intervals[static_cast<std::size_t>(indices.second)],
                        coefficient);
                }
                if (all_exact) {
                    std::vector<BigInt> moment_row(
                        stable.begin(), stable.begin() + maximum + 1
                    );
                    for (int degree = 0; degree <= maximum; ++degree) {
                        moment_row[static_cast<std::size_t>(degree)] +=
                            intervals[static_cast<std::size_t>(degree)].first;
                    }
                    if (chain_difference(moment_row, m) != margin) {
                        throw std::runtime_error("exact expansion/direct Chain mismatch");
                    }
                    ++exact_steps;
                } else {
                    ++interval_steps;
                }
                if (margin <= 0) throw std::runtime_error("nonpositive row-gated margin");
                const std::string label = std::string(1, family) + "_"
                    + std::to_string(rank) + " m=" + std::to_string(m);
                if (!have_minimum || margin < minimum) {
                    have_minimum = true;
                    minimum = margin;
                    minimum_label = label;
                }
            }
        }
    }

    const BigInt expected_checksum(
        "-522736234492485419094885313899428626303654826589487960349093081180846922600009362789711747"
    );
    if (row_count != 58 || steps != 870 || exact_steps != 416
        || interval_steps != 454 || minimum_label != "B_2 m=1"
        || minimum != 1046 || ledger_checksum != expected_checksum) {
        throw std::runtime_error("row-gated scope/checksum regression");
    }
    std::cout << "generated_correction_entries=" << correction_entries
              << " correction_checksum=" << correction_checksum << '\n';
    std::cout << "minimum_positive_margin_label=" << minimum_label
              << " minimum_positive_margin=" << minimum << '\n';
    std::cout << "ledger_checksum=" << ledger_checksum << '\n';
    std::cout << "families=2 ranks=58 row_gated_steps=870 exact_steps=416 "
                 "interval_steps=454 failures=0 status=PASS\n";
}

void verify_source_claims(
    const std::vector<MomentRow>& rows,
    const std::vector<BigInt>& stable,
    const std::vector<std::string>& source_paths
) {
    const auto claims = read_source_claims(source_paths);
    std::set<std::pair<char, int>> groups;
    for (const auto& [key, expected] : claims) {
        const auto [family, rank, degree] = key;
        const MomentRow* row = find_row(rows, family, rank);
        if (row == nullptr || degree < 0 || degree > row->moment_through) {
            throw std::runtime_error("B/C source claim is outside reconstructed scope");
        }
        const BigInt actual = row->moments[static_cast<std::size_t>(degree)]
            - stable[static_cast<std::size_t>(degree)];
        if (actual != expected) {
            throw std::runtime_error("B/C source claim mismatch");
        }
        groups.emplace(family, rank);
    }
    if (claims.size() != 182 || groups.size() != 14) {
        throw std::runtime_error("B/C source-claim scope regression");
    }
    std::cout << "source_claims_checked=" << claims.size()
              << " source_groups_checked=" << groups.size() << '\n';
}

}  // namespace fast_bc_bridge

int main(int argc, char** argv) {
    try {
        // The included generic engine also supplies exact-tail routines used
        // by its ordinary entry point; retain compile-time warning coverage
        // although this bridge ledger has no tail rows.
        const auto polynomial_tail_routine = &verify_polynomial_tail;
        const auto rational_tail_routine = &verify_rational_cap_tail;
        (void)polynomial_tail_routine;
        (void)rational_tail_routine;
        bool progress = false;
        std::string stable_path = "../references/oeis_A002137_stable.txt";
        std::vector<std::string> source_paths;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--stable" && index + 1 < argc) {
                stable_path = argv[++index];
            } else if (argument == "--progress") {
                progress = true;
            } else if (argument == "--source-log" && index + 1 < argc) {
                source_paths.emplace_back(argv[++index]);
            } else {
                throw std::runtime_error(
                    "usage: verifier [--stable PATH] [--source-log PATH]... [--progress]"
                );
            }
        }
        std::cout << "BC_ROW_GATED_BOUNDED_LITTLEWOOD maximum=61";
#ifdef _OPENMP
        std::cout << " parallel_threads=" << omp_get_max_threads();
#endif
        std::cout << '\n';
        const auto stable = stable_moments(63);
        fast_bc_bridge::compare_stable_table(stable_path, stable);
        const auto rows = reconstruct_bounded_littlewood_moments(61, progress);
        if (source_paths.empty()) {
            throw std::runtime_error("no B/C source logs supplied");
        }
        fast_bc_bridge::verify_source_claims(rows, stable, source_paths);
        fast_bc_bridge::verify_bridge(rows, stable);
        std::cout << "BC_ROW_GATED_BOUNDED_LITTLEWOOD VERIFICATION: ALL PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "BC_ROW_GATED_BOUNDED_LITTLEWOOD ERROR: " << error.what() << '\n';
        return 1;
    }
}
