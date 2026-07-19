#include <gmpxx.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr double LOG10_2 = 0.301029995663981195213738894724493027;

double log10_abs(const mpz_class& value) {
    if (value == 0) return -std::numeric_limits<double>::infinity();
    mpz_class absolute = value >= 0 ? value : -value;
    long exponent = 0;
    const double mantissa = mpz_get_d_2exp(&exponent, absolute.get_mpz_t());
    return std::log10(mantissa) + static_cast<double>(exponent) * LOG10_2;
}

std::vector<mpz_class> stable_moments(int maximum) {
    std::vector<mpz_class> moments(static_cast<std::size_t>(maximum + 1));
    moments[0] = 1;
    if (maximum >= 1) moments[1] = 0;
    for (int index = 1; index < maximum; ++index) {
        mpz_class next = index * moments[static_cast<std::size_t>(index)]
            + index * moments[static_cast<std::size_t>(index - 1)];
        if (index >= 2) {
            next -= (mpz_class(index) * (index - 1) / 2)
                * moments[static_cast<std::size_t>(index - 2)];
        }
        if (next <= 0) {
            std::cerr << "stable recurrence lost positivity at index "
                      << (index + 1) << '\n';
            std::exit(2);
        }
        moments[static_cast<std::size_t>(index + 1)] = std::move(next);
    }
    return moments;
}

std::vector<mpz_class> binomial_row(int n) {
    std::vector<mpz_class> row(static_cast<std::size_t>(n + 1));
    row[0] = 1;
    for (int k = 0; k < n; ++k) {
        row[static_cast<std::size_t>(k + 1)] =
            row[static_cast<std::size_t>(k)] * (n - k) / (k + 1);
    }
    return row;
}

void add_if_valid(
    mpz_class& coefficient,
    const std::vector<mpz_class>& row,
    int index,
    int factor
) {
    if (index < 0 || index >= static_cast<int>(row.size())) return;
    coefficient += factor * row[static_cast<std::size_t>(index)];
}

mpz_class pair_coefficient(
    const std::vector<mpz_class>& row,
    int first,
    int second,
    int scale
) {
    const int factor = 2 * scale;
    mpz_class coefficient = 0;
    add_if_valid(coefficient, row, first - 2, factor);
    if (second != first) add_if_valid(coefficient, row, second - 2, factor);
    add_if_valid(coefficient, row, first - 1, -factor);
    if (second != first) add_if_valid(coefficient, row, second - 1, -factor);
    return coefficient;
}

struct Interval {
    mpz_class lower;
    mpz_class upper;
};

using ExactCorrections = std::map<int, std::map<int, mpz_class>>;

constexpr const char* EXACT_MOMENT_SOURCE_SHA256 =
    "bc22efcafafe2d40b5f5430cb8f9cb046bc6b9b188da684e61dc6fdfaa27077e";
constexpr const char* EXACT_MOMENT_SOURCE_HOST = "nb1cb2f";

ExactCorrections read_exact_correction_logs(
    const std::vector<std::string>& paths,
    const std::vector<mpz_class>& stable
) {
    const std::regex full_pattern(
        R"(^\s*D_(\d+) m_(\d+) \+= O_even (-?\d+) \+ det (-?\d+); Delta_(\d+) = (-?\d+); moment_(\d+) = (\d+)$)"
    );
    const std::regex stable_window_pattern(
        R"(^\s*D_(\d+) Delta_(\d+) = (-?\d+))"
    );
    const std::regex exact_weyl_pattern(
        R"(^\s*D_(\d+) exact_Weyl m_(\d+) = (-?\d+); Delta_(\d+) = (-?\d+))"
    );
    ExactCorrections corrections;
    std::size_t full_rows = 0;
    std::size_t stable_window_rows = 0;
    std::size_t exact_weyl_rows = 0;
    std::size_t ignored_legacy_rows = 0;
    for (const std::string& path : paths) {
        std::ifstream input(path);
        if (!input) {
            std::cerr << "cannot open exact-correction log: " << path << '\n';
            std::exit(2);
        }
        int success_rows = 0;
        int host_rows = 0;
        int source_hash_rows = 0;
        bool found_full_row = false;
        bool found_usable_row = false;
        std::string line;
        while (std::getline(input, line)) {
            if (line == "__EXIT_STATUS=0"
                || line == "__PART_EXIT_STATUS=0") {
                ++success_rows;
            }
            if (line.rfind("__EXIT_STATUS=", 0) == 0
                && line != "__EXIT_STATUS=0") {
                std::cerr << "nonzero exact-correction status in " << path << '\n';
                std::exit(2);
            }
            if (line.rfind("__PART_EXIT_STATUS=", 0) == 0
                && line != "__PART_EXIT_STATUS=0") {
                std::cerr << "nonzero exact-correction part status in " << path << '\n';
                std::exit(2);
            }
            if (line == std::string("__HOST=") + EXACT_MOMENT_SOURCE_HOST) {
                ++host_rows;
            } else if (line.rfind("__HOST=", 0) == 0) {
                std::cerr << "exact-moment source is not from machine C: "
                          << path << '\n';
                std::exit(2);
            }
            if (line == std::string("__SOURCE_SHA256=")
                    + EXACT_MOMENT_SOURCE_SHA256) {
                ++source_hash_rows;
            } else if (line.rfind("__SOURCE_SHA256=", 0) == 0) {
                std::cerr << "unexpected exact-moment source hash: "
                          << path << '\n';
                std::exit(2);
            }
            std::smatch match;
            int rank = 0;
            int index = 0;
            mpz_class value;
            if (std::regex_search(line, match, exact_weyl_pattern)) {
                rank = std::stoi(match[1].str());
                const int printed_moment = std::stoi(match[2].str());
                const mpz_class moment_value(match[3].str());
                index = std::stoi(match[4].str());
                value = mpz_class(match[5].str());
                if (printed_moment != index || moment_value < 0
                    || index < 0 || index >= static_cast<int>(stable.size())
                    || stable[static_cast<std::size_t>(index)] + value
                        != moment_value) {
                    std::cerr << "malformed exact Weyl moment in "
                              << path << '\n';
                    std::exit(2);
                }
                ++exact_weyl_rows;
            } else if (std::regex_search(line, match, full_pattern)) {
                rank = std::stoi(match[1].str());
                const int printed_moment = std::stoi(match[2].str());
                const mpz_class even_overflow(match[3].str());
                const mpz_class determinant(match[4].str());
                index = std::stoi(match[5].str());
                value = mpz_class(match[6].str());
                const int moment_index = std::stoi(match[7].str());
                const mpz_class moment_value(match[8].str());
                if (printed_moment != index
                    || moment_index != index
                    || index < 0 || index >= static_cast<int>(stable.size())
                    || even_overflow + determinant != value
                    || stable[static_cast<std::size_t>(index)] + value
                        != moment_value) {
                    std::cerr << "malformed full exact correction in "
                              << path << '\n';
                    std::exit(2);
                }
                found_full_row = true;
                ++full_rows;
            } else if (std::regex_search(
                           line, match, stable_window_pattern)) {
                rank = std::stoi(match[1].str());
                index = std::stoi(match[2].str());
                value = mpz_class(match[3].str());
                if (index > 2 * rank) {
                    // These legacy rows contain only the determinant
                    // component B_j.  It equals the full correction only
                    // while the O(2 rank) overflow A_j still vanishes.
                    ++ignored_legacy_rows;
                    continue;
                }
                ++stable_window_rows;
            } else {
                continue;
            }
            auto& rank_corrections = corrections[rank];
            const auto existing = rank_corrections.find(index);
            if (existing != rank_corrections.end() && existing->second != value) {
                std::cerr << "conflicting exact correction D_" << rank
                          << " Delta_" << index << '\n';
                std::exit(2);
            }
            rank_corrections[index] = value;
            found_usable_row = true;
        }
        if (success_rows != 1) {
            std::cerr << path << " does not record exactly one zero status\n";
            std::exit(2);
        }
        if (found_full_row && (host_rows != 1 || source_hash_rows != 1)) {
            std::cerr << path
                      << " lacks unique machine-C exact-source provenance\n";
            std::exit(2);
        }
        if (!found_usable_row) {
            std::cerr << path << " contains no usable exact-correction rows\n";
            std::exit(2);
        }
    }
    if (!paths.empty()) {
        std::cout << "exact_correction_input_logs=" << paths.size()
                  << " full_rows=" << full_rows
                  << " exact_weyl_rows=" << exact_weyl_rows
                  << " stable_window_legacy_rows=" << stable_window_rows
                  << " ignored_out_of_window_legacy_rows="
                  << ignored_legacy_rows << '\n';
    }
    return corrections;
}

std::vector<mpz_class> determinant_upper_bounds(
    int maximum,
    const std::vector<mpz_class>& stable
) {
    std::vector<mpz_class> upper(static_cast<std::size_t>(maximum + 1));
    for (int index = 0; index <= maximum; ++index) {
        // The O^-(2 rank) component identity gives
        // 0 <= B_j <= p_j=s_j-A_j <= s_j for every parity.
        upper[static_cast<std::size_t>(index)] =
            stable[static_cast<std::size_t>(index)];
    }
    return upper;
}

Interval correction_interval(
    int index,
    int rank,
    const std::vector<mpz_class>& stable,
    const std::vector<mpz_class>& determinant_upper,
    const std::map<int, mpz_class>& exact_corrections
) {
    Interval proved;
    if (index < rank) {
        proved = {0, 0};
    } else if (index <= 2 * rank) {
        // The omitted even-row sum A_j(rank) vanishes in the stable O(2 rank)
        // window.  Hence delta_j=B_j>=0.  The determinant-isotypic integral
        // obeys B_j<=p_j<=s_j for every parity.
        proved = {0, determinant_upper[static_cast<std::size_t>(index)]};
    } else {
        proved = {
            -stable[static_cast<std::size_t>(index)],
            determinant_upper[static_cast<std::size_t>(index)]
        };
    }

    const auto exact = exact_corrections.find(index);
    if (exact == exact_corrections.end()) return proved;
    if (exact->second < proved.lower || exact->second > proved.upper) {
        std::cerr << "exact correction lies outside proved box: D_" << rank
                  << " Delta_" << index << '\n';
        std::exit(2);
    }
    return {exact->second, exact->second};
}

mpz_class linear_lower(const mpz_class& coefficient, const Interval& interval) {
    return coefficient >= 0
        ? coefficient * interval.lower
        : coefficient * interval.upper;
}

mpz_class product_lower(
    const mpz_class& coefficient,
    const Interval& first,
    const Interval& second,
    bool diagonal
) {
    mpz_class product_min;
    mpz_class product_max;
    if (diagonal) {
        const mpz_class lower_square = first.lower * first.lower;
        const mpz_class upper_square = first.upper * first.upper;
        product_max = std::max(lower_square, upper_square);
        if (first.lower <= 0 && first.upper >= 0) {
            product_min = 0;
        } else {
            product_min = std::min(lower_square, upper_square);
        }
    } else {
        const std::array<mpz_class, 4> products{
            first.lower * second.lower,
            first.lower * second.upper,
            first.upper * second.lower,
            first.upper * second.upper,
        };
        product_min = *std::min_element(products.begin(), products.end());
        product_max = *std::max_element(products.begin(), products.end());
    }
    return coefficient >= 0
        ? coefficient * product_min
        : coefficient * product_max;
}

void accumulate_diagonal(
    int n,
    int scale,
    int rank,
    const std::vector<mpz_class>& stable,
    const std::vector<mpz_class>& determinant_upper,
    const std::map<int, mpz_class>& exact_corrections,
    std::vector<mpz_class>& linear,
    mpz_class& stable_value,
    mpz_class& quadratic_interval_lower
) {
    const std::vector<mpz_class> row = binomial_row(n);
    const int pair_sum = n + 2;
    for (int first = 0; first <= pair_sum / 2; ++first) {
        const int second = pair_sum - first;
        const mpz_class coefficient = pair_coefficient(row, first, second, scale);
        if (coefficient == 0) continue;

        stable_value += coefficient
            * stable[static_cast<std::size_t>(first)]
            * stable[static_cast<std::size_t>(second)];

        const Interval first_interval = correction_interval(
            first, rank, stable, determinant_upper, exact_corrections);
        const Interval second_interval = correction_interval(
            second, rank, stable, determinant_upper, exact_corrections);
        quadratic_interval_lower += product_lower(
            coefficient,
            first_interval,
            second_interval,
            first == second
        );

        if (first == second) {
            linear[static_cast<std::size_t>(first)] +=
                2 * coefficient * stable[static_cast<std::size_t>(first)];
        } else {
            linear[static_cast<std::size_t>(first)] +=
                coefficient * stable[static_cast<std::size_t>(second)];
            linear[static_cast<std::size_t>(second)] +=
                coefficient * stable[static_cast<std::size_t>(first)];
        }
    }
}

struct AuditResult {
    int rank = 0;
    int chain_m = 0;
    mpz_class stable_value;
    mpz_class linear_interval_lower;
    mpz_class quadratic_interval_lower;
    mpz_class lower;
};

AuditResult audit_step(
    int rank,
    int chain_m,
    const std::vector<mpz_class>& stable,
    const std::vector<mpz_class>& determinant_upper,
    const std::map<int, mpz_class>& exact_corrections
) {
    const int maximum = 2 * chain_m + 5;
    if (maximum < 0
        || maximum >= static_cast<int>(stable.size())
        || maximum >= static_cast<int>(determinant_upper.size())) {
        std::cerr << "precomputed moment table is too short\n";
        std::exit(2);
    }

    std::vector<mpz_class> linear(stable.size());
    AuditResult result;
    result.rank = rank;
    result.chain_m = chain_m;

    accumulate_diagonal(
        2 * chain_m + 3,
        1,
        rank,
        stable,
        determinant_upper,
        exact_corrections,
        linear,
        result.stable_value,
        result.quadratic_interval_lower
    );
    accumulate_diagonal(
        2 * chain_m + 1,
        -4,
        rank,
        stable,
        determinant_upper,
        exact_corrections,
        linear,
        result.stable_value,
        result.quadratic_interval_lower
    );

    for (int index = rank; index <= maximum; ++index) {
        result.linear_interval_lower += linear_lower(
            linear[static_cast<std::size_t>(index)],
            correction_interval(
                index, rank, stable, determinant_upper, exact_corrections)
        );
    }
    result.lower = result.stable_value
        + result.linear_interval_lower
        + result.quadratic_interval_lower;
    return result;
}

struct CertificateRank {
    int rank;
    int onset;
};

struct ExactPrefixRank {
    int rank;
    int exact_through;
};

constexpr std::array<ExactPrefixRank, 21> D4_24_EXACT_PREFIXES{{
    {4, 59}, {5, 44}, {6, 40}, {7, 40}, {8, 38}, {9, 38},
    {10, 38}, {11, 38}, {12, 38}, {13, 38}, {14, 38}, {15, 38},
    {16, 38}, {17, 38}, {18, 38}, {19, 34}, {20, 30}, {21, 30},
    {22, 30}, {23, 30}, {24, 30},
}};

constexpr std::array<CertificateRank, 28> D25_52_ONSETS{{
    {25, 71}, {26, 73}, {27, 71}, {28, 73},
    {29, 73}, {30, 73}, {31, 75}, {32, 77},
    {33, 79}, {34, 79}, {35, 81}, {36, 83},
    {37, 83}, {38, 85}, {39, 85}, {40, 87},
    {41, 89}, {42, 89}, {43, 91}, {44, 93},
    {45, 93}, {46, 95}, {47, 97}, {48, 97},
    {49, 99}, {50, 99}, {51, 101}, {52, 103},
}};

}  // namespace

int main(int argc, char** argv) {
    int rank = 98;
    int chain_m = 102;
    int onset = 0;
    bool certificate_d4_24 = false;
    bool certificate_d25_52 = false;
    bool certificate_d53_63_row_gated = false;
    std::vector<std::string> exact_correction_paths;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--certificate-d4-24") {
            certificate_d4_24 = true;
        } else if (argument == "--certificate-d25-52") {
            certificate_d25_52 = true;
        } else if (argument == "--certificate-d53-63-row-gated") {
            certificate_d53_63_row_gated = true;
        } else if (argument == "--rank" && index + 1 < argc) {
            rank = std::atoi(argv[++index]);
        } else if (argument == "--chain-m" && index + 1 < argc) {
            chain_m = std::atoi(argv[++index]);
        } else if (argument == "--onset" && index + 1 < argc) {
            onset = std::atoi(argv[++index]);
        } else if (argument == "--exact-correction-log" && index + 1 < argc) {
            exact_correction_paths.emplace_back(argv[++index]);
        } else {
            std::cerr << "usage: " << argv[0]
                      << " [--rank B] [--chain-m M]"
                      << " [--onset ODD_N]"
                      << " [--exact-correction-log PATH]"
                      << " [--certificate-d4-24]"
                      << " [--certificate-d25-52]"
                      << " [--certificate-d53-63-row-gated]\n";
            return 2;
        }
    }
    if (rank < 1 || chain_m < 0 || onset < 0) {
        std::cerr << "rank and chain_m must be nonnegative\n";
        return 2;
    }
    if (onset != 0 && (onset % 2 == 0 || onset < rank - 1)) {
        std::cerr << "onset must be an odd integer at least rank-1\n";
        return 2;
    }
    if (static_cast<int>(certificate_d4_24)
            + static_cast<int>(certificate_d25_52)
            + static_cast<int>(certificate_d53_63_row_gated) > 1) {
        std::cerr << "select at most one fixed certificate mode\n";
        return 2;
    }
    const std::vector<mpz_class> source_stable =
        exact_correction_paths.empty()
            ? std::vector<mpz_class>{}
            : stable_moments(600);
    const ExactCorrections exact_corrections =
        read_exact_correction_logs(exact_correction_paths, source_stable);
    const std::map<int, mpz_class> no_exact_corrections;

    if (certificate_d4_24) {
        if (exact_correction_paths.empty()) {
            std::cerr << "D4-D24 certificate requires exact-correction logs\n";
            return 2;
        }
        constexpr int LAST_CHAIN_M = 29;
        constexpr int MAXIMUM_MOMENT = 2 * LAST_CHAIN_M + 5;
        const std::vector<mpz_class> stable =
            stable_moments(MAXIMUM_MOMENT + 1);
        const std::vector<mpz_class> determinant_upper =
            determinant_upper_bounds(MAXIMUM_MOMENT, stable);

        std::vector<std::pair<int, int>> tasks;
        for (const ExactPrefixRank& item : D4_24_EXACT_PREFIXES) {
            const auto supplied = exact_corrections.find(item.rank);
            if (supplied == exact_corrections.end()) {
                std::cerr << "missing exact corrections for D_"
                          << item.rank << '\n';
                return 2;
            }
            for (int index = item.rank; index <= item.exact_through; ++index) {
                if (supplied->second.find(index) == supplied->second.end()) {
                    std::cerr << "missing exact correction D_" << item.rank
                              << " Delta_" << index << '\n';
                    return 2;
                }
                (void)correction_interval(
                    index,
                    item.rank,
                    stable,
                    determinant_upper,
                    supplied->second
                );
            }
            const int first_chain_m = (item.rank - 4) / 2;
            for (int m = first_chain_m; m <= LAST_CHAIN_M; ++m) {
                tasks.emplace_back(item.rank, m);
            }
        }

        std::vector<AuditResult> results(tasks.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (long long index = 0;
             index < static_cast<long long>(tasks.size());
             ++index) {
            const int task_rank = tasks[static_cast<std::size_t>(index)].first;
            results[static_cast<std::size_t>(index)] = audit_step(
                task_rank,
                tasks[static_cast<std::size_t>(index)].second,
                stable,
                determinant_upper,
                exact_corrections.at(task_rank)
            );
        }

        int failures = 0;
        bool have_minimum = false;
        mpz_class minimum_lower;
        int minimum_rank = 0;
        int minimum_chain_m = 0;
        for (const AuditResult& result : results) {
            if (result.lower <= 0) ++failures;
            if (!have_minimum || result.lower < minimum_lower) {
                minimum_lower = result.lower;
                minimum_rank = result.rank;
                minimum_chain_m = result.chain_m;
                have_minimum = true;
            }
        }

        std::cout << "D4-D24 exact-prefix/two-sided-suffix GMP interval certificate\n"
                  << "arithmetic=exact_integer"
                  << " correction_window=j<rank:0,rank<=j<=2rank:[0,U_j],"
                  << "j>2rank:[-s_j,U_j]\n"
                  << "upper_bound=U_j:s_j_all_parities\n";
#ifdef _OPENMP
        std::cout << "openmp_max_threads=" << omp_get_max_threads() << '\n';
#else
        std::cout << "openmp_max_threads=1\n";
#endif

        std::size_t offset = 0;
        for (const ExactPrefixRank& item : D4_24_EXACT_PREFIXES) {
            const int first_chain_m = (item.rank - 4) / 2;
            bool have_rank_minimum = false;
            mpz_class rank_minimum;
            int rank_minimum_m = 0;
            for (int m = first_chain_m; m <= LAST_CHAIN_M; ++m, ++offset) {
                const AuditResult& result = results[offset];
                if (!have_rank_minimum || result.lower < rank_minimum) {
                    rank_minimum = result.lower;
                    rank_minimum_m = m;
                    have_rank_minimum = true;
                }
            }
            std::cout << "row=D_" << item.rank
                      << " exact_prefix=" << item.rank << ".."
                      << item.exact_through
                      << " first_chain_m=" << first_chain_m
                      << " last_chain_m=" << LAST_CHAIN_M
                      << " checks=" << (LAST_CHAIN_M - first_chain_m + 1)
                      << " minimum_lower_sign=" << sgn(rank_minimum)
                      << " minimum_lower=" << rank_minimum
                      << " minimum_lower_log10=" << std::setprecision(18)
                      << log10_abs(rank_minimum)
                      << " minimum_chain_m=" << rank_minimum_m << '\n';
        }
        std::cout << "MINIMUM lower_sign=" << sgn(minimum_lower)
                  << " lower=" << minimum_lower
                  << " lower_log10=" << std::setprecision(18)
                  << log10_abs(minimum_lower)
                  << " rank=D_" << minimum_rank
                  << " chain_m=" << minimum_chain_m << '\n'
                  << "SUMMARY ranks_checked=" << D4_24_EXACT_PREFIXES.size()
                  << " chain_steps_checked=" << results.size()
                  << " failures=" << failures << '\n'
                  << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << '\n';
        return failures == 0 ? 0 : 1;
    }

    if (certificate_d25_52) {
        constexpr int MAXIMUM_MOMENT = 105;
        const std::vector<mpz_class> stable =
            stable_moments(MAXIMUM_MOMENT + 1);
        const std::vector<mpz_class> determinant_upper =
            determinant_upper_bounds(MAXIMUM_MOMENT, stable);

        std::vector<std::pair<int, int>> tasks;
        for (const CertificateRank& item : D25_52_ONSETS) {
            if (item.onset % 2 == 0 || item.onset < 1) {
                std::cerr << "invalid onset for D_" << item.rank << '\n';
                return 2;
            }
            const int first_chain_m = (item.rank - 4) / 2;
            const int last_chain_m = (item.onset - 3) / 2;
            for (int m = first_chain_m; m <= last_chain_m; ++m) {
                tasks.emplace_back(item.rank, m);
            }
        }

        std::vector<AuditResult> results(tasks.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (long long index = 0;
             index < static_cast<long long>(tasks.size());
             ++index) {
            results[static_cast<std::size_t>(index)] = audit_step(
                tasks[static_cast<std::size_t>(index)].first,
                tasks[static_cast<std::size_t>(index)].second,
                stable,
                determinant_upper,
                no_exact_corrections
            );
        }

        int failures = 0;
        bool have_minimum = false;
        mpz_class minimum_lower;
        int minimum_rank = 0;
        int minimum_chain_m = 0;
        for (const AuditResult& result : results) {
            if (result.lower <= 0) ++failures;
            if (!have_minimum || result.lower < minimum_lower) {
                minimum_lower = result.lower;
                minimum_rank = result.rank;
                minimum_chain_m = result.chain_m;
                have_minimum = true;
            }
        }

        std::cout << "D25-D52 correction-aware GMP interval certificate\n"
                  << "arithmetic=exact_integer"
                  << " correction_window=j<rank:0,rank<=j<=2rank:[0,U_j],"
                  << "j>2rank:[-s_j,U_j]\n"
                  << "upper_bound=U_j:s_j_all_parities\n";
#ifdef _OPENMP
        std::cout << "openmp_max_threads=" << omp_get_max_threads() << '\n';
#else
        std::cout << "openmp_max_threads=1\n";
#endif

        std::size_t offset = 0;
        for (const CertificateRank& item : D25_52_ONSETS) {
            const int first_chain_m = (item.rank - 4) / 2;
            const int last_chain_m = (item.onset - 3) / 2;
            bool have_rank_minimum = false;
            mpz_class rank_minimum;
            int rank_minimum_m = 0;
            for (int m = first_chain_m; m <= last_chain_m; ++m, ++offset) {
                const AuditResult& result = results[offset];
                if (!have_rank_minimum || result.lower < rank_minimum) {
                    rank_minimum = result.lower;
                    rank_minimum_m = m;
                    have_rank_minimum = true;
                }
            }
            std::cout << "row=D_" << item.rank
                      << " onset=" << item.onset
                      << " first_chain_m=" << first_chain_m
                      << " last_chain_m=" << last_chain_m
                      << " checks=" << (last_chain_m - first_chain_m + 1)
                      << " minimum_lower_sign=" << sgn(rank_minimum)
                      << " minimum_lower=" << rank_minimum
                      << " minimum_lower_log10=" << std::setprecision(18)
                      << log10_abs(rank_minimum)
                      << " minimum_chain_m=" << rank_minimum_m << '\n';
        }
        std::cout << "MINIMUM lower_sign=" << sgn(minimum_lower)
                  << " lower=" << minimum_lower
                  << " lower_log10=" << std::setprecision(18)
                  << log10_abs(minimum_lower)
                  << " rank=D_" << minimum_rank
                  << " chain_m=" << minimum_chain_m << '\n'
                  << "SUMMARY ranks_checked=" << D25_52_ONSETS.size()
                  << " chain_steps_checked=" << results.size()
                  << " failures=" << failures << '\n'
                  << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << '\n';
        return failures == 0 ? 0 : 1;
    }

    if (certificate_d53_63_row_gated) {
        constexpr int FIRST_RANK = 53;
        constexpr int LAST_RANK = 63;
        constexpr int LAST_CHAIN_M = 29;
        constexpr int MAXIMUM_MOMENT = 2 * LAST_CHAIN_M + 5;
        const std::vector<mpz_class> stable =
            stable_moments(MAXIMUM_MOMENT + 1);
        const std::vector<mpz_class> determinant_upper =
            determinant_upper_bounds(MAXIMUM_MOMENT, stable);

        std::vector<std::pair<int, int>> tasks;
        for (int task_rank = FIRST_RANK; task_rank <= LAST_RANK; ++task_rank) {
            const int first_chain_m = (task_rank - 4) / 2;
            for (int m = first_chain_m; m <= LAST_CHAIN_M; ++m) {
                tasks.emplace_back(task_rank, m);
            }
        }
        std::vector<AuditResult> results(tasks.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (long long index = 0;
             index < static_cast<long long>(tasks.size());
             ++index) {
            results[static_cast<std::size_t>(index)] = audit_step(
                tasks[static_cast<std::size_t>(index)].first,
                tasks[static_cast<std::size_t>(index)].second,
                stable,
                determinant_upper,
                no_exact_corrections
            );
        }

        int failures = 0;
        bool have_minimum = false;
        mpz_class minimum_lower;
        int minimum_rank = 0;
        int minimum_chain_m = 0;
        for (const AuditResult& result : results) {
            if (result.lower <= 0) ++failures;
            if (!have_minimum || result.lower < minimum_lower) {
                have_minimum = true;
                minimum_lower = result.lower;
                minimum_rank = result.rank;
                minimum_chain_m = result.chain_m;
            }
        }

        std::cout << "D53-D63 source-generated row-gated GMP interval certificate\n"
                  << "arithmetic=exact_integer"
                  << " correction_window=j<rank:0,rank<=j<=2rank:[0,s_j]\n";
#ifdef _OPENMP
        std::cout << "openmp_max_threads=" << omp_get_max_threads() << '\n';
#else
        std::cout << "openmp_max_threads=1\n";
#endif
        std::size_t offset = 0;
        for (int task_rank = FIRST_RANK; task_rank <= LAST_RANK; ++task_rank) {
            const int first_chain_m = (task_rank - 4) / 2;
            bool have_rank_minimum = false;
            mpz_class rank_minimum;
            int rank_minimum_m = 0;
            for (int m = first_chain_m; m <= LAST_CHAIN_M; ++m, ++offset) {
                const AuditResult& result = results[offset];
                if (!have_rank_minimum || result.lower < rank_minimum) {
                    have_rank_minimum = true;
                    rank_minimum = result.lower;
                    rank_minimum_m = m;
                }
            }
            std::cout << "row=D_" << task_rank
                      << " first_chain_m=" << first_chain_m
                      << " last_chain_m=" << LAST_CHAIN_M
                      << " checks=" << (LAST_CHAIN_M - first_chain_m + 1)
                      << " minimum_lower_sign=" << sgn(rank_minimum)
                      << " minimum_lower=" << rank_minimum
                      << " minimum_chain_m=" << rank_minimum_m << '\n';
        }
        std::cout << "MINIMUM lower_sign=" << sgn(minimum_lower)
                  << " lower=" << minimum_lower
                  << " rank=D_" << minimum_rank
                  << " chain_m=" << minimum_chain_m << '\n'
                  << "SUMMARY ranks_checked=" << (LAST_RANK - FIRST_RANK + 1)
                  << " chain_steps_checked=" << results.size()
                  << " failures=" << failures << '\n'
                  << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << '\n';
        return failures == 0 ? 0 : 1;
    }

    if (onset != 0) {
        const int first_chain_m = (rank - 4) / 2;
        const int last_chain_m = (onset - 3) / 2;
        if (last_chain_m < first_chain_m) {
            std::cerr << "onset leaves an empty chain range\n";
            return 2;
        }
        const int maximum = 2 * last_chain_m + 5;
        const std::vector<mpz_class> stable = stable_moments(maximum + 1);
        const std::vector<mpz_class> determinant_upper =
            determinant_upper_bounds(maximum, stable);
        const auto supplied = exact_corrections.find(rank);
        const std::map<int, mpz_class>& rank_exact_corrections =
            supplied == exact_corrections.end()
                ? no_exact_corrections
                : supplied->second;

        const int checks = last_chain_m - first_chain_m + 1;
        std::vector<AuditResult> results(static_cast<std::size_t>(checks));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int offset = 0; offset < checks; ++offset) {
            results[static_cast<std::size_t>(offset)] = audit_step(
                rank,
                first_chain_m + offset,
                stable,
                determinant_upper,
                rank_exact_corrections
            );
        }

        int failures = 0;
        mpz_class minimum_lower = results.front().lower;
        int minimum_chain_m = first_chain_m;
        for (const AuditResult& result : results) {
            if (result.lower <= 0) ++failures;
            if (result.lower < minimum_lower) {
                minimum_lower = result.lower;
                minimum_chain_m = result.chain_m;
            }
        }
        std::cout << "D correction-aware GMP interval chain certificate\n"
                  << "arithmetic=exact_integer"
                  << " correction_window=j<rank:0,rank<=j<=2rank:[0,U_j],"
                  << "j>2rank:[-s_j,U_j]\n"
                  << "rank=D_" << rank
                  << " onset=" << onset
                  << " first_chain_m=" << first_chain_m
                  << " last_chain_m=" << last_chain_m
                  << " checks=" << checks
                  << " exact_corrections_available="
                  << rank_exact_corrections.size() << '\n'
                  << "MINIMUM lower_sign=" << sgn(minimum_lower)
                  << " lower=" << minimum_lower
                  << " lower_log10=" << std::setprecision(18)
                  << log10_abs(minimum_lower)
                  << " rank=D_" << rank
                  << " chain_m=" << minimum_chain_m << '\n'
                  << "SUMMARY ranks_checked=1 chain_steps_checked=" << checks
                  << " failures=" << failures << '\n'
                  << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << '\n';
        return failures == 0 ? 0 : 1;
    }

    const int maximum = 2 * chain_m + 5;
    const std::vector<mpz_class> stable = stable_moments(maximum + 1);
    const std::vector<mpz_class> determinant_upper =
        determinant_upper_bounds(maximum, stable);
    const auto supplied = exact_corrections.find(rank);
    const std::map<int, mpz_class>& rank_exact_corrections =
        supplied == exact_corrections.end()
            ? no_exact_corrections
            : supplied->second;
    const AuditResult result = audit_step(
        rank,
        chain_m,
        stable,
        determinant_upper,
        rank_exact_corrections
    );
    std::cout << "D stable-window/two-sided suffix GMP interval audit\n"
              << "rank=D_" << rank
              << " chain_m=" << chain_m
              << " odd_n=" << (2 * chain_m + 3)
              << " maximum_moment=" << maximum
              << " nonnegative_through=" << (2 * rank)
              << " exact_corrections_available="
              << rank_exact_corrections.size() << '\n'
              << std::setprecision(18)
              << "stable_sign=" << sgn(result.stable_value)
              << " stable_log10=" << log10_abs(result.stable_value) << '\n'
              << "linear_interval_sign=" << sgn(result.linear_interval_lower)
              << " linear_interval_log10=" << log10_abs(result.linear_interval_lower) << '\n'
              << "quadratic_interval_sign=" << sgn(result.quadratic_interval_lower)
              << " quadratic_interval_log10=" << log10_abs(result.quadratic_interval_lower) << '\n'
              << "lower_sign=" << sgn(result.lower)
              << " lower_log10=" << log10_abs(result.lower) << '\n'
              << "UPPER_BOUND_STATUS: uniform determinant bound U_j=s_j\n"
              << "SUMMARY ranks_checked=1 chain_steps_checked=1 failures="
              << (result.lower > 0 ? 0 : 1) << '\n'
              << "__EXIT_STATUS=" << (result.lower > 0 ? 0 : 1) << '\n';
    return result.lower > 0 ? 0 : 1;
}
