// Hybrid exact supplier for every active type-B H8--H27 frontier.
//
// For B_b with q=b+1, the bounded Littlewood identity gives
//
//   bad(q,j) = stable(j) - moment(B_b,j).
//
// Standardization injects every bad semistandard tableau into a standard
// tableau of the same shape.  The hook-length formula alone closes 296 cases.
// The remaining 41 low-rank cases are reconstructed by the bounded-Littlewood
// determinant/CRT engine.  The reverse traversals remain independent controls.

#include <array>
#include <cstddef>

// Supply the determinant engine with the exact B-only ledger required here.
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

inline constexpr std::array<RowCutoff, 6> row_cutoffs{{
    {'B', 14, 59, 57, TailMethod::directed_interval},
    {'B', 15, 61, 59, TailMethod::directed_interval},
    {'B', 16, 63, 61, TailMethod::directed_interval},
    {'B', 17, 65, 63, TailMethod::directed_interval},
    {'B', 18, 67, 65, TailMethod::directed_interval},
    {'B', 19, 69, 67, TailMethod::directed_interval},
}};

inline const RowCutoff* find_row(char family, int rank) {
    for (const RowCutoff& row : row_cutoffs) {
        if (row.family == family && row.rank == rank) return &row;
    }
    return nullptr;
}

inline constexpr int required_maximum_moment = 67;
inline constexpr std::size_t polynomial_rows = 0;
inline constexpr std::size_t rational_cap_rows = 0;
inline constexpr std::size_t directed_interval_rows = row_cutoffs.size();
inline constexpr std::size_t full_residual_pairs = 0;

}  // namespace full_q3_bcd_remaining

#define main full_q3_bounded_littlewood_main_disabled
#include "verify_full_q3_bcd_bounded_littlewood_gmp.cpp"
#undef main

#include <algorithm>
#include <cstdint>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace active_b_frontier {

using ResidueTable = std::vector<std::vector<std::uint32_t>>;

ResidueTable residues_for_prime_node(
    const std::vector<MomentRow>& rows,
    int maximum_moment,
    std::uint32_t prime,
    int node_index,
    const std::vector<std::vector<std::uint32_t>>& weights
) {
    const int degree_limit = 2 * maximum_moment;
    const MontgomeryField field(prime);
    ResidueTable residues(
        rows.size(),
        std::vector<std::uint32_t>(
            static_cast<std::size_t>(maximum_moment + 1), 0U
        )
    );
    const std::uint32_t inverse_two = field.inverse(field.from_unsigned(2U));

    int maximum_rank = 0;
    for (const MomentRow& row : rows) {
        if (row.family != 'B') {
            throw std::runtime_error("non-B row in active B frontier ledger");
        }
        maximum_rank = std::max(maximum_rank, row.rank);
    }

    const std::uint32_t node = field.from_unsigned(
        static_cast<std::uint64_t>(node_index)
    );
    const NodeData data = build_node_data(
        node, 2 * maximum_rank - 1, field, degree_limit
    );
    const auto minus_determinants = leading_principal_determinants(
        build_b_matrix(
            data.elementary_pair,
            maximum_rank,
            false,
            field,
            degree_limit
        ),
        maximum_rank,
        field,
        degree_limit
    );
    const auto plus_determinants = leading_principal_determinants(
        build_b_matrix(
            data.elementary_pair,
            maximum_rank,
            true,
            field,
            degree_limit
        ),
        maximum_rank,
        field,
        degree_limit
    );

    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        const MomentRow& row = rows[row_index];
        const Series first = multiply_series(
            data.elementary_sum,
            minus_determinants[static_cast<std::size_t>(row.rank)],
            field,
            degree_limit
        );
        const Series second = multiply_series(
            data.alternating_elementary_sum,
            plus_determinants[static_cast<std::size_t>(row.rank)],
            field,
            degree_limit
        );
        const Series generating = scale_series(
            add_series(first, second, field, degree_limit),
            inverse_two,
            field,
            degree_limit
        );
        for (int moment = 0; moment <= row.moment_through; ++moment) {
            const std::uint32_t coefficient = generating.coefficient[
                static_cast<std::size_t>(2 * moment)
            ];
            residues[row_index][static_cast<std::size_t>(moment)] =
                field.multiply(
                    weights[static_cast<std::size_t>(moment)]
                           [static_cast<std::size_t>(node_index)],
                    coefficient
                );
        }
    }
    return residues;
}

std::vector<MomentRow> reconstruct_moments(
    int maximum_moment,
    std::size_t& primes_consumed,
    int& modulus_bits
) {
    std::vector<MomentRow> rows;
    for (const auto& cutoff : full_q3_bcd_remaining::row_cutoffs) {
        rows.push_back(MomentRow{
            cutoff.family,
            cutoff.rank,
            cutoff.moment_through,
            std::vector<BigInt>(
                static_cast<std::size_t>(cutoff.moment_through + 1), 0
            )
        });
    }

    // A finite B moment is a subcount of its unrestricted stable moment, so
    // the stable sequence is a rigorous and much smaller CRT uniqueness bound
    // than dim(B_b)^j.
    const std::vector<BigInt> stable = stable_moments(maximum_moment);
    BigInt largest_bound = 0;
    for (const MomentRow& row : rows) {
        largest_bound = std::max(
            largest_bound,
            stable[static_cast<std::size_t>(row.moment_through)]
        );
    }

    std::vector<std::uint32_t> primes;
    BigInt planned_modulus = 1;
    for (std::uint32_t prime : descending_primes(32U)) {
        primes.push_back(prime);
        planned_modulus *= prime;
        if (planned_modulus > largest_bound) break;
    }
    if (planned_modulus <= largest_bound) {
        throw std::runtime_error("stable-bound CRT prime table is too short");
    }

    using WeightTable = std::vector<std::vector<std::uint32_t>>;
    std::vector<WeightTable> weights_by_prime(primes.size());
    const int prime_count = static_cast<int>(primes.size());
    const int node_count = maximum_moment + 1;

    // The former implementation assigned one OpenMP job per prime.  At this
    // frontier that exposed only 22 jobs and made every job traverse all 122
    // interpolation nodes serially.  The mathematical evaluations are
    // independent on the full prime-by-node grid, so precompute the small
    // Gaussian-functional tables once per prime and flatten that grid.  The
    // reduction below remains in increasing node order and the CRT merge
    // remains in the fixed descending-prime order.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int prime_index = 0; prime_index < prime_count; ++prime_index) {
        const std::size_t index = static_cast<std::size_t>(prime_index);
        const MontgomeryField field(primes[index]);
        weights_by_prime[index] = gaussian_functional_weights(
            maximum_moment, field
        );
    }

    const int task_count = prime_count * node_count;
    std::vector<ResidueTable> residues_by_task(
        static_cast<std::size_t>(task_count)
    );
    std::vector<std::string> errors(static_cast<std::size_t>(task_count));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int task_index = 0; task_index < task_count; ++task_index) {
        const int prime_index = task_index / node_count;
        const int node_index = task_index % node_count;
        const std::size_t task = static_cast<std::size_t>(task_index);
        const std::size_t prime = static_cast<std::size_t>(prime_index);
        try {
            residues_by_task[task] = residues_for_prime_node(
                rows,
                maximum_moment,
                primes[prime],
                node_index,
                weights_by_prime[prime]
            );
        } catch (const std::exception& error) {
            errors[task] = error.what();
        } catch (...) {
            errors[task] = "unknown prime-field node failure";
        }
    }
    for (std::size_t index = 0; index < errors.size(); ++index) {
        if (!errors[index].empty()) {
            const std::size_t prime_index = index
                / static_cast<std::size_t>(node_count);
            const std::size_t node_index = index
                % static_cast<std::size_t>(node_count);
            throw std::runtime_error(
                std::string("prime ") + std::to_string(primes[prime_index])
                + " node " + std::to_string(node_index)
                + " failed: " + errors[index]
            );
        }
    }

    std::vector<ResidueTable> residues_by_prime(primes.size());
    for (int prime_index = 0; prime_index < prime_count; ++prime_index) {
        const std::size_t prime = static_cast<std::size_t>(prime_index);
        const MontgomeryField field(primes[prime]);
        ResidueTable& total = residues_by_prime[prime];
        total.assign(
            rows.size(),
            std::vector<std::uint32_t>(
                static_cast<std::size_t>(maximum_moment + 1), 0U
            )
        );
        for (int node_index = 0; node_index < node_count; ++node_index) {
            const ResidueTable& contribution = residues_by_task[
                static_cast<std::size_t>(prime_index * node_count + node_index)
            ];
            for (std::size_t row_index = 0;
                 row_index < rows.size();
                 ++row_index) {
                for (int moment = 0;
                     moment <= rows[row_index].moment_through;
                     ++moment) {
                    std::uint32_t& value = total[row_index][
                        static_cast<std::size_t>(moment)
                    ];
                    value = field.add(
                        value,
                        contribution[row_index][
                            static_cast<std::size_t>(moment)
                        ]
                    );
                }
            }
        }
        for (auto& row : total) {
            for (std::uint32_t& value : row) {
                value = field.to_unsigned(value);
            }
        }
        std::cout << "ACTIVE_B_FRONTIER evaluated_prime=" << primes[prime]
                  << " prime_index=" << prime_index << '/' << prime_count
                  << std::endl;
    }

    BigInt modulus_product = 1;
    primes_consumed = 0U;
    for (std::size_t prime_index = 0;
         prime_index < primes.size();
         ++prime_index) {
        const std::uint32_t prime = primes[prime_index];
        const ResidueTable& residues = residues_by_prime[prime_index];
        const std::uint32_t modulus_residue = static_cast<std::uint32_t>(
            mpz_fdiv_ui(modulus_product.get_mpz_t(), prime)
        );
        const std::uint32_t inverse_modulus = normal_mod_inverse(
            modulus_residue, prime
        );
        for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
            MomentRow& row = rows[row_index];
            for (int moment = 0; moment <= row.moment_through; ++moment) {
                BigInt& value = row.moments[static_cast<std::size_t>(moment)];
                const std::uint32_t current = static_cast<std::uint32_t>(
                    mpz_fdiv_ui(value.get_mpz_t(), prime)
                );
                const std::uint32_t residue =
                    residues[row_index][static_cast<std::size_t>(moment)];
                const std::uint32_t difference = residue >= current
                    ? residue - current
                    : static_cast<std::uint32_t>(
                          static_cast<std::uint64_t>(residue) + prime - current
                      );
                const std::uint32_t multiplier = static_cast<std::uint32_t>(
                    (static_cast<std::uint64_t>(difference) * inverse_modulus)
                    % prime
                );
                value += modulus_product * multiplier;
            }
        }
        modulus_product *= prime;
        ++primes_consumed;
    }
    modulus_bits = static_cast<int>(
        mpz_sizeinbase(modulus_product.get_mpz_t(), 2)
    );
    if (modulus_product <= largest_bound) {
        throw std::runtime_error("reconstructed CRT modulus misses stable bound");
    }
    for (const MomentRow& row : rows) {
        for (int moment = 0; moment <= row.moment_through; ++moment) {
            const BigInt& value = row.moments[static_cast<std::size_t>(moment)];
            if (value < 0
                || value > stable[static_cast<std::size_t>(moment)]) {
                throw std::runtime_error(
                    "finite B moment lies outside the stable uniqueness bound"
                );
            }
        }
    }
    return rows;
}

struct Case {
    int h;
    int q;
    int j;
};

std::vector<Case> active_cases() {
    std::vector<Case> cases;
    for (int q = 15; q <= 17; ++q) cases.push_back(Case{8, q, 2 * q + 8});
    for (int h = 9; h <= 23; ++h) {
        const int q_high = (3 * h + 1) / 2 + 4;
        for (int q = 15; q <= q_high; ++q) {
            cases.push_back(Case{h, q, 2 * q + h});
        }
    }
    for (int h = 24; h <= 27; ++h) {
        for (int q = 15; q <= 2 * h - 7; ++q) {
            cases.push_back(Case{h, q, 2 * q + h});
        }
    }
    if (cases.size() != 337U) {
        throw std::runtime_error("active B frontier scope is not 337 cases");
    }
    return cases;
}

BigInt standard_tableau_dimension(
    const std::vector<int>& shape,
    const std::vector<BigInt>& factorials
) {
    const int row_count = static_cast<int>(shape.size());
    int size = 0;
    for (int part : shape) size += part;
    BigInt numerator = factorials[static_cast<std::size_t>(size)];
    BigInt denominator = 1;
    for (int row = 0; row < row_count; ++row) {
        denominator *= factorials[static_cast<std::size_t>(
            shape[static_cast<std::size_t>(row)] + row_count - row - 1
        )];
        for (int lower = row + 1; lower < row_count; ++lower) {
            numerator *=
                shape[static_cast<std::size_t>(row)]
                - shape[static_cast<std::size_t>(lower)]
                + lower - row;
        }
    }
    if (denominator == 0 || numerator % denominator != 0) {
        throw std::runtime_error("nonintegral hook-length dimension");
    }
    return numerator / denominator;
}

void add_partition_dimensions(
    int remaining,
    int maximum_part,
    int first_part,
    std::vector<int>& suffix,
    const std::vector<BigInt>& factorials,
    BigInt& total,
    std::size_t& shapes
) {
    if (remaining == 0) {
        std::vector<int> shape;
        shape.reserve(2U * (suffix.size() + 1U));
        shape.push_back(first_part);
        shape.push_back(first_part);
        for (int part : suffix) {
            shape.push_back(part);
            shape.push_back(part);
        }
        total += standard_tableau_dimension(shape, factorials);
        ++shapes;
        return;
    }
    for (int part = std::min(remaining, maximum_part); part >= 1; --part) {
        suffix.push_back(part);
        add_partition_dimensions(
            remaining - part,
            part,
            first_part,
            suffix,
            factorials,
            total,
            shapes
        );
        suffix.pop_back();
    }
}

// A bad B frontier tableau has column-even shape
//   lambda=(mu_1,mu_1,mu_2,mu_2,...),
// with |mu|=j and mu_1>=2q.  Standardizing its two copies of each label gives
// an injective map into the standard tableaux of shape lambda.  Summing the
// hook-length dimensions therefore gives a rigorous upper bound for bad(q,j).
BigInt hook_length_bad_upper_bound(
    const Case& frontier,
    const std::vector<BigInt>& factorials,
    std::size_t& shapes
) {
    const int excess = frontier.j - 2 * frontier.q;
    if (excess < 0) throw std::runtime_error("negative frontier excess");
    BigInt total = 0;
    shapes = 0U;
    std::vector<int> suffix;
    for (int first_excess = 0;
         first_excess <= excess;
         ++first_excess) {
        const int remainder = excess - first_excess;
        add_partition_dimensions(
            remainder,
            remainder,
            2 * frontier.q + first_excess,
            suffix,
            factorials,
            total,
            shapes
        );
    }
    return total;
}

}  // namespace active_b_frontier

int main() {
    try {
        std::cout << std::unitbuf;
        std::size_t primes_consumed = 0U;
        int modulus_bits = 0;
        const std::vector<MomentRow> rows =
            active_b_frontier::reconstruct_moments(
                full_q3_bcd_remaining::required_maximum_moment,
                primes_consumed,
                modulus_bits
            );
        constexpr int active_maximum_moment = 121;
        const std::vector<BigInt> stable = stable_moments(
            active_maximum_moment
        );
        std::vector<BigInt> factorials(
            static_cast<std::size_t>(2 * active_maximum_moment + 1), 1
        );
        for (int value = 1; value <= 2 * active_maximum_moment; ++value) {
            factorials[static_cast<std::size_t>(value)] =
                value * factorials[static_cast<std::size_t>(value - 1)];
        }
        std::map<int, const MomentRow*> by_rank;
        for (const MomentRow& row : rows) {
            if (!by_rank.emplace(row.rank, &row).second) {
                throw std::runtime_error("duplicate rank in active B ledger");
            }
        }

        int failures = 0;
        int cases_checked = 0;
        int analytic_cases = 0;
        int exact_cases = 0;
        std::size_t hook_shapes = 0U;
        BigInt minimum_margin;
        bool have_minimum = false;
        std::map<int, std::pair<int, BigInt>> summaries_by_h;
        for (const active_b_frontier::Case& frontier :
             active_b_frontier::active_cases()) {
            const BigInt& stable_value = stable[
                static_cast<std::size_t>(frontier.j)
            ];
            std::size_t case_shapes = 0U;
            const BigInt hook_bound =
                active_b_frontier::hook_length_bad_upper_bound(
                    frontier, factorials, case_shapes
                );
            hook_shapes += case_shapes;
            BigInt bad_bound = hook_bound;
            BigInt margin = stable_value - 2 * hook_bound;
            const char* method = "hook";
            if (margin > 0) {
                ++analytic_cases;
            } else {
                const int rank = frontier.q - 1;
                const auto found = by_rank.find(rank);
                if (found == by_rank.end()
                    || frontier.j > found->second->moment_through) {
                    throw std::runtime_error(
                        "exact residual frontier row is absent from ledger"
                    );
                }
                const BigInt finite = found->second->moments[
                    static_cast<std::size_t>(frontier.j)
                ];
                bad_bound = stable_value - finite;
                margin = stable_value - 2 * bad_bound;
                method = "determinant";
                ++exact_cases;
            }
            const bool pass = bad_bound >= 0 && margin > 0;
            if (!pass) ++failures;
            if (!have_minimum || margin < minimum_margin) {
                minimum_margin = margin;
                have_minimum = true;
            }
            auto& h_summary = summaries_by_h[frontier.h];
            if (h_summary.first == 0 || margin < h_summary.second) {
                h_summary.second = margin;
            }
            ++h_summary.first;
            ++cases_checked;
            if (!pass) {
                std::cout << "ACTIVE_B_FRONTIER_FAILURE h=" << frontier.h
                          << " q=" << frontier.q
                          << " rank=" << frontier.q - 1
                          << " j=" << frontier.j
                          << " method=" << method
                          << " stable=" << stable_value
                          << " bad_bound=" << bad_bound
                          << " margin=" << margin << '\n';
            }
        }
        if (cases_checked != 337 || analytic_cases != 296 || exact_cases != 41) {
            ++failures;
        }
        for (const auto& [h, summary] : summaries_by_h) {
            std::cout << "ACTIVE_B_FRONTIER_H_SUMMARY h=" << h
                      << " cases=" << summary.first
                      << " minimum_margin=" << summary.second << '\n';
        }
        std::cout << "ACTIVE_B_FRONTIER_HYBRID rows=" << rows.size()
                  << " cases=" << cases_checked
                  << " analytic_cases=" << analytic_cases
                  << " exact_cases=" << exact_cases
                  << " active_maximum_moment=" << active_maximum_moment
                  << " determinant_maximum_moment="
                  << full_q3_bcd_remaining::required_maximum_moment
                  << " hook_shapes=" << hook_shapes
                  << " primes=" << primes_consumed
                  << " modulus_bits=" << modulus_bits
                  << " minimum_margin=" << minimum_margin
                  << " failures=" << failures << '\n';
        if (failures == 0) {
            std::cout << "ACTIVE_B_FRONTIER_HYBRID "
                         "VERIFICATION: ALL PASS\n";
        }
        std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << '\n';
        return failures == 0 ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "ACTIVE_B_FRONTIER_BOUNDED_LITTLEWOOD FAILURE: "
                  << error.what() << '\n';
        std::cout << "__EXIT_STATUS=1\n";
        return 1;
    }
}
