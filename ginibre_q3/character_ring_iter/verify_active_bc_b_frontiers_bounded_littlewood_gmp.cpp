// Exact determinant/CRT supplier for every active type-B H8--H27 frontier.
//
// For B_b with q=b+1, the bounded Littlewood identity gives
//
//   bad(q,j) = stable(j) - moment(B_b,j).
//
// Thus a single finite-moment reconstruction certifies the half-stable
// inequality for all 337 cases which previously required separate reverse-
// Pieri traversals.  The reverse traversals remain available as independent
// controls, but are not inputs to this verifier.

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

inline constexpr std::array<RowCutoff, 33> row_cutoffs{{
    {'B', 14, 59, 58, TailMethod::directed_interval},
    {'B', 15, 61, 60, TailMethod::directed_interval},
    {'B', 16, 63, 62, TailMethod::directed_interval},
    {'B', 17, 65, 64, TailMethod::directed_interval},
    {'B', 18, 67, 66, TailMethod::directed_interval},
    {'B', 19, 69, 68, TailMethod::directed_interval},
    {'B', 20, 71, 70, TailMethod::directed_interval},
    {'B', 21, 73, 72, TailMethod::directed_interval},
    {'B', 22, 75, 74, TailMethod::directed_interval},
    {'B', 23, 77, 76, TailMethod::directed_interval},
    {'B', 24, 79, 78, TailMethod::directed_interval},
    {'B', 25, 81, 80, TailMethod::directed_interval},
    {'B', 26, 83, 82, TailMethod::directed_interval},
    {'B', 27, 85, 84, TailMethod::directed_interval},
    {'B', 28, 87, 86, TailMethod::directed_interval},
    {'B', 29, 89, 88, TailMethod::directed_interval},
    {'B', 30, 91, 90, TailMethod::directed_interval},
    {'B', 31, 93, 92, TailMethod::directed_interval},
    {'B', 32, 95, 94, TailMethod::directed_interval},
    {'B', 33, 97, 96, TailMethod::directed_interval},
    {'B', 34, 99, 98, TailMethod::directed_interval},
    {'B', 35, 101, 100, TailMethod::directed_interval},
    {'B', 36, 103, 102, TailMethod::directed_interval},
    {'B', 37, 105, 104, TailMethod::directed_interval},
    {'B', 38, 107, 106, TailMethod::directed_interval},
    {'B', 39, 109, 108, TailMethod::directed_interval},
    {'B', 40, 111, 110, TailMethod::directed_interval},
    {'B', 41, 113, 111, TailMethod::directed_interval},
    {'B', 42, 115, 113, TailMethod::directed_interval},
    {'B', 43, 117, 115, TailMethod::directed_interval},
    {'B', 44, 119, 117, TailMethod::directed_interval},
    {'B', 45, 121, 119, TailMethod::directed_interval},
    {'B', 46, 123, 121, TailMethod::directed_interval},
}};

inline const RowCutoff* find_row(char family, int rank) {
    for (const RowCutoff& row : row_cutoffs) {
        if (row.family == family && row.rank == rank) return &row;
    }
    return nullptr;
}

inline constexpr int required_maximum_moment = 121;
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

ResidueTable residues_for_prime(
    const std::vector<MomentRow>& rows,
    int maximum_moment,
    std::uint32_t prime
) {
    const int degree_limit = 2 * maximum_moment;
    const MontgomeryField field(prime);
    const auto weights = gaussian_functional_weights(maximum_moment, field);
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

    for (int node_index = 0; node_index <= maximum_moment; ++node_index) {
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
                    field.add(
                        residues[row_index][static_cast<std::size_t>(moment)],
                        field.multiply(
                            weights[static_cast<std::size_t>(moment)]
                                   [static_cast<std::size_t>(node_index)],
                            coefficient
                        )
                    );
            }
        }
    }

    for (auto& row : residues) {
        for (std::uint32_t& value : row) value = field.to_unsigned(value);
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

    std::vector<ResidueTable> residues_by_prime(primes.size());
    std::vector<std::string> errors(primes.size());
    const int prime_count = static_cast<int>(primes.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int prime_index = 0; prime_index < prime_count; ++prime_index) {
        const std::size_t index = static_cast<std::size_t>(prime_index);
        try {
            residues_by_prime[index] = residues_for_prime(
                rows, maximum_moment, primes[index]
            );
#ifdef _OPENMP
#pragma omp critical(active_b_frontier_progress)
#endif
            std::cout << "ACTIVE_B_FRONTIER evaluated_prime=" << primes[index]
                      << " prime_index=" << prime_index << '/' << prime_count
                      << std::endl;
        } catch (const std::exception& error) {
            errors[index] = error.what();
        } catch (...) {
            errors[index] = "unknown prime-field failure";
        }
    }
    for (std::size_t index = 0; index < errors.size(); ++index) {
        if (!errors[index].empty()) {
            throw std::runtime_error(
                std::string("prime ") + std::to_string(primes[index])
                + " failed: " + errors[index]
            );
        }
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
        const std::vector<BigInt> stable = stable_moments(
            full_q3_bcd_remaining::required_maximum_moment
        );
        std::map<int, const MomentRow*> by_rank;
        for (const MomentRow& row : rows) {
            if (!by_rank.emplace(row.rank, &row).second) {
                throw std::runtime_error("duplicate rank in active B ledger");
            }
        }

        int failures = 0;
        int cases_checked = 0;
        BigInt minimum_margin;
        bool have_minimum = false;
        std::map<int, std::pair<int, BigInt>> summaries_by_h;
        for (const active_b_frontier::Case& frontier :
             active_b_frontier::active_cases()) {
            const int rank = frontier.q - 1;
            const auto found = by_rank.find(rank);
            if (found == by_rank.end()
                || frontier.j > found->second->moment_through) {
                throw std::runtime_error("active frontier row is absent from ledger");
            }
            const BigInt finite = found->second->moments[
                static_cast<std::size_t>(frontier.j)
            ];
            const BigInt& stable_value = stable[
                static_cast<std::size_t>(frontier.j)
            ];
            const BigInt bad = stable_value - finite;
            const BigInt margin = stable_value - 2 * bad;
            const bool pass = bad >= 0 && margin > 0;
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
                          << " rank=" << rank
                          << " j=" << frontier.j
                          << " finite=" << finite
                          << " stable=" << stable_value
                          << " bad=" << bad
                          << " margin=" << margin << '\n';
            }
        }
        if (cases_checked != 337) ++failures;
        for (const auto& [h, summary] : summaries_by_h) {
            std::cout << "ACTIVE_B_FRONTIER_H_SUMMARY h=" << h
                      << " cases=" << summary.first
                      << " minimum_margin=" << summary.second << '\n';
        }
        std::cout << "ACTIVE_B_FRONTIER_BOUNDED_LITTLEWOOD rows=" << rows.size()
                  << " cases=" << cases_checked
                  << " maximum_moment="
                  << full_q3_bcd_remaining::required_maximum_moment
                  << " primes=" << primes_consumed
                  << " modulus_bits=" << modulus_bits
                  << " minimum_margin=" << minimum_margin
                  << " failures=" << failures << '\n';
        if (failures == 0) {
            std::cout << "ACTIVE_B_FRONTIER_BOUNDED_LITTLEWOOD "
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
