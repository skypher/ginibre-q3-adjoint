// Referee-side full-scale oracle for the H25--H28 B-frontier totals.
//
// The live fleet counts column-even horizontal-two-strip paths by a backward
// state traversal.  Here the bounded Littlewood determinants independently
// reconstruct the finite B_b adjoint moments.  The exact correction identity
// then gives
//
//     bad(q,j) = stable(j) - moment(B_{q-1},j).
//
// This file deliberately reuses the separately source-audited determinant
// implementation from Part III, but supplies a frontier-specific row ledger,
// a sharper stable-moment CRT bound, and a B-only evaluation path.

#include <array>
#include <cstddef>

// Replace only the data header seen by the included implementation.  No live
// or replay-pinned source is modified.
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

// Ranks 14--40 need the H28 value 2b+30.  Ranks 41--46 occur only through
// H27 at the top of their range and need 2b+29.
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
#include "../character_ring_iter/verify_full_q3_bcd_bounded_littlewood_gmp.cpp"
#undef main

#include <algorithm>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace referee_oracle {

using ResidueTable = std::vector<std::vector<std::uint32_t>>;

ResidueTable b_moment_residues_for_prime(
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

    int b_maximum = 0;
    for (const MomentRow& row : rows) {
        if (row.family != 'B') {
            throw std::runtime_error("non-B row in B-frontier oracle");
        }
        b_maximum = std::max(b_maximum, row.rank);
    }
    const int largest_f_index = 2 * b_maximum - 1;

    for (int node_index = 0; node_index <= maximum_moment; ++node_index) {
        const std::uint32_t node = field.from_unsigned(
            static_cast<std::uint64_t>(node_index)
        );
        const NodeData data = build_node_data(
            node, largest_f_index, field, degree_limit
        );
        const auto b_minus = leading_principal_determinants(
            build_b_matrix(
                data.elementary_pair,
                b_maximum,
                false,
                field,
                degree_limit
            ),
            b_maximum,
            field,
            degree_limit
        );
        const auto b_plus = leading_principal_determinants(
            build_b_matrix(
                data.elementary_pair,
                b_maximum,
                true,
                field,
                degree_limit
            ),
            b_maximum,
            field,
            degree_limit
        );

        for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
            const MomentRow& row = rows[row_index];
            const Series first = multiply_series(
                data.elementary_sum,
                b_minus[static_cast<std::size_t>(row.rank)],
                field,
                degree_limit
            );
            const Series second = multiply_series(
                data.alternating_elementary_sum,
                b_plus[static_cast<std::size_t>(row.rank)],
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

    // The finite B moment is the unrestricted stable sum minus a nonnegative
    // excluded-width sum.  Hence the stable moment is a valid, much sharper
    // reconstruction bound than dim(G)^j for this oracle.
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
            residues_by_prime[index] = b_moment_residues_for_prime(
                rows, maximum_moment, primes[index]
            );
#ifdef _OPENMP
#pragma omp critical(full_q3_frontier_determinant_progress)
#endif
            std::cout << "FRONTIER_DETERMINANT evaluated_prime="
                      << primes[index]
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
                    "finite B moment lies outside the exact stable bound"
                );
            }
        }
    }
    return rows;
}

using Key = std::tuple<int, int, int>;  // h, q, j

std::map<Key, BigInt> read_expected(const std::string& directory) {
    const std::array<std::pair<int, const char*>, 4> files{{
        {25, "bc_twentyfifth_frontier_gmp_certificate.log"},
        {26, "bc_twentysixth_frontier_gmp_certificate.log"},
        {27, "bc_twentyseventh_frontier_gmp_certificate.log"},
        {28, "bc_twentyeighth_frontier_gmp_certificate.log"},
    }};
    const std::regex pattern(
        R"(^B_BAD_COUNT q=([0-9]+) j=([0-9]+) value=([0-9]+)$)"
    );
    std::map<Key, BigInt> expected;
    for (const auto& [h, name] : files) {
        std::ifstream input(directory + "/" + name);
        if (!input) throw std::runtime_error(std::string("cannot open ") + name);
        std::string line;
        while (std::getline(input, line)) {
            std::smatch match;
            if (!std::regex_match(line, match, pattern)) continue;
            const int q = std::stoi(match[1].str());
            const int j = std::stoi(match[2].str());
            const BigInt value(match[3].str(), 10);
            const Key key{h, q, j};
            if (!expected.emplace(key, value).second) {
                throw std::runtime_error("duplicate accepted B bad-count row");
            }
        }
    }
    if (expected.size() != 120U) {
        throw std::runtime_error(
            "accepted B bad-count scope changed: "
            + std::to_string(expected.size())
        );
    }
    return expected;
}

}  // namespace referee_oracle

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: verify_frontier_bounded_littlewood_oracle "
                "CERTIFICATES_POST_M29_DIRECTORY"
            );
        }
        std::cout << std::unitbuf;
        std::size_t primes_consumed = 0U;
        int modulus_bits = 0;
        const std::vector<MomentRow> rows =
            referee_oracle::reconstruct_moments(
                full_q3_bcd_remaining::required_maximum_moment,
                primes_consumed,
                modulus_bits
            );
        const std::vector<BigInt> stable = stable_moments(
            full_q3_bcd_remaining::required_maximum_moment
        );
        std::map<int, const MomentRow*> by_rank;
        for (const MomentRow& row : rows) by_rank.emplace(row.rank, &row);
        const std::map<referee_oracle::Key, BigInt> expected =
            referee_oracle::read_expected(argv[1]);

        int failures = 0;
        int cases = 0;
        for (const auto& [key, accepted] : expected) {
            const auto [h, q, j] = key;
            const int rank = q - 1;
            const auto found = by_rank.find(rank);
            if (found == by_rank.end()
                || j > found->second->moment_through) {
                throw std::runtime_error("frontier row is absent from oracle ledger");
            }
            const BigInt finite =
                found->second->moments[static_cast<std::size_t>(j)];
            const BigInt bad = stable[static_cast<std::size_t>(j)] - finite;
            if (bad != accepted) ++failures;
            ++cases;
            std::cout << "B_DETERMINANT_BAD_COUNT h=" << h
                      << " q=" << q
                      << " rank=" << rank
                      << " j=" << j
                      << " finite=" << finite
                      << " stable=" << stable[static_cast<std::size_t>(j)]
                      << " bad=" << bad
                      << " accepted=" << accepted
                      << " match=" << (bad == accepted ? 1 : 0) << '\n';
        }
        if (cases != 120) ++failures;
        std::cout << "FRONTIER_BOUNDED_LITTLEWOOD_ORACLE rows=" << rows.size()
                  << " cases=" << cases
                  << " maximum_moment="
                  << full_q3_bcd_remaining::required_maximum_moment
                  << " primes=" << primes_consumed
                  << " modulus_bits=" << modulus_bits
                  << " failures=" << failures << '\n';
        if (failures == 0) {
            std::cout
                << "FRONTIER_BOUNDED_LITTLEWOOD_ORACLE VERIFICATION: ALL PASS\n";
        }
        std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << '\n';
        return failures == 0 ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "FRONTIER_BOUNDED_LITTLEWOOD_ORACLE FAILURE: "
                  << error.what() << '\n';
        std::cout << "__EXIT_STATUS=1\n";
        return 1;
    }
}
