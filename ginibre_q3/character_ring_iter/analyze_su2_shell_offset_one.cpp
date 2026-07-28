#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Matrix = std::vector<std::vector<int>>;
using IntegerMatrix = std::vector<std::vector<Integer>>;

int parse_positive(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error("bound must be a positive integer");
    }
    return static_cast<int>(value);
}

bool fuses_half(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= source + label
        && source + target + label <= 2 * level;
}

std::vector<Integer> multiply_row(
    const std::vector<Integer>& state,
    const Matrix& matrix
) {
    std::vector<Integer> next(matrix.size());
    for (std::size_t source = 0; source < matrix.size(); ++source) {
        for (std::size_t target = 0; target < matrix.size(); ++target) {
            next[target] += state[source] * matrix[source][target];
        }
    }
    return next;
}

IntegerMatrix multiply_right_sparse(
    const IntegerMatrix& left,
    const Matrix& right
) {
    IntegerMatrix result(
        left.size(),
        std::vector<Integer>(right.size())
    );
    for (std::size_t source = 0; source < left.size(); ++source) {
        for (std::size_t middle = 0; middle < right.size(); ++middle) {
            if (left[source][middle] == 0) {
                continue;
            }
            for (std::size_t target = 0;
                 target < right.size();
                 ++target) {
                if (right[middle][target] != 0) {
                    result[source][target] +=
                        left[source][middle] * right[middle][target];
                }
            }
        }
    }
    return result;
}

Integer binomial_integer(int top, int bottom) {
    if (bottom < 0 || top < bottom) {
        return 0;
    }
    Integer result = 1;
    for (int index = 1; index <= bottom; ++index) {
        result *= top - bottom + index;
        result /= index;
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2 && argc != 3) {
            throw std::runtime_error("usage: MAXIMUM_LABEL [--table]");
        }
        const int maximum_label = parse_positive(argv[1]);
        const bool print_table =
            argc == 3 && std::string(argv[2]) == "--table";
        if (argc == 3 && !print_table) {
            throw std::runtime_error("the only optional flag is --table");
        }
        std::uint64_t parameters = 0U;
        std::uint64_t shell_entries = 0U;
        std::uint64_t negative_total = 0U;
        std::uint64_t negative_even_packet = 0U;
        std::uint64_t negative_odd_packet = 0U;
        std::uint64_t negative_odd_main_band = 0U;
        std::uint64_t negative_odd_main_nonfull_even = 0U;
        bool printed_total = false;
        bool printed_even = false;
        bool printed_odd = false;
        bool printed_odd_main = false;
        std::map<std::pair<int, int>, std::uint64_t>
            negative_odd_main_profiles;
        std::map<int, std::map<int, Integer>> exact_profile_values;
        std::map<int, std::map<std::pair<int, int>, Integer>>
            full_profile_values;
        std::map<int, std::map<int, Integer>>
            full_profile_univariate_values;

        for (int label = 1; label <= maximum_label; ++label) {
            ++parameters;
            const int level = 2 * label + 2;
            const int paired = label + 1;
            const int quotient_size = paired + 1;
            const int center = paired;
            Matrix plus(
                static_cast<std::size_t>(quotient_size),
                std::vector<int>(
                    static_cast<std::size_t>(quotient_size)
                )
            );
            Matrix minus = plus;
            for (int source = 0; source < paired; ++source) {
                for (int target = 0; target < paired; ++target) {
                    const int same = fuses_half(
                        level,
                        label,
                        source,
                        target
                    ) ? 1 : 0;
                    const int crossed = fuses_half(
                        level,
                        label,
                        source,
                        level - target
                    ) ? 1 : 0;
                    plus[static_cast<std::size_t>(source)]
                        [static_cast<std::size_t>(target)] =
                            same + crossed;
                    minus[static_cast<std::size_t>(source)]
                         [static_cast<std::size_t>(target)] =
                            same - crossed;
                }
                const int joins = fuses_half(
                    level,
                    label,
                    source,
                    level / 2
                ) ? 1 : 0;
                plus[static_cast<std::size_t>(source)]
                    [static_cast<std::size_t>(center)] = joins;
                plus[static_cast<std::size_t>(center)]
                    [static_cast<std::size_t>(source)] = 2 * joins;
            }
            plus[static_cast<std::size_t>(center)]
                [static_cast<std::size_t>(center)] = 1;

            Matrix crossing = plus;
            for (int source = 0; source < quotient_size; ++source) {
                for (int target = 0; target < quotient_size; ++target) {
                    crossing[static_cast<std::size_t>(source)]
                            [static_cast<std::size_t>(target)] -=
                        minus[static_cast<std::size_t>(source)]
                             [static_cast<std::size_t>(target)];
                }
            }
            for (int source = 0; source < paired; ++source) {
                for (int target = 0; target < paired; ++target) {
                    const int expected =
                        target == label - source
                            || target == label + 1 - source
                        ? 1
                        : 0;
                    if (
                        minus[static_cast<std::size_t>(source)]
                             [static_cast<std::size_t>(target)]
                        != expected
                    ) {
                        throw std::runtime_error(
                            "offset-one odd quotient is not JN_1"
                        );
                    }
                }
            }

            std::vector<std::vector<Integer>> plus_rows(
                6U,
                std::vector<Integer>(
                    static_cast<std::size_t>(quotient_size)
                )
            );
            plus_rows[0][0] = 1;
            for (int power = 1; power <= 5; ++power) {
                plus_rows[static_cast<std::size_t>(power)] =
                    multiply_row(
                        plus_rows[static_cast<std::size_t>(power - 1)],
                        plus
                    );
            }
            std::vector<std::vector<Integer>> prefix(
                6U,
                std::vector<Integer>(static_cast<std::size_t>(paired))
            );
            for (int power = 0; power <= 5; ++power) {
                for (int target = 0; target < paired; ++target) {
                    for (int source = 0;
                         source < quotient_size;
                         ++source) {
                        prefix[static_cast<std::size_t>(power)]
                              [static_cast<std::size_t>(target)] +=
                            plus_rows[static_cast<std::size_t>(power)]
                                     [static_cast<std::size_t>(source)]
                            * crossing[
                                static_cast<std::size_t>(source)
                              ][static_cast<std::size_t>(target)];
                    }
                }
            }
            if (print_table && label == maximum_label) {
                for (int target = 0; target < paired; ++target) {
                    std::cout
                        << "OFFSET_ONE_PREFIX"
                        << " Q=" << label
                        << " V=" << target;
                    for (int power = 1; power <= 5; ++power) {
                        std::cout
                            << " d" << power << '='
                            << prefix[
                                static_cast<std::size_t>(power)
                            ][static_cast<std::size_t>(target)];
                    }
                    std::cout << '\n';
                }
            }

            std::vector<IntegerMatrix> minus_powers(
                5U,
                IntegerMatrix(
                    static_cast<std::size_t>(paired),
                    std::vector<Integer>(
                        static_cast<std::size_t>(paired)
                    )
                )
            );
            for (int vertex = 0; vertex < paired; ++vertex) {
                minus_powers[0][static_cast<std::size_t>(vertex)]
                               [static_cast<std::size_t>(vertex)] = 1;
            }
            Matrix active_minus(
                static_cast<std::size_t>(paired),
                std::vector<int>(static_cast<std::size_t>(paired))
            );
            for (int source = 0; source < paired; ++source) {
                for (int target = 0; target < paired; ++target) {
                    active_minus[static_cast<std::size_t>(source)]
                                [static_cast<std::size_t>(target)] =
                        minus[static_cast<std::size_t>(source)]
                             [static_cast<std::size_t>(target)];
                }
            }
            for (int power = 1; power <= 4; ++power) {
                minus_powers[static_cast<std::size_t>(power)] =
                    multiply_right_sparse(
                        minus_powers[
                            static_cast<std::size_t>(power - 1)
                        ],
                        active_minus
                    );
            }

            const Integer f4 = 2 * label + 1;
            const Integer f5 =
                (
                    5 * Integer(label) * label + 5 * label + 2
                ) / 2
                - (
                    label >= 5
                    ? Integer(label - 3) * (label - 4) / 2
                    : Integer(0)
                );
            std::vector<std::vector<Integer>> even_column(
                static_cast<std::size_t>(paired),
                std::vector<Integer>(static_cast<std::size_t>(paired))
            );
            std::vector<std::vector<Integer>> odd_column = even_column;
            for (int crossing_target = 0;
                 crossing_target < paired;
                 ++crossing_target) {
                for (int target = 0; target < paired; ++target) {
                    even_column[
                        static_cast<std::size_t>(crossing_target)
                    ][static_cast<std::size_t>(target)] =
                        f4 * (
                            prefix[1][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[4][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                            + prefix[3][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[2][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                            + prefix[5][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[0][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                        )
                        - f5 * (
                            prefix[2][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[2][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                            + prefix[4][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[0][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                        );
                    odd_column[
                        static_cast<std::size_t>(crossing_target)
                    ][static_cast<std::size_t>(target)] =
                        f4 * (
                            prefix[2][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[3][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                            + prefix[4][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[1][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                        )
                        - f5 * (
                            prefix[1][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[3][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                            + prefix[3][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[1][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                        );
                }
            }

            std::vector<Integer> even_suffix(
                static_cast<std::size_t>(paired)
            );
            std::vector<Integer> odd_suffix(
                static_cast<std::size_t>(paired)
            );
            for (int rho = paired - 1; rho >= 0; --rho) {
                for (int target = 0; target < paired; ++target) {
                    ++shell_entries;
                    even_suffix[static_cast<std::size_t>(target)] +=
                        even_column[static_cast<std::size_t>(rho)]
                                   [static_cast<std::size_t>(target)];
                    odd_suffix[static_cast<std::size_t>(target)] +=
                        odd_column[static_cast<std::size_t>(rho)]
                                  [static_cast<std::size_t>(target)];
                    const Integer total =
                        even_suffix[static_cast<std::size_t>(target)]
                        + odd_suffix[static_cast<std::size_t>(target)];
                    if (label >= 7) {
                        const int even_offset = rho - target;
                        const int odd_offset =
                            rho - (label - target);
                        const int reflected_distance = label - target;
                        if (
                            odd_offset == 2
                            && even_offset >= -2
                            && even_offset <= 1
                        ) {
                            int base = 0;
                            if (even_offset == -2) {
                                base = 6;
                            } else if (
                                even_offset == -1
                                || even_offset == 0
                            ) {
                                base = 5;
                            } else {
                                base = 4;
                            }
                            exact_profile_values[even_offset][
                                target - base
                            ] = total;
                        }
                        if (
                            odd_offset == 1
                            && even_offset <= -3
                        ) {
                            const int slack =
                                label - 2 * reflected_distance - 4;
                            if (reflected_distance >= 2) {
                                full_profile_values[0][{
                                    reflected_distance - 2,
                                    slack
                                }] = total;
                            } else if (reflected_distance == 1) {
                                full_profile_univariate_values[1][
                                    slack - 1
                                ] = total;
                            } else {
                                full_profile_univariate_values[2][
                                    slack - 3
                                ] = total;
                            }
                        }
                        if (
                            odd_offset == 2
                            && even_offset <= -3
                        ) {
                            const int slack =
                                label - 2 * reflected_distance - 5;
                            if (reflected_distance >= 2) {
                                full_profile_values[3][{
                                    reflected_distance - 2,
                                    slack
                                }] = total;
                            } else if (reflected_distance == 1) {
                                full_profile_univariate_values[5][
                                    slack
                                ] = total;
                            } else {
                                full_profile_univariate_values[4][
                                    slack - 2
                                ] = total;
                            }
                        }
                    }
                    if (
                        even_suffix[static_cast<std::size_t>(target)] < 0
                    ) {
                        ++negative_even_packet;
                        if (!printed_even) {
                            printed_even = true;
                            std::cout
                                << "FIRST_NEGATIVE_OFFSET_ONE_EVEN_PACKET"
                                << " Q=" << label
                                << " rho=" << rho
                                << " target=" << target
                                << " value="
                                << even_suffix[
                                    static_cast<std::size_t>(target)
                                ] << '\n';
                        }
                    }
                    if (
                        odd_suffix[static_cast<std::size_t>(target)] < 0
                    ) {
                        ++negative_odd_packet;
                        if (label >= 7) {
                            ++negative_odd_main_band;
                            if (rho > target - 2) {
                                ++negative_odd_main_nonfull_even;
                            }
                            const int even_offset = std::clamp(
                                rho - target,
                                -3,
                                4
                            );
                            const int odd_offset = std::clamp(
                                rho - (label - target),
                                -3,
                                4
                            );
                            ++negative_odd_main_profiles[
                                {even_offset, odd_offset}
                            ];
                            if (!printed_odd_main) {
                                printed_odd_main = true;
                                std::cout
                                    << "FIRST_NEGATIVE_OFFSET_ONE_ODD_MAIN"
                                    << " Q=" << label
                                    << " rho=" << rho
                                    << " target=" << target
                                    << " even="
                                    << even_suffix[
                                        static_cast<std::size_t>(target)
                                    ]
                                    << " odd="
                                    << odd_suffix[
                                        static_cast<std::size_t>(target)
                                    ]
                                    << " total=" << total
                                    << " even_cut_offset="
                                    << rho - target
                                    << " odd_cut_offset="
                                    << rho - (label - target)
                                    << '\n';
                            }
                        }
                        if (!printed_odd) {
                            printed_odd = true;
                            std::cout
                                << "FIRST_NEGATIVE_OFFSET_ONE_ODD_PACKET"
                                << " Q=" << label
                                << " rho=" << rho
                                << " target=" << target
                                << " value="
                                << odd_suffix[
                                    static_cast<std::size_t>(target)
                                ] << '\n';
                        }
                    }
                    if (total < 0) {
                        ++negative_total;
                        if (!printed_total) {
                            printed_total = true;
                            std::cout
                                << "FIRST_NEGATIVE_OFFSET_ONE_TOTAL"
                                << " Q=" << label
                                << " rho=" << rho
                                << " target=" << target
                                << " value=" << total << '\n';
                        }
                    }
                }
            }
        }

        for (const auto& [profile, count]
             : negative_odd_main_profiles) {
            std::cout
                << "OFFSET_ONE_NEGATIVE_ODD_PROFILE"
                << " even_cut_offset=" << profile.first
                << " odd_cut_offset=" << profile.second
                << " count=" << count << '\n';
        }
        for (const auto& [even_offset, indexed_values]
             : exact_profile_values) {
            std::vector<Integer> differences;
            for (const auto& [index, value] : indexed_values) {
                if (
                    index != static_cast<int>(differences.size())
                    || differences.size() >= 12U
                ) {
                    continue;
                }
                differences.push_back(value);
            }
            if (differences.size() < 7U) {
                throw std::runtime_error(
                    "incomplete exact-profile Newton grid"
                );
            }
            int order = 0;
            std::size_t nonzero_terms = 0U;
            Integer minimum_coefficient = 0;
            bool initialized_coefficient = false;
            while (!differences.empty()) {
                const Integer coefficient = differences.front();
                if (coefficient != 0) {
                    ++nonzero_terms;
                    if (
                        !initialized_coefficient
                        || coefficient < minimum_coefficient
                    ) {
                        minimum_coefficient = coefficient;
                        initialized_coefficient = true;
                    }
                    if (coefficient < 0 || order > 5) {
                        throw std::runtime_error(
                            "invalid exact-profile Newton certificate"
                        );
                    }
                }
                std::vector<Integer> next;
                next.reserve(
                    differences.size() > 0U
                        ? differences.size() - 1U
                        : 0U
                );
                for (std::size_t index = 1U;
                     index < differences.size();
                     ++index) {
                    next.push_back(
                        differences[index] - differences[index - 1U]
                    );
                }
                differences = std::move(next);
                ++order;
            }
            std::cout
                << "OFFSET_ONE_NEWTON_CERTIFICATE"
                << " even_cut_offset=" << even_offset
                << " odd_cut_offset=2"
                << " nonzero_terms=" << nonzero_terms
                << " minimum_coefficient=" << minimum_coefficient
                << " maximum_degree=5"
                << " result=PASS_NONNEGATIVE_NEWTON_EXPANSION\n";
        }
        for (const auto& [profile, values]
             : full_profile_values) {
            constexpr int maximum_order = 6;
            std::size_t nonzero_terms = 0U;
            Integer minimum_coefficient = 0;
            bool initialized_coefficient = false;
            for (int first_order = 0;
                 first_order <= maximum_order;
                 ++first_order) {
                for (int second_order = 0;
                     second_order <= maximum_order;
                     ++second_order) {
                    Integer coefficient = 0;
                    bool complete = true;
                    for (int first = 0;
                         first <= first_order;
                         ++first) {
                        for (int second = 0;
                             second <= second_order;
                             ++second) {
                            const auto found = values.find({
                                first,
                                second
                            });
                            if (found == values.end()) {
                                complete = false;
                                continue;
                            }
                            Integer term =
                                binomial_integer(first_order, first)
                                * binomial_integer(second_order, second)
                                * found->second;
                            if (
                                (
                                    first_order - first
                                    + second_order - second
                                ) % 2 == 0
                            ) {
                                coefficient += term;
                            } else {
                                coefficient -= term;
                            }
                        }
                    }
                    if (!complete) {
                        throw std::runtime_error(
                            "incomplete full-profile Newton grid"
                        );
                    }
                    if (coefficient != 0) {
                        ++nonzero_terms;
                        if (
                            !initialized_coefficient
                            || coefficient < minimum_coefficient
                        ) {
                            minimum_coefficient = coefficient;
                            initialized_coefficient = true;
                        }
                        if (
                            coefficient < 0
                            || first_order + second_order > 5
                        ) {
                            throw std::runtime_error(
                                "invalid full-profile Newton certificate"
                            );
                        }
                    }
                }
            }
            std::cout
                << "OFFSET_ONE_FULL_NEWTON_CERTIFICATE"
                << " profile=" << profile
                << " nonzero_terms=" << nonzero_terms
                << " minimum_coefficient=" << minimum_coefficient
                << " maximum_total_degree=5"
                << " result=PASS_NONNEGATIVE_NEWTON_EXPANSION\n";
        }
        for (const auto& [profile, indexed_values]
             : full_profile_univariate_values) {
            std::vector<Integer> differences;
            for (const auto& [index, value] : indexed_values) {
                if (
                    index != static_cast<int>(differences.size())
                    || differences.size() >= 8U
                ) {
                    continue;
                }
                differences.push_back(value);
            }
            if (differences.size() < 7U) {
                throw std::runtime_error(
                    "incomplete boundary-profile Newton grid"
                );
            }
            int order = 0;
            std::size_t nonzero_terms = 0U;
            Integer minimum_coefficient = 0;
            bool initialized_coefficient = false;
            while (!differences.empty()) {
                const Integer coefficient = differences.front();
                if (coefficient != 0) {
                    ++nonzero_terms;
                    if (
                        !initialized_coefficient
                        || coefficient < minimum_coefficient
                    ) {
                        minimum_coefficient = coefficient;
                        initialized_coefficient = true;
                    }
                    if (coefficient < 0 || order > 5) {
                        throw std::runtime_error(
                            "invalid boundary-profile Newton certificate"
                        );
                    }
                }
                std::vector<Integer> next;
                for (std::size_t index = 1U;
                     index < differences.size();
                     ++index) {
                    next.push_back(
                        differences[index] - differences[index - 1U]
                    );
                }
                differences = std::move(next);
                ++order;
            }
            std::cout
                << "OFFSET_ONE_FULL_NEWTON_CERTIFICATE"
                << " profile=" << profile
                << " nonzero_terms=" << nonzero_terms
                << " minimum_coefficient=" << minimum_coefficient
                << " maximum_total_degree=5"
                << " result=PASS_NONNEGATIVE_NEWTON_EXPANSION\n";
        }
        std::cout
            << "SU2_SHELL_OFFSET_ONE"
            << " maximum_label=" << maximum_label
            << " parameters=" << parameters
            << " shell_entries=" << shell_entries
            << " negative_even_packet=" << negative_even_packet
            << " negative_odd_packet=" << negative_odd_packet
            << " negative_odd_main_band=" << negative_odd_main_band
            << " negative_odd_main_nonfull_even="
            << negative_odd_main_nonfull_even
            << " negative_total=" << negative_total
            << " result="
            << (
                negative_total == 0U
                    ? "PASS_OFFSET_ONE_SUFFIX_DISCOVERY"
                    : "FAIL_OFFSET_ONE_SUFFIX_CANDIDATE"
            )
            << '\n';
        return negative_total == 0U ? 0 : 1;
    } catch (const std::exception& exception) {
        std::cerr << "error: " << exception.what() << '\n';
        return 2;
    }
}
