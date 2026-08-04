#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <deque>
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
    if (
        end == text
        || *end != '\0'
        || value <= 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error("argument must be a positive integer");
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
    std::vector<Integer> result(matrix.size());
    for (std::size_t source = 0; source < matrix.size(); ++source) {
        if (state[source] == 0) {
            continue;
        }
        for (std::size_t target = 0; target < matrix.size(); ++target) {
            if (matrix[source][target] != 0) {
                result[target] += state[source] * matrix[source][target];
            }
        }
    }
    return result;
}

IntegerMatrix multiply_right(
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
            for (std::size_t target = 0; target < right.size(); ++target) {
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
        if (argc != 3 && argc != 4) {
            throw std::runtime_error(
                "usage: MAXIMUM_Q OFFSET [--certificate|--residue-tails|"
                "--h2-residue-tails|--h2-family]"
            );
        }
        const int maximum_q = parse_positive(argv[1]);
        const int offset = parse_positive(argv[2]);
        const bool certificate =
            argc == 4 && std::string(argv[3]) == "--certificate";
        const bool residue_tails =
            argc == 4 && std::string(argv[3]) == "--residue-tails";
        const bool h2_residue_tails =
            argc == 4 && std::string(argv[3]) == "--h2-residue-tails";
        const bool h2_family =
            argc == 4 && std::string(argv[3]) == "--h2-family";
        if (
            argc == 4
            && !certificate
            && !residue_tails
            && !h2_residue_tails
            && !h2_family
        ) {
            throw std::runtime_error(
                "the only optional flag is --certificate, --residue-tails, "
                "--h2-residue-tails, or --h2-family"
            );
        }
        const int certificate_maximum_q =
            offset <= 7 ? 100 : (offset <= 9 ? 120 : 140);
        if (
            certificate
            && (
                offset < 3
                || offset > 11
                || maximum_q != certificate_maximum_q
            )
        ) {
            throw std::runtime_error(
                "certificate mode requires the canonical bound for OFFSET"
            );
        }
        if (offset < 1) {
            throw std::runtime_error("offset must be positive");
        }
        const int gap = offset - 1;
        const int even_radius = 2 * gap;
        const int reflected_shift = gap / 2;
        const int odd_lower =
            gap % 2 == 0
                ? -3 * gap / 2
                : (1 - 3 * gap) / 2;
        const int odd_upper =
            gap % 2 == 0
                ? 3 * gap / 2
                : (1 + 3 * gap) / 2;
        const int stable_boundary = 4 * gap + 1;

        std::uint64_t entries = 0U;
        std::uint64_t residue_tail_entries = 0U;
        std::uint64_t negative_residue_tails = 0U;
        bool printed_negative_residue_tail = false;
        std::uint64_t h2_residue_tail_entries = 0U;
        std::uint64_t negative_h2_residue_tails = 0U;
        std::uint64_t negative_h2_atoms = 0U;
        int maximum_h2_payment_span = 0;
        bool has_maximum_h2_payment_span = false;
        int maximum_h2_payment_q = -1;
        int maximum_h2_payment_target = -1;
        int maximum_h2_payment_residue = -1;
        int maximum_h2_payment_debit = -1;
        int maximum_h2_payment_credit = -1;
        std::uint64_t negative_even = 0U;
        std::uint64_t negative_odd = 0U;
        std::uint64_t negative_total = 0U;
        int maximum_even_radius = 0;
        int maximum_odd_radius = 0;
        bool printed_even = false;
        bool printed_odd = false;
        bool printed_total = false;
        std::map<std::pair<int, int>, std::uint64_t>
            negative_even_profiles;
        std::map<std::pair<int, int>, std::uint64_t>
            negative_odd_profiles;
        std::map<std::string, std::map<int, Integer>>
            classification_ray_values;
        std::map<std::string, std::map<std::pair<int, int>, Integer>>
            classification_cone_values;
        std::map<std::string, std::pair<std::uint64_t, Integer>>
            classification_box_values;

        for (int q = 1; q <= maximum_q; ++q) {
            const int level = 2 * q + offset;
            const int paired = (level + 1) / 2;
            const bool has_center = level % 2 == 0;
            const int quotient_size = paired + (has_center ? 1 : 0);
            const int center = has_center ? paired : -1;
            Matrix plus(
                static_cast<std::size_t>(quotient_size),
                std::vector<int>(static_cast<std::size_t>(quotient_size))
            );
            Matrix minus(
                static_cast<std::size_t>(paired),
                std::vector<int>(static_cast<std::size_t>(paired))
            );
            for (int source = 0; source < paired; ++source) {
                for (int target = 0; target < paired; ++target) {
                    const int same = fuses_half(
                        level,
                        q,
                        source,
                        target
                    ) ? 1 : 0;
                    const int crossed = fuses_half(
                        level,
                        q,
                        source,
                        level - target
                    ) ? 1 : 0;
                    plus[static_cast<std::size_t>(source)]
                        [static_cast<std::size_t>(target)] =
                            same + crossed;
                    minus[static_cast<std::size_t>(source)]
                         [static_cast<std::size_t>(target)] =
                            same - crossed;
                    if (
                        minus[static_cast<std::size_t>(source)]
                             [static_cast<std::size_t>(target)] < 0
                    ) {
                        throw std::runtime_error(
                            "negative odd quotient entry"
                        );
                    }
                }
                if (has_center) {
                    const int joins = fuses_half(
                        level,
                        q,
                        source,
                        level / 2
                    ) ? 1 : 0;
                    plus[static_cast<std::size_t>(source)]
                        [static_cast<std::size_t>(center)] = joins;
                    plus[static_cast<std::size_t>(center)]
                        [static_cast<std::size_t>(source)] = 2 * joins;
                }
            }
            if (has_center) {
                plus[static_cast<std::size_t>(center)]
                    [static_cast<std::size_t>(center)] = 1;
            }

            Matrix crossing(
                static_cast<std::size_t>(quotient_size),
                std::vector<int>(static_cast<std::size_t>(paired))
            );
            for (int source = 0; source < quotient_size; ++source) {
                for (int target = 0; target < paired; ++target) {
                    crossing[static_cast<std::size_t>(source)]
                            [static_cast<std::size_t>(target)] =
                        plus[static_cast<std::size_t>(source)]
                            [static_cast<std::size_t>(target)]
                        - (
                            source < paired
                                ? minus[
                                    static_cast<std::size_t>(source)
                                  ][static_cast<std::size_t>(target)]
                                : 0
                        );
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

            std::vector<Integer> minus_root(
                static_cast<std::size_t>(paired)
            );
            minus_root[0] = 1;
            Integer minus_four = 0;
            Integer minus_five = 0;
            for (int power = 1; power <= 5; ++power) {
                minus_root = multiply_row(minus_root, minus);
                if (power == 4) {
                    minus_four = minus_root[0];
                } else if (power == 5) {
                    minus_five = minus_root[0];
                }
            }
            const Integer f4 = (
                plus_rows[4][0] + minus_four
            ) / 2;
            const Integer f5 = (
                plus_rows[5][0] + minus_five
            ) / 2;

            std::vector<IntegerMatrix> minus_powers(
                5U,
                IntegerMatrix(
                    static_cast<std::size_t>(paired),
                    std::vector<Integer>(static_cast<std::size_t>(paired))
                )
            );
            for (int vertex = 0; vertex < paired; ++vertex) {
                minus_powers[0][static_cast<std::size_t>(vertex)]
                               [static_cast<std::size_t>(vertex)] = 1;
            }
            for (int power = 1; power <= 4; ++power) {
                minus_powers[static_cast<std::size_t>(power)] =
                    multiply_right(
                        minus_powers[
                            static_cast<std::size_t>(power - 1)
                        ],
                        minus
                    );
            }

            IntegerMatrix even_column(
                static_cast<std::size_t>(paired),
                std::vector<Integer>(static_cast<std::size_t>(paired))
            );
            IntegerMatrix odd_column = even_column;
            IntegerMatrix h2_column = even_column;
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
                    h2_column[
                        static_cast<std::size_t>(crossing_target)
                    ][static_cast<std::size_t>(target)] =
                        f4 * (
                            prefix[3][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[2][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                            + prefix[2][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[3][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                            + prefix[1][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[4][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                        ) - f5 * (
                            prefix[2][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[2][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                            + prefix[1][
                                static_cast<std::size_t>(crossing_target)
                            ] * minus_powers[3][
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(target)]
                        );
                    if (
                        even_column[
                            static_cast<std::size_t>(crossing_target)
                        ][static_cast<std::size_t>(target)] != 0
                    ) {
                        maximum_even_radius = std::max(
                            maximum_even_radius,
                            std::abs(crossing_target - target)
                        );
                    }
                    if (
                        odd_column[
                            static_cast<std::size_t>(crossing_target)
                        ][static_cast<std::size_t>(target)] != 0
                    ) {
                        const int reflected = paired - 1 - target;
                        maximum_odd_radius = std::max(
                            maximum_odd_radius,
                            std::abs(crossing_target - reflected)
                        );
                    }
                }
            }

            if ((h2_residue_tails || h2_family) && q >= 7 && offset >= 12) {
                struct Reserve {
                    int crossing = 0;
                    Integer amount = 0;
                };
                for (int target = 0; target < paired; ++target) {
                    for (int residue = 0; residue < 4; ++residue) {
                        Integer tail = 0;
                        std::deque<Reserve> reserve;
                        int crossing = paired - 1;
                        while (crossing >= 0 && crossing % 4 != residue) {
                            --crossing;
                        }
                        for (; crossing >= 0; crossing -= 4) {
                            const Integer& value = h2_column[
                                static_cast<std::size_t>(crossing)
                            ][static_cast<std::size_t>(target)];
                            tail += value;
                            ++h2_residue_tail_entries;
                            if (tail < 0) {
                                ++negative_h2_residue_tails;
                            }
                            if (value >= 0) {
                                if (value != 0) {
                                    reserve.push_back(Reserve{crossing, value});
                                }
                                continue;
                            }
                            ++negative_h2_atoms;
                            Integer debt = -value;
                            while (debt != 0) {
                                if (reserve.empty()) {
                                    throw std::runtime_error(
                                        "negative H2 residue atom lacks reserve"
                                    );
                                }
                                Reserve& nearest = reserve.back();
                                const Integer payment = std::min(
                                    nearest.amount, debt
                                );
                                nearest.amount -= payment;
                                debt -= payment;
                                const int span = nearest.crossing - crossing;
                                if (span > maximum_h2_payment_span) {
                                    maximum_h2_payment_span = span;
                                    has_maximum_h2_payment_span = true;
                                    maximum_h2_payment_q = q;
                                    maximum_h2_payment_target = target;
                                    maximum_h2_payment_residue = residue;
                                    maximum_h2_payment_debit = crossing;
                                    maximum_h2_payment_credit = nearest.crossing;
                                }
                                if (nearest.amount == 0) {
                                    reserve.pop_back();
                                }
                            }
                        }
                    }
                }
            }

            if (
                h2_family
                && offset == 12
                && q >= 15
                && q % 4 == 3
            ) {
                const int target = q + 5;
                const int debit = 12;
                const int credit = q - 3;
                if (target >= paired || credit >= paired) {
                    throw std::runtime_error("H2 family target is out of range");
                }
                Integer tail = 0;
                for (int crossing = paired - 1;
                     crossing >= debit;
                     --crossing) {
                    if (crossing % 4 == 0) {
                        tail += h2_column[
                            static_cast<std::size_t>(crossing)
                        ][static_cast<std::size_t>(target)];
                    }
                }
                std::cout
                    << "H2_D11_FAMILY"
                    << " Q=" << q
                    << " x=" << target
                    << " f4=" << f4
                    << " f5=" << f5
                    << " debit=" << debit
                    << " debit_value=" << h2_column[
                        static_cast<std::size_t>(debit)
                    ][static_cast<std::size_t>(target)]
                    << " debit_d=("
                    << prefix[1][static_cast<std::size_t>(debit)] << ','
                    << prefix[2][static_cast<std::size_t>(debit)] << ','
                    << prefix[3][static_cast<std::size_t>(debit)] << ')'
                    << " debit_p=("
                    << minus_powers[2U][static_cast<std::size_t>(debit)]
                       [static_cast<std::size_t>(target)] << ','
                    << minus_powers[3U][static_cast<std::size_t>(debit)]
                       [static_cast<std::size_t>(target)] << ','
                    << minus_powers[4U][static_cast<std::size_t>(debit)]
                       [static_cast<std::size_t>(target)] << ')'
                    << " credit=" << credit
                    << " credit_value=" << h2_column[
                        static_cast<std::size_t>(credit)
                    ][static_cast<std::size_t>(target)]
                    << " credit_d=("
                    << prefix[1][static_cast<std::size_t>(credit)] << ','
                    << prefix[2][static_cast<std::size_t>(credit)] << ','
                    << prefix[3][static_cast<std::size_t>(credit)] << ')'
                    << " credit_p=("
                    << minus_powers[2U][static_cast<std::size_t>(credit)]
                       [static_cast<std::size_t>(target)] << ','
                    << minus_powers[3U][static_cast<std::size_t>(credit)]
                       [static_cast<std::size_t>(target)] << ','
                    << minus_powers[4U][static_cast<std::size_t>(credit)]
                       [static_cast<std::size_t>(target)] << ')'
                    << " tail_at_debit=" << tail << '\n';
            }

            std::vector<Integer> even_suffix(
                static_cast<std::size_t>(paired)
            );
            std::vector<Integer> odd_suffix(
                static_cast<std::size_t>(paired)
            );
            std::array<std::vector<Integer>, 4U> residue_suffix{};
            for (std::vector<Integer>& suffix : residue_suffix) {
                suffix.assign(static_cast<std::size_t>(paired), 0);
            }
            for (int rho = paired - 1; rho >= 0; --rho) {
                for (int target = 0; target < paired; ++target) {
                    ++entries;
                    even_suffix[static_cast<std::size_t>(target)] +=
                        even_column[static_cast<std::size_t>(rho)]
                                   [static_cast<std::size_t>(target)];
                    odd_suffix[static_cast<std::size_t>(target)] +=
                        odd_column[static_cast<std::size_t>(rho)]
                                  [static_cast<std::size_t>(target)];
                    const Integer total =
                        even_suffix[static_cast<std::size_t>(target)]
                        + odd_suffix[static_cast<std::size_t>(target)];
                    std::vector<Integer>& current_residue_tail =
                        residue_suffix[static_cast<std::size_t>(rho % 4)];
                    current_residue_tail[static_cast<std::size_t>(target)] +=
                        even_column[static_cast<std::size_t>(rho)]
                                   [static_cast<std::size_t>(target)]
                        + odd_column[static_cast<std::size_t>(rho)]
                                  [static_cast<std::size_t>(target)];
                    if (
                        residue_tails
                        && q >= 7
                        && offset >= 11
                    ) {
                        ++residue_tail_entries;
                        if (
                            current_residue_tail[
                                static_cast<std::size_t>(target)
                            ] < 0
                        ) {
                            ++negative_residue_tails;
                            if (!printed_negative_residue_tail) {
                                printed_negative_residue_tail = true;
                                std::cout
                                    << "FIRST_NEGATIVE_FIXED_OFFSET_RESIDUE_TAIL"
                                    << " Q=" << q
                                    << " offset=" << offset
                                    << " rho=" << rho
                                    << " target=" << target
                                    << " residue=" << rho % 4
                                    << " value="
                                    << current_residue_tail[
                                        static_cast<std::size_t>(target)
                                    ] << '\n';
                            }
                        }
                    }
                    const int reflected = paired - 1 - target;
                    const std::pair<int, int> profile{
                        rho - target,
                        rho - reflected
                    };
                    if (
                        (offset >= 3 && offset <= 11)
                        && q >= 7
                    ) {
                        const int x = target;
                        const int y = q - target;
                        const int a = profile.first;
                        const int b = profile.second;
                        const auto record_box =
                            [&classification_box_values](
                                const std::string& key,
                                const Integer& value
                            ) {
                                auto& record =
                                    classification_box_values[key];
                                if (
                                    record.first == 0U
                                    || value < record.second
                                ) {
                                    record.second = value;
                                }
                                ++record.first;
                            };
                        const auto record_rectangular =
                            [
                                &classification_ray_values,
                                &classification_cone_values,
                                &record_box,
                                x,
                                y,
                                stable_boundary
                            ](
                                const std::string& key,
                                const Integer& value
                            ) {
                                if (
                                    x < stable_boundary
                                    && y < stable_boundary
                                ) {
                                    record_box(key, value);
                                } else if (x < stable_boundary) {
                                    classification_ray_values[
                                        key + "_x" + std::to_string(x)
                                    ][y - stable_boundary] = value;
                                } else if (y < stable_boundary) {
                                    classification_ray_values[
                                        key + "_y" + std::to_string(y)
                                    ][x - stable_boundary] = value;
                                } else {
                                    classification_cone_values[key][{
                                        x - stable_boundary,
                                        y - stable_boundary
                                    }] = value;
                                }
                            };
                        if (
                            a >= -even_radius
                            && a <= even_radius
                        ) {
                            record_rectangular(
                                "E_a" + std::to_string(a),
                                even_suffix[
                                    static_cast<std::size_t>(target)
                                ]
                            );
                        }
                        if (b == odd_lower) {
                            record_rectangular(
                                "O_b" + std::to_string(b),
                                odd_suffix[
                                    static_cast<std::size_t>(target)
                                ]
                            );
                        }
                        if (
                            b >= odd_lower + 1
                            && b <= odd_upper
                        ) {
                            const std::string base_key =
                                "T_b" + std::to_string(b);
                            if (a <= -even_radius) {
                                const int slack =
                                    -even_radius - a;
                                if (y < stable_boundary) {
                                    if (x < stable_boundary) {
                                        record_box(
                                            base_key + "_afull",
                                            total
                                        );
                                    } else {
                                        classification_ray_values[
                                            base_key + "_afull_y"
                                            + std::to_string(y)
                                        ][
                                            x - std::max(
                                                stable_boundary,
                                                y + reflected_shift
                                                + b + even_radius
                                            )
                                        ] = total;
                                    }
                                } else {
                                    classification_cone_values[
                                        base_key + "_afull"
                                    ][{
                                        y - stable_boundary,
                                        slack
                                    }] = total;
                                }
                            } else if (a <= even_radius) {
                                const int difference =
                                    reflected_shift + b - a;
                                const int stable_y = std::max(
                                    stable_boundary,
                                    stable_boundary - difference
                                );
                                if (y < stable_y) {
                                    record_box(
                                        base_key + "_a"
                                        + std::to_string(a),
                                        total
                                    );
                                } else {
                                    classification_ray_values[
                                        base_key + "_a"
                                        + std::to_string(a)
                                    ][y - stable_y] = total;
                                }
                            } else {
                                const int slack =
                                    a - even_radius - 1;
                                if (x < stable_boundary) {
                                    if (y < stable_boundary) {
                                        record_box(
                                            base_key + "_aempty",
                                            total
                                        );
                                    } else {
                                        classification_ray_values[
                                            base_key + "_aempty_x"
                                            + std::to_string(x)
                                        ][
                                            y - std::max(
                                                stable_boundary,
                                                x + even_radius + 1
                                                - reflected_shift - b
                                            )
                                        ] = total;
                                    }
                                } else {
                                    classification_cone_values[
                                        base_key + "_aempty"
                                    ][{
                                        x - stable_boundary,
                                        slack
                                    }] = total;
                                }
                            }
                        }
                    }
                    if (
                        even_suffix[static_cast<std::size_t>(target)] < 0
                    ) {
                        ++negative_even;
                        ++negative_even_profiles[profile];
                        if (!printed_even) {
                            printed_even = true;
                            std::cout
                                << "FIRST_NEGATIVE_FIXED_OFFSET_EVEN"
                                << " Q=" << q
                                << " offset=" << offset
                                << " rho=" << rho
                                << " target=" << target
                                << " profile=(" << profile.first
                                << ',' << profile.second << ')'
                                << " value="
                                << even_suffix[
                                    static_cast<std::size_t>(target)
                                ] << '\n';
                        }
                    }
                    if (
                        odd_suffix[static_cast<std::size_t>(target)] < 0
                    ) {
                        ++negative_odd;
                        ++negative_odd_profiles[profile];
                        if (!printed_odd) {
                            printed_odd = true;
                            std::cout
                                << "FIRST_NEGATIVE_FIXED_OFFSET_ODD"
                                << " Q=" << q
                                << " offset=" << offset
                                << " rho=" << rho
                                << " target=" << target
                                << " profile=(" << profile.first
                                << ',' << profile.second << ')'
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
                                << "FIRST_NEGATIVE_FIXED_OFFSET_TOTAL"
                                << " Q=" << q
                                << " offset=" << offset
                                << " rho=" << rho
                                << " target=" << target
                                << " profile=(" << profile.first
                                << ',' << profile.second << ')'
                                << " value=" << total << '\n';
                        }
                    }
                }
            }
        }

        if (!certificate && !h2_family) {
            for (const auto& [profile, count]
                 : negative_even_profiles) {
                std::cout
                    << "FIXED_OFFSET_NEGATIVE_EVEN_PROFILE"
                    << " a=" << profile.first
                    << " b=" << profile.second
                    << " count=" << count << '\n';
            }
            for (const auto& [profile, count]
                 : negative_odd_profiles) {
                std::cout
                    << "FIXED_OFFSET_NEGATIVE_ODD_PROFILE"
                    << " a=" << profile.first
                    << " b=" << profile.second
                    << " count=" << count << '\n';
            }
        } else {
            struct BandSummary {
                int minimum_a = 0;
                int maximum_a = 0;
                std::uint64_t profiles = 0U;
                std::uint64_t entries = 0U;
            };
            std::map<int, BandSummary> bands;
            for (const auto& [profile, count]
                 : negative_odd_profiles) {
                auto& band = bands[profile.second];
                if (band.profiles == 0U) {
                    band.minimum_a = profile.first;
                    band.maximum_a = profile.first;
                } else {
                    band.minimum_a = std::min(
                        band.minimum_a,
                        profile.first
                    );
                    band.maximum_a = std::max(
                        band.maximum_a,
                        profile.first
                    );
                }
                ++band.profiles;
                band.entries += count;
            }
            for (const auto& [b, band] : bands) {
                std::cout
                    << "FIXED_OFFSET_NEGATIVE_ODD_BAND"
                    << " b=" << b
                    << " minimum_a=" << band.minimum_a
                    << " maximum_a=" << band.maximum_a
                    << " profiles=" << band.profiles
                    << " entries=" << band.entries << '\n';
            }
        }
        if (
            offset >= 3
            && (offset <= 10 || certificate)
        ) {
            std::size_t ray_certificates = 0U;
            Integer ray_minimum = 0;
            bool initialized_ray = false;
            for (const auto& [key, values]
                 : classification_ray_values) {
                std::vector<Integer> coefficients;
                for (int order = 0; order <= 9; ++order) {
                    Integer coefficient = 0;
                    for (int index = 0; index <= order; ++index) {
                        const auto found = values.find(index);
                        if (found == values.end()) {
                            throw std::runtime_error(
                                "incomplete fixed-offset ray: " + key
                            );
                        }
                        const Integer term =
                            binomial_integer(order, index)
                            * found->second;
                        if ((order - index) % 2 == 0) {
                            coefficient += term;
                        } else {
                            coefficient -= term;
                        }
                    }
                    if (
                        (order <= 8 && coefficient < 0)
                        || (order == 9 && coefficient != 0)
                    ) {
                        std::cerr
                            << "FAILED_FIXED_OFFSET_RAY"
                            << " key=" << key
                            << " order=" << order
                            << " coefficient=" << coefficient;
                        for (int sample = 0; sample <= 9; ++sample) {
                            const auto found = values.find(sample);
                            if (found != values.end()) {
                                std::cerr
                                    << " v" << sample << '='
                                    << found->second;
                            }
                        }
                        std::cerr << '\n';
                        throw std::runtime_error(
                            "invalid fixed-offset ray certificate"
                        );
                    }
                    if (order <= 8) {
                        coefficients.push_back(coefficient);
                        if (
                            !initialized_ray
                            || coefficient < ray_minimum
                        ) {
                            ray_minimum = coefficient;
                            initialized_ray = true;
                        }
                    }
                }
                for (const auto& [index, value] : values) {
                    Integer reconstructed = 0;
                    for (
                        std::size_t order = 0U;
                        order < coefficients.size();
                        ++order
                    ) {
                        reconstructed +=
                            coefficients[order]
                            * binomial_integer(
                                index,
                                static_cast<int>(order)
                            );
                    }
                    if (reconstructed != value) {
                        throw std::runtime_error(
                            "fixed-offset ray reconstruction mismatch: "
                            + key
                        );
                    }
                }
                ++ray_certificates;
            }

            std::size_t cone_certificates = 0U;
            Integer cone_minimum = 0;
            bool initialized_cone = false;
            for (const auto& [key, values]
                 : classification_cone_values) {
                std::map<std::pair<int, int>, Integer> coefficients;
                for (int first_order = 0;
                     first_order <= 9;
                     ++first_order) {
                    for (int second_order = 0;
                         second_order <= 9 - first_order;
                         ++second_order) {
                        Integer coefficient = 0;
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
                                    throw std::runtime_error(
                                        "incomplete fixed-offset cone: "
                                        + key
                                    );
                                }
                                Integer term =
                                    binomial_integer(
                                        first_order,
                                        first
                                    )
                                    * binomial_integer(
                                        second_order,
                                        second
                                    )
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
                        const int total_order =
                            first_order + second_order;
                        if (
                            (total_order <= 8 && coefficient < 0)
                            || (
                                total_order == 9
                                && coefficient != 0
                            )
                        ) {
                            std::cerr
                                << "FAILED_FIXED_OFFSET_CONE"
                                << " key=" << key
                                << " first_order=" << first_order
                                << " second_order=" << second_order
                                << " coefficient=" << coefficient
                                << '\n';
                            throw std::runtime_error(
                                "invalid fixed-offset cone certificate"
                            );
                        }
                        if (total_order <= 8) {
                            coefficients[{
                                first_order,
                                second_order
                            }] = coefficient;
                            if (
                                !initialized_cone
                                || coefficient < cone_minimum
                            ) {
                                cone_minimum = coefficient;
                                initialized_cone = true;
                            }
                        }
                    }
                }
                for (const auto& [index, value] : values) {
                    Integer reconstructed = 0;
                    for (const auto& [order, coefficient]
                         : coefficients) {
                        reconstructed +=
                            coefficient
                            * binomial_integer(
                                index.first,
                                order.first
                            )
                            * binomial_integer(
                                index.second,
                                order.second
                            );
                    }
                    if (reconstructed != value) {
                        throw std::runtime_error(
                            "fixed-offset cone reconstruction mismatch: "
                            + key
                        );
                    }
                }
                ++cone_certificates;
            }

            std::uint64_t box_entries = 0U;
            Integer box_minimum = 0;
            bool initialized_box = false;
            for (const auto& [key, record]
                 : classification_box_values) {
                if (record.second < 0) {
                    throw std::runtime_error(
                        "negative fixed-offset boundary box: " + key
                    );
                }
                box_entries += record.first;
                if (
                    !initialized_box
                    || record.second < box_minimum
                ) {
                    box_minimum = record.second;
                    initialized_box = true;
                }
            }
            bool valid_census = false;
            if (offset == 3) {
                valid_census =
                        box_entries == 906U
                        && ray_certificates == 320U
                        && cone_certificates == 22U
                        && maximum_even_radius == 4
                        && maximum_odd_radius == 3
                        && negative_even == 0U
                        && negative_odd == 5682U
                        && negative_total == 0U
                        && negative_odd_profiles.size() == 482U
                        && entries == 358950U;
            } else if (offset == 4) {
                valid_census =
                        box_entries == 3289U
                        && ray_certificates == 662U
                        && cone_certificates == 32U
                        && maximum_even_radius == 6
                        && maximum_odd_radius == 5
                        && negative_even == 0U
                        && negative_odd == 8174U
                        && negative_total == 0U
                        && negative_odd_profiles.size() == 588U
                        && entries == 358950U;
            } else if (offset == 5) {
                valid_census =
                    box_entries == 8032U
                    && ray_certificates == 1158U
                    && cone_certificates == 42U
                    && maximum_even_radius == 8
                    && maximum_odd_radius == 6
                    && negative_even == 0U
                    && negative_odd == 11228U
                    && negative_total == 0U
                    && negative_odd_profiles.size() == 802U
                    && entries == 369550U;
            } else if (offset == 6) {
                valid_census =
                    box_entries == 15336U
                    && ray_certificates == 1754U
                    && cone_certificates == 52U
                    && maximum_even_radius == 10
                    && maximum_odd_radius == 8
                    && negative_even == 0U
                    && negative_odd == 13518U
                    && negative_total == 0U
                    && negative_odd_profiles.size() == 913U
                    && entries == 369550U;
            } else if (offset == 7) {
                valid_census =
                    box_entries == 26745U
                    && ray_certificates == 2518U
                    && cone_certificates == 62U
                    && maximum_even_radius == 12
                    && maximum_odd_radius == 9
                    && negative_even == 0U
                    && negative_odd == 16636U
                    && negative_total == 0U
                    && negative_odd_profiles.size() == 1103U
                    && entries == 380350U;
            } else if (offset == 8) {
                valid_census =
                    box_entries == 41584U
                    && ray_certificates == 3368U
                    && cone_certificates == 72U
                    && maximum_even_radius == 14
                    && maximum_odd_radius == 11
                    && negative_even == 0U
                    && negative_odd == 27056U
                    && negative_total == 0U
                    && negative_odd_profiles.size() == 1503U
                    && entries == 643220U;
            } else if (offset == 9) {
                valid_census =
                    box_entries == 62498U
                    && ray_certificates == 4400U
                    && cone_certificates == 82U
                    && maximum_even_radius == 16
                    && maximum_odd_radius == 12
                    && negative_even == 0U
                    && negative_odd == 31586U
                    && negative_total == 0U
                    && negative_odd_profiles.size() == 1730U
                    && entries == 658820U;
            } else if (offset == 10) {
                valid_census =
                    box_entries == 87509U
                    && ray_certificates == 5504U
                    && cone_certificates == 92U
                    && maximum_even_radius == 18
                    && maximum_odd_radius == 14
                    && negative_even == 0U
                    && negative_odd == 47160U
                    && negative_total == 0U
                    && negative_odd_profiles.size() == 2225U
                    && entries == 1026690U;
            } else if (offset == 11) {
                valid_census =
                    box_entries == 120772U
                    && ray_certificates == 6804U
                    && cone_certificates == 102U
                    && maximum_even_radius == 20
                    && maximum_odd_radius == 15
                    && negative_even == 0U
                    && negative_odd == 53250U
                    && negative_total == 0U
                    && negative_odd_profiles.size() == 2481U
                    && entries == 1047970U;
            }
            if (certificate && !valid_census) {
                throw std::runtime_error(
                    "incomplete fixed-offset certificate census"
                );
            }
            std::cout
                << "FIXED_OFFSET_CERTIFICATE"
                << " offset=" << offset
                << " boundary_box_entries=" << box_entries
                << " boundary_box_minimum=" << box_minimum
                << " ray_certificates=" << ray_certificates
                << " ray_minimum_coefficient=" << ray_minimum
                << " cone_certificates=" << cone_certificates
                << " cone_minimum_coefficient=" << cone_minimum
                << " maximum_total_degree=8"
                << " result=PASS_NONNEGATIVE_NEWTON_EXPANSIONS\n";
        }
        std::cout
            << "SU2_SHELL_FIXED_OFFSET"
            << " maximum_q=" << maximum_q
            << " offset=" << offset
            << " entries=" << entries
            << " maximum_even_radius=" << maximum_even_radius
            << " maximum_odd_radius=" << maximum_odd_radius
            << " negative_even=" << negative_even
            << " negative_odd=" << negative_odd
            << " negative_total=" << negative_total
            << " even_profiles=" << negative_even_profiles.size()
            << " odd_profiles=" << negative_odd_profiles.size()
            << " residue_tail_entries=" << residue_tail_entries
            << " negative_residue_tails=" << negative_residue_tails
            << " h2_residue_tail_entries=" << h2_residue_tail_entries
            << " negative_h2_residue_tails=" << negative_h2_residue_tails
            << " negative_h2_atoms=" << negative_h2_atoms
            << " maximum_h2_payment_span=" << maximum_h2_payment_span
            << " maximum_h2_payment_witness=";
        if (has_maximum_h2_payment_span) {
            std::cout
                << "Q:" << maximum_h2_payment_q
                << ",target:" << maximum_h2_payment_target
                << ",residue:" << maximum_h2_payment_residue
                << ",debit:" << maximum_h2_payment_debit
                << ",credit:" << maximum_h2_payment_credit;
        } else {
            std::cout << "none";
        }
        std::cout
            << " result="
            << (
                negative_total == 0U && negative_residue_tails == 0U
                    ? "PASS_FIXED_OFFSET_DISCOVERY"
                    : "FAIL_FIXED_OFFSET_CANDIDATE"
            ) << '\n';
        return negative_total == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
