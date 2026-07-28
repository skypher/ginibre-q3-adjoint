#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Matrix = std::vector<std::vector<int>>;

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

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error("usage: MAXIMUM_HALF_LEVEL");
        }
        const int maximum_level = parse_positive(argv[1]);
        std::uint64_t parameters = 0U;
        std::uint64_t shell_rows = 0U;
        std::uint64_t entries = 0U;
        std::uint64_t negative_k1 = 0U;
        std::uint64_t negative_k2 = 0U;
        std::uint64_t negative_terminal_k1 = 0U;
        std::uint64_t negative_terminal_k2 = 0U;
        std::uint64_t negative_shifted_k2 = 0U;
        std::uint64_t negative_terminal_shifted_k2 = 0U;
        std::uint64_t negative_smoothed_terminal_k2 = 0U;
        std::uint64_t negative_smoothed_single_column_k2 = 0U;
        std::uint64_t negative_shifted_single_column_k2 = 0U;
        std::uint64_t negative_nearest_shifted_single_column_k2 = 0U;
        std::uint64_t smoothed_single_column_sign_recrossings = 0U;
        std::uint64_t smoothed_single_column_rows = 0U;
        std::uint64_t negative_smoothed_adjacent_column_pairs = 0U;
        std::uint64_t negative_smoothed_innermost_columns = 0U;
        std::uint64_t negative_smoothed_boundary_corrections = 0U;
        int maximum_smoothed_single_column_sign_changes = 0;
        bool printed_k1 = false;
        bool printed_k2 = false;
        bool printed_shifted_k2 = false;
        bool printed_smoothed_k2 = false;
        bool printed_smoothed_single_column_k2 = false;
        bool printed_shifted_single_column_k2 = false;
        bool printed_smoothed_adjacent_column_pair = false;
        bool printed_smoothed_innermost_column = false;
        bool printed_smoothed_boundary_correction = false;

        for (int level = 3; level <= maximum_level; ++level) {
            for (int label = 1; 2 * label < level; ++label) {
                ++parameters;
                const int paired = (level + 1) / 2;
                const bool has_center = (level % 2) == 0;
                const int quotient_size = paired + (has_center ? 1 : 0);
                const int center = has_center ? paired : -1;
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
                    if (has_center) {
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
                }
                if (has_center) {
                    plus[static_cast<std::size_t>(center)]
                        [static_cast<std::size_t>(center)] = 1;
                }
                Matrix crossing = plus;
                for (int source = 0;
                     source < quotient_size;
                     ++source) {
                    for (int target = 0;
                         target < quotient_size;
                         ++target) {
                        crossing[static_cast<std::size_t>(source)]
                                [static_cast<std::size_t>(target)] -=
                            minus[static_cast<std::size_t>(source)]
                                 [static_cast<std::size_t>(target)];
                    }
                }

                std::vector<Integer> plus_row(
                    static_cast<std::size_t>(quotient_size)
                );
                std::vector<Integer> minus_row = plus_row;
                plus_row[0] = 1;
                minus_row[0] = 1;
                Integer plus_four = 0;
                Integer plus_five = 0;
                Integer minus_four = 0;
                Integer minus_five = 0;
                for (int power = 1; power <= 5; ++power) {
                    plus_row = multiply_row(plus_row, plus);
                    minus_row = multiply_row(minus_row, minus);
                    if (power == 4) {
                        plus_four = plus_row[0];
                        minus_four = minus_row[0];
                    } else if (power == 5) {
                        plus_five = plus_row[0];
                        minus_five = minus_row[0];
                    }
                }
                const Integer f4 = (plus_four + minus_four) / 2;
                const Integer f5 = (plus_five + minus_five) / 2;
                std::vector<std::vector<Integer>>
                    smoothed_single_column_values(
                        static_cast<std::size_t>(paired),
                        std::vector<Integer>(
                            static_cast<std::size_t>(quotient_size)
                        )
                    );

                for (int rho = 0; rho < paired; ++rho) {
                    ++shell_rows;
                    const int block_size = 2 * quotient_size;
                    Matrix block(
                        static_cast<std::size_t>(block_size),
                        std::vector<int>(
                            static_cast<std::size_t>(block_size)
                        )
                    );
                    for (int source = 0;
                         source < quotient_size;
                         ++source) {
                        for (int target = 0;
                             target < quotient_size;
                             ++target) {
                            block[static_cast<std::size_t>(source)]
                                 [static_cast<std::size_t>(target)] =
                                plus[static_cast<std::size_t>(source)]
                                    [static_cast<std::size_t>(target)];
                            block[static_cast<std::size_t>(
                                quotient_size + source
                            )][static_cast<std::size_t>(
                                quotient_size + target
                            )] =
                                minus[static_cast<std::size_t>(source)]
                                     [static_cast<std::size_t>(target)];
                            if (target >= rho && target < paired) {
                                block[static_cast<std::size_t>(source)]
                                     [static_cast<std::size_t>(
                                         quotient_size + target
                                     )] =
                                    crossing[
                                        static_cast<std::size_t>(source)
                                      ][static_cast<std::size_t>(target)];
                            }
                        }
                    }
                    std::vector<std::vector<Integer>> rows(
                        7U,
                        std::vector<Integer>(
                            static_cast<std::size_t>(block_size)
                        )
                    );
                    rows[0][0] = 1;
                    for (int power = 1; power <= 6; ++power) {
                        rows[static_cast<std::size_t>(power)] =
                            multiply_row(
                                rows[static_cast<std::size_t>(power - 1)],
                                block
                            );
                    }
                    for (int target = 0; target < block_size; ++target) {
                        ++entries;
                        const Integer k1 =
                            rows[3][static_cast<std::size_t>(target)]
                            - rows[2][static_cast<std::size_t>(target)];
                        const Integer k2 =
                            f4 * rows[5][static_cast<std::size_t>(target)]
                            - f5 * rows[4][static_cast<std::size_t>(target)];
                        const Integer shifted_k2 =
                            f4 * rows[6][static_cast<std::size_t>(target)]
                            - f5 * rows[5][static_cast<std::size_t>(target)];
                        const bool terminal =
                            target >= quotient_size;
                        if (k1 < 0) {
                            ++negative_k1;
                            if (terminal) {
                                ++negative_terminal_k1;
                            }
                            if (!printed_k1) {
                                printed_k1 = true;
                                std::cout
                                    << "FIRST_NEGATIVE_BLOCK_K1"
                                    << " K=" << level
                                    << " Q=" << label
                                    << " rho=" << rho
                                    << " target=" << target
                                    << " terminal=" << terminal
                                    << " value=" << k1 << '\n';
                            }
                        }
                        if (k2 < 0) {
                            ++negative_k2;
                            if (terminal) {
                                ++negative_terminal_k2;
                            }
                            if (!printed_k2) {
                                printed_k2 = true;
                                std::cout
                                    << "FIRST_NEGATIVE_BLOCK_K2"
                                    << " K=" << level
                                    << " Q=" << label
                                    << " rho=" << rho
                                    << " target=" << target
                                    << " terminal=" << terminal
                                    << " value=" << k2 << '\n';
                            }
                        }
                        if (shifted_k2 < 0) {
                            ++negative_shifted_k2;
                            if (terminal) {
                                ++negative_terminal_shifted_k2;
                            }
                            if (!printed_shifted_k2) {
                                printed_shifted_k2 = true;
                                std::cout
                                    << "FIRST_NEGATIVE_BLOCK_SHIFTED_K2"
                                    << " K=" << level
                                    << " Q=" << label
                                    << " rho=" << rho
                                    << " target=" << target
                                    << " terminal=" << terminal
                                    << " value=" << shifted_k2 << '\n';
                            }
                        }
                    }
                    for (int target = 0;
                         target < quotient_size;
                         ++target) {
                        Integer smoothed = 0;
                        for (int source = 0;
                             source < quotient_size;
                             ++source) {
                            const std::size_t block_source =
                                static_cast<std::size_t>(
                                    quotient_size + source
                                );
                            const Integer unshifted =
                                f4 * rows[5][block_source]
                                - f5 * rows[4][block_source];
                            smoothed += unshifted
                                * minus[
                                    static_cast<std::size_t>(source)
                                  ][static_cast<std::size_t>(target)];
                        }
                        Integer folded_reserve = 0;
                        for (int source = rho;
                             source < paired;
                             ++source) {
                            const Integer plus_kernel =
                                f4 * plus_row[
                                    static_cast<std::size_t>(source)
                                ]
                                - f5 * rows[4][
                                    static_cast<std::size_t>(source)
                                ];
                            folded_reserve += plus_kernel
                                * minus[
                                    static_cast<std::size_t>(source)
                                  ][static_cast<std::size_t>(target)];
                        }
                        const Integer boundary_correction =
                            smoothed - folded_reserve;
                        if (boundary_correction < 0) {
                            ++negative_smoothed_boundary_corrections;
                            if (
                                rho > 0
                                && !printed_smoothed_boundary_correction
                            ) {
                                printed_smoothed_boundary_correction = true;
                                std::cout
                                    << "FIRST_NEGATIVE_SMOOTHED_BOUNDARY_CORRECTION"
                                    << " K=" << level
                                    << " Q=" << label
                                    << " rho=" << rho
                                    << " target=" << target
                                    << " value=" << boundary_correction
                                    << " reserve=" << folded_reserve
                                    << " total=" << smoothed << '\n';
                            }
                        }
                        if (smoothed < 0) {
                            ++negative_smoothed_terminal_k2;
                            if (!printed_smoothed_k2) {
                                printed_smoothed_k2 = true;
                                std::cout
                                    << "FIRST_NEGATIVE_SMOOTHED_BLOCK_K2"
                                    << " K=" << level
                                    << " Q=" << label
                                    << " rho=" << rho
                                    << " target=" << target
                                    << " value=" << smoothed << '\n';
                            }
                        }
                    }

                    Matrix single_column_block = block;
                    for (int source = 0;
                         source < quotient_size;
                         ++source) {
                        for (int target = 0;
                             target < paired;
                             ++target) {
                            if (target != rho) {
                                single_column_block[
                                    static_cast<std::size_t>(source)
                                ][static_cast<std::size_t>(
                                    quotient_size + target
                                )] = 0;
                            }
                        }
                    }
                    std::vector<std::vector<Integer>> single_column_rows(
                        7U,
                        std::vector<Integer>(
                            static_cast<std::size_t>(block_size)
                        )
                    );
                    single_column_rows[0][0] = 1;
                    for (int power = 1; power <= 6; ++power) {
                        single_column_rows[
                            static_cast<std::size_t>(power)
                        ] = multiply_row(
                            single_column_rows[
                                static_cast<std::size_t>(power - 1)
                            ],
                            single_column_block
                        );
                    }
                    for (int target = 0;
                         target < quotient_size;
                         ++target) {
                        const std::size_t terminal_target =
                            static_cast<std::size_t>(
                                quotient_size + target
                            );
                        const Integer shifted =
                            f4 * single_column_rows[6][terminal_target]
                            - f5 * single_column_rows[5][terminal_target];
                        if (shifted < 0) {
                            ++negative_shifted_single_column_k2;
                            if (level == 2 * label + 1) {
                                ++negative_nearest_shifted_single_column_k2;
                            }
                            if (!printed_shifted_single_column_k2) {
                                printed_shifted_single_column_k2 = true;
                                std::cout
                                    << "FIRST_NEGATIVE_SHIFTED_SINGLE_COLUMN_K2"
                                    << " K=" << level
                                    << " Q=" << label
                                    << " rho=" << rho
                                    << " target=" << target
                                    << " value=" << shifted << '\n';
                            }
                        }
                        Integer smoothed = 0;
                        for (int source = 0;
                             source < quotient_size;
                             ++source) {
                            const std::size_t block_source =
                                static_cast<std::size_t>(
                                    quotient_size + source
                                );
                            const Integer unshifted =
                                f4 * single_column_rows[5][block_source]
                                - f5
                                    * single_column_rows[4][block_source];
                            smoothed += unshifted
                                * minus[
                                    static_cast<std::size_t>(source)
                                  ][static_cast<std::size_t>(target)];
                        }
                        smoothed_single_column_values[
                            static_cast<std::size_t>(rho)
                        ][static_cast<std::size_t>(target)] = smoothed;
                        if (smoothed < 0) {
                            ++negative_smoothed_single_column_k2;
                            if (!printed_smoothed_single_column_k2) {
                                printed_smoothed_single_column_k2 = true;
                                std::cout
                                    << "FIRST_NEGATIVE_SMOOTHED_SINGLE_COLUMN_K2"
                                    << " K=" << level
                                    << " Q=" << label
                                    << " rho=" << rho
                                    << " target=" << target
                                    << " value=" << smoothed << '\n';
                            }
                        }
                    }
                }
                for (int target = 0;
                     target < quotient_size;
                     ++target) {
                    ++smoothed_single_column_rows;
                    int previous_sign = 0;
                    int sign_changes = 0;
                    for (int rho = paired - 1; rho >= 0; --rho) {
                        const Integer& value =
                            smoothed_single_column_values[
                                static_cast<std::size_t>(rho)
                            ][static_cast<std::size_t>(target)];
                        const int sign = value > 0 ? 1 : (value < 0 ? -1 : 0);
                        if (sign == 0) {
                            continue;
                        }
                        if (previous_sign != 0 && sign != previous_sign) {
                            ++sign_changes;
                        }
                        previous_sign = sign;
                    }
                    if (sign_changes > maximum_smoothed_single_column_sign_changes) {
                        maximum_smoothed_single_column_sign_changes =
                            sign_changes;
                    }
                    if (sign_changes > 1) {
                        ++smoothed_single_column_sign_recrossings;
                    }
                    const Integer& innermost =
                        smoothed_single_column_values[
                            static_cast<std::size_t>(paired - 1)
                        ][static_cast<std::size_t>(target)];
                    if (innermost < 0) {
                        ++negative_smoothed_innermost_columns;
                        if (!printed_smoothed_innermost_column) {
                            printed_smoothed_innermost_column = true;
                            std::cout
                                << "FIRST_NEGATIVE_SMOOTHED_INNERMOST_COLUMN"
                                << " K=" << level
                                << " Q=" << label
                                << " target=" << target
                                << " value=" << innermost << '\n';
                        }
                    }
                    for (int rho = 0; rho + 1 < paired; ++rho) {
                        const Integer pair =
                            smoothed_single_column_values[
                                static_cast<std::size_t>(rho)
                            ][static_cast<std::size_t>(target)]
                            + smoothed_single_column_values[
                                static_cast<std::size_t>(rho + 1)
                            ][static_cast<std::size_t>(target)];
                        if (pair < 0) {
                            ++negative_smoothed_adjacent_column_pairs;
                            if (!printed_smoothed_adjacent_column_pair) {
                                printed_smoothed_adjacent_column_pair = true;
                                std::cout
                                    << "FIRST_NEGATIVE_SMOOTHED_ADJACENT_PAIR"
                                    << " K=" << level
                                    << " Q=" << label
                                    << " rho=" << rho
                                    << " target=" << target
                                    << " value=" << pair << '\n';
                            }
                        }
                    }
                }
            }
        }
        std::cout
            << "SU2_SHELL_BLOCK_KERNEL"
            << " maximum_half_level=" << maximum_level
            << " parameters=" << parameters
            << " shell_rows=" << shell_rows
            << " entries=" << entries
            << " negative_k1=" << negative_k1
            << " negative_terminal_k1=" << negative_terminal_k1
            << " negative_k2=" << negative_k2
            << " negative_terminal_k2=" << negative_terminal_k2
            << " negative_shifted_k2=" << negative_shifted_k2
            << " negative_terminal_shifted_k2="
            << negative_terminal_shifted_k2
            << " negative_smoothed_terminal_k2="
            << negative_smoothed_terminal_k2
            << " negative_smoothed_single_column_k2="
            << negative_smoothed_single_column_k2
            << " negative_shifted_single_column_k2="
            << negative_shifted_single_column_k2
            << " negative_nearest_shifted_single_column_k2="
            << negative_nearest_shifted_single_column_k2
            << " smoothed_single_column_rows="
            << smoothed_single_column_rows
            << " smoothed_single_column_sign_recrossings="
            << smoothed_single_column_sign_recrossings
            << " maximum_smoothed_single_column_sign_changes="
            << maximum_smoothed_single_column_sign_changes
            << " negative_smoothed_innermost_columns="
            << negative_smoothed_innermost_columns
            << " negative_smoothed_adjacent_column_pairs="
            << negative_smoothed_adjacent_column_pairs
            << " negative_smoothed_boundary_corrections="
            << negative_smoothed_boundary_corrections
            << " result="
            << (
                negative_k1 == 0U
                    && negative_shifted_k2 == 0U
                    && negative_smoothed_terminal_k2 == 0U
                    ? "PASS_SMOOTHED_BLOCK_KERNEL_DISCOVERY"
                    : "FAIL_SMOOTHED_BLOCK_KERNEL_CANDIDATE"
            )
            << '\n';
        return negative_k1 == 0U
                && negative_shifted_k2 == 0U
                && negative_smoothed_terminal_k2 == 0U
            ? 0
            : 1;
    } catch (const std::exception& exception) {
        std::cerr << "error: " << exception.what() << '\n';
        return 2;
    }
}
