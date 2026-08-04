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
    const std::vector<Integer>& row,
    const Matrix& matrix
) {
    std::vector<Integer> result(matrix.size());
    for (std::size_t source = 0U; source < matrix.size(); ++source) {
        if (row[source] == 0) {
            continue;
        }
        for (std::size_t target = 0U; target < matrix.size(); ++target) {
            if (matrix[source][target] != 0) {
                result[target] += row[source] * matrix[source][target];
            }
        }
    }
    return result;
}

IntegerMatrix multiply_right(const IntegerMatrix& left, const Matrix& right) {
    IntegerMatrix result(
        left.size(),
        std::vector<Integer>(right.size())
    );
    for (std::size_t source = 0U; source < left.size(); ++source) {
        for (std::size_t middle = 0U; middle < right.size(); ++middle) {
            if (left[source][middle] == 0) {
                continue;
            }
            for (std::size_t target = 0U; target < right.size(); ++target) {
                if (right[middle][target] != 0) {
                    result[source][target] +=
                        left[source][middle] * right[middle][target];
                }
            }
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error("usage: MAXIMUM_K");
        }
        const int maximum_k = parse_positive(argv[1]);
        std::uint64_t parameters = 0U;
        std::uint64_t rows = 0U;
        std::uint64_t groups = 0U;
        std::uint64_t negative_groups = 0U;
        Integer minimum = 0;
        bool initialized_minimum = false;
        bool negative_found = false;
        int witness_level = 0;
        int witness_label = 0;
        int witness_rho = 0;
        int witness_target = 0;
        int witness_source = 0;
        Integer witness_value = 0;

        for (int level = 3; level <= maximum_k; ++level) {
            for (int label = 1; 2 * label < level; ++label) {
                ++parameters;
                const int paired = (level + 1) / 2;
                const bool has_center = level % 2 == 0;
                const int plus_size = paired + (has_center ? 1 : 0);
                const int center = has_center ? paired : -1;
                Matrix plus(
                    static_cast<std::size_t>(plus_size),
                    std::vector<int>(static_cast<std::size_t>(plus_size))
                );
                Matrix minus(
                    static_cast<std::size_t>(paired),
                    std::vector<int>(static_cast<std::size_t>(paired))
                );
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
                            [static_cast<std::size_t>(target)] = same + crossed;
                        minus[static_cast<std::size_t>(source)]
                             [static_cast<std::size_t>(target)] = same - crossed;
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
                Matrix crossing(
                    static_cast<std::size_t>(plus_size),
                    std::vector<int>(static_cast<std::size_t>(paired))
                );
                for (int source = 0; source < plus_size; ++source) {
                    for (int target = 0; target < paired; ++target) {
                        crossing[static_cast<std::size_t>(source)]
                                [static_cast<std::size_t>(target)] =
                            plus[static_cast<std::size_t>(source)]
                                [static_cast<std::size_t>(target)]
                            - (source < paired
                                ? minus[static_cast<std::size_t>(source)]
                                       [static_cast<std::size_t>(target)]
                                : 0);
                    }
                }

                std::vector<std::vector<Integer>> plus_rows(
                    6U,
                    std::vector<Integer>(static_cast<std::size_t>(plus_size))
                );
                plus_rows[0][0] = 1;
                for (int power = 1; power <= 5; ++power) {
                    plus_rows[static_cast<std::size_t>(power)] = multiply_row(
                        plus_rows[static_cast<std::size_t>(power - 1)],
                        plus
                    );
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
                const Integer f4 = (plus_rows[4][0] + minus_four) / 2;
                const Integer f5 = (plus_rows[5][0] + minus_five) / 2;

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
                            minus_powers[static_cast<std::size_t>(power - 1)],
                            minus
                        );
                }

                for (int rho = 0; rho < paired; ++rho) {
                    for (int target = 0; target < paired; ++target) {
                        Integer direct = 0;
                        Integer swapped = 0;
                        for (int crossing_target = rho;
                             crossing_target < paired;
                             ++crossing_target) {
                            Integer d1 = 0;
                            Integer d2 = 0;
                            Integer d3 = 0;
                            Integer d4 = 0;
                            Integer d5 = 0;
                            for (int source = 0; source < plus_size; ++source) {
                                const int weight = crossing[
                                    static_cast<std::size_t>(source)
                                ][static_cast<std::size_t>(crossing_target)];
                                d1 += plus_rows[1][static_cast<std::size_t>(source)] * weight;
                                d2 += plus_rows[2][static_cast<std::size_t>(source)] * weight;
                                d3 += plus_rows[3][static_cast<std::size_t>(source)] * weight;
                                d4 += plus_rows[4][static_cast<std::size_t>(source)] * weight;
                                d5 += plus_rows[5][static_cast<std::size_t>(source)] * weight;
                            }
                            direct += f4 * (
                                d1 * minus_powers[4][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)]
                                + d2 * minus_powers[3][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)]
                                + d3 * minus_powers[2][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)]
                                + d4 * minus_powers[1][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)]
                                + d5 * minus_powers[0][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)]
                            ) - f5 * (
                                d1 * minus_powers[3][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)]
                                + d2 * minus_powers[2][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)]
                                + d3 * minus_powers[1][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)]
                                + d4 * minus_powers[0][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)]
                            );
                        }
                        for (int source = 0; source < plus_size; ++source) {
                            Integer tail0 = 0;
                            Integer tail1 = 0;
                            Integer tail2 = 0;
                            Integer tail3 = 0;
                            Integer tail4 = 0;
                            for (int crossing_target = rho;
                                 crossing_target < paired;
                                 ++crossing_target) {
                                const int weight = crossing[
                                    static_cast<std::size_t>(source)
                                ][static_cast<std::size_t>(crossing_target)];
                                tail0 += weight * minus_powers[0][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)];
                                tail1 += weight * minus_powers[1][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)];
                                tail2 += weight * minus_powers[2][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)];
                                tail3 += weight * minus_powers[3][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)];
                                tail4 += weight * minus_powers[4][static_cast<std::size_t>(crossing_target)][static_cast<std::size_t>(target)];
                            }
                            const Integer group = f4 * (
                                plus_rows[1][static_cast<std::size_t>(source)] * tail4
                                + plus_rows[2][static_cast<std::size_t>(source)] * tail3
                                + plus_rows[3][static_cast<std::size_t>(source)] * tail2
                                + plus_rows[4][static_cast<std::size_t>(source)] * tail1
                                + plus_rows[5][static_cast<std::size_t>(source)] * tail0
                            ) - f5 * (
                                plus_rows[1][static_cast<std::size_t>(source)] * tail3
                                + plus_rows[2][static_cast<std::size_t>(source)] * tail2
                                + plus_rows[3][static_cast<std::size_t>(source)] * tail1
                                + plus_rows[4][static_cast<std::size_t>(source)] * tail0
                            );
                            swapped += group;
                            ++groups;
                            if (!initialized_minimum || group < minimum) {
                                initialized_minimum = true;
                                minimum = group;
                            }
                            if (group < 0) {
                                ++negative_groups;
                                if (!negative_found) {
                                    negative_found = true;
                                    witness_level = level;
                                    witness_label = label;
                                    witness_rho = rho;
                                    witness_target = target;
                                    witness_source = source;
                                    witness_value = group;
                                }
                            }
                        }
                        ++rows;
                        if (direct != swapped) {
                            throw std::runtime_error("tail-swap identity mismatch");
                        }
                        if (negative_found) {
                            std::cout
                                << "FIRST_NEGATIVE_TAIL_SWAPPED_GROUP"
                                << " K=" << witness_level
                                << " Q=" << witness_label
                                << " rho=" << witness_rho
                                << " target=" << witness_target
                                << " source=" << witness_source
                                << " value=" << witness_value
                                << " tail_swap=PASS_EXACT_IDENTITY\n";
                            return EXIT_SUCCESS;
                        }
                    }
                }
            }
        }
        std::cout
            << "SU2_SHELL_TAIL_SWAP"
            << " maximum_K=" << maximum_k
            << " parameters=" << parameters
            << " rows=" << rows
            << " groups=" << groups
            << " negative_groups=" << negative_groups
            << " minimum=" << minimum
            << " result=PASS_NO_NEGATIVE_GROUPS_BOUNDED_DIAGNOSTIC\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
