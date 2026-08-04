#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Matrix = std::vector<std::vector<Integer>>;

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(
            std::string(name) + " must be a positive integer"
        );
    }
    return static_cast<int>(value);
}

bool fuses(int level, int first, int second, int output) {
    return std::abs(first - second) <= output
        && output <= std::min(
            first + second, 2 * level - first - second
        )
        && ((first + second + output) & 1) == 0;
}

Matrix fusion_matrix(int level, int label) {
    const int width = level + 1;
    Matrix result(
        static_cast<std::size_t>(width),
        std::vector<Integer>(static_cast<std::size_t>(width))
    );
    for (int row = 0; row <= level; ++row) {
        for (int column = 0; column <= level; ++column) {
            if (fuses(level, label, row, column)) {
                result[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(column)] = 1;
            }
        }
    }
    return result;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const std::size_t width = left.size();
    Matrix result(width, std::vector<Integer>(width));
    for (std::size_t row = 0U; row < width; ++row) {
        for (std::size_t middle = 0U; middle < width; ++middle) {
            if (left[row][middle] == 0) {
                continue;
            }
            for (std::size_t column = 0U; column < width; ++column) {
                if (right[middle][column] != 0) {
                    result[row][column]
                        += left[row][middle] * right[middle][column];
                }
            }
        }
    }
    return result;
}

bool in_square_support(int label, int coordinate) {
    return coordinate >= 0 && coordinate <= 2 * label
        && (coordinate & 1) == 0;
}

void apply_additive_compound(
    int level,
    int label,
    std::vector<Integer>& state
) {
    const int width = level + 1;
    const auto index = [width](int first, int second) {
        return static_cast<std::size_t>(first * width + second);
    };
    std::vector<Integer> next(
        static_cast<std::size_t>(width * width)
    );
    for (int first = 0; first <= level; ++first) {
        for (int second = 0; second <= level; ++second) {
            const Integer& coefficient = state[index(first, second)];
            if (coefficient == 0) {
                continue;
            }
            for (int output = 0; output <= level; ++output) {
                if (fuses(level, label, first, output)) {
                    next[index(output, second)] += coefficient;
                }
                if (fuses(level, label, second, output)) {
                    next[index(first, output)] += coefficient;
                }
            }
        }
    }
    state.swap(next);
}

void apply_square(
    int level,
    int label,
    std::vector<Integer>& state
) {
    apply_additive_compound(level, label, state);
    apply_additive_compound(level, label, state);
}

int verify_pair_power_two_formula(int maximum_level) {
    std::size_t parameter_rows = 0U;
    std::size_t coefficient_rows = 0U;
    for (int level = 6; level <= maximum_level; level += 2) {
        const int width = level + 1;
        const auto index = [width](int first, int second) {
            return static_cast<std::size_t>(first * width + second);
        };
        for (int label = 2; 2 * label < level; label += 2) {
            ++parameter_rows;
            const Matrix first = fusion_matrix(level, label);
            const Matrix second = multiply(first, first);
            const Matrix fourth = multiply(second, second);
            for (int source = 1; source <= level; ++source) {
                std::vector<Integer> state(
                    static_cast<std::size_t>(width * width)
                );
                state[index(source, 0)] = 1;
                state[index(0, source)] = -1;
                apply_square(level, label, state);
                apply_square(level, label, state);
                for (int target = 1; target <= level; ++target) {
                    if ((target & 1) != (source & 1)) {
                        continue;
                    }
                    ++coefficient_rows;
                    Integer formula
                        = fourth[static_cast<std::size_t>(target)]
                                 [static_cast<std::size_t>(source)]
                        + 6 * second[static_cast<std::size_t>(target)]
                                    [static_cast<std::size_t>(source)]
                        + 4 * first[static_cast<std::size_t>(target)]
                                   [static_cast<std::size_t>(source)];
                    if (target == source) {
                        formula += label + 1;
                    }
                    if (target == label) {
                        formula -= 4 * second[static_cast<std::size_t>(source)]
                                               [static_cast<std::size_t>(label)];
                    }
                    if (in_square_support(label, source)
                        && in_square_support(label, target)) {
                        formula -= 6;
                    }
                    if (source == label) {
                        formula -= 4 * second[static_cast<std::size_t>(target)]
                                               [static_cast<std::size_t>(label)];
                    }
                    const Integer& direct = state[index(target, 0)];
                    if (direct != formula || direct < 0) {
                        std::cout
                            << "SU2_TERMINAL_BOUNDARY_PAIR_POWER_TWO"
                            << " result=FAIL"
                            << " level=" << level
                            << " label=" << label
                            << " source=" << source
                            << " target=" << target
                            << " direct=" << direct
                            << " formula=" << formula
                            << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
        }
    }
    std::cout
        << "SU2_TERMINAL_BOUNDARY_PAIR_POWER_TWO"
        << " maximum_level=" << maximum_level
        << " parameter_rows=" << parameter_rows
        << " coefficient_rows=" << coefficient_rows
        << " result=PASS_IDENTITY\n";
    return EXIT_SUCCESS;
}

int verify_pair_power_three_formula(int maximum_level) {
    std::size_t parameter_rows = 0U;
    std::size_t coefficient_rows = 0U;
    for (int level = 6; level <= maximum_level; level += 2) {
        const int width = level + 1;
        const auto index = [width](int first, int second) {
            return static_cast<std::size_t>(first * width + second);
        };
        for (int label = 2; 2 * label < level; label += 2) {
            ++parameter_rows;
            const Matrix first = fusion_matrix(level, label);
            const Matrix second = multiply(first, first);
            const Matrix third = multiply(second, first);
            const Matrix fourth = multiply(second, second);
            const Matrix sixth = multiply(fourth, second);
            const Integer return_two =
                second[static_cast<std::size_t>(label)]
                      [static_cast<std::size_t>(label)];
            const Integer return_three =
                third[static_cast<std::size_t>(label)]
                     [static_cast<std::size_t>(label)];
            const Integer return_four =
                fourth[static_cast<std::size_t>(label)]
                      [static_cast<std::size_t>(label)];
            for (int source = 1; source <= level; ++source) {
                std::vector<Integer> state(
                    static_cast<std::size_t>(width * width)
                );
                state[index(source, 0)] = 1;
                state[index(0, source)] = -1;
                for (int step = 0; step < 3; ++step) {
                    apply_square(level, label, state);
                }
                for (int target = 1; target <= level; ++target) {
                    if ((target & 1) != (source & 1)) {
                        continue;
                    }
                    ++coefficient_rows;
                    Integer formula
                        = sixth[static_cast<std::size_t>(target)]
                                [static_cast<std::size_t>(source)]
                        + 15 * fourth[static_cast<std::size_t>(target)]
                                      [static_cast<std::size_t>(source)]
                        + 20 * third[static_cast<std::size_t>(target)]
                                      [static_cast<std::size_t>(source)]
                        + 15 * return_two
                            * second[static_cast<std::size_t>(target)]
                                    [static_cast<std::size_t>(source)]
                        + 6 * return_three
                            * first[static_cast<std::size_t>(target)]
                                   [static_cast<std::size_t>(source)];
                    if (target == source) {
                        formula += return_four;
                    }
                    if (target == label) {
                        formula -= 6
                            * fourth[static_cast<std::size_t>(source)]
                                    [static_cast<std::size_t>(label)];
                    }
                    if (in_square_support(label, target)) {
                        formula -= 15
                            * third[static_cast<std::size_t>(source)]
                                   [static_cast<std::size_t>(label)];
                    }
                    formula -= 20
                        * second[static_cast<std::size_t>(source)]
                                [static_cast<std::size_t>(label)]
                        * second[static_cast<std::size_t>(target)]
                                [static_cast<std::size_t>(label)];
                    if (in_square_support(label, source)) {
                        formula -= 15
                            * third[static_cast<std::size_t>(target)]
                                   [static_cast<std::size_t>(label)];
                    }
                    if (source == label) {
                        formula -= 6
                            * fourth[static_cast<std::size_t>(target)]
                                    [static_cast<std::size_t>(label)];
                    }
                    const Integer& direct = state[index(target, 0)];
                    if (direct != formula || direct < 0) {
                        std::cout
                            << "SU2_TERMINAL_BOUNDARY_PAIR_POWER_THREE"
                            << " result=FAIL"
                            << " level=" << level
                            << " label=" << label
                            << " source=" << source
                            << " target=" << target
                            << " direct=" << direct
                            << " formula=" << formula
                            << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
        }
    }
    std::cout
        << "SU2_TERMINAL_BOUNDARY_PAIR_POWER_THREE"
        << " maximum_level=" << maximum_level
        << " parameter_rows=" << parameter_rows
        << " coefficient_rows=" << coefficient_rows
        << " result=PASS_IDENTITY\n";
    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3
            && std::string(argv[1]) == "--verify-pair-power-two-formula") {
            const int maximum_level =
                parse_positive(argv[2], "maximum level");
            if (maximum_level < 6) {
                throw std::runtime_error("require maximum level>=6");
            }
            return verify_pair_power_two_formula(maximum_level);
        }
        if (argc == 3
            && std::string(argv[1]) == "--verify-pair-power-three-formula") {
            const int maximum_level =
                parse_positive(argv[2], "maximum level");
            if (maximum_level < 6) {
                throw std::runtime_error("require maximum level>=6");
            }
            return verify_pair_power_three_formula(maximum_level);
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_terminal_boundary_block "
                "MAXIMUM_LEVEL MAXIMUM_PAIR_POWER | "
                "--verify-pair-power-two-formula MAXIMUM_LEVEL | "
                "--verify-pair-power-three-formula MAXIMUM_LEVEL"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_power =
            parse_positive(argv[2], "maximum pair power");
        if (maximum_level < 6) {
            throw std::runtime_error("require maximum level>=6");
        }

        std::size_t parameter_rows = 0;
        std::size_t matrix_rows = 0;
        std::size_t coefficient_rows = 0;
        for (int level = 6; level <= maximum_level; level += 2) {
            const int width = level + 1;
            const auto index = [width](int first, int second) {
                return static_cast<std::size_t>(
                    first * width + second
                );
            };
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameter_rows;
                for (int source = 1; source <= level; ++source) {
                    std::vector<Integer> state(
                        static_cast<std::size_t>(width * width)
                    );
                    state[index(source, 0)] = 1;
                    state[index(0, source)] = -1;
                    for (int power = 1;
                         power <= maximum_power; ++power) {
                        apply_square(level, label, state);
                        ++matrix_rows;
                        for (int target = 1;
                             target <= level; ++target) {
                            if ((target & 1) != (source & 1)) {
                                continue;
                            }
                            ++coefficient_rows;
                            const Integer& coefficient =
                                state[index(target, 0)];
                            if (coefficient < 0) {
                                std::cout
                                    << "SU2_TERMINAL_BOUNDARY_BLOCK"
                                    << " result=FAIL"
                                    << " level=" << level
                                    << " label=" << label
                                    << " pair_power=" << power
                                    << " source=" << source
                                    << " target=" << target
                                    << " coefficient=" << coefficient
                                    << '\n';
                                return EXIT_FAILURE;
                            }
                            if (
                                state[index(0, target)]
                                != -coefficient
                            ) {
                                throw std::runtime_error(
                                    "antisymmetry invariant failed"
                                );
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_TERMINAL_BOUNDARY_BLOCK"
            << " maximum_level=" << maximum_level
            << " maximum_pair_power=" << maximum_power
            << " parameter_rows=" << parameter_rows
            << " matrix_rows=" << matrix_rows
            << " coefficient_rows=" << coefficient_rows
            << " negative_coefficients=0"
            << " result=PASS_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_BOUNDARY_BLOCK FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
