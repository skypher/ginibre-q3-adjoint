#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Vector = std::vector<Integer>;
using Matrix = std::vector<Vector>;

int parse_positive(const char* text, const char* name) {
    const std::string value{text};
    std::size_t used = 0U;
    const long long parsed = std::stoll(value, &used);
    if (used != value.size() || parsed <= 0LL) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const std::size_t dimension = left.size();
    Matrix result(dimension, Vector(dimension, 0));
    for (std::size_t row = 0U; row < dimension; ++row) {
        for (std::size_t middle = 0U; middle < dimension; ++middle) {
            if (left[row][middle] == 0) {
                continue;
            }
            for (std::size_t column = 0U; column < dimension; ++column) {
                result[row][column] += left[row][middle]
                    * right[middle][column];
            }
        }
    }
    return result;
}

Matrix fusion_matrix(int level, int factor) {
    Matrix result(
        static_cast<std::size_t>(level + 1),
        Vector(static_cast<std::size_t>(level + 1), 0)
    );
    for (int source = 0; source <= level; ++source) {
        const int lower = std::abs(source - factor);
        const int upper = std::min(
            source + factor,
            2 * level - source - factor
        );
        for (int target = lower; target <= upper; target += 2) {
            result[static_cast<std::size_t>(target)]
                  [static_cast<std::size_t>(source)] = 1;
        }
    }
    return result;
}

Integer minor(
    const Matrix& matrix,
    int first_row,
    int second_row,
    int first_column,
    int second_column
) {
    return matrix[static_cast<std::size_t>(first_row)]
                 [static_cast<std::size_t>(first_column)]
            * matrix[static_cast<std::size_t>(second_row)]
                    [static_cast<std::size_t>(second_column)]
        - matrix[static_cast<std::size_t>(first_row)]
                 [static_cast<std::size_t>(second_column)]
            * matrix[static_cast<std::size_t>(second_row)]
                    [static_cast<std::size_t>(first_column)];
}

Integer reflected_sum(
    const Matrix& power,
    int factor,
    int target,
    int rho,
    bool outer_band
) {
    const int level = static_cast<int>(power.size()) - 1;
    Integer result = 0;
    for (int first = 0; first <= level; ++first) {
        for (int second = first + 1; second <= level; ++second) {
            const bool in_window = first >= rho && second <= level - rho;
            if ((outer_band && in_window) || (!outer_band && !in_window)) {
                continue;
            }
            result += minor(power, 0, factor, first, second)
                * minor(power, 0, target, first, second);
        }
    }
    return result;
}

void replay_obstructions() {
    const Matrix fusion = fusion_matrix(6, 2);
    const Matrix square = multiply(fusion, fusion);
    const Matrix fourth = multiply(square, square);
    const Matrix cube = multiply(square, fusion);
    const Matrix sixth = multiply(cube, cube);
    const Integer left = minor(square, 0, 2, 2, 4);
    const Integer right = minor(square, 0, 4, 2, 4);
    const Integer half_power_two = minor(fourth, 0, 2, 0, 4);
    const Integer central_window = reflected_sum(square, 2, 4, 1, false);
    const Integer half_power_three = minor(sixth, 0, 2, 0, 4);
    const Integer outer_band = reflected_sum(cube, 2, 4, 1, true);
    if (left != -1 || right != 1 || half_power_two != 6
        || central_window != -1 || half_power_three != 35
        || outer_band != -1) {
        throw std::runtime_error("finite boundary obstruction replay mismatch");
    }
    std::cout
        << "SU2_FINITE_BOUNDARY_CAUCHY_BINET_OBSTRUCTIONS"
        << " channel=(-1,1;6)"
        << " central_window=-1"
        << " outer_band=(-1;35)"
        << " result=PASS_EXACT\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2
            && std::string(argv[1]) == "--replay-obstructions") {
            replay_obstructions();
            return EXIT_SUCCESS;
        }
        const bool window_mode = argc == 4
            && std::string(argv[1]) == "--reflected-window";
        const bool outer_mode = argc == 4
            && std::string(argv[1]) == "--reflected-outer-band";
        const bool cumulative_mode = window_mode || outer_mode;
        if ((!cumulative_mode && argc != 3)
            || (cumulative_mode && argc != 4)) {
            throw std::runtime_error(
                "usage: probe_su2_finite_boundary_cauchy_binet "
                "MAXIMUM_LEVEL MAXIMUM_HALF_POWER"
                " | --reflected-window MAXIMUM_LEVEL MAXIMUM_HALF_POWER"
                " | --reflected-outer-band MAXIMUM_LEVEL MAXIMUM_HALF_POWER"
                " | --replay-obstructions"
            );
        }
        const int offset = cumulative_mode ? 1 : 0;
        const int maximum_level = parse_positive(
            argv[1 + offset],
            "maximum level"
        );
        const int maximum_half_power = parse_positive(
            argv[2 + offset],
            "maximum half power"
        );
        std::uint64_t decompositions = 0U;
        std::uint64_t term_checks = 0U;
        for (int level = 1; level <= maximum_level; ++level) {
            for (int factor = 1; factor <= level; ++factor) {
                const Matrix fusion = fusion_matrix(level, factor);
                Matrix power = fusion;
                for (int half_power = 1;
                     half_power <= maximum_half_power;
                     ++half_power) {
                    const Matrix square = multiply(power, power);
                    for (int target = 1; target <= level; ++target) {
                        const Integer exact = minor(
                            square,
                            0,
                            factor,
                            0,
                            target
                        );
                        Integer reconstructed = 0;
                        Matrix contributions(
                            static_cast<std::size_t>(level + 1),
                            Vector(static_cast<std::size_t>(level + 1), 0)
                        );
                        for (int first = 0; first <= level; ++first) {
                            for (int second = first + 1;
                                 second <= level;
                                 ++second) {
                                const Integer left = minor(
                                    power,
                                    0,
                                    factor,
                                    first,
                                    second
                                );
                                const Integer right = minor(
                                    power,
                                    0,
                                    target,
                                    first,
                                    second
                                );
                                const Integer contribution = left * right;
                                contributions[static_cast<std::size_t>(first)]
                                             [static_cast<std::size_t>(second)] =
                                    contribution;
                                reconstructed += contribution;
                                ++term_checks;
                                if (!cumulative_mode && contribution < 0) {
                                    std::cout
                                        << "FINITE_BOUNDARY_CAUCHY_BINET_TERM_OBSTRUCTION"
                                        << " level=" << level
                                        << " factor=" << factor
                                        << " half_power=" << half_power
                                        << " target=" << target
                                        << " intermediate=(" << first << ','
                                        << second << ')'
                                        << " left=" << left
                                        << " right=" << right
                                        << " contribution=" << contribution
                                        << " anchored_minor=" << exact
                                        << '\n';
                                    return EXIT_FAILURE;
                                }
                            }
                        }
                        if (reconstructed != exact) {
                            throw std::runtime_error(
                                "Cauchy--Binet reconstruction mismatch"
                            );
                        }
                        if (cumulative_mode) {
                            for (int rho = 0; 2 * rho <= level; ++rho) {
                                Integer window = 0;
                                for (int first = 0;
                                     first <= level;
                                     ++first) {
                                    for (int second = first + 1;
                                         second <= level;
                                         ++second) {
                                        const bool in_window = first >= rho
                                            && second <= level - rho;
                                        if ((window_mode && !in_window)
                                            || (outer_mode && in_window)) {
                                            continue;
                                        }
                                        window += contributions[
                                            static_cast<std::size_t>(first)
                                        ][static_cast<std::size_t>(second)];
                                    }
                                }
                                if (window < 0) {
                                    std::cout
                                        << (window_mode
                                            ? "FINITE_BOUNDARY_CAUCHY_BINET_WINDOW_OBSTRUCTION"
                                            : "FINITE_BOUNDARY_CAUCHY_BINET_OUTER_BAND_OBSTRUCTION")
                                        << " level=" << level
                                        << " factor=" << factor
                                        << " half_power=" << half_power
                                        << " target=" << target
                                        << " rho=" << rho
                                        << " window=" << window
                                        << " anchored_minor=" << exact
                                        << '\n';
                                    return EXIT_FAILURE;
                                }
                            }
                        }
                        ++decompositions;
                    }
                    power = multiply(power, fusion);
                }
            }
        }
        std::cout
            << (window_mode
                ? "SU2_FINITE_BOUNDARY_CAUCHY_BINET_REFLECTED_WINDOW"
                : (outer_mode
                    ? "SU2_FINITE_BOUNDARY_CAUCHY_BINET_REFLECTED_OUTER_BAND"
                    : "SU2_FINITE_BOUNDARY_CAUCHY_BINET"))
            << " maximum_level=" << maximum_level
            << " maximum_half_power=" << maximum_half_power
            << " decompositions=" << decompositions
            << " term_checks=" << term_checks
            << " result="
            << (window_mode
                ? "PASS_WINDOW"
                : (outer_mode ? "PASS_OUTER_BAND" : "PASS_TERMWISE"))
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
