#include <boost/multiprecision/cpp_int.hpp>
#include <z3++.h>

#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Pattern = std::array<int, 6>;

int parse_power(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value < 2
        || value > 20
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error("power must lie in [2,20]");
    }
    return static_cast<int>(value);
}

Integer binomial_nonnegative(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    if (k > n - k) {
        k = n - k;
    }
    Integer result = 1;
    for (int index = 1; index <= k; ++index) {
        result *= n - k + index;
        result /= index;
    }
    return result;
}

Integer generalized_binomial(int upper, int lower) {
    Integer result = 1;
    for (int offset = 0; offset < lower; ++offset) {
        result *= upper - offset;
    }
    for (int divisor = 2; divisor <= lower; ++divisor) {
        result /= divisor;
    }
    return result;
}

Integer coefficient_formula(
    int power,
    int q,
    int coordinate,
    int quotient
) {
    const int width = 2 * q + 1;
    const int centered_degree = power * q + coordinate;
    Integer result = 0;
    for (int index = 0; index <= quotient; ++index) {
        Integer term =
            binomial_nonnegative(power, index)
            * generalized_binomial(
                centered_degree - index * width + power - 1,
                power - 1
            );
        result += index % 2 == 0 ? term : -term;
    }
    return result;
}

Integer chamber_value(
    int power,
    int q,
    int coordinate,
    int category
) {
    return
        category == power
        ? Integer{0}
        : coefficient_formula(
            power,
            q,
            coordinate,
            category
        );
}

Integer determinant(
    int power,
    const Pattern& pattern,
    int q_variable,
    int x_variable,
    int a
) {
    const int q = q_variable + 1;
    const int x = x_variable + 1;
    const int y = x + 2 * a + 1;
    const Integer vx_minus =
        chamber_value(power, q, x - 1, pattern[0]);
    const Integer vx =
        chamber_value(power, q, x, pattern[1]);
    const Integer vx_plus =
        chamber_value(power, q, x + 1, pattern[2]);
    const Integer vy =
        chamber_value(power, q, y, pattern[3]);
    const Integer vy_plus =
        chamber_value(power, q, y + 1, pattern[4]);
    const Integer vy_plus_two =
        chamber_value(power, q, y + 2, pattern[5]);
    return
        (vx - vy) * (vx - vy_plus_two)
        - (vx_plus - vy_plus) * (vx_minus - vy_plus);
}

std::vector<int> next_categories(int category, int power) {
    if (category == power) {
        return {power};
    }
    std::vector<int> result{category};
    if (category + 1 < power) {
        result.push_back(category + 1);
    }
    result.push_back(power);
    return result;
}

void add_chamber(
    z3::solver& solver,
    const z3::expr& q,
    const z3::expr& coordinate,
    int power,
    int category
) {
    if (category == power) {
        solver.add(coordinate >= power * q + 1);
        return;
    }
    const z3::expr width = 2 * q + 1;
    const z3::expr centered_degree = power * q + coordinate;
    solver.add(category * width <= centered_degree);
    solver.add(
        centered_degree <= (category + 1) * width - 1
    );
    solver.add(coordinate <= power * q);
}

bool integer_feasible(int power, const Pattern& pattern) {
    z3::context context;
    z3::solver solver(context, "QF_LIA");
    const z3::expr q = context.int_const("q");
    const z3::expr x = context.int_const("x");
    const z3::expr a = context.int_const("a");
    const z3::expr y = x + 2 * a + 1;
    solver.add(q >= 1);
    solver.add(x >= 1);
    solver.add(a >= 0);
    solver.add(y <= power * q);
    add_chamber(solver, q, x - 1, power, pattern[0]);
    add_chamber(solver, q, x, power, pattern[1]);
    add_chamber(solver, q, x + 1, power, pattern[2]);
    add_chamber(solver, q, y, power, pattern[3]);
    add_chamber(solver, q, y + 1, power, pattern[4]);
    add_chamber(solver, q, y + 2, power, pattern[5]);
    return solver.check() == z3::sat;
}

std::size_t flat_index(
    int first,
    int second,
    int third,
    int width
) {
    return static_cast<std::size_t>(
        (first * width + second) * width + third
    );
}

void transform_axis(
    std::vector<Integer>& grid,
    int width,
    int axis
) {
    for (int first = 0; first < width; ++first) {
        for (int second = 0; second < width; ++second) {
            std::vector<Integer> line(
                static_cast<std::size_t>(width)
            );
            for (int position = 0; position < width; ++position) {
                const int coordinates[3] = {
                    axis == 0 ? position : first,
                    axis == 1
                        ? position
                        : (axis == 0 ? first : second),
                    axis == 2 ? position : second
                };
                line[static_cast<std::size_t>(position)] =
                    grid[flat_index(
                        coordinates[0],
                        coordinates[1],
                        coordinates[2],
                        width
                    )];
            }
            for (int order = 0; order < width; ++order) {
                const int coordinates[3] = {
                    axis == 0 ? order : first,
                    axis == 1
                        ? order
                        : (axis == 0 ? first : second),
                    axis == 2 ? order : second
                };
                grid[flat_index(
                    coordinates[0],
                    coordinates[1],
                    coordinates[2],
                    width
                )] = line[0];
                for (int index = 0;
                     index + 1 < width - order;
                     ++index) {
                    line[static_cast<std::size_t>(index)] =
                        line[static_cast<std::size_t>(index + 1)]
                        - line[static_cast<std::size_t>(index)];
                }
            }
        }
    }
}

struct CertificateResult {
    std::size_t nonzero = 0;
    std::size_t negative = 0;
    std::size_t degree_violations = 0;
    std::array<int, 3> first_negative_order{0, 0, 0};
    Integer first_negative_value = 0;
};

CertificateResult certify_pattern(
    int power,
    const Pattern& pattern
) {
    const int degree_bound = 2 * (power - 1);
    const int width = degree_bound + 3;
    std::vector<Integer> grid(
        static_cast<std::size_t>(width * width * width)
    );
    for (int q = 0; q < width; ++q) {
        for (int x = 0; x < width; ++x) {
            for (int a = 0; a < width; ++a) {
                grid[flat_index(q, x, a, width)] =
                    determinant(power, pattern, q, x, a);
            }
        }
    }
    transform_axis(grid, width, 0);
    transform_axis(grid, width, 1);
    transform_axis(grid, width, 2);

    CertificateResult result;
    for (int q_order = 0; q_order < width; ++q_order) {
        for (int x_order = 0; x_order < width; ++x_order) {
            for (int a_order = 0; a_order < width; ++a_order) {
                const Integer& value = grid[flat_index(
                    q_order,
                    x_order,
                    a_order,
                    width
                )];
                if (value != 0) {
                    ++result.nonzero;
                }
                if (
                    q_order + x_order + a_order > degree_bound
                    && value != 0
                ) {
                    ++result.degree_violations;
                }
                if (value < 0) {
                    ++result.negative;
                    if (result.first_negative_value == 0) {
                        result.first_negative_order = {
                            q_order,
                            x_order,
                            a_order
                        };
                        result.first_negative_value = value;
                    }
                }
            }
        }
    }
    return result;
}

std::string show_pattern(const Pattern& pattern) {
    return
        "(" + std::to_string(pattern[0])
        + "," + std::to_string(pattern[1])
        + "," + std::to_string(pattern[2])
        + ";" + std::to_string(pattern[3])
        + "," + std::to_string(pattern[4])
        + "," + std::to_string(pattern[5]) + ")";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_power_tp2_newton_chambers power"
            );
        }
        const int power = parse_power(argv[1]);
        std::size_t patterns = 0;
        std::size_t feasible_patterns = 0;
        std::size_t globally_certified_patterns = 0;
        std::size_t negative_coefficient_patterns = 0;
        std::size_t degree_violations = 0;
        std::string first_failed_pattern;
        std::array<int, 3> first_failed_order{0, 0, 0};
        Integer first_failed_value = 0;

        for (int x_minus_category = 0;
             x_minus_category < power;
             ++x_minus_category) {
            for (
                const int x_category :
                next_categories(x_minus_category, power)
            ) {
                if (x_category == power) {
                    continue;
                }
                for (
                    const int x_plus_category :
                    next_categories(x_category, power)
                ) {
                    for (int y_category = 0;
                         y_category < power;
                         ++y_category) {
                        for (
                            const int y_plus_category :
                            next_categories(y_category, power)
                        ) {
                            for (
                                const int y_plus_two_category :
                                next_categories(
                                    y_plus_category,
                                    power
                                )
                            ) {
                                ++patterns;
                                const Pattern pattern{
                                    x_minus_category,
                                    x_category,
                                    x_plus_category,
                                    y_category,
                                    y_plus_category,
                                    y_plus_two_category
                                };
                                if (
                                    !integer_feasible(power, pattern)
                                ) {
                                    continue;
                                }
                                ++feasible_patterns;
                                const CertificateResult result =
                                    certify_pattern(power, pattern);
                                degree_violations +=
                                    result.degree_violations;
                                if (
                                    result.negative == 0
                                    && result.degree_violations == 0
                                ) {
                                    ++globally_certified_patterns;
                                } else {
                                    ++negative_coefficient_patterns;
                                    if (first_failed_pattern.empty()) {
                                        first_failed_pattern =
                                            show_pattern(pattern);
                                        first_failed_order =
                                            result
                                                .first_negative_order;
                                        first_failed_value =
                                            result
                                                .first_negative_value;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_POWER_TP2_NEWTON_CHAMBERS"
            << " power=" << power
            << " patterns=" << patterns
            << " feasible_patterns=" << feasible_patterns
            << " globally_certified_patterns="
            << globally_certified_patterns
            << " negative_coefficient_patterns="
            << negative_coefficient_patterns
            << " degree_violations=" << degree_violations
            << " first_failed_pattern="
            << (
                first_failed_pattern.empty()
                    ? "{}"
                    : first_failed_pattern
            )
            << " first_failed_order=("
            << first_failed_order[0] << ','
            << first_failed_order[1] << ','
            << first_failed_order[2] << ')'
            << " first_failed_value=" << first_failed_value
            << " result="
            << (
                globally_certified_patterns == feasible_patterns
                && degree_violations == 0
                    ? "PASS_GLOBAL_NEWTON_CERTIFICATE"
                    : "PARTIAL"
            )
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
