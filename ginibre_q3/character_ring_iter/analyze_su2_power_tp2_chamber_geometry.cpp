#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>
#include <z3++.h>

#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using Pattern = std::array<int, 6>;

struct Constraint {
    std::array<long long, 3> normal;
    long long lower;
};

struct AffineCoordinate {
    std::array<long long, 3> coefficients;
    long long constant;
};

std::string show_pattern(const Pattern& pattern);

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

void append_chamber_constraints(
    std::vector<Constraint>& constraints,
    int power,
    int category,
    const AffineCoordinate& coordinate
) {
    if (category == power) {
        Constraint support;
        support.normal = coordinate.coefficients;
        support.normal[0] -= power;
        support.lower = 1 - coordinate.constant;
        constraints.push_back(support);
        return;
    }

    Constraint lower;
    lower.normal = coordinate.coefficients;
    lower.normal[0] += power - 2 * category;
    lower.lower = category - coordinate.constant;
    constraints.push_back(lower);

    Constraint upper;
    for (int index = 0; index < 3; ++index) {
        upper.normal[static_cast<std::size_t>(index)] =
            -coordinate.coefficients[static_cast<std::size_t>(index)];
    }
    upper.normal[0] += 2 * category + 2 - power;
    upper.lower = coordinate.constant - category;
    constraints.push_back(upper);

    Constraint support_upper;
    for (int index = 0; index < 3; ++index) {
        support_upper.normal[static_cast<std::size_t>(index)] =
            -coordinate.coefficients[static_cast<std::size_t>(index)];
    }
    support_upper.normal[0] += power;
    support_upper.lower = coordinate.constant;
    constraints.push_back(support_upper);
}

std::vector<Constraint> make_constraints(
    int power,
    const Pattern& pattern
) {
    std::vector<Constraint> constraints{
        Constraint{{1, 0, 0}, 1},
        Constraint{{0, 1, 0}, 1},
        Constraint{{0, 0, 1}, 0},
        Constraint{{power, -1, -2}, 1}
    };
    const AffineCoordinate x_minus{{0, 1, 0}, -1};
    const AffineCoordinate x{{0, 1, 0}, 0};
    const AffineCoordinate x_plus{{0, 1, 0}, 1};
    const AffineCoordinate y{{0, 1, 2}, 1};
    const AffineCoordinate y_plus{{0, 1, 2}, 2};
    const AffineCoordinate y_plus_two{{0, 1, 2}, 3};
    append_chamber_constraints(
        constraints,
        power,
        pattern[0],
        x_minus
    );
    append_chamber_constraints(
        constraints,
        power,
        pattern[1],
        x
    );
    append_chamber_constraints(
        constraints,
        power,
        pattern[2],
        x_plus
    );
    append_chamber_constraints(
        constraints,
        power,
        pattern[3],
        y
    );
    append_chamber_constraints(
        constraints,
        power,
        pattern[4],
        y_plus
    );
    append_chamber_constraints(
        constraints,
        power,
        pattern[5],
        y_plus_two
    );
    return constraints;
}

Integer determinant(
    const std::array<std::array<long long, 3>, 3>& matrix
) {
    return
        Integer{matrix[0][0]}
            * (
                Integer{matrix[1][1]} * matrix[2][2]
                - Integer{matrix[1][2]} * matrix[2][1]
            )
        - Integer{matrix[0][1]}
            * (
                Integer{matrix[1][0]} * matrix[2][2]
                - Integer{matrix[1][2]} * matrix[2][0]
            )
        + Integer{matrix[0][2]}
            * (
                Integer{matrix[1][0]} * matrix[2][1]
                - Integer{matrix[1][1]} * matrix[2][0]
            );
}

bool feasible_point(
    const std::vector<Constraint>& constraints,
    const std::array<Rational, 3>& point
) {
    for (const Constraint& constraint : constraints) {
        Rational value = 0;
        for (int index = 0; index < 3; ++index) {
            value +=
                constraint.normal[static_cast<std::size_t>(index)]
                * point[static_cast<std::size_t>(index)];
        }
        if (value < constraint.lower) {
            return false;
        }
    }
    return true;
}

std::string rational_key(const Rational& value) {
    return
        value.numerator().convert_to<std::string>()
        + "/" + value.denominator().convert_to<std::string>();
}

std::string point_key(const std::array<Rational, 3>& point) {
    return
        rational_key(point[0]) + ","
        + rational_key(point[1]) + ","
        + rational_key(point[2]);
}

std::vector<std::array<Rational, 3>> vertices(
    const std::vector<Constraint>& constraints
) {
    std::vector<std::array<Rational, 3>> result;
    std::set<std::string> seen;
    for (std::size_t first = 0; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1;
             second < constraints.size();
             ++second) {
            for (std::size_t third = second + 1;
                 third < constraints.size();
                 ++third) {
                const std::array<std::size_t, 3> selected{
                    first,
                    second,
                    third
                };
                std::array<std::array<long long, 3>, 3> matrix{};
                for (int row = 0; row < 3; ++row) {
                    matrix[static_cast<std::size_t>(row)] =
                        constraints[selected[
                            static_cast<std::size_t>(row)
                        ]].normal;
                }
                const Integer denominator = determinant(matrix);
                if (denominator == 0) {
                    continue;
                }
                std::array<Rational, 3> point{
                    Rational{0},
                    Rational{0},
                    Rational{0}
                };
                for (int column = 0; column < 3; ++column) {
                    auto replaced = matrix;
                    for (int row = 0; row < 3; ++row) {
                        replaced[static_cast<std::size_t>(row)]
                            [static_cast<std::size_t>(column)] =
                            constraints[selected[
                                static_cast<std::size_t>(row)
                            ]].lower;
                    }
                    Integer numerator = determinant(replaced);
                    Integer positive_denominator = denominator;
                    if (positive_denominator < 0) {
                        numerator = -numerator;
                        positive_denominator =
                            -positive_denominator;
                    }
                    point[static_cast<std::size_t>(column)] =
                        Rational{numerator, positive_denominator};
                }
                if (!feasible_point(constraints, point)) {
                    continue;
                }
                const std::string key = point_key(point);
                if (seen.insert(key).second) {
                    result.push_back(point);
                }
            }
        }
    }
    return result;
}

long long dot(
    const std::array<long long, 3>& left,
    const std::array<long long, 3>& right
) {
    return
        left[0] * right[0]
        + left[1] * right[1]
        + left[2] * right[2];
}

std::array<long long, 3> cross(
    const std::array<long long, 3>& left,
    const std::array<long long, 3>& right
) {
    return {
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0]
    };
}

std::array<long long, 3> normalize_ray(
    std::array<long long, 3> ray
) {
    long long divisor = 0;
    for (const long long coordinate : ray) {
        divisor = std::gcd(divisor, std::llabs(coordinate));
    }
    if (divisor > 1) {
        for (long long& coordinate : ray) {
            coordinate /= divisor;
        }
    }
    return ray;
}

std::vector<std::array<long long, 3>> rays(
    const std::vector<Constraint>& constraints
) {
    std::vector<std::array<long long, 3>> result;
    std::set<std::array<long long, 3>> seen;
    for (std::size_t first = 0; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1;
             second < constraints.size();
             ++second) {
            std::array<long long, 3> ray = cross(
                constraints[first].normal,
                constraints[second].normal
            );
            if (ray == std::array<long long, 3>{0, 0, 0}) {
                continue;
            }
            bool positive = true;
            bool negative = true;
            for (const Constraint& constraint : constraints) {
                const long long product =
                    dot(constraint.normal, ray);
                positive = positive && product >= 0;
                negative = negative && product <= 0;
            }
            if (!positive && !negative) {
                continue;
            }
            if (negative) {
                for (long long& coordinate : ray) {
                    coordinate = -coordinate;
                }
            }
            ray = normalize_ray(ray);
            if (seen.insert(ray).second) {
                result.push_back(ray);
            }
        }
    }
    return result;
}

bool integer_feasible(
    int power,
    const Pattern& pattern,
    const std::vector<Constraint>& constraints
) {
    z3::context context;
    z3::solver stored_solver(context, "QF_LIA");
    z3::solver direct_solver(context, "QF_LIA");
    const z3::expr q = context.int_const("q");
    const z3::expr x = context.int_const("x");
    const z3::expr a = context.int_const("a");
    const z3::expr y = x + 2 * a + 1;
    const std::array<z3::expr, 3> variables{q, x, a};
    for (const Constraint& constraint : constraints) {
        z3::expr expression = context.int_val(0);
        for (int index = 0; index < 3; ++index) {
            expression =
                expression
                + static_cast<int>(
                    constraint.normal[
                        static_cast<std::size_t>(index)
                    ]
                )
                    * variables[static_cast<std::size_t>(index)];
        }
        stored_solver.add(
            expression >= static_cast<int>(constraint.lower)
        );
    }
    direct_solver.add(q >= 1);
    direct_solver.add(x >= 1);
    direct_solver.add(a >= 0);
    direct_solver.add(y <= power * q);
    const std::array<z3::expr, 6> coordinates{
        x - 1,
        x,
        x + 1,
        y,
        y + 1,
        y + 2
    };
    for (int index = 0; index < 6; ++index) {
        const int category =
            pattern[static_cast<std::size_t>(index)];
        const z3::expr& coordinate =
            coordinates[static_cast<std::size_t>(index)];
        if (category == power) {
            direct_solver.add(coordinate >= power * q + 1);
        } else {
            const z3::expr width = 2 * q + 1;
            const z3::expr centered_degree =
                power * q + coordinate;
            direct_solver.add(category * width <= centered_degree);
            direct_solver.add(
                centered_degree
                    <= (category + 1) * width - 1
            );
            direct_solver.add(coordinate <= power * q);
        }
    }
    const z3::check_result stored = stored_solver.check();
    const z3::check_result direct = direct_solver.check();
    if (stored != direct) {
        throw std::runtime_error(
            "stored/direct chamber feasibility mismatch for "
            + show_pattern(pattern)
        );
    }
    return direct == z3::sat;
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
                "usage: analyze_su2_power_tp2_chamber_geometry power"
            );
        }
        const int power = parse_power(argv[1]);
        std::size_t patterns = 0;
        std::size_t feasible_patterns = 0;
        std::map<std::pair<std::size_t, std::size_t>, std::size_t>
            geometry_counts;
        std::string maximum_pattern;
        std::size_t maximum_vertices = 0;
        std::size_t maximum_rays = 0;

        for (int x_minus = 0; x_minus < power; ++x_minus) {
            for (
                const int x_category :
                next_categories(x_minus, power)
            ) {
                if (x_category == power) {
                    continue;
                }
                for (
                    const int x_plus :
                    next_categories(x_category, power)
                ) {
                    for (int y_category = 0;
                         y_category < power;
                         ++y_category) {
                        for (
                            const int y_plus :
                            next_categories(y_category, power)
                        ) {
                            for (
                                const int y_plus_two :
                                next_categories(y_plus, power)
                            ) {
                                ++patterns;
                                const Pattern pattern{
                                    x_minus,
                                    x_category,
                                    x_plus,
                                    y_category,
                                    y_plus,
                                    y_plus_two
                                };
                                const std::vector<Constraint>
                                    constraints =
                                        make_constraints(
                                            power,
                                            pattern
                                        );
                                if (
                                    !integer_feasible(
                                        power,
                                        pattern,
                                        constraints
                                    )
                                ) {
                                    continue;
                                }
                                ++feasible_patterns;
                                const auto pattern_vertices =
                                    vertices(constraints);
                                const auto pattern_rays =
                                    rays(constraints);
                                ++geometry_counts[{
                                    pattern_vertices.size(),
                                    pattern_rays.size()
                                }];
                                if (
                                    pattern_vertices.size()
                                        + pattern_rays.size()
                                    > maximum_vertices
                                        + maximum_rays
                                ) {
                                    maximum_vertices =
                                        pattern_vertices.size();
                                    maximum_rays =
                                        pattern_rays.size();
                                    maximum_pattern =
                                        show_pattern(pattern);
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_POWER_TP2_CHAMBER_GEOMETRY"
            << " power=" << power
            << " patterns=" << patterns
            << " feasible_patterns=" << feasible_patterns
            << " geometry_counts={";
        bool first = true;
        for (const auto& [geometry, count] : geometry_counts) {
            if (!first) {
                std::cout << ',';
            }
            first = false;
            std::cout
                << geometry.first << 'v' << geometry.second
                << "r:" << count;
        }
        std::cout
            << "}"
            << " maximum_pattern=" << maximum_pattern
            << " maximum_vertices=" << maximum_vertices
            << " maximum_rays=" << maximum_rays
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
