#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/*
 * Exact unbounded certificate for the adjacent minors of N_Q^m.
 *
 * On each bounded-composition chamber the minor is a polynomial
 * P(Q,x,a) of degree at most 2(m-1).  The substitution
 *
 *   r=1/Q,  xi=x/Q,  alpha=a/Q
 *
 * turns the chamber into a compact rational polytope and P into its
 * degree-2(m-1) homogenization.  The code enumerates every real chamber
 * having r>0, covers its closure by rational simplices, and computes all
 * simplex Bernstein coefficients exactly over cpp_int rationals.
 * Nonnegative coefficients prove the polynomial nonnegative on the whole
 * real chamber, hence on every required integer/parity point.
 *
 * A requested subdivision depth is only an exact refinement mechanism.
 * Reaching the depth with a negative coefficient returns failure; no cap,
 * timeout, sample, or numerical approximation is accepted as proof.
 */

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using Exponent = std::array<int, 3>;
using Point = std::array<Rational, 3>;
using Pattern = std::array<int, 6>;

struct Polynomial {
    std::map<Exponent, Rational> coefficients;
};

using ValueTable = std::array<std::vector<Polynomial>, 6>;

struct Constraint {
    std::array<long long, 3> normal;
    long long lower;
};

struct AffineCoordinate {
    std::array<long long, 3> coefficients;
    long long constant;
};

struct CompactConstraint {
    std::array<long long, 3> coefficients;
    long long constant;
};

struct SimplexResult {
    bool certified = false;
    std::size_t nodes = 0;
    std::size_t leaves = 0;
    std::size_t splits = 0;
    std::size_t negative_leaves = 0;
    int maximum_depth = 0;
    Rational minimum_coefficient = 0;
    bool minimum_initialized = false;
};

int parse_bounded_integer(
    const char* text,
    int minimum,
    int maximum,
    const std::string& name
) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value < minimum
        || value > maximum
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error(
            name + " must lie in ["
            + std::to_string(minimum) + ","
            + std::to_string(maximum) + "]"
        );
    }
    return static_cast<int>(value);
}

Polynomial constant_polynomial(const Rational& value) {
    Polynomial result;
    if (value != 0) {
        result.coefficients[{0, 0, 0}] = value;
    }
    return result;
}

Polynomial variable_polynomial(int coordinate) {
    Polynomial result;
    Exponent exponent{0, 0, 0};
    exponent[static_cast<std::size_t>(coordinate)] = 1;
    result.coefficients[exponent] = 1;
    return result;
}

void add_to(
    Polynomial& target,
    const Polynomial& source,
    const Rational& scale = Rational{1}
) {
    for (const auto& [exponent, coefficient] : source.coefficients) {
        Rational& destination = target.coefficients[exponent];
        destination += scale * coefficient;
        if (destination == 0) {
            target.coefficients.erase(exponent);
        }
    }
}

Polynomial add(
    const Polynomial& left,
    const Polynomial& right
) {
    Polynomial result = left;
    add_to(result, right);
    return result;
}

Polynomial subtract(
    const Polynomial& left,
    const Polynomial& right
) {
    Polynomial result = left;
    add_to(result, right, Rational{-1});
    return result;
}

Polynomial multiply(
    const Polynomial& left,
    const Polynomial& right
) {
    Polynomial result;
    for (
        const auto& [left_exponent, left_coefficient] :
        left.coefficients
    ) {
        for (
            const auto& [right_exponent, right_coefficient] :
            right.coefficients
        ) {
            Exponent exponent{};
            for (int coordinate = 0; coordinate < 3; ++coordinate) {
                exponent[static_cast<std::size_t>(coordinate)] =
                    left_exponent[
                        static_cast<std::size_t>(coordinate)
                    ]
                    + right_exponent[
                        static_cast<std::size_t>(coordinate)
                    ];
            }
            result.coefficients[exponent] +=
                left_coefficient * right_coefficient;
        }
    }
    return result;
}

Polynomial scale(
    const Polynomial& polynomial,
    const Rational& factor
) {
    Polynomial result;
    add_to(result, polynomial, factor);
    return result;
}

Integer factorial(int value) {
    Integer result = 1;
    for (int factor = 2; factor <= value; ++factor) {
        result *= factor;
    }
    return result;
}

Integer binomial_integer(int upper, int lower) {
    if (lower < 0 || lower > upper) {
        return 0;
    }
    lower = std::min(lower, upper - lower);
    Integer result = 1;
    for (int index = 1; index <= lower; ++index) {
        result *= upper - lower + index;
        result /= index;
    }
    return result;
}

Polynomial choose_polynomial(
    const Polynomial& upper,
    int lower
) {
    Polynomial result = constant_polynomial(1);
    for (int offset = 0; offset < lower; ++offset) {
        result = multiply(
            result,
            subtract(
                upper,
                constant_polynomial(offset)
            )
        );
    }
    return scale(result, Rational{1, factorial(lower)});
}

Polynomial coefficient_formula(
    int power_value,
    const Polynomial& q,
    const Polynomial& coordinate,
    int quotient
) {
    const Polynomial width = add(scale(q, 2), constant_polynomial(1));
    const Polynomial centered_degree =
        add(scale(q, power_value), coordinate);
    Polynomial result;
    for (int index = 0; index <= quotient; ++index) {
        const Polynomial upper = add(
            subtract(
                centered_degree,
                scale(width, index)
            ),
            constant_polynomial(power_value - 1)
        );
        Polynomial term = scale(
            choose_polynomial(upper, power_value - 1),
            binomial_integer(power_value, index)
        );
        add_to(
            result,
            term,
            index % 2 == 0 ? Rational{1} : Rational{-1}
        );
    }
    return result;
}

Polynomial chamber_value(
    int power_value,
    const Polynomial& q,
    const Polynomial& coordinate,
    int category
) {
    return
        category == power_value
            ? Polynomial{}
            : coefficient_formula(
                power_value,
                q,
                coordinate,
                category
            );
}

ValueTable make_value_table(
    int power_value,
    bool wall_case
) {
    const Polynomial q = variable_polynomial(0);
    const Polynomial x = variable_polynomial(1);
    const Polynomial a = variable_polynomial(2);
    const Polynomial y =
        add(add(x, scale(a, 2)), constant_polynomial(1));
    const std::array<Polynomial, 6> coordinates =
        wall_case
            ? std::array<Polynomial, 6>{
                x,
                x,
                add(x, constant_polynomial(1)),
                y,
                add(y, constant_polynomial(1)),
                add(y, constant_polynomial(2))
            }
            : std::array<Polynomial, 6>{
                subtract(x, constant_polynomial(1)),
                x,
                add(x, constant_polynomial(1)),
                y,
                add(y, constant_polynomial(1)),
                add(y, constant_polynomial(2))
            };
    ValueTable result;
    for (int coordinate = 0; coordinate < 6; ++coordinate) {
        for (int category = 0; category <= power_value; ++category) {
            result[static_cast<std::size_t>(coordinate)].push_back(
                chamber_value(
                    power_value,
                    q,
                    coordinates[static_cast<std::size_t>(coordinate)],
                    category
                )
            );
        }
    }
    return result;
}

Polynomial determinant_polynomial(
    const ValueTable& values,
    const Pattern& pattern,
    bool wall_case
) {
    const Polynomial& vx_minus =
        values[wall_case ? 2U : 0U][
            static_cast<std::size_t>(
                wall_case ? pattern[2] : pattern[0]
            )
        ];
    const Polynomial& vx =
        values[1][static_cast<std::size_t>(pattern[1])];
    const Polynomial& vx_plus =
        values[2][static_cast<std::size_t>(pattern[2])];
    const Polynomial& vy =
        values[3][static_cast<std::size_t>(pattern[3])];
    const Polynomial& vy_plus =
        values[4][static_cast<std::size_t>(pattern[4])];
    const Polynomial& vy_plus_two =
        values[5][static_cast<std::size_t>(pattern[5])];
    return subtract(
        multiply(
            subtract(vx, vy),
            subtract(vx, vy_plus_two)
        ),
        multiply(
            subtract(vx_plus, vy_plus),
            subtract(vx_minus, vy_plus)
        )
    );
}

Polynomial homogenize(
    const Polynomial& polynomial,
    int degree
) {
    Polynomial result;
    for (const auto& [exponent, coefficient] : polynomial.coefficients) {
        const int total_degree =
            exponent[0] + exponent[1] + exponent[2];
        if (total_degree > degree) {
            throw std::runtime_error(
                "determinant exceeds claimed degree"
            );
        }
        const Exponent compact_exponent{
            degree - total_degree,
            exponent[1],
            exponent[2]
        };
        result.coefficients[compact_exponent] += coefficient;
    }
    return result;
}

std::vector<int> next_categories(int category, int power_value) {
    if (category == power_value) {
        return {power_value};
    }
    std::vector<int> result{category};
    if (category + 1 < power_value) {
        result.push_back(category + 1);
    }
    result.push_back(power_value);
    return result;
}

void append_chamber_constraints(
    std::vector<Constraint>& constraints,
    int power_value,
    int category,
    const AffineCoordinate& coordinate
) {
    if (category == power_value) {
        Constraint support;
        support.normal = coordinate.coefficients;
        support.normal[0] -= power_value;
        support.lower = 1 - coordinate.constant;
        constraints.push_back(support);
        return;
    }

    Constraint lower;
    lower.normal = coordinate.coefficients;
    lower.normal[0] += power_value - 2 * category;
    lower.lower = category - coordinate.constant;
    constraints.push_back(lower);

    Constraint upper;
    for (int index = 0; index < 3; ++index) {
        upper.normal[static_cast<std::size_t>(index)] =
            -coordinate.coefficients[static_cast<std::size_t>(index)];
    }
    upper.normal[0] += 2 * category + 2 - power_value;
    upper.lower = coordinate.constant - category;
    constraints.push_back(upper);

    Constraint support_upper;
    for (int index = 0; index < 3; ++index) {
        support_upper.normal[static_cast<std::size_t>(index)] =
            -coordinate.coefficients[static_cast<std::size_t>(index)];
    }
    support_upper.normal[0] += power_value;
    support_upper.lower = coordinate.constant;
    constraints.push_back(support_upper);
}

std::vector<Constraint> make_constraints(
    int power_value,
    const Pattern& pattern,
    bool wall_case
) {
    std::vector<Constraint> constraints;
    constraints.push_back(Constraint{{1, 0, 0}, 1});
    if (wall_case) {
        constraints.push_back(Constraint{{0, 1, 0}, 0});
        constraints.push_back(Constraint{{0, -1, 0}, 0});
    } else {
        constraints.push_back(Constraint{{0, 1, 0}, 1});
    }
    constraints.push_back(Constraint{{0, 0, 1}, 0});
    constraints.push_back(
        Constraint{{power_value, -1, -2}, 1}
    );
    const std::array<AffineCoordinate, 6> coordinates =
        wall_case
            ? std::array<AffineCoordinate, 6>{
                AffineCoordinate{{0, 1, 0}, 0},
                AffineCoordinate{{0, 1, 0}, 0},
                AffineCoordinate{{0, 1, 0}, 1},
                AffineCoordinate{{0, 1, 2}, 1},
                AffineCoordinate{{0, 1, 2}, 2},
                AffineCoordinate{{0, 1, 2}, 3}
            }
            : std::array<AffineCoordinate, 6>{
                AffineCoordinate{{0, 1, 0}, -1},
                AffineCoordinate{{0, 1, 0}, 0},
                AffineCoordinate{{0, 1, 0}, 1},
                AffineCoordinate{{0, 1, 2}, 1},
                AffineCoordinate{{0, 1, 2}, 2},
                AffineCoordinate{{0, 1, 2}, 3}
            };
    for (int index = 0; index < 6; ++index) {
        append_chamber_constraints(
            constraints,
            power_value,
            pattern[static_cast<std::size_t>(index)],
            coordinates[static_cast<std::size_t>(index)]
        );
    }
    return constraints;
}

std::vector<CompactConstraint> compact_constraints(
    const std::vector<Constraint>& constraints
) {
    std::vector<CompactConstraint> result;
    result.push_back(
        CompactConstraint{{1, 0, 0}, 0}
    );
    for (const Constraint& constraint : constraints) {
        result.push_back(CompactConstraint{
            {
                -constraint.lower,
                constraint.normal[1],
                constraint.normal[2]
            },
            constraint.normal[0]
        });
    }
    return result;
}

long long determinant3(
    const std::array<std::array<long long, 3>, 3>& matrix
) {
    return
        matrix[0][0]
            * (
                matrix[1][1] * matrix[2][2]
                - matrix[1][2] * matrix[2][1]
            )
        - matrix[0][1]
            * (
                matrix[1][0] * matrix[2][2]
                - matrix[1][2] * matrix[2][0]
            )
        + matrix[0][2]
            * (
                matrix[1][0] * matrix[2][1]
                - matrix[1][1] * matrix[2][0]
            );
}

bool feasible_point(
    const std::vector<CompactConstraint>& constraints,
    const Point& point
) {
    for (const CompactConstraint& constraint : constraints) {
        Rational value{constraint.constant};
        for (int coordinate = 0; coordinate < 3; ++coordinate) {
            value +=
                constraint.coefficients[
                    static_cast<std::size_t>(coordinate)
                ]
                * point[static_cast<std::size_t>(coordinate)];
        }
        if (value < 0) {
            return false;
        }
    }
    return true;
}

bool satisfies_constraints(
    const std::vector<Constraint>& constraints,
    int q,
    int x,
    int a
) {
    const std::array<long long, 3> point{
        q,
        x,
        a
    };
    for (const Constraint& constraint : constraints) {
        long long value = 0;
        for (int coordinate = 0; coordinate < 3; ++coordinate) {
            value +=
                constraint.normal[
                    static_cast<std::size_t>(coordinate)
                ]
                * point[static_cast<std::size_t>(coordinate)];
        }
        if (value < constraint.lower) {
            return false;
        }
    }
    return true;
}

bool certify_finite_integer_empty(
    const std::vector<Constraint>& constraints,
    const std::vector<Point>& points,
    int power_value,
    bool wall_case
) {
    Rational minimum_r = points.front()[0];
    for (const Point& point : points) {
        minimum_r = std::min(minimum_r, point[0]);
    }
    if (minimum_r <= 0) {
        return false;
    }
    const Rational q_bound = Rational{1} / minimum_r;
    const Integer upper_integer =
        q_bound.numerator() / q_bound.denominator();
    if (
        upper_integer
        > Integer{
            std::numeric_limits<int>::max() / power_value
        }
    ) {
        return false;
    }
    const int upper_q = upper_integer.convert_to<int>();
    for (int q = 1; q <= upper_q; ++q) {
        const int first_x = wall_case ? 0 : 1;
        const int last_x = wall_case ? 0 : power_value * q - 1;
        for (int x = first_x; x <= last_x; ++x) {
            const int last_a =
                (power_value * q - x - 1) / 2;
            for (int a = 0; a <= last_a; ++a) {
                if (satisfies_constraints(constraints, q, x, a)) {
                    return false;
                }
            }
        }
    }
    return true;
}

std::string rational_key(const Rational& value) {
    return
        value.numerator().convert_to<std::string>()
        + "/" + value.denominator().convert_to<std::string>();
}

std::string point_key(const Point& point) {
    return
        rational_key(point[0]) + ","
        + rational_key(point[1]) + ","
        + rational_key(point[2]);
}

std::vector<Point> vertices(
    const std::vector<CompactConstraint>& constraints
) {
    std::vector<Point> result;
    std::set<std::string> seen;
    for (std::size_t first = 0; first < constraints.size(); ++first) {
        for (
            std::size_t second = first + 1;
            second < constraints.size();
            ++second
        ) {
            for (
                std::size_t third = second + 1;
                third < constraints.size();
                ++third
            ) {
                const std::array<std::size_t, 3> selected{
                    first,
                    second,
                    third
                };
                std::array<std::array<long long, 3>, 3> matrix{};
                std::array<long long, 3> right{};
                for (int row = 0; row < 3; ++row) {
                    matrix[static_cast<std::size_t>(row)] =
                        constraints[
                            selected[static_cast<std::size_t>(row)]
                        ].coefficients;
                    right[static_cast<std::size_t>(row)] =
                        -constraints[
                            selected[static_cast<std::size_t>(row)]
                        ].constant;
                }
                const long long denominator = determinant3(matrix);
                if (denominator == 0) {
                    continue;
                }
                std::array<long long, 3> numerators{};
                long long positive_denominator = denominator;
                for (int column = 0; column < 3; ++column) {
                    auto replaced = matrix;
                    for (int row = 0; row < 3; ++row) {
                        replaced[static_cast<std::size_t>(row)]
                            [static_cast<std::size_t>(column)] =
                            right[static_cast<std::size_t>(row)];
                    }
                    numerators[static_cast<std::size_t>(column)] =
                        determinant3(replaced);
                }
                if (positive_denominator < 0) {
                    positive_denominator = -positive_denominator;
                    for (long long& numerator : numerators) {
                        numerator = -numerator;
                    }
                }
                bool feasible = true;
                for (const CompactConstraint& constraint : constraints) {
                    long long value =
                        constraint.constant * positive_denominator;
                    for (int coordinate = 0; coordinate < 3;
                         ++coordinate) {
                        value +=
                            constraint.coefficients[
                                static_cast<std::size_t>(coordinate)
                            ]
                            * numerators[
                                static_cast<std::size_t>(coordinate)
                            ];
                    }
                    feasible = feasible && value >= 0;
                }
                if (!feasible) {
                    continue;
                }
                Point point{};
                for (int coordinate = 0; coordinate < 3; ++coordinate) {
                    point[static_cast<std::size_t>(coordinate)] =
                        Rational{
                            Integer{
                                numerators[
                                    static_cast<std::size_t>(coordinate)
                                ]
                            },
                            Integer{positive_denominator}
                        };
                }
                if (seen.insert(point_key(point)).second) {
                    result.push_back(point);
                }
            }
        }
    }
    return result;
}

int point_rank(const std::vector<Point>& points) {
    if (points.size() <= 1U) {
        return 0;
    }
    std::vector<std::array<Rational, 3>> rows;
    for (std::size_t index = 1; index < points.size(); ++index) {
        std::array<Rational, 3> row{};
        for (int coordinate = 0; coordinate < 3; ++coordinate) {
            row[static_cast<std::size_t>(coordinate)] =
                points[index][static_cast<std::size_t>(coordinate)]
                - points[0][static_cast<std::size_t>(coordinate)];
        }
        rows.push_back(row);
    }
    int rank = 0;
    for (int column = 0; column < 3 && rank < 3; ++column) {
        std::size_t pivot = static_cast<std::size_t>(rank);
        while (
            pivot < rows.size()
            && rows[pivot][static_cast<std::size_t>(column)] == 0
        ) {
            ++pivot;
        }
        if (pivot == rows.size()) {
            continue;
        }
        std::swap(rows[static_cast<std::size_t>(rank)], rows[pivot]);
        const Rational divisor =
            rows[static_cast<std::size_t>(rank)]
                [static_cast<std::size_t>(column)];
        for (int entry = column; entry < 3; ++entry) {
            rows[static_cast<std::size_t>(rank)]
                [static_cast<std::size_t>(entry)] /= divisor;
        }
        for (std::size_t row = 0; row < rows.size(); ++row) {
            if (row == static_cast<std::size_t>(rank)) {
                continue;
            }
            const Rational factor =
                rows[row][static_cast<std::size_t>(column)];
            for (int entry = column; entry < 3; ++entry) {
                rows[row][static_cast<std::size_t>(entry)] -=
                    factor
                    * rows[static_cast<std::size_t>(rank)]
                        [static_cast<std::size_t>(entry)];
            }
        }
        ++rank;
    }
    return rank;
}

Polynomial compose_on_simplex(
    const Polynomial& polynomial,
    const std::vector<Point>& simplex
) {
    const int dimension =
        static_cast<int>(simplex.size()) - 1;
    std::array<Polynomial, 3> coordinates;
    for (int coordinate = 0; coordinate < 3; ++coordinate) {
        coordinates[static_cast<std::size_t>(coordinate)] =
            constant_polynomial(
                simplex[0][static_cast<std::size_t>(coordinate)]
            );
        for (int variable = 0; variable < dimension; ++variable) {
            add_to(
                coordinates[static_cast<std::size_t>(coordinate)],
                variable_polynomial(variable),
                simplex[static_cast<std::size_t>(variable + 1)]
                    [static_cast<std::size_t>(coordinate)]
                - simplex[0][static_cast<std::size_t>(coordinate)]
            );
        }
    }
    std::array<std::vector<Polynomial>, 3> coordinate_powers;
    std::array<int, 3> maximum_exponents{0, 0, 0};
    for (const auto& [exponent, coefficient] : polynomial.coefficients) {
        (void)coefficient;
        for (int coordinate = 0; coordinate < 3; ++coordinate) {
            maximum_exponents[
                static_cast<std::size_t>(coordinate)
            ] = std::max(
                maximum_exponents[
                    static_cast<std::size_t>(coordinate)
                ],
                exponent[static_cast<std::size_t>(coordinate)]
            );
        }
    }
    for (int coordinate = 0; coordinate < 3; ++coordinate) {
        coordinate_powers[static_cast<std::size_t>(coordinate)]
            .push_back(constant_polynomial(1));
        for (
            int exponent = 1;
            exponent
                <= maximum_exponents[
                    static_cast<std::size_t>(coordinate)
                ];
            ++exponent
        ) {
            coordinate_powers[static_cast<std::size_t>(coordinate)]
                .push_back(
                    multiply(
                        coordinate_powers[
                            static_cast<std::size_t>(coordinate)
                        ].back(),
                        coordinates[
                            static_cast<std::size_t>(coordinate)
                        ]
                    )
                );
        }
    }
    Polynomial result;
    for (const auto& [exponent, coefficient] : polynomial.coefficients) {
        Polynomial term = multiply(
            coordinate_powers[0][
                static_cast<std::size_t>(exponent[0])
            ],
            coordinate_powers[1][
                static_cast<std::size_t>(exponent[1])
            ]
        );
        term = multiply(
            term,
            coordinate_powers[2][
                static_cast<std::size_t>(exponent[2])
            ]
        );
        add_to(result, term, coefficient);
    }
    return result;
}

std::vector<Rational> bernstein_coefficients(
    const Polynomial& polynomial,
    int dimension,
    int degree
) {
    const int width = degree + 1;
    const std::size_t cube_size =
        static_cast<std::size_t>(width)
        * static_cast<std::size_t>(width)
        * static_cast<std::size_t>(width);
    const auto cube_index = [width](const Exponent& exponent) {
        return static_cast<std::size_t>(
            (exponent[0] * width + exponent[1]) * width
            + exponent[2]
        );
    };

    std::vector<std::vector<Integer>> binomials(
        static_cast<std::size_t>(width),
        std::vector<Integer>(static_cast<std::size_t>(width))
    );
    for (int upper = 0; upper <= degree; ++upper) {
        for (int lower = 0; lower <= upper; ++lower) {
            binomials[static_cast<std::size_t>(upper)]
                [static_cast<std::size_t>(lower)] =
                binomial_integer(upper, lower);
        }
    }

    std::vector<Rational> grid(cube_size);
    for (const auto& [beta, coefficient] : polynomial.coefficients) {
        const int beta_total = beta[0] + beta[1] + beta[2];
        if (beta_total > degree) {
            throw std::runtime_error(
                "composed polynomial degree overflow"
            );
        }
        Integer numerator = 1;
        for (int coordinate = 0; coordinate < dimension; ++coordinate) {
            numerator *= factorial(
                beta[static_cast<std::size_t>(coordinate)]
            );
        }
        for (int coordinate = dimension; coordinate < 3; ++coordinate) {
            if (beta[static_cast<std::size_t>(coordinate)] != 0) {
                throw std::runtime_error(
                    "inactive simplex variable present"
                );
            }
        }
        Integer denominator = 1;
        for (int offset = 0; offset < beta_total; ++offset) {
            denominator *= degree - offset;
        }
        grid[cube_index(beta)] +=
            coefficient * Rational{numerator, denominator};
    }

    for (int axis = 0; axis < dimension; ++axis) {
        std::vector<Rational> transformed(cube_size);
        for (int first = 0; first <= degree; ++first) {
            for (int second = 0; second <= degree; ++second) {
                for (int third = 0; third <= degree; ++third) {
                    Exponent alpha{first, second, third};
                    Rational value = 0;
                    for (
                        int beta_coordinate = 0;
                        beta_coordinate
                            <= alpha[static_cast<std::size_t>(axis)];
                        ++beta_coordinate
                    ) {
                        Exponent beta = alpha;
                        beta[static_cast<std::size_t>(axis)] =
                            beta_coordinate;
                        value +=
                            binomials[
                                static_cast<std::size_t>(
                                    alpha[
                                        static_cast<std::size_t>(axis)
                                    ]
                                )
                            ][static_cast<std::size_t>(beta_coordinate)]
                            * grid[cube_index(beta)];
                    }
                    transformed[cube_index(alpha)] = value;
                }
            }
        }
        grid = std::move(transformed);
    }

    std::vector<Rational> result;
    const int first_limit = dimension >= 1 ? degree : 0;
    for (int first = 0; first <= first_limit; ++first) {
        const int second_limit =
            dimension >= 2 ? degree - first : 0;
        for (int second = 0; second <= second_limit; ++second) {
            const int third_limit =
                dimension >= 3
                    ? degree - first - second
                    : 0;
            for (int third = 0; third <= third_limit; ++third) {
                result.push_back(
                    grid[cube_index({first, second, third})]
                );
            }
        }
    }
    return result;
}

std::vector<Rational> reference_bernstein_coefficients(
    const Polynomial& polynomial,
    int dimension,
    int degree
) {
    std::vector<Rational> result;
    const int first_limit = dimension >= 1 ? degree : 0;
    for (int first = 0; first <= first_limit; ++first) {
        const int second_limit =
            dimension >= 2 ? degree - first : 0;
        for (int second = 0; second <= second_limit; ++second) {
            const int third_limit =
                dimension >= 3
                    ? degree - first - second
                    : 0;
            for (int third = 0; third <= third_limit; ++third) {
                const Exponent alpha{first, second, third};
                Rational value = 0;
                for (
                    const auto& [beta, coefficient] :
                    polynomial.coefficients
                ) {
                    const int beta_total =
                        beta[0] + beta[1] + beta[2];
                    bool admissible = beta_total <= degree;
                    Integer numerator = 1;
                    for (
                        int coordinate = 0;
                        coordinate < dimension;
                        ++coordinate
                    ) {
                        const int beta_coordinate =
                            beta[static_cast<std::size_t>(coordinate)];
                        const int alpha_coordinate =
                            alpha[static_cast<std::size_t>(coordinate)];
                        admissible =
                            admissible
                            && beta_coordinate <= alpha_coordinate;
                        for (
                            int offset = 0;
                            offset < beta_coordinate;
                            ++offset
                        ) {
                            numerator *= alpha_coordinate - offset;
                        }
                    }
                    for (int coordinate = dimension; coordinate < 3;
                         ++coordinate) {
                        admissible =
                            admissible
                            && beta[
                                static_cast<std::size_t>(coordinate)
                            ] == 0;
                    }
                    if (!admissible) {
                        continue;
                    }
                    Integer denominator = 1;
                    for (int offset = 0; offset < beta_total; ++offset) {
                        denominator *= degree - offset;
                    }
                    value += coefficient
                        * Rational{numerator, denominator};
                }
                result.push_back(value);
            }
        }
    }
    return result;
}

void verify_bernstein_transform() {
    constexpr int degree = 4;
    for (int dimension = 0; dimension <= 3; ++dimension) {
        Polynomial test;
        for (int first = 0; first <= degree; ++first) {
            for (int second = 0; second <= degree - first; ++second) {
                for (
                    int third = 0;
                    third <= degree - first - second;
                    ++third
                ) {
                    const Exponent exponent{first, second, third};
                    bool active = true;
                    for (int coordinate = dimension; coordinate < 3;
                         ++coordinate) {
                        active =
                            active
                            && exponent[
                                static_cast<std::size_t>(coordinate)
                            ] == 0;
                    }
                    if (!active) {
                        continue;
                    }
                    test.coefficients[exponent] =
                        1 + 2 * first + 3 * second + 5 * third;
                }
            }
        }
        if (
            bernstein_coefficients(test, dimension, degree)
            != reference_bernstein_coefficients(
                test,
                dimension,
                degree
            )
        ) {
            throw std::runtime_error(
                "separable Bernstein transform replay failed"
            );
        }
    }
}

Rational squared_distance(
    const Point& left,
    const Point& right
) {
    Rational result = 0;
    for (int coordinate = 0; coordinate < 3; ++coordinate) {
        const Rational difference =
            left[static_cast<std::size_t>(coordinate)]
            - right[static_cast<std::size_t>(coordinate)];
        result += difference * difference;
    }
    return result;
}

std::pair<std::size_t, std::size_t> longest_edge(
    const std::vector<Point>& simplex
) {
    std::pair<std::size_t, std::size_t> result{0U, 1U};
    Rational maximum = -1;
    for (std::size_t first = 0; first < simplex.size(); ++first) {
        for (
            std::size_t second = first + 1;
            second < simplex.size();
            ++second
        ) {
            const Rational distance =
                squared_distance(simplex[first], simplex[second]);
            if (distance > maximum) {
                maximum = distance;
                result = {first, second};
            }
        }
    }
    return result;
}

void certify_simplex_recursive(
    const Polynomial& polynomial,
    const std::vector<Point>& simplex,
    int degree,
    int depth,
    int maximum_depth,
    SimplexResult& result
) {
    ++result.nodes;
    result.maximum_depth = std::max(result.maximum_depth, depth);
    const int dimension =
        static_cast<int>(simplex.size()) - 1;
    const Polynomial composed =
        compose_on_simplex(polynomial, simplex);
    const std::vector<Rational> coefficients =
        bernstein_coefficients(composed, dimension, degree);
    bool nonnegative = true;
    for (const Rational& coefficient : coefficients) {
        if (
            !result.minimum_initialized
            || coefficient < result.minimum_coefficient
        ) {
            result.minimum_coefficient = coefficient;
            result.minimum_initialized = true;
        }
        nonnegative = nonnegative && coefficient >= 0;
    }
    if (nonnegative) {
        ++result.leaves;
        return;
    }
    if (depth == maximum_depth || dimension == 0) {
        ++result.negative_leaves;
        return;
    }

    ++result.splits;
    const auto [first, second] = longest_edge(simplex);
    Point midpoint{};
    for (int coordinate = 0; coordinate < 3; ++coordinate) {
        midpoint[static_cast<std::size_t>(coordinate)] =
            (
                simplex[first][static_cast<std::size_t>(coordinate)]
                + simplex[second][static_cast<std::size_t>(coordinate)]
            ) / 2;
    }
    std::vector<Point> left = simplex;
    std::vector<Point> right = simplex;
    left[first] = midpoint;
    right[second] = midpoint;
    certify_simplex_recursive(
        polynomial,
        left,
        degree,
        depth + 1,
        maximum_depth,
        result
    );
    certify_simplex_recursive(
        polynomial,
        right,
        degree,
        depth + 1,
        maximum_depth,
        result
    );
}

bool affinely_independent(const std::vector<Point>& simplex) {
    return point_rank(simplex)
        == static_cast<int>(simplex.size()) - 1;
}

std::vector<std::vector<Point>> covering_simplices(
    const std::vector<Point>& points,
    const std::vector<CompactConstraint>& constraints
) {
    const int dimension = point_rank(points);
    std::vector<std::vector<Point>> result;
    if (dimension == 0) {
        result.push_back({points.front()});
        return result;
    }
    if (dimension <= 2) {
        for (std::size_t first = 0; first < points.size(); ++first) {
            for (
                std::size_t second = first + 1;
                second < points.size();
                ++second
            ) {
                if (dimension == 1) {
                    std::vector<Point> simplex{
                        points[first],
                        points[second]
                    };
                    if (affinely_independent(simplex)) {
                        result.push_back(std::move(simplex));
                    }
                    continue;
                }
                for (
                    std::size_t third = second + 1;
                    third < points.size();
                    ++third
                ) {
                    std::vector<Point> simplex{
                        points[first],
                        points[second],
                        points[third]
                    };
                    if (affinely_independent(simplex)) {
                        result.push_back(std::move(simplex));
                    }
                }
            }
        }
        return result;
    }

    Point center{};
    for (const Point& point : points) {
        for (int coordinate = 0; coordinate < 3; ++coordinate) {
            center[static_cast<std::size_t>(coordinate)] +=
                point[static_cast<std::size_t>(coordinate)];
        }
    }
    for (Rational& coordinate : center) {
        coordinate /= static_cast<int>(points.size());
    }

    std::set<std::string> seen;
    for (const CompactConstraint& constraint : constraints) {
        std::vector<Point> active;
        for (const Point& point : points) {
            Rational value{constraint.constant};
            for (int coordinate = 0; coordinate < 3; ++coordinate) {
                value +=
                    constraint.coefficients[
                        static_cast<std::size_t>(coordinate)
                    ]
                    * point[static_cast<std::size_t>(coordinate)];
            }
            if (value == 0) {
                active.push_back(point);
            }
        }
        for (std::size_t first = 0; first < active.size(); ++first) {
            for (
                std::size_t second = first + 1;
                second < active.size();
                ++second
            ) {
                for (
                    std::size_t third = second + 1;
                    third < active.size();
                    ++third
                ) {
                    std::vector<Point> simplex{
                        center,
                        active[first],
                        active[second],
                        active[third]
                    };
                    if (!affinely_independent(simplex)) {
                        continue;
                    }
                    std::array<std::string, 3> keys{
                        point_key(active[first]),
                        point_key(active[second]),
                        point_key(active[third])
                    };
                    std::sort(keys.begin(), keys.end());
                    const std::string key =
                        keys[0] + "|" + keys[1] + "|" + keys[2];
                    if (seen.insert(key).second) {
                        result.push_back(std::move(simplex));
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
        std::cout << std::unitbuf;
        verify_bernstein_transform();
        if (argc < 2 || argc > 5) {
            throw std::runtime_error(
                "usage: prove_su2_power_tp2_compact_bernstein "
                "power [maximum_depth] [feasible_pattern_index] "
                "[wall]"
            );
        }
        if (
            argc == 5
            && std::string{argv[4]} != "wall"
        ) {
            throw std::runtime_error(
                "the optional domain must be wall"
            );
        }
        const bool wall_case = argc == 5;
        const int power_value = parse_bounded_integer(
            argv[1],
            2,
            20,
            "power"
        );
        const int maximum_depth =
            argc >= 3
                ? parse_bounded_integer(
                    argv[2],
                    0,
                    30,
                    "maximum_depth"
                )
                : 0;
        const int selected_pattern =
            argc >= 4
                ? parse_bounded_integer(
                    argv[3],
                    -1,
                    1000000,
                    "feasible_pattern_index"
                )
                : -1;
        const int degree = 2 * (power_value - 1);
        const ValueTable value_table =
            make_value_table(power_value, wall_case);

        std::size_t patterns = 0;
        std::size_t feasible_patterns = 0;
        std::size_t analyzed_patterns = 0;
        std::size_t certified_patterns = 0;
        std::size_t integer_empty_patterns = 0;
        std::size_t compact_vertices = 0;
        std::size_t covering_simplex_count = 0;
        SimplexResult total;
        std::string first_failed_pattern;

        for (int x_minus = 0; x_minus < power_value; ++x_minus) {
            for (
                const int x_category :
                next_categories(x_minus, power_value)
            ) {
                if (x_category == power_value) {
                    continue;
                }
                if (wall_case && x_category != x_minus) {
                    continue;
                }
                for (
                    const int x_plus :
                    next_categories(x_category, power_value)
                ) {
                    for (int y_category = 0;
                         y_category < power_value;
                         ++y_category) {
                        for (
                            const int y_plus :
                            next_categories(y_category, power_value)
                        ) {
                            for (
                                const int y_plus_two :
                                next_categories(y_plus, power_value)
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
                                    chamber_constraints =
                                        make_constraints(
                                            power_value,
                                            pattern,
                                            wall_case
                                        );
                                const std::vector<CompactConstraint>
                                    compact =
                                        compact_constraints(
                                            chamber_constraints
                                        );
                                const std::vector<Point> points =
                                    vertices(compact);
                                bool has_finite_point = false;
                                for (const Point& point : points) {
                                    has_finite_point =
                                        has_finite_point
                                        || point[0] > 0;
                                }
                                if (!has_finite_point) {
                                    continue;
                                }
                                const int pattern_index =
                                    static_cast<int>(
                                        feasible_patterns
                                    );
                                ++feasible_patterns;
                                if (
                                    selected_pattern >= 0
                                    && pattern_index
                                        != selected_pattern
                                ) {
                                    continue;
                                }
                                ++analyzed_patterns;

                                if (
                                    !feasible_point(
                                        compact,
                                        points.front()
                                    )
                                ) {
                                    throw std::runtime_error(
                                        "vertex feasibility replay failed"
                                    );
                                }
                                compact_vertices += points.size();
                                const std::vector<std::vector<Point>>
                                    simplices =
                                        covering_simplices(
                                            points,
                                            compact
                                        );
                                if (simplices.empty()) {
                                    throw std::runtime_error(
                                        "compact chamber has no cover"
                                    );
                                }
                                covering_simplex_count +=
                                    simplices.size();

                                const Polynomial determinant =
                                    determinant_polynomial(
                                        value_table,
                                        pattern,
                                        wall_case
                                    );
                                const Polynomial homogeneous =
                                    homogenize(determinant, degree);
                                SimplexResult chamber_result;
                                for (
                                    const std::vector<Point>& simplex :
                                    simplices
                                ) {
                                    certify_simplex_recursive(
                                        homogeneous,
                                        simplex,
                                        degree,
                                        0,
                                        maximum_depth,
                                        chamber_result
                                    );
                                }
                                const bool integer_empty =
                                    chamber_result.negative_leaves != 0
                                    && certify_finite_integer_empty(
                                        chamber_constraints,
                                        points,
                                        power_value,
                                        wall_case
                                    );
                                chamber_result.certified =
                                    chamber_result.negative_leaves == 0
                                    || integer_empty;
                                if (integer_empty) {
                                    ++integer_empty_patterns;
                                }
                                if (chamber_result.certified) {
                                    ++certified_patterns;
                                } else if (first_failed_pattern.empty()) {
                                    first_failed_pattern =
                                        show_pattern(pattern);
                                }
                                total.nodes += chamber_result.nodes;
                                total.leaves += chamber_result.leaves;
                                total.splits += chamber_result.splits;
                                total.negative_leaves +=
                                    chamber_result.negative_leaves;
                                total.maximum_depth = std::max(
                                    total.maximum_depth,
                                    chamber_result.maximum_depth
                                );
                                if (
                                    chamber_result.minimum_initialized
                                    && (
                                        !total.minimum_initialized
                                        || chamber_result
                                                .minimum_coefficient
                                            < total.minimum_coefficient
                                    )
                                ) {
                                    total.minimum_coefficient =
                                        chamber_result
                                            .minimum_coefficient;
                                    total.minimum_initialized = true;
                                }

                                std::cout
                                    << "SU2_POWER_TP2_COMPACT_BERNSTEIN"
                                    << " power=" << power_value
                                    << " domain="
                                    << (
                                        wall_case
                                            ? "wall"
                                            : "interior"
                                    )
                                    << " feasible_pattern_index="
                                    << pattern_index
                                    << " pattern="
                                    << show_pattern(pattern)
                                    << " dimension="
                                    << point_rank(points)
                                    << " vertices=" << points.size()
                                    << " simplices="
                                    << simplices.size()
                                    << " nodes="
                                    << chamber_result.nodes
                                    << " negative_leaves="
                                    << chamber_result.negative_leaves
                                    << " result="
                                    << (
                                        chamber_result.certified
                                            ? (
                                                integer_empty
                                                    ? "PASS_INTEGER_EMPTY_CHAMBER"
                                                    : "PASS_REAL_CHAMBER"
                                            )
                                            : "PARTIAL"
                                    )
                                    << '\n';
                            }
                        }
                    }
                }
            }
        }

        const bool all_selected =
            selected_pattern < 0 || analyzed_patterns == 1;
        const bool proved =
            all_selected
            && certified_patterns == analyzed_patterns
            && analyzed_patterns > 0;
        std::cout
            << "SU2_POWER_TP2_COMPACT_BERNSTEIN_SUMMARY"
            << " power=" << power_value
            << " domain="
            << (wall_case ? "wall" : "interior")
            << " maximum_depth=" << maximum_depth
            << " patterns=" << patterns
            << " feasible_patterns=" << feasible_patterns
            << " analyzed_patterns=" << analyzed_patterns
            << " certified_patterns=" << certified_patterns
            << " integer_empty_patterns="
            << integer_empty_patterns
            << " compact_vertices=" << compact_vertices
            << " covering_simplices=" << covering_simplex_count
            << " nodes=" << total.nodes
            << " leaves=" << total.leaves
            << " splits=" << total.splits
            << " negative_leaves=" << total.negative_leaves
            << " maximum_reached_depth=" << total.maximum_depth
            << " minimum_bernstein_coefficient="
            << (
                total.minimum_initialized
                    ? rational_key(total.minimum_coefficient)
                    : "{}"
            )
            << " first_failed_pattern="
            << (
                first_failed_pattern.empty()
                    ? "{}"
                    : first_failed_pattern
            )
            << " result="
            << (
                proved
                    ? (
                        integer_empty_patterns == 0
                            ? "PASS_REAL_COMPACT_CERTIFICATE"
                            : "PASS_INTEGER_COMPACT_CERTIFICATE"
                    )
                    : "PARTIAL"
            )
            << '\n';
        return proved ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
