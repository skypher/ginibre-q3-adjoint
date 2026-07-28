#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>
#include <z3++.h>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using Exponent = std::array<int, 3>;

class Polynomial {
public:
    Polynomial() = default;

    explicit Polynomial(const Rational& value) {
        if (value != 0) {
            terms_[Exponent{0, 0, 0}] = value;
        }
    }

    static Polynomial variable(int index) {
        Polynomial result;
        Exponent exponent{0, 0, 0};
        exponent[static_cast<std::size_t>(index)] = 1;
        result.terms_[exponent] = Rational(1);
        return result;
    }

    Polynomial& operator+=(const Polynomial& other) {
        for (const auto& [exponent, coefficient] : other.terms_) {
            terms_[exponent] += coefficient;
            if (terms_[exponent] == 0) {
                terms_.erase(exponent);
            }
        }
        return *this;
    }

    Polynomial& operator-=(const Polynomial& other) {
        for (const auto& [exponent, coefficient] : other.terms_) {
            terms_[exponent] -= coefficient;
            if (terms_[exponent] == 0) {
                terms_.erase(exponent);
            }
        }
        return *this;
    }

    Polynomial& operator*=(const Polynomial& other) {
        std::map<Exponent, Rational> product;
        for (const auto& [left_exponent, left_coefficient] : terms_) {
            for (const auto& [right_exponent, right_coefficient]
                 : other.terms_) {
                Exponent exponent{0, 0, 0};
                for (std::size_t index = 0U;
                     index < exponent.size();
                     ++index) {
                    exponent[index] =
                        left_exponent[index] + right_exponent[index];
                }
                product[exponent] +=
                    left_coefficient * right_coefficient;
            }
        }
        terms_.clear();
        for (const auto& [exponent, coefficient] : product) {
            if (coefficient != 0) {
                terms_[exponent] = coefficient;
            }
        }
        return *this;
    }

    const std::map<Exponent, Rational>& terms() const {
        return terms_;
    }

    bool operator<(const Polynomial& other) const {
        return terms_ < other.terms_;
    }

private:
    std::map<Exponent, Rational> terms_;
};

Polynomial operator+(Polynomial left, const Polynomial& right) {
    left += right;
    return left;
}

Polynomial operator-(Polynomial left, const Polynomial& right) {
    left -= right;
    return left;
}

Polynomial operator*(Polynomial left, const Polynomial& right) {
    left *= right;
    return left;
}

Polynomial constant(long value) {
    return Polynomial(Rational(value));
}

Polynomial scale(const Polynomial& value, long factor) {
    return value * Polynomial(Rational(factor));
}

Polynomial binomial(const Polynomial& top, int order) {
    Polynomial result = constant(1);
    Integer factorial = 1;
    for (int index = 0; index < order; ++index) {
        result *= top - constant(index);
        factorial *= index + 1;
    }
    return result * Polynomial(Rational(1, factorial));
}

Polynomial power(const Polynomial& value, int exponent) {
    Polynomial result = constant(1);
    for (int index = 0; index < exponent; ++index) {
        result *= value;
    }
    return result;
}

Polynomial substitute(
    const Polynomial& polynomial,
    const std::array<Polynomial, 3>& values
) {
    Polynomial result;
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        Polynomial term(coefficient);
        for (std::size_t index = 0U;
             index < exponent.size();
             ++index) {
            term *= power(values[index], exponent[index]);
        }
        result += term;
    }
    return result;
}

struct Hinge {
    std::string name;
    Polynomial top;
    int order;
};

struct Chamber {
    std::uint64_t mask;
    std::vector<Polynomial> constraints;
    Polynomial margin;
};

Polynomial selected_binomial(
    const std::vector<Hinge>& hinges,
    std::uint64_t mask,
    std::size_t index
) {
    if ((mask & (std::uint64_t{1} << index)) == 0U) {
        return Polynomial();
    }
    return binomial(hinges[index].top, hinges[index].order);
}

Polynomial multiplicity_six(
    const std::vector<Hinge>& hinges,
    std::uint64_t mask
) {
    return selected_binomial(hinges, mask, 0)
        - scale(selected_binomial(hinges, mask, 1), 6)
        + scale(selected_binomial(hinges, mask, 2), 15);
}

Polynomial multiplicity_seven(
    const std::vector<Hinge>& hinges,
    std::uint64_t mask
) {
    return selected_binomial(hinges, mask, 6)
        - scale(selected_binomial(hinges, mask, 7), 7)
        + scale(selected_binomial(hinges, mask, 8), 21)
        - scale(selected_binomial(hinges, mask, 9), 35);
}

Chamber make_chamber(
    std::uint64_t mask,
    const std::vector<Hinge>& hinges
) {
    const Polynomial q = Polynomial::variable(0);
    const Polynomial h = Polynomial::variable(1);
    const Polynomial y = Polynomial::variable(2);

    Chamber chamber{
        mask,
        {
            q - constant(1),
            h - constant(1),
            scale(q, 5) - constant(2) - h,
            y,
            scale(q, 2) + constant(1) + h - y
        },
        Polynomial()
    };
    for (std::size_t index = 0U; index < hinges.size(); ++index) {
        chamber.constraints.push_back(
            (mask & (std::uint64_t{1} << index)) != 0U
                ? hinges[index].top - constant(hinges[index].order)
                : constant(hinges[index].order - 1)
                    - hinges[index].top
        );
    }

    const Polynomial u6 =
        multiplicity_six(hinges, mask)
        - selected_binomial(hinges, mask, 3)
        + scale(selected_binomial(hinges, mask, 4), 6)
        + selected_binomial(hinges, mask, 5);
    const Polynomial u7 =
        multiplicity_seven(hinges, mask)
        - selected_binomial(hinges, mask, 10)
        + scale(selected_binomial(hinges, mask, 11), 7)
        - scale(selected_binomial(hinges, mask, 12), 21)
        + selected_binomial(hinges, mask, 13)
        - scale(selected_binomial(hinges, mask, 14), 7)
        - selected_binomial(hinges, mask, 15);

    const Polynomial stable_f6 =
        binomial(scale(q, 6) + constant(4), 4)
        - scale(binomial(scale(q, 4) + constant(3), 4), 6)
        + scale(binomial(scale(q, 2) + constant(2), 4), 15);
    const Polynomial f6 =
        stable_f6
        - selected_binomial(hinges, mask, 16)
        + (
            (mask & (std::uint64_t{1} << 16)) != 0U
                ? binomial(scale(q, 2) - scale(h, 2), 4)
                : Polynomial()
        );

    const Polynomial stable_f7 =
        binomial(scale(q, 7) + constant(5), 5)
        - scale(binomial(scale(q, 5) + constant(4), 5), 7)
        + scale(binomial(scale(q, 3) + constant(3), 5), 21)
        - scale(selected_binomial(hinges, mask, 21), 35);
    const Polynomial f7 =
        stable_f7
        - selected_binomial(hinges, mask, 17)
        + scale(selected_binomial(hinges, mask, 18), 7)
        + selected_binomial(hinges, mask, 19)
        - scale(selected_binomial(hinges, mask, 20), 7);

    chamber.margin = f6 * u7 - f7 * u6;
    return chamber;
}

z3::expr rational_value(z3::context& context, const Rational& value) {
    std::ostringstream stream;
    stream << value.numerator();
    if (value.denominator() != 1) {
        stream << '/' << value.denominator();
    }
    return context.real_val(stream.str().c_str());
}

Integer numeral_integer(const z3::expr& value) {
    std::string text;
    if (!value.is_numeral(text)) {
        throw std::runtime_error("Z3 model value is not a numeral");
    }
    return Integer(text);
}

Rational model_rational(const z3::expr& value) {
    if (!value.is_numeral()) {
        throw std::runtime_error("Z3 model value is not rational");
    }
    return Rational(
        numeral_integer(value.numerator()),
        numeral_integer(value.denominator())
    );
}

z3::expr affine_integer_value(
    z3::context& context,
    const Polynomial& polynomial,
    const std::array<z3::expr, 3>& variables
) {
    z3::expr result = context.int_val(0);
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        if (coefficient.denominator() != 1) {
            throw std::runtime_error(
                "affine chamber constraint is not integral"
            );
        }
        int total_degree = 0;
        std::size_t variable_index = 0U;
        for (std::size_t index = 0U;
             index < exponent.size();
             ++index) {
            total_degree += exponent[index];
            if (exponent[index] == 1) {
                variable_index = index;
            } else if (exponent[index] != 0) {
                throw std::runtime_error(
                    "chamber constraint is not affine"
                );
            }
        }
        if (total_degree > 1) {
            throw std::runtime_error(
                "chamber constraint is not affine"
            );
        }
        std::ostringstream stream;
        stream << coefficient.numerator();
        const z3::expr coefficient_value =
            context.int_val(stream.str().c_str());
        result = result + (
            total_degree == 0
                ? coefficient_value
                : coefficient_value * variables[variable_index]
        );
    }
    return result;
}

std::vector<Polynomial> irredundant_constraints(
    const Chamber& chamber
) {
    z3::context context;
    const std::array<z3::expr, 3> variables{
        context.int_const("Q"),
        context.int_const("H"),
        context.int_const("Y")
    };
    std::vector<std::size_t> retained(chamber.constraints.size());
    std::iota(retained.begin(), retained.end(), 0U);
    for (std::size_t position = 0U;
         position < retained.size();) {
        z3::solver solver(context);
        for (std::size_t other = 0U;
             other < retained.size();
             ++other) {
            if (other != position) {
                solver.add(
                    affine_integer_value(
                        context,
                        chamber.constraints[retained[other]],
                        variables
                    ) >= 0
                );
            }
        }
        solver.add(
            affine_integer_value(
                context,
                chamber.constraints[retained[position]],
                variables
            ) < 0
        );
        if (solver.check() == z3::unsat) {
            retained.erase(
                retained.begin()
                    + static_cast<std::ptrdiff_t>(position)
            );
        } else {
            ++position;
        }
    }
    std::vector<Polynomial> result;
    result.reserve(retained.size());
    for (const std::size_t index : retained) {
        result.push_back(chamber.constraints[index]);
    }
    return result;
}

void enumerate_products_recursive(
    const std::vector<Polynomial>& constraints,
    int remaining_degree,
    std::size_t minimum_index,
    const Polynomial& current,
    std::map<Polynomial, Polynomial>& products
) {
    products.emplace(current, current);
    if (remaining_degree == 0) {
        return;
    }
    for (std::size_t index = minimum_index;
         index < constraints.size();
         ++index) {
        enumerate_products_recursive(
            constraints,
            remaining_degree - 1,
            index,
            current * constraints[index],
            products
        );
    }
}

std::array<Rational, 4> affine_coefficients(
    const Polynomial& polynomial
) {
    std::array<Rational, 4> result{
        Rational(0),
        Rational(0),
        Rational(0),
        Rational(0)
    };
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        int total_degree = 0;
        std::size_t variable_index = 0U;
        for (std::size_t index = 0U;
             index < exponent.size();
             ++index) {
            total_degree += exponent[index];
            if (exponent[index] == 1) {
                variable_index = index;
            } else if (exponent[index] != 0) {
                throw std::runtime_error(
                    "constraint is not affine"
                );
            }
        }
        if (total_degree == 0) {
            result[0] += coefficient;
        } else if (total_degree == 1) {
            result[variable_index + 1U] += coefficient;
        } else {
            throw std::runtime_error("constraint is not affine");
        }
    }
    return result;
}

bool nonnegative_newton_coefficients(
    const Polynomial& polynomial,
    std::size_t& positive_coefficients
) {
    constexpr int maximum_degree = 9;
    std::array<std::array<Integer, maximum_degree + 1U>,
               maximum_degree + 1U> stirling{};
    stirling[0][0] = 1;
    for (int degree = 1; degree <= maximum_degree; ++degree) {
        for (int parts = 1; parts <= degree; ++parts) {
            stirling[static_cast<std::size_t>(degree)]
                     [static_cast<std::size_t>(parts)] =
                stirling[static_cast<std::size_t>(degree - 1)]
                         [static_cast<std::size_t>(parts - 1)]
                + Integer(parts)
                    * stirling[static_cast<std::size_t>(degree - 1)]
                              [static_cast<std::size_t>(parts)];
        }
    }
    std::array<Integer, maximum_degree + 1U> factorial{};
    factorial[0] = 1;
    for (int degree = 1; degree <= maximum_degree; ++degree) {
        factorial[static_cast<std::size_t>(degree)] =
            Integer(degree)
            * factorial[static_cast<std::size_t>(degree - 1)];
    }

    std::map<Exponent, Rational> newton;
    for (const auto& [powers, coefficient] : polynomial.terms()) {
        if (
            powers[0] > maximum_degree
            || powers[1] > maximum_degree
            || powers[2] > maximum_degree
        ) {
            throw std::runtime_error(
                "Newton conversion exceeds degree table"
            );
        }
        for (int first = 0; first <= powers[0]; ++first) {
            for (int second = 0; second <= powers[1]; ++second) {
                for (int third = 0; third <= powers[2]; ++third) {
                    const Integer multiplier =
                        stirling[
                            static_cast<std::size_t>(powers[0])
                        ][static_cast<std::size_t>(first)]
                        * factorial[static_cast<std::size_t>(first)]
                        * stirling[
                            static_cast<std::size_t>(powers[1])
                        ][static_cast<std::size_t>(second)]
                        * factorial[static_cast<std::size_t>(second)]
                        * stirling[
                            static_cast<std::size_t>(powers[2])
                        ][static_cast<std::size_t>(third)]
                        * factorial[static_cast<std::size_t>(third)];
                    if (multiplier != 0) {
                        newton[Exponent{first, second, third}] +=
                            coefficient * Rational(multiplier);
                    }
                }
            }
        }
    }
    positive_coefficients = 0U;
    for (const auto& [exponent, coefficient] : newton) {
        static_cast<void>(exponent);
        if (coefficient < 0) {
            return false;
        }
        if (coefficient > 0) {
            ++positive_coefficients;
        }
    }
    return true;
}

bool direct_basis_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints,
    bool emit = true,
    int zero_coordinate_count = 0
) {
    if (constraints.size() < 3U) {
        return false;
    }
    const Polynomial x0 = Polynomial::variable(0);
    const Polynomial x1 = Polynomial::variable(1);
    const Polynomial x2 = Polynomial::variable(2);
    const std::array<Polynomial, 3> basis_variables{x0, x1, x2};
    for (std::size_t first = 0U;
         first + 2U < constraints.size();
         ++first) {
        if (zero_coordinate_count >= 1 && first != 0U) {
            continue;
        }
        for (std::size_t second = first + 1U;
             second + 1U < constraints.size();
             ++second) {
            if (zero_coordinate_count >= 2 && second != 1U) {
                continue;
            }
            for (std::size_t third = second + 1U;
                 third < constraints.size();
                 ++third) {
                if (zero_coordinate_count >= 3 && third != 2U) {
                    continue;
                }
                const std::array<std::size_t, 3> indices{
                    first,
                    second,
                    third
                };
                std::array<std::array<Rational, 7>, 3> rows{};
                for (std::size_t row = 0U; row < 3U; ++row) {
                    const auto affine = affine_coefficients(
                        constraints[indices[row]]
                    );
                    for (std::size_t column = 0U;
                         column < 3U;
                         ++column) {
                        rows[row][column] = affine[column + 1U];
                    }
                    rows[row][3U] = -affine[0];
                    for (std::size_t column = 0U;
                         column < 3U;
                         ++column) {
                        rows[row][column + 4U] =
                            row == column ? Rational(1) : Rational(0);
                    }
                }
                bool invertible = true;
                for (std::size_t column = 0U;
                     column < 3U;
                     ++column) {
                    std::size_t pivot = column;
                    while (
                        pivot < 3U
                        && rows[pivot][column] == 0
                    ) {
                        ++pivot;
                    }
                    if (pivot == 3U) {
                        invertible = false;
                        break;
                    }
                    if (pivot != column) {
                        std::swap(rows[pivot], rows[column]);
                    }
                    const Rational pivot_value = rows[column][column];
                    for (std::size_t entry = column;
                         entry < 7U;
                         ++entry) {
                        rows[column][entry] /= pivot_value;
                    }
                    for (std::size_t row = 0U; row < 3U; ++row) {
                        if (row == column) {
                            continue;
                        }
                        const Rational factor = rows[row][column];
                        for (std::size_t entry = column;
                             entry < 7U;
                             ++entry) {
                            rows[row][entry] -=
                                factor * rows[column][entry];
                        }
                    }
                }
                if (!invertible) {
                    continue;
                }
                std::array<Polynomial, 3> original_variables{
                    Polynomial(),
                    Polynomial(),
                    Polynomial()
                };
                for (std::size_t variable = 0U;
                     variable < 3U;
                     ++variable) {
                    original_variables[variable] =
                        Polynomial(rows[variable][3U]);
                    for (std::size_t basis = 0U;
                         basis < 3U;
                         ++basis) {
                        original_variables[variable] +=
                            basis_variables[basis]
                            * Polynomial(
                                rows[variable][basis + 4U]
                            );
                    }
                }
                Polynomial pulled_back = substitute(
                    chamber.margin,
                    original_variables
                );
                if (zero_coordinate_count > 0) {
                    pulled_back = substitute(
                        pulled_back,
                        std::array<Polynomial, 3>{
                            zero_coordinate_count >= 1
                                ? constant(0) : x0,
                            zero_coordinate_count >= 2
                                ? constant(0) : x1,
                            zero_coordinate_count >= 3
                                ? constant(0) : x2
                        }
                    );
                }
                bool nonnegative = true;
                std::size_t positive = 0U;
                for (const auto& [exponent, coefficient]
                     : pulled_back.terms()) {
                    static_cast<void>(exponent);
                    if (coefficient < 0) {
                        nonnegative = false;
                        break;
                    }
                    if (coefficient > 0) {
                        ++positive;
                    }
                }
                std::string basis_name = "monomial";
                if (!nonnegative) {
                    bool integral_slacks = true;
                    for (const std::size_t index : indices) {
                        const auto affine = affine_coefficients(
                            constraints[index]
                        );
                        for (const Rational& coefficient : affine) {
                            if (coefficient.denominator() != 1) {
                                integral_slacks = false;
                            }
                        }
                    }
                    if (
                        integral_slacks
                        && nonnegative_newton_coefficients(
                            pulled_back,
                            positive
                        )
                    ) {
                        nonnegative = true;
                        basis_name = "newton";
                    }
                }
                if (nonnegative) {
                    if (emit) {
                        std::cout
                            << "SU2_K3_INTERMEDIATE_DIRECT_BASIS"
                            << " mask=" << chamber.mask
                            << " facets=(" << first
                            << ',' << second
                        << ',' << third << ')'
                        << " fixed_coordinates="
                        << zero_coordinate_count
                            << " basis=" << basis_name
                            << " positive_coefficients=" << positive
                            << " result=PASS_EXACT_IDENTITY"
                            << std::endl;
                    }
                    return true;
                }
            }
        }
    }
    return false;
}

Rational evaluate(
    const Polynomial& polynomial,
    const std::array<Integer, 3>& values
) {
    Rational result(0);
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        Rational term = coefficient;
        for (std::size_t index = 0U;
             index < exponent.size();
             ++index) {
            Integer factor = 1;
            for (int power_index = 0;
                 power_index < exponent[index];
                 ++power_index) {
                factor *= values[index];
            }
            term *= Rational(factor);
        }
        result += term;
    }
    return result;
}

std::uint64_t check_integer_points_through_q(
    const Chamber& chamber,
    const std::vector<Hinge>& hinges,
    int maximum_q
) {
    std::uint64_t points = 0U;
    for (int q_value = 1; q_value <= maximum_q; ++q_value) {
        for (int h_value = 1;
             h_value <= 5 * q_value - 2;
             ++h_value) {
            for (int y_value = 0;
                 y_value <= 2 * q_value + 1 + h_value;
                 ++y_value) {
                const std::array<Integer, 3> values{
                    Integer(q_value),
                    Integer(h_value),
                    Integer(y_value)
                };
                std::uint64_t mask = 0U;
                for (std::size_t index = 0U;
                     index < hinges.size();
                     ++index) {
                    if (
                        evaluate(hinges[index].top, values)
                        >= hinges[index].order
                    ) {
                        mask |= std::uint64_t{1} << index;
                    }
                }
                if (mask != chamber.mask) {
                    continue;
                }
                ++points;
                if (evaluate(chamber.margin, values) < 0) {
                    throw std::runtime_error(
                        "finite chamber check has a negative point"
                    );
                }
            }
        }
    }
    return points;
}

bool bounded_integer_certificate(
    const Chamber& chamber,
    const std::vector<Hinge>& hinges,
    bool emit = true
) {
    z3::context context;
    z3::solver solver(context);
    std::vector<z3::expr> multipliers;
    multipliers.reserve(chamber.constraints.size());
    for (std::size_t index = 0U;
         index < chamber.constraints.size();
         ++index) {
        const std::string name =
            "bound_" + std::to_string(chamber.mask)
            + '_' + std::to_string(index);
        multipliers.push_back(context.real_const(name.c_str()));
        solver.add(multipliers.back() >= 0);
    }
    for (std::size_t variable = 0U; variable < 3U; ++variable) {
        z3::expr coefficient_sum = context.real_val(0);
        for (std::size_t index = 0U;
             index < chamber.constraints.size();
             ++index) {
            const auto affine = affine_coefficients(
                chamber.constraints[index]
            );
            coefficient_sum = coefficient_sum
                + rational_value(context, affine[variable + 1U])
                    * multipliers[index];
        }
        solver.add(
            coefficient_sum
            == context.real_val(variable == 0U ? -1 : 0)
        );
    }
    if (solver.check() != z3::sat) {
        return false;
    }

    const z3::model model = solver.get_model();
    Polynomial reconstructed;
    std::size_t nonzero_multipliers = 0U;
    for (std::size_t index = 0U;
         index < chamber.constraints.size();
         ++index) {
        const Rational multiplier = model_rational(
            model.eval(multipliers[index], true)
        );
        if (multiplier < 0) {
            throw std::runtime_error(
                "negative chamber-bound multiplier"
            );
        }
        if (multiplier != 0) {
            ++nonzero_multipliers;
            reconstructed +=
                chamber.constraints[index] * Polynomial(multiplier);
        }
    }
    const auto affine = affine_coefficients(reconstructed);
    if (
        affine[1] != -1
        || affine[2] != 0
        || affine[3] != 0
    ) {
        throw std::runtime_error(
            "chamber-bound model fails exact slope replay"
        );
    }
    const Rational bound = affine[0];
    if (bound < 1) {
        throw std::runtime_error(
            "feasible chamber received a bound below one"
        );
    }
    const Integer maximum_q_integer =
        bound.numerator() / bound.denominator();
    if (
        maximum_q_integer
        > Integer(std::numeric_limits<int>::max())
    ) {
        throw std::runtime_error(
            "finite chamber bound exceeds int range"
        );
    }
    const int maximum_q = maximum_q_integer.convert_to<int>();

    const std::uint64_t points = check_integer_points_through_q(
        chamber,
        hinges,
        maximum_q
    );
    if (points == 0U) {
        throw std::runtime_error(
            "bounded feasible chamber has no enumerated point"
        );
    }
    if (emit) {
        std::cout
            << "SU2_K3_INTERMEDIATE_BOUNDED"
            << " mask=" << chamber.mask
            << " rational_Q_bound=" << bound
            << " integer_Q_maximum=" << maximum_q
            << " nonzero_multipliers=" << nonzero_multipliers
            << " points=" << points
            << " result=PASS_EXACT_FINITE"
            << std::endl;
    }
    return true;
}

bool q_threshold_certificate(
    const Chamber& chamber,
    const std::vector<Hinge>& hinges,
    const std::vector<Polynomial>& base_constraints
) {
    const Polynomial q = Polynomial::variable(0);
    for (int threshold = 2; threshold <= 100; ++threshold) {
        std::vector<Polynomial> tail_constraints = base_constraints;
        tail_constraints.push_back(q - constant(threshold));
        if (!direct_basis_certificate(chamber, tail_constraints)) {
            continue;
        }
        const std::uint64_t lower_points =
            check_integer_points_through_q(
                chamber,
                hinges,
                threshold - 1
            );
        std::cout
            << "SU2_K3_INTERMEDIATE_Q_THRESHOLD"
            << " mask=" << chamber.mask
            << " threshold=" << threshold
            << " lower_points=" << lower_points
            << " result=PASS_EXACT_DICHOTOMY"
            << std::endl;
        return true;
    }
    return false;
}

bool integer_feasible_constraints(
    const std::vector<Polynomial>& constraints
);

bool slack_strip_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& base_constraints
) {
    for (std::size_t facet = 0U;
         facet < base_constraints.size();
         ++facet) {
        for (int threshold = 1; threshold <= 30; ++threshold) {
            std::vector<Polynomial> tail_constraints =
                base_constraints;
            tail_constraints.push_back(
                base_constraints[facet] - constant(threshold)
            );
            if (
                !direct_basis_certificate(
                    chamber,
                    tail_constraints,
                    false
                )
                && integer_feasible_constraints(tail_constraints)
            ) {
                continue;
            }
            bool boundary_passed = true;
            for (int value = 0; value < threshold; ++value) {
                std::vector<Polynomial> plane_constraints;
                plane_constraints.reserve(
                    base_constraints.size() + 1U
                );
                plane_constraints.push_back(
                    base_constraints[facet] - constant(value)
                );
                for (std::size_t index = 0U;
                     index < base_constraints.size();
                     ++index) {
                    if (index != facet) {
                        plane_constraints.push_back(
                            base_constraints[index]
                        );
                    }
                }
                if (
                    !direct_basis_certificate(
                        chamber,
                        plane_constraints,
                        false,
                        1
                    )
                ) {
                    std::vector<Polynomial> exact_plane =
                        base_constraints;
                    exact_plane.push_back(
                        base_constraints[facet] - constant(value)
                    );
                    exact_plane.push_back(
                        constant(value) - base_constraints[facet]
                    );
                    if (
                        integer_feasible_constraints(exact_plane)
                    ) {
                        boundary_passed = false;
                        break;
                    }
                }
            }
            if (!boundary_passed) {
                continue;
            }
            std::cout
                << "SU2_K3_INTERMEDIATE_SLACK_STRIP"
                << " mask=" << chamber.mask
                << " facet=" << facet
                << " threshold=" << threshold
                << " boundary_planes=" << threshold
                << " result=PASS_EXACT_INTEGER_PARTITION"
                << std::endl;
            return true;
        }
    }
    return false;
}

bool constant_sum_compositions_recursive(
    const Chamber& chamber,
    const std::vector<Polynomial>& base_constraints,
    const std::vector<std::size_t>& subset,
    const std::vector<long>& weights,
    long remaining,
    std::size_t position,
    std::vector<long>& values,
    std::uint64_t& slices
) {
    if (position + 1U == subset.size()) {
        if ((remaining % weights[position]) != 0L) {
            return true;
        }
        values[position] = remaining / weights[position];
        ++slices;
        const int fixed_count =
            static_cast<int>(subset.size() - 1U);
        std::vector<Polynomial> slice_constraints;
        slice_constraints.reserve(
            base_constraints.size()
            + static_cast<std::size_t>(fixed_count)
        );
        for (int index = 0; index < fixed_count; ++index) {
            slice_constraints.push_back(
                base_constraints[
                    subset[static_cast<std::size_t>(index)]
                ] - constant(values[static_cast<std::size_t>(index)])
            );
        }
        for (std::size_t index = 0U;
             index < base_constraints.size();
             ++index) {
            bool fixed = false;
            for (int fixed_index = 0;
                 fixed_index < fixed_count;
                 ++fixed_index) {
                if (
                    subset[static_cast<std::size_t>(fixed_index)]
                    == index
                ) {
                    fixed = true;
                }
            }
            if (!fixed) {
                slice_constraints.push_back(base_constraints[index]);
            }
        }
        if (
            direct_basis_certificate(
                chamber,
                slice_constraints,
                false,
                fixed_count
            )
        ) {
            return true;
        }
        std::vector<Polynomial> exact_slice = base_constraints;
        for (int index = 0; index < fixed_count; ++index) {
            const Polynomial equality =
                base_constraints[
                    subset[static_cast<std::size_t>(index)]
                ] - constant(values[static_cast<std::size_t>(index)]);
            exact_slice.push_back(equality);
            exact_slice.push_back(
                constant(0) - equality
            );
        }
        return !integer_feasible_constraints(exact_slice);
    }
    for (
        long value = 0;
        value * weights[position] <= remaining;
        ++value
    ) {
        values[position] = value;
        if (
            !constant_sum_compositions_recursive(
                chamber,
                base_constraints,
                subset,
                weights,
                remaining - weights[position] * value,
                position + 1U,
                values,
                slices
            )
        ) {
            return false;
        }
    }
    return true;
}

bool constant_sum_partition_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& base_constraints
) {
    std::uint64_t combination_count = 1U;
    for (std::size_t index = 0U;
         index < base_constraints.size();
         ++index) {
        combination_count *= 4U;
    }
    for (std::uint64_t code = 1U;
         code < combination_count;
         ++code) {
        std::uint64_t digits = code;
        std::vector<std::size_t> subset;
        std::vector<long> weights;
        Polynomial sum;
        long weight_gcd = 0L;
        for (std::size_t index = 0U;
             index < base_constraints.size();
             ++index) {
            const long weight =
                static_cast<long>(digits % 4U);
            digits /= 4U;
            if (weight != 0L) {
                subset.push_back(index);
                weights.push_back(weight);
                weight_gcd = std::gcd(weight_gcd, weight);
                sum += scale(base_constraints[index], weight);
            }
        }
        if (
            subset.size() < 2U
            || subset.size() > 4U
            || weight_gcd != 1L
        ) {
            continue;
        }
        const auto affine = affine_coefficients(sum);
        if (
            affine[1] != 0
            || affine[2] != 0
            || affine[3] != 0
            || affine[0].denominator() != 1
        ) {
            continue;
        }
        const Integer constant_value = affine[0].numerator();
        if (constant_value < 0 || constant_value > 30) {
            continue;
        }
        const long total = constant_value.convert_to<long>();
        std::vector<long> values(subset.size(), 0L);
        std::uint64_t slices = 0U;
        if (
            constant_sum_compositions_recursive(
                chamber,
                base_constraints,
                subset,
                weights,
                total,
                0U,
                values,
                slices
            )
        ) {
            std::cout
                << "SU2_K3_INTERMEDIATE_CONSTANT_SUM"
                << " mask=" << chamber.mask
                << " subset_size=" << subset.size()
                << " maximum_weight="
                << *std::max_element(weights.begin(), weights.end())
                << " total=" << total
                << " slices=" << slices
                << " result=PASS_EXACT_INTEGER_PARTITION"
                << std::endl;
            return true;
        }
    }
    return false;
}

bool pair_cut_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& base_constraints
) {
    for (std::size_t first = 0U;
         first < base_constraints.size();
         ++first) {
        for (std::size_t second = first + 1U;
             second < base_constraints.size();
             ++second) {
            std::vector<Polynomial> first_branch = base_constraints;
            first_branch.push_back(
                base_constraints[first] - base_constraints[second]
            );
            if (!direct_basis_certificate(chamber, first_branch)) {
                continue;
            }
            std::vector<Polynomial> second_branch = base_constraints;
            second_branch.push_back(
                base_constraints[second] - base_constraints[first]
                - constant(1)
            );
            if (!direct_basis_certificate(chamber, second_branch)) {
                continue;
            }
            std::cout
                << "SU2_K3_INTERMEDIATE_PAIR_CUT"
                << " mask=" << chamber.mask
                << " pair=(" << first << ',' << second << ')'
                << " result=PASS_EXACT_INTEGER_DICHOTOMY"
                << std::endl;
            return true;
        }
    }
    return false;
}

bool minimum_slack_partition_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& base_constraints
) {
    for (std::size_t minimum_index = 0U;
         minimum_index < base_constraints.size();
         ++minimum_index) {
        std::vector<Polynomial> branch = base_constraints;
        for (std::size_t other = 0U;
             other < base_constraints.size();
             ++other) {
            if (other == minimum_index) {
                continue;
            }
            branch.push_back(
                base_constraints[other]
                - base_constraints[minimum_index]
                - constant(other < minimum_index ? 1 : 0)
            );
        }
        if (!direct_basis_certificate(chamber, branch, false)) {
            return false;
        }
    }
    std::cout
        << "SU2_K3_INTERMEDIATE_MINIMUM_SLACK"
        << " mask=" << chamber.mask
        << " branches=" << base_constraints.size()
        << " result=PASS_EXACT_INTEGER_PARTITION"
        << std::endl;
    return true;
}

long floor_divide(long numerator, long denominator) {
    long quotient = numerator / denominator;
    const long remainder = numerator % denominator;
    if (remainder != 0 && numerator < 0) {
        --quotient;
    }
    return quotient;
}

std::vector<Polynomial> rank_one_cg_cuts(
    const std::vector<Polynomial>& base_constraints,
    int maximum_denominator
) {
    std::vector<std::array<long, 4>> affine_constraints;
    affine_constraints.reserve(base_constraints.size());
    for (const Polynomial& constraint : base_constraints) {
        const auto affine = affine_coefficients(constraint);
        std::array<long, 4> integral{};
        for (std::size_t index = 0U; index < 4U; ++index) {
            if (affine[index].denominator() != 1) {
                throw std::runtime_error(
                    "CG source constraint is not integral"
                );
            }
            integral[index] =
                affine[index].numerator().convert_to<long>();
        }
        affine_constraints.push_back(integral);
    }

    std::set<Polynomial> known(
        base_constraints.begin(),
        base_constraints.end()
    );
    std::vector<Polynomial> cuts;
    const Polynomial q = Polynomial::variable(0);
    const Polynomial h = Polynomial::variable(1);
    const Polynomial y = Polynomial::variable(2);
    const std::array<Polynomial, 3> variables{q, h, y};
    for (int denominator = 2;
         denominator <= maximum_denominator;
         ++denominator) {
        std::uint64_t combination_count = 1U;
        for (std::size_t index = 0U;
             index < base_constraints.size();
             ++index) {
            combination_count *=
                static_cast<std::uint64_t>(denominator);
        }
        for (std::uint64_t code = 1U;
             code < combination_count;
             ++code) {
            std::uint64_t digits = code;
            std::array<long, 4> sum{0L, 0L, 0L, 0L};
            for (std::size_t constraint = 0U;
                 constraint < base_constraints.size();
                 ++constraint) {
                const long multiplier = static_cast<long>(
                    digits
                    % static_cast<std::uint64_t>(denominator)
                );
                digits /=
                    static_cast<std::uint64_t>(denominator);
                for (std::size_t index = 0U; index < 4U; ++index) {
                    sum[index] +=
                        multiplier
                        * affine_constraints[constraint][index];
                }
            }
            if (
                (sum[1] % denominator) != 0L
                || (sum[2] % denominator) != 0L
                || (sum[3] % denominator) != 0L
            ) {
                continue;
            }
            Polynomial cut = constant(
                floor_divide(sum[0], denominator)
            );
            for (std::size_t variable = 0U;
                 variable < 3U;
                 ++variable) {
                cut += scale(
                    variables[variable],
                    sum[variable + 1U] / denominator
                );
            }
            if (cut.terms().empty()) {
                continue;
            }
            if (known.insert(cut).second) {
                cuts.push_back(cut);
            }
        }
    }
    return cuts;
}

bool cg_direct_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& base_constraints
) {
    constexpr int maximum_denominator = 3;
    const std::vector<Polynomial> cuts =
        rank_one_cg_cuts(
            base_constraints,
            maximum_denominator
        );
    std::vector<Polynomial> enriched = base_constraints;
    enriched.insert(enriched.end(), cuts.begin(), cuts.end());
    if (!direct_basis_certificate(chamber, enriched, false)) {
        return false;
    }
    std::cout
        << "SU2_K3_INTERMEDIATE_CG_DIRECT"
        << " mask=" << chamber.mask
        << " base_facets=" << base_constraints.size()
        << " maximum_denominator=" << maximum_denominator
        << " cg_cuts=" << cuts.size()
        << " result=PASS_EXACT_LATTICE_IDENTITY"
        << std::endl;
    return true;
}

bool integer_feasible_constraints(
    const std::vector<Polynomial>& constraints
) {
    z3::context context;
    z3::solver solver(context);
    const std::array<z3::expr, 3> variables{
        context.int_const("Q"),
        context.int_const("H"),
        context.int_const("Y")
    };
    for (const Polynomial& constraint : constraints) {
        solver.add(
            affine_integer_value(context, constraint, variables) >= 0
        );
    }
    return solver.check() == z3::sat;
}

struct CutTreeStatistics {
    std::uint64_t nodes = 0U;
    std::uint64_t direct_leaves = 0U;
    std::uint64_t empty_leaves = 0U;
};

bool pair_cut_tree_recursive(
    const Chamber& chamber,
    const std::vector<Polynomial>& base_constraints,
    const std::vector<Polynomial>& current_constraints,
    int remaining_depth,
    std::vector<bool>& used_pairs,
    CutTreeStatistics& statistics
) {
    ++statistics.nodes;
    if (
        direct_basis_certificate(
            chamber,
            current_constraints,
            false
        )
    ) {
        ++statistics.direct_leaves;
        return true;
    }
    if (remaining_depth == 0) {
        if (!integer_feasible_constraints(current_constraints)) {
            ++statistics.empty_leaves;
            return true;
        }
        return false;
    }

    std::size_t pair_index = 0U;
    for (std::size_t first = 0U;
         first < base_constraints.size();
         ++first) {
        for (std::size_t second = first + 1U;
             second < base_constraints.size();
             ++second, ++pair_index) {
            if (used_pairs[pair_index]) {
                continue;
            }
            used_pairs[pair_index] = true;
            std::vector<Polynomial> first_branch =
                current_constraints;
            first_branch.push_back(
                base_constraints[first] - base_constraints[second]
            );
            CutTreeStatistics first_statistics = statistics;
            if (
                pair_cut_tree_recursive(
                    chamber,
                    base_constraints,
                    first_branch,
                    remaining_depth - 1,
                    used_pairs,
                    first_statistics
                )
            ) {
                std::vector<Polynomial> second_branch =
                    current_constraints;
                second_branch.push_back(
                    base_constraints[second]
                    - base_constraints[first]
                    - constant(1)
                );
                CutTreeStatistics second_statistics =
                    first_statistics;
                if (
                    pair_cut_tree_recursive(
                        chamber,
                        base_constraints,
                        second_branch,
                        remaining_depth - 1,
                        used_pairs,
                        second_statistics
                    )
                ) {
                    statistics = second_statistics;
                    used_pairs[pair_index] = false;
                    return true;
                }
            }
            used_pairs[pair_index] = false;
        }
    }
    return false;
}

bool pair_cut_tree_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& base_constraints,
    int maximum_depth
) {
    const std::size_t pair_count =
        base_constraints.size()
            * (base_constraints.size() - 1U) / 2U;
    std::vector<bool> used_pairs(pair_count, false);
    CutTreeStatistics statistics;
    if (
        !pair_cut_tree_recursive(
            chamber,
            base_constraints,
            base_constraints,
            maximum_depth,
            used_pairs,
            statistics
        )
    ) {
        return false;
    }
    std::cout
        << "SU2_K3_INTERMEDIATE_PAIR_TREE"
        << " mask=" << chamber.mask
        << " depth=" << maximum_depth
        << " nodes=" << statistics.nodes
        << " direct_leaves=" << statistics.direct_leaves
        << " empty_leaves=" << statistics.empty_leaves
        << " result=PASS_EXACT_INTEGER_PARTITION"
        << std::endl;
    return true;
}

bool handelman_feasible(
    const Chamber& chamber,
    const std::vector<Hinge>& hinges,
    int maximum_degree
) {
    const std::vector<Polynomial> constraints =
        irredundant_constraints(chamber);
    if (direct_basis_certificate(chamber, constraints)) {
        return true;
    }
    if (bounded_integer_certificate(chamber, hinges)) {
        return true;
    }
    if (pair_cut_certificate(chamber, constraints)) {
        return true;
    }
    if (constant_sum_partition_certificate(chamber, constraints)) {
        return true;
    }
    if (minimum_slack_partition_certificate(chamber, constraints)) {
        return true;
    }
    if (cg_direct_certificate(chamber, constraints)) {
        return true;
    }
    if (slack_strip_certificate(chamber, constraints)) {
        return true;
    }
    if (pair_cut_tree_certificate(chamber, constraints, 2)) {
        return true;
    }
    if (q_threshold_certificate(chamber, hinges, constraints)) {
        return true;
    }
    std::map<Polynomial, Polynomial> unique_products;
    enumerate_products_recursive(
        constraints,
        maximum_degree,
        0U,
        constant(1),
        unique_products
    );
    std::vector<Polynomial> products;
    products.reserve(unique_products.size());
    for (const auto& [key, value] : unique_products) {
        static_cast<void>(key);
        products.push_back(value);
    }

    std::cout
        << "SU2_K3_INTERMEDIATE_PRODUCTS"
        << " mask=" << chamber.mask
        << " facets=" << constraints.size()
        << " degree=" << maximum_degree
        << " products=" << products.size()
        << " result=BUILT"
        << std::endl;

    z3::context context;
    z3::solver solver(context);
    std::vector<z3::expr> coefficients;
    coefficients.reserve(products.size());
    for (std::size_t index = 0U; index < products.size(); ++index) {
        const std::string name =
            "c_" + std::to_string(chamber.mask)
            + '_' + std::to_string(index);
        coefficients.push_back(context.real_const(name.c_str()));
        solver.add(coefficients.back() >= 0);
    }

    for (int q_power = 0; q_power <= maximum_degree; ++q_power) {
        for (int h_power = 0;
             h_power <= maximum_degree - q_power;
             ++h_power) {
            for (int y_power = 0;
                 y_power <= maximum_degree - q_power - h_power;
                 ++y_power) {
                const Exponent exponent{q_power, h_power, y_power};
                z3::expr left = context.real_val(0);
                for (std::size_t index = 0U;
                     index < products.size();
                     ++index) {
                    const auto position =
                        products[index].terms().find(exponent);
                    if (position != products[index].terms().end()) {
                        left = left
                            + rational_value(context, position->second)
                                * coefficients[index];
                    }
                }
                const auto target_position =
                    chamber.margin.terms().find(exponent);
                const Rational target =
                    target_position == chamber.margin.terms().end()
                    ? Rational(0)
                    : target_position->second;
                solver.add(left == rational_value(context, target));
            }
        }
    }

    const z3::check_result result = solver.check();
    std::size_t nonzero_coefficients = 0U;
    if (result == z3::sat) {
        const z3::model model = solver.get_model();
        Polynomial reconstructed;
        for (std::size_t index = 0U; index < products.size(); ++index) {
            const Rational coefficient = model_rational(
                model.eval(coefficients[index], true)
            );
            if (coefficient < 0) {
                throw std::runtime_error(
                    "negative Handelman model coefficient"
                );
            }
            if (coefficient != 0) {
                ++nonzero_coefficients;
                reconstructed +=
                    products[index] * Polynomial(coefficient);
            }
        }
        if (reconstructed.terms() != chamber.margin.terms()) {
            throw std::runtime_error(
                "Handelman model fails exact polynomial replay"
            );
        }
    }

    std::cout
        << "SU2_K3_INTERMEDIATE_HANDELMAN"
        << " mask=" << chamber.mask
        << " facets=" << constraints.size()
        << " products=" << products.size()
        << " nonzero_coefficients=" << nonzero_coefficients
        << " result="
        << (
            result == z3::sat
                ? "PASS_EXACT_IDENTITY"
                : result == z3::unsat ? "UNSAT" : "UNKNOWN"
        )
        << std::endl;
    return result == z3::sat;
}

std::vector<std::uint64_t> feasible_masks(
    const std::vector<Hinge>& hinges
) {
    z3::context context;
    z3::solver solver(context);
    const z3::expr q = context.int_const("Q");
    const z3::expr h = context.int_const("H");
    const z3::expr y = context.int_const("Y");
    solver.add(q >= 1);
    solver.add(h >= 1);
    solver.add(h <= 5 * q - 2);
    solver.add(y >= 0);
    solver.add(y <= 2 * q + 1 + h);
    const std::array<z3::expr, 3> variables{q, h, y};

    std::vector<std::uint64_t> masks;
    while (solver.check() == z3::sat) {
        const z3::model model = solver.get_model();
        std::uint64_t mask = 0U;
        z3::expr block = context.bool_val(false);
        for (std::size_t index = 0U; index < hinges.size(); ++index) {
            const z3::expr top = affine_integer_value(
                context,
                hinges[index].top,
                variables
            );
            const z3::expr active_expression =
                top >= hinges[index].order;
            const bool active =
                model.eval(active_expression, true).bool_value()
                == Z3_L_TRUE;
            if (active) {
                mask |= std::uint64_t{1} << index;
                block = block || !active_expression;
            } else {
                block = block || active_expression;
            }
        }
        masks.push_back(mask);
        solver.add(block);
    }
    std::sort(masks.begin(), masks.end());
    const auto unique_end = std::unique(masks.begin(), masks.end());
    if (unique_end != masks.end()) {
        throw std::runtime_error(
            "activation-mask blocking produced a duplicate"
        );
    }
    return masks;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Polynomial q = Polynomial::variable(0);
        const Polynomial h = Polynomial::variable(1);
        const Polynomial y = Polynomial::variable(2);
        const std::vector<Hinge> hinges{
            {"u6_m1", scale(q, 6) - y + constant(4), 4},
            {"u6_m2", scale(q, 4) - y + constant(3), 4},
            {"u6_m3", scale(q, 2) - y + constant(2), 4},
            {
                "u6_image1",
                scale(q, 2) - scale(h, 2) + y + constant(1),
                4
            },
            {"u6_image2", y - scale(h, 2), 4},
            {
                "u6_image3",
                scale(q, 2) - scale(h, 2) - y,
                4
            },
            {"u7_m1", scale(q, 7) - y + constant(5), 5},
            {"u7_m2", scale(q, 5) - y + constant(4), 5},
            {"u7_m3", scale(q, 3) - y + constant(3), 5},
            {"u7_m4", q - y + constant(2), 5},
            {
                "u7_image1",
                scale(q, 3) - scale(h, 2) + y + constant(2),
                5
            },
            {
                "u7_image2",
                q - scale(h, 2) + y + constant(1),
                5
            },
            {
                "u7_image3",
                y - q - scale(h, 2),
                5
            },
            {
                "u7_image4",
                scale(q, 3) - scale(h, 2) - y + constant(1),
                5
            },
            {"u7_image5", q - scale(h, 2) - y, 5},
            {
                "u7_image6",
                y - q - scale(h, 4) - constant(2),
                5
            },
            {
                "f6_image1",
                scale(q, 2) - scale(h, 2) + constant(1),
                4
            },
            {
                "f7_image1",
                scale(q, 3) - scale(h, 2) + constant(2),
                5
            },
            {
                "f7_image2",
                q - scale(h, 2) + constant(1),
                5
            },
            {
                "f7_image3",
                scale(q, 3) - scale(h, 2) + constant(1),
                5
            },
            {"f7_image4", q - scale(h, 2), 5},
            {"f7_m4", q + constant(2), 5}
        };
        if (hinges.size() > 63U) {
            throw std::runtime_error("activation mask exceeds uint64");
        }

        bool patterns_only = false;
        bool direct_only = false;
        bool elementary_only = false;
        bool threshold_only = false;
        bool pair_only = false;
        bool tree_only = false;
        bool minimum_only = false;
        bool cg_only = false;
        bool strip_only = false;
        bool sum_only = false;
        bool describe_only = false;
        std::uint64_t selected_mask = 0U;
        bool has_selected_mask = false;
        int maximum_degree = 9;
        if (argc == 2 && std::string(argv[1]) == "--patterns") {
            patterns_only = true;
        } else if (
            argc == 2
            && std::string(argv[1]) == "--direct-only"
        ) {
            direct_only = true;
        } else if (
            argc == 2
            && std::string(argv[1]) == "--elementary-only"
        ) {
            elementary_only = true;
        } else if (
            argc == 2
            && std::string(argv[1]) == "--threshold-only"
        ) {
            threshold_only = true;
        } else if (
            argc == 2
            && std::string(argv[1]) == "--pair-only"
        ) {
            pair_only = true;
        } else if (
            argc == 2
            && std::string(argv[1]) == "--tree-only"
        ) {
            tree_only = true;
        } else if (
            argc == 2
            && std::string(argv[1]) == "--minimum-only"
        ) {
            minimum_only = true;
        } else if (
            argc == 2
            && std::string(argv[1]) == "--cg-only"
        ) {
            cg_only = true;
        } else if (
            argc == 2
            && std::string(argv[1]) == "--strip-only"
        ) {
            strip_only = true;
        } else if (
            argc == 2
            && std::string(argv[1]) == "--sum-only"
        ) {
            sum_only = true;
        } else if (
            (argc == 3 || argc == 5)
            && std::string(argv[1]) == "--mask"
        ) {
            selected_mask = std::stoull(argv[2]);
            has_selected_mask = true;
            if (argc == 5) {
                if (std::string(argv[3]) != "--degree") {
                    throw std::runtime_error("expected --degree");
                }
                maximum_degree = std::stoi(argv[4]);
            }
        } else if (
            argc == 3
            && std::string(argv[1]) == "--describe"
        ) {
            selected_mask = std::stoull(argv[2]);
            has_selected_mask = true;
            describe_only = true;
        } else if (argc != 1) {
            throw std::runtime_error(
                "usage: prove_su2_k3_intermediate [--patterns] "
                "[--direct-only] [--elementary-only] "
                "[--threshold-only] [--pair-only] [--tree-only] "
                "[--minimum-only] [--cg-only] [--strip-only] "
                "[--sum-only] "
                "[--mask MASK [--degree DEGREE]] "
                "[--describe MASK]"
            );
        }
        if (maximum_degree < 9 || maximum_degree > 12) {
            throw std::runtime_error(
                "Handelman degree must lie from 9 through 12"
            );
        }

        const std::vector<std::uint64_t> masks = feasible_masks(hinges);
        if (patterns_only) {
            std::cout
                << "SU2_K3_INTERMEDIATE_PATTERNS"
                << " hinges=" << hinges.size()
                << " feasible_chambers=" << masks.size()
                << " result=PASS_EXACT_ENUMERATION\n";
            return EXIT_SUCCESS;
        }

        std::size_t attempted = 0U;
        std::size_t certified = 0U;
        for (const std::uint64_t mask : masks) {
            if (has_selected_mask && mask != selected_mask) {
                continue;
            }
            ++attempted;
            const Chamber chamber = make_chamber(mask, hinges);
            const std::vector<Polynomial> constraints =
                irredundant_constraints(chamber);
            if (describe_only) {
                for (std::size_t index = 0U;
                     index < constraints.size();
                     ++index) {
                    const auto affine =
                        affine_coefficients(constraints[index]);
                    std::cout
                        << "SU2_K3_INTERMEDIATE_FACET"
                        << " mask=" << mask
                        << " index=" << index
                        << " affine=("
                        << affine[0] << ','
                        << affine[1] << ','
                        << affine[2] << ','
                        << affine[3] << ")\n";
                }
                ++certified;
                continue;
            }
            std::cout
                << "SU2_K3_INTERMEDIATE_CHAMBER"
                << " mask=" << mask
                << " result=FEASIBLE"
                << std::endl;
            const bool passed = direct_only
                ? direct_basis_certificate(
                    chamber,
                    constraints
                )
                : elementary_only
                    ? (
                        direct_basis_certificate(
                            chamber,
                            constraints
                        )
                        || bounded_integer_certificate(chamber, hinges)
                    )
                    : threshold_only
                        ? (
                            direct_basis_certificate(
                                chamber,
                                constraints
                            )
                            || bounded_integer_certificate(
                                chamber,
                                hinges
                            )
                            || q_threshold_certificate(
                                chamber,
                                hinges,
                                constraints
                            )
                        )
                    : pair_only
                        ? (
                            direct_basis_certificate(
                                chamber,
                                constraints
                            )
                            || bounded_integer_certificate(
                                chamber,
                                hinges
                            )
                            || pair_cut_certificate(
                                chamber,
                                constraints
                            )
                        )
                    : tree_only
                        ? (
                            direct_basis_certificate(
                                chamber,
                                constraints
                            )
                            || bounded_integer_certificate(
                                chamber,
                                hinges
                            )
                            || pair_cut_certificate(
                                chamber,
                                constraints
                            )
                            || pair_cut_tree_certificate(
                                chamber,
                                constraints,
                                2
                            )
                        )
                    : minimum_only
                        ? (
                            direct_basis_certificate(
                                chamber,
                                constraints
                            )
                            || bounded_integer_certificate(
                                chamber,
                                hinges
                            )
                            || pair_cut_certificate(
                                chamber,
                                constraints
                            )
                            || minimum_slack_partition_certificate(
                                chamber,
                                constraints
                            )
                        )
                    : cg_only
                        ? (
                            direct_basis_certificate(
                                chamber,
                                constraints
                            )
                            || bounded_integer_certificate(
                                chamber,
                                hinges
                            )
                            || pair_cut_certificate(
                                chamber,
                                constraints
                            )
                            || cg_direct_certificate(
                                chamber,
                                constraints
                            )
                        )
                    : strip_only
                        ? (
                            direct_basis_certificate(
                                chamber,
                                constraints
                            )
                            || bounded_integer_certificate(
                                chamber,
                                hinges
                            )
                            || pair_cut_certificate(
                                chamber,
                                constraints
                            )
                            || slack_strip_certificate(
                                chamber,
                                constraints
                            )
                        )
                    : sum_only
                        ? (
                            direct_basis_certificate(
                                chamber,
                                constraints
                            )
                            || bounded_integer_certificate(
                                chamber,
                                hinges
                            )
                            || pair_cut_certificate(
                                chamber,
                                constraints
                            )
                            || constant_sum_partition_certificate(
                                chamber,
                                constraints
                            )
                        )
                    : handelman_feasible(
                        chamber,
                        hinges,
                        maximum_degree
                    );
            if (passed) {
                ++certified;
            }
        }
        if (has_selected_mask && attempted == 0U) {
            throw std::runtime_error(
                "selected mask is not an exact feasible chamber"
            );
        }
        if (describe_only) {
            std::cout
                << "SU2_K3_INTERMEDIATE_DESCRIPTION"
                << " mask=" << selected_mask
                << " result=PASS_EXACT_FACETS\n";
            return EXIT_SUCCESS;
        }
        if (direct_only) {
            std::cout
                << "SU2_K3_INTERMEDIATE_DIRECT_SCAN"
                << " feasible_chambers=" << masks.size()
                << " certified_chambers=" << certified
                << " result=PASS_EXACT_SCAN\n";
            return EXIT_SUCCESS;
        }
        if (elementary_only) {
            std::cout
                << "SU2_K3_INTERMEDIATE_ELEMENTARY_SCAN"
                << " feasible_chambers=" << masks.size()
                << " certified_chambers=" << certified
                << " result=PASS_EXACT_SCAN\n";
            return EXIT_SUCCESS;
        }
        if (threshold_only) {
            std::cout
                << "SU2_K3_INTERMEDIATE_THRESHOLD_SCAN"
                << " feasible_chambers=" << masks.size()
                << " certified_chambers=" << certified
                << " result=PASS_EXACT_SCAN\n";
            return EXIT_SUCCESS;
        }
        if (pair_only) {
            std::cout
                << "SU2_K3_INTERMEDIATE_PAIR_SCAN"
                << " feasible_chambers=" << masks.size()
                << " certified_chambers=" << certified
                << " result=PASS_EXACT_SCAN\n";
            return EXIT_SUCCESS;
        }
        if (tree_only) {
            std::cout
                << "SU2_K3_INTERMEDIATE_TREE_SCAN"
                << " feasible_chambers=" << masks.size()
                << " certified_chambers=" << certified
                << " result=PASS_EXACT_SCAN\n";
            return EXIT_SUCCESS;
        }
        if (minimum_only) {
            std::cout
                << "SU2_K3_INTERMEDIATE_MINIMUM_SCAN"
                << " feasible_chambers=" << masks.size()
                << " certified_chambers=" << certified
                << " result=PASS_EXACT_SCAN\n";
            return EXIT_SUCCESS;
        }
        if (cg_only) {
            std::cout
                << "SU2_K3_INTERMEDIATE_CG_SCAN"
                << " feasible_chambers=" << masks.size()
                << " certified_chambers=" << certified
                << " result=PASS_EXACT_SCAN\n";
            return EXIT_SUCCESS;
        }
        if (strip_only) {
            std::cout
                << "SU2_K3_INTERMEDIATE_STRIP_SCAN"
                << " feasible_chambers=" << masks.size()
                << " certified_chambers=" << certified
                << " result=PASS_EXACT_SCAN\n";
            return EXIT_SUCCESS;
        }
        if (sum_only) {
            std::cout
                << "SU2_K3_INTERMEDIATE_SUM_SCAN"
                << " feasible_chambers=" << masks.size()
                << " certified_chambers=" << certified
                << " result=PASS_EXACT_SCAN\n";
            return EXIT_SUCCESS;
        }
        std::cout
            << "SU2_K3_INTERMEDIATE"
            << " feasible_chambers=" << masks.size()
            << " attempted_chambers=" << attempted
            << " certified_chambers=" << certified
            << " degree=" << maximum_degree
            << " result="
            << (
                attempted == certified
                    ? "PASS_EXACT_CERTIFICATE"
                    : "INCOMPLETE"
            )
            << '\n';
        return attempted == certified ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_K3_INTERMEDIATE FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
