#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
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

    static Polynomial variable(std::size_t index) {
        Polynomial result;
        Exponent exponent{0, 0, 0};
        exponent[index] = 1;
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
                for (std::size_t index = 0; index < 3U; ++index) {
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

Polynomial power(Polynomial value, int exponent) {
    Polynomial result = constant(1);
    while (exponent > 0) {
        if ((exponent & 1) != 0) {
            result *= value;
        }
        exponent /= 2;
        if (exponent != 0) {
            value *= value;
        }
    }
    return result;
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

Polynomial substitute(
    const Polynomial& polynomial,
    const std::array<Polynomial, 3>& values
) {
    Polynomial result;
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        Polynomial term(coefficient);
        for (std::size_t index = 0; index < 3U; ++index) {
            term *= power(values[index], exponent[index]);
        }
        result += term;
    }
    return result;
}

Rational evaluate(
    const Polynomial& polynomial,
    const std::array<int, 3>& values
) {
    Rational result(0);
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        Rational term = coefficient;
        for (std::size_t index = 0; index < 3U; ++index) {
            Integer factor = 1;
            for (int power_index = 0;
                 power_index < exponent[index];
                 ++power_index) {
                factor *= values[index];
            }
            term *= factor;
        }
        result += term;
    }
    return result;
}

using Affine = std::array<Rational, 4>;

Affine affine_coefficients(const Polynomial& polynomial) {
    Affine result{Rational(0), Rational(0), Rational(0), Rational(0)};
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        const int degree = exponent[0] + exponent[1] + exponent[2];
        if (degree > 1) {
            throw std::runtime_error("constraint is not affine");
        }
        if (degree == 0) {
            result[0] += coefficient;
        } else {
            for (std::size_t index = 0; index < 3U; ++index) {
                if (exponent[index] == 1) {
                    result[index + 1U] += coefficient;
                }
            }
        }
    }
    return result;
}

struct Hinge {
    Polynomial slack;
};

struct Term {
    std::size_t hinge;
    long coefficient;
};

struct Formula {
    std::vector<Hinge> hinges;
    std::vector<Term> u8;
    std::vector<Term> u9;
    std::vector<Term> f8;
    std::vector<Term> f9;
};

long binomial_long(int top, int bottom) {
    if (bottom < 0 || bottom > top) {
        return 0;
    }
    long result = 1;
    for (int index = 1; index <= bottom; ++index) {
        result = result * (top - bottom + index) / index;
    }
    return result;
}

void append_block(
    Formula& formula,
    std::vector<Term>& terms,
    const Polynomial& q,
    const Polynomial& image,
    int tensor_power,
    int maximum_image,
    long outer_sign
) {
    for (int index = 0; index <= maximum_image; ++index) {
        const Polynomial slack =
            scale(q, tensor_power - 2 * index)
            - image
            - constant(index);
        const std::size_t hinge = formula.hinges.size();
        formula.hinges.push_back(Hinge{slack});
        const long alternating = index % 2 == 0 ? 1 : -1;
        terms.push_back(Term{
            hinge,
            outer_sign
                * alternating
                * binomial_long(tensor_power, index)
        });
    }
}

Formula make_formula() {
    const Polynomial q = Polynomial::variable(0);
    const Polynomial h = Polynomial::variable(1);
    const Polynomial y = Polynomial::variable(2);
    const Polynomial period =
        scale(q, 4) + scale(h, 2) + constant(4);
    Formula formula;

    const auto append_eight = [&](
        std::vector<Term>& terms,
        const Polynomial& endpoint
    ) {
        append_block(formula, terms, q, endpoint, 8, 3, 1);
        append_block(
            formula,
            terms,
            q,
            period - constant(1) - endpoint,
            8,
            2,
            -1
        );
        append_block(
            formula,
            terms,
            q,
            period + endpoint,
            8,
            1,
            1
        );
        append_block(
            formula,
            terms,
            q,
            scale(period, 2) - constant(1) - endpoint,
            8,
            0,
            -1
        );
    };
    const auto append_nine = [&](
        std::vector<Term>& terms,
        const Polynomial& endpoint
    ) {
        append_block(formula, terms, q, endpoint, 9, 4, 1);
        append_block(
            formula,
            terms,
            q,
            period - constant(1) - endpoint,
            9,
            3,
            -1
        );
        append_block(
            formula,
            terms,
            q,
            period + endpoint,
            9,
            2,
            1
        );
        append_block(
            formula,
            terms,
            q,
            scale(period, 2) - constant(1) - endpoint,
            9,
            1,
            -1
        );
        append_block(
            formula,
            terms,
            q,
            scale(period, 2) + endpoint,
            9,
            0,
            1
        );
    };

    append_eight(formula.u8, y);
    append_nine(formula.u9, y);
    append_eight(formula.f8, constant(0));
    append_nine(formula.f9, constant(0));
    if (formula.hinges.size() != 50U) {
        throw std::runtime_error("unexpected C4 hinge count");
    }
    return formula;
}

Polynomial selected_sum(
    const Formula& formula,
    std::uint64_t mask,
    const std::vector<Term>& terms,
    int order
) {
    Polynomial result;
    for (const Term& term : terms) {
        if ((mask & (std::uint64_t{1} << term.hinge)) != 0U) {
            result += scale(
                binomial(
                    formula.hinges[term.hinge].slack
                        + constant(order),
                    order
                ),
                term.coefficient
            );
        }
    }
    return result;
}

struct Chamber {
    std::uint64_t mask;
    std::vector<Polynomial> constraints;
    Polynomial margin;
};

Chamber make_chamber(
    const Formula& formula,
    std::uint64_t mask
) {
    const Polynomial q = Polynomial::variable(0);
    const Polynomial h = Polynomial::variable(1);
    const Polynomial y = Polynomial::variable(2);
    Chamber chamber{
        mask,
        {
            q - constant(1),
            h,
            y,
            scale(q, 2) + constant(1) + h - y
        },
        Polynomial()
    };
    for (std::size_t index = 0; index < formula.hinges.size(); ++index) {
        const Polynomial& slack = formula.hinges[index].slack;
        chamber.constraints.push_back(
            (mask & (std::uint64_t{1} << index)) != 0U
                ? slack
                : constant(-1) - slack
        );
    }
    const Polynomial u8 = selected_sum(formula, mask, formula.u8, 6);
    const Polynomial u9 = selected_sum(formula, mask, formula.u9, 7);
    const Polynomial f8 = selected_sum(formula, mask, formula.f8, 6);
    const Polynomial f9 = selected_sum(formula, mask, formula.f9, 7);
    chamber.margin = f8 * u9 - f9 * u8;
    return chamber;
}

z3::expr z3_affine(
    z3::context& context,
    const Polynomial& polynomial,
    const std::array<z3::expr, 3>& variables
) {
    const Affine affine = affine_coefficients(polynomial);
    z3::expr result = context.int_val(
        affine[0].numerator().convert_to<std::string>().c_str()
    );
    for (std::size_t index = 0; index < 3U; ++index) {
        if (affine[index + 1U].denominator() != 1) {
            throw std::runtime_error(
                "nonintegral coefficient in Presburger constraint"
            );
        }
        result = result
            + context.int_val(
                affine[index + 1U]
                    .numerator().convert_to<std::string>().c_str()
            ) * variables[index];
    }
    return result;
}

std::vector<std::uint64_t> feasible_masks(const Formula& formula) {
    z3::context context;
    z3::solver solver(context);
    const z3::expr q = context.int_const("Q");
    const z3::expr h = context.int_const("H");
    const z3::expr y = context.int_const("Y");
    const std::array<z3::expr, 3> variables{q, h, y};
    solver.add(q >= 1);
    solver.add(h >= 0);
    solver.add(y >= 0);
    solver.add(y <= 2 * q + 1 + h);

    std::vector<std::uint64_t> masks;
    while (solver.check() == z3::sat) {
        const z3::model model = solver.get_model();
        std::uint64_t mask = 0U;
        z3::expr block = context.bool_val(false);
        for (std::size_t index = 0; index < formula.hinges.size(); ++index) {
            const z3::expr slack = z3_affine(
                context,
                formula.hinges[index].slack,
                variables
            );
            const z3::expr active = slack >= 0;
            const bool value =
                model.eval(active, true).bool_value() == Z3_L_TRUE;
            if (value) {
                mask |= std::uint64_t{1} << index;
            }
            block = block || (value ? !active : active);
        }
        masks.push_back(mask);
        solver.add(block);
    }
    std::sort(masks.begin(), masks.end());
    if (std::unique(masks.begin(), masks.end()) != masks.end()) {
        throw std::runtime_error("duplicate feasible activation mask");
    }
    return masks;
}

std::vector<Polynomial> irredundant_constraints(
    const Chamber& chamber
) {
    std::vector<Polynomial> constraints = chamber.constraints;
    for (std::size_t position = 0; position < constraints.size();) {
        z3::context context;
        z3::solver solver(context);
        const std::array<z3::expr, 3> variables{
            context.int_const("Q"),
            context.int_const("H"),
            context.int_const("Y")
        };
        for (std::size_t other = 0; other < constraints.size(); ++other) {
            if (other != position) {
                solver.add(
                    z3_affine(context, constraints[other], variables) >= 0
                );
            }
        }
        solver.add(
            z3_affine(context, constraints[position], variables) < 0
        );
        if (solver.check() == z3::unsat) {
            constraints.erase(
                constraints.begin()
                    + static_cast<std::ptrdiff_t>(position)
            );
        } else {
            ++position;
        }
    }
    return constraints;
}

Rational determinant(
    const std::array<std::array<Rational, 3>, 3>& matrix
) {
    return
        matrix[0][0]
            * (matrix[1][1] * matrix[2][2]
               - matrix[1][2] * matrix[2][1])
        - matrix[0][1]
            * (matrix[1][0] * matrix[2][2]
               - matrix[1][2] * matrix[2][0])
        + matrix[0][2]
            * (matrix[1][0] * matrix[2][1]
               - matrix[1][1] * matrix[2][0]);
}

std::array<Polynomial, 3> inverse_facet_map(
    const std::array<Polynomial, 3>& facets
) {
    std::array<Affine, 3> affine{
        affine_coefficients(facets[0]),
        affine_coefficients(facets[1]),
        affine_coefficients(facets[2])
    };
    std::array<std::array<Rational, 3>, 3> matrix{};
    for (std::size_t row = 0; row < 3U; ++row) {
        for (std::size_t column = 0; column < 3U; ++column) {
            matrix[row][column] = affine[row][column + 1U];
        }
    }
    const Rational det = determinant(matrix);
    if (det == 0) {
        throw std::runtime_error("singular facet map");
    }
    const Rational a = matrix[0][0];
    const Rational b = matrix[0][1];
    const Rational c = matrix[0][2];
    const Rational d = matrix[1][0];
    const Rational e = matrix[1][1];
    const Rational f = matrix[1][2];
    const Rational g_entry = matrix[2][0];
    const Rational h = matrix[2][1];
    const Rational i = matrix[2][2];
    const std::array<std::array<Rational, 3>, 3> inverse{{
        {{(e * i - f * h) / det,
          (c * h - b * i) / det,
          (b * f - c * e) / det}},
        {{(f * g_entry - d * i) / det,
          (a * i - c * g_entry) / det,
          (c * d - a * f) / det}},
        {{(d * h - e * g_entry) / det,
          (b * g_entry - a * h) / det,
          (a * e - b * d) / det}}
    }};

    const std::array<Polynomial, 3> g{
        Polynomial::variable(0),
        Polynomial::variable(1),
        Polynomial::variable(2)
    };
    std::array<Polynomial, 3> result;
    for (std::size_t row = 0; row < 3U; ++row) {
        for (std::size_t column = 0; column < 3U; ++column) {
            result[row] +=
                (g[column] - Polynomial(affine[column][0]))
                * Polynomial(inverse[row][column]);
        }
    }
    return result;
}

bool nonnegative_monomial(const Polynomial& polynomial) {
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        static_cast<void>(exponent);
        if (coefficient < 0) {
            return false;
        }
    }
    return true;
}

bool nonnegative_newton(const Polynomial& polynomial) {
    Exponent degree{0, 0, 0};
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        static_cast<void>(coefficient);
        for (std::size_t index = 0; index < 3U; ++index) {
            degree[index] = std::max(degree[index], exponent[index]);
        }
    }
    const std::size_t size0 = static_cast<std::size_t>(degree[0] + 1);
    const std::size_t size1 = static_cast<std::size_t>(degree[1] + 1);
    const std::size_t size2 = static_cast<std::size_t>(degree[2] + 1);
    const auto offset = [size1, size2](
        std::size_t first,
        std::size_t second,
        std::size_t third
    ) {
        return (first * size1 + second) * size2 + third;
    };
    std::vector<Rational> table(size0 * size1 * size2);
    for (std::size_t first = 0; first < size0; ++first) {
        for (std::size_t second = 0; second < size1; ++second) {
            for (std::size_t third = 0; third < size2; ++third) {
                table[offset(first, second, third)] = evaluate(
                    polynomial,
                    std::array<int, 3>{
                        static_cast<int>(first),
                        static_cast<int>(second),
                        static_cast<int>(third)
                    }
                );
            }
        }
    }
    for (std::size_t second = 0; second < size1; ++second) {
        for (std::size_t third = 0; third < size2; ++third) {
            for (std::size_t order = 0; order < size0; ++order) {
                for (std::size_t first = size0 - 1U;
                     first > order;
                     --first) {
                    table[offset(first, second, third)] -=
                        table[offset(first - 1U, second, third)];
                }
            }
        }
    }
    for (std::size_t first = 0; first < size0; ++first) {
        for (std::size_t third = 0; third < size2; ++third) {
            for (std::size_t order = 0; order < size1; ++order) {
                for (std::size_t second = size1 - 1U;
                     second > order;
                     --second) {
                    table[offset(first, second, third)] -=
                        table[offset(first, second - 1U, third)];
                }
            }
        }
    }
    for (std::size_t first = 0; first < size0; ++first) {
        for (std::size_t second = 0; second < size1; ++second) {
            for (std::size_t order = 0; order < size2; ++order) {
                for (std::size_t third = size2 - 1U;
                     third > order;
                     --third) {
                    table[offset(first, second, third)] -=
                        table[offset(first, second, third - 1U)];
                }
            }
        }
    }
    return std::all_of(
        table.begin(),
        table.end(),
        [](const Rational& coefficient) {
            return coefficient >= 0;
        }
    );
}

bool direct_facet_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints,
    const std::vector<bool>& forced_zero = {}
) {
    for (std::size_t first = 0; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < constraints.size();
             ++second) {
            for (std::size_t third = second + 1U;
                 third < constraints.size();
                 ++third) {
                const std::array<Polynomial, 3> facets{
                    constraints[first],
                    constraints[second],
                    constraints[third]
                };
                std::array<Polynomial, 3> inverse;
                try {
                    inverse = inverse_facet_map(facets);
                } catch (const std::runtime_error&) {
                    continue;
                }
                Polynomial pulled =
                    substitute(chamber.margin, inverse);
                if (!forced_zero.empty()) {
                    const std::array<std::size_t, 3> selected{
                        first,
                        second,
                        third
                    };
                    const std::array<Polynomial, 3> coordinates{
                        forced_zero[selected[0]]
                            ? constant(0)
                            : Polynomial::variable(0),
                        forced_zero[selected[1]]
                            ? constant(0)
                            : Polynomial::variable(1),
                        forced_zero[selected[2]]
                            ? constant(0)
                            : Polynomial::variable(2)
                    };
                    pulled = substitute(pulled, coordinates);
                }
                const bool monomial = nonnegative_monomial(pulled);
                const bool newton =
                    !monomial && nonnegative_newton(pulled);
                if (monomial || newton) {
                    std::cout
                        << "SU2_K4_INTERMEDIATE_DIRECT"
                        << " mask=" << chamber.mask
                        << " facets=(" << first
                        << ',' << second
                        << ',' << third << ')'
                        << " basis="
                        << (monomial ? "monomial" : "newton")
                        << " result=PASS_EXACT_IDENTITY"
                        << std::endl;
                    return true;
                }
            }
        }
    }
    return false;
}

bool nonnegative_basis(const Polynomial& polynomial) {
    return nonnegative_monomial(polynomial)
        || nonnegative_newton(polynomial);
}

bool direct_branch_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
);

bool constraints_imply_sum(
    const std::vector<Polynomial>& constraints,
    const Polynomial& first,
    const Polynomial& second,
    int minimum_sum
) {
    z3::context context;
    z3::solver solver(context);
    const std::array<z3::expr, 3> variables{
        context.int_const("Q"),
        context.int_const("H"),
        context.int_const("Y")
    };
    for (const Polynomial& constraint : constraints) {
        solver.add(z3_affine(context, constraint, variables) >= 0);
    }
    solver.add(
        z3_affine(context, first + second, variables) < minimum_sum
    );
    return solver.check() == z3::unsat;
}

bool sum_cone_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    for (std::size_t first = 0; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < constraints.size();
             ++second) {
            int minimum_sum = 0;
            for (int candidate = 1; candidate <= 20; ++candidate) {
                if (
                    constraints_imply_sum(
                        constraints,
                        constraints[first],
                        constraints[second],
                        candidate
                    )
                ) {
                    minimum_sum = candidate;
                } else {
                    break;
                }
            }
            if (minimum_sum == 0) {
                continue;
            }
            for (std::size_t third = 0;
                 third < constraints.size();
                 ++third) {
                if (third == first || third == second) {
                    continue;
                }
                std::array<Polynomial, 3> inverse;
                try {
                    inverse = inverse_facet_map(
                        std::array<Polynomial, 3>{
                            constraints[first],
                            constraints[second],
                            constraints[third]
                        }
                    );
                } catch (const std::runtime_error&) {
                    continue;
                }
                const Polynomial pulled =
                    substitute(chamber.margin, inverse);
                const Polynomial x = Polynomial::variable(0);
                const Polynomial y = Polynomial::variable(1);
                const Polynomial z = Polynomial::variable(2);
                if (
                    !nonnegative_basis(
                        substitute(
                            pulled,
                            std::array<Polynomial, 3>{
                                x + constant(minimum_sum),
                                y,
                                z
                            }
                        )
                    )
                ) {
                    continue;
                }
                bool passed = true;
                for (int value = 0;
                     value < minimum_sum && passed;
                     ++value) {
                    passed = nonnegative_basis(
                        substitute(
                            pulled,
                            std::array<Polynomial, 3>{
                                constant(value),
                                y + constant(minimum_sum - value),
                                z
                            }
                        )
                    );
                }
                if (!passed) {
                    continue;
                }
                std::cout
                    << "SU2_K4_INTERMEDIATE_SUM_CONE"
                    << " mask=" << chamber.mask
                    << " facets=(" << first
                    << ',' << second
                    << ',' << third << ')'
                    << " minimum_sum=" << minimum_sum
                    << " slices=" << minimum_sum + 1
                    << " result=PASS_EXACT_PARTITION"
                    << std::endl;
                return true;
            }
        }
    }
    return false;
}

std::vector<std::vector<Polynomial>> sum_cone_branches(
    const std::vector<Polynomial>& constraints,
    std::size_t first,
    std::size_t second,
    int minimum_sum
) {
    std::vector<std::vector<Polynomial>> branches;
    std::vector<Polynomial> tail = constraints;
    tail.push_back(constraints[first] - constant(minimum_sum));
    branches.push_back(std::move(tail));
    for (int value = 0; value < minimum_sum; ++value) {
        std::vector<Polynomial> slice = constraints;
        slice.push_back(constraints[first] - constant(value));
        slice.push_back(constant(value) - constraints[first]);
        slice.push_back(
            constraints[second]
                - constant(minimum_sum - value)
        );
        branches.push_back(std::move(slice));
    }
    return branches;
}

bool double_sum_cone_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    struct SumPair {
        std::size_t first;
        std::size_t second;
        int minimum;
    };
    std::vector<SumPair> pairs;
    for (std::size_t first = 0; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < constraints.size();
             ++second) {
            int minimum = 0;
            for (int candidate = 1; candidate <= 10; ++candidate) {
                if (
                    constraints_imply_sum(
                        constraints,
                        constraints[first],
                        constraints[second],
                        candidate
                    )
                ) {
                    minimum = candidate;
                } else {
                    break;
                }
            }
            if (minimum > 0) {
                pairs.push_back(SumPair{first, second, minimum});
            }
        }
    }
    for (std::size_t left = 0; left < pairs.size(); ++left) {
        for (std::size_t right = left + 1U;
             right < pairs.size();
             ++right) {
            const SumPair& first_pair = pairs[left];
            const SumPair& second_pair = pairs[right];
            if (
                first_pair.first == second_pair.first
                || first_pair.first == second_pair.second
                || first_pair.second == second_pair.first
                || first_pair.second == second_pair.second
            ) {
                continue;
            }
            const auto first_branches = sum_cone_branches(
                constraints,
                first_pair.first,
                first_pair.second,
                first_pair.minimum
            );
            bool passed = true;
            std::size_t pieces = 0U;
            for (const auto& first_branch : first_branches) {
                const auto second_branches = sum_cone_branches(
                    first_branch,
                    second_pair.first,
                    second_pair.second,
                    second_pair.minimum
                );
                for (const auto& branch : second_branches) {
                    ++pieces;
                    if (!direct_branch_certificate(chamber, branch)) {
                        passed = false;
                        break;
                    }
                }
                if (!passed) {
                    break;
                }
            }
            if (!passed) {
                continue;
            }
            std::cout
                << "SU2_K4_INTERMEDIATE_DOUBLE_SUM"
                << " mask=" << chamber.mask
                << " pairs=("
                << first_pair.first << ',' << first_pair.second
                << ';'
                << second_pair.first << ',' << second_pair.second
                << ')'
                << " minima=("
                << first_pair.minimum << ','
                << second_pair.minimum << ')'
                << " pieces=" << pieces
                << " result=PASS_EXACT_PRODUCT_PARTITION"
                << std::endl;
            return true;
        }
    }
    return false;
}

std::vector<bool> forced_zero_constraints(
    const std::vector<Polynomial>& constraints
) {
    std::vector<bool> result(constraints.size(), false);
    for (std::size_t target = 0; target < constraints.size(); ++target) {
        z3::context context;
        z3::solver solver(context);
        const std::array<z3::expr, 3> variables{
            context.int_const("Q"),
            context.int_const("H"),
            context.int_const("Y")
        };
        for (const Polynomial& constraint : constraints) {
            solver.add(z3_affine(context, constraint, variables) >= 0);
        }
        solver.add(
            z3_affine(context, constraints[target], variables) > 0
        );
        result[target] = solver.check() == z3::unsat;
    }
    return result;
}

bool integer_feasible(
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
        solver.add(z3_affine(context, constraint, variables) >= 0);
    }
    return solver.check() == z3::sat;
}

bool direct_branch_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    // Most partition leaves already expose a useful integral facet
    // coordinate system.  Try that exact identity before asking Z3 to
    // remove redundant inequalities and detect forced zero facets.
    // The slower path below is still required for empty and
    // lower-dimensional leaves.
    if (direct_facet_certificate(chamber, constraints)) {
        return true;
    }
    if (!integer_feasible(constraints)) {
        return true;
    }
    const Chamber branch{
        chamber.mask,
        constraints,
        chamber.margin
    };
    const std::vector<Polynomial> reduced =
        irredundant_constraints(branch);
    return direct_facet_certificate(
        chamber,
        reduced,
        forced_zero_constraints(reduced)
    );
}

bool constant_three_sum_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    for (std::size_t first = 0; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < constraints.size();
             ++second) {
            for (std::size_t third = second + 1U;
                 third < constraints.size();
                 ++third) {
                const Polynomial sum =
                    constraints[first]
                    + constraints[second]
                    + constraints[third];
                int total = -1;
                if (sum.terms().empty()) {
                    total = 0;
                } else if (
                    sum.terms().size() == 1U
                    && sum.terms().begin()->first
                        == Exponent{0, 0, 0}
                    && sum.terms().begin()->second.denominator() == 1
                ) {
                    const Integer numerator =
                        sum.terms().begin()->second.numerator();
                    if (numerator >= 0 && numerator <= 10) {
                        total = numerator.convert_to<int>();
                    }
                }
                if (total < 0) {
                    continue;
                }

                bool passed = true;
                std::size_t pieces = 0U;
                for (int first_value = 0;
                     first_value <= total && passed;
                     ++first_value) {
                    for (int second_value = 0;
                         first_value + second_value <= total;
                         ++second_value) {
                        ++pieces;
                        std::vector<Polynomial> branch = constraints;
                        branch.push_back(
                            constraints[first]
                                - constant(first_value)
                        );
                        branch.push_back(
                            constant(first_value)
                                - constraints[first]
                        );
                        branch.push_back(
                            constraints[second]
                                - constant(second_value)
                        );
                        branch.push_back(
                            constant(second_value)
                                - constraints[second]
                        );
                        if (
                            !direct_branch_certificate(
                                chamber,
                                branch
                            )
                        ) {
                            passed = false;
                            break;
                        }
                    }
                }
                if (!passed) {
                    continue;
                }
                std::cout
                    << "SU2_K4_INTERMEDIATE_CONSTANT_SUM"
                    << " mask=" << chamber.mask
                    << " facets=(" << first
                    << ',' << second
                    << ',' << third << ')'
                    << " total=" << total
                    << " pieces=" << pieces
                    << " result=PASS_EXACT_COMPOSITIONS"
                    << std::endl;
                return true;
            }
        }
    }
    return false;
}

bool translated_orthant_certificate(
    const Polynomial& margin,
    const std::array<Polynomial, 3>& facets,
    int minimum_sum,
    int minimum_first_two = 0
) {
    if (
        minimum_sum < 0
        || minimum_first_two < 0
        || minimum_first_two > minimum_sum
    ) {
        return false;
    }
    std::array<Polynomial, 3> inverse;
    try {
        inverse = inverse_facet_map(facets);
    } catch (const std::runtime_error&) {
        return false;
    }
    const Polynomial pulled = substitute(margin, inverse);
    const Polynomial x = Polynomial::variable(0);
    const Polynomial y = Polynomial::variable(1);
    const Polynomial z = Polynomial::variable(2);
    for (int first = 0; first <= minimum_sum; ++first) {
        for (int second = 0;
             first + second <= minimum_sum;
             ++second) {
            if (first + second < minimum_first_two) {
                continue;
            }
            const int third = minimum_sum - first - second;
            if (
                !nonnegative_basis(
                    substitute(
                        pulled,
                        std::array<Polynomial, 3>{
                            x + constant(first),
                            y + constant(second),
                            z + constant(third)
                        }
                    )
                )
            ) {
                return false;
            }
        }
    }
    return true;
}

bool equal_sum_square_cone_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    struct Pair {
        std::size_t first;
        std::size_t second;
        int minimum;
    };
    std::vector<Pair> pairs;
    for (std::size_t first = 0; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < constraints.size();
             ++second) {
            int minimum = 0;
            for (int candidate = 1; candidate <= 20; ++candidate) {
                if (
                    constraints_imply_sum(
                        constraints,
                        constraints[first],
                        constraints[second],
                        candidate
                    )
                ) {
                    minimum = candidate;
                } else {
                    break;
                }
            }
            if (minimum > 0) {
                pairs.push_back(Pair{first, second, minimum});
            }
        }
    }
    for (std::size_t left = 0; left < pairs.size(); ++left) {
        for (std::size_t right = left + 1U;
             right < pairs.size();
             ++right) {
            const Pair& left_pair = pairs[left];
            const Pair& right_pair = pairs[right];
            if (
                left_pair.first == right_pair.first
                || left_pair.first == right_pair.second
                || left_pair.second == right_pair.first
                || left_pair.second == right_pair.second
            ) {
                continue;
            }
            const Polynomial offset_polynomial =
                constraints[right_pair.first]
                + constraints[right_pair.second]
                - constraints[left_pair.first]
                - constraints[left_pair.second];
            int offset = 0;
            if (offset_polynomial.terms().empty()) {
                offset = 0;
            } else if (
                offset_polynomial.terms().size() == 1U
                && offset_polynomial.terms().begin()->first
                    == Exponent{0, 0, 0}
                && offset_polynomial.terms().begin()
                    ->second.denominator() == 1
            ) {
                const Integer numerator =
                    offset_polynomial.terms().begin()
                        ->second.numerator();
                if (numerator < -1 || numerator > 1) {
                    continue;
                }
                offset = numerator.convert_to<int>();
            } else {
                continue;
            }
            if (
                right_pair.minimum
                    != left_pair.minimum + offset
            ) {
                continue;
            }
            const std::array<std::size_t, 2> left_indices{
                left_pair.first,
                left_pair.second
            };
            const std::array<std::size_t, 2> right_indices{
                right_pair.first,
                right_pair.second
            };
            for (std::size_t left_choice = 0;
                 left_choice < 2U;
                 ++left_choice) {
                for (std::size_t right_choice = 0;
                     right_choice < 2U;
                     ++right_choice) {
                    const Polynomial& a =
                        constraints[left_indices[left_choice]];
                    const Polynomial& b =
                        constraints[left_indices[1U - left_choice]];
                    const Polynomial& c =
                        constraints[right_indices[right_choice]];
                    const Polynomial& d =
                        constraints[right_indices[1U - right_choice]];

                    // On a>=c use coordinates
                    // (b,a-c,c), whose sum is the common total.
                    const bool first_branch =
                        translated_orthant_certificate(
                            chamber.margin,
                            std::array<Polynomial, 3>{
                                b,
                                a - c,
                                c
                            },
                            left_pair.minimum,
                            offset == -1 ? 1 : 0
                        );
                    if (!first_branch) {
                        continue;
                    }

                    // The complementary integer branch is c>=a+1.
                    // Coordinates (d,c-a-1,a) have sum total-1.
                    const bool second_branch =
                        translated_orthant_certificate(
                            chamber.margin,
                            std::array<Polynomial, 3>{
                                d,
                                c - a - constant(1),
                                a
                            },
                            left_pair.minimum + offset - 1
                        );
                    if (!second_branch) {
                        continue;
                    }
                    std::cout
                        << "SU2_K4_INTERMEDIATE_EQUAL_SUM_CONE"
                        << " mask=" << chamber.mask
                        << " pairs=("
                        << left_pair.first << ','
                        << left_pair.second << ';'
                        << right_pair.first << ','
                        << right_pair.second << ')'
                        << " orientations=("
                        << left_choice << ','
                        << right_choice << ')'
                        << " minimum=" << left_pair.minimum
                        << " offset=" << offset
                        << " result=PASS_EXACT_SQUARE_TRIANGULATION"
                        << std::endl;
                    return true;
                }
            }
        }
    }
    return false;
}

bool pair_cut_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    for (std::size_t first = 0; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < constraints.size();
             ++second) {
            std::vector<Polynomial> first_branch = constraints;
            first_branch.push_back(
                constraints[first] - constraints[second]
            );
            if (!direct_branch_certificate(chamber, first_branch)) {
                continue;
            }
            std::vector<Polynomial> second_branch = constraints;
            second_branch.push_back(
                constraints[second]
                - constraints[first]
                - constant(1)
            );
            if (!direct_branch_certificate(chamber, second_branch)) {
                continue;
            }
            std::cout
                << "SU2_K4_INTERMEDIATE_PAIR_CUT"
                << " mask=" << chamber.mask
                << " pair=(" << first
                << ',' << second << ')'
                << " result=PASS_EXACT_INTEGER_DICHOTOMY"
                << std::endl;
            return true;
        }
    }
    return false;
}

bool pair_cut_tree_constraints(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints,
    int depth
) {
    if (direct_branch_certificate(chamber, constraints)) {
        return true;
    }
    if (depth == 0) {
        return false;
    }
    for (std::size_t first = 0; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < constraints.size();
             ++second) {
            std::vector<Polynomial> first_branch = constraints;
            first_branch.push_back(
                constraints[first] - constraints[second]
            );
            if (
                !pair_cut_tree_constraints(
                    chamber,
                    first_branch,
                    depth - 1
                )
            ) {
                continue;
            }
            std::vector<Polynomial> second_branch = constraints;
            second_branch.push_back(
                constraints[second]
                - constraints[first]
                - constant(1)
            );
            if (
                !pair_cut_tree_constraints(
                    chamber,
                    second_branch,
                    depth - 1
                )
            ) {
                continue;
            }
            return true;
        }
    }
    return false;
}

bool pair_cut_tree_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    if (!pair_cut_tree_constraints(chamber, constraints, 2)) {
        return false;
    }
    std::cout
        << "SU2_K4_INTERMEDIATE_PAIR_TREE"
        << " mask=" << chamber.mask
        << " depth=2"
        << " result=PASS_EXACT_INTEGER_TREE"
        << std::endl;
    return true;
}

z3::expr z3_rational(
    z3::context& context,
    const Rational& value
) {
    std::ostringstream stream;
    stream << value.numerator();
    if (value.denominator() != 1) {
        stream << '/' << value.denominator();
    }
    return context.real_val(stream.str().c_str());
}

Rational model_rational(const z3::expr& value) {
    if (!value.is_numeral()) {
        throw std::runtime_error("Z3 model value is not rational");
    }
    std::string numerator;
    std::string denominator;
    if (
        !value.numerator().is_numeral(numerator)
        || !value.denominator().is_numeral(denominator)
    ) {
        throw std::runtime_error("Z3 rational parts are not numerals");
    }
    return Rational(Integer(numerator), Integer(denominator));
}

std::optional<Rational> affine_upper_bound(
    const std::vector<Polynomial>& constraints,
    std::size_t variable
) {
    z3::context context;
    z3::solver solver(context);
    std::vector<z3::expr> multipliers;
    multipliers.reserve(constraints.size());
    for (std::size_t index = 0; index < constraints.size(); ++index) {
        const std::string name =
            "bound_" + std::to_string(variable)
            + '_' + std::to_string(index);
        multipliers.push_back(context.real_const(name.c_str()));
        solver.add(multipliers.back() >= 0);
    }
    for (std::size_t coordinate = 0; coordinate < 3U; ++coordinate) {
        z3::expr sum = context.real_val(0);
        for (std::size_t index = 0; index < constraints.size(); ++index) {
            const Affine affine =
                affine_coefficients(constraints[index]);
            sum = sum
                + z3_rational(context, affine[coordinate + 1U])
                    * multipliers[index];
        }
        solver.add(
            sum == context.real_val(coordinate == variable ? -1 : 0)
        );
    }
    if (solver.check() != z3::sat) {
        return std::nullopt;
    }
    const z3::model model = solver.get_model();
    Polynomial reconstructed;
    for (std::size_t index = 0; index < constraints.size(); ++index) {
        const Rational multiplier = model_rational(
            model.eval(multipliers[index], true)
        );
        if (multiplier < 0) {
            throw std::runtime_error("negative upper-bound multiplier");
        }
        reconstructed +=
            constraints[index] * Polynomial(multiplier);
    }
    const Affine affine = affine_coefficients(reconstructed);
    for (std::size_t coordinate = 0; coordinate < 3U; ++coordinate) {
        const Rational expected =
            coordinate == variable ? Rational(-1) : Rational(0);
        if (affine[coordinate + 1U] != expected) {
            throw std::runtime_error(
                "upper-bound model fails exact slope replay"
            );
        }
    }
    return affine[0];
}

bool bounded_integer_certificate(
    const Chamber& chamber,
    const Formula& formula,
    const std::vector<Polynomial>& constraints
) {
    const std::optional<Rational> q_bound =
        affine_upper_bound(constraints, 0U);
    const std::optional<Rational> h_bound =
        affine_upper_bound(constraints, 1U);
    if (
        !q_bound.has_value()
        || !h_bound.has_value()
        || *q_bound < 1
        || *h_bound < 0
    ) {
        return false;
    }
    const Integer maximum_q_integer =
        q_bound->numerator() / q_bound->denominator();
    const Integer maximum_h_integer =
        h_bound->numerator() / h_bound->denominator();
    if (
        maximum_q_integer > 1000
        || maximum_h_integer > 10000
    ) {
        return false;
    }
    const int maximum_q = maximum_q_integer.convert_to<int>();
    const int maximum_h = maximum_h_integer.convert_to<int>();
    std::uint64_t points = 0U;
    for (int q_value = 1; q_value <= maximum_q; ++q_value) {
        for (int h_value = 0; h_value <= maximum_h; ++h_value) {
            for (int y_value = 0;
                 y_value <= 2 * q_value + 1 + h_value;
                 ++y_value) {
                const std::array<int, 3> values{
                    q_value,
                    h_value,
                    y_value
                };
                std::uint64_t mask = 0U;
                for (std::size_t index = 0;
                     index < formula.hinges.size();
                     ++index) {
                    if (
                        evaluate(formula.hinges[index].slack, values)
                        >= 0
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
                        "bounded chamber contains a negative point"
                    );
                }
            }
        }
    }
    if (points == 0U) {
        throw std::runtime_error(
            "bounded feasible chamber has no enumerated point"
        );
    }
    std::cout
        << "SU2_K4_INTERMEDIATE_BOUNDED"
        << " mask=" << chamber.mask
        << " Q_bound=" << *q_bound
        << " H_bound=" << *h_bound
        << " points=" << points
        << " result=PASS_EXACT_FINITE"
        << std::endl;
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::optional<std::uint64_t> selected_mask;
        bool describe = false;
        if (argc == 3 && std::string(argv[1]) == "--mask") {
            selected_mask = std::stoull(argv[2]);
        } else if (argc == 3 && std::string(argv[1]) == "--describe") {
            selected_mask = std::stoull(argv[2]);
            describe = true;
        } else if (argc != 1) {
            throw std::runtime_error(
                "usage: prove_su2_k4_intermediate "
                "[--mask MASK|--describe MASK]"
            );
        }
        const Formula formula = make_formula();
        const std::vector<std::uint64_t> masks =
            selected_mask.has_value()
                ? std::vector<std::uint64_t>{*selected_mask}
                : feasible_masks(formula);
        std::size_t certified = 0U;
        for (std::size_t position = 0; position < masks.size(); ++position) {
            const Chamber chamber =
                make_chamber(formula, masks[position]);
            if (
                selected_mask.has_value()
                && !integer_feasible(chamber.constraints)
            ) {
                throw std::runtime_error(
                    "selected activation mask is infeasible"
                );
            }
            const std::vector<Polynomial> domain_constraints(
                chamber.constraints.begin(),
                chamber.constraints.begin() + 4
            );
            bool passed = direct_facet_certificate(
                chamber,
                domain_constraints
            );
            if (!passed) {
                const std::vector<Polynomial> constraints =
                    irredundant_constraints(chamber);
                if (describe) {
                    for (std::size_t index = 0;
                         index < constraints.size();
                         ++index) {
                        const Affine affine =
                            affine_coefficients(constraints[index]);
                        std::cout
                            << "SU2_K4_INTERMEDIATE_FACET"
                            << " mask=" << chamber.mask
                            << " index=" << index
                            << " affine=("
                            << affine[0] << ','
                            << affine[1] << ','
                            << affine[2] << ','
                            << affine[3] << ")\n";
                    }
                    return EXIT_SUCCESS;
                }
                const std::vector<bool> forced_zero =
                    forced_zero_constraints(constraints);
                passed = direct_facet_certificate(
                    chamber,
                    constraints,
                    forced_zero
                );
                if (!passed) {
                    passed = constant_three_sum_certificate(
                        chamber,
                        constraints
                    );
                }
                if (!passed) {
                    passed = bounded_integer_certificate(
                        chamber,
                        formula,
                        constraints
                    );
                }
                if (!passed) {
                    passed = sum_cone_certificate(
                        chamber,
                        constraints
                    );
                }
                if (!passed) {
                    passed = equal_sum_square_cone_certificate(
                        chamber,
                        constraints
                    );
                }
                if (!passed) {
                    passed = double_sum_cone_certificate(
                        chamber,
                        constraints
                    );
                }
                if (!passed) {
                    passed = pair_cut_certificate(
                        chamber,
                        constraints
                    );
                }
                if (!passed) {
                    passed = pair_cut_tree_certificate(
                        chamber,
                        constraints
                    );
                }
            }
            if (passed) {
                ++certified;
            } else {
                std::cout
                    << "SU2_K4_INTERMEDIATE_UNRESOLVED"
                    << " mask=" << chamber.mask
                    << " result=NEEDS_STRONGER_CERTIFICATE"
                    << std::endl;
            }
            if ((position + 1U) % 25U == 0U) {
                std::cerr
                    << "SU2_K4_INTERMEDIATE"
                    << " progress=" << position + 1U
                    << '/' << masks.size()
                    << " certified=" << certified
                    << std::endl;
            }
        }
        const bool complete = selected_mask.has_value()
            ? certified == 1U
            : certified == masks.size();
        std::cout
            << "SU2_K4_INTERMEDIATE"
            << " hinges=" << formula.hinges.size()
            << " feasible_chambers=" << masks.size()
            << " certified_chambers=" << certified
            << " result="
            << (complete ? "PASS_EXACT_CERTIFICATE" : "INCOMPLETE")
            << '\n';
        return complete ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_K4_INTERMEDIATE FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
