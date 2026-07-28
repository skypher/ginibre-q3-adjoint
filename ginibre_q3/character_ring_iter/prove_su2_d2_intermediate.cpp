#include <array>
#include <cstdlib>
#include <iostream>
#include <map>
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
                for (std::size_t index = 0; index < exponent.size();
                     ++index) {
                    exponent[index] =
                        left_exponent[index] + right_exponent[index];
                }
                product[exponent] += left_coefficient * right_coefficient;
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

Polynomial scale(
    const Polynomial& value,
    long numerator,
    long denominator = 1
) {
    return value * Polynomial(Rational(numerator, denominator));
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
        for (std::size_t index = 0; index < exponent.size(); ++index) {
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
    std::vector<Polynomial> constraints;
    Polynomial margin;
};

Polynomial selected_binomial(
    const Hinge& hinge,
    unsigned mask,
    std::size_t index
) {
    if ((mask & (1U << index)) == 0U) {
        return Polynomial();
    }
    return binomial(hinge.top, hinge.order);
}

Chamber make_chamber(
    unsigned mask,
    const std::vector<Hinge>& hinges
) {
    const Polynomial q = Polynomial::variable(0);
    const Polynomial h = Polynomial::variable(1);
    const Polynomial y = Polynomial::variable(2);

    Chamber chamber;
    chamber.constraints = {
        q - constant(1),
        h,
        scale(q, 3) - constant(2) - h,
        y - constant(1),
        scale(q, 2) + constant(1) + h - y
    };
    for (std::size_t index = 0; index < hinges.size(); ++index) {
        if ((mask & (1U << index)) != 0U) {
            chamber.constraints.push_back(
                hinges[index].top - constant(hinges[index].order)
            );
        } else {
            chamber.constraints.push_back(
                constant(hinges[index].order - 1) - hinges[index].top
            );
        }
    }

    const Polynomial m4 =
        selected_binomial(hinges[0], mask, 0)
        - scale(selected_binomial(hinges[1], mask, 1), 4);
    const Polynomial m5 =
        binomial(scale(q, 5) - y + constant(3), 3)
        - scale(selected_binomial(hinges[2], mask, 2), 5)
        + scale(selected_binomial(hinges[3], mask, 3), 10);
    const Polynomial u4 =
        m4 - selected_binomial(hinges[4], mask, 4);
    const Polynomial u5 =
        m5
        - selected_binomial(hinges[5], mask, 5)
        + scale(selected_binomial(hinges[6], mask, 6), 5)
        + selected_binomial(hinges[7], mask, 7);
    const Polynomial stable_f5 = scale(
        scale(q * q, 5) + scale(q, 5) + constant(2),
        1,
        2
    );
    const Polynomial f5 =
        stable_f5 - selected_binomial(hinges[8], mask, 8);
    chamber.margin = (scale(q, 2) + constant(1)) * u5 - f5 * u4;
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

z3::expr polynomial_value(
    z3::context& context,
    const Polynomial& polynomial,
    const std::array<z3::expr, 3>& variables
) {
    z3::expr result = context.real_val(0);
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        z3::expr term = rational_value(context, coefficient);
        for (std::size_t index = 0; index < exponent.size(); ++index) {
            for (int power = 0; power < exponent[index]; ++power) {
                term = term * variables[index];
            }
        }
        result = result + term;
    }
    return result;
}

bool integer_feasible(
    z3::context& context,
    const std::vector<Polynomial>& constraints
) {
    const std::array<z3::expr, 3> variables{
        context.int_const("Q"),
        context.int_const("H"),
        context.int_const("Y")
    };
    z3::solver solver(context);
    for (const Polynomial& constraint : constraints) {
        solver.add(polynomial_value(context, constraint, variables) >= 0);
    }
    return solver.check() == z3::sat;
}

void certify_nonnegative_coefficients(
    const std::string& name,
    const Polynomial& polynomial
) {
    std::size_t positive = 0U;
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        static_cast<void>(exponent);
        if (coefficient < 0) {
            throw std::runtime_error(
                name + " has a negative coefficient"
            );
        }
        if (coefficient > 0) {
            ++positive;
        }
    }
    std::cout
        << "SU2_D2_DIRECT_CONE"
        << " name=" << name
        << " positive_coefficients=" << positive
        << " result=PASS_COEFFICIENTS"
        << std::endl;
}

void certify_chamber_53_cone(
    const std::string& name,
    const Chamber& chamber,
    const std::array<Polynomial, 3>& values
) {
    for (const Polynomial& constraint : chamber.constraints) {
        const Polynomial pulled_back = substitute(constraint, values);
        for (const auto& [exponent, coefficient]
             : pulled_back.terms()) {
            static_cast<void>(exponent);
            if (coefficient < 0) {
                throw std::runtime_error(
                    name + " leaves chamber 53"
                );
            }
        }
    }
    certify_nonnegative_coefficients(
        name,
        substitute(chamber.margin, values)
    );
}

void certify_chamber_53(const Chamber& chamber) {
    const Polynomial x = Polynomial::variable(0);
    const Polynomial v = Polynomial::variable(1);
    const Polynomial z = Polynomial::variable(2);

    certify_chamber_53_cone(
        "mask53_a_zero",
        chamber,
        std::array<Polynomial, 3>{
            scale(x, 2) + v + constant(2),
            x + v,
            scale(
                scale(x, 2) + v + constant(2),
                2
            )
        }
    );
    certify_chamber_53_cone(
        "mask53_b_zero",
        chamber,
        std::array<Polynomial, 3>{
            scale(x, 2) + z + constant(2),
            scale(x, 3) + z + constant(1),
            scale(x, 6) + scale(z, 2) + constant(5)
        }
    );
    certify_chamber_53_cone(
        "mask53_b_above_a",
        chamber,
        std::array<Polynomial, 3>{
            x + scale(v, 2) + z + constant(2),
            x + v + z,
            scale(x, 3) + scale(v, 4)
                + scale(z, 2) + constant(5)
        }
    );
    certify_chamber_53_cone(
        "mask53_a_above_b",
        chamber,
        std::array<Polynomial, 3>{
            x + scale(v, 2) + z + constant(3),
            x + scale(v, 3) + z + constant(2),
            scale(x, 3) + scale(v, 6)
                + scale(z, 2) + constant(8)
        }
    );
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

bool handelman_feasible(
    z3::context& context,
    const Chamber& chamber,
    unsigned chamber_index,
    int maximum_degree
) {
    if (chamber_index == 53U) {
        certify_chamber_53(chamber);
        std::cout
            << "SU2_D2_HANDELMAN"
            << " chamber=53"
            << " products=0"
            << " nonzero_coefficients=0"
            << " result=PASS_EXACT_DIRECT_CONES"
            << std::endl;
        return true;
    }
    std::map<Polynomial, Polynomial> unique_products;
    enumerate_products_recursive(
        chamber.constraints,
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

    z3::solver solver(context);
    std::vector<z3::expr> coefficients;
    coefficients.reserve(products.size());
    for (std::size_t index = 0; index < products.size(); ++index) {
        const std::string name =
            "c_" + std::to_string(chamber_index)
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
                for (std::size_t index = 0;
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
        for (std::size_t index = 0; index < products.size(); ++index) {
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
        << "SU2_D2_HANDELMAN"
        << " chamber=" << chamber_index
        << " products=" << products.size()
        << " nonzero_coefficients=" << nonzero_coefficients
        << " result="
        << (result == z3::sat
            ? "PASS_EXACT_IDENTITY"
            : result == z3::unsat ? "UNSAT" : "UNKNOWN")
        << std::endl;
    return result == z3::sat;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Polynomial q = Polynomial::variable(0);
        const Polynomial h = Polynomial::variable(1);
        const Polynomial y = Polynomial::variable(2);
        const std::vector<Hinge> hinges{
            {"m4_first", scale(q, 4) - y + constant(2), 2},
            {"m4_second", scale(q, 2) - y + constant(1), 2},
            {"m5_second", scale(q, 3) - y + constant(2), 3},
            {"m5_third", q - y + constant(1), 3},
            {"u4_reflection", y - scale(h, 2) - constant(1), 2},
            {"u5_lower_first", q + y - scale(h, 2), 3},
            {
                "u5_lower_second",
                y - q - scale(h, 2) - constant(1),
                3
            },
            {
                "u5_upper",
                q - scale(h, 2) - y - constant(1),
                3
            },
            {"f5_reflection", q - scale(h, 2) - constant(1), 2}
        };
        z3::context context;
        bool patterns_only = false;
        int selected_chamber = -1;
        int maximum_degree = 4;
        if (argc == 2 && std::string(argv[1]) == "--patterns") {
            patterns_only = true;
        } else if (argc == 3
                   && std::string(argv[1]) == "--chamber") {
            selected_chamber = std::stoi(argv[2]);
            if (selected_chamber < 0
                || selected_chamber >= (1 << hinges.size())) {
                throw std::runtime_error("chamber mask is out of range");
            }
        } else if (argc == 5
                   && std::string(argv[1]) == "--chamber"
                   && std::string(argv[3]) == "--degree") {
            selected_chamber = std::stoi(argv[2]);
            maximum_degree = std::stoi(argv[4]);
            if (selected_chamber < 0
                || selected_chamber >= (1 << hinges.size())
                || maximum_degree < 4
                || maximum_degree > 8) {
                throw std::runtime_error(
                    "invalid chamber mask or Handelman degree"
                );
            }
        } else if (argc != 1) {
            throw std::runtime_error(
                "usage: prove_su2_d2_intermediate [--patterns] "
                "[--chamber MASK [--degree DEGREE]]"
            );
        }
        std::size_t feasible = 0U;
        std::size_t certified = 0U;
        for (unsigned mask = 0U; mask < (1U << hinges.size()); ++mask) {
            if (selected_chamber >= 0
                && mask != static_cast<unsigned>(selected_chamber)) {
                continue;
            }
            const Chamber chamber = make_chamber(mask, hinges);
            if (!integer_feasible(context, chamber.constraints)) {
                continue;
            }
            ++feasible;
            std::cout
                << "SU2_D2_INTERMEDIATE_CHAMBER"
                << " mask=" << mask
                << " result=FEASIBLE"
                << std::endl;
            if (patterns_only) {
                continue;
            }
            if (handelman_feasible(
                    context,
                    chamber,
                    mask,
                    maximum_degree
                )) {
                ++certified;
            }
        }
        if (patterns_only) {
            std::cout
                << "SU2_D2_INTERMEDIATE_PATTERNS"
                << " feasible_chambers=" << feasible
                << " result=PASS_ENUMERATION\n";
            return EXIT_SUCCESS;
        }
        std::cout
            << "SU2_D2_INTERMEDIATE"
            << " feasible_chambers=" << feasible
            << " certified_chambers=" << certified
            << " result="
            << (feasible == certified
                ? "PASS_HANDELMAN"
                : "INCOMPLETE")
            << '\n';
        return feasible == certified ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_D2_INTERMEDIATE FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
