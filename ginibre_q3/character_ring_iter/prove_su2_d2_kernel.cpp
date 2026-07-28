#include <array>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

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

Polynomial multiplicity_four(
    const Polynomial& q,
    const Polynomial& y,
    bool second_image
) {
    Polynomial result =
        binomial(scale(q, 4) - y + constant(2), 2);
    if (second_image) {
        result -= scale(
            binomial(scale(q, 2) - y + constant(1), 2),
            4
        );
    }
    return result;
}

Polynomial multiplicity_five(
    const Polynomial& q,
    const Polynomial& y,
    bool second_image,
    bool third_image
) {
    Polynomial result =
        binomial(scale(q, 5) - y + constant(3), 3);
    if (second_image) {
        result -= scale(
            binomial(scale(q, 3) - y + constant(2), 3),
            5
        );
    }
    if (third_image) {
        result += scale(
            binomial(q - y + constant(1), 3),
            10
        );
    }
    return result;
}

Polynomial stable_margin(
    const Polynomial& q,
    const Polynomial& y,
    bool four_second,
    bool five_second,
    bool five_third
) {
    const Polynomial f4 = scale(q, 2) + constant(1);
    const Polynomial f5 = scale(
        scale(q * q, 5) + scale(q, 5) + constant(2),
        1,
        2
    );
    return f4
            * multiplicity_five(
                q, y, five_second, five_third
            )
        - f5 * multiplicity_four(q, y, four_second);
}

Polynomial nearest_wall_margin(
    const Polynomial& q,
    const Polynomial& y,
    bool four_second,
    bool four_reflection,
    bool five_third,
    bool five_upper_reflection,
    bool five_lower_reflection
) {
    Polynomial fourth =
        multiplicity_four(q, y, four_second);
    if (four_reflection) {
        fourth -= binomial(y - constant(1), 2);
    }
    Polynomial fifth =
        multiplicity_five(q, y, true, five_third)
        - binomial(q + y, 3);
    if (five_upper_reflection) {
        fifth += scale(
            binomial(y - q - constant(1), 3),
            5
        );
    }
    if (five_lower_reflection) {
        fifth += binomial(q - y - constant(1), 3);
    }
    const Polynomial f4 = scale(q, 2) + constant(1);
    const Polynomial stable_f5 = scale(
        scale(q * q, 5) + scale(q, 5) + constant(2),
        1,
        2
    );
    const Polynomial finite_f5 =
        stable_f5 - binomial(q - constant(1), 2);
    return f4 * fifth - finite_f5 * fourth;
}

Integer integer_binomial(int top, int order) {
    if (top < order) {
        return 0;
    }
    Integer result = 1;
    for (int index = 0; index < order; ++index) {
        result *= top - index;
        result /= index + 1;
    }
    return result;
}

Integer nearest_wall_margin_integer(int q, int y) {
    const Integer fourth =
        integer_binomial(4 * q - y + 2, 2)
        - 4 * integer_binomial(2 * q - y + 1, 2)
        - integer_binomial(y - 1, 2);
    const Integer fifth =
        integer_binomial(5 * q - y + 3, 3)
        - 5 * integer_binomial(3 * q - y + 2, 3)
        + 10 * integer_binomial(q - y + 1, 3)
        - integer_binomial(q + y, 3)
        + 5 * integer_binomial(y - q - 1, 3)
        + integer_binomial(q - y - 1, 3);
    const Integer f4 = 2 * q + 1;
    const Integer f5 =
        (5 * Integer(q) * q + 5 * q + 2) / 2
        - integer_binomial(q - 1, 2);
    return f4 * fifth - f5 * fourth;
}

void certify(
    const std::string& name,
    const Polynomial& polynomial
) {
    std::size_t positive = 0U;
    std::size_t zero = 0U;
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        static_cast<void>(exponent);
        if (coefficient < 0) {
            throw std::runtime_error(
                name + " has a negative coefficient"
            );
        }
        if (coefficient == 0) {
            ++zero;
        } else {
            ++positive;
        }
    }
    std::cout
        << "SU2_D2_STABLE_CHAMBER"
        << " name=" << name
        << " positive_coefficients=" << positive
        << " zero_coefficients=" << zero
        << " result=PASS_COEFFICIENTS\n";
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        std::cout
            << "  coefficient=" << coefficient
            << " exponent=(" << exponent[0]
            << ',' << exponent[1]
            << ',' << exponent[2] << ")\n";
    }
}

}  // namespace

int main() {
    try {
        const Polynomial x = Polynomial::variable(0);
        const Polynomial y = Polynomial::variable(1);

        certify(
            "Y_le_Q_minus_2",
            stable_margin(
                x + y + constant(2),
                x,
                true,
                true,
                true
            )
        );
        certify(
            "Q_minus_1_le_Y_le_2Q_minus_1_left_edge",
            stable_margin(
                x + constant(1),
                x,
                true,
                true,
                false
            )
        );
        certify(
            "Q_minus_1_le_Y_le_2Q_minus_1_interior",
            stable_margin(
                x + y + constant(1),
                scale(x, 2) + y + constant(1),
                true,
                true,
                false
            )
        );
        certify(
            "2Q_le_Y_le_3Q_minus_1",
            stable_margin(
                x + y + constant(1),
                scale(x, 3) + scale(y, 2) + constant(2),
                false,
                true,
                false
            )
        );
        certify(
            "3Q_le_Y_le_4Q_left_edge",
            stable_margin(
                x + constant(1),
                scale(x, 3) + constant(3),
                false,
                false,
                false
            )
        );
        certify(
            "3Q_le_Y_le_4Q_interior",
            stable_margin(
                x + y + constant(1),
                scale(x, 4) + scale(y, 3) + constant(4),
                false,
                false,
                false
            )
        );
        Integer nearest_small_minimum = 0;
        bool nearest_small_initialized = false;
        int nearest_small_q = 0;
        int nearest_small_y = 0;
        for (int q_value = 1; q_value <= 6; ++q_value) {
            for (int y_value = 0;
                 y_value <= 2 * q_value + 1;
                 ++y_value) {
                const Integer margin =
                    nearest_wall_margin_integer(q_value, y_value);
                if (!nearest_small_initialized
                    || margin < nearest_small_minimum) {
                    nearest_small_initialized = true;
                    nearest_small_minimum = margin;
                    nearest_small_q = q_value;
                    nearest_small_y = y_value;
                }
                if (margin < 0) {
                    throw std::runtime_error(
                        "nearest-wall small case is negative"
                    );
                }
            }
        }
        std::cout
            << "SU2_D2_NEAREST_WALL_SMALL"
            << " maximum_Q=6"
            << " minimum=" << nearest_small_minimum
            << " witness=(" << nearest_small_q
            << ',' << nearest_small_y << ')'
            << " result=PASS_EXACT\n";

        const Polynomial large_q = x + constant(7);
        for (int y_value = 0; y_value <= 2; ++y_value) {
            certify(
                "nearest_Y_" + std::to_string(y_value),
                nearest_wall_margin(
                    large_q,
                    constant(y_value),
                    true,
                    false,
                    true,
                    false,
                    true
                )
            );
        }
        certify(
            "nearest_3_le_Y_le_Q_minus_4",
            nearest_wall_margin(
                x + y + constant(7),
                x + constant(3),
                true,
                true,
                true,
                false,
                true
            )
        );
        for (int offset = -3; offset <= -2; ++offset) {
            certify(
                "nearest_Y_eq_Q"
                    + std::to_string(offset),
                nearest_wall_margin(
                    large_q,
                    large_q + constant(offset),
                    true,
                    true,
                    true,
                    false,
                    false
                )
            );
        }
        for (int offset = -1; offset <= 3; ++offset) {
            certify(
                "nearest_Y_eq_Q_plus_"
                    + std::to_string(offset),
                nearest_wall_margin(
                    large_q,
                    large_q + constant(offset),
                    true,
                    true,
                    false,
                    false,
                    false
                )
            );
        }
        certify(
            "nearest_Q_plus_4_le_Y_le_2Q_minus_1",
            nearest_wall_margin(
                x + y + constant(5),
                scale(x, 2) + y + constant(9),
                true,
                true,
                false,
                true,
                false
            )
        );
        for (int offset = 0; offset <= 1; ++offset) {
            certify(
                "nearest_Y_eq_2Q_plus_"
                    + std::to_string(offset),
                nearest_wall_margin(
                    large_q,
                    scale(large_q, 2) + constant(offset),
                    false,
                    true,
                    false,
                    true,
                    false
                )
            );
        }
        std::cout
            << "SU2_D2_STABLE_KERNEL result=PASS_EXACT\n";
        std::cout
            << "SU2_D2_NEAREST_WALL_KERNEL result=PASS_EXACT\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_D2_STABLE_KERNEL FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
