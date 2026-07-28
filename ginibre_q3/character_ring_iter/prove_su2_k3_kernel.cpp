#include <algorithm>
#include <array>
#include <cstdint>
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
using Exponent = std::array<int, 2>;

class Polynomial {
public:
    Polynomial() = default;

    explicit Polynomial(const Rational& value) {
        if (value != 0) {
            terms_[Exponent{0, 0}] = value;
        }
    }

    static Polynomial variable(int index) {
        Polynomial result;
        Exponent exponent{0, 0};
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
                const Exponent exponent{
                    left_exponent[0] + right_exponent[0],
                    left_exponent[1] + right_exponent[1]
                };
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
    const std::array<Polynomial, 2>& values
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

Polynomial multiplicity_six(
    const Polynomial& q,
    const Polynomial& y,
    bool second,
    bool third
) {
    Polynomial result =
        binomial(scale(q, 6) - y + constant(4), 4);
    if (second) {
        result -= scale(
            binomial(scale(q, 4) - y + constant(3), 4),
            6
        );
    }
    if (third) {
        result += scale(
            binomial(scale(q, 2) - y + constant(2), 4),
            15
        );
    }
    return result;
}

Polynomial multiplicity_seven(
    const Polynomial& q,
    const Polynomial& y,
    bool second,
    bool third,
    bool fourth
) {
    Polynomial result =
        binomial(scale(q, 7) - y + constant(5), 5);
    if (second) {
        result -= scale(
            binomial(scale(q, 5) - y + constant(4), 5),
            7
        );
    }
    if (third) {
        result += scale(
            binomial(scale(q, 3) - y + constant(3), 5),
            21
        );
    }
    if (fourth) {
        result -= scale(
            binomial(q - y + constant(2), 5),
            35
        );
    }
    return result;
}

Polynomial stable_margin(
    const Polynomial& q,
    const Polynomial& y,
    bool six_first,
    bool six_second,
    bool six_third,
    bool seven_first,
    bool seven_second,
    bool seven_third,
    bool seven_fourth
) {
    const Polynomial f6 =
        multiplicity_six(q, constant(0), true, true);
    const Polynomial f7 =
        multiplicity_seven(q, constant(0), true, true, true);
    const Polynomial u6 = six_first
        ? multiplicity_six(q, y, six_second, six_third)
        : Polynomial();
    const Polynomial u7 = seven_first
        ? multiplicity_seven(
            q,
            y,
            seven_second,
            seven_third,
            seven_fourth
        )
        : Polynomial();
    return f6 * u7 - f7 * u6;
}

Polynomial nearest_six(
    const Polynomial& q,
    const Polynomial& y,
    bool ordinary_third,
    bool upper_y,
    bool lower_image
) {
    Polynomial result =
        multiplicity_six(q, y, true, ordinary_third)
        - binomial(scale(q, 2) + y + constant(1), 4);
    if (upper_y) {
        result += scale(binomial(y, 4), 6);
    }
    if (lower_image) {
        result += binomial(scale(q, 2) - y, 4);
    }
    return result;
}

Polynomial nearest_seven(
    const Polynomial& q,
    const Polynomial& y,
    bool ordinary_fourth,
    bool upper_image,
    bool lower_image,
    bool far_upper_image
) {
    Polynomial result =
        multiplicity_seven(q, y, true, true, ordinary_fourth)
        - binomial(scale(q, 3) + y + constant(2), 5)
        + scale(binomial(q + y + constant(1), 5), 7)
        + binomial(scale(q, 3) - y + constant(1), 5);
    if (upper_image) {
        result -= scale(binomial(y - q, 5), 21);
    }
    if (lower_image) {
        result -= scale(binomial(q - y, 5), 7);
    }
    if (far_upper_image) {
        result -= binomial(y - q - constant(2), 5);
    }
    return result;
}

Polynomial nearest_margin(
    const Polynomial& q,
    const Polynomial& y,
    bool six_ordinary_third,
    bool six_upper_y,
    bool six_lower_image,
    bool seven_ordinary_fourth,
    bool seven_upper_image,
    bool seven_lower_image,
    bool seven_far_upper_image
) {
    const Polynomial f6 = nearest_six(
        q,
        constant(0),
        true,
        false,
        true
    );
    const Polynomial f7 = nearest_seven(
        q,
        constant(0),
        true,
        false,
        true,
        false
    );
    return f6
            * nearest_seven(
                q,
                y,
                seven_ordinary_fourth,
                seven_upper_image,
                seven_lower_image,
                seven_far_upper_image
            )
        - f7
            * nearest_six(
                q,
                y,
                six_ordinary_third,
                six_upper_y,
                six_lower_image
            );
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

Integer multiplicity_six_integer(int q, int y) {
    return integer_binomial(6 * q - y + 4, 4)
        - 6 * integer_binomial(4 * q - y + 3, 4)
        + 15 * integer_binomial(2 * q - y + 2, 4);
}

Integer multiplicity_seven_integer(int q, int y) {
    return integer_binomial(7 * q - y + 5, 5)
        - 7 * integer_binomial(5 * q - y + 4, 5)
        + 21 * integer_binomial(3 * q - y + 3, 5)
        - 35 * integer_binomial(q - y + 2, 5);
}

Integer nearest_six_integer(int q, int y) {
    return multiplicity_six_integer(q, y)
        - integer_binomial(2 * q + y + 1, 4)
        + 6 * integer_binomial(y, 4)
        + integer_binomial(2 * q - y, 4);
}

Integer nearest_seven_integer(int q, int y) {
    return multiplicity_seven_integer(q, y)
        - integer_binomial(3 * q + y + 2, 5)
        + 7 * integer_binomial(q + y + 1, 5)
        - 21 * integer_binomial(y - q, 5)
        + integer_binomial(3 * q - y + 1, 5)
        - 7 * integer_binomial(q - y, 5)
        - integer_binomial(y - q - 2, 5);
}

bool fuses(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= std::min(
            source + label,
            2 * level - source - label
        )
        && ((source + label + target) & 1) == 0;
}

std::vector<Integer> multiply(
    int level,
    int label,
    const std::vector<Integer>& state
) {
    std::vector<Integer> next(static_cast<std::size_t>(level + 1));
    for (int source = 0; source <= level; ++source) {
        const Integer& coefficient =
            state[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        for (int target = 0; target <= level; ++target) {
            if (fuses(level, label, source, target)) {
                next[static_cast<std::size_t>(target)] += coefficient;
            }
        }
    }
    return next;
}

std::uint64_t audit_formula_case(
    int half_level,
    int half_label,
    bool nearest
) {
    const int level = 2 * half_level;
    const int label = 2 * half_label;
    std::vector<Integer> state(static_cast<std::size_t>(level + 1));
    state[0] = 1;
    std::uint64_t comparisons = 0U;
    for (int power_index = 1; power_index <= 7; ++power_index) {
        state = multiply(level, label, state);
        if (power_index != 6 && power_index != 7) {
            continue;
        }
        for (int y_value = 0;
             y_value <= half_level;
             ++y_value) {
            const Integer expected = power_index == 6
                ? (
                    nearest
                        ? nearest_six_integer(half_label, y_value)
                        : multiplicity_six_integer(
                            half_label,
                            y_value
                        )
                )
                : (
                    nearest
                        ? nearest_seven_integer(half_label, y_value)
                        : multiplicity_seven_integer(
                            half_label,
                            y_value
                        )
                );
            if (
                state[static_cast<std::size_t>(2 * y_value)]
                != expected
            ) {
                std::cout
                    << "SU2_K3_FORMULA_MISMATCH"
                    << " K=" << half_level
                    << " Q=" << half_label
                    << " nearest=" << (nearest ? 1 : 0)
                    << " power=" << power_index
                    << " Y=" << y_value
                    << " recurrence="
                    << state[static_cast<std::size_t>(2 * y_value)]
                    << " formula=" << expected << '\n';
                throw std::runtime_error(
                    "finite recurrence disagrees with formula"
                );
            }
            ++comparisons;
        }
    }
    return comparisons;
}

void certify(
    const std::string& name,
    const Polynomial& polynomial
) {
    std::size_t positive = 0U;
    std::size_t negative = 0U;
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        if (coefficient < 0) {
            ++negative;
            std::cout
                << "SU2_K3_NEGATIVE_COEFFICIENT"
                << " name=" << name
                << " exponent=(" << exponent[0]
                << ',' << exponent[1] << ')'
                << " coefficient=" << coefficient << '\n';
        }
        if (coefficient > 0) {
            ++positive;
        }
    }
    if (negative != 0U) {
        throw std::runtime_error(
            name + " has a negative coefficient"
        );
    }
    std::cout
        << "SU2_K3_CHAMBER"
        << " name=" << name
        << " positive_coefficients=" << positive
        << " result=PASS_COEFFICIENTS\n";
}

bool has_nonnegative_coefficients(const Polynomial& polynomial) {
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        static_cast<void>(exponent);
        if (coefficient < 0) {
            return false;
        }
    }
    return true;
}

void certify_sum_cone(
    const std::string& name,
    const Polynomial& polynomial,
    int minimum_sum
) {
    const Polynomial x = Polynomial::variable(0);
    const Polynomial y = Polynomial::variable(1);
    for (int x_threshold = 0; x_threshold <= 30; ++x_threshold) {
        for (int y_threshold = 0; y_threshold <= 30; ++y_threshold) {
            const Polynomial interior = substitute(
                polynomial,
                std::array<Polynomial, 2>{
                    x + constant(x_threshold),
                    y + constant(y_threshold)
                }
            );
            if (!has_nonnegative_coefficients(interior)) {
                continue;
            }
            bool valid = true;
            for (int x_value = 0;
                 x_value < x_threshold && valid;
                 ++x_value) {
                const int y_start =
                    std::max(0, minimum_sum - x_value);
                const Polynomial vertical = substitute(
                    polynomial,
                    std::array<Polynomial, 2>{
                        constant(x_value),
                        y + constant(y_start)
                    }
                );
                valid = has_nonnegative_coefficients(vertical);
            }
            for (int y_value = 0;
                 y_value < y_threshold && valid;
                 ++y_value) {
                const int x_start = std::max(
                    x_threshold,
                    minimum_sum - y_value
                );
                const Polynomial horizontal = substitute(
                    polynomial,
                    std::array<Polynomial, 2>{
                        x + constant(x_start),
                        constant(y_value)
                    }
                );
                valid = has_nonnegative_coefficients(horizontal);
            }
            if (!valid) {
                continue;
            }
            std::cout
                << "SU2_K3_SUM_CONE"
                << " name=" << name
                << " minimum_sum=" << minimum_sum
                << " x_threshold=" << x_threshold
                << " y_threshold=" << y_threshold
                << " result=PASS_COEFFICIENTS\n";
            return;
        }
    }
    throw std::runtime_error(
        name + " has no sum-cone certificate within threshold 30"
    );
}

}  // namespace

int main() {
    try {
        for (int q = 1; q <= 2; ++q) {
            for (int y = 0; y <= 7 * q; ++y) {
                const Integer margin =
                    multiplicity_six_integer(q, 0)
                        * multiplicity_seven_integer(q, y)
                    - multiplicity_seven_integer(q, 0)
                        * multiplicity_six_integer(q, y);
                if (margin < 0) {
                    throw std::runtime_error(
                        "small stable chamber is negative"
                    );
                }
            }
        }
        std::cout
            << "SU2_K3_STABLE_SMALL maximum_Q=2 result=PASS_EXACT\n";

        const Polynomial x = Polynomial::variable(0);
        const Polynomial y = Polynomial::variable(1);

        certify(
            "0_le_Y_le_Q_minus_3",
            stable_margin(
                x + y + constant(3),
                x,
                true,
                true,
                true,
                true,
                true,
                true,
                true
            )
        );
        certify_sum_cone(
            "Q_minus_2_le_Y_le_2Q_minus_2",
            stable_margin(
                x + y,
                scale(x, 2) + y - constant(2),
                true,
                true,
                true,
                true,
                true,
                true,
                false
            ),
            3
        );
        certify(
            "2Q_minus_1_le_Y_le_3Q_minus_2",
            stable_margin(
                x + y + constant(1),
                scale(x, 3) + scale(y, 2) + constant(1),
                true,
                true,
                false,
                true,
                true,
                true,
                false
            )
        );
        certify_sum_cone(
            "3Q_minus_1_le_Y_le_4Q_minus_1",
            stable_margin(
                x + y,
                scale(x, 4) + scale(y, 3) - constant(1),
                true,
                true,
                false,
                true,
                true,
                false,
                false
            ),
            3
        );
        certify(
            "4Q_le_Y_le_5Q_minus_1",
            stable_margin(
                x + y + constant(1),
                scale(x, 5) + scale(y, 4) + constant(4),
                true,
                false,
                false,
                true,
                true,
                false,
                false
            )
        );
        certify_sum_cone(
            "5Q_le_Y_le_6Q",
            stable_margin(
                x + y,
                scale(x, 6) + scale(y, 5),
                true,
                false,
                false,
                true,
                false,
                false,
                false
            ),
            3
        );
        certify(
            "6Q_plus_1_le_Y_le_7Q",
            stable_margin(
                x + y + constant(1),
                scale(x, 7) + scale(y, 6) + constant(7),
                false,
                false,
                false,
                true,
                false,
                false,
                false
            )
        );
        std::cout << "SU2_K3_STABLE_KERNEL result=PASS_EXACT\n";

        for (int q = 1; q <= 10; ++q) {
            for (int y_value = 0;
                 y_value <= 2 * q + 1;
                 ++y_value) {
                const Integer margin =
                    nearest_six_integer(q, 0)
                        * nearest_seven_integer(q, y_value)
                    - nearest_seven_integer(q, 0)
                        * nearest_six_integer(q, y_value);
                if (margin < 0) {
                    throw std::runtime_error(
                        "small nearest-wall chamber is negative"
                    );
                }
            }
        }
        std::cout
            << "SU2_K3_NEAREST_SMALL maximum_Q=10 result=PASS_EXACT\n";

        const Polynomial large_q = x + constant(11);
        for (int y_value = 0; y_value <= 3; ++y_value) {
            certify(
                "nearest_Y_" + std::to_string(y_value),
                nearest_margin(
                    large_q,
                    constant(y_value),
                    true,
                    false,
                    true,
                    true,
                    false,
                    true,
                    false
                )
            );
        }
        certify(
            "nearest_4_le_Y_le_Q_minus_5",
            nearest_margin(
                x + y + constant(11),
                x + constant(4),
                true,
                true,
                true,
                true,
                false,
                true,
                false
            )
        );
        for (int offset = -4; offset <= -3; ++offset) {
            certify(
                "nearest_Y_eq_Q_" + std::to_string(offset),
                nearest_margin(
                    large_q,
                    large_q + constant(offset),
                    true,
                    true,
                    true,
                    true,
                    false,
                    false,
                    false
                )
            );
        }
        for (int offset = -2; offset <= 4; ++offset) {
            certify(
                "nearest_Y_eq_Q_plus_" + std::to_string(offset),
                nearest_margin(
                    large_q,
                    large_q + constant(offset),
                    true,
                    true,
                    true,
                    false,
                    false,
                    false,
                    false
                )
            );
        }
        for (int offset = 5; offset <= 6; ++offset) {
            certify(
                "nearest_Y_eq_Q_plus_" + std::to_string(offset),
                nearest_margin(
                    large_q,
                    large_q + constant(offset),
                    true,
                    true,
                    true,
                    false,
                    true,
                    false,
                    false
                )
            );
        }
        certify(
            "nearest_Q_plus_7_le_Y_le_2Q_minus_4",
            nearest_margin(
                x + y + constant(11),
                scale(x, 2) + y + constant(18),
                true,
                true,
                true,
                false,
                true,
                false,
                true
            )
        );
        for (int offset = -3; offset <= -2; ++offset) {
            certify(
                "nearest_Y_eq_2Q_" + std::to_string(offset),
                nearest_margin(
                    large_q,
                    scale(large_q, 2) + constant(offset),
                    true,
                    true,
                    false,
                    false,
                    true,
                    false,
                    true
                )
            );
        }
        for (int offset = -1; offset <= 1; ++offset) {
            certify(
                "nearest_Y_eq_2Q_plus_" + std::to_string(offset),
                nearest_margin(
                    large_q,
                    scale(large_q, 2) + constant(offset),
                    false,
                    true,
                    false,
                    false,
                    true,
                    false,
                    true
                )
            );
        }
        std::cout << "SU2_K3_NEAREST_KERNEL result=PASS_EXACT\n";

        std::uint64_t recurrence_comparisons = 0U;
        for (int q_value = 1; q_value <= 40; ++q_value) {
            recurrence_comparisons += audit_formula_case(
                2 * q_value + 1,
                q_value,
                true
            );
            recurrence_comparisons += audit_formula_case(
                7 * q_value,
                q_value,
                false
            );
        }
        std::cout
            << "SU2_K3_FORMULA_AUDIT"
            << " maximum_Q=40"
            << " comparisons=" << recurrence_comparisons
            << " result=PASS_EXACT\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_K3_KERNEL FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
