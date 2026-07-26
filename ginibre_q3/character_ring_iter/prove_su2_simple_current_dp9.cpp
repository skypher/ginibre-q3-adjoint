#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using Exponent = std::pair<int, int>;

class Polynomial {
public:
    Polynomial() = default;

    explicit Polynomial(const Rational& value) {
        if (value != 0) {
            terms_[{0, 0}] = value;
        }
    }

    static Polynomial variable(int index) {
        Polynomial result;
        result.terms_[index == 0 ? Exponent{1, 0} : Exponent{0, 1}]
            = Rational(1);
        return result;
    }

    Polynomial& operator+=(const Polynomial& other) {
        for (const auto& [power, coefficient] : other.terms_) {
            terms_[power] += coefficient;
            if (terms_[power] == 0) {
                terms_.erase(power);
            }
        }
        return *this;
    }

    Polynomial& operator-=(const Polynomial& other) {
        for (const auto& [power, coefficient] : other.terms_) {
            terms_[power] -= coefficient;
            if (terms_[power] == 0) {
                terms_.erase(power);
            }
        }
        return *this;
    }

    Polynomial& operator*=(const Polynomial& other) {
        std::map<Exponent, Rational> product;
        for (const auto& [left_power, left_coefficient] : terms_) {
            for (const auto& [right_power, right_coefficient]
                 : other.terms_) {
                const Exponent power{
                    left_power.first + right_power.first,
                    left_power.second + right_power.second
                };
                product[power] += left_coefficient * right_coefficient;
            }
        }
        terms_.clear();
        for (const auto& [power, coefficient] : product) {
            if (coefficient != 0) {
                terms_[power] = coefficient;
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

Polynomial scale(const Polynomial& value, long numerator, long denominator = 1) {
    return value * Polynomial(Rational(numerator, denominator));
}

Polynomial power(Polynomial base, int exponent) {
    Polynomial result = constant(1);
    while (exponent > 0) {
        if ((exponent & 1) != 0) {
            result *= base;
        }
        exponent /= 2;
        if (exponent != 0) {
            base *= base;
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

Polynomial dp9(
    const Polynomial& q_half,
    const Polynomial& wall_distance,
    bool upper_half
) {
    const Polynomial& q = q_half;
    const Polynomial& c = wall_distance;

    const Polynomial g4 = q + c + constant(1);
    const Polynomial g5 =
        binomial(scale(q, 2) + c + constant(2), 2)
        - scale(binomial(c + constant(1), 2), 5);
    const Polynomial g6 =
        binomial(scale(q, 3) + c + constant(3), 3)
        - scale(binomial(q + c + constant(2), 3), 6);

    const Polynomial f4 = scale(q, 2) + constant(1);
    const Polynomial f5_stable =
        scale(
            scale(power(q, 2), 5) + scale(q, 5) + constant(2),
            1,
            2
        );
    const Polynomial f6_stable =
        scale(power(q, 3), 4)
        + scale(power(q, 2), 6)
        + scale(q, 4)
        + constant(1);
    const Polynomial f7_stable =
        binomial(scale(q, 7) + constant(5), 5)
        - scale(binomial(scale(q, 5) + constant(4), 5), 7)
        + scale(binomial(scale(q, 3) + constant(3), 5), 21)
        - scale(binomial(q + constant(2), 5), 35);

    Polynomial f5 = f5_stable;
    const Polynomial f6 =
        f6_stable - binomial(scale(c, 2) + constant(2), 3);
    Polynomial f7 =
        f7_stable
        - binomial(q + scale(c, 2) + constant(3), 4);
    if (upper_half) {
        f5 -= binomial(scale(c, 2) - q + constant(1), 2);
        f7 += scale(
            binomial(scale(c, 2) - q + constant(2), 4),
            7
        );
    }

    const Polynomial delta3 = g4 * f6 - f7;
    const Polynomial delta5 = g6 * f4 - g5 * f5;
    return scale(delta3, 2) + scale(delta5, 3);
}

Polynomial reverse_nine_band_five_six(
    const Polynomial& q_half,
    const Polynomial& wall_distance
) {
    const Polynomial& q = q_half;
    const Polynomial& a = wall_distance;
    const Polynomial g6 = binomial(a + constant(3), 3);
    const Polynomial g7 = binomial(q + a + constant(4), 4);
    const Polynomial g8 =
        binomial(scale(q, 2) + a + constant(5), 5)
        - scale(binomial(a + constant(4), 5), 8);
    const Polynomial g10 =
        binomial(scale(q, 4) + a + constant(7), 7)
        - scale(
            binomial(scale(q, 2) + a + constant(6), 7),
            10
        )
        + scale(binomial(a + constant(5), 7), 45);
    return scale(g6 * (scale(q, 2) + constant(1)), 126)
        + scale(g8, 36) - scale(g7, 84) + g10;
}

Polynomial reverse_nine_band_four_five(
    const Polynomial& q_half,
    const Polynomial& wall_distance
) {
    const Polynomial& q = q_half;
    const Polynomial& a = wall_distance;
    const Polynomial g5 = binomial(a + constant(2), 2);
    const Polynomial g6 = binomial(q + a + constant(3), 3);
    const Polynomial g7 =
        binomial(scale(q, 2) + a + constant(4), 4)
        - scale(binomial(a + constant(3), 4), 7);
    const Polynomial g8 =
        binomial(scale(q, 3) + a + constant(5), 5)
        - scale(
            binomial(q + a + constant(4), 5),
            8
        );
    const Polynomial g10 =
        binomial(scale(q, 5) + a + constant(7), 7)
        - scale(
            binomial(scale(q, 3) + a + constant(6), 7),
            10
        )
        + scale(
            binomial(q + a + constant(5), 7),
            45
        );
    const Polynomial f4 = scale(q, 2) + constant(1);
    const Polynomial f5 =
        scale(
            scale(power(q, 2), 5) + scale(q, 5) + constant(2),
            1,
            2
        );
    return scale(g6 * f4 - g5 * f5, 126)
        + scale(g8, 36) - scale(g7, 84) + g10;
}

Polynomial reverse_nine_band_three_four_lower(
    const Polynomial& q_half,
    const Polynomial& wall_distance
) {
    const Polynomial& q = q_half;
    const Polynomial& b = wall_distance;
    const Polynomial g4 = b + constant(1);
    const Polynomial g5 =
        binomial(q + b + constant(2), 2);
    const Polynomial g6 =
        binomial(scale(q, 2) + b + constant(3), 3)
        - scale(binomial(b + constant(2), 3), 6);
    const Polynomial g7 =
        binomial(scale(q, 3) + b + constant(4), 4)
        - scale(
            binomial(q + b + constant(3), 4),
            7
        );
    const Polynomial g8 =
        binomial(scale(q, 4) + b + constant(5), 5)
        - scale(
            binomial(scale(q, 2) + b + constant(4), 5),
            8
        )
        + scale(binomial(b + constant(3), 5), 28);
    // The omitted next affine image of g_10 is nonnegative.
    const Polynomial g10_lower =
        binomial(scale(q, 6) + b + constant(7), 7)
        - scale(
            binomial(scale(q, 4) + b + constant(6), 7),
            10
        )
        + scale(
            binomial(scale(q, 2) + b + constant(5), 7),
            45
        )
        - scale(binomial(b + constant(4), 7), 120);
    const Polynomial f4 = scale(q, 2) + constant(1);
    const Polynomial f5 =
        scale(
            scale(power(q, 2), 5) + scale(q, 5) + constant(2),
            1,
            2
        );
    const Polynomial f6 =
        scale(power(q, 3), 4)
        + scale(power(q, 2), 6)
        + scale(q, 4)
        + constant(1);
    return scale(g4 * f6, 84)
        + scale(g6 * f4 - g5 * f5, 126)
        + scale(g8, 36) - scale(g7, 84) + g10_lower;
}

Polynomial reverse_nine_strip_lower(
    const Polynomial& q_half,
    const Polynomial& wall_distance,
    bool include_f7_image,
    bool include_g7_image
) {
    const Polynomial& q = q_half;
    const Polynomial& c = wall_distance;
    const Polynomial g3 = constant(1);
    const Polynomial g4 = q + c + constant(1);
    const Polynomial g5 =
        binomial(scale(q, 2) + c + constant(2), 2)
        - scale(binomial(c + constant(1), 2), 5);
    const Polynomial g6 =
        binomial(scale(q, 3) + c + constant(3), 3)
        - scale(binomial(q + c + constant(2), 3), 6);
    Polynomial g7 =
        binomial(scale(q, 4) + c + constant(4), 4)
        - scale(
            binomial(scale(q, 2) + c + constant(3), 4),
            7
        )
        + scale(binomial(c + constant(2), 4), 21);
    if (include_g7_image) {
        g7 += binomial(
            scale(c, 3) - scale(q, 2) + constant(2),
            4
        );
    }
    // Omit the nonnegative next affine image of g_8.
    const Polynomial g8_lower =
        binomial(scale(q, 5) + c + constant(5), 5)
        - scale(
            binomial(scale(q, 3) + c + constant(4), 5),
            8
        )
        + scale(
            binomial(q + c + constant(3), 5),
            28
        );
    // Omit the nonnegative next affine image of g_10.
    const Polynomial g10_lower =
        binomial(scale(q, 7) + c + constant(7), 7)
        - scale(
            binomial(scale(q, 5) + c + constant(6), 7),
            10
        )
        + scale(
            binomial(scale(q, 3) + c + constant(5), 7),
            45
        )
        - scale(
            binomial(q + c + constant(4), 7),
            120
        );

    const Polynomial f4 = scale(q, 2) + constant(1);
    const Polynomial f5_upper =
        scale(
            scale(power(q, 2), 5) + scale(q, 5) + constant(2),
            1,
            2
        );
    const Polynomial f6 =
        scale(power(q, 3), 4)
        + scale(power(q, 2), 6)
        + scale(q, 4)
        + constant(1)
        - binomial(scale(c, 2) + constant(2), 3);
    Polynomial f7_upper =
        binomial(scale(q, 7) + constant(5), 5)
        - scale(binomial(scale(q, 5) + constant(4), 5), 7)
        + scale(binomial(scale(q, 3) + constant(3), 5), 21)
        - scale(binomial(q + constant(2), 5), 35);
    if (include_f7_image) {
        f7_upper += scale(
            binomial(scale(c, 2) - q + constant(2), 4),
            7
        );
    }

    return scale(g4 * f6 - g3 * f7_upper, 36)
        + scale(g6 * f4 - g5 * f5_upper, 126)
        + scale(g4 * f6, 48)
        + scale(g8_lower, 36) - scale(g7, 84) + g10_lower;
}

Polynomial delta_five_strip(
    const Polynomial& q_half,
    const Polynomial& wall_distance,
    bool include_wall_correction
) {
    const Polynomial& q = q_half;
    const Polynomial& c = wall_distance;
    const Polynomial g5 =
        binomial(scale(q, 2) + c + constant(2), 2)
        - scale(binomial(c + constant(1), 2), 5);
    const Polynomial g6 =
        binomial(scale(q, 3) + c + constant(3), 3)
        - scale(binomial(q + c + constant(2), 3), 6);
    const Polynomial f4 = scale(q, 2) + constant(1);
    Polynomial f5 =
        scale(
            scale(power(q, 2), 5) + scale(q, 5) + constant(2),
            1,
            2
        );
    if (include_wall_correction) {
        f5 -= binomial(scale(c, 2) - q + constant(1), 2);
    }
    return g6 * f4 - g5 * f5;
}

Polynomial eleven_core_strip_lower(
    const Polynomial& q_half,
    const Polynomial& wall_distance,
    bool include_f7_image,
    bool include_g7_image,
    bool include_g8_image,
    bool include_f5_image,
    bool include_f9_u_image,
    bool include_f9_second_image
) {
    const Polynomial& q = q_half;
    const Polynomial& c = wall_distance;
    const Polynomial g3 = constant(1);
    const Polynomial g4 = q + c + constant(1);
    const Polynomial g5 =
        binomial(scale(q, 2) + c + constant(2), 2)
        - scale(binomial(c + constant(1), 2), 5);
    const Polynomial g6 =
        binomial(scale(q, 3) + c + constant(3), 3)
        - scale(binomial(q + c + constant(2), 3), 6);
    Polynomial g7 =
        binomial(scale(q, 4) + c + constant(4), 4)
        - scale(
            binomial(scale(q, 2) + c + constant(3), 4),
            7
        )
        + scale(binomial(c + constant(2), 4), 21);
    if (include_g7_image) {
        g7 += binomial(
            scale(c, 3) - scale(q, 2) + constant(2),
            4
        );
    }
    Polynomial g8_lower =
        binomial(scale(q, 5) + c + constant(5), 5)
        - scale(
            binomial(scale(q, 3) + c + constant(4), 5),
            8
        )
        + scale(
            binomial(q + c + constant(3), 5),
            28
        );
    if (include_g8_image) {
        g8_lower += binomial(
            scale(c, 3) - q + constant(3),
            5
        );
    }

    const Polynomial f4 = scale(q, 2) + constant(1);
    Polynomial f5_upper =
        scale(
            scale(power(q, 2), 5) + scale(q, 5) + constant(2),
            1,
            2
        );
    if (include_f5_image) {
        f5_upper -= binomial(
            scale(c, 2) - q + constant(1),
            2
        );
    }
    const Polynomial f6 =
        scale(power(q, 3), 4)
        + scale(power(q, 2), 6)
        + scale(q, 4)
        + constant(1)
        - binomial(scale(c, 2) + constant(2), 3);
    Polynomial f7_upper =
        binomial(scale(q, 7) + constant(5), 5)
        - scale(binomial(scale(q, 5) + constant(4), 5), 7)
        + scale(binomial(scale(q, 3) + constant(3), 5), 21)
        - scale(binomial(q + constant(2), 5), 35);
    f7_upper -= binomial(
        q + scale(c, 2) + constant(3),
        4
    );
    if (include_f7_image) {
        f7_upper += scale(
            binomial(scale(c, 2) - q + constant(2), 4),
            7
        );
    }
    const Polynomial f8_stable =
        binomial(scale(q, 8) + constant(6), 6)
        - scale(binomial(scale(q, 6) + constant(5), 6), 8)
        + scale(binomial(scale(q, 4) + constant(4), 6), 28)
        - scale(binomial(scale(q, 2) + constant(3), 6), 56);
    const Polynomial f8 =
        f8_stable
        - binomial(
            scale(q, 2) + scale(c, 2) + constant(4),
            5
        )
        + scale(binomial(scale(c, 2) + constant(3), 5), 8);
    Polynomial f9_upper =
        binomial(scale(q, 9) + constant(7), 7)
        - scale(binomial(scale(q, 7) + constant(6), 7), 9)
        + scale(binomial(scale(q, 5) + constant(5), 7), 36)
        - scale(binomial(scale(q, 3) + constant(4), 7), 84)
        + scale(binomial(q + constant(3), 7), 126);
    f9_upper -= binomial(
        scale(q, 3) + scale(c, 2) + constant(5),
        6
    );
    f9_upper += scale(
        binomial(q + scale(c, 2) + constant(4), 6),
        9
    );
    if (include_f9_u_image) {
        f9_upper -= scale(
            binomial(scale(c, 2) - q + constant(3), 6),
            36
        );
    }
    if (include_f9_second_image) {
        f9_upper -= binomial(
            scale(c, 4) - scale(q, 3) + constant(3),
            6
        );
    }

    const Polynomial delta7 =
        g8_lower * f4 - g7 * f5_upper;
    const Polynomial delta5 =
        g6 * f6 - g5 * f7_upper;
    const Polynomial delta3 =
        g4 * f8 - g3 * f9_upper;
    return scale(delta7, 10)
        + scale(delta5, 14)
        + scale(delta3, 5);
}

long integer_binomial(int top, int bottom) {
    if (bottom < 0 || bottom > top) {
        return 0;
    }
    long result = 1;
    for (int index = 1; index <= bottom; ++index) {
        result = result * (top - bottom + index) / index;
    }
    return result;
}

Polynomial wall_primary_band(
    const Polynomial& q_half,
    const Polynomial& wall_distance,
    int power_value,
    int band
) {
    Polynomial result = constant(0);
    for (int image = 0;
         power_value - band - 2 * image >= 0;
         ++image) {
        const Polynomial top =
            scale(q_half, power_value - band - 2 * image)
            + wall_distance
            + constant(power_value - 3 - image);
        const long coefficient =
            integer_binomial(power_value, image);
        result += scale(
            binomial(top, power_value - 3),
            (image & 1) == 0 ? coefficient : -coefficient
        );
    }
    return result;
}

Polynomial eleven_core_far_band_lower(
    const Polynomial& q_half,
    const Polynomial& wall_distance,
    int band
) {
    const Polynomial g5 = wall_primary_band(
        q_half, wall_distance, 5, band
    );
    const Polynomial g6 = wall_primary_band(
        q_half, wall_distance, 6, band
    );
    const Polynomial g7 = wall_primary_band(
        q_half, wall_distance, 7, band
    );
    const Polynomial g8 = wall_primary_band(
        q_half, wall_distance, 8, band
    );
    const Polynomial f4 = scale(q_half, 2) + constant(1);
    const Polynomial f5 =
        scale(
            scale(power(q_half, 2), 5)
            + scale(q_half, 5)
            + constant(2),
            1,
            2
        );
    const Polynomial f6 =
        scale(power(q_half, 3), 4)
        + scale(power(q_half, 2), 6)
        + scale(q_half, 4)
        + constant(1);
    const Polynomial f7_stable =
        binomial(scale(q_half, 7) + constant(5), 5)
        - scale(
            binomial(scale(q_half, 5) + constant(4), 5),
            7
        )
        + scale(
            binomial(scale(q_half, 3) + constant(3), 5),
            21
        )
        - scale(binomial(q_half + constant(2), 5), 35);
    const Polynomial delta7 = g8 * f4 - g7 * f5;
    const Polynomial delta5 = g6 * f6 - g5 * f7_stable;
    return scale(delta7, 10) + scale(delta5, 14);
}

Polynomial eleven_delta_seven_strip_lower(
    const Polynomial& q,
    const Polynomial& c,
    bool include_g7_image,
    bool include_g8_image
) {
    Polynomial g7 =
        binomial(scale(q, 4) + c + constant(4), 4)
        - scale(
            binomial(scale(q, 2) + c + constant(3), 4),
            7
        )
        + scale(binomial(c + constant(2), 4), 21);
    if (include_g7_image) {
        g7 += binomial(
            scale(c, 3) - scale(q, 2) + constant(2),
            4
        );
    }
    Polynomial g8_lower =
        binomial(scale(q, 5) + c + constant(5), 5)
        - scale(
            binomial(scale(q, 3) + c + constant(4), 5),
            8
        )
        + scale(binomial(q + c + constant(3), 5), 28);
    if (include_g8_image) {
        g8_lower += binomial(
            scale(c, 3) - q + constant(3),
            5
        );
    }
    const Polynomial f4 = scale(q, 2) + constant(1);
    const Polynomial f5_upper =
        scale(
            scale(power(q, 2), 5) + scale(q, 5) + constant(2),
            1,
            2
        );
    return g8_lower * f4 - g7 * f5_upper;
}

Polynomial eleven_delta_five_strip(
    const Polynomial& q,
    const Polynomial& c,
    bool include_f7_image
) {
    const Polynomial g5 =
        binomial(scale(q, 2) + c + constant(2), 2)
        - scale(binomial(c + constant(1), 2), 5);
    const Polynomial g6 =
        binomial(scale(q, 3) + c + constant(3), 3)
        - scale(binomial(q + c + constant(2), 3), 6);
    const Polynomial f6 =
        scale(power(q, 3), 4)
        + scale(power(q, 2), 6)
        + scale(q, 4)
        + constant(1)
        - binomial(scale(c, 2) + constant(2), 3);
    Polynomial f7 =
        binomial(scale(q, 7) + constant(5), 5)
        - scale(binomial(scale(q, 5) + constant(4), 5), 7)
        + scale(binomial(scale(q, 3) + constant(3), 5), 21)
        - scale(binomial(q + constant(2), 5), 35)
        - binomial(q + scale(c, 2) + constant(3), 4);
    if (include_f7_image) {
        f7 += scale(
            binomial(scale(c, 2) - q + constant(2), 4),
            7
        );
    }
    return g6 * f6 - g5 * f7;
}

Polynomial eleven_delta_far_band(
    const Polynomial& q,
    const Polynomial& a,
    int band,
    bool seventh
) {
    const Polynomial f4 = scale(q, 2) + constant(1);
    const Polynomial f5 =
        scale(
            scale(power(q, 2), 5) + scale(q, 5) + constant(2),
            1,
            2
        );
    const Polynomial f6 =
        scale(power(q, 3), 4)
        + scale(power(q, 2), 6)
        + scale(q, 4)
        + constant(1);
    const Polynomial f7 =
        binomial(scale(q, 7) + constant(5), 5)
        - scale(binomial(scale(q, 5) + constant(4), 5), 7)
        + scale(binomial(scale(q, 3) + constant(3), 5), 21)
        - scale(binomial(q + constant(2), 5), 35);
    if (seventh) {
        return wall_primary_band(q, a, 8, band) * f4
            - wall_primary_band(q, a, 7, band) * f5;
    }
    return wall_primary_band(q, a, 6, band) * f6
        - wall_primary_band(q, a, 5, band) * f7;
}

Rational evaluate(
    const Polynomial& value,
    const Integer& first,
    const Integer& second
) {
    Rational result(0);
    for (const auto& [exponent, coefficient] : value.terms()) {
        Integer monomial = 1;
        for (int index = 0; index < exponent.first; ++index) {
            monomial *= first;
        }
        for (int index = 0; index < exponent.second; ++index) {
            monomial *= second;
        }
        result += coefficient * monomial;
    }
    return result;
}

bool fuses(int level, int first, int second, int output) {
    return std::abs(first - second) <= output
        && output <= std::min(
            first + second,
            2 * level - first - second
        )
        && ((first + second + output) & 1) == 0;
}

void multiply(int level, int label, std::vector<Integer>& state) {
    std::vector<Integer> next(static_cast<std::size_t>(level + 1));
    for (int source = 0; source <= level; ++source) {
        if (state[static_cast<std::size_t>(source)] == 0) {
            continue;
        }
        for (int output = 0; output <= level; ++output) {
            if (fuses(level, label, source, output)) {
                next[static_cast<std::size_t>(output)]
                    += state[static_cast<std::size_t>(source)];
            }
        }
    }
    state.swap(next);
}

Integer direct_dp9(int level, int label) {
    std::vector<Integer> state(static_cast<std::size_t>(level + 1));
    std::vector<Integer> closed(8U);
    std::vector<Integer> wall(7U);
    state[0] = 1;
    for (int power_index = 0; power_index <= 7; ++power_index) {
        closed[static_cast<std::size_t>(power_index)] = state[0];
        if (power_index <= 6) {
            wall[static_cast<std::size_t>(power_index)] =
                state[static_cast<std::size_t>(level)];
        }
        if (power_index != 7) {
            multiply(level, label, state);
        }
    }
    const Integer delta3 =
        wall[4] * closed[6] - wall[3] * closed[7];
    const Integer delta5 =
        wall[6] * closed[4] - wall[5] * closed[5];
    return 2 * delta3 + 3 * delta5;
}

Integer direct_reverse_nine(int level, int label) {
    std::vector<Integer> state(static_cast<std::size_t>(level + 1));
    std::vector<Integer> closed(11U);
    std::vector<Integer> wall(11U);
    state[0] = 1;
    for (int power_index = 0; power_index <= 10; ++power_index) {
        closed[static_cast<std::size_t>(power_index)] = state[0];
        wall[static_cast<std::size_t>(power_index)] =
            state[static_cast<std::size_t>(level)];
        if (power_index != 10) {
            multiply(level, label, state);
        }
    }
    const Integer delta5 =
        wall[6] * closed[4] - wall[5] * closed[5];
    return 84 * wall[4] * closed[6]
        - 36 * wall[3] * closed[7]
        + 126 * delta5
        + 36 * wall[8] - 84 * wall[7] + wall[10];
}

bool report(const char* name, const Polynomial& value) {
    bool nonnegative = true;
    std::cout << name << " terms=" << value.terms().size() << '\n';
    for (const auto& [power, coefficient] : value.terms()) {
        std::cout
            << "  x^" << power.first
            << " y^" << power.second
            << " coefficient=" << coefficient << '\n';
        nonnegative = nonnegative && coefficient >= 0;
    }
    std::cout << name << " coefficientwise_nonnegative="
              << (nonnegative ? "true" : "false") << '\n';
    return nonnegative;
}

bool audit_recurrence(
    const Polynomial& lower,
    const Polynomial& upper,
    const Polynomial& reverse_low,
    const Polynomial& reverse_unit,
    const Polynomial& reverse_middle,
    const Polynomial& reverse_high,
    const Polynomial& reverse_five_six,
    const Polynomial& reverse_at_five,
    const Polynomial& reverse_four_five,
    const Polynomial& reverse_three_four,
    const Polynomial& reverse_at_three,
    int maximum_level
) {
    std::size_t rows = 0U;
    for (int level = 4; level <= maximum_level; level += 2) {
        for (int label = 2; 2 * label < level; label += 2) {
            if (!(2 * label < level && level < 3 * label)) {
                continue;
            }
            const int q = label / 2;
            const int k = level / 2;
            const int c = 3 * q - k;
            Rational formula;
            if (2 * c <= q) {
                formula = evaluate(lower, c, q - 2 * c);
            } else {
                formula = evaluate(upper, q - c, 2 * c - q);
            }
            const Integer direct = direct_dp9(level, label);
            ++rows;
            if (formula != Rational(direct)) {
                std::cout
                    << "DP9_RECURRENCE_MISMATCH"
                    << " level=" << level
                    << " label=" << label
                    << " formula=" << formula
                    << " direct=" << direct << '\n';
                return false;
            }

            Rational reverse_bound;
            const int u = 2 * c - q;
            const int r = q - c;
            if (u <= 0) {
                reverse_bound = evaluate(
                    reverse_low, c - 1, q - 2 * c
                );
            } else if (u == 1) {
                reverse_bound = evaluate(reverse_unit, c - 2, 0);
            } else if (u - r <= 1) {
                reverse_bound = evaluate(
                    reverse_middle, u - 2, r - u + 1
                );
            } else {
                reverse_bound = evaluate(
                    reverse_high, r - 1, u - r - 2
                );
            }
            const Integer reverse_direct =
                direct_reverse_nine(level, label);
            if (reverse_bound < 0
                || reverse_bound > Rational(reverse_direct)) {
                std::cout
                    << "REVERSE9_LOWER_BOUND_MISMATCH"
                    << " level=" << level
                    << " label=" << label
                    << " lower=" << reverse_bound
                    << " direct=" << reverse_direct << '\n';
                return false;
            }
        }
    }
    std::cout
        << "DP9_RECURRENCE_AUDIT"
        << " maximum_level=" << maximum_level
        << " rows=" << rows
        << " dp9_mismatches=0"
        << " reverse_lower_mismatches=0"
        << " result=PASS_DISCOVERY\n";

    std::size_t far_rows = 0U;
    for (int level = 4; level <= maximum_level; level += 2) {
        for (int label = 2; 2 * label < level; label += 2) {
            const int q = label / 2;
            const int k = level / 2;
            if (k < 3 * q || k > 6 * q) {
                continue;
            }
            Rational lower_bound;
            if (k > 5 * q) {
                const int a = 6 * q - k;
                lower_bound = evaluate(
                    reverse_five_six, a, q - a - 1
                );
            } else if (k == 5 * q) {
                lower_bound = evaluate(reverse_at_five, q - 1, 0);
            } else if (k > 4 * q) {
                const int a = 5 * q - k;
                lower_bound = evaluate(
                    reverse_four_five, a, q - a - 1
                );
            } else if (k > 3 * q) {
                const int b = 4 * q - k;
                lower_bound = evaluate(
                    reverse_three_four, b, q - b - 1
                );
            } else {
                lower_bound = evaluate(reverse_at_three, q - 1, 0);
            }
            const Integer direct = direct_reverse_nine(level, label);
            ++far_rows;
            if (lower_bound < 0
                || lower_bound > Rational(direct)) {
                std::cout
                    << "REVERSE9_FAR_LOWER_BOUND_MISMATCH"
                    << " level=" << level
                    << " label=" << label
                    << " lower=" << lower_bound
                    << " direct=" << direct << '\n';
                return false;
            }
        }
    }
    std::cout
        << "REVERSE9_FAR_RECURRENCE_AUDIT"
        << " maximum_level=" << maximum_level
        << " rows=" << far_rows
        << " lower_mismatches=0 result=PASS_DISCOVERY\n";
    return true;
}

}  // namespace

int main() {
    try {
        const Polynomial x = Polynomial::variable(0);
        const Polynomial y = Polynomial::variable(1);

        // Lower chamber: Q=2x+y, c=x, where x>=1 and y>=0.
        const Polynomial lower = dp9(
            scale(x, 2) + y,
            x,
            false
        );

        // Upper chamber: Q=2x+y, c=x+y, where x>=1 and y>=0.
        const Polynomial upper = dp9(
            scale(x, 2) + y,
            x + y,
            true
        );

        // Reverse (9,1), 5q<k<=6q:
        // a=6Q-K>=0 and Q=a+y+1.
        const Polynomial reverse_five_six =
            reverse_nine_band_five_six(
                x + y + constant(1),
                x
            );
        const Polynomial boundary_q = x + constant(1);
        const Polynomial reverse_at_five =
            reverse_nine_band_five_six(boundary_q, boundary_q)
            - scale(
                scale(
                    scale(power(boundary_q, 2), 5)
                    + scale(boundary_q, 5)
                    + constant(2),
                    1,
                    2
                ),
                126
            );
        const Polynomial reverse_four_five =
            reverse_nine_band_four_five(
                x + y + constant(1),
                x
            );
        const Polynomial reverse_three_four_lower =
            reverse_nine_band_three_four_lower(
                x + y + constant(1),
                x
            );
        const Polynomial boundary_three_q = x + constant(1);
        const Polynomial boundary_three_f7 = scale(
            scale(power(boundary_three_q, 4), 153)
            + scale(power(boundary_three_q, 3), 302)
            + scale(power(boundary_three_q, 2), 255)
            + scale(boundary_three_q, 106)
            + constant(24),
            1,
            24
        );
        const Polynomial reverse_at_three =
            reverse_nine_band_three_four_lower(
                boundary_three_q,
                boundary_three_q
            )
            + binomial(boundary_three_q + constant(5), 7)
            - scale(boundary_three_f7, 36);
        const Polynomial reverse_strip_low =
            reverse_nine_strip_lower(
                scale(x + constant(1), 2) + y,
                x + constant(1),
                false,
                false
            );
        const Polynomial reverse_strip_unit =
            reverse_nine_strip_lower(
                scale(x, 2) + constant(3),
                x + constant(2),
                false,
                false
            );
        const Polynomial reverse_strip_middle =
            reverse_nine_strip_lower(
                scale(x, 3) + scale(y, 2) + constant(4),
                scale(x, 2) + y + constant(3),
                true,
                false
            );
        const Polynomial reverse_strip_high =
            reverse_nine_strip_lower(
                scale(x, 3) + y + constant(5),
                scale(x, 2) + y + constant(4),
                true,
                true
            );
        const Polynomial delta_five_low =
            delta_five_strip(
                scale(x + constant(1), 2) + y,
                x + constant(1),
                false
            );
        const Polynomial delta_five_unit =
            delta_five_strip(
                scale(x, 2) + constant(3),
                x + constant(2),
                true
            );
        const Polynomial delta_five_high =
            delta_five_strip(
                scale(x, 2) + y + constant(4),
                x + y + constant(3),
                true
            );
        const Polynomial eleven_core_low =
            eleven_core_strip_lower(
                scale(x + constant(1), 2) + y,
                x + constant(1),
                false,
                false,
                false,
                false,
                false,
                false
            );
        const Polynomial eleven_core_unit =
            eleven_core_strip_lower(
                scale(x, 2) + constant(3),
                x + constant(2),
                false,
                false,
                true,
                true,
                false,
                false
            );
        const Polynomial eleven_core_middle =
            eleven_core_strip_lower(
                scale(x, 3) + scale(y, 2) + constant(4),
                scale(x, 2) + y + constant(3),
                true,
                false,
                true,
                true,
                false,
                false
            );
        const Polynomial eleven_core_high =
            eleven_core_strip_lower(
                scale(x, 3) + y + constant(5),
                scale(x, 2) + y + constant(4),
                true,
                true,
                true,
                true,
                false,
                false
            );
        const Polynomial eleven_core_boundary_three =
            eleven_core_strip_lower(
                x + constant(1),
                constant(0),
                false,
                false,
                false,
                false,
                false,
                false
            );
        const Polynomial eleven_core_band_three_four =
            eleven_core_far_band_lower(
                x + y + constant(1),
                x,
                4
            );
        const Polynomial eleven_core_band_four_five =
            eleven_core_far_band_lower(
                x + y + constant(1),
                x,
                5
            );
        const Polynomial eleven_core_band_five_six =
            eleven_core_far_band_lower(
                x + y + constant(1),
                x,
                6
            );
        const Polynomial eleven_core_band_six_seven =
            eleven_core_far_band_lower(
                x + y + constant(1),
                x,
                7
            );
        const Polynomial eleven_core_band_seven_eight =
            eleven_core_far_band_lower(
                x + y + constant(1),
                x,
                8
            );
        const Polynomial eleven_delta7_low =
            eleven_delta_seven_strip_lower(
                scale(x + constant(1), 2) + y,
                x + constant(1),
                false,
                false
            );
        const Polynomial eleven_delta7_unit =
            eleven_delta_seven_strip_lower(
                scale(x, 2) + constant(3),
                x + constant(2),
                false,
                false
            );
        const Polynomial eleven_delta7_middle =
            eleven_delta_seven_strip_lower(
                scale(x, 3) + scale(y, 2) + constant(4),
                scale(x, 2) + y + constant(3),
                false,
                false
            );
        const Polynomial eleven_delta7_high =
            eleven_delta_seven_strip_lower(
                scale(x, 3) + y + constant(5),
                scale(x, 2) + y + constant(4),
                true,
                true
            );
        const Polynomial eleven_delta5_unit =
            eleven_delta_five_strip(
                scale(x, 2) + constant(3),
                x + constant(2),
                false
            );
        const Polynomial eleven_delta5_middle =
            eleven_delta_five_strip(
                scale(x, 3) + scale(y, 2) + constant(4),
                scale(x, 2) + y + constant(3),
                true
            );
        const Polynomial eleven_delta5_high =
            eleven_delta_five_strip(
                scale(x, 3) + y + constant(5),
                scale(x, 2) + y + constant(4),
                true
            );
        const Polynomial eleven_delta5_low_c1_zero =
            eleven_delta_five_strip(
                constant(2),
                constant(1),
                false
            );
        const Polynomial eleven_delta5_low_c1_tail =
            eleven_delta_five_strip(
                x + constant(3),
                constant(1),
                false
            );
        const Polynomial eleven_delta5_low_cge2 =
            eleven_delta_five_strip(
                scale(x + constant(2), 2) + y,
                x + constant(2),
                false
            );
        const Polynomial eleven_delta7_boundary_three =
            eleven_delta_seven_strip_lower(
                x + constant(1),
                constant(0),
                false,
                false
            );
        const Polynomial eleven_delta5_boundary_three =
            eleven_delta_five_strip(
                x + constant(1),
                constant(0),
                false
            );
        const Polynomial far_q = x + y + constant(1);
        const Polynomial eleven_delta7_band_three_four =
            eleven_delta_far_band(far_q, x, 4, true);
        const Polynomial eleven_delta7_band_four_five =
            eleven_delta_far_band(far_q, x, 5, true);
        const Polynomial eleven_delta7_band_five_six =
            eleven_delta_far_band(far_q, x, 6, true);
        const Polynomial eleven_delta7_band_six_seven =
            eleven_delta_far_band(far_q, x, 7, true);
        const Polynomial eleven_delta7_band_seven_eight =
            eleven_delta_far_band(far_q, x, 8, true);
        const Polynomial eleven_delta5_band_three_four =
            eleven_delta_far_band(far_q, x, 4, false);
        const Polynomial eleven_delta5_band_four_five =
            eleven_delta_far_band(far_q, x, 5, false);
        const Polynomial eleven_delta5_band_five_six =
            eleven_delta_far_band(far_q, x, 6, false);
        const Polynomial eleven_delta5_band_six_seven =
            eleven_delta_far_band(far_q, x, 7, false);
        const Polynomial eleven_delta5_band_seven_eight =
            eleven_delta_far_band(far_q, x, 8, false);

        const bool pass =
            report("DP9_LOWER", lower)
            && report("DP9_UPPER", upper)
            && report("REVERSE9_BAND_5_6", reverse_five_six)
            && report("REVERSE9_BOUNDARY_5", reverse_at_five)
            && report("REVERSE9_BAND_4_5", reverse_four_five)
            && report(
                "REVERSE9_BAND_3_4_LOWER",
                reverse_three_four_lower
            )
            && report("REVERSE9_BOUNDARY_3", reverse_at_three)
            && report("REVERSE9_STRIP_LOW", reverse_strip_low)
            && report("REVERSE9_STRIP_UNIT", reverse_strip_unit)
            && report("REVERSE9_STRIP_MIDDLE", reverse_strip_middle)
            && report("REVERSE9_STRIP_HIGH", reverse_strip_high)
            && report("DELTA5_STRIP_LOW", delta_five_low)
            && report("DELTA5_STRIP_UNIT", delta_five_unit)
            && report("DELTA5_STRIP_HIGH", delta_five_high)
            && report("ELEVEN_CORE_STRIP_LOW", eleven_core_low)
            && report("ELEVEN_CORE_STRIP_UNIT", eleven_core_unit)
            && report("ELEVEN_CORE_STRIP_MIDDLE", eleven_core_middle)
            && report("ELEVEN_CORE_STRIP_HIGH", eleven_core_high)
            && report(
                "ELEVEN_CORE_BOUNDARY_3",
                eleven_core_boundary_three
            )
            && report(
                "ELEVEN_CORE_BAND_3_4",
                eleven_core_band_three_four
            )
            && report(
                "ELEVEN_CORE_BAND_4_5",
                eleven_core_band_four_five
            )
            && report(
                "ELEVEN_CORE_BAND_5_6",
                eleven_core_band_five_six
            )
            && report(
                "ELEVEN_CORE_BAND_6_7",
                eleven_core_band_six_seven
            )
            && report(
                "ELEVEN_CORE_BAND_7_8",
                eleven_core_band_seven_eight
            )
            && report("ELEVEN_DELTA7_LOW", eleven_delta7_low)
            && report("ELEVEN_DELTA7_UNIT", eleven_delta7_unit)
            && report("ELEVEN_DELTA7_MIDDLE", eleven_delta7_middle)
            && report("ELEVEN_DELTA7_HIGH", eleven_delta7_high)
            && report(
                "ELEVEN_DELTA5_LOW_C1_ZERO",
                eleven_delta5_low_c1_zero
            )
            && report(
                "ELEVEN_DELTA5_LOW_C1_TAIL",
                eleven_delta5_low_c1_tail
            )
            && report(
                "ELEVEN_DELTA5_LOW_CGE2",
                eleven_delta5_low_cge2
            )
            && report("ELEVEN_DELTA5_UNIT", eleven_delta5_unit)
            && report("ELEVEN_DELTA5_MIDDLE", eleven_delta5_middle)
            && report("ELEVEN_DELTA5_HIGH", eleven_delta5_high)
            && report(
                "ELEVEN_DELTA7_BOUNDARY_3",
                eleven_delta7_boundary_three
            )
            && report(
                "ELEVEN_DELTA5_BOUNDARY_3",
                eleven_delta5_boundary_three
            )
            && report(
                "ELEVEN_DELTA7_BAND_3_4",
                eleven_delta7_band_three_four
            )
            && report(
                "ELEVEN_DELTA7_BAND_4_5",
                eleven_delta7_band_four_five
            )
            && report(
                "ELEVEN_DELTA7_BAND_5_6",
                eleven_delta7_band_five_six
            )
            && report(
                "ELEVEN_DELTA7_BAND_6_7",
                eleven_delta7_band_six_seven
            )
            && report(
                "ELEVEN_DELTA7_BAND_7_8",
                eleven_delta7_band_seven_eight
            )
            && report(
                "ELEVEN_DELTA5_BAND_3_4",
                eleven_delta5_band_three_four
            )
            && report(
                "ELEVEN_DELTA5_BAND_4_5",
                eleven_delta5_band_four_five
            )
            && report(
                "ELEVEN_DELTA5_BAND_5_6",
                eleven_delta5_band_five_six
            )
            && report(
                "ELEVEN_DELTA5_BAND_6_7",
                eleven_delta5_band_six_seven
            )
            && report(
                "ELEVEN_DELTA5_BAND_7_8",
                eleven_delta5_band_seven_eight
            )
            && audit_recurrence(
                lower,
                upper,
                reverse_strip_low,
                reverse_strip_unit,
                reverse_strip_middle,
                reverse_strip_high,
                reverse_five_six,
                reverse_at_five,
                reverse_four_five,
                reverse_three_four_lower,
                reverse_at_three,
                100
            );
        std::cout << "SU2_SIMPLE_CURRENT_DP9 result="
                  << (pass ? "PASS_EXACT" : "FAIL") << '\n';
        return pass ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SIMPLE_CURRENT_DP9 FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
