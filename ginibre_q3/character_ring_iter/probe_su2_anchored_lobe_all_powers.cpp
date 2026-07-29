#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;

class Polynomial {
public:
    Polynomial() = default;

    explicit Polynomial(const Rational& constant) {
        if (constant != 0) {
            coefficients_.push_back(constant);
        }
    }

    explicit Polynomial(std::vector<Rational> coefficients)
        : coefficients_(std::move(coefficients)) {
        trim();
    }

    const std::vector<Rational>& coefficients() const {
        return coefficients_;
    }

    Polynomial derivative() const {
        if (coefficients_.size() < 2) {
            return Polynomial{};
        }
        std::vector<Rational> result(coefficients_.size() - 1);
        for (std::size_t degree = 1; degree < coefficients_.size();
             ++degree) {
            result[degree - 1] =
                Rational{static_cast<long long>(degree)}
                * coefficients_[degree];
        }
        return Polynomial{std::move(result)};
    }

    Polynomial& operator+=(const Polynomial& other) {
        coefficients_.resize(
            std::max(coefficients_.size(), other.coefficients_.size())
        );
        for (std::size_t degree = 0;
             degree < other.coefficients_.size();
             ++degree) {
            coefficients_[degree] += other.coefficients_[degree];
        }
        trim();
        return *this;
    }

    Polynomial& operator-=(const Polynomial& other) {
        coefficients_.resize(
            std::max(coefficients_.size(), other.coefficients_.size())
        );
        for (std::size_t degree = 0;
             degree < other.coefficients_.size();
             ++degree) {
            coefficients_[degree] -= other.coefficients_[degree];
        }
        trim();
        return *this;
    }

private:
    void trim() {
        while (!coefficients_.empty() && coefficients_.back() == 0) {
            coefficients_.pop_back();
        }
    }

    std::vector<Rational> coefficients_;
};

Polynomial operator+(Polynomial left, const Polynomial& right) {
    left += right;
    return left;
}

Polynomial operator-(Polynomial left, const Polynomial& right) {
    left -= right;
    return left;
}

Polynomial operator*(const Polynomial& left, const Polynomial& right) {
    if (left.coefficients().empty() || right.coefficients().empty()) {
        return Polynomial{};
    }
    std::vector<Rational> result(
        left.coefficients().size() + right.coefficients().size() - 1
    );
    for (std::size_t i = 0; i < left.coefficients().size(); ++i) {
        for (std::size_t j = 0; j < right.coefficients().size(); ++j) {
            result[i + j] +=
                left.coefficients()[i] * right.coefficients()[j];
        }
    }
    return Polynomial{std::move(result)};
}

Polynomial operator*(const Rational& scalar, Polynomial value) {
    std::vector<Rational> result = value.coefficients();
    for (Rational& coefficient : result) {
        coefficient *= scalar;
    }
    return Polynomial{std::move(result)};
}

Polynomial one_plus_z_power(int exponent) {
    Polynomial result{Rational{1}};
    const Polynomial one_plus_z{
        std::vector<Rational>{Rational{1}, Rational{1}}
    };
    for (int index = 0; index < exponent; ++index) {
        result = result * one_plus_z;
    }
    return result;
}

struct RationalPolynomial {
    Polynomial numerator;
    int denominator_power = 0;
};

RationalPolynomial align_denominator(
    RationalPolynomial value,
    int denominator_power
) {
    if (value.denominator_power < denominator_power) {
        value.numerator =
            value.numerator
            * one_plus_z_power(
                denominator_power - value.denominator_power
            );
        value.denominator_power = denominator_power;
    }
    return value;
}

RationalPolynomial operator+(
    RationalPolynomial left,
    RationalPolynomial right
) {
    const int denominator_power =
        std::max(left.denominator_power, right.denominator_power);
    left = align_denominator(std::move(left), denominator_power);
    right = align_denominator(std::move(right), denominator_power);
    return RationalPolynomial{
        left.numerator + right.numerator,
        denominator_power
    };
}

RationalPolynomial operator-(
    RationalPolynomial left,
    RationalPolynomial right
) {
    const int denominator_power =
        std::max(left.denominator_power, right.denominator_power);
    left = align_denominator(std::move(left), denominator_power);
    right = align_denominator(std::move(right), denominator_power);
    return RationalPolynomial{
        left.numerator - right.numerator,
        denominator_power
    };
}

RationalPolynomial operator*(
    const RationalPolynomial& left,
    const RationalPolynomial& right
) {
    return RationalPolynomial{
        left.numerator * right.numerator,
        left.denominator_power + right.denominator_power
    };
}

RationalPolynomial operator*(
    const Rational& scalar,
    RationalPolynomial value
) {
    value.numerator = scalar * value.numerator;
    return value;
}

Polynomial second_phi_derivative(const Polynomial& polynomial) {
    const Polynomial one_plus_z{
        std::vector<Rational>{Rational{1}, Rational{1}}
    };
    const Polynomial z{
        std::vector<Rational>{Rational{0}, Rational{1}}
    };
    const Polynomial one_plus_three_z{
        std::vector<Rational>{Rational{1}, Rational{3}}
    };
    return Rational{4}
            * z
            * one_plus_z
            * one_plus_z
            * polynomial.derivative().derivative()
        + Rational{2}
            * one_plus_z
            * one_plus_three_z
            * polynomial.derivative();
}

std::vector<Polynomial> cosecant_power_sums(int n, int max_power) {
    std::vector<Polynomial> sums(
        static_cast<std::size_t>(max_power + 1)
    );
    sums[0] = Polynomial{Rational{n}};
    if (max_power == 0) {
        return sums;
    }
    const Integer n_integer{n};
    const Rational n_squared{n_integer * n_integer};
    sums[1] = Polynomial{
        std::vector<Rational>{n_squared, n_squared}
    };
    for (int power = 2; power <= max_power; ++power) {
        const Rational r{2 * power - 2};
        sums[static_cast<std::size_t>(power)] =
            Rational{1} / (r * (r + Rational{1}))
            * (
                n_squared
                    * second_phi_derivative(
                        sums[static_cast<std::size_t>(power - 1)]
                    )
                + r * r * sums[static_cast<std::size_t>(power - 1)]
            );
    }
    return sums;
}

RationalPolynomial first_block_determinant(
    int q_half,
    int target_half,
    int power
) {
    const int n = 2 * q_half + 1;
    const int max_h = q_half + target_half;
    const std::vector<Polynomial> c =
        cosecant_power_sums(n, power);

    std::vector<std::vector<RationalPolynomial>> f(
        static_cast<std::size_t>(power + 1),
        std::vector<RationalPolynomial>(
            static_cast<std::size_t>(max_h + 1)
        )
    );
    for (int a = 0; a <= max_h; ++a) {
        if (a == 0) {
            f[0][0] = RationalPolynomial{
                Polynomial{Rational{n}},
                0
            };
        } else if (a == n) {
            f[0][static_cast<std::size_t>(a)] = RationalPolynomial{
                Polynomial{
                    std::vector<Rational>{Rational{-n}, Rational{n}}
                },
                1
            };
        }
    }

    for (int p = 1; p <= power; ++p) {
        f[static_cast<std::size_t>(p)][0] =
            RationalPolynomial{c[static_cast<std::size_t>(p)], 0};
        if (max_h >= 1) {
            f[static_cast<std::size_t>(p)][1] =
                RationalPolynomial{
                    c[static_cast<std::size_t>(p)]
                    - Rational{2}
                        * c[static_cast<std::size_t>(p - 1)],
                    0
                };
        }
        for (int a = 1; a < max_h; ++a) {
            f[static_cast<std::size_t>(p)]
             [static_cast<std::size_t>(a + 1)] =
                Rational{2}
                    * f[static_cast<std::size_t>(p)]
                       [static_cast<std::size_t>(a)]
                - f[static_cast<std::size_t>(p)]
                   [static_cast<std::size_t>(a - 1)]
                - Rational{4}
                    * f[static_cast<std::size_t>(p - 1)]
                       [static_cast<std::size_t>(a)];
        }
    }

    std::vector<RationalPolynomial> m(
        static_cast<std::size_t>(max_h + 1)
    );
    RationalPolynomial running{
        c[static_cast<std::size_t>(power)],
        0
    };
    m[0] = running;
    for (int h = 1; h <= max_h; ++h) {
        running =
            running
            + Rational{2}
                * f[static_cast<std::size_t>(power)]
                   [static_cast<std::size_t>(h)];
        m[static_cast<std::size_t>(h)] = running;
    }

    RationalPolynomial interval_sum;
    const int lower = std::abs(q_half - target_half);
    for (int h = lower; h <= max_h; ++h) {
        interval_sum =
            interval_sum + m[static_cast<std::size_t>(h)];
    }
    return RationalPolynomial{
               c[static_cast<std::size_t>(power)],
               0
           }
               * interval_sum
        - m[static_cast<std::size_t>(q_half)]
              * m[static_cast<std::size_t>(target_half)];
}

long long parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0;
    const long long parsed = std::stoll(value, &consumed);
    if (consumed != value.size() || parsed <= 0) {
        throw std::invalid_argument(name + " must be positive");
    }
    return parsed;
}

}  // namespace

int main(int argc, char** argv) {
    int max_q_half = 12;
    int max_power = 10;
    try {
        if (argc >= 2) {
            max_q_half = static_cast<int>(
                parse_positive(argv[1], "max_q_half")
            );
        }
        if (argc >= 3) {
            max_power = static_cast<int>(
                parse_positive(argv[2], "max_power")
            );
        }
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_anchored_lobe_all_powers"
                " [max_q_half] [max_power]"
            );
        }
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }

    std::size_t cases = 0;
    for (int power = 1; power <= max_power; ++power) {
        for (int q_half = 1; q_half <= max_q_half; ++q_half) {
            for (int target_half = 0;
                 target_half <= 2 * q_half;
                 ++target_half) {
                const RationalPolynomial determinant =
                    first_block_determinant(
                        q_half,
                        target_half,
                        power
                    );
                ++cases;
                const auto& coefficients =
                    determinant.numerator.coefficients();
                for (std::size_t degree = 0;
                     degree < coefficients.size();
                     ++degree) {
                    if (coefficients[degree] < 0) {
                        std::cout
                            << "SU2_ANCHORED_LOBE_ALL_POWERS"
                            << " counterexample"
                            << " power=" << power
                            << " q_half=" << q_half
                            << " target_half=" << target_half
                            << " degree=" << degree
                            << " coefficient=" << coefficients[degree]
                            << " denominator_power="
                            << determinant.denominator_power
                            << '\n';
                        return EXIT_SUCCESS;
                    }
                }
            }
        }
    }

    std::cout
        << "SU2_ANCHORED_LOBE_ALL_POWERS"
        << " cases=" << cases
        << " max_q_half=" << max_q_half
        << " max_power=" << max_power
        << " negative_coefficients=0"
        << " result=NO_NEGATIVE_Z_COEFFICIENT"
        << '\n';
    return EXIT_SUCCESS;
}
