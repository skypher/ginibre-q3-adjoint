#include <array>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using Exponent = std::array<int, 2>;

class Polynomial {
public:
    Polynomial() = default;

    explicit Polynomial(const Rational& constant) {
        if (constant != 0) {
            coefficients_[Exponent{0, 0}] = constant;
        }
    }

    static Polynomial variable(int index) {
        Polynomial result;
        Exponent exponent{0, 0};
        exponent[static_cast<std::size_t>(index)] = 1;
        result.coefficients_[exponent] = 1;
        return result;
    }

    Polynomial& operator+=(const Polynomial& other) {
        for (const auto& [exponent, coefficient] : other.coefficients_) {
            coefficients_[exponent] += coefficient;
            if (coefficients_[exponent] == 0) {
                coefficients_.erase(exponent);
            }
        }
        return *this;
    }

    Polynomial& operator-=(const Polynomial& other) {
        return *this += Rational{-1} * other;
    }

    const std::map<Exponent, Rational>& coefficients() const {
        return coefficients_;
    }

private:
    std::map<Exponent, Rational> coefficients_;

    friend Polynomial operator*(
        const Polynomial& left,
        const Polynomial& right
    );
    friend Polynomial operator*(const Rational& scalar, Polynomial value);
};

Polynomial operator+(Polynomial left, const Polynomial& right) {
    left += right;
    return left;
}

Polynomial operator-(Polynomial left, const Polynomial& right) {
    left -= right;
    return left;
}

Polynomial operator*(
    const Polynomial& left,
    const Polynomial& right
) {
    Polynomial result;
    for (const auto& [left_exponent, left_coefficient] :
         left.coefficients_) {
        for (const auto& [right_exponent, right_coefficient] :
             right.coefficients_) {
            const Exponent exponent{
                left_exponent[0] + right_exponent[0],
                left_exponent[1] + right_exponent[1]
            };
            result.coefficients_[exponent] +=
                left_coefficient * right_coefficient;
            if (result.coefficients_[exponent] == 0) {
                result.coefficients_.erase(exponent);
            }
        }
    }
    return result;
}

Polynomial operator*(const Rational& scalar, Polynomial value) {
    if (scalar == 0) {
        return Polynomial{};
    }
    for (auto& [exponent, coefficient] : value.coefficients_) {
        static_cast<void>(exponent);
        coefficient *= scalar;
    }
    return value;
}

bool report(const std::string& name, const Polynomial& polynomial) {
    std::size_t negative = 0;
    std::cout << "SU2_ANCHORED_LOBE_FIRST_BLOCK"
              << " chamber=" << name
              << " coefficients={";
    bool first = true;
    for (const auto& [exponent, coefficient] :
         polynomial.coefficients()) {
        if (coefficient < 0) {
            ++negative;
        }
        if (!first) {
            std::cout << ',';
        }
        first = false;
        std::cout << '(' << exponent[0] << ',' << exponent[1]
                  << "):" << coefficient;
    }
    std::cout << "}"
              << " negative_coefficients=" << negative << '\n';
    return negative == 0;
}

}  // namespace

int main() {
    const Polynomial x = Polynomial::variable(0);
    const Polynomial y = Polynomial::variable(1);
    const Polynomial one{Rational{1}};

    // Low chamber: S=x, Q=x+y, so 0<=S<=Q.
    const Polynomial low_s = x;
    const Polynomial low_q = x + y;
    const Polynomial low_n = Rational{2} * low_q + one;
    const Polynomial low_z =
        low_n * (
            low_n
            - Rational{1, 3}
                * (Rational{2} * low_s + one)
        );
    const Polynomial low_constant =
        low_z
        - Rational{1, 2}
            * (low_n * low_n - one);

    // High chamber: B=x=S-Q-1 and Q=B+1+y, so Q+1<=S<=2Q.
    const Polynomial high_b = x;
    const Polynomial high_q = x + y + one;
    const Polynomial high_s = high_q + high_b + one;
    const Polynomial high_n = Rational{2} * high_q + one;
    const Polynomial high_t =
        Rational{1, 3}
        * high_b
        * (high_b + one)
        * (high_b + Polynomial{Rational{2}});
    const Polynomial high_common =
        high_q
        * (high_q + one)
        * (
            Rational{2} * high_s + one
            - Rational{1, 3} * high_n
        );
    const Polynomial high_z = high_common - Rational{2} * high_t;
    const Polynomial high_constant =
        high_q
            * (high_q + one)
            * (
                high_n
                    * (
                        Rational{2} * high_s + one
                        - Rational{1, 3} * high_n
                    )
                - Rational{2} * high_s * (high_s + one)
            )
        + Rational{2} * high_n * high_t;

    bool ok = true;
    ok = report("low_z", low_z) && ok;
    ok = report("low_constant", low_constant) && ok;
    ok = report("high_z", high_z) && ok;
    ok = report("high_constant", high_constant) && ok;
    std::cout
        << "SU2_ANCHORED_LOBE_FIRST_BLOCK"
        << " result="
        << (ok ? "PASS_EXACT_POLYNOMIAL_CERTIFICATE" : "FAIL")
        << '\n';
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
