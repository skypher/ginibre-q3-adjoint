#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_rational;

namespace {

struct Interval {
    cpp_rational lower;
    cpp_rational upper;
};

Interval operator+(const Interval& left, const Interval& right) {
    return {left.lower + right.lower, left.upper + right.upper};
}

Interval operator-(const Interval& left, const Interval& right) {
    return {left.lower - right.upper, left.upper - right.lower};
}

Interval operator*(const Interval& left, const Interval& right) {
    const cpp_rational products[]{
        left.lower * right.lower,
        left.lower * right.upper,
        left.upper * right.lower,
        left.upper * right.upper
    };
    return {
        *std::min_element(std::begin(products), std::end(products)),
        *std::max_element(std::begin(products), std::end(products))
    };
}

Interval point(const cpp_rational& value) {
    return {value, value};
}

Interval power(Interval base, int exponent) {
    Interval answer = point(1);
    while (exponent > 0) {
        if ((exponent & 1) != 0) {
            answer = answer * base;
        }
        exponent /= 2;
        if (exponent > 0) {
            base = base * base;
        }
    }
    return answer;
}

cpp_rational polynomial(const cpp_rational& value) {
    return value * value * value + value * value - 2 * value - 1;
}

Interval weight(const Interval& value) {
    return power(point(2) - value, 2)
        * power(point(2) + value, 5);
}

void require_sign_change(
    const Interval& interval,
    const char* name
) {
    const cpp_rational left = polynomial(interval.lower);
    const cpp_rational right = polynomial(interval.upper);
    if (left == 0 || right == 0 || (left < 0) == (right < 0)) {
        throw std::runtime_error(
            std::string("missing sign change for ") + name
        );
    }
}

}  // namespace

int main() {
    try {
        const cpp_rational scale = 1000;
        const Interval x{cpp_rational(1246) / scale,
                         cpp_rational(1247) / scale};
        const Interval y{-cpp_rational(446) / scale,
                         -cpp_rational(445) / scale};
        const Interval z{-cpp_rational(1802) / scale,
                         -cpp_rational(1801) / scale};
        require_sign_change(x, "x");
        require_sign_change(y, "y");
        require_sign_change(z, "z");

        const Interval same_sign
            = power(y - z, 2) * weight(y)
            * (point(2) + y * z);
        const Interval opposite_sign
            = power(x - z, 2) * weight(x)
            * (point(2) + x * z);
        const Interval bracket = same_sign + opposite_sign;
        if (bracket.upper >= 0) {
            throw std::runtime_error(
                "interval enclosure does not prove the negative margin"
            );
        }
        if ((point(2) + x * z).upper >= 0) {
            throw std::runtime_error(
                "interval enclosure does not classify the negative pair"
            );
        }
        if ((point(2) + y * z).lower <= 0) {
            throw std::runtime_error(
                "interval enclosure does not classify the positive pair"
            );
        }
        std::cout
            << "PASS polynomial=t^3+t^2-2t-1"
            << " bracket_lower=" << bracket.lower
            << " bracket_upper=" << bracket.upper << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
