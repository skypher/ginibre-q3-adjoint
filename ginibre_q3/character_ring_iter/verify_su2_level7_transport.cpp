#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_rational;

namespace {

struct Interval {
    cpp_rational lower;
    cpp_rational upper;
};

cpp_rational polynomial(const cpp_rational& value) {
    return value * value * value - 3 * value + 1;
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

cpp_rational square(const cpp_rational& value) {
    return value * value;
}

}  // namespace

int main() {
    try {
        const cpp_rational scale = 1000;
        const Interval a{cpp_rational(1532) / scale,
                         cpp_rational(1533) / scale};
        const Interval b{cpp_rational(347) / scale,
                         cpp_rational(348) / scale};
        const Interval d{-cpp_rational(1880) / scale,
                         -cpp_rational(1879) / scale};
        require_sign_change(a, "a");
        require_sign_change(b, "b");
        require_sign_change(d, "d");

        const cpp_rational negative_kernel_upper
            = a.upper * (-d.lower) - 2;
        const cpp_rational bd_kernel_lower
            = 2 + b.upper * d.lower;
        const cpp_rational ab_kernel_lower
            = 2 + a.lower * b.lower;
        if (2 + a.upper * d.lower >= 0
            || bd_kernel_lower <= negative_kernel_upper
            || ab_kernel_lower <= negative_kernel_upper) {
            throw std::runtime_error("kernel classification failed");
        }

        const cpp_rational m_ge_n_ratio_lower
            = square(b.lower - d.upper)
            * bd_kernel_lower
            * (4 - square(b.upper))
            / (
                square(a.upper - d.lower)
                * negative_kernel_upper
                * (4 - square(a.lower))
            );
        const cpp_rational n_ge_m_ratio_lower
            = square(a.lower - b.upper)
            * ab_kernel_lower
            * (4 - square(b.upper))
            / (
                square(a.upper - d.lower)
                * negative_kernel_upper
                * (4 - square(d.upper))
            );
        if (m_ge_n_ratio_lower <= cpp_rational(3) / 2) {
            throw std::runtime_error("m>=n ratio bound failed");
        }
        if (n_ge_m_ratio_lower <= cpp_rational(14) / 5) {
            throw std::runtime_error("n>=m ratio bound failed");
        }
        std::cout
            << "PASS polynomial=t^3-3t+1"
            << " m_ge_n_ratio_lower=" << m_ge_n_ratio_lower
            << " n_ge_m_ratio_lower=" << n_ge_m_ratio_lower
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
