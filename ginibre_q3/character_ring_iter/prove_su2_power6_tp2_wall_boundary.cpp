#include <boost/multiprecision/cpp_int.hpp>

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;

Integer binomial(int n, int k) {
    if (k < 0 || n < k) {
        return 0;
    }
    if (k > n - k) {
        k = n - k;
    }
    Integer result = 1;
    for (int index = 1; index <= k; ++index) {
        result *= n - k + index;
        result /= index;
    }
    return result;
}

Integer coefficient(int q, int coordinate) {
    constexpr int power = 6;
    if (coordinate < 0 || coordinate > power * q) {
        return 0;
    }
    const int width = 2 * q + 1;
    const int centered_degree = power * q + coordinate;
    Integer result = 0;
    for (
        int index = 0;
        index * width <= centered_degree;
        ++index
    ) {
        Integer term =
            binomial(power, index)
            * binomial(
                centered_degree - index * width + power - 1,
                power - 1
            );
        result += index % 2 == 0 ? term : -term;
    }
    return result;
}

Integer determinant(int q) {
    const int y = 2 * q + 3;
    const Integer v0 = coefficient(q, 0);
    const Integer v1 = coefficient(q, 1);
    const Integer vy = coefficient(q, y);
    const Integer vy_plus = coefficient(q, y + 1);
    const Integer vy_plus_two = coefficient(q, y + 2);
    return
        (v0 - vy) * (v0 - vy_plus_two)
        - (v1 - vy_plus) * (v1 - vy_plus);
}

}  // namespace

int main() {
    constexpr int degree_bound = 10;
    constexpr int extra_differences = 2;
    constexpr int minimum_q = 2;
    std::vector<Integer> values;
    for (
        int q = minimum_q;
        q <= minimum_q + degree_bound + extra_differences;
        ++q
    ) {
        values.push_back(determinant(q));
    }

    std::vector<Integer> newton_coefficients;
    while (!values.empty()) {
        newton_coefficients.push_back(values.front());
        std::vector<Integer> differences;
        for (std::size_t index = 1; index < values.size(); ++index) {
            differences.push_back(values[index] - values[index - 1]);
        }
        values = std::move(differences);
    }

    bool nonnegative = true;
    bool degree_verified = true;
    for (int order = 0;
         order < static_cast<int>(newton_coefficients.size());
         ++order) {
        const Integer& coefficient_value =
            newton_coefficients[static_cast<std::size_t>(order)];
        if (order <= degree_bound && coefficient_value < 0) {
            nonnegative = false;
        }
        if (order > degree_bound && coefficient_value != 0) {
            degree_verified = false;
        }
    }

    std::cout
        << "SU2_POWER6_TP2_WALL_BOUNDARY"
        << " minimum_q=" << minimum_q
        << " newton_coefficients={";
    for (std::size_t index = 0;
         index <= static_cast<std::size_t>(degree_bound);
         ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << newton_coefficients[index];
    }
    std::cout
        << "}"
        << " nonnegative=" << (nonnegative ? 1 : 0)
        << " degree_verified=" << (degree_verified ? 1 : 0)
        << " result="
        << (
            nonnegative && degree_verified
                ? "PASS_EXACT_NEWTON_CERTIFICATE"
                : "FAIL"
        )
        << '\n';
    return
        nonnegative && degree_verified
        ? EXIT_SUCCESS
        : EXIT_FAILURE;
}
