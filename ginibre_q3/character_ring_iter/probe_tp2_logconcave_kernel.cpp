#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>

namespace {

struct Polynomial {
    std::array<std::int64_t, 3U> coefficient{};
};

Polynomial add(const Polynomial& left, const Polynomial& right) {
    Polynomial result;
    for (std::size_t index = 0U; index < result.coefficient.size(); ++index) {
        result.coefficient[index] = left.coefficient[index]
            + right.coefficient[index];
    }
    return result;
}

Polynomial scale(const Polynomial& value, std::int64_t scalar) {
    Polynomial result;
    for (std::size_t index = 0U; index < result.coefficient.size(); ++index) {
        result.coefficient[index] = scalar * value.coefficient[index];
    }
    return result;
}

Polynomial multiply(const Polynomial& left, const Polynomial& right) {
    std::array<std::int64_t, 5U> raw{};
    for (std::size_t first = 0U; first < 3U; ++first) {
        for (std::size_t second = 0U; second < 3U; ++second) {
            raw[first + second] += left.coefficient[first]
                * right.coefficient[second];
        }
    }
    // x^3=5x^2-6x+1 and x^4=19x^2-29x+5.
    Polynomial result;
    result.coefficient[0U] = raw[0U] + raw[3U] + 5 * raw[4U];
    result.coefficient[1U] = raw[1U] - 6 * raw[3U] - 29 * raw[4U];
    result.coefficient[2U] = raw[2U] + 5 * raw[3U] + 19 * raw[4U];
    return result;
}

Polynomial square(const Polynomial& value) {
    return multiply(value, value);
}

Polynomial subtract_scalar(const Polynomial& value, std::int64_t scalar) {
    Polynomial result = value;
    result.coefficient[0U] -= scalar;
    return result;
}

Polynomial kernel(const Polynomial& left_t, const Polynomial& right_t) {
    const Polynomial left_x = subtract_scalar(left_t, 2);
    const Polynomial right_x = subtract_scalar(right_t, 2);
    Polynomial result = multiply(square(left_x), square(right_x));
    result = add(result, scale(square(left_x), -2));
    result = add(result, multiply(left_x, right_x));
    result = add(result, scale(square(right_x), -2));
    result.coefficient[0U] += 6;
    return result;
}

}  // namespace

int main() {
    // The three ordered roots of T^3-5T^2+6T-1 are
    // x, x^2-4x+4, and 1+3x-x^2, with x=4cos^2(pi/7).
    const std::array<Polynomial, 3U> nodes{
        Polynomial{{0, 1, 0}},
        Polynomial{{4, -4, 1}},
        Polynomial{{1, 3, -1}}
    };
    const std::array<std::int64_t, 3U> weights{{1, 2, 3}};
    Polynomial four_energy;
    for (std::size_t left = 0U; left < nodes.size(); ++left) {
        for (std::size_t right = 0U; right < nodes.size(); ++right) {
            const Polynomial difference = add(
                nodes[left], scale(nodes[right], -1)
            );
            const Polynomial term = scale(
                multiply(square(difference), kernel(nodes[left], nodes[right])),
                weights[left] * weights[right]
            );
            four_energy = add(four_energy, term);
        }
    }
    const std::array<std::int64_t, 3U> expected{{78, -130, 32}};
    if (four_energy.coefficient != expected) {
        std::cerr << "TP2_LOGCONCAVE_KERNEL result=IDENTITY_MISMATCH\n";
        return 1;
    }
    std::cout << "TP2_LOGCONCAVE_KERNEL result=COUNTEREXAMPLE"
              << " rank=3 weights=[1,2,3]"
              << " four_energy=78-130x+32x^2"
              << " x=4cos^2(pi/7)"
              << " minimal_polynomial=x^3-5x^2+6x-1\n";
    return 0;
}
