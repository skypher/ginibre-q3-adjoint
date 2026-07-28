#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Vector = std::vector<Integer>;

Integer binomial_integer(int top, int bottom) {
    if (bottom < 0 || top < bottom) {
        return 0;
    }
    Integer result = 1;
    for (int index = 1; index <= bottom; ++index) {
        result *= top - bottom + index;
        result /= index;
    }
    return result;
}

Integer ordinary_multiplicity(int power, int q, int y) {
    if (y < 0 || y > power * q) {
        return 0;
    }
    Integer result = 0;
    for (int image = 0; image <= power; ++image) {
        const Integer term =
            binomial_integer(power, image)
            * binomial_integer(
                (power - 2 * image) * q - y
                    + power - 2 - image,
                power - 2
            );
        if (image % 2 == 0) {
            result += term;
        } else {
            result -= term;
        }
    }
    return result;
}

Integer folded_multiplicity(
    int power,
    int K,
    int Q,
    int y
) {
    const int period = 2 * K + 2;
    Integer result = 0;
    for (int image = 0; image <= power; ++image) {
        const int direct = image * period + y;
        const int reflected = (image + 1) * period - 1 - y;
        result += ordinary_multiplicity(power, Q, direct);
        result -= ordinary_multiplicity(power, Q, reflected);
        if (
            direct > power * Q
            && reflected > power * Q
        ) {
            break;
        }
    }
    return result;
}

bool adjacent(int K, int Q, int source, int target) {
    return std::abs(source - Q) <= target
        && target <= std::min(
            source + Q,
            2 * K - source - Q
        );
}

Vector multiply(int K, int Q, const Vector& state) {
    Vector next(static_cast<std::size_t>(K + 1));
    for (int source = 0; source <= K; ++source) {
        const Integer& coefficient =
            state[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        for (int target = 0; target <= K; ++target) {
            if (adjacent(K, Q, source, target)) {
                next[static_cast<std::size_t>(target)] += coefficient;
            }
        }
    }
    return next;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error("usage: MAXIMUM_Q");
        }
        const int maximum_q = std::stoi(argv[1]);
        if (maximum_q < 1) {
            throw std::runtime_error("MAXIMUM_Q must be positive");
        }

        std::uint64_t parameters = 0U;
        std::uint64_t coordinates = 0U;
        for (int Q = 1; Q <= maximum_q; ++Q) {
            for (int h = 0; h <= 8 * Q; ++h) {
                const int K = 2 * Q + 1 + h;
                Vector state(static_cast<std::size_t>(K + 1));
                state[0] = 1;
                for (int power = 1; power <= 9; ++power) {
                    state = multiply(K, Q, state);
                    if (power != 8 && power != 9) {
                        continue;
                    }
                    for (int y = 0; y <= K; ++y) {
                        const Integer expected =
                            folded_multiplicity(power, K, Q, y);
                        const Integer& actual =
                            state[static_cast<std::size_t>(y)];
                        if (actual != expected) {
                            std::cerr
                                << "SU2_K4_FORMULA_MISMATCH"
                                << " K=" << K
                                << " Q=" << Q
                                << " h=" << h
                                << " power=" << power
                                << " y=" << y
                                << " recurrence=" << actual
                                << " formula=" << expected << '\n';
                            throw std::runtime_error(
                                "finite recurrence disagrees with formula"
                            );
                        }
                        ++coordinates;
                    }
                }
                ++parameters;
            }
        }
        std::cout
            << "SU2_K4_INTERMEDIATE_FORMULA_AUDIT"
            << " maximum_Q=" << maximum_q
            << " parameters=" << parameters
            << " coordinates=" << coordinates
            << " result=PASS_EXACT\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_K4_INTERMEDIATE_FORMULA_AUDIT FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
