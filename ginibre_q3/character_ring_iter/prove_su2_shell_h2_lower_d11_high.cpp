// Exact finite-row certificate for 4d<y<=5d at d=11.

#include <array>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

int positive_modulo(int value, int modulus) {
    return (value % modulus + modulus) % modulus;
}

std::vector<Integer> multiply_by_character(
    const std::vector<Integer>& state,
    int d
) {
    std::vector<Integer> result(state.size());
    for (std::size_t source = 0U; source < state.size(); ++source) {
        if (state[source] == 0) continue;
        const int label = static_cast<int>(source);
        for (int target = std::abs(label - d);
             target <= label + d;
             target += 2) {
            result[static_cast<std::size_t>(target)] += state[source];
        }
    }
    return result;
}

struct Quadratic {
    Integer square = 0;
    Integer linear = 0;
    Integer constant = 0;

    Integer evaluate(int q) const {
        return square * q * q + linear * q + constant;
    }
};

Quadratic current(int q_residue, int y, int rail) {
    constexpr int d = 11;
    std::vector<Integer> profile(static_cast<std::size_t>(y + 4 * d + 1));
    profile[static_cast<std::size_t>(y)] = 1;
    profile = multiply_by_character(profile, d);
    profile = multiply_by_character(profile, d);
    const std::vector<Integer> p2 = profile;
    profile = multiply_by_character(profile, d);
    const std::vector<Integer> p3 = profile;
    profile = multiply_by_character(profile, d);
    const std::vector<Integer> p4 = profile;
    const int rho = d + 1 + positive_modulo(rail - d - 1, 4);
    const int even_residue = positive_modulo(2 * q_residue + d - 2 * rail, 8);
    const int odd_residue = positive_modulo(2 * rail, 8);
    Quadratic result;
    for (int label = 0; label < static_cast<int>(profile.size()); ++label) {
        if (label % 8 == even_residue) {
            const int g = label <= d ? d + label : 2 * label;
            const Integer h = Integer(label) * label
                + Integer(2 * d + 4) * label
                + Integer(3 * d * d + 6 * d + 2);
            const Integer p2_value = p2[static_cast<std::size_t>(label)];
            result.square += (2 * g + 4) * p2_value;
            result.linear += (
                8 * d * d + 2 * d * g + 4 * g + 16 * d + 8 - 2 * h
            ) * p2_value;
            result.constant -= (h + Integer(2 * d * d + 3 * d) * g)
                * p2_value;
            result.linear += 4 * p4[static_cast<std::size_t>(label)];
            result.constant += 2 * p4[static_cast<std::size_t>(label)];
        }
        if (label % 8 == odd_residue && label >= 2 * rho) {
            const Integer p3_value = p3[static_cast<std::size_t>(label)];
            result.square -= 4 * p3_value;
            result.linear += 4 * (label - 2 * d - 2) * p3_value;
            result.constant += (2 * label + 4 * d * d + 4 * d) * p3_value;
        }
    }
    return result;
}

}  // namespace

int main() {
    try {
        std::size_t cases = 0U;
        for (int y = 45; y <= 55; y += 2) {
            const int first_q = y + 34;
            for (int q_residue = 0; q_residue < 4; ++q_residue) {
                const int q0 = first_q + positive_modulo(q_residue - first_q, 4);
                for (int rail = 0; rail < 4; ++rail) {
                    const Quadratic value = current(q_residue, y, rail);
                    const Integer quadratic = 16 * value.square;
                    const Integer linear = 8 * q0 * value.square
                        + 4 * value.linear;
                    const Integer constant = value.evaluate(q0);
                    ++cases;
                    if (quadratic < 0 || linear < 0 || constant < 0) {
                        std::cout
                            << "SU2_SHELL_H2_LOWER_D11_HIGH_FAILURE"
                            << " y=" << y
                            << " q_residue=" << q_residue
                            << " rail=" << rail
                            << " coefficients=(" << quadratic << ','
                            << linear << ',' << constant << ")\n";
                        return EXIT_FAILURE;
                    }
                }
            }
        }
        std::cout
            << "SU2_SHELL_H2_LOWER_D11_HIGH"
            << " cases=" << cases
            << " certified=" << cases
            << " result=PASS_EXACT_FINITE_ROW_CERTIFICATE\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
