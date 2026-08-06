// Exact probe for the separated-lower reflected crossing weights.
//
// With K=2Q+d+1 and W=ell-2V, the crossing-tail threshold is
// T=K-Q-V=(d+W)/2+1.  This source evaluates d_1,d_2,d_3 directly from
// Lemma 5A8H28P1 and prints their dependence on Q for small exact grids.

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

bool fuses(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= source + label
        && source + target + label <= 2 * level;
}

std::vector<Integer> multiply(
    const std::vector<Integer>& state,
    const std::vector<std::vector<int>>& matrix
) {
    std::vector<Integer> result(state.size());
    for (std::size_t source = 0U; source < state.size(); ++source) {
        if (state[source] == 0) {
            continue;
        }
        for (std::size_t target = 0U; target < state.size(); ++target) {
            if (matrix[source][target] != 0) {
                result[target] += state[source] * matrix[source][target];
            }
        }
    }
    return result;
}

Integer crossing_weight(int q, int d, int w, int power) {
    if ((d + w) % 2 != 0) {
        throw std::runtime_error("incompatible W parity");
    }
    const int level = 2 * q + d + 1;
    const int threshold = (d + w) / 2 + 1;
    std::vector<std::vector<int>> fusion(
        static_cast<std::size_t>(level + 1),
        std::vector<int>(static_cast<std::size_t>(level + 1))
    );
    for (int source = 0; source <= level; ++source) {
        for (int target = 0; target <= level; ++target) {
            fusion[static_cast<std::size_t>(source)]
                  [static_cast<std::size_t>(target)] =
                fuses(level, q, source, target) ? 1 : 0;
        }
    }
    std::vector<Integer> row(static_cast<std::size_t>(level + 1));
    row[0] = 1;
    for (int step = 0; step < power; ++step) {
        row = multiply(row, fusion);
    }
    Integer tail = 0;
    const int half_lower = (level - 1) / 2;
    for (int u = threshold; u <= half_lower; ++u) {
        tail += row[static_cast<std::size_t>(u)]
            + row[static_cast<std::size_t>(level - u)];
    }
    if (level % 2 == 0 && level / 2 >= threshold) {
        tail += row[static_cast<std::size_t>(level / 2)];
    }
    return 2 * tail;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4) {
            throw std::runtime_error("usage: Q D W");
        }
        const int q = std::atoi(argv[1]);
        const int d = std::atoi(argv[2]);
        const int w = std::atoi(argv[3]);
        if (q <= 3 * d || d < 0 || w < 0 || w > q + d - 1) {
            throw std::runtime_error("invalid input");
        }
        const Integer d1 = crossing_weight(q, d, w, 1);
        const Integer d2 = crossing_weight(q, d, w, 2);
        const Integer d3 = crossing_weight(q, d, w, 3);
        const Integer expected_d1 = 2;
        const Integer expected_d2 = 4 * q - d - w
            - std::max(0, w - d);
        const Integer expected_d3 = 4 * Integer(q) * q
            + Integer(4 * d + 8) * q
            - Integer(w) * w
            - Integer(2 * d + 4) * w
            - Integer(3 * d * d + 6 * d + 2);
        if (d1 != expected_d1 || d2 != expected_d2 || d3 != expected_d3) {
            throw std::runtime_error("separated-lower crossing formula mismatch");
        }
        std::cout
            << "SU2_SHELL_LOWER_CROSSING"
            << " Q=" << q
            << " d=" << d
            << " W=" << w
            << " d1=" << d1
            << " d2=" << d2
            << " d3=" << d3
            << " result=PASS_FORMULA_AUDIT"
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
