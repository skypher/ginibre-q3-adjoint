#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Matrix = std::vector<std::vector<int>>;

void apply_additive_compound(
    const Matrix& adjacency,
    std::vector<Integer>& state
) {
    const int width = static_cast<int>(adjacency.size());
    const auto index = [width](int first, int second) {
        return static_cast<std::size_t>(first * width + second);
    };
    std::vector<Integer> next(static_cast<std::size_t>(width * width));
    for (int first = 0; first < width; ++first) {
        for (int second = 0; second < width; ++second) {
            const Integer& coefficient = state[index(first, second)];
            if (coefficient == 0) {
                continue;
            }
            for (int output = 0; output < width; ++output) {
                if (adjacency[static_cast<std::size_t>(first)]
                             [static_cast<std::size_t>(output)] != 0) {
                    next[index(output, second)] += coefficient;
                }
                if (adjacency[static_cast<std::size_t>(second)]
                             [static_cast<std::size_t>(output)] != 0) {
                    next[index(first, output)] += coefficient;
                }
            }
        }
    }
    state.swap(next);
}

int replay_r4_obstruction() {
    const Matrix adjacency{
        {0, 1, 0, 0, 0},
        {1, 1, 1, 1, 0},
        {0, 1, 0, 0, 1},
        {0, 1, 0, 1, 0},
        {0, 0, 1, 0, 0},
    };
    constexpr int width = 5;
    const auto index = [](int first, int second) {
        return static_cast<std::size_t>(first * width + second);
    };
    std::vector<Integer> state(static_cast<std::size_t>(width * width));
    state[index(3, 0)] = 1;
    state[index(0, 3)] = -1;
    Integer coefficient = 0;
    for (int pair_power = 1; pair_power <= 4; ++pair_power) {
        apply_additive_compound(adjacency, state);
        apply_additive_compound(adjacency, state);
        coefficient = state[index(4, 0)];
        std::cout
            << "LEAF_LOOP_BOUNDARY_R4"
            << " pair_power=" << pair_power
            << " source=3"
            << " target=4"
            << " coefficient=" << coefficient
            << '\n';
    }
    if (coefficient != -61) {
        throw std::runtime_error("fixed leaf-loop obstruction mismatch");
    }
    std::cout
        << "LEAF_LOOP_BOUNDARY_R4"
        << " result=COUNTEREXAMPLE"
        << " vertices=5"
        << " root=0"
        << " distinguished_loop_vertex=1"
        << " final_coefficient=-61\n";
    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2 || std::string(argv[1]) != "--replay-r4-obstruction") {
            throw std::runtime_error(
                "usage: probe_leaf_loop_boundary --replay-r4-obstruction"
            );
        }
        return replay_r4_obstruction();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
