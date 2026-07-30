#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#define main analyze_su2_aim_trapezoid_coefficients_embedded_main
#include "analyze_su2_aim_trapezoid_coefficients.cpp"
#undef main

int main(int argc, char** argv) {
    try {
        if (argc != 5) {
            throw std::invalid_argument(
                "usage: probe_su2_aim_stable_potential "
                "half_level shell factor radius");
        }
        const int half_level
            = positive_argument(argv[1], "half_level");
        const int shell = positive_argument(argv[2], "shell");
        const int factor = positive_argument(argv[3], "factor");
        const int radius = std::stoi(argv[4]);
        if (half_level < 2 || shell > 2 || factor > half_level / 2
            || radius < 0 || radius >= 2 * factor) {
            throw std::invalid_argument("parameters out of range");
        }
        const int period = 2 * half_level + 2;
        std::vector<long long> potential(
            static_cast<std::size_t>(4 * period),
            0);
        for (int index = 0; index < 4 * period; ++index) {
            long long value = boundary_potential(
                index,
                half_level,
                shell,
                factor,
                radius);
            for (int label = 0; label <= half_level; ++label) {
                value -= stable_coefficient(factor, radius, label)
                    * selected_slope(
                        index,
                        half_level,
                        shell,
                        label);
            }
            potential[static_cast<std::size_t>(index)] = value;
        }
        std::cout << "K=" << half_level
                  << " L=" << period
                  << " n=" << shell
                  << " q=" << factor
                  << " r=" << radius << '\n';
        for (int block = 0; block < 4; ++block) {
            std::cout << "block=" << block << " values=[";
            for (int residue = 0; residue < period; ++residue) {
                if (residue != 0) {
                    std::cout << ',';
                }
                std::cout << potential[static_cast<std::size_t>(
                    block * period + residue)];
            }
            std::cout << "]\n";
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
