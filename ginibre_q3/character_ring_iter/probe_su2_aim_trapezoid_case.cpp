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
        if (argc != 8) {
            throw std::invalid_argument(
                "usage: probe_su2_aim_trapezoid_case "
                "half_level shell factor radius width height shift");
        }
        const int half_level
            = positive_argument(argv[1], "half_level");
        const int shell = positive_argument(argv[2], "shell");
        const int factor = positive_argument(argv[3], "factor");
        const int radius = std::stoi(argv[4]);
        const int width = positive_argument(argv[5], "width");
        const int height = positive_argument(argv[6], "height");
        const int shift = positive_argument(argv[7], "shift");
        const int period = 2 * half_level + 2;
        const int reach = 2 * factor - radius;
        const int horizon = 4 * period;
        if (half_level < 2 || shell > 2 || factor > half_level / 2
            || radius < 0 || reach <= 0 || width >= period
            || height >= period || shift < width
            || shift >= width + period
            || (shift - width) % 2 != 0) {
            throw std::invalid_argument("parameters out of range");
        }

        std::vector<long long> tau(
            static_cast<std::size_t>(horizon + 1),
            0);
        for (int left = 0; left < width; ++left) {
            for (int right = 0; right < height; ++right) {
                ++tau[static_cast<std::size_t>(left + right)];
            }
        }
        std::cout << "tau=[";
        for (int index = 0; index < width + height - 1; ++index) {
            if (index != 0) {
                std::cout << ',';
            }
            std::cout << tau[static_cast<std::size_t>(index)];
        }
        std::cout << "]\n";

        std::vector<long long> deltas(
            static_cast<std::size_t>(half_level + 1),
            0);
        for (int label = 0; label <= half_level; ++label) {
            for (int index = 0; index <= horizon; ++index) {
                const long long value = selected_slope(
                    index,
                    half_level,
                    shell,
                    label);
                if (index >= shift) {
                    deltas[static_cast<std::size_t>(label)]
                        += value
                           * tau[static_cast<std::size_t>(
                               index - shift)];
                }
                deltas[static_cast<std::size_t>(label)]
                    -= value * tau[static_cast<std::size_t>(index)];
            }
        }
        std::vector<long long> demands(
            static_cast<std::size_t>(half_level + 2),
            0);
        for (int label = half_level; label >= 0; --label) {
            demands[static_cast<std::size_t>(label)]
                = demands[static_cast<std::size_t>(label + 1)]
                  + deltas[static_cast<std::size_t>(label)];
        }
        std::cout << "fold=[";
        for (int label = 0; label <= half_level; ++label) {
            if (label != 0) {
                std::cout << ',';
            }
            const long long initial
                = lower_pair_count(width, height, label)
                  - lower_pair_count(
                      width,
                      height,
                      label - shift);
            const long long reflected
                = lower_pair_count(
                      width,
                      height,
                      period - label - 1)
                  - lower_pair_count(
                      width,
                      height,
                      period - label - 1 - shift);
            std::cout << demands[static_cast<std::size_t>(label)]
                         + initial - reflected;
        }
        std::cout << "] reserve=[";
        for (int label = 0; label <= half_level; ++label) {
            if (label != 0) {
                std::cout << ',';
            }
            std::cout << demands[static_cast<std::size_t>(label)];
        }
        std::cout << "]\n";

        long long supply_sum = 0;
        std::cout << "supplies=[";
        for (int offset = -reach; offset < reach; ++offset) {
            if (offset != -reach) {
                std::cout << ',';
            }
            const int threshold = shell * period + offset;
            const long long supply
                = lower_pair_count(width, height, threshold)
                  - lower_pair_count(
                      width,
                      height,
                      threshold - shift);
            supply_sum += supply;
            std::cout << supply;
        }
        std::cout << "] sum=" << supply_sum << '\n';

        long long demand_sum = 0;
        std::cout << "tokens=[";
        bool first = true;
        int token_count = 0;
        for (int label = 0; token_count < 2 * reach; ++label) {
            const long long weight
                = stable_weight(factor, radius, label);
            for (long long copy = 0; copy < weight; ++copy) {
                if (!first) {
                    std::cout << ',';
                }
                first = false;
                ++token_count;
                const long long demand = label <= half_level
                    ? demands[static_cast<std::size_t>(label)]
                    : 0;
                demand_sum += demand;
                std::cout << label << ':' << demand;
            }
        }
        std::cout << "] sum=" << demand_sum
                  << " margin=" << supply_sum - demand_sum << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
