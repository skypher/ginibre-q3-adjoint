#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int parse(const char* text) {
    return std::stoi(std::string(text));
}

long long diagonal_count(int left, int right, int total) {
    const int begin = std::max(0, total - right + 1);
    const int end = std::min(left - 1, total);
    return std::max(0, end - begin + 1);
}

long long value_at(const std::vector<long long>& values, int index) {
    return index < 0
            || index >= static_cast<int>(values.size())
        ? 0
        : values[static_cast<std::size_t>(index)];
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5) {
            throw std::invalid_argument(
                "usage: probe_su2_aim_three_box_boundary_sequence "
                "half_level label left padding");
        }
        const int half_level = parse(argv[1]);
        const int label = parse(argv[2]);
        const int left = parse(argv[3]);
        const int padding = parse(argv[4]);
        const int period = 2 * half_level + 2;
        const int shift = left + 2 * padding;
        const int horizon = 4 * period;
        std::vector<long long> base(
            static_cast<std::size_t>(horizon + 1),
            0);
        for (int index = 0; index <= horizon; ++index) {
            base[static_cast<std::size_t>(index)]
                = diagonal_count(left, shift, index - 1);
        }
        long long left_prefix = 0;
        long long right_prefix = 0;
        std::cout << "left=[";
        for (int layer = 0; layer < period - 1; ++layer) {
            long long reserve = 0;
            for (int wall = 1; wall <= 3; ++wall) {
                reserve += value_at(
                               base,
                               wall * period + label - layer)
                           - value_at(
                               base,
                               (wall + 1) * period - label - 1
                                   - layer);
            }
            const long long increment
                = value_at(
                      base,
                      period - label - 1 - layer)
                  - reserve;
            left_prefix += increment;
            if (layer != 0) {
                std::cout << ',';
            }
            std::cout << increment;
        }
        std::cout << "]\nleft_prefix=[";
        left_prefix = 0;
        for (int layer = 0; layer < period - 1; ++layer) {
            long long reserve = 0;
            for (int wall = 1; wall <= 3; ++wall) {
                reserve += value_at(
                               base,
                               wall * period + label - layer)
                           - value_at(
                               base,
                               (wall + 1) * period - label - 1
                                   - layer);
            }
            left_prefix += value_at(
                               base,
                               period - label - 1 - layer)
                           - reserve;
            if (layer != 0) {
                std::cout << ',';
            }
            std::cout << left_prefix;
        }
        std::cout << "]\nright_prefix=[";
        for (int layer = 0; layer < period - 1; ++layer) {
            long long reserve = 0;
            for (int wall = 1; wall <= 3; ++wall) {
                reserve += value_at(
                               base,
                               wall * period + label - layer)
                           - value_at(
                               base,
                               (wall + 1) * period - label - 1
                                   - layer);
            }
            right_prefix += value_at(
                                base,
                                period + label - 1 - layer)
                            - reserve;
            if (layer != 0) {
                std::cout << ',';
            }
            std::cout << right_prefix;
        }
        std::cout << "]\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
