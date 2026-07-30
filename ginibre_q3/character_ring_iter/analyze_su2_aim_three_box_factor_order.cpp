#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int positive_argument(const char* text) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed < 2
        || parsed > std::numeric_limits<int>::max()) {
        throw std::invalid_argument("maximum_half_level must be at least 2");
    }
    return static_cast<int>(parsed);
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

bool increments_nonnegative(
    int period,
    int label,
    int first,
    int second,
    int layers) {
    const int horizon = 4 * period;
    std::vector<long long> base(
        static_cast<std::size_t>(horizon + 1),
        0);
    for (int index = 0; index <= horizon; ++index) {
        base[static_cast<std::size_t>(index)]
            = diagonal_count(first, second, index - 1);
    }
    for (int layer = 0; layer < layers; ++layer) {
        long long margin = value_at(
            base,
            period - label - 1 - layer);
        for (int wall = 1; wall <= 3; ++wall) {
            margin -= value_at(
                          base,
                          wall * period + label - layer)
                - value_at(
                          base,
                          (wall + 1) * period - label - 1 - layer);
        }
        if (margin < 0) {
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_half_level
            = argc == 2 ? positive_argument(argv[1]) : 100;
        if (argc > 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim_three_box_factor_order "
                "[maximum_half_level]");
        }
        std::uint64_t rows = 0U;
        std::uint64_t shortest_successes = 0U;
        std::uint64_t any_successes = 0U;
        for (int half_level = 2;
             half_level <= maximum_half_level;
             ++half_level) {
            const int period = 2 * half_level + 2;
            const int maximum_label = 2 * (half_level / 2);
            for (int left = 1; left < period; ++left) {
                for (int width = 1; width < period; ++width) {
                    for (int padding = 0;
                         padding <= half_level;
                         ++padding) {
                        const int shift = left + 2 * padding;
                        if (left + width + shift - 1
                            <= 2 * period) {
                            continue;
                        }
                        const std::array<int, 3> lengths{
                            left,
                            width,
                            shift};
                        for (int label = 1;
                             label <= maximum_label;
                             ++label) {
                            ++rows;
                            bool any = false;
                            bool shortest = false;
                            const int minimum_length = *std::min_element(
                                lengths.begin(),
                                lengths.end());
                            for (int factor = 0; factor < 3; ++factor) {
                                const int first
                                    = lengths[static_cast<std::size_t>(
                                        (factor + 1) % 3)];
                                const int second
                                    = lengths[static_cast<std::size_t>(
                                        (factor + 2) % 3)];
                                const int layers
                                    = lengths[static_cast<std::size_t>(
                                        factor)];
                                const bool success
                                    = increments_nonnegative(
                                        period,
                                        label,
                                        first,
                                        second,
                                        layers);
                                any = any || success;
                                shortest = shortest
                                    || (success
                                        && layers == minimum_length);
                            }
                            any_successes += any ? 1U : 0U;
                            shortest_successes += shortest ? 1U : 0U;
                        }
                    }
                }
            }
        }
        std::cout << "SU2_AIM_THREE_BOX_FACTOR_ORDER"
                  << " maximum_half_level=" << maximum_half_level
                  << " rows=" << rows
                  << " shortest_successes=" << shortest_successes
                  << " any_successes=" << any_successes << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
