#include <algorithm>
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

long long value_at(const std::vector<long long>& values, int index) {
    return index < 0
            || index >= static_cast<int>(values.size())
        ? 0
        : values[static_cast<std::size_t>(index)];
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_half_level
            = argc == 2 ? positive_argument(argv[1]) : 100;
        if (argc > 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim_three_box_label_shape "
                "[maximum_half_level]");
        }
        std::uint64_t rows = 0U;
        std::uint64_t decreasing_rows = 0U;
        std::uint64_t increasing_rows = 0U;
        std::uint64_t concave_rows = 0U;
        std::uint64_t convex_rows = 0U;
        std::uint64_t left_minimum_rows = 0U;
        std::uint64_t right_minimum_rows = 0U;
        std::uint64_t interior_minimum_rows = 0U;
        std::uint64_t negative_rows = 0U;
        std::uint64_t band_one_slope_comparisons = 0U;
        std::uint64_t band_one_slope_failures = 0U;
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
                        const int endpoint
                            = left + width + shift - 2;
                        if (endpoint + 1 <= 2 * period) {
                            continue;
                        }
                        std::vector<long long> profile(
                            static_cast<std::size_t>(endpoint + 2),
                            0);
                        for (int first = 0; first < left; ++first) {
                            for (int second = 0;
                                 second < width;
                                 ++second) {
                                const int base = first + second + 1;
                                for (int third = 0;
                                     third < shift;
                                     ++third) {
                                    ++profile[static_cast<std::size_t>(
                                        base + third)];
                                }
                            }
                        }
                        std::vector<long long> margins;
                        for (int label = 1;
                             label <= maximum_label;
                             ++label) {
                            long long reserve = 0;
                            for (int wall = 1; wall <= 3; ++wall) {
                                reserve += value_at(
                                               profile,
                                               wall * period + label)
                                           - value_at(
                                               profile,
                                               (wall + 1) * period
                                                   - label - 1);
                            }
                            margins.push_back(
                                value_at(
                                    profile,
                                    period - label - 1)
                                - reserve);
                            const int excess
                                = endpoint + 1 - 2 * period;
                            if (excess <= period) {
                                for (int offset = 0;
                                     offset < 2 * label + 1;
                                     ++offset) {
                                    const int early
                                        = excess - label + offset;
                                    const int late
                                        = period - label - 1 + offset;
                                    const long long early_slope
                                        = value_at(profile, early + 1)
                                          - value_at(profile, early);
                                    const long long late_slope
                                        = value_at(profile, late + 1)
                                          - value_at(profile, late);
                                    ++band_one_slope_comparisons;
                                    band_one_slope_failures
                                        += early_slope < late_slope
                                        ? 1U
                                        : 0U;
                                }
                            }
                        }
                        ++rows;
                        bool decreasing = true;
                        bool increasing = true;
                        bool concave = true;
                        bool convex = true;
                        for (std::size_t index = 1;
                             index < margins.size();
                             ++index) {
                            decreasing = decreasing
                                && margins[index] <= margins[index - 1U];
                            increasing = increasing
                                && margins[index] >= margins[index - 1U];
                            if (index >= 2U) {
                                const long long second_difference
                                    = margins[index]
                                    - 2 * margins[index - 1U]
                                    + margins[index - 2U];
                                concave = concave
                                    && second_difference <= 0;
                                convex = convex
                                    && second_difference >= 0;
                            }
                        }
                        decreasing_rows += decreasing ? 1U : 0U;
                        increasing_rows += increasing ? 1U : 0U;
                        concave_rows += concave ? 1U : 0U;
                        convex_rows += convex ? 1U : 0U;
                        const auto minimum = std::min_element(
                            margins.begin(),
                            margins.end());
                        if (minimum != margins.begin()
                            && minimum + 1 != margins.end()) {
                            ++interior_minimum_rows;
                        }
                        left_minimum_rows += minimum == margins.begin()
                            ? 1U
                            : 0U;
                        right_minimum_rows += minimum + 1 == margins.end()
                            ? 1U
                            : 0U;
                        if (*minimum < 0) {
                            ++negative_rows;
                        }
                    }
                }
            }
        }
        std::cout << "SU2_AIM_THREE_BOX_LABEL_SHAPE"
                  << " maximum_half_level=" << maximum_half_level
                  << " rows=" << rows
                  << " decreasing_rows=" << decreasing_rows
                  << " increasing_rows=" << increasing_rows
                  << " concave_rows=" << concave_rows
                  << " convex_rows=" << convex_rows
                  << " left_minimum_rows=" << left_minimum_rows
                  << " right_minimum_rows=" << right_minimum_rows
                  << " interior_minimum_rows=" << interior_minimum_rows
                  << " negative_rows=" << negative_rows << '\n';
        std::cout << " band_one_slope_comparisons="
                  << band_one_slope_comparisons
                  << " band_one_slope_failures="
                  << band_one_slope_failures << '\n';
        return negative_rows == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
