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

int positive_argument(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0
        || parsed > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(std::string(name) + " must be positive");
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

struct Witness {
    int half_level = -1;
    int label = -1;
    int left = -1;
    int padding = -1;
    int first_width = -1;
    int minimizing_width = -1;
    long long value = 0;
};

std::string render(const Witness& witness) {
    return "half_level=" + std::to_string(witness.half_level)
        + " label=" + std::to_string(witness.label)
        + " left=" + std::to_string(witness.left)
        + " padding=" + std::to_string(witness.padding)
        + " first_width=" + std::to_string(witness.first_width)
        + " minimizing_width="
        + std::to_string(witness.minimizing_width)
        + " value=" + std::to_string(witness.value);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_half_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 100;
        if (argc > 2 || maximum_half_level < 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim_three_box_boundary_shape "
                "[maximum_half_level]");
        }

        std::uint64_t rows = 0U;
        std::uint64_t hard_prefixes = 0U;
        std::uint64_t negative_left = 0U;
        std::uint64_t negative_right = 0U;
        std::uint64_t interior_left_minima = 0U;
        std::uint64_t interior_right_minima = 0U;
        std::uint64_t positive_interior_left_minima = 0U;
        std::uint64_t left_increment_recrossings = 0U;
        std::uint64_t right_increment_recrossings = 0U;
        std::uint64_t left_multiple_recrossing_rows = 0U;
        std::uint64_t right_multiple_recrossing_rows = 0U;
        Witness first_negative_left;
        Witness first_negative_right;
        Witness first_interior_left;
        Witness first_interior_right;
        Witness first_positive_interior_left;
        Witness first_left_recrossing;
        Witness first_right_recrossing;

        for (int half_level = 2;
             half_level <= maximum_half_level;
             ++half_level) {
            const int period = 2 * half_level + 2;
            const int horizon = 4 * period;
            for (int left = 1; left < period; ++left) {
                for (int padding = 0;
                     padding <= half_level;
                     ++padding) {
                    const int shift = left + 2 * padding;
                    std::vector<long long> base(
                        static_cast<std::size_t>(horizon + 1),
                        0);
                    for (int index = 0; index <= horizon; ++index) {
                        base[static_cast<std::size_t>(index)]
                            = diagonal_count(
                                left,
                                shift,
                                index - 1);
                    }
                    for (int label = 1;
                         label <= 2 * (half_level / 2);
                         ++label) {
                        ++rows;
                        std::vector<long long> left_prefix(
                            static_cast<std::size_t>(period),
                            0);
                        std::vector<long long> right_prefix(
                            static_cast<std::size_t>(period),
                            0);
                        long long left_sum = 0;
                        long long right_sum = 0;
                        int previous_left_sign = 0;
                        int previous_right_sign = 0;
                        int left_recrossings = 0;
                        int right_recrossings = 0;
                        for (int layer = 0;
                             layer < period - 1;
                             ++layer) {
                            long long reserve = 0;
                            for (int wall = 1; wall <= 3; ++wall) {
                                reserve += value_at(
                                               base,
                                               wall * period + label
                                                   - layer)
                                           - value_at(
                                               base,
                                               (wall + 1) * period
                                                   - label - 1 - layer);
                            }
                            const long long left_increment
                                = value_at(
                                      base,
                                      period - label - 1 - layer)
                                  - reserve;
                            const long long right_increment
                                = value_at(
                                      base,
                                      period + label - 1 - layer)
                                  - reserve;
                            left_sum += left_increment;
                            right_sum += right_increment;
                            left_prefix[
                                static_cast<std::size_t>(layer + 1)]
                                = left_sum;
                            right_prefix[
                                static_cast<std::size_t>(layer + 1)]
                                = right_sum;
                            const int left_sign
                                = left_increment < 0
                                ? -1
                                : left_increment > 0 ? 1 : 0;
                            const int right_sign
                                = right_increment < 0
                                ? -1
                                : right_increment > 0 ? 1 : 0;
                            if (previous_left_sign < 0
                                && left_sign > 0) {
                                ++left_increment_recrossings;
                                ++left_recrossings;
                                if (
                                    first_left_recrossing.half_level
                                    < 0) {
                                    first_left_recrossing = {
                                        half_level,
                                        label,
                                        left,
                                        padding,
                                        1,
                                        layer + 1,
                                        left_increment};
                                }
                            }
                            if (previous_right_sign < 0
                                && right_sign > 0) {
                                ++right_increment_recrossings;
                                ++right_recrossings;
                                if (
                                    first_right_recrossing.half_level
                                    < 0) {
                                    first_right_recrossing = {
                                        half_level,
                                        label,
                                        left,
                                        padding,
                                        1,
                                        layer + 1,
                                        right_increment};
                                }
                            }
                            if (left_sign != 0) {
                                previous_left_sign = left_sign;
                            }
                            if (right_sign != 0) {
                                previous_right_sign = right_sign;
                            }
                        }
                        if (left_recrossings > 1) {
                            ++left_multiple_recrossing_rows;
                        }
                        if (right_recrossings > 1) {
                            ++right_multiple_recrossing_rows;
                        }

                        const int first_width = std::max(
                            1,
                            2 * period - 2 * left
                                - 2 * padding + 2);
                        if (first_width >= period) {
                            continue;
                        }
                        ++hard_prefixes;
                        long long minimum_left
                            = std::numeric_limits<long long>::max();
                        long long minimum_right
                            = std::numeric_limits<long long>::max();
                        int minimizing_left = -1;
                        int minimizing_right = -1;
                        for (int width = first_width;
                             width < period;
                             ++width) {
                            const long long left_value
                                = left_prefix[
                                    static_cast<std::size_t>(width)];
                            const long long right_value
                                = right_prefix[
                                    static_cast<std::size_t>(width)];
                            if (left_value < minimum_left) {
                                minimum_left = left_value;
                                minimizing_left = width;
                            }
                            if (right_value < minimum_right) {
                                minimum_right = right_value;
                                minimizing_right = width;
                            }
                        }
                        if (minimum_left < 0) {
                            ++negative_left;
                            if (first_negative_left.half_level < 0) {
                                first_negative_left = {
                                    half_level,
                                    label,
                                    left,
                                    padding,
                                    first_width,
                                    minimizing_left,
                                    minimum_left};
                            }
                        }
                        if (minimum_right < 0) {
                            ++negative_right;
                            if (first_negative_right.half_level < 0) {
                                first_negative_right = {
                                    half_level,
                                    label,
                                    left,
                                    padding,
                                    first_width,
                                    minimizing_right,
                                    minimum_right};
                            }
                        }
                        if (minimizing_left != first_width
                            && minimizing_left != period - 1) {
                            ++interior_left_minima;
                            if (first_interior_left.half_level < 0) {
                                first_interior_left = {
                                    half_level,
                                    label,
                                    left,
                                    padding,
                                    first_width,
                                    minimizing_left,
                                    minimum_left};
                            }
                            if (minimum_left > 0) {
                                ++positive_interior_left_minima;
                                if (
                                    first_positive_interior_left.half_level
                                    < 0) {
                                    first_positive_interior_left = {
                                        half_level,
                                        label,
                                        left,
                                        padding,
                                        first_width,
                                        minimizing_left,
                                        minimum_left};
                                }
                            }
                        }
                        if (minimizing_right != first_width
                            && minimizing_right != period - 1) {
                            ++interior_right_minima;
                            if (first_interior_right.half_level < 0) {
                                first_interior_right = {
                                    half_level,
                                    label,
                                    left,
                                    padding,
                                    first_width,
                                    minimizing_right,
                                    minimum_right};
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_AIM_THREE_BOX_BOUNDARY_SHAPE"
            << " maximum_half_level=" << maximum_half_level
            << " rows=" << rows
            << " hard_prefixes=" << hard_prefixes
            << " negative_left=" << negative_left
            << " negative_right=" << negative_right
            << " interior_left_minima=" << interior_left_minima
            << " interior_right_minima=" << interior_right_minima
            << " positive_interior_left_minima="
            << positive_interior_left_minima
            << " left_increment_recrossings="
            << left_increment_recrossings
            << " right_increment_recrossings="
            << right_increment_recrossings
            << " left_multiple_recrossing_rows="
            << left_multiple_recrossing_rows
            << " right_multiple_recrossing_rows="
            << right_multiple_recrossing_rows
            << '\n'
            << "FIRST_NEGATIVE_LEFT "
            << render(first_negative_left) << '\n'
            << "FIRST_NEGATIVE_RIGHT "
            << render(first_negative_right) << '\n'
            << "FIRST_INTERIOR_LEFT "
            << render(first_interior_left) << '\n'
            << "FIRST_INTERIOR_RIGHT "
            << render(first_interior_right) << '\n'
            << "FIRST_POSITIVE_INTERIOR_LEFT "
            << render(first_positive_interior_left) << '\n'
            << "FIRST_LEFT_RECROSSING "
            << render(first_left_recrossing) << '\n'
            << "FIRST_RIGHT_RECROSSING "
            << render(first_right_recrossing) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
