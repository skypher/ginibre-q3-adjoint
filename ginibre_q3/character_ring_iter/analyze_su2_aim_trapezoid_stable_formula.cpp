#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#define main analyze_su2_aim_trapezoid_coefficients_embedded_main
#include "analyze_su2_aim_trapezoid_coefficients.cpp"
#undef main

namespace {

struct StableFailure {
    int half_level = -1;
    int shell = -1;
    int reach = -1;
    int width = -1;
    int height = -1;
    int padding = -1;
    long long value = 0;
};

std::string render_stable_failure(const StableFailure& failure) {
    return "half_level=" + std::to_string(failure.half_level)
        + " shell=" + std::to_string(failure.shell)
        + " reach=" + std::to_string(failure.reach)
        + " width=" + std::to_string(failure.width)
        + " height=" + std::to_string(failure.height)
        + " padding=" + std::to_string(failure.padding)
        + " value=" + std::to_string(failure.value);
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
        const int maximum_half_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 40;
        if (argc > 2 || maximum_half_level < 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim_trapezoid_stable_formula "
                "[maximum_half_level]");
        }

        std::uint64_t shapes = 0U;
        std::uint64_t reserve_differences = 0U;
        std::uint64_t negative_reserve_differences = 0U;
        std::uint64_t payments = 0U;
        std::uint64_t failures = 0U;
        std::uint64_t positive_only_failures = 0U;
        std::uint64_t maximum_reserve_bound_failures = 0U;
        std::uint64_t double_prefix_bound_failures = 0U;
        std::uint64_t pairwise_reserve_failures = 0U;
        std::uint64_t left_reserve_failures = 0U;
        std::uint64_t right_reserve_failures = 0U;
        std::uint64_t rising_arc_checks = 0U;
        std::uint64_t rising_arc_failures = 0U;
        std::uint64_t crossing_arc_checks = 0U;
        std::uint64_t crossing_arc_failures = 0U;
        long long minimum_margin = std::numeric_limits<long long>::max();
        long long minimum_positive_margin
            = std::numeric_limits<long long>::max();
        long long minimum_nontrivial_pair_margin
            = std::numeric_limits<long long>::max();
        StableFailure first_reserve_failure;
        StableFailure first_payment_failure;
        StableFailure first_positive_only_failure;
        StableFailure first_maximum_reserve_bound_failure;
        StableFailure first_double_prefix_bound_failure;
        StableFailure first_pairwise_reserve_failure;
        StableFailure first_left_reserve_failure;
        StableFailure first_right_reserve_failure;
        StableFailure first_rising_arc_failure;
        StableFailure first_crossing_arc_failure;
        StableFailure minimum_positive_witness;
        StableFailure minimum_nontrivial_pair_witness;

        for (int half_level = 2;
             half_level <= maximum_half_level;
             ++half_level) {
            const int period = 2 * half_level + 2;
            const int maximum_reach
                = 2 * (half_level / 2);
            const int horizon = 4 * period;
            for (int width = 1; width < period; ++width) {
                for (int height = 1; height < period; ++height) {
                    for (int padding = 0;
                         padding <= half_level;
                         ++padding) {
                        ++shapes;
                        const int shift = width + 2 * padding;
                        std::vector<long long> crossing(
                            static_cast<std::size_t>(horizon + 1),
                            0);
                        std::vector<long long> crossing_prefix(
                            static_cast<std::size_t>(horizon + 2),
                            0);
                        for (int index = 0;
                             index <= horizon;
                             ++index) {
                            crossing[
                                static_cast<std::size_t>(index)]
                                = lower_pair_count(
                                      width,
                                      height,
                                      index)
                                  - lower_pair_count(
                                      width,
                                      height,
                                      index - shift);
                            crossing_prefix[
                                static_cast<std::size_t>(index + 1)]
                                = crossing_prefix[
                                      static_cast<std::size_t>(index)]
                                  + crossing[
                                        static_cast<std::size_t>(
                                            index)];
                        }
                        for (int shell = 1; shell <= 2; ++shell) {
                            std::vector<long long> reserve(
                                static_cast<std::size_t>(
                                    half_level + 2),
                                0);
                            std::vector<long long> positive_reserve(
                                static_cast<std::size_t>(
                                    half_level + 2),
                                0);
                            for (int label = 0;
                                 label <= half_level;
                                 ++label) {
                                for (int wall = shell;
                                     wall * period + label <= horizon;
                                     ++wall) {
                                    reserve[
                                        static_cast<std::size_t>(label)]
                                        += value_at(
                                               crossing,
                                               wall * period + label)
                                           - value_at(
                                               crossing,
                                               (wall + 1) * period
                                                   - label - 1);
                                    positive_reserve[
                                        static_cast<std::size_t>(label)]
                                        += value_at(
                                            crossing,
                                            wall * period + label);
                                }
                            }
                            std::vector<long long> reserve_prefix(
                                static_cast<std::size_t>(
                                    half_level + 2),
                                0);
                            for (int label = 0;
                                 label <= half_level;
                                 ++label) {
                                reserve_prefix[
                                    static_cast<std::size_t>(label + 1)]
                                    = reserve_prefix[
                                          static_cast<std::size_t>(
                                              label)]
                                      + reserve[
                                            static_cast<std::size_t>(
                                                label)];
                                const long long difference
                                    = reserve[
                                          static_cast<std::size_t>(
                                              label)]
                                      - reserve[
                                            static_cast<std::size_t>(
                                                label + 1)];
                                ++reserve_differences;
                                if (difference >= 0) {
                                    continue;
                                }
                                ++negative_reserve_differences;
                                if (
                                    first_reserve_failure.half_level
                                    < 0) {
                                    first_reserve_failure = {
                                        half_level,
                                        shell,
                                        label,
                                        width,
                                        height,
                                        padding,
                                        difference};
                                }
                            }
                            for (int reach = 1;
                                 reach <= maximum_reach;
                                 ++reach) {
                                const int wall = shell * period;
                                const int center_sum
                                    = width + height + shift - 1;
                                for (int label = 1;
                                     label <= reach;
                                     ++label) {
                                    const long long pair_margin
                                        = value_at(
                                              crossing,
                                              wall - label)
                                          + value_at(
                                              crossing,
                                              wall + label - 1)
                                          - 2
                                                * value_at(
                                                    reserve,
                                                    label);
                                    const long long left_margin
                                        = value_at(
                                              crossing,
                                              wall - label)
                                          - value_at(reserve, label);
                                    const long long right_margin
                                        = value_at(
                                              crossing,
                                              wall + label - 1)
                                          - value_at(reserve, label);
                                    if (shell == 1
                                        && 2 * period < center_sum
                                        && center_sum <= 3 * period
                                        && value_at(
                                               crossing,
                                               period + label)
                                            > value_at(
                                                crossing,
                                                2 * period
                                                    - label - 1)) {
                                        const int offset
                                            = center_sum - 2 * period;
                                        const long long arc_margin
                                            = value_at(
                                                  crossing,
                                                  offset + label + 1)
                                              - value_at(
                                                    crossing,
                                                    offset - label)
                                              - value_at(
                                                    crossing,
                                                    period + label)
                                              + value_at(
                                                    crossing,
                                                    period - label);
                                        ++rising_arc_checks;
                                        if (arc_margin < 0) {
                                            ++rising_arc_failures;
                                            if (
                                                first_rising_arc_failure
                                                    .half_level
                                                < 0) {
                                                first_rising_arc_failure = {
                                                    half_level,
                                                    shell,
                                                    label,
                                                    width,
                                                    height,
                                                    padding,
                                                    arc_margin};
                                            }
                                        }
                                    }
                                    if (shell == 1
                                        && 2 * period < center_sum
                                        && center_sum <= 3 * period
                                        && value_at(
                                               crossing,
                                               period + label)
                                            <= value_at(
                                                crossing,
                                                2 * period
                                                    - label - 1)) {
                                        const int offset
                                            = center_sum - 2 * period;
                                        const long long arc_margin
                                            = value_at(
                                                  crossing,
                                                  offset - label)
                                              - value_at(
                                                    crossing,
                                                    offset - period
                                                        + label + 1)
                                              - value_at(
                                                    crossing,
                                                    offset + label + 1)
                                              + value_at(
                                                    crossing,
                                                    period + label);
                                        ++crossing_arc_checks;
                                        if (arc_margin < 0) {
                                            ++crossing_arc_failures;
                                            if (
                                                first_crossing_arc_failure
                                                    .half_level
                                                < 0) {
                                                first_crossing_arc_failure = {
                                                    half_level,
                                                    shell,
                                                    label,
                                                    width,
                                                    height,
                                                    padding,
                                                    arc_margin};
                                            }
                                        }
                                    }
                                    if (value_at(reserve, label) > 0
                                        && pair_margin
                                            < minimum_nontrivial_pair_margin) {
                                        minimum_nontrivial_pair_margin
                                            = pair_margin;
                                        minimum_nontrivial_pair_witness = {
                                            half_level,
                                            shell,
                                            label,
                                            width,
                                            height,
                                            padding,
                                            pair_margin};
                                    }
                                    if (pair_margin >= 0) {
                                        if (left_margin < 0) {
                                            ++left_reserve_failures;
                                            if (
                                                first_left_reserve_failure
                                                    .half_level
                                                < 0) {
                                                first_left_reserve_failure = {
                                                    half_level,
                                                    shell,
                                                    label,
                                                    width,
                                                    height,
                                                    padding,
                                                    left_margin};
                                            }
                                        }
                                        if (right_margin < 0) {
                                            ++right_reserve_failures;
                                            if (
                                                first_right_reserve_failure
                                                    .half_level
                                                < 0) {
                                                first_right_reserve_failure = {
                                                    half_level,
                                                    shell,
                                                    label,
                                                    width,
                                                    height,
                                                    padding,
                                                    right_margin};
                                            }
                                        }
                                        continue;
                                    }
                                    ++pairwise_reserve_failures;
                                    if (
                                        first_pairwise_reserve_failure
                                            .half_level
                                        < 0) {
                                        first_pairwise_reserve_failure = {
                                            half_level,
                                            shell,
                                            label,
                                            width,
                                            height,
                                            padding,
                                            pair_margin};
                                    }
                                }
                                const long long supply
                                    = crossing_prefix[
                                          static_cast<std::size_t>(
                                              wall + reach)]
                                      - crossing_prefix[
                                            static_cast<std::size_t>(
                                                wall - reach)];
                                long long load = 0;
                                long long positive_load = 0;
                                if (reach == 1) {
                                    load = value_at(reserve, 3)
                                        + value_at(reserve, 4);
                                    positive_load = value_at(
                                                            positive_reserve,
                                                            3)
                                        + value_at(positive_reserve, 4);
                                } else if (reach % 2 == 0) {
                                    load = 2
                                            * (reserve_prefix[
                                                   static_cast<
                                                       std::size_t>(
                                                       reach)]
                                               - reserve_prefix[1])
                                        + value_at(reserve, reach)
                                        + value_at(
                                            reserve,
                                            reach + 1);
                                    for (int label = 1;
                                         label < reach;
                                         ++label) {
                                        positive_load += 2 * value_at(
                                            positive_reserve,
                                            label);
                                    }
                                    positive_load += value_at(
                                                             positive_reserve,
                                                             reach)
                                        + value_at(
                                            positive_reserve,
                                            reach + 1);
                                } else {
                                    load = 2
                                            * (reserve_prefix[
                                                   static_cast<
                                                       std::size_t>(
                                                       reach)]
                                               - reserve_prefix[2])
                                        + value_at(reserve, reach)
                                        + value_at(
                                            reserve,
                                            reach + 1)
                                        + value_at(
                                            reserve,
                                            reach + 2)
                                        + value_at(
                                            reserve,
                                            reach + 3);
                                    for (int label = 2;
                                         label < reach;
                                         ++label) {
                                        positive_load += 2 * value_at(
                                            positive_reserve,
                                            label);
                                    }
                                    for (int label = reach;
                                         label <= reach + 3;
                                         ++label) {
                                        positive_load += value_at(
                                            positive_reserve,
                                            label);
                                    }
                                }
                                const long long margin = supply - load;
                                const long long maximum_reserve_margin
                                    = supply
                                      - 2 * reach
                                            * value_at(reserve, 1);
                                const long long double_prefix_margin
                                    = supply
                                      - 2
                                            * (reserve_prefix[
                                                   static_cast<
                                                       std::size_t>(
                                                       reach + 1)]
                                               - reserve_prefix[1]);
                                const long long positive_margin
                                    = supply - positive_load;
                                ++payments;
                                minimum_margin
                                    = std::min(minimum_margin, margin);
                                if (margin > 0
                                    && margin
                                        < minimum_positive_margin) {
                                    minimum_positive_margin = margin;
                                    minimum_positive_witness = {
                                        half_level,
                                        shell,
                                        reach,
                                        width,
                                        height,
                                        padding,
                                        margin};
                                }
                                if (margin >= 0) {
                                    if (maximum_reserve_margin < 0) {
                                        ++maximum_reserve_bound_failures;
                                        if (
                                            first_maximum_reserve_bound_failure
                                                .half_level
                                            < 0) {
                                            first_maximum_reserve_bound_failure
                                                = {
                                                    half_level,
                                                    shell,
                                                    reach,
                                                    width,
                                                    height,
                                                    padding,
                                                    maximum_reserve_margin};
                                        }
                                    }
                                    if (double_prefix_margin < 0) {
                                        ++double_prefix_bound_failures;
                                        if (
                                            first_double_prefix_bound_failure
                                                .half_level
                                            < 0) {
                                            first_double_prefix_bound_failure
                                                = {
                                                    half_level,
                                                    shell,
                                                    reach,
                                                    width,
                                                    height,
                                                    padding,
                                                    double_prefix_margin};
                                        }
                                    }
                                    if (positive_margin < 0) {
                                        ++positive_only_failures;
                                        if (
                                            first_positive_only_failure
                                                .half_level
                                            < 0) {
                                            first_positive_only_failure = {
                                                half_level,
                                                shell,
                                                reach,
                                                width,
                                                height,
                                                padding,
                                                positive_margin};
                                        }
                                    }
                                    continue;
                                }
                                ++failures;
                                if (
                                    first_payment_failure.half_level
                                    < 0) {
                                    first_payment_failure = {
                                        half_level,
                                        shell,
                                        reach,
                                        width,
                                        height,
                                        padding,
                                        margin};
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_AIM_TRAPEZOID_STABLE_FORMULA"
            << " maximum_half_level=" << maximum_half_level
            << " shapes=" << shapes
            << " reserve_differences=" << reserve_differences
            << " negative_reserve_differences="
            << negative_reserve_differences
            << " payments=" << payments
            << " failures=" << failures
            << " positive_only_failures="
            << positive_only_failures
            << " maximum_reserve_bound_failures="
            << maximum_reserve_bound_failures
            << " double_prefix_bound_failures="
            << double_prefix_bound_failures
            << " pairwise_reserve_failures="
            << pairwise_reserve_failures
            << " left_reserve_failures="
            << left_reserve_failures
            << " right_reserve_failures="
            << right_reserve_failures
            << " rising_arc_checks=" << rising_arc_checks
            << " rising_arc_failures=" << rising_arc_failures
            << " crossing_arc_checks=" << crossing_arc_checks
            << " crossing_arc_failures=" << crossing_arc_failures
            << " minimum_margin=" << minimum_margin
            << " minimum_positive_margin="
            << minimum_positive_margin
            << " minimum_nontrivial_pair_margin="
            << minimum_nontrivial_pair_margin
            << '\n'
            << "FIRST_RESERVE_FAILURE "
            << render_stable_failure(first_reserve_failure)
            << '\n'
            << "FIRST_PAYMENT_FAILURE "
            << render_stable_failure(first_payment_failure)
            << '\n'
            << "FIRST_POSITIVE_ONLY_FAILURE "
            << render_stable_failure(first_positive_only_failure)
            << '\n'
            << "FIRST_MAXIMUM_RESERVE_BOUND_FAILURE "
            << render_stable_failure(
                first_maximum_reserve_bound_failure)
            << '\n'
            << "FIRST_DOUBLE_PREFIX_BOUND_FAILURE "
            << render_stable_failure(
                first_double_prefix_bound_failure)
            << '\n'
            << "FIRST_PAIRWISE_RESERVE_FAILURE "
            << render_stable_failure(
                first_pairwise_reserve_failure)
            << '\n'
            << "FIRST_LEFT_RESERVE_FAILURE "
            << render_stable_failure(first_left_reserve_failure)
            << '\n'
            << "FIRST_RIGHT_RESERVE_FAILURE "
            << render_stable_failure(first_right_reserve_failure)
            << '\n'
            << "FIRST_RISING_ARC_FAILURE "
            << render_stable_failure(first_rising_arc_failure)
            << '\n'
            << "FIRST_CROSSING_ARC_FAILURE "
            << render_stable_failure(first_crossing_arc_failure)
            << '\n'
            << "MINIMUM_POSITIVE "
            << render_stable_failure(minimum_positive_witness)
            << '\n'
            << "MINIMUM_NONTRIVIAL_PAIR "
            << render_stable_failure(
                minimum_nontrivial_pair_witness)
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
