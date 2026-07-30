#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#define main analyze_su2_aim_square_root_diagonals_embedded_main
#include "analyze_su2_aim_square_root_diagonals.cpp"
#undef main

namespace {

using Grid = std::vector<Vector>;

long long rectangle_sum(
    const Grid& prefix,
    int row_begin,
    int row_end,
    int column_begin,
    int column_end) {
    const auto value = [&](int row, int column) {
        if (row < 0 || column < 0) {
            return 0LL;
        }
        return prefix[static_cast<std::size_t>(row)]
                     [static_cast<std::size_t>(column)];
    };
    return value(row_end, column_end)
        - value(row_begin - 1, column_end)
        - value(row_end, column_begin - 1)
        + value(row_begin - 1, column_begin - 1);
}

struct IntervalFailure {
    int level = -1;
    int shell = -1;
    int factor = -1;
    int radius = -1;
    int outer_begin = -1;
    int outer_end = -1;
    int inner_begin = -1;
    int inner_end = -1;
    long long value = 0;
};

std::string render_interval_failure(
    const IntervalFailure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " shell=" + std::to_string(failure.shell)
        + " factor=" + std::to_string(failure.factor)
        + " radius=" + std::to_string(failure.radius)
        + " outer=[" + std::to_string(failure.outer_begin)
        + ',' + std::to_string(failure.outer_end) + ']'
        + " inner=[" + std::to_string(failure.inner_begin)
        + ',' + std::to_string(failure.inner_end) + ']'
        + " value=" + std::to_string(failure.value);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 20;
        const int maximum_shell = argc >= 3
            ? positive_argument(argv[2], "maximum_shell")
            : 3;
        const int maximum_index = argc >= 4
            ? positive_argument(argv[3], "maximum_index")
            : 80;
        if (argc > 4 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim_nested_intervals "
                "[maximum_half_level] [maximum_shell] "
                "[maximum_index]");
        }

        std::uint64_t payments = 0U;
        std::uint64_t nested_interval_pairs = 0U;
        std::uint64_t failures = 0U;
        std::uint64_t separated_period_blocks = 0U;
        std::uint64_t separated_period_block_failures = 0U;
        std::uint64_t width_period_increments = 0U;
        std::uint64_t width_period_increment_failures = 0U;
        std::uint64_t point_period_interactions = 0U;
        std::uint64_t point_period_interaction_failures = 0U;
        std::uint64_t negative_cores = 0U;
        std::uint64_t negative_left_extensions = 0U;
        std::uint64_t negative_right_extensions = 0U;
        std::uint64_t negative_left_one_sided = 0U;
        std::uint64_t negative_right_one_sided = 0U;
        IntervalFailure first_failure;
        IntervalFailure first_period_block_failure;
        IntervalFailure first_width_increment_failure;
        IntervalFailure first_point_period_failure;

        for (int level = 2; level <= maximum_level; ++level) {
            const int period = 2 * level + 2;
            for (int factor = 1;
                 factor <= level / 2;
                 ++factor) {
                const Matrix transform
                    = reserve_transform(level, factor);
                for (int radius = 0; radius <= level; ++radius) {
                    for (int shell = 1;
                         shell <= maximum_shell;
                         ++shell) {
                        ++payments;
                        Vector eta_prefix(
                            static_cast<std::size_t>(
                                2 * maximum_index + 1),
                            0);
                        long long cumulative = 0;
                        for (int index = 0;
                             index <= 2 * maximum_index;
                             ++index) {
                            cumulative += payment_coefficient(
                                index,
                                level,
                                shell,
                                factor,
                                radius,
                                transform[
                                    static_cast<std::size_t>(
                                        radius)]);
                            eta_prefix[
                                static_cast<std::size_t>(index)]
                                = cumulative;
                        }

                        Grid kernel(
                            static_cast<std::size_t>(
                                maximum_index + 1),
                            Vector(
                                static_cast<std::size_t>(
                                    maximum_index + 1),
                                0));
                        for (int left = 0;
                             left <= maximum_index;
                             ++left) {
                            for (int right = 0;
                                 right <= maximum_index;
                                 ++right) {
                                const int lower
                                    = std::abs(left - right);
                                const int upper = left + right;
                                kernel[
                                    static_cast<std::size_t>(left)]
                                    [static_cast<std::size_t>(right)]
                                    = eta_prefix[
                                          static_cast<std::size_t>(
                                              upper)]
                                      - (lower == 0
                                             ? 0
                                             : eta_prefix[
                                                   static_cast<
                                                       std::size_t>(
                                                       lower - 1)]);
                            }
                        }
                        Grid prefix = kernel;
                        for (int row = 0;
                             row <= maximum_index;
                             ++row) {
                            for (int column = 0;
                                 column <= maximum_index;
                                 ++column) {
                                if (row > 0) {
                                    prefix[
                                        static_cast<std::size_t>(row)]
                                        [static_cast<std::size_t>(
                                            column)]
                                        += prefix[
                                            static_cast<std::size_t>(
                                                row - 1)]
                                            [static_cast<std::size_t>(
                                                column)];
                                }
                                if (column > 0) {
                                    prefix[
                                        static_cast<std::size_t>(row)]
                                        [static_cast<std::size_t>(
                                            column)]
                                        += prefix[
                                            static_cast<std::size_t>(
                                                row)]
                                            [static_cast<std::size_t>(
                                                column - 1)];
                                }
                                if (row > 0 && column > 0) {
                                    prefix[
                                        static_cast<std::size_t>(row)]
                                        [static_cast<std::size_t>(
                                            column)]
                                        -= prefix[
                                            static_cast<std::size_t>(
                                                row - 1)]
                                            [static_cast<std::size_t>(
                                                column - 1)];
                                }
                            }
                        }

                        for (int point = 0;
                             point <= maximum_index;
                             ++point) {
                            for (int block_begin = 0;
                                 block_begin + period - 1
                                     <= maximum_index;
                                 ++block_begin) {
                                const long long interaction
                                    = rectangle_sum(
                                        prefix,
                                        point,
                                        point,
                                        block_begin,
                                        block_begin + period - 1);
                                ++point_period_interactions;
                                if (interaction >= 0) {
                                    continue;
                                }
                                ++point_period_interaction_failures;
                                if (
                                    first_point_period_failure.level
                                    < 0) {
                                    first_point_period_failure = {
                                        level,
                                        shell,
                                        factor,
                                        radius,
                                        point,
                                        point,
                                        block_begin,
                                        block_begin + period - 1,
                                        interaction};
                                }
                            }
                        }

                        for (int inner_begin = 0;
                             inner_begin <= maximum_index;
                             ++inner_begin) {
                            for (int inner_end = inner_begin;
                                 inner_end <= maximum_index;
                                 ++inner_end) {
                                long long best_left = 0;
                                int outer_begin = inner_begin;
                                for (int candidate = inner_begin - 1;
                                     candidate >= 0;
                                     --candidate) {
                                    const long long extension
                                        = rectangle_sum(
                                            prefix,
                                            candidate,
                                            inner_begin - 1,
                                            inner_begin,
                                            inner_end);
                                    if (extension < best_left) {
                                        best_left = extension;
                                        outer_begin = candidate;
                                    }
                                }
                                long long best_right = 0;
                                int outer_end = inner_end;
                                for (int candidate = inner_end + 1;
                                     candidate <= maximum_index;
                                     ++candidate) {
                                    const long long extension
                                        = rectangle_sum(
                                            prefix,
                                            inner_end + 1,
                                            candidate,
                                            inner_begin,
                                            inner_end);
                                    if (extension < best_right) {
                                        best_right = extension;
                                        outer_end = candidate;
                                    }
                                }
                                const long long core
                                    = rectangle_sum(
                                        prefix,
                                        inner_begin,
                                        inner_end,
                                        inner_begin,
                                        inner_end);
                                const long long value
                                    = core + best_left + best_right;
                                if (core < 0) {
                                    ++negative_cores;
                                }
                                if (best_left < 0) {
                                    ++negative_left_extensions;
                                }
                                if (best_right < 0) {
                                    ++negative_right_extensions;
                                }
                                if (core + best_left < 0) {
                                    ++negative_left_one_sided;
                                }
                                if (core + best_right < 0) {
                                    ++negative_right_one_sided;
                                }
                                ++nested_interval_pairs;
                                if (value < 0) {
                                    ++failures;
                                    if (first_failure.level < 0) {
                                        first_failure = {
                                            level,
                                            shell,
                                            factor,
                                            radius,
                                            outer_begin,
                                            outer_end,
                                            inner_begin,
                                            inner_end,
                                            value};
                                    }
                                }
                                for (int block_begin = 0;
                                     block_begin + period
                                         <= inner_begin;
                                     ++block_begin) {
                                    const long long block_value
                                        = rectangle_sum(
                                            prefix,
                                            block_begin,
                                            block_begin + period - 1,
                                            inner_begin,
                                            inner_end);
                                    ++separated_period_blocks;
                                    if (block_value >= 0) {
                                        continue;
                                    }
                                    ++separated_period_block_failures;
                                    if (
                                        first_period_block_failure.level
                                        < 0) {
                                        first_period_block_failure = {
                                            level,
                                            shell,
                                            factor,
                                            radius,
                                            block_begin,
                                            block_begin + period - 1,
                                            inner_begin,
                                            inner_end,
                                            block_value};
                                    }
                                }
                                for (int block_end
                                         = inner_end + period;
                                     block_end <= maximum_index;
                                     ++block_end) {
                                    const long long block_value
                                        = rectangle_sum(
                                            prefix,
                                            block_end - period + 1,
                                            block_end,
                                            inner_begin,
                                            inner_end);
                                    ++separated_period_blocks;
                                    if (block_value >= 0) {
                                        continue;
                                    }
                                    ++separated_period_block_failures;
                                    if (
                                        first_period_block_failure.level
                                        < 0) {
                                        first_period_block_failure = {
                                            level,
                                            shell,
                                            factor,
                                            radius,
                                            block_end - period + 1,
                                            block_end,
                                            inner_begin,
                                            inner_end,
                                            block_value};
                                    }
                                }
                                if (
                                    inner_end - inner_begin + 1
                                    > period) {
                                    const int reduced_inner_end
                                        = inner_end - period;
                                    long long best_width_left = 0;
                                    int width_outer_begin
                                        = inner_begin;
                                    for (int candidate
                                             = inner_begin - 1;
                                         candidate >= 0;
                                         --candidate) {
                                        const long long extension
                                            = rectangle_sum(
                                                prefix,
                                                candidate,
                                                inner_begin - 1,
                                                reduced_inner_end + 1,
                                                inner_end);
                                        if (
                                            extension
                                            < best_width_left) {
                                            best_width_left
                                                = extension;
                                            width_outer_begin
                                                = candidate;
                                        }
                                    }
                                    for (int right_gap = 0;
                                         inner_end + right_gap
                                             <= maximum_index;
                                         ++right_gap) {
                                        const int reduced_outer_end
                                            = reduced_inner_end
                                              + right_gap;
                                        const int extended_outer_end
                                            = inner_end + right_gap;
                                        const long long increment
                                            = best_width_left
                                              + rectangle_sum(
                                                    prefix,
                                                    inner_begin,
                                                    reduced_outer_end,
                                                    reduced_inner_end
                                                        + 1,
                                                    inner_end)
                                              + rectangle_sum(
                                                    prefix,
                                                    reduced_outer_end
                                                        + 1,
                                                    extended_outer_end,
                                                    inner_begin,
                                                    inner_end);
                                        ++width_period_increments;
                                        if (increment >= 0) {
                                            continue;
                                        }
                                        ++width_period_increment_failures;
                                        if (
                                            first_width_increment_failure
                                                .level
                                            < 0) {
                                            first_width_increment_failure
                                                = {
                                                    level,
                                                    shell,
                                                    factor,
                                                    radius,
                                                    width_outer_begin,
                                                    extended_outer_end,
                                                    inner_begin,
                                                    inner_end,
                                                    increment};
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_AIM_NESTED_INTERVALS"
            << " maximum_level=" << 2 * maximum_level
            << " maximum_shell=" << maximum_shell
            << " maximum_index=" << maximum_index
            << " payments=" << payments
            << " optimized_nested_interval_pairs="
            << nested_interval_pairs
            << " failures=" << failures
            << " separated_period_blocks="
            << separated_period_blocks
            << " separated_period_block_failures="
            << separated_period_block_failures
            << " width_period_increments="
            << width_period_increments
            << " width_period_increment_failures="
            << width_period_increment_failures
            << " point_period_interactions="
            << point_period_interactions
            << " point_period_interaction_failures="
            << point_period_interaction_failures
            << " negative_cores=" << negative_cores
            << " negative_left_extensions="
            << negative_left_extensions
            << " negative_right_extensions="
            << negative_right_extensions
            << " negative_left_one_sided="
            << negative_left_one_sided
            << " negative_right_one_sided="
            << negative_right_one_sided
            << '\n'
            << "FIRST_FAILURE "
            << render_interval_failure(first_failure)
            << '\n'
            << "FIRST_SEPARATED_PERIOD_BLOCK_FAILURE "
            << render_interval_failure(
                first_period_block_failure)
            << '\n'
            << "FIRST_WIDTH_PERIOD_INCREMENT_FAILURE "
            << render_interval_failure(
                first_width_increment_failure)
            << '\n'
            << "FIRST_POINT_PERIOD_INTERACTION_FAILURE "
            << render_interval_failure(
                first_point_period_failure)
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
