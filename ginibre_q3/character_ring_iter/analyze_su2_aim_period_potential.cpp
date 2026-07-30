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

struct PotentialFailure {
    int level = -1;
    int shell = -1;
    int factor = -1;
    int radius = -1;
    int index = -1;
    long long value = 0;
    Vector potential;
};

std::string render_potential_failure(
    const PotentialFailure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " shell=" + std::to_string(failure.shell)
        + " factor=" + std::to_string(failure.factor)
        + " radius=" + std::to_string(failure.radius)
        + " index=" + std::to_string(failure.index)
        + " value=" + std::to_string(failure.value)
        + " potential=" + render(failure.potential);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 100;
        const int maximum_shell = argc >= 3
            ? positive_argument(argv[2], "maximum_shell")
            : 6;
        if (argc > 3 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim_period_potential "
                "[maximum_half_level] [maximum_shell]");
        }

        std::uint64_t payments = 0U;
        std::uint64_t potential_checks = 0U;
        std::uint64_t potential_failures = 0U;
        std::uint64_t period_monotonicity_checks = 0U;
        std::uint64_t period_monotonicity_failures = 0U;
        std::uint64_t initial_split_checks = 0U;
        std::uint64_t initial_split_failures = 0U;
        int maximum_last_period_defect = -1;
        PotentialFailure first_potential;
        PotentialFailure first_period_monotonicity;
        PotentialFailure first_initial_split;

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
                        const int horizon
                            = (shell + 5) * period;
                        Vector potential(
                            static_cast<std::size_t>(horizon + 1),
                            0);
                        long long cumulative = 0;
                        for (int index = 0;
                             index <= horizon;
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
                            potential[
                                static_cast<std::size_t>(index)]
                                = cumulative;
                            ++potential_checks;
                            if (cumulative >= 0) {
                                continue;
                            }
                            ++potential_failures;
                            if (first_potential.level < 0) {
                                first_potential = {
                                    level,
                                    shell,
                                    factor,
                                    radius,
                                    index,
                                    cumulative,
                                    potential};
                            }
                        }
                        for (int index = 0;
                             index + period <= horizon;
                             ++index) {
                            const long long margin
                                = potential[
                                      static_cast<std::size_t>(
                                          index + period)]
                                  - potential[
                                      static_cast<std::size_t>(
                                          index)];
                            ++period_monotonicity_checks;
                            if (margin < 0) {
                                ++period_monotonicity_failures;
                                if (
                                    first_period_monotonicity.level
                                    < 0) {
                                    first_period_monotonicity = {
                                        level,
                                        shell,
                                        factor,
                                        radius,
                                        index,
                                        margin,
                                        potential};
                                }
                            }
                            if (margin != 0) {
                                maximum_last_period_defect = std::max(
                                    maximum_last_period_defect,
                                    index - shell * period);
                            }
                        }
                        Vector potential_prefix(
                            static_cast<std::size_t>(period + 1),
                            0);
                        long long first_period_sum = 0;
                        for (int index = 0;
                             index < period;
                             ++index) {
                            first_period_sum += potential[
                                static_cast<std::size_t>(index)];
                            potential_prefix[
                                static_cast<std::size_t>(index + 1)]
                                = first_period_sum;
                        }
                        for (int split = 0;
                             split < period;
                             ++split) {
                            const long long margin
                                = first_period_sum
                                  - potential_prefix[
                                      static_cast<std::size_t>(
                                          split)]
                                  - potential_prefix[
                                      static_cast<std::size_t>(
                                          period - split - 1)];
                            ++initial_split_checks;
                            if (margin >= 0) {
                                continue;
                            }
                            ++initial_split_failures;
                            if (first_initial_split.level < 0) {
                                first_initial_split = {
                                    level,
                                    shell,
                                    factor,
                                    radius,
                                    split,
                                    margin,
                                    potential};
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_AIM_PERIOD_POTENTIAL"
            << " maximum_level=" << 2 * maximum_level
            << " maximum_shell=" << maximum_shell
            << " payments=" << payments
            << " potential_checks=" << potential_checks
            << " potential_failures=" << potential_failures
            << " period_monotonicity_checks="
            << period_monotonicity_checks
            << " period_monotonicity_failures="
            << period_monotonicity_failures
            << " initial_split_checks=" << initial_split_checks
            << " initial_split_failures=" << initial_split_failures
            << " maximum_last_period_defect="
            << maximum_last_period_defect
            << '\n'
            << "FIRST_POTENTIAL_FAILURE "
            << render_potential_failure(first_potential)
            << '\n'
            << "FIRST_PERIOD_MONOTONICITY_FAILURE "
            << render_potential_failure(
                first_period_monotonicity)
            << '\n'
            << "FIRST_INITIAL_SPLIT_FAILURE "
            << render_potential_failure(first_initial_split)
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
