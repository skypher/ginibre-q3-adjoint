#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Vector = std::vector<cpp_int>;

struct Failure {
    int level = -1;
    int factor = -1;
    int half_power = -1;
    int left_radius = -1;
    int right_radius = -1;
    int cutoff = -1;
    cpp_int value = 0;
};

int parse_positive(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0LL
        || parsed > static_cast<long long>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    }
    return static_cast<int>(parsed);
}

Vector fusion_step(const Vector& current, int level, int factor) {
    Vector next(static_cast<std::size_t>(level + 1), 0);
    for (int source = 0; source <= level; ++source) {
        const cpp_int& coefficient = current[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        const int lower = std::abs(source - factor);
        const int upper = std::min(source + factor, 2 * level - source - factor);
        for (int target = lower; target <= upper; ++target) {
            next[static_cast<std::size_t>(target)] += coefficient;
        }
    }
    return next;
}

cpp_int wedge(
    const std::vector<Vector>& columns,
    int radius,
    int left,
    int right
) {
    return columns[0][static_cast<std::size_t>(left)]
            * columns[static_cast<std::size_t>(radius)][
                static_cast<std::size_t>(right)]
        - columns[0][static_cast<std::size_t>(right)]
            * columns[static_cast<std::size_t>(radius)][
                static_cast<std::size_t>(left)];
}

cpp_int reflected_window(
    const std::vector<Vector>& columns,
    int left_radius,
    int right_radius,
    int cutoff,
    int level
) {
    cpp_int total = 0;
    for (int left = cutoff; left <= level - cutoff; ++left) {
        for (int right = left + 1; right <= level - cutoff; ++right) {
            total += wedge(columns, left_radius, left, right)
                * wedge(columns, right_radius, left, right);
        }
    }
    return total;
}

cpp_int lower_tail(
    const std::vector<Vector>& columns,
    int left_radius,
    int right_radius,
    int cutoff,
    int level
) {
    cpp_int total = 0;
    for (int left = cutoff; left <= level; ++left) {
        for (int right = left + 1; right <= level; ++right) {
            total += wedge(columns, left_radius, left, right)
                * wedge(columns, right_radius, left, right);
        }
    }
    return total;
}

void print_failure(const char* label, const Failure& failure) {
    std::cout
        << ' ' << label << "_level=" << 2 * failure.level
        << ' ' << label << "_factor=" << 2 * failure.factor
        << ' ' << label << "_half_power=" << failure.half_power
        << ' ' << label << "_left_radius=" << 2 * failure.left_radius
        << ' ' << label << "_right_radius=" << 2 * failure.right_radius
        << ' ' << label << "_cutoff=" << failure.cutoff
        << ' ' << label << "_value=" << failure.value;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? parse_positive(argv[1], "maximum_half_level")
            : 20;
        const int maximum_half_power = argc >= 3
            ? parse_positive(argv[2], "maximum_half_power")
            : 15;
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_finite_half_power_wedge_windows "
                "[maximum_half_level [maximum_half_power]]");
        }

        std::uint64_t cases = 0U;
        std::uint64_t window_checks = 0U;
        std::uint64_t negative_windows = 0U;
        std::uint64_t anchored_window_checks = 0U;
        std::uint64_t negative_anchored_windows = 0U;
        std::uint64_t lower_tail_checks = 0U;
        std::uint64_t negative_lower_tails = 0U;
        std::uint64_t anchored_lower_tail_checks = 0U;
        std::uint64_t negative_anchored_lower_tails = 0U;
        std::uint64_t full_gram_checks = 0U;
        std::uint64_t negative_full_gram = 0U;
        Failure first_negative_window;
        Failure first_negative_anchored_window;
        Failure first_negative_lower_tail;
        Failure first_negative_anchored_lower_tail;

        for (int level = 1; level <= maximum_level; ++level) {
            for (int factor = 1; factor <= level; ++factor) {
                Vector profile(static_cast<std::size_t>(level + 1), 0);
                profile[0] = 1;
                for (int half_power = 1;
                     half_power <= maximum_half_power;
                     ++half_power) {
                    profile = fusion_step(profile, level, factor);
                    std::vector<Vector> columns;
                    columns.reserve(static_cast<std::size_t>(level + 1));
                    for (int radius = 0; radius <= level; ++radius) {
                        columns.push_back(
                            fusion_step(profile, level, radius));
                    }

                    for (int left_radius = 1;
                         left_radius <= level;
                         ++left_radius) {
                        for (int right_radius = 1;
                             right_radius <= level;
                             ++right_radius) {
                            ++cases;
                            for (int cutoff = 0;
                                 cutoff < level;
                                 ++cutoff) {
                                const cpp_int value = lower_tail(
                                    columns,
                                    left_radius,
                                    right_radius,
                                    cutoff,
                                    level);
                                ++lower_tail_checks;
                                if (value < 0) {
                                    ++negative_lower_tails;
                                    if (
                                        first_negative_lower_tail.level
                                        < 0) {
                                        first_negative_lower_tail = {
                                            level,
                                            factor,
                                            half_power,
                                            left_radius,
                                            right_radius,
                                            cutoff,
                                            value};
                                    }
                                }
                                if (left_radius == factor) {
                                    ++anchored_lower_tail_checks;
                                    if (value < 0) {
                                        ++negative_anchored_lower_tails;
                                        if (
                                            first_negative_anchored_lower_tail
                                                .level
                                            < 0) {
                                            first_negative_anchored_lower_tail
                                                = {
                                                    level,
                                                    factor,
                                                    half_power,
                                                    left_radius,
                                                    right_radius,
                                                    cutoff,
                                                    value};
                                        }
                                    }
                                }
                            }
                            for (int cutoff = 0;
                                 cutoff <= level / 2;
                                 ++cutoff) {
                                const cpp_int value = reflected_window(
                                    columns,
                                    left_radius,
                                    right_radius,
                                    cutoff,
                                    level);
                                ++window_checks;
                                if (cutoff == 0) {
                                    ++full_gram_checks;
                                    if (value < 0) {
                                        ++negative_full_gram;
                                    }
                                }
                                if (value < 0) {
                                    ++negative_windows;
                                    if (first_negative_window.level < 0) {
                                        first_negative_window = {
                                            level,
                                            factor,
                                            half_power,
                                            left_radius,
                                            right_radius,
                                            cutoff,
                                            value};
                                    }
                                }
                                if (left_radius == factor) {
                                    ++anchored_window_checks;
                                    if (value < 0) {
                                        ++negative_anchored_windows;
                                        if (
                                            first_negative_anchored_window
                                                .level
                                            < 0) {
                                            first_negative_anchored_window = {
                                                level,
                                                factor,
                                                half_power,
                                                left_radius,
                                                right_radius,
                                                cutoff,
                                                value};
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
            << "SU2_FINITE_HALF_POWER_WEDGE_WINDOWS"
            << " maximum_level=" << 2 * maximum_level
            << " maximum_half_power=" << maximum_half_power
            << " cases=" << cases
            << " window_checks=" << window_checks
            << " negative_windows=" << negative_windows
            << " anchored_window_checks=" << anchored_window_checks
            << " negative_anchored_windows="
            << negative_anchored_windows
            << " lower_tail_checks=" << lower_tail_checks
            << " negative_lower_tails=" << negative_lower_tails
            << " anchored_lower_tail_checks="
            << anchored_lower_tail_checks
            << " negative_anchored_lower_tails="
            << negative_anchored_lower_tails
            << " full_gram_checks=" << full_gram_checks
            << " negative_full_gram=" << negative_full_gram;
        if (first_negative_window.level >= 0) {
            print_failure("first_negative_window", first_negative_window);
        }
        if (first_negative_anchored_window.level >= 0) {
            print_failure(
                "first_negative_anchored_window",
                first_negative_anchored_window);
        }
        if (first_negative_lower_tail.level >= 0) {
            print_failure(
                "first_negative_lower_tail",
                first_negative_lower_tail);
        }
        if (first_negative_anchored_lower_tail.level >= 0) {
            print_failure(
                "first_negative_anchored_lower_tail",
                first_negative_anchored_lower_tail);
        }
        std::cout
            << " full_result="
            << (negative_full_gram == 0U
                    ? "NO_NEGATIVE_FULL_GRAM"
                    : "NEGATIVE_FULL_GRAM")
            << " window_result="
            << (negative_windows == 0U
                    ? "NO_NEGATIVE_REFLECTED_WINDOW"
                    : "NEGATIVE_REFLECTED_WINDOW")
            << " lower_tail_result="
            << (negative_lower_tails == 0U
                    ? "NO_NEGATIVE_LOWER_TAIL"
                    : "NEGATIVE_LOWER_TAIL")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
