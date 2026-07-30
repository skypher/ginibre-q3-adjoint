#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using Integer = std::int64_t;
using Vector = std::vector<Integer>;

std::uint64_t positive_argument(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const unsigned long long parsed = std::stoull(value, &consumed, 10);
    if (consumed != value.size() || parsed == 0U) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<std::uint64_t>(parsed);
}

std::uint64_t splitmix64(std::uint64_t& state) {
    state += UINT64_C(0x9e3779b97f4a7c15);
    std::uint64_t value = state;
    value = (value ^ (value >> 30U))
        * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U))
        * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

Vector signed_profile(
    int level,
    int maximum_coefficient,
    std::uint64_t sample) {
    std::uint64_t state
        = sample ^ UINT64_C(0xd1b54a32d192ed03);
    const int mode = static_cast<int>(splitmix64(state) % 4U);
    Vector profile(static_cast<std::size_t>(level + 1), 0);
    for (int index = 0; index <= level; ++index) {
        const Integer magnitude = static_cast<Integer>(
            splitmix64(state)
            % static_cast<std::uint64_t>(maximum_coefficient + 1));
        int sign = 1;
        if (mode == 0) {
            sign = (splitmix64(state) & 1U) == 0U ? -1 : 1;
        } else if (mode == 1) {
            sign = index == static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(level + 1))
                ? -1
                : 1;
        } else if (mode == 2) {
            sign = index % 2 == 0 ? 1 : -1;
        } else {
            const int cut = static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(level + 1));
            sign = index < cut ? 1 : -1;
        }
        profile[static_cast<std::size_t>(index)]
            = static_cast<Integer>(sign) * magnitude;
    }
    return profile;
}

Vector fusion_square(const Vector& profile, int level) {
    Vector square(static_cast<std::size_t>(level + 1), 0);
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            const Integer weight
                = profile[static_cast<std::size_t>(left)]
                  * profile[static_cast<std::size_t>(right)];
            if (weight == 0) {
                continue;
            }
            const int lower = std::abs(left - right);
            const int upper
                = std::min(left + right, 2 * level - left - right);
            for (int target = lower; target <= upper; ++target) {
                square[static_cast<std::size_t>(target)] += weight;
            }
        }
    }
    return square;
}

Vector fusion_step(const Vector& profile, int level, int factor) {
    Vector output(static_cast<std::size_t>(level + 1), 0);
    for (int source = 0; source <= level; ++source) {
        const Integer coefficient
            = profile[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        const int lower = std::abs(source - factor);
        const int upper
            = std::min(source + factor, 2 * level - source - factor);
        for (int target = lower; target <= upper; ++target) {
            output[static_cast<std::size_t>(target)] += coefficient;
        }
    }
    return output;
}

bool nonnegative_log_concave(const Vector& vector) {
    bool support_started = false;
    bool support_ended = false;
    for (std::size_t index = 0U; index < vector.size(); ++index) {
        if (vector[index] < 0) {
            return false;
        }
        if (vector[index] > 0) {
            if (support_ended) {
                return false;
            }
            support_started = true;
        } else if (support_started) {
            support_ended = true;
        }
        if (index > 0U && index + 1U < vector.size()
            && vector[index] * vector[index]
                < vector[index - 1U] * vector[index + 1U]) {
            return false;
        }
    }
    return support_started && vector[0] > 0;
}

bool boundary_admissible(const Vector& vector) {
    const std::size_t level = vector.size() - 1U;
    for (std::size_t radius = 0U; radius <= level; ++radius) {
        if (vector[0] * vector[level - radius]
            < vector[radius] * vector[level]) {
            return false;
        }
    }
    return true;
}

std::tuple<int, int, Integer> first_bad_current(
    const Vector& vector,
    int level) {
    for (int radius = 0; radius <= level; ++radius) {
        for (int target = 0; target <= level; ++target) {
            Integer translated = 0;
            const int lower = std::abs(radius - target);
            const int upper
                = std::min(
                    radius + target,
                    2 * level - radius - target);
            for (int label = lower; label <= upper; ++label) {
                translated += vector[static_cast<std::size_t>(label)];
            }
            const Integer current
                = vector[0] * translated
                  - vector[static_cast<std::size_t>(radius)]
                    * vector[static_cast<std::size_t>(target)];
            if (current < 0) {
                return {radius, target, current};
            }
        }
    }
    return {-1, -1, 0};
}

std::string render(const Vector& vector) {
    std::string text = "[";
    for (std::size_t index = 0U; index < vector.size(); ++index) {
        if (index != 0U) {
            text += ',';
        }
        text += std::to_string(vector[index]);
    }
    return text + ']';
}

struct Failure {
    std::uint64_t sample = std::numeric_limits<std::uint64_t>::max();
    int level = -1;
    int radius = -1;
    int target = -1;
    Integer value = 0;
    Vector root;
    Vector square;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::uint64_t samples = argc >= 2
            ? positive_argument(argv[1], "samples")
            : UINT64_C(1000000);
        const int maximum_level = argc >= 3
            ? static_cast<int>(positive_argument(argv[2], "maximum_level"))
            : 16;
        const int maximum_coefficient = argc >= 4
            ? static_cast<int>(
                positive_argument(argv[3], "maximum_coefficient"))
            : 6;
        if (argc > 4 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: probe_su2_signed_square_boundary_stress "
                "[samples] [maximum_half_level] [maximum_coefficient]");
        }

        std::uint64_t nonnegative_squares = 0U;
        std::uint64_t log_concave_squares = 0U;
        std::uint64_t boundary_squares = 0U;
        std::uint64_t current_failures = 0U;
        std::uint64_t insertion_checks = 0U;
        std::uint64_t insertion_failures = 0U;
        std::uint64_t full_insertion_failures = 0U;
        Failure first;
        Failure first_insertion;
        Failure first_full_insertion;
        int first_insertion_factor = -1;
        int first_full_insertion_factor = -1;

#pragma omp parallel for schedule(static) \
    reduction(+:nonnegative_squares,log_concave_squares,boundary_squares,current_failures,insertion_checks,insertion_failures,full_insertion_failures)
        for (std::int64_t signed_sample = 0;
             signed_sample < static_cast<std::int64_t>(samples);
             ++signed_sample) {
            const std::uint64_t sample
                = static_cast<std::uint64_t>(signed_sample);
            std::uint64_t state
                = sample ^ UINT64_C(0xa0761d6478bd642f);
            const int level = 2 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_level - 1));
            const Vector root
                = signed_profile(level, maximum_coefficient, sample);
            const Vector square = fusion_square(root, level);
            if (std::ranges::any_of(
                    square, [](Integer value) { return value < 0; })) {
                continue;
            }
            ++nonnegative_squares;
            if (!nonnegative_log_concave(square)) {
                continue;
            }
            ++log_concave_squares;
            if (!boundary_admissible(square)) {
                continue;
            }
            ++boundary_squares;
            const auto [radius, target, value]
                = first_bad_current(square, level);
            if (radius < 0) {
                for (int factor = 1; factor <= level; ++factor) {
                    ++insertion_checks;
                    const Vector output = fusion_step(
                        fusion_step(square, level, factor),
                        level,
                        factor);
                    const auto [
                        full_radius,
                        full_target,
                        full_value] = first_bad_current(output, level);
                    if (full_radius >= 0) {
                        ++full_insertion_failures;
#pragma omp critical
                        {
                            if (sample < first_full_insertion.sample) {
                                first_full_insertion = {
                                    sample,
                                    level,
                                    full_radius,
                                    full_target,
                                    full_value,
                                    root,
                                    square};
                                first_full_insertion_factor = factor;
                            }
                        }
                    }
                    if (boundary_admissible(output)) {
                        continue;
                    }
                    ++insertion_failures;
                    int bad_radius = -1;
                    Integer bad_value = 0;
                    for (int insertion_radius = 0;
                         insertion_radius <= level;
                         ++insertion_radius) {
                        const Integer boundary
                            = output[0]
                                * output[static_cast<std::size_t>(
                                    level - insertion_radius)]
                              - output[
                                    static_cast<std::size_t>(
                                        insertion_radius)]
                                * output[static_cast<std::size_t>(level)];
                        if (boundary < 0) {
                            bad_radius = insertion_radius;
                            bad_value = boundary;
                            break;
                        }
                    }
#pragma omp critical
                    {
                        if (sample < first_insertion.sample) {
                            first_insertion = {
                                sample,
                                level,
                                bad_radius,
                                level,
                                bad_value,
                                root,
                                square};
                            first_insertion_factor = factor;
                        }
                    }
                }
            } else {
                ++current_failures;
#pragma omp critical
                {
                    if (sample < first.sample) {
                        first = {
                            sample,
                            level,
                            radius,
                            target,
                            value,
                            root,
                            square};
                    }
                }
            }
        }

        std::cout
            << "SU2_SIGNED_SQUARE_BOUNDARY_STRESS"
            << " samples=" << samples
            << " maximum_level=" << 2 * maximum_level
            << " maximum_coefficient=" << maximum_coefficient
            << " threads=" << omp_get_max_threads()
            << " nonnegative_squares=" << nonnegative_squares
            << " log_concave_squares=" << log_concave_squares
            << " boundary_squares=" << boundary_squares
            << " current_failures=" << current_failures
            << " insertion_checks=" << insertion_checks
            << " insertion_failures=" << insertion_failures
            << " full_insertion_failures=" << full_insertion_failures
            << " first_failure=("
            << (first.level < 0 ? -1 : 2 * first.level) << ','
            << (first.radius < 0 ? -1 : 2 * first.radius) << ','
            << (first.target < 0 ? -1 : 2 * first.target) << ','
            << first.value << ')'
            << " first_root=" << render(first.root)
            << " first_square=" << render(first.square)
            << " first_insertion=("
            << (first_insertion.level < 0
                    ? -1
                    : 2 * first_insertion.level) << ','
            << (first_insertion_factor < 0
                    ? -1
                    : 2 * first_insertion_factor) << ','
            << (first_insertion.radius < 0
                    ? -1
                    : 2 * first_insertion.radius) << ','
            << first_insertion.value << ')'
            << " first_insertion_root="
            << render(first_insertion.root)
            << " first_insertion_input="
            << render(first_insertion.square)
            << " first_full_insertion=("
            << (first_full_insertion.level < 0
                    ? -1
                    : 2 * first_full_insertion.level) << ','
            << (first_full_insertion_factor < 0
                    ? -1
                    : 2 * first_full_insertion_factor) << ','
            << (first_full_insertion.radius < 0
                    ? -1
                    : 2 * first_full_insertion.radius) << ','
            << (first_full_insertion.target < 0
                    ? -1
                    : 2 * first_full_insertion.target) << ','
            << first_full_insertion.value << ')'
            << " first_full_insertion_root="
            << render(first_full_insertion.root)
            << " first_full_insertion_input="
            << render(first_full_insertion.square)
            << " result="
            << (current_failures == 0U
                    && insertion_failures == 0U
                    && full_insertion_failures == 0U
                    ? "NO_PSD_CURRENT_OR_INSERTION_FAILURE"
                    : "PSD_CURRENT_OR_INSERTION_FAILURE")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
