#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Vector = std::vector<cpp_int>;

std::uint64_t parse_positive(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const unsigned long long parsed = std::stoull(value, &consumed, 10);
    if (consumed != value.size() || parsed == 0ULL) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<std::uint64_t>(parsed);
}

std::uint64_t splitmix64(std::uint64_t& state) {
    state += UINT64_C(0x9e3779b97f4a7c15);
    std::uint64_t value = state;
    value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

Vector profile(std::uint64_t sample, int maximum_level) {
    std::uint64_t state
        = sample ^ UINT64_C(0xd1b54a32d192ed03);
    const int level = 2 + static_cast<int>(
        splitmix64(state) % static_cast<std::uint64_t>(maximum_level - 1));
    Vector values(static_cast<std::size_t>(level + 1), 0);
    const int support_left = static_cast<int>(
        splitmix64(state) % static_cast<std::uint64_t>(level + 1));
    const int maximum_width = level - support_left;
    const int support_right = support_left + static_cast<int>(
        splitmix64(state)
        % static_cast<std::uint64_t>(maximum_width + 1));
    values[static_cast<std::size_t>(support_left)]
        = 1U + splitmix64(state) % 50U;
    if (support_left == support_right) {
        return values;
    }
    values[static_cast<std::size_t>(support_left + 1)]
        = 1U + splitmix64(state) % 250U;
    constexpr std::uint64_t coefficient_cap = UINT64_C(1000000);
    for (int index = support_left + 1; index < support_right; ++index) {
        const cpp_int quotient
            = values[static_cast<std::size_t>(index)]
                * values[static_cast<std::size_t>(index)]
                / values[static_cast<std::size_t>(index - 1)];
        const std::uint64_t maximum
            = quotient > coefficient_cap
            ? coefficient_cap
            : quotient.convert_to<std::uint64_t>();
        if (maximum == 0U) {
            break;
        }
        values[static_cast<std::size_t>(index + 1)]
            = 1U + splitmix64(state) % maximum;
    }
    return values;
}

Vector fusion_square(const Vector& values) {
    const int level = static_cast<int>(values.size()) - 1;
    Vector difference(static_cast<std::size_t>(level + 2), 0);
    for (int left = 0; left <= level; ++left) {
        if (values[static_cast<std::size_t>(left)] == 0) {
            continue;
        }
        for (int right = 0; right <= level; ++right) {
            if (values[static_cast<std::size_t>(right)] == 0) {
                continue;
            }
            const int lower = std::abs(left - right);
            const int upper
                = std::min(left + right, 2 * level - left - right);
            const cpp_int weight
                = values[static_cast<std::size_t>(left)]
                    * values[static_cast<std::size_t>(right)];
            difference[static_cast<std::size_t>(lower)] += weight;
            difference[static_cast<std::size_t>(upper + 1)] -= weight;
        }
    }
    Vector square(static_cast<std::size_t>(level + 1), 0);
    cpp_int running = 0;
    for (int index = 0; index <= level; ++index) {
        running += difference[static_cast<std::size_t>(index)];
        square[static_cast<std::size_t>(index)] = running;
    }
    return square;
}

Vector fusion_step(const Vector& values, int factor) {
    const int level = static_cast<int>(values.size()) - 1;
    Vector difference(static_cast<std::size_t>(level + 2), 0);
    for (int source = 0; source <= level; ++source) {
        const cpp_int& coefficient = values[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        const int lower = std::abs(source - factor);
        const int upper
            = std::min(source + factor, 2 * level - source - factor);
        difference[static_cast<std::size_t>(lower)] += coefficient;
        difference[static_cast<std::size_t>(upper + 1)] -= coefficient;
    }
    Vector output(static_cast<std::size_t>(level + 1), 0);
    cpp_int running = 0;
    for (int index = 0; index <= level; ++index) {
        running += difference[static_cast<std::size_t>(index)];
        output[static_cast<std::size_t>(index)] = running;
    }
    return output;
}

bool boundary_admissible(const Vector& square, int& bad_radius) {
    const int level = static_cast<int>(square.size()) - 1;
    for (int radius = 0; radius <= level; ++radius) {
        if (square[0] * square[static_cast<std::size_t>(level - radius)]
            < square[static_cast<std::size_t>(radius)]
                * square[static_cast<std::size_t>(level)]) {
            bad_radius = radius;
            return false;
        }
    }
    return true;
}

bool full_current_cone(
    const Vector& values,
    int* bad_radius = nullptr,
    int* bad_target = nullptr,
    cpp_int* bad_value = nullptr) {
    if (values[0] <= 0) {
        return false;
    }
    const int level = static_cast<int>(values.size()) - 1;
    Vector prefix(static_cast<std::size_t>(level + 2), 0);
    for (int index = 0; index <= level; ++index) {
        prefix[static_cast<std::size_t>(index + 1)]
            = prefix[static_cast<std::size_t>(index)]
                + values[static_cast<std::size_t>(index)];
    }
    for (int radius = 0; radius <= level; ++radius) {
        for (int target = 0; target <= level; ++target) {
            const int lower = std::abs(radius - target);
            const int upper = std::min(
                radius + target, 2 * level - radius - target);
            const cpp_int translated
                = prefix[static_cast<std::size_t>(upper + 1)]
                    - prefix[static_cast<std::size_t>(lower)];
            if (values[0] * translated
                < values[static_cast<std::size_t>(radius)]
                    * values[static_cast<std::size_t>(target)]) {
                if (bad_radius != nullptr) {
                    *bad_radius = radius;
                }
                if (bad_target != nullptr) {
                    *bad_target = target;
                }
                if (bad_value != nullptr) {
                    *bad_value
                        = values[0] * translated
                            - values[static_cast<std::size_t>(radius)]
                                * values[static_cast<std::size_t>(target)];
                }
                return false;
            }
        }
    }
    return true;
}

bool log_concave(const Vector& values) {
    bool support_started = false;
    bool support_ended = false;
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (values[index] > 0) {
            if (support_ended) {
                return false;
            }
            support_started = true;
        } else if (support_started) {
            support_ended = true;
        }
        if (index > 0U && index + 1U < values.size()
            && values[index] * values[index]
                < values[index - 1U] * values[index + 1U]) {
            return false;
        }
    }
    return true;
}

cpp_int minimum_internal_log_concavity_margin(const Vector& values) {
    std::size_t left = 0U;
    while (left < values.size() && values[left] == 0) {
        ++left;
    }
    std::size_t right = values.size();
    while (right > left && values[right - 1U] == 0) {
        --right;
    }
    if (right <= left + 2U) {
        return 0;
    }
    cpp_int minimum
        = values[left + 1U] * values[left + 1U]
            - values[left] * values[left + 2U];
    for (std::size_t index = left + 2U; index + 1U < right; ++index) {
        const cpp_int margin
            = values[index] * values[index]
                - values[index - 1U] * values[index + 1U];
        if (margin < minimum) {
            minimum = margin;
        }
    }
    return minimum;
}

std::string render(const Vector& values) {
    std::string text = "[";
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            text += ',';
        }
        text += values[index].convert_to<std::string>();
    }
    return text + ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::uint64_t samples = argc >= 2
            ? parse_positive(argv[1], "samples")
            : UINT64_C(1000000);
        const std::uint64_t parsed_level = argc >= 3
            ? parse_positive(argv[2], "maximum_half_level")
            : UINT64_C(30);
        if (argc > 3 || parsed_level < 2U
            || parsed_level
                > static_cast<std::uint64_t>(
                    std::numeric_limits<int>::max())) {
            throw std::invalid_argument(
                "usage: probe_su2_log_concave_square_boundary "
                "[samples] [maximum_half_level>=2]");
        }
        const int maximum_level = static_cast<int>(parsed_level);

        std::atomic<std::uint64_t> first_sample(samples);
        std::atomic<std::uint64_t> first_insertion_sample(samples);
        std::atomic<std::uint64_t> first_cone_boundary_sample(samples);
        std::atomic<std::uint64_t> first_boundary_reduction_sample(samples);
        std::uint64_t checked_currents = 0U;
        std::uint64_t negative_currents = 0U;
        std::uint64_t log_concave_square_profiles = 0U;
        std::uint64_t log_concave_square_currents = 0U;
        std::uint64_t negative_log_concave_square_currents = 0U;
        std::uint64_t admissible_insertions = 0U;
        std::uint64_t negative_insertions = 0U;
        std::uint64_t log_concave_cone_insertions = 0U;
        std::uint64_t negative_log_concave_cone_boundaries = 0U;
        std::uint64_t boundary_reduction_failures = 0U;
        int first_radius = -1;
        cpp_int first_value = 0;
        cpp_int first_profile_minimum_margin = 0;
        cpp_int first_square_minimum_margin = 0;
        bool first_square_log_concave = false;
        Vector first_profile;
        Vector first_square;
        int first_insertion_factor = -1;
        int first_insertion_radius = -1;
        Vector first_insertion_input;
        Vector first_insertion_output;
        int first_cone_boundary_factor = -1;
        int first_cone_boundary_radius = -1;
        Vector first_cone_boundary_input;
        Vector first_cone_boundary_output;
        int first_boundary_reduction_radius = -1;
        int first_boundary_reduction_target = -1;
        cpp_int first_boundary_reduction_value = 0;
        Vector first_boundary_reduction_square;

#pragma omp parallel for schedule(static) reduction(+:checked_currents,negative_currents,log_concave_square_profiles,log_concave_square_currents,negative_log_concave_square_currents,admissible_insertions,negative_insertions,log_concave_cone_insertions,negative_log_concave_cone_boundaries,boundary_reduction_failures)
        for (std::uint64_t sample = 0U; sample < samples; ++sample) {
            const Vector values = profile(sample, maximum_level);
            const Vector square = fusion_square(values);
            const bool square_is_log_concave = log_concave(square);
            if (square_is_log_concave) {
                ++log_concave_square_profiles;
            }
            const int level = static_cast<int>(values.size()) - 1;
            std::uint64_t factor_state
                = sample ^ UINT64_C(0xa0761d6478bd642f);
            const int factor = 1 + static_cast<int>(
                splitmix64(factor_state)
                % static_cast<std::uint64_t>(level));
            if (full_current_cone(values)) {
                ++log_concave_cone_insertions;
                const Vector cone_output = fusion_step(
                    fusion_step(values, factor), factor);
                int cone_bad_radius = -1;
                if (!boundary_admissible(
                        cone_output, cone_bad_radius)) {
                    ++negative_log_concave_cone_boundaries;
                    std::uint64_t known
                        = first_cone_boundary_sample.load();
                    while (
                        sample < known
                        && !first_cone_boundary_sample
                                .compare_exchange_weak(known, sample)) {
                    }
                    if (sample == first_cone_boundary_sample.load()) {
#pragma omp critical
                        {
                            if (
                                sample
                                == first_cone_boundary_sample.load()) {
                                first_cone_boundary_factor = factor;
                                first_cone_boundary_radius
                                    = cone_bad_radius;
                                first_cone_boundary_input = values;
                                first_cone_boundary_output = cone_output;
                            }
                        }
                    }
                }
            }
            int input_bad_radius = -1;
            if (boundary_admissible(square, input_bad_radius)) {
                ++admissible_insertions;
                int reduction_bad_radius = -1;
                int reduction_bad_target = -1;
                cpp_int reduction_bad_value = 0;
                if (!full_current_cone(
                        square,
                        &reduction_bad_radius,
                        &reduction_bad_target,
                        &reduction_bad_value)) {
                    ++boundary_reduction_failures;
                    std::uint64_t known
                        = first_boundary_reduction_sample.load();
                    while (
                        sample < known
                        && !first_boundary_reduction_sample
                                .compare_exchange_weak(known, sample)) {
                    }
                    if (sample == first_boundary_reduction_sample.load()) {
#pragma omp critical
                        {
                            if (
                                sample
                                == first_boundary_reduction_sample.load()) {
                                first_boundary_reduction_radius
                                    = reduction_bad_radius;
                                first_boundary_reduction_target
                                    = reduction_bad_target;
                                first_boundary_reduction_value
                                    = reduction_bad_value;
                                first_boundary_reduction_square = square;
                            }
                        }
                    }
                }
                const Vector inserted_profile
                    = fusion_step(values, factor);
                const Vector inserted_square
                    = fusion_square(inserted_profile);
                int output_bad_radius = -1;
                if (!boundary_admissible(
                        inserted_square, output_bad_radius)) {
                    ++negative_insertions;
                    std::uint64_t known = first_insertion_sample.load();
                    while (
                        sample < known
                        && !first_insertion_sample.compare_exchange_weak(
                            known, sample)) {
                    }
                    if (sample == first_insertion_sample.load()) {
#pragma omp critical
                        {
                            if (sample == first_insertion_sample.load()) {
                                first_insertion_factor = factor;
                                first_insertion_radius = output_bad_radius;
                                first_insertion_input = square;
                                first_insertion_output = inserted_square;
                            }
                        }
                    }
                }
            }
            for (int radius = 0; radius <= level; ++radius) {
                ++checked_currents;
                if (square_is_log_concave) {
                    ++log_concave_square_currents;
                }
                const cpp_int current
                    = square[0]
                        * square[static_cast<std::size_t>(level - radius)]
                    - square[static_cast<std::size_t>(radius)]
                        * square[static_cast<std::size_t>(level)];
                if (current >= 0) {
                    continue;
                }
                ++negative_currents;
                if (square_is_log_concave) {
                    ++negative_log_concave_square_currents;
                }
                std::uint64_t known = first_sample.load();
                while (sample < known
                       && !first_sample.compare_exchange_weak(known, sample)) {
                }
                if (sample == first_sample.load()) {
#pragma omp critical
                    {
                        if (sample == first_sample.load()) {
                            first_radius = radius;
                            first_value = current;
                            first_profile_minimum_margin
                                = minimum_internal_log_concavity_margin(
                                    values);
                            first_square_minimum_margin
                                = minimum_internal_log_concavity_margin(
                                    square);
                            first_square_log_concave
                                = square_is_log_concave;
                            first_profile = values;
                            first_square = square;
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_LOG_CONCAVE_SQUARE_BOUNDARY"
            << " samples=" << samples
            << " maximum_level=" << 2 * maximum_level
            << " checked_currents=" << checked_currents
            << " negative_currents=" << negative_currents
            << " log_concave_square_profiles="
            << log_concave_square_profiles
            << " log_concave_square_currents="
            << log_concave_square_currents
            << " negative_log_concave_square_currents="
            << negative_log_concave_square_currents
            << " admissible_insertions=" << admissible_insertions
            << " negative_insertions=" << negative_insertions
            << " log_concave_cone_insertions="
            << log_concave_cone_insertions
            << " negative_log_concave_cone_boundaries="
            << negative_log_concave_cone_boundaries
            << " boundary_reduction_failures="
            << boundary_reduction_failures
            << " first_sample="
            << (first_sample.load() == samples
                    ? std::string("-1")
                    : std::to_string(first_sample.load()))
            << " first_radius="
            << (first_radius < 0 ? -1 : 2 * first_radius)
            << " first_value=" << first_value
            << " first_profile_minimum_log_concavity_margin="
            << first_profile_minimum_margin
            << " first_square_minimum_log_concavity_margin="
            << first_square_minimum_margin
            << " first_square_log_concave="
            << (first_square_log_concave ? "yes" : "no")
            << " first_profile=" << render(first_profile)
            << " first_square=" << render(first_square)
            << " first_insertion_sample="
            << (first_insertion_sample.load() == samples
                    ? std::string("-1")
                    : std::to_string(first_insertion_sample.load()))
            << " first_insertion_factor="
            << (first_insertion_factor < 0
                    ? -1
                    : 2 * first_insertion_factor)
            << " first_insertion_radius="
            << (first_insertion_radius < 0
                    ? -1
                    : 2 * first_insertion_radius)
            << " first_insertion_input="
            << render(first_insertion_input)
            << " first_insertion_output="
            << render(first_insertion_output)
            << " first_cone_boundary_sample="
            << (first_cone_boundary_sample.load() == samples
                    ? std::string("-1")
                    : std::to_string(first_cone_boundary_sample.load()))
            << " first_cone_boundary_factor="
            << (first_cone_boundary_factor < 0
                    ? -1
                    : 2 * first_cone_boundary_factor)
            << " first_cone_boundary_radius="
            << (first_cone_boundary_radius < 0
                    ? -1
                    : 2 * first_cone_boundary_radius)
            << " first_cone_boundary_input="
            << render(first_cone_boundary_input)
            << " first_cone_boundary_output="
            << render(first_cone_boundary_output)
            << " first_boundary_reduction_sample="
            << (first_boundary_reduction_sample.load() == samples
                    ? std::string("-1")
                    : std::to_string(
                        first_boundary_reduction_sample.load()))
            << " first_boundary_reduction_radius="
            << (first_boundary_reduction_radius < 0
                    ? -1
                    : 2 * first_boundary_reduction_radius)
            << " first_boundary_reduction_target="
            << (first_boundary_reduction_target < 0
                    ? -1
                    : 2 * first_boundary_reduction_target)
            << " first_boundary_reduction_value="
            << first_boundary_reduction_value
            << " first_boundary_reduction_square="
            << render(first_boundary_reduction_square)
            << " result="
            << (negative_insertions == 0U
                    ? "NO_BOUNDARY_INSERTION_FAILURE"
                    : "BOUNDARY_INSERTION_FAILURE")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
