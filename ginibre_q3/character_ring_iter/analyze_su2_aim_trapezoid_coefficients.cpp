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

bool selected_slope(
    int index,
    int half_level,
    int shell,
    int label) {
    if (index < 0) {
        return false;
    }
    const int period = 2 * half_level + 2;
    if (label == half_level) {
        return index >= shell * period + half_level
            && index % period == half_level;
    }
    if (index < shell * period + label) {
        return false;
    }
    const int residue = index % period;
    return residue == label || residue == period - label - 2;
}

int boundary_potential(
    int index,
    int half_level,
    int shell,
    int factor,
    int radius) {
    const int period = 2 * half_level + 2;
    const int reach = 2 * factor - radius;
    if (reach <= 0 || index < shell * period - reach) {
        return 0;
    }
    return std::min(2 * reach, index - (shell * period - reach) + 1);
}

long long trapezoid_sum(
    const std::vector<long long>& values,
    int width,
    int height,
    int shift) {
    long long result = 0;
    for (int left = 0; left < width; ++left) {
        for (int right = 0; right < height; ++right) {
            result += values[static_cast<std::size_t>(
                shift + left + right)];
        }
    }
    return result;
}

long long simplex_pairs(long long extent) {
    return extent <= 1 ? 0 : extent * (extent - 1) / 2;
}

long long lower_pair_count(int width, int height, int cut) {
    return simplex_pairs(static_cast<long long>(cut) + 1)
        - simplex_pairs(static_cast<long long>(cut - width) + 1)
        - simplex_pairs(static_cast<long long>(cut - height) + 1)
        + simplex_pairs(
            static_cast<long long>(cut - width - height) + 1);
}

long long stable_weight(
    int factor,
    int radius,
    int label) {
    if (radius < factor) {
        if (radius + 1 <= label
            && label <= 2 * factor - radius - 1) {
            return 2;
        }
        if (2 * factor - radius <= label
            && label <= 2 * factor + radius + 1) {
            return 1;
        }
    } else if (radius < 2 * factor
               && 3 * radius - 2 * factor + 2 <= label
               && label <= 2 * factor + radius + 1) {
        return 1;
    }
    return 0;
}

long long stable_coefficient(
    int factor,
    int radius,
    int label) {
    long long result = 0;
    for (int source = 0; source <= label; ++source) {
        result += stable_weight(factor, radius, source);
    }
    return result;
}

struct Witness {
    int half_level = -1;
    int shell = -1;
    int factor = -1;
    int radius = -1;
    int width = -1;
    int height = -1;
    int shift = -1;
    int label = -1;
    long long value = 0;
};

std::string render(const Witness& witness) {
    return "half_level=" + std::to_string(witness.half_level)
        + " shell=" + std::to_string(witness.shell)
        + " factor=" + std::to_string(witness.factor)
        + " radius=" + std::to_string(witness.radius)
        + " width=" + std::to_string(witness.width)
        + " height=" + std::to_string(witness.height)
        + " shift=" + std::to_string(witness.shift)
        + " label=" + std::to_string(witness.label)
        + " value=" + std::to_string(witness.value);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_half_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 20;
        if (argc > 2 || maximum_half_level < 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim_trapezoid_coefficients "
                "[maximum_half_level]");
        }

        std::uint64_t cases = 0U;
        std::uint64_t negative_reserve_counts = 0U;
        std::uint64_t mixed_reserve_counts = 0U;
        std::uint64_t stable_failures = 0U;
        std::uint64_t box_failures = 0U;
        std::uint64_t point_interval_checks = 0U;
        std::uint64_t point_interval_stable_failures = 0U;
        std::uint64_t point_interval_box_failures = 0U;
        std::uint64_t ascending_atom_failures = 0U;
        std::uint64_t descending_atom_failures = 0U;
        std::uint64_t sorted_atom_failures = 0U;
        long long minimum_stable_margin
            = std::numeric_limits<long long>::max();
        long long minimum_box_margin
            = std::numeric_limits<long long>::max();
        Witness first_negative;
        Witness first_mixed;
        Witness first_stable_failure;
        Witness first_box_failure;
        Witness first_point_interval_stable_failure;
        Witness first_point_interval_box_failure;
        Witness first_ascending_atom_failure;
        Witness first_descending_atom_failure;
        Witness first_sorted_atom_failure;

        for (int half_level = 2;
             half_level <= maximum_half_level;
             ++half_level) {
            const int period = 2 * half_level + 2;
            const int horizon = 4 * period;
            for (int shell = 1; shell <= 2; ++shell) {
                for (int factor = 1;
                     factor <= half_level / 2;
                     ++factor) {
                    for (int radius = 0;
                         radius < 2 * factor;
                         ++radius) {
                        const int reach = 2 * factor - radius;
                        std::vector<int> stable_tokens;
                        for (int label = 0;
                             static_cast<int>(stable_tokens.size())
                                 < 2 * reach;
                             ++label) {
                            const long long weight = stable_weight(
                                factor,
                                radius,
                                label);
                            for (long long copy = 0;
                                 copy < weight;
                                 ++copy) {
                                stable_tokens.push_back(label);
                            }
                        }
                        if (static_cast<int>(stable_tokens.size())
                            != 2 * reach) {
                            throw std::logic_error(
                                "stable token count mismatch");
                        }
                        std::vector<long long> boundary(
                            static_cast<std::size_t>(horizon + 1),
                            0);
                        for (int index = 0; index <= horizon; ++index) {
                            boundary[static_cast<std::size_t>(index)]
                                = boundary_potential(
                                    index,
                                    half_level,
                                    shell,
                                    factor,
                                    radius);
                        }
                        std::vector<std::vector<long long>> reserve(
                            static_cast<std::size_t>(half_level + 1),
                            std::vector<long long>(
                                static_cast<std::size_t>(horizon + 1),
                                0));
                        for (int label = 0;
                             label <= half_level;
                             ++label) {
                            for (int index = 0;
                                 index <= horizon;
                                 ++index) {
                                reserve[
                                    static_cast<std::size_t>(label)]
                                    [static_cast<std::size_t>(index)]
                                    = selected_slope(
                                        index,
                                        half_level,
                                        shell,
                                        label);
                            }
                        }

                        for (int width = 1; width < period; ++width) {
                            for (int shift = width;
                                 shift < width + period;
                                 shift += 2) {
                                for (int initial = 0;
                                     initial < period - 1;
                                     ++initial) {
                                    ++point_interval_checks;
                                    long long boundary_margin = 0;
                                    std::vector<long long> deltas(
                                        static_cast<std::size_t>(
                                            half_level + 1),
                                        0);
                                    for (int offset = 0;
                                         offset < width;
                                         ++offset) {
                                        const int old_index
                                            = initial + offset;
                                        const int new_index
                                            = old_index + shift;
                                        boundary_margin += boundary[
                                            static_cast<std::size_t>(
                                                new_index)]
                                            - boundary[
                                                static_cast<std::size_t>(
                                                    old_index)];
                                        for (int label = 0;
                                             label <= half_level;
                                             ++label) {
                                            deltas[
                                                static_cast<std::size_t>(
                                                    label)]
                                                += reserve[
                                                      static_cast<
                                                          std::size_t>(
                                                          label)]
                                                      [static_cast<
                                                          std::size_t>(
                                                          new_index)]
                                                  - reserve[
                                                      static_cast<
                                                          std::size_t>(
                                                          label)]
                                                      [static_cast<
                                                          std::size_t>(
                                                          old_index)];
                                        }
                                    }
                                    long long stable_margin
                                        = boundary_margin;
                                    long long box_margin
                                        = boundary_margin;
                                    long long suffix_delta = 0;
                                    std::vector<long long> suffix_deltas(
                                        static_cast<std::size_t>(
                                            half_level + 2),
                                        0);
                                    for (int label = half_level;
                                         label >= 0;
                                         --label) {
                                        const long long delta
                                            = deltas[
                                                static_cast<std::size_t>(
                                                    label)];
                                        stable_margin
                                            -= stable_coefficient(
                                                   factor,
                                                   radius,
                                                   label)
                                                * delta;
                                        suffix_delta += delta;
                                        suffix_deltas[
                                            static_cast<std::size_t>(
                                                label)] = suffix_delta;
                                        box_margin
                                            -= stable_weight(
                                                   factor,
                                                   radius,
                                                   label)
                                                * std::max(
                                                    0LL,
                                                    suffix_delta);
                                    }
                                    if (stable_margin < 0) {
                                        ++point_interval_stable_failures;
                                        if (
                                            first_point_interval_stable_failure
                                                .half_level
                                            < 0) {
                                            first_point_interval_stable_failure
                                                = {
                                                    half_level,
                                                    shell,
                                                    factor,
                                                    radius,
                                                    width,
                                                    initial,
                                                    shift,
                                                    -1,
                                                    stable_margin};
                                        }
                                    }
                                    if (box_margin < 0) {
                                        ++point_interval_box_failures;
                                        if (
                                            first_point_interval_box_failure
                                                .half_level
                                            < 0) {
                                            first_point_interval_box_failure
                                                = {
                                                    half_level,
                                                    shell,
                                                    factor,
                                                    radius,
                                                    width,
                                                    initial,
                                                    shift,
                                                    -1,
                                                    box_margin};
                                        }
                                    }
                                }
                            }
                            for (int height = 1;
                                 height < period;
                                 ++height) {
                                const long long boundary_initial
                                    = trapezoid_sum(
                                        boundary,
                                        width,
                                        height,
                                        0);
                                std::vector<long long> reserve_initial(
                                    static_cast<std::size_t>(
                                        half_level + 1),
                                    0);
                                for (int label = 0;
                                     label <= half_level;
                                     ++label) {
                                    reserve_initial[
                                        static_cast<std::size_t>(label)]
                                        = trapezoid_sum(
                                            reserve[
                                                static_cast<std::size_t>(
                                                    label)],
                                            width,
                                            height,
                                            0);
                                }
                                for (int shift = width;
                                     shift < width + period;
                                     shift += 2) {
                                    ++cases;
                                    const long long boundary_margin
                                        = trapezoid_sum(
                                              boundary,
                                              width,
                                              height,
                                              shift)
                                          - boundary_initial;
                                    long long stable_margin
                                        = boundary_margin;
                                    long long box_margin
                                        = boundary_margin;
                                    bool has_positive = false;
                                    bool has_negative = false;
                                    int mixed_label = -1;
                                    long long mixed_value = 0;
                                    std::vector<long long> deltas(
                                        static_cast<std::size_t>(
                                            half_level + 1),
                                        0);
                                    for (int label = 0;
                                         label <= half_level;
                                         ++label) {
                                        const long long delta
                                            = trapezoid_sum(
                                                  reserve[
                                                      static_cast<
                                                          std::size_t>(
                                                          label)],
                                                  width,
                                                  height,
                                                  shift)
                                              - reserve_initial[
                                                    static_cast<
                                                        std::size_t>(
                                                        label)];
                                        deltas[
                                            static_cast<std::size_t>(
                                                label)] = delta;
                                        has_positive
                                            = has_positive || delta > 0;
                                        has_negative
                                            = has_negative || delta < 0;
                                        if (delta < 0) {
                                            ++negative_reserve_counts;
                                            if (
                                                first_negative.half_level
                                                < 0) {
                                                first_negative = {
                                                    half_level,
                                                    shell,
                                                    factor,
                                                    radius,
                                                    width,
                                                    height,
                                                    shift,
                                                    label,
                                                    delta};
                                            }
                                        }
                                        if (delta < 0
                                            && mixed_label < 0) {
                                            mixed_label = label;
                                            mixed_value = delta;
                                        }
                                    }
                                    long long suffix_delta = 0;
                                    std::vector<long long> suffix_deltas(
                                        static_cast<std::size_t>(
                                            half_level + 2),
                                        0);
                                    for (int label = half_level;
                                         label >= 0;
                                         --label) {
                                        const long long delta
                                            = deltas[
                                                static_cast<std::size_t>(
                                                    label)];
                                        stable_margin
                                            -= stable_coefficient(
                                                   factor,
                                                   radius,
                                                   label)
                                                * delta;
                                        suffix_delta += delta;
                                        suffix_deltas[
                                            static_cast<std::size_t>(
                                                label)] = suffix_delta;
                                        box_margin
                                            -= stable_weight(
                                                   factor,
                                                   radius,
                                                   label)
                                                * std::max(
                                                    0LL,
                                                    suffix_delta);
                                    }
                                    for (int token = 0;
                                         token < 2 * reach;
                                         ++token) {
                                        const int label
                                            = stable_tokens[
                                                static_cast<std::size_t>(
                                                    token)];
                                        const long long reserve_delta
                                            = label <= half_level
                                            ? suffix_deltas[
                                                  static_cast<std::size_t>(
                                                      label)]
                                            : 0;
                                        const auto atom_margin =
                                            [&](int threshold_offset) {
                                                const int threshold
                                                    = shell * period
                                                      + threshold_offset;
                                                return lower_pair_count(
                                                           width,
                                                           height,
                                                           threshold)
                                                    - lower_pair_count(
                                                          width,
                                                          height,
                                                          threshold
                                                              - shift)
                                                    - reserve_delta;
                                            };
                                        const long long ascending_margin
                                            = atom_margin(
                                                -reach + token);
                                        if (ascending_margin < 0) {
                                            ++ascending_atom_failures;
                                            if (
                                                first_ascending_atom_failure
                                                    .half_level
                                                < 0) {
                                                first_ascending_atom_failure
                                                    = {
                                                        half_level,
                                                        shell,
                                                        factor,
                                                        radius,
                                                        width,
                                                        height,
                                                        shift,
                                                        label,
                                                        ascending_margin};
                                            }
                                        }
                                        const long long descending_margin
                                            = atom_margin(
                                                reach - 1 - token);
                                        if (descending_margin < 0) {
                                            ++descending_atom_failures;
                                            if (
                                                first_descending_atom_failure
                                                    .half_level
                                                < 0) {
                                                first_descending_atom_failure
                                                    = {
                                                        half_level,
                                                        shell,
                                                        factor,
                                                        radius,
                                                        width,
                                                        height,
                                                        shift,
                                                        label,
                                                        descending_margin};
                                            }
                                        }
                                    }
                                    std::vector<long long> supplies;
                                    std::vector<long long> demands;
                                    supplies.reserve(
                                        static_cast<std::size_t>(
                                            2 * reach));
                                    demands.reserve(
                                        static_cast<std::size_t>(
                                            2 * reach));
                                    for (int token = 0;
                                         token < 2 * reach;
                                         ++token) {
                                        const int threshold
                                            = shell * period - reach
                                              + token;
                                        supplies.push_back(
                                            lower_pair_count(
                                                width,
                                                height,
                                                threshold)
                                            - lower_pair_count(
                                                width,
                                                height,
                                                threshold - shift));
                                        const int label
                                            = stable_tokens[
                                                static_cast<std::size_t>(
                                                    token)];
                                        demands.push_back(
                                            label <= half_level
                                            ? suffix_deltas[
                                                  static_cast<std::size_t>(
                                                      label)]
                                            : 0);
                                    }
                                    std::sort(
                                        supplies.begin(),
                                        supplies.end());
                                    std::sort(
                                        demands.begin(),
                                        demands.end());
                                    for (int token = 0;
                                         token < 2 * reach;
                                         ++token) {
                                        const long long margin
                                            = supplies[
                                                  static_cast<std::size_t>(
                                                      token)]
                                              - demands[
                                                    static_cast<
                                                        std::size_t>(
                                                        token)];
                                        if (margin >= 0) {
                                            continue;
                                        }
                                        ++sorted_atom_failures;
                                        if (
                                            first_sorted_atom_failure
                                                .half_level
                                            < 0) {
                                            first_sorted_atom_failure = {
                                                half_level,
                                                shell,
                                                factor,
                                                radius,
                                                width,
                                                height,
                                                shift,
                                                token,
                                                margin};
                                        }
                                        break;
                                    }
                                    if (has_positive && has_negative) {
                                        ++mixed_reserve_counts;
                                        if (first_mixed.half_level < 0) {
                                            first_mixed = {
                                                half_level,
                                                shell,
                                                factor,
                                                radius,
                                                width,
                                                height,
                                                shift,
                                                mixed_label,
                                                mixed_value};
                                        }
                                    }
                                    minimum_stable_margin = std::min(
                                        minimum_stable_margin,
                                        stable_margin);
                                    minimum_box_margin = std::min(
                                        minimum_box_margin,
                                        box_margin);
                                    if (stable_margin < 0) {
                                        ++stable_failures;
                                        if (
                                            first_stable_failure.half_level
                                            < 0) {
                                            first_stable_failure = {
                                                half_level,
                                                shell,
                                                factor,
                                                radius,
                                                width,
                                                height,
                                                shift,
                                                -1,
                                                stable_margin};
                                        }
                                    }
                                    if (box_margin < 0) {
                                        ++box_failures;
                                        if (
                                            first_box_failure.half_level
                                            < 0) {
                                            first_box_failure = {
                                                half_level,
                                                shell,
                                                factor,
                                                radius,
                                                width,
                                                height,
                                                shift,
                                                -1,
                                                box_margin};
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
            << "SU2_AIM_TRAPEZOID_COEFFICIENTS"
            << " maximum_half_level=" << maximum_half_level
            << " cases=" << cases
            << " negative_reserve_counts="
            << negative_reserve_counts
            << " mixed_reserve_counts=" << mixed_reserve_counts
            << " stable_failures=" << stable_failures
            << " box_failures=" << box_failures
            << " point_interval_checks=" << point_interval_checks
            << " point_interval_stable_failures="
            << point_interval_stable_failures
            << " point_interval_box_failures="
            << point_interval_box_failures
            << " ascending_atom_failures="
            << ascending_atom_failures
            << " descending_atom_failures="
            << descending_atom_failures
            << " sorted_atom_failures="
            << sorted_atom_failures
            << " minimum_stable_margin=" << minimum_stable_margin
            << " minimum_box_margin=" << minimum_box_margin
            << '\n'
            << "FIRST_NEGATIVE " << render(first_negative) << '\n'
            << "FIRST_MIXED " << render(first_mixed) << '\n'
            << "FIRST_STABLE_FAILURE "
            << render(first_stable_failure) << '\n'
            << "FIRST_BOX_FAILURE "
            << render(first_box_failure) << '\n'
            << "FIRST_POINT_INTERVAL_STABLE_FAILURE "
            << render(first_point_interval_stable_failure) << '\n'
            << "FIRST_POINT_INTERVAL_BOX_FAILURE "
            << render(first_point_interval_box_failure) << '\n'
            << "FIRST_ASCENDING_ATOM_FAILURE "
            << render(first_ascending_atom_failure) << '\n'
            << "FIRST_DESCENDING_ATOM_FAILURE "
            << render(first_descending_atom_failure) << '\n'
            << "FIRST_SORTED_ATOM_FAILURE "
            << render(first_sorted_atom_failure) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
