#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value <= 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error(
            std::string(name) + " must be a positive integer"
        );
    }
    return static_cast<int>(value);
}

std::vector<Integer> multiply_character(
    const std::vector<Integer>& state,
    int label
) {
    std::vector<Integer> next(state.size() + static_cast<std::size_t>(label));
    for (std::size_t source = 0; source < state.size(); ++source) {
        if (state[source] == 0) {
            continue;
        }
        const int source_label = static_cast<int>(source);
        const int lower = std::abs(source_label - label);
        const int upper = source_label + label;
        for (int target = lower; target <= upper; target += 2) {
            next[static_cast<std::size_t>(target)] += state[source];
        }
    }
    return next;
}

Integer coefficient_at(
    const std::vector<Integer>& state,
    std::size_t index
) {
    return index < state.size() ? state[index] : Integer(0);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: probe_su2_ordinary_anchored_kernel "
                "maximum_even_label maximum_half_power"
            );
        }
        const int maximum_label =
            parse_positive(argv[1], "maximum_even_label");
        const int maximum_half_power =
            parse_positive(argv[2], "maximum_half_power");

        std::size_t kernels = 0;
        std::size_t coordinates = 0;
        std::size_t negative = 0;
        std::size_t ratio_minors = 0;
        std::size_t negative_ratio_minors = 0;
        std::size_t negative_ratio_minors_after_base = 0;
        std::size_t even_log_concavity_minors = 0;
        std::size_t negative_even_log_concavity_minors = 0;
        std::size_t ratio_sign_changes = 0;
        std::size_t negative_to_positive_changes = 0;
        int maximum_ratio_sign_changes = 0;
        int witness_label = 0;
        int witness_power = 0;
        std::size_t witness_coordinate = 0;
        Integer minimum = 0;
        int ratio_witness_label = 0;
        int ratio_witness_power = 0;
        std::size_t ratio_witness_coordinate = 0;
        Integer minimum_ratio_minor = 0;

        for (int label = 2; label <= maximum_label; label += 2) {
            std::vector<Integer> even{Integer(1)};
            std::vector<Integer> odd =
                multiply_character(even, label);
            for (
                int power = 1;
                power <= maximum_half_power;
                ++power
            ) {
                even = multiply_character(odd, label);
                odd = multiply_character(even, label);
                const Integer even_return = coefficient_at(even, 0U);
                const Integer odd_return = coefficient_at(odd, 0U);
                ++kernels;
                for (
                    std::size_t coordinate = 1;
                    coordinate + 1 < even.size();
                    ++coordinate
                ) {
                    if ((coordinate & 1U) != 0U) {
                        continue;
                    }
                    const Integer log_concavity_minor =
                        coefficient_at(even, coordinate)
                            * coefficient_at(even, coordinate)
                        - coefficient_at(even, coordinate - 2U)
                            * coefficient_at(even, coordinate + 2U);
                    ++even_log_concavity_minors;
                    if (log_concavity_minor < 0) {
                        ++negative_even_log_concavity_minors;
                    }
                }
                Integer previous_even = 0;
                Integer previous_odd = 0;
                bool have_previous = false;
                int previous_ratio_sign = 0;
                int kernel_ratio_sign_changes = 0;
                for (
                    std::size_t coordinate = 0;
                    coordinate < odd.size();
                    ++coordinate
                ) {
                    const Integer value =
                        even_return * coefficient_at(odd, coordinate)
                        - odd_return * coefficient_at(even, coordinate);
                    ++coordinates;
                    if (value < minimum) {
                        minimum = value;
                        witness_label = label;
                        witness_power = power;
                        witness_coordinate = coordinate;
                    }
                    if (value < 0) {
                        ++negative;
                    }
                    if ((coordinate & 1U) == 0U) {
                        const Integer even_value =
                            coefficient_at(even, coordinate);
                        const Integer odd_value =
                            coefficient_at(odd, coordinate);
                        if (
                            have_previous
                            && previous_even > 0
                            && even_value > 0
                        ) {
                            const Integer ratio_minor =
                                odd_value * previous_even
                                - previous_odd * even_value;
                            ++ratio_minors;
                            if (ratio_minor < minimum_ratio_minor) {
                                minimum_ratio_minor = ratio_minor;
                                ratio_witness_label = label;
                                ratio_witness_power = power;
                                ratio_witness_coordinate = coordinate;
                            }
                            if (ratio_minor < 0) {
                                ++negative_ratio_minors;
                                if (power >= 2) {
                                    ++negative_ratio_minors_after_base;
                                }
                            }
                            const int ratio_sign =
                                ratio_minor > 0
                                    ? 1
                                    : (ratio_minor < 0 ? -1 : 0);
                            if (
                                ratio_sign != 0
                                && previous_ratio_sign != 0
                                && ratio_sign != previous_ratio_sign
                            ) {
                                ++ratio_sign_changes;
                                ++kernel_ratio_sign_changes;
                                if (
                                    previous_ratio_sign < 0
                                    && ratio_sign > 0
                                ) {
                                    ++negative_to_positive_changes;
                                }
                            }
                            if (ratio_sign != 0) {
                                previous_ratio_sign = ratio_sign;
                            }
                        }
                        previous_even = even_value;
                        previous_odd = odd_value;
                        have_previous = true;
                    }
                }
                maximum_ratio_sign_changes = std::max(
                    maximum_ratio_sign_changes,
                    kernel_ratio_sign_changes
                );
            }
        }

        std::cout
            << "SU2_ORDINARY_ANCHORED_KERNEL_PROBE"
            << " maximum_even_label=" << maximum_label
            << " maximum_half_power=" << maximum_half_power
            << " kernels=" << kernels
            << " coordinates=" << coordinates
            << " negative_coordinates=" << negative
            << " ratio_minors=" << ratio_minors
            << " negative_ratio_minors=" << negative_ratio_minors
            << " negative_ratio_minors_after_base="
            << negative_ratio_minors_after_base
            << " even_log_concavity_minors="
            << even_log_concavity_minors
            << " negative_even_log_concavity_minors="
            << negative_even_log_concavity_minors
            << " ratio_sign_changes=" << ratio_sign_changes
            << " negative_to_positive_changes="
            << negative_to_positive_changes
            << " maximum_ratio_sign_changes="
            << maximum_ratio_sign_changes
            << " minimum=" << minimum
            << " witness={label=" << witness_label
            << " half_power=" << witness_power
            << " coordinate=" << witness_coordinate << '}'
            << " minimum_ratio_minor=" << minimum_ratio_minor
            << " ratio_witness={label=" << ratio_witness_label
            << " half_power=" << ratio_witness_power
            << " coordinate=" << ratio_witness_coordinate << '}'
            << " ratio_result="
            << (
                negative_ratio_minors == 0
                    ? "PASS_BOUNDED_EXACT_DIAGNOSTIC"
                    : "COUNTEREXAMPLE"
            )
            << " result="
            << (
                negative == 0
                    ? "PASS_BOUNDED_EXACT_DIAGNOSTIC"
                    : "COUNTEREXAMPLE"
            )
            << '\n';
        return negative == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
