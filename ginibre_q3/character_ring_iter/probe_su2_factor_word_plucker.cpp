#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
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
using Matrix = std::vector<Vector>;

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

Vector fusion_step(const Vector& input, int level, int factor) {
    Vector output(static_cast<std::size_t>(level + 1), 0);
    for (int source = 0; source <= level; ++source) {
        const cpp_int& coefficient
            = input[static_cast<std::size_t>(source)];
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

Matrix multiplication_matrix(const Vector& values) {
    const int level = static_cast<int>(values.size()) - 1;
    Matrix matrix(values.size(), Vector(values.size(), 0));
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            const int lower = std::abs(left - right);
            const int upper
                = std::min(left + right, 2 * level - left - right);
            for (int label = lower; label <= upper; ++label) {
                matrix[static_cast<std::size_t>(left)]
                      [static_cast<std::size_t>(right)]
                    += values[static_cast<std::size_t>(label)];
            }
        }
    }
    return matrix;
}

Matrix factor_step(const Matrix& input, int level, int factor) {
    Matrix output(
        static_cast<std::size_t>(level + 1),
        Vector(static_cast<std::size_t>(level + 1), 0));
    for (int column = 0; column <= level; ++column) {
        Vector input_column(static_cast<std::size_t>(level + 1), 0);
        for (int row = 0; row <= level; ++row) {
            input_column[static_cast<std::size_t>(row)]
                = input[static_cast<std::size_t>(row)]
                       [static_cast<std::size_t>(column)];
        }
        const Vector output_column
            = fusion_step(input_column, level, factor);
        for (int row = 0; row <= level; ++row) {
            output[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(column)]
                = output_column[static_cast<std::size_t>(row)];
        }
    }
    return output;
}

Vector triangular_cyclic_step(
    const Vector& input,
    int level,
    int factor) {
    const int period = 2 * level + 2;
    const int reduced = std::min(factor, level - factor);
    Vector output(static_cast<std::size_t>(period), 0);
    for (int source = 0; source < period; ++source) {
        if (input[static_cast<std::size_t>(source)] == 0) {
            continue;
        }
        for (int shift = -2 * reduced;
             shift <= 2 * reduced;
             ++shift) {
            int target = (source + shift) % period;
            if (target < 0) {
                target += period;
            }
            output[static_cast<std::size_t>(target)]
                += input[static_cast<std::size_t>(source)]
                  * (2 * reduced + 1 - std::abs(shift));
        }
    }
    return output;
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

std::string render(const std::vector<int>& values) {
    std::string text = "[";
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            text += ',';
        }
        text += std::to_string(values[index]);
    }
    return text + ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::uint64_t samples = argc >= 2
            ? positive_argument(argv[1], "samples")
            : UINT64_C(100000);
        const int maximum_level = argc >= 3
            ? static_cast<int>(positive_argument(argv[2], "maximum_level"))
            : 24;
        const int maximum_length = argc >= 4
            ? static_cast<int>(positive_argument(argv[3], "maximum_length"))
            : 20;
        if (argc > 4 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: probe_su2_factor_word_plucker "
                "[samples] [maximum_half_level] [maximum_word_length]");
        }

        std::uint64_t column_pairs = 0U;
        std::uint64_t signed_minors = 0U;
        std::uint64_t sign_conflicts = 0U;
        std::uint64_t current_coordinates = 0U;
        std::uint64_t current_failures = 0U;
        std::uint64_t identity_checks = 0U;
        std::uint64_t identity_failures = 0U;
        std::uint64_t word_matrix_checks = 0U;
        std::uint64_t word_matrix_failures = 0U;
        std::uint64_t semigroup_checks = 0U;
        std::uint64_t semigroup_failures = 0U;
        std::uint64_t cyclic_kernel_checks = 0U;
        std::uint64_t cyclic_kernel_failures = 0U;
        std::uint64_t cyclic_gradient_checks = 0U;
        std::uint64_t cyclic_gradient_failures = 0U;
        std::uint64_t cyclic_current_checks = 0U;
        std::uint64_t cyclic_current_failures = 0U;
        std::uint64_t gap_group_checks = 0U;
        std::uint64_t gap_group_failures = 0U;
        std::uint64_t sum_group_checks = 0U;
        std::uint64_t sum_group_failures = 0U;
        std::uint64_t cyclic_gap_group_checks = 0U;
        std::uint64_t cyclic_gap_group_failures = 0U;
        std::uint64_t cyclic_gap_prefix_checks = 0U;
        std::uint64_t cyclic_gap_prefix_failures = 0U;
        std::uint64_t cyclic_gap_suffix_checks = 0U;
        std::uint64_t cyclic_gap_suffix_failures = 0U;
        int first_level = -1;
        int first_left_column = -1;
        int first_right_column = -1;
        int first_positive_radius = -1;
        int first_negative_radius = -1;
        Vector first_profile;
        std::vector<int> first_word;
        int first_gap_failure_level = -1;
        int first_gap_failure_radius = -1;
        int first_gap_failure_target = -1;
        int first_gap_failure_group = -1;
        cpp_int first_gap_failure_value = 0;
        Vector first_gap_failure_profile;
        std::vector<int> first_gap_failure_word;
        int first_cyclic_gap_failure_level = -1;
        int first_cyclic_gap_failure_radius = -1;
        int first_cyclic_gap_failure_target = -1;
        int first_cyclic_gap_failure_group = -1;
        cpp_int first_cyclic_gap_failure_value = 0;
        Vector first_cyclic_gap_failure_profile;
        std::vector<int> first_cyclic_gap_failure_word;
        int first_cyclic_prefix_failure_level = -1;
        int first_cyclic_prefix_failure_radius = -1;
        int first_cyclic_prefix_failure_target = -1;
        int first_cyclic_prefix_failure_group = -1;
        cpp_int first_cyclic_prefix_failure_value = 0;
        Vector first_cyclic_prefix_failure_profile;
        std::vector<int> first_cyclic_prefix_failure_word;

        for (std::uint64_t sample = 0U; sample < samples; ++sample) {
            std::uint64_t state
                = sample ^ UINT64_C(0x243f6a8885a308d3);
            const int level = 2 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_level - 1));
            const int length = 1 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_length));
            const bool audit_semigroup = sample < UINT64_C(1000);
            Vector profile(static_cast<std::size_t>(level + 1), 0);
            profile[0] = 1;
            Matrix word_matrix;
            if (audit_semigroup) {
                word_matrix.assign(
                    static_cast<std::size_t>(level + 1),
                    Vector(static_cast<std::size_t>(level + 1), 0));
                for (int index = 0; index <= level; ++index) {
                    word_matrix[static_cast<std::size_t>(index)]
                               [static_cast<std::size_t>(index)] = 1;
                }
            }
            Vector cyclic_profile;
            if (audit_semigroup) {
                cyclic_profile.assign(
                    static_cast<std::size_t>(2 * level + 2), 0);
                cyclic_profile[0] = 1;
            }
            std::vector<int> word;
            word.reserve(static_cast<std::size_t>(length));
            for (int position = 0; position < length; ++position) {
                const int factor = 1 + static_cast<int>(
                    splitmix64(state)
                    % static_cast<std::uint64_t>(level));
                word.push_back(factor);
                profile = fusion_step(profile, level, factor);
                if (audit_semigroup) {
                    word_matrix = factor_step(word_matrix, level, factor);
                    cyclic_profile = triangular_cyclic_step(
                        cyclic_profile,
                        level,
                        factor);
                }
            }
            const Matrix matrix = multiplication_matrix(profile);
            if (audit_semigroup) {
                for (int row = 0; row <= level; ++row) {
                    for (int column = 0; column <= level; ++column) {
                        ++word_matrix_checks;
                        if (word_matrix[static_cast<std::size_t>(row)]
                                       [static_cast<std::size_t>(column)]
                            != matrix[static_cast<std::size_t>(row)]
                                     [static_cast<std::size_t>(column)]) {
                            ++word_matrix_failures;
                        }
                    }
                }
            }

            for (int left = 0; left < level; ++left) {
                for (int right = left + 1; right <= level; ++right) {
                    ++column_pairs;
                    int positive_radius = -1;
                    int negative_radius = -1;
                    for (int radius = 0; radius <= level; ++radius) {
                        const cpp_int minor
                            = profile[static_cast<std::size_t>(left)]
                                * matrix[
                                    static_cast<std::size_t>(radius)]
                                    [static_cast<std::size_t>(right)]
                              - profile[static_cast<std::size_t>(right)]
                                * matrix[
                                    static_cast<std::size_t>(radius)]
                                    [static_cast<std::size_t>(left)];
                        if (minor > 0) {
                            ++signed_minors;
                            positive_radius = radius;
                        } else if (minor < 0) {
                            ++signed_minors;
                            negative_radius = radius;
                        }
                    }
                    if (positive_radius >= 0 && negative_radius >= 0) {
                        ++sign_conflicts;
                        if (first_level < 0) {
                            first_level = level;
                            first_left_column = left;
                            first_right_column = right;
                            first_positive_radius = positive_radius;
                            first_negative_radius = negative_radius;
                            first_profile = profile;
                            first_word = word;
                        }
                    }
                }
            }

            Vector square(static_cast<std::size_t>(level + 1), 0);
            for (int factor = 0; factor <= level; ++factor) {
                const Vector translated
                    = fusion_step(profile, level, factor);
                for (int target = 0; target <= level; ++target) {
                    square[static_cast<std::size_t>(target)]
                        += profile[static_cast<std::size_t>(factor)]
                          * translated[static_cast<std::size_t>(target)];
                }
            }
            const Matrix square_matrix = multiplication_matrix(square);
            if (audit_semigroup) {
                const int period = 2 * level + 2;
                for (int target = 0; target <= level; ++target) {
                    ++cyclic_gradient_checks;
                    if (cyclic_profile[static_cast<std::size_t>(target)]
                            - cyclic_profile[
                                static_cast<std::size_t>(target + 1)]
                        != square[static_cast<std::size_t>(target)]) {
                        ++cyclic_gradient_failures;
                    }
                    for (int radius = 0; radius <= level; ++radius) {
                        int difference = target - radius;
                        if (difference < 0) {
                            difference += period;
                        }
                        const int reflected = target + radius + 1;
                        ++cyclic_kernel_checks;
                        if (cyclic_profile[
                                static_cast<std::size_t>(difference)]
                                - cyclic_profile[
                                    static_cast<std::size_t>(reflected)]
                            != square_matrix[
                                static_cast<std::size_t>(target)]
                                [static_cast<std::size_t>(radius)]) {
                            ++cyclic_kernel_failures;
                        }
                    }
                }
            }
            for (int radius = 0; radius <= level; ++radius) {
                for (int target = 0; target <= level; ++target) {
                    ++current_coordinates;
                    const cpp_int current
                        = square[0]
                            * square_matrix[
                                static_cast<std::size_t>(radius)]
                                [static_cast<std::size_t>(target)]
                          - square[static_cast<std::size_t>(radius)]
                            * square[static_cast<std::size_t>(target)];
                    if (current < 0) {
                        ++current_failures;
                    }
                    cpp_int plucker_sum = 0;
                    cpp_int semigroup_sum = 0;
                    Vector gap_payments(
                        static_cast<std::size_t>(level + 1), 0);
                    Vector sum_payments(
                        static_cast<std::size_t>(2 * level + 1), 0);
                    Vector cyclic_gap_payments(
                        static_cast<std::size_t>((level + 1) / 2 + 1), 0);
                    for (int left = 0; left < level; ++left) {
                        for (int right = left + 1;
                             right <= level;
                             ++right) {
                            const cpp_int radius_minor
                                = profile[static_cast<std::size_t>(left)]
                                    * matrix[
                                        static_cast<std::size_t>(radius)]
                                        [static_cast<std::size_t>(right)]
                                  - profile[
                                        static_cast<std::size_t>(right)]
                                    * matrix[
                                        static_cast<std::size_t>(radius)]
                                        [static_cast<std::size_t>(left)];
                            const cpp_int target_minor
                                = profile[static_cast<std::size_t>(left)]
                                    * matrix[
                                        static_cast<std::size_t>(target)]
                                        [static_cast<std::size_t>(right)]
                                  - profile[
                                        static_cast<std::size_t>(right)]
                                    * matrix[
                                        static_cast<std::size_t>(target)]
                                        [static_cast<std::size_t>(left)];
                            plucker_sum
                                += radius_minor * target_minor;
                            if (audit_semigroup) {
                                const cpp_int word_radius_minor
                                    = word_matrix[
                                          static_cast<std::size_t>(left)][0]
                                        * word_matrix[
                                            static_cast<std::size_t>(right)]
                                            [static_cast<std::size_t>(radius)]
                                      - word_matrix[
                                            static_cast<std::size_t>(right)][0]
                                        * word_matrix[
                                            static_cast<std::size_t>(left)]
                                            [static_cast<std::size_t>(radius)];
                                const cpp_int word_target_minor
                                    = word_matrix[
                                          static_cast<std::size_t>(left)][0]
                                        * word_matrix[
                                            static_cast<std::size_t>(right)]
                                            [static_cast<std::size_t>(target)]
                                      - word_matrix[
                                            static_cast<std::size_t>(right)][0]
                                        * word_matrix[
                                            static_cast<std::size_t>(left)]
                                            [static_cast<std::size_t>(target)];
                                semigroup_sum
                                    += word_radius_minor * word_target_minor;
                            }
                            gap_payments[static_cast<std::size_t>(
                                right - left)]
                                += radius_minor * target_minor;
                            sum_payments[static_cast<std::size_t>(
                                right + left)]
                                += radius_minor * target_minor;
                            const int gap = right - left;
                            const int cyclic_gap
                                = std::min(gap, level + 1 - gap);
                            cyclic_gap_payments[
                                static_cast<std::size_t>(cyclic_gap)]
                                += radius_minor * target_minor;
                        }
                    }
                    ++identity_checks;
                    if (plucker_sum != current) {
                        ++identity_failures;
                    }
                    if (audit_semigroup) {
                        ++semigroup_checks;
                        if (semigroup_sum != current) {
                            ++semigroup_failures;
                        }
                        const int period = 2 * level + 2;
                        cpp_int cyclic_interval = 0;
                        for (int shift = -radius;
                             shift <= radius;
                             ++shift) {
                            int index = (target + shift) % period;
                            if (index < 0) {
                                index += period;
                            }
                            const int next = (index + 1) % period;
                            cyclic_interval
                                += cyclic_profile[
                                       static_cast<std::size_t>(index)]
                                 - cyclic_profile[
                                       static_cast<std::size_t>(next)];
                        }
                        const cpp_int cyclic_current
                            = square[0] * cyclic_interval
                              - square[static_cast<std::size_t>(radius)]
                                * square[static_cast<std::size_t>(target)];
                        ++cyclic_current_checks;
                        if (cyclic_current != current) {
                            ++cyclic_current_failures;
                        }
                    }
                    for (int gap = 1; gap <= level; ++gap) {
                        ++gap_group_checks;
                        const cpp_int& payment
                            = gap_payments[static_cast<std::size_t>(gap)];
                        if (payment < 0) {
                            ++gap_group_failures;
                            if (first_gap_failure_level < 0) {
                                first_gap_failure_level = level;
                                first_gap_failure_radius = radius;
                                first_gap_failure_target = target;
                                first_gap_failure_group = gap;
                                first_gap_failure_value = payment;
                                first_gap_failure_profile = profile;
                                first_gap_failure_word = word;
                            }
                        }
                    }
                    for (int sum = 1; sum < 2 * level; ++sum) {
                        ++sum_group_checks;
                        if (sum_payments[static_cast<std::size_t>(sum)]
                            < 0) {
                            ++sum_group_failures;
                        }
                    }
                    for (std::size_t cyclic_gap = 1U;
                         cyclic_gap < cyclic_gap_payments.size();
                         ++cyclic_gap) {
                        ++cyclic_gap_group_checks;
                        const cpp_int& payment
                            = cyclic_gap_payments[cyclic_gap];
                        if (payment < 0) {
                            ++cyclic_gap_group_failures;
                            if (first_cyclic_gap_failure_level < 0) {
                                first_cyclic_gap_failure_level = level;
                                first_cyclic_gap_failure_radius = radius;
                                first_cyclic_gap_failure_target = target;
                                first_cyclic_gap_failure_group
                                    = static_cast<int>(cyclic_gap);
                                first_cyclic_gap_failure_value = payment;
                                first_cyclic_gap_failure_profile = profile;
                                first_cyclic_gap_failure_word = word;
                            }
                        }
                    }
                    cpp_int cyclic_prefix = 0;
                    for (std::size_t cyclic_gap = 1U;
                         cyclic_gap < cyclic_gap_payments.size();
                         ++cyclic_gap) {
                        cyclic_prefix += cyclic_gap_payments[cyclic_gap];
                        ++cyclic_gap_prefix_checks;
                        if (cyclic_prefix < 0) {
                            ++cyclic_gap_prefix_failures;
                            if (first_cyclic_prefix_failure_level < 0) {
                                first_cyclic_prefix_failure_level = level;
                                first_cyclic_prefix_failure_radius = radius;
                                first_cyclic_prefix_failure_target = target;
                                first_cyclic_prefix_failure_group
                                    = static_cast<int>(cyclic_gap);
                                first_cyclic_prefix_failure_value
                                    = cyclic_prefix;
                                first_cyclic_prefix_failure_profile
                                    = profile;
                                first_cyclic_prefix_failure_word = word;
                            }
                        }
                    }
                    cpp_int cyclic_suffix = 0;
                    for (std::size_t cyclic_gap
                            = cyclic_gap_payments.size();
                         cyclic_gap-- > 1U;) {
                        cyclic_suffix += cyclic_gap_payments[cyclic_gap];
                        ++cyclic_gap_suffix_checks;
                        if (cyclic_suffix < 0) {
                            ++cyclic_gap_suffix_failures;
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_FACTOR_WORD_PLUCKER"
            << " samples=" << samples
            << " maximum_level=" << 2 * maximum_level
            << " maximum_word_length=" << maximum_length
            << " column_pairs=" << column_pairs
            << " signed_minors=" << signed_minors
            << " sign_conflicts=" << sign_conflicts
            << " current_coordinates=" << current_coordinates
            << " current_failures=" << current_failures
            << " identity_checks=" << identity_checks
            << " identity_failures=" << identity_failures
            << " word_matrix_checks=" << word_matrix_checks
            << " word_matrix_failures=" << word_matrix_failures
            << " semigroup_checks=" << semigroup_checks
            << " semigroup_failures=" << semigroup_failures
            << " cyclic_kernel_checks=" << cyclic_kernel_checks
            << " cyclic_kernel_failures=" << cyclic_kernel_failures
            << " cyclic_gradient_checks=" << cyclic_gradient_checks
            << " cyclic_gradient_failures=" << cyclic_gradient_failures
            << " cyclic_current_checks=" << cyclic_current_checks
            << " cyclic_current_failures=" << cyclic_current_failures
            << " gap_group_checks=" << gap_group_checks
            << " gap_group_failures=" << gap_group_failures
            << " sum_group_checks=" << sum_group_checks
            << " sum_group_failures=" << sum_group_failures
            << " cyclic_gap_group_checks="
            << cyclic_gap_group_checks
            << " cyclic_gap_group_failures="
            << cyclic_gap_group_failures
            << " cyclic_gap_prefix_checks="
            << cyclic_gap_prefix_checks
            << " cyclic_gap_prefix_failures="
            << cyclic_gap_prefix_failures
            << " cyclic_gap_suffix_checks="
            << cyclic_gap_suffix_checks
            << " cyclic_gap_suffix_failures="
            << cyclic_gap_suffix_failures
            << " first_conflict_level="
            << (first_level < 0 ? -1 : 2 * first_level)
            << " first_columns=[" << first_left_column << ','
            << first_right_column << ']'
            << " first_radii=[" << first_positive_radius << ','
            << first_negative_radius << ']'
            << " first_word=" << render(first_word)
            << " first_profile=" << render(first_profile)
            << " first_gap_failure_level="
            << (first_gap_failure_level < 0
                    ? -1
                    : 2 * first_gap_failure_level)
            << " first_gap_failure_indices=["
            << first_gap_failure_radius << ','
            << first_gap_failure_target << ','
            << first_gap_failure_group << ']'
            << " first_gap_failure_value="
            << first_gap_failure_value
            << " first_gap_failure_word="
            << render(first_gap_failure_word)
            << " first_gap_failure_profile="
            << render(first_gap_failure_profile)
            << " first_cyclic_gap_failure_level="
            << (first_cyclic_gap_failure_level < 0
                    ? -1
                    : 2 * first_cyclic_gap_failure_level)
            << " first_cyclic_gap_failure_indices=["
            << first_cyclic_gap_failure_radius << ','
            << first_cyclic_gap_failure_target << ','
            << first_cyclic_gap_failure_group << ']'
            << " first_cyclic_gap_failure_value="
            << first_cyclic_gap_failure_value
            << " first_cyclic_gap_failure_word="
            << render(first_cyclic_gap_failure_word)
            << " first_cyclic_gap_failure_profile="
            << render(first_cyclic_gap_failure_profile)
            << " first_cyclic_prefix_failure_level="
            << (first_cyclic_prefix_failure_level < 0
                    ? -1
                    : 2 * first_cyclic_prefix_failure_level)
            << " first_cyclic_prefix_failure_indices=["
            << first_cyclic_prefix_failure_radius << ','
            << first_cyclic_prefix_failure_target << ','
            << first_cyclic_prefix_failure_group << ']'
            << " first_cyclic_prefix_failure_value="
            << first_cyclic_prefix_failure_value
            << " first_cyclic_prefix_failure_word="
            << render(first_cyclic_prefix_failure_word)
            << " first_cyclic_prefix_failure_profile="
            << render(first_cyclic_prefix_failure_profile)
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
