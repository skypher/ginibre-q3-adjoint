#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cmath>
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

Vector line_triangular_step(
    const Vector& input,
    int input_radius,
    int factor) {
    const int output_radius = input_radius + 2 * factor;
    Vector output(static_cast<std::size_t>(2 * output_radius + 1), 0);
    for (int source = -input_radius; source <= input_radius; ++source) {
        const cpp_int& coefficient = input[static_cast<std::size_t>(
            source + input_radius)];
        if (coefficient == 0) {
            continue;
        }
        for (int shift = -2 * factor; shift <= 2 * factor; ++shift) {
            output[static_cast<std::size_t>(
                source + shift + output_radius)]
                += coefficient * (2 * factor + 1 - std::abs(shift));
        }
    }
    return output;
}

cpp_int line_value(const Vector& values, int radius, int index) {
    if (index < -radius || index > radius) {
        return 0;
    }
    return values[static_cast<std::size_t>(index + radius)];
}

cpp_int line_gradient(const Vector& values, int radius, int index) {
    return line_value(values, radius, index)
        - line_value(values, radius, index + 1);
}

std::string render(const Vector& values) {
    std::string result = "[";
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            result += ',';
        }
        result += values[index].convert_to<std::string>();
    }
    return result + ']';
}

std::string render(const std::vector<int>& values) {
    std::string result = "[";
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            result += ',';
        }
        result += std::to_string(values[index]);
    }
    return result + ']';
}

struct Failure {
    int level = -1;
    int shell_left = -1;
    int shell_right = -1;
    int radius = -1;
    int target = -1;
    cpp_int value = 0;
    std::vector<int> word;
    Vector shell_a;
    Vector shell_b;
};

void record(
    Failure& failure,
    int level,
    int shell_left,
    int shell_right,
    int radius,
    int target,
    const cpp_int& value,
    const std::vector<int>& word,
    const Vector& shell_a,
    const Vector& shell_b) {
    if (failure.level >= 0) {
        return;
    }
    failure = {
        level,
        shell_left,
        shell_right,
        radius,
        target,
        value,
        word,
        shell_a,
        shell_b};
}

std::string render_failure(const Failure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " shells=[" + std::to_string(failure.shell_left)
        + "," + std::to_string(failure.shell_right) + "]"
        + " indices=[" + std::to_string(failure.radius)
        + "," + std::to_string(failure.target) + "]"
        + " value=" + failure.value.convert_to<std::string>()
        + " word=" + render(failure.word)
        + " shell_a=" + render(failure.shell_a)
        + " shell_b=" + render(failure.shell_b);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::uint64_t samples = argc >= 2
            ? positive_argument(argv[1], "samples")
            : UINT64_C(10000);
        const int maximum_level = argc >= 3
            ? static_cast<int>(positive_argument(argv[2], "maximum_level"))
            : 24;
        const int maximum_length = argc >= 4
            ? static_cast<int>(positive_argument(argv[3], "maximum_length"))
            : 20;
        if (argc > 4 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: probe_su2_affine_image_shells "
                "[samples] [maximum_half_level] [maximum_word_length]");
        }

        std::uint64_t shell_profiles = 0U;
        std::uint64_t reconstruction_checks = 0U;
        std::uint64_t reconstruction_failures = 0U;
        std::uint64_t shell_coordinate_checks = 0U;
        std::uint64_t shell_coordinate_failures = 0U;
        std::uint64_t shell_current_checks = 0U;
        std::uint64_t shell_current_failures = 0U;
        std::uint64_t cross_current_checks = 0U;
        std::uint64_t cross_current_failures = 0U;
        std::uint64_t suffix_profiles = 0U;
        std::uint64_t suffix_coordinate_checks = 0U;
        std::uint64_t suffix_coordinate_failures = 0U;
        std::uint64_t suffix_current_checks = 0U;
        std::uint64_t suffix_current_failures = 0U;
        std::uint64_t suffix_increment_checks = 0U;
        std::uint64_t suffix_increment_failures = 0U;
        std::uint64_t suffix_spectral_checks = 0U;
        std::uint64_t suffix_spectral_failures = 0U;
        std::uint64_t proper_suffix_monotonicity_checks = 0U;
        std::uint64_t proper_suffix_monotonicity_failures = 0U;
        std::uint64_t proper_suffix_log_concavity_checks = 0U;
        std::uint64_t proper_suffix_log_concavity_failures = 0U;
        std::uint64_t line_dimension_ratio_checks = 0U;
        std::uint64_t line_dimension_ratio_failures = 0U;
        int first_line_dimension_ratio_level = -1;
        int first_line_dimension_ratio_index = -1;
        cpp_int first_line_dimension_ratio_value = 0;
        long double first_negative_eigenvalue = 0.0L;
        int first_negative_eigenvalue_level = -1;
        int first_negative_eigenvalue_shell = -1;
        int first_negative_eigenvalue_mode = -1;
        Failure first_shell_coordinate;
        Failure first_shell_current;
        Failure first_cross_current;
        Failure first_suffix_coordinate;
        Failure first_suffix_current;
        Failure first_suffix_increment;

        for (std::uint64_t sample = 0U; sample < samples; ++sample) {
            std::uint64_t state
                = sample ^ UINT64_C(0xa4093822299f31d0);
            const int level = 2 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_level - 1));
            const int length = 1 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_length));
            Vector profile(static_cast<std::size_t>(level + 1), 0);
            profile[0] = 1;
            Vector line{cpp_int(1)};
            int line_radius = 0;
            std::vector<int> word;
            word.reserve(static_cast<std::size_t>(length));
            for (int position = 0; position < length; ++position) {
                const int factor = 1 + static_cast<int>(
                    splitmix64(state)
                    % static_cast<std::uint64_t>(level));
                word.push_back(factor);
                profile = fusion_step(profile, level, factor);
                const int reduced = std::min(factor, level - factor);
                line = line_triangular_step(line, line_radius, reduced);
                line_radius += 2 * reduced;
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

            const int period = 2 * level + 2;
            for (int index = 0; index < line_radius; ++index) {
                const cpp_int margin
                    = (2 * index + 3)
                        * line_gradient(line, line_radius, index)
                      - (2 * index + 1)
                        * line_gradient(line, line_radius, index + 1);
                ++line_dimension_ratio_checks;
                if (margin < 0) {
                    ++line_dimension_ratio_failures;
                    if (first_line_dimension_ratio_level < 0) {
                        first_line_dimension_ratio_level = level;
                        first_line_dimension_ratio_index = index;
                        first_line_dimension_ratio_value = margin;
                    }
                }
            }
            const int shell_count = line_radius / period + 1;
            std::vector<Vector> shells;
            std::vector<Matrix> shell_matrices;
            Vector reconstructed(static_cast<std::size_t>(level + 1), 0);
            for (int shell = 0; shell < shell_count; ++shell) {
                Vector values(static_cast<std::size_t>(level + 1), 0);
                for (int radius = 0; radius <= level; ++radius) {
                    values[static_cast<std::size_t>(radius)]
                        = line_gradient(
                            line,
                            line_radius,
                            shell * period + radius)
                          - line_gradient(
                                line,
                                line_radius,
                                (shell + 1) * period - radius - 1);
                    reconstructed[static_cast<std::size_t>(radius)]
                        += values[static_cast<std::size_t>(radius)];
                    ++shell_coordinate_checks;
                    if (values[static_cast<std::size_t>(radius)] < 0) {
                        ++shell_coordinate_failures;
                        record(
                            first_shell_coordinate,
                            level,
                            shell,
                            shell,
                            radius,
                            radius,
                            values[static_cast<std::size_t>(radius)],
                            word,
                            values,
                            values);
                    }
                }
                ++shell_profiles;
                shell_matrices.push_back(multiplication_matrix(values));
                shells.push_back(std::move(values));
            }

            for (int radius = 0; radius <= level; ++radius) {
                ++reconstruction_checks;
                if (reconstructed[static_cast<std::size_t>(radius)]
                    != square[static_cast<std::size_t>(radius)]) {
                    ++reconstruction_failures;
                }
            }

            for (std::size_t shell = 0U; shell < shells.size(); ++shell) {
                const Vector& values = shells[shell];
                const Matrix& matrix = shell_matrices[shell];
                for (int radius = 0; radius <= level; ++radius) {
                    for (int target = 0; target <= level; ++target) {
                        const cpp_int current
                            = values[0]
                                * matrix[static_cast<std::size_t>(radius)]
                                        [static_cast<std::size_t>(target)]
                              - values[static_cast<std::size_t>(radius)]
                                * values[static_cast<std::size_t>(target)];
                        ++shell_current_checks;
                        if (current < 0) {
                            ++shell_current_failures;
                            record(
                                first_shell_current,
                                level,
                                static_cast<int>(shell),
                                static_cast<int>(shell),
                                radius,
                                target,
                                current,
                                word,
                                values,
                                values);
                        }
                    }
                }
            }

            for (std::size_t left = 0U; left < shells.size(); ++left) {
                for (std::size_t right = left + 1U;
                     right < shells.size();
                     ++right) {
                    const Vector& a = shells[left];
                    const Vector& b = shells[right];
                    const Matrix& matrix_a = shell_matrices[left];
                    const Matrix& matrix_b = shell_matrices[right];
                    for (int radius = 0; radius <= level; ++radius) {
                        for (int target = 0; target <= level; ++target) {
                            const cpp_int current
                                = a[0]
                                    * matrix_b[
                                        static_cast<std::size_t>(radius)]
                                        [static_cast<std::size_t>(target)]
                                  + b[0]
                                    * matrix_a[
                                        static_cast<std::size_t>(radius)]
                                        [static_cast<std::size_t>(target)]
                                  - a[static_cast<std::size_t>(radius)]
                                    * b[static_cast<std::size_t>(target)]
                                  - b[static_cast<std::size_t>(radius)]
                                    * a[static_cast<std::size_t>(target)];
                            ++cross_current_checks;
                            if (current < 0) {
                                ++cross_current_failures;
                                record(
                                    first_cross_current,
                                    level,
                                    static_cast<int>(left),
                                    static_cast<int>(right),
                                    radius,
                                    target,
                                    current,
                                    word,
                                    a,
                                    b);
                            }
                        }
                    }
                }
            }

            Vector suffix(static_cast<std::size_t>(level + 1), 0);
            Matrix outer_current(
                static_cast<std::size_t>(level + 1),
                Vector(static_cast<std::size_t>(level + 1), 0));
            for (std::size_t shell = shells.size(); shell-- > 0U;) {
                for (int radius = 0; radius <= level; ++radius) {
                    suffix[static_cast<std::size_t>(radius)]
                        += shells[shell][static_cast<std::size_t>(radius)];
                    ++suffix_coordinate_checks;
                    if (suffix[static_cast<std::size_t>(radius)] < 0) {
                        ++suffix_coordinate_failures;
                        record(
                            first_suffix_coordinate,
                            level,
                            static_cast<int>(shell),
                            static_cast<int>(shells.size() - 1U),
                            radius,
                            radius,
                            suffix[static_cast<std::size_t>(radius)],
                            word,
                            suffix,
                            suffix);
                    }
                }
                ++suffix_profiles;
                if (shell > 0U) {
                    for (int radius = 0; radius < level; ++radius) {
                        ++proper_suffix_monotonicity_checks;
                        if (suffix[static_cast<std::size_t>(radius)]
                            < suffix[
                                static_cast<std::size_t>(radius + 1)]) {
                            ++proper_suffix_monotonicity_failures;
                        }
                    }
                    for (int radius = 1; radius < level; ++radius) {
                        ++proper_suffix_log_concavity_checks;
                        if (suffix[static_cast<std::size_t>(radius)]
                                * suffix[static_cast<std::size_t>(radius)]
                            < suffix[
                                  static_cast<std::size_t>(radius - 1)]
                                * suffix[
                                  static_cast<std::size_t>(radius + 1)]) {
                            ++proper_suffix_log_concavity_failures;
                        }
                    }
                }
                const long double pi = std::acos(-1.0L);
                for (int mode = 1; mode <= level + 1; ++mode) {
                    const long double theta
                        = static_cast<long double>(mode) * pi
                          / static_cast<long double>(2 * level + 2);
                    long double eigenvalue = 0.0L;
                    long double scale = 0.0L;
                    for (int radius = 0; radius <= level; ++radius) {
                        const long double character
                            = std::sin(
                                static_cast<long double>(2 * radius + 1)
                                * theta)
                              / std::sin(theta);
                        const long double coefficient
                            = suffix[static_cast<std::size_t>(radius)]
                                  .convert_to<long double>();
                        eigenvalue += coefficient * character;
                        scale += std::abs(coefficient * character);
                    }
                    ++suffix_spectral_checks;
                    if (eigenvalue
                        < -1.0e-12L * std::max(1.0L, scale)) {
                        ++suffix_spectral_failures;
                        if (first_negative_eigenvalue_level < 0) {
                            first_negative_eigenvalue = eigenvalue;
                            first_negative_eigenvalue_level = level;
                            first_negative_eigenvalue_shell
                                = static_cast<int>(shell);
                            first_negative_eigenvalue_mode = mode;
                        }
                    }
                }
                const Matrix suffix_matrix
                    = multiplication_matrix(suffix);
                Matrix current_matrix(
                    static_cast<std::size_t>(level + 1),
                    Vector(static_cast<std::size_t>(level + 1), 0));
                for (int radius = 0; radius <= level; ++radius) {
                    for (int target = 0; target <= level; ++target) {
                        const cpp_int current
                            = suffix[0]
                                * suffix_matrix[
                                    static_cast<std::size_t>(radius)]
                                    [static_cast<std::size_t>(target)]
                              - suffix[static_cast<std::size_t>(radius)]
                                * suffix[static_cast<std::size_t>(target)];
                        current_matrix[
                            static_cast<std::size_t>(radius)]
                            [static_cast<std::size_t>(target)] = current;
                        ++suffix_current_checks;
                        if (current < 0) {
                            ++suffix_current_failures;
                            record(
                                first_suffix_current,
                                level,
                                static_cast<int>(shell),
                                static_cast<int>(shells.size() - 1U),
                                radius,
                                target,
                                current,
                                word,
                                suffix,
                                suffix);
                        }
                        const cpp_int increment
                            = current
                              - outer_current[
                                  static_cast<std::size_t>(radius)]
                                  [static_cast<std::size_t>(target)];
                        ++suffix_increment_checks;
                        if (increment < 0) {
                            ++suffix_increment_failures;
                            record(
                                first_suffix_increment,
                                level,
                                static_cast<int>(shell),
                                static_cast<int>(shells.size() - 1U),
                                radius,
                                target,
                                increment,
                                word,
                                suffix,
                                shells[shell]);
                        }
                    }
                }
                outer_current = std::move(current_matrix);
            }
        }

        std::cout
            << "SU2_AFFINE_IMAGE_SHELLS"
            << " samples=" << samples
            << " maximum_level=" << 2 * maximum_level
            << " maximum_word_length=" << maximum_length
            << " shell_profiles=" << shell_profiles
            << " reconstruction_checks=" << reconstruction_checks
            << " reconstruction_failures=" << reconstruction_failures
            << " shell_coordinate_checks=" << shell_coordinate_checks
            << " shell_coordinate_failures=" << shell_coordinate_failures
            << " shell_current_checks=" << shell_current_checks
            << " shell_current_failures=" << shell_current_failures
            << " cross_current_checks=" << cross_current_checks
            << " cross_current_failures=" << cross_current_failures
            << " suffix_profiles=" << suffix_profiles
            << " suffix_coordinate_checks=" << suffix_coordinate_checks
            << " suffix_coordinate_failures="
            << suffix_coordinate_failures
            << " suffix_current_checks=" << suffix_current_checks
            << " suffix_current_failures=" << suffix_current_failures
            << " suffix_increment_checks=" << suffix_increment_checks
            << " suffix_increment_failures="
            << suffix_increment_failures
            << " suffix_spectral_checks=" << suffix_spectral_checks
            << " suffix_spectral_failures=" << suffix_spectral_failures
            << " proper_suffix_monotonicity_checks="
            << proper_suffix_monotonicity_checks
            << " proper_suffix_monotonicity_failures="
            << proper_suffix_monotonicity_failures
            << " proper_suffix_log_concavity_checks="
            << proper_suffix_log_concavity_checks
            << " proper_suffix_log_concavity_failures="
            << proper_suffix_log_concavity_failures
            << " line_dimension_ratio_checks="
            << line_dimension_ratio_checks
            << " line_dimension_ratio_failures="
            << line_dimension_ratio_failures
            << '\n'
            << "FIRST_SHELL_COORDINATE_FAILURE "
            << render_failure(first_shell_coordinate) << '\n'
            << "FIRST_SHELL_CURRENT_FAILURE "
            << render_failure(first_shell_current) << '\n'
            << "FIRST_CROSS_CURRENT_FAILURE "
            << render_failure(first_cross_current) << '\n'
            << "FIRST_SUFFIX_COORDINATE_FAILURE "
            << render_failure(first_suffix_coordinate) << '\n'
            << "FIRST_SUFFIX_CURRENT_FAILURE "
            << render_failure(first_suffix_current) << '\n'
            << "FIRST_SUFFIX_INCREMENT_FAILURE "
            << render_failure(first_suffix_increment) << '\n'
            << "FIRST_SUFFIX_SPECTRAL_FAILURE"
            << " level="
            << (first_negative_eigenvalue_level < 0
                    ? -1
                    : 2 * first_negative_eigenvalue_level)
            << " shell=" << first_negative_eigenvalue_shell
            << " mode=" << first_negative_eigenvalue_mode
            << " eigenvalue="
            << static_cast<double>(first_negative_eigenvalue)
            << '\n'
            << "FIRST_LINE_DIMENSION_RATIO_FAILURE"
            << " level="
            << (first_line_dimension_ratio_level < 0
                    ? -1
                    : 2 * first_line_dimension_ratio_level)
            << " index=" << first_line_dimension_ratio_index
            << " value=" << first_line_dimension_ratio_value
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
