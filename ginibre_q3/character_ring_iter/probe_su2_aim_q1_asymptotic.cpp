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

Vector triangular_step(const Vector& input) {
    Vector output(input.size() + 4U, 0);
    constexpr int weights[5] = {1, 2, 3, 2, 1};
    for (std::size_t source = 0U; source < input.size(); ++source) {
        for (std::size_t shift = 0U; shift < 5U; ++shift) {
            output[source + shift] += input[source] * weights[shift];
        }
    }
    return output;
}

cpp_int line_at(const Vector& line, int radius, int index) {
    if (index < -radius || index > radius) {
        return 0;
    }
    return line[static_cast<std::size_t>(index + radius)];
}

cpp_int multiplicity(const Vector& line, int radius, int index) {
    if (index < 0) {
        return 0;
    }
    return line_at(line, radius, index)
        - line_at(line, radius, index + 1);
}

cpp_int slope(const Vector& line, int radius, int index) {
    return multiplicity(line, radius, index)
        - multiplicity(line, radius, index + 1);
}

Vector image_suffix(
    const Vector& line,
    int radius,
    int level,
    int shell) {
    const int period = 2 * level + 2;
    Vector suffix(static_cast<std::size_t>(level + 1), 0);
    for (int image = shell;
         image * period <= radius + level + 1;
         ++image) {
        for (int label = 0; label <= level; ++label) {
            suffix[static_cast<std::size_t>(label)]
                += multiplicity(line, radius, image * period + label)
                  - multiplicity(
                        line,
                        radius,
                        (image + 1) * period - label - 1);
        }
    }
    return suffix;
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

struct Failure {
    int length = -1;
    int level = -1;
    int shell = -1;
    int label = -1;
    cpp_int value = 0;
    cpp_int left = 0;
    cpp_int right = 0;
    Vector suffix;
};

void record(
    Failure& failure,
    int length,
    int level,
    int shell,
    int label,
    const cpp_int& value,
    const cpp_int& left,
    const cpp_int& right,
    const Vector& suffix = {}) {
    if (failure.length >= 0) {
        return;
    }
    failure = {
        length,
        level,
        shell,
        label,
        value,
        left,
        right,
        suffix};
}

std::string render_failure(const Failure& failure) {
    return "length=" + std::to_string(failure.length)
        + " level="
        + std::to_string(failure.level < 0 ? -1 : 2 * failure.level)
        + " shell=" + std::to_string(failure.shell)
        + " label=" + std::to_string(failure.label)
        + " value=" + failure.value.convert_to<std::string>()
        + " left=" + failure.left.convert_to<std::string>()
        + " right=" + failure.right.convert_to<std::string>()
        + " suffix=" + render(failure.suffix);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_length = argc >= 2
            ? positive_argument(argv[1], "maximum_length")
            : 1000;
        const int maximum_level = argc >= 3
            ? positive_argument(argv[2], "maximum_half_level")
            : 8;
        if (argc > 3 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: probe_su2_aim_q1_asymptotic "
                "[maximum_length] [maximum_half_level]");
        }

        std::uint64_t lengths = 0U;
        std::uint64_t pair_checks = 0U;
        std::uint64_t pair_failures = 0U;
        std::uint64_t suffix_profiles = 0U;
        std::uint64_t coordinate_checks = 0U;
        std::uint64_t coordinate_failures = 0U;
        std::uint64_t monotonicity_checks = 0U;
        std::uint64_t monotonicity_failures = 0U;
        std::uint64_t generator_payment_checks = 0U;
        std::uint64_t generator_payment_failures = 0U;
        Failure first_pair;
        Failure first_coordinate;
        Failure first_monotonicity;
        Failure first_generator_payment;

        Vector line{cpp_int(1)};
        int radius = 0;
        for (int length = 1; length <= maximum_length; ++length) {
            line = triangular_step(line);
            radius += 2;
            ++lengths;
            for (int level = 2; level <= maximum_level; ++level) {
                const int period = 2 * level + 2;
                for (int image = 1;
                     image * period <= radius + level + 1;
                     ++image) {
                    for (int label = 0; label < level; ++label) {
                        const int left_index = image * period + label;
                        const int right_index
                            = (image + 1) * period - label - 2;
                        const cpp_int left
                            = slope(line, radius, left_index);
                        const cpp_int right
                            = slope(line, radius, right_index);
                        const cpp_int margin = left + right;
                        ++pair_checks;
                        if (margin < 0) {
                            ++pair_failures;
                            if (first_pair.length < 0) {
                                record(
                                    first_pair,
                                    length,
                                    level,
                                    image,
                                    label,
                                    margin,
                                    left,
                                    right,
                                    image_suffix(
                                        line,
                                        radius,
                                        level,
                                        image));
                            }
                        }
                    }
                }

                const int shell_count = radius / period + 1;
                for (int shell = 1; shell < shell_count; ++shell) {
                    const Vector suffix
                        = image_suffix(line, radius, level, shell);
                    ++suffix_profiles;
                    for (int label = 0; label <= level; ++label) {
                        ++coordinate_checks;
                        if (suffix[static_cast<std::size_t>(label)] < 0) {
                            ++coordinate_failures;
                            record(
                                first_coordinate,
                                length,
                                level,
                                shell,
                                label,
                                suffix[static_cast<std::size_t>(label)],
                                0,
                                0,
                                suffix);
                        }
                    }
                    for (int label = 0; label < level; ++label) {
                        const cpp_int margin
                            = suffix[static_cast<std::size_t>(label)]
                              - suffix[
                                  static_cast<std::size_t>(label + 1)];
                        ++monotonicity_checks;
                        if (margin < 0) {
                            ++monotonicity_failures;
                            record(
                                first_monotonicity,
                                length,
                                level,
                                shell,
                                label,
                                margin,
                                suffix[static_cast<std::size_t>(label)],
                                suffix[
                                    static_cast<std::size_t>(label + 1)],
                                suffix);
                        }
                    }
                }
                if (level == 2) {
                    for (int shell = 1;
                         shell * period - 2 <= radius;
                         ++shell) {
                        const Vector reserve
                            = image_suffix(
                                line,
                                radius,
                                level,
                                shell);
                        cpp_int boundary = 0;
                        for (int index = shell * period - 2;
                             index <= shell * period + 1;
                             ++index) {
                            boundary += multiplicity(
                                line,
                                radius,
                                index);
                        }
                        const cpp_int load
                            = 2 * (
                                reserve[1] - reserve[2])
                              + 2 * reserve[2];
                        const cpp_int payment = boundary - load;
                        ++generator_payment_checks;
                        if (payment < 0) {
                            ++generator_payment_failures;
                            record(
                                first_generator_payment,
                                length,
                                level,
                                shell,
                                0,
                                payment,
                                boundary,
                                load,
                                reserve);
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_AIM_Q1_ASYMPTOTIC"
            << " maximum_length=" << maximum_length
            << " maximum_level=" << 2 * maximum_level
            << " lengths=" << lengths
            << " pair_checks=" << pair_checks
            << " pair_failures=" << pair_failures
            << " suffix_profiles=" << suffix_profiles
            << " coordinate_checks=" << coordinate_checks
            << " coordinate_failures=" << coordinate_failures
            << " monotonicity_checks=" << monotonicity_checks
            << " monotonicity_failures=" << monotonicity_failures
            << " generator_payment_checks=" << generator_payment_checks
            << " generator_payment_failures="
            << generator_payment_failures
            << '\n'
            << "FIRST_PAIR_FAILURE "
            << render_failure(first_pair) << '\n'
            << "FIRST_COORDINATE_FAILURE "
            << render_failure(first_coordinate) << '\n'
            << "FIRST_MONOTONICITY_FAILURE "
            << render_failure(first_monotonicity) << '\n'
            << "FIRST_GENERATOR_PAYMENT_FAILURE "
            << render_failure(first_generator_payment) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
