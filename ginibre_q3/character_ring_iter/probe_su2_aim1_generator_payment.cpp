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
using SmallVector = std::vector<long long>;
using SmallMatrix = std::vector<SmallVector>;

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

Vector triangular_step(
    const Vector& input,
    int input_radius,
    int factor) {
    const int output_radius = input_radius + 2 * factor;
    Vector output(static_cast<std::size_t>(2 * output_radius + 1), 0);
    for (int source = -input_radius; source <= input_radius; ++source) {
        const cpp_int& coefficient = input[static_cast<std::size_t>(
            source + input_radius)];
        for (int shift = -2 * factor; shift <= 2 * factor; ++shift) {
            output[static_cast<std::size_t>(
                source + shift + output_radius)]
                += coefficient * (2 * factor + 1 - std::abs(shift));
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

Vector reserve_coordinates(
    const Vector& line,
    int radius,
    int level,
    int shell) {
    const int period = 2 * level + 2;
    Vector reserve(static_cast<std::size_t>(level + 1), 0);
    for (int label = 0; label < level; ++label) {
        for (int image = shell;
             image * period + label <= radius;
             ++image) {
            reserve[static_cast<std::size_t>(label)]
                += slope(line, radius, image * period + label)
                  + slope(
                        line,
                        radius,
                        image * period + period - label - 2);
        }
    }
    for (int image = shell;
         image * period + level <= radius;
         ++image) {
        reserve[static_cast<std::size_t>(level)]
            += slope(line, radius, image * period + level);
    }
    return reserve;
}

SmallMatrix fusion_matrix(int level, int factor) {
    SmallMatrix matrix(
        static_cast<std::size_t>(level + 1),
        SmallVector(static_cast<std::size_t>(level + 1), 0));
    for (int source = 0; source <= level; ++source) {
        const int lower = std::abs(source - factor);
        const int upper
            = std::min(source + factor, 2 * level - source - factor);
        for (int target = lower; target <= upper; ++target) {
            matrix[static_cast<std::size_t>(target)]
                  [static_cast<std::size_t>(source)] = 1;
        }
    }
    return matrix;
}

SmallMatrix multiply(
    const SmallMatrix& left,
    const SmallMatrix& right) {
    SmallMatrix output(
        left.size(),
        SmallVector(left.size(), 0));
    for (std::size_t row = 0U; row < left.size(); ++row) {
        for (std::size_t middle = 0U; middle < left.size(); ++middle) {
            for (std::size_t column = 0U;
                 column < left.size();
                 ++column) {
                output[row][column]
                    += left[row][middle] * right[middle][column];
            }
        }
    }
    return output;
}

SmallMatrix reserve_transform(int level, int factor) {
    const SmallMatrix fusion = fusion_matrix(level, factor);
    const SmallMatrix square = multiply(fusion, fusion);
    SmallMatrix transform(
        static_cast<std::size_t>(level + 1),
        SmallVector(static_cast<std::size_t>(level + 1), 0));
    for (int source = 0; source <= level; ++source) {
        SmallVector step(static_cast<std::size_t>(level + 1), 0);
        for (int index = 0; index <= source; ++index) {
            step[static_cast<std::size_t>(index)] = 1;
        }
        SmallVector image(static_cast<std::size_t>(level + 1), 0);
        for (int row = 0; row <= level; ++row) {
            for (int column = 0; column <= level; ++column) {
                image[static_cast<std::size_t>(row)]
                    += square[static_cast<std::size_t>(row)]
                              [static_cast<std::size_t>(column)]
                       * step[static_cast<std::size_t>(column)];
            }
        }
        for (int row = 0; row < level; ++row) {
            transform[static_cast<std::size_t>(row)]
                     [static_cast<std::size_t>(source)]
                = image[static_cast<std::size_t>(row)]
                  - image[static_cast<std::size_t>(row + 1)];
        }
        transform[static_cast<std::size_t>(level)]
                 [static_cast<std::size_t>(source)]
            = image[static_cast<std::size_t>(level)];
    }
    return transform;
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

std::string render(const SmallVector& values) {
    std::string result = "[";
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            result += ',';
        }
        result += std::to_string(values[index]);
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
    int factor = -1;
    int shell = -1;
    int label = -1;
    cpp_int value = 0;
    cpp_int boundary = 0;
    std::vector<int> word;
    Vector reserve;
    SmallVector transform_row;
};

void record(
    Failure& failure,
    int level,
    int factor,
    int shell,
    int label,
    const cpp_int& value,
    const cpp_int& boundary,
    const std::vector<int>& word,
    const Vector& reserve,
    const SmallVector& transform_row) {
    if (failure.level >= 0) {
        return;
    }
    failure = {
        level,
        factor,
        shell,
        label,
        value,
        boundary,
        word,
        reserve,
        transform_row};
}

std::string render_failure(const Failure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " factor=" + std::to_string(failure.factor)
        + " shell=" + std::to_string(failure.shell)
        + " label=" + std::to_string(failure.label)
        + " value=" + failure.value.convert_to<std::string>()
        + " boundary=" + failure.boundary.convert_to<std::string>()
        + " word=" + render(failure.word)
        + " reserve=" + render(failure.reserve)
        + " transform_row=" + render(failure.transform_row);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::uint64_t samples = argc >= 2
            ? positive_argument(argv[1], "samples")
            : UINT64_C(100000);
        const int maximum_level = argc >= 3
            ? static_cast<int>(positive_argument(
                argv[2],
                "maximum_half_level"))
            : 12;
        const int maximum_length = argc >= 4
            ? static_cast<int>(positive_argument(
                argv[3],
                "maximum_word_length"))
            : 30;
        if (argc > 4 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: probe_su2_aim1_generator_payment "
                "[samples] [maximum_half_level] [maximum_word_length]");
        }

        std::uint64_t payment_checks = 0U;
        std::uint64_t payment_failures = 0U;
        std::uint64_t identity_checks = 0U;
        std::uint64_t identity_failures = 0U;
        Failure first_payment;
        Failure first_identity;

        for (std::uint64_t sample = 0U; sample < samples; ++sample) {
            std::uint64_t state
                = sample ^ UINT64_C(0x243f6a8885a308d3);
            const int level = 2 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_level - 1));
            const int length = 1 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_length));
            Vector line{cpp_int(1)};
            int radius = 0;
            std::vector<int> word;
            word.reserve(static_cast<std::size_t>(length));
            for (int position = 0; position < length; ++position) {
                const int factor = 1 + static_cast<int>(
                    splitmix64(state)
                    % static_cast<std::uint64_t>(level / 2));
                word.push_back(factor);
                line = triangular_step(line, radius, factor);
                radius += 2 * factor;
            }

            const int period = 2 * level + 2;
            for (int factor = 1; factor <= level / 2; ++factor) {
                const SmallMatrix transform
                    = reserve_transform(level, factor);
                const Vector next_line
                    = triangular_step(line, radius, factor);
                const int next_radius = radius + 2 * factor;
                for (int shell = 1;
                     shell * period - 2 * factor <= radius;
                     ++shell) {
                    const Vector reserve = reserve_coordinates(
                        line,
                        radius,
                        level,
                        shell);
                    const Vector next_reserve = reserve_coordinates(
                        next_line,
                        next_radius,
                        level,
                        shell);
                    for (int label = 0; label <= level; ++label) {
                        cpp_int boundary = 0;
                        if (label < 2 * factor) {
                            for (int index
                                     = shell * period
                                       - 2 * factor + label;
                                 index
                                     <= shell * period
                                        + 2 * factor - label - 1;
                                 ++index) {
                                boundary += multiplicity(
                                    line,
                                    radius,
                                    index);
                            }
                        }
                        cpp_int exact = boundary;
                        cpp_int payment = boundary;
                        for (int source = 0; source <= level; ++source) {
                            const long long coefficient
                                = transform[
                                    static_cast<std::size_t>(label)]
                                    [static_cast<std::size_t>(source)];
                            exact += coefficient
                                * reserve[static_cast<std::size_t>(source)];
                            if (coefficient < 0) {
                                payment += coefficient
                                    * reserve[
                                        static_cast<std::size_t>(source)];
                            }
                        }
                        ++identity_checks;
                        if (exact
                            != next_reserve[
                                static_cast<std::size_t>(label)]) {
                            ++identity_failures;
                            record(
                                first_identity,
                                level,
                                factor,
                                shell,
                                label,
                                exact
                                  - next_reserve[
                                      static_cast<std::size_t>(label)],
                                boundary,
                                word,
                                reserve,
                                transform[
                                    static_cast<std::size_t>(label)]);
                        }
                        ++payment_checks;
                        if (payment < 0) {
                            ++payment_failures;
                            record(
                                first_payment,
                                level,
                                factor,
                                shell,
                                label,
                                payment,
                                boundary,
                                word,
                                reserve,
                                transform[
                                    static_cast<std::size_t>(label)]);
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_AIM1_GENERATOR_PAYMENT"
            << " samples=" << samples
            << " maximum_level=" << 2 * maximum_level
            << " maximum_word_length=" << maximum_length
            << " identity_checks=" << identity_checks
            << " identity_failures=" << identity_failures
            << " payment_checks=" << payment_checks
            << " payment_failures=" << payment_failures
            << '\n'
            << "FIRST_IDENTITY_FAILURE "
            << render_failure(first_identity) << '\n'
            << "FIRST_PAYMENT_FAILURE "
            << render_failure(first_payment) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
