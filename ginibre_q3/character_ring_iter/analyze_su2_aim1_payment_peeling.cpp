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

using Vector = std::vector<long long>;
using Matrix = std::vector<Vector>;

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

bool selected_slope(int index, int level, int label) {
    if (index < 0) {
        return false;
    }
    const int period = 2 * level + 2;
    if (label == level) {
        return index >= period + level
            && index % period == level;
    }
    if (index < period + label) {
        return false;
    }
    const int residue = index % period;
    return residue == label || residue == period - label - 2;
}

long long reserve_coefficient(int index, int level, int label) {
    return static_cast<long long>(selected_slope(index, level, label))
        - static_cast<long long>(
            selected_slope(index - 1, level, label));
}

Matrix fusion_matrix(int level, int factor) {
    Matrix matrix(
        static_cast<std::size_t>(level + 1),
        Vector(static_cast<std::size_t>(level + 1), 0));
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

Matrix matrix_product(const Matrix& left, const Matrix& right) {
    Matrix output(
        left.size(),
        Vector(left.size(), 0));
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

Matrix reserve_transform(int level, int factor) {
    const Matrix fusion = fusion_matrix(level, factor);
    const Matrix square = matrix_product(fusion, fusion);
    Matrix transform(
        static_cast<std::size_t>(level + 1),
        Vector(static_cast<std::size_t>(level + 1), 0));
    for (int source = 0; source <= level; ++source) {
        for (int row = 0; row <= level; ++row) {
            long long value = 0;
            for (int column = 0; column <= source; ++column) {
                value += square[static_cast<std::size_t>(row)]
                               [static_cast<std::size_t>(column)];
            }
            if (row < level) {
                for (int column = 0; column <= source; ++column) {
                    value -= square[static_cast<std::size_t>(row + 1)]
                                   [static_cast<std::size_t>(column)];
                }
            }
            transform[static_cast<std::size_t>(row)]
                     [static_cast<std::size_t>(source)] = value;
        }
    }
    return transform;
}

long long payment_coefficient(
    int index,
    int level,
    int next_factor,
    int label,
    const Vector& transform_row) {
    const int period = 2 * level + 2;
    long long value = label < 2 * next_factor
            && period - 2 * next_factor + label <= index
            && index <= period + 2 * next_factor - label - 1
        ? 1
        : 0;
    for (int source = 0; source <= level; ++source) {
        const long long coefficient
            = transform_row[static_cast<std::size_t>(source)];
        if (coefficient < 0) {
            value += coefficient
                * reserve_coefficient(index, level, source);
        }
    }
    return value;
}

long long line_payment_coefficient(
    int index,
    int level,
    int next_factor,
    int label,
    const Vector& transform_row) {
    if (index == 0) {
        return 2 * payment_coefficient(
            0,
            level,
            next_factor,
            label,
            transform_row);
    }
    const int absolute = std::abs(index);
    return payment_coefficient(
               absolute,
               level,
               next_factor,
               label,
               transform_row)
        - payment_coefficient(
               absolute - 1,
               level,
               next_factor,
               label,
               transform_row);
}

Vector multiply_payment(
    const Vector& input,
    int old_factor,
    int horizon) {
    Vector output(static_cast<std::size_t>(horizon + 1), 0);
    const int maximum_character = 2 * old_factor;
    for (int source = 0;
         source < static_cast<int>(input.size());
         ++source) {
        const long long coefficient
            = input[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        for (int character = 0;
             character <= maximum_character;
             ++character) {
            const int lower = std::abs(source - character);
            const int upper = std::min(source + character, horizon);
            for (int target = lower; target <= upper; ++target) {
                output[static_cast<std::size_t>(target)]
                    += coefficient;
            }
        }
    }
    return output;
}

struct Failure {
    int level = -1;
    int old_factor = -1;
    int next_factor = -1;
    int label = -1;
    int index = -1;
    long long value = 0;
    Vector transform_row;
};

void record(
    Failure& failure,
    int level,
    int old_factor,
    int next_factor,
    int label,
    int index,
    long long value,
    const Vector& transform_row) {
    if (failure.level >= 0) {
        return;
    }
    failure = {
        level,
        old_factor,
        next_factor,
        label,
        index,
        value,
        transform_row};
}

std::string render(const Vector& values) {
    std::string result = "[";
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            result += ',';
        }
        result += std::to_string(values[index]);
    }
    return result + ']';
}

std::string render_failure(const Failure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " old_factor=" + std::to_string(failure.old_factor)
        + " next_factor=" + std::to_string(failure.next_factor)
        + " label=" + std::to_string(failure.label)
        + " index=" + std::to_string(failure.index)
        + " value=" + std::to_string(failure.value)
        + " transform_row=" + render(failure.transform_row);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 30;
        if (argc > 2 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim1_payment_peeling "
                "[maximum_half_level]");
        }

        std::uint64_t payments = 0U;
        std::uint64_t coefficient_checks = 0U;
        std::uint64_t coefficient_failures = 0U;
        std::uint64_t line_checks = 0U;
        std::uint64_t line_failures = 0U;
        std::uint64_t normalized_checks = 0U;
        std::uint64_t normalized_failures = 0U;
        Failure first;
        Failure first_line;
        Failure first_normalized;

        for (int level = 2; level <= maximum_level; ++level) {
            const int period = 2 * level + 2;
            const int horizon = 6 * period;
            for (int next_factor = 1;
                 next_factor <= level / 2;
                 ++next_factor) {
                const Matrix transform
                    = reserve_transform(level, next_factor);
                for (int label = 0; label <= level; ++label) {
                    Vector payment(
                        static_cast<std::size_t>(horizon + level + 1),
                        0);
                    for (int index = 0;
                         index < static_cast<int>(payment.size());
                         ++index) {
                        payment[static_cast<std::size_t>(index)]
                            = payment_coefficient(
                                index,
                                level,
                                next_factor,
                                label,
                                transform[
                                    static_cast<std::size_t>(label)]);
                    }
                    for (int old_factor = 1;
                         old_factor <= level / 2;
                         ++old_factor) {
                        ++payments;
                        const Vector product = multiply_payment(
                            payment,
                            old_factor,
                            horizon);
                        for (int index = 0; index <= horizon; ++index) {
                            const long long value
                                = product[static_cast<std::size_t>(index)];
                            ++coefficient_checks;
                            if (value < 0) {
                                ++coefficient_failures;
                                record(
                                    first,
                                    level,
                                    old_factor,
                                    next_factor,
                                    label,
                                    index,
                                    value,
                                    transform[
                                        static_cast<std::size_t>(label)]);
                            }
                        }
                        long long normalized_cumulative = 0;
                        for (int index = 0; index <= horizon; ++index) {
                            normalized_cumulative
                                += (2 * index + 1)
                                  * product[
                                      static_cast<std::size_t>(index)];
                            ++normalized_checks;
                            if (normalized_cumulative < 0) {
                                ++normalized_failures;
                                record(
                                    first_normalized,
                                    level,
                                    old_factor,
                                    next_factor,
                                    label,
                                    index,
                                    normalized_cumulative,
                                    transform[
                                        static_cast<std::size_t>(label)]);
                            }
                        }
                        for (int index = 0; index <= horizon; ++index) {
                            long long value = 0;
                            for (int shift = -2 * old_factor;
                                 shift <= 2 * old_factor;
                                 ++shift) {
                                value += (
                                    2 * old_factor + 1 - std::abs(shift))
                                    * line_payment_coefficient(
                                        index - shift,
                                        level,
                                        next_factor,
                                        label,
                                        transform[
                                            static_cast<std::size_t>(
                                                label)]);
                            }
                            ++line_checks;
                            if (value < 0) {
                                ++line_failures;
                                record(
                                    first_line,
                                    level,
                                    old_factor,
                                    next_factor,
                                    label,
                                    index,
                                    value,
                                    transform[
                                        static_cast<std::size_t>(label)]);
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_AIM1_PAYMENT_PEELING"
            << " maximum_level=" << 2 * maximum_level
            << " payments=" << payments
            << " coefficient_checks=" << coefficient_checks
            << " coefficient_failures=" << coefficient_failures
            << " line_checks=" << line_checks
            << " line_failures=" << line_failures
            << " normalized_checks=" << normalized_checks
            << " normalized_failures=" << normalized_failures
            << '\n'
            << "FIRST_COEFFICIENT_FAILURE "
            << render_failure(first) << '\n'
            << "FIRST_LINE_FAILURE "
            << render_failure(first_line) << '\n'
            << "FIRST_NORMALIZED_FAILURE "
            << render_failure(first_normalized) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
