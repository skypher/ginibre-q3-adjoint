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

Matrix fusion_matrix(int half_level, int factor) {
    Matrix matrix(
        static_cast<std::size_t>(half_level + 1),
        Vector(static_cast<std::size_t>(half_level + 1), 0));
    for (int source = 0; source <= half_level; ++source) {
        const int lower = std::abs(source - factor);
        const int upper = std::min(
            source + factor,
            2 * half_level - source - factor);
        for (int target = lower; target <= upper; ++target) {
            matrix[static_cast<std::size_t>(target)]
                  [static_cast<std::size_t>(source)] = 1;
        }
    }
    return matrix;
}

Matrix matrix_product(const Matrix& left, const Matrix& right) {
    Matrix output(left.size(), Vector(left.size(), 0));
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

Matrix reserve_transform(int half_level, int factor) {
    const Matrix fusion = fusion_matrix(half_level, factor);
    const Matrix square = matrix_product(fusion, fusion);
    Matrix transform(
        static_cast<std::size_t>(half_level + 1),
        Vector(static_cast<std::size_t>(half_level + 1), 0));
    for (int row = 0; row <= half_level; ++row) {
        for (int source = 0; source <= half_level; ++source) {
            long long value = 0;
            for (int column = 0; column <= source; ++column) {
                value += square[static_cast<std::size_t>(row)]
                               [static_cast<std::size_t>(column)];
                if (row < half_level) {
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

std::string render(const Vector& values) {
    std::string output = "[";
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            output += ',';
        }
        output += std::to_string(values[index]);
    }
    return output + ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_half_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 30;
        const int display_half_level = argc >= 3
            ? positive_argument(argv[2], "display_half_level")
            : 8;
        if (argc > 3 || maximum_half_level < 2
            || display_half_level > maximum_half_level) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim1_payment_weights "
                "[maximum_half_level] [display_half_level]");
        }

        std::uint64_t rows = 0U;
        std::uint64_t negative_rows = 0U;
        std::uint64_t monotonicity_failures = 0U;
        std::uint64_t support_failures = 0U;
        std::uint64_t stable_dominance_failures = 0U;
        for (int half_level = 2;
             half_level <= maximum_half_level;
             ++half_level) {
            for (int factor = 1;
                 factor <= half_level / 2;
                 ++factor) {
                const Matrix transform
                    = reserve_transform(half_level, factor);
                for (int row = 0; row <= half_level; ++row) {
                    ++rows;
                    Vector load(
                        static_cast<std::size_t>(half_level + 1),
                        0);
                    bool has_negative = false;
                    for (int source = 0;
                         source <= half_level;
                         ++source) {
                        const long long value = transform[
                            static_cast<std::size_t>(row)]
                            [static_cast<std::size_t>(source)];
                        load[static_cast<std::size_t>(source)]
                            = std::max(0LL, -value);
                        has_negative = has_negative || value < 0;
                    }
                    if (has_negative) {
                        ++negative_rows;
                    }
                    if (has_negative && row >= 2 * factor) {
                        ++support_failures;
                    }

                    Vector suffix_weights(
                        static_cast<std::size_t>(half_level + 1),
                        0);
                    suffix_weights[0] = load[0];
                    for (int source = 1;
                         source <= half_level;
                         ++source) {
                        suffix_weights[
                            static_cast<std::size_t>(source)]
                            = load[static_cast<std::size_t>(source)]
                              - load[
                                  static_cast<std::size_t>(source - 1)];
                        if (suffix_weights[
                                static_cast<std::size_t>(source)] < 0) {
                            ++monotonicity_failures;
                        }
                    }
                    for (int source = 0;
                         source <= half_level;
                         ++source) {
                        long long stable_weight = 0;
                        if (row < factor) {
                            if (source >= row + 1
                                && source
                                    <= 2 * factor - row - 1) {
                                stable_weight = 2;
                            } else if (
                                source >= 2 * factor - row
                                && source
                                    <= 2 * factor + row + 1) {
                                stable_weight = 1;
                            }
                        } else if (row < 2 * factor
                                   && source
                                       >= 3 * row
                                          - 2 * factor + 2
                                   && source
                                       <= 2 * factor + row + 1) {
                            stable_weight = 1;
                        }
                        if (suffix_weights[
                                static_cast<std::size_t>(source)]
                            > stable_weight) {
                            ++stable_dominance_failures;
                        }
                    }

                    if (half_level == display_half_level
                        && has_negative) {
                        std::cout
                            << "DETAIL"
                            << " half_level=" << half_level
                            << " factor=" << factor
                            << " row=" << row
                            << " transform="
                            << render(transform[
                                static_cast<std::size_t>(row)])
                            << " margin_load=" << render(load)
                            << " suffix_weights="
                            << render(suffix_weights)
                            << '\n';
                    }
                }
            }
        }
        std::cout
            << "SU2_AIM1_PAYMENT_WEIGHTS"
            << " maximum_half_level=" << maximum_half_level
            << " rows=" << rows
            << " negative_rows=" << negative_rows
            << " monotonicity_failures=" << monotonicity_failures
            << " support_failures=" << support_failures
            << " stable_dominance_failures="
            << stable_dominance_failures
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
