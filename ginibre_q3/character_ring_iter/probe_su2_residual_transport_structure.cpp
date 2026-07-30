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

int positive_argument(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0LL
        || parsed > static_cast<long long>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
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

Vector multiply(const Matrix& matrix, const Vector& vector) {
    Vector output(matrix.size(), 0);
    for (std::size_t row = 0U; row < matrix.size(); ++row) {
        for (std::size_t column = 0U; column < matrix.size(); ++column) {
            output[row] += matrix[row][column] * vector[column];
        }
    }
    return output;
}

bool log_concave(const Vector& vector) {
    bool support_started = false;
    bool support_ended = false;
    for (std::size_t index = 0U; index < vector.size(); ++index) {
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
    return support_started;
}

bool full_cone(const Vector& vector, const std::vector<Matrix>& fusion) {
    if (vector[0] <= 0) {
        return false;
    }
    for (std::size_t radius = 0U; radius < fusion.size(); ++radius) {
        const Vector translated = multiply(fusion[radius], vector);
        for (std::size_t target = 0U; target < vector.size(); ++target) {
            if (vector[0] * translated[target]
                < vector[radius] * vector[target]) {
                return false;
            }
        }
    }
    return true;
}

bool reflected_ratio_monotone(const Vector& vector) {
    const std::size_t level = vector.size() - 1U;
    for (std::size_t index = 0U; index < level; ++index) {
        const cpp_int left
            = vector[level - index] * vector[index + 1U];
        const cpp_int right
            = vector[level - index - 1U] * vector[index];
        if (left > right) {
            return false;
        }
    }
    return true;
}

unsigned reflected_ratio_sign_changes(const Vector& vector) {
    const std::size_t level = vector.size() - 1U;
    int previous_sign = 0;
    unsigned changes = 0U;
    for (std::size_t index = 0U; index < level; ++index) {
        const cpp_int difference
            = vector[level - index - 1U] * vector[index]
              - vector[level - index] * vector[index + 1U];
        const int sign = difference > 0 ? 1 : difference < 0 ? -1 : 0;
        if (sign != 0 && previous_sign != 0 && sign != previous_sign) {
            ++changes;
        }
        if (sign != 0) {
            previous_sign = sign;
        }
    }
    return changes;
}

bool boundary_cone(const Vector& vector) {
    const std::size_t level = vector.size() - 1U;
    for (std::size_t radius = 0U; radius <= level; ++radius) {
        if (vector[0] * vector[level - radius]
            < vector[radius] * vector[level]) {
            return false;
        }
    }
    return true;
}

std::string render(const Vector& vector) {
    std::string text = "[";
    for (std::size_t index = 0U; index < vector.size(); ++index) {
        if (index != 0U) {
            text += ',';
        }
        text += vector[index].convert_to<std::string>();
    }
    return text + ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 6;
        const int maximum_coefficient = argc >= 3
            ? positive_argument(argv[2], "maximum_coefficient")
            : 3;
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_residual_transport_structure "
                "[maximum_half_level] [maximum_coefficient]");
        }

        std::uint64_t log_concave_full_profiles = 0U;
        std::uint64_t nonmonotone_reflected_ratios = 0U;
        std::vector<std::uint64_t> sign_change_counts(
            static_cast<std::size_t>(maximum_level + 1), 0U);
        std::uint64_t transforms = 0U;
        std::uint64_t boundary_failures = 0U;
        Vector first_nonmonotone;
        unsigned first_nonmonotone_changes = 0U;
        Vector first_failure_input;
        Vector first_failure_output;
        int first_failure_factor = -1;

        for (int level = 1; level <= maximum_level; ++level) {
            std::vector<Matrix> fusion;
            for (int factor = 0; factor <= level; ++factor) {
                fusion.push_back(fusion_matrix(level, factor));
            }
            const std::uint64_t base
                = static_cast<std::uint64_t>(maximum_coefficient + 1);
            std::uint64_t count = 1U;
            for (int index = 0; index <= level; ++index) {
                count *= base;
            }
            for (std::uint64_t code = 1U; code < count; ++code) {
                std::uint64_t remainder = code;
                Vector input(static_cast<std::size_t>(level + 1), 0);
                for (int index = 0; index <= level; ++index) {
                    input[static_cast<std::size_t>(index)]
                        = remainder % base;
                    remainder /= base;
                }
                if (!log_concave(input) || !full_cone(input, fusion)) {
                    continue;
                }
                ++log_concave_full_profiles;
                const unsigned changes
                    = reflected_ratio_sign_changes(input);
                ++sign_change_counts[changes];
                if (!reflected_ratio_monotone(input)) {
                    ++nonmonotone_reflected_ratios;
                    if (first_nonmonotone.empty()) {
                        first_nonmonotone = input;
                        first_nonmonotone_changes = changes;
                    }
                }
                for (int factor = 1; factor <= level; ++factor) {
                    ++transforms;
                    Vector output = multiply(
                        fusion[static_cast<std::size_t>(factor)],
                        multiply(
                            fusion[static_cast<std::size_t>(factor)],
                            input));
                    if (boundary_cone(output)) {
                        continue;
                    }
                    ++boundary_failures;
                    if (first_failure_input.empty()) {
                        first_failure_input = input;
                        first_failure_output = output;
                        first_failure_factor = factor;
                    }
                }
            }
        }

        std::cout << "SU2_RESIDUAL_TRANSPORT_STRUCTURE"
                  << " maximum_level=" << 2 * maximum_level
                  << " maximum_coefficient=" << maximum_coefficient
                  << " log_concave_full_profiles="
                  << log_concave_full_profiles
                  << " nonmonotone_reflected_ratios="
                  << nonmonotone_reflected_ratios
                  << " sign_change_counts=[";
        for (std::size_t index = 0U;
             index < sign_change_counts.size();
             ++index) {
            if (index != 0U) {
                std::cout << ',';
            }
            std::cout << sign_change_counts[index];
        }
        std::cout << "] transforms=" << transforms
                  << " boundary_failures=" << boundary_failures
                  << " first_nonmonotone=" << render(first_nonmonotone)
                  << " first_nonmonotone_changes="
                  << first_nonmonotone_changes
                  << " first_failure_factor=" << first_failure_factor
                  << " first_failure_input="
                  << render(first_failure_input)
                  << " first_failure_output="
                  << render(first_failure_output) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
