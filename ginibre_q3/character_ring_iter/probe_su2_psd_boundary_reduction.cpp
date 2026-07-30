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
    return support_started;
}

Matrix multiplication_matrix(const Vector& values) {
    const int level = static_cast<int>(values.size()) - 1;
    Matrix matrix(
        values.size(), Vector(values.size(), 0));
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

cpp_int determinant(Matrix matrix) {
    const std::size_t size = matrix.size();
    if (size == 0U) {
        return 1;
    }
    cpp_int previous = 1;
    int sign = 1;
    for (std::size_t pivot_index = 0U;
         pivot_index + 1U < size;
         ++pivot_index) {
        std::size_t pivot_row = pivot_index;
        while (pivot_row < size
               && matrix[pivot_row][pivot_index] == 0) {
            ++pivot_row;
        }
        if (pivot_row == size) {
            return 0;
        }
        if (pivot_row != pivot_index) {
            std::swap(matrix[pivot_row], matrix[pivot_index]);
            sign = -sign;
        }
        const cpp_int pivot = matrix[pivot_index][pivot_index];
        for (std::size_t row = pivot_index + 1U; row < size; ++row) {
            for (std::size_t column = pivot_index + 1U;
                 column < size;
                 ++column) {
                matrix[row][column]
                    = (matrix[row][column] * pivot
                       - matrix[row][pivot_index]
                            * matrix[pivot_index][column])
                      / previous;
            }
        }
        previous = pivot;
    }
    return sign * matrix[size - 1U][size - 1U];
}

bool positive_semidefinite(const Matrix& matrix) {
    const std::size_t size = matrix.size();
    if (size >= 63U) {
        throw std::overflow_error("principal-minor mask overflow");
    }
    const std::uint64_t masks = UINT64_C(1) << size;
    for (std::uint64_t mask = 1U; mask < masks; ++mask) {
        std::vector<std::size_t> indices;
        for (std::size_t index = 0U; index < size; ++index) {
            if ((mask & (UINT64_C(1) << index)) != 0U) {
                indices.push_back(index);
            }
        }
        Matrix principal(
            indices.size(), Vector(indices.size(), 0));
        for (std::size_t row = 0U; row < indices.size(); ++row) {
            for (std::size_t column = 0U;
                 column < indices.size();
                 ++column) {
                principal[row][column]
                    = matrix[indices[row]][indices[column]];
            }
        }
        if (determinant(principal) < 0) {
            return false;
        }
    }
    return true;
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
        const int maximum_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 6;
        const int maximum_coefficient = argc >= 3
            ? positive_argument(argv[2], "maximum_coefficient")
            : 3;
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_psd_boundary_reduction "
                "[maximum_half_level] [maximum_coefficient]");
        }

        std::uint64_t profiles = 0U;
        std::uint64_t shape_boundary_profiles = 0U;
        std::uint64_t negative_current_profiles = 0U;
        std::uint64_t psd_negative_current_profiles = 0U;
        std::uint64_t diagonal_boundary_profiles = 0U;
        std::uint64_t diagonal_boundary_failures = 0U;
        int first_level = -1;
        int first_radius = -1;
        int first_target = -1;
        cpp_int first_value = 0;
        Vector first_profile;
        int first_diagonal_level = -1;
        int first_diagonal_radius = -1;
        int first_diagonal_target = -1;
        cpp_int first_diagonal_value = 0;
        Vector first_diagonal_profile;

        for (int level = 1; level <= maximum_level; ++level) {
            const std::uint64_t base
                = static_cast<std::uint64_t>(maximum_coefficient + 1);
            std::uint64_t count = 1U;
            for (int index = 0; index <= level; ++index) {
                count *= base;
            }
            for (std::uint64_t code = 1U; code < count; ++code) {
                ++profiles;
                std::uint64_t remainder = code;
                Vector values(static_cast<std::size_t>(level + 1), 0);
                for (int index = 0; index <= level; ++index) {
                    values[static_cast<std::size_t>(index)]
                        = remainder % base;
                    remainder /= base;
                }
                if (values[0] <= 0 || !log_concave(values)) {
                    continue;
                }
                bool boundary = true;
                for (int radius = 0; radius <= level; ++radius) {
                    if (values[0]
                            * values[static_cast<std::size_t>(
                                level - radius)]
                        < values[static_cast<std::size_t>(radius)]
                            * values[static_cast<std::size_t>(level)]) {
                        boundary = false;
                        break;
                    }
                }
                if (!boundary) {
                    continue;
                }
                ++shape_boundary_profiles;
                const Matrix matrix = multiplication_matrix(values);
                bool diagonal = true;
                for (int radius = 0; radius <= level; ++radius) {
                    if (values[0]
                            * matrix[static_cast<std::size_t>(radius)]
                                    [static_cast<std::size_t>(radius)]
                        < values[static_cast<std::size_t>(radius)]
                            * values[static_cast<std::size_t>(radius)]) {
                        diagonal = false;
                        break;
                    }
                }
                if (diagonal) {
                    ++diagonal_boundary_profiles;
                }
                int bad_radius = -1;
                int bad_target = -1;
                cpp_int bad_value = 0;
                for (int radius = 0;
                     radius <= level && bad_radius < 0;
                     ++radius) {
                    for (int target = 0; target <= level; ++target) {
                        const cpp_int current
                            = values[0]
                                * matrix[static_cast<std::size_t>(radius)]
                                        [static_cast<std::size_t>(target)]
                            - values[static_cast<std::size_t>(radius)]
                                * values[static_cast<std::size_t>(target)];
                        if (current < 0) {
                            bad_radius = radius;
                            bad_target = target;
                            bad_value = current;
                            break;
                        }
                    }
                }
                if (bad_radius < 0) {
                    continue;
                }
                ++negative_current_profiles;
                if (diagonal) {
                    ++diagonal_boundary_failures;
                    if (first_diagonal_level < 0) {
                        first_diagonal_level = level;
                        first_diagonal_radius = bad_radius;
                        first_diagonal_target = bad_target;
                        first_diagonal_value = bad_value;
                        first_diagonal_profile = values;
                    }
                }
                if (!positive_semidefinite(matrix)) {
                    continue;
                }
                ++psd_negative_current_profiles;
                if (first_level < 0) {
                    first_level = level;
                    first_radius = bad_radius;
                    first_target = bad_target;
                    first_value = bad_value;
                    first_profile = values;
                }
            }
        }

        std::cout
            << "SU2_PSD_BOUNDARY_REDUCTION"
            << " maximum_level=" << 2 * maximum_level
            << " maximum_coefficient=" << maximum_coefficient
            << " profiles=" << profiles
            << " shape_boundary_profiles=" << shape_boundary_profiles
            << " negative_current_profiles=" << negative_current_profiles
            << " psd_negative_current_profiles="
            << psd_negative_current_profiles
            << " diagonal_boundary_profiles="
            << diagonal_boundary_profiles
            << " diagonal_boundary_failures="
            << diagonal_boundary_failures
            << " first_failure=("
            << (first_level < 0 ? -1 : 2 * first_level) << ','
            << (first_radius < 0 ? -1 : 2 * first_radius) << ','
            << (first_target < 0 ? -1 : 2 * first_target) << ','
            << first_value << ')'
            << " first_profile=" << render(first_profile)
            << " first_diagonal_failure=("
            << (first_diagonal_level < 0
                    ? -1
                    : 2 * first_diagonal_level) << ','
            << (first_diagonal_radius < 0
                    ? -1
                    : 2 * first_diagonal_radius) << ','
            << (first_diagonal_target < 0
                    ? -1
                    : 2 * first_diagonal_target) << ','
            << first_diagonal_value << ')'
            << " first_diagonal_profile="
            << render(first_diagonal_profile)
            << " result="
            << (psd_negative_current_profiles == 0U
                    ? "NO_PSD_BOUNDARY_REDUCTION_FAILURE"
                    : "PSD_BOUNDARY_REDUCTION_FAILURE")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
