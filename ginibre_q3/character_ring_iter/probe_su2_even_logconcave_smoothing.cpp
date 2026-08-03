#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Vector = std::vector<Integer>;
using Matrix = std::vector<Vector>;

Matrix fusion_matrix(int level, int factor) {
    const int dimension = level + 1;
    Matrix result(
        static_cast<std::size_t>(dimension),
        Vector(static_cast<std::size_t>(dimension), 0)
    );
    for (int source = 0; source <= level; ++source) {
        const int lower = std::abs(source - factor);
        const int upper = std::min(
            source + factor, 2 * level - source - factor
        );
        for (int target = lower; target <= upper; target += 2) {
            result[static_cast<std::size_t>(target)]
                  [static_cast<std::size_t>(source)] = 1;
        }
    }
    return result;
}

Vector apply_matrix(const Matrix& matrix, const Vector& input) {
    Vector output(matrix.size(), 0);
    for (std::size_t row = 0U; row < matrix.size(); ++row) {
        for (std::size_t column = 0U; column < input.size(); ++column) {
            output[row] += matrix[row][column] * input[column];
        }
    }
    return output;
}

Matrix multiplication_matrix(
    const Vector& profile,
    const std::vector<Matrix>& fusion
) {
    Matrix output(profile.size(), Vector(profile.size(), 0));
    for (std::size_t column = 0U; column < profile.size(); ++column) {
        const Vector image = apply_matrix(fusion[column], profile);
        for (std::size_t row = 0U; row < profile.size(); ++row) {
            output[row][column] = image[row];
        }
    }
    return output;
}

Integer determinant(Matrix matrix) {
    if (matrix.empty()) {
        return 1;
    }
    Integer previous = 1;
    int sign = 1;
    for (std::size_t pivot_index = 0U;
         pivot_index + 1U < matrix.size();
         ++pivot_index) {
        std::size_t pivot_row = pivot_index;
        while (pivot_row < matrix.size()
               && matrix[pivot_row][pivot_index] == 0) {
            ++pivot_row;
        }
        if (pivot_row == matrix.size()) {
            return 0;
        }
        if (pivot_row != pivot_index) {
            std::swap(matrix[pivot_row], matrix[pivot_index]);
            sign = -sign;
        }
        const Integer pivot = matrix[pivot_index][pivot_index];
        for (std::size_t row = pivot_index + 1U; row < matrix.size(); ++row) {
            for (std::size_t column = pivot_index + 1U;
                 column < matrix.size();
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
    return sign * matrix.back().back();
}

struct PrincipalObstruction {
    unsigned order = 0U;
    Integer value = 0;
};

PrincipalObstruction first_negative_principal_minor(const Matrix& matrix) {
    const std::uint64_t masks = UINT64_C(1) << matrix.size();
    for (unsigned order = 1U; order <= matrix.size(); ++order) {
        for (std::uint64_t mask = 1U; mask < masks; ++mask) {
            if (std::popcount(mask) != static_cast<int>(order)) {
                continue;
            }
            std::vector<std::size_t> indices;
            for (std::size_t index = 0U; index < matrix.size(); ++index) {
                if ((mask & (UINT64_C(1) << index)) != 0U) {
                    indices.push_back(index);
                }
            }
            Matrix principal(indices.size(), Vector(indices.size(), 0));
            for (std::size_t row = 0U; row < indices.size(); ++row) {
                for (std::size_t column = 0U;
                     column < indices.size();
                     ++column) {
                    principal[row][column]
                        = matrix[indices[row]][indices[column]];
                }
            }
            const Integer value = determinant(principal);
            if (value < 0) {
                return {order, value};
            }
        }
    }
    return {};
}

Vector boundary_margins(const Vector& profile) {
    const int level = static_cast<int>(profile.size()) - 1;
    Vector margins(profile.size(), 0);
    for (int radius = 0; radius <= level; ++radius) {
        margins[static_cast<std::size_t>(radius)]
            = profile[0U]
                * profile[static_cast<std::size_t>(level - radius)]
            - profile[static_cast<std::size_t>(radius)] * profile.back();
    }
    return margins;
}

bool compressed_log_concave(const Vector& profile) {
    for (std::size_t index = 2U; index + 2U < profile.size(); index += 2U) {
        if (profile[index] * profile[index]
            < profile[index - 2U] * profile[index + 2U]) {
            return false;
        }
    }
    return true;
}

bool all_nonnegative(const Vector& values) {
    for (const Integer& value : values) {
        if (value < 0) {
            return false;
        }
    }
    return true;
}

void print_vector(const Vector& values) {
    std::cout << '[';
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << ']';
}

}  // namespace

int main() {
    constexpr int level = 10;
    const Matrix fundamental = fusion_matrix(level, 1);
    std::vector<Matrix> fusion;
    fusion.reserve(static_cast<std::size_t>(level + 1));
    for (int factor = 0; factor <= level; ++factor) {
        fusion.push_back(fusion_matrix(level, factor));
    }
    const Vector profile{
        100000000, 0, 10000000000LL, 0, 990000000000LL, 0,
        990000000000LL, 0, 980100000000LL, 0, 980100
    };
    const Vector source_margins = boundary_margins(profile);
    const Vector smoothed = apply_matrix(
        fundamental, apply_matrix(fundamental, profile)
    );
    const Vector smoothed_margins = boundary_margins(smoothed);
    const PrincipalObstruction obstruction
        = first_negative_principal_minor(
            multiplication_matrix(profile, fusion)
        );
    if (!compressed_log_concave(profile) || !all_nonnegative(source_margins)
        || smoothed_margins[8U] >= 0 || obstruction.order == 0U) {
        std::cerr << "SU2_EVEN_LOGCONCAVE_SMOOTHING replay mismatch\n";
        return EXIT_FAILURE;
    }
    std::cout << "SU2_EVEN_LOGCONCAVE_SMOOTHING"
              << " source=";
    print_vector(profile);
    std::cout << " source_margins=";
    print_vector(source_margins);
    std::cout << " smoothed=";
    print_vector(smoothed);
    std::cout << " smoothed_margins=";
    print_vector(smoothed_margins);
    std::cout << " psd_obstruction_order=" << obstruction.order
              << " psd_obstruction_value=" << obstruction.value;
    std::cout << " result=PASS_EXACT\n";
    return EXIT_SUCCESS;
}
