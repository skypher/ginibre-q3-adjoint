#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_rational;

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

Vector fusion_step(const Vector& values, int level, int factor) {
    Vector output(values.size(), 0);
    for (int source = 0; source <= level; ++source) {
        const cpp_int& coefficient
            = values[static_cast<std::size_t>(source)];
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

cpp_int determinant(Matrix matrix) {
    if (matrix.empty()) {
        return 1;
    }
    cpp_int previous = 1;
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
        const cpp_int pivot = matrix[pivot_index][pivot_index];
        for (std::size_t row = pivot_index + 1U;
             row < matrix.size();
             ++row) {
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
    cpp_int value = 0;
    std::vector<std::size_t> indices;
};

PrincipalObstruction first_negative_principal_minor(const Matrix& matrix) {
    if (matrix.size() >= 63U) {
        throw std::overflow_error("principal-minor mask overflow");
    }
    const std::uint64_t limit = UINT64_C(1) << matrix.size();
    for (unsigned order = 1U; order <= matrix.size(); ++order) {
        for (std::uint64_t mask = 1U; mask < limit; ++mask) {
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
            const cpp_int value = determinant(principal);
            if (value < 0) {
                return {order, value, indices};
            }
        }
    }
    return {};
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

bool has_nonnegative_natural_ldl(const Matrix& matrix) {
    if (matrix.size() <= 1U) {
        return true;
    }
    const std::size_t dimension = matrix.size() - 1U;
    std::vector<std::vector<cpp_rational>> lower(
        dimension,
        std::vector<cpp_rational>(dimension, 0));
    std::vector<cpp_rational> diagonal(dimension, 0);
    for (std::size_t column = 0U; column < dimension; ++column) {
        lower[column][column] = 1;
        cpp_rational pivot
            = matrix[column + 1U][column + 1U];
        for (std::size_t previous = 0U;
             previous < column;
             ++previous) {
            pivot -= lower[column][previous]
                * lower[column][previous] * diagonal[previous];
        }
        if (pivot < 0) {
            return false;
        }
        diagonal[column] = pivot;
        for (std::size_t row = column + 1U;
             row < dimension;
             ++row) {
            cpp_rational residual
                = matrix[row + 1U][column + 1U];
            for (std::size_t previous = 0U;
                 previous < column;
                 ++previous) {
                residual -= lower[row][previous]
                    * lower[column][previous] * diagonal[previous];
            }
            if (pivot == 0) {
                if (residual != 0) {
                    return false;
                }
                lower[row][column] = 0;
            } else {
                lower[row][column] = residual / pivot;
                if (lower[row][column] < 0) {
                    return false;
                }
            }
        }
    }
    return true;
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
                "usage: probe_su2_psd_boundary_obstruction "
                "[maximum_half_level] [maximum_coefficient]");
        }

        std::uint64_t negative_current_profiles = 0U;
        std::uint64_t psd_profiles = 0U;
        std::uint64_t psd_tp2_failures = 0U;
        std::uint64_t psd_insertion_checks = 0U;
        std::uint64_t psd_insertion_failures = 0U;
        std::uint64_t psd_full_insertion_coordinates = 0U;
        std::uint64_t psd_full_insertion_failures = 0U;
        std::uint64_t psd_row_identity_checks = 0U;
        std::uint64_t psd_row_identity_failures = 0U;
        std::uint64_t psd_wall_reflection_checks = 0U;
        std::uint64_t psd_wall_reflection_failures = 0U;
        std::uint64_t psd_mixed_boundary_checks = 0U;
        std::uint64_t psd_mixed_boundary_failures = 0U;
        std::uint64_t psd_residual_columns = 0U;
        std::uint64_t psd_residual_unimodality_failures = 0U;
        std::uint64_t psd_nonnegative_ldl_failures = 0U;
        std::uint64_t psd_local_payment_checks = 0U;
        std::uint64_t psd_local_payment_failures = 0U;
        std::uint64_t psd_cumulative_payment_checks = 0U;
        std::uint64_t psd_cumulative_payment_failures = 0U;
        std::uint64_t psd_paired_payment_checks = 0U;
        std::uint64_t psd_paired_payment_failures = 0U;
        std::uint64_t canonical_pair_obstructions = 0U;
        std::uint64_t canonical_pair_misses = 0U;
        std::uint64_t negative_current_coordinates = 0U;
        std::uint64_t canonical_coordinate_obstructions = 0U;
        std::uint64_t canonical_coordinate_misses = 0U;
        std::uint64_t first_defect_categories[4][3]{};
        std::uint64_t row_one_first_defects = 0U;
        std::uint64_t row_one_predecessor_obstructions = 0U;
        std::uint64_t row_one_predecessor_misses = 0U;
        Vector first_row_one_predecessor_miss;
        std::uint64_t row_one_balanced_obstructions = 0U;
        std::uint64_t row_one_balanced_misses = 0U;
        Vector first_row_one_balanced_miss;
        std::uint64_t diagonal_admissible_negative_profiles = 0U;
        std::vector<std::uint64_t> first_obstruction_order(
            static_cast<std::size_t>(maximum_level + 2), 0U);
        std::vector<std::uint64_t> diagonal_obstruction_order(
            static_cast<std::size_t>(maximum_level + 2), 0U);
        Vector first_profile;
        cpp_int first_current = 0;
        unsigned first_order = 0U;
        cpp_int first_minor = 0;
        std::vector<std::size_t> first_minor_indices;
        Vector first_diagonal_profile;
        cpp_int first_diagonal_current = 0;
        PrincipalObstruction first_diagonal_obstruction;
        Vector first_psd_tp2_profile;
        std::tuple<int, int, int, int> first_psd_tp2_indices{
            -1, -1, -1, -1};
        cpp_int first_psd_tp2_minor = 0;
        Vector first_psd_insertion_input;
        Vector first_psd_insertion_output;
        int first_psd_insertion_factor = -1;
        int first_psd_insertion_radius = -1;
        cpp_int first_psd_insertion_value = 0;
        Vector first_psd_full_insertion_input;
        Vector first_psd_full_insertion_output;
        int first_psd_full_insertion_factor = -1;
        int first_psd_full_insertion_radius = -1;
        int first_psd_full_insertion_target = -1;
        cpp_int first_psd_full_insertion_value = 0;
        Vector first_psd_local_payment_input;
        int first_psd_local_payment_factor = -1;
        int first_psd_local_payment_shell = -1;
        int first_psd_local_payment_radius = -1;
        cpp_int first_psd_local_payment_value = 0;
        Vector first_psd_cumulative_payment_input;
        int first_psd_cumulative_payment_factor = -1;
        int first_psd_cumulative_payment_shell = -1;
        int first_psd_cumulative_payment_radius = -1;
        cpp_int first_psd_cumulative_payment_value = 0;
        Vector first_psd_paired_payment_input;
        int first_psd_paired_payment_factor = -1;
        int first_psd_paired_payment_shell = -1;
        int first_psd_paired_payment_radius = -1;
        cpp_int first_psd_paired_payment_value = 0;
        Vector first_psd_mixed_boundary_input;
        int first_psd_mixed_boundary_radius = -1;
        int first_psd_mixed_boundary_depth = -1;
        cpp_int first_psd_mixed_boundary_value = 0;
        Vector first_psd_residual_unimodality_input;
        Vector first_psd_residual_unimodality_column;
        int first_psd_residual_unimodality_target = -1;
        Vector first_psd_nonnegative_ldl_failure;
        Vector first_canonical_pair_miss;
        std::tuple<int, int, int, int> first_canonical_pair_miss_data{
            -1, -1, -1, -1};

        for (int level = 1; level <= maximum_level; ++level) {
            const std::uint64_t base
                = static_cast<std::uint64_t>(maximum_coefficient + 1);
            std::uint64_t count = 1U;
            for (int index = 0; index <= level; ++index) {
                count *= base;
            }
            for (std::uint64_t code = 1U; code < count; ++code) {
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
                    boundary = boundary
                        && values[0]
                            * values[static_cast<std::size_t>(
                                level - radius)]
                            >= values[static_cast<std::size_t>(radius)]
                                * values[static_cast<std::size_t>(level)];
                }
                if (!boundary) {
                    continue;
                }

                const Matrix matrix = multiplication_matrix(values);
                const PrincipalObstruction psd_obstruction
                    = first_negative_principal_minor(matrix);
                if (psd_obstruction.order == 0U) {
                    ++psd_profiles;
                    Matrix anchored = matrix;
                    for (int row = 0; row <= level; ++row) {
                        for (int column = 0; column <= level; ++column) {
                            anchored[static_cast<std::size_t>(row)]
                                    [static_cast<std::size_t>(column)]
                                = values[0]
                                    * matrix[
                                        static_cast<std::size_t>(row)]
                                        [static_cast<std::size_t>(column)]
                                  - values[static_cast<std::size_t>(row)]
                                    * values[
                                        static_cast<std::size_t>(column)];
                        }
                    }
                    if (!has_nonnegative_natural_ldl(anchored)) {
                        ++psd_nonnegative_ldl_failures;
                        if (first_psd_nonnegative_ldl_failure.empty()) {
                            first_psd_nonnegative_ldl_failure = values;
                        }
                    }
                    Vector boundary_currents(
                        static_cast<std::size_t>(level + 1), 0);
                    for (int radius = 0; radius <= level; ++radius) {
                        boundary_currents[static_cast<std::size_t>(radius)]
                            = values[0]
                                * values[static_cast<std::size_t>(
                                    level - radius)]
                              - values[static_cast<std::size_t>(radius)]
                                * values[static_cast<std::size_t>(level)];
                    }
                    for (int radius = 0; radius <= level; ++radius) {
                        const int maximum_depth
                            = std::min(radius, level - radius);
                        for (int depth = 0;
                             depth <= maximum_depth;
                             ++depth) {
                            cpp_int boundary_block = 0;
                            for (int index = radius - depth;
                                 index <= radius + depth;
                                 ++index) {
                                boundary_block
                                    += boundary_currents[
                                        static_cast<std::size_t>(index)];
                            }
                            const cpp_int mixed_boundary
                                = values[0] * boundary_block
                                  - values[static_cast<std::size_t>(radius)]
                                    * boundary_currents[
                                        static_cast<std::size_t>(depth)];
                            ++psd_mixed_boundary_checks;
                            if (mixed_boundary < 0) {
                                ++psd_mixed_boundary_failures;
                                if (first_psd_mixed_boundary_input.empty()) {
                                    first_psd_mixed_boundary_input = values;
                                    first_psd_mixed_boundary_radius = radius;
                                    first_psd_mixed_boundary_depth = depth;
                                    first_psd_mixed_boundary_value
                                        = mixed_boundary;
                                }
                            }

                            const int wall_target = level - depth;
                            const cpp_int wall_current
                                = values[0]
                                    * matrix[
                                        static_cast<std::size_t>(radius)]
                                        [static_cast<std::size_t>(
                                            wall_target)]
                                  - values[static_cast<std::size_t>(radius)]
                                    * values[
                                        static_cast<std::size_t>(
                                            wall_target)];
                            const cpp_int no_wall_current
                                = values[0]
                                    * matrix[
                                        static_cast<std::size_t>(radius)]
                                        [static_cast<std::size_t>(depth)]
                                  - values[static_cast<std::size_t>(radius)]
                                    * values[
                                        static_cast<std::size_t>(depth)];
                            ++psd_wall_reflection_checks;
                            if (values[0] * wall_current
                                != mixed_boundary
                                    + values[
                                        static_cast<std::size_t>(level)]
                                        * no_wall_current) {
                                ++psd_wall_reflection_failures;
                            }
                        }
                    }
                    for (int target = 0; target <= level; ++target) {
                        ++psd_residual_columns;
                        Vector residual(
                            static_cast<std::size_t>(level + 1), 0);
                        for (int radius = 0; radius <= level; ++radius) {
                            residual[static_cast<std::size_t>(radius)]
                                = values[0]
                                    * matrix[
                                        static_cast<std::size_t>(radius)]
                                        [static_cast<std::size_t>(target)]
                                  - values[
                                        static_cast<std::size_t>(radius)]
                                    * values[
                                        static_cast<std::size_t>(target)];
                        }
                        bool decreasing_seen = false;
                        bool unimodal = true;
                        for (int radius = 0; radius < level; ++radius) {
                            if (values[
                                    static_cast<std::size_t>(radius)]
                                    == 0
                                || values[
                                    static_cast<std::size_t>(radius + 1)]
                                    == 0) {
                                continue;
                            }
                            const cpp_int ratio_difference
                                = residual[
                                    static_cast<std::size_t>(radius + 1)]
                                    * values[
                                        static_cast<std::size_t>(radius)]
                                  - residual[
                                        static_cast<std::size_t>(radius)]
                                    * values[
                                        static_cast<std::size_t>(
                                            radius + 1)];
                            if (ratio_difference < 0) {
                                decreasing_seen = true;
                            } else if (
                                ratio_difference > 0
                                && decreasing_seen) {
                                unimodal = false;
                                break;
                            }
                        }
                        if (!unimodal) {
                            ++psd_residual_unimodality_failures;
                            if (first_psd_residual_unimodality_input.empty()) {
                                first_psd_residual_unimodality_input = values;
                                first_psd_residual_unimodality_column
                                    = residual;
                                first_psd_residual_unimodality_target
                                    = target;
                            }
                        }
                    }
                    for (int factor = 1; factor <= level; ++factor) {
                        ++psd_insertion_checks;
                        const int width = 2 * std::min(
                            factor, level - factor);
                        cpp_int prefix = 0;
                        cpp_int suffix = 0;
                        for (int shell = 0; shell <= width; ++shell) {
                            prefix += values[
                                static_cast<std::size_t>(shell)];
                            suffix += values[static_cast<std::size_t>(
                                level - shell)];
                        }
                        Vector cumulative_payment(
                            static_cast<std::size_t>(level + 1), 0);
                        Matrix local_payments(
                            static_cast<std::size_t>(width + 1),
                            Vector(static_cast<std::size_t>(level + 1), 0));
                        for (int shell = 0; shell <= width; ++shell) {
                            for (int radius = 0;
                                 radius <= level;
                                 ++radius) {
                                ++psd_local_payment_checks;
                                const cpp_int top_current
                                    = values[0]
                                        * matrix[
                                            static_cast<std::size_t>(
                                                level - shell)]
                                            [static_cast<std::size_t>(
                                                radius)]
                                      - values[static_cast<std::size_t>(
                                            level - shell)]
                                        * values[static_cast<std::size_t>(
                                            radius)];
                                const cpp_int bottom_current
                                    = values[0]
                                        * matrix[
                                            static_cast<std::size_t>(
                                                shell)]
                                            [static_cast<std::size_t>(
                                                radius)]
                                      - values[
                                            static_cast<std::size_t>(
                                                shell)]
                                        * values[static_cast<std::size_t>(
                                            radius)];
                                const cpp_int local_payment
                                    = prefix * top_current
                                      - suffix * bottom_current;
                                local_payments[
                                    static_cast<std::size_t>(shell)]
                                    [static_cast<std::size_t>(radius)]
                                    = local_payment;
                                cumulative_payment[
                                    static_cast<std::size_t>(radius)]
                                    += local_payment;
                                ++psd_cumulative_payment_checks;
                                if (cumulative_payment[
                                        static_cast<std::size_t>(radius)]
                                    < 0) {
                                    ++psd_cumulative_payment_failures;
                                    if (
                                        first_psd_cumulative_payment_input
                                            .empty()) {
                                        first_psd_cumulative_payment_input
                                            = values;
                                        first_psd_cumulative_payment_factor
                                            = factor;
                                        first_psd_cumulative_payment_shell
                                            = shell;
                                        first_psd_cumulative_payment_radius
                                            = radius;
                                        first_psd_cumulative_payment_value
                                            = cumulative_payment[
                                                static_cast<std::size_t>(
                                                    radius)];
                                    }
                                }
                                if (local_payment >= 0) {
                                    continue;
                                }
                                ++psd_local_payment_failures;
                                if (first_psd_local_payment_input.empty()) {
                                    first_psd_local_payment_input = values;
                                    first_psd_local_payment_factor = factor;
                                    first_psd_local_payment_shell = shell;
                                    first_psd_local_payment_radius = radius;
                                    first_psd_local_payment_value
                                        = local_payment;
                                }
                            }
                        }
                        for (int shell = 0;
                             2 * shell <= width;
                             ++shell) {
                            const int reflected_shell = width - shell;
                            for (int radius = 0;
                                 radius <= level;
                                 ++radius) {
                                ++psd_paired_payment_checks;
                                cpp_int paired_payment
                                    = local_payments[
                                        static_cast<std::size_t>(shell)]
                                        [static_cast<std::size_t>(radius)];
                                if (reflected_shell != shell) {
                                    paired_payment
                                        += local_payments[
                                            static_cast<std::size_t>(
                                                reflected_shell)]
                                            [static_cast<std::size_t>(
                                                radius)];
                                }
                                if (paired_payment >= 0) {
                                    continue;
                                }
                                ++psd_paired_payment_failures;
                                if (first_psd_paired_payment_input.empty()) {
                                    first_psd_paired_payment_input = values;
                                    first_psd_paired_payment_factor = factor;
                                    first_psd_paired_payment_shell = shell;
                                    first_psd_paired_payment_radius = radius;
                                    first_psd_paired_payment_value
                                        = paired_payment;
                                }
                            }
                        }
                        const Vector output = fusion_step(
                            fusion_step(values, level, factor),
                            level,
                            factor);
                        const Matrix output_matrix
                            = multiplication_matrix(output);
                        for (int radius = 0;
                             radius <= level;
                             ++radius) {
                            for (int target = 0;
                                 target <= level;
                                 ++target) {
                                ++psd_full_insertion_coordinates;
                                const cpp_int current
                                    = output[0]
                                        * output_matrix[
                                            static_cast<std::size_t>(
                                                radius)]
                                            [static_cast<std::size_t>(
                                                target)]
                                      - output[
                                            static_cast<std::size_t>(
                                                radius)]
                                        * output[
                                            static_cast<std::size_t>(
                                                target)];
                                if (current >= 0) {
                                    continue;
                                }
                                ++psd_full_insertion_failures;
                                if (first_psd_full_insertion_input.empty()) {
                                    first_psd_full_insertion_input = values;
                                    first_psd_full_insertion_output = output;
                                    first_psd_full_insertion_factor = factor;
                                    first_psd_full_insertion_radius = radius;
                                    first_psd_full_insertion_target = target;
                                    first_psd_full_insertion_value = current;
                                }
                            }
                        }
                        for (int radius = 0;
                             radius <= level;
                             ++radius) {
                            ++psd_row_identity_checks;
                            cpp_int image = 0;
                            for (int shell = 0;
                                 shell <= width;
                                 ++shell) {
                                image += prefix
                                    * matrix[
                                        static_cast<std::size_t>(
                                            radius)]
                                        [static_cast<std::size_t>(
                                            level - shell)];
                                image -= suffix
                                    * matrix[
                                        static_cast<std::size_t>(
                                            radius)]
                                        [static_cast<std::size_t>(
                                            shell)];
                            }
                            const cpp_int boundary_current
                                = output[0]
                                    * output[static_cast<std::size_t>(
                                        level - radius)]
                                  - output[
                                        static_cast<std::size_t>(radius)]
                                    * output[
                                        static_cast<std::size_t>(level)];
                            if (image != boundary_current) {
                                ++psd_row_identity_failures;
                            }
                        }
                        for (int radius = 0; radius <= level; ++radius) {
                            const cpp_int boundary_current
                                = output[0]
                                    * output[static_cast<std::size_t>(
                                        level - radius)]
                                  - output[static_cast<std::size_t>(radius)]
                                    * output[
                                        static_cast<std::size_t>(level)];
                            if (boundary_current >= 0) {
                                continue;
                            }
                            ++psd_insertion_failures;
                            if (first_psd_insertion_input.empty()) {
                                first_psd_insertion_input = values;
                                first_psd_insertion_output = output;
                                first_psd_insertion_factor = factor;
                                first_psd_insertion_radius = radius;
                                first_psd_insertion_value
                                    = boundary_current;
                            }
                            break;
                        }
                    }
                    bool tp2 = true;
                    for (int row1 = 0;
                         row1 <= level && tp2;
                         ++row1) {
                        for (int row2 = row1 + 1;
                             row2 <= level && tp2;
                             ++row2) {
                            for (int column1 = 0;
                                 column1 <= level && tp2;
                                 ++column1) {
                                for (int column2 = column1 + 1;
                                     column2 <= level;
                                     ++column2) {
                                    const cpp_int minor
                                        = matrix[static_cast<std::size_t>(
                                            row1)]
                                                [static_cast<std::size_t>(
                                                    column1)]
                                            * matrix[
                                                static_cast<std::size_t>(
                                                    row2)]
                                                [static_cast<std::size_t>(
                                                    column2)]
                                        - matrix[
                                            static_cast<std::size_t>(
                                                row1)]
                                                [static_cast<std::size_t>(
                                                    column2)]
                                            * matrix[
                                                static_cast<std::size_t>(
                                                    row2)]
                                                [static_cast<std::size_t>(
                                                    column1)];
                                    if (minor >= 0) {
                                        continue;
                                    }
                                    tp2 = false;
                                    ++psd_tp2_failures;
                                    if (first_psd_tp2_profile.empty()) {
                                        first_psd_tp2_profile = values;
                                        first_psd_tp2_indices = {
                                            row1,
                                            row2,
                                            column1,
                                            column2};
                                        first_psd_tp2_minor = minor;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
                bool diagonal_admissible = true;
                for (int radius = 0; radius <= level; ++radius) {
                    diagonal_admissible = diagonal_admissible
                        && values[0]
                            * matrix[static_cast<std::size_t>(radius)]
                                    [static_cast<std::size_t>(radius)]
                            >= values[static_cast<std::size_t>(radius)]
                                * values[static_cast<std::size_t>(radius)];
                }
                cpp_int bad_current = 0;
                int bad_radius = -1;
                int bad_target = -1;
                for (int radius = 0; radius <= level; ++radius) {
                    for (int target = 0; target <= level; ++target) {
                        const cpp_int current
                            = values[0]
                                * matrix[static_cast<std::size_t>(radius)]
                                        [static_cast<std::size_t>(target)]
                            - values[static_cast<std::size_t>(radius)]
                                * values[static_cast<std::size_t>(target)];
                        if (current >= 0) {
                            continue;
                        }
                        ++negative_current_coordinates;
                        if (bad_current == 0) {
                            bad_current = current;
                            bad_radius = radius;
                            bad_target = target;
                        }
                        int pair_left = std::abs(radius - target);
                        int pair_right = std::min(radius, target);
                        if (pair_left == pair_right) {
                            pair_left = 0;
                            pair_right = std::max(radius, target);
                        } else if (pair_left > pair_right) {
                            std::swap(pair_left, pair_right);
                        }
                        const cpp_int pair_minor
                            = matrix[static_cast<std::size_t>(pair_left)]
                                    [static_cast<std::size_t>(pair_left)]
                                * matrix[
                                    static_cast<std::size_t>(pair_right)]
                                    [static_cast<std::size_t>(pair_right)]
                              - matrix[
                                    static_cast<std::size_t>(pair_left)]
                                    [static_cast<std::size_t>(pair_right)]
                                * matrix[
                                    static_cast<std::size_t>(pair_right)]
                                    [static_cast<std::size_t>(pair_left)];
                        if (pair_minor < 0) {
                            ++canonical_coordinate_obstructions;
                        } else {
                            ++canonical_coordinate_misses;
                        }
                    }
                }
                if (bad_current == 0) {
                    continue;
                }
                ++negative_current_profiles;
                int canonical_left = std::abs(bad_radius - bad_target);
                int canonical_right = std::min(bad_radius, bad_target);
                if (canonical_left == canonical_right) {
                    canonical_left = 0;
                    canonical_right = std::max(bad_radius, bad_target);
                } else if (canonical_left > canonical_right) {
                    std::swap(canonical_left, canonical_right);
                }
                bool canonical_negative = false;
                if (canonical_left != canonical_right) {
                    const cpp_int canonical_minor
                        = matrix[static_cast<std::size_t>(canonical_left)]
                                [static_cast<std::size_t>(canonical_left)]
                            * matrix[
                                static_cast<std::size_t>(canonical_right)]
                                [static_cast<std::size_t>(canonical_right)]
                          - matrix[
                                static_cast<std::size_t>(canonical_left)]
                                [static_cast<std::size_t>(canonical_right)]
                            * matrix[
                                static_cast<std::size_t>(canonical_right)]
                                [static_cast<std::size_t>(canonical_left)];
                    canonical_negative = canonical_minor < 0;
                }
                if (canonical_negative) {
                    ++canonical_pair_obstructions;
                } else {
                    ++canonical_pair_misses;
                    if (first_canonical_pair_miss.empty()) {
                        first_canonical_pair_miss = values;
                        first_canonical_pair_miss_data = {
                            bad_radius,
                            bad_target,
                            canonical_left,
                            canonical_right};
                    }
                }
                const int difference = bad_target - bad_radius;
                const int rotation_category
                    = bad_target == bad_radius ? 0
                    : difference == bad_radius ? 1
                    : difference < bad_radius ? 2
                    : 3;
                const int wall_category
                    = bad_radius + bad_target < level ? 0
                    : bad_radius + bad_target == level ? 1
                    : 2;
                ++first_defect_categories[rotation_category][wall_category];
                if (bad_radius == 1 && bad_target >= 2) {
                    ++row_one_first_defects;
                    std::vector<std::size_t> predecessor_indices;
                    if (bad_target == 2) {
                        predecessor_indices = {
                            0U, static_cast<std::size_t>(bad_target)};
                    } else {
                        predecessor_indices = {
                            0U,
                            1U,
                            static_cast<std::size_t>(bad_target - 1)};
                    }
                    Matrix predecessor(
                        predecessor_indices.size(),
                        Vector(predecessor_indices.size(), 0));
                    for (std::size_t row = 0U;
                         row < predecessor_indices.size();
                         ++row) {
                        for (std::size_t column = 0U;
                             column < predecessor_indices.size();
                             ++column) {
                            predecessor[row][column]
                                = matrix[predecessor_indices[row]]
                                        [predecessor_indices[column]];
                        }
                    }
                    if (determinant(predecessor) < 0) {
                        ++row_one_predecessor_obstructions;
                    } else {
                        ++row_one_predecessor_misses;
                        if (first_row_one_predecessor_miss.empty()) {
                            first_row_one_predecessor_miss = values;
                        }
                    }
                    const int balanced_left = (bad_target - 1) / 2;
                    const int balanced_right
                        = bad_target - balanced_left;
                    std::vector<std::size_t> balanced_indices{0U};
                    if (balanced_left != 0) {
                        balanced_indices.push_back(
                            static_cast<std::size_t>(balanced_left));
                    }
                    balanced_indices.push_back(
                        static_cast<std::size_t>(balanced_right));
                    Matrix balanced(
                        balanced_indices.size(),
                        Vector(balanced_indices.size(), 0));
                    for (std::size_t row = 0U;
                         row < balanced_indices.size();
                         ++row) {
                        for (std::size_t column = 0U;
                             column < balanced_indices.size();
                             ++column) {
                            balanced[row][column]
                                = matrix[balanced_indices[row]]
                                        [balanced_indices[column]];
                        }
                    }
                    if (determinant(balanced) < 0) {
                        ++row_one_balanced_obstructions;
                    } else {
                        ++row_one_balanced_misses;
                        if (first_row_one_balanced_miss.empty()) {
                            first_row_one_balanced_miss = values;
                        }
                    }
                }
                const PrincipalObstruction& obstruction = psd_obstruction;
                ++first_obstruction_order[obstruction.order];
                if (first_profile.empty()) {
                    first_profile = values;
                    first_current = bad_current;
                    first_order = obstruction.order;
                    first_minor = obstruction.value;
                    first_minor_indices = obstruction.indices;
                }
                if (diagonal_admissible) {
                    ++diagonal_admissible_negative_profiles;
                    ++diagonal_obstruction_order[obstruction.order];
                    if (first_diagonal_profile.empty()) {
                        first_diagonal_profile = values;
                        first_diagonal_current = bad_current;
                        first_diagonal_obstruction = obstruction;
                    }
                }
            }
        }

        std::cout << "SU2_PSD_BOUNDARY_OBSTRUCTION"
                  << " maximum_level=" << 2 * maximum_level
                  << " maximum_coefficient=" << maximum_coefficient
                  << " negative_current_profiles="
                  << negative_current_profiles
                  << " psd_profiles=" << psd_profiles
                  << " psd_tp2_failures=" << psd_tp2_failures
                  << " psd_insertion_checks=" << psd_insertion_checks
                  << " psd_insertion_failures="
                  << psd_insertion_failures
                  << " psd_full_insertion_coordinates="
                  << psd_full_insertion_coordinates
                  << " psd_full_insertion_failures="
                  << psd_full_insertion_failures
                  << " psd_row_identity_checks="
                  << psd_row_identity_checks
                  << " psd_row_identity_failures="
                  << psd_row_identity_failures
                  << " psd_wall_reflection_checks="
                  << psd_wall_reflection_checks
                  << " psd_wall_reflection_failures="
                  << psd_wall_reflection_failures
                  << " psd_mixed_boundary_checks="
                  << psd_mixed_boundary_checks
                  << " psd_mixed_boundary_failures="
                  << psd_mixed_boundary_failures
                  << " psd_residual_columns="
                  << psd_residual_columns
                  << " psd_residual_unimodality_failures="
                  << psd_residual_unimodality_failures
                  << " psd_nonnegative_ldl_failures="
                  << psd_nonnegative_ldl_failures
                  << " psd_local_payment_checks="
                  << psd_local_payment_checks
                  << " psd_local_payment_failures="
                  << psd_local_payment_failures
                  << " psd_cumulative_payment_checks="
                  << psd_cumulative_payment_checks
                  << " psd_cumulative_payment_failures="
                  << psd_cumulative_payment_failures
                  << " psd_paired_payment_checks="
                  << psd_paired_payment_checks
                  << " psd_paired_payment_failures="
                  << psd_paired_payment_failures
                  << " canonical_pair_obstructions="
                  << canonical_pair_obstructions
                  << " canonical_pair_misses=" << canonical_pair_misses
                  << " negative_current_coordinates="
                  << negative_current_coordinates
                  << " canonical_coordinate_obstructions="
                  << canonical_coordinate_obstructions
                  << " canonical_coordinate_misses="
                  << canonical_coordinate_misses
                  << " row_one_first_defects="
                  << row_one_first_defects
                  << " row_one_predecessor_obstructions="
                  << row_one_predecessor_obstructions
                  << " row_one_predecessor_misses="
                  << row_one_predecessor_misses
                  << " row_one_balanced_obstructions="
                  << row_one_balanced_obstructions
                  << " row_one_balanced_misses="
                  << row_one_balanced_misses
                  << " diagonal_admissible_negative_profiles="
                  << diagonal_admissible_negative_profiles
                  << " obstruction_orders=[";
        for (std::size_t order = 0U;
             order < first_obstruction_order.size();
             ++order) {
            if (order != 0U) {
                std::cout << ',';
            }
            std::cout << first_obstruction_order[order];
        }
        std::cout << "] diagonal_obstruction_orders=[";
        for (std::size_t order = 0U;
             order < diagonal_obstruction_order.size();
             ++order) {
            if (order != 0U) {
                std::cout << ',';
            }
            std::cout << diagonal_obstruction_order[order];
        }
        std::cout << "] first_profile=" << render(first_profile)
                  << " first_current=" << first_current
                  << " first_order=" << first_order
                  << " first_minor=" << first_minor
                  << " first_minor_indices=[";
        for (std::size_t index = 0U;
             index < first_minor_indices.size();
             ++index) {
            if (index != 0U) {
                std::cout << ',';
            }
            std::cout << first_minor_indices[index];
        }
        std::cout << "] first_diagonal_profile="
                  << render(first_diagonal_profile)
                  << " first_diagonal_current=" << first_diagonal_current
                  << " first_diagonal_minor="
                  << first_diagonal_obstruction.value
                  << " first_diagonal_minor_indices=[";
        for (std::size_t index = 0U;
             index < first_diagonal_obstruction.indices.size();
             ++index) {
            if (index != 0U) {
                std::cout << ',';
            }
            std::cout << first_diagonal_obstruction.indices[index];
        }
        std::cout << "]\n";
        const auto [tp2_row1, tp2_row2, tp2_column1, tp2_column2]
            = first_psd_tp2_indices;
        std::cout << "FIRST_PSD_TP2_FAILURE"
                  << " profile=" << render(first_psd_tp2_profile)
                  << " indices=[" << tp2_row1 << ',' << tp2_row2 << ','
                  << tp2_column1 << ',' << tp2_column2 << ']'
                  << " minor=" << first_psd_tp2_minor << '\n';
        std::cout << "FIRST_PSD_INSERTION_FAILURE"
                  << " factor=" << first_psd_insertion_factor
                  << " radius=" << first_psd_insertion_radius
                  << " value=" << first_psd_insertion_value
                  << " input=" << render(first_psd_insertion_input)
                  << " output=" << render(first_psd_insertion_output)
                  << '\n';
        std::cout << "FIRST_PSD_FULL_INSERTION_FAILURE"
                  << " factor=" << first_psd_full_insertion_factor
                  << " radius=" << first_psd_full_insertion_radius
                  << " target=" << first_psd_full_insertion_target
                  << " value=" << first_psd_full_insertion_value
                  << " input="
                  << render(first_psd_full_insertion_input)
                  << " output="
                  << render(first_psd_full_insertion_output) << '\n';
        std::cout << "FIRST_PSD_LOCAL_PAYMENT_FAILURE"
                  << " factor=" << first_psd_local_payment_factor
                  << " shell=" << first_psd_local_payment_shell
                  << " radius=" << first_psd_local_payment_radius
                  << " value=" << first_psd_local_payment_value
                  << " input="
                  << render(first_psd_local_payment_input) << '\n';
        std::cout << "FIRST_PSD_CUMULATIVE_PAYMENT_FAILURE"
                  << " factor=" << first_psd_cumulative_payment_factor
                  << " shell=" << first_psd_cumulative_payment_shell
                  << " radius="
                  << first_psd_cumulative_payment_radius
                  << " value="
                  << first_psd_cumulative_payment_value
                  << " input="
                  << render(first_psd_cumulative_payment_input) << '\n';
        std::cout << "FIRST_PSD_PAIRED_PAYMENT_FAILURE"
                  << " factor=" << first_psd_paired_payment_factor
                  << " shell=" << first_psd_paired_payment_shell
                  << " radius=" << first_psd_paired_payment_radius
                  << " value=" << first_psd_paired_payment_value
                  << " input="
                  << render(first_psd_paired_payment_input) << '\n';
        std::cout << "FIRST_PSD_MIXED_BOUNDARY_FAILURE"
                  << " radius=" << first_psd_mixed_boundary_radius
                  << " depth=" << first_psd_mixed_boundary_depth
                  << " value=" << first_psd_mixed_boundary_value
                  << " input="
                  << render(first_psd_mixed_boundary_input) << '\n';
        std::cout << "FIRST_PSD_RESIDUAL_UNIMODALITY_FAILURE"
                  << " target="
                  << first_psd_residual_unimodality_target
                  << " input="
                  << render(first_psd_residual_unimodality_input)
                  << " residual="
                  << render(first_psd_residual_unimodality_column) << '\n';
        std::cout << "FIRST_PSD_NONNEGATIVE_LDL_FAILURE"
                  << " input="
                  << render(first_psd_nonnegative_ldl_failure) << '\n';
        const auto [
            miss_radius,
            miss_target,
            miss_left,
            miss_right] = first_canonical_pair_miss_data;
        std::cout << "FIRST_CANONICAL_PAIR_MISS"
                  << " profile=" << render(first_canonical_pair_miss)
                  << " current_indices=[" << miss_radius << ','
                  << miss_target << ']'
                  << " canonical_indices=[" << miss_left << ','
                  << miss_right << "]\n";
        std::cout << "FIRST_DEFECT_CATEGORIES"
                  << " rows=[diagonal,equal_rotation,inner_rotation,"
                     "outer_rotation]"
                  << " columns=[below_wall,on_wall,beyond_wall]"
                  << " counts=[";
        for (int row = 0; row < 4; ++row) {
            if (row != 0) {
                std::cout << ',';
            }
            std::cout << '[';
            for (int column = 0; column < 3; ++column) {
                if (column != 0) {
                    std::cout << ',';
                }
                std::cout << first_defect_categories[row][column];
            }
            std::cout << ']';
        }
        std::cout << "]\n";
        std::cout << "FIRST_ROW_ONE_PREDECESSOR_MISS"
                  << " profile="
                  << render(first_row_one_predecessor_miss) << '\n';
        std::cout << "FIRST_ROW_ONE_BALANCED_MISS"
                  << " profile="
                  << render(first_row_one_balanced_miss) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
