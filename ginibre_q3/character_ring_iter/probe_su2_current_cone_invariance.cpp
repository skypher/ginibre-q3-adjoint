#include <boost/multiprecision/cpp_int.hpp>

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
        const int upper = std::min(source + factor, 2 * level - source - factor);
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

Vector fusion_square(
    const Vector& vector,
    const std::vector<Matrix>& fusion) {
    Vector square(vector.size(), 0);
    for (std::size_t factor = 0U; factor < vector.size(); ++factor) {
        const Vector translated = multiply(fusion[factor], vector);
        for (std::size_t target = 0U; target < vector.size(); ++target) {
            square[target] += vector[factor] * translated[target];
        }
    }
    return square;
}

bool log_concave(const Vector& vector) {
    bool positive_seen = false;
    bool support_ended = false;
    for (std::size_t index = 0U; index < vector.size(); ++index) {
        if (vector[index] > 0) {
            if (support_ended) {
                return false;
            }
            positive_seen = true;
        } else if (positive_seen) {
            support_ended = true;
        }
        if (index > 0U && index + 1U < vector.size()
            && vector[index] * vector[index]
                < vector[index - 1U] * vector[index + 1U]) {
            return false;
        }
    }
    return true;
}

bool interval_support(const Vector& vector) {
    bool support_started = false;
    bool support_ended = false;
    for (const cpp_int& value : vector) {
        if (value > 0) {
            if (support_ended) {
                return false;
            }
            support_started = true;
        } else if (support_started) {
            support_ended = true;
        }
    }
    return support_started;
}

bool current_cone(
    const Vector& vector,
    const std::vector<Matrix>& fusion,
    int& bad_radius,
    int& bad_target,
    cpp_int& bad_value) {
    if (vector[0] <= 0) {
        return false;
    }
    for (std::size_t radius = 0U; radius < fusion.size(); ++radius) {
        const Vector translated = multiply(fusion[radius], vector);
        for (std::size_t target = 0U; target < vector.size(); ++target) {
            const cpp_int value
                = vector[0] * translated[target]
                - vector[radius] * vector[target];
            if (value < 0) {
                bad_radius = static_cast<int>(radius);
                bad_target = static_cast<int>(target);
                bad_value = value;
                return false;
            }
        }
    }
    return true;
}

bool boundary_cone(
    const Vector& vector,
    const std::vector<Matrix>& fusion,
    int& bad_radius,
    cpp_int& bad_value) {
    if (vector[0] <= 0) {
        return false;
    }
    const std::size_t wall = vector.size() - 1U;
    for (std::size_t radius = 0U; radius < fusion.size(); ++radius) {
        const Vector translated = multiply(fusion[radius], vector);
        const cpp_int value
            = vector[0] * translated[wall]
            - vector[radius] * vector[wall];
        if (value < 0) {
            bad_radius = static_cast<int>(radius);
            bad_value = value;
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

struct Failure {
    int level = -1;
    int factor = -1;
    int radius = -1;
    int target = -1;
    cpp_int value = 0;
    Vector input;
    Vector output;
};

void record(
    Failure& failure,
    int level,
    int factor,
    int radius,
    int target,
    const cpp_int& value,
    const Vector& input,
    const Vector& output) {
    if (failure.level >= 0) {
        return;
    }
    failure = {level, factor, radius, target, value, input, output};
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
                "usage: probe_su2_current_cone_invariance "
                "[maximum_half_level] [maximum_coefficient]");
        }

        std::uint64_t profiles = 0U;
        std::uint64_t cone_profiles = 0U;
        std::uint64_t transforms = 0U;
        std::uint64_t failures = 0U;
        std::uint64_t lc_profiles = 0U;
        std::uint64_t lc_transforms = 0U;
        std::uint64_t lc_failures = 0U;
        std::uint64_t square_boundary_profiles = 0U;
        std::uint64_t square_boundary_transforms = 0U;
        std::uint64_t square_boundary_failures = 0U;
        std::uint64_t square_full_failures = 0U;
        std::uint64_t diagonal_admissible_profiles = 0U;
        std::uint64_t diagonal_admissible_open_checks = 0U;
        std::uint64_t diagonal_admissible_open_failures = 0U;
        std::uint64_t boundary_profiles = 0U;
        std::uint64_t boundary_transforms = 0U;
        std::uint64_t boundary_failures = 0U;
        std::uint64_t cone_to_boundary_transforms = 0U;
        std::uint64_t cone_to_boundary_failures = 0U;
        std::uint64_t log_concave_cone_to_boundary_transforms = 0U;
        std::uint64_t log_concave_cone_to_boundary_failures = 0U;
        std::uint64_t zero_anchor_interval_transforms = 0U;
        std::uint64_t zero_anchor_interval_failures = 0U;
        Failure first;
        Failure first_lc;
        Failure first_square;
        Failure first_diagonal;
        Failure first_boundary;
        Failure first_cone_to_boundary;
        Failure first_log_concave_cone_to_boundary;
        Failure first_zero_anchor_interval;

        for (int level = 1; level <= maximum_level; ++level) {
            std::vector<Matrix> fusion;
            for (int factor = 0; factor <= level; ++factor) {
                fusion.push_back(fusion_matrix(level, factor));
            }
            std::uint64_t profile_count = 1U;
            const std::uint64_t base
                = static_cast<std::uint64_t>(maximum_coefficient + 1);
            for (int index = 0; index <= level; ++index) {
                if (profile_count
                    > std::numeric_limits<std::uint64_t>::max() / base) {
                    throw std::overflow_error("profile count overflow");
                }
                profile_count *= base;
            }

            for (std::uint64_t code = 0U; code < profile_count; ++code) {
                ++profiles;
                std::uint64_t remainder = code;
                Vector input(static_cast<std::size_t>(level + 1), 0);
                for (int index = 0; index <= level; ++index) {
                    input[static_cast<std::size_t>(index)] = remainder % base;
                    remainder /= base;
                }
                bool diagonal_admissible
                    = input[0] > 0 && log_concave(input);
                for (int diagonal = 1;
                     2 * diagonal <= level && diagonal_admissible;
                     ++diagonal) {
                    const Vector translated = multiply(
                        fusion[static_cast<std::size_t>(diagonal)], input);
                    if (input[0]
                            * translated[static_cast<std::size_t>(diagonal)]
                        < input[static_cast<std::size_t>(diagonal)]
                            * input[static_cast<std::size_t>(diagonal)]) {
                        diagonal_admissible = false;
                    }
                }
                if (diagonal_admissible) {
                    ++diagonal_admissible_profiles;
                    for (int left = 1; left <= level; ++left) {
                        const Vector translated = multiply(
                            fusion[static_cast<std::size_t>(left)], input);
                        for (int right = 1;
                             left + right <= level;
                             ++right) {
                            ++diagonal_admissible_open_checks;
                            const cpp_int open_current
                                = input[0]
                                    * translated[
                                        static_cast<std::size_t>(right)]
                                - input[static_cast<std::size_t>(left)]
                                    * input[static_cast<std::size_t>(right)];
                            if (open_current >= 0) {
                                continue;
                            }
                            ++diagonal_admissible_open_failures;
                            record(
                                first_diagonal,
                                level,
                                0,
                                left,
                                right,
                                open_current,
                                input,
                                input);
                        }
                    }
                }
                if (input[0] == 0 && interval_support(input)) {
                    for (int factor = 1; factor <= level; ++factor) {
                        const Vector output = multiply(
                            fusion[static_cast<std::size_t>(factor)],
                            multiply(
                                fusion[static_cast<std::size_t>(factor)],
                                input));
                        ++zero_anchor_interval_transforms;
                        for (int radius_index = 0;
                             radius_index <= level;
                             ++radius_index) {
                            const cpp_int reflected_current
                                = output[static_cast<std::size_t>(level)]
                                    * output[
                                        static_cast<std::size_t>(
                                            radius_index)]
                                - output[0]
                                    * output[
                                        static_cast<std::size_t>(
                                            level - radius_index)];
                            if (reflected_current >= 0) {
                                continue;
                            }
                            ++zero_anchor_interval_failures;
                            record(
                                first_zero_anchor_interval,
                                level,
                                factor,
                                radius_index,
                                level,
                                reflected_current,
                                input,
                                output);
                        }
                    }
                }
                int radius = -1;
                int target = -1;
                cpp_int value = 0;
                if (boundary_cone(input, fusion, radius, value)) {
                    ++boundary_profiles;
                    for (int factor = 1; factor <= level; ++factor) {
                        const Vector output = multiply(
                            fusion[static_cast<std::size_t>(factor)],
                            multiply(
                                fusion[static_cast<std::size_t>(factor)],
                                input));
                        ++boundary_transforms;
                        radius = -1;
                        value = 0;
                        if (boundary_cone(
                                output, fusion, radius, value)) {
                            continue;
                        }
                        ++boundary_failures;
                        record(
                            first_boundary,
                            level,
                            factor,
                            radius,
                            level,
                            value,
                            input,
                            output);
                    }
                }
                radius = -1;
                target = -1;
                value = 0;
                if (current_cone(input, fusion, radius, target, value)) {
                    ++cone_profiles;
                    const bool lc = log_concave(input);
                    if (lc) {
                        ++lc_profiles;
                    }
                    for (int factor = 1; factor <= level; ++factor) {
                        const Vector output = multiply(
                            fusion[static_cast<std::size_t>(factor)],
                            multiply(
                                fusion[static_cast<std::size_t>(factor)], input));
                        ++transforms;
                        if (lc) {
                            ++lc_transforms;
                        }
                        ++cone_to_boundary_transforms;
                        if (lc) {
                            ++log_concave_cone_to_boundary_transforms;
                        }
                        int boundary_radius = -1;
                        cpp_int boundary_value = 0;
                        if (!boundary_cone(
                                output,
                                fusion,
                                boundary_radius,
                                boundary_value)) {
                            ++cone_to_boundary_failures;
                            record(
                                first_cone_to_boundary,
                                level,
                                factor,
                                boundary_radius,
                                level,
                                boundary_value,
                                input,
                                output);
                            if (lc) {
                                ++log_concave_cone_to_boundary_failures;
                                record(
                                    first_log_concave_cone_to_boundary,
                                    level,
                                    factor,
                                    boundary_radius,
                                    level,
                                    boundary_value,
                                    input,
                                    output);
                            }
                        }
                        radius = -1;
                        target = -1;
                        value = 0;
                        if (current_cone(output, fusion, radius, target, value)) {
                            continue;
                        }
                        ++failures;
                        record(
                            first,
                            level,
                            factor,
                            radius,
                            target,
                            value,
                            input,
                            output);
                        if (lc) {
                            ++lc_failures;
                            record(
                                first_lc,
                                level,
                                factor,
                                radius,
                                target,
                                value,
                                input,
                                output);
                        }
                    }
                }

                const Vector square = fusion_square(input, fusion);
                radius = -1;
                target = -1;
                value = 0;
                if (!log_concave(square)
                    || !boundary_cone(square, fusion, radius, value)) {
                    continue;
                }
                ++square_boundary_profiles;
                for (int factor = 1; factor <= level; ++factor) {
                    const Vector output = multiply(
                        fusion[static_cast<std::size_t>(factor)],
                        multiply(fusion[static_cast<std::size_t>(factor)], square));
                    ++square_boundary_transforms;
                    radius = -1;
                    target = -1;
                    value = 0;
                    if (!boundary_cone(output, fusion, radius, value)) {
                        ++square_boundary_failures;
                    }
                    radius = -1;
                    target = -1;
                    value = 0;
                    if (!current_cone(output, fusion, radius, target, value)) {
                        ++square_full_failures;
                        record(
                            first_square,
                            level,
                            factor,
                            radius,
                            target,
                            value,
                            square,
                            output);
                    }
                }
            }
        }

        const auto doubled = [](int value) {
            return value < 0 ? -1 : 2 * value;
        };
        std::cout
            << "SU2_CURRENT_CONE_INVARIANCE"
            << " maximum_level=" << 2 * maximum_level
            << " maximum_coefficient=" << maximum_coefficient
            << " profiles=" << profiles
            << " cone_profiles=" << cone_profiles
            << " transforms=" << transforms
            << " failures=" << failures
            << " log_concave_cone_profiles=" << lc_profiles
            << " log_concave_transforms=" << lc_transforms
            << " log_concave_failures=" << lc_failures
            << " square_log_concave_boundary_profiles="
            << square_boundary_profiles
            << " square_log_concave_boundary_transforms="
            << square_boundary_transforms
            << " square_log_concave_boundary_failures="
            << square_boundary_failures
            << " square_log_concave_full_failures="
            << square_full_failures
            << " diagonal_admissible_profiles="
            << diagonal_admissible_profiles
            << " diagonal_admissible_open_checks="
            << diagonal_admissible_open_checks
            << " diagonal_admissible_open_failures="
            << diagonal_admissible_open_failures
            << " boundary_profiles=" << boundary_profiles
            << " boundary_transforms=" << boundary_transforms
            << " boundary_failures=" << boundary_failures
            << " cone_to_boundary_transforms="
            << cone_to_boundary_transforms
            << " cone_to_boundary_failures="
            << cone_to_boundary_failures
            << " log_concave_cone_to_boundary_transforms="
            << log_concave_cone_to_boundary_transforms
            << " log_concave_cone_to_boundary_failures="
            << log_concave_cone_to_boundary_failures
            << " zero_anchor_interval_transforms="
            << zero_anchor_interval_transforms
            << " zero_anchor_interval_failures="
            << zero_anchor_interval_failures
            << " first_failure=("
            << doubled(first.level) << ','
            << doubled(first.factor) << ','
            << doubled(first.radius) << ','
            << doubled(first.target) << ','
            << first.value << ')'
            << " first_input=" << render(first.input)
            << " first_output=" << render(first.output)
            << " first_log_concave_failure=("
            << doubled(first_lc.level) << ','
            << doubled(first_lc.factor) << ','
            << doubled(first_lc.radius) << ','
            << doubled(first_lc.target) << ','
            << first_lc.value << ')'
            << " first_log_concave_input=" << render(first_lc.input)
            << " first_log_concave_output=" << render(first_lc.output)
            << " first_square_failure=("
            << doubled(first_square.level) << ','
            << doubled(first_square.factor) << ','
            << doubled(first_square.radius) << ','
            << doubled(first_square.target) << ','
            << first_square.value << ')'
            << " first_square_input=" << render(first_square.input)
            << " first_square_output=" << render(first_square.output)
            << " first_diagonal_failure=("
            << doubled(first_diagonal.level) << ','
            << doubled(first_diagonal.radius) << ','
            << doubled(first_diagonal.target) << ','
            << first_diagonal.value << ')'
            << " first_diagonal_input=" << render(first_diagonal.input)
            << " first_boundary_failure=("
            << doubled(first_boundary.level) << ','
            << doubled(first_boundary.factor) << ','
            << doubled(first_boundary.radius) << ','
            << first_boundary.value << ')'
            << " first_boundary_input=" << render(first_boundary.input)
            << " first_boundary_output=" << render(first_boundary.output)
            << " first_cone_to_boundary_failure=("
            << doubled(first_cone_to_boundary.level) << ','
            << doubled(first_cone_to_boundary.factor) << ','
            << doubled(first_cone_to_boundary.radius) << ','
            << first_cone_to_boundary.value << ')'
            << " first_cone_to_boundary_input="
            << render(first_cone_to_boundary.input)
            << " first_cone_to_boundary_output="
            << render(first_cone_to_boundary.output)
            << " first_log_concave_cone_to_boundary_failure=("
            << doubled(first_log_concave_cone_to_boundary.level) << ','
            << doubled(first_log_concave_cone_to_boundary.factor) << ','
            << doubled(first_log_concave_cone_to_boundary.radius) << ','
            << first_log_concave_cone_to_boundary.value << ')'
            << " first_log_concave_cone_to_boundary_input="
            << render(first_log_concave_cone_to_boundary.input)
            << " first_log_concave_cone_to_boundary_output="
            << render(first_log_concave_cone_to_boundary.output)
            << " first_zero_anchor_interval_failure=("
            << doubled(first_zero_anchor_interval.level) << ','
            << doubled(first_zero_anchor_interval.factor) << ','
            << doubled(first_zero_anchor_interval.radius) << ','
            << first_zero_anchor_interval.value << ')'
            << " first_zero_anchor_interval_input="
            << render(first_zero_anchor_interval.input)
            << " first_zero_anchor_interval_output="
            << render(first_zero_anchor_interval.output)
            << " result="
            << (square_boundary_failures == 0U
                    && square_full_failures == 0U
                    ? "NO_SQUARE_BOUNDARY_OR_FULL_FAILURE"
                    : "SQUARE_BOUNDARY_OR_FULL_FAILURE")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
