#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

#include <algorithm>
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

Vector fusion_square(
    const Vector& profile,
    const std::vector<Matrix>& fusion) {
    Vector square(profile.size(), 0);
    for (std::size_t factor = 0U; factor < profile.size(); ++factor) {
        if (profile[factor] == 0) {
            continue;
        }
        const Vector translated = multiply(fusion[factor], profile);
        for (std::size_t target = 0U; target < profile.size(); ++target) {
            square[target] += profile[factor] * translated[target];
        }
    }
    return square;
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

int first_bad_boundary(const Vector& vector) {
    const std::size_t level = vector.size() - 1U;
    for (std::size_t radius = 0U; radius <= level; ++radius) {
        if (vector[0] * vector[level - radius]
            < vector[radius] * vector[level]) {
            return static_cast<int>(radius);
        }
    }
    return -1;
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
    std::uint64_t code = 0U;
    int factor = -1;
    int radius = -1;
    cpp_int value = 0;
    Vector profile;
    Vector square;
    Vector output;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 8;
        const int maximum_coefficient = argc >= 3
            ? positive_argument(argv[2], "maximum_coefficient")
            : 4;
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_square_boundary_insertion_extended "
                "[maximum_half_level] [maximum_coefficient]");
        }

        std::uint64_t profiles = 0U;
        std::uint64_t admissible_squares = 0U;
        std::uint64_t insertion_checks = 0U;
        std::uint64_t boundary_failures = 0U;
        Failure first_failure;

        for (int level = 1; level <= maximum_level; ++level) {
            std::vector<Matrix> fusion;
            for (int factor = 0; factor <= level; ++factor) {
                fusion.push_back(fusion_matrix(level, factor));
            }
            const std::uint64_t base
                = static_cast<std::uint64_t>(maximum_coefficient + 1);
            std::uint64_t profile_count = 1U;
            for (int index = 0; index <= level; ++index) {
                if (profile_count
                    > static_cast<std::uint64_t>(
                        std::numeric_limits<std::int64_t>::max())
                        / base) {
                    throw std::overflow_error("profile count overflow");
                }
                profile_count *= base;
            }
            profiles += profile_count - 1U;

#pragma omp parallel for schedule(dynamic, 256) \
    reduction(+:admissible_squares,insertion_checks,boundary_failures)
            for (std::int64_t signed_code = 1;
                 signed_code < static_cast<std::int64_t>(profile_count);
                 ++signed_code) {
                const std::uint64_t code
                    = static_cast<std::uint64_t>(signed_code);
                std::uint64_t remainder = code;
                Vector profile(static_cast<std::size_t>(level + 1), 0);
                for (int index = 0; index <= level; ++index) {
                    profile[static_cast<std::size_t>(index)]
                        = remainder % base;
                    remainder /= base;
                }
                const Vector square = fusion_square(profile, fusion);
                if (!log_concave(square)
                    || first_bad_boundary(square) >= 0) {
                    continue;
                }
                ++admissible_squares;
                for (int factor = 1; factor <= level; ++factor) {
                    ++insertion_checks;
                    const Vector output = multiply(
                        fusion[static_cast<std::size_t>(factor)],
                        multiply(
                            fusion[static_cast<std::size_t>(factor)],
                            square));
                    const int radius = first_bad_boundary(output);
                    if (radius < 0) {
                        continue;
                    }
                    ++boundary_failures;
                    const cpp_int value
                        = output[0]
                            * output[static_cast<std::size_t>(
                                level - radius)]
                          - output[static_cast<std::size_t>(radius)]
                            * output[static_cast<std::size_t>(level)];
#pragma omp critical
                    {
                        const auto candidate = std::tuple{
                            level, code, factor, radius};
                        const auto current = std::tuple{
                            first_failure.level,
                            first_failure.code,
                            first_failure.factor,
                            first_failure.radius};
                        if (first_failure.level < 0 || candidate < current) {
                            first_failure = {
                                level,
                                code,
                                factor,
                                radius,
                                value,
                                profile,
                                square,
                                output};
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_SQUARE_BOUNDARY_INSERTION_EXTENDED"
            << " maximum_level=" << 2 * maximum_level
            << " maximum_coefficient=" << maximum_coefficient
            << " threads=" << omp_get_max_threads()
            << " profiles=" << profiles
            << " admissible_squares=" << admissible_squares
            << " insertion_checks=" << insertion_checks
            << " boundary_failures=" << boundary_failures
            << " first_failure=("
            << (first_failure.level < 0 ? -1 : 2 * first_failure.level)
            << ','
            << (first_failure.factor < 0 ? -1 : 2 * first_failure.factor)
            << ','
            << (first_failure.radius < 0 ? -1 : 2 * first_failure.radius)
            << ',' << first_failure.value << ')'
            << " first_profile=" << render(first_failure.profile)
            << " first_square=" << render(first_failure.square)
            << " first_output=" << render(first_failure.output)
            << " result="
            << (boundary_failures == 0U
                    ? "NO_SQUARE_BOUNDARY_INSERTION_FAILURE"
                    : "SQUARE_BOUNDARY_INSERTION_FAILURE")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
