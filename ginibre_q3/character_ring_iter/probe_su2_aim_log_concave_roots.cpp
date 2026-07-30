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
using Matrix = std::vector<std::vector<long long>>;

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

bool log_concave_interval(const Vector& values) {
    bool seen = false;
    bool ended = false;
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (values[index] > 0) {
            if (ended) {
                return false;
            }
            seen = true;
        } else if (seen) {
            ended = true;
        }
        if (index > 0U && index + 1U < values.size()
            && values[index] * values[index]
                < values[index - 1U] * values[index + 1U]) {
            return false;
        }
    }
    return seen;
}

Vector ordinary_square(const Vector& root) {
    const int maximum = static_cast<int>(root.size()) - 1;
    Vector square(static_cast<std::size_t>(2 * maximum + 1), 0);
    for (int left = 0; left <= maximum; ++left) {
        for (int right = 0; right <= maximum; ++right) {
            const cpp_int weight
                = root[static_cast<std::size_t>(left)]
                  * root[static_cast<std::size_t>(right)];
            if (weight == 0) {
                continue;
            }
            for (int target = std::abs(left - right);
                 target <= left + right;
                 ++target) {
                square[static_cast<std::size_t>(target)] += weight;
            }
        }
    }
    return square;
}

cpp_int at(const Vector& values, int index) {
    if (index < 0 || index >= static_cast<int>(values.size())) {
        return 0;
    }
    return values[static_cast<std::size_t>(index)];
}

Vector image_suffix(
    const Vector& square,
    int level,
    int shell) {
    const int period = 2 * level + 2;
    Vector suffix(static_cast<std::size_t>(level + 1), 0);
    const int maximum = static_cast<int>(square.size()) - 1;
    for (int image = shell; image * period <= maximum + level + 1; ++image) {
        for (int radius = 0; radius <= level; ++radius) {
            suffix[static_cast<std::size_t>(radius)]
                += at(square, image * period + radius)
                  - at(square, (image + 1) * period - radius - 1);
        }
    }
    return suffix;
}

Matrix fusion_matrix(int level, int factor) {
    Matrix matrix(
        static_cast<std::size_t>(level + 1),
        std::vector<long long>(
            static_cast<std::size_t>(level + 1),
            0));
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
        std::vector<long long>(left.size(), 0));
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
        std::vector<long long>(
            static_cast<std::size_t>(level + 1),
            0));
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
    int level = -1;
    int shell = -1;
    int radius = -1;
    cpp_int value = 0;
    Vector root;
    Vector square;
    Vector suffix;
};

struct PaymentFailure {
    int level = -1;
    int shell = -1;
    int factor = -1;
    int radius = -1;
    cpp_int value = 0;
    cpp_int boundary = 0;
    Vector root;
    Vector square;
    Vector reserve;
};

void record(
    Failure& failure,
    int level,
    int shell,
    int radius,
    const cpp_int& value,
    const Vector& root,
    const Vector& square,
    const Vector& suffix) {
    if (failure.level >= 0) {
        return;
    }
    failure = {level, shell, radius, value, root, square, suffix};
}

std::string render_failure(const Failure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " shell=" + std::to_string(failure.shell)
        + " radius=" + std::to_string(failure.radius)
        + " value=" + failure.value.convert_to<std::string>()
        + " root=" + render(failure.root)
        + " square=" + render(failure.square)
        + " suffix=" + render(failure.suffix);
}

std::string render_failure(const PaymentFailure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " shell=" + std::to_string(failure.shell)
        + " factor=" + std::to_string(failure.factor)
        + " radius=" + std::to_string(failure.radius)
        + " value=" + failure.value.convert_to<std::string>()
        + " boundary=" + failure.boundary.convert_to<std::string>()
        + " root=" + render(failure.root)
        + " square=" + render(failure.square)
        + " reserve=" + render(failure.reserve);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_index = argc >= 2
            ? positive_argument(argv[1], "maximum_index")
            : 8;
        const int maximum_coefficient = argc >= 3
            ? positive_argument(argv[2], "maximum_coefficient")
            : 5;
        const int maximum_level = argc >= 4
            ? positive_argument(argv[3], "maximum_half_level")
            : 8;
        if (argc > 4 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: probe_su2_aim_log_concave_roots "
                "[maximum_index] [maximum_coefficient] "
                "[maximum_half_level]");
        }

        const std::uint64_t base
            = static_cast<std::uint64_t>(maximum_coefficient + 1);
        std::uint64_t profile_count = 1U;
        for (int index = 0; index <= maximum_index; ++index) {
            if (profile_count
                > std::numeric_limits<std::uint64_t>::max() / base) {
                throw std::overflow_error("profile count overflow");
            }
            profile_count *= base;
        }

        std::uint64_t profiles = 0U;
        std::uint64_t log_concave_roots = 0U;
        std::uint64_t suffix_profiles = 0U;
        std::uint64_t coordinate_checks = 0U;
        std::uint64_t coordinate_failures = 0U;
        std::uint64_t monotonicity_checks = 0U;
        std::uint64_t monotonicity_failures = 0U;
        std::uint64_t log_concavity_checks = 0U;
        std::uint64_t log_concavity_failures = 0U;
        std::uint64_t image_pair_margin_checks = 0U;
        std::uint64_t image_pair_margin_failures = 0U;
        std::uint64_t payment_checks = 0U;
        std::uint64_t payment_failures = 0U;
        Failure first_image_pair_margin;
        Failure first_coordinate;
        Failure first_monotonicity;
        Failure first_log_concavity;
        PaymentFailure first_payment;

        for (std::uint64_t code = 1U; code < profile_count; ++code) {
            ++profiles;
            std::uint64_t remainder = code;
            Vector root(static_cast<std::size_t>(maximum_index + 1), 0);
            for (int index = 0; index <= maximum_index; ++index) {
                root[static_cast<std::size_t>(index)] = remainder % base;
                remainder /= base;
            }
            if (!log_concave_interval(root)) {
                continue;
            }
            ++log_concave_roots;
            const Vector square = ordinary_square(root);
            const int maximum_square
                = static_cast<int>(square.size()) - 1;
            for (int level = 2; level <= maximum_level; ++level) {
                const int period = 2 * level + 2;
                std::vector<Matrix> transforms(
                    static_cast<std::size_t>(level / 2 + 1));
                for (int factor = 1;
                     factor <= level / 2;
                     ++factor) {
                    transforms[static_cast<std::size_t>(factor)]
                        = reserve_transform(level, factor);
                }
                for (int shell = 1;
                     shell * period <= maximum_square + level + 1;
                     ++shell) {
                    const Vector suffix
                        = image_suffix(square, level, shell);
                    ++suffix_profiles;
                    for (int radius = 0; radius <= level; ++radius) {
                        ++coordinate_checks;
                        if (suffix[static_cast<std::size_t>(radius)] < 0) {
                            ++coordinate_failures;
                            record(
                                first_coordinate,
                                level,
                                shell,
                                radius,
                                suffix[static_cast<std::size_t>(radius)],
                                root,
                                square,
                                suffix);
                        }
                    }
                    for (int radius = 0; radius < level; ++radius) {
                        ++monotonicity_checks;
                        const cpp_int margin
                            = suffix[static_cast<std::size_t>(radius)]
                              - suffix[
                                  static_cast<std::size_t>(radius + 1)];
                        if (margin < 0) {
                            ++monotonicity_failures;
                            record(
                                first_monotonicity,
                                level,
                                shell,
                                radius,
                                margin,
                                root,
                                square,
                                suffix);
                        }
                        for (int image = shell;
                             image * period <= maximum_square + level + 1;
                             ++image) {
                            const int left = image * period + radius;
                            const int right
                                = (image + 1) * period - radius - 2;
                            const cpp_int image_margin
                                = at(square, left)
                                  - at(square, left + 1)
                                  + at(square, right)
                                  - at(square, right + 1);
                            ++image_pair_margin_checks;
                            if (image_margin < 0) {
                                ++image_pair_margin_failures;
                                record(
                                    first_image_pair_margin,
                                    level,
                                    image,
                                    radius,
                                    image_margin,
                                    root,
                                    square,
                                    suffix);
                            }
                        }
                    }
                    for (int radius = 1; radius < level; ++radius) {
                        ++log_concavity_checks;
                        const cpp_int margin
                            = suffix[static_cast<std::size_t>(radius)]
                                * suffix[static_cast<std::size_t>(radius)]
                              - suffix[
                                    static_cast<std::size_t>(radius - 1)]
                                * suffix[
                                    static_cast<std::size_t>(radius + 1)];
                        if (margin < 0) {
                            ++log_concavity_failures;
                            record(
                                first_log_concavity,
                                level,
                                shell,
                                radius,
                                margin,
                                root,
                                square,
                                suffix);
                        }
                    }
                    Vector reserve(
                        static_cast<std::size_t>(level + 1),
                        0);
                    for (int radius = 0; radius < level; ++radius) {
                        reserve[static_cast<std::size_t>(radius)]
                            = suffix[static_cast<std::size_t>(radius)]
                              - suffix[
                                  static_cast<std::size_t>(radius + 1)];
                    }
                    reserve[static_cast<std::size_t>(level)]
                        = suffix[static_cast<std::size_t>(level)];
                    if (std::any_of(
                            reserve.begin(),
                            reserve.end(),
                            [](const cpp_int& value) {
                                return value < 0;
                            })) {
                        continue;
                    }
                    for (int factor = 1;
                         factor <= level / 2;
                         ++factor) {
                        const Matrix& transform
                            = transforms[
                                static_cast<std::size_t>(factor)];
                        for (int radius = 0;
                             radius <= level;
                             ++radius) {
                            cpp_int boundary = 0;
                            if (radius < 2 * factor) {
                                for (int index
                                         = shell * period
                                           - 2 * factor + radius;
                                     index
                                         <= shell * period
                                            + 2 * factor - radius - 1;
                                     ++index) {
                                    boundary += at(square, index);
                                }
                            }
                            cpp_int payment = boundary;
                            for (int source = 0;
                                 source <= level;
                                 ++source) {
                                const long long coefficient
                                    = transform[
                                        static_cast<std::size_t>(
                                            radius)]
                                        [static_cast<std::size_t>(
                                            source)];
                                if (coefficient < 0) {
                                    payment += coefficient
                                        * reserve[
                                            static_cast<std::size_t>(
                                                source)];
                                }
                            }
                            ++payment_checks;
                            if (payment < 0) {
                                ++payment_failures;
                                if (first_payment.level < 0) {
                                    first_payment = {
                                        level,
                                        shell,
                                        factor,
                                        radius,
                                        payment,
                                        boundary,
                                        root,
                                        square,
                                        reserve};
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_AIM_LOG_CONCAVE_ROOTS"
            << " maximum_index=" << maximum_index
            << " maximum_coefficient=" << maximum_coefficient
            << " maximum_level=" << 2 * maximum_level
            << " profiles=" << profiles
            << " log_concave_roots=" << log_concave_roots
            << " suffix_profiles=" << suffix_profiles
            << " coordinate_checks=" << coordinate_checks
            << " coordinate_failures=" << coordinate_failures
            << " monotonicity_checks=" << monotonicity_checks
            << " monotonicity_failures=" << monotonicity_failures
            << " log_concavity_checks=" << log_concavity_checks
            << " log_concavity_failures=" << log_concavity_failures
            << " image_pair_margin_checks=" << image_pair_margin_checks
            << " image_pair_margin_failures="
            << image_pair_margin_failures
            << " payment_checks=" << payment_checks
            << " payment_failures=" << payment_failures
            << '\n'
            << "FIRST_COORDINATE_FAILURE "
            << render_failure(first_coordinate) << '\n'
            << "FIRST_MONOTONICITY_FAILURE "
            << render_failure(first_monotonicity) << '\n'
            << "FIRST_LOG_CONCAVITY_FAILURE "
            << render_failure(first_log_concavity) << '\n'
            << "FIRST_IMAGE_PAIR_MARGIN_FAILURE "
            << render_failure(first_image_pair_margin) << '\n'
            << "FIRST_PAYMENT_FAILURE "
            << render_failure(first_payment) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
