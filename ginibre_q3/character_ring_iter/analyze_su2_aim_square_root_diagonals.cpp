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

bool selected_slope(
    int index,
    int level,
    int shell,
    int label) {
    if (index < 0) {
        return false;
    }
    const int period = 2 * level + 2;
    if (label == level) {
        return index >= shell * period + level
            && index % period == level;
    }
    if (index < shell * period + label) {
        return false;
    }
    const int residue = index % period;
    return residue == label || residue == period - label - 2;
}

long long reserve_coefficient(
    int index,
    int level,
    int shell,
    int label) {
    return static_cast<long long>(
               selected_slope(index, level, shell, label))
        - static_cast<long long>(
               selected_slope(index - 1, level, shell, label));
}

long long payment_coefficient(
    int index,
    int level,
    int shell,
    int factor,
    int radius,
    const Vector& transform_row) {
    const int period = 2 * level + 2;
    long long value = radius < 2 * factor
            && shell * period - 2 * factor + radius <= index
            && index
                <= shell * period + 2 * factor - radius - 1
        ? 1
        : 0;
    for (int source = 0; source <= level; ++source) {
        const long long coefficient
            = transform_row[static_cast<std::size_t>(source)];
        if (coefficient < 0) {
            value += coefficient
                * reserve_coefficient(
                    index,
                    level,
                    shell,
                    source);
        }
    }
    return value;
}

struct Failure {
    int level = -1;
    int shell = -1;
    int factor = -1;
    int radius = -1;
    int total = -1;
    int pair_index = -1;
    long long value = 0;
    Vector coefficients;
    Vector reserve_suffixes;
};

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
        + " shell=" + std::to_string(failure.shell)
        + " factor=" + std::to_string(failure.factor)
        + " radius=" + std::to_string(failure.radius)
        + " total=" + std::to_string(failure.total)
        + " pair_index=" + std::to_string(failure.pair_index)
        + " value=" + std::to_string(failure.value)
        + " coefficients=" + render(failure.coefficients)
        + " reserve_suffixes=" + render(failure.reserve_suffixes);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 60;
        const int maximum_shell = argc >= 3
            ? positive_argument(argv[2], "maximum_shell")
            : 4;
        const int maximum_periods = argc >= 4
            ? positive_argument(argv[3], "maximum_periods")
            : 12;
        if (argc > 4 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim_square_root_diagonals "
                "[maximum_half_level] [maximum_shell] "
                "[maximum_periods]");
        }

        std::uint64_t payments = 0U;
        std::uint64_t diagonals = 0U;
        std::uint64_t suffix_checks = 0U;
        std::uint64_t suffix_failures = 0U;
        Failure first_failure;

        for (int level = 2; level <= maximum_level; ++level) {
            const int period = 2 * level + 2;
            const int maximum_total = maximum_periods * period;
            for (int factor = 1;
                 factor <= level / 2;
                 ++factor) {
                const Matrix transform
                    = reserve_transform(level, factor);
                for (int radius = 0; radius <= level; ++radius) {
                    for (int shell = 1;
                         shell <= maximum_shell;
                         ++shell) {
                        ++payments;
                        Vector eta(
                            static_cast<std::size_t>(
                                maximum_total + 1),
                            0);
                        for (int index = 0;
                             index <= maximum_total;
                             ++index) {
                            eta[static_cast<std::size_t>(index)]
                                = payment_coefficient(
                                    index,
                                    level,
                                    shell,
                                    factor,
                                    radius,
                                    transform[
                                        static_cast<std::size_t>(
                                            radius)]);
                        }
                        Vector eta_prefix(
                            eta.size(),
                            0);
                        long long cumulative = 0;
                        for (int index = 0;
                             index <= maximum_total;
                             ++index) {
                            cumulative += eta[
                                static_cast<std::size_t>(index)];
                            eta_prefix[
                                static_cast<std::size_t>(index)]
                                = cumulative;
                        }
                        for (int total = 0;
                             total <= maximum_total;
                             ++total) {
                            ++diagonals;
                            const int middle = total / 2;
                            Vector coefficients(
                                static_cast<std::size_t>(
                                    middle + 1),
                                0);
                            for (int left = 0;
                                 left <= middle;
                                 ++left) {
                                const int right = total - left;
                                const int gap = right - left;
                                const long long kernel
                                    = eta_prefix[
                                          static_cast<std::size_t>(
                                              total)]
                                      - (gap == 0
                                             ? 0
                                             : eta_prefix[
                                                   static_cast<
                                                       std::size_t>(
                                                       gap - 1)]);
                                coefficients[
                                    static_cast<std::size_t>(left)]
                                    = left == right
                                    ? kernel
                                    : 2 * kernel;
                            }
                            long long suffix = 0;
                            for (int left = middle;
                                 left >= 0;
                                 --left) {
                                suffix += coefficients[
                                    static_cast<std::size_t>(left)];
                                ++suffix_checks;
                                if (suffix >= 0) {
                                    continue;
                                }
                                ++suffix_failures;
                                if (first_failure.level < 0) {
                                    Vector reserve_suffixes(
                                        static_cast<std::size_t>(
                                            level + 1),
                                        0);
                                    for (int source = 0;
                                         source <= level;
                                         ++source) {
                                        long long reserve_suffix = 0;
                                        for (int pair = left;
                                             pair <= middle;
                                             ++pair) {
                                            const int right
                                                = total - pair;
                                            const int gap
                                                = right - pair;
                                            long long kernel = 0;
                                            for (int index = gap;
                                                 index <= total;
                                                 ++index) {
                                                kernel
                                                    += reserve_coefficient(
                                                        index,
                                                        level,
                                                        shell,
                                                        source);
                                            }
                                            reserve_suffix
                                                += pair == right
                                                ? kernel
                                                : 2 * kernel;
                                        }
                                        reserve_suffixes[
                                            static_cast<std::size_t>(
                                                source)]
                                            = reserve_suffix;
                                    }
                                    first_failure = {
                                        level,
                                        shell,
                                        factor,
                                        radius,
                                        total,
                                        left,
                                        suffix,
                                        coefficients,
                                        reserve_suffixes};
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_AIM_SQUARE_ROOT_DIAGONALS"
            << " maximum_level=" << 2 * maximum_level
            << " maximum_shell=" << maximum_shell
            << " maximum_periods=" << maximum_periods
            << " adjustment=none"
            << " payments=" << payments
            << " diagonals=" << diagonals
            << " suffix_checks=" << suffix_checks
            << " suffix_failures=" << suffix_failures
            << '\n'
            << "FIRST_SUFFIX_FAILURE "
            << render_failure(first_failure)
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
