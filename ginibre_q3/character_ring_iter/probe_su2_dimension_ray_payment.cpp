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

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::multiprecision::cpp_rational;
using Vector = std::vector<Integer>;
using Matrix = std::vector<Vector>;

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long parsed = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || parsed <= 0
        || parsed > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

std::uint64_t splitmix64(std::uint64_t& state) {
    state += UINT64_C(0x9e3779b97f4a7c15);
    std::uint64_t value = state;
    value = (value ^ (value >> 30U))
        * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U))
        * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

Vector triangular_step(const Vector& input, int radius, int factor) {
    const int next_radius = radius + 2 * factor;
    Vector output(static_cast<std::size_t>(2 * next_radius + 1));
    for (int source = -radius; source <= radius; ++source) {
        const Integer& weight = input[static_cast<std::size_t>(source + radius)];
        if (weight == 0) {
            continue;
        }
        for (int shift = -2 * factor; shift <= 2 * factor; ++shift) {
            output[static_cast<std::size_t>(source + shift + next_radius)] +=
                weight * (2 * factor + 1 - std::abs(shift));
        }
    }
    return output;
}

Integer line_value(const Vector& values, int radius, int index) {
    if (index < -radius || index > radius) {
        return 0;
    }
    return values[static_cast<std::size_t>(index + radius)];
}

Integer line_gradient(const Vector& values, int radius, int index) {
    return line_value(values, radius, index)
        - line_value(values, radius, index + 1);
}

Matrix multiplication_matrix(const Vector& values) {
    const int level = static_cast<int>(values.size()) - 1;
    Matrix result(values.size(), Vector(values.size()));
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            const int lower = std::abs(left - right);
            const int upper = std::min(
                left + right,
                2 * level - left - right
            );
            for (int label = lower; label <= upper; ++label) {
                result[static_cast<std::size_t>(left)]
                      [static_cast<std::size_t>(right)] +=
                    values[static_cast<std::size_t>(label)];
            }
        }
    }
    return result;
}

std::string render(const std::vector<int>& word) {
    std::string result{"["};
    for (std::size_t index = 0U; index < word.size(); ++index) {
        if (index != 0U) {
            result += ',';
        }
        result += std::to_string(word[index]);
    }
    return result + ']';
}

Integer current(
    const Vector& left,
    const Matrix& left_matrix,
    const Vector& right,
    const Matrix& right_matrix,
    int row,
    int column
) {
    return left[0] * right_matrix[static_cast<std::size_t>(row)]
                                   [static_cast<std::size_t>(column)]
        + right[0] * left_matrix[static_cast<std::size_t>(row)]
                                 [static_cast<std::size_t>(column)]
        - left[static_cast<std::size_t>(row)]
            * right[static_cast<std::size_t>(column)]
        - right[static_cast<std::size_t>(row)]
            * left[static_cast<std::size_t>(column)];
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2
            && std::string(argv[1]) == "--replay-ray-obstruction") {
            constexpr int level = 7;
            Vector ray(static_cast<std::size_t>(level + 1));
            ray[0] = 1;
            ray[1] = 3;
            const Matrix matrix = multiplication_matrix(ray);
            const Integer cross = current(ray, matrix, ray, matrix, 1, 1);
            if (cross != -10) {
                throw std::runtime_error("dimension-ray replay mismatch");
            }
            std::cout
                << "SU2_DIMENSION_RAY_OBSTRUCTION"
                << " level=7 ray=1 indices=(1,1)"
                << " anchored_current=" << cross / 2
                << " cross_current=" << cross
                << " result=PASS_EXACT\n";
            return EXIT_SUCCESS;
        }
        if (argc != 4) {
            throw std::runtime_error(
                "usage: probe_su2_dimension_ray_payment "
                "SAMPLES MAXIMUM_LEVEL MAXIMUM_LENGTH"
            );
        }
        const int samples = parse_positive(argv[1], "samples");
        const int maximum_level = parse_positive(argv[2], "maximum level");
        const int maximum_length = parse_positive(argv[3], "maximum length");
        if (maximum_level < 2) {
            throw std::runtime_error("maximum level must be at least two");
        }

        std::uint64_t profiles = 0U;
        std::uint64_t ratio_checks = 0U;
        std::uint64_t ratio_failures = 0U;
        std::uint64_t ray_cross_checks = 0U;
        std::uint64_t ray_cross_failures = 0U;
        std::uint64_t ray_current_checks = 0U;
        std::uint64_t ray_current_failures = 0U;
        std::uint64_t central_checks = 0U;
        std::uint64_t central_failures = 0U;
        bool printed_ratio = false;
        bool printed_cross = false;
        bool printed_ray = false;
        bool printed_central = false;

        for (int sample = 0; sample < samples; ++sample) {
            std::uint64_t state = static_cast<std::uint64_t>(sample)
                ^ UINT64_C(0x243f6a8885a308d3);
            const int level = 2 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_level - 1)
            );
            const int length = 1 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_length)
            );
            Vector line{Integer(1)};
            int radius = 0;
            std::vector<int> word;
            word.reserve(static_cast<std::size_t>(length));
            for (int position = 0; position < length; ++position) {
                const int factor = 1 + static_cast<int>(
                    splitmix64(state)
                    % static_cast<std::uint64_t>(level)
                );
                word.push_back(factor);
                const int reduced = std::min(factor, level - factor);
                line = triangular_step(line, radius, reduced);
                radius += 2 * reduced;
            }

            const int period = 2 * level + 2;
            const int shell_count = radius / period + 1;
            Vector central(static_cast<std::size_t>(level + 1));
            Vector proper(static_cast<std::size_t>(level + 1));
            for (int shell = 0; shell < shell_count; ++shell) {
                for (int index = 0; index <= level; ++index) {
                    const Integer value = line_gradient(
                        line,
                        radius,
                        shell * period + index
                    ) - line_gradient(
                        line,
                        radius,
                        (shell + 1) * period - index - 1
                    );
                    if (shell == 0) {
                        central[static_cast<std::size_t>(index)] = value;
                    } else {
                        proper[static_cast<std::size_t>(index)] += value;
                    }
                }
            }

            std::vector<Rational> coefficients(
                static_cast<std::size_t>(level + 1)
            );
            for (int cutoff = 0; cutoff <= level; ++cutoff) {
                const Integer numerator = proper[static_cast<std::size_t>(cutoff)]
                    * (2 * cutoff + 3)
                    - (cutoff == level
                           ? Integer(0)
                           : proper[static_cast<std::size_t>(cutoff + 1)]
                                 * (2 * cutoff + 1));
                const Integer denominator = (2 * cutoff + 1)
                    * (2 * cutoff + 3);
                ++ratio_checks;
                if (numerator < 0) {
                    ++ratio_failures;
                    if (!printed_ratio) {
                        printed_ratio = true;
                        std::cout
                            << "FIRST_DIMENSION_RAY_RATIO_FAILURE"
                            << " level=" << level
                            << " cutoff=" << cutoff
                            << " margin=" << numerator
                            << " word=" << render(word) << '\n';
                    }
                }
                coefficients[static_cast<std::size_t>(cutoff)] =
                    Rational(numerator) / Rational(denominator);
            }

            std::vector<Vector> rays;
            std::vector<Matrix> ray_matrices;
            rays.reserve(static_cast<std::size_t>(level + 1));
            ray_matrices.reserve(static_cast<std::size_t>(level + 1));
            for (int cutoff = 0; cutoff <= level; ++cutoff) {
                Vector ray(static_cast<std::size_t>(level + 1));
                for (int index = 0; index <= cutoff; ++index) {
                    ray[static_cast<std::size_t>(index)] = 2 * index + 1;
                }
                ray_matrices.push_back(multiplication_matrix(ray));
                rays.push_back(std::move(ray));
            }

            const Matrix central_matrix = multiplication_matrix(central);
            for (int left = 0; left <= level; ++left) {
                for (int right = left; right <= level; ++right) {
                    for (int row = 0; row <= level; ++row) {
                        for (int column = 0; column <= level; ++column) {
                            const Integer cross = current(
                                rays[static_cast<std::size_t>(left)],
                                ray_matrices[static_cast<std::size_t>(left)],
                                rays[static_cast<std::size_t>(right)],
                                ray_matrices[static_cast<std::size_t>(right)],
                                row,
                                column
                            );
                            ++ray_cross_checks;
                            if (cross < 0) {
                                ++ray_cross_failures;
                                if (!printed_cross) {
                                    printed_cross = true;
                                    std::cout
                                        << "FIRST_DIMENSION_RAY_CROSS_FAILURE"
                                        << " level=" << level
                                        << " rays=(" << left << ',' << right
                                        << ") indices=(" << row << ',' << column
                                        << ") value=" << cross << '\n';
                                }
                            }
                        }
                    }
                }
            }

            for (int row = 0; row <= level; ++row) {
                for (int column = 0; column <= level; ++column) {
                    Rational lower = current(
                        central,
                        central_matrix,
                        central,
                        central_matrix,
                        row,
                        column
                    ) / 2;
                    for (int cutoff = 0; cutoff <= level; ++cutoff) {
                        const Rational& coefficient =
                            coefficients[static_cast<std::size_t>(cutoff)];
                        if (coefficient == 0) {
                            continue;
                        }
                        const Integer cross = current(
                            central,
                            central_matrix,
                            rays[static_cast<std::size_t>(cutoff)],
                            ray_matrices[static_cast<std::size_t>(cutoff)],
                            row,
                            column
                        );
                        const Integer self = current(
                            rays[static_cast<std::size_t>(cutoff)],
                            ray_matrices[static_cast<std::size_t>(cutoff)],
                            rays[static_cast<std::size_t>(cutoff)],
                            ray_matrices[static_cast<std::size_t>(cutoff)],
                            row,
                            column
                        ) / 2;
                        ++ray_current_checks;
                        if (self < 0) {
                            ++ray_current_failures;
                            if (!printed_ray) {
                                printed_ray = true;
                                std::cout
                                    << "FIRST_DIMENSION_RAY_CURRENT_FAILURE"
                                    << " level=" << level
                                    << " ray=" << cutoff
                                    << " indices=(" << row << ',' << column
                                    << ") value=" << self << '\n';
                            }
                        }
                        lower += coefficient * Rational(cross)
                            + coefficient * coefficient * Rational(self);
                    }
                    ++central_checks;
                    if (lower < 0) {
                        ++central_failures;
                        if (!printed_central) {
                            printed_central = true;
                            std::cout
                                << "FIRST_DIMENSION_RAY_PAYMENT_FAILURE"
                                << " level=" << level
                                << " indices=(" << row << ',' << column
                                << ") value=" << lower
                                << " word=" << render(word) << '\n';
                        }
                    }
                }
            }
            ++profiles;
        }

        std::cout
            << "SU2_DIMENSION_RAY_PAYMENT"
            << " samples=" << samples
            << " maximum_level=" << maximum_level
            << " maximum_length=" << maximum_length
            << " profiles=" << profiles
            << " ratio_checks=" << ratio_checks
            << " ratio_failures=" << ratio_failures
            << " ray_cross_checks=" << ray_cross_checks
            << " ray_cross_failures=" << ray_cross_failures
            << " ray_current_checks=" << ray_current_checks
            << " ray_current_failures=" << ray_current_failures
            << " central_checks=" << central_checks
            << " central_failures=" << central_failures
            << " result="
            << (ratio_failures == 0U && ray_cross_failures == 0U
                    && ray_current_failures == 0U && central_failures == 0U
                    ? "PASS_EXACT_SAMPLE"
                    : "COUNTEREXAMPLE")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
