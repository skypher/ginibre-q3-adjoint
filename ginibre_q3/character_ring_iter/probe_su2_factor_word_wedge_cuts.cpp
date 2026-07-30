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

std::uint64_t positive_argument(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const unsigned long long parsed = std::stoull(value, &consumed, 10);
    if (consumed != value.size() || parsed == 0U) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<std::uint64_t>(parsed);
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

Vector fusion_step(const Vector& input, int level, int factor) {
    Vector output(static_cast<std::size_t>(level + 1), 0);
    for (int source = 0; source <= level; ++source) {
        const cpp_int& coefficient
            = input[static_cast<std::size_t>(source)];
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

std::string render(const std::vector<int>& values) {
    std::string result = "[";
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            result += ',';
        }
        result += std::to_string(values[index]);
    }
    return result + ']';
}

struct Failure {
    int level = -1;
    int radius = -1;
    int cut = -1;
    cpp_int value = 0;
    std::vector<int> word;
    Vector profile;
    Vector square;
};

void record(
    Failure& failure,
    int level,
    int radius,
    int cut,
    const cpp_int& value,
    const std::vector<int>& word,
    const Vector& profile,
    const Vector& square) {
    if (failure.level >= 0) {
        return;
    }
    failure = {level, radius, cut, value, word, profile, square};
}

std::string render_failure(const Failure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " radius=" + std::to_string(failure.radius)
        + " cut=" + std::to_string(failure.cut)
        + " value=" + failure.value.convert_to<std::string>()
        + " word=" + render(failure.word)
        + " profile=" + render(failure.profile)
        + " square=" + render(failure.square);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::uint64_t samples = argc >= 2
            ? positive_argument(argv[1], "samples")
            : UINT64_C(5000);
        const int maximum_level = argc >= 3
            ? static_cast<int>(positive_argument(argv[2], "maximum_level"))
            : 24;
        const int maximum_length = argc >= 4
            ? static_cast<int>(positive_argument(argv[3], "maximum_length"))
            : 20;
        if (argc > 4 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: probe_su2_factor_word_wedge_cuts "
                "[samples] [maximum_half_level] [maximum_word_length]");
        }

        std::uint64_t wedge_states = 0U;
        std::uint64_t star_checks = 0U;
        std::uint64_t star_failures = 0U;
        std::uint64_t cut_flux_checks = 0U;
        std::uint64_t cut_flux_failures = 0U;
        std::uint64_t central_window_checks = 0U;
        std::uint64_t central_window_failures = 0U;
        std::uint64_t star_prefix_checks = 0U;
        std::uint64_t star_prefix_failures = 0U;
        std::uint64_t star_suffix_checks = 0U;
        std::uint64_t star_suffix_failures = 0U;
        Failure first_star;
        Failure first_cut_flux;
        Failure first_central_window;
        Failure first_star_prefix;
        Failure first_star_suffix;

        for (std::uint64_t sample = 0U; sample < samples; ++sample) {
            std::uint64_t state
                = sample ^ UINT64_C(0x13198a2e03707344);
            const int level = 2 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_level - 1));
            const int length = 1 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_length));
            Vector profile(static_cast<std::size_t>(level + 1), 0);
            profile[0] = 1;
            std::vector<int> word;
            word.reserve(static_cast<std::size_t>(length));
            for (int position = 0; position < length; ++position) {
                const int factor = 1 + static_cast<int>(
                    splitmix64(state)
                    % static_cast<std::uint64_t>(level));
                word.push_back(factor);
                profile = fusion_step(profile, level, factor);
            }

            Vector square(static_cast<std::size_t>(level + 1), 0);
            for (int factor = 0; factor <= level; ++factor) {
                const Vector translated
                    = fusion_step(profile, level, factor);
                for (int target = 0; target <= level; ++target) {
                    square[static_cast<std::size_t>(target)]
                        += profile[static_cast<std::size_t>(factor)]
                          * translated[static_cast<std::size_t>(target)];
                }
            }
            const Matrix kernel = multiplication_matrix(square);

            for (int radius = 1; radius <= level; ++radius) {
                ++wedge_states;
                Matrix wedge(
                    static_cast<std::size_t>(level + 1),
                    Vector(static_cast<std::size_t>(level + 1), 0));
                for (int left = 0; left < level; ++left) {
                    for (int right = left + 1; right <= level; ++right) {
                        wedge[static_cast<std::size_t>(left)]
                             [static_cast<std::size_t>(right)]
                            = square[static_cast<std::size_t>(left)]
                                * kernel[static_cast<std::size_t>(right)]
                                        [static_cast<std::size_t>(radius)]
                              - square[static_cast<std::size_t>(right)]
                                * kernel[static_cast<std::size_t>(left)]
                                        [static_cast<std::size_t>(radius)];
                    }
                }

                for (int target = 1; target <= level; ++target) {
                    ++star_checks;
                    const cpp_int& value = wedge[0][
                        static_cast<std::size_t>(target)];
                    if (value < 0) {
                        ++star_failures;
                        record(
                            first_star,
                            level,
                            radius,
                            target,
                            value,
                            word,
                            profile,
                            square);
                    }
                }

                for (int cut = 0; cut < level; ++cut) {
                    cpp_int value = 0;
                    for (int left = 0; left <= cut; ++left) {
                        for (int right = cut + 1; right <= level; ++right) {
                            value += wedge[static_cast<std::size_t>(left)]
                                          [static_cast<std::size_t>(right)];
                        }
                    }
                    ++cut_flux_checks;
                    if (value < 0) {
                        ++cut_flux_failures;
                        record(
                            first_cut_flux,
                            level,
                            radius,
                            cut,
                            value,
                            word,
                            profile,
                            square);
                    }
                }

                for (int depth = 0; 2 * depth < level; ++depth) {
                    cpp_int value = 0;
                    for (int left = depth; left < level - depth; ++left) {
                        for (int right = left + 1;
                             right <= level - depth;
                             ++right) {
                            value += wedge[static_cast<std::size_t>(left)]
                                          [static_cast<std::size_t>(right)];
                        }
                    }
                    ++central_window_checks;
                    if (value < 0) {
                        ++central_window_failures;
                        record(
                            first_central_window,
                            level,
                            radius,
                            depth,
                            value,
                            word,
                            profile,
                            square);
                    }
                }

                cpp_int star_prefix = 0;
                for (int target = 1; target <= level; ++target) {
                    star_prefix += wedge[0][
                        static_cast<std::size_t>(target)];
                    ++star_prefix_checks;
                    if (star_prefix < 0) {
                        ++star_prefix_failures;
                        record(
                            first_star_prefix,
                            level,
                            radius,
                            target,
                            star_prefix,
                            word,
                            profile,
                            square);
                    }
                }

                cpp_int star_suffix = 0;
                for (int target = level; target >= 1; --target) {
                    star_suffix += wedge[0][
                        static_cast<std::size_t>(target)];
                    ++star_suffix_checks;
                    if (star_suffix < 0) {
                        ++star_suffix_failures;
                        record(
                            first_star_suffix,
                            level,
                            radius,
                            target,
                            star_suffix,
                            word,
                            profile,
                            square);
                    }
                }
            }
        }

        std::cout
            << "SU2_FACTOR_WORD_WEDGE_CUTS"
            << " samples=" << samples
            << " maximum_level=" << 2 * maximum_level
            << " maximum_word_length=" << maximum_length
            << " wedge_states=" << wedge_states
            << " star_checks=" << star_checks
            << " star_failures=" << star_failures
            << " cut_flux_checks=" << cut_flux_checks
            << " cut_flux_failures=" << cut_flux_failures
            << " central_window_checks=" << central_window_checks
            << " central_window_failures=" << central_window_failures
            << " star_prefix_checks=" << star_prefix_checks
            << " star_prefix_failures=" << star_prefix_failures
            << " star_suffix_checks=" << star_suffix_checks
            << " star_suffix_failures=" << star_suffix_failures
            << '\n'
            << "FIRST_STAR_FAILURE " << render_failure(first_star) << '\n'
            << "FIRST_CUT_FLUX_FAILURE "
            << render_failure(first_cut_flux) << '\n'
            << "FIRST_CENTRAL_WINDOW_FAILURE "
            << render_failure(first_central_window) << '\n'
            << "FIRST_STAR_PREFIX_FAILURE "
            << render_failure(first_star_prefix) << '\n'
            << "FIRST_STAR_SUFFIX_FAILURE "
            << render_failure(first_star_suffix) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
