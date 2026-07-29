#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(
            std::string(name) + " must be a positive integer"
        );
    }
    return static_cast<int>(value);
}

bool fuses(int level, int first, int second, int output) {
    return std::abs(first - second) <= output
        && output <= std::min(
            first + second, 2 * level - first - second
        )
        && ((first + second + output) & 1) == 0;
}

void apply_additive_compound(
    int level,
    int label,
    std::vector<Integer>& state
) {
    const int width = level + 1;
    const auto index = [width](int first, int second) {
        return static_cast<std::size_t>(first * width + second);
    };
    std::vector<Integer> next(
        static_cast<std::size_t>(width * width)
    );
    for (int first = 0; first <= level; ++first) {
        for (int second = 0; second <= level; ++second) {
            const Integer& coefficient = state[index(first, second)];
            if (coefficient == 0) {
                continue;
            }
            for (int output = 0; output <= level; ++output) {
                if (fuses(level, label, first, output)) {
                    next[index(output, second)] += coefficient;
                }
                if (fuses(level, label, second, output)) {
                    next[index(first, output)] += coefficient;
                }
            }
        }
    }
    state.swap(next);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_terminal_wedge_square "
                "MAXIMUM_LEVEL MAXIMUM_PAIR_STEPS"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_steps =
            parse_positive(argv[2], "maximum pair steps");
        if (maximum_level < 3) {
            throw std::runtime_error("require maximum level>=3");
        }

        std::size_t parameter_rows = 0;
        std::size_t state_rows = 0;
        std::size_t coefficient_rows = 0;
        std::size_t negative_interior_coefficients = 0;
        bool interior_witness_present = false;
        int witness_level = 0;
        int witness_label = 0;
        int witness_step = 0;
        int witness_first = 0;
        int witness_second = 0;
        Integer witness_coefficient = 0;
        for (int level = 4; level <= maximum_level; level += 2) {
            const int width = level + 1;
            const auto index = [width](int first, int second) {
                return static_cast<std::size_t>(
                    first * width + second
                );
            };
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameter_rows;
                std::vector<Integer> state(
                    static_cast<std::size_t>(width * width)
                );
                state[index(label, 0)] = 1;
                state[index(0, label)] = -1;
                apply_additive_compound(level, label, state);

                for (int step = 0; step <= maximum_steps; ++step) {
                    ++state_rows;
                    for (int first = 1; first <= level; ++first) {
                        for (int second = 0;
                             second < first; ++second) {
                            ++coefficient_rows;
                            const Integer& coefficient =
                                state[index(first, second)];
                            if (coefficient < 0) {
                                if (second == 0) {
                                    std::cout
                                        << "SU2_TERMINAL_WEDGE_SQUARE"
                                        << " boundary_result=FAIL"
                                        << " level=" << level
                                        << " label=" << label
                                        << " pair_step=" << step
                                        << " power=" << 2 * step + 1
                                        << " wedge=(" << first
                                        << ',' << second << ')'
                                        << " coefficient=" << coefficient
                                        << '\n';
                                    return EXIT_FAILURE;
                                }
                                ++negative_interior_coefficients;
                                if (!interior_witness_present) {
                                    interior_witness_present = true;
                                    witness_level = level;
                                    witness_label = label;
                                    witness_step = step;
                                    witness_first = first;
                                    witness_second = second;
                                    witness_coefficient = coefficient;
                                }
                            }
                            if (
                                state[index(second, first)]
                                != -coefficient
                            ) {
                                throw std::runtime_error(
                                    "antisymmetry invariant failed"
                                );
                            }
                        }
                    }
                    if (step != maximum_steps) {
                        apply_additive_compound(
                            level, label, state
                        );
                        apply_additive_compound(
                            level, label, state
                        );
                    }
                }
            }
        }

        std::cout
            << "SU2_TERMINAL_WEDGE_SQUARE"
            << " maximum_level=" << maximum_level
            << " maximum_pair_steps=" << maximum_steps
            << " parameter_rows=" << parameter_rows
            << " state_rows=" << state_rows
            << " coefficient_rows=" << coefficient_rows
            << " negative_boundary_coefficients=0"
            << " negative_interior_coefficients="
            << negative_interior_coefficients;
        if (interior_witness_present) {
            std::cout
                << " first_negative_interior={level="
                << witness_level
                << ",label=" << witness_label
                << ",pair_step=" << witness_step
                << ",power=" << 2 * witness_step + 1
                << ",wedge=(" << witness_first
                << ',' << witness_second << ')'
                << ",coefficient=" << witness_coefficient
                << '}';
        }
        std::cout
            << " boundary_result=PASS_DISCOVERY"
            << " full_oriented_cone="
            << (
                negative_interior_coefficients == 0
                    ? "PASS_DISCOVERY" : "FAIL"
            )
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_WEDGE_SQUARE FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
