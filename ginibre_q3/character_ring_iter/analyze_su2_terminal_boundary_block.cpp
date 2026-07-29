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

void apply_square(
    int level,
    int label,
    std::vector<Integer>& state
) {
    apply_additive_compound(level, label, state);
    apply_additive_compound(level, label, state);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_terminal_boundary_block "
                "MAXIMUM_LEVEL MAXIMUM_PAIR_POWER"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_power =
            parse_positive(argv[2], "maximum pair power");
        if (maximum_level < 6) {
            throw std::runtime_error("require maximum level>=6");
        }

        std::size_t parameter_rows = 0;
        std::size_t matrix_rows = 0;
        std::size_t coefficient_rows = 0;
        for (int level = 6; level <= maximum_level; level += 2) {
            const int width = level + 1;
            const auto index = [width](int first, int second) {
                return static_cast<std::size_t>(
                    first * width + second
                );
            };
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameter_rows;
                for (int source = 1; source <= level; ++source) {
                    std::vector<Integer> state(
                        static_cast<std::size_t>(width * width)
                    );
                    state[index(source, 0)] = 1;
                    state[index(0, source)] = -1;
                    for (int power = 1;
                         power <= maximum_power; ++power) {
                        apply_square(level, label, state);
                        ++matrix_rows;
                        for (int target = 1;
                             target <= level; ++target) {
                            if ((target & 1) != (source & 1)) {
                                continue;
                            }
                            ++coefficient_rows;
                            const Integer& coefficient =
                                state[index(target, 0)];
                            if (coefficient < 0) {
                                std::cout
                                    << "SU2_TERMINAL_BOUNDARY_BLOCK"
                                    << " result=FAIL"
                                    << " level=" << level
                                    << " label=" << label
                                    << " pair_power=" << power
                                    << " source=" << source
                                    << " target=" << target
                                    << " coefficient=" << coefficient
                                    << '\n';
                                return EXIT_FAILURE;
                            }
                            if (
                                state[index(0, target)]
                                != -coefficient
                            ) {
                                throw std::runtime_error(
                                    "antisymmetry invariant failed"
                                );
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_TERMINAL_BOUNDARY_BLOCK"
            << " maximum_level=" << maximum_level
            << " maximum_pair_power=" << maximum_power
            << " parameter_rows=" << parameter_rows
            << " matrix_rows=" << matrix_rows
            << " coefficient_rows=" << coefficient_rows
            << " negative_coefficients=0"
            << " result=PASS_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_BOUNDARY_BLOCK FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
