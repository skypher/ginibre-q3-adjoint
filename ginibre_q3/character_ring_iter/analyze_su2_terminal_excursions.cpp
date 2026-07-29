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

void retain_interior(int level, std::vector<Integer>& state) {
    const int width = level + 1;
    for (int label = 0; label <= level; ++label) {
        state[static_cast<std::size_t>(label * width)] = 0;
        state[static_cast<std::size_t>(label)] = 0;
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_terminal_excursions "
                "MAXIMUM_LEVEL MAXIMUM_INTERIOR_POWER"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_power =
            parse_positive(argv[2], "maximum interior power");
        if (maximum_level < 6) {
            throw std::runtime_error("require maximum level>=6");
        }

        std::size_t parameter_rows = 0;
        std::size_t kernel_rows = 0;
        std::size_t negative_kernel_entries = 0;
        bool negative_witness_present = false;
        int witness_level = 0;
        int witness_label = 0;
        int witness_power = 0;
        int witness_target = 0;
        Integer witness_coefficient = 0;
        for (int level = 6; level <= maximum_level; level += 2) {
            const int width = level + 1;
            const auto index = [width](int first, int second) {
                return static_cast<std::size_t>(
                    first * width + second
                );
            };
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameter_rows;
                std::vector<Integer> interior(
                    static_cast<std::size_t>(width * width)
                );
                for (int source = 2;
                     source <= 2 * label; source += 2) {
                    interior[index(source, 0)] = 1;
                    interior[index(0, source)] = -1;
                }
                apply_square(level, label, interior);
                retain_interior(level, interior);
                std::vector<Integer> cumulative(
                    static_cast<std::size_t>(width)
                );

                for (int power = 0;
                     power <= maximum_power; ++power) {
                    std::vector<Integer> next = interior;
                    apply_square(level, label, next);
                    for (int target = 2;
                         target <= level; target += 2) {
                        ++kernel_rows;
                        const Integer& coefficient =
                            next[index(target, 0)];
                        if (coefficient < 0) {
                            ++negative_kernel_entries;
                            if (!negative_witness_present) {
                                negative_witness_present = true;
                                witness_level = level;
                                witness_label = label;
                                witness_power = power;
                                witness_target = target;
                                witness_coefficient = coefficient;
                            }
                        }
                        Integer& prefix = cumulative[
                            static_cast<std::size_t>(target)
                        ];
                        prefix += coefficient;
                        if (prefix < 0) {
                            std::cout
                                << "SU2_TERMINAL_EXCURSIONS"
                                << " cumulative_result=FAIL"
                                << " level=" << level
                                << " label=" << label
                                << " maximum_interior_power=" << power
                                << " source=initial_boundary_seed"
                                << " target=" << target
                                << " cumulative=" << prefix
                                << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                    retain_interior(level, next);
                    interior.swap(next);
                }
            }
        }

        std::cout
            << "SU2_TERMINAL_EXCURSIONS"
            << " maximum_level=" << maximum_level
            << " maximum_interior_power=" << maximum_power
            << " parameter_rows=" << parameter_rows
            << " kernel_rows=" << kernel_rows
            << " negative_kernel_entries="
            << negative_kernel_entries;
        if (negative_witness_present) {
            std::cout
                << " first_negative_kernel={level="
                << witness_level
                << ",label=" << witness_label
                << ",interior_power=" << witness_power
                << ",target=" << witness_target
                << ",coefficient=" << witness_coefficient
                << '}';
        }
        std::cout
            << " cumulative_result=PASS_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_EXCURSIONS FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
