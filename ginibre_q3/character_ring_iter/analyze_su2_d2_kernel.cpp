#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
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
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(value);
}

bool fuses(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= std::min(
            source + label,
            2 * level - source - label
        )
        && ((source + label + target) & 1) == 0;
}

void multiply(
    int level,
    int label,
    std::vector<Integer>& state
) {
    std::vector<Integer> next(static_cast<std::size_t>(level + 1));
    for (int source = 0; source <= level; ++source) {
        const Integer& coefficient =
            state[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        for (int target = 0; target <= level; ++target) {
            if (fuses(level, label, source, target)) {
                next[static_cast<std::size_t>(target)] += coefficient;
            }
        }
    }
    state.swap(next);
}

std::vector<Integer> power_column(
    int level,
    int label,
    int source,
    int power
) {
    std::vector<Integer> state(static_cast<std::size_t>(level + 1));
    state[static_cast<std::size_t>(source)] = 1;
    for (int step = 0; step < power; ++step) {
        multiply(level, label, state);
    }
    return state;
}

Integer integer_binomial(int top, int order) {
    if (top < order) {
        return 0;
    }
    Integer result = 1;
    for (int index = 0; index < order; ++index) {
        result *= top - index;
        result /= index + 1;
    }
    return result;
}

Integer ordinary_multiplicity(int power, int half_label, int half_target) {
    Integer result = 0;
    for (int image = 0; image <= power; ++image) {
        const int top =
            power * half_label - half_target
            - image * (2 * half_label + 1)
            + power - 2;
        const Integer term = integer_binomial(top, power - 2);
        if ((image & 1) == 0) {
            result += integer_binomial(power, image) * term;
        } else {
            result -= integer_binomial(power, image) * term;
        }
    }
    return result;
}

Integer affine_multiplicity(
    int power,
    int half_level,
    int half_label,
    int half_target
) {
    const int period = 2 * half_level + 2;
    return ordinary_multiplicity(power, half_label, half_target)
        - ordinary_multiplicity(
            power,
            half_label,
            period - half_target - 1
        )
        + ordinary_multiplicity(
            power,
            half_label,
            period + half_target
        );
}

}  // namespace

struct PowerResult {
    bool initialized = false;
    bool failed = false;
    Integer minimum = 0;
    int level = 0;
    int label = 0;
    int vertex = 0;
};

int main(int argc, char** argv) {
    try {
        if (argc == 4 && std::string(argv[1]) == "--case") {
            const int level = parse_positive(argv[2], "level");
            const int label = parse_positive(argv[3], "label");
            if ((level & 1) != 0 || (label & 1) != 0
                || 2 * label >= level) {
                throw std::runtime_error(
                    "case requires even level/label and 2*label<level"
                );
            }
            const std::vector<Integer> fourth =
                power_column(level, label, 0, 4);
            const std::vector<Integer> fifth =
                power_column(level, label, 0, 5);
            const Integer f4 = fourth[0];
            const Integer f5 = fifth[0];
            std::cout
                << "SU2_D2_KERNEL_CASE"
                << " level=" << level
                << " label=" << label
                << " f4=" << f4
                << " f5=" << f5 << '\n';
            for (int vertex = 0; vertex <= level; vertex += 2) {
                const Integer margin =
                    f4 * fifth[static_cast<std::size_t>(vertex)]
                    - f5 * fourth[static_cast<std::size_t>(vertex)];
                if (fourth[static_cast<std::size_t>(vertex)] != 0
                    || fifth[static_cast<std::size_t>(vertex)] != 0) {
                    std::cout
                        << "  vertex=" << vertex
                        << " N4="
                        << fourth[static_cast<std::size_t>(vertex)]
                        << " N5="
                        << fifth[static_cast<std::size_t>(vertex)]
                        << " margin=" << margin << '\n';
                }
            }
            return EXIT_SUCCESS;
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_d2_kernel MAXIMUM_LEVEL MAXIMUM_POWER "
                "| --case LEVEL LABEL"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_power =
            parse_positive(argv[2], "maximum power");
        if (maximum_power < 4) {
            throw std::runtime_error(
                "maximum power must be at least four"
            );
        }
        std::uint64_t parameters = 0U;
        std::uint64_t coordinates = 0U;
        std::uint64_t formula_checks = 0U;
        std::vector<PowerResult> results(
            static_cast<std::size_t>(maximum_power + 1)
        );
        std::map<std::pair<int, int>, Integer> previous_five_margin;
        std::uint64_t monotonicity_comparisons = 0U;
        std::uint64_t monotonicity_failures = 0U;
        bool correction_initialized = false;
        bool correction_failed = false;
        Integer minimum_correction = 0;
        int correction_level = 0;
        int correction_label = 0;
        int correction_vertex = 0;
        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameters;
                const std::vector<Integer> closed_four =
                    power_column(level, label, 0, 4);
                const std::vector<Integer> closed_five =
                    power_column(level, label, 0, 5);
                const Integer f4 = closed_four[0];
                const Integer f5 = closed_five[0];
                const int half_level = level / 2;
                const int half_label = label / 2;
                for (int target = 0; target <= level; target += 2) {
                    const int half_target = target / 2;
                    const Integer formula_four = affine_multiplicity(
                        4, half_level, half_label, half_target
                    );
                    const Integer formula_five = affine_multiplicity(
                        5, half_level, half_label, half_target
                    );
                    formula_checks += 2U;
                    if (formula_four
                            != closed_four[
                                static_cast<std::size_t>(target)
                            ]
                        || formula_five
                            != closed_five[
                                static_cast<std::size_t>(target)
                            ]) {
                        throw std::runtime_error(
                            "Kac-Walton formula audit mismatch"
                        );
                    }
                }
                std::vector<Integer> previous =
                    power_column(level, label, level, 3);
                for (int power = 4; power <= maximum_power; ++power) {
                    std::vector<Integer> current = previous;
                    multiply(level, label, current);
                    PowerResult& result =
                        results[static_cast<std::size_t>(power)];
                    for (int vertex = 0; vertex <= level; vertex += 2) {
                        ++coordinates;
                        const Integer margin =
                            f4 * current[
                                static_cast<std::size_t>(vertex)
                            ]
                            - f5 * previous[
                                static_cast<std::size_t>(vertex)
                            ];
                        if (!result.initialized
                            || margin < result.minimum) {
                            result.initialized = true;
                            result.minimum = margin;
                            result.level = level;
                            result.label = label;
                            result.vertex = vertex;
                        }
                        if (margin < 0) {
                            result.failed = true;
                        }
                        if (power == 5) {
                            const int reflected_vertex = level - vertex;
                            const int half_target = reflected_vertex / 2;
                            const Integer ordinary_four =
                                ordinary_multiplicity(
                                    4, half_label, half_target
                                );
                            const Integer ordinary_five =
                                ordinary_multiplicity(
                                    5, half_label, half_target
                                );
                            const Integer ordinary_f4 =
                                ordinary_multiplicity(4, half_label, 0);
                            const Integer ordinary_f5 =
                                ordinary_multiplicity(5, half_label, 0);
                            const Integer stable_margin =
                                ordinary_f4 * ordinary_five
                                - ordinary_f5 * ordinary_four;
                            const Integer correction =
                                margin - stable_margin;
                            if (!correction_initialized
                                || correction < minimum_correction) {
                                correction_initialized = true;
                                minimum_correction = correction;
                                correction_level = level;
                                correction_label = label;
                                correction_vertex = reflected_vertex;
                            }
                            if (correction < 0) {
                                correction_failed = true;
                            }
                            const std::pair<int, int> key{
                                label, reflected_vertex
                            };
                            const auto previous_position =
                                previous_five_margin.find(key);
                            if (previous_position
                                != previous_five_margin.end()) {
                                ++monotonicity_comparisons;
                                if (margin < previous_position->second) {
                                    ++monotonicity_failures;
                                    if (monotonicity_failures == 1U) {
                                        std::cout
                                            << "SU2_D2_KERNEL"
                                            << " MONOTONICITY_COUNTEREXAMPLE"
                                            << " level=" << level
                                            << " label=" << label
                                            << " reflected_vertex="
                                            << reflected_vertex
                                            << " prior="
                                            << previous_position->second
                                            << " current=" << margin
                                            << '\n';
                                    }
                                }
                            }
                            previous_five_margin[key] = margin;
                        }
                    }
                    previous.swap(current);
                }
            }
        }
        std::cout
            << "SU2_D2_KERNEL"
            << " maximum_level=" << maximum_level
            << " maximum_power=" << maximum_power
            << " parameters=" << parameters
            << " coordinates=" << coordinates
            << " formula_checks=" << formula_checks
            << " monotonicity_comparisons="
                << monotonicity_comparisons
            << " monotonicity_failures=" << monotonicity_failures
            << " result=PASS_DISCOVERY\n";
        std::cout
            << "FIVE_STEP_WALL_CORRECTION"
            << " minimum=" << minimum_correction
            << " witness=(" << correction_level
            << ',' << correction_label
            << ',' << correction_vertex << ')'
            << " result="
            << (correction_failed
                ? "FAIL_POINTWISE"
                : "PASS_POINTWISE")
            << '\n';
        for (int power = 4; power <= maximum_power; ++power) {
            const PowerResult& result =
                results[static_cast<std::size_t>(power)];
            std::cout
                << "POWER power=" << power
                << " minimum=" << result.minimum
                << " witness=(" << result.level
                << ',' << result.label
                << ',' << result.vertex << ')'
                << " result="
                << (result.failed ? "FAIL_POINTWISE" : "PASS_POINTWISE")
                << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_D2_KERNEL FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
