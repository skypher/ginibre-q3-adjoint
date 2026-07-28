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

void multiply_minus(
    int level,
    int label,
    std::vector<Integer>& state
) {
    const int width = level + 1;
    std::vector<Integer> next(
        static_cast<std::size_t>(width * width)
    );
    const auto index = [width](int first, int second) {
        return static_cast<std::size_t>(first * width + second);
    };
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
                    next[index(first, output)] -= coefficient;
                }
            }
        }
    }
    state.swap(next);
}

void multiply_plus(
    int level,
    int label,
    std::vector<Integer>& state
) {
    const int width = level + 1;
    std::vector<Integer> next(
        static_cast<std::size_t>(width * width)
    );
    const auto index = [width](int first, int second) {
        return static_cast<std::size_t>(first * width + second);
    };
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

Integer signed_value(
    int level,
    int label,
    int first_count,
    int second_count
) {
    const int width = level + 1;
    std::vector<Integer> state(
        static_cast<std::size_t>(width * width)
    );
    state[0] = 1;
    for (int index = 0; index < first_count; ++index) {
        multiply_minus(level, label, state);
    }
    for (int index = 0; index < second_count; ++index) {
        multiply_minus(level, level - label, state);
    }
    return state[0];
}

struct GroupedBoundaryValues {
    Integer all_minus;
    Integer mixed;
};

GroupedBoundaryValues grouped_boundary_values(
    int level,
    int label,
    int first_count,
    int second_count
) {
    const int width = level + 1;
    const auto boundary = [width, level](
        const std::vector<Integer>& state
    ) -> Integer {
        return state[static_cast<std::size_t>(level * width)];
    };

    std::vector<Integer> all_minus(
        static_cast<std::size_t>(width * width)
    );
    all_minus[0] = 1;
    for (int index = 0;
         index < first_count + second_count; ++index) {
        multiply_minus(level, label, all_minus);
    }

    std::vector<Integer> mixed(
        static_cast<std::size_t>(width * width)
    );
    mixed[0] = 1;
    for (int index = 0; index < first_count; ++index) {
        multiply_minus(level, label, mixed);
    }
    for (int index = 0; index < second_count; ++index) {
        multiply_plus(level, label, mixed);
    }
    return {
        boundary(all_minus),
        boundary(mixed)
    };
}

std::vector<std::vector<Integer>> binomial_table(int maximum) {
    std::vector<std::vector<Integer>> choose(
        static_cast<std::size_t>(maximum + 1)
    );
    for (int n = 0; n <= maximum; ++n) {
        choose[static_cast<std::size_t>(n)].assign(
            static_cast<std::size_t>(n + 1), Integer{0}
        );
        choose[static_cast<std::size_t>(n)][0] = 1;
        choose[static_cast<std::size_t>(n)][
            static_cast<std::size_t>(n)
        ] = 1;
        for (int r = 1; r < n; ++r) {
            choose[static_cast<std::size_t>(n)][
                static_cast<std::size_t>(r)
            ] =
                choose[static_cast<std::size_t>(n - 1)][
                    static_cast<std::size_t>(r - 1)
                ]
                + choose[static_cast<std::size_t>(n - 1)][
                    static_cast<std::size_t>(r)
                ];
        }
    }
    return choose;
}

Integer alternating_path_value(
    int level,
    int label,
    int first_count,
    int second_count
) {
    const int total = first_count + second_count;
    const auto choose = binomial_table(total);
    std::vector<Integer> state(
        static_cast<std::size_t>(level + 1)
    );
    std::vector<Integer> next(
        static_cast<std::size_t>(level + 1)
    );
    std::vector<Integer> closed(
        static_cast<std::size_t>(total + 1)
    );
    std::vector<Integer> wall(
        static_cast<std::size_t>(total + 1)
    );
    state[0] = 1;
    closed[0] = 1;
    for (int power = 1; power <= total; ++power) {
        std::fill(next.begin(), next.end(), Integer{0});
        for (int source = 0; source <= level; ++source) {
            if (state[static_cast<std::size_t>(source)] == 0) {
                continue;
            }
            for (int output = 0; output <= level; ++output) {
                if (fuses(level, label, source, output)) {
                    next[static_cast<std::size_t>(output)]
                        += state[static_cast<std::size_t>(source)];
                }
            }
        }
        state.swap(next);
        closed[static_cast<std::size_t>(power)] = state[0];
        wall[static_cast<std::size_t>(power)] =
            state[static_cast<std::size_t>(level)];
    }

    Integer result = 0;
    for (int exponent = 0; exponent <= total; ++exponent) {
        Integer weight = 0;
        for (int odd = 1; odd <= second_count; odd += 2) {
            const int from_first = exponent - odd;
            if (from_first < 0 || from_first > first_count) {
                continue;
            }
            weight +=
                choose[static_cast<std::size_t>(second_count)][
                    static_cast<std::size_t>(odd)
                ]
                * choose[static_cast<std::size_t>(first_count)][
                    static_cast<std::size_t>(from_first)
                ];
        }
        Integer term =
            2 * weight * wall[static_cast<std::size_t>(exponent)]
            * closed[static_cast<std::size_t>(total - exponent)];
        if ((exponent & 1) == 0) {
            result += term;
        } else {
            result -= term;
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 5) {
            const int level = parse_positive(argv[1], "level");
            const int label = parse_positive(argv[2], "label");
            const int first = parse_positive(argv[3], "first count");
            const int second = parse_positive(argv[4], "second count");
            if (level < 2 || label >= level - label
                || (first & 1) == 0 || (second & 1) == 0) {
                throw std::runtime_error(
                    "require label<level-label and odd counts"
                );
            }
            const Integer value =
                signed_value(level, label, first, second);
            const GroupedBoundaryValues grouped =
                grouped_boundary_values(
                    level, label, first, second
                );
            const Integer paths = alternating_path_value(
                level, label, first, second
            );
            std::cout
                << "SU2_SIMPLE_CURRENT_PAIR_SINGLE"
                << " level=" << level
                << " labels={" << label << ',' << level - label << '}'
                << " counts={" << first << ',' << second << '}'
                << " value=" << value
                << " all_minus_boundary=" << grouped.all_minus
                << " mixed_boundary=" << grouped.mixed
                << " grouped_boundary="
                << grouped.all_minus + grouped.mixed
                << " alternating_paths=" << paths
                << " identity="
                << (value == grouped.all_minus + grouped.mixed
                        && value == paths
                        ? "PASS" : "FAIL")
                << '\n';
            return value == grouped.all_minus + grouped.mixed
                    && value == paths
                ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_simple_current_pair "
                "MAXIMUM_LEVEL MAXIMUM_ODD_COUNT"
                " | LEVEL LABEL ODD_COUNT_1 ODD_COUNT_2"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_count =
            parse_positive(argv[2], "maximum odd count");
        if (maximum_level < 2 || maximum_count > 99) {
            throw std::runtime_error(
                "require level>=2 and odd count bound<=99"
            );
        }

        bool initialized = false;
        Integer minimum = 0;
        int minimum_level = 0;
        int minimum_label = 0;
        int minimum_first = 0;
        int minimum_second = 0;
        Integer maximum = 0;
        int maximum_level_witness = 0;
        int maximum_label_witness = 0;
        int maximum_first = 0;
        int maximum_second = 0;
        bool nonbase_initialized = false;
        Integer nonbase_minimum = 0;
        int nonbase_level = 0;
        int nonbase_label = 0;
        int nonbase_first = 0;
        int nonbase_second = 0;
        bool all_minus_initialized = false;
        Integer all_minus_minimum = 0;
        int all_minus_level = 0;
        int all_minus_label = 0;
        int all_minus_first = 0;
        int all_minus_second = 0;
        bool mixed_initialized = false;
        Integer mixed_minimum = 0;
        int mixed_level = 0;
        int mixed_label = 0;
        int mixed_first = 0;
        int mixed_second = 0;
        std::size_t zeros = 0U;
        std::size_t rows = 0U;
        for (int level = 2; level <= maximum_level; ++level) {
            for (int label = 1; label < level - label; ++label) {
                for (int first = 1; first <= maximum_count; first += 2) {
                    for (int second = 1;
                         second <= maximum_count; second += 2) {
                        const Integer value = signed_value(
                            level, label, first, second
                        );
                        const GroupedBoundaryValues grouped =
                            grouped_boundary_values(
                                level, label, first, second
                            );
                        const Integer paths = alternating_path_value(
                            level, label, first, second
                        );
                        if (value
                                != grouped.all_minus + grouped.mixed
                            || value != paths) {
                            std::cout
                                << "SU2_SIMPLE_CURRENT_PAIR IDENTITY_FAIL"
                                << " level=" << level
                                << " labels={" << label << ','
                                << level - label << '}'
                                << " counts={" << first << ','
                                << second << '}'
                                << " value=" << value
                                << " all_minus_boundary="
                                << grouped.all_minus
                                << " mixed_boundary=" << grouped.mixed
                                << " grouped_boundary="
                                << grouped.all_minus + grouped.mixed
                                << " alternating_paths=" << paths << '\n';
                            return EXIT_FAILURE;
                        }
                        if (
                            grouped.all_minus < 0
                            || grouped.mixed < 0
                        ) {
                            std::cout
                                << "SU2_SIMPLE_CURRENT_PAIR"
                                << " COMPONENT_COUNTEREXAMPLE"
                                << " level=" << level
                                << " labels={" << label << ','
                                << level - label << '}'
                                << " counts={" << first << ','
                                << second << '}'
                                << " all_minus_boundary="
                                << grouped.all_minus
                                << " mixed_boundary=" << grouped.mixed
                                << '\n';
                            return EXIT_FAILURE;
                        }
                        if (!all_minus_initialized
                            || grouped.all_minus < all_minus_minimum) {
                            all_minus_initialized = true;
                            all_minus_minimum = grouped.all_minus;
                            all_minus_level = level;
                            all_minus_label = label;
                            all_minus_first = first;
                            all_minus_second = second;
                        }
                        if (!mixed_initialized
                            || grouped.mixed < mixed_minimum) {
                            mixed_initialized = true;
                            mixed_minimum = grouped.mixed;
                            mixed_level = level;
                            mixed_label = label;
                            mixed_first = first;
                            mixed_second = second;
                        }
                        ++rows;
                        if (!initialized || value < minimum) {
                            initialized = true;
                            minimum = value;
                            minimum_level = level;
                            minimum_label = label;
                            minimum_first = first;
                            minimum_second = second;
                        }
                        if (!initialized || value > maximum) {
                            maximum = value;
                            maximum_level_witness = level;
                            maximum_label_witness = label;
                            maximum_first = first;
                            maximum_second = second;
                        }
                        if (value == 0) {
                            ++zeros;
                        }
                        if (first + second > 2
                            && (!nonbase_initialized
                                || value < nonbase_minimum)) {
                            nonbase_initialized = true;
                            nonbase_minimum = value;
                            nonbase_level = level;
                            nonbase_label = label;
                            nonbase_first = first;
                            nonbase_second = second;
                        }
                        if (value < 0) {
                            std::cout
                                << "SU2_SIMPLE_CURRENT_PAIR COUNTEREXAMPLE"
                                << " level=" << level
                                << " labels={" << label << ','
                                << level - label << '}'
                                << " counts={" << first << ','
                                << second << '}'
                                << " value=" << value << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }
        std::cout
            << "SU2_SIMPLE_CURRENT_PAIR"
            << " rows=" << rows
            << " maximum_level=" << maximum_level
            << " maximum_odd_count=" << maximum_count
            << " minimum=" << minimum
            << " maximum=" << maximum
            << " zeros=" << zeros
            << " witness_level=" << minimum_level
            << " witness_labels={" << minimum_label << ','
            << minimum_level - minimum_label << '}'
            << " witness_counts={" << minimum_first << ','
            << minimum_second << '}'
            << " maximum_witness_level=" << maximum_level_witness
            << " maximum_witness_labels={"
            << maximum_label_witness << ','
            << maximum_level_witness - maximum_label_witness << '}'
            << " maximum_witness_counts={" << maximum_first << ','
            << maximum_second << '}'
            << " nonbase_minimum=" << nonbase_minimum
            << " nonbase_witness_level=" << nonbase_level
            << " nonbase_witness_labels={" << nonbase_label << ','
            << nonbase_level - nonbase_label << '}'
            << " nonbase_witness_counts={" << nonbase_first << ','
            << nonbase_second << '}'
            << " all_minus_minimum=" << all_minus_minimum
            << " all_minus_witness_level=" << all_minus_level
            << " all_minus_witness_labels={"
            << all_minus_label << ','
            << all_minus_level - all_minus_label << '}'
            << " all_minus_witness_counts={" << all_minus_first << ','
            << all_minus_second << '}'
            << " mixed_minimum=" << mixed_minimum
            << " mixed_witness_level=" << mixed_level
            << " mixed_witness_labels={" << mixed_label << ','
            << mixed_level - mixed_label << '}'
            << " mixed_witness_counts={" << mixed_first << ','
            << mixed_second << '}'
            << " grouped_identity=PASS alternating_path_identity=PASS"
            << " counterexamples=0 result=PASS_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SIMPLE_CURRENT_PAIR FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
