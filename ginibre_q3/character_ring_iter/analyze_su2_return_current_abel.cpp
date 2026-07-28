#include <algorithm>
#include <cstdint>
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

std::vector<Integer> multiply(
    int level,
    int label,
    const std::vector<Integer>& state
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
    return next;
}

Integer binomial(int n, int r) {
    if (r < 0 || r > n) {
        return 0;
    }
    r = std::min(r, n - r);
    Integer result = 1;
    for (int index = 1; index <= r; ++index) {
        result *= n - r + index;
        result /= index;
    }
    return result;
}

struct Witness {
    bool initialized = false;
    Integer value = 0;
    int level = 0;
    int label = 0;
    int prefix = 0;
    int truncation = 0;
    int degree = 0;

    void observe(
        const Integer& candidate,
        int candidate_level,
        int candidate_label,
        int candidate_prefix,
        int candidate_truncation,
        int candidate_degree
    ) {
        if (!initialized || candidate < value) {
            initialized = true;
            value = candidate;
            level = candidate_level;
            label = candidate_label;
            prefix = candidate_prefix;
            truncation = candidate_truncation;
            degree = candidate_degree;
        }
    }
};

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_return_current_abel "
                "MAXIMUM_LEVEL MAXIMUM_PREFIX"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_prefix =
            parse_positive(argv[2], "maximum prefix");
        if (maximum_prefix < 4) {
            throw std::runtime_error(
                "maximum prefix must be at least four"
            );
        }

        std::uint64_t parameters = 0U;
        std::uint64_t rows = 0U;
        std::uint64_t cumulative_coefficients = 0U;
        std::uint64_t cumulative_negatives = 0U;
        std::uint64_t second_cumulative_negatives = 0U;
        std::uint64_t first_passage_steps = 0U;
        std::uint64_t first_passage_decreases = 0U;
        Witness cumulative_minimum;
        Witness second_cumulative_minimum;
        Witness first_passage_minimum;

        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameters;
                const int maximum_degree = 2 * maximum_prefix + 2;
                std::vector<Integer> returns(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                std::vector<Integer> endpoints(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                std::vector<Integer> state(
                    static_cast<std::size_t>(level + 1)
                );
                state[0] = 1;
                returns[0] = 1;
                for (int power = 1;
                     power <= maximum_degree;
                     ++power) {
                    state = multiply(level, label, state);
                    returns[static_cast<std::size_t>(power)] = state[0];
                    endpoints[static_cast<std::size_t>(power)] =
                        state[static_cast<std::size_t>(level)];
                }

                // G=FE, so division by F (whose constant term is one)
                // gives the exact first-passage series E.
                std::vector<Integer> first_passage(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                for (int degree = 0;
                     degree <= maximum_degree;
                     ++degree) {
                    Integer value =
                        endpoints[static_cast<std::size_t>(degree)];
                    for (int index = 1; index <= degree; ++index) {
                        value -= returns[static_cast<std::size_t>(index)]
                            * first_passage[
                                static_cast<std::size_t>(degree - index)
                            ];
                    }
                    first_passage[static_cast<std::size_t>(degree)] =
                        value;
                    if (degree == 0) {
                        continue;
                    }
                    ++first_passage_steps;
                    const Integer difference =
                        first_passage[static_cast<std::size_t>(degree)]
                        - first_passage[
                            static_cast<std::size_t>(degree - 1)
                        ];
                    first_passage_minimum.observe(
                        difference,
                        level,
                        label,
                        0,
                        0,
                        degree
                    );
                    if (difference < 0) {
                        ++first_passage_decreases;
                    }
                }

                for (int prefix = 4;
                     prefix <= maximum_prefix;
                     ++prefix) {
                    const int n = 2 * prefix + 2;
                    std::vector<Integer> current_polynomial(
                        static_cast<std::size_t>(n + 1)
                    );
                    for (int truncation = 0;
                         truncation <= prefix - 1;
                         ++truncation) {
                        const int even = 2 * truncation;
                        const int odd = even + 1;
                        const Integer weight = binomial(n - 1, even);
                        current_polynomial[
                            static_cast<std::size_t>(even)
                        ] += weight
                            * returns[static_cast<std::size_t>(even)];
                        current_polynomial[
                            static_cast<std::size_t>(odd)
                        ] -= weight
                            * returns[static_cast<std::size_t>(odd)];
                        if (truncation < 2) {
                            continue;
                        }
                        ++rows;
                        Integer cumulative = 0;
                        Integer second_cumulative = 0;
                        for (int degree = 0; degree <= n; ++degree) {
                            Integer coefficient = 0;
                            for (int index = 0;
                                 index <= degree;
                                 ++index) {
                                coefficient += current_polynomial[
                                    static_cast<std::size_t>(index)
                                ] * returns[
                                    static_cast<std::size_t>(
                                        degree - index
                                    )
                                ];
                            }
                            cumulative += coefficient;
                            second_cumulative += cumulative;
                            ++cumulative_coefficients;
                            cumulative_minimum.observe(
                                cumulative,
                                level,
                                label,
                                prefix,
                                truncation,
                                degree
                            );
                            if (cumulative < 0) {
                                ++cumulative_negatives;
                                if (cumulative_negatives == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_CUMULATIVE"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " degree=" << degree
                                        << " value=" << cumulative
                                        << '\n';
                                }
                            }
                            second_cumulative_minimum.observe(
                                second_cumulative,
                                level,
                                label,
                                prefix,
                                truncation,
                                degree
                            );
                            if (second_cumulative < 0) {
                                ++second_cumulative_negatives;
                                if (second_cumulative_negatives == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_SECOND_CUMULATIVE"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " degree=" << degree
                                        << " value=" << second_cumulative
                                        << '\n';
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_RETURN_CURRENT_ABEL"
            << " maximum_level=" << maximum_level
            << " maximum_prefix=" << maximum_prefix
            << " parameters=" << parameters
            << " rows=" << rows
            << " cumulative_coefficients=" << cumulative_coefficients
            << " cumulative_negatives=" << cumulative_negatives
            << " cumulative_minimum=" << cumulative_minimum.value
            << " cumulative_witness=("
            << cumulative_minimum.level << ','
            << cumulative_minimum.label << ','
            << cumulative_minimum.prefix << ','
            << cumulative_minimum.truncation << ','
            << cumulative_minimum.degree << ')'
            << " second_cumulative_negatives="
                << second_cumulative_negatives
            << " second_cumulative_minimum="
                << second_cumulative_minimum.value
            << " second_cumulative_witness=("
            << second_cumulative_minimum.level << ','
            << second_cumulative_minimum.label << ','
            << second_cumulative_minimum.prefix << ','
            << second_cumulative_minimum.truncation << ','
            << second_cumulative_minimum.degree << ')'
            << " first_passage_steps=" << first_passage_steps
            << " first_passage_decreases=" << first_passage_decreases
            << " first_passage_minimum="
                << first_passage_minimum.value
            << " first_passage_witness=("
            << first_passage_minimum.level << ','
            << first_passage_minimum.label << ','
            << first_passage_minimum.degree << ')'
            << " result="
            << (
                second_cumulative_negatives == 0U
                    && first_passage_decreases == 0U
                ? "PASS_SECOND_ABEL_DISCOVERY"
                : "FAIL_ABEL_CANDIDATE"
            )
            << '\n';
        return second_cumulative_negatives == 0U
                && first_passage_decreases == 0U
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_RETURN_CURRENT_ABEL FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
