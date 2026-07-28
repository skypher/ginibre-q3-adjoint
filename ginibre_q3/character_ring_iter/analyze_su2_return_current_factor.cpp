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

struct Minimum {
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
                "usage: analyze_su2_return_current_factor "
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

        Minimum minimum;
        Minimum last_exit_minimum;
        std::uint64_t parameters = 0U;
        std::uint64_t rows = 0U;
        std::uint64_t coefficients = 0U;
        std::uint64_t negatives = 0U;
        std::uint64_t last_exit_negatives = 0U;

        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameters;
                const int maximum_degree = 2 * maximum_prefix + 2;
                std::vector<Integer> returns(
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
                }

                std::vector<Integer> return_square(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                for (int degree = 0;
                     degree <= maximum_degree;
                     ++degree) {
                    for (int left = 0; left <= degree; ++left) {
                        return_square[
                            static_cast<std::size_t>(degree)
                        ] += returns[static_cast<std::size_t>(left)]
                            * returns[
                                static_cast<std::size_t>(degree - left)
                            ];
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
                        for (int degree = 0; degree <= n; ++degree) {
                            Integer value = 0;
                            Integer last_exit_value = 0;
                            for (int index = 0;
                                 index <= degree;
                                 ++index) {
                                value += current_polynomial[
                                    static_cast<std::size_t>(index)
                                ] * return_square[
                                    static_cast<std::size_t>(
                                        degree - index
                                    )
                                ];
                                last_exit_value += current_polynomial[
                                    static_cast<std::size_t>(index)
                                ] * returns[
                                    static_cast<std::size_t>(
                                        degree - index
                                    )
                                ];
                            }
                            ++coefficients;
                            minimum.observe(
                                value,
                                level,
                                label,
                                prefix,
                                truncation,
                                degree
                            );
                            if (value < 0) {
                                ++negatives;
                                if (negatives == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " degree=" << degree
                                        << " value=" << value << '\n';
                                }
                            }
                            last_exit_minimum.observe(
                                last_exit_value,
                                level,
                                label,
                                prefix,
                                truncation,
                                degree
                            );
                            if (last_exit_value < 0) {
                                ++last_exit_negatives;
                                if (last_exit_negatives == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_LAST_EXIT"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " degree=" << degree
                                        << " value=" << last_exit_value
                                        << '\n';
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_RETURN_CURRENT_FACTOR"
            << " maximum_level=" << maximum_level
            << " maximum_prefix=" << maximum_prefix
            << " parameters=" << parameters
            << " rows=" << rows
            << " coefficients=" << coefficients
            << " negatives=" << negatives
            << " last_exit_negatives=" << last_exit_negatives
            << " minimum=" << minimum.value
            << " witness=("
            << minimum.level << ','
            << minimum.label << ','
            << minimum.prefix << ','
            << minimum.truncation << ','
            << minimum.degree << ')'
            << " last_exit_minimum=" << last_exit_minimum.value
            << " last_exit_witness=("
            << last_exit_minimum.level << ','
            << last_exit_minimum.label << ','
            << last_exit_minimum.prefix << ','
            << last_exit_minimum.truncation << ','
            << last_exit_minimum.degree << ')'
            << " result="
            << (negatives == 0U
                ? "PASS_COEFFICIENTWISE_DISCOVERY"
                : "FAIL_COEFFICIENTWISE")
            << '\n';
        return negatives == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_RETURN_CURRENT_FACTOR FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
