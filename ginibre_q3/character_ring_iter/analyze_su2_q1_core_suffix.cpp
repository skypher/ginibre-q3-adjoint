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

std::vector<Integer> multiply_q1(
    const std::vector<Integer>& state
) {
    const int half_level = static_cast<int>(state.size()) - 1;
    std::vector<Integer> next(state.size());
    next[0] = state[1];
    next[static_cast<std::size_t>(half_level)] =
        state[static_cast<std::size_t>(half_level - 1)];
    for (int vertex = 1; vertex < half_level; ++vertex) {
        next[static_cast<std::size_t>(vertex)] =
            state[static_cast<std::size_t>(vertex - 1)]
            + state[static_cast<std::size_t>(vertex)]
            + state[static_cast<std::size_t>(vertex + 1)];
    }
    return next;
}

std::vector<Integer> divide_by_fibonacci(
    const std::vector<Integer>& series
) {
    std::vector<Integer> quotient(series.size());
    for (std::size_t degree = 0U;
         degree < series.size();
         ++degree) {
        quotient[degree] = series[degree];
        if (degree >= 1U) {
            quotient[degree] += quotient[degree - 1U];
        }
        if (degree >= 2U) {
            quotient[degree] += quotient[degree - 2U];
        }
    }
    return quotient;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_q1_core_suffix "
                "MAXIMUM_HALF_LEVEL MAXIMUM_PREFIX"
            );
        }
        const int maximum_half_level =
            parse_positive(argv[1], "maximum half-level");
        const int maximum_prefix =
            parse_positive(argv[2], "maximum prefix");
        if (maximum_half_level < 3 || maximum_prefix < 4) {
            throw std::runtime_error(
                "half-level must be at least three and prefix at least four"
            );
        }

        std::uint64_t parameters = 0U;
        std::uint64_t rows = 0U;
        std::uint64_t identities = 0U;
        std::uint64_t residual_coordinates = 0U;
        std::uint64_t residual_negatives = 0U;
        std::uint64_t suffix_coordinates = 0U;
        std::uint64_t suffix_negatives = 0U;
        std::uint64_t hard_suffix_coordinates = 0U;
        std::uint64_t hard_suffix_negatives = 0U;
        std::uint64_t negative_targets = 0U;
        bool printed_first_suffix_negative = false;
        bool printed_first_hard_suffix_negative = false;

        const int maximum_degree = 2 * maximum_prefix + 2;
        for (int half_level = 3;
             half_level <= maximum_half_level;
             ++half_level) {
            ++parameters;
            std::vector<Integer> returns(
                static_cast<std::size_t>(maximum_degree + 1)
            );
            std::vector<Integer> endpoints(
                static_cast<std::size_t>(maximum_degree + 1)
            );
            std::vector<Integer> state(
                static_cast<std::size_t>(half_level + 1)
            );
            state[0] = 1;
            returns[0] = 1;
            for (int degree = 1;
                 degree <= maximum_degree;
                 ++degree) {
                state = multiply_q1(state);
                returns[static_cast<std::size_t>(degree)] = state[0];
                endpoints[static_cast<std::size_t>(degree)] =
                    state[static_cast<std::size_t>(half_level)];
            }

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
                first_passage[static_cast<std::size_t>(degree)] = value;
            }

            std::vector<Integer> denominator(
                static_cast<std::size_t>(maximum_degree + 1)
            );
            denominator[0] = 1;
            for (int copy = 0; copy < half_level - 1; ++copy) {
                std::vector<Integer> product(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                for (int degree = 0;
                     degree <= maximum_degree;
                     ++degree) {
                    product[static_cast<std::size_t>(degree)] =
                        denominator[static_cast<std::size_t>(degree)];
                    if (degree >= 1) {
                        product[static_cast<std::size_t>(degree)] -=
                            denominator[
                                static_cast<std::size_t>(degree - 1)
                            ];
                    }
                    if (degree >= 2) {
                        product[static_cast<std::size_t>(degree)] -=
                            denominator[
                                static_cast<std::size_t>(degree - 2)
                            ];
                    }
                }
                denominator = std::move(product);
            }

            std::vector<Integer> residual(
                static_cast<std::size_t>(maximum_degree + 1)
            );
            for (int degree = 0;
                 degree <= maximum_degree;
                 ++degree) {
                for (int left = 0; left <= degree; ++left) {
                    residual[static_cast<std::size_t>(degree)] +=
                        first_passage[static_cast<std::size_t>(left)]
                        * denominator[
                            static_cast<std::size_t>(degree - left)
                        ];
                }
                ++residual_coordinates;
                if (residual[static_cast<std::size_t>(degree)] < 0) {
                    ++residual_negatives;
                }
            }
            if (
                residual[static_cast<std::size_t>(half_level)] <= 0
            ) {
                throw std::runtime_error(
                    "incorrect minimal residual coefficient"
                );
            }

            for (int prefix = 4;
                 prefix <= maximum_prefix;
                 ++prefix) {
                const int n = 2 * prefix + 2;
                std::vector<Integer> current(
                    static_cast<std::size_t>(n + 1)
                );
                Integer direct_target = 0;
                for (int truncation = 0;
                     truncation < prefix;
                     ++truncation) {
                    const int even = 2 * truncation;
                    const int odd = even + 1;
                    const Integer weight = binomial(n - 1, even);
                    current[static_cast<std::size_t>(even)] +=
                        weight * returns[static_cast<std::size_t>(even)];
                    current[static_cast<std::size_t>(odd)] -=
                        weight * returns[static_cast<std::size_t>(odd)];
                    direct_target += weight * (
                        returns[static_cast<std::size_t>(even)]
                            * endpoints[
                                static_cast<std::size_t>(n - even)
                            ]
                        - returns[static_cast<std::size_t>(odd)]
                            * endpoints[
                                static_cast<std::size_t>(n - odd)
                            ]
                    );

                    std::vector<Integer> h(
                        static_cast<std::size_t>(n + 1)
                    );
                    for (int degree = 0; degree <= n; ++degree) {
                        for (int index = 0; index <= degree; ++index) {
                            h[static_cast<std::size_t>(degree)] +=
                                current[static_cast<std::size_t>(index)]
                                * returns[
                                    static_cast<std::size_t>(
                                        degree - index
                                    )
                                ];
                        }
                    }
                    for (int copy = 0;
                         copy < half_level - 1;
                         ++copy) {
                        h = divide_by_fibonacci(h);
                    }

                    Integer suffix = 0;
                    for (int residual_degree = n;
                         residual_degree >= half_level;
                         --residual_degree) {
                        suffix += residual[
                            static_cast<std::size_t>(residual_degree)
                        ] * h[
                            static_cast<std::size_t>(
                                n - residual_degree
                            )
                        ];
                        ++suffix_coordinates;
                        if (truncation >= 4) {
                            ++hard_suffix_coordinates;
                        }
                        if (suffix < 0) {
                            ++suffix_negatives;
                            if (!printed_first_suffix_negative) {
                                printed_first_suffix_negative = true;
                                std::cout
                                    << "FIRST_NEGATIVE_CORE_SUFFIX"
                                    << " half_level=" << half_level
                                    << " prefix=" << prefix
                                    << " truncation=" << truncation
                                    << " residual_degree="
                                        << residual_degree
                                    << " value=" << suffix << '\n';
                            }
                            if (truncation >= 4) {
                                ++hard_suffix_negatives;
                                if (
                                    !printed_first_hard_suffix_negative
                                ) {
                                    printed_first_hard_suffix_negative =
                                        true;
                                    std::cout
                                        << "FIRST_HARD_NEGATIVE_CORE_SUFFIX"
                                        << " half_level=" << half_level
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " residual_degree="
                                            << residual_degree
                                        << " value=" << suffix << '\n';
                                }
                            }
                        }
                    }
                    ++rows;
                    ++identities;
                    if (suffix != direct_target) {
                        throw std::runtime_error(
                            "core convolution disagrees with endpoint current"
                        );
                    }
                    if (direct_target < 0) {
                        ++negative_targets;
                    }
                }
            }
        }

        std::cout
            << "SU2_Q1_CORE_SUFFIX"
            << " maximum_half_level=" << maximum_half_level
            << " maximum_prefix=" << maximum_prefix
            << " parameters=" << parameters
            << " rows=" << rows
            << " identities=" << identities
            << " residual_coordinates=" << residual_coordinates
            << " residual_negatives=" << residual_negatives
            << " suffix_coordinates=" << suffix_coordinates
            << " suffix_negatives=" << suffix_negatives
            << " hard_suffix_coordinates=" << hard_suffix_coordinates
            << " hard_suffix_negatives=" << hard_suffix_negatives
            << " negative_targets=" << negative_targets
            << " result="
                << (
                    residual_negatives == 0U
                        && hard_suffix_negatives == 0U
                        && negative_targets == 0U
                        ? "PASS_Q1_CORE_SUFFIX_DISCOVERY"
                        : "FAIL_Q1_CORE_SUFFIX"
                )
            << '\n';
        return residual_negatives == 0U
                && hard_suffix_negatives == 0U
                && negative_targets == 0U
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
