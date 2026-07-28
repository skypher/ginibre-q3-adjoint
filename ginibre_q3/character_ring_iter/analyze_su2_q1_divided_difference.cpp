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

Integer falling(int n, int r) {
    Integer result = 1;
    for (int index = 0; index < r; ++index) {
        result *= n - index;
    }
    return result;
}

Integer factorial(int n) {
    return falling(n, n);
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

Integer power_of_four(int exponent) {
    Integer result = 1;
    for (int index = 0; index < exponent; ++index) {
        result *= 4;
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_q1_divided_difference "
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

        std::uint64_t rows = 0U;
        std::uint64_t bernstein_coordinates = 0U;
        std::uint64_t negative_bernstein_coordinates = 0U;
        std::uint64_t hard_negative_bernstein_coordinates = 0U;
        std::uint64_t negative_targets = 0U;
        std::uint64_t bivariate_coordinates = 0U;
        std::uint64_t negative_bivariate_coordinates = 0U;
        bool printed_first_bernstein_negative = false;
        bool printed_first_hard_bernstein_negative = false;
        bool printed_first_target_negative = false;
        bool printed_first_bivariate_negative = false;

        for (int half_level = 3;
             half_level <= maximum_half_level;
             ++half_level) {
            const int maximum_degree = 2 * maximum_prefix + 2;
            std::vector<Integer> returns(
                static_cast<std::size_t>(maximum_degree + 1)
            );
            std::vector<Integer> endpoints(
                static_cast<std::size_t>(maximum_degree + 1)
            );
            std::vector<std::vector<Integer>> powers(
                static_cast<std::size_t>(maximum_degree + 1),
                std::vector<Integer>(
                    static_cast<std::size_t>(half_level + 1)
                )
            );
            std::vector<Integer> state(
                static_cast<std::size_t>(half_level + 1)
            );
            state[0] = 1;
            powers[0] = state;
            returns[0] = 1;
            for (int degree = 1; degree <= maximum_degree; ++degree) {
                state = multiply_q1(state);
                powers[static_cast<std::size_t>(degree)] = state;
                returns[static_cast<std::size_t>(degree)] = state[0];
                endpoints[static_cast<std::size_t>(degree)] =
                    state[static_cast<std::size_t>(half_level)];
            }

            for (int prefix = 4; prefix <= maximum_prefix; ++prefix) {
                const int m = 2 * prefix + 1;
                for (int truncation = 0;
                     truncation < prefix;
                     ++truncation) {
                    ++rows;
                    std::vector<Integer> polynomial(
                        static_cast<std::size_t>(m + 2)
                    );
                    Integer target = 0;
                    for (int slice = 0;
                         slice <= truncation;
                         ++slice) {
                        const int even = 2 * slice;
                        const int first_degree = m - even + 1;
                        const int second_degree = m - even;
                        const Integer weight = binomial(m, even);
                        polynomial[
                            static_cast<std::size_t>(first_degree)
                        ] += weight
                            * returns[static_cast<std::size_t>(even)];
                        polynomial[
                            static_cast<std::size_t>(second_degree)
                        ] -= weight
                            * returns[
                                static_cast<std::size_t>(even + 1)
                            ];
                        target += weight * (
                            endpoints[
                                static_cast<std::size_t>(first_degree)
                            ] * returns[static_cast<std::size_t>(even)]
                            - endpoints[
                                static_cast<std::size_t>(second_degree)
                            ] * returns[
                                static_cast<std::size_t>(even + 1)
                            ]
                        );
                    }
                    if (target < 0) {
                        ++negative_targets;
                        if (!printed_first_target_negative) {
                            printed_first_target_negative = true;
                            std::cout
                                << "FIRST_NEGATIVE_TARGET"
                                << " half_level=" << half_level
                                << " prefix=" << prefix
                                << " truncation=" << truncation
                                << " value=" << target << '\n';
                        }
                    }

                    if (truncation >= 4) {
                        for (int first_target = 0;
                             first_target <= half_level;
                             ++first_target) {
                            for (int second_target = 0;
                                 second_target <= half_level;
                                 ++second_target) {
                                Integer bivariate = 0;
                                for (int slice = 0;
                                     slice <= truncation;
                                     ++slice) {
                                    const int even = 2 * slice;
                                    const int first_degree = m - even;
                                    const Integer weight =
                                        binomial(m, even);
                                    bivariate += weight * (
                                        powers[
                                            static_cast<std::size_t>(
                                                first_degree + 1
                                            )
                                        ][static_cast<std::size_t>(
                                            first_target
                                        )] * powers[
                                            static_cast<std::size_t>(even)
                                        ][static_cast<std::size_t>(
                                            second_target
                                        )]
                                        - powers[
                                            static_cast<std::size_t>(
                                                first_degree
                                            )
                                        ][static_cast<std::size_t>(
                                            first_target
                                        )] * powers[
                                            static_cast<std::size_t>(
                                                even + 1
                                            )
                                        ][static_cast<std::size_t>(
                                            second_target
                                        )]
                                    );
                                }
                                ++bivariate_coordinates;
                                if (bivariate < 0) {
                                    ++negative_bivariate_coordinates;
                                    if (
                                        !printed_first_bivariate_negative
                                    ) {
                                        printed_first_bivariate_negative =
                                            true;
                                        std::cout
                                            << "FIRST_NEGATIVE_BIVARIATE"
                                            << " half_level=" << half_level
                                            << " prefix=" << prefix
                                            << " truncation="
                                                << truncation
                                            << " first_target="
                                                << first_target
                                            << " second_target="
                                                << second_target
                                            << " value=" << bivariate
                                            << '\n';
                                    }
                                }
                            }
                        }
                    }

                    const int derivative_degree =
                        m + 1 - half_level;
                    if (derivative_degree < 0) {
                        continue;
                    }
                    std::vector<Integer> derivative(
                        static_cast<std::size_t>(
                            derivative_degree + 1
                        )
                    );
                    for (int power = half_level;
                         power <= m + 1;
                         ++power) {
                        derivative[
                            static_cast<std::size_t>(
                                power - half_level
                            )
                        ] = polynomial[
                            static_cast<std::size_t>(power)
                        ] * falling(power, half_level);
                    }

                    // Substitute x=4u-1.  The Q=1 spectrum lies in
                    // [-1,3].
                    std::vector<Integer> unit_power(
                        static_cast<std::size_t>(
                            derivative_degree + 1
                        )
                    );
                    for (int source_power = 0;
                         source_power <= derivative_degree;
                         ++source_power) {
                        for (int unit_degree = 0;
                             unit_degree <= source_power;
                             ++unit_degree) {
                            Integer term = derivative[
                                static_cast<std::size_t>(source_power)
                            ] * binomial(source_power, unit_degree)
                                * power_of_four(unit_degree);
                            if (
                                ((source_power - unit_degree) & 1)
                                != 0
                            ) {
                                term = -term;
                            }
                            unit_power[
                                static_cast<std::size_t>(unit_degree)
                            ] += term;
                        }
                    }

                    // For p(u)=sum_i a_i u^i of degree d, the
                    // degree-d Bernstein coefficient beta_k times d! is
                    // sum_{i<=k} a_i (k)_i (d-i)!.
                    for (int index = 0;
                         index <= derivative_degree;
                         ++index) {
                        Integer scaled_bernstein = 0;
                        for (int unit_degree = 0;
                             unit_degree <= index;
                             ++unit_degree) {
                            scaled_bernstein += unit_power[
                                static_cast<std::size_t>(unit_degree)
                            ] * falling(index, unit_degree)
                                * factorial(
                                    derivative_degree - unit_degree
                                );
                        }
                        ++bernstein_coordinates;
                        if (scaled_bernstein < 0) {
                            ++negative_bernstein_coordinates;
                            if (truncation >= 4) {
                                ++hard_negative_bernstein_coordinates;
                                if (
                                    !printed_first_hard_bernstein_negative
                                ) {
                                    printed_first_hard_bernstein_negative =
                                        true;
                                    std::cout
                                        << "FIRST_HARD_NEGATIVE_BERNSTEIN"
                                        << " half_level=" << half_level
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " derivative_order="
                                            << half_level
                                        << " index=" << index
                                        << " degree=" << derivative_degree
                                        << " value=" << scaled_bernstein
                                        << '\n';
                                }
                            }
                            if (!printed_first_bernstein_negative) {
                                printed_first_bernstein_negative = true;
                                std::cout
                                    << "FIRST_NEGATIVE_BERNSTEIN"
                                    << " half_level=" << half_level
                                    << " prefix=" << prefix
                                    << " truncation=" << truncation
                                    << " derivative_order="
                                        << half_level
                                    << " index=" << index
                                    << " degree=" << derivative_degree
                                    << " value=" << scaled_bernstein
                                    << '\n';
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_Q1_DIVIDED_DIFFERENCE"
            << " maximum_half_level=" << maximum_half_level
            << " maximum_prefix=" << maximum_prefix
            << " rows=" << rows
            << " bernstein_coordinates=" << bernstein_coordinates
            << " negative_bernstein_coordinates="
                << negative_bernstein_coordinates
            << " hard_negative_bernstein_coordinates="
                << hard_negative_bernstein_coordinates
            << " negative_targets=" << negative_targets
            << " bivariate_coordinates=" << bivariate_coordinates
            << " negative_bivariate_coordinates="
                << negative_bivariate_coordinates
            << " result="
                << (
                    negative_targets == 0U
                        ? "PASS_Q1_TARGET_DISCOVERY"
                        : "FAIL_Q1_TARGET"
                )
            << '\n';
        return negative_targets == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
