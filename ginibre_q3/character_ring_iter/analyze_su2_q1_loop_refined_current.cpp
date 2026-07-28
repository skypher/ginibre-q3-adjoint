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
using LoopPolynomial = std::vector<Integer>;
using EndpointPolynomials = std::vector<LoopPolynomial>;
using PowerTable = std::vector<EndpointPolynomials>;

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

PowerTable build_powers(int half_level, int maximum_degree) {
    PowerTable powers(
        static_cast<std::size_t>(maximum_degree + 1),
        EndpointPolynomials(
            static_cast<std::size_t>(half_level + 1),
            LoopPolynomial(
                static_cast<std::size_t>(maximum_degree + 1)
            )
        )
    );
    powers[0][0][0] = 1;
    for (int degree = 0; degree < maximum_degree; ++degree) {
        for (int vertex = 0; vertex <= half_level; ++vertex) {
            for (int loops = 0; loops <= degree; ++loops) {
                const Integer& value =
                    powers[static_cast<std::size_t>(degree)]
                          [static_cast<std::size_t>(vertex)]
                          [static_cast<std::size_t>(loops)];
                if (value == 0) {
                    continue;
                }
                if (vertex > 0) {
                    powers[static_cast<std::size_t>(degree + 1)]
                          [static_cast<std::size_t>(vertex - 1)]
                          [static_cast<std::size_t>(loops)] += value;
                }
                if (vertex < half_level) {
                    powers[static_cast<std::size_t>(degree + 1)]
                          [static_cast<std::size_t>(vertex + 1)]
                          [static_cast<std::size_t>(loops)] += value;
                }
                if (0 < vertex && vertex < half_level) {
                    powers[static_cast<std::size_t>(degree + 1)]
                          [static_cast<std::size_t>(vertex)]
                          [static_cast<std::size_t>(loops + 1)] += value;
                }
            }
        }
    }
    return powers;
}

void add_product(
    LoopPolynomial& target,
    const LoopPolynomial& left,
    int left_degree,
    const LoopPolynomial& right,
    int right_degree,
    const Integer& weight
) {
    for (int left_loops = 0;
         left_loops <= left_degree;
         ++left_loops) {
        if (left[static_cast<std::size_t>(left_loops)] == 0) {
            continue;
        }
        for (int right_loops = 0;
             right_loops <= right_degree;
             ++right_loops) {
            if (right[static_cast<std::size_t>(right_loops)] == 0) {
                continue;
            }
            target[
                static_cast<std::size_t>(left_loops + right_loops)
            ] += weight
                * left[static_cast<std::size_t>(left_loops)]
                * right[static_cast<std::size_t>(right_loops)];
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_q1_loop_refined_current "
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
        std::uint64_t coordinates = 0U;
        std::uint64_t negative_coordinates = 0U;
        std::uint64_t hard_coordinates = 0U;
        std::uint64_t hard_negative_coordinates = 0U;
        std::uint64_t negative_targets = 0U;
        bool printed_first_negative = false;
        bool printed_first_hard_negative = false;

        const int maximum_degree = 2 * maximum_prefix + 2;
        for (int half_level = 3;
             half_level <= maximum_half_level;
             ++half_level) {
            const PowerTable powers =
                build_powers(half_level, maximum_degree);
            for (int prefix = 4;
                 prefix <= maximum_prefix;
                 ++prefix) {
                const int n = 2 * prefix + 2;
                LoopPolynomial current(
                    static_cast<std::size_t>(n + 1)
                );
                for (int truncation = 0;
                     truncation < prefix;
                     ++truncation) {
                    const int even = 2 * truncation;
                    const int positive_degree = n - even;
                    const int negative_degree = positive_degree - 1;
                    const Integer weight = binomial(n - 1, even);
                    add_product(
                        current,
                        powers[static_cast<std::size_t>(positive_degree)]
                              [static_cast<std::size_t>(half_level)],
                        positive_degree,
                        powers[static_cast<std::size_t>(even)][0],
                        even,
                        weight
                    );
                    add_product(
                        current,
                        powers[static_cast<std::size_t>(negative_degree)]
                              [static_cast<std::size_t>(half_level)],
                        negative_degree,
                        powers[static_cast<std::size_t>(even + 1)][0],
                        even + 1,
                        -weight
                    );

                    ++rows;
                    Integer target = 0;
                    for (int loops = 0; loops <= n; ++loops) {
                        const Integer& value =
                            current[static_cast<std::size_t>(loops)];
                        target += value;
                        ++coordinates;
                        if (truncation >= 4) {
                            ++hard_coordinates;
                        }
                        if (value < 0) {
                            ++negative_coordinates;
                            if (truncation >= 4) {
                                ++hard_negative_coordinates;
                                if (!printed_first_hard_negative) {
                                    printed_first_hard_negative = true;
                                    std::cout
                                        << "FIRST_HARD_NEGATIVE_LOOP_COEFFICIENT"
                                        << " half_level=" << half_level
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " loops=" << loops
                                        << " nonloop_edges=" << n - loops
                                        << " value=" << value << '\n';
                                }
                            }
                            if (!printed_first_negative) {
                                printed_first_negative = true;
                                std::cout
                                    << "FIRST_NEGATIVE_LOOP_COEFFICIENT"
                                    << " half_level=" << half_level
                                    << " prefix=" << prefix
                                    << " truncation=" << truncation
                                    << " loops=" << loops
                                    << " nonloop_edges=" << n - loops
                                    << " value=" << value << '\n';
                            }
                        }
                    }
                    if (target < 0) {
                        ++negative_targets;
                    }
                }
            }
        }

        std::cout
            << "SU2_Q1_LOOP_REFINED_CURRENT"
            << " maximum_half_level=" << maximum_half_level
            << " maximum_prefix=" << maximum_prefix
            << " rows=" << rows
            << " coordinates=" << coordinates
            << " negative_coordinates=" << negative_coordinates
            << " hard_coordinates=" << hard_coordinates
            << " hard_negative_coordinates="
                << hard_negative_coordinates
            << " negative_targets=" << negative_targets
            << " result="
                << (
                    negative_coordinates == 0U
                        ? "PASS_LOOP_REFINED_CURRENT_DISCOVERY"
                        : "FAIL_LOOP_REFINED_CURRENT"
                )
            << '\n';
        return negative_coordinates == 0U
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
