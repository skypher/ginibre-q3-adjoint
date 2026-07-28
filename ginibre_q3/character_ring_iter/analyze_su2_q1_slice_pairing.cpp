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

std::vector<Integer> multiply_fully_looped_path(
    const std::vector<Integer>& state
) {
    std::vector<Integer> next(state.size());
    if (state.size() == 1U) {
        next[0] = state[0];
        return next;
    }
    next[0] = state[0] + state[1];
    const std::size_t last = state.size() - 1U;
    next[last] = state[last - 1U] + state[last];
    for (std::size_t vertex = 1U; vertex < last; ++vertex) {
        next[vertex] =
            state[vertex - 1U] + state[vertex] + state[vertex + 1U];
    }
    return next;
}

Integer power_of_two(int exponent) {
    Integer result = 1;
    for (int index = 0; index < exponent; ++index) {
        result *= 2;
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 4) {
            throw std::runtime_error(
                "usage: analyze_su2_q1_slice_pairing "
                "MAXIMUM_HALF_LEVEL MAXIMUM_PREFIX [--frontier]"
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
        const bool print_frontier =
            argc == 4 && std::string(argv[3]) == "--frontier";
        if (argc == 4 && !print_frontier) {
            throw std::runtime_error(
                "the optional third argument must be --frontier"
            );
        }

        std::uint64_t profiles = 0U;
        std::uint64_t slice_coordinates = 0U;
        std::uint64_t negative_slices = 0U;
        std::uint64_t negative_even_level_slices = 0U;
        std::uint64_t negative_odd_level_slices = 0U;
        std::uint64_t negative_above_three_slices = 0U;
        std::uint64_t sign_recrossings = 0U;
        std::uint64_t prefix_coordinates = 0U;
        std::uint64_t negative_prefixes = 0U;
        std::uint64_t mirror_pairs = 0U;
        std::uint64_t negative_mirror_pairs = 0U;
        std::uint64_t terminal_negatives = 0U;
        std::uint64_t return_ratio_minors = 0U;
        std::uint64_t negative_return_ratio_minors = 0U;
        std::uint64_t endpoint_ratio_minors = 0U;
        std::uint64_t negative_endpoint_ratio_minors = 0U;
        std::uint64_t wrong_parity_endpoint_ratio_minors = 0U;
        std::uint64_t positive_above_three_endpoint_ratio_minors = 0U;
        std::uint64_t negative_level_three_endpoint_ratio_minors = 0U;
        std::uint64_t level_three_identities = 0U;
        std::uint64_t failed_level_three_identities = 0U;
        std::uint64_t incidence_identities = 0U;
        std::uint64_t failed_incidence_identities = 0U;
        std::uint64_t fully_looped_log_concavity_minors = 0U;
        std::uint64_t negative_fully_looped_log_concavity_minors = 0U;
        std::uint64_t fully_looped_turan_reserve_steps = 0U;
        std::uint64_t negative_fully_looped_turan_reserve_steps = 0U;
        std::uint64_t incidence_ratio_sandwiches = 0U;
        std::uint64_t wrong_above_three_incidence_ratio_sandwiches = 0U;
        std::uint64_t negative_even_above_three_incidence_ratio_sandwiches =
            0U;
        std::uint64_t positive_odd_above_three_incidence_ratio_sandwiches =
            0U;
        std::uint64_t negative_above_three_incidence_ratio_sandwiches = 0U;
        int maximum_sign_changes = 0;
        bool printed_first_mirror_negative = false;
        bool printed_first_recrossing = false;
        bool printed_first_endpoint_ratio_negative = false;
        bool printed_first_negative_slice = false;
        bool printed_first_wrong_incidence_ratio_sandwich = false;
        bool printed_first_negative_incidence_ratio_sandwich = false;

        const int maximum_degree = 2 * maximum_prefix + 2;
        for (int half_level = 3;
             half_level <= maximum_half_level;
             ++half_level) {
            std::vector<Integer> returns(
                static_cast<std::size_t>(maximum_degree + 1)
            );
            std::vector<Integer> endpoints(
                static_cast<std::size_t>(maximum_degree + 1)
            );
            std::vector<Integer> near_returns(
                static_cast<std::size_t>(maximum_degree + 1)
            );
            std::vector<Integer> near_endpoints(
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
                near_returns[static_cast<std::size_t>(degree)] = state[1];
                endpoints[static_cast<std::size_t>(degree)] =
                    state[static_cast<std::size_t>(half_level)];
                near_endpoints[static_cast<std::size_t>(degree)] =
                    state[
                        static_cast<std::size_t>(half_level - 1)
                    ];
            }
            std::vector<Integer> fully_looped_endpoints(
                static_cast<std::size_t>(maximum_degree + 1)
            );
            std::vector<Integer> fully_looped_state(
                static_cast<std::size_t>(half_level)
            );
            fully_looped_state[0] = 1;
            fully_looped_endpoints[0] =
                fully_looped_state[
                    static_cast<std::size_t>(half_level - 1)
                ];
            for (int degree = 1;
                 degree <= maximum_degree;
                 ++degree) {
                fully_looped_state =
                    multiply_fully_looped_path(fully_looped_state);
                fully_looped_endpoints[
                    static_cast<std::size_t>(degree)
                ] = fully_looped_state[
                    static_cast<std::size_t>(half_level - 1)
                ];
            }
            for (int degree = 1;
                 degree <= maximum_degree;
                 ++degree) {
                ++incidence_identities;
                if (
                    endpoints[static_cast<std::size_t>(degree)]
                        + endpoints[
                            static_cast<std::size_t>(degree - 1)
                        ]
                    != fully_looped_endpoints[
                        static_cast<std::size_t>(degree - 1)
                    ]
                ) {
                    ++failed_incidence_identities;
                }
            }
            for (int degree = 1;
                 degree + 1 <= maximum_degree;
                 ++degree) {
                if (
                    fully_looped_endpoints[
                        static_cast<std::size_t>(degree - 1)
                    ] == 0
                ) {
                    continue;
                }
                const Integer minor =
                    fully_looped_endpoints[
                        static_cast<std::size_t>(degree)
                    ] * fully_looped_endpoints[
                        static_cast<std::size_t>(degree)
                    ]
                    - fully_looped_endpoints[
                        static_cast<std::size_t>(degree - 1)
                    ] * fully_looped_endpoints[
                        static_cast<std::size_t>(degree + 1)
                    ];
                ++fully_looped_log_concavity_minors;
                if (minor < 0) {
                    ++negative_fully_looped_log_concavity_minors;
                }
            }
            for (int degree = half_level > 3
                    ? half_level
                    : maximum_degree + 1;
                 degree + 2 <= maximum_degree;
                 ++degree) {
                const Integer previous =
                    fully_looped_endpoints[
                        static_cast<std::size_t>(degree - 1)
                    ];
                const Integer current =
                    fully_looped_endpoints[
                        static_cast<std::size_t>(degree)
                    ];
                const Integer next =
                    fully_looped_endpoints[
                        static_cast<std::size_t>(degree + 1)
                    ];
                const Integer after_next =
                    fully_looped_endpoints[
                        static_cast<std::size_t>(degree + 2)
                    ];
                const Integer current_turan =
                    current * current - previous * next;
                const Integer next_turan =
                    next * next - current * after_next;
                const Integer reserve_step =
                    next * next_turan * (previous + current)
                    - current * current_turan * (next + after_next);
                ++fully_looped_turan_reserve_steps;
                if (reserve_step < 0) {
                    ++negative_fully_looped_turan_reserve_steps;
                }
            }
            for (int degree = 1;
                 degree <= maximum_degree;
                 ++degree) {
                if (
                    endpoints[static_cast<std::size_t>(degree - 1)] == 0
                    || fully_looped_endpoints[
                        static_cast<std::size_t>(degree - 1)
                    ] == 0
                ) {
                    continue;
                }
                const Integer minor =
                    endpoints[static_cast<std::size_t>(degree)]
                        * fully_looped_endpoints[
                            static_cast<std::size_t>(degree - 1)
                        ]
                    - endpoints[static_cast<std::size_t>(degree - 1)]
                        * fully_looped_endpoints[
                            static_cast<std::size_t>(degree)
                        ];
                ++incidence_ratio_sandwiches;
                if (
                    half_level > 3
                    && (
                        (degree % 2 == 0 && minor < 0)
                        || (degree % 2 != 0 && minor > 0)
                    )
                ) {
                    ++wrong_above_three_incidence_ratio_sandwiches;
                    if (!printed_first_wrong_incidence_ratio_sandwich) {
                        printed_first_wrong_incidence_ratio_sandwich = true;
                        std::cout
                            << "FIRST_WRONG_INCIDENCE_RATIO_SANDWICH"
                            << " half_level=" << half_level
                            << " degree=" << degree
                            << " value=" << minor << '\n';
                    }
                }
                if (
                    half_level > 3 && degree % 2 == 0 && minor < 0
                ) {
                    ++negative_even_above_three_incidence_ratio_sandwiches;
                }
                if (
                    half_level > 3 && degree % 2 != 0 && minor > 0
                ) {
                    ++positive_odd_above_three_incidence_ratio_sandwiches;
                }
                if (half_level > 3 && minor < 0) {
                    ++negative_above_three_incidence_ratio_sandwiches;
                    if (!printed_first_negative_incidence_ratio_sandwich) {
                        printed_first_negative_incidence_ratio_sandwich =
                            true;
                        std::cout
                            << "FIRST_NEGATIVE_INCIDENCE_RATIO_SANDWICH"
                            << " half_level=" << half_level
                            << " degree=" << degree
                            << " value=" << minor << '\n';
                    }
                }
            }
            for (int even = 2;
                 even + 2 <= maximum_degree;
                 even += 2) {
                const Integer original_minor =
                    endpoints[static_cast<std::size_t>(even + 1)]
                        * endpoints[static_cast<std::size_t>(even)]
                    - endpoints[static_cast<std::size_t>(even + 2)]
                        * endpoints[static_cast<std::size_t>(even - 1)];
                const Integer reduced_minor =
                    endpoints[static_cast<std::size_t>(even + 1)]
                        * fully_looped_endpoints[
                            static_cast<std::size_t>(even - 1)
                        ]
                    - endpoints[static_cast<std::size_t>(even - 1)]
                        * fully_looped_endpoints[
                            static_cast<std::size_t>(even + 1)
                        ];
                ++incidence_identities;
                if (original_minor != reduced_minor) {
                    ++failed_incidence_identities;
                }
            }
            if (half_level == 3) {
                for (int ratio_index = 1;
                     2 * ratio_index + 1 <= maximum_degree;
                     ++ratio_index) {
                    const int even = 2 * ratio_index;
                    const Integer c =
                        returns[static_cast<std::size_t>(even)]
                        + endpoints[static_cast<std::size_t>(even)];
                    const Integer p =
                        2 * returns[
                            static_cast<std::size_t>(even + 1)
                        ];
                    const Integer previous_p =
                        ratio_index == 1
                            ? Integer(0)
                            : 2 * returns[
                                static_cast<std::size_t>(even - 1)
                            ];
                    ++level_three_identities;
                    if (c * c - 1 != p * previous_p) {
                        ++failed_level_three_identities;
                    }
                }
            }
            for (int slice = 1;
                 2 * slice + 3 <= maximum_degree;
                 ++slice) {
                const int even = 2 * slice;
                const Integer minor =
                    returns[static_cast<std::size_t>(even)]
                        * returns[static_cast<std::size_t>(even + 3)]
                    - returns[static_cast<std::size_t>(even + 1)]
                        * returns[static_cast<std::size_t>(even + 2)];
                ++return_ratio_minors;
                if (minor < 0) {
                    ++negative_return_ratio_minors;
                }
            }
            for (int ratio_index = 1;
                 2 * ratio_index + 2 <= maximum_degree;
                 ++ratio_index) {
                const int even = 2 * ratio_index;
                if (
                    endpoints[static_cast<std::size_t>(even - 1)] == 0
                    || endpoints[static_cast<std::size_t>(even)] == 0
                ) {
                    continue;
                }
                const Integer minor =
                    endpoints[static_cast<std::size_t>(even + 2)]
                        * endpoints[
                            static_cast<std::size_t>(even - 1)
                        ]
                    - endpoints[static_cast<std::size_t>(even + 1)]
                        * endpoints[static_cast<std::size_t>(even)];
                ++endpoint_ratio_minors;
                if (minor < 0) {
                    ++negative_endpoint_ratio_minors;
                    if (!printed_first_endpoint_ratio_negative) {
                        printed_first_endpoint_ratio_negative = true;
                        std::cout
                            << "FIRST_NEGATIVE_ENDPOINT_RATIO_MINOR"
                            << " half_level=" << half_level
                            << " ratio_index=" << ratio_index
                            << " value=" << minor << '\n';
                    }
                }
                if (
                    (half_level % 2 == 1 && minor < 0)
                    || (half_level % 2 == 0 && minor > 0)
                ) {
                    ++wrong_parity_endpoint_ratio_minors;
                }
                if (half_level > 3 && minor > 0) {
                    ++positive_above_three_endpoint_ratio_minors;
                }
                if (half_level == 3 && minor < 0) {
                    ++negative_level_three_endpoint_ratio_minors;
                }
            }

            for (int prefix = 4;
                 prefix <= maximum_prefix;
                 ++prefix) {
                ++profiles;
                const int m = 2 * prefix + 1;
                const int n = m + 1;
                std::vector<Integer> weighted_slices(
                    static_cast<std::size_t>(prefix)
                );
                Integer cumulative = 0;
                int previous_sign = 0;
                int sign_changes = 0;
                int first_negative_slice = -1;
                for (int slice = 0;
                     slice < prefix;
                     ++slice) {
                    const int even = 2 * slice;
                    const Integer unweighted =
                        returns[static_cast<std::size_t>(even)]
                            * endpoints[
                                static_cast<std::size_t>(n - even)
                            ]
                        - returns[
                            static_cast<std::size_t>(even + 1)
                        ] * endpoints[
                            static_cast<std::size_t>(n - even - 1)
                        ];
                    const Integer weighted =
                        binomial(m, even) * unweighted;
                    weighted_slices[static_cast<std::size_t>(slice)] =
                        weighted;
                    ++slice_coordinates;
                    if (weighted < 0) {
                        ++negative_slices;
                        if (first_negative_slice < 0) {
                            first_negative_slice = slice;
                        }
                        if (!printed_first_negative_slice) {
                            printed_first_negative_slice = true;
                            std::cout
                                << "FIRST_NEGATIVE_SLICE"
                                << " half_level=" << half_level
                                << " prefix=" << prefix
                                << " slice=" << slice
                                << " value=" << weighted << '\n';
                        }
                        if (half_level % 2 == 0) {
                            ++negative_even_level_slices;
                        } else {
                            ++negative_odd_level_slices;
                        }
                        if (half_level > 3) {
                            ++negative_above_three_slices;
                        }
                    }
                    const int sign =
                        weighted > 0 ? 1 : (weighted < 0 ? -1 : 0);
                    if (sign != 0) {
                        if (previous_sign != 0
                            && sign != previous_sign) {
                            ++sign_changes;
                            if (sign_changes >= 2) {
                                ++sign_recrossings;
                                if (!printed_first_recrossing) {
                                    printed_first_recrossing = true;
                                    std::cout
                                        << "FIRST_SLICE_RECROSSING"
                                        << " half_level=" << half_level
                                        << " prefix=" << prefix
                                        << " slice=" << slice
                                        << " value=" << weighted << '\n';
                                }
                            }
                        }
                        previous_sign = sign;
                    }
                    cumulative += weighted;
                    ++prefix_coordinates;
                    if (cumulative < 0) {
                        ++negative_prefixes;
                    }
                }
                maximum_sign_changes =
                    std::max(maximum_sign_changes, sign_changes);
                if (cumulative < 0) {
                    ++terminal_negatives;
                }
                if (half_level == 3) {
                    Integer unrestricted_difference = 0;
                    for (int second_steps = 0;
                         second_steps <= m;
                         ++second_steps) {
                        unrestricted_difference +=
                            binomial(m, second_steps) * (
                                near_endpoints[
                                    static_cast<std::size_t>(
                                        m - second_steps
                                    )
                                ] * returns[
                                    static_cast<std::size_t>(second_steps)
                                ]
                                - endpoints[
                                    static_cast<std::size_t>(
                                        m - second_steps
                                    )
                                ] * near_returns[
                                    static_cast<std::size_t>(second_steps)
                                ]
                            );
                    }
                    Integer signed_difference = 0;
                    for (int second_steps = 0;
                         second_steps <= n;
                         ++second_steps) {
                        Integer term =
                            binomial(n, second_steps)
                            * endpoints[
                                static_cast<std::size_t>(
                                    n - second_steps
                                )
                            ] * returns[
                                static_cast<std::size_t>(second_steps)
                            ];
                        if (second_steps % 2 != 0) {
                            term = -term;
                        }
                        signed_difference += term;
                    }
                    const Integer signed_formula =
                        power_of_two(3 * n / 2 - 4)
                        - power_of_two(n - 3);
                    ++level_three_identities;
                    if (
                        2 * cumulative
                            != unrestricted_difference + signed_difference
                        || unrestricted_difference < 0
                        || signed_difference != signed_formula
                        || signed_formula <= 0
                    ) {
                        ++failed_level_three_identities;
                    }
                }
                if (print_frontier && prefix == maximum_prefix) {
                    std::cout
                        << "SLICE_FRONTIER"
                        << " half_level=" << half_level
                        << " prefix=" << prefix
                        << " first_negative_slice="
                            << first_negative_slice
                        << " terminal=" << cumulative << '\n';
                }

                for (int lower = 0;
                     lower < prefix - 1 - lower;
                     ++lower) {
                    const int upper = prefix - 1 - lower;
                    const Integer paired =
                        weighted_slices[
                            static_cast<std::size_t>(lower)
                        ] + weighted_slices[
                            static_cast<std::size_t>(upper)
                        ];
                    ++mirror_pairs;
                    if (paired < 0) {
                        ++negative_mirror_pairs;
                        if (!printed_first_mirror_negative) {
                            printed_first_mirror_negative = true;
                            std::cout
                                << "FIRST_NEGATIVE_MIRROR_PAIR"
                                << " half_level=" << half_level
                                << " prefix=" << prefix
                                << " lower_slice=" << lower
                                << " upper_slice=" << upper
                                << " value=" << paired << '\n';
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_Q1_SLICE_PAIRING"
            << " maximum_half_level=" << maximum_half_level
            << " maximum_prefix=" << maximum_prefix
            << " profiles=" << profiles
            << " slice_coordinates=" << slice_coordinates
            << " negative_slices=" << negative_slices
            << " negative_even_level_slices="
                << negative_even_level_slices
            << " negative_odd_level_slices="
                << negative_odd_level_slices
            << " negative_above_three_slices="
                << negative_above_three_slices
            << " sign_recrossings=" << sign_recrossings
            << " maximum_sign_changes=" << maximum_sign_changes
            << " prefix_coordinates=" << prefix_coordinates
            << " negative_prefixes=" << negative_prefixes
            << " mirror_pairs=" << mirror_pairs
            << " negative_mirror_pairs=" << negative_mirror_pairs
            << " terminal_negatives=" << terminal_negatives
            << " return_ratio_minors=" << return_ratio_minors
            << " negative_return_ratio_minors="
                << negative_return_ratio_minors
            << " endpoint_ratio_minors=" << endpoint_ratio_minors
            << " negative_endpoint_ratio_minors="
                << negative_endpoint_ratio_minors
            << " wrong_parity_endpoint_ratio_minors="
                << wrong_parity_endpoint_ratio_minors
            << " positive_above_three_endpoint_ratio_minors="
                << positive_above_three_endpoint_ratio_minors
            << " negative_level_three_endpoint_ratio_minors="
                << negative_level_three_endpoint_ratio_minors
            << " level_three_identities=" << level_three_identities
            << " failed_level_three_identities="
                << failed_level_three_identities
            << " incidence_identities=" << incidence_identities
            << " failed_incidence_identities="
                << failed_incidence_identities
            << " fully_looped_log_concavity_minors="
                << fully_looped_log_concavity_minors
            << " negative_fully_looped_log_concavity_minors="
                << negative_fully_looped_log_concavity_minors
            << " fully_looped_turan_reserve_steps="
                << fully_looped_turan_reserve_steps
            << " negative_fully_looped_turan_reserve_steps="
                << negative_fully_looped_turan_reserve_steps
            << " incidence_ratio_sandwiches="
                << incidence_ratio_sandwiches
            << " wrong_above_three_incidence_ratio_sandwiches="
                << wrong_above_three_incidence_ratio_sandwiches
            << " negative_even_above_three_incidence_ratio_sandwiches="
                << negative_even_above_three_incidence_ratio_sandwiches
            << " positive_odd_above_three_incidence_ratio_sandwiches="
                << positive_odd_above_three_incidence_ratio_sandwiches
            << " negative_above_three_incidence_ratio_sandwiches="
                << negative_above_three_incidence_ratio_sandwiches
            << " result="
                << (
                    sign_recrossings == 0U
                        && negative_prefixes == 0U
                        && negative_mirror_pairs > 0U
                        && terminal_negatives == 0U
                        && negative_above_three_slices == 0U
                        && negative_return_ratio_minors == 0U
                        && positive_above_three_endpoint_ratio_minors == 0U
                        && negative_level_three_endpoint_ratio_minors == 0U
                        && failed_level_three_identities == 0U
                        && failed_incidence_identities == 0U
                        && negative_fully_looped_log_concavity_minors == 0U
                        && negative_fully_looped_turan_reserve_steps == 0U
                        && negative_even_above_three_incidence_ratio_sandwiches
                            == 0U
                        ? "PASS_Q1_SLICE_PAIRING_DISCOVERY"
                        : "FAIL_Q1_SLICE_PAIRING"
                )
            << '\n';
        return sign_recrossings == 0U
                && negative_prefixes == 0U
                && negative_mirror_pairs > 0U
                && terminal_negatives == 0U
                && negative_above_three_slices == 0U
                && negative_return_ratio_minors == 0U
                && positive_above_three_endpoint_ratio_minors == 0U
                && negative_level_three_endpoint_ratio_minors == 0U
                && failed_level_three_identities == 0U
                && failed_incidence_identities == 0U
                && negative_fully_looped_log_concavity_minors == 0U
                && negative_fully_looped_turan_reserve_steps == 0U
                && negative_even_above_three_incidence_ratio_sandwiches
                    == 0U
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
