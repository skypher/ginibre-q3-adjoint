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
    int lower = 0;
    int upper = 0;
    int vertex = 0;

    void observe(
        const Integer& candidate,
        int candidate_level,
        int candidate_label,
        int candidate_prefix,
        int candidate_lower,
        int candidate_upper,
        int candidate_vertex
    ) {
        if (!initialized || candidate < value) {
            initialized = true;
            value = candidate;
            level = candidate_level;
            label = candidate_label;
            prefix = candidate_prefix;
            lower = candidate_lower;
            upper = candidate_upper;
            vertex = candidate_vertex;
        }
    }
};

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_reflected_prefix_pairs "
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

        Minimum pair_minimum;
        Minimum cumulative_minimum;
        Minimum post_three_minimum;
        Minimum lower_minimum;
        std::uint64_t parameters = 0U;
        std::uint64_t pairs = 0U;
        std::uint64_t coordinates = 0U;
        std::uint64_t negative_pairs = 0U;
        std::uint64_t negative_cumulative = 0U;
        std::uint64_t negative_lower = 0U;
        std::uint64_t negative_post_three = 0U;
        std::uint64_t positive_upper = 0U;
        std::uint64_t negative_upper = 0U;

        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameters;
                std::vector<std::vector<Integer>> powers(
                    static_cast<std::size_t>(2 * maximum_prefix + 3)
                );
                powers[0].assign(
                    static_cast<std::size_t>(level + 1),
                    Integer(0)
                );
                powers[0][0] = 1;
                for (int power = 1;
                     power <= 2 * maximum_prefix + 2;
                     ++power) {
                    powers[static_cast<std::size_t>(power)] =
                        multiply(
                            level,
                            label,
                            powers[static_cast<std::size_t>(power - 1)]
                        );
                }

                for (int prefix = 4;
                     prefix <= maximum_prefix;
                     ++prefix) {
                    const int total_length = 2 * prefix + 2;
                    std::vector<Integer> cumulative(
                        static_cast<std::size_t>(level + 1)
                    );
                    std::vector<Integer> post_three(
                        static_cast<std::size_t>(level + 1)
                    );
                    for (int slice = 4;
                         slice <= prefix - 1;
                         ++slice) {
                        const Integer weight =
                            binomial(total_length - 1, 2 * slice);
                        for (int vertex = 0;
                             vertex <= level;
                             vertex += 2) {
                            post_three[
                                static_cast<std::size_t>(vertex)
                            ] += weight
                                * (
                                    powers[static_cast<std::size_t>(
                                        2 * slice
                                    )][0]
                                        * powers[
                                            static_cast<std::size_t>(
                                                total_length - 2 * slice
                                            )
                                        ][static_cast<std::size_t>(
                                            vertex
                                        )]
                                    - powers[static_cast<std::size_t>(
                                        2 * slice + 1
                                    )][0]
                                        * powers[
                                            static_cast<std::size_t>(
                                                total_length
                                                - 2 * slice - 1
                                            )
                                        ][static_cast<std::size_t>(
                                            vertex
                                        )]
                                );
                            post_three_minimum.observe(
                                post_three[
                                    static_cast<std::size_t>(vertex)
                                ],
                                level,
                                label,
                                prefix,
                                3,
                                slice,
                                vertex
                            );
                            if (post_three[
                                    static_cast<std::size_t>(vertex)
                                ] < 0) {
                                ++negative_post_three;
                                if (negative_post_three == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_POST_THREE"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << slice
                                        << " vertex=" << vertex
                                        << " tail="
                                        << post_three[
                                            static_cast<std::size_t>(
                                                vertex
                                            )
                                        ] << '\n';
                                }
                            }
                        }
                    }
                    for (int upper = prefix / 2 + 1;
                         upper <= prefix - 1;
                         ++upper) {
                        const int lower = prefix - upper;
                        ++pairs;
                        for (int vertex = 0;
                             vertex <= level;
                             vertex += 2) {
                            ++coordinates;
                            const Integer lower_slice =
                                powers[static_cast<std::size_t>(
                                    2 * lower
                                )][0]
                                    * powers[static_cast<std::size_t>(
                                        2 * upper + 2
                                    )][static_cast<std::size_t>(vertex)]
                                - powers[static_cast<std::size_t>(
                                    2 * lower + 1
                                )][0]
                                    * powers[static_cast<std::size_t>(
                                        2 * upper + 1
                                    )][static_cast<std::size_t>(vertex)];
                            const Integer upper_slice =
                                powers[static_cast<std::size_t>(
                                    2 * upper
                                )][0]
                                    * powers[static_cast<std::size_t>(
                                        2 * lower + 2
                                    )][static_cast<std::size_t>(vertex)]
                                - powers[static_cast<std::size_t>(
                                    2 * upper + 1
                                )][0]
                                    * powers[static_cast<std::size_t>(
                                        2 * lower + 1
                                    )][static_cast<std::size_t>(vertex)];
                            const Integer paired =
                                (2 * lower + 1) * lower_slice
                                + (2 * upper + 1) * upper_slice;
                            cumulative[static_cast<std::size_t>(vertex)]
                                += paired;
                            pair_minimum.observe(
                                paired,
                                level,
                                label,
                                prefix,
                                lower,
                                upper,
                                vertex
                            );
                            lower_minimum.observe(
                                lower_slice,
                                level,
                                label,
                                prefix,
                                lower,
                                upper,
                                vertex
                            );
                            cumulative_minimum.observe(
                                cumulative[
                                    static_cast<std::size_t>(vertex)
                                ],
                                level,
                                label,
                                prefix,
                                lower,
                                upper,
                                vertex
                            );
                            if (paired < 0) {
                                ++negative_pairs;
                                if (negative_pairs == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_PAIR"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " lower=" << lower
                                        << " upper=" << upper
                                        << " vertex=" << vertex
                                        << " lower_slice=" << lower_slice
                                        << " upper_slice=" << upper_slice
                                        << " pair=" << paired << '\n';
                                }
                            }
                            if (lower_slice < 0) {
                                ++negative_lower;
                            }
                            if (upper_slice < 0) {
                                ++negative_upper;
                            } else if (upper_slice > 0) {
                                ++positive_upper;
                            }
                            if (cumulative[
                                    static_cast<std::size_t>(vertex)
                                ] < 0) {
                                ++negative_cumulative;
                                if (negative_cumulative == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_CUMULATIVE"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " last_lower=" << lower
                                        << " last_upper=" << upper
                                        << " vertex=" << vertex
                                        << " cumulative="
                                        << cumulative[
                                            static_cast<std::size_t>(
                                                vertex
                                            )
                                        ] << '\n';
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_REFLECTED_PREFIX_PAIRS"
            << " maximum_level=" << maximum_level
            << " maximum_prefix=" << maximum_prefix
            << " parameters=" << parameters
            << " pairs=" << pairs
            << " coordinates=" << coordinates
            << " negative_pairs=" << negative_pairs
            << " negative_cumulative=" << negative_cumulative
            << " negative_lower=" << negative_lower
            << " negative_post_three=" << negative_post_three
            << " positive_upper=" << positive_upper
            << " negative_upper=" << negative_upper
            << " pair_minimum=" << pair_minimum.value
            << " pair_witness=("
            << pair_minimum.level << ','
            << pair_minimum.label << ','
            << pair_minimum.prefix << ','
            << pair_minimum.lower << ','
            << pair_minimum.upper << ','
            << pair_minimum.vertex << ')'
            << " cumulative_minimum=" << cumulative_minimum.value
            << " cumulative_witness=("
            << cumulative_minimum.level << ','
            << cumulative_minimum.label << ','
            << cumulative_minimum.prefix << ','
            << cumulative_minimum.lower << ','
            << cumulative_minimum.upper << ','
            << cumulative_minimum.vertex << ')'
            << " post_three_minimum=" << post_three_minimum.value
            << " post_three_witness=("
            << post_three_minimum.level << ','
            << post_three_minimum.label << ','
            << post_three_minimum.prefix << ','
            << post_three_minimum.upper << ','
            << post_three_minimum.vertex << ')'
            << " lower_minimum=" << lower_minimum.value
            << " lower_witness=("
            << lower_minimum.level << ','
            << lower_minimum.label << ','
            << lower_minimum.prefix << ','
            << lower_minimum.lower << ','
            << lower_minimum.upper << ','
            << lower_minimum.vertex << ')'
            << " result="
            << (negative_cumulative == 0U && negative_lower == 0U
                ? "PASS_POINTWISE_DISCOVERY"
                : "FAIL_POINTWISE")
            << '\n';
        return negative_cumulative == 0U && negative_lower == 0U
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_REFLECTED_PREFIX_PAIRS FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
