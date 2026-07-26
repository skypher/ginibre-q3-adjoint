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

bool fuses(int level, int first, int second, int output) {
    return std::abs(first - second) <= output
        && output <= std::min(
            first + second,
            2 * level - first - second
        )
        && ((first + second + output) & 1) == 0;
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
        for (int output = 0; output <= level; ++output) {
            if (fuses(level, label, source, output)) {
                next[static_cast<std::size_t>(output)] += coefficient;
            }
        }
    }
    state.swap(next);
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

struct Moments {
    std::vector<Integer> closed;
    std::vector<Integer> wall;
};

Moments moments(int level, int label, int maximum_power) {
    Moments result{
        std::vector<Integer>(
            static_cast<std::size_t>(maximum_power + 1)
        ),
        std::vector<Integer>(
            static_cast<std::size_t>(maximum_power + 1)
        )
    };
    std::vector<Integer> state(
        static_cast<std::size_t>(level + 1)
    );
    state[0] = 1;
    for (int power = 0; power <= maximum_power; ++power) {
        result.closed[static_cast<std::size_t>(power)] = state[0];
        result.wall[static_cast<std::size_t>(power)] =
            state[static_cast<std::size_t>(level)];
        if (power != maximum_power) {
            multiply(level, label, state);
        }
    }
    return result;
}

struct PrefixResult {
    Integer packet;
    Integer core;
    std::vector<Integer> determinants;
    std::vector<Integer> partial_cores;
};

PrefixResult prefix(const Moments& value, int prefix_length) {
    const int odd_power = 2 * prefix_length + 1;
    const int total_power = odd_power + 1;
    PrefixResult result{
        value.wall[static_cast<std::size_t>(2 * prefix_length + 2)]
            + binomial(odd_power, 2)
                * (
                    value.wall[
                        static_cast<std::size_t>(2 * prefix_length)
                    ]
                    - value.wall[
                        static_cast<std::size_t>(
                            2 * prefix_length - 1
                        )
                    ]
        ),
        0,
        {},
        {}
    };
    result.determinants.reserve(
        static_cast<std::size_t>(prefix_length - 2)
    );
    for (int index = 2; index <= prefix_length - 1; ++index) {
        const int wall_power = 2 * prefix_length + 2 - 2 * index;
        const int closed_power = 2 * index;
        const Integer determinant =
            value.wall[static_cast<std::size_t>(wall_power)]
                * value.closed[static_cast<std::size_t>(closed_power)]
            - value.wall[static_cast<std::size_t>(wall_power - 1)]
                * value.closed[
                    static_cast<std::size_t>(closed_power + 1)
                ];
        const Integer weight = binomial(odd_power, 2 * index);
        const Integer scaled_weight = total_power * weight;
        if (scaled_weight
                != (total_power - 2 * index)
                    * binomial(total_power, 2 * index)
            || scaled_weight
                != (2 * index + 1)
                    * binomial(total_power, 2 * index + 1)) {
            throw std::logic_error(
                "binomial-flux identity failed"
            );
        }
        result.determinants.push_back(determinant);
        result.core += weight * determinant;
        result.partial_cores.push_back(result.core);
        result.packet += weight * determinant;
    }
    return result;
}

struct Minimum {
    bool initialized = false;
    Integer value = 0;
    int level = 0;
    int label = 0;
    int prefix_length = 0;
    int determinant_index = 0;

    void observe(
        const Integer& candidate,
        int candidate_level,
        int candidate_label,
        int candidate_prefix,
        int candidate_index
    ) {
        if (!initialized || candidate < value) {
            initialized = true;
            value = candidate;
            level = candidate_level;
            label = candidate_label;
            prefix_length = candidate_prefix;
            determinant_index = candidate_index;
        }
    }
};

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 5 && std::string(argv[1]) == "--case") {
            const int level = parse_positive(argv[2], "level");
            const int label = parse_positive(argv[3], "label");
            const int length = parse_positive(argv[4], "prefix");
            if ((level & 1) != 0 || (label & 1) != 0
                || 2 * label >= level || length < 3) {
                throw std::runtime_error(
                    "case requires even level/label, 2*label<level, "
                    "and prefix at least 3"
                );
            }
            const PrefixResult result = prefix(
                moments(level, label, 2 * length + 2),
                length
            );
            std::cout
                << "SU2_SIMPLE_CURRENT_HIERARCHY_CASE"
                << " level=" << level
                << " label=" << label
                << " prefix=" << length
                << " packet=" << result.packet
                << " core=" << result.core << '\n';
            const int odd_power = 2 * length + 1;
            for (std::size_t offset = 0U;
                 offset < result.determinants.size();
                 ++offset) {
                const int index = static_cast<int>(offset) + 2;
                std::cout
                    << "  index=" << index
                    << " weight=" << binomial(odd_power, 2 * index)
                    << " determinant=" << result.determinants[offset]
                    << " partial_core=" << result.partial_cores[offset]
                    << '\n';
            }
            return EXIT_SUCCESS;
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_simple_current_hierarchy "
                "MAXIMUM_LEVEL MAXIMUM_PREFIX | "
                "--case LEVEL EVEN_LABEL PREFIX"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_prefix =
            parse_positive(argv[2], "maximum prefix");
        if (maximum_prefix < 3) {
            throw std::runtime_error("maximum prefix must be at least 3");
        }

        Minimum packet_minimum;
        Minimum core_minimum;
        std::vector<Minimum> prefix_packet_minima(
            static_cast<std::size_t>(maximum_prefix + 1)
        );
        std::vector<Minimum> prefix_core_minima(
            static_cast<std::size_t>(maximum_prefix + 1)
        );
        std::vector<std::vector<Minimum>> determinant_minima(
            static_cast<std::size_t>(maximum_prefix + 1),
            std::vector<Minimum>(
                static_cast<std::size_t>(maximum_prefix + 1)
            )
        );
        std::vector<std::vector<Minimum>> partial_core_minima(
            static_cast<std::size_t>(maximum_prefix + 1),
            std::vector<Minimum>(
                static_cast<std::size_t>(maximum_prefix + 1)
            )
        );
        std::size_t rows = 0U;
        std::size_t partial_core_rows = 0U;
        std::size_t single_crossing_rows = 0U;
        std::size_t sign_recrossings = 0U;
        int first_recrossing_level = 0;
        int first_recrossing_label = 0;
        int first_recrossing_prefix = 0;
        int first_recrossing_index = 0;
        for (int level = 4; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                const Moments value = moments(
                    level,
                    label,
                    2 * maximum_prefix + 2
                );
                for (int length = 3;
                     length <= maximum_prefix;
                     ++length) {
                    const PrefixResult result = prefix(value, length);
                    ++rows;
                    bool negative_seen = false;
                    bool recrossing = false;
                    for (std::size_t offset = 0U;
                         offset < result.determinants.size();
                         ++offset) {
                        const Integer& determinant =
                            result.determinants[offset];
                        if (determinant < 0) {
                            negative_seen = true;
                        } else if (determinant > 0 && negative_seen) {
                            recrossing = true;
                            if (sign_recrossings == 0U) {
                                first_recrossing_level = level;
                                first_recrossing_label = label;
                                first_recrossing_prefix = length;
                                first_recrossing_index =
                                    static_cast<int>(offset) + 2;
                            }
                            ++sign_recrossings;
                        }
                    }
                    if (!recrossing) {
                        ++single_crossing_rows;
                    }
                    packet_minimum.observe(
                        result.packet, level, label, length, 0
                    );
                    core_minimum.observe(
                        result.core, level, label, length, 0
                    );
                    prefix_packet_minima[
                        static_cast<std::size_t>(length)
                    ].observe(
                        result.packet, level, label, length, 0
                    );
                    prefix_core_minima[
                        static_cast<std::size_t>(length)
                    ].observe(
                        result.core, level, label, length, 0
                    );
                    for (std::size_t offset = 0U;
                         offset < result.determinants.size();
                         ++offset) {
                        const int determinant_index =
                            static_cast<int>(offset) + 2;
                        determinant_minima[
                            static_cast<std::size_t>(length)
                        ][
                            static_cast<std::size_t>(determinant_index)
                        ].observe(
                            result.determinants[offset],
                            level,
                            label,
                            length,
                            determinant_index
                        );
                        partial_core_minima[
                            static_cast<std::size_t>(length)
                        ][
                            static_cast<std::size_t>(determinant_index)
                        ].observe(
                            result.partial_cores[offset],
                            level,
                            label,
                            length,
                            determinant_index
                        );
                        if (length >= 4) {
                            ++partial_core_rows;
                            if (result.partial_cores[offset] < 0) {
                                std::cout
                                    << "SU2_SIMPLE_CURRENT_HIERARCHY"
                                    << " PARTIAL_CORE_COUNTEREXAMPLE"
                                    << " level=" << level
                                    << " label=" << label
                                    << " prefix=" << length
                                    << " truncation="
                                    << determinant_index
                                    << " value="
                                    << result.partial_cores[offset]
                                    << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                    }
                    if (result.packet < 0) {
                        std::cout
                            << "SU2_SIMPLE_CURRENT_HIERARCHY"
                            << " COUNTEREXAMPLE level=" << level
                            << " label=" << label
                            << " prefix=" << length
                            << " packet=" << result.packet
                            << " core=" << result.core
                            << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
        }

        std::cout
            << "SU2_SIMPLE_CURRENT_HIERARCHY"
            << " rows=" << rows
            << " partial_core_rows=" << partial_core_rows
            << " single_crossing_rows=" << single_crossing_rows
            << " sign_recrossings=" << sign_recrossings
            << " first_recrossing=(" << first_recrossing_level
            << ',' << first_recrossing_label
            << ',' << first_recrossing_prefix
            << ',' << first_recrossing_index << ')'
            << " maximum_level=" << maximum_level
            << " maximum_prefix=" << maximum_prefix
            << " minimum_packet=" << packet_minimum.value
            << " packet_witness=(" << packet_minimum.level
            << ',' << packet_minimum.label
            << ',' << packet_minimum.prefix_length << ')'
            << " minimum_core=" << core_minimum.value
            << " core_witness=(" << core_minimum.level
            << ',' << core_minimum.label
            << ',' << core_minimum.prefix_length << ')'
            << " result=PASS_DISCOVERY\n";
        for (int length = 3; length <= maximum_prefix; ++length) {
            const Minimum& packet =
                prefix_packet_minima[static_cast<std::size_t>(length)];
            const Minimum& core =
                prefix_core_minima[static_cast<std::size_t>(length)];
            std::cout
                << "PREFIX length=" << length
                << " minimum_packet=" << packet.value
                << " packet_witness=(" << packet.level
                << ',' << packet.label << ')'
                << " minimum_core=" << core.value
                << " core_witness=(" << core.level
                << ',' << core.label << ")\n";
            for (int index = 2; index <= length - 1; ++index) {
                const Minimum& minimum = determinant_minima[
                    static_cast<std::size_t>(length)
                ][static_cast<std::size_t>(index)];
                const Minimum& partial = partial_core_minima[
                    static_cast<std::size_t>(length)
                ][static_cast<std::size_t>(index)];
                std::cout
                    << "  D index=" << index
                    << " minimum=" << minimum.value
                    << " witness=(" << minimum.level
                    << ',' << minimum.label << ')'
                    << " partial_core=" << partial.value
                    << " partial_witness=(" << partial.level
                    << ',' << partial.label << ")\n";
            }
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_SIMPLE_CURRENT_HIERARCHY FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
