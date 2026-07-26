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

bool fuses(int level, int first, int second, int output) {
    return std::abs(first - second) <= output
        && output <= std::min(
            first + second, 2 * level - first - second
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
}

struct Candidate {
    Integer delta_five;
    Integer paired;
    Integer reverse_core;
    Integer packet;
    Integer reverse_packet;
    Integer eleven_delta_seven;
    Integer eleven_delta_five;
    Integer eleven_delta_three;
    Integer eleven_paired;
    Integer eleven_core;
    Integer eleven_packet;
};

Candidate candidate(int level, int label) {
    std::vector<Integer> state(
        static_cast<std::size_t>(level + 1)
    );
    std::vector<Integer> closed(13U);
    std::vector<Integer> wall(13U);
    state[0] = 1;
    for (int power = 0; power <= 12; ++power) {
        closed[static_cast<std::size_t>(power)] = state[0];
        wall[static_cast<std::size_t>(power)] =
            state[static_cast<std::size_t>(level)];
        if (power != 12) {
            multiply(level, label, state);
        }
    }
    const Integer delta_three =
        wall[4] * closed[6] - wall[3] * closed[7];
    const Integer delta_five =
        wall[6] * closed[4] - wall[5] * closed[5];
    const Integer delta_seven = wall[8] - wall[7];
    const Integer eleven_delta_seven =
        wall[8] * closed[4] - wall[7] * closed[5];
    const Integer eleven_delta_five =
        wall[6] * closed[6] - wall[5] * closed[7];
    const Integer eleven_delta_three =
        wall[4] * closed[8] - wall[3] * closed[9];
    return Candidate{
        delta_five,
        2 * delta_three + 3 * delta_five,
        3 * delta_five
            + 2 * (wall[4] * closed[6] - wall[7]),
        84 * delta_three + 126 * delta_five
            + 36 * delta_seven + wall[10],
        84 * wall[4] * closed[6]
            - 36 * wall[3] * closed[7]
            + 126 * delta_five
            + 36 * wall[8] - 84 * wall[7] + wall[10],
        eleven_delta_seven,
        eleven_delta_five,
        eleven_delta_three,
        14 * eleven_delta_five + 5 * eleven_delta_three,
        10 * eleven_delta_seven
            + 14 * eleven_delta_five + 5 * eleven_delta_three,
        wall[12] + 55 * (wall[10] - wall[9])
            + 330 * eleven_delta_seven
            + 462 * eleven_delta_five
            + 165 * eleven_delta_three
    };
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--scan") {
            const int maximum_level =
                parse_positive(argv[2], "maximum level");
            bool initialized = false;
            Integer minimum = 0;
            Integer minimum_delta_five = 0;
            Integer minimum_reverse = 0;
            Integer minimum_reverse_core = 0;
            Integer minimum_eleven_delta_seven = 0;
            Integer minimum_eleven_delta_five = 0;
            Integer minimum_eleven_delta_three = 0;
            Integer minimum_eleven_packet = 0;
            Integer minimum_eleven_paired = 0;
            Integer minimum_eleven_core = 0;
            int witness_level = 0;
            int witness_label = 0;
            int delta_five_witness_level = 0;
            int delta_five_witness_label = 0;
            int reverse_witness_level = 0;
            int reverse_witness_label = 0;
            int reverse_core_witness_level = 0;
            int reverse_core_witness_label = 0;
            int eleven_seven_witness_level = 0;
            int eleven_seven_witness_label = 0;
            int eleven_five_witness_level = 0;
            int eleven_five_witness_label = 0;
            int eleven_three_witness_level = 0;
            int eleven_three_witness_label = 0;
            int eleven_packet_witness_level = 0;
            int eleven_packet_witness_label = 0;
            int eleven_paired_witness_level = 0;
            int eleven_paired_witness_label = 0;
            int eleven_core_witness_level = 0;
            int eleven_core_witness_label = 0;
            std::size_t rows = 0U;
            for (int level = 4; level <= maximum_level; level += 2) {
                for (int label = 2;
                     2 * label < level; label += 2) {
                    const Candidate result = candidate(level, label);
                    ++rows;
                    if (!initialized || result.paired < minimum) {
                        minimum = result.paired;
                        witness_level = level;
                        witness_label = label;
                    }
                    if (!initialized
                        || result.delta_five < minimum_delta_five) {
                        minimum_delta_five = result.delta_five;
                        delta_five_witness_level = level;
                        delta_five_witness_label = label;
                    }
                    if (!initialized
                        || result.reverse_packet < minimum_reverse) {
                        minimum_reverse = result.reverse_packet;
                        reverse_witness_level = level;
                        reverse_witness_label = label;
                    }
                    if (!initialized
                        || result.reverse_core < minimum_reverse_core) {
                        minimum_reverse_core = result.reverse_core;
                        reverse_core_witness_level = level;
                        reverse_core_witness_label = label;
                    }
                    if (!initialized
                        || result.eleven_delta_seven
                            < minimum_eleven_delta_seven) {
                        minimum_eleven_delta_seven =
                            result.eleven_delta_seven;
                        eleven_seven_witness_level = level;
                        eleven_seven_witness_label = label;
                    }
                    if (!initialized
                        || result.eleven_delta_five
                            < minimum_eleven_delta_five) {
                        minimum_eleven_delta_five =
                            result.eleven_delta_five;
                        eleven_five_witness_level = level;
                        eleven_five_witness_label = label;
                    }
                    if (!initialized
                        || result.eleven_delta_three
                            < minimum_eleven_delta_three) {
                        minimum_eleven_delta_three =
                            result.eleven_delta_three;
                        eleven_three_witness_level = level;
                        eleven_three_witness_label = label;
                    }
                    if (!initialized
                        || result.eleven_packet
                            < minimum_eleven_packet) {
                        minimum_eleven_packet =
                            result.eleven_packet;
                        eleven_packet_witness_level = level;
                        eleven_packet_witness_label = label;
                    }
                    if (!initialized
                        || result.eleven_paired
                            < minimum_eleven_paired) {
                        minimum_eleven_paired =
                            result.eleven_paired;
                        eleven_paired_witness_level = level;
                        eleven_paired_witness_label = label;
                    }
                    if (!initialized
                        || result.eleven_core < minimum_eleven_core) {
                        minimum_eleven_core = result.eleven_core;
                        eleven_core_witness_level = level;
                        eleven_core_witness_label = label;
                    }
                    initialized = true;
                    if (result.paired < 0 || result.packet < 0
                        || result.reverse_packet < 0
                        || result.eleven_packet < 0
                        || result.eleven_core < 0) {
                        std::cout
                            << "SU2_SIMPLE_CURRENT_MOMENTS"
                            << " COUNTEREXAMPLE level=" << level
                            << " label=" << label
                            << " paired=" << result.paired
                            << " delta_five=" << result.delta_five
                            << " reverse_core=" << result.reverse_core
                            << " half_packet=" << result.packet
                            << " reverse_half_packet="
                            << result.reverse_packet
                            << " eleven_half_packet="
                            << result.eleven_packet
                            << " eleven_paired="
                            << result.eleven_paired
                            << " eleven_core=" << result.eleven_core
                            << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
            std::cout
                << "SU2_SIMPLE_CURRENT_MOMENTS"
                << " rows=" << rows
                << " maximum_level=" << maximum_level
                << " minimum_paired=" << minimum
                << " witness_level=" << witness_level
                << " witness_label=" << witness_label
                << " minimum_delta_five=" << minimum_delta_five
                << " delta_five_witness_level="
                << delta_five_witness_level
                << " delta_five_witness_label="
                << delta_five_witness_label
                << " minimum_reverse_half_packet=" << minimum_reverse
                << " reverse_witness_level=" << reverse_witness_level
                << " reverse_witness_label=" << reverse_witness_label
                << " minimum_reverse_core=" << minimum_reverse_core
                << " reverse_core_witness_level="
                << reverse_core_witness_level
                << " reverse_core_witness_label="
                << reverse_core_witness_label
                << " minimum_eleven_delta7="
                << minimum_eleven_delta_seven
                << " eleven_delta7_witness=("
                << eleven_seven_witness_level << ','
                << eleven_seven_witness_label << ')'
                << " minimum_eleven_delta5="
                << minimum_eleven_delta_five
                << " eleven_delta5_witness=("
                << eleven_five_witness_level << ','
                << eleven_five_witness_label << ')'
                << " minimum_eleven_delta3="
                << minimum_eleven_delta_three
                << " eleven_delta3_witness=("
                << eleven_three_witness_level << ','
                << eleven_three_witness_label << ')'
                << " minimum_eleven_half_packet="
                << minimum_eleven_packet
                << " eleven_packet_witness=("
                << eleven_packet_witness_level << ','
                << eleven_packet_witness_label << ')'
                << " minimum_eleven_paired="
                << minimum_eleven_paired
                << " eleven_paired_witness=("
                << eleven_paired_witness_level << ','
                << eleven_paired_witness_label << ')'
                << " minimum_eleven_core="
                << minimum_eleven_core
                << " eleven_core_witness=("
                << eleven_core_witness_level << ','
                << eleven_core_witness_label << ')'
                << " counterexamples=0 result=PASS_DISCOVERY\n";
            return EXIT_SUCCESS;
        }
        if (argc != 4) {
            throw std::runtime_error(
                "usage: analyze_su2_simple_current_moments "
                "LEVEL EVEN_LABEL MAXIMUM_POWER | --scan MAXIMUM_LEVEL"
            );
        }
        const int level = parse_positive(argv[1], "level");
        const int label = parse_positive(argv[2], "label");
        const int maximum = parse_positive(argv[3], "maximum power");
        if ((level & 1) != 0 || (label & 1) != 0
            || 2 * label >= level || maximum < 2) {
            throw std::runtime_error(
                "require even level/label, 2*label<level, power>=2"
            );
        }

        std::vector<Integer> state(
            static_cast<std::size_t>(level + 1)
        );
        std::vector<Integer> closed(
            static_cast<std::size_t>(maximum + 1)
        );
        std::vector<Integer> wall(
            static_cast<std::size_t>(maximum + 1)
        );
        state[0] = 1;
        for (int power = 0; power <= maximum; ++power) {
            closed[static_cast<std::size_t>(power)] = state[0];
            wall[static_cast<std::size_t>(power)] =
                state[static_cast<std::size_t>(level)];
            std::cout
                << "power=" << power
                << " f=" << closed[static_cast<std::size_t>(power)]
                << " g=" << wall[static_cast<std::size_t>(power)]
                << '\n';
            if (power != maximum) {
                multiply(level, label, state);
            }
        }
        if (maximum >= 10) {
            const Integer delta_three =
                wall[4] * closed[6] - wall[3] * closed[7];
            const Integer delta_five =
                wall[6] * closed[4] - wall[5] * closed[5];
            const Integer delta_seven = wall[8] - wall[7];
            const Integer packet =
                84 * delta_three + 126 * delta_five
                + 36 * delta_seven + wall[10];
            const Integer paired =
                2 * delta_three + 3 * delta_five;
            const Integer reverse_packet =
                84 * wall[4] * closed[6]
                - 36 * wall[3] * closed[7]
                + 126 * delta_five
                + 36 * wall[8] - 84 * wall[7] + wall[10];
            const Integer reverse_core =
                3 * delta_five
                + 2 * (wall[4] * closed[6] - wall[7]);
            const bool pass =
                paired >= 0 && packet >= 0 && reverse_packet >= 0;
            std::cout
                << "B9 level=" << level
                << " label=" << label
                << " delta3=" << delta_three
                << " delta5=" << delta_five
                << " delta7=" << delta_seven
                << " paired=" << paired
                << " reverse_core=" << reverse_core
                << " half_packet=" << packet
                << " reverse_half_packet=" << reverse_packet
                << " result="
                << (pass ? "PASS_CANDIDATE" : "FAIL_CANDIDATE")
                << '\n';
            if (!pass) {
                return EXIT_FAILURE;
            }
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SIMPLE_CURRENT_MOMENTS FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
