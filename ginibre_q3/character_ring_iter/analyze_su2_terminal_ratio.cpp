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

std::vector<std::vector<Integer>> path_table(
    int level,
    int label,
    int maximum_power
) {
    std::vector<std::vector<Integer>> table(
        static_cast<std::size_t>(maximum_power + 1),
        std::vector<Integer>(static_cast<std::size_t>(level + 1))
    );
    table[0][0] = 1;
    for (int power = 1; power <= maximum_power; ++power) {
        for (int source = 0; source <= level; ++source) {
            const Integer& coefficient =
                table[static_cast<std::size_t>(power - 1)]
                     [static_cast<std::size_t>(source)];
            if (coefficient == 0) {
                continue;
            }
            for (int output = 0; output <= level; ++output) {
                if (fuses(level, label, source, output)) {
                    table[static_cast<std::size_t>(power)]
                         [static_cast<std::size_t>(output)]
                        += coefficient;
                }
            }
        }
    }
    return table;
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
        choose[static_cast<std::size_t>(n)]
              [static_cast<std::size_t>(n)] = 1;
        for (int r = 1; r < n; ++r) {
            choose[static_cast<std::size_t>(n)]
                  [static_cast<std::size_t>(r)] =
                choose[static_cast<std::size_t>(n - 1)]
                      [static_cast<std::size_t>(r - 1)]
                + choose[static_cast<std::size_t>(n - 1)]
                        [static_cast<std::size_t>(r)];
        }
    }
    return choose;
}

struct Witness {
    bool present = false;
    Integer determinant = 0;
    int level = 0;
    int label = 0;
    int length = 0;
    int first = 0;
    int second = 0;
    Integer g_first = 0;
    Integer f_first = 0;
    Integer g_second = 0;
    Integer f_second = 0;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_terminal_ratio "
                "MAXIMUM_LEVEL MAXIMUM_EVEN_LENGTH"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_length =
            parse_positive(argv[2], "maximum even length");
        if (maximum_level < 3 || maximum_length < 4
            || (maximum_length & 1) != 0) {
            throw std::runtime_error(
                "require maximum level>=3 and even length>=4"
            );
        }

        const auto choose = binomial_table(maximum_length);
        std::size_t parameter_rows = 0;
        std::size_t determinant_rows = 0;
        std::size_t negative_determinants = 0;
        std::size_t current_rows = 0;
        std::size_t negative_currents = 0;
        std::size_t outer_suffix_rows = 0;
        std::size_t negative_outer_suffixes = 0;
        Witness first_negative;
        Witness first_negative_outer_suffix;

        for (int level = 3; level <= maximum_level; ++level) {
            for (int label = 1; 2 * label < level; ++label) {
                const auto paths =
                    path_table(level, label, maximum_length);
                const auto f = [&](int power) -> const Integer& {
                    return paths[static_cast<std::size_t>(power)][0];
                };
                const auto g = [&](int power) -> const Integer& {
                    return paths[static_cast<std::size_t>(power)]
                                [static_cast<std::size_t>(level)];
                };
                ++parameter_rows;

                for (int length = 4;
                     length <= maximum_length; length += 2) {
                    Integer paired_current = 0;
                    std::vector<Integer> pair_contributions(
                        static_cast<std::size_t>(length + 1)
                    );
                    for (int first = length / 2 + 1;
                         first <= length; ++first) {
                        const int second = length - first;
                        const Integer determinant =
                            g(first) * f(second)
                            - g(second) * f(first);
                        const Integer numerator =
                            (2 * first - length)
                            * choose[static_cast<std::size_t>(length)]
                                    [static_cast<std::size_t>(first)];
                        if (numerator % length != 0) {
                            throw std::runtime_error(
                                "nonintegral paired coefficient"
                            );
                        }
                        const Integer contribution =
                            (numerator / length) * determinant;
                        pair_contributions[
                            static_cast<std::size_t>(first)
                        ] = contribution;
                        paired_current += contribution;
                        ++determinant_rows;
                        if (determinant < 0) {
                            ++negative_determinants;
                            if (!first_negative.present) {
                                first_negative = {
                                    true,
                                    determinant,
                                    level,
                                    label,
                                    length,
                                    first,
                                    second,
                                    g(first),
                                    f(first),
                                    g(second),
                                    f(second)
                                };
                            }
                        }
                    }
                    Integer outer_suffix = 0;
                    for (int first = length;
                         first > length / 2; --first) {
                        outer_suffix += pair_contributions[
                            static_cast<std::size_t>(first)
                        ];
                        ++outer_suffix_rows;
                        if (outer_suffix < 0) {
                            ++negative_outer_suffixes;
                            if (!first_negative_outer_suffix.present) {
                                const int second = length - first;
                                first_negative_outer_suffix = {
                                    true,
                                    outer_suffix,
                                    level,
                                    label,
                                    length,
                                    first,
                                    second,
                                    g(first),
                                    f(first),
                                    g(second),
                                    f(second)
                                };
                            }
                        }
                    }

                    Integer direct_current = 0;
                    for (int exponent = 0;
                         exponent < length; ++exponent) {
                        direct_current +=
                            choose[static_cast<std::size_t>(length - 1)]
                                  [static_cast<std::size_t>(exponent)]
                            * (
                                g(exponent + 1)
                                    * f(length - 1 - exponent)
                                - g(exponent)
                                    * f(length - exponent)
                            );
                    }
                    if (direct_current != paired_current) {
                        throw std::runtime_error(
                            "paired-current identity failed"
                        );
                    }
                    ++current_rows;
                    if (direct_current < 0) {
                        ++negative_currents;
                    }
                }
            }
        }

        std::cout
            << "SU2_TERMINAL_RATIO"
            << " maximum_level=" << maximum_level
            << " maximum_even_length=" << maximum_length
            << " parameter_rows=" << parameter_rows
            << " determinant_rows=" << determinant_rows
            << " negative_determinants=" << negative_determinants
            << " current_rows=" << current_rows
            << " negative_currents=" << negative_currents
            << " outer_suffix_rows=" << outer_suffix_rows
            << " negative_outer_suffixes=" << negative_outer_suffixes
            << " paired_identity=PASS";
        if (first_negative.present) {
            std::cout
                << " first_negative={level=" << first_negative.level
                << ",label=" << first_negative.label
                << ",length=" << first_negative.length
                << ",A=" << first_negative.first
                << ",B=" << first_negative.second
                << ",det=" << first_negative.determinant
                << ",g_A=" << first_negative.g_first
                << ",f_A=" << first_negative.f_first
                << ",g_B=" << first_negative.g_second
                << ",f_B=" << first_negative.f_second
                << '}';
        }
        if (first_negative_outer_suffix.present) {
            std::cout
                << " first_negative_outer_suffix={level="
                << first_negative_outer_suffix.level
                << ",label=" << first_negative_outer_suffix.label
                << ",length=" << first_negative_outer_suffix.length
                << ",cutoff_A=" << first_negative_outer_suffix.first
                << ",sum=" << first_negative_outer_suffix.determinant
                << '}';
        }
        std::cout
            << " determinant_cone="
            << (negative_determinants == 0 ? "PASS_DISCOVERY" : "FAIL")
            << " current_sign="
            << (negative_currents == 0 ? "PASS_DISCOVERY" : "FAIL")
            << " outer_suffix_cone="
            << (
                negative_outer_suffixes == 0
                    ? "PASS_DISCOVERY" : "FAIL"
            )
            << '\n';
        return negative_currents == 0
            ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_RATIO FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
