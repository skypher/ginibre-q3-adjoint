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

Integer choose_or_zero(
    const std::vector<std::vector<Integer>>& choose,
    int n,
    int r
) {
    if (r < 0 || r > n) {
        return 0;
    }
    return choose[static_cast<std::size_t>(n)]
                 [static_cast<std::size_t>(r)];
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_terminal_first_passage "
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
        std::size_t active_first_passage_rows = 0;
        std::size_t negative_kernel_rows = 0;
        std::size_t kernel_sign_pattern_rows = 0;
        std::size_t reverse_kernel_sign_changes = 0;
        std::size_t first_passage_suffix_rows = 0;
        std::size_t negative_first_passage_suffixes = 0;
        std::size_t depth_one_pair_rows = 0;
        std::size_t negative_depth_one_pairs = 0;
        bool depth_one_witness_present = false;
        int depth_one_witness_level = 0;
        int depth_one_witness_label = 0;
        int depth_one_witness_length = 0;
        Integer depth_one_witness_value = 0;
        std::size_t current_rows = 0;
        bool witness_present = false;
        int witness_level = 0;
        int witness_label = 0;
        int witness_length = 0;
        int witness_first = 0;
        Integer witness_hitting = 0;
        Integer witness_kernel = 0;
        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                const auto paths =
                    path_table(level, label, maximum_length);
                const auto f = [&](int power) -> const Integer& {
                    return paths[static_cast<std::size_t>(power)][0];
                };
                const auto g = [&](int power) -> const Integer& {
                    return paths[static_cast<std::size_t>(power)]
                                [static_cast<std::size_t>(level)];
                };
                std::vector<Integer> first_passage(
                    static_cast<std::size_t>(maximum_length + 1)
                );
                for (int power = 0;
                     power <= maximum_length; ++power) {
                    Integer value = g(power);
                    for (int first = 0; first < power; ++first) {
                        value -= first_passage[
                                     static_cast<std::size_t>(first)
                                 ] * f(power - first);
                    }
                    if (value < 0) {
                        throw std::runtime_error(
                            "negative first-passage coefficient"
                        );
                    }
                    first_passage[static_cast<std::size_t>(power)] =
                        value;
                }
                for (int power = 0;
                     power <= maximum_length; ++power) {
                    Integer reconstructed = 0;
                    for (int first = 0; first <= power; ++first) {
                        reconstructed += first_passage[
                                             static_cast<std::size_t>(
                                                 first
                                             )
                                         ] * f(power - first);
                    }
                    if (reconstructed != g(power)) {
                        throw std::runtime_error(
                            "first-passage convolution failed"
                        );
                    }
                }
                ++parameter_rows;

                for (int length = 4;
                     length <= maximum_length; length += 2) {
                    Integer decomposed_current = 0;
                    std::vector<Integer> contributions(
                        static_cast<std::size_t>(length + 1)
                    );
                    std::vector<Integer> kernels(
                        static_cast<std::size_t>(length + 1)
                    );
                    for (int first = 0;
                         first <= length; ++first) {
                        const Integer& hitting =
                            first_passage[
                                static_cast<std::size_t>(first)
                            ];
                        if (hitting == 0) {
                            continue;
                        }
                        Integer kernel = 0;
                        for (int left = 0;
                             left <= length - first; ++left) {
                            const Integer coefficient =
                                choose_or_zero(
                                    choose,
                                    length - 1,
                                    left + first - 1
                                )
                                - choose_or_zero(
                                    choose,
                                    length - 1,
                                    left + first
                                );
                            kernel += coefficient
                                * f(left)
                                * f(length - first - left);
                        }
                        ++active_first_passage_rows;
                        if (kernel < 0) {
                            ++negative_kernel_rows;
                            if (!witness_present) {
                                witness_present = true;
                                witness_level = level;
                                witness_label = label;
                                witness_length = length;
                                witness_first = first;
                                witness_hitting = hitting;
                                witness_kernel = kernel;
                            }
                        }
                        contributions[
                            static_cast<std::size_t>(first)
                        ] = hitting * kernel;
                        kernels[static_cast<std::size_t>(first)] =
                            kernel;
                        decomposed_current += contributions[
                            static_cast<std::size_t>(first)
                        ];
                    }
                    Integer suffix = 0;
                    for (int first = length;
                         first >= 0; --first) {
                        if (
                            first_passage[
                                static_cast<std::size_t>(first)
                            ] == 0
                        ) {
                            continue;
                        }
                        suffix += contributions[
                            static_cast<std::size_t>(first)
                        ];
                        ++first_passage_suffix_rows;
                        if (suffix < 0) {
                            ++negative_first_passage_suffixes;
                            std::cout
                                << "SU2_TERMINAL_FIRST_PASSAGE"
                                << " suffix_result=FAIL"
                                << " level=" << level
                                << " label=" << label
                                << " length=" << length
                                << " cutoff_first_hit=" << first
                                << " suffix=" << suffix
                                << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                    bool positive_kernel_seen = false;
                    for (int first = 0;
                         first <= length; ++first) {
                        if (
                            first_passage[
                                static_cast<std::size_t>(first)
                            ] == 0
                        ) {
                            continue;
                        }
                        const Integer& kernel =
                            kernels[static_cast<std::size_t>(first)];
                        if (kernel > 0) {
                            positive_kernel_seen = true;
                        } else if (
                            kernel < 0 && positive_kernel_seen
                        ) {
                            ++reverse_kernel_sign_changes;
                            std::cout
                                << "SU2_TERMINAL_FIRST_PASSAGE"
                                << " sign_pattern_result=FAIL"
                                << " level=" << level
                                << " label=" << label
                                << " length=" << length
                                << " first_hit=" << first
                                << " kernel=" << kernel
                                << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                    ++kernel_sign_pattern_rows;
                    const int middle = length / 2;
                    if (
                        first_passage[
                            static_cast<std::size_t>(middle - 1)
                        ] != 0
                    ) {
                        ++depth_one_pair_rows;
                        const Integer pair =
                            contributions[
                                static_cast<std::size_t>(middle - 1)
                            ]
                            + contributions[
                                static_cast<std::size_t>(middle)
                            ];
                        if (pair < 0) {
                            ++negative_depth_one_pairs;
                            if (!depth_one_witness_present) {
                                depth_one_witness_present = true;
                                depth_one_witness_level = level;
                                depth_one_witness_label = label;
                                depth_one_witness_length = length;
                                depth_one_witness_value = pair;
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
                    if (decomposed_current != direct_current) {
                        throw std::runtime_error(
                            "first-passage current identity failed"
                        );
                    }
                    if (direct_current < 0) {
                        throw std::runtime_error(
                            "negative terminal current"
                        );
                    }
                    ++current_rows;
                }
            }
        }

        std::cout
            << "SU2_TERMINAL_FIRST_PASSAGE"
            << " maximum_level=" << maximum_level
            << " maximum_even_length=" << maximum_length
            << " parameter_rows=" << parameter_rows
            << " active_first_passage_rows="
            << active_first_passage_rows
            << " negative_kernel_rows=" << negative_kernel_rows
            << " kernel_sign_pattern_rows="
            << kernel_sign_pattern_rows
            << " reverse_kernel_sign_changes="
            << reverse_kernel_sign_changes
            << " first_passage_suffix_rows="
            << first_passage_suffix_rows
            << " negative_first_passage_suffixes="
            << negative_first_passage_suffixes
            << " depth_one_pair_rows=" << depth_one_pair_rows
            << " negative_depth_one_pairs="
            << negative_depth_one_pairs
            << " current_rows=" << current_rows
            << " convolution_identity=PASS"
            << " current_identity=PASS"
            << " kernel_cone="
            << (negative_kernel_rows == 0 ? "PASS_DISCOVERY" : "FAIL");
        if (witness_present) {
            std::cout
                << " first_negative_kernel={level=" << witness_level
                << ",label=" << witness_label
                << ",length=" << witness_length
                << ",first_hit=" << witness_first
                << ",hitting_weight=" << witness_hitting
                << ",kernel=" << witness_kernel
                << '}';
        }
        if (depth_one_witness_present) {
            std::cout
                << " first_negative_depth_one_pair={level="
                << depth_one_witness_level
                << ",label=" << depth_one_witness_label
                << ",length=" << depth_one_witness_length
                << ",pair=" << depth_one_witness_value
                << '}';
        }
        std::cout
            << " suffix_result=PASS_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_FIRST_PASSAGE FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
