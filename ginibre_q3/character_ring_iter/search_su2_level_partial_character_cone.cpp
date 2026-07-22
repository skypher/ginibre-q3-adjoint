#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

template <class Function>
void for_each_output(int level, int left, int right, Function function) {
    const int lower = std::abs(left - right);
    const int upper = std::min(left + right, 2 * level - left - right);
    for (int output = lower; output <= upper; output += 2) {
        function(output);
    }
}

std::vector<std::vector<int>> multisets(int level, int maximum_size) {
    std::vector<std::vector<int>> result;
    std::vector<int> current;
    const auto visit = [&](const auto& self, int first) -> void {
        if (!current.empty()) {
            result.push_back(current);
        }
        if (current.size() == static_cast<std::size_t>(maximum_size)) {
            return;
        }
        for (int label = first; label <= level; ++label) {
            current.push_back(label);
            self(self, label);
            current.pop_back();
        }
    };
    visit(visit, 1);
    return result;
}

std::vector<cpp_int> signed_product(
    int level,
    const std::vector<int>& labels,
    std::uint64_t minus_mask
) {
    const int width = level + 1;
    std::vector<cpp_int> current(
        static_cast<std::size_t>(width * width)
    );
    current[0] = 1;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        const int label = labels[index];
        const bool minus = ((minus_mask >> index) & 1U) != 0U;
        std::vector<cpp_int> next(current.size());
        for (int left = 0; left <= level; ++left) {
            for (int right = 0; right <= level; ++right) {
                const cpp_int& value = current[static_cast<std::size_t>(
                    left * width + right
                )];
                if (value == 0) {
                    continue;
                }
                for_each_output(level, left, label, [&](int output) {
                    next[static_cast<std::size_t>(
                        output * width + right
                    )] += value;
                });
                for_each_output(level, right, label, [&](int output) {
                    next[static_cast<std::size_t>(
                        left * width + output
                    )] += minus ? -value : value;
                });
            }
        }
        current = std::move(next);
    }
    return current;
}

std::vector<cpp_int> add_plus_factor(
    int level,
    const std::vector<cpp_int>& current,
    int label
) {
    const int width = level + 1;
    if (current.size() != static_cast<std::size_t>(width * width)) {
        throw std::runtime_error("coefficient matrix has wrong size");
    }
    std::vector<cpp_int> next(current.size(), 0);
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            const cpp_int& value = current[static_cast<std::size_t>(
                left * width + right
            )];
            if (value == 0) {
                continue;
            }
            for_each_output(level, left, label, [&](int output) {
                next[static_cast<std::size_t>(
                    output * width + right
                )] += value;
            });
            for_each_output(level, right, label, [&](int output) {
                next[static_cast<std::size_t>(
                    left * width + output
                )] += value;
            });
        }
    }
    return next;
}

std::vector<cpp_int> fusion_hankel_defect(
    int level,
    const std::vector<cpp_int>& coefficient
) {
    const int width = level + 1;
    if (coefficient.size()
        != static_cast<std::size_t>(width * width)) {
        throw std::runtime_error("coefficient matrix has wrong size");
    }
    std::vector<cpp_int> defect(coefficient.size(), 0);
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            cpp_int hankel = 0;
            for_each_output(level, left, right, [&](int channel) {
                hankel += coefficient[static_cast<std::size_t>(
                    channel * width
                )];
            });
            const std::size_t entry = static_cast<std::size_t>(
                left * width + right
            );
            defect[entry] = hankel - coefficient[entry];
        }
    }
    return defect;
}

std::vector<cpp_int> subset_multiplicities(
    int level,
    const std::vector<int>& labels,
    std::uint64_t mask
) {
    std::vector<cpp_int> current(static_cast<std::size_t>(level + 1), 0);
    current[0] = 1;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        if (((mask >> index) & 1U) == 0U) {
            continue;
        }
        std::vector<cpp_int> next(current.size(), 0);
        for (int input = 0; input <= level; ++input) {
            const cpp_int& value = current[static_cast<std::size_t>(input)];
            if (value == 0) {
                continue;
            }
            for_each_output(level, input, labels[index], [&](int output) {
                next[static_cast<std::size_t>(output)] += value;
            });
        }
        current = std::move(next);
    }
    return current;
}

bool inspect_word(
    int level,
    const std::vector<int>& labels,
    std::uint64_t minus_mask,
    int& failed_target,
    cpp_int& failed_value
) {
    const int width = level + 1;
    const std::vector<cpp_int> current
        = signed_product(level, labels, minus_mask);
    for (int target = 0; target <= level; ++target) {
        const cpp_int& value = current[static_cast<std::size_t>(
            target * width
        )];
        if (value < 0) {
            failed_target = target;
            failed_value = value;
            return false;
        }
    }
    return true;
}

int popcount(std::uint64_t value) {
    int answer = 0;
    while (value != 0U) {
        answer += static_cast<int>(value & 1U);
        value >>= 1U;
    }
    return answer;
}

void print_word(const std::vector<int>& labels, std::uint64_t minus_mask) {
    std::cout << '[';
    for (std::size_t index = 0; index < labels.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << (((minus_mask >> index) & 1U) != 0U ? '-' : '+')
                  << labels[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4
            && std::string(argv[1]) == "simple-current-exterior") {
            const int maximum_level = std::stoi(argv[2]);
            const int maximum_factors = std::stoi(argv[3]);
            if (maximum_level < 2 || maximum_factors < 0
                || maximum_factors >= 63) {
                throw std::runtime_error(
                    "invalid simple-current exterior bounds"
                );
            }
            std::uint64_t words_tested = 0U;
            std::uint64_t boundary_entries_tested = 0U;
            std::uint64_t interior_entries_tested = 0U;
            std::uint64_t negative_interior_entries = 0U;
            std::uint64_t odd_orbit_entries_tested = 0U;
            bool printed_negative_interior = false;
            for (int level = 2; level <= maximum_level; ++level) {
                const int width = level + 1;
                auto words = multisets(level - 1, maximum_factors);
                words.insert(words.begin(), std::vector<int>{});
                for (const auto& labels : words) {
                    std::vector<cpp_int> exterior(
                        static_cast<std::size_t>(width * width), 0
                    );
                    exterior[static_cast<std::size_t>(level * width)] = 1;
                    const auto add_wedge = [width](
                        std::vector<cpp_int>& destination,
                        int first,
                        int second,
                        const cpp_int& value
                    ) {
                        if (first > second) {
                            destination[static_cast<std::size_t>(
                                first * width + second
                            )] += value;
                        } else if (second > first) {
                            destination[static_cast<std::size_t>(
                                second * width + first
                            )] -= value;
                        }
                    };
                    for (int label : labels) {
                        std::vector<cpp_int> next(exterior.size(), 0);
                        for (int first = 1; first <= level; ++first) {
                            for (int second = 0; second < first; ++second) {
                                const cpp_int& value = exterior[
                                    static_cast<std::size_t>(
                                        first * width + second
                                    )
                                ];
                                if (value == 0) {
                                    continue;
                                }
                                for_each_output(
                                    level, first, label, [&](int output) {
                                        add_wedge(
                                            next, output, second, value
                                        );
                                    }
                                );
                                for_each_output(
                                    level, second, label, [&](int output) {
                                        add_wedge(next, first, output, value);
                                    }
                                );
                            }
                        }
                        exterior = std::move(next);
                    }
                    const std::vector<cpp_int> defect
                        = fusion_hankel_defect(
                            level, signed_product(level, labels, 0U)
                        );
                    for (int target = 1; target <= level; ++target) {
                        const cpp_int& exterior_value = exterior[
                            static_cast<std::size_t>(target * width)
                        ];
                        const cpp_int& boundary = defect[
                            static_cast<std::size_t>(
                                target * width + level
                            )
                        ];
                        ++boundary_entries_tested;
                        if (exterior_value != boundary || boundary < 0) {
                            std::cout
                                << "SU2_SIMPLE_CURRENT_EXTERIOR result=FAIL"
                                << " level=" << level << " word=";
                            print_word(labels, 0U);
                            std::cout << " target=" << target
                                      << " exterior=" << exterior_value
                                      << " boundary=" << boundary << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                    if ((level & 1) != 0) {
                        const int orbit_width = (level + 1) / 2;
                        std::vector<cpp_int> orbit(
                            static_cast<std::size_t>(
                                orbit_width * orbit_width
                            ), 0
                        );
                        orbit[0] = 1;
                        int odd_label_parity = 0;
                        for (int label : labels) {
                            odd_label_parity ^= label & 1;
                            std::vector<cpp_int> next(orbit.size(), 0);
                            for (int left = 0; left < orbit_width; ++left) {
                                for (int right = 0;
                                     right < orbit_width; ++right) {
                                    const cpp_int& value = orbit[
                                        static_cast<std::size_t>(
                                            left * orbit_width + right
                                        )
                                    ];
                                    if (value == 0) {
                                        continue;
                                    }
                                    for_each_output(
                                        level, left, label, [&](int output) {
                                            const int folded = std::min(
                                                output, level - output
                                            );
                                            next[static_cast<std::size_t>(
                                                folded * orbit_width + right
                                            )] += value;
                                        }
                                    );
                                    for_each_output(
                                        level, right, label, [&](int output) {
                                            const int folded = std::min(
                                                output, level - output
                                            );
                                            cpp_int& destination = next[
                                                static_cast<std::size_t>(
                                                    left * orbit_width + folded
                                                )
                                            ];
                                            if ((label & 1) == 0) {
                                                destination += value;
                                            } else {
                                                destination -= value;
                                            }
                                        }
                                    );
                                }
                            }
                            orbit = std::move(next);
                        }
                        for (int folded_target = 0;
                             folded_target < orbit_width;
                             ++folded_target) {
                            const cpp_int& orbit_value = orbit[
                                static_cast<std::size_t>(
                                    folded_target * orbit_width
                                )
                            ];
                            const int lower_target = folded_target;
                            const int upper_target = level - folded_target;
                            const bool use_lower
                                = odd_label_parity
                                    == ((level - folded_target) & 1);
                            const int physical_target
                                = use_lower ? lower_target : upper_target;
                            const int zero_target
                                = use_lower ? upper_target : lower_target;
                            const cpp_int& exterior_value = exterior[
                                static_cast<std::size_t>(
                                    physical_target * width
                                )
                            ];
                            const cpp_int& zero_value = exterior[
                                static_cast<std::size_t>(
                                    zero_target * width
                                )
                            ];
                            ++odd_orbit_entries_tested;
                            if (orbit_value != exterior_value
                                || zero_value != 0) {
                                std::cout
                                    << "SU2_SIMPLE_CURRENT_EXTERIOR "
                                    << "result=FAIL_ODD_ORBIT level="
                                    << level << " word=";
                                print_word(labels, 0U);
                                std::cout
                                    << " folded_target=" << folded_target
                                    << " physical_target=" << physical_target
                                    << " orbit=" << orbit_value
                                    << " exterior=" << exterior_value
                                    << " paired_zero=" << zero_value << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                    }
                    for (int first = 1; first <= level; ++first) {
                        for (int second = 0; second < first; ++second) {
                            const cpp_int& value = exterior[
                                static_cast<std::size_t>(
                                    first * width + second
                                )
                            ];
                            const int reflected_first = level - second;
                            const int reflected_second = level - first;
                            const cpp_int& reflected = exterior[
                                static_cast<std::size_t>(
                                    reflected_first * width
                                    + reflected_second
                                )
                            ];
                            ++interior_entries_tested;
                            if (value != reflected) {
                                std::cout
                                    << "SU2_SIMPLE_CURRENT_EXTERIOR "
                                    << "result=FAIL_SECTOR level=" << level
                                    << " word=";
                                print_word(labels, 0U);
                                std::cout << " state=(" << first << ','
                                          << second << ") value=" << value
                                          << " reflected=" << reflected
                                          << '\n';
                                return EXIT_FAILURE;
                            }
                            if (value < 0) {
                                ++negative_interior_entries;
                                if (!printed_negative_interior) {
                                    std::cout
                                        << "negative_interior level="
                                        << level << " word=";
                                    print_word(labels, 0U);
                                    std::cout << " state=(" << first << ','
                                              << second << ") value="
                                              << value << '\n';
                                    printed_negative_interior = true;
                                }
                            }
                        }
                    }
                    ++words_tested;
                }
                std::cout << "progress level=" << level
                          << " words=" << words_tested
                          << " boundary_entries="
                          << boundary_entries_tested
                          << " interior_entries=" << interior_entries_tested
                          << " negative_interior_entries="
                          << negative_interior_entries
                          << " odd_orbit_entries="
                          << odd_orbit_entries_tested << '\n';
            }
            std::cout << "SU2_SIMPLE_CURRENT_EXTERIOR levels=2.."
                      << maximum_level
                      << " maximum_factors=" << maximum_factors
                      << " words=" << words_tested
                      << " boundary_entries=" << boundary_entries_tested
                      << " interior_entries=" << interior_entries_tested
                      << " negative_interior_entries="
                      << negative_interior_entries
                      << " odd_orbit_entries="
                      << odd_orbit_entries_tested
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc == 4
            && std::string(argv[1]) == "finite-c2-wall-packets") {
            const int maximum_level = std::stoi(argv[2]);
            const int maximum_factors = std::stoi(argv[3]);
            if (maximum_level < 2 || maximum_factors < 0
                || maximum_factors >= 63) {
                throw std::runtime_error(
                    "invalid finite C2 wall-packet bounds"
                );
            }
            std::uint64_t words_tested = 0U;
            std::uint64_t coefficients_tested = 0U;
            std::uint64_t negative_coefficients = 0U;
            std::uint64_t prefixes_tested = 0U;
            std::uint64_t diagonal_prefix_failures = 0U;
            std::uint64_t two_layer_packets_tested = 0U;
            bool printed_diagonal_prefix_failure = false;
            cpp_int minimum_prefix = 0;
            cpp_int minimum_two_layer_packet = 0;
            bool have_prefix = false;
            bool have_two_layer_packet = false;
            for (int level = 2; level <= maximum_level; ++level) {
                const int width = level + 1;
                auto words = multisets(level - 1, maximum_factors);
                words.insert(words.begin(), std::vector<int>{});
                for (const auto& labels : words) {
                    std::vector<cpp_int> product(
                        static_cast<std::size_t>(width * width), 0
                    );
                    for (int first = 0; first < level; ++first) {
                        const int second = level - 1 - first;
                        product[static_cast<std::size_t>(
                            first * width + second
                        )] = 1;
                    }
                    for (int label : labels) {
                        product = add_plus_factor(level, product, label);
                    }
                    const auto at = [&](int first, int second) -> cpp_int {
                        if (first < 0 || first > level
                            || second < 0 || second > level) {
                            return 0;
                        }
                        return product[static_cast<std::size_t>(
                            first * width + second
                        )];
                    };
                    std::vector<cpp_int> packet(
                        static_cast<std::size_t>(level * level), 0
                    );
                    const auto packet_index = [&](int first, int second) {
                        return static_cast<std::size_t>(
                            first * level + second
                        );
                    };
                    for (int first = 0; first < level; ++first) {
                        for (int second = 0; second <= first; ++second) {
                            cpp_int& coefficient = packet[
                                packet_index(first, second)
                            ];
                            coefficient = at(first, second)
                                + at(first + 2, second)
                                - at(first + 1, second - 1)
                                - at(first + 1, second + 1);
                            ++coefficients_tested;
                            if (coefficient < 0) {
                                ++negative_coefficients;
                            }
                        }
                    }
                    const std::vector<cpp_int> defect
                        = fusion_hankel_defect(
                            level, signed_product(level, labels, 0U)
                        );
                    for (int first = 0; first < level; ++first) {
                        const cpp_int& row = packet[
                            packet_index(first, 0)
                        ];
                        const cpp_int& boundary = defect[
                            static_cast<std::size_t>(
                                (first + 1) * width + level
                            )
                        ];
                        if (row != boundary) {
                            std::cout
                                << "SU2_FINITE_C2_WALL_PACKETS "
                                << "result=FAIL_IDENTITY level=" << level
                                << " word=";
                            print_word(labels, 0U);
                            std::cout << " weight=(" << first
                                      << ",0) coefficient=" << row
                                      << " boundary=" << boundary << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                    for (int degree = 0;
                         degree <= 2 * (level - 1); ++degree) {
                        const int first_second
                            = std::max(0, degree - (level - 1));
                        const int last_second = degree / 2;
                        cpp_int prefix = 0;
                        for (int second = first_second;
                             second <= last_second; ++second) {
                            const int first = degree - second;
                            prefix += packet[packet_index(first, second)];
                            ++prefixes_tested;
                            if (!have_prefix || prefix < minimum_prefix) {
                                minimum_prefix = prefix;
                                have_prefix = true;
                            }
                            if (prefix < 0) {
                                ++diagonal_prefix_failures;
                                if (!printed_diagonal_prefix_failure) {
                                    std::cout
                                        << "diagonal_prefix_failure level="
                                        << level << " word=";
                                    print_word(labels, 0U);
                                    std::cout << " degree=" << degree
                                              << " through_second=" << second
                                              << " prefix=" << prefix
                                              << " profile=";
                                    for (int profile_second = first_second;
                                         profile_second <= last_second;
                                         ++profile_second) {
                                        const int profile_first
                                            = degree - profile_second;
                                        std::cout << " (" << profile_first
                                                  << ',' << profile_second
                                                  << ")="
                                                  << packet[packet_index(
                                                         profile_first,
                                                         profile_second
                                                     )];
                                    }
                                    std::cout << " full_profile=";
                                    for (int profile_first = 0;
                                         profile_first < level;
                                         ++profile_first) {
                                        for (int profile_second = 0;
                                             profile_second <= profile_first;
                                             ++profile_second) {
                                            const cpp_int& value = packet[
                                                packet_index(
                                                    profile_first,
                                                    profile_second
                                                )
                                            ];
                                            if (value != 0) {
                                                std::cout
                                                    << " (" << profile_first
                                                    << ',' << profile_second
                                                    << ")=" << value;
                                            }
                                        }
                                    }
                                    std::cout << '\n';
                                    printed_diagonal_prefix_failure = true;
                                }
                            }
                        }
                    }
                    for (int first_target = 1;
                         first_target <= level; ++first_target) {
                        for (int second_target = 1;
                             second_target <= first_target;
                             ++second_target) {
                            if (first_target + second_target - 2
                                >= level) {
                                continue;
                            }
                            cpp_int two_layer_packet = 0;
                            for (int index = 0;
                                 index < second_target; ++index) {
                                two_layer_packet += packet[packet_index(
                                    first_target + second_target - 2 - index,
                                    index
                                )];
                            }
                            for (int index = 1;
                                 index < second_target; ++index) {
                                two_layer_packet += packet[packet_index(
                                    first_target + second_target - 3 - index,
                                    index - 1
                                )];
                            }
                            ++two_layer_packets_tested;
                            if (!have_two_layer_packet
                                || two_layer_packet
                                    < minimum_two_layer_packet) {
                                minimum_two_layer_packet = two_layer_packet;
                                have_two_layer_packet = true;
                            }
                            if (two_layer_packet < 0) {
                                std::cout
                                    << "SU2_FINITE_C2_WALL_PACKETS "
                                    << "result=FAIL_TWO_LAYER level="
                                    << level << " word=";
                                print_word(labels, 0U);
                                std::cout << " target=(" << first_target
                                          << ',' << second_target
                                          << ") packet=" << two_layer_packet
                                          << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                    }
                    ++words_tested;
                }
                std::cout << "progress level=" << level
                          << " words=" << words_tested
                          << " coefficients=" << coefficients_tested
                          << " negative_coefficients="
                          << negative_coefficients
                          << " prefixes=" << prefixes_tested
                          << " diagonal_prefix_failures="
                          << diagonal_prefix_failures
                          << " minimum_prefix="
                          << (have_prefix ? minimum_prefix : cpp_int{0})
                          << " two_layer_packets="
                          << two_layer_packets_tested
                          << " minimum_two_layer_packet="
                          << (have_two_layer_packet
                                  ? minimum_two_layer_packet : cpp_int{0})
                          << '\n';
            }
            std::cout << "SU2_FINITE_C2_WALL_PACKETS levels=2.."
                      << maximum_level
                      << " maximum_factors=" << maximum_factors
                      << " words=" << words_tested
                      << " coefficients=" << coefficients_tested
                      << " negative_coefficients=" << negative_coefficients
                      << " prefixes=" << prefixes_tested
                      << " diagonal_prefix_failures="
                      << diagonal_prefix_failures
                      << " minimum_prefix="
                      << (have_prefix ? minimum_prefix : cpp_int{0})
                      << " two_layer_packets=" << two_layer_packets_tested
                      << " minimum_two_layer_packet="
                      << (have_two_layer_packet
                              ? minimum_two_layer_packet : cpp_int{0})
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc == 4
            && std::string(argv[1]) == "finite-c2-wall-cone") {
            const int maximum_level = std::stoi(argv[2]);
            const int maximum_factors = std::stoi(argv[3]);
            if (maximum_level < 2 || maximum_factors < 0) {
                throw std::runtime_error(
                    "invalid finite C2 wall-cone bounds"
                );
            }
            std::uint64_t words_tested = 0U;
            std::uint64_t coefficients_tested = 0U;
            for (int level = 2; level <= maximum_level; ++level) {
                const int width = level + 1;
                auto words = multisets(level - 1, maximum_factors);
                words.insert(words.begin(), std::vector<int>{});
                for (const auto& labels : words) {
                    std::vector<cpp_int> product(
                        static_cast<std::size_t>(width * width), 0
                    );
                    for (int first = 0; first < level; ++first) {
                        const int second = level - 1 - first;
                        product[static_cast<std::size_t>(
                            first * width + second
                        )] = 1;
                    }
                    for (int label : labels) {
                        product = add_plus_factor(level, product, label);
                    }
                    const auto at = [&](int first, int second) -> cpp_int {
                        if (first < 0 || first > level
                            || second < 0 || second > level) {
                            return 0;
                        }
                        return product[static_cast<std::size_t>(
                            first * width + second
                        )];
                    };
                    const std::vector<cpp_int> defect
                        = fusion_hankel_defect(
                            level, signed_product(level, labels, 0U)
                        );
                    for (int first = 0; first < level; ++first) {
                        for (int second = 0; second <= first; ++second) {
                            const cpp_int coefficient = at(first, second)
                                + at(first + 2, second)
                                - at(first + 1, second - 1)
                                - at(first + 1, second + 1);
                            ++coefficients_tested;
                            if (second == 0) {
                                const cpp_int& boundary = defect[
                                    static_cast<std::size_t>(
                                        (first + 1) * width + level
                                    )
                                ];
                                if (coefficient != boundary) {
                                    std::cout
                                        << "SU2_FINITE_C2_WALL_CONE "
                                        << "result=FAIL_IDENTITY level="
                                        << level << " word=";
                                    print_word(labels, 0U);
                                    std::cout << " weight=(" << first << ','
                                              << second << ") coefficient="
                                              << coefficient
                                              << " boundary=" << boundary
                                              << '\n';
                                    return EXIT_FAILURE;
                                }
                            }
                            if (coefficient < 0) {
                                std::cout
                                    << "SU2_FINITE_C2_WALL_CONE "
                                    << "result=FAIL level=" << level
                                    << " word=";
                                print_word(labels, 0U);
                                std::cout << " weight=(" << first << ','
                                          << second << ") coefficient="
                                          << coefficient << '\n';
                                std::cout << "profile";
                                for (int profile_first = 0;
                                     profile_first < level;
                                     ++profile_first) {
                                    for (int profile_second = 0;
                                         profile_second <= profile_first;
                                         ++profile_second) {
                                        const cpp_int profile_coefficient
                                            = at(profile_first, profile_second)
                                            + at(profile_first + 2,
                                                 profile_second)
                                            - at(profile_first + 1,
                                                 profile_second - 1)
                                            - at(profile_first + 1,
                                                 profile_second + 1);
                                        if (profile_coefficient != 0) {
                                            std::cout << " ("
                                                      << profile_first << ','
                                                      << profile_second << ")="
                                                      << profile_coefficient;
                                        }
                                    }
                                }
                                std::cout << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                    }
                    ++words_tested;
                }
                std::cout << "progress level=" << level
                          << " words=" << words_tested
                          << " coefficients=" << coefficients_tested << '\n';
            }
            std::cout << "SU2_FINITE_C2_WALL_CONE levels=2.."
                      << maximum_level
                      << " maximum_factors=" << maximum_factors
                      << " words=" << words_tested
                      << " coefficients=" << coefficients_tested
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc == 4
            && std::string(argv[1]) == "finite-c2-wall-generator") {
            const int maximum_level = std::stoi(argv[2]);
            const int maximum_label = std::stoi(argv[3]);
            if (maximum_level < 2 || maximum_label < 1) {
                throw std::runtime_error(
                    "invalid finite C2 wall-generator bounds"
                );
            }
            std::uint64_t coefficients_tested = 0U;
            std::uint64_t nonzero_coefficients = 0U;
            for (int level = 2; level <= maximum_level; ++level) {
                const int width = level + 1;
                const int last_label = std::min(level - 1, maximum_label);
                for (int label = 1; label <= last_label; ++label) {
                    std::vector<cpp_int> product(
                        static_cast<std::size_t>(width * width), 0
                    );
                    for (int first = 0; first < level; ++first) {
                        const int second = level - 1 - first;
                        for_each_output(
                            level, first, label, [&](int output) {
                                product[static_cast<std::size_t>(
                                    output * width + second
                                )] += 1;
                            }
                        );
                        for_each_output(
                            level, second, label, [&](int output) {
                                product[static_cast<std::size_t>(
                                    first * width + output
                                )] += 1;
                            }
                        );
                    }
                    const auto at = [&](int first, int second) -> cpp_int {
                        if (first < 0 || first > level
                            || second < 0 || second > level) {
                            return 0;
                        }
                        return product[static_cast<std::size_t>(
                            first * width + second
                        )];
                    };
                    for (int first = 0; first < level; ++first) {
                        for (int second = 0; second <= first; ++second) {
                            const cpp_int coefficient = at(first, second)
                                + at(first + 2, second)
                                - at(first + 1, second - 1)
                                - at(first + 1, second + 1);
                            const cpp_int expected
                                = ((first == level - 1 && second == label)
                                        || (first == level - 1 - label
                                            && second == 0))
                                ? 1 : 0;
                            ++coefficients_tested;
                            if (coefficient != 0) {
                                ++nonzero_coefficients;
                            }
                            if (coefficient != expected) {
                                std::cout
                                    << "SU2_FINITE_C2_WALL_GENERATOR "
                                    << "result=FAIL level=" << level
                                    << " label=" << label
                                    << " weight=(" << first << ','
                                    << second << ") coefficient="
                                    << coefficient
                                    << " expected=" << expected << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                    }
                }
                std::cout << "progress level=" << level
                          << " labels=1.." << last_label
                          << " coefficients=" << coefficients_tested
                          << " nonzero=" << nonzero_coefficients << '\n';
            }
            std::cout << "SU2_FINITE_C2_WALL_GENERATOR levels=2.."
                      << maximum_level
                      << " maximum_label=" << maximum_label
                      << " coefficients=" << coefficients_tested
                      << " nonzero=" << nonzero_coefficients
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc == 4
            && std::string(argv[1]) == "finite-c2-fundamental") {
            const int maximum_level = std::stoi(argv[2]);
            const int maximum_factors = std::stoi(argv[3]);
            if (maximum_level < 1 || maximum_factors < 0) {
                throw std::runtime_error(
                    "invalid finite C2 fundamental bounds"
                );
            }
            std::uint64_t entries_tested = 0U;
            for (int level = 1; level <= maximum_level; ++level) {
                const int width = level + 1;
                std::vector<std::vector<int>> state_index(
                    static_cast<std::size_t>(level),
                    std::vector<int>(static_cast<std::size_t>(level), -1)
                );
                int state_count = 0;
                for (int first = 0; first < level; ++first) {
                    for (int second = 0; second <= first; ++second) {
                        state_index[static_cast<std::size_t>(first)]
                                   [static_cast<std::size_t>(second)]
                            = state_count;
                        ++state_count;
                    }
                }
                std::vector<std::vector<cpp_int>> walks(
                    static_cast<std::size_t>(level),
                    std::vector<cpp_int>(
                        static_cast<std::size_t>(state_count), 0
                    )
                );
                for (int boundary = 1; boundary <= level; ++boundary) {
                    const int state = state_index[
                        static_cast<std::size_t>(boundary - 1)
                    ][0];
                    walks[static_cast<std::size_t>(boundary - 1)]
                         [static_cast<std::size_t>(state)] = 1;
                }
                std::vector<int> labels;
                for (int factors = 0; factors <= maximum_factors; ++factors) {
                    const std::vector<cpp_int> defect
                        = fusion_hankel_defect(
                            level, signed_product(level, labels, 0U)
                        );
                    for (int left = 1; left <= level; ++left) {
                        const int target_state = state_index[
                            static_cast<std::size_t>(left - 1)
                        ][0];
                        for (int right = 1; right <= level; ++right) {
                            const cpp_int& direct = defect[
                                static_cast<std::size_t>(
                                    left * width + right
                                )
                            ];
                            const cpp_int& triangular = walks[
                                static_cast<std::size_t>(right - 1)
                            ][static_cast<std::size_t>(target_state)];
                            ++entries_tested;
                            if (direct != triangular || direct < 0) {
                                std::cout
                                    << "SU2_FINITE_C2_FUNDAMENTAL result=FAIL"
                                    << " level=" << level
                                    << " factors=" << factors
                                    << " left=" << left
                                    << " right=" << right
                                    << " defect=" << direct
                                    << " triangular=" << triangular << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                    }
                    if (factors == maximum_factors) {
                        break;
                    }
                    std::vector<std::vector<cpp_int>> next(
                        walks.size(),
                        std::vector<cpp_int>(
                            static_cast<std::size_t>(state_count), 0
                        )
                    );
                    for (std::size_t boundary = 0;
                         boundary < walks.size(); ++boundary) {
                        for (int first = 0; first < level; ++first) {
                            for (int second = 0; second <= first; ++second) {
                                const int source = state_index[
                                    static_cast<std::size_t>(first)
                                ][static_cast<std::size_t>(second)];
                                const cpp_int& value = walks[boundary][
                                    static_cast<std::size_t>(source)
                                ];
                                if (value == 0) {
                                    continue;
                                }
                                const auto add = [&](int next_first,
                                                     int next_second) {
                                    if (next_first < 0
                                        || next_first >= level
                                        || next_second < 0
                                        || next_second > next_first) {
                                        return;
                                    }
                                    const int destination = state_index[
                                        static_cast<std::size_t>(next_first)
                                    ][static_cast<std::size_t>(next_second)];
                                    next[boundary][static_cast<std::size_t>(
                                        destination
                                    )] += value;
                                };
                                add(first + 1, second);
                                add(first - 1, second);
                                add(first, second + 1);
                                add(first, second - 1);
                            }
                        }
                    }
                    walks = std::move(next);
                    labels.push_back(1);
                }
                std::cout << "progress level=" << level
                          << " factors=0.." << maximum_factors
                          << " entries=" << entries_tested << '\n';
            }
            std::cout << "SU2_FINITE_C2_FUNDAMENTAL levels=1.."
                      << maximum_level
                      << " factors=0.." << maximum_factors
                      << " entries=" << entries_tested
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc == 4
            && std::string(argv[1]) == "simple-current-boundary-local") {
            const int maximum_level = std::stoi(argv[2]);
            const int maximum_factors = std::stoi(argv[3]);
            if (maximum_level < 2 || maximum_factors < 0
                || maximum_factors >= 63) {
                throw std::runtime_error(
                    "invalid local simple-current boundary bounds"
                );
            }
            std::uint64_t words_tested = 0U;
            std::uint64_t allocation_pairs_tested = 0U;
            std::uint64_t entries_tested = 0U;
            std::uint64_t local_boundary_failures = 0U;
            std::uint64_t local_scbi_failures = 0U;
            std::uint64_t shell_boundary_failures = 0U;
            std::uint64_t shell_scbi_failures = 0U;
            bool printed_boundary_failure = false;
            bool printed_scbi_failure = false;
            bool printed_shell_boundary_failure = false;
            bool printed_shell_scbi_failure = false;
            for (int level = 2; level <= maximum_level; ++level) {
                auto words = multisets(level - 1, maximum_factors);
                words.insert(words.begin(), std::vector<int>{});
                for (const auto& labels : words) {
                    ++words_tested;
                    const int width = level + 1;
                    const std::uint64_t allocation_count
                        = std::uint64_t{1} << labels.size();
                    const std::uint64_t full_mask = allocation_count - 1U;
                    const std::size_t shell_count = labels.size() / 2U + 1U;
                    std::vector<std::vector<cpp_int>> shell_coefficients(
                        shell_count,
                        std::vector<cpp_int>(
                            static_cast<std::size_t>(width * width), 0
                        )
                    );
                    for (std::uint64_t mask = 0; mask < allocation_count;
                         ++mask) {
                        const std::uint64_t complement = full_mask ^ mask;
                        const std::vector<cpp_int> left
                            = subset_multiplicities(level, labels, mask);
                        const std::vector<cpp_int> right
                            = subset_multiplicities(
                                level, labels, complement
                            );
                        const std::size_t selected = static_cast<std::size_t>(
                            popcount(mask)
                        );
                        const std::size_t shell = std::min(
                            selected, labels.size() - selected
                        );
                        for (int first = 0; first <= level; ++first) {
                            for (int second = 0; second <= level; ++second) {
                                shell_coefficients[shell][
                                    static_cast<std::size_t>(
                                        first * width + second
                                    )
                                ] += left[static_cast<std::size_t>(first)]
                                    * right[static_cast<std::size_t>(second)];
                            }
                        }
                        if (mask > complement) {
                            continue;
                        }
                        std::vector<cpp_int> coefficient(
                            static_cast<std::size_t>(width * width), 0
                        );
                        for (int first = 0; first <= level; ++first) {
                            for (int second = 0; second <= level; ++second) {
                                coefficient[static_cast<std::size_t>(
                                    first * width + second
                                )] = left[static_cast<std::size_t>(first)]
                                        * right[static_cast<std::size_t>(second)]
                                    + right[static_cast<std::size_t>(first)]
                                        * left[static_cast<std::size_t>(second)];
                            }
                        }
                        const std::vector<cpp_int> defect
                            = fusion_hankel_defect(level, coefficient);
                        ++allocation_pairs_tested;
                        for (int target = 0; target <= level; ++target) {
                            if (std::binary_search(
                                    labels.begin(), labels.end(), target
                                )) {
                                continue;
                            }
                            const cpp_int& boundary_value = defect[
                                static_cast<std::size_t>(
                                    target * width + level
                                )
                            ];
                            if (boundary_value < 0) {
                                ++local_boundary_failures;
                                if (!printed_boundary_failure) {
                                    std::cout
                                        << "local_boundary_failure level="
                                        << level << " prefix=";
                                    print_word(labels, 0U);
                                    std::cout << " allocation=" << mask
                                              << " target=" << target
                                              << " value=" << boundary_value
                                              << '\n';
                                    printed_boundary_failure = true;
                                }
                            }
                            const int first_next
                                = labels.empty() ? 1 : labels.back();
                            for (int next_label = first_next;
                                 next_label < level; ++next_label) {
                                if (target == next_label) {
                                    continue;
                                }
                                cpp_int boundary = 0;
                                for_each_output(
                                    level, target, next_label, [&](int output) {
                                        boundary += defect[
                                            static_cast<std::size_t>(
                                                output * width + level
                                            )
                                        ];
                                    }
                                );
                                const cpp_int& lower = defect[
                                    static_cast<std::size_t>(
                                        target * width + level - next_label
                                    )
                                ];
                                const cpp_int& sink = defect[
                                    static_cast<std::size_t>(
                                        (level - target) * width + next_label
                                    )
                                ];
                                const cpp_int combined
                                    = boundary + lower - sink;
                                ++entries_tested;
                                if (combined < 0) {
                                    ++local_scbi_failures;
                                    if (!printed_scbi_failure) {
                                        std::cout
                                            << "local_scbi_failure level="
                                            << level << " prefix=";
                                        print_word(labels, 0U);
                                        std::cout << " allocation=" << mask
                                                  << " next=" << next_label
                                                  << " target=" << target
                                                  << " boundary=" << boundary
                                                  << " lower=" << lower
                                                  << " sink=" << sink
                                                  << " combined=" << combined
                                                  << '\n';
                                        printed_scbi_failure = true;
                                    }
                                }
                            }
                        }
                    }
                    for (std::size_t shell = 0; shell < shell_count; ++shell) {
                        const std::vector<cpp_int> defect
                            = fusion_hankel_defect(
                                level, shell_coefficients[shell]
                            );
                        for (int target = 0; target <= level; ++target) {
                            if (std::binary_search(
                                    labels.begin(), labels.end(), target
                                )) {
                                continue;
                            }
                            const cpp_int& boundary_value = defect[
                                static_cast<std::size_t>(
                                    target * width + level
                                )
                            ];
                            if (boundary_value < 0) {
                                ++shell_boundary_failures;
                                if (!printed_shell_boundary_failure) {
                                    std::cout
                                        << "shell_boundary_failure level="
                                        << level << " prefix=";
                                    print_word(labels, 0U);
                                    std::cout << " shell=" << shell
                                              << " target=" << target
                                              << " value=" << boundary_value
                                              << '\n';
                                    printed_shell_boundary_failure = true;
                                }
                            }
                            const int first_next
                                = labels.empty() ? 1 : labels.back();
                            for (int next_label = first_next;
                                 next_label < level; ++next_label) {
                                if (target == next_label) {
                                    continue;
                                }
                                cpp_int boundary = 0;
                                for_each_output(
                                    level, target, next_label, [&](int output) {
                                        boundary += defect[
                                            static_cast<std::size_t>(
                                                output * width + level
                                            )
                                        ];
                                    }
                                );
                                const cpp_int& lower = defect[
                                    static_cast<std::size_t>(
                                        target * width + level - next_label
                                    )
                                ];
                                const cpp_int& sink = defect[
                                    static_cast<std::size_t>(
                                        (level - target) * width + next_label
                                    )
                                ];
                                const cpp_int combined
                                    = boundary + lower - sink;
                                if (combined < 0) {
                                    ++shell_scbi_failures;
                                    if (!printed_shell_scbi_failure) {
                                        std::cout
                                            << "shell_scbi_failure level="
                                            << level << " prefix=";
                                        print_word(labels, 0U);
                                        std::cout << " shell=" << shell
                                                  << " next=" << next_label
                                                  << " target=" << target
                                                  << " boundary=" << boundary
                                                  << " lower=" << lower
                                                  << " sink=" << sink
                                                  << " combined=" << combined
                                                  << '\n';
                                        printed_shell_scbi_failure = true;
                                    }
                                }
                            }
                        }
                    }
                }
                std::cout << "progress level=" << level
                          << " words=" << words_tested
                          << " allocation_pairs="
                          << allocation_pairs_tested
                          << " entries=" << entries_tested
                          << " local_boundary_failures="
                          << local_boundary_failures
                          << " local_scbi_failures="
                          << local_scbi_failures
                          << " shell_boundary_failures="
                          << shell_boundary_failures
                          << " shell_scbi_failures="
                          << shell_scbi_failures << '\n';
            }
            std::cout << "SU2_SIMPLE_CURRENT_BOUNDARY_LOCAL levels=2.."
                      << maximum_level
                      << " maximum_prefix_factors=" << maximum_factors
                      << " words=" << words_tested
                      << " allocation_pairs=" << allocation_pairs_tested
                      << " entries=" << entries_tested
                      << " local_boundary_failures="
                      << local_boundary_failures
                      << " local_scbi_failures=" << local_scbi_failures
                      << " shell_boundary_failures="
                      << shell_boundary_failures
                      << " shell_scbi_failures=" << shell_scbi_failures
                      << " result="
                      << (local_boundary_failures == 0U
                              && local_scbi_failures == 0U
                              && shell_boundary_failures == 0U
                              && shell_scbi_failures == 0U
                          ? "PASS" : "FAIL")
                      << '\n';
            return local_boundary_failures == 0U
                    && local_scbi_failures == 0U
                    && shell_boundary_failures == 0U
                    && shell_scbi_failures == 0U
                ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        if (argc == 4
            && std::string(argv[1]) == "simple-current-boundary") {
            const int maximum_level = std::stoi(argv[2]);
            const int maximum_factors = std::stoi(argv[3]);
            if (maximum_level < 2 || maximum_factors < 0
                || maximum_factors >= 63) {
                throw std::runtime_error(
                    "invalid simple-current boundary bounds"
                );
            }
            std::uint64_t words_tested = 0U;
            std::uint64_t entries_tested = 0U;
            std::uint64_t boundary_only_failures = 0U;
            std::uint64_t lower_only_failures = 0U;
            std::uint64_t strict_boundary_failures = 0U;
            std::uint64_t strict_joint_failures = 0U;
            std::uint64_t fully_disjoint_joint_failures = 0U;
            bool printed_boundary_failure = false;
            bool printed_lower_failure = false;
            bool printed_strict_boundary_failure = false;
            bool printed_strict_joint_failure = false;
            bool printed_fully_disjoint_joint_failure = false;
            cpp_int minimum_combined = 0;
            bool have_combined = false;
            for (int level = 2; level <= maximum_level; ++level) {
                auto words = multisets(level - 1, maximum_factors);
                words.insert(words.begin(), std::vector<int>{});
                for (const auto& labels : words) {
                    const int width = level + 1;
                    const std::vector<cpp_int> coefficient
                        = signed_product(level, labels, 0U);
                    const std::vector<cpp_int> defect
                        = fusion_hankel_defect(level, coefficient);
                    ++words_tested;
                    const int first_next
                        = labels.empty() ? 1 : labels.back();
                    for (int next_label = first_next;
                         next_label < level; ++next_label) {
                        const bool strict_next = labels.empty()
                            || next_label > labels.back();
                        std::vector<int> extended = labels;
                        extended.push_back(next_label);
                        const std::vector<cpp_int> next_defect
                            = fusion_hankel_defect(
                                level,
                                signed_product(level, extended, 0U)
                        );
                        for (int target = 0; target <= level; ++target) {
                            if (std::binary_search(
                                    extended.begin(), extended.end(), target
                                )) {
                                continue;
                            }
                            cpp_int boundary = 0;
                            for_each_output(
                                level, target, next_label,
                                [&](int output) {
                                    boundary += defect[
                                        static_cast<std::size_t>(
                                            output * width + level
                                        )
                                    ];
                                }
                            );
                            const cpp_int& lower = defect[
                                static_cast<std::size_t>(
                                    target * width + level - next_label
                                )
                            ];
                            const cpp_int& sink = defect[
                                static_cast<std::size_t>(
                                    (level - target) * width + next_label
                                )
                            ];
                            const cpp_int combined = boundary + lower - sink;
                            const cpp_int& direct = next_defect[
                                static_cast<std::size_t>(
                                    target * width + level
                                )
                            ];
                            ++entries_tested;
                            if (combined != direct || direct < 0) {
                                std::cout
                                    << "SU2_SIMPLE_CURRENT_BOUNDARY "
                                    << "result=FAIL level=" << level
                                    << " prefix=";
                                print_word(labels, 0U);
                                std::cout << " next=" << next_label
                                          << " target=" << target
                                          << " boundary=" << boundary
                                          << " lower=" << lower
                                          << " sink=" << sink
                                          << " combined=" << combined
                                          << " direct=" << direct << '\n';
                                return EXIT_FAILURE;
                            }
                            if (!have_combined || combined < minimum_combined) {
                                minimum_combined = combined;
                                have_combined = true;
                            }
                            if (boundary < sink) {
                                ++boundary_only_failures;
                                if (!printed_boundary_failure) {
                                    std::cout
                                        << "boundary_only_failure level="
                                        << level << " prefix=";
                                    print_word(labels, 0U);
                                    std::cout << " next=" << next_label
                                              << " target=" << target
                                              << " boundary=" << boundary
                                              << " lower=" << lower
                                              << " sink=" << sink << '\n';
                                    printed_boundary_failure = true;
                                }
                                if (strict_next) {
                                    ++strict_boundary_failures;
                                    if (!printed_strict_boundary_failure) {
                                        std::cout
                                            << "strict_boundary_only_failure"
                                            << " level=" << level
                                            << " prefix=";
                                        print_word(labels, 0U);
                                        std::cout << " next=" << next_label
                                                  << " target=" << target
                                                  << " boundary=" << boundary
                                                  << " lower=" << lower
                                                  << " sink=" << sink << '\n';
                                        printed_strict_boundary_failure = true;
                                    }
                                }
                            }
                            if (lower < sink) {
                                ++lower_only_failures;
                                if (!printed_lower_failure) {
                                    std::cout
                                        << "lower_only_failure level="
                                        << level << " prefix=";
                                    print_word(labels, 0U);
                                    std::cout << " next=" << next_label
                                              << " target=" << target
                                              << " boundary=" << boundary
                                              << " lower=" << lower
                                              << " sink=" << sink << '\n';
                                    printed_lower_failure = true;
                                }
                            }
                            if (strict_next && boundary < sink
                                && lower < sink) {
                                ++strict_joint_failures;
                                if (!printed_strict_joint_failure) {
                                    std::cout << "strict_joint_failure"
                                              << " level=" << level
                                              << " prefix=";
                                    print_word(labels, 0U);
                                    std::cout << " next=" << next_label
                                              << " target=" << target
                                              << " boundary=" << boundary
                                              << " lower=" << lower
                                              << " sink=" << sink << '\n';
                                    printed_strict_joint_failure = true;
                                }
                                if (!std::binary_search(
                                        labels.begin(), labels.end(),
                                        level - target
                                    )) {
                                    ++fully_disjoint_joint_failures;
                                    if (!printed_fully_disjoint_joint_failure) {
                                        std::cout
                                            << "fully_disjoint_joint_failure"
                                            << " level=" << level
                                            << " prefix=";
                                        print_word(labels, 0U);
                                        std::cout << " next=" << next_label
                                                  << " target=" << target
                                                  << " boundary=" << boundary
                                                  << " lower=" << lower
                                                  << " sink=" << sink << '\n';
                                        printed_fully_disjoint_joint_failure
                                            = true;
                                    }
                                }
                            }
                        }
                    }
                }
                std::cout << "progress level=" << level
                          << " words=" << words_tested
                          << " entries=" << entries_tested
                          << " boundary_only_failures="
                          << boundary_only_failures
                          << " lower_only_failures="
                          << lower_only_failures
                          << " strict_boundary_only_failures="
                          << strict_boundary_failures
                          << " strict_joint_failures="
                          << strict_joint_failures
                          << " fully_disjoint_joint_failures="
                          << fully_disjoint_joint_failures << '\n';
            }
            std::cout << "SU2_SIMPLE_CURRENT_BOUNDARY levels=2.."
                      << maximum_level
                      << " maximum_prefix_factors=" << maximum_factors
                      << " words=" << words_tested
                      << " entries=" << entries_tested
                      << " boundary_only_failures="
                      << boundary_only_failures
                      << " lower_only_failures=" << lower_only_failures
                      << " strict_boundary_only_failures="
                      << strict_boundary_failures
                      << " strict_joint_failures="
                      << strict_joint_failures
                      << " fully_disjoint_joint_failures="
                      << fully_disjoint_joint_failures
                      << " minimum_combined="
                      << (have_combined ? minimum_combined : cpp_int{0})
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc == 4
            && std::string(argv[1]) == "simple-current-wall") {
            const int maximum_level = std::stoi(argv[2]);
            const int maximum_factors = std::stoi(argv[3]);
            if (maximum_level < 2 || maximum_factors < 0
                || maximum_factors >= 63) {
                throw std::runtime_error(
                    "invalid simple-current wall bounds"
                );
            }
            std::uint64_t total_tested = 0U;
            std::uint64_t total_active = 0U;
            std::uint64_t total_tight = 0U;
            std::uint64_t printed_tight = 0U;
            cpp_int minimum_slack = 0;
            bool have_slack = false;
            const std::size_t sector_count
                = static_cast<std::size_t>(maximum_factors + 1);
            std::vector<std::uint64_t> sector_tested(sector_count, 0U);
            std::vector<std::uint64_t> sector_active(sector_count, 0U);
            std::vector<std::uint64_t> sector_tight(sector_count, 0U);
            std::vector<cpp_int> sector_minimum(sector_count, 0);
            std::vector<bool> sector_has_minimum(sector_count, false);
            for (int level = 2; level <= maximum_level; ++level) {
                auto words = multisets(level - 1, maximum_factors);
                words.insert(words.begin(), std::vector<int>{});
                for (const auto& labels : words) {
                    const std::uint64_t masks
                        = std::uint64_t{1} << labels.size();
                    for (std::uint64_t mask = 0U; mask < masks; ++mask) {
                        if ((popcount(mask) & 1) != 0) {
                            continue;
                        }
                        const int width = level + 1;
                        const std::vector<cpp_int> product
                            = signed_product(level, labels, mask);
                        const cpp_int& h_00 = product[0];
                        const cpp_int& h_kk = product[
                            static_cast<std::size_t>(
                                level * width + level
                            )
                        ];
                        const cpp_int slack = h_00 - h_kk;
                        const std::size_t minus_count
                            = static_cast<std::size_t>(popcount(mask));
                        ++total_tested;
                        ++sector_tested[minus_count];
                        if (h_00 < 0 || slack < 0) {
                            std::cout
                                << "SU2_SIMPLE_CURRENT_WALL result=FAIL"
                                << " level=" << level << " remainder=";
                            print_word(labels, mask);
                            std::cout << " h_00=" << h_00
                                      << " h_kk=" << h_kk
                                      << " slack=" << slack << '\n';
                            return EXIT_FAILURE;
                        }
                        if (h_kk == 0) {
                            continue;
                        }
                        ++total_active;
                        ++sector_active[minus_count];
                        if (!have_slack || slack < minimum_slack) {
                            minimum_slack = slack;
                            have_slack = true;
                        }
                        if (!sector_has_minimum[minus_count]
                            || slack < sector_minimum[minus_count]) {
                            sector_minimum[minus_count] = slack;
                            sector_has_minimum[minus_count] = true;
                        }
                        if (slack == 0) {
                            ++total_tight;
                            ++sector_tight[minus_count];
                            if (printed_tight < 12U) {
                                std::cout << "tight_simple_current_wall"
                                          << " level=" << level
                                          << " remainder=";
                                print_word(labels, mask);
                                std::cout << " h_00=" << h_00
                                          << " h_kk=" << h_kk << '\n';
                                ++printed_tight;
                            }
                        }
                    }
                }
                std::cout << "progress level=" << level
                          << " tested=" << total_tested
                          << " active=" << total_active
                          << " tight=" << total_tight << '\n';
            }
            for (std::size_t minus_count = 0U;
                 minus_count < sector_count; ++minus_count) {
                if (sector_tested[minus_count] == 0U) {
                    continue;
                }
                std::cout << "minus_sector count=" << minus_count
                          << " tested=" << sector_tested[minus_count]
                          << " active=" << sector_active[minus_count]
                          << " tight=" << sector_tight[minus_count]
                          << " minimum_active_slack="
                          << (sector_has_minimum[minus_count]
                                  ? sector_minimum[minus_count]
                                  : cpp_int{0})
                          << '\n';
            }
            std::cout << "SU2_SIMPLE_CURRENT_WALL levels=2.."
                      << maximum_level
                      << " maximum_remainder_factors=" << maximum_factors
                      << " tested=" << total_tested
                      << " active=" << total_active
                      << " tight=" << total_tight
                      << " minimum_active_slack="
                      << (have_slack ? minimum_slack : cpp_int{0})
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: search_su2_level_partial_character_cone "
                "MAXIMUM_LEVEL MAXIMUM_FACTORS\n"
                "   or: search_su2_level_partial_character_cone "
                "simple-current-wall MAXIMUM_LEVEL MAXIMUM_FACTORS\n"
                "   or: search_su2_level_partial_character_cone "
                "simple-current-boundary MAXIMUM_LEVEL MAXIMUM_FACTORS\n"
                "   or: search_su2_level_partial_character_cone "
                "simple-current-boundary-local MAXIMUM_LEVEL "
                "MAXIMUM_FACTORS\n"
                "   or: search_su2_level_partial_character_cone "
                "simple-current-exterior MAXIMUM_LEVEL MAXIMUM_FACTORS\n"
                "   or: search_su2_level_partial_character_cone "
                "finite-c2-fundamental MAXIMUM_LEVEL MAXIMUM_FACTORS\n"
                "   or: search_su2_level_partial_character_cone "
                "finite-c2-wall-generator MAXIMUM_LEVEL MAXIMUM_LABEL\n"
                "   or: search_su2_level_partial_character_cone "
                "finite-c2-wall-packets MAXIMUM_LEVEL MAXIMUM_FACTORS\n"
                "   or: search_su2_level_partial_character_cone "
                "finite-c2-wall-cone MAXIMUM_LEVEL MAXIMUM_FACTORS"
            );
        }
        const int maximum_level = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_level < 1 || maximum_factors < 1
            || maximum_factors >= 63) {
            throw std::runtime_error("invalid search bound");
        }

        std::uint64_t total_tested = 0;
        for (int level = 1; level <= maximum_level; ++level) {
            const auto words = multisets(level, maximum_factors);
            std::uint64_t level_tested = 0;
            for (const auto& labels : words) {
                const std::uint64_t masks
                    = std::uint64_t{1} << labels.size();
                for (std::uint64_t mask = 0; mask < masks; ++mask) {
                    int failed_target = -1;
                    cpp_int failed_value = 0;
                    ++level_tested;
                    if (!inspect_word(
                            level, labels, mask,
                            failed_target, failed_value
                        )) {
                        std::cout << "FAIL level=" << level << " factors=";
                        print_word(labels, mask);
                        std::cout << " target=" << failed_target
                                  << " coefficient=" << failed_value << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
            total_tested += level_tested;
            std::cout << "SU2_LEVEL_PARTIAL_CHARACTER_CONE level=" << level
                      << " tested=" << level_tested << " PASS\n";
        }
        std::cout << "SU2_LEVEL_PARTIAL_CHARACTER_CONE total_tested="
                  << total_tested << " maximum_level=" << maximum_level
                  << " maximum_factors=" << maximum_factors << " PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
