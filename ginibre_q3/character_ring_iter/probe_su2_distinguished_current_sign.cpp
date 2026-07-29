#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value <= 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error(
            std::string(name) + " must be a positive integer"
        );
    }
    return static_cast<int>(value);
}

std::vector<Integer> multiply(
    const std::vector<Integer>& profile,
    int q
) {
    std::vector<Integer> result(
        profile.size() + static_cast<std::size_t>(q)
    );
    for (std::size_t source = 0; source < profile.size(); ++source) {
        if (profile[source] == 0) {
            continue;
        }
        const int i = static_cast<int>(source);
        for (int target = std::abs(i - q);
             target <= i + q;
             ++target) {
            result[static_cast<std::size_t>(target)] += profile[source];
        }
    }
    return result;
}

int sign(const Integer& value) {
    return value > 0 ? 1 : (value < 0 ? -1 : 0);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: probe_su2_distinguished_current_sign "
                "maximum_half_label maximum_power"
            );
        }
        const int maximum_q =
            parse_positive(argv[1], "maximum_half_label");
        const int maximum_power =
            parse_positive(argv[2], "maximum_power");

        std::size_t profiles = 0;
        std::size_t adjacent_minors = 0;
        std::size_t negative_minors = 0;
        std::size_t positive_to_negative = 0;
        std::size_t negative_to_positive = 0;
        int maximum_sign_changes = 0;
        int maximum_negative_blocks = 0;
        int first_negative_q = 0;
        int first_negative_power = 0;
        int first_negative_index = 0;
        Integer first_negative_value = 0;
        std::vector<int> negative_profiles_by_power(
            static_cast<std::size_t>(maximum_power + 1)
        );
        std::vector<int> current_decrease_profiles_by_power(
            static_cast<std::size_t>(maximum_power + 1)
        );
        std::size_t current_comparisons = 0;
        std::size_t current_decreases = 0;

        for (int q = 1; q <= maximum_q; ++q) {
            std::vector<Integer> profile{Integer(1)};
            for (int power = 1; power <= maximum_power; ++power) {
                profile = multiply(profile, q);
                const std::vector<Integer> next = multiply(profile, q);
                int previous_sign = 0;
                int sign_changes = 0;
                int negative_blocks = 0;
                bool inside_negative_block = false;
                bool profile_has_negative = false;
                for (std::size_t index = 1;
                     index < profile.size();
                     ++index) {
                    if (profile[index - 1U] == 0 || profile[index] == 0) {
                        continue;
                    }
                    const Integer minor =
                        next[index] * profile[index - 1U]
                        - next[index - 1U] * profile[index];
                    ++adjacent_minors;
                    const int current_sign = sign(minor);
                    if (current_sign < 0) {
                        ++negative_minors;
                        profile_has_negative = true;
                        if (first_negative_q == 0) {
                            first_negative_q = q;
                            first_negative_power = power;
                            first_negative_index =
                                static_cast<int>(index);
                            first_negative_value = minor;
                        }
                        if (!inside_negative_block) {
                            ++negative_blocks;
                            inside_negative_block = true;
                        }
                    } else if (current_sign > 0) {
                        inside_negative_block = false;
                    }
                    if (
                        current_sign != 0
                        && previous_sign != 0
                        && current_sign != previous_sign
                    ) {
                        ++sign_changes;
                        if (previous_sign > 0) {
                            ++positive_to_negative;
                        } else {
                            ++negative_to_positive;
                        }
                    }
                    if (current_sign != 0) {
                        previous_sign = current_sign;
                    }
                }
                maximum_sign_changes =
                    std::max(maximum_sign_changes, sign_changes);
                maximum_negative_blocks =
                    std::max(maximum_negative_blocks, negative_blocks);
                if (profile_has_negative) {
                    ++negative_profiles_by_power[
                        static_cast<std::size_t>(power)
                    ];
                }
                bool profile_has_current_decrease = false;
                const int support =
                    static_cast<int>(profile.size()) - 1;
                for (int gap = 1; gap <= support + q; ++gap) {
                    Integer previous = 0;
                    bool have_previous = false;
                    for (int left = 0; left <= support; ++left) {
                        const int right = left + gap;
                        const Integer profile_right =
                            right <= support
                                ? profile[static_cast<std::size_t>(right)]
                                : Integer(0);
                        const Integer next_left =
                            next[static_cast<std::size_t>(left)];
                        const Integer next_right =
                            right < static_cast<int>(next.size())
                                ? next[static_cast<std::size_t>(right)]
                                : Integer(0);
                        const Integer current =
                            profile[static_cast<std::size_t>(left)]
                                * next_right
                            - profile_right * next_left;
                        if (have_previous) {
                            ++current_comparisons;
                            if (current < previous) {
                                ++current_decreases;
                                profile_has_current_decrease = true;
                            }
                        }
                        previous = current;
                        have_previous = true;
                    }
                }
                if (profile_has_current_decrease) {
                    ++current_decrease_profiles_by_power[
                        static_cast<std::size_t>(power)
                    ];
                }
                ++profiles;
            }
        }

        std::cout
            << "SU2_DISTINGUISHED_CURRENT_SIGN"
            << " maximum_half_label=" << maximum_q
            << " maximum_power=" << maximum_power
            << " profiles=" << profiles
            << " adjacent_minors=" << adjacent_minors
            << " negative_minors=" << negative_minors
            << " positive_to_negative=" << positive_to_negative
            << " negative_to_positive=" << negative_to_positive
            << " maximum_sign_changes=" << maximum_sign_changes
            << " maximum_negative_blocks=" << maximum_negative_blocks
            << " current_comparisons=" << current_comparisons
            << " current_decreases=" << current_decreases
            << " first_negative={q=" << first_negative_q
            << " power=" << first_negative_power
            << " index=" << first_negative_index
            << " value=" << first_negative_value << "}";
        std::cout << " negative_profiles_by_power={";
        bool first = true;
        for (int power = 1; power <= maximum_power; ++power) {
            const int count = negative_profiles_by_power[
                static_cast<std::size_t>(power)
            ];
            if (count == 0) {
                continue;
            }
            if (!first) {
                std::cout << ',';
            }
            first = false;
            std::cout << power << ':' << count;
        }
        std::cout << "} current_decrease_profiles_by_power={";
        first = true;
        for (int power = 1; power <= maximum_power; ++power) {
            const int count = current_decrease_profiles_by_power[
                static_cast<std::size_t>(power)
            ];
            if (count == 0) {
                continue;
            }
            if (!first) {
                std::cout << ',';
            }
            first = false;
            std::cout << power << ':' << count;
        }
        std::cout << "}\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
