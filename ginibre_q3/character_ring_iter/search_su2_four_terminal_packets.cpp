#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Character = std::vector<cpp_int>;

template <class Function>
void for_each_output(
    const int level,
    const int left,
    const int right,
    Function function
) {
    const int maximum = std::min(left + right, 2 * level - left - right);
    for (int output = std::abs(left - right);
         output <= maximum;
         output += 2) {
        function(output);
    }
}

Character tensor_multiplicities(
    const int level,
    const std::vector<int>& labels
) {
    Character current(static_cast<std::size_t>(level + 1));
    Character next(static_cast<std::size_t>(level + 1));
    current[0] = 1;
    for (const int label : labels) {
        std::fill(next.begin(), next.end(), cpp_int(0));
        for (int source = 0; source <= level; ++source) {
            const cpp_int value = current[static_cast<std::size_t>(source)];
            if (value == 0) {
                continue;
            }
            for_each_output(level, source, label, [&](const int output) {
                next[static_cast<std::size_t>(output)] += value;
            });
        }
        current.swap(next);
    }
    return current;
}

cpp_int transition_multiplicity(
    const int level,
    const std::vector<int>& even_word,
    const unsigned long long even_mask,
    const int start,
    const int end,
    const std::vector<int>& extra_labels = {}
) {
    std::vector<int> labels;
    labels.reserve(1U + even_word.size() + extra_labels.size());
    labels.push_back(start);
    for (std::size_t index = 0; index < even_word.size(); ++index) {
        if (((even_mask >> index) & 1ULL) != 0ULL) {
            labels.push_back(even_word[index]);
        }
    }
    labels.insert(labels.end(), extra_labels.begin(), extra_labels.end());
    return tensor_multiplicities(level, labels)
        [static_cast<std::size_t>(end)];
}

cpp_int invariant_multiplicity(
    const int level,
    const std::vector<int>& even_word,
    const unsigned long long even_mask,
    const std::vector<int>& odd_labels,
    const unsigned int odd_mask
) {
    std::vector<int> labels;
    labels.reserve(even_word.size() + odd_labels.size());
    for (std::size_t index = 0; index < even_word.size(); ++index) {
        if (((even_mask >> index) & 1ULL) != 0ULL) {
            labels.push_back(even_word[index]);
        }
    }
    for (std::size_t index = 0; index < odd_labels.size(); ++index) {
        if (((odd_mask >> index) & 1U) != 0U) {
            labels.push_back(odd_labels[index]);
        }
    }
    return tensor_multiplicities(level, labels)[0];
}

void print_word(const std::vector<int>& word) {
    std::cout << '[';
    for (std::size_t index = 0; index < word.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << word[index];
    }
    std::cout << ']';
}

void print_selected(
    const std::vector<int>& word,
    const unsigned long long mask
) {
    std::cout << '[';
    bool first = true;
    for (std::size_t index = 0; index < word.size(); ++index) {
        if (((mask >> index) & 1ULL) == 0ULL) {
            continue;
        }
        if (!first) {
            std::cout << ',';
        }
        first = false;
        std::cout << word[index];
    }
    std::cout << ']';
}

struct Packet {
    cpp_int total = 0;
    cpp_int zero_or_four = 0;
    cpp_int positive_pair = 0;
    cpp_int crossed_pair = 0;
};

Packet packet(
    const int level,
    const std::vector<int>& even_word,
    const unsigned long long even_mask,
    const int q,
    const int r,
    const int a
) {
    // Odd terminals are ordered as q,r,a,1.  The first two carry minus signs
    // after the central parity toggle; a and 1 carry plus signs.
    const std::vector<int> odd_labels{q, r, a, 1};
    constexpr unsigned int full_odd_mask = 15U;
    const unsigned long long full_even_mask = even_word.empty()
        ? 0ULL
        : (1ULL << even_word.size()) - 1ULL;
    Packet answer;
    for (unsigned int odd_mask = 0U; odd_mask <= full_odd_mask; ++odd_mask) {
        const int odd_count = __builtin_popcount(odd_mask);
        if ((odd_count & 1) != 0) {
            continue;
        }
        const cpp_int left = invariant_multiplicity(
            level, even_word, even_mask, odd_labels, odd_mask
        );
        const cpp_int right = invariant_multiplicity(
            level,
            even_word,
            full_even_mask ^ even_mask,
            odd_labels,
            full_odd_mask ^ odd_mask
        );
        const cpp_int value = left * right;
        const bool negative = (__builtin_popcount(odd_mask & 3U) & 1) != 0;
        answer.total += negative ? -value : value;
        if (odd_count == 0 || odd_count == 4) {
            answer.zero_or_four += value;
        } else if (negative) {
            answer.crossed_pair += value;
        } else {
            answer.positive_pair += value;
        }
    }
    return answer;
}

cpp_int oriented_carrier_packet(
    const int level,
    const std::vector<int>& even_word,
    const unsigned long long first_mask,
    const int q,
    const int r,
    const int a
) {
    const unsigned long long full_mask = even_word.empty()
        ? 0ULL
        : (1ULL << even_word.size()) - 1ULL;
    const unsigned long long second_mask = full_mask ^ first_mask;
    const cpp_int first_positive = transition_multiplicity(
        level, even_word, second_mask, 0, 0
    ) * transition_multiplicity(
        level, even_word, first_mask, 0, a, {q, r, 1}
    );
    const cpp_int second_positive = transition_multiplicity(
        level, even_word, first_mask, 1, a
    ) * transition_multiplicity(
        level, even_word, second_mask, 0, 0, {q, r}
    );
    const cpp_int first_crossed = transition_multiplicity(
        level, even_word, first_mask, r, a
    ) * transition_multiplicity(
        level, even_word, second_mask, q, 1
    );
    const cpp_int second_crossed = transition_multiplicity(
        level, even_word, first_mask, q, a
    ) * transition_multiplicity(
        level, even_word, second_mask, r, 1
    );
    return first_positive + second_positive
        - first_crossed - second_crossed;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4 && argc != 5) {
            throw std::runtime_error(
                "usage: search_su2_four_terminal_packets "
                "MINIMUM_LEVEL MAXIMUM_LEVEL MAXIMUM_EVEN_FACTORS "
                "[--largest-pair|--largest-block]"
            );
        }
        const int minimum_level = std::stoi(argv[1]);
        const int maximum_level = std::stoi(argv[2]);
        const int maximum_factors = std::stoi(argv[3]);
        const bool largest_pair_mode = argc == 5
            && std::string(argv[4]) == "--largest-pair";
        const bool largest_block_mode = argc == 5
            && std::string(argv[4]) == "--largest-block";
        if (argc == 5 && !largest_pair_mode && !largest_block_mode) {
            throw std::runtime_error("unknown mode");
        }
        if (minimum_level < 2 || maximum_level < minimum_level
            || maximum_factors < 0
            || maximum_factors >= 63) {
            throw std::runtime_error("invalid search bound");
        }

        long long words = 0;
        long long packets = 0;
        long long carrier_identities = 0;
        bool found_failure = false;
        for (int level = minimum_level;
             level <= maximum_level && !found_failure;
             ++level) {
            std::vector<int> word;
            const auto visit = [&](const auto& self, const int next_even) -> void {
                if (found_failure) {
                    return;
                }
                ++words;
                const unsigned long long allocation_count =
                    1ULL << word.size();
                for (int q = 1; q <= level && !found_failure; q += 2) {
                    for (int r = q; r <= level && !found_failure; r += 2) {
                        for (int a = 1; a <= level && !found_failure; a += 2) {
                            if (a == 1 || q == 1 || r == 1
                                || q == a || r == a) {
                                continue;
                            }
                            if ((largest_pair_mode || largest_block_mode)
                                && word.empty()) {
                                continue;
                            }
                            std::size_t floating_start = word.size();
                            if (largest_pair_mode && !word.empty()) {
                                floating_start = word.size() - 1U;
                            } else if (largest_block_mode && !word.empty()) {
                                floating_start = word.size() - 1U;
                                while (floating_start > 0U
                                    && word[floating_start - 1U]
                                        == word.back()) {
                                    --floating_start;
                                }
                            }
                            const unsigned long long tested_allocations =
                                (largest_pair_mode || largest_block_mode)
                                    ? 1ULL << floating_start
                                    : allocation_count;
                            const std::size_t floating_count =
                                word.size() - floating_start;
                            const unsigned long long floating_allocations =
                                1ULL << floating_count;
                            for (unsigned long long mask = 0ULL;
                                 mask < tested_allocations;
                                 ++mask) {
                                ++packets;
                                Packet value = packet(
                                    level, word, mask, q, r, a
                                );
                                const unsigned long long full_mask =
                                    allocation_count - 1ULL;
                                const cpp_int carrier_pair =
                                    oriented_carrier_packet(
                                        level, word, mask, q, r, a
                                    ) + oriented_carrier_packet(
                                        level,
                                        word,
                                        full_mask ^ mask,
                                        q,
                                        r,
                                        a
                                    );
                                if (value.total != carrier_pair) {
                                    throw std::runtime_error(
                                        "carrier-switch identity failed"
                                    );
                                }
                                ++carrier_identities;
                                if (largest_pair_mode || largest_block_mode) {
                                    value = Packet{};
                                    for (unsigned long long floating_mask = 0ULL;
                                         floating_mask < floating_allocations;
                                         ++floating_mask) {
                                        const Packet member = packet(
                                            level,
                                            word,
                                            mask | (floating_mask
                                                << floating_start),
                                            q,
                                            r,
                                            a
                                        );
                                        const unsigned long long member_mask =
                                            mask | (floating_mask
                                                << floating_start);
                                        const cpp_int member_carrier =
                                            oriented_carrier_packet(
                                                level,
                                                word,
                                                member_mask,
                                                q,
                                                r,
                                                a
                                            ) + oriented_carrier_packet(
                                                level,
                                                word,
                                                full_mask ^ member_mask,
                                                q,
                                                r,
                                                a
                                            );
                                        if (member.total != member_carrier) {
                                            throw std::runtime_error(
                                                "carrier-switch orbit identity failed"
                                            );
                                        }
                                        ++carrier_identities;
                                        value.total += member.total;
                                        value.zero_or_four +=
                                            member.zero_or_four;
                                        value.positive_pair +=
                                            member.positive_pair;
                                        value.crossed_pair +=
                                            member.crossed_pair;
                                    }
                                }
                                if (value.total >= 0) {
                                    continue;
                                }
                                found_failure = true;
                                std::cout << ((largest_pair_mode
                                    || largest_block_mode)
                                    ? "LARGEST_ORBIT_FAIL level="
                                    : "FIXED_ALLOCATION_FAIL level=")
                                          << level << " even_word=";
                                print_word(word);
                                std::cout << " first_even=";
                                print_selected(word, mask);
                                std::cout << " q=" << q << " r=" << r
                                          << " a=" << a
                                          << " zero_or_four="
                                          << value.zero_or_four
                                          << " positive_pair="
                                          << value.positive_pair
                                          << " crossed_pair="
                                          << value.crossed_pair
                                          << " total=" << value.total << '\n';
                                if ((largest_pair_mode || largest_block_mode)
                                    && !word.empty()) {
                                    int first_working_suffix = -1;
                                    for (std::size_t suffix = 1U;
                                         suffix <= word.size();
                                         ++suffix) {
                                        const std::size_t prefix_count =
                                            word.size() - suffix;
                                        const unsigned long long prefix_masks =
                                            1ULL << prefix_count;
                                        const unsigned long long suffix_masks =
                                            1ULL << suffix;
                                        bool all_nonnegative = true;
                                        for (unsigned long long prefix_mask = 0ULL;
                                             prefix_mask < prefix_masks;
                                             ++prefix_mask) {
                                            cpp_int orbit = 0;
                                            for (unsigned long long suffix_mask = 0ULL;
                                                 suffix_mask < suffix_masks;
                                                 ++suffix_mask) {
                                                orbit += packet(
                                                    level,
                                                    word,
                                                    prefix_mask
                                                        | (suffix_mask
                                                            << prefix_count),
                                                    q,
                                                    r,
                                                    a
                                                ).total;
                                            }
                                            if (orbit < 0) {
                                                all_nonnegative = false;
                                                break;
                                            }
                                        }
                                        if (all_nonnegative) {
                                            first_working_suffix =
                                                static_cast<int>(suffix);
                                            break;
                                        }
                                    }
                                    std::cout << "FIRST_WORKING_SUFFIX factors="
                                              << first_working_suffix << '\n';
                                }
                                if (!largest_pair_mode && !largest_block_mode) {
                                    cpp_int aggregate = 0;
                                    std::cout << "ALLOCATION_PROFILE";
                                    for (unsigned long long profile_mask = 0ULL;
                                         profile_mask < allocation_count;
                                         ++profile_mask) {
                                        const Packet profile = packet(
                                            level, word, profile_mask, q, r, a
                                        );
                                        aggregate += profile.total;
                                        std::cout << ' ' << profile_mask << ':'
                                                  << profile.total;
                                    }
                                    std::cout << " aggregate=" << aggregate
                                              << '\n';
                                }
                                break;
                            }
                        }
                    }
                }
                if (found_failure
                    || word.size()
                        == static_cast<std::size_t>(maximum_factors)) {
                    return;
                }
                for (int label = next_even; label <= level; label += 2) {
                    word.push_back(label);
                    self(self, label);
                    word.pop_back();
                    if (found_failure) {
                        return;
                    }
                }
            };
            visit(visit, 2);
        }
        if (found_failure) {
            return EXIT_SUCCESS;
        }
        std::cout << ((largest_pair_mode || largest_block_mode)
                          ? "SU2_FOUR_TERMINAL_LARGEST_ORBITS PASS words="
                          : "SU2_FOUR_TERMINAL_FIXED_ALLOCATIONS PASS words=")
                  << words << " packets=" << packets
                  << " carrier_identities=" << carrier_identities
                  << " minimum_level=" << minimum_level
                  << " maximum_level=" << maximum_level
                  << " maximum_even_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
