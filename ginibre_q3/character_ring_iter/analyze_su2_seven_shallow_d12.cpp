#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#define SU2_SEVEN_PARITY_ORBITS_NO_MAIN
#include "analyze_su2_seven_parity_orbits.cpp"
#undef SU2_SEVEN_PARITY_ORBITS_NO_MAIN

namespace {

struct SmallRankWitness {
    bool initialized = false;
    std::int64_t value = 0;
    int level = 0;
    int orbit_index = 0;
    unsigned int selected_mask = 0U;
    int endpoint_position = 0;
    Labels labels;
    std::int64_t local = 0;
    std::int64_t pair = 0;
    std::int64_t positive = 0;
    std::int64_t negative = 0;
    std::array<std::int64_t, 5> selected_profile{};
    std::array<std::int64_t, 5> complement_profile{};
    int active_rank_one = 0;
    int active_rank_two = 0;
};

struct SmallRankSummary {
    std::uint64_t cases = 0;
    SmallRankWitness direct;
    SmallRankWitness ceiling;
    SmallRankWitness local;
};

void consider_small_rank(
    SmallRankWitness& witness,
    std::int64_t value,
    int level,
    int orbit_index,
    unsigned int selected_mask,
    int endpoint_position,
    const Labels& labels,
    std::int64_t local,
    std::int64_t pair,
    std::int64_t positive,
    std::int64_t negative,
    const std::array<std::int64_t, 5>& selected_profile,
    const std::array<std::int64_t, 5>& complement_profile,
    int active_rank_one,
    int active_rank_two
);

void print_small_rank_witness(
    const char* name,
    const SmallRankWitness& witness
);

std::array<SmallRankWitness, 4> double_neighbor_margin;
std::array<SmallRankWitness, 4> double_neighbor_w_payment;
std::array<SmallRankWitness, 4> double_neighbor_tail_payment;
std::array<std::array<SmallRankWitness, 4>, 7>
    neighbor_local_with_payment;
std::array<std::array<SmallRankWitness, 4>, 7>
    neighbor_pair_without_payment;
std::array<std::array<SmallRankWitness, 4>, 7>
    neighbor_local_without_payment;
std::array<std::array<SmallRankWitness, 4>, 7>
    neighbor_pair_with_payment;
std::array<std::array<SmallRankWitness, 4>, 7>
    neighbor_local_total_payment;
std::array<std::array<SmallRankSummary, 4>, 2>
    neighbor_zero_selected_by_sign;
std::array<std::array<SmallRankWitness, 4>, 2>
    neighbor_zero_selected_supply;
std::array<std::array<SmallRankWitness, 4>, 2>
    neighbor_zero_selected_local;
std::array<std::array<SmallRankWitness, 4>, 2>
    neighbor_zero_selected_pair;
SmallRankWitness rank_one_switch_margin;
std::array<SmallRankWitness, 3> rank_one_switch_component_margin;
SmallRankWitness rank_one_switch_bad_maximum;
std::array<SmallRankWitness, 3> rank_one_switch_component_bad_maximum;
std::array<std::array<SmallRankWitness, 2>, 3>
    rank_one_switch_endpoint_margin;
std::uint64_t rank_one_switch_cases = 0;

void scan_no_two_selected_profile(int maximum_level) {
    std::int64_t minimum = std::numeric_limits<std::int64_t>::max();
    int minimum_level = 0;
    int minimum_label = 0;
    int minimum_lower = 0;
    int minimum_upper = 0;
    std::array<std::int64_t, 5> minimum_profile{};
    std::array<std::int64_t, 3> no_two_weight_minimum{
        std::numeric_limits<std::int64_t>::max(),
        std::numeric_limits<std::int64_t>::max(),
        std::numeric_limits<std::int64_t>::max()
    };
    std::array<std::array<std::int64_t, 5>, 3>
        no_two_weight_profile{};
    std::array<std::int64_t, 3> one_two_minimum{
        std::numeric_limits<std::int64_t>::max(),
        std::numeric_limits<std::int64_t>::max(),
        std::numeric_limits<std::int64_t>::max()
    };
    std::array<std::array<std::int64_t, 5>, 3>
        one_two_profile{};
    std::array<std::array<std::int64_t, 3>, 3>
        factorized_minimum{};
    for (auto& row : factorized_minimum) {
        row.fill(std::numeric_limits<std::int64_t>::max());
    }
    std::uint64_t cases = 0;
    for (int level = 7; level <= maximum_level; level += 2) {
        for (int label = 4; label <= level - 3; label += 2) {
            for (int lower = 0; lower <= label; lower += 2) {
                for (int upper = label;
                     upper <= level - 1; upper += 2) {
                    if ((upper - lower) / 2 + 1 < 4) {
                        continue;
                    }
                    ++cases;
                    std::array<std::int64_t, 5> profile{};
                    for (int intermediate = lower;
                         intermediate <= upper;
                         intermediate += 2) {
                        for (std::size_t output = 0U;
                             output < profile.size(); ++output) {
                            profile[output] += multiplicity(
                                Labels{intermediate, label},
                                static_cast<int>(2U * output),
                                level
                            );
                        }
                    }
                    const std::int64_t weighted =
                        2 * profile[0]
                        + 5 * profile[1]
                        + 5 * profile[2]
                        + 3 * profile[3]
                        + profile[4];
                    const std::array<std::int64_t, 3> weights{
                        2 * profile[0]
                            + 6 * profile[1]
                            + 7 * profile[2]
                            + 4 * profile[3],
                        2 * profile[0]
                            + 5 * profile[1]
                            + 6 * profile[2]
                            + 4 * profile[3],
                        weighted
                    };
                    for (std::size_t index = 0U;
                         index < weights.size(); ++index) {
                        if (weights[index]
                            < no_two_weight_minimum[index]) {
                            no_two_weight_minimum[index] =
                                weights[index];
                            no_two_weight_profile[index] = profile;
                        }
                    }
                    if (weighted < minimum) {
                        minimum = weighted;
                        minimum_level = level;
                        minimum_label = label;
                        minimum_lower = lower;
                        minimum_upper = upper;
                        minimum_profile = profile;
                    }
                }
            }
        }
    }
    for (int level = 5; level <= maximum_level; level += 2) {
        std::array<
            std::vector<std::array<std::int64_t, 5>>,
            3
        > selected_profiles;
        std::array<
            std::vector<std::array<std::int64_t, 5>>,
            3
        > complement_profiles;
        for (int first = 2; first <= level - 3; first += 2) {
            for (int second = first;
                 second <= level - 3; second += 2) {
                for (int third = second;
                     third <= level - 3; third += 2) {
                    const Labels triple{first, second, third};
                    const int copies_two =
                        static_cast<int>(std::count(
                            triple.begin(),
                            triple.end(),
                            2
                        ));
                    if (copies_two <= 2
                        && invariant(triple, level) == 1) {
                        std::array<std::int64_t, 5> profile{};
                        for (std::size_t output = 0U;
                             output < profile.size(); ++output) {
                            profile[output] = multiplicity(
                                triple,
                                static_cast<int>(2U * output),
                                level
                            );
                        }
                        selected_profiles[
                            static_cast<std::size_t>(copies_two)
                        ].push_back(profile);
                    }
                    Labels block = triple;
                    block.push_back(level - 1);
                    if (copies_two <= 2
                        && invariant(block, level) == 2) {
                        std::array<std::int64_t, 5> profile{};
                        for (std::size_t output = 0U;
                             output < profile.size(); ++output) {
                            profile[output] = multiplicity(
                                block,
                                static_cast<int>(2U * output),
                                level
                            );
                        }
                        complement_profiles[
                            static_cast<std::size_t>(copies_two)
                        ].push_back(profile);
                    }
                }
            }
        }
        for (auto& profiles : selected_profiles) {
            std::sort(profiles.begin(), profiles.end());
            profiles.erase(
                std::unique(profiles.begin(), profiles.end()),
                profiles.end()
            );
        }
        for (auto& profiles : complement_profiles) {
            std::sort(profiles.begin(), profiles.end());
            profiles.erase(
                std::unique(profiles.begin(), profiles.end()),
                profiles.end()
            );
        }
        for (std::size_t selected_twos = 0U;
             selected_twos < selected_profiles.size();
             ++selected_twos) {
            for (std::size_t complement_twos = 0U;
                 complement_twos < complement_profiles.size();
                 ++complement_twos) {
                if (selected_twos + complement_twos > 2U) {
                    continue;
                }
                for (const auto& selected :
                     selected_profiles[selected_twos]) {
                    for (const auto& complement :
                         complement_profiles[complement_twos]) {
                        std::int64_t local = 0;
                        for (std::size_t output = 0U;
                             output < selected.size(); ++output) {
                            local += selected[output]
                                * complement[output];
                        }
                        factorized_minimum[selected_twos][
                            complement_twos
                        ] = std::min(
                            factorized_minimum[selected_twos][
                                complement_twos
                            ],
                            local
                        );
                    }
                }
            }
        }
    }
    for (int level = 7; level <= maximum_level; level += 2) {
        for (int middle = 4; middle <= level - 3; middle += 2) {
            for (int delta = -1; delta <= 1; ++delta) {
                const int last = middle + 2 * delta;
                if (last < 4 || last > level - 3) {
                    continue;
                }
                std::array<std::int64_t, 5> profile{};
                for (std::size_t output = 0U;
                     output < profile.size(); ++output) {
                    profile[output] = multiplicity(
                        Labels{2, middle, last},
                        static_cast<int>(2U * output),
                        level
                    );
                }
                const std::int64_t weighted =
                    2 * profile[0]
                    + 5 * profile[1]
                    + 6 * profile[2]
                    + 4 * profile[3];
                const std::size_t index =
                    static_cast<std::size_t>(delta + 1);
                if (weighted < one_two_minimum[index]) {
                    one_two_minimum[index] = weighted;
                    one_two_profile[index] = profile;
                }
            }
        }
    }
    std::cout
        << "NO_TWO_SELECTED_PROFILE cases=" << cases
        << " minimum=" << minimum
        << " level=" << minimum_level
        << " label=" << minimum_label
        << " interval=[" << minimum_lower << ',' << minimum_upper
        << "] profile=[";
    for (std::size_t index = 0U;
         index < minimum_profile.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << minimum_profile[index];
    }
    std::cout << "]\n";
    for (std::size_t complement_twos = 0U;
         complement_twos < no_two_weight_minimum.size();
         ++complement_twos) {
        std::cout
            << "NO_TWO_SELECTED_PROFILE complement_twos="
            << complement_twos
            << " minimum="
            << no_two_weight_minimum[complement_twos]
            << " profile=[";
        for (std::size_t output = 0U;
             output
                 < no_two_weight_profile[complement_twos].size();
             ++output) {
            if (output != 0U) {
                std::cout << ',';
            }
            std::cout
                << no_two_weight_profile[complement_twos][output];
        }
        std::cout << "]\n";
    }
    for (int delta = -1; delta <= 1; ++delta) {
        const std::size_t index =
            static_cast<std::size_t>(delta + 1);
        std::cout
            << "NO_TWO_SELECTED_PROFILE one_two_delta=" << delta
            << " minimum=" << one_two_minimum[index]
            << " profile=[";
        for (std::size_t output = 0U;
             output < one_two_profile[index].size(); ++output) {
            if (output != 0U) {
                std::cout << ',';
            }
            std::cout << one_two_profile[index][output];
        }
        std::cout << "]\n";
    }
    for (std::size_t selected_twos = 0U;
         selected_twos < factorized_minimum.size();
         ++selected_twos) {
        for (std::size_t complement_twos = 0U;
             complement_twos
                 < factorized_minimum[selected_twos].size();
             ++complement_twos) {
            if (selected_twos + complement_twos > 2U) {
                continue;
            }
            std::cout
                << "NO_TWO_SELECTED_PROFILE factorized_twos="
                << selected_twos << ',' << complement_twos
                << " minimum="
                << factorized_minimum[selected_twos][
                    complement_twos
                ]
                << "\n";
        }
    }
    if (minimum < 41) {
        throw std::runtime_error(
            "counterexample to no-two selected profile"
        );
    }
    for (const std::int64_t value : one_two_minimum) {
        if (value < 42) {
            throw std::runtime_error(
                "counterexample to one-two selected profile"
            );
        }
    }
    std::cout << "NO_TWO_SELECTED_PROFILE PASS\n";
}

void scan_plus_neighbor_rank_two(
    int maximum_level
) {
    SmallRankWitness direct_margin;
    SmallRankWitness supply_margin;
    SmallRankWitness completed_direct_margin;
    SmallRankWitness completed_supply_margin;
    SmallRankWitness no_triple_supply_margin;
    SmallRankWitness no_triple_local_margin;
    std::array<std::array<SmallRankWitness, 3>, 3>
        no_triple_local_by_twos{};
    SmallRankWitness low_rank_triple_completed_supply_margin;
    SmallRankWitness nontrivial_supply_margin;
    SmallRankWitness all_two_direct_margin;
    std::uint64_t cases = 0;
    std::uint64_t supply_below_forty = 0;
    std::uint64_t completed_supply_below_forty = 0;
    std::uint64_t nontrivial_supply_below_forty = 0;

    for (int level = 5; level <= maximum_level; level += 2) {
        Labels even;
        for (int label = 2; label < level - 1; label += 2) {
            even.push_back(label);
        }
        const std::vector<Labels> minus_lists =
            sorted_lists(even, 6);
        for (const Labels& minus : minus_lists) {
            Labels labels = minus;
            labels.push_back(level - 1);

            std::int64_t pair = 0;
            std::int64_t completed_pair = 0;
            const int copies_two = static_cast<int>(std::count(
                labels.begin(),
                labels.begin() + 6,
                2
            ));
            int triple_two_complement_rank = -1;
            if (copies_two >= 3) {
                Labels triple_two_rest;
                for (std::size_t index = 3U;
                     index < labels.size(); ++index) {
                    triple_two_rest.push_back(labels[index]);
                }
                triple_two_complement_rank =
                    static_cast<int>(invariant(
                        triple_two_rest,
                        level
                    ));
            }
            for (std::size_t first = 0U; first < 6U; ++first) {
                for (std::size_t second = first + 1U;
                     second < 6U; ++second) {
                    if (labels[first] != labels[second]) {
                        continue;
                    }
                    Labels rest;
                    for (std::size_t index = 0U;
                         index < labels.size(); ++index) {
                        if (index != first && index != second) {
                            rest.push_back(labels[index]);
                        }
                    }
                    const std::int64_t retained = local_channels(
                        Labels{rest[0], rest[1]},
                        Labels{rest[2], rest[3], rest[4]},
                        level
                    );
                    pair += retained;
                    completed_pair += (
                        labels[first] == 2 && copies_two >= 3
                    ) ? invariant(rest, level) : retained;
                }
            }

            std::int64_t positive = 0;
            std::int64_t negative = 0;
            std::int64_t maximum = 0;
            std::vector<std::pair<unsigned int, std::int64_t>>
                negative_terms;
            for (unsigned int mask = 0U;
                 mask < (1U << 7U); ++mask) {
                if (popcount(mask) != 3) {
                    continue;
                }
                const std::int64_t term =
                    invariant(subset(labels, mask, true), level)
                    * invariant(subset(labels, mask, false), level);
                int minus_in_cut = 0;
                for (int position = 0; position < 6; ++position) {
                    minus_in_cut += static_cast<int>(
                        (mask >> position) & 1U
                    );
                }
                if ((minus_in_cut & 1) == 0) {
                    positive += term;
                } else {
                    negative += term;
                    maximum = std::max(maximum, term);
                    if (term > 0) {
                        negative_terms.emplace_back(mask, term);
                    }
                }
            }
            if (maximum != 2) {
                continue;
            }
            for (const auto& [mask, term] : negative_terms) {
                if (term != maximum) {
                    continue;
                }
                ++cases;
                const Labels selected = subset(labels, mask, true);
                const Labels complement =
                    subset(labels, mask, false);
                const std::int64_t local =
                    local_channels(selected, complement, level);
                std::array<std::int64_t, 5> selected_profile{};
                std::array<std::int64_t, 5> complement_profile{};
                for (std::size_t output = 0U;
                     output < selected_profile.size(); ++output) {
                    selected_profile[output] = multiplicity(
                        selected, static_cast<int>(2U * output), level
                    );
                    complement_profile[output] = multiplicity(
                        complement,
                        static_cast<int>(2U * output),
                        level
                    );
                }
                int rank_one = 0;
                int rank_two = 0;
                for (const auto& [other_mask, other_term]
                     : negative_terms) {
                    static_cast<void>(other_mask);
                    rank_one += static_cast<int>(other_term == 1);
                    rank_two += static_cast<int>(other_term == 2);
                }
                const bool all_two =
                    selected[0] == 2
                    && selected[1] == 2
                    && selected[2] == 2;
                if (local + pair < 40) {
                    ++supply_below_forty;
                    nontrivial_supply_below_forty +=
                        static_cast<std::uint64_t>(!all_two);
                    if (supply_below_forty <= 64U) {
                        std::cout
                            << "PLUS_NEIGHBOR_D2 low_supply"
                            << " level=" << level
                            << " mask=" << mask
                            << " all_two=" << (all_two ? 1 : 0)
                            << " labels=[";
                        for (std::size_t index = 0U;
                             index < labels.size(); ++index) {
                            if (index != 0U) {
                                std::cout << ',';
                            }
                            std::cout << labels[index];
                        }
                        std::cout
                            << "] L=" << local
                            << " R=" << pair
                            << " U=" << positive
                            << " T=" << negative
                            << '\n';
                    }
                }
                completed_supply_below_forty +=
                    static_cast<std::uint64_t>(
                        local + completed_pair < 40
                    );
                consider_small_rank(
                    direct_margin,
                    local + pair + positive - negative,
                    level,
                    6,
                    mask,
                    6,
                    labels,
                    local,
                    pair,
                    positive,
                    negative,
                    selected_profile,
                    complement_profile,
                    rank_one,
                    rank_two
                );
                consider_small_rank(
                    completed_direct_margin,
                    local + completed_pair + positive - negative,
                    level,
                    6,
                    mask,
                    6,
                    labels,
                    local,
                    completed_pair,
                    positive,
                    negative,
                    selected_profile,
                    complement_profile,
                    rank_one,
                    rank_two
                );
                consider_small_rank(
                    supply_margin,
                    local + pair,
                    level,
                    6,
                    mask,
                    6,
                    labels,
                    local,
                    pair,
                    positive,
                    negative,
                    selected_profile,
                    complement_profile,
                    rank_one,
                    rank_two
                );
                consider_small_rank(
                    completed_supply_margin,
                    local + completed_pair,
                    level,
                    6,
                    mask,
                    6,
                    labels,
                    local,
                    completed_pair,
                    positive,
                    negative,
                    selected_profile,
                    complement_profile,
                    rank_one,
                    rank_two
                );
                if (copies_two <= 2) {
                    const int selected_twos =
                        static_cast<int>(std::count(
                            selected.begin(),
                            selected.end(),
                            2
                        ));
                    const int complement_twos =
                        copies_two - selected_twos;
                    consider_small_rank(
                        no_triple_supply_margin,
                        local + pair,
                        level,
                        6,
                        mask,
                        6,
                        labels,
                        local,
                        pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        rank_one,
                        rank_two
                    );
                    consider_small_rank(
                        no_triple_local_by_twos[
                            static_cast<std::size_t>(selected_twos)
                        ][static_cast<std::size_t>(complement_twos)],
                        local,
                        level,
                        6,
                        mask,
                        6,
                        labels,
                        local,
                        pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        rank_one,
                        rank_two
                    );
                    consider_small_rank(
                        no_triple_local_margin,
                        local,
                        level,
                        6,
                        mask,
                        6,
                        labels,
                        local,
                        pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        rank_one,
                        rank_two
                    );
                } else if (triple_two_complement_rank <= 1) {
                    consider_small_rank(
                        low_rank_triple_completed_supply_margin,
                        local + completed_pair,
                        level,
                        6,
                        mask,
                        6,
                        labels,
                        local,
                        completed_pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        rank_one,
                        rank_two
                    );
                }
                if (!all_two) {
                    consider_small_rank(
                        nontrivial_supply_margin,
                        local + pair,
                        level,
                        6,
                        mask,
                        6,
                        labels,
                        local,
                        pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        rank_one,
                        rank_two
                    );
                } else {
                    consider_small_rank(
                        all_two_direct_margin,
                        local + pair + positive - negative,
                        level,
                        6,
                        mask,
                        6,
                        labels,
                        local,
                        pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        rank_one,
                        rank_two
                    );
                }
            }
        }
        std::cout
            << "PLUS_NEIGHBOR_D2 progress_level=" << level
            << " cases=" << cases
            << " supply_below_40=" << supply_below_forty
            << " completed_supply_below_40="
            << completed_supply_below_forty
            << " nontrivial_below_40="
            << nontrivial_supply_below_forty
            << '\n' << std::flush;
    }

    std::cout
        << "PLUS_NEIGHBOR_D2 cases=" << cases
        << " supply_below_40=" << supply_below_forty
        << " completed_supply_below_40="
        << completed_supply_below_forty
        << " nontrivial_below_40="
        << nontrivial_supply_below_forty;
    print_small_rank_witness("direct_margin", direct_margin);
    print_small_rank_witness("supply", supply_margin);
    print_small_rank_witness(
        "completed_direct_margin",
        completed_direct_margin
    );
    print_small_rank_witness(
        "completed_supply",
        completed_supply_margin
    );
    print_small_rank_witness(
        "no_triple_supply",
        no_triple_supply_margin
    );
    print_small_rank_witness(
        "no_triple_local",
        no_triple_local_margin
    );
    for (std::size_t selected_twos = 0U;
         selected_twos < no_triple_local_by_twos.size();
         ++selected_twos) {
        for (std::size_t complement_twos = 0U;
             complement_twos
                 < no_triple_local_by_twos[selected_twos].size();
             ++complement_twos) {
            if (selected_twos + complement_twos > 2U) {
                continue;
            }
            const std::string name =
                "local_twos_"
                + std::to_string(selected_twos)
                + "_"
                + std::to_string(complement_twos);
            print_small_rank_witness(
                name.c_str(),
                no_triple_local_by_twos[selected_twos][
                    complement_twos
                ]
            );
        }
    }
    print_small_rank_witness(
        "low_rank_triple_completed_supply",
        low_rank_triple_completed_supply_margin
    );
    print_small_rank_witness(
        "nontrivial_supply", nontrivial_supply_margin
    );
    print_small_rank_witness(
        "all_two_direct_margin", all_two_direct_margin
    );
    std::cout << '\n';
    if (
        !completed_supply_margin.initialized
        || completed_supply_margin.value < 40
        || completed_supply_below_forty != 0U
    ) {
        throw std::runtime_error(
            "bounded counterexample to completed plus-neighbor supply"
        );
    }
    std::cout << "PLUS_NEIGHBOR_D2 PASS\n";
}

void consider_small_rank(
    SmallRankWitness& witness,
    std::int64_t value,
    int level,
    int orbit_index,
    unsigned int selected_mask,
    int endpoint_position,
    const Labels& labels,
    std::int64_t local,
    std::int64_t pair,
    std::int64_t positive,
    std::int64_t negative,
    const std::array<std::int64_t, 5>& selected_profile,
    const std::array<std::int64_t, 5>& complement_profile,
    int active_rank_one,
    int active_rank_two
) {
    if (!witness.initialized || value < witness.value) {
        witness = SmallRankWitness{
            true,
            value,
            level,
            orbit_index,
            selected_mask,
            endpoint_position,
            labels,
            local,
            pair,
            positive,
            negative,
            selected_profile,
            complement_profile,
            active_rank_one,
            active_rank_two
        };
    }
}

int endpoint_kind(int label, int level) {
    if (label == 1) {
        return 0;
    }
    if (label == level - 1) {
        return 1;
    }
    if (label == level) {
        return 2;
    }
    return -1;
}

void scan_small_rank_level(
    int orbit_index,
    int level,
    std::array<
        std::array<
            std::array<
                std::array<SmallRankSummary, 2>,
                3
            >,
            2
        >,
        8
    >& summaries,
    SmallRankWitness& lemma_b16,
    std::uint64_t& lemma_b16_cases,
    std::array<std::array<SmallRankSummary, 4>, 7>& orbit6_neighbor
) {
    const Orbit& orbit =
        orbits[static_cast<std::size_t>(orbit_index)];
    Labels even;
    Labels odd;
    for (int label = 1; label <= level; ++label) {
        ((label & 1) == 0 ? even : odd).push_back(label);
    }

    const int even_minus =
        orbit.minus_count - orbit.odd_minus_count;
    const int odd_plus =
        orbit.odd_count - orbit.odd_minus_count;
    const int plus_count = 7 - orbit.minus_count;
    const int even_plus = plus_count - odd_plus;
    const std::vector<Labels> minus_even =
        sorted_lists(even, even_minus);
    const std::vector<Labels> minus_odd =
        sorted_lists(odd, orbit.odd_minus_count);
    const std::vector<Labels> plus_even =
        sorted_lists(even, even_plus);
    const std::vector<Labels> plus_odd =
        sorted_lists(odd, odd_plus);

    std::vector<Labels> minus_lists;
    std::vector<Labels> plus_lists;
    for (const Labels& first : minus_even) {
        for (const Labels& second : minus_odd) {
            minus_lists.push_back(merge_labels(first, second));
        }
    }
    for (const Labels& first : plus_even) {
        for (const Labels& second : plus_odd) {
            plus_lists.push_back(merge_labels(first, second));
        }
    }

    constexpr unsigned int full_mask = (1U << 7U) - 1U;
    for (const Labels& minus : minus_lists) {
        for (const Labels& plus : plus_lists) {
            if (!disjoint_support(minus, plus)) {
                continue;
            }
            Labels labels = minus;
            labels.insert(labels.end(), plus.begin(), plus.end());

            std::int64_t pair = 0;
            for (std::size_t first = 0U;
                 first < labels.size(); ++first) {
                for (std::size_t second = first + 1U;
                     second < labels.size(); ++second) {
                    if (labels[first] != labels[second]) {
                        continue;
                    }
                    Labels rest;
                    for (std::size_t index = 0U;
                         index < labels.size(); ++index) {
                        if (index != first && index != second) {
                            rest.push_back(labels[index]);
                        }
                    }
                    pair += local_channels(
                        Labels{rest[0], rest[1]},
                        Labels{rest[2], rest[3], rest[4]},
                        level
                    );
                }
            }

            std::int64_t positive = 0;
            std::int64_t negative = 0;
            std::int64_t maximum = 0;
            int active_rank_one = 0;
            int active_rank_two = 0;
            std::vector<std::pair<unsigned int, std::int64_t>>
                negative_terms;
            for (unsigned int mask = 0U;
                 mask <= full_mask; ++mask) {
                if (popcount(mask) != 3) {
                    continue;
                }
                const std::int64_t term =
                    invariant(subset(labels, mask, true), level)
                    * invariant(subset(labels, mask, false), level);
                int minus_in_cut = 0;
                for (int index = 0;
                     index < orbit.minus_count; ++index) {
                    minus_in_cut += static_cast<int>(
                        (mask >> index) & 1U
                    );
                }
                if ((minus_in_cut & 1) == 0) {
                    positive += term;
                    continue;
                }
                negative += term;
                if (term > 0) {
                    negative_terms.emplace_back(mask, term);
                    active_rank_one += static_cast<int>(term == 1);
                    active_rank_two += static_cast<int>(term == 2);
                    maximum = std::max(maximum, term);
                }
            }
            if (maximum < 1 || maximum > 2) {
                continue;
            }

            for (const auto& [mask, term] : negative_terms) {
                if (term != maximum) {
                    continue;
                }
                const Labels selected = subset(labels, mask, true);
                const Labels complement = subset(labels, mask, false);
                const std::int64_t local =
                    local_channels(selected, complement, level);
                std::array<std::int64_t, 5> selected_profile{};
                std::array<std::int64_t, 5> complement_profile{};
                for (std::size_t output = 0U;
                     output < selected_profile.size(); ++output) {
                    selected_profile[output] = multiplicity(
                        selected, static_cast<int>(2U * output), level
                    );
                    complement_profile[output] = multiplicity(
                        complement,
                        static_cast<int>(2U * output),
                        level
                    );
                }
                int shallow_count = 0;
                int top_position = -1;
                for (int position = 0; position < 7; ++position) {
                    const int label =
                        labels[static_cast<std::size_t>(position)];
                    shallow_count += static_cast<int>(
                        std::min(label, level - label) + 1 <= 2
                    );
                    if (label == level) {
                        top_position = position;
                    }
                }
                const int selected_minus_count = [&]() {
                    int result = 0;
                    for (int position = 0;
                         position < orbit.minus_count;
                         ++position) {
                        result += static_cast<int>(
                            ((mask >> position) & 1U) != 0U
                        );
                    }
                    return result;
                }();
                if (maximum == 1) {
                    std::array<std::int64_t, 3>
                        bad_by_intersection{};
                    std::int64_t bad_switches = 0;
                    for (const auto& [
                             other_mask,
                             other_term
                         ] : negative_terms) {
                        if (other_mask == mask) {
                            continue;
                        }
                        const unsigned int switched_mask =
                            mask ^ other_mask;
                        const std::int64_t switched_term =
                            invariant(
                                subset(
                                    labels,
                                    switched_mask,
                                    true
                                ),
                                level
                            )
                            * invariant(
                                subset(
                                    labels,
                                    switched_mask,
                                    false
                                ),
                                level
                            );
                        if (other_term != 1
                            || switched_term != 0) {
                            continue;
                        }
                        ++bad_switches;
                        const int intersection = popcount(
                            mask & other_mask
                        );
                        if (intersection < 0
                            || intersection > 2) {
                            throw std::runtime_error(
                                "invalid distinct triple "
                                "intersection"
                            );
                        }
                        ++bad_by_intersection[
                            static_cast<std::size_t>(
                                intersection
                            )
                        ];
                    }
                    ++rank_one_switch_cases;
                    consider_small_rank(
                        rank_one_switch_margin,
                        local - 1 - bad_switches,
                        level,
                        orbit_index,
                        mask,
                        top_position,
                        labels,
                        local,
                        pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        active_rank_one,
                        active_rank_two
                    );
                    consider_small_rank(
                        rank_one_switch_bad_maximum,
                        -bad_switches,
                        level,
                        orbit_index,
                        mask,
                        top_position,
                        labels,
                        local,
                        pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        active_rank_one,
                        active_rank_two
                    );
                    for (std::size_t intersection = 0U;
                         intersection
                             < bad_by_intersection.size();
                         ++intersection) {
                        consider_small_rank(
                            rank_one_switch_component_margin[
                                intersection
                            ],
                            local - 1
                                - bad_by_intersection[
                                    intersection
                                ],
                            level,
                            orbit_index,
                            mask,
                            top_position,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                        consider_small_rank(
                            rank_one_switch_component_bad_maximum[
                                intersection
                            ],
                            -bad_by_intersection[
                                intersection
                            ],
                            level,
                            orbit_index,
                            mask,
                            top_position,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                    }
                    for (int position = 0;
                         position < 7;
                         ++position) {
                        const int kind = endpoint_kind(
                            labels[
                                static_cast<std::size_t>(position)
                            ],
                            level
                        );
                        if (kind < 0) {
                            continue;
                        }
                        const std::size_t side =
                            ((mask >> position) & 1U) != 0U
                                ? 0U
                                : 1U;
                        consider_small_rank(
                            rank_one_switch_endpoint_margin[
                                static_cast<std::size_t>(kind)
                            ][side],
                            local - 1 - bad_switches,
                            level,
                            orbit_index,
                            mask,
                            position,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                    }
                }
                const bool top_is_selected =
                    top_position >= 0
                    && ((mask >> top_position) & 1U) != 0U;
                const bool lemma_b16_packet =
                    maximum == 2
                    && shallow_count == 1
                    && top_is_selected
                    && (
                        (
                            orbit_index == 3
                            && (
                                (
                                    (level & 1) == 1
                                    && top_position
                                        >= orbit.minus_count
                                    && selected_minus_count == 1
                                )
                                || (
                                    (level & 1) == 0
                                    && top_position
                                        < orbit.minus_count
                                    && selected_minus_count == 3
                                )
                            )
                        )
                        || (
                            orbit_index == 6
                            && (level & 1) == 0
                        )
                        || (
                            orbit_index == 7
                            && (level & 1) == 0
                        )
                    );
                if (lemma_b16_packet) {
                    ++lemma_b16_cases;
                    consider_small_rank(
                        lemma_b16,
                        local + pair - negative,
                        level,
                        orbit_index,
                        mask,
                        top_position,
                        labels,
                        local,
                        pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        active_rank_one,
                        active_rank_two
                    );
                }
                if (orbit_index == 6
                    && (level & 1) == 1
                    && maximum == 2) {
                    int total_neighbor = 0;
                    int selected_neighbor = 0;
                    for (int position = 0;
                         position < 7; ++position) {
                        if (labels[
                                static_cast<std::size_t>(position)
                            ] != level - 1) {
                            continue;
                        }
                        ++total_neighbor;
                        selected_neighbor += static_cast<int>(
                            ((mask >> position) & 1U) != 0U
                        );
                    }
                    if (total_neighbor > 0) {
                        std::int64_t negative_with_neighbor = 0;
                        std::int64_t negative_without_neighbor = 0;
                        for (const auto& [negative_mask, weight]
                             : negative_terms) {
                            bool contains_neighbor = false;
                            for (int position = 0;
                                 position < 6; ++position) {
                                contains_neighbor =
                                    contains_neighbor
                                    || (
                                        labels[static_cast<std::size_t>(
                                            position
                                        )] == level - 1
                                        && (
                                            (
                                                negative_mask
                                                >> position
                                            ) & 1U
                                        ) != 0U
                                    );
                            }
                            (contains_neighbor
                                 ? negative_with_neighbor
                                 : negative_without_neighbor) += weight;
                        }
                        SmallRankSummary& neighbor_summary =
                            orbit6_neighbor[
                                static_cast<std::size_t>(total_neighbor)
                            ][static_cast<std::size_t>(
                                selected_neighbor
                            )];
                        ++neighbor_summary.cases;
                        consider_small_rank(
                            neighbor_summary.direct,
                            local + pair - negative,
                            level,
                            orbit_index,
                            mask,
                            0,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                        const std::size_t total_index =
                            static_cast<std::size_t>(total_neighbor);
                        const std::size_t selected_index =
                            static_cast<std::size_t>(
                                selected_neighbor
                            );
                        consider_small_rank(
                            neighbor_local_with_payment[total_index]
                                [selected_index],
                            local - negative_with_neighbor,
                            level,
                            orbit_index,
                            mask,
                            0,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                        consider_small_rank(
                            neighbor_pair_without_payment[total_index]
                                [selected_index],
                            pair - negative_without_neighbor,
                            level,
                            orbit_index,
                            mask,
                            0,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                        consider_small_rank(
                            neighbor_local_without_payment[total_index]
                                [selected_index],
                            local - negative_without_neighbor,
                            level,
                            orbit_index,
                            mask,
                            0,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                        consider_small_rank(
                            neighbor_pair_with_payment[total_index]
                                [selected_index],
                            pair - negative_with_neighbor,
                            level,
                            orbit_index,
                            mask,
                            0,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                        consider_small_rank(
                            neighbor_local_total_payment[total_index]
                                [selected_index],
                            local - negative,
                            level,
                            orbit_index,
                            mask,
                            0,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                        if (total_neighbor == 1
                            && selected_neighbor == 0) {
                            int neighbor_position = -1;
                            for (int position = 0;
                                 position < 7; ++position) {
                                if (labels[
                                        static_cast<std::size_t>(
                                            position
                                        )
                                    ] == level - 1) {
                                    neighbor_position = position;
                                }
                            }
                            if (neighbor_position < 0) {
                                throw std::runtime_error(
                                    "missing sole neighbor position"
                                );
                            }
                            const Labels selected_labels =
                                subset(labels, mask, true);
                            int category = 0;
                            for (std::size_t first = 0U;
                                 first < selected_labels.size();
                                 ++first) {
                                for (std::size_t second = first + 1U;
                                     second < selected_labels.size();
                                     ++second) {
                                    if (selected_labels[first]
                                        != selected_labels[second]) {
                                        continue;
                                    }
                                    int third = 0;
                                    for (std::size_t index = 0U;
                                         index < selected_labels.size();
                                         ++index) {
                                        if (index != first
                                            && index != second) {
                                            third =
                                                selected_labels[index];
                                        }
                                    }
                                    category = third == 2
                                        ? 1
                                        : (third == 4 ? 2 : 3);
                                }
                            }
                            SmallRankSummary& sign_summary =
                                neighbor_zero_selected_by_sign[
                                    static_cast<std::size_t>(
                                        neighbor_position < 6 ? 0 : 1
                                    )
                                ][static_cast<std::size_t>(category)];
                            const std::size_t sign_index =
                                static_cast<std::size_t>(
                                    neighbor_position < 6 ? 0 : 1
                                );
                            const std::size_t category_index =
                                static_cast<std::size_t>(category);
                            ++sign_summary.cases;
                            consider_small_rank(
                                sign_summary.direct,
                                local + pair - negative,
                                level,
                                orbit_index,
                                mask,
                                neighbor_position,
                                labels,
                                local,
                                pair,
                                positive,
                                negative,
                                selected_profile,
                                complement_profile,
                                active_rank_one,
                                active_rank_two
                            );
                            consider_small_rank(
                                neighbor_zero_selected_supply[sign_index]
                                    [category_index],
                                local + pair,
                                level,
                                orbit_index,
                                mask,
                                neighbor_position,
                                labels,
                                local,
                                pair,
                                positive,
                                negative,
                                selected_profile,
                                complement_profile,
                                active_rank_one,
                                active_rank_two
                            );
                            consider_small_rank(
                                neighbor_zero_selected_local[sign_index]
                                    [category_index],
                                local,
                                level,
                                orbit_index,
                                mask,
                                neighbor_position,
                                labels,
                                local,
                                pair,
                                positive,
                                negative,
                                selected_profile,
                                complement_profile,
                                active_rank_one,
                                active_rank_two
                            );
                            consider_small_rank(
                                neighbor_zero_selected_pair[sign_index]
                                    [category_index],
                                pair,
                                level,
                                orbit_index,
                                mask,
                                neighbor_position,
                                labels,
                                local,
                                pair,
                                positive,
                                negative,
                                selected_profile,
                                complement_profile,
                                active_rank_one,
                                active_rank_two
                            );
                            consider_small_rank(
                                sign_summary.ceiling,
                                38 + pair - negative,
                                level,
                                orbit_index,
                                mask,
                                neighbor_position,
                                labels,
                                local,
                                pair,
                                positive,
                                negative,
                                selected_profile,
                                complement_profile,
                                active_rank_one,
                                active_rank_two
                            );
                            consider_small_rank(
                                sign_summary.local,
                                40 - negative,
                                level,
                                orbit_index,
                                mask,
                                neighbor_position,
                                labels,
                                local,
                                pair,
                                positive,
                                negative,
                                selected_profile,
                                complement_profile,
                                active_rank_one,
                                active_rank_two
                            );
                        }
                    }
                    if (total_neighbor == 2
                        && selected_neighbor == 2) {
                        std::vector<int> q_positions;
                        for (int position = 0; position < 6;
                             ++position) {
                            if (labels[
                                    static_cast<std::size_t>(position)
                                ] != level - 1) {
                                q_positions.push_back(position);
                            }
                        }
                        if (q_positions.size() != 4U) {
                            throw std::runtime_error(
                                "invalid double-neighbor Q block"
                            );
                        }
                        int copies_two = 0;
                        for (const int position : q_positions) {
                            copies_two += static_cast<int>(
                                labels[
                                    static_cast<std::size_t>(position)
                                ] == 2
                            );
                        }
                        std::int64_t edge_weight = 0;
                        for (std::size_t first = 0U;
                             first < q_positions.size(); ++first) {
                            for (std::size_t second = first + 1U;
                                 second < q_positions.size();
                                 ++second) {
                                const int first_position =
                                    q_positions[first];
                                const int second_position =
                                    q_positions[second];
                                const int sum =
                                    labels[static_cast<std::size_t>(
                                        first_position
                                    )]
                                    + labels[static_cast<std::size_t>(
                                        second_position
                                    )];
                                if (sum != level - 1
                                    && sum != level + 1) {
                                    continue;
                                }
                                Labels rest{
                                    level - 1,
                                    labels[6]
                                };
                                for (const int position : q_positions) {
                                    if (position != first_position
                                        && position
                                            != second_position) {
                                        rest.push_back(labels[
                                            static_cast<std::size_t>(
                                                position
                                            )
                                        ]);
                                    }
                                }
                                edge_weight += invariant(rest, level);
                            }
                        }
                        std::int64_t no_neighbor = 0;
                        for (std::size_t omitted = 0U;
                             omitted < q_positions.size(); ++omitted) {
                            Labels cut;
                            Labels rest{
                                level - 1,
                                level - 1,
                                labels[6]
                            };
                            for (std::size_t index = 0U;
                                 index < q_positions.size(); ++index) {
                                const int label = labels[
                                    static_cast<std::size_t>(
                                        q_positions[index]
                                    )
                                ];
                                (index == omitted ? rest : cut)
                                    .push_back(label);
                            }
                            no_neighbor += invariant(cut, level)
                                * invariant(rest, level);
                        }
                        const std::size_t n_index =
                            static_cast<std::size_t>(copies_two - 1);
                        consider_small_rank(
                            double_neighbor_margin[n_index],
                            local - negative,
                            level,
                            orbit_index,
                            mask,
                            0,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                        consider_small_rank(
                            double_neighbor_w_payment[n_index],
                            complement_profile[1] - edge_weight,
                            level,
                            orbit_index,
                            mask,
                            0,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                        consider_small_rank(
                            double_neighbor_tail_payment[n_index],
                            complement_profile[2]
                                - (
                                    2 * copies_two
                                    + no_neighbor
                                    - 2
                                ),
                            level,
                            orbit_index,
                            mask,
                            0,
                            labels,
                            local,
                            pair,
                            positive,
                            negative,
                            selected_profile,
                            complement_profile,
                            active_rank_one,
                            active_rank_two
                        );
                    }
                }
                for (int position = 0; position < 7; ++position) {
                    const int kind = endpoint_kind(
                        labels[static_cast<std::size_t>(position)],
                        level
                    );
                    if (kind < 0) {
                        continue;
                    }
                    const int side =
                        ((mask >> position) & 1U) != 0U ? 0 : 1;
                    SmallRankSummary& summary =
                        summaries[
                            static_cast<std::size_t>(orbit_index)
                        ][static_cast<std::size_t>(maximum - 1)]
                         [static_cast<std::size_t>(kind)]
                         [static_cast<std::size_t>(side)];
                    ++summary.cases;
                    const std::int64_t direct =
                        local + pair + positive - negative;
                    consider_small_rank(
                        summary.direct,
                        direct,
                        level,
                        orbit_index,
                        mask,
                        position,
                        labels,
                        local,
                        pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        active_rank_one,
                        active_rank_two
                    );
                    consider_small_rank(
                        summary.ceiling,
                        local + pair + positive - 20 * maximum,
                        level,
                        orbit_index,
                        mask,
                        position,
                        labels,
                        local,
                        pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        active_rank_one,
                        active_rank_two
                    );
                    consider_small_rank(
                        summary.local,
                        local,
                        level,
                        orbit_index,
                        mask,
                        position,
                        labels,
                        local,
                        pair,
                        positive,
                        negative,
                        selected_profile,
                        complement_profile,
                        active_rank_one,
                        active_rank_two
                    );
                }
            }
        }
    }
}

void print_small_rank_witness(
    const char* name,
    const SmallRankWitness& witness
) {
    std::cout << ' ' << name << '=';
    if (!witness.initialized) {
        std::cout << "none";
        return;
    }
    std::cout << witness.value
              << "@k" << witness.level
              << ":orbit" << witness.orbit_index
              << ":mask" << witness.selected_mask
              << ":position" << witness.endpoint_position
              << ":labels";
    print_labels(witness.labels);
    std::cout << ":(L,P,U,T)=("
              << witness.local << ','
              << witness.pair << ','
              << witness.positive << ','
              << witness.negative << ")"
              << ":ranks=("
              << witness.active_rank_one << ','
              << witness.active_rank_two << ")"
              << ":A[";
    for (std::size_t index = 0U;
         index < witness.selected_profile.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << witness.selected_profile[index];
    }
    std::cout << "]:B[";
    for (std::size_t index = 0U;
         index < witness.complement_profile.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << witness.complement_profile[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (
            argc == 12
            && std::string(argv[1]) == "--switch-word"
        ) {
            const int level = parse_positive(argv[2]);
            const int minus_count = parse_positive(argv[3]);
            const unsigned int selected_mask =
                static_cast<unsigned int>(parse_positive(argv[4]));
            if (minus_count > 6
                || popcount(selected_mask) != 3) {
                throw std::runtime_error(
                    "invalid switch-word sign count or mask"
                );
            }
            Labels labels;
            for (int index = 0; index < 7; ++index) {
                labels.push_back(parse_positive(argv[index + 5]));
            }
            const std::int64_t selected_rank =
                invariant(
                    subset(labels, selected_mask, true),
                    level
                )
                * invariant(
                    subset(labels, selected_mask, false),
                    level
                );
            const std::int64_t local = local_channels(
                subset(labels, selected_mask, true),
                subset(labels, selected_mask, false),
                level
            );
            std::int64_t bad = 0;
            constexpr unsigned int full_mask =
                (1U << 7U) - 1U;
            std::cout
                << "switch_word selected_rank=" << selected_rank
                << " local=" << local << '\n';
            for (unsigned int mask = 0U;
                 mask <= full_mask; ++mask) {
                if (popcount(mask) != 3
                    || mask == selected_mask) {
                    continue;
                }
                int minus_in_cut = 0;
                for (int position = 0;
                     position < minus_count;
                     ++position) {
                    minus_in_cut += static_cast<int>(
                        ((mask >> position) & 1U) != 0U
                    );
                }
                if ((minus_in_cut & 1) == 0) {
                    continue;
                }
                const std::int64_t rank =
                    invariant(subset(labels, mask, true), level)
                    * invariant(
                        subset(labels, mask, false),
                        level
                    );
                const unsigned int switched_mask =
                    selected_mask ^ mask;
                const std::int64_t switched =
                    invariant(
                        subset(labels, switched_mask, true),
                        level
                    )
                    * invariant(
                        subset(labels, switched_mask, false),
                        level
                    );
                if (rank == 0) {
                    continue;
                }
                const bool is_bad =
                    rank == 1 && switched == 0;
                bad += static_cast<std::int64_t>(is_bad);
                std::cout
                    << " D=" << mask
                    << " intersection="
                    << popcount(mask & selected_mask)
                    << " rank=" << rank
                    << " switched=" << switched
                    << " bad=" << (is_bad ? 1 : 0)
                    << '\n';
            }
            std::cout
                << "bad=" << bad
                << " margin=" << (local - 1 - bad)
                << '\n';
            return EXIT_SUCCESS;
        }
        if (
            argc == 3
            && std::string(argv[1]) == "--plus-neighbor"
        ) {
            scan_plus_neighbor_rank_two(parse_positive(argv[2]));
            return EXIT_SUCCESS;
        }
        if (
            argc == 3
            && std::string(argv[1]) == "--no-two-profile"
        ) {
            scan_no_two_selected_profile(parse_positive(argv[2]));
            return EXIT_SUCCESS;
        }
        if (argc == 10 && std::string(argv[1]) == "--word") {
            const int level = parse_positive(argv[2]);
            Labels labels;
            for (int index = 0; index < 7; ++index) {
                labels.push_back(parse_positive(argv[index + 3]));
            }
            std::int64_t pair = 0;
            std::int64_t completed_pair = 0;
            for (std::size_t first = 0U;
                 first < labels.size(); ++first) {
                for (std::size_t second = first + 1U;
                     second < labels.size(); ++second) {
                    const bool same_sign =
                        (first < 6U) == (second < 6U);
                    if (!same_sign
                        || labels[first] != labels[second]) {
                        continue;
                    }
                    Labels rest;
                    for (std::size_t index = 0U;
                         index < labels.size(); ++index) {
                        if (index != first && index != second) {
                            rest.push_back(labels[index]);
                        }
                    }
                    const std::int64_t retained = local_channels(
                        Labels{rest[0], rest[1]},
                        Labels{rest[2], rest[3], rest[4]},
                        level
                    );
                    pair += retained;
                    int same_sign_copies = 0;
                    for (std::size_t index = 0U;
                         index < labels.size(); ++index) {
                        same_sign_copies += static_cast<int>(
                            ((index < 6U) == (first < 6U))
                            && labels[index] == 2
                        );
                    }
                    completed_pair += (
                        labels[first] == 2
                        && same_sign_copies >= 3
                    ) ? invariant(rest, level) : retained;
                }
            }
            std::int64_t positive = 0;
            std::int64_t negative = 0;
            std::int64_t maximum = 0;
            constexpr unsigned int selected_mask = 7U;
            for (unsigned int mask = 0U;
                 mask < (1U << 7U); ++mask) {
                if (popcount(mask) != 3) {
                    continue;
                }
                const std::int64_t term =
                    invariant(subset(labels, mask, true), level)
                    * invariant(subset(labels, mask, false), level);
                int minus_in_cut = 0;
                for (int position = 0; position < 6; ++position) {
                    minus_in_cut += static_cast<int>(
                        (mask >> position) & 1U
                    );
                }
                if ((minus_in_cut & 1) == 0) {
                    positive += term;
                } else {
                    negative += term;
                    maximum = std::max(maximum, term);
                }
            }
            const Labels selected =
                subset(labels, selected_mask, true);
            const Labels complement =
                subset(labels, selected_mask, false);
            const std::int64_t local =
                local_channels(selected, complement, level);
            std::cout
                << "word level=" << level
                << " selected_rank="
                << invariant(selected, level)
                    * invariant(complement, level)
                << " maximum=" << maximum
                << " L=" << local
                << " R=" << pair
                << " R₂=" << completed_pair
                << " U=" << positive
                << " T=" << negative
                << " margin=" << (
                    local + pair + positive - negative
                )
                << " margin₂=" << (
                    local + completed_pair + positive - negative
                )
                << '\n';
            return EXIT_SUCCESS;
        }
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_seven_shallow_d12 maximum_level "
                "| --plus-neighbor maximum_level "
                "| --no-two-profile maximum_level "
                "| --switch-word k minus_count mask "
                "l0 l1 l2 l3 l4 l5 l6 "
                "| --word k l0 l1 l2 l3 l4 l5 l6"
            );
        }
        const int maximum_level = parse_positive(argv[1]);
        std::array<
            std::array<
                std::array<
                    std::array<SmallRankSummary, 2>,
                    3
                >,
                2
            >,
            8
        > summaries{};
        SmallRankWitness lemma_b16;
        std::uint64_t lemma_b16_cases = 0;
        std::array<std::array<SmallRankSummary, 4>, 7>
            orbit6_neighbor{};
        for (int level = 4; level <= maximum_level; ++level) {
            for (int orbit_index = 1;
                 orbit_index < 8; ++orbit_index) {
                scan_small_rank_level(
                    orbit_index,
                    level,
                    summaries,
                    lemma_b16,
                    lemma_b16_cases,
                    orbit6_neighbor
                );
            }
            std::cout << "SU2_SEVEN_SHALLOW_D12 progress="
                      << (level - 3) << '/'
                      << (maximum_level - 3)
                      << " level=" << level << '\n' << std::flush;
        }

        for (int orbit_index = 1;
             orbit_index < 8; ++orbit_index) {
            for (int rank = 1; rank <= 2; ++rank) {
                for (int kind = 0; kind < 3; ++kind) {
                    for (int side = 0; side < 2; ++side) {
                        const SmallRankSummary& summary =
                            summaries[
                                static_cast<std::size_t>(orbit_index)
                            ][static_cast<std::size_t>(rank - 1)]
                             [static_cast<std::size_t>(kind)]
                             [static_cast<std::size_t>(side)];
                        if (summary.cases == 0U) {
                            continue;
                        }
                        std::cout
                            << "orbit=" << orbit_index
                            << " rank=" << rank
                            << " kind=" << kind
                            << " side="
                            << (side == 0 ? "selected" : "complement")
                            << " cases=" << summary.cases;
                        print_small_rank_witness(
                            "direct", summary.direct
                        );
                        print_small_rank_witness(
                            "ceiling", summary.ceiling
                        );
                        print_small_rank_witness(
                            "local", summary.local
                        );
                        std::cout << '\n';
                    }
                }
            }
        }
        std::cout << "lemma=5A7B16 cases=" << lemma_b16_cases;
        print_small_rank_witness("margin", lemma_b16);
        std::cout << '\n';
        if (lemma_b16.initialized && lemma_b16.value < 0) {
            throw std::runtime_error(
                "bounded counterexample to Lemma 5A7B16"
            );
        }
        std::cout
            << "rank_one_switch cases="
            << rank_one_switch_cases;
        print_small_rank_witness(
            "total_margin",
            rank_one_switch_margin
        );
        print_small_rank_witness(
            "negative_bad_max",
            rank_one_switch_bad_maximum
        );
        for (std::size_t intersection = 0U;
             intersection
                 < rank_one_switch_component_margin.size();
             ++intersection) {
            const std::string name =
                "intersection_"
                + std::to_string(intersection)
                + "_margin";
            print_small_rank_witness(
                name.c_str(),
                rank_one_switch_component_margin[intersection]
            );
            const std::string maximum_name =
                "intersection_"
                + std::to_string(intersection)
                + "_negative_bad_max";
            print_small_rank_witness(
                maximum_name.c_str(),
                rank_one_switch_component_bad_maximum[
                    intersection
                ]
            );
        }
        std::cout << '\n';
        if (
            rank_one_switch_margin.initialized
            && rank_one_switch_margin.value < 0
        ) {
            throw std::runtime_error(
                "bounded counterexample to rank-one "
                "bad-switch payment"
            );
        }
        for (std::size_t kind = 0U;
             kind < rank_one_switch_endpoint_margin.size();
             ++kind) {
            for (std::size_t side = 0U;
                 side
                     < rank_one_switch_endpoint_margin[kind].size();
                 ++side) {
                std::cout
                    << "rank_one_switch_endpoint kind=" << kind
                    << " side="
                    << (side == 0U ? "selected" : "complement");
                print_small_rank_witness(
                    "margin",
                    rank_one_switch_endpoint_margin[kind][side]
                );
                std::cout << '\n';
            }
        }
        for (std::size_t total = 1U;
             total < orbit6_neighbor.size(); ++total) {
            for (std::size_t selected = 0U;
                 selected < orbit6_neighbor[total].size();
                 ++selected) {
                const SmallRankSummary& summary =
                    orbit6_neighbor[total][selected];
                if (summary.cases == 0U) {
                    continue;
                }
                std::cout
                    << "orbit6_neighbor total=" << total
                    << " selected=" << selected
                    << " cases=" << summary.cases;
                print_small_rank_witness(
                    "local_pair_margin", summary.direct
                );
                print_small_rank_witness(
                    "local_with_margin",
                    neighbor_local_with_payment[total][selected]
                );
                print_small_rank_witness(
                    "pair_without_margin",
                    neighbor_pair_without_payment[total][selected]
                );
                print_small_rank_witness(
                    "local_without_margin",
                    neighbor_local_without_payment[total][selected]
                );
                print_small_rank_witness(
                    "pair_with_margin",
                    neighbor_pair_with_payment[total][selected]
                );
                print_small_rank_witness(
                    "local_total_margin",
                    neighbor_local_total_payment[total][selected]
                );
                std::cout << '\n';
                const bool lemma_b17 =
                    total == 2U && selected == 2U;
                const bool lemma_b18 =
                    total == 3U && selected == 2U;
                const bool lemma_b19 =
                    total == 1U && selected == 1U;
                const bool lemma_b20 =
                    total == 2U && selected == 1U;
                if (
                    (lemma_b17 || lemma_b20)
                    && neighbor_local_total_payment[total][selected]
                        .value < 0
                ) {
                    throw std::runtime_error(
                        "bounded counterexample to a direct "
                        "rank-two neighbor lemma"
                    );
                }
                if (
                    (lemma_b18 || lemma_b19)
                    && summary.direct.value < 0
                ) {
                    throw std::runtime_error(
                        "bounded counterexample to an equality-paid "
                        "rank-two neighbor lemma"
                    );
                }
            }
        }
        for (std::size_t index = 0U;
             index < double_neighbor_margin.size(); ++index) {
            if (!double_neighbor_margin[index].initialized) {
                continue;
            }
            std::cout << "double_neighbor_d2 N=" << (index + 1U);
            print_small_rank_witness(
                "local_margin", double_neighbor_margin[index]
            );
            print_small_rank_witness(
                "W_payment", double_neighbor_w_payment[index]
            );
            print_small_rank_witness(
                "tail_payment",
                double_neighbor_tail_payment[index]
            );
            std::cout << '\n';
        }
        for (std::size_t sign = 0U;
             sign < neighbor_zero_selected_by_sign.size(); ++sign) {
            for (std::size_t category = 0U;
                 category
                     < neighbor_zero_selected_by_sign[sign].size();
                 ++category) {
                const SmallRankSummary& summary =
                    neighbor_zero_selected_by_sign[sign][category];
                if (summary.cases == 0U) {
                    continue;
                }
                std::cout
                    << "neighbor_zero_selected sign="
                    << (sign == 0U ? "minus" : "plus")
                    << " category=" << category
                    << " cases=" << summary.cases;
                print_small_rank_witness(
                    "direct_margin", summary.direct
                );
                print_small_rank_witness(
                    "bound38_margin", summary.ceiling
                );
                print_small_rank_witness(
                    "raw40_deficit", summary.local
                );
                print_small_rank_witness(
                    "supply",
                    neighbor_zero_selected_supply[sign][category]
                );
                print_small_rank_witness(
                    "local",
                    neighbor_zero_selected_local[sign][category]
                );
                print_small_rank_witness(
                    "pair",
                    neighbor_zero_selected_pair[sign][category]
                );
                std::cout << '\n';
            }
        }
        std::cout << "SU2_SEVEN_SHALLOW_D12 PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SEVEN_SHALLOW_D12 FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
