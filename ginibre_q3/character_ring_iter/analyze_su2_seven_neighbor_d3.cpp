#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#define main analyze_su2_four_minus_shallow_main
#include "analyze_su2_four_minus_shallow.cpp"
#undef main

namespace {

struct NeighborWitness {
    bool initialized = false;
    std::int64_t value = 0;
    int level = 0;
    Labels labels;
    std::int64_t local = 0;
    std::int64_t pair = 0;
    std::int64_t positive = 0;
    std::int64_t negative = 0;
    std::array<std::int64_t, 20> negative_terms{};
    std::array<std::int64_t, 15> positive_terms{};
    std::array<std::int64_t, 29> auxiliary{};
};

struct NeighborSummary {
    std::uint64_t cases = 0;
    NeighborWitness exact_margin;
    NeighborWitness local_pair_margin;
    NeighborWitness local_pair_ceiling;
    NeighborWitness local_margin;
    NeighborWitness pair_margin;
    NeighborWitness positive_margin;
    NeighborWitness negative;
    NeighborWitness edge_to_output_two;
    NeighborWitness output_four_residual;
    std::array<NeighborWitness, 4> local_demand_by_n;
    std::array<NeighborWitness, 4> maximum_edges_by_n;
    std::array<std::array<NeighborWitness, 7>, 2>
        sole_ceiling_by_selected_pair_and_cross;
    std::array<std::array<NeighborWitness, 7>, 2>
        sole_local_by_selected_pair_and_cross;
};

void consider_neighbor(
    NeighborWitness& witness,
    std::int64_t value,
    int level,
    const Labels& labels,
    std::int64_t local,
    std::int64_t pair,
    std::int64_t positive,
    std::int64_t negative,
    const std::array<std::int64_t, 20>& negative_terms,
    const std::array<std::int64_t, 15>& positive_terms,
    const std::array<std::int64_t, 29>& auxiliary
) {
    if (!witness.initialized || value < witness.value) {
        witness = NeighborWitness{
            true,
            value,
            level,
            labels,
            local,
            pair,
            positive,
            negative,
            negative_terms,
            positive_terms,
            auxiliary
        };
    }
}

std::array<std::int64_t, 2> pair_reservoir_exact(
    const Labels& labels,
    int level
) {
    std::array<std::int64_t, 2> result{};
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
            result[first == 0U && second == 1U ? 0U : 1U]
                += invariant(rest, level);
        }
    }
    return result;
}

std::int64_t selected_local_exact(
    const Labels& labels,
    int level
) {
    const Labels selected{labels[0], labels[1], labels[2]};
    const Labels complement{
        labels[3], labels[4], labels[5], labels[6]
    };
    std::int64_t result = 0;
    for (int low = 0; low <= 8; low += 2) {
        result += multiplicity(selected, low, level)
            * multiplicity(complement, low, level);
    }
    for (int reflected_low = 1;
         reflected_low <= 9; reflected_low += 2) {
        const int high = level - reflected_low;
        if (high <= 8) {
            continue;
        }
        result += multiplicity(selected, high, level)
            * multiplicity(complement, high, level);
    }
    return result;
}

void cut_terms(
    const Labels& labels,
    int level,
    std::array<std::int64_t, 20>& negative_terms,
    std::array<std::int64_t, 15>& positive_terms,
    std::int64_t& negative,
    std::int64_t& positive
) {
    negative = 0;
    positive = 0;
    std::size_t negative_index = 0U;
    std::size_t positive_index = 0U;
    for (unsigned int first = 0U; first < 6U; ++first) {
        for (unsigned int second = first + 1U;
             second < 6U; ++second) {
            for (unsigned int third = second + 1U;
                 third < 6U; ++third) {
                const unsigned int mask =
                    (1U << first)
                    | (1U << second)
                    | (1U << third);
                const std::int64_t term =
                    invariant(subset(labels, mask, true), level)
                    * invariant(subset(labels, mask, false), level);
                negative_terms[negative_index] = term;
                ++negative_index;
                negative += term;
            }
        }
    }
    for (unsigned int first = 0U; first < 6U; ++first) {
        for (unsigned int second = first + 1U;
             second < 6U; ++second) {
            const unsigned int mask =
                (1U << 6U)
                | (1U << first)
                | (1U << second);
            const std::int64_t term =
                invariant(subset(labels, mask, true), level)
                * invariant(subset(labels, mask, false), level);
            positive_terms[positive_index] = term;
            ++positive_index;
            positive += term;
        }
    }
}

void print_neighbor_witness(
    const char* name,
    const NeighborWitness& witness
) {
    std::cout << ' ' << name << '=';
    if (!witness.initialized) {
        std::cout << "none";
        return;
    }
    std::cout << witness.value << "@k" << witness.level << ":L";
    print_labels(witness.labels);
    std::cout << ":(L,Peq,U,T)=("
              << witness.local << ','
              << witness.pair << ','
              << witness.positive << ','
              << witness.negative << "):neg[";
    for (std::size_t index = 0U;
         index < witness.negative_terms.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << witness.negative_terms[index];
    }
    std::cout << "]:pos[";
    for (std::size_t index = 0U;
         index < witness.positive_terms.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << witness.positive_terms[index];
    }
    std::cout << "]:aux[";
    for (std::size_t index = 0U;
         index < witness.auxiliary.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << witness.auxiliary[index];
    }
    std::cout << ']';
}

void record_case(
    NeighborSummary& summary,
    int level,
    const Labels& labels,
    bool doubled
) {
    std::array<std::int64_t, 20> negative_terms{};
    std::array<std::int64_t, 15> positive_terms{};
    std::int64_t negative = 0;
    std::int64_t positive = 0;
    cut_terms(
        labels,
        level,
        negative_terms,
        positive_terms,
        negative,
        positive
    );
    if (negative_terms[0] != 3) {
        return;
    }
    for (const std::int64_t term : negative_terms) {
        if (term > 3) {
            return;
        }
    }
    const std::int64_t local =
        selected_local_exact(labels, level);
    const std::array<std::int64_t, 2> pair_parts =
        pair_reservoir_exact(labels, level);
    const std::int64_t pair = pair_parts[0] + pair_parts[1];
    const Labels complement{
        labels[3], labels[4], labels[5], labels[6]
    };
    std::array<std::int64_t, 29> auxiliary{};
    auxiliary[0] = multiplicity(complement, 2, level);
    auxiliary[1] = multiplicity(complement, 4, level);
    std::size_t negative_index = 0U;
    for (unsigned int first = 0U; first < 6U; ++first) {
        for (unsigned int second = first + 1U;
             second < 6U; ++second) {
            for (unsigned int third = second + 1U;
                 third < 6U; ++third) {
                int endpoint_count = 0;
                if (doubled) {
                    endpoint_count =
                        static_cast<int>(first == 1U || first == 2U)
                        + static_cast<int>(
                            second == 1U || second == 2U
                        )
                        + static_cast<int>(
                            third == 1U || third == 2U
                        );
                    auxiliary[
                        endpoint_count == 2
                            ? 2U
                            : (endpoint_count == 1 ? 3U : 4U)
                    ] += negative_terms[negative_index];
                } else {
                    const bool contains =
                        first == 2U || second == 2U || third == 2U;
                    auxiliary[contains ? 2U : 3U]
                        += negative_terms[negative_index];
                }
                ++negative_index;
            }
        }
    }
    auxiliary[5] = pair_parts[0];
    auxiliary[6] = pair_parts[1];
    int copies_two = 0;
    for (std::size_t index : {0U, 3U, 4U, 5U}) {
        copies_two += static_cast<int>(labels[index] == 2);
    }
    auxiliary[7] = copies_two;
    int edge_count = 0;
    constexpr std::array<std::size_t, 4> q_indices{{
        0U, 3U, 4U, 5U
    }};
    for (std::size_t first = 0U;
         first < q_indices.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < q_indices.size(); ++second) {
            const int sum = labels[q_indices[first]]
                + labels[q_indices[second]];
            edge_count += static_cast<int>(
                sum == level - 1 || sum == level + 1
            );
        }
    }
    auxiliary[8] = edge_count;
    int cross_equalities = 0;
    for (std::size_t selected_index = 0U;
         selected_index < 2U; ++selected_index) {
        for (std::size_t complement_index = 3U;
             complement_index < 6U; ++complement_index) {
            cross_equalities += static_cast<int>(
                labels[selected_index] == labels[complement_index]
            );
        }
    }
    const Labels selected{labels[0], labels[1], labels[2]};
    for (int index = 0; index < 5; ++index) {
        const int low = 2 * index;
        const int high = level - (2 * index + 1);
        auxiliary[static_cast<std::size_t>(9 + index)] =
            multiplicity(selected, low, level);
        auxiliary[static_cast<std::size_t>(14 + index)] =
            high > 8 ? multiplicity(selected, high, level) : 0;
        auxiliary[static_cast<std::size_t>(19 + index)] =
            multiplicity(complement, low, level);
        auxiliary[static_cast<std::size_t>(24 + index)] =
            high > 8 ? multiplicity(complement, high, level) : 0;
    }
    ++summary.cases;
    consider_neighbor(
        summary.exact_margin,
        local + pair + positive - negative,
        level,
        labels,
        local,
        pair,
        positive,
        negative,
        negative_terms,
        positive_terms,
        auxiliary
    );
    consider_neighbor(
        summary.local_pair_margin,
        local + pair - negative,
        level,
        labels,
        local,
        pair,
        positive,
        negative,
        negative_terms,
        positive_terms,
        auxiliary
    );
    consider_neighbor(
        summary.local_pair_ceiling,
        local + pair - 60,
        level,
        labels,
        local,
        pair,
        positive,
        negative,
        negative_terms,
        positive_terms,
        auxiliary
    );
    consider_neighbor(
        summary.local_margin,
        local - negative,
        level,
        labels,
        local,
        pair,
        positive,
        negative,
        negative_terms,
        positive_terms,
        auxiliary
    );
    consider_neighbor(
        summary.pair_margin,
        pair - negative,
        level,
        labels,
        local,
        pair,
        positive,
        negative,
        negative_terms,
        positive_terms,
        auxiliary
    );
    consider_neighbor(
        summary.positive_margin,
        positive - negative,
        level,
        labels,
        local,
        pair,
        positive,
        negative,
        negative_terms,
        positive_terms,
        auxiliary
    );
    consider_neighbor(
        summary.negative,
        -negative,
        level,
        labels,
        local,
        pair,
        positive,
        negative,
        negative_terms,
        positive_terms,
        auxiliary
    );
    if (doubled && copies_two >= 1 && copies_two <= 4) {
        const std::int64_t one_endpoint_demand = auxiliary[3];
        const std::int64_t no_endpoint_demand = auxiliary[4];
        consider_neighbor(
            summary.edge_to_output_two,
            auxiliary[0] - one_endpoint_demand / 2,
            level,
            labels,
            local,
            pair,
            positive,
            negative,
            negative_terms,
            positive_terms,
            auxiliary
        );
        consider_neighbor(
            summary.output_four_residual,
            auxiliary[1]
                - 3 * (copies_two - 1)
                - no_endpoint_demand,
            level,
            labels,
            local,
            pair,
            positive,
            negative,
            negative_terms,
            positive_terms,
            auxiliary
        );
        consider_neighbor(
            summary.local_demand_by_n[
                static_cast<std::size_t>(copies_two - 1)
            ],
            local - negative,
            level,
            labels,
            local,
            pair,
            positive,
            negative,
            negative_terms,
            positive_terms,
            auxiliary
        );
        consider_neighbor(
            summary.maximum_edges_by_n[
                static_cast<std::size_t>(copies_two - 1)
            ],
            -edge_count,
            level,
            labels,
            local,
            pair,
            positive,
            negative,
            negative_terms,
            positive_terms,
            auxiliary
        );
    } else if (!doubled) {
        const std::size_t selected_pair =
            labels[0] == labels[1] ? 1U : 0U;
        const std::size_t cross =
            static_cast<std::size_t>(cross_equalities);
        consider_neighbor(
            summary.sole_ceiling_by_selected_pair_and_cross[
                selected_pair
            ][cross],
            local + pair - 60,
            level,
            labels,
            local,
            pair,
            positive,
            negative,
            negative_terms,
            positive_terms,
            auxiliary
        );
        consider_neighbor(
            summary.sole_local_by_selected_pair_and_cross[
                selected_pair
            ][cross],
            local,
            level,
            labels,
            local,
            pair,
            positive,
            negative,
            negative_terms,
            positive_terms,
            auxiliary
        );
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_seven_neighbor_d3 maximum_level"
            );
        }
        const int maximum_level = parse_positive(argv[1]);
        NeighborSummary sole;
        NeighborSummary doubled;
        for (int level = 5; level <= maximum_level; level += 2) {
            Labels deep;
            for (int label = 2; label <= level - 2; label += 2) {
                deep.push_back(label);
            }
            const std::vector<Labels> selected_pairs =
                sorted_lists(deep, 2);
            const std::vector<Labels> complements =
                sorted_lists(deep, 3);
            for (const Labels& complement_minus : complements) {
                for (const int plus : deep) {
                    if (std::binary_search(
                            complement_minus.begin(),
                            complement_minus.end(),
                            plus
                        )) {
                        continue;
                    }
                    Labels complement = complement_minus;
                    complement.push_back(plus);
                    if (invariant(complement, level) != 3) {
                        continue;
                    }
                    for (const Labels& selected_pair : selected_pairs) {
                        if (plus == selected_pair[0]
                            || plus == selected_pair[1]) {
                            continue;
                        }
                        Labels labels{
                            selected_pair[0],
                            selected_pair[1],
                            level - 1,
                            complement_minus[0],
                            complement_minus[1],
                            complement_minus[2],
                            plus
                        };
                        if (invariant(
                                Labels{
                                    labels[0], labels[1], labels[2]
                                },
                                level
                            ) == 1) {
                            record_case(sole, level, labels, false);
                        }
                    }
                    if (plus == 2) {
                        continue;
                    }
                    Labels labels{
                        2,
                        level - 1,
                        level - 1,
                        complement_minus[0],
                        complement_minus[1],
                        complement_minus[2],
                        plus
                    };
                    record_case(doubled, level, labels, true);
                }
            }
        }
        std::cout
            << "SU2_SEVEN_NEIGHBOR_D3 maximum_level="
            << maximum_level << '\n';
        for (const auto& [name, summary] :
             std::array<std::pair<const char*, NeighborSummary*>, 2>{{
                 {"sole", &sole},
                 {"double", &doubled}
             }}) {
            std::cout << "type=" << name
                      << " cases=" << summary->cases;
            print_neighbor_witness(
                "exact_margin", summary->exact_margin
            );
            print_neighbor_witness(
                "local_pair_margin",
                summary->local_pair_margin
            );
            print_neighbor_witness(
                "local_pair_ceiling",
                summary->local_pair_ceiling
            );
            print_neighbor_witness(
                "local_margin", summary->local_margin
            );
            print_neighbor_witness(
                "pair_margin", summary->pair_margin
            );
            print_neighbor_witness(
                "positive_margin", summary->positive_margin
            );
            print_neighbor_witness(
                "max_negative", summary->negative
            );
            if (name == std::string("double")) {
                print_neighbor_witness(
                    "edge_to_output_two",
                    summary->edge_to_output_two
                );
                print_neighbor_witness(
                    "output_four_residual",
                    summary->output_four_residual
                );
                for (std::size_t index = 0U;
                     index < summary->local_demand_by_n.size();
                     ++index) {
                    const std::string row =
                        "local_demand_N"
                        + std::to_string(index + 1U);
                    print_neighbor_witness(
                        row.c_str(),
                        summary->local_demand_by_n[index]
                    );
                    const std::string edge_row =
                        "maximum_edges_N"
                        + std::to_string(index + 1U);
                    print_neighbor_witness(
                        edge_row.c_str(),
                        summary->maximum_edges_by_n[index]
                    );
                }
            } else {
                for (std::size_t selected_pair = 0U;
                     selected_pair
                         < summary
                             ->sole_ceiling_by_selected_pair_and_cross
                             .size();
                     ++selected_pair) {
                    for (std::size_t cross = 0U;
                         cross
                             < summary
                                 ->sole_ceiling_by_selected_pair_and_cross[
                                     selected_pair
                                 ].size();
                         ++cross) {
                        const NeighborWitness& witness =
                            summary
                                ->sole_ceiling_by_selected_pair_and_cross[
                                    selected_pair
                                ][cross];
                        if (!witness.initialized) {
                            continue;
                        }
                        const std::string row =
                            "ceiling_selected_pair"
                            + std::to_string(selected_pair)
                            + "_cross"
                            + std::to_string(cross);
                        print_neighbor_witness(
                            row.c_str(), witness
                        );
                        const std::string local_row =
                            "local_selected_pair"
                            + std::to_string(selected_pair)
                            + "_cross"
                            + std::to_string(cross);
                        print_neighbor_witness(
                            local_row.c_str(),
                            summary
                                ->sole_local_by_selected_pair_and_cross[
                                    selected_pair
                                ][cross]
                        );
                    }
                }
            }
            std::cout << '\n';
        }
        std::cout << "SU2_SEVEN_NEIGHBOR_D3 PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SEVEN_NEIGHBOR_D3 FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
