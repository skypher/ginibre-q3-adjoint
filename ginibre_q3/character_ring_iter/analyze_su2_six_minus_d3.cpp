#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#define main analyze_su2_four_minus_shallow_main
#include "analyze_su2_four_minus_shallow.cpp"
#undef main

namespace {

struct D3Witness {
    bool initialized = false;
    std::int64_t value = 0;
    int level = 0;
    Labels minus;
    int plus = 0;
    Labels selected;
    std::int64_t local = 0;
    std::int64_t pair = 0;
    std::int64_t positive = 0;
    std::int64_t negative = 0;
};

struct D3Stratum {
    std::uint64_t selected_cuts = 0;
    D3Witness exact_margin;
    D3Witness ceiling_margin;
    D3Witness local_ceiling;
    D3Witness anchor_pair_ceiling;
    D3Witness no_pair_local_ceiling;
    D3Witness sole_pair_ceiling;
    D3Witness equality_pair_ceiling;
    D3Witness pair_ceiling;
    D3Witness positive_ceiling;
    D3Witness local_margin;
    D3Witness demand_defect;
};

constexpr int stratum_count = 6;

struct EqualityPattern {
    std::uint64_t selected_cuts = 0;
    D3Witness pair_ceiling;
};

int selected_stratum(const Labels& selected) {
    if (selected[0] != selected[1]
        && selected[1] != selected[2]) {
        return 0;
    }
    if (selected[0] == selected[2]) {
        return 1;
    }
    int singleton = selected[0];
    if (selected[0] == selected[1]) {
        singleton = selected[2];
    } else if (selected[1] == selected[2]) {
        singleton = selected[0];
    }
    if (singleton == 2) {
        return 2;
    }
    if (singleton == 4) {
        return 3;
    }
    if (singleton == 6) {
        return 4;
    }
    return 5;
}

void consider_d3(
    D3Witness& witness,
    std::int64_t value,
    int level,
    const Labels& minus,
    int plus,
    const Labels& selected,
    std::int64_t local,
    std::int64_t pair,
    std::int64_t positive,
    std::int64_t negative
) {
    if (!witness.initialized || value < witness.value) {
        witness = D3Witness{
            true,
            value,
            level,
            minus,
            plus,
            selected,
            local,
            pair,
            positive,
            negative
        };
    }
}

void print_d3_witness(
    const char* name,
    const D3Witness& witness
) {
    std::cout << ' ' << name << '=';
    if (!witness.initialized) {
        std::cout << "none";
        return;
    }
    std::cout << witness.value << "@k" << witness.level << ":M";
    print_labels(witness.minus);
    std::cout << ":P[" << witness.plus << "]:A";
    print_labels(witness.selected);
    std::cout << ":(L,Peq,U,T)=("
              << witness.local << ','
              << witness.pair << ','
              << witness.positive << ','
              << witness.negative << ')';
}

std::int64_t low_local(
    const Labels& selected,
    const Labels& complement,
    int level
) {
    std::int64_t result = 0;
    for (int output = 0; output <= 6; output += 2) {
        result += multiplicity(selected, output, level)
            * multiplicity(complement, output, level);
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_six_minus_d3 maximum_level"
            );
        }
        const int maximum_level = parse_positive(argv[1]);
        if (maximum_level < 6) {
            throw std::runtime_error(
                "maximum level must be at least six"
            );
        }

        std::array<D3Stratum, stratum_count> strata{};
        std::array<
            std::array<std::array<EqualityPattern, 4>, 10>,
            4
        > equality_patterns{};
        std::uint64_t words = 0;
        std::uint64_t distinct_short_checks = 0;
        std::uint64_t distinct_cross_checks = 0;
        for (int level = 6; level <= maximum_level; ++level) {
            Labels alphabet;
            for (int label = 2; label <= level - 2; ++label) {
                alphabet.push_back(label);
            }
            const std::vector<Labels> minus_lists =
                sorted_lists(alphabet, 6);
            for (const Labels& minus : minus_lists) {
                for (const int plus : alphabet) {
                    if (std::binary_search(
                            minus.begin(), minus.end(), plus
                        )) {
                        continue;
                    }
                    int total = plus;
                    for (const int label : minus) {
                        total += label;
                    }
                    if ((total & 1) != 0) {
                        continue;
                    }
                    ++words;
                    Labels labels = minus;
                    labels.push_back(plus);

                    std::int64_t pair = 0;
                    int equality_pairs = 0;
                    for (std::size_t first = 0U;
                         first < minus.size(); ++first) {
                        for (std::size_t second = first + 1U;
                             second < minus.size(); ++second) {
                            if (minus[first] != minus[second]) {
                                continue;
                            }
                            ++equality_pairs;
                            Labels rest;
                            for (std::size_t index = 0U;
                                 index < labels.size(); ++index) {
                                if (index != first
                                    && index != second) {
                                    rest.push_back(labels[index]);
                                }
                            }
                            pair += invariant(rest, level);
                        }
                    }

                    std::int64_t positive = 0;
                    for (unsigned int first = 0U; first < 6U; ++first) {
                        for (unsigned int second = first + 1U;
                             second < 6U; ++second) {
                            const unsigned int mask =
                                (1U << 6U)
                                | (1U << first)
                                | (1U << second);
                            positive +=
                                invariant(subset(labels, mask, true), level)
                                * invariant(
                                    subset(labels, mask, false), level
                                );
                        }
                    }

                    std::int64_t negative = 0;
                    std::int64_t maximum_cut = 0;
                    std::vector<unsigned int> maximal_masks;
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
                                    invariant(
                                        subset(labels, mask, true), level
                                    )
                                    * invariant(
                                        subset(labels, mask, false), level
                                    );
                                negative += term;
                                if (term > maximum_cut) {
                                    maximum_cut = term;
                                    maximal_masks.clear();
                                    maximal_masks.push_back(mask);
                                } else if (term == maximum_cut) {
                                    maximal_masks.push_back(mask);
                                }
                            }
                        }
                    }
                    if (maximum_cut != 3) {
                        continue;
                    }

                    for (const unsigned int mask : maximal_masks) {
                        const Labels selected =
                            subset(labels, mask, true);
                        const Labels complement =
                            subset(labels, mask, false);
                        const std::int64_t local =
                            low_local(selected, complement, level);
                        int selected_pairs = 0;
                        int complement_pairs = 0;
                        int cross_pairs = 0;
                        for (std::size_t first = 0U;
                             first < selected.size(); ++first) {
                            for (std::size_t second = first + 1U;
                                 second < selected.size(); ++second) {
                                selected_pairs += static_cast<int>(
                                    selected[first] == selected[second]
                                );
                            }
                            for (std::size_t second = 0U;
                                 second + 1U < complement.size(); ++second) {
                                cross_pairs += static_cast<int>(
                                    selected[first] == complement[second]
                                );
                            }
                        }
                        for (std::size_t first = 0U;
                             first + 1U < complement.size(); ++first) {
                            for (std::size_t second = first + 1U;
                                 second + 1U < complement.size(); ++second) {
                                complement_pairs += static_cast<int>(
                                    complement[first] == complement[second]
                                );
                            }
                        }
                        EqualityPattern& equality_pattern =
                            equality_patterns[
                                static_cast<std::size_t>(selected_pairs)
                            ][static_cast<std::size_t>(cross_pairs)]
                            [static_cast<std::size_t>(complement_pairs)];
                        ++equality_pattern.selected_cuts;
                        consider_d3(
                            equality_pattern.pair_ceiling,
                            local + pair - 60,
                            level, minus, plus, selected,
                            local, pair, positive, negative
                        );
                        if (equality_pairs == 0
                            && selected_stratum(selected) == 0) {
                            const std::int64_t b2 =
                                multiplicity(complement, 2, level);
                            const std::int64_t b4 =
                                multiplicity(complement, 4, level);
                            const std::int64_t b6 =
                                multiplicity(complement, 6, level);
                            const std::int64_t base_payment =
                                3 + 2 * b2 + 3 * b4 + 2 * b6;
                            if (base_payment < 60) {
                                ++distinct_short_checks;
                                if (multiplicity(selected, 6, level) < 4) {
                                    throw std::runtime_error(
                                        "distinct short band lacks "
                                        "four output-six paths"
                                    );
                                }
                            }
                        }
                        if (selected_stratum(selected) == 0
                            && local < 60) {
                            ++distinct_cross_checks;
                            if (local + pair < 60) {
                                throw std::runtime_error(
                                    "distinct equal-pair payment failed"
                                );
                            }
                        }
                        if (local + pair < 60) {
                            throw std::runtime_error(
                                "rank-three equal-pair payment failed"
                            );
                        }
                        std::int64_t anchor_pair = 0;
                        if (selected[0] == selected[1]) {
                            anchor_pair = multiplicity(
                                complement, selected[2], level
                            );
                        } else if (selected[1] == selected[2]) {
                            anchor_pair = multiplicity(
                                complement, selected[0], level
                            );
                        }
                        D3Stratum& stratum = strata[
                            static_cast<std::size_t>(
                                selected_stratum(selected)
                            )
                        ];
                        ++stratum.selected_cuts;
                        consider_d3(
                            stratum.exact_margin,
                            local + pair + positive - negative,
                            level, minus, plus, selected,
                            local, pair, positive, negative
                        );
                        consider_d3(
                            stratum.ceiling_margin,
                            local + pair + positive - 60,
                            level, minus, plus, selected,
                            local, pair, positive, negative
                        );
                        consider_d3(
                            stratum.local_ceiling,
                            local - 60,
                            level, minus, plus, selected,
                            local, pair, positive, negative
                        );
                        if (anchor_pair > 0) {
                            consider_d3(
                                stratum.anchor_pair_ceiling,
                                local + anchor_pair - 60,
                                level, minus, plus, selected,
                                local, anchor_pair, positive, negative
                            );
                        }
                        if (equality_pairs == 0) {
                            consider_d3(
                                stratum.no_pair_local_ceiling,
                                local - 60,
                                level, minus, plus, selected,
                                local, pair, positive, negative
                            );
                        }
                        if (equality_pairs == 1 && anchor_pair > 0) {
                            consider_d3(
                                stratum.sole_pair_ceiling,
                                local + anchor_pair - 60,
                                level, minus, plus, selected,
                                local, anchor_pair, positive, negative
                            );
                        }
                        if (equality_pairs > 0) {
                            consider_d3(
                                stratum.equality_pair_ceiling,
                                local + pair - 60,
                                level, minus, plus, selected,
                                local, pair, positive, negative
                            );
                        }
                        consider_d3(
                            stratum.pair_ceiling,
                            local + pair - 60,
                            level, minus, plus, selected,
                            local, pair, positive, negative
                        );
                        consider_d3(
                            stratum.positive_ceiling,
                            local + positive - 60,
                            level, minus, plus, selected,
                            local, pair, positive, negative
                        );
                        consider_d3(
                            stratum.local_margin,
                            local + pair - negative,
                            level, minus, plus, selected,
                            local, pair, positive, negative
                        );
                        consider_d3(
                            stratum.demand_defect,
                            60 - negative,
                            level, minus, plus, selected,
                            local, pair, positive, negative
                        );
                    }
                }
            }
        }

        constexpr std::array<const char*, stratum_count> names{{
            "distinct", "all_equal", "c2", "c4", "c6", "c8plus"
        }};
        std::cout << "SU2_SIX_MINUS_D3 maximum_level="
                  << maximum_level << " words=" << words
                  << " distinct_short_checks="
                  << distinct_short_checks
                  << " distinct_cross_checks="
                  << distinct_cross_checks << '\n';
        for (std::size_t index = 0U;
             index < strata.size(); ++index) {
            const D3Stratum& stratum = strata[index];
            std::cout << "type=" << names[index]
                      << " selected_cuts=" << stratum.selected_cuts;
            print_d3_witness("exact_margin", stratum.exact_margin);
            print_d3_witness("ceiling_margin", stratum.ceiling_margin);
            print_d3_witness("local_ceiling", stratum.local_ceiling);
            print_d3_witness(
                "anchor_pair_ceiling",
                stratum.anchor_pair_ceiling
            );
            print_d3_witness(
                "no_pair_local_ceiling",
                stratum.no_pair_local_ceiling
            );
            print_d3_witness(
                "sole_pair_ceiling",
                stratum.sole_pair_ceiling
            );
            print_d3_witness(
                "equality_pair_ceiling",
                stratum.equality_pair_ceiling
            );
            print_d3_witness("pair_ceiling", stratum.pair_ceiling);
            print_d3_witness(
                "positive_ceiling", stratum.positive_ceiling
            );
            print_d3_witness("local_margin", stratum.local_margin);
            print_d3_witness("demand_defect", stratum.demand_defect);
            std::cout << '\n';
        }
        for (std::size_t selected_pairs = 0U;
             selected_pairs < equality_patterns.size();
             ++selected_pairs) {
            for (std::size_t cross_pairs = 0U;
                 cross_pairs
                     < equality_patterns[selected_pairs].size();
                 ++cross_pairs) {
                for (std::size_t complement_pairs = 0U;
                     complement_pairs
                         < equality_patterns[selected_pairs][cross_pairs]
                               .size();
                     ++complement_pairs) {
                    const EqualityPattern& pattern =
                        equality_patterns[selected_pairs][cross_pairs]
                            [complement_pairs];
                    if (pattern.selected_cuts == 0U) {
                        continue;
                    }
                    std::cout
                        << "equality_pattern selected_pairs="
                        << selected_pairs
                        << " cross_pairs=" << cross_pairs
                        << " complement_pairs=" << complement_pairs
                        << " selected_cuts=" << pattern.selected_cuts;
                    print_d3_witness(
                        "pair_ceiling", pattern.pair_ceiling
                    );
                    std::cout << '\n';
                }
            }
        }
        std::cout << "SU2_SIX_MINUS_D3 PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SIX_MINUS_D3 FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
