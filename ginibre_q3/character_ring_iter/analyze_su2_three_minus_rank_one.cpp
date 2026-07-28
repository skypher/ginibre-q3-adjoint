#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <set>
#include <vector>

namespace {

constexpr unsigned selected_mask = 0b0000111U;
constexpr unsigned minus_mask = 0b0001111U;

bool fusion(int k, int first, int second, int output) {
    return output >= 0
        && output <= k
        && ((first + second + output) & 1) == 0
        && output >= std::abs(first - second)
        && output <= first + second
        && first + second + output <= 2 * k;
}

int triple_output(
    int k,
    int first,
    int second,
    int third,
    int output
) {
    int result = 0;
    for (int intermediate = 0; intermediate <= k; ++intermediate) {
        if (fusion(k, first, second, intermediate)
            && fusion(k, intermediate, third, output)) {
            ++result;
        }
    }
    return result;
}

int fourfold_invariant(
    int k,
    int first,
    int second,
    int third,
    int fourth
) {
    int result = 0;
    for (int intermediate = 0; intermediate <= k; ++intermediate) {
        if (fusion(k, first, second, intermediate)
            && fusion(k, third, fourth, intermediate)) {
            ++result;
        }
    }
    return result;
}

std::array<int, 3> triple_positions(unsigned mask) {
    std::array<int, 3> result{};
    int next = 0;
    for (int position = 0; position < 7; ++position) {
        if (((mask >> position) & 1U) != 0U) {
            result[static_cast<std::size_t>(next)] = position;
            ++next;
        }
    }
    return result;
}

std::array<int, 4> complement_positions(unsigned mask) {
    std::array<int, 4> result{};
    int next = 0;
    for (int position = 0; position < 7; ++position) {
        if (((mask >> position) & 1U) == 0U) {
            result[static_cast<std::size_t>(next)] = position;
            ++next;
        }
    }
    return result;
}

int triple_cut_weight(
    int k,
    const std::array<int, 7>& labels,
    unsigned mask
) {
    const std::array<int, 3> cut = triple_positions(mask);
    const std::array<int, 4> rest = complement_positions(mask);
    return static_cast<int>(fusion(
        k,
        labels[static_cast<std::size_t>(cut[0])],
        labels[static_cast<std::size_t>(cut[1])],
        labels[static_cast<std::size_t>(cut[2])]
    )) * fourfold_invariant(
        k,
        labels[static_cast<std::size_t>(rest[0])],
        labels[static_cast<std::size_t>(rest[1])],
        labels[static_cast<std::size_t>(rest[2])],
        labels[static_cast<std::size_t>(rest[3])]
    );
}

bool switched_weight_is_zero(
    int k,
    const std::array<int, 7>& labels,
    unsigned candidate
) {
    const unsigned difference = selected_mask ^ candidate;
    const int size = std::popcount(difference);
    if (size == 6) {
        return true;
    }
    if (size == 4) {
        return triple_cut_weight(
            k,
            labels,
            (~difference) & 0b1111111U
        ) == 0;
    }
    if (size == 2) {
        int first = -1;
        int second = -1;
        for (int position = 0; position < 7; ++position) {
            if (((difference >> position) & 1U) == 0U) {
                continue;
            }
            if (first < 0) {
                first = position;
            } else {
                second = position;
            }
        }
        return labels[static_cast<std::size_t>(first)]
            != labels[static_cast<std::size_t>(second)];
    }
    return false;
}

int bad_switch_count(
    int k,
    const std::array<int, 7>& labels
) {
    int result = 0;
    for (unsigned candidate = 0U;
         candidate < (1U << 7U);
         ++candidate) {
        if (std::popcount(candidate) != 3
            || candidate == selected_mask
            || (std::popcount(candidate & minus_mask) & 1) == 0) {
            continue;
        }
        if (triple_cut_weight(k, labels, candidate) == 1
            && switched_weight_is_zero(k, labels, candidate)) {
            ++result;
        }
    }
    return result;
}

std::set<unsigned> bad_switch_masks(
    int k,
    const std::array<int, 7>& labels
) {
    std::set<unsigned> result;
    for (unsigned candidate = 0U;
         candidate < (1U << 7U);
         ++candidate) {
        if (std::popcount(candidate) != 3
            || candidate == selected_mask
            || (std::popcount(candidate & minus_mask) & 1) == 0) {
            continue;
        }
        if (triple_cut_weight(k, labels, candidate) == 1
            && switched_weight_is_zero(k, labels, candidate)) {
            result.insert(candidate);
        }
    }
    return result;
}

int selected_local(
    int k,
    const std::array<int, 7>& labels
) {
    int result = 0;
    for (int output = 0; output <= k; ++output) {
        const int selected = triple_output(
            k, labels[0], labels[1], labels[2], output
        );
        const int reflected_complement = triple_output(
            k,
            k - labels[4],
            labels[5],
            labels[6],
            output
        );
        result += selected * reflected_complement;
    }
    return result;
}

int selected_endpoint_local(
    int k,
    const std::array<int, 7>& labels
) {
    int result = 0;
    for (int offset = 0; offset <= 8; offset += 2) {
        if (offset <= k) {
            result += triple_output(
                k, labels[0], labels[1], labels[2], offset
            ) * triple_output(
                k,
                k - labels[4],
                labels[5],
                labels[6],
                offset
            );
        }
        const int high = k - offset;
        if (high > 8) {
            result += triple_output(
                k, labels[0], labels[1], labels[2], high
            ) * triple_output(
                k,
                k - labels[4],
                labels[5],
                labels[6],
                high
            );
        }
    }
    return result;
}

bool selected_is_maximal_rank_one(
    int k,
    const std::array<int, 7>& labels
) {
    if (triple_cut_weight(k, labels, selected_mask) != 1) {
        return false;
    }
    for (unsigned candidate = 0U;
         candidate < (1U << 7U);
         ++candidate) {
        if (std::popcount(candidate) != 3
            || (std::popcount(candidate & minus_mask) & 1) == 0) {
            continue;
        }
        if (triple_cut_weight(k, labels, candidate) > 1) {
            return false;
        }
    }
    return true;
}

std::set<unsigned> bad_switch_masks_for(
    int k,
    const std::array<int, 7>& labels,
    unsigned chosen_selected_mask
) {
    std::set<unsigned> result;
    for (unsigned candidate = 0U;
         candidate < (1U << 7U);
         ++candidate) {
        if (std::popcount(candidate) != 3
            || candidate == chosen_selected_mask
            || (std::popcount(candidate & minus_mask) & 1) == 0) {
            continue;
        }
        const unsigned difference =
            chosen_selected_mask ^ candidate;
        bool switched_zero = false;
        const int difference_size = std::popcount(difference);
        if (difference_size == 6) {
            switched_zero = true;
        } else if (difference_size == 4) {
            switched_zero = triple_cut_weight(
                k,
                labels,
                (~difference) & 0b1111111U
            ) == 0;
        } else if (difference_size == 2) {
            int first = -1;
            int second = -1;
            for (int position = 0; position < 7; ++position) {
                if (((difference >> position) & 1U) == 0U) {
                    continue;
                }
                if (first < 0) {
                    first = position;
                } else {
                    second = position;
                }
            }
            switched_zero =
                labels[static_cast<std::size_t>(first)]
                != labels[static_cast<std::size_t>(second)];
        }
        if (triple_cut_weight(k, labels, candidate) == 1
            && switched_zero) {
            result.insert(candidate);
        }
    }
    return result;
}

bool chosen_is_maximal_rank_one(
    int k,
    const std::array<int, 7>& labels,
    unsigned chosen_selected_mask
) {
    if (triple_cut_weight(k, labels, chosen_selected_mask) != 1) {
        return false;
    }
    for (unsigned candidate = 0U;
         candidate < (1U << 7U);
         ++candidate) {
        if (std::popcount(candidate) != 3
            || (std::popcount(candidate & minus_mask) & 1) == 0) {
            continue;
        }
        if (triple_cut_weight(k, labels, candidate) > 1) {
            return false;
        }
    }
    return true;
}

int fourfold_output(
    int k,
    int first,
    int second,
    int third,
    int fourth,
    int output
) {
    int result = 0;
    for (int left = 0; left <= k; ++left) {
        if (!fusion(k, first, second, left)) {
            continue;
        }
        for (int right = 0; right <= k; ++right) {
            if (fusion(k, third, fourth, right)
                && fusion(k, left, right, output)) {
                ++result;
            }
        }
    }
    return result;
}

int selected_fundamental_local(
    int k,
    const std::array<int, 7>& labels
) {
    int result = 0;
    for (int offset = 0; offset <= 8; offset += 2) {
        if (offset <= k) {
            result += triple_output(
                k, labels[0], labels[1], labels[2], offset
            ) * fourfold_output(
                k,
                labels[3],
                labels[4],
                labels[5],
                labels[6],
                offset
            );
        }
        const int high = k - offset;
        if (high > 8) {
            result += triple_output(
                k, labels[0], labels[1], labels[2], high
            ) * fourfold_output(
                k,
                labels[3],
                labels[4],
                labels[5],
                labels[6],
                high
            );
        }
    }
    return result;
}

int selected_one_minus_two_plus_complement_top_local(
    int k,
    const std::array<int, 7>& labels
) {
    int result = 0;
    for (int offset = 0; offset <= 8; offset += 2) {
        if (offset <= k) {
            result += triple_output(
                k, labels[0], labels[4], labels[5], offset
            ) * triple_output(
                k,
                k - labels[2],
                labels[3],
                labels[6],
                offset
            );
        }
        const int high = k - offset;
        if (high > 8) {
            result += triple_output(
                k, labels[0], labels[4], labels[5], high
            ) * triple_output(
                k,
                k - labels[2],
                labels[3],
                labels[6],
                high
            );
        }
    }
    return result;
}

int selected_one_minus_two_plus_reflected_local(
    int k,
    const std::array<int, 7>& labels
) {
    int result = 0;
    for (int offset = 0; offset <= 8; offset += 2) {
        if (offset <= k) {
            result += triple_output(
                k, labels[0], labels[4], labels[5], offset
            ) * fourfold_output(
                k,
                labels[1],
                labels[2],
                labels[3],
                labels[6],
                offset
            );
        }
        const int high = k - offset;
        if (high > 8) {
            result += triple_output(
                k, labels[0], labels[4], labels[5], high
            ) * fourfold_output(
                k,
                labels[1],
                labels[2],
                labels[3],
                labels[6],
                high
            );
        }
    }
    return result;
}

}  // namespace

int main() {
    int maximum_bad = -1;
    std::array<int, 7> maximum_labels{};
    int maximum_local = -1;
    int cases = 0;
    std::array<int, 3> family_maximum{-1, -1, -1};
    std::array<std::set<unsigned>, 3> family_bad_masks{};
    for (int k = 8; k <= 200; k += 2) {
        for (int a = 4; 2 * a <= k; a += 2) {
            const std::array<std::array<int, 3>, 3> plus_families{{
                {1, k - 1, 2},
                {1, 1, k - 2},
                {k - 1, k - 1, k - 2},
            }};
            for (std::size_t family = 0U;
                 family < plus_families.size();
                 ++family) {
                const auto& plus = plus_families[family];
                const std::array<int, 7> labels{
                    a, k - a, k, k,
                    plus[0], plus[1], plus[2]
                };
                bool disjoint = true;
                for (int minus = 0; minus < 4; ++minus) {
                    for (int positive = 4; positive < 7; ++positive) {
                        disjoint = disjoint
                            && labels[static_cast<std::size_t>(minus)]
                                != labels[
                                    static_cast<std::size_t>(positive)
                                ];
                    }
                }
                if (!disjoint) {
                    continue;
                }
                ++cases;
                const int local = selected_local(k, labels);
                const int bad = bad_switch_count(k, labels);
                const std::set<unsigned> bad_masks =
                    bad_switch_masks(k, labels);
                family_maximum[family] = std::max(
                    family_maximum[family], bad
                );
                family_bad_masks[family].insert(
                    bad_masks.begin(), bad_masks.end()
                );
                if (local != 4) {
                    std::cerr << "unexpected local=" << local << '\n';
                    return 1;
                }
                if (bad > maximum_bad) {
                    maximum_bad = bad;
                    maximum_local = local;
                    maximum_labels = labels;
                }
            }
        }
    }
    std::cout
        << "SU2_THREE_MINUS_COMPLEMENT_TOP_L4 cases=" << cases
        << " maximum_bad=" << maximum_bad
        << " local=" << maximum_local
        << " labels=";
    for (std::size_t index = 0U;
         index < maximum_labels.size();
         ++index) {
        std::cout << (index == 0U ? "" : ",") << maximum_labels[index];
    }
    std::cout << " result="
        << (maximum_bad <= 3 ? "PASS" : "FAIL")
        << '\n';
    for (std::size_t family = 0U;
         family < family_maximum.size();
         ++family) {
        std::cout << "family=" << family
            << " maximum_bad=" << family_maximum[family]
            << " possible_bad_masks=";
        bool first = true;
        for (const unsigned mask : family_bad_masks[family]) {
            std::cout << (first ? "" : ",") << mask;
            first = false;
        }
        std::cout << '\n';
    }

    std::array<int, 16> local_maximum_bad{};
    local_maximum_bad.fill(-1);
    std::array<std::array<int, 7>, 16> local_witness{};
    std::set<unsigned> general_bad_masks;
    std::set<std::vector<unsigned>> general_bad_sets;
    int overall_maximum_bad = -1;
    int overall_maximum_local = -1;
    std::array<int, 7> overall_witness{};
    std::set<unsigned> overall_bad_masks;
    std::int64_t general_cases = 0;
    for (int k = 4; k <= 30; k += 2) {
        for (int first = 2; first <= k; first += 2) {
            for (int second = first;
                 second <= k;
                 second += 2) {
                for (int third = second;
                     third <= k;
                     third += 2) {
                    for (int first_plus = 1;
                         first_plus <= k;
                         first_plus += 2) {
                        for (int second_plus = first_plus;
                             second_plus <= k;
                             second_plus += 2) {
                            for (int even_plus = 2;
                                 even_plus <= k;
                                 even_plus += 2) {
                                const std::array<int, 7> labels{
                                    first,
                                    second,
                                    third,
                                    k,
                                    first_plus,
                                    second_plus,
                                    even_plus
                                };
                                bool disjoint = true;
                                for (int minus = 0;
                                     minus < 4;
                                     ++minus) {
                                    for (int positive = 4;
                                         positive < 7;
                                         ++positive) {
                                        disjoint = disjoint
                                            && labels[
                                                static_cast<std::size_t>(
                                                    minus
                                                )
                                            ] != labels[
                                                static_cast<std::size_t>(
                                                    positive
                                                )
                                            ];
                                    }
                                }
                                if (!disjoint
                                    || !selected_is_maximal_rank_one(
                                        k, labels
                                    )) {
                                    continue;
                                }
                                ++general_cases;
                                const int local =
                                    selected_endpoint_local(k, labels);
                                const int bad =
                                    bad_switch_count(k, labels);
                                if (bad > overall_maximum_bad) {
                                    overall_maximum_bad = bad;
                                    overall_maximum_local = local;
                                    overall_witness = labels;
                                    overall_bad_masks =
                                        bad_switch_masks(k, labels);
                                }
                                if (local < 0) {
                                    std::cerr
                                        << "unexpected endpoint local="
                                        << local << '\n';
                                    return 1;
                                }
                                if (local >= 16) {
                                    continue;
                                }
                                const std::set<unsigned> bad_masks =
                                    bad_switch_masks(k, labels);
                                general_bad_masks.insert(
                                    bad_masks.begin(), bad_masks.end()
                                );
                                general_bad_sets.emplace(
                                    bad_masks.begin(), bad_masks.end()
                                );
                                if (bad > local_maximum_bad[
                                        static_cast<std::size_t>(local)
                                    ]) {
                                    local_maximum_bad[
                                        static_cast<std::size_t>(local)
                                    ] = bad;
                                    local_witness[
                                        static_cast<std::size_t>(local)
                                    ] = labels;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout
        << "SU2_THREE_MINUS_COMPLEMENT_TOP_GENERAL levels=4..30"
        << " cases=" << general_cases
        << " overall_maximum_bad=" << overall_maximum_bad
        << " overall_local=" << overall_maximum_local
        << " overall_labels=";
    for (std::size_t index = 0U;
         index < overall_witness.size();
         ++index) {
        std::cout << (index == 0U ? "" : ",")
            << overall_witness[index];
    }
    std::cout << " overall_masks=";
    bool first_overall = true;
    for (const unsigned mask : overall_bad_masks) {
        std::cout << (first_overall ? "" : ",") << mask;
        first_overall = false;
    }
    std::cout << '\n';
    for (std::size_t local = 0U;
         local < local_maximum_bad.size();
         ++local) {
        if (local_maximum_bad[local] < 0) {
            continue;
        }
        std::cout << "local=" << local
            << " maximum_bad=" << local_maximum_bad[local]
            << " labels=";
        for (std::size_t index = 0U;
             index < local_witness[local].size();
             ++index) {
            std::cout << (index == 0U ? "" : ",")
                << local_witness[local][index];
        }
        std::cout << '\n';
    }
    std::cout << "general_possible_bad_masks=";
    bool first_mask = true;
    for (const unsigned mask : general_bad_masks) {
        std::cout << (first_mask ? "" : ",") << mask;
        first_mask = false;
    }
    std::cout << '\n';
    std::cout << "general_bad_sets=";
    bool first_set = true;
    for (const auto& masks : general_bad_sets) {
        std::cout << (first_set ? "" : ";") << "{";
        for (std::size_t index = 0U;
             index < masks.size();
             ++index) {
            std::cout << (index == 0U ? "" : ",") << masks[index];
        }
        std::cout << "}";
        first_set = false;
    }
    std::cout << '\n';

    std::array<int, 16> fundamental_local_maximum_bad{};
    fundamental_local_maximum_bad.fill(-1);
    std::array<std::array<int, 7>, 16> fundamental_witness{};
    std::set<std::vector<unsigned>> fundamental_bad_sets;
    std::int64_t fundamental_cases = 0;
    int fundamental_overall_maximum_bad = -1;
    int fundamental_overall_local = -1;
    std::array<int, 7> fundamental_overall_witness{};
    for (int k = 4; k <= 24; k += 2) {
        for (int first = 2; first <= k; first += 2) {
            for (int second = first;
                 second <= k;
                 second += 2) {
                for (int third = second;
                     third <= k;
                     third += 2) {
                    for (int fourth = 2;
                         fourth <= k;
                         fourth += 2) {
                        for (int odd_plus = 1;
                             odd_plus <= k;
                             odd_plus += 2) {
                            for (int even_plus = 2;
                                 even_plus <= k;
                                 even_plus += 2) {
                                if (even_plus == first
                                    || even_plus == second
                                    || even_plus == third
                                    || even_plus == fourth) {
                                    continue;
                                }
                                const std::array<int, 7> labels{
                                    first,
                                    second,
                                    third,
                                    fourth,
                                    1,
                                    odd_plus,
                                    even_plus
                                };
                                if (!selected_is_maximal_rank_one(
                                        k, labels
                                    )) {
                                    continue;
                                }
                                ++fundamental_cases;
                                const int local =
                                    selected_fundamental_local(
                                        k, labels
                                    );
                                const int bad =
                                    bad_switch_count(k, labels);
                                if (bad
                                    > fundamental_overall_maximum_bad) {
                                    fundamental_overall_maximum_bad = bad;
                                    fundamental_overall_local = local;
                                    fundamental_overall_witness = labels;
                                }
                                if (local < 0) {
                                    std::cerr
                                        << "unexpected fundamental local="
                                        << local << '\n';
                                    return 1;
                                }
                                if (local >= 16) {
                                    continue;
                                }
                                const std::set<unsigned> bad_masks =
                                    bad_switch_masks(k, labels);
                                fundamental_bad_sets.emplace(
                                    bad_masks.begin(), bad_masks.end()
                                );
                                if (
                                    bad > fundamental_local_maximum_bad[
                                        static_cast<std::size_t>(local)
                                    ]
                                ) {
                                    fundamental_local_maximum_bad[
                                        static_cast<std::size_t>(local)
                                    ] = bad;
                                    fundamental_witness[
                                        static_cast<std::size_t>(local)
                                    ] = labels;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout
        << "SU2_THREE_MINUS_COMPLEMENT_FUNDAMENTAL levels=4..24"
        << " cases=" << fundamental_cases
        << " overall_maximum_bad="
        << fundamental_overall_maximum_bad
        << " overall_local=" << fundamental_overall_local
        << " overall_labels=";
    for (std::size_t index = 0U;
         index < fundamental_overall_witness.size();
         ++index) {
        std::cout << (index == 0U ? "" : ",")
            << fundamental_overall_witness[index];
    }
    std::cout << '\n';
    for (std::size_t local = 0U;
         local < fundamental_local_maximum_bad.size();
         ++local) {
        if (fundamental_local_maximum_bad[local] < 0) {
            continue;
        }
        std::cout << "fundamental_local=" << local
            << " maximum_bad="
            << fundamental_local_maximum_bad[local]
            << " labels=";
        for (std::size_t index = 0U;
             index < fundamental_witness[local].size();
             ++index) {
            std::cout << (index == 0U ? "" : ",")
                << fundamental_witness[local][index];
        }
        std::cout << '\n';
    }
    std::cout << "fundamental_bad_sets=";
    first_set = true;
    for (const auto& masks : fundamental_bad_sets) {
        std::cout << (first_set ? "" : ";") << "{";
        for (std::size_t index = 0U;
             index < masks.size();
             ++index) {
            std::cout << (index == 0U ? "" : ",") << masks[index];
        }
        std::cout << "}";
        first_set = false;
    }
    std::cout << '\n';

    constexpr unsigned one_minus_two_plus_selected = 0b0110001U;
    std::array<int, 16> corrected_top_local_maximum_bad{};
    corrected_top_local_maximum_bad.fill(-1);
    std::array<std::array<int, 7>, 16> corrected_top_witness{};
    std::set<std::vector<unsigned>> corrected_top_bad_sets;
    std::int64_t corrected_top_cases = 0;
    for (int k = 4; k <= 24; k += 2) {
        for (int selected_minus = 2;
             selected_minus <= k;
             selected_minus += 2) {
            for (int second_minus = 2;
                 second_minus <= k;
                 second_minus += 2) {
                for (int third_minus = second_minus;
                     third_minus <= k;
                     third_minus += 2) {
                    for (int first_plus = 1;
                         first_plus <= k;
                         first_plus += 2) {
                        for (int second_plus = first_plus;
                             second_plus <= k;
                             second_plus += 2) {
                            for (int even_plus = 2;
                                 even_plus <= k;
                                 even_plus += 2) {
                                if (even_plus == selected_minus
                                    || even_plus == k
                                    || even_plus == second_minus
                                    || even_plus == third_minus) {
                                    continue;
                                }
                                const std::array<int, 7> labels{
                                    selected_minus,
                                    k,
                                    second_minus,
                                    third_minus,
                                    first_plus,
                                    second_plus,
                                    even_plus
                                };
                                if (!chosen_is_maximal_rank_one(
                                        k,
                                        labels,
                                        one_minus_two_plus_selected
                                    )) {
                                    continue;
                                }
                                ++corrected_top_cases;
                                const int local =
                                    selected_one_minus_two_plus_complement_top_local(
                                        k, labels
                                    );
                                if (local < 0) {
                                    std::cerr
                                        << "unexpected corrected top local="
                                        << local << '\n';
                                    return 1;
                                }
                                if (local >= 16) {
                                    continue;
                                }
                                const std::set<unsigned> bad_masks =
                                    bad_switch_masks_for(
                                        k,
                                        labels,
                                        one_minus_two_plus_selected
                                    );
                                const int bad =
                                    static_cast<int>(bad_masks.size());
                                corrected_top_bad_sets.emplace(
                                    bad_masks.begin(), bad_masks.end()
                                );
                                if (
                                    bad > corrected_top_local_maximum_bad[
                                        static_cast<std::size_t>(local)
                                    ]
                                ) {
                                    corrected_top_local_maximum_bad[
                                        static_cast<std::size_t>(local)
                                    ] = bad;
                                    corrected_top_witness[
                                        static_cast<std::size_t>(local)
                                    ] = labels;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout
        << "SU2_CORRECTED_ONE_MINUS_TWO_PLUS_COMPLEMENT_TOP"
        << " levels=4..24 cases=" << corrected_top_cases << '\n';
    for (std::size_t local = 0U;
         local < corrected_top_local_maximum_bad.size();
         ++local) {
        if (corrected_top_local_maximum_bad[local] < 0) {
            continue;
        }
        std::cout << "corrected_top_local=" << local
            << " maximum_bad="
            << corrected_top_local_maximum_bad[local]
            << " labels=";
        for (std::size_t index = 0U;
             index < corrected_top_witness[local].size();
             ++index) {
            std::cout << (index == 0U ? "" : ",")
                << corrected_top_witness[local][index];
        }
        std::cout << '\n';
    }
    std::cout << "corrected_top_bad_sets=";
    first_set = true;
    for (const auto& masks : corrected_top_bad_sets) {
        std::cout << (first_set ? "" : ";") << "{";
        for (std::size_t index = 0U;
             index < masks.size();
             ++index) {
            std::cout << (index == 0U ? "" : ",") << masks[index];
        }
        std::cout << "}";
        first_set = false;
    }
    std::cout << '\n';

    std::array<int, 16> corrected_reflected_local_maximum_bad{};
    corrected_reflected_local_maximum_bad.fill(-1);
    std::array<std::array<int, 7>, 16> corrected_reflected_witness{};
    std::set<std::vector<unsigned>> corrected_reflected_bad_sets;
    std::int64_t corrected_reflected_cases = 0;
    for (int k = 4; k <= 22; k += 2) {
        for (int selected_minus = 2;
             selected_minus <= k;
             selected_minus += 2) {
            for (int first_minus = 2;
                 first_minus <= k;
                 first_minus += 2) {
                for (int second_minus = first_minus;
                     second_minus <= k;
                     second_minus += 2) {
                    for (int third_minus = second_minus;
                         third_minus <= k;
                         third_minus += 2) {
                        for (int second_plus = 1;
                             second_plus <= k;
                             second_plus += 2) {
                            for (int even_plus = 2;
                                 even_plus <= k;
                                 even_plus += 2) {
                                if (even_plus == selected_minus
                                    || even_plus == first_minus
                                    || even_plus == second_minus
                                    || even_plus == third_minus) {
                                    continue;
                                }
                                const std::array<int, 7> labels{
                                    selected_minus,
                                    first_minus,
                                    second_minus,
                                    third_minus,
                                    k - 1,
                                    second_plus,
                                    even_plus
                                };
                                if (!chosen_is_maximal_rank_one(
                                        k,
                                        labels,
                                        one_minus_two_plus_selected
                                    )) {
                                    continue;
                                }
                                ++corrected_reflected_cases;
                                const int local =
                                    selected_one_minus_two_plus_reflected_local(
                                        k, labels
                                    );
                                if (local < 0) {
                                    std::cerr
                                        << "unexpected corrected reflected local="
                                        << local << '\n';
                                    return 1;
                                }
                                if (local >= 16) {
                                    continue;
                                }
                                const std::set<unsigned> bad_masks =
                                    bad_switch_masks_for(
                                        k,
                                        labels,
                                        one_minus_two_plus_selected
                                    );
                                const int bad =
                                    static_cast<int>(bad_masks.size());
                                corrected_reflected_bad_sets.emplace(
                                    bad_masks.begin(), bad_masks.end()
                                );
                                if (
                                    bad
                                    > corrected_reflected_local_maximum_bad[
                                        static_cast<std::size_t>(local)
                                    ]
                                ) {
                                    corrected_reflected_local_maximum_bad[
                                        static_cast<std::size_t>(local)
                                    ] = bad;
                                    corrected_reflected_witness[
                                        static_cast<std::size_t>(local)
                                    ] = labels;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout
        << "SU2_CORRECTED_SELECTED_REFLECTED_FUNDAMENTAL"
        << " levels=4..22 cases=" << corrected_reflected_cases
        << '\n';
    for (std::size_t local = 0U;
         local < corrected_reflected_local_maximum_bad.size();
         ++local) {
        if (corrected_reflected_local_maximum_bad[local] < 0) {
            continue;
        }
        std::cout << "corrected_reflected_local=" << local
            << " maximum_bad="
            << corrected_reflected_local_maximum_bad[local]
            << " labels=";
        for (std::size_t index = 0U;
             index < corrected_reflected_witness[local].size();
             ++index) {
            std::cout << (index == 0U ? "" : ",")
                << corrected_reflected_witness[local][index];
        }
        std::cout << '\n';
    }
    std::cout << "corrected_reflected_bad_sets=";
    first_set = true;
    for (const auto& masks : corrected_reflected_bad_sets) {
        std::cout << (first_set ? "" : ";") << "{";
        for (std::size_t index = 0U;
             index < masks.size();
             ++index) {
            std::cout << (index == 0U ? "" : ",") << masks[index];
        }
        std::cout << "}";
        first_set = false;
    }
    std::cout << '\n';

    std::array<int, 20> corrected_neighbor_local_maximum_bad{};
    corrected_neighbor_local_maximum_bad.fill(-1);
    std::array<std::array<int, 7>, 20> corrected_neighbor_witness{};
    std::set<std::vector<unsigned>> corrected_neighbor_bad_sets;
    std::int64_t corrected_neighbor_cases = 0;
    for (int k = 3; k <= 23; k += 2) {
        for (int selected_minus = 2;
             selected_minus <= k;
             selected_minus += 2) {
            for (int second_minus = 2;
                 second_minus <= k;
                 second_minus += 2) {
                for (int third_minus = second_minus;
                     third_minus <= k;
                     third_minus += 2) {
                    for (int first_plus = 1;
                         first_plus <= k;
                         first_plus += 2) {
                        for (int second_plus = first_plus;
                             second_plus <= k;
                             second_plus += 2) {
                            for (int even_plus = 2;
                                 even_plus <= k;
                                 even_plus += 2) {
                                if (even_plus == selected_minus
                                    || even_plus == k - 1
                                    || even_plus == second_minus
                                    || even_plus == third_minus) {
                                    continue;
                                }
                                const std::array<int, 7> labels{
                                    selected_minus,
                                    k - 1,
                                    second_minus,
                                    third_minus,
                                    first_plus,
                                    second_plus,
                                    even_plus
                                };
                                if (!chosen_is_maximal_rank_one(
                                        k,
                                        labels,
                                        one_minus_two_plus_selected
                                    )) {
                                    continue;
                                }
                                ++corrected_neighbor_cases;
                                const int local =
                                    selected_one_minus_two_plus_reflected_local(
                                        k, labels
                                    );
                                if (local < 0) {
                                    std::cerr
                                        << "unexpected corrected neighbor local="
                                        << local << '\n';
                                    return 1;
                                }
                                if (
                                    local
                                    >= static_cast<int>(
                                        corrected_neighbor_local_maximum_bad
                                            .size()
                                    )
                                ) {
                                    continue;
                                }
                                const std::set<unsigned> bad_masks =
                                    bad_switch_masks_for(
                                        k,
                                        labels,
                                        one_minus_two_plus_selected
                                    );
                                const int bad =
                                    static_cast<int>(bad_masks.size());
                                corrected_neighbor_bad_sets.emplace(
                                    bad_masks.begin(), bad_masks.end()
                                );
                                if (
                                    bad
                                    > corrected_neighbor_local_maximum_bad[
                                        static_cast<std::size_t>(local)
                                    ]
                                ) {
                                    corrected_neighbor_local_maximum_bad[
                                        static_cast<std::size_t>(local)
                                    ] = bad;
                                    corrected_neighbor_witness[
                                        static_cast<std::size_t>(local)
                                    ] = labels;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout
        << "SU2_CORRECTED_COMPLEMENT_REFLECTED_NEIGHBOR"
        << " levels=3..23 cases=" << corrected_neighbor_cases
        << '\n';
    for (std::size_t local = 0U;
         local < corrected_neighbor_local_maximum_bad.size();
         ++local) {
        if (corrected_neighbor_local_maximum_bad[local] < 0) {
            continue;
        }
        std::cout << "corrected_neighbor_local=" << local
            << " maximum_bad="
            << corrected_neighbor_local_maximum_bad[local]
            << " labels=";
        for (std::size_t index = 0U;
             index < corrected_neighbor_witness[local].size();
             ++index) {
            std::cout << (index == 0U ? "" : ",")
                << corrected_neighbor_witness[local][index];
        }
        std::cout << '\n';
    }
    std::cout << "corrected_neighbor_bad_sets=";
    first_set = true;
    for (const auto& masks : corrected_neighbor_bad_sets) {
        std::cout << (first_set ? "" : ";") << "{";
        for (std::size_t index = 0U;
             index < masks.size();
             ++index) {
            std::cout << (index == 0U ? "" : ",") << masks[index];
        }
        std::cout << "}";
        first_set = false;
    }
    std::cout << '\n';
    return maximum_bad <= 3 ? 0 : 1;
}
