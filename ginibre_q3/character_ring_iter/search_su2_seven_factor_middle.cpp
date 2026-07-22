#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Witness {
    bool initialized = false;
    std::int64_t value = 0;
    std::vector<int> minus;
    std::vector<int> plus;
    std::int64_t invariant = 0;
    std::int64_t negative_middle = 0;
    std::int64_t positive_middle = 0;
    std::int64_t pair_term = 0;
};

std::int64_t invariant_multiplicity(const std::vector<int>& labels) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    std::vector<std::int64_t> current(static_cast<std::size_t>(total + 1), 0);
    std::vector<std::int64_t> next(static_cast<std::size_t>(total + 1), 0);
    current[0] = 1;
    int current_maximum = 0;
    for (int label : labels) {
        std::fill(next.begin(), next.end(), 0);
        for (int source = 0; source <= current_maximum; ++source) {
            const std::int64_t multiplicity = current[static_cast<std::size_t>(source)];
            if (multiplicity == 0) {
                continue;
            }
            for (int output = std::abs(source - label);
                 output <= source + label; output += 2) {
                next[static_cast<std::size_t>(output)] += multiplicity;
            }
        }
        current_maximum += label;
        current.swap(next);
    }
    return current[0];
}

void combinations_rec(
    int maximum_label,
    int remaining,
    int first,
    std::vector<int>& current,
    std::vector<std::vector<int>>& output
) {
    if (remaining == 0) {
        output.push_back(current);
        return;
    }
    for (int label = first; label <= maximum_label; ++label) {
        current.push_back(label);
        combinations_rec(maximum_label, remaining - 1, label, current, output);
        current.pop_back();
    }
}

std::vector<std::vector<int>> combinations(int maximum_label, int size) {
    std::vector<std::vector<int>> output;
    std::vector<int> current;
    combinations_rec(maximum_label, size, 1, current, output);
    return output;
}

bool supports_are_disjoint(
    const std::vector<int>& minus,
    const std::vector<int>& plus
) {
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < minus.size() && j < plus.size()) {
        if (minus[i] == plus[j]) {
            return false;
        }
        if (minus[i] < plus[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    return true;
}

bool all_labels_are_distinct(const std::vector<int>& labels) {
    std::vector<int> sorted = labels;
    std::sort(sorted.begin(), sorted.end());
    return std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
}

int popcount(unsigned int value) {
    int count = 0;
    while (value != 0U) {
        count += static_cast<int>(value & 1U);
        value >>= 1U;
    }
    return count;
}

std::vector<int> selected(
    const std::vector<int>& labels,
    unsigned int mask,
    bool take_mask
) {
    std::vector<int> result;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        const bool present = ((mask >> i) & 1U) != 0U;
        if (present == take_mask) {
            result.push_back(labels[i]);
        }
    }
    return result;
}

std::int64_t q3_subset(
    const std::vector<int>& labels,
    std::size_t minus_count
) {
    std::int64_t answer = 0;
    constexpr unsigned int full_mask = (1U << 7U) - 1U;
    for (unsigned int mask = 0; mask <= full_mask; ++mask) {
        const std::int64_t term =
            invariant_multiplicity(selected(labels, mask, true))
            * invariant_multiplicity(selected(labels, mask, false));
        int parity = 0;
        for (std::size_t i = 0; i < minus_count; ++i) {
            parity ^= static_cast<int>((mask >> i) & 1U);
        }
        answer += parity == 0 ? term : -term;
    }
    return answer;
}

void consider(
    Witness& witness,
    std::int64_t value,
    const std::vector<int>& minus,
    const std::vector<int>& plus,
    std::int64_t invariant,
    std::int64_t negative_middle,
    std::int64_t positive_middle,
    std::int64_t pair_term
) {
    if (!witness.initialized || value < witness.value) {
        witness = Witness{
            true, value, minus, plus, invariant,
            negative_middle, positive_middle, pair_term
        };
    }
}

void print_vector(const std::vector<int>& values) {
    std::cout << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0U) {
            std::cout << ',';
        }
        std::cout << values[i];
    }
    std::cout << ']';
}

void print_witness(const std::string& name, const Witness& witness) {
    std::cout << name << '=' << witness.value << " minus=";
    print_vector(witness.minus);
    std::cout << " plus=";
    print_vector(witness.plus);
    std::cout << " invariant=" << witness.invariant
              << " negative_middle=" << witness.negative_middle
              << " positive_middle=" << witness.positive_middle
              << " pair_term=" << witness.pair_term << '\n';
}

int parse_positive(const char* text) {
    const long value = std::strtol(text, nullptr, 10);
    if (value <= 0 || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error("maximum label must be positive");
    }
    return static_cast<int>(value);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: search_su2_seven_factor_middle maximum_label"
            );
        }
        const int maximum_label = parse_positive(argv[1]);
        Witness raw_domination;
        Witness signed_middle_domination;
        Witness exact_half_witness;
        Witness signed_distinct_witness;
        Witness pair_domination_witness;
        Witness active_exact_witness;
        std::uint64_t cases = 0;

        for (int minus_count = 2; minus_count <= 6; minus_count += 2) {
            const auto minus_lists = combinations(maximum_label, minus_count);
            const auto plus_lists = combinations(maximum_label, 7 - minus_count);
            for (const auto& minus : minus_lists) {
                for (const auto& plus : plus_lists) {
                    if (!supports_are_disjoint(minus, plus)) {
                        continue;
                    }
                    ++cases;
                    std::vector<int> labels = minus;
                    labels.insert(labels.end(), plus.begin(), plus.end());
                    const std::int64_t invariant = invariant_multiplicity(labels);

                    std::int64_t pair_term = 0;
                    for (std::size_t i = 0; i < labels.size(); ++i) {
                        for (std::size_t j = i + 1U; j < labels.size(); ++j) {
                            if (labels[i] != labels[j]) {
                                continue;
                            }
                            std::vector<int> rest;
                            for (std::size_t t = 0; t < labels.size(); ++t) {
                                if (t != i && t != j) {
                                    rest.push_back(labels[t]);
                                }
                            }
                            const bool opposite_signs =
                                (i < minus.size()) != (j < minus.size());
                            const std::int64_t term = invariant_multiplicity(rest);
                            pair_term += opposite_signs ? -term : term;
                        }
                    }

                    std::int64_t negative_middle = 0;
                    std::int64_t positive_middle = 0;
                    std::int64_t all_middle = 0;
                    constexpr unsigned int full_mask = (1U << 7U) - 1U;
                    for (unsigned int mask = 0; mask <= full_mask; ++mask) {
                        if (popcount(mask) != 3) {
                            continue;
                        }
                        const std::int64_t term =
                            invariant_multiplicity(selected(labels, mask, true))
                            * invariant_multiplicity(selected(labels, mask, false));
                        all_middle += term;
                        int parity = 0;
                        for (std::size_t i = 0; i < minus.size(); ++i) {
                            parity ^= static_cast<int>((mask >> i) & 1U);
                        }
                        if (parity == 0) {
                            positive_middle += term;
                        } else {
                            negative_middle += term;
                        }
                    }

                    const std::int64_t exact_half = invariant + pair_term
                        + positive_middle - negative_middle;
                    if (cases <= 1024U
                        && 2 * exact_half != q3_subset(labels, minus.size())) {
                        throw std::runtime_error(
                            "seven-factor middle decomposition disagrees with subset sum"
                        );
                    }
                    consider(
                        raw_domination, invariant - all_middle, minus, plus,
                        invariant, negative_middle, positive_middle, pair_term
                    );
                    consider(
                        signed_middle_domination, invariant - negative_middle,
                        minus, plus, invariant, negative_middle,
                        positive_middle, pair_term
                    );
                    if (all_labels_are_distinct(labels)) {
                        consider(
                            signed_distinct_witness,
                            invariant - negative_middle, minus, plus,
                            invariant, negative_middle,
                            positive_middle, pair_term
                        );
                    }
                    consider(
                        pair_domination_witness,
                        invariant + pair_term - negative_middle,
                        minus, plus, invariant, negative_middle,
                        positive_middle, pair_term
                    );
                    consider(
                        exact_half_witness, exact_half, minus, plus,
                        invariant, negative_middle, positive_middle, pair_term
                    );
                    if (invariant != 0 || negative_middle != 0
                        || positive_middle != 0 || pair_term != 0) {
                        consider(
                            active_exact_witness, exact_half, minus, plus,
                            invariant, negative_middle,
                            positive_middle, pair_term
                        );
                    }
                }
            }
        }

        std::cout << "SU2_SEVEN_FACTOR_MIDDLE cases=" << cases
                  << " maximum_label=" << maximum_label << '\n';
        print_witness("raw_domination_minimum", raw_domination);
        print_witness("signed_domination_minimum", signed_middle_domination);
        print_witness("signed_distinct_minimum", signed_distinct_witness);
        print_witness("pair_domination_minimum", pair_domination_witness);
        print_witness("exact_half_minimum", exact_half_witness);
        print_witness("active_exact_half_minimum", active_exact_witness);
        std::cout << "SU2_SEVEN_FACTOR_MIDDLE "
                  << (exact_half_witness.value < 0 ? "FAIL" : "PASS") << '\n';
        return exact_half_witness.value < 0 ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SEVEN_FACTOR_MIDDLE FAILURE: " << error.what() << '\n';
        return 1;
    }
}
