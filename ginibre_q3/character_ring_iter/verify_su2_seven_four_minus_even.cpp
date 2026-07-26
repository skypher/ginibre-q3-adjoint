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

using Labels = std::vector<int>;

struct LocalWitness {
    bool initialized = false;
    std::int64_t margin = 0;
    Labels minus;
    Labels plus;
    int level = -1;
};

struct Totals {
    std::uint64_t cases = 0;
    std::uint64_t active_cases = 0;
    std::uint64_t raw_bound_failures = 0;
    std::uint64_t raw_failures_with_simple_current = 0;
    std::uint64_t raw_failures_without_simple_current = 0;
    std::uint64_t raw_failures_all_cap_at_least_three = 0;
    std::uint64_t finite_deep_large_cut_cases = 0;
    std::int64_t minimum_margin = std::numeric_limits<std::int64_t>::max();
    std::int64_t minimum_no_k_margin =
        std::numeric_limits<std::int64_t>::max();
    std::int64_t minimum_deep_margin =
        std::numeric_limits<std::int64_t>::max();
    std::int64_t minimum_half = std::numeric_limits<std::int64_t>::max();
    Labels minimum_minus;
    Labels minimum_plus;
    int minimum_level = -1;
    Labels minimum_half_minus;
    Labels minimum_half_plus;
    int minimum_half_level = -1;
    Labels minimum_no_k_minus;
    Labels minimum_no_k_plus;
    int minimum_no_k_level = -1;
    Labels minimum_deep_minus;
    Labels minimum_deep_plus;
    int minimum_deep_level = -1;
    std::array<LocalWitness, 6> deep_small_local{};
    std::array<LocalWitness, 6> deep_small_worst_local{};
};

int parse_nonnegative(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value < 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string(name) + " must be nonnegative");
    }
    return static_cast<int>(value);
}

void sorted_lists_rec(
    const Labels& alphabet,
    int remaining,
    std::size_t first,
    Labels& current,
    std::vector<Labels>& output
) {
    if (remaining == 0) {
        output.push_back(current);
        return;
    }
    for (std::size_t i = first; i < alphabet.size(); ++i) {
        current.push_back(alphabet[i]);
        sorted_lists_rec(alphabet, remaining - 1, i, current, output);
        current.pop_back();
    }
}

std::vector<Labels> sorted_lists(const Labels& alphabet, int size) {
    std::vector<Labels> output;
    Labels current;
    sorted_lists_rec(alphabet, size, 0U, current, output);
    return output;
}

bool disjoint_support(const Labels& first, const Labels& second) {
    std::size_t i = 0U;
    std::size_t j = 0U;
    while (i < first.size() && j < second.size()) {
        if (first[i] == second[j]) {
            return false;
        }
        if (first[i] < second[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    return true;
}

std::int64_t multiplicity(const Labels& labels, int output, int level) {
    int bound = level;
    if (level < 0) {
        bound = 0;
        for (const int label : labels) {
            bound += label;
        }
    }
    std::vector<std::int64_t> current(
        static_cast<std::size_t>(bound + 1), 0
    );
    std::vector<std::int64_t> next(
        static_cast<std::size_t>(bound + 1), 0
    );
    current[0] = 1;
    int current_maximum = 0;
    for (const int label : labels) {
        std::fill(next.begin(), next.end(), 0);
        const int source_bound =
            level < 0 ? current_maximum : std::min(current_maximum, level);
        for (int source = 0; source <= source_bound; ++source) {
            const std::int64_t multiplicity =
                current[static_cast<std::size_t>(source)];
            if (multiplicity == 0) {
                continue;
            }
            const int upper = level < 0
                ? source + label
                : std::min(source + label, 2 * level - source - label);
            for (int destination = std::abs(source - label);
                 destination <= upper; destination += 2) {
                next[static_cast<std::size_t>(destination)] += multiplicity;
            }
        }
        current_maximum = level < 0
            ? current_maximum + label
            : std::min(level, current_maximum + label);
        current.swap(next);
    }
    if (output < 0 || output > bound) {
        return 0;
    }
    return current[static_cast<std::size_t>(output)];
}

std::int64_t invariant(const Labels& labels, int level) {
    return multiplicity(labels, 0, level);
}

Labels subset(const Labels& labels, unsigned int mask, bool selected) {
    Labels output;
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        const bool is_selected = ((mask >> i) & 1U) != 0U;
        if (is_selected == selected) {
            output.push_back(labels[i]);
        }
    }
    return output;
}

int popcount(unsigned int value) {
    int result = 0;
    while (value != 0U) {
        result += static_cast<int>(value & 1U);
        value >>= 1U;
    }
    return result;
}

void print_labels(const Labels& labels) {
    std::cout << '[';
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        if (i != 0U) {
            std::cout << ',';
        }
        std::cout << labels[i];
    }
    std::cout << ']';
}

void check_case(
    const Labels& minus,
    const Labels& plus,
    int level,
    bool require_raw_bound,
    Totals& totals
) {
    Labels labels = minus;
    labels.insert(labels.end(), plus.begin(), plus.end());
    const std::int64_t sevenfold = invariant(labels, level);
    std::int64_t negative_middle = 0;
    std::int64_t positive_middle = 0;
    std::int64_t maximum_cut = 0;
    int active_cuts = 0;
    std::vector<unsigned int> maximum_masks;
    constexpr unsigned int full_mask = (1U << 7U) - 1U;
    for (unsigned int mask = 0U; mask <= full_mask; ++mask) {
        if (popcount(mask) != 3) {
            continue;
        }
        int minus_in_cut = 0;
        for (std::size_t i = 0U; i < 4U; ++i) {
            minus_in_cut += static_cast<int>((mask >> i) & 1U);
        }
        const std::int64_t term =
            invariant(subset(labels, mask, true), level)
            * invariant(subset(labels, mask, false), level);
        if ((minus_in_cut % 2) == 0) {
            positive_middle += term;
        } else {
            negative_middle += term;
            if (term > maximum_cut) {
                maximum_cut = term;
                maximum_masks.clear();
                maximum_masks.push_back(mask);
            } else if (term != 0 && term == maximum_cut) {
                maximum_masks.push_back(mask);
            }
            active_cuts += term != 0 ? 1 : 0;
        }
    }

    std::int64_t pair_term = 0;
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        for (std::size_t j = i + 1U; j < labels.size(); ++j) {
            if (labels[i] != labels[j]) {
                continue;
            }
            Labels rest;
            for (std::size_t t = 0U; t < labels.size(); ++t) {
                if (t != i && t != j) {
                    rest.push_back(labels[t]);
                }
            }
            pair_term += invariant(rest, level);
        }
    }
    const std::int64_t exact_half =
        sevenfold + pair_term + positive_middle - negative_middle;

    ++totals.cases;
    if (exact_half < totals.minimum_half) {
        totals.minimum_half = exact_half;
        totals.minimum_half_minus = minus;
        totals.minimum_half_plus = plus;
        totals.minimum_half_level = level;
    }
    if (exact_half < 0) {
        throw std::runtime_error("negative exact half-contraction");
    }
    if (maximum_cut == 0) {
        if (negative_middle != 0) {
            throw std::runtime_error("zero maximum with nonzero demand");
        }
        return;
    }
    ++totals.active_cases;
    if (active_cuts > 16) {
        throw std::runtime_error("more than sixteen negative cuts");
    }
    if (negative_middle > 16 * maximum_cut) {
        throw std::runtime_error("maximal-cut bound failed");
    }
    const std::int64_t margin = sevenfold - 16 * maximum_cut;
    const bool all_cap_at_least_three = level >= 0
        && std::all_of(
            labels.begin(), labels.end(), [level](const int label) {
                return std::min(label, level - label) >= 2;
            }
        );
    if (all_cap_at_least_three && maximum_cut >= 6) {
        ++totals.finite_deep_large_cut_cases;
        if (margin < 0) {
            throw std::runtime_error("finite deep large-cut bound failed");
        }
    }
    if (all_cap_at_least_three && maximum_cut >= 1 && maximum_cut <= 5) {
        std::int64_t best_local = 0;
        std::int64_t worst_local = std::numeric_limits<std::int64_t>::max();
        constexpr std::array<int, 4> low_outputs{0, 2, 4, 6};
        for (const unsigned int mask : maximum_masks) {
            const Labels first = subset(labels, mask, true);
            const Labels second = subset(labels, mask, false);
            std::int64_t local = 0;
            for (const int output : low_outputs) {
                local += multiplicity(first, output, level)
                    * multiplicity(second, output, level);
            }
            best_local = std::max(best_local, local);
            worst_local = std::min(worst_local, local);
        }
        const std::int64_t local_margin =
            best_local - 16 * maximum_cut;
        LocalWitness& witness = totals.deep_small_local[
            static_cast<std::size_t>(maximum_cut)
        ];
        if (!witness.initialized || local_margin < witness.margin) {
            witness = LocalWitness{
                true, local_margin, minus, plus, level
            };
        }
        const std::int64_t worst_local_margin =
            worst_local - 16 * maximum_cut;
        if (worst_local_margin < 0) {
            throw std::runtime_error(
                "finite deep small-cut local bound failed"
            );
        }
        LocalWitness& worst_witness = totals.deep_small_worst_local[
            static_cast<std::size_t>(maximum_cut)
        ];
        if (!worst_witness.initialized
            || worst_local_margin < worst_witness.margin) {
            worst_witness = LocalWitness{
                true, worst_local_margin, minus, plus, level
            };
        }
    }
    if (margin < totals.minimum_margin) {
        totals.minimum_margin = margin;
        totals.minimum_minus = minus;
        totals.minimum_plus = plus;
        totals.minimum_level = level;
    }
    if (margin < 0) {
        ++totals.raw_bound_failures;
        const bool has_simple_current = level >= 0
            && std::find(labels.begin(), labels.end(), level) != labels.end();
        if (has_simple_current) {
            ++totals.raw_failures_with_simple_current;
        } else {
            ++totals.raw_failures_without_simple_current;
            if (margin < totals.minimum_no_k_margin) {
                totals.minimum_no_k_margin = margin;
                totals.minimum_no_k_minus = minus;
                totals.minimum_no_k_plus = plus;
                totals.minimum_no_k_level = level;
            }
        }
        const bool has_cap_at_most_two = level >= 0
            && std::any_of(
                labels.begin(), labels.end(), [level](const int label) {
                    return std::min(label, level - label) <= 1;
                }
            );
        if (!has_cap_at_most_two) {
            ++totals.raw_failures_all_cap_at_least_three;
            if (margin < totals.minimum_deep_margin) {
                totals.minimum_deep_margin = margin;
                totals.minimum_deep_minus = minus;
                totals.minimum_deep_plus = plus;
                totals.minimum_deep_level = level;
            }
        }
        if (!require_raw_bound) {
            return;
        }
        std::cout << "COUNTEREXAMPLE level=" << level << " minus=";
        print_labels(minus);
        std::cout << " plus=";
        print_labels(plus);
        std::cout << " N=" << sevenfold
                  << " T_minus=" << negative_middle
                  << " d=" << maximum_cut
                  << " active=" << active_cuts << '\n';
        throw std::runtime_error("N >= 16d failed");
    }
}

void enumerate_box(
    const Labels& alphabet,
    int level,
    bool require_raw_bound,
    Totals& totals
) {
    const std::vector<Labels> minus_lists = sorted_lists(alphabet, 4);
    const std::vector<Labels> plus_lists = sorted_lists(alphabet, 3);
    for (const Labels& minus : minus_lists) {
        for (const Labels& plus : plus_lists) {
            if (disjoint_support(minus, plus)) {
                check_case(minus, plus, level, require_raw_bound, totals);
            }
        }
    }
}

void print_totals(const char* name, const Totals& totals) {
    std::cout << name << "_cases=" << totals.cases
              << ' ' << name << "_active=" << totals.active_cases;
    if (totals.active_cases != 0U) {
        std::cout << ' ' << name << "_minimum_N_minus_16d="
                  << totals.minimum_margin
                  << ' ' << name << "_minimum_level="
                  << totals.minimum_level
                  << ' ' << name << "_minimum_minus=";
        print_labels(totals.minimum_minus);
        std::cout << ' ' << name << "_minimum_plus=";
        print_labels(totals.minimum_plus);
    }
    std::cout << ' ' << name << "_raw_bound_failures="
              << totals.raw_bound_failures
              << ' ' << name << "_raw_failures_with_k="
              << totals.raw_failures_with_simple_current
              << ' ' << name << "_raw_failures_without_k="
              << totals.raw_failures_without_simple_current
              << ' ' << name << "_raw_failures_all_cap_ge_3="
              << totals.raw_failures_all_cap_at_least_three
              << ' ' << name << "_finite_deep_d_ge_6="
              << totals.finite_deep_large_cut_cases;
    if (totals.raw_failures_without_simple_current != 0U) {
        std::cout << ' ' << name << "_minimum_no_k_margin="
                  << totals.minimum_no_k_margin
                  << ' ' << name << "_minimum_no_k_level="
                  << totals.minimum_no_k_level
                  << ' ' << name << "_minimum_no_k_minus=";
        print_labels(totals.minimum_no_k_minus);
        std::cout << ' ' << name << "_minimum_no_k_plus=";
        print_labels(totals.minimum_no_k_plus);
    }
    if (totals.raw_failures_all_cap_at_least_three != 0U) {
        std::cout << ' ' << name << "_minimum_deep_margin="
                  << totals.minimum_deep_margin
                  << ' ' << name << "_minimum_deep_level="
                  << totals.minimum_deep_level
                  << ' ' << name << "_minimum_deep_minus=";
        print_labels(totals.minimum_deep_minus);
        std::cout << ' ' << name << "_minimum_deep_plus=";
        print_labels(totals.minimum_deep_plus);
    }
    if (totals.cases != 0U) {
        std::cout << ' ' << name << "_minimum_half="
                  << totals.minimum_half
                  << ' ' << name << "_minimum_half_level="
                  << totals.minimum_half_level
                  << ' ' << name << "_minimum_half_minus=";
        print_labels(totals.minimum_half_minus);
        std::cout << ' ' << name << "_minimum_half_plus=";
        print_labels(totals.minimum_half_plus);
    }
    for (std::size_t d = 1U; d <= 5U; ++d) {
        const LocalWitness& witness = totals.deep_small_local[d];
        if (!witness.initialized) {
            continue;
        }
        std::cout << ' ' << name << "_deep_d" << d
                  << "_best_local_margin=" << witness.margin
                  << ' ' << name << "_deep_d" << d
                  << "_level=" << witness.level
                  << ' ' << name << "_deep_d" << d << "_minus=";
        print_labels(witness.minus);
        std::cout << ' ' << name << "_deep_d" << d << "_plus=";
        print_labels(witness.plus);
        const LocalWitness& worst_witness =
            totals.deep_small_worst_local[d];
        std::cout << ' ' << name << "_deep_d" << d
                  << "_worst_local_margin=" << worst_witness.margin;
    }
    std::cout << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 4) {
            throw std::runtime_error(
                "usage: verify_su2_seven_four_minus_even "
                "ordinary_maximum_even_label finite_half_box_maximum_level "
                "[finite_full_maximum_level]"
            );
        }
        const int ordinary_maximum = parse_nonnegative(
            argv[1], "ordinary maximum"
        );
        const int finite_maximum = parse_nonnegative(
            argv[2], "finite maximum"
        );
        const int finite_full_maximum = argc == 4
            ? parse_nonnegative(argv[3], "finite full maximum")
            : 0;
        if ((ordinary_maximum % 2) != 0) {
            throw std::runtime_error("ordinary maximum must be even");
        }

        Labels ordinary_alphabet;
        for (int label = 2; label <= ordinary_maximum; label += 2) {
            ordinary_alphabet.push_back(label);
        }
        Totals ordinary;
        enumerate_box(ordinary_alphabet, -1, true, ordinary);

        Totals finite;
        for (int level = 3; level <= finite_maximum; ++level) {
            const int half_box = (level - 3) / 2;
            Labels finite_alphabet;
            for (int label = 2; label <= half_box; label += 2) {
                finite_alphabet.push_back(label);
            }
            enumerate_box(finite_alphabet, level, true, finite);
        }

        Totals finite_full;
        for (int level = 3; level <= finite_full_maximum; ++level) {
            Labels finite_alphabet;
            for (int label = 2; label <= level; label += 2) {
                finite_alphabet.push_back(label);
            }
            enumerate_box(finite_alphabet, level, false, finite_full);
        }

        std::cout << "SU2_SEVEN_FOUR_MINUS_EVEN "
                  << "ordinary_maximum=" << ordinary_maximum
                  << " finite_maximum_level=" << finite_maximum
                  << " finite_full_maximum_level="
                  << finite_full_maximum << '\n';
        print_totals("ordinary", ordinary);
        print_totals("finite_half_box", finite);
        print_totals("finite_full_discovery", finite_full);
        std::cout << "SU2_SEVEN_FOUR_MINUS_EVEN PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SEVEN_FOUR_MINUS_EVEN FAILURE: "
                  << error.what() << '\n';
        return 1;
    }
}
