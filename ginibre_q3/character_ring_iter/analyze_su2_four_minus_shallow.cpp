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

struct Witness {
    bool initialized = false;
    std::int64_t value = 0;
    int level = 0;
    Labels minus;
    Labels plus;
    std::int64_t invariant = 0;
    std::int64_t pair = 0;
    std::int64_t positive = 0;
    std::int64_t negative = 0;
    std::int64_t maximum_cut = 0;
};

struct Stratum {
    std::uint64_t cases = 0;
    std::uint64_t active_cases = 0;
    Witness raw;
    Witness pair;
    Witness active_pair;
    Witness local_pair;
    Witness local_pair_low;
    Witness exact;
    int maximum_ratio = 0;
};

int parse_positive(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error("maximum level must be positive");
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
    std::vector<std::int64_t> current(
        static_cast<std::size_t>(level + 1), 0
    );
    std::vector<std::int64_t> next(
        static_cast<std::size_t>(level + 1), 0
    );
    current[0] = 1;
    int current_maximum = 0;
    for (const int label : labels) {
        std::fill(next.begin(), next.end(), 0);
        for (int source = 0; source <= current_maximum; ++source) {
            const std::int64_t multiplicity =
                current[static_cast<std::size_t>(source)];
            if (multiplicity == 0) {
                continue;
            }
            const int upper = std::min(
                source + label, 2 * level - source - label
            );
            for (int destination = std::abs(source - label);
                 destination <= upper; destination += 2) {
                next[static_cast<std::size_t>(destination)] += multiplicity;
            }
        }
        current_maximum = std::min(level, current_maximum + label);
        current.swap(next);
    }
    if (output < 0 || output > level) {
        return 0;
    }
    return current[static_cast<std::size_t>(output)];
}

std::int64_t invariant(const Labels& labels, int level) {
    return multiplicity(labels, 0, level);
}

int popcount(unsigned int value) {
    int result = 0;
    while (value != 0U) {
        result += static_cast<int>(value & 1U);
        value >>= 1U;
    }
    return result;
}

Labels subset(const Labels& labels, unsigned int mask, bool selected) {
    Labels result;
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        const bool present = ((mask >> i) & 1U) != 0U;
        if (present == selected) {
            result.push_back(labels[i]);
        }
    }
    return result;
}

void consider(
    Witness& witness,
    std::int64_t value,
    int level,
    const Labels& minus,
    const Labels& plus,
    std::int64_t sevenfold,
    std::int64_t pair,
    std::int64_t positive,
    std::int64_t negative,
    std::int64_t maximum_cut
) {
    if (!witness.initialized || value < witness.value) {
        witness = Witness{
            true, value, level, minus, plus, sevenfold, pair,
            positive, negative, maximum_cut
        };
    }
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

void print_witness(const char* name, const Witness& witness) {
    std::cout << ' ' << name << '=';
    if (!witness.initialized) {
        std::cout << "none";
        return;
    }
    std::cout << witness.value
              << "@k" << witness.level << ":M";
    print_labels(witness.minus);
    std::cout << ":P";
    print_labels(witness.plus);
    std::cout << ":(N,P,T+,T-,d)=("
              << witness.invariant << ',' << witness.pair << ','
              << witness.positive << ',' << witness.negative << ','
              << witness.maximum_cut << ')';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_four_minus_shallow maximum_level"
            );
        }
        const int maximum_level = parse_positive(argv[1]);
        std::array<std::array<Stratum, 5>, 2> minus_shallow{};
        std::array<std::array<Stratum, 4>, 2> plus_shallow{};

        for (int level = 4; level <= maximum_level; ++level) {
            const int shallow = (level % 2) == 0 ? level : level - 1;
            Labels alphabet;
            for (int label = 2; label <= level; label += 2) {
                alphabet.push_back(label);
            }
            const std::vector<Labels> minus_lists =
                sorted_lists(alphabet, 4);
            const std::vector<Labels> plus_lists =
                sorted_lists(alphabet, 3);
            for (const Labels& minus : minus_lists) {
                const int minus_count = static_cast<int>(
                    std::count(minus.begin(), minus.end(), shallow)
                );
                for (const Labels& plus : plus_lists) {
                    if (!disjoint_support(minus, plus)) {
                        continue;
                    }
                    const int plus_count = static_cast<int>(
                        std::count(plus.begin(), plus.end(), shallow)
                    );
                    if (minus_count == 0 && plus_count == 0) {
                        continue;
                    }

                    Labels labels = minus;
                    labels.insert(labels.end(), plus.begin(), plus.end());
                    const std::int64_t sevenfold =
                        invariant(labels, level);
                    std::int64_t pair = 0;
                    std::int64_t pair_low = 0;
                    std::vector<int> probe_outputs;
                    for (int output = 0; output <= 8; output += 2) {
                        probe_outputs.push_back(output);
                        const int reflected = shallow - output;
                        if (reflected >= 0
                            && std::find(
                                probe_outputs.begin(),
                                probe_outputs.end(),
                                reflected
                            ) == probe_outputs.end()) {
                            probe_outputs.push_back(reflected);
                        }
                    }
                    for (std::size_t i = 0U; i < labels.size(); ++i) {
                        for (std::size_t j = i + 1U;
                             j < labels.size(); ++j) {
                            if (labels[i] != labels[j]) {
                                continue;
                            }
                            Labels rest;
                            for (std::size_t t = 0U;
                                 t < labels.size(); ++t) {
                                if (t != i && t != j) {
                                    rest.push_back(labels[t]);
                                }
                            }
                            pair += invariant(rest, level);
                            const Labels rest_pair{
                                rest[0], rest[1]
                            };
                            const Labels rest_triple{
                                rest[2], rest[3], rest[4]
                            };
                            for (const int output : probe_outputs) {
                                pair_low += multiplicity(
                                    rest_pair, output, level
                                ) * multiplicity(
                                    rest_triple, output, level
                                );
                            }
                        }
                    }

                    std::int64_t positive = 0;
                    std::int64_t negative = 0;
                    std::int64_t maximum_cut = 0;
                    std::vector<unsigned int> maximum_masks;
                    constexpr unsigned int full_mask =
                        (1U << 7U) - 1U;
                    for (unsigned int mask = 0U;
                         mask <= full_mask; ++mask) {
                        if (popcount(mask) != 3) {
                            continue;
                        }
                        const std::int64_t term =
                            invariant(subset(labels, mask, true), level)
                            * invariant(
                                subset(labels, mask, false), level
                            );
                        int minus_in_cut = 0;
                        for (std::size_t i = 0U; i < 4U; ++i) {
                            minus_in_cut += static_cast<int>(
                                (mask >> i) & 1U
                            );
                        }
                        if ((minus_in_cut % 2) == 0) {
                            positive += term;
                        } else {
                            negative += term;
                            if (term > maximum_cut) {
                                maximum_cut = term;
                                maximum_masks.clear();
                                maximum_masks.push_back(mask);
                            } else if (
                                term != 0 && term == maximum_cut
                            ) {
                                maximum_masks.push_back(mask);
                            }
                        }
                    }

                    std::int64_t best_local = 0;
                    for (const unsigned int mask : maximum_masks) {
                        const Labels first =
                            subset(labels, mask, true);
                        const Labels second =
                            subset(labels, mask, false);
                        std::int64_t local = 0;
                        for (const int output : probe_outputs) {
                            local += multiplicity(
                                first, output, level
                            ) * multiplicity(
                                second, output, level
                            );
                        }
                        best_local = std::max(best_local, local);
                    }

                    const std::size_t parity =
                        static_cast<std::size_t>(level % 2);
                    Stratum& stratum = minus_count > 0
                        ? minus_shallow[parity][
                            static_cast<std::size_t>(minus_count)]
                        : plus_shallow[parity][
                            static_cast<std::size_t>(plus_count)];
                    ++stratum.cases;
                    consider(
                        stratum.raw, sevenfold - negative, level,
                        minus, plus, sevenfold, pair, positive,
                        negative, maximum_cut
                    );
                    consider(
                        stratum.pair,
                        sevenfold + pair - negative,
                        level, minus, plus, sevenfold, pair, positive,
                        negative, maximum_cut
                    );
                    if (negative > 0) {
                        ++stratum.active_cases;
                        consider(
                            stratum.active_pair,
                            sevenfold + pair - negative,
                            level, minus, plus, sevenfold, pair, positive,
                            negative, maximum_cut
                        );
                        consider(
                            stratum.local_pair,
                            best_local + pair - negative,
                            level, minus, plus, sevenfold, pair, positive,
                            negative, maximum_cut
                        );
                        consider(
                            stratum.local_pair_low,
                            best_local + pair_low - negative,
                            level, minus, plus, sevenfold, pair_low,
                            positive, negative, maximum_cut
                        );
                        const int ratio = static_cast<int>(
                            (negative + maximum_cut - 1) / maximum_cut
                        );
                        stratum.maximum_ratio =
                            std::max(stratum.maximum_ratio, ratio);
                    }
                    consider(
                        stratum.exact,
                        sevenfold + pair + positive - negative,
                        level, minus, plus, sevenfold, pair, positive,
                        negative, maximum_cut
                    );
                }
            }
        }

        std::cout << "SU2_FOUR_MINUS_SHALLOW maximum_level="
                  << maximum_level << '\n';
        for (std::size_t parity = 0U; parity < 2U; ++parity) {
            for (std::size_t count = 1U; count <= 4U; ++count) {
                const Stratum& stratum = minus_shallow[parity][count];
                if (stratum.cases == 0U) {
                    continue;
                }
                std::cout << "parity=" << parity
                          << " sign=minus count=" << count
                          << " cases=" << stratum.cases
                          << " active=" << stratum.active_cases;
                print_witness("raw", stratum.raw);
                print_witness("pair", stratum.pair);
                print_witness("active_pair", stratum.active_pair);
                print_witness("local_pair", stratum.local_pair);
                print_witness("local_pair_low", stratum.local_pair_low);
                print_witness("exact", stratum.exact);
                std::cout << " max_ceil_T_over_d="
                          << stratum.maximum_ratio;
                std::cout << '\n';
            }
            for (std::size_t count = 1U; count <= 3U; ++count) {
                const Stratum& stratum = plus_shallow[parity][count];
                if (stratum.cases == 0U) {
                    continue;
                }
                std::cout << "parity=" << parity
                          << " sign=plus count=" << count
                          << " cases=" << stratum.cases
                          << " active=" << stratum.active_cases;
                print_witness("raw", stratum.raw);
                print_witness("pair", stratum.pair);
                print_witness("active_pair", stratum.active_pair);
                print_witness("local_pair", stratum.local_pair);
                print_witness("local_pair_low", stratum.local_pair_low);
                print_witness("exact", stratum.exact);
                std::cout << " max_ceil_T_over_d="
                          << stratum.maximum_ratio;
                std::cout << '\n';
            }
        }
        std::cout << "SU2_FOUR_MINUS_SHALLOW PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_FOUR_MINUS_SHALLOW FAILURE: "
                  << error.what() << '\n';
        return 1;
    }
}
