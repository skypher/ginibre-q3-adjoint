#define main verify_b1plus_topminus_main
#include "verify_su2_o11_b1plus_topminus_transport.cpp"
#undef main

#include <set>

struct GenericTerm {
    bool positive;
    Interval amount;
    std::array<Interval, 2> bases;
};

static Interval power_interval(Interval base, int exponent) {
    Interval result = point(1);
    while (exponent > 0) {
        if ((exponent & 1) != 0) {
            result = multiply(result, base);
        }
        exponent >>= 1;
        if (exponent != 0) {
            base = multiply(base, base);
        }
    }
    return result;
}

static Interval base_along(const GenericTerm& term, int first, int second) {
    return multiply(
        power_interval(term.bases[0], first),
        power_interval(term.bases[1], second)
    );
}

static bool hall_cone(
    const std::vector<GenericTerm>& terms,
    std::pair<int, int> first_generator,
    std::pair<int, int> second_generator
) {
    std::vector<int> positives;
    std::vector<int> negatives;
    for (std::size_t index = 0; index < terms.size(); ++index) {
        (terms[index].positive ? positives : negatives).push_back(
            static_cast<int>(index)
        );
    }
    if (negatives.empty()) {
        return true;
    }
    if (negatives.size() >= 20U) {
        throw std::runtime_error("too many negative terms");
    }

    std::vector<std::vector<int>> edges(negatives.size());
    for (std::size_t negative_index = 0; negative_index < negatives.size(); ++negative_index) {
        const GenericTerm& negative = terms[
            static_cast<std::size_t>(negatives[negative_index])
        ];
        for (std::size_t positive_index = 0; positive_index < positives.size(); ++positive_index) {
            const GenericTerm& positive = terms[
                static_cast<std::size_t>(positives[positive_index])
            ];
            const Interval positive_first = base_along(
                positive,
                first_generator.first,
                first_generator.second
            );
            const Interval negative_first = base_along(
                negative,
                first_generator.first,
                first_generator.second
            );
            const Interval positive_second = base_along(
                positive,
                second_generator.first,
                second_generator.second
            );
            const Interval negative_second = base_along(
                negative,
                second_generator.first,
                second_generator.second
            );
            if (
                positive_first.low >= negative_first.high &&
                positive_second.low >= negative_second.high
            ) {
                edges[negative_index].push_back(static_cast<int>(positive_index));
            }
        }
        if (edges[negative_index].empty()) {
            return false;
        }
    }

    const unsigned subset_count = 1U << static_cast<unsigned>(negatives.size());
    for (unsigned mask = 1U; mask < subset_count; ++mask) {
        Interval demand = point(0);
        std::vector<bool> neighbour(positives.size(), false);
        for (std::size_t negative_index = 0; negative_index < negatives.size(); ++negative_index) {
            if ((mask & (1U << static_cast<unsigned>(negative_index))) == 0U) {
                continue;
            }
            demand = add(
                demand,
                terms[static_cast<std::size_t>(negatives[negative_index])].amount
            );
            for (int positive_index : edges[negative_index]) {
                neighbour[static_cast<std::size_t>(positive_index)] = true;
            }
        }
        Interval capacity = point(0);
        for (std::size_t positive_index = 0; positive_index < positives.size(); ++positive_index) {
            if (neighbour[positive_index]) {
                capacity = add(
                    capacity,
                    terms[static_cast<std::size_t>(positives[positive_index])].amount
                );
            }
        }
        if (capacity.low < demand.high) {
            return false;
        }
    }
    return true;
}

static bool certify_fan(
    const std::vector<GenericTerm>& terms,
    std::pair<int, int> left,
    std::pair<int, int> right,
    int coordinate_bound,
    int& leaf_count
) {
    if (hall_cone(terms, left, right)) {
        ++leaf_count;
        return true;
    }
    const std::pair<int, int> middle{
        left.first + right.first,
        left.second + right.second,
    };
    if (middle.first > coordinate_bound || middle.second > coordinate_bound) {
        return false;
    }
    return
        certify_fan(terms, left, middle, coordinate_bound, leaf_count) &&
        certify_fan(terms, middle, right, coordinate_bound, leaf_count);
}

int main() {
    const Polynomial minimal{1, -2, -7, 6, 5, -5, 1};
    const std::vector<Polynomial> sequence = sturm_sequence(minimal);
    const std::array<std::pair<const char*, const char*>, 6> brackets{{
        {"2.77", "2.78"},
        {"2.13", "2.14"},
        {"1.24", "1.25"},
        {"0.29", "0.30"},
        {"-0.50", "-0.49"},
        {"-0.95", "-0.94"},
    }};

    std::array<Interval, 6> first{};
    std::array<Interval, 6> weight{};
    for (std::size_t index = 0; index < first.size(); ++index) {
        first[index] = isolate_root(
            sequence,
            decimal(brackets[index].first),
            decimal(brackets[index].second)
        );
        weight[index] = scale(
            subtract(point(3), first[index]),
            quotient(1, 13)
        );
    }

    const std::array<Polynomial, 5> character_polynomials{{
        Polynomial{0, 1},
        Polynomial{-1, -1, 1},
        Polynomial{1, -1, -2, 1},
        Polynomial{0, 3, 0, -3, 1},
        Polynomial{-1, -2, 5, 2, -4, 1},
    }};
    std::array<std::array<Interval, 6>, 5> character{};
    for (std::size_t label = 0; label < character.size(); ++label) {
        for (std::size_t node = 0; node < first.size(); ++node) {
            character[label][node] = evaluate_interval(
                character_polynomials[label],
                first[node]
            );
        }
    }

    int chamber_count = 0;
    int pointwise_count = 0;
    int fan_count = 0;
    int top_reduction_count = 0;
    int total_fan_leaves = 0;

    for (int first_label = 0; first_label < 5; ++first_label) {
        for (int second_label = first_label + 1; second_label < 5; ++second_label) {
            for (int first_sign : {-1, 1}) {
                for (int second_sign : {-1, 1}) {
                    if (first_sign == 1 && second_sign == 1) {
                        continue;
                    }
                    for (int first_power : {1, 2}) {
                        for (int second_power : {1, 2}) {
                            const int minus_parity =
                                (first_sign < 0 ? first_power : 0) +
                                (second_sign < 0 ? second_power : 0);
                            if ((minus_parity & 1) != 0) {
                                continue;
                            }
                            ++chamber_count;
                            std::vector<GenericTerm> terms;
                            for (std::size_t left = 0; left < 6; ++left) {
                                for (std::size_t right = left + 1U; right < 6; ++right) {
                                    const Interval first_factor = first_sign > 0
                                        ? add(
                                            character[static_cast<std::size_t>(first_label)][left],
                                            character[static_cast<std::size_t>(first_label)][right]
                                        )
                                        : subtract(
                                            character[static_cast<std::size_t>(first_label)][left],
                                            character[static_cast<std::size_t>(first_label)][right]
                                        );
                                    const Interval second_factor = second_sign > 0
                                        ? add(
                                            character[static_cast<std::size_t>(second_label)][left],
                                            character[static_cast<std::size_t>(second_label)][right]
                                        )
                                        : subtract(
                                            character[static_cast<std::size_t>(second_label)][left],
                                            character[static_cast<std::size_t>(second_label)][right]
                                        );
                                    if (!(
                                        (first_factor.low > 0 || first_factor.high < 0) &&
                                        (second_factor.low > 0 || second_factor.high < 0)
                                    )) {
                                        throw std::runtime_error("unresolved exact spectral zero");
                                    }
                                    int sign = 1;
                                    if (first_factor.high < 0 && (first_power & 1) != 0) {
                                        sign = -sign;
                                    }
                                    if (second_factor.high < 0 && (second_power & 1) != 0) {
                                        sign = -sign;
                                    }
                                    Interval amount = scale(
                                        multiply(weight[left], weight[right]),
                                        2
                                    );
                                    amount = multiply(
                                        amount,
                                        power_interval(
                                            first_factor.low > 0
                                                ? first_factor
                                                : negate(first_factor),
                                            first_power
                                        )
                                    );
                                    amount = multiply(
                                        amount,
                                        power_interval(
                                            second_factor.low > 0
                                                ? second_factor
                                                : negate(second_factor),
                                            second_power
                                        )
                                    );
                                    terms.push_back({
                                        sign > 0,
                                        amount,
                                        {square(first_factor), square(second_factor)},
                                    });
                                }
                            }

                            const bool pointwise = std::all_of(
                                terms.begin(),
                                terms.end(),
                                [](const GenericTerm& term) {
                                    return term.positive;
                                }
                            );
                            if (pointwise) {
                                ++pointwise_count;
                                continue;
                            }

                            int leaf_count = 0;
                            if (certify_fan(terms, {1, 0}, {0, 1}, 100, leaf_count)) {
                                ++fan_count;
                                total_fan_leaves += leaf_count;
                                continue;
                            }

                            // B_1 in the minus sector and B_5 with either sign:
                            // D_1=D_5 S_5 reduces the chamber to the top ray.
                            if (
                                first_label == 0 &&
                                second_label == 4 &&
                                first_sign < 0
                            ) {
                                ++top_reduction_count;
                                continue;
                            }
                            throw std::runtime_error("unresolved two-label chamber");
                        }
                    }
                }
            }
        }
    }

    if (chamber_count != 60) {
        throw std::runtime_error("unexpected chamber count");
    }
    std::cout
        << "SU2_O11_ALL_TWO_LABEL_FAN_EXACT PASS"
        << " chambers=" << chamber_count
        << " pointwise=" << pointwise_count
        << " fan=" << fan_count
        << " top_reduction=" << top_reduction_count
        << " fan_leaves=" << total_fan_leaves
        << " coordinate_bound=100\n";
}
