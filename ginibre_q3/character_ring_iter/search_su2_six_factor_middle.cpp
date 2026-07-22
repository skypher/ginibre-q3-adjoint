#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <map>
#include <queue>
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
    for (std::size_t index = 0; index < labels.size(); ++index) {
        const bool present = ((mask >> index) & 1U) != 0U;
        if (present == take_mask) {
            result.push_back(labels[index]);
        }
    }
    return result;
}

std::int64_t q3_subset(
    const std::vector<int>& labels,
    std::size_t minus_count
) {
    std::int64_t answer = 0;
    constexpr unsigned int full_mask = (1U << 6U) - 1U;
    for (unsigned int mask = 0; mask <= full_mask; ++mask) {
        const std::int64_t term =
            invariant_multiplicity(selected(labels, mask, true))
            * invariant_multiplicity(selected(labels, mask, false));
        int minus_parity = 0;
        for (std::size_t i = 0; i < minus_count; ++i) {
            minus_parity ^= static_cast<int>((mask >> i) & 1U);
        }
        answer += minus_parity == 0 ? term : -term;
    }
    return answer;
}

using EdgeMatrix = std::vector<std::vector<int>>;

std::vector<EdgeMatrix> admissible_triple_graphs(const std::vector<int>& labels) {
    std::vector<EdgeMatrix> graphs;
    constexpr unsigned int full_mask = (1U << 6U) - 1U;
    for (unsigned int mask = 0; mask <= full_mask; ++mask) {
        if (popcount(mask) != 3 || mask > (full_mask ^ mask)) {
            continue;
        }
        if (invariant_multiplicity(selected(labels, mask, true)) == 0
            || invariant_multiplicity(selected(labels, mask, false)) == 0) {
            continue;
        }
        EdgeMatrix graph(6U, std::vector<int>(6U, 0));
        for (bool take_mask : {false, true}) {
            std::vector<std::size_t> vertices;
            for (std::size_t i = 0; i < labels.size(); ++i) {
                if ((((mask >> i) & 1U) != 0U) == take_mask) {
                    vertices.push_back(i);
                }
            }
            const std::size_t i = vertices[0];
            const std::size_t j = vertices[1];
            const std::size_t k = vertices[2];
            graph[i][j] = (labels[i] + labels[j] - labels[k]) / 2;
            graph[i][k] = (labels[i] + labels[k] - labels[j]) / 2;
            graph[j][k] = (labels[j] + labels[k] - labels[i]) / 2;
        }
        graphs.push_back(std::move(graph));
    }
    return graphs;
}

bool has_separating_vertex_order(const std::vector<int>& labels) {
    const auto graphs = admissible_triple_graphs(labels);
    std::vector<std::size_t> order{0U, 1U, 2U, 3U, 4U, 5U};
    do {
        std::vector<std::size_t> position(6U, 0U);
        for (std::size_t rank = 0; rank < order.size(); ++rank) {
            position[order[rank]] = rank;
        }
        std::set<std::vector<int>> leading_exponents;
        bool separated = true;
        for (const auto& graph : graphs) {
            std::vector<int> exponent(6U, 0);
            for (std::size_t i = 0; i < 6U; ++i) {
                for (std::size_t j = i + 1U; j < 6U; ++j) {
                    if (position[i] < position[j]) {
                        exponent[i] += graph[i][j];
                    } else {
                        exponent[j] += graph[i][j];
                    }
                }
            }
            if (!leading_exponents.insert(exponent).second) {
                separated = false;
                break;
            }
        }
        if (separated) {
            return true;
        }
    } while (std::next_permutation(order.begin(), order.end()));
    return false;
}

std::vector<unsigned int> nontrivial_splits() {
    std::vector<unsigned int> splits;
    constexpr unsigned int full_mask = (1U << 6U) - 1U;
    for (unsigned int mask = 1U; mask < full_mask; ++mask) {
        const unsigned int complement = full_mask ^ mask;
        if ((mask & 1U) == 0U || popcount(mask) < 2
            || popcount(complement) < 2) {
            continue;
        }
        splits.push_back(mask);
    }
    return splits;
}

bool compatible_splits(unsigned int first, unsigned int second) {
    constexpr unsigned int full_mask = (1U << 6U) - 1U;
    const unsigned int first_complement = full_mask ^ first;
    const unsigned int second_complement = full_mask ^ second;
    return (first & second) == 0U
        || (first & second_complement) == 0U
        || (first_complement & second) == 0U
        || (first_complement & second_complement) == 0U;
}

int cut_weight(const EdgeMatrix& graph, unsigned int split) {
    int weight = 0;
    for (std::size_t i = 0; i < 6U; ++i) {
        for (std::size_t j = i + 1U; j < 6U; ++j) {
            const bool first_side = ((split >> i) & 1U) != 0U;
            const bool second_side = ((split >> j) & 1U) != 0U;
            if (first_side != second_side) {
                weight += graph[i][j];
            }
        }
    }
    return weight;
}

bool has_separating_tree(const std::vector<int>& labels) {
    const auto graphs = admissible_triple_graphs(labels);
    const auto splits = nontrivial_splits();
    for (std::size_t i = 0; i < splits.size(); ++i) {
        for (std::size_t j = i + 1U; j < splits.size(); ++j) {
            if (!compatible_splits(splits[i], splits[j])) {
                continue;
            }
            for (std::size_t k = j + 1U; k < splits.size(); ++k) {
                if (!compatible_splits(splits[i], splits[k])
                    || !compatible_splits(splits[j], splits[k])) {
                    continue;
                }
                std::set<std::vector<int>> values;
                bool separated = true;
                for (const auto& graph : graphs) {
                    const std::vector<int> value{
                        cut_weight(graph, splits[i]),
                        cut_weight(graph, splits[j]),
                        cut_weight(graph, splits[k])
                    };
                    if (!values.insert(value).second) {
                        separated = false;
                        break;
                    }
                }
                if (separated) {
                    return true;
                }
            }
        }
    }
    return graphs.size() <= 1U;
}

std::vector<int> flatten_graph(const EdgeMatrix& graph) {
    std::vector<int> flattened;
    for (std::size_t i = 0; i < 6U; ++i) {
        for (std::size_t j = i + 1U; j < 6U; ++j) {
            flattened.push_back(graph[i][j]);
        }
    }
    return flattened;
}

EdgeMatrix uncross_graph(
    EdgeMatrix graph,
    const std::vector<std::size_t>& order,
    bool use_nested_resolution
) {
    std::vector<std::size_t> vertex_at_position(6U, 0U);
    for (std::size_t position = 0; position < 6U; ++position) {
        vertex_at_position[position] = order[position];
    }
    while (true) {
        bool changed = false;
        for (std::size_t a = 0; a < 6U && !changed; ++a) {
            for (std::size_t b = a + 1U; b < 6U && !changed; ++b) {
                for (std::size_t c = b + 1U; c < 6U && !changed; ++c) {
                    for (std::size_t d = c + 1U; d < 6U; ++d) {
                        const std::size_t va = vertex_at_position[a];
                        const std::size_t vb = vertex_at_position[b];
                        const std::size_t vc = vertex_at_position[c];
                        const std::size_t vd = vertex_at_position[d];
                        const auto get_edge = [&graph](std::size_t x, std::size_t y) {
                            return graph[std::min(x, y)][std::max(x, y)];
                        };
                        if (get_edge(va, vc) == 0 || get_edge(vb, vd) == 0) {
                            continue;
                        }
                        --graph[std::min(va, vc)][std::max(va, vc)];
                        --graph[std::min(vb, vd)][std::max(vb, vd)];
                        if (use_nested_resolution) {
                            ++graph[std::min(va, vd)][std::max(va, vd)];
                            ++graph[std::min(vb, vc)][std::max(vb, vc)];
                        } else {
                            ++graph[std::min(va, vb)][std::max(va, vb)];
                            ++graph[std::min(vc, vd)][std::max(vc, vd)];
                        }
                        changed = true;
                        break;
                    }
                }
            }
        }
        if (!changed) {
            return graph;
        }
    }
}

bool has_injective_uncrossing(const std::vector<int>& labels) {
    const auto graphs = admissible_triple_graphs(labels);
    std::vector<std::size_t> order{0U, 1U, 2U, 3U, 4U, 5U};
    do {
        for (bool nested : {false, true}) {
            std::set<std::vector<int>> outputs;
            bool injective = true;
            for (const auto& graph : graphs) {
                if (!outputs.insert(
                        flatten_graph(uncross_graph(graph, order, nested))
                    ).second) {
                    injective = false;
                    break;
                }
            }
            if (injective) {
                return true;
            }
        }
    } while (std::next_permutation(order.begin(), order.end()));
    return graphs.size() <= 1U;
}

std::pair<bool, std::array<std::size_t, 4>> first_crossing(
    const EdgeMatrix& graph,
    const std::vector<std::size_t>& order
) {
    for (std::size_t a = 0; a < 6U; ++a) {
        for (std::size_t b = a + 1U; b < 6U; ++b) {
            for (std::size_t c = b + 1U; c < 6U; ++c) {
                for (std::size_t d = c + 1U; d < 6U; ++d) {
                    const std::size_t va = order[a];
                    const std::size_t vb = order[b];
                    const std::size_t vc = order[c];
                    const std::size_t vd = order[d];
                    if (graph[std::min(va, vc)][std::max(va, vc)] != 0
                        && graph[std::min(vb, vd)][std::max(vb, vd)] != 0) {
                        return {true, {va, vb, vc, vd}};
                    }
                }
            }
        }
    }
    return {false, {0U, 0U, 0U, 0U}};
}

std::set<std::vector<int>> all_uncrossings(
    const EdgeMatrix& initial,
    const std::vector<std::size_t>& order
) {
    std::set<std::vector<int>> seen;
    std::set<std::vector<int>> outputs;
    std::queue<EdgeMatrix> pending;
    seen.insert(flatten_graph(initial));
    pending.push(initial);
    while (!pending.empty()) {
        EdgeMatrix graph = std::move(pending.front());
        pending.pop();
        const auto [found, crossing] = first_crossing(graph, order);
        if (!found) {
            outputs.insert(flatten_graph(graph));
            continue;
        }
        const auto [a, b, c, d] = crossing;
        for (bool nested : {false, true}) {
            EdgeMatrix next = graph;
            --next[std::min(a, c)][std::max(a, c)];
            --next[std::min(b, d)][std::max(b, d)];
            if (nested) {
                ++next[std::min(a, d)][std::max(a, d)];
                ++next[std::min(b, c)][std::max(b, c)];
            } else {
                ++next[std::min(a, b)][std::max(a, b)];
                ++next[std::min(c, d)][std::max(c, d)];
            }
            const auto key = flatten_graph(next);
            if (seen.insert(key).second) {
                pending.push(std::move(next));
            }
        }
    }
    return outputs;
}

bool augment_matching(
    std::size_t source,
    const std::vector<std::vector<std::size_t>>& adjacency,
    std::vector<int>& output_match,
    std::vector<bool>& visited
) {
    for (std::size_t output : adjacency[source]) {
        if (visited[output]) {
            continue;
        }
        visited[output] = true;
        if (output_match[output] < 0
            || augment_matching(
                static_cast<std::size_t>(output_match[output]),
                adjacency, output_match, visited
            )) {
            output_match[output] = static_cast<int>(source);
            return true;
        }
    }
    return false;
}

bool support_has_full_matching(
    const std::vector<EdgeMatrix>& graphs,
    const std::vector<std::size_t>& order
) {
    std::map<std::vector<int>, std::size_t> output_index;
    std::vector<std::vector<std::vector<int>>> supports;
    for (const auto& graph : graphs) {
        const auto output_set = all_uncrossings(graph, order);
        supports.emplace_back(output_set.begin(), output_set.end());
        for (const auto& output : output_set) {
            if (output_index.find(output) == output_index.end()) {
                output_index.emplace(output, output_index.size());
            }
        }
    }
    std::vector<std::vector<std::size_t>> adjacency(graphs.size());
    for (std::size_t source = 0; source < supports.size(); ++source) {
        for (const auto& output : supports[source]) {
            adjacency[source].push_back(output_index.at(output));
        }
    }
    std::vector<int> output_match(output_index.size(), -1);
    for (std::size_t source = 0; source < graphs.size(); ++source) {
        std::vector<bool> visited(output_index.size(), false);
        if (!augment_matching(source, adjacency, output_match, visited)) {
            return false;
        }
    }
    return true;
}

bool has_uncrossing_support_matching(const std::vector<int>& labels) {
    const auto graphs = admissible_triple_graphs(labels);
    std::vector<std::size_t> order{0U, 1U, 2U, 3U, 4U, 5U};
    return support_has_full_matching(graphs, order);
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
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << ']';
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
                "usage: search_su2_six_factor_middle maximum_label"
            );
        }
        const int maximum_label = parse_positive(argv[1]);
        Witness core_witness;
        Witness exact_witness;
        std::uint64_t cases = 0;
        std::uint64_t separation_tests = 0;
        std::vector<int> separation_failure;
        std::vector<int> tree_separation_failure;
        std::vector<int> uncrossing_failure;
        std::vector<int> support_matching_failure;
        std::array<std::int64_t, 10> maximum_triples_by_small_dimension{};
        std::array<std::vector<int>, 10> small_dimension_witness{};

        for (int minus_count = 2; minus_count <= 6; minus_count += 2) {
            const auto minus_lists = combinations(maximum_label, minus_count);
            const auto plus_lists = combinations(maximum_label, 6 - minus_count);
            for (const auto& minus : minus_lists) {
                for (const auto& plus : plus_lists) {
                    if (!supports_are_disjoint(minus, plus)) {
                        continue;
                    }
                    ++cases;
                    std::vector<int> labels = minus;
                    labels.insert(labels.end(), plus.begin(), plus.end());
                    const std::int64_t invariant = invariant_multiplicity(labels);

                    if (minus_count == 6) {
                        ++separation_tests;
                        if (separation_failure.empty()
                            && !has_separating_vertex_order(labels)) {
                            separation_failure = labels;
                        }
                        if (tree_separation_failure.empty()
                            && !has_separating_tree(labels)) {
                            tree_separation_failure = labels;
                        }
                        if (uncrossing_failure.empty()
                            && !has_injective_uncrossing(labels)) {
                            uncrossing_failure = labels;
                        }
                        if (support_matching_failure.empty()
                            && !has_uncrossing_support_matching(labels)) {
                            support_matching_failure = labels;
                        }
                        const auto graphs = admissible_triple_graphs(labels);
                        if (invariant >= 0 && invariant < 10
                            && static_cast<std::int64_t>(graphs.size())
                                > maximum_triples_by_small_dimension[
                                    static_cast<std::size_t>(invariant)]) {
                            maximum_triples_by_small_dimension[
                                static_cast<std::size_t>(invariant)] =
                                static_cast<std::int64_t>(graphs.size());
                            small_dimension_witness[
                                static_cast<std::size_t>(invariant)] = labels;
                        }
                    }

                    std::int64_t pair_term = 0;
                    for (std::size_t i = 0; i < labels.size(); ++i) {
                        for (std::size_t j = i + 1; j < labels.size(); ++j) {
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
                    constexpr unsigned int full_mask = (1U << 6U) - 1U;
                    for (unsigned int mask = 0; mask <= full_mask; ++mask) {
                        if (popcount(mask) != 3) {
                            continue;
                        }
                        const unsigned int complement = full_mask ^ mask;
                        if (mask > complement) {
                            continue;
                        }
                        const std::int64_t term =
                            invariant_multiplicity(selected(labels, mask, true))
                            * invariant_multiplicity(selected(labels, mask, false));
                        int minus_parity = 0;
                        for (std::size_t i = 0; i < minus.size(); ++i) {
                            minus_parity ^= static_cast<int>((mask >> i) & 1U);
                        }
                        if (minus_parity == 0) {
                            positive_middle += term;
                        } else {
                            negative_middle += term;
                        }
                    }

                    const std::int64_t core = invariant - negative_middle;
                    const std::int64_t exact_half = invariant + pair_term
                        + positive_middle - negative_middle;
                    if (2 * exact_half != q3_subset(labels, minus.size())) {
                        throw std::runtime_error(
                            "six-factor middle decomposition disagrees with subset sum"
                        );
                    }
                    if (minus_count == 6) {
                        const std::int64_t triple_splits =
                            static_cast<std::int64_t>(
                                admissible_triple_graphs(labels).size()
                            );
                        if (invariant < triple_splits) {
                            throw std::runtime_error(
                                "sixfold invariant does not dominate triple splits"
                            );
                        }
                    }
                    consider(
                        core_witness, core, minus, plus, invariant,
                        negative_middle, positive_middle, pair_term
                    );
                    consider(
                        exact_witness, exact_half, minus, plus, invariant,
                        negative_middle, positive_middle, pair_term
                    );
                }
            }
        }

        std::cout << "SU2_SIX_FACTOR_MIDDLE cases=" << cases
                  << " maximum_label=" << maximum_label << '\n';
        std::cout << "leading_order_tests=" << separation_tests
                  << " result="
                  << (separation_failure.empty() ? "PASS" : "FAIL");
        if (!separation_failure.empty()) {
            std::cout << " labels=";
            print_vector(separation_failure);
        }
        std::cout << '\n';
        std::cout << "small_dimension_table";
        for (std::size_t dimension = 0; dimension < 10U; ++dimension) {
            std::cout << " N" << dimension << "="
                      << maximum_triples_by_small_dimension[dimension];
            if (!small_dimension_witness[dimension].empty()) {
                std::cout << '@';
                print_vector(small_dimension_witness[dimension]);
            }
        }
        std::cout << '\n';
        std::cout << "uncrossing_support_matching_tests=" << separation_tests
                  << " result="
                  << (support_matching_failure.empty() ? "PASS" : "FAIL");
        if (!support_matching_failure.empty()) {
            std::cout << " labels=";
            print_vector(support_matching_failure);
        }
        std::cout << '\n';
        std::cout << "uncrossing_tests=" << separation_tests
                  << " result="
                  << (uncrossing_failure.empty() ? "PASS" : "FAIL");
        if (!uncrossing_failure.empty()) {
            std::cout << " labels=";
            print_vector(uncrossing_failure);
        }
        std::cout << '\n';
        std::cout << "tree_valuation_tests=" << separation_tests
                  << " result="
                  << (tree_separation_failure.empty() ? "PASS" : "FAIL");
        if (!tree_separation_failure.empty()) {
            std::cout << " labels=";
            print_vector(tree_separation_failure);
        }
        std::cout << '\n';
        std::cout << "core_minimum=" << core_witness.value << " minus=";
        print_vector(core_witness.minus);
        std::cout << " plus=";
        print_vector(core_witness.plus);
        std::cout << " invariant=" << core_witness.invariant
                  << " negative_middle=" << core_witness.negative_middle
                  << " positive_middle=" << core_witness.positive_middle
                  << " pair_term=" << core_witness.pair_term << '\n';
        std::cout << "exact_half_minimum=" << exact_witness.value << " minus=";
        print_vector(exact_witness.minus);
        std::cout << " plus=";
        print_vector(exact_witness.plus);
        std::cout << " invariant=" << exact_witness.invariant
                  << " negative_middle=" << exact_witness.negative_middle
                  << " positive_middle=" << exact_witness.positive_middle
                  << " pair_term=" << exact_witness.pair_term << '\n';
        std::cout << "SU2_SIX_FACTOR_MIDDLE "
                  << (exact_witness.value < 0 ? "FAIL" : "PASS") << '\n';
        return exact_witness.value < 0 ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SIX_FACTOR_MIDDLE FAILURE: " << error.what() << '\n';
        return 1;
    }
}
