#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Graph = std::vector<std::vector<int>>;
using Key = std::vector<int>;

int popcount(unsigned int value) {
    int answer = 0;
    while (value != 0U) {
        answer += static_cast<int>(value & 1U);
        value >>= 1U;
    }
    return answer;
}

bool edges_cross(
    std::size_t first_left,
    std::size_t first_right,
    std::size_t second_left,
    std::size_t second_right
) {
    return (first_left < second_left && second_left < first_right
            && first_right < second_right)
        || (second_left < first_left && first_left < second_right
            && second_right < first_right);
}

bool can_add_edge(const Graph& graph, std::size_t left, std::size_t right) {
    for (std::size_t first = 0U; first < graph.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < graph.size(); ++second) {
            if (graph[first][second] != 0
                && edges_cross(left, right, first, second)) {
                return false;
            }
        }
    }
    return true;
}

void enumerate_graphs_rec(
    std::vector<int>& remaining,
    Graph& graph,
    std::vector<Graph>& output
);

void distribute_vertex_degree(
    std::size_t vertex,
    std::size_t neighbor,
    int needed,
    std::vector<int>& remaining,
    Graph& graph,
    std::vector<Graph>& output
) {
    if (neighbor == graph.size()) {
        if (needed == 0) {
            const int saved = remaining[vertex];
            remaining[vertex] = 0;
            enumerate_graphs_rec(remaining, graph, output);
            remaining[vertex] = saved;
        }
        return;
    }
    int available = 0;
    for (std::size_t next = neighbor; next < graph.size(); ++next) {
        available += remaining[next];
    }
    if (available < needed) {
        return;
    }
    const int maximum = std::min(needed, remaining[neighbor]);
    for (int count = 0; count <= maximum; ++count) {
        if (count != 0 && !can_add_edge(graph, vertex, neighbor)) {
            continue;
        }
        graph[vertex][neighbor] = count;
        remaining[neighbor] -= count;
        distribute_vertex_degree(
            vertex, neighbor + 1U, needed - count,
            remaining, graph, output
        );
        remaining[neighbor] += count;
        graph[vertex][neighbor] = 0;
    }
}

void enumerate_graphs_rec(
    std::vector<int>& remaining,
    Graph& graph,
    std::vector<Graph>& output
) {
    std::size_t vertex = 0U;
    while (vertex < graph.size() && remaining[vertex] == 0) {
        ++vertex;
    }
    if (vertex == graph.size()) {
        output.push_back(graph);
        return;
    }
    distribute_vertex_degree(
        vertex, vertex + 1U, remaining[vertex],
        remaining, graph, output
    );
}

std::vector<Graph> noncrossing_graphs(const std::vector<int>& degrees) {
    const int total = [&degrees]() {
        int sum = 0;
        for (int degree : degrees) {
            sum += degree;
        }
        return sum;
    }();
    if ((total & 1) != 0) {
        return {};
    }
    std::vector<int> remaining = degrees;
    Graph graph(degrees.size(), std::vector<int>(degrees.size(), 0));
    std::vector<Graph> output;
    enumerate_graphs_rec(remaining, graph, output);
    return output;
}

Graph graph_sum(const Graph& first, const Graph& second) {
    Graph result = first;
    for (std::size_t i = 0U; i < result.size(); ++i) {
        for (std::size_t j = i + 1U; j < result.size(); ++j) {
            result[i][j] += second[i][j];
        }
    }
    return result;
}

bool has_crossing(const Graph& graph) {
    for (std::size_t a = 0U; a < graph.size(); ++a) {
        for (std::size_t b = a + 1U; b < graph.size(); ++b) {
            for (std::size_t c = b + 1U; c < graph.size(); ++c) {
                for (std::size_t d = c + 1U; d < graph.size(); ++d) {
                    if (graph[a][c] != 0 && graph[b][d] != 0) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

Key flatten(const Graph& graph) {
    Key key;
    key.reserve(graph.size() * (graph.size() - 1U) / 2U);
    for (std::size_t i = 0U; i < graph.size(); ++i) {
        for (std::size_t j = i + 1U; j < graph.size(); ++j) {
            key.push_back(graph[i][j]);
        }
    }
    return key;
}

struct Counts {
    std::uint64_t odd_sources = 0U;
    std::uint64_t odd_zero_sources = 0U;
    std::uint64_t positive_sources = 0U;
    std::uint64_t positive_zero_sources = 0U;
    std::size_t odd_rank = 0U;
    bool fiber_parity = true;
    std::uint64_t fiber_odd = 0U;
    std::uint64_t fiber_even = 0U;
    bool xor_payment = true;
    unsigned int xor_mask = 0U;
    std::uint64_t xor_demand = 0U;
    std::uint64_t xor_capacity = 0U;
    std::uint64_t order_count = 0U;
    std::uint64_t rank_sum = 0U;
    std::uint64_t odd_pair_collisions = 0U;
    std::uint64_t pair_collision_sum = 0U;
    std::uint64_t component_count = 0U;
    std::uint64_t minimum_component_states = 0U;
    std::uint64_t minimum_component_relations = 0U;
    std::uint64_t minimum_component_positive = 0U;
    std::uint64_t minimum_component_slack = 0U;
    std::uint64_t active_component_count = 0U;
    std::uint64_t loaded_component_count = 0U;
    Key alpha_witness;
};

struct ToricFiber {
    std::array<std::uint64_t, 2U> parity_counts{};
    std::vector<unsigned int> odd_cuts;
    std::vector<unsigned int> nonempty_even_cuts;
};

struct CutFiberData {
    std::vector<unsigned int> odd_cuts;
    std::vector<unsigned int> positive_cuts;
};

struct AlphaAverage {
    std::uint64_t odd_occurrences = 0U;
    std::uint64_t positive_occurrences = 0U;
    std::uint64_t hit_orders = 0U;
};

using WeightedPairCounts = std::map<
    std::pair<unsigned int, std::uint64_t>, std::uint64_t
>;

struct FiberCounts {
    std::uint64_t odd = 0U;
    std::uint64_t positive = 0U;
};

using GraphCache = std::map<std::vector<int>, std::vector<Graph>>;
using ToricCache = std::map<std::vector<int>, std::vector<Key>>;

const std::vector<Graph>& cached_graphs(
    const std::vector<int>& degrees,
    GraphCache& cache
) {
    const auto found = cache.find(degrees);
    if (found != cache.end()) {
        return found->second;
    }
    return cache.emplace(degrees, noncrossing_graphs(degrees)).first->second;
}

Counts analyze_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    GraphCache& cache
) {
    const std::size_t factors = labels.size();
    const unsigned int subset_limit = 1U << factors;
    std::set<Key> odd_nonzero_outputs;
    Counts counts;
    for (unsigned int subset = 0U; subset < subset_limit; ++subset) {
        std::vector<int> left(factors + 1U, 0);
        std::vector<int> right(factors + 1U, 0);
        for (std::size_t i = 0U; i < factors; ++i) {
            if (((subset >> i) & 1U) != 0U) {
                left[i] = labels[i];
            } else {
                right[i] = labels[i];
            }
        }
        right[factors] = target;
        const std::vector<Graph>& left_graphs = cached_graphs(left, cache);
        const std::vector<Graph>& right_graphs = cached_graphs(right, cache);
        const bool odd = (popcount(subset & minus_mask) & 1) != 0;
        for (const Graph& first : left_graphs) {
            for (const Graph& second : right_graphs) {
                const Graph product = graph_sum(first, second);
                const bool zero = has_crossing(product);
                if (odd) {
                    ++counts.odd_sources;
                    if (zero) {
                        ++counts.odd_zero_sources;
                    } else {
                        odd_nonzero_outputs.insert(flatten(product));
                    }
                } else if (subset != 0U) {
                    ++counts.positive_sources;
                    if (zero) {
                        ++counts.positive_zero_sources;
                    }
                }
            }
        }
    }
    counts.odd_rank = odd_nonzero_outputs.size();
    return counts;
}

bool valid_toric_top(
    const std::vector<int>& degrees,
    const Key& top
) {
    int top_total = 0;
    int degree_total = 0;
    for (std::size_t i = 0U; i < degrees.size(); ++i) {
        top_total += top[i];
        degree_total += degrees[i];
    }
    if (2 * top_total != degree_total) {
        return false;
    }
    int top_before = 0;
    int bottom_through = 0;
    for (std::size_t i = 0U; i < degrees.size(); ++i) {
        bottom_through += degrees[i] - top[i];
        if (bottom_through > top_before) {
            return false;
        }
        top_before += top[i];
    }
    return true;
}

void enumerate_toric_tops_rec(
    const std::vector<int>& degrees,
    std::size_t vertex,
    Key& top,
    std::vector<Key>& output
) {
    if (vertex == degrees.size()) {
        if (valid_toric_top(degrees, top)) {
            output.push_back(top);
        }
        return;
    }
    for (int count = 0; count <= degrees[vertex]; ++count) {
        top[vertex] = count;
        enumerate_toric_tops_rec(degrees, vertex + 1U, top, output);
    }
}

std::vector<Key> toric_tops(const std::vector<int>& degrees) {
    Key top(degrees.size(), 0);
    std::vector<Key> output;
    enumerate_toric_tops_rec(degrees, 0U, top, output);
    return output;
}

const std::vector<Key>& cached_toric_tops(
    const std::vector<int>& degrees,
    ToricCache& cache
) {
    const auto found = cache.find(degrees);
    if (found != cache.end()) {
        return found->second;
    }
    return cache.emplace(degrees, toric_tops(degrees)).first->second;
}

Counts analyze_gt_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache,
    const std::vector<std::size_t>& order,
    std::map<Key, AlphaAverage>* alpha_average = nullptr,
    std::map<Key, FiberCounts>* state_fibers = nullptr,
    WeightedPairCounts* weighted_pair_counts = nullptr,
    std::map<unsigned int, std::uint64_t>* channel_capacity = nullptr,
    std::vector<std::vector<unsigned int>>* odd_fibers = nullptr,
    std::vector<CutFiberData>* cut_fibers = nullptr
) {
    const std::size_t factors = labels.size();
    if (order.size() != factors + 1U) {
        throw std::runtime_error("GT vertex order has the wrong size");
    }
    std::vector<std::size_t> position(factors + 1U, 0U);
    std::vector<bool> seen(factors + 1U, false);
    for (std::size_t i = 0U; i < order.size(); ++i) {
        if (order[i] >= order.size()) {
            throw std::runtime_error("GT vertex order is out of range");
        }
        if (seen[order[i]]) {
            throw std::runtime_error("GT vertex order has a duplicate");
        }
        seen[order[i]] = true;
        position[order[i]] = i;
    }
    const unsigned int subset_limit = 1U << factors;
    std::set<Key> odd_outputs;
    std::map<Key, ToricFiber> fibers;
    std::map<unsigned int, std::uint64_t> positive_capacity;
    Counts counts;
    for (unsigned int subset = 0U; subset < subset_limit; ++subset) {
        std::vector<int> left(factors + 1U, 0);
        std::vector<int> right(factors + 1U, 0);
        for (std::size_t i = 0U; i < factors; ++i) {
            if (((subset >> i) & 1U) != 0U) {
                left[position[i]] = labels[i];
            } else {
                right[position[i]] = labels[i];
            }
        }
        right[position[factors]] = target;
        const std::vector<Key>& left_tops = cached_toric_tops(left, cache);
        const std::vector<Key>& right_tops = cached_toric_tops(right, cache);
        const bool odd = (popcount(subset & minus_mask) & 1) != 0;
        for (const Key& first : left_tops) {
            for (const Key& second : right_tops) {
                Key product(first.size(), 0);
                for (std::size_t i = 0U; i < product.size(); ++i) {
                    product[i] = first[i] + second[i];
                }
                if (odd) {
                    ++counts.odd_sources;
                    odd_outputs.insert(product);
                    ToricFiber& fiber = fibers[std::move(product)];
                    ++fiber.parity_counts[1U];
                    fiber.odd_cuts.push_back(subset);
                } else if (subset != 0U) {
                    ++counts.positive_sources;
                    ToricFiber& fiber = fibers[std::move(product)];
                    ++fiber.parity_counts[0U];
                    fiber.nonempty_even_cuts.push_back(subset);
                    ++positive_capacity[subset];
                    if (channel_capacity != nullptr) {
                        ++(*channel_capacity)[subset];
                    }
                } else {
                    ++fibers[std::move(product)].parity_counts[0U];
                }
            }
        }
    }
    counts.odd_rank = odd_outputs.size();
    std::map<unsigned int, std::uint64_t> used_capacity;
    std::map<unsigned int, std::uint64_t> xor_demand;
    for (const auto& [key, fiber] : fibers) {
        if (alpha_average != nullptr || state_fibers != nullptr) {
            Key raw_top(key.size(), 0);
            for (std::size_t vertex = 0U;
                 vertex < position.size(); ++vertex) {
                raw_top[vertex] = key[position[vertex]];
            }
            if (alpha_average != nullptr) {
                AlphaAverage& average = (*alpha_average)[raw_top];
                average.odd_occurrences += static_cast<std::uint64_t>(
                    fiber.odd_cuts.size()
                );
                average.positive_occurrences += static_cast<std::uint64_t>(
                    fiber.nonempty_even_cuts.size()
                );
                if (!fiber.odd_cuts.empty()) {
                    ++average.hit_orders;
                }
            }
            if (state_fibers != nullptr) {
                state_fibers->emplace(
                    std::move(raw_top),
                    FiberCounts{
                        static_cast<std::uint64_t>(fiber.odd_cuts.size()),
                        static_cast<std::uint64_t>(
                            fiber.nonempty_even_cuts.size()
                        )
                    }
                );
            }
        }
        const std::uint64_t odd_in_fiber
            = static_cast<std::uint64_t>(fiber.odd_cuts.size());
        if (cut_fibers != nullptr
            && (!fiber.odd_cuts.empty()
                || !fiber.nonempty_even_cuts.empty())) {
            cut_fibers->push_back(CutFiberData{
                fiber.odd_cuts, fiber.nonempty_even_cuts
            });
        }
        if (odd_fibers != nullptr && odd_in_fiber >= 2U) {
            odd_fibers->push_back(fiber.odd_cuts);
        }
        if (odd_in_fiber >= 2U) {
            counts.odd_pair_collisions
                += odd_in_fiber * (odd_in_fiber - 1U) / 2U;
            if (weighted_pair_counts != nullptr) {
                for (std::size_t first = 0U;
                     first < fiber.odd_cuts.size(); ++first) {
                    for (std::size_t second = first + 1U;
                         second < fiber.odd_cuts.size(); ++second) {
                        ++(*weighted_pair_counts)[{
                            fiber.odd_cuts[first] ^ fiber.odd_cuts[second],
                            odd_in_fiber
                        }];
                    }
                }
            }
        }
        if (fiber.parity_counts[1U] > fiber.parity_counts[0U]) {
            counts.fiber_parity = false;
            counts.fiber_even = fiber.parity_counts[0U];
            counts.fiber_odd = fiber.parity_counts[1U];
        }
        if (!fiber.odd_cuts.empty()) {
            const unsigned int anchor = fiber.odd_cuts.front();
            const std::size_t local_payments = std::min(
                fiber.odd_cuts.size() - 1U,
                fiber.nonempty_even_cuts.size()
            );
            for (std::size_t i = 0U; i < local_payments; ++i) {
                ++used_capacity[fiber.nonempty_even_cuts[i]];
            }
            for (std::size_t i = 1U + local_payments;
                 i < fiber.odd_cuts.size(); ++i) {
                ++xor_demand[anchor ^ fiber.odd_cuts[i]];
            }
        }
    }
    for (const auto& [mask, demand] : xor_demand) {
        const std::uint64_t total_capacity = positive_capacity[mask];
        const std::uint64_t local_use = used_capacity[mask];
        if (local_use > total_capacity) {
            throw std::runtime_error("local GT use exceeds channel capacity");
        }
        const std::uint64_t capacity = total_capacity - local_use;
        if (demand > capacity) {
            counts.xor_payment = false;
            counts.xor_mask = mask;
            counts.xor_demand = demand;
            counts.xor_capacity = capacity;
            break;
        }
    }
    return counts;
}

Counts analyze_gt_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::vector<std::size_t> order(labels.size() + 1U, 0U);
    for (std::size_t i = 0U; i < order.size(); ++i) {
        order[i] = i;
    }
    return analyze_gt_case(labels, target, minus_mask, cache, order);
}

Counts analyze_gt_zigzag_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::vector<std::size_t> order;
    order.reserve(labels.size() + 1U);
    for (std::size_t parity = 0U; parity < 2U; ++parity) {
        for (std::size_t i = parity; i < labels.size() + 1U; i += 2U) {
            order.push_back(i);
        }
    }
    return analyze_gt_case(labels, target, minus_mask, cache, order);
}

Counts analyze_gt_label_rounds_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::map<int, std::vector<std::size_t>> label_groups;
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        label_groups[labels[i]].push_back(i);
    }
    const bool target_label_present = label_groups.contains(target);
    if (target_label_present) {
        label_groups[target].push_back(labels.size());
    }
    std::size_t maximum_group = 0U;
    for (const auto& [label, group] : label_groups) {
        (void)label;
        maximum_group = std::max(maximum_group, group.size());
    }
    std::vector<std::size_t> order;
    order.reserve(labels.size() + 1U);
    for (std::size_t round = 0U; round < maximum_group; ++round) {
        for (const auto& [label, group] : label_groups) {
            (void)label;
            if (round < group.size()) {
                order.push_back(group[round]);
            }
        }
    }
    if (!target_label_present) {
        order.push_back(labels.size());
    }
    return analyze_gt_case(labels, target, minus_mask, cache, order);
}

Counts analyze_gt_xor_rounds_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::map<int, std::vector<std::size_t>> label_groups;
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        label_groups[labels[i]].push_back(i);
    }
    std::size_t maximum_group = 0U;
    for (const auto& [label, group] : label_groups) {
        (void)label;
        maximum_group = std::max(maximum_group, group.size());
    }
    std::vector<std::size_t> order;
    order.reserve(labels.size() + 1U);
    for (std::size_t round = 0U; round < maximum_group; ++round) {
        for (const auto& [label, group] : label_groups) {
            (void)label;
            if (round < group.size()) {
                order.push_back(group[round]);
            }
        }
    }
    order.push_back(labels.size());
    return analyze_gt_case(labels, target, minus_mask, cache, order);
}

Counts analyze_gt_deficit_xor_label_rounds_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    Counts counts = analyze_gt_label_rounds_case(
        labels, target, minus_mask, cache
    );
    if (counts.odd_sources <= counts.positive_sources) {
        counts.xor_payment = true;
    }
    return counts;
}

Counts analyze_gt_xor_insert_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::map<int, std::vector<std::size_t>> label_groups;
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        label_groups[labels[i]].push_back(i);
    }
    std::size_t maximum_group = 0U;
    for (const auto& [label, group] : label_groups) {
        (void)label;
        maximum_group = std::max(maximum_group, group.size());
    }
    std::vector<std::size_t> remainder_order;
    remainder_order.reserve(labels.size());
    for (std::size_t round = 0U; round < maximum_group; ++round) {
        for (const auto& [label, group] : label_groups) {
            (void)label;
            if (round < group.size()) {
                remainder_order.push_back(group[round]);
            }
        }
    }
    Counts last;
    for (std::size_t gap = 0U; gap <= remainder_order.size(); ++gap) {
        std::vector<std::size_t> order = remainder_order;
        order.insert(
            order.begin() + static_cast<std::ptrdiff_t>(gap), labels.size()
        );
        last = analyze_gt_case(labels, target, minus_mask, cache, order);
        if (gap == 0U && last.odd_sources <= last.positive_sources) {
            last.xor_payment = true;
            return last;
        }
        if (last.xor_payment) {
            return last;
        }
    }
    return last;
}

Counts analyze_gt_scalar_balanced_case(
    const std::vector<int>& labels,
    unsigned int minus_mask,
    ToricCache& cache
) {
    const std::size_t factors = labels.size();
    std::map<int, std::vector<std::size_t>> label_groups;
    for (std::size_t i = 0U; i < factors; ++i) {
        label_groups[labels[i]].push_back(i);
    }
    std::map<int, std::size_t> used;
    std::vector<std::size_t> order;
    order.reserve(factors + 1U);
    for (std::size_t step = 0U; step < factors; ++step) {
        bool selected = false;
        int selected_label = 0;
        std::int64_t selected_deficit = 0;
        for (const auto& [label, group] : label_groups) {
            const std::size_t already_used = used[label];
            if (already_used == group.size()) {
                continue;
            }
            const std::int64_t deficit
                = static_cast<std::int64_t>(step + 1U)
                    * static_cast<std::int64_t>(group.size())
                - static_cast<std::int64_t>(already_used)
                    * static_cast<std::int64_t>(factors);
            if (!selected || deficit > selected_deficit
                || (deficit == selected_deficit
                    && group.size()
                        < label_groups.at(selected_label).size())) {
                selected = true;
                selected_label = label;
                selected_deficit = deficit;
            }
        }
        if (!selected) {
            throw std::runtime_error("balanced label schedule is incomplete");
        }
        const std::size_t selected_index = used[selected_label];
        order.push_back(label_groups.at(selected_label)[selected_index]);
        ++used[selected_label];
    }
    order.push_back(factors);
    return analyze_gt_case(labels, 0, minus_mask, cache, order);
}

Counts analyze_gt_scalar_best_case(
    const std::vector<int>& labels,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::vector<int> degree_order = labels;
    std::sort(degree_order.begin(), degree_order.end());
    Counts last;
    do {
        std::map<int, std::size_t> occurrence;
        std::map<int, std::vector<std::size_t>> original_indices;
        for (std::size_t i = 0U; i < labels.size(); ++i) {
            original_indices[labels[i]].push_back(i);
        }
        std::vector<std::size_t> order;
        order.reserve(labels.size() + 1U);
        for (int degree : degree_order) {
            const std::size_t index = occurrence[degree];
            order.push_back(original_indices.at(degree)[index]);
            ++occurrence[degree];
        }
        order.push_back(labels.size());
        last = analyze_gt_case(labels, 0, minus_mask, cache, order);
        if (last.fiber_parity) {
            return last;
        }
    } while (std::next_permutation(degree_order.begin(), degree_order.end()));
    return last;
}

bool support_disjoint(
    const std::vector<int>& labels,
    unsigned int minus_mask
) {
    std::map<int, int> signs;
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        const int sign = ((minus_mask >> i) & 1U) != 0U ? -1 : 1;
        const auto [position, inserted] = signs.emplace(labels[i], sign);
        if (!inserted && position->second != sign) {
            return false;
        }
    }
    return true;
}

std::vector<int> signed_vertex_classes(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask
) {
    std::vector<std::pair<int, int>> keys;
    keys.reserve(labels.size() + 1U);
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        const int sign = ((minus_mask >> i) & 1U) != 0U ? -1 : 1;
        keys.emplace_back(labels[i], sign);
    }
    const int target_sign = (popcount(minus_mask) & 1) != 0 ? -1 : 1;
    keys.emplace_back(target, target_sign);
    std::map<std::pair<int, int>, int> identifiers;
    for (const std::pair<int, int>& key : keys) {
        if (!identifiers.contains(key)) {
            identifiers.emplace(
                key, static_cast<int>(identifiers.size())
            );
        }
    }
    std::vector<int> classes;
    classes.reserve(keys.size());
    for (const std::pair<int, int>& key : keys) {
        classes.push_back(identifiers.at(key));
    }
    return classes;
}

std::vector<std::size_t> representative_vertex_order(
    const std::vector<int>& vertex_classes,
    const std::vector<int>& class_order
) {
    if (vertex_classes.size() != class_order.size()) {
        throw std::runtime_error("class order has the wrong size");
    }
    int maximum_class = -1;
    for (int vertex_class : vertex_classes) {
        maximum_class = std::max(maximum_class, vertex_class);
    }
    std::vector<std::vector<std::size_t>> members(
        static_cast<std::size_t>(maximum_class + 1)
    );
    for (std::size_t i = 0U; i < vertex_classes.size(); ++i) {
        members[static_cast<std::size_t>(vertex_classes[i])].push_back(i);
    }
    std::vector<std::size_t> used(members.size(), 0U);
    std::vector<std::size_t> order;
    order.reserve(class_order.size());
    for (int vertex_class : class_order) {
        if (vertex_class < 0
            || static_cast<std::size_t>(vertex_class) >= members.size()) {
            throw std::runtime_error("class order is out of range");
        }
        const std::size_t class_index
            = static_cast<std::size_t>(vertex_class);
        if (used[class_index] >= members[class_index].size()) {
            throw std::runtime_error("class order has excess multiplicity");
        }
        order.push_back(members[class_index][used[class_index]]);
        ++used[class_index];
    }
    return order;
}

Counts analyze_gt_best_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::vector<std::size_t> order(labels.size() + 1U, 0U);
    for (std::size_t i = 0U; i < order.size(); ++i) {
        order[i] = i;
    }
    Counts best;
    bool initialized = false;
    do {
        const Counts counts = analyze_gt_case(
            labels, target, minus_mask, cache, order
        );
        if (!initialized || counts.odd_rank > best.odd_rank) {
            best = counts;
            initialized = true;
        }
        const std::uint64_t kernel = counts.odd_sources
            - static_cast<std::uint64_t>(counts.odd_rank);
        if (kernel <= counts.positive_sources) {
            return counts;
        }
    } while (std::next_permutation(order.begin(), order.end()));
    return best;
}

Counts analyze_gt_fiber_best_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::vector<std::size_t> order(labels.size() + 1U, 0U);
    for (std::size_t i = 0U; i < order.size(); ++i) {
        order[i] = i;
    }
    Counts last;
    do {
        last = analyze_gt_case(labels, target, minus_mask, cache, order);
        if (last.fiber_parity) {
            return last;
        }
    } while (std::next_permutation(order.begin(), order.end()));
    return last;
}

Counts analyze_gt_deficit_fiber_best_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    Counts first = analyze_gt_case(labels, target, minus_mask, cache);
    if (first.odd_sources <= first.positive_sources) {
        first.fiber_parity = true;
        return first;
    }
    return analyze_gt_fiber_best_case(labels, target, minus_mask, cache);
}

Counts analyze_gt_xor_best_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::vector<std::size_t> order(labels.size() + 1U, 0U);
    for (std::size_t i = 0U; i < order.size(); ++i) {
        order[i] = i;
    }
    Counts last;
    do {
        last = analyze_gt_case(labels, target, minus_mask, cache, order);
        if (last.xor_payment) {
            return last;
        }
    } while (std::next_permutation(order.begin(), order.end()));
    return last;
}

Counts analyze_gt_monotone_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::vector<std::size_t> order(labels.size() + 1U, 0U);
    for (std::size_t i = 0U; i < order.size(); ++i) {
        order[i] = i;
    }
    Counts first = analyze_gt_case(
        labels, target, minus_mask, cache, order
    );
    if (first.odd_sources <= first.positive_sources) {
        first.xor_payment = true;
        return first;
    }
    const std::size_t required_rank = static_cast<std::size_t>(
        first.odd_sources - first.positive_sources
    );
    const std::vector<int> vertex_classes = signed_vertex_classes(
        labels, target, minus_mask
    );
    std::vector<int> class_order = vertex_classes;
    std::sort(class_order.begin(), class_order.end());
    std::map<std::vector<int>, std::size_t> landscape;
    do {
        order = representative_vertex_order(vertex_classes, class_order);
        const Counts counts = analyze_gt_case(
            labels, target, minus_mask, cache, order
        );
        landscape.emplace(class_order, counts.odd_rank);
    } while (std::next_permutation(class_order.begin(), class_order.end()));
    std::set<std::vector<int>> reachable;
    std::queue<std::vector<int>> frontier;
    for (const auto& [candidate, rank] : landscape) {
        if (rank >= required_rank) {
            reachable.insert(candidate);
            frontier.push(candidate);
        }
    }
    while (!frontier.empty()) {
        const std::vector<int> current = frontier.front();
        frontier.pop();
        const std::size_t current_rank = landscape.at(current);
        for (std::size_t i = 0U; i + 1U < current.size(); ++i) {
            if (current[i] == current[i + 1U]) {
                continue;
            }
            std::vector<int> predecessor = current;
            std::swap(predecessor[i], predecessor[i + 1U]);
            if (reachable.contains(predecessor)
                || landscape.at(predecessor) > current_rank) {
                continue;
            }
            reachable.insert(predecessor);
            frontier.push(predecessor);
        }
    }
    first.xor_payment = reachable.size() == landscape.size();
    if (!first.xor_payment) {
        first.xor_demand = static_cast<std::uint64_t>(
            landscape.size() - reachable.size()
        );
        first.xor_capacity = static_cast<std::uint64_t>(required_rank);
    }
    return first;
}

Counts analyze_gt_average_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache,
    bool require_pair_bound
) {
    std::vector<std::size_t> order(labels.size() + 1U, 0U);
    for (std::size_t i = 0U; i < order.size(); ++i) {
        order[i] = i;
    }
    Counts first = analyze_gt_case(
        labels, target, minus_mask, cache, order
    );
    if (first.odd_sources <= first.positive_sources) {
        first.xor_payment = true;
        return first;
    }
    const std::vector<int> vertex_classes = signed_vertex_classes(
        labels, target, minus_mask
    );
    std::vector<int> class_order = vertex_classes;
    std::sort(class_order.begin(), class_order.end());
    do {
        order = representative_vertex_order(vertex_classes, class_order);
        const Counts counts = analyze_gt_case(
            labels, target, minus_mask, cache, order
        );
        if (first.rank_sum
            > std::numeric_limits<std::uint64_t>::max()
                - static_cast<std::uint64_t>(counts.odd_rank)) {
            throw std::runtime_error("GT average rank sum overflow");
        }
        first.rank_sum += static_cast<std::uint64_t>(counts.odd_rank);
        if (first.pair_collision_sum
            > std::numeric_limits<std::uint64_t>::max()
                - counts.odd_pair_collisions) {
            throw std::runtime_error("GT average pair sum overflow");
        }
        first.pair_collision_sum += counts.odd_pair_collisions;
        ++first.order_count;
    } while (std::next_permutation(class_order.begin(), class_order.end()));
    const std::uint64_t required_rank
        = first.odd_sources - first.positive_sources;
    if (required_rank != 0U
        && first.order_count
            > std::numeric_limits<std::uint64_t>::max() / required_rank) {
        throw std::runtime_error("GT average threshold overflow");
    }
    const std::uint64_t required_sum = first.order_count * required_rank;
    if (first.positive_sources != 0U
        && first.order_count > std::numeric_limits<std::uint64_t>::max()
            / first.positive_sources) {
        throw std::runtime_error("GT average positive sum overflow");
    }
    const std::uint64_t positive_sum
        = first.order_count * first.positive_sources;
    first.xor_payment = require_pair_bound
        ? first.pair_collision_sum <= positive_sum
        : first.rank_sum >= required_sum;
    if (!first.xor_payment) {
        first.xor_demand = require_pair_bound
            ? first.pair_collision_sum : required_sum;
        first.xor_capacity = require_pair_bound
            ? positive_sum : first.rank_sum;
    }
    return first;
}

Counts analyze_gt_weighted_pair_average_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::vector<std::size_t> order(labels.size() + 1U, 0U);
    std::iota(order.begin(), order.end(), 0U);
    Counts first = analyze_gt_case(
        labels, target, minus_mask, cache, order
    );
    if (first.odd_sources <= first.positive_sources) {
        first.xor_payment = true;
        return first;
    }

    const std::vector<int> vertex_classes = signed_vertex_classes(
        labels, target, minus_mask
    );
    std::vector<int> class_order = vertex_classes;
    std::sort(class_order.begin(), class_order.end());
    WeightedPairCounts weighted_pairs;
    std::map<unsigned int, std::uint64_t> capacities;
    do {
        order = representative_vertex_order(vertex_classes, class_order);
        (void)analyze_gt_case(
            labels, target, minus_mask, cache, order, nullptr, nullptr,
            &weighted_pairs, &capacities
        );
        ++first.order_count;
    } while (std::next_permutation(class_order.begin(), class_order.end()));

    std::uint64_t maximum_odd = 1U;
    for (const auto& [key, count] : weighted_pairs) {
        (void)count;
        maximum_odd = std::max(maximum_odd, key.second);
    }
    boost::multiprecision::cpp_int common_denominator = 1;
    for (std::uint64_t divisor = 2U;
         divisor <= maximum_odd; ++divisor) {
        const std::uint64_t remainder
            = (common_denominator % divisor).convert_to<std::uint64_t>();
        const std::uint64_t common = std::gcd(divisor, remainder);
        common_denominator
            = common_denominator / common * divisor;
    }

    std::map<unsigned int, boost::multiprecision::cpp_int> demands;
    for (const auto& [key, count] : weighted_pairs) {
        const unsigned int mask = key.first;
        const std::uint64_t odd = key.second;
        demands[mask] += boost::multiprecision::cpp_int(count) * 2
            * (common_denominator / odd);
    }
    first.xor_payment = true;
    for (const auto& [mask, demand] : demands) {
        const boost::multiprecision::cpp_int capacity
            = boost::multiprecision::cpp_int(capacities[mask])
                * common_denominator;
        if (demand > capacity) {
            first.xor_payment = false;
            first.xor_mask = mask;
            break;
        }
    }
    return first;
}

struct TreePairEdge {
    std::size_t fiber = 0U;
    unsigned int left = 0U;
    unsigned int right = 0U;
    unsigned int channel = 0U;
};

struct PositiveResource {
    std::size_t fiber = 0U;
    unsigned int channel = 0U;
};

std::vector<std::size_t> graphic_pair_circuit(
    const std::vector<TreePairEdge>& edges,
    const std::vector<bool>& selected,
    std::size_t candidate
) {
    const TreePairEdge& edge = edges[candidate];
    std::map<unsigned int, std::pair<unsigned int, std::size_t>> parent;
    std::queue<unsigned int> frontier;
    parent.emplace(edge.left, std::pair{edge.left, edges.size()});
    frontier.push(edge.left);
    while (!frontier.empty() && !parent.contains(edge.right)) {
        const unsigned int current = frontier.front();
        frontier.pop();
        for (std::size_t i = 0U; i < edges.size(); ++i) {
            if (!selected[i] || edges[i].fiber != edge.fiber) {
                continue;
            }
            unsigned int neighbor = 0U;
            if (edges[i].left == current) {
                neighbor = edges[i].right;
            } else if (edges[i].right == current) {
                neighbor = edges[i].left;
            } else {
                continue;
            }
            if (!parent.contains(neighbor)) {
                parent.emplace(neighbor, std::pair{current, i});
                frontier.push(neighbor);
            }
        }
    }
    std::vector<std::size_t> circuit;
    if (!parent.contains(edge.right)) {
        return circuit;
    }
    unsigned int current = edge.right;
    while (current != edge.left) {
        const auto [previous, path_edge] = parent.at(current);
        circuit.push_back(path_edge);
        current = previous;
    }
    return circuit;
}

std::vector<std::ptrdiff_t> transversal_pair_matching(
    const std::vector<std::vector<std::size_t>>& allowed,
    const std::vector<bool>& selected,
    std::size_t resource_count
) {
    std::vector<std::ptrdiff_t> resource_match(resource_count, -1);
    for (std::size_t element = 0U; element < selected.size(); ++element) {
        if (!selected[element]) {
            continue;
        }
        std::vector<bool> seen(resource_count, false);
        std::function<bool(std::size_t)> augment
            = [&](std::size_t current) -> bool {
                for (std::size_t resource : allowed[current]) {
                    if (seen[resource]) {
                        continue;
                    }
                    seen[resource] = true;
                    if (resource_match[resource] < 0
                        || augment(static_cast<std::size_t>(
                            resource_match[resource]
                        ))) {
                        resource_match[resource]
                            = static_cast<std::ptrdiff_t>(current);
                        return true;
                    }
                }
                return false;
            };
        if (!augment(element)) {
            throw std::runtime_error(
                "selected tree pairs have no transversal matching"
            );
        }
    }
    return resource_match;
}

std::pair<bool, std::vector<std::size_t>> transversal_pair_circuit(
    const std::vector<std::vector<std::size_t>>& allowed,
    const std::vector<std::ptrdiff_t>& resource_match,
    std::size_t candidate
) {
    std::vector<bool> seen_resource(resource_match.size(), false);
    std::vector<bool> seen_element(allowed.size(), false);
    std::queue<std::size_t> frontier;
    for (std::size_t resource : allowed[candidate]) {
        if (!seen_resource[resource]) {
            seen_resource[resource] = true;
            frontier.push(resource);
        }
    }
    while (!frontier.empty()) {
        const std::size_t resource = frontier.front();
        frontier.pop();
        if (resource_match[resource] < 0) {
            return {true, {}};
        }
        const std::size_t element
            = static_cast<std::size_t>(resource_match[resource]);
        if (seen_element[element]) {
            continue;
        }
        seen_element[element] = true;
        for (std::size_t next : allowed[element]) {
            if (!seen_resource[next]) {
                seen_resource[next] = true;
                frontier.push(next);
            }
        }
    }
    std::vector<std::size_t> circuit;
    for (std::size_t element = 0U; element < seen_element.size(); ++element) {
        if (seen_element[element]) {
            circuit.push_back(element);
        }
    }
    return {false, circuit};
}

bool graphic_pair_addable(
    const std::vector<TreePairEdge>& edges,
    const std::vector<bool>& selected,
    std::size_t candidate,
    std::size_t removed
) {
    const TreePairEdge& edge = edges[candidate];
    std::set<unsigned int> reached{edge.left};
    std::queue<unsigned int> frontier;
    frontier.push(edge.left);
    while (!frontier.empty()) {
        const unsigned int current = frontier.front();
        frontier.pop();
        for (std::size_t i = 0U; i < edges.size(); ++i) {
            if (!selected[i] || i == removed
                || edges[i].fiber != edge.fiber) {
                continue;
            }
            unsigned int neighbor = 0U;
            if (edges[i].left == current) {
                neighbor = edges[i].right;
            } else if (edges[i].right == current) {
                neighbor = edges[i].left;
            } else {
                continue;
            }
            if (reached.insert(neighbor).second) {
                frontier.push(neighbor);
            }
        }
    }
    return !reached.contains(edge.right);
}

Counts analyze_gt_tree_pair_average_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::vector<std::size_t> order(labels.size() + 1U, 0U);
    std::iota(order.begin(), order.end(), 0U);
    Counts first = analyze_gt_case(
        labels, target, minus_mask, cache, order
    );
    if (first.odd_sources <= first.positive_sources) {
        first.xor_payment = true;
        return first;
    }

    const std::vector<int> vertex_classes = signed_vertex_classes(
        labels, target, minus_mask
    );
    std::vector<int> class_order = vertex_classes;
    std::sort(class_order.begin(), class_order.end());
    std::vector<std::vector<unsigned int>> odd_fibers;
    std::map<unsigned int, std::uint64_t> capacities;
    do {
        order = representative_vertex_order(vertex_classes, class_order);
        (void)analyze_gt_case(
            labels, target, minus_mask, cache, order, nullptr, nullptr,
            nullptr, &capacities, &odd_fibers
        );
        ++first.order_count;
    } while (std::next_permutation(class_order.begin(), class_order.end()));

    std::vector<TreePairEdge> edges;
    std::uint64_t target_rank = 0U;
    for (std::size_t fiber = 0U; fiber < odd_fibers.size(); ++fiber) {
        const std::vector<unsigned int>& cuts = odd_fibers[fiber];
        target_rank += static_cast<std::uint64_t>(cuts.size() - 1U);
        for (std::size_t left = 0U; left < cuts.size(); ++left) {
            for (std::size_t right = left + 1U;
                 right < cuts.size(); ++right) {
                edges.push_back(TreePairEdge{
                    fiber, cuts[left], cuts[right], cuts[left] ^ cuts[right]
                });
            }
        }
    }

    std::vector<bool> selected(edges.size(), false);
    std::map<unsigned int, std::uint64_t> used;
    std::uint64_t cardinality = 0U;
    while (cardinality < target_rank) {
        const std::ptrdiff_t unvisited = -2;
        std::vector<std::ptrdiff_t> parent(edges.size(), unvisited);
        std::queue<std::size_t> frontier;
        for (std::size_t edge = 0U; edge < edges.size(); ++edge) {
            if (!selected[edge]
                && graphic_pair_addable(
                    edges, selected, edge, edges.size()
                )) {
                parent[edge] = -1;
                frontier.push(edge);
            }
        }

        std::size_t sink = edges.size();
        while (!frontier.empty() && sink == edges.size()) {
            const std::size_t current = frontier.front();
            frontier.pop();
            if (!selected[current]) {
                const unsigned int channel = edges[current].channel;
                if (used[channel] < capacities[channel]) {
                    sink = current;
                    break;
                }
                for (std::size_t edge = 0U; edge < edges.size(); ++edge) {
                    if (selected[edge] && parent[edge] == unvisited
                        && edges[edge].channel == channel) {
                        parent[edge] = static_cast<std::ptrdiff_t>(current);
                        frontier.push(edge);
                    }
                }
            } else {
                for (std::size_t edge = 0U; edge < edges.size(); ++edge) {
                    if (!selected[edge] && parent[edge] == unvisited
                        && graphic_pair_addable(
                            edges, selected, edge, current
                        )) {
                        parent[edge] = static_cast<std::ptrdiff_t>(current);
                        frontier.push(edge);
                    }
                }
            }
        }
        if (sink == edges.size()) {
            break;
        }
        std::ptrdiff_t current = static_cast<std::ptrdiff_t>(sink);
        while (current >= 0) {
            const std::size_t edge = static_cast<std::size_t>(current);
            selected[edge] = !selected[edge];
            current = parent[edge];
        }
        used.clear();
        cardinality = 0U;
        for (std::size_t edge = 0U; edge < edges.size(); ++edge) {
            if (selected[edge]) {
                ++used[edges[edge].channel];
                ++cardinality;
            }
        }
    }
    first.xor_payment = cardinality == target_rank;
    if (!first.xor_payment) {
        first.xor_demand = target_rank;
        first.xor_capacity = cardinality;
    }
    return first;
}

Counts analyze_gt_transversal_pair_average_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache,
    bool boolean_channels = false,
    bool universal_channels = false,
    bool fiber_boolean_channels = false
) {
    std::vector<std::size_t> order(labels.size() + 1U, 0U);
    std::iota(order.begin(), order.end(), 0U);
    Counts first = analyze_gt_case(
        labels, target, minus_mask, cache, order
    );
    if (first.odd_sources <= first.positive_sources) {
        first.xor_payment = true;
        return first;
    }

    const std::vector<int> vertex_classes = signed_vertex_classes(
        labels, target, minus_mask
    );
    std::vector<int> class_order = vertex_classes;
    std::sort(class_order.begin(), class_order.end());
    std::vector<CutFiberData> cut_fibers;
    do {
        order = representative_vertex_order(vertex_classes, class_order);
        (void)analyze_gt_case(
            labels, target, minus_mask, cache, order, nullptr, nullptr,
            nullptr, nullptr, nullptr, &cut_fibers
        );
        ++first.order_count;
    } while (std::next_permutation(class_order.begin(), class_order.end()));

    std::vector<PositiveResource> resources;
    std::map<std::size_t, std::vector<std::size_t>> resources_by_fiber;
    std::map<unsigned int, std::vector<std::size_t>> resources_by_channel;
    for (std::size_t fiber = 0U; fiber < cut_fibers.size(); ++fiber) {
        for (unsigned int channel : cut_fibers[fiber].positive_cuts) {
            const std::size_t resource = resources.size();
            resources.push_back(PositiveResource{fiber, channel});
            resources_by_fiber[fiber].push_back(resource);
            resources_by_channel[channel].push_back(resource);
        }
    }

    std::vector<TreePairEdge> edges;
    std::uint64_t target_rank = 0U;
    for (std::size_t fiber = 0U; fiber < cut_fibers.size(); ++fiber) {
        const std::vector<unsigned int>& cuts = cut_fibers[fiber].odd_cuts;
        if (cuts.size() < 2U) {
            continue;
        }
        target_rank += static_cast<std::uint64_t>(cuts.size() - 1U);
        for (std::size_t left = 0U; left < cuts.size(); ++left) {
            for (std::size_t right = left + 1U;
                 right < cuts.size(); ++right) {
                edges.push_back(TreePairEdge{
                    fiber, cuts[left], cuts[right], cuts[left] ^ cuts[right]
                });
            }
        }
    }

    std::vector<std::vector<std::size_t>> allowed(edges.size());
    for (std::size_t edge = 0U; edge < edges.size(); ++edge) {
        allowed[edge] = resources_by_fiber[edges[edge].fiber];
        std::vector<unsigned int> channels{edges[edge].channel};
        if (universal_channels) {
            for (const auto& [channel, channel_resources]
                 : resources_by_channel) {
                (void)channel_resources;
                channels.push_back(channel);
            }
        } else if (fiber_boolean_channels) {
            const unsigned int universe
                = (1U << labels.size()) - 1U;
            std::vector<unsigned int> atoms{universe};
            for (unsigned int cut : cut_fibers[edges[edge].fiber].odd_cuts) {
                std::vector<unsigned int> refined;
                for (unsigned int atom : atoms) {
                    if ((atom & cut) != 0U) {
                        refined.push_back(atom & cut);
                    }
                    if ((atom & ~cut & universe) != 0U) {
                        refined.push_back(atom & ~cut & universe);
                    }
                }
                atoms = std::move(refined);
            }
            if (atoms.size() >= std::numeric_limits<unsigned int>::digits) {
                throw std::runtime_error(
                    "too many fibre Boolean atoms"
                );
            }
            const unsigned int atom_limit = 1U << atoms.size();
            for (unsigned int atom_mask = 1U;
                 atom_mask < atom_limit; ++atom_mask) {
                unsigned int channel = 0U;
                for (std::size_t atom = 0U; atom < atoms.size(); ++atom) {
                    if (((atom_mask >> atom) & 1U) != 0U) {
                        channel |= atoms[atom];
                    }
                }
                channels.push_back(channel);
            }
        } else if (boolean_channels) {
            const unsigned int universe
                = (1U << labels.size()) - 1U;
            const std::array<unsigned int, 4U> atoms{
                edges[edge].left & edges[edge].right,
                edges[edge].left & ~edges[edge].right,
                edges[edge].right & ~edges[edge].left,
                universe & ~(edges[edge].left | edges[edge].right)
            };
            for (unsigned int atom_mask = 1U;
                 atom_mask < 16U; ++atom_mask) {
                unsigned int channel = 0U;
                for (std::size_t atom = 0U; atom < atoms.size(); ++atom) {
                    if (((atom_mask >> atom) & 1U) != 0U) {
                        channel |= atoms[atom];
                    }
                }
                channels.push_back(channel);
            }
        }
        for (unsigned int channel : channels) {
            if (channel == 0U
                || (popcount(channel & minus_mask) & 1) != 0) {
                continue;
            }
            const std::vector<std::size_t>& channel_resources
                = resources_by_channel[channel];
            allowed[edge].insert(
                allowed[edge].end(),
                channel_resources.begin(), channel_resources.end()
            );
        }
        std::sort(allowed[edge].begin(), allowed[edge].end());
        allowed[edge].erase(
            std::unique(allowed[edge].begin(), allowed[edge].end()),
            allowed[edge].end()
        );
    }

    std::vector<bool> selected(edges.size(), false);
    std::uint64_t cardinality = 0U;
    while (cardinality < target_rank) {
        const std::vector<std::ptrdiff_t> resource_match
            = transversal_pair_matching(
                allowed, selected, resources.size()
            );
        std::vector<std::vector<std::size_t>> graphic_circuits(edges.size());
        std::vector<std::vector<std::size_t>> transversal_circuits(
            edges.size()
        );
        std::vector<bool> transversal_addable(edges.size(), false);
        for (std::size_t edge = 0U; edge < edges.size(); ++edge) {
            if (selected[edge]) {
                continue;
            }
            graphic_circuits[edge]
                = graphic_pair_circuit(edges, selected, edge);
            auto [addable, circuit] = transversal_pair_circuit(
                allowed, resource_match, edge
            );
            transversal_addable[edge] = addable;
            transversal_circuits[edge] = std::move(circuit);
        }

        const std::ptrdiff_t unvisited = -2;
        std::vector<std::ptrdiff_t> parent(edges.size(), unvisited);
        std::queue<std::size_t> frontier;
        for (std::size_t edge = 0U; edge < edges.size(); ++edge) {
            if (!selected[edge] && graphic_circuits[edge].empty()) {
                parent[edge] = -1;
                frontier.push(edge);
            }
        }
        std::size_t sink = edges.size();
        while (!frontier.empty() && sink == edges.size()) {
            const std::size_t current = frontier.front();
            frontier.pop();
            if (!selected[current]) {
                if (transversal_addable[current]) {
                    sink = current;
                    break;
                }
                for (std::size_t next : transversal_circuits[current]) {
                    if (parent[next] == unvisited) {
                        parent[next] = static_cast<std::ptrdiff_t>(current);
                        frontier.push(next);
                    }
                }
            } else {
                for (std::size_t next = 0U; next < edges.size(); ++next) {
                    if (selected[next] || parent[next] != unvisited) {
                        continue;
                    }
                    if (std::find(
                            graphic_circuits[next].begin(),
                            graphic_circuits[next].end(), current
                        ) != graphic_circuits[next].end()) {
                        parent[next] = static_cast<std::ptrdiff_t>(current);
                        frontier.push(next);
                    }
                }
            }
        }
        if (sink == edges.size()) {
            break;
        }
        std::ptrdiff_t current = static_cast<std::ptrdiff_t>(sink);
        while (current >= 0) {
            const std::size_t edge = static_cast<std::size_t>(current);
            selected[edge] = !selected[edge];
            current = parent[edge];
        }
        ++cardinality;
    }

    first.xor_payment = cardinality == target_rank;
    if (!first.xor_payment) {
        first.xor_demand = target_rank;
        first.xor_capacity = cardinality;
    }
    return first;
}

Counts analyze_gt_alpha_average_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    std::vector<std::size_t> order(labels.size() + 1U, 0U);
    for (std::size_t i = 0U; i < order.size(); ++i) {
        order[i] = i;
    }
    Counts first = analyze_gt_case(
        labels, target, minus_mask, cache, order
    );
    if (first.odd_sources <= first.positive_sources) {
        first.xor_payment = true;
        return first;
    }
    std::map<Key, AlphaAverage> averages;
    do {
        (void)analyze_gt_case(
            labels, target, minus_mask, cache, order, &averages
        );
        ++first.order_count;
    } while (std::next_permutation(order.begin(), order.end()));
    first.xor_payment = true;
    for (const auto& [raw_top, average] : averages) {
        (void)raw_top;
        if (average.odd_occurrences < average.hit_orders) {
            throw std::runtime_error("GT alpha average count underflow");
        }
        const std::uint64_t relations
            = average.odd_occurrences - average.hit_orders;
        first.rank_sum += average.hit_orders;
        if (relations > average.positive_occurrences) {
            first.xor_payment = false;
            first.xor_demand = relations;
            first.xor_capacity = average.positive_occurrences;
            first.alpha_witness = raw_top;
            break;
        }
    }
    return first;
}

using StateKey = std::pair<std::vector<std::size_t>, Key>;

Key bender_knuth_raw_neighbor(
    const std::vector<int>& degrees,
    const StateKey& state,
    std::size_t first
) {
    const std::vector<std::size_t>& order = state.first;
    const Key& raw_top = state.second;
    if (first + 1U >= order.size()) {
        throw std::runtime_error("GT state swap is out of range");
    }
    int height = 0;
    for (std::size_t position = 0U; position < first; ++position) {
        const std::size_t vertex = order[position];
        height += 2 * raw_top[vertex] - degrees[vertex];
    }
    const std::size_t first_vertex = order[first];
    const std::size_t second_vertex = order[first + 1U];
    const int first_bottom
        = degrees[first_vertex] - raw_top[first_vertex];
    const int second_bottom
        = degrees[second_vertex] - raw_top[second_vertex];
    const int transferred = std::max(
        0, first_bottom + second_bottom - height
    );
    Key neighbor = raw_top;
    neighbor[first_vertex] -= transferred;
    neighbor[second_vertex] += transferred;
    return neighbor;
}

Counts analyze_gt_bk_component_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache,
    bool print_components = false
) {
    std::vector<std::size_t> order(labels.size() + 1U, 0U);
    for (std::size_t i = 0U; i < order.size(); ++i) {
        order[i] = i;
    }
    Counts first = analyze_gt_case(
        labels, target, minus_mask, cache, order
    );
    if (first.odd_sources <= first.positive_sources) {
        first.xor_payment = true;
        return first;
    }
    std::map<StateKey, FiberCounts> states;
    do {
        std::map<Key, FiberCounts> fibers;
        (void)analyze_gt_case(
            labels, target, minus_mask, cache, order, nullptr, &fibers
        );
        for (const auto& [raw_top, counts] : fibers) {
            states.emplace(StateKey{order, raw_top}, counts);
        }
        ++first.order_count;
    } while (std::next_permutation(order.begin(), order.end()));
    first.rank_sum = static_cast<std::uint64_t>(states.size());

    std::vector<int> degrees = labels;
    degrees.push_back(target);
    std::set<StateKey> visited;
    first.xor_payment = true;
    for (const auto& [initial, initial_counts] : states) {
        (void)initial_counts;
        if (visited.contains(initial)) {
            continue;
        }
        std::queue<StateKey> frontier;
        visited.insert(initial);
        frontier.push(initial);
        std::uint64_t relations = 0U;
        std::uint64_t positive = 0U;
        std::uint64_t odd_occurrences = 0U;
        std::uint64_t active_states = 0U;
        std::uint64_t component_states = 0U;
        while (!frontier.empty()) {
            const StateKey current = frontier.front();
            frontier.pop();
            ++component_states;
            const FiberCounts counts = states.at(current);
            relations += counts.odd
                - static_cast<std::uint64_t>(counts.odd != 0U);
            odd_occurrences += counts.odd;
            active_states += counts.odd != 0U ? 1U : 0U;
            positive += counts.positive;
            for (std::size_t i = 0U; i + 1U < current.first.size(); ++i) {
                StateKey neighbor{
                    current.first,
                    bender_knuth_raw_neighbor(degrees, current, i)
                };
                std::swap(neighbor.first[i], neighbor.first[i + 1U]);
                if (!states.contains(neighbor)) {
                    throw std::runtime_error(
                        "Bender--Knuth neighbor is missing"
                    );
                }
                if (visited.insert(neighbor).second) {
                    frontier.push(std::move(neighbor));
                }
            }
        }
        ++first.component_count;
        if (print_components) {
            std::cout << "component=" << first.component_count
                      << " states=" << component_states
                      << " odd=" << odd_occurrences
                      << " active_states=" << active_states
                      << " relations=" << relations
                      << " positive=" << positive
                      << " slack="
                      << (positive >= relations
                          ? positive - relations : 0U)
                      << '\n';
        }
        if (odd_occurrences != 0U || positive != 0U) {
            ++first.active_component_count;
        }
        if (relations != 0U || positive != 0U) {
            ++first.loaded_component_count;
        }
        const std::uint64_t slack
            = positive >= relations ? positive - relations : 0U;
        if (first.component_count == 1U
            || slack < first.minimum_component_slack) {
            first.minimum_component_states = component_states;
            first.minimum_component_relations = relations;
            first.minimum_component_positive = positive;
            first.minimum_component_slack = slack;
        }
        if (relations > positive) {
            first.xor_payment = false;
            first.xor_demand = relations;
            first.xor_capacity = positive;
            first.alpha_witness = initial.second;
            break;
        }
    }
    return first;
}

Counts analyze_gt_bk_single_active_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    Counts counts = analyze_gt_bk_component_case(
        labels, target, minus_mask, cache
    );
    counts.xor_payment = counts.active_component_count <= 1U;
    if (!counts.xor_payment) {
        counts.xor_demand = counts.active_component_count;
        counts.xor_capacity = 1U;
    }
    return counts;
}

Counts analyze_gt_bk_single_loaded_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    ToricCache& cache
) {
    Counts counts = analyze_gt_bk_component_case(
        labels, target, minus_mask, cache
    );
    counts.xor_payment = counts.loaded_component_count <= 1U;
    if (!counts.xor_payment) {
        counts.xor_demand = counts.loaded_component_count;
        counts.xor_capacity = 1U;
    }
    return counts;
}

void print_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    const Counts& counts
) {
    std::cout << "labels=[";
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        if (i != 0U) {
            std::cout << ',';
        }
        std::cout << labels[i];
    }
    std::cout << "] target=" << target
              << " minus_mask=" << minus_mask
              << " odd_sources=" << counts.odd_sources
              << " odd_rank=" << counts.odd_rank
              << " degenerate_kernel="
              << counts.odd_sources
                    - static_cast<std::uint64_t>(counts.odd_rank)
              << " positive_sources=" << counts.positive_sources
              << " odd_zero=" << counts.odd_zero_sources
              << " positive_zero=" << counts.positive_zero_sources
              << " fiber_parity="
              << (counts.fiber_parity ? "PASS" : "FAIL")
              << " xor_payment="
              << (counts.xor_payment ? "PASS" : "FAIL");
    if (!counts.fiber_parity) {
        std::cout << " fiber_odd=" << counts.fiber_odd
                  << " fiber_even=" << counts.fiber_even;
    }
    if (!counts.xor_payment) {
        std::cout << " xor_mask=" << counts.xor_mask
                  << " xor_demand=" << counts.xor_demand
                  << " xor_capacity=" << counts.xor_capacity;
    }
    if (counts.order_count != 0U) {
        std::cout << " orders=" << counts.order_count
                  << " rank_sum=" << counts.rank_sum;
    }
    if (counts.pair_collision_sum != 0U) {
        std::cout << " pair_collision_sum="
                  << counts.pair_collision_sum;
    }
    if (!counts.alpha_witness.empty()) {
        std::cout << " alpha_witness=[";
        for (std::size_t i = 0U; i < counts.alpha_witness.size(); ++i) {
            if (i != 0U) {
                std::cout << ',';
            }
            std::cout << counts.alpha_witness[i];
        }
        std::cout << ']';
    }
    if (counts.component_count != 0U) {
        std::cout << " components=" << counts.component_count
                  << " active_components="
                  << counts.active_component_count
                  << " loaded_components="
                  << counts.loaded_component_count
                  << " minimum_component_states="
                  << counts.minimum_component_states
                  << " minimum_component_relations="
                  << counts.minimum_component_relations
                  << " minimum_component_positive="
                  << counts.minimum_component_positive
                  << " minimum_component_slack="
                  << counts.minimum_component_slack;
    }
}

bool increment_labels(std::vector<int>& labels, int maximum_label) {
    for (std::size_t reverse = 0U; reverse < labels.size(); ++reverse) {
        const std::size_t index = labels.size() - 1U - reverse;
        if (labels[index] < maximum_label) {
            const int next = labels[index] + 1;
            for (std::size_t j = index; j < labels.size(); ++j) {
                labels[j] = next;
            }
            return true;
        }
    }
    return false;
}

int parse_nonnegative(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value < 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc >= 6 && std::string(argv[1]) == "gt-average-case") {
            const std::string method = argv[2];
            const int target = parse_nonnegative(argv[3], "target");
            const int parsed_mask = parse_nonnegative(argv[4], "minus mask");
            std::vector<int> labels;
            for (int argument = 5; argument < argc; ++argument) {
                labels.push_back(parse_nonnegative(argv[argument], "label"));
            }
            if (labels.empty()
                || labels.size() >= std::numeric_limits<unsigned int>::digits
                || parsed_mask >= static_cast<int>(1U << labels.size())) {
                throw std::runtime_error("invalid gt-average-case data");
            }
            const unsigned int minus_mask
                = static_cast<unsigned int>(parsed_mask);
            ToricCache cache;
            Counts counts;
            if (method == "rank") {
                counts = analyze_gt_average_case(
                    labels, target, minus_mask, cache, false
                );
            } else if (method == "weighted") {
                counts = analyze_gt_weighted_pair_average_case(
                    labels, target, minus_mask, cache
                );
            } else if (method == "tree") {
                counts = analyze_gt_tree_pair_average_case(
                    labels, target, minus_mask, cache
                );
            } else if (method == "transversal" || method == "boolean"
                       || method == "universal"
                       || method == "fiber-boolean") {
                counts = analyze_gt_transversal_pair_average_case(
                    labels, target, minus_mask, cache,
                    method == "boolean", method == "universal",
                    method == "fiber-boolean"
                );
            } else {
                throw std::runtime_error(
                    "gt-average-case method must be rank, weighted, tree, "
                    "transversal, boolean, fiber-boolean, or universal"
                );
            }
            print_case(labels, target, minus_mask, counts);
            std::cout << " method=" << method << " result="
                      << (counts.xor_payment ? "PASS" : "FAIL") << '\n';
            return counts.xor_payment ? 0 : 1;
        }
        if (argc >= 5 && std::string(argv[1]) == "gt-bk-case") {
            const int target = parse_nonnegative(argv[2], "target");
            const int parsed_mask = parse_nonnegative(argv[3], "minus mask");
            std::vector<int> labels;
            for (int argument = 4; argument < argc; ++argument) {
                labels.push_back(parse_nonnegative(argv[argument], "label"));
            }
            if (labels.empty()
                || labels.size() >= std::numeric_limits<unsigned int>::digits
                || parsed_mask >= static_cast<int>(1U << labels.size())) {
                throw std::runtime_error("invalid gt-bk-case data");
            }
            const unsigned int minus_mask
                = static_cast<unsigned int>(parsed_mask);
            ToricCache cache;
            const Counts counts = analyze_gt_bk_component_case(
                labels, target, minus_mask, cache, true
            );
            print_case(labels, target, minus_mask, counts);
            std::cout << " result="
                      << (counts.xor_payment ? "PASS" : "FAIL") << '\n';
            return counts.xor_payment ? 0 : 1;
        }
        if (argc >= 5 && std::string(argv[1]) == "gt-case") {
            const int target = parse_nonnegative(argv[2], "target");
            const int parsed_mask = parse_nonnegative(argv[3], "minus mask");
            std::vector<int> labels;
            for (int argument = 4; argument < argc; ++argument) {
                labels.push_back(parse_nonnegative(argv[argument], "label"));
            }
            if (labels.empty()
                || labels.size() >= std::numeric_limits<unsigned int>::digits
                || parsed_mask >= static_cast<int>(1U << labels.size())) {
                throw std::runtime_error("invalid gt-case data");
            }
            const unsigned int minus_mask
                = static_cast<unsigned int>(parsed_mask);
            std::vector<std::size_t> order(labels.size() + 1U, 0U);
            for (std::size_t i = 0U; i < order.size(); ++i) {
                order[i] = i;
            }
            ToricCache cache;
            Counts best;
            std::vector<std::size_t> best_order = order;
            std::vector<std::size_t> first_fiber_order;
            std::vector<std::size_t> first_xor_order;
            std::map<
                std::vector<std::size_t>,
                std::pair<std::size_t, bool>
            > order_landscape;
            bool initialized = false;
            std::uint64_t orders = 0U;
            do {
                const Counts counts = analyze_gt_case(
                    labels, target, minus_mask, cache, order
                );
                ++orders;
                order_landscape.emplace(
                    order,
                    std::pair<std::size_t, bool>{
                        counts.odd_rank, counts.xor_payment
                    }
                );
                if (!initialized || counts.odd_rank > best.odd_rank) {
                    best = counts;
                    best_order = order;
                    initialized = true;
                }
                if (first_fiber_order.empty() && counts.fiber_parity) {
                    first_fiber_order = order;
                }
                if (first_xor_order.empty() && counts.xor_payment) {
                    first_xor_order = order;
                }
            } while (std::next_permutation(order.begin(), order.end()));
            print_case(labels, target, minus_mask, best);
            std::cout << " orders=" << orders << " best_order=[";
            for (std::size_t i = 0U; i < best_order.size(); ++i) {
                if (i != 0U) {
                    std::cout << ',';
                }
                std::cout << best_order[i];
            }
            const std::uint64_t kernel = best.odd_sources
                - static_cast<std::uint64_t>(best.odd_rank);
            std::cout << "] result="
                      << (kernel <= best.positive_sources ? "PASS" : "FAIL")
                      << " fiber_result="
                      << (!first_fiber_order.empty() ? "PASS" : "FAIL")
                      << " first_fiber_order=[";
            for (std::size_t i = 0U; i < first_fiber_order.size(); ++i) {
                if (i != 0U) {
                    std::cout << ',';
                }
                std::cout << first_fiber_order[i];
            }
            std::cout << "]\n";
            std::cout << "xor_result="
                      << (!first_xor_order.empty() ? "PASS" : "FAIL")
                      << " first_xor_order=[";
            for (std::size_t i = 0U; i < first_xor_order.size(); ++i) {
                if (i != 0U) {
                    std::cout << ',';
                }
                std::cout << first_xor_order[i];
            }
            std::cout << "]\n";
            std::uint64_t failing_local_maxima = 0U;
            std::vector<std::size_t> first_failing_local_maximum;
            if (best.odd_sources > best.positive_sources) {
                for (const auto& [candidate, data] : order_landscape) {
                    if (data.second) {
                        continue;
                    }
                    bool strict_ascent = false;
                    for (std::size_t i = 0U; i + 1U < candidate.size(); ++i) {
                        std::vector<std::size_t> neighbor = candidate;
                        std::swap(neighbor[i], neighbor[i + 1U]);
                        if (order_landscape.at(neighbor).first > data.first) {
                            strict_ascent = true;
                            break;
                        }
                    }
                    if (!strict_ascent) {
                        ++failing_local_maxima;
                        if (first_failing_local_maximum.empty()) {
                            first_failing_local_maximum = candidate;
                        }
                    }
                }
            }
            std::cout << "strict_ascent_result="
                      << (failing_local_maxima == 0U ? "PASS" : "FAIL")
                      << " failing_local_maxima=" << failing_local_maxima
                      << " first_failing_local_maximum=[";
            for (std::size_t i = 0U;
                 i < first_failing_local_maximum.size(); ++i) {
                if (i != 0U) {
                    std::cout << ',';
                }
                std::cout << first_failing_local_maximum[i];
            }
            std::cout << "]\n";
            std::set<std::vector<std::size_t>> monotone_reachable;
            std::map<
                std::vector<std::size_t>,
                std::vector<std::size_t>
            > monotone_next;
            std::queue<std::vector<std::size_t>> frontier;
            for (const auto& [candidate, data] : order_landscape) {
                if (data.second) {
                    monotone_reachable.insert(candidate);
                    frontier.push(candidate);
                }
            }
            while (!frontier.empty()) {
                const std::vector<std::size_t> current = frontier.front();
                frontier.pop();
                const std::size_t current_rank
                    = order_landscape.at(current).first;
                for (std::size_t i = 0U; i + 1U < current.size(); ++i) {
                    std::vector<std::size_t> predecessor = current;
                    std::swap(predecessor[i], predecessor[i + 1U]);
                    if (monotone_reachable.contains(predecessor)
                        || order_landscape.at(predecessor).first
                            > current_rank) {
                        continue;
                    }
                    monotone_reachable.insert(predecessor);
                    monotone_next.emplace(predecessor, current);
                    frontier.push(predecessor);
                }
            }
            std::uint64_t monotone_traps = 0U;
            std::vector<std::size_t> first_monotone_trap;
            if (best.odd_sources > best.positive_sources) {
                for (const auto& [candidate, data] : order_landscape) {
                    if (!data.second
                        && !monotone_reachable.contains(candidate)) {
                        ++monotone_traps;
                        if (first_monotone_trap.empty()) {
                            first_monotone_trap = candidate;
                        }
                    }
                }
            }
            std::cout << "hxc_monotone_reachability_result="
                      << (monotone_traps == 0U ? "PASS" : "FAIL")
                      << " monotone_traps=" << monotone_traps
                      << " first_monotone_trap=[";
            for (std::size_t i = 0U; i < first_monotone_trap.size(); ++i) {
                if (i != 0U) {
                    std::cout << ',';
                }
                std::cout << first_monotone_trap[i];
            }
            std::cout << "]\n";
            if (!first_failing_local_maximum.empty()
                && monotone_reachable.contains(
                    first_failing_local_maximum
                )) {
                std::vector<std::size_t> path_order
                    = first_failing_local_maximum;
                std::cout << "hxc_plateau_escape=";
                bool first_step = true;
                while (!order_landscape.at(path_order).second) {
                    const std::vector<std::size_t>& next
                        = monotone_next.at(path_order);
                    std::size_t swap_index = 0U;
                    while (swap_index + 1U < path_order.size()
                           && (path_order[swap_index] != next[swap_index + 1U]
                               || path_order[swap_index + 1U]
                                   != next[swap_index])) {
                        ++swap_index;
                    }
                    if (!first_step) {
                        std::cout << ',';
                    }
                    std::cout << order_landscape.at(path_order).first
                              << ":s" << swap_index;
                    first_step = false;
                    path_order = next;
                }
                std::cout << ',' << order_landscape.at(path_order).first
                          << ":PASS\n";
            }
            const std::size_t required_rank
                = best.odd_sources > best.positive_sources
                    ? static_cast<std::size_t>(
                        best.odd_sources - best.positive_sources
                    )
                    : 0U;
            const auto rank_reachable_to = [&order_landscape](
                std::size_t threshold
            ) {
                std::set<std::vector<std::size_t>> reachable;
                std::queue<std::vector<std::size_t>> rank_frontier;
                for (const auto& [candidate, data] : order_landscape) {
                    if (data.first >= threshold) {
                        reachable.insert(candidate);
                        rank_frontier.push(candidate);
                    }
                }
                while (!rank_frontier.empty()) {
                    const std::vector<std::size_t> current
                        = rank_frontier.front();
                    rank_frontier.pop();
                    const std::size_t current_rank
                        = order_landscape.at(current).first;
                    for (std::size_t i = 0U;
                         i + 1U < current.size(); ++i) {
                        std::vector<std::size_t> predecessor = current;
                        std::swap(predecessor[i], predecessor[i + 1U]);
                        if (reachable.contains(predecessor)
                            || order_landscape.at(predecessor).first
                                > current_rank) {
                            continue;
                        }
                        reachable.insert(predecessor);
                        rank_frontier.push(predecessor);
                    }
                }
                return reachable;
            };
            const std::set<std::vector<std::size_t>> rank_reachable
                = rank_reachable_to(required_rank);
            std::uint64_t rank_monotone_traps = 0U;
            std::vector<std::size_t> first_rank_monotone_trap;
            for (const auto& [candidate, data] : order_landscape) {
                if (data.first < required_rank
                    && !rank_reachable.contains(candidate)) {
                    ++rank_monotone_traps;
                    if (first_rank_monotone_trap.empty()) {
                        first_rank_monotone_trap = candidate;
                    }
                }
            }
            std::cout << "rank_threshold=" << required_rank
                      << " rank_monotone_reachability_result="
                      << (rank_monotone_traps == 0U ? "PASS" : "FAIL")
                      << " rank_monotone_traps=" << rank_monotone_traps
                      << " first_rank_monotone_trap=[";
            for (std::size_t i = 0U;
                 i < first_rank_monotone_trap.size(); ++i) {
                if (i != 0U) {
                    std::cout << ',';
                }
                std::cout << first_rank_monotone_trap[i];
            }
            std::cout << "]\n";
            std::map<std::size_t, std::uint64_t> rank_histogram;
            for (const auto& [candidate, data] : order_landscape) {
                (void)candidate;
                ++rank_histogram[data.first];
            }
            std::cout << "rank_histogram=[";
            bool first_rank_entry = true;
            for (const auto& [rank, count] : rank_histogram) {
                if (!first_rank_entry) {
                    std::cout << ',';
                }
                std::cout << rank << ':' << count;
                first_rank_entry = false;
            }
            std::cout << "]\n";
            const std::set<std::vector<std::size_t>> maximum_reachable
                = rank_reachable_to(best.odd_rank);
            std::uint64_t maximum_monotone_traps = 0U;
            std::vector<std::size_t> first_maximum_monotone_trap;
            for (const auto& [candidate, data] : order_landscape) {
                if (data.first < best.odd_rank
                    && !maximum_reachable.contains(candidate)) {
                    ++maximum_monotone_traps;
                    if (first_maximum_monotone_trap.empty()) {
                        first_maximum_monotone_trap = candidate;
                    }
                }
            }
            std::cout << "maximum_rank=" << best.odd_rank
                      << " maximum_monotone_reachability_result="
                      << (maximum_monotone_traps == 0U ? "PASS" : "FAIL")
                      << " maximum_monotone_traps="
                      << maximum_monotone_traps
                      << " first_maximum_monotone_trap=[";
            for (std::size_t i = 0U;
                 i < first_maximum_monotone_trap.size(); ++i) {
                if (i != 0U) {
                    std::cout << ',';
                }
                std::cout << first_maximum_monotone_trap[i];
            }
            std::cout << "]\n";
            return kernel <= best.positive_sources ? 0 : 1;
        }
        if (argc != 6) {
            throw std::runtime_error(
                "usage: analyze_su2_plucker_sr_parity MODEL MAXIMUM_LABEL "
                "MAXIMUM_FACTORS MAXIMUM_TARGET MAXIMUM_TOTAL_DEGREE; "
                "MODEL is sr, gt, gt-zigzag, gt-best, gt-fiber, "
                "gt-fiber-best, gt-sd-best, gt-one-minus, or "
                "gt-deficit-fiber-best, "
                "gt-label-rounds, gt-scalar-sorted, or gt-scalar-balanced; "
                "gt-scalar-best, gt-xor-best; or "
                "gt-xor-rounds; or "
                "gt-xor-label-rounds; or "
                "gt-deficit-xor-label-rounds; or "
                "gt-deficit-xor-insert; or "
                "gt-average, gt-average-pairs, "
                "gt-average-weighted-pairs, gt-average-tree-pairs, or "
                "gt-average-transversal-pairs, gt-average-boolean-pairs, "
                "or gt-alpha-average; or "
                "gt-bk-components; or "
                "gt-bk-single-active; or "
                "gt-bk-single-loaded; or "
                "gt-monotone; or "
                "gt-case TARGET "
                "MINUS_MASK LABEL..."
            );
        }
        const std::string model = argv[1];
        if (model != "sr" && model != "gt" && model != "gt-zigzag"
            && model != "gt-best" && model != "gt-fiber"
            && model != "gt-fiber-best" && model != "gt-sd-best"
            && model != "gt-deficit-fiber-best"
            && model != "gt-one-minus" && model != "gt-label-rounds"
            && model != "gt-scalar-sorted"
            && model != "gt-scalar-balanced" && model != "gt-scalar-best"
            && model != "gt-xor-best" && model != "gt-xor-rounds"
            && model != "gt-xor-label-rounds"
            && model != "gt-deficit-xor-label-rounds"
            && model != "gt-deficit-xor-insert"
            && model != "gt-average"
            && model != "gt-average-pairs"
            && model != "gt-average-weighted-pairs"
            && model != "gt-average-tree-pairs"
            && model != "gt-average-transversal-pairs"
            && model != "gt-average-boolean-pairs"
            && model != "gt-alpha-average"
            && model != "gt-bk-components"
            && model != "gt-bk-single-active"
            && model != "gt-bk-single-loaded"
            && model != "gt-monotone") {
            throw std::runtime_error(
                "model must be sr, gt, gt-zigzag, gt-best, gt-fiber, or "
                "gt-fiber-best, gt-sd-best, gt-one-minus, or "
                "gt-deficit-fiber-best, "
                "gt-label-rounds, gt-scalar-sorted, or gt-scalar-balanced"
                ", gt-scalar-best, gt-xor-best, gt-xor-rounds, "
                "gt-xor-label-rounds, gt-deficit-xor-label-rounds"
                ", gt-deficit-xor-insert"
                ", gt-average"
                ", gt-average-pairs"
                ", gt-average-weighted-pairs"
                ", gt-average-tree-pairs"
                ", gt-average-transversal-pairs"
                ", gt-average-boolean-pairs"
                ", gt-alpha-average"
                ", gt-bk-components"
                ", gt-bk-single-active"
                ", gt-bk-single-loaded"
                ", gt-monotone"
            );
        }
        const int maximum_label = parse_nonnegative(argv[2], "maximum label");
        const int maximum_factors = parse_nonnegative(argv[3], "maximum factors");
        const int maximum_target = parse_nonnegative(argv[4], "maximum target");
        const int maximum_total_degree
            = parse_nonnegative(argv[5], "maximum total degree");
        if (maximum_label == 0 || maximum_factors == 0
            || maximum_factors >= static_cast<int>(
                std::numeric_limits<unsigned int>::digits
            )) {
            throw std::runtime_error("search bounds are out of range");
        }

        GraphCache cache;
        ToricCache toric_cache;
        std::uint64_t cases = 0U;
        for (int factors = 1; factors <= maximum_factors; ++factors) {
            std::vector<int> labels(static_cast<std::size_t>(factors), 1);
            bool more = true;
            while (more) {
                int label_sum = 0;
                for (int label : labels) {
                    label_sum += label;
                }
                for (int target = 0; target <= maximum_target; ++target) {
                    if ((model == "gt-scalar-sorted"
                         || model == "gt-scalar-balanced"
                         || model == "gt-scalar-best")
                        && target != 0) {
                        continue;
                    }
                    if (label_sum + target > maximum_total_degree
                        || ((label_sum + target) & 1) != 0) {
                        continue;
                    }
                    const unsigned int minus_limit = 1U << factors;
                    for (unsigned int minus_mask = 1U;
                         minus_mask < minus_limit; ++minus_mask) {
                        if (model == "gt-one-minus"
                            && popcount(minus_mask) != 1) {
                            continue;
                        }
                        if ((model == "gt-scalar-sorted"
                             || model == "gt-scalar-balanced"
                             || model == "gt-scalar-best")
                            && (popcount(minus_mask) & 1) != 0) {
                            continue;
                        }
                        if ((model == "gt-label-rounds"
                             || model == "gt-sd-best"
                             || model == "gt-deficit-fiber-best"
                             || model == "gt-one-minus"
                             || model == "gt-scalar-sorted"
                             || model == "gt-scalar-balanced"
                             || model == "gt-scalar-best"
                             || model == "gt-xor-best"
                             || model == "gt-xor-rounds"
                             || model == "gt-xor-label-rounds"
                             || model == "gt-deficit-xor-label-rounds"
                             || model == "gt-deficit-xor-insert"
                             || model == "gt-average"
                             || model == "gt-average-pairs"
                             || model == "gt-average-weighted-pairs"
                             || model == "gt-average-tree-pairs"
                             || model == "gt-average-transversal-pairs"
                             || model == "gt-average-boolean-pairs"
                             || model == "gt-alpha-average"
                             || model == "gt-bk-components"
                             || model == "gt-bk-single-active"
                             || model == "gt-bk-single-loaded"
                             || model == "gt-monotone")
                            && !support_disjoint(labels, minus_mask)) {
                            continue;
                        }
                        const Counts counts = model == "sr"
                            ? analyze_case(labels, target, minus_mask, cache)
                            : (model == "gt" || model == "gt-fiber"
                               || model == "gt-one-minus"
                               || model == "gt-scalar-sorted"
                                ? analyze_gt_case(
                                    labels, target, minus_mask, toric_cache
                                )
                                : (model == "gt-zigzag"
                                    ? analyze_gt_zigzag_case(
                                        labels, target, minus_mask, toric_cache
                                    )
                                    : (model == "gt-best"
                                        ? analyze_gt_best_case(
                                            labels, target, minus_mask,
                                            toric_cache
                                        )
                                        : (model == "gt-fiber-best"
                                           || model == "gt-sd-best"
                                            ? analyze_gt_fiber_best_case(
                                                labels, target, minus_mask,
                                                toric_cache
                                            )
                                            : (model
                                                == "gt-deficit-fiber-best"
                                                ? analyze_gt_deficit_fiber_best_case(
                                                    labels, target, minus_mask,
                                                    toric_cache
                                                )
                                            : (model == "gt-scalar-balanced"
                                                ? analyze_gt_scalar_balanced_case(
                                                    labels, minus_mask,
                                                    toric_cache
                                                )
                                                : (model == "gt-scalar-best"
                                                    ? analyze_gt_scalar_best_case(
                                                        labels, minus_mask,
                                                        toric_cache
                                                    )
                                                    : (model == "gt-xor-best"
                                                        ? analyze_gt_xor_best_case(
                                                            labels, target,
                                                            minus_mask,
                                                            toric_cache
                                                        )
                                                        : (model
                                                            == "gt-xor-rounds"
                                                            ? analyze_gt_xor_rounds_case(
                                                                labels, target,
                                                                minus_mask,
                                                                toric_cache
                                                            )
                                                        : (model
                                                            == "gt-deficit-xor-label-rounds"
                                                            ? analyze_gt_deficit_xor_label_rounds_case(
                                                                labels, target,
                                                                minus_mask,
                                                                toric_cache
                                                            )
                                                        : (model == "gt-monotone"
                                                            ? analyze_gt_monotone_case(
                                                                labels, target,
                                                                minus_mask,
                                                                toric_cache
                                                            )
                                                        : (model == "gt-average"
                                                            || model == "gt-average-pairs"
                                                            || model == "gt-average-weighted-pairs"
                                                            || model == "gt-average-tree-pairs"
                                                            || model == "gt-average-transversal-pairs"
                                                            || model == "gt-average-boolean-pairs"
                                                            ? (model == "gt-average-weighted-pairs"
                                                                ? analyze_gt_weighted_pair_average_case(
                                                                    labels, target,
                                                                    minus_mask,
                                                                    toric_cache
                                                                )
                                                                : (model == "gt-average-tree-pairs"
                                                                    ? analyze_gt_tree_pair_average_case(
                                                                        labels, target,
                                                                        minus_mask,
                                                                        toric_cache
                                                                    )
                                                                : (model == "gt-average-transversal-pairs"
                                                                    || model == "gt-average-boolean-pairs"
                                                                    ? analyze_gt_transversal_pair_average_case(
                                                                        labels, target,
                                                                        minus_mask,
                                                                        toric_cache,
                                                                        model == "gt-average-boolean-pairs"
                                                                    )
                                                                : analyze_gt_average_case(
                                                                    labels, target,
                                                                    minus_mask,
                                                                    toric_cache,
                                                                    model == "gt-average-pairs"
                                                                ))))
                                                        : (model == "gt-alpha-average"
                                                            ? analyze_gt_alpha_average_case(
                                                                labels, target,
                                                                minus_mask,
                                                                toric_cache
                                                            )
                                                        : (model == "gt-bk-components"
                                                            ? analyze_gt_bk_component_case(
                                                                labels, target,
                                                                minus_mask,
                                                                toric_cache
                                                            )
                                                        : (model == "gt-bk-single-active"
                                                            ? analyze_gt_bk_single_active_case(
                                                                labels, target,
                                                                minus_mask,
                                                                toric_cache
                                                            )
                                                        : (model == "gt-bk-single-loaded"
                                                            ? analyze_gt_bk_single_loaded_case(
                                                                labels, target,
                                                                minus_mask,
                                                                toric_cache
                                                            )
                                                        : (model
                                                            == "gt-deficit-xor-insert"
                                                            ? analyze_gt_xor_insert_case(
                                                                labels, target,
                                                                minus_mask,
                                                                toric_cache
                                                            )
                                                        : analyze_gt_label_rounds_case(
                                                            labels, target,
                                                            minus_mask,
                                                            toric_cache
                                                        ))))))))))))))))));
                        ++cases;
                        const std::uint64_t kernel = counts.odd_sources
                            - static_cast<std::uint64_t>(counts.odd_rank);
                        const bool failure = model == "gt-fiber"
                                || model == "gt-fiber-best"
                                || model == "gt-sd-best"
                                || model == "gt-deficit-fiber-best"
                                || model == "gt-label-rounds"
                                || model == "gt-one-minus"
                                || model == "gt-scalar-sorted"
                                || model == "gt-scalar-balanced"
                                || model == "gt-scalar-best"
                            ? !counts.fiber_parity
                            : model == "gt-xor-best"
                                || model == "gt-xor-rounds"
                                || model == "gt-xor-label-rounds"
                                || model == "gt-deficit-xor-label-rounds"
                                || model == "gt-deficit-xor-insert"
                                || model == "gt-average"
                                || model == "gt-average-pairs"
                                || model == "gt-average-weighted-pairs"
                                || model == "gt-average-tree-pairs"
                                || model == "gt-average-transversal-pairs"
                                || model == "gt-average-boolean-pairs"
                                || model == "gt-alpha-average"
                                || model == "gt-bk-components"
                                || model == "gt-bk-single-active"
                                || model == "gt-bk-single-loaded"
                                || model == "gt-monotone"
                                ? !counts.xor_payment
                            : kernel > counts.positive_sources;
                        if (failure) {
                            print_case(labels, target, minus_mask, counts);
                            std::cout << " result=FAIL\n";
                            return 1;
                        }
                    }
                }
                more = increment_labels(labels, maximum_label);
            }
            std::cout << "progress factors=" << factors
                      << " cases=" << cases
                      << " cached_degrees="
                      << (model == "sr" ? cache.size() : toric_cache.size())
                      << '\n' << std::flush;
        }
        std::cout << "SU2_PLUCKER_"
                  << (model == "sr" ? "SR"
                      : (model == "gt" ? "GT"
                         : (model == "gt-zigzag" ? "GT_ZIGZAG"
                            : (model == "gt-best" ? "GT_BEST"
                               : (model == "gt-fiber" ? "GT_FIBER"
                                  : (model == "gt-fiber-best"
                                        ? "GT_FIBER_BEST"
                                    : (model == "gt-sd-best"
                                        ? "GT_SD_BEST"
                                      : (model == "gt-deficit-fiber-best"
                                            ? "GT_DEFICIT_FIBER_BEST"
                                      : (model == "gt-one-minus"
                                            ? "GT_ONE_MINUS"
                                        : (model == "gt-scalar-sorted"
                                            ? "GT_SCALAR_SORTED"
                                          : (model == "gt-scalar-balanced"
                                                ? "GT_SCALAR_BALANCED"
                                            : (model == "gt-scalar-best"
                                                ? "GT_SCALAR_BEST"
                                              : (model == "gt-xor-best"
                                                    ? "GT_XOR_BEST"
                                                : (model == "gt-xor-rounds"
                                                    ? "GT_XOR_ROUNDS"
                                                : (model
                                                    == "gt-xor-label-rounds"
                                                    ? "GT_XOR_LABEL_ROUNDS"
                                                : (model
                                                    == "gt-deficit-xor-label-rounds"
                                                    ? "GT_DEFICIT_XOR_LABEL_ROUNDS"
                                                : (model == "gt-monotone"
                                                    ? "GT_RANK_MONOTONE"
                                                : (model == "gt-average"
                                                    ? "GT_AVERAGE"
                                                : (model == "gt-average-pairs"
                                                    ? "GT_AVERAGE_PAIRS"
                                                : (model == "gt-average-weighted-pairs"
                                                    ? "GT_AVERAGE_WEIGHTED_PAIRS"
                                                : (model == "gt-average-tree-pairs"
                                                    ? "GT_AVERAGE_TREE_PAIRS"
                                                : (model == "gt-average-transversal-pairs"
                                                    ? "GT_AVERAGE_TRANSVERSAL_PAIRS"
                                                : (model == "gt-average-boolean-pairs"
                                                    ? "GT_AVERAGE_BOOLEAN_PAIRS"
                                                : (model == "gt-alpha-average"
                                                    ? "GT_ALPHA_AVERAGE"
                                                : (model == "gt-bk-components"
                                                    ? "GT_BK_COMPONENTS"
                                                : (model == "gt-bk-single-active"
                                                    ? "GT_BK_SINGLE_ACTIVE"
                                                : (model == "gt-bk-single-loaded"
                                                    ? "GT_BK_SINGLE_LOADED"
                                                : (model
                                                    == "gt-deficit-xor-insert"
                                                    ? "GT_DEFICIT_XOR_INSERT"
                                                    : "GT_LABEL_ROUNDS"))))))))))))))))))))))))))))
                  << "_PARITY cases=" << cases
                  << " cached_degrees="
                  << (model == "sr" ? cache.size() : toric_cache.size())
                  << " result=PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
