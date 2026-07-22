#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

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
};

struct ToricFiber {
    std::array<std::uint64_t, 2U> parity_counts{};
    std::vector<unsigned int> odd_cuts;
    std::vector<unsigned int> nonempty_even_cuts;
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
    const std::vector<std::size_t>& order
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
        (void)key;
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
    std::map<
        std::vector<std::size_t>,
        std::pair<std::size_t, bool>
    > landscape;
    do {
        const Counts counts = analyze_gt_case(
            labels, target, minus_mask, cache, order
        );
        landscape.emplace(
            order,
            std::pair<std::size_t, bool>{
                counts.odd_rank, counts.xor_payment
            }
        );
    } while (std::next_permutation(order.begin(), order.end()));
    std::set<std::vector<std::size_t>> reachable;
    std::queue<std::vector<std::size_t>> frontier;
    for (const auto& [candidate, data] : landscape) {
        if (data.second) {
            reachable.insert(candidate);
            frontier.push(candidate);
        }
    }
    while (!frontier.empty()) {
        const std::vector<std::size_t> current = frontier.front();
        frontier.pop();
        const std::size_t current_rank = landscape.at(current).first;
        for (std::size_t i = 0U; i + 1U < current.size(); ++i) {
            std::vector<std::size_t> predecessor = current;
            std::swap(predecessor[i], predecessor[i + 1U]);
            if (reachable.contains(predecessor)
                || landscape.at(predecessor).first > current_rank) {
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
        first.xor_capacity = 0U;
    }
    return first;
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
            std::cout << "monotone_reachability_result="
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
                std::cout << "plateau_escape=";
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
            && model != "gt-monotone") {
            throw std::runtime_error(
                "model must be sr, gt, gt-zigzag, gt-best, gt-fiber, or "
                "gt-fiber-best, gt-sd-best, gt-one-minus, or "
                "gt-deficit-fiber-best, "
                "gt-label-rounds, gt-scalar-sorted, or gt-scalar-balanced"
                ", gt-scalar-best, gt-xor-best, gt-xor-rounds, "
                "gt-xor-label-rounds, gt-deficit-xor-label-rounds"
                ", gt-deficit-xor-insert"
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
                                                        )))))))))))));
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
                                                    ? "GT_MONOTONE"
                                                : (model
                                                    == "gt-deficit-xor-insert"
                                                    ? "GT_DEFICIT_XOR_INSERT"
                                                    : "GT_LABEL_ROUNDS"))))))))))))))))))
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
