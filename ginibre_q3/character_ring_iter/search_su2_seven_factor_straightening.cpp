#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

namespace {

constexpr std::size_t vertex_count = 7U;
constexpr std::int64_t modulus = 1'000'000'007LL;
using Graph = std::vector<std::vector<int>>;
using Key = std::vector<int>;
using Expansion = std::map<Key, std::int64_t>;

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

int popcount(unsigned int value) {
    int answer = 0;
    while (value != 0U) {
        answer += static_cast<int>(value & 1U);
        value >>= 1U;
    }
    return answer;
}

Key flatten(const Graph& graph) {
    Key key;
    key.reserve(vertex_count * (vertex_count - 1U) / 2U);
    for (std::size_t i = 0; i < vertex_count; ++i) {
        for (std::size_t j = i + 1U; j < vertex_count; ++j) {
            key.push_back(graph[i][j]);
        }
    }
    return key;
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
    for (std::size_t i = 0; i < vertex_count; ++i) {
        for (std::size_t j = i + 1U; j < vertex_count; ++j) {
            if (graph[i][j] != 0 && edges_cross(left, right, i, j)) {
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
    if (neighbor == vertex_count) {
        if (needed == 0) {
            const int saved = remaining[vertex];
            remaining[vertex] = 0;
            enumerate_graphs_rec(remaining, graph, output);
            remaining[vertex] = saved;
        }
        return;
    }
    int available = 0;
    for (std::size_t j = neighbor; j < vertex_count; ++j) {
        available += remaining[j];
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
    std::size_t vertex = 0;
    while (vertex < vertex_count && remaining[vertex] == 0) {
        ++vertex;
    }
    if (vertex == vertex_count) {
        output.push_back(graph);
        return;
    }
    distribute_vertex_degree(
        vertex, vertex + 1U, remaining[vertex], remaining, graph, output
    );
}

std::vector<Graph> noncrossing_graphs(
    const std::vector<int>& labels,
    unsigned int mask
) {
    std::vector<int> remaining(vertex_count, 0);
    for (std::size_t i = 0; i < vertex_count; ++i) {
        if (((mask >> i) & 1U) != 0U) {
            remaining[i] = labels[i];
        }
    }
    Graph graph(vertex_count, std::vector<int>(vertex_count, 0));
    std::vector<Graph> output;
    enumerate_graphs_rec(remaining, graph, output);
    return output;
}

std::pair<bool, std::array<std::size_t, 4>> first_crossing(
    const Graph& graph
) {
    for (std::size_t a = 0; a < vertex_count; ++a) {
        for (std::size_t b = a + 1U; b < vertex_count; ++b) {
            for (std::size_t c = b + 1U; c < vertex_count; ++c) {
                for (std::size_t d = c + 1U; d < vertex_count; ++d) {
                    if (graph[a][c] != 0 && graph[b][d] != 0) {
                        return {true, {a, b, c, d}};
                    }
                }
            }
        }
    }
    return {false, {0U, 0U, 0U, 0U}};
}

void add_expansion(Expansion& target, const Expansion& source) {
    for (const auto& [key, coefficient] : source) {
        std::int64_t& value = target[key];
        value += coefficient;
        if (value >= modulus) {
            value -= modulus;
        }
    }
}

const Expansion& straighten(
    const Graph& graph,
    std::map<Key, Expansion>& memo
) {
    const Key key = flatten(graph);
    const auto found = memo.find(key);
    if (found != memo.end()) {
        return found->second;
    }
    const auto [has_crossing, crossing] = first_crossing(graph);
    if (!has_crossing) {
        Expansion expansion;
        expansion.emplace(key, 1);
        return memo.emplace(key, std::move(expansion)).first->second;
    }
    const auto [a, b, c, d] = crossing;
    Graph separated = graph;
    --separated[a][c];
    --separated[b][d];
    ++separated[a][b];
    ++separated[c][d];
    Graph nested = graph;
    --nested[a][c];
    --nested[b][d];
    ++nested[a][d];
    ++nested[b][c];
    Expansion expansion;
    add_expansion(expansion, straighten(separated, memo));
    add_expansion(expansion, straighten(nested, memo));
    return memo.emplace(key, std::move(expansion)).first->second;
}

Graph graph_sum(const Graph& first, const Graph& second) {
    Graph result = first;
    for (std::size_t i = 0; i < vertex_count; ++i) {
        for (std::size_t j = i + 1U; j < vertex_count; ++j) {
            result[i][j] += second[i][j];
        }
    }
    return result;
}

bool compatible_splits(unsigned int first, unsigned int second) {
    constexpr unsigned int full_mask = (1U << vertex_count) - 1U;
    const unsigned int first_complement = full_mask ^ first;
    const unsigned int second_complement = full_mask ^ second;
    return (first & second) == 0U
        || (first & second_complement) == 0U
        || (first_complement & second) == 0U
        || (first_complement & second_complement) == 0U;
}

const std::vector<std::array<unsigned int, 4>>& tree_split_sets() {
    static const std::vector<std::array<unsigned int, 4>> trees = [] {
        constexpr unsigned int full_mask = (1U << vertex_count) - 1U;
        std::vector<unsigned int> splits;
        for (unsigned int mask = 1U; mask < full_mask; ++mask) {
            if ((mask & 1U) == 0U || popcount(mask) < 2
                || popcount(full_mask ^ mask) < 2) {
                continue;
            }
            splits.push_back(mask);
        }
        std::vector<std::array<unsigned int, 4>> result;
        for (std::size_t a = 0; a < splits.size(); ++a) {
            for (std::size_t b = a + 1U; b < splits.size(); ++b) {
                if (!compatible_splits(splits[a], splits[b])) {
                    continue;
                }
                for (std::size_t c = b + 1U; c < splits.size(); ++c) {
                    if (!compatible_splits(splits[a], splits[c])
                        || !compatible_splits(splits[b], splits[c])) {
                        continue;
                    }
                    for (std::size_t d = c + 1U; d < splits.size(); ++d) {
                        if (compatible_splits(splits[a], splits[d])
                            && compatible_splits(splits[b], splits[d])
                            && compatible_splits(splits[c], splits[d])) {
                            result.push_back({
                                splits[a], splits[b], splits[c], splits[d]
                            });
                        }
                    }
                }
            }
        }
        return result;
    }();
    return trees;
}

int cut_weight(const Graph& graph, unsigned int split) {
    int answer = 0;
    for (std::size_t i = 0; i < vertex_count; ++i) {
        for (std::size_t j = i + 1U; j < vertex_count; ++j) {
            if (((split >> i) & 1U) != ((split >> j) & 1U)) {
                answer += graph[i][j];
            }
        }
    }
    return answer;
}

std::int64_t tree_valuation_collisions(
    const std::vector<Graph>& sources,
    std::int64_t sufficient_bound
) {
    std::int64_t best = static_cast<std::int64_t>(sources.size());
    for (const auto& tree : tree_split_sets()) {
        std::set<std::array<int, 4>> values;
        for (const Graph& source : sources) {
            values.insert({
                cut_weight(source, tree[0]),
                cut_weight(source, tree[1]),
                cut_weight(source, tree[2]),
                cut_weight(source, tree[3])
            });
        }
        const std::int64_t collisions =
            static_cast<std::int64_t>(sources.size())
            - static_cast<std::int64_t>(values.size());
        best = std::min(best, collisions);
        if (best <= sufficient_bound) {
            break;
        }
    }
    return best;
}

std::int64_t modular_power(std::int64_t base, std::int64_t exponent) {
    std::int64_t answer = 1;
    while (exponent != 0) {
        if ((exponent & 1LL) != 0) {
            answer = (answer * base) % modulus;
        }
        base = (base * base) % modulus;
        exponent >>= 1LL;
    }
    return answer;
}

std::size_t modular_rank(std::vector<std::vector<std::int64_t>> matrix) {
    if (matrix.empty()) {
        return 0U;
    }
    const std::size_t rows = matrix.size();
    const std::size_t columns = matrix.front().size();
    std::size_t rank = 0;
    for (std::size_t column = 0; column < columns && rank < rows; ++column) {
        std::size_t pivot = rank;
        while (pivot < rows && matrix[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == rows) {
            continue;
        }
        std::swap(matrix[rank], matrix[pivot]);
        const std::int64_t inverse = modular_power(
            matrix[rank][column], modulus - 2
        );
        for (std::size_t j = column; j < columns; ++j) {
            matrix[rank][j] = (matrix[rank][j] * inverse) % modulus;
        }
        for (std::size_t i = 0; i < rows; ++i) {
            if (i == rank || matrix[i][column] == 0) {
                continue;
            }
            const std::int64_t factor = matrix[i][column];
            for (std::size_t j = column; j < columns; ++j) {
                matrix[i][j] = (
                    matrix[i][j] - factor * matrix[rank][j]
                ) % modulus;
                if (matrix[i][j] < 0) {
                    matrix[i][j] += modulus;
                }
            }
        }
        ++rank;
    }
    return rank;
}

std::vector<std::vector<std::int64_t>> modular_left_kernel(
    const std::vector<std::vector<std::int64_t>>& matrix
) {
    const std::size_t source_count = matrix.size();
    if (source_count == 0U) {
        return {};
    }
    const std::size_t target_count = matrix.front().size();
    std::vector<std::vector<std::int64_t>> transpose(
        target_count, std::vector<std::int64_t>(source_count, 0)
    );
    for (std::size_t source = 0; source < source_count; ++source) {
        for (std::size_t target = 0; target < target_count; ++target) {
            transpose[target][source] = matrix[source][target];
        }
    }
    std::vector<std::size_t> pivots;
    std::size_t row = 0U;
    for (std::size_t column = 0U;
         column < source_count && row < target_count; ++column) {
        std::size_t pivot = row;
        while (pivot < target_count && transpose[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == target_count) {
            continue;
        }
        std::swap(transpose[row], transpose[pivot]);
        const std::int64_t inverse = modular_power(
            transpose[row][column], modulus - 2
        );
        for (std::size_t j = column; j < source_count; ++j) {
            transpose[row][j] = (transpose[row][j] * inverse) % modulus;
        }
        for (std::size_t i = 0U; i < target_count; ++i) {
            if (i == row || transpose[i][column] == 0) {
                continue;
            }
            const std::int64_t factor = transpose[i][column];
            for (std::size_t j = column; j < source_count; ++j) {
                transpose[i][j] = (
                    transpose[i][j] - factor * transpose[row][j]
                ) % modulus;
                if (transpose[i][j] < 0) {
                    transpose[i][j] += modulus;
                }
            }
        }
        pivots.push_back(column);
        ++row;
    }
    std::vector<bool> is_pivot(source_count, false);
    for (std::size_t pivot : pivots) {
        is_pivot[pivot] = true;
    }
    std::vector<std::vector<std::int64_t>> basis;
    for (std::size_t free_column = 0U; free_column < source_count;
         ++free_column) {
        if (is_pivot[free_column]) {
            continue;
        }
        std::vector<std::int64_t> vector(source_count, 0);
        vector[free_column] = 1;
        for (std::size_t pivot_row = 0U; pivot_row < pivots.size();
             ++pivot_row) {
            const std::int64_t coefficient =
                transpose[pivot_row][free_column];
            vector[pivots[pivot_row]] =
                coefficient == 0 ? 0 : modulus - coefficient;
        }
        basis.push_back(std::move(vector));
    }
    return basis;
}

std::int64_t modular_binomial(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    k = std::min(k, n - k);
    std::int64_t answer = 1;
    for (int i = 1; i <= k; ++i) {
        answer = (answer * static_cast<std::int64_t>(n - k + i)) % modulus;
        answer = (answer * modular_power(i, modulus - 2)) % modulus;
    }
    return answer;
}

Expansion expand_graph_coordinates(const Graph& graph) {
    Expansion current;
    current.emplace(Key(vertex_count, 0), 1);
    for (std::size_t i = 0U; i < vertex_count; ++i) {
        for (std::size_t j = i + 1U; j < vertex_count; ++j) {
            const int multiplicity = graph[i][j];
            if (multiplicity == 0) {
                continue;
            }
            Expansion next;
            for (const auto& [key, coefficient] : current) {
                for (int chosen = 0; chosen <= multiplicity; ++chosen) {
                    Key output = key;
                    output[i] += chosen;
                    output[j] += multiplicity - chosen;
                    std::int64_t term = modular_binomial(
                        multiplicity, chosen
                    );
                    if ((chosen & 1) != 0) {
                        term = modulus - term;
                    }
                    std::int64_t& value = next[output];
                    value = (
                        value + (coefficient * term) % modulus
                    ) % modulus;
                }
            }
            current = std::move(next);
        }
    }
    return current;
}

Expansion singlet_projection(
    const Graph& graph,
    const std::vector<int>& labels,
    std::size_t first,
    std::size_t second,
    int pair_code
) {
    const int degree = labels[first];
    Expansion answer;
    for (const auto& [key, coefficient] : expand_graph_coordinates(graph)) {
        if (key[first] + key[second] != degree) {
            continue;
        }
        const int first_y_degree = key[first];
        std::int64_t weight = modular_power(
            modular_binomial(degree, first_y_degree), modulus - 2
        );
        if ((first_y_degree & 1) != 0) {
            weight = modulus - weight;
        }
        Key output;
        output.reserve(vertex_count - 1U);
        output.push_back(pair_code);
        for (std::size_t i = 0U; i < vertex_count; ++i) {
            if (i != first && i != second) {
                output.push_back(key[i]);
            }
        }
        std::int64_t& value = answer[output];
        value = (value + coefficient * weight) % modulus;
    }
    return answer;
}

std::int64_t residue_rank_on_kernel(
    const std::vector<int>& labels,
    std::size_t minus_count,
    const std::vector<Graph>& sources,
    const std::vector<unsigned int>& source_masks,
    const std::vector<std::vector<std::int64_t>>& source_matrix
) {
    const auto kernel = modular_left_kernel(source_matrix);
    if (kernel.empty()) {
        return 0;
    }
    std::vector<Expansion> residue_rows(sources.size());
    for (std::size_t first = 0U; first < vertex_count; ++first) {
        for (std::size_t second = first + 1U; second < vertex_count; ++second) {
            const bool same_sign = (first < minus_count) == (second < minus_count);
            if (!same_sign || labels[first] != labels[second]) {
                continue;
            }
            const int pair_code = static_cast<int>(
                first * vertex_count + second
            );
            for (std::size_t source = 0U; source < sources.size(); ++source) {
                const bool contains_first =
                    ((source_masks[source] >> first) & 1U) != 0U;
                const bool contains_second =
                    ((source_masks[source] >> second) & 1U) != 0U;
                if (!contains_first || contains_second) {
                    continue;
                }
                const Expansion projection = singlet_projection(
                    sources[source], labels, first, second, pair_code
                );
                add_expansion(residue_rows[source], projection);
            }
        }
    }
    std::map<Key, std::size_t> columns;
    for (const Expansion& row : residue_rows) {
        for (const auto& [key, coefficient] : row) {
            (void)coefficient;
            if (columns.find(key) == columns.end()) {
                columns.emplace(key, columns.size());
            }
        }
    }
    std::vector<std::vector<std::int64_t>> image(
        kernel.size(), std::vector<std::int64_t>(columns.size(), 0)
    );
    for (std::size_t relation = 0U; relation < kernel.size(); ++relation) {
        for (std::size_t source = 0U; source < sources.size(); ++source) {
            if (kernel[relation][source] == 0) {
                continue;
            }
            for (const auto& [key, coefficient] : residue_rows[source]) {
                std::int64_t& value = image[relation][columns.at(key)];
                value = (
                    value + kernel[relation][source] * coefficient
                ) % modulus;
            }
        }
    }
    return static_cast<std::int64_t>(modular_rank(std::move(image)));
}

void print_graph(const Graph& graph) {
    bool first = true;
    std::cout << '{';
    for (std::size_t i = 0U; i < vertex_count; ++i) {
        for (std::size_t j = i + 1U; j < vertex_count; ++j) {
            if (graph[i][j] == 0) {
                continue;
            }
            if (!first) {
                std::cout << ',';
            }
            first = false;
            std::cout << i << '-' << j << ':' << graph[i][j];
        }
    }
    std::cout << '}';
}

std::int64_t leading_collisions_for_order(
    const std::vector<Graph>& sources,
    const std::array<std::size_t, vertex_count>& order
) {
    std::map<Key, Expansion> memo;
    std::set<Key> smallest_terms;
    std::set<Key> largest_terms;
    for (const Graph& source : sources) {
        Graph permuted(vertex_count, std::vector<int>(vertex_count, 0));
        for (std::size_t i = 0U; i < vertex_count; ++i) {
            for (std::size_t j = i + 1U; j < vertex_count; ++j) {
                const std::size_t old_i = order[i];
                const std::size_t old_j = order[j];
                permuted[i][j] = source[std::min(old_i, old_j)]
                    [std::max(old_i, old_j)];
            }
        }
        const Expansion& expansion = straighten(permuted, memo);
        if (!expansion.empty()) {
            smallest_terms.insert(expansion.begin()->first);
            largest_terms.insert(expansion.rbegin()->first);
        }
    }
    return static_cast<std::int64_t>(sources.size())
        - static_cast<std::int64_t>(
            std::max(smallest_terms.size(), largest_terms.size())
        );
}

void print_order(const std::array<std::size_t, vertex_count>& order) {
    std::cout << '[';
    for (std::size_t i = 0U; i < vertex_count; ++i) {
        if (i != 0U) {
            std::cout << ',';
        }
        std::cout << order[i];
    }
    std::cout << ']';
}

std::pair<std::int64_t, std::array<std::size_t, vertex_count>>
best_circular_order_collisions(const std::vector<Graph>& sources) {
    std::array<std::size_t, vertex_count> order{};
    for (std::size_t i = 0U; i < vertex_count; ++i) {
        order[i] = i;
    }
    std::int64_t best = static_cast<std::int64_t>(sources.size());
    auto best_order = order;
    do {
        const std::int64_t collisions = leading_collisions_for_order(
            sources, order
        );
        if (collisions < best) {
            best = collisions;
            best_order = order;
            if (best == 0) {
                break;
            }
        }
    } while (std::next_permutation(order.begin() + 1, order.end()));
    return {best, best_order};
}

struct CaseResult {
    std::int64_t sources = 0;
    std::int64_t invariant = 0;
    std::int64_t rank = 0;
    std::int64_t kernel = 0;
    std::int64_t pair_term = 0;
    std::int64_t minimum_leading_collisions = 0;
    std::int64_t tree_collisions = 0;
    std::int64_t residue_rank = 0;
};

CaseResult analyze_case(
    const std::vector<int>& labels,
    std::size_t minus_count,
    bool verbose = false,
    bool evaluate_tree = true,
    bool evaluate_residue = false
) {
    constexpr unsigned int full_mask = (1U << vertex_count) - 1U;
    std::vector<Graph> sources;
    std::vector<unsigned int> source_masks;
    for (unsigned int mask = 0; mask <= full_mask; ++mask) {
        if (popcount(mask) != 3) {
            continue;
        }
        int parity = 0;
        for (std::size_t i = 0; i < minus_count; ++i) {
            parity ^= static_cast<int>((mask >> i) & 1U);
        }
        if (parity == 0) {
            continue;
        }
        const auto triple = noncrossing_graphs(labels, mask);
        const auto complement = noncrossing_graphs(labels, full_mask ^ mask);
        for (const Graph& first : triple) {
            for (const Graph& second : complement) {
                sources.push_back(graph_sum(first, second));
                source_masks.push_back(mask);
            }
        }
    }

    std::int64_t pair_term = 0;
    for (std::size_t i = 0; i < vertex_count; ++i) {
        for (std::size_t j = i + 1U; j < vertex_count; ++j) {
            if (labels[i] != labels[j]) {
                continue;
            }
            unsigned int rest_mask = full_mask;
            rest_mask ^= 1U << i;
            rest_mask ^= 1U << j;
            pair_term += static_cast<std::int64_t>(
                noncrossing_graphs(labels, rest_mask).size()
            );
        }
    }

    std::map<Key, Expansion> memo;
    std::vector<Expansion> expansions;
    std::map<Key, std::size_t> columns;
    std::map<Key, std::size_t> smallest_terms;
    std::map<Key, std::size_t> largest_terms;
    for (const Graph& source : sources) {
        expansions.push_back(straighten(source, memo));
        if (!expansions.back().empty()) {
            smallest_terms.emplace(expansions.back().begin()->first, 0U);
            largest_terms.emplace(expansions.back().rbegin()->first, 0U);
        }
        for (const auto& [key, coefficient] : expansions.back()) {
            (void)coefficient;
            if (columns.find(key) == columns.end()) {
                columns.emplace(key, columns.size());
            }
        }
    }
    std::vector<std::vector<std::int64_t>> matrix(
        expansions.size(), std::vector<std::int64_t>(columns.size(), 0)
    );
    for (std::size_t i = 0; i < expansions.size(); ++i) {
        for (const auto& [key, coefficient] : expansions[i]) {
            matrix[i][columns.at(key)] = coefficient;
        }
    }
    const std::int64_t rank = static_cast<std::int64_t>(modular_rank(matrix));
    const std::int64_t invariant = invariant_multiplicity(labels);
    const std::int64_t valuation_collisions = evaluate_tree
        ? tree_valuation_collisions(sources, pair_term)
        : -1;
    const std::int64_t residue_rank = (verbose || evaluate_residue)
        ? residue_rank_on_kernel(
            labels, minus_count, sources, source_masks, matrix
        )
        : -1;
    if (verbose) {
        const auto kernel_basis = modular_left_kernel(matrix);
        std::cout << "sources=" << sources.size() << " rank=" << rank
                  << " kernel=" << kernel_basis.size()
                  << " invariant=" << invariant
                  << " pair_term=" << pair_term << '\n';
        std::cout << "residue_rank_on_kernel=" << residue_rank << '\n';
        std::array<std::size_t, vertex_count> label_order{};
        for (std::size_t i = 0U; i < vertex_count; ++i) {
            label_order[i] = i;
        }
        std::stable_sort(
            label_order.begin(), label_order.end(),
            [&labels](std::size_t first, std::size_t second) {
                return labels[first] < labels[second];
            }
        );
        std::cout << "label_order=";
        print_order(label_order);
        std::cout << " label_order_collisions="
                  << leading_collisions_for_order(sources, label_order) << '\n';
        std::reverse(label_order.begin(), label_order.end());
        std::cout << "reverse_label_order=";
        print_order(label_order);
        std::cout << " reverse_label_order_collisions="
                  << leading_collisions_for_order(sources, label_order) << '\n';
        const auto [best_circular_collisions, best_circular_order] =
            best_circular_order_collisions(sources);
        std::cout << "best_circular_order=";
        print_order(best_circular_order);
        std::cout << " best_circular_order_collisions="
                  << best_circular_collisions << '\n';
        for (std::size_t i = 0U; i < sources.size(); ++i) {
            std::cout << "source " << i << " mask=" << source_masks[i]
                      << " graph=";
            print_graph(sources[i]);
            std::cout << '\n';
        }
        for (std::size_t relation = 0U; relation < kernel_basis.size();
             ++relation) {
            std::cout << "relation " << relation << ':';
            for (std::size_t source = 0U; source < sources.size(); ++source) {
                std::int64_t coefficient = kernel_basis[relation][source];
                if (coefficient == 0) {
                    continue;
                }
                if (coefficient > modulus / 2) {
                    coefficient -= modulus;
                }
                std::cout << ' ' << source << ':' << coefficient;
            }
            std::cout << '\n';
        }
    }
    return CaseResult{
        static_cast<std::int64_t>(sources.size()), invariant, rank,
        static_cast<std::int64_t>(sources.size()) - rank, pair_term,
        static_cast<std::int64_t>(sources.size())
            - static_cast<std::int64_t>(
                std::max(smallest_terms.size(), largest_terms.size())
            ),
        valuation_collisions, residue_rank
    };
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
    for (int label : minus) {
        if (std::binary_search(plus.begin(), plus.end(), label)) {
            return false;
        }
    }
    return true;
}

bool has_repeated_fundamental(
    const std::vector<int>& minus,
    const std::vector<int>& plus
) {
    return static_cast<int>(std::count(minus.begin(), minus.end(), 1)) >= 2
        || static_cast<int>(std::count(plus.begin(), plus.end(), 1)) >= 2;
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

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--two-odd-edge") {
            const long parsed = std::strtol(argv[2], nullptr, 10);
            if (parsed < 3 || parsed > std::numeric_limits<int>::max()) {
                throw std::runtime_error("invalid two-odd-edge maximum label");
            }
            const int maximum_label = static_cast<int>(parsed);
            long long cases = 0;
            std::int64_t maximum_kernel = 0;
            for (int q = 3; q <= maximum_label; q += 2) {
                for (int r = q; r <= maximum_label; r += 2) {
                    for (int a = 3; a <= maximum_label; a += 2) {
                        if (a == q || a == r) {
                            continue;
                        }
                        for (int first = 2; first <= maximum_label; first += 2) {
                            for (int second = first;
                                 second <= maximum_label;
                                 second += 2) {
                                for (int third = second;
                                     third <= maximum_label;
                                     third += 2) {
                                    const std::vector<int> labels{
                                        q, r, a, 1, first, second, third
                                    };
                                    const CaseResult result = analyze_case(
                                        labels, 2U, false, false, true
                                    );
                                    ++cases;
                                    maximum_kernel = std::max(
                                        maximum_kernel, result.kernel
                                    );
                                    if (result.residue_rank != result.kernel) {
                                        std::cout
                                            << "TWO_ODD_EDGE_RESIDUE_FAIL labels=";
                                        print_vector(labels);
                                        std::cout << " sources=" << result.sources
                                                  << " rank=" << result.rank
                                                  << " kernel=" << result.kernel
                                                  << " residue_rank="
                                                  << result.residue_rank
                                                  << " pair_term="
                                                  << result.pair_term << '\n';
                                        return EXIT_SUCCESS;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            std::cout << "SU2_SEVEN_FACTOR_TWO_ODD_EDGE_RESIDUE PASS cases="
                      << cases << " maximum_kernel=" << maximum_kernel
                      << " maximum_label=" << maximum_label << '\n';
            return EXIT_SUCCESS;
        }
        if (argc == 10 && std::string(argv[1]) == "--case") {
            const long parsed_minus_count = std::strtol(argv[2], nullptr, 10);
            if (parsed_minus_count < 0
                || parsed_minus_count > static_cast<long>(vertex_count)) {
                throw std::runtime_error("invalid case minus count");
            }
            std::vector<int> labels;
            for (int i = 3; i < 10; ++i) {
                const long label = std::strtol(argv[i], nullptr, 10);
                if (label <= 0 || label > std::numeric_limits<int>::max()) {
                    throw std::runtime_error("invalid case label");
                }
                labels.push_back(static_cast<int>(label));
            }
            (void)analyze_case(
                labels, static_cast<std::size_t>(parsed_minus_count), true
            );
            return 0;
        }
        bool evaluate_tree = true;
        const char* maximum_label_argument = nullptr;
        if (argc == 2) {
            maximum_label_argument = argv[1];
        } else if (argc == 3 && std::string(argv[1]) == "--fast") {
            evaluate_tree = false;
            maximum_label_argument = argv[2];
        } else {
            throw std::runtime_error(
                "usage: search_su2_seven_factor_straightening MAXIMUM_LABEL"
                " | --fast MAXIMUM_LABEL"
                " | --case MINUS_COUNT LABEL0 ... LABEL6"
                " | --two-odd-edge MAXIMUM_LABEL"
            );
        }
        const long parsed = std::strtol(maximum_label_argument, nullptr, 10);
        if (parsed <= 0 || parsed > std::numeric_limits<int>::max()) {
            throw std::runtime_error("invalid maximum label");
        }
        const int maximum_label = static_cast<int>(parsed);
        struct SignedCase {
            std::vector<int> minus;
            std::vector<int> plus;
        };
        std::vector<SignedCase> cases;
        for (int minus_count = 2; minus_count <= 6; minus_count += 2) {
            const auto minus_lists = combinations(maximum_label, minus_count);
            const auto plus_lists = combinations(maximum_label, 7 - minus_count);
            for (const auto& minus : minus_lists) {
                for (const auto& plus : plus_lists) {
                    if (supports_are_disjoint(minus, plus)) {
                        cases.push_back({minus, plus});
                    }
                }
            }
        }

        bool core_failed = false;
        bool tree_criterion_failed = false;
        std::int64_t worst_relation_margin = std::numeric_limits<std::int64_t>::max();
        std::int64_t worst_count_margin = std::numeric_limits<std::int64_t>::max();
        std::int64_t active_relation_margin = std::numeric_limits<std::int64_t>::max();
        std::int64_t excess_relation_margin = std::numeric_limits<std::int64_t>::max();
        std::int64_t worst_tree_margin = std::numeric_limits<std::int64_t>::max();
        std::int64_t maximum_kernel_without_repeated_fundamental = -1;
        SignedCase worst_relation_case;
        SignedCase worst_count_case;
        SignedCase active_relation_case;
        SignedCase excess_relation_case;
        SignedCase worst_tree_case;
        CaseResult worst_relation_result;
        CaseResult worst_count_result;
        CaseResult active_relation_result;
        CaseResult excess_relation_result;
        CaseResult worst_tree_result;
        SignedCase maximum_no_fundamental_case;
        CaseResult maximum_no_fundamental_result;

#pragma omp parallel
        {
            bool local_core_failed = false;
            bool local_tree_criterion_failed = false;
            std::int64_t local_relation_margin =
                std::numeric_limits<std::int64_t>::max();
            std::int64_t local_count_margin =
                std::numeric_limits<std::int64_t>::max();
            std::int64_t local_active_margin =
                std::numeric_limits<std::int64_t>::max();
            std::int64_t local_excess_margin =
                std::numeric_limits<std::int64_t>::max();
            std::int64_t local_tree_margin =
                std::numeric_limits<std::int64_t>::max();
            std::int64_t local_maximum_no_fundamental_kernel = -1;
            SignedCase local_relation_case;
            SignedCase local_count_case;
            SignedCase local_active_case;
            SignedCase local_excess_case;
            SignedCase local_tree_case;
            CaseResult local_relation_result;
            CaseResult local_count_result;
            CaseResult local_active_result;
            CaseResult local_excess_result;
            CaseResult local_tree_result;
            SignedCase local_maximum_no_fundamental_case;
            CaseResult local_maximum_no_fundamental_result;

#pragma omp for schedule(dynamic)
            for (std::size_t index = 0; index < cases.size(); ++index) {
                std::vector<int> labels = cases[index].minus;
                labels.insert(
                    labels.end(), cases[index].plus.begin(), cases[index].plus.end()
                );
                const CaseResult result = analyze_case(
                    labels, cases[index].minus.size(), false, evaluate_tree
                );
                const std::int64_t relation_margin =
                    result.pair_term - result.kernel;
                const std::int64_t count_margin =
                    result.invariant + result.pair_term - result.sources;
                const std::int64_t tree_margin =
                    result.pair_term - result.tree_collisions;
                if (relation_margin < local_relation_margin) {
                    local_relation_margin = relation_margin;
                    local_relation_case = cases[index];
                    local_relation_result = result;
                }
                if (count_margin < local_count_margin) {
                    local_count_margin = count_margin;
                    local_count_case = cases[index];
                    local_count_result = result;
                }
                if (result.sources > 0 && relation_margin < local_active_margin) {
                    local_active_margin = relation_margin;
                    local_active_case = cases[index];
                    local_active_result = result;
                }
                if (result.sources > result.invariant
                    && relation_margin < local_excess_margin) {
                    local_excess_margin = relation_margin;
                    local_excess_case = cases[index];
                    local_excess_result = result;
                }
                if (tree_margin < local_tree_margin) {
                    local_tree_margin = tree_margin;
                    local_tree_case = cases[index];
                    local_tree_result = result;
                }
                if (!has_repeated_fundamental(
                        cases[index].minus, cases[index].plus
                    )
                    && result.kernel > local_maximum_no_fundamental_kernel) {
                    local_maximum_no_fundamental_kernel = result.kernel;
                    local_maximum_no_fundamental_case = cases[index];
                    local_maximum_no_fundamental_result = result;
                }
                if (relation_margin < 0 || count_margin < 0
                    || result.rank > result.invariant) {
                    local_core_failed = true;
                }
                if (result.tree_collisions >= 0
                    && result.tree_collisions > result.pair_term) {
                    local_tree_criterion_failed = true;
                }
            }

#pragma omp critical
            {
                core_failed = core_failed || local_core_failed;
                tree_criterion_failed =
                    tree_criterion_failed || local_tree_criterion_failed;
                if (local_relation_margin < worst_relation_margin) {
                    worst_relation_margin = local_relation_margin;
                    worst_relation_case = std::move(local_relation_case);
                    worst_relation_result = local_relation_result;
                }
                if (local_count_margin < worst_count_margin) {
                    worst_count_margin = local_count_margin;
                    worst_count_case = std::move(local_count_case);
                    worst_count_result = local_count_result;
                }
                if (local_active_margin < active_relation_margin) {
                    active_relation_margin = local_active_margin;
                    active_relation_case = std::move(local_active_case);
                    active_relation_result = local_active_result;
                }
                if (local_excess_margin < excess_relation_margin) {
                    excess_relation_margin = local_excess_margin;
                    excess_relation_case = std::move(local_excess_case);
                    excess_relation_result = local_excess_result;
                }
                if (local_tree_margin < worst_tree_margin) {
                    worst_tree_margin = local_tree_margin;
                    worst_tree_case = std::move(local_tree_case);
                    worst_tree_result = local_tree_result;
                }
                if (local_maximum_no_fundamental_kernel
                    > maximum_kernel_without_repeated_fundamental) {
                    maximum_kernel_without_repeated_fundamental =
                        local_maximum_no_fundamental_kernel;
                    maximum_no_fundamental_case =
                        std::move(local_maximum_no_fundamental_case);
                    maximum_no_fundamental_result =
                        local_maximum_no_fundamental_result;
                }
            }
        }

        std::cout << "SU2_SEVEN_FACTOR_STRAIGHTENING cases=" << cases.size()
                  << " maximum_label=" << maximum_label << '\n';
        std::cout << "minimum_relation_margin=" << worst_relation_margin
                  << " minus=";
        print_vector(worst_relation_case.minus);
        std::cout << " plus=";
        print_vector(worst_relation_case.plus);
        std::cout << " sources=" << worst_relation_result.sources
                  << " rank=" << worst_relation_result.rank
                  << " kernel=" << worst_relation_result.kernel
                  << " pair_term=" << worst_relation_result.pair_term << '\n';
        std::cout << "minimum_count_margin=" << worst_count_margin
                  << " minus=";
        print_vector(worst_count_case.minus);
        std::cout << " plus=";
        print_vector(worst_count_case.plus);
        std::cout << " sources=" << worst_count_result.sources
                  << " invariant=" << worst_count_result.invariant
                  << " pair_term=" << worst_count_result.pair_term << '\n';
        if (active_relation_margin != std::numeric_limits<std::int64_t>::max()) {
            std::cout << "active_minimum_relation_margin="
                      << active_relation_margin << " minus=";
            print_vector(active_relation_case.minus);
            std::cout << " plus=";
            print_vector(active_relation_case.plus);
            std::cout << " sources=" << active_relation_result.sources
                      << " rank=" << active_relation_result.rank
                      << " kernel=" << active_relation_result.kernel
                      << " invariant=" << active_relation_result.invariant
                      << " pair_term=" << active_relation_result.pair_term
                      << " minimum_leading_collisions="
                      << active_relation_result.minimum_leading_collisions
                      << " tree_collisions="
                      << active_relation_result.tree_collisions
                      << " residue_rank="
                      << active_relation_result.residue_rank << '\n';
        }
        if (excess_relation_margin != std::numeric_limits<std::int64_t>::max()) {
            std::cout << "excess_minimum_relation_margin="
                      << excess_relation_margin << " minus=";
            print_vector(excess_relation_case.minus);
            std::cout << " plus=";
            print_vector(excess_relation_case.plus);
            std::cout << " sources=" << excess_relation_result.sources
                      << " rank=" << excess_relation_result.rank
                      << " kernel=" << excess_relation_result.kernel
                      << " invariant=" << excess_relation_result.invariant
                      << " pair_term=" << excess_relation_result.pair_term
                      << " minimum_leading_collisions="
                      << excess_relation_result.minimum_leading_collisions
                      << " tree_collisions="
                      << excess_relation_result.tree_collisions
                      << " residue_rank="
                      << excess_relation_result.residue_rank << '\n';
        }
        if (evaluate_tree) {
            std::cout << "minimum_tree_margin=" << worst_tree_margin
                      << " minus=";
            print_vector(worst_tree_case.minus);
            std::cout << " plus=";
            print_vector(worst_tree_case.plus);
            std::cout << " sources=" << worst_tree_result.sources
                      << " kernel=" << worst_tree_result.kernel
                      << " pair_term=" << worst_tree_result.pair_term
                      << " tree_collisions="
                      << worst_tree_result.tree_collisions << '\n';
        }
        std::cout << "maximum_kernel_without_repeated_fundamental="
                  << maximum_kernel_without_repeated_fundamental << " minus=";
        print_vector(maximum_no_fundamental_case.minus);
        std::cout << " plus=";
        print_vector(maximum_no_fundamental_case.plus);
        std::cout << " sources=" << maximum_no_fundamental_result.sources
                  << " rank=" << maximum_no_fundamental_result.rank
                  << " pair_term=" << maximum_no_fundamental_result.pair_term
                  << '\n';
        std::cout << "SU2_SEVEN_FACTOR_STRAIGHTENING_CORE "
                  << (core_failed ? "FAIL" : "PASS") << '\n';
        std::cout << "SU2_SEVEN_FACTOR_SINGLE_TREE_CRITERION "
                  << (evaluate_tree
                        ? (tree_criterion_failed ? "FALSE" : "PASS")
                        : "SKIPPED")
                  << '\n';
        return core_failed ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SEVEN_FACTOR_STRAIGHTENING FAILURE: "
                  << error.what() << '\n';
        return 1;
    }
}
