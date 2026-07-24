#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using Key = std::vector<int>;
using Graph = std::vector<std::vector<int>>;
using Expansion = std::map<Key, Integer>;
using ChannelKey = std::tuple<unsigned int, Key, Key>;

struct Edge {
    std::size_t first = 0U;
    std::size_t second = 0U;
};

int popcount(unsigned int value) {
    int answer = 0;
    while (value != 0U) {
        answer += static_cast<int>(value & 1U);
        value >>= 1U;
    }
    return answer;
}

bool valid_top(const std::vector<int>& degrees, const Key& top) {
    int degree_total = 0;
    int top_total = 0;
    int top_before = 0;
    int bottom_through = 0;
    for (std::size_t vertex = 0U; vertex < degrees.size(); ++vertex) {
        if (top[vertex] < 0 || top[vertex] > degrees[vertex]) {
            return false;
        }
        degree_total += degrees[vertex];
        top_total += top[vertex];
        bottom_through += degrees[vertex] - top[vertex];
        if (bottom_through > top_before) {
            return false;
        }
        top_before += top[vertex];
    }
    return 2 * top_total == degree_total;
}

void enumerate_tops(
    const std::vector<int>& degrees,
    std::size_t vertex,
    Key& top,
    std::vector<Key>& output
) {
    if (vertex == degrees.size()) {
        if (valid_top(degrees, top)) {
            output.push_back(top);
        }
        return;
    }
    for (int value = 0; value <= degrees[vertex]; ++value) {
        top[vertex] = value;
        enumerate_tops(degrees, vertex + 1U, top, output);
    }
}

std::vector<Key> tops(const std::vector<int>& degrees) {
    Key top(degrees.size(), 0);
    std::vector<Key> output;
    enumerate_tops(degrees, 0U, top, output);
    return output;
}

Key flatten(const Graph& graph) {
    Key key;
    key.reserve(graph.size() * (graph.size() - 1U) / 2U);
    for (std::size_t left = 0U; left < graph.size(); ++left) {
        for (std::size_t right = left + 1U;
             right < graph.size(); ++right) {
            key.push_back(graph[left][right]);
        }
    }
    return key;
}

Graph unflatten(const Key& key, std::size_t vertices) {
    Graph graph(vertices, std::vector<int>(vertices, 0));
    std::size_t index = 0U;
    for (std::size_t left = 0U; left < vertices; ++left) {
        for (std::size_t right = left + 1U;
             right < vertices; ++right) {
            if (index >= key.size()) {
                throw std::runtime_error("short graph key");
            }
            graph[left][right] = key[index++];
        }
    }
    if (index != key.size()) {
        throw std::runtime_error("long graph key");
    }
    return graph;
}

Graph standard_graph(
    const std::vector<int>& degrees,
    const Key& top
) {
    std::vector<std::size_t> top_tokens;
    std::vector<std::size_t> bottom_tokens;
    for (std::size_t vertex = 0U; vertex < degrees.size(); ++vertex) {
        for (int copy = 0; copy < top[vertex]; ++copy) {
            top_tokens.push_back(vertex);
        }
        for (int copy = top[vertex]; copy < degrees[vertex]; ++copy) {
            bottom_tokens.push_back(vertex);
        }
    }
    if (top_tokens.size() != bottom_tokens.size()) {
        throw std::runtime_error("unbalanced standard graph");
    }
    Graph graph(degrees.size(), std::vector<int>(degrees.size(), 0));
    for (std::size_t token = 0U; token < top_tokens.size(); ++token) {
        const std::size_t left = top_tokens[token];
        const std::size_t right = bottom_tokens[token];
        if (left >= right) {
            throw std::runtime_error("nonstandard token pair");
        }
        ++graph[left][right];
    }
    return graph;
}

std::vector<Graph> ordered_standard_graphs(
    const std::vector<int>& degrees,
    const std::vector<std::size_t>& order
) {
    if (order.size() != degrees.size()) {
        throw std::runtime_error("ordered basis has wrong vertex count");
    }
    std::vector<bool> seen(order.size(), false);
    std::vector<int> ordered_degrees(order.size(), 0);
    for (std::size_t position = 0U; position < order.size(); ++position) {
        if (order[position] >= order.size() || seen[order[position]]) {
            throw std::runtime_error("ordered basis is not a permutation");
        }
        seen[order[position]] = true;
        ordered_degrees[position] = degrees[order[position]];
    }
    std::vector<Graph> answer;
    for (const Key& top : tops(ordered_degrees)) {
        const Graph ordered_graph = standard_graph(ordered_degrees, top);
        Graph graph(degrees.size(), std::vector<int>(degrees.size(), 0));
        for (std::size_t first = 0U; first < order.size(); ++first) {
            for (std::size_t second = first + 1U;
                 second < order.size(); ++second) {
                const std::size_t original_first
                    = std::min(order[first], order[second]);
                const std::size_t original_second
                    = std::max(order[first], order[second]);
                graph[original_first][original_second]
                    += ordered_graph[first][second];
            }
        }
        answer.push_back(std::move(graph));
    }
    return answer;
}

void add_scaled(
    Expansion& target,
    const Expansion& source,
    int sign
) {
    for (const auto& [key, coefficient] : source) {
        target[key] += sign * coefficient;
        if (target[key] == 0) {
            target.erase(key);
        }
    }
}

const Expansion& straighten(
    const Key& key,
    std::size_t vertices,
    std::map<Key, Expansion>& cache
) {
    const auto found = cache.find(key);
    if (found != cache.end()) {
        return found->second;
    }
    Graph graph = unflatten(key, vertices);
    for (std::size_t a = 0U; a < vertices; ++a) {
        for (std::size_t b = a + 1U; b < vertices; ++b) {
            for (std::size_t c = b + 1U; c < vertices; ++c) {
                for (std::size_t d = c + 1U; d < vertices; ++d) {
                    if (graph[a][d] == 0 || graph[b][c] == 0) {
                        continue;
                    }
                    --graph[a][d];
                    --graph[b][c];

                    Graph first = graph;
                    ++first[a][c];
                    ++first[b][d];
                    const Expansion first_expansion = straighten(
                        flatten(first), vertices, cache
                    );

                    Graph second = graph;
                    ++second[a][b];
                    ++second[c][d];
                    const Expansion second_expansion = straighten(
                        flatten(second), vertices, cache
                    );

                    Expansion answer;
                    add_scaled(answer, first_expansion, 1);
                    add_scaled(answer, second_expansion, -1);
                    return cache.emplace(key, std::move(answer)).first->second;
                }
            }
        }
    }
    return cache.emplace(key, Expansion{{key, 1}}).first->second;
}

Graph add_graphs(const Graph& first, const Graph& second) {
    if (first.size() != second.size()) {
        throw std::runtime_error("graph size mismatch");
    }
    Graph answer = first;
    for (std::size_t left = 0U; left < answer.size(); ++left) {
        for (std::size_t right = left + 1U;
             right < answer.size(); ++right) {
            answer[left][right] += second[left][right];
        }
    }
    return answer;
}

std::vector<Edge> edge_copies(
    const Graph& graph,
    unsigned int mask,
    bool inside
) {
    std::vector<Edge> edges;
    for (std::size_t first = 0U; first < graph.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < graph.size(); ++second) {
            const bool first_inside = first < 32U
                && ((mask >> first) & 1U) != 0U;
            const bool second_inside = second < 32U
                && ((mask >> second) & 1U) != 0U;
            if (first_inside != second_inside) {
                continue;
            }
            if (first_inside != inside) {
                continue;
            }
            for (int copy = 0; copy < graph[first][second]; ++copy) {
                edges.push_back({first, second});
            }
        }
    }
    return edges;
}

bool remove_edge(Graph& graph, const Edge& edge) {
    if (graph[edge.first][edge.second] == 0) {
        return false;
    }
    --graph[edge.first][edge.second];
    return true;
}

bool add_oriented_edge(
    Graph& graph,
    std::size_t first,
    std::size_t second,
    Integer& coefficient
) {
    if (first == second) {
        coefficient = 0;
        return false;
    }
    if (first > second) {
        std::swap(first, second);
        coefficient = -coefficient;
    }
    ++graph[first][second];
    return true;
}

Graph restrict_graph(const Graph& graph, unsigned int mask) {
    Graph answer(graph.size(), std::vector<int>(graph.size(), 0));
    for (std::size_t first = 0U; first < graph.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < graph.size(); ++second) {
            if (first < 32U && second < 32U
                && ((mask >> first) & 1U) != 0U
                && ((mask >> second) & 1U) != 0U) {
                answer[first][second] = graph[first][second];
            }
        }
    }
    return answer;
}

bool respects_cut(
    const Graph& graph,
    unsigned int mask,
    std::size_t factors
) {
    for (std::size_t first = 0U; first < graph.size(); ++first) {
        const bool first_inside = first < factors
            && ((mask >> first) & 1U) != 0U;
        for (std::size_t second = first + 1U;
             second < graph.size(); ++second) {
            const bool second_inside = second < factors
                && ((mask >> second) & 1U) != 0U;
            if (graph[first][second] != 0
                && first_inside != second_inside) {
                return false;
            }
        }
    }
    return true;
}

std::size_t exact_rank(std::vector<std::vector<Integer>> matrix) {
    if (matrix.empty()) {
        return 0U;
    }
    const std::size_t rows = matrix.size();
    const std::size_t columns = matrix.front().size();
    std::vector<std::vector<Rational>> rational(
        rows, std::vector<Rational>(columns)
    );
    for (std::size_t row = 0U; row < rows; ++row) {
        if (matrix[row].size() != columns) {
            throw std::runtime_error("ragged rank matrix");
        }
        for (std::size_t column = 0U; column < columns; ++column) {
            rational[row][column] = Rational(matrix[row][column]);
        }
    }
    std::size_t rank = 0U;
    for (std::size_t column = 0U;
         column < columns && rank < rows; ++column) {
        std::size_t pivot = rank;
        while (pivot < rows && rational[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == rows) {
            continue;
        }
        std::swap(rational[rank], rational[pivot]);
        const Rational pivot_value = rational[rank][column];
        for (std::size_t entry = column; entry < columns; ++entry) {
            rational[rank][entry] /= pivot_value;
        }
        for (std::size_t row = 0U; row < rows; ++row) {
            if (row == rank || rational[row][column] == 0) {
                continue;
            }
            const Rational factor = rational[row][column];
            for (std::size_t entry = column; entry < columns; ++entry) {
                rational[row][entry] -= factor * rational[rank][entry];
            }
        }
        ++rank;
    }
    return rank;
}

std::vector<std::vector<Rational>> exact_nullspace(
    const std::vector<std::vector<Integer>>& matrix
) {
    if (matrix.empty()) {
        return {};
    }
    const std::size_t rows = matrix.size();
    const std::size_t columns = matrix.front().size();
    std::vector<std::vector<Rational>> reduced(
        rows, std::vector<Rational>(columns)
    );
    for (std::size_t row = 0U; row < rows; ++row) {
        if (matrix[row].size() != columns) {
            throw std::runtime_error("ragged nullspace matrix");
        }
        for (std::size_t column = 0U; column < columns; ++column) {
            reduced[row][column] = Rational(matrix[row][column]);
        }
    }
    std::vector<std::size_t> pivots;
    std::size_t rank = 0U;
    for (std::size_t column = 0U;
         column < columns && rank < rows; ++column) {
        std::size_t pivot = rank;
        while (pivot < rows && reduced[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == rows) {
            continue;
        }
        std::swap(reduced[rank], reduced[pivot]);
        const Rational pivot_value = reduced[rank][column];
        for (std::size_t entry = column; entry < columns; ++entry) {
            reduced[rank][entry] /= pivot_value;
        }
        for (std::size_t row = 0U; row < rows; ++row) {
            if (row == rank || reduced[row][column] == 0) {
                continue;
            }
            const Rational factor = reduced[row][column];
            for (std::size_t entry = column; entry < columns; ++entry) {
                reduced[row][entry] -= factor * reduced[rank][entry];
            }
        }
        pivots.push_back(column);
        ++rank;
    }
    std::vector<bool> is_pivot(columns, false);
    for (std::size_t pivot : pivots) {
        is_pivot[pivot] = true;
    }
    std::vector<std::vector<Rational>> basis;
    for (std::size_t free = 0U; free < columns; ++free) {
        if (is_pivot[free]) {
            continue;
        }
        std::vector<Rational> vector(columns, Rational(0));
        vector[free] = Rational(1);
        for (std::size_t row = 0U; row < pivots.size(); ++row) {
            vector[pivots[row]] = -reduced[row][free];
        }
        basis.push_back(std::move(vector));
    }
    return basis;
}

std::pair<bool, std::vector<Rational>> solve_left(
    const std::vector<std::vector<Integer>>& rows,
    const std::vector<Rational>& target
) {
    const std::size_t variables = rows.size();
    const std::size_t equations = target.size();
    if (variables != 0U && rows.front().size() != equations) {
        throw std::runtime_error("left solve dimension mismatch");
    }
    std::vector<std::vector<Rational>> augmented(
        equations, std::vector<Rational>(variables + 1U, Rational(0))
    );
    for (std::size_t equation = 0U; equation < equations; ++equation) {
        for (std::size_t variable = 0U;
             variable < variables; ++variable) {
            if (rows[variable].size() != equations) {
                throw std::runtime_error("ragged left solve matrix");
            }
            augmented[equation][variable] = rows[variable][equation];
        }
        augmented[equation][variables] = target[equation];
    }
    std::vector<std::size_t> pivots;
    std::size_t rank = 0U;
    for (std::size_t variable = 0U;
         variable < variables && rank < equations; ++variable) {
        std::size_t pivot = rank;
        while (pivot < equations && augmented[pivot][variable] == 0) {
            ++pivot;
        }
        if (pivot == equations) {
            continue;
        }
        std::swap(augmented[rank], augmented[pivot]);
        const Rational pivot_value = augmented[rank][variable];
        for (std::size_t entry = variable;
             entry <= variables; ++entry) {
            augmented[rank][entry] /= pivot_value;
        }
        for (std::size_t equation = 0U;
             equation < equations; ++equation) {
            if (equation == rank || augmented[equation][variable] == 0) {
                continue;
            }
            const Rational factor = augmented[equation][variable];
            for (std::size_t entry = variable;
                 entry <= variables; ++entry) {
                augmented[equation][entry]
                    -= factor * augmented[rank][entry];
            }
        }
        pivots.push_back(variable);
        ++rank;
    }
    for (std::size_t equation = rank; equation < equations; ++equation) {
        bool zero = true;
        for (std::size_t variable = 0U;
             variable < variables; ++variable) {
            zero = zero && augmented[equation][variable] == 0;
        }
        if (zero && augmented[equation][variables] != 0) {
            return {false, {}};
        }
    }
    std::vector<Rational> solution(variables, Rational(0));
    for (std::size_t row = 0U; row < pivots.size(); ++row) {
        solution[pivots[row]] = augmented[row][variables];
    }
    return {true, std::move(solution)};
}

std::vector<Integer> integerize(const std::vector<Rational>& vector) {
    Integer denominator_product = 1;
    for (const Rational& coefficient : vector) {
        denominator_product *= coefficient.denominator();
    }
    std::vector<Integer> answer(vector.size(), 0);
    for (std::size_t entry = 0U; entry < vector.size(); ++entry) {
        answer[entry] = vector[entry].numerator()
            * (denominator_product / vector[entry].denominator());
    }
    return answer;
}

std::vector<std::vector<Integer>> column_matrix(
    const std::vector<std::vector<Integer>>& columns,
    std::size_t rows
);

std::vector<std::vector<Integer>> apply_row_map(
    const std::vector<std::vector<Integer>>& rows,
    const std::vector<std::vector<Integer>>& linear_map
) {
    if (linear_map.empty()) {
        if (rows.empty()) {
            return {};
        }
        if (!rows.front().empty()) {
            throw std::runtime_error("empty row map has nonempty domain");
        }
        return std::vector<std::vector<Integer>>(rows.size());
    }
    const std::size_t domain = linear_map.size();
    const std::size_t codomain = linear_map.front().size();
    for (const auto& image : linear_map) {
        if (image.size() != codomain) {
            throw std::runtime_error("ragged row map");
        }
    }
    std::vector<std::vector<Integer>> answer(
        rows.size(), std::vector<Integer>(codomain, 0)
    );
    for (std::size_t row = 0U; row < rows.size(); ++row) {
        if (rows[row].size() != domain) {
            throw std::runtime_error("row-map domain mismatch");
        }
        for (std::size_t source = 0U; source < domain; ++source) {
            if (rows[row][source] == 0) {
                continue;
            }
            for (std::size_t target = 0U; target < codomain; ++target) {
                answer[row][target]
                    += rows[row][source] * linear_map[source][target];
            }
        }
    }
    return answer;
}

std::vector<std::vector<Integer>> span_intersection(
    const std::vector<std::vector<Integer>>& first,
    const std::vector<std::vector<Integer>>& second,
    std::size_t ambient_dimension
) {
    if (first.empty() || second.empty()) {
        return {};
    }
    std::vector<std::vector<Integer>> combined = first;
    combined.insert(combined.end(), second.begin(), second.end());
    const auto relations = exact_nullspace(
        column_matrix(combined, ambient_dimension)
    );
    std::vector<std::vector<Integer>> basis;
    std::size_t basis_rank = 0U;
    for (const std::vector<Rational>& relation : relations) {
        std::vector<Rational> target(ambient_dimension, Rational(0));
        for (std::size_t source = 0U; source < first.size(); ++source) {
            for (std::size_t coordinate = 0U;
                 coordinate < ambient_dimension; ++coordinate) {
                target[coordinate]
                    += relation[source] * Rational(first[source][coordinate]);
            }
        }
        if (std::all_of(
                target.begin(), target.end(),
                [](const Rational& value) { return value == 0; }
            )) {
            continue;
        }
        std::vector<Integer> integer_target = integerize(target);
        auto trial = basis;
        trial.push_back(integer_target);
        const std::size_t next_rank = exact_rank(trial);
        if (next_rank != basis_rank) {
            basis.push_back(std::move(integer_target));
            basis_rank = next_rank;
        }
    }
    return basis;
}

Graph delete_graph_vertex(const Graph& graph, std::size_t deleted) {
    if (deleted >= graph.size()) {
        throw std::runtime_error("deleted graph vertex is out of range");
    }
    Graph answer(
        graph.size() - 1U,
        std::vector<int>(graph.size() - 1U, 0)
    );
    for (std::size_t first = 0U; first < graph.size(); ++first) {
        if (first == deleted) {
            continue;
        }
        const std::size_t new_first
            = first < deleted ? first : first - 1U;
        for (std::size_t second = first + 1U;
             second < graph.size(); ++second) {
            if (second == deleted) {
                continue;
            }
            const std::size_t new_second
                = second < deleted ? second : second - 1U;
            answer[new_first][new_second] = graph[first][second];
        }
    }
    return answer;
}

struct Result {
    std::uint64_t odd_sources = 0U;
    std::uint64_t positive_sources = 0U;
    std::size_t codomain_dimension = 0U;
    std::size_t rank = 0U;
    std::size_t straightening_states = 0U;
    struct PairIntersection {
        unsigned int first = 0U;
        unsigned int second = 0U;
        std::size_t dimension = 0U;
        std::uint64_t xor_channel_dimension = 0U;
    };
    std::vector<PairIntersection> pair_intersections;
    std::vector<std::pair<unsigned int, std::size_t>> column_labels;
    std::vector<Graph> column_graphs;
    struct PositiveSource {
        unsigned int mask = 0U;
        std::size_t basis = 0U;
        Graph graph;
    };
    std::vector<PositiveSource> positive_source_graphs;
    std::vector<std::vector<Rational>> nullspace;
    std::size_t raw_collision_relations = 0U;
    std::size_t pair_intersection_relations = 0U;
    std::size_t edge_exchange_relations = 0U;
    std::size_t cut_boundary_relations = 0U;
    std::size_t cut_boundary_rank = 0U;
    std::vector<std::size_t> cut_order_ranks;
    std::size_t cut_order_combined_rank = 0U;
    std::size_t section_defect_relations = 0U;
    std::vector<std::size_t> section_trial_ranks;
    std::vector<std::size_t> section_individual_ranks;
    struct SectionCriterion {
        std::size_t rank = 0U;
        std::size_t kernel = 0U;
        std::size_t carrier_intersection = 0U;
        bool pass = false;
    };
    std::vector<SectionCriterion> section_criteria;
    std::size_t candidate_relation_rank = 0U;
    bool candidate_relations_valid = true;
    struct FundamentalLiftTerm {
        std::array<std::size_t, 3> triple{};
        Graph coefficient;
        Rational scalar;
        struct Owner {
            unsigned int mask = 0U;
            std::size_t basis = 0U;
            std::size_t face = 0U;
            Rational ratio;
        };
        std::vector<Owner> owners;
    };
    std::vector<std::vector<FundamentalLiftTerm>> fundamental_lifts;
    bool fundamental_lifts_valid = true;
    bool fundamental_owned_lifts_valid = true;
    bool fundamental_combined_valid = true;
    std::size_t fundamental_kernel_quotient_rank = 0U;
    std::size_t fundamental_candidate_quotient_rank = 0U;
    std::size_t fundamental_covered_rank = 0U;
    std::vector<Integer> fundamental_missing_relation;
    std::size_t positive_internal_kernel_rank = 0U;
    std::size_t optimal_section_rank = 0U;
    bool optimal_section_criterion = true;
    std::size_t allocation_dimension = 0U;
    std::size_t carrier_dimension = 0U;
    std::size_t f0_dimension = 0U;
    std::size_t f0_kernel_dimension = 0U;
    std::size_t negative_stripped_rank = 0U;
    std::size_t positive_stripped_rank = 0U;
    std::size_t positive_global_rank = 0U;
    std::size_t carrier_overlap_dimension = 0U;
    bool ambient_kernel_criterion = true;
    std::size_t positive_allocation_dimension = 0U;
    std::size_t unshared_negative_carrier_dimension = 0U;
    std::size_t residual_negative_demand = 0U;
    std::size_t residual_positive_supply = 0U;
    bool optimal_section_cover_feasible = true;
    std::size_t optimal_defect_rank = 0U;
    std::size_t optimal_defect_kernel = 0U;
    std::size_t optimal_relation_rank = 0U;
    std::size_t merged_row_kernel_dimension = 0U;
    bool deletion_contraction_valid = true;
    std::size_t projected_allocation_dimension = 0U;
    std::size_t projected_positive_allocation_dimension = 0U;
    std::size_t projected_unshared_carrier_dimension = 0U;
    std::size_t projected_positive_rank = 0U;
    std::size_t negative_lower_branch_dimension = 0U;
    std::size_t positive_lower_branch_dimension = 0U;
    std::size_t negative_lift_obstruction_rank = 0U;
    std::size_t positive_lift_obstruction_rank = 0U;
    bool enriched_deletion_contraction_valid = true;
    bool projected_carrier_criterion = true;
};

std::vector<std::vector<Integer>> column_matrix(
    const std::vector<std::vector<Integer>>& columns,
    std::size_t rows
) {
    std::vector<std::vector<Integer>> matrix(
        rows, std::vector<Integer>(columns.size(), 0)
    );
    for (std::size_t column = 0U; column < columns.size(); ++column) {
        if (columns[column].size() != rows) {
            throw std::runtime_error("column dimension mismatch");
        }
        for (std::size_t row = 0U; row < rows; ++row) {
            matrix[row][column] = columns[column][row];
        }
    }
    return matrix;
}

Result analyze_case(
    const std::vector<int>& labels,
    int target,
    unsigned int minus_mask,
    bool diagnostics = true,
    bool relations = false,
    bool fundamental = false
) {
    const std::size_t factors = labels.size();
    std::vector<int> degrees = labels;
    degrees.push_back(target);
    const std::vector<Key> global_tops = tops(degrees);
    std::map<Key, std::size_t> basis_index;
    for (const Key& top : global_tops) {
        const Key graph_key = flatten(standard_graph(degrees, top));
        basis_index.emplace(graph_key, basis_index.size());
    }
    if (basis_index.size() != global_tops.size()) {
        throw std::runtime_error("duplicate global standard graph");
    }

    std::vector<std::vector<Integer>> columns;
    std::map<unsigned int, std::vector<std::vector<Integer>>> odd_columns;
    std::map<unsigned int, std::vector<std::size_t>> odd_channel_sources;
    std::map<ChannelKey, std::size_t> odd_source_indices;
    std::map<Key, std::vector<std::size_t>> raw_odd_sources;
    std::map<unsigned int, std::uint64_t> channel_dimensions;
    std::map<Key, Expansion> straightening_cache;
    const unsigned int subset_limit = 1U << factors;
    Result result;
    result.codomain_dimension = basis_index.size();
    for (unsigned int subset = 0U; subset < subset_limit; ++subset) {
        std::vector<int> left_degrees(factors + 1U, 0);
        std::vector<int> right_degrees(factors + 1U, 0);
        for (std::size_t factor = 0U; factor < factors; ++factor) {
            if (((subset >> factor) & 1U) != 0U) {
                left_degrees[factor] = labels[factor];
            } else {
                right_degrees[factor] = labels[factor];
            }
        }
        right_degrees[factors] = target;
        const std::vector<Key> left_tops = tops(left_degrees);
        const std::vector<Key> right_tops = tops(right_degrees);
        const std::uint64_t channel_dimension
            = static_cast<std::uint64_t>(left_tops.size())
                * static_cast<std::uint64_t>(right_tops.size());
        const bool odd = (popcount(subset & minus_mask) & 1) != 0;
        channel_dimensions[subset] = channel_dimension;
        if (odd) {
            result.odd_sources += channel_dimension;
        } else if (subset != 0U) {
            result.positive_sources += channel_dimension;
        }
        if (!odd) {
            if ((diagnostics || relations || fundamental) && subset != 0U) {
                std::size_t basis = 0U;
                for (const Key& left_top : left_tops) {
                    const Graph left_graph = standard_graph(
                        left_degrees, left_top
                    );
                    for (const Key& right_top : right_tops) {
                        const Graph right_graph = standard_graph(
                            right_degrees, right_top
                        );
                        result.positive_source_graphs.push_back(
                            Result::PositiveSource{
                                subset,
                                basis++,
                                add_graphs(left_graph, right_graph)
                            }
                        );
                    }
                }
            }
            continue;
        }
        for (const Key& left_top : left_tops) {
            const Graph left_graph = standard_graph(
                left_degrees, left_top
            );
            for (const Key& right_top : right_tops) {
                const Graph right_graph = standard_graph(
                    right_degrees, right_top
                );
                const Key product_key = flatten(
                    add_graphs(left_graph, right_graph)
                );
                const Expansion& expansion = straighten(
                    product_key, degrees.size(), straightening_cache
                );
                std::vector<Integer> column(basis_index.size(), 0);
                for (const auto& [key, coefficient] : expansion) {
                    const auto row = basis_index.find(key);
                    if (row == basis_index.end()) {
                        throw std::runtime_error(
                            "straightening left the standard basis"
                        );
                    }
                    column[row->second] += coefficient;
                }
                odd_columns[subset].push_back(column);
                const std::size_t global_source = columns.size();
                const ChannelKey source_key{
                    subset, flatten(left_graph), flatten(right_graph)
                };
                if (!odd_source_indices.emplace(
                        source_key, global_source
                    ).second) {
                    throw std::runtime_error("duplicate odd source basis key");
                }
                odd_channel_sources[subset].push_back(global_source);
                raw_odd_sources[product_key].push_back(global_source);
                result.column_labels.emplace_back(
                    subset, odd_columns[subset].size() - 1U
                );
                if (diagnostics || fundamental) {
                    result.column_graphs.push_back(
                        add_graphs(left_graph, right_graph)
                    );
                }
                columns.push_back(std::move(column));
            }
        }
    }

    const std::vector<std::vector<Integer>> full_matrix
        = column_matrix(columns, basis_index.size());
    result.rank = exact_rank(full_matrix);
    std::vector<std::vector<Integer>> candidate_relations;
    if (diagnostics || relations || fundamental) {
        for (const auto& [raw_graph, sources] : raw_odd_sources) {
            (void)raw_graph;
            for (std::size_t index = 1U; index < sources.size(); ++index) {
                std::vector<Integer> relation(columns.size(), 0);
                relation[sources.front()] = -1;
                relation[sources[index]] = 1;
                candidate_relations.push_back(std::move(relation));
                ++result.raw_collision_relations;
            }
        }
        for (auto first = odd_columns.begin();
             first != odd_columns.end(); ++first) {
            for (auto second = std::next(first);
                 second != odd_columns.end(); ++second) {
                std::vector<std::vector<Integer>> combined = first->second;
                combined.insert(
                    combined.end(), second->second.begin(), second->second.end()
                );
                const auto local_relations = exact_nullspace(
                    column_matrix(combined, basis_index.size())
                );
                std::vector<std::size_t> global_sources
                    = odd_channel_sources.at(first->first);
                const auto& second_sources
                    = odd_channel_sources.at(second->first);
                global_sources.insert(
                    global_sources.end(),
                    second_sources.begin(), second_sources.end()
                );
                for (const std::vector<Rational>& local : local_relations) {
                    Integer denominator_product = 1;
                    for (const Rational& coefficient : local) {
                        denominator_product *= coefficient.denominator();
                    }
                    std::vector<Integer> relation(columns.size(), 0);
                    for (std::size_t source = 0U;
                         source < local.size(); ++source) {
                        relation[global_sources[source]]
                            = local[source].numerator()
                                * (denominator_product
                                   / local[source].denominator());
                    }
                    candidate_relations.push_back(std::move(relation));
                    ++result.pair_intersection_relations;
                }
            }
        }

        const auto lift_to_odd_subset = [
            &odd_source_indices,
            &straightening_cache,
            factors,
            minus_mask,
            source_count = columns.size()
        ](
            const Graph& graph,
            unsigned int subset,
            std::vector<Integer>& lift
        ) {
            if ((popcount(subset & minus_mask) & 1) == 0
                || !respects_cut(graph, subset, factors)) {
                return false;
            }
            const Graph left = restrict_graph(graph, subset);
            Graph right = graph;
            for (std::size_t first = 0U;
                 first < graph.size(); ++first) {
                for (std::size_t second = first + 1U;
                     second < graph.size(); ++second) {
                    if (first < factors && second < factors
                        && ((subset >> first) & 1U) != 0U
                        && ((subset >> second) & 1U) != 0U) {
                        right[first][second] = 0;
                    }
                }
            }
            const Key left_key = flatten(left);
            const Key right_key = flatten(right);
            const Expansion& left_expansion = straighten(
                left_key, graph.size(), straightening_cache
            );
            const Expansion& right_expansion = straighten(
                right_key, graph.size(), straightening_cache
            );
            std::vector<Integer> answer(source_count, 0);
            for (const auto& [left_graph_key, left_coefficient]
                 : left_expansion) {
                for (const auto& [right_graph_key, right_coefficient]
                     : right_expansion) {
                    const auto found = odd_source_indices.find(
                        ChannelKey{
                            subset, left_graph_key, right_graph_key
                        }
                    );
                    if (found == odd_source_indices.end()) {
                        throw std::runtime_error(
                            "odd lift left the channel basis"
                        );
                    }
                    answer[found->second]
                        += left_coefficient * right_coefficient;
                }
            }
            lift = std::move(answer);
            return true;
        };

        const auto lift_to_odd = [
            &lift_to_odd_subset,
            subset_limit
        ](
            const Graph& graph,
            std::vector<Integer>& lift
        ) {
            for (unsigned int subset = 0U;
                 subset < subset_limit; ++subset) {
                if (lift_to_odd_subset(graph, subset, lift)) {
                    return true;
                }
            }
            return false;
        };

        std::vector<std::vector<Integer>> cut_boundary_relations;
        if (popcount(minus_mask) == 1) {
            std::size_t fundamental = 0U;
            while (((minus_mask >> fundamental) & 1U) == 0U) {
                ++fundamental;
            }
            if (labels[fundamental] == 1) {
                const auto build_cut_boundary = [
                    &lift_to_odd_subset,
                    factors,
                    fundamental
                ](const std::vector<Result::PositiveSource>& sources) {
                    std::vector<std::vector<Integer>> relations;
                    for (const Result::PositiveSource& source : sources) {
                        Graph stripped = source.graph;
                        std::size_t neighbor = stripped.size();
                        for (std::size_t vertex = 0U;
                             vertex < stripped.size(); ++vertex) {
                            if (vertex == fundamental) {
                                continue;
                            }
                            const std::size_t first
                                = std::min(vertex, fundamental);
                            const std::size_t second
                                = std::max(vertex, fundamental);
                            if (stripped[first][second] == 0) {
                                continue;
                            }
                            if (stripped[first][second] != 1
                                || neighbor != stripped.size()) {
                                throw std::runtime_error(
                                    "positive fundamental edge is not unique"
                                );
                            }
                            neighbor = vertex;
                            stripped[first][second] = 0;
                        }
                        if (neighbor == stripped.size()) {
                            throw std::runtime_error(
                                "positive fundamental edge is missing"
                            );
                        }
                        std::vector<bool> component(
                            stripped.size(), false
                        );
                        std::vector<std::size_t> queue{neighbor};
                        component[neighbor] = true;
                        for (std::size_t head = 0U;
                             head < queue.size(); ++head) {
                            const std::size_t vertex = queue[head];
                            for (std::size_t other = 0U;
                                 other < stripped.size(); ++other) {
                                const std::size_t first
                                    = std::min(vertex, other);
                                const std::size_t second
                                    = std::max(vertex, other);
                                if (first == second
                                    || stripped[first][second] == 0
                                    || component[other]) {
                                    continue;
                                }
                                component[other] = true;
                                queue.push_back(other);
                            }
                        }
                        if (component[factors]) {
                            continue;
                        }
                        unsigned int component_mask = 0U;
                        for (std::size_t vertex = 0U;
                             vertex < factors; ++vertex) {
                            if (component[vertex]) {
                                component_mask |= 1U << vertex;
                            }
                        }
                        if ((component_mask & source.mask) != 0U) {
                            throw std::runtime_error(
                                "positive cut meets neighbor component"
                            );
                        }
                        const unsigned int base_mask
                            = (1U << fundamental) | component_mask;
                        const unsigned int extended_mask
                            = base_mask | source.mask;
                        std::vector<Integer> base_lift;
                        std::vector<Integer> extended_lift;
                        if (!lift_to_odd_subset(
                                source.graph, base_mask, base_lift
                            )
                            || !lift_to_odd_subset(
                                source.graph,
                                extended_mask,
                                extended_lift
                            )) {
                            throw std::runtime_error(
                                "cut-boundary odd lift failed"
                            );
                        }
                        for (std::size_t entry = 0U;
                             entry < base_lift.size(); ++entry) {
                            extended_lift[entry] -= base_lift[entry];
                        }
                        if (std::any_of(
                                extended_lift.begin(),
                                extended_lift.end(),
                                [](const Integer& value) {
                                    return value != 0;
                                }
                            )) {
                            relations.push_back(std::move(extended_lift));
                        }
                    }
                    return relations;
                };
                cut_boundary_relations = build_cut_boundary(
                    result.positive_source_graphs
                );
                result.cut_boundary_relations
                    = cut_boundary_relations.size();
                candidate_relations.insert(
                    candidate_relations.end(),
                    cut_boundary_relations.begin(),
                    cut_boundary_relations.end()
                );

                if (diagnostics) {
                    std::vector<std::vector<std::size_t>> trial_orders;
                    std::vector<std::vector<Integer>> all_order_relations;
                    const std::size_t vertices = degrees.size();
                    for (std::size_t shift = 0U;
                         shift < vertices; ++shift) {
                    std::vector<std::size_t> order(vertices, 0U);
                    for (std::size_t position = 0U;
                         position < vertices; ++position) {
                        order[position] = (position + shift) % vertices;
                    }
                    trial_orders.push_back(std::move(order));
                    }
                    for (std::size_t target_position = 0U;
                         target_position < vertices; ++target_position) {
                        for (std::size_t fundamental_position = 0U;
                             fundamental_position < vertices;
                             ++fundamental_position) {
                            if (target_position == fundamental_position) {
                                continue;
                            }
                            std::vector<std::size_t> order;
                            order.reserve(vertices);
                            std::size_t next_other = 0U;
                            for (std::size_t position = 0U;
                                 position < vertices; ++position) {
                                if (position == target_position) {
                                    order.push_back(factors);
                                } else if (position
                                           == fundamental_position) {
                                    order.push_back(fundamental);
                                } else {
                                    while (next_other == fundamental
                                           || next_other == factors) {
                                        ++next_other;
                                    }
                                    order.push_back(next_other++);
                                }
                            }
                            trial_orders.push_back(std::move(order));
                        }
                    }
                    for (const std::vector<std::size_t>& order
                         : trial_orders) {
                        std::vector<Result::PositiveSource> ordered_sources;
                        for (unsigned int subset = 1U;
                             subset < subset_limit; ++subset) {
                            if ((popcount(subset & minus_mask) & 1) != 0) {
                                continue;
                            }
                            std::vector<int> left_degrees(
                                factors + 1U, 0
                            );
                            std::vector<int> right_degrees(
                                factors + 1U, 0
                            );
                            for (std::size_t factor = 0U;
                                 factor < factors; ++factor) {
                                if (((subset >> factor) & 1U) != 0U) {
                                    left_degrees[factor] = labels[factor];
                                } else {
                                    right_degrees[factor] = labels[factor];
                                }
                            }
                            right_degrees[factors] = target;
                            std::size_t basis = 0U;
                            const auto left_graphs = ordered_standard_graphs(
                                left_degrees, order
                            );
                            const auto right_graphs
                                = ordered_standard_graphs(
                                    right_degrees, order
                                );
                            for (const Graph& left : left_graphs) {
                                for (const Graph& right : right_graphs) {
                                    ordered_sources.push_back(
                                        Result::PositiveSource{
                                            subset,
                                            basis++,
                                            add_graphs(left, right)
                                        }
                                    );
                                }
                            }
                        }
                        auto ordered_relations
                            = build_cut_boundary(ordered_sources);
                        result.cut_order_ranks.push_back(exact_rank(
                            ordered_relations
                        ));
                        all_order_relations.insert(
                            all_order_relations.end(),
                            ordered_relations.begin(),
                            ordered_relations.end()
                        );
                    }
                    result.cut_order_combined_rank
                        = exact_rank(std::move(all_order_relations));
                }
            }
        }
        result.cut_boundary_rank = exact_rank(cut_boundary_relations);

        const auto add_edge_exchange = [
            &candidate_relations,
            &lift_to_odd,
            &result,
            source_count = columns.size()
        ](
            const Graph& positive_graph,
            const Edge& ab,
            const Edge& cd,
            const Edge& ef
        ) {
            Graph residual = positive_graph;
            if (!remove_edge(residual, ab)
                || !remove_edge(residual, cd)
                || !remove_edge(residual, ef)) {
                throw std::runtime_error("edge-exchange removal failed");
            }
            const std::size_t a = ab.first;
            const std::size_t b = ab.second;
            const std::size_t c = cd.first;
            const std::size_t d = cd.second;
            const std::size_t e = ef.first;
            const std::size_t f = ef.second;
            const std::array<std::array<std::pair<std::size_t, std::size_t>, 3>, 4>
                terms{{
                    {{{a, c}, {b, d}, {e, f}}},
                    {{{a, d}, {b, c}, {e, f}}},
                    {{{a, e}, {b, f}, {c, d}}},
                    {{{a, f}, {b, e}, {c, d}}}
                }};
            const std::array<int, 4> signs{{1, -1, -1, 1}};
            std::vector<Integer> relation(source_count, 0);
            for (std::size_t term_index = 0U;
                 term_index < terms.size(); ++term_index) {
                Graph term = residual;
                Integer coefficient = signs[term_index];
                for (const auto& [first, second] : terms[term_index]) {
                    if (!add_oriented_edge(
                            term, first, second, coefficient
                        )) {
                        break;
                    }
                }
                if (coefficient == 0) {
                    continue;
                }
                std::vector<Integer> term_lift;
                if (!lift_to_odd(term, term_lift)) {
                    return;
                }
                for (std::size_t source = 0U;
                     source < source_count; ++source) {
                    relation[source] += coefficient * term_lift[source];
                }
            }
            if (std::any_of(
                    relation.begin(), relation.end(),
                    [](const Integer& value) { return value != 0; }
                )) {
                candidate_relations.push_back(std::move(relation));
                ++result.edge_exchange_relations;
            }
        };

        for (const Result::PositiveSource& source
             : result.positive_source_graphs) {
            const std::vector<Edge> left_edges = edge_copies(
                source.graph, source.mask, true
            );
            const std::vector<Edge> right_edges = edge_copies(
                source.graph, source.mask, false
            );
            for (const Edge& first : left_edges) {
                for (std::size_t second = 0U;
                     second < right_edges.size(); ++second) {
                    for (std::size_t third = second + 1U;
                         third < right_edges.size(); ++third) {
                        add_edge_exchange(
                            source.graph, first,
                            right_edges[second], right_edges[third]
                        );
                    }
                }
            }
            for (const Edge& first : right_edges) {
                for (std::size_t second = 0U;
                     second < left_edges.size(); ++second) {
                    for (std::size_t third = second + 1U;
                         third < left_edges.size(); ++third) {
                        add_edge_exchange(
                            source.graph, first,
                            left_edges[second], left_edges[third]
                        );
                    }
                }
            }
        }
        for (const std::vector<Integer>& relation : candidate_relations) {
            for (std::size_t row = 0U; row < full_matrix.size(); ++row) {
                Integer value = 0;
                for (std::size_t source = 0U;
                     source < relation.size(); ++source) {
                    value += full_matrix[row][source] * relation[source];
                }
                if (value != 0) {
                    result.candidate_relations_valid = false;
                }
            }
        }
        result.candidate_relation_rank = exact_rank(candidate_relations);
    }
    if (diagnostics || fundamental) {
        result.nullspace = exact_nullspace(full_matrix);
        if (popcount(minus_mask) == 1) {
            std::size_t fundamental = 0U;
            while (((minus_mask >> fundamental) & 1U) == 0U) {
                ++fundamental;
            }
            if (labels[fundamental] == 1) {
                using F0Key = std::pair<std::size_t, Key>;
                struct F1Element {
                    std::array<std::size_t, 3> triple{};
                    Graph coefficient;
                };
                std::vector<F0Key> f0_basis;
                std::map<F0Key, std::size_t> f0_indices;
                for (std::size_t neighbor = 0U;
                     neighbor < degrees.size(); ++neighbor) {
                    if (neighbor == fundamental || degrees[neighbor] == 0) {
                        continue;
                    }
                    std::vector<int> residual = degrees;
                    residual[fundamental] = 0;
                    --residual[neighbor];
                    for (const Key& top : tops(residual)) {
                        const Key graph_key = flatten(
                            standard_graph(residual, top)
                        );
                        const F0Key key{neighbor, graph_key};
                        f0_indices.emplace(key, f0_basis.size());
                        f0_basis.push_back(key);
                    }
                }
                result.f0_dimension = f0_basis.size();
                std::vector<F1Element> f1_basis;
                for (std::size_t first = 0U;
                     first < degrees.size(); ++first) {
                    if (first == fundamental) {
                        continue;
                    }
                    for (std::size_t second = first + 1U;
                         second < degrees.size(); ++second) {
                        if (second == fundamental) {
                            continue;
                        }
                        for (std::size_t third = second + 1U;
                             third < degrees.size(); ++third) {
                            if (third == fundamental) {
                                continue;
                            }
                            std::vector<int> residual = degrees;
                            residual[fundamental] = 0;
                            --residual[first];
                            --residual[second];
                            --residual[third];
                            if (std::any_of(
                                    residual.begin(), residual.end(),
                                    [](int value) { return value < 0; }
                                )) {
                                continue;
                            }
                            for (const Key& top : tops(residual)) {
                                f1_basis.push_back(F1Element{
                                    {first, second, third},
                                    standard_graph(residual, top)
                                });
                            }
                        }
                    }
                }
                std::vector<std::vector<Integer>> differential(
                    f1_basis.size(),
                    std::vector<Integer>(f0_basis.size(), 0)
                );
                std::vector<std::array<std::vector<Integer>, 3>> faces(
                    f1_basis.size()
                );
                for (auto& source_faces : faces) {
                    for (auto& face : source_faces) {
                        face.assign(f0_basis.size(), 0);
                    }
                }
                const auto add_f0_term = [
                    &differential,
                    &faces,
                    &f0_indices,
                    &straightening_cache
                ](
                    std::size_t source,
                    std::size_t face,
                    const Graph& coefficient,
                    std::size_t edge_first,
                    std::size_t edge_second,
                    std::size_t output_neighbor,
                    int scalar
                ) {
                    Graph product = coefficient;
                    Integer orientation = scalar;
                    add_oriented_edge(
                        product, edge_first, edge_second, orientation
                    );
                    const Key product_key = flatten(product);
                    for (const auto& [graph_key, expansion_coefficient]
                         : straighten(
                               product_key, product.size(),
                               straightening_cache
                           )) {
                        const auto found = f0_indices.find(
                            F0Key{output_neighbor, graph_key}
                        );
                        if (found == f0_indices.end()) {
                            throw std::runtime_error(
                                "fundamental differential left F0"
                            );
                        }
                        const Integer value
                            = orientation * expansion_coefficient;
                        differential[source][found->second] += value;
                        faces[source][face][found->second] += value;
                    }
                };
                for (std::size_t source = 0U;
                     source < f1_basis.size(); ++source) {
                    const auto [first, second, third]
                        = f1_basis[source].triple;
                    add_f0_term(
                        source, 0U, f1_basis[source].coefficient,
                        second, third, first, 1
                    );
                    add_f0_term(
                        source, 1U, f1_basis[source].coefficient,
                        first, third, second, -1
                    );
                    add_f0_term(
                        source, 2U, f1_basis[source].coefficient,
                        first, second, third, 1
                    );
                }
                const auto strip_to_f0 = [
                    &degrees,
                    fundamental,
                    &f0_indices,
                    &straightening_cache
                ](const Graph& graph) {
                    std::vector<Integer> answer(f0_indices.size(), 0);
                    Graph stripped = graph;
                    std::size_t neighbor = degrees.size();
                    for (std::size_t vertex = 0U;
                         vertex < degrees.size(); ++vertex) {
                        if (vertex == fundamental) {
                            continue;
                        }
                        const std::size_t left
                            = std::min(vertex, fundamental);
                        const std::size_t right
                            = std::max(vertex, fundamental);
                        if (stripped[left][right] != 0) {
                            if (stripped[left][right] != 1
                                || neighbor != degrees.size()) {
                                throw std::runtime_error(
                                    "fundamental edge is not unique"
                                );
                            }
                            neighbor = vertex;
                            stripped[left][right] = 0;
                        }
                    }
                    if (neighbor == degrees.size()) {
                        throw std::runtime_error(
                            "fundamental edge is missing"
                        );
                    }
                    const Integer orientation
                        = fundamental < neighbor ? 1 : -1;
                    const Key stripped_key = flatten(stripped);
                    for (const auto& [graph_key, expansion_coefficient]
                         : straighten(
                               stripped_key, stripped.size(),
                               straightening_cache
                           )) {
                        const auto found = f0_indices.find(
                            F0Key{neighbor, graph_key}
                        );
                        if (found == f0_indices.end()) {
                            throw std::runtime_error(
                                "fundamental stripping left F0"
                            );
                        }
                        answer[found->second]
                            += orientation * expansion_coefficient;
                    }
                    return answer;
                };
                std::vector<std::vector<Integer>> iota(
                    columns.size(),
                    std::vector<Integer>(f0_basis.size(), 0)
                );
                for (std::size_t source = 0U;
                     source < result.column_graphs.size(); ++source) {
                    iota[source] = strip_to_f0(
                        result.column_graphs[source]
                    );
                }
                std::vector<std::vector<Integer>> positive_iota;
                positive_iota.reserve(result.positive_source_graphs.size());
                for (const Result::PositiveSource& source
                     : result.positive_source_graphs) {
                    positive_iota.push_back(strip_to_f0(source.graph));
                }
                const std::vector<std::vector<Integer>> section_base_relations
                    = candidate_relations;
                std::size_t section_begin = candidate_relations.size();
                const std::size_t iota_rank = exact_rank(iota);
                result.negative_stripped_rank = iota_rank;
                result.positive_stripped_rank = exact_rank(positive_iota);
                const std::size_t allocation_dimension
                    = columns.size() - iota_rank;
                const std::size_t carrier_dimension
                    = iota_rank - result.rank;
                result.allocation_dimension = allocation_dimension;
                result.carrier_dimension = carrier_dimension;
                const auto record_section_criterion = [
                    &result,
                    &iota,
                    iota_rank,
                    allocation_dimension,
                    carrier_dimension
                ](const std::vector<std::vector<Integer>>& defects) {
                    const std::size_t defect_rank = exact_rank(defects);
                    auto combined = defects;
                    combined.insert(
                        combined.end(), iota.begin(), iota.end()
                    );
                    const std::size_t intersection = defect_rank + iota_rank
                        - exact_rank(std::move(combined));
                    const std::size_t defect_kernel
                        = defects.size() - defect_rank;
                    result.section_criteria.push_back(
                        Result::SectionCriterion{
                            defect_rank,
                            defect_kernel,
                            intersection,
                            intersection == carrier_dimension
                                && defect_kernel >= allocation_dimension
                        }
                    );
                };
                if (diagnostics) {
                std::vector<std::vector<Integer>> section(
                    basis_index.size(),
                    std::vector<Integer>(f0_basis.size(), 0)
                );
                for (const auto& [graph_key, row] : basis_index) {
                    section[row] = strip_to_f0(
                        unflatten(graph_key, degrees.size())
                    );
                }
                std::vector<std::vector<Integer>> positive_defects;
                positive_defects.reserve(
                    result.positive_source_graphs.size()
                );
                for (std::size_t positive = 0U;
                     positive < result.positive_source_graphs.size();
                     ++positive) {
                    std::vector<Integer> defect = positive_iota[positive];
                    const Key raw_key = flatten(
                        result.positive_source_graphs[positive].graph
                    );
                    for (const auto& [graph_key, coefficient]
                         : straighten(
                               raw_key, degrees.size(), straightening_cache
                           )) {
                        const auto found = basis_index.find(graph_key);
                        if (found == basis_index.end()) {
                            throw std::runtime_error(
                                "positive defect left the global basis"
                            );
                        }
                        for (std::size_t column = 0U;
                             column < f0_basis.size(); ++column) {
                            defect[column]
                                -= coefficient * section[found->second][column];
                        }
                    }
                    positive_defects.push_back(std::move(defect));
                }
                record_section_criterion(positive_defects);
                std::vector<std::vector<Integer>> defect_intersection_columns
                    = positive_defects;
                defect_intersection_columns.insert(
                    defect_intersection_columns.end(),
                    iota.begin(), iota.end()
                );
                const auto defect_intersection = exact_nullspace(
                    column_matrix(
                        defect_intersection_columns, f0_basis.size()
                    )
                );
                for (const std::vector<Rational>& solution
                     : defect_intersection) {
                    std::vector<Rational> relation(
                        columns.size(), Rational(0)
                    );
                    for (std::size_t source = 0U;
                         source < columns.size(); ++source) {
                        relation[source] = solution[
                            positive_defects.size() + source
                        ];
                    }
                    if (std::any_of(
                            relation.begin(), relation.end(),
                            [](const Rational& value) { return value != 0; }
                        )) {
                        candidate_relations.push_back(integerize(relation));
                        ++result.section_defect_relations;
                    }
                }
                result.section_trial_ranks.push_back(
                    exact_rank(candidate_relations)
                );
                {
                    auto individual = section_base_relations;
                    individual.insert(
                        individual.end(),
                        candidate_relations.begin()
                            + static_cast<std::ptrdiff_t>(section_begin),
                        candidate_relations.end()
                    );
                    result.section_individual_ranks.push_back(
                        exact_rank(std::move(individual))
                    );
                }
                section_begin = candidate_relations.size();
                std::vector<std::vector<Integer>> phi_rows(
                    f0_basis.size(),
                    std::vector<Integer>(basis_index.size(), 0)
                );
                for (std::size_t source = 0U;
                     source < f0_basis.size(); ++source) {
                    Graph product = unflatten(
                        f0_basis[source].second, degrees.size()
                    );
                    Integer orientation = 1;
                    add_oriented_edge(
                        product, fundamental, f0_basis[source].first,
                        orientation
                    );
                    const Key product_key = flatten(product);
                    for (const auto& [graph_key, coefficient]
                         : straighten(
                               product_key, degrees.size(),
                               straightening_cache
                           )) {
                        const auto found = basis_index.find(graph_key);
                        if (found == basis_index.end()) {
                            throw std::runtime_error(
                                "F0 multiplication left the global basis"
                            );
                        }
                        phi_rows[source][found->second]
                            += orientation * coefficient;
                    }
                }
                std::vector<std::vector<Rational>> gaussian_section(
                    basis_index.size(),
                    std::vector<Rational>(f0_basis.size(), Rational(0))
                );
                for (std::size_t row = 0U;
                     row < basis_index.size(); ++row) {
                    std::vector<Rational> target(
                        basis_index.size(), Rational(0)
                    );
                    target[row] = Rational(1);
                    const auto [solved, lift] = solve_left(phi_rows, target);
                    if (!solved) {
                        throw std::runtime_error("phi is not onto");
                    }
                    gaussian_section[row] = lift;
                }
                std::vector<std::vector<Integer>> gaussian_defects;
                gaussian_defects.reserve(
                    result.positive_source_graphs.size()
                );
                for (std::size_t positive = 0U;
                     positive < result.positive_source_graphs.size();
                     ++positive) {
                    std::vector<Rational> defect(
                        f0_basis.size(), Rational(0)
                    );
                    for (std::size_t column = 0U;
                         column < f0_basis.size(); ++column) {
                        defect[column] = Rational(
                            positive_iota[positive][column]
                        );
                    }
                    const Key raw_key = flatten(
                        result.positive_source_graphs[positive].graph
                    );
                    for (const auto& [graph_key, coefficient]
                         : straighten(
                               raw_key, degrees.size(), straightening_cache
                           )) {
                        const auto found = basis_index.find(graph_key);
                        if (found == basis_index.end()) {
                            throw std::runtime_error(
                                "Gaussian defect left the global basis"
                            );
                        }
                        for (std::size_t column = 0U;
                             column < f0_basis.size(); ++column) {
                            defect[column] -= Rational(coefficient)
                                * gaussian_section[found->second][column];
                        }
                    }
                    gaussian_defects.push_back(integerize(defect));
                }
                record_section_criterion(gaussian_defects);
                std::vector<std::vector<Integer>> gaussian_intersection_columns
                    = gaussian_defects;
                gaussian_intersection_columns.insert(
                    gaussian_intersection_columns.end(),
                    iota.begin(), iota.end()
                );
                const auto gaussian_intersection = exact_nullspace(
                    column_matrix(
                        gaussian_intersection_columns, f0_basis.size()
                    )
                );
                for (const std::vector<Rational>& solution
                     : gaussian_intersection) {
                    std::vector<Rational> relation(
                        columns.size(), Rational(0)
                    );
                    for (std::size_t source = 0U;
                         source < columns.size(); ++source) {
                        relation[source] = solution[
                            gaussian_defects.size() + source
                        ];
                    }
                    if (std::any_of(
                            relation.begin(), relation.end(),
                            [](const Rational& value) { return value != 0; }
                        )) {
                        candidate_relations.push_back(integerize(relation));
                        ++result.section_defect_relations;
                    }
                }
                result.section_trial_ranks.push_back(
                    exact_rank(candidate_relations)
                );
                {
                    auto individual = section_base_relations;
                    individual.insert(
                        individual.end(),
                        candidate_relations.begin()
                            + static_cast<std::ptrdiff_t>(section_begin),
                        candidate_relations.end()
                    );
                    result.section_individual_ranks.push_back(
                        exact_rank(std::move(individual))
                    );
                }
                section_begin = candidate_relations.size();
                const auto append_rational_section = [
                    &result,
                    &positive_iota,
                    &basis_index,
                    &straightening_cache,
                    &degrees,
                    &iota,
                    &columns,
                    &candidate_relations,
                    &record_section_criterion,
                    f0_dimension = f0_basis.size()
                ](const std::vector<std::vector<Rational>>& trial_section) {
                    std::vector<std::vector<Integer>> defects;
                    defects.reserve(result.positive_source_graphs.size());
                    for (std::size_t positive = 0U;
                         positive < result.positive_source_graphs.size();
                         ++positive) {
                        std::vector<Rational> defect(
                            f0_dimension, Rational(0)
                        );
                        for (std::size_t column = 0U;
                             column < f0_dimension; ++column) {
                            defect[column] = Rational(
                                positive_iota[positive][column]
                            );
                        }
                        const Key raw_key = flatten(
                            result.positive_source_graphs[positive].graph
                        );
                        for (const auto& [graph_key, coefficient]
                             : straighten(
                                   raw_key, degrees.size(),
                                   straightening_cache
                               )) {
                            const auto found = basis_index.find(graph_key);
                            if (found == basis_index.end()) {
                                throw std::runtime_error(
                                    "trial defect left the global basis"
                                );
                            }
                            for (std::size_t column = 0U;
                                 column < f0_dimension; ++column) {
                                defect[column] -= Rational(coefficient)
                                    * trial_section[found->second][column];
                            }
                        }
                        defects.push_back(integerize(defect));
                    }
                    record_section_criterion(defects);
                    std::vector<std::vector<Integer>> intersection_columns
                        = defects;
                    intersection_columns.insert(
                        intersection_columns.end(), iota.begin(), iota.end()
                    );
                    const auto intersection = exact_nullspace(
                        column_matrix(intersection_columns, f0_dimension)
                    );
                    for (const std::vector<Rational>& solution : intersection) {
                        std::vector<Rational> relation(
                            columns.size(), Rational(0)
                        );
                        for (std::size_t source = 0U;
                             source < columns.size(); ++source) {
                            relation[source]
                                = solution[defects.size() + source];
                        }
                        if (std::any_of(
                                relation.begin(), relation.end(),
                                [](const Rational& value) {
                                    return value != 0;
                                }
                            )) {
                            candidate_relations.push_back(
                                integerize(relation)
                            );
                            ++result.section_defect_relations;
                        }
                    }
                };
                if (!differential.empty()) {
                    for (std::size_t trial = 0U; trial < 4U; ++trial) {
                        auto trial_section = gaussian_section;
                        for (std::size_t row = 0U;
                             row < trial_section.size(); ++row) {
                            const std::size_t kernel_row
                                = (row * 131U + trial * 197U + 17U)
                                    % differential.size();
                            const int sign
                                = ((row + trial) & 1U) == 0U ? 1 : -1;
                            for (std::size_t column = 0U;
                                 column < f0_basis.size(); ++column) {
                                trial_section[row][column]
                                    += Rational(
                                        sign * differential[kernel_row][column]
                                    );
                            }
                        }
                        append_rational_section(trial_section);
                        result.section_trial_ranks.push_back(
                            exact_rank(candidate_relations)
                        );
                        auto individual = section_base_relations;
                        individual.insert(
                            individual.end(),
                            candidate_relations.begin()
                                + static_cast<std::ptrdiff_t>(section_begin),
                            candidate_relations.end()
                        );
                        result.section_individual_ranks.push_back(
                            exact_rank(std::move(individual))
                        );
                        section_begin = candidate_relations.size();
                    }
                }
                }
                result.candidate_relation_rank = exact_rank(
                    candidate_relations
                );
                std::vector<std::vector<
                    Result::FundamentalLiftTerm::Owner
                >> f1_owners(f1_basis.size());
                for (std::size_t source = 0U;
                     source < f1_basis.size(); ++source) {
                    for (std::size_t face = 0U; face < 3U; ++face) {
                        std::vector<Rational> face_target(
                            f0_basis.size(), Rational(0)
                        );
                        for (std::size_t column = 0U;
                             column < f0_basis.size(); ++column) {
                            face_target[column]
                                = Rational(faces[source][face][column]);
                        }
                        const auto [face_solved, face_lift]
                            = solve_left(positive_iota, face_target);
                        if (!face_solved) {
                            continue;
                        }
                        for (std::size_t positive = 0U;
                             positive < face_lift.size(); ++positive) {
                            if (face_lift[positive] != 0) {
                                f1_owners[source].push_back(
                                    Result::FundamentalLiftTerm::Owner{
                                        result.positive_source_graphs[
                                            positive
                                        ].mask,
                                        result.positive_source_graphs[
                                            positive
                                        ].basis,
                                        face,
                                        face_lift[positive]
                                    }
                                );
                            }
                        }
                    }
                }
                std::vector<std::size_t> owned_sources;
                std::vector<std::vector<Integer>> owned_differential;
                for (std::size_t source = 0U;
                     source < f1_basis.size(); ++source) {
                    if (!f1_owners[source].empty()) {
                        owned_sources.push_back(source);
                        owned_differential.push_back(differential[source]);
                    }
                }
                std::vector<std::vector<Rational>> kernel_targets;
                kernel_targets.reserve(result.nullspace.size());
                for (const std::vector<Rational>& relation
                     : result.nullspace) {
                    std::vector<Rational> target(
                        f0_basis.size(), Rational(0)
                    );
                    for (std::size_t source = 0U;
                         source < relation.size(); ++source) {
                        for (std::size_t column = 0U;
                             column < f0_basis.size(); ++column) {
                            target[column]
                                += relation[source] * iota[source][column];
                        }
                    }
                    kernel_targets.push_back(std::move(target));
                }
                std::vector<std::vector<Integer>> positive_global;
                positive_global.reserve(
                    result.positive_source_graphs.size()
                );
                for (const Result::PositiveSource& source
                     : result.positive_source_graphs) {
                    std::vector<Integer> row(basis_index.size(), 0);
                    const Key raw_key = flatten(source.graph);
                    for (const auto& [graph_key, coefficient]
                         : straighten(
                               raw_key, degrees.size(), straightening_cache
                           )) {
                        const auto found = basis_index.find(graph_key);
                        if (found == basis_index.end()) {
                            throw std::runtime_error(
                                "positive multiplication left global basis"
                            );
                        }
                        row[found->second] += coefficient;
                    }
                    positive_global.push_back(std::move(row));
                }
                const auto positive_global_kernel = exact_nullspace(
                    column_matrix(positive_global, basis_index.size())
                );
                result.positive_global_rank = exact_rank(positive_global);
                std::vector<std::vector<Integer>> positive_internal_rows;
                for (const std::vector<Rational>& relation
                     : positive_global_kernel) {
                    std::vector<Rational> target(
                        f0_basis.size(), Rational(0)
                    );
                    for (std::size_t positive = 0U;
                         positive < relation.size(); ++positive) {
                        for (std::size_t column = 0U;
                             column < f0_basis.size(); ++column) {
                            target[column] += relation[positive]
                                * Rational(positive_iota[positive][column]);
                        }
                    }
                    if (std::any_of(
                            target.begin(), target.end(),
                            [](const Rational& value) { return value != 0; }
                        )) {
                        positive_internal_rows.push_back(
                            integerize(target)
                        );
                    }
                }
                result.positive_internal_kernel_rank = exact_rank(
                    positive_internal_rows
                );
                std::vector<std::vector<Integer>> optimal_section_rows
                    = positive_internal_rows;
                for (const std::vector<Rational>& target : kernel_targets) {
                    optimal_section_rows.push_back(integerize(target));
                }
                result.optimal_section_rank = exact_rank(
                    std::move(optimal_section_rows)
                );
                result.carrier_overlap_dimension
                    = result.carrier_dimension
                        + result.positive_internal_kernel_rank
                        - result.optimal_section_rank;
                result.positive_allocation_dimension
                    = result.positive_sources
                        - result.positive_stripped_rank;
                result.unshared_negative_carrier_dimension
                    = result.carrier_dimension
                        - result.carrier_overlap_dimension;
                result.residual_negative_demand
                    = result.allocation_dimension
                        + result.unshared_negative_carrier_dimension;
                result.residual_positive_supply
                    = result.positive_allocation_dimension
                        + result.positive_global_rank;
                result.optimal_section_cover_feasible
                    = result.unshared_negative_carrier_dimension
                        <= result.positive_global_rank;
                result.optimal_section_criterion
                    = result.positive_sources >= result.optimal_section_rank
                    && result.positive_sources
                            - result.optimal_section_rank
                        >= allocation_dimension;
                if (result.optimal_section_criterion
                    != (result.residual_negative_demand
                        <= result.residual_positive_supply)) {
                    throw std::runtime_error(
                        "residual capacity identity failed"
                    );
                }
                if (result.f0_dimension < result.codomain_dimension) {
                    throw std::runtime_error("F0 is smaller than its image");
                }
                result.f0_kernel_dimension
                    = result.f0_dimension - result.codomain_dimension;
                std::vector<int> merged_degrees;
                merged_degrees.reserve(degrees.size() - 1U);
                for (std::size_t vertex = 0U;
                     vertex < degrees.size(); ++vertex) {
                    if (vertex == fundamental) {
                        continue;
                    }
                    merged_degrees.push_back(
                        vertex + 1U == degrees.size()
                            ? degrees[vertex] + 1 : degrees[vertex]
                    );
                }
                const std::size_t merged_target
                    = merged_degrees.size() - 1U;
                std::vector<F0Key> merged_f0_basis;
                std::map<F0Key, std::size_t> merged_f0_indices;
                for (std::size_t neighbor = 0U;
                     neighbor < merged_degrees.size(); ++neighbor) {
                    if (neighbor == merged_target
                        || merged_degrees[neighbor] == 0) {
                        continue;
                    }
                    std::vector<int> residual = merged_degrees;
                    --residual[merged_target];
                    --residual[neighbor];
                    for (const Key& top : tops(residual)) {
                        const F0Key key{
                            neighbor,
                            flatten(standard_graph(residual, top))
                        };
                        merged_f0_indices.emplace(
                            key, merged_f0_basis.size()
                        );
                        merged_f0_basis.push_back(key);
                    }
                }
                const std::size_t merged_f0_dimension
                    = merged_f0_basis.size();
                const std::size_t merged_codomain_dimension
                    = tops(merged_degrees).size();
                if (merged_f0_dimension < merged_codomain_dimension) {
                    throw std::runtime_error(
                        "merged row space is smaller than its image"
                    );
                }
                result.merged_row_kernel_dimension
                    = merged_f0_dimension - merged_codomain_dimension;
                result.deletion_contraction_valid
                    = result.merged_row_kernel_dimension
                        == result.f0_kernel_dimension;

                std::map<Key, Expansion> merged_straightening_cache;
                std::vector<std::vector<Integer>> projection(
                    f0_basis.size(),
                    std::vector<Integer>(merged_f0_dimension, 0)
                );
                const std::size_t target_vertex = degrees.size() - 1U;
                for (std::size_t source = 0U;
                     source < f0_basis.size(); ++source) {
                    const std::size_t neighbor = f0_basis[source].first;
                    if (neighbor == target_vertex) {
                        continue;
                    }
                    const std::size_t merged_neighbor
                        = neighbor < fundamental
                            ? neighbor : neighbor - 1U;
                    const Graph graph = delete_graph_vertex(
                        unflatten(
                            f0_basis[source].second, degrees.size()
                        ),
                        fundamental
                    );
                    const Key graph_key_input = flatten(graph);
                    for (const auto& [graph_key, coefficient]
                         : straighten(
                               graph_key_input,
                               merged_degrees.size(),
                               merged_straightening_cache
                           )) {
                        const auto found = merged_f0_indices.find(
                            F0Key{merged_neighbor, graph_key}
                        );
                        if (found == merged_f0_indices.end()) {
                            throw std::runtime_error(
                                "deletion projection left merged F0"
                            );
                        }
                        projection[source][found->second] += coefficient;
                    }
                }
                if (exact_rank(projection) != merged_f0_dimension) {
                    throw std::runtime_error(
                        "deletion projection is not onto merged F0"
                    );
                }

                const auto projected_kernel = apply_row_map(
                    differential, projection
                );
                const std::size_t differential_rank
                    = exact_rank(differential);
                const std::size_t projected_kernel_rank
                    = exact_rank(projected_kernel);
                if (differential_rank != result.f0_kernel_dimension
                    || projected_kernel_rank
                        != result.merged_row_kernel_dimension
                    || differential_rank != projected_kernel_rank) {
                    throw std::runtime_error(
                        "deletion projection failed kernel isomorphism"
                    );
                }
                const auto projected_iota = apply_row_map(iota, projection);
                const auto projected_positive_iota = apply_row_map(
                    positive_iota, projection
                );
                const std::size_t projected_negative_rank
                    = exact_rank(projected_iota);
                const std::size_t projected_positive_rank
                    = exact_rank(projected_positive_iota);
                const auto projected_negative_carrier = span_intersection(
                    projected_iota,
                    projected_kernel,
                    merged_f0_dimension
                );
                const auto projected_positive_carrier = span_intersection(
                    projected_positive_iota,
                    projected_kernel,
                    merged_f0_dimension
                );
                const std::size_t projected_negative_carrier_dimension
                    = exact_rank(projected_negative_carrier);
                const std::size_t projected_positive_carrier_dimension
                    = exact_rank(projected_positive_carrier);
                const auto projected_carrier_overlap = span_intersection(
                    projected_negative_carrier,
                    projected_positive_carrier,
                    merged_f0_dimension
                );
                const std::size_t projected_carrier_overlap_dimension
                    = exact_rank(projected_carrier_overlap);

                std::vector<std::vector<Integer>>
                    negative_carrier_rows;
                negative_carrier_rows.reserve(kernel_targets.size());
                for (const std::vector<Rational>& row : kernel_targets) {
                    if (std::any_of(
                            row.begin(), row.end(),
                            [](const Rational& value) {
                                return value != 0;
                            }
                        )) {
                        negative_carrier_rows.push_back(integerize(row));
                    }
                }
                const auto projected_realized_negative_carrier
                    = apply_row_map(negative_carrier_rows, projection);
                const auto projected_realized_positive_carrier
                    = apply_row_map(positive_internal_rows, projection);
                if (exact_rank(projected_realized_negative_carrier)
                        != result.carrier_dimension
                    || exact_rank(projected_realized_positive_carrier)
                        != result.positive_internal_kernel_rank) {
                    throw std::runtime_error(
                        "realized carrier lost rank under deletion"
                    );
                }
                auto negative_containment = projected_negative_carrier;
                negative_containment.insert(
                    negative_containment.end(),
                    projected_realized_negative_carrier.begin(),
                    projected_realized_negative_carrier.end()
                );
                auto positive_containment = projected_positive_carrier;
                positive_containment.insert(
                    positive_containment.end(),
                    projected_realized_positive_carrier.begin(),
                    projected_realized_positive_carrier.end()
                );
                if (exact_rank(negative_containment)
                        != projected_negative_carrier_dimension
                    || exact_rank(positive_containment)
                        != projected_positive_carrier_dimension) {
                    throw std::runtime_error(
                        "realized carrier left projected carrier space"
                    );
                }

                result.negative_lower_branch_dimension
                    = iota_rank - projected_negative_rank;
                result.positive_lower_branch_dimension
                    = result.positive_stripped_rank
                        - projected_positive_rank;
                result.projected_allocation_dimension
                    = columns.size() - projected_negative_rank;
                result.projected_positive_allocation_dimension
                    = result.positive_source_graphs.size()
                        - projected_positive_rank;
                result.projected_unshared_carrier_dimension
                    = projected_negative_carrier_dimension
                        - projected_carrier_overlap_dimension;
                result.projected_positive_rank
                    = projected_positive_rank
                        - projected_positive_carrier_dimension;
                if (projected_negative_carrier_dimension
                        < result.carrier_dimension
                    || projected_positive_carrier_dimension
                        < result.positive_internal_kernel_rank) {
                    throw std::runtime_error(
                        "negative lift-obstruction rank"
                    );
                }
                result.negative_lift_obstruction_rank
                    = projected_negative_carrier_dimension
                        - result.carrier_dimension;
                result.positive_lift_obstruction_rank
                    = projected_positive_carrier_dimension
                        - result.positive_internal_kernel_rank;

                result.enriched_deletion_contraction_valid
                    = result.allocation_dimension
                        + result.negative_lower_branch_dimension
                            == result.projected_allocation_dimension
                    && result.positive_allocation_dimension
                        + result.positive_lower_branch_dimension
                            == result.projected_positive_allocation_dimension
                    && result.positive_global_rank
                        == result.projected_positive_rank
                            + result.positive_lower_branch_dimension
                            + result.positive_lift_obstruction_rank
                    && result.unshared_negative_carrier_dimension
                        <= result.projected_unshared_carrier_dimension
                            + result.positive_lift_obstruction_rank;
                if (!result.enriched_deletion_contraction_valid) {
                    throw std::runtime_error(
                        "enriched deletion-contraction recurrence failed"
                    );
                }
                result.projected_carrier_criterion
                    = result.projected_unshared_carrier_dimension
                            <= result.projected_positive_rank
                        && result.projected_allocation_dimension
                                + result.projected_unshared_carrier_dimension
                            <= result.projected_positive_allocation_dimension
                                + result.projected_positive_rank;
                result.ambient_kernel_criterion
                    = result.positive_sources
                            >= result.f0_kernel_dimension
                    && result.positive_sources
                            - result.f0_kernel_dimension
                        >= allocation_dimension;

                if (result.optimal_section_cover_feasible) {
                    std::vector<std::vector<Integer>> domain_basis;
                    std::vector<std::vector<Rational>> desired_outputs;
                    domain_basis.reserve(result.positive_source_graphs.size());
                    desired_outputs.reserve(
                        result.positive_source_graphs.size()
                    );
                    for (const std::vector<Rational>& relation
                         : positive_global_kernel) {
                        std::vector<Integer> integer_relation
                            = integerize(relation);
                        std::vector<Rational> target(
                            f0_basis.size(), Rational(0)
                        );
                        for (std::size_t positive = 0U;
                             positive < integer_relation.size(); ++positive) {
                            for (std::size_t column = 0U;
                                 column < f0_basis.size(); ++column) {
                                target[column]
                                    += Rational(integer_relation[positive])
                                        * Rational(
                                            positive_iota[positive][column]
                                        );
                            }
                        }
                        domain_basis.push_back(
                            std::move(integer_relation)
                        );
                        desired_outputs.push_back(std::move(target));
                    }
                    std::vector<std::vector<Integer>> independent_mu_rows;
                    std::vector<std::size_t> complement_sources;
                    std::size_t mu_rank = 0U;
                    for (std::size_t positive = 0U;
                         positive < positive_global.size(); ++positive) {
                        auto trial = independent_mu_rows;
                        trial.push_back(positive_global[positive]);
                        const std::size_t next_rank = exact_rank(trial);
                        if (next_rank == mu_rank) {
                            continue;
                        }
                        independent_mu_rows.push_back(
                            positive_global[positive]
                        );
                        complement_sources.push_back(positive);
                        mu_rank = next_rank;
                    }
                    if (mu_rank != result.positive_global_rank) {
                        throw std::runtime_error(
                            "positive complement rank mismatch"
                        );
                    }
                    std::vector<std::vector<Integer>> carrier_span
                        = positive_internal_rows;
                    std::vector<std::vector<Rational>> carrier_complement;
                    std::size_t carrier_span_rank = exact_rank(carrier_span);
                    for (const std::vector<Rational>& target
                         : kernel_targets) {
                        std::vector<Integer> integer_target
                            = integerize(target);
                        auto trial = carrier_span;
                        trial.push_back(integer_target);
                        const std::size_t next_rank = exact_rank(trial);
                        if (next_rank == carrier_span_rank) {
                            continue;
                        }
                        carrier_span.push_back(std::move(integer_target));
                        carrier_complement.push_back(target);
                        carrier_span_rank = next_rank;
                    }
                    if (carrier_complement.size()
                        != result.unshared_negative_carrier_dimension) {
                        throw std::runtime_error(
                            "negative carrier quotient rank mismatch"
                        );
                    }
                    for (std::size_t complement = 0U;
                         complement < complement_sources.size();
                         ++complement) {
                        std::vector<Integer> coordinate(
                            result.positive_source_graphs.size(), 0
                        );
                        coordinate[complement_sources[complement]] = 1;
                        domain_basis.push_back(std::move(coordinate));
                        if (complement < carrier_complement.size()) {
                            desired_outputs.push_back(
                                carrier_complement[complement]
                            );
                        } else {
                            desired_outputs.emplace_back(
                                f0_basis.size(), Rational(0)
                            );
                        }
                    }
                    if (domain_basis.size()
                            != result.positive_source_graphs.size()
                        || desired_outputs.size() != domain_basis.size()) {
                        throw std::runtime_error(
                            "optimal defect basis is incomplete"
                        );
                    }
                    std::vector<std::vector<Integer>> optimal_defects;
                    optimal_defects.reserve(domain_basis.size());
                    for (std::size_t positive = 0U;
                         positive < domain_basis.size(); ++positive) {
                        std::vector<Rational> unit(
                            domain_basis.size(), Rational(0)
                        );
                        unit[positive] = Rational(1);
                        const auto [solved, coordinates]
                            = solve_left(domain_basis, unit);
                        if (!solved) {
                            throw std::runtime_error(
                                "optimal defect domain basis is singular"
                            );
                        }
                        std::vector<Rational> target(
                            f0_basis.size(), Rational(0)
                        );
                        for (std::size_t basis = 0U;
                             basis < coordinates.size(); ++basis) {
                            for (std::size_t column = 0U;
                                 column < f0_basis.size(); ++column) {
                                target[column] += coordinates[basis]
                                    * desired_outputs[basis][column];
                            }
                        }
                        optimal_defects.push_back(integerize(target));
                    }
                    result.optimal_defect_rank = exact_rank(optimal_defects);
                    result.optimal_defect_kernel
                        = optimal_defects.size()
                            - result.optimal_defect_rank;
                    if (result.optimal_defect_rank
                            != result.optimal_section_rank
                        || result.optimal_defect_kernel
                            < result.allocation_dimension) {
                        throw std::runtime_error(
                            "constructed optimal defect has wrong capacity"
                        );
                    }

                    const auto allocation_relations = exact_nullspace(
                        column_matrix(iota, f0_basis.size())
                    );
                    for (const std::vector<Rational>& relation
                         : allocation_relations) {
                        candidate_relations.push_back(
                            integerize(relation)
                        );
                    }
                    std::vector<std::vector<Integer>> intersection_columns
                        = optimal_defects;
                    intersection_columns.insert(
                        intersection_columns.end(), iota.begin(), iota.end()
                    );
                    const auto optimal_intersection = exact_nullspace(
                        column_matrix(
                            intersection_columns, f0_basis.size()
                        )
                    );
                    for (const std::vector<Rational>& solution
                         : optimal_intersection) {
                        std::vector<Rational> relation(
                            columns.size(), Rational(0)
                        );
                        for (std::size_t source = 0U;
                             source < columns.size(); ++source) {
                            relation[source] = solution[
                                optimal_defects.size() + source
                            ];
                        }
                        if (std::any_of(
                                relation.begin(), relation.end(),
                                [](const Rational& value) {
                                    return value != 0;
                                }
                            )) {
                            candidate_relations.push_back(
                                integerize(relation)
                            );
                        }
                    }
                    result.optimal_relation_rank = exact_rank(
                        candidate_relations
                    );
                    for (const std::vector<Integer>& relation
                         : candidate_relations) {
                        for (std::size_t row = 0U;
                             row < full_matrix.size(); ++row) {
                            Integer value = 0;
                            for (std::size_t source = 0U;
                                 source < relation.size(); ++source) {
                                value += full_matrix[row][source]
                                    * relation[source];
                            }
                            if (value != 0) {
                                throw std::runtime_error(
                                    "optimal defect produced invalid relation"
                                );
                            }
                        }
                    }
                    result.candidate_relation_rank
                        = result.optimal_relation_rank;
                }
                const std::size_t owned_rank = exact_rank(
                    owned_differential
                );
                std::vector<std::vector<Integer>> owned_and_kernel
                    = owned_differential;
                for (const std::vector<Rational>& target : kernel_targets) {
                    owned_and_kernel.push_back(integerize(target));
                }
                result.fundamental_kernel_quotient_rank
                    = exact_rank(std::move(owned_and_kernel)) - owned_rank;
                std::vector<std::vector<Integer>> owned_and_candidates
                    = owned_differential;
                for (const std::vector<Integer>& relation
                     : candidate_relations) {
                    std::vector<Integer> target(f0_basis.size(), 0);
                    for (std::size_t source = 0U;
                         source < relation.size(); ++source) {
                        for (std::size_t column = 0U;
                             column < f0_basis.size(); ++column) {
                            target[column]
                                += relation[source] * iota[source][column];
                        }
                    }
                    owned_and_candidates.push_back(std::move(target));
                }
                result.fundamental_candidate_quotient_rank
                    = exact_rank(std::move(owned_and_candidates)) - owned_rank;
                std::vector<std::vector<Integer>> kernel_integer_relations;
                std::vector<std::vector<Integer>> carrier_columns;
                for (const std::vector<Rational>& relation
                     : result.nullspace) {
                    std::vector<Integer> integer_relation
                        = integerize(relation);
                    std::vector<Integer> target(f0_basis.size(), 0);
                    for (std::size_t source = 0U;
                         source < integer_relation.size(); ++source) {
                        for (std::size_t column = 0U;
                             column < f0_basis.size(); ++column) {
                            target[column] += integer_relation[source]
                                * iota[source][column];
                        }
                    }
                    kernel_integer_relations.push_back(
                        std::move(integer_relation)
                    );
                    carrier_columns.push_back(std::move(target));
                }
                carrier_columns.insert(
                    carrier_columns.end(), owned_differential.begin(),
                    owned_differential.end()
                );
                std::vector<std::vector<Integer>> covered_relations
                    = candidate_relations;
                const auto relative_solutions = exact_nullspace(
                    column_matrix(carrier_columns, f0_basis.size())
                );
                for (const std::vector<Rational>& solution
                     : relative_solutions) {
                    std::vector<Rational> relation(
                        columns.size(), Rational(0)
                    );
                    for (std::size_t kernel = 0U;
                         kernel < kernel_integer_relations.size(); ++kernel) {
                        for (std::size_t source = 0U;
                             source < columns.size(); ++source) {
                            relation[source] += solution[kernel]
                                * Rational(
                                    kernel_integer_relations[kernel][source]
                                );
                        }
                    }
                    if (std::any_of(
                            relation.begin(), relation.end(),
                            [](const Rational& value) { return value != 0; }
                        )) {
                        covered_relations.push_back(integerize(relation));
                    }
                }
                result.fundamental_covered_rank = exact_rank(
                    covered_relations
                );
                result.fundamental_combined_valid
                    = result.fundamental_covered_rank
                        == result.nullspace.size();
                if (!result.fundamental_combined_valid) {
                    std::size_t covered_rank = result.fundamental_covered_rank;
                    for (const std::vector<Integer>& relation
                         : kernel_integer_relations) {
                        auto augmented = covered_relations;
                        augmented.push_back(relation);
                        const std::size_t augmented_rank
                            = exact_rank(std::move(augmented));
                        if (augmented_rank > covered_rank) {
                            result.fundamental_missing_relation = relation;
                            break;
                        }
                    }
                }
                for (const std::vector<Rational>& target : kernel_targets) {
                    const auto [owned_solved, owned_lift]
                        = solve_left(owned_differential, target);
                    std::vector<Rational> lift(
                        f1_basis.size(), Rational(0)
                    );
                    bool solved = owned_solved;
                    if (owned_solved) {
                        for (std::size_t source = 0U;
                             source < owned_sources.size(); ++source) {
                            lift[owned_sources[source]] = owned_lift[source];
                        }
                    } else {
                        result.fundamental_owned_lifts_valid = false;
                        const auto fallback = solve_left(
                            differential, target
                        );
                        solved = fallback.first;
                        if (solved) {
                            lift = fallback.second;
                        }
                    }
                    if (!solved) {
                        result.fundamental_lifts_valid = false;
                        result.fundamental_lifts.emplace_back();
                        continue;
                    }
                    std::vector<Result::FundamentalLiftTerm> terms;
                    for (std::size_t source = 0U;
                         source < lift.size(); ++source) {
                        if (lift[source] == 0) {
                            continue;
                        }
                        Result::FundamentalLiftTerm term;
                        term.triple = f1_basis[source].triple;
                        term.coefficient = f1_basis[source].coefficient;
                        term.scalar = lift[source];
                        term.owners = f1_owners[source];
                        terms.push_back(std::move(term));
                    }
                    result.fundamental_lifts.push_back(std::move(terms));
                }
            }
        }
        if (diagnostics) {
            for (auto first = odd_columns.begin();
                 first != odd_columns.end(); ++first) {
                for (auto second = std::next(first);
                     second != odd_columns.end(); ++second) {
                    std::vector<std::vector<Integer>> combined = first->second;
                    combined.insert(
                        combined.end(), second->second.begin(),
                        second->second.end()
                    );
                    const std::size_t union_rank = exact_rank(
                        column_matrix(combined, basis_index.size())
                    );
                    const std::size_t intersection = first->second.size()
                        + second->second.size() - union_rank;
                    result.pair_intersections.push_back(
                        Result::PairIntersection{
                            first->first,
                            second->first,
                            intersection,
                            channel_dimensions[first->first ^ second->first]
                        }
                    );
                }
            }
        }
    }
    result.straightening_states = straightening_cache.size();
    return result;
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

void print_vector(const std::vector<int>& values) {
    std::cout << '[';
    for (std::size_t i = 0U; i < values.size(); ++i) {
        if (i != 0U) {
            std::cout << ',';
        }
        std::cout << values[i];
    }
    std::cout << ']';
}

void print_graph(const Graph& graph) {
    bool first = true;
    std::cout << '{';
    for (std::size_t left = 0U; left < graph.size(); ++left) {
        for (std::size_t right = left + 1U;
             right < graph.size(); ++right) {
            if (graph[left][right] == 0) {
                continue;
            }
            if (!first) {
                std::cout << ',';
            }
            first = false;
            std::cout << left << '-' << right;
            if (graph[left][right] != 1) {
                std::cout << '^' << graph[left][right];
            }
        }
    }
    std::cout << '}';
}

bool increment_labels(std::vector<int>& labels, int maximum_label) {
    for (std::size_t offset = 0U; offset < labels.size(); ++offset) {
        const std::size_t position = labels.size() - 1U - offset;
        if (labels[position] >= maximum_label) {
            continue;
        }
        const int next = labels[position] + 1;
        for (std::size_t i = position; i < labels.size(); ++i) {
            labels[i] = next;
        }
        return true;
    }
    return false;
}

bool support_disjoint(
    const std::vector<int>& labels,
    unsigned int minus_mask
) {
    for (std::size_t first = 0U; first < labels.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < labels.size(); ++second) {
            if (labels[first] == labels[second]
                && ((minus_mask >> first) & 1U)
                    != ((minus_mask >> second) & 1U)) {
                return false;
            }
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 6
            && (std::string(argv[1]) == "sweep"
                || std::string(argv[1]) == "relation-sweep"
                || std::string(argv[1]) == "fundamental-sweep")) {
            const bool relation_sweep
                = std::string(argv[1]) == "relation-sweep";
            const bool fundamental_sweep
                = std::string(argv[1]) == "fundamental-sweep";
            const bool capacity_verbose
                = std::getenv("SU2_PLUCKER_CAPACITY_VERBOSE") != nullptr;
            const int maximum_label
                = parse_nonnegative(argv[2], "maximum label");
            const int maximum_factors
                = parse_nonnegative(argv[3], "maximum factors");
            const int maximum_target
                = parse_nonnegative(argv[4], "maximum target");
            const int maximum_total
                = parse_nonnegative(argv[5], "maximum total");
            if (maximum_label <= 0 || maximum_factors <= 0
                || maximum_factors
                    >= static_cast<int>(
                        std::numeric_limits<unsigned int>::digits
                    )) {
                throw std::runtime_error("sweep bounds are out of range");
            }
            std::uint64_t cases = 0U;
            std::uint64_t deficit_cases = 0U;
            std::uint64_t nonzero_kernels = 0U;
            std::uint64_t maximum_kernel = 0U;
            std::uint64_t minimum_slack
                = std::numeric_limits<std::uint64_t>::max();
            std::uint64_t minimum_residual_slack
                = std::numeric_limits<std::uint64_t>::max();
            std::uint64_t minimum_cover_slack
                = std::numeric_limits<std::uint64_t>::max();
            std::uint64_t projected_cases = 0U;
            std::uint64_t projected_failures = 0U;
            std::uint64_t printed_projected_failures = 0U;
            std::uint64_t minimum_projected_cover_slack
                = std::numeric_limits<std::uint64_t>::max();
            std::uint64_t minimum_projected_capacity_slack
                = std::numeric_limits<std::uint64_t>::max();
            std::uint64_t maximum_projected_cover_deficit = 0U;
            std::uint64_t maximum_projected_capacity_deficit = 0U;
            std::uint64_t printed_nonzero = 0U;
            std::uint64_t printed_tight = 0U;
            for (int factors = 1; factors <= maximum_factors; ++factors) {
                std::vector<int> labels(
                    static_cast<std::size_t>(factors), 1
                );
                bool more = true;
                while (more) {
                    int label_total = 0;
                    for (int label : labels) {
                        label_total += label;
                    }
                    for (int target = 0; target <= maximum_target; ++target) {
                        if (label_total + target > maximum_total
                            || ((label_total + target) & 1) != 0) {
                            continue;
                        }
                        const unsigned int mask_limit = 1U << factors;
                        for (unsigned int minus_mask = 1U;
                             minus_mask < mask_limit; ++minus_mask) {
                            if (!support_disjoint(labels, minus_mask)) {
                                continue;
                            }
                            if (fundamental_sweep) {
                                std::size_t negative = 0U;
                                while (((minus_mask >> negative) & 1U) == 0U) {
                                    ++negative;
                                }
                                if (popcount(minus_mask) != 1
                                    || labels[negative] != 1) {
                                    continue;
                                }
                            }
                            const Result result = analyze_case(
                                labels, target, minus_mask, false,
                                relation_sweep, fundamental_sweep
                            );
                            ++cases;
                            if (fundamental_sweep && target > 0) {
                                ++projected_cases;
                                if (!result.enriched_deletion_contraction_valid) {
                                    std::cout
                                        << "SU2_PLUCKER_PROJECTED_CARRIER "
                                        << "result=FAIL labels=";
                                    print_vector(labels);
                                    std::cout
                                        << " target=" << target
                                        << " minus_mask=" << minus_mask
                                        << " hat_alpha="
                                        << result.projected_allocation_dimension
                                        << " hat_beta="
                                        << result.projected_positive_allocation_dimension
                                        << " hat_c="
                                        << result.projected_unshared_carrier_dimension
                                        << " bar_r="
                                        << result.projected_positive_rank
                                        << " A0="
                                        << result.negative_lower_branch_dimension
                                        << " B0="
                                        << result.positive_lower_branch_dimension
                                        << " theta_minus="
                                        << result.negative_lift_obstruction_rank
                                        << " theta_plus="
                                        << result.positive_lift_obstruction_rank
                                        << " recurrences="
                                        << (result.enriched_deletion_contraction_valid
                                                ? "PASS" : "FAIL")
                                        << " inequalities="
                                        << (result.projected_carrier_criterion
                                                ? "PASS" : "FAIL")
                                        << '\n';
                                    return 1;
                                }
                                const std::size_t projected_demand
                                    = result.projected_unshared_carrier_dimension;
                                const std::size_t projected_supply
                                    = result.projected_positive_rank;
                                const std::size_t projected_capacity_demand
                                    = result.projected_allocation_dimension
                                        + projected_demand;
                                const std::size_t projected_capacity_supply
                                    = result.projected_positive_allocation_dimension
                                        + projected_supply;
                                if (projected_demand <= projected_supply
                                    && projected_capacity_demand
                                        <= projected_capacity_supply) {
                                    minimum_projected_cover_slack = std::min(
                                        minimum_projected_cover_slack,
                                        static_cast<std::uint64_t>(
                                            projected_supply
                                                - projected_demand
                                        )
                                    );
                                    minimum_projected_capacity_slack = std::min(
                                        minimum_projected_capacity_slack,
                                        static_cast<std::uint64_t>(
                                            projected_capacity_supply
                                                - projected_capacity_demand
                                        )
                                    );
                                } else {
                                    ++projected_failures;
                                    if (projected_demand > projected_supply) {
                                        maximum_projected_cover_deficit
                                            = std::max(
                                                maximum_projected_cover_deficit,
                                                static_cast<std::uint64_t>(
                                                    projected_demand
                                                        - projected_supply
                                                )
                                            );
                                    }
                                    if (projected_capacity_demand
                                            > projected_capacity_supply) {
                                        maximum_projected_capacity_deficit
                                            = std::max(
                                                maximum_projected_capacity_deficit,
                                                static_cast<std::uint64_t>(
                                                    projected_capacity_demand
                                                        - projected_capacity_supply
                                                )
                                            );
                                    }
                                    if (printed_projected_failures < 12U) {
                                        std::cout
                                            << "projected_failure labels=";
                                        print_vector(labels);
                                        std::cout
                                            << " target=" << target
                                            << " minus_mask=" << minus_mask
                                            << " hat_alpha="
                                            << result.projected_allocation_dimension
                                            << " hat_beta="
                                            << result.projected_positive_allocation_dimension
                                            << " hat_c="
                                            << result.projected_unshared_carrier_dimension
                                            << " bar_r="
                                            << result.projected_positive_rank
                                            << " A0="
                                            << result.negative_lower_branch_dimension
                                            << " B0="
                                            << result.positive_lower_branch_dimension
                                            << " theta_minus="
                                            << result.negative_lift_obstruction_rank
                                            << " theta_plus="
                                            << result.positive_lift_obstruction_rank
                                            << '\n';
                                        ++printed_projected_failures;
                                    }
                                }
                            }
                            if (result.odd_sources
                                <= result.positive_sources) {
                                continue;
                            }
                            ++deficit_cases;
                            const std::uint64_t kernel
                                = result.odd_sources
                                    - static_cast<std::uint64_t>(result.rank);
                            if (fundamental_sweep
                                && (!result.fundamental_lifts_valid
                                    || !result.fundamental_combined_valid
                                    || !result.optimal_section_criterion
                                    || !result.optimal_section_cover_feasible
                                    || !result.deletion_contraction_valid
                                    || !result.enriched_deletion_contraction_valid)) {
                                std::cout
                                    << "SU2_PLUCKER_FUNDAMENTAL_LIFT "
                                    << "result=FAIL labels=";
                                print_vector(labels);
                                std::cout << " target=" << target
                                          << " minus_mask=" << minus_mask
                                          << " kernel=" << kernel
                                          << " lift="
                                          << (result.fundamental_lifts_valid
                                                  ? "PASS" : "FAIL")
                                          << " owned="
                                          << (result.fundamental_owned_lifts_valid
                                                  ? "PASS" : "FAIL")
                                          << " candidate_rank="
                                          << result.candidate_relation_rank
                                          << " kernel_quotient_rank="
                                          << result.fundamental_kernel_quotient_rank
                                          << " candidate_quotient_rank="
                                          << result.fundamental_candidate_quotient_rank
                                          << " covered_rank="
                                          << result.fundamental_covered_rank
                                          << " allocation="
                                          << result.allocation_dimension
                                          << " carrier="
                                          << result.carrier_dimension
                                          << " optimal_section_rank="
                                          << result.optimal_section_rank
                                          << " optimal="
                                          << (result.optimal_section_criterion
                                                  ? "PASS" : "FAIL")
                                          << " cover="
                                          << (result.optimal_section_cover_feasible
                                                  ? "PASS" : "FAIL")
                                          << " optimal_defect_rank="
                                          << result.optimal_defect_rank
                                          << " optimal_defect_kernel="
                                          << result.optimal_defect_kernel
                                          << " optimal_relation_rank="
                                          << result.optimal_relation_rank
                                          << " merged_row_kernel="
                                          << result.merged_row_kernel_dimension
                                          << " deletion_contraction="
                                          << (result.deletion_contraction_valid
                                                  ? "PASS" : "FAIL")
                                          << " projected="
                                          << (result.projected_carrier_criterion
                                                  ? "PASS" : "FAIL")
                                          << '\n';
                                return 1;
                            }
                            if (fundamental_sweep) {
                                minimum_residual_slack = std::min(
                                    minimum_residual_slack,
                                    static_cast<std::uint64_t>(
                                        result.residual_positive_supply
                                            - result.residual_negative_demand
                                    )
                                );
                                minimum_cover_slack = std::min(
                                    minimum_cover_slack,
                                    static_cast<std::uint64_t>(
                                        result.positive_global_rank
                                            - result.unshared_negative_carrier_dimension
                                    )
                                );
                            }
                            if (fundamental_sweep && capacity_verbose) {
                                std::cout << "capacity labels=";
                                print_vector(labels);
                                std::cout << " target=" << target
                                          << " odd=" << result.odd_sources
                                          << " positive="
                                          << result.positive_sources
                                          << " allocation="
                                          << result.allocation_dimension
                                          << " negative_carrier="
                                          << result.carrier_dimension
                                          << " positive_carrier="
                                          << result.positive_internal_kernel_rank
                                          << " carrier_sum="
                                          << result.optimal_section_rank
                                          << " carrier_overlap="
                                          << result.carrier_overlap_dimension
                                          << " f0=" << result.f0_dimension
                                          << " ker_f0="
                                          << result.f0_kernel_dimension
                                          << " rank_iota_minus="
                                          << result.negative_stripped_rank
                                          << " rank_iota_plus="
                                          << result.positive_stripped_rank
                                          << " rank_mu_plus="
                                          << result.positive_global_rank
                                          << " positive_allocation="
                                          << result.positive_allocation_dimension
                                          << " unshared_carrier="
                                          << result.unshared_negative_carrier_dimension
                                          << " residual_demand="
                                          << result.residual_negative_demand
                                          << " residual_supply="
                                          << result.residual_positive_supply
                                          << " hat_alpha="
                                          << result.projected_allocation_dimension
                                          << " hat_beta="
                                          << result.projected_positive_allocation_dimension
                                          << " hat_c="
                                          << result.projected_unshared_carrier_dimension
                                          << " bar_r="
                                          << result.projected_positive_rank
                                          << " A0="
                                          << result.negative_lower_branch_dimension
                                          << " B0="
                                          << result.positive_lower_branch_dimension
                                          << " theta_minus="
                                          << result.negative_lift_obstruction_rank
                                          << " theta_plus="
                                          << result.positive_lift_obstruction_rank
                                          << " ambient="
                                          << (result.ambient_kernel_criterion
                                                  ? "PASS" : "FAIL")
                                          << '\n';
                            }
                            if (relation_sweep
                                && (!result.candidate_relations_valid
                                    || result.candidate_relation_rank
                                        != kernel)) {
                                std::cout
                                    << "SU2_PLUCKER_RELATION_KERNEL "
                                    << "result=FAIL labels=";
                                print_vector(labels);
                                std::cout << " target=" << target
                                          << " minus_mask=" << minus_mask
                                          << " kernel=" << kernel
                                          << " candidate_rank="
                                          << result.candidate_relation_rank
                                          << " raw_relations="
                                          << result.raw_collision_relations
                                          << " pair_relations="
                                          << result.pair_intersection_relations
                                          << " edge_relations="
                                          << result.edge_exchange_relations
                                          << " cut_boundary_relations="
                                          << result.cut_boundary_relations
                                          << " cut_boundary_rank="
                                          << result.cut_boundary_rank
                                          << " section_defect_relations="
                                          << result.section_defect_relations
                                          << " valid="
                                          << (result.candidate_relations_valid
                                                  ? "PASS" : "FAIL")
                                          << '\n';
                                return 1;
                            }
                            if (kernel > result.positive_sources) {
                                std::cout << "SU2_PLUCKER_ORIGINAL_KERNEL "
                                          << "result=FAIL labels=";
                                print_vector(labels);
                                std::cout << " target=" << target
                                          << " minus_mask=" << minus_mask
                                          << " odd=" << result.odd_sources
                                          << " rank=" << result.rank
                                          << " kernel=" << kernel
                                          << " positive="
                                          << result.positive_sources << '\n';
                                return 1;
                            }
                            if (kernel != 0U) {
                                ++nonzero_kernels;
                                if (printed_nonzero < 12U) {
                                    std::cout << "nonzero_kernel labels=";
                                    print_vector(labels);
                                    std::cout << " target=" << target
                                              << " minus_mask=" << minus_mask
                                              << " codomain="
                                              << result.codomain_dimension
                                              << " odd=" << result.odd_sources
                                              << " rank=" << result.rank
                                              << " kernel=" << kernel
                                              << " positive="
                                              << result.positive_sources
                                              << '\n';
                                    ++printed_nonzero;
                                }
                            }
                            if (kernel != 0U
                                && kernel == result.positive_sources
                                && printed_tight < 12U) {
                                std::cout << "tight_kernel labels=";
                                print_vector(labels);
                                std::cout << " target=" << target
                                          << " minus_mask=" << minus_mask
                                          << " codomain="
                                          << result.codomain_dimension
                                          << " odd=" << result.odd_sources
                                          << " rank=" << result.rank
                                          << " kernel=" << kernel
                                          << " positive="
                                          << result.positive_sources
                                          << '\n';
                                ++printed_tight;
                            }
                            maximum_kernel = std::max(maximum_kernel, kernel);
                            minimum_slack = std::min(
                                minimum_slack,
                                result.positive_sources - kernel
                            );
                        }
                    }
                    more = increment_labels(labels, maximum_label);
                }
                std::cout << "progress factors=" << factors
                          << " cases=" << cases
                          << " deficit_cases=" << deficit_cases
                          << " nonzero_kernels=" << nonzero_kernels
                          << '\n' << std::flush;
            }
            std::cout << "SU2_PLUCKER_ORIGINAL_KERNEL cases=" << cases
                      << " deficit_cases=" << deficit_cases
                      << " nonzero_kernels=" << nonzero_kernels
                      << " maximum_kernel=" << maximum_kernel
                      << " minimum_slack="
                      << (deficit_cases == 0U ? 0U : minimum_slack)
                      << " result=PASS\n";
            if (relation_sweep) {
                std::cout << "SU2_PLUCKER_RELATION_KERNEL cases="
                          << cases << " deficit_cases=" << deficit_cases
                          << " result=PASS\n";
            }
            if (fundamental_sweep) {
                std::cout << "SU2_PLUCKER_FUNDAMENTAL_LIFT cases="
                          << cases << " deficit_cases=" << deficit_cases
                          << " minimum_residual_slack="
                          << (minimum_residual_slack
                                      == std::numeric_limits<
                                          std::uint64_t>::max()
                                  ? 0U : minimum_residual_slack)
                          << " minimum_cover_slack="
                          << (minimum_cover_slack
                                      == std::numeric_limits<
                                          std::uint64_t>::max()
                                  ? 0U : minimum_cover_slack)
                          << " result=PASS\n";
                std::cout << "SU2_PLUCKER_PROJECTED_CARRIER cases="
                          << projected_cases
                          << " failures=" << projected_failures
                          << " minimum_cover_slack="
                          << (minimum_projected_cover_slack
                                      == std::numeric_limits<
                                          std::uint64_t>::max()
                                  ? 0U : minimum_projected_cover_slack)
                          << " minimum_capacity_slack="
                          << (minimum_projected_capacity_slack
                                      == std::numeric_limits<
                                          std::uint64_t>::max()
                                  ? 0U : minimum_projected_capacity_slack)
                          << " maximum_cover_deficit="
                          << maximum_projected_cover_deficit
                          << " maximum_capacity_deficit="
                          << maximum_projected_capacity_deficit
                          << " recurrences=PASS inequalities="
                          << (projected_failures == 0U ? "PASS" : "FAIL")
                          << '\n';
            }
            return 0;
        }
        if (argc < 4) {
            throw std::runtime_error(
                "usage: analyze_su2_plucker_kernel TARGET MINUS_MASK LABEL..."
            );
        }
        const int target = parse_nonnegative(argv[1], "target");
        const int parsed_mask = parse_nonnegative(argv[2], "minus mask");
        std::vector<int> labels;
        for (int argument = 3; argument < argc; ++argument) {
            labels.push_back(parse_nonnegative(argv[argument], "label"));
        }
        if (labels.size() >= std::numeric_limits<unsigned int>::digits
            || parsed_mask >= static_cast<int>(1U << labels.size())) {
            throw std::runtime_error("case bounds are out of range");
        }
        const unsigned int minus_mask
            = static_cast<unsigned int>(parsed_mask);
        const Result result = analyze_case(labels, target, minus_mask);
        const std::uint64_t kernel
            = result.odd_sources - static_cast<std::uint64_t>(result.rank);
        std::cout << "labels=";
        print_vector(labels);
        std::cout << " target=" << target
                  << " minus_mask=" << minus_mask
                  << " codomain_dimension=" << result.codomain_dimension
                  << " odd_sources=" << result.odd_sources
                  << " original_rank=" << result.rank
                  << " original_kernel=" << kernel
                  << " positive_sources=" << result.positive_sources
                  << " straightening_states="
                  << result.straightening_states
                  << " raw_collision_relations="
                  << result.raw_collision_relations
                  << " pair_intersection_relations="
                  << result.pair_intersection_relations
                  << " edge_exchange_relations="
                  << result.edge_exchange_relations
                  << " cut_boundary_relations="
                  << result.cut_boundary_relations
                  << " cut_boundary_rank="
                  << result.cut_boundary_rank
                  << " cut_order_best_rank="
                  << (result.cut_order_ranks.empty()
                          ? 0U
                          : *std::max_element(
                                result.cut_order_ranks.begin(),
                                result.cut_order_ranks.end()
                            ))
                  << " cut_order_combined_rank="
                  << result.cut_order_combined_rank
                  << " section_defect_relations="
                  << result.section_defect_relations
                  << " candidate_relation_rank="
                  << result.candidate_relation_rank
                  << " candidate_relations_valid="
                  << (result.candidate_relations_valid ? "PASS" : "FAIL")
                  << " fundamental_lifts_valid="
                  << (result.fundamental_lifts_valid ? "PASS" : "FAIL")
                  << " fundamental_owned_lifts_valid="
                  << (result.fundamental_owned_lifts_valid
                          ? "PASS" : "FAIL")
                  << " fundamental_combined_valid="
                  << (result.fundamental_combined_valid ? "PASS" : "FAIL")
                  << " fundamental_kernel_quotient_rank="
                  << result.fundamental_kernel_quotient_rank
                  << " fundamental_candidate_quotient_rank="
                  << result.fundamental_candidate_quotient_rank
                  << " fundamental_covered_rank="
                  << result.fundamental_covered_rank
                  << " positive_internal_kernel_rank="
                  << result.positive_internal_kernel_rank
                  << " allocation_dimension="
                  << result.allocation_dimension
                  << " carrier_dimension="
                  << result.carrier_dimension
                  << " f0_dimension=" << result.f0_dimension
                  << " f0_kernel_dimension="
                  << result.f0_kernel_dimension
                  << " negative_stripped_rank="
                  << result.negative_stripped_rank
                  << " positive_stripped_rank="
                  << result.positive_stripped_rank
                  << " positive_global_rank="
                  << result.positive_global_rank
                  << " carrier_overlap_dimension="
                  << result.carrier_overlap_dimension
                  << " positive_allocation_dimension="
                  << result.positive_allocation_dimension
                  << " unshared_negative_carrier_dimension="
                  << result.unshared_negative_carrier_dimension
                  << " residual_negative_demand="
                  << result.residual_negative_demand
                  << " residual_positive_supply="
                  << result.residual_positive_supply
                  << " optimal_section_cover_feasible="
                  << (result.optimal_section_cover_feasible
                          ? "PASS" : "FAIL")
                  << " optimal_defect_rank="
                  << result.optimal_defect_rank
                  << " optimal_defect_kernel="
                  << result.optimal_defect_kernel
                  << " optimal_relation_rank="
                  << result.optimal_relation_rank
                  << " merged_row_kernel_dimension="
                  << result.merged_row_kernel_dimension
                  << " deletion_contraction_valid="
                  << (result.deletion_contraction_valid ? "PASS" : "FAIL")
                  << " projected_allocation_dimension="
                  << result.projected_allocation_dimension
                  << " projected_positive_allocation_dimension="
                  << result.projected_positive_allocation_dimension
                  << " projected_unshared_carrier_dimension="
                  << result.projected_unshared_carrier_dimension
                  << " projected_positive_rank="
                  << result.projected_positive_rank
                  << " negative_lower_branch_dimension="
                  << result.negative_lower_branch_dimension
                  << " positive_lower_branch_dimension="
                  << result.positive_lower_branch_dimension
                  << " negative_lift_obstruction_rank="
                  << result.negative_lift_obstruction_rank
                  << " positive_lift_obstruction_rank="
                  << result.positive_lift_obstruction_rank
                  << " enriched_deletion_contraction_valid="
                  << (result.enriched_deletion_contraction_valid
                          ? "PASS" : "FAIL")
                  << " projected_carrier_criterion="
                  << (result.projected_carrier_criterion ? "PASS" : "FAIL")
                  << " ambient_kernel_criterion="
                  << (result.ambient_kernel_criterion ? "PASS" : "FAIL")
                  << " optimal_section_rank="
                  << result.optimal_section_rank
                  << " optimal_section_criterion="
                  << (result.optimal_section_criterion ? "PASS" : "FAIL")
                  << " result="
                  << (kernel <= result.positive_sources ? "PASS" : "FAIL")
                  << '\n';
        if (!result.section_trial_ranks.empty()) {
            std::cout << "cut_order_ranks=";
            for (std::size_t trial = 0U;
                 trial < result.cut_order_ranks.size(); ++trial) {
                if (trial != 0U) {
                    std::cout << ',';
                }
                std::cout << result.cut_order_ranks[trial];
            }
            std::cout << '\n';
            std::cout << "section_trial_ranks=";
            for (std::size_t trial = 0U;
                 trial < result.section_trial_ranks.size(); ++trial) {
                if (trial != 0U) {
                    std::cout << ',';
                }
                std::cout << result.section_trial_ranks[trial];
            }
            std::cout << '\n';
            std::cout << "section_individual_ranks=";
            for (std::size_t trial = 0U;
                 trial < result.section_individual_ranks.size(); ++trial) {
                if (trial != 0U) {
                    std::cout << ',';
                }
                std::cout << result.section_individual_ranks[trial];
            }
            std::cout << '\n';
            std::cout << "section_criteria=";
            for (std::size_t trial = 0U;
                 trial < result.section_criteria.size(); ++trial) {
                if (trial != 0U) {
                    std::cout << ',';
                }
                const Result::SectionCriterion& criterion
                    = result.section_criteria[trial];
                std::cout << criterion.rank << ':' << criterion.kernel
                          << ':' << criterion.carrier_intersection << ':'
                          << (criterion.pass ? "PASS" : "FAIL");
            }
            std::cout << '\n';
        }
        for (const Result::PairIntersection& pair
             : result.pair_intersections) {
            std::cout << "pair first=" << pair.first
                      << " second=" << pair.second
                      << " intersection=" << pair.dimension
                      << " xor=" << (pair.first ^ pair.second)
                      << " xor_channel_dimension="
                      << pair.xor_channel_dimension << '\n';
        }
        for (std::size_t relation = 0U;
             relation < result.nullspace.size(); ++relation) {
            std::cout << "null_relation=" << relation;
            for (std::size_t column = 0U;
                 column < result.nullspace[relation].size(); ++column) {
                const Rational& coefficient
                    = result.nullspace[relation][column];
                if (coefficient == 0) {
                    continue;
                }
                std::cout << " mask=" << result.column_labels[column].first
                          << " basis=" << result.column_labels[column].second
                          << " coefficient=" << coefficient.numerator();
                if (coefficient.denominator() != 1) {
                    std::cout << '/' << coefficient.denominator();
                }
            }
            std::cout << '\n';
        }
        if (std::getenv("SU2_PLUCKER_VERBOSE") != nullptr) {
            for (std::size_t column = 0U;
                 column < result.column_graphs.size(); ++column) {
                std::cout << "odd_source mask="
                          << result.column_labels[column].first
                          << " basis="
                          << result.column_labels[column].second
                          << " graph=";
                print_graph(result.column_graphs[column]);
                std::cout << '\n';
            }
            for (const Result::PositiveSource& source
                 : result.positive_source_graphs) {
                std::cout << "positive_source mask=" << source.mask
                          << " basis=" << source.basis
                          << " graph=";
                print_graph(source.graph);
                std::cout << '\n';
            }
            for (std::size_t relation = 0U;
                 relation < result.fundamental_lifts.size(); ++relation) {
                std::cout << "fundamental_lift relation=" << relation
                          << " terms="
                          << result.fundamental_lifts[relation].size()
                          << '\n';
                for (const Result::FundamentalLiftTerm& term
                     : result.fundamental_lifts[relation]) {
                    std::cout << "  triple=" << term.triple[0] << ','
                              << term.triple[1] << ',' << term.triple[2]
                              << " coefficient=" << term.scalar.numerator();
                    if (term.scalar.denominator() != 1) {
                        std::cout << '/' << term.scalar.denominator();
                    }
                    std::cout << " graph=";
                    print_graph(term.coefficient);
                    std::cout << " owners=";
                    if (term.owners.empty()) {
                        std::cout << "none";
                    }
                    for (std::size_t owner = 0U;
                         owner < term.owners.size(); ++owner) {
                        if (owner != 0U) {
                            std::cout << ';';
                        }
                        const auto& assignment = term.owners[owner];
                        std::cout << assignment.mask << ':'
                                  << assignment.basis << ":f"
                                  << assignment.face << ':'
                                  << assignment.ratio.numerator();
                        if (assignment.ratio.denominator() != 1) {
                            std::cout << '/'
                                      << assignment.ratio.denominator();
                        }
                    }
                    std::cout << '\n';
                }
            }
            if (!result.fundamental_missing_relation.empty()) {
                std::cout << "fundamental_missing_relation";
                for (std::size_t source = 0U;
                     source < result.fundamental_missing_relation.size();
                     ++source) {
                    if (result.fundamental_missing_relation[source] != 0) {
                        std::cout << " mask="
                                  << result.column_labels[source].first
                                  << " basis="
                                  << result.column_labels[source].second
                                  << " coefficient="
                                  << result.fundamental_missing_relation[
                                         source
                                     ];
                    }
                }
                std::cout << '\n';
            }
        }
        return kernel <= result.positive_sources ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
