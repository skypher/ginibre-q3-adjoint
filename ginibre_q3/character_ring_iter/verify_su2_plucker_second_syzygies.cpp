#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <omp.h>

namespace {

constexpr std::int64_t modulus = 1'000'000'007LL;
using Graph = std::vector<std::vector<int>>;
using Key = std::vector<int>;
using Expansion = std::map<Key, std::int64_t>;
using Support = std::set<Key>;
using SparseVector = std::map<std::size_t, std::int64_t>;

int popcount(unsigned int value) {
    int answer = 0;
    while (value != 0U) {
        answer += static_cast<int>(value & 1U);
        value >>= 1U;
    }
    return answer;
}

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

std::int64_t normalized(std::int64_t value) {
    value %= modulus;
    if (value < 0) {
        value += modulus;
    }
    return value;
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

std::int64_t deterministic_mask_weight(
    unsigned int mask,
    std::size_t trial
) {
    std::uint64_t value = static_cast<std::uint64_t>(mask)
        + 0x9e3779b97f4a7c15ULL
            * static_cast<std::uint64_t>(trial + 1U);
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    value ^= value >> 31U;
    return static_cast<std::int64_t>(
        value % static_cast<std::uint64_t>(modulus - 1)
    ) + 1;
}

std::size_t modular_rank(std::vector<std::vector<std::int64_t>> matrix) {
    if (matrix.empty()) {
        return 0U;
    }
    const std::size_t rows = matrix.size();
    const std::size_t columns = matrix.front().size();
    std::size_t rank = 0U;
    for (std::size_t column = 0U; column < columns && rank < rows; ++column) {
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
        for (std::size_t i = 0U; i < rows; ++i) {
            if (i == rank || matrix[i][column] == 0) {
                continue;
            }
            const std::int64_t factor = matrix[i][column];
            for (std::size_t j = column; j < columns; ++j) {
                matrix[i][j] = normalized(
                    matrix[i][j] - factor * matrix[rank][j]
                );
            }
        }
        ++rank;
    }
    return rank;
}

std::vector<std::int64_t> binary_vertex_factor_coefficients(
    const Graph& graph,
    std::size_t vertex
) {
    if (vertex >= graph.size()) {
        throw std::runtime_error("binary-factor vertex out of range");
    }
    std::vector<std::int64_t> coefficients{1};
    for (std::size_t neighbor = 0U; neighbor < graph.size(); ++neighbor) {
        if (neighbor == vertex) {
            continue;
        }
        const std::size_t left = std::min(vertex, neighbor);
        const std::size_t right = std::max(vertex, neighbor);
        const std::int64_t slope
            = static_cast<std::int64_t>(neighbor + 1U);
        for (int copy = 0; copy < graph[left][right]; ++copy) {
            std::vector<std::int64_t> next(coefficients.size() + 1U, 0);
            for (std::size_t power = 0U;
                 power < coefficients.size(); ++power) {
                next[power] = normalized(next[power] + coefficients[power]);
                next[power + 1U] = normalized(
                    next[power + 1U] - slope * coefficients[power]
                );
            }
            coefficients = std::move(next);
        }
    }
    return coefficients;
}

std::size_t binary_vertex_factor_rank(
    const std::vector<Graph>& graphs,
    std::size_t vertex
) {
    std::vector<std::vector<std::int64_t>> rows;
    rows.reserve(graphs.size());
    for (const Graph& graph : graphs) {
        rows.push_back(binary_vertex_factor_coefficients(graph, vertex));
    }
    return modular_rank(std::move(rows));
}

std::size_t binary_terminal_tensor_rank(
    const std::vector<Graph>& graphs,
    std::size_t first_vertex,
    std::size_t second_vertex
) {
    std::vector<std::vector<std::int64_t>> rows;
    rows.reserve(graphs.size());
    for (const Graph& graph : graphs) {
        const std::vector<std::int64_t> first
            = binary_vertex_factor_coefficients(graph, first_vertex);
        const std::vector<std::int64_t> second
            = binary_vertex_factor_coefficients(graph, second_vertex);
        std::vector<std::int64_t> tensor(
            first.size() * second.size(), 0
        );
        for (std::size_t i = 0U; i < first.size(); ++i) {
            for (std::size_t j = 0U; j < second.size(); ++j) {
                tensor[i * second.size() + j]
                    = (first[i] * second[j]) % modulus;
            }
        }
        rows.push_back(std::move(tensor));
    }
    return modular_rank(std::move(rows));
}

bool binary_terminal_jet_separates(
    const std::vector<Graph>& graphs,
    std::size_t jet_vertex,
    std::size_t remaining_vertex,
    std::size_t root_vertex
) {
    std::map<int, std::vector<Graph>> vanishing_order_groups;
    for (const Graph& graph : graphs) {
        if (jet_vertex >= graph.size() || remaining_vertex >= graph.size()
            || root_vertex >= graph.size() || root_vertex == jet_vertex) {
            throw std::runtime_error("binary jet vertex out of range");
        }
        const std::size_t left = std::min(jet_vertex, root_vertex);
        const std::size_t right = std::max(jet_vertex, root_vertex);
        vanishing_order_groups[graph[left][right]].push_back(graph);
    }
    for (const auto& [order, group] : vanishing_order_groups) {
        (void)order;
        if (binary_vertex_factor_rank(group, remaining_vertex)
            != group.size()) {
            return false;
        }
    }
    return true;
}

bool binary_terminal_double_valuation_separates(
    const std::vector<Graph>& graphs,
    std::size_t first_terminal,
    std::size_t second_terminal,
    std::size_t first_root
) {
    std::map<int, std::vector<Graph>> groups;
    for (const Graph& graph : graphs) {
        const std::size_t left = std::min(first_terminal, first_root);
        const std::size_t right = std::max(first_terminal, first_root);
        groups[graph[left][right]].push_back(graph);
    }
    for (const auto& [order, group] : groups) {
        (void)order;
        bool group_separated = group.size() <= 1U;
        for (std::size_t second_root = 2U;
             second_root < 6U && !group_separated; ++second_root) {
            std::set<int> orders;
            for (const Graph& graph : group) {
                const std::size_t left
                    = std::min(second_terminal, second_root);
                const std::size_t right
                    = std::max(second_terminal, second_root);
                orders.insert(graph[left][right]);
            }
            group_separated = orders.size() == group.size();
        }
        if (!group_separated) {
            return false;
        }
    }
    return true;
}

bool binary_gcd_hall_condition(
    const std::vector<Graph>& graphs,
    std::size_t terminal
) {
    if (graphs.size() <= 1U) {
        return true;
    }
    if (graphs.size() >= std::numeric_limits<unsigned int>::digits) {
        throw std::runtime_error("too many rows for binary Hall subsets");
    }
    int degree = 0;
    for (std::size_t root = 2U; root < 6U; ++root) {
        const std::size_t left = std::min(terminal, root);
        const std::size_t right = std::max(terminal, root);
        degree += graphs.front()[left][right];
    }
    const unsigned int subset_limit = 1U << graphs.size();
    for (unsigned int subset = 1U; subset < subset_limit; ++subset) {
        int gcd_degree = 0;
        for (std::size_t root = 2U; root < 6U; ++root) {
            int minimum = degree;
            for (std::size_t row = 0U; row < graphs.size(); ++row) {
                if (((subset >> row) & 1U) == 0U) {
                    continue;
                }
                const std::size_t left = std::min(terminal, root);
                const std::size_t right = std::max(terminal, root);
                minimum = std::min(minimum, graphs[row][left][right]);
            }
            gcd_degree += minimum;
        }
        if (popcount(subset) > degree - gcd_degree + 1) {
            return false;
        }
    }
    return true;
}

bool binary_terminal_hall_jet_separates(
    const std::vector<Graph>& graphs,
    std::size_t jet_vertex,
    std::size_t remaining_vertex,
    std::size_t root_vertex
) {
    std::map<int, std::vector<Graph>> groups;
    for (const Graph& graph : graphs) {
        const std::size_t left = std::min(jet_vertex, root_vertex);
        const std::size_t right = std::max(jet_vertex, root_vertex);
        groups[graph[left][right]].push_back(graph);
    }
    for (const auto& [order, group] : groups) {
        (void)order;
        if (!binary_gcd_hall_condition(group, remaining_vertex)) {
            return false;
        }
    }
    return true;
}

bool binary_lovett_vstar_shape(
    const std::vector<Graph>& graphs,
    std::size_t terminal
) {
    if (graphs.empty()) {
        return true;
    }
    std::array<int, 4U> common{};
    common.fill(std::numeric_limits<int>::max());
    for (const Graph& graph : graphs) {
        for (std::size_t root = 2U; root < 6U; ++root) {
            const std::size_t left = std::min(terminal, root);
            const std::size_t right = std::max(terminal, root);
            common[root - 2U] = std::min(
                common[root - 2U], graph[left][right]
            );
        }
    }
    for (std::size_t distinguished = 0U;
         distinguished < common.size(); ++distinguished) {
        bool valid = true;
        for (const Graph& graph : graphs) {
            for (std::size_t coordinate = 0U;
                 coordinate < common.size(); ++coordinate) {
                if (coordinate == distinguished) {
                    continue;
                }
                const std::size_t root = coordinate + 2U;
                const std::size_t left = std::min(terminal, root);
                const std::size_t right = std::max(terminal, root);
                const int residual
                    = graph[left][right] - common[coordinate];
                if (residual < 0 || residual > 1) {
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                break;
            }
        }
        if (valid) {
            return true;
        }
    }
    return false;
}

bool binary_terminal_lovett_jet_separates(
    const std::vector<Graph>& graphs,
    std::size_t jet_vertex,
    std::size_t remaining_vertex,
    std::size_t root_vertex
) {
    std::map<int, std::vector<Graph>> groups;
    for (const Graph& graph : graphs) {
        const std::size_t left = std::min(jet_vertex, root_vertex);
        const std::size_t right = std::max(jet_vertex, root_vertex);
        groups[graph[left][right]].push_back(graph);
    }
    for (const auto& [order, group] : groups) {
        (void)order;
        if (!binary_gcd_hall_condition(group, remaining_vertex)
            || !binary_lovett_vstar_shape(group, remaining_vertex)) {
            return false;
        }
    }
    return true;
}

bool exponent_rows_merge_certificate(
    const std::vector<std::vector<int>>& rows,
    std::pair<std::size_t, std::size_t>* witness = nullptr
) {
    if (rows.size() <= 1U) {
        if (witness != nullptr) {
            *witness = {0U, 1U};
        }
        return true;
    }
    const std::size_t coordinates = rows.front().size();
    const auto specialized_rank = [&rows]() {
        std::vector<std::vector<std::int64_t>> coefficient_rows;
        coefficient_rows.reserve(rows.size());
        for (const std::vector<int>& row : rows) {
            std::vector<std::int64_t> coefficients{1};
            for (std::size_t coordinate = 0U;
                 coordinate < row.size(); ++coordinate) {
                const std::int64_t slope
                    = static_cast<std::int64_t>(coordinate + 1U);
                for (int copy = 0; copy < row[coordinate]; ++copy) {
                    std::vector<std::int64_t> next(
                        coefficients.size() + 1U, 0
                    );
                    for (std::size_t power = 0U;
                         power < coefficients.size(); ++power) {
                        next[power] = normalized(
                            next[power] + coefficients[power]
                        );
                        next[power + 1U] = normalized(
                            next[power + 1U]
                            - slope * coefficients[power]
                        );
                    }
                    coefficients = std::move(next);
                }
            }
            coefficient_rows.push_back(std::move(coefficients));
        }
        return modular_rank(std::move(coefficient_rows));
    };
    if (coordinates <= 3U) {
        return specialized_rank() == rows.size();
    }
    for (std::size_t first = 0U; first < coordinates; ++first) {
        for (std::size_t second = first + 1U; second < coordinates;
             ++second) {
            std::vector<std::vector<int>> merged;
            merged.reserve(rows.size());
            for (const std::vector<int>& row : rows) {
                std::vector<int> next;
                next.reserve(coordinates - 1U);
                for (std::size_t coordinate = 0U;
                     coordinate < coordinates; ++coordinate) {
                    if (coordinate == second) {
                        continue;
                    }
                    next.push_back(
                        coordinate == first
                            ? row[first] + row[second]
                            : row[coordinate]
                    );
                }
                merged.push_back(std::move(next));
            }
            if (exponent_rows_merge_certificate(merged)) {
                if (witness != nullptr) {
                    *witness = {first, second};
                }
                return true;
            }
        }
    }
    return false;
}

bool exponent_rows_gcd_hall_condition(
    const std::vector<std::vector<int>>& rows
) {
    if (rows.size() <= 1U) {
        return true;
    }
    if (rows.size() >= std::numeric_limits<unsigned int>::digits) {
        throw std::runtime_error("too many exponent rows for Hall subsets");
    }
    const int degree = std::accumulate(
        rows.front().begin(), rows.front().end(), 0
    );
    for (const std::vector<int>& row : rows) {
        if (std::accumulate(row.begin(), row.end(), 0) != degree) {
            throw std::runtime_error("unequal exponent-row degrees");
        }
    }
    const unsigned int subset_limit = 1U << rows.size();
    for (unsigned int subset = 1U; subset < subset_limit; ++subset) {
        int gcd_degree = 0;
        for (std::size_t coordinate = 0U;
             coordinate < rows.front().size(); ++coordinate) {
            int minimum = degree;
            for (std::size_t row = 0U; row < rows.size(); ++row) {
                if (((subset >> row) & 1U) != 0U) {
                    minimum = std::min(minimum, rows[row][coordinate]);
                }
            }
            gcd_degree += minimum;
        }
        if (popcount(subset) > degree - gcd_degree + 1) {
            return false;
        }
    }
    return true;
}

bool exponent_rows_lovett_vstar_shape(
    const std::vector<std::vector<int>>& rows
) {
    if (rows.empty()) {
        return true;
    }
    std::vector<int> common(
        rows.front().size(), std::numeric_limits<int>::max()
    );
    for (const std::vector<int>& row : rows) {
        for (std::size_t coordinate = 0U;
             coordinate < row.size(); ++coordinate) {
            common[coordinate] = std::min(
                common[coordinate], row[coordinate]
            );
        }
    }
    for (std::size_t distinguished = 0U;
         distinguished < common.size(); ++distinguished) {
        bool valid = true;
        for (const std::vector<int>& row : rows) {
            for (std::size_t coordinate = 0U;
                 coordinate < row.size(); ++coordinate) {
                if (coordinate == distinguished) {
                    continue;
                }
                const int residual
                    = row[coordinate] - common[coordinate];
                if (residual < 0 || residual > 1) {
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                break;
            }
        }
        if (valid) {
            return true;
        }
    }
    return false;
}

bool exponent_rows_v2_shape(
    const std::vector<std::vector<int>>& rows
) {
    if (rows.empty()) {
        return true;
    }
    std::vector<int> common(
        rows.front().size(), std::numeric_limits<int>::max()
    );
    for (const std::vector<int>& row : rows) {
        for (std::size_t coordinate = 0U;
             coordinate < row.size(); ++coordinate) {
            common[coordinate] = std::min(
                common[coordinate], row[coordinate]
            );
        }
    }
    std::size_t binary_coordinates = 0U;
    for (std::size_t coordinate = 0U;
         coordinate < common.size(); ++coordinate) {
        bool binary = true;
        for (const std::vector<int>& row : rows) {
            const int residual = row[coordinate] - common[coordinate];
            if (residual < 0 || residual > 1) {
                binary = false;
                break;
            }
        }
        if (binary) {
            ++binary_coordinates;
        }
    }
    return binary_coordinates + 2U >= common.size();
}

bool exponent_rows_lovett_merge_certificate(
    const std::vector<std::vector<int>>& rows,
    std::pair<std::size_t, std::size_t>* witness = nullptr
) {
    if (rows.size() <= 1U) {
        if (witness != nullptr) {
            *witness = {0U, 1U};
        }
        return true;
    }
    const std::size_t coordinates = rows.front().size();
    for (std::size_t first = 0U; first < coordinates; ++first) {
        for (std::size_t second = first + 1U; second < coordinates;
             ++second) {
            std::vector<std::vector<int>> merged;
            merged.reserve(rows.size());
            for (const std::vector<int>& row : rows) {
                std::vector<int> next;
                next.reserve(coordinates - 1U);
                for (std::size_t coordinate = 0U;
                     coordinate < coordinates; ++coordinate) {
                    if (coordinate == second) {
                        continue;
                    }
                    next.push_back(
                        coordinate == first
                            ? row[first] + row[second]
                            : row[coordinate]
                    );
                }
                merged.push_back(std::move(next));
            }
            if (exponent_rows_lovett_vstar_shape(merged)
                && exponent_rows_gcd_hall_condition(merged)) {
                if (witness != nullptr) {
                    *witness = {first, second};
                }
                return true;
            }
        }
    }
    return false;
}

bool exponent_rows_v2_merge_certificate(
    const std::vector<std::vector<int>>& rows
) {
    if (rows.size() <= 1U) {
        return true;
    }
    const std::size_t coordinates = rows.front().size();
    for (std::size_t first = 0U; first < coordinates; ++first) {
        for (std::size_t second = first + 1U; second < coordinates;
             ++second) {
            std::vector<std::vector<int>> merged;
            merged.reserve(rows.size());
            for (const std::vector<int>& row : rows) {
                std::vector<int> next;
                next.reserve(coordinates - 1U);
                for (std::size_t coordinate = 0U;
                     coordinate < coordinates; ++coordinate) {
                    if (coordinate == second) {
                        continue;
                    }
                    next.push_back(
                        coordinate == first
                            ? row[first] + row[second]
                            : row[coordinate]
                    );
                }
                merged.push_back(std::move(next));
            }
            if (exponent_rows_v2_shape(merged)
                && exponent_rows_gcd_hall_condition(merged)) {
                return true;
            }
        }
    }
    return false;
}

bool binary_merge_certifies(
    const std::vector<Graph>& graphs,
    std::size_t terminal,
    std::pair<std::size_t, std::size_t>* witness = nullptr
) {
    std::vector<std::vector<int>> rows;
    rows.reserve(graphs.size());
    for (const Graph& graph : graphs) {
        std::vector<int> row;
        row.reserve(4U);
        for (std::size_t root = 2U; root < 6U; ++root) {
            const std::size_t left = std::min(terminal, root);
            const std::size_t right = std::max(terminal, root);
            row.push_back(graph[left][right]);
        }
        rows.push_back(std::move(row));
    }
    return exponent_rows_merge_certificate(rows, witness);
}

bool binary_lovett_merge_certifies(
    const std::vector<Graph>& graphs,
    std::size_t terminal,
    std::pair<std::size_t, std::size_t>* witness = nullptr
) {
    std::vector<std::vector<int>> rows;
    rows.reserve(graphs.size());
    for (const Graph& graph : graphs) {
        std::vector<int> row;
        row.reserve(4U);
        for (std::size_t root = 2U; root < 6U; ++root) {
            const std::size_t left = std::min(terminal, root);
            const std::size_t right = std::max(terminal, root);
            row.push_back(graph[left][right]);
        }
        rows.push_back(std::move(row));
    }
    return exponent_rows_lovett_merge_certificate(rows, witness);
}

bool binary_v2_merge_certifies(
    const std::vector<Graph>& graphs,
    std::size_t terminal
) {
    std::vector<std::vector<int>> rows;
    rows.reserve(graphs.size());
    for (const Graph& graph : graphs) {
        std::vector<int> row;
        row.reserve(4U);
        for (std::size_t root = 2U; root < 6U; ++root) {
            const std::size_t left = std::min(terminal, root);
            const std::size_t right = std::max(terminal, root);
            row.push_back(graph[left][right]);
        }
        rows.push_back(std::move(row));
    }
    return exponent_rows_v2_merge_certificate(rows);
}

bool binary_v2_certifies(
    const std::vector<Graph>& graphs,
    std::size_t terminal
) {
    std::vector<std::vector<int>> rows;
    rows.reserve(graphs.size());
    for (const Graph& graph : graphs) {
        std::vector<int> row;
        row.reserve(4U);
        for (std::size_t root = 2U; root < 6U; ++root) {
            const std::size_t left = std::min(terminal, root);
            const std::size_t right = std::max(terminal, root);
            row.push_back(graph[left][right]);
        }
        rows.push_back(std::move(row));
    }
    return exponent_rows_v2_shape(rows)
        && exponent_rows_gcd_hall_condition(rows);
}

bool binary_terminal_merge_jet_separates(
    const std::vector<Graph>& graphs,
    std::size_t jet_vertex,
    std::size_t remaining_vertex,
    std::size_t root_vertex,
    std::vector<int>* witness_signature = nullptr
) {
    if (witness_signature != nullptr) {
        witness_signature->clear();
    }
    std::map<int, std::vector<Graph>> groups;
    for (const Graph& graph : graphs) {
        const std::size_t left = std::min(jet_vertex, root_vertex);
        const std::size_t right = std::max(jet_vertex, root_vertex);
        groups[graph[left][right]].push_back(graph);
    }
    for (const auto& [order, group] : groups) {
        std::pair<std::size_t, std::size_t> merger{0U, 1U};
        if (!binary_merge_certifies(group, remaining_vertex, &merger)) {
            if (witness_signature != nullptr) {
                witness_signature->clear();
            }
            return false;
        }
        if (witness_signature != nullptr) {
            witness_signature->push_back(order);
            witness_signature->push_back(static_cast<int>(group.size()));
            witness_signature->push_back(static_cast<int>(
                merger.first * 4U + merger.second
            ));
        }
    }
    return true;
}

bool binary_terminal_lovett_merge_jet_separates(
    const std::vector<Graph>& graphs,
    std::size_t jet_vertex,
    std::size_t remaining_vertex,
    std::size_t root_vertex
) {
    std::map<int, std::vector<Graph>> groups;
    for (const Graph& graph : graphs) {
        const std::size_t left = std::min(jet_vertex, root_vertex);
        const std::size_t right = std::max(jet_vertex, root_vertex);
        groups[graph[left][right]].push_back(graph);
    }
    for (const auto& [order, group] : groups) {
        (void)order;
        if (!binary_lovett_merge_certifies(group, remaining_vertex)) {
            return false;
        }
    }
    return true;
}

bool binary_terminal_v2_merge_jet_separates(
    const std::vector<Graph>& graphs,
    std::size_t jet_vertex,
    std::size_t remaining_vertex,
    std::size_t root_vertex
) {
    std::map<int, std::vector<Graph>> groups;
    for (const Graph& graph : graphs) {
        const std::size_t left = std::min(jet_vertex, root_vertex);
        const std::size_t right = std::max(jet_vertex, root_vertex);
        groups[graph[left][right]].push_back(graph);
    }
    for (const auto& [order, group] : groups) {
        (void)order;
        if (!binary_v2_merge_certifies(group, remaining_vertex)) {
            return false;
        }
    }
    return true;
}

bool binary_terminal_v2_jet_separates(
    const std::vector<Graph>& graphs,
    std::size_t jet_vertex,
    std::size_t remaining_vertex,
    std::size_t root_vertex
) {
    std::map<int, std::vector<Graph>> groups;
    for (const Graph& graph : graphs) {
        const std::size_t left = std::min(jet_vertex, root_vertex);
        const std::size_t right = std::max(jet_vertex, root_vertex);
        groups[graph[left][right]].push_back(graph);
    }
    for (const auto& [order, group] : groups) {
        (void)order;
        if (!binary_v2_certifies(group, remaining_vertex)) {
            return false;
        }
    }
    return true;
}

std::vector<std::vector<std::int64_t>> modular_left_kernel(
    const std::vector<std::vector<std::int64_t>>& matrix
) {
    if (matrix.empty()) {
        return {};
    }
    const std::size_t source_count = matrix.size();
    const std::size_t target_count = matrix.front().size();
    std::vector<std::vector<std::int64_t>> transpose(
        target_count, std::vector<std::int64_t>(source_count, 0)
    );
    for (std::size_t source = 0U; source < source_count; ++source) {
        for (std::size_t target = 0U; target < target_count; ++target) {
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
                transpose[i][j] = normalized(
                    transpose[i][j] - factor * transpose[row][j]
                );
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

std::pair<bool, std::vector<std::int64_t>> solve_left_combination(
    const std::vector<std::vector<std::int64_t>>& rows,
    const std::vector<std::int64_t>& target
) {
    if (rows.empty()) {
        const bool zero = std::all_of(
            target.begin(), target.end(),
            [](std::int64_t value) { return value == 0; }
        );
        return {zero, {}};
    }
    const std::size_t variables = rows.size();
    const std::size_t equations = rows.front().size();
    if (target.size() != equations) {
        throw std::runtime_error("left-combination target dimension mismatch");
    }
    std::vector<std::vector<std::int64_t>> augmented(
        equations, std::vector<std::int64_t>(variables + 1U, 0)
    );
    for (std::size_t equation = 0U; equation < equations; ++equation) {
        for (std::size_t variable = 0U; variable < variables; ++variable) {
            augmented[equation][variable] = rows[variable][equation];
        }
        augmented[equation][variables] = target[equation];
    }
    std::vector<std::size_t> pivot_columns;
    std::size_t pivot_row = 0U;
    for (std::size_t column = 0U;
         column < variables && pivot_row < equations; ++column) {
        std::size_t pivot = pivot_row;
        while (pivot < equations && augmented[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == equations) {
            continue;
        }
        std::swap(augmented[pivot_row], augmented[pivot]);
        const std::int64_t inverse = modular_power(
            augmented[pivot_row][column], modulus - 2
        );
        for (std::size_t j = column; j <= variables; ++j) {
            augmented[pivot_row][j] =
                (augmented[pivot_row][j] * inverse) % modulus;
        }
        for (std::size_t row = 0U; row < equations; ++row) {
            if (row == pivot_row || augmented[row][column] == 0) {
                continue;
            }
            const std::int64_t factor = augmented[row][column];
            for (std::size_t j = column; j <= variables; ++j) {
                augmented[row][j] = normalized(
                    augmented[row][j]
                    - factor * augmented[pivot_row][j]
                );
            }
        }
        pivot_columns.push_back(column);
        ++pivot_row;
    }
    for (std::size_t row = pivot_row; row < equations; ++row) {
        bool zero = true;
        for (std::size_t column = 0U; column < variables; ++column) {
            if (augmented[row][column] != 0) {
                zero = false;
                break;
            }
        }
        if (zero && augmented[row][variables] != 0) {
            return {false, {}};
        }
    }
    std::vector<std::int64_t> solution(variables, 0);
    for (std::size_t row = 0U; row < pivot_columns.size(); ++row) {
        solution[pivot_columns[row]] = augmented[row][variables];
    }
    return {true, std::move(solution)};
}

Key flatten(const Graph& graph) {
    Key key;
    const std::size_t vertices = graph.size();
    key.reserve(vertices * (vertices - 1U) / 2U);
    for (std::size_t i = 0U; i < vertices; ++i) {
        for (std::size_t j = i + 1U; j < vertices; ++j) {
            key.push_back(graph[i][j]);
        }
    }
    return key;
}

Graph unflatten(const Key& key, std::size_t vertices) {
    Graph graph(vertices, std::vector<int>(vertices, 0));
    std::size_t position = 0U;
    for (std::size_t i = 0U; i < vertices; ++i) {
        for (std::size_t j = i + 1U; j < vertices; ++j) {
            graph[i][j] = key[position];
            ++position;
        }
    }
    return graph;
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
    const std::size_t vertices = graph.size();
    for (std::size_t i = 0U; i < vertices; ++i) {
        for (std::size_t j = i + 1U; j < vertices; ++j) {
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
    const std::size_t vertices = graph.size();
    if (neighbor == vertices) {
        if (needed == 0) {
            const int saved = remaining[vertex];
            remaining[vertex] = 0;
            enumerate_graphs_rec(remaining, graph, output);
            remaining[vertex] = saved;
        }
        return;
    }
    int available = 0;
    for (std::size_t j = neighbor; j < vertices; ++j) {
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
    const std::size_t vertices = graph.size();
    std::size_t vertex = 0U;
    while (vertex < vertices && remaining[vertex] == 0) {
        ++vertex;
    }
    if (vertex == vertices) {
        output.push_back(graph);
        return;
    }
    distribute_vertex_degree(
        vertex, vertex + 1U, remaining[vertex], remaining, graph, output
    );
}

std::vector<Graph> noncrossing_graphs(const std::vector<int>& degrees) {
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

const Expansion& straighten(
    const Graph& graph,
    std::map<Key, Expansion>& memo
);

const Support& straighten_support(
    const Graph& graph,
    std::map<Key, Support>& memo
);

int cut_weight(const Graph& graph, unsigned int split);

Graph permute_graph(
    const Graph& graph,
    const std::vector<std::size_t>& order
) {
    const std::size_t vertices = graph.size();
    if (order.size() != vertices) {
        throw std::runtime_error("graph permutation has wrong size");
    }
    Graph result(vertices, std::vector<int>(vertices, 0));
    for (std::size_t i = 0U; i < vertices; ++i) {
        for (std::size_t j = i + 1U; j < vertices; ++j) {
            const std::size_t old_i = order[i];
            const std::size_t old_j = order[j];
            result[i][j] = graph[std::min(old_i, old_j)]
                [std::max(old_i, old_j)];
        }
    }
    return result;
}

Graph swap_graph_vertices(
    const Graph& graph,
    std::size_t first,
    std::size_t second
) {
    std::vector<std::size_t> order(graph.size(), 0U);
    for (std::size_t i = 0U; i < graph.size(); ++i) {
        order[i] = i;
    }
    std::swap(order[first], order[second]);
    return permute_graph(graph, order);
}

Graph merge_graph_vertices(
    const Graph& graph,
    std::size_t retained,
    std::size_t removed
) {
    if (retained == removed || retained >= graph.size()
        || removed >= graph.size()) {
        throw std::runtime_error("invalid graph merge vertices");
    }
    Graph result(
        graph.size() - 1U,
        std::vector<int>(graph.size() - 1U, 0)
    );
    const auto compressed = [removed](std::size_t vertex) {
        return vertex < removed ? vertex : vertex - 1U;
    };
    for (std::size_t i = 0U; i < graph.size(); ++i) {
        if (i == removed) {
            continue;
        }
        for (std::size_t j = i + 1U; j < graph.size(); ++j) {
            if (j == removed) {
                continue;
            }
            result[compressed(i)][compressed(j)] += graph[i][j];
        }
    }
    const std::size_t retained_new = compressed(retained);
    for (std::size_t vertex = 0U; vertex < graph.size(); ++vertex) {
        if (vertex == retained || vertex == removed) {
            continue;
        }
        const std::size_t left = std::min(removed, vertex);
        const std::size_t right = std::max(removed, vertex);
        const std::size_t vertex_new = compressed(vertex);
        result[std::min(retained_new, vertex_new)]
              [std::max(retained_new, vertex_new)] += graph[left][right];
    }
    return result;
}

std::pair<std::size_t, std::size_t> extremal_collision_counts(
    const std::vector<Graph>& graphs,
    const std::vector<std::size_t>& order
) {
    std::map<Key, Support> memo;
    std::set<Key> least;
    std::set<Key> greatest;
    for (const Graph& graph : graphs) {
        const Graph permuted = permute_graph(graph, order);
        const Support& support = straighten_support(permuted, memo);
        if (support.empty()) {
            throw std::runtime_error("nonzero bracket monomial vanished");
        }
        least.insert(*support.begin());
        greatest.insert(*support.rbegin());
    }
    return {
        graphs.size() - least.size(),
        graphs.size() - greatest.size()
    };
}

std::size_t leading_collision_count(
    const std::vector<Graph>& graphs,
    const std::vector<std::size_t>& order
) {
    const auto [least_collisions, greatest_collisions]
        = extremal_collision_counts(graphs, order);
    return std::min(least_collisions, greatest_collisions);
}

std::pair<std::size_t, std::vector<std::size_t>>
best_leading_collision_count(const std::vector<Graph>& graphs) {
    if (graphs.empty()) {
        return {0U, {}};
    }
    const std::size_t vertices = graphs.front().size();
    std::vector<std::size_t> order(vertices, 0U);
    for (std::size_t i = 0U; i < vertices; ++i) {
        order[i] = i;
    }
    std::size_t best = graphs.size();
    std::vector<std::size_t> best_order = order;
    do {
        const std::size_t collisions = leading_collision_count(graphs, order);
        if (collisions < best) {
            best = collisions;
            best_order = order;
            if (best == 0U) {
                break;
            }
        }
    } while (std::next_permutation(order.begin() + 1, order.end()));
    return {best, best_order};
}

bool compatible_splits(
    unsigned int first,
    unsigned int second,
    unsigned int full_mask
) {
    const unsigned int first_complement = full_mask ^ first;
    const unsigned int second_complement = full_mask ^ second;
    return (first & second) == 0U
        || (first & second_complement) == 0U
        || (first_complement & second) == 0U
        || (first_complement & second_complement) == 0U;
}

int cut_weight(const Graph& graph, unsigned int split) {
    int weight = 0;
    for (std::size_t i = 0U; i < graph.size(); ++i) {
        for (std::size_t j = i + 1U; j < graph.size(); ++j) {
            if (((split >> i) & 1U) != ((split >> j) & 1U)) {
                weight += graph[i][j];
            }
        }
    }
    return weight;
}

std::pair<std::size_t, std::vector<unsigned int>>
best_tree_collision_count(const std::vector<Graph>& graphs) {
    if (graphs.empty()) {
        return {0U, {}};
    }
    const std::size_t vertices = graphs.front().size();
    if (vertices < 3U) {
        return {0U, {}};
    }
    const unsigned int full_mask = (1U << vertices) - 1U;
    std::vector<unsigned int> splits;
    for (unsigned int mask = 1U; mask < full_mask; ++mask) {
        if ((mask & 1U) == 0U || popcount(mask) < 2
            || popcount(full_mask ^ mask) < 2) {
            continue;
        }
        splits.push_back(mask);
    }
    std::size_t best = graphs.size();
    std::vector<unsigned int> best_tree;
    std::vector<unsigned int> chosen;
    const std::size_t target_size = vertices - 3U;
    const auto search = [&](const auto& self, std::size_t start) -> bool {
        if (chosen.size() == target_size) {
            std::set<std::vector<int>> signatures;
            for (const Graph& graph : graphs) {
                std::vector<int> signature;
                signature.reserve(chosen.size());
                for (const unsigned int split : chosen) {
                    signature.push_back(cut_weight(graph, split));
                }
                signatures.insert(std::move(signature));
            }
            const std::size_t collisions = graphs.size()
                - signatures.size();
            if (collisions < best) {
                best = collisions;
                best_tree = chosen;
            }
            return best == 0U;
        }
        for (std::size_t index = start; index < splits.size(); ++index) {
            bool compatible = true;
            for (const unsigned int split : chosen) {
                if (!compatible_splits(split, splits[index], full_mask)) {
                    compatible = false;
                    break;
                }
            }
            if (!compatible) {
                continue;
            }
            chosen.push_back(splits[index]);
            if (self(self, index + 1U)) {
                return true;
            }
            chosen.pop_back();
        }
        return false;
    };
    (void)search(search, 0U);
    return {best, best_tree};
}

std::vector<int> raw_leading_signature(
    const Graph& graph,
    const std::vector<std::size_t>& order
) {
    const std::size_t vertices = graph.size();
    if (order.size() != vertices) {
        throw std::runtime_error("raw leading order has wrong size");
    }
    std::vector<std::size_t> rank(vertices, 0U);
    for (std::size_t position = 0U; position < vertices; ++position) {
        rank[order[position]] = position;
    }
    std::vector<int> x_degree(vertices, 0);
    for (std::size_t i = 0U; i < vertices; ++i) {
        for (std::size_t j = i + 1U; j < vertices; ++j) {
            if (rank[i] < rank[j]) {
                x_degree[i] += graph[i][j];
            } else {
                x_degree[j] += graph[i][j];
            }
        }
    }
    return x_degree;
}

std::size_t raw_leading_collision_count(
    const std::vector<Graph>& graphs,
    const std::vector<std::size_t>& order
) {
    std::set<std::vector<int>> signatures;
    for (const Graph& graph : graphs) {
        signatures.insert(raw_leading_signature(graph, order));
    }
    return graphs.size() - signatures.size();
}

std::pair<std::size_t, std::vector<std::size_t>>
best_raw_leading_collision_count(const std::vector<Graph>& graphs) {
    if (graphs.empty()) {
        return {0U, {}};
    }
    const std::size_t vertices = graphs.front().size();
    std::vector<std::size_t> order(vertices, 0U);
    for (std::size_t i = 0U; i < vertices; ++i) {
        order[i] = i;
    }
    std::size_t best = graphs.size();
    std::vector<std::size_t> best_order = order;
    do {
        const std::size_t collisions = raw_leading_collision_count(
            graphs, order
        );
        if (collisions < best) {
            best = collisions;
            best_order = order;
            if (best == 0U) {
                break;
            }
        }
    } while (std::next_permutation(order.begin(), order.end()));
    return {best, best_order};
}

std::pair<std::size_t, std::vector<std::size_t>>
best_paired_raw_leading_collision_count(
    const std::vector<std::pair<Graph, Graph>>& graph_pairs
) {
    if (graph_pairs.empty()) {
        return {0U, {}};
    }
    const std::size_t vertices = graph_pairs.front().first.size();
    std::vector<std::size_t> order(vertices, 0U);
    for (std::size_t i = 0U; i < vertices; ++i) {
        order[i] = i;
    }
    std::size_t best = graph_pairs.size();
    std::vector<std::size_t> best_order = order;
    do {
        std::set<std::vector<int>> signatures;
        for (const auto& [first, second] : graph_pairs) {
            const std::vector<int> first_signature = raw_leading_signature(
                first, order
            );
            const std::vector<int> second_signature = raw_leading_signature(
                second, order
            );
            const auto less_in_term_order = [&order](
                const std::vector<int>& left,
                const std::vector<int>& right
            ) {
                for (std::size_t vertex : order) {
                    if (left[vertex] != right[vertex]) {
                        return left[vertex] < right[vertex];
                    }
                }
                return false;
            };
            signatures.insert(
                less_in_term_order(first_signature, second_signature)
                    ? second_signature : first_signature
            );
        }
        const std::size_t collisions = graph_pairs.size() - signatures.size();
        if (collisions < best) {
            best = collisions;
            best_order = order;
            if (best == 0U) {
                break;
            }
        }
    } while (std::next_permutation(order.begin(), order.end()));
    return {best, best_order};
}

using GraphFamilyKey = std::vector<Key>;

bool divisor_separates_graph_monomials_memo(
    const std::vector<Graph>& graphs,
    std::map<GraphFamilyKey, bool>& memo
) {
    GraphFamilyKey family_key;
    family_key.reserve(graphs.size());
    for (const Graph& graph : graphs) {
        family_key.push_back(flatten(graph));
    }
    std::sort(family_key.begin(), family_key.end());
    const auto memoized = memo.find(family_key);
    if (memoized != memo.end()) {
        return memoized->second;
    }
    if (graphs.size() <= 1U) {
        memo.emplace(std::move(family_key), true);
        return true;
    }
    std::set<Key> distinct;
    for (const Graph& graph : graphs) {
        distinct.insert(flatten(graph));
    }
    if (distinct.size() != graphs.size()) {
        memo.emplace(std::move(family_key), false);
        return false;
    }
    const std::size_t vertices = graphs.front().size();
    if (vertices <= 3U) {
        memo.emplace(std::move(family_key), true);
        return true;
    }
    for (std::size_t first = 0U; first < vertices; ++first) {
        for (std::size_t second = first + 1U; second < vertices; ++second) {
            std::map<int, std::vector<Graph>> groups;
            for (const Graph& graph : graphs) {
                Graph divided = graph;
                divided[first][second] = 0;
                groups[graph[first][second]].push_back(
                    merge_graph_vertices(divided, first, second)
                );
            }
            bool separates = true;
            for (const auto& [exponent, group] : groups) {
                (void)exponent;
                if (!divisor_separates_graph_monomials_memo(group, memo)) {
                    separates = false;
                    break;
                }
            }
            if (separates) {
                memo.emplace(std::move(family_key), true);
                return true;
            }
        }
    }
    memo.emplace(std::move(family_key), false);
    return false;
}

bool divisor_separates_graph_monomials(const std::vector<Graph>& graphs) {
    std::map<GraphFamilyKey, bool> memo;
    return divisor_separates_graph_monomials_memo(graphs, memo);
}

std::pair<bool, std::pair<std::size_t, std::size_t>>
first_divisor_separating_edge(const std::vector<Graph>& graphs) {
    if (graphs.size() <= 1U) {
        return {true, {0U, 0U}};
    }
    const std::size_t vertices = graphs.front().size();
    std::map<GraphFamilyKey, bool> memo;
    for (std::size_t first = 0U; first < vertices; ++first) {
        for (std::size_t second = first + 1U; second < vertices; ++second) {
            std::map<int, std::vector<Graph>> groups;
            for (const Graph& graph : graphs) {
                Graph divided = graph;
                divided[first][second] = 0;
                groups[graph[first][second]].push_back(
                    merge_graph_vertices(divided, first, second)
                );
            }
            bool separates = true;
            for (const auto& [exponent, group] : groups) {
                (void)exponent;
                if (!divisor_separates_graph_monomials_memo(group, memo)) {
                    separates = false;
                    break;
                }
            }
            if (separates) {
                return {true, {first, second}};
            }
        }
    }
    return {false, {0U, 0U}};
}

std::vector<Graph> noncrossing_graphs(
    const std::vector<int>& degrees,
    unsigned int mask
) {
    std::vector<int> restricted(degrees.size(), 0);
    for (std::size_t i = 0U; i < degrees.size(); ++i) {
        if (((mask >> i) & 1U) != 0U) {
            restricted[i] = degrees[i];
        }
    }
    return noncrossing_graphs(restricted);
}

std::pair<bool, std::array<std::size_t, 4>> first_crossing(
    const Graph& graph
) {
    const std::size_t vertices = graph.size();
    for (std::size_t a = 0U; a < vertices; ++a) {
        for (std::size_t b = a + 1U; b < vertices; ++b) {
            for (std::size_t c = b + 1U; c < vertices; ++c) {
                for (std::size_t d = c + 1U; d < vertices; ++d) {
                    if (graph[a][c] != 0 && graph[b][d] != 0) {
                        return {true, {a, b, c, d}};
                    }
                }
            }
        }
    }
    return {false, {0U, 0U, 0U, 0U}};
}

Graph separated_normal_form(Graph graph) {
    while (true) {
        const auto [has_crossing, crossing] = first_crossing(graph);
        if (!has_crossing) {
            return graph;
        }
        const auto [a, b, c, d] = crossing;
        --graph[a][c];
        --graph[b][d];
        ++graph[a][b];
        ++graph[c][d];
    }
}

void add_scaled(
    Expansion& target,
    const Expansion& source,
    std::int64_t scale
) {
    for (const auto& [key, coefficient] : source) {
        std::int64_t& value = target[key];
        value = normalized(value + scale * coefficient);
        if (value == 0) {
            target.erase(key);
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
    add_scaled(expansion, straighten(separated, memo), 1);
    add_scaled(expansion, straighten(nested, memo), 1);
    return memo.emplace(key, std::move(expansion)).first->second;
}

std::size_t expansion_rows_rank(const std::vector<Expansion>& rows) {
    std::map<Key, std::size_t> columns;
    for (const Expansion& row : rows) {
        for (const auto& [key, coefficient] : row) {
            (void)coefficient;
            columns.emplace(key, columns.size());
        }
    }
    std::vector<std::vector<std::int64_t>> matrix(
        rows.size(), std::vector<std::int64_t>(columns.size(), 0)
    );
    for (std::size_t row = 0U; row < rows.size(); ++row) {
        for (const auto& [key, coefficient] : rows[row]) {
            matrix[row][columns.at(key)] = coefficient;
        }
    }
    return modular_rank(std::move(matrix));
}

const Support& straighten_support(
    const Graph& graph,
    std::map<Key, Support>& memo
) {
    const Key key = flatten(graph);
    const auto found = memo.find(key);
    if (found != memo.end()) {
        return found->second;
    }
    const auto [has_crossing, crossing] = first_crossing(graph);
    if (!has_crossing) {
        Support support;
        support.insert(key);
        return memo.emplace(key, std::move(support)).first->second;
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
    Support support = straighten_support(separated, memo);
    const Support& nested_support = straighten_support(nested, memo);
    support.insert(nested_support.begin(), nested_support.end());
    return memo.emplace(key, std::move(support)).first->second;
}

bool subtract_degree(
    const std::vector<int>& degree,
    const std::vector<int>& base,
    std::vector<int>& result
) {
    result.resize(degree.size());
    for (std::size_t i = 0U; i < degree.size(); ++i) {
        if (degree[i] < base[i]) {
            return false;
        }
        result[i] = degree[i] - base[i];
    }
    return true;
}

std::pair<std::array<std::size_t, 3>, int> oriented_triple(
    std::size_t first,
    std::size_t second,
    std::size_t third
) {
    std::array<std::size_t, 3> values{first, second, third};
    int sign = 1;
    for (std::size_t i = 0U; i < 3U; ++i) {
        for (std::size_t j = i + 1U; j < 3U; ++j) {
            if (values[i] > values[j]) {
                std::swap(values[i], values[j]);
                sign = -sign;
            }
        }
    }
    return {values, sign};
}

Key module_key(
    const std::array<std::size_t, 3>& triple,
    const Graph& graph
) {
    Key key;
    key.reserve(3U + graph.size() * (graph.size() - 1U) / 2U);
    for (std::size_t index : triple) {
        key.push_back(static_cast<int>(index));
    }
    const Key graph_key = flatten(graph);
    key.insert(key.end(), graph_key.begin(), graph_key.end());
    return key;
}

Key module_key(std::size_t vertex, const Graph& graph) {
    Key key;
    key.reserve(1U + graph.size() * (graph.size() - 1U) / 2U);
    key.push_back(static_cast<int>(vertex));
    const Key graph_key = flatten(graph);
    key.insert(key.end(), graph_key.begin(), graph_key.end());
    return key;
}

struct F1BasisElement {
    std::array<std::size_t, 3> triple{};
    Graph coefficient;
};

void add_f1_term(
    SparseVector& vector,
    const Graph& coefficient,
    std::size_t edge_first,
    std::size_t edge_second,
    std::size_t triple_first,
    std::size_t triple_second,
    std::size_t triple_third,
    std::int64_t scalar,
    const std::map<Key, std::size_t>& f1_indices,
    std::map<Key, Expansion>& memo
) {
    if (edge_first > edge_second) {
        std::swap(edge_first, edge_second);
        scalar = -scalar;
    }
    const auto [triple, triple_sign] = oriented_triple(
        triple_first, triple_second, triple_third
    );
    scalar *= triple_sign;
    Graph product = coefficient;
    ++product[edge_first][edge_second];
    for (const auto& [graph_key, expansion_coefficient] :
         straighten(product, memo)) {
        const Graph standard = unflatten(graph_key, product.size());
        const auto found = f1_indices.find(module_key(triple, standard));
        if (found == f1_indices.end()) {
            throw std::runtime_error("second-syzygy term left F1 multidegree");
        }
        std::int64_t& value = vector[found->second];
        value = normalized(value + scalar * expansion_coefficient);
        if (value == 0) {
            vector.erase(found->second);
        }
    }
}

struct DegreeResult {
    bool pass = true;
    std::size_t f1_dimension = 0U;
    std::size_t kernel_dimension = 0U;
    std::size_t generator_count = 0U;
    std::size_t generator_rank = 0U;
};

void print_degree(const std::vector<int>& degree);
void print_graph(const Graph& graph);

DegreeResult analyze_degree(const std::vector<int>& degree) {
    const std::size_t vertices = degree.size();
    std::vector<F1BasisElement> f1_basis;
    std::map<Key, std::size_t> f1_indices;
    for (std::size_t i = 0U; i < vertices; ++i) {
        for (std::size_t j = i + 1U; j < vertices; ++j) {
            for (std::size_t k = j + 1U; k < vertices; ++k) {
                std::vector<int> base(vertices, 0);
                ++base[i];
                ++base[j];
                ++base[k];
                std::vector<int> residual;
                if (!subtract_degree(degree, base, residual)) {
                    continue;
                }
                for (Graph graph : noncrossing_graphs(residual)) {
                    const std::array<std::size_t, 3> triple{i, j, k};
                    const std::size_t index = f1_basis.size();
                    f1_indices.emplace(module_key(triple, graph), index);
                    f1_basis.push_back({triple, std::move(graph)});
                }
            }
        }
    }
    if (f1_basis.empty()) {
        return {};
    }

    std::map<Key, Expansion> memo;
    std::vector<std::map<Key, std::int64_t>> differential_sparse(
        f1_basis.size()
    );
    std::map<Key, std::size_t> f0_indices;
    const auto add_differential_term = [&memo, &f0_indices](
        std::map<Key, std::int64_t>& row,
        const Graph& coefficient,
        std::size_t edge_first,
        std::size_t edge_second,
        std::size_t output_vertex,
        std::int64_t scalar
    ) {
        Graph product = coefficient;
        ++product[edge_first][edge_second];
        for (const auto& [graph_key, expansion_coefficient] :
             straighten(product, memo)) {
            const Graph standard = unflatten(graph_key, product.size());
            const Key key = module_key(output_vertex, standard);
            f0_indices.emplace(key, f0_indices.size());
            std::int64_t& value = row[key];
            value = normalized(value + scalar * expansion_coefficient);
            if (value == 0) {
                row.erase(key);
            }
        }
    };
    for (std::size_t source = 0U; source < f1_basis.size(); ++source) {
        const auto [i, j, k] = f1_basis[source].triple;
        add_differential_term(
            differential_sparse[source], f1_basis[source].coefficient,
            j, k, i, 1
        );
        add_differential_term(
            differential_sparse[source], f1_basis[source].coefficient,
            i, k, j, -1
        );
        add_differential_term(
            differential_sparse[source], f1_basis[source].coefficient,
            i, j, k, 1
        );
    }
    std::vector<std::vector<std::int64_t>> differential(
        f1_basis.size(), std::vector<std::int64_t>(f0_indices.size(), 0)
    );
    for (std::size_t source = 0U; source < f1_basis.size(); ++source) {
        for (const auto& [key, coefficient] : differential_sparse[source]) {
            differential[source][f0_indices.at(key)] = coefficient;
        }
    }
    const std::size_t differential_rank = modular_rank(differential);
    const std::size_t kernel_dimension = f1_basis.size() - differential_rank;

    std::vector<SparseVector> generators;
    for (std::size_t anchor = 0U; anchor < vertices; ++anchor) {
        std::vector<std::size_t> others;
        for (std::size_t vertex = 0U; vertex < vertices; ++vertex) {
            if (vertex != anchor) {
                others.push_back(vertex);
            }
        }
        for (std::size_t x = 0U; x < others.size(); ++x) {
            for (std::size_t y = x + 1U; y < others.size(); ++y) {
                for (std::size_t z = y + 1U; z < others.size(); ++z) {
                    const std::size_t b = others[x];
                    const std::size_t c = others[y];
                    const std::size_t d = others[z];
                    std::vector<int> base(vertices, 0);
                    base[anchor] = 2;
                    ++base[b];
                    ++base[c];
                    ++base[d];
                    std::vector<int> residual;
                    if (!subtract_degree(degree, base, residual)) {
                        continue;
                    }
                    for (const Graph& coefficient : noncrossing_graphs(residual)) {
                        SparseVector generator;
                        add_f1_term(
                            generator, coefficient, anchor, b,
                            anchor, c, d, 1, f1_indices, memo
                        );
                        add_f1_term(
                            generator, coefficient, anchor, c,
                            anchor, b, d, -1, f1_indices, memo
                        );
                        add_f1_term(
                            generator, coefficient, anchor, d,
                            anchor, b, c, 1, f1_indices, memo
                        );
                        generators.push_back(std::move(generator));
                    }
                }
            }
        }
    }
    if (vertices >= 5U) {
        for (std::size_t a = 0U; a < vertices; ++a) {
            for (std::size_t b = a + 1U; b < vertices; ++b) {
                for (std::size_t c = b + 1U; c < vertices; ++c) {
                    for (std::size_t d = c + 1U; d < vertices; ++d) {
                        for (std::size_t e = d + 1U; e < vertices; ++e) {
                            const std::array<std::size_t, 5> set{a, b, c, d, e};
                            std::vector<int> base(vertices, 0);
                            for (std::size_t vertex : set) {
                                ++base[vertex];
                            }
                            std::vector<int> residual;
                            if (!subtract_degree(degree, base, residual)) {
                                continue;
                            }
                            for (const Graph& coefficient :
                                 noncrossing_graphs(residual)) {
                                for (std::size_t distinguished = 0U;
                                     distinguished < 5U; ++distinguished) {
                                    SparseVector generator;
                                    for (std::size_t x = 0U; x < 5U; ++x) {
                                        if (x == distinguished) {
                                            continue;
                                        }
                                        std::array<std::size_t, 3> remaining{};
                                        std::size_t position = 0U;
                                        for (std::size_t z = 0U; z < 5U; ++z) {
                                            if (z != x && z != distinguished) {
                                                remaining[position] = set[z];
                                                ++position;
                                            }
                                        }
                                        const int sign = (x % 2U == 0U) ? 1 : -1;
                                        add_f1_term(
                                            generator, coefficient,
                                            std::min(set[x], set[distinguished]),
                                            std::max(set[x], set[distinguished]),
                                            remaining[0],
                                            remaining[1], remaining[2], sign,
                                            f1_indices, memo
                                        );
                                    }
                                    generators.push_back(std::move(generator));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    bool composition_zero = true;
    for (const SparseVector& generator : generators) {
        std::vector<std::int64_t> image(f0_indices.size(), 0);
        for (const auto& [source, coefficient] : generator) {
            for (std::size_t target = 0U; target < f0_indices.size(); ++target) {
                image[target] = normalized(
                    image[target] + coefficient * differential[source][target]
                );
            }
        }
        if (std::any_of(
                image.begin(), image.end(),
                [](std::int64_t value) { return value != 0; }
            )) {
            composition_zero = false;
            break;
        }
    }
    std::vector<std::vector<std::int64_t>> generator_matrix(
        generators.size(), std::vector<std::int64_t>(f1_basis.size(), 0)
    );
    for (std::size_t row = 0U; row < generators.size(); ++row) {
        for (const auto& [column, coefficient] : generators[row]) {
            generator_matrix[row][column] = coefficient;
        }
    }
    const std::size_t generator_rank = modular_rank(std::move(generator_matrix));
    if (generator_rank != kernel_dimension
        && std::getenv("PLUCKER_VERBOSE") != nullptr) {
        const auto kernel = modular_left_kernel(differential);
#pragma omp critical
        {
            std::cout << "KERNEL_BASIS degree=";
            print_degree(degree);
            std::cout << '\n';
            for (std::size_t relation = 0U; relation < kernel.size();
                 ++relation) {
                std::cout << "relation " << relation << ':';
                for (std::size_t source = 0U; source < kernel[relation].size();
                     ++source) {
                    std::int64_t coefficient = kernel[relation][source];
                    if (coefficient == 0) {
                        continue;
                    }
                    if (coefficient > modulus / 2) {
                        coefficient -= modulus;
                    }
                    const auto [i, j, k] = f1_basis[source].triple;
                    std::cout << " (" << coefficient << ")p";
                    bool printed_edge = false;
                    for (std::size_t u = 0U; u < vertices; ++u) {
                        for (std::size_t v = u + 1U; v < vertices; ++v) {
                            if (f1_basis[source].coefficient[u][v] != 0) {
                                std::cout << u << v;
                                printed_edge = true;
                            }
                        }
                    }
                    if (!printed_edge) {
                        std::cout << "1";
                    }
                    std::cout << "e" << i << j << k;
                }
                std::cout << '\n';
            }
        }
    }
    return DegreeResult{
        composition_zero && generator_rank == kernel_dimension,
        f1_basis.size(), kernel_dimension, generators.size(), generator_rank
    };
}

struct F0BasisElement {
    std::size_t vertex = 0U;
    Graph coefficient;
};

std::size_t compressed_vertex(std::size_t vertex, std::size_t removed) {
    if (vertex == removed) {
        throw std::runtime_error("cannot compress removed vertex");
    }
    return vertex < removed ? vertex : vertex - 1U;
}

std::size_t expanded_vertex(std::size_t vertex, std::size_t inserted) {
    return vertex < inserted ? vertex : vertex + 1U;
}

Graph strip_vertex(const Graph& graph, std::size_t removed) {
    Graph result(
        graph.size() - 1U,
        std::vector<int>(graph.size() - 1U, 0)
    );
    for (std::size_t i = 0U; i < graph.size(); ++i) {
        if (i == removed) {
            continue;
        }
        for (std::size_t j = i + 1U; j < graph.size(); ++j) {
            if (j == removed) {
                continue;
            }
            result[compressed_vertex(i, removed)][compressed_vertex(j, removed)] =
                graph[i][j];
        }
    }
    return result;
}

Graph insert_vertex_and_edge(
    const Graph& coefficient,
    std::size_t inserted,
    std::size_t neighbor
) {
    Graph result(
        coefficient.size() + 1U,
        std::vector<int>(coefficient.size() + 1U, 0)
    );
    for (std::size_t i = 0U; i < coefficient.size(); ++i) {
        for (std::size_t j = i + 1U; j < coefficient.size(); ++j) {
            result[expanded_vertex(i, inserted)][expanded_vertex(j, inserted)] =
                coefficient[i][j];
        }
    }
    const std::size_t old_neighbor = expanded_vertex(neighbor, inserted);
    ++result[std::min(inserted, old_neighbor)][std::max(inserted, old_neighbor)];
    return result;
}

struct TauEvenResult {
    std::size_t domain_dimension = 0U;
    std::size_t global_kernel = 0U;
    std::size_t iota_kernel = 0U;
    std::size_t carrier_kernel = 0U;
    std::int64_t fundamental_pair_reservoir = 0;
    std::size_t merged_raw_collisions = 0U;
    std::size_t merged_standard_collisions = 0U;
    std::vector<std::size_t> merged_standard_order;
    bool divisor_separated = false;
    std::array<std::size_t, 2> sigma_domain{};
    std::array<std::size_t, 2> sigma_global_kernel{};
    std::array<std::size_t, 2> sigma_iota_kernel{};
    bool sigma_valid = false;
    bool commuting_rank = false;
    bool separation_diagnostics_valid = false;
};

std::vector<Graph> tau_even_representatives(
    const std::vector<int>& labels
) {
    constexpr std::size_t a_vertex = 2U;
    constexpr std::size_t f_vertex = 3U;
    const std::size_t full_vertices = labels.size();
    if (full_vertices < 5U
        || full_vertices >= std::numeric_limits<unsigned int>::digits
        || labels[a_vertex] != 1
        || labels[f_vertex] != 1) {
        throw std::runtime_error("invalid tau-even allocation labels");
    }
    const unsigned int full_mask = (1U << full_vertices) - 1U;
    std::vector<Graph> representatives;
    for (std::size_t side = 0U; side < 2U; ++side) {
        for (std::size_t even = 4U; even < full_vertices; ++even) {
            const unsigned int mask = (1U << side) | (1U << a_vertex)
                | (1U << even);
            const auto triple = noncrossing_graphs(labels, mask);
            const auto complement = noncrossing_graphs(
                labels, full_mask ^ mask
            );
            for (const Graph& left : triple) {
                for (const Graph& right : complement) {
                    representatives.push_back(graph_sum(left, right));
                }
            }
        }
    }
    return representatives;
}

std::vector<Graph> tau_even_merged_graphs(const std::vector<int>& labels) {
    constexpr std::size_t a_vertex = 2U;
    constexpr std::size_t f_vertex = 3U;
    const std::vector<Graph> representatives
        = tau_even_representatives(labels);
    std::vector<Graph> merged_graphs;
    merged_graphs.reserve(representatives.size());
    for (const Graph& representative : representatives) {
        merged_graphs.push_back(merge_graph_vertices(
            representative, a_vertex, f_vertex
        ));
    }
    return merged_graphs;
}

struct EqualTerminalDiagonalResult {
    std::size_t plus_domain = 0U;
    std::size_t plus_rank = 0U;
    std::size_t plus_terminal_tensor_rank = 0U;
    bool plus_diagonal_v2 = true;
    std::size_t minus_domain = 0U;
    std::size_t minus_rank = 0U;
    std::size_t minus_terminal_tensor_rank = 0U;
    std::size_t minus_wronskian_rank = 0U;
    bool antisymmetric_divisible = true;
};

EqualTerminalDiagonalResult analyze_equal_terminal_diagonal(
    const std::vector<int>& labels
) {
    if (labels.size() != 7U || labels[0U] != labels[1U]
        || labels[2U] != 1 || labels[3U] != 1) {
        throw std::runtime_error("invalid equal-terminal diagonal labels");
    }
    const std::vector<Graph> graphs = tau_even_merged_graphs(labels);
    std::map<Key, std::size_t> lookup;
    for (std::size_t source = 0U; source < graphs.size(); ++source) {
        if (!lookup.emplace(flatten(graphs[source]), source).second) {
            throw std::runtime_error("duplicate equal-terminal source graph");
        }
    }

    std::map<Key, Expansion> six_vertex_memo;
    std::map<Key, Expansion> five_vertex_memo;
    std::vector<Expansion> plus_rows;
    std::vector<Expansion> minus_rows;
    std::vector<std::vector<std::int64_t>> plus_terminal_rows;
    std::vector<std::vector<std::int64_t>> minus_terminal_rows;
    std::vector<std::vector<std::int64_t>> minus_wronskian_rows;
    std::vector<std::vector<int>> plus_diagonal_exponents;
    bool divisible = true;
    for (std::size_t source = 0U; source < graphs.size(); ++source) {
        const Graph swapped = swap_graph_vertices(graphs[source], 0U, 1U);
        const auto found = lookup.find(flatten(swapped));
        if (found == lookup.end()) {
            throw std::runtime_error("equal-terminal source orbit is incomplete");
        }
        const std::size_t partner = found->second;
        if (source > partner) {
            continue;
        }

        const Graph plus_diagonal = merge_graph_vertices(
            graphs[source], 0U, 1U
        );
        plus_rows.push_back(straighten(plus_diagonal, five_vertex_memo));
        std::vector<int> diagonal_exponents;
        diagonal_exponents.reserve(4U);
        for (std::size_t root = 1U; root < plus_diagonal.size(); ++root) {
            diagonal_exponents.push_back(plus_diagonal[0U][root]);
        }
        plus_diagonal_exponents.push_back(std::move(diagonal_exponents));

        const std::vector<std::int64_t> first_terminal
            = binary_vertex_factor_coefficients(graphs[source], 0U);
        const std::vector<std::int64_t> second_terminal
            = binary_vertex_factor_coefficients(graphs[source], 1U);
        if (first_terminal.size() != second_terminal.size()) {
            throw std::runtime_error("equal-terminal factor degrees differ");
        }
        std::vector<std::int64_t> plus_terminal(
            first_terminal.size() * second_terminal.size(), 0
        );
        std::vector<std::int64_t> minus_terminal(
            first_terminal.size() * second_terminal.size(), 0
        );
        for (std::size_t first = 0U; first < first_terminal.size(); ++first) {
            for (std::size_t second = 0U;
                 second < second_terminal.size(); ++second) {
                const std::size_t position
                    = first * second_terminal.size() + second;
                const std::int64_t forward = normalized(
                    first_terminal[first] * second_terminal[second]
                );
                const std::int64_t backward = normalized(
                    second_terminal[first] * first_terminal[second]
                );
                plus_terminal[position] = normalized(forward + backward);
                minus_terminal[position] = normalized(forward - backward);
            }
        }
        plus_terminal_rows.push_back(std::move(plus_terminal));

        if (source == partner) {
            continue;
        }
        minus_terminal_rows.push_back(std::move(minus_terminal));
        std::vector<std::int64_t> wronskian(
            2U * first_terminal.size() - 3U, 0
        );
        for (std::size_t first = 0U; first < first_terminal.size(); ++first) {
            for (std::size_t second = 0U;
                 second < second_terminal.size(); ++second) {
                if (first + second == 0U) {
                    continue;
                }
                const std::int64_t difference
                    = static_cast<std::int64_t>(second)
                    - static_cast<std::int64_t>(first);
                const std::int64_t product = normalized(
                    first_terminal[first] * second_terminal[second]
                );
                const std::size_t power = first + second - 1U;
                wronskian[power] = normalized(
                    wronskian[power]
                    + normalized(difference) * product
                );
            }
        }
        minus_wronskian_rows.push_back(std::move(wronskian));
        Expansion antisymmetric;
        add_scaled(
            antisymmetric,
            straighten(graphs[source], six_vertex_memo), 1
        );
        add_scaled(
            antisymmetric,
            straighten(graphs[partner], six_vertex_memo), -1
        );
        Expansion divided_diagonal;
        for (const auto& [key, coefficient] : antisymmetric) {
            Graph quotient = unflatten(key, 6U);
            if (quotient[0U][1U] == 0) {
                divisible = false;
                continue;
            }
            --quotient[0U][1U];
            if (quotient[0U][1U] != 0) {
                // After division by one p_(qr), every term with another
                // p_(qr) factor vanishes on the terminal diagonal.
                continue;
            }
            const Graph diagonal = merge_graph_vertices(
                quotient, 0U, 1U
            );
            add_scaled(
                divided_diagonal,
                straighten(diagonal, five_vertex_memo), coefficient
            );
        }
        Expansion expected_wronskian;
        std::array<int, 4U> first_exponents{};
        std::array<int, 4U> second_exponents{};
        Graph nonterminal(5U, std::vector<int>(5U, 0));
        for (std::size_t root = 0U; root < 4U; ++root) {
            first_exponents[root] = graphs[source][0U][root + 2U];
            second_exponents[root] = graphs[source][1U][root + 2U];
        }
        for (std::size_t first_root = 0U; first_root < 4U; ++first_root) {
            for (std::size_t second_root = first_root + 1U;
                 second_root < 4U; ++second_root) {
                nonterminal[first_root + 1U][second_root + 1U]
                    = graphs[source][first_root + 2U][second_root + 2U];
            }
        }
        for (std::size_t first_root = 0U; first_root < 4U; ++first_root) {
            for (std::size_t second_root = 0U; second_root < 4U;
                 ++second_root) {
                if (first_root == second_root
                    || first_exponents[first_root] == 0
                    || second_exponents[second_root] == 0) {
                    continue;
                }
                Graph term = nonterminal;
                for (std::size_t root = 0U; root < 4U; ++root) {
                    term[0U][root + 1U]
                        = first_exponents[root] + second_exponents[root]
                        - static_cast<int>(root == first_root)
                        - static_cast<int>(root == second_root);
                }
                ++term[std::min(first_root, second_root) + 1U]
                      [std::max(first_root, second_root) + 1U];
                const std::int64_t orientation
                    = first_root < second_root ? 1 : -1;
                const std::int64_t scale = orientation
                    * static_cast<std::int64_t>(first_exponents[first_root])
                    * static_cast<std::int64_t>(second_exponents[second_root]);
                add_scaled(
                    expected_wronskian,
                    straighten(term, five_vertex_memo), scale
                );
            }
        }
        Expansion negated_expected;
        add_scaled(negated_expected, expected_wronskian, -1);
        Expansion scaled_divided;
        add_scaled(
            scaled_divided, divided_diagonal,
            static_cast<std::int64_t>(labels[0U])
        );
        if (scaled_divided != expected_wronskian
            && scaled_divided != negated_expected) {
            std::cerr << "Wronskian mismatch divided_terms="
                      << divided_diagonal.size()
                      << " expected_terms=" << expected_wronskian.size()
                      << '\n';
            std::set<Key> mismatch_keys;
            for (const auto& [key, coefficient] : divided_diagonal) {
                (void)coefficient;
                mismatch_keys.insert(key);
            }
            for (const auto& [key, coefficient] : expected_wronskian) {
                (void)coefficient;
                mismatch_keys.insert(key);
            }
            std::size_t printed = 0U;
            for (const Key& key : mismatch_keys) {
                const auto divided_found = divided_diagonal.find(key);
                const auto expected_found = expected_wronskian.find(key);
                const std::int64_t divided_coefficient
                    = divided_found == divided_diagonal.end()
                    ? 0 : divided_found->second;
                const std::int64_t expected_coefficient
                    = expected_found == expected_wronskian.end()
                    ? 0 : expected_found->second;
                if (normalized(
                        static_cast<std::int64_t>(labels[0U])
                            * divided_coefficient
                    ) == expected_coefficient) {
                    continue;
                }
                std::cerr << "  divided=" << divided_coefficient
                          << " expected=" << expected_coefficient
                          << " key=";
                for (int entry : key) {
                    std::cerr << entry << ',';
                }
                std::cerr << '\n';
                ++printed;
                if (printed == 10U) {
                    break;
                }
            }
            throw std::runtime_error(
                "divided antisymmetric row disagrees with Wronskian identity"
            );
        }
        minus_rows.push_back(std::move(divided_diagonal));
    }
    const bool plus_diagonal_v2
        = (exponent_rows_v2_shape(plus_diagonal_exponents)
            && exponent_rows_gcd_hall_condition(plus_diagonal_exponents))
        || exponent_rows_v2_merge_certificate(plus_diagonal_exponents);
    return EqualTerminalDiagonalResult{
        plus_rows.size(), expansion_rows_rank(plus_rows),
        modular_rank(std::move(plus_terminal_rows)),
        plus_diagonal_v2,
        minus_rows.size(), expansion_rows_rank(minus_rows),
        modular_rank(std::move(minus_terminal_rows)),
        modular_rank(std::move(minus_wronskian_rows)), divisible
    };
}

using TerminalExponentPair = std::array<int, 8U>;

struct EqualTerminalFactorResult {
    std::size_t plus_domain = 0U;
    std::size_t plus_tensor_rank = 0U;
    std::size_t minus_domain = 0U;
    std::size_t minus_wronskian_rank = 0U;
    bool strict_min_v2 = false;
    std::size_t component_forms_domain = 0U;
    std::size_t component_forms_rank = 0U;
    bool component_forms_v2 = false;
    bool unordered_pairs_distinct = false;
    bool plus_leading_v2 = false;
    bool minus_leading_v2 = false;
};

EqualTerminalFactorResult analyze_equal_terminal_factors(
    const std::vector<int>& labels
) {
    if (labels.size() != 7U || labels[0U] != labels[1U]
        || labels[2U] != 1 || labels[3U] != 1) {
        throw std::runtime_error("invalid equal-terminal factor labels");
    }
    const std::vector<Graph> graphs = tau_even_merged_graphs(labels);
    std::map<Key, std::size_t> lookup;
    for (std::size_t source = 0U; source < graphs.size(); ++source) {
        if (!lookup.emplace(flatten(graphs[source]), source).second) {
            throw std::runtime_error("duplicate equal-terminal factor source");
        }
    }
    std::vector<std::vector<std::int64_t>> plus_rows;
    std::vector<std::vector<std::int64_t>> wronskian_rows;
    std::map<std::array<int, 4U>, std::vector<std::int64_t>>
        component_coefficient_rows;
    std::vector<std::pair<std::array<int, 4U>, std::array<int, 4U>>>
        exponent_pairs;
    std::vector<bool> antisymmetric_orbits;
    for (std::size_t source = 0U; source < graphs.size(); ++source) {
        const Graph swapped = swap_graph_vertices(graphs[source], 0U, 1U);
        const auto found = lookup.find(flatten(swapped));
        if (found == lookup.end()) {
            throw std::runtime_error("incomplete equal-terminal factor orbit");
        }
        const std::size_t partner = found->second;
        if (source > partner) {
            continue;
        }
        const std::vector<std::int64_t> first_terminal
            = binary_vertex_factor_coefficients(graphs[source], 0U);
        const std::vector<std::int64_t> second_terminal
            = binary_vertex_factor_coefficients(graphs[source], 1U);
        std::array<int, 4U> first_exponents{};
        std::array<int, 4U> second_exponents{};
        for (std::size_t root = 0U; root < 4U; ++root) {
            first_exponents[root] = graphs[source][0U][root + 2U];
            second_exponents[root] = graphs[source][1U][root + 2U];
        }
        exponent_pairs.push_back({first_exponents, second_exponents});
        antisymmetric_orbits.push_back(source != partner);
        component_coefficient_rows.emplace(
            first_exponents, first_terminal
        );
        component_coefficient_rows.emplace(
            second_exponents, second_terminal
        );
        std::vector<std::int64_t> plus(
            first_terminal.size() * second_terminal.size(), 0
        );
        for (std::size_t first = 0U; first < first_terminal.size(); ++first) {
            for (std::size_t second = 0U;
                 second < second_terminal.size(); ++second) {
                const std::size_t position
                    = first * second_terminal.size() + second;
                plus[position] = normalized(
                    first_terminal[first] * second_terminal[second]
                    + second_terminal[first] * first_terminal[second]
                );
            }
        }
        plus_rows.push_back(std::move(plus));
        if (source == partner) {
            continue;
        }
        std::vector<std::int64_t> wronskian(
            2U * first_terminal.size() - 3U, 0
        );
        for (std::size_t first = 0U; first < first_terminal.size(); ++first) {
            for (std::size_t second = 0U;
                 second < second_terminal.size(); ++second) {
                if (first + second == 0U) {
                    continue;
                }
                const std::int64_t difference
                    = static_cast<std::int64_t>(second)
                    - static_cast<std::int64_t>(first);
                const std::int64_t product = normalized(
                    first_terminal[first] * second_terminal[second]
                );
                const std::size_t power = first + second - 1U;
                wronskian[power] = normalized(
                    wronskian[power]
                    + normalized(difference) * product
                );
            }
        }
        wronskian_rows.push_back(std::move(wronskian));
    }
    std::set<std::array<int, 4U>> component_forms;
    std::set<std::pair<std::array<int, 4U>, std::array<int, 4U>>>
        unordered_pairs;
    for (const auto& [first, second] : exponent_pairs) {
        component_forms.insert(first);
        component_forms.insert(second);
        unordered_pairs.insert(first <= second
            ? std::make_pair(first, second)
            : std::make_pair(second, first));
    }
    std::vector<std::vector<int>> component_rows;
    component_rows.reserve(component_forms.size());
    for (const std::array<int, 4U>& row : component_forms) {
        component_rows.emplace_back(row.begin(), row.end());
    }
    const bool component_forms_v2
        = (exponent_rows_v2_shape(component_rows)
            && exponent_rows_gcd_hall_condition(component_rows))
        || exponent_rows_v2_merge_certificate(component_rows);
    const bool unordered_pairs_distinct
        = unordered_pairs.size() == exponent_pairs.size();
    std::vector<std::vector<std::int64_t>> component_coefficients;
    component_coefficients.reserve(component_coefficient_rows.size());
    for (const auto& [exponents, coefficients] : component_coefficient_rows) {
        (void)exponents;
        component_coefficients.push_back(coefficients);
    }
    const std::size_t component_forms_rank
        = modular_rank(std::move(component_coefficients));

    const auto leading_v2_certificate = [&](bool antisymmetric) {
        for (std::size_t root = 0U; root < 4U; ++root) {
            using SupportRow = std::vector<std::array<int, 4U>>;
            std::map<int, std::vector<SupportRow>> groups;
            for (std::size_t orbit = 0U; orbit < exponent_pairs.size();
                 ++orbit) {
                if (antisymmetric && !antisymmetric_orbits[orbit]) {
                    continue;
                }
                const auto& [first, second] = exponent_pairs[orbit];
                SupportRow support;
                if (first[root] < second[root]) {
                    support.push_back(second);
                } else if (second[root] < first[root]) {
                    support.push_back(first);
                } else if (first == second) {
                    if (!antisymmetric) {
                        support.push_back(first);
                    }
                } else {
                    support.push_back(first);
                    support.push_back(second);
                }
                if (!support.empty()) {
                    groups[std::min(first[root], second[root])].push_back(
                        std::move(support)
                    );
                }
            }
            bool root_passes = true;
            for (const auto& [order, support_rows] : groups) {
                (void)order;
                std::set<std::array<int, 4U>> forms;
                for (const SupportRow& support : support_rows) {
                    forms.insert(support.begin(), support.end());
                }
                std::vector<std::vector<int>> rows;
                rows.reserve(forms.size());
                for (const std::array<int, 4U>& form : forms) {
                    rows.emplace_back(form.begin(), form.end());
                }
                const bool forms_independent
                    = (exponent_rows_v2_shape(rows)
                        && exponent_rows_gcd_hall_condition(rows))
                    || exponent_rows_v2_merge_certificate(rows);
                if (!forms_independent) {
                    root_passes = false;
                    break;
                }
                std::vector<bool> live(support_rows.size(), true);
                std::size_t remaining = support_rows.size();
                while (remaining != 0U) {
                    std::map<std::array<int, 4U>, std::size_t> counts;
                    for (std::size_t row = 0U; row < support_rows.size();
                         ++row) {
                        if (!live[row]) {
                            continue;
                        }
                        for (const std::array<int, 4U>& form
                             : support_rows[row]) {
                            ++counts[form];
                        }
                    }
                    std::size_t removable = support_rows.size();
                    for (std::size_t row = 0U; row < support_rows.size();
                         ++row) {
                        if (!live[row]) {
                            continue;
                        }
                        for (const std::array<int, 4U>& form
                             : support_rows[row]) {
                            if (counts[form] == 1U) {
                                removable = row;
                                break;
                            }
                        }
                        if (removable != support_rows.size()) {
                            break;
                        }
                    }
                    if (removable == support_rows.size()) {
                        root_passes = false;
                        break;
                    }
                    live[removable] = false;
                    --remaining;
                }
                if (!root_passes) {
                    break;
                }
            }
            if (root_passes) {
                return true;
            }
        }
        return false;
    };
    const bool plus_leading_v2 = leading_v2_certificate(false);
    const bool minus_leading_v2 = leading_v2_certificate(true);

    bool strict_min_v2 = false;
    for (std::size_t jet_root = 0U; jet_root < 4U && !strict_min_v2;
         ++jet_root) {
        bool strict = true;
        std::map<int, std::vector<std::vector<int>>> groups;
        for (const auto& [first, second] : exponent_pairs) {
            if (first[jet_root] == second[jet_root]) {
                strict = false;
                break;
            }
            const std::array<int, 4U>& leading_form
                = first[jet_root] < second[jet_root] ? second : first;
            groups[std::min(first[jet_root], second[jet_root])].push_back(
                std::vector<int>(leading_form.begin(), leading_form.end())
            );
        }
        if (!strict) {
            continue;
        }
        bool groups_pass = true;
        for (const auto& [order, rows] : groups) {
            (void)order;
            const bool direct = exponent_rows_v2_shape(rows)
                && exponent_rows_gcd_hall_condition(rows);
            if (!direct && !exponent_rows_v2_merge_certificate(rows)) {
                groups_pass = false;
                break;
            }
        }
        strict_min_v2 = groups_pass;
    }
    return EqualTerminalFactorResult{
        plus_rows.size(), modular_rank(std::move(plus_rows)),
        wronskian_rows.size(), modular_rank(std::move(wronskian_rows)),
        strict_min_v2, component_coefficient_rows.size(),
        component_forms_rank, component_forms_v2, unordered_pairs_distinct,
        plus_leading_v2, minus_leading_v2
    };
}

std::multiset<TerminalExponentPair> explicit_terminal_source_rows(
    const std::vector<int>& labels
) {
    if (labels.size() != 7U || labels[2U] != 1 || labels[3U] != 1) {
        throw std::runtime_error("invalid explicit terminal-source labels");
    }
    const std::array<int, 2U> terminal{labels[0U], labels[1U]};
    const std::array<int, 3U> even{labels[4U], labels[5U], labels[6U]};
    std::multiset<TerminalExponentPair> rows;
    for (std::size_t selected_terminal = 0U;
         selected_terminal < terminal.size(); ++selected_terminal) {
        const int selected_degree = terminal[selected_terminal];
        const int complement_degree = terminal[1U - selected_terminal];
        for (std::size_t selected = 0U; selected < even.size(); ++selected) {
            std::array<int, 4U> selected_exponents{};
            if (even[selected] == selected_degree - 1) {
                selected_exponents[0U] = 1;
                selected_exponents[selected + 1U] = selected_degree - 1;
            } else if (even[selected] == selected_degree + 1) {
                selected_exponents[selected + 1U] = selected_degree;
            } else {
                continue;
            }
            std::array<std::size_t, 2U> remaining{};
            std::size_t position = 0U;
            for (std::size_t index = 0U; index < even.size(); ++index) {
                if (index != selected) {
                    remaining[position] = index;
                    ++position;
                }
            }
            const int x = even[remaining[0U]];
            const int y = even[remaining[1U]];
            for (std::size_t channel = 0U; channel < 2U; ++channel) {
                std::array<int, 4U> complement_exponents{};
                bool active = false;
                if (channel == 0U) {
                    active = y - x <= complement_degree - 1
                        && x + y >= complement_degree - 1;
                    if (active) {
                        complement_exponents[0U] = 1;
                        complement_exponents[remaining[0U] + 1U]
                            = (complement_degree - 1 + x - y) / 2;
                        complement_exponents[remaining[1U] + 1U]
                            = (complement_degree - 1 - x + y) / 2;
                    }
                } else if (y - x == complement_degree + 1) {
                    active = true;
                    complement_exponents[remaining[1U] + 1U]
                        = complement_degree;
                } else {
                    active = y - x <= complement_degree - 1
                        && x + y >= complement_degree + 1;
                    if (active) {
                        complement_exponents[remaining[0U] + 1U]
                            = (complement_degree + x - y - 1) / 2;
                        complement_exponents[remaining[1U] + 1U]
                            = (complement_degree - x + y + 1) / 2;
                    }
                }
                if (!active) {
                    continue;
                }
                TerminalExponentPair row{};
                for (std::size_t coordinate = 0U; coordinate < 4U;
                     ++coordinate) {
                    row[(selected_terminal == 0U ? 0U : 4U) + coordinate]
                        = selected_exponents[coordinate];
                    row[(selected_terminal == 0U ? 4U : 0U) + coordinate]
                        = complement_exponents[coordinate];
                }
                rows.insert(row);
            }
        }
    }
    return rows;
}

std::multiset<TerminalExponentPair> generated_terminal_source_rows(
    const std::vector<int>& labels
) {
    const std::vector<Graph> graphs = tau_even_merged_graphs(labels);
    std::multiset<TerminalExponentPair> rows;
    for (const Graph& graph : graphs) {
        TerminalExponentPair row{};
        if (graph[0U][1U] != 0) {
            throw std::runtime_error("unexpected terminal-terminal source edge");
        }
        for (std::size_t root = 2U; root < 6U; ++root) {
            row[root - 2U] = graph[0U][root];
            row[4U + root - 2U] = graph[1U][root];
        }
        rows.insert(row);
    }
    return rows;
}

TauEvenResult analyze_tau_even_allocation(
    const std::vector<int>& labels,
    bool include_separation_diagnostics = true
) {
    constexpr std::size_t a_vertex = 2U;
    constexpr std::size_t f_vertex = 3U;
    const std::size_t full_vertices = labels.size();
    if (full_vertices < 5U
        || full_vertices >= std::numeric_limits<unsigned int>::digits
        || labels[a_vertex] != 1
        || labels[f_vertex] != 1) {
        throw std::runtime_error("invalid tau-even allocation labels");
    }
    std::vector<Graph> representatives = tau_even_representatives(labels);

    std::map<Key, Expansion> global_memo;
    std::map<Key, Expansion> reduced_memo;
    std::vector<Expansion> global_rows;
    std::vector<Expansion> iota_rows;
    std::map<Key, std::size_t> global_columns;
    std::map<Key, std::size_t> iota_columns;
    for (const Graph& representative : representatives) {
        const Graph partner = swap_graph_vertices(
            representative, a_vertex, f_vertex
        );
        Expansion global_row;
        add_scaled(
            global_row, straighten(representative, global_memo), 1
        );
        add_scaled(global_row, straighten(partner, global_memo), 1);
        for (const auto& [key, coefficient] : global_row) {
            (void)coefficient;
            global_columns.emplace(key, global_columns.size());
        }
        global_rows.push_back(std::move(global_row));

        Expansion iota_row;
        for (const Graph* component : {&representative, &partner}) {
            std::size_t neighbor = full_vertices;
            for (std::size_t vertex = 0U; vertex < full_vertices; ++vertex) {
                if (vertex == f_vertex) {
                    continue;
                }
                const std::size_t left = std::min(vertex, f_vertex);
                const std::size_t right = std::max(vertex, f_vertex);
                if ((*component)[left][right] != 0) {
                    if ((*component)[left][right] != 1
                        || neighbor != full_vertices) {
                        throw std::runtime_error(
                            "tau-even fundamental edge is not unique"
                        );
                    }
                    neighbor = vertex;
                }
            }
            if (neighbor == full_vertices) {
                throw std::runtime_error(
                    "tau-even fundamental edge is missing"
                );
            }
            const Graph stripped = strip_vertex(*component, f_vertex);
            const std::size_t reduced_neighbor = compressed_vertex(
                neighbor, f_vertex
            );
            const std::int64_t orientation = neighbor < f_vertex ? -1 : 1;
            for (const auto& [graph_key, coefficient] :
                 straighten(stripped, reduced_memo)) {
                const Graph standard = unflatten(
                    graph_key, full_vertices - 1U
                );
                const Key key = module_key(reduced_neighbor, standard);
                std::int64_t& value = iota_row[key];
                value = normalized(
                    value + orientation * coefficient
                );
                if (value == 0) {
                    iota_row.erase(key);
                }
            }
        }
        for (const auto& [key, coefficient] : iota_row) {
            (void)coefficient;
            iota_columns.emplace(key, iota_columns.size());
        }
        iota_rows.push_back(std::move(iota_row));
    }

    std::vector<std::vector<std::int64_t>> global_matrix(
        global_rows.size(),
        std::vector<std::int64_t>(global_columns.size(), 0)
    );
    std::vector<std::vector<std::int64_t>> iota_matrix(
        iota_rows.size(), std::vector<std::int64_t>(iota_columns.size(), 0)
    );
    for (std::size_t row = 0U; row < global_rows.size(); ++row) {
        for (const auto& [key, coefficient] : global_rows[row]) {
            global_matrix[row][global_columns.at(key)] = coefficient;
        }
        for (const auto& [key, coefficient] : iota_rows[row]) {
            iota_matrix[row][iota_columns.at(key)] = coefficient;
        }
    }
    std::array<std::size_t, 2> sigma_domain{};
    std::array<std::size_t, 2> sigma_global_kernel{};
    std::array<std::size_t, 2> sigma_iota_kernel{};
    bool sigma_valid = labels[0] == labels[1];
    if (sigma_valid) {
        std::map<Key, std::size_t> representative_lookup;
        for (std::size_t source = 0U; source < representatives.size();
             ++source) {
            if (!representative_lookup.emplace(
                    flatten(representatives[source]), source
                ).second) {
                sigma_valid = false;
                break;
            }
        }
        std::vector<std::size_t> sigma_source(representatives.size(), 0U);
        if (sigma_valid) {
            for (std::size_t source = 0U; source < representatives.size();
                 ++source) {
                const Graph swapped = swap_graph_vertices(
                    representatives[source], 0U, 1U
                );
                const auto found = representative_lookup.find(flatten(swapped));
                if (found == representative_lookup.end()) {
                    sigma_valid = false;
                    break;
                }
                sigma_source[source] = found->second;
            }
        }
        if (sigma_valid) {
            for (std::size_t parity = 0U; parity < 2U; ++parity) {
                std::vector<std::vector<std::int64_t>> parity_global;
                std::vector<std::vector<std::int64_t>> parity_iota;
                for (std::size_t source = 0U;
                     source < representatives.size(); ++source) {
                    const std::size_t partner = sigma_source[source];
                    if (source > partner || (parity == 1U && source == partner)) {
                        continue;
                    }
                    std::vector<std::int64_t> global_row(
                        global_columns.size(), 0
                    );
                    std::vector<std::int64_t> iota_row(
                        iota_columns.size(), 0
                    );
                    const std::int64_t sign = parity == 0U ? 1 : -1;
                    for (std::size_t column = 0U;
                         column < global_columns.size(); ++column) {
                        global_row[column] = normalized(
                            global_matrix[source][column]
                            + sign * global_matrix[partner][column]
                        );
                    }
                    for (std::size_t column = 0U;
                         column < iota_columns.size(); ++column) {
                        iota_row[column] = normalized(
                            iota_matrix[source][column]
                            + sign * iota_matrix[partner][column]
                        );
                    }
                    parity_global.push_back(std::move(global_row));
                    parity_iota.push_back(std::move(iota_row));
                }
                sigma_domain[parity] = parity_global.size();
                sigma_global_kernel[parity] = sigma_domain[parity]
                    - modular_rank(std::move(parity_global));
                sigma_iota_kernel[parity] = sigma_domain[parity]
                    - modular_rank(std::move(parity_iota));
            }
        }
    }
    const std::size_t global_rank = modular_rank(std::move(global_matrix));
    const std::size_t iota_rank = modular_rank(std::move(iota_matrix));
    TauEvenResult result;
    result.domain_dimension = representatives.size();
    result.global_kernel = representatives.size() - global_rank;
    result.iota_kernel = representatives.size() - iota_rank;
    result.carrier_kernel = result.global_kernel - result.iota_kernel;
    std::vector<int> without_fundamental_pair;
    without_fundamental_pair.reserve(labels.size() - 2U);
    for (std::size_t vertex = 0U; vertex < labels.size(); ++vertex) {
        if (vertex != a_vertex && vertex != f_vertex) {
            without_fundamental_pair.push_back(labels[vertex]);
        }
    }
    result.fundamental_pair_reservoir = invariant_multiplicity(
        without_fundamental_pair
    );
    result.sigma_domain = sigma_domain;
    result.sigma_global_kernel = sigma_global_kernel;
    result.sigma_iota_kernel = sigma_iota_kernel;
    result.sigma_valid = sigma_valid;
    result.commuting_rank = global_rank <= iota_rank;
    if (include_separation_diagnostics) {
        const std::vector<Graph> merged_graphs = tau_even_merged_graphs(labels);
        const auto [merged_raw_collisions, merged_raw_order]
            = best_raw_leading_collision_count(merged_graphs);
        (void)merged_raw_order;
        const auto [merged_standard_collisions, merged_standard_order]
            = best_leading_collision_count(merged_graphs);
        result.merged_raw_collisions = merged_raw_collisions;
        result.merged_standard_collisions = merged_standard_collisions;
        result.merged_standard_order = merged_standard_order;
        result.divisor_separated = divisor_separates_graph_monomials(
            merged_graphs
        );
        result.separation_diagnostics_valid = true;
    }
    return result;
}

struct FullBackgroundAllocationResult {
    std::size_t domain_dimension = 0U;
    std::size_t global_kernel = 0U;
    std::size_t iota_kernel = 0U;
    std::size_t carrier_kernel = 0U;
    std::size_t weighted_stacked_kernel = 0U;
    std::size_t weighted_power = 0U;
    std::array<std::size_t, 8> weighted_trial_stacked_kernel{};
    std::array<std::size_t, 8> hashed_trial_stacked_kernel{};
    bool weighted_hash = false;
    bool weighted_payment_rank = false;
    std::array<std::size_t, 2> sigma_domain{};
    std::array<std::size_t, 2> sigma_global_kernel{};
    std::array<std::size_t, 2> sigma_iota_kernel{};
    std::array<std::size_t, 2> sigma_weighted_stacked_kernel{};
    std::array<bool, 2> sigma_weighted_payment_rank{};
    bool sigma_valid = false;
};

FullBackgroundAllocationResult analyze_full_background_allocation(
    const std::vector<int>& labels
) {
    constexpr std::size_t fundamental = 3U;
    const std::size_t full_vertices = labels.size();
    if (full_vertices < 7U
        || full_vertices >= std::numeric_limits<unsigned int>::digits
        || labels[fundamental] != 1) {
        throw std::runtime_error("invalid full-background allocation labels");
    }
    std::vector<int> reduced_degree;
    reduced_degree.reserve(full_vertices - 1U);
    for (std::size_t vertex = 0U; vertex < full_vertices; ++vertex) {
        if (vertex != fundamental) {
            reduced_degree.push_back(labels[vertex]);
        }
    }
    const std::size_t reduced_vertices = reduced_degree.size();
    std::map<Key, std::size_t> f0_indices;
    for (std::size_t vertex = 0U; vertex < reduced_vertices; ++vertex) {
        std::vector<int> base(reduced_vertices, 0);
        ++base[vertex];
        std::vector<int> residual;
        if (!subtract_degree(reduced_degree, base, residual)) {
            continue;
        }
        for (Graph graph : noncrossing_graphs(residual)) {
            f0_indices.emplace(
                module_key(vertex, graph), f0_indices.size()
            );
        }
    }

    std::vector<Graph> sources;
    std::vector<unsigned int> source_masks;
    for (std::size_t first = 0U; first < full_vertices; ++first) {
        for (std::size_t second = first + 1U; second < full_vertices;
             ++second) {
            for (std::size_t third = second + 1U; third < full_vertices;
                 ++third) {
                const unsigned int mask = (1U << first) | (1U << second)
                    | (1U << third);
                const int minus_parity = static_cast<int>(mask & 1U)
                    ^ static_cast<int>((mask >> 1U) & 1U);
                if (minus_parity == 0) {
                    continue;
                }
                const unsigned int full_mask
                    = (1U << full_vertices) - 1U;
                const auto left = noncrossing_graphs(labels, mask);
                const auto right = noncrossing_graphs(
                    labels, full_mask ^ mask
                );
                for (const Graph& left_graph : left) {
                    for (const Graph& right_graph : right) {
                        sources.push_back(graph_sum(left_graph, right_graph));
                        source_masks.push_back(mask);
                    }
                }
            }
        }
    }

    std::map<Key, Expansion> reduced_memo;
    std::vector<std::vector<std::int64_t>> iota_matrix(
        sources.size(),
        std::vector<std::int64_t>(f0_indices.size(), 0)
    );
    for (std::size_t source = 0U; source < sources.size(); ++source) {
        std::size_t old_neighbor = full_vertices;
        for (std::size_t vertex = 0U; vertex < full_vertices; ++vertex) {
            if (vertex == fundamental) {
                continue;
            }
            const std::size_t left = std::min(vertex, fundamental);
            const std::size_t right = std::max(vertex, fundamental);
            if (sources[source][left][right] != 0) {
                if (sources[source][left][right] != 1
                    || old_neighbor != full_vertices) {
                    throw std::runtime_error(
                        "full-background fundamental edge is not unique"
                    );
                }
                old_neighbor = vertex;
            }
        }
        if (old_neighbor == full_vertices) {
            throw std::runtime_error(
                "full-background fundamental edge is missing"
            );
        }
        const Graph stripped = strip_vertex(sources[source], fundamental);
        const std::size_t neighbor = compressed_vertex(
            old_neighbor, fundamental
        );
        const std::int64_t orientation
            = old_neighbor < fundamental ? -1 : 1;
        for (const auto& [graph_key, coefficient] :
             straighten(stripped, reduced_memo)) {
            const Graph standard = unflatten(graph_key, reduced_vertices);
            const auto found = f0_indices.find(module_key(neighbor, standard));
            if (found == f0_indices.end()) {
                throw std::runtime_error(
                    "full-background source left F0 multidegree"
                );
            }
            iota_matrix[source][found->second] = normalized(
                orientation * coefficient
            );
        }
    }

    std::map<Key, Expansion> full_memo;
    std::vector<Expansion> global_rows;
    std::map<Key, std::size_t> global_columns;
    global_rows.reserve(sources.size());
    for (const Graph& source : sources) {
        global_rows.push_back(straighten(source, full_memo));
        for (const auto& [key, coefficient] : global_rows.back()) {
            (void)coefficient;
            global_columns.emplace(key, global_columns.size());
        }
    }
    std::vector<std::vector<std::int64_t>> global_matrix(
        sources.size(),
        std::vector<std::int64_t>(global_columns.size(), 0)
    );
    for (std::size_t source = 0U; source < sources.size(); ++source) {
        for (const auto& [key, coefficient] : global_rows[source]) {
            global_matrix[source][global_columns.at(key)] = coefficient;
        }
    }

    constexpr std::size_t a_vertex = 2U;
    std::vector<int> retained_map(full_vertices, -1);
    std::size_t retained_count = 0U;
    for (std::size_t vertex = 0U; vertex < full_vertices; ++vertex) {
        if (vertex != a_vertex && vertex != fundamental) {
            retained_map[vertex] = static_cast<int>(retained_count);
            ++retained_count;
        }
    }
    std::map<Key, Expansion> residue_memo;
    std::vector<Expansion> residue_rows(sources.size());
    std::map<Key, std::size_t> residue_columns;
    for (std::size_t source = 0U; source < sources.size(); ++source) {
        std::size_t a_neighbor = full_vertices;
        std::size_t f_neighbor = full_vertices;
        for (std::size_t vertex = 0U; vertex < full_vertices; ++vertex) {
            if (vertex == a_vertex || vertex == fundamental) {
                continue;
            }
            const std::size_t a_left = std::min(a_vertex, vertex);
            const std::size_t a_right = std::max(a_vertex, vertex);
            const std::size_t f_left = std::min(fundamental, vertex);
            const std::size_t f_right = std::max(fundamental, vertex);
            if (sources[source][a_left][a_right] != 0) {
                a_neighbor = vertex;
            }
            if (sources[source][f_left][f_right] != 0) {
                f_neighbor = vertex;
            }
        }
        if (a_neighbor == full_vertices || f_neighbor == full_vertices
            || a_neighbor == f_neighbor) {
            continue;
        }
        Graph residue(
            retained_count, std::vector<int>(retained_count, 0)
        );
        for (std::size_t first = 0U; first < full_vertices; ++first) {
            if (first == a_vertex || first == fundamental) {
                continue;
            }
            for (std::size_t second = first + 1U; second < full_vertices;
                 ++second) {
                if (second == a_vertex || second == fundamental) {
                    continue;
                }
                residue[static_cast<std::size_t>(retained_map[first])]
                    [static_cast<std::size_t>(retained_map[second])]
                    = sources[source][first][second];
            }
        }
        const int a_orientation = a_vertex < a_neighbor ? 1 : -1;
        const int f_orientation = fundamental < f_neighbor ? 1 : -1;
        const int output_orientation = a_neighbor < f_neighbor ? 1 : -1;
        const std::int64_t scale = static_cast<std::int64_t>(
            a_orientation * f_orientation * output_orientation
        );
        const std::size_t first = static_cast<std::size_t>(
            retained_map[std::min(a_neighbor, f_neighbor)]
        );
        const std::size_t second = static_cast<std::size_t>(
            retained_map[std::max(a_neighbor, f_neighbor)]
        );
        ++residue[first][second];
        add_scaled(
            residue_rows[source], straighten(residue, residue_memo), scale
        );
        for (const auto& [key, coefficient] : residue_rows[source]) {
            (void)coefficient;
            residue_columns.emplace(key, residue_columns.size());
        }
    }
    FullBackgroundAllocationResult result;
    result.domain_dimension = sources.size();
    const std::size_t global_rank = modular_rank(global_matrix);
    const std::size_t iota_rank = modular_rank(iota_matrix);
    result.global_kernel = sources.size() - global_rank;
    result.iota_kernel = sources.size() - iota_rank;
    result.carrier_kernel = result.global_kernel - result.iota_kernel;
    std::size_t weighted_stacked_rank = 0U;
    std::vector<std::vector<std::int64_t>> weighted_stacked_matrix;
    for (std::size_t trial = 0U;
         trial < result.weighted_trial_stacked_kernel.size(); ++trial) {
        std::vector<std::vector<std::int64_t>> trial_matrix(
            sources.size(),
            std::vector<std::int64_t>(
                global_columns.size() + residue_columns.size(), 0
            )
        );
        for (std::size_t source = 0U; source < sources.size(); ++source) {
            for (std::size_t column = 0U; column < global_columns.size();
                 ++column) {
                trial_matrix[source][column] = global_matrix[source][column];
            }
            const std::int64_t weight = modular_power(
                static_cast<std::int64_t>(source_masks[source] + 1U),
                static_cast<std::int64_t>(trial + 1U)
            );
            for (const auto& [key, coefficient] : residue_rows[source]) {
                trial_matrix[source][
                    global_columns.size() + residue_columns.at(key)
                ] = (weight * coefficient) % modulus;
            }
        }
        const std::size_t trial_rank = modular_rank(trial_matrix);
        result.weighted_trial_stacked_kernel[trial]
            = sources.size() - trial_rank;
        if (trial == 0U || trial_rank > weighted_stacked_rank) {
            weighted_stacked_rank = trial_rank;
            weighted_stacked_matrix = std::move(trial_matrix);
            result.weighted_power = trial + 1U;
        }
    }
    for (std::size_t trial = 0U;
         trial < result.hashed_trial_stacked_kernel.size(); ++trial) {
        std::vector<std::vector<std::int64_t>> trial_matrix(
            sources.size(),
            std::vector<std::int64_t>(
                global_columns.size() + residue_columns.size(), 0
            )
        );
        for (std::size_t source = 0U; source < sources.size(); ++source) {
            for (std::size_t column = 0U; column < global_columns.size();
                 ++column) {
                trial_matrix[source][column] = global_matrix[source][column];
            }
            const std::int64_t weight = deterministic_mask_weight(
                source_masks[source], trial
            );
            for (const auto& [key, coefficient] : residue_rows[source]) {
                trial_matrix[source][
                    global_columns.size() + residue_columns.at(key)
                ] = (weight * coefficient) % modulus;
            }
        }
        const std::size_t trial_rank = modular_rank(trial_matrix);
        result.hashed_trial_stacked_kernel[trial]
            = sources.size() - trial_rank;
        if (trial_rank > weighted_stacked_rank) {
            weighted_stacked_rank = trial_rank;
            weighted_stacked_matrix = std::move(trial_matrix);
            result.weighted_power = trial + 1U;
            result.weighted_hash = true;
        }
    }
    result.weighted_stacked_kernel = sources.size() - weighted_stacked_rank;
    result.weighted_payment_rank = weighted_stacked_rank >= iota_rank;
    result.sigma_valid = labels[0U] == labels[1U];
    if (!result.sigma_valid) {
        return result;
    }
    std::map<std::pair<unsigned int, Key>, std::size_t> source_lookup;
    for (std::size_t source = 0U; source < sources.size(); ++source) {
        if (!source_lookup.emplace(
                std::make_pair(source_masks[source], flatten(sources[source])),
                source
            ).second) {
            result.sigma_valid = false;
            return result;
        }
    }
    std::vector<std::size_t> sigma_source(sources.size(), 0U);
    for (std::size_t source = 0U; source < sources.size(); ++source) {
        unsigned int swapped_mask = source_masks[source];
        const unsigned int first_bit = swapped_mask & 1U;
        const unsigned int second_bit = (swapped_mask >> 1U) & 1U;
        swapped_mask &= ~3U;
        swapped_mask |= first_bit << 1U;
        swapped_mask |= second_bit;
        const Graph swapped_graph = swap_graph_vertices(
            sources[source], 0U, 1U
        );
        const auto found = source_lookup.find(std::make_pair(
            swapped_mask, flatten(swapped_graph)
        ));
        if (found == source_lookup.end()) {
            result.sigma_valid = false;
            return result;
        }
        sigma_source[source] = found->second;
    }
    for (std::size_t parity = 0U; parity < 2U; ++parity) {
        std::vector<std::vector<std::int64_t>> parity_global;
        std::vector<std::vector<std::int64_t>> parity_iota;
        std::vector<std::vector<std::int64_t>> parity_weighted_stacked;
        for (std::size_t source = 0U; source < sources.size(); ++source) {
            const std::size_t partner = sigma_source[source];
            if (source > partner || (parity == 1U && source == partner)) {
                continue;
            }
            const std::int64_t sign = parity == 0U ? 1 : -1;
            const auto combine_rows = [source, partner, sign](
                const std::vector<std::vector<std::int64_t>>& matrix
            ) {
                const std::size_t columns = matrix.empty()
                    ? 0U : matrix.front().size();
                std::vector<std::int64_t> row(columns, 0);
                for (std::size_t column = 0U; column < columns; ++column) {
                    row[column] = normalized(
                        matrix[source][column]
                        + sign * matrix[partner][column]
                    );
                }
                return row;
            };
            parity_global.push_back(combine_rows(global_matrix));
            parity_iota.push_back(combine_rows(iota_matrix));
            parity_weighted_stacked.push_back(
                combine_rows(weighted_stacked_matrix)
            );
        }
        result.sigma_domain[parity] = parity_iota.size();
        const std::size_t parity_global_rank = modular_rank(
            std::move(parity_global)
        );
        const std::size_t parity_iota_rank = modular_rank(parity_iota);
        const std::size_t parity_weighted_stacked_rank = modular_rank(
            std::move(parity_weighted_stacked)
        );
        result.sigma_global_kernel[parity] = result.sigma_domain[parity]
            - parity_global_rank;
        result.sigma_iota_kernel[parity] = result.sigma_domain[parity]
            - parity_iota_rank;
        result.sigma_weighted_stacked_kernel[parity]
            = result.sigma_domain[parity] - parity_weighted_stacked_rank;
        result.sigma_weighted_payment_rank[parity]
            = parity_weighted_stacked_rank >= parity_iota_rank;
    }
    return result;
}

struct SevenLiftResult {
    std::size_t source_kernel = 0U;
    std::size_t iota_kernel = 0U;
    std::size_t carrier_kernel = 0U;
    std::int64_t fundamental_pair_reservoir = 0;
    std::size_t split_residue_rank = 0U;
    std::size_t split_allocation_rank = 0U;
    std::size_t linear_residue_rank = 0U;
    std::size_t linear_kernel_iota_rank = 0U;
    std::size_t minus_residue_rank = 0U;
    std::size_t side_global_defect = 0U;
    std::size_t q_stacked_kernel_dimension = 0U;
    std::size_t q_stacked_iota_rank = 0U;
    std::size_t tau_even_best_raw_collisions = 0U;
    bool commuting_square = false;
};

SevenLiftResult analyze_seven_factor_lift(
    const std::vector<int>& labels,
    bool report
) {
    constexpr std::size_t full_vertices = 7U;
    constexpr std::size_t fundamental = 3U;
    constexpr unsigned int full_mask = (1U << full_vertices) - 1U;
    if (labels.size() != full_vertices || labels[fundamental] != 1) {
        throw std::runtime_error("invalid seven-factor lift labels");
    }
    std::vector<int> reduced_degree;
    reduced_degree.reserve(full_vertices - 1U);
    for (std::size_t vertex = 0U; vertex < full_vertices; ++vertex) {
        if (vertex != fundamental) {
            reduced_degree.push_back(labels[vertex]);
        }
    }
    const std::size_t vertices = reduced_degree.size();

    std::vector<F0BasisElement> f0_basis;
    std::map<Key, std::size_t> f0_indices;
    for (std::size_t i = 0U; i < vertices; ++i) {
        std::vector<int> base(vertices, 0);
        ++base[i];
        std::vector<int> residual;
        if (!subtract_degree(reduced_degree, base, residual)) {
            continue;
        }
        for (Graph graph : noncrossing_graphs(residual)) {
            const std::size_t index = f0_basis.size();
            f0_indices.emplace(module_key(i, graph), index);
            f0_basis.push_back({i, std::move(graph)});
        }
    }

    std::vector<F1BasisElement> f1_basis;
    for (std::size_t i = 0U; i < vertices; ++i) {
        for (std::size_t j = i + 1U; j < vertices; ++j) {
            for (std::size_t k = j + 1U; k < vertices; ++k) {
                std::vector<int> base(vertices, 0);
                ++base[i];
                ++base[j];
                ++base[k];
                std::vector<int> residual;
                if (!subtract_degree(reduced_degree, base, residual)) {
                    continue;
                }
                for (Graph graph : noncrossing_graphs(residual)) {
                    const std::array<std::size_t, 3> triple{i, j, k};
                    f1_basis.push_back({triple, std::move(graph)});
                }
            }
        }
    }

    std::map<Key, Expansion> reduced_memo;
    std::vector<std::vector<std::int64_t>> differential(
        f1_basis.size(), std::vector<std::int64_t>(f0_basis.size(), 0)
    );
    const auto add_differential_term = [
        &reduced_memo, &f0_indices, &differential
    ](
        std::size_t source,
        const Graph& coefficient,
        std::size_t edge_first,
        std::size_t edge_second,
        std::size_t output_vertex,
        std::int64_t scalar
    ) {
        Graph product = coefficient;
        ++product[edge_first][edge_second];
        for (const auto& [graph_key, expansion_coefficient] :
             straighten(product, reduced_memo)) {
            const Graph standard = unflatten(graph_key, product.size());
            const auto found = f0_indices.find(module_key(output_vertex, standard));
            if (found == f0_indices.end()) {
                throw std::runtime_error("differential left F0 multidegree");
            }
            std::int64_t& value = differential[source][found->second];
            value = normalized(value + scalar * expansion_coefficient);
        }
    };
    for (std::size_t source = 0U; source < f1_basis.size(); ++source) {
        const auto [i, j, k] = f1_basis[source].triple;
        add_differential_term(
            source, f1_basis[source].coefficient, j, k, i, 1
        );
        add_differential_term(
            source, f1_basis[source].coefficient, i, k, j, -1
        );
        add_differential_term(
            source, f1_basis[source].coefficient, i, j, k, 1
        );
    }

    std::map<Key, Expansion> full_memo;
    std::vector<Expansion> phi_expansions(f0_basis.size());
    for (std::size_t source = 0U; source < f0_basis.size(); ++source) {
        const std::size_t old_neighbor = expanded_vertex(
            f0_basis[source].vertex, fundamental
        );
        const std::int64_t orientation = old_neighbor < fundamental ? -1 : 1;
        const Graph full_graph = insert_vertex_and_edge(
            f0_basis[source].coefficient,
            fundamental,
            f0_basis[source].vertex
        );
        add_scaled(
            phi_expansions[source], straighten(full_graph, full_memo), orientation
        );
    }

    std::vector<Graph> sources;
    std::vector<unsigned int> source_masks;
    for (unsigned int mask = 0U; mask <= full_mask; ++mask) {
        if (popcount(mask) != 3) {
            continue;
        }
        const int minus_parity = static_cast<int>(mask & 1U)
            ^ static_cast<int>((mask >> 1U) & 1U);
        if (minus_parity == 0) {
            continue;
        }
        const auto first = noncrossing_graphs(labels, mask);
        const auto second = noncrossing_graphs(labels, full_mask ^ mask);
        for (const Graph& left : first) {
            for (const Graph& right : second) {
                sources.push_back(graph_sum(left, right));
                source_masks.push_back(mask);
            }
        }
    }
    std::vector<Expansion> source_expansions;
    std::vector<SparseVector> iota_rows;
    std::map<Key, std::size_t> global_indices;
    bool commuting_square = true;
    for (const Graph& source : sources) {
        source_expansions.push_back(straighten(source, full_memo));
        for (const auto& [key, coefficient] : source_expansions.back()) {
            (void)coefficient;
            global_indices.emplace(key, global_indices.size());
        }
        std::size_t old_neighbor = full_vertices;
        for (std::size_t vertex = 0U; vertex < full_vertices; ++vertex) {
            if (vertex == fundamental) {
                continue;
            }
            const std::size_t left = std::min(vertex, fundamental);
            const std::size_t right = std::max(vertex, fundamental);
            if (source[left][right] != 0) {
                if (source[left][right] != 1 || old_neighbor != full_vertices) {
                    throw std::runtime_error("fundamental edge is not unique");
                }
                old_neighbor = vertex;
            }
        }
        if (old_neighbor == full_vertices) {
            throw std::runtime_error("fundamental edge is missing");
        }
        Graph stripped = strip_vertex(source, fundamental);
        const std::size_t neighbor = compressed_vertex(old_neighbor, fundamental);
        const std::int64_t orientation = old_neighbor < fundamental ? -1 : 1;
        SparseVector iota;
        for (const auto& [graph_key, coefficient] :
             straighten(stripped, reduced_memo)) {
            const Graph standard = unflatten(graph_key, vertices);
            const auto found = f0_indices.find(module_key(neighbor, standard));
            if (found == f0_indices.end()) {
                throw std::runtime_error("stripped source left F0 multidegree");
            }
            iota[found->second] = normalized(orientation * coefficient);
        }
        Expansion recomposed;
        for (const auto& [f0_source, coefficient] : iota) {
            add_scaled(recomposed, phi_expansions[f0_source], coefficient);
        }
        if (recomposed != source_expansions.back()) {
            commuting_square = false;
        }
        iota_rows.push_back(std::move(iota));
    }

    std::vector<std::vector<std::int64_t>> source_matrix(
        sources.size(), std::vector<std::int64_t>(global_indices.size(), 0)
    );
    std::vector<std::vector<std::int64_t>> iota_matrix(
        sources.size(), std::vector<std::int64_t>(f0_basis.size(), 0)
    );
    for (std::size_t source = 0U; source < sources.size(); ++source) {
        for (const auto& [key, coefficient] : source_expansions[source]) {
            source_matrix[source][global_indices.at(key)] = coefficient;
        }
        for (const auto& [column, coefficient] : iota_rows[source]) {
            iota_matrix[source][column] = coefficient;
        }
    }
    const auto source_kernel = modular_left_kernel(source_matrix);
    const auto iota_kernel_basis = modular_left_kernel(iota_matrix);
    const std::size_t iota_rank = modular_rank(iota_matrix);
    const std::size_t iota_kernel = sources.size() - iota_rank;
    std::vector<int> without_fundamental_pair;
    for (std::size_t vertex = 0U; vertex < labels.size(); ++vertex) {
        if (vertex != 2U && vertex != fundamental) {
            without_fundamental_pair.push_back(labels[vertex]);
        }
    }
    const std::int64_t fundamental_pair_reservoir = labels[2] == 1
        ? invariant_multiplicity(without_fundamental_pair) : 0;
    std::size_t split_residue_rank = 0U;
    std::size_t split_allocation_rank = 0U;
    std::array<std::size_t, 8> weighted_trial_ranks{};
    std::array<std::size_t, 6> component_weight_ranks{};
    std::array<std::size_t, 2> side_source_counts{};
    std::array<std::size_t, 2> side_global_ranks{};
    std::array<std::size_t, 2> side_augmented_ranks{};
    std::size_t q_stacked_kernel_dimension = 0U;
    std::size_t q_stacked_iota_rank = 0U;
    std::array<std::size_t, 2> tau_projection_ranks{};
    std::array<std::size_t, 2> tau_projection_iota_ranks{};
    std::size_t tau_even_dimension = 0U;
    std::size_t tau_even_global_kernel = 0U;
    std::size_t tau_even_iota_kernel = 0U;
    std::size_t tau_even_raw_collisions = 0U;
    std::size_t tau_even_best_raw_collisions = 0U;
    std::vector<std::size_t> tau_even_best_raw_order;
    bool tau_source_involution = true;
    std::size_t linear_kernel_iota_rank = 0U;
    if (labels[2] == 1) {
        constexpr std::size_t a_vertex = 2U;
        std::vector<int> retained_map(full_vertices, -1);
        std::size_t retained_count = 0U;
        for (std::size_t vertex = 0U; vertex < full_vertices; ++vertex) {
            if (vertex != a_vertex && vertex != fundamental) {
                retained_map[vertex] = static_cast<int>(retained_count);
                ++retained_count;
            }
        }
        std::map<Key, Expansion> residue_memo;
        std::vector<Expansion> residue_rows(sources.size());
        std::map<Key, std::size_t> residue_columns;
        for (std::size_t source_index = 0U; source_index < sources.size();
             ++source_index) {
            const Graph& source = sources[source_index];
            std::size_t a_neighbor = full_vertices;
            std::size_t f_neighbor = full_vertices;
            for (std::size_t vertex = 0U; vertex < full_vertices; ++vertex) {
                if (vertex == a_vertex || vertex == fundamental) {
                    continue;
                }
                const std::size_t a_left = std::min(a_vertex, vertex);
                const std::size_t a_right = std::max(a_vertex, vertex);
                const std::size_t f_left = std::min(fundamental, vertex);
                const std::size_t f_right = std::max(fundamental, vertex);
                if (source[a_left][a_right] != 0) {
                    a_neighbor = vertex;
                }
                if (source[f_left][f_right] != 0) {
                    f_neighbor = vertex;
                }
            }
            if (a_neighbor == full_vertices || f_neighbor == full_vertices
                || a_neighbor == f_neighbor) {
                continue;
            }
            Graph residue(retained_count, std::vector<int>(retained_count, 0));
            for (std::size_t i = 0U; i < full_vertices; ++i) {
                if (i == a_vertex || i == fundamental) {
                    continue;
                }
                for (std::size_t j = i + 1U; j < full_vertices; ++j) {
                    if (j == a_vertex || j == fundamental) {
                        continue;
                    }
                    residue[static_cast<std::size_t>(retained_map[i])]
                        [static_cast<std::size_t>(retained_map[j])] = source[i][j];
                }
            }
            const int a_orientation = a_vertex < a_neighbor ? 1 : -1;
            const int f_orientation = fundamental < f_neighbor ? 1 : -1;
            const int output_orientation = a_neighbor < f_neighbor ? 1 : -1;
            const std::int64_t contraction_scale = static_cast<std::int64_t>(
                a_orientation * f_orientation * output_orientation
            );
            const std::size_t first = static_cast<std::size_t>(
                retained_map[std::min(a_neighbor, f_neighbor)]
            );
            const std::size_t second = static_cast<std::size_t>(
                retained_map[std::max(a_neighbor, f_neighbor)]
            );
            ++residue[first][second];
            add_scaled(
                residue_rows[source_index],
                straighten(residue, residue_memo), contraction_scale
            );
            for (const auto& [key, coefficient] : residue_rows[source_index]) {
                (void)coefficient;
                residue_columns.emplace(key, residue_columns.size());
            }
        }
        std::vector<std::vector<std::int64_t>> residue_matrix(
            sources.size(), std::vector<std::int64_t>(residue_columns.size(), 0)
        );
        for (std::size_t source = 0U; source < sources.size(); ++source) {
            for (const auto& [key, coefficient] : residue_rows[source]) {
                residue_matrix[source][residue_columns.at(key)] = coefficient;
            }
        }
        const auto residue_image_rank = [&residue_matrix](
            const std::vector<std::vector<std::int64_t>>& domain,
            const std::vector<std::int64_t>& weights
        ) {
            const std::size_t columns = residue_matrix.empty()
                ? 0U : residue_matrix.front().size();
            std::vector<std::vector<std::int64_t>> image(
                domain.size(), std::vector<std::int64_t>(columns, 0)
            );
            for (std::size_t row = 0U; row < domain.size(); ++row) {
                for (std::size_t source = 0U; source < residue_matrix.size();
                     ++source) {
                    if (domain[row][source] == 0) {
                        continue;
                    }
                    for (std::size_t column = 0U; column < columns; ++column) {
                        const std::int64_t weighted = (
                            domain[row][source] * weights[source]
                        ) % modulus;
                        image[row][column] = normalized(
                            image[row][column]
                            + weighted * residue_matrix[source][column]
                        );
                    }
                }
            }
            return modular_rank(std::move(image));
        };
        for (int trial = 0; trial < 8; ++trial) {
            std::vector<std::int64_t> weights(sources.size(), 0);
            for (std::size_t source = 0U; source < sources.size(); ++source) {
                weights[source] = modular_power(
                    static_cast<std::int64_t>(source_masks[source] + 1U),
                    static_cast<std::int64_t>(trial + 1)
                );
            }
            const std::size_t trial_rank = residue_image_rank(
                source_kernel, weights
            );
            weighted_trial_ranks[static_cast<std::size_t>(trial)] = trial_rank;
            split_residue_rank = std::max(
                split_residue_rank,
                trial_rank
            );
            split_allocation_rank = std::max(
                split_allocation_rank,
                residue_image_rank(iota_kernel_basis, weights)
            );
        }
        for (std::size_t mode = 0U; mode < component_weight_ranks.size();
             ++mode) {
            std::vector<std::int64_t> weights(sources.size(), 0);
            for (std::size_t source = 0U; source < sources.size(); ++source) {
                const unsigned int mask = source_masks[source];
                const std::int64_t minus_weight = static_cast<std::int64_t>(
                    mask & 3U
                );
                const std::int64_t fundamental_weight = static_cast<std::int64_t>(
                    mask & 12U
                );
                const std::int64_t even_weight = static_cast<std::int64_t>(
                    mask & 112U
                );
                if (mode == 0U) {
                    weights[source] = minus_weight;
                } else if (mode == 1U) {
                    weights[source] = fundamental_weight;
                } else if (mode == 2U) {
                    weights[source] = even_weight;
                } else if (mode == 3U) {
                    weights[source] = minus_weight + fundamental_weight;
                } else if (mode == 4U) {
                    weights[source] = minus_weight + even_weight;
                } else {
                    weights[source] = fundamental_weight + even_weight;
                }
            }
            component_weight_ranks[mode] = residue_image_rank(
                source_kernel, weights
            );
        }
        for (std::size_t side = 0U; side < 2U; ++side) {
            std::vector<std::vector<std::int64_t>> side_matrix;
            std::vector<std::vector<std::int64_t>> side_global_matrix;
            for (std::size_t source = 0U; source < sources.size(); ++source) {
                if (((source_masks[source] >> side) & 1U) == 0U) {
                    continue;
                }
                std::vector<std::int64_t> row(
                    global_indices.size() + residue_columns.size(), 0
                );
                for (std::size_t column = 0U; column < global_indices.size();
                     ++column) {
                    row[column] = source_matrix[source][column];
                }
                side_global_matrix.push_back(source_matrix[source]);
                for (std::size_t column = 0U; column < residue_columns.size();
                     ++column) {
                    row[global_indices.size() + column] =
                        residue_matrix[source][column];
                }
                side_matrix.push_back(std::move(row));
            }
            side_source_counts[side] = side_matrix.size();
            side_global_ranks[side] = modular_rank(
                std::move(side_global_matrix)
            );
            side_augmented_ranks[side] = modular_rank(std::move(side_matrix));
        }
        std::vector<std::vector<std::int64_t>> q_stacked_matrix(
            sources.size(),
            std::vector<std::int64_t>(
                global_indices.size() + residue_columns.size(), 0
            )
        );
        for (std::size_t source = 0U; source < sources.size(); ++source) {
            for (std::size_t column = 0U; column < global_indices.size();
                 ++column) {
                q_stacked_matrix[source][column] = source_matrix[source][column];
            }
            if ((source_masks[source] & 1U) != 0U) {
                for (std::size_t column = 0U;
                     column < residue_columns.size(); ++column) {
                    q_stacked_matrix[source][global_indices.size() + column] =
                        residue_matrix[source][column];
                }
            }
        }
        const auto q_stacked_kernel = modular_left_kernel(q_stacked_matrix);
        q_stacked_kernel_dimension = q_stacked_kernel.size();
        std::vector<std::vector<std::int64_t>> q_stacked_iota(
            q_stacked_kernel.size(),
            std::vector<std::int64_t>(f0_basis.size(), 0)
        );
        for (std::size_t row = 0U; row < q_stacked_kernel.size(); ++row) {
            for (std::size_t source = 0U; source < sources.size(); ++source) {
                for (std::size_t column = 0U; column < f0_basis.size();
                     ++column) {
                    q_stacked_iota[row][column] = normalized(
                        q_stacked_iota[row][column]
                        + q_stacked_kernel[row][source]
                            * iota_matrix[source][column]
                    );
                }
            }
        }
        q_stacked_iota_rank = modular_rank(q_stacked_iota);
        if (report && std::getenv("IAB_VERBOSE") != nullptr) {
            for (std::size_t row = 0U; row < q_stacked_kernel.size(); ++row) {
                std::cout << "IAB_KERNEL_RELATION " << row << ':';
                for (std::size_t source = 0U; source < sources.size();
                     ++source) {
                    std::int64_t coefficient = q_stacked_kernel[row][source];
                    if (coefficient == 0) {
                        continue;
                    }
                    if (coefficient > modulus / 2) {
                        coefficient -= modulus;
                    }
                    std::cout << " coefficient=" << coefficient
                              << " mask=" << source_masks[source]
                              << " graph=";
                    print_graph(sources[source]);
                }
                std::cout << '\n';
            }
        }

        std::map<std::pair<unsigned int, Key>, std::size_t> source_lookup;
        for (std::size_t source = 0U; source < sources.size(); ++source) {
            source_lookup.emplace(
                std::make_pair(source_masks[source], flatten(sources[source])),
                source
            );
        }
        std::vector<std::size_t> tau_source(sources.size(), 0U);
        for (std::size_t source = 0U; source < sources.size(); ++source) {
            unsigned int swapped_mask = source_masks[source];
            const unsigned int a_bit = (swapped_mask >> 2U) & 1U;
            const unsigned int f_bit = (swapped_mask >> 3U) & 1U;
            swapped_mask &= ~(12U);
            swapped_mask |= a_bit << 3U;
            swapped_mask |= f_bit << 2U;
            const Graph swapped_graph = swap_graph_vertices(
                sources[source], 2U, 3U
            );
            const auto found = source_lookup.find(std::make_pair(
                swapped_mask, flatten(swapped_graph)
            ));
            if (found == source_lookup.end()) {
                tau_source_involution = false;
                break;
            }
            tau_source[source] = found->second;
        }
        if (tau_source_involution) {
            std::vector<std::vector<std::int64_t>> tau_even_global;
            std::vector<std::vector<std::int64_t>> tau_even_iota;
            std::set<std::vector<int>> tau_even_signatures;
            std::vector<std::pair<Graph, Graph>> tau_even_graph_pairs;
            std::vector<std::size_t> identity_order(full_vertices, 0U);
            for (std::size_t vertex = 0U; vertex < full_vertices; ++vertex) {
                identity_order[vertex] = vertex;
            }
            for (std::size_t source = 0U; source < sources.size(); ++source) {
                if (source >= tau_source[source]) {
                    continue;
                }
                std::vector<std::int64_t> global_row(global_indices.size(), 0);
                std::vector<std::int64_t> iota_row(f0_basis.size(), 0);
                for (std::size_t column = 0U; column < global_indices.size();
                     ++column) {
                    global_row[column] = normalized(
                        source_matrix[source][column]
                        + source_matrix[tau_source[source]][column]
                    );
                }
                for (std::size_t column = 0U; column < f0_basis.size();
                     ++column) {
                    iota_row[column] = normalized(
                        iota_matrix[source][column]
                        + iota_matrix[tau_source[source]][column]
                    );
                }
                tau_even_global.push_back(std::move(global_row));
                tau_even_iota.push_back(std::move(iota_row));
                const std::vector<int> first_signature = raw_leading_signature(
                    sources[source], identity_order
                );
                const std::vector<int> second_signature = raw_leading_signature(
                    sources[tau_source[source]], identity_order
                );
                tau_even_signatures.insert(
                    std::max(first_signature, second_signature)
                );
                tau_even_graph_pairs.emplace_back(
                    sources[source], sources[tau_source[source]]
                );
            }
            tau_even_dimension = tau_even_global.size();
            tau_even_global_kernel = tau_even_dimension
                - modular_rank(tau_even_global);
            tau_even_iota_kernel = tau_even_dimension
                - modular_rank(tau_even_iota);
            tau_even_raw_collisions = tau_even_dimension
                - tau_even_signatures.size();
            const auto [best_collisions, best_order] =
                best_paired_raw_leading_collision_count(
                    tau_even_graph_pairs
                );
            tau_even_best_raw_collisions = best_collisions;
            tau_even_best_raw_order = best_order;
            for (std::size_t parity = 0U; parity < 2U; ++parity) {
                std::vector<std::vector<std::int64_t>> projected(
                    q_stacked_kernel.size(),
                    std::vector<std::int64_t>(sources.size(), 0)
                );
                for (std::size_t row = 0U; row < q_stacked_kernel.size();
                     ++row) {
                    for (std::size_t source = 0U; source < sources.size();
                         ++source) {
                        const std::int64_t tau_value = q_stacked_kernel[row]
                            [tau_source[source]];
                        projected[row][source] = parity == 0U
                            ? normalized(
                                q_stacked_kernel[row][source] + tau_value
                              )
                            : normalized(
                                q_stacked_kernel[row][source] - tau_value
                              );
                    }
                }
                tau_projection_ranks[parity] = modular_rank(projected);
                std::vector<std::vector<std::int64_t>> projected_iota(
                    projected.size(),
                    std::vector<std::int64_t>(f0_basis.size(), 0)
                );
                for (std::size_t row = 0U; row < projected.size(); ++row) {
                    for (std::size_t source = 0U; source < sources.size();
                         ++source) {
                        for (std::size_t column = 0U;
                             column < f0_basis.size(); ++column) {
                            projected_iota[row][column] = normalized(
                                projected_iota[row][column]
                                + projected[row][source]
                                    * iota_matrix[source][column]
                            );
                        }
                    }
                }
                tau_projection_iota_ranks[parity] = modular_rank(
                    std::move(projected_iota)
                );
            }
        }
        std::vector<std::vector<std::int64_t>> stacked_matrix(
            sources.size(),
            std::vector<std::int64_t>(
                global_indices.size() + residue_columns.size(), 0
            )
        );
        for (std::size_t source = 0U; source < sources.size(); ++source) {
            for (std::size_t column = 0U; column < global_indices.size();
                 ++column) {
                stacked_matrix[source][column] = source_matrix[source][column];
            }
            const std::int64_t weight = static_cast<std::int64_t>(
                source_masks[source] + 1U
            );
            for (std::size_t column = 0U; column < residue_columns.size();
                 ++column) {
                stacked_matrix[source][global_indices.size() + column] = (
                    weight * residue_matrix[source][column]
                ) % modulus;
            }
        }
        const auto stacked_kernel = modular_left_kernel(stacked_matrix);
        std::vector<std::vector<std::int64_t>> stacked_kernel_iota(
            stacked_kernel.size(), std::vector<std::int64_t>(f0_basis.size(), 0)
        );
        for (std::size_t row = 0U; row < stacked_kernel.size(); ++row) {
            for (std::size_t source = 0U; source < sources.size(); ++source) {
                if (stacked_kernel[row][source] == 0) {
                    continue;
                }
                for (std::size_t column = 0U; column < f0_basis.size();
                     ++column) {
                    stacked_kernel_iota[row][column] = normalized(
                        stacked_kernel_iota[row][column]
                        + stacked_kernel[row][source] * iota_matrix[source][column]
                    );
                }
            }
        }
        linear_kernel_iota_rank = modular_rank(std::move(stacked_kernel_iota));
    }
    const SevenLiftResult result{
        source_kernel.size(), iota_kernel,
        source_kernel.size() - iota_kernel,
        fundamental_pair_reservoir, split_residue_rank,
        split_allocation_rank, weighted_trial_ranks[0],
        linear_kernel_iota_rank, component_weight_ranks[0],
        (side_source_counts[0] - side_global_ranks[0])
            + (side_source_counts[1] - side_global_ranks[1]),
        q_stacked_kernel_dimension,
        q_stacked_iota_rank,
        tau_even_best_raw_collisions,
        commuting_square
    };
    if (report) {
        std::cout << "SEVEN_FACTOR_PLUCKER_LIFT labels=[";
        for (std::size_t i = 0U; i < labels.size(); ++i) {
            if (i != 0U) {
                std::cout << ',';
            }
            std::cout << labels[i];
        }
        std::cout << "] sources=" << sources.size()
                  << " global_rank=" << modular_rank(source_matrix)
                  << " kernel=" << source_kernel.size()
                  << " iota_rank=" << iota_rank
                  << " iota_kernel=" << iota_kernel
                  << " carrier_kernel=" << source_kernel.size() - iota_kernel
                  << " af_pair=" << fundamental_pair_reservoir
                  << " weighted_residue_rank=" << split_residue_rank
                  << " weighted_allocation_rank=" << split_allocation_rank
                  << " linear_kernel_iota_rank=" << linear_kernel_iota_rank
                  << " f0=" << f0_basis.size()
                  << " f1=" << f1_basis.size()
                  << " d_rank=" << modular_rank(differential)
                  << " square=" << (commuting_square ? "PASS" : "FAIL")
                  << '\n';
        if (labels[2] == 1) {
            std::cout << "WEIGHTED_TRIAL_RANKS";
            for (std::size_t trial = 0U; trial < weighted_trial_ranks.size();
                 ++trial) {
                std::cout << ' ' << trial + 1U << ':'
                          << weighted_trial_ranks[trial];
            }
            std::cout << '\n';
            std::cout << "COMPONENT_WEIGHT_RANKS minus:"
                      << component_weight_ranks[0]
                      << " fundamental:" << component_weight_ranks[1]
                      << " even:" << component_weight_ranks[2]
                      << " minus+fundamental:" << component_weight_ranks[3]
                      << " minus+even:" << component_weight_ranks[4]
                      << " fundamental+even:" << component_weight_ranks[5]
                      << '\n';
            std::cout << "SIDE_AUGMENTED_RANKS q:"
                      << side_augmented_ranks[0] << '/'
                      << side_source_counts[0] << " r:"
                      << side_augmented_ranks[1] << '/'
                      << side_source_counts[1] << '\n';
            std::cout << "SIDE_GLOBAL_RANKS q:"
                      << side_global_ranks[0] << '/'
                      << side_source_counts[0] << " r:"
                      << side_global_ranks[1] << '/'
                      << side_source_counts[1] << '\n';
            std::cout << "IAB_TAU_DIAGNOSTIC kernel="
                      << q_stacked_kernel_dimension
                      << " iota_rank=" << q_stacked_iota_rank
                      << " tau="
                      << (tau_source_involution ? "PASS" : "FAIL")
                      << " plus_rank=" << tau_projection_ranks[0]
                      << " plus_iota_rank="
                      << tau_projection_iota_ranks[0]
                      << " minus_rank=" << tau_projection_ranks[1]
                      << " minus_iota_rank="
                      << tau_projection_iota_ranks[1] << '\n';
            std::cout << "TAU_EVEN_SOURCE dimension=" << tau_even_dimension
                      << " global_kernel=" << tau_even_global_kernel
                      << " iota_kernel=" << tau_even_iota_kernel
                      << " raw_collisions=" << tau_even_raw_collisions
                      << " best_raw_collisions="
                      << tau_even_best_raw_collisions << " best_order=[";
            for (std::size_t place = 0U;
                 place < tau_even_best_raw_order.size(); ++place) {
                if (place != 0U) {
                    std::cout << ',';
                }
                std::cout << tau_even_best_raw_order[place];
            }
            std::cout << "]\n";
            for (std::size_t side = 0U; side < 2U; ++side) {
                std::vector<Graph> side_sources;
                for (std::size_t source = 0U; source < sources.size();
                     ++source) {
                    if (((source_masks[source] >> side) & 1U) != 0U) {
                        side_sources.push_back(sources[source]);
                    }
                }
                const auto [collisions, order] =
                    best_leading_collision_count(side_sources);
                const std::vector<std::size_t> canonical_order = side == 0U
                    ? std::vector<std::size_t>{0U, 2U, 4U, 5U, 6U, 1U, 3U}
                    : std::vector<std::size_t>{0U, 3U, 1U, 2U, 4U, 5U, 6U};
                const std::size_t canonical_collisions =
                    leading_collision_count(side_sources, canonical_order);
                std::cout << "SIDE_LEADING side="
                          << (side == 0U ? 'q' : 'r')
                          << " collisions=" << collisions
                          << " canonical_collisions=" << canonical_collisions
                          << " order=[";
                for (std::size_t position = 0U; position < order.size();
                     ++position) {
                    if (position != 0U) {
                        std::cout << ',';
                    }
                    std::cout << order[position];
                }
                std::cout << "]\n";

                const std::vector<std::size_t> normalized_order = side == 0U
                    ? std::vector<std::size_t>{0U, 2U, 4U, 5U, 6U, 1U, 3U}
                    : std::vector<std::size_t>{1U, 2U, 4U, 5U, 6U, 0U, 3U};
                std::vector<int> normalized_labels(full_vertices, 0);
                for (std::size_t position = 0U; position < full_vertices;
                     ++position) {
                    normalized_labels[position] = labels[
                        normalized_order[position]
                    ];
                }
                std::vector<Graph> normalized_sources;
                for (std::size_t fundamental_position : {1U, 6U}) {
                    for (std::size_t even_position = 2U;
                         even_position <= 4U; ++even_position) {
                        const unsigned int mask = 1U
                            | (1U << fundamental_position)
                            | (1U << even_position);
                        const auto triple = noncrossing_graphs(
                            normalized_labels, mask
                        );
                        const auto complement = noncrossing_graphs(
                            normalized_labels, full_mask ^ mask
                        );
                        for (const Graph& left : triple) {
                            for (const Graph& right : complement) {
                                normalized_sources.push_back(
                                    graph_sum(left, right)
                                );
                            }
                        }
                    }
                }
                std::vector<std::size_t> identity(full_vertices, 0U);
                for (std::size_t position = 0U; position < full_vertices;
                     ++position) {
                    identity[position] = position;
                }
                std::cout << "SIDE_NORMALIZED_LEADING side="
                          << (side == 0U ? 'q' : 'r')
                          << " sources=" << normalized_sources.size()
                          << " collisions="
                          << leading_collision_count(
                                 normalized_sources, identity
                             )
                          << '\n';
                const auto [tree_collisions, tree] =
                    best_tree_collision_count(normalized_sources);
                std::cout << "SIDE_NORMALIZED_TREE side="
                          << (side == 0U ? 'q' : 'r')
                          << " collisions=" << tree_collisions
                          << " splits=";
                for (std::size_t index = 0U; index < tree.size(); ++index) {
                    if (index != 0U) {
                        std::cout << ',';
                    }
                    std::cout << tree[index];
                }
                std::cout << '\n';
                const auto [raw_collisions, raw_order] =
                    best_raw_leading_collision_count(normalized_sources);
                std::cout << "SIDE_NORMALIZED_RAW side="
                          << (side == 0U ? 'q' : 'r')
                          << " collisions=" << raw_collisions
                          << " identity_collisions="
                          << raw_leading_collision_count(
                                 normalized_sources, identity
                             )
                          << " order=[";
                for (std::size_t place = 0U; place < raw_order.size();
                     ++place) {
                    if (place != 0U) {
                        std::cout << ',';
                    }
                    std::cout << raw_order[place];
                }
                std::cout << "]\n";
                if (std::getenv("OSI_VERBOSE") != nullptr) {
                    std::map<Key, Expansion> normalized_memo;
                    for (std::size_t source = 0U;
                         source < normalized_sources.size(); ++source) {
                        const Expansion& expansion = straighten(
                            normalized_sources[source], normalized_memo
                        );
                        std::cout << "  normalized_source " << source
                                  << " terms=" << expansion.size()
                                  << " graph=";
                        print_graph(normalized_sources[source]);
                        std::cout << " least=";
                        print_graph(unflatten(
                            expansion.begin()->first, full_vertices
                        ));
                        std::cout << " greatest=";
                        print_graph(unflatten(
                            expansion.rbegin()->first, full_vertices
                        ));
                        std::cout << '\n';
                        if (std::getenv("OSI_FULL_VERBOSE") != nullptr) {
                            for (const auto& [key, coefficient] : expansion) {
                                std::cout << "    term coefficient="
                                          << coefficient << " graph=";
                                print_graph(unflatten(key, full_vertices));
                                std::cout << '\n';
                            }
                        }
                    }
                }
                if (std::getenv("OSI_BASIS_SEARCH") != nullptr) {
                    std::vector<std::size_t> basis_order(full_vertices, 0U);
                    basis_order[0] = side;
                    std::size_t position = 1U;
                    for (std::size_t vertex = 0U; vertex < full_vertices;
                         ++vertex) {
                        if (vertex != side) {
                            basis_order[position] = vertex;
                            ++position;
                        }
                    }
                    std::size_t best_basis_collisions = sources.size();
                    std::vector<std::size_t> best_basis_order = basis_order;
                    do {
                        std::vector<int> basis_labels(full_vertices, 0);
                        std::size_t a_position = full_vertices;
                        std::size_t f_position = full_vertices;
                        std::vector<std::size_t> even_positions;
                        for (std::size_t place = 0U; place < full_vertices;
                             ++place) {
                            const std::size_t vertex = basis_order[place];
                            basis_labels[place] = labels[vertex];
                            if (vertex == 2U) {
                                a_position = place;
                            } else if (vertex == 3U) {
                                f_position = place;
                            } else if (vertex >= 4U) {
                                even_positions.push_back(place);
                            }
                        }
                        std::vector<Graph> basis_sources;
                        for (std::size_t fundamental_position :
                             {a_position, f_position}) {
                            for (std::size_t even_position : even_positions) {
                                const unsigned int mask = 1U
                                    | (1U << fundamental_position)
                                    | (1U << even_position);
                                const auto triple = noncrossing_graphs(
                                    basis_labels, mask
                                );
                                const auto complement = noncrossing_graphs(
                                    basis_labels, full_mask ^ mask
                                );
                                for (const Graph& left : triple) {
                                    for (const Graph& right : complement) {
                                        basis_sources.push_back(
                                            graph_sum(left, right)
                                        );
                                    }
                                }
                            }
                        }
                        const std::size_t basis_collisions =
                            leading_collision_count(basis_sources, identity);
                        if (basis_collisions < best_basis_collisions) {
                            best_basis_collisions = basis_collisions;
                            best_basis_order = basis_order;
                            if (basis_collisions == 0U) {
                                break;
                            }
                        }
                    } while (std::next_permutation(
                        basis_order.begin() + 1, basis_order.end()
                    ));
                    std::cout << "SIDE_BASIS_SEARCH side="
                              << (side == 0U ? 'q' : 'r')
                              << " collisions=" << best_basis_collisions
                              << " order=[";
                    for (std::size_t place = 0U;
                         place < best_basis_order.size(); ++place) {
                        if (place != 0U) {
                            std::cout << ',';
                        }
                        std::cout << best_basis_order[place];
                    }
                    std::cout << "]\n";
                }
            }
        }
    }
    if (!commuting_square || !report) {
        return result;
    }
    std::vector<std::vector<std::int64_t>> kernel_targets;
    kernel_targets.reserve(source_kernel.size());
    for (std::size_t relation = 0U; relation < source_kernel.size();
         ++relation) {
        std::vector<std::int64_t> target(f0_basis.size(), 0);
        for (std::size_t source = 0U; source < sources.size(); ++source) {
            const std::int64_t kernel_coefficient = source_kernel[relation][source];
            if (kernel_coefficient == 0) {
                continue;
            }
            for (const auto& [column, coefficient] : iota_rows[source]) {
                target[column] = normalized(
                    target[column] + kernel_coefficient * coefficient
                );
            }
        }
        kernel_targets.push_back(std::move(target));
    }
    std::vector<std::size_t> carrier_relations;
    std::vector<std::vector<std::int64_t>> carrier_targets;
    std::size_t carrier_rank = 0U;
    for (std::size_t relation = 0U; relation < kernel_targets.size();
         ++relation) {
        auto candidate = carrier_targets;
        candidate.push_back(kernel_targets[relation]);
        const std::size_t candidate_rank = modular_rank(std::move(candidate));
        if (candidate_rank > carrier_rank) {
            carrier_rank = candidate_rank;
            carrier_relations.push_back(relation);
            carrier_targets.push_back(kernel_targets[relation]);
        }
    }
    if (carrier_rank != source_kernel.size() - iota_kernel) {
        throw std::runtime_error("carrier quotient dimension mismatch");
    }
    std::vector<std::size_t> a_columns;
    for (std::size_t column = 0U; column < f0_basis.size(); ++column) {
        if (f0_basis[column].vertex == 2U) {
            a_columns.push_back(column);
        }
    }
    std::vector<std::vector<std::int64_t>> a_projection(
        carrier_targets.size(),
        std::vector<std::int64_t>(a_columns.size(), 0)
    );
    for (std::size_t row = 0U; row < carrier_targets.size(); ++row) {
        for (std::size_t column = 0U; column < a_columns.size(); ++column) {
            a_projection[row][column] = carrier_targets[row][a_columns[column]];
        }
    }
    const std::size_t a_projection_rank = modular_rank(std::move(a_projection));
    std::cout << "CARRIER_PROJECTION carrier_rank=" << carrier_rank
              << " a_projection_rank=" << a_projection_rank
              << " a_columns=" << a_columns.size() << '\n';
    std::vector<std::vector<std::int64_t>> two_terminal_differential;
    std::vector<std::vector<std::int64_t>> a_terminal_differential;
    for (std::size_t source = 0U; source < f1_basis.size(); ++source) {
        int terminal_count = 0;
        for (std::size_t vertex : f1_basis[source].triple) {
            if (vertex < 3U) {
                ++terminal_count;
            }
        }
        if (terminal_count >= 2) {
            two_terminal_differential.push_back(differential[source]);
        }
        if (std::find(
                f1_basis[source].triple.begin(),
                f1_basis[source].triple.end(), 2U
            ) != f1_basis[source].triple.end()) {
            a_terminal_differential.push_back(differential[source]);
        }
    }
    for (std::size_t carrier = 0U; carrier < carrier_relations.size();
         ++carrier) {
        const std::size_t relation = carrier_relations[carrier];
        const std::vector<std::int64_t>& target = kernel_targets[relation];
        const auto [solved, lift] = solve_left_combination(differential, target);
        if (!solved) {
            std::cout << "carrier_lift " << carrier << " source_relation="
                      << relation << " FAIL\n";
            continue;
        }
        const auto [two_terminal_solved, two_terminal_lift] =
            solve_left_combination(two_terminal_differential, target);
        const auto [a_terminal_solved, a_terminal_lift] =
            solve_left_combination(a_terminal_differential, target);
        std::map<unsigned int, std::size_t> mask_counts;
        std::size_t terms = 0U;
        for (std::size_t source = 0U; source < lift.size(); ++source) {
            if (lift[source] == 0) {
                continue;
            }
            ++terms;
            unsigned int terminal_mask = 0U;
            for (std::size_t vertex : f1_basis[source].triple) {
                if (vertex < 3U) {
                    terminal_mask |= 1U << vertex;
                }
            }
            ++mask_counts[terminal_mask];
        }
        std::cout << "carrier_lift " << carrier << " source_relation="
                  << relation << " PASS terms=" << terms
                  << " two_terminal="
                  << (two_terminal_solved ? "PASS" : "FAIL");
        if (two_terminal_solved) {
            const std::size_t two_terminal_terms = static_cast<std::size_t>(
                std::count_if(
                    two_terminal_lift.begin(), two_terminal_lift.end(),
                    [](std::int64_t value) { return value != 0; }
                )
            );
            std::cout << ':' << two_terminal_terms;
        }
        std::cout << " a_terminal="
                  << (a_terminal_solved ? "PASS" : "FAIL");
        if (a_terminal_solved) {
            const std::size_t a_terminal_terms = static_cast<std::size_t>(
                std::count_if(
                    a_terminal_lift.begin(), a_terminal_lift.end(),
                    [](std::int64_t value) { return value != 0; }
                )
            );
            std::cout << ':' << a_terminal_terms;
        }
        std::cout
                  << " terminal_masks=";
        bool first = true;
        for (const auto& [mask, count] : mask_counts) {
            if (!first) {
                std::cout << ',';
            }
            first = false;
            std::cout << mask << ':' << count;
        }
        std::cout << '\n';
        if (std::getenv("PLUCKER_VERBOSE") != nullptr) {
            for (std::size_t source = 0U; source < lift.size(); ++source) {
                if (lift[source] == 0) {
                    continue;
                }
                std::int64_t coefficient = lift[source];
                if (coefficient > modulus / 2) {
                    coefficient -= modulus;
                }
                const auto [i, j, k] = f1_basis[source].triple;
                std::cout << "  term coefficient=" << coefficient
                          << " triple=" << i << ',' << j << ',' << k
                          << " graph=";
                print_graph(f1_basis[source].coefficient);
                std::cout << '\n';
            }
        }
    }
    return result;
}

void enumerate_degrees_rec(
    std::size_t position,
    int remaining_total,
    int maximum_coordinate,
    std::vector<int>& degree,
    std::vector<std::vector<int>>& output
) {
    if (position == degree.size()) {
        int total = 0;
        for (int value : degree) {
            total += value;
        }
        if (total >= 3 && (total & 1) != 0) {
            output.push_back(degree);
        }
        return;
    }
    const int maximum = std::min(maximum_coordinate, remaining_total);
    for (int value = 0; value <= maximum; ++value) {
        degree[position] = value;
        enumerate_degrees_rec(
            position + 1U, remaining_total - value,
            maximum_coordinate, degree, output
        );
    }
}

void enumerate_even_backgrounds_rec(
    std::size_t position,
    int minimum_label,
    int maximum_label,
    std::vector<int>& background,
    std::vector<std::vector<int>>& output
) {
    if (position == background.size()) {
        output.push_back(background);
        return;
    }
    for (int label = minimum_label; label <= maximum_label; label += 2) {
        background[position] = label;
        enumerate_even_backgrounds_rec(
            position + 1U, label, maximum_label, background, output
        );
    }
}

void print_degree(const std::vector<int>& degree) {
    std::cout << '[';
    for (std::size_t i = 0U; i < degree.size(); ++i) {
        if (i != 0U) {
            std::cout << ',';
        }
        std::cout << degree[i];
    }
    std::cout << ']';
}

void print_graph(const Graph& graph) {
    std::cout << '{';
    bool first = true;
    for (std::size_t i = 0U; i < graph.size(); ++i) {
        for (std::size_t j = i + 1U; j < graph.size(); ++j) {
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

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 7 && std::string(argv[1]) == "--tau-even-v2-case") {
            std::vector<int> labels{0, 0, 1, 1, 0, 0, 0};
            const std::array<std::size_t, 5U> positions{
                0U, 1U, 4U, 5U, 6U
            };
            for (std::size_t argument = 2U; argument < 7U; ++argument) {
                const long parsed = std::strtol(
                    argv[static_cast<int>(argument)], nullptr, 10
                );
                if (parsed < 1 || parsed > std::numeric_limits<int>::max()) {
                    throw std::runtime_error("invalid tau-even V2-case label");
                }
                labels[positions[argument - 2U]] = static_cast<int>(parsed);
            }
            const std::vector<Graph> graphs = tau_even_merged_graphs(labels);
            std::cout << "TAU_EVEN_V2_CASE labels=";
            print_degree(labels);
            std::cout << " sources=" << graphs.size() << " choices=";
            bool any = false;
            for (std::size_t terminal = 0U; terminal < 2U; ++terminal) {
                for (std::size_t root = 2U; root < 6U; ++root) {
                    const bool direct = binary_terminal_v2_jet_separates(
                        graphs, terminal, 1U - terminal, root
                    );
                    const bool merged
                        = binary_terminal_v2_merge_jet_separates(
                            graphs, terminal, 1U - terminal, root
                        );
                    if (terminal != 0U || root != 2U) {
                        std::cout << ',';
                    }
                    std::cout << (terminal == 0U ? 'q' : 'r')
                              << (root - 2U) << ':' << direct << '/'
                              << merged;
                    any = any || direct || merged;
                }
            }
            std::cout << " any=" << any << '\n';
            return any ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        if (argc == 6
            && std::string(argv[1]) == "--tau-even-equal-factor-case") {
            std::vector<int> labels{0, 0, 1, 1, 0, 0, 0};
            for (std::size_t index = 0U; index < 4U; ++index) {
                const long parsed = std::strtol(argv[index + 2U], nullptr, 10);
                if (parsed < 1 || parsed > std::numeric_limits<int>::max()) {
                    throw std::runtime_error(
                        "invalid equal-terminal factor-case label"
                    );
                }
                if (index == 0U) {
                    labels[0U] = static_cast<int>(parsed);
                    labels[1U] = static_cast<int>(parsed);
                } else {
                    labels[index + 3U] = static_cast<int>(parsed);
                }
            }
            const std::vector<Graph> graphs = tau_even_merged_graphs(labels);
            std::map<Key, std::size_t> lookup;
            for (std::size_t source = 0U; source < graphs.size(); ++source) {
                lookup.emplace(flatten(graphs[source]), source);
            }
            for (std::size_t source = 0U; source < graphs.size(); ++source) {
                const Graph swapped
                    = swap_graph_vertices(graphs[source], 0U, 1U);
                const std::size_t partner = lookup.at(flatten(swapped));
                if (source > partner) {
                    continue;
                }
                std::cout << "TAU_EVEN_EQUAL_FACTOR_ROW source=" << source
                          << " partner=" << partner << " first=[";
                for (std::size_t root = 0U; root < 4U; ++root) {
                    if (root != 0U) {
                        std::cout << ',';
                    }
                    std::cout << graphs[source][0U][root + 2U];
                }
                std::cout << "] second=[";
                for (std::size_t root = 0U; root < 4U; ++root) {
                    if (root != 0U) {
                        std::cout << ',';
                    }
                    std::cout << graphs[source][1U][root + 2U];
                }
                std::cout << "]\n";
            }
            const EqualTerminalFactorResult result
                = analyze_equal_terminal_factors(labels);
            std::cout << "TAU_EVEN_EQUAL_FACTOR_CASE labels=";
            print_degree(labels);
            std::cout << " plus=" << result.plus_tensor_rank << '/'
                      << result.plus_domain << " wronskian="
                      << result.minus_wronskian_rank << '/'
                      << result.minus_domain << " strict_min_v2="
                      << result.strict_min_v2 << " components="
                      << result.component_forms_rank << '/'
                      << result.component_forms_domain
                      << " component_forms_v2="
                      << result.component_forms_v2
                      << " unordered_pairs_distinct="
                      << result.unordered_pairs_distinct
                      << " plus_leading_v2=" << result.plus_leading_v2
                      << " minus_leading_v2=" << result.minus_leading_v2
                      << '\n';
            return EXIT_SUCCESS;
        }
        if (argc == 3
            && std::string(argv[1]) == "--tau-even-equal-factor-sweep") {
            const long parsed_maximum = std::strtol(argv[2], nullptr, 10);
            if (parsed_maximum < 3
                || parsed_maximum > std::numeric_limits<int>::max()) {
                throw std::runtime_error(
                    "invalid equal-terminal factor-sweep maximum"
                );
            }
            const int maximum = static_cast<int>(parsed_maximum);
            std::vector<std::vector<int>> cases;
            for (int q = 3; q <= maximum; q += 2) {
                for (int first = 2; first <= maximum; first += 2) {
                    for (int second = first; second <= maximum; second += 2) {
                        for (int third = second; third <= maximum; third += 2) {
                            cases.push_back({
                                q, q, 1, 1, first, second, third
                            });
                        }
                    }
                }
            }
            std::atomic<std::size_t> nonempty{0U};
            std::atomic<std::size_t> plus_residual{0U};
            std::atomic<std::size_t> strict_min_v2_full{0U};
            std::atomic<std::size_t> q3_wronskian_residual{0U};
            std::atomic<std::size_t> special_wronskian_residual{0U};
            std::vector<std::vector<int>> strict_min_v2_residuals;
            std::vector<std::vector<int>> unexpected_plus;
            std::vector<std::vector<int>> unexpected_wronskian;
#pragma omp parallel for schedule(dynamic)
            for (std::size_t index = 0U; index < cases.size(); ++index) {
                const EqualTerminalFactorResult result
                    = analyze_equal_terminal_factors(cases[index]);
                if (result.plus_domain == 0U) {
                    continue;
                }
                nonempty.fetch_add(1U, std::memory_order_relaxed);
                if (result.strict_min_v2) {
                    strict_min_v2_full.fetch_add(1U, std::memory_order_relaxed);
                } else {
#pragma omp critical
                    strict_min_v2_residuals.push_back(cases[index]);
                }
                const bool special = cases[index][4U] == 2
                    && cases[index][5U] == cases[index][0U] + 1
                    && cases[index][6U] == cases[index][0U] + 1;
                if (result.plus_tensor_rank < result.plus_domain) {
                    if (special) {
                        plus_residual.fetch_add(1U, std::memory_order_relaxed);
                    } else {
#pragma omp critical
                        unexpected_plus.push_back(cases[index]);
                    }
                }
                if (result.minus_wronskian_rank < result.minus_domain) {
                    if (cases[index][0U] == 3) {
                        q3_wronskian_residual.fetch_add(
                            1U, std::memory_order_relaxed
                        );
                    } else if (special) {
                        special_wronskian_residual.fetch_add(
                            1U, std::memory_order_relaxed
                        );
                    } else {
#pragma omp critical
                        unexpected_wronskian.push_back(cases[index]);
                    }
                }
            }
            std::sort(
                strict_min_v2_residuals.begin(),
                strict_min_v2_residuals.end()
            );
            std::sort(unexpected_plus.begin(), unexpected_plus.end());
            std::sort(
                unexpected_wronskian.begin(), unexpected_wronskian.end()
            );
            for (const std::vector<int>& labels : strict_min_v2_residuals) {
                std::cout << "TAU_EVEN_EQUAL_STRICT_MIN_V2_RESIDUAL labels=";
                print_degree(labels);
                std::cout << '\n';
            }
            for (const std::vector<int>& labels : unexpected_plus) {
                std::cout << "TAU_EVEN_EQUAL_FACTOR_UNEXPECTED_PLUS labels=";
                print_degree(labels);
                std::cout << '\n';
            }
            for (const std::vector<int>& labels : unexpected_wronskian) {
                std::cout
                    << "TAU_EVEN_EQUAL_FACTOR_UNEXPECTED_WRONSKIAN labels=";
                print_degree(labels);
                std::cout << '\n';
            }
            std::cout
                << "TAU_EVEN_EQUAL_FACTOR_SWEEP_DONE maximum=" << maximum
                << " candidate_cases=" << cases.size()
                << " nonempty=" << nonempty.load()
                << " plus_residual=" << plus_residual.load()
                << " strict_min_v2_full=" << strict_min_v2_full.load()
                << " q3_wronskian_residual="
                << q3_wronskian_residual.load()
                << " special_wronskian_residual="
                << special_wronskian_residual.load()
                << " unexpected_plus=" << unexpected_plus.size()
                << " unexpected_wronskian="
                << unexpected_wronskian.size()
                << " threads=" << omp_get_max_threads() << '\n';
            return unexpected_plus.empty() && unexpected_wronskian.empty()
                ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        if (argc == 3
            && std::string(argv[1]) == "--tau-even-source-formula-sweep") {
            const long parsed_maximum = std::strtol(argv[2], nullptr, 10);
            if (parsed_maximum < 5
                || parsed_maximum > std::numeric_limits<int>::max()) {
                throw std::runtime_error(
                    "invalid terminal-source formula-sweep maximum"
                );
            }
            const int maximum = static_cast<int>(parsed_maximum);
            std::vector<std::vector<int>> cases;
            for (int q = 3; q <= maximum; q += 2) {
                for (int r = q; r <= maximum; r += 2) {
                    for (int first = 2; first <= maximum; first += 2) {
                        for (int second = first; second <= maximum;
                             second += 2) {
                            for (int third = second; third <= maximum;
                                 third += 2) {
                                cases.push_back({
                                    q, r, 1, 1, first, second, third
                                });
                            }
                        }
                    }
                }
            }
            std::atomic<std::size_t> nonempty{0U};
            std::vector<std::vector<int>> failures;
#pragma omp parallel for schedule(dynamic)
            for (std::size_t index = 0U; index < cases.size(); ++index) {
                const std::multiset<TerminalExponentPair> explicit_rows
                    = explicit_terminal_source_rows(cases[index]);
                const std::multiset<TerminalExponentPair> generated_rows
                    = generated_terminal_source_rows(cases[index]);
                if (!generated_rows.empty()) {
                    nonempty.fetch_add(1U, std::memory_order_relaxed);
                }
                if (explicit_rows != generated_rows) {
#pragma omp critical
                    failures.push_back(cases[index]);
                }
            }
            for (const std::vector<int>& labels : failures) {
                std::cout << "TAU_EVEN_SOURCE_FORMULA_FAILURE labels=";
                print_degree(labels);
                std::cout << " explicit_rows="
                          << explicit_terminal_source_rows(labels).size()
                          << " generated_rows="
                          << generated_terminal_source_rows(labels).size()
                          << '\n';
            }
            std::cout
                << "TAU_EVEN_SOURCE_FORMULA_SWEEP_DONE maximum=" << maximum
                << " cases=" << cases.size()
                << " nonempty=" << nonempty.load()
                << " failures=" << failures.size()
                << " threads=" << omp_get_max_threads() << '\n';
            return failures.empty() ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        if (argc == 6
            && std::string(argv[1]) == "--tau-even-equal-diagonal") {
            std::vector<int> labels{0, 0, 1, 1, 0, 0, 0};
            const long parsed_q = std::strtol(argv[2], nullptr, 10);
            if (parsed_q < 3 || (parsed_q & 1L) == 0
                || parsed_q > std::numeric_limits<int>::max()) {
                throw std::runtime_error("invalid equal-terminal odd label");
            }
            labels[0U] = static_cast<int>(parsed_q);
            labels[1U] = static_cast<int>(parsed_q);
            for (std::size_t argument = 3U; argument < 6U; ++argument) {
                const long parsed = std::strtol(
                    argv[static_cast<int>(argument)], nullptr, 10
                );
                if (parsed < 2 || (parsed & 1L) != 0
                    || parsed > std::numeric_limits<int>::max()) {
                    throw std::runtime_error("invalid equal-terminal even label");
                }
                labels[argument + 1U] = static_cast<int>(parsed);
            }
            if (labels[4U] > labels[5U] || labels[5U] > labels[6U]) {
                throw std::runtime_error("equal-terminal labels are not ordered");
            }
            const EqualTerminalDiagonalResult result
                = analyze_equal_terminal_diagonal(labels);
            std::cout << "TAU_EVEN_EQUAL_DIAGONAL labels=";
            print_degree(labels);
            std::cout << " plus=" << result.plus_rank << '/'
                      << result.plus_domain
                      << " plus_terminal_tensor="
                      << result.plus_terminal_tensor_rank << '/'
                      << result.plus_domain
                      << " plus_diagonal_v2=" << result.plus_diagonal_v2
                      << " minus=" << result.minus_rank << '/'
                      << result.minus_domain
                      << " minus_terminal_tensor="
                      << result.minus_terminal_tensor_rank << '/'
                      << result.minus_domain
                      << " minus_wronskian="
                      << result.minus_wronskian_rank << '/'
                      << result.minus_domain
                      << " antisymmetric_divisible="
                      << result.antisymmetric_divisible << '\n';
            return result.antisymmetric_divisible
                    && result.plus_rank == result.plus_domain
                ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        if (argc == 3
            && std::string(argv[1]) == "--tau-even-equal-diagonal-sweep") {
            const long parsed_maximum = std::strtol(argv[2], nullptr, 10);
            if (parsed_maximum < 3
                || parsed_maximum > std::numeric_limits<int>::max()) {
                throw std::runtime_error(
                    "invalid equal-terminal diagonal-sweep maximum"
                );
            }
            const int maximum = static_cast<int>(parsed_maximum);
            std::vector<std::vector<int>> cases;
            for (int q = 3; q <= maximum; q += 2) {
                for (int first = 2; first <= maximum; first += 2) {
                    for (int second = first; second <= maximum; second += 2) {
                        for (int third = second; third <= maximum; third += 2) {
                            cases.push_back({
                                q, q, 1, 1, first, second, third
                            });
                        }
                    }
                }
            }
            struct EqualFailure {
                std::vector<int> labels;
                EqualTerminalDiagonalResult result;
            };
            std::atomic<std::size_t> nonempty{0U};
            std::atomic<std::size_t> plus_full{0U};
            std::atomic<std::size_t> minus_expected{0U};
            std::vector<EqualFailure> failures;
            std::vector<EqualFailure> wronskian_residuals;
#pragma omp parallel for schedule(dynamic)
            for (std::size_t index = 0U; index < cases.size(); ++index) {
                const EqualTerminalDiagonalResult result
                    = analyze_equal_terminal_diagonal(cases[index]);
                if (result.plus_domain == 0U) {
                    continue;
                }
                nonempty.fetch_add(1U, std::memory_order_relaxed);
                const bool symmetric_residual_family
                    = cases[index][4U] == 2
                    && cases[index][5U] == cases[index][0U] + 1
                    && cases[index][6U] == cases[index][0U] + 1;
                const bool plus_pass = result.antisymmetric_divisible
                    && (result.plus_terminal_tensor_rank == result.plus_domain
                        || symmetric_residual_family);
                const bool cubic = cases[index]
                    == std::vector<int>{3, 3, 1, 1, 2, 2, 2};
                const std::size_t expected_minus_kernel = cubic ? 1U : 0U;
                const bool minus_pass = cubic
                    ? (result.minus_domain >= result.minus_rank
                        && result.minus_domain - result.minus_rank
                            == expected_minus_kernel)
                    : (result.minus_terminal_tensor_rank == result.minus_domain
                        || result.minus_rank == result.minus_domain);
                if (plus_pass) {
                    plus_full.fetch_add(1U, std::memory_order_relaxed);
                }
                if (minus_pass) {
                    minus_expected.fetch_add(1U, std::memory_order_relaxed);
                }
                if (!plus_pass || !minus_pass) {
#pragma omp critical
                    failures.push_back(EqualFailure{cases[index], result});
                }
                if (result.minus_wronskian_rank < result.minus_domain) {
#pragma omp critical
                    wronskian_residuals.push_back(
                        EqualFailure{cases[index], result}
                    );
                }
            }
            for (const EqualFailure& failure : failures) {
                std::cout << "TAU_EVEN_EQUAL_DIAGONAL_FAILURE labels=";
                print_degree(failure.labels);
                std::cout << " plus=" << failure.result.plus_rank << '/'
                          << failure.result.plus_domain
                          << " plus_terminal_tensor="
                          << failure.result.plus_terminal_tensor_rank << '/'
                          << failure.result.plus_domain
                          << " plus_diagonal_v2="
                          << failure.result.plus_diagonal_v2
                          << " minus=" << failure.result.minus_rank << '/'
                          << failure.result.minus_domain
                          << " minus_terminal_tensor="
                          << failure.result.minus_terminal_tensor_rank << '/'
                          << failure.result.minus_domain
                          << " minus_wronskian="
                          << failure.result.minus_wronskian_rank << '/'
                          << failure.result.minus_domain
                          << " antisymmetric_divisible="
                          << failure.result.antisymmetric_divisible << '\n';
            }
            for (const EqualFailure& residual : wronskian_residuals) {
                std::cout << "TAU_EVEN_EQUAL_WRONSKIAN_RESIDUAL labels=";
                print_degree(residual.labels);
                std::cout << " rank=" << residual.result.minus_wronskian_rank
                          << '/' << residual.result.minus_domain << '\n';
            }
            std::cout
                << "TAU_EVEN_EQUAL_DIAGONAL_SWEEP_DONE maximum=" << maximum
                << " candidate_cases=" << cases.size()
                << " nonempty=" << nonempty.load()
                << " plus_full=" << plus_full.load()
                << " minus_expected=" << minus_expected.load()
                << " failures=" << failures.size()
                << " wronskian_residuals=" << wronskian_residuals.size()
                << " threads=" << omp_get_max_threads() << '\n';
            return failures.empty() ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        if (argc == 3
            && std::string(argv[1]) == "--tau-even-v2-residual-sweep") {
            const long parsed_maximum = std::strtol(argv[2], nullptr, 10);
            if (parsed_maximum < 4
                || parsed_maximum > std::numeric_limits<int>::max()) {
                throw std::runtime_error(
                    "invalid tau-even-v2-residual-sweep maximum"
                );
            }
            const int maximum = static_cast<int>(parsed_maximum);
            std::size_t nonempty = 0U;
            std::size_t exact_merge_full = 0U;
            std::size_t v2_union_full = 0U;
            std::size_t exact_residual = 0U;
            std::size_t v2_symbolic_uncovered = 0U;
            std::array<std::size_t, 8U> v2_choice_counts{};
            std::array<std::size_t, 8U> first_v2_choice_counts{};
            const bool collect_v2_choices
                = std::getenv("PLUCKER_V2_CHOICE_SUMMARY") != nullptr;
            for (int q = 3; q <= maximum; q += 2) {
                for (int r = q + 2; r <= maximum; r += 2) {
                    for (int first = 2; first <= maximum; first += 2) {
                        for (int second = first; second <= maximum;
                             second += 2) {
                            for (int third = second; third <= maximum;
                                 third += 2) {
                                const std::vector<int> labels{
                                    q, r, 1, 1, first, second, third
                                };
                                const std::vector<Graph> graphs
                                    = tau_even_merged_graphs(labels);
                                if (graphs.empty()) {
                                    continue;
                                }
                                ++nonempty;
                                bool exact_merge = false;
                                bool v2_union = false;
                                std::array<bool, 8U> v2_choices{};
                                for (std::size_t root = 2U;
                                     root < 6U; ++root) {
                                    exact_merge = exact_merge
                                        || binary_terminal_merge_jet_separates(
                                               graphs, 0U, 1U, root
                                           )
                                        || binary_terminal_merge_jet_separates(
                                               graphs, 1U, 0U, root
                                           );
                                    const std::size_t q_choice = root - 2U;
                                    const std::size_t r_choice = root + 2U;
                                    const bool q_v2
                                        = binary_terminal_v2_jet_separates(
                                              graphs, 0U, 1U, root
                                          )
                                        || binary_terminal_v2_merge_jet_separates(
                                               graphs, 0U, 1U, root
                                           );
                                    const bool r_v2
                                        = binary_terminal_v2_jet_separates(
                                              graphs, 1U, 0U, root
                                          )
                                        || binary_terminal_v2_merge_jet_separates(
                                               graphs, 1U, 0U, root
                                           );
                                    v2_choices[q_choice] = q_v2;
                                    v2_choices[r_choice] = r_v2;
                                    v2_union = v2_union || q_v2 || r_v2;
                                }
                                bool recorded_first_v2 = false;
                                for (std::size_t choice = 0U;
                                     choice < v2_choices.size(); ++choice) {
                                    if (!v2_choices[choice]) {
                                        continue;
                                    }
                                    ++v2_choice_counts[choice];
                                    if (!recorded_first_v2) {
                                        recorded_first_v2 = true;
                                        ++first_v2_choice_counts[choice];
                                        if (choice != 0U
                                            && std::getenv(
                                                   "PLUCKER_V2_NONDEFAULT_VERBOSE"
                                               ) != nullptr) {
                                            std::cout
                                                << "TAU_EVEN_V2_CERT labels=";
                                            print_degree(labels);
                                            std::cout << " choice=" << choice
                                                      << '\n';
                                        }
                                    }
                                }
                                if (exact_merge) {
                                    ++exact_merge_full;
                                }
                                if (v2_union) {
                                    ++v2_union_full;
                                }
                                if (exact_merge && !v2_union) {
                                    ++exact_residual;
                                    std::cout
                                        << "TAU_EVEN_V2_RESIDUAL labels=";
                                    print_degree(labels);
                                    std::cout << " sources=" << graphs.size()
                                              << '\n';
                                }
                                const int adjacent_neighbors
                                    = static_cast<int>(first == q + 1)
                                    + static_cast<int>(second == q + 1)
                                    + static_cast<int>(third == q + 1);
                                const bool circuit_family
                                    = first == 2 && second == q + 1
                                    && third == r + 1;
                                const bool adjacent_two_neighbor_family
                                    = r == q + 2
                                    && adjacent_neighbors >= 2;
                                const bool cubic_repeated_neighbor_family
                                    = q == 3 && first == second
                                    && second == third
                                    && (first == r - 1 || first == r + 1);
                                const bool cubic_special_family
                                    = labels
                                      == std::vector<int>{
                                             3, 5, 1, 1, 4, 4, 8
                                         };
                                const bool symbolic_family
                                    = circuit_family
                                    || adjacent_two_neighbor_family
                                    || cubic_repeated_neighbor_family
                                    || cubic_special_family;
                                if (!v2_union && !symbolic_family) {
                                    ++v2_symbolic_uncovered;
                                    std::cout
                                        << "TAU_EVEN_V2_UNCOVERED labels=";
                                    print_degree(labels);
                                    std::cout << " sources=" << graphs.size()
                                              << '\n';
                                }
                            }
                        }
                    }
                }
            }
            std::cout
                << "TAU_EVEN_V2_RESIDUAL_DONE maximum=" << maximum
                << " nonempty=" << nonempty
                << " exact_merge_full=" << exact_merge_full
                << " v2_union_full=" << v2_union_full
                << " exact_residual=" << exact_residual
                << " v2_symbolic_uncovered="
                << v2_symbolic_uncovered << '\n';
            if (collect_v2_choices) {
                std::cout
                    << "TAU_EVEN_V2_CHOICES_qv_qe1_qe2_qe3_rv_re1_re2_re3=";
                for (std::size_t choice = 0U;
                     choice < v2_choice_counts.size(); ++choice) {
                    if (choice != 0U) {
                        std::cout << ',';
                    }
                    std::cout << v2_choice_counts[choice];
                }
                std::cout
                    << " first_qv_qe1_qe2_qe3_rv_re1_re2_re3=";
                for (std::size_t choice = 0U;
                     choice < first_v2_choice_counts.size(); ++choice) {
                    if (choice != 0U) {
                        std::cout << ',';
                    }
                    std::cout << first_v2_choice_counts[choice];
                }
                std::cout << '\n';
            }
            return EXIT_SUCCESS;
        }
        if (argc == 3
            && std::string(argv[1]) == "--tau-even-binary-rank-sweep") {
            const long parsed_maximum = std::strtol(argv[2], nullptr, 10);
            if (parsed_maximum < 4
                || parsed_maximum > std::numeric_limits<int>::max()) {
                throw std::runtime_error(
                    "invalid tau-even-binary-rank-sweep maximum"
                );
            }
            const int maximum = static_cast<int>(parsed_maximum);
            std::map<int, std::array<std::size_t, 3>> by_q;
            std::size_t nonempty = 0U;
            std::size_t q_full = 0U;
            std::size_t r_full = 0U;
            std::size_t either_full = 0U;
            std::size_t tensor_full = 0U;
            std::size_t jet_full = 0U;
            std::size_t v_jet_full = 0U;
            std::size_t hall_jet_full = 0U;
            std::size_t lovett_jet_full = 0U;
            std::size_t merge_jet_full = 0U;
            std::size_t lovett_merge_jet_full = 0U;
            std::size_t merge_lovett_choice_mismatches = 0U;
            std::size_t v2_merge_jet_full = 0U;
            std::size_t merge_v2_choice_mismatches = 0U;
            std::size_t v2_jet_full = 0U;
            std::size_t v2_or_merge_jet_full = 0U;
            std::size_t hall_rank_choice_mismatches = 0U;
            std::size_t double_valuation_full = 0U;
            std::size_t classification_failures = 0U;
            std::array<std::size_t, 8U> jet_choice_counts{};
            std::array<std::size_t, 8U> merge_choice_counts{};
            std::array<std::size_t, 8U> first_merge_choice_counts{};
            std::map<std::vector<int>, std::size_t> merge_pattern_counts;
            const bool collect_merge_patterns
                = std::getenv("PLUCKER_MERGE_PATTERN_SUMMARY") != nullptr;
            for (int q = 3; q <= maximum; q += 2) {
                for (int r = q + 2; r <= maximum; r += 2) {
                    for (int first = 2; first <= maximum; first += 2) {
                        for (int second = first; second <= maximum;
                             second += 2) {
                            for (int third = second; third <= maximum;
                                 third += 2) {
                                const std::vector<int> labels{
                                    q, r, 1, 1, first, second, third
                                };
                                const std::vector<Graph> graphs
                                    = tau_even_merged_graphs(labels);
                                if (graphs.empty()) {
                                    continue;
                                }
                                ++nonempty;
                                const std::size_t q_rank
                                    = binary_vertex_factor_rank(graphs, 0U);
                                const std::size_t r_rank
                                    = binary_vertex_factor_rank(graphs, 1U);
                                const std::size_t tensor_rank
                                    = binary_terminal_tensor_rank(
                                        graphs, 0U, 1U
                                    );
                                bool jet_separated = false;
                                bool v_jet_separated = false;
                                bool hall_jet_separated = false;
                                bool lovett_jet_separated = false;
                                bool merge_jet_separated = false;
                                bool lovett_merge_jet_separated = false;
                                bool v2_merge_jet_separated = false;
                                bool v2_jet_separated = false;
                                bool double_valuation_separated = false;
                                std::array<bool, 8U> merge_choices{};
                                std::array<std::vector<int>, 8U>
                                    merge_witnesses;
                                for (std::size_t root = 2U;
                                     root < 6U; ++root) {
                                    const bool q_jet
                                        = binary_terminal_jet_separates(
                                              graphs, 0U, 1U, root
                                          );
                                    const bool r_jet
                                        = binary_terminal_jet_separates(
                                              graphs, 1U, 0U, root
                                          );
                                    const bool q_hall
                                        = binary_terminal_hall_jet_separates(
                                              graphs, 0U, 1U, root
                                          );
                                    const bool r_hall
                                        = binary_terminal_hall_jet_separates(
                                              graphs, 1U, 0U, root
                                          );
                                    lovett_jet_separated
                                        = lovett_jet_separated
                                        || binary_terminal_lovett_jet_separates(
                                               graphs, 0U, 1U, root
                                           )
                                        || binary_terminal_lovett_jet_separates(
                                               graphs, 1U, 0U, root
                                           );
                                    const std::size_t q_choice = root - 2U;
                                    const std::size_t r_choice = root + 2U;
                                    merge_choices[q_choice]
                                        = binary_terminal_merge_jet_separates(
                                              graphs, 0U, 1U, root,
                                              &merge_witnesses[q_choice]
                                          );
                                    merge_choices[r_choice]
                                        = binary_terminal_merge_jet_separates(
                                              graphs, 1U, 0U, root,
                                              &merge_witnesses[r_choice]
                                          );
                                    const bool q_lovett_merge
                                        = binary_terminal_lovett_merge_jet_separates(
                                              graphs, 0U, 1U, root
                                          );
                                    const bool r_lovett_merge
                                        = binary_terminal_lovett_merge_jet_separates(
                                              graphs, 1U, 0U, root
                                          );
                                    const bool q_v2_merge
                                        = binary_terminal_v2_merge_jet_separates(
                                              graphs, 0U, 1U, root
                                          );
                                    const bool r_v2_merge
                                        = binary_terminal_v2_merge_jet_separates(
                                              graphs, 1U, 0U, root
                                          );
                                    const bool q_v2
                                        = binary_terminal_v2_jet_separates(
                                              graphs, 0U, 1U, root
                                          );
                                    const bool r_v2
                                        = binary_terminal_v2_jet_separates(
                                              graphs, 1U, 0U, root
                                          );
                                    if (q_lovett_merge
                                            != merge_choices[q_choice]
                                        || r_lovett_merge
                                            != merge_choices[r_choice]) {
                                        ++merge_lovett_choice_mismatches;
                                    }
                                    if (q_v2_merge
                                            != merge_choices[q_choice]
                                        || r_v2_merge
                                            != merge_choices[r_choice]) {
                                        ++merge_v2_choice_mismatches;
                                    }
                                    lovett_merge_jet_separated
                                        = lovett_merge_jet_separated
                                        || q_lovett_merge || r_lovett_merge;
                                    v2_merge_jet_separated
                                        = v2_merge_jet_separated
                                        || q_v2_merge || r_v2_merge;
                                    v2_jet_separated
                                        = v2_jet_separated || q_v2 || r_v2;
                                    merge_jet_separated
                                        = merge_jet_separated
                                        || merge_choices[q_choice]
                                        || merge_choices[r_choice];
                                    if (q_hall != q_jet || r_hall != r_jet) {
                                        ++hall_rank_choice_mismatches;
                                    }
                                    hall_jet_separated = hall_jet_separated
                                        || q_hall || r_hall;
                                    if (q_jet) {
                                        ++jet_choice_counts[root - 2U];
                                    }
                                    if (r_jet) {
                                        ++jet_choice_counts[4U + root - 2U];
                                    }
                                    if (root == 2U) {
                                        v_jet_separated = q_jet || r_jet;
                                    }
                                    if (!jet_separated) {
                                        jet_separated = q_jet || r_jet;
                                    }
                                    if (!double_valuation_separated) {
                                        double_valuation_separated
                                            = binary_terminal_double_valuation_separates(
                                                  graphs, 0U, 1U, root
                                              )
                                            || binary_terminal_double_valuation_separates(
                                                  graphs, 1U, 0U, root
                                              );
                                    }
                                }
                                bool recorded_first_merge = false;
                                for (std::size_t choice = 0U;
                                     choice < merge_choices.size(); ++choice) {
                                    if (!merge_choices[choice]) {
                                        continue;
                                    }
                                    ++merge_choice_counts[choice];
                                    if (recorded_first_merge) {
                                        continue;
                                    }
                                    recorded_first_merge = true;
                                    ++first_merge_choice_counts[choice];
                                    const std::vector<int>& witness
                                        = merge_witnesses[choice];
                                    if (collect_merge_patterns) {
                                        std::vector<int> pattern{
                                            static_cast<int>(choice)
                                        };
                                        for (std::size_t item = 0U;
                                             item < witness.size();
                                             item += 3U) {
                                            pattern.push_back(
                                                witness[item + 1U]
                                            );
                                            pattern.push_back(
                                                witness[item + 2U]
                                            );
                                        }
                                        ++merge_pattern_counts[pattern];
                                    }
                                    const bool verbose_merge_certificate
                                        = std::getenv(
                                              "PLUCKER_MERGE_CERT_VERBOSE"
                                          ) != nullptr
                                        || (choice != 0U
                                            && std::getenv(
                                                   "PLUCKER_MERGE_NONDEFAULT_VERBOSE"
                                               ) != nullptr);
                                    if (verbose_merge_certificate) {
                                        std::cout
                                            << "TAU_EVEN_MERGE_CERT labels=";
                                        print_degree(labels);
                                        std::cout << " choice=" << choice
                                                  << " witness=";
                                        print_degree(witness);
                                        std::cout << '\n';
                                    }
                                }
                                const bool q_is_full = q_rank == graphs.size();
                                const bool r_is_full = r_rank == graphs.size();
                                if (q_is_full) {
                                    ++q_full;
                                }
                                if (r_is_full) {
                                    ++r_full;
                                }
                                if (q_is_full || r_is_full) {
                                    ++either_full;
                                }
                                if (tensor_rank == graphs.size()) {
                                    ++tensor_full;
                                }
                                if (jet_separated) {
                                    ++jet_full;
                                }
                                if (v_jet_separated) {
                                    ++v_jet_full;
                                } else if (std::getenv(
                                               "PLUCKER_V_JET_VERBOSE"
                                           ) != nullptr) {
                                    std::cout
                                        << "TAU_EVEN_V_JET_DEFECT labels=";
                                    print_degree(labels);
                                    std::cout << " sources=" << graphs.size()
                                              << '\n';
                                }
                                if (hall_jet_separated) {
                                    ++hall_jet_full;
                                }
                                if (lovett_jet_separated) {
                                    ++lovett_jet_full;
                                } else if (std::getenv(
                                               "PLUCKER_LOVETT_VERBOSE"
                                           ) != nullptr) {
                                    std::cout
                                        << "TAU_EVEN_LOVETT_JET_DEFECT labels=";
                                    print_degree(labels);
                                    std::cout << " sources=" << graphs.size()
                                              << '\n';
                                }
                                if (merge_jet_separated) {
                                    ++merge_jet_full;
                                } else if (std::getenv(
                                               "PLUCKER_MERGE_JET_VERBOSE"
                                           ) != nullptr) {
                                    std::cout
                                        << "TAU_EVEN_MERGE_JET_DEFECT labels=";
                                    print_degree(labels);
                                    std::cout << " sources=" << graphs.size()
                                              << '\n';
                                }
                                if (lovett_merge_jet_separated) {
                                    ++lovett_merge_jet_full;
                                }
                                if (v2_merge_jet_separated) {
                                    ++v2_merge_jet_full;
                                }
                                if (v2_jet_separated) {
                                    ++v2_jet_full;
                                }
                                if (v2_jet_separated
                                    || v2_merge_jet_separated) {
                                    ++v2_or_merge_jet_full;
                                }
                                if (merge_jet_separated
                                    && !v2_jet_separated
                                    && !v2_merge_jet_separated
                                    && std::getenv(
                                           "PLUCKER_V2_RESIDUAL_VERBOSE"
                                       ) != nullptr) {
                                    std::cout
                                        << "TAU_EVEN_V2_RESIDUAL labels=";
                                    print_degree(labels);
                                    std::cout << " sources=" << graphs.size()
                                              << '\n';
                                }
                                if (double_valuation_separated) {
                                    ++double_valuation_full;
                                }
                                const bool exceptional_family
                                    = first == 2 && second == q + 1
                                    && third == r + 1;
                                const std::size_t expected_exceptional_sources
                                    = 4U
                                    + (r == q + 2 ? 1U : 0U)
                                    + (q == 3 ? 2U : 0U);
                                const bool classification_valid
                                    = (tensor_rank != graphs.size())
                                          == exceptional_family
                                    && (!exceptional_family
                                        || (graphs.size()
                                                == expected_exceptional_sources
                                            && tensor_rank + 1U
                                                == graphs.size()));
                                if (!classification_valid) {
                                    ++classification_failures;
                                    std::cout
                                        << "TAU_EVEN_BINARY_CLASSIFICATION_FAIL labels=";
                                    print_degree(labels);
                                    std::cout << " sources=" << graphs.size()
                                              << " tensor_rank="
                                              << tensor_rank << '\n';
                                }
                                std::array<std::size_t, 3>& counts = by_q[q];
                                ++counts[0];
                                if (q_is_full) {
                                    ++counts[1];
                                }
                                if (q_is_full || r_is_full) {
                                    ++counts[2];
                                }
                                if (std::getenv("PLUCKER_BINARY_RANK_VERBOSE")
                                        != nullptr
                                    && (tensor_rank != graphs.size()
                                        || !jet_separated)) {
                                    std::cout
                                        << "TAU_EVEN_BINARY_DIAGNOSTIC labels=";
                                    print_degree(labels);
                                    std::cout << " sources=" << graphs.size()
                                              << " q_rank=" << q_rank
                                              << " r_rank=" << r_rank
                                              << " tensor_rank="
                                              << tensor_rank
                                              << " jet_separated="
                                              << jet_separated
                                              << " double_valuation="
                                              << double_valuation_separated
                                              << '\n';
                                }
                            }
                        }
                    }
                }
            }
            for (const auto& [q, counts] : by_q) {
                std::cout << "TAU_EVEN_BINARY_RANK_Q q=" << q
                          << " cases=" << counts[0]
                          << " q_full=" << counts[1]
                          << " either_full=" << counts[2] << '\n';
            }
            std::cout << "TAU_EVEN_BINARY_RANK_DONE maximum=" << maximum
                      << " nonempty=" << nonempty
                      << " q_full=" << q_full << " r_full=" << r_full
                      << " either_full=" << either_full
                      << " tensor_full=" << tensor_full
                      << " jet_full=" << jet_full
                      << " v_jet_full=" << v_jet_full
                      << " hall_jet_full=" << hall_jet_full
                      << " lovett_jet_full=" << lovett_jet_full
                      << " merge_jet_full=" << merge_jet_full
                      << " lovett_merge_jet_full="
                      << lovett_merge_jet_full
                      << " merge_lovett_choice_mismatches="
                      << merge_lovett_choice_mismatches
                      << " v2_merge_jet_full=" << v2_merge_jet_full
                      << " merge_v2_choice_mismatches="
                      << merge_v2_choice_mismatches
                      << " v2_jet_full=" << v2_jet_full
                      << " v2_or_merge_jet_full="
                      << v2_or_merge_jet_full
                      << " hall_rank_choice_mismatches="
                      << hall_rank_choice_mismatches
                      << " double_valuation_full="
                      << double_valuation_full
                      << " classification_failures="
                      << classification_failures
                      << " jet_choices_qv_qe1_qe2_qe3_rv_re1_re2_re3=";
            for (std::size_t choice = 0U;
                 choice < jet_choice_counts.size(); ++choice) {
                if (choice != 0U) {
                    std::cout << ',';
                }
                std::cout << jet_choice_counts[choice];
            }
            std::cout << '\n';
            std::cout
                << "TAU_EVEN_MERGE_CHOICES_qv_qe1_qe2_qe3_rv_re1_re2_re3=";
            for (std::size_t choice = 0U;
                 choice < merge_choice_counts.size(); ++choice) {
                if (choice != 0U) {
                    std::cout << ',';
                }
                std::cout << merge_choice_counts[choice];
            }
            std::cout
                << " first_qv_qe1_qe2_qe3_rv_re1_re2_re3=";
            for (std::size_t choice = 0U;
                 choice < first_merge_choice_counts.size(); ++choice) {
                if (choice != 0U) {
                    std::cout << ',';
                }
                std::cout << first_merge_choice_counts[choice];
            }
            std::cout << '\n';
            if (collect_merge_patterns) {
                for (const auto& [pattern, count] : merge_pattern_counts) {
                    std::cout << "TAU_EVEN_MERGE_PATTERN signature=";
                    print_degree(pattern);
                    std::cout << " cases=" << count << '\n';
                }
            }
            return EXIT_SUCCESS;
        }
        if (argc == 3
            && std::string(argv[1]) == "--tau-even-order-cover") {
            const long parsed_maximum = std::strtol(argv[2], nullptr, 10);
            if (parsed_maximum < 4
                || parsed_maximum > std::numeric_limits<int>::max()) {
                throw std::runtime_error(
                    "invalid tau-even-order-cover maximum"
                );
            }
            const int maximum = static_cast<int>(parsed_maximum);
            std::vector<std::vector<int>> labels_by_case;
            std::vector<std::vector<Graph>> graphs_by_case;
            for (int q = 3; q <= maximum; q += 2) {
                for (int r = q + 2; r <= maximum; r += 2) {
                    for (int first = 2; first <= maximum; first += 2) {
                        for (int second = first; second <= maximum;
                             second += 2) {
                            for (int third = second; third <= maximum;
                                 third += 2) {
                                std::vector<int> labels{
                                    q, r, 1, 1, first, second, third
                                };
                                std::vector<Graph> graphs
                                    = tau_even_merged_graphs(labels);
                                if (!graphs.empty()) {
                                    labels_by_case.push_back(std::move(labels));
                                    graphs_by_case.push_back(std::move(graphs));
                                }
                            }
                        }
                    }
                }
            }
            std::vector<std::vector<std::size_t>> orders;
            std::vector<std::size_t> order{0U, 1U, 2U, 3U, 4U, 5U};
            do {
                orders.push_back(order);
            } while (std::next_permutation(order.begin() + 1, order.end()));
            std::vector<std::vector<std::uint8_t>> covers(
                orders.size(),
                std::vector<std::uint8_t>(graphs_by_case.size(), 0U)
            );
#pragma omp parallel for schedule(dynamic)
            for (std::size_t order_index = 0U;
                 order_index < orders.size(); ++order_index) {
                for (std::size_t case_index = 0U;
                     case_index < graphs_by_case.size(); ++case_index) {
                    if (leading_collision_count(
                            graphs_by_case[case_index], orders[order_index]
                        ) == 0U) {
                        covers[order_index][case_index] = 1U;
                    }
                }
            }
            std::vector<std::uint8_t> uncovered(graphs_by_case.size(), 1U);
            std::size_t remaining = graphs_by_case.size();
            std::size_t charts = 0U;
            while (remaining != 0U) {
                std::size_t best_order = orders.size();
                std::size_t best_count = 0U;
                for (std::size_t order_index = 0U;
                     order_index < orders.size(); ++order_index) {
                    std::size_t count = 0U;
                    for (std::size_t case_index = 0U;
                         case_index < graphs_by_case.size(); ++case_index) {
                        if (uncovered[case_index] != 0U
                            && covers[order_index][case_index] != 0U) {
                            ++count;
                        }
                    }
                    if (count > best_count) {
                        best_count = count;
                        best_order = order_index;
                    }
                }
                if (best_count == 0U) {
                    break;
                }
                ++charts;
                std::cout << "TAU_EVEN_ORDER_COVER_CHART index=" << charts
                          << " order=[";
                for (std::size_t position = 0U;
                     position < orders[best_order].size(); ++position) {
                    if (position != 0U) {
                        std::cout << ',';
                    }
                    std::cout << orders[best_order][position];
                }
                std::cout << "] newly_covered=" << best_count << '\n';
                for (std::size_t case_index = 0U;
                     case_index < graphs_by_case.size(); ++case_index) {
                    if (uncovered[case_index] != 0U
                        && covers[best_order][case_index] != 0U) {
                        uncovered[case_index] = 0U;
                        --remaining;
                    }
                }
            }
            for (std::size_t case_index = 0U;
                 case_index < labels_by_case.size(); ++case_index) {
                if (uncovered[case_index] != 0U) {
                    std::cout << "TAU_EVEN_ORDER_COVER_UNCOVERED labels=";
                    print_degree(labels_by_case[case_index]);
                    std::cout << " domain="
                              << graphs_by_case[case_index].size() << '\n';
                }
            }
            std::cout << "TAU_EVEN_ORDER_COVER_DONE maximum=" << maximum
                      << " nonempty_cases=" << graphs_by_case.size()
                      << " charts=" << charts
                      << " uncovered=" << remaining
                      << " threads=" << omp_get_max_threads() << '\n';
            return EXIT_SUCCESS;
        }
        if (argc == 3
            && std::string(argv[1]) == "--tau-even-divisor-sweep") {
            const long parsed_maximum = std::strtol(argv[2], nullptr, 10);
            if (parsed_maximum < 4
                || parsed_maximum > std::numeric_limits<int>::max()) {
                throw std::runtime_error(
                    "invalid tau-even-divisor-sweep maximum"
                );
            }
            const int maximum = static_cast<int>(parsed_maximum);
            std::vector<std::vector<int>> cases;
            for (int q = 3; q <= maximum; q += 2) {
                for (int r = q + 2; r <= maximum; r += 2) {
                    for (int first = 2; first <= maximum; first += 2) {
                        for (int second = first; second <= maximum;
                             second += 2) {
                            for (int third = second; third <= maximum;
                                 third += 2) {
                                cases.push_back({
                                    q, r, 1, 1, first, second, third
                                });
                            }
                        }
                    }
                }
            }
            struct DivisorFailure {
                std::vector<int> labels;
                std::size_t domain = 0U;
                std::size_t standard_collisions = 0U;
                std::vector<std::size_t> standard_order;
            };
            std::atomic<std::size_t> checked{0U};
            std::atomic<std::size_t> nonempty{0U};
            std::array<std::atomic<std::size_t>, 36U> root_counts{};
            std::array<std::atomic<std::size_t>, 3U> chart_counts{};
            std::array<std::atomic<std::size_t>, 4U>
                partition_extremal_counts{};
            std::vector<DivisorFailure> failures;
            std::vector<DivisorFailure> raw_failures;
            std::vector<DivisorFailure> tree_failures;
            std::vector<DivisorFailure> nonterminal_root_cases;
            std::vector<DivisorFailure> three_chart_failures;
            std::vector<DivisorFailure> second_chart_cases;
            std::vector<DivisorFailure> third_chart_cases;
            std::vector<DivisorFailure> chart_classification_failures;
            std::vector<DivisorFailure> chosen_chart_nonstandard_cases;
            std::vector<DivisorFailure> partition_chart_failures;
            std::vector<DivisorFailure> partition_raw_collisions;
            std::vector<DivisorFailure> partition_least_only_cases;
            std::vector<DivisorFailure> separated_normal_form_failures;
            const bool audit_raw
                = std::getenv("PLUCKER_RAW_AUDIT") != nullptr;
            const bool audit_tree
                = std::getenv("PLUCKER_TREE_AUDIT") != nullptr;
            const bool audit_three_charts
                = std::getenv("PLUCKER_THREE_CHART_AUDIT") != nullptr;
#pragma omp parallel for schedule(dynamic)
            for (std::size_t index = 0U; index < cases.size(); ++index) {
                const std::vector<Graph> graphs
                    = tau_even_merged_graphs(cases[index]);
                if (!graphs.empty()) {
                    nonempty.fetch_add(1U, std::memory_order_relaxed);
                }
                const auto [divisor_separated, root]
                    = first_divisor_separating_edge(graphs);
                if (divisor_separated && graphs.size() > 1U) {
                    root_counts[root.first * 6U + root.second].fetch_add(
                        1U, std::memory_order_relaxed
                    );
                    if (std::getenv("PLUCKER_ROOT_VERBOSE") != nullptr
                        && root != std::pair<std::size_t, std::size_t>{0U, 1U}) {
#pragma omp critical
                        {
                            nonterminal_root_cases.push_back(DivisorFailure{
                                cases[index], graphs.size(), 0U,
                                {root.first, root.second}
                            });
                        }
                    }
                }
                if (!divisor_separated) {
                    const auto [collisions, order]
                        = best_leading_collision_count(graphs);
#pragma omp critical
                    {
                        failures.push_back(DivisorFailure{
                            cases[index], graphs.size(), collisions, order
                        });
                    }
                }
                if (audit_raw && !graphs.empty()) {
                    const auto [collisions, order]
                        = best_raw_leading_collision_count(graphs);
                    if (collisions != 0U) {
#pragma omp critical
                        {
                            raw_failures.push_back(DivisorFailure{
                                cases[index], graphs.size(), collisions, order
                            });
                        }
                    }
                }
                if (audit_tree && !graphs.empty()) {
                    const auto [collisions, tree]
                        = best_tree_collision_count(graphs);
                    if (collisions != 0U) {
                        std::vector<std::size_t> encoded_tree;
                        encoded_tree.reserve(tree.size());
                        for (const unsigned int split : tree) {
                            encoded_tree.push_back(
                                static_cast<std::size_t>(split)
                            );
                        }
#pragma omp critical
                        {
                            tree_failures.push_back(DivisorFailure{
                                cases[index], graphs.size(), collisions,
                                std::move(encoded_tree)
                            });
                        }
                    }
                }
                if (audit_three_charts && !graphs.empty()) {
                    const std::array<std::vector<std::size_t>, 3> charts{
                        std::vector<std::size_t>{0U, 3U, 4U, 1U, 5U, 2U},
                        std::vector<std::size_t>{0U, 1U, 5U, 3U, 2U, 4U},
                        std::vector<std::size_t>{0U, 5U, 1U, 3U, 4U, 2U}
                    };
                    const bool cubic_exception = cases[index]
                        == std::vector<int>{3, 5, 1, 1, 4, 4, 4};
                    const bool double_lower
                        = cases[index][4] == cases[index][0] - 1
                        && cases[index][5] == cases[index][0] - 1;
                    const bool lower_wall_exception = double_lower
                        && cases[index][6] == cases[index][0] - 1
                        && cases[index][1] == 2 * cases[index][0] - 1;
                    const std::size_t partition_chart = cubic_exception
                        ? 2U
                        : (double_lower && !lower_wall_exception ? 1U : 0U);
                    const std::size_t partition_raw_collision_count
                        = raw_leading_collision_count(
                            graphs, charts[partition_chart]
                        );
                    if (partition_raw_collision_count != 0U) {
#pragma omp critical
                        {
                            partition_raw_collisions.push_back(
                                DivisorFailure{
                                    cases[index], graphs.size(),
                                    partition_raw_collision_count, {}
                                }
                            );
                        }
                    }
                    const auto [least_collisions, greatest_collisions]
                        = extremal_collision_counts(
                            graphs, charts[partition_chart]
                        );
                    bool separated_normal_form_matches = true;
                    std::map<Key, Support> partition_support_memo;
                    for (const Graph& graph : graphs) {
                        const Graph permuted = permute_graph(
                            graph, charts[partition_chart]
                        );
                        const Support& support = straighten_support(
                            permuted, partition_support_memo
                        );
                        if (flatten(separated_normal_form(permuted))
                            != *support.rbegin()) {
                            separated_normal_form_matches = false;
                            break;
                        }
                    }
                    if (!separated_normal_form_matches) {
#pragma omp critical
                        {
                            separated_normal_form_failures.push_back(
                                DivisorFailure{
                                    cases[index], graphs.size(),
                                    partition_chart, {}
                                }
                            );
                        }
                    }
                    const std::size_t extremal_kind = least_collisions == 0U
                        ? (greatest_collisions == 0U ? 0U : 1U)
                        : (greatest_collisions == 0U ? 2U : 3U);
                    partition_extremal_counts[extremal_kind].fetch_add(
                        1U, std::memory_order_relaxed
                    );
                    if (extremal_kind == 1U) {
#pragma omp critical
                        {
                            partition_least_only_cases.push_back(
                                DivisorFailure{
                                    cases[index], graphs.size(), 0U, {}
                                }
                            );
                        }
                    }
                    if (least_collisions != 0U
                        && greatest_collisions != 0U) {
#pragma omp critical
                        {
                            partition_chart_failures.push_back(
                                DivisorFailure{
                                    cases[index], graphs.size(),
                                    partition_chart, {}
                                }
                            );
                        }
                    }
                    bool covered = false;
                    for (std::size_t chart_index = 0U;
                         chart_index < charts.size(); ++chart_index) {
                        if (leading_collision_count(
                                graphs, charts[chart_index]
                            ) == 0U) {
                            covered = true;
                            bool all_noncrossing = true;
                            for (const Graph& graph : graphs) {
                                if (first_crossing(permute_graph(
                                        graph, charts[chart_index]
                                    )).first) {
                                    all_noncrossing = false;
                                    break;
                                }
                            }
                            if (!all_noncrossing) {
#pragma omp critical
                                {
                                    chosen_chart_nonstandard_cases.push_back(
                                        DivisorFailure{
                                            cases[index], graphs.size(),
                                            chart_index, {}
                                        }
                                    );
                                }
                            }
                            chart_counts[chart_index].fetch_add(
                                1U, std::memory_order_relaxed
                            );
                            if (chart_index == 1U
                                && std::getenv(
                                    "PLUCKER_THREE_CHART_VERBOSE"
                                ) != nullptr) {
#pragma omp critical
                                {
                                    second_chart_cases.push_back(
                                        DivisorFailure{
                                            cases[index], graphs.size(), 0U, {}
                                        }
                                    );
                                }
                            }
                            if (chart_index == 1U
                                && !(cases[index][4] == cases[index][0] - 1
                                    && cases[index][5]
                                        == cases[index][0] - 1)) {
#pragma omp critical
                                {
                                    chart_classification_failures.push_back(
                                        DivisorFailure{
                                            cases[index], graphs.size(),
                                            chart_index, {}
                                        }
                                    );
                                }
                            }
                            if (chart_index == 2U
                                && std::getenv(
                                    "PLUCKER_THREE_CHART_VERBOSE"
                                ) != nullptr) {
#pragma omp critical
                                {
                                    third_chart_cases.push_back(
                                        DivisorFailure{
                                            cases[index], graphs.size(), 0U, {}
                                        }
                                    );
                                }
                            }
                            if (chart_index == 2U
                                && cases[index]
                                    != std::vector<int>{3, 5, 1, 1, 4, 4, 4}) {
#pragma omp critical
                                {
                                    chart_classification_failures.push_back(
                                        DivisorFailure{
                                            cases[index], graphs.size(),
                                            chart_index, {}
                                        }
                                    );
                                }
                            }
                            break;
                        }
                    }
                    if (!covered) {
#pragma omp critical
                        {
                            three_chart_failures.push_back(DivisorFailure{
                                cases[index], graphs.size(), 0U, {}
                            });
                        }
                    }
                }
                const std::size_t completed = checked.fetch_add(
                    1U, std::memory_order_relaxed
                ) + 1U;
                if (cases.size() > 10000U && completed % 5000U == 0U) {
#pragma omp critical
                    {
                        std::cout << "TAU_EVEN_DIVISOR_PROGRESS completed="
                                  << completed << '/' << cases.size() << '\n';
                        std::cout.flush();
                    }
                }
            }
            std::sort(
                failures.begin(), failures.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                raw_failures.begin(), raw_failures.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                tree_failures.begin(), tree_failures.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                nonterminal_root_cases.begin(), nonterminal_root_cases.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                three_chart_failures.begin(), three_chart_failures.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                second_chart_cases.begin(), second_chart_cases.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                third_chart_cases.begin(), third_chart_cases.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                chart_classification_failures.begin(),
                chart_classification_failures.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                chosen_chart_nonstandard_cases.begin(),
                chosen_chart_nonstandard_cases.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                partition_chart_failures.begin(),
                partition_chart_failures.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                partition_raw_collisions.begin(),
                partition_raw_collisions.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                partition_least_only_cases.begin(),
                partition_least_only_cases.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            std::sort(
                separated_normal_form_failures.begin(),
                separated_normal_form_failures.end(),
                [](const DivisorFailure& first,
                   const DivisorFailure& second) {
                    return first.labels < second.labels;
                }
            );
            for (const DivisorFailure& failure : failures) {
                std::cout << "TAU_EVEN_DIVISOR_UNSEPARATED labels=";
                print_degree(failure.labels);
                std::cout << " domain=" << failure.domain
                          << " standard_collisions="
                          << failure.standard_collisions
                          << " standard_order=[";
                for (std::size_t position = 0U;
                     position < failure.standard_order.size(); ++position) {
                    if (position != 0U) {
                        std::cout << ',';
                    }
                    std::cout << failure.standard_order[position];
                }
                std::cout << "]\n";
            }
            for (const DivisorFailure& failure : raw_failures) {
                std::cout << "TAU_EVEN_RAW_UNSEPARATED labels=";
                print_degree(failure.labels);
                std::cout << " domain=" << failure.domain
                          << " best_collisions="
                          << failure.standard_collisions
                          << " best_order=[";
                for (std::size_t position = 0U;
                     position < failure.standard_order.size(); ++position) {
                    if (position != 0U) {
                        std::cout << ',';
                    }
                    std::cout << failure.standard_order[position];
                }
                std::cout << "]\n";
            }
            for (const DivisorFailure& failure : tree_failures) {
                std::cout << "TAU_EVEN_TREE_UNSEPARATED labels=";
                print_degree(failure.labels);
                std::cout << " domain=" << failure.domain
                          << " best_collisions="
                          << failure.standard_collisions << '\n';
            }
            for (const DivisorFailure& root_case : nonterminal_root_cases) {
                std::cout << "TAU_EVEN_DIVISOR_NONTERMINAL_ROOT labels=";
                print_degree(root_case.labels);
                std::cout << " domain=" << root_case.domain
                          << " root=" << root_case.standard_order[0]
                          << '-' << root_case.standard_order[1] << '\n';
            }
            for (const DivisorFailure& failure : three_chart_failures) {
                std::cout << "TAU_EVEN_THREE_CHART_UNCOVERED labels=";
                print_degree(failure.labels);
                std::cout << " domain=" << failure.domain << '\n';
            }
            for (const DivisorFailure& third_case : third_chart_cases) {
                std::cout << "TAU_EVEN_THIRD_CHART_CASE labels=";
                print_degree(third_case.labels);
                std::cout << " domain=" << third_case.domain << '\n';
            }
            for (const DivisorFailure& second_case : second_chart_cases) {
                std::cout << "TAU_EVEN_SECOND_CHART_CASE labels=";
                print_degree(second_case.labels);
                std::cout << " domain=" << second_case.domain << '\n';
            }
            for (const DivisorFailure& failure :
                 chart_classification_failures) {
                std::cout << "TAU_EVEN_CHART_CLASSIFICATION_FAIL labels=";
                print_degree(failure.labels);
                std::cout << " domain=" << failure.domain
                          << " chart_index="
                          << failure.standard_collisions << '\n';
            }
            if (std::getenv("PLUCKER_THREE_CHART_VERBOSE") != nullptr) {
                for (const DivisorFailure& failure :
                     chosen_chart_nonstandard_cases) {
                    std::cout << "TAU_EVEN_CHOSEN_CHART_NONSTANDARD labels=";
                    print_degree(failure.labels);
                    std::cout << " domain=" << failure.domain
                              << " chart_index="
                              << failure.standard_collisions << '\n';
                }
            }
            for (const DivisorFailure& failure : partition_chart_failures) {
                std::cout << "TAU_EVEN_PARTITION_CHART_FAIL labels=";
                print_degree(failure.labels);
                std::cout << " domain=" << failure.domain
                          << " chart_index="
                          << failure.standard_collisions << '\n';
            }
            if (std::getenv("PLUCKER_THREE_CHART_VERBOSE") != nullptr) {
                for (const DivisorFailure& failure :
                     partition_raw_collisions) {
                    std::cout << "TAU_EVEN_PARTITION_RAW_COLLISION labels=";
                    print_degree(failure.labels);
                    std::cout << " domain=" << failure.domain
                              << " collisions="
                              << failure.standard_collisions << '\n';
                }
            }
            std::size_t maximum_partition_raw_collisions = 0U;
            for (const DivisorFailure& failure : partition_raw_collisions) {
                maximum_partition_raw_collisions = std::max(
                    maximum_partition_raw_collisions,
                    failure.standard_collisions
                );
            }
            for (const DivisorFailure& least_case :
                 partition_least_only_cases) {
                std::cout << "TAU_EVEN_PARTITION_LEAST_ONLY labels=";
                print_degree(least_case.labels);
                std::cout << " domain=" << least_case.domain << '\n';
            }
            if (std::getenv("PLUCKER_THREE_CHART_VERBOSE") != nullptr) {
                for (const DivisorFailure& failure :
                     separated_normal_form_failures) {
                    std::cout << "TAU_EVEN_SEPARATED_NORMAL_FORM_FAIL labels=";
                    print_degree(failure.labels);
                    std::cout << " domain=" << failure.domain
                              << " chart_index="
                              << failure.standard_collisions << '\n';
                }
            }
            std::cout << "TAU_EVEN_DIVISOR_ROOTS";
            for (std::size_t first = 0U; first < 6U; ++first) {
                for (std::size_t second = first + 1U; second < 6U; ++second) {
                    const std::size_t count
                        = root_counts[first * 6U + second].load();
                    if (count != 0U) {
                        std::cout << ' ' << first << '-' << second
                                  << ':' << count;
                    }
                }
            }
            std::cout << '\n';
            std::cout << "TAU_EVEN_DIVISOR_SWEEP_PASS maximum=" << maximum
                      << " cases=" << checked.load()
                      << " nonempty=" << nonempty.load()
                      << " unseparated=" << failures.size()
                      << " raw_unseparated=" << raw_failures.size()
                      << " tree_unseparated=" << tree_failures.size()
                      << " three_chart_uncovered="
                      << three_chart_failures.size()
                      << " chart_counts=" << chart_counts[0].load() << ','
                      << chart_counts[1].load() << ','
                      << chart_counts[2].load()
                      << " chart_classification_failures="
                      << chart_classification_failures.size()
                      << " chosen_chart_nonstandard="
                      << chosen_chart_nonstandard_cases.size()
                      << " partition_chart_failures="
                      << partition_chart_failures.size()
                      << " partition_raw_collision_cases="
                      << partition_raw_collisions.size()
                      << " partition_raw_collision_max="
                      << maximum_partition_raw_collisions
                      << " partition_extrema_both_least_greatest_neither="
                      << partition_extremal_counts[0].load() << ','
                      << partition_extremal_counts[1].load() << ','
                      << partition_extremal_counts[2].load() << ','
                      << partition_extremal_counts[3].load()
                      << " separated_normal_form_failures="
                      << separated_normal_form_failures.size()
                      << " threads=" << omp_get_max_threads() << '\n';
            return EXIT_SUCCESS;
        }
        if (argc == 3 && std::string(argv[1]) == "--tau-even-sweep") {
            const long parsed_maximum = std::strtol(argv[2], nullptr, 10);
            if (parsed_maximum < 3
                || parsed_maximum > std::numeric_limits<int>::max()) {
                throw std::runtime_error("invalid tau-even-sweep maximum");
            }
            const int maximum = static_cast<int>(parsed_maximum);
            std::vector<std::vector<int>> cases;
            for (int q = 3; q <= maximum; q += 2) {
                for (int r = q; r <= maximum; r += 2) {
                    for (int first = 2; first <= maximum; first += 2) {
                        for (int second = first; second <= maximum;
                             second += 2) {
                            for (int third = second; third <= maximum;
                                 third += 2) {
                                cases.push_back({
                                    q, r, 1, 1, first, second, third
                                });
                            }
                        }
                    }
                }
            }
            std::atomic<bool> pass{true};
            std::atomic<std::size_t> checked{0U};
            std::atomic<std::size_t> unseparated{0U};
            std::atomic<std::size_t> unseparated_unequal{0U};
            std::vector<std::pair<std::vector<int>, TauEvenResult>>
                kernel_cases;
            std::vector<std::vector<int>> unseparated_cases;
            std::map<std::vector<std::size_t>, std::size_t>
                unequal_standard_orders;
            std::vector<std::vector<int>> unequal_standard_failures;
            std::vector<int> failure_labels;
            TauEvenResult failure_result;
#pragma omp parallel for schedule(dynamic)
            for (std::size_t index = 0U; index < cases.size(); ++index) {
                if (!pass.load(std::memory_order_relaxed)) {
                    continue;
                }
                const TauEvenResult result = analyze_tau_even_allocation(
                    cases[index]
                );
                const std::size_t completed = checked.fetch_add(
                    1U, std::memory_order_relaxed
                ) + 1U;
                if (cases.size() > 1000U && completed % 250U == 0U) {
#pragma omp critical
                    {
                        std::cout << "TAU_EVEN_PROGRESS completed="
                                  << completed << '/' << cases.size() << '\n';
                        std::cout.flush();
                    }
                }
                if (!result.divisor_separated) {
                    unseparated.fetch_add(1U, std::memory_order_relaxed);
                    if (cases[index][0] != cases[index][1]) {
                        unseparated_unequal.fetch_add(
                            1U, std::memory_order_relaxed
                        );
                    }
#pragma omp critical
                    {
                        unseparated_cases.push_back(cases[index]);
                    }
                }
                if (result.global_kernel != 0U) {
#pragma omp critical
                    {
                        kernel_cases.emplace_back(cases[index], result);
                    }
                }
                if (cases[index][0] != cases[index][1]) {
#pragma omp critical
                    {
                        if (result.merged_standard_collisions == 0U) {
                            ++unequal_standard_orders[
                                result.merged_standard_order
                            ];
                        } else {
                            unequal_standard_failures.push_back(cases[index]);
                        }
                    }
                }
                if (!result.commuting_rank
                    || result.global_kernel != result.iota_kernel) {
                    bool expected = true;
                    if (pass.compare_exchange_strong(expected, false)) {
#pragma omp critical
                        {
                            failure_labels = cases[index];
                            failure_result = result;
                        }
                    }
                }
            }
            if (!pass.load()) {
                std::cout << "TAU_EVEN_SWEEP_FAIL labels=";
                print_degree(failure_labels);
                std::cout << " domain=" << failure_result.domain_dimension
                          << " global_kernel=" << failure_result.global_kernel
                          << " iota_kernel=" << failure_result.iota_kernel
                          << " merged_raw_collisions="
                          << failure_result.merged_raw_collisions
                          << " merged_standard_collisions="
                          << failure_result.merged_standard_collisions
                          << " divisor_separated="
                          << failure_result.divisor_separated
                          << " sigma_valid=" << failure_result.sigma_valid
                          << " sigma_plus_kernel="
                          << failure_result.sigma_global_kernel[0]
                          << " sigma_minus_kernel="
                          << failure_result.sigma_global_kernel[1]
                          << " commuting_rank="
                          << failure_result.commuting_rank << '\n';
                return EXIT_SUCCESS;
            }
            std::sort(
                kernel_cases.begin(), kernel_cases.end(),
                [](const auto& first, const auto& second) {
                    return first.first < second.first;
                }
            );
            std::sort(unseparated_cases.begin(), unseparated_cases.end());
            for (const auto& [labels, result] : kernel_cases) {
                std::cout << "TAU_EVEN_KERNEL labels=";
                print_degree(labels);
                std::cout << " domain=" << result.domain_dimension
                          << " kernel=" << result.global_kernel
                          << " merged_raw_collisions="
                          << result.merged_raw_collisions
                          << " merged_standard_collisions="
                          << result.merged_standard_collisions
                          << " divisor_separated="
                          << result.divisor_separated
                          << " sigma_plus_kernel="
                          << result.sigma_global_kernel[0]
                          << " sigma_minus_kernel="
                          << result.sigma_global_kernel[1] << '\n';
            }
            if (std::getenv("PLUCKER_VERBOSE") != nullptr) {
                for (const std::vector<int>& labels : unseparated_cases) {
                    std::cout << "TAU_EVEN_UNSEPARATED labels=";
                    print_degree(labels);
                    std::cout << '\n';
                }
            }
            std::cout << "TAU_EVEN_UNEQUAL_STANDARD_ORDERS";
            for (const auto& [order, count] : unequal_standard_orders) {
                std::cout << " [";
                for (std::size_t position = 0U; position < order.size();
                     ++position) {
                    if (position != 0U) {
                        std::cout << ',';
                    }
                    std::cout << order[position];
                }
                std::cout << "]:" << count;
            }
            std::cout << " failures=" << unequal_standard_failures.size()
                      << '\n';
            std::cout << "TAU_EVEN_SWEEP_PASS maximum=" << maximum
                      << " cases=" << checked.load()
                      << " kernel_cases=" << kernel_cases.size()
                      << " unseparated=" << unseparated.load()
                      << " unseparated_unequal="
                      << unseparated_unequal.load()
                      << " threads=" << omp_get_max_threads() << '\n';
            return EXIT_SUCCESS;
        }
        if (argc >= 7
            && std::string(argv[1]) == "--tau-even-background") {
            std::vector<int> labels;
            labels.reserve(static_cast<std::size_t>(argc - 1));
            for (int argument = 2; argument < argc; ++argument) {
                const long parsed = std::strtol(argv[argument], nullptr, 10);
                if (parsed < 1
                    || parsed > std::numeric_limits<int>::max()) {
                    throw std::runtime_error(
                        "invalid tau-even-background label"
                    );
                }
                if (argument == 4) {
                    labels.push_back(1);
                    labels.push_back(1);
                }
                labels.push_back(static_cast<int>(parsed));
            }
            if (labels.size()
                    >= std::numeric_limits<unsigned int>::digits
                || labels.size() < 7U || (labels[0] & 1) == 0
                || (labels[1] & 1) == 0 || labels[0] < 3
                || labels[1] < 3) {
                throw std::runtime_error(
                    "invalid tau-even-background terminal labels"
                );
            }
            for (std::size_t even = 4U; even < labels.size(); ++even) {
                if ((labels[even] & 1) != 0) {
                    throw std::runtime_error(
                        "tau-even-background labels must be positive even"
                    );
                }
            }
            const TauEvenResult result = analyze_tau_even_allocation(
                labels, false
            );
            const FullBackgroundAllocationResult full_allocation
                = analyze_full_background_allocation(labels);
            const bool naive_kernel_equality
                = result.global_kernel == result.iota_kernel;
            const bool total_payment_capacity = result.global_kernel
                <= full_allocation.iota_kernel;
            const bool sector_payment_capacity = labels[0] != labels[1]
                || (result.sigma_valid && full_allocation.sigma_valid
                    && result.sigma_global_kernel[0]
                        <= full_allocation.sigma_iota_kernel[0]
                    && result.sigma_global_kernel[1]
                        <= full_allocation.sigma_iota_kernel[1]);
            const bool weighted_sector_payment = labels[0] != labels[1]
                || (full_allocation.sigma_valid
                    && full_allocation.sigma_weighted_payment_rank[0]
                    && full_allocation.sigma_weighted_payment_rank[1]);
            const bool pass = result.commuting_rank
                && full_allocation.weighted_payment_rank
                && weighted_sector_payment;
            std::cout << "TAU_EVEN_BACKGROUND_WEIGHTED_PAYMENT_"
                      << (pass ? "PASS" : "FAIL") << " labels=";
            print_degree(labels);
            std::cout << " even_factors=" << labels.size() - 4U
                      << " domain=" << result.domain_dimension
                      << " global_kernel=" << result.global_kernel
                      << " iota_kernel=" << result.iota_kernel
                      << " carrier_kernel=" << result.carrier_kernel
                      << " af_pair="
                      << result.fundamental_pair_reservoir
                      << " full_allocation_domain="
                      << full_allocation.domain_dimension
                      << " full_allocation_kernel="
                      << full_allocation.iota_kernel
                      << " full_global_kernel="
                      << full_allocation.global_kernel
                      << " full_carrier_kernel="
                      << full_allocation.carrier_kernel
                      << " weighted_stacked_kernel="
                      << full_allocation.weighted_stacked_kernel
                      << " weighted_power="
                      << full_allocation.weighted_power
                      << " weighted_family="
                      << (full_allocation.weighted_hash ? "hash" : "power")
                      << " weighted_payment_rank="
                      << full_allocation.weighted_payment_rank
                      << " naive_kernel_equality="
                      << naive_kernel_equality
                      << " one_sided_capacity="
                      << total_payment_capacity << ','
                      << sector_payment_capacity
                      << " sigma_valid=" << result.sigma_valid
                      << " sigma_plus=" << result.sigma_global_kernel[0]
                      << '/' << result.sigma_iota_kernel[0]
                      << " sigma_minus=" << result.sigma_global_kernel[1]
                      << '/' << result.sigma_iota_kernel[1]
                      << " full_sigma_allocation="
                      << full_allocation.sigma_iota_kernel[0] << ','
                      << full_allocation.sigma_iota_kernel[1]
                      << " full_sigma_weighted_stacked="
                      << full_allocation.sigma_weighted_stacked_kernel[0]
                      << ','
                      << full_allocation.sigma_weighted_stacked_kernel[1]
                      << " commuting_rank=" << result.commuting_rank
                      << '\n';
            std::cout << "TAU_EVEN_BACKGROUND_WEIGHTED_TRIAL_KERNELS";
            for (std::size_t trial = 0U;
                 trial < full_allocation
                             .weighted_trial_stacked_kernel.size();
                 ++trial) {
                std::cout << ' ' << trial + 1U << ':'
                          << full_allocation
                                 .weighted_trial_stacked_kernel[trial];
            }
            std::cout << '\n';
            std::cout << "TAU_EVEN_BACKGROUND_HASHED_TRIAL_KERNELS";
            for (std::size_t trial = 0U;
                 trial < full_allocation.hashed_trial_stacked_kernel.size();
                 ++trial) {
                std::cout << ' ' << trial + 1U << ':'
                          << full_allocation
                                 .hashed_trial_stacked_kernel[trial];
            }
            std::cout << '\n';
            return pass ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        if (argc == 4
            && std::string(argv[1]) == "--tau-even-background-sweep") {
            const long parsed_maximum = std::strtol(argv[2], nullptr, 10);
            const long parsed_factors = std::strtol(argv[3], nullptr, 10);
            const long maximum_factors
                = static_cast<long>(
                    std::numeric_limits<unsigned int>::digits
                ) - 5L;
            if (parsed_maximum < 3
                || parsed_maximum > std::numeric_limits<int>::max()
                || parsed_factors < 3 || parsed_factors > maximum_factors) {
                throw std::runtime_error(
                    "invalid tau-even-background-sweep range"
                );
            }
            const int maximum = static_cast<int>(parsed_maximum);
            const std::size_t even_factors
                = static_cast<std::size_t>(parsed_factors);
            std::vector<std::vector<int>> backgrounds;
            std::vector<int> background(even_factors, 2);
            enumerate_even_backgrounds_rec(
                0U, 2, maximum, background, backgrounds
            );
            std::vector<std::vector<int>> cases;
            for (int q = 3; q <= maximum; q += 2) {
                for (int r = q; r <= maximum; r += 2) {
                    for (const std::vector<int>& evens : backgrounds) {
                        std::vector<int> labels{q, r, 1, 1};
                        labels.insert(labels.end(), evens.begin(), evens.end());
                        cases.push_back(std::move(labels));
                    }
                }
            }
            std::atomic<bool> pass{true};
            std::atomic<std::size_t> checked{0U};
            std::atomic<std::size_t> nonempty{0U};
            std::atomic<std::size_t> kernel_cases{0U};
            std::atomic<std::size_t> naive_equality_failures{0U};
            std::atomic<std::size_t> maximum_domain{0U};
            std::atomic<std::size_t> maximum_kernel{0U};
            std::vector<int> failure_labels;
            TauEvenResult failure_result;
            FullBackgroundAllocationResult failure_full_allocation;
#pragma omp parallel for schedule(dynamic)
            for (std::size_t index = 0U; index < cases.size(); ++index) {
                if (!pass.load(std::memory_order_relaxed)) {
                    continue;
                }
                const TauEvenResult result = analyze_tau_even_allocation(
                    cases[index], false
                );
                const FullBackgroundAllocationResult full_allocation
                    = analyze_full_background_allocation(cases[index]);
                checked.fetch_add(1U, std::memory_order_relaxed);
                if (result.domain_dimension != 0U) {
                    nonempty.fetch_add(1U, std::memory_order_relaxed);
                }
                if (result.global_kernel != 0U) {
                    kernel_cases.fetch_add(1U, std::memory_order_relaxed);
                }
                if (result.global_kernel != result.iota_kernel) {
                    naive_equality_failures.fetch_add(
                        1U, std::memory_order_relaxed
                    );
                }
                std::size_t observed = maximum_domain.load(
                    std::memory_order_relaxed
                );
                while (observed < result.domain_dimension
                       && !maximum_domain.compare_exchange_weak(
                           observed, result.domain_dimension,
                           std::memory_order_relaxed
                       )) {
                }
                observed = maximum_kernel.load(std::memory_order_relaxed);
                while (observed < result.global_kernel
                       && !maximum_kernel.compare_exchange_weak(
                           observed, result.global_kernel,
                           std::memory_order_relaxed
                       )) {
                }
                const bool weighted_sector_payment
                    = cases[index][0] != cases[index][1]
                    || (full_allocation.sigma_valid
                        && full_allocation.sigma_weighted_payment_rank[0]
                        && full_allocation.sigma_weighted_payment_rank[1]);
                if (!result.commuting_rank
                    || !full_allocation.weighted_payment_rank
                    || !weighted_sector_payment) {
                    bool expected = true;
                    if (pass.compare_exchange_strong(expected, false)) {
#pragma omp critical
                        {
                            failure_labels = cases[index];
                            failure_result = result;
                            failure_full_allocation = full_allocation;
                        }
                    }
                }
            }
            if (!pass.load()) {
                std::cout << "TAU_EVEN_BACKGROUND_SWEEP_FAIL labels=";
                print_degree(failure_labels);
                std::cout << " domain=" << failure_result.domain_dimension
                          << " global_kernel="
                          << failure_result.global_kernel
                          << " iota_kernel=" << failure_result.iota_kernel
                          << " carrier_kernel="
                          << failure_result.carrier_kernel
                          << " af_pair="
                          << failure_result.fundamental_pair_reservoir
                          << " full_allocation_kernel="
                          << failure_full_allocation.iota_kernel
                          << " full_global_kernel="
                          << failure_full_allocation.global_kernel
                          << " full_carrier_kernel="
                          << failure_full_allocation.carrier_kernel
                          << " weighted_stacked_kernel="
                          << failure_full_allocation.weighted_stacked_kernel
                          << " weighted_power="
                          << failure_full_allocation.weighted_power
                          << " weighted_family="
                          << (failure_full_allocation.weighted_hash
                                  ? "hash" : "power")
                          << " weighted_payment_rank="
                          << failure_full_allocation.weighted_payment_rank
                          << " sigma_valid=" << failure_result.sigma_valid
                          << " sigma_plus="
                          << failure_result.sigma_global_kernel[0] << '/'
                          << failure_result.sigma_iota_kernel[0]
                          << " sigma_minus="
                          << failure_result.sigma_global_kernel[1] << '/'
                          << failure_result.sigma_iota_kernel[1]
                          << " full_sigma_allocation="
                          << failure_full_allocation.sigma_iota_kernel[0]
                          << ','
                          << failure_full_allocation.sigma_iota_kernel[1]
                          << " full_sigma_weighted_stacked="
                          << failure_full_allocation
                                 .sigma_weighted_stacked_kernel[0]
                          << ','
                          << failure_full_allocation
                                 .sigma_weighted_stacked_kernel[1]
                          << " commuting_rank="
                          << failure_result.commuting_rank
                          << " checked=" << checked.load() << '\n';
                return EXIT_FAILURE;
            }
            std::cout
                << "TAU_EVEN_BACKGROUND_WEIGHTED_PAYMENT_SWEEP_PASS maximum="
                      << maximum << " even_factors=" << even_factors
                      << " cases=" << checked.load()
                      << " nonempty=" << nonempty.load()
                      << " kernel_cases=" << kernel_cases.load()
                      << " naive_equality_failures="
                      << naive_equality_failures.load()
                      << " maximum_domain=" << maximum_domain.load()
                      << " maximum_kernel=" << maximum_kernel.load()
                      << " threads=" << omp_get_max_threads() << '\n';
            return EXIT_SUCCESS;
        }
        if (argc == 7 && std::string(argv[1]) == "--tau-even") {
            std::vector<int> labels{0, 0, 1, 1, 0, 0, 0};
            const std::array<std::size_t, 5> positions{0U, 1U, 4U, 5U, 6U};
            for (std::size_t argument = 2U; argument < 7U; ++argument) {
                const long parsed = std::strtol(
                    argv[static_cast<int>(argument)], nullptr, 10
                );
                if (parsed < 1 || parsed > std::numeric_limits<int>::max()) {
                    throw std::runtime_error("invalid tau-even label");
                }
                labels[positions[argument - 2U]] = static_cast<int>(parsed);
            }
            const TauEvenResult result = analyze_tau_even_allocation(labels);
            std::cout << "TAU_EVEN labels=";
            print_degree(labels);
            std::cout << " domain=" << result.domain_dimension
                      << " global_kernel=" << result.global_kernel
                      << " iota_kernel=" << result.iota_kernel
                      << " merged_raw_collisions="
                      << result.merged_raw_collisions
                      << " merged_standard_collisions="
                      << result.merged_standard_collisions
                      << " merged_standard_order=[";
            for (std::size_t place = 0U;
                 place < result.merged_standard_order.size(); ++place) {
                if (place != 0U) {
                    std::cout << ',';
                }
                std::cout << result.merged_standard_order[place];
            }
            std::cout << ']'
                      << " divisor_separated=" << result.divisor_separated
                      << " sigma_valid=" << result.sigma_valid
                      << " sigma_plus=" << result.sigma_global_kernel[0]
                      << '/' << result.sigma_iota_kernel[0]
                      << " sigma_minus=" << result.sigma_global_kernel[1]
                      << '/' << result.sigma_iota_kernel[1]
                      << " commuting_rank=" << result.commuting_rank << '\n';
            if (std::getenv("PLUCKER_VERBOSE") != nullptr) {
                const std::vector<Graph> merged_graphs
                    = tau_even_merged_graphs(labels);
                for (std::size_t index = 0U;
                     index < merged_graphs.size(); ++index) {
                    std::cout << "TAU_EVEN_MERGED_GRAPH index=" << index
                              << ' ';
                    print_graph(merged_graphs[index]);
                    std::cout << '\n';
                }
            }
            if (std::getenv("PLUCKER_EXTREMA") != nullptr) {
                const std::vector<Graph> merged_graphs
                    = tau_even_merged_graphs(labels);
                const bool cubic_exception = labels
                    == std::vector<int>{3, 5, 1, 1, 4, 4, 4};
                const bool double_lower = labels[4] == labels[0] - 1
                    && labels[5] == labels[0] - 1;
                const bool lower_wall_exception = double_lower
                    && labels[6] == labels[0] - 1
                    && labels[1] == 2 * labels[0] - 1;
                const std::array<std::vector<std::size_t>, 3> charts{
                    std::vector<std::size_t>{0U, 3U, 4U, 1U, 5U, 2U},
                    std::vector<std::size_t>{0U, 1U, 5U, 3U, 2U, 4U},
                    std::vector<std::size_t>{0U, 5U, 1U, 3U, 4U, 2U}
                };
                const std::size_t chart_index = cubic_exception
                    ? 2U
                    : (double_lower && !lower_wall_exception ? 1U : 0U);
                std::map<Key, Support> support_memo;
                for (std::size_t index = 0U;
                     index < merged_graphs.size(); ++index) {
                    const Graph permuted = permute_graph(
                        merged_graphs[index], charts[chart_index]
                    );
                    const Support& support = straighten_support(
                        permuted, support_memo
                    );
                    std::cout << "TAU_EVEN_EXTREMA index=" << index
                              << " chart=" << chart_index << " least=";
                    print_graph(unflatten(*support.begin(), 6U));
                    std::cout << " greatest=";
                    print_graph(unflatten(*support.rbegin(), 6U));
                    std::cout << '\n';
                }
            }
            return EXIT_SUCCESS;
        }
        if (argc == 3 && std::string(argv[1]) == "--seven-lift-sweep") {
            const long parsed_maximum = std::strtol(argv[2], nullptr, 10);
            if (parsed_maximum < 3
                || parsed_maximum > std::numeric_limits<int>::max()) {
                throw std::runtime_error("invalid seven-lift-sweep maximum");
            }
            const int maximum = static_cast<int>(parsed_maximum);
            std::vector<std::vector<int>> cases;
            for (int q = 3; q <= maximum; q += 2) {
                for (int r = q; r <= maximum; r += 2) {
                    for (int a = 1; a <= maximum; a += 2) {
                        if (a != 1 && (a == q || a == r)) {
                            continue;
                        }
                        for (int first = 2; first <= maximum; first += 2) {
                            for (int second = first; second <= maximum;
                                 second += 2) {
                                for (int third = second; third <= maximum;
                                     third += 2) {
                                    cases.push_back({
                                        q, r, a, 1, first, second, third
                                    });
                                }
                            }
                        }
                    }
                }
            }
            std::atomic<bool> pass{true};
            std::atomic<std::size_t> checked{0U};
            std::atomic<std::size_t> nonzero_carrier{0U};
            std::atomic<std::size_t> maximum_carrier{0U};
            std::vector<int> failure_labels;
            SevenLiftResult failure_result;
#pragma omp parallel for schedule(dynamic)
            for (std::size_t index = 0U; index < cases.size(); ++index) {
                if (!pass.load(std::memory_order_relaxed)) {
                    continue;
                }
                const SevenLiftResult result = analyze_seven_factor_lift(
                    cases[index], false
                );
                checked.fetch_add(1U, std::memory_order_relaxed);
                if (result.carrier_kernel != 0U) {
                    nonzero_carrier.fetch_add(1U, std::memory_order_relaxed);
                    if (std::getenv("PLUCKER_VERBOSE") != nullptr) {
#pragma omp critical
                        {
                            std::cout << "carrier_case labels=[";
                            for (std::size_t i = 0U; i < cases[index].size(); ++i) {
                                if (i != 0U) {
                                    std::cout << ',';
                                }
                                std::cout << cases[index][i];
                            }
                            std::cout << "] kernel=" << result.source_kernel
                                      << " allocation=" << result.iota_kernel
                                      << " carrier=" << result.carrier_kernel
                                      << " af_pair="
                                      << result.fundamental_pair_reservoir
                                      << " q_stacked_kernel="
                                      << result.q_stacked_kernel_dimension
                                      << '\n';
                        }
                    }
                }
                if (result.iota_kernel != 0U
                    && std::getenv("PLUCKER_VERBOSE") != nullptr) {
#pragma omp critical
                    {
                        std::cout << "allocation_case labels=[";
                        for (std::size_t i = 0U; i < cases[index].size(); ++i) {
                            if (i != 0U) {
                                std::cout << ',';
                            }
                            std::cout << cases[index][i];
                        }
                        std::cout << "] allocation=" << result.iota_kernel
                                  << " q_stacked_kernel="
                                  << result.q_stacked_kernel_dimension
                                  << '\n';
                    }
                }
                std::size_t observed = maximum_carrier.load(
                    std::memory_order_relaxed
                );
                while (observed < result.carrier_kernel
                       && !maximum_carrier.compare_exchange_weak(
                           observed, result.carrier_kernel,
                           std::memory_order_relaxed
                       )) {
                }
                const bool bounded = result.fundamental_pair_reservoir >= 0
                    && result.carrier_kernel
                        <= static_cast<std::size_t>(
                            result.fundamental_pair_reservoir
                        );
                const bool residue_detects = result.split_residue_rank
                    >= result.carrier_kernel;
                const bool linear_residue_detects = result.linear_residue_rank
                    >= result.carrier_kernel;
                const bool minus_residue_detects = result.minus_residue_rank
                    >= result.carrier_kernel;
                const bool side_injective = result.side_global_defect == 0U;
                const bool q_stacked_allocation =
                    result.q_stacked_iota_rank == 0U;
                if (!result.commuting_square || !bounded || !residue_detects
                    || !linear_residue_detects || !minus_residue_detects
                    || !side_injective || !q_stacked_allocation) {
                    bool expected = true;
                    if (pass.compare_exchange_strong(expected, false)) {
#pragma omp critical
                        {
                            failure_labels = cases[index];
                            failure_result = result;
                        }
                    }
                }
            }
            if (!pass.load()) {
                std::cout << "SEVEN_LIFT_SWEEP_FAIL labels=[";
                for (std::size_t i = 0U; i < failure_labels.size(); ++i) {
                    if (i != 0U) {
                        std::cout << ',';
                    }
                    std::cout << failure_labels[i];
                }
                std::cout << "] carrier=" << failure_result.carrier_kernel
                          << " af_pair="
                          << failure_result.fundamental_pair_reservoir
                          << " weighted_residue="
                          << failure_result.split_residue_rank
                          << " linear_residue="
                          << failure_result.linear_residue_rank
                          << " kernel_iota_rank="
                          << failure_result.linear_kernel_iota_rank
                          << " minus_residue="
                          << failure_result.minus_residue_rank
                          << " side_defect="
                          << failure_result.side_global_defect
                          << " q_stacked_iota="
                          << failure_result.q_stacked_iota_rank
                          << " tau_even_best_raw_collisions="
                          << failure_result.tau_even_best_raw_collisions
                          << " square=" << failure_result.commuting_square
                          << '\n';
                return EXIT_SUCCESS;
            }
            std::cout << "SEVEN_LIFT_SWEEP_PASS maximum=" << maximum
                      << " cases=" << checked.load()
                      << " nonzero_carrier=" << nonzero_carrier.load()
                      << " maximum_carrier=" << maximum_carrier.load()
                      << " threads=" << omp_get_max_threads() << '\n';
            return EXIT_SUCCESS;
        }
        if (argc == 8 && std::string(argv[1]) == "--seven-lift") {
            std::vector<int> labels;
            labels.reserve(7U);
            for (int argument = 2; argument < 8; ++argument) {
                const long parsed = std::strtol(argv[argument], nullptr, 10);
                if (parsed < 1 || parsed > std::numeric_limits<int>::max()) {
                    throw std::runtime_error("invalid seven-lift label");
                }
                labels.push_back(static_cast<int>(parsed));
                if (argument == 4) {
                    labels.push_back(1);
                }
            }
            if ((labels[0] & 1) == 0 || (labels[1] & 1) == 0
                || (labels[2] & 1) == 0
                || (labels[4] & 1) != 0 || (labels[5] & 1) != 0
                || (labels[6] & 1) != 0) {
                throw std::runtime_error("seven-lift parity mismatch");
            }
            (void)analyze_seven_factor_lift(labels, true);
            return EXIT_SUCCESS;
        }
        if (argc != 4) {
            throw std::runtime_error(
                "usage: verify_su2_plucker_second_syzygies"
                " VERTICES MAXIMUM_COORDINATE MAXIMUM_TOTAL"
                " | --seven-lift Q R A E1 E2 E3"
                " | --seven-lift-sweep MAXIMUM_LABEL"
                " | --tau-even Q R E1 E2 E3"
                " | --tau-even-sweep MAXIMUM_LABEL"
                " | --tau-even-background Q R E1 E2 E3 [E4 ...]"
                " | --tau-even-background-sweep MAXIMUM_LABEL EVEN_FACTORS"
            );
        }
        const long parsed_vertices = std::strtol(argv[1], nullptr, 10);
        const long parsed_coordinate = std::strtol(argv[2], nullptr, 10);
        const long parsed_total = std::strtol(argv[3], nullptr, 10);
        if (parsed_vertices < 3 || parsed_vertices > 9
            || parsed_coordinate < 1
            || parsed_coordinate > std::numeric_limits<int>::max()
            || parsed_total < 3
            || parsed_total > std::numeric_limits<int>::max()) {
            throw std::runtime_error("invalid argument range");
        }
        const std::size_t vertices = static_cast<std::size_t>(parsed_vertices);
        const int maximum_coordinate = static_cast<int>(parsed_coordinate);
        const int maximum_total = static_cast<int>(parsed_total);
        std::vector<std::vector<int>> degrees;
        std::vector<int> degree(vertices, 0);
        enumerate_degrees_rec(
            0U, maximum_total, maximum_coordinate, degree, degrees
        );

        std::atomic<bool> pass{true};
        std::atomic<std::size_t> checked{0U};
        std::atomic<std::size_t> nonzero_kernels{0U};
        std::atomic<std::size_t> maximum_kernel{0U};
        std::vector<int> failure_degree;
        DegreeResult failure_result;
#pragma omp parallel for schedule(dynamic)
        for (std::size_t index = 0U; index < degrees.size(); ++index) {
            if (!pass.load(std::memory_order_relaxed)) {
                continue;
            }
            const DegreeResult result = analyze_degree(degrees[index]);
            checked.fetch_add(1U, std::memory_order_relaxed);
            if (result.kernel_dimension != 0U) {
                nonzero_kernels.fetch_add(1U, std::memory_order_relaxed);
            }
            std::size_t observed = maximum_kernel.load(std::memory_order_relaxed);
            while (observed < result.kernel_dimension
                   && !maximum_kernel.compare_exchange_weak(
                       observed, result.kernel_dimension,
                       std::memory_order_relaxed
                   )) {
            }
            if (!result.pass) {
                bool expected = true;
                if (pass.compare_exchange_strong(expected, false)) {
#pragma omp critical
                    {
                        failure_degree = degrees[index];
                        failure_result = result;
                    }
                }
            }
        }
        if (!pass.load()) {
            std::cout << "PLUCKER_SECOND_SYZYGY_FAIL degree=";
            print_degree(failure_degree);
            std::cout << " f1=" << failure_result.f1_dimension
                      << " kernel=" << failure_result.kernel_dimension
                      << " generators=" << failure_result.generator_count
                      << " generator_rank=" << failure_result.generator_rank
                      << '\n';
            return EXIT_SUCCESS;
        }
        std::cout << "PLUCKER_SECOND_SYZYGY_PASS vertices=" << vertices
                  << " maximum_coordinate=" << maximum_coordinate
                  << " maximum_total=" << maximum_total
                  << " checked=" << checked.load()
                  << " nonzero_kernels=" << nonzero_kernels.load()
                  << " maximum_kernel=" << maximum_kernel.load()
                  << " threads=" << omp_get_max_threads() << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
