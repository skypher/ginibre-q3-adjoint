#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
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
    std::vector<std::vector<Rational>> nullspace;
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
    bool diagnostics = true
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
                result.column_labels.emplace_back(
                    subset, odd_columns[subset].size() - 1U
                );
                columns.push_back(std::move(column));
            }
        }
    }

    const std::vector<std::vector<Integer>> full_matrix
        = column_matrix(columns, basis_index.size());
    result.rank = exact_rank(full_matrix);
    if (diagnostics) {
        result.nullspace = exact_nullspace(full_matrix);
        for (auto first = odd_columns.begin(); first != odd_columns.end(); ++first) {
            for (auto second = std::next(first);
                 second != odd_columns.end(); ++second) {
                std::vector<std::vector<Integer>> combined = first->second;
                combined.insert(
                    combined.end(), second->second.begin(), second->second.end()
                );
                const std::size_t union_rank = exact_rank(
                    column_matrix(combined, basis_index.size())
                );
                const std::size_t intersection
                    = first->second.size() + second->second.size() - union_rank;
                result.pair_intersections.push_back(Result::PairIntersection{
                    first->first,
                    second->first,
                    intersection,
                    channel_dimensions[first->first ^ second->first]
                });
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
        if (argc == 6 && std::string(argv[1]) == "sweep") {
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
                            const Result result = analyze_case(
                                labels, target, minus_mask, false
                            );
                            ++cases;
                            if (result.odd_sources
                                <= result.positive_sources) {
                                continue;
                            }
                            ++deficit_cases;
                            const std::uint64_t kernel
                                = result.odd_sources
                                    - static_cast<std::uint64_t>(result.rank);
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
                  << " result="
                  << (kernel <= result.positive_sources ? "PASS" : "FAIL")
                  << '\n';
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
        return kernel <= result.positive_sources ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
