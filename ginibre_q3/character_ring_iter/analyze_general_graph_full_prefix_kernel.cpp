#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Matrix = std::vector<std::vector<int>>;

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(value);
}

Integer binomial(int n, int r) {
    if (r < 0 || r > n) {
        return 0;
    }
    r = std::min(r, n - r);
    Integer result = 1;
    for (int index = 1; index <= r; ++index) {
        result *= n - r + index;
        result /= index;
    }
    return result;
}

std::vector<Integer> multiply(
    const Matrix& adjacency,
    const std::vector<Integer>& state
) {
    const int vertices = static_cast<int>(adjacency.size());
    std::vector<Integer> next(static_cast<std::size_t>(vertices));
    for (int source = 0; source < vertices; ++source) {
        for (int target = 0; target < vertices; ++target) {
            if (adjacency[static_cast<std::size_t>(source)]
                         [static_cast<std::size_t>(target)] != 0) {
                next[static_cast<std::size_t>(target)] +=
                    state[static_cast<std::size_t>(source)];
            }
        }
    }
    return next;
}

bool connected_from_root(const Matrix& adjacency) {
    const int vertices = static_cast<int>(adjacency.size());
    std::vector<bool> reached(static_cast<std::size_t>(vertices), false);
    reached[0] = true;
    bool changed = true;
    while (changed) {
        changed = false;
        for (int source = 0; source < vertices; ++source) {
            if (!reached[static_cast<std::size_t>(source)]) {
                continue;
            }
            for (int target = 0; target < vertices; ++target) {
                if (adjacency[static_cast<std::size_t>(source)]
                             [static_cast<std::size_t>(target)] != 0
                    && !reached[static_cast<std::size_t>(target)]) {
                    reached[static_cast<std::size_t>(target)] = true;
                    changed = true;
                }
            }
        }
    }
    return std::all_of(
        reached.begin(),
        reached.end(),
        [](bool value) { return value; }
    );
}

bool chordal(const Matrix& adjacency) {
    const int vertices = static_cast<int>(adjacency.size());
    std::vector<int> weight(static_cast<std::size_t>(vertices), 0);
    std::vector<int> number(static_cast<std::size_t>(vertices), 0);
    std::vector<bool> selected(
        static_cast<std::size_t>(vertices),
        false
    );
    for (int step = 0; step < vertices; ++step) {
        int choice = -1;
        for (int vertex = 0; vertex < vertices; ++vertex) {
            if (!selected[static_cast<std::size_t>(vertex)]
                && (
                    choice < 0
                    || weight[static_cast<std::size_t>(vertex)]
                        > weight[static_cast<std::size_t>(choice)]
                )) {
                choice = vertex;
            }
        }
        selected[static_cast<std::size_t>(choice)] = true;
        number[static_cast<std::size_t>(choice)] = vertices - step;
        for (int neighbor = 0; neighbor < vertices; ++neighbor) {
            if (!selected[static_cast<std::size_t>(neighbor)]
                && adjacency[static_cast<std::size_t>(choice)]
                            [static_cast<std::size_t>(neighbor)] != 0) {
                ++weight[static_cast<std::size_t>(neighbor)];
            }
        }
    }
    for (int vertex = 0; vertex < vertices; ++vertex) {
        std::vector<int> higher_neighbors;
        for (int neighbor = 0; neighbor < vertices; ++neighbor) {
            if (neighbor != vertex
                && adjacency[static_cast<std::size_t>(vertex)]
                            [static_cast<std::size_t>(neighbor)] != 0
                && number[static_cast<std::size_t>(neighbor)]
                    > number[static_cast<std::size_t>(vertex)]) {
                higher_neighbors.push_back(neighbor);
            }
        }
        for (std::size_t first = 0U;
             first < higher_neighbors.size();
             ++first) {
            for (std::size_t second = first + 1U;
                 second < higher_neighbors.size();
                 ++second) {
                if (adjacency[
                        static_cast<std::size_t>(
                            higher_neighbors[first]
                        )
                    ][static_cast<std::size_t>(
                        higher_neighbors[second]
                    )] == 0) {
                    return false;
                }
            }
        }
    }
    return true;
}

bool outside_in_peo(const Matrix& adjacency) {
    const int vertices = static_cast<int>(adjacency.size());
    std::vector<int> order;
    order.reserve(static_cast<std::size_t>(vertices));
    int lower = 0;
    int upper = vertices - 1;
    while (lower <= upper) {
        order.push_back(lower);
        if (lower != upper) {
            order.push_back(upper);
        }
        ++lower;
        --upper;
    }
    std::vector<bool> removed(
        static_cast<std::size_t>(vertices),
        false
    );
    for (const int vertex : order) {
        std::vector<int> remaining_neighbors;
        for (int neighbor = 0; neighbor < vertices; ++neighbor) {
            if (neighbor != vertex
                && !removed[static_cast<std::size_t>(neighbor)]
                && adjacency[static_cast<std::size_t>(vertex)]
                            [static_cast<std::size_t>(neighbor)] != 0) {
                remaining_neighbors.push_back(neighbor);
            }
        }
        for (std::size_t first = 0U;
             first < remaining_neighbors.size();
             ++first) {
            for (std::size_t second = first + 1U;
                 second < remaining_neighbors.size();
                 ++second) {
                if (adjacency[
                        static_cast<std::size_t>(
                            remaining_neighbors[first]
                        )
                    ][static_cast<std::size_t>(
                        remaining_neighbors[second]
                    )] == 0) {
                    return false;
                }
            }
        }
        removed[static_cast<std::size_t>(vertex)] = true;
    }
    return true;
}

void print_graph(const Matrix& adjacency) {
    std::cout << " adjacency=";
    for (const auto& row : adjacency) {
        std::cout << '[';
        for (const int entry : row) {
            std::cout << entry;
        }
        std::cout << ']';
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 4) {
            throw std::runtime_error(
                "usage: analyze_general_graph_full_prefix_kernel "
                "MAXIMUM_VERTICES MAXIMUM_PREFIX "
                "[--root-leaf|--reflected-leaves|--reflected-roots|"
                "--reflected-convex-leaves|--reflected-chordal-leaves|"
                "--reflected-outside-in-leaves|"
                "--reflected-outside-in-loop-leaves|"
                "--reflected-loop-leaves]"
            );
        }
        const bool root_leaf =
            argc == 4 && std::string(argv[3]) == "--root-leaf";
        const bool reflected_leaves =
            argc == 4 && std::string(argv[3]) == "--reflected-leaves";
        const bool reflected_roots =
            argc == 4 && std::string(argv[3]) == "--reflected-roots";
        const bool reflected_convex_leaves =
            argc == 4
            && std::string(argv[3]) == "--reflected-convex-leaves";
        const bool reflected_chordal_leaves =
            argc == 4
            && std::string(argv[3]) == "--reflected-chordal-leaves";
        const bool reflected_outside_in_leaves =
            argc == 4
            && std::string(argv[3]) == "--reflected-outside-in-leaves";
        const bool reflected_outside_in_loop_leaves =
            argc == 4
            && std::string(argv[3])
                == "--reflected-outside-in-loop-leaves";
        const bool reflected_loop_leaves =
            argc == 4
            && std::string(argv[3]) == "--reflected-loop-leaves";
        if (argc == 4 && !root_leaf && !reflected_leaves
            && !reflected_roots && !reflected_convex_leaves
            && !reflected_chordal_leaves
            && !reflected_outside_in_leaves
            && !reflected_outside_in_loop_leaves
            && !reflected_loop_leaves) {
            throw std::runtime_error(
                "the optional argument is --root-leaf or "
                "--reflected-leaves, --reflected-roots, or "
                "--reflected-convex-leaves, or "
                "--reflected-chordal-leaves, or "
                "--reflected-outside-in-leaves, or "
                "--reflected-outside-in-loop-leaves, or "
                "--reflected-loop-leaves"
            );
        }
        const int maximum_vertices =
            parse_positive(argv[1], "maximum vertices");
        const int maximum_prefix =
            parse_positive(argv[2], "maximum prefix");
        if (
            maximum_vertices
                > (
                    reflected_leaves || reflected_roots
                        || reflected_convex_leaves
                        || reflected_chordal_leaves
                        || reflected_outside_in_leaves
                        || reflected_outside_in_loop_leaves
                        || reflected_loop_leaves
                    ? (
                        reflected_convex_leaves
                            || reflected_outside_in_loop_leaves
                            || reflected_loop_leaves
                        ? 9
                        : 8
                    )
                    : 6
                )
        ) {
            throw std::runtime_error(
                "exhaustive mode supports at most six general vertices "
                "or eight reflected vertices; the reflected convex-leaf "
                "and reflected outside-in-loop-leaf modes support nine"
            );
        }
        if (maximum_prefix < 4) {
            throw std::runtime_error(
                "maximum prefix must be at least four"
            );
        }

        std::uint64_t graphs = 0U;
        std::uint64_t connected_graphs = 0U;
        std::uint64_t coordinates = 0U;
        for (int vertices = 1;
             vertices <= maximum_vertices;
             ++vertices) {
            using Edge = std::pair<int, int>;
            std::vector<std::vector<Edge>> edge_orbits;
            std::map<Edge, int> edge_orbit_ids;
            for (int row = 0; row < vertices; ++row) {
                for (int column = row;
                     column < vertices;
                     ++column) {
                    const Edge edge{row, column};
                    if (edge_orbit_ids.find(edge)
                        != edge_orbit_ids.end()) {
                        continue;
                    }
                    const int orbit_id =
                        static_cast<int>(edge_orbits.size());
                    edge_orbits.push_back({edge});
                    edge_orbit_ids.emplace(edge, orbit_id);
                    if (
                        reflected_leaves || reflected_roots
                            || reflected_convex_leaves
                            || reflected_chordal_leaves
                            || reflected_outside_in_leaves
                            || reflected_outside_in_loop_leaves
                            || reflected_loop_leaves
                    ) {
                        const int reflected_row =
                            vertices - 1 - column;
                        const int reflected_column =
                            vertices - 1 - row;
                        const Edge reflected_edge{
                            reflected_row,
                            reflected_column
                        };
                        if (reflected_edge != edge) {
                            edge_orbits.back().push_back(reflected_edge);
                            edge_orbit_ids.emplace(
                                reflected_edge,
                                orbit_id
                            );
                        }
                    }
                }
            }
            const int bits = static_cast<int>(edge_orbits.size());
            if (bits >= 63) {
                throw std::runtime_error("graph mask exceeds uint64");
            }
            const std::uint64_t masks = std::uint64_t{1} << bits;
            for (std::uint64_t mask = 0U; mask < masks; ++mask) {
                ++graphs;
                if (reflected_leaves || reflected_convex_leaves
                    || reflected_chordal_leaves
                    || reflected_outside_in_leaves
                    || reflected_outside_in_loop_leaves
                    || reflected_loop_leaves) {
                    const int root_loop_orbit =
                        edge_orbit_ids.at(Edge{0, 0});
                    if (
                        (
                            mask
                            & (
                                std::uint64_t{1}
                                << root_loop_orbit
                            )
                        ) != 0U
                    ) {
                        continue;
                    }
                    int root_degree = 0;
                    for (int column = 1;
                         column < vertices;
                         ++column) {
                        const int orbit =
                            edge_orbit_ids.at(Edge{0, column});
                        if (
                            (
                                mask
                                & (std::uint64_t{1} << orbit)
                            ) != 0U
                        ) {
                            ++root_degree;
                        }
                    }
                    if (root_degree != 1) {
                        continue;
                    }
                }
                Matrix adjacency(
                    static_cast<std::size_t>(vertices),
                    std::vector<int>(
                        static_cast<std::size_t>(vertices),
                        0
                    )
                );
                for (int bit = 0; bit < bits; ++bit) {
                    const int entry =
                        (mask & (std::uint64_t{1} << bit)) == 0U
                            ? 0
                            : 1;
                    for (const Edge& edge :
                         edge_orbits[static_cast<std::size_t>(bit)]) {
                        adjacency[static_cast<std::size_t>(edge.first)]
                                 [static_cast<std::size_t>(edge.second)] =
                            entry;
                        adjacency[static_cast<std::size_t>(edge.second)]
                                 [static_cast<std::size_t>(edge.first)] =
                            entry;
                    }
                }
                if (!connected_from_root(adjacency)) {
                    continue;
                }
                if (
                    root_leaf || reflected_leaves
                        || reflected_convex_leaves
                        || reflected_chordal_leaves
                        || reflected_outside_in_leaves
                        || reflected_outside_in_loop_leaves
                        || reflected_loop_leaves
                ) {
                    int root_degree = 0;
                    for (int vertex = 0;
                         vertex < vertices;
                         ++vertex) {
                        root_degree +=
                            adjacency[0][static_cast<std::size_t>(vertex)];
                    }
                    if (adjacency[0][0] != 0 || root_degree != 1) {
                        continue;
                    }
                }
                if (
                    reflected_leaves || reflected_roots
                        || reflected_convex_leaves
                        || reflected_chordal_leaves
                        || reflected_outside_in_leaves
                        || reflected_outside_in_loop_leaves
                        || reflected_loop_leaves
                ) {
                    bool reflected = true;
                    for (int row = 0; row < vertices; ++row) {
                        for (int column = 0;
                             column < vertices;
                             ++column) {
                            if (
                                adjacency[
                                    static_cast<std::size_t>(row)
                                ][static_cast<std::size_t>(column)]
                                != adjacency[
                                    static_cast<std::size_t>(
                                        vertices - 1 - row
                                    )
                                ][static_cast<std::size_t>(
                                    vertices - 1 - column
                                )]
                            ) {
                                reflected = false;
                            }
                        }
                    }
                    if (!reflected) {
                        continue;
                    }
                }
                if (reflected_convex_leaves) {
                    bool row_convex = true;
                    for (int row = 0; row < vertices; ++row) {
                        int first = vertices;
                        int last = -1;
                        for (int column = 0;
                             column < vertices;
                             ++column) {
                            if (
                                adjacency[
                                    static_cast<std::size_t>(row)
                                ][static_cast<std::size_t>(column)] != 0
                            ) {
                                first = std::min(first, column);
                                last = std::max(last, column);
                            }
                        }
                        for (int column = first;
                             column <= last;
                             ++column) {
                            if (
                                adjacency[
                                    static_cast<std::size_t>(row)
                                ][static_cast<std::size_t>(column)] == 0
                            ) {
                                row_convex = false;
                            }
                        }
                    }
                    if (!row_convex) {
                        continue;
                    }
                }
                if (reflected_chordal_leaves && !chordal(adjacency)) {
                    continue;
                }
                if (reflected_outside_in_leaves
                    && !outside_in_peo(adjacency)) {
                    continue;
                }
                if (reflected_outside_in_loop_leaves) {
                    if (!outside_in_peo(adjacency)) {
                        continue;
                    }
                    int root_neighbor = -1;
                    for (int vertex = 1;
                         vertex < vertices;
                         ++vertex) {
                        if (adjacency[0][
                                static_cast<std::size_t>(vertex)
                            ] != 0) {
                            root_neighbor = vertex;
                            break;
                        }
                    }
                    if (root_neighbor < 0
                        || adjacency[
                            static_cast<std::size_t>(root_neighbor)
                        ][static_cast<std::size_t>(root_neighbor)] == 0) {
                        continue;
                    }
                }
                if (reflected_loop_leaves) {
                    int root_neighbor = -1;
                    for (int vertex = 1;
                         vertex < vertices;
                         ++vertex) {
                        if (adjacency[0][
                                static_cast<std::size_t>(vertex)
                            ] != 0) {
                            root_neighbor = vertex;
                            break;
                        }
                    }
                    if (root_neighbor < 0
                        || adjacency[
                            static_cast<std::size_t>(root_neighbor)
                        ][static_cast<std::size_t>(root_neighbor)] == 0) {
                        continue;
                    }
                }
                ++connected_graphs;
                std::vector<std::vector<Integer>> powers(
                    static_cast<std::size_t>(2 * maximum_prefix + 3)
                );
                powers[0].assign(
                    static_cast<std::size_t>(vertices),
                    Integer(0)
                );
                powers[0][0] = 1;
                for (int power = 1;
                     power <= 2 * maximum_prefix + 2;
                     ++power) {
                    powers[static_cast<std::size_t>(power)] =
                        multiply(
                            adjacency,
                            powers[static_cast<std::size_t>(power - 1)]
                        );
                }
                for (int prefix = 4;
                     prefix <= maximum_prefix;
                     ++prefix) {
                    const int total_length = 2 * prefix + 2;
                    std::vector<Integer> kernel(
                        static_cast<std::size_t>(vertices)
                    );
                    for (int truncation = 0;
                         truncation <= prefix - 1;
                         ++truncation) {
                        const int even_power = 2 * truncation;
                        const int odd_power = even_power + 1;
                        const Integer weight =
                            binomial(total_length - 1, even_power);
                        const Integer& even_return =
                            powers[
                                static_cast<std::size_t>(even_power)
                            ][0];
                        const Integer& odd_return =
                            powers[
                                static_cast<std::size_t>(odd_power)
                            ][0];
                        const int first_vertex =
                            reflected_leaves || reflected_roots
                                    || reflected_convex_leaves
                                    || reflected_chordal_leaves
                                    || reflected_outside_in_leaves
                                    || reflected_outside_in_loop_leaves
                                    || reflected_loop_leaves
                                ? vertices - 1
                                : 0;
                        for (int vertex = first_vertex;
                             vertex < vertices;
                             ++vertex) {
                            kernel[static_cast<std::size_t>(vertex)] +=
                                weight
                                * (
                                    even_return
                                        * powers[
                                            static_cast<std::size_t>(
                                                total_length - even_power
                                            )
                                        ][static_cast<std::size_t>(vertex)]
                                    - odd_return
                                        * powers[
                                            static_cast<std::size_t>(
                                                total_length - odd_power
                                            )
                                        ][static_cast<std::size_t>(vertex)]
                                );
                        }
                        if (truncation < 2) {
                            continue;
                        }
                        for (int vertex = 0;
                             vertex < vertices;
                             ++vertex) {
                            ++coordinates;
                            const Integer& value =
                                kernel[static_cast<std::size_t>(vertex)];
                            if (value >= 0) {
                                continue;
                            }
                            std::cout
                                << "GENERAL_GRAPH_FULL_PREFIX_KERNEL"
                                << " counterexample"
                                << " vertices=" << vertices
                                << " mask=" << mask
                                << " prefix=" << prefix
                                << " truncation=" << truncation
                                << " vertex=" << vertex
                                << " value=" << value;
                            print_graph(adjacency);
                            std::cout << " result=FAIL_UNIVERSAL\n";
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }
        std::cout
            << "GENERAL_GRAPH_FULL_PREFIX_KERNEL"
            << " maximum_vertices=" << maximum_vertices
            << " maximum_prefix=" << maximum_prefix
            << " root_leaf=" << (root_leaf ? 1 : 0)
            << " reflected_leaves=" << (reflected_leaves ? 1 : 0)
            << " reflected_roots=" << (reflected_roots ? 1 : 0)
            << " reflected_convex_leaves="
                << (reflected_convex_leaves ? 1 : 0)
            << " reflected_chordal_leaves="
                << (reflected_chordal_leaves ? 1 : 0)
            << " reflected_outside_in_leaves="
                << (reflected_outside_in_leaves ? 1 : 0)
            << " reflected_outside_in_loop_leaves="
                << (reflected_outside_in_loop_leaves ? 1 : 0)
            << " reflected_loop_leaves="
                << (reflected_loop_leaves ? 1 : 0)
            << " graphs=" << graphs
            << " connected_graphs=" << connected_graphs
            << " coordinates=" << coordinates
            << " result=PASS_EXHAUSTIVE_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "GENERAL_GRAPH_FULL_PREFIX_KERNEL FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
