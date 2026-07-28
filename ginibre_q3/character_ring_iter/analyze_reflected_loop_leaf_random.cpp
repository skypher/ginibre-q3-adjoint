#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
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

bool connected(const Matrix& adjacency) {
    const int vertices = static_cast<int>(adjacency.size());
    std::vector<bool> reached(
        static_cast<std::size_t>(vertices),
        false
    );
    reached[0] = true;
    for (int pass = 0; pass < vertices; ++pass) {
        for (int source = 0; source < vertices; ++source) {
            if (!reached[static_cast<std::size_t>(source)]) {
                continue;
            }
            for (int target = 0; target < vertices; ++target) {
                if (adjacency[static_cast<std::size_t>(source)]
                             [static_cast<std::size_t>(target)] != 0) {
                    reached[static_cast<std::size_t>(target)] = true;
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

bool outside_in_peo(const Matrix& adjacency) {
    const int vertices = static_cast<int>(adjacency.size());
    std::vector<bool> removed(
        static_cast<std::size_t>(vertices),
        false
    );
    int lower = 0;
    int upper = vertices - 1;
    while (lower <= upper) {
        for (int side = 0; side < 2; ++side) {
            if (side == 1 && lower == upper) {
                continue;
            }
            const int vertex = side == 0 ? lower : upper;
            std::vector<int> neighbors;
            for (int candidate = 0;
                 candidate < vertices;
                 ++candidate) {
                if (candidate != vertex
                    && !removed[static_cast<std::size_t>(candidate)]
                    && adjacency[static_cast<std::size_t>(vertex)]
                                [static_cast<std::size_t>(candidate)] != 0) {
                    neighbors.push_back(candidate);
                }
            }
            for (std::size_t first = 0U;
                 first < neighbors.size();
                 ++first) {
                for (std::size_t second = first + 1U;
                     second < neighbors.size();
                     ++second) {
                    if (adjacency[
                            static_cast<std::size_t>(neighbors[first])
                        ][static_cast<std::size_t>(neighbors[second])] == 0) {
                        return false;
                    }
                }
            }
            removed[static_cast<std::size_t>(vertex)] = true;
        }
        ++lower;
        --upper;
    }
    return true;
}

void row_convex_closure(Matrix& adjacency) {
    const int vertices = static_cast<int>(adjacency.size());
    bool changed = true;
    while (changed) {
        changed = false;
        for (int row = 0; row < vertices; ++row) {
            int first = vertices;
            int last = -1;
            for (int column = 0; column < vertices; ++column) {
                if (adjacency[static_cast<std::size_t>(row)]
                             [static_cast<std::size_t>(column)] != 0) {
                    first = std::min(first, column);
                    last = std::max(last, column);
                }
            }
            for (int column = first; column <= last; ++column) {
                const int reflected_row = vertices - 1 - column;
                const int reflected_column = vertices - 1 - row;
                if (adjacency[static_cast<std::size_t>(row)]
                             [static_cast<std::size_t>(column)] == 0) {
                    changed = true;
                }
                adjacency[static_cast<std::size_t>(row)]
                         [static_cast<std::size_t>(column)] = 1;
                adjacency[static_cast<std::size_t>(column)]
                         [static_cast<std::size_t>(row)] = 1;
                adjacency[static_cast<std::size_t>(reflected_row)]
                         [static_cast<std::size_t>(reflected_column)] = 1;
                adjacency[static_cast<std::size_t>(reflected_column)]
                         [static_cast<std::size_t>(reflected_row)] = 1;
            }
        }
    }
}

bool row_convex(const Matrix& adjacency) {
    const int vertices = static_cast<int>(adjacency.size());
    for (int row = 0; row < vertices; ++row) {
        int first = vertices;
        int last = -1;
        for (int column = 0; column < vertices; ++column) {
            if (adjacency[static_cast<std::size_t>(row)]
                         [static_cast<std::size_t>(column)] != 0) {
                first = std::min(first, column);
                last = std::max(last, column);
            }
        }
        for (int column = first; column <= last; ++column) {
            if (adjacency[static_cast<std::size_t>(row)]
                         [static_cast<std::size_t>(column)] == 0) {
                return false;
            }
        }
    }
    return true;
}

std::vector<Integer> multiply(
    const Matrix& adjacency,
    const std::vector<Integer>& state
) {
    const int vertices = static_cast<int>(adjacency.size());
    std::vector<Integer> next(static_cast<std::size_t>(vertices));
    for (int source = 0; source < vertices; ++source) {
        if (state[static_cast<std::size_t>(source)] == 0) {
            continue;
        }
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
        if (argc != 4 && argc != 5) {
            throw std::runtime_error(
                "usage: analyze_reflected_loop_leaf_random "
                "MAXIMUM_VERTICES MAXIMUM_PREFIX SAMPLES_PER_SIZE "
                "[--outside-in|--constructed-outside-in|"
                "--constructed-row-convex]"
            );
        }
        const bool require_outside_in =
            argc == 5 && std::string(argv[4]) == "--outside-in";
        const bool construct_outside_in =
            argc == 5
            && std::string(argv[4]) == "--constructed-outside-in";
        const bool construct_row_convex =
            argc == 5
            && std::string(argv[4]) == "--constructed-row-convex";
        if (argc == 5 && !require_outside_in
            && !construct_outside_in && !construct_row_convex) {
            throw std::runtime_error(
                "the optional argument is --outside-in or "
                "--constructed-outside-in or --constructed-row-convex"
            );
        }
        const int maximum_vertices =
            parse_positive(argv[1], "maximum vertices");
        const int maximum_prefix =
            parse_positive(argv[2], "maximum prefix");
        const int samples_per_size =
            parse_positive(argv[3], "samples per size");
        if (maximum_vertices < 6 || maximum_prefix < 4) {
            throw std::runtime_error(
                "maximum vertices must be at least six and "
                "maximum prefix at least four"
            );
        }

        std::mt19937_64 generator(0x5a17c0deULL);
        std::uint64_t attempts = 0U;
        std::uint64_t connected_graphs = 0U;
        std::uint64_t prefix_rows = 0U;

        for (int vertices = 6;
             vertices <= maximum_vertices;
             ++vertices) {
            for (int sample = 0; sample < samples_per_size; ++sample) {
                ++attempts;
                Matrix adjacency(
                    static_cast<std::size_t>(vertices),
                    std::vector<int>(
                        static_cast<std::size_t>(vertices),
                        0
                    )
                );
                const int maximum_neighbor = (vertices - 2) / 2;
                const int root_neighbor =
                    1 + static_cast<int>(
                        generator()
                        % static_cast<std::uint64_t>(maximum_neighbor)
                    );
                const int reflected_neighbor =
                    vertices - 1 - root_neighbor;
                const int density_tenths =
                    1 + static_cast<int>(generator() % 9U);
                if (construct_outside_in) {
                    const int inner_lower = (vertices - 1) / 2;
                    const int inner_upper = vertices / 2;
                    adjacency[
                        static_cast<std::size_t>(inner_lower)
                    ][static_cast<std::size_t>(inner_lower)] =
                        static_cast<int>(generator() & 1U);
                    adjacency[
                        static_cast<std::size_t>(inner_upper)
                    ][static_cast<std::size_t>(inner_upper)] =
                        adjacency[
                            static_cast<std::size_t>(inner_lower)
                        ][static_cast<std::size_t>(inner_lower)];
                    if (inner_lower != inner_upper) {
                        adjacency[
                            static_cast<std::size_t>(inner_lower)
                        ][static_cast<std::size_t>(inner_upper)] = 1;
                        adjacency[
                            static_cast<std::size_t>(inner_upper)
                        ][static_cast<std::size_t>(inner_lower)] = 1;
                    }
                    for (int row = inner_lower - 1;
                         row >= 1;
                         --row) {
                        const int reflected_row = vertices - 1 - row;
                        const int first_inner = row + 1;
                        const int last_inner = reflected_row - 1;
                        const int anchor =
                            first_inner
                            + static_cast<int>(
                                generator()
                                % static_cast<std::uint64_t>(
                                    last_inner - first_inner + 1
                                )
                            );
                        std::vector<int> clique{anchor};
                        for (int candidate = first_inner;
                             candidate <= last_inner;
                             ++candidate) {
                            if (candidate == anchor
                                || static_cast<int>(generator() % 10U)
                                    >= density_tenths) {
                                continue;
                            }
                            bool joins_clique = true;
                            for (const int member : clique) {
                                if (adjacency[
                                        static_cast<std::size_t>(
                                            candidate
                                        )
                                    ][static_cast<std::size_t>(
                                        member
                                    )] == 0) {
                                    joins_clique = false;
                                }
                            }
                            if (joins_clique) {
                                clique.push_back(candidate);
                            }
                        }
                        for (const int neighbor : clique) {
                            const int reflected_neighbor_of_row =
                                vertices - 1 - neighbor;
                            adjacency[
                                static_cast<std::size_t>(row)
                            ][static_cast<std::size_t>(neighbor)] = 1;
                            adjacency[
                                static_cast<std::size_t>(neighbor)
                            ][static_cast<std::size_t>(row)] = 1;
                            adjacency[
                                static_cast<std::size_t>(reflected_row)
                            ][static_cast<std::size_t>(
                                reflected_neighbor_of_row
                            )] = 1;
                            adjacency[
                                static_cast<std::size_t>(
                                    reflected_neighbor_of_row
                                )
                            ][static_cast<std::size_t>(reflected_row)] = 1;
                        }
                        const int loop =
                            static_cast<int>(generator() & 1U);
                        adjacency[static_cast<std::size_t>(row)]
                                 [static_cast<std::size_t>(row)] = loop;
                        adjacency[
                            static_cast<std::size_t>(reflected_row)
                        ][static_cast<std::size_t>(reflected_row)] = loop;
                    }
                } else {
                    for (int row = 1; row < vertices - 1; ++row) {
                        for (int column = row;
                             column < vertices - 1;
                             ++column) {
                            const int reflected_row =
                                vertices - 1 - column;
                            const int reflected_column =
                                vertices - 1 - row;
                            if (row > reflected_row
                                || (
                                    row == reflected_row
                                    && column > reflected_column
                                )) {
                                continue;
                            }
                            const int entry =
                                static_cast<int>(generator() % 10U)
                                    < density_tenths
                                ? 1
                                : 0;
                            adjacency[static_cast<std::size_t>(row)]
                                     [static_cast<std::size_t>(column)] =
                                entry;
                            adjacency[static_cast<std::size_t>(column)]
                                     [static_cast<std::size_t>(row)] =
                                entry;
                            adjacency[
                                static_cast<std::size_t>(reflected_row)
                            ][static_cast<std::size_t>(reflected_column)] =
                                entry;
                            adjacency[
                                static_cast<std::size_t>(reflected_column)
                            ][static_cast<std::size_t>(reflected_row)] =
                                entry;
                        }
                    }
                }
                adjacency[0][static_cast<std::size_t>(root_neighbor)] = 1;
                adjacency[static_cast<std::size_t>(root_neighbor)][0] = 1;
                adjacency[static_cast<std::size_t>(vertices - 1)]
                         [static_cast<std::size_t>(reflected_neighbor)] = 1;
                adjacency[static_cast<std::size_t>(reflected_neighbor)]
                         [static_cast<std::size_t>(vertices - 1)] = 1;
                adjacency[static_cast<std::size_t>(root_neighbor)]
                         [static_cast<std::size_t>(root_neighbor)] = 1;
                adjacency[static_cast<std::size_t>(reflected_neighbor)]
                         [static_cast<std::size_t>(reflected_neighbor)] = 1;
                if (construct_row_convex) {
                    row_convex_closure(adjacency);
                    if (!row_convex(adjacency)) {
                        throw std::runtime_error(
                            "row-convex closure failed"
                        );
                    }
                }
                if (!connected(adjacency)) {
                    continue;
                }
                if (construct_outside_in
                    && !outside_in_peo(adjacency)) {
                    throw std::runtime_error(
                        "constructed graph violates outside-in PEO"
                    );
                }
                if (require_outside_in
                    && !outside_in_peo(adjacency)) {
                    continue;
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
                    const int n = 2 * prefix + 2;
                    Integer current = 0;
                    for (int truncation = 0;
                         truncation <= prefix - 1;
                         ++truncation) {
                        const int even = 2 * truncation;
                        const int odd = even + 1;
                        current += binomial(n - 1, even)
                            * (
                                powers[static_cast<std::size_t>(even)][0]
                                    * powers[static_cast<std::size_t>(
                                        n - even
                                    )][static_cast<std::size_t>(
                                        vertices - 1
                                    )]
                                - powers[static_cast<std::size_t>(odd)][0]
                                    * powers[static_cast<std::size_t>(
                                        n - odd
                                    )][static_cast<std::size_t>(
                                        vertices - 1
                                    )]
                            );
                        if (truncation < 2) {
                            continue;
                        }
                        ++prefix_rows;
                        if (current < 0) {
                            std::cout
                                << "REFLECTED_LOOP_LEAF_RANDOM"
                                << " counterexample"
                                << " vertices=" << vertices
                                << " sample=" << sample
                                << " root_neighbor=" << root_neighbor
                                << " density_tenths=" << density_tenths
                                << " prefix=" << prefix
                                << " truncation=" << truncation
                                << " value=" << current;
                            print_graph(adjacency);
                            std::cout << " result=FAIL_RANDOM\n";
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }

        std::cout
            << "REFLECTED_LOOP_LEAF_RANDOM"
            << " maximum_vertices=" << maximum_vertices
            << " maximum_prefix=" << maximum_prefix
            << " samples_per_size=" << samples_per_size
            << " require_outside_in="
                << (require_outside_in ? 1 : 0)
            << " construct_outside_in="
                << (construct_outside_in ? 1 : 0)
            << " construct_row_convex="
                << (construct_row_convex ? 1 : 0)
            << " attempts=" << attempts
            << " connected_graphs=" << connected_graphs
            << " prefix_rows=" << prefix_rows
            << " result=PASS_RANDOM_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "REFLECTED_LOOP_LEAF_RANDOM FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
