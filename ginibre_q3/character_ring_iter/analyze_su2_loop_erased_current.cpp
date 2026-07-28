#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using Integer = boost::multiprecision::cpp_int;
using Series = std::vector<Integer>;

bool adjacent(int K, int Q, int x, int y) {
    return std::abs(x - Q) <= y
        && y <= std::min(x + Q, 2 * K - x - Q);
}

Integer binomial(int n, int r) {
    if (r < 0 || r > n) {
        return 0;
    }
    r = std::min(r, n - r);
    Integer result = 1;
    for (int i = 1; i <= r; ++i) {
        result *= n - r + i;
        result /= i;
    }
    return result;
}

Series product(const Series& a, const Series& b) {
    Series result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; i + j < a.size(); ++j) {
            result[i + j] += a[i] * b[j];
        }
    }
    return result;
}

Series returns(
    int K,
    int Q,
    int base,
    const std::vector<bool>& allowed,
    int maximum_degree
) {
    Series result(static_cast<std::size_t>(maximum_degree + 1));
    Series state(static_cast<std::size_t>(K + 1));
    state[static_cast<std::size_t>(base)] = 1;
    result[0] = 1;
    for (int degree = 1; degree <= maximum_degree; ++degree) {
        Series next(static_cast<std::size_t>(K + 1));
        for (int x = 1; x <= K; ++x) {
            if (!allowed[static_cast<std::size_t>(x)]
                || state[static_cast<std::size_t>(x)] == 0) {
                continue;
            }
            for (int y = 1; y <= K; ++y) {
                if (allowed[static_cast<std::size_t>(y)]
                    && adjacent(K, Q, x, y)) {
                    next[static_cast<std::size_t>(y)]
                        += state[static_cast<std::size_t>(x)];
                }
            }
        }
        state = std::move(next);
        result[static_cast<std::size_t>(degree)] =
            state[static_cast<std::size_t>(base)];
    }
    return result;
}

Series full_returns(int K, int Q, int maximum_degree) {
    Series result(static_cast<std::size_t>(maximum_degree + 1));
    Series state(static_cast<std::size_t>(K + 1));
    state[0] = 1;
    result[0] = 1;
    for (int degree = 1; degree <= maximum_degree; ++degree) {
        Series next(static_cast<std::size_t>(K + 1));
        for (int x = 0; x <= K; ++x) {
            for (int y = 0; y <= K; ++y) {
                if (adjacent(K, Q, x, y)) {
                    next[static_cast<std::size_t>(y)]
                        += state[static_cast<std::size_t>(x)];
                }
            }
        }
        state = std::move(next);
        result[static_cast<std::size_t>(degree)] = state[0];
    }
    return result;
}

Series full_endpoints(int K, int Q, int maximum_degree) {
    Series result(static_cast<std::size_t>(maximum_degree + 1));
    Series state(static_cast<std::size_t>(K + 1));
    state[0] = 1;
    for (int degree = 0; degree <= maximum_degree; ++degree) {
        result[static_cast<std::size_t>(degree)] =
            state[static_cast<std::size_t>(K)];
        Series next(static_cast<std::size_t>(K + 1));
        for (int x = 0; x <= K; ++x) {
            for (int y = 0; y <= K; ++y) {
                if (adjacent(K, Q, x, y)) {
                    next[static_cast<std::size_t>(y)]
                        += state[static_cast<std::size_t>(x)];
                }
            }
        }
        state = std::move(next);
    }
    return result;
}

int main(int argc, char** argv) {
    try {
    if (argc != 5 && argc != 6) {
        throw std::runtime_error("usage: K Q j t [--verbose]");
    }
    const bool verbose =
        argc == 6 && std::string(argv[5]) == "--verbose";
    if (argc == 6 && !verbose) {
        throw std::runtime_error("the optional argument is --verbose");
    }
    const int K = std::stoi(argv[1]);
    const int Q = std::stoi(argv[2]);
    const int j = std::stoi(argv[3]);
    const int t = std::stoi(argv[4]);
    if (Q < 1 || K <= 2 * Q || j < 4 || t < 0 || t >= j) {
        throw std::runtime_error("invalid parameters");
    }
    const int n = 2 * j + 2;
    const int m = n - 1;
    const Series f = full_returns(K, Q, n);
    const Series g = full_endpoints(K, Q, n);
    Series p(static_cast<std::size_t>(n + 1));
    for (int s = 0; s <= t; ++s) {
        p[static_cast<std::size_t>(2 * s)] +=
            binomial(m, 2 * s) * f[static_cast<std::size_t>(2 * s)];
        p[static_cast<std::size_t>(2 * s + 1)] -=
            binomial(m, 2 * s) * f[static_cast<std::size_t>(2 * s + 1)];
    }
    const Series h = product(p, f);
    Integer direct_target = 0;
    for (int degree = 0; degree <= n; ++degree) {
        direct_target +=
            p[static_cast<std::size_t>(degree)]
            * g[static_cast<std::size_t>(n - degree)];
    }

    const int start = Q;
    const int target = K;
    std::vector<bool> visited(static_cast<std::size_t>(K + 1), false);
    std::vector<int> path{start};
    visited[static_cast<std::size_t>(start)] = true;
    std::uint64_t skeletons = 0;
    std::uint64_t negative = 0;
    Integer minimum = 0;
    bool initialized = false;
    std::vector<int> minimum_path;
    Integer total = 0;
    Series path_sum(static_cast<std::size_t>(n + 1));
    std::vector<Integer> first_edge_totals(
        static_cast<std::size_t>(K + 1)
    );

    std::function<void(int)> dfs = [&](int vertex) {
        if (vertex == target) {
            ++skeletons;
            std::vector<bool> allowed(
                static_cast<std::size_t>(K + 1),
                false
            );
            for (int x = 1; x <= K; ++x) {
                allowed[static_cast<std::size_t>(x)] = true;
            }
            Series weight(static_cast<std::size_t>(n + 1));
            weight[0] = 1;
            for (const int x : path) {
                weight = product(
                    weight,
                    returns(K, Q, x, allowed, n)
                );
                allowed[static_cast<std::size_t>(x)] = false;
            }
            const int shift = static_cast<int>(path.size());
            for (int degree = 0; degree + shift <= n; ++degree) {
                path_sum[static_cast<std::size_t>(degree + shift)]
                    += weight[static_cast<std::size_t>(degree)];
            }
            Integer margin = 0;
            for (int degree = 0; degree + shift <= n; ++degree) {
                margin += weight[static_cast<std::size_t>(degree)]
                    * h[static_cast<std::size_t>(n - shift - degree)];
            }
            total += margin;
            first_edge_totals[
                static_cast<std::size_t>(path[1])
            ] += margin;
            if (verbose) {
                std::cout << "SKELETON path=";
                for (const int x : path) {
                    std::cout << x << ',';
                }
                std::cout << " margin=" << margin << '\n';
            }
            if (!initialized || margin < minimum) {
                initialized = true;
                minimum = margin;
                minimum_path = path;
            }
            if (margin < 0) {
                ++negative;
                if (negative == 1) {
                    std::cout << "FIRST_NEGATIVE skeleton=";
                    for (const int x : path) {
                        std::cout << x << ',';
                    }
                    std::cout << " margin=" << margin << '\n';
                }
            }
            return;
        }
        for (int next = 1; next <= K; ++next) {
            if (!visited[static_cast<std::size_t>(next)]
                && adjacent(K, Q, vertex, next)) {
                visited[static_cast<std::size_t>(next)] = true;
                path.push_back(next);
                dfs(next);
                path.pop_back();
                visited[static_cast<std::size_t>(next)] = false;
            }
        }
    };
    dfs(start);
    Series interior_state(static_cast<std::size_t>(K + 1));
    Series direct_e(static_cast<std::size_t>(n + 1));
    interior_state[static_cast<std::size_t>(start)] = 1;
    for (int degree = 0; degree + 1 <= n; ++degree) {
        direct_e[static_cast<std::size_t>(degree + 1)] =
            interior_state[static_cast<std::size_t>(target)];
        Series next(static_cast<std::size_t>(K + 1));
        for (int x = 1; x <= K; ++x) {
            for (int y = 1; y <= K; ++y) {
                if (adjacent(K, Q, x, y)) {
                    next[static_cast<std::size_t>(y)]
                        += interior_state[static_cast<std::size_t>(x)];
                }
            }
        }
        interior_state = std::move(next);
    }
    for (int degree = 0; degree <= n; ++degree) {
        if (path_sum[static_cast<std::size_t>(degree)]
            != direct_e[static_cast<std::size_t>(degree)]) {
            throw std::runtime_error(
                "loop-erased path-sum identity failed at degree "
                + std::to_string(degree)
            );
        }
    }
    const Series reconstructed_g = product(f, direct_e);
    for (int degree = 0; degree <= n; ++degree) {
        if (reconstructed_g[static_cast<std::size_t>(degree)]
            != g[static_cast<std::size_t>(degree)]) {
            throw std::runtime_error(
                "last-exit identity failed at degree "
                + std::to_string(degree)
            );
        }
    }
    if (total != direct_target) {
        throw std::runtime_error(
            "loop-erased current sum does not equal direct current"
        );
    }
    std::uint64_t negative_first_edge_groups = 0U;
    Integer minimum_first_edge_group = 0;
    bool initialized_first_edge_group = false;
    for (int neighbor = 1; neighbor <= K; ++neighbor) {
        const Integer& margin =
            first_edge_totals[static_cast<std::size_t>(neighbor)];
        if (margin == 0) {
            continue;
        }
        if (
            !initialized_first_edge_group
            || margin < minimum_first_edge_group
        ) {
            initialized_first_edge_group = true;
            minimum_first_edge_group = margin;
        }
        if (margin < 0) {
            ++negative_first_edge_groups;
        }
        if (verbose) {
            std::cout << "FIRST_EDGE_GROUP neighbor=" << neighbor
                      << " margin=" << margin << '\n';
        }
    }
    std::cout << "LOOP_ERASED_CURRENT"
              << " K=" << K
              << " Q=" << Q
              << " j=" << j
              << " t=" << t
              << " skeletons=" << skeletons
              << " negative=" << negative
              << " minimum=" << minimum
              << " minimum_skeleton=";
    for (const int vertex : minimum_path) {
        std::cout << vertex << ',';
    }
    std::cout
              << " negative_first_edge_groups="
              << negative_first_edge_groups
              << " minimum_first_edge_group="
              << minimum_first_edge_group
              << " total=" << total
              << " direct=" << direct_target
              << " result="
              << (
                    negative == 0
                        ? "PASS_EXACT_DECOMPOSITION_SKELETONWISE"
                        : "PASS_EXACT_DECOMPOSITION_WITH_SKELETON_NEGATIVES"
                 )
              << '\n';
    return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_LOOP_ERASED_CURRENT FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
