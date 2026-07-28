#include <algorithm>
#include <cstdlib>
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

std::vector<std::vector<int>> neighbors(int K, int Q) {
    std::vector<std::vector<int>> result(
        static_cast<std::size_t>(K + 1)
    );
    for (int x = 0; x <= K; ++x) {
        for (int y = 0; y <= K; ++y) {
            if (adjacent(K, Q, x, y)) {
                result[static_cast<std::size_t>(x)].push_back(y);
            }
        }
    }
    return result;
}

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error("usage: MAXIMUM_K MAXIMUM_PREFIX");
        }
        const int maximum_k = std::stoi(argv[1]);
        const int maximum_prefix = std::stoi(argv[2]);
        if (maximum_k < 3 || maximum_prefix < 5) {
            throw std::runtime_error("bounds are too small");
        }
        const int maximum_degree = 2 * maximum_prefix + 2;
        std::uint64_t parameters = 0;
        std::uint64_t rows = 0;
        std::uint64_t identities = 0;
        std::uint64_t groups = 0;
        std::uint64_t negative_groups = 0;
        std::uint64_t nonspecial_negative_groups = 0;
        std::uint64_t special_positive_groups = 0;
        std::uint64_t outside_in_prefixes = 0;
        std::uint64_t negative_outside_in_prefixes = 0;
        Integer minimum_nonspecial = 0;
        bool initialized_nonspecial = false;
        Integer minimum_total = 0;
        bool initialized_total = false;

        for (int K = 3; K <= maximum_k; ++K) {
            for (int Q = 1; K > 2 * Q; ++Q) {
                ++parameters;
                const auto graph = neighbors(K, Q);
                Series f(static_cast<std::size_t>(maximum_degree + 1));
                Series g(static_cast<std::size_t>(maximum_degree + 1));
                Series state(static_cast<std::size_t>(K + 1));
                state[0] = 1;
                f[0] = 1;
                for (int degree = 1;
                     degree <= maximum_degree;
                     ++degree) {
                    Series next(static_cast<std::size_t>(K + 1));
                    for (int x = 0; x <= K; ++x) {
                        for (const int y :
                             graph[static_cast<std::size_t>(x)]) {
                            next[static_cast<std::size_t>(y)]
                                += state[static_cast<std::size_t>(x)];
                        }
                    }
                    state = std::move(next);
                    f[static_cast<std::size_t>(degree)] = state[0];
                    g[static_cast<std::size_t>(degree)] =
                        state[static_cast<std::size_t>(K)];
                }

                Series positive_return(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                state.assign(static_cast<std::size_t>(K + 1), Integer(0));
                state[static_cast<std::size_t>(Q)] = 1;
                positive_return[0] = 1;
                for (int degree = 1;
                     degree <= maximum_degree;
                     ++degree) {
                    Series next(static_cast<std::size_t>(K + 1));
                    for (int x = 1; x <= K; ++x) {
                        for (const int y :
                             graph[static_cast<std::size_t>(x)]) {
                            if (y >= 1) {
                                next[static_cast<std::size_t>(y)]
                                    += state[static_cast<std::size_t>(x)];
                            }
                        }
                    }
                    state = std::move(next);
                    positive_return[static_cast<std::size_t>(degree)] =
                        state[static_cast<std::size_t>(Q)];
                }

                std::vector<Series> escape(
                    static_cast<std::size_t>(K + 1),
                    Series(static_cast<std::size_t>(maximum_degree + 1))
                );
                state.assign(static_cast<std::size_t>(K + 1), Integer(0));
                state[static_cast<std::size_t>(K)] = 1;
                for (int degree = 0;
                     degree <= maximum_degree;
                     ++degree) {
                    for (int y = 1; y <= K; ++y) {
                        if (y != Q) {
                            escape[static_cast<std::size_t>(y)]
                                  [static_cast<std::size_t>(degree)] =
                                state[static_cast<std::size_t>(y)];
                        }
                    }
                    Series next(static_cast<std::size_t>(K + 1));
                    for (int x = 1; x <= K; ++x) {
                        if (x == Q) {
                            continue;
                        }
                        for (const int y :
                             graph[static_cast<std::size_t>(x)]) {
                            if (y >= 1 && y != Q) {
                                next[static_cast<std::size_t>(y)]
                                    += state[static_cast<std::size_t>(x)];
                            }
                        }
                    }
                    state = std::move(next);
                }

                for (int j = 5; j <= maximum_prefix; ++j) {
                    const int n = 2 * j + 2;
                    Series p(static_cast<std::size_t>(maximum_degree + 1));
                    for (int t = 0; t <= j - 1; ++t) {
                        p[static_cast<std::size_t>(2 * t)] +=
                            binomial(n - 1, 2 * t)
                            * f[static_cast<std::size_t>(2 * t)];
                        p[static_cast<std::size_t>(2 * t + 1)] -=
                            binomial(n - 1, 2 * t)
                            * f[static_cast<std::size_t>(2 * t + 1)];
                        if (t < 4) {
                            continue;
                        }
                        ++rows;
                        const Series common =
                            product(product(p, f), positive_return);
                        Integer total = 0;
                        std::vector<Integer> distance_buckets(
                            static_cast<std::size_t>(K + 1)
                        );
                        int maximum_distance = 0;
                        for (const int y :
                             graph[static_cast<std::size_t>(Q)]) {
                            if (y == 0 || y == Q) {
                                continue;
                            }
                            ++groups;
                            Integer margin = 0;
                            for (int degree = 0;
                                 degree <= n - 2;
                                 ++degree) {
                                margin += common[
                                    static_cast<std::size_t>(degree)
                                ] * escape[static_cast<std::size_t>(y)][
                                    static_cast<std::size_t>(
                                        n - 2 - degree
                                    )];
                            }
                            total += margin;
                            const int distance = std::abs(y - (K - Q));
                            distance_buckets[
                                static_cast<std::size_t>(distance)
                            ] += margin;
                            maximum_distance =
                                std::max(maximum_distance, distance);
                            const bool special = y == K - Q;
                            if (!special
                                && (
                                    !initialized_nonspecial
                                    || margin < minimum_nonspecial
                                )) {
                                initialized_nonspecial = true;
                                minimum_nonspecial = margin;
                            }
                            if (margin < 0) {
                                ++negative_groups;
                                if (!special) {
                                    ++nonspecial_negative_groups;
                                    if (nonspecial_negative_groups == 1) {
                                        std::cout
                                            << "FIRST_NEGATIVE_NONSPECIAL"
                                            << " K=" << K
                                            << " Q=" << Q
                                            << " j=" << j
                                            << " t=" << t
                                            << " y=" << y
                                            << " margin=" << margin
                                            << '\n';
                                    }
                                }
                            } else if (special) {
                                ++special_positive_groups;
                            }
                        }
                        Integer outside_in = 0;
                        for (int distance = maximum_distance;
                             distance >= 0;
                             --distance) {
                            outside_in += distance_buckets[
                                static_cast<std::size_t>(distance)
                            ];
                            ++outside_in_prefixes;
                            if (outside_in < 0) {
                                ++negative_outside_in_prefixes;
                                if (negative_outside_in_prefixes == 1) {
                                    std::cout
                                        << "FIRST_NEGATIVE_OUTSIDE_IN"
                                        << " K=" << K
                                        << " Q=" << Q
                                        << " j=" << j
                                        << " t=" << t
                                        << " distance=" << distance
                                        << " margin=" << outside_in
                                        << '\n';
                                }
                            }
                        }
                        Integer direct = 0;
                        for (int degree = 0; degree <= n; ++degree) {
                            direct += p[static_cast<std::size_t>(degree)]
                                * g[static_cast<std::size_t>(
                                    n - degree
                                )];
                        }
                        if (total != direct) {
                            throw std::runtime_error(
                                "last-exit group identity failed"
                            );
                        }
                        ++identities;
                        if (!initialized_total || total < minimum_total) {
                            initialized_total = true;
                            minimum_total = total;
                        }
                        if (total < 0) {
                            throw std::runtime_error(
                                "negative full-prefix total"
                            );
                        }
                    }
                }
            }
        }
        std::cout
            << "SU2_LAST_EXIT_GROUP_SCAN"
            << " maximum_K=" << maximum_k
            << " maximum_prefix=" << maximum_prefix
            << " parameters=" << parameters
            << " rows=" << rows
            << " identities=" << identities
            << " groups=" << groups
            << " negative_groups=" << negative_groups
            << " nonspecial_negative_groups="
            << nonspecial_negative_groups
            << " special_positive_groups=" << special_positive_groups
            << " outside_in_prefixes=" << outside_in_prefixes
            << " negative_outside_in_prefixes="
            << negative_outside_in_prefixes
            << " minimum_nonspecial=" << minimum_nonspecial
            << " minimum_total=" << minimum_total
            << " result="
            << (
                negative_outside_in_prefixes == 0
                    ? "PASS_OUTSIDE_IN_PREFIXES"
                    : "FAIL_OUTSIDE_IN_PREFIXES"
            )
            << '\n';
        return negative_outside_in_prefixes == 0
            ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "SU2_LAST_EXIT_GROUP_SCAN FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
