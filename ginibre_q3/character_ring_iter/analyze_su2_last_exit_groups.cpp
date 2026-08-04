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
        std::uint64_t special_negative_groups = 0;
        std::uint64_t special_positive_groups = 0;
        std::uint64_t outside_in_prefixes = 0;
        std::uint64_t negative_outside_in_prefixes = 0;
        std::uint64_t renewal_current_coefficients = 0;
        std::uint64_t negative_renewal_current_coefficients = 0;
        std::uint64_t twice_return_current_coefficients = 0;
        std::uint64_t negative_twice_return_current_coefficients = 0;
        std::uint64_t green_minor_identities = 0;
        std::uint64_t outside_in_core_identities = 0;
        std::uint64_t outside_in_core_coefficients = 0;
        std::uint64_t negative_outside_in_core_coefficients = 0;
        std::uint64_t negative_higher_q_core_coefficients = 0;
        std::uint64_t negative_higher_q_twice_return_coefficients = 0;
        std::vector<std::uint64_t> interval_piece_tests(6, 0);
        std::vector<std::uint64_t> negative_interval_piece_tests(6, 0);
        std::uint64_t simultaneous_negative_upper_pieces = 0;
        std::uint64_t simultaneous_negative_ray_pairs = 0;
        std::uint64_t strict_simultaneous_negative_upper_pieces = 0;
        std::uint64_t strict_negative_left_rays = 0;
        std::uint64_t strict_deficit_reserve_tests = 0;
        std::uint64_t negative_strict_deficit_reserves = 0;
        std::uint64_t local_triple_reserve_tests = 0;
        std::uint64_t negative_local_triple_reserves = 0;
        std::uint64_t nested_radius_reserve_tests = 0;
        std::uint64_t negative_nested_radius_reserves = 0;
        std::uint64_t odd_boundary_current_coefficients = 0;
        std::uint64_t negative_odd_boundary_current_coefficients = 0;
        std::uint64_t negative_higher_q_odd_boundary_coefficients = 0;
        std::uint64_t t4_kernel_coordinates = 0;
        std::uint64_t negative_t4_kernel_coordinates = 0;
        std::uint64_t negative_higher_q_t4_kernel_coordinates = 0;
        std::uint64_t full_t4_kernel_coordinates = 0;
        std::uint64_t negative_full_t4_kernel_coordinates = 0;
        std::uint64_t t4_induction_coordinates = 0;
        std::uint64_t negative_t4_induction_coordinates = 0;
        std::uint64_t t4_base_difference_coordinates = 0;
        std::uint64_t negative_t4_base_difference_coordinates = 0;
        std::uint64_t nonzero_t4_ninth_difference_coordinates = 0;
        std::uint64_t t4_c4_coordinates = 0;
        std::uint64_t negative_t4_c4_coordinates = 0;
        std::uint64_t t4_grouped_block_coordinates = 0;
        std::uint64_t negative_t4_grouped_block_coordinates = 0;
        std::uint64_t negative_t4_low_component_coordinates = 0;
        std::uint64_t negative_t4_weak_group_coordinates = 0;
        int last_negative_left_ray_K = -1;
        int last_negative_left_ray_Q = -1;
        int last_negative_left_ray_distance = -1;
        const char* const interval_piece_names[6] = {
            "LOWER_LEFT",
            "UPPER_LEFT",
            "UPPER_RIGHT",
            "LEFT_RAY",
            "OUTER_FLANKS",
            "UPPER_PAIR"
        };
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
                std::vector<Series> green(
                    static_cast<std::size_t>(K + 1),
                    Series(static_cast<std::size_t>(maximum_degree + 1))
                );
                Series state(static_cast<std::size_t>(K + 1));
                state[0] = 1;
                f[0] = 1;
                green[0][0] = 1;
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
                    for (int x = 0; x <= K; ++x) {
                        green[static_cast<std::size_t>(x)]
                             [static_cast<std::size_t>(degree)] =
                            state[static_cast<std::size_t>(x)];
                    }
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

                Series nonempty_returns = f;
                nonempty_returns[0] = 0;
                std::vector<Series> first_passage(
                    static_cast<std::size_t>(K + 1),
                    Series(static_cast<std::size_t>(maximum_degree + 1))
                );
                for (int x = 0; x <= K; ++x) {
                    for (int degree = 0;
                         degree <= maximum_degree;
                         ++degree) {
                        Integer value =
                            green[static_cast<std::size_t>(x)]
                                 [static_cast<std::size_t>(degree)];
                        for (int return_degree = 1;
                             return_degree <= degree;
                             ++return_degree) {
                            value -= f[
                                static_cast<std::size_t>(return_degree)
                            ] * first_passage[
                                static_cast<std::size_t>(x)
                            ][static_cast<std::size_t>(
                                degree - return_degree
                            )];
                        }
                        if (value < 0) {
                            throw std::runtime_error(
                                "negative first-passage coefficient"
                            );
                        }
                        first_passage[static_cast<std::size_t>(x)]
                                     [static_cast<std::size_t>(degree)] =
                            value;
                    }
                }

                int maximum_channel_distance = 0;
                std::vector<Series> direct_buckets(
                    static_cast<std::size_t>(K + 1),
                    Series(static_cast<std::size_t>(maximum_degree + 1))
                );
                std::vector<Series> reflected_buckets = direct_buckets;
                std::vector<Series> escape_buckets = direct_buckets;
                for (const int y :
                     graph[static_cast<std::size_t>(Q)]) {
                    if (y == 0 || y == Q) {
                        continue;
                    }
                    const Series left =
                        product(nonempty_returns, escape[
                            static_cast<std::size_t>(y)
                        ]);
                    const Series reflected = product(
                        nonempty_returns,
                        green[static_cast<std::size_t>(K - y)]
                    );
                    const Series correction = product(
                        green[static_cast<std::size_t>(y)],
                        green[static_cast<std::size_t>(K - Q)]
                    );
                    for (int degree = 0;
                         degree <= maximum_degree;
                         ++degree) {
                        Integer right =
                            reflected[static_cast<std::size_t>(degree)];
                        if (degree >= 1) {
                            right -= correction[
                                static_cast<std::size_t>(degree - 1)
                            ];
                        }
                        if (
                            left[static_cast<std::size_t>(degree)]
                            != right
                        ) {
                            throw std::runtime_error(
                                "two-vertex Green minor identity failed"
                            );
                        }
                    }
                    ++green_minor_identities;

                    const int distance = std::abs(y - (K - Q));
                    maximum_channel_distance =
                        std::max(maximum_channel_distance, distance);
                    for (int degree = 0;
                         degree <= maximum_degree;
                         ++degree) {
                        direct_buckets[
                            static_cast<std::size_t>(distance)
                        ][static_cast<std::size_t>(degree)] +=
                            first_passage[static_cast<std::size_t>(y)]
                                         [static_cast<std::size_t>(degree)];
                        reflected_buckets[
                            static_cast<std::size_t>(distance)
                        ][static_cast<std::size_t>(degree)] +=
                            first_passage[
                                static_cast<std::size_t>(K - y)
                            ][static_cast<std::size_t>(degree)];
                        escape_buckets[
                            static_cast<std::size_t>(distance)
                        ][static_cast<std::size_t>(degree)] +=
                            escape[static_cast<std::size_t>(y)]
                                  [static_cast<std::size_t>(degree)];
                    }
                }

                Series direct_tail(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                Series reflected_tail = direct_tail;
                Series escape_tail = direct_tail;
                const Series twice_returns = product(f, f);
                for (int distance = maximum_channel_distance;
                     distance >= 0;
                     --distance) {
                    for (int degree = 0;
                         degree <= maximum_degree;
                         ++degree) {
                        direct_tail[static_cast<std::size_t>(degree)] +=
                            direct_buckets[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(degree)];
                        reflected_tail[static_cast<std::size_t>(degree)] +=
                            reflected_buckets[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(degree)];
                        escape_tail[static_cast<std::size_t>(degree)] +=
                            escape_buckets[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(degree)];
                    }
                    Series core = product(
                        first_passage[static_cast<std::size_t>(Q)],
                        reflected_tail
                    );
                    const Series core_correction = product(
                        first_passage[
                            static_cast<std::size_t>(K - Q)
                        ],
                        direct_tail
                    );
                    for (int degree = 0;
                         degree <= maximum_degree;
                         ++degree) {
                        core[static_cast<std::size_t>(degree)] -=
                            core_correction[
                                static_cast<std::size_t>(degree)
                            ];
                        ++outside_in_core_coefficients;
                        if (core[static_cast<std::size_t>(degree)] < 0) {
                            ++negative_outside_in_core_coefficients;
                            if (
                                negative_outside_in_core_coefficients
                                == 1
                            ) {
                                std::cout
                                    << "FIRST_NEGATIVE_OUTSIDE_IN_CORE"
                                    << " K=" << K
                                    << " Q=" << Q
                                    << " distance=" << distance
                                    << " degree=" << degree
                                    << " value=" << core[
                                        static_cast<std::size_t>(degree)
                                    ]
                                    << '\n';
                            }
                            if (Q >= 2) {
                                ++negative_higher_q_core_coefficients;
                                if (
                                    negative_higher_q_core_coefficients
                                    == 1
                                ) {
                                    std::cout
                                        << "FIRST_NEGATIVE_HIGHER_Q_CORE"
                                        << " K=" << K
                                        << " Q=" << Q
                                        << " distance=" << distance
                                        << " degree=" << degree
                                        << " value=" << core[
                                            static_cast<std::size_t>(degree)
                                        ]
                                        << '\n';
                                }
                            }
                        }
                    }
                    const Series reconstructed =
                        product(twice_returns, core);
                    const Series killed = product(
                        green[static_cast<std::size_t>(Q)],
                        escape_tail
                    );
                    for (int degree = 0;
                         degree <= maximum_degree;
                         ++degree) {
                        if (
                            reconstructed[
                                static_cast<std::size_t>(degree)
                            ] != killed[static_cast<std::size_t>(degree)]
                        ) {
                            throw std::runtime_error(
                                "outside-in Green core identity failed"
                            );
                        }
                    }
                    ++outside_in_core_identities;
                }

                std::vector<Series> endpoint_powers(
                    11,
                    Series(static_cast<std::size_t>(K + 1))
                );
                endpoint_powers[0][static_cast<std::size_t>(K)] = 1;
                for (int power = 1; power <= 10; ++power) {
                    for (int x = 0; x <= K; ++x) {
                        for (const int y :
                             graph[static_cast<std::size_t>(x)]) {
                            endpoint_powers[
                                static_cast<std::size_t>(power)
                            ][static_cast<std::size_t>(y)] +=
                                endpoint_powers[
                                    static_cast<std::size_t>(power - 1)
                                ][static_cast<std::size_t>(x)];
                        }
                    }
                }
                for (int x = 0; x <= K; ++x) {
                    const std::size_t index =
                        static_cast<std::size_t>(K - x);
                    const Integer c4 =
                        f[8] * endpoint_powers[9][index]
                        - f[9] * endpoint_powers[8][index];
                    ++t4_c4_coordinates;
                    if (c4 < 0) {
                        ++negative_t4_c4_coordinates;
                        if (negative_t4_c4_coordinates == 1) {
                            std::cout
                                << "FIRST_NEGATIVE_T4_C4"
                                << " K=" << K
                                << " Q=" << Q
                                << " x=" << x
                                << " value=" << c4 << '\n';
                        }
                    }
                    const Integer component_2_r0 =
                        f[4] * endpoint_powers[8][index]
                        - f[5] * endpoint_powers[7][index];
                    const Integer component_3_r0 =
                        f[6] * endpoint_powers[6][index]
                        - f[7] * endpoint_powers[5][index];
                    const Integer component_4_r0 =
                        f[8] * endpoint_powers[4][index]
                        - f[9] * endpoint_powers[3][index];
                    const Integer component_2_r1 =
                        f[4] * endpoint_powers[10][index]
                        - f[5] * endpoint_powers[9][index];
                    const Integer component_3_r1 =
                        f[6] * endpoint_powers[8][index]
                        - f[7] * endpoint_powers[7][index];
                    const Integer component_4_r1 =
                        f[8] * endpoint_powers[6][index]
                        - f[9] * endpoint_powers[5][index];
                    const Integer component_3_r2 =
                        f[6] * endpoint_powers[10][index]
                        - f[7] * endpoint_powers[9][index];
                    const Integer component_4_r2 =
                        f[8] * endpoint_powers[8][index]
                        - f[9] * endpoint_powers[7][index];
                    const Integer grouped_blocks[3] = {
                        330 * component_2_r0
                            + 462 * component_3_r0
                            + 165 * component_4_r0,
                        385 * component_2_r1
                            + 1254 * component_3_r1
                            + 1122 * component_4_r1,
                        2035 * component_3_r2
                            + 4026 * component_4_r2
                    };
                    for (int order = 0; order <= 2; ++order) {
                        ++t4_grouped_block_coordinates;
                        if (grouped_blocks[order] < 0) {
                            ++negative_t4_grouped_block_coordinates;
                            if (
                                negative_t4_grouped_block_coordinates
                                == 1
                            ) {
                                std::cout
                                    << "FIRST_NEGATIVE_T4_GROUPED_BLOCK"
                                    << " K=" << K
                                    << " Q=" << Q
                                    << " order=" << order
                                    << " x=" << x
                                    << " value=" << grouped_blocks[order]
                                    << '\n';
                            }
                        }
                    }
                    if (component_4_r0 < 0) {
                        ++negative_t4_low_component_coordinates;
                        if (
                            negative_t4_low_component_coordinates
                            == 1
                        ) {
                            std::cout
                                << "FIRST_NEGATIVE_T4_LOW_COMPONENT"
                                << " K=" << K
                                << " Q=" << Q
                                << " order=0 s=4"
                                << " x=" << x
                                << " value=" << component_4_r0 << '\n';
                        }
                    }
                    const Integer weak_group =
                        1254 * component_3_r1
                        + 1122 * component_4_r1;
                    if (weak_group < 0) {
                        ++negative_t4_weak_group_coordinates;
                        if (negative_t4_weak_group_coordinates == 1) {
                            std::cout
                                << "FIRST_NEGATIVE_T4_WEAK_GROUP"
                                << " K=" << K
                                << " Q=" << Q
                                << " order=1"
                                << " x=" << x
                                << " value=" << weak_group << '\n';
                        }
                    }
                }

                Series previous_full_t4_kernel;
                std::vector<Series> full_t4_kernels;
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
                        if (t == 4) {
                            Series t4_kernel(
                                static_cast<std::size_t>(K + 1)
                            );
                            for (int x = 0; x <= K; ++x) {
                                const std::size_t index =
                                    static_cast<std::size_t>(x);
                                t4_kernel[index] =
                                    endpoint_powers[9][index]
                                    + binomial(n - 1, 2)
                                        * (
                                            f[2] * endpoint_powers[7][index]
                                            - f[3]
                                                * endpoint_powers[6][index]
                                        )
                                    + binomial(n - 1, 4)
                                        * (
                                            f[4] * endpoint_powers[5][index]
                                            - f[5]
                                                * endpoint_powers[4][index]
                                        )
                                    + binomial(n - 1, 6)
                                        * (
                                            f[6] * endpoint_powers[3][index]
                                            - f[7]
                                                * endpoint_powers[2][index]
                                        )
                                    + binomial(n - 1, 8)
                                        * (
                                            f[8] * endpoint_powers[1][index]
                                            - f[9]
                                                * endpoint_powers[0][index]
                                        );
                            }
                            for (int step = 0; step < 3; ++step) {
                                Series next(
                                    static_cast<std::size_t>(K + 1)
                                );
                                for (int x = 0; x <= K; ++x) {
                                    for (const int y :
                                         graph[
                                             static_cast<std::size_t>(x)
                                         ]) {
                                        next[static_cast<std::size_t>(y)] +=
                                            t4_kernel[
                                                static_cast<std::size_t>(x)
                                            ];
                                    }
                                }
                                t4_kernel = std::move(next);
                            }
                            for (int x = 0; x <= K; ++x) {
                                ++t4_kernel_coordinates;
                                if (
                                    t4_kernel[static_cast<std::size_t>(x)]
                                    < 0
                                ) {
                                    ++negative_t4_kernel_coordinates;
                                    if (
                                        negative_t4_kernel_coordinates == 1
                                    ) {
                                        std::cout
                                            << "FIRST_NEGATIVE_T4_KERNEL"
                                            << " K=" << K
                                            << " Q=" << Q
                                            << " j=" << j
                                            << " x=" << x
                                            << " value=" << t4_kernel[
                                                static_cast<std::size_t>(x)
                                            ]
                                            << '\n';
                                    }
                                    if (Q >= 2) {
                                        ++negative_higher_q_t4_kernel_coordinates;
                                        if (
                                            negative_higher_q_t4_kernel_coordinates
                                            == 1
                                        ) {
                                            std::cout
                                                << "FIRST_NEGATIVE_HIGHER_Q_T4_KERNEL"
                                                << " K=" << K
                                                << " Q=" << Q
                                                << " j=" << j
                                                << " x=" << x
                                                << " value=" << t4_kernel[
                                                    static_cast<std::size_t>(
                                                        x
                                                    )
                                                ]
                                                << '\n';
                                        }
                                    }
                                }
                            }
                            Series full_t4_kernel = t4_kernel;
                            for (int step = 3;
                                 step < n - 9;
                                 ++step) {
                                Series next(
                                    static_cast<std::size_t>(K + 1)
                                );
                                for (int x = 0; x <= K; ++x) {
                                    for (const int y :
                                         graph[
                                             static_cast<std::size_t>(x)
                                         ]) {
                                        next[static_cast<std::size_t>(y)] +=
                                            full_t4_kernel[
                                                static_cast<std::size_t>(x)
                                            ];
                                    }
                                }
                                full_t4_kernel = std::move(next);
                            }
                            for (int x = 0; x <= K; ++x) {
                                ++full_t4_kernel_coordinates;
                                if (
                                    full_t4_kernel[
                                        static_cast<std::size_t>(x)
                                    ] < 0
                                ) {
                                    ++negative_full_t4_kernel_coordinates;
                                    if (
                                        negative_full_t4_kernel_coordinates
                                        == 1
                                    ) {
                                        std::cout
                                            << "FIRST_NEGATIVE_FULL_T4_KERNEL"
                                            << " K=" << K
                                            << " Q=" << Q
                                            << " j=" << j
                                            << " x=" << x
                                            << " value=" << full_t4_kernel[
                                                static_cast<std::size_t>(x)
                                            ]
                                            << '\n';
                                    }
                                }
                            }
                            if (!previous_full_t4_kernel.empty()) {
                                Series advanced =
                                    previous_full_t4_kernel;
                                for (int step = 0; step < 2; ++step) {
                                    Series next(
                                        static_cast<std::size_t>(K + 1)
                                    );
                                    for (int x = 0; x <= K; ++x) {
                                        for (const int y :
                                             graph[
                                                 static_cast<std::size_t>(x)
                                             ]) {
                                            next[
                                                static_cast<std::size_t>(y)
                                            ] += advanced[
                                                static_cast<std::size_t>(x)
                                            ];
                                        }
                                    }
                                    advanced = std::move(next);
                                }
                                for (int x = 0; x <= K; ++x) {
                                    ++t4_induction_coordinates;
                                    const Integer increment =
                                        full_t4_kernel[
                                            static_cast<std::size_t>(x)
                                        ] - advanced[
                                            static_cast<std::size_t>(x)
                                        ];
                                    if (increment < 0) {
                                        ++negative_t4_induction_coordinates;
                                        if (
                                            negative_t4_induction_coordinates
                                            == 1
                                        ) {
                                            std::cout
                                                << "FIRST_NEGATIVE_T4_INDUCTION"
                                                << " K=" << K
                                                << " Q=" << Q
                                                << " j=" << j
                                                << " x=" << x
                                                << " value=" << increment
                                                << '\n';
                                        }
                                    }
                                }
                            }
                            previous_full_t4_kernel = full_t4_kernel;
                            full_t4_kernels.push_back(full_t4_kernel);
                        }
                        ++rows;
                        const Series renewal_current =
                            product(p, nonempty_returns);
                        const Series twice_return_current =
                            product(product(p, f), f);
                        Series odd_boundary = green[
                            static_cast<std::size_t>(Q)
                        ];
                        for (int degree = 0;
                             degree <= maximum_degree;
                             ++degree) {
                            odd_boundary[
                                static_cast<std::size_t>(degree)
                            ] -= green[
                                static_cast<std::size_t>(K - Q)
                            ][static_cast<std::size_t>(degree)];
                        }
                        const Series odd_boundary_current =
                            product(p, odd_boundary);
                        for (int degree = 0; degree <= n; ++degree) {
                            ++renewal_current_coefficients;
                            if (renewal_current[
                                    static_cast<std::size_t>(degree)
                                ] < 0) {
                                ++negative_renewal_current_coefficients;
                                if (
                                    negative_renewal_current_coefficients
                                    == 1
                                ) {
                                    std::cout
                                        << "FIRST_NEGATIVE_RENEWAL_CURRENT"
                                        << " K=" << K
                                        << " Q=" << Q
                                        << " j=" << j
                                        << " t=" << t
                                        << " degree=" << degree
                                        << " value=" << renewal_current[
                                            static_cast<std::size_t>(degree)
                                        ]
                                        << '\n';
                                }
                            }
                            ++twice_return_current_coefficients;
                            if (twice_return_current[
                                    static_cast<std::size_t>(degree)
                                ] < 0) {
                                ++negative_twice_return_current_coefficients;
                                if (
                                    negative_twice_return_current_coefficients
                                    == 1
                                ) {
                                    std::cout
                                        << "FIRST_NEGATIVE_TWICE_RETURN_CURRENT"
                                        << " K=" << K
                                        << " Q=" << Q
                                        << " j=" << j
                                        << " t=" << t
                                        << " degree=" << degree
                                        << " value=" << twice_return_current[
                                            static_cast<std::size_t>(degree)
                                        ]
                                        << '\n';
                                }
                                if (Q >= 2) {
                                    ++negative_higher_q_twice_return_coefficients;
                                    if (
                                        negative_higher_q_twice_return_coefficients
                                        == 1
                                    ) {
                                        std::cout
                                            << "FIRST_NEGATIVE_HIGHER_Q_TWICE_RETURN"
                                            << " K=" << K
                                            << " Q=" << Q
                                            << " j=" << j
                                            << " t=" << t
                                            << " degree=" << degree
                                            << " value="
                                            << twice_return_current[
                                                static_cast<std::size_t>(
                                                    degree
                                                )
                                            ]
                                            << '\n';
                                    }
                                }
                            }
                            ++odd_boundary_current_coefficients;
                            if (odd_boundary_current[
                                    static_cast<std::size_t>(degree)
                                ] < 0) {
                                ++negative_odd_boundary_current_coefficients;
                                if (
                                    negative_odd_boundary_current_coefficients
                                    == 1
                                ) {
                                    std::cout
                                        << "FIRST_NEGATIVE_ODD_BOUNDARY_CURRENT"
                                        << " K=" << K
                                        << " Q=" << Q
                                        << " j=" << j
                                        << " t=" << t
                                        << " degree=" << degree
                                        << " value="
                                        << odd_boundary_current[
                                            static_cast<std::size_t>(degree)
                                        ]
                                        << '\n';
                                }
                                if (Q >= 2) {
                                    ++negative_higher_q_odd_boundary_coefficients;
                                    if (
                                        negative_higher_q_odd_boundary_coefficients
                                        == 1
                                    ) {
                                        std::cout
                                            << "FIRST_NEGATIVE_HIGHER_Q_ODD_BOUNDARY_CURRENT"
                                            << " K=" << K
                                            << " Q=" << Q
                                            << " j=" << j
                                            << " t=" << t
                                            << " degree=" << degree
                                            << " value="
                                            << odd_boundary_current[
                                                static_cast<std::size_t>(
                                                    degree
                                                )
                                            ]
                                            << '\n';
                                    }
                                }
                            }
                        }
                        const Series common =
                            product(product(p, f), positive_return);
                        Integer total = 0;
                        std::vector<Integer> distance_buckets(
                            static_cast<std::size_t>(K + 1)
                        );
                        std::vector<Integer> channel_margins(
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
                            channel_margins[
                                static_cast<std::size_t>(y)
                            ] = margin;
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
                                } else {
                                    ++special_negative_groups;
                                    if (special_negative_groups == 1) {
                                        std::cout
                                            << "FIRST_NEGATIVE_SPECIAL"
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

                            Integer lower_left = 0;
                            Integer upper_left = 0;
                            Integer upper_right = 0;
                            bool has_lower_left = false;
                            bool has_upper_left = false;
                            bool has_upper_right = false;
                            for (const int y :
                                 graph[static_cast<std::size_t>(Q)]) {
                                if (
                                    y == 0
                                    || y == Q
                                    || std::abs(y - (K - Q)) < distance
                                ) {
                                    continue;
                                }
                                const Integer& margin = channel_margins[
                                    static_cast<std::size_t>(y)
                                ];
                                if (y < Q) {
                                    lower_left += margin;
                                    has_lower_left = true;
                                } else if (y <= K - Q - distance) {
                                    upper_left += margin;
                                    has_upper_left = true;
                                } else {
                                    upper_right += margin;
                                    has_upper_right = true;
                                }
                            }
                            const Integer pieces[6] = {
                                lower_left,
                                upper_left,
                                upper_right,
                                lower_left + upper_left,
                                lower_left + upper_right,
                                upper_left + upper_right
                            };
                            const bool present[6] = {
                                has_lower_left,
                                has_upper_left,
                                has_upper_right,
                                has_lower_left || has_upper_left,
                                has_lower_left || has_upper_right,
                                has_upper_left || has_upper_right
                            };
                            for (std::size_t piece = 0;
                                 piece < 6;
                                 ++piece) {
                                if (!present[piece]) {
                                    continue;
                                }
                                ++interval_piece_tests[piece];
                                if (pieces[piece] < 0) {
                                    ++negative_interval_piece_tests[piece];
                                    if (
                                        negative_interval_piece_tests[piece]
                                        == 1
                                    ) {
                                        std::cout
                                            << "FIRST_NEGATIVE_INTERVAL_PIECE"
                                            << " piece="
                                            << interval_piece_names[piece]
                                            << " K=" << K
                                            << " Q=" << Q
                                            << " j=" << j
                                            << " t=" << t
                                            << " distance=" << distance
                                            << " value=" << pieces[piece]
                                            << '\n';
                                    }
                                }
                            }
                            if (
                                has_upper_left
                                && has_upper_right
                                && upper_left < 0
                                && upper_right < 0
                            ) {
                                ++simultaneous_negative_upper_pieces;
                                if (simultaneous_negative_upper_pieces == 1) {
                                    std::cout
                                        << "FIRST_SIMULTANEOUS_NEGATIVE_UPPER"
                                        << " K=" << K
                                        << " Q=" << Q
                                        << " j=" << j
                                        << " t=" << t
                                        << " distance=" << distance
                                        << " lower_left=" << lower_left
                                        << " upper_left=" << upper_left
                                        << " upper_right=" << upper_right
                                        << " total="
                                        << lower_left
                                            + upper_left
                                            + upper_right
                                        << '\n';
                                }
                                if (distance >= 1) {
                                    ++strict_simultaneous_negative_upper_pieces;
                                    if (
                                        strict_simultaneous_negative_upper_pieces
                                        == 1
                                    ) {
                                        std::cout
                                            << "FIRST_STRICT_SIMULTANEOUS_NEGATIVE_UPPER"
                                            << " K=" << K
                                            << " Q=" << Q
                                            << " j=" << j
                                            << " t=" << t
                                            << " distance=" << distance
                                            << " lower_left=" << lower_left
                                            << " upper_left=" << upper_left
                                            << " upper_right=" << upper_right
                                            << " total="
                                            << lower_left
                                                + upper_left
                                                + upper_right
                                            << '\n';
                                    }
                                }
                            }
                            if (
                                present[3]
                                && present[5]
                                && pieces[3] < 0
                                && pieces[5] < 0
                            ) {
                                ++simultaneous_negative_ray_pairs;
                                if (simultaneous_negative_ray_pairs == 1) {
                                    std::cout
                                        << "FIRST_SIMULTANEOUS_NEGATIVE_RAYS"
                                        << " K=" << K
                                        << " Q=" << Q
                                        << " j=" << j
                                        << " t=" << t
                                        << " distance=" << distance
                                        << " left_ray=" << pieces[3]
                                        << " upper_pair=" << pieces[5]
                                    << '\n';
                                }
                            }
                            if (
                                present[3]
                                && pieces[3] < 0
                                && (
                                    K != last_negative_left_ray_K
                                    || Q != last_negative_left_ray_Q
                                    || distance
                                        != last_negative_left_ray_distance
                                )
                            ) {
                                last_negative_left_ray_K = K;
                                last_negative_left_ray_Q = Q;
                                last_negative_left_ray_distance = distance;
                                std::cout
                                    << "NEGATIVE_LEFT_RAY_FAMILY"
                                    << " K=" << K
                                    << " Q=" << Q
                                    << " j=" << j
                                    << " t=" << t
                                    << " distance=" << distance
                                    << " left_ray=" << pieces[3]
                                    << " upper_right=" << upper_right
                                    << " total="
                                    << lower_left
                                        + upper_left
                                        + upper_right
                                    << '\n';
                            }
                            if (
                                distance >= 1
                                && present[3]
                                && pieces[3] < 0
                            ) {
                                ++strict_negative_left_rays;
                                if (strict_negative_left_rays == 1) {
                                    std::cout
                                        << "FIRST_STRICT_NEGATIVE_LEFT_RAY"
                                        << " K=" << K
                                        << " Q=" << Q
                                        << " j=" << j
                                        << " t=" << t
                                        << " distance=" << distance
                                        << " value=" << pieces[3]
                                    << '\n';
                                }
                            }
                            if (distance >= 1) {
                                ++strict_deficit_reserve_tests;
                                Integer deficit_reserve = lower_left;
                                if (upper_left < 0) {
                                    deficit_reserve += upper_left;
                                }
                                if (upper_right < 0) {
                                    deficit_reserve += upper_right;
                                }
                                if (deficit_reserve < 0) {
                                    ++negative_strict_deficit_reserves;
                                    if (
                                        negative_strict_deficit_reserves
                                        == 1
                                    ) {
                                        std::cout
                                            << "FIRST_NEGATIVE_STRICT_DEFICIT_RESERVE"
                                            << " K=" << K
                                            << " Q=" << Q
                                            << " j=" << j
                                            << " t=" << t
                                            << " distance=" << distance
                                            << " lower_left=" << lower_left
                                            << " upper_left=" << upper_left
                                            << " upper_right=" << upper_right
                                            << " reserve=" << deficit_reserve
                                            << '\n';
                                    }
                                }
                                for (int radius = 1;
                                     radius < Q;
                                     ++radius) {
                                    const int lower = Q - radius;
                                    const int upper_left_label = Q + radius;
                                    const int upper_right_label =
                                        K - Q + radius;
                                    const bool use_upper_left =
                                        upper_left_label <= 2 * Q
                                        && std::abs(
                                            upper_left_label - (K - Q)
                                        ) >= distance
                                        && upper_left_label
                                            <= K - Q - distance;
                                    const bool use_upper_right =
                                        upper_right_label <= 2 * Q
                                        && std::abs(
                                            upper_right_label - (K - Q)
                                        ) >= distance
                                        && upper_right_label
                                            > K - Q - distance;
                                    if (
                                        !use_upper_left
                                        && !use_upper_right
                                    ) {
                                        continue;
                                    }
                                    if (
                                        std::abs(lower - (K - Q))
                                        < distance
                                    ) {
                                        throw std::runtime_error(
                                            "local lower partner is absent"
                                        );
                                    }
                                    Integer local_reserve =
                                        channel_margins[
                                            static_cast<std::size_t>(lower)
                                        ];
                                    if (
                                        use_upper_left
                                        && channel_margins[
                                            static_cast<std::size_t>(
                                                upper_left_label
                                            )
                                        ] < 0
                                    ) {
                                        local_reserve += channel_margins[
                                            static_cast<std::size_t>(
                                                upper_left_label
                                            )
                                        ];
                                    }
                                    if (
                                        use_upper_right
                                        && channel_margins[
                                            static_cast<std::size_t>(
                                                upper_right_label
                                            )
                                        ] < 0
                                    ) {
                                        local_reserve += channel_margins[
                                            static_cast<std::size_t>(
                                                upper_right_label
                                            )
                                        ];
                                    }
                                    ++local_triple_reserve_tests;
                                    if (local_reserve < 0) {
                                        ++negative_local_triple_reserves;
                                        if (
                                            negative_local_triple_reserves
                                            == 1
                                        ) {
                                            std::cout
                                                << "FIRST_NEGATIVE_LOCAL_TRIPLE_RESERVE"
                                                << " K=" << K
                                                << " Q=" << Q
                                                << " j=" << j
                                                << " t=" << t
                                                << " distance=" << distance
                                                << " radius=" << radius
                                                << " lower="
                                                << channel_margins[
                                                    static_cast<std::size_t>(
                                                        lower
                                                    )
                                                ]
                                                << " upper_left="
                                                << (
                                                    use_upper_left
                                                        ? channel_margins[
                                                            static_cast<
                                                                std::size_t
                                                            >(
                                                                upper_left_label
                                                            )
                                                        ] : Integer(0)
                                                )
                                                << " upper_right="
                                                << (
                                                    use_upper_right
                                                        ? channel_margins[
                                                            static_cast<
                                                                std::size_t
                                                            >(
                                                                upper_right_label
                                                            )
                                                        ] : Integer(0)
                                                )
                                                << " reserve="
                                                << local_reserve
                                                << '\n';
                                        }
                                    }
                                }
                                std::vector<Integer> radius_reserves(
                                    static_cast<std::size_t>(Q + 1)
                                );
                                std::vector<bool> radius_present(
                                    static_cast<std::size_t>(Q + 1),
                                    false
                                );
                                for (int radius = 1;
                                     radius <= Q;
                                     ++radius) {
                                    const int lower = Q - radius;
                                    const bool lower_present =
                                        lower >= 1
                                        &&
                                        std::abs(lower - (K - Q))
                                        >= distance;
                                    Integer& reserve = radius_reserves[
                                        static_cast<std::size_t>(radius)
                                    ];
                                    if (lower_present) {
                                        radius_present[
                                            static_cast<std::size_t>(radius)
                                        ] = true;
                                        reserve = channel_margins[
                                            static_cast<std::size_t>(lower)
                                        ];
                                    }
                                    const int upper_left_label = Q + radius;
                                    if (
                                        upper_left_label <= 2 * Q
                                        && upper_left_label
                                            <= K - Q - distance
                                        && channel_margins[
                                            static_cast<std::size_t>(
                                                upper_left_label
                                            )
                                        ] < 0
                                    ) {
                                        reserve += channel_margins[
                                            static_cast<std::size_t>(
                                                upper_left_label
                                            )
                                        ];
                                        radius_present[
                                            static_cast<std::size_t>(radius)
                                        ] = true;
                                    }
                                    const int upper_right_label =
                                        K - Q + radius;
                                    if (
                                        upper_right_label <= 2 * Q
                                        && radius >= distance
                                        && channel_margins[
                                            static_cast<std::size_t>(
                                                upper_right_label
                                            )
                                        ] < 0
                                    ) {
                                        reserve += channel_margins[
                                            static_cast<std::size_t>(
                                                upper_right_label
                                            )
                                        ];
                                        radius_present[
                                            static_cast<std::size_t>(radius)
                                        ] = true;
                                    }
                                }
                                Integer nested_reserve = 0;
                                for (int radius = Q;
                                     radius >= 1;
                                     --radius) {
                                    if (!radius_present[
                                            static_cast<std::size_t>(radius)
                                        ]) {
                                        continue;
                                    }
                                    nested_reserve += radius_reserves[
                                        static_cast<std::size_t>(radius)
                                    ];
                                    if (radius == Q) {
                                        continue;
                                    }
                                    ++nested_radius_reserve_tests;
                                    if (nested_reserve < 0) {
                                        ++negative_nested_radius_reserves;
                                        if (
                                            negative_nested_radius_reserves
                                            == 1
                                        ) {
                                            std::cout
                                                << "FIRST_NEGATIVE_NESTED_RADIUS_RESERVE"
                                                << " K=" << K
                                                << " Q=" << Q
                                                << " j=" << j
                                                << " t=" << t
                                                << " distance=" << distance
                                                << " radius=" << radius
                                                << " value="
                                                << nested_reserve
                                                << '\n';
                                        }
                                    }
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
                std::vector<Series> t4_differences = full_t4_kernels;
                for (int order = 0;
                     order <= 8 && !t4_differences.empty();
                     ++order) {
                    for (int x = 0; x <= K; ++x) {
                        ++t4_base_difference_coordinates;
                        const Integer& value = t4_differences[0][
                            static_cast<std::size_t>(x)
                        ];
                        if (value < 0) {
                            ++negative_t4_base_difference_coordinates;
                            if (
                                negative_t4_base_difference_coordinates
                                == 1
                            ) {
                                std::cout
                                    << "FIRST_NEGATIVE_T4_BASE_DIFFERENCE"
                                    << " K=" << K
                                    << " Q=" << Q
                                    << " order=" << order
                                    << " x=" << x
                                    << " value=" << value
                                    << '\n';
                            }
                        }
                    }
                    std::vector<Series> next_differences;
                    for (std::size_t index = 0;
                         index + 1 < t4_differences.size();
                         ++index) {
                        Series advanced = t4_differences[index];
                        for (int step = 0; step < 2; ++step) {
                            Series next(
                                static_cast<std::size_t>(K + 1)
                            );
                            for (int x = 0; x <= K; ++x) {
                                for (const int y :
                                     graph[static_cast<std::size_t>(x)]) {
                                    next[static_cast<std::size_t>(y)] +=
                                        advanced[
                                            static_cast<std::size_t>(x)
                                        ];
                                }
                            }
                            advanced = std::move(next);
                        }
                        Series difference =
                            t4_differences[index + 1];
                        for (int x = 0; x <= K; ++x) {
                            difference[static_cast<std::size_t>(x)] -=
                                advanced[static_cast<std::size_t>(x)];
                        }
                        next_differences.push_back(
                            std::move(difference)
                        );
                    }
                    t4_differences = std::move(next_differences);
                }
                if (full_t4_kernels.size() >= 10) {
                    std::vector<Series> ninth = full_t4_kernels;
                    for (int order = 0; order < 9; ++order) {
                        std::vector<Series> next_differences;
                        for (std::size_t index = 0;
                             index + 1 < ninth.size();
                             ++index) {
                            Series advanced = ninth[index];
                            for (int step = 0; step < 2; ++step) {
                                Series next(
                                    static_cast<std::size_t>(K + 1)
                                );
                                for (int x = 0; x <= K; ++x) {
                                    for (const int y :
                                         graph[
                                             static_cast<std::size_t>(x)
                                         ]) {
                                        next[
                                            static_cast<std::size_t>(y)
                                        ] += advanced[
                                            static_cast<std::size_t>(x)
                                        ];
                                    }
                                }
                                advanced = std::move(next);
                            }
                            Series difference = ninth[index + 1];
                            for (int x = 0; x <= K; ++x) {
                                difference[
                                    static_cast<std::size_t>(x)
                                ] -= advanced[
                                    static_cast<std::size_t>(x)
                                ];
                            }
                            next_differences.push_back(
                                std::move(difference)
                            );
                        }
                        ninth = std::move(next_differences);
                    }
                    for (const Series& difference : ninth) {
                        for (const Integer& value : difference) {
                            if (value != 0) {
                                ++nonzero_t4_ninth_difference_coordinates;
                            }
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
            << " special_negative_groups=" << special_negative_groups
            << " special_positive_groups=" << special_positive_groups
            << " outside_in_prefixes=" << outside_in_prefixes
            << " negative_outside_in_prefixes="
            << negative_outside_in_prefixes
            << " renewal_current_coefficients="
            << renewal_current_coefficients
            << " negative_renewal_current_coefficients="
            << negative_renewal_current_coefficients
            << " twice_return_current_coefficients="
            << twice_return_current_coefficients
            << " negative_twice_return_current_coefficients="
            << negative_twice_return_current_coefficients
            << " green_minor_identities=" << green_minor_identities
            << " outside_in_core_identities="
            << outside_in_core_identities
            << " outside_in_core_coefficients="
            << outside_in_core_coefficients
            << " negative_outside_in_core_coefficients="
            << negative_outside_in_core_coefficients
            << " negative_higher_q_core_coefficients="
            << negative_higher_q_core_coefficients
            << " negative_higher_q_twice_return_coefficients="
            << negative_higher_q_twice_return_coefficients
            << " interval_piece_tests="
            << interval_piece_tests[0] << ','
            << interval_piece_tests[1] << ','
            << interval_piece_tests[2] << ','
            << interval_piece_tests[3] << ','
            << interval_piece_tests[4] << ','
            << interval_piece_tests[5]
            << " negative_interval_piece_tests="
            << negative_interval_piece_tests[0] << ','
            << negative_interval_piece_tests[1] << ','
            << negative_interval_piece_tests[2] << ','
            << negative_interval_piece_tests[3] << ','
            << negative_interval_piece_tests[4] << ','
            << negative_interval_piece_tests[5]
            << " simultaneous_negative_upper_pieces="
            << simultaneous_negative_upper_pieces
            << " simultaneous_negative_ray_pairs="
            << simultaneous_negative_ray_pairs
            << " strict_simultaneous_negative_upper_pieces="
            << strict_simultaneous_negative_upper_pieces
            << " strict_negative_left_rays="
            << strict_negative_left_rays
            << " strict_deficit_reserve_tests="
            << strict_deficit_reserve_tests
            << " negative_strict_deficit_reserves="
            << negative_strict_deficit_reserves
            << " local_triple_reserve_tests="
            << local_triple_reserve_tests
            << " negative_local_triple_reserves="
            << negative_local_triple_reserves
            << " nested_radius_reserve_tests="
            << nested_radius_reserve_tests
            << " negative_nested_radius_reserves="
            << negative_nested_radius_reserves
            << " odd_boundary_current_coefficients="
            << odd_boundary_current_coefficients
            << " negative_odd_boundary_current_coefficients="
            << negative_odd_boundary_current_coefficients
            << " negative_higher_q_odd_boundary_coefficients="
            << negative_higher_q_odd_boundary_coefficients
            << " t4_kernel_coordinates=" << t4_kernel_coordinates
            << " negative_t4_kernel_coordinates="
            << negative_t4_kernel_coordinates
            << " negative_higher_q_t4_kernel_coordinates="
            << negative_higher_q_t4_kernel_coordinates
            << " full_t4_kernel_coordinates="
            << full_t4_kernel_coordinates
            << " negative_full_t4_kernel_coordinates="
            << negative_full_t4_kernel_coordinates
            << " t4_induction_coordinates="
            << t4_induction_coordinates
            << " negative_t4_induction_coordinates="
            << negative_t4_induction_coordinates
            << " t4_base_difference_coordinates="
            << t4_base_difference_coordinates
            << " negative_t4_base_difference_coordinates="
            << negative_t4_base_difference_coordinates
            << " nonzero_t4_ninth_difference_coordinates="
            << nonzero_t4_ninth_difference_coordinates
            << " t4_c4_coordinates=" << t4_c4_coordinates
            << " negative_t4_c4_coordinates="
            << negative_t4_c4_coordinates
            << " t4_grouped_block_coordinates="
            << t4_grouped_block_coordinates
            << " negative_t4_grouped_block_coordinates="
            << negative_t4_grouped_block_coordinates
            << " negative_t4_low_component_coordinates="
            << negative_t4_low_component_coordinates
            << " negative_t4_weak_group_coordinates="
            << negative_t4_weak_group_coordinates
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
