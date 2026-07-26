#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Minus = std::array<int, 2>;
using Plus = std::array<int, 5>;

template <class F>
void outputs(const int k, const int a, const int b, F f) {
    const int upper = std::min(a + b, 2 * k - a - b);
    for (int c = std::abs(a - b); c <= upper; c += 2) {
        f(c);
    }
}

bool contains_output(const int k, const int a, const int b, const int target) {
    bool found = false;
    outputs(k, a, b, [&](const int c) {
        if (c == target) {
            found = true;
        }
    });
    return found;
}

std::vector<int> fusion_outputs(const int k, const int a, const int b) {
    std::vector<int> result;
    outputs(k, a, b, [&](const int c) {
        result.push_back(c);
    });
    return result;
}

int intersection_size(const std::vector<int>& lhs, const std::vector<int>& rhs) {
    int count = 0;
    for (const int x : lhs) {
        count += std::binary_search(rhs.begin(), rhs.end(), x) ? 1 : 0;
    }
    return count;
}

int interval_pair_profile(
    const int k,
    const std::vector<int>& lhs,
    const std::vector<int>& rhs,
    const int target
) {
    int count = 0;
    for (const int x : lhs) {
        for (const int y : rhs) {
            count += contains_output(k, x, y, target) ? 1 : 0;
        }
    }
    return count;
}

std::int64_t mult(const int k, const std::vector<int>& labels) {
    std::vector<std::int64_t> cur(static_cast<std::size_t>(k + 1));
    std::vector<std::int64_t> next(static_cast<std::size_t>(k + 1));
    cur[0] = 1;
    for (const int p : labels) {
        std::fill(next.begin(), next.end(), 0);
        for (int x = 0; x <= k; ++x) {
            const auto value = cur[static_cast<std::size_t>(x)];
            if (value == 0) {
                continue;
            }
            outputs(k, x, p, [&](const int y) {
                next[static_cast<std::size_t>(y)] += value;
            });
        }
        cur.swap(next);
    }
    return cur[0];
}

std::vector<std::int64_t> decomposition(
    const int k,
    const std::vector<int>& labels
) {
    std::vector<std::int64_t> cur(static_cast<std::size_t>(k + 1));
    std::vector<std::int64_t> next(static_cast<std::size_t>(k + 1));
    cur[0] = 1;
    for (const int p : labels) {
        std::fill(next.begin(), next.end(), 0);
        for (int x = 0; x <= k; ++x) {
            const auto value = cur[static_cast<std::size_t>(x)];
            if (value == 0) {
                continue;
            }
            outputs(k, x, p, [&](const int y) {
                next[static_cast<std::size_t>(y)] += value;
            });
        }
        cur.swap(next);
    }
    return cur;
}

bool disjoint(const Minus& minus, const Plus& plus) {
    for (const int q : minus) {
        for (const int p : plus) {
            if (q == p) {
                return false;
            }
        }
    }
    return true;
}

int cap(const int k, const int x) {
    return std::min(x, k - x) + 1;
}

struct Stats {
    std::int64_t n = 0;
    std::int64_t peq = 0;
    std::int64_t peq_minus = 0;
    std::int64_t peq_plus = 0;
    std::int64_t u1 = 0;
    std::int64_t u2 = 0;
    std::int64_t t = 0;
    int max_cut = 0;
    int active = 0;
};

Stats statistics(const int k, const Minus& minus, const Plus& plus) {
    Stats s;
    s.n = mult(k, {
        minus[0], minus[1], plus[0], plus[1], plus[2], plus[3], plus[4]
    });
    if (minus[0] == minus[1]) {
        s.peq_minus =
            mult(k, {plus[0], plus[1], plus[2], plus[3], plus[4]});
        s.peq += s.peq_minus;
    }
    for (std::size_t i = 0; i < plus.size(); ++i) {
        std::vector<int> rest;
        for (std::size_t j = 0; j < plus.size(); ++j) {
            if (j != i) {
                rest.push_back(plus[j]);
            }
        }
        s.u1 += mult(k, {minus[0], minus[1], plus[i]}) * mult(k, rest);
    }
    for (std::size_t i = 0; i < plus.size(); ++i) {
        for (std::size_t j = i + 1; j < plus.size(); ++j) {
            std::vector<int> rest;
            std::vector<int> positive_rest{minus[0], minus[1]};
            for (std::size_t h = 0; h < plus.size(); ++h) {
                if (h == i || h == j) {
                    continue;
                }
                rest.push_back(plus[h]);
                positive_rest.push_back(plus[h]);
            }
            if (plus[i] == plus[j]) {
                const std::int64_t contribution = mult(k, positive_rest);
                s.peq_plus += contribution;
                s.peq += contribution;
            }
            s.u2 += mult(k, {minus[0], minus[1], plus[i], plus[j]})
                * mult(k, rest);
            for (std::size_t o = 0; o < minus.size(); ++o) {
                const auto triple = mult(k, {minus[o], plus[i], plus[j]});
                if (triple == 0) {
                    continue;
                }
                auto block = rest;
                block.push_back(minus[1U - o]);
                const auto value = triple * mult(k, block);
                if (value == 0) {
                    continue;
                }
                ++s.active;
                s.t += value;
                s.max_cut = std::max<int>(s.max_cut, static_cast<int>(value));
            }
        }
    }
    return s;
}

std::int64_t first_two_equal_plus_reservoir(
    const int k,
    const Minus& minus,
    const Plus& plus
) {
    const std::vector<int> minus_outputs =
        fusion_outputs(k, minus[0], minus[1]);
    std::int64_t result = 0;
    for (std::size_t i = 0; i < plus.size(); ++i) {
        for (std::size_t j = i + 1; j < plus.size(); ++j) {
            if (plus[i] != plus[j]) {
                continue;
            }
            std::vector<int> rest;
            for (std::size_t h = 0; h < plus.size(); ++h) {
                if (h != i && h != j) {
                    rest.push_back(plus[h]);
                }
            }
            const std::vector<std::int64_t> profile =
                decomposition(k, rest);
            const std::size_t count =
                std::min<std::size_t>(2U, minus_outputs.size());
            for (std::size_t output = 0; output < count; ++output) {
                result += profile[
                    static_cast<std::size_t>(minus_outputs[output])
                ];
            }
        }
    }
    return result;
}

struct GraphStats {
    std::array<int, 2> edge_count{};
    std::array<int, 2> weight_two{};
    std::array<int, 2> mask{};
    std::array<int, 2> weight_two_mask{};
};

int incident_vertex_mask(const int edge_mask) {
    int vertices = 0;
    int bit = 0;
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j, ++bit) {
            if ((edge_mask & (1 << bit)) != 0) {
                vertices |= 1 << i;
                vertices |= 1 << j;
            }
        }
    }
    return vertices;
}

GraphStats graph_stats(const int k, const Minus& minus, const Plus& plus) {
    GraphStats g;
    int bit = 0;
    for (std::size_t i = 0; i < plus.size(); ++i) {
        for (std::size_t j = i + 1; j < plus.size(); ++j, ++bit) {
            std::vector<int> rest;
            for (std::size_t h = 0; h < plus.size(); ++h) {
                if (h != i && h != j) {
                    rest.push_back(plus[h]);
                }
            }
            for (std::size_t o = 0; o < minus.size(); ++o) {
                const int r = minus[o];
                const int s = minus[1U - o];
                if (!contains_output(k, plus[i], plus[j], r)) {
                    continue;
                }
                const auto rank = mult(k, {s, rest[0], rest[1], rest[2]});
                if (rank == 0 || rank > 2) {
                    continue;
                }
                ++g.edge_count[o];
                if (rank == 2) {
                    ++g.weight_two[o];
                    g.weight_two_mask[o] |= 1 << bit;
                }
                g.mask[o] |= 1 << bit;
            }
        }
    }
    return g;
}

std::int64_t best_partial_cut_bound(
    const int k,
    const Minus& minus,
    const Plus& plus,
    const int max_output
) {
    std::int64_t best = 0;
    for (std::size_t i = 0; i < plus.size(); ++i) {
        for (std::size_t j = i + 1; j < plus.size(); ++j) {
            std::vector<int> rest;
            for (std::size_t h = 0; h < plus.size(); ++h) {
                if (h != i && h != j) {
                    rest.push_back(plus[h]);
                }
            }
            for (std::size_t o = 0; o < minus.size(); ++o) {
                const int r = minus[o];
                if (!contains_output(k, plus[i], plus[j], r)) {
                    continue;
                }
                auto block = rest;
                block.push_back(minus[1U - o]);
                const auto active_profile =
                    decomposition(k, {r, plus[i], plus[j]});
                const auto complement_profile = decomposition(k, block);
                if (complement_profile[0] == 0
                    || complement_profile[0] > 2) {
                    continue;
                }
                std::int64_t bound = 0;
                const int last = std::min(k, max_output);
                for (int output = 0; output <= last; ++output) {
                    bound += active_profile[static_cast<std::size_t>(output)]
                        * complement_profile[
                            static_cast<std::size_t>(output)
                        ];
                }
                best = std::max(best, bound);
            }
        }
    }
    return best;
}

std::int64_t partial_cut_bound(
    const int k,
    const Minus& minus,
    const Plus& plus,
    const int orientation,
    const int first,
    const int second,
    const int max_output
) {
    if (!contains_output(
            k,
            plus[static_cast<std::size_t>(first)],
            plus[static_cast<std::size_t>(second)],
            minus[static_cast<std::size_t>(orientation)]
        )) {
        return 0;
    }
    std::vector<int> rest;
    for (int index = 0; index < 5; ++index) {
        if (index != first && index != second) {
            rest.push_back(plus[static_cast<std::size_t>(index)]);
        }
    }
    const auto active_profile = decomposition(
        k,
        {
            minus[static_cast<std::size_t>(orientation)],
            plus[static_cast<std::size_t>(first)],
            plus[static_cast<std::size_t>(second)]
        }
    );
    rest.push_back(minus[static_cast<std::size_t>(1 - orientation)]);
    const auto complement_profile = decomposition(k, rest);
    std::int64_t result = 0;
    const int last = std::min(k, max_output);
    for (int output = 0; output <= last; ++output) {
        result += active_profile[static_cast<std::size_t>(output)]
            * complement_profile[static_cast<std::size_t>(output)];
    }
    return result;
}

std::int64_t worst_partial_cut_bound(
    const int k,
    const Minus& minus,
    const Plus& plus,
    const int max_output,
    const int required_rank
) {
    bool set = false;
    std::int64_t worst = 0;
    for (std::size_t i = 0; i < plus.size(); ++i) {
        for (std::size_t j = i + 1; j < plus.size(); ++j) {
            std::vector<int> rest;
            for (std::size_t h = 0; h < plus.size(); ++h) {
                if (h != i && h != j) {
                    rest.push_back(plus[h]);
                }
            }
            for (std::size_t o = 0; o < minus.size(); ++o) {
                const int r = minus[o];
                if (!contains_output(k, plus[i], plus[j], r)) {
                    continue;
                }
                auto block = rest;
                block.push_back(minus[1U - o]);
                const auto active_profile =
                    decomposition(k, {r, plus[i], plus[j]});
                const auto complement_profile = decomposition(k, block);
                if (complement_profile[0] == 0
                    || complement_profile[0] > 2) {
                    continue;
                }
                if (required_rank > 0
                    && complement_profile[0] != required_rank) {
                    continue;
                }
                std::int64_t bound = 0;
                const int last = std::min(k, max_output);
                for (int output = 0; output <= last; ++output) {
                    bound += active_profile[static_cast<std::size_t>(output)]
                        * complement_profile[
                            static_cast<std::size_t>(output)
                        ];
                }
                if (!set || bound < worst) {
                    set = true;
                    worst = bound;
                }
            }
        }
    }
    return worst;
}

std::int64_t sum_partial_cut_bounds(
    const int k,
    const Minus& minus,
    const Plus& plus,
    const int max_output
) {
    std::int64_t sum = 0;
    for (std::size_t i = 0; i < plus.size(); ++i) {
        for (std::size_t j = i + 1; j < plus.size(); ++j) {
            std::vector<int> rest;
            for (std::size_t h = 0; h < plus.size(); ++h) {
                if (h != i && h != j) {
                    rest.push_back(plus[h]);
                }
            }
            for (std::size_t o = 0; o < minus.size(); ++o) {
                const int r = minus[o];
                if (!contains_output(k, plus[i], plus[j], r)) {
                    continue;
                }
                auto block = rest;
                block.push_back(minus[1U - o]);
                const auto active_profile =
                    decomposition(k, {r, plus[i], plus[j]});
                const auto complement_profile = decomposition(k, block);
                if (complement_profile[0] == 0
                    || complement_profile[0] > 2) {
                    continue;
                }
                const int last = std::min(k, max_output);
                for (int output = 0; output <= last; ++output) {
                    sum += active_profile[static_cast<std::size_t>(output)]
                        * complement_profile[
                            static_cast<std::size_t>(output)
                        ];
                }
            }
        }
    }
    return sum;
}

std::int64_t best_max_weighted_incidence_bound(
    const int k,
    const Minus& minus,
    const Plus& plus,
    const int max_output
) {
    struct Cut {
        int i = 0;
        int j = 0;
        int weight = 0;
        std::int64_t bound = 0;
    };
    std::vector<Cut> cuts;
    std::array<int, 5> weighted_degree{};
    for (std::size_t i = 0; i < plus.size(); ++i) {
        for (std::size_t j = i + 1; j < plus.size(); ++j) {
            std::vector<int> rest;
            for (std::size_t h = 0; h < plus.size(); ++h) {
                if (h != i && h != j) {
                    rest.push_back(plus[h]);
                }
            }
            for (std::size_t o = 0; o < minus.size(); ++o) {
                const int r = minus[o];
                if (!contains_output(k, plus[i], plus[j], r)) {
                    continue;
                }
                auto block = rest;
                block.push_back(minus[1U - o]);
                const auto active_profile =
                    decomposition(k, {r, plus[i], plus[j]});
                const auto complement_profile = decomposition(k, block);
                const auto rank = complement_profile[0];
                if (rank == 0 || rank > 2) {
                    continue;
                }
                std::int64_t bound = 0;
                const int last = std::min(k, max_output);
                for (int output = 0; output <= last; ++output) {
                    bound += active_profile[static_cast<std::size_t>(output)]
                        * complement_profile[
                            static_cast<std::size_t>(output)
                        ];
                }
                const int weight = static_cast<int>(rank);
                cuts.push_back({
                    static_cast<int>(i), static_cast<int>(j), weight, bound
                });
                weighted_degree[i] += weight;
                weighted_degree[j] += weight;
            }
        }
    }
    int best_score = -1;
    std::int64_t best_bound = 0;
    for (const Cut& cut : cuts) {
        const int score =
            weighted_degree[static_cast<std::size_t>(cut.i)]
            + weighted_degree[static_cast<std::size_t>(cut.j)]
            - cut.weight;
        if (score > best_score) {
            best_score = score;
            best_bound = cut.bound;
        } else if (score == best_score) {
            best_bound = std::max(best_bound, cut.bound);
        }
    }
    return best_bound;
}

std::int64_t best_light_color_bound(
    const int k,
    const Minus& minus,
    const Plus& plus,
    const int max_output
) {
    std::array<int, 2> color_weight{};
    std::array<std::int64_t, 2> best{};
    for (std::size_t i = 0; i < plus.size(); ++i) {
        for (std::size_t j = i + 1; j < plus.size(); ++j) {
            std::vector<int> rest;
            for (std::size_t h = 0; h < plus.size(); ++h) {
                if (h != i && h != j) {
                    rest.push_back(plus[h]);
                }
            }
            for (std::size_t o = 0; o < minus.size(); ++o) {
                const int r = minus[o];
                if (!contains_output(k, plus[i], plus[j], r)) {
                    continue;
                }
                auto block = rest;
                block.push_back(minus[1U - o]);
                const auto active_profile =
                    decomposition(k, {r, plus[i], plus[j]});
                const auto complement_profile = decomposition(k, block);
                const auto rank = complement_profile[0];
                if (rank == 0 || rank > 2) {
                    continue;
                }
                color_weight[o] += static_cast<int>(rank);
                std::int64_t bound = 0;
                const int last = std::min(k, max_output);
                for (int output = 0; output <= last; ++output) {
                    bound += active_profile[static_cast<std::size_t>(output)]
                        * complement_profile[
                            static_cast<std::size_t>(output)
                        ];
                }
                best[o] = std::max(best[o], bound);
            }
        }
    }
    if (color_weight[0] < color_weight[1]) {
        return best[0];
    }
    if (color_weight[1] < color_weight[0]) {
        return best[1];
    }
    return std::max(best[0], best[1]);
}

std::pair<std::int64_t, std::int64_t>
rank_two_incidence_one_bounds(
    const int k,
    const Minus& minus,
    const Plus& plus,
    const int max_output
) {
    bool set = false;
    std::int64_t worst = 0;
    std::int64_t best = 0;
    for (std::size_t i = 0; i < plus.size(); ++i) {
        for (std::size_t j = i + 1; j < plus.size(); ++j) {
            const int incidence =
                (cap(k, plus[i]) == 2 ? 1 : 0)
                + (cap(k, plus[j]) == 2 ? 1 : 0);
            if (incidence != 1) {
                continue;
            }
            std::vector<int> rest;
            for (std::size_t h = 0; h < plus.size(); ++h) {
                if (h != i && h != j) {
                    rest.push_back(plus[h]);
                }
            }
            for (std::size_t o = 0; o < minus.size(); ++o) {
                const int r = minus[o];
                if (!contains_output(k, plus[i], plus[j], r)) {
                    continue;
                }
                auto block = rest;
                block.push_back(minus[1U - o]);
                const auto active_profile =
                    decomposition(k, {r, plus[i], plus[j]});
                const auto complement_profile = decomposition(k, block);
                if (complement_profile[0] != 2) {
                    continue;
                }
                std::int64_t bound = 0;
                const int last = std::min(k, max_output);
                for (int output = 0; output <= last; ++output) {
                    bound += active_profile[
                        static_cast<std::size_t>(output)
                    ] * complement_profile[
                        static_cast<std::size_t>(output)
                    ];
                }
                if (!set) {
                    set = true;
                    worst = bound;
                    best = bound;
                } else {
                    worst = std::min(worst, bound);
                    best = std::max(best, bound);
                }
            }
        }
    }
    return {worst, best};
}

void print_case(
    const char* tag,
    const int k,
    const Minus& minus,
    const Plus& plus,
    const Stats& s,
    const GraphStats& g
) {
    std::cout << tag << " k=" << k << " minus=[" << minus[0] << ','
              << minus[1] << "] plus=[";
    for (std::size_t i = 0; i < plus.size(); ++i) {
        std::cout << (i == 0 ? "" : ",") << plus[i];
    }
    std::cout << "] c=(" << g.edge_count[0] << ',' << g.edge_count[1]
              << ") h=(" << g.weight_two[0] << ',' << g.weight_two[1]
              << ") mask=(" << g.mask[0] << ',' << g.mask[1]
              << ") hmask=(" << g.weight_two_mask[0] << ','
              << g.weight_two_mask[1] << ")"
              << " d=" << s.max_cut << " N=" << s.n
              << " P=" << s.peq << "(M" << s.peq_minus
              << ",P" << s.peq_plus << ')'
              << " U1=" << s.u1 << " U2=" << s.u2 << " T=" << s.t
              << " raw=" << s.n + s.peq - s.t
              << " paid=" << s.n + s.peq + s.u1 + s.u2 - s.t << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 6 && std::string(argv[1]) == "--random-direct") {
        const bool shallow_mode = std::string(argv[2]) == "shallow";
        if (!shallow_mode && std::string(argv[2]) != "deep") {
            std::cerr << "MODE must be shallow or deep\n";
            return EXIT_FAILURE;
        }
        const std::uint64_t seed = static_cast<std::uint64_t>(
            std::strtoull(argv[3], nullptr, 10)
        );
        const std::int64_t trials = std::atoll(argv[4]);
        const int max_k = std::atoi(argv[5]);
        std::mt19937_64 generator(seed);
        std::int64_t tested = 0;
        bool minimum_set = false;
        std::int64_t minimum = 0;
        Minus minimum_minus{};
        Plus minimum_plus{};
        int minimum_k = 0;
        Stats minimum_stats{};
        for (std::int64_t trial = 0; trial < trials; ++trial) {
            const int k = 4 + static_cast<int>(
                generator() % static_cast<std::uint64_t>(max_k - 3)
            );
            Minus minus{
                1 + static_cast<int>(
                    generator() % static_cast<std::uint64_t>(k)
                ),
                1 + static_cast<int>(
                    generator() % static_cast<std::uint64_t>(k)
                )
            };
            std::sort(minus.begin(), minus.end());
            const bool is_shallow =
                std::min(cap(k, minus[0]), cap(k, minus[1])) <= 2;
            if (is_shallow != shallow_mode) {
                continue;
            }
            Plus plus{};
            for (int& label : plus) {
                label = 1 + static_cast<int>(
                    generator() % static_cast<std::uint64_t>(k)
                );
            }
            std::sort(plus.begin(), plus.end());
            if (!disjoint(minus, plus)) {
                continue;
            }
            const Stats s = statistics(k, minus, plus);
            if (s.max_cut == 0 || s.max_cut > 2) {
                continue;
            }
            ++tested;
            const int cutoff = shallow_mode ? 4 : 8;
            const std::int64_t direct =
                worst_partial_cut_bound(
                    k, minus, plus, cutoff, s.max_cut
                )
                + first_two_equal_plus_reservoir(k, minus, plus)
                + s.u1 + s.u2 + (minus[0] == minus[1] ? 1 : 0)
                - s.t;
            if (!minimum_set || direct < minimum) {
                minimum_set = true;
                minimum = direct;
                minimum_minus = minus;
                minimum_plus = plus;
                minimum_k = k;
                minimum_stats = s;
            }
            if (direct < 0) {
                break;
            }
        }
        std::cout << "random-direct mode="
                  << (shallow_mode ? "shallow" : "deep")
                  << " seed=" << seed
                  << " trials=" << trials
                  << " tested=" << tested;
        if (minimum_set) {
            std::cout << " min_direct=" << minimum << ' ';
            print_case(
                "minimum",
                minimum_k,
                minimum_minus,
                minimum_plus,
                minimum_stats,
                graph_stats(minimum_k, minimum_minus, minimum_plus)
            );
        } else {
            std::cout << " none\n";
        }
        return minimum_set && minimum < 0 ? EXIT_FAILURE : EXIT_SUCCESS;
    }

    if (argc == 10 && std::string(argv[1]) == "--case") {
        const int k = std::atoi(argv[2]);
        const Minus minus{std::atoi(argv[3]), std::atoi(argv[4])};
        const Plus plus{
            std::atoi(argv[5]),
            std::atoi(argv[6]),
            std::atoi(argv[7]),
            std::atoi(argv[8]),
            std::atoi(argv[9])
        };
        const Stats s = statistics(k, minus, plus);
        const std::int64_t peq_first_two =
            first_two_equal_plus_reservoir(k, minus, plus);
        std::cout << "case k=" << k
                  << " minus=[" << minus[0] << ',' << minus[1] << ']'
                  << " plus=[" << plus[0] << ',' << plus[1] << ','
                  << plus[2] << ',' << plus[3] << ',' << plus[4] << ']'
                  << " N=" << s.n
                  << " P=" << s.peq
                  << " Pminus=" << s.peq_minus
                  << " Pplus=" << s.peq_plus
                  << " Plow2=" << peq_first_two
                  << " U1=" << s.u1
                  << " U2=" << s.u2
                  << " T=" << s.t
                  << " signed=" << s.n + s.peq + s.u1 + s.u2 - s.t
                  << " d=" << s.max_cut << '\n';
        for (int orientation = 0; orientation < 2; ++orientation) {
            for (int first = 0; first < 5; ++first) {
                for (int second = first + 1; second < 5; ++second) {
                    const std::int64_t local2 = partial_cut_bound(
                        k, minus, plus, orientation, first, second, 2
                    );
                    const std::int64_t local4 = partial_cut_bound(
                        k, minus, plus, orientation, first, second, 4
                    );
                    const std::int64_t local6 = partial_cut_bound(
                        k, minus, plus, orientation, first, second, 6
                    );
                    const std::int64_t local8 = partial_cut_bound(
                        k, minus, plus, orientation, first, second, 8
                    );
                    if (local2 != 0 || local8 != 0) {
                        std::cout << "  orientation=" << orientation
                                  << " pair=(" << first << ',' << second
                                  << ") local2=" << local2
                                  << " local4=" << local4
                                  << " local6=" << local6
                                  << " local8=" << local8 << '\n';
                    }
                }
            }
        }
        return EXIT_SUCCESS;
    }

    if (argc == 4 && std::string(argv[1]) == "--interval-band3") {
        const int min_k = std::atoi(argv[2]);
        const int max_k = std::atoi(argv[3]);
        for (int k = min_k; k <= max_k; ++k) {
            std::array<std::pair<int, int>, 3> minima{{
                {0, 0},
                {1000000, 1000000},
                {1000000, 1000000}
            }};
            std::array<std::array<int, 4>, 3> witness{};
            for (int a = 0; a <= k; ++a)
            for (int b = a; b <= k; ++b) {
                const auto lhs = fusion_outputs(k, a, b);
                if (lhs.size() < 3) {
                    continue;
                }
                for (int c = 0; c <= k; ++c)
                for (int d = c; d <= k; ++d) {
                    const auto rhs = fusion_outputs(k, c, d);
                    if (rhs.size() < 3) {
                        continue;
                    }
                    const int overlap = intersection_size(lhs, rhs);
                    if (overlap < 1 || overlap > 2) {
                        continue;
                    }
                    const int b2 = interval_pair_profile(k, lhs, rhs, 2);
                    const int b4 = interval_pair_profile(k, lhs, rhs, 4);
                    if (std::pair{b2, b4}
                        < minima[static_cast<std::size_t>(overlap)]) {
                        minima[static_cast<std::size_t>(overlap)] = {b2, b4};
                        witness[static_cast<std::size_t>(overlap)] = {a, b, c, d};
                    }
                }
            }
            std::cout << "level " << k;
            for (int overlap = 1; overlap <= 2; ++overlap) {
                const auto [b2, b4] =
                    minima[static_cast<std::size_t>(overlap)];
                const auto labels = witness[static_cast<std::size_t>(overlap)];
                std::cout << " d=" << overlap << ":(" << b2 << ',' << b4
                          << ")@[" << labels[0] << ',' << labels[1] << "]|["
                          << labels[2] << ',' << labels[3] << ']';
            }
            std::cout << '\n';
        }
        return EXIT_SUCCESS;
    }

    if (argc == 4
        && (std::string(argv[1]) == "--no-j-onecolor-minima"
            || std::string(argv[1]) == "--no-j-onecolor-dense-endpoint")) {
        const bool dense_only =
            std::string(argv[1]) == "--no-j-onecolor-dense-endpoint";
        const int min_k = std::atoi(argv[2]);
        const int max_k = std::atoi(argv[3]);
        for (int k = min_k; k <= max_k; ++k) {
            bool set = false;
            Minus best_minus{};
            Plus best_plus{};
            Stats best_stats{};
            GraphStats best_graph{};
            std::int64_t best_margin = 0;
            std::int64_t best_peq_margin = 0;
            std::int64_t best_np = 0;
            bool np_set = false;
            Minus np_minus{};
            Plus np_plus{};
            Stats np_stats{};
            GraphStats np_graph{};
            bool peq_set = false;
            Minus peq_minus{};
            Plus peq_plus{};
            Stats peq_stats{};
            GraphStats peq_graph{};
            std::int64_t endpoint_cases = 0;
            std::int64_t equality_cases = 0;
            int greatest_depth = 0;
            Minus deepest_minus{};
            Plus deepest_plus{};
            Stats deepest_stats{};
            GraphStats deepest_graph{};
            std::map<std::int64_t, std::int64_t> small_margins;
            std::map<std::tuple<int, int, int, int>,
                     std::pair<std::int64_t, std::int64_t>> strata;
            std::map<std::tuple<int, int, int, int>,
                     std::pair<std::int64_t, std::int64_t>> cap_two_strata;
            for (int q = 1; q <= k; ++q) {
                for (int a = q; a <= k; ++a) {
                    const Minus minus{q, a};
                    if (std::min({q, a, k - q, k - a}) < 2) {
                        continue;
                    }
                    for (int p1 = 1; p1 <= k - 1; ++p1)
                    for (int p2 = p1; p2 <= k - 1; ++p2)
                    for (int p3 = p2; p3 <= k - 1; ++p3)
                    for (int p4 = p3; p4 <= k - 1; ++p4)
                    for (int p5 = p4; p5 <= k - 1; ++p5) {
                        const Plus plus{p1, p2, p3, p4, p5};
                        if (!disjoint(minus, plus)) {
                            continue;
                        }
                        const auto s = statistics(k, minus, plus);
                        if (s.max_cut == 0 || s.max_cut > 2) {
                            continue;
                        }
                        const auto g = graph_stats(k, minus, plus);
                        const bool one_color =
                            (g.edge_count[0] > 0) != (g.edge_count[1] > 0);
                        if (!one_color) {
                            continue;
                        }
                        const std::size_t color =
                            g.edge_count[0] > 0 ? 0U : 1U;
                        const std::size_t inactive = 1U - color;
                        const int c = g.edge_count[color];
                        const int h = g.weight_two[color];
                        if (dense_only
                            && !((h == 0 && c >= 8) || (h > 0 && c + h >= 12))) {
                            continue;
                        }
                        const int vertices = incident_vertex_mask(g.mask[color]);
                        bool low_incident = false;
                        int min_incident_cap = k + 1;
                        for (int i = 0; i < 5; ++i) {
                            if ((vertices & (1 << i)) != 0) {
                                min_incident_cap = std::min(
                                    min_incident_cap,
                                    cap(k, plus[static_cast<std::size_t>(i)])
                                );
                            }
                        }
                        low_incident = min_incident_cap <= 3;
                        if (dense_only) {
                            std::array<int, 5> low_positions{};
                            int low_count = 0;
                            for (int i = 0; i < 5; ++i) {
                                if (cap(k, plus[static_cast<std::size_t>(i)]) == 2) {
                                    low_positions[static_cast<std::size_t>(low_count)] =
                                        i;
                                    ++low_count;
                                }
                            }
                            if (low_count < 2) {
                                continue;
                            }
                            if (low_count == 2) {
                                const int maximal_mask = s.max_cut == 2
                                    ? g.weight_two_mask[color]
                                    : g.mask[color];
                                const int low_vertex_mask =
                                    (1 << low_positions[0])
                                    | (1 << low_positions[1]);
                                if ((incident_vertex_mask(maximal_mask)
                                     & low_vertex_mask) != 0) {
                                    continue;
                                }
                            }
                        } else if (
                            cap(k, minus[static_cast<std::size_t>(inactive)]) != 3
                            && !low_incident
                        ) {
                            continue;
                        }
                        ++endpoint_cases;
                        int case_depth = 0;
                        for (const int label : minus) {
                            case_depth = std::max(case_depth, cap(k, label));
                        }
                        for (const int label : plus) {
                            case_depth = std::max(case_depth, cap(k, label));
                        }
                        if (case_depth > greatest_depth) {
                            greatest_depth = case_depth;
                            deepest_minus = minus;
                            deepest_plus = plus;
                            deepest_stats = s;
                            deepest_graph = g;
                        }
                        const std::int64_t margin = s.n - s.t;
                        if (dense_only) {
                            int ones = 0;
                            int jones = 0;
                            for (const int label : plus) {
                                ones += label == 1 ? 1 : 0;
                                jones += label == k - 1 ? 1 : 0;
                            }
                            const auto key = std::make_tuple(ones, jones, c, h);
                            const std::int64_t raw_margin = s.n + s.peq - s.t;
                            auto [it, inserted] = cap_two_strata.emplace(
                                key,
                                std::pair<std::int64_t, std::int64_t>{
                                    0, raw_margin
                                }
                            );
                            ++it->second.first;
                            if (!inserted) {
                                it->second.second =
                                    std::min(it->second.second, raw_margin);
                            }
                        }
                        const std::int64_t peq_margin = s.peq - s.t;
                        if (!np_set || s.n + s.peq < best_np) {
                            np_set = true;
                            best_np = s.n + s.peq;
                            np_minus = minus;
                            np_plus = plus;
                            np_stats = s;
                            np_graph = g;
                        }
                        if (!peq_set || peq_margin < best_peq_margin) {
                            peq_set = true;
                            best_peq_margin = peq_margin;
                            peq_minus = minus;
                            peq_plus = plus;
                            peq_stats = s;
                            peq_graph = g;
                        }
                        const auto stratum_key = std::make_tuple(
                            c, h,
                            cap(k, minus[static_cast<std::size_t>(inactive)]) == 3,
                            std::min(min_incident_cap, 4)
                        );
                        auto [it, inserted] = strata.emplace(
                            stratum_key,
                            std::pair<std::int64_t, std::int64_t>{0, margin}
                        );
                        ++it->second.first;
                        if (!inserted) {
                            it->second.second = std::min(it->second.second, margin);
                        }
                        if (margin == 0) {
                            ++equality_cases;
                        }
                        if (s.t >= 8 && margin <= 10) {
                            ++small_margins[margin];
                            if (k == max_k) {
                                print_case("near-endpoint", k, minus, plus, s, g);
                            }
                        }
                        if (dense_only && k == max_k && s.n <= 19) {
                            print_case("small-N-endpoint", k, minus, plus, s, g);
                        }
                        if (!set || margin < best_margin) {
                            set = true;
                            best_margin = margin;
                            best_minus = minus;
                            best_plus = plus;
                            best_stats = s;
                            best_graph = g;
                        }
                    }
                }
            }
            std::cout << "level " << k << " endpoint_cases=" << endpoint_cases
                      << " equality_cases=" << equality_cases;
            std::cout << " dense_small_margins={";
            bool first = true;
            for (const auto& [margin, count] : small_margins) {
                std::cout << (first ? "" : ",") << margin << ':' << count;
                first = false;
            }
            std::cout << '}';
            if (set) {
                std::cout << " min_N_minus_T=" << best_margin << ' ';
                print_case("best", k, best_minus, best_plus,
                           best_stats, best_graph);
            } else {
                std::cout << " none\n";
            }
            if (k == max_k) {
                std::cout << "  greatest_cap=" << greatest_depth << ' ';
                if (greatest_depth > 0) {
                    print_case("deepest", k, deepest_minus, deepest_plus,
                               deepest_stats, deepest_graph);
                } else {
                    std::cout << "none\n";
                }
                std::cout << "  min_Peq_minus_T=" << best_peq_margin << ' ';
                if (peq_set) {
                    print_case("peq-worst", k, peq_minus, peq_plus,
                               peq_stats, peq_graph);
                } else {
                    std::cout << "none\n";
                }
                std::cout << "  min_N_plus_Peq=" << best_np << ' ';
                if (np_set) {
                    print_case("np-worst", k, np_minus, np_plus,
                               np_stats, np_graph);
                } else {
                    std::cout << "none\n";
                }
                for (const auto& [key, value] : strata) {
                    const auto [c, h, low_minus, min_plus_cap] = key;
                    std::cout << "  stratum c=" << c << " h=" << h
                              << " inactive_cap3=" << low_minus
                              << " min_incident_cap=" << min_plus_cap
                              << " count=" << value.first
                              << " min_N_minus_T=" << value.second << '\n';
                }
                for (const auto& [key, value] : cap_two_strata) {
                    const auto [ones, jones, c, h] = key;
                    std::cout << "  cap2-stratum n1=" << ones
                              << " nJ1=" << jones
                              << " c=" << c << " h=" << h
                              << " count=" << value.first
                              << " min_raw=" << value.second << '\n';
                }
            }
        }
        return EXIT_SUCCESS;
    }

    if (argc == 4 && std::string(argv[1]) == "--no-j-onecolor-summary") {
        const int min_k = std::atoi(argv[2]);
        const int max_k = std::atoi(argv[3]);
        for (int k = min_k; k <= max_k; ++k) {
            bool set = false;
            Minus best_minus{};
            Plus best_plus{};
            Stats best_stats{};
            GraphStats best_graph{};
            std::int64_t best_margin = 0;
            std::int64_t total = 0;
            for (int q = 1; q <= k; ++q) {
                for (int a = q; a <= k; ++a) {
                    const Minus minus{q, a};
                    if (std::min({q, a, k - q, k - a}) < 2) {
                        continue;
                    }
                    for (int p1 = 1; p1 <= k - 1; ++p1)
                    for (int p2 = p1; p2 <= k - 1; ++p2)
                    for (int p3 = p2; p3 <= k - 1; ++p3)
                    for (int p4 = p3; p4 <= k - 1; ++p4)
                    for (int p5 = p4; p5 <= k - 1; ++p5) {
                        const Plus plus{p1, p2, p3, p4, p5};
                        if (!disjoint(minus, plus)) {
                            continue;
                        }
                        const auto s = statistics(k, minus, plus);
                        if (s.max_cut == 0 || s.max_cut > 2) {
                            continue;
                        }
                        const auto g = graph_stats(k, minus, plus);
                        const bool one_color =
                            (g.edge_count[0] > 0) != (g.edge_count[1] > 0);
                        if (!one_color) {
                            continue;
                        }
                        const std::size_t color =
                            g.edge_count[0] > 0 ? 0U : 1U;
                        const int c = g.edge_count[color];
                        const int h = g.weight_two[color];
                        if (c + h <= s.n) {
                            continue;
                        }
                        ++total;
                        const std::int64_t margin = s.n + s.peq - (c + h);
                        if (!set || margin < best_margin) {
                            set = true;
                            best_margin = margin;
                            best_minus = minus;
                            best_plus = plus;
                            best_stats = s;
                            best_graph = g;
                        }
                    }
                }
            }
            std::cout << "level " << k << " candidates=" << total;
            if (set) {
                std::cout << " min_margin=" << best_margin << ' ';
                print_case("best", k, best_minus, best_plus,
                           best_stats, best_graph);
            } else {
                std::cout << " none\n";
            }
        }
        return EXIT_SUCCESS;
    }

    if (argc == 4
        && (std::string(argv[1]) == "--no-j-mixed-summary"
            || std::string(argv[1]) == "--no-j-mixed-minima"
            || std::string(argv[1]) == "--no-j-mixed-dense")) {
        const bool only_raw_deficits =
            std::string(argv[1]) == "--no-j-mixed-summary";
        const bool dense_only =
            std::string(argv[1]) == "--no-j-mixed-dense";
        const int min_k = std::atoi(argv[2]);
        const int max_k = std::atoi(argv[3]);
        for (int k = min_k; k <= max_k; ++k) {
            bool set = false;
            Minus best_minus{};
            Plus best_plus{};
            Stats best_stats{};
            GraphStats best_graph{};
            std::int64_t best_paid_margin = 0;
            std::int64_t best_raw_margin = 0;
            Minus raw_minus{};
            Plus raw_plus{};
            Stats raw_stats{};
            GraphStats raw_graph{};
            bool raw_set = false;
            std::int64_t best_double_margin = 0;
            bool double_set = false;
            Minus double_minus{};
            Plus double_plus{};
            Stats double_stats{};
            GraphStats double_graph{};
            std::int64_t best_n_margin = 0;
            bool n_set = false;
            Minus n_minus{};
            Plus n_plus{};
            Stats n_stats{};
            GraphStats n_graph{};
            constexpr std::array<int, 5> partial_outputs{4, 6, 8, 10, 12};
            std::array<std::int64_t, partial_outputs.size()> partial_violations{};
            std::array<std::int64_t, partial_outputs.size()> partial_min_margin{};
            std::array<bool, partial_outputs.size()> partial_set{};
            std::int64_t all_cut_violations = 0;
            std::int64_t all_cut_min_margin = 0;
            bool all_cut_set = false;
            std::int64_t maximal_cut_violations = 0;
            std::int64_t maximal_cut_min_margin = 0;
            bool maximal_cut_set = false;
            std::int64_t average_violations = 0;
            std::int64_t average_min_margin = 0;
            bool average_set = false;
            std::int64_t incidence_violations = 0;
            std::int64_t incidence_min_margin = 0;
            bool incidence_set = false;
            std::int64_t light_color_violations = 0;
            std::int64_t light_color_min_margin = 0;
            bool light_color_set = false;
            std::int64_t max_demand = 0;
            std::array<std::int64_t, 6> cap_two_max_demand{};
            std::array<std::int64_t, 6> cap_two_rank_one_max{};
            std::array<std::int64_t, 6> cap_two_min_margin{};
            std::array<bool, 6> cap_two_set{};
            std::array<std::array<std::int64_t, 3>, 2>
                cap_two_rank_two_class_max{};
            std::array<std::int64_t, 2> n2_t1_worst_violations{};
            std::array<std::int64_t, 2> n2_t1_best_violations{};
            std::array<std::int64_t, 2> n2_t1_worst_margin{};
            std::array<std::int64_t, 2> n2_t1_best_margin{};
            std::array<bool, 2> n2_t1_set{};
            std::array<std::int64_t, 2> n2_opposite_missing_endpoint{};
            std::set<std::tuple<int, int, int, int>> mixed_patterns;
            std::int64_t raw_deficits = 0;
            std::int64_t paid_deficits = 0;
            for (int q = 1; q <= k; ++q) {
                for (int a = q; a <= k; ++a) {
                    const Minus minus{q, a};
                    if (std::min({q, a, k - q, k - a}) < 2) {
                        continue;
                    }
                    for (int p1 = 1; p1 <= k - 1; ++p1)
                    for (int p2 = p1; p2 <= k - 1; ++p2)
                    for (int p3 = p2; p3 <= k - 1; ++p3)
                    for (int p4 = p3; p4 <= k - 1; ++p4)
                    for (int p5 = p4; p5 <= k - 1; ++p5) {
                        const Plus plus{p1, p2, p3, p4, p5};
                        if (!disjoint(minus, plus)) {
                            continue;
                        }
                        const auto s = statistics(k, minus, plus);
                        if (s.max_cut == 0 || s.max_cut > 2) {
                            continue;
                        }
                        const auto g = graph_stats(k, minus, plus);
                        if (minus[0] == minus[1]
                            || g.edge_count[0] == 0
                            || g.edge_count[1] == 0) {
                            continue;
                        }
                        const int c_total =
                            g.edge_count[0] + g.edge_count[1];
                        const int h_total =
                            g.weight_two[0] + g.weight_two[1];
                        if (dense_only
                            && !((h_total == 0 && c_total >= 8)
                                 || (h_total > 0
                                     && c_total + h_total >= 12))) {
                            continue;
                        }
                        if (only_raw_deficits && s.t <= s.n + s.peq) {
                            continue;
                        }
                        if (dense_only) {
                            max_demand = std::max(max_demand, s.t);
                            int cap_two_count_for_case = 0;
                            for (const int label : plus) {
                                cap_two_count_for_case +=
                                    cap(k, label) == 2 ? 1 : 0;
                            }
                            if (cap_two_count_for_case == 2
                                && s.max_cut == 2) {
                                int cap_mask = 0;
                                int low_count = 0;
                                int high_count = 0;
                                for (std::size_t vertex = 0;
                                     vertex < plus.size(); ++vertex) {
                                    if (plus[vertex] == 1) {
                                        ++low_count;
                                        cap_mask |= 1
                                            << static_cast<int>(vertex);
                                    } else if (plus[vertex] == k - 1) {
                                        ++high_count;
                                        cap_mask |= 1
                                            << static_cast<int>(vertex);
                                    }
                                }
                                std::array<bool, 3> has_incidence{};
                                bool has_low_incidence_one = false;
                                bool has_high_incidence_one = false;
                                int edge_bit = 0;
                                for (int i = 0; i < 5; ++i) {
                                    for (int j = i + 1; j < 5;
                                         ++j, ++edge_bit) {
                                        const int bit = 1 << edge_bit;
                                        if (((g.weight_two_mask[0]
                                              | g.weight_two_mask[1])
                                             & bit) == 0) {
                                            continue;
                                        }
                                        const int incidence =
                                            ((cap_mask >> i) & 1)
                                            + ((cap_mask >> j) & 1);
                                        has_incidence[
                                            static_cast<std::size_t>(
                                                incidence
                                            )
                                        ] = true;
                                        if (incidence == 1) {
                                            has_low_incidence_one =
                                                has_low_incidence_one
                                                || (((cap_mask >> i) & 1)
                                                    && plus[
                                                        static_cast<
                                                            std::size_t
                                                        >(i)
                                                    ] == 1)
                                                || (((cap_mask >> j) & 1)
                                                    && plus[
                                                        static_cast<
                                                            std::size_t
                                                        >(j)
                                                    ] == 1);
                                            has_high_incidence_one =
                                                has_high_incidence_one
                                                || (((cap_mask >> i) & 1)
                                                    && plus[
                                                        static_cast<
                                                            std::size_t
                                                        >(i)
                                                    ] == k - 1)
                                                || (((cap_mask >> j) & 1)
                                                    && plus[
                                                        static_cast<
                                                            std::size_t
                                                        >(j)
                                                    ] == k - 1);
                                        }
                                    }
                                }
                                const bool opposite =
                                    low_count == 1 && high_count == 1;
                                const std::size_t endpoint_type =
                                    opposite ? 1U : 0U;
                                if (opposite) {
                                    n2_opposite_missing_endpoint[0] +=
                                        has_low_incidence_one ? 0 : 1;
                                    n2_opposite_missing_endpoint[1] +=
                                        has_high_incidence_one ? 0 : 1;
                                }
                                const auto [worst_t1, best_t1] =
                                    rank_two_incidence_one_bounds(
                                        k, minus, plus, 4
                                    );
                                const std::int64_t worst_t1_margin =
                                    worst_t1 - s.t;
                                const std::int64_t best_t1_margin =
                                    best_t1 - s.t;
                                n2_t1_worst_violations[endpoint_type] +=
                                    worst_t1_margin < 0 ? 1 : 0;
                                n2_t1_best_violations[endpoint_type] +=
                                    best_t1_margin < 0 ? 1 : 0;
                                if (!n2_t1_set[endpoint_type]
                                    || worst_t1_margin
                                        < n2_t1_worst_margin[
                                            endpoint_type]) {
                                    n2_t1_worst_margin[endpoint_type] =
                                        worst_t1_margin;
                                }
                                if (!n2_t1_set[endpoint_type]
                                    || best_t1_margin
                                        < n2_t1_best_margin[
                                            endpoint_type]) {
                                    n2_t1_best_margin[endpoint_type] =
                                        best_t1_margin;
                                }
                                n2_t1_set[endpoint_type] = true;
                                const std::array<int, 3> priority =
                                    opposite
                                    ? std::array<int, 3>{1, 0, 2}
                                    : std::array<int, 3>{0, 1, 2};
                                for (const int incidence : priority) {
                                    if (!has_incidence[
                                            static_cast<std::size_t>(
                                                incidence
                                            )
                                        ]) {
                                        continue;
                                    }
                                    auto& maximum =
                                        cap_two_rank_two_class_max[
                                            opposite ? 1U : 0U
                                        ][static_cast<std::size_t>(
                                            incidence
                                        )];
                                    maximum = std::max(maximum, s.t);
                                    break;
                                }
                            }
                            mixed_patterns.emplace(
                                g.mask[0], g.mask[1],
                                g.weight_two_mask[0],
                                g.weight_two_mask[1]
                            );
                            const std::int64_t worst_bound =
                                worst_partial_cut_bound(k, minus, plus, 4, 0);
                            const std::int64_t worst_margin =
                                worst_bound - s.t;
                            all_cut_violations += worst_margin < 0 ? 1 : 0;
                            if (!all_cut_set
                                || worst_margin < all_cut_min_margin) {
                                all_cut_set = true;
                                all_cut_min_margin = worst_margin;
                            }
                            const std::int64_t worst_maximal_bound =
                                worst_partial_cut_bound(
                                    k, minus, plus, 4, s.max_cut
                                );
                            const std::int64_t maximal_margin =
                                worst_maximal_bound - s.t;
                            maximal_cut_violations +=
                                maximal_margin < 0 ? 1 : 0;
                            if (!maximal_cut_set
                                || maximal_margin
                                    < maximal_cut_min_margin) {
                                maximal_cut_set = true;
                                maximal_cut_min_margin = maximal_margin;
                            }
                            const std::int64_t bound_sum =
                                sum_partial_cut_bounds(k, minus, plus, 4);
                            const std::int64_t average_margin =
                                bound_sum - s.active * s.t;
                            average_violations +=
                                average_margin < 0 ? 1 : 0;
                            if (!average_set
                                || average_margin < average_min_margin) {
                                average_set = true;
                                average_min_margin = average_margin;
                            }
                            const std::int64_t incidence_bound =
                                best_max_weighted_incidence_bound(
                                    k, minus, plus, 4
                                );
                            const std::int64_t incidence_margin =
                                incidence_bound - s.t;
                            incidence_violations +=
                                incidence_margin < 0 ? 1 : 0;
                            if (!incidence_set
                                || incidence_margin
                                    < incidence_min_margin) {
                                incidence_set = true;
                                incidence_min_margin = incidence_margin;
                            }
                            const std::int64_t light_color_bound =
                                best_light_color_bound(
                                    k, minus, plus, 4
                                );
                            const std::int64_t light_color_margin =
                                light_color_bound - s.t;
                            light_color_violations +=
                                light_color_margin < 0 ? 1 : 0;
                            if (!light_color_set
                                || light_color_margin
                                    < light_color_min_margin) {
                                light_color_set = true;
                                light_color_min_margin = light_color_margin;
                            }
                            for (std::size_t b = 0;
                                 b < partial_outputs.size(); ++b) {
                                const std::int64_t bound =
                                    best_partial_cut_bound(
                                        k, minus, plus, partial_outputs[b]
                                    );
                                const std::int64_t margin = bound - s.t;
                                partial_violations[b] += margin < 0 ? 1 : 0;
                                if (!partial_set[b]
                                    || margin < partial_min_margin[b]) {
                                    partial_set[b] = true;
                                    partial_min_margin[b] = margin;
                                }
                                if (b == 0U) {
                                    int cap_two_count = 0;
                                    for (const int label : plus) {
                                        cap_two_count +=
                                            cap(k, label) == 2 ? 1 : 0;
                                    }
                                    const std::size_t index =
                                        static_cast<std::size_t>(
                                            cap_two_count
                                        );
                                    cap_two_max_demand[index] = std::max(
                                        cap_two_max_demand[index], s.t
                                    );
                                    if (h_total == 0) {
                                        cap_two_rank_one_max[index] =
                                            std::max(
                                                cap_two_rank_one_max[index],
                                                s.t
                                            );
                                    }
                                    if (!cap_two_set[index]
                                        || margin
                                            < cap_two_min_margin[index]) {
                                        cap_two_set[index] = true;
                                        cap_two_min_margin[index] = margin;
                                    }
                                }
                            }
                        }
                        ++raw_deficits;
                        const std::int64_t paid_margin =
                            s.n + s.peq + s.u1 + s.u2 - s.t;
                        const std::int64_t raw_margin = s.n + s.peq - s.t;
                        const std::int64_t double_margin =
                            s.n + s.peq - 2 * s.t;
                        const std::int64_t n_margin = s.n - s.t;
                        if (!n_set || n_margin < best_n_margin) {
                            n_set = true;
                            best_n_margin = n_margin;
                            n_minus = minus;
                            n_plus = plus;
                            n_stats = s;
                            n_graph = g;
                        }
                        if (!double_set || double_margin < best_double_margin) {
                            double_set = true;
                            best_double_margin = double_margin;
                            double_minus = minus;
                            double_plus = plus;
                            double_stats = s;
                            double_graph = g;
                        }
                        if (!raw_set || raw_margin < best_raw_margin) {
                            raw_set = true;
                            best_raw_margin = raw_margin;
                            raw_minus = minus;
                            raw_plus = plus;
                            raw_stats = s;
                            raw_graph = g;
                        }
                        paid_deficits += paid_margin < 0 ? 1 : 0;
                        if (!set || paid_margin < best_paid_margin) {
                            set = true;
                            best_paid_margin = paid_margin;
                            best_minus = minus;
                            best_plus = plus;
                            best_stats = s;
                            best_graph = g;
                        }
                        if (dense_only && k == max_k
                            && s.n + s.peq <= 40) {
                            print_case("small-NP-mixed", k, minus, plus, s, g);
                        }
                    }
                }
            }
            std::cout << "level " << k
                      << " raw_deficits=" << raw_deficits
                      << " paid_deficits=" << paid_deficits;
            if (set) {
                std::cout << " min_paid_margin=" << best_paid_margin << ' ';
                print_case("best-mixed", k, best_minus, best_plus,
                           best_stats, best_graph);
                std::cout << "  min_raw_margin=" << best_raw_margin << ' ';
                print_case("best-raw-mixed", k, raw_minus, raw_plus,
                           raw_stats, raw_graph);
                std::cout << "  min_double_margin=" << best_double_margin
                          << ' ';
                print_case("best-double-mixed", k, double_minus, double_plus,
                           double_stats, double_graph);
                std::cout << "  min_N_minus_T=" << best_n_margin << ' ';
                print_case("best-N-mixed", k, n_minus, n_plus,
                           n_stats, n_graph);
            } else {
                std::cout << " none\n";
            }
            if (dense_only) {
                std::cout << "  partial";
                for (std::size_t b = 0; b < partial_outputs.size(); ++b) {
                    std::cout << " j" << partial_outputs[b]
                              << "=(violations:" << partial_violations[b]
                              << ",min:" << partial_min_margin[b] << ')';
                }
                std::cout << " all_cuts_j4=(violations:"
                          << all_cut_violations
                          << ",min:" << all_cut_min_margin
                          << ") maximal_cuts_j4=(violations:"
                          << maximal_cut_violations
                          << ",min:" << maximal_cut_min_margin
                          << ") average_j4=(violations:"
                          << average_violations
                          << ",min:" << average_min_margin
                          << ") max_weighted_incidence_j4=(violations:"
                          << incidence_violations
                          << ",min:" << incidence_min_margin
                          << ") light_color_j4=(violations:"
                          << light_color_violations
                          << ",min:" << light_color_min_margin
                          << ") max_demand=" << max_demand
                          << " patterns=" << mixed_patterns.size() << '\n';
                std::cout << "  cap_two";
                for (std::size_t count = 0;
                     count < cap_two_set.size(); ++count) {
                    if (cap_two_set[count]) {
                        std::cout << " n" << count
                                  << "=(maxT:"
                                  << cap_two_max_demand[count]
                                  << ",rank1:"
                                  << cap_two_rank_one_max[count]
                                  << ",min:" << cap_two_min_margin[count]
                                  << ')';
                    }
                }
                std::cout << "\n  n2_rank2"
                          << " equal(t0,t1,t2)=("
                          << cap_two_rank_two_class_max[0][0] << ','
                          << cap_two_rank_two_class_max[0][1] << ','
                          << cap_two_rank_two_class_max[0][2] << ')'
                          << " opposite(t0,t1,t2)=("
                          << cap_two_rank_two_class_max[1][0] << ','
                          << cap_two_rank_two_class_max[1][1] << ','
                          << cap_two_rank_two_class_max[1][2] << ')'
                          << " t1_equal=(worst_viol:"
                          << n2_t1_worst_violations[0]
                          << ",worst_min:" << n2_t1_worst_margin[0]
                          << ",best_viol:" << n2_t1_best_violations[0]
                          << ",best_min:" << n2_t1_best_margin[0] << ')'
                          << " t1_opposite=(worst_viol:"
                          << n2_t1_worst_violations[1]
                          << ",worst_min:" << n2_t1_worst_margin[1]
                          << ",best_viol:" << n2_t1_best_violations[1]
                          << ",best_min:" << n2_t1_best_margin[1] << ')'
                          << " opposite_missing_t1=(low:"
                          << n2_opposite_missing_endpoint[0]
                          << ",high:" << n2_opposite_missing_endpoint[1]
                          << ")\n";
            }
        }
        return EXIT_SUCCESS;
    }

    if (argc == 4 && std::string(argv[1]) == "--lone-j-endpoint") {
        const int min_k = std::atoi(argv[2]);
        const int max_k = std::atoi(argv[3]);
        for (int k = min_k; k <= max_k; ++k) {
            const std::array<int, 4> minus_labels{3, 4, k - 4, k - 3};
            const std::array<int, 4> plus_x_labels{1, 2, k - 2, k - 1};
            std::map<std::tuple<int, int, int, int, int>, int> summary;
            std::set<Minus> minus_pairs;
            for (const int q : minus_labels) {
                if (q < 1 || q > k || std::min(q, k - q) < 2) {
                    continue;
                }
                for (const int a0 : {q, k - q}) {
                    const Minus minus{std::min(q, a0), std::max(q, a0)};
                    if (std::min({minus[0], minus[1], k - minus[0],
                                  k - minus[1]}) < 2) {
                        continue;
                    }
                    if (!minus_pairs.insert(minus).second) {
                        continue;
                    }
                    for (int i0 = 0; i0 < 4; ++i0)
                    for (int i1 = i0; i1 < 4; ++i1)
                    for (int i2 = i1; i2 < 4; ++i2)
                    for (int i3 = i2; i3 < 4; ++i3) {
                        Plus plus{
                            plus_x_labels[static_cast<std::size_t>(i0)],
                            plus_x_labels[static_cast<std::size_t>(i1)],
                            plus_x_labels[static_cast<std::size_t>(i2)],
                            plus_x_labels[static_cast<std::size_t>(i3)],
                            k
                        };
                        std::sort(plus.begin(), plus.end());
                        if (!disjoint(minus, plus)) {
                            continue;
                        }
                        const auto s = statistics(k, minus, plus);
                        if (s.max_cut == 0 || s.max_cut > 2) {
                            continue;
                        }
                        const auto g = graph_stats(k, minus, plus);
                        const int c = g.edge_count[0] + g.edge_count[1];
                        const int h = g.weight_two[0] + g.weight_two[1];
                        const auto key = std::make_tuple(
                            minus[0], minus[1], s.n + s.peq - s.t, c, h
                        );
                        ++summary[key];
                    }
                }
            }
            std::cout << "level " << k << '\n';
            for (const auto& [key, count] : summary) {
                const auto [q, a, raw, c, h] = key;
                std::cout << "  minus=[" << q << ',' << a << "] raw="
                          << raw << " c=" << c << " h=" << h
                          << " count=" << count << '\n';
            }
        }
        return EXIT_SUCCESS;
    }

    if (argc == 4
        && (std::string(argv[1]) == "--shallow-minus"
            || std::string(argv[1]) == "--deep-minus")) {
        const bool shallow_mode =
            std::string(argv[1]) == "--shallow-minus";
        const std::string mode_name =
            shallow_mode ? "shallow-minus" : "deep-minus";
        const int min_k = std::atoi(argv[2]);
        const int max_k = std::atoi(argv[3]);
        for (int k = min_k; k <= max_k; ++k) {
            std::int64_t cases = 0;
            std::int64_t active_cases = 0;
            std::int64_t failures = 0;
            bool signed_set = false;
            bool raw_set = false;
            bool frontier_set = false;
            bool local_set = false;
            bool equal_local_set = false;
            std::int64_t best_signed = 0;
            std::int64_t best_raw = 0;
            std::int64_t best_frontier = 0;
            std::int64_t best_local_four = 0;
            std::int64_t best_local_eight = 0;
            std::int64_t best_local_four_pair = 0;
            std::int64_t best_local_eight_pair = 0;
            std::int64_t best_equal_local = 0;
            std::int64_t best_equal_local_eight = 0;
            Minus signed_minus{};
            Minus raw_minus{};
            Minus frontier_minus{};
            Minus local_four_minus{};
            Minus local_eight_minus{};
            Plus signed_plus{};
            Plus raw_plus{};
            Plus frontier_plus{};
            Plus local_four_plus{};
            Plus local_eight_plus{};
            Minus local_four_pair_minus{};
            Minus local_eight_pair_minus{};
            Plus local_four_pair_plus{};
            Plus local_eight_pair_plus{};
            Minus equal_local_minus{};
            Plus equal_local_plus{};
            Minus equal_local_eight_minus{};
            Plus equal_local_eight_plus{};
            Stats signed_stats{};
            Stats raw_stats{};
            Stats frontier_stats{};
            Stats local_four_stats{};
            Stats local_eight_stats{};
            Stats local_four_pair_stats{};
            Stats local_eight_pair_stats{};
            Stats equal_local_stats{};
            Stats equal_local_eight_stats{};
            std::map<std::tuple<int, int, int>,
                     std::pair<std::int64_t, std::int64_t>> strata;
            for (int q = 1; q <= k; ++q) {
                for (int a = q; a <= k; ++a) {
                    const Minus minus{q, a};
                    const int q_cap = cap(k, q);
                    const int a_cap = cap(k, a);
                    if ((shallow_mode && std::min(q_cap, a_cap) > 2)
                        || (!shallow_mode
                            && std::min(q_cap, a_cap) <= 2)) {
                        continue;
                    }
                    for (int p1 = 1; p1 <= k; ++p1)
                    for (int p2 = p1; p2 <= k; ++p2)
                    for (int p3 = p2; p3 <= k; ++p3)
                    for (int p4 = p3; p4 <= k; ++p4)
                    for (int p5 = p4; p5 <= k; ++p5) {
                        const Plus plus{p1, p2, p3, p4, p5};
                        if (!disjoint(minus, plus)) {
                            continue;
                        }
                        ++cases;
                        const Stats s = statistics(k, minus, plus);
                        if (s.max_cut == 0 || s.max_cut > 2) {
                            continue;
                        }
                        ++active_cases;
                        const std::int64_t raw = s.n + s.peq - s.t;
                        const std::int64_t signed_margin =
                            raw + s.u1 + s.u2;
                        if (signed_margin < 0) {
                            ++failures;
                        }
                        if (!signed_set || signed_margin < best_signed) {
                            signed_set = true;
                            best_signed = signed_margin;
                            signed_minus = minus;
                            signed_plus = plus;
                            signed_stats = s;
                        }
                        if (!raw_set || raw < best_raw) {
                            raw_set = true;
                            best_raw = raw;
                            raw_minus = minus;
                            raw_plus = plus;
                            raw_stats = s;
                        }
                        if (!(q == 1 && a == 1)
                            && (!frontier_set
                                || signed_margin < best_frontier)) {
                            frontier_set = true;
                            best_frontier = signed_margin;
                            frontier_minus = minus;
                            frontier_plus = plus;
                            frontier_stats = s;
                        }
                        if (q != a) {
                            const std::int64_t positive = s.u1 + s.u2;
                            const std::int64_t peq_first_two =
                                first_two_equal_plus_reservoir(
                                    k, minus, plus
                                );
                            const std::int64_t local_four =
                                best_partial_cut_bound(k, minus, plus, 4)
                                + positive - s.t;
                            const std::int64_t local_eight =
                                best_partial_cut_bound(k, minus, plus, 8)
                                + positive - s.t;
                            const std::int64_t local_four_worst =
                                worst_partial_cut_bound(
                                    k, minus, plus, 4, s.max_cut
                                )
                                + positive - s.t;
                            const std::int64_t local_eight_worst =
                                worst_partial_cut_bound(
                                    k, minus, plus, 8, s.max_cut
                                )
                                + positive - s.t;
                            const std::int64_t local_four_pair =
                                local_four_worst + peq_first_two;
                            const std::int64_t local_eight_pair =
                                local_eight_worst + peq_first_two;
                            if (!local_set
                                || local_four < best_local_four) {
                                best_local_four = local_four;
                                local_four_minus = minus;
                                local_four_plus = plus;
                                local_four_stats = s;
                            }
                            if (!local_set
                                || local_eight < best_local_eight) {
                                best_local_eight = local_eight;
                                local_eight_minus = minus;
                                local_eight_plus = plus;
                                local_eight_stats = s;
                            }
                            if (!local_set
                                || local_four_pair
                                    < best_local_four_pair) {
                                best_local_four_pair = local_four_pair;
                                local_four_pair_minus = minus;
                                local_four_pair_plus = plus;
                                local_four_pair_stats = s;
                            }
                            if (!local_set
                                || local_eight_pair
                                    < best_local_eight_pair) {
                                best_local_eight_pair = local_eight_pair;
                                local_eight_pair_minus = minus;
                                local_eight_pair_plus = plus;
                                local_eight_pair_stats = s;
                            }
                            local_set = true;
                        }
                        if (q == a && q != 1) {
                            const std::int64_t peq_first_two =
                                first_two_equal_plus_reservoir(
                                    k, minus, plus
                                );
                            const std::int64_t equal_local =
                                worst_partial_cut_bound(
                                    k, minus, plus, 4, s.max_cut
                                )
                                + peq_first_two + s.u1 + s.u2
                                + 1 - s.t;
                            const std::int64_t equal_local_eight =
                                worst_partial_cut_bound(
                                    k, minus, plus, 8, s.max_cut
                                )
                                + peq_first_two + s.u1 + s.u2
                                + 1 - s.t;
                            if (!equal_local_set
                                || equal_local < best_equal_local) {
                                best_equal_local = equal_local;
                                equal_local_minus = minus;
                                equal_local_plus = plus;
                                equal_local_stats = s;
                            }
                            if (!equal_local_set
                                || equal_local_eight
                                    < best_equal_local_eight) {
                                best_equal_local_eight =
                                    equal_local_eight;
                                equal_local_eight_minus = minus;
                                equal_local_eight_plus = plus;
                                equal_local_eight_stats = s;
                            }
                            equal_local_set = true;
                        }
                        const auto key = std::make_tuple(
                            std::min(q_cap, 3),
                            std::min(a_cap, 3),
                            s.max_cut
                        );
                        auto [it, inserted] = strata.emplace(
                            key,
                            std::pair<std::int64_t, std::int64_t>{
                                0, signed_margin
                            }
                        );
                        ++it->second.first;
                        if (!inserted) {
                            it->second.second = std::min(
                                it->second.second, signed_margin
                            );
                        }
                    }
                }
            }
            std::cout << mode_name << " level=" << k
                      << " cases=" << cases
                      << " active_rank_le_2=" << active_cases
                      << " signed_failures=" << failures;
            if (signed_set) {
                std::cout << " min_signed=" << best_signed << ' ';
                print_case(
                    "signed-worst", k, signed_minus, signed_plus,
                    signed_stats, graph_stats(k, signed_minus, signed_plus)
                );
                std::cout << "  min_raw=" << best_raw << ' ';
                print_case(
                    "raw-worst", k, raw_minus, raw_plus,
                    raw_stats, graph_stats(k, raw_minus, raw_plus)
                );
                if (frontier_set) {
                    std::cout << "  min_non_q1_equal=" << best_frontier
                              << ' ';
                    print_case(
                        "frontier-worst", k, frontier_minus, frontier_plus,
                        frontier_stats,
                        graph_stats(k, frontier_minus, frontier_plus)
                    );
                }
                if (local_set) {
                    std::cout << "  distinct_min_local4_plus_U="
                              << best_local_four << ' ';
                    print_case(
                        "local4-worst", k, local_four_minus,
                        local_four_plus, local_four_stats,
                        graph_stats(k, local_four_minus, local_four_plus)
                    );
                    std::cout << "  distinct_min_local8_plus_U="
                              << best_local_eight << ' ';
                    print_case(
                        "local8-worst", k, local_eight_minus,
                        local_eight_plus, local_eight_stats,
                        graph_stats(k, local_eight_minus, local_eight_plus)
                    );
                    std::cout
                        << "  distinct_min_worstmax_local4_plus_Pplus_U="
                              << best_local_four_pair << ' ';
                    print_case(
                        "worstmax-local2-pair", k, local_four_pair_minus,
                        local_four_pair_plus, local_four_pair_stats,
                        graph_stats(
                            k, local_four_pair_minus, local_four_pair_plus
                        )
                    );
                    std::cout << "  distinct_min_local8_plus_Pplus_U="
                              << best_local_eight_pair << ' ';
                    print_case(
                        "local8-pair-worst", k, local_eight_pair_minus,
                        local_eight_pair_plus, local_eight_pair_stats,
                        graph_stats(
                            k, local_eight_pair_minus, local_eight_pair_plus
                        )
                    );
                }
                if (equal_local_set) {
                    std::cout
                        << "  equal_non_q1_min_worstmax_local4_plus_Plow_U="
                        << best_equal_local << ' ';
                    print_case(
                        "equal-local-worst", k, equal_local_minus,
                        equal_local_plus, equal_local_stats,
                        graph_stats(k, equal_local_minus, equal_local_plus)
                    );
                    std::cout
                        << "  equal_non_q1_min_worstmax_local8_plus_Plow_U="
                        << best_equal_local_eight << ' ';
                    print_case(
                        "equal-local8-worst", k,
                        equal_local_eight_minus,
                        equal_local_eight_plus,
                        equal_local_eight_stats,
                        graph_stats(
                            k,
                            equal_local_eight_minus,
                            equal_local_eight_plus
                        )
                    );
                }
            } else {
                std::cout << " none\n";
            }
            for (const auto& [key, value] : strata) {
                const auto [q_cap, a_cap, rank] = key;
                std::cout << "  shallow-stratum qcap=" << q_cap
                          << " acap=" << a_cap
                          << " d=" << rank
                          << " count=" << value.first
                          << " min_signed=" << value.second << '\n';
            }
        }
        return EXIT_SUCCESS;
    }

    if (argc != 3) {
        std::cerr << "usage: analyze_su2_d12_graph_residual MIN_K MAX_K\n"
                  << "   or: analyze_su2_d12_graph_residual"
                  << " --random-direct MODE SEED TRIALS MAX_K\n"
                  << "   or: analyze_su2_d12_graph_residual"
                  << " --case K Q A P0 P1 P2 P3 P4\n"
                  << "   or: analyze_su2_d12_graph_residual"
                  << " --no-j-onecolor-minima MIN_K MAX_K\n"
                  << "   or: analyze_su2_d12_graph_residual"
                  << " --no-j-onecolor-dense-endpoint MIN_K MAX_K\n"
                  << "   or: analyze_su2_d12_graph_residual"
                  << " --lone-j-endpoint MIN_K MAX_K\n"
                  << "   or: analyze_su2_d12_graph_residual"
                  << " --shallow-minus MIN_K MAX_K\n"
                  << "   or: analyze_su2_d12_graph_residual"
                  << " --deep-minus MIN_K MAX_K\n";
        return EXIT_FAILURE;
    }
    const int min_k = std::atoi(argv[1]);
    const int max_k = std::atoi(argv[2]);
    std::map<std::tuple<int, int, int, int, int>, std::array<int, 8>> counts;
    std::set<std::tuple<int, int, int, int, int>> printed;
    for (int k = min_k; k <= max_k; ++k) {
        for (int q = 1; q <= k; ++q) {
            for (int a = q; a <= k; ++a) {
                const Minus minus{q, a};
                if (std::min({q, a, k - q, k - a}) < 2) {
                    continue;
                }
                for (int p1 = 1; p1 <= k - 1; ++p1)
                for (int p2 = p1; p2 <= k - 1; ++p2)
                for (int p3 = p2; p3 <= k - 1; ++p3)
                for (int p4 = p3; p4 <= k - 1; ++p4)
                for (int p5 = p4; p5 <= k - 1; ++p5) {
                    const Plus plus{p1, p2, p3, p4, p5};
                    if (!disjoint(minus, plus)) {
                        continue;
                    }
                    const auto s = statistics(k, minus, plus);
                    if (s.max_cut == 0 || s.max_cut > 2 || s.t <= s.n) {
                        continue;
                    }
                    const auto g = graph_stats(k, minus, plus);
                    const bool one_color =
                        (g.edge_count[0] > 0) != (g.edge_count[1] > 0);
                    if (!one_color) {
                        continue;
                    }
                    const std::size_t color =
                        g.edge_count[0] > 0 ? 0U : 1U;
                    const int c = g.edge_count[color];
                    const int h = g.weight_two[color];
                    int low_cap = 0;
                    for (const int p : plus) {
                        if (cap(k, p) <= 3) {
                            ++low_cap;
                        }
                    }
                    const auto key = std::make_tuple(
                        c, h, s.n + s.peq - s.t < 0 ? 1 : 0,
                        s.u1 + s.u2 > 0 ? 1 : 0, low_cap
                    );
                    ++counts[key][0];
                    if (printed.insert(key).second) {
                        print_case("one-color", k, minus, plus, s, g);
                    }
                }
            }
        }
        std::cerr << "finished k=" << k << '\n';
    }
    for (const auto& [key, value] : counts) {
        const auto [c, h, raw_negative, has_u, low_cap] = key;
        std::cout << "summary c=" << c << " h=" << h
                  << " rawneg=" << raw_negative << " hasU=" << has_u
                  << " lowcap=" << low_cap << " count=" << value[0] << '\n';
    }
    return EXIT_SUCCESS;
}
