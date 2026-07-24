#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
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
        s.peq += mult(k, {plus[0], plus[1], plus[2], plus[3], plus[4]});
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
                s.peq += mult(k, positive_rest);
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

struct GraphStats {
    std::array<int, 2> edge_count{};
    std::array<int, 2> weight_two{};
    std::array<int, 2> mask{};
};

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
                }
                g.mask[o] |= 1 << bit;
            }
        }
    }
    return g;
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
              << ") mask=(" << g.mask[0] << ',' << g.mask[1] << ")"
              << " d=" << s.max_cut << " N=" << s.n << " P=" << s.peq
              << " U1=" << s.u1 << " U2=" << s.u2 << " T=" << s.t
              << " raw=" << s.n + s.peq - s.t
              << " paid=" << s.n + s.peq + s.u1 + s.u2 - s.t << '\n';
}

}  // namespace

int main(int argc, char** argv) {
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
                        const int color = g.edge_count[0] > 0 ? 0 : 1;
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

    if (argc != 3) {
        std::cerr << "usage: analyze_su2_d12_graph_residual MIN_K MAX_K\n"
                  << "   or: analyze_su2_d12_graph_residual"
                  << " --lone-j-endpoint MIN_K MAX_K\n";
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
                    const int color = g.edge_count[0] > 0 ? 0 : 1;
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
