#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Minus = std::array<int, 2>;
using Plus = std::array<int, 5>;

template <class F>
void outputs(int k, int a, int b, F f) {
    const int hi = std::min(a + b, 2 * k - a - b);
    for (int c = std::abs(a - b); c <= hi; c += 2) f(c);
}

std::int64_t mult(int k, std::vector<int> labels) {
    std::array<std::int64_t, 128> a{}, b{};
    a[0] = 1;
    for (int p : labels) {
        b.fill(0);
        for (int x = 0; x <= k; ++x) {
            const auto source = static_cast<std::size_t>(x);
            if (!a[source]) continue;
            outputs(k, x, p, [&](int y) {
                b[static_cast<std::size_t>(y)] += a[source];
            });
        }
        a = b;
    }
    return a[0];
}

bool disjoint(const Minus& m, const Plus& p) {
    for (int x : m) for (int y : p) if (x == y) return false;
    return true;
}

struct Stats {
    std::int64_t n{}, peq{}, u1{}, u2{}, t{};
    int c{};
    std::int64_t d{};
};

Stats stats(int k, const Minus& minus, const Plus& plus) {
    Stats s;
    s.n = mult(k, {
        minus[0], minus[1], plus[0], plus[1], plus[2], plus[3], plus[4]
    });
    if (minus[0] == minus[1]) {
        s.peq += mult(k, {plus[0], plus[1], plus[2], plus[3], plus[4]});
    }
    for (std::size_t i = 0; i < plus.size(); ++i) {
        std::vector<int> r;
        for (std::size_t j = 0; j < plus.size(); ++j) {
            if (j != i) r.push_back(plus[j]);
        }
        s.u1 += mult(k, {minus[0], minus[1], plus[i]}) * mult(k, r);
    }
    for (std::size_t i = 0; i < plus.size(); ++i) {
      for (std::size_t j = i + 1U; j < plus.size(); ++j) {
        std::vector<int> r, b{minus[0], minus[1]};
        for (std::size_t h = 0; h < plus.size(); ++h) if (h != i && h != j) {
            r.push_back(plus[h]);
            b.push_back(plus[h]);
        }
        if (plus[i] == plus[j]) s.peq += mult(k, b);
        s.u2 += mult(k, {minus[0], minus[1], plus[i], plus[j]})
              * mult(k, r);
        for (std::size_t o = 0; o < minus.size(); ++o) {
            const auto x = mult(k, {minus[o], plus[i], plus[j]});
            if (!x) continue;
            std::vector<int> z{minus[1U-o]};
            z.insert(z.end(), r.begin(), r.end());
            const auto v = x * mult(k, z);
            if (v) {
                ++s.c;
                s.t += v;
                s.d = std::max(s.d, v);
            }
        }
      }
    }
    return s;
}

std::string key(const Minus& m, const Plus& p) {
    std::ostringstream out;
    out << m[0] << ',' << m[1] << ':';
    for (int x : p) out << x << ',';
    return out.str();
}

struct Case {
    int k{}, placement{}, top_t{}, high{};
    Minus minus{};
    Plus plus{};
    Stats s{};
};

int high_cut_count(int k, int placement, const Minus& minus, const Plus& plus) {
    int total = 0;
    for (std::size_t i = 0; i < plus.size(); ++i) {
      for (std::size_t j = i + 1U; j < plus.size(); ++j) {
        for (std::size_t o = 0; o < minus.size(); ++o) {
            const bool distinguished =
                placement == 0
                    ? minus[o] == k
                    : plus[i] == k || plus[j] == k;
            if (!distinguished) continue;
            if (!mult(k, {minus[o], plus[i], plus[j]})) continue;
            std::vector<int> complement{minus[1U-o]};
            for (std::size_t h = 0; h < plus.size(); ++h) if (h != i && h != j) {
                complement.push_back(plus[h]);
            }
            if (mult(k, complement)) ++total;
        }
      }
    }
    return total;
}

int coefficient(int k, int a, int b, int c) {
    int value = 0;
    outputs(k, a, b, [&](int x) { value += x == c; });
    return value;
}

std::vector<int> erase_once(std::vector<int> values, int x) {
    const auto it = std::find(values.begin(), values.end(), x);
    if (it == values.end()) return {};
    values.erase(it);
    return values;
}

std::array<int, 4> top_edge_formula(const Case& z) {
    const int k = z.k;
    if (z.top_t != 1) return {-1, -1, -1, -1};
    if (z.placement == 0) {
        const int a = z.minus[0] == k ? z.minus[1] : z.minus[0];
        auto r = erase_once(
            erase_once(std::vector<int>(z.plus.begin(), z.plus.end()), 1),
            k-1
        );
        if (r.size() != 3) return {-2, -2, -2, -2};
        int w = 0, e1 = 0, eh = 0, e0 = 0, e2 = 0, rho = 0;
        for (std::size_t i = 0; i < r.size(); ++i) {
            const std::size_t j = (i + 1U) % 3U;
            const std::size_t l = (i + 2U) % 3U;
            e1 += std::abs(a-r[i]) == 1 && std::abs(r[j]-r[l]) == 1;
            eh += (a+r[i] == k-1 || a+r[i] == k+1)
               && (r[j]+r[l] == k-1 || r[j]+r[l] == k+1);
            if (r[i] == 2) e2 += coefficient(k, r[j], r[l], a);
        }
        auto kr = r;
        kr.push_back(k);
        e0 = (a == k-2) * static_cast<int>(mult(k, kr));
        for (std::size_t i = 0; i < r.size(); ++i) {
          for (std::size_t j = i+1U; j < r.size(); ++j) {
            if (r[i]+r[j] != k) continue;
            ++rho;
            const std::size_t l = 3U-i-j;
            w += a+r[l] == k-2 || a+r[l] == k+2;
            w += 2 * (a+r[l] == k);
          }
        }
        return {3+w+e1+eh+e0+e2, rho, w, e1+eh+e0+e2};
    }
    if (std::find(z.minus.begin(), z.minus.end(), 1) == z.minus.end()) {
        return {-1, -1, -1, -1};
    }
    const int a = z.minus[0] == 1 ? z.minus[1] : z.minus[0];
    const int y = k-a;
    auto r = erase_once(
        erase_once(std::vector<int>(z.plus.begin(), z.plus.end()), k),
        k-1
    );
    if (r.size() != 3) return {-2, -2, -2, -2};
    int f = 0, n = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
        const std::size_t j = (i + 1U) % 3U;
        const std::size_t l = (i + 2U) % 3U;
        n += r[i] == y;
        if (r[i] == k-2) f += coefficient(k, r[j], r[l], y);
        f += (r[i] == y) * (
            coefficient(k, r[j], r[l], k-2)
            + coefficient(k, r[j], r[l], k)
        );
        f += (a+r[i] == k-1 || a+r[i] == k+1)
           * coefficient(k, r[j], r[l], k-1);
        if (r[i] == 2) f += coefficient(k, r[j], r[l], a);
    }
    for (std::size_t i = 0; i < r.size(); ++i) {
      for (std::size_t j = i+1U; j < r.size(); ++j) {
        const std::size_t l = 3U-i-j;
        f += std::abs(r[i]-r[j]) == 1 && std::abs(r[l]-a) == 1;
      }
    }
    return {3+f, n, 0, f};
}

}  // namespace

int main(int argc, char** argv) {
    const int max_k = argc > 1 ? std::stoi(argv[1]) : 18;
    std::vector<Case> cases;

    for (int k = 4; k <= max_k; ++k) {
        std::map<std::string, Case> level;

        // Placement 0: k is a minus label and {u,k-u} is the selected pair.
        for (int u = 1; u <= k / 2; ++u) {
            for (int a = 1; a < k; ++a) {
                for (int r1 = 1; r1 <= k; ++r1)
                for (int r2 = r1; r2 <= k; ++r2)
                for (int r3 = r2; r3 <= k; ++r3) {
                    Minus m{a, k};
                    Plus p{u, k-u, r1, r2, r3};
                    std::sort(m.begin(), m.end());
                    std::sort(p.begin(), p.end());
                    if (!disjoint(m, p)) continue;
                    bool all_even = true;
                    for (int z : m) all_even = all_even && z % 2 == 0;
                    for (int z : p) all_even = all_even && z % 2 == 0;
                    if (all_even) continue;
                    if (mult(k, {a, r1, r2, r3}) != 3) continue;
                    const auto s = stats(k, m, p);
                    if (s.d != 3) continue;
                    auto z = Case{
                        k, 0, u, high_cut_count(k, 0, m, p), m, p, s
                    };
                    level.emplace(key(m, p), z);
                }
            }
        }

        // Placement 1: k is a plus label, q and x=k-q complete the top triple.
        for (int q = 1; q <= k-1; ++q) {
            const int x = k-q;
            for (int a = 1; a <= k; ++a) {
                Minus m{q, a};
                for (int r1 = 1; r1 <= k; ++r1)
                for (int r2 = r1; r2 <= k; ++r2)
                for (int r3 = r2; r3 <= k; ++r3) {
                    Plus p{k, x, r1, r2, r3};
                    std::sort(m.begin(), m.end());
                    std::sort(p.begin(), p.end());
                    if (!disjoint(m, p)) continue;
                    bool all_even = true;
                    for (int z : m) all_even = all_even && z % 2 == 0;
                    for (int z : p) all_even = all_even && z % 2 == 0;
                    if (all_even) continue;
                    if (mult(k, {a, r1, r2, r3}) != 3) continue;
                    const auto s = stats(k, m, p);
                    if (s.d != 3) continue;
                    auto z = Case{
                        k, 1, std::min(q, x), high_cut_count(k, 1, m, p),
                        m, p, s
                    };
                    auto [it, inserted] = level.emplace(key(m, p), z);
                    if (!inserted) it->second.top_t = std::max(it->second.top_t, z.top_t);
                }
            }
        }

        for (const auto& [_, z] : level) cases.push_back(z);
        std::cerr << "finished k=" << k << " cases=" << level.size() << '\n';
    }

    std::sort(cases.begin(), cases.end(), [](const Case& x, const Case& y) {
        const auto xm = x.s.n + x.s.peq - x.s.t;
        const auto ym = y.s.n + y.s.peq - y.s.t;
        return std::tie(xm, x.k, x.placement, x.top_t, x.minus, x.plus)
             < std::tie(ym, y.k, y.placement, y.top_t, y.minus, y.plus);
    });

    std::map<std::tuple<int,int,int>, std::array<std::int64_t, 6>> summary;
    std::map<std::pair<int,int>, int> top_edge_maximum;
    std::map<std::pair<int,int>, std::array<int, 2>> top_edge_parts;
    for (const auto& z : cases) {
        const auto [formula, parameter, first_part, second_part] =
            top_edge_formula(z);
        if (formula >= 0) {
            if (formula != z.s.t) {
                std::cerr << "top-edge formula mismatch k=" << z.k
                          << " formula=" << formula << " T=" << z.s.t << '\n';
                return 2;
            }
            top_edge_maximum[{z.placement, parameter}] = std::max(
                top_edge_maximum[{z.placement, parameter}], formula
            );
            auto& parts = top_edge_parts[{z.placement, parameter}];
            parts[0] = std::max(parts[0], first_part);
            parts[1] = std::max(parts[1], second_part);
            if (z.placement == 0 && parameter == 1 && second_part > 4) {
                std::cout << "top-edge-exception k=" << z.k
                          << " minus=[" << z.minus[0] << ',' << z.minus[1]
                          << "] plus=[";
                for (std::size_t i = 0; i < z.plus.size(); ++i) {
                    std::cout << (i ? "," : "") << z.plus[i];
                }
                std::cout << "] W=" << first_part << " E=" << second_part
                          << " T=" << formula << '\n';
            }
        }
        const auto raw = z.s.n + z.s.peq - z.s.t;
        const auto paid = raw + z.s.u1 + z.s.u2;
        auto& a = summary[{z.placement, std::min(z.top_t, 3), z.high}];
        if (a[0] == 0) {
            a = {
                1, raw, paid, z.s.n, z.s.t,
                z.s.peq + z.s.u1 + z.s.u2
            };
        } else {
            ++a[0];
            a[1] = std::min(a[1], raw);
            a[2] = std::min(a[2], paid);
            a[3] = std::min(a[3], z.s.n);
            a[4] = std::max(a[4], z.s.t);
            a[5] = std::min(a[5], z.s.peq + z.s.u1 + z.s.u2);
        }
    }
    for (const auto& [ix, a] : summary) {
        std::cout << "summary placement=" << std::get<0>(ix)
                  << " t=" << std::get<1>(ix)
                  << " high=" << std::get<2>(ix)
                  << " count=" << a[0] << " rawmin=" << a[1]
                  << " paidmin=" << a[2] << " Nmin=" << a[3]
                  << " Tmax=" << a[4] << " supplymin=" << a[5] << '\n';
    }
    for (const auto& [ix, value] : top_edge_maximum) {
        std::cout << "top-edge placement=" << ix.first
                  << " parameter=" << ix.second << " Tmax=" << value
                  << " firstmax=" << top_edge_parts[ix][0]
                  << " secondmax=" << top_edge_parts[ix][1] << '\n';
    }

    auto crude_cases = cases;
    std::sort(crude_cases.begin(), crude_cases.end(), [](const Case& x, const Case& y) {
        const auto xd = x.s.c + 2 * x.high - x.s.n - x.s.peq;
        const auto yd = y.s.c + 2 * y.high - y.s.n - y.s.peq;
        return std::tie(yd, x.k, x.placement, x.top_t, x.minus, x.plus)
             < std::tie(xd, y.k, y.placement, y.top_t, y.minus, y.plus);
    });
    int crude_printed = 0;
    for (const auto& z : crude_cases) {
        const auto deficit = z.s.c + 2 * z.high - z.s.n - z.s.peq;
        if (deficit <= 0 || crude_printed >= 400) break;
        std::cout << "crude k=" << z.k << " placement=" << z.placement
                  << " top_t=" << z.top_t << " high=" << z.high
                  << " minus=[" << z.minus[0] << ',' << z.minus[1]
                  << "] plus=[";
        for (std::size_t i = 0; i < z.plus.size(); ++i) {
            std::cout << (i ? "," : "") << z.plus[i];
        }
        std::cout << "] c=" << z.s.c << " d=" << z.s.d
                  << " N=" << z.s.n << " P=" << z.s.peq
                  << " U1=" << z.s.u1 << " U2=" << z.s.u2
                  << " T=" << z.s.t << " deficit=" << deficit << '\n';
        ++crude_printed;
    }

    auto repetition_cases = cases;
    std::sort(
        repetition_cases.begin(), repetition_cases.end(),
        [](const Case& x, const Case& y) {
            const auto xs = x.s.peq + x.s.u1 + x.s.u2;
            const auto ys = y.s.peq + y.s.u1 + y.s.u2;
            return std::tie(xs, x.s.n, x.s.t, x.k, x.placement, x.minus, x.plus)
                 < std::tie(ys, y.s.n, y.s.t, y.k, y.placement, y.minus, y.plus);
        }
    );
    int repetition_printed = 0;
    for (const auto& z : repetition_cases) {
        if ((z.high < 6 && !(z.top_t == 2 && z.high >= 3))
            || repetition_printed >= 240) continue;
        std::cout << "repeat k=" << z.k << " placement=" << z.placement
                  << " top_t=" << z.top_t << " high=" << z.high
                  << " minus=[" << z.minus[0] << ',' << z.minus[1]
                  << "] plus=[";
        for (std::size_t i = 0; i < z.plus.size(); ++i) {
            std::cout << (i ? "," : "") << z.plus[i];
        }
        std::cout << "] c=" << z.s.c << " N=" << z.s.n
                  << " P=" << z.s.peq << " U1=" << z.s.u1
                  << " U2=" << z.s.u2 << " T=" << z.s.t << '\n';
        ++repetition_printed;
    }

    int printed = 0;
    for (const auto& z : cases) {
        const auto raw = z.s.n + z.s.peq - z.s.t;
        if (raw >= 0 && printed >= 240) continue;
        std::cout << "case k=" << z.k << " placement=" << z.placement
                  << " top_t=" << z.top_t
                  << " high=" << z.high
                  << " minus=[" << z.minus[0] << ',' << z.minus[1]
                  << "] plus=[";
        for (std::size_t i = 0; i < z.plus.size(); ++i) {
            std::cout << (i ? "," : "") << z.plus[i];
        }
        std::cout << "] c=" << z.s.c << " d=" << z.s.d
                  << " N=" << z.s.n << " P=" << z.s.peq
                  << " U1=" << z.s.u1 << " U2=" << z.s.u2
                  << " T=" << z.s.t << " raw=" << raw
                  << " paid=" << raw + z.s.u1 + z.s.u2 << '\n';
        ++printed;
    }
}
