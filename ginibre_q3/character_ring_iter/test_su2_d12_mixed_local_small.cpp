#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

template <class F>
void outputs(const int k, const int a, const int b, F f) {
    const int upper = std::min(a + b, 2 * k - a - b);
    for (int c = std::abs(a - b); c <= upper; c += 2) {
        f(c);
    }
}

bool contains_output(const int k, const int a, const int b, const int target) {
    return std::abs(a - b) <= target
        && target <= std::min(a + b, 2 * k - a - b)
        && ((a + b + target) % 2 == 0);
}

std::vector<std::int64_t> decomposition(
    const int k,
    const std::vector<int>& labels
) {
    std::vector<std::int64_t> cur(static_cast<std::size_t>(k + 1));
    std::vector<std::int64_t> next(static_cast<std::size_t>(k + 1));
    cur[0] = 1;
    for (const int label : labels) {
        std::fill(next.begin(), next.end(), 0);
        for (int x = 0; x <= k; ++x) {
            const auto value = cur[static_cast<std::size_t>(x)];
            if (value == 0) {
                continue;
            }
            outputs(k, x, label, [&](const int y) {
                next[static_cast<std::size_t>(y)] += value;
            });
        }
        cur.swap(next);
    }
    return cur;
}

}  // namespace

int main() {
    std::int64_t checked = 0;
    for (int k = 6; k <= 8; ++k) {
        for (int q = 2; q <= k - 2; ++q) {
            for (int a = q + 1; a <= k - 2; ++a) {
                for (int p0 = 1; p0 <= k - 1; ++p0)
                for (int p1 = p0; p1 <= k - 1; ++p1)
                for (int p2 = p1; p2 <= k - 1; ++p2)
                for (int p3 = p2; p3 <= k - 1; ++p3)
                for (int p4 = p3; p4 <= k - 1; ++p4) {
                    const std::array<int, 5> plus{p0, p1, p2, p3, p4};
                    bool disjoint = true;
                    for (const int p : plus) {
                        disjoint = disjoint && p != q && p != a;
                    }
                    if (!disjoint) {
                        continue;
                    }

                    int c = 0;
                    int h = 0;
                    std::array<int, 2> colors{};
                    std::int64_t best_local = 0;
                    bool valid = true;
                    for (std::size_t i = 0; i < plus.size(); ++i) {
                        for (std::size_t j = i + 1; j < plus.size(); ++j) {
                            std::vector<int> rest;
                            for (std::size_t u = 0; u < plus.size(); ++u) {
                                if (u != i && u != j) {
                                    rest.push_back(plus[u]);
                                }
                            }
                            for (std::size_t orientation = 0;
                                 orientation < 2U; ++orientation) {
                                const int r = orientation == 0U ? q : a;
                                const int s = orientation == 0U ? a : q;
                                if (!contains_output(
                                        k, plus[i], plus[j], r
                                    )) {
                                    continue;
                                }
                                auto block = rest;
                                block.push_back(s);
                                const auto complement =
                                    decomposition(k, block);
                                const auto rank = complement[0];
                                if (rank == 0) {
                                    continue;
                                }
                                if (rank > 2) {
                                    valid = false;
                                    continue;
                                }

                                ++c;
                                ++colors[orientation];
                                h += rank == 2 ? 1 : 0;
                                const auto active_profile = decomposition(
                                    k, {r, plus[i], plus[j]}
                                );
                                std::int64_t local = 0;
                                for (int output = 0;
                                     output <= std::min(k, 4); ++output) {
                                    local += active_profile[
                                        static_cast<std::size_t>(output)
                                    ] * complement[
                                        static_cast<std::size_t>(output)
                                    ];
                                }
                                best_local = std::max(best_local, local);
                            }
                        }
                    }

                    if (!valid || colors[0] == 0 || colors[1] == 0) {
                        continue;
                    }
                    const int demand = c + h;
                    const bool dense =
                        (h == 0 && c >= 8) || (h > 0 && demand >= 12);
                    if (!dense) {
                        continue;
                    }
                    ++checked;
                    if (best_local < demand) {
                        std::cerr
                            << "SU2_D12_MIXED_LOCAL_SMALL FAIL"
                            << " k=" << k
                            << " minus=[" << q << ',' << a << ']'
                            << " plus=[" << p0 << ',' << p1 << ',' << p2
                            << ',' << p3 << ',' << p4 << ']'
                            << " local=" << best_local
                            << " demand=" << demand << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
        }
    }

    std::cout
        << "SU2_D12_MIXED_LOCAL_SMALL levels=6..8"
        << " cases=" << checked
        << " channels=0,2,4 result=PASS\n";
    return EXIT_SUCCESS;
}
