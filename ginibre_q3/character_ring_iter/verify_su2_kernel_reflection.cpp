#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

cpp_int binom(const int n, const int r) {
    if (r < 0 || r > n) {
        return 0;
    }
    const int rr = std::min(r, n - r);
    cpp_int value = 1;
    for (int i = 1; i <= rr; ++i) {
        value *= n - rr + i;
        value /= i;
    }
    return value;
}

cpp_int factorial(const int n) {
    cpp_int value = 1;
    for (int i = 2; i <= n; ++i) {
        value *= i;
    }
    return value;
}

cpp_int unrestricted(const int r, const int u) {
    if (std::abs(u) > r || ((r + u) & 1) != 0) {
        return 0;
    }
    return binom(r, (r + u) / 2);
}

cpp_int reflected_path(const int k, const int r, const int a) {
    const int h = k + 2;
    cpp_int value = 0;
    const int bound = r + 2;
    for (int j = -bound; j <= bound; ++j) {
        value += unrestricted(r, a + 2 * j * h);
        value -= unrestricted(r, a + 2 + 2 * j * h);
    }
    return value;
}

std::vector<std::vector<cpp_int>> adjacency_paths(
    const int k,
    const int max_r
) {
    std::vector<std::vector<cpp_int>> paths(
        static_cast<std::size_t>(max_r + 1),
        std::vector<cpp_int>(static_cast<std::size_t>(k + 1), 0)
    );
    paths[0][0] = 1;
    for (int r = 0; r < max_r; ++r) {
        for (int a = 0; a <= k; ++a) {
            if (a > 0) {
                paths[static_cast<std::size_t>(r + 1)]
                     [static_cast<std::size_t>(a - 1)]
                    += paths[static_cast<std::size_t>(r)]
                            [static_cast<std::size_t>(a)];
            }
            if (a < k) {
                paths[static_cast<std::size_t>(r + 1)]
                     [static_cast<std::size_t>(a + 1)]
                    += paths[static_cast<std::size_t>(r)]
                            [static_cast<std::size_t>(a)];
            }
        }
    }
    return paths;
}

cpp_int q_convolution(const int n, const int u, const int v) {
    cpp_int value = 0;
    for (int r = 0; r <= n; ++r) {
        value += binom(n, r) * unrestricted(r, u)
               * unrestricted(n - r, v);
    }
    return value;
}

cpp_int q_closed(const int n, const int u, const int v) {
    const int au = std::abs(u);
    const int av = std::abs(v);
    const int slack_numerator = n - au - av;
    if (slack_numerator < 0 || (slack_numerator & 1) != 0) {
        return 0;
    }
    const int t = slack_numerator / 2;
    const cpp_int nf = factorial(n);
    return nf * nf
         / (factorial(t) * factorial(au + t) * factorial(av + t)
            * factorial(au + av + t));
}

cpp_int l_direct(
    const int n,
    const int a,
    const int b,
    const std::vector<std::vector<cpp_int>>& paths
) {
    cpp_int value = 0;
    for (int r = 0; r <= n; ++r) {
        value += binom(n, r)
               * paths[static_cast<std::size_t>(r)]
                      [static_cast<std::size_t>(a)]
               * paths[static_cast<std::size_t>(n - r)]
                      [static_cast<std::size_t>(b)];
    }
    return value;
}

cpp_int l_reflected(const int k, const int n, const int a, const int b) {
    const int h = k + 2;
    const int bound = n + 2;
    cpp_int value = 0;
    for (int j = -bound; j <= bound; ++j) {
        for (int ell = -bound; ell <= bound; ++ell) {
            const int x0 = a + 2 * j * h;
            const int x1 = a + 2 + 2 * j * h;
            const int y0 = b + 2 * ell * h;
            const int y1 = b + 2 + 2 * ell * h;
            value += q_closed(n, x0, y0);
            value -= q_closed(n, x1, y0);
            value -= q_closed(n, x0, y1);
            value += q_closed(n, x1, y1);
        }
    }
    return value;
}

cpp_int l_stable(const int n, const int a, const int b) {
    const int slack_numerator = n - a - b;
    if (slack_numerator < 0 || (slack_numerator & 1) != 0) {
        return 0;
    }
    const int s = slack_numerator / 2;
    return factorial(n) * (a + 1) * (b + 1) * factorial(n + 2)
         / (factorial(s) * factorial(a + s + 1)
            * factorial(b + s + 1) * factorial(a + b + s + 2));
}

[[noreturn]] void mismatch(
    const std::string& kind,
    const int k,
    const int n,
    const int a,
    const int b,
    const cpp_int& lhs,
    const cpp_int& rhs
) {
    std::cerr << "FAIL kind=" << kind << " k=" << k << " n=" << n
              << " a=" << a << " b=" << b << " lhs=" << lhs
              << " rhs=" << rhs << '\n';
    throw std::runtime_error("exact verification mismatch");
}

}  // namespace

int main() {
    constexpr int max_k = 12;
    constexpr int max_n = 24;
    std::uint64_t path_checks = 0;
    std::uint64_t q_checks = 0;
    std::uint64_t reflection_checks = 0;
    std::uint64_t stable_checks = 0;

    for (int k = 0; k <= max_k; ++k) {
        const auto paths = adjacency_paths(k, max_n);
        for (int r = 0; r <= max_n; ++r) {
            for (int a = 0; a <= k; ++a) {
                const cpp_int reflected = reflected_path(k, r, a);
                const cpp_int direct
                    = paths[static_cast<std::size_t>(r)]
                           [static_cast<std::size_t>(a)];
                ++path_checks;
                if (reflected != direct) {
                    mismatch(
                        "one-dimensional-reflection",
                        k,
                        r,
                        a,
                        -1,
                        reflected,
                        direct
                    );
                }
            }
        }
    }

    for (int n = 0; n <= max_n; ++n) {
        for (int u = -max_n - 2; u <= max_n + 2; ++u) {
            for (int v = -max_n - 2; v <= max_n + 2; ++v) {
                const cpp_int convolution = q_convolution(n, u, v);
                const cpp_int closed = q_closed(n, u, v);
                ++q_checks;
                if (convolution != closed) {
                    mismatch(
                        "Q-closed-form",
                        -1,
                        n,
                        u,
                        v,
                        convolution,
                        closed
                    );
                }
            }
        }
    }

    for (int k = 0; k <= max_k; ++k) {
        const auto paths = adjacency_paths(k, max_n);
        for (int n = 0; n <= max_n; ++n) {
            for (int a = 0; a <= k; ++a) {
                for (int b = 0; b <= k; ++b) {
                    const cpp_int direct = l_direct(n, a, b, paths);
                    const cpp_int reflected = l_reflected(k, n, a, b);
                    ++reflection_checks;
                    if (direct != reflected) {
                        mismatch(
                            "two-dimensional-reflection",
                            k,
                            n,
                            a,
                            b,
                            direct,
                            reflected
                        );
                    }
                    if (k >= n) {
                        const cpp_int stable = l_stable(n, a, b);
                        ++stable_checks;
                        if (direct != stable) {
                            mismatch(
                                "stable-collapse",
                                k,
                                n,
                                a,
                                b,
                                direct,
                                stable
                            );
                        }
                    }
                }
            }
        }
    }

    std::cout << "PASS exact_cpp_int"
              << " max_k=" << max_k << " max_N=" << max_n
              << " path_checks=" << path_checks
              << " Q_checks=" << q_checks
              << " reflection_checks=" << reflection_checks
              << " stable_checks=" << stable_checks << '\n';
    return 0;
}
