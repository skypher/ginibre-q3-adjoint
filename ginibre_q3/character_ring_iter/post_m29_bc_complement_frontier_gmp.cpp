#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>

namespace {

constexpr double LOG10_2 = 0.301029995663981195213738894724493027;

struct Check {
    int q = 0;
    int window_bound = 0;
    int frontier = -1;
    int lhs_bits = 0;
    int rhs_bits = 0;
    int margin_bits = 0;
    double ok_ratio_log10 = -std::numeric_limits<double>::infinity();
    int next_k = -1;
    double next_ratio_log10 = std::numeric_limits<double>::infinity();
    bool ok = false;
};

double log10_abs_mpz(const mpz_class& value) {
    if (sgn(value) == 0) return -std::numeric_limits<double>::infinity();
    mpz_class abs_value = value >= 0 ? value : -value;
    long exponent = 0;
    const double mantissa = mpz_get_d_2exp(&exponent, abs_value.get_mpz_t());
    return std::log10(mantissa) + static_cast<double>(exponent) * LOG10_2;
}

std::string display(double value) {
    if (std::isinf(value)) return value > 0 ? "+inf" : "-inf";
    if (std::isnan(value)) return "nan";
    std::ostringstream out;
    out << std::setprecision(18) << value;
    return out.str();
}

int bit_length(const mpz_class& value) {
    if (sgn(value) <= 0) return 0;
    return static_cast<int>(mpz_sizeinbase(value.get_mpz_t(), 2));
}

int residual_window_bound(int q) {
    if (q == 15) return 912;
    if (16 <= q && q <= 17) return 401;
    if (q == 18) return 287;
    if (19 <= q && q <= 20) return 283;
    if (21 <= q && q <= 22) return 303;
    if (23 <= q && q <= 24) return 337;
    if (25 <= q && q <= 28) return 423;
    return 30899 - q;
}

std::vector<mpz_class> stable_moments(int max_index, bool progress) {
    std::vector<mpz_class> moments(static_cast<std::size_t>(max_index + 1));
    moments[0] = 1;
    if (max_index >= 1) moments[1] = 0;
    for (int index = 1; index < max_index; ++index) {
        mpz_class value = index * moments[static_cast<std::size_t>(index)]
            + index * moments[static_cast<std::size_t>(index - 1)];
        if (index >= 2) {
            value -= (mpz_class(index) * (index - 1) / 2)
                * moments[static_cast<std::size_t>(index - 2)];
        }
        if (sgn(value) <= 0) {
            std::cerr << "stable moment recurrence lost positivity at index "
                      << (index + 1) << "\n";
            std::cout << "__EXIT_STATUS=1\n";
            std::fflush(stdout);
            std::exit(1);
        }
        moments[static_cast<std::size_t>(index + 1)] = std::move(value);
        if (progress && ((index + 1) == max_index || (index + 1) % 10 == 0)) {
            std::cout << "stable_progress index=" << (index + 1)
                      << "/" << max_index << std::endl;
            std::fflush(stdout);
        }
    }
    return moments;
}

std::vector<mpz_class> cumulative_counts(const std::vector<mpz_class>& moments) {
    std::vector<mpz_class> cumulative(moments.size());
    mpz_class running = 0;
    for (std::size_t k = 0; k < moments.size(); ++k) {
        running += mpz_class(k + 1) * moments[k];
        cumulative[k] = running;
    }
    return cumulative;
}

}  // namespace

int main(int argc, char** argv) {
    int q_min = 15;
    int q_max = 218;
    int k_cap = 218;
    bool progress = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--q-min" && i + 1 < argc) {
            q_min = std::atoi(argv[++i]);
        } else if (arg == "--q-max" && i + 1 < argc) {
            q_max = std::atoi(argv[++i]);
        } else if (arg == "--k-cap" && i + 1 < argc) {
            k_cap = std::atoi(argv[++i]);
        } else if (arg == "--progress") {
            progress = true;
        } else {
            std::cerr
                << "usage: " << argv[0]
                << " [--q-min N] [--q-max N] [--k-cap N] [--progress]\n";
            std::cout << "__EXIT_STATUS=2\n";
            return 2;
        }
    }

    if (q_min < 15 || q_max < q_min || q_max > 218 || k_cap < 1) {
        std::cerr << "invalid complement-frontier verifier arguments\n";
        std::cout << "__EXIT_STATUS=2\n";
        return 2;
    }

    char hostname[256] = {};
    if (gethostname(hostname, sizeof(hostname) - 1) != 0) {
        hostname[0] = '\0';
    }

    std::cout << "B/C residual exact complement-frontier GMP verifier\n"
              << "VERIFIES exact frontier for window_bound(q) * M_K < 2^(2q)\n"
              << "where M_K=sum_{d=0}^K (d+1)*s_d and "
              << "s_{n+1}=n*s_n+n*s_{n-1}-binom(n,2)*s_{n-2}\n"
              << "window_bound(q) is the residual j-q bound from "
              << "cor:post29-bc-window-size-strata, with 30899-q for q>=29\n"
              << "host=" << (hostname[0] == '\0' ? "unknown" : hostname) << "\n"
              << "q_min=" << q_min
              << " q_max=" << q_max
              << " k_cap=" << k_cap
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;

    const std::vector<mpz_class> moments = stable_moments(k_cap, progress);
    const std::vector<mpz_class> cumulative = cumulative_counts(moments);

    std::vector<Check> checks(static_cast<std::size_t>(q_max - q_min + 1));
    int failures = 0;
    int completed = 0;

#pragma omp parallel for schedule(dynamic, 1) reduction(+:failures)
    for (int q = q_min; q <= q_max; ++q) {
        const int slot = q - q_min;
        const int window_bound = residual_window_bound(q);
        mpz_class rhs;
        mpz_ui_pow_ui(rhs.get_mpz_t(), 2UL, static_cast<unsigned long>(2 * q));

        int frontier = -1;
        mpz_class frontier_lhs = 0;
        for (int k = 0; k <= k_cap; ++k) {
            const mpz_class lhs = window_bound * cumulative[static_cast<std::size_t>(k)];
            if (lhs < rhs) {
                frontier = k;
                frontier_lhs = lhs;
            } else {
                break;
            }
        }

        Check check;
        check.q = q;
        check.window_bound = window_bound;
        check.frontier = frontier;
        check.ok = frontier >= 0 && frontier < k_cap && frontier <= q - 2;
        if (!check.ok) ++failures;

        if (frontier >= 0) {
            const mpz_class margin = rhs - frontier_lhs;
            check.lhs_bits = bit_length(frontier_lhs);
            check.rhs_bits = bit_length(rhs);
            check.margin_bits = bit_length(margin);
            check.ok_ratio_log10 = log10_abs_mpz(frontier_lhs) - log10_abs_mpz(rhs);
        }

        if (frontier + 1 <= k_cap) {
            check.next_k = frontier + 1;
            const mpz_class next_lhs =
                window_bound * cumulative[static_cast<std::size_t>(frontier + 1)];
            check.next_ratio_log10 = log10_abs_mpz(next_lhs) - log10_abs_mpz(rhs);
            if (next_lhs < rhs) {
                check.ok = false;
                ++failures;
            }
        } else {
            check.ok = false;
            ++failures;
        }

        checks[static_cast<std::size_t>(slot)] = check;

        if (progress) {
#pragma omp critical
            {
                ++completed;
                if (completed == static_cast<int>(checks.size()) || completed % 25 == 0) {
                    std::cout << "progress completed=" << completed
                              << "/" << checks.size()
                              << " last_q=" << q
                              << " failures=" << failures << std::endl;
                    std::fflush(stdout);
                }
            }
        }
    }

    Check worst_ok = checks.front();
    Check closest_next = checks.front();
    int min_frontier = std::numeric_limits<int>::max();
    int max_frontier = std::numeric_limits<int>::min();
    int max_frontier_minus_q = std::numeric_limits<int>::min();
    for (const Check& check : checks) {
        if (check.ok_ratio_log10 > worst_ok.ok_ratio_log10) worst_ok = check;
        if (check.next_ratio_log10 < closest_next.next_ratio_log10) closest_next = check;
        min_frontier = std::min(min_frontier, check.frontier);
        max_frontier = std::max(max_frontier, check.frontier);
        max_frontier_minus_q = std::max(max_frontier_minus_q, check.frontier - check.q);
    }

    for (const Check& check : checks) {
        std::cout << "frontier q=" << check.q
                  << " window_bound=" << check.window_bound
                  << " K=" << check.frontier
                  << " lhs_bits=" << check.lhs_bits
                  << " rhs_bits=" << check.rhs_bits
                  << " margin_bits=" << check.margin_bits
                  << " log10_lhs_over_rhs=" << display(check.ok_ratio_log10)
                  << " next_K=" << check.next_k
                  << " next_log10_lhs_over_rhs=" << display(check.next_ratio_log10)
                  << " status=" << (check.ok ? "OK" : "FAIL")
                  << std::endl;
    }

    const int anchors[] = {9, 10, 15, 39, 40, 83};
    for (int k : anchors) {
        if (k <= k_cap) {
            std::cout << "anchor_value M_" << k << "=" << cumulative[static_cast<std::size_t>(k)] << std::endl;
        }
    }

    std::cout << "SUMMARY failures=" << failures
              << " q_checked=" << checks.size()
              << " min_frontier=" << min_frontier
              << " max_frontier=" << max_frontier
              << " max_frontier_minus_q=" << max_frontier_minus_q
              << " worst_ok_q=" << worst_ok.q
              << " worst_ok_K=" << worst_ok.frontier
              << " worst_ok_log10_lhs_over_rhs=" << display(worst_ok.ok_ratio_log10)
              << " closest_next_q=" << closest_next.q
              << " closest_next_K=" << closest_next.next_k
              << " closest_next_log10_lhs_over_rhs=" << display(closest_next.next_ratio_log10)
              << std::endl;
    std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << std::endl;

    return failures == 0 ? 0 : 1;
}
