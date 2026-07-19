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
    int k = 0;
    int lhs_bits = 0;
    int rhs_bits = 0;
    int margin_bits = 0;
    double ratio_log10 = -std::numeric_limits<double>::infinity();
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

int ceil_div_3(int q) {
    return (q + 2) / 3;
}

int bit_length(const mpz_class& value) {
    if (sgn(value) <= 0) return 0;
    return static_cast<int>(mpz_sizeinbase(value.get_mpz_t(), 2));
}

}  // namespace

int main(int argc, char** argv) {
    int q_min = 86;
    int q_max = 218;
    int offset = 10;
    int row_bound = 30899;
    bool progress = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--q-min" && i + 1 < argc) {
            q_min = std::atoi(argv[++i]);
        } else if (arg == "--q-max" && i + 1 < argc) {
            q_max = std::atoi(argv[++i]);
        } else if (arg == "--offset" && i + 1 < argc) {
            offset = std::atoi(argv[++i]);
        } else if (arg == "--row-bound" && i + 1 < argc) {
            row_bound = std::atoi(argv[++i]);
        } else if (arg == "--progress") {
            progress = true;
        } else {
            std::cerr
                << "usage: " << argv[0]
                << " [--q-min N] [--q-max N] [--offset N]"
                << " [--row-bound N] [--progress]\n";
            std::cout << "__EXIT_STATUS=2\n";
            return 2;
        }
    }

    if (q_min < 1 || q_max < q_min || offset < 0 || row_bound < 1) {
        std::cerr << "invalid complement-packing verifier arguments\n";
        std::cout << "__EXIT_STATUS=2\n";
        return 2;
    }

    const int k_min = ceil_div_3(q_min) + offset;
    const int k_max = ceil_div_3(q_max) + offset;
    char hostname[256] = {};
    if (gethostname(hostname, sizeof(hostname) - 1) != 0) {
        hostname[0] = '\0';
    }

    std::cout << "B/C residual exact complement-packing GMP verifier\n"
              << "VERIFIES exact inequality: row_bound * M_{ceil(q/3)+offset} < 2^(2q)\n"
              << "where M_K=sum_{d=0}^K (d+1)*s_d and "
              << "s_{n+1}=n*s_n+n*s_{n-1}-binom(n,2)*s_{n-2}\n"
              << "host=" << (hostname[0] == '\0' ? "unknown" : hostname) << "\n"
              << "q_min=" << q_min
              << " q_max=" << q_max
              << " offset=" << offset
              << " row_bound=" << row_bound
              << " k_min=" << k_min
              << " k_max=" << k_max
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;

    const std::vector<mpz_class> moments = stable_moments(k_max, progress);
    const std::vector<mpz_class> cumulative = cumulative_counts(moments);

    std::vector<Check> checks(static_cast<std::size_t>(q_max - q_min + 1));
    int completed = 0;
    int failures = 0;

#pragma omp parallel for schedule(dynamic, 1) reduction(+:failures)
    for (int q = q_min; q <= q_max; ++q) {
        const int slot = q - q_min;
        const int k = ceil_div_3(q) + offset;
        mpz_class lhs = row_bound * cumulative[static_cast<std::size_t>(k)];
        mpz_class rhs;
        mpz_ui_pow_ui(rhs.get_mpz_t(), 2UL, static_cast<unsigned long>(2 * q));
        const bool ok = lhs < rhs;
        if (!ok) ++failures;
        const mpz_class margin = rhs - lhs;
        Check check;
        check.q = q;
        check.k = k;
        check.lhs_bits = bit_length(lhs);
        check.rhs_bits = bit_length(rhs);
        check.margin_bits = ok ? bit_length(margin) : 0;
        check.ratio_log10 = log10_abs_mpz(lhs) - log10_abs_mpz(rhs);
        check.ok = ok;
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

    Check worst = checks.front();
    for (const Check& check : checks) {
        if (check.ratio_log10 > worst.ratio_log10) worst = check;
    }

    for (const Check& check : checks) {
        std::cout << "check q=" << check.q
                  << " K=" << check.k
                  << " lhs_bits=" << check.lhs_bits
                  << " rhs_bits=" << check.rhs_bits
                  << " margin_bits=" << check.margin_bits
                  << " log10_lhs_over_rhs=" << display(check.ratio_log10)
                  << " status=" << (check.ok ? "OK" : "FAIL")
                  << std::endl;
    }

    std::cout << "anchor_value M_39=" << cumulative[39] << std::endl;
    std::cout << "anchor_value M_40=" << cumulative[40] << std::endl;
    std::cout << "anchor_value M_83=" << cumulative[83] << std::endl;

    std::cout << "SUMMARY failures=" << failures
              << " q_checked=" << checks.size()
              << " k_min=" << k_min
              << " k_max=" << k_max
              << " worst_q=" << worst.q
              << " worst_K=" << worst.k
              << " worst_log10_lhs_over_rhs=" << display(worst.ratio_log10)
              << " worst_lhs_bits=" << worst.lhs_bits
              << " worst_rhs_bits=" << worst.rhs_bits
              << " worst_margin_bits=" << worst.margin_bits
              << std::endl;
    std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << std::endl;

    return failures == 0 ? 0 : 1;
}
