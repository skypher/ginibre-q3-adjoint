#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double LOG10_2_DOUBLE = 0.301029995663981195213738894724493027;

struct Row {
    char family = 'B';
    int rank = 0;
};

struct ChainResult {
    int chain_m = 0;
    int odd_n = 0;
    int frontier_boundary = 0;
    int active_rows = 0;
    int lower_sign = 0;
    double lower_log10 = -std::numeric_limits<double>::infinity();
    double stable_log10 = -std::numeric_limits<double>::infinity();
    double positive_linear_log10 = -std::numeric_limits<double>::infinity();
    double negative_quadratic_log10 = -std::numeric_limits<double>::infinity();
    double uniform_lambda_log10 = std::numeric_limits<double>::quiet_NaN();
    double relative_log10 = std::numeric_limits<double>::quiet_NaN();
};

double log10_abs_mpz(const mpz_class& value) {
    if (sgn(value) == 0) return -std::numeric_limits<double>::infinity();
    mpz_class abs_value = value >= 0 ? value : -value;
    long exponent = 0;
    const double mantissa = mpz_get_d_2exp(&exponent, abs_value.get_mpz_t());
    return std::log10(mantissa) + static_cast<double>(exponent) * LOG10_2_DOUBLE;
}

double log10_sum(double first, double second) {
    if (!std::isfinite(first)) return second;
    if (!std::isfinite(second)) return first;
    const double base = std::max(first, second);
    return base + std::log10(
        std::pow(10.0, first - base) + std::pow(10.0, second - base)
    );
}

double uniform_lambda_threshold_log10(
    const mpz_class& stable_value,
    const mpz_class& positive_linear_abs,
    const mpz_class& negative_quadratic_abs
) {
    if (sgn(stable_value) <= 0) return std::numeric_limits<double>::quiet_NaN();
    const bool have_linear = sgn(positive_linear_abs) > 0;
    const bool have_quadratic = sgn(negative_quadratic_abs) > 0;
    if (!have_linear && !have_quadratic) {
        return std::numeric_limits<double>::infinity();
    }
    const double log_stable = log10_abs_mpz(stable_value);
    if (!have_quadratic) {
        return log_stable - log10_abs_mpz(positive_linear_abs);
    }
    if (!have_linear) {
        return 0.5 * (log_stable - log10_abs_mpz(negative_quadratic_abs));
    }

    const double log_linear = log10_abs_mpz(positive_linear_abs);
    const double log_quadratic = log10_abs_mpz(negative_quadratic_abs);
    const double log_sqrt = 0.5 * log10_sum(
        2.0 * log_linear,
        std::log10(4.0) + log_quadratic + log_stable
    );
    const double log_denom = log10_sum(log_linear, log_sqrt);
    return std::log10(2.0) + log_stable - log_denom;
}

std::vector<mpz_class> stable_moments(int max_index) {
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
            std::exit(1);
        }
        moments[static_cast<std::size_t>(index + 1)] = std::move(value);
    }
    return moments;
}

std::vector<mpz_class> stable_chain_values(int max_n) {
    std::vector<mpz_class> values(static_cast<std::size_t>(max_n + 1));
    values[0] = 2;
    for (int n = 0; n < max_n; ++n) {
        const long n1 = n;
        const long n2 = n1 * (n - 1);
        const long n3 = n2 * (n - 2);
        const long n4 = n3 * (n - 3);
        mpz_class numerator = (4 * n1 + 2) * values[static_cast<std::size_t>(n)];
        if (n >= 1) {
            numerator += (4 * n1 - 3 * n2)
                * values[static_cast<std::size_t>(n - 1)];
        }
        if (n >= 2) {
            numerator += (n3 - 6 * n2)
                * values[static_cast<std::size_t>(n - 2)];
        }
        if (n >= 3) {
            numerator += 4 * n3 * values[static_cast<std::size_t>(n - 3)];
        }
        if (n >= 4) {
            numerator -= n4 * values[static_cast<std::size_t>(n - 4)];
        }
        if (!mpz_even_p(numerator.get_mpz_t())) {
            std::cerr << "stable chain recurrence lost integrality at n=" << n << "\n";
            std::exit(1);
        }
        values[static_cast<std::size_t>(n + 1)] = numerator / 2;
        if (sgn(values[static_cast<std::size_t>(n + 1)]) <= 0) {
            std::cerr << "stable chain recurrence lost positivity at n=" << (n + 1) << "\n";
            std::exit(1);
        }
    }
    return values;
}

long long a_sign_polynomial(int n_value, int k) {
    const long long pair_sum = static_cast<long long>(n_value) + 2;
    const long long difference = 2LL * k - pair_sum;
    return difference * difference - pair_sum;
}

mpz_class a_term(
    int n_value,
    int k,
    const std::vector<mpz_class>& stable
) {
    const int pair_sum = n_value + 2;
    const int reflected = pair_sum - k;
    const int t = std::min(k, reflected);
    mpz_class coefficient;
    if (t == 0) {
        coefficient = 1;
    } else if (t == 1) {
        coefficient = n_value - 2;
    } else {
        mpz_bin_uiui(
            coefficient.get_mpz_t(),
            static_cast<unsigned long>(n_value),
            static_cast<unsigned long>(t - 2)
        );
        coefficient *= static_cast<long>(a_sign_polynomial(n_value, k));
        const long denominator = static_cast<long>(t) * (t - 1);
        if (!mpz_divisible_ui_p(
                coefficient.get_mpz_t(),
                static_cast<unsigned long>(denominator))) {
            std::cerr << "nonintegral A coefficient at n=" << n_value
                      << " k=" << k << "\n";
            std::exit(1);
        }
        coefficient /= denominator;
    }
    return coefficient * stable[static_cast<std::size_t>(k)]
        * stable[static_cast<std::size_t>(reflected)];
}

template <typename Predicate>
mpz_class central_negative_sum(
    int n_value,
    const std::vector<mpz_class>& stable,
    Predicate include
) {
    const int pair_sum = n_value + 2;
    const int radius = static_cast<int>(std::sqrt(static_cast<double>(pair_sum))) + 2;
    const int center = pair_sum / 2;
    const int low = std::max(0, center - radius);
    const int high = std::min(pair_sum, center + radius + 1);
    mpz_class sum = 0;
    for (int k = low; k <= high; ++k) {
        if (a_sign_polynomial(n_value, k) < 0 && include(k, pair_sum - k)) {
            sum += a_term(n_value, k, stable);
        }
    }
    if (sgn(sum) > 0) {
        std::cerr << "negative A sum has wrong sign at n=" << n_value << "\n";
        std::exit(1);
    }
    return sum;
}

mpz_class central_negative_sum_all(
    int n_value,
    const std::vector<mpz_class>& stable
) {
    const int pair_sum = n_value + 2;
    const int center = pair_sum / 2;
    const int radius = static_cast<int>(std::sqrt(static_cast<double>(pair_sum))) + 2;
    const int low = std::max(0, center - radius);
    mpz_class sum = 0;
    for (int k = low; k <= center; ++k) {
        if (a_sign_polynomial(n_value, k) >= 0) continue;
        mpz_class term = a_term(n_value, k, stable);
        sum += (2 * k == pair_sum) ? term : 2 * term;
    }
    if (sgn(sum) > 0) {
        std::cerr << "negative A sum has wrong sign at n=" << n_value << "\n";
        std::exit(1);
    }
    return sum;
}

bool central_support_is_active(int n_value, int frontier_boundary) {
    const int pair_sum = n_value + 2;
    const int center = pair_sum / 2;
    const int radius = static_cast<int>(std::sqrt(static_cast<double>(pair_sum))) + 2;
    const int low = std::max(0, center - radius);
    const int high = std::min(pair_sum, center + radius + 1);
    for (int k = low; k <= high; ++k) {
        if (a_sign_polynomial(n_value, k) < 0
            && (k < frontier_boundary || pair_sum - k < frontier_boundary)) {
            return false;
        }
    }
    return true;
}

mpz_class positive_boundary_sum(
    int n_value,
    int frontier_boundary,
    const std::vector<mpz_class>& stable
) {
    if (frontier_boundary <= 0 || frontier_boundary > n_value + 2) {
        std::cerr << "invalid positive-boundary range at n=" << n_value << "\n";
        std::exit(1);
    }
    mpz_class sum = 0;
    mpz_class binomial = 1;
    for (int k = 0; k < frontier_boundary; ++k) {
        if (a_sign_polynomial(n_value, k) > 0) {
            mpz_class coefficient;
            if (k == 0) {
                coefficient = 1;
            } else if (k == 1) {
                coefficient = n_value - 2;
            } else {
                coefficient = binomial * static_cast<long>(a_sign_polynomial(n_value, k));
                const long denominator = static_cast<long>(k) * (k - 1);
                if (!mpz_divisible_ui_p(
                        coefficient.get_mpz_t(),
                        static_cast<unsigned long>(denominator))) {
                    std::cerr << "nonintegral boundary coefficient at n=" << n_value
                              << " k=" << k << "\n";
                    std::exit(1);
                }
                coefficient /= denominator;
            }
            sum += coefficient * stable[static_cast<std::size_t>(k)]
                * stable[static_cast<std::size_t>(n_value + 2 - k)];
        }
        if (k >= 2 && k + 1 < frontier_boundary) {
            binomial *= n_value - k + 2;
            binomial /= k - 1;
        }
    }
    if (sgn(sum) < 0) {
        std::cerr << "positive boundary sum has wrong sign at n=" << n_value << "\n";
        std::exit(1);
    }
    return sum;
}

mpz_class binomial_or_zero(int n_value, int k) {
    if (k < 0 || k > n_value) return 0;
    mpz_class value;
    mpz_bin_uiui(
        value.get_mpz_t(),
        static_cast<unsigned long>(n_value),
        static_cast<unsigned long>(k)
    );
    return value;
}

void verify_algebraic_reduction(
    int max_n,
    const std::vector<mpz_class>& stable,
    const std::vector<mpz_class>& stable_chain
) {
    for (int n_value = 0; n_value <= max_n; ++n_value) {
        const int pair_sum = n_value + 2;
        mpz_class direct_sum = 0;
        for (int k = 0; k <= pair_sum; ++k) {
            const mpz_class coefficient =
                binomial_or_zero(n_value, k - 2)
                + binomial_or_zero(n_value, k)
                - 2 * binomial_or_zero(n_value, k - 1);
            const mpz_class direct_term =
                coefficient * stable[static_cast<std::size_t>(k)]
                * stable[static_cast<std::size_t>(pair_sum - k)];
            if (direct_term != a_term(n_value, k, stable)) {
                std::cerr << "closed A coefficient mismatch at n=" << n_value
                          << " k=" << k << "\n";
                std::exit(1);
            }
            direct_sum += direct_term;
        }
        if (direct_sum != stable_chain[static_cast<std::size_t>(n_value)]) {
            std::cerr << "stable chain recurrence mismatch at n=" << n_value << "\n";
            std::exit(1);
        }
    }
}

ChainResult verify_chain_m_exact(
    int chain_m,
    int frontier_boundary,
    int active_rows,
    const std::vector<mpz_class>& stable,
    const std::vector<mpz_class>& stable_chain,
    int lambda_num,
    int lambda_den
) {
    const int n_first = 2 * chain_m + 3;
    const int n_second = 2 * chain_m + 1;
    const int pair_sum_second = n_second + 2;
    if (2 * frontier_boundary > pair_sum_second + 1) {
        std::cerr << "overlapping frontier boundary at m=" << chain_m << "\n";
        std::exit(1);
    }

    const mpz_class stable_value = stable_chain[static_cast<std::size_t>(n_first)]
        - 4 * stable_chain[static_cast<std::size_t>(n_second)];
    if (sgn(stable_value) <= 0) {
        std::cerr << "stable chain value lost positivity at m=" << chain_m << "\n";
        std::exit(1);
    }

    const mpz_class negative_first_all =
        central_negative_sum_all(n_first, stable);
    const mpz_class negative_second_all =
        central_negative_sum_all(n_second, stable);
    const mpz_class positive_first_boundary =
        positive_boundary_sum(n_first, frontier_boundary, stable);
    const mpz_class positive_second_boundary =
        positive_boundary_sum(n_second, frontier_boundary, stable);

    const mpz_class positive_first_active =
        stable_chain[static_cast<std::size_t>(n_first)]
        - negative_first_all - positive_first_boundary;
    const mpz_class negative_second_active = central_support_is_active(
        n_second,
        frontier_boundary
    ) ? negative_second_all : central_negative_sum(
            n_second,
            stable,
            [frontier_boundary](int k, int) { return k >= frontier_boundary; }
        );
    if (sgn(positive_first_active) < 0 || sgn(negative_second_active) > 0) {
        std::cerr << "linear majorant component has wrong sign at m=" << chain_m << "\n";
        std::exit(1);
    }
    const mpz_class positive_linear_abs =
        2 * positive_first_active - 8 * negative_second_active;

    const mpz_class negative_first_restricted = central_support_is_active(
        n_first,
        frontier_boundary
    ) ? negative_first_all : central_negative_sum(
            n_first,
            stable,
            [frontier_boundary](int k, int reflected) {
                return k >= frontier_boundary && reflected >= frontier_boundary;
            }
        );
    const mpz_class positive_second_restricted =
        stable_chain[static_cast<std::size_t>(n_second)]
        - negative_second_all - 2 * positive_second_boundary;
    if (sgn(positive_second_restricted) < 0) {
        std::cerr << "quadratic majorant component has wrong sign at m=" << chain_m << "\n";
        std::exit(1);
    }
    const mpz_class negative_quadratic =
        negative_first_restricted - 4 * positive_second_restricted;
    const mpz_class negative_quadratic_abs = -negative_quadratic;
    if (sgn(positive_linear_abs) < 0 || sgn(negative_quadratic_abs) < 0) {
        std::cerr << "penalty has wrong sign at m=" << chain_m << "\n";
        std::exit(1);
    }
    const mpz_class lower_scaled =
        mpz_class(lambda_den) * lambda_den * stable_value
        - mpz_class(lambda_num) * lambda_den * positive_linear_abs
        - mpz_class(lambda_num) * lambda_num * negative_quadratic_abs;

    ChainResult result;
    result.chain_m = chain_m;
    result.odd_n = 2 * chain_m + 3;
    result.frontier_boundary = frontier_boundary;
    result.active_rows = active_rows;
    result.lower_sign = sgn(lower_scaled);
    result.lower_log10 =
        log10_abs_mpz(lower_scaled) - 2.0 * std::log10(static_cast<double>(lambda_den));
    result.stable_log10 = log10_abs_mpz(stable_value);
    result.positive_linear_log10 = log10_abs_mpz(positive_linear_abs);
    result.negative_quadratic_log10 = log10_abs_mpz(negative_quadratic_abs);
    result.uniform_lambda_log10 = uniform_lambda_threshold_log10(
        stable_value,
        positive_linear_abs,
        negative_quadratic_abs
    );
    if (result.lower_sign > 0 && sgn(stable_value) > 0) {
        result.relative_log10 = result.lower_log10 - result.stable_log10;
    }
    return result;
}

std::vector<Row> make_rows(
    bool include_b,
    bool include_c,
    int b_rank_lo,
    int b_rank_hi,
    int c_rank_lo,
    int c_rank_hi
) {
    std::vector<Row> rows;
    if (include_b) {
        for (int rank = b_rank_lo; rank <= b_rank_hi; ++rank) rows.push_back({'B', rank});
    }
    if (include_c) {
        for (int rank = c_rank_lo; rank <= c_rank_hi; ++rank) rows.push_back({'C', rank});
    }
    return rows;
}

int boundary_for(const Row& row) {
    return 2 * row.rank + 2;
}

int row_index(const Row& row, int c_rank_hi) {
    if (row.family == 'B') return row.rank;
    return c_rank_hi + 1 + row.rank;
}

std::vector<int> read_bc_onsets(
    const std::string& path,
    const std::vector<Row>& rows,
    int b_rank_hi,
    int c_rank_hi
) {
    const int size = c_rank_hi + 1 + c_rank_hi + 1;
    std::vector<int> onsets(static_cast<std::size_t>(size), 0);
    std::ifstream input(path);
    if (!input) {
        std::cerr << "could not open onset log: " << path << "\n";
        std::exit(2);
    }

    std::string line;
    while (std::getline(input, line)) {
        if (!(line.rfind("B_", 0) == 0 || line.rfind("C_", 0) == 0)) continue;
        std::istringstream in(line);
        std::string label;
        int onset = 0;
        in >> label >> onset;
        if (!in || label.size() < 3) continue;
        const char family = label[0];
        const int rank = std::atoi(label.c_str() + 2);
        if (family == 'B' && rank >= 0 && rank <= b_rank_hi) {
            onsets[static_cast<std::size_t>(rank)] = onset;
        } else if (family == 'C' && rank >= 0 && rank <= c_rank_hi) {
            onsets[static_cast<std::size_t>(c_rank_hi + 1 + rank)] = onset;
        }
    }

    for (const Row& row : rows) {
        const int onset = onsets[static_cast<std::size_t>(row_index(row, c_rank_hi))];
        if (onset <= 0 || onset % 2 == 0) {
            std::cerr << "missing odd onset for " << row.family << "_" << row.rank << "\n";
            std::exit(2);
        }
    }
    return onsets;
}

int frontier_for_chain_m(
    int chain_m,
    const std::vector<Row>& rows,
    const std::vector<int>& bridge_max,
    int c_rank_hi,
    int& active_rows
) {
    int frontier = 0;
    active_rows = 0;
    for (const Row& row : rows) {
        if (bridge_max[static_cast<std::size_t>(row_index(row, c_rank_hi))] >= chain_m) {
            ++active_rows;
            const int boundary = boundary_for(row);
            if (frontier == 0 || boundary < frontier) frontier = boundary;
        }
    }
    return frontier;
}

}  // namespace

int main(int argc, char** argv) {
    int b_rank_lo = 14;
    int b_rank_hi = 217;
    int c_rank_lo = 20;
    int c_rank_hi = 217;
    int chain_lo = 31;
    int chain_hi_override = 0;
    int progress_step = 50;
    int lambda_num = 1;
    int lambda_den = 1;
    bool include_b = true;
    bool include_c = true;
    bool progress = false;
    std::string onset_log;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--b-rank-lo" && i + 1 < argc) {
            b_rank_lo = std::atoi(argv[++i]);
        } else if (arg == "--b-rank-hi" && i + 1 < argc) {
            b_rank_hi = std::atoi(argv[++i]);
        } else if (arg == "--c-rank-lo" && i + 1 < argc) {
            c_rank_lo = std::atoi(argv[++i]);
        } else if (arg == "--c-rank-hi" && i + 1 < argc) {
            c_rank_hi = std::atoi(argv[++i]);
        } else if (arg == "--only-b") {
            include_b = true;
            include_c = false;
        } else if (arg == "--only-c") {
            include_b = false;
            include_c = true;
        } else if (arg == "--chain-lo" && i + 1 < argc) {
            chain_lo = std::atoi(argv[++i]);
        } else if (arg == "--chain-hi" && i + 1 < argc) {
            chain_hi_override = std::atoi(argv[++i]);
        } else if (arg == "--onset-log" && i + 1 < argc) {
            onset_log = argv[++i];
        } else if (arg == "--progress-step" && i + 1 < argc) {
            progress_step = std::max(1, std::atoi(argv[++i]));
        } else if (arg == "--lambda-num" && i + 1 < argc) {
            lambda_num = std::atoi(argv[++i]);
        } else if (arg == "--lambda-den" && i + 1 < argc) {
            lambda_den = std::atoi(argv[++i]);
        } else if (arg == "--progress") {
            progress = true;
        } else {
            std::cerr
                << "usage: " << argv[0]
                << " --onset-log PATH [--b-rank-lo N] [--b-rank-hi N]"
                << " [--c-rank-lo N] [--c-rank-hi N] [--only-b] [--only-c]"
                << " [--chain-lo M] [--chain-hi M] [--progress]"
                << " [--progress-step N] [--lambda-num A] [--lambda-den B]\n";
            return 2;
        }
    }

    if (onset_log.empty()
        || (!include_b && !include_c)
        || b_rank_lo < 1
        || c_rank_lo < 1
        || b_rank_hi < b_rank_lo
        || c_rank_hi < c_rank_lo
        || chain_lo < 0
        || lambda_num < 0
        || lambda_den <= 0
        || lambda_num > lambda_den) {
        std::cerr << "invalid B/C bridge frontier GMP arguments\n";
        return 2;
    }

    const std::vector<Row> rows =
        make_rows(include_b, include_c, b_rank_lo, b_rank_hi, c_rank_lo, c_rank_hi);
    const std::vector<int> onsets = read_bc_onsets(onset_log, rows, b_rank_hi, c_rank_hi);
    const int bridge_size = c_rank_hi + 1 + c_rank_hi + 1;
    std::vector<int> bridge_max(static_cast<std::size_t>(bridge_size), 0);
    int chain_hi = 0;
    for (const Row& row : rows) {
        const int onset = onsets[static_cast<std::size_t>(row_index(row, c_rank_hi))];
        const int max_m = (onset - 3) / 2;
        bridge_max[static_cast<std::size_t>(row_index(row, c_rank_hi))] = max_m;
        chain_hi = std::max(chain_hi, max_m);
    }
    if (chain_hi_override > 0) chain_hi = std::min(chain_hi, chain_hi_override);
    if (chain_hi < chain_lo) {
        std::cerr << "empty chain interval\n";
        return 2;
    }

    const int max_index = 2 * chain_hi + 5;
    std::cout << "B/C interval stable-envelope bridge frontier GMP verifier\n"
              << "b_rank_lo=" << b_rank_lo << " b_rank_hi=" << b_rank_hi
              << " c_rank_lo=" << c_rank_lo << " c_rank_hi=" << c_rank_hi
              << " include_b=" << include_b << " include_c=" << include_c
              << " chain_lo=" << chain_lo << " chain_hi=" << chain_hi
              << " max_index=" << max_index
              << " lambda=" << lambda_num << "/" << lambda_den
              << " onset_log=" << onset_log
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;
    std::cout << "SIGN_MODEL corrections_delta<=0 active from boundary=2*rank+2; "
              << "frontier lower bound uses |delta_j|<=lambda*s_j with a separable "
              << "positive-linear majorant and exact negative quadratic coefficients"
              << std::endl;

    std::cout << "stable_moments_start max_index=" << max_index << std::endl;
    const std::vector<mpz_class> stable = stable_moments(max_index);
    std::cout << "stable_moments_done max_index=" << max_index << std::endl;
    const std::vector<mpz_class> stable_chain = stable_chain_values(2 * chain_hi + 3);
    const int algebraic_crosscheck_max_n = std::min(128, 2 * chain_hi + 3);
    verify_algebraic_reduction(algebraic_crosscheck_max_n, stable, stable_chain);
    std::cout << "algebraic_reduction_crosscheck_done max_n="
              << algebraic_crosscheck_max_n << std::endl;

    const int total = chain_hi - chain_lo + 1;
    std::vector<ChainResult> results(static_cast<std::size_t>(total));
    int completed = 0;
#pragma omp parallel for schedule(dynamic, 1)
    for (int offset = 0; offset < total; ++offset) {
        const int chain_m = chain_lo + offset;
        int active_rows = 0;
        const int frontier = frontier_for_chain_m(
            chain_m,
            rows,
            bridge_max,
            c_rank_hi,
            active_rows
        );
        if (frontier != 0) {
            results[static_cast<std::size_t>(offset)] = verify_chain_m_exact(
                chain_m,
                frontier,
                active_rows,
                stable,
                stable_chain,
                lambda_num,
                lambda_den
            );
        } else {
            results[static_cast<std::size_t>(offset)].chain_m = chain_m;
            results[static_cast<std::size_t>(offset)].odd_n = 2 * chain_m + 3;
        }

        int done_now = 0;
#pragma omp atomic capture
        done_now = ++completed;
        if (progress && (done_now == total || done_now % progress_step == 0)) {
#pragma omp critical
            {
                std::cout << "progress completed=" << done_now << "/" << total
                          << " last_m=" << chain_m << std::endl;
                std::fflush(stdout);
            }
        }
    }

    int failures = 0;
    bool have_worst = false;
    ChainResult worst;
    bool have_min_lambda = false;
    ChainResult min_lambda;
    for (const ChainResult& result : results) {
        if (result.frontier_boundary == 0) continue;
        if (result.lower_sign <= 0) ++failures;
        if (result.lower_sign > 0) {
            if (!have_worst || result.relative_log10 < worst.relative_log10) {
                worst = result;
                have_worst = true;
            }
        }
        if (std::isfinite(result.uniform_lambda_log10)
            && (!have_min_lambda || result.uniform_lambda_log10 < min_lambda.uniform_lambda_log10)) {
            min_lambda = result;
            have_min_lambda = true;
        }
    }

    std::cout << std::setprecision(18)
              << "worst_positive\tm\todd_n\tfrontier_boundary\tfrontier_rank"
              << "\tactive_rows\tlower_log10\tstable_log10\trelative_log10\n";
    if (have_worst) {
        std::cout << "worst_positive\t" << worst.chain_m
                  << '\t' << worst.odd_n
                  << '\t' << worst.frontier_boundary
                  << '\t' << ((worst.frontier_boundary - 2) / 2)
                  << '\t' << worst.active_rows
                  << '\t' << worst.lower_log10
                  << '\t' << worst.stable_log10
                  << '\t' << worst.relative_log10
                  << '\n';
    }

    std::cout << "sample\tm\todd_n\tfrontier_boundary\tfrontier_rank"
              << "\tactive_rows\tlower_sign\tlower_log10\tstable_log10"
              << "\tpositive_linear_upper_log10\tnegative_quadratic_log10"
              << "\tuniform_lambda_log10\trelative_log10\n";
    const int sample_stride = std::max(1, total / 10);
    for (int offset = 0; offset < total; offset += sample_stride) {
        const ChainResult& result = results[static_cast<std::size_t>(offset)];
        if (result.frontier_boundary == 0) continue;
        std::cout << "sample\t" << result.chain_m
                  << '\t' << result.odd_n
                  << '\t' << result.frontier_boundary
                  << '\t' << ((result.frontier_boundary - 2) / 2)
                  << '\t' << result.active_rows
                  << '\t' << result.lower_sign
                  << '\t' << result.lower_log10
                  << '\t' << result.stable_log10
                  << '\t' << result.positive_linear_log10
                  << '\t' << result.negative_quadratic_log10
                  << '\t' << result.uniform_lambda_log10
                  << '\t' << result.relative_log10
                  << '\n';
    }

    std::cout << "SUMMARY failures=" << failures;
    if (have_worst) {
        std::cout << " worst_m=" << worst.chain_m
                  << " worst_odd_n=" << worst.odd_n
                  << " worst_frontier_boundary=" << worst.frontier_boundary
                  << " worst_frontier_rank=" << ((worst.frontier_boundary - 2) / 2)
                  << " worst_relative_log10=" << worst.relative_log10
                  << " worst_lower_log10=" << worst.lower_log10;
    }
    if (have_min_lambda) {
        std::cout << " min_uniform_lambda_m=" << min_lambda.chain_m
                  << " min_uniform_lambda_odd_n=" << min_lambda.odd_n
                  << " min_uniform_lambda_frontier_boundary=" << min_lambda.frontier_boundary
                  << " min_uniform_lambda_rank=" << ((min_lambda.frontier_boundary - 2) / 2)
                  << " min_uniform_lambda_log10=" << min_lambda.uniform_lambda_log10;
    }
    std::cout << std::endl;
    return failures == 0 ? 0 : 1;
}
