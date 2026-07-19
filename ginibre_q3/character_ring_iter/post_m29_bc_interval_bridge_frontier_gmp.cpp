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

std::vector<mpz_class> binom_row(int n_value) {
    std::vector<mpz_class> row(static_cast<std::size_t>(n_value + 1));
    row[0] = 1;
    for (int k = 0; k < n_value; ++k) {
        row[static_cast<std::size_t>(k + 1)] =
            row[static_cast<std::size_t>(k)] * (n_value - k) / (k + 1);
    }
    return row;
}

void add_binom_if_valid(
    mpz_class& coeff,
    const std::vector<mpz_class>& row,
    int k,
    int factor
) {
    if (k < 0 || k >= static_cast<int>(row.size())) return;
    if (factor == 1) {
        coeff += row[static_cast<std::size_t>(k)];
    } else if (factor == -1) {
        coeff -= row[static_cast<std::size_t>(k)];
    } else {
        coeff += factor * row[static_cast<std::size_t>(k)];
    }
}

mpz_class pair_coefficient(
    const std::vector<mpz_class>& row,
    int first,
    int second,
    int scale
) {
    const int factor = 2 * scale;
    mpz_class coeff = 0;
    const int pos_ks[2] = {first - 2, second - 2};
    for (int slot = 0; slot < 2; ++slot) {
        if (slot == 1 && pos_ks[1] == pos_ks[0]) continue;
        add_binom_if_valid(coeff, row, pos_ks[slot], factor);
    }
    const int neg_ks[2] = {first - 1, second - 1};
    for (int slot = 0; slot < 2; ++slot) {
        if (slot == 1 && neg_ks[1] == neg_ks[0]) continue;
        add_binom_if_valid(coeff, row, neg_ks[slot], -factor);
    }
    return coeff;
}

void accumulate_diagonal(
    int n_value,
    int scale,
    int frontier_boundary,
    const std::vector<mpz_class>& stable,
    std::vector<mpz_class>& linear_coeffs,
    mpz_class& stable_value,
    mpz_class& negative_quadratic
) {
    const std::vector<mpz_class> row = binom_row(n_value);
    const int pair_sum = n_value + 2;
    for (int first = 0; first <= pair_sum / 2; ++first) {
        const int second = pair_sum - first;
        mpz_class coeff = pair_coefficient(row, first, second, scale);
        if (sgn(coeff) == 0) continue;

        mpz_class pair_value =
            coeff * stable[static_cast<std::size_t>(first)]
            * stable[static_cast<std::size_t>(second)];
        stable_value += pair_value;
        if (sgn(coeff) < 0 && first >= frontier_boundary && second >= frontier_boundary) {
            negative_quadratic += pair_value;
        }

        if (first == second) {
            linear_coeffs[static_cast<std::size_t>(first)] +=
                2 * coeff * stable[static_cast<std::size_t>(first)];
        } else {
            linear_coeffs[static_cast<std::size_t>(first)] +=
                coeff * stable[static_cast<std::size_t>(second)];
            linear_coeffs[static_cast<std::size_t>(second)] +=
                coeff * stable[static_cast<std::size_t>(first)];
        }
    }
}

ChainResult verify_chain_m_exact(
    int chain_m,
    int frontier_boundary,
    int active_rows,
    const std::vector<mpz_class>& stable,
    int lambda_num,
    int lambda_den
) {
    const int max_index = 2 * chain_m + 5;
    std::vector<mpz_class> linear_coeffs(static_cast<std::size_t>(max_index + 1));
    mpz_class stable_value = 0;
    mpz_class negative_quadratic = 0;
    accumulate_diagonal(
        2 * chain_m + 3,
        1,
        frontier_boundary,
        stable,
        linear_coeffs,
        stable_value,
        negative_quadratic
    );
    accumulate_diagonal(
        2 * chain_m + 1,
        -4,
        frontier_boundary,
        stable,
        linear_coeffs,
        stable_value,
        negative_quadratic
    );

    mpz_class negative_linear = 0;
    for (int index = frontier_boundary; index <= max_index; ++index) {
        const mpz_class& coeff = linear_coeffs[static_cast<std::size_t>(index)];
        if (sgn(coeff) > 0) {
            negative_linear -= coeff * stable[static_cast<std::size_t>(index)];
        }
    }
    const mpz_class positive_linear_abs = -negative_linear;
    const mpz_class negative_quadratic_abs = -negative_quadratic;
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
              << "frontier lower bound uses |delta_j|<=lambda*s_j with positive linear "
              << "and negative quadratic coefficients"
              << std::endl;

    std::cout << "stable_moments_start max_index=" << max_index << std::endl;
    const std::vector<mpz_class> stable = stable_moments(max_index);
    std::cout << "stable_moments_done max_index=" << max_index << std::endl;

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
              << "\tpositive_linear_log10\tnegative_quadratic_log10"
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
