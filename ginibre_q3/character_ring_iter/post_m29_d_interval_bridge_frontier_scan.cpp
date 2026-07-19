#include <omp.h>

#include <algorithm>
#include <cctype>
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

constexpr long double NEG_INF = -std::numeric_limits<long double>::infinity();
constexpr long double LOG10_2 = 0.301029995663981195213738894724493027L;
constexpr long double LOG10_8 = 0.903089986991943585641216684173479084L;

struct SignedLog {
    int sign = 0;
    long double log_abs = NEG_INF;
};

struct FrontierResult {
    int chain_m = 0;
    int odd_n = 0;
    int frontier_rank = 0;
    int active_rows = 0;
    SignedLog stable;
    SignedLog negative_linear;
    SignedLog negative_quadratic;
    SignedLog lower;
};

bool is_zero(const SignedLog& value) {
    return value.sign == 0;
}

void add_signed_log(SignedLog& acc, int sign, long double log_abs) {
    if (sign == 0 || !std::isfinite(log_abs)) return;
    if (acc.sign == 0) {
        acc.sign = sign > 0 ? 1 : -1;
        acc.log_abs = log_abs;
        return;
    }
    const int term_sign = sign > 0 ? 1 : -1;
    const long double hi = std::max(acc.log_abs, log_abs);
    const long double lo = std::min(acc.log_abs, log_abs);
    if (acc.sign == term_sign) {
        if (hi - lo > 80.0L) {
            acc.log_abs = hi;
        } else {
            acc.log_abs = hi + std::log10(1.0L + std::pow(10.0L, lo - hi));
        }
        return;
    }
    if (std::fabs(acc.log_abs - log_abs) < 1.0e-18L) {
        acc.sign = 0;
        acc.log_abs = NEG_INF;
        return;
    }
    const int hi_sign = (acc.log_abs > log_abs) ? acc.sign : term_sign;
    if (hi - lo > 80.0L) {
        acc.sign = hi_sign;
        acc.log_abs = hi;
    } else {
        const long double diff = 1.0L - std::pow(10.0L, lo - hi);
        if (!(diff > 0.0L)) {
            acc.sign = 0;
            acc.log_abs = NEG_INF;
        } else {
            acc.sign = hi_sign;
            acc.log_abs = hi + std::log10(diff);
        }
    }
}

SignedLog add_values(const SignedLog& first, const SignedLog& second) {
    SignedLog out = first;
    add_signed_log(out, second.sign, second.log_abs);
    return out;
}

std::string display_signed_log(const SignedLog& value) {
    if (value.sign == 0) return "0";
    std::ostringstream out;
    out << (value.sign > 0 ? "+" : "-") << std::setprecision(18)
        << static_cast<double>(value.log_abs);
    return out.str();
}

std::vector<long double> log10_factorials(int max_n) {
    std::vector<long double> out(static_cast<std::size_t>(max_n + 1), 0.0L);
    for (int n = 2; n <= max_n; ++n) {
        out[static_cast<std::size_t>(n)] =
            out[static_cast<std::size_t>(n - 1)]
            + std::log10(static_cast<long double>(n));
    }
    return out;
}

long double log10_binom(
    const std::vector<long double>& log_fact,
    int n,
    int k
) {
    if (k < 0 || k > n) return NEG_INF;
    return log_fact[static_cast<std::size_t>(n)]
        - log_fact[static_cast<std::size_t>(k)]
        - log_fact[static_cast<std::size_t>(n - k)];
}

std::vector<SignedLog> stable_moment_logs(int max_index) {
    std::vector<SignedLog> moments(static_cast<std::size_t>(max_index + 1));
    moments[0] = SignedLog{1, 0.0L};
    if (max_index >= 1) moments[1] = SignedLog{0, NEG_INF};
    for (int index = 1; index < max_index; ++index) {
        SignedLog next;
        if (!is_zero(moments[static_cast<std::size_t>(index)])) {
            add_signed_log(
                next,
                moments[static_cast<std::size_t>(index)].sign,
                std::log10(static_cast<long double>(index))
                    + moments[static_cast<std::size_t>(index)].log_abs
            );
        }
        if (!is_zero(moments[static_cast<std::size_t>(index - 1)])) {
            add_signed_log(
                next,
                moments[static_cast<std::size_t>(index - 1)].sign,
                std::log10(static_cast<long double>(index))
                    + moments[static_cast<std::size_t>(index - 1)].log_abs
            );
        }
        if (index >= 2 && !is_zero(moments[static_cast<std::size_t>(index - 2)])) {
            const long double factor =
                std::log10(static_cast<long double>(index))
                + std::log10(static_cast<long double>(index - 1))
                - LOG10_2;
            add_signed_log(
                next,
                -moments[static_cast<std::size_t>(index - 2)].sign,
                factor + moments[static_cast<std::size_t>(index - 2)].log_abs
            );
        }
        if (next.sign <= 0) {
            std::cerr << "stable moment recurrence lost positivity at index "
                      << (index + 1) << "\n";
            std::exit(1);
        }
        moments[static_cast<std::size_t>(index + 1)] = next;
    }
    return moments;
}

SignedLog q3_pair_coefficient(
    int n_value,
    int first,
    int second,
    int scale_sign,
    long double scale_log_abs,
    const std::vector<long double>& log_fact
) {
    SignedLog coeff;
    const int pos_ks[2] = {first - 2, second - 2};
    for (int slot = 0; slot < 2; ++slot) {
        if (slot == 1 && pos_ks[1] == pos_ks[0]) continue;
        const int k = pos_ks[slot];
        if (k >= 0 && k <= n_value) {
            add_signed_log(
                coeff,
                scale_sign,
                scale_log_abs + log10_binom(log_fact, n_value, k)
            );
        }
    }
    const int neg_ks[2] = {first - 1, second - 1};
    for (int slot = 0; slot < 2; ++slot) {
        if (slot == 1 && neg_ks[1] == neg_ks[0]) continue;
        const int k = neg_ks[slot];
        if (k >= 0 && k <= n_value) {
            add_signed_log(
                coeff,
                -scale_sign,
                scale_log_abs + log10_binom(log_fact, n_value, k)
            );
        }
    }
    return coeff;
}

void accumulate_diagonal(
    int n_value,
    int scale_sign,
    long double scale_log_abs,
    int frontier_rank,
    const std::vector<long double>& log_fact,
    const std::vector<SignedLog>& stable_moments,
    std::vector<SignedLog>& linear_coeffs,
    SignedLog& stable_value,
    SignedLog& negative_quadratic
) {
    const int pair_sum = n_value + 2;
    for (int first = 0; first <= pair_sum / 2; ++first) {
        const int second = pair_sum - first;
        const SignedLog coeff = q3_pair_coefficient(
            n_value,
            first,
            second,
            scale_sign,
            scale_log_abs,
            log_fact
        );
        if (coeff.sign == 0) continue;
        const SignedLog& m_first = stable_moments[static_cast<std::size_t>(first)];
        const SignedLog& m_second = stable_moments[static_cast<std::size_t>(second)];
        if (m_first.sign == 0 || m_second.sign == 0) continue;
        const long double pair_log =
            coeff.log_abs + m_first.log_abs + m_second.log_abs;
        add_signed_log(stable_value, coeff.sign, pair_log);
        if (coeff.sign < 0 && first >= frontier_rank && second >= frontier_rank) {
            add_signed_log(negative_quadratic, coeff.sign, pair_log);
        }
        if (first == second) {
            add_signed_log(
                linear_coeffs[static_cast<std::size_t>(first)],
                coeff.sign,
                coeff.log_abs + LOG10_2 + m_first.log_abs
            );
        } else {
            add_signed_log(
                linear_coeffs[static_cast<std::size_t>(first)],
                coeff.sign,
                coeff.log_abs + m_second.log_abs
            );
            add_signed_log(
                linear_coeffs[static_cast<std::size_t>(second)],
                coeff.sign,
                coeff.log_abs + m_first.log_abs
            );
        }
    }
}

FrontierResult scan_chain_m(
    int chain_m,
    int frontier_rank,
    int active_rows,
    const std::vector<long double>& log_fact,
    const std::vector<SignedLog>& stable_moments
) {
    const int max_index = 2 * chain_m + 5;
    std::vector<SignedLog> linear_coeffs(static_cast<std::size_t>(max_index + 1));
    FrontierResult result;
    result.chain_m = chain_m;
    result.odd_n = 2 * chain_m + 3;
    result.frontier_rank = frontier_rank;
    result.active_rows = active_rows;
    accumulate_diagonal(
        2 * chain_m + 3,
        1,
        LOG10_2,
        frontier_rank,
        log_fact,
        stable_moments,
        linear_coeffs,
        result.stable,
        result.negative_quadratic
    );
    accumulate_diagonal(
        2 * chain_m + 1,
        -1,
        LOG10_8,
        frontier_rank,
        log_fact,
        stable_moments,
        linear_coeffs,
        result.stable,
        result.negative_quadratic
    );
    for (int index = frontier_rank; index <= max_index; ++index) {
        const SignedLog& coeff = linear_coeffs[static_cast<std::size_t>(index)];
        const SignedLog& moment = stable_moments[static_cast<std::size_t>(index)];
        if (coeff.sign < 0 && moment.sign > 0) {
            add_signed_log(
                result.negative_linear,
                coeff.sign,
                coeff.log_abs + moment.log_abs
            );
        }
    }
    result.lower = result.stable;
    result.lower = add_values(result.lower, result.negative_linear);
    result.lower = add_values(result.lower, result.negative_quadratic);
    return result;
}

std::vector<int> read_d_onsets(
    const std::string& path,
    int rank_lo,
    int rank_hi
) {
    std::vector<int> onsets(static_cast<std::size_t>(rank_hi + 1), 0);
    std::ifstream input(path);
    if (!input) {
        std::cerr << "could not open onset log: " << path << "\n";
        std::exit(2);
    }
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("D_", 0) != 0) continue;
        std::istringstream in(line);
        std::string label;
        int onset = 0;
        in >> label >> onset;
        if (!in || label.size() < 3) continue;
        const int rank = std::atoi(label.c_str() + 2);
        if (rank >= rank_lo && rank <= rank_hi) {
            onsets[static_cast<std::size_t>(rank)] = onset;
        }
    }
    for (int rank = rank_lo; rank <= rank_hi; ++rank) {
        const int onset = onsets[static_cast<std::size_t>(rank)];
        if (onset <= 0 || onset % 2 == 0) {
            std::cerr << "missing odd onset for D_" << rank << "\n";
            std::exit(2);
        }
    }
    return onsets;
}

int frontier_for_chain_m(
    int chain_m,
    int rank_lo,
    int rank_hi,
    const std::vector<int>& bridge_max,
    int& active_rows
) {
    int frontier = 0;
    active_rows = 0;
    for (int rank = rank_lo; rank <= rank_hi; ++rank) {
        if (bridge_max[static_cast<std::size_t>(rank)] >= chain_m) {
            ++active_rows;
            if (frontier == 0 || rank < frontier) frontier = rank;
        }
    }
    return frontier;
}

}  // namespace

int main(int argc, char** argv) {
    int rank_lo = 98;
    int rank_hi = 279;
    int chain_lo = 31;
    int chain_hi_override = 0;
    int progress_step = 250;
    bool progress = false;
    std::string onset_log;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--rank-lo" && i + 1 < argc) {
            rank_lo = std::atoi(argv[++i]);
        } else if (arg == "--rank-hi" && i + 1 < argc) {
            rank_hi = std::atoi(argv[++i]);
        } else if (arg == "--chain-lo" && i + 1 < argc) {
            chain_lo = std::atoi(argv[++i]);
        } else if (arg == "--chain-hi" && i + 1 < argc) {
            chain_hi_override = std::atoi(argv[++i]);
        } else if (arg == "--onset-log" && i + 1 < argc) {
            onset_log = argv[++i];
        } else if (arg == "--progress-step" && i + 1 < argc) {
            progress_step = std::max(1, std::atoi(argv[++i]));
        } else if (arg == "--progress") {
            progress = true;
        } else {
            std::cerr
                << "usage: " << argv[0]
                << " --onset-log PATH [--rank-lo N] [--rank-hi N]"
                << " [--chain-lo M] [--chain-hi M] [--progress]"
                << " [--progress-step N]\n";
            return 2;
        }
    }
    if (onset_log.empty() || rank_lo < 1 || rank_hi < rank_lo || chain_lo < 0) {
        std::cerr << "invalid D bridge frontier scan arguments\n";
        return 2;
    }

    const std::vector<int> onsets = read_d_onsets(onset_log, rank_lo, rank_hi);
    std::vector<int> bridge_max(static_cast<std::size_t>(rank_hi + 1), 0);
    int chain_hi = 0;
    for (int rank = rank_lo; rank <= rank_hi; ++rank) {
        bridge_max[static_cast<std::size_t>(rank)] =
            (onsets[static_cast<std::size_t>(rank)] - 3) / 2;
        chain_hi = std::max(chain_hi, bridge_max[static_cast<std::size_t>(rank)]);
    }
    if (chain_hi_override > 0) chain_hi = std::min(chain_hi, chain_hi_override);
    if (chain_hi < chain_lo) {
        std::cerr << "empty chain interval\n";
        return 2;
    }

    const int max_index = 2 * chain_hi + 5;
    std::cout << "D interval stable-envelope bridge frontier diagnostic\n"
              << "rank_lo=" << rank_lo << " rank_hi=" << rank_hi
              << " chain_lo=" << chain_lo << " chain_hi=" << chain_hi
              << " max_index=" << max_index
              << " onset_log=" << onset_log
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;

    const std::vector<long double> log_fact = log10_factorials(max_index);
    const std::vector<SignedLog> stable_moments = stable_moment_logs(max_index);

    const int total = chain_hi - chain_lo + 1;
    std::vector<FrontierResult> results(static_cast<std::size_t>(total));
    int completed = 0;
#pragma omp parallel for schedule(dynamic, 8)
    for (int offset = 0; offset < total; ++offset) {
        const int chain_m = chain_lo + offset;
        int active_rows = 0;
        const int frontier = frontier_for_chain_m(
            chain_m,
            rank_lo,
            rank_hi,
            bridge_max,
            active_rows
        );
        if (frontier != 0) {
            results[static_cast<std::size_t>(offset)] = scan_chain_m(
                chain_m,
                frontier,
                active_rows,
                log_fact,
                stable_moments
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
            }
        }
    }

    int failures = 0;
    FrontierResult worst;
    bool have_worst = false;
    for (const FrontierResult& result : results) {
        if (result.frontier_rank == 0) continue;
        if (result.lower.sign <= 0) ++failures;
        if (result.lower.sign > 0) {
            if (!have_worst) {
                worst = result;
                have_worst = true;
            } else {
                const long double worst_relative =
                    worst.lower.log_abs - worst.stable.log_abs;
                const long double result_relative =
                    result.lower.log_abs - result.stable.log_abs;
                if (result_relative < worst_relative) worst = result;
            }
        }
    }

    std::cout << std::setprecision(18)
              << "worst_positive\tm\todd_n\tfrontier\tactive_rows"
              << "\tlower_log10\tstable_log10\trelative_log10"
              << "\tnegative_linear_log10\tnegative_quadratic_log10\n";
    if (have_worst) {
        std::cout << "worst_positive\t" << worst.chain_m
                  << '\t' << worst.odd_n
                  << "\tD_" << worst.frontier_rank
                  << '\t' << worst.active_rows
                  << '\t' << static_cast<double>(worst.lower.log_abs)
                  << '\t' << static_cast<double>(worst.stable.log_abs)
                  << '\t' << static_cast<double>(worst.lower.log_abs - worst.stable.log_abs)
                  << '\t' << display_signed_log(worst.negative_linear)
                  << '\t' << display_signed_log(worst.negative_quadratic)
                  << '\n';
    }

    std::cout << "sample\tm\todd_n\tfrontier\tactive_rows"
              << "\tlower_log10\tstable_log10\trelative_log10"
              << "\tnegative_linear_log10\tnegative_quadratic_log10\n";
    const int sample_stride = std::max(1, total / 20);
    for (int offset = 0; offset < total; offset += sample_stride) {
        const FrontierResult& result = results[static_cast<std::size_t>(offset)];
        if (result.frontier_rank == 0) continue;
        std::cout << "sample\t" << result.chain_m
                  << '\t' << result.odd_n
                  << "\tD_" << result.frontier_rank
                  << '\t' << result.active_rows
                  << '\t' << display_signed_log(result.lower)
                  << '\t' << display_signed_log(result.stable)
                  << '\t';
        if (result.lower.sign > 0 && result.stable.sign > 0) {
            std::cout << static_cast<double>(result.lower.log_abs - result.stable.log_abs);
        } else {
            std::cout << "nan";
        }
        std::cout << '\t' << display_signed_log(result.negative_linear)
                  << '\t' << display_signed_log(result.negative_quadratic)
                  << '\n';
    }

    std::cout << "SUMMARY failures=" << failures;
    if (have_worst) {
        std::cout << " worst_m=" << worst.chain_m
                  << " worst_odd_n=" << worst.odd_n
                  << " worst_frontier=D_" << worst.frontier_rank
                  << " worst_relative_log10="
                  << static_cast<double>(worst.lower.log_abs - worst.stable.log_abs)
                  << " worst_lower_log10=" << static_cast<double>(worst.lower.log_abs);
    }
    std::cout << std::endl;
    return failures == 0 ? 0 : 1;
}
