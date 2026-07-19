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

struct ChainResult {
    int chain_m = 0;
    int odd_n = 0;
    int frontier_rank = 0;
    int active_rows = 0;
    int lower_sign = 0;
    double lower_log10 = -std::numeric_limits<double>::infinity();
    double stable_log10 = -std::numeric_limits<double>::infinity();
    double relative_log10 = std::numeric_limits<double>::quiet_NaN();
};

double log10_abs_mpz(const mpz_class& value) {
    if (sgn(value) == 0) return -std::numeric_limits<double>::infinity();
    mpz_class abs_value = value >= 0 ? value : -value;
    long exponent = 0;
    const double mantissa = mpz_get_d_2exp(&exponent, abs_value.get_mpz_t());
    return std::log10(mantissa) + static_cast<double>(exponent) * LOG10_2_DOUBLE;
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
    int frontier_rank,
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
        if (sgn(coeff) < 0 && first >= frontier_rank && second >= frontier_rank) {
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
    int frontier_rank,
    int active_rows,
    const std::vector<mpz_class>& stable
) {
    const int max_index = 2 * chain_m + 5;
    std::vector<mpz_class> linear_coeffs(static_cast<std::size_t>(max_index + 1));
    mpz_class stable_value = 0;
    mpz_class negative_quadratic = 0;
    accumulate_diagonal(
        2 * chain_m + 3,
        1,
        frontier_rank,
        stable,
        linear_coeffs,
        stable_value,
        negative_quadratic
    );
    accumulate_diagonal(
        2 * chain_m + 1,
        -4,
        frontier_rank,
        stable,
        linear_coeffs,
        stable_value,
        negative_quadratic
    );

    mpz_class negative_linear = 0;
    for (int index = frontier_rank; index <= max_index; ++index) {
        const mpz_class& coeff = linear_coeffs[static_cast<std::size_t>(index)];
        if (sgn(coeff) < 0) {
            negative_linear += coeff * stable[static_cast<std::size_t>(index)];
        }
    }
    const mpz_class lower = stable_value + negative_linear + negative_quadratic;

    ChainResult result;
    result.chain_m = chain_m;
    result.odd_n = 2 * chain_m + 3;
    result.frontier_rank = frontier_rank;
    result.active_rows = active_rows;
    result.lower_sign = sgn(lower);
    result.lower_log10 = log10_abs_mpz(lower);
    result.stable_log10 = log10_abs_mpz(stable_value);
    if (result.lower_sign > 0 && sgn(stable_value) > 0) {
        result.relative_log10 = result.lower_log10 - result.stable_log10;
    }
    return result;
}

std::vector<int> read_d_onsets(
    const std::string& path,
    int rank_lo,
    int rank_hi
) {
    std::vector<int> onsets(static_cast<std::size_t>(rank_hi + 1), 0);
    std::vector<bool> seen(static_cast<std::size_t>(rank_hi + 1), false);
    std::ifstream input(path);
    if (!input) {
        std::cerr << "could not open onset log: " << path << "\n";
        std::exit(2);
    }
    std::string line;
    while (std::getline(input, line)) {
        std::string label;
        int onset = 0;
        if (line.rfind("D_", 0) == 0) {
            std::istringstream in(line);
            in >> label >> onset;
            if (!in) continue;
        } else if (line.rfind("row=D_", 0) == 0) {
            const std::size_t rank_end = line.find(' ', 6);
            const std::size_t onset_begin = line.find("n=", rank_end);
            if (rank_end == std::string::npos || onset_begin == std::string::npos) continue;
            label = line.substr(4, rank_end - 4);
            onset = std::atoi(line.c_str() + onset_begin + 2);
        } else {
            continue;
        }
        if (label.size() < 3) continue;
        const int rank = std::atoi(label.c_str() + 2);
        if (rank >= rank_lo && rank <= rank_hi) {
            if (seen[static_cast<std::size_t>(rank)]) {
                std::cerr << "duplicate onset for D_" << rank << "\n";
                std::exit(2);
            }
            seen[static_cast<std::size_t>(rank)] = true;
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
    int progress_step = 10;
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
        std::cerr << "invalid D bridge frontier GMP arguments\n";
        return 2;
    }

    const std::vector<int> onsets = read_d_onsets(onset_log, rank_lo, rank_hi);
    std::vector<int> bridge_max(static_cast<std::size_t>(rank_hi + 1), 0);
    int chain_hi = 0;
    int minimum_a_free_slack = std::numeric_limits<int>::max();
    for (int rank = rank_lo; rank <= rank_hi; ++rank) {
        const int a_free_slack =
            2 * rank - onsets[static_cast<std::size_t>(rank)] - 2;
        if (a_free_slack < 0) {
            std::cerr << "onset for D_" << rank
                      << " leaves the A-free moment window: onset="
                      << onsets[static_cast<std::size_t>(rank)]
                      << " maximum_moment="
                      << onsets[static_cast<std::size_t>(rank)] + 2
                      << " window_endpoint=" << 2 * rank << "\n";
            return 2;
        }
        minimum_a_free_slack = std::min(minimum_a_free_slack, a_free_slack);
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
    std::cout << "D interval stable-envelope bridge frontier GMP verifier\n"
              << "rank_lo=" << rank_lo << " rank_hi=" << rank_hi
              << " chain_lo=" << chain_lo << " chain_hi=" << chain_hi
              << " max_index=" << max_index
              << " minimum_a_free_slack=" << minimum_a_free_slack
              << " onset_log=" << onset_log
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;

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
            rank_lo,
            rank_hi,
            bridge_max,
            active_rows
        );
        if (frontier != 0) {
            results[static_cast<std::size_t>(offset)] = verify_chain_m_exact(
                chain_m,
                frontier,
                active_rows,
                stable
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
    for (const ChainResult& result : results) {
        if (result.frontier_rank == 0) continue;
        if (result.lower_sign <= 0) ++failures;
        if (result.lower_sign > 0) {
            if (!have_worst || result.relative_log10 < worst.relative_log10) {
                worst = result;
                have_worst = true;
            }
        }
    }

    std::cout << std::setprecision(18)
              << "worst_positive\tm\todd_n\tfrontier\tactive_rows"
              << "\tlower_log10\tstable_log10\trelative_log10\n";
    if (have_worst) {
        std::cout << "worst_positive\t" << worst.chain_m
                  << '\t' << worst.odd_n
                  << "\tD_" << worst.frontier_rank
                  << '\t' << worst.active_rows
                  << '\t' << worst.lower_log10
                  << '\t' << worst.stable_log10
                  << '\t' << worst.relative_log10
                  << '\n';
    }

    std::cout << "sample\tm\todd_n\tfrontier\tactive_rows"
              << "\tlower_sign\tlower_log10\tstable_log10\trelative_log10\n";
    const int sample_stride = std::max(1, total / 10);
    for (int offset = 0; offset < total; offset += sample_stride) {
        const ChainResult& result = results[static_cast<std::size_t>(offset)];
        if (result.frontier_rank == 0) continue;
        std::cout << "sample\t" << result.chain_m
                  << '\t' << result.odd_n
                  << "\tD_" << result.frontier_rank
                  << '\t' << result.active_rows
                  << '\t' << result.lower_sign
                  << '\t' << result.lower_log10
                  << '\t' << result.stable_log10
                  << '\t' << result.relative_log10
                  << '\n';
    }

    std::cout << "SUMMARY failures=" << failures;
    if (have_worst) {
        std::cout << " worst_m=" << worst.chain_m
                  << " worst_odd_n=" << worst.odd_n
                  << " worst_frontier=D_" << worst.frontier_rank
                  << " worst_relative_log10=" << worst.relative_log10
                  << " worst_lower_log10=" << worst.lower_log10;
    }
    std::cout << std::endl;
    return failures == 0 ? 0 : 1;
}
