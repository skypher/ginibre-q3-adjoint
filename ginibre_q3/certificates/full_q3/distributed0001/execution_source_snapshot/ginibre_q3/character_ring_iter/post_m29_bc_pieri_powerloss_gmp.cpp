#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double LOG10_2 = 0.301029995663981195213738894724493027;

struct Row {
    char family = 'B';
    int rank = 0;
};

struct Worst {
    Row row;
    int onset = 0;
    int boundary = 0;
    int q = 0;
    int moment = 0;
    double ratio_log10 = -std::numeric_limits<double>::infinity();
    double margin_log10 = std::numeric_limits<double>::infinity();
};

int row_index(const Row& row, int c_rank_hi) {
    if (row.family == 'B') return row.rank;
    return c_rank_hi + 1 + row.rank;
}

std::vector<Row> make_rows(int b_rank_lo, int b_rank_hi, int c_rank_lo, int c_rank_hi) {
    std::vector<Row> rows;
    for (int rank = b_rank_lo; rank <= b_rank_hi; ++rank) rows.push_back({'B', rank});
    for (int rank = c_rank_lo; rank <= c_rank_hi; ++rank) rows.push_back({'C', rank});
    return rows;
}

std::vector<int> read_onsets(
    const std::string& path,
    const std::vector<Row>& rows,
    int b_rank_hi,
    int c_rank_hi
) {
    std::vector<int> onsets(static_cast<std::size_t>(c_rank_hi + 1 + c_rank_hi + 1), 0);
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
        if (family == 'B' && 0 <= rank && rank <= b_rank_hi) {
            onsets[static_cast<std::size_t>(rank)] = onset;
        } else if (family == 'C' && 0 <= rank && rank <= c_rank_hi) {
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
            std::exit(1);
        }
        moments[static_cast<std::size_t>(index + 1)] = std::move(value);
        if (progress && ((index + 1) == max_index || (index + 1) % 1000 == 0)) {
            std::cout << "stable_progress index=" << (index + 1)
                      << "/" << max_index << std::endl;
            std::fflush(stdout);
        }
    }
    return moments;
}

mpz_class binom_initial(int n, int k) {
    if (k < 0 || k > n) return 0;
    mpz_class out;
    mpz_bin_uiui(out.get_mpz_t(), static_cast<unsigned long>(n), static_cast<unsigned long>(k));
    return out;
}

void advance_binom_same_k(mpz_class& binom, int n, int k) {
    binom *= (n + 1);
    binom /= (n + 1 - k);
}

void update_worst_ratio(Worst& worst, const Worst& candidate) {
    if (candidate.ratio_log10 > worst.ratio_log10) worst = candidate;
}

void update_worst_margin(Worst& worst, const Worst& candidate) {
    if (candidate.margin_log10 < worst.margin_log10) worst = candidate;
}

}  // namespace

int main(int argc, char** argv) {
    int b_rank_lo = 14;
    int b_rank_hi = 217;
    int c_rank_lo = 20;
    int c_rank_hi = 217;
    int loss_q_mult = 2;
    int loss_add = 0;
    bool progress = false;
    std::string onset_log;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--onset-log" && i + 1 < argc) {
            onset_log = argv[++i];
        } else if (arg == "--loss-q-mult" && i + 1 < argc) {
            loss_q_mult = std::atoi(argv[++i]);
        } else if (arg == "--loss-add" && i + 1 < argc) {
            loss_add = std::atoi(argv[++i]);
        } else if (arg == "--progress") {
            progress = true;
        } else {
            std::cerr
                << "usage: " << argv[0]
                << " --onset-log PATH [--loss-q-mult A] [--loss-add B] [--progress]\n";
            return 2;
        }
    }

    if (onset_log.empty() || loss_q_mult < 0 || loss_add < 0) {
        std::cerr << "invalid B/C power-loss GMP verifier arguments\n";
        return 2;
    }

    const std::vector<Row> rows = make_rows(b_rank_lo, b_rank_hi, c_rank_lo, c_rank_hi);
    const std::vector<int> onsets = read_onsets(onset_log, rows, b_rank_hi, c_rank_hi);
    int max_index = 0;
    std::uint64_t total_checks = 0;
    for (const Row& row : rows) {
        const int onset = onsets[static_cast<std::size_t>(row_index(row, c_rank_hi))];
        max_index = std::max(max_index, onset + 2);
        const int boundary = 2 * row.rank + 2;
        total_checks += static_cast<std::uint64_t>(onset + 2 - boundary + 1);
    }

    std::cout << "B/C residual Pieri-compression power-loss arithmetic GMP verifier\n"
              << "VERIFIES exact inequality: 2^(A*q+B+1)*binom(j,q)*s_{j-q} <= s_j\n"
              << "A=" << loss_q_mult << " B=" << loss_add << " q=rank+1\n"
              << "This is only the arithmetic implication of a possible branch-loss supplier.\n"
              << "rows=" << rows.size()
              << " total_checks=" << total_checks
              << " max_index=" << max_index
              << " onset_log=" << onset_log
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;

    const std::vector<mpz_class> stable = stable_moments(max_index, progress);
    std::vector<Worst> worst_ratios(rows.size());
    std::vector<Worst> worst_margins(rows.size());
    std::vector<int> row_failures(rows.size(), 0);
    int completed = 0;

#pragma omp parallel for schedule(dynamic, 1)
    for (int row_slot = 0; row_slot < static_cast<int>(rows.size()); ++row_slot) {
        const Row row = rows[static_cast<std::size_t>(row_slot)];
        const int onset = onsets[static_cast<std::size_t>(row_index(row, c_rank_hi))];
        const int boundary = 2 * row.rank + 2;
        const int q = row.rank + 1;
        const int exponent = loss_q_mult * q + loss_add + 1;
        mpz_class pow2_factor;
        mpz_ui_pow_ui(pow2_factor.get_mpz_t(), 2UL, static_cast<unsigned long>(exponent));
        mpz_class binom = binom_initial(boundary, q);
        Worst row_worst_ratio;
        Worst row_worst_margin;
        row_worst_ratio.row = row;
        row_worst_margin.row = row;
        row_worst_ratio.onset = onset;
        row_worst_margin.onset = onset;
        row_worst_ratio.boundary = boundary;
        row_worst_margin.boundary = boundary;
        row_worst_ratio.q = q;
        row_worst_margin.q = q;
        int failures = 0;

        for (int moment = boundary; moment <= onset + 2; ++moment) {
            const mpz_class bound =
                pow2_factor * binom * stable[static_cast<std::size_t>(moment - q)];
            const mpz_class margin = stable[static_cast<std::size_t>(moment)] - bound;
            const double ratio_log10 =
                log10_abs_mpz(bound)
                - log10_abs_mpz(stable[static_cast<std::size_t>(moment)]);
            const double margin_log10 = log10_abs_mpz(margin);
            Worst candidate{row, onset, boundary, q, moment, ratio_log10, margin_log10};
            update_worst_ratio(row_worst_ratio, candidate);
            update_worst_margin(row_worst_margin, candidate);
            if (sgn(margin) <= 0) {
                ++failures;
#pragma omp critical
                {
                    std::cout << "FAIL row=" << row.family << "_" << row.rank
                              << " onset=" << onset
                              << " boundary=" << boundary
                              << " q=" << q
                              << " moment=" << moment
                              << " ratio_log10=" << display(ratio_log10)
                              << " margin_sign=" << sgn(margin)
                              << " margin_log10=" << display(margin_log10)
                              << std::endl;
                    std::fflush(stdout);
                }
            }
            if (moment < onset + 2) advance_binom_same_k(binom, moment, q);
        }

        worst_ratios[static_cast<std::size_t>(row_slot)] = row_worst_ratio;
        worst_margins[static_cast<std::size_t>(row_slot)] = row_worst_margin;
        row_failures[static_cast<std::size_t>(row_slot)] = failures;
        int done_now = 0;
#pragma omp atomic capture
        done_now = ++completed;
        if (progress && (done_now == static_cast<int>(rows.size()) || done_now % 25 == 0)) {
#pragma omp critical
            {
                std::cout << "progress completed=" << done_now
                          << "/" << rows.size()
                          << " last=" << row.family << "_" << row.rank
                          << " failures=" << failures
                          << std::endl;
                std::fflush(stdout);
            }
        }
    }

    Worst worst_ratio;
    Worst worst_margin;
    int failures = 0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        update_worst_ratio(worst_ratio, worst_ratios[i]);
        update_worst_margin(worst_margin, worst_margins[i]);
        failures += row_failures[i];
    }

    std::cout << "worst_ratio"
              << "\trow=" << worst_ratio.row.family << "_" << worst_ratio.row.rank
              << "\tonset=" << worst_ratio.onset
              << "\tboundary=" << worst_ratio.boundary
              << "\tq=" << worst_ratio.q
              << "\tmoment=" << worst_ratio.moment
              << "\tlog10_bound_over_stable=" << display(worst_ratio.ratio_log10)
              << "\n";
    std::cout << "worst_margin"
              << "\trow=" << worst_margin.row.family << "_" << worst_margin.row.rank
              << "\tonset=" << worst_margin.onset
              << "\tboundary=" << worst_margin.boundary
              << "\tq=" << worst_margin.q
              << "\tmoment=" << worst_margin.moment
              << "\tmargin_log10=" << display(worst_margin.margin_log10)
              << "\n";
    std::cout << "SUMMARY failures=" << failures
              << " rows=" << rows.size()
              << " total_checks=" << total_checks
              << " max_index=" << max_index
              << " loss_q_mult=" << loss_q_mult
              << " loss_add=" << loss_add
              << " worst_ratio_log10=" << display(worst_ratio.ratio_log10)
              << " worst_margin_log10=" << display(worst_margin.margin_log10)
              << std::endl;
    std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << std::endl;
    return failures == 0 ? 0 : 1;
}
