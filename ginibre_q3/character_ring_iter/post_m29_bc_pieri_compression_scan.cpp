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

constexpr long double NEG_INF = -std::numeric_limits<long double>::infinity();
constexpr long double LOG10_2 = 0.301029995663981195213738894724493027L;
constexpr long double LOG10_HALF = -LOG10_2;

struct Row {
    char family = 'B';
    int rank = 0;
};

struct Worst {
    long double log10_bound = NEG_INF;
    Row row;
    int onset = 0;
    int boundary = 0;
    int q = 0;
    int moment = 0;
};

struct SignedLog {
    int sign = 0;
    long double log_abs = NEG_INF;
};

void add_signed_log(SignedLog& acc, int sign, long double log_abs) {
    if (sign == 0 || !std::isfinite(log_abs)) return;
    const int term_sign = sign > 0 ? 1 : -1;
    if (acc.sign == 0) {
        acc.sign = term_sign;
        acc.log_abs = log_abs;
        return;
    }
    const long double hi = std::max(acc.log_abs, log_abs);
    const long double lo = std::min(acc.log_abs, log_abs);
    if (acc.sign == term_sign) {
        acc.log_abs = hi - lo > 80.0L
            ? hi
            : hi + std::log10(1.0L + std::pow(10.0L, lo - hi));
        return;
    }
    if (std::fabs(acc.log_abs - log_abs) < 1.0e-18L) {
        acc.sign = 0;
        acc.log_abs = NEG_INF;
        return;
    }
    const int hi_sign = acc.log_abs > log_abs ? acc.sign : term_sign;
    if (hi - lo > 80.0L) {
        acc.sign = hi_sign;
        acc.log_abs = hi;
        return;
    }
    const long double diff = 1.0L - std::pow(10.0L, lo - hi);
    if (!(diff > 0.0L)) {
        acc.sign = 0;
        acc.log_abs = NEG_INF;
        return;
    }
    acc.sign = hi_sign;
    acc.log_abs = hi + std::log10(diff);
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

long double log10_binom(const std::vector<long double>& log_fact, int n, int k) {
    if (k < 0 || k > n) return NEG_INF;
    return log_fact[static_cast<std::size_t>(n)]
        - log_fact[static_cast<std::size_t>(k)]
        - log_fact[static_cast<std::size_t>(n - k)];
}

std::vector<long double> stable_moment_logs(int max_index) {
    std::vector<SignedLog> moments(static_cast<std::size_t>(max_index + 1));
    moments[0] = SignedLog{1, 0.0L};
    if (max_index >= 1) moments[1] = SignedLog{0, NEG_INF};
    for (int index = 1; index < max_index; ++index) {
        SignedLog next;
        if (moments[static_cast<std::size_t>(index)].sign != 0) {
            add_signed_log(
                next,
                moments[static_cast<std::size_t>(index)].sign,
                std::log10(static_cast<long double>(index))
                    + moments[static_cast<std::size_t>(index)].log_abs
            );
        }
        if (moments[static_cast<std::size_t>(index - 1)].sign != 0) {
            add_signed_log(
                next,
                moments[static_cast<std::size_t>(index - 1)].sign,
                std::log10(static_cast<long double>(index))
                    + moments[static_cast<std::size_t>(index - 1)].log_abs
            );
        }
        if (index >= 2 && moments[static_cast<std::size_t>(index - 2)].sign != 0) {
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
    std::vector<long double> logs(static_cast<std::size_t>(max_index + 1), NEG_INF);
    for (int index = 0; index <= max_index; ++index) {
        logs[static_cast<std::size_t>(index)] =
            moments[static_cast<std::size_t>(index)].log_abs;
    }
    return logs;
}

std::vector<Row> make_rows(int b_rank_lo, int b_rank_hi, int c_rank_lo, int c_rank_hi) {
    std::vector<Row> rows;
    for (int rank = b_rank_lo; rank <= b_rank_hi; ++rank) rows.push_back({'B', rank});
    for (int rank = c_rank_lo; rank <= c_rank_hi; ++rank) rows.push_back({'C', rank});
    return rows;
}

int row_index(const Row& row, int c_rank_hi) {
    if (row.family == 'B') return row.rank;
    return c_rank_hi + 1 + row.rank;
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

void update_worst(Worst& worst, const Worst& candidate) {
    if (candidate.log10_bound > worst.log10_bound) worst = candidate;
}

std::string display(long double value) {
    if (std::isinf(value)) return value > 0 ? "+inf" : "-inf";
    std::ostringstream out;
    out << std::setprecision(18) << static_cast<double>(value);
    return out.str();
}

}  // namespace

int main(int argc, char** argv) {
    int b_rank_lo = 14;
    int b_rank_hi = 217;
    int c_rank_lo = 20;
    int c_rank_hi = 217;
    bool progress = false;
    std::string onset_log;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--onset-log" && i + 1 < argc) {
            onset_log = argv[++i];
        } else if (arg == "--progress") {
            progress = true;
        } else if (arg == "--b-rank-lo" && i + 1 < argc) {
            b_rank_lo = std::atoi(argv[++i]);
        } else if (arg == "--b-rank-hi" && i + 1 < argc) {
            b_rank_hi = std::atoi(argv[++i]);
        } else if (arg == "--c-rank-lo" && i + 1 < argc) {
            c_rank_lo = std::atoi(argv[++i]);
        } else if (arg == "--c-rank-hi" && i + 1 < argc) {
            c_rank_hi = std::atoi(argv[++i]);
        } else {
            std::cerr
                << "usage: " << argv[0]
                << " --onset-log PATH [--progress]"
                << " [--b-rank-lo N] [--b-rank-hi N]"
                << " [--c-rank-lo N] [--c-rank-hi N]\n";
            return 2;
        }
    }

    if (onset_log.empty()
        || b_rank_lo < 1
        || c_rank_lo < 1
        || b_rank_hi < b_rank_lo
        || c_rank_hi < c_rank_lo) {
        std::cerr << "invalid B/C Pieri compression scan arguments\n";
        return 2;
    }

    const std::vector<Row> rows = make_rows(b_rank_lo, b_rank_hi, c_rank_lo, c_rank_hi);
    const std::vector<int> onsets = read_onsets(onset_log, rows, b_rank_hi, c_rank_hi);
    int max_index = 0;
    for (const Row& row : rows) {
        max_index = std::max(
            max_index,
            onsets[static_cast<std::size_t>(row_index(row, c_rank_hi))] + 2
        );
    }

    std::cout << "B/C Pieri compression-bound diagnostic\n"
              << "DIAGNOSTIC_ONLY candidate bound is "
              << "|delta_j| <= binom(j,q) stable_{j-q}; "
              << "q_half=rank+1, q_full=2*rank+2\n"
              << "b_rank_lo=" << b_rank_lo << " b_rank_hi=" << b_rank_hi
              << " c_rank_lo=" << c_rank_lo << " c_rank_hi=" << c_rank_hi
              << " rows=" << rows.size()
              << " max_index=" << max_index
              << " onset_log=" << onset_log
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;

    const std::vector<long double> log_fact = log10_factorials(max_index);
    const std::vector<long double> stable_logs = stable_moment_logs(max_index);

    std::vector<Worst> half_worsts(rows.size());
    std::vector<Worst> full_worsts(rows.size());
    int completed = 0;
#pragma omp parallel for schedule(dynamic, 1)
    for (int row_index_value = 0; row_index_value < static_cast<int>(rows.size()); ++row_index_value) {
        const Row row = rows[static_cast<std::size_t>(row_index_value)];
        const int onset = onsets[static_cast<std::size_t>(row_index(row, c_rank_hi))];
        const int boundary = 2 * row.rank + 2;
        const int q_half = row.rank + 1;
        const int q_full = boundary;
        Worst half;
        Worst full;
        half.row = row;
        half.onset = onset;
        half.boundary = boundary;
        half.q = q_half;
        full.row = row;
        full.onset = onset;
        full.boundary = boundary;
        full.q = q_full;

        for (int moment = boundary; moment <= onset + 2; ++moment) {
            if (moment >= q_half) {
                const long double value =
                    log10_binom(log_fact, moment, q_half)
                    + stable_logs[static_cast<std::size_t>(moment - q_half)]
                    - stable_logs[static_cast<std::size_t>(moment)];
                update_worst(half, Worst{value, row, onset, boundary, q_half, moment});
            }
            if (moment >= q_full) {
                const long double value =
                    log10_binom(log_fact, moment, q_full)
                    + stable_logs[static_cast<std::size_t>(moment - q_full)]
                    - stable_logs[static_cast<std::size_t>(moment)];
                update_worst(full, Worst{value, row, onset, boundary, q_full, moment});
            }
        }

        half_worsts[static_cast<std::size_t>(row_index_value)] = half;
        full_worsts[static_cast<std::size_t>(row_index_value)] = full;
        int done_now = 0;
#pragma omp atomic capture
        done_now = ++completed;
        if (progress && (done_now == static_cast<int>(rows.size()) || done_now % 50 == 0)) {
#pragma omp critical
            {
                std::cout << "progress completed=" << done_now
                          << "/" << rows.size()
                          << " last=" << row.family << "_" << row.rank
                          << std::endl;
                std::fflush(stdout);
            }
        }
    }

    Worst worst_half;
    Worst worst_full;
    for (const Worst& item : half_worsts) update_worst(worst_half, item);
    for (const Worst& item : full_worsts) update_worst(worst_full, item);

    auto print_worst = [](const std::string& label, const Worst& worst) {
        std::cout << label
                  << "\trow=" << worst.row.family << "_" << worst.row.rank
                  << "\tonset=" << worst.onset
                  << "\tboundary=" << worst.boundary
                  << "\tq=" << worst.q
                  << "\tmoment=" << worst.moment
                  << "\tlog10_bound=" << display(worst.log10_bound)
                  << "\n";
    };

    print_worst("worst_half", worst_half);
    print_worst("worst_full", worst_full);
    std::cout << "SUMMARY worst_half_log10=" << display(worst_half.log10_bound)
              << " worst_half_below_half=" << (worst_half.log10_bound < LOG10_HALF)
              << " worst_full_log10=" << display(worst_full.log10_bound)
              << " worst_full_below_half=" << (worst_full.log10_bound < LOG10_HALF)
              << std::endl;
    return 0;
}
