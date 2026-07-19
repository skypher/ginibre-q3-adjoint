#define main classical_boundary_certificate_legacy_main
#include "classical_boundary_certificate.cpp"
#undef main

#include "full_q3_bcd_low_tail_data.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

struct FullQ3ResidualRow {
    char family;
    int rank;
    int tail_onset;
    int stable_through;
    int residual_through;
    std::vector<BigInt> moments;
};

struct FullQ3PieriMomentTask {
    int moment;
    std::vector<std::size_t> active_rows;
    std::vector<std::vector<Partition>> target_sets;
};

std::vector<std::vector<BigInt>> run_full_q3_reverse_tasks(
    const std::vector<FullQ3PieriMomentTask>& tasks,
    bool progress,
    const std::string& family_label
) {
    std::vector<std::vector<BigInt>> results(tasks.size());
    int completed = 0;
#ifdef _OPENMP
    const int thread_count = omp_get_max_threads();
#else
    const int thread_count = 1;
#endif
    std::cout << "FULL_Q3_RESIDUAL family=" << family_label
              << " moment_tasks=" << tasks.size()
              << " parallel_threads=" << thread_count << '\n';
    const int task_count = static_cast<int>(tasks.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int task_index = 0; task_index < task_count; ++task_index) {
        const FullQ3PieriMomentTask& task =
            tasks[static_cast<std::size_t>(task_index)];
        results[static_cast<std::size_t>(task_index)] =
            pieri_e2_coefficient_sums_reverse(
                task.moment,
                task.target_sets,
                false,
                "full-Q3 " + family_label + " m_" + std::to_string(task.moment)
            );
        if (progress) {
#ifdef _OPENMP
#pragma omp critical(full_q3_residual_progress)
#endif
            {
                ++completed;
                std::cout << "FULL_Q3_RESIDUAL family=" << family_label
                          << " completed_task=" << completed << '/' << task_count
                          << " moment=" << task.moment << std::endl;
            }
        }
    }
    return results;
}

std::vector<BigInt> regenerate_stable_classical_moments(int maximum) {
    if (maximum < 0) throw std::runtime_error("negative stable-moment maximum");
    std::vector<BigInt> moments(static_cast<std::size_t>(maximum + 1));
    moments[0] = 1;
    if (maximum >= 1) moments[1] = 0;
    for (int n = 1; n < maximum; ++n) {
        moments[static_cast<std::size_t>(n + 1)] =
            n * moments[static_cast<std::size_t>(n)]
            + n * moments[static_cast<std::size_t>(n - 1)];
        if (n >= 2) {
            moments[static_cast<std::size_t>(n + 1)] -=
                binom_int(n, 2) * moments[static_cast<std::size_t>(n - 2)];
        }
    }
    return moments;
}

BigInt full_q3_hierarchy_value(
    const std::vector<BigInt>& moments,
    int a,
    int n
) {
    if (a < 0 || n < 0) throw std::runtime_error("negative full-Q3 index");
    const int total_degree = 2 * a + n;
    if (total_degree >= static_cast<int>(moments.size())) {
        throw std::runtime_error("full-Q3 value exceeds moment range");
    }
    BigInt answer = 0;
    for (int j = 0; j <= 2 * a; ++j) {
        for (int k = 0; k <= n; ++k) {
            BigInt term =
                binom_int(2 * a, j) * binom_int(n, k)
                * moments[static_cast<std::size_t>(2 * a - j + k)]
                * moments[static_cast<std::size_t>(j + n - k)];
            if (j % 2 == 0) {
                answer += term;
            } else {
                answer -= term;
            }
        }
    }
    return answer;
}

int first_odd_strictly_above(int value) {
    int answer = value + 1;
    if (answer % 2 == 0) ++answer;
    return answer;
}

void run_full_q3_b18_d31_residual_certificate(
    const std::string& stable_path,
    bool progress
) {
    constexpr int maximum_moment = 71;
    std::vector<BigInt> stable = read_stable_moments(stable_path, maximum_moment);
    const std::vector<BigInt> regenerated =
        regenerate_stable_classical_moments(maximum_moment);
    if (stable != regenerated) {
        throw std::runtime_error(
            "vendored stable classical moments fail the exact recurrence"
        );
    }

    std::vector<FullQ3ResidualRow> b_rows;
    for (std::size_t offset = 0;
         offset < full_q3_bcd_low_tail::b_onsets.size();
         ++offset) {
        const int rank =
            full_q3_bcd_low_tail::b_first_rank + static_cast<int>(offset);
        const int onset = full_q3_bcd_low_tail::b_onsets[offset];
        b_rows.push_back(FullQ3ResidualRow{
            'B', rank, onset, 2 * rank + 1, onset - 2, stable
        });
    }

    std::vector<FullQ3ResidualRow> d_rows;
    for (std::size_t offset = 0;
         offset < full_q3_bcd_low_tail::d_onsets.size();
         ++offset) {
        const int rank =
            full_q3_bcd_low_tail::d_first_rank + static_cast<int>(offset);
        const int onset = full_q3_bcd_low_tail::d_onsets[offset];
        d_rows.push_back(FullQ3ResidualRow{
            'D', rank, onset, rank - 1, onset - 2, stable
        });
    }

    std::cout << "full-Q3 exact finite residual certificate"
              << " B_18..B_21 and D_31..D_70\n";
    std::cout << "stable_moments_regenerated=" << stable.size() << '\n';

    std::vector<FullQ3PieriMomentTask> b_tasks;
    for (int moment = 0; moment <= maximum_moment; ++moment) {
        FullQ3PieriMomentTask task{moment, {}, {}};
        for (std::size_t row_index = 0; row_index < b_rows.size(); ++row_index) {
            const FullQ3ResidualRow& row = b_rows[row_index];
            if (moment < 2 * row.rank + 2 || moment > row.residual_through) continue;
            task.active_rows.push_back(row_index);
            task.target_sets.push_back(b_boundary_targets(row.rank, moment));
            if (task.target_sets.back().empty()) {
                throw std::runtime_error("empty active B boundary target set");
            }
        }
        if (!task.active_rows.empty()) b_tasks.push_back(std::move(task));
    }
    const std::vector<std::vector<BigInt>> b_sums =
        run_full_q3_reverse_tasks(b_tasks, progress, "B");
    for (std::size_t task_index = 0; task_index < b_tasks.size(); ++task_index) {
        const FullQ3PieriMomentTask& task = b_tasks[task_index];
        const std::vector<BigInt>& sums = b_sums[task_index];
        if (sums.size() != task.active_rows.size()) {
            throw std::runtime_error("B reverse-Pieri color ledger mismatch");
        }
        for (std::size_t color = 0; color < task.active_rows.size(); ++color) {
            FullQ3ResidualRow& row = b_rows[task.active_rows[color]];
            if (sums[color] <= 0) {
                throw std::runtime_error("nonpositive B boundary coefficient sum");
            }
            row.moments[static_cast<std::size_t>(task.moment)] -= sums[color];
            if (row.moments[static_cast<std::size_t>(task.moment)] < 0) {
                throw std::runtime_error("negative B invariant-tensor moment");
            }
            std::cout << "B_" << row.rank << " m_" << task.moment
                      << " correction=-" << sums[color]
                      << " moment="
                      << row.moments[static_cast<std::size_t>(task.moment)]
                      << '\n';
        }
    }

    std::vector<FullQ3PieriMomentTask> d_tasks;
    for (int moment = 0; moment <= maximum_moment; ++moment) {
        FullQ3PieriMomentTask task{moment, {}, {}};
        for (std::size_t row_index = 0; row_index < d_rows.size(); ++row_index) {
            const FullQ3ResidualRow& row = d_rows[row_index];
            if (moment < row.rank || moment > row.residual_through) continue;
            if (moment > 2 * row.rank) {
                throw std::runtime_error(
                    "D residual crossed the O-even stability boundary"
                );
            }
            task.active_rows.push_back(row_index);
            task.target_sets.push_back(
                d_determinant_targets(row.rank, moment - row.rank)
            );
            if (task.target_sets.back().empty()) {
                throw std::runtime_error("empty active D determinant target set");
            }
        }
        if (!task.active_rows.empty()) d_tasks.push_back(std::move(task));
    }
    const std::vector<std::vector<BigInt>> d_sums =
        run_full_q3_reverse_tasks(d_tasks, progress, "D");
    for (std::size_t task_index = 0; task_index < d_tasks.size(); ++task_index) {
        const FullQ3PieriMomentTask& task = d_tasks[task_index];
        const std::vector<BigInt>& sums = d_sums[task_index];
        if (sums.size() != task.active_rows.size()) {
            throw std::runtime_error("D reverse-Pieri color ledger mismatch");
        }
        for (std::size_t color = 0; color < task.active_rows.size(); ++color) {
            FullQ3ResidualRow& row = d_rows[task.active_rows[color]];
            const int depth = task.moment - row.rank;
            if (depth <= 4) {
                const BigInt closed =
                    d_determinant_correction(row.rank, task.moment, false);
                if (sums[color] != closed) {
                    throw std::runtime_error(
                        "D closed-form/reverse-Pieri correction mismatch"
                    );
                }
            }
            if (sums[color] <= 0) {
                throw std::runtime_error("nonpositive D determinant correction");
            }
            row.moments[static_cast<std::size_t>(task.moment)] += sums[color];
            std::cout << "D_" << row.rank << " m_" << task.moment
                      << " determinant_correction=" << sums[color]
                      << " moment="
                      << row.moments[static_cast<std::size_t>(task.moment)]
                      << '\n';
        }
    }

    std::size_t b_checked = 0;
    std::size_t d_checked = 0;
    bool have_minimum = false;
    BigInt minimum;
    char minimum_family = 0;
    int minimum_rank = 0;
    int minimum_a = 0;
    int minimum_n = 0;

    const auto verify_rows = [&](const std::vector<FullQ3ResidualRow>& rows) {
        for (const FullQ3ResidualRow& row : rows) {
            const int first_total = first_odd_strictly_above(row.stable_through);
            std::size_t row_checked = 0;
            for (int total = first_total; total < row.tail_onset; total += 2) {
                if (total > row.residual_through) {
                    throw std::runtime_error("residual total-degree ledger overflow");
                }
                for (int a = 2; 2 * a < total; ++a) {
                    const int n = total - 2 * a;
                    if (n <= 0 || n % 2 == 0) {
                        throw std::runtime_error("invalid residual hierarchy index");
                    }
                    const BigInt value = full_q3_hierarchy_value(row.moments, a, n);
                    if (value <= 0) {
                        throw std::runtime_error(
                            std::string("nonpositive full-Q3 residual at ")
                            + row.family + "_" + std::to_string(row.rank)
                            + ", a=" + std::to_string(a)
                            + ", n=" + std::to_string(n)
                        );
                    }
                    if (!have_minimum || value < minimum) {
                        minimum = value;
                        minimum_family = row.family;
                        minimum_rank = row.rank;
                        minimum_a = a;
                        minimum_n = n;
                        have_minimum = true;
                    }
                    ++row_checked;
                }
            }
            if (row.family == 'B') {
                b_checked += row_checked;
            } else if (row.family == 'D') {
                d_checked += row_checked;
            } else {
                throw std::runtime_error("unknown full-Q3 residual family");
            }
            std::cout << "FULL_Q3_RESIDUAL row=" << row.family << '_' << row.rank
                      << " stable_through=" << row.stable_through
                      << " tail_onset=" << row.tail_onset
                      << " pairs=" << row_checked << '\n';
        }
    };

    verify_rows(b_rows);
    verify_rows(d_rows);
    if (b_checked != 137U || d_checked != 4732U) {
        throw std::runtime_error("full-Q3 residual scope count mismatch");
    }
    if (!have_minimum) {
        throw std::runtime_error("empty full-Q3 residual certificate");
    }
    std::cout << "FULL_Q3_RESIDUAL B_pairs=" << b_checked
              << " D_pairs=" << d_checked
              << " total_pairs=" << (b_checked + d_checked) << '\n';
    std::cout << "FULL_Q3_RESIDUAL minimum="
              << minimum_family << '_' << minimum_rank
              << " a=" << minimum_a << " n=" << minimum_n
              << " value=" << minimum << '\n';
    std::cout << "FULL_Q3_B18_D31_RESIDUAL VERIFICATION: ALL PASS\n";
}

int main(int argc, char** argv) {
    try {
        bool progress = false;
        std::string stable_path = "references/oeis_A002137_stable.txt";
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--progress") {
                progress = true;
            } else if (argument == "--stable-moments" && index + 1 < argc) {
                stable_path = argv[++index];
            } else {
                std::cerr << "usage: " << argv[0]
                          << " [--progress] [--stable-moments PATH]\n";
                return EXIT_FAILURE;
            }
        }
        run_full_q3_b18_d31_residual_certificate(stable_path, progress);
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "FULL_Q3_B18_D31_RESIDUAL VERIFICATION: FAIL: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

