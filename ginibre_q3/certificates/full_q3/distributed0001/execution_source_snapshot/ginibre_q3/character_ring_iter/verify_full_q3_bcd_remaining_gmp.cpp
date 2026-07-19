#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wpedantic"
#endif
#define main classical_boundary_certificate_legacy_main
#include "classical_boundary_certificate.cpp"
#undef main
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "full_q3_bcd_remaining_data.hpp"

#include <array>
#include <initializer_list>
#include <limits>

namespace {

enum class MomentOperation { set, add, subtract };

struct RemainingRow {
    char family;
    int rank;
    int stable_through;
    int tail_onset;
    int moment_through;
    full_q3_bcd_remaining::TailMethod tail_method;
    std::vector<BigInt> moments;
};

struct MomentAssignment {
    std::size_t row;
    MomentOperation operation;
};

struct MomentTask {
    int moment;
    std::string label;
    std::vector<MomentAssignment> assignments;
    std::vector<std::vector<Partition>> targets;
};

std::vector<Partition> b_kept_targets(int rank, int moment) {
    const int matrix_size = 2 * rank + 1;
    std::vector<Partition> targets;
    for (const Partition& half : integer_partitions_cached(moment)) {
        if (static_cast<int>(half.size()) > matrix_size) continue;
        Partition target;
        target.reserve(half.size());
        for (int part : half) target.push_back(2 * part);
        targets.push_back(std::move(target));
    }
    return targets;
}

std::vector<Partition> c_kept_targets(int rank, int moment) {
    std::vector<Partition> targets;
    for (const Partition& paired : integer_partitions_cached(moment)) {
        if (static_cast<int>(paired.size()) > rank) continue;
        Partition target;
        target.reserve(2 * paired.size());
        for (int part : paired) {
            target.push_back(part);
            target.push_back(part);
        }
        targets.push_back(conjugate_partition(target));
    }
    return targets;
}

std::vector<BigInt> regenerate_stable_moments(int maximum) {
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

BigInt hierarchy_value(
    const std::vector<BigInt>& moments,
    int a,
    int n
) {
    const int total = 2 * a + n;
    if (a < 0 || n < 0 || total >= static_cast<int>(moments.size())) {
        throw std::runtime_error("invalid full-Q3 hierarchy index");
    }
    BigInt answer = 0;
    for (int j = 0; j <= 2 * a; ++j) {
        for (int k = 0; k <= n; ++k) {
            const BigInt term =
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

int first_odd_above(int value) {
    int answer = value + 1;
    if (answer % 2 == 0) ++answer;
    return answer;
}

mpq_class rational(unsigned long numerator, unsigned long denominator) {
    if (denominator == 0UL) {
        throw std::runtime_error("zero exact-rational denominator");
    }
    mpq_class answer(numerator);
    answer /= denominator;
    answer.canonicalize();
    return answer;
}

mpq_class rational_power(mpq_class base, int exponent) {
    if (exponent < 0) {
        throw std::runtime_error("negative exact-rational exponent");
    }
    mpq_class answer(1);
    while (exponent > 0) {
        if (exponent % 2 != 0) answer *= base;
        exponent /= 2;
        if (exponent > 0) base *= base;
    }
    return answer;
}

BigInt factorial_exact(int value) {
    if (value < 0) throw std::runtime_error("negative factorial argument");
    BigInt answer = 1;
    for (int factor = 2; factor <= value; ++factor) answer *= factor;
    return answer;
}

BigInt power_of_two_exact(int exponent) {
    if (exponent < 0) throw std::runtime_error("negative power-of-two exponent");
    BigInt answer;
    mpz_ui_pow_ui(
        answer.get_mpz_t(),
        2UL,
        static_cast<unsigned long>(exponent)
    );
    return answer;
}

BigInt degree_factor_exact(std::initializer_list<int> degrees) {
    BigInt answer = 1;
    for (int degree : degrees) answer *= factorial_exact(degree);
    return answer;
}

int negative_endpoint(char family, int rank) {
    if (family == 'B' || family == 'C') return rank;
    if (family == 'D') return rank % 2 == 0 ? rank : rank - 2;
    throw std::runtime_error("unknown classical family in negative endpoint");
}

mpq_class exact_tail_margin(
    const mpq_class& probability,
    const mpq_class& push,
    int endpoint,
    int onset
) {
    if (probability <= 0 || probability > 1) {
        throw std::runtime_error("invalid exact tail probability");
    }
    if (push <= 2 * endpoint) {
        throw std::runtime_error("exact tail push does not exceed negative cap");
    }
    return probability * rational_power(push, onset)
        - rational_power(mpq_class(2 * endpoint), onset);
}

struct PolynomialTailSpec {
    char family;
    int rank;
    int power;
    unsigned long threshold_numerator;
    unsigned long threshold_denominator;
    unsigned long central_numerator;
    unsigned long central_denominator;
    int onset;
};

constexpr std::array<PolynomialTailSpec, 6> polynomial_tail_specs{{
    {'B', 3, 8, 227, 20, 23, 20, 37},
    {'B', 4, 10, 331, 20, 5, 4, 41},
    {'B', 5, 10, 75, 4, 6, 5, 55},
    {'C', 3, 8, 57, 5, 6, 5, 37},
    {'D', 4, 9, 57, 4, 23, 20, 47},
    {'D', 5, 7, 131, 10, 6, 5, 33},
}};

const PolynomialTailSpec* find_polynomial_tail(char family, int rank) {
    for (const PolynomialTailSpec& spec : polynomial_tail_specs) {
        if (spec.family == family && spec.rank == rank) return &spec;
    }
    return nullptr;
}

void verify_polynomial_tail(const RemainingRow& row) {
    const PolynomialTailSpec* spec = find_polynomial_tail(row.family, row.rank);
    if (spec == nullptr || spec->onset != row.tail_onset) {
        throw std::runtime_error("polynomial-tail ledger mismatch");
    }
    const int k = spec->power;
    const int required_moment = 4 * k + 2;
    if (required_moment > row.moment_through
        || required_moment >= static_cast<int>(row.moments.size())) {
        throw std::runtime_error("polynomial tail lacks required exact moments");
    }
    const mpq_class threshold = rational(
        spec->threshold_numerator, spec->threshold_denominator
    );
    const mpq_class central = rational(
        spec->central_numerator, spec->central_denominator
    );
    const mpq_class a_value =
        mpq_class(row.moments[static_cast<std::size_t>(2 * k + 1)])
        - threshold
            * mpq_class(row.moments[static_cast<std::size_t>(2 * k)]);
    const mpq_class b_value =
        mpq_class(row.moments[static_cast<std::size_t>(4 * k + 2)])
        - 2 * threshold
            * mpq_class(row.moments[static_cast<std::size_t>(4 * k + 1)])
        + threshold * threshold
            * mpq_class(row.moments[static_cast<std::size_t>(4 * k)]);
    if (a_value <= 0 || b_value <= 0 || central <= 1) {
        throw std::runtime_error("invalid polynomial one-sided-tail inputs");
    }
    const mpq_class probability =
        (1 - 1 / (central * central)) * a_value * a_value / b_value;
    const mpq_class push = threshold - central;
    const mpq_class margin = exact_tail_margin(
        probability,
        push,
        negative_endpoint(row.family, row.rank),
        row.tail_onset
    );
    if (margin <= 0) {
        throw std::runtime_error("polynomial-tail comparison is not positive");
    }
    std::cout << "FULL_Q3_EXACT_TAIL row=" << row.family << '_' << row.rank
              << " method=polynomial onset=" << row.tail_onset
              << " power=" << k
              << " margin_numerator=" << margin.get_num()
              << " margin_denominator=" << margin.get_den() << '\n';
}

void verify_rational_cap_tail(const RemainingRow& row) {
    const char family = row.family;
    const int rank = row.rank;
    int dimension = 0;
    int kappa = 0;
    int endpoint = 0;
    mpq_class radius;
    mpq_class central;
    mpq_class push;
    BigInt degree_factor;
    BigInt weyl_order;
    mpq_class cap_probability;

    if ((family == 'B' || family == 'C') && rank == 2) {
        dimension = 10;
        endpoint = 2;
        radius = rational(1, 1);
        central = rational(6, 5);
        push = rational(39, 5);
        const mpq_class width = rational(1, 4);
        const mpq_class cap_integral =
            rational_power(width, 10)
            * (rational(1, 21) - rational(2, 25) + rational(1, 21));
        cap_probability =
            rational(1, 512) * rational_power(width, 4) * cap_integral;
    } else if (family == 'B' && rank == 6) {
        dimension = 78;
        kappa = 11;
        endpoint = 6;
        radius = rational(1053, 50);
        central = rational(3, 2);
        push = rational(1386, 25);
        degree_factor = degree_factor_exact({2, 4, 6, 8, 10, 12});
        weyl_order = power_of_two_exact(6) * factorial_exact(6);
        const int shape = dimension / 2;
        const mpq_class density = 1 / (
            rational_power(rational(44, 7), rank / 2)
            * factorial_exact(shape)
        );
        cap_probability =
            rational(1, 1) / weyl_order
            * rational_power(1 - radius / 24, 2)
            * degree_factor * density
            * rational_power(radius / (2 * kappa), shape)
            / power_of_two_exact(rank);
        if (radius >= 24) {
            throw std::runtime_error("B6 product-sine cap left its domain");
        }
    } else if (family == 'D' && rank == 6) {
        dimension = 66;
        kappa = 10;
        endpoint = 6;
        radius = rational(371, 20);
        central = rational(7, 5);
        push = rational(921, 20);
        degree_factor = degree_factor_exact({2, 4, 6, 8, 10, 6});
        weyl_order = power_of_two_exact(5) * factorial_exact(6);
        const int shape = dimension / 2;
        const mpq_class density = 1 / (
            rational_power(rational(44, 7), rank / 2)
            * factorial_exact(shape)
        );
        cap_probability =
            rational(1, 1) / weyl_order
            * rational_power(1 - radius / 24, 2)
            * degree_factor * density
            * rational_power(radius / (2 * kappa), shape);
        if (radius >= 24) {
            throw std::runtime_error("D6 product-sine cap left its domain");
        }
    } else if (family == 'D' && rank == 7) {
        dimension = 91;
        kappa = 12;
        endpoint = 5;
        radius = rational(18, 1);
        central = rational(8, 5);
        push = rational(357, 5);
        degree_factor = degree_factor_exact({2, 4, 6, 8, 10, 12, 7});
        weyl_order = power_of_two_exact(6) * factorial_exact(7);
        const int half_index = (dimension + 1) / 2;
        const mpq_class density =
            rational_power(rational(4, 1), half_index)
            * factorial_exact(half_index)
            / (
                12 * rational_power(rational(22, 7), 4)
                * factorial_exact(2 * half_index)
            );
        cap_probability =
            rational(1, 1) / weyl_order
            * rational_power(1 - radius / 24, 2)
            * degree_factor * density
            * rational_power(radius / (2 * kappa), (dimension - 1) / 2)
            * rational(433, 500);
        if (radius >= 24
            || rational_power(rational(433, 500), 2)
                >= radius / (2 * kappa)
            || rational_power(rational(3, 2), 2) <= 2) {
            throw std::runtime_error("D7 rational cap bound is invalid");
        }
    } else if (family == 'D' && rank == 9) {
        dimension = 153;
        kappa = 16;
        endpoint = 7;
        radius = rational(763, 10);
        central = rational(8, 5);
        push = rational(751, 10);
        degree_factor = degree_factor_exact({2, 4, 6, 8, 10, 12, 14, 16, 9});
        weyl_order = power_of_two_exact(8) * factorial_exact(9);
        const int half_index = (dimension + 1) / 2;
        const mpq_class density =
            rational(7, 5)
            * rational_power(rational(4, 1), half_index)
            * factorial_exact(half_index)
            / (
                rational_power(rational(44, 7), 5)
                * factorial_exact(2 * half_index)
            );
        cap_probability =
            rational(1, 1) / weyl_order
            * rational_power(1 - 2 * radius / (24 * kappa), kappa)
            * degree_factor * density
            * rational_power(radius / (2 * kappa), (dimension - 1) / 2)
            * rational(193, 125);
        if (radius >= 12 * kappa
            || rational_power(rational(193, 125), 2)
                >= radius / (2 * kappa)
            || rational_power(rational(10, 7), 2) <= 2) {
            throw std::runtime_error("D9 rational cap bound is invalid");
        }
    } else {
        throw std::runtime_error("unknown exact rational-cap row");
    }

    if (push != dimension - radius - central
        || endpoint != negative_endpoint(family, rank)) {
        throw std::runtime_error("rational-cap geometry ledger mismatch");
    }
    if (kappa > 0 && radius >= 9 * kappa) {
        throw std::runtime_error("rational cap leaves the eigenangle chart");
    }
    const mpq_class central_probability = 1 - 1 / (central * central);
    const mpq_class probability = cap_probability * central_probability;
    const mpq_class margin = exact_tail_margin(
        probability, push, endpoint, row.tail_onset
    );
    if (margin <= 0) {
        throw std::runtime_error("rational-cap tail comparison is not positive");
    }
    std::cout << "FULL_Q3_EXACT_TAIL row=" << family << '_' << rank
              << " method=rational_cap onset=" << row.tail_onset
              << " margin_numerator=" << margin.get_num()
              << " margin_denominator=" << margin.get_den() << '\n';
}

void add_bc_task(
    std::vector<MomentTask>& tasks,
    std::vector<RemainingRow>& rows,
    char family,
    int moment
) {
    MomentTask task{moment, std::string(1, family), {}, {}};
    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        const RemainingRow& row = rows[row_index];
        if (row.family != family || moment <= row.stable_through
            || moment > row.moment_through) {
            continue;
        }
        std::vector<Partition> kept = family == 'B'
            ? b_kept_targets(row.rank, moment)
            : c_kept_targets(row.rank, moment);
        std::vector<Partition> overflow = family == 'B'
            ? b_boundary_targets(row.rank, moment)
            : c_boundary_targets(row.rank, moment);
        if (kept.empty() && overflow.empty()) {
            throw std::runtime_error("empty B/C partition ledger");
        }
        const bool use_kept = kept.size() <= overflow.size();
        task.assignments.push_back(MomentAssignment{
            row_index,
            use_kept ? MomentOperation::set : MomentOperation::subtract
        });
        task.targets.push_back(use_kept ? std::move(kept) : std::move(overflow));
    }
    if (!task.assignments.empty()) tasks.push_back(std::move(task));
}

void add_d_even_task(
    std::vector<MomentTask>& tasks,
    std::vector<RemainingRow>& rows,
    int moment
) {
    MomentTask task{moment, "D-even", {}, {}};
    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        const RemainingRow& row = rows[row_index];
        if (row.family != 'D' || moment <= 2 * row.rank
            || moment > row.moment_through) {
            continue;
        }
        std::vector<Partition> kept = o_even_kept_targets(row.rank, moment);
        std::vector<Partition> overflow = o_even_boundary_targets(row.rank, moment);
        if (kept.empty() && overflow.empty()) {
            throw std::runtime_error("empty D-even partition ledger");
        }
        const bool use_kept = kept.size() <= overflow.size();
        task.assignments.push_back(MomentAssignment{
            row_index,
            use_kept ? MomentOperation::set : MomentOperation::subtract
        });
        task.targets.push_back(use_kept ? std::move(kept) : std::move(overflow));
    }
    if (!task.assignments.empty()) tasks.push_back(std::move(task));
}

void add_d_determinant_task(
    std::vector<MomentTask>& tasks,
    const std::vector<RemainingRow>& rows,
    int moment
) {
    MomentTask task{moment, "D-det", {}, {}};
    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        const RemainingRow& row = rows[row_index];
        if (row.family != 'D' || moment < row.rank
            || moment > row.moment_through) {
            continue;
        }
        std::vector<Partition> targets =
            d_determinant_targets(row.rank, moment - row.rank);
        if (targets.empty()) {
            throw std::runtime_error("empty D determinant partition ledger");
        }
        task.assignments.push_back(
            MomentAssignment{row_index, MomentOperation::add}
        );
        task.targets.push_back(std::move(targets));
    }
    if (!task.assignments.empty()) tasks.push_back(std::move(task));
}

void run_remaining_certificate(
    const std::string& stable_path,
    int maximum,
    bool progress
) {
    if (maximum < 7
        || maximum > full_q3_bcd_remaining::required_maximum_moment) {
        throw std::runtime_error(
            "remaining-moment maximum exceeds the certified row ledger"
        );
    }
    const bool full_certificate =
        maximum == full_q3_bcd_remaining::required_maximum_moment;
    const std::vector<BigInt> stable = read_stable_moments(stable_path, maximum);
    if (stable != regenerate_stable_moments(maximum)) {
        throw std::runtime_error("stable-moment recurrence mismatch");
    }

    std::vector<RemainingRow> rows;
    std::size_t polynomial_count = 0;
    std::size_t rational_cap_count = 0;
    std::size_t directed_interval_count = 0;
    int largest_moment = 0;
    int expected_b_rank = 2;
    int expected_c_rank = 2;
    int expected_d_rank = 4;
    for (const full_q3_bcd_remaining::RowCutoff& cutoff :
         full_q3_bcd_remaining::row_cutoffs) {
        if ((cutoff.family == 'B' && cutoff.rank != expected_b_rank++)
            || (cutoff.family == 'C' && cutoff.rank != expected_c_rank++)
            || (cutoff.family == 'D' && cutoff.rank != expected_d_rank++)
            || (cutoff.family != 'B' && cutoff.family != 'C'
                && cutoff.family != 'D')) {
            throw std::runtime_error("remaining-row rank sequence mismatch");
        }
        const int stable_through = cutoff.family == 'D'
            ? cutoff.rank - 1
            : 2 * cutoff.rank + 1;
        if (cutoff.tail_onset % 2 == 0
            || cutoff.tail_onset <= stable_through
            || cutoff.moment_through < cutoff.tail_onset - 2) {
            throw std::runtime_error("invalid remaining-row coverage ledger");
        }
        largest_moment = std::max(largest_moment, cutoff.moment_through);
        switch (cutoff.tail_method) {
            case full_q3_bcd_remaining::TailMethod::polynomial:
                ++polynomial_count;
                break;
            case full_q3_bcd_remaining::TailMethod::rational_cap:
                ++rational_cap_count;
                break;
            case full_q3_bcd_remaining::TailMethod::directed_interval:
                ++directed_interval_count;
                break;
        }
        rows.push_back(RemainingRow{
            cutoff.family,
            cutoff.rank,
            stable_through,
            cutoff.tail_onset,
            cutoff.moment_through,
            cutoff.tail_method,
            stable
        });
    }
    if (rows.size() != full_q3_bcd_remaining::row_cutoffs.size()
        || polynomial_count != full_q3_bcd_remaining::polynomial_rows
        || rational_cap_count != full_q3_bcd_remaining::rational_cap_rows
        || directed_interval_count
            != full_q3_bcd_remaining::directed_interval_rows
        || largest_moment
            != full_q3_bcd_remaining::required_maximum_moment
        || expected_b_rank != 18 || expected_c_rank != 29
        || expected_d_rank != 31) {
        throw std::runtime_error("remaining-row ledger mismatch");
    }

    int completed = 0;
#ifdef _OPENMP
    const int threads = omp_get_max_threads();
#else
    const int threads = 1;
#endif
    std::cout << "FULL_Q3_REMAINING rows=" << rows.size()
              << " maximum=" << maximum
              << " full_certificate=" << (full_certificate ? 1 : 0)
              << " parallel_threads=" << threads
              << " maximum_concurrent_tasks=4" << std::endl;
    for (int moment = 2; moment <= maximum; ++moment) {
        std::vector<MomentTask> tasks;
        add_bc_task(tasks, rows, 'B', moment);
        add_bc_task(tasks, rows, 'C', moment);
        add_d_even_task(tasks, rows, moment);
        add_d_determinant_task(tasks, rows, moment);
        std::vector<std::vector<BigInt>> results(tasks.size());
        const int task_count = static_cast<int>(tasks.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int task_index = 0; task_index < task_count; ++task_index) {
            const MomentTask& task = tasks[static_cast<std::size_t>(task_index)];
            results[static_cast<std::size_t>(task_index)] =
                pieri_e2_coefficient_sums_reverse(
                    task.moment,
                    task.targets,
                    false,
                    "remaining " + task.label + " m_"
                        + std::to_string(task.moment)
                );
        }

        for (std::size_t task_index = 0; task_index < tasks.size(); ++task_index) {
            const MomentTask& task = tasks[task_index];
            const std::vector<BigInt>& sums = results[task_index];
            if (sums.size() != task.assignments.size()) {
                throw std::runtime_error("reverse-Pieri assignment ledger mismatch");
            }
            for (std::size_t color = 0; color < sums.size(); ++color) {
                RemainingRow& row = rows[task.assignments[color].row];
                BigInt& value =
                    row.moments[static_cast<std::size_t>(task.moment)];
                switch (task.assignments[color].operation) {
                    case MomentOperation::set:
                        value = sums[color];
                        break;
                    case MomentOperation::add:
                        value += sums[color];
                        break;
                    case MomentOperation::subtract:
                        value -= sums[color];
                        break;
                }
            }
            ++completed;
            if (progress) {
                std::cout << "FULL_Q3_REMAINING completed_task=" << completed
                          << " label=" << task.label
                          << " moment=" << task.moment << std::endl;
            }
        }
    }

    const RemainingRow* b2 = nullptr;
    const RemainingRow* c2 = nullptr;
    for (const RemainingRow& row : rows) {
        if (row.moments[0] != 1 || row.moments[1] != 0 || row.moments[2] != 1) {
            throw std::runtime_error("invalid invariant-tensor moment prefix");
        }
        const int checked_moment = std::min(maximum, row.moment_through);
        for (int moment = 0; moment <= checked_moment; ++moment) {
            const BigInt& value = row.moments[static_cast<std::size_t>(moment)];
            if (value < 0) throw std::runtime_error("negative finite-rank moment");
            if (moment <= row.stable_through
                && value != stable[static_cast<std::size_t>(moment)]) {
                throw std::runtime_error("finite-rank stable-prefix mismatch");
            }
            std::cout << "FULL_Q3_MOMENT row=" << row.family << '_' << row.rank
                      << " degree=" << moment << " value=" << value << '\n';
        }
        if (row.family == 'B' && row.rank == 2) b2 = &row;
        if (row.family == 'C' && row.rank == 2) c2 = &row;
    }
    if (b2 == nullptr || c2 == nullptr
        || b2->moment_through != c2->moment_through) {
        throw std::runtime_error("B2/C2 Lie-isomorphism moment mismatch");
    }
    for (int moment = 0;
         moment <= std::min(maximum, b2->moment_through);
         ++moment) {
        if (b2->moments[static_cast<std::size_t>(moment)]
            != c2->moments[static_cast<std::size_t>(moment)]) {
            throw std::runtime_error("B2/C2 Lie-isomorphism moment mismatch");
        }
    }

    std::size_t polynomial_checked = 0;
    std::size_t rational_cap_checked = 0;
    for (const RemainingRow& row : rows) {
        if (row.tail_method
            == full_q3_bcd_remaining::TailMethod::polynomial) {
            if (maximum >= row.moment_through) {
                verify_polynomial_tail(row);
                ++polynomial_checked;
            }
        } else if (row.tail_method
                   == full_q3_bcd_remaining::TailMethod::rational_cap) {
            verify_rational_cap_tail(row);
            ++rational_cap_checked;
        }
    }
    if (rational_cap_checked != full_q3_bcd_remaining::rational_cap_rows
        || (full_certificate
            && polynomial_checked
                != full_q3_bcd_remaining::polynomial_rows)) {
        throw std::runtime_error("exact-tail replay ledger mismatch");
    }

    std::size_t checked = 0;
    bool have_minimum = false;
    BigInt minimum;
    char minimum_family = 0;
    int minimum_rank = 0;
    int minimum_a = 0;
    int minimum_n = 0;
    for (const RemainingRow& row : rows) {
        std::size_t row_checked = 0;
        for (int total = first_odd_above(row.stable_through);
             total <= std::min(maximum, row.tail_onset - 2);
             total += 2) {
            for (int a = 2; 2 * a < total; ++a) {
                const int n = total - 2 * a;
                const BigInt value = hierarchy_value(row.moments, a, n);
                if (value <= 0) {
                    throw std::runtime_error(
                        std::string("nonpositive hierarchy value at ")
                        + row.family + '_' + std::to_string(row.rank)
                        + " a=" + std::to_string(a)
                        + " n=" + std::to_string(n)
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
                ++checked;
                ++row_checked;
            }
        }
        std::cout << "FULL_Q3_REMAINING_ROW row=" << row.family << '_'
                  << row.rank << " stable_through=" << row.stable_through
                  << " checked_through="
                  << std::min(maximum, row.tail_onset - 2)
                  << " moment_through="
                  << std::min(maximum, row.moment_through)
                  << " tail_onset=" << row.tail_onset
                  << " pairs=" << row_checked << '\n';
    }
    if (!have_minimum || checked == 0U) {
        throw std::runtime_error("empty remaining hierarchy audit");
    }
    if (full_certificate
        && checked != full_q3_bcd_remaining::full_residual_pairs) {
        throw std::runtime_error("remaining hierarchy scope-count mismatch");
    }
    std::cout << "FULL_Q3_REMAINING pairs=" << checked
              << " minimum=" << minimum_family << '_' << minimum_rank
              << " a=" << minimum_a << " n=" << minimum_n
              << " value=" << minimum << '\n';
    if (full_certificate) {
        std::cout << "FULL_Q3_BCD_REMAINING VERIFICATION: ALL PASS\n";
    } else {
        std::cout << "FULL_Q3_BCD_REMAINING DIAGNOSTIC PREFIX: PASS\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        bool progress = false;
        int maximum = full_q3_bcd_remaining::required_maximum_moment;
        std::string stable_path = "references/oeis_A002137_stable.txt";
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--progress") {
                progress = true;
            } else if (argument == "--maximum" && index + 1 < argc) {
                maximum = std::stoi(argv[++index]);
            } else if (argument == "--stable-moments" && index + 1 < argc) {
                stable_path = argv[++index];
            } else {
                std::cerr << "usage: " << argv[0]
                          << " [--progress] [--maximum N]"
                          << " [--stable-moments PATH]\n";
                return EXIT_FAILURE;
            }
        }
        run_remaining_certificate(stable_path, maximum, progress);
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "FULL_Q3_BCD_REMAINING VERIFICATION: FAIL: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
