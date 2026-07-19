#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr long double LOG10_E = 0.434294481903251827651128918916605082L;
constexpr long double LOG10_2 = 0.301029995663981195213738894724493027L;

struct RankResult {
    int rank = 0;
    int onset = 0;
    int slack_target = 0;
    long double onset_margin = -std::numeric_limits<long double>::infinity();
    long double slack_margin = -std::numeric_limits<long double>::infinity();
    long double radius = 0.0L;
    long double c_push = 0.0L;
};

struct FormulaParams {
    long double radius_slope = 40.0L;
    long double radius_offset = -300.0L;
    long double radius_span_slope = 40.0L;
    long double radius_span_offset = -1180.0L;
    long double radius_unit = 0.0000295386650000L;
};

int next_odd_at_least(int n) {
    return (n % 2 == 0) ? n + 1 : n;
}

long double log10_factorial(int n) {
    long double total = 0.0L;
    for (int k = 2; k <= n; ++k) total += std::log10(static_cast<long double>(k));
    return total;
}

long double log10_half_integer_density_factor(int rank, int dimension) {
    if (rank % 2 != 1 || dimension % 2 != 1) {
        return -std::numeric_limits<long double>::infinity();
    }
    const int half_index = (dimension + 1) / 2;
    const int pi_power = (rank + 1) / 2;
    return std::log10(7.0L / 5.0L)
        + static_cast<long double>(half_index) * std::log10(4.0L)
        + log10_factorial(half_index)
        - static_cast<long double>(pi_power) * std::log10(44.0L / 7.0L)
        - log10_factorial(2 * half_index);
}

long double log10_density_factor(int rank, int dimension) {
    if (rank % 2 == 1 && dimension % 2 == 1) {
        return log10_half_integer_density_factor(rank, dimension);
    }
    if (rank % 2 == 0 && dimension % 2 == 0) {
        return -static_cast<long double>(rank / 2) * std::log10(44.0L / 7.0L)
            - log10_factorial(dimension / 2);
    }
    return -std::numeric_limits<long double>::infinity();
}

long double log10_degree_factor(int rank) {
    long double total = 0.0L;
    for (int j = 1; j <= rank; ++j) total += log10_factorial(2 * j);
    return total;
}

long double log10_weyl_order(int rank) {
    return static_cast<long double>(rank) * LOG10_2 + log10_factorial(rank);
}

long double sine_exponent(long double radius, int kappa) {
    const long double k = static_cast<long double>(kappa);
    const long double q = 2.0L * radius / (39.0L * k);
    if (!(q > 0.0L && q < 1.0L)) return std::numeric_limits<long double>::infinity();
    return radius / 12.0L
        + radius * radius / (1440.0L * k)
        + radius * radius * radius / (k * k * 90720.0L * (1.0L - q));
}

long double d_generated_radius(int rank, const FormulaParams& params) {
    const long double radius_min = params.radius_slope * rank + params.radius_offset;
    const long double radius_span =
        params.radius_span_slope * rank + params.radius_span_offset;
    return radius_min + params.radius_unit * radius_span;
}

long double d_interval_margin(
    int rank,
    int target_n,
    const FormulaParams& params,
    long double& c_push_out
) {
    const int dimension = rank * (2 * rank - 1);
    const int kappa = 2 * rank - 2;
    const int c_value = (rank % 2 == 0) ? rank : rank - 2;

    constexpr long double central_t = 1.64578692310241648L;
    constexpr long double quartic_r = 3.06396669928381904L;
    constexpr long double quartic_s = 0.317259520411341256L;
    constexpr long double cheb_cutoff = 80.0L;
    constexpr long double log10_neg = -15.0L;

    const long double radius = d_generated_radius(rank, params);
    const long double x_floor = static_cast<long double>(dimension) - radius;
    const long double c_push = x_floor - central_t;
    c_push_out = c_push;
    if (!(c_push > 2.0L * c_value)) {
        return -std::numeric_limits<long double>::infinity();
    }

    const long double base_adjust =
        LOG10_2 + log10_factorial(rank) - log10_factorial(2 * rank);
    const long double base_log =
        -log10_weyl_order(rank)
        + log10_degree_factor(rank)
        + log10_density_factor(rank, dimension)
        + base_adjust;

    const long double sine_log = -sine_exponent(radius, kappa) * LOG10_E;
    const long double radial_log =
        ((static_cast<long double>(dimension - 1) / 2.0L) + 0.5L)
        * std::log10(radius / (2.0L * kappa));
    const long double ball_log = base_log + sine_log + radial_log;

    const long double quartic_moment =
        6.0L - 2.0L * (quartic_r + quartic_s)
        + (quartic_r * quartic_r + quartic_s * quartic_s
           + 4.0L * quartic_r * quartic_s)
        + quartic_r * quartic_r * quartic_s * quartic_s;
    const long double denom =
        (central_t + quartic_r) * (central_t + quartic_r)
        * (central_t + quartic_s) * (central_t + quartic_s);
    if (!(denom > 0.0L)) return -std::numeric_limits<long double>::infinity();

    const long double tail_upper =
        (x_floor + central_t) * (x_floor + central_t) * quartic_moment / denom;
    const long double push_weight = x_floor * x_floor + 1.0L - tail_upper;
    const long double deriv =
        2.0L * x_floor - 2.0L * quartic_moment * (x_floor + central_t) / denom;
    const long double second = 2.0L - 2.0L * quartic_moment / denom;
    if (!(push_weight > 0.0L && deriv > 0.0L && second > 0.0L)) {
        return -std::numeric_limits<long double>::infinity();
    }

    const long double weighted_log = ball_log + std::log10(push_weight);
    const long double target_log =
        std::max(
            log10_neg + static_cast<long double>(target_n) * std::log10(2.0L * c_value),
            LOG10_2 + static_cast<long double>(target_n) * std::log10(cheb_cutoff))
        + LOG10_2 - static_cast<long double>(target_n) * std::log10(c_push);
    return weighted_log - target_log;
}

RankResult scan_rank(int rank, int slack, const FormulaParams& params) {
    int lo = 63;
    int hi = next_odd_at_least(std::max(201, rank * rank));
    long double c_push = 0.0L;
    while (!(d_interval_margin(rank, hi, params, c_push) > 0.0L)) {
        hi = next_odd_at_least(2 * hi + 1);
        if (hi > 1000000) {
            RankResult failed;
            failed.rank = rank;
            failed.onset = 0;
            failed.radius = d_generated_radius(rank, params);
            failed.c_push = c_push;
            return failed;
        }
    }
    lo = next_odd_at_least(lo);
    if (d_interval_margin(rank, lo, params, c_push) > 0.0L) {
        hi = lo;
    } else {
        int false_n = lo;
        int true_n = hi;
        while (false_n + 2 < true_n) {
            int mid = false_n + ((true_n - false_n) / 4) * 2;
            if (mid <= false_n) mid = false_n + 2;
            if (mid >= true_n) mid = true_n - 2;
            if (d_interval_margin(rank, mid, params, c_push) > 0.0L) {
                true_n = mid;
            } else {
                false_n = mid;
            }
        }
        hi = true_n;
    }

    RankResult out;
    out.rank = rank;
    out.onset = hi;
    out.slack_target = next_odd_at_least(hi + slack);
    out.radius = d_generated_radius(rank, params);
    out.onset_margin = d_interval_margin(rank, out.onset, params, out.c_push);
    out.slack_margin = d_interval_margin(rank, out.slack_target, params, c_push);
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    int rank_lo = 98;
    int rank_hi = 279;
    int slack = 32;
    bool progress = false;
    FormulaParams params;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--rank-lo" && i + 1 < argc) {
            rank_lo = std::atoi(argv[++i]);
        } else if (arg == "--rank-hi" && i + 1 < argc) {
            rank_hi = std::atoi(argv[++i]);
        } else if (arg == "--slack" && i + 1 < argc) {
            slack = std::atoi(argv[++i]);
        } else if (arg == "--radius-slope" && i + 1 < argc) {
            params.radius_slope = std::strtold(argv[++i], nullptr);
        } else if (arg == "--radius-offset" && i + 1 < argc) {
            params.radius_offset = std::strtold(argv[++i], nullptr);
        } else if (arg == "--radius-span-slope" && i + 1 < argc) {
            params.radius_span_slope = std::strtold(argv[++i], nullptr);
        } else if (arg == "--radius-span-offset" && i + 1 < argc) {
            params.radius_span_offset = std::strtold(argv[++i], nullptr);
        } else if (arg == "--radius-unit" && i + 1 < argc) {
            params.radius_unit = std::strtold(argv[++i], nullptr);
        } else if (arg == "--progress") {
            progress = true;
        } else {
            std::cerr << "usage: " << argv[0]
                      << " [--rank-lo N] [--rank-hi N] [--slack N]"
                      << " [--radius-slope X] [--radius-offset X]"
                      << " [--radius-span-slope X] [--radius-span-offset X]"
                      << " [--radius-unit X] [--progress]\n";
            return 2;
        }
    }
    if (rank_lo < 48 || rank_hi < rank_lo) {
        std::cerr << "invalid D-rank interval\n";
        return 2;
    }

    std::vector<RankResult> results(static_cast<std::size_t>(rank_hi - rank_lo + 1));
    std::cout << "D interval fixed-template direct-tail diagnostic\n"
              << "rank_lo=" << rank_lo << " rank_hi=" << rank_hi
              << " slack=" << slack
              << " radius_slope=" << static_cast<double>(params.radius_slope)
              << " radius_offset=" << static_cast<double>(params.radius_offset)
              << " radius_span_slope=" << static_cast<double>(params.radius_span_slope)
              << " radius_span_offset=" << static_cast<double>(params.radius_span_offset)
              << " radius_unit=" << static_cast<double>(params.radius_unit)
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;

#pragma omp parallel for schedule(dynamic, 1)
    for (int rank = rank_lo; rank <= rank_hi; ++rank) {
        RankResult result = scan_rank(rank, slack, params);
        results[static_cast<std::size_t>(rank - rank_lo)] = result;
        if (progress) {
#pragma omp critical
            {
                std::cout << "rank_done D_" << rank
                          << " onset=" << result.onset
                          << " margin=" << static_cast<double>(result.onset_margin)
                          << " slack_margin=" << static_cast<double>(result.slack_margin)
                          << std::endl;
            }
        }
    }

    RankResult worst_onset = results.front();
    RankResult worst_slack = results.front();
    int failures = 0;
    for (const RankResult& result : results) {
        if (result.onset == 0) {
            ++failures;
            continue;
        }
        if (result.onset_margin < worst_onset.onset_margin) worst_onset = result;
        if (result.slack_margin < worst_slack.slack_margin) worst_slack = result;
    }

    std::cout << std::setprecision(18)
              << "rank\tonset\tmargin_log10\tslack_target\tslack_margin_log10"
              << "\tradius\tc_push\n";
    for (const RankResult& result : results) {
        std::cout << "D_" << result.rank << '\t'
                  << result.onset << '\t'
                  << static_cast<double>(result.onset_margin) << '\t'
                  << result.slack_target << '\t'
                  << static_cast<double>(result.slack_margin) << '\t'
                  << static_cast<double>(result.radius) << '\t'
                  << static_cast<double>(result.c_push) << '\n';
    }
    std::cout << "SUMMARY failures=" << failures
              << " worst_onset=D_" << worst_onset.rank
              << " onset=" << worst_onset.onset
              << " margin_log10=" << static_cast<double>(worst_onset.onset_margin)
              << " worst_slack=D_" << worst_slack.rank
              << " target=" << worst_slack.slack_target
              << " margin_log10=" << static_cast<double>(worst_slack.slack_margin)
              << std::endl;
    return failures == 0 ? 0 : 1;
}
