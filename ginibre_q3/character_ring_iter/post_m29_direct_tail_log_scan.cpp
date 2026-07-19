#include <omp.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr long double LOG10_E = 0.434294481903251827651128918916605082L;
constexpr long double LOG10_2 = 0.301029995663981195213738894724493027L;

struct ChebBound {
    int degree;
    int scale;
    long double cutoff;
    long double log10_neg;
    int majorant_degree;
};

struct ScanConfig {
    std::string label;
    int rank;
    int dimension;
    int kappa;
    int c_value;
    int target_n;
    long double radius_min;
    long double radius_max;
    long double radius_step;
    long double central_min;
    long double central_max;
    long double central_step;
    long double q_min;
    long double q_max;
    long double q_step;
    std::vector<ChebBound> cheb_bounds;
    long double base_adjust_log10 = 0.0L;
};

struct Candidate {
    long double margin = -std::numeric_limits<long double>::infinity();
    long double radius = 0.0L;
    long double central_t = 0.0L;
    long double quartic_r = 0.0L;
    long double quartic_s = 0.0L;
    long double c_push = 0.0L;
    ChebBound bound{};
};

long double log10_factorial(int n) {
    long double total = 0.0L;
    for (int k = 2; k <= n; ++k) total += std::log10(static_cast<long double>(k));
    return total;
}

long double log10_half_integer_density_factor(int rank, int dimension) {
    // Matches the conservative rational lower bound in
    // classical_tail_constants.py:half_integer_density_factor.
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

[[maybe_unused]] Candidate scan_config(const ScanConfig& cfg, bool print_progress) {
    const long double base_log =
        -log10_weyl_order(cfg.rank)
        + log10_degree_factor(cfg.rank)
        + log10_density_factor(cfg.rank, cfg.dimension)
        + cfg.base_adjust_log10;

    const std::int64_t radius_count =
        static_cast<std::int64_t>(std::floor((cfg.radius_max - cfg.radius_min) / cfg.radius_step)) + 1;
    const std::int64_t central_count =
        static_cast<std::int64_t>(std::floor((cfg.central_max - cfg.central_min) / cfg.central_step)) + 1;
    const std::int64_t q_count =
        static_cast<std::int64_t>(std::floor((cfg.q_max - cfg.q_min) / cfg.q_step)) + 1;
    const std::int64_t total_outer = radius_count * central_count;

    Candidate global_best;
    std::cout << "case " << cfg.label
              << " target_n=" << cfg.target_n
              << " OpenMP threads=" << omp_get_max_threads()
              << " outer=" << total_outer
              << " quartic_pairs=" << (q_count * q_count)
              << std::endl;

#pragma omp parallel
    {
        Candidate local_best;
#pragma omp for schedule(dynamic, 16)
        for (std::int64_t outer = 0; outer < total_outer; ++outer) {
            const std::int64_t radius_index = outer / central_count;
            const std::int64_t central_index = outer % central_count;
            const long double radius = cfg.radius_min + cfg.radius_step * radius_index;
            const long double central_t = cfg.central_min + cfg.central_step * central_index;
            const long double q_bound = 2.0L * radius / (39.0L * cfg.kappa);
            if (!(q_bound > 0.0L && q_bound < 1.0L)) continue;

            const long double x_floor = static_cast<long double>(cfg.dimension) - radius;
            const long double c_push = x_floor - central_t;
            if (!(c_push > 2.0L * cfg.c_value)) continue;

            const long double sine_log = -sine_exponent(radius, cfg.kappa) * LOG10_E;
            const long double radial_log =
                ((static_cast<long double>(cfg.dimension - 1) / 2.0L) + 0.5L)
                * std::log10(radius / (2.0L * cfg.kappa));
            const long double ball_log = base_log + sine_log + radial_log;

            for (std::int64_t qr_index = 0; qr_index < q_count; ++qr_index) {
                const long double quartic_r = cfg.q_min + cfg.q_step * qr_index;
                for (std::int64_t qs_index = 0; qs_index < q_count; ++qs_index) {
                    const long double quartic_s = cfg.q_min + cfg.q_step * qs_index;
                    const long double quartic_moment =
                        6.0L - 2.0L * (quartic_r + quartic_s)
                        + (quartic_r * quartic_r + quartic_s * quartic_s
                           + 4.0L * quartic_r * quartic_s)
                        + quartic_r * quartic_r * quartic_s * quartic_s;
                    const long double denom =
                        (central_t + quartic_r) * (central_t + quartic_r)
                        * (central_t + quartic_s) * (central_t + quartic_s);
                    if (!(denom > 0.0L)) continue;
                    const long double tail_upper =
                        (x_floor + central_t) * (x_floor + central_t)
                        * quartic_moment / denom;
                    const long double push_weight = x_floor * x_floor + 1.0L - tail_upper;
                    const long double deriv =
                        2.0L * x_floor
                        - 2.0L * quartic_moment * (x_floor + central_t) / denom;
                    const long double second = 2.0L - 2.0L * quartic_moment / denom;
                    if (!(push_weight > 0.0L && deriv > 0.0L && second > 0.0L)) continue;

                    const long double weighted_log = ball_log + std::log10(push_weight);
                    for (const ChebBound& bound : cfg.cheb_bounds) {
                        const long double target_log =
                            std::max(
                                bound.log10_neg
                                    + static_cast<long double>(cfg.target_n)
                                        * std::log10(2.0L * cfg.c_value),
                                LOG10_2
                                    + static_cast<long double>(cfg.target_n)
                                        * std::log10(bound.cutoff))
                            + LOG10_2
                            - static_cast<long double>(cfg.target_n) * std::log10(c_push);
                        const long double margin = weighted_log - target_log;
                        if (margin > local_best.margin) {
                            local_best.margin = margin;
                            local_best.radius = radius;
                            local_best.central_t = central_t;
                            local_best.quartic_r = quartic_r;
                            local_best.quartic_s = quartic_s;
                            local_best.c_push = c_push;
                            local_best.bound = bound;
                        }
                    }
                }
            }
            if (print_progress && outer % 10000 == 0) {
#pragma omp critical
                {
                    std::cout << "progress outer=" << outer << "/" << total_outer
                              << " local_best=" << static_cast<double>(local_best.margin)
                              << std::endl;
                }
            }
        }
#pragma omp critical
        {
            if (local_best.margin > global_best.margin) global_best = local_best;
        }
    }
    return global_best;
}

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

long double unit_from_u64(std::uint64_t x) {
    return static_cast<long double>(x >> 11)
        / static_cast<long double>(1ULL << 53);
}

long double sample_range(std::uint64_t seed, long double lo, long double hi) {
    return lo + (hi - lo) * unit_from_u64(seed);
}

Candidate sample_config(const ScanConfig& cfg, std::uint64_t samples, bool print_progress) {
    const long double base_log =
        -log10_weyl_order(cfg.rank)
        + log10_degree_factor(cfg.rank)
        + log10_density_factor(cfg.rank, cfg.dimension)
        + cfg.base_adjust_log10;

    Candidate global_best;
    std::cout << "case " << cfg.label
              << " target_n=" << cfg.target_n
              << " OpenMP threads=" << omp_get_max_threads()
              << " samples=" << samples
              << std::endl;

#pragma omp parallel
    {
        Candidate local_best;
#pragma omp for schedule(dynamic, 4096)
        for (std::int64_t signed_index = 0;
             signed_index < static_cast<std::int64_t>(samples);
             ++signed_index) {
            const std::uint64_t index = static_cast<std::uint64_t>(signed_index);
            const long double radius =
                sample_range(splitmix64(4 * index + 1), cfg.radius_min, cfg.radius_max);
            const long double central_t =
                sample_range(splitmix64(4 * index + 2), cfg.central_min, cfg.central_max);
            const long double quartic_r =
                sample_range(splitmix64(4 * index + 3), cfg.q_min, cfg.q_max);
            const long double quartic_s =
                sample_range(splitmix64(4 * index + 4), cfg.q_min, cfg.q_max);

            const long double q_bound = 2.0L * radius / (39.0L * cfg.kappa);
            if (!(q_bound > 0.0L && q_bound < 1.0L)) continue;
            const long double x_floor = static_cast<long double>(cfg.dimension) - radius;
            const long double c_push = x_floor - central_t;
            if (!(c_push > 2.0L * cfg.c_value)) continue;

            const long double sine_log = -sine_exponent(radius, cfg.kappa) * LOG10_E;
            const long double radial_log =
                ((static_cast<long double>(cfg.dimension - 1) / 2.0L) + 0.5L)
                * std::log10(radius / (2.0L * cfg.kappa));
            const long double ball_log = base_log + sine_log + radial_log;

            const long double quartic_moment =
                6.0L - 2.0L * (quartic_r + quartic_s)
                + (quartic_r * quartic_r + quartic_s * quartic_s
                   + 4.0L * quartic_r * quartic_s)
                + quartic_r * quartic_r * quartic_s * quartic_s;
            const long double denom =
                (central_t + quartic_r) * (central_t + quartic_r)
                * (central_t + quartic_s) * (central_t + quartic_s);
            if (!(denom > 0.0L)) continue;
            const long double tail_upper =
                (x_floor + central_t) * (x_floor + central_t)
                * quartic_moment / denom;
            const long double push_weight = x_floor * x_floor + 1.0L - tail_upper;
            const long double deriv =
                2.0L * x_floor
                - 2.0L * quartic_moment * (x_floor + central_t) / denom;
            const long double second = 2.0L - 2.0L * quartic_moment / denom;
            if (!(push_weight > 0.0L && deriv > 0.0L && second > 0.0L)) continue;

            const long double weighted_log = ball_log + std::log10(push_weight);
            for (const ChebBound& bound : cfg.cheb_bounds) {
                const long double target_log =
                    std::max(
                        bound.log10_neg
                            + static_cast<long double>(cfg.target_n)
                                * std::log10(2.0L * cfg.c_value),
                        LOG10_2
                            + static_cast<long double>(cfg.target_n)
                                * std::log10(bound.cutoff))
                    + LOG10_2
                    - static_cast<long double>(cfg.target_n) * std::log10(c_push);
                const long double margin = weighted_log - target_log;
                if (margin > local_best.margin) {
                    local_best.margin = margin;
                    local_best.radius = radius;
                    local_best.central_t = central_t;
                    local_best.quartic_r = quartic_r;
                    local_best.quartic_s = quartic_s;
                    local_best.c_push = c_push;
                    local_best.bound = bound;
                }
            }

            if (print_progress && index % 1000000ULL == 0ULL) {
#pragma omp critical
                {
                    std::cout << "sample " << index << "/" << samples
                              << " local_best=" << static_cast<double>(local_best.margin)
                              << std::endl;
                }
            }
        }
#pragma omp critical
        {
            if (local_best.margin > global_best.margin) global_best = local_best;
        }
    }
    return global_best;
}

void print_best(const std::string& label, const Candidate& best) {
    std::cout << std::setprecision(18)
              << "best " << label
              << " margin_log10=" << static_cast<double>(best.margin)
              << " radius=" << static_cast<double>(best.radius)
              << " central_t=" << static_cast<double>(best.central_t)
              << " quartic_r=" << static_cast<double>(best.quartic_r)
              << " quartic_s=" << static_cast<double>(best.quartic_s)
              << " c_push=" << static_cast<double>(best.c_push)
              << " cheb_degree=" << best.bound.degree
              << " cheb_scale=" << best.bound.scale
              << " cutoff=" << static_cast<double>(best.bound.cutoff)
              << " majorant_degree=" << best.bound.majorant_degree
              << " log10_neg=" << static_cast<double>(best.bound.log10_neg)
              << std::endl;
    std::cout << (best.margin > 0.0L ? "RESULT CANDIDATE" : "RESULT NO_CANDIDATE")
              << std::endl;
}

int next_odd_at_least(int n) {
    return (n % 2 == 0) ? n + 1 : n;
}

int prev_odd_at_most(int n) {
    return (n % 2 == 0) ? n - 1 : n;
}

int binary_search_onset(
    ScanConfig cfg,
    int lo_input,
    int hi_input,
    std::uint64_t samples,
    bool progress
) {
    const int lo = next_odd_at_least(lo_input);
    const int hi = prev_odd_at_most(hi_input);
    if (lo <= 0 || hi < lo) {
        std::cerr << "invalid binary-search range after odd normalization: "
                  << lo_input << ".." << hi_input << std::endl;
        return 2;
    }

    std::cout << "binary_search case=" << cfg.label
              << " lo=" << lo
              << " hi=" << hi
              << " samples_per_probe=" << samples
              << std::endl;

    auto run_probe = [&](int target_n) {
        cfg.target_n = target_n;
        std::cout << "binary_probe target_n=" << target_n << std::endl;
        Candidate best = sample_config(cfg, samples, progress);
        print_best(cfg.label, best);
        return best;
    };

    Candidate lo_best = run_probe(lo);
    if (lo_best.margin > 0.0L) {
        std::cout << "BINARY_ONSET target_n=" << lo << std::endl;
        print_best(cfg.label + "-binary-onset", lo_best);
        return 0;
    }

    Candidate hi_best = run_probe(hi);
    if (!(hi_best.margin > 0.0L)) {
        std::cout << "RESULT NO_ONSET_IN_RANGE" << std::endl;
        return 1;
    }

    int false_n = lo;
    int true_n = hi;
    Candidate true_best = hi_best;
    while (false_n + 2 < true_n) {
        int mid = false_n + ((true_n - false_n) / 4) * 2;
        if (mid <= false_n) mid = false_n + 2;
        if (mid >= true_n) mid = true_n - 2;
        Candidate mid_best = run_probe(mid);
        if (mid_best.margin > 0.0L) {
            true_n = mid;
            true_best = mid_best;
        } else {
            false_n = mid;
        }
        std::cout << "binary_interval false_below_or_at=" << false_n
                  << " true_at_or_above=" << true_n
                  << std::endl;
    }

    std::cout << "BINARY_ONSET target_n=" << true_n << std::endl;
    print_best(cfg.label + "-binary-onset", true_best);
    return 0;
}

ScanConfig c17_config() {
    // The Chebyshev constants are exact-replay outputs compressed to log10
    // form from the local C_17 delta window Delta_36..Delta_39.  The scan is
    // diagnostic; any positive candidate must still be replayed exactly.
    return {
        "C17-post-m29-direct-tail",
        17,
        17 * (2 * 17 + 1),
        18,
        17,
        63,
        150.0L,
        450.0L,
        0.05L,
        0.05L,
        8.0L,
        0.01L,
        0.05L,
        16.0L,
        0.01L,
        {
            {9, 14500, 22.5L, -11.030059343790285L, 36},
        },
    };
}

ScanConfig b13_exact54_config() {
    // Diagnostic B_13 scan using the exact local delta window
    // Delta_28..Delta_54.  A positive candidate at target_n=83 would combine
    // with the finite Chain bridge through n=81 and avoid the active
    // Delta_55..Delta_62 computations.  Any positive candidate must still be
    // replayed exactly in classical_tail_constants.py.
    return {
        "B13-post-m29-direct-tail-exact54",
        13,
        13 * (2 * 13 + 1),
        25,
        13,
        83,
        80.0L,
        320.0L,
        0.05L,
        0.05L,
        8.0L,
        0.01L,
        0.05L,
        16.0L,
        0.01L,
        {
            {5, 6000, 18.5L, -6.695190058371046860L, 20},
            {6, 6000, 18.5L, -7.613237175329700790L, 24},
            {7, 8000, 18.5L, -8.404111992361151806L, 28},
            {8, 8000, 18.5L, -9.257272424776203934L, 32},
            {9, 10000, 18.5L, -9.855985406262590232L, 36},
            {10, 10000, 18.5L, -10.672026549881664437L, 40},
            {11, 10000, 18.5L, -11.231280535439637447L, 44},
            {12, 12000, 18.5L, -11.909867817497868714L, 48},
            {13, 12000, 18.5L, -12.502699299334224747L, 52},
        },
    };
}

ScanConfig b13_exact57_config() {
    ScanConfig cfg = b13_exact54_config();
    cfg.label = "B13-post-m29-direct-tail-exact57";
    cfg.cheb_bounds.push_back({14, 12000, 18.5L, -12.749889658256577718L, 56});
    cfg.cheb_bounds.push_back({14, 14000, 18.5L, -13.009845176748612516L, 56});
    cfg.cheb_bounds.push_back({14, 12000, 20.5L, -13.691379993675340643L, 56});
    cfg.cheb_bounds.push_back({14, 14000, 20.5L, -13.917154612706397643L, 56});
    cfg.cheb_bounds.push_back({14, 12000, 22.5L, -14.566922635202487868L, 56});
    cfg.cheb_bounds.push_back({14, 14000, 22.5L, -14.762539205748723248L, 56});
    cfg.cheb_bounds.push_back({14, 12000, 24.5L, -15.385756505330945743L, 56});
    cfg.cheb_bounds.push_back({14, 14000, 24.5L, -15.554549889854243361L, 56});
    return cfg;
}

ScanConfig b14_exact39_config() {
    // Diagnostic B_14 scan using the exact local delta window
    // Delta_30..Delta_39 from m30_BC13_18_pair_range_exact39.log.  A positive
    // candidate at target_n=63 would close the B_14 all-later tail from the
    // first post-m=29 hit, after exact rational replay.
    return {
        "B14-post-m29-direct-tail-exact39",
        14,
        14 * (2 * 14 + 1),
        27,
        14,
        63,
        100.0L,
        380.0L,
        0.05L,
        0.05L,
        8.0L,
        0.01L,
        0.05L,
        18.0L,
        0.01L,
        {
            {5, 10000, 18.0L, -6.101033912526418L, 20},
            {6, 10000, 18.0L, -7.138595699550926L, 24},
            {7, 10000, 18.0L, -8.128723706361448L, 28},
            {8, 10000, 18.0L, -9.039596833937026L, 32},
            {9, 10000, 18.0L, -9.752027785102982L, 36},
            {10, 12000, 18.0L, -10.416901268595012L, 40},
            {11, 12000, 18.0L, -11.043698000493166L, 44},
            {12, 14000, 18.0L, -11.645466295440215L, 48},
            {13, 14000, 18.0L, -12.206013672695605L, 52},
            {14, 16000, 18.0L, -12.758498169651077L, 56},
            {5, 10000, 22.0L, -6.805543607377039L, 20},
            {6, 10000, 22.0L, -7.984007334301111L, 24},
            {7, 10000, 22.0L, -9.115037280252331L, 28},
            {8, 10000, 22.0L, -10.166812346955439L, 32},
            {9, 10000, 22.0L, -11.020145237248698L, 36},
            {10, 12000, 22.0L, -11.775693890093137L, 40},
            {11, 12000, 22.0L, -12.538369884141105L, 44},
            {12, 14000, 22.0L, -13.221769411123262L, 48},
            {13, 14000, 22.0L, -13.913675381352249L, 52},
            {14, 16000, 22.0L, -14.540158395205495L, 56},
            {5, 10000, 26.0L, -7.414791723763578L, 20},
            {6, 10000, 26.0L, -8.715105074138094L, 24},
            {7, 10000, 26.0L, -9.967984643397614L, 28},
            {8, 10000, 26.0L, -11.141609333407203L, 32},
            {9, 10000, 26.0L, -12.116791847006951L, 36},
            {10, 12000, 26.0L, -12.955492509511359L, 40},
            {11, 12000, 26.0L, -13.836148365501145L, 44},
            {12, 14000, 26.0L, -14.595254973152066L, 48},
            {13, 14000, 26.0L, -15.401618073550125L, 52},
            {14, 16000, 26.0L, -16.097412576522174L, 56},
            {5, 10000, 30.0L, -7.952455949737946L, 20},
            {6, 10000, 30.0L, -9.360302145347525L, 24},
            {7, 10000, 30.0L, -10.720714559808989L, 28},
            {8, 10000, 30.0L, -12.001872095020218L, 32},
            {9, 10000, 30.0L, -13.084587453821584L, 36},
            {10, 12000, 30.0L, -14.000054800010702L, 40},
            {11, 12000, 30.0L, -14.985166885050432L, 44},
            {12, 14000, 30.0L, -15.814810276168160L, 48},
            {13, 14000, 30.0L, -16.722802985150892L, 52},
            {14, 16000, 30.0L, -17.483704371319860L, 56},
            {5, 10000, 35.0L, -8.547330845035674L, 20},
            {6, 10000, 35.0L, -10.074152019717062L, 24},
            {7, 10000, 35.0L, -11.553539413240202L, 28},
            {8, 10000, 35.0L, -12.953671927513028L, 32},
            {9, 10000, 35.0L, -14.155362265376013L, 36},
            {10, 12000, 35.0L, -15.159183685819158L, 40},
            {11, 12000, 35.0L, -16.260208659439741L, 44},
            {12, 14000, 35.0L, -17.171720969759107L, 48},
            {13, 14000, 35.0L, -18.192789569874421L, 52},
            {14, 16000, 35.0L, -19.029831783845950L, 56},
            {5, 10000, 42.0L, -9.271208401482426L, 20},
            {6, 10000, 42.0L, -10.942805087456556L, 24},
            {7, 10000, 42.0L, -12.566967992269639L, 28},
            {8, 10000, 42.0L, -14.111876017832387L, 32},
            {9, 10000, 42.0L, -15.458341866985265L, 36},
            {10, 12000, 42.0L, -16.574102894873349L, 40},
            {11, 12000, 42.0L, -17.816619789399340L, 44},
            {12, 14000, 42.0L, -18.832789732550452L, 48},
            {13, 14000, 42.0L, -19.992280729565039L, 52},
            {14, 16000, 42.0L, -20.927456411445178L, 56},
        },
    };
}

ScanConfig c18_exact39_config() {
    // Diagnostic C_18 direct-tail scan using the exact local delta window
    // Delta_38..Delta_39 from m30_BC13_18_pair_range_exact39.log.  A positive
    // candidate at target_n=63 would close C_18 without new delta jobs, after
    // exact replay.
    return {
        "C18-post-m29-direct-tail-exact39",
        18,
        18 * (2 * 18 + 1),
        19,
        18,
        63,
        120.0L,
        365.0L,
        0.05L,
        0.05L,
        8.0L,
        0.01L,
        0.05L,
        16.0L,
        0.01L,
        {
            {5, 11000, 20.5L, -7.029657855760262919L, 20},
            {6, 13000, 20.5L, -7.980463122222140271L, 24},
            {7, 15000, 20.5L, -8.848280514062807356L, 28},
            {8, 15000, 20.5L, -9.687217502391845869L, 32},
            {9, 18000, 20.5L, -10.467358077194163002L, 36},
            {5, 11000, 22.5L, -7.379462899673150389L, 20},
            {6, 13000, 22.5L, -8.390129251429463864L, 24},
            {7, 13000, 22.5L, -9.319080928118474105L, 28},
            {8, 15000, 22.5L, -10.220900254147082364L, 32},
            {9, 18000, 22.5L, -11.048303064729054768L, 36},
            {5, 11000, 25.0L, -7.781335792605759138L, 20},
            {6, 13000, 25.0L, -8.861575115516878043L, 24},
            {7, 13000, 25.0L, -9.869101102887270827L, 28},
            {8, 15000, 25.0L, -10.836019353192511971L, 32},
            {9, 18000, 25.0L, -11.719296695714589873L, 36},
            {5, 11000, 30.0L, -8.490976596434499868L, 20},
            {6, 13000, 30.0L, -9.696013978803513567L, 24},
            {7, 13000, 30.0L, -10.842613110055083325L, 28},
            {8, 15000, 30.0L, -11.927091472703111208L, 32},
            {9, 18000, 30.0L, -12.912948921413303083L, 36},
        },
    };
}

ScanConfig c18_exact47_config() {
    // Diagnostic C_18 direct-tail scan using the exact local delta window
    // Delta_38..Delta_47.  The goal is to test whether the direct tail can
    // start before the m=38 Chain bridge step.
    return {
        "C18-post-m29-direct-tail-exact47",
        18,
        18 * (2 * 18 + 1),
        19,
        18,
        79,
        120.0L,
        365.0L,
        0.05L,
        0.05L,
        8.0L,
        0.01L,
        0.05L,
        16.0L,
        0.01L,
        {
            {5, 11000, 20.5L, -7.029657855760262919L, 20},
            {6, 13000, 20.5L, -7.980463122222140271L, 24},
            {7, 15000, 20.5L, -8.848280514062807356L, 28},
            {8, 15000, 20.5L, -9.687217502391845869L, 32},
            {9, 18000, 20.5L, -10.467358077194163002L, 36},
            {10, 18000, 20.5L, -11.136346641454281325L, 40},
            {11, 22000, 20.5L, -11.839459236318006674L, 44},
            {5, 11000, 22.5L, -7.379462899673150389L, 20},
            {6, 13000, 22.5L, -8.390129251429463864L, 24},
            {7, 13000, 22.5L, -9.319080928118474105L, 28},
            {8, 15000, 22.5L, -10.220900254147082364L, 32},
            {9, 18000, 22.5L, -11.048303064729054768L, 36},
            {10, 18000, 22.5L, -11.781841072048607089L, 40},
            {11, 22000, 22.5L, -12.521130149230231154L, 44},
            {5, 11000, 25.0L, -7.781335792605759138L, 20},
            {6, 13000, 25.0L, -8.861575115516878043L, 24},
            {7, 13000, 25.0L, -9.869101102887270827L, 28},
            {8, 15000, 25.0L, -10.836019353192511971L, 32},
            {9, 18000, 25.0L, -11.719296695714589873L, 36},
            {10, 18000, 25.0L, -12.527389550921412820L, 40},
            {11, 22000, 25.0L, -13.310366269709135167L, 44},
            {5, 11000, 30.0L, -8.490976596434499868L, 20},
            {6, 13000, 30.0L, -9.696013978803513567L, 24},
            {7, 13000, 30.0L, -10.842613110055083325L, 28},
            {8, 15000, 30.0L, -11.927091472703111208L, 32},
            {9, 18000, 30.0L, -12.912948921413303083L, 36},
            {10, 18000, 30.0L, -13.853669801697762409L, 40},
            {11, 22000, 30.0L, -14.719154343835811005L, 44},
            {5, 11000, 35.0L, -9.104254417906034291L, 20},
            {6, 13000, 35.0L, -10.418966523378784927L, 24},
            {7, 13000, 35.0L, -11.686057745392908203L, 28},
            {8, 15000, 35.0L, -12.874606142392963193L, 32},
            {9, 18000, 35.0L, -13.952881616312993174L, 36},
            {10, 18000, 35.0L, -15.009150573808526019L, 40},
            {11, 22000, 35.0L, -15.951192249754356567L, 44},
        },
    };
}

ScanConfig c19_stable_config() {
    // Diagnostic C_19 direct-tail scan using only stable moments below the
    // first C_19 boundary Delta_40.  A positive candidate at target_n=63
    // would close C_19 without new exact-delta jobs, after exact replay.
    return {
        "C19-post-m29-direct-tail-stable",
        19,
        19 * (2 * 19 + 1),
        20,
        19,
        63,
        140.0L,
        410.0L,
        0.05L,
        0.05L,
        8.0L,
        0.01L,
        0.05L,
        16.0L,
        0.01L,
        {
            {5, 11000, 20.5L, -7.004532338304401L, 20},
            {6, 13000, 20.5L, -7.971050748417397L, 24},
            {7, 15000, 20.5L, -8.860416764339966L, 28},
            {8, 17000, 20.5L, -9.685764456356026L, 32},
            {9, 19000, 20.5L, -10.456995634037611L, 36},
            {5, 11000, 22.5L, -7.358675856993777L, 20},
            {6, 13000, 22.5L, -8.386533246257557L, 24},
            {7, 15000, 22.5L, -9.334785274509734L, 28},
            {8, 17000, 22.5L, -10.216787147820710L, 32},
            {9, 19000, 22.5L, -11.042633105062820L, 36},
            {5, 11000, 25.0L, -7.765053333809909L, 20},
            {6, 13000, 25.0L, -8.864069471181729L, 24},
            {7, 15000, 25.0L, -9.880813700910267L, 28},
            {8, 17000, 25.0L, -10.828866344127490L, 32},
            {9, 19000, 25.0L, -11.718519341820894L, 36},
            {5, 11000, 27.5L, -8.137537413444818L, 20},
            {6, 13000, 27.5L, -9.302460708637888L, 24},
            {7, 15000, 27.5L, -10.382805496933294L, 28},
            {8, 17000, 27.5L, -11.392333839317587L, 32},
            {9, 19000, 27.5L, -12.341498729441724L, 36},
            {5, 11000, 30.0L, -8.481467855376750L, 20},
            {6, 13000, 30.0L, -9.707791933622644L, 24},
            {7, 15000, 30.0L, -10.847522485298512L, 28},
            {8, 17000, 30.0L, -11.914570370991015L, 32},
            {9, 19000, 30.0L, -12.919521487268980L, 36},
            {5, 11000, 35.0L, -9.099465085193827L, 20},
            {6, 13000, 35.0L, -10.437340005073587L, 24},
            {7, 15000, 35.0L, -11.685269610415418L, 28},
            {8, 17000, 35.0L, -12.857391450327569L, 32},
            {9, 19000, 35.0L, -13.964495872004647L, 36},
        },
    };
}

ScanConfig c19_exact40_config() {
    // Diagnostic C_19 direct-tail scan using the exact Delta_40 boundary
    // correction.  The degree-10 majorants are exact-replay outputs with
    // that correction included.  A positive candidate must still be replayed
    // exactly in classical_tail_constants.py.
    ScanConfig cfg = c19_stable_config();
    cfg.label = "C19-post-m29-direct-tail-exact40";
    cfg.cheb_bounds.push_back({10, 21000, 20.5L, -11.181799601893502L, 40});
    cfg.cheb_bounds.push_back({10, 21000, 22.5L, -11.820180902726378L, 40});
    cfg.cheb_bounds.push_back({10, 21000, 25.0L, -12.557804984810929L, 40});
    cfg.cheb_bounds.push_back({10, 21000, 27.5L, -13.238475823987372L, 40});
    cfg.cheb_bounds.push_back({10, 21000, 30.0L, -13.870670664659343L, 40});
    cfg.cheb_bounds.push_back({10, 21000, 35.0L, -15.015062711217738L, 40});
    return cfg;
}

ScanConfig d12_exact51_config() {
    // Diagnostic D_12 direct-tail scan using the exact local determinant
    // corrections through Delta_51.  The base adjustment replaces the B/C
    // Weyl order and degree product by the type-D values.
    ScanConfig cfg{
        "D12-post-m29-direct-tail-exact51",
        12,
        12 * (2 * 12 - 1),
        22,
        12,
        79,
        100.0L,
        220.0L,
        0.05L,
        0.05L,
        8.0L,
        0.01L,
        0.05L,
        16.0L,
        0.01L,
        {
            {10, 7000, 16.5L, -9.869019015118463L, 40},
            {10, 8000, 16.5L, -9.961439339215385L, 40},
            {10, 9000, 16.5L, -9.709417669649355L, 40},
            {11, 8000, 16.5L, -10.573586382824606L, 44},
            {11, 9000, 16.5L, -10.484440035332206L, 44},
            {12, 8000, 16.5L, -10.863478513158098L, 48},
            {12, 9000, 16.5L, -11.165189067298880L, 48},
            {12, 10000, 16.5L, -10.957275729346270L, 48},
            {13, 9000, 16.5L, -11.592298244788275L, 52},
            {13, 10000, 16.5L, -11.660333751791756L, 52},
            {14, 9000, 16.5L, -11.695574182705789L, 56},
            {14, 10000, 16.5L, -12.214574062261708L, 56},
            {15, 10000, 16.5L, -12.447453255427288L, 60},
            {15, 12000, 16.5L, -12.488343736796258L, 60},
            {16, 12000, 16.5L, -13.131803224493600L, 64},
            {16, 14000, 16.5L, -12.582595137927873L, 64},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d14_exact49_config() {
    // Diagnostic D_14 direct-tail scan using the exact local determinant
    // corrections through Delta_49.  A positive candidate at target_n=325
    // would combine with the stable-envelope Chain bridge through n=323.
    ScanConfig cfg{
        "D14-post-m29-direct-tail-exact49",
        14,
        14 * (2 * 14 - 1),
        26,
        14,
        325,
        120.0L,
        340.0L,
        0.05L,
        0.05L,
        8.0L,
        0.01L,
        0.05L,
        16.0L,
        0.01L,
        {
            {10, 12000, 60.0L, -19.274977883217459862L, 40},
            {11, 12000, 60.0L, -20.963062411247534556L, 44},
            {12, 12000, 60.0L, -22.442664146282567117L, 48},
            {13, 14000, 60.0L, -23.618862205152424849L, 52},
            {14, 14000, 60.0L, -25.063439010089737735L, 56},
            {15, 14000, 60.0L, -26.159563280098737437L, 60},
            {16, 16000, 60.0L, -27.523953388769683670L, 64},
            {17, 16000, 60.0L, -28.611520818063127081L, 68},
            {18, 18000, 60.0L, -29.848114472145880427L, 72},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d15_exact48_config() {
    // Diagnostic D_15 direct-tail scan using the accepted determinant
    // corrections through Delta_48.  A positive candidate at target_n=525
    // would combine with the stable-envelope Chain bridge through n=523.
    ScanConfig cfg{
        "D15-post-m29-direct-tail-exact48",
        15,
        15 * (2 * 15 - 1),
        28,
        15,
        525,
        130.0L,
        390.0L,
        0.05L,
        0.05L,
        8.0L,
        0.01L,
        0.05L,
        16.0L,
        0.01L,
        {
            {8, 10000, 70.0L, -17.597134115735272530L, 32},
            {9, 12000, 70.0L, -19.101883109140032957L, 36},
            {10, 12000, 70.0L, -20.847429661724291350L, 40},
            {11, 14000, 70.0L, -22.247784273530285759L, 44},
            {12, 14000, 70.0L, -23.860896183447891872L, 48},
            {13, 16000, 70.0L, -25.233798871472558289L, 52},
            {8, 10000, 80.0L, -18.534413718800749393L, 32},
            {9, 12000, 80.0L, -20.143491308338326462L, 36},
            {10, 12000, 80.0L, -22.004772105277936589L, 40},
            {11, 14000, 80.0L, -23.505737812587454982L, 44},
            {12, 14000, 80.0L, -25.233209135146623225L, 48},
            {13, 16000, 80.0L, -26.703221324897555178L, 52},
            {8, 10000, 90.0L, -19.375216850754142683L, 32},
            {9, 12000, 90.0L, -21.079168293520496036L, 36},
            {10, 12000, 90.0L, -23.044413199924804303L, 40},
            {11, 14000, 90.0L, -24.637241139323108996L, 44},
            {12, 14000, 90.0L, -26.467576400676435355L, 48},
            {13, 16000, 90.0L, -28.026595384857500903L, 52},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d16_exact45_config() {
    // Diagnostic D_16 direct-tail scan using exact determinant corrections
    // through Delta_45.  A positive candidate would combine with the
    // stable-envelope Chain bridge up to the preceding odd exponent.
    ScanConfig cfg{
        "D16-post-m29-direct-tail-exact45",
        16,
        16 * (2 * 16 - 1),
        30,
        16,
        725,
        120.0L,
        440.0L,
        0.05L,
        0.05L,
        8.0L,
        0.01L,
        0.05L,
        18.0L,
        0.01L,
        {
            {8, 10000, 80.0L, -18.354907296276095963L, 32},
            {8, 12000, 80.0L, -18.392687792079101428L, 32},
            {8, 14000, 80.0L, -17.965999812146435488L, 32},
            {8, 16000, 80.0L, -17.490275879779076718L, 32},
            {9, 10000, 80.0L, -19.811380222591239431L, 36},
            {9, 12000, 80.0L, -20.200386888432234400L, 36},
            {9, 14000, 80.0L, -19.995794354493057199L, 36},
            {9, 16000, 80.0L, -19.504874787795159818L, 36},
            {10, 10000, 80.0L, -21.130469631262300823L, 40},
            {10, 12000, 80.0L, -21.674847414761387654L, 40},
            {10, 14000, 80.0L, -21.881783053202653377L, 40},
            {10, 16000, 80.0L, -21.486988097329856373L, 40},
            {8, 10000, 90.0L, -19.193168672695051669L, 32},
            {8, 12000, 90.0L, -19.222798157132046981L, 32},
            {8, 14000, 90.0L, -18.788192645904970846L, 32},
            {8, 16000, 90.0L, -18.304773713828936366L, 32},
            {9, 10000, 90.0L, -20.754424271062561047L, 36},
            {9, 12000, 90.0L, -21.134261049116801701L, 36},
            {9, 14000, 90.0L, -20.920761292471397041L, 36},
            {9, 16000, 90.0L, -20.421184851101244817L, 36},
            {10, 10000, 90.0L, -22.178296351786002560L, 40},
            {10, 12000, 90.0L, -22.712485371077576701L, 40},
            {10, 14000, 90.0L, -22.909524095400811916L, 40},
            {10, 16000, 90.0L, -22.505110389892180933L, 40},
            {8, 10000, 100.0L, -19.953198959218241271L, 32},
            {8, 12000, 100.0L, -19.976195173053994836L, 32},
            {8, 14000, 100.0L, -19.535127238054286636L, 32},
            {8, 16000, 100.0L, -19.045409512904299731L, 32},
            {9, 10000, 100.0L, -21.609458343401144020L, 36},
            {9, 12000, 100.0L, -21.981832692028987708L, 36},
            {9, 14000, 100.0L, -21.761062708639386187L, 36},
            {9, 16000, 100.0L, -21.254400125061039262L, 36},
            {10, 10000, 100.0L, -23.128334209939993116L, 40},
            {10, 12000, 100.0L, -23.654231640980000861L, 40},
            {10, 14000, 100.0L, -23.843192335587460207L, 40},
            {10, 16000, 100.0L, -23.430905138736406457L, 40},
            {8, 10000, 110.0L, -20.649196948021696585L, 32},
            {8, 12000, 110.0L, -20.666691747943971791L, 32},
            {8, 14000, 110.0L, -20.220250918913009741L, 32},
            {8, 16000, 110.0L, -19.725283888856921521L, 32},
            {9, 10000, 110.0L, -22.392456080805033025L, 36},
            {9, 12000, 110.0L, -22.758641338780208230L, 36},
            {9, 14000, 110.0L, -22.531826849605437246L, 36},
            {9, 16000, 110.0L, -22.019258798007740552L, 36},
            {10, 10000, 110.0L, -23.998331695944301600L, 40},
            {10, 12000, 110.0L, -24.517352359592479161L, 40},
            {10, 14000, 110.0L, -24.699596936660867641L, 40},
            {10, 16000, 110.0L, -24.280748108677173036L, 40},
            {8, 10000, 120.0L, -21.291792087704067171L, 32},
            {8, 12000, 120.0L, -21.304652040604196372L, 32},
            {8, 14000, 120.0L, -20.853675290653725938L, 32},
            {8, 16000, 120.0L, -20.354267793879472492L, 32},
            {9, 10000, 120.0L, -23.115375612947701711L, 36},
            {9, 12000, 120.0L, -23.476346668022969766L, 36},
            {9, 14000, 120.0L, -23.244429267813757178L, 36},
            {9, 16000, 120.0L, -22.726865691158110394L, 36},
            {10, 10000, 120.0L, -24.801575620547282597L, 40},
            {10, 12000, 120.0L, -25.314802725417763440L, 40},
            {10, 14000, 120.0L, -25.491377401336777098L, 40},
            {10, 16000, 120.0L, -25.066977989955361750L, 40},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d17_exact44_config() {
    // Diagnostic D_17 direct-tail scan using exact determinant corrections
    // through Delta_44.  A positive candidate combines with the
    // stable-envelope Chain bridge up to the preceding odd exponent.
    ScanConfig cfg{
        "D17-post-m29-direct-tail-exact44",
        17,
        17 * (2 * 17 - 1),
        32,
        17,
        287,
        200.0L,
        270.0L,
        0.05L,
        0.50L,
        1.50L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -18.362513018228378002L, 32},
            {9, 14000, 80.0L, -20.131534254571019215L, 36},
            {10, 16000, 80.0L, -21.784345118255504923L, 40},
            {11, 16000, 80.0L, -23.357809375390360583L, 44},
            {8, 12000, 90.0L, -19.191289419468069915L, 32},
            {9, 14000, 90.0L, -21.055860944708608494L, 36},
            {10, 16000, 90.0L, -22.802661719641180544L, 40},
            {11, 16000, 90.0L, -24.477957636914603766L, 44},
            {8, 12000, 100.0L, -19.942418951120048173L, 32},
            {9, 14000, 100.0L, -21.894326476541669249L, 36},
            {10, 16000, 100.0L, -23.727173598727322851L, 40},
            {11, 16000, 100.0L, -25.494920703909360303L, 44},
            {8, 12000, 110.0L, -20.629979267416901696L, 32},
            {9, 14000, 110.0L, -22.662390146071237454L, 36},
            {10, 16000, 110.0L, -24.574657287025879400L, 40},
            {11, 16000, 110.0L, -26.427152761037772508L, 44},
            {8, 12000, 120.0L, -21.264515799079283987L, 32},
            {9, 14000, 120.0L, -23.371655336346907891L, 36},
            {10, 16000, 120.0L, -25.357728611159734386L, 40},
            {11, 16000, 120.0L, -27.288531217585012992L, 44},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d18_exact43_config() {
    // Diagnostic D_18 direct-tail scan using exact determinant corrections
    // through Delta_43.  A positive candidate combines with the
    // stable-envelope Chain bridge up to the preceding odd exponent.
    ScanConfig cfg{
        "D18-post-m29-direct-tail-exact43",
        18,
        18 * (2 * 18 - 1),
        34,
        18,
        365,
        290.0L,
        330.0L,
        0.05L,
        0.50L,
        1.50L,
        0.01L,
        0.20L,
        5.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -18.094572166203064761L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d19_exact44_config() {
    // Diagnostic D_19 direct-tail scan using exact determinant corrections
    // through Delta_44.  A positive candidate combines with the
    // stable-envelope Chain bridge up to the preceding odd exponent.
    ScanConfig cfg{
        "D19-post-m29-direct-tail-exact44",
        19,
        19 * (2 * 19 - 1),
        36,
        19,
        405,
        320.0L,
        380.0L,
        0.05L,
        0.50L,
        1.50L,
        0.01L,
        0.20L,
        5.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -17.794644479668351023L, 32},
            {9, 14000, 80.0L, -19.537416490104534041L, 36},
            {10, 16000, 80.0L, -21.214837266556195562L, 40},
            {8, 12000, 90.0L, -18.621349670596307461L, 32},
            {8, 12000, 100.0L, -19.368888451333196593L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d28_exact35_config() {
    // Diagnostic D_28 direct-tail scan using exact determinant corrections
    // through Delta_35.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D28-post-m29-direct-tail-exact35",
        28,
        28 * (2 * 28 - 1),
        54,
        28,
        421,
        760.0L,
        1040.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -16.328976478657791756L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d29_exact35_config() {
    // Diagnostic D_29 direct-tail scan using exact determinant corrections
    // through Delta_35.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D29-post-m29-direct-tail-exact35",
        29,
        29 * (2 * 29 - 1),
        56,
        29,
        433,
        820.0L,
        1080.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -16.247072859292647640L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d30_exact35_config() {
    // Diagnostic D_30 direct-tail scan using exact determinant corrections
    // through Delta_35.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D30-post-m29-direct-tail-exact35",
        30,
        30 * (2 * 30 - 1),
        58,
        30,
        449,
        880.0L,
        1120.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -16.173526254152999250L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d31_exact35_config() {
    // Diagnostic D_31 direct-tail scan using exact determinant corrections
    // through Delta_35.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D31-post-m29-direct-tail-exact35",
        31,
        31 * (2 * 31 - 1),
        60,
        29,
        481,
        920.0L,
        1220.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -16.107229314344294043L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d32_exact35_config() {
    // Diagnostic D_32 direct-tail scan using exact determinant corrections
    // through Delta_35.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D32-post-m29-direct-tail-exact35",
        32,
        32 * (2 * 32 - 1),
        62,
        32,
        497,
        980.0L,
        1205.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -16.047253398809909007L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d33_exact35_config() {
    // Diagnostic D_33 direct-tail scan using exact determinant corrections
    // through Delta_35.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D33-post-m29-direct-tail-exact35",
        33,
        33 * (2 * 33 - 1),
        64,
        31,
        513,
        1020.0L,
        1240.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.992814855710117387L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d34_exact35_config() {
    // Diagnostic D_34 direct-tail scan using exact determinant corrections
    // through Delta_35.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D34-post-m29-direct-tail-exact35",
        34,
        34 * (2 * 34 - 1),
        66,
        34,
        545,
        1060.0L,
        1300.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.943248529192711727L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d35_exact35_config() {
    // Diagnostic D_35 direct-tail scan using exact determinant correction
    // through Delta_35.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D35-post-m29-direct-tail-exact35",
        35,
        35 * (2 * 35 - 1),
        68,
        33,
        577,
        1100.0L,
        1320.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.897986772779222504L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d36_stable35_config() {
    // Diagnostic D_36 direct-tail scan using stable determinant envelopes
    // from the first D_36 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D36-post-m29-direct-tail-stable35",
        36,
        36 * (2 * 36 - 1),
        70,
        36,
        609,
        1140.0L,
        1360.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.856542695385082799L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d37_stable35_config() {
    // Diagnostic D_37 direct-tail scan using stable determinant envelopes
    // from the first D_37 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D37-post-m29-direct-tail-stable35",
        37,
        37 * (2 * 37 - 1),
        72,
        35,
        625,
        1180.0L,
        1500.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.818496688171853488L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d38_stable35_config() {
    // Diagnostic D_38 direct-tail scan using stable determinant envelopes
    // from the first D_38 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D38-post-m29-direct-tail-stable35",
        38,
        38 * (2 * 38 - 1),
        74,
        38,
        657,
        1220.0L,
        1600.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.783485515049097713L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d39_stable35_config() {
    // Diagnostic D_39 direct-tail scan using stable determinant envelopes
    // from the first D_39 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D39-post-m29-direct-tail-stable35",
        39,
        39 * (2 * 39 - 1),
        76,
        37,
        673,
        1260.0L,
        1660.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.751193421807068945L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d40_stable35_config() {
    // Diagnostic D_40 direct-tail scan using stable determinant envelopes
    // from the first D_40 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D40-post-m29-direct-tail-stable35",
        40,
        40 * (2 * 40 - 1),
        78,
        40,
        705,
        1300.0L,
        1720.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.721344846262262362L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d41_stable35_config() {
    // Diagnostic D_41 direct-tail scan using stable determinant envelopes
    // from the first D_41 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D41-post-m29-direct-tail-stable35",
        41,
        41 * (2 * 41 - 1),
        80,
        39,
        737,
        1340.0L,
        1800.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.693698406882747677L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d42_stable35_config() {
    // Diagnostic D_42 direct-tail scan using stable determinant envelopes
    // from the first D_42 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D42-post-m29-direct-tail-stable35",
        42,
        42 * (2 * 42 - 1),
        82,
        42,
        769,
        1380.0L,
        1880.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.668041918921885536L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d43_stable35_config() {
    // Diagnostic D_43 direct-tail scan using stable determinant envelopes
    // from the first D_43 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D43-post-m29-direct-tail-stable35",
        43,
        43 * (2 * 43 - 1),
        84,
        41,
        801,
        1420.0L,
        1960.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.644188241382080461L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d44_stable35_config() {
    // Diagnostic D_44 direct-tail scan using stable determinant envelopes
    // from the first D_44 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D44-post-m29-direct-tail-stable35",
        44,
        44 * (2 * 44 - 1),
        86,
        44,
        833,
        1460.0L,
        2040.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.621971799639412574L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d45_stable35_config() {
    // Diagnostic D_45 direct-tail scan using stable determinant envelopes
    // from the first D_45 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D45-post-m29-direct-tail-stable35",
        45,
        45 * (2 * 45 - 1),
        88,
        43,
        881,
        1500.0L,
        2120.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.601245660528103301L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d46_stable35_config() {
    // Diagnostic D_46 direct-tail scan using stable determinant envelopes
    // from the first D_46 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D46-post-m29-direct-tail-stable35",
        46,
        46 * (2 * 46 - 1),
        90,
        46,
        913,
        1540.0L,
        2200.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.581879061475227242L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d47_stable35_config() {
    // Diagnostic D_47 direct-tail scan using stable determinant envelopes
    // from the first D_47 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D47-post-m29-direct-tail-stable35",
        47,
        47 * (2 * 47 - 1),
        92,
        45,
        977,
        1580.0L,
        2280.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.563755314628980167L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d48_stable35_config() {
    // Diagnostic D_48 direct-tail scan using stable determinant envelopes
    // from the first D_48 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D48-post-m29-direct-tail-stable35",
        48,
        48 * (2 * 48 - 1),
        94,
        48,
        1009,
        1620.0L,
        2360.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d49_stable35_config() {
    // Diagnostic D_49 direct-tail scan using stable determinant envelopes
    // from the first D_49 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D49-post-m29-direct-tail-stable35",
        49,
        49 * (2 * 49 - 1),
        96,
        47,
        1073,
        1660.0L,
        2440.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d50_stable35_config() {
    // Diagnostic D_50 direct-tail scan using stable determinant envelopes
    // from the first D_50 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D50-post-m29-direct-tail-stable35",
        50,
        50 * (2 * 50 - 1),
        98,
        50,
        1093,
        1700.0L,
        2520.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d51_stable35_config() {
    // Diagnostic D_51 direct-tail scan using stable determinant envelopes
    // from the first D_51 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D51-post-m29-direct-tail-stable35",
        51,
        51 * (2 * 51 - 1),
        100,
        49,
        1155,
        1740.0L,
        2600.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d52_stable35_config() {
    // Diagnostic D_52 direct-tail scan using stable determinant envelopes
    // from the first D_52 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D52-post-m29-direct-tail-stable35",
        52,
        52 * (2 * 52 - 1),
        102,
        52,
        1191,
        1780.0L,
        2680.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d53_stable35_config() {
    // Diagnostic D_53 direct-tail scan using stable determinant envelopes
    // from the first D_53 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D53-post-m29-direct-tail-stable35",
        53,
        53 * (2 * 53 - 1),
        104,
        51,
        1255,
        1820.0L,
        2760.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d54_stable35_config() {
    // Diagnostic D_54 direct-tail scan using stable determinant envelopes
    // from the first D_54 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D54-post-m29-direct-tail-stable35",
        54,
        54 * (2 * 54 - 1),
        106,
        54,
        1293,
        1860.0L,
        2840.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d55_stable35_config() {
    // Diagnostic D_55 direct-tail scan using stable determinant envelopes
    // from the first D_55 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D55-post-m29-direct-tail-stable35",
        55,
        55 * (2 * 55 - 1),
        108,
        53,
        1359,
        1900.0L,
        2920.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d56_stable35_config() {
    // Diagnostic D_56 direct-tail scan using stable determinant envelopes
    // from the first D_56 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D56-post-m29-direct-tail-stable35",
        56,
        56 * (2 * 56 - 1),
        110,
        56,
        1399,
        1940.0L,
        3000.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d57_stable35_config() {
    // Diagnostic D_57 direct-tail scan using stable determinant envelopes
    // from the first D_57 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D57-post-m29-direct-tail-stable35",
        57,
        57 * (2 * 57 - 1),
        112,
        55,
        1469,
        1980.0L,
        3080.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d58_stable35_config() {
    // Diagnostic D_58 direct-tail scan using stable determinant envelopes
    // from the first D_58 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D58-post-m29-direct-tail-stable35",
        58,
        58 * (2 * 58 - 1),
        114,
        58,
        1511,
        2020.0L,
        3160.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d59_stable35_config() {
    // Diagnostic D_59 direct-tail scan using stable determinant envelopes
    // from the first D_59 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D59-post-m29-direct-tail-stable35",
        59,
        59 * (2 * 59 - 1),
        116,
        57,
        1585,
        2060.0L,
        3240.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d60_stable35_config() {
    // Diagnostic D_60 direct-tail scan using stable determinant envelopes
    // from the first D_60 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D60-post-m29-direct-tail-stable35",
        60,
        60 * (2 * 60 - 1),
        118,
        60,
        1627,
        2100.0L,
        3320.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d61_stable35_config() {
    // Diagnostic D_61 direct-tail scan using stable determinant envelopes
    // from the first D_61 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D61-post-m29-direct-tail-stable35",
        61,
        61 * (2 * 61 - 1),
        120,
        59,
        1703,
        2140.0L,
        3400.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d62_stable35_config() {
    // Diagnostic D_62 direct-tail scan using stable determinant envelopes
    // from the first D_62 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D62-post-m29-direct-tail-stable35",
        62,
        62 * (2 * 62 - 1),
        122,
        62,
        1753,
        2180.0L,
        3480.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d63_stable35_config() {
    // Diagnostic D_63 direct-tail scan using stable determinant envelopes
    // from the first D_63 boundary.  The common sine-product range requires
    // 2*radius/(39*kappa)<1, enforced by the scanner.
    ScanConfig cfg{
        "D63-post-m29-direct-tail-stable35",
        63,
        63 * (2 * 63 - 1),
        124,
        61,
        1833,
        2220.0L,
        3560.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

ScanConfig d_stable35_generated_config(int rank) {
    if (rank < 48) {
        std::cerr << "--d-stable35-rank is calibrated for D-ranks >= 48" << std::endl;
        std::exit(2);
    }
    const int dimension = rank * (2 * rank - 1);
    const int kappa = 2 * rank - 2;
    const int c_value = (rank % 2 == 0) ? rank : rank - 2;
    ScanConfig cfg{
        "D" + std::to_string(rank) + "-post-m29-direct-tail-stable35-generated",
        rank,
        dimension,
        kappa,
        c_value,
        next_odd_at_least(30 * rank),
        40.0L * rank - 300.0L,
        80.0L * rank - 1480.0L,
        0.05L,
        0.50L,
        2.00L,
        0.01L,
        0.20L,
        6.00L,
        0.01L,
        {
            {8, 12000, 80.0L, -15.0L, 32},
        },
    };
    cfg.base_adjust_log10 =
        LOG10_2 + log10_factorial(cfg.rank) - log10_factorial(2 * cfg.rank);
    return cfg;
}

}  // namespace

int main(int argc, char** argv) {
    bool progress = false;
    bool binary_search = false;
    std::string case_name = "c17";
    std::uint64_t samples = 5000000ULL;
    int target_n_override = 0;
    int target_lo_override = 0;
    int target_hi_override = 0;
    int d_stable35_rank = 0;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--progress") {
            progress = true;
        } else if (arg == "--binary-search") {
            binary_search = true;
        } else if (arg == "--c17") {
            case_name = "c17";
        } else if (arg == "--b13-exact54") {
            case_name = "b13-exact54";
        } else if (arg == "--b13-exact57") {
            case_name = "b13-exact57";
        } else if (arg == "--b14-exact39") {
            case_name = "b14-exact39";
        } else if (arg == "--c18-exact39") {
            case_name = "c18-exact39";
        } else if (arg == "--c18-exact47") {
            case_name = "c18-exact47";
        } else if (arg == "--c19-stable") {
            case_name = "c19-stable";
        } else if (arg == "--c19-exact40") {
            case_name = "c19-exact40";
        } else if (arg == "--d12-exact51") {
            case_name = "d12-exact51";
        } else if (arg == "--d14-exact49") {
            case_name = "d14-exact49";
        } else if (arg == "--d15-exact48") {
            case_name = "d15-exact48";
        } else if (arg == "--d16-exact45") {
            case_name = "d16-exact45";
        } else if (arg == "--d17-exact44") {
            case_name = "d17-exact44";
        } else if (arg == "--d18-exact43") {
            case_name = "d18-exact43";
        } else if (arg == "--d19-exact44") {
            case_name = "d19-exact44";
        } else if (arg == "--d28-exact35") {
            case_name = "d28-exact35";
        } else if (arg == "--d29-exact35") {
            case_name = "d29-exact35";
        } else if (arg == "--d30-exact35") {
            case_name = "d30-exact35";
        } else if (arg == "--d31-exact35") {
            case_name = "d31-exact35";
        } else if (arg == "--d32-exact35") {
            case_name = "d32-exact35";
        } else if (arg == "--d33-exact35") {
            case_name = "d33-exact35";
        } else if (arg == "--d34-exact35") {
            case_name = "d34-exact35";
        } else if (arg == "--d35-exact35") {
            case_name = "d35-exact35";
        } else if (arg == "--d36-stable35") {
            case_name = "d36-stable35";
        } else if (arg == "--d37-stable35") {
            case_name = "d37-stable35";
        } else if (arg == "--d38-stable35") {
            case_name = "d38-stable35";
        } else if (arg == "--d39-stable35") {
            case_name = "d39-stable35";
        } else if (arg == "--d40-stable35") {
            case_name = "d40-stable35";
        } else if (arg == "--d41-stable35") {
            case_name = "d41-stable35";
        } else if (arg == "--d42-stable35") {
            case_name = "d42-stable35";
        } else if (arg == "--d43-stable35") {
            case_name = "d43-stable35";
        } else if (arg == "--d44-stable35") {
            case_name = "d44-stable35";
        } else if (arg == "--d45-stable35") {
            case_name = "d45-stable35";
        } else if (arg == "--d46-stable35") {
            case_name = "d46-stable35";
        } else if (arg == "--d47-stable35") {
            case_name = "d47-stable35";
        } else if (arg == "--d48-stable35") {
            case_name = "d48-stable35";
        } else if (arg == "--d49-stable35") {
            case_name = "d49-stable35";
        } else if (arg == "--d50-stable35") {
            case_name = "d50-stable35";
        } else if (arg == "--d51-stable35") {
            case_name = "d51-stable35";
        } else if (arg == "--d52-stable35") {
            case_name = "d52-stable35";
        } else if (arg == "--d53-stable35") {
            case_name = "d53-stable35";
        } else if (arg == "--d54-stable35") {
            case_name = "d54-stable35";
        } else if (arg == "--d55-stable35") {
            case_name = "d55-stable35";
        } else if (arg == "--d56-stable35") {
            case_name = "d56-stable35";
        } else if (arg == "--d57-stable35") {
            case_name = "d57-stable35";
        } else if (arg == "--d58-stable35") {
            case_name = "d58-stable35";
        } else if (arg == "--d59-stable35") {
            case_name = "d59-stable35";
        } else if (arg == "--d60-stable35") {
            case_name = "d60-stable35";
        } else if (arg == "--d61-stable35") {
            case_name = "d61-stable35";
        } else if (arg == "--d62-stable35") {
            case_name = "d62-stable35";
        } else if (arg == "--d63-stable35") {
            case_name = "d63-stable35";
        } else if (arg == "--d-stable35-rank") {
            if (i + 1 >= argc) {
                std::cerr << "missing --d-stable35-rank value" << std::endl;
                return 2;
            }
            case_name = "d-stable35-rank";
            d_stable35_rank = std::atoi(argv[++i]);
        } else if (arg == "--samples") {
            if (i + 1 >= argc) {
                std::cerr << "missing --samples value" << std::endl;
                return 2;
            }
            samples = std::strtoull(argv[++i], nullptr, 10);
        } else if (arg == "--target-n") {
            if (i + 1 >= argc) {
                std::cerr << "missing --target-n value" << std::endl;
                return 2;
            }
            target_n_override = std::atoi(argv[++i]);
        } else if (arg == "--target-lo") {
            if (i + 1 >= argc) {
                std::cerr << "missing --target-lo value" << std::endl;
                return 2;
            }
            target_lo_override = std::atoi(argv[++i]);
        } else if (arg == "--target-hi") {
            if (i + 1 >= argc) {
                std::cerr << "missing --target-hi value" << std::endl;
                return 2;
            }
            target_hi_override = std::atoi(argv[++i]);
        } else {
            std::cerr << "unknown argument: " << arg << std::endl;
            return 2;
        }
    }

    ScanConfig cfg = c17_config();
    if (case_name == "b13-exact54") cfg = b13_exact54_config();
    if (case_name == "b13-exact57") cfg = b13_exact57_config();
    if (case_name == "b14-exact39") cfg = b14_exact39_config();
    if (case_name == "c18-exact39") cfg = c18_exact39_config();
    if (case_name == "c18-exact47") cfg = c18_exact47_config();
    if (case_name == "c19-stable") cfg = c19_stable_config();
    if (case_name == "c19-exact40") cfg = c19_exact40_config();
    if (case_name == "d12-exact51") cfg = d12_exact51_config();
    if (case_name == "d14-exact49") cfg = d14_exact49_config();
    if (case_name == "d15-exact48") cfg = d15_exact48_config();
    if (case_name == "d16-exact45") cfg = d16_exact45_config();
    if (case_name == "d17-exact44") cfg = d17_exact44_config();
    if (case_name == "d18-exact43") cfg = d18_exact43_config();
    if (case_name == "d19-exact44") cfg = d19_exact44_config();
    if (case_name == "d28-exact35") cfg = d28_exact35_config();
    if (case_name == "d29-exact35") cfg = d29_exact35_config();
    if (case_name == "d30-exact35") cfg = d30_exact35_config();
    if (case_name == "d31-exact35") cfg = d31_exact35_config();
    if (case_name == "d32-exact35") cfg = d32_exact35_config();
    if (case_name == "d33-exact35") cfg = d33_exact35_config();
    if (case_name == "d34-exact35") cfg = d34_exact35_config();
    if (case_name == "d35-exact35") cfg = d35_exact35_config();
    if (case_name == "d36-stable35") cfg = d36_stable35_config();
    if (case_name == "d37-stable35") cfg = d37_stable35_config();
    if (case_name == "d38-stable35") cfg = d38_stable35_config();
    if (case_name == "d39-stable35") cfg = d39_stable35_config();
    if (case_name == "d40-stable35") cfg = d40_stable35_config();
    if (case_name == "d41-stable35") cfg = d41_stable35_config();
    if (case_name == "d42-stable35") cfg = d42_stable35_config();
    if (case_name == "d43-stable35") cfg = d43_stable35_config();
    if (case_name == "d44-stable35") cfg = d44_stable35_config();
    if (case_name == "d45-stable35") cfg = d45_stable35_config();
    if (case_name == "d46-stable35") cfg = d46_stable35_config();
    if (case_name == "d47-stable35") cfg = d47_stable35_config();
    if (case_name == "d48-stable35") cfg = d48_stable35_config();
    if (case_name == "d49-stable35") cfg = d49_stable35_config();
    if (case_name == "d50-stable35") cfg = d50_stable35_config();
    if (case_name == "d51-stable35") cfg = d51_stable35_config();
    if (case_name == "d52-stable35") cfg = d52_stable35_config();
    if (case_name == "d53-stable35") cfg = d53_stable35_config();
    if (case_name == "d54-stable35") cfg = d54_stable35_config();
    if (case_name == "d55-stable35") cfg = d55_stable35_config();
    if (case_name == "d56-stable35") cfg = d56_stable35_config();
    if (case_name == "d57-stable35") cfg = d57_stable35_config();
    if (case_name == "d58-stable35") cfg = d58_stable35_config();
    if (case_name == "d59-stable35") cfg = d59_stable35_config();
    if (case_name == "d60-stable35") cfg = d60_stable35_config();
    if (case_name == "d61-stable35") cfg = d61_stable35_config();
    if (case_name == "d62-stable35") cfg = d62_stable35_config();
    if (case_name == "d63-stable35") cfg = d63_stable35_config();
    if (case_name == "d-stable35-rank") cfg = d_stable35_generated_config(d_stable35_rank);
    if (binary_search) {
        const int lo =
            target_lo_override > 0
                ? target_lo_override
                : (target_n_override > 0 ? target_n_override : cfg.target_n);
        const int hi = target_hi_override > 0 ? target_hi_override : cfg.target_n;
        return binary_search_onset(cfg, lo, hi, samples, progress);
    }
    if (target_n_override > 0) cfg.target_n = target_n_override;
    const Candidate best = sample_config(cfg, samples, progress);
    print_best(cfg.label, best);
    return best.margin > 0.0L ? 0 : 1;
}
