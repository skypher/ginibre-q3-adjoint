#if __has_include(<mpfr.h>)
#include <mpfr.h>
#else
#include <gmp.h>
extern "C" {
typedef long mpfr_prec_t;
typedef long mpfr_exp_t;
typedef int mpfr_sign_t;
typedef enum {
    MPFR_RNDN = 0,
    MPFR_RNDZ = 1,
    MPFR_RNDU = 2,
    MPFR_RNDD = 3,
    MPFR_RNDA = 4,
    MPFR_RNDF = 5,
    MPFR_RNDNA = -1
} mpfr_rnd_t;
typedef struct {
    mpfr_prec_t _mpfr_prec;
    mpfr_sign_t _mpfr_sign;
    mpfr_exp_t _mpfr_exp;
    mp_limb_t* _mpfr_d;
} __mpfr_struct;
typedef __mpfr_struct mpfr_t[1];
typedef __mpfr_struct* mpfr_ptr;
typedef const __mpfr_struct* mpfr_srcptr;

void mpfr_init2(mpfr_ptr, mpfr_prec_t);
void mpfr_clear(mpfr_ptr);
void mpfr_set_zero(mpfr_ptr, int);
int mpfr_set(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
void mpfr_set_prec(mpfr_ptr, mpfr_prec_t);
void mpfr_swap(mpfr_ptr, mpfr_ptr);
int mpfr_set_str(mpfr_ptr, const char*, int, mpfr_rnd_t);
int mpfr_asprintf(char**, const char*, ...);
void mpfr_free_str(char*);
int mpfr_set_ui(mpfr_ptr, unsigned long, mpfr_rnd_t);
int mpfr_set_si(mpfr_ptr, long, mpfr_rnd_t);
int mpfr_log10(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_log(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_div_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_add(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_sub(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_mul_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_sqr(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_mul(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_div(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_ui_div(mpfr_ptr, unsigned long, mpfr_srcptr, mpfr_rnd_t);
int mpfr_add_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_sub_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_ui_sub(mpfr_ptr, unsigned long, mpfr_srcptr, mpfr_rnd_t);
int mpfr_neg(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_cmp(mpfr_srcptr, mpfr_srcptr);
int mpfr_cmp_ui(mpfr_srcptr, unsigned long);
int mpfr_sgn(mpfr_srcptr);
int mpfr_pow_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
const char* mpfr_get_version(void);
}
#endif
#include <omp.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

mpfr_prec_t g_prec = 384;

struct Real {
    mpfr_t v;

    Real() {
        mpfr_init2(v, g_prec);
        mpfr_set_zero(v, 0);
    }

    Real(const Real& other) {
        mpfr_init2(v, g_prec);
        mpfr_set(v, other.v, MPFR_RNDN);
    }

    Real& operator=(const Real& other) {
        if (this != &other) {
            mpfr_set_prec(v, g_prec);
            mpfr_set(v, other.v, MPFR_RNDN);
        }
        return *this;
    }

    Real(Real&& other) noexcept {
        mpfr_init2(v, g_prec);
        mpfr_swap(v, other.v);
    }

    Real& operator=(Real&& other) noexcept {
        if (this != &other) {
            mpfr_swap(v, other.v);
        }
        return *this;
    }

    ~Real() {
        mpfr_clear(v);
    }
};

struct FactorialLogs {
    std::vector<Real> lo;
    std::vector<Real> hi;

    explicit FactorialLogs(int max_n)
        : lo(static_cast<std::size_t>(max_n + 1)),
          hi(static_cast<std::size_t>(max_n + 1)) {}
};

struct Constants {
    Real log10_2_lo;
    Real log10_2_hi;
    Real log10_4_lo;
    Real log10_44_7_hi;
    Real log10_7_5_lo;
    Real log10_80_hi;
    Real log10_e_hi;
    Real central_t_lo;
    Real central_t_hi;
    Real quartic_r_lo;
    Real quartic_r_hi;
    Real quartic_s_lo;
    Real quartic_s_hi;
    Real quartic_moment_hi;
    Real quartic_denom_lo;
};

struct Params {
    int b_rank_lo = 14;
    int b_rank_hi = 217;
    int c_rank_lo = 20;
    int c_rank_hi = 217;
    int b_radius_slope = 22;
    int c_radius_slope = 19;
    int b_radius_offset = 0;
    int c_radius_offset = 0;
    int search_lo = 63;
    int search_hi = 40001;
    int slack = 32;
    int log10_neg_power = 15;
    bool include_b = true;
    bool include_c = true;
    bool progress = false;
};

struct Row {
    char family = 'B';
    int rank = 0;
};

struct RowResult {
    char family = 'B';
    int rank = 0;
    int onset = 0;
    int slack_target = 0;
    int radius = 0;
    int ok = 0;
    std::string reason;
    Real c_push_lower;
    Real onset_margin_lower;
    Real slack_margin_lower;
};

[[noreturn]] void die(const std::string& message, int code = 2) {
    std::cerr << message << "\n";
    std::exit(code);
}

int next_odd_at_least(int n) {
    return (n % 2 == 0) ? n + 1 : n;
}

int prev_odd_at_most(int n) {
    return (n % 2 == 0) ? n - 1 : n;
}

void set_decimal(Real& out, const char* value, mpfr_rnd_t rnd) {
    if (mpfr_set_str(out.v, value, 10, rnd) != 0) {
        die(std::string("invalid decimal constant: ") + value);
    }
}

std::string format_mpfr(const mpfr_t value, int digits = 24) {
    char fmt[32];
    std::snprintf(fmt, sizeof(fmt), "%%.%dRg", digits);
    char* raw = nullptr;
    if (mpfr_asprintf(&raw, fmt, value) < 0 || raw == nullptr) {
        return "<format-error>";
    }
    std::string out(raw);
    mpfr_free_str(raw);
    return out;
}

void log10_ui(Real& out, unsigned long value, mpfr_rnd_t rnd) {
    Real tmp;
    mpfr_set_ui(tmp.v, value, rnd);
    mpfr_log10(out.v, tmp.v, rnd);
}

void log10_rational(Real& out, unsigned long num, unsigned long den, mpfr_rnd_t rnd) {
    Real tmp;
    mpfr_set_ui(tmp.v, num, rnd);
    mpfr_div_ui(tmp.v, tmp.v, den, rnd);
    mpfr_log10(out.v, tmp.v, rnd);
}

FactorialLogs precompute_factorial_logs(int max_n, bool progress) {
    FactorialLogs facts(max_n);
    mpfr_set_zero(facts.lo[0].v, 0);
    mpfr_set_zero(facts.hi[0].v, 0);
    for (int k = 1; k <= max_n; ++k) {
        Real lower_log;
        Real upper_log;
        log10_ui(lower_log, static_cast<unsigned long>(k), MPFR_RNDD);
        log10_ui(upper_log, static_cast<unsigned long>(k), MPFR_RNDU);
        mpfr_add(
            facts.lo[static_cast<std::size_t>(k)].v,
            facts.lo[static_cast<std::size_t>(k - 1)].v,
            lower_log.v,
            MPFR_RNDD
        );
        mpfr_add(
            facts.hi[static_cast<std::size_t>(k)].v,
            facts.hi[static_cast<std::size_t>(k - 1)].v,
            upper_log.v,
            MPFR_RNDU
        );
        if (progress && (k == max_n || k % 10000 == 0)) {
            std::cout << "factorial_progress k=" << k << "/" << max_n << std::endl;
            std::fflush(stdout);
        }
    }
    return facts;
}

void compute_quartic_constants(Constants& constants) {
    set_decimal(constants.central_t_lo, "1.64578692310241648", MPFR_RNDD);
    set_decimal(constants.central_t_hi, "1.64578692310241648", MPFR_RNDU);
    set_decimal(constants.quartic_r_lo, "3.06396669928381904", MPFR_RNDD);
    set_decimal(constants.quartic_r_hi, "3.06396669928381904", MPFR_RNDU);
    set_decimal(constants.quartic_s_lo, "0.317259520411341256", MPFR_RNDD);
    set_decimal(constants.quartic_s_hi, "0.317259520411341256", MPFR_RNDU);

    Real sum_lo;
    Real neg_term_hi;
    mpfr_add(sum_lo.v, constants.quartic_r_lo.v, constants.quartic_s_lo.v, MPFR_RNDD);
    mpfr_mul_ui(neg_term_hi.v, sum_lo.v, 2, MPFR_RNDD);
    mpfr_neg(neg_term_hi.v, neg_term_hi.v, MPFR_RNDU);

    Real r2_hi;
    Real s2_hi;
    Real rs_hi;
    Real four_rs_hi;
    Real positive_quad_hi;
    mpfr_sqr(r2_hi.v, constants.quartic_r_hi.v, MPFR_RNDU);
    mpfr_sqr(s2_hi.v, constants.quartic_s_hi.v, MPFR_RNDU);
    mpfr_mul(rs_hi.v, constants.quartic_r_hi.v, constants.quartic_s_hi.v, MPFR_RNDU);
    mpfr_mul_ui(four_rs_hi.v, rs_hi.v, 4, MPFR_RNDU);
    mpfr_add(positive_quad_hi.v, r2_hi.v, s2_hi.v, MPFR_RNDU);
    mpfr_add(positive_quad_hi.v, positive_quad_hi.v, four_rs_hi.v, MPFR_RNDU);

    Real r2s2_hi;
    mpfr_mul(r2s2_hi.v, r2_hi.v, s2_hi.v, MPFR_RNDU);

    mpfr_set_ui(constants.quartic_moment_hi.v, 6, MPFR_RNDU);
    mpfr_add(constants.quartic_moment_hi.v, constants.quartic_moment_hi.v, neg_term_hi.v, MPFR_RNDU);
    mpfr_add(constants.quartic_moment_hi.v, constants.quartic_moment_hi.v, positive_quad_hi.v, MPFR_RNDU);
    mpfr_add(constants.quartic_moment_hi.v, constants.quartic_moment_hi.v, r2s2_hi.v, MPFR_RNDU);

    Real tr_lo;
    Real ts_lo;
    Real tr2_lo;
    Real ts2_lo;
    mpfr_add(tr_lo.v, constants.central_t_lo.v, constants.quartic_r_lo.v, MPFR_RNDD);
    mpfr_add(ts_lo.v, constants.central_t_lo.v, constants.quartic_s_lo.v, MPFR_RNDD);
    mpfr_sqr(tr2_lo.v, tr_lo.v, MPFR_RNDD);
    mpfr_sqr(ts2_lo.v, ts_lo.v, MPFR_RNDD);
    mpfr_mul(constants.quartic_denom_lo.v, tr2_lo.v, ts2_lo.v, MPFR_RNDD);
}

Constants make_constants() {
    Constants constants;
    log10_ui(constants.log10_2_lo, 2, MPFR_RNDD);
    log10_ui(constants.log10_2_hi, 2, MPFR_RNDU);
    mpfr_mul_ui(constants.log10_4_lo.v, constants.log10_2_lo.v, 2, MPFR_RNDD);
    log10_rational(constants.log10_44_7_hi, 44, 7, MPFR_RNDU);
    log10_rational(constants.log10_7_5_lo, 7, 5, MPFR_RNDD);
    log10_ui(constants.log10_80_hi, 80, MPFR_RNDU);

    Real ten;
    Real ln10_lo;
    mpfr_set_ui(ten.v, 10, MPFR_RNDD);
    mpfr_log(ln10_lo.v, ten.v, MPFR_RNDD);
    mpfr_ui_div(constants.log10_e_hi.v, 1, ln10_lo.v, MPFR_RNDU);

    compute_quartic_constants(constants);
    return constants;
}

int dimension_for(const Row& row) {
    return row.rank * (2 * row.rank + 1);
}

int kappa_for(const Row& row) {
    if (row.family == 'B') return 2 * row.rank - 1;
    if (row.family == 'C') return row.rank + 1;
    die("unsupported family");
}

int radius_for(const Row& row, const Params& params) {
    if (row.family == 'B') return params.b_radius_slope * row.rank + params.b_radius_offset;
    if (row.family == 'C') return params.c_radius_slope * row.rank + params.c_radius_offset;
    die("unsupported family");
}

void density_factor_lower(
    Real& out,
    int rank,
    int dimension,
    const Constants& constants,
    const FactorialLogs& facts
) {
    mpfr_set_zero(out.v, 0);
    if (rank % 2 == 1 && dimension % 2 == 1) {
        const int half_index = (dimension + 1) / 2;
        const int pi_power = (rank + 1) / 2;
        Real tmp;
        mpfr_set(out.v, constants.log10_7_5_lo.v, MPFR_RNDD);
        mpfr_mul_ui(tmp.v, constants.log10_4_lo.v, static_cast<unsigned long>(half_index), MPFR_RNDD);
        mpfr_add(out.v, out.v, tmp.v, MPFR_RNDD);
        mpfr_add(out.v, out.v, facts.lo[static_cast<std::size_t>(half_index)].v, MPFR_RNDD);
        mpfr_mul_ui(tmp.v, constants.log10_44_7_hi.v, static_cast<unsigned long>(pi_power), MPFR_RNDU);
        mpfr_sub(out.v, out.v, tmp.v, MPFR_RNDD);
        mpfr_sub(out.v, out.v, facts.hi[static_cast<std::size_t>(2 * half_index)].v, MPFR_RNDD);
        return;
    }
    if (rank % 2 == 0 && dimension % 2 == 0) {
        Real tmp;
        mpfr_mul_ui(tmp.v, constants.log10_44_7_hi.v, static_cast<unsigned long>(rank / 2), MPFR_RNDU);
        mpfr_sub(out.v, out.v, tmp.v, MPFR_RNDD);
        mpfr_sub(out.v, out.v, facts.hi[static_cast<std::size_t>(dimension / 2)].v, MPFR_RNDD);
        return;
    }
    die("invalid B/C density parity");
}

bool sine_exponent_upper(Real& out, int radius, int kappa, std::string& reason) {
    Real term1;
    Real term2;
    Real term3;
    Real numerator;
    Real q_hi;
    Real one_minus_q_lo;
    Real denom_lo;

    mpfr_set_ui(term1.v, static_cast<unsigned long>(radius), MPFR_RNDU);
    mpfr_div_ui(term1.v, term1.v, 12, MPFR_RNDU);

    mpfr_set_ui(numerator.v, static_cast<unsigned long>(radius), MPFR_RNDU);
    mpfr_sqr(numerator.v, numerator.v, MPFR_RNDU);
    mpfr_set(term2.v, numerator.v, MPFR_RNDU);
    mpfr_div_ui(term2.v, term2.v, static_cast<unsigned long>(1440 * kappa), MPFR_RNDU);

    mpfr_set_ui(q_hi.v, static_cast<unsigned long>(2 * radius), MPFR_RNDU);
    mpfr_div_ui(q_hi.v, q_hi.v, static_cast<unsigned long>(39 * kappa), MPFR_RNDU);
    mpfr_ui_sub(one_minus_q_lo.v, 1, q_hi.v, MPFR_RNDD);
    if (mpfr_sgn(one_minus_q_lo.v) <= 0) {
        reason = "sine q upper bound is not below one";
        return false;
    }

    const unsigned long denom_int =
        static_cast<unsigned long>(kappa)
        * static_cast<unsigned long>(kappa)
        * 90720UL;
    mpfr_mul_ui(denom_lo.v, one_minus_q_lo.v, denom_int, MPFR_RNDD);
    if (mpfr_sgn(denom_lo.v) <= 0) {
        reason = "sine denominator lower bound is not positive";
        return false;
    }

    mpfr_set_ui(numerator.v, static_cast<unsigned long>(radius), MPFR_RNDU);
    mpfr_pow_ui(numerator.v, numerator.v, 3, MPFR_RNDU);
    mpfr_div(term3.v, numerator.v, denom_lo.v, MPFR_RNDU);

    mpfr_add(out.v, term1.v, term2.v, MPFR_RNDU);
    mpfr_add(out.v, out.v, term3.v, MPFR_RNDU);
    return true;
}

void base_log_lower(
    Real& out,
    int rank,
    int dimension,
    const Constants& constants,
    const FactorialLogs& facts
) {
    Real weyl_hi;
    Real degree_lo;
    Real density_lo;

    mpfr_mul_ui(weyl_hi.v, constants.log10_2_hi.v, static_cast<unsigned long>(rank), MPFR_RNDU);
    mpfr_add(weyl_hi.v, weyl_hi.v, facts.hi[static_cast<std::size_t>(rank)].v, MPFR_RNDU);

    mpfr_set_zero(degree_lo.v, 0);
    for (int j = 1; j <= rank; ++j) {
        mpfr_add(degree_lo.v, degree_lo.v, facts.lo[static_cast<std::size_t>(2 * j)].v, MPFR_RNDD);
    }

    density_factor_lower(density_lo, rank, dimension, constants, facts);

    mpfr_set_zero(out.v, 0);
    mpfr_sub(out.v, out.v, weyl_hi.v, MPFR_RNDD);
    mpfr_add(out.v, out.v, degree_lo.v, MPFR_RNDD);
    mpfr_add(out.v, out.v, density_lo.v, MPFR_RNDD);
}

bool push_weight_lower(
    Real& out,
    int x_floor,
    const Constants& constants,
    std::string& reason
) {
    Real x_plus_t_hi;
    Real x_plus_t_sq_hi;
    Real numerator_hi;
    Real tail_upper;
    Real x_sq_plus_one;
    Real deriv_term_hi;
    Real deriv_lower;
    Real second_term_hi;
    Real second_lower;

    mpfr_set_si(x_plus_t_hi.v, x_floor, MPFR_RNDU);
    mpfr_add(x_plus_t_hi.v, x_plus_t_hi.v, constants.central_t_hi.v, MPFR_RNDU);
    mpfr_sqr(x_plus_t_sq_hi.v, x_plus_t_hi.v, MPFR_RNDU);
    mpfr_mul(numerator_hi.v, x_plus_t_sq_hi.v, constants.quartic_moment_hi.v, MPFR_RNDU);
    mpfr_div(tail_upper.v, numerator_hi.v, constants.quartic_denom_lo.v, MPFR_RNDU);

    mpfr_set_si(x_sq_plus_one.v, x_floor, MPFR_RNDD);
    mpfr_sqr(x_sq_plus_one.v, x_sq_plus_one.v, MPFR_RNDD);
    mpfr_add_ui(x_sq_plus_one.v, x_sq_plus_one.v, 1, MPFR_RNDD);
    mpfr_sub(out.v, x_sq_plus_one.v, tail_upper.v, MPFR_RNDD);

    mpfr_mul(deriv_term_hi.v, constants.quartic_moment_hi.v, x_plus_t_hi.v, MPFR_RNDU);
    mpfr_mul_ui(deriv_term_hi.v, deriv_term_hi.v, 2, MPFR_RNDU);
    mpfr_div(deriv_term_hi.v, deriv_term_hi.v, constants.quartic_denom_lo.v, MPFR_RNDU);
    mpfr_set_si(deriv_lower.v, 2 * x_floor, MPFR_RNDD);
    mpfr_sub(deriv_lower.v, deriv_lower.v, deriv_term_hi.v, MPFR_RNDD);

    mpfr_mul_ui(second_term_hi.v, constants.quartic_moment_hi.v, 2, MPFR_RNDU);
    mpfr_div(second_term_hi.v, second_term_hi.v, constants.quartic_denom_lo.v, MPFR_RNDU);
    mpfr_set_ui(second_lower.v, 2, MPFR_RNDD);
    mpfr_sub(second_lower.v, second_lower.v, second_term_hi.v, MPFR_RNDD);

    if (mpfr_sgn(out.v) <= 0) {
        reason = "push weight lower bound is not positive";
        return false;
    }
    if (mpfr_sgn(deriv_lower.v) <= 0) {
        reason = "push weight derivative lower bound is not positive";
        return false;
    }
    if (mpfr_sgn(second_lower.v) <= 0) {
        reason = "push weight second derivative lower bound is not positive";
        return false;
    }
    return true;
}

bool margin_lower(
    Real& out,
    Real& c_push_lower_out,
    int& radius_out,
    const Row& row,
    int target_n,
    const Params& params,
    const Constants& constants,
    const FactorialLogs& facts,
    std::string& reason
) {
    const int rank = row.rank;
    const int dimension = dimension_for(row);
    const int kappa = kappa_for(row);
    const int c_value = rank;
    const int radius = radius_for(row, params);
    radius_out = radius;
    if (radius <= 0) {
        reason = "nonpositive radius";
        return false;
    }
    const int x_floor = dimension - radius;
    if (x_floor <= 0) {
        reason = "nonpositive x floor";
        return false;
    }

    mpfr_set_si(c_push_lower_out.v, x_floor, MPFR_RNDD);
    mpfr_sub(c_push_lower_out.v, c_push_lower_out.v, constants.central_t_hi.v, MPFR_RNDD);
    Real two_c;
    mpfr_set_ui(two_c.v, static_cast<unsigned long>(2 * c_value), MPFR_RNDD);
    if (mpfr_cmp(c_push_lower_out.v, two_c.v) <= 0) {
        reason = "c_push lower bound is not above 2c";
        return false;
    }
    if (mpfr_cmp_ui(c_push_lower_out.v, 80) <= 0) {
        reason = "c_push lower bound is not above cutoff 80";
        return false;
    }

    Real base_lo;
    Real sine_exp_hi;
    Real sine_log_lo;
    Real radial_ratio_lo;
    Real radial_log_lo;
    Real radial_scaled_lo;
    Real ball_log_lo;
    Real push_lo;
    Real log_push_lo;
    Real weighted_log_lo;
    Real target_a_hi;
    Real target_b_hi;
    Real target_hi;
    Real log_2c_hi;
    Real log_cpush_lo;
    Real tmp;

    base_log_lower(base_lo, rank, dimension, constants, facts);

    if (!sine_exponent_upper(sine_exp_hi, radius, kappa, reason)) {
        return false;
    }
    mpfr_mul(sine_log_lo.v, sine_exp_hi.v, constants.log10_e_hi.v, MPFR_RNDU);
    mpfr_neg(sine_log_lo.v, sine_log_lo.v, MPFR_RNDD);

    mpfr_set_ui(radial_ratio_lo.v, static_cast<unsigned long>(radius), MPFR_RNDD);
    mpfr_div_ui(radial_ratio_lo.v, radial_ratio_lo.v, static_cast<unsigned long>(2 * kappa), MPFR_RNDD);
    if (mpfr_cmp_ui(radial_ratio_lo.v, 1) <= 0) {
        reason = "radial ratio lower bound is not above one";
        return false;
    }
    mpfr_log10(radial_log_lo.v, radial_ratio_lo.v, MPFR_RNDD);
    mpfr_mul_ui(radial_scaled_lo.v, radial_log_lo.v, static_cast<unsigned long>(dimension), MPFR_RNDD);
    mpfr_div_ui(radial_scaled_lo.v, radial_scaled_lo.v, 2, MPFR_RNDD);

    mpfr_add(ball_log_lo.v, base_lo.v, sine_log_lo.v, MPFR_RNDD);
    mpfr_add(ball_log_lo.v, ball_log_lo.v, radial_scaled_lo.v, MPFR_RNDD);

    if (!push_weight_lower(push_lo, x_floor, constants, reason)) {
        return false;
    }
    mpfr_log10(log_push_lo.v, push_lo.v, MPFR_RNDD);
    mpfr_add(weighted_log_lo.v, ball_log_lo.v, log_push_lo.v, MPFR_RNDD);

    log10_ui(log_2c_hi, static_cast<unsigned long>(2 * c_value), MPFR_RNDU);
    mpfr_mul_ui(target_a_hi.v, log_2c_hi.v, static_cast<unsigned long>(target_n), MPFR_RNDU);
    mpfr_sub_ui(target_a_hi.v, target_a_hi.v, static_cast<unsigned long>(params.log10_neg_power), MPFR_RNDU);

    mpfr_mul_ui(target_b_hi.v, constants.log10_80_hi.v, static_cast<unsigned long>(target_n), MPFR_RNDU);
    mpfr_add(target_b_hi.v, target_b_hi.v, constants.log10_2_hi.v, MPFR_RNDU);

    if (mpfr_cmp(target_a_hi.v, target_b_hi.v) >= 0) {
        mpfr_set(target_hi.v, target_a_hi.v, MPFR_RNDU);
    } else {
        mpfr_set(target_hi.v, target_b_hi.v, MPFR_RNDU);
    }
    mpfr_add(target_hi.v, target_hi.v, constants.log10_2_hi.v, MPFR_RNDU);

    mpfr_log10(log_cpush_lo.v, c_push_lower_out.v, MPFR_RNDD);
    mpfr_mul_ui(tmp.v, log_cpush_lo.v, static_cast<unsigned long>(target_n), MPFR_RNDD);
    mpfr_sub(target_hi.v, target_hi.v, tmp.v, MPFR_RNDU);

    mpfr_sub(out.v, weighted_log_lo.v, target_hi.v, MPFR_RNDD);
    return true;
}

bool margin_is_positive(
    const Row& row,
    int target_n,
    const Params& params,
    const Constants& constants,
    const FactorialLogs& facts,
    Real& margin,
    Real& c_push,
    int& radius,
    std::string& reason
) {
    if (!margin_lower(margin, c_push, radius, row, target_n, params, constants, facts, reason)) {
        return false;
    }
    return mpfr_sgn(margin.v) > 0;
}

bool find_onset(
    int& onset,
    Real& onset_margin,
    Real& c_push,
    int& radius,
    const Row& row,
    const Params& params,
    const Constants& constants,
    const FactorialLogs& facts,
    std::string& reason
) {
    const int lo = next_odd_at_least(params.search_lo);
    const int hi = prev_odd_at_most(params.search_hi);
    if (lo <= 0 || hi < lo) {
        reason = "invalid search interval";
        return false;
    }

    Real lo_margin;
    Real lo_c_push;
    int lo_radius = 0;
    std::string lo_reason;
    if (margin_is_positive(row, lo, params, constants, facts, lo_margin, lo_c_push, lo_radius, lo_reason)) {
        onset = lo;
        onset_margin = lo_margin;
        c_push = lo_c_push;
        radius = lo_radius;
        return true;
    }

    Real hi_margin;
    Real hi_c_push;
    int hi_radius = 0;
    if (!margin_is_positive(row, hi, params, constants, facts, hi_margin, hi_c_push, hi_radius, reason)) {
        if (reason.empty()) reason = "no positive margin at search_hi";
        return false;
    }

    int false_k = (lo - 1) / 2;
    int true_k = (hi - 1) / 2;
    Real true_margin = hi_margin;
    Real true_c_push = hi_c_push;
    int true_radius = hi_radius;

    while (false_k + 1 < true_k) {
        int mid_k = (false_k + true_k) / 2;
        if (mid_k <= false_k) ++mid_k;
        if (mid_k >= true_k) --mid_k;
        const int mid = 2 * mid_k + 1;
        Real mid_margin;
        Real mid_c_push;
        int mid_radius = 0;
        std::string mid_reason;
        if (margin_is_positive(row, mid, params, constants, facts, mid_margin, mid_c_push, mid_radius, mid_reason)) {
            true_k = mid_k;
            true_margin = mid_margin;
            true_c_push = mid_c_push;
            true_radius = mid_radius;
        } else {
            false_k = mid_k;
        }
    }

    onset = 2 * true_k + 1;
    onset_margin = true_margin;
    c_push = true_c_push;
    radius = true_radius;
    return true;
}

RowResult verify_row(
    const Row& row,
    const Params& params,
    const Constants& constants,
    const FactorialLogs& facts
) {
    RowResult result;
    result.family = row.family;
    result.rank = row.rank;

    if (!find_onset(
            result.onset,
            result.onset_margin_lower,
            result.c_push_lower,
            result.radius,
            row,
            params,
            constants,
            facts,
            result.reason
        )) {
        return result;
    }

    result.slack_target = result.onset + params.slack;
    if (result.slack_target % 2 == 0) ++result.slack_target;
    Real slack_c_push_lower;
    int slack_radius = 0;
    if (!margin_lower(
            result.slack_margin_lower,
            slack_c_push_lower,
            slack_radius,
            row,
            result.slack_target,
            params,
            constants,
            facts,
            result.reason
        )) {
        return result;
    }
    if (mpfr_sgn(result.slack_margin_lower.v) <= 0) {
        result.reason = "slack margin lower bound is not positive";
        return result;
    }
    result.ok = 1;
    return result;
}

std::vector<Row> make_rows(const Params& params) {
    std::vector<Row> rows;
    if (params.include_b) {
        for (int rank = params.b_rank_lo; rank <= params.b_rank_hi; ++rank) {
            rows.push_back({'B', rank});
        }
    }
    if (params.include_c) {
        for (int rank = params.c_rank_lo; rank <= params.c_rank_hi; ++rank) {
            rows.push_back({'C', rank});
        }
    }
    return rows;
}

void print_usage(const char* argv0) {
    std::cerr
        << "usage: " << argv0
        << " [--b-rank-lo N] [--b-rank-hi N] [--c-rank-lo N] [--c-rank-hi N]"
        << " [--only-b] [--only-c]"
        << " [--b-radius-slope N] [--c-radius-slope N]"
        << " [--b-radius-offset N] [--c-radius-offset N]"
        << " [--search-lo N] [--search-hi N] [--slack N]"
        << " [--log10-neg-power N] [--precision BITS] [--progress]\n";
}

}  // namespace

int main(int argc, char** argv) {
    Params params;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--b-rank-lo" && i + 1 < argc) {
            params.b_rank_lo = std::atoi(argv[++i]);
        } else if (arg == "--b-rank-hi" && i + 1 < argc) {
            params.b_rank_hi = std::atoi(argv[++i]);
        } else if (arg == "--c-rank-lo" && i + 1 < argc) {
            params.c_rank_lo = std::atoi(argv[++i]);
        } else if (arg == "--c-rank-hi" && i + 1 < argc) {
            params.c_rank_hi = std::atoi(argv[++i]);
        } else if (arg == "--only-b") {
            params.include_b = true;
            params.include_c = false;
        } else if (arg == "--only-c") {
            params.include_b = false;
            params.include_c = true;
        } else if (arg == "--b-radius-slope" && i + 1 < argc) {
            params.b_radius_slope = std::atoi(argv[++i]);
        } else if (arg == "--c-radius-slope" && i + 1 < argc) {
            params.c_radius_slope = std::atoi(argv[++i]);
        } else if (arg == "--b-radius-offset" && i + 1 < argc) {
            params.b_radius_offset = std::atoi(argv[++i]);
        } else if (arg == "--c-radius-offset" && i + 1 < argc) {
            params.c_radius_offset = std::atoi(argv[++i]);
        } else if (arg == "--search-lo" && i + 1 < argc) {
            params.search_lo = std::atoi(argv[++i]);
        } else if (arg == "--search-hi" && i + 1 < argc) {
            params.search_hi = std::atoi(argv[++i]);
        } else if (arg == "--slack" && i + 1 < argc) {
            params.slack = std::atoi(argv[++i]);
        } else if (arg == "--log10-neg-power" && i + 1 < argc) {
            params.log10_neg_power = std::atoi(argv[++i]);
        } else if (arg == "--precision" && i + 1 < argc) {
            g_prec = static_cast<mpfr_prec_t>(std::max(128, std::atoi(argv[++i])));
        } else if (arg == "--progress") {
            params.progress = true;
        } else {
            print_usage(argv[0]);
            return 2;
        }
    }

    if ((!params.include_b && !params.include_c)
        || params.b_rank_lo < 1
        || params.c_rank_lo < 1
        || params.b_rank_hi < params.b_rank_lo
        || params.c_rank_hi < params.c_rank_lo
        || params.search_lo <= 0
        || params.search_hi < params.search_lo
        || params.slack < 0
        || params.log10_neg_power < 0) {
        print_usage(argv[0]);
        return 2;
    }

    const std::vector<Row> rows = make_rows(params);
    int max_dimension = 0;
    int max_rank = 0;
    for (const Row& row : rows) {
        max_dimension = std::max(max_dimension, dimension_for(row));
        max_rank = std::max(max_rank, row.rank);
    }
    const int max_fact = std::max(max_dimension + 1, 2 * max_rank);

    std::cout << "B/C interval fixed-template direct-tail MPFR onset verifier\n"
              << "b_rank_lo=" << params.b_rank_lo
              << " b_rank_hi=" << params.b_rank_hi
              << " c_rank_lo=" << params.c_rank_lo
              << " c_rank_hi=" << params.c_rank_hi
              << " include_b=" << params.include_b
              << " include_c=" << params.include_c
              << " b_radius_slope=" << params.b_radius_slope
              << " c_radius_slope=" << params.c_radius_slope
              << " b_radius_offset=" << params.b_radius_offset
              << " c_radius_offset=" << params.c_radius_offset
              << " search_lo=" << params.search_lo
              << " search_hi=" << params.search_hi
              << " slack=" << params.slack
              << " cheb_degree=8 cheb_scale=12000 cutoff=80 majorant_degree=32"
              << " log10_neg=-" << params.log10_neg_power
              << " mpfr_precision=" << g_prec
              << " mpfr_version=" << mpfr_get_version()
              << " max_fact=" << max_fact
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;
    std::cout << "CONDITIONAL_INPUT "
              << "requires uniform B/C Chebyshev negative-tail bound "
              << "degree=8 scale=12000 cutoff=80 majorant_degree=32 "
              << "log10_negative_upper<=-" << params.log10_neg_power
              << std::endl;

    std::cout << "factorial_log_precompute_start max_fact=" << max_fact << std::endl;
    const FactorialLogs facts = precompute_factorial_logs(max_fact, params.progress);
    std::cout << "factorial_log_precompute_done max_fact=" << max_fact << std::endl;

    const Constants constants = make_constants();
    std::vector<RowResult> results(rows.size());
    int completed = 0;

#pragma omp parallel for schedule(dynamic, 1)
    for (int index = 0; index < static_cast<int>(rows.size()); ++index) {
        const Row row = rows[static_cast<std::size_t>(index)];
        RowResult result = verify_row(row, params, constants, facts);
        results[static_cast<std::size_t>(index)] = result;

        int done_now = 0;
#pragma omp atomic capture
        done_now = ++completed;
        if (params.progress) {
#pragma omp critical
            {
                std::cout << "row_done " << result.family << "_" << result.rank
                          << " onset=" << result.onset
                          << " ok=" << result.ok
                          << " margin_lower_log10="
                          << format_mpfr(result.onset_margin_lower.v, 18)
                          << " slack_margin_lower_log10="
                          << format_mpfr(result.slack_margin_lower.v, 18)
                          << " completed=" << done_now << "/" << rows.size();
                if (!result.ok) {
                    std::cout << " reason=\"" << result.reason << "\"";
                }
                std::cout << std::endl;
                std::fflush(stdout);
            }
        }
    }

    int failures = 0;
    bool have_worst_onset = false;
    bool have_worst_slack = false;
    RowResult worst_onset;
    RowResult worst_slack;
    for (const RowResult& result : results) {
        if (!result.ok) {
            ++failures;
            continue;
        }
        if (!have_worst_onset
            || mpfr_cmp(result.onset_margin_lower.v, worst_onset.onset_margin_lower.v) < 0) {
            worst_onset = result;
            have_worst_onset = true;
        }
        if (!have_worst_slack
            || mpfr_cmp(result.slack_margin_lower.v, worst_slack.slack_margin_lower.v) < 0) {
            worst_slack = result;
            have_worst_slack = true;
        }
    }

    std::cout << "row\tonset\tmargin_lower_log10\tslack_target"
              << "\tslack_margin_lower_log10\tradius\tc_push_lower\n";
    for (const RowResult& result : results) {
        std::cout << result.family << "_" << result.rank << '\t'
                  << result.onset << '\t'
                  << format_mpfr(result.onset_margin_lower.v, 24) << '\t'
                  << result.slack_target << '\t'
                  << format_mpfr(result.slack_margin_lower.v, 24) << '\t'
                  << result.radius << '\t'
                  << format_mpfr(result.c_push_lower.v, 24);
        if (!result.ok) {
            std::cout << "\tFAIL\t" << result.reason;
        }
        std::cout << '\n';
    }

    std::cout << "SUMMARY failures=" << failures;
    if (have_worst_onset) {
        std::cout << " worst_onset=" << worst_onset.family << "_" << worst_onset.rank
                  << " onset=" << worst_onset.onset
                  << " margin_lower_log10="
                  << format_mpfr(worst_onset.onset_margin_lower.v, 24);
    }
    if (have_worst_slack) {
        std::cout << " worst_slack=" << worst_slack.family << "_" << worst_slack.rank
                  << " target=" << worst_slack.slack_target
                  << " margin_lower_log10="
                  << format_mpfr(worst_slack.slack_margin_lower.v, 24);
    }
    std::cout << std::endl;
    std::fflush(stdout);

    return failures == 0 ? 0 : 1;
}
