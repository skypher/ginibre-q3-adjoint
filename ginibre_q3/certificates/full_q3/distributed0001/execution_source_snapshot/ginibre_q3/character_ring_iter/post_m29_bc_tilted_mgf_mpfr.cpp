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
void mpfr_swap(mpfr_ptr, mpfr_ptr);
int mpfr_set_ui(mpfr_ptr, unsigned long, mpfr_rnd_t);
int mpfr_add(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_sub(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_mul(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_mul_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_div(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_div_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_sqr(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_sqrt(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_log(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_exp(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_neg(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_cmp(mpfr_srcptr, mpfr_srcptr);
int mpfr_sgn(mpfr_srcptr);
int mpfr_asprintf(char**, const char*, ...);
void mpfr_free_str(char*);
const char* mpfr_get_version(void);
}
#endif

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr int B_MIN = 62;
constexpr int B_MAX = 123;
constexpr int C_MIN = 62;
constexpr int C_MAX = 217;
constexpr unsigned long R_NUM = 20001;
constexpr unsigned long R_DEN = 10000;
constexpr unsigned long ALPHA_NUM = 979;
constexpr unsigned long ALPHA_DEN = 250;
constexpr unsigned long H_NUM = 601;
constexpr unsigned long H_DEN = 500;
constexpr mpfr_prec_t PRECISION = 384;

struct Real {
    mpfr_t value;

    Real() {
        mpfr_init2(value, PRECISION);
        mpfr_set_zero(value, 0);
    }

    Real(const Real& other) {
        mpfr_init2(value, PRECISION);
        mpfr_set(value, other.value, MPFR_RNDN);
    }

    Real& operator=(const Real& other) {
        if (this != &other) mpfr_set(value, other.value, MPFR_RNDN);
        return *this;
    }

    Real(Real&& other) noexcept {
        mpfr_init2(value, PRECISION);
        mpfr_swap(value, other.value);
    }

    Real& operator=(Real&& other) noexcept {
        if (this != &other) mpfr_swap(value, other.value);
        return *this;
    }

    ~Real() { mpfr_clear(value); }
};

struct RowMargins {
    int rank = 0;
    Real domain;
    Real square_domain;
    Real tilted_mass;
    Real tilted_over_twice_tau;
    Real push_over_twice_negative;
    Real push_over_twice_central;
};

[[noreturn]] void die(const std::string& message) {
    std::cerr << message << '\n';
    std::exit(2);
}

std::string format_real(const Real& x, int digits = 28) {
    char format[32];
    std::snprintf(format, sizeof(format), "%%.%dRg", digits);
    char* raw = nullptr;
    if (mpfr_asprintf(&raw, format, x.value) < 0 || raw == nullptr) {
        return "<format-error>";
    }
    std::string answer(raw);
    mpfr_free_str(raw);
    return answer;
}

void set_rational(
    Real& out,
    unsigned long numerator,
    unsigned long denominator,
    mpfr_rnd_t rounding
) {
    mpfr_set_ui(out.value, numerator, rounding);
    mpfr_div_ui(out.value, out.value, denominator, rounding);
}

void log_ui(Real& out, unsigned long value, mpfr_rnd_t rounding) {
    mpfr_set_ui(out.value, value, rounding);
    mpfr_log(out.value, out.value, rounding);
}

void log_rational(
    Real& out,
    unsigned long numerator,
    unsigned long denominator,
    mpfr_rnd_t rounding
) {
    set_rational(out, numerator, denominator, rounding);
    mpfr_log(out.value, out.value, rounding);
}

std::vector<Real> factorial_log_lowers(int maximum) {
    std::vector<Real> answer(static_cast<std::size_t>(maximum + 1));
    for (int k = 1; k <= maximum; ++k) {
        Real term;
        log_ui(term, static_cast<unsigned long>(k), MPFR_RNDD);
        mpfr_add(
            answer[static_cast<std::size_t>(k)].value,
            answer[static_cast<std::size_t>(k - 1)].value,
            term.value,
            MPFR_RNDD
        );
    }
    return answer;
}

// CJL Proposition 3.3 and Lemma 3.7 give
//   exp(v^2/2)(1-beta_b(v)) <= E exp(v T_1)
//                            <= exp(v^2/2)/(1-beta_b(v)),
// where d is the total matrix dimension and
//   beta_d(v)=exp(v^2/(d+1)) v^d /
//             ((1-v^2/(d+1)^2) d!).
// This routine returns an upward-rounded beta_d(v), given an upper bound on v.
void beta_upper(
    Real& out,
    int dimension,
    const Real& v_upper,
    const std::vector<Real>& factorial_logs
) {
    const unsigned long d = static_cast<unsigned long>(dimension);

    Real v_square_upper;
    mpfr_sqr(v_square_upper.value, v_upper.value, MPFR_RNDU);

    Real first_upper;
    mpfr_div_ui(first_upper.value, v_square_upper.value, d + 1, MPFR_RNDU);

    Real quotient_upper;
    mpfr_div_ui(quotient_upper.value, v_square_upper.value, d + 1, MPFR_RNDU);
    mpfr_div_ui(quotient_upper.value, quotient_upper.value, d + 1, MPFR_RNDU);

    Real one;
    mpfr_set_ui(one.value, 1, MPFR_RNDD);
    Real denominator_lower;
    mpfr_sub(
        denominator_lower.value,
        one.value,
        quotient_upper.value,
        MPFR_RNDD
    );
    if (mpfr_sgn(denominator_lower.value) <= 0) {
        die("CJL beta denominator is not positive");
    }
    Real minus_log_denominator_upper;
    mpfr_log(
        minus_log_denominator_upper.value,
        denominator_lower.value,
        MPFR_RNDD
    );
    mpfr_neg(
        minus_log_denominator_upper.value,
        minus_log_denominator_upper.value,
        MPFR_RNDU
    );

    Real d_log_v_upper;
    mpfr_log(d_log_v_upper.value, v_upper.value, MPFR_RNDU);
    mpfr_mul_ui(d_log_v_upper.value, d_log_v_upper.value, d, MPFR_RNDU);

    Real log_beta_upper;
    mpfr_add(
        log_beta_upper.value,
        first_upper.value,
        minus_log_denominator_upper.value,
        MPFR_RNDU
    );
    mpfr_add(
        log_beta_upper.value,
        log_beta_upper.value,
        d_log_v_upper.value,
        MPFR_RNDU
    );
    mpfr_sub(
        log_beta_upper.value,
        log_beta_upper.value,
        factorial_logs[static_cast<std::size_t>(d)].value,
        MPFR_RNDU
    );
    mpfr_exp(out.value, log_beta_upper.value, MPFR_RNDU);
}

void one_minus_beta_lower(Real& out, const Real& beta) {
    Real one;
    mpfr_set_ui(one.value, 1, MPFR_RNDD);
    mpfr_sub(out.value, one.value, beta.value, MPFR_RNDD);
    if (mpfr_sgn(out.value) <= 0) die("CJL beta upper bound reached one");
}

void reciprocal_upper(Real& out, const Real& positive_lower) {
    Real one;
    mpfr_set_ui(one.value, 1, MPFR_RNDU);
    mpfr_div(out.value, one.value, positive_lower.value, MPFR_RNDU);
}

struct Geometry {
    Real sqrt_b_lower;
    Real sqrt_b_upper;
    Real a_lower;
    Real a_upper;
    Real t_upper;
    Real s_upper;
    Real plus_upper;
};

Geometry geometry_for(int rank) {
    Geometry geometry;

    Real b;
    mpfr_set_ui(b.value, static_cast<unsigned long>(rank), MPFR_RNDN);
    mpfr_sqrt(geometry.sqrt_b_lower.value, b.value, MPFR_RNDD);
    mpfr_sqrt(geometry.sqrt_b_upper.value, b.value, MPFR_RNDU);

    Real alpha_lower;
    Real alpha_upper;
    set_rational(alpha_lower, ALPHA_NUM, ALPHA_DEN, MPFR_RNDD);
    set_rational(alpha_upper, ALPHA_NUM, ALPHA_DEN, MPFR_RNDU);
    mpfr_mul(
        geometry.a_lower.value,
        alpha_lower.value,
        geometry.sqrt_b_lower.value,
        MPFR_RNDD
    );
    mpfr_mul(
        geometry.a_upper.value,
        alpha_upper.value,
        geometry.sqrt_b_upper.value,
        MPFR_RNDU
    );

    Real twice_rb_upper;
    set_rational(
        twice_rb_upper,
        2 * R_NUM * static_cast<unsigned long>(rank),
        R_DEN,
        MPFR_RNDU
    );
    Real three_a_upper;
    mpfr_mul_ui(three_a_upper.value, geometry.a_upper.value, 3, MPFR_RNDU);
    Real threshold_square_upper;
    mpfr_add(
        threshold_square_upper.value,
        twice_rb_upper.value,
        three_a_upper.value,
        MPFR_RNDU
    );
    mpfr_sqrt(
        geometry.t_upper.value,
        threshold_square_upper.value,
        MPFR_RNDU
    );

    Real h_upper;
    set_rational(h_upper, H_NUM, H_DEN, MPFR_RNDU);
    mpfr_add(
        geometry.s_upper.value,
        geometry.t_upper.value,
        h_upper.value,
        MPFR_RNDU
    );
    mpfr_add(
        geometry.plus_upper.value,
        geometry.s_upper.value,
        h_upper.value,
        MPFR_RNDU
    );
    return geometry;
}

void domain_margin_lower(
    Real& out,
    int pair_count,
    const Real& source_upper
) {
    Real five_fourths;
    set_rational(five_fourths, 5, 4, MPFR_RNDU);
    Real negative_five_fourths;
    mpfr_neg(
        negative_five_fourths.value,
        five_fourths.value,
        MPFR_RNDD
    );
    Real exponential_lower;
    mpfr_exp(
        exponential_lower.value,
        negative_five_fourths.value,
        MPFR_RNDD
    );
    Real source_radius_lower;
    mpfr_mul_ui(
        source_radius_lower.value,
        exponential_lower.value,
        static_cast<unsigned long>(2 * pair_count),
        MPFR_RNDD
    );
    mpfr_sub(
        out.value,
        source_radius_lower.value,
        source_upper.value,
        MPFR_RNDD
    );
}

// Lower bound for the one-sided defining-trace tail obtained by tilting at
// s=t+h and retaining the interval (s-h,s+h).  The returned logarithm is
// log V_b, where P{T_1 >= t} >= V_b.
void tilted_tail_log_lower(
    Real& log_tail,
    Real& tilted_mass,
    int dimension,
    const Geometry& geometry,
    const std::vector<Real>& factorial_logs
) {
    Real beta_s_upper;
    Real beta_plus_upper;
    Real beta_minus_upper;
    beta_upper(beta_s_upper, dimension, geometry.s_upper, factorial_logs);
    beta_upper(beta_plus_upper, dimension, geometry.plus_upper, factorial_logs);
    beta_upper(beta_minus_upper, dimension, geometry.t_upper, factorial_logs);

    Real one_minus_s_lower;
    Real one_minus_plus_lower;
    Real one_minus_minus_lower;
    one_minus_beta_lower(one_minus_s_lower, beta_s_upper);
    one_minus_beta_lower(one_minus_plus_lower, beta_plus_upper);
    one_minus_beta_lower(one_minus_minus_lower, beta_minus_upper);

    Real reciprocal_s_upper;
    Real reciprocal_plus_upper;
    Real reciprocal_minus_upper;
    reciprocal_upper(reciprocal_s_upper, one_minus_s_lower);
    reciprocal_upper(reciprocal_plus_upper, one_minus_plus_lower);
    reciprocal_upper(reciprocal_minus_upper, one_minus_minus_lower);

    Real h_upper;
    set_rational(h_upper, H_NUM, H_DEN, MPFR_RNDU);
    Real exponential_upper;
    // For an upper bound on exp(-h^2/2), round h^2/2 downward before
    // negating and exponentiating upward.
    Real h_lower;
    set_rational(h_lower, H_NUM, H_DEN, MPFR_RNDD);
    Real half_h_square_lower;
    mpfr_sqr(half_h_square_lower.value, h_lower.value, MPFR_RNDD);
    mpfr_div_ui(
        half_h_square_lower.value,
        half_h_square_lower.value,
        2,
        MPFR_RNDD
    );
    Real negative_half_h_square_lower;
    mpfr_neg(
        negative_half_h_square_lower.value,
        half_h_square_lower.value,
        MPFR_RNDU
    );
    mpfr_exp(exponential_upper.value, negative_half_h_square_lower.value, MPFR_RNDU);

    Real reciprocal_sum_upper;
    mpfr_add(
        reciprocal_sum_upper.value,
        reciprocal_plus_upper.value,
        reciprocal_minus_upper.value,
        MPFR_RNDU
    );
    Real lost_mass_upper;
    mpfr_mul(
        lost_mass_upper.value,
        exponential_upper.value,
        reciprocal_s_upper.value,
        MPFR_RNDU
    );
    mpfr_mul(
        lost_mass_upper.value,
        lost_mass_upper.value,
        reciprocal_sum_upper.value,
        MPFR_RNDU
    );

    Real one;
    mpfr_set_ui(one.value, 1, MPFR_RNDD);
    mpfr_sub(
        tilted_mass.value,
        one.value,
        lost_mass_upper.value,
        MPFR_RNDD
    );
    if (mpfr_sgn(tilted_mass.value) <= 0) {
        die("tilted interval mass lower bound is not positive");
    }

    Real log_one_minus_s_lower;
    Real log_mass_lower;
    mpfr_log(
        log_one_minus_s_lower.value,
        one_minus_s_lower.value,
        MPFR_RNDD
    );
    mpfr_log(log_mass_lower.value, tilted_mass.value, MPFR_RNDD);

    Real s_square_upper;
    mpfr_sqr(s_square_upper.value, geometry.s_upper.value, MPFR_RNDU);
    mpfr_div_ui(s_square_upper.value, s_square_upper.value, 2, MPFR_RNDU);
    Real s_h_upper;
    mpfr_mul(
        s_h_upper.value,
        geometry.s_upper.value,
        h_upper.value,
        MPFR_RNDU
    );
    Real exponent_upper;
    mpfr_add(
        exponent_upper.value,
        s_square_upper.value,
        s_h_upper.value,
        MPFR_RNDU
    );
    Real negative_exponent_lower;
    mpfr_neg(negative_exponent_lower.value, exponent_upper.value, MPFR_RNDD);

    mpfr_add(
        log_tail.value,
        log_one_minus_s_lower.value,
        log_mass_lower.value,
        MPFR_RNDD
    );
    mpfr_add(
        log_tail.value,
        log_tail.value,
        negative_exponent_lower.value,
        MPFR_RNDD
    );
}

void square_trace_tail_log_upper(
    Real& out,
    const Geometry& geometry,
    const Real& log4_upper
) {
    Real one;
    mpfr_set_ui(one.value, 1, MPFR_RNDU);
    Real a_minus_one_lower;
    mpfr_sub(
        a_minus_one_lower.value,
        geometry.a_lower.value,
        one.value,
        MPFR_RNDD
    );
    if (mpfr_sgn(a_minus_one_lower.value) <= 0) die("a must exceed one");
    Real square_lower;
    mpfr_sqr(square_lower.value, a_minus_one_lower.value, MPFR_RNDD);
    mpfr_div_ui(square_lower.value, square_lower.value, 4, MPFR_RNDD);
    mpfr_sub(out.value, log4_upper.value, square_lower.value, MPFR_RNDU);
}

// Rains' p=2 decomposition gives, for U Haar in Sp(2b),
//   -Tr(U^2)-1 = -(Y_0+Y_1),
// where Y_0 and Y_1 are independent centered full traces in
// O^-(2 ceil((b+1)/2)) and O^+(2 ceil(b/2)+1).  CJL's MGF upper
// sandwich and Chernoff at v=(a-1)/2 therefore give
//   P{-Tr(U^2)>=a}
//     <= exp(-(a-1)^2/4)/((1-beta_{d0}(v))(1-beta_{d1}(v))).
void c_square_trace_tail_log_upper(
    Real& out,
    Real& domain,
    int rank,
    const Geometry& geometry,
    const std::vector<Real>& factorial_logs
) {
    const int m0 = (rank + 2) / 2;
    const int m1 = (rank + 1) / 2;
    const int pair_count0 = m0 - 1;
    const int pair_count1 = m1;
    const int minimum_pair_count =
        pair_count0 < pair_count1 ? pair_count0 : pair_count1;
    const int dimension0 = 2 * m0;
    const int dimension1 = 2 * m1 + 1;

    Real one_upper;
    mpfr_set_ui(one_upper.value, 1, MPFR_RNDU);
    Real a_minus_one_upper;
    mpfr_sub(
        a_minus_one_upper.value,
        geometry.a_upper.value,
        one_upper.value,
        MPFR_RNDU
    );
    Real v_upper;
    mpfr_div_ui(v_upper.value, a_minus_one_upper.value, 2, MPFR_RNDU);
    domain_margin_lower(domain, minimum_pair_count, v_upper);

    Real beta0_upper;
    Real beta1_upper;
    beta_upper(beta0_upper, dimension0, v_upper, factorial_logs);
    beta_upper(beta1_upper, dimension1, v_upper, factorial_logs);
    Real one_minus_beta0_lower;
    Real one_minus_beta1_lower;
    one_minus_beta_lower(one_minus_beta0_lower, beta0_upper);
    one_minus_beta_lower(one_minus_beta1_lower, beta1_upper);

    Real one_lower;
    mpfr_set_ui(one_lower.value, 1, MPFR_RNDD);
    Real a_minus_one_lower;
    mpfr_sub(
        a_minus_one_lower.value,
        geometry.a_lower.value,
        one_lower.value,
        MPFR_RNDD
    );
    if (mpfr_sgn(a_minus_one_lower.value) <= 0) die("a must exceed one");
    Real gaussian_exponent_lower;
    mpfr_sqr(
        gaussian_exponent_lower.value,
        a_minus_one_lower.value,
        MPFR_RNDD
    );
    mpfr_div_ui(
        gaussian_exponent_lower.value,
        gaussian_exponent_lower.value,
        4,
        MPFR_RNDD
    );

    Real log_one_minus0_lower;
    Real log_one_minus1_lower;
    mpfr_log(
        log_one_minus0_lower.value,
        one_minus_beta0_lower.value,
        MPFR_RNDD
    );
    mpfr_log(
        log_one_minus1_lower.value,
        one_minus_beta1_lower.value,
        MPFR_RNDD
    );
    mpfr_neg(out.value, gaussian_exponent_lower.value, MPFR_RNDU);
    mpfr_sub(
        out.value,
        out.value,
        log_one_minus0_lower.value,
        MPFR_RNDU
    );
    mpfr_sub(
        out.value,
        out.value,
        log_one_minus1_lower.value,
        MPFR_RNDU
    );
}

void positive_push_log_lower(
    Real& out,
    int rank,
    const Real& tilted_tail_log_lower,
    const Real& log2_upper
) {
    // log((r b)^2 (1-(alpha sqrt(b))^(-2)) V_b/2).
    Real log_rb_lower;
    log_rational(
        log_rb_lower,
        R_NUM * static_cast<unsigned long>(rank),
        R_DEN,
        MPFR_RNDD
    );
    mpfr_mul_ui(log_rb_lower.value, log_rb_lower.value, 2, MPFR_RNDD);

    const unsigned long alpha_square_b_num =
        ALPHA_NUM * ALPHA_NUM * static_cast<unsigned long>(rank);
    const unsigned long alpha_square_den = ALPHA_DEN * ALPHA_DEN;
    if (alpha_square_b_num <= alpha_square_den) die("alpha^2 b must exceed one");
    Real log_variance_factor_lower;
    log_rational(
        log_variance_factor_lower,
        alpha_square_b_num - alpha_square_den,
        alpha_square_b_num,
        MPFR_RNDD
    );

    mpfr_add(
        out.value,
        log_rb_lower.value,
        log_variance_factor_lower.value,
        MPFR_RNDD
    );
    mpfr_add(
        out.value,
        out.value,
        tilted_tail_log_lower.value,
        MPFR_RNDD
    );
    mpfr_sub(out.value, out.value, log2_upper.value, MPFR_RNDD);
}

void negative_push_term_log_upper(
    Real& out,
    int rank,
    const Real& tau_tail_log_upper
) {
    const unsigned long b = static_cast<unsigned long>(rank);
    Real leading_upper;
    log_ui(leading_upper, 2 * (b * b + 1), MPFR_RNDU);

    Real ratio_log_upper;
    log_rational(ratio_log_upper, 2 * R_DEN, R_NUM, MPFR_RNDU);
    mpfr_mul_ui(
        ratio_log_upper.value,
        ratio_log_upper.value,
        static_cast<unsigned long>(2 * rank + 29),
        MPFR_RNDU
    );

    mpfr_add(
        out.value,
        leading_upper.value,
        tau_tail_log_upper.value,
        MPFR_RNDU
    );
    mpfr_add(out.value, out.value, ratio_log_upper.value, MPFR_RNDU);
}

void central_push_term_log_upper(
    Real& out,
    int rank,
    const Geometry& geometry,
    const Real& log2_upper
) {
    // eta/r = a/(r b), with a=alpha sqrt(b).
    Real rb_lower;
    set_rational(
        rb_lower,
        R_NUM * static_cast<unsigned long>(rank),
        R_DEN,
        MPFR_RNDD
    );
    Real ratio_upper;
    mpfr_div(
        ratio_upper.value,
        geometry.a_upper.value,
        rb_lower.value,
        MPFR_RNDU
    );
    Real one;
    mpfr_set_ui(one.value, 1, MPFR_RNDD);
    if (mpfr_cmp(ratio_upper.value, one.value) >= 0) die("eta/r must be below one");

    Real ratio_log_upper;
    mpfr_log(ratio_log_upper.value, ratio_upper.value, MPFR_RNDU);
    mpfr_mul_ui(
        ratio_log_upper.value,
        ratio_log_upper.value,
        static_cast<unsigned long>(2 * rank + 29),
        MPFR_RNDU
    );
    mpfr_add(out.value, log2_upper.value, ratio_log_upper.value, MPFR_RNDU);
}

RowMargins verify_b_row(
    int rank,
    const std::vector<Real>& factorial_logs,
    const Real& log2_upper,
    const Real& log4_upper
) {
    const Geometry geometry = geometry_for(rank);
    Real tilted_log_lower;
    Real tau_log_upper;
    Real push_log_lower;
    Real negative_log_upper;
    Real central_log_upper;

    RowMargins result;
    result.rank = rank;
    domain_margin_lower(result.domain, rank, geometry.plus_upper);
    mpfr_set_ui(result.square_domain.value, 1, MPFR_RNDD);
    tilted_tail_log_lower(
        tilted_log_lower,
        result.tilted_mass,
        2 * rank + 1,
        geometry,
        factorial_logs
    );
    square_trace_tail_log_upper(tau_log_upper, geometry, log4_upper);
    positive_push_log_lower(
        push_log_lower,
        rank,
        tilted_log_lower,
        log2_upper
    );
    negative_push_term_log_upper(negative_log_upper, rank, tau_log_upper);
    central_push_term_log_upper(
        central_log_upper,
        rank,
        geometry,
        log2_upper
    );

    // V_b >= 2 U_b.
    mpfr_sub(
        result.tilted_over_twice_tau.value,
        tilted_log_lower.value,
        log2_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.tilted_over_twice_tau.value,
        result.tilted_over_twice_tau.value,
        tau_log_upper.value,
        MPFR_RNDD
    );

    // The positive push lower bound dominates twice each RHS contribution.
    mpfr_sub(
        result.push_over_twice_negative.value,
        push_log_lower.value,
        log2_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.push_over_twice_negative.value,
        result.push_over_twice_negative.value,
        negative_log_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.push_over_twice_central.value,
        push_log_lower.value,
        log2_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.push_over_twice_central.value,
        result.push_over_twice_central.value,
        central_log_upper.value,
        MPFR_RNDD
    );
    return result;
}

RowMargins verify_c_row(
    int rank,
    const std::vector<Real>& factorial_logs,
    const Real& log2_upper
) {
    const Geometry geometry = geometry_for(rank);
    Real tilted_log_lower;
    Real tau_log_upper;
    Real push_log_lower;
    Real negative_log_upper;
    Real central_log_upper;

    RowMargins result;
    result.rank = rank;
    domain_margin_lower(result.domain, rank, geometry.plus_upper);
    tilted_tail_log_lower(
        tilted_log_lower,
        result.tilted_mass,
        2 * rank,
        geometry,
        factorial_logs
    );
    c_square_trace_tail_log_upper(
        tau_log_upper,
        result.square_domain,
        rank,
        geometry,
        factorial_logs
    );
    positive_push_log_lower(
        push_log_lower,
        rank,
        tilted_log_lower,
        log2_upper
    );
    negative_push_term_log_upper(negative_log_upper, rank, tau_log_upper);
    central_push_term_log_upper(
        central_log_upper,
        rank,
        geometry,
        log2_upper
    );

    mpfr_sub(
        result.tilted_over_twice_tau.value,
        tilted_log_lower.value,
        log2_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.tilted_over_twice_tau.value,
        result.tilted_over_twice_tau.value,
        tau_log_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.push_over_twice_negative.value,
        push_log_lower.value,
        log2_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.push_over_twice_negative.value,
        result.push_over_twice_negative.value,
        negative_log_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.push_over_twice_central.value,
        push_log_lower.value,
        log2_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.push_over_twice_central.value,
        result.push_over_twice_central.value,
        central_log_upper.value,
        MPFR_RNDD
    );
    return result;
}

bool positive(const Real& x) { return mpfr_sgn(x.value) > 0; }

}  // namespace

int main() {
    std::cout << "B/C tilted-MGF post-m29 trace certificate\n";
    std::cout << "mpfr_version=" << mpfr_get_version()
              << " precision_bits=" << PRECISION << '\n';
    std::cout << "source_mgf=CJL_2022_Proposition_3.3_and_Lemma_3.7\n";
    std::cout << "source_square_B="
                 "Rains_2003_Theorem_1.5_equation_23_p2"
                 "+CJL_2022_Lemma_2.2_m1\n";
    std::cout << "source_square_C="
                 "Rains_2003_Theorem_1.5_equation_22_p2_and_equation_36"
                 "+CJL_2022_Proposition_3.3_and_Lemma_3.7\n";
    std::cout << "r=" << R_NUM << '/' << R_DEN
              << " alpha=" << ALPHA_NUM << '/' << ALPHA_DEN
              << " h=" << H_NUM << '/' << H_DEN
              << " tail_start=2b+29 ranks=B_" << B_MIN << "..B_" << B_MAX
              << ",C_" << C_MIN << "..C_" << C_MAX
              << '\n';

    Real log2_upper;
    log_ui(log2_upper, 2, MPFR_RNDU);
    Real log4_upper;
    mpfr_mul_ui(log4_upper.value, log2_upper.value, 2, MPFR_RNDU);
    const auto factorial_logs = factorial_log_lowers(2 * C_MAX + 1);

    int failures = 0;
    std::vector<RowMargins> b_rows;
    b_rows.reserve(static_cast<std::size_t>(B_MAX - B_MIN + 1));
    for (int rank = B_MIN; rank <= B_MAX; ++rank) {
        RowMargins row = verify_b_row(
            rank,
            factorial_logs,
            log2_upper,
            log4_upper
        );
        const bool ok =
            positive(row.domain)
            && positive(row.tilted_mass)
            && positive(row.tilted_over_twice_tau)
            && positive(row.push_over_twice_negative)
            && positive(row.push_over_twice_central);
        if (!ok) ++failures;
        std::cout << "ROW B_" << rank
                  << " domain=" << format_real(row.domain)
                  << " tilted_mass=" << format_real(row.tilted_mass)
                  << " tilted_over_2tau="
                  << format_real(row.tilted_over_twice_tau)
                  << " push_over_2negative="
                  << format_real(row.push_over_twice_negative)
                  << " push_over_2central="
                  << format_real(row.push_over_twice_central)
                  << " ok=" << (ok ? 1 : 0) << '\n';
        b_rows.push_back(std::move(row));
    }

    auto print_minimum = [](
        const char* family,
        const std::vector<RowMargins>& rows,
        const char* name,
        auto member
    ) {
        const RowMargins* minimum = &rows.front();
        for (const RowMargins& row : rows) {
            if (mpfr_cmp((row.*member).value, (minimum->*member).value) < 0) {
                minimum = &row;
            }
        }
        std::cout << "MINIMUM " << name << ' ' << family << '_' << minimum->rank
                  << " value=" << format_real(minimum->*member) << '\n';
    };

    print_minimum(
        "B", b_rows, "domain", &RowMargins::domain
    );
    print_minimum(
        "B", b_rows, "tilted_mass", &RowMargins::tilted_mass
    );
    print_minimum(
        "B",
        b_rows,
        "tilted_over_2tau",
        &RowMargins::tilted_over_twice_tau
    );
    print_minimum(
        "B",
        b_rows,
        "push_over_2negative",
        &RowMargins::push_over_twice_negative
    );
    print_minimum(
        "B",
        b_rows,
        "push_over_2central",
        &RowMargins::push_over_twice_central
    );

    std::vector<RowMargins> c_rows;
    c_rows.reserve(static_cast<std::size_t>(C_MAX - C_MIN + 1));
    for (int rank = C_MIN; rank <= C_MAX; ++rank) {
        RowMargins row = verify_c_row(
            rank,
            factorial_logs,
            log2_upper
        );
        const bool ok =
            positive(row.domain)
            && positive(row.square_domain)
            && positive(row.tilted_mass)
            && positive(row.tilted_over_twice_tau)
            && positive(row.push_over_twice_negative)
            && positive(row.push_over_twice_central);
        if (!ok) ++failures;
        std::cout << "ROW C_" << rank
                  << " domain=" << format_real(row.domain)
                  << " square_domain=" << format_real(row.square_domain)
                  << " tilted_mass=" << format_real(row.tilted_mass)
                  << " tilted_over_2tau="
                  << format_real(row.tilted_over_twice_tau)
                  << " push_over_2negative="
                  << format_real(row.push_over_twice_negative)
                  << " push_over_2central="
                  << format_real(row.push_over_twice_central)
                  << " ok=" << (ok ? 1 : 0) << '\n';
        c_rows.push_back(std::move(row));
    }

    print_minimum(
        "C", c_rows, "domain", &RowMargins::domain
    );
    print_minimum(
        "C", c_rows, "square_domain", &RowMargins::square_domain
    );
    print_minimum(
        "C", c_rows, "tilted_mass", &RowMargins::tilted_mass
    );
    print_minimum(
        "C",
        c_rows,
        "tilted_over_2tau",
        &RowMargins::tilted_over_twice_tau
    );
    print_minimum(
        "C",
        c_rows,
        "push_over_2negative",
        &RowMargins::push_over_twice_negative
    );
    print_minimum(
        "C",
        c_rows,
        "push_over_2central",
        &RowMargins::push_over_twice_central
    );

    std::cout << "SUMMARY B_rows_checked=" << b_rows.size()
              << " C_rows_checked=" << c_rows.size()
              << " rows_checked=" << (b_rows.size() + c_rows.size())
              << " failures=" << failures << '\n';
    std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << '\n';
    return failures == 0 ? 0 : 1;
}
