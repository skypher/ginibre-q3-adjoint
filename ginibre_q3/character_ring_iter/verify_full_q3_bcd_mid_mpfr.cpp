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
int mpfr_sqrt(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_log(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_exp(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_lngamma(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_neg(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_cmp(mpfr_srcptr, mpfr_srcptr);
int mpfr_sgn(mpfr_srcptr);
int mpfr_asprintf(char**, const char*, ...);
void mpfr_free_str(char*);
const char* mpfr_get_version(void);
}
#endif

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "full_q3_bcd_low_tail_data.hpp"

namespace {

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
        if (this != &other) {
            mpfr_set(value, other.value, MPFR_RNDN);
        }
        return *this;
    }

    Real(Real&& other) noexcept {
        mpfr_init2(value, PRECISION);
        mpfr_swap(value, other.value);
    }

    Real& operator=(Real&& other) noexcept {
        if (this != &other) {
            mpfr_swap(value, other.value);
        }
        return *this;
    }

    ~Real() { mpfr_clear(value); }
};

struct Interval {
    Real lo;
    Real hi;
};

struct FredholmBound {
    Real denominator_lower;
    Real a_upper;
    Real b_upper;
};

struct Parameters {
    char family;
    int first_rank;
    int last_rank;
    unsigned long central_numerator;
    unsigned long central_denominator;
    unsigned long q_numerator;
    unsigned long q_denominator;
    unsigned long h_numerator;
    unsigned long h_denominator;
    bool q_times_square_root;
};

struct RowResult {
    Real margin_lower;
    Real rho_lower;
    Real trace_minus_square_log_lower;
    Real max_a_upper;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

Interval rational(unsigned long numerator, unsigned long denominator = 1) {
    Interval answer;
    mpfr_set_ui(answer.lo.value, numerator, MPFR_RNDD);
    mpfr_div_ui(answer.lo.value, answer.lo.value, denominator, MPFR_RNDD);
    mpfr_set_ui(answer.hi.value, numerator, MPFR_RNDU);
    mpfr_div_ui(answer.hi.value, answer.hi.value, denominator, MPFR_RNDU);
    return answer;
}

Interval add(const Interval& left, const Interval& right) {
    Interval answer;
    mpfr_add(answer.lo.value, left.lo.value, right.lo.value, MPFR_RNDD);
    mpfr_add(answer.hi.value, left.hi.value, right.hi.value, MPFR_RNDU);
    return answer;
}

Interval negate(const Interval& value) {
    Interval answer;
    mpfr_neg(answer.lo.value, value.hi.value, MPFR_RNDD);
    mpfr_neg(answer.hi.value, value.lo.value, MPFR_RNDU);
    return answer;
}

Interval subtract(const Interval& left, const Interval& right) {
    return add(left, negate(right));
}

void assign_min(Real& target, const Real& candidate) {
    if (mpfr_cmp(candidate.value, target.value) < 0) {
        mpfr_set(target.value, candidate.value, MPFR_RNDN);
    }
}

void assign_max(Real& target, const Real& candidate) {
    if (mpfr_cmp(candidate.value, target.value) > 0) {
        mpfr_set(target.value, candidate.value, MPFR_RNDN);
    }
}

Interval multiply(const Interval& left, const Interval& right) {
    std::array<Real, 4> lower;
    std::array<Real, 4> upper;
    const std::array<mpfr_srcptr, 2> left_values{left.lo.value, left.hi.value};
    const std::array<mpfr_srcptr, 2> right_values{right.lo.value, right.hi.value};
    std::size_t index = 0;
    for (mpfr_srcptr left_value : left_values) {
        for (mpfr_srcptr right_value : right_values) {
            mpfr_mul(lower[index].value, left_value, right_value, MPFR_RNDD);
            mpfr_mul(upper[index].value, left_value, right_value, MPFR_RNDU);
            ++index;
        }
    }
    Interval answer;
    answer.lo = lower[0];
    answer.hi = upper[0];
    for (std::size_t position = 1; position < lower.size(); ++position) {
        assign_min(answer.lo, lower[position]);
        assign_max(answer.hi, upper[position]);
    }
    return answer;
}

Interval reciprocal_positive(const Interval& value) {
    if (mpfr_sgn(value.lo.value) <= 0) fail("nonpositive interval divisor");
    Interval answer;
    const Interval one = rational(1);
    mpfr_div(answer.lo.value, one.lo.value, value.hi.value, MPFR_RNDD);
    mpfr_div(answer.hi.value, one.hi.value, value.lo.value, MPFR_RNDU);
    return answer;
}

Interval square_root(const Interval& value) {
    if (mpfr_sgn(value.lo.value) < 0) fail("square root of a negative interval");
    Interval answer;
    mpfr_sqrt(answer.lo.value, value.lo.value, MPFR_RNDD);
    mpfr_sqrt(answer.hi.value, value.hi.value, MPFR_RNDU);
    return answer;
}

Interval logarithm(const Interval& value) {
    if (mpfr_sgn(value.lo.value) <= 0) fail("logarithm of a nonpositive interval");
    Interval answer;
    mpfr_log(answer.lo.value, value.lo.value, MPFR_RNDD);
    mpfr_log(answer.hi.value, value.hi.value, MPFR_RNDU);
    return answer;
}

Interval exponential(const Interval& value) {
    Interval answer;
    mpfr_exp(answer.lo.value, value.lo.value, MPFR_RNDD);
    mpfr_exp(answer.hi.value, value.hi.value, MPFR_RNDU);
    return answer;
}

Interval scale(const Interval& value, unsigned long numerator, unsigned long denominator = 1) {
    return multiply(value, rational(numerator, denominator));
}

std::string format(const Real& value) {
    char* raw = nullptr;
    if (mpfr_asprintf(&raw, "%.32Rg", value.value) < 0 || raw == nullptr) {
        return "<format-error>";
    }
    std::string answer(raw);
    mpfr_free_str(raw);
    return answer;
}

Real log_factorial_lower(int dimension) {
    Real argument;
    mpfr_set_ui(argument.value, static_cast<unsigned long>(dimension + 1), MPFR_RNDD);
    Real answer;
    mpfr_lngamma(answer.value, argument.value, MPFR_RNDD);
    return answer;
}

FredholmBound fredholm_bound(int dimension, const Real& argument_upper) {
    if (dimension < 1 || mpfr_sgn(argument_upper.value) <= 0) {
        fail("invalid Fredholm dimension or argument");
    }
    const unsigned long d = static_cast<unsigned long>(dimension);

    Real d_plus_one;
    mpfr_set_ui(d_plus_one.value, d + 1, MPFR_RNDD);

    Real argument_square_upper;
    mpfr_mul(
        argument_square_upper.value,
        argument_upper.value,
        argument_upper.value,
        MPFR_RNDU
    );

    Real first_upper;
    mpfr_div(first_upper.value, argument_square_upper.value, d_plus_one.value, MPFR_RNDU);

    Real log_argument_upper;
    mpfr_log(log_argument_upper.value, argument_upper.value, MPFR_RNDU);
    mpfr_mul_ui(log_argument_upper.value, log_argument_upper.value, d, MPFR_RNDU);

    const Real factorial_lower = log_factorial_lower(dimension);

    Real linear_ratio_upper;
    mpfr_div(linear_ratio_upper.value, argument_upper.value, d_plus_one.value, MPFR_RNDU);
    Real one;
    mpfr_set_ui(one.value, 1, MPFR_RNDD);
    Real linear_denominator_lower;
    mpfr_sub(
        linear_denominator_lower.value,
        one.value,
        linear_ratio_upper.value,
        MPFR_RNDD
    );
    if (mpfr_sgn(linear_denominator_lower.value) <= 0) {
        fail("Fredholm linear denominator is nonpositive");
    }

    Real quadratic_ratio_upper;
    mpfr_div(
        quadratic_ratio_upper.value,
        argument_square_upper.value,
        d_plus_one.value,
        MPFR_RNDU
    );
    mpfr_div(
        quadratic_ratio_upper.value,
        quadratic_ratio_upper.value,
        d_plus_one.value,
        MPFR_RNDU
    );
    Real quadratic_denominator_lower;
    mpfr_sub(
        quadratic_denominator_lower.value,
        one.value,
        quadratic_ratio_upper.value,
        MPFR_RNDD
    );
    if (mpfr_sgn(quadratic_denominator_lower.value) <= 0) {
        fail("Fredholm quadratic denominator is nonpositive");
    }

    Real minus_log_linear_upper;
    mpfr_log(
        minus_log_linear_upper.value,
        linear_denominator_lower.value,
        MPFR_RNDD
    );
    mpfr_neg(minus_log_linear_upper.value, minus_log_linear_upper.value, MPFR_RNDU);

    Real minus_log_quadratic_upper;
    mpfr_log(
        minus_log_quadratic_upper.value,
        quadratic_denominator_lower.value,
        MPFR_RNDD
    );
    mpfr_neg(
        minus_log_quadratic_upper.value,
        minus_log_quadratic_upper.value,
        MPFR_RNDU
    );

    Real common_log_upper;
    mpfr_add(common_log_upper.value, first_upper.value, log_argument_upper.value, MPFR_RNDU);
    mpfr_sub(common_log_upper.value, common_log_upper.value, factorial_lower.value, MPFR_RNDU);

    Real log_b_upper;
    mpfr_add(
        log_b_upper.value,
        common_log_upper.value,
        minus_log_quadratic_upper.value,
        MPFR_RNDU
    );

    Real log_a_upper;
    mpfr_add(log_a_upper.value, common_log_upper.value, minus_log_linear_upper.value, MPFR_RNDU);
    Real half_minus_log_quadratic_upper;
    mpfr_div_ui(
        half_minus_log_quadratic_upper.value,
        minus_log_quadratic_upper.value,
        2,
        MPFR_RNDU
    );
    mpfr_add(log_a_upper.value, log_a_upper.value, half_minus_log_quadratic_upper.value, MPFR_RNDU);

    FredholmBound answer;
    mpfr_exp(answer.a_upper.value, log_a_upper.value, MPFR_RNDU);
    mpfr_exp(answer.b_upper.value, log_b_upper.value, MPFR_RNDU);
    if (mpfr_cmp(answer.a_upper.value, one.value) >= 0) {
        fail("Fredholm trace-norm envelope reached one");
    }
    mpfr_sub(
        answer.denominator_lower.value,
        one.value,
        answer.b_upper.value,
        MPFR_RNDD
    );
    if (mpfr_sgn(answer.denominator_lower.value) <= 0) {
        fail("Fredholm determinant denominator is nonpositive");
    }
    return answer;
}

void update_max_a(Real& maximum, const FredholmBound& bound) {
    if (mpfr_cmp(bound.a_upper.value, maximum.value) > 0) {
        mpfr_set(maximum.value, bound.a_upper.value, MPFR_RNDN);
    }
}

std::pair<int, int> square_dimensions(char family, int rank) {
    if (family == 'C') {
        const int m0 = (rank + 2) / 2;
        const int m1 = (rank + 1) / 2;
        return {2 * m0, 2 * m1 + 1};
    }
    if (family == 'D') {
        const int m0 = (rank + 1) / 2;
        const int m1 = rank / 2;
        return {2 * m0, 2 * m1 + 1};
    }
    return {0, 0};
}

RowResult evaluate_row(
    const Parameters& parameters,
    int rank,
    int total_degree_override = 0
) {
    const char family = parameters.family;
    const int endpoint =
        family != 'D' || rank % 2 == 0 ? rank : rank - 2;
    const int stable_tail_degree =
        family == 'B' || family == 'C'
            ? 2 * rank + 3
            : (rank % 2 == 0 ? rank + 1 : rank);
    const int total_degree =
        total_degree_override > 0 ? total_degree_override : stable_tail_degree;

    const Interval rank_value = rational(static_cast<unsigned long>(rank));
    const Interval endpoint_value = rational(static_cast<unsigned long>(endpoint));
    const Interval total_value = rational(static_cast<unsigned long>(total_degree));
    const Interval r = rational(2000001, 1000000);
    const Interval central_cutoff = rational(
        parameters.central_numerator,
        parameters.central_denominator
    );
    const Interval q_coefficient = rational(
        parameters.q_numerator,
        parameters.q_denominator
    );
    const Interval h = rational(parameters.h_numerator, parameters.h_denominator);

    const Interval c = multiply(r, endpoint_value);
    const Interval q_scale =
        parameters.q_times_square_root ? square_root(rank_value) : rank_value;
    const Interval q = multiply(q_coefficient, q_scale);
    Interval threshold_square = add(scale(c, 2), scale(central_cutoff, 2));
    threshold_square = add(threshold_square, q);
    const Interval t = square_root(threshold_square);
    // The full defining character is centered in all three families.  In
    // type B the fixed +1 eigenvalue cancels the -1 mean of the random
    // eigenvalue contribution in the Fredholm formula.
    const Interval s = add(t, h);
    const Interval s_minus_h = subtract(s, h);
    const Interval s_plus_h = add(s, h);
    const Interval event_upper = add(t, scale(h, 2));
    if (mpfr_sgn(s_minus_h.lo.value) <= 0) fail("nonpositive trace tilt argument");

    const int trace_dimension = family == 'B' ? 2 * rank + 1 : 2 * rank;
    FredholmBound trace_t = fredholm_bound(trace_dimension, s_minus_h.hi);
    FredholmBound trace_s = fredholm_bound(trace_dimension, s.hi);
    FredholmBound trace_high = fredholm_bound(trace_dimension, s_plus_h.hi);

    Real max_a;
    update_max_a(max_a, trace_t);
    update_max_a(max_a, trace_s);
    update_max_a(max_a, trace_high);

    // With s=t+h, both tilted exit losses equal exp(-h^2/2).
    const Interval h_square = multiply(h, h);
    const Interval minus_half_h_square = negate(scale(h_square, 1, 2));

    Interval d_s;
    d_s.lo = trace_s.denominator_lower;
    d_s.hi = rational(1).hi;
    Interval d_t;
    d_t.lo = trace_t.denominator_lower;
    d_t.hi = rational(1).hi;
    Interval d_high;
    d_high.lo = trace_high.denominator_lower;
    d_high.hi = rational(1).hi;

    const Interval reciprocal_sum = add(reciprocal_positive(d_t), reciprocal_positive(d_high));
    const Interval rho_loss = multiply(
        multiply(exponential(minus_half_h_square), reciprocal_positive(d_s)),
        reciprocal_sum
    );
    const Interval rho = subtract(rational(1), rho_loss);
    if (mpfr_sgn(rho.lo.value) <= 0) {
        fail(
            std::string("signed tilted trace mass is nonpositive at ") + family + "_" +
            std::to_string(rank)
        );
    }

    // The two one-sided Fredholm estimates have the same numerical lower
    // bound even when the trace law is not symmetric, so their disjoint
    // union supplies a factor two in every family.
    Interval log_p1 = logarithm(rational(2));
    log_p1 = add(log_p1, logarithm(d_s));
    log_p1 = add(log_p1, logarithm(rho));
    log_p1 = add(log_p1, scale(multiply(s, s), 1, 2));
    log_p1 = subtract(log_p1, multiply(s, event_upper));

    const Interval q_minus_one = subtract(q, rational(1));
    if (mpfr_sgn(q_minus_one.lo.value) <= 0) fail("square-trace cutoff is not above one");
    const Interval vq = scale(q_minus_one, 1, 2);
    Interval log_square_tail = negate(scale(multiply(q_minus_one, q_minus_one), 1, 4));

    const auto [square_d0, square_d1] = square_dimensions(family, rank);
    if (family != 'B') {
        FredholmBound square_q0 = fredholm_bound(square_d0, vq.hi);
        FredholmBound square_q1 = fredholm_bound(square_d1, vq.hi);
        update_max_a(max_a, square_q0);
        update_max_a(max_a, square_q1);

        Interval square_denominator0;
        square_denominator0.lo = square_q0.denominator_lower;
        square_denominator0.hi = rational(1).hi;
        Interval square_denominator1;
        square_denominator1.lo = square_q1.denominator_lower;
        square_denominator1.hi = rational(1).hi;
        log_square_tail = subtract(log_square_tail, logarithm(square_denominator0));
        log_square_tail = subtract(log_square_tail, logarithm(square_denominator1));
    }

    // A rigorous lower bound for log(p1-U).
    Real p1_lower;
    mpfr_exp(p1_lower.value, log_p1.lo.value, MPFR_RNDD);
    Real square_tail_upper;
    mpfr_exp(square_tail_upper.value, log_square_tail.hi.value, MPFR_RNDU);
    Real trace_minus_square_lower;
    mpfr_sub(
        trace_minus_square_lower.value,
        p1_lower.value,
        square_tail_upper.value,
        MPFR_RNDD
    );
    if (mpfr_sgn(trace_minus_square_lower.value) <= 0) {
        fail("defining-trace lower bound did not exceed square-trace tail");
    }
    Real log_trace_minus_square_lower;
    mpfr_log(
        log_trace_minus_square_lower.value,
        trace_minus_square_lower.value,
        MPFR_RNDD
    );

    const Interval central_probability = subtract(
        rational(1),
        reciprocal_positive(multiply(central_cutoff, central_cutoff))
    );
    Real log_p_lower;
    const Interval log_central = logarithm(central_probability);
    mpfr_add(
        log_p_lower.value,
        log_central.lo.value,
        log_trace_minus_square_lower.value,
        MPFR_RNDD
    );

    // v_m=(-1+sqrt(1+8L))/4.
    const Interval discriminant = add(rational(1), scale(total_value, 8));
    const Interval vm = scale(subtract(square_root(discriminant), rational(1)), 1, 4);
    if (mpfr_sgn(vm.lo.value) <= 0) fail("moment tilt is nonpositive");

    Interval moment_log = subtract(logarithm(total_value), rational(1));
    moment_log = subtract(moment_log, logarithm(vm));
    moment_log = scale(moment_log, static_cast<unsigned long>(total_degree));
    moment_log = add(moment_log, vm);
    moment_log = add(moment_log, multiply(vm, vm));

    if (family != 'B') {
        FredholmBound square_m0 = fredholm_bound(square_d0, vm.hi);
        FredholmBound square_m1 = fredholm_bound(square_d1, vm.hi);
        update_max_a(max_a, square_m0);
        update_max_a(max_a, square_m1);

        Interval square_denominator0;
        square_denominator0.lo = square_m0.denominator_lower;
        square_denominator0.hi = rational(1).hi;
        Interval square_denominator1;
        square_denominator1.lo = square_m1.denominator_lower;
        square_denominator1.hi = rational(1).hi;
        moment_log = subtract(moment_log, logarithm(square_denominator0));
        moment_log = subtract(moment_log, logarithm(square_denominator1));
    }

    const Interval log_c = logarithm(c);
    Real positive_log_lower;
    mpfr_mul_ui(
        positive_log_lower.value,
        log_c.lo.value,
        static_cast<unsigned long>(total_degree),
        MPFR_RNDD
    );
    mpfr_add(
        positive_log_lower.value,
        positive_log_lower.value,
        log_p_lower.value,
        MPFR_RNDD
    );

    RowResult answer;
    mpfr_sub(
        answer.margin_lower.value,
        positive_log_lower.value,
        moment_log.hi.value,
        MPFR_RNDD
    );
    if (mpfr_sgn(answer.margin_lower.value) <= 0) {
        fail("full-Q3 trace-tail margin is nonpositive");
    }
    answer.rho_lower = rho.lo;
    mpfr_sub(
        answer.trace_minus_square_log_lower.value,
        log_p1.lo.value,
        log_square_tail.hi.value,
        MPFR_RNDD
    );
    answer.max_a_upper = max_a;
    return answer;
}

}  // namespace

int main() {
    try {
        const std::array<Parameters, 3> parameter_sets{{
            {'B', 22, 295, 3, 2, 21, 5, 13, 10, true},
            {'C', 29, 295, 3, 2, 4, 1, 13, 10, true},
            {'D', 71, 297, 3, 2, 18, 5, 6, 5, true},
        }};

        bool have_minimum = false;
        Real minimum_margin;
        Real minimum_rho;
        Real minimum_trace_separation;
        Real maximum_a;
        std::string minimum_margin_label;
        std::string minimum_rho_label;
        std::string minimum_trace_label;
        std::string maximum_a_label;
        int checked = 0;

        for (const Parameters& parameters : parameter_sets) {
            for (int rank = parameters.first_rank; rank <= parameters.last_rank; ++rank) {
                RowResult result = evaluate_row(parameters, rank);
                const std::string label = std::string(1, parameters.family) + "_" +
                                          std::to_string(rank);
                if (!have_minimum ||
                    mpfr_cmp(result.margin_lower.value, minimum_margin.value) < 0) {
                    minimum_margin = result.margin_lower;
                    minimum_margin_label = label;
                }
                if (!have_minimum || mpfr_cmp(result.rho_lower.value, minimum_rho.value) < 0) {
                    minimum_rho = result.rho_lower;
                    minimum_rho_label = label;
                }
                if (!have_minimum ||
                    mpfr_cmp(
                        result.trace_minus_square_log_lower.value,
                        minimum_trace_separation.value
                    ) < 0) {
                    minimum_trace_separation = result.trace_minus_square_log_lower;
                    minimum_trace_label = label;
                }
                if (!have_minimum || mpfr_cmp(result.max_a_upper.value, maximum_a.value) > 0) {
                    maximum_a = result.max_a_upper;
                    maximum_a_label = label;
                }
                have_minimum = true;
                ++checked;
            }
        }

        const int b_rows = 295 - 22 + 1;
        const int c_rows = 295 - 29 + 1;
        const int d_rows = 297 - 71 + 1;
        const int expected = b_rows + c_rows + d_rows;
        if (!have_minimum || checked != expected) fail("finite rank ledger mismatch");

        const auto& d_low_onsets = full_q3_bcd_low_tail::d_onsets;
        const auto& b_low_onsets = full_q3_bcd_low_tail::b_onsets;

        bool have_low_minimum = false;
        Real low_minimum_margin;
        Real low_minimum_rho;
        Real low_minimum_trace_separation;
        Real low_maximum_a;
        std::string low_minimum_margin_label;
        std::string low_minimum_rho_label;
        std::string low_minimum_trace_label;
        std::string low_maximum_a_label;
        int low_checked = 0;

        const auto record_low = [&](char family, int rank, int onset, const RowResult& result) {
            const std::string label = std::string(1, family) + "_" +
                                      std::to_string(rank) + "@" +
                                      std::to_string(onset);
            if (!have_low_minimum ||
                mpfr_cmp(result.margin_lower.value, low_minimum_margin.value) < 0) {
                low_minimum_margin = result.margin_lower;
                low_minimum_margin_label = label;
            }
            if (!have_low_minimum ||
                mpfr_cmp(result.rho_lower.value, low_minimum_rho.value) < 0) {
                low_minimum_rho = result.rho_lower;
                low_minimum_rho_label = label;
            }
            if (!have_low_minimum ||
                mpfr_cmp(
                    result.trace_minus_square_log_lower.value,
                    low_minimum_trace_separation.value
                ) < 0) {
                low_minimum_trace_separation = result.trace_minus_square_log_lower;
                low_minimum_trace_label = label;
            }
            if (!have_low_minimum ||
                mpfr_cmp(result.max_a_upper.value, low_maximum_a.value) > 0) {
                low_maximum_a = result.max_a_upper;
                low_maximum_a_label = label;
            }
            have_low_minimum = true;
            ++low_checked;
        };

        for (int rank = 18; rank <= 21; ++rank) {
            const bool first = rank == 18;
            const Parameters parameters{
                'B', rank, rank,
                first ? 13UL : 3UL,
                first ? 10UL : 2UL,
                first ? 6UL : 1UL,
                first ? 5UL : 1UL,
                first ? 5UL : 6UL,
                first ? 4UL : 5UL,
                false,
            };
            const int onset = b_low_onsets[static_cast<std::size_t>(
                rank - full_q3_bcd_low_tail::b_first_rank
            )];
            record_low('B', rank, onset, evaluate_row(parameters, rank, onset));
        }

        for (int rank = 31; rank <= 70; ++rank) {
            unsigned long q_numerator = 1;
            unsigned long q_denominator = 2;
            if (rank <= 38) {
                q_numerator = 7;
                q_denominator = 10;
            } else if (rank <= 52) {
                q_numerator = 3;
                q_denominator = 5;
            }
            const Parameters parameters{
                'D', rank, rank, 3, 2,
                q_numerator, q_denominator, 6, 5, false,
            };
            const int onset = d_low_onsets[static_cast<std::size_t>(
                rank - full_q3_bcd_low_tail::d_first_rank
            )];
            record_low('D', rank, onset, evaluate_row(parameters, rank, onset));
        }

        if (!have_low_minimum || low_checked != 44) fail("low-tail ledger mismatch");

        std::cout << "BCD full-Q3 finite middle-rank directed-MPFR certificate\n";
        std::cout << "precision_bits=" << PRECISION
                  << " mpfr_version=" << mpfr_get_version() << '\n';
        std::cout << "rows_checked=" << checked
                  << " B_rows=" << b_rows
                  << " C_rows=" << c_rows
                  << " D_rows=" << d_rows << '\n';
        std::cout << "minimum_margin=" << minimum_margin_label << " "
                  << format(minimum_margin) << '\n';
        std::cout << "minimum_rho=" << minimum_rho_label << " "
                  << format(minimum_rho) << '\n';
        std::cout << "minimum_log_trace_over_square=" << minimum_trace_label << " "
                  << format(minimum_trace_separation) << '\n';
        std::cout << "maximum_fredholm_a=" << maximum_a_label << " "
                  << format(maximum_a) << '\n';
        std::cout << "low_tail_rows_checked=" << low_checked
                  << " B_low_rows=4 D_low_rows=40\n";
        std::cout << "low_tail_minimum_margin=" << low_minimum_margin_label << " "
                  << format(low_minimum_margin) << '\n';
        std::cout << "low_tail_minimum_rho=" << low_minimum_rho_label << " "
                  << format(low_minimum_rho) << '\n';
        std::cout << "low_tail_minimum_log_trace_over_square="
                  << low_minimum_trace_label << " "
                  << format(low_minimum_trace_separation) << '\n';
        std::cout << "low_tail_maximum_fredholm_a=" << low_maximum_a_label << " "
                  << format(low_maximum_a) << '\n';
        std::cout << "BCD_MID_FULL_Q3 VERIFICATION: ALL PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "BCD finite middle-rank full-Q3 verification failure: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
