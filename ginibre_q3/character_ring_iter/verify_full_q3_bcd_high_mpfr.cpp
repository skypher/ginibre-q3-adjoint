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
int mpfr_log(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_exp(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_sqrt(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
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
#include <stdexcept>
#include <string>

namespace {

constexpr mpfr_prec_t PRECISION = 384;
constexpr unsigned long CUTOFF = 296;

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

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

Real rational(
    unsigned long numerator,
    unsigned long denominator,
    mpfr_rnd_t rounding
) {
    Real answer;
    mpfr_set_ui(answer.value, numerator, rounding);
    mpfr_div_ui(answer.value, answer.value, denominator, rounding);
    return answer;
}

Real integer(unsigned long value) {
    Real answer;
    mpfr_set_ui(answer.value, value, MPFR_RNDN);
    return answer;
}

std::string format(const Real& value) {
    char* raw = nullptr;
    if (mpfr_asprintf(&raw, "%.36Rg", value.value) < 0 || raw == nullptr) {
        return "<format-error>";
    }
    std::string answer(raw);
    mpfr_free_str(raw);
    return answer;
}

Real exponential_one(mpfr_rnd_t rounding) {
    Real one = integer(1);
    Real answer;
    mpfr_exp(answer.value, one.value, rounding);
    return answer;
}

Real cutoff_base_lower() {
    // 4 e r C / (2 C + 3), with r=2000001/1000000.
    Real answer = exponential_one(MPFR_RNDD);
    mpfr_mul_ui(answer.value, answer.value, 4, MPFR_RNDD);
    const Real r = rational(2000001, 1000000, MPFR_RNDD);
    mpfr_mul(answer.value, answer.value, r.value, MPFR_RNDD);
    mpfr_mul_ui(answer.value, answer.value, CUTOFF, MPFR_RNDD);
    mpfr_div_ui(answer.value, answer.value, 2 * CUTOFF + 3, MPFR_RNDD);
    return answer;
}

Real high_rank_margin_lower() {
    // H(C)=log(3/32)-A*C-20+(C+1)log(base), A=4571/2000.
    const Real three_over_32 = rational(3, 32, MPFR_RNDD);
    Real answer;
    mpfr_log(answer.value, three_over_32.value, MPFR_RNDD);

    const Real ac_upper = rational(4571 * CUTOFF, 2000, MPFR_RNDU);
    mpfr_sub(answer.value, answer.value, ac_upper.value, MPFR_RNDD);
    const Real twenty = integer(20);
    mpfr_sub(answer.value, answer.value, twenty.value, MPFR_RNDD);

    const Real base = cutoff_base_lower();
    Real log_base;
    mpfr_log(log_base.value, base.value, MPFR_RNDD);
    mpfr_mul_ui(log_base.value, log_base.value, CUTOFF + 1, MPFR_RNDD);
    mpfr_add(answer.value, answer.value, log_base.value, MPFR_RNDD);
    return answer;
}

Real high_rank_derivative_lower() {
    // H'(C)=-A+log(base)+3(C+1)/(C(2C+3)).
    const Real a_upper = rational(4571, 2000, MPFR_RNDU);
    Real answer;
    mpfr_set_zero(answer.value, 0);
    mpfr_sub(answer.value, answer.value, a_upper.value, MPFR_RNDD);

    const Real base = cutoff_base_lower();
    Real log_base;
    mpfr_log(log_base.value, base.value, MPFR_RNDD);
    mpfr_add(answer.value, answer.value, log_base.value, MPFR_RNDD);

    const unsigned long numerator = 3 * (CUTOFF + 1);
    const unsigned long denominator = CUTOFF * (2 * CUTOFF + 3);
    const Real correction = rational(numerator, denominator, MPFR_RNDD);
    mpfr_add(answer.value, answer.value, correction.value, MPFR_RNDD);
    return answer;
}

Real high_rank_log_component_lower() {
    // log(4 e r C / (2 C + 3)) - A at the cutoff.  The logarithmic
    // factor is increasing in C, so positivity here is the fail-closed
    // monotonicity input used in the paper for every C >= CUTOFF.
    const Real a_upper = rational(4571, 2000, MPFR_RNDU);
    const Real base = cutoff_base_lower();
    Real answer;
    mpfr_log(answer.value, base.value, MPFR_RNDD);
    mpfr_sub(answer.value, answer.value, a_upper.value, MPFR_RNDD);
    return answer;
}

Real power_four_11(mpfr_rnd_t rounding) {
    Real answer = integer(1);
    const Real four = integer(4);
    for (int index = 0; index < 11; ++index) {
        mpfr_mul(answer.value, answer.value, four.value, rounding);
    }
    return answer;
}

Real exp_four_thirds_upper() {
    const Real exponent = rational(4, 3, MPFR_RNDU);
    Real answer;
    mpfr_exp(answer.value, exponent.value, MPFR_RNDU);
    return answer;
}

Real fredholm_b11_upper() {
    // b_11(4)=exp(4/3)/(1-1/9) * 4^11/11!.
    Real answer = exp_four_thirds_upper();
    const Real power = power_four_11(MPFR_RNDU);
    mpfr_mul(answer.value, answer.value, power.value, MPFR_RNDU);
    mpfr_div_ui(answer.value, answer.value, 39916800, MPFR_RNDU);
    const Real denominator = rational(8, 9, MPFR_RNDD);
    mpfr_div(answer.value, answer.value, denominator.value, MPFR_RNDU);
    return answer;
}

Real fredholm_a11_upper() {
    // a_11(4)=exp(4/3)4^11/11! /[(1-1/3)sqrt(1-1/9)].
    Real answer = exp_four_thirds_upper();
    const Real power = power_four_11(MPFR_RNDU);
    mpfr_mul(answer.value, answer.value, power.value, MPFR_RNDU);
    mpfr_div_ui(answer.value, answer.value, 39916800, MPFR_RNDU);

    const Real linear_denominator = rational(2, 3, MPFR_RNDD);
    mpfr_div(answer.value, answer.value, linear_denominator.value, MPFR_RNDU);

    const Real square_argument = rational(8, 9, MPFR_RNDD);
    Real square_root_lower;
    mpfr_sqrt(square_root_lower.value, square_argument.value, MPFR_RNDD);
    mpfr_div(answer.value, answer.value, square_root_lower.value, MPFR_RNDU);
    return answer;
}

void require_greater(const Real& value, unsigned long numerator, unsigned long denominator) {
    const Real threshold = rational(numerator, denominator, MPFR_RNDU);
    if (mpfr_cmp(value.value, threshold.value) <= 0) {
        fail("directed lower bound did not exceed its threshold");
    }
}

void require_less(const Real& value, unsigned long numerator, unsigned long denominator) {
    const Real threshold = rational(numerator, denominator, MPFR_RNDD);
    if (mpfr_cmp(value.value, threshold.value) >= 0) {
        fail("directed upper bound did not lie below its threshold");
    }
}

}  // namespace

int main() {
    try {
        const Real base = cutoff_base_lower();
        const Real margin = high_rank_margin_lower();
        const Real derivative = high_rank_derivative_lower();
        const Real log_component = high_rank_log_component_lower();
        const Real fredholm_a = fredholm_a11_upper();
        const Real fredholm_b = fredholm_b11_upper();

        require_greater(base, 1, 1);
        require_greater(margin, 8, 1);
        require_greater(derivative, 1, 10);
        require_greater(log_component, 9, 100);
        require_less(fredholm_a, 1, 1);
        require_less(fredholm_b, 1, 2);

        std::cout << "BCD full-Q3 high-rank directed-MPFR certificate\n";
        std::cout << "precision_bits=" << PRECISION
                  << " mpfr_version=" << mpfr_get_version()
                  << " cutoff=" << CUTOFF << '\n';
        std::cout << "base_lower=" << format(base) << '\n';
        std::cout << "H_cutoff_lower=" << format(margin) << '\n';
        std::cout << "H_derivative_cutoff_lower=" << format(derivative) << '\n';
        std::cout << "log_base_minus_A_cutoff_lower="
                  << format(log_component) << '\n';
        std::cout << "fredholm_a11_v4_upper=" << format(fredholm_a) << '\n';
        std::cout << "fredholm_b11_v4_upper=" << format(fredholm_b) << '\n';
        std::cout << "BCD_HIGH_FULL_Q3 VERIFICATION: ALL PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "BCD high-rank full-Q3 verification failure: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
