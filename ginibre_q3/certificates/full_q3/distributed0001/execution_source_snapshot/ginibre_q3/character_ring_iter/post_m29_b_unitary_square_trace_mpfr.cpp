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
int mpfr_mul_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_div_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_sqr(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_log(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_neg(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_cmp(mpfr_srcptr, mpfr_srcptr);
int mpfr_sgn(mpfr_srcptr);
int mpfr_asprintf(char**, const char*, ...);
void mpfr_free_str(char*);
const char* mpfr_get_version(void);
}
#endif

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr int B_MIN = 124;
constexpr int B_MAX = 217;
constexpr unsigned long R_NUM = 2001;
constexpr unsigned long R_DEN = 1000;
constexpr unsigned long ETA_NUM = 7;
constexpr unsigned long ETA_DEN = 20;
constexpr unsigned long L_NUM = 1263;  // 2r + 3 eta = 1263/250.
constexpr unsigned long L_DEN = 250;
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
    Real gaussian_over_tv;
    Real gaussian_over_tau;
    Real push_over_twice_negative_tail;
    Real push_over_twice_central_tail;
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

// Lower bound for the two-sided standard-Gaussian tail used in the paper:
// G(x) = exp(-x/2) sqrt(x)/(sqrt(2)(x+1)), x=(2r+3eta)b.
void gaussian_log_lower(Real& out, int rank, const Real& log2_upper) {
    const unsigned long x_num =
        L_NUM * static_cast<unsigned long>(rank);

    Real x_upper;
    set_rational(x_upper, x_num, L_DEN, MPFR_RNDU);
    Real negative_half_x;
    mpfr_div_ui(negative_half_x.value, x_upper.value, 2, MPFR_RNDU);
    mpfr_neg(negative_half_x.value, negative_half_x.value, MPFR_RNDD);

    Real log_x_lower;
    log_rational(log_x_lower, x_num, L_DEN, MPFR_RNDD);
    mpfr_div_ui(log_x_lower.value, log_x_lower.value, 2, MPFR_RNDD);

    Real x_plus_one_upper;
    set_rational(
        x_plus_one_upper,
        x_num + L_DEN,
        L_DEN,
        MPFR_RNDU
    );
    Real log_x_plus_one_upper;
    mpfr_log(
        log_x_plus_one_upper.value,
        x_plus_one_upper.value,
        MPFR_RNDU
    );

    Real half_log2_upper;
    mpfr_div_ui(half_log2_upper.value, log2_upper.value, 2, MPFR_RNDU);

    mpfr_set(out.value, negative_half_x.value, MPFR_RNDD);
    mpfr_add(out.value, out.value, log_x_lower.value, MPFR_RNDD);
    mpfr_sub(out.value, out.value, log_x_plus_one_upper.value, MPFR_RNDD);
    mpfr_sub(out.value, out.value, half_log2_upper.value, MPFR_RNDD);
}

// Upper bound from CJL Theorem 1.4, using a downward-rounded exact
// logarithmic sum for (2b)!:
// E_b <= 5 (log(2b))^(1/4) / sqrt((2b)!).
void one_trace_error_log_upper(
    Real& out,
    int rank,
    const std::vector<Real>& factorial_logs
) {
    Real log5_upper;
    log_ui(log5_upper, 5, MPFR_RNDU);

    Real log_2b_upper;
    log_ui(log_2b_upper, static_cast<unsigned long>(2 * rank), MPFR_RNDU);
    Real log_log_2b_upper;
    mpfr_log(log_log_2b_upper.value, log_2b_upper.value, MPFR_RNDU);
    mpfr_div_ui(log_log_2b_upper.value, log_log_2b_upper.value, 4, MPFR_RNDU);

    Real half_factorial_log_lower;
    mpfr_div_ui(
        half_factorial_log_lower.value,
        factorial_logs[static_cast<std::size_t>(2 * rank)].value,
        2,
        MPFR_RNDD
    );

    mpfr_add(out.value, log5_upper.value, log_log_2b_upper.value, MPFR_RNDU);
    mpfr_sub(out.value, out.value, half_factorial_log_lower.value, MPFR_RNDU);
}

// Rains (23) at p=2 gives tau_B-1 = 2 Re Tr U(b).  CJL Lemma 2.2
// therefore gives U_b = 4 exp(-(eta*b-1)^2/4).
void unitary_square_tail_log_upper(
    Real& out,
    int rank,
    const Real& log4_upper
) {
    const unsigned long t_num =
        ETA_NUM * static_cast<unsigned long>(rank) - ETA_DEN;
    Real t_lower;
    set_rational(t_lower, t_num, ETA_DEN, MPFR_RNDD);
    Real t_square_lower;
    mpfr_sqr(t_square_lower.value, t_lower.value, MPFR_RNDD);
    mpfr_div_ui(t_square_lower.value, t_square_lower.value, 4, MPFR_RNDD);

    mpfr_sub(out.value, log4_upper.value, t_square_lower.value, MPFR_RNDU);
}

void positive_push_log_lower(
    Real& out,
    int rank,
    const Real& gaussian_lower,
    const Real& log2_upper
) {
    // log((r b)^2 (1-(eta b)^(-2)) G_b/2).
    Real log_rb_lower;
    log_rational(
        log_rb_lower,
        R_NUM * static_cast<unsigned long>(rank),
        R_DEN,
        MPFR_RNDD
    );
    mpfr_mul_ui(log_rb_lower.value, log_rb_lower.value, 2, MPFR_RNDD);

    const unsigned long eta_sq_b_sq_num =
        ETA_NUM * ETA_NUM * static_cast<unsigned long>(rank)
        * static_cast<unsigned long>(rank);
    const unsigned long eta_sq_den = ETA_DEN * ETA_DEN;
    if (eta_sq_b_sq_num <= eta_sq_den) die("eta*b must exceed one");
    Real log_den_lower;
    log_rational(
        log_den_lower,
        eta_sq_b_sq_num - eta_sq_den,
        eta_sq_b_sq_num,
        MPFR_RNDD
    );

    mpfr_add(out.value, log_rb_lower.value, log_den_lower.value, MPFR_RNDD);
    mpfr_add(out.value, out.value, gaussian_lower.value, MPFR_RNDD);
    mpfr_sub(out.value, out.value, log2_upper.value, MPFR_RNDD);
}

void negative_push_term_log_upper(
    Real& out,
    int rank,
    const Real& unitary_tail_upper
) {
    // log(2(b^2+1) U_b (2/r)^(2b+1)).
    Real leading_upper;
    const unsigned long b = static_cast<unsigned long>(rank);
    log_ui(leading_upper, 2 * (b * b + 1), MPFR_RNDU);

    Real ratio_log_upper;
    log_rational(ratio_log_upper, 2 * R_DEN, R_NUM, MPFR_RNDU);
    mpfr_mul_ui(
        ratio_log_upper.value,
        ratio_log_upper.value,
        static_cast<unsigned long>(2 * rank + 1),
        MPFR_RNDU
    );

    mpfr_add(out.value, leading_upper.value, unitary_tail_upper.value, MPFR_RNDU);
    mpfr_add(out.value, out.value, ratio_log_upper.value, MPFR_RNDU);
}

void central_push_term_log_upper(
    Real& out,
    int rank,
    const Real& log2_upper
) {
    // log(2 (eta/r)^(2b+1)).
    Real ratio_log_upper;
    log_rational(
        ratio_log_upper,
        ETA_NUM * R_DEN,
        ETA_DEN * R_NUM,
        MPFR_RNDU
    );
    mpfr_mul_ui(
        ratio_log_upper.value,
        ratio_log_upper.value,
        static_cast<unsigned long>(2 * rank + 1),
        MPFR_RNDU
    );
    mpfr_add(out.value, log2_upper.value, ratio_log_upper.value, MPFR_RNDU);
}

RowMargins verify_row(
    int rank,
    const std::vector<Real>& factorial_logs,
    const Real& log2_upper,
    const Real& log4_upper
) {
    Real gaussian_lower;
    Real error_upper;
    Real tau_upper;
    Real push_lower;
    Real negative_term_upper;
    Real central_term_upper;

    gaussian_log_lower(gaussian_lower, rank, log2_upper);
    one_trace_error_log_upper(error_upper, rank, factorial_logs);
    unitary_square_tail_log_upper(tau_upper, rank, log4_upper);
    positive_push_log_lower(push_lower, rank, gaussian_lower, log2_upper);
    negative_push_term_log_upper(negative_term_upper, rank, tau_upper);
    central_push_term_log_upper(central_term_upper, rank, log2_upper);

    RowMargins result;
    result.rank = rank;

    // E_b <= G_b/4 and U_b <= G_b/4.
    mpfr_sub(
        result.gaussian_over_tv.value,
        gaussian_lower.value,
        log4_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.gaussian_over_tv.value,
        result.gaussian_over_tv.value,
        error_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.gaussian_over_tau.value,
        gaussian_lower.value,
        log4_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.gaussian_over_tau.value,
        result.gaussian_over_tau.value,
        tau_upper.value,
        MPFR_RNDD
    );

    // The push lower bound dominates twice each of the two negative terms.
    mpfr_sub(
        result.push_over_twice_negative_tail.value,
        push_lower.value,
        log2_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.push_over_twice_negative_tail.value,
        result.push_over_twice_negative_tail.value,
        negative_term_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.push_over_twice_central_tail.value,
        push_lower.value,
        log2_upper.value,
        MPFR_RNDD
    );
    mpfr_sub(
        result.push_over_twice_central_tail.value,
        result.push_over_twice_central_tail.value,
        central_term_upper.value,
        MPFR_RNDD
    );
    return result;
}

bool positive(const Real& x) { return mpfr_sgn(x.value) > 0; }

}  // namespace

int main() {
    std::cout << "B-unitary-square post-m29 trace certificate\n";
    std::cout << "mpfr_version=" << mpfr_get_version()
              << " precision_bits=" << PRECISION << '\n';
    std::cout << "source_identity=Rains_2003_Theorem_1.5_equation_23_p2\n";
    std::cout << "source_concentration=CJL_2022_Lemma_2.2_m1\n";
    std::cout << "source_positive_tail=CJL_2022_Theorem_1.4\n";
    std::cout << "r=" << R_NUM << '/' << R_DEN
              << " eta=" << ETA_NUM << '/' << ETA_DEN
              << " ranks=B_" << B_MIN << "..B_" << B_MAX << '\n';

    Real log2_upper;
    log_ui(log2_upper, 2, MPFR_RNDU);
    Real log4_upper;
    mpfr_mul_ui(log4_upper.value, log2_upper.value, 2, MPFR_RNDU);
    const auto factorial_logs = factorial_log_lowers(2 * B_MAX);

    int failures = 0;
    std::vector<RowMargins> rows;
    rows.reserve(static_cast<std::size_t>(B_MAX - B_MIN + 1));
    for (int rank = B_MIN; rank <= B_MAX; ++rank) {
        RowMargins row = verify_row(
            rank,
            factorial_logs,
            log2_upper,
            log4_upper
        );
        const bool ok =
            positive(row.gaussian_over_tv)
            && positive(row.gaussian_over_tau)
            && positive(row.push_over_twice_negative_tail)
            && positive(row.push_over_twice_central_tail);
        if (!ok) ++failures;
        std::cout << "ROW B_" << rank
                  << " gaussian_over_tv=" << format_real(row.gaussian_over_tv)
                  << " gaussian_over_tau=" << format_real(row.gaussian_over_tau)
                  << " push_over_2negative="
                  << format_real(row.push_over_twice_negative_tail)
                  << " push_over_2central="
                  << format_real(row.push_over_twice_central_tail)
                  << " ok=" << (ok ? 1 : 0) << '\n';
        rows.push_back(std::move(row));
    }

    auto print_minimum = [&](const char* name, auto member) {
        const RowMargins* minimum = &rows.front();
        for (const RowMargins& row : rows) {
            if (mpfr_cmp((row.*member).value, (minimum->*member).value) < 0) {
                minimum = &row;
            }
        }
        std::cout << "MINIMUM " << name << " B_" << minimum->rank
                  << " value=" << format_real(minimum->*member) << '\n';
    };

    print_minimum("gaussian_over_tv", &RowMargins::gaussian_over_tv);
    print_minimum("gaussian_over_tau", &RowMargins::gaussian_over_tau);
    print_minimum(
        "push_over_2negative",
        &RowMargins::push_over_twice_negative_tail
    );
    print_minimum(
        "push_over_2central",
        &RowMargins::push_over_twice_central_tail
    );
    std::cout << "SUMMARY rows_checked=" << rows.size()
              << " failures=" << failures << '\n';
    std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << '\n';
    return failures == 0 ? 0 : 1;
}
