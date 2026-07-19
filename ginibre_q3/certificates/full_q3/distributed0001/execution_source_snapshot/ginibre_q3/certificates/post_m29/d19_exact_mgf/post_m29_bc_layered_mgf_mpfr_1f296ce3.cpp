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
int mpfr_set_z(mpfr_ptr, mpz_srcptr, mpfr_rnd_t);
int mpfr_add(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_sub(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_mul(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_mul_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_div(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_div_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_ui_div(mpfr_ptr, unsigned long, mpfr_srcptr, mpfr_rnd_t);
int mpfr_sqr(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_sqrt(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_pow_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_exp(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_log(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_neg(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_cmp(mpfr_srcptr, mpfr_srcptr);
int mpfr_cmp_ui(mpfr_srcptr, unsigned long);
int mpfr_sgn(mpfr_srcptr);
int mpfr_asprintf(char**, const char*, ...);
void mpfr_free_str(char*);
const char* mpfr_get_version(void);
}
#endif

#include <gmpxx.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr mpfr_prec_t PRECISION = 384;
constexpr int LAYERS = 8;

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

struct Interval {
    Real lo;
    Real hi;
};

[[noreturn]] void die(const std::string& message) {
    std::cerr << "ERROR " << message << '\n';
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

std::string format_optional_real(
    const Real& x,
    bool initialized,
    int digits = 28
) {
    return initialized ? format_real(x, digits) : "not_run";
}

Interval exact_ui(unsigned long value) {
    Interval out;
    mpfr_set_ui(out.lo.value, value, MPFR_RNDD);
    mpfr_set_ui(out.hi.value, value, MPFR_RNDU);
    return out;
}

Interval exact_z(const mpz_class& value) {
    Interval out;
    mpfr_set_z(out.lo.value, value.get_mpz_t(), MPFR_RNDD);
    mpfr_set_z(out.hi.value, value.get_mpz_t(), MPFR_RNDU);
    return out;
}

Interval rational(unsigned long numerator, unsigned long denominator) {
    Interval out = exact_ui(numerator);
    mpfr_div_ui(out.lo.value, out.lo.value, denominator, MPFR_RNDD);
    mpfr_div_ui(out.hi.value, out.hi.value, denominator, MPFR_RNDU);
    return out;
}

Interval add(const Interval& a, const Interval& b) {
    Interval out;
    mpfr_add(out.lo.value, a.lo.value, b.lo.value, MPFR_RNDD);
    mpfr_add(out.hi.value, a.hi.value, b.hi.value, MPFR_RNDU);
    return out;
}

Interval sub(const Interval& a, const Interval& b) {
    Interval out;
    mpfr_sub(out.lo.value, a.lo.value, b.hi.value, MPFR_RNDD);
    mpfr_sub(out.hi.value, a.hi.value, b.lo.value, MPFR_RNDU);
    return out;
}

void assign_min(Real& out, const Real& candidate) {
    if (mpfr_cmp(candidate.value, out.value) < 0) {
        mpfr_set(out.value, candidate.value, MPFR_RNDD);
    }
}

void assign_max(Real& out, const Real& candidate) {
    if (mpfr_cmp(candidate.value, out.value) > 0) {
        mpfr_set(out.value, candidate.value, MPFR_RNDU);
    }
}

Interval mul(const Interval& a, const Interval& b) {
    std::array<Real, 4> lower;
    std::array<Real, 4> upper;
    mpfr_mul(lower[0].value, a.lo.value, b.lo.value, MPFR_RNDD);
    mpfr_mul(lower[1].value, a.lo.value, b.hi.value, MPFR_RNDD);
    mpfr_mul(lower[2].value, a.hi.value, b.lo.value, MPFR_RNDD);
    mpfr_mul(lower[3].value, a.hi.value, b.hi.value, MPFR_RNDD);
    mpfr_mul(upper[0].value, a.lo.value, b.lo.value, MPFR_RNDU);
    mpfr_mul(upper[1].value, a.lo.value, b.hi.value, MPFR_RNDU);
    mpfr_mul(upper[2].value, a.hi.value, b.lo.value, MPFR_RNDU);
    mpfr_mul(upper[3].value, a.hi.value, b.hi.value, MPFR_RNDU);

    Interval out;
    mpfr_set(out.lo.value, lower[0].value, MPFR_RNDD);
    mpfr_set(out.hi.value, upper[0].value, MPFR_RNDU);
    for (int i = 1; i < 4; ++i) {
        assign_min(out.lo, lower[static_cast<std::size_t>(i)]);
        assign_max(out.hi, upper[static_cast<std::size_t>(i)]);
    }
    return out;
}

Interval neg(const Interval& a) {
    Interval out;
    mpfr_neg(out.lo.value, a.hi.value, MPFR_RNDD);
    mpfr_neg(out.hi.value, a.lo.value, MPFR_RNDU);
    return out;
}

Interval reciprocal_positive(const Interval& a) {
    if (mpfr_sgn(a.lo.value) <= 0) die("reciprocal denominator is not positive");
    Interval out;
    mpfr_ui_div(out.lo.value, 1, a.hi.value, MPFR_RNDD);
    mpfr_ui_div(out.hi.value, 1, a.lo.value, MPFR_RNDU);
    return out;
}

Interval div_positive(const Interval& a, const Interval& b) {
    return mul(a, reciprocal_positive(b));
}

Interval square(const Interval& a) {
    if (mpfr_sgn(a.lo.value) >= 0) {
        Interval out;
        mpfr_sqr(out.lo.value, a.lo.value, MPFR_RNDD);
        mpfr_sqr(out.hi.value, a.hi.value, MPFR_RNDU);
        return out;
    }
    if (mpfr_sgn(a.hi.value) <= 0) return square(neg(a));
    Interval out;
    mpfr_set_zero(out.lo.value, 0);
    Real left;
    Real right;
    mpfr_sqr(left.value, a.lo.value, MPFR_RNDU);
    mpfr_sqr(right.value, a.hi.value, MPFR_RNDU);
    if (mpfr_cmp(left.value, right.value) >= 0) {
        mpfr_set(out.hi.value, left.value, MPFR_RNDU);
    } else {
        mpfr_set(out.hi.value, right.value, MPFR_RNDU);
    }
    return out;
}

Interval pow_positive(const Interval& a, unsigned long exponent) {
    if (mpfr_sgn(a.lo.value) < 0) die("positive integer power received negative interval");
    Interval out;
    mpfr_pow_ui(out.lo.value, a.lo.value, exponent, MPFR_RNDD);
    mpfr_pow_ui(out.hi.value, a.hi.value, exponent, MPFR_RNDU);
    return out;
}

Interval sqrt_positive(const Interval& a) {
    if (mpfr_sgn(a.lo.value) < 0) die("square root received negative interval");
    Interval out;
    mpfr_sqrt(out.lo.value, a.lo.value, MPFR_RNDD);
    mpfr_sqrt(out.hi.value, a.hi.value, MPFR_RNDU);
    return out;
}

Interval exp_interval(const Interval& a) {
    Interval out;
    mpfr_exp(out.lo.value, a.lo.value, MPFR_RNDD);
    mpfr_exp(out.hi.value, a.hi.value, MPFR_RNDU);
    return out;
}

Interval log_positive(const Interval& a) {
    if (mpfr_sgn(a.lo.value) <= 0) die("logarithm received nonpositive interval");
    Interval out;
    mpfr_log(out.lo.value, a.lo.value, MPFR_RNDD);
    mpfr_log(out.hi.value, a.hi.value, MPFR_RNDU);
    return out;
}

Interval scale_ui(const Interval& a, unsigned long value) {
    return mul(a, exact_ui(value));
}

bool certainly_positive(const Interval& a) {
    return mpfr_sgn(a.lo.value) > 0;
}

bool certainly_less_than_one(const Interval& a) {
    return mpfr_cmp_ui(a.hi.value, 1) < 0;
}

Interval one_minus(const Interval& a) {
    return sub(exact_ui(1), a);
}

unsigned long binomial_ui(unsigned n, unsigned k) {
    if (k > n) return 0;
    if (k > n - k) k = n - k;
    unsigned long answer = 1;
    for (unsigned i = 1; i <= k; ++i) {
        answer = answer * (n - k + i) / i;
    }
    return answer;
}

std::vector<mpz_class> stable_push_moments(int maximum) {
    std::vector<mpq_class> gaussian(static_cast<std::size_t>(maximum + 1));
    gaussian[0] = 1;
    for (int k = 0; k < maximum; ++k) {
        mpq_class previous = (k == 0) ? mpq_class(0) : gaussian[static_cast<std::size_t>(k - 1)];
        gaussian[static_cast<std::size_t>(k + 1)] =
            (-gaussian[static_cast<std::size_t>(k)] + previous) / (k + 1);
    }

    std::vector<mpz_class> factorial(static_cast<std::size_t>(maximum + 1));
    factorial[0] = 1;
    for (int k = 1; k <= maximum; ++k) {
        factorial[static_cast<std::size_t>(k)] =
            factorial[static_cast<std::size_t>(k - 1)] * k;
    }

    std::vector<mpz_class> moments(static_cast<std::size_t>(maximum + 1));
    for (int q = 0; q <= maximum; ++q) {
        mpq_class coefficient = 0;
        for (int j = 0; j <= q; ++j) {
            const mpq_class pole = mpq_class((j + 2) * (j + 1), 2) + 1;
            coefficient += pole * gaussian[static_cast<std::size_t>(q - j)];
        }
        coefficient *= factorial[static_cast<std::size_t>(q)];
        coefficient.canonicalize();
        if (coefficient.get_den() != 1) die("stable push moment is not integral");
        moments[static_cast<std::size_t>(q)] = coefficient.get_num();
        if (moments[static_cast<std::size_t>(q)] <= 0) {
            die("stable push moment is not positive");
        }
    }
    return moments;
}

std::vector<mpz_class> stable_adjoint_moments(int maximum) {
    std::vector<mpz_class> factorials(static_cast<std::size_t>(maximum + 1));
    factorials[0] = 1;
    for (int index = 1; index <= maximum; ++index) {
        factorials[static_cast<std::size_t>(index)] =
            factorials[static_cast<std::size_t>(index - 1)] * index;
    }

    std::vector<mpz_class> moments(static_cast<std::size_t>(maximum + 1));
    moments[0] = 1;
    for (int n = 0; n < maximum; ++n) {
        mpz_class next = 0;
        for (int index = 1; index <= n + 1; ++index) {
            mpz_class cumulant;
            if (index == 1) {
                cumulant = 0;
            } else if (index == 2) {
                cumulant = 1;
            } else {
                cumulant = factorials[static_cast<std::size_t>(index - 1)] / 2;
            }
            next += mpz_class(binomial_ui(
                static_cast<unsigned>(n),
                static_cast<unsigned>(index - 1)
            )) * cumulant * moments[static_cast<std::size_t>(n + 1 - index)];
        }
        moments[static_cast<std::size_t>(n + 1)] = std::move(next);
    }
    return moments;
}

std::vector<mpz_class> push_moments_from_adjoint(
    const std::vector<mpz_class>& adjoint_moments,
    int maximum
) {
    if (static_cast<int>(adjoint_moments.size()) <= maximum + 2) {
        die("insufficient adjoint moments for push-moment reconstruction");
    }
    std::vector<mpz_class> push_moments(static_cast<std::size_t>(maximum + 1));
    for (int exponent = 0; exponent <= maximum; ++exponent) {
        mpz_class sum = 0;
        for (int index = 0; index <= exponent; ++index) {
            const mpz_class coefficient = binomial_ui(
                static_cast<unsigned>(exponent),
                static_cast<unsigned>(index)
            );
            sum += coefficient * (
                adjoint_moments[static_cast<std::size_t>(index + 2)]
                    * adjoint_moments[static_cast<std::size_t>(exponent - index)]
                - adjoint_moments[static_cast<std::size_t>(index + 1)]
                    * adjoint_moments[static_cast<std::size_t>(exponent - index + 1)]
            );
        }
        push_moments[static_cast<std::size_t>(exponent)] = 2 * sum;
        if (push_moments[static_cast<std::size_t>(exponent)] <= 0) {
            die("reconstructed push moment is not positive");
        }
    }
    return push_moments;
}

// CJL Lemma 3.1 reduces the SO(2b+1) defining-trace MGF to a finite
// Toeplitz--Hankel determinant.  The entries use I_k(2v), whose positive
// series is summed here with an explicit geometric remainder.  All interval
// endpoints are rounded outwards.
Interval bessel_i_twice_argument(
    int order,
    const Interval& argument,
    const std::vector<mpz_class>& factorials
) {
    if (order < 0) die("negative Bessel order");
    if (!certainly_positive(argument)) die("nonpositive Bessel argument");
    if (static_cast<std::size_t>(order) >= factorials.size()) {
        die("Bessel factorial table is too short");
    }

    Interval term = div_positive(
        pow_positive(argument, static_cast<unsigned long>(order)),
        exact_z(factorials[static_cast<std::size_t>(order)])
    );
    Interval sum = term;
    const Interval argument_sq = square(argument);
    constexpr int SERIES_TERMS = 180;
    for (int m = 0; m < SERIES_TERMS; ++m) {
        const unsigned long denominator =
            static_cast<unsigned long>(m + 1)
            * static_cast<unsigned long>(m + order + 1);
        term = div_positive(mul(term, argument_sq), exact_ui(denominator));
        sum = add(sum, term);
    }

    const unsigned long next_denominator =
        static_cast<unsigned long>(SERIES_TERMS + 1)
        * static_cast<unsigned long>(SERIES_TERMS + order + 1);
    const Interval next_ratio = div_positive(argument_sq, exact_ui(next_denominator));
    if (!certainly_less_than_one(next_ratio)) die("Bessel tail ratio reached one");
    const Interval first_omitted = mul(term, next_ratio);
    const Interval tail = div_positive(first_omitted, one_minus(next_ratio));
    mpfr_add(sum.hi.value, sum.hi.value, tail.hi.value, MPFR_RNDU);
    return sum;
}

Interval log_so_odd_defining_mgf(
    int rank,
    const Interval& argument,
    int sign,
    const std::vector<mpz_class>& factorials,
    Real& minimum_pivot
) {
    if (rank < 1) die("SO odd rank is not positive");
    if (sign != 1 && sign != -1) die("invalid MGF sign");
    std::vector<Interval> bessel(static_cast<std::size_t>(2 * rank));
    for (int order = 0; order < 2 * rank; ++order) {
        bessel[static_cast<std::size_t>(order)] =
            bessel_i_twice_argument(order, argument, factorials);
    }

    std::vector<std::vector<Interval>> matrix(
        static_cast<std::size_t>(rank),
        std::vector<Interval>(static_cast<std::size_t>(rank))
    );
    for (int row = 0; row < rank; ++row) {
        for (int column = 0; column < rank; ++column) {
            const Interval& toeplitz = bessel[static_cast<std::size_t>(std::abs(row - column))];
            const Interval& hankel = bessel[static_cast<std::size_t>(row + column + 1)];
            matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] =
                (sign == 1) ? sub(toeplitz, hankel) : add(toeplitz, hankel);
        }
    }

    // Both matrices are Gram matrices for a strictly positive Jacobi weight.
    // Interval Cholesky therefore also certifies positivity of every pivot.
    std::vector<std::vector<Interval>> cholesky(
        static_cast<std::size_t>(rank),
        std::vector<Interval>(static_cast<std::size_t>(rank))
    );
    Interval log_determinant = exact_ui(0);
    for (int row = 0; row < rank; ++row) {
        for (int column = 0; column < row; ++column) {
            Interval numerator = matrix[static_cast<std::size_t>(row)]
                                          [static_cast<std::size_t>(column)];
            for (int k = 0; k < column; ++k) {
                numerator = sub(
                    numerator,
                    mul(
                        cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(k)],
                        cholesky[static_cast<std::size_t>(column)][static_cast<std::size_t>(k)]
                    )
                );
            }
            cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] =
                div_positive(
                    numerator,
                    cholesky[static_cast<std::size_t>(column)][static_cast<std::size_t>(column)]
                );
        }
        Interval pivot = matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
        for (int k = 0; k < row; ++k) {
            pivot = sub(
                pivot,
                square(cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(k)])
            );
        }
        if (!certainly_positive(pivot)) die("SO odd interval Cholesky pivot failed");
        if (mpfr_cmp(pivot.lo.value, minimum_pivot.value) < 0) {
            mpfr_set(minimum_pivot.value, pivot.lo.value, MPFR_RNDD);
        }
        cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)] =
            sqrt_positive(pivot);
        log_determinant = add(log_determinant, log_positive(pivot));
    }

    return (sign == 1)
        ? add(argument, log_determinant)
        : sub(log_determinant, argument);
}

// CJL Lemma 3.1, the E_n^{++} case.  This is the defining-trace MGF of
// O^-(2n+2): the fixed +1 and -1 eigenvalues cancel, and the random part has
// n conjugate pairs.  For psi(x)=exp(2vx), the matrix is
// I_{|j-k|}(2v)-I_{j+k+2}(2v), 0<=j,k<n.
[[maybe_unused]] Interval log_o_even_minus_defining_mgf(
    int total_dimension,
    const Interval& argument,
    const std::vector<mpz_class>& factorials,
    Real& minimum_pivot
) {
    if (total_dimension < 4 || total_dimension % 2 != 0) {
        die("O-even-minus dimension must be even and at least four");
    }
    const int pairs = total_dimension / 2 - 1;
    std::vector<Interval> bessel(static_cast<std::size_t>(2 * pairs + 1));
    for (int order = 0; order <= 2 * pairs; ++order) {
        bessel[static_cast<std::size_t>(order)] =
            bessel_i_twice_argument(order, argument, factorials);
    }

    std::vector<std::vector<Interval>> matrix(
        static_cast<std::size_t>(pairs),
        std::vector<Interval>(static_cast<std::size_t>(pairs))
    );
    for (int row = 0; row < pairs; ++row) {
        for (int column = 0; column < pairs; ++column) {
            matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] = sub(
                bessel[static_cast<std::size_t>(std::abs(row - column))],
                bessel[static_cast<std::size_t>(row + column + 2)]
            );
        }
    }

    std::vector<std::vector<Interval>> cholesky(
        static_cast<std::size_t>(pairs),
        std::vector<Interval>(static_cast<std::size_t>(pairs))
    );
    Interval log_determinant = exact_ui(0);
    for (int row = 0; row < pairs; ++row) {
        for (int column = 0; column < row; ++column) {
            Interval numerator = matrix[static_cast<std::size_t>(row)]
                                       [static_cast<std::size_t>(column)];
            for (int k = 0; k < column; ++k) {
                numerator = sub(
                    numerator,
                    mul(
                        cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(k)],
                        cholesky[static_cast<std::size_t>(column)][static_cast<std::size_t>(k)]
                    )
                );
            }
            cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] =
                div_positive(
                    numerator,
                    cholesky[static_cast<std::size_t>(column)][static_cast<std::size_t>(column)]
                );
        }
        Interval pivot = matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
        for (int k = 0; k < row; ++k) {
            pivot = sub(
                pivot,
                square(cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(k)])
            );
        }
        if (!certainly_positive(pivot)) die("O-even-minus interval Cholesky pivot failed");
        if (mpfr_cmp(pivot.lo.value, minimum_pivot.value) < 0) {
            mpfr_set(minimum_pivot.value, pivot.lo.value, MPFR_RNDD);
        }
        cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)] =
            sqrt_positive(pivot);
        log_determinant = add(log_determinant, log_positive(pivot));
    }
    return log_determinant;
}

// CJL Lemma 3.1, the E_n^{--} case, normalized at v=0.  This is the
// defining-trace MGF of SO(2n)=O^+(2n).  The factor 1/2 compensates for the
// doubled (0,0) entry in the Toeplitz--Hankel determinant.
Interval log_so_even_defining_mgf(
    int rank,
    const Interval& argument,
    const std::vector<mpz_class>& factorials,
    Real& minimum_pivot
) {
    if (rank < 2) die("SO even rank must be at least two");
    std::vector<Interval> bessel(static_cast<std::size_t>(2 * rank - 1));
    for (int order = 0; order <= 2 * rank - 2; ++order) {
        bessel[static_cast<std::size_t>(order)] =
            bessel_i_twice_argument(order, argument, factorials);
    }

    std::vector<std::vector<Interval>> matrix(
        static_cast<std::size_t>(rank),
        std::vector<Interval>(static_cast<std::size_t>(rank))
    );
    for (int row = 0; row < rank; ++row) {
        for (int column = 0; column < rank; ++column) {
            matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] = add(
                bessel[static_cast<std::size_t>(std::abs(row - column))],
                bessel[static_cast<std::size_t>(row + column)]
            );
        }
    }

    std::vector<std::vector<Interval>> cholesky(
        static_cast<std::size_t>(rank),
        std::vector<Interval>(static_cast<std::size_t>(rank))
    );
    Interval log_determinant = neg(log_positive(exact_ui(2)));
    for (int row = 0; row < rank; ++row) {
        for (int column = 0; column < row; ++column) {
            Interval numerator = matrix[static_cast<std::size_t>(row)]
                                       [static_cast<std::size_t>(column)];
            for (int k = 0; k < column; ++k) {
                numerator = sub(
                    numerator,
                    mul(
                        cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(k)],
                        cholesky[static_cast<std::size_t>(column)][static_cast<std::size_t>(k)]
                    )
                );
            }
            cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] =
                div_positive(
                    numerator,
                    cholesky[static_cast<std::size_t>(column)][static_cast<std::size_t>(column)]
                );
        }
        Interval pivot = matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
        for (int k = 0; k < row; ++k) {
            pivot = sub(
                pivot,
                square(cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(k)])
            );
        }
        if (!certainly_positive(pivot)) die("SO-even interval Cholesky pivot failed");
        if (mpfr_cmp(pivot.lo.value, minimum_pivot.value) < 0) {
            mpfr_set(minimum_pivot.value, pivot.lo.value, MPFR_RNDD);
        }
        cholesky[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)] =
            sqrt_positive(pivot);
        log_determinant = add(log_determinant, log_positive(pivot));
    }
    return log_determinant;
}

struct FredholmBounds {
    Interval beta;
    Interval trace_norm;
    Interval determinant_lower;
    Interval determinant_upper;
};

FredholmBounds fredholm_bounds(
    int dimension,
    const Interval& argument,
    const mpz_class& factorial
) {
    if (!certainly_positive(argument)) die("Fredholm argument is not positive");
    const Interval d_plus_one = exact_ui(static_cast<unsigned long>(dimension + 1));
    if (mpfr_cmp(argument.hi.value, d_plus_one.lo.value) >= 0) {
        die("Fredholm argument reached d+1");
    }

    const Interval argument_sq = square(argument);
    const Interval exp_factor = exp_interval(div_positive(argument_sq, d_plus_one));
    const Interval power = pow_positive(argument, static_cast<unsigned long>(dimension));
    const Interval factorial_interval = exact_z(factorial);
    const Interval ratio = div_positive(argument, d_plus_one);
    const Interval ratio_sq = square(ratio);
    const Interval one_minus_ratio_sq = one_minus(ratio_sq);
    if (!certainly_positive(one_minus_ratio_sq)) die("Fredholm quadratic denominator failed");

    FredholmBounds out;
    out.beta = div_positive(
        mul(exp_factor, power),
        mul(factorial_interval, one_minus_ratio_sq)
    );

    const Interval one_minus_ratio = one_minus(ratio);
    if (!certainly_positive(one_minus_ratio)) die("Fredholm linear denominator failed");
    out.trace_norm = div_positive(
        mul(exp_factor, power),
        mul(
            factorial_interval,
            mul(one_minus_ratio, sqrt_positive(one_minus_ratio_sq))
        )
    );

    if (!certainly_less_than_one(out.beta)) {
        die(
            "Fredholm beta reached one: dimension=" + std::to_string(dimension)
            + " argument_hi=" + format_real(argument.hi, 20)
            + " beta_hi=" + format_real(out.beta.hi, 20)
        );
    }
    if (!certainly_less_than_one(out.trace_norm)) {
        die(
            "Fredholm trace norm reached one: dimension=" + std::to_string(dimension)
            + " argument_hi=" + format_real(argument.hi, 20)
            + " trace_norm_hi=" + format_real(out.trace_norm.hi, 20)
        );
    }
    out.determinant_lower = one_minus(out.beta);
    out.determinant_upper = reciprocal_positive(out.determinant_lower);
    return out;
}

struct RationalParameter {
    unsigned long numerator;
    unsigned long denominator;
};

struct ParameterSet {
    const char* label;
    char family;
    int rank_low;
    int rank_high;
    RationalParameter r;
    RationalParameter square_alpha;
    RationalParameter negative_alpha;
    RationalParameter character_cutoff;
    RationalParameter lower_gap;
    RationalParameter upper_start;
    RationalParameter upper_step;
    RationalParameter lower_lambda_fraction;
    RationalParameter upper_lambda_fraction;
    RationalParameter square_lambda_fraction;
    RationalParameter negative_lambda_fraction;
    int polynomial_k;
    int polynomial_m;
    bool search_from_d_first_hit;
};

using DAdjointMomentMap = std::map<std::pair<int, int>, mpz_class>;

void load_d_adjoint_moment_source_log(
    const std::string& path,
    DAdjointMomentMap& moments
) {
    std::ifstream input(path);
    if (!input) die("cannot open D adjoint-moment source log: " + path);

    std::string line;
    int loaded = 0;
    bool success = false;
    bool have_host = false;
    bool have_source_hash = false;
    while (std::getline(input, line)) {
        if (line == "__EXIT_STATUS=0"
            || line == "__PART_EXIT_STATUS=0") {
            success = true;
        }
        if (line.rfind("__HOST=", 0) == 0 && line.size() > 7) {
            have_host = true;
        }
        if (line.rfind("__SOURCE_SHA256=", 0) == 0 && line.size() == 80) {
            have_source_hash = true;
        }
        const std::size_t rank_marker = line.find("D_");
        if (rank_marker == std::string::npos) continue;
        const std::size_t index_marker = line.find(" m_", rank_marker + 2);
        const std::size_t value_marker = line.find("; moment_", index_marker);
        if (index_marker == std::string::npos
            || value_marker == std::string::npos) {
            continue;
        }
        const std::size_t equals_marker = line.find(" = ", value_marker);
        if (equals_marker == std::string::npos) continue;

        const int rank = std::stoi(line.substr(
            rank_marker + 2,
            index_marker - (rank_marker + 2)
        ));
        const int index = std::stoi(line.substr(
            index_marker + 3,
            value_marker - (index_marker + 3)
        ));
        const std::size_t value_index_marker = value_marker + 9;
        const int printed_index = std::stoi(line.substr(
            value_index_marker,
            equals_marker - value_index_marker
        ));
        if (printed_index != index) {
            die("mismatched moment index in D adjoint-moment source log");
        }

        mpz_class value;
        const std::string digits = line.substr(equals_marker + 3);
        if (mpz_set_str(value.get_mpz_t(), digits.c_str(), 10) != 0) {
            die("invalid integer in D adjoint-moment source log");
        }
        const std::pair<int, int> key{rank, index};
        auto found = moments.find(key);
        if (found != moments.end() && found->second != value) {
            die("conflicting D adjoint moment in source logs");
        }
        if (found == moments.end()) {
            moments.emplace(key, std::move(value));
            ++loaded;
        }
    }
    if (!success) {
        die("D adjoint-moment source log lacks __EXIT_STATUS=0: " + path);
    }
    if (!have_host || !have_source_hash) {
        die("D adjoint-moment source log lacks machine/source provenance: " + path);
    }
    if (loaded == 0) die("D adjoint-moment source log supplied no moment rows");
}

const std::array<ParameterSet, 5> PARAMETER_SETS{{
    {
        "B16-special", 'B', 16, 16,
        {200001, 100000}, {45546, 10000}, {4284, 1000}, {1277, 1000},
        {856, 1000}, {1847, 1000}, {34, 1000}, {1, 1}, {567, 1000},
        {1, 1}, {101331, 100000}, 5, 10, false
    },
    {
        "B17-special", 'B', 17, 17,
        {20446, 10000}, {44739, 10000}, {42902, 10000}, {1341, 1000},
        {972, 1000}, {15366, 10000}, {332, 10000}, {1, 1}, {7291, 10000},
        {1, 1}, {101331, 100000}, 5, 10, false
    },
    {
        "B18-61-uniform", 'B', 18, 61,
        {20674, 10000}, {441721, 100000}, {429465, 100000}, {13501, 10000},
        {103311, 100000}, {139822, 100000}, {31205, 1000000}, {1, 1}, {1, 1},
        {1, 1}, {101331, 100000}, 5, 8, false
    },
    {
        "C20-61-uniform", 'C', 20, 61,
        {226172, 100000}, {460576, 100000}, {438971, 100000}, {138552, 100000},
        {105713, 100000}, {135104, 100000}, {30177, 1000000}, {1, 1}, {1, 1},
        {75841, 100000}, {79538, 100000}, 5, 8, false
    },
    {
        "D53-295-short-onset", 'D', 53, 295,
        {20674, 10000}, {441721, 100000}, {429465, 100000}, {13501, 10000},
        {103311, 100000}, {139822, 100000}, {31205, 1000000}, {1, 1}, {1, 1},
        {1, 1}, {101331, 100000}, 5, 8, true
    }
}};

Interval parameter(const RationalParameter& value) {
    return rational(value.numerator, value.denominator);
}

struct ExactLayeredTail {
    Interval probability;
    Interval first_mass;
    Interval monotonicity;
};

ExactLayeredTail exact_so_even_layered_tail(
    int rank,
    const Interval& saddle,
    const Interval& lower_gap,
    const Interval& upper_start,
    const Interval& upper_step,
    const std::vector<mpz_class>& factorials,
    Real& minimum_pivot
) {
    const Interval log_mgf_s = log_so_even_defining_mgf(
        rank, saddle, factorials, minimum_pivot
    );
    const Interval lower_argument = sub(saddle, lower_gap);
    const Interval log_mgf_lower = log_so_even_defining_mgf(
        rank, lower_argument, factorials, minimum_pivot
    );
    const Interval lower_exponent = add(
        sub(mul(lower_gap, saddle), square(lower_gap)),
        sub(log_mgf_lower, log_mgf_s)
    );
    const Interval lower_tail = exp_interval(lower_exponent);

    std::array<Interval, LAYERS> cumulative_mass;
    std::array<Interval, LAYERS> weights;
    for (int i = 0; i < LAYERS; ++i) {
        const Interval upper_gap = add(
            upper_start,
            scale_ui(upper_step, static_cast<unsigned long>(i))
        );
        const Interval upper_argument = add(saddle, upper_gap);
        const Interval log_mgf_upper = log_so_even_defining_mgf(
            rank, upper_argument, factorials, minimum_pivot
        );
        const Interval upper_exponent = add(
            neg(add(mul(upper_gap, saddle), square(upper_gap))),
            sub(log_mgf_upper, log_mgf_s)
        );
        const Interval upper_tail = exp_interval(upper_exponent);
        cumulative_mass[static_cast<std::size_t>(i)] = sub(
            sub(exact_ui(1), lower_tail),
            upper_tail
        );
        if (!certainly_positive(cumulative_mass[static_cast<std::size_t>(i)])) {
            die("exact SO-even layered tilted mass failed");
        }
        weights[static_cast<std::size_t>(i)] = exp_interval(
            neg(mul(saddle, upper_gap))
        );
    }

    Interval monotonicity = sub(cumulative_mass[1], cumulative_mass[0]);
    for (int i = 1; i + 1 < LAYERS; ++i) {
        const Interval difference = sub(
            cumulative_mass[static_cast<std::size_t>(i + 1)],
            cumulative_mass[static_cast<std::size_t>(i)]
        );
        if (mpfr_cmp(difference.lo.value, monotonicity.lo.value) < 0) {
            mpfr_set(monotonicity.lo.value, difference.lo.value, MPFR_RNDD);
        }
    }
    if (!certainly_positive(monotonicity)) {
        die("exact SO-even layer monotonicity failed");
    }

    Interval layer_sum = mul(weights[LAYERS - 1], cumulative_mass[LAYERS - 1]);
    for (int i = 0; i + 1 < LAYERS; ++i) {
        const Interval weight_difference = sub(
            weights[static_cast<std::size_t>(i)],
            weights[static_cast<std::size_t>(i + 1)]
        );
        if (!certainly_positive(weight_difference)) {
            die("exact SO-even layer weights are not decreasing");
        }
        layer_sum = add(
            layer_sum,
            mul(weight_difference, cumulative_mass[static_cast<std::size_t>(i)])
        );
    }
    if (!certainly_positive(layer_sum)) die("exact SO-even layer sum failed");

    ExactLayeredTail out;
    out.probability = mul(
        exp_interval(sub(log_mgf_s, square(saddle))),
        layer_sum
    );
    out.first_mass = cumulative_mass[0];
    out.monotonicity = monotonicity;
    return out;
}

struct RowMargins {
    std::string label;
    int rank = 0;
    int exponent = 0;
    Real geometry;
    Real fredholm;
    Real exact_defining_pivot;
    Real exact_square_pivot;
    Real first_layer_mass;
    Real layer_monotonicity;
    Real positive_difference;
    Real polynomial;
    Real final_margin;
};

void update_min(Real& current, const Real& candidate, bool& initialized) {
    if (!initialized || mpfr_cmp(candidate.value, current.value) < 0) {
        mpfr_set(current.value, candidate.value, MPFR_RNDD);
        initialized = true;
    }
}

Interval square_tail_upper(
    char family,
    int rank,
    const Interval& cutoff,
    const Interval& lambda,
    const std::vector<mpz_class>& factorials,
    Real& fredholm_margin,
    Real& exact_square_pivot,
    bool exact_d_square
) {
    Interval exponent = add(
        neg(mul(lambda, sub(cutoff, exact_ui(1)))),
        square(lambda)
    );
    Interval answer = exp_interval(exponent);
    if (family == 'B') return answer;
    if (family != 'C' && family != 'D') die("unknown square-trace family");

    const int m0 = family == 'D' ? (rank + 1) / 2 : (rank + 2) / 2;
    const int m1 = family == 'D' ? rank / 2 : (rank + 1) / 2;
    const int d0 = 2 * m0;
    const int d1 = 2 * m1 + 1;
    if (family == 'D' && exact_d_square) {
        // Rains (22) gives O^+(2m0) and O^-(2m1+1), with forced
        // eigenvalues omitted.  Restoring the odd component's fixed -1 gives
        // tau_D - 1 = Y_even,+ + Y_odd,-.
        const Interval log_first = log_so_even_defining_mgf(
            m0,
            lambda,
            factorials,
            exact_square_pivot
        );
        const Interval log_second = log_so_odd_defining_mgf(
            m1,
            lambda,
            -1,
            factorials,
            exact_square_pivot
        );
        return exp_interval(add(
            neg(mul(lambda, sub(cutoff, exact_ui(1)))),
            add(log_first, log_second)
        ));
    }
    for (int dimension : {d0, d1}) {
        FredholmBounds bounds = fredholm_bounds(
            dimension,
            lambda,
            factorials[static_cast<std::size_t>(dimension)]
        );
        Interval trace_margin = one_minus(bounds.trace_norm);
        if (mpfr_cmp(trace_margin.lo.value, fredholm_margin.value) < 0) {
            mpfr_set(fredholm_margin.value, trace_margin.lo.value, MPFR_RNDD);
        }
        answer = mul(answer, bounds.determinant_upper);
    }
    return answer;
}

RowMargins verify_row(
    const ParameterSet& parameters,
    int rank,
    const std::vector<mpz_class>& factorials,
    const std::vector<mpz_class>& push_moments,
    bool enforce_d_a_free = true,
    bool exact_d_square = false,
    bool exact_d_defining = false,
    bool exact_d_push_moments = false
) {
    RowMargins margins;
    margins.label = parameters.label;
    margins.rank = rank;
    mpfr_set_ui(margins.fredholm.value, 1, MPFR_RNDD);
    mpfr_set_ui(margins.exact_defining_pivot.value, 1, MPFR_RNDD);
    mpfr_set_ui(margins.exact_square_pivot.value, 1, MPFR_RNDD);

    const Interval b = exact_ui(static_cast<unsigned long>(rank));
    const Interval sqrt_b = sqrt_positive(b);
    const Interval r = parameter(parameters.r);
    const Interval square_alpha = parameter(parameters.square_alpha);
    const Interval negative_alpha = parameter(parameters.negative_alpha);
    const Interval character_cutoff = parameter(parameters.character_cutoff);
    const Interval lower_gap = parameter(parameters.lower_gap);
    const Interval upper_start = parameter(parameters.upper_start);
    const Interval upper_step = parameter(parameters.upper_step);
    const Interval lower_fraction = parameter(parameters.lower_lambda_fraction);
    const Interval upper_fraction = parameter(parameters.upper_lambda_fraction);
    const Interval square_fraction = parameter(parameters.square_lambda_fraction);
    const Interval negative_fraction = parameter(parameters.negative_lambda_fraction);

    const Interval c = mul(r, b);
    const Interval q = mul(square_alpha, sqrt_b);
    const Interval a = mul(negative_alpha, sqrt_b);
    const Interval two_b = scale_ui(b, 2);
    const int endpoint = parameters.family == 'D' && rank % 2 != 0
        ? rank - 2
        : rank;
    const Interval two_endpoint = scale_ui(
        exact_ui(static_cast<unsigned long>(endpoint)),
        2
    );

    Interval geometry_margin = sub(c, two_b);
    Interval a_above_one = sub(a, exact_ui(1));
    Interval a_below_support = sub(two_endpoint, a);
    Interval d_above_one = sub(character_cutoff, exact_ui(1));
    for (const Interval* candidate : {&a_above_one, &a_below_support, &d_above_one}) {
        if (mpfr_cmp(candidate->lo.value, geometry_margin.lo.value) < 0) {
            mpfr_set(geometry_margin.lo.value, candidate->lo.value, MPFR_RNDD);
        }
    }
    if (!certainly_positive(geometry_margin)) die(margins.label + " geometry failed");
    mpfr_set(margins.geometry.value, geometry_margin.lo.value, MPFR_RNDD);

    const Interval threshold = sqrt_positive(add(
        add(scale_ui(c, 2), scale_ui(character_cutoff, 2)),
        q
    ));
    const Interval saddle = add(threshold, lower_gap);
    const int defining_dimension = (parameters.family == 'B') ? 2 * rank + 1 : 2 * rank;

    Interval defining_tail;
    if (parameters.family == 'D' && exact_d_defining) {
        ExactLayeredTail exact_tail = exact_so_even_layered_tail(
            rank,
            saddle,
            lower_gap,
            upper_start,
            upper_step,
            factorials,
            margins.exact_defining_pivot
        );
        defining_tail = exact_tail.probability;
        mpfr_set(
            margins.first_layer_mass.value,
            exact_tail.first_mass.lo.value,
            MPFR_RNDD
        );
        mpfr_set(
            margins.layer_monotonicity.value,
            exact_tail.monotonicity.lo.value,
            MPFR_RNDD
        );
    } else {
    FredholmBounds saddle_bounds = fredholm_bounds(
        defining_dimension,
        saddle,
        factorials[static_cast<std::size_t>(defining_dimension)]
    );
    Interval saddle_trace_margin = one_minus(saddle_bounds.trace_norm);
    mpfr_set(margins.fredholm.value, saddle_trace_margin.lo.value, MPFR_RNDD);
    const Interval saddle_det_lower = saddle_bounds.determinant_lower;

    const Interval lower_lambda = mul(lower_fraction, lower_gap);
    const Interval lower_argument = sub(saddle, lower_lambda);
    FredholmBounds lower_bounds = fredholm_bounds(
        defining_dimension,
        lower_argument,
        factorials[static_cast<std::size_t>(defining_dimension)]
    );
    Interval lower_trace_margin = one_minus(lower_bounds.trace_norm);
    if (mpfr_cmp(lower_trace_margin.lo.value, margins.fredholm.value) < 0) {
        mpfr_set(margins.fredholm.value, lower_trace_margin.lo.value, MPFR_RNDD);
    }

    Interval lower_exponent = sub(
        div_positive(square(lower_lambda), exact_ui(2)),
        mul(lower_lambda, lower_gap)
    );
    Interval lower_tail = div_positive(
        exp_interval(lower_exponent),
        mul(saddle_det_lower, lower_bounds.determinant_lower)
    );

    std::array<Interval, LAYERS> upper_gaps;
    std::array<Interval, LAYERS> cumulative_mass;
    std::array<Interval, LAYERS> weights;
    for (int i = 0; i < LAYERS; ++i) {
        upper_gaps[static_cast<std::size_t>(i)] = add(
            upper_start,
            scale_ui(upper_step, static_cast<unsigned long>(i))
        );
        const Interval lambda = mul(
            upper_fraction,
            upper_gaps[static_cast<std::size_t>(i)]
        );
        const Interval upper_argument = add(saddle, lambda);
        FredholmBounds upper_bounds = fredholm_bounds(
            defining_dimension,
            upper_argument,
            factorials[static_cast<std::size_t>(defining_dimension)]
        );
        Interval upper_trace_margin = one_minus(upper_bounds.trace_norm);
        if (mpfr_cmp(upper_trace_margin.lo.value, margins.fredholm.value) < 0) {
            mpfr_set(margins.fredholm.value, upper_trace_margin.lo.value, MPFR_RNDD);
        }

        const Interval upper_exponent = sub(
            div_positive(square(lambda), exact_ui(2)),
            mul(lambda, upper_gaps[static_cast<std::size_t>(i)])
        );
        const Interval upper_tail = div_positive(
            exp_interval(upper_exponent),
            mul(saddle_det_lower, upper_bounds.determinant_lower)
        );
        cumulative_mass[static_cast<std::size_t>(i)] = sub(
            sub(exact_ui(1), lower_tail),
            upper_tail
        );
        if (!certainly_positive(cumulative_mass[static_cast<std::size_t>(i)])) {
            die(margins.label + " layered tilted mass failed");
        }
        weights[static_cast<std::size_t>(i)] = exp_interval(
            neg(mul(saddle, upper_gaps[static_cast<std::size_t>(i)]))
        );
    }
    mpfr_set(
        margins.first_layer_mass.value,
        cumulative_mass[0].lo.value,
        MPFR_RNDD
    );

    Interval layer_monotonicity = sub(cumulative_mass[1], cumulative_mass[0]);
    for (int i = 1; i + 1 < LAYERS; ++i) {
        Interval difference = sub(
            cumulative_mass[static_cast<std::size_t>(i + 1)],
            cumulative_mass[static_cast<std::size_t>(i)]
        );
        if (mpfr_cmp(difference.lo.value, layer_monotonicity.lo.value) < 0) {
            mpfr_set(layer_monotonicity.lo.value, difference.lo.value, MPFR_RNDD);
        }
    }
    if (!certainly_positive(layer_monotonicity)) die(margins.label + " layer monotonicity failed");
    mpfr_set(margins.layer_monotonicity.value, layer_monotonicity.lo.value, MPFR_RNDD);

    Interval layer_sum = mul(weights[LAYERS - 1], cumulative_mass[LAYERS - 1]);
    for (int i = 0; i + 1 < LAYERS; ++i) {
        Interval weight_difference = sub(
            weights[static_cast<std::size_t>(i)],
            weights[static_cast<std::size_t>(i + 1)]
        );
        if (!certainly_positive(weight_difference)) die("layer weights are not decreasing");
        layer_sum = add(
            layer_sum,
            mul(weight_difference, cumulative_mass[static_cast<std::size_t>(i)])
        );
    }
    if (!certainly_positive(layer_sum)) die("layer sum is not positive");

    defining_tail = mul(
        saddle_det_lower,
        mul(exp_interval(neg(div_positive(square(saddle), exact_ui(2)))), layer_sum)
    );
    }

    const Interval q_minus_one = sub(q, exact_ui(1));
    const Interval square_lambda = mul(
        square_fraction,
        div_positive(q_minus_one, exact_ui(2))
    );
    const Interval square_tail = square_tail_upper(
        parameters.family,
        rank,
        q,
        square_lambda,
        factorials,
        margins.fredholm,
        margins.exact_square_pivot,
        exact_d_square
    );
    const Interval positive_probability = sub(scale_ui(defining_tail, 2), square_tail);
    if (!certainly_positive(positive_probability)) die(margins.label + " positive probability failed");
    mpfr_set(margins.positive_difference.value, positive_probability.lo.value, MPFR_RNDD);

    const Interval chebyshev_factor = one_minus(reciprocal_positive(square(character_cutoff)));
    if (!certainly_positive(chebyshev_factor)) die("character cutoff factor failed");
    const Interval push_lower = mul(square(c), mul(chebyshev_factor, positive_probability));

    const Interval a_minus_one = sub(a, exact_ui(1));
    const Interval negative_lambda = mul(
        negative_fraction,
        div_positive(a_minus_one, exact_ui(2))
    );
    const Interval negative_square_tail = square_tail_upper(
        parameters.family,
        rank,
        a,
        negative_lambda,
        factorials,
        margins.fredholm,
        margins.exact_square_pivot,
        exact_d_square
    );
    Interval negative_mass = mul(
        rational(8, 1),
        mul(
            exp_interval(neg(exact_ui(2))),
            div_positive(negative_square_tail, square(negative_lambda))
        )
    );

    const int k = parameters.polynomial_k;
    const int m = parameters.polynomial_m;
    if (parameters.family == 'D'
        && !exact_d_push_moments
        && 2 * k + 2 * m + 2 >= rank) {
        die(margins.label + " stable push-moment degree exceeds the D window");
    }
    Interval polynomial_expectation = exact_ui(0);
    for (int j = 0; j <= 2 * m; ++j) {
        Interval term = mul(
            exact_ui(binomial_ui(static_cast<unsigned>(2 * m), static_cast<unsigned>(j))),
            div_positive(
                exact_z(push_moments[static_cast<std::size_t>(2 * k + j)]),
                pow_positive(a, static_cast<unsigned long>(j))
            )
        );
        polynomial_expectation = (j % 2 == 0)
            ? add(polynomial_expectation, term)
            : sub(polynomial_expectation, term);
    }
    if (!certainly_positive(polynomial_expectation)) die("one-sided polynomial expectation failed");
    mpfr_set(margins.polynomial.value, polynomial_expectation.lo.value, MPFR_RNDD);

    const int d_first_hit = std::max(
        63,
        rank % 2 == 0 ? rank - 1 : rank - 2
    );
    const Interval support_ratio = div_positive(two_b, c);
    const Interval central_coefficient = div_positive(
        polynomial_expectation,
        mul(
            pow_positive(exact_ui(4), static_cast<unsigned long>(m)),
            pow_positive(a, static_cast<unsigned long>(2 * k))
        )
    );
    const Interval central_ratio = div_positive(a, c);
    int exponent = parameters.search_from_d_first_hit ? d_first_hit : 2 * rank + 29;
    Interval log_margin;
    for (;;) {
        const Interval negative_rhs = mul(
            negative_mass,
            pow_positive(support_ratio, static_cast<unsigned long>(exponent))
        );
        const Interval central_rhs = mul(
            central_coefficient,
            pow_positive(central_ratio, static_cast<unsigned long>(exponent))
        );
        const Interval rhs = add(negative_rhs, central_rhs);
        log_margin = sub(log_positive(push_lower), log_positive(rhs));
        if (certainly_positive(log_margin)) break;
        if (!parameters.search_from_d_first_hit || exponent >= 100000) {
            die(margins.label + " final pushforward comparison failed at rank " + std::to_string(rank));
        }
        exponent += 2;
    }
    if (enforce_d_a_free && parameters.family == 'D' && exponent + 2 > 2 * rank) {
        die(
            margins.label + " direct onset leaves the A-free D moment window: rank="
            + std::to_string(rank) + " onset=" + std::to_string(exponent)
            + " excess=" + std::to_string(exponent + 2 - 2 * rank)
        );
    }
    margins.exponent = exponent;
    mpfr_set(margins.final_margin.value, log_margin.lo.value, MPFR_RNDD);
    return margins;
}

struct ExactRowMargins {
    int rank = 0;
    Real geometry;
    Real cholesky_pivot;
    Real first_layer_mass;
    Real layer_monotonicity;
    Real positive_difference;
    Real polynomial;
    Real final_margin;
};

ExactLayeredTail exact_so_odd_layered_tail(
    int rank,
    int sign,
    const Interval& saddle,
    const Interval& lower_gap,
    const Interval& upper_start,
    const Interval& upper_step,
    const std::vector<mpz_class>& factorials,
    Real& minimum_pivot
) {
    const Interval log_mgf_s = log_so_odd_defining_mgf(
        rank, saddle, sign, factorials, minimum_pivot
    );
    const Interval lower_argument = sub(saddle, lower_gap);
    const Interval log_mgf_lower = log_so_odd_defining_mgf(
        rank, lower_argument, sign, factorials, minimum_pivot
    );
    const Interval lower_exponent = add(
        sub(mul(lower_gap, saddle), square(lower_gap)),
        sub(log_mgf_lower, log_mgf_s)
    );
    const Interval lower_tail = exp_interval(lower_exponent);

    std::array<Interval, LAYERS> cumulative_mass;
    std::array<Interval, LAYERS> weights;
    for (int i = 0; i < LAYERS; ++i) {
        const Interval upper_gap = add(
            upper_start,
            scale_ui(upper_step, static_cast<unsigned long>(i))
        );
        const Interval upper_argument = add(saddle, upper_gap);
        const Interval log_mgf_upper = log_so_odd_defining_mgf(
            rank, upper_argument, sign, factorials, minimum_pivot
        );
        const Interval upper_exponent = add(
            neg(add(mul(upper_gap, saddle), square(upper_gap))),
            sub(log_mgf_upper, log_mgf_s)
        );
        const Interval upper_tail = exp_interval(upper_exponent);
        cumulative_mass[static_cast<std::size_t>(i)] = sub(
            sub(exact_ui(1), lower_tail),
            upper_tail
        );
        if (!certainly_positive(cumulative_mass[static_cast<std::size_t>(i)])) {
            die("exact SO odd layered tilted mass failed");
        }
        weights[static_cast<std::size_t>(i)] = exp_interval(
            neg(mul(saddle, upper_gap))
        );
    }

    Interval monotonicity = sub(cumulative_mass[1], cumulative_mass[0]);
    for (int i = 1; i + 1 < LAYERS; ++i) {
        const Interval difference = sub(
            cumulative_mass[static_cast<std::size_t>(i + 1)],
            cumulative_mass[static_cast<std::size_t>(i)]
        );
        if (mpfr_cmp(difference.lo.value, monotonicity.lo.value) < 0) {
            mpfr_set(monotonicity.lo.value, difference.lo.value, MPFR_RNDD);
        }
    }
    if (!certainly_positive(monotonicity)) {
        die("exact SO odd layer monotonicity failed");
    }

    Interval layer_sum = mul(weights[LAYERS - 1], cumulative_mass[LAYERS - 1]);
    for (int i = 0; i + 1 < LAYERS; ++i) {
        const Interval weight_difference = sub(
            weights[static_cast<std::size_t>(i)],
            weights[static_cast<std::size_t>(i + 1)]
        );
        if (!certainly_positive(weight_difference)) {
            die("exact SO odd layer weights are not decreasing");
        }
        layer_sum = add(
            layer_sum,
            mul(weight_difference, cumulative_mass[static_cast<std::size_t>(i)])
        );
    }
    if (!certainly_positive(layer_sum)) die("exact SO odd layer sum failed");

    ExactLayeredTail out;
    out.probability = mul(
        exp_interval(sub(log_mgf_s, square(saddle))),
        layer_sum
    );
    out.first_mass = cumulative_mass[0];
    out.monotonicity = monotonicity;
    return out;
}

ExactRowMargins verify_exact_low_b_row(
    int rank,
    const std::vector<mpz_class>& factorials,
    const std::vector<mpz_class>& push_moments
) {
    if (rank != 14 && rank != 15) die("unexpected exact low-B rank");
    ExactRowMargins margins;
    margins.rank = rank;
    mpfr_set_ui(margins.cholesky_pivot.value, 1, MPFR_RNDD);

    const Interval b = exact_ui(static_cast<unsigned long>(rank));
    const Interval sqrt_b = sqrt_positive(b);
    const Interval r = rational(16, 5);
    const Interval square_alpha = rational(529, 100);
    const Interval negative_alpha = rational(213, 50);
    const Interval character_cutoff = rational(143, 100);
    const Interval lower_gap = rational(109, 100);
    const Interval upper_start = rational(13, 10);
    const Interval upper_step = rational(31, 1000);
    const Interval negative_fraction = rational(509, 500);

    const Interval c = mul(r, b);
    const Interval q = mul(square_alpha, sqrt_b);
    const Interval a = mul(negative_alpha, sqrt_b);
    const Interval two_b = scale_ui(b, 2);
    Interval geometry = sub(c, two_b);
    const Interval a_above_one = sub(a, exact_ui(1));
    const Interval a_below_support = sub(two_b, a);
    const Interval cutoff_above_one = sub(character_cutoff, exact_ui(1));
    for (const Interval* candidate : std::array<const Interval*, 3>{
             &a_above_one,
             &a_below_support,
             &cutoff_above_one
         }) {
        if (mpfr_cmp(candidate->lo.value, geometry.lo.value) < 0) {
            mpfr_set(geometry.lo.value, candidate->lo.value, MPFR_RNDD);
        }
    }
    if (!certainly_positive(geometry)) die("exact low-B geometry failed");
    mpfr_set(margins.geometry.value, geometry.lo.value, MPFR_RNDD);

    const Interval threshold = sqrt_positive(add(
        add(scale_ui(c, 2), scale_ui(character_cutoff, 2)),
        q
    ));
    const Interval saddle = add(threshold, lower_gap);
    ExactLayeredTail positive_tail = exact_so_odd_layered_tail(
        rank, 1, saddle, lower_gap, upper_start, upper_step, factorials,
        margins.cholesky_pivot
    );
    ExactLayeredTail negative_tail = exact_so_odd_layered_tail(
        rank, -1, saddle, lower_gap, upper_start, upper_step, factorials,
        margins.cholesky_pivot
    );
    Interval first_mass = positive_tail.first_mass;
    if (mpfr_cmp(negative_tail.first_mass.lo.value, first_mass.lo.value) < 0) {
        mpfr_set(first_mass.lo.value, negative_tail.first_mass.lo.value, MPFR_RNDD);
    }
    mpfr_set(margins.first_layer_mass.value, first_mass.lo.value, MPFR_RNDD);
    Interval layer_monotonicity = positive_tail.monotonicity;
    if (mpfr_cmp(
            negative_tail.monotonicity.lo.value,
            layer_monotonicity.lo.value
        ) < 0) {
        mpfr_set(
            layer_monotonicity.lo.value,
            negative_tail.monotonicity.lo.value,
            MPFR_RNDD
        );
    }
    mpfr_set(
        margins.layer_monotonicity.value,
        layer_monotonicity.lo.value,
        MPFR_RNDD
    );

    const Interval defining_tail = add(
        positive_tail.probability,
        negative_tail.probability
    );
    const Interval q_minus_one = sub(q, exact_ui(1));
    const Interval square_lambda = div_positive(q_minus_one, exact_ui(2));
    const Interval square_tail = exp_interval(add(
        neg(mul(square_lambda, q_minus_one)),
        square(square_lambda)
    ));
    const Interval positive_probability = sub(defining_tail, square_tail);
    if (!certainly_positive(positive_probability)) {
        die("exact low-B positive probability failed");
    }
    mpfr_set(
        margins.positive_difference.value,
        positive_probability.lo.value,
        MPFR_RNDD
    );
    const Interval chebyshev_factor = one_minus(reciprocal_positive(square(character_cutoff)));
    const Interval push_lower = mul(
        square(c),
        mul(chebyshev_factor, positive_probability)
    );

    const Interval a_minus_one = sub(a, exact_ui(1));
    const Interval negative_lambda = mul(
        negative_fraction,
        div_positive(a_minus_one, exact_ui(2))
    );
    const Interval negative_square_tail = exp_interval(add(
        neg(mul(negative_lambda, a_minus_one)),
        square(negative_lambda)
    ));
    const Interval negative_mass = mul(
        exact_ui(8),
        mul(
            exp_interval(neg(exact_ui(2))),
            div_positive(negative_square_tail, square(negative_lambda))
        )
    );

    constexpr int k = 5;
    const int m = rank - 6;
    Interval polynomial_expectation = exact_ui(0);
    for (int j = 0; j <= 2 * m; ++j) {
        const Interval term = mul(
            exact_ui(binomial_ui(static_cast<unsigned>(2 * m), static_cast<unsigned>(j))),
            div_positive(
                exact_z(push_moments[static_cast<std::size_t>(2 * k + j)]),
                pow_positive(a, static_cast<unsigned long>(j))
            )
        );
        polynomial_expectation = (j % 2 == 0)
            ? add(polynomial_expectation, term)
            : sub(polynomial_expectation, term);
    }
    if (!certainly_positive(polynomial_expectation)) {
        die("exact low-B one-sided polynomial expectation failed");
    }
    mpfr_set(margins.polynomial.value, polynomial_expectation.lo.value, MPFR_RNDD);

    const int exponent = 2 * rank + 29;
    const Interval negative_rhs = mul(
        negative_mass,
        pow_positive(div_positive(two_b, c), static_cast<unsigned long>(exponent))
    );
    const Interval central_coefficient = div_positive(
        polynomial_expectation,
        mul(
            pow_positive(exact_ui(4), static_cast<unsigned long>(m)),
            pow_positive(a, static_cast<unsigned long>(2 * k))
        )
    );
    const Interval central_rhs = mul(
        central_coefficient,
        pow_positive(div_positive(a, c), static_cast<unsigned long>(exponent))
    );
    const Interval rhs = add(negative_rhs, central_rhs);
    const Interval log_margin = sub(log_positive(push_lower), log_positive(rhs));
    if (!certainly_positive(log_margin)) {
        die("exact low-B final pushforward comparison failed");
    }
    mpfr_set(margins.final_margin.value, log_margin.lo.value, MPFR_RNDD);
    return margins;
}

}  // namespace

int main(int argc, char** argv) {
    int d_rank_low = 53;
    int d_rank_high = 295;
    bool d_only = false;
    bool allow_d_shallow_wall = false;
    bool exact_d_square = false;
    bool exact_d_defining = false;
    bool d_low_certificate = false;
    bool d12_24_exact_moment_certificate = false;
    std::vector<std::string> d_adjoint_moment_source_logs;
    int d_polynomial_k = 5;
    int d_polynomial_m = 8;
    RationalParameter d_square_alpha{441721, 100000};
    RationalParameter d_square_lambda_fraction{1, 1};
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--d-only") {
            d_only = true;
        } else if (argument == "--allow-d-shallow-wall") {
            allow_d_shallow_wall = true;
        } else if (argument == "--exact-d-square-mgf") {
            exact_d_square = true;
        } else if (argument == "--exact-d-defining-mgf") {
            exact_d_defining = true;
        } else if (argument == "--d-low-certificate") {
            d_low_certificate = true;
        } else if (argument == "--d12-24-exact-moment-certificate") {
            d12_24_exact_moment_certificate = true;
        } else if (argument == "--d-rank-low" && index + 1 < argc) {
            d_rank_low = std::stoi(argv[++index]);
        } else if (argument == "--d-rank-high" && index + 1 < argc) {
            d_rank_high = std::stoi(argv[++index]);
        } else if (argument == "--d-square-alpha" && index + 2 < argc) {
            d_square_alpha = {
                std::stoul(argv[++index]),
                std::stoul(argv[++index])
            };
        } else if (argument == "--d-square-lambda-fraction" && index + 2 < argc) {
            d_square_lambda_fraction = {
                std::stoul(argv[++index]),
                std::stoul(argv[++index])
            };
        } else if (argument == "--d-polynomial-k" && index + 1 < argc) {
            d_polynomial_k = std::stoi(argv[++index]);
        } else if (argument == "--d-polynomial-m" && index + 1 < argc) {
            d_polynomial_m = std::stoi(argv[++index]);
        } else if (argument == "--d-adjoint-moment-source-log"
                   && index + 1 < argc) {
            d_adjoint_moment_source_logs.emplace_back(argv[++index]);
        } else {
            die("unknown or incomplete command-line argument: " + argument);
        }
    }
    if (d_low_certificate) {
        d_only = true;
        d_rank_low = 25;
        d_rank_high = 52;
        allow_d_shallow_wall = true;
        exact_d_square = true;
        exact_d_defining = true;
        d_square_alpha = {441721, 100000};
        d_square_lambda_fraction = {1, 1};
    }
    if (d12_24_exact_moment_certificate) {
        d_only = true;
        d_rank_low = 12;
        d_rank_high = 24;
        allow_d_shallow_wall = true;
        exact_d_square = true;
        exact_d_defining = true;
        d_square_lambda_fraction = {1, 1};
        if (d_adjoint_moment_source_logs.empty()) {
            die("the D12-D24 exact-moment certificate requires adjoint-moment source logs");
        }
    }
    if (d_rank_low < 4 || d_rank_high < d_rank_low || d_rank_high > 295) {
        die("invalid D diagnostic rank range");
    }
    if (d_square_alpha.denominator == 0
        || d_square_lambda_fraction.denominator == 0) {
        die("zero denominator in D diagnostic parameter");
    }
    if (d_polynomial_k < 0 || d_polynomial_m < 0
        || 2 * d_polynomial_k + 2 * d_polynomial_m > 32) {
        die("invalid D diagnostic polynomial degrees");
    }

    constexpr int MAXIMUM = 2 * 295 + 2;
    std::vector<mpz_class> factorials(static_cast<std::size_t>(MAXIMUM + 1));
    factorials[0] = 1;
    for (int k = 1; k <= MAXIMUM; ++k) {
        factorials[static_cast<std::size_t>(k)] =
            factorials[static_cast<std::size_t>(k - 1)] * k;
    }
    const std::vector<mpz_class> push_moments = stable_push_moments(32);
    const std::vector<mpz_class> stable_adjoint = stable_adjoint_moments(34);
    if (push_moments_from_adjoint(stable_adjoint, 32) != push_moments) {
        die("stable adjoint-to-push moment reconstruction failed");
    }
    DAdjointMomentMap d_adjoint_moments;
    for (const std::string& path : d_adjoint_moment_source_logs) {
        load_d_adjoint_moment_source_log(path, d_adjoint_moments);
    }

    std::cout << "Post-m=29 layered B/C/D MGF certificate\n"
              << "mpfr_precision=" << PRECISION
              << " mpfr_version=" << mpfr_get_version()
              << " layers=" << LAYERS << '\n';
    std::cout << "exact_D_adjoint_moment_source_logs="
              << d_adjoint_moment_source_logs.size()
              << " exact_D_adjoint_moment_rows_loaded=" << d_adjoint_moments.size()
              << '\n';
    std::cout << "stable_push_moments";
    for (int q = 10; q <= 30; ++q) {
        std::cout << " K_" << q << '=' << push_moments[static_cast<std::size_t>(q)];
    }
    std::cout << '\n';

    Real minimum_geometry;
    Real minimum_fredholm;
    Real minimum_cholesky_pivot;
    Real minimum_exact_d_defining_pivot;
    Real minimum_exact_d_square_pivot;
    Real minimum_first_layer;
    Real minimum_layer_monotonicity;
    Real minimum_positive_difference;
    Real minimum_polynomial;
    Real minimum_final_margin;
    bool have_geometry = false;
    bool have_fredholm = false;
    bool have_cholesky_pivot = false;
    bool have_exact_d_defining_pivot = false;
    bool have_exact_d_square_pivot = false;
    bool have_first_layer = false;
    bool have_layer_monotonicity = false;
    bool have_positive_difference = false;
    bool have_polynomial = false;
    bool have_final_margin = false;
    std::string final_row;
    int rows = 0;
    int b_rows = 0;
    int c_rows = 0;
    int d_rows = 0;

    if (!d_only) {
      std::cout << "parameter_set=B14-15-exact-determinant family=B ranks=14..15"
                << " polynomial_k=5 polynomial_m=rank-6\n";
      for (int rank : {14, 15}) {
        ExactRowMargins row = verify_exact_low_b_row(rank, factorials, push_moments);
        ++rows;
        ++b_rows;
        update_min(minimum_geometry, row.geometry, have_geometry);
        update_min(
            minimum_cholesky_pivot,
            row.cholesky_pivot,
            have_cholesky_pivot
        );
        update_min(minimum_first_layer, row.first_layer_mass, have_first_layer);
        update_min(
            minimum_layer_monotonicity,
            row.layer_monotonicity,
            have_layer_monotonicity
        );
        update_min(
            minimum_positive_difference,
            row.positive_difference,
            have_positive_difference
        );
        update_min(minimum_polynomial, row.polynomial, have_polynomial);
        if (!have_final_margin
            || mpfr_cmp(row.final_margin.value, minimum_final_margin.value) < 0) {
            mpfr_set(minimum_final_margin.value, row.final_margin.value, MPFR_RNDD);
            have_final_margin = true;
            final_row = "B_" + std::to_string(rank);
        }
          std::cout << "row=B_" << rank
                    << " n=" << 2 * rank + 29
                    << " exact_toeplitz_hankel=1"
                    << " cholesky_pivot_lower=" << format_real(row.cholesky_pivot, 20)
                    << " final_log_margin_lower=" << format_real(row.final_margin, 20)
                    << '\n';
      }
    }

    for (const ParameterSet& stored_parameters : PARAMETER_SETS) {
        if (d_only && stored_parameters.family != 'D') continue;
        ParameterSet parameters = stored_parameters;
        if (parameters.family == 'D') {
            parameters.square_alpha = d_square_alpha;
            parameters.square_lambda_fraction = d_square_lambda_fraction;
            parameters.polynomial_k = d_polynomial_k;
            parameters.polynomial_m = d_polynomial_m;
        }
        const int rank_low = parameters.family == 'D'
            ? d_rank_low
            : parameters.rank_low;
        const int rank_high = parameters.family == 'D'
            ? d_rank_high
            : parameters.rank_high;
        std::cout << "parameter_set="
                  << (d12_24_exact_moment_certificate
                      ? "D12-24-exact-push-moment-onset"
                      : (d_low_certificate
                          ? "D25-52-exact-MGF-onset"
                          : parameters.label))
                  << " family=" << parameters.family
                  << " ranks=" << rank_low << ".." << rank_high;
        if (d12_24_exact_moment_certificate) {
            std::cout << " polynomial_schedule="
                      << "D12-16:(4,6),D17:(5,5),D18:(4,7),D19:(5,6),"
                      << "D20-21:(5,7),D22-23:(5,8),D24:(5,9)"
                      << " square_alpha_schedule="
                      << "D12-14:19/4,D15-16:9/2,D17-22:441721/100000,"
                      << "D23:17/4,D24:441721/100000"
                      << " defining_schedule=D12-24:exact"
                      << " square_schedule=D12-24:exact"
                      << " push_moment_schedule=D12-24:exact_from_adjoint_moments";
        } else if (d_low_certificate) {
            std::cout << " polynomial_schedule="
                      << "D25:(5,6),D26:(4,7),D27-28:(5,7),D29-52:(5,8)"
                      << " defining_schedule=D25-28:exact,D29-52:Fredholm"
                      << " square_schedule=D25-52:exact";
        } else {
            std::cout << " polynomial_k=" << parameters.polynomial_k
                      << " polynomial_m=" << parameters.polynomial_m;
        }
        std::cout << '\n';
        for (int rank = rank_low; rank <= rank_high; ++rank) {
            ParameterSet rank_parameters = parameters;
            if (d12_24_exact_moment_certificate) {
                if (rank <= 16) {
                    rank_parameters.polynomial_k = 4;
                    rank_parameters.polynomial_m = 6;
                    if (rank <= 14) {
                        rank_parameters.square_alpha = {19, 4};
                    } else {
                        rank_parameters.square_alpha = {9, 2};
                    }
                } else if (rank == 17) {
                    rank_parameters.polynomial_k = 5;
                    rank_parameters.polynomial_m = 5;
                    rank_parameters.square_alpha = {441721, 100000};
                } else if (rank == 18) {
                    rank_parameters.polynomial_k = 4;
                    rank_parameters.polynomial_m = 7;
                    rank_parameters.square_alpha = {441721, 100000};
                } else if (rank == 19) {
                    rank_parameters.polynomial_k = 5;
                    rank_parameters.polynomial_m = 6;
                    rank_parameters.square_alpha = {441721, 100000};
                } else if (rank <= 21) {
                    rank_parameters.polynomial_k = 5;
                    rank_parameters.polynomial_m = 7;
                    rank_parameters.square_alpha = {441721, 100000};
                } else if (rank <= 23) {
                    rank_parameters.polynomial_k = 5;
                    rank_parameters.polynomial_m = 8;
                    rank_parameters.square_alpha = rank == 23
                        ? RationalParameter{17, 4}
                        : RationalParameter{441721, 100000};
                } else {
                    rank_parameters.polynomial_k = 5;
                    rank_parameters.polynomial_m = 9;
                    rank_parameters.square_alpha = {441721, 100000};
                }
            } else if (d_low_certificate) {
                if (rank == 25) {
                    rank_parameters.polynomial_k = 5;
                    rank_parameters.polynomial_m = 6;
                } else if (rank == 26) {
                    rank_parameters.polynomial_k = 4;
                    rank_parameters.polynomial_m = 7;
                } else if (rank <= 28) {
                    rank_parameters.polynomial_k = 5;
                    rank_parameters.polynomial_m = 7;
                } else {
                    rank_parameters.polynomial_k = 5;
                    rank_parameters.polynomial_m = 8;
                }
            }
            const bool rank_exact_d_defining = exact_d_defining
                && (!d_low_certificate || rank <= 28);
            std::vector<mpz_class> rank_push_moments = push_moments;
            const bool rank_exact_d_push_moments =
                parameters.family == 'D' && !d_adjoint_moment_source_logs.empty();
            if (rank_exact_d_push_moments) {
                const int required_push_moment =
                    2 * rank_parameters.polynomial_k
                    + 2 * rank_parameters.polynomial_m;
                const int required_adjoint_moment = required_push_moment + 2;
                std::vector<mpz_class> rank_adjoint_moments = stable_adjoint;
                for (int moment = rank; moment <= required_adjoint_moment; ++moment) {
                    auto found = d_adjoint_moments.find({rank, moment});
                    if (found == d_adjoint_moments.end()) {
                        die(
                            "missing exact D adjoint moment: rank="
                            + std::to_string(rank)
                            + " moment=" + std::to_string(moment)
                        );
                    }
                    if (found->second < 0) {
                        die("negative exact D adjoint moment");
                    }
                    rank_adjoint_moments[static_cast<std::size_t>(moment)] =
                        found->second;
                }
                rank_push_moments = push_moments_from_adjoint(
                    rank_adjoint_moments,
                    required_push_moment
                );
            }
            RowMargins row = verify_row(
                rank_parameters,
                rank,
                factorials,
                rank_push_moments,
                !allow_d_shallow_wall,
                exact_d_square,
                rank_exact_d_defining,
                rank_exact_d_push_moments
            );
            if (d12_24_exact_moment_certificate && row.exponent != 63) {
                die(
                    "D12-D24 exact-moment onset is not 63: rank="
                    + std::to_string(rank)
                    + " onset=" + std::to_string(row.exponent)
                );
            }
            ++rows;
            if (parameters.family == 'B') {
                ++b_rows;
            } else if (parameters.family == 'C') {
                ++c_rows;
            } else {
                ++d_rows;
            }
            update_min(minimum_geometry, row.geometry, have_geometry);
            const bool row_uses_fredholm = parameters.family != 'D'
                || !exact_d_square
                || !rank_exact_d_defining;
            if (row_uses_fredholm) {
                update_min(minimum_fredholm, row.fredholm, have_fredholm);
            }
            if (parameters.family == 'D' && rank_exact_d_defining) {
                update_min(
                    minimum_exact_d_defining_pivot,
                    row.exact_defining_pivot,
                    have_exact_d_defining_pivot
                );
            }
            if (parameters.family == 'D' && exact_d_square) {
                update_min(
                    minimum_exact_d_square_pivot,
                    row.exact_square_pivot,
                    have_exact_d_square_pivot
                );
            }
            update_min(minimum_first_layer, row.first_layer_mass, have_first_layer);
            update_min(
                minimum_layer_monotonicity,
                row.layer_monotonicity,
                have_layer_monotonicity
            );
            update_min(
                minimum_positive_difference,
                row.positive_difference,
                have_positive_difference
            );
            update_min(minimum_polynomial, row.polynomial, have_polynomial);
            if (!have_final_margin
                || mpfr_cmp(row.final_margin.value, minimum_final_margin.value) < 0) {
                mpfr_set(minimum_final_margin.value, row.final_margin.value, MPFR_RNDD);
                have_final_margin = true;
                final_row = std::string(1, parameters.family) + "_" + std::to_string(rank);
            }
            std::cout << "row=" << parameters.family << '_' << rank
                      << " n=" << row.exponent
                      << " polynomial_k=" << rank_parameters.polynomial_k
                      << " polynomial_m=" << rank_parameters.polynomial_m
                      << (rank_exact_d_push_moments
                          ? " push_moments=reconstructed_from_exact_adjoint_source"
                          : " push_moments=stable")
                      << (parameters.family == 'D'
                          ? " A_free_slack=" + std::to_string(2 * rank - row.exponent - 2)
                          : std::string())
                      << (parameters.family == 'D' && rank_exact_d_defining
                          ? " exact_defining_pivot_lower="
                              + format_real(row.exact_defining_pivot, 20)
                          : std::string())
                      << (parameters.family == 'D' && exact_d_square
                          ? " exact_square_pivot_lower="
                              + format_real(row.exact_square_pivot, 20)
                          : std::string())
                      << " final_log_margin_lower=" << format_real(row.final_margin, 20)
                      << '\n';
        }
    }

    std::cout << "MINIMUM geometry=" << format_real(minimum_geometry)
              << " fredholm_trace_margin="
              << format_optional_real(minimum_fredholm, have_fredholm)
              << " exact_cholesky_pivot="
              << format_optional_real(
                     minimum_cholesky_pivot,
                     have_cholesky_pivot
                 )
              << " exact_D_defining_pivot="
              << format_optional_real(
                     minimum_exact_d_defining_pivot,
                     have_exact_d_defining_pivot
                 )
              << " exact_D_square_pivot="
              << format_optional_real(
                     minimum_exact_d_square_pivot,
                     have_exact_d_square_pivot
                 )
              << " first_layer_mass=" << format_real(minimum_first_layer)
              << " layer_monotonicity=" << format_real(minimum_layer_monotonicity)
              << " positive_difference=" << format_real(minimum_positive_difference)
              << " polynomial_expectation=" << format_real(minimum_polynomial)
              << " final_log_margin=" << format_real(minimum_final_margin)
              << " final_row=" << final_row << '\n';
    std::cout << "SUMMARY rows_checked=" << rows
              << " B_rows_checked=" << b_rows
              << " C_rows_checked=" << c_rows
              << " D_rows_checked=" << d_rows
              << " failures=0\n";
    std::cout << "__EXIT_STATUS=0\n";
    return 0;
}
