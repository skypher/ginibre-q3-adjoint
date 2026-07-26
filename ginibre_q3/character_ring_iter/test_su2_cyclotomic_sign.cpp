// Independent regression for the exact cyclotomic sign boundary used by the
// finite SU(2)_k AM-GM verifier.
//
// The test checks, for n=3..32:
//   1. the outward interval for zeta=2cos(pi/n) contains an independently
//      evaluated 4096-bit value;
//   2. interval Horner evaluation of the exact minimal polynomial contains 0;
//   3. Field::sign agrees with independent 4096-bit direct evaluation on 200
//      deterministic nonzero field elements.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include "su2_cyclotomic.h"

namespace {

void mul_interval(mpfr_ptr out_lo, mpfr_ptr out_hi,
                  mpfr_srcptr a_lo, mpfr_srcptr a_hi,
                  mpfr_srcptr b_lo, mpfr_srcptr b_hi) {
    const mpfr_prec_t prec = mpfr_get_prec(out_lo);
    mpfr_t down[4], up[4];
    for (int i = 0; i < 4; ++i) {
        mpfr_init2(down[i], prec);
        mpfr_init2(up[i], prec);
    }
    mpfr_mul(down[0], a_lo, b_lo, MPFR_RNDD);
    mpfr_mul(down[1], a_lo, b_hi, MPFR_RNDD);
    mpfr_mul(down[2], a_hi, b_lo, MPFR_RNDD);
    mpfr_mul(down[3], a_hi, b_hi, MPFR_RNDD);
    mpfr_set(out_lo, down[0], MPFR_RNDD);
    for (int i = 1; i < 4; ++i)
        if (mpfr_less_p(down[i], out_lo)) mpfr_set(out_lo, down[i], MPFR_RNDD);

    mpfr_mul(up[0], a_lo, b_lo, MPFR_RNDU);
    mpfr_mul(up[1], a_lo, b_hi, MPFR_RNDU);
    mpfr_mul(up[2], a_hi, b_lo, MPFR_RNDU);
    mpfr_mul(up[3], a_hi, b_hi, MPFR_RNDU);
    mpfr_set(out_hi, up[0], MPFR_RNDU);
    for (int i = 1; i < 4; ++i)
        if (mpfr_greater_p(up[i], out_hi)) mpfr_set(out_hi, up[i], MPFR_RNDU);

    for (int i = 0; i < 4; ++i) {
        mpfr_clear(down[i]);
        mpfr_clear(up[i]);
    }
}

void eval_interval(const cyclo::Poly& p, mpfr_srcptr zlo, mpfr_srcptr zhi,
                   mpfr_ptr out_lo, mpfr_ptr out_hi) {
    const mpfr_prec_t prec = mpfr_get_prec(out_lo);
    mpfr_set_zero(out_lo, 1);
    mpfr_set_zero(out_hi, 1);
    mpfr_t next_lo, next_hi;
    mpfr_init2(next_lo, prec);
    mpfr_init2(next_hi, prec);
    for (std::size_t i = p.size(); i-- > 0;) {
        mul_interval(next_lo, next_hi, out_lo, out_hi, zlo, zhi);
        mpfr_add_z(next_lo, next_lo, p[i].get_mpz_t(), MPFR_RNDD);
        mpfr_add_z(next_hi, next_hi, p[i].get_mpz_t(), MPFR_RNDU);
        mpfr_set(out_lo, next_lo, MPFR_RNDD);
        mpfr_set(out_hi, next_hi, MPFR_RNDU);
    }
    mpfr_clear(next_lo);
    mpfr_clear(next_hi);
}

void direct_root(mpfr_ptr z, int n) {
    mpfr_const_pi(z, MPFR_RNDN);
    mpfr_div_si(z, z, n, MPFR_RNDN);
    mpfr_cos(z, z, MPFR_RNDN);
    mpfr_mul_ui(z, z, 2, MPFR_RNDN);
}

// The pre-audit construction, retained only as a negative control: it rounds
// cos(theta_lo) downward and cos(theta_hi) upward before swapping.  Those are
// inward directions after the monotonicity reversal.
void legacy_inward_root(mpfr_ptr zlo, mpfr_ptr zhi, int n) {
    mpfr_const_pi(zlo, MPFR_RNDD);
    mpfr_div_si(zlo, zlo, n, MPFR_RNDD);
    mpfr_cos(zlo, zlo, MPFR_RNDD);
    mpfr_mul_ui(zlo, zlo, 2, MPFR_RNDD);
    mpfr_const_pi(zhi, MPFR_RNDU);
    mpfr_div_si(zhi, zhi, n, MPFR_RNDU);
    mpfr_cos(zhi, zhi, MPFR_RNDU);
    mpfr_mul_ui(zhi, zhi, 2, MPFR_RNDU);
    if (mpfr_greater_p(zlo, zhi)) mpfr_swap(zlo, zhi);
}

int direct_sign(const cyclo::Field::Elt& a, int n) {
    constexpr mpfr_prec_t prec = 4096;
    mpfr_t z, value;
    mpfr_init2(z, prec);
    mpfr_init2(value, prec);
    direct_root(z, n);
    mpfr_set_zero(value, 1);
    for (std::size_t i = a.size(); i-- > 0;) {
        mpfr_mul(value, value, z, MPFR_RNDN);
        mpfr_add_z(value, value, a[i].get_mpz_t(), MPFR_RNDN);
    }
    const int result = mpfr_sgn(value);
    mpfr_clear(z);
    mpfr_clear(value);
    return result;
}

std::uint64_t next_state(std::uint64_t& state) {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

}  // namespace

int main() {
    constexpr mpfr_prec_t interval_prec = 512;
    constexpr mpfr_prec_t reference_prec = 4096;
    std::uint64_t checked_elements = 0;
    int rejected_legacy_intervals = 0;

    for (int n = 3; n <= 32; ++n) {
        cyclo::Field field;
        field.init(n);

        mpfr_t zlo, zhi, zref, psi_lo, psi_hi;
        mpfr_init2(zlo, interval_prec);
        mpfr_init2(zhi, interval_prec);
        mpfr_init2(zref, reference_prec);
        mpfr_init2(psi_lo, interval_prec);
        mpfr_init2(psi_hi, interval_prec);
        field.rootInterval(zlo, zhi);
        direct_root(zref, n);

        if (mpfr_less_p(zref, zlo) || mpfr_greater_p(zref, zhi)) {
            std::fprintf(stderr, "root enclosure missed reference at n=%d\n", n);
            return 1;
        }
        legacy_inward_root(zlo, zhi, n);
        if (mpfr_less_p(zref, zlo) || mpfr_greater_p(zref, zhi))
            ++rejected_legacy_intervals;
        field.rootInterval(zlo, zhi);
        eval_interval(field.psi, zlo, zhi, psi_lo, psi_hi);
        if (mpfr_sgn(psi_lo) > 0 || mpfr_sgn(psi_hi) < 0) {
            std::fprintf(stderr, "minimal polynomial interval missed zero at n=%d\n", n);
            return 1;
        }
        if (field.sign(field.reduce(field.psi)) != 0) {
            std::fprintf(stderr, "minimal polynomial did not reduce to zero at n=%d\n", n);
            return 1;
        }

        std::uint64_t state = UINT64_C(0x9e3779b97f4a7c15) ^
                              static_cast<std::uint64_t>(n);
        for (int trial = 0; trial < 200; ++trial) {
            cyclo::Field::Elt a(static_cast<std::size_t>(field.D), mpz_class(0));
            for (int i = 0; i < field.D; ++i) {
                const long coefficient =
                    static_cast<long>(next_state(state) % UINT64_C(2000001)) - 1000000L;
                a[static_cast<std::size_t>(i)] = coefficient;
            }
            if (a.back() == 0) a.back() = 1;
            cyclo::trim(a);
            const int exact = field.sign(a);
            const int reference = direct_sign(a, n);
            if (reference == 0 || exact != reference) {
                std::fprintf(stderr,
                             "sign disagreement at n=%d trial=%d exact=%d reference=%d\n",
                             n, trial, exact, reference);
                return 1;
            }
            ++checked_elements;
        }

        mpfr_clear(zlo);
        mpfr_clear(zhi);
        mpfr_clear(zref);
        mpfr_clear(psi_lo);
        mpfr_clear(psi_hi);
    }

    if (rejected_legacy_intervals == 0) {
        std::fprintf(stderr, "negative control failed to reject the legacy root interval\n");
        return 1;
    }
    std::printf("SU2_CYCLOTOMIC_SIGN levels=30 elements=%llu legacy_rejected=%d result=PASS\n",
                static_cast<unsigned long long>(checked_elements),
                rejected_legacy_intervals);
    return 0;
}
