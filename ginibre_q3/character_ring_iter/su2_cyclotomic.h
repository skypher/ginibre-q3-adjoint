// Exact arithmetic in Z[zeta], zeta = 2cos(pi/n), for deciding tied
// certificate rows.
//
// Strict interval arithmetic cannot certify a row that holds with exact
// equality, and the proved-level calibration shows such ties are the entire
// engine residue at levels three and four.  Every spectral value in play is
// an algebraic integer in Z[zeta]: node values are Dickson polynomials
// D_m(zeta) (D_m(2cos t) = 2cos(mt)), characters are Chebyshev values
// S_a(D_{j+1}(zeta)), bases and their squares stay in the ring, and a
// domination row  prod_p lambda_p^{u_p} >= lambda_x^d  with integer weights
// is the sign of one ring element.  Zero means the row holds with equality,
// which weighted AM-GM permits.
//
// The minimal polynomial psi of zeta is derived exactly: Phi_{2n} is computed
// by integer cyclotomic division, and Phi_{2n}(z) = z^{D} Psi(z + 1/z) is
// inverted by integer back-substitution whose residual must vanish
// identically.  Nothing in the construction rests on floating point; interval
// evaluation is used only to determine the sign of provably nonzero elements,
// at precision raised until the enclosure excludes zero.

#pragma once

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

namespace cyclo {

using Poly = std::vector<mpz_class>;   // coefficient i of x^i; no trailing zeros

inline void trim(Poly& p) {
    while (!p.empty() && p.back() == 0) p.pop_back();
}

inline Poly mul(const Poly& a, const Poly& b) {
    if (a.empty() || b.empty()) return {};
    Poly r(a.size() + b.size() - 1, mpz_class(0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j)
            r[i + j] += a[i] * b[j];
    trim(r);
    return r;
}

// Exact division of integer polynomials; aborts if not exact.
inline Poly divExact(const Poly& num, const Poly& den) {
    Poly r = num, q;
    if (den.empty()) std::abort();
    q.assign(num.size() >= den.size() ? num.size() - den.size() + 1 : 0, mpz_class(0));
    while (r.size() >= den.size() && !r.empty()) {
        if (r.back() % den.back() != 0) std::abort();
        mpz_class c = r.back() / den.back();
        const std::size_t shift = r.size() - den.size();
        q[shift] = c;
        for (std::size_t i = 0; i < den.size(); ++i) r[shift + i] -= c * den[i];
        trim(r);
    }
    if (!r.empty()) std::abort();
    trim(q);
    return q;
}

// Cyclotomic polynomial Phi_m by the product formula, exactly.
inline Poly cyclotomic(int m) {
    Poly num(static_cast<std::size_t>(m) + 1, mpz_class(0));
    num[0] = -1;
    num[static_cast<std::size_t>(m)] = 1;                       // z^m - 1
    for (int d = 1; d < m; ++d)
        if (m % d == 0) num = divExact(num, cyclotomic(d));
    return num;
}

// The exact field data for one spectrum.
struct Field {
    int n = 0;                 // zeta = 2cos(pi/n)
    int D = 0;                 // degree of psi
    Poly psi;                  // monic minimal polynomial of zeta

    // An element is a Poly of degree < D, representing its value at zeta.
    using Elt = Poly;

    void init(int n_) {
        n = n_;
        Poly phi = cyclotomic(2 * n);
        D = static_cast<int>(phi.size()) - 1;
        if (D % 2 != 0) std::abort();
        D /= 2;
        // Invert Phi(z) = z^D Psi(z + 1/z): peel a_i from the top and
        // subtract a_i * z^{D-i} (z^2+1)^i; the residual must vanish.
        Poly residual = phi;
        psi.assign(static_cast<std::size_t>(D) + 1, mpz_class(0));
        Poly zsq1 = {mpz_class(1), mpz_class(0), mpz_class(1)};   // z^2 + 1
        for (int i = D; i >= 0; --i) {
            const std::size_t top = static_cast<std::size_t>(D + i);
            mpz_class a = (top < residual.size()) ? residual[top] : mpz_class(0);
            psi[static_cast<std::size_t>(i)] = a;
            if (a == 0) continue;
            Poly term = {mpz_class(1)};
            for (int e = 0; e < i; ++e) term = mul(term, zsq1);
            Poly shifted(static_cast<std::size_t>(D - i), mpz_class(0));
            shifted.insert(shifted.end(), term.begin(), term.end());
            if (shifted.size() > residual.size()) residual.resize(shifted.size(), mpz_class(0));
            for (std::size_t j = 0; j < shifted.size(); ++j) residual[j] -= a * shifted[j];
            trim(residual);
        }
        if (!residual.empty()) {
            std::fprintf(stderr, "cyclotomic: Psi extraction residual nonzero for n=%d\n", n);
            std::abort();
        }
        if (psi.back() != 1) {
            std::fprintf(stderr, "cyclotomic: psi not monic for n=%d\n", n);
            std::abort();
        }
    }

    Elt fromInt(long v) const { Elt e; if (v) e.push_back(mpz_class(v)); return e; }

    Elt reduce(Poly p) const {
        while (static_cast<int>(p.size()) > D) {
            const mpz_class c = p.back();          // psi is monic
            const std::size_t shift = p.size() - psi.size();
            for (std::size_t i = 0; i < psi.size(); ++i) p[shift + i] -= c * psi[i];
            trim(p);
        }
        return p;
    }

    Elt add(const Elt& a, const Elt& b) const {
        Elt r = a;
        if (b.size() > r.size()) r.resize(b.size(), mpz_class(0));
        for (std::size_t i = 0; i < b.size(); ++i) r[i] += b[i];
        trim(r);
        return r;
    }
    Elt sub(const Elt& a, const Elt& b) const {
        Elt r = a;
        if (b.size() > r.size()) r.resize(b.size(), mpz_class(0));
        for (std::size_t i = 0; i < b.size(); ++i) r[i] -= b[i];
        trim(r);
        return r;
    }
    Elt mulE(const Elt& a, const Elt& b) const { return reduce(mul(a, b)); }
    Elt powE(Elt base, unsigned long e) const {
        Elt r = fromInt(1);
        while (e) {
            if (e & 1UL) r = mulE(r, base);
            base = mulE(base, base);
            e >>= 1;
        }
        return r;
    }
    static bool isZero(const Elt& a) { return a.empty(); }

    // Dickson value D_m(zeta) = 2cos(m pi / n), by the recurrence.
    Elt dickson(int m) const {
        Elt d0 = fromInt(2), d1 = {mpz_class(0), mpz_class(1)};   // zeta itself
        if (m == 0) return d0;
        for (int i = 1; i < m; ++i) {
            Elt d2 = sub(mulE(d1, {mpz_class(0), mpz_class(1)}), d0);
            d0 = d1;
            d1 = d2;
        }
        return d1;
    }

    // Chebyshev value S_a(y) with S_0 = 1, S_1 = y: this is
    // U_a(cos t) when y = 2cos t.
    Elt chebyshevS(int a, const Elt& y) const {
        Elt s0 = fromInt(1), s1 = y;
        if (a == 0) return s0;
        for (int i = 1; i < a; ++i) {
            Elt s2 = sub(mulE(y, s1), s0);
            s0 = s1;
            s1 = s2;
        }
        return s1;
    }

    // Outward enclosure of zeta = 2 cos(pi/n).  Since cos is decreasing on
    // [0,pi/3] for n>=3, the lower zeta endpoint comes from the upper angle
    // endpoint and the upper zeta endpoint comes from the lower angle
    // endpoint.  The rounding directions must follow that reversal; computing
    // cos(theta_lo) downward and cos(theta_hi) upward and merely swapping the
    // results is not an outward enclosure.
    void rootInterval(mpfr_ptr zlo, mpfr_ptr zhi) const {
        if (n < 3 || mpfr_get_prec(zlo) != mpfr_get_prec(zhi)) std::abort();
        const mpfr_prec_t prec = mpfr_get_prec(zlo);
        mpfr_t theta_lo, theta_hi;
        mpfr_init2(theta_lo, prec);
        mpfr_init2(theta_hi, prec);
        mpfr_const_pi(theta_lo, MPFR_RNDD);
        mpfr_div_si(theta_lo, theta_lo, n, MPFR_RNDD);
        mpfr_const_pi(theta_hi, MPFR_RNDU);
        mpfr_div_si(theta_hi, theta_hi, n, MPFR_RNDU);
        mpfr_cos(zlo, theta_hi, MPFR_RNDD);
        mpfr_mul_ui(zlo, zlo, 2, MPFR_RNDD);
        mpfr_cos(zhi, theta_lo, MPFR_RNDU);
        mpfr_mul_ui(zhi, zhi, 2, MPFR_RNDU);
        mpfr_clear(theta_lo);
        mpfr_clear(theta_hi);
        if (mpfr_greater_p(zlo, zhi)) std::abort();
    }

    // Sign of a nonzero element's value at zeta, by interval evaluation at
    // rising precision.  Returns 0 only for the zero element.
    int sign(const Elt& a) const {
        if (isZero(a)) return 0;
        // The needed precision scales with the element's height: a
        // denominator-5000 escalation raises lambda to the 5000th power and
        // its coefficients to tens of thousands of bits, and resolving a
        // near-cancelling difference requires precision above that height.
        // A nonzero algebraic number of degree D and height H cannot be
        // smaller than roughly H^-D, so D*height plus slack always suffices.
        std::size_t hbits = 0;
        for (const auto& c : a)
            hbits = std::max(hbits, mpz_sizeinbase(c.get_mpz_t(), 2));
        const mpfr_prec_t cap = static_cast<mpfr_prec_t>(
            static_cast<std::size_t>(D) * (hbits + 64) + 4096);
        for (mpfr_prec_t prec = 512; prec <= cap; prec *= 2) {
            mpfr_t zlo, zhi, plo, phi_, t;
            mpfr_init2(zlo, prec); mpfr_init2(zhi, prec);
            mpfr_init2(plo, prec); mpfr_init2(phi_, prec); mpfr_init2(t, prec);
            rootInterval(zlo, zhi);
            // Horner in interval arithmetic over [zlo, zhi] (both positive
            // for n >= 3, but keep it general with four-corner products).
            mpfr_set_ui(plo, 0, MPFR_RNDD);
            mpfr_set_ui(phi_, 0, MPFR_RNDU);
            for (std::size_t i = a.size(); i-- > 0;) {
                // [plo,phi] = [plo,phi] * [zlo,zhi] + a[i]
                mpfr_t c00, c01, c10, c11;
                mpfr_init2(c00, prec); mpfr_init2(c01, prec);
                mpfr_init2(c10, prec); mpfr_init2(c11, prec);
                mpfr_mul(c00, plo, zlo, MPFR_RNDD);
                mpfr_mul(c01, plo, zhi, MPFR_RNDD);
                mpfr_mul(c10, phi_, zlo, MPFR_RNDD);
                mpfr_mul(c11, phi_, zhi, MPFR_RNDD);
                mpfr_min(t, c00, c01, MPFR_RNDD);
                mpfr_min(t, t, c10, MPFR_RNDD);
                mpfr_min(t, t, c11, MPFR_RNDD);
                mpfr_t nlo; mpfr_init2(nlo, prec); mpfr_set(nlo, t, MPFR_RNDD);
                mpfr_mul(c00, plo, zlo, MPFR_RNDU);
                mpfr_mul(c01, plo, zhi, MPFR_RNDU);
                mpfr_mul(c10, phi_, zlo, MPFR_RNDU);
                mpfr_mul(c11, phi_, zhi, MPFR_RNDU);
                mpfr_max(t, c00, c01, MPFR_RNDU);
                mpfr_max(t, t, c10, MPFR_RNDU);
                mpfr_max(t, t, c11, MPFR_RNDU);
                mpfr_set(phi_, t, MPFR_RNDU);
                mpfr_set(plo, nlo, MPFR_RNDD);
                mpfr_clear(nlo);
                mpfr_clear(c00); mpfr_clear(c01); mpfr_clear(c10); mpfr_clear(c11);
                mpfr_add_z(plo, plo, a[i].get_mpz_t(), MPFR_RNDD);
                mpfr_add_z(phi_, phi_, a[i].get_mpz_t(), MPFR_RNDU);
            }
            const int slo = mpfr_sgn(plo), shi = mpfr_sgn(phi_);
            mpfr_clear(zlo); mpfr_clear(zhi); mpfr_clear(plo); mpfr_clear(phi_); mpfr_clear(t);
            if (slo > 0) return 1;
            if (shi < 0) return -1;
        }
        std::fprintf(stderr, "cyclotomic: sign undecided at height-scaled cap (element nonzero)\n");
        std::abort();
    }
};

}  // namespace cyclo
