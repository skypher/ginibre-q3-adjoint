// Rank-parameterised AM-GM certificate stack for the odd simple-current orbit
// rings, in C++.
//
// This is a port of verify_su2_orbit_amgm_stack.py.  The Python version is
// correct but slow: a full rank-six pass takes about eight minutes and a full
// rank-seven pass would take roughly fifteen hours, which makes measurement
// the bottleneck rather than the mathematics.
//
// Trust boundary is unchanged from the Python.  Long double proposes, via a
// linear program; nothing it produces is accepted.  Every certificate is
// re-checked in MPFR interval arithmetic with directed rounding, on spectral
// data rebuilt independently in interval mode, and the structural conditions
// are checked exactly over the rationals.
//
// Build (see Makefile.research target verify_su2_orbit_amgm_stack):
//   g++ -O2 -std=c++20 -Wall -Wextra -Wpedantic -Wconversion
//       -Wsign-conversion -Wshadow -Werror -pthread
//       verify_su2_orbit_amgm_stack.cpp -l:libmpfr.so.6 -lgmp
//       -o verify_su2_orbit_amgm_stack
// Run:
//   ./verify_su2_orbit_amgm_stack --rank 6 [--decompose] [--threads N]

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

#include <gmp.h>
#include <mpfr.h>

namespace {

constexpr mpfr_prec_t kPrec = 512;
constexpr long double kHallTol = 1.0e-12L;

// ---------------------------------------------------------------- intervals

// Minimal outward-rounded interval.  Every operation rounds the lower bound
// down and the upper bound up, so an enclosure is never optimistic.
struct Iv {
    mpfr_t lo, hi;
    Iv() { mpfr_init2(lo, kPrec); mpfr_init2(hi, kPrec); mpfr_set_zero(lo, 1); mpfr_set_zero(hi, 1); }
    Iv(const Iv& o) { mpfr_init2(lo, kPrec); mpfr_init2(hi, kPrec); mpfr_set(lo, o.lo, MPFR_RNDD); mpfr_set(hi, o.hi, MPFR_RNDU); }
    Iv& operator=(const Iv& o) { if (this != &o) { mpfr_set(lo, o.lo, MPFR_RNDD); mpfr_set(hi, o.hi, MPFR_RNDU); } return *this; }
    ~Iv() { mpfr_clear(lo); mpfr_clear(hi); }
    void set_ui(unsigned long v) { mpfr_set_ui(lo, v, MPFR_RNDD); mpfr_set_ui(hi, v, MPFR_RNDU); }
};

void iv_add(Iv& r, const Iv& a, const Iv& b) {
    mpfr_add(r.lo, a.lo, b.lo, MPFR_RNDD);
    mpfr_add(r.hi, a.hi, b.hi, MPFR_RNDU);
}
void iv_sub(Iv& r, const Iv& a, const Iv& b) {
    mpfr_sub(r.lo, a.lo, b.hi, MPFR_RNDD);
    mpfr_sub(r.hi, a.hi, b.lo, MPFR_RNDU);
}

// General product; the operands here are not sign-constant so all four
// corner products are considered.
void iv_mul(Iv& r, const Iv& a, const Iv& b) {
    mpfr_t c[4];
    for (int i = 0; i < 4; ++i) mpfr_init2(c[i], kPrec);
    mpfr_mul(c[0], a.lo, b.lo, MPFR_RNDD);
    mpfr_mul(c[1], a.lo, b.hi, MPFR_RNDD);
    mpfr_mul(c[2], a.hi, b.lo, MPFR_RNDD);
    mpfr_mul(c[3], a.hi, b.hi, MPFR_RNDD);
    mpfr_t lo; mpfr_init2(lo, kPrec); mpfr_set(lo, c[0], MPFR_RNDD);
    for (int i = 1; i < 4; ++i) if (mpfr_less_p(c[i], lo)) mpfr_set(lo, c[i], MPFR_RNDD);
    mpfr_mul(c[0], a.lo, b.lo, MPFR_RNDU);
    mpfr_mul(c[1], a.lo, b.hi, MPFR_RNDU);
    mpfr_mul(c[2], a.hi, b.lo, MPFR_RNDU);
    mpfr_mul(c[3], a.hi, b.hi, MPFR_RNDU);
    mpfr_t hi; mpfr_init2(hi, kPrec); mpfr_set(hi, c[0], MPFR_RNDU);
    for (int i = 1; i < 4; ++i) if (mpfr_greater_p(c[i], hi)) mpfr_set(hi, c[i], MPFR_RNDU);
    mpfr_set(r.lo, lo, MPFR_RNDD);
    mpfr_set(r.hi, hi, MPFR_RNDU);
    mpfr_clear(lo); mpfr_clear(hi);
    for (int i = 0; i < 4; ++i) mpfr_clear(c[i]);
}

// log of a strictly positive interval.
bool iv_log(Iv& r, const Iv& a) {
    if (mpfr_sgn(a.lo) <= 0) return false;
    mpfr_log(r.lo, a.lo, MPFR_RNDD);
    mpfr_log(r.hi, a.hi, MPFR_RNDU);
    return true;
}

// r += a * (num/den), with num/den a nonnegative rational.
void iv_addmul_q(Iv& r, const Iv& a, long num, long den) {
    Iv t;
    mpfr_mul_si(t.lo, a.lo, num, MPFR_RNDD);
    mpfr_div_si(t.lo, t.lo, den, MPFR_RNDD);
    mpfr_mul_si(t.hi, a.hi, num, MPFR_RNDU);
    mpfr_div_si(t.hi, t.hi, den, MPFR_RNDU);
    Iv s; iv_add(s, r, t); r = s;
}

// ---------------------------------------------------------------- spectral

struct Model {
    int rank = 0, labels = 0, order = 0;
    std::vector<long double> w;                       // node weights
    std::vector<std::vector<long double>> v;          // character table
    std::vector<Iv> iw;
    std::vector<std::vector<Iv>> iv_tab;
};

void build_model(Model& m, int rank) {
    m.rank = rank;
    m.labels = rank - 1;
    m.order = 2 * rank + 1;
    const long double pi = acosl(-1.0L);
    m.w.assign(static_cast<std::size_t>(rank), 0.0L);
    m.v.assign(static_cast<std::size_t>(rank), std::vector<long double>(static_cast<std::size_t>(m.labels), 0.0L));
    for (int t = 1; t <= rank; ++t) {
        const long double a = pi * static_cast<long double>(t) / static_cast<long double>(m.order);
        const std::size_t i = static_cast<std::size_t>(t - 1);
        m.w[i] = 4.0L * sinl(a) * sinl(a) / static_cast<long double>(m.order);
        for (int j = 1; j <= m.labels; ++j)
            m.v[i][static_cast<std::size_t>(j - 1)] =
                sinl(static_cast<long double>(2 * j + 1) * a) / sinl(a);
    }
    // Independent interval rebuild: nodes as sin(pi t / n) with interval pi.
    m.iw.resize(static_cast<std::size_t>(rank));
    m.iv_tab.assign(static_cast<std::size_t>(rank), std::vector<Iv>(static_cast<std::size_t>(m.labels)));
    Iv ivpi;
    mpfr_const_pi(ivpi.lo, MPFR_RNDD);
    mpfr_const_pi(ivpi.hi, MPFR_RNDU);
    for (int t = 1; t <= rank; ++t) {
        const std::size_t i = static_cast<std::size_t>(t - 1);
        Iv a;
        mpfr_mul_si(a.lo, ivpi.lo, t, MPFR_RNDD);
        mpfr_div_si(a.lo, a.lo, m.order, MPFR_RNDD);
        mpfr_mul_si(a.hi, ivpi.hi, t, MPFR_RNDU);
        mpfr_div_si(a.hi, a.hi, m.order, MPFR_RNDU);
        Iv s;                                  // sin is increasing on (0, pi/2]
        mpfr_sin(s.lo, a.lo, MPFR_RNDD);
        mpfr_sin(s.hi, a.hi, MPFR_RNDU);
        if (mpfr_cmp_ui(a.hi, 0) > 0 && mpfr_greater_p(a.hi, ivpi.hi)) { /* unreachable */ }
        // For t up to rank, a <= pi*rank/(2rank+1) < pi/2, so sin is monotone.
        Iv s2; iv_mul(s2, s, s);
        mpfr_mul_ui(m.iw[i].lo, s2.lo, 4, MPFR_RNDD);
        mpfr_div_si(m.iw[i].lo, m.iw[i].lo, m.order, MPFR_RNDD);
        mpfr_mul_ui(m.iw[i].hi, s2.hi, 4, MPFR_RNDU);
        mpfr_div_si(m.iw[i].hi, m.iw[i].hi, m.order, MPFR_RNDU);
        for (int j = 1; j <= m.labels; ++j) {
            Iv na;                            // (2j+1) a
            mpfr_mul_si(na.lo, a.lo, 2 * j + 1, MPFR_RNDD);
            mpfr_mul_si(na.hi, a.hi, 2 * j + 1, MPFR_RNDU);
            Iv sn;                            // sin over a possibly wide range
            mpfr_t c0, c1;
            mpfr_init2(c0, kPrec); mpfr_init2(c1, kPrec);
            mpfr_sin(c0, na.lo, MPFR_RNDD);
            mpfr_sin(c1, na.hi, MPFR_RNDD);
            mpfr_set(sn.lo, mpfr_less_p(c0, c1) ? c0 : c1, MPFR_RNDD);
            mpfr_sin(c0, na.lo, MPFR_RNDU);
            mpfr_sin(c1, na.hi, MPFR_RNDU);
            mpfr_set(sn.hi, mpfr_greater_p(c0, c1) ? c0 : c1, MPFR_RNDU);
            mpfr_clear(c0); mpfr_clear(c1);
            // Widen by the possible turning point: the enclosure above is only
            // valid when sin is monotone across [na.lo, na.hi].  The interval
            // is extremely narrow here (pi to 512 bits), so verify monotonicity
            // by checking the cosine keeps one sign; otherwise widen to [-1,1]
            // scaled, which simply fails the certificate rather than lying.
            mpfr_t clo, chi;
            mpfr_init2(clo, kPrec); mpfr_init2(chi, kPrec);
            mpfr_cos(clo, na.lo, MPFR_RNDD);
            mpfr_cos(chi, na.hi, MPFR_RNDU);
            if (mpfr_sgn(clo) * mpfr_sgn(chi) < 0) {
                mpfr_set_si(sn.lo, -1, MPFR_RNDD);
                mpfr_set_si(sn.hi, 1, MPFR_RNDU);
            }
            mpfr_clear(clo); mpfr_clear(chi);
            // chi_(2j) = sin((2j+1)a)/sin(a)
            Iv q;
            mpfr_div(q.lo, sn.lo, s.hi, MPFR_RNDD);
            mpfr_div(q.hi, sn.hi, s.lo, MPFR_RNDU);
            if (mpfr_sgn(sn.lo) < 0) {            // negative numerator: swap divisors
                mpfr_div(q.lo, sn.lo, s.lo, MPFR_RNDD);
                mpfr_div(q.hi, sn.hi, s.hi, MPFR_RNDU);
                if (mpfr_sgn(sn.hi) > 0) {        // straddles zero
                    mpfr_div(q.lo, sn.lo, s.lo, MPFR_RNDD);
                    mpfr_div(q.hi, sn.hi, s.lo, MPFR_RNDU);
                }
            }
            m.iv_tab[i][static_cast<std::size_t>(j - 1)] = q;
        }
    }
}

// ---------------------------------------------------------------- terms

struct Term {
    int sign = 1;
    long double coeff = 0.0L;
    std::vector<long double> lam;
    Iv icoeff;
    std::vector<Iv> ilam;
};

// Build the off-diagonal term list for a chamber, in both arithmetics.
std::vector<Term> chamber_terms(const Model& m, const std::vector<int>& signs,
                                const std::vector<int>& powers) {
    std::vector<Term> out;
    const std::size_t L = static_cast<std::size_t>(m.labels);
    for (int i = 0; i < m.rank; ++i) {
        for (int j = i + 1; j < m.rank; ++j) {
            Term t;
            t.sign = 1;
            t.coeff = 2.0L * m.w[static_cast<std::size_t>(i)] * m.w[static_cast<std::size_t>(j)];
            t.lam.assign(L, 1.0L);
            t.ilam.assign(L, Iv());
            for (std::size_t l = 0; l < L; ++l) t.ilam[l].set_ui(1);
            Iv acc;
            iv_mul(acc, m.iw[static_cast<std::size_t>(i)], m.iw[static_cast<std::size_t>(j)]);
            mpfr_mul_ui(acc.lo, acc.lo, 2, MPFR_RNDD);
            mpfr_mul_ui(acc.hi, acc.hi, 2, MPFR_RNDU);
            bool degenerate = false;
            for (std::size_t l = 0; l < L; ++l) {
                if (signs[l] == 0) continue;
                const long double b = m.v[static_cast<std::size_t>(i)][l] +
                    static_cast<long double>(signs[l]) * m.v[static_cast<std::size_t>(j)][l];
                if (fabsl(b) < 1.0e-30L) { degenerate = true; break; }
                Iv ib;
                if (signs[l] > 0) iv_add(ib, m.iv_tab[static_cast<std::size_t>(i)][l], m.iv_tab[static_cast<std::size_t>(j)][l]);
                else iv_sub(ib, m.iv_tab[static_cast<std::size_t>(i)][l], m.iv_tab[static_cast<std::size_t>(j)][l]);
                if (mpfr_sgn(ib.lo) <= 0 && mpfr_sgn(ib.hi) >= 0) { degenerate = true; break; }
                const int q = powers[l];
                if (b < 0 && (q % 2 == 1)) t.sign = -t.sign;
                long double ab = fabsl(b);
                for (int e = 0; e < q; ++e) t.coeff *= ab;
                Iv iab = ib;
                if (mpfr_sgn(ib.hi) < 0) { mpfr_neg(iab.lo, ib.hi, MPFR_RNDD); mpfr_neg(iab.hi, ib.lo, MPFR_RNDU); }
                for (int e = 0; e < q; ++e) { Iv tmp; iv_mul(tmp, acc, iab); acc = tmp; }
                t.lam[l] = b * b;
                Iv sq; iv_mul(sq, ib, ib);
                t.ilam[l] = sq;
            }
            if (!degenerate) { t.icoeff = acc; out.push_back(t); }
        }
    }
    return out;
}

// ---------------------------------------------------------------- simplex

// Two-phase dense simplex, maximising c.x subject to A x <= b, x >= 0.
// Dantzig's rule with a Bland fallback and an iteration cap.  Returns false if
// infeasible, unbounded, or abandoned.  The caller re-validates the result.
bool simplex(const std::vector<std::vector<double>>& A, const std::vector<double>& b,
             const std::vector<double>& c, int n_real, std::vector<double>& out) {
    const int m = static_cast<int>(A.size());
    if (m == 0) return false;
    const int n = static_cast<int>(A[0].size());
    const int W = n + m + m + 1;
    std::vector<std::vector<double>> T(static_cast<std::size_t>(m + 1), std::vector<double>(static_cast<std::size_t>(W), 0.0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) T[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        T[static_cast<std::size_t>(i)][static_cast<std::size_t>(n + i)] = 1.0;
        T[static_cast<std::size_t>(i)][static_cast<std::size_t>(W - 1)] = b[static_cast<std::size_t>(i)];
        if (T[static_cast<std::size_t>(i)][static_cast<std::size_t>(W - 1)] < 0.0)
            for (int j = 0; j <= n + m; ++j) {
                if (j < n + m) T[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = -T[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            }
        if (b[static_cast<std::size_t>(i)] < 0.0) T[static_cast<std::size_t>(i)][static_cast<std::size_t>(W - 1)] = -b[static_cast<std::size_t>(i)];
        T[static_cast<std::size_t>(i)][static_cast<std::size_t>(n + m + i)] = 1.0;
    }
    std::vector<int> basis(static_cast<std::size_t>(m));
    for (int i = 0; i < m; ++i) basis[static_cast<std::size_t>(i)] = n + m + i;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < W; ++j)
            T[static_cast<std::size_t>(m)][static_cast<std::size_t>(j)] -= 0.0;
    // phase I objective: minimise sum of artificials
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < W; ++j)
            T[static_cast<std::size_t>(m)][static_cast<std::size_t>(j)] -= T[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    for (int i = 0; i < m; ++i) T[static_cast<std::size_t>(m)][static_cast<std::size_t>(n + m + i)] += 1.0;

    const double tol = 1e-9;
    auto pivot = [&](int limit) -> bool {
        const long cap = 200L * (W + m);
        const long bland_after = 4L * (W + m);
        for (long it = 1;; ++it) {
            if (it > cap) return false;
            int col = -1;
            if (it < bland_after) {
                double best = -tol;
                for (int j = 0; j < limit; ++j)
                    if (T[static_cast<std::size_t>(m)][static_cast<std::size_t>(j)] < best) { best = T[static_cast<std::size_t>(m)][static_cast<std::size_t>(j)]; col = j; }
            } else {
                for (int j = 0; j < limit; ++j)
                    if (T[static_cast<std::size_t>(m)][static_cast<std::size_t>(j)] < -tol) { col = j; break; }
            }
            if (col < 0) return true;
            int row = -1; double best = 0.0;
            for (int i = 0; i < m; ++i) {
                const double a = T[static_cast<std::size_t>(i)][static_cast<std::size_t>(col)];
                if (a > tol) {
                    const double r = T[static_cast<std::size_t>(i)][static_cast<std::size_t>(W - 1)] / a;
                    if (row < 0 || r < best - tol || (std::fabs(r - best) <= tol && basis[static_cast<std::size_t>(i)] < basis[static_cast<std::size_t>(row)])) { best = r; row = i; }
                }
            }
            if (row < 0) return false;
            const double p = T[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)];
            for (int j = 0; j < W; ++j) T[static_cast<std::size_t>(row)][static_cast<std::size_t>(j)] /= p;
            for (int i = 0; i <= m; ++i) {
                if (i == row) continue;
                const double f = T[static_cast<std::size_t>(i)][static_cast<std::size_t>(col)];
                if (f != 0.0)
                    for (int j = 0; j < W; ++j) T[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] -= f * T[static_cast<std::size_t>(row)][static_cast<std::size_t>(j)];
            }
            basis[static_cast<std::size_t>(row)] = col;
        }
    };

    if (!pivot(W - 1)) return false;
    if (T[static_cast<std::size_t>(m)][static_cast<std::size_t>(W - 1)] > tol) return false;   // infeasible

    // drive artificials out of the basis
    for (int i = 0; i < m; ++i) {
        if (basis[static_cast<std::size_t>(i)] >= n + m) {
            int piv = -1;
            for (int j = 0; j < n + m; ++j)
                if (std::fabs(T[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]) > tol) { piv = j; break; }
            if (piv < 0) continue;
            const double p = T[static_cast<std::size_t>(i)][static_cast<std::size_t>(piv)];
            for (int j = 0; j < W; ++j) T[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] /= p;
            for (int r = 0; r <= m; ++r) {
                if (r == i) continue;
                const double f = T[static_cast<std::size_t>(r)][static_cast<std::size_t>(piv)];
                if (f != 0.0)
                    for (int j = 0; j < W; ++j) T[static_cast<std::size_t>(r)][static_cast<std::size_t>(j)] -= f * T[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            }
            basis[static_cast<std::size_t>(i)] = piv;
        }
    }

    // phase II objective
    for (int j = 0; j < W; ++j) T[static_cast<std::size_t>(m)][static_cast<std::size_t>(j)] = 0.0;
    for (int j = 0; j < n; ++j) T[static_cast<std::size_t>(m)][static_cast<std::size_t>(j)] = -c[static_cast<std::size_t>(j)];
    for (int i = 0; i < m; ++i) {
        const int bi = basis[static_cast<std::size_t>(i)];
        if (bi < n + m && T[static_cast<std::size_t>(m)][static_cast<std::size_t>(bi)] != 0.0) {
            const double f = T[static_cast<std::size_t>(m)][static_cast<std::size_t>(bi)];
            for (int j = 0; j < W; ++j) T[static_cast<std::size_t>(m)][static_cast<std::size_t>(j)] -= f * T[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }
    }
    for (int i = 0; i < m; ++i)
        for (int j = n + m; j < n + m + m; ++j) T[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = 0.0;
    if (!pivot(n + m)) return false;

    out.assign(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < m; ++i) {
        const int bi = basis[static_cast<std::size_t>(i)];
        if (bi < n_real) out[static_cast<std::size_t>(bi)] = T[static_cast<std::size_t>(i)][static_cast<std::size_t>(W - 1)];
    }
    return true;
}

// Verify one rational allocation in interval arithmetic.  `lam` and `coef` are
// the regime's interval data, already carrying the minimum-exponent shift.
bool verify_alloc(const std::vector<std::vector<Iv>>& lam, const std::vector<Iv>& coef,
                  const std::vector<int>& pos, const std::vector<int>& neg,
                  const std::vector<std::vector<long>>& units, long denom, std::size_t K) {
    // structural conditions, exact over the integers
    for (std::size_t xi = 0; xi < neg.size(); ++xi) {
        long s = 0;
        for (std::size_t pi = 0; pi < pos.size(); ++pi) s += units[pi][xi];
        if (s != denom) return false;
    }
    for (std::size_t pi = 0; pi < pos.size(); ++pi) {
        long s = 0;
        for (std::size_t xi = 0; xi < neg.size(); ++xi) s += units[pi][xi];
        if (s > denom) return false;
    }
    // logarithms of the interval data
    std::vector<std::vector<Iv>> llam(lam.size());
    for (std::size_t t = 0; t < lam.size(); ++t) {
        llam[t].resize(K);
        for (std::size_t k = 0; k < K; ++k)
            if (!iv_log(llam[t][k], lam[t][k])) return false;
    }
    std::vector<Iv> lcoef(coef.size());
    for (std::size_t t = 0; t < coef.size(); ++t)
        if (!iv_log(lcoef[t], coef[t])) return false;

    for (std::size_t xi = 0; xi < neg.size(); ++xi) {
        for (std::size_t k = 0; k < K; ++k) {
            Iv acc;
            for (std::size_t pi = 0; pi < pos.size(); ++pi)
                if (units[pi][xi])
                    iv_addmul_q(acc, llam[static_cast<std::size_t>(pos[pi])][k], units[pi][xi], denom);
            Iv slack; iv_sub(slack, acc, llam[static_cast<std::size_t>(neg[xi])][k]);
            if (mpfr_sgn(slack.lo) <= 0) return false;
        }
        Iv acc;
        for (std::size_t pi = 0; pi < pos.size(); ++pi)
            if (units[pi][xi])
                iv_addmul_q(acc, lcoef[static_cast<std::size_t>(pos[pi])], units[pi][xi], denom);
        Iv slack; iv_sub(slack, acc, lcoef[static_cast<std::size_t>(neg[xi])]);
        if (mpfr_sgn(slack.lo) <= 0) return false;
    }
    return true;
}


// ------------------------------------------------------- nodes and recursion

// A node of the decomposition: terms over K free coordinates, in both
// arithmetics.  The float side proposes, the interval side decides.
struct Node {
    std::vector<int> sign;
    std::vector<long double> coeff;
    std::vector<std::vector<long double>> lam;   // [term][k]
    std::vector<Iv> icoeff;
    std::vector<std::vector<Iv>> ilam;
    std::size_t k = 0;
};

struct Stats {
    std::uint64_t amgm = 0, leaves = 0, nodes = 0, cap = 400000;
    bool root_lp_ok = false;   // did a root allocation exist at all
    double root_margin = 0.0;  // and with what worst-case slack
    std::array<std::uint64_t, 5> den{};
};

// A zero-dimensional node is one exponent point: decide it outright.
bool leaf_nonneg(const Node& nd) {
    Iv total;
    for (std::size_t t = 0; t < nd.sign.size(); ++t) {
        Iv nxt;
        if (nd.sign[t] > 0) iv_add(nxt, total, nd.icoeff[t]);
        else iv_sub(nxt, total, nd.icoeff[t]);
        total = nxt;
    }
    return mpfr_sgn(total.lo) >= 0;
}

// Try one capacitated allocation at this node.
bool amgm_node(const Node& nd, Stats& st, bool* proposal_existed = nullptr,
               double* margin_out = nullptr) {
    std::vector<int> pos, neg;
    for (std::size_t t = 0; t < nd.sign.size(); ++t)
        (nd.sign[t] < 0 ? neg : pos).push_back(static_cast<int>(t));
    if (neg.empty()) return true;                    // nothing to dominate
    if (pos.empty()) return false;

    const std::size_t P = pos.size(), X = neg.size(), K = nd.k;
    const std::size_t nv = P * X + 1, dj = P * X;
    auto aidx = [&](std::size_t pi, std::size_t xi) { return pi * X + xi; };

    std::vector<std::vector<double>> A;
    std::vector<double> bb;
    for (std::size_t xi = 0; xi < X; ++xi) {
        std::vector<double> r(nv, 0.0);
        for (std::size_t pi = 0; pi < P; ++pi) r[aidx(pi, xi)] = 1.0;
        A.push_back(r); bb.push_back(1.0);
        for (double& q : r) q = -q;
        A.push_back(r); bb.push_back(-1.0);
    }
    for (std::size_t pi = 0; pi < P; ++pi) {
        std::vector<double> r(nv, 0.0);
        for (std::size_t xi = 0; xi < X; ++xi) r[aidx(pi, xi)] = 1.0;
        A.push_back(r); bb.push_back(1.0);
    }
    for (std::size_t xi = 0; xi < X; ++xi) {
        for (std::size_t kk = 0; kk < K; ++kk) {
            std::vector<double> r(nv, 0.0);
            for (std::size_t pi = 0; pi < P; ++pi)
                r[aidx(pi, xi)] = -static_cast<double>(logl(nd.lam[static_cast<std::size_t>(pos[pi])][kk]));
            r[dj] = 1.0;
            A.push_back(r);
            bb.push_back(-static_cast<double>(logl(nd.lam[static_cast<std::size_t>(neg[xi])][kk])));
        }
        std::vector<double> r(nv, 0.0);
        for (std::size_t pi = 0; pi < P; ++pi)
            r[aidx(pi, xi)] = -static_cast<double>(logl(nd.coeff[static_cast<std::size_t>(pos[pi])]));
        r[dj] = 1.0;
        A.push_back(r);
        bb.push_back(-static_cast<double>(logl(nd.coeff[static_cast<std::size_t>(neg[xi])])));
    }
    std::vector<double> c(nv, 0.0); c[dj] = 1.0;
    std::vector<double> sol;
    if (!simplex(A, bb, c, static_cast<int>(nv), sol)) return false;

    // Never trust the solver: discard a proposal violating its own constraints.
    for (double q : sol) if (q < -1e-7) return false;
    for (std::size_t i = 0; i < A.size(); ++i) {
        double acc = 0.0;
        for (std::size_t j = 0; j < nv; ++j) acc += A[i][j] * sol[j];
        if (acc - bb[i] > 1e-6) return false;
    }
    if (sol[dj] <= 0.0) return false;
    if (proposal_existed) *proposal_existed = true;   // an allocation exists
    if (margin_out) *margin_out = sol[dj];            // its worst-case slack

    // Rounding, done properly.
    //
    // Round-to-nearest with a crude fixup was losing certificates that the
    // linear program had already found: at rank six 158 of 326 open regimes,
    // and at rank seven 3,631 of 6,344, had a feasible root allocation whose
    // certificate died in the conversion to a rational.  A finer denominator
    // does not rescue them, which points at the rounding rule rather than the
    // resolution.
    //
    // So start from the floor and hand out the remaining units one at a time,
    // each to whichever positive term maximises the resulting minimum slack.
    // Capacity is tracked across columns, since spending a positive term on one
    // negative leaves less of it for the next.
    std::vector<double> G;                    // per positive, per constraint
    const std::size_t NC = K + 1;             // K domination rows plus the constant row
    G.assign(P * NC, 0.0);
    for (std::size_t pi = 0; pi < P; ++pi) {
        for (std::size_t kk = 0; kk < K; ++kk)
            G[pi * NC + kk] = static_cast<double>(logl(nd.lam[static_cast<std::size_t>(pos[pi])][kk]));
        G[pi * NC + K] = static_cast<double>(logl(nd.coeff[static_cast<std::size_t>(pos[pi])]));
    }
    std::vector<double> H(X * NC, 0.0);
    for (std::size_t xi = 0; xi < X; ++xi) {
        for (std::size_t kk = 0; kk < K; ++kk)
            H[xi * NC + kk] = static_cast<double>(logl(nd.lam[static_cast<std::size_t>(neg[xi])][kk]));
        H[xi * NC + K] = static_cast<double>(logl(nd.coeff[static_cast<std::size_t>(neg[xi])]));
    }

    static const long ladder[5] = {100, 200, 500, 1000, 5000};
    for (int di = 0; di < 5; ++di) {
        const long d = ladder[di];
        std::vector<long> cap_left(P, d);
        std::vector<std::vector<long>> units(P, std::vector<long>(X, 0));
        bool ok = true;
        for (std::size_t xi = 0; xi < X && ok; ++xi) {
            long assigned = 0;
            for (std::size_t pi = 0; pi < P; ++pi) {
                long u = static_cast<long>(std::floor(sol[aidx(pi, xi)] * static_cast<double>(d)));
                if (u < 0) u = 0;
                if (u > cap_left[pi]) u = cap_left[pi];
                units[pi][xi] = u;
                assigned += u;
            }
            // Hand out the remainder greedily by resulting minimum slack.
            std::vector<double> acc(NC, 0.0);
            for (std::size_t pi = 0; pi < P; ++pi)
                if (units[pi][xi])
                    for (std::size_t cc = 0; cc < NC; ++cc)
                        acc[cc] += static_cast<double>(units[pi][xi]) * G[pi * NC + cc];
            while (assigned < d) {
                int best = -1;
                double best_score = 0.0;
                for (std::size_t pi = 0; pi < P; ++pi) {
                    if (units[pi][xi] >= cap_left[pi]) continue;
                    double worst = 0.0;
                    bool first = true;
                    for (std::size_t cc = 0; cc < NC; ++cc) {
                        const double v = (acc[cc] + G[pi * NC + cc]) / static_cast<double>(d) - H[xi * NC + cc];
                        if (first || v < worst) { worst = v; first = false; }
                    }
                    if (best < 0 || worst > best_score) { best_score = worst; best = static_cast<int>(pi); }
                }
                if (best < 0) { ok = false; break; }      // capacity exhausted
                const std::size_t bp = static_cast<std::size_t>(best);
                units[bp][xi] += 1;
                for (std::size_t cc = 0; cc < NC; ++cc) acc[cc] += G[bp * NC + cc];
                ++assigned;
            }
            if (!ok) break;
            for (std::size_t pi = 0; pi < P; ++pi) cap_left[pi] -= units[pi][xi];
        }
        if (!ok) continue;
        if (verify_alloc(nd.ilam, nd.icoeff, pos, neg, units, d, K)) {
            ++st.den[static_cast<std::size_t>(di)];
            return true;
        }
    }
    return false;
}

// How deep to split, by dimension.  A uniform cap is wrong: the recorded O11
// final chamber needed a threshold of 75 before its tail closed, with 73,150
// finite lattice checks below it, so a ray can need tens of splits.  Depth is
// cheap in low dimension, where one split costs one leaf, and expensive in high
// dimension, where the face subtree branches.  The node budget is the real
// ceiling.
int depth_for(std::size_t k) {
    switch (k) {
        case 0: return 0;
        case 1: return 240;
        case 2: return 40;
        case 3: return 14;
        case 4: return 9;
        default: return 6;
    }
}

// Close a node by allocation, or split it and recurse.
//
// For a free coordinate l the index set splits exactly:
//     {u_l >= 0} = {u_l = 0} disjoint-union {u_l >= 1}.
// The face drops l; the tail absorbs one factor of lambda into every
// coefficient.  Both are strictly simpler and their union is the parent, so the
// split is sound by construction.  Acceptance never widens: a node closes only
// by an interval-verified allocation or an interval evaluation.
bool certify_node(const Node& nd, int depth, int max_depth, Stats& st) {
    if (++st.nodes > st.cap) return false;
    if (nd.k == 0) { ++st.leaves; return leaf_nonneg(nd); }
    if (amgm_node(nd, st, depth == 0 ? &st.root_lp_ok : nullptr,
                  depth == 0 ? &st.root_margin : nullptr)) { ++st.amgm; return true; }
    if (depth >= max_depth) return false;

    // Split where the negative terms hold the largest advantage: that is the
    // obstruction the allocation could not pay for.
    std::size_t best = 0;
    long double best_score = -1.0L;
    for (std::size_t kk = 0; kk < nd.k; ++kk) {
        long double worst = 0.0L;
        for (std::size_t t = 0; t < nd.sign.size(); ++t)
            if (nd.sign[t] < 0) worst = std::max(worst, nd.lam[t][kk]);
        if (worst > best_score) { best_score = worst; best = kk; }
    }

    // Face u_best = 0: drop the coordinate.
    Node face;
    face.sign = nd.sign;
    face.coeff = nd.coeff;
    face.icoeff = nd.icoeff;
    face.k = nd.k - 1;
    face.lam.resize(nd.sign.size());
    face.ilam.resize(nd.sign.size());
    for (std::size_t t = 0; t < nd.sign.size(); ++t)
        for (std::size_t kk = 0; kk < nd.k; ++kk)
            if (kk != best) { face.lam[t].push_back(nd.lam[t][kk]); face.ilam[t].push_back(nd.ilam[t][kk]); }
    if (!certify_node(face, 0, depth_for(face.k), st)) return false;

    // Tail u_best >= 1: absorb one factor of lambda.
    Node tail;
    tail.sign = nd.sign;
    tail.lam = nd.lam;
    tail.ilam = nd.ilam;
    tail.k = nd.k;
    tail.coeff.resize(nd.sign.size());
    tail.icoeff.resize(nd.sign.size());
    for (std::size_t t = 0; t < nd.sign.size(); ++t) {
        tail.coeff[t] = nd.coeff[t] * nd.lam[t][best];
        iv_mul(tail.icoeff[t], nd.icoeff[t], nd.ilam[t][best]);
    }
    return certify_node(tail, depth + 1, max_depth, st);
}

}  // namespace

// ---------------------------------------------------------------- driver

int main(int argc, char** argv) {
    int rank = 6, threads = static_cast<int>(std::thread::hardware_concurrency());
    bool decompose = false;
    int soundness = 0;      // sample this many certified regimes
    bool control = false;   // check the verifier rejects corrupted allocations
    bool dump = false;      // print one line per unresolved regime
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--rank") && i + 1 < argc) rank = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--threads") && i + 1 < argc) threads = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--decompose")) decompose = true;
        else if (!std::strcmp(argv[i], "--soundness") && i + 1 < argc) soundness = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--control")) control = true;
        else if (!std::strcmp(argv[i], "--dump-open")) dump = true;
    }
    if (rank < 3) { std::fprintf(stderr, "rank must be at least 3\n"); return 2; }
    if (threads < 1) threads = 1;

    Model model;
    build_model(model, rank);
    const int L = model.labels;

    // Cross-check the two spectral models against each other.  They are built
    // independently -- long double trigonometry versus 512-bit interval
    // trigonometry -- so agreement is real evidence that neither is wrong.
    //
    // The comparison must be at long double precision, not interval precision.
    // The enclosures here are tight to about 1e-154 while the long double value
    // carries its own error of about 1e-19, so the float value legitimately
    // sits outside the interval; demanding containment would always fail.
    {
        int checked = 0;
        const long double tol = 1.0e-15L;
        auto agrees = [&](const Iv& iv, long double f) {
            const long double lo = mpfr_get_ld(iv.lo, MPFR_RNDD);
            const long double hi = mpfr_get_ld(iv.hi, MPFR_RNDU);
            const long double mid = 0.5L * (lo + hi);
            return fabsl(mid - f) <= tol * (1.0L + fabsl(f));
        };
        for (int i = 0; i < rank; ++i) {
            const std::size_t si = static_cast<std::size_t>(i);
            if (!agrees(model.iw[si], model.w[si])) {
                std::fprintf(stderr, "spectral cross-check failed: weight %d\n", i);
                return 3;
            }
            ++checked;
            for (int l = 0; l < L; ++l) {
                const std::size_t sl = static_cast<std::size_t>(l);
                if (!agrees(model.iv_tab[si][sl], model.v[si][sl])) {
                    std::fprintf(stderr, "spectral cross-check failed: chi node %d label %d\n", i, l);
                    return 3;
                }
                ++checked;
            }
        }
        std::printf("  spectral_crosscheck %d values agree\n", checked);
    }

    std::atomic<int> next{0};
    std::mutex mu;
    std::mutex dump_mu;
    std::uint64_t direct = 0, residual = 0, pointwise = 0, sep = 0;
    std::uint64_t certified_total = 0;
    std::array<std::uint64_t, 5> den_used{};
    std::uint64_t sound_regimes = 0, sound_points = 0, sound_viol = 0;
    std::uint64_t ctl_tried = 0, ctl_refused = 0;
    std::uint64_t flat_total = 0, split_total = 0;
    std::uint64_t node_total = 0, budget_hit = 0;

    const int supports = static_cast<int>(std::llround(std::pow(3.0, L)));

    auto worker = [&]() {
        std::uint64_t ld = 0, lr = 0, lp = 0, ls = 0, lcert = 0;
        std::array<std::uint64_t, 5> lden{};
        std::uint64_t lsound_regimes = 0, lsound_points = 0, lsound_viol = 0;
        std::uint64_t lflat = 0, lsplit = 0, lnodes = 0, lbudget = 0;
        std::uint64_t lctl_tried = 0, lctl_refused = 0;
        for (;;) {
            const int support = next.fetch_add(1);
            if (support >= supports) break;
            std::vector<int> signs(static_cast<std::size_t>(L), 0);
            int z = support, active = 0, minus = 0;
            for (int l = 0; l < L; ++l) {
                const int d = z % 3; z /= 3;
                signs[static_cast<std::size_t>(l)] = (d == 0) ? 0 : (d == 1 ? 1 : -1);
                if (d) { ++active; if (d == 2) ++minus; }
            }
            if (!active || !minus) continue;
            std::vector<int> act;
            for (int l = 0; l < L; ++l) if (signs[static_cast<std::size_t>(l)]) act.push_back(l);

            for (int parity = 0; parity < (1 << active); ++parity) {
                std::vector<int> powers(static_cast<std::size_t>(L), 0);
                int mp = 0;
                for (std::size_t i = 0; i < act.size(); ++i) {
                    const int l = act[i];
                    const bool odd = ((parity >> i) & 1) != 0;
                    powers[static_cast<std::size_t>(l)] = odd ? 1 : 2;
                    if (signs[static_cast<std::size_t>(l)] < 0 && odd) mp ^= 1;
                }
                if (mp) continue;
                std::vector<Term> terms = chamber_terms(model, signs, powers);
                if (terms.empty()) continue;
                bool anyneg = false;
                for (const Term& t : terms) if (t.sign < 0) anyneg = true;
                if (!anyneg) { lp += static_cast<std::uint64_t>(1) << act.size(); continue; }

                for (int res = 0; res < (1 << active); ++res) {
                    std::vector<int> freev;
                    for (std::size_t i = 0; i < act.size(); ++i)
                        if ((res >> i) & 1) freev.push_back(act[i]);
                    // Hall transport in long double, matching the census.
                    std::vector<long double> mag(terms.size());
                    for (std::size_t t = 0; t < terms.size(); ++t) {
                        long double v = terms[t].coeff;
                        for (int l : freev) v *= terms[t].lam[static_cast<std::size_t>(l)];
                        mag[t] = v;
                    }
                    std::vector<int> neg, pos;
                    for (std::size_t t = 0; t < terms.size(); ++t)
                        (terms[t].sign < 0 ? neg : pos).push_back(static_cast<int>(t));
                    auto ge = [&](long double a, long double b) {
                        return a + kHallTol * (1.0L + fabsl(a) + fabsl(b)) >= b;
                    };
                    std::vector<std::uint32_t> edge(neg.size(), 0);
                    for (std::size_t x = 0; x < neg.size(); ++x)
                        for (std::size_t p = 0; p < pos.size(); ++p) {
                            bool ok = true;
                            for (int l : freev)
                                if (!ge(terms[static_cast<std::size_t>(pos[p])].lam[static_cast<std::size_t>(l)],
                                        terms[static_cast<std::size_t>(neg[x])].lam[static_cast<std::size_t>(l)])) { ok = false; break; }
                            if (ok) edge[x] |= (1U << p);
                        }
                    bool hall = true;
                    const std::uint64_t lim = std::uint64_t{1} << neg.size();
                    for (std::uint64_t s = 1; s < lim && hall; ++s) {
                        long double d = 0.0L, cap = 0.0L;
                        std::uint32_t nb = 0;
                        for (std::size_t i = 0; i < neg.size(); ++i)
                            if ((s >> i) & 1U) { d += mag[static_cast<std::size_t>(neg[i])]; nb |= edge[i]; }
                        for (std::size_t p = 0; p < pos.size(); ++p)
                            if ((nb >> p) & 1U) cap += mag[static_cast<std::size_t>(pos[p])];
                        if (!ge(cap, d)) hall = false;
                    }
                    if (hall) { ++ld; continue; }
                    ++lr;

                    // Build the root node and certify it.  With --decompose
                    // off the depth is zero, so this is exactly the flat
                    // allocation stage and must reproduce its numbers.
                    const std::size_t K = freev.size();
                    Node root;
                    root.k = K;
                    root.sign.resize(terms.size());
                    root.coeff.resize(terms.size());
                    root.icoeff.resize(terms.size());
                    root.lam.assign(terms.size(), {});
                    root.ilam.assign(terms.size(), {});
                    for (std::size_t t = 0; t < terms.size(); ++t) {
                        root.sign[t] = terms[t].sign;
                        long double cc = terms[t].coeff;
                        Iv icc = terms[t].icoeff;
                        for (std::size_t kk = 0; kk < K; ++kk) {
                            const std::size_t l = static_cast<std::size_t>(freev[kk]);
                            root.lam[t].push_back(terms[t].lam[l]);
                            root.ilam[t].push_back(terms[t].ilam[l]);
                            cc *= terms[t].lam[l];
                            Iv tmp; iv_mul(tmp, icc, terms[t].ilam[l]); icc = tmp;
                        }
                        root.coeff[t] = cc;
                        root.icoeff[t] = icc;
                    }
                    Stats st;
                    const int md = decompose ? depth_for(K) : 0;
                    const bool certified = certify_node(root, 0, md, st);
                    if (!certified) {
                        // Distinguish "no certificate exists at this depth"
                        // from "ran out of node budget", since only the second
                        // is fixable by spending more.
                        if (st.nodes > st.cap) ++lbudget;
                        if (dump) {
                            std::size_t nneg = 0;
                            for (int sg : root.sign) if (sg < 0) ++nneg;
                            std::lock_guard<std::mutex> dg(dump_mu);
                            std::printf("OPEN support=%d parity=%d k=%zu neg=%zu pos=%zu"
                                        " root_lp=%s margin=%.3e nodes=%llu\n",
                                        support, parity, K, nneg, root.sign.size() - nneg,
                                        st.root_lp_ok ? "feasible" : "infeasible", st.root_margin,
                                        static_cast<unsigned long long>(st.nodes));
                        }
                        continue;
                    }
                    ++lcert;
                    if (st.amgm == 1 && st.leaves == 0) ++lflat; else ++lsplit;
                    lnodes += st.nodes;
                    for (std::size_t di = 0; di < 5; ++di) lden[di] += st.den[di];

                    // Soundness: the certificate claims the signed sum is
                    // nonnegative across the whole regime.  Evaluate the actual
                    // sum at pseudo-random exponent points, independently of the
                    // certificate logic.
                    if (soundness && lsound_regimes < static_cast<std::uint64_t>(soundness)) {
                        ++lsound_regimes;
                        std::uint64_t rng = static_cast<std::uint64_t>(0x9e3779b97f4a7c15ULL);
                        rng ^= static_cast<std::uint64_t>(support) << 20;
                        rng ^= static_cast<std::uint64_t>(parity) << 8;
                        rng ^= static_cast<std::uint64_t>(res);
                        for (int trial = 0; trial < 40; ++trial) {
                            std::vector<int> expo(K);
                            for (std::size_t kk = 0; kk < K; ++kk) {
                                rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
                                expo[kk] = 1 + static_cast<int>(rng % 12);
                            }
                            long double tot = 0.0L;
                            for (std::size_t t = 0; t < terms.size(); ++t) {
                                long double v = terms[t].coeff;
                                for (std::size_t kk = 0; kk < K; ++kk)
                                    for (int q = 0; q < expo[kk]; ++q)
                                        v *= terms[t].lam[static_cast<std::size_t>(freev[kk])];
                                tot += static_cast<long double>(terms[t].sign) * v;
                            }
                            ++lsound_points;
                            if (tot < -1.0e-18L * (1.0L + fabsl(tot))) ++lsound_viol;
                        }
                    }

                    // Negative control: corrupt a verified allocation at the
                    // root and require the verifier to reject it.
                    if (control && lctl_tried < 200 && st.amgm >= 1 && K > 0) {
                        std::vector<int> cpos, cneg;
                        for (std::size_t t = 0; t < root.sign.size(); ++t)
                            (root.sign[t] < 0 ? cneg : cpos).push_back(static_cast<int>(t));
                        if (!cneg.empty() && !cpos.empty()) {
                            const long d = 100;
                            std::vector<std::vector<long>> u(cpos.size(), std::vector<long>(cneg.size(), 0));
                            for (std::size_t xi = 0; xi < cneg.size(); ++xi) u[0][xi] = d;
                            ++lctl_tried;                              // capacity blown
                            if (!verify_alloc(root.ilam, root.icoeff, cpos, cneg, u, d, K)) ++lctl_refused;
                            for (std::size_t xi = 0; xi < cneg.size(); ++xi) u[0][xi] = d + 1;
                            ++lctl_tried;                              // normalisation broken
                            if (!verify_alloc(root.ilam, root.icoeff, cpos, cneg, u, d, K)) ++lctl_refused;
                        }
                    }
                }
            }
        }
        std::lock_guard<std::mutex> g(mu);
        direct += ld; residual += lr; pointwise += lp; sep += ls;
        certified_total += lcert;
        for (std::size_t i = 0; i < 5; ++i) den_used[i] += lden[i];
        sound_regimes += lsound_regimes; sound_points += lsound_points; sound_viol += lsound_viol;
        ctl_tried += lctl_tried; ctl_refused += lctl_refused;
        flat_total += lflat; split_total += lsplit;
        node_total += lnodes; budget_hit += lbudget;
    };

    std::vector<std::thread> pool;
    for (int i = 0; i < threads; ++i) pool.emplace_back(worker);
    for (auto& t : pool) t.join();

    std::printf("SU2_ORBIT_AMGM_CPP rank=%d level=%d threads=%d\n", rank, 2 * rank - 1, threads);
    std::printf("  direct_hall      %llu\n", static_cast<unsigned long long>(direct));
    std::printf("  residual         %llu\n", static_cast<unsigned long long>(residual));
    std::printf("  pointwise        %llu\n", static_cast<unsigned long long>(pointwise));
    std::printf("  certified        %llu\n", static_cast<unsigned long long>(certified_total));
    std::printf("  by_flat_amgm     %llu\n", static_cast<unsigned long long>(flat_total));
    std::printf("  by_split         %llu\n", static_cast<unsigned long long>(split_total));
    std::printf("  split_nodes      %llu\n", static_cast<unsigned long long>(node_total));
    std::printf("  budget_exhausted %llu\n", static_cast<unsigned long long>(budget_hit));
    static const int ladder_out[5] = {100, 200, 500, 1000, 5000};
    std::printf("  denominators    ");
    for (std::size_t i = 0; i < 5; ++i)
        if (den_used[i]) std::printf(" %d:%llu", ladder_out[i], static_cast<unsigned long long>(den_used[i]));
    std::printf("\n");
    if (sound_points)
        std::printf("  soundness        regimes=%llu points=%llu violations=%llu\n",
                    static_cast<unsigned long long>(sound_regimes),
                    static_cast<unsigned long long>(sound_points),
                    static_cast<unsigned long long>(sound_viol));
    if (ctl_tried)
        std::printf("  negative_control corrupted=%llu refused=%llu\n",
                    static_cast<unsigned long long>(ctl_tried),
                    static_cast<unsigned long long>(ctl_refused));
    if (residual)
        std::printf("  coverage         %.1f%% of residual regimes\n",
                    100.0 * static_cast<double>(certified_total) / static_cast<double>(residual));
    return 0;
}
