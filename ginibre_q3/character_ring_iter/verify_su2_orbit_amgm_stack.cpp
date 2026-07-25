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

#include "su2_cyclotomic.h"

namespace {

constexpr mpfr_prec_t kPrec = 512;
constexpr long double kHallTol = 1.0e-12L;
bool g_total_mass_sweep = false;
// Full fusion-ring mode: labels 1..k at all k+1 Verlinde nodes, instead of the
// odd-level orbit ring.  The two-node spectral structure is identical, so the
// whole certificate pipeline applies unchanged; the orbit reduction was an
// optimisation, never a prerequisite.  This is what opens the even levels.
bool g_full_level = false;
// In full-level mode every regime is taken through the LP + interval pipeline;
// the long double Hall pre-filter is only a speed statistic, and a closed
// level must not rest on a tolerance-based test.
bool g_certify_all = false;
// Fusion level for exact integer leaf evaluation.  In full-level mode this
// is the level itself with labels 1..k; the orbit ring embeds as the even
// labels 2,4,... at level 2*rank-1.
int g_leaf_level = 0;
int g_label_stride = 1;
// Exact cyclotomic field data for tied-row escalation.  A certificate row the
// interval test cannot decide is settled by the sign of one element of
// Z[2cos(pi/n)]; zero means the row holds with equality, which weighted AM-GM
// permits.  This is what closes the exact-tie residue.
cyclo::Field g_field;
std::vector<cyclo::Field::Elt> g_exact_y;                    // per node
std::vector<std::vector<cyclo::Field::Elt>> g_exact_chi;     // per node, label

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

// Sound enclosure of sin over an interval that may contain an extremum.
// The endpoint enclosure is only valid where sin is monotone; when the cosine
// changes sign across the interval an extremum lies inside, and the affected
// bound is widened to the true extreme value +-1.
void iv_sin_enclose(Iv& r, const Iv& a) {
    mpfr_t c0, c1;
    mpfr_init2(c0, kPrec); mpfr_init2(c1, kPrec);
    mpfr_sin(c0, a.lo, MPFR_RNDD);
    mpfr_sin(c1, a.hi, MPFR_RNDD);
    mpfr_set(r.lo, mpfr_less_p(c0, c1) ? c0 : c1, MPFR_RNDD);
    mpfr_sin(c0, a.lo, MPFR_RNDU);
    mpfr_sin(c1, a.hi, MPFR_RNDU);
    mpfr_set(r.hi, mpfr_greater_p(c0, c1) ? c0 : c1, MPFR_RNDU);
    mpfr_cos(c0, a.lo, MPFR_RNDD);
    mpfr_cos(c1, a.hi, MPFR_RNDU);
    const int s0 = mpfr_sgn(c0), s1 = mpfr_sgn(c1);
    if (s0 >= 0 && s1 <= 0) mpfr_set_si(r.hi, 1, MPFR_RNDU);    // max inside
    if (s0 <= 0 && s1 >= 0) mpfr_set_si(r.lo, -1, MPFR_RNDD);   // min inside
    mpfr_clear(c0); mpfr_clear(c1);
}

// x / d for a divisor interval with d.lo > 0.
void iv_div_pos(Iv& r, const Iv& x, const Iv& d) {
    mpfr_div(r.lo, x.lo, mpfr_sgn(x.lo) >= 0 ? d.hi : d.lo, MPFR_RNDD);
    mpfr_div(r.hi, x.hi, mpfr_sgn(x.hi) >= 0 ? d.lo : d.hi, MPFR_RNDU);
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

// Full fusion-ring spectrum of SU(2)_k, in the normalisation of the level 1-4
// note: node j = 0..k carries weight w_j = 2/(k+2) sin^2((j+1)pi/(k+2)) and
// character values chi_a(j) = sin((a+1)(j+1)pi/(k+2)) / sin((j+1)pi/(k+2)).
// Validated end to end against exact integer fusion DP by
// verify_su2_level_spectral_dp.py before this builder was written.
void build_model_level(Model& m, int k) {
    m.rank = k + 1;            // node count; the field name predates this mode
    m.labels = k;
    m.order = k + 2;
    const long double pi = acosl(-1.0L);
    m.w.assign(static_cast<std::size_t>(m.rank), 0.0L);
    m.v.assign(static_cast<std::size_t>(m.rank),
               std::vector<long double>(static_cast<std::size_t>(m.labels), 0.0L));
    for (int j = 0; j <= k; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        const long double ang = pi * static_cast<long double>(j + 1) / static_cast<long double>(m.order);
        m.w[sj] = 2.0L * sinl(ang) * sinl(ang) / static_cast<long double>(m.order);
        for (int a = 1; a <= k; ++a)
            m.v[sj][static_cast<std::size_t>(a - 1)] =
                sinl(static_cast<long double>(a + 1) * ang) / sinl(ang);
    }
    m.iw.resize(static_cast<std::size_t>(m.rank));
    m.iv_tab.assign(static_cast<std::size_t>(m.rank),
                    std::vector<Iv>(static_cast<std::size_t>(m.labels)));
    Iv ivpi;
    mpfr_const_pi(ivpi.lo, MPFR_RNDD);
    mpfr_const_pi(ivpi.hi, MPFR_RNDU);
    for (int j = 0; j <= k; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        Iv ang;
        mpfr_mul_si(ang.lo, ivpi.lo, j + 1, MPFR_RNDD);
        mpfr_div_si(ang.lo, ang.lo, m.order, MPFR_RNDD);
        mpfr_mul_si(ang.hi, ivpi.hi, j + 1, MPFR_RNDU);
        mpfr_div_si(ang.hi, ang.hi, m.order, MPFR_RNDU);
        // The base angle lies in (0, pi), where sin may pass its maximum, so
        // the extremum-aware enclosure is required here, unlike the orbit case.
        Iv s;
        iv_sin_enclose(s, ang);
        if (mpfr_sgn(s.lo) <= 0) {
            std::fprintf(stderr, "level model: sin enclosure not positive at node %d\n", j);
            std::abort();
        }
        Iv s2; iv_mul(s2, s, s);
        mpfr_mul_ui(m.iw[sj].lo, s2.lo, 2, MPFR_RNDD);
        mpfr_div_si(m.iw[sj].lo, m.iw[sj].lo, m.order, MPFR_RNDD);
        mpfr_mul_ui(m.iw[sj].hi, s2.hi, 2, MPFR_RNDU);
        mpfr_div_si(m.iw[sj].hi, m.iw[sj].hi, m.order, MPFR_RNDU);
        for (int a = 1; a <= k; ++a) {
            Iv na;
            mpfr_mul_si(na.lo, ang.lo, a + 1, MPFR_RNDD);
            mpfr_mul_si(na.hi, ang.hi, a + 1, MPFR_RNDU);
            Iv sn;
            iv_sin_enclose(sn, na);
            iv_div_pos(m.iv_tab[sj][static_cast<std::size_t>(a - 1)], sn, s);
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
    // Exact counterparts: the coefficient numerator (the positive rational
    // scalar common to every term cancels in row comparisons) and the exact
    // lambda per label.
    cyclo::Field::Elt ecoeff;
    std::vector<cyclo::Field::Elt> elam;
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
            // Exact numerator: (4 - y_i^2)(4 - y_j^2) = (2n)^2 w_i w_j up to
            // the positive rational scalar shared by every term.
            const cyclo::Field::Elt four = g_field.fromInt(4);
            cyclo::Field::Elt eacc = g_field.mulE(
                g_field.sub(four, g_field.mulE(g_exact_y[static_cast<std::size_t>(i)], g_exact_y[static_cast<std::size_t>(i)])),
                g_field.sub(four, g_field.mulE(g_exact_y[static_cast<std::size_t>(j)], g_exact_y[static_cast<std::size_t>(j)])));
            t.elam.assign(L, g_field.fromInt(1));
            bool degenerate = false;
            for (std::size_t l = 0; l < L; ++l) {
                if (signs[l] == 0) continue;
                const long double b = m.v[static_cast<std::size_t>(i)][l] +
                    static_cast<long double>(signs[l]) * m.v[static_cast<std::size_t>(j)][l];
                Iv ib;
                if (signs[l] > 0) iv_add(ib, m.iv_tab[static_cast<std::size_t>(i)][l], m.iv_tab[static_cast<std::size_t>(j)][l]);
                else iv_sub(ib, m.iv_tab[static_cast<std::size_t>(i)][l], m.iv_tab[static_cast<std::size_t>(j)][l]);
                const bool iv_zero = (mpfr_sgn(ib.lo) <= 0 && mpfr_sgn(ib.hi) >= 0);
                // Full-level spectra have genuine interior character zeros:
                // chi_a(j) = 0 exactly whenever (a+1)(j+1) is a multiple of
                // k+2.  The float side computes sin(m*pi) as roundoff of order
                // 1e-19 rather than 0, so its zero test must sit at roundoff
                // scale.  A genuine nonzero base for these small-height
                // algebraic values is at least of order 1e-2, five orders
                // above the threshold, so the two cannot be confused.
                const bool fl_zero = (fabsl(b) < 1.0e-12L);
                // A base is dropped only when both arithmetics call it zero;
                // any disagreement means a model bug and must stop the run.
                if (iv_zero != fl_zero) {
                    std::fprintf(stderr, "degenerate-base disagreement at pair (%d,%d) label %zu\n", i, j, l);
                    std::abort();
                }
                if (iv_zero) { degenerate = true; break; }
                // The interval is authoritative for the sign; the float must
                // agree or the run stops.
                const bool iv_negative = (mpfr_sgn(ib.hi) < 0);
                if (iv_negative != (b < 0)) {
                    std::fprintf(stderr, "sign disagreement at pair (%d,%d) label %zu\n", i, j, l);
                    std::abort();
                }
                cyclo::Field::Elt eb;
                if (signs[l] > 0)
                    eb = g_field.add(g_exact_chi[static_cast<std::size_t>(i)][l], g_exact_chi[static_cast<std::size_t>(j)][l]);
                else
                    eb = g_field.sub(g_exact_chi[static_cast<std::size_t>(i)][l], g_exact_chi[static_cast<std::size_t>(j)][l]);
                const int q = powers[l];
                if (iv_negative && (q % 2 == 1)) t.sign = -t.sign;
                cyclo::Field::Elt eabs = iv_negative ? g_field.sub(g_field.fromInt(0), eb) : eb;
                for (int e = 0; e < q; ++e) eacc = g_field.mulE(eacc, eabs);
                t.elam[l] = g_field.mulE(eb, eb);
                long double ab = fabsl(b);
                for (int e = 0; e < q; ++e) t.coeff *= ab;
                Iv iab = ib;
                if (mpfr_sgn(ib.hi) < 0) { mpfr_neg(iab.lo, ib.hi, MPFR_RNDD); mpfr_neg(iab.hi, ib.lo, MPFR_RNDU); }
                for (int e = 0; e < q; ++e) { Iv tmp; iv_mul(tmp, acc, iab); acc = tmp; }
                t.lam[l] = b * b;
                Iv sq; iv_mul(sq, ib, ib);
                t.ilam[l] = sq;
            }
            if (!degenerate) { t.icoeff = acc; t.ecoeff = eacc; out.push_back(t); }
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
// Total allocated mass S = sn/sd >= 1.  Weighted AM-GM needs only the weights
// inside the geometric mean to sum to one, so with beta = alpha/S the
// domination rows are unchanged and the constant row gains a free log S.
// S = 1 is the ordinary convex certificate.
// Exact escalation of one row: prod_p base_p^{u_p} >= base_x^denom, decided in
// Z[zeta].  Returns true when the exact sign is nonnegative (equality allowed).
bool exact_row_ok(const std::vector<cyclo::Field::Elt>& base,
                  const std::vector<int>& pos, int x,
                  const std::vector<std::vector<long>>& units, std::size_t xi,
                  long denom) {
    cyclo::Field::Elt lhs = g_field.fromInt(1);
    for (std::size_t pi = 0; pi < pos.size(); ++pi) {
        const long u = units[pi][xi];
        if (u > 0)
            lhs = g_field.mulE(lhs, g_field.powE(base[static_cast<std::size_t>(pos[pi])],
                                                 static_cast<unsigned long>(u)));
    }
    cyclo::Field::Elt rhs = g_field.powE(base[static_cast<std::size_t>(x)],
                                         static_cast<unsigned long>(denom));
    return g_field.sign(g_field.sub(lhs, rhs)) >= 0;
}

bool verify_alloc(const std::vector<std::vector<Iv>>& lam, const std::vector<Iv>& coef,
                  const std::vector<int>& pos, const std::vector<int>& neg,
                  const std::vector<std::vector<long>>& units, long denom, std::size_t K,
                  long sn = 1, long sd = 1,
                  const std::vector<std::vector<cyclo::Field::Elt>>* elam = nullptr,
                  const std::vector<cyclo::Field::Elt>* ecoef = nullptr) {
    // structural conditions, exact over the integers
    for (std::size_t xi = 0; xi < neg.size(); ++xi) {
        long s = 0;
        for (std::size_t pi = 0; pi < pos.size(); ++pi) s += units[pi][xi];
        if (s != denom) return false;
    }
    // capacity: sum_x S * beta[p][x] <= 1, i.e. sn * sum_x u <= sd * denom
    for (std::size_t pi = 0; pi < pos.size(); ++pi) {
        long s = 0;
        for (std::size_t xi = 0; xi < neg.size(); ++xi) s += units[pi][xi];
        if (sn * s > sd * denom) return false;
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
            if (mpfr_sgn(slack.lo) <= 0) {
                // The interval cannot certify a row that holds with exact
                // equality.  Escalate to the exact field: the row is the sign
                // of one element of Z[zeta], and zero passes.
                if (elam == nullptr) return false;
                std::vector<cyclo::Field::Elt> base(elam->size());
                for (std::size_t t = 0; t < elam->size(); ++t) base[t] = (*elam)[t][k];
                if (!exact_row_ok(base, pos, neg[xi], units, xi, denom)) return false;
            }
        }
        Iv acc;
        for (std::size_t pi = 0; pi < pos.size(); ++pi)
            if (units[pi][xi])
                iv_addmul_q(acc, lcoef[static_cast<std::size_t>(pos[pi])], units[pi][xi], denom);
        if (sn != sd) {                       // add log S, enclosed rigorously
            Iv sv;
            mpfr_set_si(sv.lo, sn, MPFR_RNDD); mpfr_div_si(sv.lo, sv.lo, sd, MPFR_RNDD);
            mpfr_set_si(sv.hi, sn, MPFR_RNDU); mpfr_div_si(sv.hi, sv.hi, sd, MPFR_RNDU);
            Iv ls;
            if (!iv_log(ls, sv)) return false;
            Iv t; iv_add(t, acc, ls); acc = t;
        }
        Iv slack; iv_sub(slack, acc, lcoef[static_cast<std::size_t>(neg[xi])]);
        if (mpfr_sgn(slack.lo) <= 0) {
            // Exact escalation of the constant row; the shared positive
            // rational scalar cancels because the weights sum to denom on
            // both sides.  Only valid without the total-mass shift.
            if (ecoef == nullptr || sn != sd) return false;
            if (!exact_row_ok(*ecoef, pos, neg[xi], units, xi, denom)) return false;
        }
    }
    return true;
}


// ------------------------------------------------------- nodes and recursion
struct Node;
struct Stats;
bool amgm_node_S(const Node&, Stats&, long, long, bool*, double*);

// A node of the decomposition: terms over K free coordinates, in both
// arithmetics.  The float side proposes, the interval side decides.
struct Node {
    std::vector<int> sign;
    std::vector<long double> coeff;
    std::vector<std::vector<long double>> lam;   // [term][k]
    std::vector<Iv> icoeff;
    std::vector<std::vector<Iv>> ilam;
    std::size_t k = 0;
    // Exact-leaf bookkeeping: which label each free column refers to, the
    // current minimum exponent of every label, and the chamber's signs.  A
    // face split fixes a label at its current exponent; a tail split raises
    // one.  At a zero-dimensional leaf these data reconstruct the word
    // exactly, so the corner can be computed as an integer.
    std::vector<int> col_label;
    std::vector<int> expo;
    const std::vector<int>* chamber_signs = nullptr;
    std::vector<cyclo::Field::Elt> ecoeff;
    std::vector<std::vector<cyclo::Field::Elt>> elam;    // [term][k]
};

struct Stats {
    std::uint64_t amgm = 0, leaves = 0, nodes = 0, cap = 400000;
    std::uint64_t exact_leaves = 0;
    std::uint64_t abel = 0;
    bool root_lp_ok = false;   // did a root allocation exist at all
    double root_margin = 0.0;  // and with what worst-case slack
    std::array<std::uint64_t, 5> den{};
};

// Exact integer corner of the word given by (signs, exponents), by dynamic
// programming over coefficient matrices in the doubled fusion ring at level
// g_leaf_level.  This is the same computation as the repository's integral
// fusion leaf evaluations: an exact integer, so J = 0 certifies as
// nonnegative, which no strict interval test can do.
//
// A strictly negative return would be a counterexample to GKS2* itself and is
// reported loudly rather than treated as an ordinary failure.
bool leaf_exact_nonneg(const std::vector<int>& signs, const std::vector<int>& expo) {
    const int k = g_leaf_level;
    const int dim = k + 1;
    std::vector<mpz_t> grid(static_cast<std::size_t>(dim * dim));
    std::vector<mpz_t> next_grid(static_cast<std::size_t>(dim * dim));
    for (auto& z : grid) mpz_init(z);
    for (auto& z : next_grid) mpz_init(z);
    mpz_set_ui(grid[0], 1);
    auto idx = [&](int c, int d) { return static_cast<std::size_t>(c * dim + d); };
    for (std::size_t l = 0; l < signs.size(); ++l) {
        if (signs[l] == 0) continue;
        const int a = g_label_stride * (static_cast<int>(l) + 1);
        const int eps = signs[l];
        for (int e = 0; e < expo[l]; ++e) {
            for (auto& z : next_grid) mpz_set_ui(z, 0);
            for (int c = 0; c < dim; ++c)
                for (int d = 0; d < dim; ++d) {
                    if (mpz_sgn(grid[idx(c, d)]) == 0) continue;
                    const int chi = std::min(a + c, 2 * k - a - c);
                    for (int c2 = std::abs(a - c); c2 <= chi; c2 += 2)
                        mpz_add(next_grid[idx(c2, d)], next_grid[idx(c2, d)], grid[idx(c, d)]);
                    const int dhi = std::min(a + d, 2 * k - a - d);
                    for (int d2 = std::abs(a - d); d2 <= dhi; d2 += 2) {
                        if (eps > 0) mpz_add(next_grid[idx(c, d2)], next_grid[idx(c, d2)], grid[idx(c, d)]);
                        else mpz_sub(next_grid[idx(c, d2)], next_grid[idx(c, d2)], grid[idx(c, d)]);
                    }
                }
            std::swap(grid, next_grid);
        }
    }
    const int sgn = mpz_sgn(grid[0]);
    if (sgn < 0) {
        std::fprintf(stderr, "COUNTEREXAMPLE CANDIDATE: exact corner negative at a leaf\n");
        for (std::size_t l = 0; l < signs.size(); ++l)
            if (signs[l])
                std::fprintf(stderr, "  label %d sign %+d exponent %d\n",
                             g_label_stride * (static_cast<int>(l) + 1), signs[l], expo[l]);
    }
    for (auto& z : grid) mpz_clear(z);
    for (auto& z : next_grid) mpz_clear(z);
    return sgn >= 0;
}

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
    // Total-mass values to sweep.  S = 1 first, so behaviour is unchanged
    // wherever the ordinary certificate already works.
    // Off by default.  The sweep is mathematically sound but empirically
    // inert here: total positive supply is the number of positive terms and
    // demand at total mass S is S times the number of negatives, so S > 1 needs
    // #pos/#neg >= S, and 80% of root-infeasible regimes have that ratio below
    // 1.5.  Enabling it costs five times the linear-program solves for no
    // measured gain.  Retained behind a flag because the generalisation is
    // correct and may apply in a setting with more positive mass.
    static const long SN[5] = {1, 3, 2, 3, 4};
    static const long SD[5] = {1, 2, 1, 1, 1};
    const int nsweep = g_total_mass_sweep ? 5 : 1;
    for (int si = 0; si < nsweep; ++si) {
        if (amgm_node_S(nd, st, SN[si], SD[si], si == 0 ? proposal_existed : nullptr,
                        si == 0 ? margin_out : nullptr)) return true;
    }
    return false;
}

bool amgm_node_S(const Node& nd, Stats& st, long sn, long sd,
                 bool* proposal_existed, double* margin_out) {
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
        A.push_back(r); bb.push_back(static_cast<double>(sd) / static_cast<double>(sn));
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
        bb.push_back(-static_cast<double>(logl(nd.coeff[static_cast<std::size_t>(neg[xi])]))
                     + std::log(static_cast<double>(sn) / static_cast<double>(sd)));
    }
    std::vector<double> c(nv, 0.0); c[dj] = 1.0;
    std::vector<double> sol;
    if (!simplex(A, bb, c, static_cast<int>(nv), sol)) return false;
    // The allocation polytope's optimum is usually a face, not a point, and
    // which vertex the solver lands on decides whether rational rounding can
    // reproduce its tie structure.  On failure, retry from other optimal-face
    // vertices reached by small deterministic objective perturbations; the
    // verifier still decides every candidate, so this only varies the
    // proposer.
    const int kJitter = 24;

    // Never trust the solver: discard a proposal violating its own constraints.
    for (double q : sol) if (q < -1e-7) return false;
    for (std::size_t i = 0; i < A.size(); ++i) {
        double acc = 0.0;
        for (std::size_t j = 0; j < nv; ++j) acc += A[i][j] * sol[j];
        if (acc - bb[i] > 1e-6) return false;
    }
    // A zero-margin vertex is admissible: rows that hold with exact equality
    // are decided by the exact escalation in verification, so the proposal
    // gate must not demand strict slack.  Tied rays -- one negative whose
    // lambda exactly equals the best positive's -- are certified this way.
    if (sol[dj] < -1.0e-9) return false;
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
        H[xi * NC + K] = static_cast<double>(logl(nd.coeff[static_cast<std::size_t>(neg[xi])]))
                         - std::log(static_cast<double>(sn) / static_cast<double>(sd));
    }

    static const long ladder[5] = {100, 200, 500, 1000, 5000};
    for (int jitter = 0; jitter <= kJitter; ++jitter) {
    if (jitter > 0) {
        std::vector<double> cj = c;
        std::uint64_t rng = 0x243f6a8885a308d3ULL ^ (static_cast<std::uint64_t>(jitter) << 32);
        for (std::size_t v = 0; v + 1 < nv; ++v) {
            rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
            cj[v] = 1.0e-7 * static_cast<double>(rng % 1000);
        }
        std::vector<double> sj;
        if (!simplex(A, bb, cj, static_cast<int>(nv), sj)) continue;
        bool jfeas = true;
        for (double q : sj) if (q < -1e-7) { jfeas = false; break; }
        if (jfeas)
            for (std::size_t i = 0; i < A.size() && jfeas; ++i) {
                double acc2 = 0.0;
                for (std::size_t jx = 0; jx < nv; ++jx) acc2 += A[i][jx] * sj[jx];
                if (acc2 - bb[i] > 1e-6) jfeas = false;
            }
        if (!jfeas || sj[dj] < -1.0e-9) continue;
        sol = sj;
    }
    for (int di = 0; di < 5; ++di) {
        const long d = ladder[di];
        std::vector<long> cap_left(P, (sd * d) / sn);
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
        if (verify_alloc(nd.ilam, nd.icoeff, pos, neg, units, d, K, sn, sd,
                         nd.elam.empty() ? nullptr : &nd.elam,
                         nd.ecoeff.empty() ? nullptr : &nd.ecoeff)) {
            ++st.den[static_cast<std::size_t>(di)];
            return true;
        }
    }
    }
    return false;
}

// For a regime whose root allocation is infeasible, find which constraint
// family is responsible by re-solving with each dropped in turn.
//
//   (D) coordinatewise domination in each free coordinate
//   (K) the single constant condition on the coefficients
//
// If dropping (K) restores feasibility the obstruction is one of magnitude:
// the positives cannot pay the negative's constant even though they dominate
// its growth.  If dropping (D) restores it the obstruction is directional.
// If neither does, both bind and the regime is far from certifiable.
int diagnose_infeasible(const Node& nd) {
    std::vector<int> pos, neg;
    for (std::size_t t = 0; t < nd.sign.size(); ++t)
        (nd.sign[t] < 0 ? neg : pos).push_back(static_cast<int>(t));
    if (neg.empty() || pos.empty()) return 0;
    const std::size_t P = pos.size(), X = neg.size(), K = nd.k;
    const std::size_t nv = P * X + 1, dj = P * X;
    auto aidx = [&](std::size_t pi, std::size_t xi) { return pi * X + xi; };

    auto solve = [&](bool useD, bool useK) -> bool {
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
            if (useD)
                for (std::size_t kk = 0; kk < K; ++kk) {
                    std::vector<double> r(nv, 0.0);
                    for (std::size_t pi = 0; pi < P; ++pi)
                        r[aidx(pi, xi)] = -static_cast<double>(logl(nd.lam[static_cast<std::size_t>(pos[pi])][kk]));
                    r[dj] = 1.0;
                    A.push_back(r);
                    bb.push_back(-static_cast<double>(logl(nd.lam[static_cast<std::size_t>(neg[xi])][kk])));
                }
            if (useK) {
                std::vector<double> r(nv, 0.0);
                for (std::size_t pi = 0; pi < P; ++pi)
                    r[aidx(pi, xi)] = -static_cast<double>(logl(nd.coeff[static_cast<std::size_t>(pos[pi])]));
                r[dj] = 1.0;
                A.push_back(r);
                bb.push_back(-static_cast<double>(logl(nd.coeff[static_cast<std::size_t>(neg[xi])])));
            }
        }
        std::vector<double> c(nv, 0.0); c[dj] = 1.0;
        std::vector<double> sol;
        if (!simplex(A, bb, c, static_cast<int>(nv), sol)) return false;
        return sol[dj] > 0.0;
    };

    const bool dropK = solve(true, false);     // domination alone
    const bool dropD = solve(false, true);     // constant alone
    if (dropK && dropD) return 3;              // each alone fine, together not
    if (dropK) return 1;                       // constant condition blocks
    if (dropD) return 2;                       // domination blocks
    return 4;                                  // both block individually
}

// Exact Abel chain certificate.
//
// Group the node's terms by exact equality of their lambda vectors in
// Z[zeta].  A group's signed coefficient sum sigma_g is exact and linear --
// precisely what an AM-GM allocation cannot see, since a geometric mean of
// several tied positives is strictly below their sum.  If the groups form a
// chain under exact coordinatewise dominance, then for every exponent vector
// u >= 0 the group values x_g = prod lambda^u are ordered the same way, and
// Abel summation gives
//
//     sum_g sigma_g x_g = sum_g S_g (x_g - x_{g+1}) + S_G x_G,
//
// with every bracket nonnegative.  So exact nonnegativity of every prefix sum
// S_g certifies the node for all exponents simultaneously -- no tail bound,
// no finite leaves, any dimension.  The r -> infinity and r = 0 limits force
// the first and last prefix conditions, so on rays this is close to
// necessary as well.
bool abel_chain_certificate(const Node& nd) {
    if (nd.ecoeff.empty() || nd.elam.empty()) return false;
    const std::size_t T = nd.sign.size();
    // Group by exact lambda-vector equality.
    std::vector<std::vector<std::size_t>> groups;
    for (std::size_t t = 0; t < T; ++t) {
        bool placed = false;
        for (auto& g : groups) {
            bool same = true;
            for (std::size_t kk = 0; kk < nd.k && same; ++kk)
                if (!cyclo::Field::isZero(g_field.sub(nd.elam[t][kk], nd.elam[g[0]][kk])))
                    same = false;
            if (same) { g.push_back(t); placed = true; break; }
        }
        if (!placed) groups.push_back({t});
    }
    // Exact coordinatewise dominance must totally order the groups.
    const std::size_t G = groups.size();
    auto dominates = [&](std::size_t a, std::size_t b) {
        for (std::size_t kk = 0; kk < nd.k; ++kk)
            if (g_field.sign(g_field.sub(nd.elam[groups[a][0]][kk],
                                         nd.elam[groups[b][0]][kk])) < 0)
                return false;
        return true;
    };
    std::vector<std::size_t> order(G);
    for (std::size_t i = 0; i < G; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return dominates(a, b) && !dominates(b, a);
    });
    for (std::size_t i = 0; i + 1 < G; ++i)
        if (!dominates(order[i], order[i + 1])) return false;   // not a chain
    // Exact prefix sums of the signed coefficient numerators; the positive
    // rational scalar shared by every term cancels in each sum's sign.
    cyclo::Field::Elt prefix = g_field.fromInt(0);
    for (std::size_t i = 0; i < G; ++i) {
        for (std::size_t t : groups[order[i]]) {
            if (nd.sign[t] > 0) prefix = g_field.add(prefix, nd.ecoeff[t]);
            else prefix = g_field.sub(prefix, nd.ecoeff[t]);
        }
        if (g_field.sign(prefix) < 0) return false;
    }
    return true;
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
    if (nd.k == 0) {
        ++st.leaves;
        // Interval evaluation is the fast path; where it cannot decide --
        // exact zeros have enclosures straddling zero -- the integer DP is
        // the decider, not a fallback heuristic.
        if (leaf_nonneg(nd)) return true;
        if (nd.chamber_signs != nullptr && g_leaf_level > 0) {
            ++st.exact_leaves;
            return leaf_exact_nonneg(*nd.chamber_signs, nd.expo);
        }
        return false;
    }
    if (amgm_node(nd, st, depth == 0 ? &st.root_lp_ok : nullptr,
                  depth == 0 ? &st.root_margin : nullptr)) { ++st.amgm; return true; }
    if (abel_chain_certificate(nd)) { ++st.abel; return true; }
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
    face.expo = nd.expo;
    face.chamber_signs = nd.chamber_signs;
    face.col_label = nd.col_label;
    face.ecoeff = nd.ecoeff;
    face.elam.resize(nd.elam.size());
    for (std::size_t t = 0; t < nd.elam.size(); ++t)
        for (std::size_t kk = 0; kk < nd.k; ++kk)
            if (kk != best) face.elam[t].push_back(nd.elam[t][kk]);
    if (best < face.col_label.size()) face.col_label.erase(face.col_label.begin() + static_cast<std::ptrdiff_t>(best));
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
    tail.col_label = nd.col_label;
    tail.expo = nd.expo;
    tail.chamber_signs = nd.chamber_signs;
    tail.elam = nd.elam;
    tail.ecoeff.resize(nd.ecoeff.size());
    for (std::size_t t = 0; t < nd.ecoeff.size(); ++t)
        tail.ecoeff[t] = g_field.mulE(nd.ecoeff[t], nd.elam[t][best]);
    if (best < nd.col_label.size()) tail.expo[static_cast<std::size_t>(nd.col_label[best])] += 1;
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
    int level = 0;
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
        else if (!std::strcmp(argv[i], "--total-mass")) g_total_mass_sweep = true;
        else if (!std::strcmp(argv[i], "--level") && i + 1 < argc) level = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--certify-all")) g_certify_all = true;
    }
    if (level > 0) {
        if (level < 2) { std::fprintf(stderr, "level must be at least 2\n"); return 2; }
        g_full_level = true;
        // A closed level cannot rest on the tolerance-based Hall pre-filter,
        // and at higher levels its bitmasks overflow anyway.
        g_certify_all = true;
    } else if (rank < 3) { std::fprintf(stderr, "rank must be at least 3\n"); return 2; }
    if (threads < 1) threads = 1;

    Model model;
    if (g_full_level) build_model_level(model, level);
    else build_model(model, rank);
    g_leaf_level = g_full_level ? level : (2 * rank - 1);
    g_label_stride = g_full_level ? 1 : 2;
    g_field.init(g_full_level ? level + 2 : 2 * rank + 1);
    {
        // Exact spectral tables: node values as Dickson polynomials of zeta,
        // characters as Chebyshev values of those.  Cross-checked against the
        // long double model so all three spectral models must agree.
        const int nnodes = model.rank;
        g_exact_y.assign(static_cast<std::size_t>(nnodes), {});
        g_exact_chi.assign(static_cast<std::size_t>(nnodes), {});
        for (int j = 0; j < nnodes; ++j) {
            const int mstep = g_full_level ? (j + 1) : (j + 1);
            g_exact_y[static_cast<std::size_t>(j)] = g_field.dickson(mstep);
            g_exact_chi[static_cast<std::size_t>(j)].resize(static_cast<std::size_t>(model.labels));
            for (int l = 0; l < model.labels; ++l) {
                const int a = g_label_stride * (l + 1);
                g_exact_chi[static_cast<std::size_t>(j)][static_cast<std::size_t>(l)] =
                    g_field.chebyshevS(a, g_exact_y[static_cast<std::size_t>(j)]);
            }
        }
        int agree = 0;
        for (int j = 0; j < nnodes; ++j)
            for (int l = 0; l < model.labels; ++l) {
                // Evaluate the exact element crudely in double via Horner at
                // zeta ~ 2cos(pi/n); roundoff is far below the 1e-6 gate for
                // these small-height values, and the exact/interval agreement
                // is enforced separately by the tied-row escalation itself.
                const double z = 2.0 * std::cos(3.141592653589793238 / static_cast<double>(g_field.n));
                const auto& e = g_exact_chi[static_cast<std::size_t>(j)][static_cast<std::size_t>(l)];
                double v = 0.0;
                for (std::size_t i = e.size(); i-- > 0;) v = v * z + e[i].get_d();
                if (std::fabs(v - static_cast<double>(model.v[static_cast<std::size_t>(j)][static_cast<std::size_t>(l)])) > 1e-6) {
                    std::fprintf(stderr, "exact spectral cross-check failed: node %d label %d\n", j, l);
                    return 3;
                }
                ++agree;
            }
        std::printf("  exact_crosscheck %d values agree (field degree %d)\n", agree, g_field.D);
    }
    const int L = model.labels;
    const int nodes = model.rank;

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
        for (int i = 0; i < nodes; ++i) {
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
    std::uint64_t direct = 0, residual = 0, pointwise = 0, sep = 0, parity_zero = 0;
    std::uint64_t certified_total = 0;
    std::array<std::uint64_t, 5> den_used{};
    std::uint64_t sound_regimes = 0, sound_points = 0, sound_viol = 0;
    std::uint64_t ctl_tried = 0, ctl_refused = 0;
    std::uint64_t flat_total = 0, split_total = 0;
    std::uint64_t node_total = 0, budget_hit = 0, exact_leaf_total = 0, abel_total = 0;

    const int supports = static_cast<int>(std::llround(std::pow(3.0, L)));

    auto worker = [&]() {
        std::uint64_t ld = 0, lr = 0, lp = 0, ls = 0, lcert = 0, lpz = 0;
        std::array<std::uint64_t, 5> lden{};
        std::uint64_t lsound_regimes = 0, lsound_points = 0, lsound_viol = 0;
        std::uint64_t lflat = 0, lsplit = 0, lnodes = 0, lbudget = 0, lexact = 0, label_ = 0;
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
                if (g_full_level) {
                    // Parity selection rule (proved in
                    // SU2_PARITY_SELECTION_RULE_2026_07_24.md): the corner
                    // vanishes identically unless the total degree
                    // sum_a a*p_a is even.  Residual increments change it by
                    // 2a, so the parity is fixed per chamber, and odd-parity
                    // chambers are exactly zero rather than certifiable.
                    int deg = 0;
                    for (std::size_t i = 0; i < act.size(); ++i) {
                        const int l = act[i];
                        if (((l + 1) & 1) != 0) deg += powers[static_cast<std::size_t>(l)];
                    }
                    if (deg & 1) { lpz += static_cast<std::uint64_t>(1) << act.size(); continue; }
                }
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
                    if (!g_certify_all && pos.size() <= 31 && neg.size() <= 22)
                    for (std::size_t x = 0; x < neg.size(); ++x)
                        for (std::size_t p = 0; p < pos.size(); ++p) {
                            bool ok = true;
                            for (int l : freev)
                                if (!ge(terms[static_cast<std::size_t>(pos[p])].lam[static_cast<std::size_t>(l)],
                                        terms[static_cast<std::size_t>(neg[x])].lam[static_cast<std::size_t>(l)])) { ok = false; break; }
                            if (ok) edge[x] |= (1U << p);
                        }
                    bool hall = false;
                    if (!g_certify_all && pos.size() <= 31 && neg.size() <= 22) {
                    hall = true;
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
                    }
                    if (hall) { ++ld; continue; }
                    ++lr;

                    // Build the root node and certify it.  With --decompose
                    // off the depth is zero, so this is exactly the flat
                    // allocation stage and must reproduce its numbers.
                    const std::size_t K = freev.size();
                    Node root;
                    root.k = K;
                    root.col_label = freev;
                    root.expo = powers;
                    root.chamber_signs = &signs;
                    root.ecoeff.resize(terms.size());
                    root.elam.assign(terms.size(), {});
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
                        cyclo::Field::Elt ec = terms[t].ecoeff;
                        for (std::size_t kk = 0; kk < K; ++kk) {
                            const std::size_t l = static_cast<std::size_t>(freev[kk]);
                            root.elam[t].push_back(terms[t].elam[l]);
                            ec = g_field.mulE(ec, terms[t].elam[l]);
                        }
                        root.ecoeff[t] = ec;
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
                                        " root_lp=%s margin=%.3e block=%d nodes=%llu\n",
                                        support, parity, K, nneg, root.sign.size() - nneg,
                                        st.root_lp_ok ? "feasible" : "infeasible", st.root_margin,
                                        st.root_lp_ok ? 0 : diagnose_infeasible(root),
                                        static_cast<unsigned long long>(st.nodes));
                        }
                        continue;
                    }
                    ++lcert;
                    if (st.amgm == 1 && st.leaves == 0) ++lflat; else ++lsplit;
                    lnodes += st.nodes;
                    lexact += st.exact_leaves;
                    label_ += st.abel;
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
                            long double tot = 0.0L, mass = 0.0L;
                            for (std::size_t t = 0; t < terms.size(); ++t) {
                                long double v = terms[t].coeff;
                                for (std::size_t kk = 0; kk < K; ++kk)
                                    for (int q = 0; q < expo[kk]; ++q)
                                        v *= terms[t].lam[static_cast<std::size_t>(freev[kk])];
                                tot += static_cast<long double>(terms[t].sign) * v;
                                mass += fabsl(v);
                            }
                            ++lsound_points;
                            // Exact-zero points are certified now, so the
                            // sampler visits sums that cancel exactly; its
                            // roundoff scales with the cancelling mass, not
                            // with the tiny result.
                            if (tot < -1.0e-16L * (1.0L + mass)) ++lsound_viol;
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
                            // With a single negative, giving it all of one
                            // positive's capacity is a legitimate allocation,
                            // not a corruption, so that control applies only
                            // when at least two negatives compete.
                            if (cneg.size() >= 2) {
                                for (std::size_t xi = 0; xi < cneg.size(); ++xi) u[0][xi] = d;
                                ++lctl_tried;                          // capacity blown
                                if (!verify_alloc(root.ilam, root.icoeff, cpos, cneg, u, d, K)) ++lctl_refused;
                            }
                            for (std::size_t xi = 0; xi < cneg.size(); ++xi) u[0][xi] = d + 1;
                            ++lctl_tried;                              // normalisation broken
                            if (!verify_alloc(root.ilam, root.icoeff, cpos, cneg, u, d, K)) ++lctl_refused;
                        }
                    }
                }
            }
        }
        std::lock_guard<std::mutex> g(mu);
        direct += ld; residual += lr; pointwise += lp; sep += ls; parity_zero += lpz;
        certified_total += lcert;
        for (std::size_t i = 0; i < 5; ++i) den_used[i] += lden[i];
        sound_regimes += lsound_regimes; sound_points += lsound_points; sound_viol += lsound_viol;
        ctl_tried += lctl_tried; ctl_refused += lctl_refused;
        flat_total += lflat; split_total += lsplit;
        node_total += lnodes; budget_hit += lbudget; exact_leaf_total += lexact; abel_total += label_;
    };

    std::vector<std::thread> pool;
    for (int i = 0; i < threads; ++i) pool.emplace_back(worker);
    for (auto& t : pool) t.join();

    if (g_full_level)
        std::printf("SU2_LEVEL_AMGM_CPP level=%d labels=%d nodes=%d threads=%d certify_all=1\n",
                    level, L, nodes, threads);
    else
        std::printf("SU2_ORBIT_AMGM_CPP rank=%d level=%d threads=%d\n", rank, 2 * rank - 1, threads);
    if (g_full_level)
        std::printf("  parity_zero      %llu\n", static_cast<unsigned long long>(parity_zero));
    std::printf("  direct_hall      %llu\n", static_cast<unsigned long long>(direct));
    std::printf("  residual         %llu\n", static_cast<unsigned long long>(residual));
    std::printf("  pointwise        %llu\n", static_cast<unsigned long long>(pointwise));
    std::printf("  certified        %llu\n", static_cast<unsigned long long>(certified_total));
    std::printf("  by_flat_amgm     %llu\n", static_cast<unsigned long long>(flat_total));
    std::printf("  by_split         %llu\n", static_cast<unsigned long long>(split_total));
    std::printf("  split_nodes      %llu\n", static_cast<unsigned long long>(node_total));
    std::printf("  budget_exhausted %llu\n", static_cast<unsigned long long>(budget_hit));
    std::printf("  exact_leaves     %llu\n", static_cast<unsigned long long>(exact_leaf_total));
    std::printf("  abel_chain       %llu\n", static_cast<unsigned long long>(abel_total));
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
    if (g_full_level) {
        const std::uint64_t nonzero = residual + pointwise;
        const std::uint64_t closed = certified_total + pointwise;
        std::printf("  regimes_nonzero  %llu\n", static_cast<unsigned long long>(nonzero));
        std::printf("  open             %llu\n", static_cast<unsigned long long>(nonzero - closed));
        if (nonzero)
            std::printf("  coverage         %.1f%% of nonzero regimes\n",
                        100.0 * static_cast<double>(closed) / static_cast<double>(nonzero));
    } else if (residual)
        std::printf("  coverage         %.1f%% of residual regimes\n",
                    100.0 * static_cast<double>(certified_total) / static_cast<double>(residual));
    return 0;
}
