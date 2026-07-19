#if __has_include(<mpfr.h>)
#include <mpfr.h>
#else
#include <gmp.h>
extern "C" {
typedef long mpfr_prec_t;
typedef long mpfr_exp_t;
typedef int mpfr_sign_t;
typedef enum { MPFR_RNDN = 0, MPFR_RNDZ = 1, MPFR_RNDU = 2, MPFR_RNDD = 3,
               MPFR_RNDA = 4, MPFR_RNDF = 5, MPFR_RNDNA = -1 } mpfr_rnd_t;
typedef struct { mpfr_prec_t _mpfr_prec; mpfr_sign_t _mpfr_sign; mpfr_exp_t _mpfr_exp;
                 mp_limb_t* _mpfr_d; } __mpfr_struct;
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
int mpfr_div(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_div_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
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
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr mpfr_prec_t PRECISION = 384;
constexpr int CUTOFF = 296;

struct Real {
  mpfr_t value;
  Real() { mpfr_init2(value, PRECISION); mpfr_set_zero(value, 0); }
  Real(const Real& other) { mpfr_init2(value, PRECISION); mpfr_set(value, other.value, MPFR_RNDN); }
  Real& operator=(const Real& other) {
    if (this != &other) mpfr_set(value, other.value, MPFR_RNDN);
    return *this;
  }
  Real(Real&& other) noexcept { mpfr_init2(value, PRECISION); mpfr_swap(value, other.value); }
  Real& operator=(Real&& other) noexcept {
    if (this != &other) mpfr_swap(value, other.value);
    return *this;
  }
  ~Real() { mpfr_clear(value); }
};

struct Interval { Real lo, hi; };

Interval rational(unsigned long numerator, unsigned long denominator = 1) {
  Interval out;
  mpfr_set_ui(out.lo.value, numerator, MPFR_RNDD);
  mpfr_div_ui(out.lo.value, out.lo.value, denominator, MPFR_RNDD);
  mpfr_set_ui(out.hi.value, numerator, MPFR_RNDU);
  mpfr_div_ui(out.hi.value, out.hi.value, denominator, MPFR_RNDU);
  return out;
}

Interval add(const Interval& a, const Interval& b) {
  Interval out;
  mpfr_add(out.lo.value, a.lo.value, b.lo.value, MPFR_RNDD);
  mpfr_add(out.hi.value, a.hi.value, b.hi.value, MPFR_RNDU);
  return out;
}

Interval neg(const Interval& a) {
  Interval out;
  mpfr_neg(out.lo.value, a.hi.value, MPFR_RNDD);
  mpfr_neg(out.hi.value, a.lo.value, MPFR_RNDU);
  return out;
}

Interval sub(const Interval& a, const Interval& b) { return add(a, neg(b)); }

void assign_min(Real& target, const Real& candidate) {
  if (mpfr_cmp(candidate.value, target.value) < 0) mpfr_set(target.value, candidate.value, MPFR_RNDN);
}
void assign_max(Real& target, const Real& candidate) {
  if (mpfr_cmp(candidate.value, target.value) > 0) mpfr_set(target.value, candidate.value, MPFR_RNDN);
}

Interval mul(const Interval& a, const Interval& b) {
  Real lower[4], upper[4];
  const mpfr_srcptr av[2] = {a.lo.value, a.hi.value};
  const mpfr_srcptr bv[2] = {b.lo.value, b.hi.value};
  int k = 0;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j, ++k) {
      mpfr_mul(lower[k].value, av[i], bv[j], MPFR_RNDD);
      mpfr_mul(upper[k].value, av[i], bv[j], MPFR_RNDU);
    }
  Interval out;
  out.lo = lower[0]; out.hi = upper[0];
  for (int i = 1; i < 4; ++i) { assign_min(out.lo, lower[i]); assign_max(out.hi, upper[i]); }
  return out;
}

Interval reciprocal_positive(const Interval& a) {
  if (mpfr_sgn(a.lo.value) <= 0) throw std::runtime_error("nonpositive interval divisor");
  const Interval one = rational(1);
  Interval out;
  mpfr_div(out.lo.value, one.lo.value, a.hi.value, MPFR_RNDD);
  mpfr_div(out.hi.value, one.hi.value, a.lo.value, MPFR_RNDU);
  return out;
}

Interval div(const Interval& a, const Interval& b) { return mul(a, reciprocal_positive(b)); }

Interval log_interval(const Interval& a) {
  if (mpfr_sgn(a.lo.value) <= 0) throw std::runtime_error("log of nonpositive interval");
  Interval out;
  mpfr_log(out.lo.value, a.lo.value, MPFR_RNDD);
  mpfr_log(out.hi.value, a.hi.value, MPFR_RNDU);
  return out;
}

Interval exp_interval(const Interval& a) {
  Interval out;
  mpfr_exp(out.lo.value, a.lo.value, MPFR_RNDD);
  mpfr_exp(out.hi.value, a.hi.value, MPFR_RNDU);
  return out;
}

Interval lngamma_interval(const Interval& a) {
  Interval out;
  mpfr_lngamma(out.lo.value, a.lo.value, MPFR_RNDD);
  mpfr_lngamma(out.hi.value, a.hi.value, MPFR_RNDU);
  return out;
}

Interval square(const Interval& a) { return mul(a, a); }

std::string format(const Real& value) {
  char* raw = nullptr;
  if (mpfr_asprintf(&raw, "%.28Rg", value.value) < 0 || raw == nullptr) return "<format-error>";
  std::string out(raw);
  mpfr_free_str(raw);
  return out;
}

Interval one_trace_log(int n) {
  const Interval log_2n = log_interval(rational(2UL * n));
  return sub(
      add(log_interval(rational(5)), div(log_interval(log_2n), rational(4))),
      mul(rational(n), sub(log_2n, rational(1))));
}

Interval unitary_trace_log(int n) {
  const Interval log_n = log_interval(rational(n));
  return sub(
      add(add(log_interval(rational(19)), div(log_n, rational(4))),
          div(log_interval(log_n), rational(2))),
      lngamma_interval(rational(n + 2)));
}

Interval logadd(const Interval& first, const Interval& second) {
  return log_interval(add(exp_interval(first), exp_interval(second)));
}

Interval logsub_positive(const Interval& first, const Interval& second) {
  const Interval difference = sub(exp_interval(first), exp_interval(second));
  if (mpfr_sgn(difference.lo.value) <= 0)
    throw std::runtime_error("interval log subtraction is not provably positive");
  return log_interval(difference);
}

struct Checks {
  Interval gaussian, positive_tv, tau, unitary_vs_orth, r0, r1, r2;
  Interval gaussian_derivative, positive_tv_derivative, square_r0_derivative,
      square_r1_derivative, r2_derivative;
};

Checks evaluate(int c_value) {
  const Interval c = rational(c_value);
  const Interval r = rational(2000001, 1000000);
  const Interval eta = rational(359, 2000);
  const Interval a = rational(4571, 2000);
  const Interval l = rational(2269251, 500000);
  const Interval log2 = log_interval(rational(2));
  const Interval log_half = neg(log2);
  const Interval eta_c_minus_one = sub(mul(eta, c), rational(1));
  const Interval gaussian_tau = neg(div(square(eta_c_minus_one), rational(4)));
  const Interval tau = logadd(gaussian_tau,
                              add(log2, one_trace_log(c_value / 2)));
  const Interval den = sub(rational(1),
                           reciprocal_positive(mul(square(eta), square(c))));

  Checks out;
  out.gaussian = sub(
      sub(add(mul(sub(a, div(l, rational(2))), c),
              div(log_interval(mul(l, c)), rational(2))),
          log_interval(add(mul(l, c), rational(1)))),
      add(div(log2, rational(2)), log2));
  out.positive_tv = add(one_trace_log(c_value), mul(a, c));
  out.tau = tau;
  out.unitary_vs_orth = sub(unitary_trace_log(c_value), one_trace_log(c_value / 2));
  out.r0 = sub(add(mul(a, c), tau), log_half);

  const Interval r1_inside = div(
      mul(mul(rational(4), add(rational(1), reciprocal_positive(square(c)))),
          mul(div(r, rational(2)), square(div(r, rational(2))))),
      mul(square(r), den));
  out.r1 = sub(
      add(add(add(log_interval(r1_inside), mul(a, c)), tau),
          mul(c, log_interval(div(rational(2), r)))),
      log_half);

  const Interval r_over_eta = div(r, eta);
  const Interval r2_inside = div(
      mul(rational(4), mul(r_over_eta, square(r_over_eta))),
      mul(mul(square(r), square(c)), den));
  out.r2 = sub(
      sub(log_interval(r2_inside), mul(sub(log_interval(r_over_eta), a), c)),
      log_half);

  out.gaussian_derivative = sub(sub(a, div(l, rational(2))),
                                reciprocal_positive(mul(rational(2), c)));
  out.positive_tv_derivative = add(
      sub(a, log_interval(mul(rational(2), c))),
      reciprocal_positive(mul(mul(rational(4), c),
                              log_interval(mul(rational(2), c)))));
  out.square_r0_derivative = sub(a, div(mul(eta, eta_c_minus_one), rational(2)));
  out.square_r1_derivative = add(out.square_r0_derivative,
                                 log_interval(div(rational(2), r)));
  out.r2_derivative = sub(a, log_interval(r_over_eta));
  return out;
}

bool positive(const Interval& value) { return mpfr_sgn(value.lo.value) > 0; }
bool nonnegative(const Interval& value) { return mpfr_sgn(value.lo.value) >= 0; }
bool negative(const Interval& value) { return mpfr_sgn(value.hi.value) < 0; }
bool nonpositive(const Interval& value) { return mpfr_sgn(value.hi.value) <= 0; }

bool passes(const Checks& c) {
  return nonnegative(c.gaussian) && nonpositive(c.positive_tv) &&
         nonpositive(c.unitary_vs_orth) && nonpositive(c.r0) &&
         nonpositive(c.r1) && nonpositive(c.r2) &&
         positive(c.gaussian_derivative) && negative(c.r2_derivative);
}

Interval first_hit_margin(int c_value, int odd_n) {
  // Directed-interval version of finite_trace_pushforward_margin with the
  // source-audited first-hit parameters r=20001/10000, eta=1861/10000.
  const Interval c = rational(c_value);
  const Interval odd = rational(odd_n);
  const Interval r = rational(20001, 10000);
  const Interval eta = rational(1861, 10000);
  const Interval ell = rational(9117, 2000);  // 2r+3eta
  const Interval log2 = log_interval(rational(2));
  const Interval eta_c_minus_one = sub(mul(eta, c), rational(1));
  const Interval den = sub(rational(1), reciprocal_positive(mul(square(eta), square(c))));
  if (mpfr_sgn(den.lo.value) <= 0)
    throw std::runtime_error("first-hit Chebyshev denominator is not positive");

  const Interval x_square = mul(ell, c);
  const Interval gaussian_log = sub(
      sub(add(neg(div(x_square, rational(2))),
              div(log_interval(x_square), rational(2))),
          log_interval(add(x_square, rational(1)))),
      div(log2, rational(2)));
  const Interval positive_tv_log = one_trace_log(c_value);
  const Interval p1_log = logsub_positive(gaussian_log, positive_tv_log);

  const Interval square_gaussian_log = neg(div(square(eta_c_minus_one), rational(4)));
  const Interval square_tv_log = add(log2, one_trace_log(c_value / 2));
  const Interval p2_log = logadd(square_gaussian_log, square_tv_log);
  const Interval positive_gap_log = logsub_positive(p1_log, p2_log);

  const Interval left_log = add(
      log_interval(mul(mul(square(r), square(c)), den)), positive_gap_log);
  const Interval term1_log = add(
      add(log_interval(mul(rational(2), add(square(c), rational(1)))), p2_log),
      mul(odd, log_interval(div(rational(2), r))));
  const Interval term2_log = add(
      log2, mul(odd, log_interval(div(eta, r))));
  return sub(left_log, logadd(term1_log, term2_log));
}

Interval rains_correction(int c_value) {
  const int n_value = std::max(2, c_value / 2);
  const Interval n = rational(n_value);
  const Interval exponent = neg(add(
      mul(rational(5, 2), log_interval(n)), div(n, rational(2))));
  return exp_interval(exponent);
}

Interval chernoff_tau_log(int c_value, const Interval& eta) {
  const Interval c = rational(c_value);
  const Interval threshold = sub(mul(eta, c), rational(1));
  const Interval denominator =
      mul(rational(16), add(rational(1), rains_correction(c_value)));
  return sub(log_interval(rational(2)), div(square(threshold), denominator));
}

void require_chernoff_window(int c_value, const Interval& eta) {
  const int n_value = std::max(2, c_value / 2);
  const Interval threshold = sub(mul(eta, rational(c_value)), rational(1));
  if (mpfr_sgn(threshold.lo.value) <= 0)
    throw std::runtime_error("Chernoff threshold is not positive");
  const Interval window = mul(
      mul(mul(rational(4), exp_interval(neg(rational(3, 2)))),
          add(rational(1), rains_correction(c_value))),
      rational(n_value));
  if (mpfr_sgn(sub(window, threshold).lo.value) <= 0)
    throw std::runtime_error("Chernoff threshold is outside the source window");
}

Interval chernoff_margin(int c_value, int odd_n, const Interval& r,
                         const Interval& eta) {
  if (c_value < 132) throw std::runtime_error("CJL one-trace source range failed");
  require_chernoff_window(c_value, eta);
  const Interval c = rational(c_value);
  const Interval odd = rational(odd_n);
  const Interval ell = add(mul(rational(2), r), mul(rational(3), eta));
  const Interval log2 = log_interval(rational(2));
  const Interval den = sub(rational(1), reciprocal_positive(mul(square(eta), square(c))));
  if (mpfr_sgn(den.lo.value) <= 0)
    throw std::runtime_error("Chernoff Chebyshev denominator is not positive");

  const Interval x_square = mul(ell, c);
  const Interval gaussian_log = sub(
      sub(add(neg(div(x_square, rational(2))),
              div(log_interval(x_square), rational(2))),
          log_interval(add(x_square, rational(1)))),
      div(log2, rational(2)));
  const Interval p1_log = logsub_positive(gaussian_log, one_trace_log(c_value));
  const Interval p2_log = chernoff_tau_log(c_value, eta);
  const Interval positive_gap_log = logsub_positive(p1_log, p2_log);
  const Interval left_log = add(
      log_interval(mul(mul(square(r), square(c)), den)), positive_gap_log);
  const Interval term1_log = add(
      add(log_interval(mul(rational(2), add(square(c), rational(1)))), p2_log),
      mul(odd, log_interval(div(rational(2), r))));
  const Interval term2_log = add(
      log2, mul(odd, log_interval(div(eta, r))));
  return sub(left_log, logadd(term1_log, term2_log));
}

struct FirstHitRow {
  char family;
  int rank;
};

void print_interval(const char* label, const Interval& value);

std::vector<FirstHitRow> first_hit_rows() {
  std::vector<FirstHitRow> rows;
  for (char family : {'B', 'C'}) {
    rows.push_back({family, 276});
    for (int rank = 278; rank <= 295; ++rank) rows.push_back({family, rank});
  }
  rows.push_back({'D', 276});
  rows.push_back({'D', 278});
  for (int rank = 280; rank <= 295; ++rank) rows.push_back({'D', rank});
  rows.push_back({'D', 297});
  return rows;
}

void verify_first_hit_rows() {
  const std::vector<FirstHitRow> rows = first_hit_rows();
  if (rows.size() != 57) throw std::runtime_error("first-hit row count is not 57");
  Interval minimum;
  std::string minimum_label;
  bool have_minimum = false;
  for (const FirstHitRow& row : rows) {
    const int c_value = row.family == 'D' && (row.rank & 1) ? row.rank - 2 : row.rank;
    const int odd_n = row.family == 'D'
                          ? ((row.rank & 1) ? row.rank - 2 : row.rank - 1)
                          : 2 * row.rank + 1;
    const Interval margin = first_hit_margin(c_value, odd_n);
    if (!positive(margin))
      throw std::runtime_error("first-hit interval margin failed at " +
                               std::string(1, row.family) + "_" +
                               std::to_string(row.rank));
    if (!have_minimum || mpfr_cmp(margin.lo.value, minimum.lo.value) < 0) {
      minimum = margin;
      minimum_label = std::string(1, row.family) + "_" + std::to_string(row.rank);
      have_minimum = true;
    }
  }
  if (!have_minimum || minimum_label != "D_281")
    throw std::runtime_error("first-hit minimum row mismatch");
  std::cout << "first_hit_rows=57 minimum_row=" << minimum_label << '\n';
  print_interval("first_hit_minimum_log_margin", minimum);
}

void verify_chernoff_rows() {
  Interval minimum;
  std::string minimum_label;
  bool have_minimum = false;
  int checked = 0;
  for (char family : {'B', 'C'}) {
    for (int rank = 218; rank <= 277; ++rank) {
      if (rank == 276) continue;
      if (rank > 275 && rank != 277) continue;
      Interval r, eta;
      if (rank == 219) {
        r = rational(99, 50);
        eta = rational(56, 125);
      } else if ((rank & 1) == 0 && rank <= 266) {
        r = rational(2001, 1000);
        eta = rational(9, 20);
      } else {
        r = rational(10001, 5000);
        eta = rational(56, 125);
      }
      const Interval margin = chernoff_margin(rank, 2 * rank + 1, r, eta);
      if (!positive(margin))
        throw std::runtime_error("Chernoff interval margin failed at " +
                                 std::string(1, family) + "_" +
                                 std::to_string(rank));
      const std::string label = std::string(1, family) + "_" + std::to_string(rank);
      if (!have_minimum || mpfr_cmp(margin.lo.value, minimum.lo.value) < 0) {
        minimum = margin;
        minimum_label = label;
        have_minimum = true;
      }
      ++checked;
    }
  }
  if (checked != 118 || !have_minimum || minimum_label != "B_219")
    throw std::runtime_error("Chernoff row ledger mismatch");
  std::cout << "chernoff_rows=118 minimum_row=" << minimum_label << '\n';
  print_interval("chernoff_minimum_log_margin", minimum);
}

void print_interval(const char* label, const Interval& value) {
  std::cout << label << " in [" << format(value.lo) << ", " << format(value.hi) << "]\n";
}
}  // namespace

int main() {
  try {
    int first_bad = -1, last_bad = -1;
    for (int c = 248; c < CUTOFF; ++c) {
      if (!passes(evaluate(c))) {
        if (first_bad < 0) first_bad = c;
        last_bad = c;
      }
    }
    const Checks cutoff = evaluate(CUTOFF);
    if (!passes(cutoff)) throw std::runtime_error("cutoff row failed");
    if (!negative(cutoff.positive_tv_derivative) ||
        !negative(cutoff.square_r0_derivative) ||
        !negative(cutoff.square_r1_derivative))
      throw std::runtime_error("cutoff derivative check failed");
    for (int c = CUTOFF; c <= 10000; ++c)
      if (!passes(evaluate(c)))
        throw std::runtime_error("integer replay failed at C=" + std::to_string(c));
    verify_first_hit_rows();
    verify_chernoff_rows();

    std::cout << "Rains-square trace cutoff directed-MPFR certificate\n";
    std::cout << "precision_bits=" << PRECISION << " mpfr_version=" << mpfr_get_version()
              << " cutoff=" << CUTOFF << '\n';
    print_interval("gaussian_minus_required", cutoff.gaussian);
    print_interval("positive_tv_log_plus_A_C", cutoff.positive_tv);
    print_interval("unitary_projection_minus_half_orth", cutoff.unitary_vs_orth);
    print_interval("log_R0_minus_log_half", cutoff.r0);
    print_interval("log_R1_minus_log_half", cutoff.r1);
    print_interval("log_R2_minus_log_half", cutoff.r2);
    print_interval("gaussian_derivative", cutoff.gaussian_derivative);
    print_interval("log_R2_derivative", cutoff.r2_derivative);
    std::cout << "pre_cutoff_bad_range=" << first_bad << ".." << last_bad << '\n';
    std::cout << "integer replay OK for 296 <= C <= 10000\n";
    std::cout << "SUMMARY failures=0 directed_rounding=1\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "Rains-square MPFR certificate failure: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
