#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Polynomial = std::vector<mpz_class>;

struct ExpectedRow {
  int rank;
  int moment;
  const char* delta;
};

const ExpectedRow kExpected[] = {
    {6, 39, "64535719854392029544773326762011774756925778563684823690382903928732901693552607059275"},
    {7, 41, "81983250005769227285498956748387927182179586638327455931621250522805530412550262461119689354176"},
    {8, 41, "13776730806844857067459723475572806115606561762537749689970379383403774177954665948696550730822960"},
    {9, 43, "516741713808588211599059495979604717390970752707423594032209037987274754148154730954160916056953963450368"},
    {10, 43, "3618232844122202584832213951704951536993240226634987294481850293593859749135402230512808610277609410270817"},
    {11, 45, "33399774449175848055011934469656070323910964218709116241112806497510841469515062121126196910235665905718776943228"},
    {12, 45, "60567314448754365018106205779950325082006402099956453971960484870790769677108940276606446856248401017048182106706"},
    {13, 47, "363239305623168787675421494307053710143238482046635190860778910231693491256558565289022469174629997636148032811272314736"},
    {14, 47, "412017453922855017028903546300156146196642966415729406173097639746151298064337858733695094575907116122702045750587973615"},
    {15, 49, "2454353316530929119018984347889491575535943771589993590419065745590973632964833686716148957357041276645873773505319638231954239"},
    {16, 49, "2493119549736314583261850567832642447348332533596249440198341404887202602996116002199771825733723134557950156991608933086839975"},
    {17, 51, "16910085949149500902700908171041557824229854633053201817574276758325371717389977268687410623596430698509784622587525716043509684517756"},
    {18, 51, "16928441762086713834411725710317433241255361224654439216646167863265666371228166615581072054953650607845481890940623126734747743109775"},
};

[[noreturn]] void fail(const std::string& message) {
  throw std::runtime_error(message);
}

mpz_class binomial(unsigned n, unsigned k) {
  if (k > n) return 0;
  mpz_class out;
  mpz_bin_uiui(out.get_mpz_t(), n, k);
  return out;
}

Polynomial polynomial_add(const Polynomial& left, const Polynomial& right,
                          int right_sign = 1) {
  Polynomial out(std::max(left.size(), right.size()));
  for (std::size_t i = 0; i < left.size(); ++i) out[i] += left[i];
  for (std::size_t i = 0; i < right.size(); ++i) out[i] += right_sign * right[i];
  while (out.size() > 1 && out.back() == 0) out.pop_back();
  return out;
}

Polynomial polynomial_multiply(const Polynomial& left, const Polynomial& right) {
  Polynomial out(left.size() + right.size() - 1);
  for (std::size_t i = 0; i < left.size(); ++i) {
    for (std::size_t j = 0; j < right.size(); ++j) out[i + j] += left[i] * right[j];
  }
  while (out.size() > 1 && out.back() == 0) out.pop_back();
  return out;
}

Polynomial polynomial_power(Polynomial base, unsigned exponent) {
  Polynomial out{1};
  while (exponent != 0) {
    if (exponent & 1U) out = polynomial_multiply(out, base);
    exponent >>= 1U;
    if (exponent != 0) base = polynomial_multiply(base, base);
  }
  return out;
}

Polynomial polynomial_shift(const Polynomial& polynomial, const mpz_class& shift) {
  // Horner composition polynomial(x + shift).
  Polynomial out{0};
  const Polynomial x_plus_shift{shift, 1};
  for (auto it = polynomial.rbegin(); it != polynomial.rend(); ++it) {
    out = polynomial_multiply(out, x_plus_shift);
    out[0] += *it;
  }
  return out;
}

void require_positive_coefficients(const Polynomial& polynomial,
                                   const std::string& label) {
  if (polynomial.empty()) fail(label + " is empty");
  for (const mpz_class& coefficient : polynomial) {
    if (coefficient <= 0) fail(label + " has a nonpositive shifted coefficient");
  }
}

mpq_class rational_power(mpq_class base, unsigned exponent) {
  mpq_class out = 1;
  while (exponent != 0) {
    if (exponent & 1U) out *= base;
    exponent >>= 1U;
    if (exponent != 0) base *= base;
  }
  out.canonicalize();
  return out;
}

std::vector<mpz_class> factorials(int maximum) {
  std::vector<mpz_class> out(maximum + 1, 1);
  for (int i = 1; i <= maximum; ++i) out[i] = out[i - 1] * i;
  return out;
}

mpz_class hook_dimension(const std::vector<int>& partition,
                         const std::vector<mpz_class>& facts) {
  const int length = static_cast<int>(partition.size());
  int size = 0;
  for (int part : partition) size += part;
  mpz_class numerator = facts[size];
  mpz_class denominator = 1;
  for (int i = 0; i < length; ++i) {
    denominator *= facts[partition[i] + length - i - 1];
  }
  for (int i = 0; i < length; ++i) {
    for (int j = i + 1; j < length; ++j) {
      numerator *= partition[i] - partition[j] + j - i;
    }
  }
  if (!mpz_divisible_p(numerator.get_mpz_t(), denominator.get_mpz_t())) {
    fail("hook-dimension division was not exact");
  }
  return numerator / denominator;
}

void sum_partitions(int remaining, int max_part, int max_length,
                    std::vector<int>& partition,
                    const std::vector<mpz_class>& facts, mpz_class& total) {
  if (remaining == 0) {
    const mpz_class dimension = hook_dimension(partition, facts);
    total += dimension * dimension;
    return;
  }
  if (max_length == 0) return;
  const int first_max = std::min(max_part, remaining);
  for (int first = first_max; first >= 1; --first) {
    partition.push_back(first);
    sum_partitions(remaining - first, first, max_length - 1, partition, facts,
                   total);
    partition.pop_back();
  }
}

std::vector<mpz_class> trace_moments(int rank, int maximum) {
  const auto facts = factorials(maximum);
  std::vector<mpz_class> out(maximum + 1);
  out[0] = 1;
  std::vector<int> partition;
  partition.reserve(rank);
  for (int degree = 1; degree <= maximum; ++degree) {
    sum_partitions(degree, degree, rank, partition, facts, out[degree]);
  }
  return out;
}

std::vector<mpz_class> adjoint_moments(const std::vector<mpz_class>& trace) {
  std::vector<mpz_class> out(trace.size());
  for (unsigned r = 0; r < trace.size(); ++r) {
    mpz_class value = 0;
    for (unsigned a = 0; a <= r; ++a) {
      mpz_class term = binomial(r, a) * trace[a];
      if ((r - a) & 1U)
        value -= term;
      else
        value += term;
    }
    out[r] = value;
  }
  return out;
}

mpq_class shifted_moment(const std::vector<mpz_class>& moments, int degree,
                         const mpq_class& center) {
  mpq_class total = 0;
  mpq_class power = 1;
  const mpq_class minus_center = -center;
  // Accumulate in ascending t while maintaining (-center)^(degree-t).
  for (int t = 0; t <= degree; ++t) {
    const mpq_class center_power = rational_power(minus_center, degree - t);
    total += mpq_class(binomial(degree, t) * moments[t]) * center_power;
  }
  total.canonicalize();
  return total;
}

void verify_endpoint_constants() {
  const mpq_class delta_b_64(mpz_class(244495424),
                             mpz_class(25) * 1065343 * 1065343);
  const mpq_class delta_c_64(mpz_class(471756427),
                             mpz_class("178410035553750000"));
  if (delta_b_64 * delta_c_64 * rational_power(mpq_class(3, 2), 89) <= 72)
    fail("Proposition 64 endpoint inequality failed");

  const mpq_class eta_b_65(mpz_class("122923229137548665407"),
                           mpz_class(65536) * 1467426013 * 1467426013);
  const mpq_class eta_c_65(mpz_class("88255484124721"),
                           mpz_class("992717442773183102976"));

  // Reconstruct the N>=24 band witness from factorial trace moments.
  const mpz_class coefficients[] = {
      mpz_class("925503648"), mpz_class("348262482"),
      mpz_class("-142832255"), mpz_class("-40966535"),
      mpz_class("8723312"), mpz_class("880405"),
      mpz_class("-242416"), mpz_class("14686"), mpz_class("-274")};
  const auto stable = factorials(24);
  const mpq_class center(7, 2);
  mpq_class band_integral = 0;
  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 9; ++j) {
      const int degree = i + j;
      const mpq_class h = shifted_moment(stable, degree, center) -
                          shifted_moment(stable, degree + 2, center);
      band_integral += mpq_class(coefficients[i] * coefficients[j]) * h;
    }
  }
  const mpq_class expected_band_integral(
      mpz_class("122923229137548665407"), mpz_class(65536));
  if (band_integral != expected_band_integral)
    fail("Proposition 65 band-witness reconstruction failed");
  if (band_integral / mpq_class(mpz_class(1467426013) * 1467426013) != eta_b_65)
    fail("Proposition 65 eta_B identity failed");

  mpz_class five12;
  mpz_ui_pow_ui(five12.get_mpz_t(), 5, 12);
  const mpz_class pz_numerator = stable[12] - five12;
  mpq_class reconstructed_eta_c(pz_numerator * pz_numerator, stable[24]);
  reconstructed_eta_c.canonicalize();
  if (reconstructed_eta_c != eta_c_65)
    fail("Proposition 65 eta_C identity failed");
  if (eta_b_65 * eta_c_65 * rational_power(mpq_class(3, 2), 71) <= 200)
    fail("Proposition 65 endpoint inequality failed");
}

void verify_prop66() {
  for (int rank = 24; rank <= 27; ++rank) {
    const int first = rank == 24 ? 57 : rank == 25 ? 59 : rank == 26 ? 61 : 65;
    mpq_class minimum;
    bool have_minimum = false;
    for (int k = first; k < 70; k += 2) {
      const int K = k + 2;
      mpq_class tau(9 * binomial(K, rank + 1), factorials(rank + 1)[rank + 1]);
      const mpq_class bracket = mpq_class(7, 18) - tau * (4 * k + 10) -
                                tau * tau * (2 * k + 5);
      if (!have_minimum || bracket < minimum) {
        minimum = bracket;
        have_minimum = true;
      }
    }
    if (!have_minimum || minimum <= mpq_class(3, 10))
      fail("Proposition 66 row N=" + std::to_string(rank) + " failed");
  }
}

void verify_quadratic_transition_window() {
  // For Proposition \ref{lem:su-quadratic-window}, put theta=121/176.
  // The two error envelopes are decreasing once N>=25.  After shifting
  // N=s+25, coefficientwise positivity of denominator-theta*numerator
  // proves both ratio comparisons for every integer s>=0.
  const mpq_class theta(121, 176);
  const Polynomial N{0, 1};
  const Polynomial N_plus_1{1, 1};
  const Polynomial four_N_squared_plus_110 =
      polynomial_add(polynomial_multiply(Polynomial{4}, polynomial_power(N, 2)),
                     Polynomial{110});
  const Polynomial four_N_plus_1_squared_plus_110 =
      polynomial_add(
          polynomial_multiply(Polynomial{4}, polynomial_power(N_plus_1, 2)),
          Polynomial{110});
  const Polynomial first = polynomial_add(
      polynomial_multiply(four_N_squared_plus_110,
                          Polynomial{theta.get_den()}),
      polynomial_multiply(four_N_plus_1_squared_plus_110,
                          Polynomial{theta.get_num()}),
      -1);
  require_positive_coefficients(polynomial_shift(first, 25),
                                "quadratic-window first monotonicity");

  const Polynomial two_N_squared_plus_55 =
      polynomial_add(polynomial_multiply(Polynomial{2}, polynomial_power(N, 2)),
                     Polynomial{55});
  const Polynomial two_N_plus_1_squared_plus_55 =
      polynomial_add(
          polynomial_multiply(Polynomial{2}, polynomial_power(N_plus_1, 2)),
          Polynomial{55});
  const mpq_class theta_squared = theta * theta;
  const Polynomial second = polynomial_add(
      polynomial_multiply(two_N_squared_plus_55,
                          Polynomial{theta_squared.get_den()}),
      polynomial_multiply(two_N_plus_1_squared_plus_55,
                          Polynomial{theta_squared.get_num()}),
      -1);
  require_positive_coefficients(polynomial_shift(second, 25),
                                "quadratic-window second monotonicity");

  const mpq_class N25 = 25;
  const mpq_class error =
      9 * rational_power(theta, 26) *
          (4 * N25 * N25 / 11 + 10) +
      81 * rational_power(theta, 52) *
          (2 * N25 * N25 / 11 + 5);
  if (error >= mpq_class(129, 1000) || error >= mpq_class(7, 18))
    fail("quadratic-window N=25 error margin failed");
}

void verify_uniform_transition_strip(
    const std::vector<std::vector<mpz_class>>& trace) {
  // The only finite part of the strip theorem: all ranks below the uniform
  // analytic threshold, with exact hook-length moments.
  for (int rank = 6; rank <= 18; ++rank) {
    const auto moments = adjoint_moments(trace[rank]);
    for (int k = rank - 1; k <= rank + 33; ++k) {
      const mpz_class delta =
          moments[k + 2] * moments[k] - moments[k + 1] * moments[k + 1];
      if (delta <= 0) {
        fail("uniform transition strip failed at N=" + std::to_string(rank) +
             ", k=" + std::to_string(k));
      }
    }
  }

  // For N >= 19 and -1 <= d <= 33 the relative-error bound is worst at
  // d=33.  Its two error terms decrease with N.  Verify the exact N=19
  // margin and coefficientwise positivity of the two shifted polynomials
  // proving those monotonicities.
  constexpr int rank = 19;
  constexpr int offset = 33;
  constexpr int k = rank + offset;
  constexpr int K = k + 2;
  const auto facts = factorials(rank + 1);
  const mpq_class tau(9 * binomial(K, rank + 1), facts[rank + 1]);
  const mpq_class bracket = mpq_class(7, 18) - tau * (4 * k + 10) -
                            tau * tau * (2 * k + 5);
  if (bracket <= mpq_class(1, 10))
    fail("uniform transition strip base margin is not greater than 1/10");

  const Polynomial N_plus_2{2, 1};
  const Polynomial N_plus_36{36, 1};
  const Polynomial four_N_plus_142{142, 4};
  const Polynomial four_N_plus_146{146, 4};
  const Polynomial two_N_plus_71{71, 2};
  const Polynomial two_N_plus_73{73, 2};
  const Polynomial first = polynomial_add(
      polynomial_multiply(polynomial_power(N_plus_2, 2), four_N_plus_142),
      polynomial_multiply(N_plus_36, four_N_plus_146), -1);
  const Polynomial second = polynomial_add(
      polynomial_multiply(polynomial_power(N_plus_2, 4), two_N_plus_71),
      polynomial_multiply(polynomial_power(N_plus_36, 2), two_N_plus_73), -1);
  require_positive_coefficients(polynomial_shift(first, 19),
                                "uniform strip first monotonicity polynomial");
  require_positive_coefficients(polynomial_shift(second, 19),
                                "uniform strip second monotonicity polynomial");
}

void verify_prop67(const std::vector<std::vector<mpz_class>>& trace) {
  const mpz_class coefficients[] = {mpz_class("959324176457"),
                                    mpz_class("276017837985"),
                                    mpz_class("-58083597219"),
                                    mpz_class("-11664458272"),
                                    mpz_class("1230335979")};
  const mpq_class center(43, 10), radius(7, 5);
  const mpq_class S(mpz_class("935204207737834"), 625);
  const mpq_class band_target(mpz_class("44525786465344542780431067563007"),
                              mpz_class("10000000000"));
  const mpq_class delta_b = band_target / (mpq_class(49, 25) * S * S);
  const mpq_class expected_delta_b(
      mpz_class("44525786465344542780431067563007"),
      mpz_class("43884276324717505323728973379833856"));
  if (delta_b != expected_delta_b) fail("Proposition 67 delta_B identity failed");

  std::pair<mpq_class, int> minimum;
  bool have_minimum = false;
  for (int rank = 6; rank <= 9; ++rank) {
    mpq_class value = 0;
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) {
        const int n = i + j;
        const mpq_class h = radius * radius * shifted_moment(trace[rank], n, center) -
                            shifted_moment(trace[rank], n + 2, center);
        value += mpq_class(coefficients[i] * coefficients[j]) * h;
      }
    }
    const mpq_class excess = value - band_target;
    if (excess < 0) fail("Proposition 67 band positivity failed");
    if (!have_minimum || excess < minimum.first) {
      minimum = {excess, rank};
      have_minimum = true;
    }
  }
  if (minimum.first != mpq_class(mpz_class("1513726621221888441")) ||
      minimum.second != 9)
    fail("Proposition 67 band minimum mismatch");

  const mpq_class delta_c(mpz_class(471756427),
                          mpz_class("178410035553750000"));
  const auto stable_facts = factorials(30);
  mpz_class six15;
  mpz_ui_pow_ui(six15.get_mpz_t(), 6, 15);
  const mpz_class stable_pz_numerator = stable_facts[15] - six15;
  mpq_class reconstructed_delta_c(stable_pz_numerator * stable_pz_numerator,
                                  stable_facts[30]);
  reconstructed_delta_c.canonicalize();
  if (reconstructed_delta_c != delta_c)
    fail("Proposition 67 stable Paley-Zygmund identity failed");
  std::pair<mpq_class, int> tail_minimum;
  have_minimum = false;
  for (int rank = 6; rank <= 14; ++rank) {
    if (trace[rank][15] <= six15)
      fail("Proposition 67 Paley-Zygmund threshold is outside its domain");
    const mpq_class ratio = mpq_class(1) - mpq_class(six15, trace[rank][15]);
    const mpq_class tail = ratio * ratio * mpq_class(trace[rank][15] * trace[rank][15],
                                                     trace[rank][30]);
    const mpq_class excess = tail - delta_c;
    if (excess < 0) fail("Proposition 67 tail positivity failed");
    if (!have_minimum || excess < tail_minimum.first) {
      tail_minimum = {excess, rank};
      have_minimum = true;
    }
  }
  const mpq_class expected_tail(
      mpz_class("2602599137945210797524667371658649"),
      mpz_class("56336839730247629366109305276233265935886250000"));
  if (tail_minimum.first != expected_tail || tail_minimum.second != 14)
    fail("Proposition 67 tail minimum mismatch");
  if (delta_b * delta_c * rational_power(mpq_class(19, 10), 52) *
          rational_power(mpq_class(3, 10), 2) <=
      72)
    fail("Proposition 67 endpoint inequality failed");
}

void verify_prop68(const std::vector<std::vector<mpz_class>>& trace) {
  for (const auto& expected : kExpected) {
    const auto moments = adjoint_moments(trace[expected.rank]);
    int first = expected.rank + 33;
    if ((first & 1) == 0) ++first;
    bool have_minimum = false;
    mpz_class minimum;
    int minimum_k = -1;
    for (int k = first; k < 53; k += 2) {
      const mpz_class delta = moments[k + 2] * moments[k] - moments[k + 1] * moments[k + 1];
      if (!have_minimum || delta < minimum) {
        minimum = delta;
        minimum_k = k;
        have_minimum = true;
      }
    }
    if (!have_minimum || minimum_k != expected.moment || minimum != mpz_class(expected.delta))
      fail("Proposition 68 row N=" + std::to_string(expected.rank) + " mismatch");
  }
}

}  // namespace

int main() {
  try {
    verify_endpoint_constants();
    verify_prop66();
    verify_quadratic_transition_window();

    std::vector<std::vector<mpz_class>> trace(19);
#pragma omp parallel for schedule(dynamic, 1)
    for (int rank = 6; rank <= 18; ++rank) {
      trace[rank] = trace_moments(rank, 53);
    }
    verify_uniform_transition_strip(trace);
    verify_prop67(trace);
    verify_prop68(trace);

    for (int rank = 6; rank <= 18; ++rank) {
      std::cout << "SU(N) row N=" << rank << " exact moments through 53: ok\n";
    }
    std::cout << "Uniform transition strip N-1 <= k <= N+33: ok\n";
    std::cout << "Quadratic transition window N>=25: ok\n";
    std::cout << "All SU(N) repair certificates verified with exact GMP arithmetic.\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "SU(N) GMP certificate failure: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
