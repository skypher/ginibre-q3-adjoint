#include <iostream>
#include <stdexcept>
#include <vector>

#include <gmpxx.h>

namespace {

constexpr int H = 29;
constexpr int C_ABSORPTION_Q = 74;

mpz_class binom(int n, int k) {
  mpz_class out;
  if (k < 0 || k > n) return out;
  mpz_bin_uiui(out.get_mpz_t(), static_cast<unsigned long>(n),
               static_cast<unsigned long>(k));
  return out;
}

mpz_class pow_ui(int base, int exponent) {
  mpz_class out;
  mpz_ui_pow_ui(out.get_mpz_t(), static_cast<unsigned long>(base),
                static_cast<unsigned long>(exponent));
  return out;
}

mpz_class pow2(int exponent) {
  return pow_ui(2, exponent);
}

std::vector<mpz_class> stable_moments(int max_n) {
  std::vector<mpz_class> s(max_n + 1);
  s[0] = 1;
  if (max_n >= 1) s[1] = 0;
  for (int n = 1; n < max_n; ++n) {
    mpz_class next = n * s[n] + n * s[n - 1];
    if (n >= 2) next -= (mpz_class(n) * (n - 1) / 2) * s[n - 2];
    if (next < 0) throw std::runtime_error("stable recurrence became negative");
    s[n + 1] = next;
  }
  return s;
}

mpz_class c_twentyninth_pause_polynomial(int q) {
  const int j = 2 * q + H;
  mpz_class out = 0;
  for (int pauses = H % 2; pauses <= H; pauses += 2) {
    out += pow2(H + pauses) * binom(j, pauses);
  }
  return out;
}

mpz_class odd_double_factorial(int odd_n) {
  mpz_class out = 1;
  for (int k = 1; k <= odd_n; k += 2) out *= k;
  return out;
}

mpz_class c_twentyninth_fpf_pause_bound(int q) {
  const int j = 2 * q + H;
  mpz_class out = 0;
  for (int pauses = H % 2; pauses <= H; pauses += 2) {
    out += binom(j, pauses) * odd_double_factorial(j + pauses - 1);
  }
  return out;
}

}  // namespace

int main() {
  std::cout << std::unitbuf;
  int failures = 0;

  std::vector<int> c_fpf_q_cases;
  for (int q = 21; q <= H + 5; ++q) c_fpf_q_cases.push_back(q);
  std::vector<int> c_q_cases;
  for (int q = H + 6; q < C_ABSORPTION_Q; ++q) c_q_cases.push_back(q);
  const std::vector<mpz_class> s = stable_moments(300);

  std::cout << "BC_TWENTYNINTH_C_FRONTIER_GMP"
            << " c_fpf_cases=" << c_fpf_q_cases.size()
            << " c_arithmetic_cases=" << c_q_cases.size()
            << "\n";

  {
    const int q = C_ABSORPTION_Q;
    const mpz_class p = c_twentyninth_pause_polynomial(q);
    const mpz_class central = binom(2 * q + H, q);
    if (p > central) ++failures;
    std::cout << "C_PAUSE_ABSORPTION_BASE q=" << q
              << " P=" << p
              << " central=" << central
              << " margin=" << (central - p) << "\n";

    const mpz_class ratio_left = mpz_class(2 * q + H + 2) * (2 * q + H + 1);
    const mpz_class ratio_right = 2 * mpz_class(2 * q + 2) * (2 * q + 1);
    if (ratio_left >= ratio_right) ++failures;
    std::cout << "C_PAUSE_RATIO_BASE q=" << q
              << " left=" << ratio_left
              << " right=" << ratio_right
              << " margin=" << (ratio_right - ratio_left) << "\n";

    const mpz_class central_ratio_left =
        mpz_class(2 * q + H + 2) * (2 * q + H + 1);
    const mpz_class central_ratio_right =
        3 * mpz_class(q + 1) * (q + H + 1);
    if (central_ratio_left <= central_ratio_right) ++failures;
    std::cout << "C_CENTRAL_RATIO_BASE q=" << q
              << " left=" << central_ratio_left
              << " right=" << central_ratio_right
              << " margin=" << (central_ratio_left - central_ratio_right)
              << "\n";
  }

  for (int q : c_fpf_q_cases) {
    const int j = 2 * q + H;
    const mpz_class bound = c_twentyninth_fpf_pause_bound(q);
    const mpz_class stable = s[j];
    const mpz_class margin = stable - 2 * bound;
    if (margin < 0) {
      ++failures;
      std::cout << "C_FPF_FAIL half_stable q=" << q
                << " j=" << j
                << " bound=" << bound
                << " stable=" << stable << "\n";
    }
    std::cout << "C_FPF_BOUND q=" << q << " j=" << j
              << " value=" << bound << "\n";
    std::cout << "C_FPF_STABLE q=" << q << " j=" << j
              << " value=" << stable << "\n";
    std::cout << "C_FPF_MARGIN q=" << q << " j=" << j
              << " stable_minus_2bound=" << margin << "\n";
  }

  for (int q : c_q_cases) {
    const int j = 2 * q + H;
    const mpz_class p = c_twentyninth_pause_polynomial(q);
    const mpz_class bound = pow2(2 * q + 1) * p * s[q + H];
    const mpz_class margin = s[j] - bound;
    if (margin < 0) {
      ++failures;
      std::cout << "C_FAIL twentyninth_pause q=" << q
                << " j=" << j
                << " pause_polynomial=" << p
                << " stable=" << s[j] << "\n";
    }
    std::cout << "C_PAUSE_POLYNOMIAL q=" << q << " j=" << j
              << " value=" << p << "\n";
    std::cout << "C_BOUND q=" << q << " j=" << j
              << " value=" << bound << "\n";
    std::cout << "C_STABLE q=" << q << " j=" << j
              << " value=" << s[j] << "\n";
    std::cout << "C_MARGIN q=" << q << " j=" << j
              << " stable_minus_bound=" << margin << "\n";
  }

  std::cout << "SUMMARY failures=" << failures << "\n";
  std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << "\n";
  return failures == 0 ? 0 : 1;
}
