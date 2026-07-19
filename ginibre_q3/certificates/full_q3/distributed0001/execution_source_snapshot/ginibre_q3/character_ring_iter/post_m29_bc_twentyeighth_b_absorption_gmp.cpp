#include <iostream>
#include <stdexcept>
#include <vector>

#include <gmpxx.h>

namespace {

constexpr int H = 28;
constexpr int LOWER_BOXES = 2 * H;
constexpr int ONE_LOWER_CHOICES = 2 * LOWER_BOXES;
constexpr int TWO_LOWER_CHOICES = LOWER_BOXES * (LOWER_BOXES + 1) / 2;

mpz_class binom(int n, int k) {
  mpz_class out;
  if (k < 0 || k > n) return out;
  mpz_bin_uiui(out.get_mpz_t(), static_cast<unsigned long>(n),
               static_cast<unsigned long>(k));
  return out;
}

mpz_class factorial(int n) {
  mpz_class out;
  mpz_fac_ui(out.get_mpz_t(), static_cast<unsigned long>(n));
  return out;
}

mpz_class pow_ui(int base, int exponent) {
  mpz_class out;
  mpz_ui_pow_ui(out.get_mpz_t(), static_cast<unsigned long>(base),
                static_cast<unsigned long>(exponent));
  return out;
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

mpq_class lower_box_rational_with_floor(int quotient_floor) {
  if (quotient_floor <= 0) {
    throw std::runtime_error("quotient floor must be positive");
  }

  mpq_class out(0);
  for (int a = 0; a <= LOWER_BOXES; ++a) {
    for (int c = 0; c <= LOWER_BOXES; ++c) {
      if (a + 2 * c > LOWER_BOXES) continue;
      mpq_class term(factorial(LOWER_BOXES));
      term /= factorial(a);
      term /= factorial(c);
      term /= pow_ui(quotient_floor, LOWER_BOXES - a - c);
      term *= pow_ui(3, LOWER_BOXES - a - c);
      term *= pow_ui(TWO_LOWER_CHOICES, c);
      term /= pow_ui(ONE_LOWER_CHOICES, LOWER_BOXES - a);
      out += term;
    }
  }
  return out;
}

mpq_class lower_box_rational_for_q(int q) {
  return lower_box_rational_with_floor(2 * q - H + 1);
}

void print_rational(const char* tag, const mpq_class& value) {
  std::cout << tag
            << " numerator=" << value.get_num()
            << " denominator=" << value.get_den() << "\n";
}

}  // namespace

int main() {
  std::cout << std::unitbuf;
  int failures = 0;

  std::cout << "BC_TWENTYEIGHTH_B_ABSORPTION_GMP"
            << " h=" << H
            << " lower_boxes=" << LOWER_BOXES
            << " one_lower_choices=" << ONE_LOWER_CHOICES
            << " two_lower_choices=" << TWO_LOWER_CHOICES << "\n";

  const std::vector<mpz_class> s = stable_moments(150);

  for (int q = 42; q <= 49; ++q) {
    const int n = 2 * q + H;
    const int quotient_floor = 2 * q - H + 1;
    const mpq_class ratio = lower_box_rational_for_q(q);
    const mpz_class left_num =
        2 * ratio.get_num() * binom(n, LOWER_BOXES)
        * pow_ui(ONE_LOWER_CHOICES, LOWER_BOXES) * pow_ui(3, 2 * q - H);
    const mpz_class right_num = s[n] * ratio.get_den();
    if (left_num > right_num) ++failures;
    std::cout << "B_VARIABLE_ABSORPTION q=" << q
              << " j=" << n
              << " quotient_floor=" << quotient_floor
              << " numerator=" << ratio.get_num()
              << " denominator=" << ratio.get_den()
              << " left_num=" << left_num
              << " right_num=" << right_num
              << " margin_num=" << (right_num - left_num) << "\n";
  }

  const mpq_class fixed_ratio = lower_box_rational_with_floor(31 - H);
  const mpz_class fixed_num = fixed_ratio.get_num();
  const mpz_class fixed_den = fixed_ratio.get_den();
  const mpz_class fixed_ceil = (fixed_num + fixed_den - 1) / fixed_den;
  print_rational("B_FIXED_LOWER_BOX_RATIONAL", fixed_ratio);
  std::cout << "B_FIXED_LOWER_BOX_CEIL value=" << fixed_ceil << "\n";

  {
    const int q = 50;
    const mpz_class left = 4 * fixed_ceil * binom(2 * q + H, LOWER_BOXES)
        * pow_ui(ONE_LOWER_CHOICES, LOWER_BOXES) * pow_ui(3, 2 * q - H);
    const mpz_class right = factorial(2 * q + H - 1);
    if (left > right) ++failures;
    std::cout << "B_FIXED_ABSORPTION_BASE q=" << q
              << " left=" << left
              << " right=" << right
              << " margin=" << (right - left) << "\n";

    const mpz_class ratio_left =
        9 * mpz_class(2 * q + H + 2) * (2 * q + H + 1);
    const mpz_class ratio_right =
        30 * mpz_class(2 * q - H + 2) * (2 * q - H + 1);
    if (ratio_left >= ratio_right) ++failures;
    std::cout << "B_FIXED_ABSORPTION_RATIO_BASE q=" << q
              << " left=" << ratio_left
              << " right=" << ratio_right
              << " margin=" << (ratio_right - ratio_left) << "\n";
  }

  std::cout << "SUMMARY failures=" << failures << "\n";
  std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << "\n";
  return failures == 0 ? 0 : 1;
}
