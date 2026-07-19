#include <gmpxx.h>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr int kMaxM = 218;
constexpr int kMaxProfile = 2 * kMaxM;

mpz_class binom_ui(unsigned long n, unsigned long k) {
  mpz_class out;
  mpz_bin_uiui(out.get_mpz_t(), n, k);
  return out;
}

std::string prefix_digits(const mpz_class& x, std::size_t n = 24) {
  std::string s = x.get_str();
  if (s.size() <= n) return s;
  return s.substr(0, n) + "...";
}

}  // namespace

int main() {
  std::vector<mpz_class> fact(kMaxProfile + 1);
  fact[0] = 1;
  for (int i = 1; i <= kMaxProfile; ++i) fact[i] = fact[i - 1] * i;

  std::vector<mpz_class> C(kMaxProfile + 1);
  C[0] = 1;
  if (kMaxProfile >= 1) C[1] = 0;
  for (int c = 2; c <= kMaxProfile; ++c) {
    mpz_class total = 0;
    for (int s = 2; s <= c; ++s) {
      mpz_class cyc = 1;
      if (s != 2) cyc = fact[s - 1] / 2;
      total += binom_ui(c - 1, s - 1) * cyc * C[c - s];
    }
    C[c] = total;
  }

  std::vector<std::vector<mpz_class>> N(
      kMaxProfile + 1, std::vector<mpz_class>(kMaxProfile + 1));
  for (int c = 0; c <= kMaxProfile; ++c) N[0][c] = C[c];
  for (int b = 1; b <= kMaxProfile; ++b) {
    if (b % 2 == 1) continue;
    for (int c = 0; c <= kMaxProfile; ++c) {
      mpz_class total = 0;
      for (int k = 0; k <= c; ++k) {
        total += binom_ui(c, k) * fact[k] * N[b - 2][c - k];
      }
      N[b][c] = total * (b - 1);
    }
  }

  int failures = 0;
  int checked_m = 0;
  int global_min_b = std::numeric_limits<int>::max();
  int global_min_active = std::numeric_limits<int>::max();
  int global_min_D = std::numeric_limits<int>::max();
  int global_max_closed_active = 0;
  int active20_failures = 0;
  mpz_class active20_max = 0;
  int active20_arg_b = -1;
  int active20_arg_c = -1;
  mpz_class first_dense_count = N[0][21];
  mpz_class m29_budget = mpz_class(1) << 58;
  std::vector<int> min_active_by_m(kMaxM + 1, -1);
  std::vector<int> min_D_by_m(kMaxM + 1, -1);

  std::cout << "PROFILE_GMP max_m=" << kMaxM
            << " max_profile=" << kMaxProfile << "\n";

  for (int m = 29; m <= kMaxM; ++m) {
    ++checked_m;
    mpz_class budget = mpz_class(1) << (2 * m);
    int bad_count = 0;
    int min_b = std::numeric_limits<int>::max();
    int min_active = std::numeric_limits<int>::max();
    int min_D = std::numeric_limits<int>::max();
    int max_closed_active = 0;
    mpz_class max_ratio_num = 0;
    int max_b = -1;
    int max_c = -1;

    for (int b = 0; b <= 2 * m; ++b) {
      for (int c = 0; c + b <= 2 * m; ++c) {
        const mpz_class& count = N[b][c];
        const int active = b + c;
        const int D = b + 2 * c;
        if (active <= 20 && count > active20_max) {
          active20_max = count;
          active20_arg_b = b;
          active20_arg_c = c;
        }
        if (count > budget) {
          ++bad_count;
          if (active <= 20) ++active20_failures;
          min_b = std::min(min_b, b);
          min_active = std::min(min_active, active);
          min_D = std::min(min_D, D);
          if (count > max_ratio_num) {
            max_ratio_num = count;
            max_b = b;
            max_c = c;
          }
        } else {
          max_closed_active = std::max(max_closed_active, active);
        }
      }
    }

    if (bad_count == 0) {
      failures++;
      std::cout << "m=" << m << " ERROR no_bad_profiles\n";
      continue;
    }

    global_min_b = std::min(global_min_b, min_b);
    global_min_active = std::min(global_min_active, min_active);
    global_min_D = std::min(global_min_D, min_D);
    global_max_closed_active = std::max(global_max_closed_active, max_closed_active);
    min_active_by_m[m] = min_active;
    min_D_by_m[m] = min_D;

    if (m <= 40 || m % 20 == 0 || m == kMaxM) {
      std::cout << "m=" << m << " bad_profiles=" << bad_count
                << " min_b=" << min_b
                << " min_active=" << min_active
                << " min_D=" << min_D
                << " max_closed_active=" << max_closed_active
                << " max_profile=(" << max_b << "," << max_c << ")"
                << " max_count_prefix=" << prefix_digits(max_ratio_num)
                << " budget_bits=" << 2 * m << "\n";
    }
  }

  if (active20_failures != 0) ++failures;
  if (active20_max > m29_budget) ++failures;
  if (!(first_dense_count > m29_budget)) ++failures;
  std::cout << "ACTIVE20 max_profile=(" << active20_arg_b << ","
            << active20_arg_c << ")"
            << " max_count=" << active20_max.get_str()
            << " budget_2^58=" << m29_budget.get_str()
            << " active20_failures=" << active20_failures << "\n";
  std::cout << "FIRST_DENSE profile=(0,21)"
            << " count=" << first_dense_count.get_str()
            << " exceeds_2^58=" << (first_dense_count > m29_budget ? "yes" : "no")
            << "\n";
  int interval_start = 29;
  while (interval_start <= kMaxM) {
    const int value = min_active_by_m[interval_start];
    int interval_end = interval_start;
    while (interval_end + 1 <= kMaxM &&
           min_active_by_m[interval_end + 1] == value) {
      ++interval_end;
    }
    std::cout << "ACTIVE_FRONTIER m=" << interval_start;
    if (interval_end != interval_start) std::cout << ".." << interval_end;
    std::cout << " min_dense_active=" << value
              << " closes_active_le=" << (value - 1) << "\n";
    interval_start = interval_end + 1;
  }
  interval_start = 29;
  while (interval_start <= kMaxM) {
    const int value = min_D_by_m[interval_start];
    int interval_end = interval_start;
    while (interval_end + 1 <= kMaxM &&
           min_D_by_m[interval_end + 1] == value) {
      ++interval_end;
    }
    std::cout << "DEGREE_FRONTIER m=" << interval_start;
    if (interval_end != interval_start) std::cout << ".." << interval_end;
    std::cout << " min_dense_D=" << value
              << " closes_D_le=" << (value - 1) << "\n";
    interval_start = interval_end + 1;
  }
  std::cout << "SUMMARY checked_m=" << checked_m
            << " failures=" << failures
            << " global_min_bad_b=" << global_min_b
            << " global_min_bad_active=" << global_min_active
            << " global_min_bad_D=" << global_min_D
            << " global_max_closed_active=" << global_max_closed_active
            << "\n";

  return failures == 0 ? 0 : 1;
}
