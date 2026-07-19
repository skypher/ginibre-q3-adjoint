#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <omp.h>

namespace {

constexpr int H = 22;
constexpr int LOWER_BOXES = 2 * H;
constexpr int ONE_LOWER_CHOICES = 2 * LOWER_BOXES;
constexpr int TWO_LOWER_CHOICES = LOWER_BOXES * (LOWER_BOXES + 1) / 2;
constexpr int QUOTIENT_FLOOR = 31 - H;
constexpr int B_ABSORPTION_Q = 38;
constexpr int C_ABSORPTION_Q = 56;

struct Case {
  int q;
  int j;
};

std::string encode(const std::vector<int>& p) {
  std::string out;
  out.reserve(p.size());
  for (int x : p) out.push_back(static_cast<char>(x));
  return out;
}

std::vector<int> decode(const std::string& key) {
  std::vector<int> out;
  out.reserve(key.size());
  for (unsigned char c : key) out.push_back(static_cast<int>(c));
  return out;
}

void trim(std::vector<int>& p) {
  while (!p.empty() && p.back() == 0) p.pop_back();
}

bool is_partition(const std::vector<int>& p) {
  for (int x : p) {
    if (x < 0) return false;
  }
  for (std::size_t i = 1; i < p.size(); ++i) {
    if (p[i - 1] < p[i]) return false;
  }
  return true;
}

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

mpz_class pow2(int exponent) {
  return pow_ui(2, exponent);
}

void generate_mu_partitions(
    int remaining,
    int max_part,
    int width_wall,
    std::vector<int>& prefix,
    std::vector<std::vector<int>>& out
) {
  if (remaining == 0) {
    if (!prefix.empty() && prefix.front() >= width_wall) out.push_back(prefix);
    return;
  }
  for (int part = std::min(remaining, max_part); part >= 1; --part) {
    prefix.push_back(part);
    generate_mu_partitions(remaining - part, part, width_wall, prefix, out);
    prefix.pop_back();
  }
}

std::set<std::string> previous_horizontal_two_strip_shapes(const std::string& key) {
  const std::vector<int> mu = decode(key);
  std::set<std::string> out;
  const int rows = static_cast<int>(mu.size());

  for (int r = 0; r < rows; ++r) {
    if (mu[r] < 2) continue;
    std::vector<int> nu = mu;
    nu[r] -= 2;
    trim(nu);
    if (is_partition(nu)) out.insert(encode(nu));
  }

  for (int r = 0; r < rows; ++r) {
    if (mu[r] < 1) continue;
    for (int s = r + 1; s < rows; ++s) {
      if (mu[s] < 1) continue;
      if (mu[r] == mu[s]) continue;
      std::vector<int> nu = mu;
      --nu[r];
      --nu[s];
      trim(nu);
      if (is_partition(nu)) out.insert(encode(nu));
    }
  }

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

std::map<std::string, mpz_class> reverse_one_step(
    const std::map<std::string, mpz_class>& current
) {
  std::vector<std::pair<std::string, mpz_class>> items;
  items.reserve(current.size());
  for (const auto& kv : current) items.push_back(kv);

  const int threads = std::max(1, omp_get_max_threads());
  std::vector<std::map<std::string, mpz_class>> local(threads);

#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
#pragma omp for schedule(dynamic, 64)
    for (std::size_t i = 0; i < items.size(); ++i) {
      const auto& [shape, count] = items[i];
      const std::set<std::string> prevs = previous_horizontal_two_strip_shapes(shape);
      for (const auto& prev : prevs) local[tid][prev] += count;
    }
  }

  std::map<std::string, mpz_class> next;
  for (const auto& bucket : local) {
    for (const auto& [shape, count] : bucket) next[shape] += count;
  }
  return next;
}

mpz_class compute_b_width_bad_count(const Case& c) {
  const int width_wall = 2 * c.q;
  std::vector<std::vector<int>> mu_partitions;
  std::vector<int> prefix;
  generate_mu_partitions(c.j, c.j, width_wall, prefix, mu_partitions);

  std::map<std::string, mpz_class> current;
  for (const auto& mu : mu_partitions) {
    std::vector<int> lambda;
    lambda.reserve(2 * mu.size());
    for (int part : mu) {
      lambda.push_back(part);
      lambda.push_back(part);
    }
    current[encode(lambda)] += 1;
  }

  std::cout << "B_CASE q=" << c.q
            << " j=" << c.j
            << " width_wall=" << width_wall
            << " terminal_shapes=" << current.size() << "\n";

  for (int step = c.j; step >= 1; --step) {
    current = reverse_one_step(current);
    if (step % 5 == 0 || step <= 2 || step == c.j) {
      std::cout << "B_REVERSE_STEP q=" << c.q
                << " after_removing_to=" << (step - 1)
                << " states=" << current.size() << "\n";
    }
  }

  return current.count("") ? current.at("") : 0;
}

mpz_class c_twentysecond_pause_polynomial(int q) {
  const int j = 2 * q + H;
  mpz_class out = 0;
  for (int pauses = 0; pauses <= H; pauses += 2) {
    out += pow2(H + pauses) * binom(j, pauses);
  }
  return out;
}

mpz_class odd_double_factorial(int odd_n) {
  mpz_class out = 1;
  for (int k = 1; k <= odd_n; k += 2) out *= k;
  return out;
}

mpz_class c_twentysecond_fpf_pause_bound(int q) {
  const int j = 2 * q + H;
  mpz_class out = 0;
  for (int pauses = 0; pauses <= H; pauses += 2) {
    out += binom(j, pauses) * odd_double_factorial(j + pauses - 1);
  }
  return out;
}

mpq_class b_twentysecond_lower_box_rational() {
  mpq_class out(0);
  for (int a = 0; a <= LOWER_BOXES; ++a) {
    for (int c = 0; c <= LOWER_BOXES; ++c) {
      if (a + 2 * c > LOWER_BOXES) continue;
      mpq_class term(factorial(LOWER_BOXES));
      term /= factorial(a);
      term /= factorial(c);
      term /= pow_ui(QUOTIENT_FLOOR, LOWER_BOXES - a - c);
      term *= pow_ui(3, LOWER_BOXES - a - c);
      term *= pow_ui(TWO_LOWER_CHOICES, c);
      term /= pow_ui(ONE_LOWER_CHOICES, LOWER_BOXES - a);
      out += term;
    }
  }
  return out;
}

}  // namespace

int main() {
  std::cout << std::unitbuf;
  int failures = 0;
  const std::vector<Case> b_cases = {
      {15, 52},
      {16, 54},
      {17, 56},
      {18, 58},
      {19, 60},
      {20, 62},
      {21, 64},
      {22, 66},
      {23, 68},
      {24, 70},
      {25, 72},
      {26, 74},
      {27, 76},
      {28, 78},
      {29, 80},
      {30, 82},
      {31, 84},
      {32, 86},
      {33, 88},
      {34, 90},
      {35, 92},
      {36, 94},
      {37, 96},
  };
  const std::vector<int> c_fpf_q_cases = {21, 22, 23, 24, 25, 26, 27};
  const std::vector<int> c_q_cases = {28, 29, 30, 31, 32, 33, 34, 35,
                                      36, 37, 38, 39, 40, 41, 42, 43, 44,
                                      45, 46, 47, 48, 49, 50, 51, 52, 53,
                                      54, 55};
  const std::vector<mpz_class> s = stable_moments(280);

  std::cout << "BC_TWENTYSECOND_FRONTIER_GMP b_cases=" << b_cases.size()
            << " c_fpf_cases=" << c_fpf_q_cases.size()
            << " c_arithmetic_cases=" << c_q_cases.size()
            << " omp_max_threads=" << omp_get_max_threads() << "\n";

  const mpq_class lower_box_ratio = b_twentysecond_lower_box_rational();
  const mpz_class lower_box_num = lower_box_ratio.get_num();
  const mpz_class lower_box_den = lower_box_ratio.get_den();
  const mpz_class lower_box_ceil =
      (lower_box_num + lower_box_den - 1) / lower_box_den;
  std::cout << "B_LOWER_BOX_RATIONAL numerator=" << lower_box_num
            << " denominator=" << lower_box_den
            << " ceil=" << lower_box_ceil << "\n";

  {
    const int q = B_ABSORPTION_Q;
    const mpz_class left = 4 * lower_box_ceil * binom(2 * q + H, LOWER_BOXES)
        * pow_ui(ONE_LOWER_CHOICES, LOWER_BOXES) * pow_ui(3, 2 * q - H);
    const mpz_class right = factorial(2 * q + H - 1);
    if (left > right) ++failures;
    std::cout << "B_ABSORPTION_BASE q=" << q
              << " left=" << left
              << " right=" << right
              << " margin=" << (right - left) << "\n";

    const mpz_class ratio_left =
        9 * mpz_class(2 * q + H + 2) * (2 * q + H + 1);
    const mpz_class ratio_right =
        30 * mpz_class(2 * q - H + 2) * (2 * q - H + 1);
    if (ratio_left >= ratio_right) ++failures;
    std::cout << "B_ABSORPTION_RATIO_BASE q=" << q
              << " left=" << ratio_left
              << " right=" << ratio_right
              << " margin=" << (ratio_right - ratio_left) << "\n";
  }

  {
    const int q = C_ABSORPTION_Q;
    const mpz_class p = c_twentysecond_pause_polynomial(q);
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

  for (const Case& c : b_cases) {
    const mpz_class bad_count = compute_b_width_bad_count(c);
    const mpz_class stable = s[c.j];
    const mpz_class margin = stable - 2 * bad_count;
    if (margin < 0) {
      ++failures;
      std::cout << "B_FAIL half_stable q=" << c.q
                << " j=" << c.j
                << " bad_count=" << bad_count
                << " stable=" << stable << "\n";
    }
    std::cout << "B_BAD_COUNT q=" << c.q << " j=" << c.j
              << " value=" << bad_count << "\n";
    std::cout << "B_STABLE q=" << c.q << " j=" << c.j
              << " value=" << stable << "\n";
    std::cout << "B_MARGIN q=" << c.q << " j=" << c.j
              << " stable_minus_2bad=" << margin << "\n";
  }

  for (int q : c_fpf_q_cases) {
    const int j = 2 * q + H;
    const mpz_class bound = c_twentysecond_fpf_pause_bound(q);
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
    const mpz_class p = c_twentysecond_pause_polynomial(q);
    const mpz_class bound = pow2(2 * q + 1) * p * s[q + H];
    const mpz_class margin = s[j] - bound;
    if (margin < 0) {
      ++failures;
      std::cout << "C_FAIL twentysecond_pause q=" << q
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
