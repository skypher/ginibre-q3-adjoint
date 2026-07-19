#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gmpxx.h>
#include <omp.h>

namespace {

constexpr int kQ = 15;
constexpr int kJ = 37;
constexpr int kWidthWall = 2 * kQ;

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

void generate_mu_partitions(
    int remaining,
    int max_part,
    std::vector<int>& prefix,
    std::vector<std::vector<int>>& out
) {
  if (remaining == 0) {
    if (!prefix.empty() && prefix.front() >= kWidthWall) out.push_back(prefix);
    return;
  }
  for (int part = std::min(remaining, max_part); part >= 1; --part) {
    prefix.push_back(part);
    generate_mu_partitions(remaining - part, part, prefix, out);
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

}  // namespace

int main() {
  int failures = 0;
  std::cout << "B14_J37_BADSHAPE_GMP q=" << kQ
            << " j=" << kJ
            << " width_wall=" << kWidthWall
            << " omp_max_threads=" << omp_get_max_threads() << "\n";

  std::vector<std::vector<int>> mu_partitions;
  std::vector<int> prefix;
  generate_mu_partitions(kJ, kJ, prefix, mu_partitions);

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

  std::cout << "TERMINAL_SHAPES column_even_width_bad=" << current.size() << "\n";

  for (int step = kJ; step >= 1; --step) {
    std::map<std::string, mpz_class> next;
    for (const auto& [shape, count] : current) {
      const std::set<std::string> prevs = previous_horizontal_two_strip_shapes(shape);
      for (const auto& prev : prevs) next[prev] += count;
    }
    current.swap(next);
    if (step % 5 == 0 || step <= 2 || step == kJ) {
      std::cout << "REVERSE_STEP after_removing_to=" << (step - 1)
                << " states=" << current.size() << "\n";
    }
  }

  const mpz_class bad_count = current.count("") ? current.at("") : 0;
  const std::vector<mpz_class> s = stable_moments(kJ);
  const mpz_class stable = s[kJ];
  const mpz_class doubled_bad = 2 * bad_count;

  const mpz_class expected_bad("6170018716859356895648705996");
  const mpz_class expected_stable("991021752128621887393916677761923440691808");
  if (bad_count != expected_bad) {
    ++failures;
    std::cout << "FAIL bad_count expected=" << expected_bad
              << " actual=" << bad_count << "\n";
  }
  if (stable != expected_stable) {
    ++failures;
    std::cout << "FAIL stable_s37 expected=" << expected_stable
              << " actual=" << stable << "\n";
  }
  if (doubled_bad > stable) {
    ++failures;
    std::cout << "FAIL half_stable_bound 2*bad_count=" << doubled_bad
              << " stable=" << stable << "\n";
  }

  std::cout << "BAD_COUNT " << bad_count << "\n";
  std::cout << "STABLE_S37 " << stable << "\n";
  std::cout << "MARGIN stable_minus_2bad " << (stable - doubled_bad) << "\n";
  std::cout << "SUMMARY failures=" << failures << "\n";
  std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << "\n";
  return failures == 0 ? 0 : 1;
}
