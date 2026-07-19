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

mpz_class compute_bad_count(const Case& c) {
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

  std::cout << "CASE q=" << c.q
            << " j=" << c.j
            << " width_wall=" << width_wall
            << " terminal_shapes=" << current.size() << "\n";

  for (int step = c.j; step >= 1; --step) {
    current = reverse_one_step(current);
    if (step % 5 == 0 || step <= 2 || step == c.j) {
      std::cout << "REVERSE_STEP q=" << c.q
                << " after_removing_to=" << (step - 1)
                << " states=" << current.size() << "\n";
    }
  }

  return current.count("") ? current.at("") : 0;
}

}  // namespace

int main() {
  int failures = 0;
  const std::vector<Case> cases = {
      {15, 38},
      {16, 40},
      {17, 42},
  };
  const std::vector<mpz_class> s = stable_moments(42);

  std::cout << "B_EIGHTH_FRONTIER_GMP cases=" << cases.size()
            << " omp_max_threads=" << omp_get_max_threads() << "\n";

  for (const Case& c : cases) {
    const mpz_class bad_count = compute_bad_count(c);
    const mpz_class stable = s[c.j];
    const mpz_class margin = stable - 2 * bad_count;
    if (margin < 0) {
      ++failures;
      std::cout << "FAIL half_stable q=" << c.q
                << " j=" << c.j
                << " bad_count=" << bad_count
                << " stable=" << stable << "\n";
    }
    std::cout << "BAD_COUNT q=" << c.q << " j=" << c.j
              << " value=" << bad_count << "\n";
    std::cout << "STABLE q=" << c.q << " j=" << c.j
              << " value=" << stable << "\n";
    std::cout << "MARGIN q=" << c.q << " j=" << c.j
              << " stable_minus_2bad=" << margin << "\n";
  }

  std::cout << "SUMMARY failures=" << failures << "\n";
  std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << "\n";
  return failures == 0 ? 0 : 1;
}
