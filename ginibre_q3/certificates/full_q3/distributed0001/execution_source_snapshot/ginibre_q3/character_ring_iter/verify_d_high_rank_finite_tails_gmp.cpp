#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
using Polynomial = std::vector<mpq_class>;

struct Parameters {
  int rank, exact_through, dimension, kappa, c_value, exp_terms;
  int target_exponent, bridge_max_chain;
  mpq_class radius, quartic_s, central_t, quartic_r;
  std::string parity;
};

struct RowResult {
  bool ok = false;
  std::string error;
  int minimum_chain = -1;
  mpz_class minimum_margin;
  mpq_class tail_margin;
  double seconds = 0;
};

[[noreturn]] void fail(const std::string& message) { throw std::runtime_error(message); }

mpq_class parse_rational(const std::string& text) {
  const auto slash = text.find('/');
  if (slash == std::string::npos) return mpq_class(mpz_class(text));
  mpq_class out(mpz_class(text.substr(0, slash)), mpz_class(text.substr(slash + 1)));
  out.canonicalize();
  return out;
}

mpq_class qpow(mpq_class base, unsigned exponent) {
  mpq_class out = 1;
  while (exponent) {
    if (exponent & 1U) out *= base;
    exponent >>= 1U;
    if (exponent) base *= base;
  }
  out.canonicalize();
  return out;
}

mpz_class zpow(const mpz_class& base, unsigned exponent) {
  mpz_class out;
  mpz_pow_ui(out.get_mpz_t(), base.get_mpz_t(), exponent);
  return out;
}

mpz_class factorial(unsigned n) {
  mpz_class out;
  mpz_fac_ui(out.get_mpz_t(), n);
  return out;
}

std::vector<mpz_class> binomial_row(int n) {
  std::vector<mpz_class> row(n + 1);
  row[0] = 1;
  for (int k = 0; k < n; ++k) {
    row[k + 1] = row[k] * (n - k);
    if (!mpz_divisible_ui_p(row[k + 1].get_mpz_t(), k + 1))
      fail("inexact binomial recurrence");
    row[k + 1] /= k + 1;
  }
  return row;
}

std::vector<mpz_class> stable_moments(int maximum) {
  std::vector<mpz_class> moments(maximum + 1);
  moments[0] = 1;
  if (maximum >= 1) moments[1] = 0;
  for (int index = 1; index < maximum; ++index) {
    mpz_class value = index * moments[index] + index * moments[index - 1];
    if (index >= 2) value -= (mpz_class(index) * (index - 1) / 2) * moments[index - 2];
    moments[index + 1] = value;
  }
  return moments;
}

mpz_class q3_integer(const std::vector<mpz_class>& moments, int n,
                     const std::vector<mpz_class>& row) {
  mpz_class total = 0;
  for (int k = 0; k <= n; ++k)
    total += row[k] * (moments[k + 2] * moments[n - k] -
                       moments[k + 1] * moments[n - k + 1]);
  return 2 * total;
}

void add_linear(std::vector<mpz_class>& coefficients,
                const std::vector<mpz_class>& moments, int n, int scale,
                const std::vector<mpz_class>& row) {
  for (int k = 0; k <= n; ++k) {
    const mpz_class c = 2 * scale * row[k];
    int first = k + 2, second = n - k;
    coefficients[first] += c * moments[second];
    coefficients[second] += c * moments[first];
    first = k + 1;
    second = n - k + 1;
    coefficients[first] -= c * moments[second];
    coefficients[second] -= c * moments[first];
  }
}

mpz_class moment_bound(int index, const std::map<int, mpz_class>& deltas,
                       const std::vector<mpz_class>& stable) {
  const auto found = deltas.find(index);
  return found == deltas.end() ? stable[index] : abs(found->second);
}

mpz_class quadratic_negative_bound(int n, int scale,
                                   const std::vector<mpz_class>& row, int rank,
                                   const std::map<int, mpz_class>& deltas,
                                   const std::vector<mpz_class>& stable) {
  const int sum = n + 2;
  mpz_class total = 0;
  for (int first = 0; first <= sum / 2; ++first) {
    const int second = sum - first;
    mpz_class coefficient = 0;
    auto add = [&](int k, int sign, bool positive) {
      if (k < 0 || k > n) return;
      int left = positive ? k + 2 : k + 1;
      int right = positive ? n - k : n - k + 1;
      if (left > right) std::swap(left, right);
      if (left == first && right == second) coefficient += sign * 2 * scale * row[k];
    };
    add(first - 2, 1, true);
    if (second != first) add(second - 2, 1, true);
    add(first - 1, -1, false);
    if (second != first) add(second - 1, -1, false);
    if (coefficient < 0 && first >= rank && second >= rank)
      total += coefficient * moment_bound(first, deltas, stable) *
               moment_bound(second, deltas, stable);
  }
  return total;
}

void add_quadratic_stable_buckets(int n, int scale,
                                  const std::vector<mpz_class>& row,
                                  const std::vector<mpz_class>& stable,
                                  std::vector<mpz_class>& buckets) {
  const int sum = n + 2;
  for (int first = 0; first <= sum / 2; ++first) {
    const int second = sum - first;
    mpz_class coefficient = 0;
    auto add = [&](int k, int sign, bool positive) {
      if (k < 0 || k > n) return;
      int left = positive ? k + 2 : k + 1;
      int right = positive ? n - k : n - k + 1;
      if (left > right) std::swap(left, right);
      if (left == first && right == second) coefficient += sign * 2 * scale * row[k];
    };
    add(first - 2, 1, true);
    if (second != first) add(second - 2, 1, true);
    add(first - 1, -1, false);
    if (second != first) add(second - 1, -1, false);
    if (coefficient < 0) {
      mpz_class contribution = coefficient * stable[first];
      contribution *= stable[second];
      buckets[first] += contribution;
    }
  }
}

Polynomial poly_add(const Polynomial& a, const Polynomial& b) {
  Polynomial out(std::max(a.size(), b.size()));
  for (size_t i = 0; i < a.size(); ++i) out[i] += a[i];
  for (size_t i = 0; i < b.size(); ++i) out[i] += b[i];
  return out;
}

Polynomial poly_scale(Polynomial a, const mpq_class& factor) {
  for (auto& coefficient : a) coefficient *= factor;
  return a;
}

Polynomial poly_multiply(const Polynomial& a, const Polynomial& b) {
  Polynomial out(a.size() + b.size() - 1);
  for (size_t i = 0; i < a.size(); ++i)
    for (size_t j = 0; j < b.size(); ++j) out[i + j] += a[i] * b[j];
  return out;
}

Polynomial chebyshev_t(int degree) {
  if (degree == 0) return {mpq_class(1)};
  if (degree == 1) return {mpq_class(0), mpq_class(1)};
  Polynomial previous{mpq_class(1)}, current{mpq_class(0), mpq_class(1)};
  for (int index = 2; index <= degree; ++index) {
    Polynomial shifted(current.size() + 1);
    for (size_t j = 0; j < current.size(); ++j) shifted[j + 1] = 2 * current[j];
    Polynomial next = poly_add(shifted, poly_scale(previous, -1));
    previous = std::move(current);
    current = std::move(next);
  }
  return current;
}

Polynomial compose(const Polynomial& polynomial, const Polynomial& argument) {
  Polynomial out{mpq_class(0)}, power{mpq_class(1)};
  for (const auto& coefficient : polynomial) {
    out = poly_add(out, poly_scale(power, coefficient));
    power = poly_multiply(power, argument);
  }
  return out;
}

mpq_class evaluate(const Polynomial& polynomial, const mpq_class& value) {
  mpq_class out = 0;
  for (auto it = polynomial.rbegin(); it != polynomial.rend(); ++it) out = out * value + *it;
  return out;
}

mpq_class exp_upper_bound(const mpq_class& value, int terms) {
  if (value < 0 || terms < 0) fail("invalid exponential-series input");
  const mpz_class numerator = value.get_num();
  const mpz_class denominator_base = value.get_den();
  const mpz_class denominator = zpow(denominator_base, terms) * factorial(terms);
  mpz_class component = denominator;
  mpz_class total_numerator = component;
  for (int index = 0; index < terms; ++index) {
    component *= numerator;
    const mpz_class divisor = denominator_base * (index + 1);
    if (!mpz_divisible_p(component.get_mpz_t(), divisor.get_mpz_t()))
      fail("exponential common-denominator recurrence was not exact");
    mpz_divexact(component.get_mpz_t(), component.get_mpz_t(), divisor.get_mpz_t());
    total_numerator += component;
  }
  const mpq_class partial_sum(total_numerator, denominator);
  const mpq_class next_term(component * numerator,
                            denominator * denominator_base * (terms + 1));
  const mpq_class ratio = value / (terms + 2);
  if (ratio >= 1) fail("Taylor tail ratio is not below one");
  return partial_sum + next_term / (1 - ratio);
}

std::vector<Parameters> read_parameters(const std::string& path, int low, int high) {
  std::ifstream input(path);
  if (!input) fail("cannot open parameter table: " + path);
  std::vector<Parameters> out;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    Parameters p;
    std::string radius, quartic_s, central_t, quartic_r;
    std::istringstream row(line);
    if (!(row >> p.rank >> p.exact_through >> p.dimension >> p.kappa >> p.c_value >>
          radius >> p.exp_terms >> quartic_s >> p.parity >> central_t >> quartic_r >>
          p.target_exponent >> p.bridge_max_chain))
      fail("malformed parameter row: " + line);
    p.radius = parse_rational(radius);
    p.quartic_s = parse_rational(quartic_s);
    p.central_t = parse_rational(central_t);
    p.quartic_r = parse_rational(quartic_r);
    if (p.rank >= low && p.rank <= high) out.push_back(std::move(p));
  }
  if (out.empty()) fail("selected rank range has no parameter rows");
  for (size_t i = 1; i < out.size(); ++i)
    if (out[i].rank != out[i - 1].rank + 1) fail("parameter table has a rank gap");
  return out;
}

std::map<int, std::map<int, mpz_class>> read_delta_logs(
    const std::vector<std::string>& paths) {
  const std::regex pattern(R"(^\s*D_(\d+) Delta_(\d+) = (-?\d+))");
  std::map<int, std::map<int, mpz_class>> out;
  for (const auto& path : paths) {
    std::ifstream input(path);
    if (!input) fail("cannot open determinant-delta log: " + path);
    bool success = false;
    std::string line;
    while (std::getline(input, line)) {
      if (line.find("__EXIT_STATUS=0") != std::string::npos) success = true;
      std::smatch match;
      if (!std::regex_search(line, match, pattern)) continue;
      const int rank = std::stoi(match[1].str()), index = std::stoi(match[2].str());
      const mpz_class value(match[3].str());
      if (value < 0) fail("negative determinant delta in " + path);
      auto& slot = out[rank][index];
      if (slot != 0 && slot != value) fail("conflicting determinant deltas");
      slot = value;
    }
    if (!success) fail(path + " does not record __EXIT_STATUS=0");
  }
  return out;
}

void verify_stable_table(const std::string& path, const std::vector<mpz_class>& stable) {
  std::ifstream input(path);
  if (!input) fail("cannot open stable-moment table: " + path);
  std::vector<bool> seen(101);
  int index;
  std::string value;
  while (input >> index >> value) {
    if (index < 0 || index > 100) continue;
    if (stable[index] != mpz_class(value))
      fail("stable formula disagrees with table at index " + std::to_string(index));
    seen[index] = true;
  }
  if (std::find(seen.begin(), seen.end(), false) != seen.end())
    fail("stable-moment table is incomplete through index 100");
}

std::vector<std::pair<int, mpz_class>> precompute_stable_bridges(
    const std::vector<Parameters>& parameters, const std::vector<mpz_class>& stable) {
  std::vector<std::pair<int, mpz_class>> minima(parameters.size(), {-1, 0});
  int maximum_chain = 0;
  bool have_stable_rows = false;
  for (const auto& p : parameters) {
    if (p.exact_through < p.rank) {
      have_stable_rows = true;
      maximum_chain = std::max(maximum_chain, p.bridge_max_chain);
    }
  }
  if (!have_stable_rows) return minima;

  std::atomic<bool> failed{false};
  std::string failure_message;
#pragma omp parallel
  {
    std::vector<std::pair<int, mpz_class>> local(parameters.size(), {-1, 0});
#pragma omp for schedule(dynamic, 1)
    for (int chain = 31; chain <= maximum_chain; ++chain) {
      if (failed.load(std::memory_order_relaxed)) continue;
      try {
        const int n1 = 2 * chain + 3, n2 = 2 * chain + 1;
        const auto row1 = binomial_row(n1), row2 = binomial_row(n2);
        const mpz_class stable_diff = q3_integer(stable, n1, row1) -
                                      4 * q3_integer(stable, n2, row2);
        std::vector<mpz_class> linear(2 * chain + 6);
        add_linear(linear, stable, n1, 1, row1);
        add_linear(linear, stable, n2, -4, row2);
        std::vector<mpz_class> buckets(linear.size());
        for (size_t index = 0; index < linear.size(); ++index) {
          if (linear[index] < 0) {
            const mpz_class contribution = linear[index] * stable[index];
            buckets[index] += contribution;
          }
        }
        add_quadratic_stable_buckets(n1, 1, row1, stable, buckets);
        add_quadratic_stable_buckets(n2, -4, row2, stable, buckets);
        for (int index = static_cast<int>(buckets.size()) - 2; index >= 0; --index)
          buckets[index] += buckets[index + 1];

        for (size_t index = 0; index < parameters.size(); ++index) {
          const auto& p = parameters[index];
          if (p.exact_through >= p.rank || chain > p.bridge_max_chain) continue;
          const mpz_class margin =
              p.rank < static_cast<int>(buckets.size())
                  ? stable_diff + buckets[p.rank]
                  : stable_diff;
          if (margin <= 0)
            fail("D" + std::to_string(p.rank) + " stable bridge failed at m=" +
                 std::to_string(chain));
          if (local[index].first < 0 || margin < local[index].second)
            local[index] = {chain, margin};
        }
      } catch (const std::exception& error) {
        failed.store(true, std::memory_order_relaxed);
#pragma omp critical(d_bridge_failure)
        {
          if (failure_message.empty()) failure_message = error.what();
        }
      }
    }
#pragma omp critical(d_bridge_minima)
    {
      for (size_t index = 0; index < parameters.size(); ++index) {
        if (local[index].first < 0) continue;
        if (minima[index].first < 0 || local[index].second < minima[index].second)
          minima[index] = std::move(local[index]);
      }
    }
  }
  if (failed) fail(failure_message);
  return minima;
}

RowResult verify_row(const Parameters& p, const std::map<int, mpz_class>& supplied,
                     const std::vector<mpz_class>& stable,
                     const std::pair<int, mpz_class>* precomputed_bridge) {
  RowResult result;
  const auto start = Clock::now();
  try {
    std::map<int, mpz_class> deltas;
    for (const auto& [index, value] : supplied)
      if (index <= p.exact_through) deltas[index] = value;
    for (int index = p.rank; index <= p.exact_through; ++index)
      if (!deltas.count(index))
        fail("D" + std::to_string(p.rank) + " is missing Delta_" + std::to_string(index));

    bool have_minimum = false;
    if (precomputed_bridge != nullptr) {
      if (precomputed_bridge->first < 0) fail("missing precomputed stable bridge");
      result.minimum_chain = precomputed_bridge->first;
      result.minimum_margin = precomputed_bridge->second;
      have_minimum = true;
    }
    for (int chain = 31; precomputed_bridge == nullptr && chain <= p.bridge_max_chain;
         ++chain) {
      const int n1 = 2 * chain + 3, n2 = 2 * chain + 1;
      const auto row1 = binomial_row(n1), row2 = binomial_row(n2);
      const mpz_class stable_diff = q3_integer(stable, n1, row1) -
                                    4 * q3_integer(stable, n2, row2);
      std::vector<mpz_class> linear(2 * chain + 6);
      add_linear(linear, stable, n1, 1, row1);
      add_linear(linear, stable, n2, -4, row2);
      mpz_class linear_bound = 0;
      for (int index = p.rank; index < static_cast<int>(linear.size()); ++index) {
        if (linear[index] >= 0) continue;
        const auto found = deltas.find(index);
        linear_bound += linear[index] *
                        (found == deltas.end() ? stable[index] : found->second);
      }
      const mpz_class quadratic_bound =
          quadratic_negative_bound(n1, 1, row1, p.rank, deltas, stable) +
          quadratic_negative_bound(n2, -4, row2, p.rank, deltas, stable);
      const mpz_class margin = stable_diff + linear_bound + quadratic_bound;
      if (margin <= 0)
        fail("D" + std::to_string(p.rank) + " bridge failed at m=" + std::to_string(chain));
      if (!have_minimum || margin < result.minimum_margin) {
        result.minimum_margin = margin;
        result.minimum_chain = chain;
        have_minimum = true;
      }
    }

    const mpq_class x_floor = mpq_class(p.dimension) - p.radius;
    const mpq_class c_push = x_floor - p.central_t;
    if (c_push <= 0) fail("nonpositive pushforward cutoff");
    const mpq_class q_bound = 2 * p.radius / (39 * p.kappa);
    if (q_bound >= 1) fail("sine-series q bound is not below one");
    const mpq_class sine_exponent =
        p.radius / 12 + p.radius * p.radius / (1440 * p.kappa) +
        qpow(p.radius, 3) / (mpq_class(p.kappa) * p.kappa * 90720 * (1 - q_bound));
    const mpq_class sine_factor = 1 / exp_upper_bound(sine_exponent, p.exp_terms);

    const mpq_class quartic_moment =
        6 - 2 * (p.quartic_r + p.quartic_s) + p.quartic_r * p.quartic_r +
        p.quartic_s * p.quartic_s + 4 * p.quartic_r * p.quartic_s +
        p.quartic_r * p.quartic_r * p.quartic_s * p.quartic_s;
    const mpq_class quartic_denominator =
        qpow(p.central_t + p.quartic_r, 2) * qpow(p.central_t + p.quartic_s, 2);
    const mpq_class tail_weight_upper =
        qpow(x_floor + p.central_t, 2) * quartic_moment / quartic_denominator;
    const mpq_class pushforward_weight = x_floor * x_floor + 1 - tail_weight_upper;
    const mpq_class pushforward_derivative =
        2 * x_floor - 2 * quartic_moment * (x_floor + p.central_t) / quartic_denominator;
    const mpq_class pushforward_second = 2 - 2 * quartic_moment / quartic_denominator;
    if (pushforward_weight <= 0 || pushforward_derivative <= 0 || pushforward_second <= 0)
      fail("quartic pushforward positivity/monotonicity failed");

    const mpz_class weyl_order = zpow(mpz_class(2), p.rank - 1) * factorial(p.rank);
    mpz_class degree_factor = factorial(p.rank);
    for (int degree = 2; degree <= 2 * (p.rank - 1); degree += 2)
      degree_factor *= factorial(degree);
    mpq_class density_factor, radial_factor;
    if (p.parity == "even") {
      if (p.dimension & 1) fail("even density selected for odd dimension");
      density_factor = 1 / (qpow(mpq_class(44, 7), p.rank / 2) *
                            mpq_class(factorial(p.dimension / 2)));
      radial_factor = qpow(p.radius / (2 * p.kappa), p.dimension / 2);
    } else if (p.parity == "odd") {
      if (!(p.dimension & 1)) fail("odd density selected for even dimension");
      const int half = (p.dimension + 1) / 2, pi_power = (p.rank + 1) / 2;
      density_factor = mpq_class(7, 5) * zpow(mpz_class(4), half) * factorial(half) /
                       (qpow(mpq_class(44, 7), pi_power) * factorial(2 * half));
      radial_factor = qpow(p.radius / (2 * p.kappa), half);
    } else {
      fail("unknown density parity");
    }
    const mpq_class weighted_tail =
        pushforward_weight * sine_factor * degree_factor * density_factor * radial_factor /
        weyl_order;

    constexpr int chebyshev_degree = 8;
    std::vector<mpz_class> moments(stable.begin(), stable.begin() + 35);
    for (const auto& [index, value] : deltas)
      if (index < static_cast<int>(moments.size())) moments[index] += value;
    const mpq_class support = 2 * p.dimension, scale = 12000, cutoff = 80;
    const mpq_class q0 = cutoff * (cutoff + support);
    const Polynomial chebyshev = chebyshev_t(chebyshev_degree);
    const Polynomial argument{1, -2 * support / scale, 2 / scale};
    const Polynomial composed = compose(chebyshev, argument);
    const Polynomial majorant = poly_multiply(composed, composed);
    const mpq_class denominator = qpow(evaluate(chebyshev, 1 + 2 * q0 / scale), 2);
    mpq_class majorant_moment = 0;
    for (size_t index = 0; index < majorant.size(); ++index)
      majorant_moment += majorant[index] *
                         q3_integer(moments, index, binomial_row(index));
    const mpq_class negative_tail = majorant_moment / denominator;
    if (negative_tail <= 0) fail("Chebyshev negative-tail bound is not positive");
    const mpq_class cutoff_power = qpow(c_push, p.target_exponent);
    const mpq_class target =
        negative_tail * zpow(mpz_class(2 * p.c_value), p.target_exponent) / cutoff_power +
        2 * zpow(mpz_class(80), p.target_exponent) / cutoff_power;
    result.tail_margin = weighted_tail - target;
    if (result.tail_margin <= 0) fail("sharpened direct-tail margin is not positive");
    result.ok = true;
  } catch (const std::exception& error) {
    result.error = error.what();
  }
  result.seconds = std::chrono::duration<double>(Clock::now() - start).count();
  return result;
}

size_t bits(const mpz_class& value) {
  if (value == 0) return 0;
  const mpz_class magnitude = abs(value);
  return mpz_sizeinbase(magnitude.get_mpz_t(), 2);
}
}  // namespace

int main(int argc, char** argv) {
  try {
    std::string parameters_path, stable_path;
    std::vector<std::string> delta_paths;
    int low = 22, high = 97, threads = 0;
    bool print_exact = false;
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      auto value = [&]() -> std::string {
        if (i + 1 >= argc) fail(arg + " requires a value");
        return argv[++i];
      };
      if (arg == "--parameters") parameters_path = value();
      else if (arg == "--stable-moments") stable_path = value();
      else if (arg == "--delta-log") delta_paths.push_back(value());
      else if (arg == "--rank-lo") low = std::stoi(value());
      else if (arg == "--rank-hi") high = std::stoi(value());
      else if (arg == "--threads") threads = std::stoi(value());
      else if (arg == "--print-exact") print_exact = true;
      else fail("unknown argument: " + arg);
    }
    if (parameters_path.empty() || stable_path.empty())
      fail("required: --parameters FILE --stable-moments FILE");
    if (threads > 0) omp_set_num_threads(threads);
    const auto start = Clock::now();
    const auto parameters = read_parameters(parameters_path, low, high);
    const auto deltas = read_delta_logs(delta_paths);
    int maximum = 100;
    for (const auto& p : parameters) maximum = std::max(maximum, 2 * p.bridge_max_chain + 5);
    const auto stable = stable_moments(maximum);
    verify_stable_table(stable_path, stable);
    const auto precomputed_bridges = precompute_stable_bridges(parameters, stable);

    std::vector<RowResult> results(parameters.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t index = 0; index < parameters.size(); ++index) {
      static const std::map<int, mpz_class> empty;
      const auto found = deltas.find(parameters[index].rank);
      const auto* precomputed = parameters[index].exact_through < parameters[index].rank
                                    ? &precomputed_bridges[index]
                                    : nullptr;
      results[index] = verify_row(parameters[index],
                                  found == deltas.end() ? empty : found->second, stable,
                                  precomputed);
    }
    bool ok = true;
    for (size_t index = 0; index < parameters.size(); ++index) {
      const auto& p = parameters[index];
      const auto& result = results[index];
      if (!result.ok) {
        ok = false;
        std::cerr << "D" << p.rank << " FAILED: " << result.error << '\n';
      } else {
        std::cout << "D" << p.rank << " exact finite bridge/tail: ok; min_m="
                  << result.minimum_chain << "; bridge_bits=" << bits(result.minimum_margin)
                  << "; tail_num_bits=" << bits(result.tail_margin.get_num())
                  << "; elapsed=" << result.seconds << " s\n";
        if (print_exact) {
          std::cout << "D_" << p.rank << " minimum bridge margin m="
                    << result.minimum_chain << ": " << result.minimum_margin << '\n';
          std::cout << "D_" << p.rank << " tail margin = " << result.tail_margin << '\n';
        }
      }
    }
    if (!ok) return EXIT_FAILURE;
    const double elapsed = std::chrono::duration<double>(Clock::now() - start).count();
    std::cout << "All selected D high-rank finite-bridge/tail certificates verified "
              << "with exact GMP arithmetic; rows=" << parameters.size()
              << ", OpenMP threads=" << omp_get_max_threads() << ", elapsed=" << elapsed
              << " s.\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "D high-rank GMP certificate failure: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
