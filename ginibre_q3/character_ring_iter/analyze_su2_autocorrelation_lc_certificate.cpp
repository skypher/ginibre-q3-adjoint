#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>
#include <z3++.h>

namespace {

using Pair = std::array<int, 2>;
using Quartet = std::array<int, 4>;
using QuadraticTerms = std::map<Pair, long long>;
using QuarticTerms = std::map<Quartet, long long>;
using Exponent = std::vector<int>;
using Polynomial = std::map<Exponent, long long>;
using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using BigPolynomial = std::map<Exponent, Integer>;
using RationalPolynomial = std::map<Exponent, Rational>;

int parse_positive(const char* text, const std::string& name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value, &consumed, 10);
  if (consumed != value.size() || parsed <= 0) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  return static_cast<int>(parsed);
}

Pair pair_monomial(const int first, const int second) {
  Pair result{first, second};
  std::sort(result.begin(), result.end(), std::greater<int>());
  return result;
}

Quartet product_monomial(const Pair& first, const Pair& second) {
  Quartet result{first[0], first[1], second[0], second[1]};
  std::sort(result.begin(), result.end(), std::greater<int>());
  return result;
}

QuadraticTerms autocorrelation(const int support, const int shift) {
  QuadraticTerms result;
  if (shift > 2 * support) {
    return result;
  }
  for (int index = -support; index + shift <= support; ++index) {
    ++result[pair_monomial(std::abs(index), std::abs(index + shift))];
  }
  return result;
}

QuadraticTerms character_coefficient(const int support, const int label) {
  QuadraticTerms result = autocorrelation(support, label);
  for (const auto& [term, coefficient] :
       autocorrelation(support, label + 1)) {
    result[term] -= coefficient;
  }
  return result;
}

void add_product(QuarticTerms& result, const QuadraticTerms& left,
                 const QuadraticTerms& right, const long long sign) {
  for (const auto& [left_term, left_coefficient] : left) {
    for (const auto& [right_term, right_coefficient] : right) {
      result[product_monomial(left_term, right_term)] +=
          sign * left_coefficient * right_coefficient;
    }
  }
}

QuarticTerms radial_polynomial(const int support, const int antidiagonal,
                               const int depth) {
  const QuadraticTerms c_zero = character_coefficient(support, 0);
  QuadraticTerms top_shells =
      character_coefficient(support, antidiagonal + 1);
  for (const auto& [term, coefficient] :
       character_coefficient(support, antidiagonal + 2)) {
    top_shells[term] += coefficient;
  }
  QuarticTerms result;
  add_product(result, c_zero, top_shells, 1);
  add_product(result, character_coefficient(support, depth),
              character_coefficient(support, antidiagonal - depth), 1);
  add_product(result, character_coefficient(support, depth + 1),
              character_coefficient(support, antidiagonal - depth + 1),
              -1);
  return result;
}

QuadraticTerms wedge_polynomial(const int support, const int label,
                                const int left, const int right) {
  QuadraticTerms result;
  if (left <= support) {
    for (int source = std::abs(right - label);
         source <= right + label; ++source) {
      if (source <= support) {
        ++result[pair_monomial(left, source)];
      }
    }
  }
  if (right <= support) {
    for (int source = std::abs(left - label);
         source <= left + label; ++source) {
      if (source <= support) {
        --result[pair_monomial(right, source)];
      }
    }
  }
  return result;
}

QuarticTerms gap_prefix_polynomial(const int support, const int first_label,
                                   const int second_label,
                                   const int maximum_gap) {
  QuarticTerms result;
  const int maximum_right =
      support + std::min(first_label, second_label);
  for (int left = 0; left <= support; ++left) {
    for (int right = left + 1; right <= maximum_right; ++right) {
      if (right - left > maximum_gap) {
        continue;
      }
      add_product(
          result,
          wedge_polynomial(support, first_label, left, right),
          wedge_polynomial(support, second_label, left, right),
          1);
    }
  }
  return result;
}

Polynomial direct_polynomial(const QuarticTerms& input,
                             const int support) {
  Polynomial result;
  for (const auto& [term, coefficient] : input) {
    Exponent exponent(static_cast<std::size_t>(support + 1), 0);
    for (const int index : term) {
      ++exponent[static_cast<std::size_t>(index)];
    }
    result[std::move(exponent)] += coefficient;
  }
  return result;
}

Polynomial constant_polynomial(const int variables, const long long value) {
  Polynomial result;
  if (value != 0) {
    result.emplace(Exponent(static_cast<std::size_t>(variables), 0), value);
  }
  return result;
}

Polynomial variable_polynomial(const int variables, const int variable) {
  Exponent exponent(static_cast<std::size_t>(variables), 0);
  exponent[static_cast<std::size_t>(variable)] = 1;
  return Polynomial{{std::move(exponent), 1}};
}

void add_scaled(Polynomial& target, const Polynomial& source,
                const long long scale) {
  for (const auto& [exponent, coefficient] : source) {
    target[exponent] += scale * coefficient;
  }
}

Polynomial multiply(const Polynomial& left, const Polynomial& right) {
  Polynomial result;
  for (const auto& [left_exponent, left_coefficient] : left) {
    for (const auto& [right_exponent, right_coefficient] : right) {
      Exponent exponent = left_exponent;
      for (std::size_t index = 0; index < exponent.size(); ++index) {
        exponent[index] += right_exponent[index];
      }
      result[exponent] += left_coefficient * right_coefficient;
    }
  }
  return result;
}

Polynomial substitute_differences(const QuarticTerms& input,
                                  const int support) {
  const int variables = support + 1;
  std::vector<Polynomial> half(static_cast<std::size_t>(variables));
  for (int index = 0; index <= support; ++index) {
    for (int difference = index; difference <= support; ++difference) {
      add_scaled(half[static_cast<std::size_t>(index)],
                 variable_polynomial(variables, difference), 1);
    }
  }
  Polynomial result;
  for (const auto& [term, coefficient] : input) {
    Polynomial expanded = constant_polynomial(variables, coefficient);
    for (const int index : term) {
      expanded =
          multiply(expanded, half[static_cast<std::size_t>(index)]);
    }
    add_scaled(result, expanded, 1);
  }
  return result;
}

std::vector<Polynomial> log_concavity_generators(const int support) {
  const int variables = support + 1;
  std::vector<Polynomial> half(static_cast<std::size_t>(variables));
  for (int index = 0; index <= support; ++index) {
    for (int difference = index; difference <= support; ++difference) {
      add_scaled(half[static_cast<std::size_t>(index)],
                 variable_polynomial(variables, difference), 1);
    }
  }
  std::vector<Polynomial> result;
  for (int index = 1; index < support; ++index) {
    Polynomial generator = multiply(
        half[static_cast<std::size_t>(index)],
        half[static_cast<std::size_t>(index)]);
    add_scaled(
        generator,
        multiply(half[static_cast<std::size_t>(index - 1)],
                 half[static_cast<std::size_t>(index + 1)]),
        -1);
    result.push_back(std::move(generator));
  }
  return result;
}

void enumerate_exponents(const int variables, const int degree,
                         const int index, Exponent& exponent,
                         std::vector<Exponent>& result) {
  if (index + 1 == variables) {
    exponent[static_cast<std::size_t>(index)] = degree;
    result.push_back(exponent);
    return;
  }
  for (int value = 0; value <= degree; ++value) {
    exponent[static_cast<std::size_t>(index)] = value;
    enumerate_exponents(variables, degree - value, index + 1, exponent,
                        result);
  }
}

std::vector<Exponent> homogeneous_exponents(const int variables,
                                            const int degree) {
  std::vector<Exponent> result;
  Exponent exponent(static_cast<std::size_t>(variables), 0);
  enumerate_exponents(variables, degree, 0, exponent, result);
  return result;
}

Polynomial monomial_polynomial(const Exponent& exponent) {
  return Polynomial{{exponent, 1}};
}

BigPolynomial big_constant(const Integer& value) {
  BigPolynomial result;
  if (value != 0) {
    result.emplace(Exponent(3U, 0), value);
  }
  return result;
}

BigPolynomial big_variable(const int variable) {
  Exponent exponent(3U, 0);
  exponent[static_cast<std::size_t>(variable)] = 1;
  return BigPolynomial{{std::move(exponent), 1}};
}

void big_add_scaled(BigPolynomial& target, const BigPolynomial& source,
                    const Integer& scale) {
  for (const auto& [exponent, coefficient] : source) {
    target[exponent] += scale * coefficient;
  }
}

BigPolynomial big_multiply(const BigPolynomial& left,
                           const BigPolynomial& right) {
  BigPolynomial result;
  for (const auto& [left_exponent, left_coefficient] : left) {
    for (const auto& [right_exponent, right_coefficient] : right) {
      Exponent exponent = left_exponent;
      for (std::size_t index = 0; index < exponent.size(); ++index) {
        exponent[index] += right_exponent[index];
      }
      result[exponent] += left_coefficient * right_coefficient;
    }
  }
  return result;
}

BigPolynomial big_power(BigPolynomial base, int exponent) {
  BigPolynomial result = big_constant(1);
  while (exponent > 0) {
    if (exponent % 2 == 1) {
      result = big_multiply(result, base);
    }
    exponent /= 2;
    if (exponent > 0) {
      base = big_multiply(base, base);
    }
  }
  return result;
}

BigPolynomial support_three_parameterization(const Polynomial& target) {
  const BigPolynomial one = big_constant(1);
  const BigPolynomial t = big_variable(0);
  const BigPolynomial x = big_variable(1);
  const BigPolynomial y = big_variable(2);
  BigPolynomial one_plus_t = one;
  big_add_scaled(one_plus_t, t, 1);
  BigPolynomial one_plus_yt = one;
  big_add_scaled(one_plus_yt, big_multiply(y, t), 1);
  const BigPolynomial b =
      big_multiply(big_multiply(y, t), one_plus_t);
  const BigPolynomial a =
      big_multiply(big_multiply(x, b), one_plus_yt);
  const std::array<BigPolynomial, 4> root{a, b, t, one};

  BigPolynomial result;
  for (const auto& [exponent, coefficient] : target) {
    BigPolynomial expanded = big_constant(coefficient);
    for (std::size_t variable = 0; variable < exponent.size(); ++variable) {
      expanded = big_multiply(
          expanded,
          big_power(root[variable], exponent[variable]));
    }
    big_add_scaled(result, expanded, 1);
  }
  return result;
}

BigPolynomial support_three_ratio_parameterization(
    const Polynomial& target) {
  const BigPolynomial one = big_constant(1);
  const BigPolynomial u = big_variable(0);
  const BigPolynomial v = big_variable(1);
  const BigPolynomial w = big_variable(2);
  BigPolynomial one_minus_u = one;
  big_add_scaled(one_minus_u, u, -1);
  const BigPolynomial uv = big_multiply(u, v);
  BigPolynomial one_minus_uv = one;
  big_add_scaled(one_minus_uv, uv, -1);
  const BigPolynomial uvw = big_multiply(uv, w);
  BigPolynomial one_minus_uvw = one;
  big_add_scaled(one_minus_uvw, uvw, -1);

  const BigPolynomial a = one_minus_u;
  const BigPolynomial b = big_multiply(u, one_minus_uv);
  const BigPolynomial c = big_multiply(
      big_multiply(big_power(u, 2), v), one_minus_uvw);
  const BigPolynomial d = big_multiply(
      big_multiply(big_power(u, 3), big_power(v, 2)), w);
  const std::array<BigPolynomial, 4> root{a, b, c, d};

  BigPolynomial result;
  for (const auto& [exponent, coefficient] : target) {
    BigPolynomial expanded = big_constant(coefficient);
    for (std::size_t variable = 0; variable < exponent.size(); ++variable) {
      expanded = big_multiply(
          expanded,
          big_power(root[variable], exponent[variable]));
    }
    big_add_scaled(result, expanded, 1);
  }
  return result;
}

BigPolynomial support_three_ordered_ratio_difference_parameterization(
    const Polynomial& target) {
  const BigPolynomial one = big_constant(1);
  const BigPolynomial first_slack = big_variable(0);
  const BigPolynomial second_slack = big_variable(1);
  const BigPolynomial terminal_ratio = big_variable(2);
  BigPolynomial middle_ratio = second_slack;
  big_add_scaled(middle_ratio, terminal_ratio, 1);
  BigPolynomial initial_ratio = first_slack;
  big_add_scaled(initial_ratio, middle_ratio, 1);
  const std::array<BigPolynomial, 4> root{
      one,
      initial_ratio,
      big_multiply(initial_ratio, middle_ratio),
      big_multiply(
          big_multiply(initial_ratio, middle_ratio),
          terminal_ratio)};

  BigPolynomial result;
  for (const auto& [exponent, coefficient] : target) {
    BigPolynomial expanded = big_constant(coefficient);
    for (std::size_t variable = 0; variable < exponent.size(); ++variable) {
      expanded = big_multiply(
          expanded,
          big_power(root[variable], exponent[variable]));
    }
    big_add_scaled(result, expanded, 1);
  }
  return result;
}

BigPolynomial nd_constant(const int variables, const Integer& value) {
  BigPolynomial result;
  if (value != 0) {
    result.emplace(Exponent(static_cast<std::size_t>(variables), 0), value);
  }
  return result;
}

BigPolynomial nd_variable(const int variables, const int variable) {
  Exponent exponent(static_cast<std::size_t>(variables), 0);
  exponent[static_cast<std::size_t>(variable)] = 1;
  return BigPolynomial{{std::move(exponent), 1}};
}

BigPolynomial nd_power(BigPolynomial base, int exponent) {
  const int variables =
      static_cast<int>(base.begin()->first.size());
  BigPolynomial result = nd_constant(variables, 1);
  while (exponent > 0) {
    if (exponent % 2 == 1) {
      result = big_multiply(result, base);
    }
    exponent /= 2;
    if (exponent > 0) {
      base = big_multiply(base, base);
    }
  }
  return result;
}

BigPolynomial ordered_ratio_difference_parameterization(
    const Polynomial& target, const int support) {
  if (support <= 0) {
    throw std::invalid_argument(
        "ordered-ratio parameterization needs positive support");
  }
  const int variables = support;
  const BigPolynomial one = nd_constant(variables, 1);
  std::vector<BigPolynomial> ratio(
      static_cast<std::size_t>(support),
      nd_constant(variables, 0));
  for (int index = 0; index < support; ++index) {
    for (int slack = index; slack < support; ++slack) {
      big_add_scaled(
          ratio[static_cast<std::size_t>(index)],
          nd_variable(variables, slack),
          1);
    }
  }
  std::vector<BigPolynomial> root(
      static_cast<std::size_t>(support + 1), one);
  for (int index = 1; index <= support; ++index) {
    root[static_cast<std::size_t>(index)] =
        big_multiply(
            root[static_cast<std::size_t>(index - 1)],
            ratio[static_cast<std::size_t>(index - 1)]);
  }

  BigPolynomial result;
  for (const auto& [exponent, coefficient] : target) {
    BigPolynomial expanded = nd_constant(variables, coefficient);
    for (std::size_t variable = 0; variable < exponent.size();
         ++variable) {
      expanded = big_multiply(
          expanded,
          nd_power(root[variable], exponent[variable]));
    }
    big_add_scaled(result, expanded, 1);
  }
  return result;
}

BigPolynomial ratio_cube_parameterization(const Polynomial& target,
                                          const int support) {
  const int variables = support;
  const BigPolynomial one = nd_constant(variables, 1);
  std::vector<BigPolynomial> tail(
      static_cast<std::size_t>(support + 1), one);
  BigPolynomial ratio = one;
  for (int index = 1; index <= support; ++index) {
    ratio = big_multiply(
        ratio, nd_variable(variables, index - 1));
    tail[static_cast<std::size_t>(index)] = big_multiply(
        tail[static_cast<std::size_t>(index - 1)], ratio);
  }
  std::vector<BigPolynomial> root(
      static_cast<std::size_t>(support + 1));
  for (int index = 0; index < support; ++index) {
    root[static_cast<std::size_t>(index)] =
        tail[static_cast<std::size_t>(index)];
    big_add_scaled(
        root[static_cast<std::size_t>(index)],
        tail[static_cast<std::size_t>(index + 1)], -1);
  }
  root[static_cast<std::size_t>(support)] =
      tail[static_cast<std::size_t>(support)];

  BigPolynomial result;
  for (const auto& [exponent, coefficient] : target) {
    BigPolynomial expanded = nd_constant(variables, coefficient);
    for (std::size_t variable = 0; variable < exponent.size(); ++variable) {
      expanded = big_multiply(
          expanded,
          nd_power(root[variable], exponent[variable]));
    }
    big_add_scaled(result, expanded, 1);
  }
  return result;
}

BigPolynomial root_ratio_parameterization(const Polynomial& target,
                                           const int support) {
  const int variables = support;
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial slope = nd_variable(variables, 0);
  std::vector<BigPolynomial> root(
      static_cast<std::size_t>(support + 1), one);
  for (int index = 1; index <= support; ++index) {
    root[static_cast<std::size_t>(index)] =
        nd_power(slope, index);
    for (int ratio = 2; ratio <= index; ++ratio) {
      root[static_cast<std::size_t>(index)] = big_multiply(
          root[static_cast<std::size_t>(index)],
          nd_power(
              nd_variable(variables, ratio - 1),
              index - ratio + 1));
    }
  }

  BigPolynomial result;
  for (const auto& [exponent, coefficient] : target) {
    BigPolynomial expanded = nd_constant(variables, coefficient);
    for (std::size_t variable = 0; variable < exponent.size(); ++variable) {
      expanded = big_multiply(
          expanded,
          nd_power(root[variable], exponent[variable]));
    }
    big_add_scaled(result, expanded, 1);
  }
  return result;
}

Integer binomial(const int n, const int k) {
  if (k < 0 || k > n) {
    return 0;
  }
  Integer result = 1;
  for (int index = 1; index <= std::min(k, n - k); ++index) {
    result *= n - index + 1;
    result /= index;
  }
  return result;
}

BigPolynomial compactify_positive_variable(
    const BigPolynomial& input, const std::size_t variable) {
  if (input.empty() || input.begin()->first.empty()) {
    throw std::invalid_argument(
        "positive-variable compactification needs a variable");
  }
  if (variable >= input.begin()->first.size()) {
    throw std::invalid_argument(
        "positive-variable compactification index out of range");
  }
  int degree = 0;
  for (const auto& [exponent, coefficient] : input) {
    static_cast<void>(coefficient);
    degree = std::max(degree, exponent[variable]);
  }
  BigPolynomial result;
  for (const auto& [exponent, coefficient] : input) {
    const int complement = degree - exponent[variable];
    for (int power = 0; power <= complement; ++power) {
      Exponent transformed = exponent;
      transformed[variable] += power;
      const Integer sign = power % 2 == 0 ? 1 : -1;
      result[std::move(transformed)] +=
          coefficient * sign * binomial(complement, power);
    }
  }
  return result;
}

BigPolynomial half_shift_positive_variable(
    const BigPolynomial& input, const std::size_t variable) {
  if (input.empty() || input.begin()->first.empty()) {
    throw std::invalid_argument(
        "half-shift needs a nonempty polynomial");
  }
  if (variable >= input.begin()->first.size()) {
    throw std::invalid_argument("half-shift index out of range");
  }
  int degree = 0;
  for (const auto& [exponent, coefficient] : input) {
    if (coefficient != 0) {
      degree = std::max(degree, exponent[variable]);
    }
  }
  BigPolynomial result;
  for (const auto& [exponent, coefficient] : input) {
    if (coefficient == 0) {
      continue;
    }
    Integer denominator_clear = 1;
    for (int power = exponent[variable]; power < degree; ++power) {
      denominator_clear *= 2;
    }
    for (int power = 0; power <= exponent[variable]; ++power) {
      Exponent transformed = exponent;
      transformed[variable] = power;
      result[std::move(transformed)] +=
          coefficient * denominator_clear
          * binomial(exponent[variable], power);
    }
  }
  return result;
}

BigPolynomial reflect_unit_variable(
    const BigPolynomial& input, const std::size_t variable) {
  if (input.empty() || input.begin()->first.empty()) {
    throw std::invalid_argument(
        "unit reflection needs a nonempty polynomial");
  }
  if (variable >= input.begin()->first.size()) {
    throw std::invalid_argument("unit reflection index out of range");
  }
  BigPolynomial result;
  for (const auto& [exponent, coefficient] : input) {
    if (coefficient == 0) {
      continue;
    }
    for (int power = 0; power <= exponent[variable]; ++power) {
      Exponent transformed = exponent;
      transformed[variable] = power;
      const Integer sign = power % 2 == 0 ? 1 : -1;
      result[std::move(transformed)] +=
          coefficient * sign
          * binomial(exponent[variable], power);
    }
  }
  return result;
}

BigPolynomial compactify_first_positive_variable(
    const BigPolynomial& input) {
  return compactify_positive_variable(input, 0U);
}

BigPolynomial remove_common_monomial_factor(
    const BigPolynomial& input) {
  if (input.empty()) {
    return input;
  }
  Exponent common = input.begin()->first;
  for (const auto& [exponent, coefficient] : input) {
    if (coefficient == 0) {
      continue;
    }
    for (std::size_t variable = 0; variable < common.size(); ++variable) {
      common[variable] = std::min(common[variable], exponent[variable]);
    }
  }
  BigPolynomial result;
  for (const auto& [exponent, coefficient] : input) {
    if (coefficient == 0) {
      continue;
    }
    Exponent reduced = exponent;
    for (std::size_t variable = 0; variable < common.size(); ++variable) {
      reduced[variable] -= common[variable];
    }
    result[std::move(reduced)] += coefficient;
  }
  return result;
}

std::vector<RationalPolynomial> outer_bernstein_coefficients(
    const BigPolynomial& input, const int variables) {
  if (variables < 1) {
    throw std::invalid_argument(
        "outer Bernstein conversion needs at least one variable");
  }
  const std::size_t outer =
      static_cast<std::size_t>(variables - 1);
  int degree = 0;
  for (const auto& [exponent, coefficient] : input) {
    static_cast<void>(coefficient);
    degree = std::max(degree, exponent[outer]);
  }
  std::vector<RationalPolynomial> result(
      static_cast<std::size_t>(degree + 1));
  for (const auto& [exponent, coefficient] : input) {
    const int power = exponent[outer];
    Exponent reduced(exponent.begin(), exponent.end() - 1);
    for (int index = power; index <= degree; ++index) {
      result[static_cast<std::size_t>(index)][reduced] +=
          Rational(
              coefficient * binomial(index, power),
              binomial(degree, power));
    }
  }
  return result;
}

std::vector<RationalPolynomial> elevate_bernstein_coefficients(
    const std::vector<RationalPolynomial>& input,
    const int target_degree) {
  const int source_degree = static_cast<int>(input.size()) - 1;
  if (source_degree < 0 || target_degree < source_degree) {
    throw std::invalid_argument(
        "invalid Bernstein degree elevation");
  }
  std::vector<RationalPolynomial> result(
      static_cast<std::size_t>(target_degree + 1));
  for (int index = 0; index <= target_degree; ++index) {
    for (int source = 0; source <= source_degree; ++source) {
      const int remainder = index - source;
      if (
          remainder < 0
          || remainder > target_degree - source_degree
      ) {
        continue;
      }
      const Rational scale(
          binomial(source_degree, source)
              * binomial(
                  target_degree - source_degree,
                  remainder),
          binomial(target_degree, index));
      for (const auto& [exponent, coefficient] :
           input[static_cast<std::size_t>(source)]) {
        result[static_cast<std::size_t>(index)][exponent] +=
            scale * coefficient;
      }
    }
  }
  return result;
}

std::vector<RationalPolynomial> outer_power_coefficients(
    const BigPolynomial& input, const int variables) {
  if (variables < 1) {
    throw std::invalid_argument(
        "outer power decomposition needs at least one variable");
  }
  const std::size_t outer =
      static_cast<std::size_t>(variables - 1);
  int degree = 0;
  for (const auto& [exponent, coefficient] : input) {
    static_cast<void>(coefficient);
    degree = std::max(degree, exponent[outer]);
  }
  std::vector<RationalPolynomial> result(
      static_cast<std::size_t>(degree + 1));
  for (const auto& [exponent, coefficient] : input) {
    Exponent reduced(exponent.begin(), exponent.end() - 1);
    result[static_cast<std::size_t>(exponent[outer])][reduced] +=
        Rational(coefficient);
  }
  return result;
}

RationalPolynomial rational_polynomial(
    const BigPolynomial& input) {
  RationalPolynomial result;
  for (const auto& [exponent, coefficient] : input) {
    result[exponent] += Rational(coefficient);
  }
  return result;
}

Rational evaluate_polynomial(
    const RationalPolynomial& input,
    const std::vector<Rational>& point) {
  Rational result = 0;
  for (const auto& [exponent, coefficient] : input) {
    if (exponent.size() != point.size()) {
      throw std::invalid_argument(
          "polynomial evaluation dimension mismatch");
    }
    Rational term = coefficient;
    for (std::size_t variable = 0; variable < exponent.size();
         ++variable) {
      for (int power = 0; power < exponent[variable]; ++power) {
        term *= point[variable];
      }
    }
    result += term;
  }
  return result;
}

Integer evaluate_append_power(
    const Polynomial& input, const int append_variable,
    const int append_power, const std::vector<Integer>& old_profile) {
  if (
      append_variable < 0
      || static_cast<std::size_t>(append_variable)
             != old_profile.size()
  ) {
    throw std::invalid_argument(
        "append-power evaluation dimension mismatch");
  }
  Integer result = 0;
  for (const auto& [exponent, coefficient] : input) {
    if (
        exponent[static_cast<std::size_t>(append_variable)]
        != append_power
    ) {
      continue;
    }
    Integer term = coefficient;
    for (int variable = 0; variable < append_variable; ++variable) {
      for (
          int power = 0;
          power < exponent[static_cast<std::size_t>(variable)];
          ++power
      ) {
        term *= old_profile[static_cast<std::size_t>(variable)];
      }
    }
    result += term;
  }
  return result;
}

Integer evaluate_coordinate_power(
    const Polynomial& input, const int selected_variable,
    const int selected_power,
    const std::vector<Integer>& remaining_profile) {
  if (
      input.empty()
      || selected_variable < 0
      || static_cast<std::size_t>(selected_variable)
             >= input.begin()->first.size()
      || remaining_profile.size() + 1U
             != input.begin()->first.size()
  ) {
    throw std::invalid_argument(
        "coordinate-power evaluation dimension mismatch");
  }
  Integer result = 0;
  for (const auto& [exponent, coefficient] : input) {
    if (
        exponent[static_cast<std::size_t>(selected_variable)]
        != selected_power
    ) {
      continue;
    }
    Integer term = coefficient;
    std::size_t remaining_variable = 0U;
    for (std::size_t variable = 0; variable < exponent.size();
         ++variable) {
      if (static_cast<int>(variable) == selected_variable) {
        continue;
      }
      for (int power = 0; power < exponent[variable]; ++power) {
        term *= remaining_profile[remaining_variable];
      }
      ++remaining_variable;
    }
    result += term;
  }
  return result;
}

bool equal_polynomials(const RationalPolynomial& left,
                       const RationalPolynomial& right) {
  for (const auto& [exponent, coefficient] : left) {
    const auto iterator = right.find(exponent);
    const Rational other =
        iterator == right.end() ? Rational(0) : iterator->second;
    if (coefficient != other) {
      return false;
    }
  }
  for (const auto& [exponent, coefficient] : right) {
    const auto iterator = left.find(exponent);
    const Rational other =
        iterator == left.end() ? Rational(0) : iterator->second;
    if (coefficient != other) {
      return false;
    }
  }
  return true;
}

std::pair<int, Rational> leading_at_last_one(
    const RationalPolynomial& polynomial) {
  std::map<int, Rational> univariate;
  for (const auto& [exponent, coefficient] : polynomial) {
    univariate[exponent.front()] += coefficient;
  }
  for (auto iterator = univariate.rbegin();
       iterator != univariate.rend(); ++iterator) {
    if (iterator->second != 0) {
      return *iterator;
    }
  }
  return {-1, Rational(0)};
}

int replay_gap_prefix_append_low_coefficients() {
  constexpr int support = 3;
  const Polynomial target = direct_polynomial(
      gap_prefix_polynomial(support, 1, 2, 1), support);
  const std::vector<Integer> linear_profile{5, 10, 18};
  const std::vector<Integer> quadratic_profile{1, 2, 4};
  const Integer linear =
      evaluate_append_power(target, support, 1, linear_profile);
  const Integer constant =
      evaluate_append_power(target, support, 0, quadratic_profile);
  const Integer quadratic_linear =
      evaluate_append_power(target, support, 1, quadratic_profile);
  const Integer quadratic =
      evaluate_append_power(target, support, 2, quadratic_profile);
  const Integer cubic =
      evaluate_append_power(target, support, 3, quadratic_profile);
  const Integer quartic =
      evaluate_append_power(target, support, 4, quadratic_profile);
  const Integer append = 3;
  const Integer low_truncation =
      constant
      + append * quadratic_linear
      + append * append * quadratic;
  const Integer full =
      low_truncation
      + append * append * append * cubic
      + append * append * append * append * quartic;
  const std::vector<Integer> cubic_profile{5, 15, 18};
  const Integer cubic_constant =
      evaluate_append_power(target, support, 0, cubic_profile);
  const Integer cubic_linear =
      evaluate_append_power(target, support, 1, cubic_profile);
  const Integer cubic_quadratic =
      evaluate_append_power(target, support, 2, cubic_profile);
  const Integer cubic_cubic =
      evaluate_append_power(target, support, 3, cubic_profile);
  const Integer cubic_quartic =
      evaluate_append_power(target, support, 4, cubic_profile);
  const Integer cubic_append = 15;
  const Integer cubic_truncation =
      cubic_constant
      + cubic_append * cubic_linear
      + cubic_append * cubic_append * cubic_quadratic
      + cubic_append * cubic_append * cubic_append * cubic_cubic;
  const Integer cubic_full =
      cubic_truncation
      + cubic_append * cubic_append * cubic_append * cubic_append
            * cubic_quartic;
  const Polynomial unbounded_target = direct_polynomial(
      gap_prefix_polynomial(support, 1, 6, 1), support);
  const std::vector<Integer> unbounded_profile{1, 3, 1};
  std::array<Integer, 5> unbounded_coefficients{};
  for (int power = 0; power <= 4; ++power) {
    unbounded_coefficients[static_cast<std::size_t>(power)] =
        evaluate_append_power(
            unbounded_target, support, power, unbounded_profile);
  }
  const Integer unbounded_value =
      unbounded_coefficients[0]
      + unbounded_coefficients[1]
      + unbounded_coefficients[2]
      + unbounded_coefficients[3]
      + unbounded_coefficients[4];
  if (
      linear != -8803
      || constant != 538
      || quadratic_linear != -81
      || quadratic != -38
      || cubic != 9
      || quartic != 2
      || low_truncation != -47
      || full != 358
      || cubic_constant != 200667
      || cubic_linear != -15358
      || cubic_quadratic != -663
      || cubic_cubic != 53
      || cubic_quartic != 2
      || cubic_truncation != -3
      || cubic_full != 101247
      || unbounded_coefficients
             != std::array<Integer, 5>{0, 1, -3, 0, 1}
      || unbounded_value != -1
  ) {
    throw std::runtime_error(
        "gap-prefix append coefficient replay mismatch");
  }
  std::cout
      << "SU2_GAP_PREFIX_APPEND_LOW_COEFFICIENT_OBSTRUCTIONS"
      << " support=3 R=1 S=2 D=1"
      << " C1_at_5_10_18=-8803"
      << " coefficients_at_1_2_4=(538,-81,-38,9,2)"
      << " append=3 low_truncation=-47 full=358"
      << " coefficients_at_5_15_18=(200667,-15358,-663,53,2)"
      << " append=15 cubic_truncation=-3 full=101247"
      << " unbounded_append_profile=(1,3,1,1)"
      << " unbounded_labels=(1,6,1)"
      << " unbounded_coefficients=(0,1,-3,0,1)"
      << " unbounded_value=-1"
      << " result=PASS_EXACT"
      << '\n';
  return EXIT_SUCCESS;
}

int replay_gap_prefix_wall_cubic_reserve() {
  {
    constexpr int support = 2;
    const Polynomial target = direct_polynomial(
        gap_prefix_polynomial(support, 1, 2, 1), support);
    const std::vector<Integer> tail{1, 1};
    std::array<Integer, 5> coefficient{};
    for (int power = 0; power <= 4; ++power) {
      coefficient[static_cast<std::size_t>(power)] =
          evaluate_coordinate_power(
              target, 0, power, tail);
    }
    const Integer wall = 1;
    const Integer linear =
        coefficient[0] + wall * coefficient[1];
    const Integer full =
        linear
        + wall * wall * coefficient[2]
        + wall * wall * wall * coefficient[3];
    if (
        coefficient != std::array<Integer, 5>{3, -4, 2, 2, 0}
        || linear != -1
        || full != 3
    ) {
      throw std::runtime_error(
          "wall linear-reserve replay mismatch");
    }
  }

  constexpr int support = 4;
  const Polynomial target = direct_polynomial(
      gap_prefix_polynomial(support, 1, 2, 1), support);
  const std::vector<Integer> tail{3, 3, 2, 1};
  std::array<Integer, 5> coefficient{};
  for (int power = 0; power <= 4; ++power) {
    coefficient[static_cast<std::size_t>(power)] =
        evaluate_coordinate_power(
            target, 0, power, tail);
  }
  const Integer wall = 2;
  const Integer quadratic =
      coefficient[0]
      + wall * coefficient[1]
      + wall * wall * coefficient[2];
  const Integer full =
      quadratic + wall * wall * wall * coefficient[3];
  const Integer without_constant = full - coefficient[0];
  const Integer without_quadratic =
      full - wall * wall * coefficient[2];
  if (
      coefficient != std::array<Integer, 5>{120, -121, 30, 8, 0}
      || quadratic != -2
      || full != 62
      || without_constant != -58
      || without_quadratic != -58
  ) {
    throw std::runtime_error(
        "wall cubic-reserve replay mismatch");
  }
  std::cout
      << "SU2_GAP_PREFIX_WALL_CUBIC_RESERVE"
      << " labels=(1,2,1)"
      << " linear_profile=(1,1,1)"
      << " linear_coefficients=(3,-4,2,2,0)"
      << " linear_truncation=-1 linear_full=3"
      << " cubic_profile=(2,3,3,2,1)"
      << " cubic_coefficients=(120,-121,30,8,0)"
      << " quadratic_truncation=-2 cubic_full=62"
      << " without_constant=-58 without_quadratic=-58"
      << " result=PASS_EXACT"
      << '\n';
  return EXIT_SUCCESS;
}

int replay_wall_121_current_normal_form() {
  constexpr int minimum_support = 2;
  constexpr int maximum_support = 12;
  for (int support = minimum_support;
       support <= maximum_support; ++support) {
    const int variables = support + 1;
    const Polynomial target = direct_polynomial(
        gap_prefix_polynomial(support, 1, 2, 1), support);
    std::array<BigPolynomial, 4> actual;
    for (const auto& [exponent, coefficient] : target) {
      if (coefficient == 0 || exponent[0] > 3) {
        continue;
      }
      Exponent tail_exponent = exponent;
      const int power = tail_exponent[0];
      tail_exponent[0] = 0;
      actual[static_cast<std::size_t>(power)]
            [std::move(tail_exponent)] += coefficient;
    }

    std::vector<BigPolynomial> p(
        static_cast<std::size_t>(support + 4),
        nd_constant(variables, 0));
    for (int index = 1; index <= support; ++index) {
      p[static_cast<std::size_t>(index)] =
          nd_variable(variables, index);
    }
    const auto polynomial_power = [variables](
        const BigPolynomial& base, const int exponent) {
      if (base.empty()) {
        return nd_constant(variables, exponent == 0 ? 1 : 0);
      }
      return nd_power(base, exponent);
    };
    const BigPolynomial& a = p[1];
    const BigPolynomial& b = p[2];
    const BigPolynomial& c = p[3];
    const BigPolynomial& d = p[4];
    std::array<BigPolynomial, 4> expected;
    big_add_scaled(
        expected[1],
        big_multiply(polynomial_power(a, 2), b),
        -3);
    big_add_scaled(
        expected[1],
        big_multiply(a, polynomial_power(b, 2)),
        -2);
    big_add_scaled(
        expected[1],
        big_multiply(big_multiply(a, b), c),
        -1);
    big_add_scaled(
        expected[1],
        big_multiply(big_multiply(a, b), d),
        -1);
    big_add_scaled(
        expected[1],
        big_multiply(big_multiply(b, c), d),
        -1);
    big_add_scaled(expected[1], polynomial_power(b, 3), 1);
    big_add_scaled(
        expected[1],
        big_multiply(a, polynomial_power(c, 2)),
        1);
    big_add_scaled(expected[1], polynomial_power(c, 3), 1);
    big_add_scaled(expected[2], polynomial_power(a, 2), 1);
    big_add_scaled(expected[2], polynomial_power(b, 2), 1);
    big_add_scaled(
        expected[2], big_multiply(a, c), 1);
    big_add_scaled(
        expected[2], big_multiply(b, c), 1);
    big_add_scaled(expected[3], a, 1);
    big_add_scaled(expected[3], b, 1);
    big_add_scaled(expected[3], c, 1);

    std::vector<BigPolynomial> g(
        static_cast<std::size_t>(support + 2),
        nd_constant(variables, 0));
    for (int index = 1; index <= support + 1; ++index) {
      g[static_cast<std::size_t>(index)] =
          polynomial_power(
              p[static_cast<std::size_t>(index)], 2);
      big_add_scaled(
          g[static_cast<std::size_t>(index)],
          big_multiply(
              p[static_cast<std::size_t>(index - 1)],
              p[static_cast<std::size_t>(index + 1)]),
          -1);
    }
    std::vector<BigPolynomial> h(
        static_cast<std::size_t>(support + 3),
        nd_constant(variables, 0));
    for (int index = 1; index <= support + 2; ++index) {
      h[static_cast<std::size_t>(index)] =
          big_multiply(
              p[static_cast<std::size_t>(index)],
              p[static_cast<std::size_t>(index - 1)]);
      if (index >= 2) {
        big_add_scaled(
            h[static_cast<std::size_t>(index)],
            big_multiply(
                p[static_cast<std::size_t>(index - 2)],
                p[static_cast<std::size_t>(index + 1)]),
          -1);
      }
    }
    for (int index = 1; index <= support; ++index) {
      BigPolynomial identity = big_multiply(
          p[static_cast<std::size_t>(index)],
          h[static_cast<std::size_t>(index)]);
      big_add_scaled(
          identity,
          big_multiply(
              p[static_cast<std::size_t>(index - 1)],
              g[static_cast<std::size_t>(index)]),
          -1);
      big_add_scaled(
          identity,
          big_multiply(
              p[static_cast<std::size_t>(index + 1)],
              g[static_cast<std::size_t>(index - 1)]),
          -1);
      if (std::any_of(
              identity.begin(), identity.end(),
              [](const auto& term) { return term.second != 0; })) {
        throw std::runtime_error(
            "wall (1,2,1) feasible-Jacobi identity mismatch at support "
            + std::to_string(support)
            + " and index " + std::to_string(index));
      }
    }
    std::vector<BigPolynomial> k(
        static_cast<std::size_t>(support + 2),
        nd_constant(variables, 0));
    for (int index = 1; index <= support + 1; ++index) {
      k[static_cast<std::size_t>(index)] =
          g[static_cast<std::size_t>(index)];
      big_add_scaled(
          k[static_cast<std::size_t>(index)],
          h[static_cast<std::size_t>(index)],
          1);
      big_add_scaled(
          k[static_cast<std::size_t>(index)],
          h[static_cast<std::size_t>(index + 1)],
          1);
    }
    big_add_scaled(
        expected[0],
        big_multiply(polynomial_power(a, 3), b),
        1);
    BigPolynomial paired_current;
    for (int index = 1; index <= support; ++index) {
      BigPolynomial delta_g =
          g[static_cast<std::size_t>(index)];
      big_add_scaled(
          delta_g,
          g[static_cast<std::size_t>(index + 1)],
          -1);
      BigPolynomial delta_k =
          k[static_cast<std::size_t>(index)];
      big_add_scaled(
          delta_k,
          k[static_cast<std::size_t>(index + 1)],
          -1);
      big_add_scaled(
          paired_current, big_multiply(delta_g, delta_k), 1);
    }
    big_add_scaled(expected[0], paired_current, 1);

    const auto laurent_ratio = [variables, support](
        const int ratio_index, const bool reciprocal) {
      if (ratio_index < 2 || ratio_index > support) {
        return nd_constant(variables, 0);
      }
      Exponent exponent(static_cast<std::size_t>(variables), 0);
      const int numerator =
          reciprocal ? ratio_index - 1 : ratio_index;
      const int denominator =
          reciprocal ? ratio_index : ratio_index - 1;
      ++exponent[static_cast<std::size_t>(numerator)];
      --exponent[static_cast<std::size_t>(denominator)];
      return BigPolynomial{{std::move(exponent), 1}};
    };
    BigPolynomial band_form =
        big_multiply(g[1], g[1]);
    for (int index = 2; index <= support; ++index) {
      BigPolynomial diagonal = nd_constant(variables, 2);
      big_add_scaled(
          diagonal, laurent_ratio(index, true), 1);
      big_add_scaled(
          diagonal, laurent_ratio(index + 2, false), 1);
      big_add_scaled(
          band_form,
          big_multiply(
              diagonal,
              big_multiply(
                  g[static_cast<std::size_t>(index)],
                  g[static_cast<std::size_t>(index)])),
          1);
    }
    if (support >= 2) {
      BigPolynomial adjacent = nd_constant(variables, -2);
      big_add_scaled(adjacent, laurent_ratio(3, false), 1);
      big_add_scaled(adjacent, laurent_ratio(4, false), -1);
      big_add_scaled(
          band_form,
          big_multiply(
              adjacent, big_multiply(g[1], g[2])),
          1);
    }
    for (int index = 2; index < support; ++index) {
      BigPolynomial adjacent = nd_constant(variables, -2);
      big_add_scaled(
          adjacent, laurent_ratio(index + 1, true), 1);
      big_add_scaled(
          adjacent, laurent_ratio(index, true), -1);
      big_add_scaled(
          adjacent, laurent_ratio(index + 2, false), 1);
      big_add_scaled(
          adjacent, laurent_ratio(index + 3, false), -1);
      big_add_scaled(
          band_form,
          big_multiply(
              adjacent,
              big_multiply(
                  g[static_cast<std::size_t>(index)],
                  g[static_cast<std::size_t>(index + 1)])),
          1);
    }
    for (int index = 1; index + 2 <= support; ++index) {
      BigPolynomial second_neighbor =
          laurent_ratio(index + 2, false);
      big_add_scaled(
          second_neighbor,
          laurent_ratio(index + 2, true),
          1);
      big_add_scaled(
          band_form,
          big_multiply(
              second_neighbor,
              big_multiply(
                  g[static_cast<std::size_t>(index)],
                  g[static_cast<std::size_t>(index + 2)])),
          -1);
    }

    const auto canonicalize = [](
        const BigPolynomial& polynomial) {
      BigPolynomial result;
      for (const auto& [exponent, coefficient] : polynomial) {
        if (coefficient != 0) {
          result.emplace(exponent, coefficient);
        }
      }
      return result;
    };
    for (std::size_t power = 0; power < actual.size(); ++power) {
      if (
          canonicalize(actual[power])
          != canonicalize(expected[power])
      ) {
        throw std::runtime_error(
            "wall (1,2,1) current normal-form mismatch at support "
            + std::to_string(support)
            + " and power " + std::to_string(power));
      }
    }
    if (canonicalize(band_form) != canonicalize(paired_current)) {
      throw std::runtime_error(
          "wall (1,2,1) ratio-band identity mismatch at support "
          + std::to_string(support));
    }
  }
  const std::array<Integer, 2> obstruction_u{7, 9};
  const std::array<Integer, 2> obstruction_w{0, 12};
  Integer obstruction_u_norm = 0;
  Integer obstruction_w_norm = 0;
  Integer obstruction_pairing = 0;
  for (std::size_t index = 0; index < obstruction_u.size(); ++index) {
    obstruction_u_norm += obstruction_u[index] * obstruction_u[index];
    obstruction_w_norm += obstruction_w[index] * obstruction_w[index];
    obstruction_pairing +=
        obstruction_u[index]
        * (obstruction_u[index] + obstruction_w[index]);
  }
  if (
      obstruction_u_norm != 130
      || obstruction_w_norm != 144
      || obstruction_pairing != 238
  ) {
    throw std::runtime_error(
        "wall (1,2,1) norm-contraction obstruction mismatch");
  }
  const std::array<Integer, 3> remainder_u{1, -1, 1};
  const std::array<Integer, 3> remainder_w{-1, 1, 1};
  Integer remainder_square = 0;
  Integer remainder_cross = 0;
  for (std::size_t index = 0; index < remainder_u.size(); ++index) {
    remainder_square += remainder_u[index] * remainder_u[index];
    remainder_cross += remainder_u[index] * remainder_w[index];
  }
  if (remainder_square != 3 || remainder_cross != -1) {
    throw std::runtime_error(
        "wall (1,2,1) signed-Jacobi remainder mismatch");
  }
  const int geometric_variables = 1;
  const BigPolynomial geometric_one =
      nd_constant(geometric_variables, 1);
  const BigPolynomial geometric_r =
      nd_variable(geometric_variables, 0);
  const BigPolynomial geometric_r2 =
      nd_power(geometric_r, 2);
  const BigPolynomial geometric_r3 =
      nd_power(geometric_r, 3);
  BigPolynomial geometric_a = geometric_one;
  big_add_scaled(geometric_a, geometric_r2, 2);
  big_add_scaled(geometric_a, geometric_r3, 1);
  BigPolynomial geometric_b =
      nd_constant(geometric_variables, 0);
  big_add_scaled(geometric_b, geometric_r, 3);
  big_add_scaled(geometric_b, geometric_r2, 2);
  BigPolynomial geometric_c = geometric_one;
  big_add_scaled(geometric_c, geometric_r, 1);
  big_add_scaled(geometric_c, geometric_r2, 1);
  BigPolynomial geometric_d =
      big_multiply(geometric_a, geometric_a);
  big_add_scaled(
      geometric_d,
      big_multiply(geometric_b, geometric_c),
      3);
  BigPolynomial geometric_one_plus_r = geometric_one;
  big_add_scaled(geometric_one_plus_r, geometric_r, 1);
  BigPolynomial geometric_n =
      big_multiply(geometric_a, geometric_b);
  big_add_scaled(
      geometric_n,
      big_multiply(geometric_c, geometric_one_plus_r),
      9);
  BigPolynomial geometric_certificate =
      big_multiply(
          geometric_c,
          big_multiply(geometric_n, geometric_n));
  for (auto& [exponent, coefficient] : geometric_certificate) {
    static_cast<void>(exponent);
    coefficient *= 3;
  }
  big_add_scaled(
      geometric_certificate,
      big_multiply(
          geometric_a,
          big_multiply(geometric_n, geometric_d)),
      4);
  big_add_scaled(
      geometric_certificate,
      big_multiply(
          geometric_b,
          big_multiply(geometric_d, geometric_d)),
      -4);
  const std::array<Integer, 11> geometric_expected{
      279, 1773, 5148, 9135, 11061, 9576,
      6057, 2835, 972, 225, 27};
  std::array<Integer, 11> geometric_actual{};
  for (const auto& [exponent, coefficient] : geometric_certificate) {
    if (coefficient == 0) {
      continue;
    }
    if (exponent.size() != 1U
        || exponent[0] < 0
        || exponent[0] >= static_cast<int>(geometric_actual.size())) {
      throw std::runtime_error(
          "wall (1,2,1) geometric certificate degree mismatch");
    }
    geometric_actual[static_cast<std::size_t>(exponent[0])] =
        coefficient;
  }
  if (geometric_actual != geometric_expected) {
    throw std::runtime_error(
        "wall (1,2,1) geometric certificate mismatch");
  }
  const std::array<Integer, 8> cutoff_profile{
      0, 8, 8, 8, 6, 4, 2, 1};
  constexpr int cutoff_support = 7;
  std::array<Integer, 9> cutoff_g{};
  std::array<Integer, 10> cutoff_h{};
  std::array<Integer, 9> cutoff_k{};
  for (int index = 1; index <= cutoff_support; ++index) {
    const Integer next =
        index < cutoff_support
            ? cutoff_profile[static_cast<std::size_t>(index + 1)]
            : Integer{0};
    cutoff_g[static_cast<std::size_t>(index)] =
        cutoff_profile[static_cast<std::size_t>(index)]
            * cutoff_profile[static_cast<std::size_t>(index)]
        - cutoff_profile[static_cast<std::size_t>(index - 1)] * next;
  }
  for (int index = 1; index <= cutoff_support + 1; ++index) {
    const Integer current =
        index <= cutoff_support
            ? cutoff_profile[static_cast<std::size_t>(index)]
            : Integer{0};
    const Integer previous =
        cutoff_profile[static_cast<std::size_t>(index - 1)];
    const Integer left_two =
        index >= 2
            ? cutoff_profile[static_cast<std::size_t>(index - 2)]
            : Integer{0};
    const Integer next =
        index < cutoff_support
            ? cutoff_profile[static_cast<std::size_t>(index + 1)]
            : Integer{0};
    cutoff_h[static_cast<std::size_t>(index)] =
        current * previous - left_two * next;
  }
  for (int index = 1; index <= cutoff_support; ++index) {
    cutoff_k[static_cast<std::size_t>(index)] =
        cutoff_g[static_cast<std::size_t>(index)]
        + cutoff_h[static_cast<std::size_t>(index)]
        + cutoff_h[static_cast<std::size_t>(index + 1)];
  }
  Integer cutoff_full = 0;
  Integer cutoff_suffix = 0;
  for (int index = 1; index <= cutoff_support; ++index) {
    const Integer product =
        (cutoff_g[static_cast<std::size_t>(index)]
         - cutoff_g[static_cast<std::size_t>(index + 1)])
        * (cutoff_k[static_cast<std::size_t>(index)]
           - cutoff_k[static_cast<std::size_t>(index + 1)]);
    cutoff_full += product;
    if (index >= 2) {
      cutoff_suffix += product;
    }
  }
  if (
      cutoff_g
          != std::array<Integer, 9>{0, 64, 0, 16, 4, 4, 0, 1, 0}
      || cutoff_k
          != std::array<Integer, 9>{
              0, 128, 80, 48, 28, 14, 4, 3, 0}
      || cutoff_suffix != -230
      || cutoff_full != 2842
  ) {
    throw std::runtime_error(
        "wall (1,2,1) cutoff-current obstruction mismatch");
  }
  const auto potential_increment_coefficients = [](
      const std::array<Integer, 5>& profile,
      const int cutoff) {
    constexpr int support = 4;
    std::array<Integer, 6> g{};
    std::array<Integer, 7> h{};
    std::array<Integer, 6> k{};
    for (int index = 1; index <= support; ++index) {
      const Integer next =
          index < support
              ? profile[static_cast<std::size_t>(index + 1)]
              : Integer{0};
      g[static_cast<std::size_t>(index)] =
          profile[static_cast<std::size_t>(index)]
              * profile[static_cast<std::size_t>(index)]
          - profile[static_cast<std::size_t>(index - 1)] * next;
    }
    for (int index = 1; index <= support + 1; ++index) {
      const Integer current =
          index <= support
              ? profile[static_cast<std::size_t>(index)]
              : Integer{0};
      const Integer previous =
          profile[static_cast<std::size_t>(index - 1)];
      const Integer left_two =
          index >= 2
              ? profile[static_cast<std::size_t>(index - 2)]
              : Integer{0};
      const Integer next =
          index < support
              ? profile[static_cast<std::size_t>(index + 1)]
              : Integer{0};
      h[static_cast<std::size_t>(index)] =
          current * previous - left_two * next;
    }
    for (int index = 1; index <= support; ++index) {
      k[static_cast<std::size_t>(index)] =
          g[static_cast<std::size_t>(index)]
          + h[static_cast<std::size_t>(index)]
          + h[static_cast<std::size_t>(index + 1)];
    }
    const Integer base =
        (g[static_cast<std::size_t>(cutoff)]
         - g[static_cast<std::size_t>(cutoff + 1)])
        * (k[static_cast<std::size_t>(cutoff)]
           - k[static_cast<std::size_t>(cutoff + 1)]);
    const Integer coefficient_a =
        g[static_cast<std::size_t>(cutoff)]
              * k[static_cast<std::size_t>(cutoff + 1)]
        - g[static_cast<std::size_t>(cutoff - 1)]
              * k[static_cast<std::size_t>(cutoff)];
    const Integer coefficient_b =
        g[static_cast<std::size_t>(cutoff + 1)]
              * k[static_cast<std::size_t>(cutoff)]
        - g[static_cast<std::size_t>(cutoff)]
              * k[static_cast<std::size_t>(cutoff - 1)];
    const Integer coefficient_c =
        g[static_cast<std::size_t>(cutoff + 1)]
              * k[static_cast<std::size_t>(cutoff + 1)]
        - g[static_cast<std::size_t>(cutoff)]
              * k[static_cast<std::size_t>(cutoff)];
    return std::array<Integer, 4>{
        base, coefficient_a, coefficient_b, coefficient_c};
  };
  const std::array<std::array<Integer, 5>, 4> potential_profiles{
      std::array<Integer, 5>{0, 2, 3, 3, 2},
      std::array<Integer, 5>{0, 3, 3, 2, 1},
      std::array<Integer, 5>{0, 3, 3, 3, 1},
      std::array<Integer, 5>{0, 3, 3, 3, 2}};
  const std::array<int, 4> potential_cutoffs{3, 3, 2, 3};
  const std::array<std::array<Integer, 4>, 4> potential_expected{
      std::array<Integer, 4>{-4, -12, 14, -2},
      std::array<Integer, 4>{0, -15, -9, -3},
      std::array<Integer, 4>{0, -135, 90, 90},
      std::array<Integer, 4>{-2, 30, 12, 4}};
  std::array<std::array<Integer, 4>, 4> potential_actual{};
  for (std::size_t index = 0; index < potential_profiles.size();
       ++index) {
    potential_actual[index] = potential_increment_coefficients(
        potential_profiles[index], potential_cutoffs[index]);
  }
  const std::array<Integer, 4> farkas_weights{90, 1044, 10, 603};
  std::array<Integer, 4> farkas_sum{};
  for (std::size_t inequality = 0;
       inequality < potential_actual.size(); ++inequality) {
    for (std::size_t coefficient = 0;
         coefficient < farkas_sum.size(); ++coefficient) {
      farkas_sum[coefficient] +=
          farkas_weights[inequality]
          * potential_actual[inequality][coefficient];
    }
  }
  if (
      potential_actual != potential_expected
      || farkas_sum != std::array<Integer, 4>{-1566, 0, 0, 0}
  ) {
    throw std::runtime_error(
        "wall (1,2,1) constant-potential obstruction mismatch");
  }
  std::cout
      << "SU2_WALL_121_CURRENT_NORMAL_FORM"
      << " supports=2..12"
      << " C1_terms=8"
      << " C2_terms=4"
      << " C3_terms=3"
      << " C0_boundary=p_1^3*p_2"
      << " current_suffixes=(g_i,k_i)"
      << " jacobi_identity_supports=2..12"
      << " ratio_band_supports=2..12"
      << " norm_contraction_profile=(0,4,3)"
      << " norm_u=130 norm_w=144 current_pairing=238"
      << " remainder_profile=(0,1,1,1)"
      << " remainder_square=3 remainder_cross=-1"
      << " geometric_certificate_degree=10"
      << " geometric_certificate_coefficients="
      << "(279,1773,5148,9135,11061,9576,"
      << "6057,2835,972,225,27)"
      << " cutoff_profile=(0,8,8,8,6,4,2,1)"
      << " cutoff=2 cutoff_suffix=-230 cutoff_full=2842"
      << " constant_potential_core=4"
      << " farkas_weights=(90,1044,10,603)"
      << " farkas_sum=(-1566,0,0,0)"
      << " result=PASS_EXACT"
      << '\n';
  return EXIT_SUCCESS;
}

int replay_wall_121_saturated_recurrence_obstruction() {
  constexpr int cutoff = 3;
  const std::array<Rational, 8> p{
      Rational(0), Rational(1), Rational(1, 2), Rational(1, 8),
      Rational(1, 32), Rational(1, 128), Rational(1, 512),
      Rational(0)};
  for (int index = 1; index <= 5; ++index) {
    if (p[static_cast<std::size_t>(index)]
            * p[static_cast<std::size_t>(index)]
        < p[static_cast<std::size_t>(index - 1)]
              * p[static_cast<std::size_t>(index + 1)]) {
      throw std::runtime_error(
          "wall (1,2,1) saturated recurrence profile is not log concave");
    }
  }
  if (p[4] * p[6] != p[5] * p[5]) {
    throw std::runtime_error(
        "wall (1,2,1) outer saturation identity mismatch");
  }

  std::array<Rational, 6> g{};
  std::array<Rational, 7> h{};
  std::array<Rational, 6> k{};
  for (int index = 1; index <= cutoff + 1; ++index) {
    g[static_cast<std::size_t>(index)] =
        p[static_cast<std::size_t>(index)]
            * p[static_cast<std::size_t>(index)]
        - p[static_cast<std::size_t>(index - 1)]
              * p[static_cast<std::size_t>(index + 1)];
  }
  for (int index = 1; index <= cutoff + 2; ++index) {
    h[static_cast<std::size_t>(index)] =
        p[static_cast<std::size_t>(index)]
            * p[static_cast<std::size_t>(index - 1)]
        - (index >= 2
               ? p[static_cast<std::size_t>(index - 2)]
                     * p[static_cast<std::size_t>(index + 1)]
               : Rational(0));
  }
  for (int index = 1; index <= cutoff + 1; ++index) {
    k[static_cast<std::size_t>(index)] =
        g[static_cast<std::size_t>(index)]
        + h[static_cast<std::size_t>(index)]
        + h[static_cast<std::size_t>(index + 1)];
  }
  if (g != std::array<Rational, 6>{
               Rational(0), Rational(1), Rational(1, 8),
               Rational(0), Rational(0), Rational(0)}
      || k != std::array<Rational, 6>{
               Rational(0), Rational(3, 2), Rational(21, 32),
               Rational(1, 32), Rational(0), Rational(0)}) {
    throw std::runtime_error(
        "wall (1,2,1) saturated recurrence current mismatch");
  }

  std::array<Rational, 4> q{};
  Rational partial = p[1] * p[1] * p[1] * p[2];
  for (int index = 1; index <= cutoff; ++index) {
    partial +=
        (g[static_cast<std::size_t>(index)]
         - g[static_cast<std::size_t>(index + 1)])
        * (k[static_cast<std::size_t>(index)]
           - k[static_cast<std::size_t>(index + 1)]);
    q[static_cast<std::size_t>(index)] =
        partial
        + g[static_cast<std::size_t>(index)]
              * k[static_cast<std::size_t>(index + 1)]
        + g[static_cast<std::size_t>(index + 1)]
              * k[static_cast<std::size_t>(index)];
  }
  const std::array<Rational, 4> expected{
      Rational(0), Rational(533, 256), Rational(169, 128),
      Rational(337, 256)};
  if (q != expected
      || q[3] - q[2] != Rational(-1, 256)
      || q[3] - std::min(q[1], q[2]) != Rational(-1, 256)) {
    throw std::runtime_error(
        "wall (1,2,1) saturated recurrence obstruction mismatch");
  }

  std::cout
      << "SU2_WALL_121_SATURATED_RECURRENCE_OBSTRUCTION"
      << " profile=(0,1,1/2,1/8,1/32,1/128,1/512)"
      << " cutoff=3"
      << " saturation=p_4*p_6=p_5^2"
      << " currents_g=(1,1/8,0,0)"
      << " currents_k=(3/2,21/32,1/32,0)"
      << " Q=(533/256,169/128,337/256)"
      << " one_step_gap=-1/256"
      << " two_step_min_gap=-1/256"
      << " result=PASS_EXACT"
      << '\n';
  return EXIT_SUCCESS;
}

int replay_support_three_gap_prefix_elevation_obstruction() {
  constexpr int support = 3;
  const Polynomial target = direct_polynomial(
      gap_prefix_polynomial(support, 1, 2, 1), support);
  const std::vector<RationalPolynomial> base =
      outer_bernstein_coefficients(
          root_ratio_parameterization(target, support),
          support);
  const std::array<std::pair<int, Rational>, 5> expected{
      std::pair<int, Rational>{8, Rational(2)},
      std::pair<int, Rational>{9, Rational(1, 4)},
      std::pair<int, Rational>{10, Rational(-1, 3)},
      std::pair<int, Rational>{11, Rational(1, 4)},
      std::pair<int, Rational>{12, Rational(2)}
  };
  if (base.size() != expected.size()) {
    throw std::runtime_error(
        "gap-prefix elevation obstruction degree mismatch");
  }
  for (std::size_t index = 0; index < base.size(); ++index) {
    if (leading_at_last_one(base[index]) != expected[index]) {
      throw std::runtime_error(
          "gap-prefix elevation obstruction leading-term mismatch");
    }
  }
  for (int degree = 4; degree <= 40; ++degree) {
    const std::vector<RationalPolynomial> elevated =
        elevate_bernstein_coefficients(base, degree);
    const std::pair<int, Rational> leading =
        leading_at_last_one(elevated[2]);
    const Rational expected_coefficient(
        -4, degree * (degree - 1));
    if (
        leading.first != 10
        || leading.second != expected_coefficient
    ) {
      throw std::runtime_error(
          "gap-prefix elevated leading-term mismatch");
    }
  }
  const std::vector<RationalPolynomial> degree_twelve =
      elevate_bernstein_coefficients(base, 12);
  const Rational value = evaluate_polynomial(
      degree_twelve[2], {Rational(10), Rational(1)});
  if (value != Rational(-268504190, 11)) {
    throw std::runtime_error(
        "gap-prefix elevation obstruction evaluation mismatch");
  }
  std::cout
      << "SU2_GAP_PREFIX_ELEVATION_OBSTRUCTION"
      << " support=3 R=1 S=2 D=1"
      << " old_root=(1,b,b^2)"
      << " base_leading=(2b^8,b^9/4,-b^10/3,b^11/4,2b^12)"
      << " elevated_index=2"
      << " elevated_leading=-4b^10/(N(N-1))"
      << " N12_b10=-268504190/11"
      << " result=PASS_EXACT\n";
  return EXIT_SUCCESS;
}

int replay_root_outer_bernstein_obstruction() {
  const Polynomial target = substitute_differences(
      radial_polynomial(2, 1, 0), 2);
  const std::vector<RationalPolynomial> bernstein =
      outer_bernstein_coefficients(
          root_ratio_parameterization(target, 2), 2);
  const Polynomial previous_target = substitute_differences(
      radial_polynomial(1, 1, 0), 1);
  const RationalPolynomial previous = rational_polynomial(
      root_ratio_parameterization(previous_target, 1));
  const RationalPolynomial expected{
      {Exponent{1}, Rational(2)},
      {Exponent{2}, Rational(3)},
      {Exponent{3}, Rational(1)},
      {Exponent{4}, Rational(-1, 2)},
      {Exponent{5}, Rational(1, 3)},
      {Exponent{6}, Rational(-1, 6)}};
  if (
      bernstein.size() != 5U
      || !equal_polynomials(bernstein.front(), previous)
      || !equal_polynomials(bernstein[2], expected)
  ) {
    throw std::runtime_error(
        "outer Bernstein obstruction replay mismatch");
  }
  const Rational value =
      evaluate_polynomial(bernstein[2], {Rational(10)});
  if (value != Rational(-411040, 3)) {
    throw std::runtime_error(
        "outer Bernstein obstruction evaluation mismatch");
  }
  std::cout
      << "SU2_AUTOCORRELATION_ROOT_OUTER_BERNSTEIN_OBSTRUCTION"
      << " support=2 A=1 L=0"
      << " base_matches_support_one=1"
      << " outer_index=2"
      << " polynomial=2b+3b^2+b^3-b^4/2+b^5/3-b^6/6"
      << " B2_at_b10=-411040/3"
      << " result=PASS_EXACT\n";
  return EXIT_SUCCESS;
}

int replay_root_outer_first_variation() {
  const Polynomial target = substitute_differences(
      radial_polynomial(3, 2, 0), 3);
  const std::vector<RationalPolynomial> powers =
      outer_power_coefficients(
          root_ratio_parameterization(target, 3), 3);
  const std::vector<Rational> point{Rational(2), Rational(1)};
  const Rational base = evaluate_polynomial(powers[0], point);
  const Rational parameter_first =
      evaluate_polynomial(powers[1], point);
  const Rational outer_scale = 8;
  const Rational first_variation =
      parameter_first / outer_scale;
  if (
      base != 652 || parameter_first != -144
      || first_variation != -18
  ) {
    throw std::runtime_error(
        "outer first-variation replay mismatch");
  }
  std::cout
      << "SU2_AUTOCORRELATION_ROOT_OUTER_FIRST_VARIATION"
      << " root=[1,2,4]"
      << " A=2 L=0"
      << " base=652"
      << " outer_scale=8"
      << " parameter_first_variation=-144"
      << " Q1=-18"
      << " result=PASS_EXACT\n";
  return EXIT_SUCCESS;
}

std::vector<Exponent> grid_indices(const std::vector<int>& degrees) {
  std::vector<Exponent> result{Exponent(degrees.size(), 0)};
  for (std::size_t variable = 0; variable < degrees.size(); ++variable) {
    std::vector<Exponent> expanded;
    expanded.reserve(
        result.size()
        * static_cast<std::size_t>(degrees[variable] + 1));
    for (const Exponent& prefix : result) {
      for (int value = 0; value <= degrees[variable]; ++value) {
        Exponent index = prefix;
        index[variable] = value;
        expanded.push_back(std::move(index));
      }
    }
    result = std::move(expanded);
  }
  return result;
}

std::size_t nd_grid_index(const std::vector<int>& degrees,
                          const Exponent& index) {
  std::size_t result = 0U;
  for (std::size_t variable = 0; variable < degrees.size(); ++variable) {
    result =
        result * static_cast<std::size_t>(degrees[variable] + 1)
        + static_cast<std::size_t>(index[variable]);
  }
  return result;
}

struct NDBernsteinGrid {
  std::vector<int> degrees;
  std::vector<Rational> values;
};

NDBernsteinGrid nd_bernstein_grid(const BigPolynomial& input,
                                  const int variables) {
  NDBernsteinGrid result;
  result.degrees.assign(static_cast<std::size_t>(variables), 0);
  for (const auto& [exponent, coefficient] : input) {
    static_cast<void>(coefficient);
    for (int variable = 0; variable < variables; ++variable) {
      result.degrees[static_cast<std::size_t>(variable)] = std::max(
          result.degrees[static_cast<std::size_t>(variable)],
          exponent[static_cast<std::size_t>(variable)]);
    }
  }
  const std::vector<Exponent> indices = grid_indices(result.degrees);
  result.values.assign(indices.size(), Rational(0));
  for (const auto& [exponent, coefficient] : input) {
    result.values[nd_grid_index(result.degrees, exponent)] = coefficient;
  }

  for (int variable = 0; variable < variables; ++variable) {
    std::vector<Rational> transformed(result.values.size(), Rational(0));
    for (const Exponent& base : indices) {
      if (base[static_cast<std::size_t>(variable)] != 0) {
        continue;
      }
      const int degree =
          result.degrees[static_cast<std::size_t>(variable)];
      for (int output = 0; output <= degree; ++output) {
        Rational value = 0;
        for (int power = 0; power <= output; ++power) {
          Exponent source = base;
          source[static_cast<std::size_t>(variable)] = power;
          value +=
              result.values[nd_grid_index(result.degrees, source)]
              * Rational(
                    binomial(output, power),
                    binomial(degree, power));
        }
        Exponent destination = base;
        destination[static_cast<std::size_t>(variable)] = output;
        transformed[nd_grid_index(result.degrees, destination)] = value;
      }
    }
    result.values = std::move(transformed);
  }
  return result;
}

std::pair<NDBernsteinGrid, NDBernsteinGrid> nd_split_grid(
    const NDBernsteinGrid& input, const int dimension) {
  NDBernsteinGrid left{input.degrees, input.values};
  NDBernsteinGrid right{input.degrees, input.values};
  const std::vector<Exponent> indices = grid_indices(input.degrees);
  const int degree =
      input.degrees[static_cast<std::size_t>(dimension)];
  for (const Exponent& base : indices) {
    if (base[static_cast<std::size_t>(dimension)] != 0) {
      continue;
    }
    std::vector<std::vector<Rational>> triangle(
        static_cast<std::size_t>(degree + 1));
    triangle[0].resize(static_cast<std::size_t>(degree + 1));
    for (int index = 0; index <= degree; ++index) {
      Exponent source = base;
      source[static_cast<std::size_t>(dimension)] = index;
      triangle[0][static_cast<std::size_t>(index)] =
          input.values[nd_grid_index(input.degrees, source)];
    }
    for (int level = 1; level <= degree; ++level) {
      triangle[static_cast<std::size_t>(level)].resize(
          static_cast<std::size_t>(degree - level + 1));
      for (int index = 0; index <= degree - level; ++index) {
        triangle[static_cast<std::size_t>(level)]
                [static_cast<std::size_t>(index)] =
            (
                triangle[static_cast<std::size_t>(level - 1)]
                        [static_cast<std::size_t>(index)]
                + triangle[static_cast<std::size_t>(level - 1)]
                          [static_cast<std::size_t>(index + 1)]
            ) / 2;
      }
    }
    for (int index = 0; index <= degree; ++index) {
      Exponent destination = base;
      destination[static_cast<std::size_t>(dimension)] = index;
      left.values[nd_grid_index(left.degrees, destination)] =
          triangle[static_cast<std::size_t>(index)][0];
      right.values[nd_grid_index(right.degrees, destination)] =
          triangle[static_cast<std::size_t>(degree - index)]
                  [static_cast<std::size_t>(index)];
    }
  }
  return {std::move(left), std::move(right)};
}

struct NDSubdivisionResult {
  std::size_t nodes = 0U;
  std::size_t leaves = 0U;
  std::size_t unresolved = 0U;
  int maximum_depth = 0;
};

void nd_certify_subdivision(const NDBernsteinGrid& grid, const int depth,
                            const int depth_limit,
                            NDSubdivisionResult& result) {
  ++result.nodes;
  result.maximum_depth = std::max(result.maximum_depth, depth);
  if (std::all_of(
          grid.values.begin(), grid.values.end(),
          [](const Rational& value) { return value >= 0; })) {
    ++result.leaves;
    return;
  }
  if (depth >= depth_limit) {
    ++result.unresolved;
    return;
  }
  int dimension = depth % static_cast<int>(grid.degrees.size());
  for (std::size_t offset = 0; offset < grid.degrees.size(); ++offset) {
    const int candidate =
        (dimension + static_cast<int>(offset))
        % static_cast<int>(grid.degrees.size());
    if (grid.degrees[static_cast<std::size_t>(candidate)] > 0) {
      dimension = candidate;
      break;
    }
  }
  auto [left, right] = nd_split_grid(grid, dimension);
  nd_certify_subdivision(left, depth + 1, depth_limit, result);
  nd_certify_subdivision(right, depth + 1, depth_limit, result);
}

struct NDCornerSubdivisionResult {
  NDSubdivisionResult counts;
  std::size_t corner_leaves = 0U;
};

void nd_certify_subdivision_with_corner(
    const NDBernsteinGrid& grid, const int depth,
    const int depth_limit, const Rational& corner_upper,
    const std::vector<Rational>& lower,
    const std::vector<Rational>& upper,
    NDCornerSubdivisionResult& result) {
  ++result.counts.nodes;
  result.counts.maximum_depth =
      std::max(result.counts.maximum_depth, depth);
  if (std::all_of(
          grid.values.begin(), grid.values.end(),
          [](const Rational& value) { return value >= 0; })) {
    ++result.counts.leaves;
    return;
  }
  if (std::all_of(
          upper.begin(), upper.end(),
          [&corner_upper](const Rational& value) {
            return value <= corner_upper;
          })) {
    ++result.corner_leaves;
    return;
  }
  if (depth >= depth_limit) {
    ++result.counts.unresolved;
    return;
  }
  const int dimension =
      depth % static_cast<int>(grid.degrees.size());
  auto [left, right] = nd_split_grid(grid, dimension);
  std::vector<Rational> left_upper = upper;
  std::vector<Rational> right_lower = lower;
  const Rational middle =
      (lower[static_cast<std::size_t>(dimension)]
       + upper[static_cast<std::size_t>(dimension)]) / 2;
  left_upper[static_cast<std::size_t>(dimension)] = middle;
  right_lower[static_cast<std::size_t>(dimension)] = middle;
  nd_certify_subdivision_with_corner(
      left, depth + 1, depth_limit, corner_upper,
      lower, left_upper, result);
  nd_certify_subdivision_with_corner(
      right, depth + 1, depth_limit, corner_upper,
      right_lower, upper, result);
}

int replay_wall_121_q2_floor() {
  constexpr int variables = 3;
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial r = nd_variable(variables, 0);
  const BigPolynomial u = nd_variable(variables, 1);
  const BigPolynomial v = nd_variable(variables, 2);
  BigPolynomial one_minus_u = one;
  big_add_scaled(one_minus_u, u, -1);
  BigPolynomial one_minus_v = one;
  big_add_scaled(one_minus_v, v, -1);
  const BigPolynomial u2v =
      big_multiply(nd_power(u, 2), v);
  BigPolynomial one_minus_2u_plus_u2v = one;
  big_add_scaled(one_minus_2u_plus_u2v, u, -2);
  big_add_scaled(one_minus_2u_plus_u2v, u2v, 1);
  const BigPolynomial uv = big_multiply(u, v);
  BigPolynomial one_minus_uv = one;
  big_add_scaled(one_minus_uv, uv, -1);

  BigPolynomial reduced = one;
  big_add_scaled(
      reduced,
      big_multiply(r, one_minus_u),
      -2);
  big_add_scaled(
      reduced,
      big_multiply(
          nd_power(r, 2),
          one_minus_2u_plus_u2v),
      1);
  big_add_scaled(
      reduced,
      big_multiply(
          nd_power(r, 3),
          nd_power(one_minus_u, 2)),
      2);
  big_add_scaled(
      reduced,
      big_multiply(
          big_multiply(
              big_multiply(nd_power(r, 4), u),
              one_minus_u),
          one_minus_uv),
      2);
  big_add_scaled(
      reduced,
      big_multiply(
          big_multiply(
              big_multiply(nd_power(r, 6), nd_power(u, 3)),
              one_minus_v),
          one_minus_uv),
      1);
  big_add_scaled(
      reduced,
      big_multiply(
          big_multiply(
              nd_power(r, 7), nd_power(u, 4)),
          nd_power(one_minus_v, 2)),
      1);
  big_add_scaled(
      reduced,
      big_multiply(
          big_multiply(
              big_multiply(
                  nd_power(r, 8), nd_power(u, 5)),
              v),
          nd_power(one_minus_v, 2)),
      1);

  const BigPolynomial s = big_multiply(r, u);
  const BigPolynomial t = big_multiply(s, v);
  const BigPolynomial g1 = one;
  const BigPolynomial g2 =
      big_multiply(r, [&]() {
        BigPolynomial value = r;
        big_add_scaled(value, s, -1);
        return value;
      }());
  const BigPolynomial g3 =
      big_multiply(
          big_multiply(big_multiply(r, r), s), [&]() {
        BigPolynomial value = s;
        big_add_scaled(value, t, -1);
        return value;
      }());
  const BigPolynomial h2 = r;
  const BigPolynomial h3 =
      big_multiply(big_multiply(r, s), [&]() {
        BigPolynomial value = r;
        big_add_scaled(value, t, -1);
        return value;
      }());
  const BigPolynomial h4 =
      big_multiply(
          big_multiply(
              big_multiply(big_multiply(r, r), s), t),
          [&]() {
            BigPolynomial value = s;
            big_add_scaled(value, t, -1);
            return value;
          }());
  BigPolynomial k1 = one;
  big_add_scaled(k1, r, 1);
  BigPolynomial k2 = g2;
  big_add_scaled(k2, h2, 1);
  big_add_scaled(k2, h3, 1);
  BigPolynomial k3 = g3;
  big_add_scaled(k3, h3, 1);
  big_add_scaled(k3, h4, 1);
  BigPolynomial direct = r;
  BigPolynomial g1_minus_g2 = g1;
  big_add_scaled(g1_minus_g2, g2, -1);
  BigPolynomial k1_minus_k2 = k1;
  big_add_scaled(k1_minus_k2, k2, -1);
  big_add_scaled(
      direct,
      big_multiply(g1_minus_g2, k1_minus_k2),
      1);
  big_add_scaled(direct, big_multiply(g2, k2), 1);
  big_add_scaled(direct, big_multiply(g3, k3), 1);
  big_add_scaled(direct, one, -1);
  const auto canonicalize = [](
      const BigPolynomial& polynomial) {
    BigPolynomial result;
    for (const auto& [exponent, coefficient] : polynomial) {
      if (coefficient != 0) {
        result.emplace(exponent, coefficient);
      }
    }
    return result;
  };
  if (
      canonicalize(direct)
      != canonicalize(big_multiply(r, reduced))
  ) {
    throw std::runtime_error(
        "wall (1,2,1) Q2 reduced polynomial mismatch");
  }
  const Rational witness_reduced = evaluate_polynomial(
      rational_polynomial(reduced),
      {Rational(2), Rational(1), Rational(255, 256)});
  const Rational witness_q2 = 1 + 2 * witness_reduced;
  const Rational witness_floor = 3;
  const Rational witness_deficit = witness_floor - witness_q2;
  if (
      witness_reduced != Rational(64959, 65536)
      || witness_q2 != Rational(97727, 32768)
      || witness_deficit != Rational(577, 32768)
  ) {
    throw std::runtime_error(
        "wall (1,2,1) Q2 separated-floor witness mismatch");
  }

  const BigPolynomial c = big_multiply(r, s);
  const BigPolynomial d = big_multiply(c, t);
  BigPolynomial one_plus_r = one;
  big_add_scaled(one_plus_r, r, 1);
  BigPolynomial critical_a = one;
  big_add_scaled(critical_a, nd_power(r, 2), 1);
  big_add_scaled(
      critical_a, big_multiply(one_plus_r, c), 1);
  BigPolynomial critical_c = one_plus_r;
  big_add_scaled(critical_c, c, 1);
  BigPolynomial critical_b;
  big_add_scaled(critical_b, r, 3);
  big_add_scaled(critical_b, nd_power(r, 2), 2);
  big_add_scaled(critical_b, big_multiply(r, c), 1);
  big_add_scaled(critical_b, big_multiply(r, d), 1);
  big_add_scaled(
      critical_b, big_multiply(big_multiply(r, c), d), 1);
  big_add_scaled(critical_b, nd_power(r, 3), -1);
  big_add_scaled(critical_b, nd_power(c, 2), -1);
  big_add_scaled(critical_b, nd_power(c, 3), -1);
  BigPolynomial ratio_drop_b;
  big_add_scaled(ratio_drop_b, r, 3);
  big_add_scaled(ratio_drop_b, nd_power(r, 2), 2);
  big_add_scaled(
      ratio_drop_b,
      big_multiply(nd_power(r, 3), one_minus_u),
      -1);
  BigPolynomial second_drop_coefficient =
      big_multiply(nd_power(r, 4), nd_power(u, 2));
  big_add_scaled(
      second_drop_coefficient,
      big_multiply(nd_power(r, 6), nd_power(u, 3)),
      1);
  big_add_scaled(
      ratio_drop_b,
      big_multiply(second_drop_coefficient, one_minus_v),
      -1);
  if (canonicalize(critical_b) != canonicalize(ratio_drop_b)) {
    throw std::runtime_error(
        "wall (1,2,1) coupled branch identity mismatch");
  }
  BigPolynomial critical_delta = nd_power(critical_a, 2);
  big_add_scaled(
      critical_delta,
      big_multiply(critical_b, critical_c),
      3);
  BigPolynomial q2_payment = one;
  big_add_scaled(q2_payment, big_multiply(r, reduced), 1);
  BigPolynomial critical_n =
      big_multiply(critical_a, critical_b);
  big_add_scaled(
      critical_n,
      big_multiply(critical_c, q2_payment),
      9);
  BigPolynomial coupled_certificate =
      big_multiply(critical_c, nd_power(critical_n, 2));
  for (auto& [exponent, coefficient] : coupled_certificate) {
    static_cast<void>(exponent);
    coefficient *= 3;
  }
  big_add_scaled(
      coupled_certificate,
      big_multiply(
          big_multiply(critical_a, critical_n),
          critical_delta),
      4);
  big_add_scaled(
      coupled_certificate,
      big_multiply(critical_b, nd_power(critical_delta, 2)),
      -4);
  const BigPolynomial compact_coupled = canonicalize(
      compactify_positive_variable(coupled_certificate, 0U));
  const NDBernsteinGrid coupled_grid =
      nd_bernstein_grid(compact_coupled, variables);
  const std::size_t coupled_negative = static_cast<std::size_t>(
      std::count_if(
          coupled_grid.values.begin(), coupled_grid.values.end(),
          [](const Rational& value) { return value < 0; }));

  auto [coupled_lower_half, coupled_upper_half] =
      nd_split_grid(coupled_grid, 0);
  auto [coupled_middle_quarter, coupled_top_quarter] =
      nd_split_grid(coupled_upper_half, 0);
  static_cast<void>(coupled_top_quarter);
  NDSubdivisionResult moderate_coupled_subdivision;
  nd_certify_subdivision(
      coupled_lower_half, 0, 24, moderate_coupled_subdivision);
  nd_certify_subdivision(
      coupled_middle_quarter, 0, 24,
      moderate_coupled_subdivision);
  std::cerr << "progress wall_121_q2 moderate_complete\n" << std::flush;

  const BigPolynomial blowup_t = nd_variable(variables, 0);
  const BigPolynomial blowup_alpha = nd_variable(variables, 1);
  const BigPolynomial blowup_beta = nd_variable(variables, 2);
  BigPolynomial blowup_w = one;
  big_add_scaled(blowup_w, blowup_t, -1);
  BigPolynomial blowup_u_numerator = nd_power(blowup_w, 2);
  big_add_scaled(
      blowup_u_numerator,
      big_multiply(
          big_multiply(blowup_t, blowup_w), blowup_alpha),
      -2);
  big_add_scaled(
      blowup_u_numerator,
      big_multiply(nd_power(blowup_t, 2), blowup_alpha),
      -3);
  BigPolynomial blowup_t2_plus_u = nd_power(blowup_t, 2);
  big_add_scaled(blowup_t2_plus_u, blowup_u_numerator, 1);
  const BigPolynomial blowup_denominator = big_multiply(
      nd_power(blowup_u_numerator, 2), blowup_t2_plus_u);
  BigPolynomial two_w_plus_three_t = blowup_w;
  for (auto& [exponent, coefficient] : two_w_plus_three_t) {
    static_cast<void>(exponent);
    coefficient *= 2;
  }
  big_add_scaled(two_w_plus_three_t, blowup_t, 3);
  BigPolynomial one_minus_alpha = one;
  big_add_scaled(one_minus_alpha, blowup_alpha, -1);
  const BigPolynomial blowup_drop_numerator = big_multiply(
      big_multiply(
          big_multiply(
              big_multiply(
                  blowup_beta, nd_power(blowup_t, 4)),
              blowup_w),
          two_w_plus_three_t),
      one_minus_alpha);
  BigPolynomial blowup_v_numerator = blowup_denominator;
  big_add_scaled(
      blowup_v_numerator, blowup_drop_numerator, -1);
  const int coupled_u_degree = coupled_grid.degrees[1];
  const int coupled_v_degree = coupled_grid.degrees[2];
  BigPolynomial blowup_polynomial;
  for (const auto& [exponent, coefficient] : compact_coupled) {
    BigPolynomial term = nd_constant(variables, coefficient);
    term = big_multiply(
        term,
        nd_power(
            blowup_w,
            exponent[0]
                + 2 * (coupled_u_degree - exponent[1])));
    term = big_multiply(
        term, nd_power(blowup_u_numerator, exponent[1]));
    term = big_multiply(
        term,
        nd_power(
            blowup_denominator,
            coupled_v_degree - exponent[2]));
    term = big_multiply(
        term, nd_power(blowup_v_numerator, exponent[2]));
    big_add_scaled(blowup_polynomial, term, 1);
  }
  blowup_polynomial = canonicalize(blowup_polynomial);
  int blowup_t_degree = 0;
  for (const auto& [exponent, coefficient] : blowup_polynomial) {
    static_cast<void>(coefficient);
    blowup_t_degree = std::max(blowup_t_degree, exponent[0]);
  }
  BigPolynomial scaled_blowup;
  for (const auto& [exponent, coefficient] : blowup_polynomial) {
    Integer scale = 1;
    for (int power = exponent[0]; power < blowup_t_degree; ++power) {
      scale *= 4;
    }
    scaled_blowup[exponent] += coefficient * scale;
  }
  scaled_blowup = canonicalize(scaled_blowup);
  constexpr int blowup_t_factor = 10;
  BigPolynomial blowup_quotient;
  for (const auto& [exponent, coefficient] : scaled_blowup) {
    if (exponent[0] < blowup_t_factor) {
      throw std::runtime_error(
          "wall (1,2,1) coupled blowup t-factor mismatch");
    }
    Exponent quotient_exponent = exponent;
    quotient_exponent[0] -= blowup_t_factor;
    blowup_quotient[quotient_exponent] += coefficient;
  }
  blowup_quotient = canonicalize(blowup_quotient);
  std::cerr << "progress wall_121_q2 blowup_constructed\n"
            << std::flush;
  const Exponent leading_alpha2{0, 2, 0};
  const Exponent leading_alpha_beta{1, 1, 1};
  const Exponent leading_beta{2, 0, 1};
  const Exponent leading_alpha{2, 1, 0};
  const Exponent leading_constant{4, 0, 0};
  const Integer leading_scale = blowup_quotient.at(leading_beta);
  if (
      leading_scale <= 0
      || blowup_quotient.at(leading_alpha2) != 64 * leading_scale
      || blowup_quotient.at(leading_alpha_beta) != 8 * leading_scale
      || blowup_quotient.at(leading_alpha) != -leading_scale
      || 128 * blowup_quotient.at(leading_constant)
          != 3 * leading_scale
  ) {
    throw std::runtime_error(
        "wall (1,2,1) coupled blowup Newton face mismatch");
  }
  struct CornerBudget {
    Rational alpha_reserve{0};
    Rational cross_reserve{0};
    Rational beta_reserve{0};
    Rational linear_alpha{0};
    Rational constant_reserve{0};
    Rational discriminant_margin{0};
    std::size_t negative_monomials = 0U;
    std::size_t unclassified_negative = 0U;
  };
  const auto compute_corner_budget = [
      &blowup_quotient, &leading_scale
  ](const Rational& upper) {
    const auto corner_power = [&upper](const int exponent) {
      Rational result(1);
      for (int power = 0; power < exponent; ++power) {
        result *= upper;
      }
      return result;
    };
    Rational alpha_loss(0);
    Rational cross_loss(0);
    Rational beta_loss(0);
    Rational linear_alpha(0);
    Rational constant_loss(0);
    CornerBudget budget;
    for (const auto& [exponent, coefficient] : blowup_quotient) {
      if (coefficient >= 0) {
        continue;
      }
      ++budget.negative_monomials;
      const int t_exponent = exponent[0];
      const int alpha_exponent = exponent[1];
      const int beta_exponent = exponent[2];
      const Rational magnitude(-coefficient, leading_scale);
      if (
          alpha_exponent == 1 && beta_exponent == 0
          && t_exponent >= 2
      ) {
        linear_alpha += magnitude * corner_power(t_exponent - 2);
      } else if (
          alpha_exponent >= 1 && beta_exponent >= 1
          && t_exponent >= 1
      ) {
        cross_loss += magnitude * corner_power(
            t_exponent - 1 + alpha_exponent - 1
            + beta_exponent - 1);
      } else if (alpha_exponent >= 2) {
        alpha_loss += magnitude * corner_power(
            t_exponent + alpha_exponent - 2 + beta_exponent);
      } else if (beta_exponent >= 1 && t_exponent >= 2) {
        beta_loss += magnitude * corner_power(
            t_exponent - 2 + alpha_exponent + beta_exponent - 1);
      } else if (
          alpha_exponent == 0 && beta_exponent == 0
          && t_exponent >= 4
      ) {
        constant_loss += magnitude * corner_power(t_exponent - 4);
      } else {
        ++budget.unclassified_negative;
      }
    }
    budget.alpha_reserve = Rational(64) - alpha_loss;
    budget.cross_reserve = Rational(8) - cross_loss;
    budget.beta_reserve = Rational(1) - beta_loss;
    budget.linear_alpha = linear_alpha;
    budget.constant_reserve = Rational(3, 128) - constant_loss;
    budget.discriminant_margin =
        4 * budget.alpha_reserve * budget.constant_reserve
        - budget.linear_alpha * budget.linear_alpha;
    return budget;
  };
  int corner_power_two = 4;
  Rational corner_upper(1, 16);
  CornerBudget corner_budget =
      compute_corner_budget(corner_upper);
  for (; corner_power_two <= 64; ++corner_power_two) {
    Integer denominator = 1;
    denominator <<= corner_power_two;
    corner_upper = Rational(Integer(1), denominator);
    corner_budget = compute_corner_budget(corner_upper);
    if (
        corner_budget.unclassified_negative == 0U
        && corner_budget.alpha_reserve >= 0
        && corner_budget.cross_reserve >= 0
        && corner_budget.beta_reserve >= 0
        && corner_budget.constant_reserve >= 0
        && corner_budget.discriminant_margin >= 0
    ) {
      break;
    }
  }
  if (corner_power_two > 64) {
    throw std::runtime_error(
        "wall (1,2,1) coupled blowup corner budget failed");
  }
  std::cerr << "progress wall_121_q2 corner_budget_complete\n"
            << std::flush;
  const NDBernsteinGrid blowup_grid =
      nd_bernstein_grid(blowup_quotient, variables);
  const std::size_t blowup_negative = static_cast<std::size_t>(
      std::count_if(
          blowup_grid.values.begin(), blowup_grid.values.end(),
          [](const Rational& value) { return value < 0; }));
  NDCornerSubdivisionResult blowup_subdivision;
  const std::vector<Rational> blowup_lower(
      static_cast<std::size_t>(variables), Rational(0));
  const std::vector<Rational> blowup_upper(
      static_cast<std::size_t>(variables), Rational(1));
  const int blowup_depth_limit = 3 * corner_power_two + 18;
  std::cerr << "progress wall_121_q2 blowup_subdivision_start\n"
            << std::flush;
  nd_certify_subdivision_with_corner(
      blowup_grid, 0, blowup_depth_limit, corner_upper,
      blowup_lower, blowup_upper, blowup_subdivision);
  const BigPolynomial compact = canonicalize(
      compactify_positive_variable(reduced, 0U));
  const NDBernsteinGrid grid =
      nd_bernstein_grid(compact, variables);
  auto [lower_half, upper_half] = nd_split_grid(grid, 0);
  static_cast<void>(upper_half);
  NDSubdivisionResult result;
  nd_certify_subdivision(lower_half, 0, 18, result);
  if (
      grid.degrees != std::vector<int>{8, 5, 3}
      || coupled_grid.degrees != std::vector<int>{24, 13, 6}
      || coupled_grid.values.size() != 2450U
      || coupled_negative != 253U
      || moderate_coupled_subdivision.nodes != 14U
      || moderate_coupled_subdivision.leaves != 8U
      || moderate_coupled_subdivision.unresolved != 0U
      || moderate_coupled_subdivision.maximum_depth != 6
      || blowup_grid.degrees != std::vector<int>{76, 27, 6}
      || blowup_quotient.size() != 10861U
      || blowup_grid.values.size() != 15092U
      || blowup_negative != 75U
      || corner_budget.negative_monomials != 5396U
      || corner_budget.unclassified_negative != 0U
      || corner_power_two != 5
      || corner_budget.alpha_reserve < 0
      || corner_budget.cross_reserve < 0
      || corner_budget.beta_reserve < 0
      || corner_budget.constant_reserve < 0
      || corner_budget.discriminant_margin < 0
      || blowup_subdivision.counts.nodes != 611U
      || blowup_subdivision.counts.leaves != 305U
      || blowup_subdivision.corner_leaves != 1U
      || blowup_subdivision.counts.unresolved != 0U
      || blowup_subdivision.counts.maximum_depth != 32
      || lower_half.values.size() != 216U
      || result.nodes != 1U
      || result.leaves != 1U
      || result.unresolved != 0U
      || result.maximum_depth != 0
  ) {
    throw std::runtime_error(
        "wall (1,2,1) Q2 Bernstein certificate mismatch");
  }
  std::cout
      << "SU2_WALL_121_Q2_FLOOR"
      << " degrees=(" << grid.degrees[0]
      << "," << grid.degrees[1]
      << "," << grid.degrees[2] << ")"
      << " bernstein_coefficients=216"
      << " bernstein_leaf=[0,1/2]x[0,1]^2"
      << " analytic_tail=R>=1"
      << " separated_floor_witness_Q2=97727/32768"
      << " separated_floor=3"
      << " separated_floor_deficit=577/32768"
      << " coupled_degrees=(24,13,6)"
      << " coupled_coefficients=2450"
      << " coupled_initial_negative=253"
      << " coupled_branch_identity=1"
      << " coupled_moderate_nodes=14"
      << " coupled_moderate_leaves=8"
      << " coupled_moderate_unresolved=0"
      << " coupled_moderate_depth=6"
      << " coupled_blowup_degrees=(76,27,6)"
      << " coupled_blowup_terms=10861"
      << " coupled_blowup_coefficients=15092"
      << " coupled_blowup_initial_negative=75"
      << " coupled_corner=1/32"
      << " coupled_corner_negative_monomials=5396"
      << " coupled_corner_unclassified=0"
      << " coupled_corner_reserves_nonnegative=1"
      << " coupled_corner_discriminant_nonnegative=1"
      << " coupled_blowup_nodes=611"
      << " coupled_blowup_leaves=305"
      << " coupled_blowup_corner_leaves=1"
      << " coupled_blowup_unresolved=0"
      << " coupled_blowup_depth=32"
      << " coupled_payment=PASS_EXACT"
      << " result=PASS_EXACT"
      << '\n';
  return EXIT_SUCCESS;
}

int replay_wall_121_small_ratio_payment() {
  const Rational r(1, 3);
  const Rational tail_prefactor =
      2 * (1 + r + r * r)
      / (1 - r * r * r * r);
  const Rational current_linear_reserve(2, 9);
  const Rational tail_linear_debit(1, 27);
  const Rational current_margin =
      current_linear_reserve - tail_linear_debit;
  const Rational c2_max(34, 27);
  const Rational c3_max(13, 9);
  const Rational b_over_r_max(4);
  const Rational critical_half_min(7, 4);
  const Rational b_max(4, 3);
  const Rational demand_max =
      (c2_max + c3_max) / 4;
  if (
      tail_prefactor != Rational(117, 40)
      || !(tail_prefactor < 3)
      || current_margin != Rational(5, 27)
      || c2_max != Rational(34, 27)
      || c3_max != Rational(13, 9)
      || b_over_r_max != 4
      || !(critical_half_min > b_max)
      || demand_max != Rational(73, 108)
      || !(demand_max < 1)
  ) {
    throw std::runtime_error(
        "wall (1,2,1) small-ratio payment mismatch");
  }
  std::cout
      << "SU2_WALL_121_SMALL_RATIO_PAYMENT"
      << " ratio_max=1/3"
      << " tail_prefactor=117/40<3"
      << " current_margin=(5/27)R"
      << " C2_max=34/27"
      << " C3_max=13/9"
      << " B_max=(4R)<=4/3"
      << " critical_half_min=7/4"
      << " demand_max=73/108<1"
      << " result=PASS_EXACT\n";
  return EXIT_SUCCESS;
}

bool z3_negative_on_unit_cube(const BigPolynomial& polynomial,
                              std::string& model_text) {
  if (polynomial.empty()) {
    return false;
  }
  const int variables =
      static_cast<int>(polynomial.begin()->first.size());
  z3::context context;
  z3::solver solver(context);
  std::vector<z3::expr> coordinate;
  coordinate.reserve(static_cast<std::size_t>(variables));
  for (int variable = 0; variable < variables; ++variable) {
    coordinate.push_back(
        context.real_const(("x_" + std::to_string(variable)).c_str()));
    solver.add(coordinate.back() >= 0);
    solver.add(coordinate.back() <= 1);
  }
  z3::expr value = context.int_val(0);
  for (const auto& [exponent, coefficient] : polynomial) {
    z3::expr term = context.int_val(coefficient.str().c_str());
    for (int variable = 0; variable < variables; ++variable) {
      for (
          int power = 0;
          power < exponent[static_cast<std::size_t>(variable)];
          ++power
      ) {
        term = term * coordinate[static_cast<std::size_t>(variable)];
      }
    }
    value = value + term;
  }
  solver.add(value < 0);
  const z3::check_result result = solver.check();
  if (result == z3::unknown) {
    throw std::runtime_error(
        "Z3 returned unknown on gap-prefix polynomial");
  }
  if (result == z3::sat) {
    model_text = solver.get_model().to_string();
    return true;
  }
  return false;
}

int analyze_support_two_gap_prefix_certificate() {
  constexpr int support = 2;
  constexpr int maximum_relevant_label = 2 * support;
  constexpr int subdivision_depth = 24;
  std::size_t cases = 0U;
  std::size_t coefficients = 0U;
  std::size_t initial_negative = 0U;
  std::size_t nodes = 0U;
  std::size_t leaves = 0U;
  std::size_t unresolved = 0U;
  std::size_t nlsat_unsat = 0U;
  for (int first_label = 1;
       first_label <= maximum_relevant_label;
       ++first_label) {
    for (int second_label = 1;
         second_label <= maximum_relevant_label;
         ++second_label) {
      if (first_label == second_label) {
        continue;
      }
      const int maximum_gap =
          support + std::min(first_label, second_label);
      for (int gap = 1; gap <= maximum_gap; ++gap) {
        const Polynomial target = direct_polynomial(
            gap_prefix_polynomial(
                support, first_label, second_label, gap),
            support);
        const BigPolynomial root_parameterized =
            root_ratio_parameterization(target, support);
        const BigPolynomial parameterized =
            compactify_first_positive_variable(root_parameterized);
        const NDBernsteinGrid grid =
            nd_bernstein_grid(parameterized, support);
        NDSubdivisionResult subdivision;
        nd_certify_subdivision(
            grid, 0, subdivision_depth, subdivision);
        ++cases;
        coefficients += grid.values.size();
        initial_negative += static_cast<std::size_t>(std::count_if(
            grid.values.begin(), grid.values.end(),
            [](const Rational& value) { return value < 0; }));
        nodes += subdivision.nodes;
        leaves += subdivision.leaves;
        unresolved += subdivision.unresolved;
        if (subdivision.unresolved != 0U) {
          std::string model_text;
          if (z3_negative_on_unit_cube(parameterized, model_text)) {
            throw std::runtime_error(
                "support-two gap-prefix counterexample at labels "
                + std::to_string(first_label) + ","
                + std::to_string(second_label)
                + " and gap " + std::to_string(gap)
                + ": " + model_text);
          }
          ++nlsat_unsat;
        }
      }
    }
  }
  if (
      cases != 44U
      || coefficients != 1980U
      || initial_negative != 58U
      || nodes != 768U
      || leaves != 392U
      || unresolved != 14U
      || nlsat_unsat != 14U
  ) {
    throw std::runtime_error(
        "support-two gap-prefix certificate count mismatch");
  }
  std::cout
      << "SU2_SUPPORT_TWO_GAP_PREFIX_CERTIFICATE"
      << " cases=" << cases
      << " coefficients=" << coefficients
      << " initial_negative=" << initial_negative
      << " subdivision_nodes=" << nodes
      << " subdivision_leaves=" << leaves
      << " subdivision_unresolved=" << unresolved
      << " nlsat_unsat=" << nlsat_unsat
      << " result=PASS_EXACT"
      << '\n';
  return EXIT_SUCCESS;
}

int analyze_support_three_gap_prefix_certificate_part(
    const int selected_first_label) {
  constexpr int support = 3;
  constexpr int maximum_relevant_label = 2 * support;
  constexpr int subdivision_depth = 24;
  if (
      selected_first_label < 1
      || selected_first_label > maximum_relevant_label
  ) {
    throw std::invalid_argument(
        "support-three first label must lie in [1,6]");
  }
  std::size_t cases = 0U;
  std::size_t coefficients = 0U;
  std::size_t initial_negative = 0U;
  std::size_t nodes = 0U;
  std::size_t leaves = 0U;
  std::size_t unresolved = 0U;
  std::size_t nlsat_unsat = 0U;
  for (int second_label = 1;
       second_label <= maximum_relevant_label;
       ++second_label) {
    if (selected_first_label == second_label) {
      continue;
    }
    const int maximum_gap =
        support + std::min(selected_first_label, second_label);
    for (int gap = 1; gap <= maximum_gap; ++gap) {
      const Polynomial target = direct_polynomial(
          gap_prefix_polynomial(
              support, selected_first_label, second_label, gap),
          support);
      BigPolynomial parameterized =
          support_three_ordered_ratio_difference_parameterization(target);
      for (std::size_t variable = 0; variable < 3U; ++variable) {
        parameterized =
            compactify_positive_variable(parameterized, variable);
      }
      const BigPolynomial quotient =
          remove_common_monomial_factor(parameterized);
      const NDBernsteinGrid grid =
          nd_bernstein_grid(quotient, support);
      NDSubdivisionResult subdivision;
      nd_certify_subdivision(
          grid, 0, subdivision_depth, subdivision);
      ++cases;
      coefficients += grid.values.size();
      initial_negative += static_cast<std::size_t>(std::count_if(
          grid.values.begin(), grid.values.end(),
          [](const Rational& value) { return value < 0; }));
      nodes += subdivision.nodes;
      leaves += subdivision.leaves;
      unresolved += subdivision.unresolved;
      if (subdivision.unresolved != 0U) {
        std::string model_text;
        const bool quotient_negative =
            z3_negative_on_unit_cube(quotient, model_text);
        if (
            quotient_negative
            && z3_negative_on_unit_cube(parameterized, model_text)
        ) {
          throw std::runtime_error(
              "support-three gap-prefix counterexample at labels "
              + std::to_string(selected_first_label) + ","
              + std::to_string(second_label)
              + " and gap " + std::to_string(gap)
              + ": " + model_text);
        }
        ++nlsat_unsat;
      }
    }
    std::cout
        << "SU2_SUPPORT_THREE_GAP_PREFIX_PROGRESS"
        << " first_label=" << selected_first_label
        << " second_label=" << second_label
        << " cases=" << cases
        << " unresolved=" << unresolved
        << '\n';
    std::cout.flush();
  }
  const std::array<std::size_t, 6> expected_cases{
      20U, 24U, 27U, 29U, 30U, 30U};
  const std::array<std::size_t, 6> expected_coefficients{
      11520U, 13995U, 15750U, 16920U, 17505U, 17190U};
  const std::array<std::size_t, 6> expected_initial_negative{
      788U, 633U, 531U, 321U, 135U, 0U};
  const std::array<std::size_t, 6> expected_nodes{
      212U, 190U, 151U, 155U, 114U, 30U};
  const std::array<std::size_t, 6> expected_leaves{
      116U, 107U, 89U, 92U, 72U, 30U};
  const std::size_t partition =
      static_cast<std::size_t>(selected_first_label - 1);
  if (
      cases != expected_cases[partition]
      || coefficients != expected_coefficients[partition]
      || initial_negative != expected_initial_negative[partition]
      || nodes != expected_nodes[partition]
      || leaves != expected_leaves[partition]
      || unresolved != 0U
      || nlsat_unsat != 0U
  ) {
    throw std::runtime_error(
        "support-three gap-prefix partition count mismatch");
  }
  std::cout
      << "SU2_SUPPORT_THREE_GAP_PREFIX_CERTIFICATE_PART"
      << " first_label=" << selected_first_label
      << " cases=" << cases
      << " coefficients=" << coefficients
      << " initial_negative=" << initial_negative
      << " subdivision_nodes=" << nodes
      << " subdivision_leaves=" << leaves
      << " subdivision_unresolved=" << unresolved
      << " nlsat_unsat=" << nlsat_unsat
      << " result=PASS_EXACT"
      << '\n';
  return EXIT_SUCCESS;
}

int analyze_support_four_exceptional_gap_prefix_certificate() {
  constexpr int support = 4;
  constexpr int subdivision_depth = 32;
  const Polynomial target = direct_polynomial(
      gap_prefix_polynomial(support, 1, 2, 1), support);
  BigPolynomial parameterized =
      ordered_ratio_difference_parameterization(target, support);
  for (std::size_t variable = 0; variable < 4U; ++variable) {
    parameterized =
        compactify_positive_variable(parameterized, variable);
  }
  const BigPolynomial quotient =
      remove_common_monomial_factor(parameterized);
  const NDBernsteinGrid grid =
      nd_bernstein_grid(quotient, support);
  NDSubdivisionResult subdivision;
  nd_certify_subdivision(
      grid, 0, subdivision_depth, subdivision);
  const std::size_t initial_negative =
      static_cast<std::size_t>(std::count_if(
          grid.values.begin(), grid.values.end(),
          [](const Rational& value) { return value < 0; }));
  if (
      grid.values.size() != 9945U
      || initial_negative != 766U
      || subdivision.nodes != 55U
      || subdivision.leaves != 28U
      || subdivision.unresolved != 0U
      || subdivision.maximum_depth != 7
  ) {
    throw std::runtime_error(
        "support-four exceptional gap-prefix count mismatch");
  }
  std::cout
      << "SU2_SUPPORT_FOUR_EXCEPTIONAL_GAP_PREFIX_CERTIFICATE"
      << " labels=(1,2,1)"
      << " coefficients=" << grid.values.size()
      << " initial_negative=" << initial_negative
      << " subdivision_nodes=" << subdivision.nodes
      << " subdivision_leaves=" << subdivision.leaves
      << " subdivision_unresolved=" << subdivision.unresolved
      << " maximum_depth=" << subdivision.maximum_depth
      << " result=PASS_EXACT"
      << '\n';
  return EXIT_SUCCESS;
}

int analyze_support_four_161_saturated_blowup_certificate() {
  constexpr int variables = 3;
  constexpr int subdivision_depth = 16;
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial bounded_t = nd_variable(variables, 0);
  const BigPolynomial bounded_u = nd_variable(variables, 1);
  const BigPolynomial bounded_v = nd_variable(variables, 2);
  BigPolynomial one_minus_u = one;
  big_add_scaled(one_minus_u, bounded_u, -1);
  BigPolynomial one_minus_v = one;
  big_add_scaled(one_minus_v, bounded_v, -1);
  BigPolynomial one_minus_uv = one;
  big_add_scaled(
      one_minus_uv, big_multiply(bounded_u, bounded_v), -1);

  BigPolynomial reduced;
  big_add_scaled(
      reduced, big_multiply(one_minus_u, one_minus_uv), 1);
  BigPolynomial t_minus_square = bounded_t;
  big_add_scaled(t_minus_square, nd_power(bounded_t, 2), -1);
  big_add_scaled(
      reduced,
      big_multiply(
          big_multiply(one_minus_u, bounded_v),
          t_minus_square),
      1);
  big_add_scaled(
      reduced,
      big_multiply(
          nd_power(bounded_t, 3), nd_power(one_minus_v, 2)),
      1);
  BigPolynomial fourth_bracket = one;
  big_add_scaled(fourth_bracket, bounded_v, -2);
  big_add_scaled(
      fourth_bracket,
      big_multiply(bounded_u, nd_power(bounded_v, 2)),
      1);
  big_add_scaled(
      reduced,
      big_multiply(
          big_multiply(nd_power(bounded_t, 4), bounded_v),
          fourth_bracket),
      1);
  big_add_scaled(
      reduced,
      big_multiply(
          big_multiply(
              nd_power(bounded_t, 5), nd_power(bounded_v, 2)),
          one_minus_v),
      -2);
  big_add_scaled(
      reduced,
      big_multiply(
          nd_power(bounded_t, 6), nd_power(bounded_v, 3)),
      1);
  big_add_scaled(
      reduced,
      big_multiply(
          nd_power(bounded_t, 7), nd_power(bounded_v, 4)),
      2);

  const auto canonicalize = [](
      const BigPolynomial& polynomial) {
    BigPolynomial result;
    for (const auto& [exponent, coefficient] : polynomial) {
      if (coefficient != 0) {
        result.emplace(exponent, coefficient);
      }
    }
    return result;
  };
  const BigPolynomial reduced_canonical = canonicalize(reduced);

  constexpr int support = 4;
  const Polynomial fusion_target = direct_polynomial(
      gap_prefix_polynomial(support, 1, 6, 1), support);
  std::size_t p_zero_terms = 0U;
  bool p_zero_matches = true;
  BigPolynomial normalized_saturated;
  for (const auto& [exponent, coefficient] : fusion_target) {
    if (coefficient == 0) {
      continue;
    }
    p_zero_matches = p_zero_matches && exponent[0] <= 1;
    if (exponent[0] == 1) {
      ++p_zero_terms;
      p_zero_matches =
          p_zero_matches
          && exponent == Exponent({1, 1, 1, 0, 1})
          && coefficient == -1;
    }
    const int t_exponent =
        -exponent[0] + exponent[2]
        + 2 * exponent[3] + 3 * exponent[4] - 5;
    const int u_exponent =
        exponent[0] - exponent[2]
        - exponent[3] - exponent[4] + 4;
    normalized_saturated[
        Exponent{t_exponent, u_exponent, exponent[4]}] +=
        coefficient;
  }
  normalized_saturated = canonicalize(normalized_saturated);
  if (
      !p_zero_matches
      || p_zero_terms != 1U
      || std::any_of(
          normalized_saturated.begin(),
          normalized_saturated.end(),
          [](const auto& term) {
            return std::any_of(
                term.first.begin(), term.first.end(),
                [](const int exponent) { return exponent < 0; });
          })
      || normalized_saturated != reduced_canonical
  ) {
    throw std::runtime_error(
        "support-four (1,6,1) saturated identity mismatch");
  }

  BigPolynomial local =
      reflect_unit_variable(reduced, 1U);
  local = reflect_unit_variable(local, 2U);
  local = canonicalize(local);
  const BigPolynomial local_a = nd_variable(variables, 1);
  const BigPolynomial local_b = nd_variable(variables, 2);
  BigPolynomial local_c = one;
  big_add_scaled(local_c, local_b, -1);
  BigPolynomial first_bracket = local_a;
  big_add_scaled(first_bracket, local_b, 1);
  big_add_scaled(
      first_bracket, big_multiply(local_a, local_b), -1);
  BigPolynomial local_bracket = one;
  big_add_scaled(local_bracket, bounded_t, -1);
  big_add_scaled(
      local_bracket,
      big_multiply(
          nd_power(bounded_t, 3), nd_power(local_c, 2)),
      -1);
  BigPolynomial square_base = local_b;
  big_add_scaled(
      square_base,
      big_multiply(
          nd_power(bounded_t, 2), nd_power(local_c, 2)),
      -1);
  BigPolynomial local_decomposition;
  big_add_scaled(
      local_decomposition,
      big_multiply(local_a, first_bracket),
      1);
  big_add_scaled(
      local_decomposition,
      big_multiply(
          big_multiply(
              big_multiply(local_a, bounded_t), local_c),
          local_bracket),
      1);
  big_add_scaled(
      local_decomposition,
      big_multiply(
          nd_power(bounded_t, 3), nd_power(square_base, 2)),
      1);
  big_add_scaled(
      local_decomposition,
      big_multiply(
          big_multiply(
              nd_power(bounded_t, 4), nd_power(local_b, 2)),
          local_c),
      1);
  big_add_scaled(
      local_decomposition,
      big_multiply(nd_power(bounded_t, 6), nd_power(local_c, 3)),
      1);
  big_add_scaled(
      local_decomposition,
      big_multiply(nd_power(bounded_t, 7), nd_power(local_c, 4)),
      1);
  if (canonicalize(local_decomposition) != local) {
    throw std::runtime_error(
        "support-four (1,6,1) local decomposition mismatch");
  }

  const BigPolynomial s = nd_variable(variables, 0);
  const BigPolynomial w = nd_variable(variables, 1);
  const BigPolynomial u = nd_variable(variables, 2);
  BigPolynomial t = s;
  big_add_scaled(t, w, 1);
  BigPolynomial blowup_one_minus_u = one;
  big_add_scaled(blowup_one_minus_u, u, -1);
  BigPolynomial blowup;
  BigPolynomial term = t;
  big_add_scaled(term, big_multiply(u, w), -1);
  big_add_scaled(
      blowup, big_multiply(blowup_one_minus_u, term), 1);
  term = big_multiply(t, w);
  big_add_scaled(
      term, big_multiply(nd_power(t, 2), w), -1);
  big_add_scaled(
      blowup, big_multiply(blowup_one_minus_u, term), 1);
  big_add_scaled(blowup, nd_power(t, 4), 1);
  big_add_scaled(
      blowup, big_multiply(nd_power(t, 3), w), -2);
  big_add_scaled(
      blowup,
      big_multiply(nd_power(t, 2), nd_power(w, 2)),
      1);
  big_add_scaled(
      blowup, big_multiply(nd_power(t, 4), w), 1);
  big_add_scaled(
      blowup,
      big_multiply(nd_power(t, 3), nd_power(w, 2)),
      -2);
  big_add_scaled(
      blowup,
      big_multiply(
          big_multiply(u, nd_power(t, 2)), nd_power(w, 3)),
      1);
  big_add_scaled(
      blowup,
      big_multiply(nd_power(t, 4), nd_power(w, 2)),
      -2);
  big_add_scaled(
      blowup,
      big_multiply(nd_power(t, 3), nd_power(w, 3)),
      2);
  big_add_scaled(
      blowup,
      big_multiply(nd_power(t, 4), nd_power(w, 3)),
      1);
  big_add_scaled(
      blowup,
      big_multiply(nd_power(t, 4), nd_power(w, 4)),
      2);

  struct Certificate {
    std::size_t coefficients = 0U;
    std::size_t initial_negative = 0U;
    NDSubdivisionResult subdivision;
  };
  const auto certify = [&canonicalize](
      const BigPolynomial& polynomial) {
    const BigPolynomial canonical = canonicalize(polynomial);
    const BigPolynomial quotient =
        remove_common_monomial_factor(canonical);
    const NDBernsteinGrid grid =
        nd_bernstein_grid(quotient, variables);
    Certificate result;
    result.coefficients = grid.values.size();
    result.initial_negative =
        static_cast<std::size_t>(std::count_if(
            grid.values.begin(), grid.values.end(),
            [](const Rational& value) { return value < 0; }));
    nd_certify_subdivision(
        grid, 0, subdivision_depth, result.subdivision);
    return result;
  };

  const Certificate bounded_upper = certify(
      half_shift_positive_variable(local, 0U));
  BigPolynomial s_half =
      half_shift_positive_variable(blowup, 0U);
  s_half = compactify_positive_variable(s_half, 0U);
  s_half = compactify_positive_variable(s_half, 1U);
  const Certificate s_half_result = certify(s_half);
  BigPolynomial w_half =
      half_shift_positive_variable(blowup, 1U);
  w_half = compactify_positive_variable(w_half, 0U);
  w_half = compactify_positive_variable(w_half, 1U);
  const Certificate w_half_result = certify(w_half);

  if (
      bounded_upper.coefficients != 120U
      || bounded_upper.initial_negative != 11U
      || bounded_upper.subdivision.nodes != 13U
      || bounded_upper.subdivision.leaves != 7U
      || bounded_upper.subdivision.unresolved != 0U
      || bounded_upper.subdivision.maximum_depth != 3
      || s_half_result.coefficients != 135U
      || s_half_result.initial_negative != 25U
      || s_half_result.subdivision.nodes != 7U
      || s_half_result.subdivision.leaves != 4U
      || s_half_result.subdivision.unresolved != 0U
      || s_half_result.subdivision.maximum_depth != 2
      || w_half_result.coefficients != 135U
      || w_half_result.initial_negative != 3U
      || w_half_result.subdivision.nodes != 5U
      || w_half_result.subdivision.leaves != 3U
      || w_half_result.subdivision.unresolved != 0U
      || w_half_result.subdivision.maximum_depth != 2
  ) {
    throw std::runtime_error(
        "support-four (1,6,1) chart certificate mismatch");
  }
  std::cout
      << "SU2_SUPPORT_FOUR_161_SATURATED_BLOWUP_CERTIFICATE"
      << " p_zero_coefficient=-p_1*p_2*p_4"
      << " reduced_terms=" << reduced_canonical.size()
      << " local_identity=PASS_EXACT"
      << " bounded_upper_coefficients="
      << bounded_upper.coefficients
      << " bounded_upper_initial_negative="
      << bounded_upper.initial_negative
      << " bounded_upper_nodes="
      << bounded_upper.subdivision.nodes
      << " bounded_upper_leaves="
      << bounded_upper.subdivision.leaves
      << " bounded_upper_unresolved="
      << bounded_upper.subdivision.unresolved
      << " s_half_coefficients=" << s_half_result.coefficients
      << " s_half_initial_negative="
      << s_half_result.initial_negative
      << " s_half_nodes=" << s_half_result.subdivision.nodes
      << " s_half_leaves=" << s_half_result.subdivision.leaves
      << " s_half_unresolved="
      << s_half_result.subdivision.unresolved
      << " w_half_coefficients=" << w_half_result.coefficients
      << " w_half_initial_negative="
      << w_half_result.initial_negative
      << " w_half_nodes=" << w_half_result.subdivision.nodes
      << " w_half_leaves=" << w_half_result.subdivision.leaves
      << " w_half_unresolved="
      << w_half_result.subdivision.unresolved
      << " result=PASS_EXACT"
      << '\n';
  return EXIT_SUCCESS;
}

struct BernsteinResult {
  std::size_t coefficients = 0U;
  std::size_t negative = 0U;
  Rational minimum = 0;
  Exponent minimum_index;
  std::array<int, 3> degrees{0, 0, 0};
};

struct BernsteinGrid {
  std::array<int, 3> degrees{0, 0, 0};
  std::vector<Rational> values;
};

std::size_t grid_index(const std::array<int, 3>& degrees, const int first,
                       const int second, const int third) {
  return (
      static_cast<std::size_t>(first)
          * static_cast<std::size_t>(degrees[1] + 1)
      + static_cast<std::size_t>(second)
  ) * static_cast<std::size_t>(degrees[2] + 1)
      + static_cast<std::size_t>(third);
}

BernsteinGrid bernstein_grid(const BigPolynomial& input) {
  BernsteinGrid result;
  for (const auto& [exponent, coefficient] : input) {
    static_cast<void>(coefficient);
    for (std::size_t variable = 0; variable < 3U; ++variable) {
      result.degrees[variable] =
          std::max(result.degrees[variable], exponent[variable]);
    }
  }
  const std::size_t size =
      static_cast<std::size_t>(result.degrees[0] + 1)
      * static_cast<std::size_t>(result.degrees[1] + 1)
      * static_cast<std::size_t>(result.degrees[2] + 1);
  result.values.assign(size, Rational(0));
  for (int first = 0; first <= result.degrees[0]; ++first) {
    for (int second = 0; second <= result.degrees[1]; ++second) {
      for (int third = 0; third <= result.degrees[2]; ++third) {
        const Exponent index{first, second, third};
        Rational value = 0;
        for (const auto& [exponent, coefficient] : input) {
          bool eligible = true;
          Rational factor = coefficient;
          for (std::size_t variable = 0; variable < 3U; ++variable) {
            if (exponent[variable] > index[variable]) {
              eligible = false;
              break;
            }
            factor *= Rational(
                binomial(index[variable], exponent[variable]),
                binomial(result.degrees[variable], exponent[variable]));
          }
          if (eligible) {
            value += factor;
          }
        }
        result.values[grid_index(
            result.degrees, first, second, third)] = value;
      }
    }
  }
  return result;
}

BernsteinResult bernstein_cube(const BigPolynomial& input) {
  const BernsteinGrid grid = bernstein_grid(input);
  BernsteinResult result;
  result.degrees = grid.degrees;
  bool initialized = false;
  for (int first = 0; first <= grid.degrees[0]; ++first) {
    for (int second = 0; second <= grid.degrees[1]; ++second) {
      for (int third = 0; third <= grid.degrees[2]; ++third) {
        const Exponent index{first, second, third};
        const Rational& value =
            grid.values[grid_index(grid.degrees, first, second, third)];
        ++result.coefficients;
        if (value < 0) {
          ++result.negative;
        }
        if (!initialized || value < result.minimum) {
          initialized = true;
          result.minimum = value;
          result.minimum_index = index;
        }
      }
    }
  }
  return result;
}

std::pair<BernsteinGrid, BernsteinGrid> split_grid(
    const BernsteinGrid& input, const int dimension) {
  BernsteinGrid left{input.degrees, input.values};
  BernsteinGrid right{input.degrees, input.values};
  const int degree = input.degrees[static_cast<std::size_t>(dimension)];
  const int first_max = dimension == 0 ? 0 : input.degrees[0];
  const int second_max = dimension == 1 ? 0 : input.degrees[1];
  const int third_max = dimension == 2 ? 0 : input.degrees[2];
  for (int first_fixed = 0; first_fixed <= first_max; ++first_fixed) {
    for (int second_fixed = 0; second_fixed <= second_max;
         ++second_fixed) {
      for (int third_fixed = 0; third_fixed <= third_max; ++third_fixed) {
        std::vector<std::vector<Rational>> triangle(
            static_cast<std::size_t>(degree + 1));
        triangle[0].resize(static_cast<std::size_t>(degree + 1));
        for (int index = 0; index <= degree; ++index) {
          const int first = dimension == 0 ? index : first_fixed;
          const int second = dimension == 1 ? index : second_fixed;
          const int third = dimension == 2 ? index : third_fixed;
          triangle[0][static_cast<std::size_t>(index)] =
              input.values[grid_index(
                  input.degrees, first, second, third)];
        }
        for (int level = 1; level <= degree; ++level) {
          triangle[static_cast<std::size_t>(level)].resize(
              static_cast<std::size_t>(degree - level + 1));
          for (int index = 0; index <= degree - level; ++index) {
            triangle[static_cast<std::size_t>(level)]
                    [static_cast<std::size_t>(index)] =
                (
                    triangle[static_cast<std::size_t>(level - 1)]
                            [static_cast<std::size_t>(index)]
                    + triangle[static_cast<std::size_t>(level - 1)]
                              [static_cast<std::size_t>(index + 1)]
                ) / 2;
          }
        }
        for (int index = 0; index <= degree; ++index) {
          const int first = dimension == 0 ? index : first_fixed;
          const int second = dimension == 1 ? index : second_fixed;
          const int third = dimension == 2 ? index : third_fixed;
          left.values[grid_index(left.degrees, first, second, third)] =
              triangle[static_cast<std::size_t>(index)][0];
          right.values[grid_index(right.degrees, first, second, third)] =
              triangle[static_cast<std::size_t>(degree - index)]
                      [static_cast<std::size_t>(index)];
        }
      }
    }
  }
  return {std::move(left), std::move(right)};
}

struct SubdivisionResult {
  std::size_t nodes = 0U;
  std::size_t leaves = 0U;
  std::size_t unresolved = 0U;
  int maximum_depth = 0;
  bool has_first_unresolved = false;
  std::array<unsigned long long, 3> first_cell{0ULL, 0ULL, 0ULL};
  std::array<int, 3> first_splits{0, 0, 0};
};

void certify_subdivision(const BernsteinGrid& grid, const int depth,
                         const int depth_limit,
                         const std::array<unsigned long long, 3>& cell,
                         const std::array<int, 3>& splits,
                         SubdivisionResult& result) {
  ++result.nodes;
  result.maximum_depth = std::max(result.maximum_depth, depth);
  const bool nonnegative = std::all_of(
      grid.values.begin(), grid.values.end(),
      [](const Rational& value) { return value >= 0; });
  if (nonnegative) {
    ++result.leaves;
    return;
  }
  if (depth >= depth_limit) {
    ++result.unresolved;
    if (!result.has_first_unresolved) {
      result.has_first_unresolved = true;
      result.first_cell = cell;
      result.first_splits = splits;
    }
    return;
  }
  int dimension = depth % 3;
  for (int offset = 0; offset < 3; ++offset) {
    const int candidate = (dimension + offset) % 3;
    if (grid.degrees[static_cast<std::size_t>(candidate)] > 0) {
      dimension = candidate;
      break;
    }
  }
  auto [left, right] = split_grid(grid, dimension);
  std::array<unsigned long long, 3> left_cell = cell;
  std::array<unsigned long long, 3> right_cell = cell;
  std::array<int, 3> child_splits = splits;
  left_cell[static_cast<std::size_t>(dimension)] *= 2ULL;
  right_cell[static_cast<std::size_t>(dimension)] =
      right_cell[static_cast<std::size_t>(dimension)] * 2ULL + 1ULL;
  ++child_splits[static_cast<std::size_t>(dimension)];
  certify_subdivision(
      left, depth + 1, depth_limit, left_cell, child_splits, result);
  certify_subdivision(
      right, depth + 1, depth_limit, right_cell, child_splits, result);
}

int analyze_support_three_bernstein() {
  std::size_t cases = 0U;
  std::size_t coefficients = 0U;
  std::size_t negative = 0U;
  int first_antidiagonal = -1;
  int first_depth = -1;
  BernsteinResult first_result;
  std::size_t subdivision_nodes = 0U;
  std::size_t subdivision_leaves = 0U;
  std::size_t subdivision_unresolved = 0U;
  int subdivision_maximum_depth = 0;
  for (int antidiagonal = 1; antidiagonal <= 11; ++antidiagonal) {
    for (int depth = 0; 2 * depth < antidiagonal; ++depth) {
      const Polynomial target = substitute_differences(
          radial_polynomial(3, antidiagonal, depth), 3);
      const BigPolynomial parameterized =
          support_three_ratio_parameterization(target);
      const BernsteinResult result = bernstein_cube(parameterized);
      SubdivisionResult subdivision;
      certify_subdivision(
          bernstein_grid(parameterized), 0, 30,
          {0ULL, 0ULL, 0ULL}, {0, 0, 0}, subdivision);
      ++cases;
      coefficients += result.coefficients;
      negative += result.negative;
      subdivision_nodes += subdivision.nodes;
      subdivision_leaves += subdivision.leaves;
      subdivision_unresolved += subdivision.unresolved;
      subdivision_maximum_depth = std::max(
          subdivision_maximum_depth, subdivision.maximum_depth);
      if (subdivision.unresolved > 0U) {
        std::cout
            << "support_three_unresolved"
            << " antidiagonal=" << antidiagonal
            << " depth=" << depth
            << " nodes=" << subdivision.nodes
            << " leaves=" << subdivision.leaves
            << " unresolved=" << subdivision.unresolved
            << " maximum_depth=" << subdivision.maximum_depth
            << " first_cell=[" << subdivision.first_cell[0] << ','
            << subdivision.first_cell[1] << ','
            << subdivision.first_cell[2] << ']'
            << " first_splits=[" << subdivision.first_splits[0] << ','
            << subdivision.first_splits[1] << ','
            << subdivision.first_splits[2] << ']'
            << '\n';
      }
      if (result.negative > 0U && first_antidiagonal < 0) {
        first_antidiagonal = antidiagonal;
        first_depth = depth;
        first_result = result;
      }
    }
  }
  std::cout
      << "SU2_AUTOCORRELATION_SUPPORT_THREE_BERNSTEIN"
      << " cases=" << cases
      << " coefficients=" << coefficients
      << " negative=" << negative
      << " subdivision_nodes=" << subdivision_nodes
      << " subdivision_leaves=" << subdivision_leaves
      << " subdivision_unresolved=" << subdivision_unresolved
      << " subdivision_maximum_depth=" << subdivision_maximum_depth;
  if (first_antidiagonal >= 0) {
    std::cout
        << " first_antidiagonal=" << first_antidiagonal
        << " first_depth=" << first_depth
        << " first_degrees=[" << first_result.degrees[0] << ','
        << first_result.degrees[1] << ',' << first_result.degrees[2] << ']'
        << " first_minimum_index=[" << first_result.minimum_index[0] << ','
        << first_result.minimum_index[1] << ','
        << first_result.minimum_index[2] << ']'
        << " first_minimum=" << first_result.minimum.numerator()
        << '/' << first_result.minimum.denominator();
  }
  std::cout
      << " result="
      << (
          subdivision_unresolved == 0U
          ? "PASS"
          : "SUBDIVISION_INCOMPLETE"
      )
      << '\n';
  return EXIT_SUCCESS;
}

int dump_support_three_hard() {
  for (int antidiagonal = 1; antidiagonal <= 2; ++antidiagonal) {
    const Polynomial target = substitute_differences(
        radial_polynomial(3, antidiagonal, 0), 3);
    const BigPolynomial parameterized =
        support_three_parameterization(target);
    std::cout
        << "support_three_parameterized"
        << " antidiagonal=" << antidiagonal
        << " polynomial={";
    bool first = true;
    for (const auto& [exponent, coefficient] : parameterized) {
      if (coefficient == 0) {
        continue;
      }
      if (!first) {
        std::cout << ',';
      }
      first = false;
      std::cout << '[' << exponent[0] << ',' << exponent[1] << ','
                << exponent[2] << "]:" << coefficient;
    }
    std::cout << "}\n";
  }
  return EXIT_SUCCESS;
}

int analyze_ratio_cube_bernstein(const int support,
                                 const int depth_limit,
                                 const int maximum_antidiagonal) {
  std::size_t cases = 0U;
  std::size_t coefficients = 0U;
  std::size_t initial_negative = 0U;
  std::size_t nodes = 0U;
  std::size_t leaves = 0U;
  std::size_t unresolved = 0U;
  int maximum_depth = 0;
  int first_unresolved_antidiagonal = -1;
  int first_unresolved_depth = -1;
  for (int antidiagonal = 1; antidiagonal <= maximum_antidiagonal;
       ++antidiagonal) {
    for (int depth = 0; 2 * depth < antidiagonal; ++depth) {
      const Polynomial target = substitute_differences(
          radial_polynomial(support, antidiagonal, depth), support);
      const BigPolynomial parameterized =
          ratio_cube_parameterization(target, support);
      const NDBernsteinGrid grid =
          nd_bernstein_grid(parameterized, support);
      ++cases;
      coefficients += grid.values.size();
      initial_negative += static_cast<std::size_t>(std::count_if(
          grid.values.begin(), grid.values.end(),
          [](const Rational& value) { return value < 0; }));
      NDSubdivisionResult result;
      nd_certify_subdivision(grid, 0, depth_limit, result);
      nodes += result.nodes;
      leaves += result.leaves;
      unresolved += result.unresolved;
      maximum_depth = std::max(maximum_depth, result.maximum_depth);
      if (result.unresolved > 0U && first_unresolved_antidiagonal < 0) {
        first_unresolved_antidiagonal = antidiagonal;
        first_unresolved_depth = depth;
      }
      if (result.nodes > 1U) {
        std::cout
            << "ratio_cube_subdivision"
            << " support=" << support
            << " antidiagonal=" << antidiagonal
            << " depth=" << depth
            << " coefficients=" << grid.values.size()
            << " nodes=" << result.nodes
            << " leaves=" << result.leaves
            << " unresolved=" << result.unresolved
            << " maximum_depth=" << result.maximum_depth
            << '\n';
      }
    }
  }
  std::cout
      << "SU2_AUTOCORRELATION_RATIO_CUBE_BERNSTEIN"
      << " support=" << support
      << " cases=" << cases
      << " coefficients=" << coefficients
      << " initial_negative=" << initial_negative
      << " nodes=" << nodes
      << " leaves=" << leaves
      << " unresolved=" << unresolved
      << " maximum_depth=" << maximum_depth;
  if (first_unresolved_antidiagonal >= 0) {
    std::cout
        << " first_unresolved_antidiagonal="
        << first_unresolved_antidiagonal
        << " first_unresolved_depth=" << first_unresolved_depth;
  }
  std::cout
      << " result=" << (unresolved == 0U ? "PASS" : "INCOMPLETE")
      << '\n';
  return EXIT_SUCCESS;
}

struct CertificateResult {
  bool satisfiable = false;
  std::size_t generators = 0U;
  std::size_t active = 0U;
};

CertificateResult certify(const Polynomial& target, const int support) {
  const int variables = support + 1;
  const std::vector<Polynomial> curvature =
      log_concavity_generators(support);
  std::vector<Polynomial> generators;
  for (const Exponent& exponent : homogeneous_exponents(variables, 4)) {
    generators.push_back(monomial_polynomial(exponent));
  }
  const std::vector<Exponent> quadratic_monomials =
      homogeneous_exponents(variables, 2);
  for (const Polynomial& generator : curvature) {
    for (const Exponent& exponent : quadratic_monomials) {
      generators.push_back(
          multiply(generator, monomial_polynomial(exponent)));
    }
  }
  for (std::size_t first = 0; first < curvature.size(); ++first) {
    for (std::size_t second = first; second < curvature.size(); ++second) {
      generators.push_back(multiply(curvature[first], curvature[second]));
    }
  }
  z3::context context;
  z3::solver solver(context);
  std::vector<z3::expr> weights;
  weights.reserve(generators.size());
  for (std::size_t index = 0; index < generators.size(); ++index) {
    weights.push_back(
        context.real_const(("w_" + std::to_string(index)).c_str()));
    solver.add(weights.back() >= 0);
  }
  for (const Exponent& exponent : homogeneous_exponents(variables, 4)) {
    z3::expr sum = context.real_val(0);
    for (std::size_t index = 0; index < generators.size(); ++index) {
      const auto iterator = generators[index].find(exponent);
      if (iterator != generators[index].end() && iterator->second != 0) {
        sum = sum + weights[index] * context.int_val(
                                           static_cast<std::int64_t>(
                                               iterator->second));
      }
    }
    const auto target_iterator = target.find(exponent);
    const long long target_coefficient =
        target_iterator == target.end() ? 0 : target_iterator->second;
    solver.add(sum == context.int_val(
                          static_cast<std::int64_t>(target_coefficient)));
  }

  CertificateResult result;
  result.generators = generators.size();
  result.satisfiable = solver.check() == z3::sat;
  if (result.satisfiable) {
    const z3::model model = solver.get_model();
    for (const z3::expr& weight : weights) {
      const z3::expr value = model.eval(weight, true);
      if (value.to_string() != "0") {
        ++result.active;
      }
    }
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-root-outer-bernstein-obstruction"
    ) {
      return replay_root_outer_bernstein_obstruction();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-root-outer-first-variation"
    ) {
      return replay_root_outer_first_variation();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--support-two-gap-prefix-certificate"
    ) {
      return analyze_support_two_gap_prefix_certificate();
    }
    if (
        argc == 3
        && std::string{argv[1]}
               == "--support-three-gap-prefix-certificate-part"
    ) {
      return analyze_support_three_gap_prefix_certificate_part(
          parse_positive(argv[2], "first_label"));
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--support-four-exceptional-gap-prefix-certificate"
    ) {
      return analyze_support_four_exceptional_gap_prefix_certificate();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--support-four-161-saturated-blowup-certificate"
    ) {
      return analyze_support_four_161_saturated_blowup_certificate();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-gap-prefix-append-low-coefficients"
    ) {
      return replay_gap_prefix_append_low_coefficients();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-gap-prefix-wall-cubic-reserve"
    ) {
      return replay_gap_prefix_wall_cubic_reserve();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-wall-121-current-normal-form"
    ) {
      return replay_wall_121_current_normal_form();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-wall-121-saturated-recurrence-obstruction"
    ) {
      return replay_wall_121_saturated_recurrence_obstruction();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-wall-121-q2-floor"
    ) {
      return replay_wall_121_q2_floor();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-wall-121-small-ratio-payment"
    ) {
      return replay_wall_121_small_ratio_payment();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-gap-prefix-elevation-obstruction"
    ) {
      return replay_support_three_gap_prefix_elevation_obstruction();
    }
    if (
        argc == 2
        && std::string{argv[1]} == "--support-three-bernstein"
    ) {
      return analyze_support_three_bernstein();
    }
    if (
        argc == 2
        && std::string{argv[1]} == "--dump-support-three-hard"
    ) {
      return dump_support_three_hard();
    }
    if (
        argc == 2
        && std::string{argv[1]} == "--support-four-bernstein"
    ) {
      return analyze_ratio_cube_bernstein(4, 16, 15);
    }
    if (
        argc == 2
        && std::string{argv[1]} == "--support-five-bernstein"
    ) {
      return analyze_ratio_cube_bernstein(5, 20, 19);
    }
    if (
        argc == 2
        && std::string{argv[1]} == "--support-five-lower-bernstein"
    ) {
      return analyze_ratio_cube_bernstein(5, 20, 8);
    }
    const int maximum_support =
        argc >= 2 ? parse_positive(argv[1], "maximum_support") : 5;
    const std::string mode =
        argc >= 3 ? std::string{argv[2]} : std::string{};
    const bool dump_uncertified = mode == "--dump-uncertified";
    const bool dump_all = mode == "--dump-all";
    if (argc >= 3 && !dump_uncertified && !dump_all) {
      throw std::invalid_argument(
          "usage: analyze_su2_autocorrelation_lc_certificate "
          "[maximum_support] [--dump-uncertified|--dump-all]");
    }
    if (argc > 3) {
      throw std::invalid_argument(
          "usage: analyze_su2_autocorrelation_lc_certificate "
          "[maximum_support] [--dump-uncertified|--dump-all]");
    }

    std::size_t cases = 0U;
    std::size_t certified = 0U;
    std::size_t total_generators = 0U;
    std::size_t total_active = 0U;
    int first_support = -1;
    int first_antidiagonal = -1;
    int first_depth = -1;
    std::vector<std::array<int, 3>> uncertified_cases;
    for (int support = 1; support <= maximum_support; ++support) {
      for (int antidiagonal = 1; antidiagonal <= 4 * support - 1;
           ++antidiagonal) {
        for (int depth = 0; 2 * depth < antidiagonal; ++depth) {
          const Polynomial target = substitute_differences(
              radial_polynomial(support, antidiagonal, depth), support);
          const CertificateResult result = certify(target, support);
          ++cases;
          total_generators += result.generators;
          total_active += result.active;
          if (result.satisfiable) {
            ++certified;
          } else if (first_support < 0) {
            first_support = support;
            first_antidiagonal = antidiagonal;
            first_depth = depth;
          }
          if (!result.satisfiable && uncertified_cases.size() < 32U) {
            uncertified_cases.push_back(
                {support, antidiagonal, depth});
          }
          if (dump_all || (!result.satisfiable && dump_uncertified)) {
            std::cout
                << (result.satisfiable ? "certified" : "uncertified")
                << " support=" << support
                << " antidiagonal=" << antidiagonal
                << " depth=" << depth
                << " polynomial={";
            bool first = true;
            for (const auto& [exponent, coefficient] : target) {
              if (coefficient == 0) {
                continue;
              }
              if (!first) {
                std::cout << ',';
              }
              first = false;
              std::cout << '[';
              for (std::size_t index = 0; index < exponent.size(); ++index) {
                if (index != 0U) {
                  std::cout << ',';
                }
                std::cout << exponent[index];
              }
              std::cout << "]:" << coefficient;
            }
            std::cout << "}\n";
          }
        }
      }
    }
    std::cout
        << "SU2_AUTOCORRELATION_LC_CERTIFICATE"
        << " maximum_support=" << maximum_support
        << " cases=" << cases
        << " certified=" << certified
        << " total_generators=" << total_generators
        << " total_active=" << total_active;
    if (first_support >= 0) {
      std::cout
          << " first_uncertified_support=" << first_support
          << " first_uncertified_antidiagonal=" << first_antidiagonal
          << " first_uncertified_depth=" << first_depth;
      std::cout << " uncertified_cases={";
      for (std::size_t index = 0; index < uncertified_cases.size(); ++index) {
        if (index != 0U) {
          std::cout << ',';
        }
        std::cout << '[' << uncertified_cases[index][0] << ','
                  << uncertified_cases[index][1] << ','
                  << uncertified_cases[index][2] << ']';
      }
      std::cout << '}';
    }
    std::cout
        << " result="
        << (cases == certified ? "ALL_CERTIFIED_BOUNDED"
                               : "TRUNCATED_PREORDER_INSUFFICIENT")
        << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
