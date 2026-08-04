#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <set>
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

int parse_nonnegative(const char* text, const std::string& name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value, &consumed, 10);
  if (consumed != value.size() || parsed < 0) {
    throw std::invalid_argument(name + " must be a nonnegative integer");
  }
  return static_cast<int>(parsed);
}

std::vector<int> parse_csv_nonnegative(
    const char* text,
    const std::string& name
) {
  const std::string input{text};
  std::vector<int> result;
  std::size_t begin = 0U;
  while (begin < input.size()) {
    const std::size_t end = input.find(',', begin);
    const std::string part = input.substr(
        begin,
        end == std::string::npos ? std::string::npos : end - begin
    );
    if (part.empty()) {
      throw std::invalid_argument(name + " contains an empty coordinate");
    }
    std::size_t consumed = 0U;
    const long value = std::stol(part, &consumed, 10);
    if (consumed != part.size() || value < 0) {
      throw std::invalid_argument(name + " has an invalid coordinate");
    }
    result.push_back(static_cast<int>(value));
    if (end == std::string::npos) {
      break;
    }
    begin = end + 1U;
  }
  if (result.empty()) {
    throw std::invalid_argument(name + " is empty");
  }
  return result;
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

Integer integer_power(Integer base, int exponent) {
  Integer result = 1;
  while (exponent > 0) {
    if ((exponent & 1) != 0) {
      result *= base;
    }
    base *= base;
    exponent /= 2;
  }
  return result;
}

BigPolynomial restrict_to_dyadic_cell(
    const BigPolynomial& input,
    const std::vector<int>& cell,
    const std::vector<int>& splits
) {
  if (input.empty() || input.begin()->first.empty()
      || cell.size() != input.begin()->first.size()
      || splits.size() != cell.size()) {
    throw std::invalid_argument("invalid dyadic-cell restriction");
  }
  BigPolynomial result = input;
  for (std::size_t variable = 0U; variable < cell.size(); ++variable) {
    if (splits[variable] < 0) {
      throw std::invalid_argument("negative dyadic split count");
    }
    const Integer denominator = integer_power(Integer(2), splits[variable]);
    if (Integer(cell[variable]) >= denominator) {
      throw std::invalid_argument("dyadic cell lies outside the unit cube");
    }
    int degree = 0;
    for (const auto& [exponent, coefficient] : result) {
      static_cast<void>(coefficient);
      degree = std::max(degree, exponent[variable]);
    }
    BigPolynomial transformed;
    for (const auto& [exponent, coefficient] : result) {
      const int power = exponent[variable];
      const Integer scale = coefficient
          * integer_power(denominator, degree - power);
      for (int replacement = 0; replacement <= power; ++replacement) {
        Exponent output = exponent;
        output[variable] = replacement;
        transformed[std::move(output)] += scale
            * binomial(power, replacement)
            * integer_power(Integer(cell[variable]), power - replacement);
      }
    }
    result = std::move(transformed);
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

int replay_ordinary_12_complete_current_sign() {
  constexpr int minimum_support = 1;
  constexpr int maximum_support = 12;
  const auto canonical_quadratic = [](const QuadraticTerms& polynomial) {
    QuadraticTerms result;
    for (const auto& [monomial, coefficient] : polynomial) {
      if (coefficient != 0) {
        result.emplace(monomial, coefficient);
      }
    }
    return result;
  };
  const auto canonical_quartic = [](const QuarticTerms& polynomial) {
    QuarticTerms result;
    for (const auto& [monomial, coefficient] : polynomial) {
      if (coefficient != 0) {
        result.emplace(monomial, coefficient);
      }
    }
    return result;
  };
  const auto add_quadratic = [](
      QuadraticTerms& target,
      const QuadraticTerms& source,
      const long long scale) {
    for (const auto& [monomial, coefficient] : source) {
      target[monomial] += scale * coefficient;
    }
  };
  const auto ordinary_square_coefficient = [](
      const int support,
      const int label) {
    QuadraticTerms result;
    for (int first = 0; first <= support; ++first) {
      for (int second = first; second <= support; ++second) {
        if (
            std::abs(first - second) <= label
            && label <= first + second
        ) {
          result[pair_monomial(first, second)] +=
              first == second ? 1 : 2;
        }
      }
    }
    return result;
  };

  for (int support = minimum_support;
       support <= maximum_support; ++support) {
    const QuadraticTerms c_zero =
        ordinary_square_coefficient(support, 0);
    const QuadraticTerms c_one =
        ordinary_square_coefficient(support, 1);
    const QuadraticTerms c_two =
        ordinary_square_coefficient(support, 2);
    const QuadraticTerms c_three =
        ordinary_square_coefficient(support, 3);

    QuadraticTerms fusion_interval = c_one;
    add_quadratic(fusion_interval, c_two, 1);
    add_quadratic(fusion_interval, c_three, 1);
    QuarticTerms character_current;
    add_product(character_current, c_zero, fusion_interval, 1);
    add_product(character_current, c_one, c_two, -1);
    if (canonical_quartic(character_current).empty()) {
      throw std::runtime_error(
          "ordinary (1,2) complete current vanished unexpectedly at support "
          + std::to_string(support));
    }

    QuadraticTerms adjacent = c_one;
    add_quadratic(adjacent, c_zero, -1);
    QuadraticTerms adjacent_expected;
    for (int index = 0; index < support; ++index) {
      adjacent_expected[pair_monomial(index, index + 1)] += 2;
    }
    adjacent_expected[pair_monomial(0, 0)] -= 1;
    if (
        canonical_quadratic(adjacent)
        != canonical_quadratic(adjacent_expected)
    ) {
      throw std::runtime_error(
          "ordinary (1,2) alpha-threshold identity mismatch "
          "at support " + std::to_string(support));
    }

    QuadraticTerms middle_difference = c_one;
    add_quadratic(middle_difference, c_three, 1);
    add_quadratic(middle_difference, c_two, -1);
    QuadraticTerms middle_expected;
    for (int index = 2; index <= support; ++index) {
      middle_expected[pair_monomial(index, index)] += 1;
    }
    for (int index = 0; index < support; ++index) {
      middle_expected[pair_monomial(index, index + 1)] += 2;
    }
    if (support >= 2) {
      middle_expected[pair_monomial(0, 2)] -= 2;
    }
    for (int index = 0; index + 3 <= support; ++index) {
      middle_expected[pair_monomial(index, index + 3)] += 2;
    }
    if (
        canonical_quadratic(middle_difference)
        != canonical_quadratic(middle_expected)
    ) {
      throw std::runtime_error(
          "ordinary (1,2) middle-branch identity mismatch "
          "at support " + std::to_string(support));
    }

    QuadraticTerms spectral_left = c_zero;
    add_quadratic(spectral_left, c_one, 1);
    QuadraticTerms spectral_squares;
    for (int index = 0; index <= support; ++index) {
      spectral_squares[pair_monomial(index, index)] += 1;
      if (index + 1 <= support) {
        spectral_squares[pair_monomial(index + 1, index + 1)] += 1;
        spectral_squares[pair_monomial(index, index + 1)] += 2;
      }
    }
    if (
        canonical_quadratic(spectral_left)
        != canonical_quadratic(spectral_squares)
    ) {
      throw std::runtime_error(
          "ordinary (1,2) spectral-square identity mismatch "
          "at support " + std::to_string(support));
    }
  }

  std::cout
      << "SU2_ORDINARY_12_COMPLETE_CURRENT_SIGN"
      << " supports=1..12"
      << " current=c0*(c1+c2+c3)-c1*c2"
      << " alpha_cases=(alpha<=1,1<=alpha<=2,alpha>=2)"
      << " middle_identity=c1+c3-c2"
      << " threshold_identity=c1-c0=2*B1-p0^2"
      << " spectral_identity=<z,(N1+I)z>=sum_i(z_i+z_(i+1))^2"
      << " conclusion=J(1,2)>=0"
      << " result=PASS_EXACT"
      << '\n';
  return EXIT_SUCCESS;
}

int replay_ordinary_one_row_universal_obstruction() {
  const std::array<Integer, 7> root{
      Integer(1), Integer(2), Integer(2), Integer(0),
      Integer(1), Integer(0), Integer(0)};
  std::array<Integer, 13> square{};
  for (std::size_t left = 0U; left < root.size(); ++left) {
    for (std::size_t right = 0U; right < root.size(); ++right) {
      const std::size_t lower = left > right ? left - right : right - left;
      for (std::size_t target = lower; target <= left + right; ++target) {
        square[target] += root[left] * root[right];
      }
    }
  }
  const std::array<Integer, 13> expected_square{
      Integer(10), Integer(21), Integer(25), Integer(21),
      Integer(15), Integer(9), Integer(5), Integer(1), Integer(1),
      Integer(0), Integer(0), Integer(0), Integer(0)};
  const Integer current = square[0] * (square[7] + square[8] + square[9])
      - square[1] * square[8];
  if (square != expected_square || current != Integer(-1)) {
    throw std::runtime_error(
        "ordinary one-row universal-current obstruction mismatch");
  }
  std::cout
      << "SU2_ORDINARY_ONE_ROW_UNIVERSAL_OBSTRUCTION"
      << " root=(1,2,2,0,1)"
      << " target=(1,8)"
      << " square=(10,21,25,21,15,9,5,1,1)"
      << " current=-1"
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

int replay_wall_121_geometric_demand_ceiling() {
  const auto univariate = [](const std::vector<Rational>& coefficients) {
    RationalPolynomial result;
    for (std::size_t degree = 0; degree < coefficients.size(); ++degree) {
      if (coefficients[degree] != 0) {
        result[Exponent{static_cast<int>(degree)}] = coefficients[degree];
      }
    }
    return result;
  };
  const auto shifted_coefficients = [](const std::vector<Integer>& powers) {
    std::vector<Integer> result(powers.size(), 0);
    for (std::size_t power = 0; power < powers.size(); ++power) {
      for (std::size_t degree = 0; degree <= power; ++degree) {
        result[degree] +=
            powers[power]
            * binomial(static_cast<int>(power), static_cast<int>(degree));
      }
    }
    return result;
  };

  constexpr int variables = 2;
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial r = nd_variable(variables, 0);
  const BigPolynomial x = nd_variable(variables, 1);
  const BigPolynomial r3 = nd_power(r, 3);
  const BigPolynomial r4 = nd_power(r, 4);
  const BigPolynomial r5 = nd_power(r, 5);
  BigPolynomial g_zero = r3;
  big_add_scaled(g_zero, r4, 1);
  BigPolynomial g_one;
  big_add_scaled(g_one, r3, -3);
  big_add_scaled(g_one, r4, -2);
  big_add_scaled(g_one, r5, 1);
  BigPolynomial g_two = r;
  big_add_scaled(g_two, r3, 1);
  BigPolynomial g_three = one;
  big_add_scaled(g_three, r, 1);
  BigPolynomial g = g_zero;
  big_add_scaled(g, big_multiply(g_one, x), 1);
  big_add_scaled(g, big_multiply(g_two, nd_power(x, 2)), 1);
  big_add_scaled(g, big_multiply(g_three, nd_power(x, 3)), 1);

  const std::vector<RationalPolynomial> cubic_bernstein =
      outer_bernstein_coefficients(g, variables);
  const std::vector<RationalPolynomial> elevated =
      elevate_bernstein_coefficients(cubic_bernstein, 4);
  const std::vector<RationalPolynomial> expected{
      univariate({0, 0, 0, 1, 1}),
      univariate({0, 0, 0, Rational(1, 4), Rational(1, 2),
                  Rational(1, 4)}),
      univariate({0, Rational(1, 6), 0, Rational(-1, 3), 0,
                  Rational(1, 2)}),
      univariate({Rational(1, 4), Rational(3, 4), 0,
                  Rational(-3, 4), Rational(-1, 2),
                  Rational(3, 4)}),
      univariate({1, 2, 0, -1, -1, 1})};
  if (elevated.size() != expected.size()) {
    throw std::runtime_error(
        "wall (1,2,1) demand-ceiling Bernstein degree mismatch");
  }
  for (std::size_t index = 0; index < elevated.size(); ++index) {
    if (!equal_polynomials(elevated[index], expected[index])) {
      throw std::runtime_error(
          "wall (1,2,1) demand-ceiling Bernstein coefficient mismatch");
    }
  }

  const BigPolynomial one_variable_one = nd_constant(1, 1);
  const BigPolynomial one_variable_r = nd_variable(1, 0);
  BigPolynomial h = one_variable_one;
  big_add_scaled(h, one_variable_r, 3);
  big_add_scaled(h, nd_power(one_variable_r, 3), -3);
  big_add_scaled(h, nd_power(one_variable_r, 4), -2);
  big_add_scaled(h, nd_power(one_variable_r, 5), 3);
  BigPolynomial q = one_variable_one;
  big_add_scaled(q, one_variable_r, 2);
  big_add_scaled(q, nd_power(one_variable_r, 3), -1);
  big_add_scaled(q, nd_power(one_variable_r, 4), -1);
  big_add_scaled(q, nd_power(one_variable_r, 5), 1);
  const auto scalar_bernstein = [](const BigPolynomial& polynomial) {
    const std::vector<RationalPolynomial> converted =
        outer_bernstein_coefficients(polynomial, 1);
    std::vector<Rational> result;
    result.reserve(converted.size());
    for (const RationalPolynomial& coefficient : converted) {
      const auto iterator = coefficient.find(Exponent{});
      result.push_back(
          iterator == coefficient.end() ? Rational(0) : iterator->second);
      if (coefficient.size() > (iterator == coefficient.end() ? 0U : 1U)) {
        throw std::runtime_error(
            "wall (1,2,1) demand-ceiling scalar Bernstein mismatch");
      }
    }
    return result;
  };
  const std::vector<Rational> h_unit_expected{
      1, Rational(8, 5), Rational(11, 5), Rational(5, 2),
      Rational(9, 5), 2};
  const std::vector<Rational> q_unit_expected{
      1, Rational(7, 5), Rational(9, 5), Rational(21, 10), 2, 2};
  const std::vector<Integer> h_powers{1, 3, 0, -3, -2, 3};
  const std::vector<Integer> q_powers{1, 2, 0, -1, -1, 1};
  const std::vector<Integer> h_shift_expected{2, 1, 9, 19, 13, 3};
  const std::vector<Integer> q_shift_expected{2, 0, 1, 5, 4, 1};
  if (
      scalar_bernstein(h) != h_unit_expected
      || scalar_bernstein(q) != q_unit_expected
      || shifted_coefficients(h_powers) != h_shift_expected
      || shifted_coefficients(q_powers) != q_shift_expected
      || Integer(-2) * Integer(-2) - 4 * Integer(3) * Integer(1)
             != Integer(-8)
  ) {
    throw std::runtime_error(
        "wall (1,2,1) demand-ceiling positivity certificate mismatch");
  }
  constexpr int obstruction_support = 8;
  const std::array<Integer, 10> obstruction_profile{
      0, 3, 4, 5, 5, 4, 3, 2, 1, 0};
  for (int index = 1; index < obstruction_support; ++index) {
    if (
        obstruction_profile[static_cast<std::size_t>(index)]
              * obstruction_profile[static_cast<std::size_t>(index)]
        < obstruction_profile[static_cast<std::size_t>(index - 1)]
              * obstruction_profile[static_cast<std::size_t>(index + 1)]
    ) {
      throw std::runtime_error(
          "wall (1,2,1) terminal-floor profile is not log concave");
    }
  }
  std::array<Integer, 10> obstruction_g{};
  std::array<Integer, 11> obstruction_h{};
  std::array<Integer, 10> obstruction_k{};
  for (int index = 1; index <= obstruction_support + 1; ++index) {
    obstruction_g[static_cast<std::size_t>(index)] =
        obstruction_profile[static_cast<std::size_t>(index)]
              * obstruction_profile[static_cast<std::size_t>(index)]
        - obstruction_profile[static_cast<std::size_t>(index - 1)]
              * (index + 1 < static_cast<int>(obstruction_profile.size())
                     ? obstruction_profile[static_cast<std::size_t>(index + 1)]
                     : Integer(0));
  }
  for (int index = 1; index <= obstruction_support + 2; ++index) {
    const Integer current =
        index < static_cast<int>(obstruction_profile.size())
            ? obstruction_profile[static_cast<std::size_t>(index)]
            : Integer(0);
    const Integer previous =
        obstruction_profile[static_cast<std::size_t>(index - 1)];
    const Integer left_two =
        index >= 2
            ? obstruction_profile[static_cast<std::size_t>(index - 2)]
            : Integer(0);
    const Integer next =
        index + 1 < static_cast<int>(obstruction_profile.size())
            ? obstruction_profile[static_cast<std::size_t>(index + 1)]
            : Integer(0);
    obstruction_h[static_cast<std::size_t>(index)] =
        current * previous - left_two * next;
  }
  for (int index = 1; index <= obstruction_support + 1; ++index) {
    obstruction_k[static_cast<std::size_t>(index)] =
        obstruction_g[static_cast<std::size_t>(index)]
        + obstruction_h[static_cast<std::size_t>(index)]
        + obstruction_h[static_cast<std::size_t>(index + 1)];
  }
  const std::array<Integer, 10> obstruction_g_expected{
      0, 9, 1, 5, 5, 1, 1, 1, 1, 0};
  const std::array<Integer, 10> obstruction_k_expected{
      0, 21, 18, 19, 19, 8, 5, 5, 3, 0};
  if (
      obstruction_g != obstruction_g_expected
      || obstruction_k != obstruction_k_expected
  ) {
    throw std::runtime_error(
        "wall (1,2,1) terminal-floor current reconstruction mismatch");
  }
  Integer obstruction_energy = 0;
  for (int index = 1; index <= obstruction_support; ++index) {
    obstruction_energy +=
        (obstruction_g[static_cast<std::size_t>(index)]
         - obstruction_g[static_cast<std::size_t>(index + 1)])
        * (obstruction_k[static_cast<std::size_t>(index)]
           - obstruction_k[static_cast<std::size_t>(index + 1)]);
  }
  const Integer obstruction_wall =
      obstruction_profile[1] * obstruction_profile[1]
      * obstruction_profile[1] * obstruction_profile[2];
  const Integer obstruction_floor =
      obstruction_profile[1] * obstruction_profile[1]
      * obstruction_profile[1] * obstruction_profile[1]
      + obstruction_wall;
  const Integer obstruction_a = obstruction_profile[1];
  const Integer obstruction_b = obstruction_profile[2];
  const Integer obstruction_c = obstruction_profile[3];
  const Integer obstruction_d = obstruction_profile[4];
  const Integer obstruction_c1 =
      obstruction_b * obstruction_b * obstruction_b
      + obstruction_a * obstruction_c * obstruction_c
      + obstruction_c * obstruction_c * obstruction_c
      - 3 * obstruction_a * obstruction_a * obstruction_b
      - 2 * obstruction_a * obstruction_b * obstruction_b
      - obstruction_a * obstruction_b * obstruction_c
      - obstruction_a * obstruction_b * obstruction_d
      - obstruction_b * obstruction_c * obstruction_d;
  const Integer obstruction_c2 =
      obstruction_a * obstruction_a
      + obstruction_b * obstruction_b
      + (obstruction_a + obstruction_b) * obstruction_c;
  const Integer obstruction_c3 =
      obstruction_a + obstruction_b + obstruction_c;
  const Rational obstruction_endpoint(9, 4);
  const Rational obstruction_derivative =
      Rational(obstruction_c1)
      + 2 * Rational(obstruction_c2) * obstruction_endpoint
      + 3 * Rational(obstruction_c3) * obstruction_endpoint
            * obstruction_endpoint;
  if (
      obstruction_energy != 75
      || obstruction_wall + obstruction_energy != 183
      || obstruction_floor != 189
      || obstruction_c1 != -160
      || obstruction_c2 != 60
      || obstruction_c3 != 12
      || obstruction_derivative != Rational(1169, 4)
  ) {
    throw std::runtime_error(
        "wall (1,2,1) terminal geometric-floor obstruction mismatch");
  }
  std::cout
      << "SU2_WALL_121_GEOMETRIC_DEMAND_CEILING"
      << " normalized_demand_bound=1+R"
      << " legendre_endpoints=(u=0,u=1)"
      << " u0_x_degree=3"
      << " elevated_degree=4"
      << " elevated_coefficients="
      << "(R^3*(1+R),R^3*(1+R)^2/4,"
      << "R*(1-2R^2+3R^4)/6,H(R)/4,Q(R))"
      << " H_unit_bernstein=(1,8/5,11/5,5/2,9/5,2)"
      << " Q_unit_bernstein=(1,7/5,9/5,21/10,2,2)"
      << " H_shift=(2,1,9,19,13,3)"
      << " Q_shift=(2,0,1,5,4,1)"
      << " quadratic_discriminant=-8"
      << " terminal_floor_profile=(3,4,5,5,4,3,2,1)"
      << " terminal_paired_energy=75"
      << " terminal_current=183"
      << " terminal_geometric_floor=189"
      << " terminal_floor_gap=-6"
      << " terminal_critical_derivative=1169/4"
      << " result=PASS_EXACT"
      << '\n';
  return EXIT_SUCCESS;
}

int replay_wall_121_cubic_corrected_demand_ceiling() {
  const BigPolynomial one = nd_constant(1, 1);
  const BigPolynomial r = nd_variable(1, 0);
  const BigPolynomial r_two = nd_power(r, 2);
  BigPolynomial geometric_a = one;
  big_add_scaled(geometric_a, r_two, 2);
  big_add_scaled(geometric_a, nd_power(r, 3), 1);
  BigPolynomial geometric_b;
  big_add_scaled(geometric_b, r, 3);
  big_add_scaled(geometric_b, r_two, 2);
  BigPolynomial geometric_c = one;
  big_add_scaled(geometric_c, r, 1);
  big_add_scaled(geometric_c, r_two, 1);
  BigPolynomial geometric_current = one;
  big_add_scaled(geometric_current, r, 1);
  BigPolynomial geometric_h = nd_power(geometric_a, 2);
  for (auto& [exponent, coefficient] : geometric_h) {
    static_cast<void>(exponent);
    coefficient *= 4;
  }
  big_add_scaled(
      geometric_h, big_multiply(geometric_b, geometric_c), 3);
  BigPolynomial geometric_quadratic =
      big_multiply(geometric_a, geometric_current);
  for (auto& [exponent, coefficient] : geometric_quadratic) {
    static_cast<void>(exponent);
    coefficient *= 4;
  }
  big_add_scaled(geometric_quadratic, nd_power(geometric_b, 2), -1);
  BigPolynomial geometric_margin =
      big_multiply(geometric_quadratic, nd_power(geometric_h, 4));
  BigPolynomial geometric_first_correction = big_multiply(
      big_multiply(
          big_multiply(nd_power(geometric_a, 4), geometric_c),
          nd_power(geometric_b, 3)),
      geometric_h);
  big_add_scaled(geometric_margin, geometric_first_correction, 32);
  const BigPolynomial geometric_second_correction = big_multiply(
      big_multiply(
          nd_power(geometric_a, 4), nd_power(geometric_c, 2)),
      nd_power(geometric_b, 4));
  big_add_scaled(geometric_margin, geometric_second_correction, 144);
  const std::vector<Integer> geometric_expected{
      1024, 10240, 71808, 348480, 1357460, 4230852, 11049207,
      24376076, 46302282, 76760012, 112488037, 147955104,
      175724946, 189838320, 187530543, 169190696, 139027496,
      103376016, 68592960, 39925984, 20013840, 8414336, 2848128,
      733056, 132608, 14848, 768};
  std::vector<Integer> geometric_actual(geometric_expected.size(), 0);
  for (const auto& [exponent, coefficient] : geometric_margin) {
    if (coefficient == 0) {
      continue;
    }
    if (
        exponent.size() != 1U
        || exponent[0] < 0
        || exponent[0] >= static_cast<int>(geometric_actual.size())
    ) {
      throw std::runtime_error(
          "wall (1,2,1) geometric rational-payment degree mismatch");
    }
    geometric_actual[static_cast<std::size_t>(exponent[0])] = coefficient;
  }
  if (geometric_actual != geometric_expected) {
    throw std::runtime_error(
        "wall (1,2,1) geometric rational-payment expansion mismatch");
  }

  constexpr int support = 11;
  const std::array<Integer, 13> profile{
      0, 10, 14, 16, 17, 16, 14, 11, 8, 5, 3, 1, 0};
  for (int index = 1; index <= support; ++index) {
    if (
        profile[static_cast<std::size_t>(index)]
              * profile[static_cast<std::size_t>(index)]
        < profile[static_cast<std::size_t>(index - 1)]
              * profile[static_cast<std::size_t>(index + 1)]
    ) {
      throw std::runtime_error(
          "wall (1,2,1) quadratic-envelope profile is not log concave");
    }
  }

  std::array<Integer, 14> g{};
  std::array<Integer, 15> h{};
  std::array<Integer, 14> k{};
  for (int index = 1; index <= support + 1; ++index) {
    g[static_cast<std::size_t>(index)] =
        profile[static_cast<std::size_t>(index)]
              * profile[static_cast<std::size_t>(index)]
        - profile[static_cast<std::size_t>(index - 1)]
              * (index + 1 < static_cast<int>(profile.size())
                     ? profile[static_cast<std::size_t>(index + 1)]
                     : Integer(0));
  }
  for (int index = 1; index <= support + 2; ++index) {
    const Integer current =
        index < static_cast<int>(profile.size())
            ? profile[static_cast<std::size_t>(index)]
            : Integer(0);
    const Integer previous =
        index - 1 < static_cast<int>(profile.size())
            ? profile[static_cast<std::size_t>(index - 1)]
            : Integer(0);
    const Integer left_two =
        index >= 2
            ? profile[static_cast<std::size_t>(index - 2)]
            : Integer(0);
    const Integer next =
        index + 1 < static_cast<int>(profile.size())
            ? profile[static_cast<std::size_t>(index + 1)]
            : Integer(0);
    h[static_cast<std::size_t>(index)] =
        current * previous - left_two * next;
  }
  for (int index = 1; index <= support + 1; ++index) {
    k[static_cast<std::size_t>(index)] =
        g[static_cast<std::size_t>(index)]
        + h[static_cast<std::size_t>(index)]
        + h[static_cast<std::size_t>(index + 1)];
  }

  Integer energy = 0;
  for (int index = 1; index <= support; ++index) {
    energy +=
        (g[static_cast<std::size_t>(index)]
         - g[static_cast<std::size_t>(index + 1)])
        * (k[static_cast<std::size_t>(index)]
           - k[static_cast<std::size_t>(index + 1)]);
  }
  const Integer a = profile[1];
  const Integer b = profile[2];
  const Integer c = profile[3];
  const Integer d = profile[4];
  const Integer c_zero = a * a * a * b + energy;
  const Integer c_one =
      b * b * b + a * c * c + c * c * c
      - 3 * a * a * b - 2 * a * b * b
      - a * b * c - a * b * d - b * c * d;
  const Integer a_coefficient = a * a + b * b + (a + b) * c;
  const Integer b_coefficient = -c_one;
  const Integer c_coefficient = a + b + c;
  const Integer quadratic_margin =
      4 * a_coefficient * c_zero
      - b_coefficient * b_coefficient;
  const Integer endpoint_derivative_numerator =
      c_one * b * b
      + 2 * a_coefficient * a * a * b
      + 3 * c_coefficient * a * a * a * a;
  const Integer discriminant =
      18 * c_coefficient * a_coefficient * c_one * c_zero
      - 4 * a_coefficient * a_coefficient * a_coefficient * c_zero
      + a_coefficient * a_coefficient * c_one * c_one
      - 4 * c_coefficient * c_one * c_one * c_one
      - 27 * c_coefficient * c_coefficient * c_zero * c_zero;
  const Integer denominator =
      4 * a_coefficient * a_coefficient
      + 3 * b_coefficient * c_coefficient;
  const Integer corrected_margin =
      quadratic_margin
          * denominator * denominator * denominator * denominator
      + 32 * a_coefficient * a_coefficient
          * a_coefficient * a_coefficient
          * c_coefficient * b_coefficient * b_coefficient
          * b_coefficient * denominator
      + 144 * a_coefficient * a_coefficient
          * a_coefficient * a_coefficient
          * c_coefficient * c_coefficient
          * b_coefficient * b_coefficient
          * b_coefficient * b_coefficient;
  if (
      c_zero != 17618
      || a_coefficient != 680
      || b_coefficient != 7148
      || c_coefficient != 40
      || quadratic_margin != -3172944
      || endpoint_derivative_numerator != 1702992
      || discriminant != Integer("-15163796058880")
      || corrected_margin
             != Integer("228746021443515431465475112960000")
  ) {
    throw std::runtime_error(
        "wall (1,2,1) cubic-corrected demand replay mismatch");
  }

  const Integer long_ratio_denominator = 1000;
  const std::vector<Integer> long_ratio_numerators{
      1584, 1323, 1208, 1139, 1092, 1057, 1028, 1003, 981, 960,
      939, 919, 897, 874, 848, 818, 782, 734, 670, 572, 410};
  std::vector<Rational> long_ratios;
  long_ratios.reserve(long_ratio_numerators.size());
  for (const Integer& numerator : long_ratio_numerators) {
    long_ratios.emplace_back(numerator, long_ratio_denominator);
  }
  for (std::size_t index = 1; index < long_ratios.size(); ++index) {
    if (long_ratios[index - 1U] < long_ratios[index]) {
      throw std::runtime_error(
          "wall (1,2,1) rational-ceiling obstruction is not log concave");
    }
  }
  const std::size_t long_support = long_ratios.size() + 1U;
  std::vector<Rational> long_profile(long_support + 4U, Rational(0));
  long_profile[1] = 1;
  for (std::size_t index = 0; index < long_ratios.size(); ++index) {
    long_profile[index + 2U] =
        long_profile[index + 1U] * long_ratios[index];
  }
  std::vector<Rational> long_g(long_profile.size(), Rational(0));
  std::vector<Rational> long_h(long_profile.size(), Rational(0));
  std::vector<Rational> long_k(long_profile.size(), Rational(0));
  for (std::size_t index = 1; index <= long_support + 1U; ++index) {
    long_g[index] =
        long_profile[index] * long_profile[index]
        - long_profile[index - 1U] * long_profile[index + 1U];
  }
  for (std::size_t index = 1; index <= long_support + 2U; ++index) {
    long_h[index] = long_profile[index] * long_profile[index - 1U];
    if (index >= 2U) {
      long_h[index] -=
          long_profile[index - 2U] * long_profile[index + 1U];
    }
  }
  for (std::size_t index = 1; index <= long_support + 1U; ++index) {
    long_k[index] = long_g[index] + long_h[index] + long_h[index + 1U];
  }
  Rational long_c_zero =
      long_profile[1] * long_profile[1] * long_profile[1]
      * long_profile[2];
  for (std::size_t index = 1; index <= long_support; ++index) {
    long_c_zero +=
        (long_g[index] - long_g[index + 1U])
        * (long_k[index] - long_k[index + 1U]);
  }
  const Rational& long_a = long_profile[1];
  const Rational& long_b = long_profile[2];
  const Rational& long_c = long_profile[3];
  const Rational& long_d = long_profile[4];
  const Rational long_a_coefficient =
      long_a * long_a + long_b * long_b
      + (long_a + long_b) * long_c;
  const Rational long_c_coefficient = long_a + long_b + long_c;
  const Rational long_b_coefficient =
      3 * long_a * long_a * long_b
      + 2 * long_a * long_b * long_b
      + long_a * long_b * long_c + long_a * long_b * long_d
      + long_b * long_c * long_d - long_b * long_b * long_b
      - long_a * long_c * long_c - long_c * long_c * long_c;
  const Rational long_endpoint_derivative =
      -long_b_coefficient
      + 2 * long_a_coefficient * long_a * long_a / long_b
      + 3 * long_c_coefficient * long_a * long_a * long_a * long_a
            / (long_b * long_b);
  const Rational long_h_coefficient =
      4 * long_a_coefficient * long_a_coefficient
      + 3 * long_b_coefficient * long_c_coefficient;
  const Rational long_corrected_margin =
      (4 * long_a_coefficient * long_c_zero
       - long_b_coefficient * long_b_coefficient)
          * long_h_coefficient * long_h_coefficient
          * long_h_coefficient * long_h_coefficient
      + 32 * long_a_coefficient * long_a_coefficient
          * long_a_coefficient * long_a_coefficient
          * long_c_coefficient * long_b_coefficient
          * long_b_coefficient * long_b_coefficient
          * long_h_coefficient
      + 144 * long_a_coefficient * long_a_coefficient
          * long_a_coefficient * long_a_coefficient
          * long_c_coefficient * long_c_coefficient
          * long_b_coefficient * long_b_coefficient
          * long_b_coefficient * long_b_coefficient;
  const Rational long_j_coefficient =
      2 * long_a_coefficient * long_a_coefficient
      + 3 * long_b_coefficient * long_c_coefficient;
  const Rational long_k_coefficient =
      8 * long_a_coefficient * long_a_coefficient
          * long_j_coefficient
      + 3 * long_b_coefficient * long_c_coefficient
          * long_h_coefficient;
  const Rational long_n_coefficient =
      4 * long_a_coefficient * long_b_coefficient
          * long_j_coefficient;
  const Rational long_second_corrected_margin =
      (4 * long_a_coefficient * long_c_zero
       - long_b_coefficient * long_b_coefficient)
          * long_k_coefficient * long_k_coefficient
          * long_k_coefficient * long_k_coefficient
      + 4 * long_a_coefficient * long_c_coefficient
          * long_n_coefficient * long_n_coefficient
          * long_n_coefficient * long_k_coefficient
      + 9 * long_c_coefficient * long_c_coefficient
          * long_n_coefficient * long_n_coefficient
          * long_n_coefficient * long_n_coefficient;
  const Rational long_discriminant =
      18 * long_c_coefficient * long_a_coefficient
          * (-long_b_coefficient) * long_c_zero
      - 4 * long_a_coefficient * long_a_coefficient
          * long_a_coefficient * long_c_zero
      + long_a_coefficient * long_a_coefficient
          * long_b_coefficient * long_b_coefficient
      + 4 * long_c_coefficient * long_b_coefficient
          * long_b_coefficient * long_b_coefficient
      - 27 * long_c_coefficient * long_c_coefficient
          * long_c_zero * long_c_zero;
  if (
      long_support != 22U
      || long_c_zero <= 0
      || long_b_coefficient <= 0
      || long_endpoint_derivative <= 0
      || long_corrected_margin >= 0
      || long_second_corrected_margin <= 0
      || long_discriminant >= 0
  ) {
    throw std::runtime_error(
        "wall (1,2,1) rational-ceiling obstruction mismatch");
  }

  std::cout
      << "SU2_WALL_121_CUBIC_CORRECTED_DEMAND_CEILING"
      << " lower_tau=2*A*B/(4*A^2+3*B*C)"
      << " exact_gap=B^2/(4*A)-D=C*tau^3*(1+9*C*tau/(4*A))"
      << " sufficient_margin=(4*A*C0-B^2)*H^4"
      << "+32*A^4*C*B^3*H+144*A^4*C^2*B^4"
      << " H=4*A^2+3*B*C"
      << " geometric_margin_degree=26"
      << " geometric_margin_coefficients="
      << "(1024,10240,71808,348480,1357460,4230852,11049207,"
      << "24376076,46302282,76760012,112488037,147955104,"
      << "175724946,189838320,187530543,169190696,139027496,"
      << "103376016,68592960,39925984,20013840,8414336,2848128,"
      << "733056,132608,14848,768)"
      << " obstruction_profile=(0,10,14,16,17,16,14,11,8,5,3,1)"
      << " coefficients=(C0=17618,A=680,B=7148,C=40)"
      << " quadratic_margin=-3172944"
      << " critical_derivative_numerator=1702992"
      << " discriminant=-15163796058880"
      << " corrected_margin=228746021443515431465475112960000"
      << " rational_ceiling_obstruction_support=22"
      << " rational_ceiling_obstruction_ratios_per_1000="
      << "(1584,1323,1208,1139,1092,1057,1028,1003,981,960,"
      << "939,919,897,874,848,818,782,734,670,572,410)"
      << " rational_ceiling_obstruction_margin_negative=1"
      << " second_ceiling_margin_positive=1"
      << " rational_ceiling_obstruction_discriminant_negative=1"
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
  bool has_first_unresolved = false;
  std::vector<std::uint64_t> first_unresolved_cell;
  std::vector<int> first_unresolved_splits;
};

void nd_certify_subdivision_impl(
    const NDBernsteinGrid& grid,
    const int depth,
    const int depth_limit,
    const std::vector<std::uint64_t>& cell,
    const std::vector<int>& splits,
    NDSubdivisionResult& result
) {
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
    if (!result.has_first_unresolved) {
      result.has_first_unresolved = true;
      result.first_unresolved_cell = cell;
      result.first_unresolved_splits = splits;
    }
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
  std::vector<std::uint64_t> left_cell = cell;
  std::vector<std::uint64_t> right_cell = cell;
  std::vector<int> child_splits = splits;
  const std::size_t selected = static_cast<std::size_t>(dimension);
  left_cell[selected] *= UINT64_C(2);
  right_cell[selected] = right_cell[selected] * UINT64_C(2) + UINT64_C(1);
  ++child_splits[selected];
  nd_certify_subdivision_impl(
      left, depth + 1, depth_limit, left_cell, child_splits, result);
  nd_certify_subdivision_impl(
      right, depth + 1, depth_limit, right_cell, child_splits, result);
}

void nd_certify_subdivision(const NDBernsteinGrid& grid, const int depth,
                            const int depth_limit,
                            NDSubdivisionResult& result) {
  const std::size_t dimensions = grid.degrees.size();
  nd_certify_subdivision_impl(
      grid, depth, depth_limit,
      std::vector<std::uint64_t>(dimensions, UINT64_C(0)),
      std::vector<int>(dimensions, 0), result);
}

std::size_t nd_negative_count(const NDBernsteinGrid& grid) {
  return static_cast<std::size_t>(std::count_if(
      grid.values.begin(), grid.values.end(),
      [](const Rational& value) { return value < 0; }));
}

void nd_certify_subdivision_best_impl(
    const NDBernsteinGrid& grid,
    const int depth,
    const int depth_limit,
    const std::vector<std::uint64_t>& cell,
    const std::vector<int>& splits,
    NDSubdivisionResult& result
) {
  ++result.nodes;
  result.maximum_depth = std::max(result.maximum_depth, depth);
  if (nd_negative_count(grid) == 0U) {
    ++result.leaves;
    return;
  }
  if (depth >= depth_limit) {
    ++result.unresolved;
    if (!result.has_first_unresolved) {
      result.has_first_unresolved = true;
      result.first_unresolved_cell = cell;
      result.first_unresolved_splits = splits;
    }
    return;
  }

  int dimension = -1;
  std::size_t best_score = std::numeric_limits<std::size_t>::max();
  for (std::size_t candidate = 0U;
       candidate < grid.degrees.size();
       ++candidate) {
    if (grid.degrees[candidate] == 0) {
      continue;
    }
    const auto [left, right] = nd_split_grid(
        grid, static_cast<int>(candidate)
    );
    const std::size_t score = nd_negative_count(left)
        + nd_negative_count(right);
    if (score < best_score) {
      best_score = score;
      dimension = static_cast<int>(candidate);
    }
  }
  if (dimension < 0) {
    throw std::runtime_error("Bernstein grid has no splittable dimension");
  }

  auto [left, right] = nd_split_grid(grid, dimension);
  std::vector<std::uint64_t> left_cell = cell;
  std::vector<std::uint64_t> right_cell = cell;
  std::vector<int> child_splits = splits;
  const std::size_t selected = static_cast<std::size_t>(dimension);
  left_cell[selected] *= UINT64_C(2);
  right_cell[selected] = right_cell[selected] * UINT64_C(2) + UINT64_C(1);
  ++child_splits[selected];
  nd_certify_subdivision_best_impl(
      left, depth + 1, depth_limit, left_cell, child_splits, result);
  nd_certify_subdivision_best_impl(
      right, depth + 1, depth_limit, right_cell, child_splits, result);
}

void nd_certify_subdivision_best(const NDBernsteinGrid& grid,
                                 const int depth_limit,
                                 NDSubdivisionResult& result) {
  const std::size_t dimensions = grid.degrees.size();
  nd_certify_subdivision_best_impl(
      grid, 0, depth_limit,
      std::vector<std::uint64_t>(dimensions, UINT64_C(0)),
      std::vector<int>(dimensions, 0), result);
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

struct NDIntegerBernsteinGrid {
  std::vector<int> degrees;
  std::vector<Integer> values;
};

struct NDIntegerGridConversion {
  NDIntegerBernsteinGrid grid;
  std::size_t distinct_denominators = 0U;
  std::size_t common_denominator_bits = 0U;
};

Integer integer_gcd(Integer left, Integer right) {
  if (left < 0) {
    left = -left;
  }
  if (right < 0) {
    right = -right;
  }
  while (right != 0) {
    const Integer remainder = left % right;
    left = right;
    right = remainder;
  }
  return left;
}

NDIntegerGridConversion integer_bernstein_grid(
    const NDBernsteinGrid& input) {
  std::set<Integer> denominators;
  for (const Rational& value : input.values) {
    denominators.insert(value.denominator());
  }
  Integer common_denominator = 1;
  for (const Integer& denominator : denominators) {
    common_denominator =
        (common_denominator / integer_gcd(
            common_denominator, denominator))
        * denominator;
  }
  NDIntegerBernsteinGrid grid{input.degrees, {}};
  grid.values.reserve(input.values.size());
  for (const Rational& value : input.values) {
    grid.values.push_back(
        value.numerator()
        * (common_denominator / value.denominator()));
  }
  return {
      std::move(grid), denominators.size(),
      boost::multiprecision::msb(common_denominator) + 1U};
}

std::pair<NDIntegerBernsteinGrid, NDIntegerBernsteinGrid>
nd_split_integer_grid(
    const NDIntegerBernsteinGrid& input, const int dimension) {
  NDIntegerBernsteinGrid left{input.degrees, input.values};
  NDIntegerBernsteinGrid right{input.degrees, input.values};
  const std::vector<Exponent> indices = grid_indices(input.degrees);
  const int degree =
      input.degrees[static_cast<std::size_t>(dimension)];
  for (const Exponent& base : indices) {
    if (base[static_cast<std::size_t>(dimension)] != 0) {
      continue;
    }
    std::vector<std::vector<Integer>> triangle(
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
            triangle[static_cast<std::size_t>(level - 1)]
                    [static_cast<std::size_t>(index)]
            + triangle[static_cast<std::size_t>(level - 1)]
                      [static_cast<std::size_t>(index + 1)];
      }
    }
    for (int index = 0; index <= degree; ++index) {
      Exponent destination = base;
      destination[static_cast<std::size_t>(dimension)] = index;
      left.values[nd_grid_index(left.degrees, destination)] =
          triangle[static_cast<std::size_t>(index)][0]
          << (degree - index);
      right.values[nd_grid_index(right.degrees, destination)] =
          triangle[static_cast<std::size_t>(degree - index)]
                  [static_cast<std::size_t>(index)]
          << index;
    }
  }
  return {std::move(left), std::move(right)};
}

void nd_certify_integer_subdivision_with_corner(
    const NDIntegerBernsteinGrid& grid, const int depth,
    const int depth_limit, const Rational& corner_upper,
    const std::vector<Rational>& lower,
    const std::vector<Rational>& upper,
    NDCornerSubdivisionResult& result) {
  ++result.counts.nodes;
  result.counts.maximum_depth =
      std::max(result.counts.maximum_depth, depth);
  if (std::all_of(
          grid.values.begin(), grid.values.end(),
          [](const Integer& value) { return value >= 0; })) {
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
  auto [left, right] = nd_split_integer_grid(grid, dimension);
  std::vector<Rational> left_upper = upper;
  std::vector<Rational> right_lower = lower;
  const Rational middle =
      (lower[static_cast<std::size_t>(dimension)]
       + upper[static_cast<std::size_t>(dimension)]) / 2;
  left_upper[static_cast<std::size_t>(dimension)] = middle;
  right_lower[static_cast<std::size_t>(dimension)] = middle;
  nd_certify_integer_subdivision_with_corner(
      left, depth + 1, depth_limit, corner_upper,
      lower, left_upper, result);
  nd_certify_integer_subdivision_with_corner(
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

  BigPolynomial rational_h = nd_power(critical_a, 2);
  for (auto& [exponent, coefficient] : rational_h) {
    static_cast<void>(exponent);
    coefficient *= 4;
  }
  big_add_scaled(
      rational_h, big_multiply(critical_b, critical_c), 3);
  BigPolynomial rational_quadratic =
      big_multiply(critical_a, q2_payment);
  for (auto& [exponent, coefficient] : rational_quadratic) {
    static_cast<void>(exponent);
    coefficient *= 4;
  }
  big_add_scaled(
      rational_quadratic, nd_power(critical_b, 2), -1);
  BigPolynomial rational_certificate = big_multiply(
      rational_quadratic, nd_power(rational_h, 4));
  const BigPolynomial rational_first_correction = big_multiply(
      big_multiply(
          big_multiply(
              nd_power(critical_a, 4), critical_c),
          nd_power(critical_b, 3)),
      rational_h);
  big_add_scaled(
      rational_certificate, rational_first_correction, 32);
  const BigPolynomial rational_second_correction = big_multiply(
      big_multiply(
          nd_power(critical_a, 4), nd_power(critical_c, 2)),
      nd_power(critical_b, 4));
  big_add_scaled(
      rational_certificate, rational_second_correction, 144);
  const BigPolynomial compact_rational = canonicalize(
      compactify_positive_variable(rational_certificate, 0U));
  const NDBernsteinGrid rational_grid =
      nd_bernstein_grid(compact_rational, variables);
  const std::size_t rational_negative = static_cast<std::size_t>(
      std::count_if(
          rational_grid.values.begin(), rational_grid.values.end(),
          [](const Rational& value) { return value < 0; }));
  auto [rational_lower_half, rational_upper_half] =
      nd_split_grid(rational_grid, 0);
  auto [rational_middle_quarter, rational_top_quarter] =
      nd_split_grid(rational_upper_half, 0);
  static_cast<void>(rational_top_quarter);
  NDSubdivisionResult rational_moderate_subdivision;
  nd_certify_subdivision(
      rational_lower_half, 0, 24, rational_moderate_subdivision);
  nd_certify_subdivision(
      rational_middle_quarter, 0, 24,
      rational_moderate_subdivision);

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
  NDIntegerGridConversion blowup_integer_conversion =
      integer_bernstein_grid(blowup_grid);
  nd_certify_integer_subdivision_with_corner(
      blowup_integer_conversion.grid, 0,
      blowup_depth_limit, corner_upper,
      blowup_lower, blowup_upper, blowup_subdivision);

  const int rational_u_degree = rational_grid.degrees[1];
  const int rational_v_degree = rational_grid.degrees[2];
  int rational_maximum_w_power = 0;
  for (const auto& [exponent, coefficient] : compact_rational) {
    static_cast<void>(coefficient);
    rational_maximum_w_power = std::max(
        rational_maximum_w_power,
        exponent[0]
            + 2 * (rational_u_degree - exponent[1]));
  }
  const auto power_table = [&one](
      const BigPolynomial& base, const int maximum) {
    std::vector<BigPolynomial> powers;
    powers.reserve(static_cast<std::size_t>(maximum + 1));
    powers.push_back(one);
    for (int exponent = 1; exponent <= maximum; ++exponent) {
      powers.push_back(big_multiply(powers.back(), base));
    }
    return powers;
  };
  const std::vector<BigPolynomial> rational_w_powers =
      power_table(blowup_w, rational_maximum_w_power);
  const std::vector<BigPolynomial> rational_u_powers =
      power_table(blowup_u_numerator, rational_u_degree);
  const std::vector<BigPolynomial> rational_d_powers =
      power_table(blowup_denominator, rational_v_degree);
  const std::vector<BigPolynomial> rational_v_powers =
      power_table(blowup_v_numerator, rational_v_degree);
  std::map<std::pair<int, int>, BigPolynomial> rational_grouped_w;
  for (const auto& [exponent, coefficient] : compact_rational) {
    big_add_scaled(
        rational_grouped_w[{exponent[1], exponent[2]}],
        rational_w_powers[static_cast<std::size_t>(
            exponent[0]
                + 2 * (rational_u_degree - exponent[1]))],
        coefficient);
  }
  BigPolynomial rational_blowup_polynomial;
  for (const auto& [uv_exponents, w_polynomial] :
       rational_grouped_w) {
    BigPolynomial term = w_polynomial;
    term = big_multiply(
        term,
        rational_u_powers[static_cast<std::size_t>(
            uv_exponents.first)]);
    term = big_multiply(
        term,
        rational_d_powers[static_cast<std::size_t>(
            rational_v_degree - uv_exponents.second)]);
    term = big_multiply(
        term,
        rational_v_powers[static_cast<std::size_t>(
            uv_exponents.second)]);
    big_add_scaled(rational_blowup_polynomial, term, 1);
  }
  rational_blowup_polynomial =
      canonicalize(rational_blowup_polynomial);
  int rational_blowup_t_degree = 0;
  for (const auto& [exponent, coefficient] :
       rational_blowup_polynomial) {
    static_cast<void>(coefficient);
    rational_blowup_t_degree = std::max(
        rational_blowup_t_degree, exponent[0]);
  }
  BigPolynomial rational_scaled_blowup;
  for (const auto& [exponent, coefficient] :
       rational_blowup_polynomial) {
    Integer scale = 1;
    for (int power = exponent[0];
         power < rational_blowup_t_degree; ++power) {
      scale *= 4;
    }
    rational_scaled_blowup[exponent] += coefficient * scale;
  }
  rational_scaled_blowup = canonicalize(rational_scaled_blowup);
  constexpr int rational_blowup_t_factor = 14;
  BigPolynomial rational_blowup_quotient;
  for (const auto& [exponent, coefficient] :
       rational_scaled_blowup) {
    if (exponent[0] < rational_blowup_t_factor) {
      throw std::runtime_error(
          "wall (1,2,1) rational blowup t-factor mismatch");
    }
    Exponent quotient_exponent = exponent;
    quotient_exponent[0] -= rational_blowup_t_factor;
    rational_blowup_quotient[quotient_exponent] += coefficient;
  }
  rational_blowup_quotient =
      canonicalize(rational_blowup_quotient);
  const Integer rational_leading_scale =
      rational_blowup_quotient.at(leading_beta);
  if (
      rational_leading_scale <= 0
      || rational_blowup_quotient.at(leading_alpha2)
          != 64 * rational_leading_scale
      || rational_blowup_quotient.at(leading_alpha_beta)
          != 8 * rational_leading_scale
      || rational_blowup_quotient.at(leading_alpha)
          != -rational_leading_scale
      || 128 * rational_blowup_quotient.at(leading_constant)
          != 3 * rational_leading_scale
  ) {
    throw std::runtime_error(
        "wall (1,2,1) rational blowup Newton face mismatch");
  }
  const auto compute_rational_corner_budget = [
      &rational_blowup_quotient, &rational_leading_scale
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
    for (const auto& [exponent, coefficient] :
         rational_blowup_quotient) {
      if (coefficient >= 0) {
        continue;
      }
      ++budget.negative_monomials;
      const int t_exponent = exponent[0];
      const int alpha_exponent = exponent[1];
      const int beta_exponent = exponent[2];
      const Rational magnitude(
          -coefficient, rational_leading_scale);
      if (
          alpha_exponent == 1 && beta_exponent == 0
          && t_exponent >= 2
      ) {
        linear_alpha +=
            magnitude * corner_power(t_exponent - 2);
      } else if (
          alpha_exponent >= 1 && beta_exponent >= 1
          && t_exponent >= 1
      ) {
        cross_loss += magnitude * corner_power(
            t_exponent - 1 + alpha_exponent - 1
            + beta_exponent - 1);
      } else if (alpha_exponent >= 2) {
        alpha_loss += magnitude * corner_power(
            t_exponent + alpha_exponent - 2
            + beta_exponent);
      } else if (beta_exponent >= 1 && t_exponent >= 2) {
        beta_loss += magnitude * corner_power(
            t_exponent - 2 + alpha_exponent
            + beta_exponent - 1);
      } else if (
          alpha_exponent == 0 && beta_exponent == 0
          && t_exponent >= 4
      ) {
        constant_loss +=
            magnitude * corner_power(t_exponent - 4);
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
  int rational_corner_power_two = 4;
  Rational rational_corner_upper(1, 16);
  CornerBudget rational_corner_budget =
      compute_rational_corner_budget(rational_corner_upper);
  for (; rational_corner_power_two <= 64;
       ++rational_corner_power_two) {
    Integer denominator = 1;
    denominator <<= rational_corner_power_two;
    rational_corner_upper = Rational(Integer(1), denominator);
    rational_corner_budget =
        compute_rational_corner_budget(rational_corner_upper);
    if (
        rational_corner_budget.unclassified_negative == 0U
        && rational_corner_budget.alpha_reserve >= 0
        && rational_corner_budget.cross_reserve >= 0
        && rational_corner_budget.beta_reserve >= 0
        && rational_corner_budget.constant_reserve >= 0
        && rational_corner_budget.discriminant_margin >= 0
    ) {
      break;
    }
  }
  if (rational_corner_power_two > 64) {
    throw std::runtime_error(
        "wall (1,2,1) rational blowup corner budget failed");
  }
  const NDBernsteinGrid rational_blowup_grid =
      nd_bernstein_grid(rational_blowup_quotient, variables);
  const std::size_t rational_blowup_negative =
      static_cast<std::size_t>(std::count_if(
          rational_blowup_grid.values.begin(),
          rational_blowup_grid.values.end(),
          [](const Rational& value) { return value < 0; }));
  NDIntegerGridConversion rational_integer_conversion =
      integer_bernstein_grid(rational_blowup_grid);
  NDCornerSubdivisionResult rational_blowup_subdivision;
  const int rational_blowup_depth_limit =
      3 * rational_corner_power_two + 24;
  nd_certify_integer_subdivision_with_corner(
      rational_integer_conversion.grid, 0,
      rational_blowup_depth_limit, rational_corner_upper,
      blowup_lower, blowup_upper, rational_blowup_subdivision);
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
      || rational_grid.degrees != std::vector<int>{44, 22, 7}
      || rational_grid.values.size() != 8280U
      || rational_negative != 897U
      || rational_moderate_subdivision.nodes != 18U
      || rational_moderate_subdivision.leaves != 10U
      || rational_moderate_subdivision.unresolved != 0U
      || rational_moderate_subdivision.maximum_depth != 6
      || rational_maximum_w_power != 88
      || rational_grouped_w.size() != 133U
      || rational_blowup_t_degree != 130
      || rational_blowup_grid.degrees
          != std::vector<int>{116, 33, 7}
      || rational_blowup_quotient.size() != 26201U
      || rational_blowup_grid.values.size() != 31824U
      || rational_blowup_negative != 115U
      || rational_corner_budget.negative_monomials != 13070U
      || rational_corner_budget.unclassified_negative != 0U
      || rational_corner_power_two != 6
      || rational_corner_budget.alpha_reserve < 0
      || rational_corner_budget.cross_reserve < 0
      || rational_corner_budget.beta_reserve < 0
      || rational_corner_budget.constant_reserve < 0
      || rational_corner_budget.discriminant_margin < 0
      || rational_integer_conversion.distinct_denominators != 14594U
      || rational_integer_conversion.common_denominator_bits != 206U
      || rational_blowup_subdivision.counts.nodes != 987U
      || rational_blowup_subdivision.counts.leaves != 493U
      || rational_blowup_subdivision.corner_leaves != 1U
      || rational_blowup_subdivision.counts.unresolved != 0U
      || rational_blowup_subdivision.counts.maximum_depth != 38
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
      << " rational_degrees=(44,22,7)"
      << " rational_coefficients=8280"
      << " rational_initial_negative=897"
      << " rational_moderate_nodes=18"
      << " rational_moderate_leaves=10"
      << " rational_moderate_unresolved=0"
      << " rational_moderate_depth=6"
      << " rational_blowup_degrees=(116,33,7)"
      << " rational_blowup_terms=26201"
      << " rational_blowup_coefficients=31824"
      << " rational_blowup_initial_negative=115"
      << " rational_blowup_t_factor=14"
      << " rational_corner=1/64"
      << " rational_corner_negative_monomials=13070"
      << " rational_corner_unclassified=0"
      << " rational_corner_reserves_nonnegative=1"
      << " rational_corner_discriminant_nonnegative=1"
      << " rational_integer_denominators=14594"
      << " rational_integer_common_denominator_bits=206"
      << " rational_blowup_nodes=987"
      << " rational_blowup_leaves=493"
      << " rational_blowup_corner_leaves=1"
      << " rational_blowup_unresolved=0"
      << " rational_blowup_depth=38"
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

int replay_wall_121_rational_support_three() {
  constexpr int variables = 2;
  constexpr int support = 3;
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial r = nd_variable(variables, 0);
  const BigPolynomial u = nd_variable(variables, 1);
  const BigPolynomial zero;
  std::vector<BigPolynomial> profile(8, zero);
  profile[1] = one;
  profile[2] = r;
  profile[3] = big_multiply(nd_power(r, 2), u);
  std::vector<BigPolynomial> g(7, zero);
  std::vector<BigPolynomial> h(7, zero);
  std::vector<BigPolynomial> k(7, zero);
  for (int index = 1; index <= support + 1; ++index) {
    g[static_cast<std::size_t>(index)] = big_multiply(
        profile[static_cast<std::size_t>(index)],
        profile[static_cast<std::size_t>(index)]);
    big_add_scaled(
        g[static_cast<std::size_t>(index)],
        big_multiply(
            profile[static_cast<std::size_t>(index - 1)],
            profile[static_cast<std::size_t>(index + 1)]),
        -1);
  }
  for (int index = 1; index <= support + 2; ++index) {
    h[static_cast<std::size_t>(index)] = big_multiply(
        profile[static_cast<std::size_t>(index)],
        profile[static_cast<std::size_t>(index - 1)]);
    if (index >= 2) {
      big_add_scaled(
          h[static_cast<std::size_t>(index)],
          big_multiply(
              profile[static_cast<std::size_t>(index - 2)],
              profile[static_cast<std::size_t>(index + 1)]),
          -1);
    }
  }
  for (int index = 1; index <= support + 1; ++index) {
    k[static_cast<std::size_t>(index)] =
        g[static_cast<std::size_t>(index)];
    big_add_scaled(
        k[static_cast<std::size_t>(index)],
        h[static_cast<std::size_t>(index)], 1);
    big_add_scaled(
        k[static_cast<std::size_t>(index)],
        h[static_cast<std::size_t>(index + 1)], 1);
  }
  BigPolynomial c_zero = r;
  for (int index = 1; index <= support; ++index) {
    BigPolynomial g_difference =
        g[static_cast<std::size_t>(index)];
    big_add_scaled(
        g_difference,
        g[static_cast<std::size_t>(index + 1)], -1);
    BigPolynomial k_difference =
        k[static_cast<std::size_t>(index)];
    big_add_scaled(
        k_difference,
        k[static_cast<std::size_t>(index + 1)], -1);
    big_add_scaled(
        c_zero,
        big_multiply(g_difference, k_difference), 1);
  }
  const BigPolynomial c = profile[3];
  BigPolynomial a_coefficient = one;
  big_add_scaled(a_coefficient, nd_power(r, 2), 1);
  BigPolynomial one_plus_r = one;
  big_add_scaled(one_plus_r, r, 1);
  big_add_scaled(
      a_coefficient,
      big_multiply(one_plus_r, c), 1);
  BigPolynomial c_coefficient = one_plus_r;
  big_add_scaled(c_coefficient, c, 1);
  BigPolynomial b_coefficient;
  big_add_scaled(b_coefficient, r, 3);
  big_add_scaled(b_coefficient, nd_power(r, 2), 2);
  big_add_scaled(
      b_coefficient, big_multiply(r, c), 1);
  big_add_scaled(b_coefficient, nd_power(r, 3), -1);
  big_add_scaled(b_coefficient, nd_power(c, 2), -1);
  big_add_scaled(b_coefficient, nd_power(c, 3), -1);
  BigPolynomial rational_h = nd_power(a_coefficient, 2);
  for (auto& [exponent, coefficient] : rational_h) {
    static_cast<void>(exponent);
    coefficient *= 4;
  }
  big_add_scaled(
      rational_h,
      big_multiply(b_coefficient, c_coefficient), 3);
  BigPolynomial rational_quadratic =
      big_multiply(a_coefficient, c_zero);
  for (auto& [exponent, coefficient] : rational_quadratic) {
    static_cast<void>(exponent);
    coefficient *= 4;
  }
  big_add_scaled(
      rational_quadratic, nd_power(b_coefficient, 2), -1);
  BigPolynomial rational_certificate = big_multiply(
      rational_quadratic, nd_power(rational_h, 4));
  const BigPolynomial first_correction = big_multiply(
      big_multiply(
          big_multiply(
              nd_power(a_coefficient, 4), c_coefficient),
          nd_power(b_coefficient, 3)),
      rational_h);
  big_add_scaled(rational_certificate, first_correction, 32);
  const BigPolynomial second_correction = big_multiply(
      big_multiply(
          nd_power(a_coefficient, 4),
          nd_power(c_coefficient, 2)),
      nd_power(b_coefficient, 4));
  big_add_scaled(rational_certificate, second_correction, 144);
  const BigPolynomial compact =
      compactify_positive_variable(rational_certificate, 0U);
  BigPolynomial canonical;
  for (const auto& [exponent, coefficient] : compact) {
    if (coefficient != 0) {
      canonical[exponent] = coefficient;
    }
  }
  const NDBernsteinGrid grid =
      nd_bernstein_grid(canonical, variables);
  const std::size_t negative = static_cast<std::size_t>(
      std::count_if(
          grid.values.begin(), grid.values.end(),
          [](const Rational& value) { return value < 0; }));
  auto [lower_half, upper_half] = nd_split_grid(grid, 0);
  auto [middle_quarter, top_quarter] =
      nd_split_grid(upper_half, 0);
  auto [upper_eighth, top_eighth] =
      nd_split_grid(top_quarter, 0);
  static_cast<void>(top_eighth);
  NDSubdivisionResult subdivision;
  nd_certify_subdivision(
      lower_half, 0, 30, subdivision);
  nd_certify_subdivision(
      middle_quarter, 0, 30, subdivision);
  nd_certify_subdivision(
      upper_eighth, 0, 30, subdivision);
  const Rational branch_upper_at_four =
      3 + Rational(9, 4) * 4 - 16;
  const Rational branch_upper_derivative_at_four =
      Rational(9, 4) - 8;
  if (
      grid.degrees != std::vector<int>{44, 22}
      || grid.values.size() != 1035U
      || negative != 187U
      || subdivision.nodes != 29U
      || subdivision.leaves != 16U
      || subdivision.unresolved != 0U
      || subdivision.maximum_depth != 6
      || branch_upper_at_four != -4
      || branch_upper_derivative_at_four != Rational(-23, 4)
  ) {
    throw std::runtime_error(
        "wall (1,2,1) support-three rational certificate mismatch");
  }
  std::cout
      << "SU2_WALL_121_RATIONAL_SUPPORT_THREE"
      << " support<=3"
      << " branch_bound=R<4"
      << " certified_box=R<=7"
      << " degrees=(44,22)"
      << " bernstein_coefficients=1035"
      << " initial_negative=187"
      << " subdivision_nodes=29"
      << " subdivision_leaves=16"
      << " subdivision_unresolved=0"
      << " subdivision_depth=6"
      << " result=PASS_EXACT\n";
  return EXIT_SUCCESS;
}

BigPolynomial wall_121_terminal_rational_certificate(
    const std::vector<BigPolynomial>& profile, const int support) {
  if (
      support < 4
      || profile.size() < static_cast<std::size_t>(support + 4)
  ) {
    throw std::invalid_argument(
        "wall terminal rational certificate needs support at least four");
  }
  const int variables =
      static_cast<int>(profile[1].begin()->first.size());
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial zero;
  std::vector<BigPolynomial> g(
      static_cast<std::size_t>(support + 4), zero);
  std::vector<BigPolynomial> h(
      static_cast<std::size_t>(support + 4), zero);
  std::vector<BigPolynomial> k(
      static_cast<std::size_t>(support + 4), zero);
  for (int index = 1; index <= support + 1; ++index) {
    g[static_cast<std::size_t>(index)] = big_multiply(
        profile[static_cast<std::size_t>(index)],
        profile[static_cast<std::size_t>(index)]);
    big_add_scaled(
        g[static_cast<std::size_t>(index)],
        big_multiply(
            profile[static_cast<std::size_t>(index - 1)],
            profile[static_cast<std::size_t>(index + 1)]),
        -1);
  }
  for (int index = 1; index <= support + 2; ++index) {
    h[static_cast<std::size_t>(index)] = big_multiply(
        profile[static_cast<std::size_t>(index)],
        profile[static_cast<std::size_t>(index - 1)]);
    if (index >= 2) {
      big_add_scaled(
          h[static_cast<std::size_t>(index)],
          big_multiply(
              profile[static_cast<std::size_t>(index - 2)],
              profile[static_cast<std::size_t>(index + 1)]),
          -1);
    }
  }
  for (int index = 1; index <= support + 1; ++index) {
    k[static_cast<std::size_t>(index)] =
        g[static_cast<std::size_t>(index)];
    big_add_scaled(
        k[static_cast<std::size_t>(index)],
        h[static_cast<std::size_t>(index)], 1);
    big_add_scaled(
        k[static_cast<std::size_t>(index)],
        h[static_cast<std::size_t>(index + 1)], 1);
  }
  BigPolynomial c_zero = big_multiply(
      nd_power(profile[1], 3), profile[2]);
  for (int index = 1; index <= support; ++index) {
    BigPolynomial g_difference =
        g[static_cast<std::size_t>(index)];
    big_add_scaled(
        g_difference,
        g[static_cast<std::size_t>(index + 1)], -1);
    BigPolynomial k_difference =
        k[static_cast<std::size_t>(index)];
    big_add_scaled(
        k_difference,
        k[static_cast<std::size_t>(index + 1)], -1);
    big_add_scaled(
        c_zero,
        big_multiply(g_difference, k_difference), 1);
  }
  const BigPolynomial& a = profile[1];
  const BigPolynomial& b = profile[2];
  const BigPolynomial& c = profile[3];
  const BigPolynomial& d = profile[4];
  BigPolynomial a_coefficient = nd_power(a, 2);
  big_add_scaled(a_coefficient, nd_power(b, 2), 1);
  BigPolynomial a_plus_b = a;
  big_add_scaled(a_plus_b, b, 1);
  big_add_scaled(
      a_coefficient, big_multiply(a_plus_b, c), 1);
  BigPolynomial c_coefficient = a_plus_b;
  big_add_scaled(c_coefficient, c, 1);
  BigPolynomial b_coefficient;
  big_add_scaled(
      b_coefficient,
      big_multiply(big_multiply(nd_power(a, 2), b), one), 3);
  big_add_scaled(
      b_coefficient,
      big_multiply(big_multiply(a, nd_power(b, 2)), one), 2);
  big_add_scaled(
      b_coefficient, big_multiply(big_multiply(a, b), c), 1);
  big_add_scaled(
      b_coefficient, big_multiply(big_multiply(a, b), d), 1);
  big_add_scaled(
      b_coefficient, big_multiply(big_multiply(b, c), d), 1);
  big_add_scaled(b_coefficient, nd_power(b, 3), -1);
  big_add_scaled(
      b_coefficient, big_multiply(a, nd_power(c, 2)), -1);
  big_add_scaled(b_coefficient, nd_power(c, 3), -1);
  BigPolynomial rational_h = nd_power(a_coefficient, 2);
  for (auto& [exponent, coefficient] : rational_h) {
    static_cast<void>(exponent);
    coefficient *= 4;
  }
  big_add_scaled(
      rational_h,
      big_multiply(b_coefficient, c_coefficient), 3);
  BigPolynomial quadratic =
      big_multiply(a_coefficient, c_zero);
  for (auto& [exponent, coefficient] : quadratic) {
    static_cast<void>(exponent);
    coefficient *= 4;
  }
  big_add_scaled(quadratic, nd_power(b_coefficient, 2), -1);
  BigPolynomial result =
      big_multiply(quadratic, nd_power(rational_h, 4));
  const BigPolynomial first_correction = big_multiply(
      big_multiply(
          big_multiply(
              nd_power(a_coefficient, 4), c_coefficient),
          nd_power(b_coefficient, 3)),
      rational_h);
  big_add_scaled(result, first_correction, 32);
  const BigPolynomial second_correction = big_multiply(
      big_multiply(
          nd_power(a_coefficient, 4),
          nd_power(c_coefficient, 2)),
      nd_power(b_coefficient, 4));
  big_add_scaled(result, second_correction, 144);
  return result;
}

BigPolynomial wall_121_terminal_current(
    const std::vector<BigPolynomial>& profile, const int support);

// This is the cleared second fixed-point payment for the wall (1,2,1)
// cubic.  Unlike wall_121_terminal_rational_certificate(), it is only a
// diagnostic for the remaining arbitrary-support lemma: nonnegativity would
// imply the exact critical payment, but it is not equivalent to it.
BigPolynomial wall_121_terminal_second_iterate_certificate(
    const std::vector<BigPolynomial>& profile, const int support) {
  if (
      support < 1
      || profile.size() < static_cast<std::size_t>(support + 4)
  ) {
    throw std::invalid_argument(
        "wall terminal second-iterate certificate needs positive support");
  }
  const int variables =
      static_cast<int>(profile[1].begin()->first.size());
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial& a = profile[1];
  const BigPolynomial& b = profile[2];
  const BigPolynomial& c = profile[3];
  const BigPolynomial& d = profile[4];

  BigPolynomial a_coefficient = nd_power(a, 2);
  big_add_scaled(a_coefficient, nd_power(b, 2), 1);
  BigPolynomial a_plus_b = a;
  big_add_scaled(a_plus_b, b, 1);
  big_add_scaled(
      a_coefficient, big_multiply(a_plus_b, c), 1);

  BigPolynomial c_coefficient = a_plus_b;
  big_add_scaled(c_coefficient, c, 1);

  // b_coefficient is B=-C_1 in the critical cubic C y^3+A y^2-B y+C_0.
  BigPolynomial b_coefficient;
  big_add_scaled(
      b_coefficient,
      big_multiply(big_multiply(nd_power(a, 2), b), one), 3);
  big_add_scaled(
      b_coefficient,
      big_multiply(big_multiply(a, nd_power(b, 2)), one), 2);
  big_add_scaled(
      b_coefficient, big_multiply(big_multiply(a, b), c), 1);
  big_add_scaled(
      b_coefficient, big_multiply(big_multiply(a, b), d), 1);
  big_add_scaled(
      b_coefficient, big_multiply(big_multiply(b, c), d), 1);
  big_add_scaled(b_coefficient, nd_power(b, 3), -1);
  big_add_scaled(
      b_coefficient, big_multiply(a, nd_power(c, 2)), -1);
  big_add_scaled(b_coefficient, nd_power(c, 3), -1);

  const BigPolynomial c_zero =
      wall_121_terminal_current(profile, support);
  BigPolynomial j = nd_power(a_coefficient, 2);
  for (auto& [exponent, coefficient] : j) {
    static_cast<void>(exponent);
    coefficient *= 2;
  }
  big_add_scaled(
      j, big_multiply(b_coefficient, c_coefficient), 3);

  BigPolynomial h = nd_power(a_coefficient, 2);
  for (auto& [exponent, coefficient] : h) {
    static_cast<void>(exponent);
    coefficient *= 4;
  }
  big_add_scaled(
      h, big_multiply(b_coefficient, c_coefficient), 3);

  BigPolynomial k = big_multiply(nd_power(a_coefficient, 2), j);
  for (auto& [exponent, coefficient] : k) {
    static_cast<void>(exponent);
    coefficient *= 8;
  }
  big_add_scaled(
      k,
      big_multiply(
          big_multiply(b_coefficient, c_coefficient), h), 3);

  BigPolynomial n = big_multiply(
      big_multiply(a_coefficient, b_coefficient), j);
  for (auto& [exponent, coefficient] : n) {
    static_cast<void>(exponent);
    coefficient *= 4;
  }

  BigPolynomial quadratic = big_multiply(a_coefficient, c_zero);
  for (auto& [exponent, coefficient] : quadratic) {
    static_cast<void>(exponent);
    coefficient *= 4;
  }
  big_add_scaled(quadratic, nd_power(b_coefficient, 2), -1);

  BigPolynomial result = big_multiply(quadratic, nd_power(k, 4));
  const BigPolynomial first_correction = big_multiply(
      big_multiply(
          big_multiply(a_coefficient, c_coefficient), nd_power(n, 3)),
      k);
  big_add_scaled(result, first_correction, 4);
  const BigPolynomial second_correction = big_multiply(
      big_multiply(nd_power(c_coefficient, 2), nd_power(n, 4)), one);
  big_add_scaled(result, second_correction, 9);
  return result;
}

BigPolynomial wall_121_terminal_current(
    const std::vector<BigPolynomial>& profile, const int support) {
  if (
      support < 1
      || profile.size() < static_cast<std::size_t>(support + 4)
  ) {
    throw std::invalid_argument(
        "wall terminal current needs a positive finite support");
  }
  const BigPolynomial zero;
  std::vector<BigPolynomial> g(
      static_cast<std::size_t>(support + 4), zero);
  std::vector<BigPolynomial> h(
      static_cast<std::size_t>(support + 4), zero);
  std::vector<BigPolynomial> k(
      static_cast<std::size_t>(support + 4), zero);
  for (int index = 1; index <= support + 1; ++index) {
    g[static_cast<std::size_t>(index)] = big_multiply(
        profile[static_cast<std::size_t>(index)],
        profile[static_cast<std::size_t>(index)]);
    big_add_scaled(
        g[static_cast<std::size_t>(index)],
        big_multiply(
            profile[static_cast<std::size_t>(index - 1)],
            profile[static_cast<std::size_t>(index + 1)]),
        -1);
  }
  for (int index = 1; index <= support + 2; ++index) {
    h[static_cast<std::size_t>(index)] = big_multiply(
        profile[static_cast<std::size_t>(index)],
        profile[static_cast<std::size_t>(index - 1)]);
    if (index >= 2) {
      big_add_scaled(
          h[static_cast<std::size_t>(index)],
          big_multiply(
              profile[static_cast<std::size_t>(index - 2)],
              profile[static_cast<std::size_t>(index + 1)]),
          -1);
    }
  }
  for (int index = 1; index <= support + 1; ++index) {
    k[static_cast<std::size_t>(index)] =
        g[static_cast<std::size_t>(index)];
    big_add_scaled(
        k[static_cast<std::size_t>(index)],
        h[static_cast<std::size_t>(index)], 1);
    big_add_scaled(
        k[static_cast<std::size_t>(index)],
        h[static_cast<std::size_t>(index + 1)], 1);
  }
  BigPolynomial result = big_multiply(
      nd_power(profile[1], 3), profile[2]);
  for (int index = 1; index <= support; ++index) {
    BigPolynomial g_difference = g[static_cast<std::size_t>(index)];
    big_add_scaled(
        g_difference, g[static_cast<std::size_t>(index + 1)], -1);
    BigPolynomial k_difference = k[static_cast<std::size_t>(index)];
    big_add_scaled(
        k_difference, k[static_cast<std::size_t>(index + 1)], -1);
    big_add_scaled(result, big_multiply(g_difference, k_difference), 1);
  }
  return result;
}

BigPolynomial wall_121_second_iterate_support_five_compact() {
  constexpr int variables = 4;
  constexpr int support = 5;
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial r = nd_variable(variables, 0);
  const BigPolynomial u = nd_variable(variables, 1);
  const BigPolynomial v = nd_variable(variables, 2);
  const BigPolynomial z = nd_variable(variables, 3);
  const BigPolynomial zero;
  std::vector<BigPolynomial> profile(
      static_cast<std::size_t>(support + 4), zero);
  profile[1] = one;
  profile[2] = r;
  profile[3] = big_multiply(nd_power(r, 2), u);
  profile[4] = big_multiply(
      big_multiply(nd_power(r, 3), nd_power(u, 2)), v);
  profile[5] = big_multiply(
      big_multiply(
          big_multiply(nd_power(r, 4), nd_power(u, 3)),
          nd_power(v, 2)),
      z);
  const BigPolynomial target =
      wall_121_terminal_second_iterate_certificate(profile, support);
  const BigPolynomial compact_raw =
      compactify_positive_variable(target, 0U);
  BigPolynomial compact;
  for (const auto& [exponent, coefficient] : compact_raw) {
    if (coefficient != 0) {
      compact[exponent] = coefficient;
    }
  }
  return compact;
}

NDBernsteinGrid wall_121_second_iterate_support_five_grid(
    std::size_t& terms
) {
  const BigPolynomial compact =
      wall_121_second_iterate_support_five_compact();
  terms = compact.size();
  return nd_bernstein_grid(compact, 4);
}

struct Wall121LargeRatioBlowupResult {
  int maximum_w_power = 0;
  std::size_t grouped_blocks = 0U;
  int t_degree = 0;
  int t_factor = 0;
  std::vector<int> degrees;
  std::size_t terms = 0U;
  std::size_t coefficients = 0U;
  std::size_t negative_coefficients = 0U;
  std::size_t t_zero_terms = 0U;
  std::size_t t_zero_negative = 0U;
  std::size_t distinct_denominators = 0U;
  std::size_t common_denominator_bits = 0U;
  NDCornerSubdivisionResult subdivision;
};

Wall121LargeRatioBlowupResult wall_121_large_ratio_blowup(
    const BigPolynomial& compact, const NDBernsteinGrid& compact_grid) {
  constexpr int variables = 3;
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial t = nd_variable(variables, 0);
  const BigPolynomial alpha = nd_variable(variables, 1);
  const BigPolynomial beta = nd_variable(variables, 2);
  BigPolynomial w = one;
  big_add_scaled(w, t, -1);
  BigPolynomial u_numerator = nd_power(w, 2);
  big_add_scaled(
      u_numerator,
      big_multiply(big_multiply(t, w), alpha), -2);
  big_add_scaled(
      u_numerator,
      big_multiply(nd_power(t, 2), alpha), -3);
  BigPolynomial t2_plus_u = nd_power(t, 2);
  big_add_scaled(t2_plus_u, u_numerator, 1);
  const BigPolynomial denominator = big_multiply(
      nd_power(u_numerator, 2), t2_plus_u);
  BigPolynomial two_w_plus_three_t = w;
  for (auto& [exponent, coefficient] : two_w_plus_three_t) {
    static_cast<void>(exponent);
    coefficient *= 2;
  }
  big_add_scaled(two_w_plus_three_t, t, 3);
  BigPolynomial one_minus_alpha = one;
  big_add_scaled(one_minus_alpha, alpha, -1);
  const BigPolynomial drop_numerator = big_multiply(
      big_multiply(
          big_multiply(
              big_multiply(beta, nd_power(t, 4)), w),
          two_w_plus_three_t),
      one_minus_alpha);
  BigPolynomial v_numerator = denominator;
  big_add_scaled(v_numerator, drop_numerator, -1);
  const int u_degree = compact_grid.degrees[1];
  const int v_degree = compact_grid.degrees[2];
  int maximum_w_power = 0;
  for (const auto& [exponent, coefficient] : compact) {
    static_cast<void>(coefficient);
    maximum_w_power = std::max(
        maximum_w_power,
        exponent[0] + 2 * (u_degree - exponent[1]));
  }
  const auto power_table = [&one](
      const BigPolynomial& base, const int maximum) {
    std::vector<BigPolynomial> powers;
    powers.reserve(static_cast<std::size_t>(maximum + 1));
    powers.push_back(one);
    for (int exponent = 1; exponent <= maximum; ++exponent) {
      powers.push_back(big_multiply(powers.back(), base));
    }
    return powers;
  };
  const std::vector<BigPolynomial> w_powers =
      power_table(w, maximum_w_power);
  const std::vector<BigPolynomial> u_powers =
      power_table(u_numerator, u_degree);
  const std::vector<BigPolynomial> d_powers =
      power_table(denominator, v_degree);
  const std::vector<BigPolynomial> v_powers =
      power_table(v_numerator, v_degree);
  std::map<std::pair<int, int>, BigPolynomial> grouped_w;
  for (const auto& [exponent, coefficient] : compact) {
    big_add_scaled(
        grouped_w[{exponent[1], exponent[2]}],
        w_powers[static_cast<std::size_t>(
            exponent[0] + 2 * (u_degree - exponent[1]))],
        coefficient);
  }
  BigPolynomial blowup;
  for (const auto& [uv_exponents, w_polynomial] : grouped_w) {
    BigPolynomial term = w_polynomial;
    term = big_multiply(
        term,
        u_powers[static_cast<std::size_t>(uv_exponents.first)]);
    term = big_multiply(
        term,
        d_powers[static_cast<std::size_t>(
            v_degree - uv_exponents.second)]);
    term = big_multiply(
        term,
        v_powers[static_cast<std::size_t>(uv_exponents.second)]);
    big_add_scaled(blowup, term, 1);
  }
  BigPolynomial canonical_blowup;
  for (const auto& [exponent, coefficient] : blowup) {
    if (coefficient != 0) {
      canonical_blowup[exponent] = coefficient;
    }
  }
  int t_degree = 0;
  int t_factor = std::numeric_limits<int>::max();
  for (const auto& [exponent, coefficient] : canonical_blowup) {
    static_cast<void>(coefficient);
    t_degree = std::max(t_degree, exponent[0]);
    t_factor = std::min(t_factor, exponent[0]);
  }
  BigPolynomial quotient;
  for (const auto& [exponent, coefficient] : canonical_blowup) {
    Integer scale = 1;
    for (int power = exponent[0]; power < t_degree; ++power) {
      scale *= 4;
    }
    Exponent quotient_exponent = exponent;
    quotient_exponent[0] -= t_factor;
    quotient[quotient_exponent] += coefficient * scale;
  }
  BigPolynomial canonical_quotient;
  for (const auto& [exponent, coefficient] : quotient) {
    if (coefficient != 0) {
      canonical_quotient[exponent] = coefficient;
    }
  }
  const NDBernsteinGrid quotient_grid =
      nd_bernstein_grid(canonical_quotient, variables);
  const std::size_t negative = static_cast<std::size_t>(
      std::count_if(
          quotient_grid.values.begin(), quotient_grid.values.end(),
          [](const Rational& value) { return value < 0; }));
  std::size_t t_zero_terms = 0U;
  std::size_t t_zero_negative = 0U;
  for (const auto& [exponent, coefficient] : canonical_quotient) {
    if (exponent[0] == 0) {
      ++t_zero_terms;
      if (coefficient < 0) {
        ++t_zero_negative;
      }
    }
  }
  NDIntegerGridConversion conversion =
      integer_bernstein_grid(quotient_grid);
  NDCornerSubdivisionResult subdivision;
  const std::vector<Rational> lower(
      static_cast<std::size_t>(variables), Rational(0));
  const std::vector<Rational> upper(
      static_cast<std::size_t>(variables), Rational(1));
  nd_certify_integer_subdivision_with_corner(
      conversion.grid, 0, 42, Rational(0),
      lower, upper, subdivision);
  return {
      maximum_w_power, grouped_w.size(), t_degree, t_factor,
      quotient_grid.degrees, canonical_quotient.size(),
      quotient_grid.values.size(), negative,
      t_zero_terms, t_zero_negative,
      conversion.distinct_denominators,
      conversion.common_denominator_bits,
      std::move(subdivision)};
}

struct Wall121SupportFiveLargeRatioDiagnostic {
  int maximum_w_power = 0;
  std::size_t grouped_blocks = 0U;
  int t_degree = 0;
  int t_factor = 0;
  std::vector<int> degrees;
  std::size_t terms = 0U;
  std::size_t coefficients = 0U;
  std::size_t negative_coefficients = 0U;
  std::size_t t_zero_terms = 0U;
  std::size_t t_zero_negative = 0U;
};

struct Wall121SupportFiveLargeRatioQuotient {
  int maximum_w_power = 0;
  std::size_t grouped_blocks = 0U;
  int t_degree = 0;
  int t_factor = 0;
  BigPolynomial quotient;
};

Wall121SupportFiveLargeRatioQuotient
wall_121_support_five_large_ratio_quotient() {
  constexpr int variables = 4;
  const BigPolynomial compact =
      wall_121_second_iterate_support_five_compact();
  const NDBernsteinGrid compact_grid =
      nd_bernstein_grid(compact, variables);
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial t = nd_variable(variables, 0);
  const BigPolynomial alpha = nd_variable(variables, 1);
  const BigPolynomial beta = nd_variable(variables, 2);
  const BigPolynomial z = nd_variable(variables, 3);
  BigPolynomial w = one;
  big_add_scaled(w, t, -1);
  BigPolynomial u_numerator = nd_power(w, 2);
  big_add_scaled(
      u_numerator,
      big_multiply(big_multiply(t, w), alpha), -2);
  big_add_scaled(
      u_numerator,
      big_multiply(nd_power(t, 2), alpha), -3);
  BigPolynomial t2_plus_u = nd_power(t, 2);
  big_add_scaled(t2_plus_u, u_numerator, 1);
  const BigPolynomial denominator = big_multiply(
      nd_power(u_numerator, 2), t2_plus_u);
  BigPolynomial two_w_plus_three_t = w;
  for (auto& [exponent, coefficient] : two_w_plus_three_t) {
    static_cast<void>(exponent);
    coefficient *= 2;
  }
  big_add_scaled(two_w_plus_three_t, t, 3);
  BigPolynomial one_minus_alpha = one;
  big_add_scaled(one_minus_alpha, alpha, -1);
  const BigPolynomial drop_numerator = big_multiply(
      big_multiply(
          big_multiply(
              big_multiply(beta, nd_power(t, 4)), w),
          two_w_plus_three_t),
      one_minus_alpha);
  BigPolynomial v_numerator = denominator;
  big_add_scaled(v_numerator, drop_numerator, -1);

  const int u_degree = compact_grid.degrees[1];
  const int v_degree = compact_grid.degrees[2];
  const int z_degree = compact_grid.degrees[3];
  int maximum_w_power = 0;
  for (const auto& [exponent, coefficient] : compact) {
    static_cast<void>(coefficient);
    maximum_w_power = std::max(
        maximum_w_power,
        exponent[0] + 2 * (u_degree - exponent[1]));
  }
  const auto power_table = [&one](
      const BigPolynomial& base, const int maximum) {
    std::vector<BigPolynomial> powers;
    powers.reserve(static_cast<std::size_t>(maximum + 1));
    powers.push_back(one);
    for (int exponent = 1; exponent <= maximum; ++exponent) {
      powers.push_back(big_multiply(powers.back(), base));
    }
    return powers;
  };
  const std::vector<BigPolynomial> w_powers =
      power_table(w, maximum_w_power);
  const std::vector<BigPolynomial> u_powers =
      power_table(u_numerator, u_degree);
  const std::vector<BigPolynomial> d_powers =
      power_table(denominator, v_degree);
  const std::vector<BigPolynomial> v_powers =
      power_table(v_numerator, v_degree);
  const std::vector<BigPolynomial> z_powers =
      power_table(z, z_degree);
  std::map<Exponent, BigPolynomial> grouped_w;
  for (const auto& [exponent, coefficient] : compact) {
    const Exponent suffix{exponent[1], exponent[2], exponent[3]};
    big_add_scaled(
        grouped_w[suffix],
        w_powers[static_cast<std::size_t>(
            exponent[0] + 2 * (u_degree - exponent[1]))],
        coefficient);
  }
  BigPolynomial blowup;
  for (const auto& [suffix, w_polynomial] : grouped_w) {
    BigPolynomial term = w_polynomial;
    term = big_multiply(
        term, u_powers[static_cast<std::size_t>(suffix[0])]);
    term = big_multiply(
        term,
        d_powers[static_cast<std::size_t>(v_degree - suffix[1])]);
    term = big_multiply(
        term, v_powers[static_cast<std::size_t>(suffix[1])]);
    term = big_multiply(
        term, z_powers[static_cast<std::size_t>(suffix[2])]);
    big_add_scaled(blowup, term, 1);
  }
  BigPolynomial canonical_blowup;
  for (const auto& [exponent, coefficient] : blowup) {
    if (coefficient != 0) {
      canonical_blowup[exponent] = coefficient;
    }
  }
  int t_degree = 0;
  int t_factor = std::numeric_limits<int>::max();
  for (const auto& [exponent, coefficient] : canonical_blowup) {
    static_cast<void>(coefficient);
    t_degree = std::max(t_degree, exponent[0]);
    t_factor = std::min(t_factor, exponent[0]);
  }
  BigPolynomial quotient;
  for (const auto& [exponent, coefficient] : canonical_blowup) {
    Integer scale = 1;
    for (int power = exponent[0]; power < t_degree; ++power) {
      scale *= 4;
    }
    Exponent quotient_exponent = exponent;
    quotient_exponent[0] -= t_factor;
    quotient[quotient_exponent] += coefficient * scale;
  }
  BigPolynomial canonical_quotient;
  for (const auto& [exponent, coefficient] : quotient) {
    if (coefficient != 0) {
      canonical_quotient[exponent] = coefficient;
    }
  }
  return {
      maximum_w_power, grouped_w.size(), t_degree, t_factor,
      std::move(canonical_quotient)};
}

Wall121SupportFiveLargeRatioDiagnostic
wall_121_support_five_large_ratio_diagnostic() {
  constexpr int variables = 4;
  const Wall121SupportFiveLargeRatioQuotient prepared =
      wall_121_support_five_large_ratio_quotient();
  const NDBernsteinGrid quotient_grid =
      nd_bernstein_grid(prepared.quotient, variables);
  const std::size_t negative = static_cast<std::size_t>(std::count_if(
      quotient_grid.values.begin(), quotient_grid.values.end(),
      [](const Rational& value) { return value < 0; }));
  std::size_t t_zero_terms = 0U;
  std::size_t t_zero_negative = 0U;
  for (const auto& [exponent, coefficient] : prepared.quotient) {
    if (exponent[0] == 0) {
      ++t_zero_terms;
      if (coefficient < 0) {
        ++t_zero_negative;
      }
    }
  }
  return {
      prepared.maximum_w_power, prepared.grouped_blocks,
      prepared.t_degree, prepared.t_factor,
      quotient_grid.degrees, prepared.quotient.size(),
      quotient_grid.values.size(), negative,
      t_zero_terms, t_zero_negative};
}

struct Wall121SupportFiveLargeRatioFaceDiagnostic {
  std::vector<int> degrees;
  std::size_t terms = 0U;
  std::size_t coefficients = 0U;
  std::size_t negative_coefficients = 0U;
  int t_factor = -1;
  NDSubdivisionResult subdivision;
  bool lower_endpoint_positive = false;
  bool upper_endpoint_positive = false;
};

struct Wall121SupportFiveLargeRatioLeadingFace {
  int t_factor = 0;
  BigPolynomial polynomial;
  BigPolynomial truncated_blowup;
};

Wall121SupportFiveLargeRatioLeadingFace
wall_121_support_five_large_ratio_leading_face() {
  constexpr int variables = 4;
  constexpr int t_cap = 32;
  const BigPolynomial compact =
      wall_121_second_iterate_support_five_compact();
  int u_degree = 0;
  int v_degree = 0;
  int z_degree = 0;
  int maximum_w_power = 0;
  for (const auto& [exponent, coefficient] : compact) {
    static_cast<void>(coefficient);
    u_degree = std::max(u_degree, exponent[1]);
    v_degree = std::max(v_degree, exponent[2]);
    z_degree = std::max(z_degree, exponent[3]);
  }
  for (const auto& [exponent, coefficient] : compact) {
    static_cast<void>(coefficient);
    maximum_w_power = std::max(
        maximum_w_power,
        exponent[0] + 2 * (u_degree - exponent[1]));
  }
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial t = nd_variable(variables, 0);
  const BigPolynomial alpha = nd_variable(variables, 1);
  const BigPolynomial beta = nd_variable(variables, 2);
  const BigPolynomial z = nd_variable(variables, 3);
  BigPolynomial w = one;
  big_add_scaled(w, t, -1);
  BigPolynomial u_numerator = nd_power(w, 2);
  big_add_scaled(
      u_numerator,
      big_multiply(big_multiply(t, w), alpha), -2);
  big_add_scaled(
      u_numerator,
      big_multiply(nd_power(t, 2), alpha), -3);
  BigPolynomial t2_plus_u = nd_power(t, 2);
  big_add_scaled(t2_plus_u, u_numerator, 1);
  const BigPolynomial denominator = big_multiply(
      nd_power(u_numerator, 2), t2_plus_u);
  BigPolynomial two_w_plus_three_t = w;
  for (auto& [exponent, coefficient] : two_w_plus_three_t) {
    static_cast<void>(exponent);
    coefficient *= 2;
  }
  big_add_scaled(two_w_plus_three_t, t, 3);
  BigPolynomial one_minus_alpha = one;
  big_add_scaled(one_minus_alpha, alpha, -1);
  const BigPolynomial drop_numerator = big_multiply(
      big_multiply(
          big_multiply(
              big_multiply(beta, nd_power(t, 4)), w),
          two_w_plus_three_t),
      one_minus_alpha);
  BigPolynomial v_numerator = denominator;
  big_add_scaled(v_numerator, drop_numerator, -1);
  const auto truncate_t = [](const BigPolynomial& input) {
    BigPolynomial result;
    for (const auto& [exponent, coefficient] : input) {
      if (exponent[0] <= t_cap && coefficient != 0) {
        result[exponent] = coefficient;
      }
    }
    return result;
  };
  const auto multiply_truncated = [&truncate_t](
      const BigPolynomial& left, const BigPolynomial& right) {
    return truncate_t(big_multiply(left, right));
  };
  const auto power_table = [&multiply_truncated, &one](
      const BigPolynomial& base, const int maximum) {
    std::vector<BigPolynomial> powers;
    powers.reserve(static_cast<std::size_t>(maximum + 1));
    powers.push_back(one);
    for (int exponent = 1; exponent <= maximum; ++exponent) {
      powers.push_back(multiply_truncated(powers.back(), base));
    }
    return powers;
  };
  const std::vector<BigPolynomial> w_powers =
      power_table(w, maximum_w_power);
  const std::vector<BigPolynomial> u_powers =
      power_table(u_numerator, u_degree);
  const std::vector<BigPolynomial> d_powers =
      power_table(denominator, v_degree);
  const std::vector<BigPolynomial> v_powers =
      power_table(v_numerator, v_degree);
  const std::vector<BigPolynomial> z_powers =
      power_table(z, z_degree);
  std::map<Exponent, BigPolynomial> grouped_w;
  for (const auto& [exponent, coefficient] : compact) {
    const Exponent suffix{exponent[1], exponent[2], exponent[3]};
    big_add_scaled(
        grouped_w[suffix],
        w_powers[static_cast<std::size_t>(
            exponent[0] + 2 * (u_degree - exponent[1]))],
        coefficient);
  }
  BigPolynomial truncated_blowup;
  for (const auto& [suffix, w_polynomial] : grouped_w) {
    BigPolynomial term = truncate_t(w_polynomial);
    term = multiply_truncated(
        term, u_powers[static_cast<std::size_t>(suffix[0])]);
    term = multiply_truncated(
        term, d_powers[static_cast<std::size_t>(v_degree - suffix[1])]);
    term = multiply_truncated(
        term, v_powers[static_cast<std::size_t>(suffix[1])]);
    term = multiply_truncated(
        term, z_powers[static_cast<std::size_t>(suffix[2])]);
    big_add_scaled(truncated_blowup, term, 1);
  }
  int t_factor = t_cap + 1;
  for (const auto& [exponent, coefficient] : truncated_blowup) {
    if (coefficient != 0) {
      t_factor = std::min(t_factor, exponent[0]);
    }
  }
  if (t_factor > t_cap) {
    throw std::runtime_error(
        "support-five high-ratio leading face exceeds truncation cap");
  }
  BigPolynomial face;
  for (const auto& [exponent, coefficient] : truncated_blowup) {
    if (exponent[0] == t_factor && coefficient != 0) {
      const Exponent face_exponent{
          exponent[1], exponent[2], exponent[3]};
      face[face_exponent] += coefficient;
    }
  }
  BigPolynomial canonical_face;
  for (const auto& [exponent, coefficient] : face) {
    if (coefficient != 0) {
      canonical_face[exponent] = coefficient;
    }
  }
  return {
      t_factor, std::move(canonical_face), std::move(truncated_blowup)};
}

Wall121SupportFiveLargeRatioFaceDiagnostic
wall_121_support_five_large_ratio_leading_face_diagnostic(
    const int depth_limit
) {
  constexpr int variables = 3;
  const Wall121SupportFiveLargeRatioLeadingFace leading =
      wall_121_support_five_large_ratio_leading_face();
  const NDBernsteinGrid grid =
      nd_bernstein_grid(leading.polynomial, variables);
  const std::size_t negative = static_cast<std::size_t>(std::count_if(
      grid.values.begin(), grid.values.end(),
      [](const Rational& value) { return value < 0; }));
  NDSubdivisionResult subdivision;
  nd_certify_subdivision(grid, 0, depth_limit, subdivision);
  return {
      grid.degrees, leading.polynomial.size(), grid.values.size(), negative,
      leading.t_factor, std::move(subdivision)};
}

Wall121SupportFiveLargeRatioFaceDiagnostic
wall_121_support_five_large_ratio_leading_z_zero_diagnostic(
    const int depth_limit
) {
  constexpr int variables = 2;
  constexpr int t_cap = 32;
  const Wall121SupportFiveLargeRatioLeadingFace leading =
      wall_121_support_five_large_ratio_leading_face();
  int t_factor = t_cap + 1;
  for (const auto& [exponent, coefficient] : leading.truncated_blowup) {
    if (exponent[3] == 0 && coefficient != 0) {
      t_factor = std::min(t_factor, exponent[0]);
    }
  }
  if (t_factor > t_cap) {
    throw std::runtime_error(
        "support-five high-ratio z-zero face exceeds truncation cap");
  }
  BigPolynomial face;
  for (const auto& [exponent, coefficient] : leading.truncated_blowup) {
    if (exponent[0] == t_factor && exponent[3] == 0 && coefficient != 0) {
      const Exponent face_exponent{exponent[1], exponent[2]};
      face[face_exponent] += coefficient;
    }
  }
  BigPolynomial canonical_face;
  for (const auto& [exponent, coefficient] : face) {
    if (coefficient != 0) {
      canonical_face[exponent] = coefficient;
    }
  }
  const NDBernsteinGrid grid =
      nd_bernstein_grid(canonical_face, variables);
  const std::size_t negative = static_cast<std::size_t>(std::count_if(
      grid.values.begin(), grid.values.end(),
      [](const Rational& value) { return value < 0; }));
  NDSubdivisionResult subdivision;
  nd_certify_subdivision(grid, 0, depth_limit, subdivision);
  return {
      grid.degrees, canonical_face.size(), grid.values.size(), negative,
      t_factor, std::move(subdivision)};
}

Wall121SupportFiveLargeRatioFaceDiagnostic
wall_121_support_five_large_ratio_diagonal_corner_diagnostic(
    const int depth_limit
) {
  constexpr int variables = 3;
  constexpr int t_cap = 32;
  const Wall121SupportFiveLargeRatioLeadingFace leading =
      wall_121_support_five_large_ratio_leading_face();
  int t_factor = t_cap + 1;
  for (const auto& [exponent, coefficient] : leading.truncated_blowup) {
    if (coefficient != 0) {
      t_factor = std::min(t_factor, exponent[0] + exponent[3]);
    }
  }
  if (t_factor > t_cap) {
    throw std::runtime_error(
        "support-five high-ratio diagonal face exceeds truncation cap");
  }
  BigPolynomial face;
  for (const auto& [exponent, coefficient] : leading.truncated_blowup) {
    if (
        exponent[0] + exponent[3] == t_factor
        && coefficient != 0
    ) {
      const Exponent face_exponent{
          exponent[1], exponent[2], exponent[3]};
      face[face_exponent] += coefficient;
    }
  }
  BigPolynomial canonical_face;
  for (const auto& [exponent, coefficient] : face) {
    if (coefficient != 0) {
      canonical_face[exponent] = coefficient;
    }
  }
  const NDBernsteinGrid grid =
      nd_bernstein_grid(canonical_face, variables);
  const std::size_t negative = static_cast<std::size_t>(std::count_if(
      grid.values.begin(), grid.values.end(),
      [](const Rational& value) { return value < 0; }));
  NDSubdivisionResult subdivision;
  nd_certify_subdivision(grid, 0, depth_limit, subdivision);
  Wall121SupportFiveLargeRatioFaceDiagnostic result{
      grid.degrees, canonical_face.size(), grid.values.size(), negative,
      t_factor, std::move(subdivision)};
  const Exponent lower{0, 0, 0};
  const Exponent upper{0, 0, grid.degrees[2]};
  result.lower_endpoint_positive =
      grid.values[nd_grid_index(grid.degrees, lower)] > 0;
  result.upper_endpoint_positive =
      grid.values[nd_grid_index(grid.degrees, upper)] > 0;
  return result;
}

Wall121SupportFiveLargeRatioFaceDiagnostic
wall_121_support_five_large_ratio_complementary_corner_diagnostic(
    const int depth_limit
) {
  constexpr int variables = 3;
  constexpr int t_cap = 32;
  const Wall121SupportFiveLargeRatioLeadingFace leading =
      wall_121_support_five_large_ratio_leading_face();
  int total_factor = t_cap + 1;
  for (const auto& [exponent, coefficient] : leading.truncated_blowup) {
    if (coefficient != 0) {
      total_factor = std::min(
          total_factor, exponent[0] + exponent[3]);
    }
  }
  if (total_factor > t_cap) {
    throw std::runtime_error(
        "support-five high-ratio complementary face exceeds truncation cap");
  }
  BigPolynomial face;
  for (const auto& [exponent, coefficient] : leading.truncated_blowup) {
    if (
        exponent[0] + exponent[3] == total_factor
        && coefficient != 0
    ) {
      if (exponent[0] < leading.t_factor) {
        throw std::runtime_error(
            "support-five high-ratio factor lost in complementary chart");
      }
      const Exponent face_exponent{
          exponent[1], exponent[2], exponent[0] - leading.t_factor};
      face[face_exponent] += coefficient;
    }
  }
  BigPolynomial canonical_face;
  for (const auto& [exponent, coefficient] : face) {
    if (coefficient != 0) {
      canonical_face[exponent] = coefficient;
    }
  }
  const NDBernsteinGrid grid =
      nd_bernstein_grid(canonical_face, variables);
  const std::size_t negative = static_cast<std::size_t>(std::count_if(
      grid.values.begin(), grid.values.end(),
      [](const Rational& value) { return value < 0; }));
  NDSubdivisionResult subdivision;
  nd_certify_subdivision(grid, 0, depth_limit, subdivision);
  Wall121SupportFiveLargeRatioFaceDiagnostic result{
      grid.degrees, canonical_face.size(), grid.values.size(), negative,
      total_factor, std::move(subdivision)};
  const Exponent lower{0, 0, 0};
  const Exponent upper{0, 0, grid.degrees[2]};
  result.lower_endpoint_positive =
      grid.values[nd_grid_index(grid.degrees, lower)] > 0;
  result.upper_endpoint_positive =
      grid.values[nd_grid_index(grid.degrees, upper)] > 0;
  return result;
}

Wall121SupportFiveLargeRatioFaceDiagnostic
wall_121_support_five_large_ratio_t_zero_face(const int depth_limit) {
  constexpr int variables = 3;
  const Wall121SupportFiveLargeRatioQuotient prepared =
      wall_121_support_five_large_ratio_quotient();
  BigPolynomial face;
  for (const auto& [exponent, coefficient] : prepared.quotient) {
    if (exponent[0] == 0) {
      const Exponent face_exponent{
          exponent[1], exponent[2], exponent[3]};
      face[face_exponent] += coefficient;
    }
  }
  BigPolynomial canonical_face;
  for (const auto& [exponent, coefficient] : face) {
    if (coefficient != 0) {
      canonical_face[exponent] = coefficient;
    }
  }
  const NDBernsteinGrid grid = nd_bernstein_grid(canonical_face, variables);
  const std::size_t negative = static_cast<std::size_t>(std::count_if(
      grid.values.begin(), grid.values.end(),
      [](const Rational& value) { return value < 0; }));
  NDSubdivisionResult subdivision;
  nd_certify_subdivision(grid, 0, depth_limit, subdivision);
  return {
      grid.degrees, canonical_face.size(), grid.values.size(), negative,
      -1, std::move(subdivision)};
}

int replay_wall_121_rational_support_four() {
  constexpr int variables = 3;
  constexpr int support = 4;
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial r = nd_variable(variables, 0);
  const BigPolynomial u = nd_variable(variables, 1);
  const BigPolynomial v = nd_variable(variables, 2);
  const BigPolynomial zero;
  std::vector<BigPolynomial> profile(
      static_cast<std::size_t>(support + 4), zero);
  profile[1] = one;
  profile[2] = r;
  profile[3] = big_multiply(nd_power(r, 2), u);
  profile[4] = big_multiply(
      big_multiply(nd_power(r, 3), nd_power(u, 2)), v);
  const BigPolynomial certificate =
      wall_121_terminal_rational_certificate(profile, support);
  const BigPolynomial compact_raw =
      compactify_positive_variable(certificate, 0U);
  BigPolynomial compact;
  for (const auto& [exponent, coefficient] : compact_raw) {
    if (coefficient != 0) {
      compact[exponent] = coefficient;
    }
  }
  const NDBernsteinGrid grid =
      nd_bernstein_grid(compact, variables);
  const std::size_t negative = static_cast<std::size_t>(
      std::count_if(
          grid.values.begin(), grid.values.end(),
          [](const Rational& value) { return value < 0; }));
  auto [lower_half, upper_half] = nd_split_grid(grid, 0);
  auto [middle_quarter, top_quarter] =
      nd_split_grid(upper_half, 0);
  static_cast<void>(top_quarter);
  NDSubdivisionResult moderate;
  nd_certify_subdivision(lower_half, 0, 30, moderate);
  nd_certify_subdivision(middle_quarter, 0, 30, moderate);
  const Wall121LargeRatioBlowupResult large =
      wall_121_large_ratio_blowup(compact, grid);
  if (
      grid.degrees != std::vector<int>{47, 25, 8}
      || grid.values.size() != 11232U
      || negative != 1901U
      || moderate.nodes != 50U
      || moderate.leaves != 26U
      || moderate.unresolved != 0U
      || moderate.maximum_depth != 12
      || large.grouped_blocks != 153U
      || large.t_degree != 145
      || large.t_factor != 8
      || large.degrees != std::vector<int>{137, 41, 8}
      || large.terms != 40204U
      || large.coefficients != 52164U
      || large.negative_coefficients != 97U
      || large.t_zero_terms != 1U
      || large.t_zero_negative != 0U
      || large.distinct_denominators != 24575U
      || large.common_denominator_bits != 244U
      || large.subdivision.counts.nodes != 5U
      || large.subdivision.counts.leaves != 3U
      || large.subdivision.corner_leaves != 0U
      || large.subdivision.counts.unresolved != 0U
      || large.subdivision.counts.maximum_depth != 2
  ) {
    throw std::runtime_error(
        "wall (1,2,1) support-four rational certificate mismatch");
  }
  std::cout
      << "SU2_WALL_121_RATIONAL_SUPPORT_FOUR"
      << " support=4"
      << " moderate_degrees=(47,25,8)"
      << " moderate_coefficients=11232"
      << " moderate_initial_negative=1901"
      << " moderate_nodes=50"
      << " moderate_leaves=26"
      << " moderate_unresolved=0"
      << " moderate_depth=12"
      << " blowup_t_factor=8"
      << " blowup_degrees=(137,41,8)"
      << " blowup_terms=40204"
      << " blowup_coefficients=52164"
      << " blowup_initial_negative=97"
      << " blowup_t_zero_terms=1"
      << " blowup_t_zero_negative=0"
      << " blowup_integer_denominators=24575"
      << " blowup_integer_common_denominator_bits=244"
      << " blowup_nodes=5"
      << " blowup_leaves=3"
      << " blowup_corner_leaves=0"
      << " blowup_unresolved=0"
      << " blowup_depth=2"
      << " result=PASS_EXACT\n";
  return EXIT_SUCCESS;
}

int replay_wall_121_terminal_append_law() {
  constexpr int variables = 6;
  const BigPolynomial z = nd_variable(variables, 0);
  const BigPolynomial a = nd_variable(variables, 1);
  const BigPolynomial b = nd_variable(variables, 2);
  const BigPolynomial c = nd_variable(variables, 3);
  const BigPolynomial d = nd_variable(variables, 4);
  const BigPolynomial x = nd_variable(variables, 5);
  const BigPolynomial zero;
  const auto terminal_current = [&zero](
      const std::vector<BigPolynomial>& profile, const int support) {
    std::vector<BigPolynomial> g(
        static_cast<std::size_t>(support + 4), zero);
    std::vector<BigPolynomial> h(
        static_cast<std::size_t>(support + 4), zero);
    std::vector<BigPolynomial> k(
        static_cast<std::size_t>(support + 4), zero);
    for (int index = 1; index <= support + 1; ++index) {
      g[static_cast<std::size_t>(index)] = big_multiply(
          profile[static_cast<std::size_t>(index)],
          profile[static_cast<std::size_t>(index)]);
      big_add_scaled(
          g[static_cast<std::size_t>(index)],
          big_multiply(
              profile[static_cast<std::size_t>(index - 1)],
              profile[static_cast<std::size_t>(index + 1)]), -1);
    }
    for (int index = 1; index <= support + 2; ++index) {
      h[static_cast<std::size_t>(index)] = big_multiply(
          profile[static_cast<std::size_t>(index)],
          profile[static_cast<std::size_t>(index - 1)]);
      if (index >= 2) {
        big_add_scaled(
            h[static_cast<std::size_t>(index)],
            big_multiply(
                profile[static_cast<std::size_t>(index - 2)],
                profile[static_cast<std::size_t>(index + 1)]), -1);
      }
    }
    for (int index = 1; index <= support + 1; ++index) {
      k[static_cast<std::size_t>(index)] =
          g[static_cast<std::size_t>(index)];
      big_add_scaled(
          k[static_cast<std::size_t>(index)],
          h[static_cast<std::size_t>(index)], 1);
      big_add_scaled(
          k[static_cast<std::size_t>(index)],
          h[static_cast<std::size_t>(index + 1)], 1);
    }
    BigPolynomial result = big_multiply(
        nd_power(profile[1], 3), profile[2]);
    for (int index = 1; index <= support; ++index) {
      BigPolynomial g_difference =
          g[static_cast<std::size_t>(index)];
      big_add_scaled(
          g_difference,
          g[static_cast<std::size_t>(index + 1)], -1);
      BigPolynomial k_difference =
          k[static_cast<std::size_t>(index)];
      big_add_scaled(
          k_difference,
          k[static_cast<std::size_t>(index + 1)], -1);
      big_add_scaled(
          result,
          big_multiply(g_difference, k_difference), 1);
    }
    return result;
  };
  std::vector<BigPolynomial> appended(10, zero);
  appended[1] = z;
  appended[2] = a;
  appended[3] = b;
  appended[4] = c;
  appended[5] = d;
  appended[6] = x;
  std::vector<BigPolynomial> truncated = appended;
  truncated[6] = zero;
  BigPolynomial difference = terminal_current(appended, 6);
  big_add_scaled(
      difference, terminal_current(truncated, 5), -1);
  BigPolynomial linear = nd_power(d, 3);
  big_add_scaled(
      linear, big_multiply(big_multiply(a, c), d), -1);
  big_add_scaled(
      linear, big_multiply(c, nd_power(d, 2)), -4);
  big_add_scaled(
      linear, big_multiply(nd_power(c, 2), d), -2);
  big_add_scaled(linear, nd_power(c, 3), 2);
  big_add_scaled(
      linear, big_multiply(big_multiply(a, b), c), -1);
  big_add_scaled(
      linear, big_multiply(big_multiply(b, c), d), -2);
  big_add_scaled(
      linear, big_multiply(nd_power(b, 2), d), 1);
  big_add_scaled(linear, nd_power(b, 3), 1);
  BigPolynomial quadratic = nd_power(c, 2);
  for (auto& [exponent, coefficient] : quadratic) {
    static_cast<void>(exponent);
    coefficient *= 2;
  }
  big_add_scaled(quadratic, nd_power(d, 2), -2);
  big_add_scaled(quadratic, big_multiply(c, d), -2);
  big_add_scaled(quadratic, big_multiply(b, c), 1);
  BigPolynomial cubic = b;
  big_add_scaled(cubic, c, 2);
  big_add_scaled(cubic, d, 1);
  BigPolynomial expected = big_multiply(x, linear);
  big_add_scaled(
      expected, big_multiply(nd_power(x, 2), quadratic), 1);
  big_add_scaled(
      expected, big_multiply(nd_power(x, 3), cubic), 1);
  big_add_scaled(expected, nd_power(x, 4), 2);
  const auto canonicalize = [](const BigPolynomial& polynomial) {
    BigPolynomial result;
    for (const auto& [exponent, coefficient] : polynomial) {
      if (coefficient != 0) {
        result[exponent] = coefficient;
      }
    }
    return result;
  };
  if (canonicalize(difference) != canonicalize(expected)) {
    throw std::runtime_error(
        "wall (1,2,1) terminal append identity mismatch");
  }
  const Rational witness = evaluate_polynomial(
      rational_polynomial(difference),
      {Rational(7, 5), Rational(1), Rational(4, 3),
       Rational(5, 3), Rational(5, 3), Rational(4, 3)});
  if (witness != Rational(-1016, 81)) {
    throw std::runtime_error(
        "wall (1,2,1) terminal append obstruction mismatch");
  }
  std::cout
      << "SU2_WALL_121_TERMINAL_APPEND_LAW"
      << " arbitrary_depth_previous=(a,b,c,d)"
      << " append=x"
      << " degree=4"
      << " delta=x*(L+Q*x+(b+2*c+d)*x^2+2*x^3)"
      << " L=d^3-a*c*d-4*c*d^2-2*c^2*d+2*c^3-a*b*c-2*b*c*d+b^2*d+b^3"
      << " Q=-2*d^2-2*c*d+2*c^2+b*c"
      << " normalized_witness=(a=1,b=4/3,c=5/3,d=5/3,x=4/3)"
      << " witness_delta=" << witness
      << " result=PASS_EXACT\n";
  return EXIT_SUCCESS;
}

int replay_wall_121_small_ratio_payment() {
  const Rational r(2, 5);
  const Rational tail_prefactor =
      2 * (1 + r + r * r)
      / (1 - r * r * r * r);
  const Rational current_numerator_at_endpoint(18129, 390625);
  const Rational current_margin =
      r * current_numerator_at_endpoint / (1 - r * r * r * r);
  const Rational c2_max(173, 125);
  const Rational c3_max(39, 25);
  const Rational critical_half_min(7, 4);
  const Rational b_max(25214, 15625);
  const Rational demand_max =
      (c2_max + c3_max) / 4;
  if (
      tail_prefactor != Rational(650, 203)
      || !(current_numerator_at_endpoint > 0)
      || !(current_margin > 0)
      || c2_max != Rational(173, 125)
      || c3_max != Rational(39, 25)
      || !(critical_half_min > b_max)
      || demand_max != Rational(92, 125)
      || !(demand_max < 1)
  ) {
    throw std::runtime_error(
        "wall (1,2,1) small-ratio payment mismatch");
  }
  std::cout
      << "SU2_WALL_121_SMALL_RATIO_PAYMENT"
      << " ratio_max=2/5"
      << " tail_prefactor=650/203"
      << " current_numerator_at_ratio_max=18129/390625>0"
      << " current_margin>0"
      << " C2_max=173/125"
      << " C3_max=39/25"
      << " B_max=25214/15625"
      << " critical_half_min=7/4"
      << " demand_max=92/125<1"
      << " result=PASS_EXACT\n";
  return EXIT_SUCCESS;
}

int replay_wall_121_renewal_kernel() {
  constexpr int variables = 4;
  const BigPolynomial one = nd_constant(variables, 1);
  const BigPolynomial zero;
  const BigPolynomial r = nd_variable(variables, 0);
  const BigPolynomial s = nd_variable(variables, 1);
  const BigPolynomial t = nd_variable(variables, 2);
  const BigPolynomial q = nd_variable(variables, 3);
  const BigPolynomial rs = big_multiply(r, s);
  const BigPolynomial rst = big_multiply(rs, t);
  std::vector<BigPolynomial> profile(8U, zero);
  profile[1] = one;
  profile[2] = r;
  profile[3] = rs;
  profile[4] = rst;
  const BigPolynomial terminal = wall_121_terminal_current(profile, 4);
  std::vector<BigPolynomial> next_profile(8U, zero);
  next_profile[1] = one;
  next_profile[2] = s;
  next_profile[3] = big_multiply(s, t);
  next_profile[4] = big_multiply(next_profile[3], q);
  const BigPolynomial next_terminal =
      wall_121_terminal_current(next_profile, 4);

  const BigPolynomial b = r;
  const BigPolynomial c = rs;
  const BigPolynomial d = rst;
  const BigPolynomial x = big_multiply(rst, q);
  BigPolynomial linear = nd_power(d, 3);
  big_add_scaled(linear, big_multiply(c, d), -1);
  big_add_scaled(linear, big_multiply(c, nd_power(d, 2)), -4);
  big_add_scaled(linear, big_multiply(nd_power(c, 2), d), -2);
  big_add_scaled(linear, nd_power(c, 3), 2);
  big_add_scaled(linear, big_multiply(b, c), -1);
  big_add_scaled(linear, big_multiply(big_multiply(b, c), d), -2);
  big_add_scaled(linear, big_multiply(nd_power(b, 2), d), 1);
  big_add_scaled(linear, nd_power(b, 3), 1);
  BigPolynomial quadratic = nd_power(d, 2);
  for (auto& [exponent, coefficient] : quadratic) {
    static_cast<void>(exponent);
    coefficient *= -2;
  }
  big_add_scaled(quadratic, big_multiply(c, d), -2);
  big_add_scaled(quadratic, nd_power(c, 2), 2);
  big_add_scaled(quadratic, big_multiply(b, c), 1);
  BigPolynomial cubic = b;
  big_add_scaled(cubic, c, 2);
  big_add_scaled(cubic, d, 1);
  BigPolynomial delta = big_multiply(x, linear);
  big_add_scaled(delta, big_multiply(nd_power(x, 2), quadratic), 1);
  big_add_scaled(delta, big_multiply(nd_power(x, 3), cubic), 1);
  big_add_scaled(delta, nd_power(x, 4), 2);

  BigPolynomial renewal = terminal;
  big_add_scaled(renewal, delta, 1);
  big_add_scaled(
      renewal, big_multiply(nd_power(r, 4), next_terminal), -1);

  const auto add = [](BigPolynomial& target, const BigPolynomial& term,
                      const long coefficient) {
    big_add_scaled(target, term, coefficient);
  };
  BigPolynomial expected = one;
  add(expected, r, 1);
  add(expected, rs, 2);
  add(expected, rst, 1);
  add(expected, nd_power(r, 2), -2);
  add(expected, big_multiply(nd_power(r, 2), s), -2);
  add(expected, big_multiply(nd_power(rs, 2), t), 1);
  add(expected, nd_power(rs, 2), 2);
  add(expected, nd_power(r, 3), 1);
  add(expected, big_multiply(nd_power(r, 3), s), -4);
  add(expected, big_multiply(nd_power(r, 3), nd_power(s, 2)), -2);
  add(expected, big_multiply(big_multiply(nd_power(r, 3), nd_power(s, 2)), t), -2);
  add(expected, big_multiply(big_multiply(nd_power(r, 3), nd_power(s, 2)),
                             big_multiply(t, q)), -1);
  add(expected, big_multiply(big_multiply(nd_power(r, 3), nd_power(s, 2)),
                             nd_power(t, 2)), 1);
  add(expected, big_multiply(nd_power(r, 3), nd_power(s, 3)), 2);
  add(expected, big_multiply(big_multiply(nd_power(r, 3), nd_power(s, 3)),
                             nd_power(t, 3)), 1);
  add(expected, big_multiply(
                    big_multiply(nd_power(r, 3), nd_power(s, 3)),
                    big_multiply(nd_power(t, 2), q)),
      -1);
  add(expected, nd_power(r, 4), 1);
  const auto canonicalize = [](const BigPolynomial& polynomial) {
    BigPolynomial result;
    for (const auto& [exponent, coefficient] : polynomial) {
      if (coefficient != 0) {
        result.emplace(exponent, coefficient);
      }
    }
    return result;
  };
  if (
      canonicalize(renewal) != canonicalize(expected)
      || canonicalize(expected).size() != 18U
  ) {
    throw std::runtime_error(
        "wall (1,2,1) renewal-kernel identity mismatch");
  }
  constexpr int geometric_variables = 1;
  const BigPolynomial geometric_one = nd_constant(geometric_variables, 1);
  const BigPolynomial rho = nd_variable(geometric_variables, 0);
  std::vector<BigPolynomial> geometric_profile(9U);
  geometric_profile[1] = geometric_one;
  geometric_profile[2] = rho;
  geometric_profile[3] = nd_power(rho, 2);
  geometric_profile[4] = nd_power(rho, 3);
  const BigPolynomial geometric_terminal =
      wall_121_terminal_current(geometric_profile, 4);
  geometric_profile[5] = nd_power(rho, 4);
  const BigPolynomial geometric_extended =
      wall_121_terminal_current(geometric_profile, 5);
  BigPolynomial geometric_expected = geometric_one;
  big_add_scaled(geometric_expected, rho, 1);
  big_add_scaled(geometric_expected, nd_power(rho, 11), 1);
  big_add_scaled(geometric_expected, nd_power(rho, 12), 2);
  BigPolynomial geometric_extended_expected = geometric_one;
  big_add_scaled(geometric_extended_expected, rho, 1);
  big_add_scaled(geometric_extended_expected, nd_power(rho, 15), 1);
  big_add_scaled(geometric_extended_expected, nd_power(rho, 16), 2);
  BigPolynomial geometric_delta = geometric_extended;
  big_add_scaled(geometric_delta, geometric_terminal, -1);
  BigPolynomial rho_four_minus_one = nd_power(rho, 4);
  big_add_scaled(rho_four_minus_one, geometric_one, -1);
  BigPolynomial one_plus_two_rho = geometric_one;
  big_add_scaled(one_plus_two_rho, rho, 2);
  const BigPolynomial geometric_delta_expected = big_multiply(
      big_multiply(nd_power(rho, 11), rho_four_minus_one),
      one_plus_two_rho);
  BigPolynomial geometric_renewal = geometric_terminal;
  big_add_scaled(geometric_renewal, geometric_delta, 1);
  big_add_scaled(
      geometric_renewal,
      big_multiply(nd_power(rho, 4), geometric_terminal), -1);
  BigPolynomial geometric_renewal_expected = geometric_one;
  big_add_scaled(geometric_renewal_expected, rho, 1);
  BigPolynomial one_minus_rho_four = geometric_one;
  big_add_scaled(one_minus_rho_four, nd_power(rho, 4), -1);
  geometric_renewal_expected = big_multiply(
      geometric_renewal_expected, one_minus_rho_four);
  if (
      canonicalize(geometric_terminal) != canonicalize(geometric_expected)
      || canonicalize(geometric_extended)
             != canonicalize(geometric_extended_expected)
      || canonicalize(geometric_delta)
             != canonicalize(geometric_delta_expected)
      || canonicalize(geometric_renewal)
             != canonicalize(geometric_renewal_expected)
  ) {
    throw std::runtime_error(
        "wall (1,2,1) geometric renewal calibration mismatch");
  }
  for (int support = 4; support <= 12; ++support) {
    std::vector<BigPolynomial> geometric_tail(
        static_cast<std::size_t>(support + 4));
    for (int index = 1; index <= support; ++index) {
      geometric_tail[static_cast<std::size_t>(index)] =
          nd_power(rho, index - 1);
    }
    const BigPolynomial actual =
        wall_121_terminal_current(geometric_tail, support);
    BigPolynomial expected_tail = geometric_one;
    big_add_scaled(expected_tail, rho, 1);
    big_add_scaled(expected_tail, nd_power(rho, 4 * support - 5), 1);
    big_add_scaled(expected_tail, nd_power(rho, 4 * support - 4), 2);
    if (canonicalize(actual) != canonicalize(expected_tail)) {
      throw std::runtime_error(
          "wall (1,2,1) geometric terminal-tail identity mismatch");
    }
  }
  std::vector<BigPolynomial> tangent_profile(8U);
  tangent_profile[1] = nd_constant(geometric_variables, 8);
  tangent_profile[2] = geometric_one;
  const BigPolynomial tangent_terminal =
      wall_121_terminal_current(tangent_profile, 4);
  const BigPolynomial tangent_terminal_expected =
      nd_constant(geometric_variables, 4490);
  const BigPolynomial tangent_barrier_expected =
      nd_constant(geometric_variables, 4608);
  if (
      canonicalize(tangent_terminal)
          != canonicalize(tangent_terminal_expected)
      || !(4490 < 4608)
  ) {
    throw std::runtime_error(
        "wall (1,2,1) tangent-barrier obstruction mismatch");
  }
  std::vector<BigPolynomial> one_ratio_profile(8U);
  one_ratio_profile[1] = nd_constant(geometric_variables, 144);
  one_ratio_profile[2] = nd_constant(geometric_variables, 84);
  one_ratio_profile[3] = nd_constant(geometric_variables, 7);
  const BigPolynomial one_ratio_terminal =
      wall_121_terminal_current(one_ratio_profile, 3);
  const BigPolynomial one_ratio_terminal_expected =
      nd_constant(geometric_variables, 566506574);
  const BigPolynomial one_ratio_barrier_expected =
      nd_constant(geometric_variables, 573101568);
  if (
      canonicalize(one_ratio_terminal)
          != canonicalize(one_ratio_terminal_expected)
      || !(566506574 < 573101568)
  ) {
    throw std::runtime_error(
        "wall (1,2,1) one-ratio barrier obstruction mismatch");
  }
  std::vector<BigPolynomial> two_ratio_prefix(8U);
  two_ratio_prefix[1] = nd_constant(geometric_variables, 2);
  two_ratio_prefix[2] = nd_constant(geometric_variables, 2);
  two_ratio_prefix[3] = nd_constant(geometric_variables, 2);
  const BigPolynomial two_ratio_terminal =
      wall_121_terminal_current(two_ratio_prefix, 3);
  two_ratio_prefix[4] = geometric_one;
  BigPolynomial two_ratio_append =
      wall_121_terminal_current(two_ratio_prefix, 4);
  big_add_scaled(two_ratio_append, two_ratio_terminal, -1);
  if (
      canonicalize(two_ratio_append)
          != canonicalize(nd_constant(geometric_variables, -18))
  ) {
    throw std::runtime_error(
        "wall (1,2,1) two-ratio barrier obstruction mismatch");
  }
  BigPolynomial one_step = terminal;
  big_add_scaled(one_step, delta, 1);
  const RationalPolynomial one_step_rational =
      rational_polynomial(one_step);
  const std::vector<Rational> endpoint_control_state{
      Rational(1, 2), Rational(1, 2), Rational(3, 8)};
  const auto one_step_value = [&](const Rational& control) {
    std::vector<Rational> point = endpoint_control_state;
    point.push_back(control);
    return evaluate_polynomial(one_step_rational, point);
  };
  const Rational endpoint_zero = one_step_value(Rational(0));
  const Rational endpoint_upper = one_step_value(Rational(3, 8));
  const Rational interior_control = one_step_value(Rational(1, 4));
  std::vector<Rational> terminal_point = endpoint_control_state;
  terminal_point.push_back(Rational(0));
  const Rational terminal_control_value = evaluate_polynomial(
      rational_polynomial(terminal), terminal_point);
  const Rational endpoint_ratio = endpoint_control_state[0];
  const Rational endpoint_second_ratio = endpoint_control_state[1];
  const Rational endpoint_third_ratio = endpoint_control_state[2];
  const Rational endpoint_second_coordinate = endpoint_ratio;
  const Rational endpoint_third_coordinate =
      endpoint_ratio * endpoint_second_ratio;
  const Rational endpoint_fourth_coordinate =
      endpoint_third_coordinate * endpoint_third_ratio;
  const Rational endpoint_a(13, 8);
  const Rational endpoint_b(507, 256);
  const Rational endpoint_c(7, 4);
  const Rational endpoint_a_from_state =
      1 + endpoint_ratio * endpoint_ratio
      + (1 + endpoint_ratio) * endpoint_third_coordinate;
  const Rational endpoint_b_from_state =
      3 * endpoint_ratio + 2 * endpoint_ratio * endpoint_ratio
      + endpoint_ratio * endpoint_third_coordinate
      + endpoint_ratio * endpoint_fourth_coordinate
      + endpoint_ratio * endpoint_third_coordinate
            * endpoint_fourth_coordinate
      - endpoint_ratio * endpoint_ratio * endpoint_ratio
      - endpoint_third_coordinate * endpoint_third_coordinate
      - endpoint_third_coordinate * endpoint_third_coordinate
            * endpoint_third_coordinate;
  const Rational endpoint_c_from_state =
      1 + endpoint_ratio + endpoint_third_coordinate;
  const Rational endpoint_derivative =
      -endpoint_b + 2 * endpoint_a / endpoint_ratio
      + 3 * endpoint_c / (endpoint_ratio * endpoint_ratio);
  if (
      endpoint_zero != Rational(766557, 524288)
      || interior_control != Rational(196216899, 134217728)
      || endpoint_upper != Rational(3139522317, 2147483648)
      || terminal_control_value != endpoint_zero
      || endpoint_a_from_state != endpoint_a
      || endpoint_b_from_state != endpoint_b
      || endpoint_c_from_state != endpoint_c
      || endpoint_b != Rational(507, 256)
      || endpoint_derivative != Rational(6533, 256)
      || !(interior_control < endpoint_zero)
      || !(interior_control < endpoint_upper)
  ) {
    throw std::runtime_error(
        "wall (1,2,1) endpoint-control obstruction mismatch");
  }
  std::cout
      << "SU2_WALL_121_RENEWAL_KERNEL"
      << " variables=(r,s,t,q)"
      << " terms=18"
      << " q_degree=1"
      << " q_coefficient=-r^3*s^2*t*(1+s*t)"
      << " geometric_T4=1+rho+rho^11*(1+2*rho)"
      << " geometric_delta=rho^11*(rho^4-1)*(1+2*rho)"
      << " geometric_R=(1+rho)*(1-rho^4)"
      << " geometric_terminal_supports=4..12"
      << " tangent_barrier_boundary_residual=-59/2048"
      << " one_ratio_barrier_residual=-3297497/214990848"
      << " two_ratio_barrier_append=-9/8"
      << " endpoint_control_values=(q0=" << endpoint_zero
      << ",q1/4=" << interior_control
      << ",q3/8=" << endpoint_upper
      << ",terminate=" << terminal_control_value << ')'
      << " endpoint_control_critical_data=(B=" << endpoint_b
      << ",Phi_prime=" << endpoint_derivative << ')'
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
  std::vector<std::uint64_t> first_unresolved_cell;
  std::vector<int> first_unresolved_splits;
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
        first_unresolved_cell = result.first_unresolved_cell;
        first_unresolved_splits = result.first_unresolved_splits;
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
            << '\n' << std::flush;
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
        << " first_unresolved_depth=" << first_unresolved_depth
        << " first_unresolved_cell=(";
    for (std::size_t index = 0U;
         index < first_unresolved_cell.size();
         ++index) {
      if (index != 0U) {
        std::cout << ',';
      }
      std::cout << first_unresolved_cell[index];
    }
    std::cout << ") first_unresolved_splits=(";
    for (std::size_t index = 0U;
         index < first_unresolved_splits.size();
         ++index) {
      if (index != 0U) {
        std::cout << ',';
      }
      std::cout << first_unresolved_splits[index];
    }
    std::cout << ')';
  }
  std::cout
      << " result=" << (unresolved == 0U ? "PASS" : "INCOMPLETE")
      << '\n' << std::flush;
  return EXIT_SUCCESS;
}

int analyze_ratio_cube_best_subdivision(
    const int support,
    const int depth_limit,
    const int antidiagonal,
    const int depth
) {
  if (2 * depth >= antidiagonal) {
    throw std::invalid_argument("radial depth must be strictly below half");
  }
  const Polynomial target = substitute_differences(
      radial_polynomial(support, antidiagonal, depth), support);
  const BigPolynomial parameterized =
      ratio_cube_parameterization(target, support);
  const NDBernsteinGrid grid = nd_bernstein_grid(parameterized, support);
  NDSubdivisionResult result;
  nd_certify_subdivision_best(grid, depth_limit, result);
  std::cout
      << "SU2_AUTOCORRELATION_RATIO_CUBE_BEST"
      << " support=" << support
      << " antidiagonal=" << antidiagonal
      << " depth=" << depth
      << " coefficients=" << grid.values.size()
      << " initial_negative=" << nd_negative_count(grid)
      << " depth_limit=" << depth_limit
      << " nodes=" << result.nodes
      << " leaves=" << result.leaves
      << " unresolved=" << result.unresolved
      << " maximum_depth=" << result.maximum_depth;
  if (result.has_first_unresolved) {
    std::cout << " first_unresolved_cell=(";
    for (std::size_t index = 0U;
         index < result.first_unresolved_cell.size();
         ++index) {
      if (index != 0U) {
        std::cout << ',';
      }
      std::cout << result.first_unresolved_cell[index];
    }
    std::cout << ") first_unresolved_splits=(";
    for (std::size_t index = 0U;
         index < result.first_unresolved_splits.size();
         ++index) {
      if (index != 0U) {
        std::cout << ',';
      }
      std::cout << result.first_unresolved_splits[index];
    }
    std::cout << ')';
  }
  std::cout
      << " result=" << (result.unresolved == 0U ? "PASS" : "INCOMPLETE")
      << '\n' << std::flush;
  return EXIT_SUCCESS;
}

int analyze_ratio_cube_local_corner(
    const int support,
    const int antidiagonal,
    const int depth,
    const std::vector<int>& cell,
    const std::vector<int>& splits
) {
  if (2 * depth >= antidiagonal) {
    throw std::invalid_argument("radial depth must be strictly below half");
  }
  if (static_cast<int>(cell.size()) != support
      || splits.size() != cell.size()) {
    throw std::invalid_argument("cell dimension does not match support");
  }
  const Polynomial target = substitute_differences(
      radial_polynomial(support, antidiagonal, depth), support);
  BigPolynomial local = restrict_to_dyadic_cell(
      ratio_cube_parameterization(target, support), cell, splits);
  for (std::size_t variable = 0U; variable < cell.size(); ++variable) {
    local = reflect_unit_variable(local, variable);
  }
  local = remove_common_monomial_factor(local);
  std::size_t negative = 0U;
  bool has_first_negative = false;
  Exponent first_negative;
  Integer first_negative_value = 0;
  Integer negative_sum = 0;
  Integer constant = 0;
  int leading_degree = std::numeric_limits<int>::max();
  for (const auto& [exponent, coefficient] : local) {
    if (coefficient != 0) {
      int total_degree = 0;
      for (const int power : exponent) {
        total_degree += power;
      }
      leading_degree = std::min(leading_degree, total_degree);
    }
    const bool is_constant = std::all_of(
        exponent.begin(), exponent.end(),
        [](const int value) { return value == 0; });
    if (is_constant) {
      constant = coefficient;
    }
    if (coefficient < 0) {
      ++negative;
      negative_sum += coefficient;
      if (!has_first_negative) {
        has_first_negative = true;
        first_negative = exponent;
        first_negative_value = coefficient;
      }
    }
  }
  std::cout
      << "SU2_AUTOCORRELATION_RATIO_CUBE_LOCAL_CORNER"
      << " support=" << support
      << " antidiagonal=" << antidiagonal
      << " depth=" << depth
      << " terms=" << local.size()
      << " negative=" << negative
      << " constant=" << constant
      << " leading_degree=" << leading_degree
      << " coefficient_lower_bound=" << constant + negative_sum;
  if (has_first_negative) {
    std::cout << " first_negative_exponent=(";
    for (std::size_t index = 0U; index < first_negative.size(); ++index) {
      if (index != 0U) {
        std::cout << ',';
      }
      std::cout << first_negative[index];
    }
    std::cout << ") first_negative_value=" << first_negative_value;
  }
  std::cout
      << " result="
      << (negative == 0U || constant + negative_sum >= 0
              ? "PASS"
              : "INCOMPLETE")
      << '\n' << std::flush;
  if (leading_degree != std::numeric_limits<int>::max()) {
    std::cout << "SU2_AUTOCORRELATION_RATIO_CUBE_LOCAL_LEADING terms={";
    bool first = true;
    for (const auto& [exponent, coefficient] : local) {
      int total_degree = 0;
      for (const int power : exponent) {
        total_degree += power;
      }
      if (total_degree != leading_degree) {
        continue;
      }
      if (!first) {
        std::cout << ',';
      }
      first = false;
      std::cout << '(';
      for (std::size_t index = 0U; index < exponent.size(); ++index) {
        if (index != 0U) {
          std::cout << ',';
        }
        std::cout << exponent[index];
      }
      std::cout << "):" << coefficient;
    }
    std::cout << "}\n" << std::flush;
  }
  BigPolynomial outer_boundary;
  BigPolynomial outer_quotient;
  const std::size_t outer = cell.size() - 1U;
  for (const auto& [exponent, coefficient] : local) {
    Exponent reduced = exponent;
    if (reduced[outer] == 0) {
      outer_boundary[std::move(reduced)] += coefficient;
    } else {
      --reduced[outer];
      outer_quotient[std::move(reduced)] += coefficient;
    }
  }
  const auto negative_terms = [](const BigPolynomial& polynomial) {
    return static_cast<std::size_t>(std::count_if(
        polynomial.begin(), polynomial.end(),
        [](const auto& term) { return term.second < 0; }));
  };
  const std::size_t boundary_negative = negative_terms(outer_boundary);
  const std::size_t quotient_negative = negative_terms(outer_quotient);
  std::cout
      << "SU2_AUTOCORRELATION_RATIO_CUBE_LOCAL_OUTER_SPLIT"
      << " boundary_terms=" << outer_boundary.size()
      << " boundary_negative=" << boundary_negative
      << " quotient_terms=" << outer_quotient.size()
      << " quotient_negative=" << quotient_negative
      << " result="
      << (boundary_negative == 0U && quotient_negative == 0U
              ? "PASS"
              : "INCOMPLETE")
      << '\n' << std::flush;
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
               == "--replay-ordinary-12-complete-current-sign"
    ) {
      return replay_ordinary_12_complete_current_sign();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-ordinary-one-row-universal-obstruction"
    ) {
      return replay_ordinary_one_row_universal_obstruction();
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
               == "--replay-wall-121-geometric-demand-ceiling"
    ) {
      return replay_wall_121_geometric_demand_ceiling();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-wall-121-cubic-corrected-demand-ceiling"
    ) {
      return replay_wall_121_cubic_corrected_demand_ceiling();
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
               == "--replay-wall-121-rational-support-three"
    ) {
      return replay_wall_121_rational_support_three();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-wall-121-rational-support-four"
    ) {
      return replay_wall_121_rational_support_four();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--replay-wall-121-terminal-append-law"
    ) {
      return replay_wall_121_terminal_append_law();
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
               == "--replay-wall-121-renewal-kernel"
    ) {
      return replay_wall_121_renewal_kernel();
    }
    if (
        argc == 2
        && std::string{argv[1]}
               == "--support-five-rat2-high-ratio-diagnostic"
    ) {
      const Wall121SupportFiveLargeRatioDiagnostic result =
          wall_121_support_five_large_ratio_diagnostic();
      std::cout << "SU2_WALL_121_RAT2_SUPPORT_FIVE_HIGH_RATIO"
                << " maximum_w_power=" << result.maximum_w_power
                << " grouped_blocks=" << result.grouped_blocks
                << " t_degree=" << result.t_degree
                << " t_factor=" << result.t_factor
                << " degrees=(" << result.degrees[0] << ','
                << result.degrees[1] << ',' << result.degrees[2] << ','
                << result.degrees[3] << ')'
                << " terms=" << result.terms
                << " bernstein_coefficients=" << result.coefficients
                << " initial_negative=" << result.negative_coefficients
                << " t_zero_terms=" << result.t_zero_terms
                << " t_zero_negative=" << result.t_zero_negative
                << '\n';
      return EXIT_SUCCESS;
    }
    if (
        argc == 3
        && std::string{argv[1]}
               == "--support-five-rat2-high-ratio-t-zero"
    ) {
      const int depth_limit = parse_positive(argv[2], "depth limit");
      if (depth_limit > 24) {
        throw std::invalid_argument("depth limit must be at most 24");
      }
      const Wall121SupportFiveLargeRatioFaceDiagnostic result =
          wall_121_support_five_large_ratio_t_zero_face(depth_limit);
      std::cout << "SU2_WALL_121_RAT2_SUPPORT_FIVE_HIGH_RATIO_T_ZERO"
                << " degrees=(" << result.degrees[0] << ','
                << result.degrees[1] << ',' << result.degrees[2] << ')'
                << " terms=" << result.terms
                << " bernstein_coefficients=" << result.coefficients
                << " initial_negative=" << result.negative_coefficients
                << " depth_limit=" << depth_limit
                << " nodes=" << result.subdivision.nodes
                << " leaves=" << result.subdivision.leaves
                << " unresolved=" << result.subdivision.unresolved
                << " maximum_depth=" << result.subdivision.maximum_depth;
      if (result.subdivision.has_first_unresolved) {
        std::cout << " first_unresolved_cell=(";
        for (std::size_t index = 0U;
             index < result.subdivision.first_unresolved_cell.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.subdivision.first_unresolved_cell[index];
        }
        std::cout << ") first_unresolved_splits=(";
        for (std::size_t index = 0U;
             index < result.subdivision.first_unresolved_splits.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.subdivision.first_unresolved_splits[index];
        }
        std::cout << ')';
      }
      std::cout << '\n';
      return EXIT_SUCCESS;
    }
    if (
        argc == 3
        && std::string{argv[1]}
               == "--support-five-rat2-high-ratio-leading-face"
    ) {
      const int depth_limit = parse_positive(argv[2], "depth limit");
      if (depth_limit > 24) {
        throw std::invalid_argument("depth limit must be at most 24");
      }
      const Wall121SupportFiveLargeRatioFaceDiagnostic result =
          wall_121_support_five_large_ratio_leading_face_diagnostic(
              depth_limit);
      std::cout << "SU2_WALL_121_RAT2_SUPPORT_FIVE_HIGH_RATIO_LEADING_FACE"
                << " t_factor=" << result.t_factor
                << " degrees=(" << result.degrees[0] << ','
                << result.degrees[1] << ',' << result.degrees[2] << ')'
                << " terms=" << result.terms
                << " bernstein_coefficients=" << result.coefficients
                << " initial_negative=" << result.negative_coefficients
                << " depth_limit=" << depth_limit
                << " nodes=" << result.subdivision.nodes
                << " leaves=" << result.subdivision.leaves
                << " unresolved=" << result.subdivision.unresolved
                << " maximum_depth=" << result.subdivision.maximum_depth;
      if (result.subdivision.has_first_unresolved) {
        std::cout << " first_unresolved_cell=(";
        for (std::size_t index = 0U;
             index < result.subdivision.first_unresolved_cell.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.subdivision.first_unresolved_cell[index];
        }
        std::cout << ") first_unresolved_splits=(";
        for (std::size_t index = 0U;
             index < result.subdivision.first_unresolved_splits.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.subdivision.first_unresolved_splits[index];
        }
        std::cout << ')';
      }
      std::cout << '\n';
      return EXIT_SUCCESS;
    }
    if (
        argc == 3
        && std::string{argv[1]}
               == "--support-five-rat2-high-ratio-leading-z-zero"
    ) {
      const int depth_limit = parse_positive(argv[2], "depth limit");
      if (depth_limit > 24) {
        throw std::invalid_argument("depth limit must be at most 24");
      }
      const Wall121SupportFiveLargeRatioFaceDiagnostic result =
          wall_121_support_five_large_ratio_leading_z_zero_diagnostic(
              depth_limit);
      std::cout << "SU2_WALL_121_RAT2_SUPPORT_FIVE_HIGH_RATIO_LEADING_Z_ZERO"
                << " t_factor=" << result.t_factor
                << " degrees=(" << result.degrees[0] << ','
                << result.degrees[1] << ')'
                << " terms=" << result.terms
                << " bernstein_coefficients=" << result.coefficients
                << " initial_negative=" << result.negative_coefficients
                << " depth_limit=" << depth_limit
                << " nodes=" << result.subdivision.nodes
                << " leaves=" << result.subdivision.leaves
                << " unresolved=" << result.subdivision.unresolved
                << " maximum_depth=" << result.subdivision.maximum_depth;
      if (result.subdivision.has_first_unresolved) {
        std::cout << " first_unresolved_cell=(";
        for (std::size_t index = 0U;
             index < result.subdivision.first_unresolved_cell.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.subdivision.first_unresolved_cell[index];
        }
        std::cout << ") first_unresolved_splits=(";
        for (std::size_t index = 0U;
             index < result.subdivision.first_unresolved_splits.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.subdivision.first_unresolved_splits[index];
        }
        std::cout << ')';
      }
      std::cout << '\n';
      return EXIT_SUCCESS;
    }
    if (
        argc == 3
        && std::string{argv[1]}
               == "--support-five-rat2-high-ratio-diagonal-corner"
    ) {
      const int depth_limit = parse_positive(argv[2], "depth limit");
      if (depth_limit > 24) {
        throw std::invalid_argument("depth limit must be at most 24");
      }
      const Wall121SupportFiveLargeRatioFaceDiagnostic result =
          wall_121_support_five_large_ratio_diagonal_corner_diagnostic(
              depth_limit);
      std::cout << "SU2_WALL_121_RAT2_SUPPORT_FIVE_HIGH_RATIO_DIAGONAL_CORNER"
                << " t_factor=" << result.t_factor
                << " degrees=(" << result.degrees[0] << ','
                << result.degrees[1] << ',' << result.degrees[2] << ')'
                << " terms=" << result.terms
                << " bernstein_coefficients=" << result.coefficients
                << " initial_negative=" << result.negative_coefficients
                << " depth_limit=" << depth_limit
                << " nodes=" << result.subdivision.nodes
                << " leaves=" << result.subdivision.leaves
                << " unresolved=" << result.subdivision.unresolved
                << " maximum_depth=" << result.subdivision.maximum_depth
                << " lower_endpoint_positive="
                << (result.lower_endpoint_positive ? 1 : 0)
                << " upper_endpoint_positive="
                << (result.upper_endpoint_positive ? 1 : 0);
      if (result.subdivision.has_first_unresolved) {
        std::cout << " first_unresolved_cell=(";
        for (std::size_t index = 0U;
             index < result.subdivision.first_unresolved_cell.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.subdivision.first_unresolved_cell[index];
        }
        std::cout << ") first_unresolved_splits=(";
        for (std::size_t index = 0U;
             index < result.subdivision.first_unresolved_splits.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.subdivision.first_unresolved_splits[index];
        }
        std::cout << ')';
      }
      std::cout << '\n';
      return EXIT_SUCCESS;
    }
    if (
        argc == 3
        && std::string{argv[1]}
               == "--support-five-rat2-high-ratio-complementary-corner"
    ) {
      const int depth_limit = parse_positive(argv[2], "depth limit");
      if (depth_limit > 24) {
        throw std::invalid_argument("depth limit must be at most 24");
      }
      const Wall121SupportFiveLargeRatioFaceDiagnostic result =
          wall_121_support_five_large_ratio_complementary_corner_diagnostic(
              depth_limit);
      std::cout
          << "SU2_WALL_121_RAT2_SUPPORT_FIVE_HIGH_RATIO_COMPLEMENTARY_CORNER"
          << " total_factor=" << result.t_factor
          << " degrees=(" << result.degrees[0] << ','
          << result.degrees[1] << ',' << result.degrees[2] << ')'
          << " terms=" << result.terms
          << " bernstein_coefficients=" << result.coefficients
          << " initial_negative=" << result.negative_coefficients
          << " depth_limit=" << depth_limit
          << " nodes=" << result.subdivision.nodes
          << " leaves=" << result.subdivision.leaves
          << " unresolved=" << result.subdivision.unresolved
          << " maximum_depth=" << result.subdivision.maximum_depth
          << " lower_endpoint_positive="
          << (result.lower_endpoint_positive ? 1 : 0)
          << " upper_endpoint_positive="
          << (result.upper_endpoint_positive ? 1 : 0);
      if (result.subdivision.has_first_unresolved) {
        std::cout << " first_unresolved_cell=(";
        for (std::size_t index = 0U;
             index < result.subdivision.first_unresolved_cell.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.subdivision.first_unresolved_cell[index];
        }
        std::cout << ") first_unresolved_splits=(";
        for (std::size_t index = 0U;
             index < result.subdivision.first_unresolved_splits.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.subdivision.first_unresolved_splits[index];
        }
        std::cout << ')';
      }
      std::cout << '\n';
      return EXIT_SUCCESS;
    }
    if (
        argc == 3
        && std::string{argv[1]} == "--support-five-rat2-subdivide"
    ) {
      const int depth_limit = parse_positive(argv[2], "depth limit");
      if (depth_limit > 12) {
        throw std::invalid_argument("depth limit must be at most 12");
      }
      std::size_t terms = 0U;
      const NDBernsteinGrid grid =
          wall_121_second_iterate_support_five_grid(terms);
      const std::size_t negative = static_cast<std::size_t>(std::count_if(
          grid.values.begin(), grid.values.end(),
          [](const Rational& value) { return value < 0; }));
      NDSubdivisionResult result;
      nd_certify_subdivision(grid, 0, depth_limit, result);
      std::cout << "SU2_WALL_121_RAT2_SUPPORT_FIVE_SUBDIVISION"
                << " terms=" << terms
                << " degrees=(" << grid.degrees[0] << ','
                << grid.degrees[1] << ',' << grid.degrees[2] << ','
                << grid.degrees[3] << ')'
                << " bernstein_coefficients=" << grid.values.size()
                << " initial_negative=" << negative
                << " depth_limit=" << depth_limit
                << " nodes=" << result.nodes
                << " leaves=" << result.leaves
                << " unresolved=" << result.unresolved
                << " maximum_depth=" << result.maximum_depth
                << " result="
                << (result.unresolved == 0U ? "PASS_EXACT" : "INCOMPLETE")
                << '\n';
      if (result.has_first_unresolved) {
        std::cout << "SU2_WALL_121_RAT2_SUPPORT_FIVE_FIRST_UNRESOLVED"
                  << " cell=(";
        for (std::size_t index = 0U;
             index < result.first_unresolved_cell.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.first_unresolved_cell[index];
        }
        std::cout << ") splits=(";
        for (std::size_t index = 0U;
             index < result.first_unresolved_splits.size();
             ++index) {
          if (index != 0U) {
            std::cout << ',';
          }
          std::cout << result.first_unresolved_splits[index];
        }
        std::cout << ")\n";
      }
      return EXIT_SUCCESS;
    }
    if (
        argc == 2
        && std::string{argv[1]} == "--support-five-rat2-diagnostic"
    ) {
      std::size_t terms = 0U;
      const NDBernsteinGrid grid =
          wall_121_second_iterate_support_five_grid(terms);
      const std::size_t negative = static_cast<std::size_t>(std::count_if(
          grid.values.begin(), grid.values.end(),
          [](const Rational& value) { return value < 0; }));
      std::cout << "SU2_WALL_121_RAT2_SUPPORT_FIVE_DIAGNOSTIC"
                << " terms=" << terms
                << " degrees=(" << grid.degrees[0] << ','
                << grid.degrees[1] << ',' << grid.degrees[2] << ','
                << grid.degrees[3] << ')'
                << " bernstein_coefficients=" << grid.values.size()
                << " initial_negative=" << negative << '\n';
      return EXIT_SUCCESS;
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
        argc == 5
        && std::string{argv[1]} == "--ratio-cube-bernstein"
    ) {
      const int support = parse_positive(argv[2], "support");
      const int depth_limit = parse_positive(argv[3], "depth limit");
      const int maximum_antidiagonal =
          parse_positive(argv[4], "maximum antidiagonal");
      return analyze_ratio_cube_bernstein(
          support, depth_limit, maximum_antidiagonal);
    }
    if (
        argc == 6
        && std::string{argv[1]} == "--ratio-cube-best"
    ) {
      const int support = parse_positive(argv[2], "support");
      const int depth_limit = parse_positive(argv[3], "depth limit");
      const int antidiagonal = parse_positive(argv[4], "antidiagonal");
      const int depth = parse_nonnegative(argv[5], "radial depth");
      return analyze_ratio_cube_best_subdivision(
          support, depth_limit, antidiagonal, depth);
    }
    if (
        argc == 7
        && std::string{argv[1]} == "--ratio-cube-local-corner"
    ) {
      const int support = parse_positive(argv[2], "support");
      const int antidiagonal = parse_positive(argv[3], "antidiagonal");
      const int depth = parse_nonnegative(argv[4], "radial depth");
      return analyze_ratio_cube_local_corner(
          support, antidiagonal, depth,
          parse_csv_nonnegative(argv[5], "cell"),
          parse_csv_nonnegative(argv[6], "splits"));
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
