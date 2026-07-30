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
