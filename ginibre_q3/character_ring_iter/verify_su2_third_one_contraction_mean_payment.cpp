#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;

struct Field {
  Rational rational{Integer(0)};
  Rational radical{Integer(0)};
};

Field operator+(const Field& left, const Field& right) {
  return {
      left.rational + right.rational,
      left.radical + right.radical};
}

Field operator-(const Field& left, const Field& right) {
  return {
      left.rational - right.rational,
      left.radical - right.radical};
}

Field operator-(const Field& value) {
  return {-value.rational, -value.radical};
}

Field operator*(const Field& left, const Field& right) {
  return {
      left.rational * right.rational +
          Rational(Integer(10)) * left.radical * right.radical,
      left.rational * right.radical +
          left.radical * right.rational};
}

Field operator*(const Field& value, const Rational& scale) {
  return {value.rational * scale, value.radical * scale};
}

int sign(const Field& value) {
  if (value.radical == 0) {
    return value.rational > 0 ? 1 : (value.rational < 0 ? -1 : 0);
  }
  if (value.rational == 0) {
    return value.radical > 0 ? 1 : -1;
  }
  if ((value.rational > 0) == (value.radical > 0)) {
    return value.rational > 0 ? 1 : -1;
  }
  const Rational comparison =
      value.rational * value.rational -
      Rational(Integer(10)) * value.radical * value.radical;
  if (comparison == 0) {
    return 0;
  }
  if (value.rational > 0) {
    return comparison > 0 ? 1 : -1;
  }
  return comparison > 0 ? -1 : 1;
}

struct Polynomial {
  std::vector<std::vector<Field>> coefficients;

  Polynomial() : coefficients(1, std::vector<Field>(1)) {}
  explicit Polynomial(const Field& value)
      : coefficients(1, std::vector<Field>(1, value)) {}
};

void trim(Polynomial& polynomial) {
  while (polynomial.coefficients.size() > 1U) {
    bool zero = true;
    for (const Field& coefficient : polynomial.coefficients.back()) {
      zero = zero && sign(coefficient) == 0;
    }
    if (!zero) {
      break;
    }
    polynomial.coefficients.pop_back();
  }
  std::size_t columns = 1U;
  for (const auto& row : polynomial.coefficients) {
    for (std::size_t column = row.size(); column > 0U; --column) {
      if (sign(row[column - 1U]) != 0) {
        columns = std::max(columns, column);
        break;
      }
    }
  }
  for (auto& row : polynomial.coefficients) {
    row.resize(columns);
  }
}

Polynomial constant(const int value) {
  return Polynomial(Field{Rational(Integer(value)), Rational(Integer(0))});
}

Polynomial monomial(
    const Field& constant_term, const Field& first,
    const Field& second) {
  Polynomial result;
  result.coefficients.assign(2U, std::vector<Field>(2U));
  result.coefficients[0][0] = constant_term;
  result.coefficients[1][0] = first;
  result.coefficients[0][1] = second;
  trim(result);
  return result;
}

Polynomial operator+(
    const Polynomial& left, const Polynomial& right) {
  Polynomial result;
  result.coefficients.assign(
      std::max(left.coefficients.size(), right.coefficients.size()),
      std::vector<Field>(
          std::max(
              left.coefficients.front().size(),
              right.coefficients.front().size())));
  for (std::size_t i = 0; i < left.coefficients.size(); ++i) {
    for (std::size_t j = 0; j < left.coefficients[i].size(); ++j) {
      result.coefficients[i][j] =
          result.coefficients[i][j] + left.coefficients[i][j];
    }
  }
  for (std::size_t i = 0; i < right.coefficients.size(); ++i) {
    for (std::size_t j = 0; j < right.coefficients[i].size(); ++j) {
      result.coefficients[i][j] =
          result.coefficients[i][j] + right.coefficients[i][j];
    }
  }
  trim(result);
  return result;
}

Polynomial operator-(
    const Polynomial& left, const Polynomial& right) {
  Polynomial negative = right;
  for (auto& row : negative.coefficients) {
    for (Field& coefficient : row) {
      coefficient = -coefficient;
    }
  }
  return left + negative;
}

Polynomial operator*(
    const Polynomial& polynomial, const Rational& scale) {
  Polynomial result = polynomial;
  for (auto& row : result.coefficients) {
    for (Field& coefficient : row) {
      coefficient = coefficient * scale;
    }
  }
  trim(result);
  return result;
}

Polynomial operator*(
    const Polynomial& left, const Polynomial& right) {
  Polynomial result;
  result.coefficients.assign(
      left.coefficients.size() + right.coefficients.size() - 1U,
      std::vector<Field>(
          left.coefficients.front().size() +
          right.coefficients.front().size() - 1U));
  for (std::size_t i = 0; i < left.coefficients.size(); ++i) {
    for (std::size_t j = 0; j < left.coefficients[i].size(); ++j) {
      if (sign(left.coefficients[i][j]) == 0) {
        continue;
      }
      for (std::size_t k = 0; k < right.coefficients.size(); ++k) {
        for (std::size_t l = 0; l < right.coefficients[k].size();
             ++l) {
          result.coefficients[i + k][j + l] =
              result.coefficients[i + k][j + l] +
              left.coefficients[i][j] * right.coefficients[k][l];
        }
      }
    }
  }
  trim(result);
  return result;
}

Polynomial kernel(const Polynomial& x, const Polynomial& y) {
  const Polynomial difference = x - y;
  return difference * difference *
         (x + y - constant(1)) *
         (x * x + x * y + y * y -
          x * Rational(Integer(2)) -
          y * Rational(Integer(2)) - constant(2));
}

std::vector<Integer> multiply_by_square(
    const std::vector<Integer>& profile, const int label) {
  const int old_support = static_cast<int>(profile.size()) - 1;
  std::vector<Integer> result(
      static_cast<std::size_t>(old_support + label + 1));
  for (int source = 0; source <= old_support; ++source) {
    for (int shell = 0; shell <= label; ++shell) {
      for (int target = std::abs(source - shell);
           target <= source + shell; ++target) {
        result[static_cast<std::size_t>(target)] +=
            profile[static_cast<std::size_t>(source)];
      }
    }
  }
  return result;
}

Integer binomial(const std::size_t n, const std::size_t k) {
  Integer result = 1;
  const std::size_t order = std::min(k, n - k);
  for (std::size_t index = 1; index <= order; ++index) {
    result *= n - order + index;
    result /= index;
  }
  return result;
}

struct Check {
  std::size_t coefficients = 0U;
  std::size_t positive = 0U;
  std::size_t zero = 0U;
  std::size_t negative = 0U;
};

using Bernstein = std::vector<std::vector<Field>>;

Bernstein bernstein_coefficients(const Polynomial& polynomial) {
  const std::size_t first_degree =
      polynomial.coefficients.size() - 1U;
  const std::size_t second_degree =
      polynomial.coefficients.front().size() - 1U;
  Bernstein result(
      first_degree + 1U,
      std::vector<Field>(second_degree + 1U));
  for (std::size_t first = 0; first <= first_degree; ++first) {
    for (std::size_t second = 0; second <= second_degree; ++second) {
      Field coefficient;
      for (std::size_t i = 0; i <= first; ++i) {
        for (std::size_t j = 0; j <= second; ++j) {
          coefficient = coefficient +
              polynomial.coefficients[i][j] *
              Rational(
                  binomial(first, i) * binomial(second, j),
                  binomial(first_degree, i) *
                      binomial(second_degree, j));
        }
      }
      result[first][second] = coefficient;
    }
  }
  return result;
}

Check check_bernstein(const Bernstein& coefficients) {
  Check result;
  for (const auto& row : coefficients) {
    for (const Field& coefficient : row) {
      ++result.coefficients;
      const int coefficient_sign = sign(coefficient);
      if (coefficient_sign > 0) {
        ++result.positive;
      } else if (coefficient_sign < 0) {
        ++result.negative;
      } else {
        ++result.zero;
      }
    }
  }
  return result;
}

std::pair<Bernstein, Bernstein> split_first(
    const Bernstein& coefficients) {
  const std::size_t degree = coefficients.size() - 1U;
  const std::size_t columns = coefficients.front().size();
  Bernstein left(degree + 1U, std::vector<Field>(columns));
  Bernstein right(degree + 1U, std::vector<Field>(columns));
  for (std::size_t column = 0; column < columns; ++column) {
    std::vector<Field> work(degree + 1U);
    for (std::size_t index = 0; index <= degree; ++index) {
      work[index] = coefficients[index][column];
    }
    left[0][column] = work[0];
    right[degree][column] = work[degree];
    for (std::size_t level = 1; level <= degree; ++level) {
      for (std::size_t index = 0; index <= degree - level;
           ++index) {
        work[index] =
            (work[index] + work[index + 1U]) *
            Rational(Integer(1), Integer(2));
      }
      left[level][column] = work[0];
      right[degree - level][column] = work[degree - level];
    }
  }
  return {std::move(left), std::move(right)};
}

Bernstein transpose(const Bernstein& coefficients) {
  Bernstein result(
      coefficients.front().size(),
      std::vector<Field>(coefficients.size()));
  for (std::size_t first = 0; first < coefficients.size(); ++first) {
    for (std::size_t second = 0;
         second < coefficients[first].size(); ++second) {
      result[second][first] = coefficients[first][second];
    }
  }
  return result;
}

std::pair<Bernstein, Bernstein> split_second(
    const Bernstein& coefficients) {
  auto result = split_first(transpose(coefficients));
  return {transpose(result.first), transpose(result.second)};
}

struct Certificate {
  std::uint64_t visited = 0U;
  std::uint64_t leaves = 0U;
  std::uint64_t unresolved = 0U;
  int maximum_depth = 0;
};

struct Expected {
  Check root;
  std::uint64_t visited = 0U;
  std::uint64_t leaves = 0U;
  int maximum_depth = 0;
};

bool certify(
    const Bernstein& coefficients, const int depth,
    const int maximum_depth, Certificate& certificate) {
  ++certificate.visited;
  certificate.maximum_depth =
      std::max(certificate.maximum_depth, depth);
  const Check check = check_bernstein(coefficients);
  if (check.negative == 0U) {
    ++certificate.leaves;
    return true;
  }
  if (depth == maximum_depth) {
    ++certificate.unresolved;
    return false;
  }
  const auto children =
      depth % 2 == 0 ? split_first(coefficients)
                     : split_second(coefficients);
  const bool first = certify(
      children.first, depth + 1, maximum_depth, certificate);
  const bool second = certify(
      children.second, depth + 1, maximum_depth, certificate);
  return first && second;
}

bool verify(
    const std::string& name, const Polynomial& polynomial,
    const Expected& expected) {
  const Bernstein root = bernstein_coefficients(polynomial);
  const Check root_check = check_bernstein(root);
  Certificate certificate;
  const bool certified = certify(root, 0, 8, certificate);
  const bool pass =
      certified &&
      root_check.coefficients == expected.root.coefficients &&
      root_check.positive == expected.root.positive &&
      root_check.zero == expected.root.zero &&
      root_check.negative == expected.root.negative &&
      certificate.visited == expected.visited &&
      certificate.leaves == expected.leaves &&
      certificate.unresolved == 0U &&
      certificate.maximum_depth == expected.maximum_depth;
  std::cout
      << name
      << " root_coefficients=" << root_check.coefficients
      << " root_positive=" << root_check.positive
      << " root_zero=" << root_check.zero
      << " root_negative=" << root_check.negative
      << " visited=" << certificate.visited
      << " leaves=" << certificate.leaves
      << " unresolved=" << certificate.unresolved
      << " maximum_depth=" << certificate.maximum_depth
      << " pass=" << pass << '\n';
  return pass;
}

}  // namespace

int main() {
  const Field alpha{
      Rational(Integer(2), Integer(3)),
      Rational(Integer(2), Integer(3))};
  const Field lower_critical{
      Rational(Integer(2), Integer(3)),
      Rational(Integer(-1), Integer(3))};
  const Field three{Rational(Integer(3)), Rational(Integer(0))};
  const Field two{Rational(Integer(2)), Rational(Integer(0))};
  const Field minus_one{Rational(Integer(-1)), Rational(Integer(0))};
  std::uint64_t scalar_checks = 0U;
  std::uint64_t scalar_failures = 0U;
  const auto check_zero = [&](const Field& value) {
    ++scalar_checks;
    if (sign(value) != 0) {
      ++scalar_failures;
    }
  };
  const auto check_positive = [&](const Field& value) {
    ++scalar_checks;
    if (sign(value) <= 0) {
      ++scalar_failures;
    }
  };
  const auto check_integer = [&](const Integer& left,
                                 const Integer& right) {
    ++scalar_checks;
    if (left != right) {
      ++scalar_failures;
    }
  };
  check_positive(
      Field{Rational(Integer(37), Integer(13)), Rational(Integer(0))} -
      Field{Rational(Integer(48), Integer(17)), Rational(Integer(0))});
  check_zero(
      lower_critical * lower_critical * Rational(Integer(3)) -
      lower_critical * Rational(Integer(4)) -
      Field{Rational(Integer(2)), Rational(Integer(0))});
  check_positive(alpha - two);
  check_positive(three - alpha);
  check_positive(lower_critical - minus_one);
  check_positive(two - lower_critical);
  check_positive(
      Field{Rational(Integer(48), Integer(17)), Rational(Integer(0))} -
      alpha);
  check_zero(
      alpha * alpha * Rational(Integer(3)) -
      alpha * Rational(Integer(4)) -
      Field{Rational(Integer(12)), Rational(Integer(0))});
  std::vector<Integer> countercontrol{Integer(1)};
  for (int factor = 0; factor < 8; ++factor) {
    countercontrol = multiply_by_square(countercontrol, 2);
  }
  check_integer(countercontrol[0], Integer(227475));
  check_integer(countercontrol[1], Integer(625992));
  const Integer threshold_left =
      3 * countercontrol[1] - 2 * countercontrol[0];
  check_integer(
      threshold_left * threshold_left -
          40 * countercontrol[0] * countercontrol[0],
      Integer(-44792028324LL));
  const Polynomial alpha_polynomial(alpha);
  const Polynomial low_x =
      monomial(minus_one, three, Field{});
  const Polynomial low_x_left =
      monomial(
          minus_one, lower_critical - minus_one, Field{});
  const Polynomial low_x_right =
      monomial(
          lower_critical, two - lower_critical, Field{});
  const Polynomial low_y =
      monomial(minus_one, Field{}, three);
  const Polynomial high_y =
      monomial(two, Field{}, alpha - two);
  const Polynomial cross_z =
      monomial(alpha, Field{}, three - alpha);

  const auto low_payment = [&](const Polynomial& x,
                               const Polynomial& y) {
    return kernel(x, y) +
        (alpha_polynomial - x) *
            (alpha_polynomial - y) *
            Rational(Integer(16));
  };
  const auto cross_payment = [&](const Polynomial& x,
                                 const Polynomial& z) {
    return kernel(x, z) -
        (alpha_polynomial - x) *
            (z - alpha_polynomial) *
            Rational(Integer(8));
  };

  const bool low_low_pass =
      verify(
          "low_low", low_payment(low_x, low_y),
          Expected{Check{36U, 36U, 0U, 0U}, 1U, 1U, 0});
  const bool low_high_left_pass =
      verify(
          "low_high_left",
          low_payment(low_x_left, high_y),
          Expected{Check{36U, 34U, 2U, 0U}, 1U, 1U, 0});
  const bool low_high_right_pass =
      verify(
          "low_high_right",
          low_payment(low_x_right, high_y),
          Expected{Check{36U, 33U, 2U, 1U}, 7U, 4U, 3});
  const bool cross_left_pass =
      verify(
          "cross_left",
          cross_payment(low_x_left, cross_z),
          Expected{Check{36U, 34U, 2U, 0U}, 1U, 1U, 0});
  const bool cross_right_pass =
      verify(
          "cross_right",
          cross_payment(low_x_right, cross_z),
          Expected{Check{36U, 34U, 2U, 0U}, 1U, 1U, 0});
  const bool pass =
      low_low_pass && low_high_left_pass &&
      low_high_right_pass && cross_left_pass &&
      cross_right_pass && scalar_checks == 11U &&
      scalar_failures == 0U;
  std::cout
      << "scalar_checks=" << scalar_checks
      << " scalar_failures=" << scalar_failures
      << " result=" << (pass ? "PASS" : "FAIL") << '\n';
  return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
