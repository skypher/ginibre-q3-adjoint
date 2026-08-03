#include <array>
#include <iostream>
#include <map>
#include <stdexcept>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Exponent = std::array<int, 4>;
using Polynomial = std::map<Exponent, Integer>;

Polynomial constant(const long long value) {
  if (value == 0L) {
    return {};
  }
  return {{{0, 0, 0, 0}, Integer{value}}};
}

Polynomial variable(const int coordinate) {
  Exponent exponent{0, 0, 0, 0};
  exponent[static_cast<std::size_t>(coordinate)] = 1;
  return {{exponent, Integer{1}}};
}

void add_scaled(Polynomial& target, const Polynomial& source,
                const Integer& scale) {
  for (const auto& [exponent, coefficient] : source) {
    Integer& entry = target[exponent];
    entry += scale * coefficient;
    if (entry == 0) {
      target.erase(exponent);
    }
  }
}

Polynomial add(const Polynomial& left, const Polynomial& right) {
  Polynomial result = left;
  add_scaled(result, right, Integer{1});
  return result;
}

Polynomial subtract(const Polynomial& left, const Polynomial& right) {
  Polynomial result = left;
  add_scaled(result, right, Integer{-1});
  return result;
}

Polynomial multiply(const Polynomial& left, const Polynomial& right) {
  Polynomial result;
  for (const auto& [left_exponent, left_coefficient] : left) {
    for (const auto& [right_exponent, right_coefficient] : right) {
      Exponent exponent{0, 0, 0, 0};
      for (std::size_t coordinate = 0U; coordinate < exponent.size();
           ++coordinate) {
        exponent[coordinate] = left_exponent[coordinate]
                               + right_exponent[coordinate];
      }
      result[exponent] += left_coefficient * right_coefficient;
    }
  }
  return result;
}

Polynomial power(Polynomial base, int exponent) {
  if (exponent < 0) {
    throw std::invalid_argument("negative exponent");
  }
  Polynomial result = constant(1);
  while (exponent > 0) {
    if ((exponent % 2) != 0) {
      result = multiply(result, base);
    }
    base = multiply(base, base);
    exponent /= 2;
  }
  return result;
}

Polynomial scaled(const Polynomial& input, const long long coefficient) {
  Polynomial result;
  add_scaled(result, input, Integer{coefficient});
  return result;
}

Integer integer_power(const int base, const int exponent) {
  Integer result = 1;
  for (int index = 0; index < exponent; ++index) {
    result *= base;
  }
  return result;
}

Integer evaluate(const Polynomial& input, const std::array<int, 4>& point) {
  Integer result = 0;
  for (const auto& [exponent, coefficient] : input) {
    Integer monomial = coefficient;
    for (std::size_t coordinate = 0U; coordinate < point.size();
         ++coordinate) {
      monomial *= integer_power(
          point[coordinate], exponent[coordinate]);
    }
    result += monomial;
  }
  return result;
}

}  // namespace

int main() {
  // q=1+a, t=q(1+b), s=t(1+c), r=s(1+d) parameterizes
  // exactly r>=s>=t>=q>=1 by a,b,c,d>=0.
  const Polynomial one = constant(1);
  const Polynomial q = add(one, variable(0));
  const Polynomial t = multiply(q, add(one, variable(1)));
  const Polynomial s = multiply(t, add(one, variable(2)));
  const Polynomial r = multiply(s, add(one, variable(3)));
  const Polynomial c = multiply(r, s);
  const Polynomial d = multiply(c, t);
  const Polynomial x = multiply(d, q);

  Polynomial linear = power(d, 3);
  linear = subtract(linear, multiply(c, d));
  linear = subtract(linear, scaled(multiply(c, power(d, 2)), 4));
  linear = subtract(linear, scaled(multiply(power(c, 2), d), 2));
  linear = add(linear, scaled(power(c, 3), 2));
  linear = subtract(linear, multiply(r, c));
  linear = subtract(linear, scaled(multiply(multiply(r, c), d), 2));
  linear = add(linear, multiply(power(r, 2), d));
  linear = add(linear, power(r, 3));

  Polynomial quadratic = scaled(power(d, 2), -2);
  quadratic = subtract(quadratic, scaled(multiply(c, d), 2));
  quadratic = add(quadratic, scaled(power(c, 2), 2));
  quadratic = add(quadratic, multiply(r, c));

  Polynomial bracket = add(linear, multiply(quadratic, x));
  bracket = add(
      bracket,
      multiply(add(add(r, scaled(c, 2)), d), power(x, 2)));
  bracket = add(bracket, scaled(power(x, 3), 2));

  std::size_t negative = 0U;
  Integer minimum = 0;
  bool have_minimum = false;
  for (const auto& [exponent, coefficient] : bracket) {
    static_cast<void>(exponent);
    if (!have_minimum || coefficient < minimum) {
      minimum = coefficient;
      have_minimum = true;
    }
    if (coefficient < 0) {
      ++negative;
    }
  }
  Integer grid_minimum = 0;
  std::array<int, 4> grid_minimizer{0, 0, 0, 0};
  bool have_grid_minimum = false;
  for (int a = 0; a <= 4; ++a) {
    for (int b = 0; b <= 4; ++b) {
      for (int c_slack = 0; c_slack <= 4; ++c_slack) {
        for (int d_slack = 0; d_slack <= 4; ++d_slack) {
          const std::array<int, 4> point{a, b, c_slack, d_slack};
          const Integer value = evaluate(bracket, point);
          if (!have_grid_minimum || value < grid_minimum) {
            grid_minimum = value;
            grid_minimizer = point;
            have_grid_minimum = true;
          }
        }
      }
    }
  }
  const Integer root_a = 1;
  const Integer root_b = 10;
  const Integer root_c = 100;
  const Integer root_d = 150;
  const Integer root_x = 150;
  const Integer witness_linear =
      root_d * root_d * root_d - root_a * root_c * root_d
      - 4 * root_c * root_d * root_d - 2 * root_c * root_c * root_d
      + 2 * root_c * root_c * root_c - root_a * root_b * root_c
      - 2 * root_b * root_c * root_d + root_b * root_b * root_d
      + root_b * root_b * root_b;
  const Integer witness_quadratic =
      -2 * root_d * root_d - 2 * root_c * root_d
      + 2 * root_c * root_c + root_b * root_c;
  const Integer witness_bracket =
      witness_linear + witness_quadratic * root_x
      + (root_b + 2 * root_c + root_d) * root_x * root_x
      + 2 * root_x * root_x * root_x;
  const Integer witness_delta = root_x * witness_bracket;
  std::cout << "SU2_WALL_121_APPEND_Q_GE_ONE"
            << " monomials=" << bracket.size()
            << " negative_coefficients=" << negative
            << " minimum_coefficient=" << minimum
            << " grid_minimum=" << grid_minimum << "@(" << grid_minimizer[0]
            << ',' << grid_minimizer[1] << ',' << grid_minimizer[2] << ','
            << grid_minimizer[3] << ')'
            << " coefficientwise=" << (negative == 0U ? "PASS" : "FAIL")
            << " grid=" << (grid_minimum >= 0 ? "PASS" : "COUNTEREXAMPLE")
            << '\n';
  std::cout << "SU2_WALL_121_APPEND_Q_GE_ONE_COUNTEREXAMPLE"
            << " root=[1,10,100,150,150]"
            << " ratios=[10,10,3/2,1]"
            << " bracket=" << witness_bracket
            << " delta=" << witness_delta
            << " result="
            << ((witness_bracket < 0 && witness_delta < 0) ? "PASS" : "FAIL")
            << '\n';
  return grid_minimum >= 0 ? 0 : 1;
}
