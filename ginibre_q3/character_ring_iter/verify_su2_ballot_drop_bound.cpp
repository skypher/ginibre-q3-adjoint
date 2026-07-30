#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Polynomial = std::vector<Integer>;

void trim(Polynomial& polynomial) {
  while (polynomial.size() > 1U && polynomial.back() == 0) {
    polynomial.pop_back();
  }
}

Polynomial add(const Polynomial& left, const Polynomial& right) {
  Polynomial result(std::max(left.size(), right.size()));
  for (std::size_t index = 0; index < result.size(); ++index) {
    if (index < left.size()) {
      result[index] += left[index];
    }
    if (index < right.size()) {
      result[index] += right[index];
    }
  }
  trim(result);
  return result;
}

Polynomial scale(const Polynomial& polynomial, const Integer& factor) {
  Polynomial result = polynomial;
  for (Integer& coefficient : result) {
    coefficient *= factor;
  }
  trim(result);
  return result;
}

Polynomial multiply(const Polynomial& left, const Polynomial& right) {
  Polynomial result(left.size() + right.size() - 1U);
  for (std::size_t first = 0; first < left.size(); ++first) {
    for (std::size_t second = 0; second < right.size(); ++second) {
      result[first + second] += left[first] * right[second];
    }
  }
  trim(result);
  return result;
}

Polynomial power(Polynomial base, int exponent) {
  Polynomial result{Integer(1)};
  while (exponent > 0) {
    if ((exponent & 1) != 0) {
      result = multiply(result, base);
    }
    exponent /= 2;
    if (exponent > 0) {
      base = multiply(base, base);
    }
  }
  return result;
}

Integer binomial(const int size, int selection) {
  selection = std::min(selection, size - selection);
  Integer result = 1;
  for (int index = 1; index <= selection; ++index) {
    result = result * (size - selection + index) / index;
  }
  return result;
}

Integer ballot(const int factors, const int label) {
  if (label > factors) {
    return 0;
  }
  return Integer(2 * label + 1) *
         binomial(2 * factors, factors - label) /
         (factors + label + 1);
}

Polynomial translate(const Polynomial& polynomial, const int offset) {
  Polynomial result(polynomial.size());
  for (std::size_t exponent = 0; exponent < polynomial.size();
       ++exponent) {
    Integer offset_power = 1;
    for (std::size_t target = exponent + 1U; target > 0U; --target) {
      const std::size_t degree = target - 1U;
      result[degree] +=
          polynomial[exponent] *
          binomial(
              static_cast<int>(exponent),
              static_cast<int>(degree)) *
          offset_power;
      offset_power *= offset;
    }
  }
  trim(result);
  return result;
}

void render(const Polynomial& polynomial) {
  std::cout << '[';
  for (std::size_t index = 0; index < polynomial.size(); ++index) {
    if (index != 0U) {
      std::cout << ',';
    }
    std::cout << polynomial[index];
  }
  std::cout << ']';
}

}  // namespace

int main() {
  const Polynomial one{Integer(1)};
  const Polynomial b{Integer(0), Integer(1)};
  const Polynomial b_plus_one{Integer(1), Integer(1)};
  const Polynomial b_plus_two{Integer(2), Integer(1)};
  const Polynomial minimum_prefactor{
      Integer(3), Integer(7), Integer(2)};
  const Polynomial p = scale(b_plus_one, 3);
  const Polynomial q =
      scale(multiply(b, b_plus_two), 2);

  Polynomial truncated_binomial{Integer(0)};
  Polynomial falling{Integer(1)};
  constexpr int factorials[] = {1, 1, 2, 6, 24};
  for (int order = 0; order <= 4; ++order) {
    if (order > 0) {
      falling = multiply(
          falling,
          add(b, Polynomial{Integer(1 - order)}));
    }
    const Polynomial term = scale(
        multiply(
            multiply(falling, power(p, order)),
            power(q, 4 - order)),
        24 / factorials[order]);
    truncated_binomial = add(truncated_binomial, term);
  }

  const Polynomial left = multiply(
      minimum_prefactor,
      truncated_binomial);
  const Polynomial right = scale(
      multiply(
          Polynomial{Integer(3), Integer(4), Integer(2)},
          power(q, 4)),
      -96);
  const Polynomial numerator = add(left, right);
  const Polynomial shifted = translate(numerator, 9);
  const bool polynomial_pass = std::all_of(
      shifted.begin(), shifted.end(),
      [](const Integer& coefficient) { return coefficient >= 0; });
  std::uint64_t finite_checks = 0U;
  std::uint64_t finite_failures = 0U;
  for (int label = 1; label <= 8; ++label) {
    const int maximum_factors =
        2 * label * label + 4 * label;
    for (int factors = label; factors <= maximum_factors; ++factors) {
      ++finite_checks;
      if (ballot(factors, label) -
              ballot(factors, label + 1) >
          ballot(factors, 0)) {
        ++finite_failures;
      }
    }
  }
  const bool pass = polynomial_pass && finite_failures == 0U;

  std::cout << "SU2_BALLOT_DROP_BOUND numerator=";
  render(numerator);
  std::cout << " shifted_at_9=";
  render(shifted);
  std::cout
      << " finite_checks=" << finite_checks
      << " finite_failures=" << finite_failures;
  std::cout << " result=" << (pass ? "PASS" : "FAIL") << '\n';
  return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
