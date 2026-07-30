#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
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

Polynomial add(const Polynomial& left, const Polynomial& right,
               const Integer& right_scale = 1) {
  Polynomial result(std::max(left.size(), right.size()));
  for (std::size_t index = 0; index < left.size(); ++index) {
    result[index] += left[index];
  }
  for (std::size_t index = 0; index < right.size(); ++index) {
    result[index] += right_scale * right[index];
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

Polynomial linear(const int constant) {
  return {Integer(constant), Integer(1)};
}

Integer evaluate(const Polynomial& polynomial, const Integer& value) {
  Integer result = 0;
  for (std::size_t index = polynomial.size(); index > 0U; --index) {
    result = result * value + polynomial[index - 1U];
  }
  return result;
}

Polynomial divide_linear(const Polynomial& polynomial, const int root) {
  Polynomial quotient(polynomial.size() - 1U);
  Integer carry = polynomial.back();
  quotient.back() = carry;
  for (std::size_t degree = polynomial.size() - 1U;
       degree > 1U; --degree) {
    carry = polynomial[degree - 1U] + root * carry;
    quotient[degree - 2U] = carry;
  }
  if (polynomial.front() + root * carry != 0) {
    throw std::runtime_error("nonzero synthetic-division remainder");
  }
  trim(quotient);
  return quotient;
}

Polynomial ballot_numerator(const int label, const int maximum_label) {
  Polynomial result{Integer(2 * label + 1)};
  for (int index = 0; index < label; ++index) {
    result = multiply(result, linear(-index));
  }
  for (int index = label + 2; index <= maximum_label + 1; ++index) {
    result = multiply(result, linear(index));
  }
  return result;
}

int fusion_count(const int factor, const int source, const int target) {
  return std::max(
      0, std::min(factor, source + target) -
             std::abs(source - target) + 1);
}

struct Check {
  int factor = 0;
  int column = 0;
  int degree = 0;
  std::vector<int> roots;
  Integer minimum_coefficient = 0;
};

}  // namespace

int main() {
  try {
    std::vector<Check> checks;
    std::uint64_t denominator_checks = 0U;
    std::uint64_t denominator_failures = 0U;
    std::uint64_t coefficient_checks = 0U;
    std::uint64_t coefficient_failures = 0U;
    for (int factor = 2; factor <= 9; ++factor) {
      const int maximum_label = factor + 5;
      Polynomial denominator{Integer(1)};
      for (int index = 2; index <= maximum_label + 1; ++index) {
        denominator = multiply(denominator, linear(index));
      }
      ++denominator_checks;
      if (denominator.size() !=
          static_cast<std::size_t>(maximum_label + 1)) {
        ++denominator_failures;
      }
      std::vector<Polynomial> ballot;
      for (int label = 0; label <= maximum_label; ++label) {
        ballot.push_back(ballot_numerator(label, maximum_label));
      }
      std::vector<Polynomial> coefficient;
      for (int target = 0; target <= 5; ++target) {
        Polynomial value{Integer(0)};
        for (int source = 0; source <= factor + target; ++source) {
          value = add(
              value, ballot[static_cast<std::size_t>(source)],
              fusion_count(factor, source, target));
        }
        coefficient.push_back(value);
      }
      const Polynomial third_zero = add(
          multiply(
              coefficient[0],
              add(add(coefficient[3], coefficient[4]), coefficient[5])),
          multiply(coefficient[1], coefficient[4]), -1);
      const Polynomial third_one = add(
          add(
              multiply(
                  coefficient[0],
                  add(coefficient[4], coefficient[5])),
              multiply(coefficient[1], coefficient[2])),
          multiply(coefficient[2], coefficient[3]), -1);

      for (int column = 0; column < 2; ++column) {
        Polynomial quotient = column == 0 ? third_zero : third_one;
        std::vector<int> roots;
        for (int root = 0; root <= maximum_label; ++root) {
          while (quotient.size() > 1U &&
                 evaluate(quotient, root) == 0) {
            roots.push_back(root);
            quotient = divide_linear(quotient, root);
          }
        }
        Integer minimum = quotient.front();
        for (const Integer& coefficient_value : quotient) {
          ++coefficient_checks;
          minimum = std::min(minimum, coefficient_value);
          if (coefficient_value <= 0) {
            ++coefficient_failures;
          }
        }
        checks.push_back(
            {factor, column, static_cast<int>(quotient.size()) - 1,
             roots, minimum});
      }
    }

    const bool pass =
        checks.size() == 16U &&
        denominator_checks == 8U &&
        denominator_failures == 0U &&
        coefficient_checks == 318U &&
        coefficient_failures == 0U;
    std::cout
        << "SU2_FUNDAMENTAL_ONE_ARBITRARY_THIRD"
        << " polynomial_checks=" << checks.size()
        << " denominator_checks=" << denominator_checks
        << " denominator_failures=" << denominator_failures
        << " coefficient_checks=" << coefficient_checks
        << " coefficient_failures=" << coefficient_failures
        << " summaries={";
    for (std::size_t index = 0; index < checks.size(); ++index) {
      if (index != 0U) {
        std::cout << ';';
      }
      const Check& check = checks[index];
      std::cout
          << check.factor << ',' << check.column
          << ":degree=" << check.degree << ",roots=[";
      for (std::size_t root = 0; root < check.roots.size(); ++root) {
        if (root != 0U) {
          std::cout << ',';
        }
        std::cout << check.roots[root];
      }
      std::cout
          << "],minimum=" << check.minimum_coefficient;
    }
    std::cout
        << "} result=" << (pass ? "PASS" : "FAIL")
        << '\n';
    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
