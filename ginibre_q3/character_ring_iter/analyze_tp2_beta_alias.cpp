#include <boost/multiprecision/cpp_int.hpp>

#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

std::vector<cpp_int> fourier_coefficients(int minus_power, int plus_power) {
  const int degree = minus_power + plus_power;
  std::vector<cpp_int> result(static_cast<std::size_t>(2 * degree + 1));
  const int negative_degree = 2 * minus_power;
  const int positive_degree = 2 * plus_power;
  cpp_int previous = 0;
  cpp_int current = 1;
  const cpp_int linear = positive_degree - negative_degree;
  for (int index = 0; index <= 2 * degree; ++index) {
    result[static_cast<std::size_t>(index)] =
        (minus_power & 1) == 0 ? current : -current;
    if (index == 2 * degree) {
      break;
    }
    const cpp_int next_numerator = linear * current +
        (index - 2 * degree - 1) * previous;
    const int next_denominator = index + 1;
    if (next_numerator % next_denominator != 0) {
      throw std::runtime_error("nonintegral Fourier recurrence");
    }
    previous = current;
    current = next_numerator / next_denominator;
  }
  return result;
}

cpp_int coefficient(const std::vector<cpp_int>& values, int exponent) {
  const int degree = (static_cast<int>(values.size()) - 1) / 2;
  if (exponent < -degree || exponent > degree) {
    return 0;
  }
  return values[static_cast<std::size_t>(degree + exponent)];
}

cpp_int cyclic_coefficient(const std::vector<cpp_int>& values, int modulus,
                           int offset) {
  const int degree = (static_cast<int>(values.size()) - 1) / 2;
  const int bound = (degree + std::abs(offset)) / modulus + 2;
  cpp_int result = 0;
  for (int winding = -bound; winding <= bound; ++winding) {
    result += coefficient(values, offset + winding * modulus);
  }
  return result;
}

cpp_int derivative_alias(const std::vector<cpp_int>& values, int modulus,
                         int index) {
  const int degree = (static_cast<int>(values.size()) - 1) / 2;
  const int bound = (degree + std::abs(index)) / modulus + 2;
  cpp_int result = 0;
  for (int winding = -bound; winding <= bound; ++winding) {
    const int exponent = winding * modulus;
    result += exponent * (
        coefficient(values, exponent - 1 - index) +
        coefficient(values, exponent - 1 + index) -
        coefficient(values, exponent + 1 - index) -
        coefficient(values, exponent + 1 + index)
    );
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  const int maximum_rank = argc > 1 ? std::atoi(argv[1]) : 20;
  const int maximum_minus_power = argc > 2 ? std::atoi(argv[2]) : 40;
  const int maximum_plus_power = argc > 3 ? std::atoi(argv[3]) : 40;
  const int maximum_index = argc > 4 ? std::atoi(argv[4]) : 6;
  if (maximum_rank < 2 || maximum_minus_power < 4 ||
      maximum_plus_power < 2 || maximum_index < 0) {
    std::cerr << "usage: analyze_tp2_beta_alias [maximum-rank>=2] "
                 "[maximum-minus-power>=4] [maximum-plus-power>=2] "
                 "[maximum-index>=0]\n";
    return 2;
  }

  std::vector<std::size_t> negative(
      static_cast<std::size_t>(maximum_index + 1));
  std::vector<std::size_t> positive(
      static_cast<std::size_t>(maximum_index + 1));
  std::vector<bool> reported_negative(
      static_cast<std::size_t>(maximum_index + 1));
  std::vector<bool> reported_positive(
      static_cast<std::size_t>(maximum_index + 1));
  std::size_t systems = 0;
  std::size_t identity_checks = 0;
  for (int rank = 2; rank <= maximum_rank; ++rank) {
    const int modulus = 2 * rank + 1;
    for (int plus_power = 2; plus_power <= maximum_plus_power;
         ++plus_power) {
      for (int minus_power = plus_power + 2;
           minus_power <= maximum_minus_power; ++minus_power) {
        const int difference = minus_power - plus_power;
        const int sum = minus_power + plus_power + 1;
        const std::vector<cpp_int> values =
            fourier_coefficients(minus_power, plus_power);
        std::vector<cpp_int> moments(
            static_cast<std::size_t>(maximum_index + 2));
        for (int index = 0; index <= maximum_index + 1; ++index) {
          moments[static_cast<std::size_t>(index)] =
              2 * cyclic_coefficient(values, modulus, index);
        }
        for (int index = 0; index <= maximum_index; ++index) {
          const cpp_int recurrence =
              (sum + index) *
                  moments[static_cast<std::size_t>(index + 1)] +
              2 * difference * moments[static_cast<std::size_t>(index)] +
              (sum - index) * moments[static_cast<std::size_t>(
                  index == 0 ? 1 : index - 1
              )];
          const cpp_int direct = derivative_alias(values, modulus, index);
          ++identity_checks;
          if (recurrence != direct) {
            std::cout << "TP2_BETA_ALIAS result=IDENTITY_FAIL"
                      << " rank=" << rank
                      << " minus_power=" << minus_power
                      << " plus_power=" << plus_power
                      << " index=" << index
                      << " recurrence=" << recurrence
                      << " direct=" << direct << '\n';
            return 1;
          }
          const std::size_t slot = static_cast<std::size_t>(index);
          if (direct < 0) {
            ++negative[slot];
            if (!reported_negative[slot]) {
              reported_negative[slot] = true;
              std::cout << "TP2_BETA_ALIAS first_negative"
                        << " index=" << index
                        << " rank=" << rank
                        << " minus_power=" << minus_power
                        << " plus_power=" << plus_power
                        << " value=" << direct << '\n';
            }
          } else if (direct > 0) {
            ++positive[slot];
            if (!reported_positive[slot]) {
              reported_positive[slot] = true;
              std::cout << "TP2_BETA_ALIAS first_positive"
                        << " index=" << index
                        << " rank=" << rank
                        << " minus_power=" << minus_power
                        << " plus_power=" << plus_power
                        << " value=" << direct << '\n';
            }
          }
        }
        ++systems;
      }
    }
  }
  std::cout << "TP2_BETA_ALIAS maximum_rank=" << maximum_rank
            << " maximum_minus_power=" << maximum_minus_power
            << " maximum_plus_power=" << maximum_plus_power
            << " maximum_index=" << maximum_index
            << " systems=" << systems
            << " identity_checks=" << identity_checks;
  for (int index = 0; index <= maximum_index; ++index) {
    const std::size_t slot = static_cast<std::size_t>(index);
    std::cout << " negative_" << index << '=' << negative[slot]
              << " positive_" << index << '=' << positive[slot];
  }
  std::cout << " result=PASS\n";
  return 0;
}
