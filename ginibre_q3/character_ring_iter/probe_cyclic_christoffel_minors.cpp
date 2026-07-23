#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Vec = std::vector<cpp_int>;

cpp_int cyclic_weighted_moment(int modulus, int minus_pairs, int half_power) {
  const int power = 2 * half_power;
  const int degree = 2 * minus_pairs + power;
  Vec coefficients(static_cast<std::size_t>(2 * degree + 1));
  coefficients[static_cast<std::size_t>(degree)] = 1;
  int radius = 0;
  for (int factor = 0; factor < minus_pairs; ++factor) {
    Vec next(coefficients.size());
    for (int exponent = -radius; exponent <= radius; ++exponent) {
      const cpp_int& value =
          coefficients[static_cast<std::size_t>(degree + exponent)];
      next[static_cast<std::size_t>(degree + exponent)] += 2 * value;
      next[static_cast<std::size_t>(degree + exponent - 2)] -= value;
      next[static_cast<std::size_t>(degree + exponent + 2)] -= value;
    }
    coefficients = std::move(next);
    radius += 2;
  }
  for (int factor = 0; factor < power; ++factor) {
    Vec next(coefficients.size());
    for (int exponent = -radius; exponent <= radius; ++exponent) {
      const cpp_int& value =
          coefficients[static_cast<std::size_t>(degree + exponent)];
      next[static_cast<std::size_t>(degree + exponent - 1)] += value;
      next[static_cast<std::size_t>(degree + exponent + 1)] += value;
    }
    coefficients = std::move(next);
    ++radius;
  }
  cpp_int result = 0;
  for (int exponent = -degree; exponent <= degree; ++exponent) {
    if (exponent % modulus == 0) {
      result += coefficients[static_cast<std::size_t>(degree + exponent)];
    }
  }
  return result;
}

std::vector<Vec> chebyshev_polynomials(int maximum_index) {
  std::vector<Vec> result(static_cast<std::size_t>(maximum_index + 1));
  result[0] = Vec{cpp_int(2)};
  if (maximum_index == 0) {
    return result;
  }
  result[1] = Vec{cpp_int(-2), cpp_int(1)};
  for (int index = 1; index < maximum_index; ++index) {
    Vec next(result[static_cast<std::size_t>(index)].size() + 1U);
    const Vec& current = result[static_cast<std::size_t>(index)];
    const Vec& previous = result[static_cast<std::size_t>(index - 1)];
    for (std::size_t degree = 0; degree < current.size(); ++degree) {
      next[degree] -= 2 * current[degree];
      next[degree + 1U] += current[degree];
    }
    for (std::size_t degree = 0; degree < previous.size(); ++degree) {
      next[degree] -= previous[degree];
    }
    result[static_cast<std::size_t>(index + 1)] = std::move(next);
  }
  return result;
}

cpp_int evaluate_functional(const Vec& polynomial, const Vec& moments,
                            int shift) {
  cpp_int result = 0;
  for (std::size_t degree = 0; degree < polynomial.size(); ++degree) {
    result += polynomial[degree] *
              moments[degree + static_cast<std::size_t>(shift)];
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  const int maximum_rank = argc > 1 ? std::atoi(argv[1]) : 20;
  const int maximum_minus_pairs = argc > 2 ? std::atoi(argv[2]) : 30;
  const int maximum_half_power = argc > 3 ? std::atoi(argv[3]) : 30;
  const int maximum_index = argc > 4 ? std::atoi(argv[4]) : 8;
  if (maximum_rank < 2 || maximum_minus_pairs < 0 ||
      maximum_half_power < 0 || maximum_index < 1) {
    std::cerr << "usage: probe_cyclic_christoffel_minors "
                 "[maximum-rank>=2] [maximum-minus-pairs>=0] "
                 "[maximum-half-power>=0] [maximum-index>=1]\n";
    return 2;
  }

  const std::vector<Vec> polynomials =
      chebyshev_polynomials(maximum_index + 1);
  std::size_t checks = 0;
  std::size_t negative = 0;
  std::size_t target_checks = 0;
  std::size_t alternating_failures = 0;
  std::vector<bool> reported_negative(
      static_cast<std::size_t>(maximum_index + 1));
  for (int rank = 2; rank <= maximum_rank; ++rank) {
    const int modulus = 2 * rank + 1;
    for (int minus_pairs = 0; minus_pairs <= maximum_minus_pairs;
         ++minus_pairs) {
      for (int half_power = 0; half_power <= maximum_half_power;
           ++half_power) {
        Vec moments(static_cast<std::size_t>(maximum_index + 3));
        for (int degree = 0; degree <= maximum_index + 2; ++degree) {
          moments[static_cast<std::size_t>(degree)] =
              cyclic_weighted_moment(modulus, minus_pairs,
                                     half_power + degree);
        }
        for (int index = 1; index <= maximum_index; ++index) {
          const Vec& left = polynomials[static_cast<std::size_t>(index)];
          const Vec& right =
              polynomials[static_cast<std::size_t>(index + 1)];
          const cpp_int determinant =
              evaluate_functional(right, moments, 1) *
                  evaluate_functional(left, moments, 0) -
              evaluate_functional(right, moments, 0) *
                  evaluate_functional(left, moments, 1);
          ++checks;
          if (index == 2 && determinant < 0) {
            std::cout << "CYCLIC_CHRISTOFFEL result=INDEX2_FAIL"
                      << " rank=" << rank
                      << " minus_pairs=" << minus_pairs
                      << " half_power=" << half_power
                      << " determinant=" << determinant << '\n';
            return 1;
          }
          if (index == 2 && half_power >= 2 &&
              minus_pairs >= half_power + 2) {
            ++target_checks;
            if (determinant < 0) {
              std::cout << "CYCLIC_CHRISTOFFEL result=TARGET_FAIL"
                        << " rank=" << rank
                        << " minus_pairs=" << minus_pairs
                        << " half_power=" << half_power
                        << " index=" << index
                        << " determinant=" << determinant << '\n';
              return 1;
            }
          }
          if (determinant < 0) {
            ++negative;
            if (!reported_negative[static_cast<std::size_t>(index)]) {
              reported_negative[static_cast<std::size_t>(index)] = true;
              std::cout << "CYCLIC_CHRISTOFFEL first_negative"
                        << " rank=" << rank
                        << " minus_pairs=" << minus_pairs
                        << " half_power=" << half_power
                        << " difference=" << minus_pairs - half_power
                        << " index=" << index
                        << " determinant=" << determinant << '\n';
            }
          }
        }
        if (half_power >= 2 && minus_pairs >= half_power + 2) {
          for (int index = 1; index <= maximum_index + 1; ++index) {
            cpp_int value = evaluate_functional(
                polynomials[static_cast<std::size_t>(index)], moments, 0);
            if (index % 2 != 0) {
              value = -value;
            }
            if (value < 0) {
              ++alternating_failures;
            }
          }
        }
      }
    }
  }
  std::cout << "CYCLIC_CHRISTOFFEL maximum_rank=" << maximum_rank
            << " maximum_minus_pairs=" << maximum_minus_pairs
            << " maximum_half_power=" << maximum_half_power
            << " maximum_index=" << maximum_index << " checks=" << checks
            << " negative=" << negative << " target_checks=" << target_checks
            << " alternating_failures=" << alternating_failures
            << " result=PASS\n";
  return 0;
}
