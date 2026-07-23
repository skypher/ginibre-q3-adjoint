#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

cpp_int binomial(int top, int bottom) {
  if (bottom < 0 || bottom > top) {
    return 0;
  }
  bottom = std::min(bottom, top - bottom);
  cpp_int result = 1;
  for (int factor = 1; factor <= bottom; ++factor) {
    result *= top - bottom + factor;
    result /= factor;
  }
  return result;
}

cpp_int power(int base, int exponent) {
  cpp_int result = 1;
  cpp_int factor = base;
  for (int remaining = exponent; remaining > 0; remaining /= 2) {
    if (remaining % 2 != 0) {
      result *= factor;
    }
    if (remaining > 1) {
      factor *= factor;
    }
  }
  return result;
}

cpp_int falling_factorial(int top, int length) {
  if (length < 0 || length > top) {
    return 0;
  }
  cpp_int result = 1;
  for (int offset = 0; offset < length; ++offset) {
    result *= top - offset;
  }
  return result;
}

cpp_int derivative_value(int first_power, int second_power, int separation,
                         int argument, int derivative) {
  cpp_int result = 0;
  for (int chosen = 0; chosen <= second_power; ++chosen) {
    const int degree = first_power + chosen;
    if (degree < derivative) {
      continue;
    }
    cpp_int term = binomial(second_power, chosen) *
                   power(separation, second_power - chosen) *
                   falling_factorial(degree, derivative) *
                   power(argument, degree - derivative);
    if (chosen % 2 != 0) {
      term = -term;
    }
    result += term;
  }
  return result;
}

cpp_int determinant(std::vector<std::vector<cpp_int>> matrix) {
  const int size = static_cast<int>(matrix.size());
  if (size == 0) {
    return 1;
  }
  cpp_int previous = 1;
  int sign = 1;
  for (int pivot_index = 0; pivot_index + 1 < size; ++pivot_index) {
    int pivot_row = pivot_index;
    while (pivot_row < size &&
           matrix[static_cast<std::size_t>(pivot_row)]
                 [static_cast<std::size_t>(pivot_index)] == 0) {
      ++pivot_row;
    }
    if (pivot_row == size) {
      return 0;
    }
    if (pivot_row != pivot_index) {
      std::swap(matrix[static_cast<std::size_t>(pivot_row)],
                matrix[static_cast<std::size_t>(pivot_index)]);
      sign = -sign;
    }
    const cpp_int pivot =
        matrix[static_cast<std::size_t>(pivot_index)]
              [static_cast<std::size_t>(pivot_index)];
    for (int row = pivot_index + 1; row < size; ++row) {
      for (int column = pivot_index + 1; column < size; ++column) {
        cpp_int value =
            matrix[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(column)] *
                pivot -
            matrix[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(pivot_index)] *
                matrix[static_cast<std::size_t>(pivot_index)]
                      [static_cast<std::size_t>(column)];
        if (pivot_index > 0) {
          value /= previous;
        }
        matrix[static_cast<std::size_t>(row)]
              [static_cast<std::size_t>(column)] = std::move(value);
      }
    }
    previous = pivot;
  }
  cpp_int result = matrix.back().back();
  if (sign < 0) {
    result = -result;
  }
  return result;
}

std::vector<cpp_int> coefficients(int separation,
                                  const std::vector<int>& poles, int power_value,
                                  int maximum_degree) {
  std::vector<cpp_int> denominator(1, cpp_int(1));
  for (int pole : poles) {
    denominator.push_back(0);
    for (std::size_t degree = denominator.size() - 1U; degree > 0U;
         --degree) {
      denominator[degree] -= pole * denominator[degree - 1U];
    }
  }
  std::vector<cpp_int> numerator(
      static_cast<std::size_t>(maximum_degree + 1));
  for (int degree = 0;
       degree <= std::min(power_value, maximum_degree); ++degree) {
    cpp_int value = binomial(power_value, degree) * power(separation, degree);
    if ((power_value - degree) % 2 != 0) {
      value = -value;
    }
    numerator[static_cast<std::size_t>(degree)] = value;
  }
  std::vector<cpp_int> result(
      static_cast<std::size_t>(maximum_degree + 1));
  for (int degree = 0; degree <= maximum_degree; ++degree) {
    cpp_int value = numerator[static_cast<std::size_t>(degree)];
    const int maximum_offset = std::min(
        degree, static_cast<int>(denominator.size()) - 1);
    for (int offset = 1; offset <= maximum_offset; ++offset) {
      value -= denominator[static_cast<std::size_t>(offset)] *
               result[static_cast<std::size_t>(degree - offset)];
    }
    result[static_cast<std::size_t>(degree)] = std::move(value);
  }
  return result;
}

cpp_int evaluation_determinant(int separation, const std::vector<int>& poles,
                               int power_value, int q) {
  const int rank = static_cast<int>(poles.size());
  std::vector<std::vector<cpp_int>> matrix(
      static_cast<std::size_t>(rank),
      std::vector<cpp_int>(static_cast<std::size_t>(rank)));
  for (int column = 0; column < rank; ++column) {
    const int pole = poles[static_cast<std::size_t>(column)];
    const cpp_int f =
        power(pole, q) * power(separation - pole, power_value);
    matrix[0][static_cast<std::size_t>(column)] = f;
    matrix[1][static_cast<std::size_t>(column)] = pole * f;
    for (int row = 2; row < rank; ++row) {
      matrix[static_cast<std::size_t>(row)]
            [static_cast<std::size_t>(column)] = power(pole, row - 2);
    }
  }
  return determinant(std::move(matrix));
}

cpp_int vandermonde(const std::vector<int>& poles) {
  cpp_int result = 1;
  for (std::size_t first = 0; first < poles.size(); ++first) {
    for (std::size_t second = first + 1U; second < poles.size(); ++second) {
      result *= poles[second] - poles[first];
    }
  }
  return result;
}

bool check_wronskians(int separation, int rank, int power_value, int q,
                      std::size_t& checks) {
  for (int argument = 1; argument < separation; ++argument) {
    for (int s = 0; s <= rank - 2; ++s) {
      const int size = s + 2;
      std::vector<std::vector<cpp_int>> matrix(
          static_cast<std::size_t>(size),
          std::vector<cpp_int>(static_cast<std::size_t>(size)));
      for (int derivative = 0; derivative < size; ++derivative) {
        matrix[static_cast<std::size_t>(derivative)][0] = derivative_value(
            q, power_value, separation, argument, derivative);
        matrix[static_cast<std::size_t>(derivative)][1] = derivative_value(
            q + 1, power_value, separation, argument, derivative);
        for (int monomial = 0; monomial < s; ++monomial) {
          matrix[static_cast<std::size_t>(derivative)]
                [static_cast<std::size_t>(monomial + 2)] =
              falling_factorial(monomial, derivative) *
              power(argument, monomial - derivative);
        }
      }
      cpp_int factor = 1;
      for (int index = 0; index < s; ++index) {
        factor *= falling_factorial(index, index);
      }
      const cpp_int middle = derivative_value(
          q, power_value, separation, argument, s);
      cpp_int bracket = (s + 1) * middle * middle;
      if (s > 0) {
        bracket -=
            s * derivative_value(q, power_value, separation, argument, s - 1) *
            derivative_value(q, power_value, separation, argument, s + 1);
      }
      const cpp_int direct = determinant(std::move(matrix));
      ++checks;
      if (direct != factor * bracket || direct <= 0) {
        std::cout << "SEPARATED_POLE_ECT result=WRONSKIAN_FAIL"
                  << " separation=" << separation << " rank=" << rank
                  << " numerator_power=" << power_value << " q=" << q
                  << " argument=" << argument << " s=" << s
                  << " direct=" << direct << " formula=" << factor * bracket
                  << '\n';
        return false;
      }
    }
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  const int maximum_separation = argc > 1 ? std::atoi(argv[1]) : 8;
  const int maximum_power = argc > 2 ? std::atoi(argv[2]) : 8;
  const int maximum_degree = argc > 3 ? std::atoi(argv[3]) : 12;
  if (maximum_separation < 3 || maximum_separation > 20 ||
      maximum_power < 0 || maximum_degree < 2) {
    std::cerr << "usage: verify_separated_pole_turan_ect "
                 "[3<=maximum-separation<=20] [maximum-power>=0] "
                 "[maximum-degree>=2]\n";
    return 2;
  }

  std::size_t coefficient_checks = 0;
  std::size_t determinant_checks = 0;
  std::size_t wronskian_checks = 0;
  for (int separation = 3; separation <= maximum_separation; ++separation) {
    const int available = separation - 1;
    const unsigned long long subset_count = 1ULL << available;
    for (unsigned long long subset = 1; subset < subset_count; ++subset) {
      std::vector<int> poles;
      for (int pole = 1; pole <= available; ++pole) {
        if ((subset & (1ULL << (pole - 1))) != 0ULL) {
          poles.push_back(pole);
        }
      }
      const int rank = static_cast<int>(poles.size());
      for (int power_value = 0; power_value <= maximum_power; ++power_value) {
        const std::vector<cpp_int> values = coefficients(
            separation, poles, power_value, maximum_degree + 2);
        const int first_degree = std::max(0, power_value + 1 - rank);
        for (int degree = first_degree; degree <= maximum_degree; ++degree) {
          const cpp_int turan =
              values[static_cast<std::size_t>(degree + 1)] *
                  values[static_cast<std::size_t>(degree + 1)] -
              values[static_cast<std::size_t>(degree)] *
                  values[static_cast<std::size_t>(degree + 2)];
          ++coefficient_checks;
          if ((rank == 1 && turan != 0) || (rank >= 2 && turan <= 0)) {
            std::cout << "SEPARATED_POLE_ECT result=TURAN_FAIL"
                      << " separation=" << separation << " rank=" << rank
                      << " numerator_power=" << power_value
                      << " degree=" << degree << " turan=" << turan << '\n';
            return 1;
          }
          if (rank >= 2) {
            const int q = degree + rank - 1 - power_value;
            const cpp_int direct = evaluation_determinant(
                separation, poles, power_value, q);
            const cpp_int scaled = turan * vandermonde(poles);
            ++determinant_checks;
            if (direct != scaled) {
              std::cout << "SEPARATED_POLE_ECT result=DETERMINANT_FAIL"
                        << " separation=" << separation << " rank=" << rank
                        << " numerator_power=" << power_value
                        << " degree=" << degree << " direct=" << direct
                        << " scaled_turan=" << scaled << '\n';
              return 1;
            }
            if (!check_wronskians(separation, rank, power_value, q,
                                  wronskian_checks)) {
              return 1;
            }
          }
        }
      }
    }
  }
  std::cout << "SEPARATED_POLE_ECT maximum_separation=" << maximum_separation
            << " maximum_power=" << maximum_power
            << " maximum_degree=" << maximum_degree
            << " coefficient_checks=" << coefficient_checks
            << " determinant_checks=" << determinant_checks
            << " wronskian_checks=" << wronskian_checks << " result=PASS\n";
  return 0;
}
