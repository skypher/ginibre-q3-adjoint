#include <boost/multiprecision/cpp_int.hpp>

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Vector = std::vector<cpp_int>;

int residue(int value, int modulus) {
  value %= modulus;
  return value < 0 ? value + modulus : value;
}

Vector update(const Vector& values, int neighbour_sign) {
  const int modulus = static_cast<int>(values.size());
  Vector result(values.size());
  for (int index = 0; index < modulus; ++index) {
    result[static_cast<std::size_t>(index)] =
        2 * values[static_cast<std::size_t>(index)] +
        neighbour_sign *
            (values[static_cast<std::size_t>(residue(index - 1, modulus))] +
             values[static_cast<std::size_t>(residue(index + 1, modulus))]);
  }
  return result;
}

Vector endpoint_vector(int modulus, int minus_power, int plus_power) {
  Vector result(static_cast<std::size_t>(modulus));
  result[0] = 1;
  for (int exponent = 0; exponent < minus_power; ++exponent) {
    result = update(result, -1);
  }
  for (int exponent = 0; exponent < plus_power; ++exponent) {
    result = update(result, 1);
  }
  return result;
}

Vector adjacency_update(const Vector& values) {
  const int modulus = static_cast<int>(values.size());
  Vector result(values.size());
  for (int index = 0; index < modulus; ++index) {
    result[static_cast<std::size_t>(index)] =
        values[static_cast<std::size_t>(residue(index - 1, modulus))] +
        values[static_cast<std::size_t>(residue(index + 1, modulus))];
  }
  return result;
}

cpp_int pascal_cross(const Vector& left, const Vector& right) {
  const int modulus = static_cast<int>(left.size());
  const auto at = [&](const Vector& values, int index) -> const cpp_int& {
    return values[static_cast<std::size_t>(residue(index, modulus))];
  };
  return at(left, 2) * at(right, 4) +
         at(left, 4) * at(right, 2) +
         2 * at(left, 2) * at(right, 2) -
         2 * at(left, 3) * at(right, 3) -
         at(left, 1) * at(right, 3) -
         at(left, 3) * at(right, 1);
}

std::vector<cpp_int> Jacobi_moments(int modulus, int minus_power,
                                    int plus_power) {
  Vector values = endpoint_vector(modulus, minus_power, plus_power);
  std::vector<cpp_int> result(6);
  for (int degree = 0; degree <= 5; ++degree) {
    result[static_cast<std::size_t>(degree)] =
        modulus * values[0];
    values = adjacency_update(values);
  }
  return result;
}

void add_scaled(cpp_int& result, const cpp_int& left, const cpp_int& right,
                int coefficient, unsigned int scale) {
  cpp_int term = left * right;
  term <<= scale;
  result += coefficient * term;
}

cpp_int mixed_Jacobi_moment_form(const std::vector<cpp_int>& moments) {
  cpp_int result = 0;
  const auto add = [&](int first, int second, int coefficient) {
    const unsigned int scale = static_cast<unsigned int>(8 - first - second);
    add_scaled(result, moments[static_cast<std::size_t>(first)],
               moments[static_cast<std::size_t>(second)], coefficient, scale);
  };
  add(0, 2, 12);
  add(1, 1, -12);
  add(0, 4, -16);
  add(1, 3, 28);
  add(1, 5, 16);
  add(2, 2, -12);
  add(2, 4, -8);
  add(3, 3, -8);
  add(3, 5, -32);
  add(4, 4, 32);
  return result;
}

cpp_int Chebyshev_CD_form(const std::vector<cpp_int>& moments) {
  const cpp_int first_row_first =
      4 * moments[5] + 4 * moments[3] - 48 * moments[1];
  const cpp_int first_row_second =
      2 * moments[4] + 4 * moments[2] - 16 * moments[0];
  const cpp_int second_row_first = 4 * moments[4] - 12 * moments[2];
  const cpp_int second_row_second = 2 * moments[3] - 4 * moments[1];
  return -4 * (first_row_first * second_row_second -
               first_row_second * second_row_first);
}

cpp_int scaled_cross(int modulus, const cpp_int& cross) {
  cpp_int result = cpp_int(64) * modulus * modulus;
  return result * cross;
}

}  // namespace

int main(int argc, char** argv) {
  const int maximum_rank = argc > 1 ? std::atoi(argv[1]) : 20;
  const int maximum_minus_power = argc > 2 ? std::atoi(argv[2]) : 30;
  const int maximum_plus_power = argc > 3 ? std::atoi(argv[3]) : 30;
  if (maximum_rank < 1 || maximum_minus_power < 2 ||
      maximum_plus_power < 2) {
    std::cerr << "usage: verify_tp2_mixed_jacobi_moment "
                 "[maximum-rank>=1] [maximum-minus-power>=2] "
                 "[maximum-plus-power>=2]\n";
    return 2;
  }

  std::size_t systems = 0;
  for (int rank = 1; rank <= maximum_rank; ++rank) {
    const int modulus = 2 * rank + 1;
    for (int minus_power = 2; minus_power <= maximum_minus_power;
         ++minus_power) {
      for (int plus_power = 2; plus_power <= maximum_plus_power;
           ++plus_power) {
        const Vector left = endpoint_vector(modulus, minus_power + 1,
                                            plus_power);
        const Vector right = endpoint_vector(modulus, minus_power,
                                             plus_power + 1);
        const cpp_int cross = pascal_cross(left, right);
        const std::vector<cpp_int> moments =
            Jacobi_moments(modulus, minus_power, plus_power);
        const cpp_int moment_form = mixed_Jacobi_moment_form(moments);
        const cpp_int Chebyshev_form = Chebyshev_CD_form(moments);
        ++systems;
        if (moment_form != Chebyshev_form) {
          std::cout << "TP2_MIXED_JACOBI_MOMENT result=CHEBYSHEV_FAIL"
                    << " rank=" << rank
                    << " minus_power=" << minus_power
                    << " plus_power=" << plus_power
                    << " moment_form=" << moment_form
                    << " Chebyshev_form=" << Chebyshev_form << '\n';
          return 1;
        }
        if (moment_form != scaled_cross(modulus, cross)) {
          std::cout << "TP2_MIXED_JACOBI_MOMENT result=IDENTITY_FAIL"
                    << " rank=" << rank
                    << " minus_power=" << minus_power
                    << " plus_power=" << plus_power
                    << " moment_form=" << moment_form
                    << " scaled_cross=" << scaled_cross(modulus, cross)
                    << '\n';
          return 1;
        }
        if (rank <= 2 && cross != 0) {
          std::cout << "TP2_MIXED_JACOBI_MOMENT result=BOUNDARY_CYCLE_FAIL"
                    << " rank=" << rank
                    << " minus_power=" << minus_power
                    << " plus_power=" << plus_power
                    << " cross=" << cross << '\n';
          return 1;
        }
        if (cross < 0) {
          std::cout << "TP2_MIXED_JACOBI_MOMENT result=TARGET_FAIL"
                    << " rank=" << rank
                    << " minus_power=" << minus_power
                    << " plus_power=" << plus_power
                    << " cross=" << cross << '\n';
          return 1;
        }
      }
    }
  }
  std::cout << "TP2_MIXED_JACOBI_MOMENT maximum_rank=" << maximum_rank
            << " maximum_minus_power=" << maximum_minus_power
            << " maximum_plus_power=" << maximum_plus_power
            << " systems=" << systems << " result=PASS\n";
}
