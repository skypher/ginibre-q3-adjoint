#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;
using Rational = boost::rational<cpp_int>;

namespace {

using Vector = std::vector<cpp_int>;

struct JacobiData {
  std::vector<Rational> alpha;
  std::vector<Rational> beta;
};

Rational rational(cpp_int numerator, cpp_int denominator = 1) {
  return Rational(std::move(numerator), std::move(denominator));
}

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

JacobiData base_jacobi(int rank) {
  const int modulus = 2 * rank + 1;
  JacobiData result;
  result.alpha.resize(static_cast<std::size_t>(rank));
  result.beta.assign(static_cast<std::size_t>(rank), rational(0));
  result.alpha[0] = rational(-1, modulus - 1);
  for (int index = 1; index < rank; ++index) {
    const int first = modulus - 2 * index + 1;
    const int second = modulus - 2 * index - 1;
    result.alpha[static_cast<std::size_t>(index)] =
        rational(-2, cpp_int(first) * second);
  }
  if (rank >= 2) {
    result.beta[1] = rational(cpp_int(modulus) * (modulus - 3),
                             2 * cpp_int(modulus - 1) * (modulus - 1));
  }
  for (int index = 2; index < rank; ++index) {
    const int first = modulus - 2 * index - 1;
    const int second = modulus - 2 * index + 3;
    const int denominator = modulus - 2 * index + 1;
    result.beta[static_cast<std::size_t>(index)] =
        rational(cpp_int(first) * second,
                 4 * cpp_int(denominator) * denominator);
  }
  return result;
}

Rational base_endpoint_ratio(int rank, int index, int endpoint) {
  const int modulus = 2 * rank + 1;
  if (endpoint == 1) {
    if (index == 0) {
      return rational(modulus, modulus - 1);
    }
    return rational(modulus - 2 * index + 1,
                    2 * (modulus - 2 * index - 1));
  }
  if (endpoint == -1) {
    if (index == 0) {
      return rational(-(modulus - 2), modulus - 1);
    }
    return rational(-cpp_int(modulus - 2 * index - 2) *
                        (modulus - 2 * index + 1),
                    2 * cpp_int(modulus - 2 * index - 1) *
                        (modulus - 2 * index));
  }
  throw std::invalid_argument("endpoint must be plus or minus one");
}

void verify_base_endpoint_ratios(const JacobiData& data, int rank) {
  for (const int endpoint : {-1, 1}) {
    Rational previous = rational(1);
    for (int index = 0; index < rank; ++index) {
      const Rational ratio = base_endpoint_ratio(rank, index, endpoint);
      const Rational recurrence =
          rational(endpoint) - data.alpha[static_cast<std::size_t>(index)] -
          (index == 0
               ? rational(0)
               : data.beta[static_cast<std::size_t>(index)] / previous);
      if (ratio != recurrence) {
        throw std::runtime_error("base endpoint-ratio recurrence mismatch");
      }
      previous = ratio;
    }
  }
}

JacobiData christoffel_step(const JacobiData& data, int endpoint) {
  const std::size_t size = data.alpha.size();
  JacobiData result;
  result.alpha.resize(size);
  result.beta.assign(size, rational(0));
  std::vector<Rational> diagonal(size);
  diagonal[0] = data.alpha[0] - rational(endpoint);
  if (diagonal[0] == rational(0)) {
    throw std::runtime_error("singular endpoint Darboux pivot");
  }
  for (std::size_t index = 1U; index < size; ++index) {
    diagonal[index] = data.alpha[index] - rational(endpoint) -
                      data.beta[index] / diagonal[index - 1U];
    if (diagonal[index] == rational(0)) {
      throw std::runtime_error("singular endpoint Darboux pivot");
    }
    result.beta[index] =
        data.beta[index] * diagonal[index] / diagonal[index - 1U];
  }
  for (std::size_t index = 0U; index + 1U < size; ++index) {
    result.alpha[index] = rational(endpoint) + diagonal[index] +
                          data.beta[index + 1U] / diagonal[index];
  }
  result.alpha.back() = rational(endpoint) + diagonal.back();
  return result;
}

JacobiData endpoint_jacobi(int rank, int minus_power, int plus_power) {
  JacobiData result = base_jacobi(rank);
  verify_base_endpoint_ratios(result, rank);
  for (int exponent = 0; exponent < minus_power; ++exponent) {
    result = christoffel_step(result, 1);
  }
  for (int exponent = 0; exponent < plus_power; ++exponent) {
    result = christoffel_step(result, -1);
  }
  return result;
}

Rational recurrence_target(const JacobiData& data) {
  if (data.alpha.size() < 3U) {
    throw std::invalid_argument("rank must be at least three");
  }
  const Rational alpha_zero = data.alpha[0];
  const Rational alpha_one = data.alpha[1];
  const Rational beta_one = data.beta[1];
  const Rational beta_two = data.beta[2];
  const Rational u_zero = 2 * (alpha_zero * alpha_zero + beta_one) - 1;
  const Rational u_one = 2 * (alpha_zero + alpha_one);
  const Rational v_zero =
      4 * (alpha_zero * alpha_zero * alpha_zero +
           beta_one * (2 * alpha_zero + alpha_one)) -
      3 * alpha_zero;
  const Rational v_one =
      4 * (alpha_zero * alpha_zero + beta_one + beta_two +
           alpha_one * (alpha_zero + alpha_one)) -
      3;
  return v_zero * u_one - v_one * u_zero;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc == 5 && std::string(argv[1]) == "--sample") {
    const int rank = std::atoi(argv[2]);
    const int minus_power = std::atoi(argv[3]);
    const int plus_power = std::atoi(argv[4]);
    if (rank < 3 || minus_power < 0 || plus_power < 0) {
      std::cerr << "usage: verify_tp2_finite_darboux --sample "
                   "[rank>=3] [minus-power>=0] [plus-power>=0]\n";
      return 2;
    }
    const JacobiData data = endpoint_jacobi(rank, minus_power, plus_power);
    std::cout << "TP2_FINITE_DARBOUX_SAMPLE"
              << " rank=" << rank
              << " minus_power=" << minus_power
              << " plus_power=" << plus_power
              << " alpha_0=" << data.alpha[0]
              << " alpha_1=" << data.alpha[1]
              << " beta_1=" << data.beta[1]
              << " beta_2=" << data.beta[2]
              << " recurrence_target=" << recurrence_target(data) << '\n';
    return 0;
  }
  const int maximum_rank = argc > 1 ? std::atoi(argv[1]) : 20;
  const int maximum_minus_power = argc > 2 ? std::atoi(argv[2]) : 30;
  const int maximum_plus_power = argc > 3 ? std::atoi(argv[3]) : 30;
  if (maximum_rank < 3 || maximum_minus_power < 2 ||
      maximum_plus_power < 2) {
    std::cerr << "usage: verify_tp2_finite_darboux "
                 "[maximum-rank>=3] [maximum-minus-power>=2] "
                 "[maximum-plus-power>=2]\n";
    return 2;
  }

  std::size_t systems = 0U;
  for (int rank = 3; rank <= maximum_rank; ++rank) {
    const int modulus = 2 * rank + 1;
    for (int minus_power = 2; minus_power <= maximum_minus_power;
         ++minus_power) {
      for (int plus_power = 2; plus_power <= maximum_plus_power;
           ++plus_power) {
        const JacobiData data = endpoint_jacobi(
            rank, minus_power, plus_power
        );
        for (std::size_t index = 1U; index < data.beta.size(); ++index) {
          if (data.beta[index] <= rational(0)) {
            std::cout << "TP2_FINITE_DARBOUX result=BETA_FAIL"
                      << " rank=" << rank
                      << " minus_power=" << minus_power
                      << " plus_power=" << plus_power
                      << " index=" << index << '\n';
            return 1;
          }
        }
        const Rational target = recurrence_target(data);
        const Vector left = endpoint_vector(
            modulus, minus_power + 1, plus_power
        );
        const Vector right = endpoint_vector(
            modulus, minus_power, plus_power + 1
        );
        const cpp_int cross = pascal_cross(left, right);
        ++systems;
        if ((cross < 0) != (target > rational(0)) ||
            (cross == 0) != (target == rational(0))) {
          std::cout << "TP2_FINITE_DARBOUX result=SIGN_MISMATCH"
                    << " rank=" << rank
                    << " minus_power=" << minus_power
                    << " plus_power=" << plus_power
                    << " cross=" << cross
                    << " recurrence_target=" << target << '\n';
          return 1;
        }
        if (target > rational(0)) {
          std::cout << "TP2_FINITE_DARBOUX result=TARGET_FAIL"
                    << " rank=" << rank
                    << " minus_power=" << minus_power
                    << " plus_power=" << plus_power
                    << " recurrence_target=" << target << '\n';
          return 1;
        }
      }
    }
  }
  std::cout << "TP2_FINITE_DARBOUX maximum_rank=" << maximum_rank
            << " maximum_minus_power=" << maximum_minus_power
            << " maximum_plus_power=" << maximum_plus_power
            << " systems=" << systems << " result=PASS\n";
  return 0;
}
