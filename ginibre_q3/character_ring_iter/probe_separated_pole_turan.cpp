#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

cpp_int binomial(int top, int bottom) {
  if (bottom < 0 || bottom > top) {
    return 0;
  }
  bottom = std::min(bottom, top - bottom);
  cpp_int value = 1;
  for (int factor = 1; factor <= bottom; ++factor) {
    value *= top - bottom + factor;
    value /= factor;
  }
  return value;
}

std::vector<cpp_int> coefficients(int separation,
                                  const std::vector<int>& poles,
                                  int numerator_power,
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
  cpp_int separation_power = 1;
  for (int degree = 0;
       degree <= std::min(numerator_power, maximum_degree); ++degree) {
    cpp_int value = binomial(numerator_power, degree) * separation_power;
    if ((numerator_power - degree) % 2 != 0) {
      value = -value;
    }
    numerator[static_cast<std::size_t>(degree)] = value;
    separation_power *= separation;
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
    result[static_cast<std::size_t>(degree)] = value;
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  const int maximum_separation = argc > 1 ? std::atoi(argv[1]) : 14;
  const int maximum_power = argc > 2 ? std::atoi(argv[2]) : 20;
  const int maximum_degree = argc > 3 ? std::atoi(argv[3]) : 50;
  if (maximum_separation < 3 || maximum_separation > 30 ||
      maximum_power < 0 || maximum_degree < 2) {
    std::cerr << "usage: probe_separated_pole_turan "
                 "[3<=maximum-separation<=30] [maximum-power>=0] "
                 "[maximum-degree>=2]\n";
    return 2;
  }

  std::size_t systems = 0;
  std::size_t determinants = 0;
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
      for (int power = 0; power <= maximum_power; ++power) {
        const std::vector<cpp_int> values = coefficients(
            separation, poles, power, maximum_degree + 2);
        ++systems;
        const int first_degree = std::max(0, power + 1 - rank);
        for (int degree = first_degree; degree <= maximum_degree; ++degree) {
          const cpp_int determinant =
              values[static_cast<std::size_t>(degree + 1)] *
                  values[static_cast<std::size_t>(degree + 1)] -
              values[static_cast<std::size_t>(degree)] *
                  values[static_cast<std::size_t>(degree + 2)];
          ++determinants;
          if (determinant < 0) {
            std::cout << "SEPARATED_POLE_TURAN result=COUNTEREXAMPLE"
                      << " separation=" << separation
                      << " numerator_power=" << power
                      << " degree=" << degree << " poles=";
            for (int pole : poles) {
              std::cout << pole << ',';
            }
            std::cout << " determinant=" << determinant << '\n';
            return 1;
          }
        }
      }
    }
  }
  std::cout << "SEPARATED_POLE_TURAN maximum_separation="
            << maximum_separation << " maximum_power=" << maximum_power
            << " maximum_degree=" << maximum_degree
            << " systems=" << systems << " determinants=" << determinants
            << " result=PASS\n";
  return 0;
}
