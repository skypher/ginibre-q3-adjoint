#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
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

class DeterministicGenerator {
 public:
  std::uint64_t next() {
    state_ = state_ * 6364136223846793005ULL + 1442695040888963407ULL;
    return state_;
  }

  int bounded(int upper_bound) {
    return static_cast<int>(next() % static_cast<std::uint64_t>(upper_bound));
  }

 private:
  std::uint64_t state_ = 0x9e3779b97f4a7c15ULL;
};

bool check_system(int separation, const std::vector<int>& poles, int power,
                  int maximum_degree, std::size_t& determinants) {
  const std::vector<cpp_int> values =
      coefficients(separation, poles, power, maximum_degree + 2);
  const int rank = static_cast<int>(poles.size());
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
                << " numerator_power=" << power << " degree=" << degree
                << " poles=";
      for (int pole : poles) {
        std::cout << pole << ',';
      }
      std::cout << " determinant=" << determinant << '\n';
      return false;
    }
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  const int maximum_separation = argc > 1 ? std::atoi(argv[1]) : 14;
  const int maximum_power = argc > 2 ? std::atoi(argv[2]) : 20;
  const int maximum_degree = argc > 3 ? std::atoi(argv[3]) : 50;
  const int random_systems = argc > 4 ? std::atoi(argv[4]) : 0;
  const int random_maximum_separation = argc > 5 ? std::atoi(argv[5]) : 1000;
  if (maximum_separation < 3 || maximum_separation > 30 ||
      maximum_power < 0 || maximum_degree < 2 || random_systems < 0 ||
      random_maximum_separation < 3) {
    std::cerr << "usage: probe_separated_pole_turan "
                 "[3<=maximum-separation<=30] [maximum-power>=0] "
                 "[maximum-degree>=2] [random-systems>=0] "
                 "[random-maximum-separation>=3]\n";
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
      for (int power = 0; power <= maximum_power; ++power) {
        ++systems;
        if (!check_system(separation, poles, power, maximum_degree,
                          determinants)) {
          return 1;
        }
      }
    }
  }

  DeterministicGenerator generator;
  for (int trial = 0; trial < random_systems; ++trial) {
    const int separation =
        3 + generator.bounded(random_maximum_separation - 2);
    const int available = separation - 1;
    const int maximum_rank = std::min(available, 30);
    const int rank = 1 + generator.bounded(maximum_rank);
    std::vector<bool> selected(static_cast<std::size_t>(available + 1));
    std::vector<int> poles;
    poles.reserve(static_cast<std::size_t>(rank));
    while (static_cast<int>(poles.size()) < rank) {
      const int pole = 1 + generator.bounded(available);
      if (!selected[static_cast<std::size_t>(pole)]) {
        selected[static_cast<std::size_t>(pole)] = true;
        poles.push_back(pole);
      }
    }
    std::sort(poles.begin(), poles.end());
    const int power = generator.bounded(maximum_power + 1);
    ++systems;
    if (!check_system(separation, poles, power, maximum_degree,
                      determinants)) {
      return 1;
    }
  }
  std::cout << "SEPARATED_POLE_TURAN maximum_separation="
            << maximum_separation << " maximum_power=" << maximum_power
            << " maximum_degree=" << maximum_degree
            << " random_systems=" << random_systems
            << " random_maximum_separation=" << random_maximum_separation
            << " systems=" << systems << " determinants=" << determinants
            << " result=PASS\n";
  return 0;
}
