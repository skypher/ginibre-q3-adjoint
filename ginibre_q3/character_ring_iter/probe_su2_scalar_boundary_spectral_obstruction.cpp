#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

Integer value(const std::vector<Integer>& profile, const int index) {
  return index >= 0 && index < static_cast<int>(profile.size())
             ? profile[static_cast<std::size_t>(index)]
             : Integer(0);
}

Integer current(const std::vector<Integer>& profile,
                const int first, const int second) {
  Integer interval = 0;
  for (int label = std::abs(first - second);
       label <= first + second; ++label) {
    interval += value(profile, label);
  }
  return profile.front() * interval -
         value(profile, first) * value(profile, second);
}

void require(const bool condition, const char* message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

}  // namespace

int main() {
  try {
    const std::vector<Integer> profile{
        10485760, 20971520, 31457280, 34611200,
        16777216, 10240, 5};
    const int support = static_cast<int>(profile.size()) - 1;

    Integer minimum_shape =
        profile[1] * profile[1] - profile[0] * profile[2];
    for (std::size_t index = 1U; index + 1U < profile.size(); ++index) {
      const Integer margin =
          profile[index] * profile[index] -
          profile[index - 1U] * profile[index + 1U];
      minimum_shape = std::min(minimum_shape, margin);
    }

    Integer minimum_reflection =
        profile[0] * profile[5] - profile[1] * profile[6];
    for (int depth = 0; depth < support / 2; ++depth) {
      const Integer margin =
          profile[static_cast<std::size_t>(depth)] *
                  profile[static_cast<std::size_t>(
                      support - depth - 1)] -
          profile[static_cast<std::size_t>(depth + 1)] *
                  profile[static_cast<std::size_t>(
                      support - depth)];
      minimum_reflection = std::min(minimum_reflection, margin);
    }

    Integer central_value = 0;
    Integer quarter_value = 0;
    for (std::size_t index = 0; index < profile.size(); ++index) {
      central_value +=
          index % 2U == 0U ? profile[index] : -profile[index];
      quarter_value +=
          index % 4U <= 1U ? profile[index] : -profile[index];
    }

    Integer minimum_first_row = current(profile, 1, 2);
    for (int second = 2; second <= support; ++second) {
      minimum_first_row =
          std::min(minimum_first_row, current(profile, 1, second));
    }

    Integer minimum_diagonal = current(profile, 1, 1);
    for (int radius = 1; radius <= support; ++radius) {
      minimum_diagonal =
          std::min(minimum_diagonal, current(profile, radius, radius));
    }

    const Integer obstruction = current(profile, 2, 3);
    require(minimum_shape == 20971520, "shape margin mismatch");
    require(
        minimum_reflection == 107269324800,
        "reflection margin mismatch");
    require(central_value == 3127301, "central value mismatch");
    require(
        minimum_first_row == 107321753600,
        "first-row margin mismatch");
    require(
        minimum_diagonal == 725834792960,
        "diagonal margin mismatch");
    require(obstruction == -64424509440, "current mismatch");
    require(quarter_value == -17823749, "quarter value mismatch");

    std::cout
        << "SU2_SCALAR_BOUNDARY_SPECTRAL_OBSTRUCTION"
        << " minimum_shape=" << minimum_shape
        << " minimum_reflection=" << minimum_reflection
        << " central_value=" << central_value
        << " minimum_first_row=" << minimum_first_row
        << " minimum_diagonal=" << minimum_diagonal
        << " J_2_3=" << obstruction
        << " quarter_value=" << quarter_value
        << " result=PASS\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& exception) {
    std::cerr << exception.what() << '\n';
    return EXIT_FAILURE;
  }
}
