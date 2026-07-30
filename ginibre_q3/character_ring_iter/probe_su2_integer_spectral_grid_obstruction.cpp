#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;

template <typename Scalar>
Scalar density(const std::array<Integer, 7>& profile,
               const Scalar& point) {
  Scalar previous = Scalar(1);
  Scalar present = point;
  Scalar result =
      Scalar(profile[0]) * previous + Scalar(profile[1]) * present;
  for (std::size_t index = 2U; index < profile.size(); ++index) {
    const Scalar next = (point - Scalar(1)) * present - previous;
    result += Scalar(profile[index]) * next;
    previous = present;
    present = next;
  }
  return result;
}

void require(const bool condition, const char* message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

}  // namespace

int main() {
  try {
    const std::array<Integer, 7> profile{
        8, 36, 40, 40, 36, 8, 1};

    Integer minimum_shape =
        profile[1] * profile[1] - profile[0] * profile[2];
    for (std::size_t index = 1U; index + 1U < profile.size(); ++index) {
      const Integer margin =
          profile[index] * profile[index] -
          profile[index - 1U] * profile[index + 1U];
      minimum_shape = std::min(minimum_shape, margin);
    }

    const Integer current =
        profile[0] *
            (profile[1] + profile[2] + profile[3] +
             profile[4] + profile[5]) -
        profile[2] * profile[3];

    const std::array<Integer, 5> expected_grid{1, 1, 7, 1, 1021};
    for (int point = -1; point <= 3; ++point) {
      require(
          density(profile, Integer(point)) ==
              expected_grid[static_cast<std::size_t>(point + 1)],
          "integer spectral value mismatch");
    }

    const Rational rational_value =
        density(profile, Rational(Integer(5), Integer(3)));
    require(minimum_shape == 28, "shape margin mismatch");
    require(current == -320, "current mismatch");
    require(
        rational_value ==
            Rational(Integer(-18431), Integer(729)),
        "rational spectral value mismatch");

    std::cout
        << "SU2_INTEGER_SPECTRAL_GRID_OBSTRUCTION"
        << " minimum_shape=" << minimum_shape
        << " grid=-1:1,0:1,1:7,2:1,3:1021"
        << " value_5_over_3=" << rational_value
        << " J_2_3=" << current
        << " result=PASS\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& exception) {
    std::cerr << exception.what() << '\n';
    return EXIT_FAILURE;
  }
}
