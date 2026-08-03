#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Profile = std::vector<Integer>;

Profile multiply_by_square(const Profile& profile, const int label) {
  const int old_support = static_cast<int>(profile.size()) - 1;
  Profile result(static_cast<std::size_t>(old_support + label + 1));
  for (int source = 0; source <= old_support; ++source) {
    for (int shell = 0; shell <= label; ++shell) {
      for (int target = std::abs(source - shell);
           target <= source + shell; ++target) {
        result[static_cast<std::size_t>(target)] +=
            profile[static_cast<std::size_t>(source)];
      }
    }
  }
  return result;
}

Profile fundamental_square_profile(const int fundamentals) {
  Profile result{Integer(1)};
  for (int index = 0; index < fundamentals; ++index) {
    result = multiply_by_square(result, 1);
  }
  return result;
}

Integer factor_coefficient(const int first, const int second,
                           const int shell) {
  Integer result = 0;
  for (int left = 0; left <= first; ++left) {
    const int lower = std::abs(left - shell);
    const int upper = std::min(left + shell, second);
    if (lower <= upper) {
      result += upper - lower + 1;
    }
  }
  return result;
}

Profile direct_factor_profile(const int first, const int gap,
                              const int maximum_shell) {
  Profile result(static_cast<std::size_t>(maximum_shell + 1));
  const int second = first + gap;
  for (int shell = 0; shell <= maximum_shell; ++shell) {
    result[static_cast<std::size_t>(shell)] =
        factor_coefficient(first, second, shell);
  }
  return result;
}

Integer tail_factor_coefficient(const int first, const int gap,
                                const int shell) {
  Integer result = Integer(2 * shell + 1) * first - shell * shell + shell
                   + 1;
  if (gap < shell) {
    const int deficit = shell - gap;
    result -= deficit * (deficit + 1) / 2;
  }
  return result;
}

Profile tail_factor_profile(const int first, const int gap,
                            const int maximum_shell) {
  Profile result(static_cast<std::size_t>(maximum_shell + 1));
  for (int shell = 0; shell <= maximum_shell; ++shell) {
    result[static_cast<std::size_t>(shell)] =
        tail_factor_coefficient(first, gap, shell);
  }
  return result;
}

std::array<Integer, 6U> low_profile(const Profile& fundamental,
                                    const Profile& factor) {
  std::array<Integer, 6U> result{};
  const int fundamental_support = static_cast<int>(fundamental.size()) - 1;
  const int factor_support = static_cast<int>(factor.size()) - 1;
  for (int output = 0; output <= 5; ++output) {
    for (int factor_shell = 0; factor_shell <= factor_support;
         ++factor_shell) {
      for (int source = 0; source <= fundamental_support; ++source) {
        if (std::abs(factor_shell - source) <= output
            && output <= factor_shell + source) {
          result[static_cast<std::size_t>(output)] +=
              factor[static_cast<std::size_t>(factor_shell)]
              * fundamental[static_cast<std::size_t>(source)];
        }
      }
    }
  }
  return result;
}

Integer third_zero(const std::array<Integer, 6U>& profile) {
  return profile[0U] * (profile[3U] + profile[4U] + profile[5U])
         - profile[1U] * profile[4U];
}

Integer third_one(const std::array<Integer, 6U>& profile) {
  return profile[0U] * (profile[4U] + profile[5U])
         + profile[1U] * profile[2U] - profile[2U] * profile[3U];
}

void verify_direct_reconstruction() {
  struct TestCase {
    int fundamentals;
    int first;
    int gap;
  };
  constexpr std::array<TestCase, 4U> tests{{
      {0, 2, 1},
      {6, 2, 1},
      {12, 7, 3},
      {39, 15, 44},
  }};
  for (const TestCase& test : tests) {
    const int maximum_shell = test.fundamentals + 5;
    const Profile fundamental = fundamental_square_profile(test.fundamentals);
    const Profile factor =
        direct_factor_profile(test.first, test.gap, maximum_shell);
    const std::array<Integer, 6U> reconstructed =
        low_profile(fundamental, factor);
    Profile direct = multiply_by_square(fundamental, test.first);
    direct = multiply_by_square(direct, test.first + test.gap);
    for (int shell = 0; shell <= 5; ++shell) {
      if (reconstructed[static_cast<std::size_t>(shell)]
          != direct[static_cast<std::size_t>(shell)]) {
        throw std::runtime_error("direct low-profile reconstruction failed");
      }
    }
  }
}

struct Quadratic {
  Integer constant = 0;
  Integer linear = 0;
  Integer quadratic = 0;
};

Quadratic interpolate_quadratic(const Integer& value_zero,
                                const Integer& value_one,
                                const Integer& value_two) {
  const Integer second_difference = value_two - 2 * value_one + value_zero;
  if (second_difference % 2 != 0) {
    throw std::runtime_error("nonintegral quadratic interpolation");
  }
  const Integer quadratic = second_difference / 2;
  return {value_zero, value_one - value_zero - quadratic, quadratic};
}

void require_nonnegative(const Quadratic& polynomial, const int fundamentals,
                         const int gap, const char* name) {
  if (polynomial.constant < 0 || polynomial.linear < 0
      || polynomial.quadratic < 0) {
    throw std::runtime_error(
        std::string{"negative tail coefficient in "} + name
        + " at m=" + std::to_string(fundamentals)
        + ", gap=" + std::to_string(gap));
  }
}

void observe_minimum(const Integer& value, Integer& minimum,
                     bool& initialized) {
  if (!initialized || value < minimum) {
    minimum = value;
    initialized = true;
  }
}

}  // namespace

int main() {
  try {
    constexpr int maximum_fundamentals = 39;
    verify_direct_reconstruction();
    std::size_t finite_cells = 0U;
    std::size_t tail_cells = 0U;
    Integer minimum_zero = 0;
    Integer minimum_one = 0;
    bool initialized_zero = false;
    bool initialized_one = false;

    for (int fundamentals = 0; fundamentals <= maximum_fundamentals;
         ++fundamentals) {
      const int maximum_shell = fundamentals + 5;
      const int tail_start = 2 * maximum_shell;
      const Profile fundamental = fundamental_square_profile(fundamentals);

      // The finite box includes gap=maximum_shell, which represents every
      // larger gap because all shells through maximum_shell have stabilized.
      for (int first = 2; first < tail_start; ++first) {
        for (int gap = 1; gap <= maximum_shell; ++gap) {
          const Profile factor =
              direct_factor_profile(first, gap, maximum_shell);
          if (gap == maximum_shell
              && factor != direct_factor_profile(
                               first, gap + 1, maximum_shell)) {
            throw std::runtime_error("finite gap stabilization failed");
          }
          const std::array<Integer, 6U> profile =
              low_profile(fundamental, factor);
          const Integer zero = third_zero(profile);
          const Integer one = third_one(profile);
          if (zero < 0 || one < 0) {
            throw std::runtime_error("negative finite two-arbitrary radial");
          }
          observe_minimum(zero, minimum_zero, initialized_zero);
          observe_minimum(one, minimum_one, initialized_one);
          ++finite_cells;
        }
      }

      // For first>=2*maximum_shell the direct lattice count is
      // (2s+1)first-s^2+s+1-T_(max(s-gap,0)) on every needed shell s.
      // Its two radial expressions are therefore quadratics in
      // first-tail_start with exact nonnegative shifted coefficients.
      for (int gap = 1; gap <= maximum_shell; ++gap) {
        std::array<Integer, 3U> zero_values{};
        std::array<Integer, 3U> one_values{};
        for (int offset = 0; offset <= 2; ++offset) {
          const int first = tail_start + offset;
          const Profile direct =
              direct_factor_profile(first, gap, maximum_shell);
          const Profile tail =
              tail_factor_profile(first, gap, maximum_shell);
          if (direct != tail) {
            throw std::runtime_error("tail lattice-count formula failed");
          }
          const std::array<Integer, 6U> profile =
              low_profile(fundamental, tail);
          zero_values[static_cast<std::size_t>(offset)] = third_zero(profile);
          one_values[static_cast<std::size_t>(offset)] = third_one(profile);
        }
        require_nonnegative(interpolate_quadratic(
                                zero_values[0U], zero_values[1U],
                                zero_values[2U]),
                            fundamentals, gap, "K30");
        require_nonnegative(interpolate_quadratic(
                                one_values[0U], one_values[1U],
                                one_values[2U]),
                            fundamentals, gap, "K31");
        ++tail_cells;
      }
    }

    std::cout << "SU2_TWO_ARBITRARY_THIRD_FULL"
              << " maximum_fundamentals=" << maximum_fundamentals
              << " finite_cells=" << finite_cells
              << " tail_cells=" << tail_cells
              << " minimum_K30=" << minimum_zero
              << " minimum_K31=" << minimum_one
              << " result=PASS_EXACT\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
