#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
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

std::vector<cpp_int> fourier_coefficients(int minus_pairs, int half_power) {
  const int degree = minus_pairs + half_power;
  std::vector<cpp_int> result(static_cast<std::size_t>(2 * degree + 1));
  for (int selected_minus = 0; selected_minus <= 2 * minus_pairs;
       ++selected_minus) {
    for (int selected_plus = 0; selected_plus <= 2 * half_power;
         ++selected_plus) {
      cpp_int value = binomial(2 * minus_pairs, selected_minus) *
                      binomial(2 * half_power, selected_plus);
      if ((minus_pairs + selected_minus) % 2 != 0) {
        value = -value;
      }
      const int exponent = selected_minus + selected_plus - degree;
      result[static_cast<std::size_t>(degree + exponent)] += value;
    }
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

std::vector<int> winding_indices(int modulus, int degree) {
  std::vector<int> result;
  const int bound = (degree + 4) / modulus + 2;
  for (int winding = -bound; winding <= bound; ++winding) {
    bool active = false;
    for (int offset = 1; offset <= 4; ++offset) {
      if (std::abs(offset + winding * modulus) <= degree) {
        active = true;
      }
    }
    if (active) {
      result.push_back(winding);
    }
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  const int maximum_rank = argc > 1 ? std::atoi(argv[1]) : 20;
  const int maximum_minus_pairs = argc > 2 ? std::atoi(argv[2]) : 30;
  const int maximum_half_power = argc > 3 ? std::atoi(argv[3]) : 30;
  if (maximum_rank < 2 || maximum_minus_pairs < 0 ||
      maximum_half_power < 0) {
    std::cerr << "usage: probe_tp2_winding_sectors [maximum-rank>=2] "
                 "[maximum-minus-pairs>=0] [maximum-half-power>=0]\n";
    return 2;
  }

  std::size_t systems = 0;
  std::size_t sectors = 0;
  std::size_t negative_sectors = 0;
  std::size_t negative_sum_groups = 0;
  std::size_t negative_shell_groups = 0;
  std::size_t negative_shell_prefixes = 0;
  std::size_t same_sign_adjacent_coefficients = 0;
  std::size_t negative_cyclic_corrections = 0;
  bool reported = false;
  bool reported_sum = false;
  bool reported_shell = false;
  bool reported_shell_prefix = false;
  for (int rank = 5; rank <= maximum_rank; ++rank) {
    const int modulus = 2 * rank + 1;
    for (int half_power = 2; half_power <= maximum_half_power; ++half_power) {
      for (int minus_pairs = half_power + 2;
           minus_pairs <= maximum_minus_pairs; ++minus_pairs) {
        const std::vector<cpp_int> values =
            fourier_coefficients(minus_pairs, half_power);
        const int degree = minus_pairs + half_power;
        const std::vector<int> windings = winding_indices(modulus, degree);
        cpp_int total = 0;
        std::map<int, cpp_int> sum_groups;
        std::map<int, cpp_int> shell_groups;
        for (int first : windings) {
          for (int second : windings) {
            const auto c = [&](int offset, int winding) {
              return coefficient(values, offset + winding * modulus);
            };
            const cpp_int twice_sector =
                c(2, first) * c(4, second) +
                c(4, first) * c(2, second) +
                2 * c(2, first) * c(2, second) -
                2 * c(3, first) * c(3, second) -
                c(1, first) * c(3, second) -
                c(3, first) * c(1, second);
            ++sectors;
            total += twice_sector;
            sum_groups[first + second] += twice_sector;
            shell_groups[std::max(std::abs(first), std::abs(second))] +=
                twice_sector;
            if (twice_sector < 0) {
              ++negative_sectors;
              if (!reported) {
                reported = true;
                std::cout << "TP2_WINDING first_negative"
                          << " rank=" << rank
                          << " minus_pairs=" << minus_pairs
                          << " half_power=" << half_power
                          << " first=" << first << " second=" << second
                          << " twice_sector=" << twice_sector << '\n';
              }
            }
          }
        }
        for (const auto& [shell, value] : shell_groups) {
          if (value < 0) {
            ++negative_shell_groups;
            if (!reported_shell) {
              reported_shell = true;
              std::cout << "TP2_WINDING first_negative_shell_group"
                        << " rank=" << rank
                        << " minus_pairs=" << minus_pairs
                        << " half_power=" << half_power
                        << " shell=" << shell << " value=" << value
                        << '\n';
            }
          }
        }
        cpp_int shell_prefix = 0;
        for (const auto& [shell, value] : shell_groups) {
          shell_prefix += value;
          if (shell_prefix < 0) {
            ++negative_shell_prefixes;
            if (!reported_shell_prefix) {
              reported_shell_prefix = true;
              std::cout << "TP2_WINDING first_negative_shell_prefix"
                        << " rank=" << rank
                        << " minus_pairs=" << minus_pairs
                        << " half_power=" << half_power
                        << " shell=" << shell
                        << " prefix=" << shell_prefix << '\n';
            }
          }
        }
        for (const auto& [winding_sum, value] : sum_groups) {
          if (value < 0) {
            ++negative_sum_groups;
            if (!reported_sum) {
              reported_sum = true;
              std::cout << "TP2_WINDING first_negative_sum_group"
                        << " rank=" << rank
                        << " minus_pairs=" << minus_pairs
                        << " half_power=" << half_power
                        << " winding_sum=" << winding_sum
                        << " value=" << value << '\n';
            }
          }
        }
        const cpp_int c1 = cyclic_coefficient(values, modulus, 1);
        const cpp_int c2 = cyclic_coefficient(values, modulus, 2);
        const cpp_int c3 = cyclic_coefficient(values, modulus, 3);
        const cpp_int c4 = cyclic_coefficient(values, modulus, 4);
        const cpp_int target = c2 * c4 + c2 * c2 - c3 * c3 - c1 * c3;
        const cpp_int ordinary_c1 = coefficient(values, 1);
        const cpp_int ordinary_c2 = coefficient(values, 2);
        const cpp_int ordinary_c3 = coefficient(values, 3);
        const cpp_int ordinary_c4 = coefficient(values, 4);
        const cpp_int ordinary_target =
            ordinary_c2 * ordinary_c4 + ordinary_c2 * ordinary_c2 -
            ordinary_c3 * ordinary_c3 - ordinary_c1 * ordinary_c3;
        if (total != 2 * target) {
          std::cout << "TP2_WINDING result=DECOMPOSITION_FAIL"
                    << " rank=" << rank
                    << " minus_pairs=" << minus_pairs
                    << " half_power=" << half_power << '\n';
          return 1;
        }
        const std::vector<cpp_int> next_minus_values =
            fourier_coefficients(minus_pairs + 1, half_power);
        const std::vector<cpp_int> next_half_values =
            fourier_coefficients(minus_pairs, half_power + 1);
        const cpp_int minus_christoffel =
            c3 * cyclic_coefficient(next_minus_values, modulus, 2) -
            c2 * cyclic_coefficient(next_minus_values, modulus, 3);
        const cpp_int plus_christoffel =
            c2 * cyclic_coefficient(next_half_values, modulus, 3) -
            c3 * cyclic_coefficient(next_half_values, modulus, 2);
        if (target != minus_christoffel || target != plus_christoffel) {
          std::cout << "TP2_WINDING result=CHRISTOFFEL_FAIL"
                    << " rank=" << rank
                    << " minus_pairs=" << minus_pairs
                    << " half_power=" << half_power << '\n';
          return 1;
        }
        if (c2 * c3 > 0) {
          ++same_sign_adjacent_coefficients;
        }
        if (target < ordinary_target) {
          ++negative_cyclic_corrections;
          if (negative_cyclic_corrections == 1U) {
            std::cout << "TP2_WINDING first_negative_cyclic_correction"
                      << " rank=" << rank
                      << " minus_pairs=" << minus_pairs
                      << " half_power=" << half_power
                      << " cyclic=" << target
                      << " ordinary=" << ordinary_target << '\n';
          }
        }
        ++systems;
        if (total < 0) {
          std::cout << "TP2_WINDING result=TOTAL_FAIL"
                    << " rank=" << rank
                    << " minus_pairs=" << minus_pairs
                    << " half_power=" << half_power
                    << " twice_total=" << total << '\n';
          return 1;
        }
      }
    }
  }
  std::cout << "TP2_WINDING maximum_rank=" << maximum_rank
            << " maximum_minus_pairs=" << maximum_minus_pairs
            << " maximum_half_power=" << maximum_half_power
            << " systems=" << systems << " sectors=" << sectors
            << " negative_sectors=" << negative_sectors
            << " negative_sum_groups=" << negative_sum_groups
            << " negative_shell_groups=" << negative_shell_groups
            << " negative_shell_prefixes=" << negative_shell_prefixes
            << " same_sign_adjacent_coefficients="
            << same_sign_adjacent_coefficients
            << " negative_cyclic_corrections="
            << negative_cyclic_corrections
            << " result=PASS\n";
  return 0;
}
