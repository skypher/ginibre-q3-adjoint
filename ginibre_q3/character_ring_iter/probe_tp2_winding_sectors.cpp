#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

std::vector<cpp_int> fourier_coefficients(int minus_pairs, int half_power) {
  const int degree = minus_pairs + half_power;
  std::vector<cpp_int> result(static_cast<std::size_t>(2 * degree + 1));
  const int negative_degree = 2 * minus_pairs;
  const int positive_degree = 2 * half_power;
  cpp_int previous = 0;
  cpp_int current = 1;
  const cpp_int linear = positive_degree - negative_degree;
  for (int coefficient = 0; coefficient <= 2 * degree; ++coefficient) {
    result[static_cast<std::size_t>(coefficient)] =
        (minus_pairs & 1) == 0 ? current : -current;
    if (coefficient == 2 * degree) {
      break;
    }
    const cpp_int next_numerator = linear * current
        + (coefficient - 2 * degree - 1) * previous;
    const int next_denominator = coefficient + 1;
    if (next_numerator % next_denominator != 0) {
      throw std::runtime_error("nonintegral Fourier recurrence");
    }
    previous = current;
    current = next_numerator / next_denominator;
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
  std::size_t negative_shell_suffixes = 0;
  std::size_t negative_physical_image_prefixes = 0;
  std::size_t negative_antiperiodic_targets = 0;
  std::size_t negative_periodic_antiperiodic_sums = 0;
  std::size_t negative_periodic_antiperiodic_differences = 0;
  std::size_t negative_exponent_recurrence_cross_terms = 0;
  std::size_t negative_exponent_cross_box_prefixes = 0;
  std::size_t boxed_identity_checks = 0;
  std::size_t same_sign_adjacent_coefficients = 0;
  std::size_t negative_cyclic_corrections = 0;
  bool reported = false;
  bool reported_sum = false;
  bool reported_shell = false;
  bool reported_shell_prefix = false;
  bool reported_shell_suffix = false;
  bool reported_physical_image_prefix = false;
  bool reported_antiperiodic_target = false;
  bool reported_periodic_antiperiodic_sum = false;
  bool reported_periodic_antiperiodic_difference = false;
  bool reported_exponent_recurrence_cross_term = false;
  bool reported_exponent_cross_box_prefix = false;
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
        std::map<int, std::array<cpp_int, 4>> shell_coefficients;
        std::map<int, std::array<cpp_int, 4>>
            physical_image_coefficients;
        for (const int winding : windings) {
          const int shell = std::abs(winding);
          std::array<cpp_int, 4>& coefficients =
              shell_coefficients[shell];
          for (int offset = 1; offset <= 4; ++offset) {
            const cpp_int value = coefficient(
                values, offset + winding * modulus
            );
            coefficients[static_cast<std::size_t>(offset - 1)] += value;
            const int distance = std::abs(offset + winding * modulus);
            physical_image_coefficients[distance]
                [static_cast<std::size_t>(offset - 1)] += value;
          }
        }
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
        std::array<cpp_int, 4> boxed_coefficients{};
        for (const auto& [shell, value] : shell_groups) {
          shell_prefix += value;
          const auto coefficients = shell_coefficients.find(shell);
          if (coefficients == shell_coefficients.end()) {
            std::cout << "TP2_WINDING result=BOX_SHELL_MISSING"
                      << " rank=" << rank
                      << " minus_pairs=" << minus_pairs
                      << " half_power=" << half_power
                      << " shell=" << shell << '\n';
            return 1;
          }
          for (std::size_t offset = 0U;
               offset < boxed_coefficients.size(); ++offset) {
            boxed_coefficients[offset] += coefficients->second[offset];
          }
          const cpp_int twice_boxed_current = 2 * (
              boxed_coefficients[1U] * boxed_coefficients[3U]
              + boxed_coefficients[1U] * boxed_coefficients[1U]
              - boxed_coefficients[2U] * boxed_coefficients[2U]
              - boxed_coefficients[0U] * boxed_coefficients[2U]
          );
          ++boxed_identity_checks;
          if (shell_prefix != twice_boxed_current) {
            std::cout << "TP2_WINDING result=BOX_IDENTITY_FAIL"
                      << " rank=" << rank
                      << " minus_pairs=" << minus_pairs
                      << " half_power=" << half_power
                      << " shell=" << shell
                      << " sector_sum=" << shell_prefix
                      << " boxed_current=" << twice_boxed_current
                      << '\n';
            return 1;
          }
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
        cpp_int shell_before = 0;
        for (const auto& [shell, value] : shell_groups) {
          const cpp_int shell_suffix = total - shell_before;
          if (shell_suffix < 0) {
            ++negative_shell_suffixes;
            if (!reported_shell_suffix) {
              reported_shell_suffix = true;
              std::cout << "TP2_WINDING first_negative_shell_suffix"
                        << " rank=" << rank
                        << " minus_pairs=" << minus_pairs
                        << " half_power=" << half_power
                        << " shell=" << shell
                        << " suffix=" << shell_suffix << '\n';
            }
          }
          shell_before += value;
        }
        if (shell_before != total) {
          std::cout << "TP2_WINDING result=SHELL_TOTAL_MISMATCH"
                    << " rank=" << rank
                    << " minus_pairs=" << minus_pairs
                    << " half_power=" << half_power << '\n';
          return 1;
        }
        std::array<cpp_int, 4> physical_image_coefficients_sum{};
        for (const auto& [distance, coefficients]
             : physical_image_coefficients) {
          for (std::size_t offset = 0U;
               offset < physical_image_coefficients_sum.size(); ++offset) {
            physical_image_coefficients_sum[offset] += coefficients[offset];
          }
          const cpp_int twice_physical_image_current = 2 * (
              physical_image_coefficients_sum[1U]
                  * physical_image_coefficients_sum[3U]
              + physical_image_coefficients_sum[1U]
                  * physical_image_coefficients_sum[1U]
              - physical_image_coefficients_sum[2U]
                  * physical_image_coefficients_sum[2U]
              - physical_image_coefficients_sum[0U]
                  * physical_image_coefficients_sum[2U]
          );
          if (twice_physical_image_current < 0) {
            ++negative_physical_image_prefixes;
            if (!reported_physical_image_prefix) {
              reported_physical_image_prefix = true;
              std::cout << "TP2_WINDING first_negative_physical_image_prefix"
                        << " rank=" << rank
                        << " minus_pairs=" << minus_pairs
                        << " half_power=" << half_power
                        << " distance=" << distance
                        << " prefix=" << twice_physical_image_current
                        << '\n';
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
        const std::vector<cpp_int> values_after_minus =
            fourier_coefficients(minus_pairs + 1, half_power);
        const std::vector<cpp_int> values_after_plus =
            fourier_coefficients(minus_pairs, half_power + 1);
        const auto cyclic_current = [&](const std::vector<cpp_int>& input) {
          const cpp_int input_c1 = cyclic_coefficient(input, modulus, 1);
          const cpp_int input_c2 = cyclic_coefficient(input, modulus, 2);
          const cpp_int input_c3 = cyclic_coefficient(input, modulus, 3);
          const cpp_int input_c4 = cyclic_coefficient(input, modulus, 4);
          return std::array<cpp_int, 5>{
              input_c1,
              input_c2,
              input_c3,
              input_c4,
              input_c2 * input_c4 + input_c2 * input_c2
                  - input_c3 * input_c3 - input_c1 * input_c3
          };
        };
        const auto after_minus = cyclic_current(values_after_minus);
        const auto after_plus = cyclic_current(values_after_plus);
        const auto exponent_cross = [](
                                        const std::array<cpp_int, 4>& left,
                                        const std::array<cpp_int, 4>& right) {
          return left[1U] * right[3U] + left[3U] * right[1U]
                 + 2 * left[1U] * right[1U]
                 - 2 * left[2U] * right[2U]
                 - left[0U] * right[2U] - left[2U] * right[0U];
        };
        const std::array<cpp_int, 4> after_minus_coefficients{
            after_minus[0U], after_minus[1U], after_minus[2U],
            after_minus[3U]
        };
        const std::array<cpp_int, 4> after_plus_coefficients{
            after_plus[0U], after_plus[1U], after_plus[2U], after_plus[3U]
        };
        const cpp_int twice_exponent_cross = exponent_cross(
            after_minus_coefficients, after_plus_coefficients
        );
        if (16 * target
            != after_minus[4U] + after_plus[4U]
                + twice_exponent_cross) {
          std::cout << "TP2_WINDING result=EXPONENT_RECURRENCE_FAIL"
                    << " rank=" << rank
                    << " minus_pairs=" << minus_pairs
                    << " half_power=" << half_power << '\n';
          return 1;
        }
        if (twice_exponent_cross < 0) {
          ++negative_exponent_recurrence_cross_terms;
          if (!reported_exponent_recurrence_cross_term) {
            reported_exponent_recurrence_cross_term = true;
            std::cout << "TP2_WINDING first_negative_exponent_recurrence_cross"
                      << " rank=" << rank
                      << " minus_pairs=" << minus_pairs
                      << " half_power=" << half_power
                      << " twice_cross=" << twice_exponent_cross << '\n';
          }
        }
        const std::vector<int> cross_windings = winding_indices(
            modulus, degree + 1
        );
        using CrossCoefficientPair = std::pair<
            std::array<cpp_int, 4>, std::array<cpp_int, 4>
        >;
        std::map<int, CrossCoefficientPair> cross_shell_coefficients;
        for (const int winding : cross_windings) {
          CrossCoefficientPair& shell_pair = cross_shell_coefficients[
              std::abs(winding)
          ];
          for (int offset = 1; offset <= 4; ++offset) {
            const std::size_t index = static_cast<std::size_t>(offset - 1);
            shell_pair.first[index] += coefficient(
                values_after_minus, offset + winding * modulus
            );
            shell_pair.second[index] += coefficient(
                values_after_plus, offset + winding * modulus
            );
          }
        }
        std::array<cpp_int, 4> cross_box_left{};
        std::array<cpp_int, 4> cross_box_right{};
        for (const auto& [shell, shell_pair] : cross_shell_coefficients) {
          for (std::size_t index = 0U;
               index < cross_box_left.size(); ++index) {
            cross_box_left[index] += shell_pair.first[index];
            cross_box_right[index] += shell_pair.second[index];
          }
          const cpp_int cross_box_value = exponent_cross(
              cross_box_left, cross_box_right
          );
          if (cross_box_value < 0) {
            ++negative_exponent_cross_box_prefixes;
            if (!reported_exponent_cross_box_prefix) {
              reported_exponent_cross_box_prefix = true;
              std::cout << "TP2_WINDING first_negative_exponent_cross_box"
                        << " rank=" << rank
                        << " minus_pairs=" << minus_pairs
                        << " half_power=" << half_power
                        << " shell=" << shell
                        << " value=" << cross_box_value << '\n';
            }
          }
        }
        const auto antiperiodic_coefficient = [&](int offset) {
          cpp_int result = 0;
          for (const int winding : windings) {
            const cpp_int value = coefficient(
                values, offset + winding * modulus
            );
            result += (winding & 1) == 0 ? value : -value;
          }
          return result;
        };
        const cpp_int anti_c1 = antiperiodic_coefficient(1);
        const cpp_int anti_c2 = antiperiodic_coefficient(2);
        const cpp_int anti_c3 = antiperiodic_coefficient(3);
        const cpp_int anti_c4 = antiperiodic_coefficient(4);
        const cpp_int antiperiodic_target =
            anti_c2 * anti_c4 + anti_c2 * anti_c2
            - anti_c3 * anti_c3 - anti_c1 * anti_c3;
        if (antiperiodic_target < 0) {
          ++negative_antiperiodic_targets;
          if (!reported_antiperiodic_target) {
            reported_antiperiodic_target = true;
            std::cout << "TP2_WINDING first_negative_antiperiodic_target"
                      << " rank=" << rank
                      << " minus_pairs=" << minus_pairs
                      << " half_power=" << half_power
                      << " value=" << antiperiodic_target << '\n';
          }
        }
        if (target + antiperiodic_target < 0) {
          ++negative_periodic_antiperiodic_sums;
          if (!reported_periodic_antiperiodic_sum) {
            reported_periodic_antiperiodic_sum = true;
            std::cout << "TP2_WINDING first_negative_periodic_antiperiodic_sum"
                      << " rank=" << rank
                      << " minus_pairs=" << minus_pairs
                      << " half_power=" << half_power
                      << " periodic=" << target
                      << " antiperiodic=" << antiperiodic_target
                      << " sum=" << target + antiperiodic_target << '\n';
          }
        }
        if (target - antiperiodic_target < 0) {
          ++negative_periodic_antiperiodic_differences;
          if (!reported_periodic_antiperiodic_difference) {
            reported_periodic_antiperiodic_difference = true;
            std::cout
                << "TP2_WINDING first_negative_periodic_antiperiodic_difference"
                << " rank=" << rank
                << " minus_pairs=" << minus_pairs
                << " half_power=" << half_power
                << " periodic=" << target
                << " antiperiodic=" << antiperiodic_target
                << " difference=" << target - antiperiodic_target << '\n';
          }
        }
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
            << " negative_shell_suffixes=" << negative_shell_suffixes
            << " negative_physical_image_prefixes="
            << negative_physical_image_prefixes
            << " negative_antiperiodic_targets="
            << negative_antiperiodic_targets
            << " negative_periodic_antiperiodic_sums="
            << negative_periodic_antiperiodic_sums
            << " negative_periodic_antiperiodic_differences="
            << negative_periodic_antiperiodic_differences
            << " negative_exponent_recurrence_cross_terms="
            << negative_exponent_recurrence_cross_terms
            << " negative_exponent_cross_box_prefixes="
            << negative_exponent_cross_box_prefixes
            << " boxed_identity_checks=" << boxed_identity_checks
            << " same_sign_adjacent_coefficients="
            << same_sign_adjacent_coefficients
            << " negative_cyclic_corrections="
            << negative_cyclic_corrections
            << " result=PASS\n";
  return 0;
}
