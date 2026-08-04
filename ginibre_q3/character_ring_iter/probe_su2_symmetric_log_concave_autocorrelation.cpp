#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

int parse_positive(const char* text, const std::string& name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value, &consumed, 10);
  if (consumed != value.size() || parsed <= 0
      || parsed > std::numeric_limits<int>::max()) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  return static_cast<int>(parsed);
}

std::uint64_t splitmix64(std::uint64_t& state) {
  state += UINT64_C(0x9e3779b97f4a7c15);
  std::uint64_t value = state;
  value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31U);
}

Integer half_value(const std::vector<int>& half, const int index) {
  const int absolute = std::abs(index);
  return absolute < static_cast<int>(half.size())
             ? Integer(half[static_cast<std::size_t>(absolute)])
             : Integer(0);
}

Integer toeplitz_minor(const std::vector<int>& half, const int row_first,
                       const int row_second, const int column_first,
                       const int column_second) {
  return half_value(half, column_first - row_first) *
             half_value(half, column_second - row_second) -
         half_value(half, column_second - row_first) *
             half_value(half, column_first - row_second);
}

std::vector<Integer> autocorrelation(const std::vector<int>& half) {
  const int radius = static_cast<int>(half.size()) - 1;
  std::vector<Integer> result(static_cast<std::size_t>(2 * radius + 2));
  for (int shift = 0; shift <= 2 * radius; ++shift) {
    for (int index = -radius; index <= radius; ++index) {
      result[static_cast<std::size_t>(shift)] +=
          half_value(half, index) * half_value(half, index + shift);
    }
  }
  return result;
}

Integer profile_value(const std::vector<Integer>& profile, const int index) {
  return index >= 0 && index < static_cast<int>(profile.size())
             ? profile[static_cast<std::size_t>(index)]
             : Integer(0);
}

Integer current(const std::vector<Integer>& profile, const int radius,
                const int target) {
  Integer interval = 0;
  for (int label = std::abs(radius - target);
       label <= radius + target; ++label) {
    interval += profile_value(profile, label);
  }
  return profile.front() * interval -
         profile_value(profile, radius) * profile_value(profile, target);
}

std::vector<Integer> transform(const std::vector<Integer>& profile,
                               const int label) {
  std::vector<Integer> result(profile.size() +
                              static_cast<std::size_t>(label));
  for (int target = 0; target < static_cast<int>(result.size()); ++target) {
    for (int source = std::abs(target - label);
         source <= target + label; ++source) {
      result[static_cast<std::size_t>(target)] +=
          profile_value(profile, source);
    }
  }
  return result;
}

Integer wedge(const std::vector<Integer>& root,
              const std::vector<Integer>& image, const int first,
              const int second) {
  return profile_value(root, first) * profile_value(image, second) -
         profile_value(root, second) * profile_value(image, first);
}

std::string render(const std::vector<int>& profile) {
  std::string result = "[";
  for (std::size_t index = 0; index < profile.size(); ++index) {
    if (index != 0U) {
      result += ",";
    }
    result += std::to_string(profile[index]);
  }
  return result + "]";
}

struct Counters {
  std::uint64_t profiles = 0;
  std::uint64_t currents = 0;
  std::uint64_t failures = 0;
  std::uint64_t inward_minor_payments = 0;
  std::uint64_t inward_minor_payment_failures = 0;
  std::uint64_t fixed_gap_energies = 0;
  std::uint64_t fixed_gap_failures = 0;
  std::uint64_t gap_suffixes = 0;
  std::uint64_t gap_suffix_failures = 0;
  std::uint64_t cauchy_binet_coefficients = 0;
  std::uint64_t cauchy_binet_coefficient_failures = 0;
  std::uint64_t cauchy_binet_gap_suffixes = 0;
  std::uint64_t cauchy_binet_gap_suffix_failures = 0;
  std::uint64_t reflected_cauchy_binet_pairs = 0;
  std::uint64_t reflected_cauchy_binet_failures = 0;
  std::uint64_t cauchy_binet_identity_checks = 0;
  std::uint64_t cauchy_binet_identity_failures = 0;
  std::vector<int> first_half;
  std::vector<Integer> first_character;
  int first_radius = -1;
  int first_target = -1;
  Integer first_value = 0;
  std::vector<int> first_inward_minor_half;
  int first_inward_minor_radius = -1;
  int first_inward_minor_target = -1;
  int first_inward_minor_side = -1;
  Integer first_inward_minor_value = 0;
  std::vector<int> first_fixed_gap_half;
  int first_fixed_gap_radius = -1;
  int first_fixed_gap_target = -1;
  int first_fixed_gap = -1;
  Integer first_fixed_gap_value = 0;
  std::vector<int> first_gap_suffix_half;
  int first_gap_suffix_radius = -1;
  int first_gap_suffix_target = -1;
  int first_gap_suffix = -1;
  Integer first_gap_suffix_value = 0;
  std::vector<Integer> first_gap_suffix_energies;
  std::vector<int> first_cauchy_binet_half;
  int first_cauchy_binet_radius = -1;
  int first_cauchy_binet_target = -1;
  int first_cauchy_binet_i = 0;
  int first_cauchy_binet_j = 0;
  Integer first_cauchy_binet_outer = 0;
  Integer first_cauchy_binet_inner = 0;
  std::vector<int> first_cauchy_binet_gap_suffix_half;
  int first_cauchy_binet_gap_suffix_radius = -1;
  int first_cauchy_binet_gap_suffix_target = -1;
  int first_cauchy_binet_gap_suffix = -1;
  Integer first_cauchy_binet_gap_suffix_value = 0;
  std::vector<int> first_reflected_cauchy_binet_half;
  int first_reflected_cauchy_binet_radius = -1;
  int first_reflected_cauchy_binet_target = -1;
  int first_reflected_cauchy_binet_i = 0;
  int first_reflected_cauchy_binet_j = 0;
  Integer first_reflected_cauchy_binet_value = 0;
};

void inspect(const std::vector<int>& half, Counters& counters) {
  const std::vector<Integer> weights = autocorrelation(half);
  std::vector<Integer> character(weights.size());
  for (std::size_t index = 0; index < weights.size(); ++index) {
    character[index] =
        weights[index] -
        (index + 1U < weights.size() ? weights[index + 1U] : Integer(0));
  }
  while (!character.empty() && character.back() == 0) {
    character.pop_back();
  }
  std::vector<Integer> root(half.size());
  for (std::size_t index = 0; index < half.size(); ++index) {
    root[index] =
        Integer(half[index]) -
        (index + 1U < half.size() ? Integer(half[index + 1U]) : Integer(0));
  }
  while (!root.empty() && root.back() == 0) {
    root.pop_back();
  }
  ++counters.profiles;
  const int support = static_cast<int>(character.size()) - 1;
  for (int radius = 0; radius <= support; ++radius) {
    for (int target = radius; target <= support; ++target) {
      const Integer value = current(character, radius, target);
      ++counters.currents;
      if (value < 0) {
        ++counters.failures;
        if (counters.first_half.empty()) {
          counters.first_half = half;
          counters.first_character = character;
          counters.first_radius = radius;
          counters.first_target = target;
          counters.first_value = value;
        }
      }
      if (radius >= 1 && target > radius) {
        const Integer left_payment =
            character.front() *
                profile_value(weights, target - radius) -
            profile_value(character, radius) *
                profile_value(weights, target);
        const Integer right_payment =
            profile_value(character, radius) *
                profile_value(weights, target + 1) -
            character.front() *
                profile_value(weights, radius + target + 1);
        for (int side = 0; side < 2; ++side) {
          const Integer payment =
              side == 0 ? left_payment : right_payment;
          ++counters.inward_minor_payments;
          if (payment < 0) {
            ++counters.inward_minor_payment_failures;
            if (counters.first_inward_minor_half.empty()) {
              counters.first_inward_minor_half = half;
              counters.first_inward_minor_radius = radius;
              counters.first_inward_minor_target = target;
              counters.first_inward_minor_side = side;
              counters.first_inward_minor_value = payment;
            }
          }
        }
      }
    }
  }

  std::vector<std::vector<Integer>> images(
      static_cast<std::size_t>(support + 1));
  for (int label = 0; label <= support; ++label) {
    images[static_cast<std::size_t>(label)] = transform(root, label);
  }
  const int endpoint =
      static_cast<int>(root.size()) - 1 + 2 * support;
  for (int radius = 1; radius <= support; ++radius) {
    for (int target = radius + 1; target <= support; ++target) {
      std::vector<Integer> gap_energies(
          static_cast<std::size_t>(endpoint + 1));
      for (int gap = 1; gap <= endpoint; ++gap) {
        Integer energy = 0;
        for (int first = 0; first + gap <= endpoint; ++first) {
          energy +=
              wedge(root, images[static_cast<std::size_t>(radius)], first,
                    first + gap) *
              wedge(root, images[static_cast<std::size_t>(target)], first,
                    first + gap);
        }
        gap_energies[static_cast<std::size_t>(gap)] = energy;
        ++counters.fixed_gap_energies;
        if (energy < 0) {
          ++counters.fixed_gap_failures;
          if (counters.first_fixed_gap_half.empty()) {
            counters.first_fixed_gap_half = half;
            counters.first_fixed_gap_radius = radius;
            counters.first_fixed_gap_target = target;
            counters.first_fixed_gap = gap;
            counters.first_fixed_gap_value = energy;
          }
        }
      }
      Integer suffix = 0;
      for (int gap = endpoint; gap >= 1; --gap) {
        suffix += gap_energies[static_cast<std::size_t>(gap)];
        ++counters.gap_suffixes;
        if (suffix < 0) {
          ++counters.gap_suffix_failures;
          if (counters.first_gap_suffix_half.empty()) {
            counters.first_gap_suffix_half = half;
            counters.first_gap_suffix_radius = radius;
            counters.first_gap_suffix_target = target;
            counters.first_gap_suffix = gap;
            counters.first_gap_suffix_value = suffix;
            counters.first_gap_suffix_energies = gap_energies;
          }
        }
      }

      const int half_radius = static_cast<int>(half.size()) - 1;
      const int last = radius + target + 1;
      Integer cauchy_binet_total = 0;
      std::vector<Integer> cauchy_binet_gap_contributions(
          static_cast<std::size_t>(radius + 2 * half_radius + 1));
      for (int first = -half_radius; first <= radius + half_radius; ++first) {
        for (int second = first + 1; second <= radius + half_radius;
             ++second) {
          const Integer outer =
              toeplitz_minor(half, 0, radius, first, second);
          if (outer <= 0) {
            continue;
          }
          const Integer inner =
              toeplitz_minor(half, first, second, 0, target) -
              toeplitz_minor(half, first, second, -1, target) +
              toeplitz_minor(half, first, second, radius, last) -
              toeplitz_minor(half, first, second, radius + 1, last);
          cauchy_binet_total += outer * inner;
          cauchy_binet_gap_contributions[static_cast<std::size_t>(
              second - first)] += outer * inner;
          ++counters.cauchy_binet_coefficients;
          if (inner < 0) {
            ++counters.cauchy_binet_coefficient_failures;
            if (counters.first_cauchy_binet_half.empty()) {
              counters.first_cauchy_binet_half = half;
              counters.first_cauchy_binet_radius = radius;
              counters.first_cauchy_binet_target = target;
              counters.first_cauchy_binet_i = first;
              counters.first_cauchy_binet_j = second;
              counters.first_cauchy_binet_outer = outer;
              counters.first_cauchy_binet_inner = inner;
            }
          }

          const int reflected_first = radius - second;
          const int reflected_second = radius - first;
          if (first > reflected_first ||
              (first == reflected_first && second > reflected_second)) {
            continue;
          }
          const Integer reflected_outer = toeplitz_minor(
              half, 0, radius, reflected_first, reflected_second);
          const Integer reflected_inner =
              toeplitz_minor(half, reflected_first, reflected_second, 0,
                             target) -
              toeplitz_minor(half, reflected_first, reflected_second, -1,
                             target) +
              toeplitz_minor(half, reflected_first, reflected_second, radius,
                             last) -
              toeplitz_minor(half, reflected_first, reflected_second,
                             radius + 1, last);
          const Integer paired_inner =
              first == reflected_first && second == reflected_second
                  ? inner
                  : inner + reflected_inner;
          ++counters.reflected_cauchy_binet_pairs;
          if (reflected_outer != outer || paired_inner < 0) {
            ++counters.reflected_cauchy_binet_failures;
            if (counters.first_reflected_cauchy_binet_half.empty()) {
              counters.first_reflected_cauchy_binet_half = half;
              counters.first_reflected_cauchy_binet_radius = radius;
              counters.first_reflected_cauchy_binet_target = target;
              counters.first_reflected_cauchy_binet_i = first;
              counters.first_reflected_cauchy_binet_j = second;
              counters.first_reflected_cauchy_binet_value = paired_inner;
            }
          }
        }
      }
      Integer cauchy_binet_gap_suffix = 0;
      for (int gap = radius + 2 * half_radius; gap >= 1; --gap) {
        cauchy_binet_gap_suffix +=
            cauchy_binet_gap_contributions[static_cast<std::size_t>(gap)];
        ++counters.cauchy_binet_gap_suffixes;
        if (cauchy_binet_gap_suffix < 0) {
          ++counters.cauchy_binet_gap_suffix_failures;
          if (counters.first_cauchy_binet_gap_suffix_half.empty()) {
            counters.first_cauchy_binet_gap_suffix_half = half;
            counters.first_cauchy_binet_gap_suffix_radius = radius;
            counters.first_cauchy_binet_gap_suffix_target = target;
            counters.first_cauchy_binet_gap_suffix = gap;
            counters.first_cauchy_binet_gap_suffix_value =
                cauchy_binet_gap_suffix;
          }
        }
      }
      ++counters.cauchy_binet_identity_checks;
      if (cauchy_binet_total != current(character, radius, target)) {
        ++counters.cauchy_binet_identity_failures;
      }
    }
  }
}

void enumerate(const int length, const int maximum_coefficient,
               const int index, std::vector<int>& half, Counters& counters) {
  if (index == length) {
    inspect(half, counters);
    return;
  }
  const int upper =
      index == 0 ? maximum_coefficient
                 : half[static_cast<std::size_t>(index - 1)];
  for (int value = 1; value <= upper; ++value) {
    if (index >= 2) {
      const std::int64_t middle =
          half[static_cast<std::size_t>(index - 1)];
      const std::int64_t previous =
          half[static_cast<std::size_t>(index - 2)];
      if (middle * middle <
          previous * static_cast<std::int64_t>(value)) {
        continue;
      }
    }
    half[static_cast<std::size_t>(index)] = value;
    enumerate(length, maximum_coefficient, index + 1, half, counters);
  }
}

int replay_upper_orthant_obstruction() {
  const std::vector<int> half{5, 4, 3, 2, 1};
  const int radius = 1;
  const int target = 2;
  const int cutoff_first = 2;
  const int cutoff_second = 3;
  const int last = radius + target + 1;
  const int half_radius = static_cast<int>(half.size()) - 1;
  Integer value = 0;
  std::string terms = "[";
  bool first_term = true;
  for (int first = cutoff_first;
       first <= radius + half_radius;
       ++first) {
    for (int second = std::max(first + 1, cutoff_second);
         second <= radius + half_radius;
         ++second) {
      const Integer outer =
          toeplitz_minor(half, 0, radius, first, second);
      if (outer == 0) {
        continue;
      }
      const Integer inner =
          toeplitz_minor(half, first, second, 0, target) -
          toeplitz_minor(half, first, second, -1, target) +
          toeplitz_minor(half, first, second, radius, last) -
          toeplitz_minor(half, first, second, radius + 1, last);
      const Integer contribution = outer * inner;
      value += contribution;
      if (!first_term) {
        terms += ',';
      }
      terms += "(" + std::to_string(first) + ","
               + std::to_string(second) + "):"
               + contribution.convert_to<std::string>();
      first_term = false;
    }
  }
  terms += ']';
  Integer minimum_log_concavity_margin = 0;
  for (std::size_t index = 1U; index + 1U < half.size(); ++index) {
    const Integer margin =
        Integer(half[index]) * half[index] -
        Integer(half[index - 1U]) * half[index + 1U];
    if (index == 1U || margin < minimum_log_concavity_margin) {
      minimum_log_concavity_margin = margin;
    }
  }
  const bool passed =
      value == -2 && minimum_log_concavity_margin == 1;
  std::cout
      << "SU2_CAUCHY_BINET_UPPER_ORTHANT_OBSTRUCTION"
      << " half=" << render(half)
      << " radius=" << radius
      << " target=" << target
      << " cutoff_i=" << cutoff_first
      << " cutoff_j=" << cutoff_second
      << " minimum_log_concavity_margin="
      << minimum_log_concavity_margin
      << " terms=" << terms
      << " orthant_value=" << value
      << " result=" << (passed ? "PASS_EXACT" : "FAIL")
      << '\n';
  return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}

std::vector<int> random_symmetric_log_concave_half(
    std::uint64_t& state, const int maximum_length, const int denominator) {
  if (denominator > std::numeric_limits<int>::max() / 2) {
    throw std::invalid_argument("denominator is too large");
  }
  const int length = 1 + static_cast<int>(
      splitmix64(state) % static_cast<std::uint64_t>(maximum_length));
  std::vector<int> numerators(static_cast<std::size_t>(length - 1));
  for (int& numerator : numerators) {
    numerator = 1 + static_cast<int>(
        splitmix64(state) % static_cast<std::uint64_t>(denominator));
  }
  std::sort(numerators.begin(), numerators.end(), std::greater<int>());

  std::int64_t value = 1;
  for (int index = 1; index < length; ++index) {
    if (value > std::numeric_limits<int>::max() / denominator) {
      throw std::overflow_error("random half profile exceeds int range");
    }
    value *= denominator;
  }
  std::vector<int> half(static_cast<std::size_t>(length));
  half[0] = static_cast<int>(value);
  for (int index = 1; index < length; ++index) {
    value = (value / denominator)
        * numerators[static_cast<std::size_t>(index - 1)];
    half[static_cast<std::size_t>(index)] = static_cast<int>(value);
  }
  return half;
}

std::vector<int> half_from_ratio_numerators(
    const std::vector<int>& numerators, const int denominator) {
  std::int64_t value = 1;
  for (std::size_t index = 0U; index < numerators.size(); ++index) {
    if (value > std::numeric_limits<int>::max() / denominator) {
      throw std::overflow_error("ratio-grid half profile exceeds int range");
    }
    value *= denominator;
  }
  std::vector<int> half(numerators.size() + 1U);
  half[0] = static_cast<int>(value);
  for (std::size_t index = 0U; index < numerators.size(); ++index) {
    value = (value / denominator) * numerators[index];
    half[index + 1U] = static_cast<int>(value);
  }
  return half;
}

struct CurrentSearch {
  std::uint64_t profiles = 0U;
  std::uint64_t currents = 0U;
  std::uint64_t failures = 0U;
  std::vector<int> first_half;
  int first_radius = -1;
  int first_target = -1;
  Integer first_value = 0;
};

void inspect_current_only(const std::vector<int>& half, CurrentSearch& search) {
  const std::vector<Integer> weights = autocorrelation(half);
  std::vector<Integer> character(weights.size());
  for (std::size_t index = 0U; index < weights.size(); ++index) {
    character[index] = weights[index]
        - (index + 1U < weights.size() ? weights[index + 1U] : Integer(0));
  }
  while (!character.empty() && character.back() == 0) {
    character.pop_back();
  }
  ++search.profiles;
  const int support = static_cast<int>(character.size()) - 1;
  for (int radius = 0; radius <= support; ++radius) {
    for (int target = radius; target <= support; ++target) {
      const Integer value = current(character, radius, target);
      ++search.currents;
      if (value < 0) {
        ++search.failures;
        if (search.first_half.empty()) {
          search.first_half = half;
          search.first_radius = radius;
          search.first_target = target;
          search.first_value = value;
        }
      }
    }
  }
}

void print_current_search(
    const std::string& name, const CurrentSearch& search) {
  std::cout << name
            << " profiles=" << search.profiles
            << " currents=" << search.currents
            << " failures=" << search.failures;
  if (!search.first_half.empty()) {
    std::cout << " first_half=" << render(search.first_half)
              << " first_R=" << search.first_radius
              << " first_S=" << search.first_target
              << " first_value=" << search.first_value;
  }
  std::cout << " result="
            << (search.failures == 0U ? "NO_COUNTEREXAMPLE" : "FAIL") << '\n';
}

void enumerate_ratio_grid(const int remaining, const int upper,
                          const int denominator,
                          std::vector<int>& numerators,
                          CurrentSearch& search) {
  if (remaining == 0) {
    inspect_current_only(half_from_ratio_numerators(numerators, denominator),
                         search);
    return;
  }
  for (int value = upper; value >= 1; --value) {
    numerators.push_back(value);
    enumerate_ratio_grid(remaining - 1, value, denominator, numerators,
                         search);
    numerators.pop_back();
  }
}

void ratio_grid_search(const int maximum_length, const int denominator) {
  CurrentSearch search;
  for (int length = 1; length <= maximum_length; ++length) {
    std::vector<int> numerators;
    enumerate_ratio_grid(length - 1, denominator, denominator, numerators,
                         search);
  }
  std::cout << "SU2_SYMMETRIC_LOG_CONCAVE_RATIO_GRID"
            << " maximum_length=" << maximum_length
            << " denominator=" << denominator << ' ';
  print_current_search("", search);
  if (search.failures != 0U) {
    throw std::runtime_error("symmetric log-concave ratio grid failed");
  }
}

void random_current_search(
    const int samples, const int maximum_length, const int denominator) {
  CurrentSearch search;
  std::uint64_t state = UINT64_C(0xbb67ae8584caa73b);
  for (int sample = 0; sample < samples; ++sample) {
    const std::vector<int> half = random_symmetric_log_concave_half(
        state, maximum_length, denominator);
    inspect_current_only(half, search);
  }
  std::cout << "SU2_SYMMETRIC_LOG_CONCAVE_RANDOM"
            << " samples=" << samples
            << " maximum_length=" << maximum_length
            << " denominator=" << denominator << ' ';
  print_current_search("", search);
  if (search.failures != 0U) {
    throw std::runtime_error("random symmetric log-concave search failed");
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc >= 3 && std::string(argv[1]) == "--profile") {
      std::vector<int> half;
      half.reserve(static_cast<std::size_t>(argc - 2));
      for (int index = 2; index < argc; ++index) {
        half.push_back(parse_positive(argv[index], "profile entry"));
      }
      Counters counters;
      inspect(half, counters);
      std::cout << "SU2_SYMMETRIC_LOG_CONCAVE_PROFILE"
                << " half=" << render(half)
                << " currents=" << counters.currents
                << " failures=" << counters.failures;
      if (!counters.first_half.empty()) {
        std::cout << " first_radius=" << counters.first_radius
                  << " first_target=" << counters.first_target
                  << " first_value=" << counters.first_value;
      }
      std::cout << '\n';
      return EXIT_SUCCESS;
    }
    if (argc == 2 &&
        std::string(argv[1]) == "--replay-upper-orthant") {
      return replay_upper_orthant_obstruction();
    }
    if (argc == 5 && std::string(argv[1]) == "--random-current") {
      random_current_search(
          parse_positive(argv[2], "samples"),
          parse_positive(argv[3], "maximum_length"),
          parse_positive(argv[4], "denominator"));
      return EXIT_SUCCESS;
    }
    if (argc == 4 && std::string(argv[1]) == "--ratio-grid") {
      ratio_grid_search(
          parse_positive(argv[2], "maximum_length"),
          parse_positive(argv[3], "denominator"));
      return EXIT_SUCCESS;
    }
    const int maximum_length =
        argc >= 2 ? parse_positive(argv[1], "maximum_length") : 7;
    const int maximum_coefficient =
        argc >= 3 ? parse_positive(argv[2], "maximum_coefficient") : 8;
    if (argc > 3) {
      throw std::invalid_argument(
          "usage: probe_su2_symmetric_log_concave_autocorrelation "
          "[maximum_length] [maximum_coefficient]\n"
          "       probe_su2_symmetric_log_concave_autocorrelation "
          "--replay-upper-orthant\n"
          "       probe_su2_symmetric_log_concave_autocorrelation "
          "--random-current SAMPLES MAXIMUM_LENGTH DENOMINATOR\n"
          "       probe_su2_symmetric_log_concave_autocorrelation "
          "--ratio-grid MAXIMUM_LENGTH DENOMINATOR\n"
          "       probe_su2_symmetric_log_concave_autocorrelation "
          "--profile H_0 H_1 ... H_m");
    }

    Counters counters;
    for (int length = 1; length <= maximum_length; ++length) {
      std::vector<int> half(static_cast<std::size_t>(length));
      enumerate(length, maximum_coefficient, 0, half, counters);
    }

    std::cout << "SU2_SYMMETRIC_LOG_CONCAVE_AUTOCORRELATION"
              << " maximum_length=" << maximum_length
              << " maximum_coefficient=" << maximum_coefficient
              << " profiles=" << counters.profiles
              << " currents=" << counters.currents
              << " failures=" << counters.failures
              << " inward_minor_payments="
              << counters.inward_minor_payments
              << " inward_minor_payment_failures="
              << counters.inward_minor_payment_failures
              << " fixed_gap_energies=" << counters.fixed_gap_energies
              << " fixed_gap_failures=" << counters.fixed_gap_failures
              << " gap_suffixes=" << counters.gap_suffixes
              << " gap_suffix_failures=" << counters.gap_suffix_failures
              << " cauchy_binet_coefficients="
              << counters.cauchy_binet_coefficients
              << " cauchy_binet_coefficient_failures="
              << counters.cauchy_binet_coefficient_failures
              << " cauchy_binet_gap_suffixes="
              << counters.cauchy_binet_gap_suffixes
              << " cauchy_binet_gap_suffix_failures="
              << counters.cauchy_binet_gap_suffix_failures
              << " reflected_cauchy_binet_pairs="
              << counters.reflected_cauchy_binet_pairs
              << " reflected_cauchy_binet_failures="
              << counters.reflected_cauchy_binet_failures
              << " cauchy_binet_identity_checks="
              << counters.cauchy_binet_identity_checks
              << " cauchy_binet_identity_failures="
              << counters.cauchy_binet_identity_failures
              << '\n';
    if (!counters.first_half.empty()) {
      std::cout << "first_half=" << render(counters.first_half)
                << " first_radius=" << counters.first_radius
                << " first_target=" << counters.first_target
                << " first_value=" << counters.first_value
                << " first_character=[";
      for (std::size_t index = 0; index < counters.first_character.size();
           ++index) {
        if (index != 0U) {
          std::cout << ',';
        }
        std::cout << counters.first_character[index];
      }
      std::cout << "]\n";
    }
    if (!counters.first_inward_minor_half.empty()) {
      std::cout << "first_inward_minor_half="
                << render(counters.first_inward_minor_half)
                << " first_inward_minor_radius="
                << counters.first_inward_minor_radius
                << " first_inward_minor_target="
                << counters.first_inward_minor_target
                << " first_inward_minor_side="
                << counters.first_inward_minor_side
                << " first_inward_minor_value="
                << counters.first_inward_minor_value << '\n';
    }
    if (!counters.first_fixed_gap_half.empty()) {
      std::cout << "first_fixed_gap_half="
                << render(counters.first_fixed_gap_half)
                << " first_fixed_gap_radius="
                << counters.first_fixed_gap_radius
                << " first_fixed_gap_target="
                << counters.first_fixed_gap_target
                << " first_fixed_gap=" << counters.first_fixed_gap
                << " first_fixed_gap_value="
                << counters.first_fixed_gap_value << '\n';
    }
    if (!counters.first_gap_suffix_half.empty()) {
      std::cout << "first_gap_suffix_half="
                << render(counters.first_gap_suffix_half)
                << " first_gap_suffix_radius="
                << counters.first_gap_suffix_radius
                << " first_gap_suffix_target="
                << counters.first_gap_suffix_target
                << " first_gap_suffix=" << counters.first_gap_suffix
                << " first_gap_suffix_value="
                << counters.first_gap_suffix_value
                << " first_gap_suffix_nonzero_energies=[";
      bool printed = false;
      for (std::size_t gap = 1;
           gap < counters.first_gap_suffix_energies.size(); ++gap) {
        if (counters.first_gap_suffix_energies[gap] == 0) {
          continue;
        }
        if (printed) {
          std::cout << ',';
        }
        std::cout << gap << ':'
                  << counters.first_gap_suffix_energies[gap];
        printed = true;
      }
      std::cout << "]\n";
    }
    if (!counters.first_cauchy_binet_half.empty()) {
      std::cout << "first_cauchy_binet_half="
                << render(counters.first_cauchy_binet_half)
                << " first_cauchy_binet_radius="
                << counters.first_cauchy_binet_radius
                << " first_cauchy_binet_target="
                << counters.first_cauchy_binet_target
                << " first_cauchy_binet_i="
                << counters.first_cauchy_binet_i
                << " first_cauchy_binet_j="
                << counters.first_cauchy_binet_j
                << " first_cauchy_binet_outer="
                << counters.first_cauchy_binet_outer
                << " first_cauchy_binet_inner="
                << counters.first_cauchy_binet_inner << '\n';
    }
    if (!counters.first_cauchy_binet_gap_suffix_half.empty()) {
      std::cout << "first_cauchy_binet_gap_suffix_half="
                << render(counters.first_cauchy_binet_gap_suffix_half)
                << " first_cauchy_binet_gap_suffix_radius="
                << counters.first_cauchy_binet_gap_suffix_radius
                << " first_cauchy_binet_gap_suffix_target="
                << counters.first_cauchy_binet_gap_suffix_target
                << " first_cauchy_binet_gap_suffix="
                << counters.first_cauchy_binet_gap_suffix
                << " first_cauchy_binet_gap_suffix_value="
                << counters.first_cauchy_binet_gap_suffix_value << '\n';
    }
    if (!counters.first_reflected_cauchy_binet_half.empty()) {
      std::cout << "first_reflected_cauchy_binet_half="
                << render(counters.first_reflected_cauchy_binet_half)
                << " first_reflected_cauchy_binet_radius="
                << counters.first_reflected_cauchy_binet_radius
                << " first_reflected_cauchy_binet_target="
                << counters.first_reflected_cauchy_binet_target
                << " first_reflected_cauchy_binet_i="
                << counters.first_reflected_cauchy_binet_i
                << " first_reflected_cauchy_binet_j="
                << counters.first_reflected_cauchy_binet_j
                << " first_reflected_cauchy_binet_value="
                << counters.first_reflected_cauchy_binet_value << '\n';
    }
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
