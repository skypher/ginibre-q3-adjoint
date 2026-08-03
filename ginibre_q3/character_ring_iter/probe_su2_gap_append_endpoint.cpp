#include <algorithm>
#include <cstddef>
#include <cstdlib>
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
  if (consumed != value.size() || parsed <= 0L
      || parsed > std::numeric_limits<int>::max()) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  return static_cast<int>(parsed);
}

Integer value_at(const std::vector<Integer>& profile, const int index) {
  return index >= 0 && index < static_cast<int>(profile.size())
             ? profile[static_cast<std::size_t>(index)]
             : Integer(0);
}

std::vector<Integer> transform(const std::vector<Integer>& profile,
                               const int label) {
  std::vector<Integer> result(profile.size()
                              + static_cast<std::size_t>(label));
  for (int target = 0; target < static_cast<int>(result.size()); ++target) {
    for (int source = std::abs(target - label); source <= target + label;
         ++source) {
      result[static_cast<std::size_t>(target)] += value_at(profile, source);
    }
  }
  return result;
}

Integer gap_prefix(const std::vector<Integer>& profile, const int first_label,
                   const int second_label, const int maximum_gap) {
  const std::vector<Integer> first = transform(profile, first_label);
  const std::vector<Integer> second = transform(profile, second_label);
  const int endpoint = static_cast<int>(first.size()) - 1;
  Integer result = 0;
  for (int i = 0; i <= endpoint; ++i) {
    for (int j = i + 1; j <= endpoint && j - i <= maximum_gap; ++j) {
      const Integer first_wedge =
          value_at(profile, i) * value_at(first, j)
          - value_at(profile, j) * value_at(first, i);
      const Integer second_wedge =
          value_at(profile, i) * value_at(second, j)
          - value_at(profile, j) * value_at(second, i);
      result += first_wedge * second_wedge;
    }
  }
  return result;
}

Integer fourth_power(const Integer& value) {
  const Integer square = value * value;
  return square * square;
}

std::string render(const std::vector<int>& profile) {
  std::string result = "[";
  for (std::size_t index = 0U; index < profile.size(); ++index) {
    if (index != 0U) {
      result += ',';
    }
    result += std::to_string(profile[index]);
  }
  return result + ']';
}

struct Search {
  int maximum_length;
  int maximum_entry;
  int denominator;
  std::uint64_t profiles = 0U;
  std::uint64_t cases = 0U;
  std::uint64_t interior_below_endpoints = 0U;
  std::uint64_t negative_values = 0U;
  bool have_endpoint_witness = false;
  bool have_negative_witness = false;
  std::vector<int> endpoint_profile;
  std::vector<int> negative_profile;
  int endpoint_r = 0;
  int endpoint_s = 0;
  int endpoint_d = 0;
  int endpoint_sample = 0;
  int negative_r = 0;
  int negative_s = 0;
  int negative_d = 0;
  int negative_sample = 0;
  Integer endpoint_value = 0;
  Integer endpoint_zero = 0;
  Integer endpoint_top = 0;
  Integer negative_value = 0;
};

void inspect(const std::vector<int>& old, Search& search) {
  ++search.profiles;
  const int old_length = static_cast<int>(old.size());
  const int numerator = old.back() * old.back();
  const int denominator = old[static_cast<std::size_t>(old_length - 2)];
  const int bound_denominator = denominator * search.denominator;
  std::vector<Integer> zero(old.begin(), old.end());
  std::vector<Integer> endpoint;
  endpoint.reserve(old.size() + 1U);
  for (const int entry : old) {
    endpoint.push_back(Integer(entry) * denominator);
  }
  endpoint.push_back(numerator);

  for (int first_label = 1; first_label <= 2 * old_length; ++first_label) {
    for (int second_label = first_label + 1;
         second_label <= 2 * old_length; ++second_label) {
      for (int maximum_gap = 1;
           maximum_gap <= old_length + first_label; ++maximum_gap) {
        const Integer zero_value =
            gap_prefix(zero, first_label, second_label, maximum_gap);
        const Integer endpoint_value =
            gap_prefix(endpoint, first_label, second_label, maximum_gap);
        ++search.cases;
        for (int sample = 1; sample < search.denominator; ++sample) {
          std::vector<Integer> trial;
          trial.reserve(old.size() + 1U);
          for (const int entry : old) {
            trial.push_back(Integer(entry) * bound_denominator);
          }
          trial.push_back(Integer(numerator) * sample);
          const Integer trial_value =
              gap_prefix(trial, first_label, second_label, maximum_gap);
          const Integer scaled_zero =
              zero_value * fourth_power(Integer(bound_denominator));
          const Integer scaled_endpoint =
              endpoint_value * fourth_power(Integer(search.denominator));
          if (trial_value < 0) {
            ++search.negative_values;
            if (!search.have_negative_witness) {
              search.have_negative_witness = true;
              search.negative_profile = old;
              search.negative_r = first_label;
              search.negative_s = second_label;
              search.negative_d = maximum_gap;
              search.negative_sample = sample;
              search.negative_value = trial_value;
            }
          }
          if (trial_value < scaled_zero && trial_value < scaled_endpoint) {
            ++search.interior_below_endpoints;
            if (!search.have_endpoint_witness) {
              search.have_endpoint_witness = true;
              search.endpoint_profile = old;
              search.endpoint_r = first_label;
              search.endpoint_s = second_label;
              search.endpoint_d = maximum_gap;
              search.endpoint_sample = sample;
              search.endpoint_value = trial_value;
              search.endpoint_zero = scaled_zero;
              search.endpoint_top = scaled_endpoint;
            }
          }
        }
      }
    }
  }
}

void enumerate(const int target_length, const int maximum_entry,
               std::vector<int>& profile, Search& search) {
  if (static_cast<int>(profile.size()) == target_length) {
    inspect(profile, search);
    return;
  }
  const int upper = profile.empty() ? maximum_entry : profile.back();
  for (int entry = 1; entry <= upper; ++entry) {
    if (profile.size() >= 2U) {
      const Integer middle = profile.back();
      const Integer previous = profile[profile.size() - 2U];
      if (middle * middle < previous * entry) {
        continue;
      }
    }
    profile.push_back(entry);
    enumerate(target_length, maximum_entry, profile, search);
    profile.pop_back();
  }
}

int replay_flat_endpoint_obstruction() {
  const std::vector<Integer> zero{1, 1, 0};
  const std::vector<Integer> endpoint{1, 1, 1};
  const std::vector<Integer> interior{16, 16, 1};
  constexpr int first_label = 1;
  constexpr int second_label = 2;
  constexpr int maximum_gap = 1;
  const Integer zero_value =
      gap_prefix(zero, first_label, second_label, maximum_gap);
  const Integer endpoint_value =
      gap_prefix(endpoint, first_label, second_label, maximum_gap);
  const Integer interior_value =
      gap_prefix(interior, first_label, second_label, maximum_gap);
  const Integer expected_interior =
      3 * fourth_power(Integer(16)) - Integer(16) * Integer(16) * Integer(16)
      - 3 * Integer(16) * Integer(16) + 2 * Integer(16) + 2;
  const bool passed = zero_value == 3 && endpoint_value == 3
                      && interior_value == expected_interior
                      && interior_value < zero_value * fourth_power(Integer(16))
                      && interior_value
                             < endpoint_value * fourth_power(Integer(16));
  std::cout << "SU2_GAP_APPEND_FLAT_ENDPOINT_OBSTRUCTION"
            << " old=(1,1)"
            << " R=1 S=2 D=1"
            << " polynomial=3-x-3x^2+2x^3+2x^4"
            << " x=1/16"
            << " P0=" << zero_value
            << " P1=" << endpoint_value
            << " scaled_P_1_16=" << interior_value
            << " result=" << (passed ? "PASS_EXACT" : "FAIL") << '\n';
  return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1])
                         == "--replay-flat-endpoint-obstruction") {
      return replay_flat_endpoint_obstruction();
    }
    if (argc != 4) {
      throw std::invalid_argument(
          "usage: probe_su2_gap_append_endpoint "
          "<maximum-old-length> <maximum-entry> <interior-denominator>\n"
          "   or: probe_su2_gap_append_endpoint "
          "--replay-flat-endpoint-obstruction");
    }
    Search search{};
    search.maximum_length =
        parse_positive(argv[1], "maximum-old-length");
    search.maximum_entry = parse_positive(argv[2], "maximum-entry");
    search.denominator =
        parse_positive(argv[3], "interior-denominator");
    if (search.maximum_length < 2) {
      throw std::invalid_argument("maximum-old-length must be at least two");
    }
    if (search.denominator < 2) {
      throw std::invalid_argument("interior-denominator must be at least two");
    }
    for (int length = 2; length <= search.maximum_length; ++length) {
      std::vector<int> profile;
      enumerate(length, search.maximum_entry, profile, search);
    }
    std::cout << "SU2_GAP_APPEND_ENDPOINT"
              << " maximum_old_length=" << search.maximum_length
              << " maximum_entry=" << search.maximum_entry
              << " interior_denominator=" << search.denominator
              << " profiles=" << search.profiles
              << " cases=" << search.cases
              << " interior_below_endpoints="
              << search.interior_below_endpoints
              << " negative_values=" << search.negative_values;
    if (search.have_endpoint_witness) {
      std::cout << " endpoint_witness_old=" << render(search.endpoint_profile)
                << " R=" << search.endpoint_r
                << " S=" << search.endpoint_s
                << " D=" << search.endpoint_d
                << " sample=" << search.endpoint_sample
                << " trial=" << search.endpoint_value
                << " zero=" << search.endpoint_zero
                << " top=" << search.endpoint_top;
    }
    if (search.have_negative_witness) {
      std::cout << " negative_witness_old=" << render(search.negative_profile)
                << " R=" << search.negative_r
                << " S=" << search.negative_s
                << " D=" << search.negative_d
                << " sample=" << search.negative_sample
                << " trial=" << search.negative_value;
    }
    const bool passed = search.interior_below_endpoints == 0U
                        && search.negative_values == 0U;
    std::cout << " result=" << (passed ? "NO_COUNTEREXAMPLE" : "COUNTEREXAMPLE")
              << '\n';
    return passed ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
