#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
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
  if (consumed != value.size() || parsed <= 0) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  return static_cast<int>(parsed);
}

Integer at(const std::vector<Integer>& profile, const int index) {
  return index >= 0 && index < static_cast<int>(profile.size())
             ? profile[static_cast<std::size_t>(index)]
             : Integer(0);
}

Integer first_radial(const std::vector<Integer>& profile) {
  return at(profile, 0) *
             (at(profile, 1) + at(profile, 2) + at(profile, 3)) -
         at(profile, 1) * at(profile, 2);
}

std::vector<Integer> multiply_by_square(
    const std::vector<Integer>& profile, const int label) {
  const int old_support = static_cast<int>(profile.size()) - 1;
  std::vector<Integer> result(
      static_cast<std::size_t>(old_support + label + 1));
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

void enumerate_words(const int maximum_label, const int length,
                     const int depth, const int minimum_label,
                     std::vector<int>& word,
                     std::vector<std::vector<int>>& words) {
  if (depth == length) {
    words.push_back(word);
    return;
  }
  for (int label = minimum_label; label <= maximum_label; ++label) {
    word[static_cast<std::size_t>(depth)] = label;
    enumerate_words(
        maximum_label, length, depth + 1, label, word, words);
  }
}

struct Inspection {
  std::uint64_t insertion_identities = 0U;
  std::uint64_t catalan_identities = 0U;
  std::uint64_t boundary_ratio_checks = 0U;
  std::uint64_t radial_checks = 0U;
  std::uint64_t failures = 0U;
  Integer final_zero = 0;
  Integer final_one = 0;
};

Inspection inspect_word(const std::vector<int>& word) {
  Inspection inspection;
  std::vector<Integer> profile{Integer(1)};
  for (std::size_t factor = 0; factor < word.size(); ++factor) {
    const int label = word[factor];
    const std::vector<Integer> next =
        multiply_by_square(profile, label);

    Integer predicted_zero = 0;
    for (int index = 0; index <= label; ++index) {
      predicted_zero += at(profile, index);
    }
    Integer predicted_one = at(profile, 0);
    for (int index = 1; index <= label - 1; ++index) {
      predicted_one += 3 * at(profile, index);
    }
    predicted_one +=
        2 * at(profile, label) + at(profile, label + 1);
    ++inspection.insertion_identities;
    if (at(next, 0) != predicted_zero ||
        at(next, 1) != predicted_one) {
      ++inspection.failures;
    }

    ++inspection.boundary_ratio_checks;
    if (at(next, 1) < at(next, 0)) {
      ++inspection.failures;
    }
    const bool spectral_threshold_applies =
        (factor >= 1U && label >= 2) ||
        (factor >= 3U && label == 1);
    if (spectral_threshold_applies &&
        at(next, 1) < 2 * at(next, 0)) {
      ++inspection.failures;
    }

    const Integer radial = first_radial(next);
    ++inspection.radial_checks;
    if (radial < 0) {
      ++inspection.failures;
    }
    profile = next;
  }
  inspection.final_zero = at(profile, 0);
  inspection.final_one = at(profile, 1);
  if (std::all_of(
          word.begin(), word.end(),
          [](const int label) { return label == 1; })) {
    Integer catalan = 1;
    for (std::size_t index = 0; index < word.size(); ++index) {
      const Integer step = static_cast<unsigned long long>(index);
      catalan = catalan * 2 * (2 * step + 1) / (step + 2);
    }
    const Integer size = static_cast<unsigned long long>(word.size());
    const Integer next_catalan =
        catalan * 2 * (2 * size + 1) / (size + 2);
    ++inspection.catalan_identities;
    if (inspection.final_zero != catalan ||
        inspection.final_one != next_catalan - catalan) {
      ++inspection.failures;
    }
  }
  return inspection;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const int maximum_label =
        argc >= 2 ? parse_positive(argv[1], "maximum_label") : 8;
    const int maximum_factors =
        argc >= 3 ? parse_positive(argv[2], "maximum_factors") : 8;
    if (argc > 3) {
      throw std::invalid_argument(
          "usage: verify_su2_first_lower_radial_payment "
          "[maximum_label] [maximum_factors]");
    }

    std::vector<std::vector<int>> words;
    for (int length = 1; length <= maximum_factors; ++length) {
      std::vector<int> word(static_cast<std::size_t>(length));
      enumerate_words(
          maximum_label, length, 0, 1, word, words);
    }

    std::atomic<std::uint64_t> insertion_identities{0U};
    std::atomic<std::uint64_t> catalan_identities{0U};
    std::atomic<std::uint64_t> boundary_ratio_checks{0U};
    std::atomic<std::uint64_t> radial_checks{0U};
    std::atomic<std::uint64_t> failures{0U};
    std::vector<Integer> minimum_numerator(
        static_cast<std::size_t>(maximum_factors + 1));
    std::vector<Integer> minimum_denominator(
        static_cast<std::size_t>(maximum_factors + 1));
    std::vector<std::vector<int>> minimum_word(
        static_cast<std::size_t>(maximum_factors + 1));

#pragma omp parallel for schedule(dynamic)
    for (std::size_t index = 0; index < words.size(); ++index) {
      const Inspection inspection = inspect_word(words[index]);
      insertion_identities.fetch_add(
          inspection.insertion_identities, std::memory_order_relaxed);
      catalan_identities.fetch_add(
          inspection.catalan_identities, std::memory_order_relaxed);
      boundary_ratio_checks.fetch_add(
          inspection.boundary_ratio_checks, std::memory_order_relaxed);
      radial_checks.fetch_add(
          inspection.radial_checks, std::memory_order_relaxed);
      failures.fetch_add(
          inspection.failures, std::memory_order_relaxed);
      const std::size_t length = words[index].size();
#pragma omp critical
      {
        const Integer candidate_cross =
            inspection.final_one * minimum_denominator[length];
        const Integer incumbent_cross =
            minimum_numerator[length] * inspection.final_zero;
        if (minimum_word[length].empty() ||
            candidate_cross < incumbent_cross ||
            (candidate_cross == incumbent_cross &&
             words[index] < minimum_word[length])) {
          minimum_numerator[length] = inspection.final_one;
          minimum_denominator[length] = inspection.final_zero;
          minimum_word[length] = words[index];
        }
      }
    }

    const bool pass = failures.load() == 0U;
    std::cout
        << "SU2_FIRST_LOWER_RADIAL_PAYMENT"
        << " maximum_label=" << maximum_label
        << " maximum_factors=" << maximum_factors
        << " words=" << words.size()
        << " insertion_identities=" << insertion_identities.load()
        << " catalan_identities=" << catalan_identities.load()
        << " boundary_ratio_checks=" << boundary_ratio_checks.load()
        << " radial_checks=" << radial_checks.load()
        << " failures=" << failures.load()
        << " minimum_boundary_ratios={";
    for (int length = 1; length <= maximum_factors; ++length) {
      if (length != 1) {
        std::cout << ';';
      }
      const std::size_t index = static_cast<std::size_t>(length);
      std::cout
          << length << ':' << minimum_numerator[index]
          << '/' << minimum_denominator[index] << "@[";
      for (std::size_t label = 0; label < minimum_word[index].size();
           ++label) {
        if (label != 0U) {
          std::cout << ',';
        }
        std::cout << minimum_word[index][label];
      }
      std::cout << ']';
    }
    std::cout
        << "} result=" << (pass ? "PASS" : "FAIL")
        << '\n';
    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
