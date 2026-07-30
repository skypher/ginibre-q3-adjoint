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

Integer second_radial(const std::vector<Integer>& profile) {
  return at(profile, 0) *
             (at(profile, 2) + at(profile, 3) + at(profile, 4)) -
         at(profile, 1) * at(profile, 3);
}

Integer radial(const std::vector<Integer>& profile, const int antidiagonal,
               const int contraction) {
  return at(profile, 0) *
             (at(profile, antidiagonal + 1) +
              at(profile, antidiagonal + 2)) +
         at(profile, contraction) *
             at(profile, antidiagonal - contraction) -
         at(profile, contraction + 1) *
             at(profile, antidiagonal - contraction + 1);
}

bool mean_threshold(const std::vector<Integer>& profile) {
  const Integer zero = at(profile, 0);
  const Integer one = at(profile, 1);
  return 3 * one * one >= 4 * zero * one + 8 * zero * zero;
}

struct DeletionPrefixInspection {
  bool applicable = false;
  bool certified = false;
  Integer best_prefix = 0;
  Integer best_zero = 1;
};

bool deletion_prefix_exception(const std::vector<int>& word) {
  return word == std::vector<int>({1, 1, 2, 2});
}

DeletionPrefixInspection deletion_prefix_inspection(
    const std::vector<int>& word) {
  DeletionPrefixInspection result;
  const int nonfundamental = static_cast<int>(std::count_if(
      word.begin(), word.end(),
      [](const int label) { return label >= 2; }));
  if (word.size() < 4U || nonfundamental < 2) {
    return result;
  }
  result.applicable = true;
  for (std::size_t deleted = 0; deleted < word.size(); ++deleted) {
    const int label = word[deleted];
    if (label < 2) {
      continue;
    }
    std::vector<Integer> profile{Integer(1)};
    for (std::size_t index = 0; index < word.size(); ++index) {
      if (index != deleted) {
        profile = multiply_by_square(profile, word[index]);
      }
    }
    Integer prefix = 0;
    for (int shell = 0; shell <= label; ++shell) {
      prefix += at(profile, shell);
    }
    const Integer zero = at(profile, 0);
    if (result.best_prefix * zero < prefix * result.best_zero) {
      result.best_prefix = prefix;
      result.best_zero = zero;
    }
    if (19 * prefix >= 102 * zero) {
      result.certified = true;
    }
  }
  return result;
}

Integer one_factor_expected(const int label) {
  if (label == 1) {
    return 0;
  }
  if (label <= 3) {
    return 1;
  }
  return 2;
}

Integer two_factor_expected(const int first_label,
                            const int second_label) {
  const int gap = second_label - first_label;
  const int gap_class = std::min(gap, 4);
  if (first_label == 1) {
    constexpr int values[] = {2, 4, 4, 6, 8};
    return values[gap_class];
  }
  if (first_label == 2) {
    constexpr int values[] = {12, 9, 13, 15, 18};
    return values[gap_class];
  }
  if (first_label == 3) {
    constexpr int values[] = {18, 14, 18, 20, 24};
    return values[gap_class];
  }
  switch (gap_class) {
    case 0:
      return 18 * (first_label - 2);
    case 1:
      return 11 * first_label - 19;
    case 2:
      return 11 * first_label - 15;
    case 3:
      return 11 * first_label - 13;
    default:
      return 12 * (first_label - 1);
  }
}

Integer binomial(const int size, int selection) {
  selection = std::min(selection, size - selection);
  Integer result = 1;
  for (int index = 1; index <= selection; ++index) {
    result =
        result * (size - selection + index) / index;
  }
  return result;
}

Integer fundamental_coefficient(const int factors, const int label) {
  if (label > factors) {
    return 0;
  }
  return Integer(2 * label + 1) *
         binomial(2 * factors, factors - label) /
         (factors + label + 1);
}

bool three_factor_exception(const std::vector<int>& word) {
  return (word[0] == 1 && word[1] == 1) ||
         (word[0] == 1 && word[1] == 2 && word[2] == 2) ||
         (word[0] == 2 && word[1] == 2 && word[2] == 2);
}

Integer three_factor_exception_expected(const std::vector<int>& word) {
  if (word[0] == 1 && word[1] == 1) {
    constexpr int values[] = {21, 37, 32, 42, 58, 64};
    return values[std::min(word[2], 6) - 1];
  }
  if (word[0] == 1) {
    return 78;
  }
  return 216;
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

std::string render(const std::vector<int>& word) {
  std::string result{"["};
  for (std::size_t index = 0; index < word.size(); ++index) {
    if (index != 0U) {
      result += ",";
    }
    result += std::to_string(word[index]);
  }
  return result + "]";
}

struct Inspection {
  std::uint64_t insertions = 0U;
  std::uint64_t below_threshold = 0U;
  std::uint64_t second_radial_failures = 0U;
  std::uint64_t third_radial_zero_failures = 0U;
  std::uint64_t third_radial_one_failures = 0U;
  std::uint64_t formula_checks = 0U;
  std::uint64_t formula_failures = 0U;
  std::uint64_t adjacent_drop_checks = 0U;
  std::uint64_t adjacent_drop_failures = 0U;
  std::uint64_t deletion_prefix_candidates = 0U;
  std::uint64_t deletion_prefix_uncovered = 0U;
  std::uint64_t deletion_prefix_unexpected = 0U;
  Integer deletion_prefix_best = 0;
  Integer deletion_prefix_zero = 1;
  Integer final_zero = 0;
  Integer final_one = 0;
  Integer final_second_radial = 0;
  Integer final_third_radial_zero = 0;
  Integer final_third_radial_one = 0;
};

Inspection inspect_word(const std::vector<int>& word) {
  Inspection result;
  std::vector<Integer> profile{Integer(1)};
  for (const int label : word) {
    const std::vector<Integer> next = multiply_by_square(profile, label);
    ++result.insertions;
    profile = next;
  }
  result.final_zero = at(profile, 0);
  result.final_one = at(profile, 1);
  result.final_second_radial = second_radial(profile);
  result.final_third_radial_zero = radial(profile, 3, 0);
  result.final_third_radial_one = radial(profile, 3, 1);
  if (!mean_threshold(profile)) {
    ++result.below_threshold;
  }
  if (result.final_second_radial < 0) {
    ++result.second_radial_failures;
  }
  if (result.final_third_radial_zero < 0) {
    ++result.third_radial_zero_failures;
  }
  if (result.final_third_radial_one < 0) {
    ++result.third_radial_one_failures;
  }
  const int support = static_cast<int>(profile.size()) - 1;
  for (int label = 0; label <= support; ++label) {
    ++result.adjacent_drop_checks;
    if (at(profile, label) - at(profile, label + 1) >
        at(profile, 0)) {
      ++result.adjacent_drop_failures;
    }
  }
  const DeletionPrefixInspection deletion =
      deletion_prefix_inspection(word);
  if (deletion.applicable) {
    ++result.deletion_prefix_candidates;
    result.deletion_prefix_best = deletion.best_prefix;
    result.deletion_prefix_zero = deletion.best_zero;
    if (!deletion.certified) {
      ++result.deletion_prefix_uncovered;
      if (!deletion_prefix_exception(word)) {
        ++result.deletion_prefix_unexpected;
      }
    }
  }
  if (word.size() == 1U) {
    ++result.formula_checks;
    if (result.final_second_radial != one_factor_expected(word[0])) {
      ++result.formula_failures;
    }
  }
  if (word.size() == 2U) {
    ++result.formula_checks;
    if (result.final_second_radial !=
        two_factor_expected(word[0], word[1])) {
      ++result.formula_failures;
    }
  }
  if (word.size() == 3U) {
    ++result.formula_checks;
    if (three_factor_exception(word)) {
      if (result.final_second_radial !=
          three_factor_exception_expected(word)) {
        ++result.formula_failures;
      }
    } else if (2 * result.final_one < 5 * result.final_zero) {
      ++result.formula_failures;
    }
  }
  if (std::all_of(
          word.begin(), word.end(),
          [](const int label) { return label == 1; })) {
    const int factors = static_cast<int>(word.size());
    for (int label = 0; label <= 4; ++label) {
      ++result.formula_checks;
      if (at(profile, label) !=
          fundamental_coefficient(factors, label)) {
        ++result.formula_failures;
      }
    }
    const Integer count = factors;
    const Integer denominator =
        (count + 2) * (count + 2) * (count + 3) *
        (count + 4) * (count + 5);
    const Integer numerator =
        168 * count * (count - 1) * (2 * count + 1);
    ++result.formula_checks;
    if (result.final_second_radial * denominator !=
        result.final_zero * result.final_zero * numerator) {
      ++result.formula_failures;
    }
    const Integer third_denominator =
        (count + 2) * (count + 3) * (count + 4) *
        (count + 5) * (count + 6);
    const Integer third_zero_polynomial =
        (((((2 * count + 37) * count + 256) * count + 803) *
               count +
           1062) *
              count +
          360);
    const Integer third_zero_numerator =
        360 * count * (count - 1) * (count - 2) *
        third_zero_polynomial;
    ++result.formula_checks;
    if (result.final_third_radial_zero *
            third_denominator * third_denominator !=
        result.final_zero * result.final_zero *
            third_zero_numerator) {
      ++result.formula_failures;
    }
    const Integer third_one_polynomial =
        ((((((26 * count + 513) * count + 3836) * count + 13233) *
                  count +
              19808) *
                 count +
             8484) *
                count +
            720);
    const Integer third_one_numerator =
        60 * count * (count - 1) * third_one_polynomial;
    ++result.formula_checks;
    if (result.final_third_radial_one *
            third_denominator * third_denominator !=
        result.final_zero * result.final_zero *
            third_one_numerator) {
      ++result.formula_failures;
    }
  }
  if (word.size() >= 4U && word.size() <= 6U &&
      word.back() == 2 &&
      std::all_of(
          word.begin(), word.end() - 1,
          [](const int label) { return label == 1; })) {
    constexpr int expected[] = {365, 3731, 39340};
    ++result.formula_checks;
    if (result.final_second_radial !=
        expected[word.size() - 4U]) {
      ++result.formula_failures;
    }
  }
  const int fundamental_prefix =
      static_cast<int>(word.size()) - 1;
  if (!word.empty() &&
      std::all_of(
          word.begin(), word.end() - 1,
          [](const int label) { return label == 1; }) &&
      word.back() >= fundamental_prefix + 5) {
    const Integer count = fundamental_prefix;
    const Integer central_binomial =
        binomial(2 * fundamental_prefix, fundamental_prefix);
    const Integer denominator =
        (count + 1) * (count + 2) * (count + 3) *
        (count + 4) * (count + 5);
    const Integer zero_polynomial =
        ((((((((240 * count + 3768) * count + 23808) * count +
                     81648) *
                        count +
                    179952) *
                       count +
                   281400) *
                      count +
                  287904) *
                     count +
                 149280) *
                    count +
                28800);
    ++result.formula_checks;
    if (result.final_third_radial_zero *
            denominator * denominator !=
        central_binomial * central_binomial * zero_polynomial) {
      ++result.formula_failures;
    }
    const Integer one_polynomial =
        ((((((((520 * count + 9172) * count + 65144) * count +
                     240520) *
                        count +
                    501224) *
                       count +
                   610276) *
                      count +
                  442264) *
                     count +
                 175680) *
                    count +
                28800);
    ++result.formula_checks;
    if (result.final_third_radial_one *
            denominator * denominator !=
        central_binomial * central_binomial * one_polynomial) {
      ++result.formula_failures;
    }
  }
  if (deletion_prefix_exception(word)) {
    ++result.formula_checks;
    if (result.final_second_radial != 816) {
      ++result.formula_failures;
    }
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const int maximum_label =
        argc >= 2 ? parse_positive(argv[1], "maximum_label") : 12;
    const int maximum_factors =
        argc >= 3 ? parse_positive(argv[2], "maximum_factors") : 8;
    if (argc > 3) {
      throw std::invalid_argument(
          "usage: probe_su2_second_lower_radial_threshold "
          "[maximum_label] [maximum_factors]");
    }

    const std::vector<Integer> arbitrary_profile{
        Integer(1), Integer(0), Integer(3)};
    const std::vector<Integer> arbitrary_image =
        multiply_by_square(arbitrary_profile, 1);
    const Integer insertion_countercontrol =
        second_radial(arbitrary_image);

    std::vector<std::vector<int>> words;
    for (int length = 1; length <= maximum_factors; ++length) {
      std::vector<int> word(static_cast<std::size_t>(length));
      enumerate_words(
          maximum_label, length, 0, 1, word, words);
    }

    std::atomic<std::uint64_t> insertions{0U};
    std::atomic<std::uint64_t> below_threshold{0U};
    std::atomic<std::uint64_t> second_radial_failures{0U};
    std::atomic<std::uint64_t> third_radial_zero_failures{0U};
    std::atomic<std::uint64_t> third_radial_one_failures{0U};
    std::atomic<std::uint64_t> formula_checks{0U};
    std::atomic<std::uint64_t> formula_failures{0U};
    std::atomic<std::uint64_t> adjacent_drop_checks{0U};
    std::atomic<std::uint64_t> adjacent_drop_failures{0U};
    std::atomic<std::uint64_t> deletion_prefix_candidates{0U};
    std::atomic<std::uint64_t> deletion_prefix_uncovered{0U};
    std::atomic<std::uint64_t> deletion_prefix_unexpected{0U};
    std::vector<int> first_deletion_prefix_uncovered;
    Integer minimum_deletion_prefix = 0;
    Integer minimum_deletion_zero = 1;
    std::vector<int> minimum_deletion_prefix_word;
    std::vector<Integer> minimum_numerator(
        static_cast<std::size_t>(maximum_factors + 1));
    std::vector<Integer> minimum_denominator(
        static_cast<std::size_t>(maximum_factors + 1));
    std::vector<std::vector<int>> minimum_word(
        static_cast<std::size_t>(maximum_factors + 1));
    std::vector<Integer> minimum_third_zero(
        static_cast<std::size_t>(maximum_factors + 1));
    std::vector<Integer> minimum_third_one(
        static_cast<std::size_t>(maximum_factors + 1));
    std::vector<std::vector<int>> minimum_third_zero_word(
        static_cast<std::size_t>(maximum_factors + 1));
    std::vector<std::vector<int>> minimum_third_one_word(
        static_cast<std::size_t>(maximum_factors + 1));

#pragma omp parallel for schedule(dynamic)
    for (std::size_t index = 0; index < words.size(); ++index) {
      const Inspection inspection = inspect_word(words[index]);
      insertions.fetch_add(
          inspection.insertions, std::memory_order_relaxed);
      below_threshold.fetch_add(
          inspection.below_threshold, std::memory_order_relaxed);
      second_radial_failures.fetch_add(
          inspection.second_radial_failures, std::memory_order_relaxed);
      third_radial_zero_failures.fetch_add(
          inspection.third_radial_zero_failures,
          std::memory_order_relaxed);
      third_radial_one_failures.fetch_add(
          inspection.third_radial_one_failures,
          std::memory_order_relaxed);
      formula_checks.fetch_add(
          inspection.formula_checks, std::memory_order_relaxed);
      formula_failures.fetch_add(
          inspection.formula_failures, std::memory_order_relaxed);
      adjacent_drop_checks.fetch_add(
          inspection.adjacent_drop_checks, std::memory_order_relaxed);
      adjacent_drop_failures.fetch_add(
          inspection.adjacent_drop_failures, std::memory_order_relaxed);
      deletion_prefix_candidates.fetch_add(
          inspection.deletion_prefix_candidates,
          std::memory_order_relaxed);
      deletion_prefix_uncovered.fetch_add(
          inspection.deletion_prefix_uncovered,
          std::memory_order_relaxed);
      deletion_prefix_unexpected.fetch_add(
          inspection.deletion_prefix_unexpected,
          std::memory_order_relaxed);
      const std::size_t length = words[index].size();
#pragma omp critical
      {
        if (inspection.deletion_prefix_uncovered != 0U &&
            (first_deletion_prefix_uncovered.empty() ||
             words[index] < first_deletion_prefix_uncovered)) {
          first_deletion_prefix_uncovered = words[index];
        }
        if (inspection.deletion_prefix_candidates != 0U &&
            !deletion_prefix_exception(words[index]) &&
            (minimum_deletion_prefix_word.empty() ||
             inspection.deletion_prefix_best *
                     minimum_deletion_zero <
                 minimum_deletion_prefix *
                     inspection.deletion_prefix_zero)) {
          minimum_deletion_prefix = inspection.deletion_prefix_best;
          minimum_deletion_zero = inspection.deletion_prefix_zero;
          minimum_deletion_prefix_word = words[index];
        }
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
        if (minimum_third_zero_word[length].empty() ||
            inspection.final_third_radial_zero <
                minimum_third_zero[length] ||
            (inspection.final_third_radial_zero ==
                 minimum_third_zero[length] &&
             words[index] < minimum_third_zero_word[length])) {
          minimum_third_zero[length] =
              inspection.final_third_radial_zero;
          minimum_third_zero_word[length] = words[index];
        }
        if (minimum_third_one_word[length].empty() ||
            inspection.final_third_radial_one <
                minimum_third_one[length] ||
            (inspection.final_third_radial_one ==
                 minimum_third_one[length] &&
             words[index] < minimum_third_one_word[length])) {
          minimum_third_one[length] =
              inspection.final_third_radial_one;
          minimum_third_one_word[length] = words[index];
        }
      }
    }

    const bool pass =
        second_radial_failures.load() == 0U &&
        third_radial_zero_failures.load() == 0U &&
        third_radial_one_failures.load() == 0U &&
        formula_failures.load() == 0U &&
        adjacent_drop_failures.load() == 0U &&
        deletion_prefix_unexpected.load() == 0U &&
        insertion_countercontrol == -3;
    std::cout
        << "SU2_SECOND_LOWER_RADIAL_THRESHOLD"
        << " maximum_label=" << maximum_label
        << " maximum_factors=" << maximum_factors
        << " words=" << words.size()
        << " insertions=" << insertions.load()
        << " below_threshold=" << below_threshold.load()
        << " second_radial_failures=" << second_radial_failures.load()
        << " third_radial_zero_failures="
        << third_radial_zero_failures.load()
        << " third_radial_one_failures="
        << third_radial_one_failures.load()
        << " formula_checks=" << formula_checks.load()
        << " formula_failures=" << formula_failures.load()
        << " adjacent_drop_checks=" << adjacent_drop_checks.load()
        << " adjacent_drop_failures="
        << adjacent_drop_failures.load()
        << " deletion_prefix_candidates="
        << deletion_prefix_candidates.load()
        << " deletion_prefix_uncovered="
        << deletion_prefix_uncovered.load()
        << " deletion_prefix_unexpected="
        << deletion_prefix_unexpected.load()
        << " first_deletion_prefix_uncovered="
        << render(first_deletion_prefix_uncovered)
        << " minimum_deletion_prefix_ratio="
        << minimum_deletion_prefix << '/' << minimum_deletion_zero
        << '@' << render(minimum_deletion_prefix_word)
        << " insertion_countercontrol=" << insertion_countercontrol
        << " minimum_means={";
    for (int length = 1; length <= maximum_factors; ++length) {
      if (length != 1) {
        std::cout << ';';
      }
      const std::size_t index = static_cast<std::size_t>(length);
      std::cout
          << length << ':' << minimum_numerator[index]
          << '/' << minimum_denominator[index] << '@'
          << render(minimum_word[index]);
    }
    std::cout
        << "} minimum_third_zero={";
    for (int length = 1; length <= maximum_factors; ++length) {
      if (length != 1) {
        std::cout << ';';
      }
      const std::size_t index = static_cast<std::size_t>(length);
      std::cout << length << ':' << minimum_third_zero[index]
                << '@' << render(minimum_third_zero_word[index]);
    }
    std::cout << "} minimum_third_one={";
    for (int length = 1; length <= maximum_factors; ++length) {
      if (length != 1) {
        std::cout << ';';
      }
      const std::size_t index = static_cast<std::size_t>(length);
      std::cout << length << ':' << minimum_third_one[index]
                << '@' << render(minimum_third_one_word[index]);
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
