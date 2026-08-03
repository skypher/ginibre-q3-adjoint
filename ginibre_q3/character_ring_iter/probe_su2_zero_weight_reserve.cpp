#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Profile = std::vector<Integer>;

int parse_positive(const char* text, const std::string& name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value, &consumed, 10);
  if (consumed != value.size() || parsed <= 0) {
    throw std::invalid_argument(name + " must be positive");
  }
  return static_cast<int>(parsed);
}

Profile multiply_character(const Profile& profile, const int label) {
  const int old_support = static_cast<int>(profile.size()) - 1;
  Profile result(static_cast<std::size_t>(old_support + label + 1));
  for (int source = 0; source <= old_support; ++source) {
    if (profile[static_cast<std::size_t>(source)] == 0) {
      continue;
    }
    for (int target = std::abs(source - label); target <= source + label;
         target += 2) {
      result[static_cast<std::size_t>(target)] +=
          profile[static_cast<std::size_t>(source)];
    }
  }
  return result;
}

Profile square_profile(const Profile& profile) {
  const int support = static_cast<int>(profile.size()) - 1;
  Profile result(static_cast<std::size_t>(support + 1));
  for (int left = 0; left <= support; ++left) {
    if (profile[static_cast<std::size_t>(left)] == 0) {
      continue;
    }
    for (int right = 0; right <= support; ++right) {
      if (profile[static_cast<std::size_t>(right)] == 0) {
        continue;
      }
      const Integer weight = profile[static_cast<std::size_t>(left)]
                             * profile[static_cast<std::size_t>(right)];
      for (int output = std::abs(left - right); output <= left + right;
           output += 2) {
        result[static_cast<std::size_t>(output / 2)] += weight;
      }
    }
  }
  return result;
}

Profile square_full_profile(const Profile& profile) {
  const int support = static_cast<int>(profile.size()) - 1;
  Profile result(static_cast<std::size_t>(2 * support + 1));
  for (int left = 0; left <= support; ++left) {
    if (profile[static_cast<std::size_t>(left)] == 0) {
      continue;
    }
    for (int right = 0; right <= support; ++right) {
      if (profile[static_cast<std::size_t>(right)] == 0) {
        continue;
      }
      const Integer weight = profile[static_cast<std::size_t>(left)]
                             * profile[static_cast<std::size_t>(right)];
      for (int output = std::abs(left - right); output <= left + right;
           output += 2) {
        result[static_cast<std::size_t>(output)] += weight;
      }
    }
  }
  return result;
}

Profile bounded_composition_profile(const std::vector<int>& word) {
  Profile result{Integer(1)};
  for (const int label : word) {
    const int old_degree = static_cast<int>(result.size()) - 1;
    Profile next(static_cast<std::size_t>(old_degree + label + 1));
    for (int degree = 0; degree <= old_degree; ++degree) {
      for (int increment = 0; increment <= label; ++increment) {
        next[static_cast<std::size_t>(degree + increment)] +=
            result[static_cast<std::size_t>(degree)];
      }
    }
    result = std::move(next);
  }
  return result;
}

Integer total(const Profile& profile) {
  Integer result = 0;
  for (const Integer& entry : profile) {
    result += entry;
  }
  return result;
}

Integer zero_weight_total(const Profile& profile) {
  Integer result = 0;
  for (std::size_t label = 0U; label < profile.size(); label += 2U) {
    result += profile[label];
  }
  return result;
}

void print_word(const std::vector<int>& word) {
  std::cout << '[';
  for (std::size_t index = 0U; index < word.size(); ++index) {
    if (index != 0U) {
      std::cout << ',';
    }
    std::cout << word[index];
  }
  std::cout << ']';
}

void print_profile(const Profile& profile) {
  std::cout << '[';
  for (std::size_t index = 0U; index < profile.size(); ++index) {
    if (index != 0U) {
      std::cout << ',';
    }
    std::cout << profile[index];
  }
  std::cout << ']';
}

struct Witness {
  bool present = false;
  std::vector<int> word;
  int append = 0;
  Integer invariant = 0;
  Integer zero_weight = 0;
  Integer next_invariant = 0;
  Integer next_zero_weight = 0;
};

struct RatioWitness {
  bool present = false;
  std::vector<int> word;
  int lowering = 0;
  Integer margin = 0;
};

int total_label(const std::vector<int>& word) {
  int result = 0;
  for (const int label : word) {
    result += label;
  }
  return result;
}

bool has_log_concave_half_profile(const Profile& root) {
  if (root.empty() || root[0U] == 0) {
    return false;
  }
  std::size_t end = root.size() - 1U;
  while (end > 0U && root[end] == 0) {
    --end;
  }
  Integer suffix = 0;
  std::vector<Integer> half(end + 1U);
  for (std::size_t reverse = end + 1U; reverse > 0U; --reverse) {
    const std::size_t index = reverse - 1U;
    suffix += root[index];
    half[index] = suffix;
    if (root[index] == 0) {
      return false;
    }
  }
  for (std::size_t index = 1U; index + 1U < half.size(); ++index) {
    if (half[index] * half[index]
        < half[index - 1U] * half[index + 1U]) {
      return false;
    }
  }
  return true;
}

int arbitrary_root_mode(const int maximum_support, const int maximum_entry,
                        const int maximum_append,
                        const bool require_log_concavity) {
  Profile root(static_cast<std::size_t>(maximum_support + 1));
  std::size_t roots = 0U;
  std::size_t failures = 0U;
  Profile witness;
  int witness_append = 0;
  Integer witness_invariant = 0;
  Integer witness_zero_weight = 0;
  Integer witness_next_invariant = 0;
  Integer witness_next_zero_weight = 0;
  const std::function<void(int)> visit = [&](const int index) {
    if (index > maximum_support) {
      bool nonzero = false;
      for (const Integer& entry : root) {
        nonzero = nonzero || entry != 0;
      }
      if (!nonzero) {
        return;
      }
      if (require_log_concavity && !has_log_concave_half_profile(root)) {
        return;
      }
      ++roots;
      const Profile squared = square_full_profile(root);
      const Integer invariant = squared[0U];
      const Integer zero_weight = zero_weight_total(squared);
      for (int append = 1; append <= maximum_append; ++append) {
        const Profile next_root = multiply_character(root, append);
        const Profile next_squared = square_full_profile(next_root);
        const Integer next_invariant = next_squared[0U];
        const Integer next_zero_weight = zero_weight_total(next_squared);
        if (Integer(append + 1) * next_zero_weight * invariant
            < (Integer(append + 1) * zero_weight
               + Integer(append) * invariant) * next_invariant) {
          ++failures;
          if (witness.empty()) {
            witness = root;
            witness_append = append;
            witness_invariant = invariant;
            witness_zero_weight = zero_weight;
            witness_next_invariant = next_invariant;
            witness_next_zero_weight = next_zero_weight;
          }
        }
      }
      return;
    }
    for (int entry = 0; entry <= maximum_entry; ++entry) {
      root[static_cast<std::size_t>(index)] = entry;
      visit(index + 1);
    }
  };
  visit(0);
  std::cout << "SU2_ZERO_WEIGHT_FRACTIONAL_ARBITRARY_ROOT"
            << " maximum_support=" << maximum_support
            << " maximum_entry=" << maximum_entry
            << " maximum_append=" << maximum_append
            << " logconcave=" << (require_log_concavity ? 1 : 0)
            << " roots=" << roots
            << " failures=" << failures;
  if (!witness.empty()) {
    std::cout << " first_root=";
    print_profile(witness);
    std::cout << " append=" << witness_append;
    std::cout << " invariant=" << witness_invariant
              << " zero_weight=" << witness_zero_weight
              << " next_invariant=" << witness_next_invariant
              << " next_zero_weight=" << witness_next_zero_weight;
  }
  std::cout << " result="
            << (failures == 0U ? "PASS_EXACT_BOX" : "COUNTEREXAMPLE")
            << '\n';
  return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}

int geometric_root_mode(const int support, const int maximum_append) {
  Profile root(static_cast<std::size_t>(support + 1));
  Integer value = 1;
  root[static_cast<std::size_t>(support)] = 1;
  if (support > 0) {
    root[static_cast<std::size_t>(support - 1)] = 1;
  }
  for (int index = support - 2; index >= 0; --index) {
    value *= 2;
    root[static_cast<std::size_t>(index)] = value;
  }
  const Profile squared = square_full_profile(root);
  const Integer invariant = squared[0U];
  const Integer zero_weight = zero_weight_total(squared);
  int failures = 0;
  int first_append = 0;
  Integer first_next_invariant = 0;
  Integer first_next_zero_weight = 0;
  for (int append = 1; append <= maximum_append; ++append) {
    const Profile next_root = multiply_character(root, append);
    const Profile next_squared = square_full_profile(next_root);
    const Integer next_invariant = next_squared[0U];
    const Integer next_zero_weight = zero_weight_total(next_squared);
    if (Integer(append + 1) * next_zero_weight * invariant
        < (Integer(append + 1) * zero_weight
           + Integer(append) * invariant) * next_invariant) {
      ++failures;
      if (first_append == 0) {
        first_append = append;
        first_next_invariant = next_invariant;
        first_next_zero_weight = next_zero_weight;
      }
    }
  }
  std::cout << "SU2_ZERO_WEIGHT_FRACTIONAL_GEOMETRIC"
            << " support=" << support
            << " maximum_append=" << maximum_append
            << " root=";
  print_profile(root);
  std::cout << " logconcave="
            << (has_log_concave_half_profile(root) ? 1 : 0)
            << " invariant=" << invariant
            << " zero_weight=" << zero_weight
            << " failures=" << failures;
  if (first_append != 0) {
    std::cout << " first_append=" << first_append
              << " next_invariant=" << first_next_invariant
              << " next_zero_weight=" << first_next_zero_weight;
  }
  std::cout << " result="
            << (failures == 0 ? "PASS_EXACT" : "COUNTEREXAMPLE")
            << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

int parity_geometric_mode(const int support) {
  Profile root(static_cast<std::size_t>(support + 1));
  Integer value = 1;
  root[static_cast<std::size_t>(support)] = 1;
  if (support > 0) {
    root[static_cast<std::size_t>(support - 1)] = 1;
  }
  for (int index = support - 2; index >= 0; --index) {
    value *= 2;
    root[static_cast<std::size_t>(index)] = value;
  }
  const Profile squared = square_profile(root);
  const Integer invariant = squared[0U];
  const Integer zero_weight = total(squared);
  Profile appended(root.size());
  for (std::size_t index = 0U; index < root.size(); ++index) {
    appended[index] = root[index];
    if (index + 1U < root.size()) {
      appended[index] += root[index + 1U];
    }
  }
  const Profile next_squared = square_profile(appended);
  const Integer next_invariant = next_squared[0U];
  const Integer next_zero_weight = total(next_squared);
  const Integer margin =
      2 * next_zero_weight * invariant
      - (2 * zero_weight + invariant) * next_invariant;
  std::cout << "SU2_ZERO_WEIGHT_FRACTIONAL_PARITY_GEOMETRIC"
            << " support=" << support << " root=";
  print_profile(root);
  std::cout << " invariant=" << invariant
            << " zero_weight=" << zero_weight
            << " next_invariant=" << next_invariant
            << " next_zero_weight=" << next_zero_weight
            << " margin=" << margin
            << " result="
            << (margin >= 0 ? "PASS_EXACT" : "COUNTEREXAMPLE")
            << '\n';
  return margin >= 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 5 && std::string{argv[1]} == "--arbitrary-root") {
      return arbitrary_root_mode(
          parse_positive(argv[2], "maximum support"),
          parse_positive(argv[3], "maximum entry"),
          parse_positive(argv[4], "maximum append"), false);
    }
    if (argc == 5 && std::string{argv[1]} == "--logconcave-root") {
      return arbitrary_root_mode(
          parse_positive(argv[2], "maximum support"),
          parse_positive(argv[3], "maximum entry"),
          parse_positive(argv[4], "maximum append"), true);
    }
    if (argc == 4 && std::string{argv[1]} == "--geometric-root") {
      return geometric_root_mode(
          parse_positive(argv[2], "support"),
          parse_positive(argv[3], "maximum append"));
    }
    if (argc == 3 && std::string{argv[1]} == "--parity-geometric") {
      return parity_geometric_mode(parse_positive(argv[2], "support"));
    }
    bool check_append = true;
    int argument_offset = 1;
    if (argc == 4 && std::string{argv[1]} == "--reserve-only") {
      check_append = false;
      argument_offset = 2;
    } else if (argc != 3) {
      throw std::invalid_argument(
          "usage: probe_su2_zero_weight_reserve [--reserve-only] "
          "MAXIMUM_LABEL MAXIMUM_FACTORS | --arbitrary-root "
          "MAXIMUM_SUPPORT MAXIMUM_ENTRY MAXIMUM_APPEND | "
          "--logconcave-root MAXIMUM_SUPPORT MAXIMUM_ENTRY MAXIMUM_APPEND | "
          "--geometric-root SUPPORT MAXIMUM_APPEND | "
          "--parity-geometric SUPPORT");
    }
    const int maximum_label =
        parse_positive(argv[argument_offset], "maximum label");
    const int maximum_factors =
        parse_positive(argv[argument_offset + 1], "maximum factor count");

    std::size_t words = 0U;
    std::size_t reserve_failures = 0U;
    std::size_t total_label_reserve_failures = 0U;
    std::size_t lower_half_ratio_failures = 0U;
    bool total_label_minimum_initialized = false;
    Integer total_label_minimum = 0;
    Witness total_label_minimum_witness;
    bool lower_half_ratio_minimum_initialized = false;
    Integer lower_half_ratio_minimum = 0;
    RatioWitness lower_half_ratio_minimum_witness;
    std::size_t append_failures = 0U;
    std::size_t fractional_append_failures = 0U;
    Witness reserve_witness;
    Witness total_label_witness;
    Witness append_witness;
    std::vector<int> word;

    const std::function<void(const Profile&, int)> visit =
        [&](const Profile& profile, const int minimum_label) {
          if (!word.empty()) {
            const Profile squared = square_profile(profile);
            const Integer invariant = squared[0U];
            const Integer zero_weight = total(squared);
            ++words;
            if (zero_weight
                < Integer(static_cast<int>(word.size()) + 1) * invariant) {
              ++reserve_failures;
              if (!reserve_witness.present) {
                reserve_witness = {true, word, 0, invariant, zero_weight,
                                   0, 0};
              }
            }
            if (zero_weight
                < Integer(total_label(word) + 1) * invariant) {
              ++total_label_reserve_failures;
              if (!total_label_witness.present) {
                total_label_witness = {true, word, 0, invariant,
                                       zero_weight, 0, 0};
              }
            }
            const Integer total_label_margin =
                zero_weight
                - Integer(total_label(word) + 1) * invariant;
            if (!total_label_minimum_initialized
                || total_label_margin < total_label_minimum) {
              total_label_minimum_initialized = true;
              total_label_minimum = total_label_margin;
              total_label_minimum_witness = {true, word, 0, invariant,
                                             zero_weight, 0, 0};
            }
            const int label_sum = total_label(word);
            const Profile compositions = bounded_composition_profile(word);
            for (int lowering = 1; 2 * lowering <= label_sum; ++lowering) {
              const Integer ratio_margin =
                  Integer(label_sum - lowering + 1)
                      * compositions[static_cast<std::size_t>(lowering - 1)]
                  - Integer(lowering)
                      * compositions[static_cast<std::size_t>(lowering)];
              if (ratio_margin < 0) {
                ++lower_half_ratio_failures;
              }
              if (!lower_half_ratio_minimum_initialized
                  || ratio_margin < lower_half_ratio_minimum) {
                lower_half_ratio_minimum_initialized = true;
                lower_half_ratio_minimum = ratio_margin;
                lower_half_ratio_minimum_witness = {
                    true, word, lowering, ratio_margin
                };
              }
            }
            if (check_append) {
              for (int append = 1; append <= maximum_label; ++append) {
                const Profile next_profile =
                    multiply_character(profile, append);
                const Profile next_squared = square_profile(next_profile);
                const Integer next_invariant = next_squared[0U];
                const Integer next_zero_weight = total(next_squared);
                if (next_zero_weight * invariant
                    < (zero_weight + invariant) * next_invariant) {
                  ++append_failures;
                  if (!append_witness.present) {
                    append_witness = {true, word, append, invariant,
                                      zero_weight, next_invariant,
                                      next_zero_weight};
                  }
                }
                if (Integer(append + 1) * next_zero_weight * invariant
                    < (Integer(append + 1) * zero_weight
                       + Integer(append) * invariant)
                          * next_invariant) {
                  ++fractional_append_failures;
                }
              }
            }
          }
          if (static_cast<int>(word.size()) == maximum_factors) {
            return;
          }
          for (int label = minimum_label; label <= maximum_label; ++label) {
            word.push_back(label);
            visit(multiply_character(profile, label), label);
            word.pop_back();
          }
        };

    visit(Profile{Integer(1)}, 1);
    std::cout << "SU2_ZERO_WEIGHT_RESERVE"
              << " maximum_label=" << maximum_label
              << " maximum_factors=" << maximum_factors
              << " words=" << words
              << " reserve_failures=" << reserve_failures
              << " total_label_reserve_failures="
              << total_label_reserve_failures
              << " lower_half_ratio_failures="
              << lower_half_ratio_failures
              << " lower_half_ratio_minimum="
              << lower_half_ratio_minimum
              << " total_label_minimum=" << total_label_minimum
              << " append_checked=" << (check_append ? 1 : 0)
              << " append_failures=" << append_failures
              << " fractional_append_failures="
              << fractional_append_failures;
    if (reserve_witness.present) {
      std::cout << " reserve_first=";
      print_word(reserve_witness.word);
      std::cout << " invariant=" << reserve_witness.invariant
                << " zero_weight=" << reserve_witness.zero_weight;
    }
    if (total_label_witness.present) {
      std::cout << " total_label_reserve_first=";
      print_word(total_label_witness.word);
      std::cout << " invariant=" << total_label_witness.invariant
                << " zero_weight=" << total_label_witness.zero_weight;
    }
    if (total_label_minimum_witness.present) {
      std::cout << " total_label_minimum_word=";
      print_word(total_label_minimum_witness.word);
      std::cout << " invariant=" << total_label_minimum_witness.invariant
                << " zero_weight="
                << total_label_minimum_witness.zero_weight;
    }
    if (lower_half_ratio_minimum_witness.present) {
      std::cout << " lower_half_ratio_minimum_word=";
      print_word(lower_half_ratio_minimum_witness.word);
      std::cout << " lowering="
                << lower_half_ratio_minimum_witness.lowering;
    }
    if (append_witness.present) {
      std::cout << " append_first=";
      print_word(append_witness.word);
      std::cout << '+' << append_witness.append;
      std::cout << " invariant=" << append_witness.invariant
                << " zero_weight=" << append_witness.zero_weight
                << " next_invariant=" << append_witness.next_invariant
                << " next_zero_weight=" << append_witness.next_zero_weight;
    }
    std::cout << " result="
              << (reserve_failures == 0U
                      && total_label_reserve_failures == 0U
                      && lower_half_ratio_failures == 0U
                      && append_failures == 0U
                      ? "PASS_EXACT_BOX"
                      : "COUNTEREXAMPLE")
              << '\n';
    return reserve_failures == 0U
               && total_label_reserve_failures == 0U
               && lower_half_ratio_failures == 0U
               && append_failures == 0U
               ? EXIT_SUCCESS
               : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
