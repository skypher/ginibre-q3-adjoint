#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

std::vector<cpp_int> multiply_by_square(const std::vector<cpp_int>& profile,
                                        const int q) {
  const int old_support = static_cast<int>(profile.size()) - 1;
  const int new_support = old_support + 2 * q;
  std::vector<cpp_int> result(static_cast<std::size_t>(new_support + 1));
  for (int a = 0; a <= old_support; ++a) {
    for (int t = 0; t <= 2 * q; ++t) {
      const int low = std::abs(a - t);
      const int high = a + t;
      for (int b = low; b <= high; ++b) {
        result[static_cast<std::size_t>(b)] +=
            profile[static_cast<std::size_t>(a)];
      }
    }
  }
  return result;
}

std::vector<cpp_int> multiply_by_irrep(const std::vector<cpp_int>& profile,
                                       const int q) {
  const int old_support = static_cast<int>(profile.size()) - 1;
  std::vector<cpp_int> result(
      static_cast<std::size_t>(old_support + q + 1));
  for (int source = 0; source <= old_support; ++source) {
    for (int target = std::abs(source - q); target <= source + q;
         ++target) {
      result[static_cast<std::size_t>(target)] +=
          profile[static_cast<std::size_t>(source)];
    }
  }
  return result;
}

cpp_int profile_value(const std::vector<cpp_int>& profile, const int index) {
  return index >= 0 && index < static_cast<int>(profile.size())
             ? profile[static_cast<std::size_t>(index)]
             : cpp_int(0);
}

std::vector<cpp_int> transform(const std::vector<cpp_int>& profile,
                               const int label) {
  std::vector<cpp_int> result(profile.size() +
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

cpp_int wedge(const std::vector<cpp_int>& root,
              const std::vector<cpp_int>& image, const int first,
              const int second) {
  return profile_value(root, first) * profile_value(image, second) -
         profile_value(root, second) * profile_value(image, first);
}

cpp_int kernel_entry(const std::vector<cpp_int>& profile, const int a,
                     const int b) {
  const int low = std::abs(a - b);
  const int high =
      std::min(a + b, static_cast<int>(profile.size()) - 1);
  cpp_int value = 0;
  for (int t = low; t <= high; ++t) {
    value += profile[static_cast<std::size_t>(t)];
  }
  return value;
}

cpp_int b2_coefficient(const std::vector<cpp_int>& profile, const int first,
                       const int second) {
  if (second == 0) {
    return profile_value(profile, 0) *
               (profile_value(profile, first) +
                profile_value(profile, first + 1) +
                profile_value(profile, first + 2)) -
           profile_value(profile, first + 1) *
               profile_value(profile, 1);
  }
  return profile_value(profile, second) *
             (profile_value(profile, first) +
              profile_value(profile, first + 2)) -
         profile_value(profile, first + 1) *
             (profile_value(profile, second - 1) +
              profile_value(profile, second + 1));
}

struct Counters {
  std::uint64_t words = 0;
  std::uint64_t adjacent_minors = 0;
  std::uint64_t root_adjacent_minors = 0;
  std::uint64_t root_adjacent_failures = 0;
  std::uint64_t root_star_minors = 0;
  std::uint64_t root_star_failures = 0;
  std::uint64_t boundary_minors = 0;
  std::uint64_t adjacent_failures = 0;
  std::uint64_t strictly_balanced_adjacent_failures = 0;
  int maximum_adjacent_failure_balance_slack = -1;
  std::uint64_t boundary_failures = 0;
  std::uint64_t ratio_rows = 0;
  std::uint64_t ratio_multiturn_failures = 0;
  std::uint64_t far_endpoint_failures = 0;
  std::uint64_t current_curvatures = 0;
  std::uint64_t current_curvature_failures = 0;
  std::uint64_t exterior_gap_suffixes = 0;
  std::uint64_t exterior_gap_suffix_failures = 0;
  std::uint64_t exterior_pair_products = 0;
  std::uint64_t exterior_pair_product_failures = 0;
  std::uint64_t b2_triangle_rows = 0;
  std::uint64_t b2_triangle_row_failures = 0;
  std::uint64_t b2_triangle_columns = 0;
  std::uint64_t b2_triangle_column_failures = 0;
  std::vector<int> first_adjacent_word;
  std::vector<int> first_root_adjacent_word;
  std::vector<int> first_root_star_word;
  std::vector<int> first_strictly_balanced_adjacent_word;
  std::vector<int> first_boundary_word;
  std::vector<int> first_multiturn_word;
  std::vector<int> first_far_endpoint_word;
  std::vector<int> first_current_curvature_word;
  std::vector<int> first_exterior_gap_suffix_word;
  std::vector<int> first_exterior_pair_product_word;
  std::vector<int> first_b2_triangle_row_word;
  std::vector<int> first_b2_triangle_column_word;
  int first_adjacent_a = -1;
  int first_adjacent_b = -1;
  int first_boundary_a = -1;
  int first_boundary_b = -1;
  cpp_int first_adjacent_northwest = 0;
  cpp_int first_adjacent_northeast = 0;
  cpp_int first_adjacent_southwest = 0;
  cpp_int first_adjacent_southeast = 0;
  cpp_int first_adjacent_value = 0;
  int first_root_adjacent_a = -1;
  int first_root_adjacent_b = -1;
  cpp_int first_root_adjacent_value = 0;
  int first_root_star_a = -1;
  int first_root_star_b = -1;
  cpp_int first_root_star_value = 0;
  cpp_int first_strictly_balanced_adjacent_value = 0;
  cpp_int first_boundary_value = 0;
  int first_multiturn_radius = -1;
  int first_multiturn_target = -1;
  int first_far_endpoint_radius = -1;
  cpp_int first_far_endpoint_value = 0;
  int first_current_curvature_radius = -1;
  int first_current_curvature_target = -1;
  cpp_int first_current_curvature_value = 0;
  int first_exterior_gap_suffix_radius = -1;
  int first_exterior_gap_suffix_target = -1;
  int first_exterior_gap_suffix_gap = -1;
  cpp_int first_exterior_gap_suffix_value = 0;
  int first_exterior_pair_product_radius = -1;
  int first_exterior_pair_product_target = -1;
  int first_exterior_pair_product_first = -1;
  int first_exterior_pair_product_second = -1;
  cpp_int first_exterior_pair_product_value = 0;
  int first_b2_triangle_row_radius = -1;
  int first_b2_triangle_row_target = -1;
  int first_b2_triangle_row = -1;
  cpp_int first_b2_triangle_row_value = 0;
  int first_b2_triangle_column_radius = -1;
  int first_b2_triangle_column_target = -1;
  int first_b2_triangle_column = -1;
  cpp_int first_b2_triangle_column_value = 0;
};

void inspect_word(const std::vector<int>& word, const bool tp2_only,
                  Counters& counters) {
  std::vector<cpp_int> profile{cpp_int(1)};
  std::vector<cpp_int> root{cpp_int(1)};
  for (const int q : word) {
    profile = multiply_by_square(profile, q);
    root = multiply_by_irrep(root, q);
  }
  ++counters.words;

  const int support = static_cast<int>(profile.size()) - 1;
  int sum_before_maximum = 0;
  for (std::size_t index = 0; index + 1U < word.size(); ++index) {
    sum_before_maximum += word[index];
  }
  const bool strictly_balanced =
      word.back() < sum_before_maximum;
  const int endpoint = support + 2;
  const int root_support = static_cast<int>(root.size()) - 1;
  const int root_endpoint = root_support + 2;
  for (int a = 0; a < root_endpoint; ++a) {
    for (int b = 0; b < root_endpoint; ++b) {
      const cpp_int adjacent =
          kernel_entry(root, a, b) * kernel_entry(root, a + 1, b + 1) -
          kernel_entry(root, a, b + 1) * kernel_entry(root, a + 1, b);
      ++counters.root_adjacent_minors;
      if (adjacent < 0) {
        ++counters.root_adjacent_failures;
        if (counters.first_root_adjacent_word.empty()) {
          counters.first_root_adjacent_word = word;
          counters.first_root_adjacent_a = a;
          counters.first_root_adjacent_b = b;
          counters.first_root_adjacent_value = adjacent;
        }
      }
    }
  }
  for (int a = 1; a <= root_endpoint; ++a) {
    for (int b = 1; b <= root_endpoint; ++b) {
      const cpp_int star =
          kernel_entry(root, 0, 0) * kernel_entry(root, a, b) -
          kernel_entry(root, 0, b) * kernel_entry(root, a, 0);
      ++counters.root_star_minors;
      if (star < 0) {
        ++counters.root_star_failures;
        if (counters.first_root_star_word.empty()) {
          counters.first_root_star_word = word;
          counters.first_root_star_a = a;
          counters.first_root_star_b = b;
          counters.first_root_star_value = star;
        }
      }
    }
  }
  for (int a = 0; a <= endpoint; ++a) {
    for (int b = 0; b <= endpoint; ++b) {
      const cpp_int value =
          kernel_entry(profile, a, b) * kernel_entry(profile, a + 1, b + 1) -
          kernel_entry(profile, a, b + 1) * kernel_entry(profile, a + 1, b);
      ++counters.adjacent_minors;
      if (value < 0) {
        ++counters.adjacent_failures;
        counters.maximum_adjacent_failure_balance_slack =
            std::max(counters.maximum_adjacent_failure_balance_slack,
                     sum_before_maximum - word.back());
        if (strictly_balanced) {
          ++counters.strictly_balanced_adjacent_failures;
          if (counters.first_strictly_balanced_adjacent_word.empty()) {
            counters.first_strictly_balanced_adjacent_word = word;
            counters.first_strictly_balanced_adjacent_value = value;
          }
        }
        if (counters.first_adjacent_word.empty()) {
          counters.first_adjacent_word = word;
          counters.first_adjacent_a = a;
          counters.first_adjacent_b = b;
          counters.first_adjacent_northwest =
              kernel_entry(profile, a, b);
          counters.first_adjacent_northeast =
              kernel_entry(profile, a, b + 1);
          counters.first_adjacent_southwest =
              kernel_entry(profile, a + 1, b);
          counters.first_adjacent_southeast =
              kernel_entry(profile, a + 1, b + 1);
          counters.first_adjacent_value = value;
        }
      }
    }
  }

  const cpp_int c0 = profile.front();
  for (int a = 0; a <= support; ++a) {
    for (int b = a; b <= support; ++b) {
      const cpp_int value =
          c0 * kernel_entry(profile, a, b) -
          profile[static_cast<std::size_t>(a)] *
              profile[static_cast<std::size_t>(b)];
      ++counters.boundary_minors;
      if (value < 0) {
        ++counters.boundary_failures;
        if (counters.first_boundary_word.empty()) {
          counters.first_boundary_word = word;
          counters.first_boundary_a = a;
          counters.first_boundary_b = b;
          counters.first_boundary_value = value;
        }
      }
    }
  }

  if (tp2_only) {
    return;
  }

  for (int radius = 0; radius <= support; ++radius) {
    ++counters.ratio_rows;
    bool saw_negative = false;
    bool multiturn = false;
    for (int target = 0; target < support; ++target) {
      const cpp_int cross =
          profile[static_cast<std::size_t>(target)] *
              kernel_entry(profile, radius, target + 1) -
          profile[static_cast<std::size_t>(target + 1)] *
              kernel_entry(profile, radius, target);
      if (cross < 0) {
        saw_negative = true;
      } else if (cross > 0 && saw_negative) {
        multiturn = true;
        if (counters.first_multiturn_word.empty()) {
          counters.first_multiturn_word = word;
          counters.first_multiturn_radius = radius;
          counters.first_multiturn_target = target;
        }
        break;
      }
    }
    if (multiturn) {
      ++counters.ratio_multiturn_failures;
    }

    const cpp_int endpoint_value =
        c0 * kernel_entry(profile, radius, support) -
        profile[static_cast<std::size_t>(radius)] *
            profile[static_cast<std::size_t>(support)];
    if (endpoint_value < 0) {
      ++counters.far_endpoint_failures;
      if (counters.first_far_endpoint_word.empty()) {
        counters.first_far_endpoint_word = word;
        counters.first_far_endpoint_radius = radius;
        counters.first_far_endpoint_value = endpoint_value;
      }
    }
  }

  const auto current = [&](const int radius, const int target) {
    if (radius < 0 || target < 0) {
      return cpp_int(0);
    }
    return c0 * kernel_entry(profile, radius, target) -
           profile[static_cast<std::size_t>(radius)] *
               profile[static_cast<std::size_t>(target)];
  };
  for (int radius = 1; radius <= support; ++radius) {
    for (int target = 1; target <= support; ++target) {
      const cpp_int curvature =
          current(radius, target) - current(radius - 1, target) -
          current(radius, target - 1) +
          current(radius - 1, target - 1);
      ++counters.current_curvatures;
      if (curvature < 0) {
        ++counters.current_curvature_failures;
        if (counters.first_current_curvature_word.empty()) {
          counters.first_current_curvature_word = word;
          counters.first_current_curvature_radius = radius;
          counters.first_current_curvature_target = target;
          counters.first_current_curvature_value = curvature;
        }
      }
    }
  }

  for (int radius = 1; radius <= support; ++radius) {
    for (int target = radius; target <= support; ++target) {
      const int smaller = radius - 1;
      const int larger = target - 1;
      for (int row = 0; row <= smaller; ++row) {
        cpp_int row_sum = 0;
        for (int contraction = 0;
             contraction <= smaller - row; ++contraction) {
          row_sum += b2_coefficient(
              profile, smaller + larger - row - 2 * contraction, row);
        }
        ++counters.b2_triangle_rows;
        if (row_sum < 0) {
          ++counters.b2_triangle_row_failures;
          if (counters.first_b2_triangle_row_word.empty()) {
            counters.first_b2_triangle_row_word = word;
            counters.first_b2_triangle_row_radius = radius;
            counters.first_b2_triangle_row_target = target;
            counters.first_b2_triangle_row = row;
            counters.first_b2_triangle_row_value = row_sum;
          }
        }
      }
      for (int contraction = 0; contraction <= smaller; ++contraction) {
        cpp_int column_sum = 0;
        for (int row = 0; row <= smaller - contraction; ++row) {
          column_sum += b2_coefficient(
              profile, smaller + larger - row - 2 * contraction, row);
        }
        ++counters.b2_triangle_columns;
        if (column_sum < 0) {
          ++counters.b2_triangle_column_failures;
          if (counters.first_b2_triangle_column_word.empty()) {
            counters.first_b2_triangle_column_word = word;
            counters.first_b2_triangle_column_radius = radius;
            counters.first_b2_triangle_column_target = target;
            counters.first_b2_triangle_column = contraction;
            counters.first_b2_triangle_column_value = column_sum;
          }
        }
      }
    }
  }

  std::vector<std::vector<cpp_int>> images(
      static_cast<std::size_t>(support + 1));
  for (int label = 1; label <= support; ++label) {
    images[static_cast<std::size_t>(label)] = transform(root, label);
  }
  for (int radius = 1; radius <= support; ++radius) {
    for (int target = radius + 1; target <= support; ++target) {
      const int exterior_endpoint =
          static_cast<int>(root.size()) - 1 + target;
      std::vector<cpp_int> gap_contributions(
          static_cast<std::size_t>(exterior_endpoint + 1));
      for (int first = 0; first <= exterior_endpoint; ++first) {
        for (int second = first + 1; second <= exterior_endpoint; ++second) {
          const cpp_int product =
              wedge(root, images[static_cast<std::size_t>(radius)], first,
                    second) *
              wedge(root, images[static_cast<std::size_t>(target)], first,
                    second);
          gap_contributions[static_cast<std::size_t>(second - first)] +=
              product;
          ++counters.exterior_pair_products;
          if (product < 0) {
            ++counters.exterior_pair_product_failures;
            if (counters.first_exterior_pair_product_word.empty()) {
              counters.first_exterior_pair_product_word = word;
              counters.first_exterior_pair_product_radius = radius;
              counters.first_exterior_pair_product_target = target;
              counters.first_exterior_pair_product_first = first;
              counters.first_exterior_pair_product_second = second;
              counters.first_exterior_pair_product_value = product;
            }
          }
        }
      }
      cpp_int suffix = 0;
      for (int gap = exterior_endpoint; gap >= 1; --gap) {
        suffix += gap_contributions[static_cast<std::size_t>(gap)];
        ++counters.exterior_gap_suffixes;
        if (suffix < 0) {
          ++counters.exterior_gap_suffix_failures;
          if (counters.first_exterior_gap_suffix_word.empty()) {
            counters.first_exterior_gap_suffix_word = word;
            counters.first_exterior_gap_suffix_radius = radius;
            counters.first_exterior_gap_suffix_target = target;
            counters.first_exterior_gap_suffix_gap = gap;
            counters.first_exterior_gap_suffix_value = suffix;
          }
        }
      }
    }
  }
}

void enumerate_words(const int maximum_q, const int length, const int depth,
                     const int minimum_q, std::vector<int>& word,
                     const bool tp2_only, const bool distinct,
                     Counters& counters) {
  if (depth == length) {
    inspect_word(word, tp2_only, counters);
    return;
  }
  for (int q = minimum_q; q <= maximum_q; ++q) {
    word[static_cast<std::size_t>(depth)] = q;
    enumerate_words(maximum_q, length, depth + 1, distinct ? q + 1 : q,
                    word, tp2_only, distinct, counters);
  }
}

void print_word(const std::vector<int>& word) {
  std::cout << '[';
  for (std::size_t i = 0; i < word.size(); ++i) {
    if (i != 0U) {
      std::cout << ',';
    }
    std::cout << word[i];
  }
  std::cout << ']';
}

void replay_distinct_tp2_obstruction() {
  std::vector<cpp_int> profile{cpp_int(1)};
  for (const int q : {1, 2, 5}) {
    profile = multiply_by_square(profile, q);
  }
  const cpp_int northwest = kernel_entry(profile, 0, 6);
  const cpp_int northeast = kernel_entry(profile, 0, 7);
  const cpp_int southwest = kernel_entry(profile, 1, 6);
  const cpp_int southeast = kernel_entry(profile, 1, 7);
  const cpp_int determinant =
      northwest * southeast - northeast * southwest;
  if (northwest != 220 || northeast != 210 || southwest != 652 ||
      southeast != 622 || determinant != -80) {
    throw std::runtime_error("distinct TP2 obstruction replay mismatch");
  }
  std::cout
      << "SU2_ORDINARY_DISTINCT_TP2_OBSTRUCTION"
      << " word=[1,2,5]"
      << " rows=[0,1]"
      << " columns=[6,7]"
      << " matrix=[220,210;652,622]"
      << " determinant=-80"
      << " result=PASS_EXACT\n";
}

}  // namespace

int main(int argc, char** argv) {
  const std::string mode =
      argc >= 4 ? std::string(argv[3]) : std::string{};
  const bool replay_distinct_tp2 =
      mode == "--replay-distinct-tp2-obstruction";
  const int maximum_q = argc >= 2 ? std::stoi(argv[1]) : 8;
  const int length = argc >= 3 ? std::stoi(argv[2]) : 3;
  const bool tp2_only =
      mode == "--tp2-only" || mode == "--distinct-tp2";
  const bool distinct =
      mode == "--distinct-tp2";
  if (argc >= 4 && !tp2_only && !replay_distinct_tp2) {
    std::cerr
        << "third argument must be --tp2-only, --distinct-tp2, or "
           "--replay-distinct-tp2-obstruction\n";
    return 2;
  }
  if (argc > 4) {
    std::cerr << "usage: " << argv[0]
              << " [maximum_q>=1] [length>=1] "
                 "[(--tp2-only|--distinct-tp2|"
                 "--replay-distinct-tp2-obstruction)]\n";
    return 2;
  }
  if (maximum_q < 1 || length < 1) {
    std::cerr << "usage: " << argv[0] << " [maximum_q>=1] [length>=1]\n";
    return 2;
  }
  if (replay_distinct_tp2) {
    replay_distinct_tp2_obstruction();
    return 0;
  }

  Counters counters;
  std::vector<int> word(static_cast<std::size_t>(length));
  enumerate_words(maximum_q, length, 0, 1, word, tp2_only, distinct,
                  counters);

  std::cout << "SU2_ORDINARY_VARIABLE_BOX_TP2"
            << " maximum_q=" << maximum_q << " length=" << length
            << " tp2_only=" << (tp2_only ? 1 : 0)
            << " distinct=" << (distinct ? 1 : 0)
            << " words=" << counters.words
            << " root_adjacent_minors=" << counters.root_adjacent_minors
            << " root_adjacent_failures="
            << counters.root_adjacent_failures
            << " root_star_minors=" << counters.root_star_minors
            << " root_star_failures=" << counters.root_star_failures
            << " adjacent_minors=" << counters.adjacent_minors
            << " adjacent_failures=" << counters.adjacent_failures
            << " strictly_balanced_adjacent_failures="
            << counters.strictly_balanced_adjacent_failures
            << " maximum_adjacent_failure_balance_slack="
            << counters.maximum_adjacent_failure_balance_slack
            << " boundary_minors=" << counters.boundary_minors
            << " boundary_failures=" << counters.boundary_failures
            << " ratio_rows=" << counters.ratio_rows
            << " ratio_multiturn_failures="
            << counters.ratio_multiturn_failures
            << " far_endpoint_failures=" << counters.far_endpoint_failures
            << " current_curvatures=" << counters.current_curvatures
            << " current_curvature_failures="
            << counters.current_curvature_failures
            << " exterior_gap_suffixes="
            << counters.exterior_gap_suffixes
            << " exterior_gap_suffix_failures="
            << counters.exterior_gap_suffix_failures
            << " exterior_pair_products="
            << counters.exterior_pair_products
            << " exterior_pair_product_failures="
            << counters.exterior_pair_product_failures
            << " b2_triangle_rows=" << counters.b2_triangle_rows
            << " b2_triangle_row_failures="
            << counters.b2_triangle_row_failures
            << " b2_triangle_columns=" << counters.b2_triangle_columns
            << " b2_triangle_column_failures="
            << counters.b2_triangle_column_failures
            << '\n';
  if (!counters.first_root_adjacent_word.empty()) {
    std::cout << "first_root_adjacent_word=";
    print_word(counters.first_root_adjacent_word);
    std::cout << " first_root_adjacent_a=" << counters.first_root_adjacent_a
              << " first_root_adjacent_b=" << counters.first_root_adjacent_b
              << " first_root_adjacent_value="
              << counters.first_root_adjacent_value << '\n';
  }
  if (!counters.first_root_star_word.empty()) {
    std::cout << "first_root_star_word=";
    print_word(counters.first_root_star_word);
    std::cout << " first_root_star_a=" << counters.first_root_star_a
              << " first_root_star_b=" << counters.first_root_star_b
              << " first_root_star_value=" << counters.first_root_star_value
              << '\n';
  }
  if (!counters.first_adjacent_word.empty()) {
    std::cout << "first_adjacent_word=";
    print_word(counters.first_adjacent_word);
    std::cout << " first_adjacent_a=" << counters.first_adjacent_a
              << " first_adjacent_b=" << counters.first_adjacent_b
              << " first_adjacent_entries=["
              << counters.first_adjacent_northwest << ','
              << counters.first_adjacent_northeast << ';'
              << counters.first_adjacent_southwest << ','
              << counters.first_adjacent_southeast << ']'
              << " first_adjacent_value=" << counters.first_adjacent_value
              << '\n';
  }
  if (!counters.first_strictly_balanced_adjacent_word.empty()) {
    std::cout << "first_strictly_balanced_adjacent_word=";
    print_word(counters.first_strictly_balanced_adjacent_word);
    std::cout << " first_strictly_balanced_adjacent_value="
              << counters.first_strictly_balanced_adjacent_value << '\n';
  }
  if (!counters.first_boundary_word.empty()) {
    std::cout << "first_boundary_word=";
    print_word(counters.first_boundary_word);
    std::cout << " first_boundary_a=" << counters.first_boundary_a
              << " first_boundary_b=" << counters.first_boundary_b
              << " first_boundary_value=" << counters.first_boundary_value
              << '\n';
  }
  if (!counters.first_multiturn_word.empty()) {
    std::cout << "first_multiturn_word=";
    print_word(counters.first_multiturn_word);
    std::cout << " first_multiturn_radius="
              << counters.first_multiturn_radius
              << " first_multiturn_target="
              << counters.first_multiturn_target << '\n';
  }
  if (!counters.first_far_endpoint_word.empty()) {
    std::cout << "first_far_endpoint_word=";
    print_word(counters.first_far_endpoint_word);
    std::cout << " first_far_endpoint_radius="
              << counters.first_far_endpoint_radius
              << " first_far_endpoint_value="
              << counters.first_far_endpoint_value << '\n';
  }
  if (!counters.first_current_curvature_word.empty()) {
    std::cout << "first_current_curvature_word=";
    print_word(counters.first_current_curvature_word);
    std::cout << " first_current_curvature_radius="
              << counters.first_current_curvature_radius
              << " first_current_curvature_target="
              << counters.first_current_curvature_target
              << " first_current_curvature_value="
              << counters.first_current_curvature_value << '\n';
  }
  if (!counters.first_exterior_gap_suffix_word.empty()) {
    std::cout << "first_exterior_gap_suffix_word=";
    print_word(counters.first_exterior_gap_suffix_word);
    std::cout << " first_exterior_gap_suffix_radius="
              << counters.first_exterior_gap_suffix_radius
              << " first_exterior_gap_suffix_target="
              << counters.first_exterior_gap_suffix_target
              << " first_exterior_gap_suffix_gap="
              << counters.first_exterior_gap_suffix_gap
              << " first_exterior_gap_suffix_value="
              << counters.first_exterior_gap_suffix_value << '\n';
  }
  if (!counters.first_exterior_pair_product_word.empty()) {
    std::cout << "first_exterior_pair_product_word=";
    print_word(counters.first_exterior_pair_product_word);
    std::cout << " first_exterior_pair_product_radius="
              << counters.first_exterior_pair_product_radius
              << " first_exterior_pair_product_target="
              << counters.first_exterior_pair_product_target
              << " first_exterior_pair_product_first="
              << counters.first_exterior_pair_product_first
              << " first_exterior_pair_product_second="
              << counters.first_exterior_pair_product_second
              << " first_exterior_pair_product_value="
              << counters.first_exterior_pair_product_value << '\n';
  }
  if (!counters.first_b2_triangle_row_word.empty()) {
    std::cout << "first_b2_triangle_row_word=";
    print_word(counters.first_b2_triangle_row_word);
    std::cout << " first_b2_triangle_row_radius="
              << counters.first_b2_triangle_row_radius
              << " first_b2_triangle_row_target="
              << counters.first_b2_triangle_row_target
              << " first_b2_triangle_row="
              << counters.first_b2_triangle_row
              << " first_b2_triangle_row_value="
              << counters.first_b2_triangle_row_value << '\n';
  }
  if (!counters.first_b2_triangle_column_word.empty()) {
    std::cout << "first_b2_triangle_column_word=";
    print_word(counters.first_b2_triangle_column_word);
    std::cout << " first_b2_triangle_column_radius="
              << counters.first_b2_triangle_column_radius
              << " first_b2_triangle_column_target="
              << counters.first_b2_triangle_column_target
              << " first_b2_triangle_column="
              << counters.first_b2_triangle_column
              << " first_b2_triangle_column_value="
              << counters.first_b2_triangle_column_value << '\n';
  }
  return counters.boundary_failures == 0U ? 0 : 1;
}
