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

struct Counters {
  std::uint64_t words = 0;
  std::uint64_t adjacent_minors = 0;
  std::uint64_t boundary_minors = 0;
  std::uint64_t adjacent_failures = 0;
  std::uint64_t boundary_failures = 0;
  std::uint64_t ratio_rows = 0;
  std::uint64_t ratio_multiturn_failures = 0;
  std::uint64_t far_endpoint_failures = 0;
  std::uint64_t current_curvatures = 0;
  std::uint64_t current_curvature_failures = 0;
  std::vector<int> first_adjacent_word;
  std::vector<int> first_boundary_word;
  std::vector<int> first_multiturn_word;
  std::vector<int> first_far_endpoint_word;
  std::vector<int> first_current_curvature_word;
  int first_adjacent_a = -1;
  int first_adjacent_b = -1;
  int first_boundary_a = -1;
  int first_boundary_b = -1;
  cpp_int first_adjacent_northwest = 0;
  cpp_int first_adjacent_northeast = 0;
  cpp_int first_adjacent_southwest = 0;
  cpp_int first_adjacent_southeast = 0;
  cpp_int first_adjacent_value = 0;
  cpp_int first_boundary_value = 0;
  int first_multiturn_radius = -1;
  int first_multiturn_target = -1;
  int first_far_endpoint_radius = -1;
  cpp_int first_far_endpoint_value = 0;
  int first_current_curvature_radius = -1;
  int first_current_curvature_target = -1;
  cpp_int first_current_curvature_value = 0;
};

void inspect_word(const std::vector<int>& word, Counters& counters) {
  std::vector<cpp_int> profile{cpp_int(1)};
  for (const int q : word) {
    profile = multiply_by_square(profile, q);
  }
  ++counters.words;

  const int support = static_cast<int>(profile.size()) - 1;
  const int endpoint = support + 2;
  for (int a = 0; a <= endpoint; ++a) {
    for (int b = 0; b <= endpoint; ++b) {
      const cpp_int value =
          kernel_entry(profile, a, b) * kernel_entry(profile, a + 1, b + 1) -
          kernel_entry(profile, a, b + 1) * kernel_entry(profile, a + 1, b);
      ++counters.adjacent_minors;
      if (value < 0) {
        ++counters.adjacent_failures;
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
}

void enumerate_words(const int maximum_q, const int length, const int depth,
                     const int minimum_q, std::vector<int>& word,
                     Counters& counters) {
  if (depth == length) {
    inspect_word(word, counters);
    return;
  }
  for (int q = minimum_q; q <= maximum_q; ++q) {
    word[static_cast<std::size_t>(depth)] = q;
    enumerate_words(maximum_q, length, depth + 1, q, word, counters);
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

}  // namespace

int main(int argc, char** argv) {
  const int maximum_q = argc >= 2 ? std::stoi(argv[1]) : 8;
  const int length = argc >= 3 ? std::stoi(argv[2]) : 3;
  if (maximum_q < 1 || length < 1) {
    std::cerr << "usage: " << argv[0] << " [maximum_q>=1] [length>=1]\n";
    return 2;
  }

  Counters counters;
  std::vector<int> word(static_cast<std::size_t>(length));
  enumerate_words(maximum_q, length, 0, 1, word, counters);

  std::cout << "SU2_ORDINARY_VARIABLE_BOX_TP2"
            << " maximum_q=" << maximum_q << " length=" << length
            << " words=" << counters.words
            << " adjacent_minors=" << counters.adjacent_minors
            << " adjacent_failures=" << counters.adjacent_failures
            << " boundary_minors=" << counters.boundary_minors
            << " boundary_failures=" << counters.boundary_failures
            << " ratio_rows=" << counters.ratio_rows
            << " ratio_multiturn_failures="
            << counters.ratio_multiturn_failures
            << " far_endpoint_failures=" << counters.far_endpoint_failures
            << " current_curvatures=" << counters.current_curvatures
            << " current_curvature_failures="
            << counters.current_curvature_failures
            << '\n';
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
  return counters.boundary_failures == 0U ? 0 : 1;
}
