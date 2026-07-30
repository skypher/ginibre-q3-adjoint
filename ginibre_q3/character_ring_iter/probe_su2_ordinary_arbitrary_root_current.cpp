#include <algorithm>
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

Integer value(const std::vector<Integer>& profile, const int index) {
  return index >= 0 && index < static_cast<int>(profile.size())
             ? profile[static_cast<std::size_t>(index)]
             : Integer(0);
}

std::vector<Integer> square_character(const std::vector<int>& root) {
  const int support = static_cast<int>(root.size()) - 1;
  std::vector<Integer> result(static_cast<std::size_t>(2 * support + 1));
  for (int first = 0; first <= support; ++first) {
    for (int second = 0; second <= support; ++second) {
      const Integer weight =
          Integer(root[static_cast<std::size_t>(first)]) *
          Integer(root[static_cast<std::size_t>(second)]);
      for (int label = std::abs(first - second);
           label <= first + second; ++label) {
        result[static_cast<std::size_t>(label)] += weight;
      }
    }
  }
  return result;
}

Integer current(const std::vector<Integer>& square, const int radius,
                const int target) {
  Integer interval = 0;
  for (int label = std::abs(radius - target);
       label <= radius + target; ++label) {
    interval += value(square, label);
  }
  return square.front() * interval -
         value(square, radius) * value(square, target);
}

Integer b2_coefficient(const std::vector<Integer>& square, const int first,
                       const int second) {
  if (second == 0) {
    return value(square, 0) *
               (value(square, first) + value(square, first + 1) +
                value(square, first + 2)) -
           value(square, first + 1) * value(square, 1);
  }
  return value(square, second) *
             (value(square, first) + value(square, first + 2)) -
         value(square, first + 1) *
             (value(square, second - 1) + value(square, second + 1));
}

Integer b2_triangle(const std::vector<Integer>& square, const int radius,
                    const int target) {
  const int smaller = std::min(radius, target) - 1;
  const int larger = std::max(radius, target) - 1;
  Integer result = 0;
  for (int row = 0; row <= smaller; ++row) {
    for (int contraction = 0; contraction <= smaller - row;
         ++contraction) {
      result += b2_coefficient(
          square, smaller + larger - row - 2 * contraction, row);
    }
  }
  return result;
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

void replay_square_log_concave_radial_obstruction() {
  const std::vector<int> root{3, 1, 5, 3, 5, 2, 2};
  const std::vector<Integer> expected{
      77, 172, 284, 313, 339, 294, 254,
      179, 125, 68, 36, 12, 4};
  const std::vector<Integer> square = square_character(root);
  if (square != expected) {
    throw std::runtime_error("square-profile replay mismatch");
  }
  Integer minimum_margin =
      square[1] * square[1] - square[0] * square[2];
  for (std::size_t index = 1; index + 1U < square.size(); ++index) {
    const Integer margin =
        square[index] * square[index] -
        square[index - 1U] * square[index + 1U];
    minimum_margin = std::min(minimum_margin, margin);
  }
  if (minimum_margin < 0) {
    throw std::runtime_error("square profile is not log concave");
  }
  const Integer radial =
      square[0] * (square[5] + square[6]) +
      square[1] * square[3] - square[2] * square[4];
  if (radial != -244) {
    throw std::runtime_error("radial-column replay mismatch");
  }
  std::cout
      << "SU2_ORDINARY_SQUARE_LC_RADIAL_OBSTRUCTION"
      << " root=" << render(root)
      << " square=[77,172,284,313,339,294,254,179,125,68,36,12,4]"
      << " minimum_log_concavity_margin=" << minimum_margin
      << " A=4 L=1 value=" << radial
      << " result=PASS_EXACT\n";
}

struct Counters {
  std::uint64_t roots = 0;
  std::uint64_t currents = 0;
  std::uint64_t failures = 0;
  std::uint64_t log_concave_roots = 0;
  std::uint64_t log_concave_currents = 0;
  std::uint64_t log_concave_failures = 0;
  std::uint64_t interval_support_roots = 0;
  std::uint64_t interval_support_currents = 0;
  std::uint64_t interval_support_failures = 0;
  std::uint64_t b2_identity_failures = 0;
  std::uint64_t negative_b2_coefficients = 0;
  std::uint64_t square_log_concave_roots = 0;
  std::uint64_t square_log_concave_b2_columns = 0;
  std::uint64_t square_log_concave_b2_column_failures = 0;
  std::vector<int> first_root;
  std::vector<Integer> first_square;
  int first_radius = -1;
  int first_target = -1;
  Integer first_value = 0;
  std::vector<int> first_log_concave_root;
  int first_log_concave_radius = -1;
  int first_log_concave_target = -1;
  Integer first_log_concave_value = 0;
  std::vector<int> first_interval_support_root;
  int first_interval_support_radius = -1;
  int first_interval_support_target = -1;
  Integer first_interval_support_value = 0;
  std::vector<int> first_square_log_concave_b2_column_root;
  std::vector<Integer> first_square_log_concave_b2_column_square;
  int first_square_log_concave_b2_column_radius = -1;
  int first_square_log_concave_b2_column_target = -1;
  int first_square_log_concave_b2_column = -1;
  Integer first_square_log_concave_b2_column_value = 0;
};

void inspect(const std::vector<int>& root, Counters& counters) {
  const std::vector<Integer> square = square_character(root);
  ++counters.roots;
  bool log_concave = true;
  bool interval_support = true;
  bool saw_positive = false;
  bool saw_zero_after_positive = false;
  for (int index = 0; index < static_cast<int>(root.size()); ++index) {
    const int coefficient = root[static_cast<std::size_t>(index)];
    if (coefficient > 0) {
      if (saw_zero_after_positive) {
        log_concave = false;
        interval_support = false;
      }
      saw_positive = true;
    } else if (saw_positive) {
      saw_zero_after_positive = true;
    }
    if (index >= 1 && index + 1 < static_cast<int>(root.size())) {
      const std::int64_t middle = coefficient;
      if (middle * middle <
          static_cast<std::int64_t>(
              root[static_cast<std::size_t>(index - 1)]) *
              root[static_cast<std::size_t>(index + 1)]) {
        log_concave = false;
      }
    }
  }
  if (log_concave) {
    ++counters.log_concave_roots;
  }
  if (interval_support) {
    ++counters.interval_support_roots;
  }
  bool square_log_concave = true;
  bool square_saw_zero = false;
  for (int index = 0; index < static_cast<int>(square.size()); ++index) {
    if (square[static_cast<std::size_t>(index)] == 0) {
      square_saw_zero = true;
    } else if (square_saw_zero) {
      square_log_concave = false;
    }
    if (index >= 1 && index + 1 < static_cast<int>(square.size()) &&
        square[static_cast<std::size_t>(index)] *
                square[static_cast<std::size_t>(index)] <
            square[static_cast<std::size_t>(index - 1)] *
                square[static_cast<std::size_t>(index + 1)]) {
      square_log_concave = false;
    }
  }
  if (square_log_concave) {
    ++counters.square_log_concave_roots;
  }
  const int support = static_cast<int>(square.size()) - 1;
  for (int first = 0; first <= support; ++first) {
    for (int second = 0; second <= first; ++second) {
      if (b2_coefficient(square, first, second) < 0) {
        ++counters.negative_b2_coefficients;
      }
    }
  }
  for (int radius = 0; radius <= support; ++radius) {
    for (int target = radius; target <= support; ++target) {
      const Integer margin = current(square, radius, target);
      ++counters.currents;
      if (log_concave) {
        ++counters.log_concave_currents;
        if (margin < 0) {
          ++counters.log_concave_failures;
          if (counters.first_log_concave_root.empty()) {
            counters.first_log_concave_root = root;
            counters.first_log_concave_radius = radius;
            counters.first_log_concave_target = target;
            counters.first_log_concave_value = margin;
          }
        }
      }
      if (interval_support) {
        ++counters.interval_support_currents;
        if (margin < 0) {
          ++counters.interval_support_failures;
          if (counters.first_interval_support_root.empty()) {
            counters.first_interval_support_root = root;
            counters.first_interval_support_radius = radius;
            counters.first_interval_support_target = target;
            counters.first_interval_support_value = margin;
          }
        }
      }
      if (radius >= 1 &&
          b2_triangle(square, radius, target) != margin) {
        ++counters.b2_identity_failures;
      }
      if (square_log_concave && radius >= 1 && target > radius) {
        const int smaller = radius - 1;
        const int larger = target - 1;
        for (int contraction = 0; contraction <= smaller; ++contraction) {
          Integer column_sum = 0;
          for (int row = 0; row <= smaller - contraction; ++row) {
            column_sum += b2_coefficient(
                square,
                smaller + larger - row - 2 * contraction, row);
          }
          ++counters.square_log_concave_b2_columns;
          if (column_sum < 0) {
            ++counters.square_log_concave_b2_column_failures;
            if (counters.first_square_log_concave_b2_column_root.empty()) {
              counters.first_square_log_concave_b2_column_root = root;
              counters.first_square_log_concave_b2_column_square = square;
              counters.first_square_log_concave_b2_column_radius = radius;
              counters.first_square_log_concave_b2_column_target = target;
              counters.first_square_log_concave_b2_column = contraction;
              counters.first_square_log_concave_b2_column_value = column_sum;
            }
          }
        }
      }
      if (margin < 0) {
        ++counters.failures;
        if (counters.first_root.empty()) {
          counters.first_root = root;
          counters.first_square = square;
          counters.first_radius = radius;
          counters.first_target = target;
          counters.first_value = margin;
        }
      }
    }
  }
}

void enumerate(const int length, const int maximum_coefficient,
               const int index, std::vector<int>& root, Counters& counters) {
  if (index == length) {
    if (root.back() != 0) {
      inspect(root, counters);
    }
    return;
  }
  const int minimum = index + 1 == length ? 1 : 0;
  for (int coefficient = minimum; coefficient <= maximum_coefficient;
       ++coefficient) {
    root[static_cast<std::size_t>(index)] = coefficient;
    enumerate(length, maximum_coefficient, index + 1, root, counters);
  }
}

void enumerate_log_concave(const int length, const int maximum_coefficient,
                           const int maximum_shift, const int index,
                           std::vector<int>& core, Counters& counters) {
  if (index == length) {
    for (int shift = 0; shift <= maximum_shift; ++shift) {
      std::vector<int> root(static_cast<std::size_t>(shift), 0);
      root.insert(root.end(), core.begin(), core.end());
      inspect(root, counters);
    }
    return;
  }
  for (int coefficient = 1; coefficient <= maximum_coefficient;
       ++coefficient) {
    if (index >= 2) {
      const std::int64_t previous =
          core[static_cast<std::size_t>(index - 2)];
      const std::int64_t middle =
          core[static_cast<std::size_t>(index - 1)];
      if (middle * middle <
          previous * static_cast<std::int64_t>(coefficient)) {
        continue;
      }
    }
    core[static_cast<std::size_t>(index)] = coefficient;
    enumerate_log_concave(length, maximum_coefficient, maximum_shift,
                          index + 1, core, counters);
  }
}

void enumerate_interval_support(const int length,
                                const int maximum_coefficient,
                                const int maximum_shift, const int index,
                                std::vector<int>& core,
                                Counters& counters) {
  if (index == length) {
    for (int shift = 0; shift <= maximum_shift; ++shift) {
      std::vector<int> root(static_cast<std::size_t>(shift), 0);
      root.insert(root.end(), core.begin(), core.end());
      inspect(root, counters);
    }
    return;
  }
  for (int coefficient = 1; coefficient <= maximum_coefficient;
       ++coefficient) {
    core[static_cast<std::size_t>(index)] = coefficient;
    enumerate_interval_support(length, maximum_coefficient, maximum_shift,
                               index + 1, core, counters);
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const int maximum_length =
        argc >= 2 ? parse_positive(argv[1], "maximum_length") : 6;
    const int maximum_coefficient =
        argc >= 3 ? parse_positive(argv[2], "maximum_coefficient") : 5;
    const std::string mode =
        argc >= 4 ? std::string(argv[3]) : std::string("--all");
    const bool log_concave_only = mode == "--log-concave-only";
    const bool interval_support_only = mode == "--interval-support-only";
    const bool scaled_gap_family = mode == "--scaled-gap-family";
    const bool replay_square_lc_radial =
        mode == "--replay-square-log-concave-radial";
    if (mode != "--all" && !log_concave_only &&
        !interval_support_only && !scaled_gap_family &&
        !replay_square_lc_radial) {
      throw std::invalid_argument(
          "third argument must be --log-concave-only or "
          "--interval-support-only or --scaled-gap-family or "
          "--replay-square-log-concave-radial");
    }
    const int maximum_shift =
        argc >= 5 ? parse_positive(argv[4], "maximum_shift") : 5;
    if (argc > 5) {
      throw std::invalid_argument(
          "usage: probe_su2_ordinary_arbitrary_root_current "
          "[maximum_length] [maximum_coefficient] "
          "[(--log-concave-only|--interval-support-only|"
          "--scaled-gap-family|"
          "--replay-square-log-concave-radial) "
          "[maximum_shift]]");
    }

    if (replay_square_lc_radial) {
      replay_square_log_concave_radial_obstruction();
      return EXIT_SUCCESS;
    }

    Counters counters;
    if (scaled_gap_family) {
      for (int scale = 1; scale <= maximum_coefficient; ++scale) {
        inspect({scale, 2 * scale, 2 * scale, 1, scale}, counters);
      }
    } else {
      for (int length = 1; length <= maximum_length; ++length) {
        std::vector<int> root(static_cast<std::size_t>(length));
        if (log_concave_only) {
          enumerate_log_concave(length, maximum_coefficient, maximum_shift,
                                0, root, counters);
        } else if (interval_support_only) {
          enumerate_interval_support(length, maximum_coefficient,
                                     maximum_shift, 0, root, counters);
        } else {
          enumerate(length, maximum_coefficient, 0, root, counters);
        }
      }
    }

    std::cout << "SU2_ORDINARY_ARBITRARY_ROOT_CURRENT"
              << " maximum_length=" << maximum_length
              << " maximum_coefficient=" << maximum_coefficient
              << " log_concave_only=" << (log_concave_only ? 1 : 0)
              << " interval_support_only="
              << (interval_support_only ? 1 : 0)
              << " scaled_gap_family=" << (scaled_gap_family ? 1 : 0)
              << " maximum_shift=" << maximum_shift
              << " roots=" << counters.roots
              << " currents=" << counters.currents
              << " failures=" << counters.failures
              << " log_concave_roots=" << counters.log_concave_roots
              << " log_concave_currents="
              << counters.log_concave_currents
              << " log_concave_failures="
              << counters.log_concave_failures
              << " interval_support_roots="
              << counters.interval_support_roots
              << " interval_support_currents="
              << counters.interval_support_currents
              << " interval_support_failures="
              << counters.interval_support_failures
              << " negative_b2_coefficients="
              << counters.negative_b2_coefficients
              << " square_log_concave_roots="
              << counters.square_log_concave_roots
              << " square_log_concave_b2_columns="
              << counters.square_log_concave_b2_columns
              << " square_log_concave_b2_column_failures="
              << counters.square_log_concave_b2_column_failures
              << " b2_identity_failures="
              << counters.b2_identity_failures << '\n';
    if (!counters.first_root.empty()) {
      std::cout << "first_root=" << render(counters.first_root)
                << " first_square=[";
      for (std::size_t index = 0; index < counters.first_square.size();
           ++index) {
        if (index != 0U) {
          std::cout << ',';
        }
        std::cout << counters.first_square[index];
      }
      std::cout << "] first_radius=" << counters.first_radius
                << " first_target=" << counters.first_target
                << " first_value=" << counters.first_value << '\n';
    }
    if (!counters.first_log_concave_root.empty()) {
      std::cout << "first_log_concave_root="
                << render(counters.first_log_concave_root)
                << " first_log_concave_radius="
                << counters.first_log_concave_radius
                << " first_log_concave_target="
                << counters.first_log_concave_target
                << " first_log_concave_value="
                << counters.first_log_concave_value << '\n';
    }
    if (!counters.first_interval_support_root.empty()) {
      std::cout << "first_interval_support_root="
                << render(counters.first_interval_support_root)
                << " first_interval_support_radius="
                << counters.first_interval_support_radius
                << " first_interval_support_target="
                << counters.first_interval_support_target
                << " first_interval_support_value="
                << counters.first_interval_support_value << '\n';
    }
    if (!counters.first_square_log_concave_b2_column_root.empty()) {
      std::cout << "first_square_log_concave_b2_column_root="
                << render(counters.first_square_log_concave_b2_column_root)
                << " first_square_log_concave_b2_column_square=[";
      for (std::size_t index = 0;
           index < counters.first_square_log_concave_b2_column_square.size();
           ++index) {
        if (index != 0U) {
          std::cout << ',';
        }
        std::cout
            << counters.first_square_log_concave_b2_column_square[index];
      }
      std::cout << ']'
                << " first_square_log_concave_b2_column_radius="
                << counters.first_square_log_concave_b2_column_radius
                << " first_square_log_concave_b2_column_target="
                << counters.first_square_log_concave_b2_column_target
                << " first_square_log_concave_b2_column="
                << counters.first_square_log_concave_b2_column
                << " first_square_log_concave_b2_column_value="
                << counters.first_square_log_concave_b2_column_value
                << '\n';
    }
    return counters.b2_identity_failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
