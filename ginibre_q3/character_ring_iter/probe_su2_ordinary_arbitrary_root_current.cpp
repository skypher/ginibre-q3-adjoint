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

struct Counters {
  std::uint64_t roots = 0;
  std::uint64_t currents = 0;
  std::uint64_t failures = 0;
  std::uint64_t b2_identity_failures = 0;
  std::uint64_t negative_b2_coefficients = 0;
  std::vector<int> first_root;
  std::vector<Integer> first_square;
  int first_radius = -1;
  int first_target = -1;
  Integer first_value = 0;
};

void inspect(const std::vector<int>& root, Counters& counters) {
  const std::vector<Integer> square = square_character(root);
  ++counters.roots;
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
      if (radius >= 1 &&
          b2_triangle(square, radius, target) != margin) {
        ++counters.b2_identity_failures;
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

}  // namespace

int main(int argc, char** argv) {
  try {
    const int maximum_length =
        argc >= 2 ? parse_positive(argv[1], "maximum_length") : 6;
    const int maximum_coefficient =
        argc >= 3 ? parse_positive(argv[2], "maximum_coefficient") : 5;
    if (argc > 3) {
      throw std::invalid_argument(
          "usage: probe_su2_ordinary_arbitrary_root_current "
          "[maximum_length] [maximum_coefficient]");
    }

    Counters counters;
    for (int length = 1; length <= maximum_length; ++length) {
      std::vector<int> root(static_cast<std::size_t>(length));
      enumerate(length, maximum_coefficient, 0, root, counters);
    }

    std::cout << "SU2_ORDINARY_ARBITRARY_ROOT_CURRENT"
              << " maximum_length=" << maximum_length
              << " maximum_coefficient=" << maximum_coefficient
              << " roots=" << counters.roots
              << " currents=" << counters.currents
              << " failures=" << counters.failures
              << " negative_b2_coefficients="
              << counters.negative_b2_coefficients
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
    return counters.b2_identity_failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
