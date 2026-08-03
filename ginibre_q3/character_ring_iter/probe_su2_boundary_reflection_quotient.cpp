#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int parse_positive(const char* text, const std::string& name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value, &consumed, 10);
  if (consumed != value.size() || parsed <= 0L) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  return static_cast<int>(parsed);
}

bool fuses(const int level, const int factor, const int source,
           const int target) {
  return std::abs(source - factor) <= target
         && target <= std::min(source + factor,
                               2 * level - source - factor)
         && ((source + factor + target) & 1) == 0;
}

int plus_transition(const int level, const int factor, const int from_x,
                    const int from_y, const int to_x, const int to_y) {
  int result = 0;
  if (from_y == to_y && fuses(level, factor, from_x, to_x)) {
    ++result;
  }
  if (from_x == to_x && fuses(level, factor, from_y, to_y)) {
    ++result;
  }
  return result;
}

struct Witness {
  int level = 0;
  int factor = 0;
  int source_x = 0;
  int source_y = 0;
  int target_x = 0;
  int target_y = 0;
  int value = 0;
};

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 3) {
      throw std::invalid_argument(
          "usage: probe_su2_boundary_reflection_quotient MAXIMUM_LEVEL "
          "MAXIMUM_FACTOR");
    }
    const int maximum_level = parse_positive(argv[1], "MAXIMUM_LEVEL");
    const int maximum_factor = parse_positive(argv[2], "MAXIMUM_FACTOR");
    std::uint64_t entries = 0U;
    std::uint64_t negative_entries = 0U;
    bool have_witness = false;
    Witness witness;
    for (int level = 1; level <= maximum_level; ++level) {
      struct Point { int x; int y; };
      std::vector<Point> chamber;
      for (int x = 0; x <= level; ++x) {
        for (int y = 0; y <= level - x; ++y) {
          if (x + y < level) {
            chamber.push_back({x, y});
          }
        }
      }
      const int last_factor = std::min(level - 1, maximum_factor);
      for (int factor = 1; factor <= last_factor; ++factor) {
        for (const Point source : chamber) {
          for (const Point target : chamber) {
            const int quotient =
                plus_transition(level, factor, source.x, source.y,
                                target.x, target.y)
                - plus_transition(level, factor, source.x, source.y,
                                  level - target.y, level - target.x);
            ++entries;
            if (quotient < 0) {
              ++negative_entries;
              if (!have_witness) {
                have_witness = true;
                witness = {level, factor, source.x, source.y,
                           target.x, target.y, quotient};
              }
            }
          }
        }
      }
    }
    std::cout << "SU2_BOUNDARY_REFLECTION_QUOTIENT"
              << " maximum_level=" << maximum_level
              << " maximum_factor=" << maximum_factor
              << " entries=" << entries
              << " negative_entries=" << negative_entries;
    if (have_witness) {
      std::cout << " first_negative={level=" << witness.level
                << ",factor=" << witness.factor
                << ",source=(" << witness.source_x << ',' << witness.source_y
                << "),target=(" << witness.target_x << ',' << witness.target_y
                << "),value=" << witness.value << '}';
    }
    std::cout << " result="
              << (negative_entries == 0U ? "NONNEGATIVE" : "NEGATIVE")
              << '\n';
    return negative_entries == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
