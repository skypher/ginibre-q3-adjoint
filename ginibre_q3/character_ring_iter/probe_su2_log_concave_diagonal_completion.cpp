#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::int64_t value(const std::vector<int>& profile, const int index) {
  return index >= 0 && index < static_cast<int>(profile.size())
             ? profile[static_cast<std::size_t>(index)]
             : 0;
}

std::int64_t current(const std::vector<int>& profile, const int radius,
                     const int target) {
  std::int64_t interval = 0;
  for (int label = std::abs(radius - target);
       label <= radius + target; ++label) {
    interval += value(profile, label);
  }
  return value(profile, 0) * interval -
         value(profile, radius) * value(profile, target);
}

int parse_positive(const char* text, const std::string& name) {
  const std::string value_text{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value_text, &consumed, 10);
  if (consumed != value_text.size() || parsed <= 0) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  return static_cast<int>(parsed);
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
  std::uint64_t log_concave_profiles = 0;
  std::uint64_t diagonal_admissible_profiles = 0;
  std::uint64_t off_diagonal_currents = 0;
  std::uint64_t failures = 0;
  std::uint64_t schur_escape_failures = 0;
  std::uint64_t boundary_admissible_profiles = 0;
  std::uint64_t boundary_completion_currents = 0;
  std::uint64_t boundary_completion_failures = 0;
  std::uint64_t schur_boundary_admissible_profiles = 0;
  std::uint64_t schur_boundary_currents = 0;
  std::uint64_t schur_boundary_failures = 0;
  std::vector<int> first_profile;
  int first_radius = -1;
  int first_target = -1;
  std::int64_t first_value = 0;
  std::vector<int> first_schur_escape_profile;
  int first_schur_escape_radius = -1;
  int first_schur_escape_target = -1;
  std::int64_t first_schur_escape_cross = 0;
  std::int64_t first_schur_escape_determinant = 0;
  std::vector<int> first_boundary_completion_profile;
  int first_boundary_completion_radius = -1;
  int first_boundary_completion_target = -1;
  std::int64_t first_boundary_completion_value = 0;
  std::vector<int> first_schur_boundary_profile;
  int first_schur_boundary_radius = -1;
  int first_schur_boundary_target = -1;
  std::int64_t first_schur_boundary_value = 0;
};

void inspect_core(const std::vector<int>& core, const int maximum_shift,
                  Counters& counters) {
  for (int shift = 0; shift <= maximum_shift; ++shift) {
    std::vector<int> profile(static_cast<std::size_t>(shift), 0);
    profile.insert(profile.end(), core.begin(), core.end());
    ++counters.log_concave_profiles;
    const int support = static_cast<int>(profile.size()) - 1;
    bool diagonal_admissible = true;
    for (int radius = 1; radius <= support; ++radius) {
      if (current(profile, radius, radius) < 0) {
        diagonal_admissible = false;
        break;
      }
    }
    if (!diagonal_admissible) {
      continue;
    }
    ++counters.diagonal_admissible_profiles;
    bool boundary_admissible = true;
    for (int radius = 1; radius < support; ++radius) {
      if (current(profile, radius, support) < 0) {
        boundary_admissible = false;
        break;
      }
    }
    if (boundary_admissible) {
      ++counters.boundary_admissible_profiles;
      for (int radius = 1; radius <= support; ++radius) {
        for (int target = radius + 1; target < support; ++target) {
          const std::int64_t margin = current(profile, radius, target);
          ++counters.boundary_completion_currents;
          if (margin < 0) {
            ++counters.boundary_completion_failures;
            if (counters.first_boundary_completion_profile.empty()) {
              counters.first_boundary_completion_profile = profile;
              counters.first_boundary_completion_radius = radius;
              counters.first_boundary_completion_target = target;
              counters.first_boundary_completion_value = margin;
            }
          }
        }
      }

      bool schur_admissible = true;
      for (int radius = 1; radius <= support && schur_admissible; ++radius) {
        const std::int64_t radius_diagonal =
            current(profile, radius, radius);
        for (int target = radius + 1; target <= support; ++target) {
          const std::int64_t cross = current(profile, radius, target);
          const std::int64_t target_diagonal =
              current(profile, target, target);
          if (radius_diagonal * target_diagonal - cross * cross < 0) {
            schur_admissible = false;
            break;
          }
        }
      }
      if (schur_admissible) {
        ++counters.schur_boundary_admissible_profiles;
        for (int radius = 1; radius <= support; ++radius) {
          for (int target = radius + 1; target < support; ++target) {
            const std::int64_t margin = current(profile, radius, target);
            ++counters.schur_boundary_currents;
            if (margin < 0) {
              ++counters.schur_boundary_failures;
              if (counters.first_schur_boundary_profile.empty()) {
                counters.first_schur_boundary_profile = profile;
                counters.first_schur_boundary_radius = radius;
                counters.first_schur_boundary_target = target;
                counters.first_schur_boundary_value = margin;
              }
            }
          }
        }
      }
    }
    for (int radius = 1; radius <= support; ++radius) {
      for (int target = radius + 1; target <= support; ++target) {
        const std::int64_t margin = current(profile, radius, target);
        ++counters.off_diagonal_currents;
        if (margin < 0) {
          ++counters.failures;
          if (counters.first_profile.empty()) {
            counters.first_profile = profile;
            counters.first_radius = radius;
            counters.first_target = target;
            counters.first_value = margin;
          }
          const std::int64_t radius_diagonal =
              current(profile, radius, radius);
          const std::int64_t target_diagonal =
              current(profile, target, target);
          const std::int64_t schur_determinant =
              radius_diagonal * target_diagonal - margin * margin;
          if (schur_determinant >= 0) {
            ++counters.schur_escape_failures;
            if (counters.first_schur_escape_profile.empty()) {
              counters.first_schur_escape_profile = profile;
              counters.first_schur_escape_radius = radius;
              counters.first_schur_escape_target = target;
              counters.first_schur_escape_cross = margin;
              counters.first_schur_escape_determinant = schur_determinant;
            }
          }
        }
      }
    }
  }
}

void enumerate(const int length, const int maximum_coefficient,
               const int maximum_shift, const int index,
               std::vector<int>& core, Counters& counters) {
  if (index == length) {
    inspect_core(core, maximum_shift, counters);
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
    enumerate(length, maximum_coefficient, maximum_shift, index + 1, core,
              counters);
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const int maximum_length =
        argc >= 2 ? parse_positive(argv[1], "maximum_length") : 7;
    const int maximum_coefficient =
        argc >= 3 ? parse_positive(argv[2], "maximum_coefficient") : 8;
    const int maximum_shift =
        argc >= 4 ? parse_positive(argv[3], "maximum_shift") : 4;
    if (argc > 4) {
      throw std::invalid_argument(
          "usage: probe_su2_log_concave_diagonal_completion "
          "[maximum_length] [maximum_coefficient] [maximum_shift]");
    }

    Counters counters;
    for (int length = 1; length <= maximum_length; ++length) {
      std::vector<int> core(static_cast<std::size_t>(length));
      enumerate(length, maximum_coefficient, maximum_shift, 0, core,
                counters);
    }

    std::cout << "SU2_LOG_CONCAVE_DIAGONAL_COMPLETION"
              << " maximum_length=" << maximum_length
              << " maximum_coefficient=" << maximum_coefficient
              << " maximum_shift=" << maximum_shift
              << " log_concave_profiles=" << counters.log_concave_profiles
              << " diagonal_admissible_profiles="
              << counters.diagonal_admissible_profiles
              << " off_diagonal_currents=" << counters.off_diagonal_currents
              << " failures=" << counters.failures
              << " schur_escape_failures="
              << counters.schur_escape_failures
              << " boundary_admissible_profiles="
              << counters.boundary_admissible_profiles
              << " boundary_completion_currents="
              << counters.boundary_completion_currents
              << " boundary_completion_failures="
              << counters.boundary_completion_failures
              << " schur_boundary_admissible_profiles="
              << counters.schur_boundary_admissible_profiles
              << " schur_boundary_currents="
              << counters.schur_boundary_currents
              << " schur_boundary_failures="
              << counters.schur_boundary_failures << '\n';
    if (!counters.first_profile.empty()) {
      std::cout << "first_profile=" << render(counters.first_profile)
                << " first_radius=" << counters.first_radius
                << " first_target=" << counters.first_target
                << " first_value=" << counters.first_value << '\n';
    }
    if (!counters.first_schur_escape_profile.empty()) {
      std::cout << "first_schur_escape_profile="
                << render(counters.first_schur_escape_profile)
                << " first_schur_escape_radius="
                << counters.first_schur_escape_radius
                << " first_schur_escape_target="
                << counters.first_schur_escape_target
                << " first_schur_escape_cross="
                << counters.first_schur_escape_cross
                << " first_schur_escape_determinant="
                << counters.first_schur_escape_determinant << '\n';
    }
    if (!counters.first_boundary_completion_profile.empty()) {
      std::cout << "first_boundary_completion_profile="
                << render(counters.first_boundary_completion_profile)
                << " first_boundary_completion_radius="
                << counters.first_boundary_completion_radius
                << " first_boundary_completion_target="
                << counters.first_boundary_completion_target
                << " first_boundary_completion_value="
                << counters.first_boundary_completion_value << '\n';
    }
    if (!counters.first_schur_boundary_profile.empty()) {
      std::cout << "first_schur_boundary_profile="
                << render(counters.first_schur_boundary_profile)
                << " first_schur_boundary_radius="
                << counters.first_schur_boundary_radius
                << " first_schur_boundary_target="
                << counters.first_schur_boundary_target
                << " first_schur_boundary_value="
                << counters.first_schur_boundary_value << '\n';
      const std::vector<int>& profile =
          counters.first_schur_boundary_profile;
      const int support = static_cast<int>(profile.size()) - 1;
      std::int64_t minimum_diagonal = current(profile, 1, 1);
      std::int64_t minimum_boundary = current(profile, 1, support);
      std::int64_t minimum_schur_determinant =
          current(profile, 1, 1) * current(profile, 2, 2) -
          current(profile, 1, 2) * current(profile, 1, 2);
      for (int radius = 1; radius <= support; ++radius) {
        minimum_diagonal =
            std::min(minimum_diagonal, current(profile, radius, radius));
        if (radius < support) {
          minimum_boundary =
              std::min(minimum_boundary,
                       current(profile, radius, support));
        }
        for (int target = radius + 1; target <= support; ++target) {
          const std::int64_t cross = current(profile, radius, target);
          minimum_schur_determinant =
              std::min(minimum_schur_determinant,
                       current(profile, radius, radius) *
                               current(profile, target, target) -
                           cross * cross);
        }
      }
      std::int64_t pi_over_two_value = 0;
      for (int label = 0; label <= support; ++label) {
        pi_over_two_value +=
            (label % 2 == 0 ? 1 : -1) * value(profile, label);
      }
      std::cout << "first_schur_boundary_minimum_diagonal="
                << minimum_diagonal
                << " first_schur_boundary_minimum_boundary="
                << minimum_boundary
                << " first_schur_boundary_minimum_schur_determinant="
                << minimum_schur_determinant
                << " first_schur_boundary_character_at_pi_over_two="
                << pi_over_two_value << '\n';
    }
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
