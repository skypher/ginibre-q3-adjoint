#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Pair = std::pair<int, int>;
using State = std::map<Pair, cpp_int>;
using Vector = std::vector<cpp_int>;

int residue(int value, int modulus) {
  value %= modulus;
  return value < 0 ? value + modulus : value;
}

void add_wedge(State& state, int first, int second, const cpp_int& value) {
  if (first == second || value == 0) {
    return;
  }
  if (first < second) {
    state[{first, second}] += value;
  } else {
    state[{second, first}] -= value;
  }
}

State compound_laplacian_update(const State& state, int modulus, int step) {
  const std::vector<std::pair<int, int>> shifts{{0, 2}, {-step, -1},
                                                {step, -1}};
  State result;
  for (const auto& [pair, value] : state) {
    for (const auto& [first_shift, first_weight] : shifts) {
      for (const auto& [second_shift, second_weight] : shifts) {
        add_wedge(result, residue(pair.first + first_shift, modulus),
                  residue(pair.second + second_shift, modulus),
                  value * first_weight * second_weight);
      }
    }
  }
  for (auto iterator = result.begin(); iterator != result.end();) {
    if (iterator->second == 0) {
      iterator = result.erase(iterator);
    } else {
      ++iterator;
    }
  }
  return result;
}

State compound_plus_update(const State& state, int modulus) {
  const std::vector<std::pair<int, int>> shifts{{0, 2}, {-1, 1}, {1, 1}};
  State result;
  for (const auto& [pair, value] : state) {
    for (const auto& [first_shift, first_weight] : shifts) {
      for (const auto& [second_shift, second_weight] : shifts) {
        add_wedge(result, residue(pair.first + first_shift, modulus),
                  residue(pair.second + second_shift, modulus),
                  value * first_weight * second_weight);
      }
    }
  }
  for (auto iterator = result.begin(); iterator != result.end();) {
    if (iterator->second == 0) {
      iterator = result.erase(iterator);
    } else {
      ++iterator;
    }
  }
  return result;
}

State compound_binomial_update(const State& state, int modulus, int sign) {
  const std::vector<std::pair<int, int>> shifts{{0, 1}, {-1, sign}};
  State result;
  for (const auto& [pair, value] : state) {
    for (const auto& [first_shift, first_weight] : shifts) {
      for (const auto& [second_shift, second_weight] : shifts) {
        add_wedge(result, residue(pair.first + first_shift, modulus),
                  residue(pair.second + second_shift, modulus),
                  value * first_weight * second_weight);
      }
    }
  }
  for (auto iterator = result.begin(); iterator != result.end();) {
    if (iterator->second == 0) {
      iterator = result.erase(iterator);
    } else {
      ++iterator;
    }
  }
  return result;
}

State endpoint_power(State state, int modulus, int minus_power,
                     int plus_power) {
  for (int power = 0; power < minus_power; ++power) {
    state = compound_laplacian_update(state, modulus, 1);
  }
  for (int power = 0; power < plus_power; ++power) {
    state = compound_plus_update(state, modulus);
  }
  return state;
}

State add_states(State left, const State& right) {
  for (const auto& [pair, value] : right) {
    left[pair] += value;
  }
  for (auto iterator = left.begin(); iterator != left.end();) {
    if (iterator->second == 0) {
      iterator = left.erase(iterator);
    } else {
      ++iterator;
    }
  }
  return left;
}

cpp_int coefficient(const State& state, int first, int second) {
  int sign = 1;
  if (first > second) {
    std::swap(first, second);
    sign = -1;
  }
  const auto found = state.find({first, second});
  return found == state.end() ? cpp_int(0) : sign * found->second;
}

cpp_int inner_product(const State& left, const State& right) {
  cpp_int result = 0;
  for (const auto& [pair, value] : left) {
    const auto found = right.find(pair);
    if (found != right.end()) {
      result += value * found->second;
    }
  }
  return result;
}

void print_state(const char* label, const State& state) {
  std::cout << "TP2_COMPOUND " << label << '=';
  for (const auto& [pair, value] : state) {
    std::cout << ' ' << value << "*(" << pair.first << ',' << pair.second
              << ')';
  }
  std::cout << '\n';
}

Vector laplacian_update(const Vector& vector, int step) {
  const int modulus = static_cast<int>(vector.size());
  Vector result(vector.size());
  for (int index = 0; index < modulus; ++index) {
    result[static_cast<std::size_t>(index)] =
        2 * vector[static_cast<std::size_t>(index)] -
        vector[static_cast<std::size_t>(residue(index - step, modulus))] -
        vector[static_cast<std::size_t>(residue(index + step, modulus))];
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  const int maximum_rank = argc > 1 ? std::atoi(argv[1]) : 20;
  const int maximum_difference = argc > 2 ? std::atoi(argv[2]) : 30;
  const int maximum_half_power = argc > 3 ? std::atoi(argv[3]) : 30;
  const bool dump_base = argc > 4 && std::atoi(argv[4]) != 0;
  if (maximum_rank < 5 || maximum_difference < 2 ||
      maximum_half_power < 2) {
    std::cerr << "usage: probe_tp2_multiplicative_compound "
                 "[maximum-rank>=5] [maximum-difference>=2] "
                 "[maximum-half-power>=2] [dump-base=0|1]\n";
    return 2;
  }

  std::size_t systems = 0;
  std::size_t sign_conflicts = 0;
  std::size_t negative_targets = 0;
  std::size_t square_source_checks = 0;
  std::size_t identity_checks = 0;
  std::size_t positive_cauchy_binet_terms = 0;
  std::size_t positive_cauchy_binet_separation_groups = 0;
  std::size_t positive_cauchy_binet_center_groups = 0;
  bool reported_cauchy_binet = false;
  bool reported_cauchy_binet_group = false;
  bool reported_cauchy_binet_center = false;
  bool reported_conflict = false;
  for (int rank = 5; rank <= maximum_rank; ++rank) {
    const int modulus = 2 * rank + 1;
    State source;
    add_wedge(source, 0, residue(-1, modulus), -1);
    add_wedge(source, 0, 1, -1);

    State target_wedge;
    add_wedge(target_wedge, 2, 3, 1);
    State transformed_source = endpoint_power(source, modulus, 2, 1);
    State transformed_target = endpoint_power(target_wedge, modulus, 2, 1);
    if (dump_base && rank == 5) {
      print_state("base_source", transformed_source);
      print_state("base_target", transformed_target);
    }
    ++square_source_checks;
    const State square_defect =
        add_states(transformed_source, transformed_target);
    if (!square_defect.empty() && rank == 5) {
      std::cout << "TP2_COMPOUND symmetric_square_intertwiner=false\n";
    }
    if (rank == 5) {
      int intertwiner_splits = 0;
      for (int minus_left = 0; minus_left <= 4; ++minus_left) {
        for (int plus_left = 0; plus_left <= 2; ++plus_left) {
          const State left = endpoint_power(target_wedge, modulus, minus_left,
                                            plus_left);
          const State right = endpoint_power(source, modulus, 4 - minus_left,
                                             2 - plus_left);
          if (add_states(left, right).empty()) {
            ++intertwiner_splits;
            std::cout << "TP2_COMPOUND intertwiner_split minus_left="
                      << minus_left << " plus_left=" << plus_left << '\n';
          }
        }
      }
      if (intertwiner_splits == 0) {
        std::cout << "TP2_COMPOUND asymmetric_intertwiner=false\n";
      }
      int polynomial_intertwiners = 0;
      for (int minus_power = 0; minus_power <= 6; ++minus_power) {
        for (int plus_power = 0; plus_power <= 6; ++plus_power) {
          const State forward = endpoint_power(
              transformed_target, modulus, minus_power, plus_power);
          if (add_states(forward, transformed_source).empty()) {
            ++polynomial_intertwiners;
            std::cout << "TP2_COMPOUND polynomial_intertwiner direction=UV"
                      << " minus_power=" << minus_power
                      << " plus_power=" << plus_power << '\n';
          }
          const State reverse = endpoint_power(
              transformed_source, modulus, minus_power, plus_power);
          if (add_states(reverse, transformed_target).empty()) {
            ++polynomial_intertwiners;
            std::cout << "TP2_COMPOUND polynomial_intertwiner direction=VU"
                      << " minus_power=" << minus_power
                      << " plus_power=" << plus_power << '\n';
          }
        }
      }
      if (polynomial_intertwiners == 0) {
        std::cout << "TP2_COMPOUND polynomial_intertwiner=false\n";
      }
      std::size_t nonorthogonal_defect_moments = 0;
      for (int difference = 0; difference <= maximum_difference;
           ++difference) {
        State defect_moment =
            endpoint_power(square_defect, modulus, difference, 0);
        for (int half_power = 0; half_power <= maximum_half_power;
             ++half_power) {
          if (half_power > 0) {
            defect_moment =
                compound_laplacian_update(defect_moment, modulus, 2);
          }
          if (inner_product(transformed_target, defect_moment) != 0) {
            ++nonorthogonal_defect_moments;
          }
        }
      }
      std::cout << "TP2_COMPOUND nonorthogonal_defect_moments="
                << nonorthogonal_defect_moments << '\n';
    }

    std::vector<State> difference_states(
        static_cast<std::size_t>(maximum_difference + 1));
    difference_states[0] = source;
    for (int difference = 1; difference <= maximum_difference; ++difference) {
      difference_states[static_cast<std::size_t>(difference)] =
          compound_laplacian_update(
              difference_states[static_cast<std::size_t>(difference - 1)],
              modulus, 1);
    }

    std::map<Pair, int> coordinate_sign;
    for (int difference = 2; difference <= maximum_difference; ++difference) {
      State state = difference_states[static_cast<std::size_t>(difference)];
      for (int half_power = 1; half_power <= maximum_half_power; ++half_power) {
        state = compound_laplacian_update(state, modulus, 2);
        if (half_power < 2) {
          continue;
        }
        ++systems;
        const cpp_int target = -coefficient(state, 2, 3);
        Vector direct(static_cast<std::size_t>(modulus));
        direct[0] = 1;
        for (int power = 0; power < difference; ++power) {
          direct = laplacian_update(direct, 1);
        }
        for (int power = 0; power < half_power; ++power) {
          direct = laplacian_update(direct, 2);
        }
        const cpp_int direct_target =
            direct[2] * direct[4] + direct[2] * direct[2] -
            direct[3] * direct[3] - direct[1] * direct[3];
        ++identity_checks;
        if (target != direct_target) {
          std::cout << "TP2_COMPOUND result=IDENTITY_FAIL"
                    << " rank=" << rank << " difference=" << difference
                    << " half_power=" << half_power
                    << " compound=" << target
                    << " direct=" << direct_target << '\n';
          return 1;
        }
        State binomial_target = target_wedge;
        State binomial_source = source;
        const int minus_power = difference + half_power;
        for (int power = 0; power < minus_power; ++power) {
          binomial_target =
              compound_binomial_update(binomial_target, modulus, -1);
          binomial_source =
              compound_binomial_update(binomial_source, modulus, -1);
        }
        for (int power = 0; power < half_power; ++power) {
          binomial_target =
              compound_binomial_update(binomial_target, modulus, 1);
          binomial_source =
              compound_binomial_update(binomial_source, modulus, 1);
        }
        const cpp_int cauchy_binet_target =
            -inner_product(binomial_target, binomial_source);
        if (target != cauchy_binet_target) {
          std::cout << "TP2_COMPOUND result=CAUCHY_BINET_FAIL"
                    << " rank=" << rank << " difference=" << difference
                    << " half_power=" << half_power
                    << " compound=" << target
                    << " cauchy_binet=" << cauchy_binet_target << '\n';
          return 1;
        }
        std::map<int, cpp_int> cauchy_binet_separation_groups;
        std::map<int, cpp_int> cauchy_binet_center_groups;
        for (const auto& [pair, value] : binomial_target) {
          const auto found = binomial_source.find(pair);
          if (found == binomial_source.end()) {
            continue;
          }
          const cpp_int product = value * found->second;
          const int forward_separation = pair.second - pair.first;
          const int separation =
              std::min(forward_separation, modulus - forward_separation);
          cauchy_binet_separation_groups[separation] += product;
          cauchy_binet_center_groups[
              residue(pair.first + pair.second, modulus)] += product;
          if (product > 0) {
            ++positive_cauchy_binet_terms;
            if (!reported_cauchy_binet) {
              reported_cauchy_binet = true;
              std::cout << "TP2_COMPOUND first_positive_cauchy_binet_term"
                        << " rank=" << rank << " difference=" << difference
                        << " half_power=" << half_power << " pair=("
                        << pair.first << ',' << pair.second << ")"
                        << " target_minor=" << value
                        << " source_minor=" << found->second << '\n';
            }
          }
        }
        for (const auto& [center, value] : cauchy_binet_center_groups) {
          if (value > 0) {
            ++positive_cauchy_binet_center_groups;
            if (!reported_cauchy_binet_center) {
              reported_cauchy_binet_center = true;
              std::cout << "TP2_COMPOUND first_positive_cauchy_binet_center"
                        << " rank=" << rank
                        << " difference=" << difference
                        << " half_power=" << half_power
                        << " center=" << center << " value=" << value
                        << '\n';
            }
          }
        }
        for (const auto& [separation, value] :
             cauchy_binet_separation_groups) {
          if (value > 0) {
            ++positive_cauchy_binet_separation_groups;
            if (!reported_cauchy_binet_group) {
              reported_cauchy_binet_group = true;
              std::cout
                  << "TP2_COMPOUND first_positive_cauchy_binet_separation"
                  << " rank=" << rank << " difference=" << difference
                  << " half_power=" << half_power
                  << " separation=" << separation << " value=" << value
                  << '\n';
            }
          }
        }
        if (target < 0) {
          ++negative_targets;
          std::cout << "TP2_COMPOUND result=TARGET_FAIL"
                    << " rank=" << rank << " difference=" << difference
                    << " half_power=" << half_power
                    << " target=" << target << '\n';
          return 1;
        }
        for (const auto& [pair, value] : state) {
          const int sign = value > 0 ? 1 : -1;
          const auto [iterator, inserted] =
              coordinate_sign.emplace(pair, sign);
          if (!inserted && iterator->second != sign) {
            ++sign_conflicts;
            if (!reported_conflict) {
              reported_conflict = true;
              std::cout << "TP2_COMPOUND first_sign_conflict"
                        << " rank=" << rank << " difference=" << difference
                        << " half_power=" << half_power << " pair=("
                        << pair.first << ',' << pair.second << ")"
                        << " previous_sign=" << iterator->second
                        << " current_sign=" << sign << '\n';
            }
          }
        }
      }
    }
  }
  std::cout << "TP2_COMPOUND maximum_rank=" << maximum_rank
            << " maximum_difference=" << maximum_difference
            << " maximum_half_power=" << maximum_half_power
            << " systems=" << systems << " sign_conflicts=" << sign_conflicts
            << " negative_targets=" << negative_targets
            << " square_source_checks=" << square_source_checks
            << " identity_checks=" << identity_checks
            << " positive_cauchy_binet_terms="
            << positive_cauchy_binet_terms
            << " positive_cauchy_binet_separation_groups="
            << positive_cauchy_binet_separation_groups
            << " positive_cauchy_binet_center_groups="
            << positive_cauchy_binet_center_groups
            << " result=PASS\n";
  return 0;
}
