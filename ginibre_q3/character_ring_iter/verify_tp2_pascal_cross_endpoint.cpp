#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
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

Vector laplacian_update(const Vector& values) {
  const int modulus = static_cast<int>(values.size());
  Vector result(values.size());
  for (int index = 0; index < modulus; ++index) {
    result[static_cast<std::size_t>(index)] =
        2 * values[static_cast<std::size_t>(index)] -
        values[static_cast<std::size_t>(residue(index - 1, modulus))] -
        values[static_cast<std::size_t>(residue(index + 1, modulus))];
  }
  return result;
}

Vector plus_update(const Vector& values) {
  const int modulus = static_cast<int>(values.size());
  Vector result(values.size());
  for (int index = 0; index < modulus; ++index) {
    result[static_cast<std::size_t>(index)] =
        2 * values[static_cast<std::size_t>(index)] +
        values[static_cast<std::size_t>(residue(index - 1, modulus))] +
        values[static_cast<std::size_t>(residue(index + 1, modulus))];
  }
  return result;
}

Vector endpoint_vector(int modulus, int minus_power, int plus_power) {
  Vector result(static_cast<std::size_t>(modulus));
  result[0] = 1;
  for (int power = 0; power < minus_power; ++power) {
    result = laplacian_update(result);
  }
  for (int power = 0; power < plus_power; ++power) {
    result = plus_update(result);
  }
  return result;
}

State wedge_of(const Vector& left, const Vector& right) {
  if (left.size() != right.size()) {
    throw std::runtime_error("incompatible wedge vectors");
  }
  State result;
  const int modulus = static_cast<int>(left.size());
  for (int first = 0; first < modulus; ++first) {
    for (int second = first + 1; second < modulus; ++second) {
      const cpp_int value =
          left[static_cast<std::size_t>(first)] *
              right[static_cast<std::size_t>(second)] -
          left[static_cast<std::size_t>(second)] *
              right[static_cast<std::size_t>(first)];
      add_wedge(result, first, second, value);
    }
  }
  return result;
}

State compound_update(const State& state, int modulus, bool laplacian) {
  const std::vector<std::pair<int, int>> shifts = laplacian
      ? std::vector<std::pair<int, int>>{{0, 2}, {-1, -1}, {1, -1}}
      : std::vector<std::pair<int, int>>{{0, 2}, {-1, 1}, {1, 1}};
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

State compound_endpoint(State state, int modulus, int minus_power,
                        int plus_power) {
  for (int power = 0; power < minus_power; ++power) {
    state = compound_update(state, modulus, true);
  }
  for (int power = 0; power < plus_power; ++power) {
    state = compound_update(state, modulus, false);
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

State scale_state(State state, const cpp_int& scalar) {
  for (auto iterator = state.begin(); iterator != state.end();) {
    iterator->second *= scalar;
    if (iterator->second == 0) {
      iterator = state.erase(iterator);
    } else {
      ++iterator;
    }
  }
  return state;
}

State translate_state(const State& state, int modulus, int shift) {
  State result;
  for (const auto& [pair, value] : state) {
    add_wedge(result, residue(pair.first + shift, modulus),
              residue(pair.second + shift, modulus), value);
  }
  return result;
}

State halve_state(State state) {
  for (auto& [pair, value] : state) {
    static_cast<void>(pair);
    if (value % 2 != 0) {
      throw std::runtime_error("nonintegral half-state");
    }
    value /= 2;
  }
  return state;
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

cpp_int pascal_cross(const Vector& left, const Vector& right) {
  const auto at = [](const Vector& values, int index) -> const cpp_int& {
    return values[static_cast<std::size_t>(index)];
  };
  return at(left, 2) * at(right, 4) +
         at(left, 4) * at(right, 2) +
         2 * at(left, 2) * at(right, 2) -
         2 * at(left, 3) * at(right, 3) -
         at(left, 1) * at(right, 3) -
         at(left, 3) * at(right, 1);
}

}  // namespace

int main(int argc, char** argv) {
  const int maximum_rank = argc > 1 ? std::atoi(argv[1]) : 20;
  const int maximum_minus_power = argc > 2 ? std::atoi(argv[2]) : 40;
  const int maximum_plus_power = argc > 3 ? std::atoi(argv[3]) : 40;
  const bool all_pairs = argc > 4 && std::atoi(argv[4]) != 0;
  if (maximum_rank < 2 || maximum_minus_power < 4 ||
      maximum_plus_power < 2) {
    std::cerr << "usage: verify_tp2_pascal_cross_endpoint "
                 "[maximum-rank>=2] [maximum-minus-power>=4] "
                 "[maximum-plus-power>=2] [all-pairs=0|1]\n";
    return 2;
  }

  std::size_t systems = 0;
  for (int rank = 2; rank <= maximum_rank; ++rank) {
    const int modulus = 2 * rank + 1;
    Vector identity(static_cast<std::size_t>(modulus));
    identity[0] = 1;
    const Vector first_laplacian = laplacian_update(identity);
    const Vector second_laplacian = laplacian_update(first_laplacian);
    Vector reflected_identity = identity;
    for (int index = 0; index < modulus; ++index) {
      reflected_identity[static_cast<std::size_t>(index)] =
          2 * identity[static_cast<std::size_t>(index)] -
          first_laplacian[static_cast<std::size_t>(index)];
    }
    const State cross_boundary = wedge_of(reflected_identity,
                                          second_laplacian);
    State target_wedge;
    add_wedge(target_wedge, 2, 3, 1);
    const State cross_operator = add_states(
        add_states(scale_state(target_wedge, 16),
                   scale_state(compound_update(target_wedge, modulus, true),
                               -1)),
        scale_state(compound_update(target_wedge, modulus, false), -1)
    );
    const State centre_difference = add_states(
        translate_state(cross_operator, modulus, -3),
        scale_state(translate_state(cross_operator, modulus, -2), -1)
    );
    const State centre_boundary = halve_state(centre_difference);
    if (cross_boundary != centre_boundary) {
      std::cout << "TP2_PASCAL_CROSS_ENDPOINT result=CENTRE_IDENTITY_FAIL"
                << " rank=" << rank << '\n';
      return 1;
    }
    for (int plus_power = 2; plus_power <= maximum_plus_power;
         ++plus_power) {
      const int first_minus_power = all_pairs ? 2 : plus_power + 2;
      for (int minus_power = first_minus_power;
           minus_power <= maximum_minus_power; ++minus_power) {
        const Vector left = endpoint_vector(modulus, minus_power + 1,
                                            plus_power);
        const Vector right = endpoint_vector(modulus, minus_power,
                                             plus_power + 1);
        const cpp_int direct_cross = pascal_cross(left, right);
        const State propagated = compound_endpoint(
            cross_boundary, modulus, minus_power, plus_power
        );
        const cpp_int endpoint_cross = -2 * coefficient(propagated, 2, 3);
        ++systems;
        if (direct_cross != endpoint_cross) {
          std::cout << "TP2_PASCAL_CROSS_ENDPOINT result=IDENTITY_FAIL"
                    << " rank=" << rank
                    << " minus_power=" << minus_power
                    << " plus_power=" << plus_power
                    << " direct=" << direct_cross
                    << " endpoint=" << endpoint_cross << '\n';
          return 1;
        }
        if (direct_cross < 0) {
          std::cout << "TP2_PASCAL_CROSS_ENDPOINT result=TARGET_FAIL"
                    << " rank=" << rank
                    << " minus_power=" << minus_power
                    << " plus_power=" << plus_power
                    << " value=" << direct_cross << '\n';
          return 1;
        }
      }
    }
  }
  std::cout << "TP2_PASCAL_CROSS_ENDPOINT maximum_rank=" << maximum_rank
            << " maximum_minus_power=" << maximum_minus_power
            << " maximum_plus_power=" << maximum_plus_power
            << " all_pairs=" << (all_pairs ? 1 : 0)
            << " systems=" << systems << " result=PASS\n";
  return 0;
}
