#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

using boost::multiprecision::cpp_int;
using Rational = boost::rational<cpp_int>;

namespace {

struct Interval {
  Rational lower;
  Rational upper;
};

struct PairData {
  Interval endpoint_minus;
  Interval endpoint_plus;
  Interval base_weight;
  Interval signed_base_amplitude;
};

struct ApproximatePairData {
  long double endpoint_minus;
  long double endpoint_plus;
  long double base_amplitude;
};

struct IndexedPairData {
  PairData data;
  PairData coarse_data;
  ApproximatePairData approximate_data;
  int first;
  int second;
};

struct CoverState {
  std::uint64_t low = 0U;
  std::uint64_t high = 0U;

  bool operator==(const CoverState& other) const {
    return low == other.low && high == other.high;
  }
};

struct CoverStateHash {
  std::size_t operator()(const CoverState& state) const {
    const std::uint64_t mixed = state.low ^
        (state.high + 0x9e3779b97f4a7c15ULL + (state.low << 6U) +
         (state.low >> 2U));
    return static_cast<std::size_t>(mixed);
  }
};

struct PartitionMatroidInfo {
  bool valid;
  std::size_t class_count;
};

int parse_positive(const char* text, const char* name) {
  const int value = std::stoi(text);
  if (value <= 0) {
    throw std::runtime_error(std::string(name) + " must be positive");
  }
  return value;
}

Rational minimum(const Rational& left, const Rational& right) {
  return left < right ? left : right;
}

Rational maximum(const Rational& left, const Rational& right) {
  return left > right ? left : right;
}

Interval exact(long long value) {
  return Interval{Rational(value), Rational(value)};
}

Interval add(const Interval& left, const Interval& right) {
  return Interval{left.lower + right.lower, left.upper + right.upper};
}

Interval negate(const Interval& value) {
  return Interval{-value.upper, -value.lower};
}

Interval subtract(const Interval& left, const Interval& right) {
  return add(left, negate(right));
}

Interval multiply(const Interval& left, const Interval& right) {
  const Rational values[4]{
      left.lower * right.lower,
      left.lower * right.upper,
      left.upper * right.lower,
      left.upper * right.upper};
  Rational lower = values[0];
  Rational upper = values[0];
  for (int index = 1; index < 4; ++index) {
    lower = minimum(lower, values[index]);
    upper = maximum(upper, values[index]);
  }
  return Interval{lower, upper};
}

Interval square_nonnegative(const Interval& value) {
  if (value.lower < 0) {
    throw std::runtime_error("square_nonnegative received a negative interval");
  }
  return Interval{value.lower * value.lower, value.upper * value.upper};
}

Interval scale(const Interval& value, const Rational& factor) {
  return multiply(value, Interval{factor, factor});
}

Interval cosine_at_rational(const Rational& argument, int maximum_term) {
  const Rational square = argument * argument;
  Rational term(1);
  Rational sum(1);
  for (int index = 1; index <= maximum_term; ++index) {
    term *= square;
    term /= Rational(static_cast<long long>(2 * index - 1) *
                     static_cast<long long>(2 * index));
    if ((index & 1) == 0) {
      sum += term;
    } else {
      sum -= term;
    }
  }
  Rational next = term * square;
  next /= Rational(static_cast<long long>(2 * maximum_term + 1) *
                   static_cast<long long>(2 * maximum_term + 2));
  if (((maximum_term + 1) & 1) != 0) {
    next = -next;
  }
  return Interval{minimum(sum, sum + next), maximum(sum, sum + next)};
}

std::pair<Rational, Rational> pi_bounds() {
  // Archimedes' rational enclosure, sufficient for the certified small-rank
  // diagnostics below: 333/106 < pi < 355/113.
  return {Rational(333, 106), Rational(355, 113)};
}

Interval cosine(int index, int modulus, int maximum_term) {
  if (index <= 0 || 2 * index >= modulus) {
    throw std::runtime_error("cosine index is outside the acute range");
  }
  const auto [pi_lower, pi_upper] = pi_bounds();
  const Rational lower_argument =
      2 * pi_lower * Rational(index, modulus);
  const Rational upper_argument =
      2 * pi_upper * Rational(index, modulus);
  const Interval lower_value = cosine_at_rational(
      upper_argument, maximum_term
  );
  const Interval upper_value = cosine_at_rational(
      lower_argument, maximum_term
  );
  return Interval{lower_value.lower, upper_value.upper};
}

PairData make_pair(const Interval& first, const Interval& second) {
  const Interval one = exact(1);
  const Interval difference = subtract(first, second);
  if (difference.lower <= 0) {
    throw std::runtime_error("nodes are not strictly ordered");
  }
  const Interval first_square = square_nonnegative(first.lower >= 0
      ? first : negate(first));
  const Interval second_square = square_nonnegative(second.lower >= 0
      ? second : negate(second));
  const Interval product = multiply(first, second);
  const Interval h = add(
      add(exact(6), scale(product, Rational(4))),
      add(scale(add(first_square, second_square), Rational(-8)),
          scale(multiply(first_square, second_square), Rational(16))));
  const Interval common = multiply(
      square_nonnegative(difference), subtract(one, product));
  const Interval endpoint_minus = multiply(subtract(one, first),
                                            subtract(one, second));
  const Interval endpoint_plus = multiply(add(one, first), add(one, second));
  const Interval base_weight = multiply(square_nonnegative(endpoint_minus),
                                        square_nonnegative(endpoint_plus));
  return PairData{endpoint_minus, endpoint_plus, base_weight,
                  multiply(common, h)};
}

ApproximatePairData make_approximate_pair(
    int first, int second, int modulus
) {
  const long double first_angle =
      2.0L * std::numbers::pi_v<long double> *
      static_cast<long double>(first) / static_cast<long double>(modulus);
  const long double second_angle =
      2.0L * std::numbers::pi_v<long double> *
      static_cast<long double>(second) / static_cast<long double>(modulus);
  const long double first_node = std::cos(first_angle);
  const long double second_node = std::cos(second_angle);
  const long double endpoint_minus =
      (1.0L - first_node) * (1.0L - second_node);
  const long double endpoint_plus =
      (1.0L + first_node) * (1.0L + second_node);
  const long double kernel = 6.0L + 4.0L * first_node * second_node -
      8.0L * (first_node * first_node + second_node * second_node) +
      16.0L * first_node * first_node * second_node * second_node;
  const long double amplitude =
      2.0L * (first_node - second_node) * (first_node - second_node) *
      (1.0L - first_node * second_node) * kernel *
      endpoint_minus * endpoint_minus * endpoint_plus * endpoint_plus;
  return ApproximatePairData{
      endpoint_minus, endpoint_plus, std::fabs(amplitude)
  };
}

bool at_least(const Interval& credit, const Interval& load) {
  return credit.lower >= load.upper;
}

bool certified_edge(const PairData& load, const PairData& credit) {
  if (!at_least(credit.endpoint_minus, load.endpoint_minus) ||
      !at_least(credit.endpoint_plus, load.endpoint_plus)) {
    return false;
  }
  const Interval credit_value = multiply(
      credit.signed_base_amplitude, credit.base_weight);
  const Interval load_value = multiply(
      negate(load.signed_base_amplitude), load.base_weight);
  return at_least(credit_value, load_value);
}

bool certified_two_credit_edge(
    const PairData& load,
    const PairData& first_credit,
    const PairData& second_credit,
    const Rational& first_share,
    const Rational& second_share
) {
  if (!at_least(
          multiply(first_credit.endpoint_minus, second_credit.endpoint_minus),
          square_nonnegative(load.endpoint_minus)
      ) ||
      !at_least(
          multiply(first_credit.endpoint_plus, second_credit.endpoint_plus),
          square_nonnegative(load.endpoint_plus)
      )) {
    return false;
  }
  const Interval first_value = multiply(
      scale(first_credit.signed_base_amplitude, first_share),
      first_credit.base_weight
  );
  const Interval second_value = multiply(
      scale(second_credit.signed_base_amplitude, second_share),
      second_credit.base_weight
  );
  const Interval load_value = multiply(
      negate(load.signed_base_amplitude), load.base_weight
  );
  return at_least(
      scale(multiply(first_value, second_value), Rational(4)),
      square_nonnegative(load_value)
  );
}

bool possibly_two_credit_edge(
    const PairData& load,
    const PairData& first_credit,
    const PairData& second_credit,
    const Rational& first_share,
    const Rational& second_share
) {
  // This is a deliberately coarse interval rejection filter.  Returning
  // false proves that the exact edge inequalities fail; returning true is
  // followed by certified_two_credit_edge on the sharper intervals.
  const Interval endpoint_minus_product = multiply(
      first_credit.endpoint_minus, second_credit.endpoint_minus
  );
  const Interval endpoint_minus_square = multiply(
      load.endpoint_minus, load.endpoint_minus
  );
  if (endpoint_minus_product.upper < endpoint_minus_square.lower) {
    return false;
  }
  const Interval endpoint_plus_product = multiply(
      first_credit.endpoint_plus, second_credit.endpoint_plus
  );
  const Interval endpoint_plus_square = multiply(
      load.endpoint_plus, load.endpoint_plus
  );
  if (endpoint_plus_product.upper < endpoint_plus_square.lower) {
    return false;
  }
  const Interval first_value = multiply(
      scale(first_credit.signed_base_amplitude, first_share),
      first_credit.base_weight
  );
  const Interval second_value = multiply(
      scale(second_credit.signed_base_amplitude, second_share),
      second_credit.base_weight
  );
  const Interval load_value = multiply(
      negate(load.signed_base_amplitude), load.base_weight
  );
  const Interval credit_product = scale(
      multiply(first_value, second_value), Rational(4)
  );
  const Interval load_square = multiply(load_value, load_value);
  return credit_product.upper >= load_square.lower;
}

bool approximate_two_credit_edge(
    const ApproximatePairData& load,
    const ApproximatePairData& first_credit,
    const ApproximatePairData& second_credit
) {
  constexpr long double tolerance = 1.0e-12L;
  const auto at_least_approximately = [tolerance](long double left,
                                                   long double right) {
    return left + tolerance *
        (1.0L + std::fabs(left) + std::fabs(right)) >= right;
  };
  return at_least_approximately(
             first_credit.endpoint_minus * second_credit.endpoint_minus,
             load.endpoint_minus * load.endpoint_minus
         ) &&
      at_least_approximately(
          first_credit.endpoint_plus * second_credit.endpoint_plus,
          load.endpoint_plus * load.endpoint_plus
      ) &&
      at_least_approximately(
          first_credit.base_amplitude * second_credit.base_amplitude,
          load.base_amplitude * load.base_amplitude
      );
}

bool augment(int load, const std::vector<std::vector<int>>& adjacency,
             std::vector<int>& matched_credit, std::vector<bool>& seen) {
  for (const int credit : adjacency[static_cast<std::size_t>(load)]) {
    const std::size_t slot = static_cast<std::size_t>(credit);
    if (seen[slot]) {
      continue;
    }
    seen[slot] = true;
    if (matched_credit[slot] < 0 ||
        augment(matched_credit[slot], adjacency, matched_credit, seen)) {
      matched_credit[slot] = load;
      return true;
    }
  }
  return false;
}

int matching_size(const std::vector<std::vector<int>>& adjacency,
                  int credit_count) {
  std::vector<int> matched_credit(static_cast<std::size_t>(credit_count), -1);
  int matched = 0;
  for (int load = 0; load < static_cast<int>(adjacency.size()); ++load) {
    std::vector<bool> seen(static_cast<std::size_t>(credit_count), false);
    if (augment(load, adjacency, matched_credit, seen)) {
      ++matched;
    }
  }
  return matched;
}

bool two_credit_greedy_cover(
    const std::vector<std::vector<std::pair<int, int>>>& adjacency,
    int credit_count,
    unsigned char maximum_uses,
    std::vector<std::pair<int, int>>& cover
) {
  std::vector<unsigned char> used(static_cast<std::size_t>(credit_count), 0U);
  std::vector<unsigned char> assigned(adjacency.size(), 0U);
  cover.assign(adjacency.size(), {-1, -1});
  for (std::size_t depth = 0U; depth < adjacency.size(); ++depth) {
    int chosen_load = -1;
    std::size_t chosen_count = 0U;
    for (std::size_t load = 0U; load < adjacency.size(); ++load) {
      if (assigned[load] != 0U) {
        continue;
      }
      std::size_t feasible = 0U;
      for (const auto& [first, second] : adjacency[load]) {
        if (used[static_cast<std::size_t>(first)] < maximum_uses &&
            used[static_cast<std::size_t>(second)] < maximum_uses) {
          ++feasible;
        }
      }
      if (feasible == 0U) {
        return false;
      }
      if (chosen_load < 0 || feasible < chosen_count) {
        chosen_load = static_cast<int>(load);
        chosen_count = feasible;
      }
    }
    int chosen_first = -1;
    int chosen_second = -1;
    std::size_t chosen_future = 0U;
    for (const auto& [first, second]
         : adjacency[static_cast<std::size_t>(chosen_load)]) {
      const std::size_t first_slot = static_cast<std::size_t>(first);
      const std::size_t second_slot = static_cast<std::size_t>(second);
      if (used[first_slot] >= maximum_uses ||
          used[second_slot] >= maximum_uses) {
        continue;
      }
      ++used[first_slot];
      ++used[second_slot];
      std::size_t future = 0U;
      bool feasible_future = true;
      for (std::size_t load = 0U; load < adjacency.size(); ++load) {
        if (assigned[load] != 0U ||
            static_cast<int>(load) == chosen_load) {
          continue;
        }
        std::size_t options = 0U;
        for (const auto& [next_first, next_second] : adjacency[load]) {
          if (used[static_cast<std::size_t>(next_first)] < maximum_uses &&
              used[static_cast<std::size_t>(next_second)] < maximum_uses) {
            ++options;
          }
        }
        if (options == 0U) {
          feasible_future = false;
          break;
        }
        future += options;
      }
      --used[first_slot];
      --used[second_slot];
      if (feasible_future && (chosen_first < 0 || future > chosen_future)) {
        chosen_first = first;
        chosen_second = second;
        chosen_future = future;
      }
    }
    if (chosen_first < 0) {
      return false;
    }
    ++used[static_cast<std::size_t>(chosen_first)];
    ++used[static_cast<std::size_t>(chosen_second)];
    assigned[static_cast<std::size_t>(chosen_load)] = 1U;
    cover[static_cast<std::size_t>(chosen_load)] = {
        chosen_first, chosen_second
    };
  }
  return true;
}

bool two_credit_random_greedy_cover(
    const std::vector<std::vector<std::pair<int, int>>>& adjacency,
    int credit_count,
    unsigned char maximum_uses,
    std::size_t candidate_width,
    std::mt19937_64& generator,
    std::vector<std::pair<int, int>>& cover
) {
  std::vector<unsigned char> used(static_cast<std::size_t>(credit_count), 0U);
  std::vector<unsigned char> assigned(adjacency.size(), 0U);
  cover.assign(adjacency.size(), {-1, -1});
  for (std::size_t depth = 0U; depth < adjacency.size(); ++depth) {
    std::size_t minimum_feasible = std::numeric_limits<std::size_t>::max();
    std::vector<int> load_choices;
    for (std::size_t load = 0U; load < adjacency.size(); ++load) {
      if (assigned[load] != 0U) {
        continue;
      }
      std::size_t feasible = 0U;
      for (const auto& [first, second] : adjacency[load]) {
        if (used[static_cast<std::size_t>(first)] < maximum_uses &&
            used[static_cast<std::size_t>(second)] < maximum_uses) {
          ++feasible;
        }
      }
      if (feasible == 0U) {
        return false;
      }
      if (feasible < minimum_feasible) {
        minimum_feasible = feasible;
        load_choices.clear();
      }
      if (feasible == minimum_feasible) {
        load_choices.push_back(static_cast<int>(load));
      }
    }
    std::uniform_int_distribution<std::size_t> load_pick(
        0U, load_choices.size() - 1U
    );
    const int chosen_load = load_choices[load_pick(generator)];
    std::vector<std::pair<std::size_t, std::pair<int, int>>> edge_choices;
    for (const auto& [first, second]
         : adjacency[static_cast<std::size_t>(chosen_load)]) {
      const std::size_t first_slot = static_cast<std::size_t>(first);
      const std::size_t second_slot = static_cast<std::size_t>(second);
      if (used[first_slot] >= maximum_uses ||
          used[second_slot] >= maximum_uses) {
        continue;
      }
      ++used[first_slot];
      ++used[second_slot];
      std::size_t future = 0U;
      bool feasible_future = true;
      for (std::size_t load = 0U; load < adjacency.size(); ++load) {
        if (assigned[load] != 0U ||
            static_cast<int>(load) == chosen_load) {
          continue;
        }
        std::size_t options = 0U;
        for (const auto& [next_first, next_second] : adjacency[load]) {
          if (used[static_cast<std::size_t>(next_first)] < maximum_uses &&
              used[static_cast<std::size_t>(next_second)] < maximum_uses) {
            ++options;
          }
        }
        if (options == 0U) {
          feasible_future = false;
          break;
        }
        future += options;
      }
      --used[first_slot];
      --used[second_slot];
      if (!feasible_future) {
        continue;
      }
      edge_choices.emplace_back(future, std::make_pair(first, second));
    }
    if (edge_choices.empty()) {
      return false;
    }
    std::sort(
        edge_choices.begin(), edge_choices.end(),
        [](const auto& left, const auto& right) {
          return left.first > right.first;
        }
    );
    const std::size_t choice_count = std::min<std::size_t>(
        candidate_width, edge_choices.size()
    );
    std::uniform_int_distribution<std::size_t> edge_pick(0U, choice_count - 1U);
    const auto [first, second] = edge_choices[edge_pick(generator)].second;
    ++used[static_cast<std::size_t>(first)];
    ++used[static_cast<std::size_t>(second)];
    assigned[static_cast<std::size_t>(chosen_load)] = 1U;
    cover[static_cast<std::size_t>(chosen_load)] = {first, second};
  }
  return true;
}

bool two_credit_cover_search(
    const std::vector<std::vector<std::pair<int, int>>>& adjacency,
    const std::vector<int>& order,
    std::vector<unsigned char>& used,
    unsigned char maximum_uses,
    std::unordered_set<CoverState, CoverStateHash>& failed_states,
    std::vector<std::pair<int, int>>& cover,
    std::size_t depth
) {
  if (depth == order.size()) {
    return true;
  }
  std::size_t remaining_capacity = 0U;
  CoverState state;
  for (std::size_t index = 0U; index < used.size(); ++index) {
    remaining_capacity += static_cast<std::size_t>(maximum_uses - used[index]);
    const std::uint64_t value = static_cast<std::uint64_t>(used[index]);
    if (index < 32U) {
      state.low |= value << (2U * index);
    } else {
      state.high |= value << (2U * (index - 32U));
    }
  }
  const std::size_t remaining_loads = order.size() - depth;
  if (remaining_capacity < 2U * remaining_loads ||
      failed_states.contains(state)) {
    return false;
  }
  const int load = order[depth];
  for (const auto& [first, second]
       : adjacency[static_cast<std::size_t>(load)]) {
    const std::size_t first_slot = static_cast<std::size_t>(first);
    const std::size_t second_slot = static_cast<std::size_t>(second);
    if (used[first_slot] >= maximum_uses ||
        used[second_slot] >= maximum_uses) {
      continue;
    }
    ++used[first_slot];
    ++used[second_slot];
    cover[static_cast<std::size_t>(load)] = {first, second};
    if (two_credit_cover_search(
            adjacency, order, used, maximum_uses, failed_states, cover,
            depth + 1U
        )) {
      return true;
    }
    --used[first_slot];
    --used[second_slot];
    cover[static_cast<std::size_t>(load)] = {-1, -1};
  }
  failed_states.insert(state);
  return false;
}

bool two_credit_cover(
    const std::vector<std::vector<std::pair<int, int>>>& adjacency,
    int credit_count,
    unsigned char maximum_uses,
    std::vector<std::pair<int, int>>& cover
) {
  std::vector<int> order(adjacency.size());
  for (std::size_t index = 0U; index < order.size(); ++index) {
    order[index] = static_cast<int>(index);
  }
  std::sort(
      order.begin(), order.end(),
      [&adjacency](int left, int right) {
        return adjacency[static_cast<std::size_t>(left)].size() <
               adjacency[static_cast<std::size_t>(right)].size();
      }
  );
  std::vector<unsigned char> used(static_cast<std::size_t>(credit_count), 0U);
  if (used.size() > 64U) {
    throw std::runtime_error("two-credit state exceeds 128-bit memoization");
  }
  if (two_credit_greedy_cover(
          adjacency, credit_count, maximum_uses, cover
      )) {
    return true;
  }
  cover.assign(adjacency.size(), {-1, -1});
  std::unordered_set<CoverState, CoverStateHash> failed_states;
  return two_credit_cover_search(
      adjacency, order, used, maximum_uses, failed_states, cover, 0U
  );
}

bool verifies_two_credit_cover(
    const std::vector<std::vector<std::pair<int, int>>>& adjacency,
    const std::vector<std::pair<int, int>>& cover,
    int credit_count,
    unsigned char maximum_uses
) {
  if (cover.size() != adjacency.size()) {
    return false;
  }
  std::vector<unsigned char> used(static_cast<std::size_t>(credit_count), 0U);
  for (std::size_t load = 0U; load < adjacency.size(); ++load) {
    const auto selected = cover[load];
    if (std::find(
            adjacency[load].begin(), adjacency[load].end(), selected
        ) == adjacency[load].end()) {
      return false;
    }
    const std::size_t first = static_cast<std::size_t>(selected.first);
    const std::size_t second = static_cast<std::size_t>(selected.second);
    ++used[first];
    ++used[second];
    if (used[first] > maximum_uses || used[second] > maximum_uses) {
      return false;
    }
  }
  return true;
}

PartitionMatroidInfo partition_matroid_info(
    const std::vector<std::pair<int, int>>& edges,
    int credit_count
) {
  std::vector<unsigned char> neighbour(
      static_cast<std::size_t>(credit_count), 0U
  );
  for (const auto& [first, second] : edges) {
    neighbour[static_cast<std::size_t>(first)] = 1U;
    neighbour[static_cast<std::size_t>(second)] = 1U;
  }
  std::vector<int> vertices;
  for (int credit = 0; credit < credit_count; ++credit) {
    if (neighbour[static_cast<std::size_t>(credit)] != 0U) {
      vertices.push_back(credit);
    }
  }
  std::vector<std::vector<unsigned char>> eligible(
      vertices.size(), std::vector<unsigned char>(vertices.size(), 0U)
  );
  std::vector<int> local_index(static_cast<std::size_t>(credit_count), -1);
  for (std::size_t index = 0U; index < vertices.size(); ++index) {
    local_index[static_cast<std::size_t>(vertices[index])] =
        static_cast<int>(index);
  }
  for (const auto& [first, second] : edges) {
    const int first_local = local_index[static_cast<std::size_t>(first)];
    const int second_local = local_index[static_cast<std::size_t>(second)];
    eligible[static_cast<std::size_t>(first_local)]
            [static_cast<std::size_t>(second_local)] = 1U;
    eligible[static_cast<std::size_t>(second_local)]
            [static_cast<std::size_t>(first_local)] = 1U;
  }
  std::vector<int> class_index(vertices.size(), -1);
  std::size_t class_count = 0U;
  for (std::size_t seed = 0U; seed < vertices.size(); ++seed) {
    if (class_index[seed] >= 0) {
      continue;
    }
    std::vector<std::size_t> pending{seed};
    class_index[seed] = static_cast<int>(class_count);
    for (std::size_t cursor = 0U; cursor < pending.size(); ++cursor) {
      const std::size_t current = pending[cursor];
      for (std::size_t next = 0U; next < vertices.size(); ++next) {
        if (current == next || eligible[current][next] != 0U ||
            class_index[next] >= 0) {
          continue;
        }
        class_index[next] = static_cast<int>(class_count);
        pending.push_back(next);
      }
    }
    ++class_count;
  }
  for (std::size_t first = 0U; first < vertices.size(); ++first) {
    for (std::size_t second = first + 1U; second < vertices.size(); ++second) {
      const bool same_class = class_index[first] == class_index[second];
      const bool is_eligible = eligible[first][second] != 0U;
      if (same_class == is_eligible) {
        return PartitionMatroidInfo{false, class_count};
      }
    }
  }
  return PartitionMatroidInfo{true, class_count};
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const bool valid_arguments =
        argc == 2 ||
        (argc == 3 &&
         (std::string(argv[2]) == "--two-credit" ||
          std::string(argv[2]) == "--two-credit-half")) ||
        (argc == 4 && std::string(argv[2]) == "--two-credit-half" &&
         (std::string(argv[3]) == "--print-cover" ||
          std::string(argv[3]) == "--greedy-cover" ||
          std::string(argv[3]) == "--coarse-greedy-certificate" ||
          std::string(argv[3]) == "--approx-greedy-certificate" ||
          std::string(argv[3]) == "--approx-greedy-certificate-single" ||
          std::string(argv[3]) == "--approx-random-certificate-single" ||
          std::string(argv[3]) == "--analyze-graph" ||
          std::string(argv[3]) == "--dump-hyperedges" ||
          std::string(argv[3]) == "--uniform-fractional"));
    if (!valid_arguments) {
      throw std::runtime_error(
          "usage: analyze_tp2_mixed_jacobi_matching MAXIMUM_RANK "
          "[--two-credit|--two-credit-half] "
          "[--two-credit-half --print-cover|--greedy-cover|"
          "--coarse-greedy-certificate|--approx-greedy-certificate|"
          "--approx-greedy-certificate-single|"
          "--approx-random-certificate-single|"
          "--analyze-graph|"
          "--dump-hyperedges|"
          "--uniform-fractional]");
    }
    const int maximum_rank = parse_positive(argv[1], "maximum rank");
    const bool two_credit = argc >= 3;
    const bool two_credit_half = argc >= 3 &&
        std::string(argv[2]) == "--two-credit-half";
    const bool print_cover = argc == 4 &&
        std::string(argv[3]) == "--print-cover";
    const bool greedy_cover = argc == 4 &&
        std::string(argv[3]) == "--greedy-cover";
    const bool coarse_greedy_certificate = argc == 4 &&
        std::string(argv[3]) == "--coarse-greedy-certificate";
    const bool approximate_greedy_certificate = argc == 4 &&
        std::string(argv[3]) == "--approx-greedy-certificate";
    const bool approximate_greedy_certificate_single = argc == 4 &&
        std::string(argv[3]) == "--approx-greedy-certificate-single";
    const bool approximate_random_certificate_single = argc == 4 &&
        std::string(argv[3]) == "--approx-random-certificate-single";
    const bool approximate_greedy_mode = approximate_greedy_certificate ||
        approximate_greedy_certificate_single ||
        approximate_random_certificate_single;
    const bool approximate_random_mode = approximate_random_certificate_single;
    const bool analyze_graph = argc == 4 &&
        std::string(argv[3]) == "--analyze-graph";
    const bool dump_hyperedges = argc == 4 &&
        std::string(argv[3]) == "--dump-hyperedges";
    const bool uniform_fractional = argc == 4 &&
        std::string(argv[3]) == "--uniform-fractional";
    std::size_t graphs = 0U;
    std::size_t loads = 0U;
    std::size_t edges = 0U;
    std::size_t two_credit_edges = 0U;
    int first_failure_rank = -1;
    int maximum_loads = 0;
    const int minimum_rank = approximate_greedy_certificate_single
        || approximate_random_certificate_single
        ? maximum_rank : 3;
    for (int rank = minimum_rank; rank <= maximum_rank; ++rank) {
      const int modulus = 2 * rank + 1;
      std::vector<Interval> nodes;
      std::vector<Interval> coarse_nodes;
      for (int index = 1; index <= rank; ++index) {
        nodes.push_back(cosine(index, modulus, 16));
        coarse_nodes.push_back(cosine(index, modulus, 8));
      }
      std::vector<IndexedPairData> negative;
      std::vector<IndexedPairData> positive;
      for (int first = 0; first < rank; ++first) {
        for (int second = first + 1; second < rank; ++second) {
          // Every 5-divisible grid contains the n=5 pair; its H-kernel is
          // exactly zero, so it is neither a load nor a credit.
          if (modulus % 5 == 0 && first == modulus / 5 - 1 &&
              second == 2 * modulus / 5 - 1) {
            continue;
          }
          const PairData pair = make_pair(
              nodes[static_cast<std::size_t>(first)],
              nodes[static_cast<std::size_t>(second)]
          );
          const PairData coarse_pair = make_pair(
              coarse_nodes[static_cast<std::size_t>(first)],
              coarse_nodes[static_cast<std::size_t>(second)]
          );
          const ApproximatePairData approximate_pair = make_approximate_pair(
              first + 1, second + 1, modulus
          );
          if (pair.signed_base_amplitude.upper < 0) {
            negative.push_back(IndexedPairData{
                pair, coarse_pair, approximate_pair, first, second
            });
          } else if (pair.signed_base_amplitude.lower > 0) {
            positive.push_back(IndexedPairData{
                pair, coarse_pair, approximate_pair, first, second
            });
          } else {
            std::cout << "TP2_MIXED_JACOBI_MATCHING result=UNRESOLVED_SIGN"
                      << " rank=" << rank
                      << " first=" << first
                      << " second=" << second
                      << '\n';
            return EXIT_FAILURE;
          }
        }
      }
      std::vector<std::vector<int>> adjacency(negative.size());
      for (std::size_t load = 0U; load < negative.size(); ++load) {
        for (std::size_t credit = 0U; credit < positive.size(); ++credit) {
          if (certified_edge(negative[load].data, positive[credit].data)) {
            adjacency[load].push_back(static_cast<int>(credit));
            ++edges;
          }
        }
      }
      if (!two_credit) {
        const int matched = matching_size(
            adjacency, static_cast<int>(positive.size())
        );
        if (matched != static_cast<int>(negative.size()) &&
            first_failure_rank < 0) {
          first_failure_rank = rank;
          std::cout << "TP2_MIXED_JACOBI_MATCHING first_failure_rank="
                    << rank << " loads=" << negative.size()
                    << " credits=" << positive.size()
                    << " matched=" << matched << '\n';
        }
      }
      if (two_credit) {
        std::vector<std::vector<std::pair<int, int>>> two_credit_adjacency(
            negative.size()
        );
        for (std::size_t load = 0U; load < negative.size(); ++load) {
          for (std::size_t first = 0U; first < positive.size(); ++first) {
            for (std::size_t second = first + 1U;
                 second < positive.size(); ++second) {
              const Rational share = two_credit_half
                  ? Rational(1, 2) : Rational(1);
              const bool candidate_edge = approximate_greedy_mode
                  // Approximation is used only to choose a candidate cover;
                  // every selected pair is interval-certified below.
                  ? approximate_two_credit_edge(
                      negative[load].approximate_data,
                      positive[first].approximate_data,
                      positive[second].approximate_data
                  )
                  : possibly_two_credit_edge(
                      negative[load].coarse_data,
                      positive[first].coarse_data,
                      positive[second].coarse_data,
                      share, share
                  );
              if (candidate_edge &&
                  (coarse_greedy_certificate ||
                   approximate_greedy_mode ||
                   certified_two_credit_edge(
                      negative[load].data, positive[first].data,
                      positive[second].data,
                      share, share
                   ))) {
                two_credit_adjacency[load].emplace_back(
                    static_cast<int>(first), static_cast<int>(second)
                );
                ++two_credit_edges;
              }
            }
          }
        }
        std::vector<int> credit_degree(positive.size(), 0);
        for (const auto& load_edges : two_credit_adjacency) {
          for (const auto& [first, second] : load_edges) {
            ++credit_degree[static_cast<std::size_t>(first)];
            ++credit_degree[static_cast<std::size_t>(second)];
          }
        }
        if (!approximate_greedy_mode) {
          for (auto& load_edges : two_credit_adjacency) {
            std::sort(
                load_edges.begin(), load_edges.end(),
                [&credit_degree](
                    const std::pair<int, int>& left,
                    const std::pair<int, int>& right
                ) {
                  const int left_degree =
                      credit_degree[static_cast<std::size_t>(left.first)] +
                      credit_degree[static_cast<std::size_t>(left.second)];
                  const int right_degree =
                      credit_degree[static_cast<std::size_t>(right.first)] +
                      credit_degree[static_cast<std::size_t>(right.second)];
                  return left_degree < right_degree;
                }
            );
          }
        }
        if (analyze_graph && rank == maximum_rank) {
          for (std::size_t load = 0U; load < negative.size(); ++load) {
            std::vector<unsigned char> neighbour(positive.size(), 0U);
            for (const auto& [first, second] : two_credit_adjacency[load]) {
              neighbour[static_cast<std::size_t>(first)] = 1U;
              neighbour[static_cast<std::size_t>(second)] = 1U;
            }
            const std::size_t neighbour_count = static_cast<std::size_t>(
                std::count(neighbour.begin(), neighbour.end(), 1U)
            );
            const std::size_t complete_edges =
                neighbour_count * (neighbour_count - 1U) / 2U;
            const PartitionMatroidInfo partition = partition_matroid_info(
                two_credit_adjacency[load], static_cast<int>(positive.size())
            );
            std::cout << "TP2_MIXED_JACOBI_HYPERGRAPH rank=" << rank
                      << " load=" << negative[load].first + 1
                      << ',' << negative[load].second + 1
                      << " pair_edges=" << two_credit_adjacency[load].size()
                      << " neighbours=" << neighbour_count
                      << " clique="
                      << (two_credit_adjacency[load].size() == complete_edges
                              ? "yes" : "no")
                      << " partition_matroid="
                      << (partition.valid ? "yes" : "no")
                      << " classes=" << partition.class_count
                      << '\n';
          }
        }
        if (dump_hyperedges && rank == maximum_rank) {
          for (std::size_t load = 0U; load < negative.size(); ++load) {
            std::cout << "TP2_MIXED_JACOBI_HYPEREDGES rank=" << rank
                      << " load=" << negative[load].first + 1
                      << ',' << negative[load].second + 1
                      << " eligible=";
            bool first_edge = true;
            for (const auto& [first, second] : two_credit_adjacency[load]) {
              if (!first_edge) {
                std::cout << ',';
              }
              first_edge = false;
              const IndexedPairData& first_credit =
                  positive[static_cast<std::size_t>(first)];
              const IndexedPairData& second_credit =
                  positive[static_cast<std::size_t>(second)];
              std::cout << first_credit.first + 1 << ':'
                        << first_credit.second + 1 << ';'
                        << second_credit.first + 1 << ':'
                        << second_credit.second + 1;
            }
            std::cout << '\n';
          }
        }
        if (greedy_cover || coarse_greedy_certificate ||
            approximate_greedy_mode) {
          std::vector<std::pair<int, int>> cover;
          const unsigned char maximum_uses = two_credit_half
              ? static_cast<unsigned char>(2U)
              : static_cast<unsigned char>(1U);
          bool covered = false;
          if (approximate_random_mode) {
            constexpr std::uint64_t attempts = 256U;
            const std::size_t candidate_width = rank <= 22 ? 16U : 32U;
            for (std::uint64_t attempt = 0U; attempt < attempts; ++attempt) {
              std::mt19937_64 generator(
                  UINT64_C(0x9e3779b97f4a7c15) ^
                  static_cast<std::uint64_t>(rank) *
                      UINT64_C(0xbf58476d1ce4e5b9) ^ attempt
              );
              if (two_credit_random_greedy_cover(
                  two_credit_adjacency,
                  static_cast<int>(positive.size()), maximum_uses,
                  candidate_width,
                  generator, cover
              )) {
                covered = true;
                break;
              }
            }
          } else {
            covered = two_credit_greedy_cover(
                two_credit_adjacency, static_cast<int>(positive.size()),
                maximum_uses, cover
            );
          }
          if (covered && !verifies_two_credit_cover(
                  two_credit_adjacency, cover,
                  static_cast<int>(positive.size()), maximum_uses
              )) {
            throw std::runtime_error("greedy cover verification failure");
          }
          if (covered && (coarse_greedy_certificate ||
                          approximate_greedy_mode)) {
            const Rational share(1, 2);
            for (std::size_t load = 0U; load < negative.size(); ++load) {
              const auto [first, second] = cover[load];
              if (!certified_two_credit_edge(
                      negative[load].data,
                      positive[static_cast<std::size_t>(first)].data,
                      positive[static_cast<std::size_t>(second)].data,
                      share, share
                  )) {
                covered = false;
                break;
              }
            }
          }
          if (!covered && first_failure_rank < 0) {
            first_failure_rank = rank;
            std::cout << (coarse_greedy_certificate
                              ? "TP2_MIXED_JACOBI_COARSE_CANDIDATE"
                              : (approximate_random_mode
                                  ? "TP2_MIXED_JACOBI_APPROX_RANDOM_CANDIDATE"
                                  : (approximate_greedy_mode
                                  ? "TP2_MIXED_JACOBI_APPROX_CANDIDATE"
                                  : "TP2_MIXED_JACOBI_GREEDY_COVER")))
                      << " first_failure_rank=" << rank
                      << " loads=" << negative.size()
                      << " credits=" << positive.size()
                      << " two_credit_edges=" << two_credit_edges << '\n';
          }
          if (rank == maximum_rank && covered) {
            for (std::size_t load = 0U; load < negative.size(); ++load) {
              const auto [first, second] = cover[load];
              std::cout << (coarse_greedy_certificate
                                ? "TP2_MIXED_JACOBI_COARSE_CERTIFICATE"
                                : (approximate_random_mode
                                    ? "TP2_MIXED_JACOBI_APPROX_RANDOM_CERTIFICATE"
                                    : (approximate_greedy_mode
                                    ? "TP2_MIXED_JACOBI_APPROX_CERTIFICATE"
                                    : "TP2_MIXED_JACOBI_GREEDY_COVER")))
                        << " rank=" << rank
                        << " load=" << negative[load].first + 1
                        << ',' << negative[load].second + 1
                        << " credits="
                        << positive[static_cast<std::size_t>(first)].first + 1
                        << ','
                        << positive[static_cast<std::size_t>(first)].second + 1
                        << ';'
                        << positive[static_cast<std::size_t>(second)].first + 1
                        << ','
                        << positive[static_cast<std::size_t>(second)].second + 1
                        << '\n';
            }
          }
        } else if (uniform_fractional) {
          std::vector<Rational> use(positive.size(), Rational(0));
          for (const auto& load_edges : two_credit_adjacency) {
            if (load_edges.empty()) {
              throw std::runtime_error("load has no two-credit edge");
            }
            const Rational share(
                1, static_cast<int>(load_edges.size())
            );
            for (const auto& [first, second] : load_edges) {
              use[static_cast<std::size_t>(first)] += share;
              use[static_cast<std::size_t>(second)] += share;
            }
          }
          const auto maximum = std::max_element(use.begin(), use.end());
          if (*maximum > Rational(2) && first_failure_rank < 0) {
            first_failure_rank = rank;
            std::cout << "TP2_MIXED_JACOBI_UNIFORM_FRACTIONAL"
                      << " first_failure_rank=" << rank
                      << " maximum_credit_use=" << *maximum << '\n';
          }
        } else if (!analyze_graph && !dump_hyperedges) {
          std::vector<std::pair<int, int>> cover;
          const bool covered = two_credit_cover(
              two_credit_adjacency, static_cast<int>(positive.size()),
              two_credit_half ? static_cast<unsigned char>(2U)
                              : static_cast<unsigned char>(1U),
              cover
          );
          if (covered && !verifies_two_credit_cover(
                  two_credit_adjacency, cover,
                  static_cast<int>(positive.size()),
                  two_credit_half ? static_cast<unsigned char>(2U)
                                  : static_cast<unsigned char>(1U)
              )) {
            throw std::runtime_error("two-credit cover verification failure");
          }
          if (!covered && first_failure_rank < 0) {
            first_failure_rank = rank;
            std::cout << (two_credit_half
                              ? "TP2_MIXED_JACOBI_TWO_CREDIT_HALF"
                              : "TP2_MIXED_JACOBI_TWO_CREDIT")
                      << " first_failure_rank="
                      << rank << " loads=" << negative.size()
                      << " credits=" << positive.size()
                      << " two_credit_edges=" << two_credit_edges << '\n';
          }
          if (print_cover && rank == maximum_rank && covered) {
            for (std::size_t load = 0U; load < negative.size(); ++load) {
              const auto [first, second] = cover[load];
              std::cout << "TP2_MIXED_JACOBI_COVER rank=" << rank
                        << " load=" << negative[load].first + 1
                        << ',' << negative[load].second + 1
                        << " credits="
                        << positive[static_cast<std::size_t>(first)].first + 1
                        << ','
                        << positive[static_cast<std::size_t>(first)].second + 1
                        << ';'
                        << positive[static_cast<std::size_t>(second)].first + 1
                        << ','
                        << positive[static_cast<std::size_t>(second)].second + 1
                        << '\n';
            }
          }
        }
      }
      ++graphs;
      loads += negative.size();
      maximum_loads = std::max(maximum_loads, static_cast<int>(negative.size()));
    }
    std::cout << "TP2_MIXED_JACOBI_MATCHING maximum_rank=" << maximum_rank
              << " graphs=" << graphs << " loads=" << loads
              << " edges=" << edges << " maximum_loads=" << maximum_loads
              << " two_credit_edges=" << two_credit_edges
              << " first_failure_rank=" << first_failure_rank
              << " result=" << (analyze_graph || dump_hyperedges
                      ? "GRAPH_ONLY"
                      : (greedy_cover || coarse_greedy_certificate ||
                         approximate_greedy_mode
                          ? (first_failure_rank < 0
                              ? (coarse_greedy_certificate
                                  ? "COARSE_CERTIFICATE_PASS"
                                  : (approximate_random_mode
                                      ? "APPROX_RANDOM_CERTIFICATE_PASS"
                                      : (approximate_greedy_mode
                                      ? "APPROX_CERTIFICATE_PASS"
                                      : "GREEDY_PASS")))
                              : (coarse_greedy_certificate
                                  ? "COARSE_CANDIDATE_FAIL"
                                  : (approximate_random_mode
                                      ? "APPROX_RANDOM_CANDIDATE_FAIL"
                                      : (approximate_greedy_mode
                                      ? "APPROX_CANDIDATE_FAIL"
                                      : "GREEDY_FAIL"))))
                          : (first_failure_rank < 0 ? "PASS" : "FAIL")))
              << '\n';
    return analyze_graph || dump_hyperedges || first_failure_rank < 0
        ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "TP2_MIXED_JACOBI_MATCHING error=" << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
