#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
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

struct IndexedPairData {
  PairData data;
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

Interval cosine_at_rational(const Rational& argument) {
  const Rational square = argument * argument;
  Rational term(1);
  Rational sum(1);
  constexpr int maximum_term = 16;
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

Interval cosine(int index, int modulus) {
  if (index <= 0 || 2 * index >= modulus) {
    throw std::runtime_error("cosine index is outside the acute range");
  }
  const auto [pi_lower, pi_upper] = pi_bounds();
  const Rational lower_argument =
      2 * pi_lower * Rational(index, modulus);
  const Rational upper_argument =
      2 * pi_upper * Rational(index, modulus);
  const Interval lower_value = cosine_at_rational(upper_argument);
  const Interval upper_value = cosine_at_rational(lower_argument);
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
          std::string(argv[3]) == "--analyze-graph" ||
          std::string(argv[3]) == "--uniform-fractional"));
    if (!valid_arguments) {
      throw std::runtime_error(
          "usage: analyze_tp2_mixed_jacobi_matching MAXIMUM_RANK "
          "[--two-credit|--two-credit-half] "
          "[--two-credit-half --print-cover|--analyze-graph|"
          "--uniform-fractional]");
    }
    const int maximum_rank = parse_positive(argv[1], "maximum rank");
    const bool two_credit = argc >= 3;
    const bool two_credit_half = argc >= 3 &&
        std::string(argv[2]) == "--two-credit-half";
    const bool print_cover = argc == 4 &&
        std::string(argv[3]) == "--print-cover";
    const bool analyze_graph = argc == 4 &&
        std::string(argv[3]) == "--analyze-graph";
    const bool uniform_fractional = argc == 4 &&
        std::string(argv[3]) == "--uniform-fractional";
    std::size_t graphs = 0U;
    std::size_t loads = 0U;
    std::size_t edges = 0U;
    std::size_t two_credit_edges = 0U;
    int first_failure_rank = -1;
    int maximum_loads = 0;
    for (int rank = 3; rank <= maximum_rank; ++rank) {
      const int modulus = 2 * rank + 1;
      std::vector<Interval> nodes;
      for (int index = 1; index <= rank; ++index) {
        nodes.push_back(cosine(index, modulus));
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
          if (pair.signed_base_amplitude.upper < 0) {
            negative.push_back(IndexedPairData{pair, first, second});
          } else if (pair.signed_base_amplitude.lower > 0) {
            positive.push_back(IndexedPairData{pair, first, second});
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
              if (certified_two_credit_edge(
                      negative[load].data, positive[first].data,
                      positive[second].data,
                      share, share
                  )) {
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
            std::cout << "TP2_MIXED_JACOBI_HYPERGRAPH rank=" << rank
                      << " load=" << negative[load].first + 1
                      << ',' << negative[load].second + 1
                      << " pair_edges=" << two_credit_adjacency[load].size()
                      << " neighbours=" << neighbour_count
                      << " clique="
                      << (two_credit_adjacency[load].size() == complete_edges
                              ? "yes" : "no")
                      << '\n';
          }
        }
        if (uniform_fractional) {
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
        } else if (!analyze_graph) {
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
              << " result=" << (analyze_graph ? "GRAPH_ONLY"
                                                : (first_failure_rank < 0
                                                       ? "PASS" : "FAIL"))
              << '\n';
    return analyze_graph || first_failure_rank < 0 ? EXIT_SUCCESS
                                                    : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "TP2_MIXED_JACOBI_MATCHING error=" << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
