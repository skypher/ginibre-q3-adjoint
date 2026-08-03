#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
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
  const cpp_int denominator("100000000000000000000");
  const cpp_int lower_numerator("314159265358979323846");
  return {Rational(lower_numerator, denominator),
          Rational(lower_numerator + 1, denominator)};
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

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2) {
      throw std::runtime_error(
          "usage: analyze_tp2_mixed_jacobi_matching MAXIMUM_RANK");
    }
    const int maximum_rank = parse_positive(argv[1], "maximum rank");
    std::size_t graphs = 0U;
    std::size_t loads = 0U;
    std::size_t edges = 0U;
    int first_failure_rank = -1;
    int maximum_loads = 0;
    for (int rank = 3; rank <= maximum_rank; ++rank) {
      const int modulus = 2 * rank + 1;
      std::vector<Interval> nodes;
      for (int index = 1; index <= rank; ++index) {
        nodes.push_back(cosine(index, modulus));
      }
      std::vector<PairData> negative;
      std::vector<PairData> positive;
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
            negative.push_back(pair);
          } else if (pair.signed_base_amplitude.lower > 0) {
            positive.push_back(pair);
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
          if (certified_edge(negative[load], positive[credit])) {
            adjacency[load].push_back(static_cast<int>(credit));
            ++edges;
          }
        }
      }
      const int matched = matching_size(adjacency, static_cast<int>(positive.size()));
      if (matched != static_cast<int>(negative.size()) && first_failure_rank < 0) {
        first_failure_rank = rank;
        std::cout << "TP2_MIXED_JACOBI_MATCHING first_failure_rank="
                  << rank << " loads=" << negative.size()
                  << " credits=" << positive.size()
                  << " matched=" << matched << '\n';
      }
      ++graphs;
      loads += negative.size();
      maximum_loads = std::max(maximum_loads, static_cast<int>(negative.size()));
    }
    std::cout << "TP2_MIXED_JACOBI_MATCHING maximum_rank=" << maximum_rank
              << " graphs=" << graphs << " loads=" << loads
              << " edges=" << edges << " maximum_loads=" << maximum_loads
              << " first_failure_rank=" << first_failure_rank
              << " result=" << (first_failure_rank < 0 ? "PASS" : "FAIL")
              << '\n';
    return first_failure_rank < 0 ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "TP2_MIXED_JACOBI_MATCHING error=" << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
