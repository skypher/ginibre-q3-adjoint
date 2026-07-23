#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Vec = std::vector<cpp_int>;

cpp_int binomial(int top, int bottom) {
  if (bottom < 0 || bottom > top) {
    return 0;
  }
  bottom = std::min(bottom, top - bottom);
  cpp_int result = 1;
  for (int factor = 1; factor <= bottom; ++factor) {
    result *= top - bottom + factor;
    result /= factor;
  }
  return result;
}

Vec separated_pole_coefficients(int rank, int minus_pairs,
                                int maximum_index) {
  if (maximum_index < 0) {
    return {};
  }
  Vec numerator(static_cast<std::size_t>(maximum_index + 1));
  const int numerator_power = minus_pairs - 1;
  for (int degree = 0;
       degree <= std::min(numerator_power, maximum_index); ++degree) {
    cpp_int value = binomial(numerator_power, degree);
    value <<= 2 * degree;
    if ((numerator_power - degree) % 2 != 0) {
      value = -value;
    }
    numerator[static_cast<std::size_t>(degree)] = value;
  }

  Vec denominator(static_cast<std::size_t>(rank + 1));
  for (int degree = 0; degree <= rank; ++degree) {
    cpp_int value = binomial(2 * rank - degree, degree);
    if (degree % 2 != 0) {
      value = -value;
    }
    denominator[static_cast<std::size_t>(degree)] = value;
  }

  Vec coefficients(static_cast<std::size_t>(maximum_index + 1));
  for (int degree = 0; degree <= maximum_index; ++degree) {
    cpp_int value = numerator[static_cast<std::size_t>(degree)];
    for (int offset = 1; offset <= std::min(rank, degree); ++offset) {
      value -= denominator[static_cast<std::size_t>(offset)] *
               coefficients[static_cast<std::size_t>(degree - offset)];
    }
    coefficients[static_cast<std::size_t>(degree)] = value;
  }
  return coefficients;
}

cpp_int cyclic_weighted_moment(int modulus, int minus_pairs, int power) {
  const int degree = 2 * minus_pairs + power;
  Vec coefficients(static_cast<std::size_t>(2 * degree + 1));
  coefficients[static_cast<std::size_t>(degree)] = 1;
  int radius = 0;
  for (int factor = 0; factor < minus_pairs; ++factor) {
    Vec next(coefficients.size());
    for (int exponent = -radius; exponent <= radius; ++exponent) {
      const cpp_int& value =
          coefficients[static_cast<std::size_t>(degree + exponent)];
      next[static_cast<std::size_t>(degree + exponent)] += 2 * value;
      next[static_cast<std::size_t>(degree + exponent - 2)] -= value;
      next[static_cast<std::size_t>(degree + exponent + 2)] -= value;
    }
    coefficients = std::move(next);
    radius += 2;
  }
  for (int factor = 0; factor < power; ++factor) {
    Vec next(coefficients.size());
    for (int exponent = -radius; exponent <= radius; ++exponent) {
      const cpp_int& value =
          coefficients[static_cast<std::size_t>(degree + exponent)];
      next[static_cast<std::size_t>(degree + exponent - 1)] += value;
      next[static_cast<std::size_t>(degree + exponent + 1)] += value;
    }
    coefficients = std::move(next);
    ++radius;
  }
  cpp_int result = 0;
  for (int exponent = -degree; exponent <= degree; ++exponent) {
    if (exponent % modulus == 0) {
      result += coefficients[static_cast<std::size_t>(degree + exponent)];
    }
  }
  return result;
}

Vec apply_tadpole(const Vec& input, int rank, int sign_second) {
  Vec output(static_cast<std::size_t>(rank * rank));
  const auto index = [rank](int first, int second) {
    return static_cast<std::size_t>(first * rank + second);
  };
  for (int first = 0; first < rank; ++first) {
    for (int second = 0; second < rank; ++second) {
      const cpp_int& value = input[index(first, second)];
      if (value == 0) {
        continue;
      }
      if (first > 0) {
        output[index(first - 1, second)] += value;
      }
      if (first + 1 < rank) {
        output[index(first + 1, second)] += value;
      } else {
        output[index(first, second)] += value;
      }
      if (second > 0) {
        output[index(first, second - 1)] += sign_second * value;
      }
      if (second + 1 < rank) {
        output[index(first, second + 1)] += sign_second * value;
      } else {
        output[index(first, second)] += sign_second * value;
      }
    }
  }
  return output;
}

Vec apply_path(const Vec& input, int vertices, int sign_second) {
  Vec output(static_cast<std::size_t>(vertices * vertices));
  const auto index = [vertices](int first, int second) {
    return static_cast<std::size_t>(first * vertices + second);
  };
  for (int first = 0; first < vertices; ++first) {
    for (int second = 0; second < vertices; ++second) {
      const cpp_int& value = input[index(first, second)];
      if (value == 0) {
        continue;
      }
      if (first > 0) {
        output[index(first - 1, second)] += value;
      }
      if (first + 1 < vertices) {
        output[index(first + 1, second)] += value;
      }
      if (second > 0) {
        output[index(first, second - 1)] += sign_second * value;
      }
      if (second + 1 < vertices) {
        output[index(first, second + 1)] += sign_second * value;
      }
    }
  }
  return output;
}

Vec apply_tadpole_side(const Vec& input, int rank, bool first_side) {
  Vec output(static_cast<std::size_t>(rank * rank));
  const auto index = [rank](int first, int second) {
    return static_cast<std::size_t>(first * rank + second);
  };
  for (int first = 0; first < rank; ++first) {
    for (int second = 0; second < rank; ++second) {
      const cpp_int& value = input[index(first, second)];
      if (value == 0) {
        continue;
      }
      const int coordinate = first_side ? first : second;
      if (coordinate > 0) {
        const int next_first = first_side ? first - 1 : first;
        const int next_second = first_side ? second : second - 1;
        output[index(next_first, next_second)] += value;
      }
      if (coordinate + 1 < rank) {
        const int next_first = first_side ? first + 1 : first;
        const int next_second = first_side ? second : second + 1;
        output[index(next_first, next_second)] += value;
      } else {
        output[index(first, second)] += value;
      }
    }
  }
  return output;
}

Vec apply_fourth_character_side(const Vec& input, int rank, bool first_side) {
  const Vec square = apply_tadpole_side(
      apply_tadpole_side(input, rank, first_side), rank, first_side);
  const Vec fourth = apply_tadpole_side(
      apply_tadpole_side(square, rank, first_side), rank, first_side);
  Vec output(input.size());
  for (std::size_t index = 0; index < input.size(); ++index) {
    output[index] = fourth[index] - 3 * square[index] + input[index];
  }
  return output;
}

cpp_int tadpole_value(int rank, int minus_pairs, int plus_power) {
  Vec state(static_cast<std::size_t>(rank * rank));
  state[0] = 1;
  for (int step = 0; step < 2 * minus_pairs; ++step) {
    state = apply_tadpole(state, rank, -1);
  }
  for (int step = 0; step < plus_power; ++step) {
    state = apply_tadpole(state, rank, +1);
  }
  return state[0];
}

cpp_int unfolded_value(int rank, int minus_pairs, int plus_power) {
  const int vertices = 2 * rank;
  Vec state(static_cast<std::size_t>(vertices * vertices));
  const auto index = [vertices](int first, int second) {
    return static_cast<std::size_t>(first * vertices + second);
  };
  state[index(vertices - 1, 0)] = 1;
  for (int step = 0; step < plus_power; ++step) {
    state = apply_path(state, vertices, +1);
  }
  for (int step = 0; step < 2 * minus_pairs; ++step) {
    state = apply_path(state, vertices, -1);
  }
  return 2 * state[0];
}

cpp_int top_partial_value(int rank, int first_excess, int second_excess) {
  const int total_minus_pairs = 4 + first_excess + second_excess;
  const int plus_power = 4 + 2 * first_excess;
  Vec state(static_cast<std::size_t>(rank * rank));
  state[0] = 1;
  for (int step = 0; step < 2 * total_minus_pairs; ++step) {
    state = apply_tadpole(state, rank, -1);
  }
  for (int step = 0; step < plus_power; ++step) {
    state = apply_tadpole(state, rank, +1);
  }
  const Vec first = apply_fourth_character_side(state, rank, true);
  const Vec second = apply_fourth_character_side(state, rank, false);
  return first[0] + second[0];
}

cpp_int top_partial_moment_value(int rank, int first_excess,
                                 int second_excess) {
  const int modulus = 2 * rank + 1;
  const int minus_pairs = 4 + first_excess + second_excess;
  const int initial_power = 4 + 2 * first_excess;
  Vec moments(5U);
  for (int index = 0; index < 5; ++index) {
    moments[static_cast<std::size_t>(index)] = cyclic_weighted_moment(
        modulus, minus_pairs, initial_power + 2 * index);
  }
  const cpp_int base = moments[0] * moments[2] - moments[1] * moments[1];
  Vec second(3U);
  Vec fourth(3U);
  for (int index = 0; index < 3; ++index) {
    second[static_cast<std::size_t>(index)] =
        moments[static_cast<std::size_t>(index + 1)] -
        2 * moments[static_cast<std::size_t>(index)];
    fourth[static_cast<std::size_t>(index)] =
        moments[static_cast<std::size_t>(index + 2)] -
        4 * moments[static_cast<std::size_t>(index + 1)] +
        2 * moments[static_cast<std::size_t>(index)];
  }
  const cpp_int second_determinant =
      second[0] * second[2] - second[1] * second[1];
  const cpp_int fourth_determinant =
      fourth[0] * fourth[2] - fourth[1] * fourth[1];
  const cpp_int twice_result =
      2 * base + second_determinant + fourth_determinant;
  const cpp_int sixth0 = moments[3] - 6 * moments[2] + 9 * moments[1] -
                         2 * moments[0];
  const cpp_int sixth1 = moments[4] - 6 * moments[3] + 9 * moments[2] -
                         2 * moments[1];
  const cpp_int christoffel_darboux =
      sixth1 * fourth[0] - sixth0 * fourth[1];
  if (twice_result != christoffel_darboux) {
    throw std::runtime_error("Christoffel-Darboux reduction mismatch");
  }
  if (twice_result % 2 != 0) {
    throw std::runtime_error("nonintegral top-partial moment formula");
  }
  return twice_result / 2;
}

}  // namespace

int main(int argc, char** argv) {
  const int maximum_rank = argc > 1 ? std::atoi(argv[1]) : 8;
  const int maximum_minus_pairs = argc > 2 ? std::atoi(argv[2]) : 8;
  const int maximum_plus_power = argc > 3 ? std::atoi(argv[3]) : 15;
  const bool verbose = argc > 4 && std::string(argv[4]) == "verbose";
  if (maximum_rank < 1 || maximum_minus_pairs < 0 || maximum_plus_power < 0) {
    std::cerr << "usage: analyze_su2_tadpole_ray [maximum-rank>=1] "
                 "[maximum-minus-pairs>=0] [maximum-plus-power>=0]\n";
    return 2;
  }

  std::size_t checks = 0;
  std::size_t rational_checks = 0;
  std::size_t partial_checks = 0;
  for (int rank = 1; rank <= maximum_rank; ++rank) {
    for (int minus_pairs = 0; minus_pairs <= maximum_minus_pairs;
         ++minus_pairs) {
      for (int plus_power = 0; plus_power <= maximum_plus_power;
           ++plus_power) {
        const cpp_int direct = tadpole_value(rank, minus_pairs, plus_power);
        const cpp_int unfolded = plus_power % 2 == 1
                                     ? unfolded_value(rank, minus_pairs,
                                                      plus_power)
                                     : direct;
        cpp_int moment0 = 0;
        cpp_int moment2 = 0;
        cpp_int moment4 = 0;
        cpp_int turan = 2 * direct;
        if (plus_power % 2 == 1) {
          const int modulus = 2 * rank + 1;
          moment0 = cyclic_weighted_moment(modulus, minus_pairs, plus_power);
          moment2 =
              cyclic_weighted_moment(modulus, minus_pairs, plus_power + 2);
          moment4 =
              cyclic_weighted_moment(modulus, minus_pairs, plus_power + 4);
          turan = moment2 * moment2 - moment0 * moment4;
          if (minus_pairs > 0) {
            const int half_plus = (plus_power - 1) / 2;
            const int index = minus_pairs + half_plus - rank;
            const Vec coefficients = separated_pole_coefficients(
                rank, minus_pairs, index + 2);
            const auto coefficient = [&coefficients](int wanted) {
              return wanted < 0 ? cpp_int(0)
                                : coefficients[static_cast<std::size_t>(wanted)];
            };
            if (moment0 != -2 * coefficient(index) ||
                moment2 != -2 * coefficient(index + 1) ||
                moment4 != -2 * coefficient(index + 2)) {
              std::cout << "rank=" << rank
                        << " minus_pairs=" << minus_pairs
                        << " plus_power=" << plus_power
                        << " result=RATIONAL_MISMATCH\n";
              return 1;
            }
            ++rational_checks;
          }
        }
        ++checks;
        if (direct != unfolded || 2 * direct != turan || direct < 0) {
          std::cout << "rank=" << rank << " minus_pairs=" << minus_pairs
                    << " plus_power=" << plus_power << " direct=" << direct
                    << " unfolded=" << unfolded << " turan=" << turan
                    << " result=FAIL\n";
          return 1;
        }
        if (verbose && rank <= 4 && plus_power % 2 == 1 && direct != 0) {
          std::cout << "value rank=" << rank << " minus_pairs=" << minus_pairs
                    << " plus_power=" << plus_power << " moments=" << moment0
                    << ',' << moment2 << ',' << moment4 << " value=" << direct
                    << '\n';
        }
      }
    }
  }
  for (int rank = 5; rank <= maximum_rank; ++rank) {
    for (int first_excess = 0; first_excess <= maximum_minus_pairs;
         ++first_excess) {
      for (int second_excess = 0; second_excess <= maximum_minus_pairs;
           ++second_excess) {
        const cpp_int direct =
            top_partial_value(rank, first_excess, second_excess);
        const cpp_int moment =
            top_partial_moment_value(rank, first_excess, second_excess);
        ++partial_checks;
        if (direct != moment || direct < 0) {
          std::cout << "rank=" << rank << " first_excess=" << first_excess
                    << " second_excess=" << second_excess
                    << " direct=" << direct << " moment=" << moment
                    << " result=PARTIAL_FAIL\n";
          return 1;
        }
      }
    }
  }
  std::cout << "SU2_TADPOLE_RAY maximum_rank=" << maximum_rank
            << " maximum_minus_pairs=" << maximum_minus_pairs
            << " maximum_plus_power=" << maximum_plus_power
            << " checks=" << checks
            << " rational_checks=" << rational_checks
            << " partial_checks=" << partial_checks << " result=PASS\n";
  return 0;
}
