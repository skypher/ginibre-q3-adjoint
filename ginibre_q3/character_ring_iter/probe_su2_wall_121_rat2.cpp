#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;

struct Counters {
  std::size_t tails = 0U;
  std::size_t critical = 0U;
  std::size_t smallest_support = 0U;
};

Rational power(Rational base, int exponent) {
  Rational result(1);
  while (exponent > 0) {
    if ((exponent & 1) != 0) {
      result *= base;
    }
    base *= base;
    exponent /= 2;
  }
  return result;
}

Rational at(const std::vector<Rational>& values, std::size_t index) {
  return index < values.size() ? values[index] : Rational(0);
}

Rational terminal_current(const std::vector<Rational>& p) {
  Rational paired(0);
  for (std::size_t index = 1U; index < p.size(); ++index) {
    const Rational g = power(at(p, index), 2)
        - at(p, index - 1U) * at(p, index + 1U);
    const Rational next_g = power(at(p, index + 1U), 2)
        - at(p, index) * at(p, index + 2U);
    const Rational h = at(p, index) * at(p, index - 1U)
        - (index >= 2U ? at(p, index - 2U) * at(p, index + 1U)
                       : Rational(0));
    const Rational next_h = at(p, index + 1U) * at(p, index)
        - at(p, index - 1U) * at(p, index + 2U);
    const Rational next_next_h = at(p, index + 2U) * at(p, index + 1U)
        - at(p, index) * at(p, index + 3U);
    const Rational k = g + h + next_h;
    const Rational next_k = next_g + next_h + next_next_h;
    paired += (g - next_g) * (k - next_k);
  }
  return power(p[1], 3) * p[2] + paired;
}

void check_tail(
    const std::vector<Rational>& ratios,
    Counters& counters,
    const bool expect_unit_floor_payment_obstruction = false
) {
  // The wall tail is p_1=1, p_{i+1}=r_i p_i; ratios are nonincreasing.
  std::vector<Rational> p(ratios.size() + 2U, Rational(0));
  p[1] = Rational(1);
  for (std::size_t index = 0U; index < ratios.size(); ++index) {
    p[index + 2U] = p[index + 1U] * ratios[index];
  }
  if (p.size() < 5U) {
    return;
  }

  const Rational C_zero = terminal_current(p);
  if (C_zero < 1) {
    std::cout << "SU2_WALL_121_UNIT_FLOOR_COUNTEREXAMPLE ratios=[";
    for (std::size_t index = 0U; index < ratios.size(); ++index) {
      if (index != 0U) {
        std::cout << ',';
      }
      std::cout << ratios[index].numerator() << '/'
                << ratios[index].denominator();
    }
    std::cout << "] C0=" << C_zero.numerator() << '/'
              << C_zero.denominator() << '\n';
    std::exit(EXIT_FAILURE);
  }

  const Rational a = p[1];
  const Rational b = p[2];
  const Rational c = p[3];
  const Rational d = p[4];
  const Rational c_one =
      power(b, 3) + a * power(c, 2) + power(c, 3)
      - 3 * power(a, 2) * b - 2 * a * power(b, 2)
      - a * b * c - a * b * d - b * c * d;
  const Rational A = power(a, 2) + power(b, 2) + (a + b) * c;
  const Rational B = -c_one;
  const Rational C = a + b + c;
  if (B <= 0 || A <= 0 || C <= 0 || b <= 0) {
    return;
  }
  const Rational Y = power(a, 2) / b;
  const Rational derivative = -B + 2 * A * Y + 3 * C * power(Y, 2);
  if (derivative <= 0) {
    return;
  }
  ++counters.critical;

  const Rational H = 4 * power(A, 2) + 3 * B * C;
  const Rational J = 2 * power(A, 2) + 3 * B * C;
  const Rational K = 8 * power(A, 2) * J + 3 * B * C * H;
  const Rational N = 4 * A * B * J;
  const Rational unit_floor_margin =
      (4 * A - power(B, 2)) * power(K, 4)
      + 4 * A * C * power(N, 3) * K
      + 9 * power(C, 2) * power(N, 4);
  if (expect_unit_floor_payment_obstruction) {
    if (unit_floor_margin >= 0) {
      throw std::runtime_error("expected unit-floor payment obstruction");
    }
    std::cout << "SU2_WALL_121_RAT2_UNIT_FLOOR_PAYMENT_OBSTRUCTION"
              << " ratios=[";
    for (std::size_t index = 0U; index < ratios.size(); ++index) {
      if (index != 0U) {
        std::cout << ',';
      }
      std::cout << ratios[index].numerator() << '/'
                << ratios[index].denominator();
    }
    std::cout << "] C0=" << C_zero.numerator() << '/'
              << C_zero.denominator()
              << " A=" << A.numerator() << '/' << A.denominator()
              << " B=" << B.numerator() << '/' << B.denominator()
              << " C=" << C.numerator() << '/' << C.denominator()
              << " margin=" << unit_floor_margin.numerator() << '/'
              << unit_floor_margin.denominator()
              << " result=PASS_EXACT_OBSTRUCTION\n";
    std::exit(EXIT_SUCCESS);
  }
  const Rational margin =
      (4 * A * C_zero - power(B, 2)) * power(K, 4)
      + 4 * A * C * power(N, 3) * K
      + 9 * power(C, 2) * power(N, 4);
  if (margin < 0) {
    std::cout << "SU2_WALL_121_RAT2_COUNTEREXAMPLE ratios=[";
    for (std::size_t index = 0U; index < ratios.size(); ++index) {
      if (index != 0U) {
        std::cout << ',';
      }
      std::cout << ratios[index].numerator() << '/' << ratios[index].denominator();
    }
    std::cout << "] C0=" << C_zero.numerator() << '/' << C_zero.denominator()
              << " A=" << A.numerator() << '/' << A.denominator()
              << " B=" << B.numerator() << '/' << B.denominator()
              << " C=" << C.numerator() << '/' << C.denominator()
              << " margin=" << margin.numerator() << '/' << margin.denominator()
              << '\n';
    std::exit(EXIT_FAILURE);
  }
}

void enumerate(
    const std::vector<Rational>& alphabet,
    int maximum_length,
    std::size_t first,
    std::vector<Rational>& ratios,
    Counters& counters
) {
  if (ratios.size() >= 3U) {
    ++counters.tails;
    check_tail(ratios, counters);
  }
  if (static_cast<int>(ratios.size()) == maximum_length) {
    return;
  }
  for (std::size_t index = first; index < alphabet.size(); ++index) {
    ratios.push_back(alphabet[index]);
    enumerate(alphabet, maximum_length, index, ratios, counters);
    ratios.pop_back();
  }
}

int parse_positive(const char* text, const char* name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value, &consumed);
  if (consumed != value.size() || parsed < 3L || parsed > 12L) {
    throw std::invalid_argument(std::string(name) + " must lie in [3,12]");
  }
  return static_cast<int>(parsed);
}

std::uint64_t parse_positive_u64(const char* text, const char* name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const unsigned long long parsed = std::stoull(value, &consumed, 10);
  if (consumed != value.size() || parsed == 0U) {
    throw std::invalid_argument(std::string(name) + " must be positive");
  }
  return static_cast<std::uint64_t>(parsed);
}

int parse_positive_int(const char* text, const char* name) {
  const std::uint64_t parsed = parse_positive_u64(text, name);
  if (parsed > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument(std::string(name) + " is too large");
  }
  return static_cast<int>(parsed);
}

std::uint64_t splitmix64(std::uint64_t& state) {
  state += UINT64_C(0x9e3779b97f4a7c15);
  std::uint64_t value = state;
  value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31U);
}

void random_search(
    std::uint64_t samples,
    int maximum_length,
    std::uint64_t maximum_numerator,
    std::uint64_t denominator
) {
  if (maximum_length < 3) {
    throw std::invalid_argument("maximum_length must be at least three");
  }
  std::uint64_t state = UINT64_C(0x7e8b6d4c1f2a3905);
  Counters counters;
  for (std::uint64_t sample = 0U; sample < samples; ++sample) {
    const int length = 3 + static_cast<int>(
        splitmix64(state)
        % static_cast<std::uint64_t>(maximum_length - 2)
    );
    std::vector<Rational> ratios;
    ratios.reserve(static_cast<std::size_t>(length));
    for (int index = 0; index < length; ++index) {
      ratios.emplace_back(
          Integer(1) + Integer(splitmix64(state) % maximum_numerator),
          Integer(denominator)
      );
    }
    std::sort(ratios.begin(), ratios.end(), std::greater<Rational>());
    ++counters.tails;
    check_tail(ratios, counters);
  }
  std::cout << "SU2_WALL_121_RAT2_RANDOM_PASS samples=" << samples
            << " maximum_length=" << maximum_length
            << " maximum_numerator=" << maximum_numerator
            << " denominator=" << denominator
            << " seed=0x7e8b6d4c1f2a3905"
            << " tails=" << counters.tails
            << " critical=" << counters.critical << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (
        argc == 2
        && std::string(argv[1]) == "--replay-unit-floor-payment-obstruction"
    ) {
      Counters counters;
      check_tail(
          {Rational(4), Rational(4), Rational(4)}, counters, true);
      throw std::runtime_error("unit-floor obstruction replay did not exit");
    }
    if (argc == 6 && std::string(argv[1]) == "--random") {
      const std::uint64_t samples = parse_positive_u64(argv[2], "samples");
      const int maximum_length =
          parse_positive_int(argv[3], "maximum_length");
      const std::uint64_t maximum_numerator =
          parse_positive_u64(argv[4], "maximum_numerator");
      const std::uint64_t denominator =
          parse_positive_u64(argv[5], "denominator");
      random_search(samples, maximum_length, maximum_numerator, denominator);
      return EXIT_SUCCESS;
    }
    const int maximum_length = argc == 2
        ? parse_positive(argv[1], "maximum_length") : 8;
    if (argc > 2) {
      throw std::invalid_argument(
          "usage: probe_su2_wall_121_rat2 [maximum_length]"
          " | --random SAMPLES MAXIMUM_LENGTH MAXIMUM_NUMERATOR DENOMINATOR"
          " | --replay-unit-floor-payment-obstruction"
      );
    }
    // Descending order makes recursive index choices exactly the
    // nonincreasing-ratio log-concavity condition.
    const std::vector<Rational> alphabet{
        Rational(4), Rational(3), Rational(2), Rational(4, 3),
        Rational(1), Rational(2, 3), Rational(1, 2),
        Rational(1, 3), Rational(1, 4)
    };
    Counters counters;
    std::vector<Rational> ratios;
    enumerate(alphabet, maximum_length, 0U, ratios, counters);
    std::cout << "SU2_WALL_121_RAT2_GRID_PASS maximum_length="
              << maximum_length
              << " ratio_alphabet=4,3,2,4/3,1,2/3,1/2,1/3,1/4"
              << " tails=" << counters.tails
              << " critical=" << counters.critical
              << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "SU2_WALL_121_RAT2_ERROR " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
