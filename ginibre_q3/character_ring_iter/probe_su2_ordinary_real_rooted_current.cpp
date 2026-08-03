#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using Polynomial = std::vector<Rational>;

int parse_nonnegative(const char* text, const std::string& name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value, &consumed, 10);
  if (consumed != value.size() || parsed < 0) {
    throw std::invalid_argument(name + " must be a nonnegative integer");
  }
  return static_cast<int>(parsed);
}

int parse_positive(const char* text, const std::string& name) {
  const int result = parse_nonnegative(text, name);
  if (result == 0) {
    throw std::invalid_argument(name + " must be positive");
  }
  return result;
}

void trim(Polynomial& polynomial) {
  while (!polynomial.empty() && polynomial.back() == Rational(0)) {
    polynomial.pop_back();
  }
}

int degree(const Polynomial& polynomial) {
  return static_cast<int>(polynomial.size()) - 1;
}

Polynomial add(const Polynomial& left, const Polynomial& right,
               const Rational scale_right = Rational(1)) {
  Polynomial result(std::max(left.size(), right.size()), Rational(0));
  for (std::size_t index = 0U; index < left.size(); ++index) {
    result[index] += left[index];
  }
  for (std::size_t index = 0U; index < right.size(); ++index) {
    result[index] += scale_right * right[index];
  }
  trim(result);
  return result;
}

Polynomial multiply(const Polynomial& left, const Polynomial& right) {
  if (left.empty() || right.empty()) {
    return {};
  }
  Polynomial result(left.size() + right.size() - 1U, Rational(0));
  for (std::size_t first = 0U; first < left.size(); ++first) {
    for (std::size_t second = 0U; second < right.size(); ++second) {
      result[first + second] += left[first] * right[second];
    }
  }
  trim(result);
  return result;
}

Polynomial derivative(const Polynomial& polynomial) {
  if (polynomial.size() <= 1U) {
    return {};
  }
  Polynomial result(polynomial.size() - 1U, Rational(0));
  for (std::size_t index = 1U; index < polynomial.size(); ++index) {
    result[index - 1U] =
        Rational(static_cast<unsigned long>(index)) * polynomial[index];
  }
  trim(result);
  return result;
}

std::pair<Polynomial, Polynomial> divide_with_remainder(
    Polynomial numerator, const Polynomial& denominator) {
  if (denominator.empty()) {
    throw std::invalid_argument("division by the zero polynomial");
  }
  Polynomial quotient(
      static_cast<std::size_t>(std::max(0, degree(numerator) -
                                            degree(denominator) + 1)),
      Rational(0));
  while (!numerator.empty() && degree(numerator) >= degree(denominator)) {
    const int shift = degree(numerator) - degree(denominator);
    const Rational factor = numerator.back() / denominator.back();
    quotient[static_cast<std::size_t>(shift)] += factor;
    for (int index = 0; index <= degree(denominator); ++index) {
      numerator[static_cast<std::size_t>(index + shift)] -=
          factor * denominator[static_cast<std::size_t>(index)];
    }
    trim(numerator);
  }
  trim(quotient);
  return {std::move(quotient), std::move(numerator)};
}

Polynomial monic(Polynomial polynomial) {
  if (polynomial.empty()) {
    return polynomial;
  }
  const Rational leading = polynomial.back();
  for (Rational& coefficient : polynomial) {
    coefficient /= leading;
  }
  return polynomial;
}

Polynomial polynomial_gcd(Polynomial left, Polynomial right) {
  while (!right.empty()) {
    Polynomial remainder = divide_with_remainder(left, right).second;
    left = std::move(right);
    right = std::move(remainder);
  }
  return monic(std::move(left));
}

Rational evaluate(const Polynomial& polynomial, const int argument) {
  Rational value(0);
  for (auto iterator = polynomial.rbegin(); iterator != polynomial.rend();
       ++iterator) {
    value = value * Rational(argument) + *iterator;
  }
  return value;
}

int sign(const Rational& value) {
  return value > Rational(0) ? 1 : (value < Rational(0) ? -1 : 0);
}

int sturm_variations(const std::vector<Polynomial>& sequence,
                     const int argument) {
  int previous = 0;
  int variations = 0;
  for (const Polynomial& polynomial : sequence) {
    const int current = sign(evaluate(polynomial, argument));
    if (current == 0) {
      continue;
    }
    if (previous != 0 && previous != current) {
      ++variations;
    }
    previous = current;
  }
  return variations;
}

bool has_only_interior_real_roots(const Polynomial& polynomial) {
  if (degree(polynomial) <= 0) {
    return true;
  }
  if (evaluate(polynomial, -1) == Rational(0) ||
      evaluate(polynomial, 3) == Rational(0)) {
    return false;
  }
  const Polynomial gcd = polynomial_gcd(polynomial, derivative(polynomial));
  const auto [square_free, remainder] = divide_with_remainder(polynomial, gcd);
  if (!remainder.empty()) {
    throw std::runtime_error("nonexact square-free quotient");
  }
  std::vector<Polynomial> sturm{square_free, derivative(square_free)};
  while (!sturm.back().empty()) {
    Polynomial remainder_value =
        divide_with_remainder(sturm[sturm.size() - 2U], sturm.back()).second;
    for (Rational& coefficient : remainder_value) {
      coefficient = -coefficient;
    }
    trim(remainder_value);
    if (remainder_value.empty()) {
      break;
    }
    sturm.push_back(std::move(remainder_value));
  }
  const int roots = sturm_variations(sturm, -1) - sturm_variations(sturm, 3);
  return roots == degree(square_free);
}

std::vector<Polynomial> beta_basis(const int maximum_degree) {
  std::vector<Polynomial> result(static_cast<std::size_t>(maximum_degree + 1));
  result[0] = Polynomial{Rational(1)};
  if (maximum_degree == 0) {
    return result;
  }
  result[1] = Polynomial{Rational(0), Rational(1)};
  const Polynomial x_minus_one{Rational(-1), Rational(1)};
  for (int index = 1; index < maximum_degree; ++index) {
    result[static_cast<std::size_t>(index + 1)] =
        add(multiply(x_minus_one, result[static_cast<std::size_t>(index)]),
            result[static_cast<std::size_t>(index - 1)], Rational(-1));
  }
  return result;
}

Polynomial character_polynomial(const std::vector<int>& profile,
                                const std::vector<Polynomial>& basis) {
  Polynomial result;
  for (std::size_t index = 0U; index < profile.size(); ++index) {
    if (profile[index] == 0) {
      continue;
    }
    Polynomial term = basis[index];
    for (Rational& coefficient : term) {
      coefficient *= Rational(profile[index]);
    }
    result = add(result, term);
  }
  return result;
}

using Profile = std::vector<Integer>;

Integer coefficient(const Profile& profile, const int index) {
  return index >= 0 && index < static_cast<int>(profile.size())
             ? profile[static_cast<std::size_t>(index)]
             : Integer(0);
}

Profile square_character(const std::vector<int>& root) {
  const int support = static_cast<int>(root.size()) - 1;
  Profile result(static_cast<std::size_t>(2 * support + 1), Integer(0));
  for (int first = 0; first <= support; ++first) {
    for (int second = 0; second <= support; ++second) {
      const Integer weight =
          Integer(root[static_cast<std::size_t>(first)]) *
          Integer(root[static_cast<std::size_t>(second)]);
      for (int label = std::abs(first - second); label <= first + second;
           ++label) {
        result[static_cast<std::size_t>(label)] += weight;
      }
    }
  }
  return result;
}

Integer current(const Profile& square, const int radius, const int target) {
  Integer interval = 0;
  for (int label = std::abs(radius - target); label <= radius + target;
       ++label) {
    interval += coefficient(square, label);
  }
  return square.front() * interval - coefficient(square, radius) *
         coefficient(square, target);
}

std::string render(const std::vector<int>& profile) {
  std::string result = "[";
  for (std::size_t index = 0U; index < profile.size(); ++index) {
    if (index != 0U) {
      result += ',';
    }
    result += std::to_string(profile[index]);
  }
  return result + ']';
}

struct Counters {
  std::uint64_t examined = 0U;
  std::uint64_t shaped = 0U;
  std::uint64_t real_rooted = 0U;
  std::uint64_t currents = 0U;
  std::uint64_t failures = 0U;
  std::vector<int> first_failure_profile;
  int first_failure_radius = -1;
  int first_failure_target = -1;
  Integer first_failure_value = 0;
};

bool has_interval_log_concave_support(const std::vector<int>& profile) {
  int first = -1;
  int last = -1;
  for (int index = 0; index < static_cast<int>(profile.size()); ++index) {
    if (profile[static_cast<std::size_t>(index)] > 0) {
      if (first < 0) {
        first = index;
      }
      last = index;
    }
  }
  if (first < 0) {
    return false;
  }
  for (int index = first; index <= last; ++index) {
    if (profile[static_cast<std::size_t>(index)] == 0) {
      return false;
    }
  }
  for (int index = first + 1; index < last; ++index) {
    const std::int64_t middle = profile[static_cast<std::size_t>(index)];
    if (middle * middle <
        static_cast<std::int64_t>(profile[static_cast<std::size_t>(index - 1)]) *
            profile[static_cast<std::size_t>(index + 1)]) {
      return false;
    }
  }
  return true;
}

void inspect(const std::vector<int>& profile,
             const std::vector<Polynomial>& basis, Counters& counters) {
  ++counters.examined;
  if (!has_interval_log_concave_support(profile)) {
    return;
  }
  ++counters.shaped;
  if (!has_only_interior_real_roots(character_polynomial(profile, basis))) {
    return;
  }
  ++counters.real_rooted;
  const Profile square = square_character(profile);
  const int support = static_cast<int>(square.size()) - 1;
  for (int radius = 0; radius <= support; ++radius) {
    for (int target = radius; target <= support; ++target) {
      const Integer margin = current(square, radius, target);
      ++counters.currents;
      if (margin < 0) {
        ++counters.failures;
        if (counters.first_failure_profile.empty()) {
          counters.first_failure_profile = profile;
          counters.first_failure_radius = radius;
          counters.first_failure_target = target;
          counters.first_failure_value = margin;
        }
      }
    }
  }
}

void enumerate(const int maximum_degree, const int maximum_coefficient,
               const std::vector<Polynomial>& basis, const int index,
               std::vector<int>& profile, Counters& counters) {
  if (index == maximum_degree) {
    for (int value = 1; value <= maximum_coefficient; ++value) {
      profile[static_cast<std::size_t>(index)] = value;
      inspect(profile, basis, counters);
    }
    return;
  }
  for (int value = 0; value <= maximum_coefficient; ++value) {
    profile[static_cast<std::size_t>(index)] = value;
    enumerate(maximum_degree, maximum_coefficient, basis, index + 1, profile,
              counters);
  }
}

void replay_broad_real_rooted_obstruction() {
  const std::vector<int> profile{0, 1, 2, 2, 2, 0, 2};
  const std::vector<Polynomial> basis = beta_basis(6);
  if (has_interval_log_concave_support(profile)) {
    throw std::runtime_error("replay profile unexpectedly has interval support");
  }
  if (!has_only_interior_real_roots(character_polynomial(profile, basis))) {
    throw std::runtime_error("replay profile is not interior real rooted");
  }
  const Integer margin = current(square_character(profile), 1, 12);
  if (margin != -12) {
    throw std::runtime_error("real-rooted replay current mismatch");
  }
  std::cout << "SU2_ORDINARY_REAL_ROOTED_BROAD_OBSTRUCTION"
            << " profile=" << render(profile)
            << " interior_real_rooted=1"
            << " interval_log_concave=0"
            << " R=1 S=12 current=" << margin
            << " result=PASS_EXACT\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--replay-broad-obstruction") {
      replay_broad_real_rooted_obstruction();
      return EXIT_SUCCESS;
    }
    const int maximum_degree =
        argc >= 2 ? parse_nonnegative(argv[1], "maximum degree") : 6;
    const int maximum_coefficient =
        argc >= 3 ? parse_positive(argv[2], "maximum coefficient") : 4;
    Counters counters;
    const std::vector<Polynomial> basis = beta_basis(maximum_degree);
    for (int degree_value = 0; degree_value <= maximum_degree; ++degree_value) {
      std::vector<int> profile(static_cast<std::size_t>(degree_value + 1), 0);
      enumerate(degree_value, maximum_coefficient, basis, 0, profile, counters);
    }
    std::cout << "SU2_ORDINARY_REAL_ROOTED_CURRENT"
              << " maximum_degree=" << maximum_degree
              << " maximum_coefficient=" << maximum_coefficient
              << " examined=" << counters.examined
              << " shaped=" << counters.shaped
              << " real_rooted=" << counters.real_rooted
              << " currents=" << counters.currents
              << " failures=" << counters.failures;
    if (!counters.first_failure_profile.empty()) {
      std::cout << " first_profile=" << render(counters.first_failure_profile)
                << " first_R=" << counters.first_failure_radius
                << " first_S=" << counters.first_failure_target
                << " first_value=" << counters.first_failure_value;
    }
    std::cout << " result="
              << (counters.failures == 0U ? "NO_COUNTEREXAMPLE" : "FAIL")
              << '\n';
    return counters.failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
