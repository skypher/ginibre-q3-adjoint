#include <array>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;

int parse_positive(const char* text, const std::string& name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value, &consumed, 10);
  if (consumed != value.size() || parsed <= 0) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  return static_cast<int>(parsed);
}

Rational terminal_current(const Rational& r, const Rational& s,
                          const Rational& t) {
  std::array<Rational, 8> p{};
  p[1] = Rational(1);
  p[2] = r;
  p[3] = r * s;
  p[4] = r * s * t;

  std::array<Rational, 8> g{};
  std::array<Rational, 8> h{};
  std::array<Rational, 8> k{};
  for (int index = 1; index <= 5; ++index) {
    g[static_cast<std::size_t>(index)] =
        p[static_cast<std::size_t>(index)] *
            p[static_cast<std::size_t>(index)] -
        p[static_cast<std::size_t>(index - 1)] *
            p[static_cast<std::size_t>(index + 1)];
  }
  for (int index = 1; index <= 6; ++index) {
    h[static_cast<std::size_t>(index)] =
        p[static_cast<std::size_t>(index)] *
        p[static_cast<std::size_t>(index - 1)];
    if (index >= 2) {
      h[static_cast<std::size_t>(index)] -=
          p[static_cast<std::size_t>(index - 2)] *
          p[static_cast<std::size_t>(index + 1)];
    }
  }
  for (int index = 1; index <= 5; ++index) {
    k[static_cast<std::size_t>(index)] =
        g[static_cast<std::size_t>(index)] +
        h[static_cast<std::size_t>(index)] +
        h[static_cast<std::size_t>(index + 1)];
  }

  Rational result = p[1] * p[1] * p[1] * p[2];
  for (int index = 1; index <= 4; ++index) {
    result += (g[static_cast<std::size_t>(index)] -
               g[static_cast<std::size_t>(index + 1)]) *
              (k[static_cast<std::size_t>(index)] -
               k[static_cast<std::size_t>(index + 1)]);
  }
  return result;
}

Rational renewal(const Rational& r, const Rational& s, const Rational& t,
                 const Rational& q) {
  const Rational r2 = r * r;
  const Rational r3 = r2 * r;
  const Rational r4 = r2 * r2;
  const Rational s2 = s * s;
  const Rational s3 = s2 * s;
  const Rational t2 = t * t;
  const Rational t3 = t2 * t;
  const Rational polynomial =
      Rational(1) + r * (Rational(1) + Rational(2) * s + s * t) +
      r2 * (Rational(-2) - Rational(2) * s + Rational(2) * s2 + s2 * t) +
      r3 * (Rational(1) - Rational(4) * s - Rational(2) * s2 -
            Rational(2) * s2 * t + s2 * t2 + Rational(2) * s3 +
            s3 * t3) +
      r4;
  return polynomial - r3 * s2 * t * (Rational(1) + s * t) * q;
}

Rational two_term_terminal(const Rational& value) {
  const Rational square = value * value;
  return Rational(1) + value - Rational(2) * square + square * value +
         Rational(2) * square * square;
}

std::string render(const Rational& value) {
  return value.numerator().convert_to<std::string>() + "/" +
         value.denominator().convert_to<std::string>();
}

int run_grid(const int denominator) {
  std::uint64_t terminal_checks = 0U;
  std::uint64_t renewal_checks = 0U;
  for (int numerator_r = 0; numerator_r <= denominator; ++numerator_r) {
    const Rational r{Integer(numerator_r), Integer(denominator)};
    for (int numerator_s = 0; numerator_s <= numerator_r; ++numerator_s) {
      const Rational s{Integer(numerator_s), Integer(denominator)};
      for (int numerator_t = 0; numerator_t <= numerator_s; ++numerator_t) {
        const Rational t{Integer(numerator_t), Integer(denominator)};
        const Rational terminal = terminal_current(r, s, t);
        ++terminal_checks;
        const Rational barrier = two_term_terminal(t);
        if (terminal < barrier) {
          std::cout << "SU2_WALL_121_SIMPLE_BARRIER_TERMINAL_FAILURE"
                    << " r=" << render(r) << " s=" << render(s)
                    << " t=" << render(t)
                    << " value=" << render(terminal - barrier)
                    << '\n';
          return EXIT_FAILURE;
        }
        const Rational coefficient =
            r - s * s * t * (Rational(1) + s * t);
        const Rational q = coefficient < 0 ? t : Rational(0);
        const Rational margin =
            renewal(r, s, t, q) + r * r * r * r * two_term_terminal(q) -
            barrier;
        ++renewal_checks;
        if (margin < 0) {
          std::cout << "SU2_WALL_121_SIMPLE_BARRIER_RENEWAL_FAILURE"
                    << " r=" << render(r) << " s=" << render(s)
                    << " t=" << render(t) << " q=" << render(q)
                    << " value=" << render(margin) << '\n';
          return EXIT_FAILURE;
        }
      }
    }
  }
  std::cout << "SU2_WALL_121_SIMPLE_BARRIER_GRID_PASS"
            << " denominator=" << denominator
            << " terminal_checks=" << terminal_checks
            << " renewal_checks=" << renewal_checks << '\n';
  return EXIT_SUCCESS;
}

int replay_last_ratio_obstructions() {
  const Rational small{1, 64};
  const Rational small_margin =
      renewal(small, small, small, Rational(0)) +
      small * small * small * small * Rational(1) -
      (Rational(1) + small);
  const Rational small_expected{-16773119, 18014398509481984LL};
  if (small_margin != small_expected || small_margin >= 0) {
    throw std::runtime_error("linear last-ratio barrier replay mismatch");
  }

  const Rational rho{25, 32};
  const Rational quartic_margin =
      renewal(rho, rho, rho, Rational(0)) +
      rho * rho * rho * rho * two_term_terminal(Rational(0)) -
      two_term_terminal(rho);
  const Rational quartic_expected{-217245774375LL, 35184372088832LL};
  if (quartic_margin != quartic_expected || quartic_margin >= 0) {
    throw std::runtime_error("quartic last-ratio barrier replay mismatch");
  }

  std::cout << "SU2_WALL_121_LAST_RATIO_BARRIER_OBSTRUCTIONS"
            << " linear_rho=1/64 linear_margin=" << render(small_margin)
            << " quartic_rho=25/32 quartic_margin="
            << render(quartic_margin) << " result=PASS_EXACT\n";
  return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string{argv[1]} == "--replay-last-ratio") {
      return replay_last_ratio_obstructions();
    }
    const int denominator =
        argc == 2 ? parse_positive(argv[1], "denominator") : 24;
    if (argc > 2) {
      throw std::invalid_argument(
          "usage: probe_su2_wall_121_simple_barrier [denominator]");
    }
    return run_grid(denominator);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
