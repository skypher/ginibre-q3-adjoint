#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// This is a falsification probe for the remaining wall-(1,2,1) Bellman
// payment.  It does not certify a sign: every reported minimum is recomputed
// from the displayed rational grid point before it is used as a proof lead.

namespace {

int parse_positive(const char* text, const std::string& name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value, &consumed, 10);
  if (consumed != value.size() || parsed <= 0L) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  return static_cast<int>(parsed);
}

long double terminal_current(const long double r, const long double s,
                             const long double t) {
  const std::vector<long double> p{0.0L, 1.0L, r, r * s, r * s * t,
                                   0.0L, 0.0L, 0.0L};
  std::vector<long double> g(8U, 0.0L);
  std::vector<long double> h(8U, 0.0L);
  std::vector<long double> k(8U, 0.0L);
  for (int index = 1; index <= 5; ++index) {
    g[static_cast<std::size_t>(index)] =
        p[static_cast<std::size_t>(index)]
            * p[static_cast<std::size_t>(index)]
        - p[static_cast<std::size_t>(index - 1)]
            * p[static_cast<std::size_t>(index + 1)];
  }
  for (int index = 1; index <= 6; ++index) {
    h[static_cast<std::size_t>(index)] =
        p[static_cast<std::size_t>(index)]
            * p[static_cast<std::size_t>(index - 1)];
    if (index >= 2) {
      h[static_cast<std::size_t>(index)] -=
          p[static_cast<std::size_t>(index - 2)]
              * p[static_cast<std::size_t>(index + 1)];
    }
  }
  for (int index = 1; index <= 5; ++index) {
    k[static_cast<std::size_t>(index)] =
        g[static_cast<std::size_t>(index)]
        + h[static_cast<std::size_t>(index)]
        + h[static_cast<std::size_t>(index + 1)];
  }
  long double result = r;
  for (int index = 1; index <= 4; ++index) {
    result += (g[static_cast<std::size_t>(index)]
                   - g[static_cast<std::size_t>(index + 1)])
              * (k[static_cast<std::size_t>(index)]
                     - k[static_cast<std::size_t>(index + 1)]);
  }
  return result;
}

long double renewal(const long double r, const long double s,
                    const long double t, const long double q) {
  const long double p =
      1.0L + r * (1.0L + 2.0L * s + s * t)
      + r * r * (-2.0L - 2.0L * s + 2.0L * s * s + s * s * t)
      + r * r * r
            * (1.0L - 4.0L * s - 2.0L * s * s - 2.0L * s * s * t
               + s * s * t * t + 2.0L * s * s * s
               + s * s * s * t * t * t)
      + r * r * r * r;
  return p - r * r * r * s * s * t * (1.0L + s * t) * q;
}

long double append_delta(const long double r, const long double s,
                         const long double t, const long double q) {
  const long double b = r;
  const long double c = r * s;
  const long double d = c * t;
  const long double x = d * q;
  const long double linear =
      d * d * d - c * d - 4.0L * c * d * d - 2.0L * c * c * d
      + 2.0L * c * c * c - b * c - 2.0L * b * c * d + b * b * d
      + b * b * b;
  const long double quadratic =
      -2.0L * d * d - 2.0L * c * d + 2.0L * c * c + b * c;
  return x * (linear + quadratic * x + (b + 2.0L * c + d) * x * x
              + 2.0L * x * x * x);
}

long double geometric_barrier(const long double ratio) {
  if (ratio <= 1.0L) {
    return 1.0L + ratio;
  }
  return terminal_current(ratio, ratio, ratio);
}

long double saturated_barrier(const long double r, const long double s,
                              const long double t, const int depth) {
  long double value = geometric_barrier(t);
  long double first = r;
  long double second = s;
  long double third = t;
  for (int step = 0; step < depth; ++step) {
    value = std::min(
        terminal_current(first, second, third),
        renewal(first, second, third, third)
            + first * first * first * first * value);
    first = second;
    second = third;
  }
  return value;
}

long double clipped_geometric_barrier(const long double r,
                                      const long double s,
                                      const long double t) {
  return std::min(terminal_current(r, s, t), 1.0L + r);
}

long double clipped_last_barrier(const long double r, const long double s,
                                 const long double t) {
  static_cast<void>(r);
  static_cast<void>(s);
  return std::min(terminal_current(r, s, t), 1.0L + t);
}

long double clipped_tangent_barrier(const long double r, const long double s,
                                    const long double t) {
  const long double tangent =
      1.0L + t + (1.0L - 2.0L * t) * (r - s)
      + (1.0L - t * t - t * t * t * t) * (s - t);
  return std::min(terminal_current(r, s, t), tangent);
}

long double clipped_shift_barrier(const long double r, const long double s,
                                  const long double t) {
  static_cast<void>(s);
  const long double shift_compatible =
      1.0L + t + (1.0L - 2.0L * t) * (r - t);
  return std::min(terminal_current(r, s, t), shift_compatible);
}

long double clipped_downward_barrier(const long double r,
                                     const long double s,
                                     const long double t) {
  static_cast<void>(s);
  return std::min(terminal_current(r, s, t), 1.0L + 2.0L * t - r);
}

long double clipped_steep_barrier(const long double r, const long double s,
                                  const long double t) {
  static_cast<void>(s);
  const long double slope = 1.0L - 2.5L * t;
  return std::min(
      terminal_current(r, s, t), 1.0L + t + slope * (r - t));
}

long double clipped_quadratic_barrier(const long double r,
                                      const long double s,
                                      const long double t) {
  static_cast<void>(s);
  const long double slope = 1.0L - 5.0L * t * t;
  return std::min(
      terminal_current(r, s, t), 1.0L + t + slope * (r - t));
}

bool critical_demand(const long double r, const long double s,
                     const long double t, long double& demand) {
  if (r == 0.0L) {
    return false;
  }
  const long double c = r * s;
  const long double d = c * t;
  const long double a_coefficient = 1.0L + r * r + (1.0L + r) * c;
  const long double c_coefficient = 1.0L + r + c;
  const long double b_coefficient =
      3.0L * r + 2.0L * r * r + r * c + r * d + r * c * d
      - r * r * r - c * c - c * c * c;
  if (b_coefficient <= 0.0L) {
    return false;
  }
  const long double endpoint_derivative =
      -b_coefficient + 2.0L * a_coefficient / r
      + 3.0L * c_coefficient / (r * r);
  if (endpoint_derivative <= 0.0L) {
    return false;
  }
  const long double discriminant =
      a_coefficient * a_coefficient
      + 3.0L * b_coefficient * c_coefficient;
  const long double tau =
      (-a_coefficient + std::sqrt(discriminant))
      / (3.0L * c_coefficient);
  demand = a_coefficient * tau * tau
           + 2.0L * c_coefficient * tau * tau * tau;
  return true;
}

struct Minimum {
  long double value = std::numeric_limits<long double>::infinity();
  int r = 0;
  int s = 0;
  int t = 0;
  int q = 0;
};

void update(Minimum& minimum, const long double value, const int r,
            const int s, const int t, const int q) {
  if (value < minimum.value) {
    minimum.value = value;
    minimum.r = r;
    minimum.s = s;
    minimum.t = t;
    minimum.q = q;
  }
}

void print_minimum(const std::string& name, const Minimum& minimum,
                   const int denominator) {
  std::cout << ' ' << name << '=' << std::setprecision(18)
            << minimum.value << "@(" << minimum.r << '/' << denominator
            << ',' << minimum.s << '/' << denominator << ',' << minimum.t
            << '/' << denominator;
  if (minimum.q >= 0) {
    std::cout << ',' << minimum.q << '/' << denominator;
  }
  std::cout << ')';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 3 || argc > 5) {
      throw std::invalid_argument(
          "usage: probe_su2_wall_121_bellman_barrier MAX_RATIO_NUMERATOR "
          "DENOMINATOR [BELLMAN_DEPTH] "
          "[--candidate-one|--candidate-r|--candidate-saturated|"
          "--candidate-saturated-deep|--candidate-tangent|"
          "--candidate-shift|--candidate-downward|--candidate-last|"
          "--candidate-steep|--candidate-quadratic]");
    }
    const int maximum = parse_positive(argv[1], "MAX_RATIO_NUMERATOR");
    const int denominator = parse_positive(argv[2], "DENOMINATOR");
    int depth = 8;
    bool candidate_one = false;
    bool candidate_r = false;
    bool candidate_saturated = false;
    bool candidate_saturated_deep = false;
    bool candidate_tangent = false;
    bool candidate_shift = false;
    bool candidate_downward = false;
    bool candidate_last = false;
    bool candidate_steep = false;
    bool candidate_quadratic = false;
    bool candidate_clipped = false;
    if (argc >= 4) {
      if (std::string{argv[3]} == "--candidate-one") {
        candidate_one = true;
      } else if (std::string{argv[3]} == "--candidate-r") {
        candidate_r = true;
      } else if (std::string{argv[3]} == "--candidate-saturated") {
        candidate_saturated = true;
      } else if (std::string{argv[3]} == "--candidate-saturated-deep") {
        candidate_saturated_deep = true;
      } else if (std::string{argv[3]} == "--candidate-tangent") {
        candidate_tangent = true;
      } else if (std::string{argv[3]} == "--candidate-shift") {
        candidate_shift = true;
      } else if (std::string{argv[3]} == "--candidate-downward") {
        candidate_downward = true;
      } else if (std::string{argv[3]} == "--candidate-last") {
        candidate_last = true;
      } else if (std::string{argv[3]} == "--candidate-steep") {
        candidate_steep = true;
      } else if (std::string{argv[3]} == "--candidate-quadratic") {
        candidate_quadratic = true;
      } else if (std::string{argv[3]} == "--candidate-clipped") {
        candidate_clipped = true;
      } else {
        depth = parse_positive(argv[3], "BELLMAN_DEPTH");
      }
    }
    if (argc == 5) {
      if (std::string{argv[4]} == "--candidate-one") {
        candidate_one = true;
      } else if (std::string{argv[4]} == "--candidate-r") {
        candidate_r = true;
      } else if (std::string{argv[4]} == "--candidate-saturated") {
        candidate_saturated = true;
      } else if (std::string{argv[4]} == "--candidate-saturated-deep") {
        candidate_saturated_deep = true;
      } else if (std::string{argv[4]} == "--candidate-tangent") {
        candidate_tangent = true;
      } else if (std::string{argv[4]} == "--candidate-shift") {
        candidate_shift = true;
      } else if (std::string{argv[4]} == "--candidate-downward") {
        candidate_downward = true;
      } else if (std::string{argv[4]} == "--candidate-last") {
        candidate_last = true;
      } else if (std::string{argv[4]} == "--candidate-steep") {
        candidate_steep = true;
      } else if (std::string{argv[4]} == "--candidate-quadratic") {
        candidate_quadratic = true;
      } else if (std::string{argv[4]} == "--candidate-clipped") {
        candidate_clipped = true;
      } else {
        throw std::invalid_argument("unknown barrier candidate");
      }
    }
    if (static_cast<int>(candidate_one) + static_cast<int>(candidate_r)
            + static_cast<int>(candidate_saturated)
            + static_cast<int>(candidate_saturated_deep)
            + static_cast<int>(candidate_tangent)
            + static_cast<int>(candidate_shift)
            + static_cast<int>(candidate_downward)
            + static_cast<int>(candidate_last)
            + static_cast<int>(candidate_steep)
            + static_cast<int>(candidate_quadratic)
            + static_cast<int>(candidate_clipped)
        > 1) {
      throw std::invalid_argument("choose only one barrier candidate");
    }
    const auto barrier = [candidate_clipped, candidate_one, candidate_r,
                          candidate_saturated,
                          candidate_saturated_deep,
                          candidate_tangent,
                          candidate_shift,
                          candidate_downward,
                          candidate_last,
                          candidate_steep,
                          candidate_quadratic](const long double r,
                                               const long double s,
                                               const long double t) {
      if (candidate_one) {
        return 1.0L;
      }
      if (candidate_clipped) {
        return clipped_geometric_barrier(r, s, t);
      }
      if (candidate_saturated) {
        return saturated_barrier(r, s, t, 2);
      }
      if (candidate_saturated_deep) {
        return saturated_barrier(r, s, t, 8);
      }
      if (candidate_tangent) {
        return clipped_tangent_barrier(r, s, t);
      }
      if (candidate_shift) {
        return clipped_shift_barrier(r, s, t);
      }
      if (candidate_downward) {
        return clipped_downward_barrier(r, s, t);
      }
      if (candidate_last) {
        return clipped_last_barrier(r, s, t);
      }
      if (candidate_steep) {
        return clipped_steep_barrier(r, s, t);
      }
      if (candidate_quadratic) {
        return clipped_quadratic_barrier(r, s, t);
      }
      return 1.0L + (candidate_r ? r : t);
    };
    const std::size_t side = static_cast<std::size_t>(maximum + 1);
    const std::size_t entries = side * side * side;
    const auto flat_index = [side](const int r, const int s, const int t) {
      return (static_cast<std::size_t>(r) * side
              + static_cast<std::size_t>(s))
                 * side
             + static_cast<std::size_t>(t);
    };
    Minimum terminal;
    Minimum terminal_floor;
    Minimum renewal_barrier;
    Minimum critical_renewal_barrier;
    Minimum append_above_one;
    Minimum terminal_payment;
    Minimum barrier_payment;
    Minimum renewal_payment;
    Minimum critical_renewal_payment;
    Minimum critical_to_noncritical_payment;
    Minimum critical_to_nongeometric_payment;
    Minimum critical_unit_floor_payment;
    Minimum critical_unit_floor_q_le_one_payment;
    Minimum critical_unit_floor_q_lt_one_payment;
    Minimum bellman_margin;
    std::vector<long double> bellman(entries, 0.0L);
    std::vector<std::vector<int>> policies;
    std::uint64_t states = 0U;
    for (int r_index = 0; r_index <= maximum; ++r_index) {
      const long double r = static_cast<long double>(r_index)
                            / static_cast<long double>(denominator);
      for (int s_index = 0; s_index <= r_index; ++s_index) {
        const long double s = static_cast<long double>(s_index)
                              / static_cast<long double>(denominator);
        for (int t_index = 0; t_index <= s_index; ++t_index) {
          const long double t = static_cast<long double>(t_index)
                                / static_cast<long double>(denominator);
          const long double terminal_residual =
              terminal_current(r, s, t) - barrier(r, s, t);
          update(terminal, terminal_residual, r_index, s_index, t_index,
                 -1);
          update(terminal_floor, terminal_current(r, s, t) - 1.0L,
                 r_index, s_index, t_index, -1);
          bellman[flat_index(r_index, s_index, t_index)] =
              terminal_current(r, s, t);
          long double local_demand = 0.0L;
          const bool critical = critical_demand(r, s, t, local_demand);
          update(
              terminal_payment,
              terminal_current(r, s, t) - (critical ? local_demand : 0.0L),
              r_index, s_index, t_index, -1
          );
          update(
              barrier_payment,
              barrier(r, s, t) - (critical ? local_demand : 0.0L),
              r_index, s_index, t_index, -1
          );
          for (int q_index = 0; q_index <= t_index; ++q_index) {
            const long double q = static_cast<long double>(q_index)
                                  / static_cast<long double>(denominator);
            const long double barrier_residual =
                renewal(r, s, t, q)
                + std::pow(r, 4) * barrier(s, t, q) - barrier(r, s, t);
            update(renewal_barrier, barrier_residual, r_index, s_index,
                   t_index, q_index);
            if (critical) {
              update(critical_renewal_barrier, barrier_residual, r_index,
                     s_index, t_index, q_index);
            }
            long double next_demand = 0.0L;
            const bool next_critical =
                critical_demand(s, t, q, next_demand);
            const long double payment_residual =
                renewal(r, s, t, q)
                + r * r * r * r * (next_critical ? next_demand : 0.0L)
                - (critical ? local_demand : 0.0L);
            update(renewal_payment, payment_residual, r_index, s_index,
                   t_index, q_index);
            if (critical) {
              update(critical_renewal_payment, payment_residual, r_index,
                     s_index, t_index, q_index);
              if (!next_critical) {
                update(critical_to_noncritical_payment, payment_residual,
                       r_index, s_index, t_index, q_index);
              } else if (s_index != t_index || t_index != q_index) {
                update(critical_to_nongeometric_payment, payment_residual,
                       r_index, s_index, t_index, q_index);
              }
              update(critical_unit_floor_payment,
                     renewal(r, s, t, q) + r * r * r * r - local_demand,
                     r_index, s_index, t_index, q_index);
              if (q_index <= denominator) {
                update(critical_unit_floor_q_le_one_payment,
                       renewal(r, s, t, q) + r * r * r * r - local_demand,
                       r_index, s_index, t_index, q_index);
              }
              if (q_index < denominator) {
                update(critical_unit_floor_q_lt_one_payment,
                       renewal(r, s, t, q) + r * r * r * r - local_demand,
                       r_index, s_index, t_index, q_index);
              }
            }
            if (q_index >= denominator) {
              update(append_above_one, append_delta(r, s, t, q), r_index,
                     s_index, t_index, q_index);
            }
            ++states;
          }
        }
      }
    }
    for (int iteration = 0; iteration < depth; ++iteration) {
      std::vector<long double> next = bellman;
      std::vector<int> policy(entries, -1);
      for (int r_index = 0; r_index <= maximum; ++r_index) {
        const long double r = static_cast<long double>(r_index)
                              / static_cast<long double>(denominator);
        const long double r_four = r * r * r * r;
        for (int s_index = 0; s_index <= r_index; ++s_index) {
          const long double s = static_cast<long double>(s_index)
                                / static_cast<long double>(denominator);
          for (int t_index = 0; t_index <= s_index; ++t_index) {
            const long double t = static_cast<long double>(t_index)
                                  / static_cast<long double>(denominator);
            long double best =
                bellman[flat_index(r_index, s_index, t_index)];
            for (int q_index = 0; q_index <= t_index; ++q_index) {
              const long double q = static_cast<long double>(q_index)
                                    / static_cast<long double>(denominator);
              best = std::min(
                  best,
                  renewal(r, s, t, q)
                      + r_four
                            * bellman[flat_index(s_index, t_index, q_index)]
              );
              const long double candidate =
                  renewal(r, s, t, q)
                  + r_four
                        * bellman[flat_index(s_index, t_index, q_index)];
              if (candidate < next[flat_index(r_index, s_index, t_index)]) {
                next[flat_index(r_index, s_index, t_index)] = candidate;
                policy[flat_index(r_index, s_index, t_index)] = q_index;
              }
            }
            next[flat_index(r_index, s_index, t_index)] = best;
          }
        }
      }
      bellman = std::move(next);
      policies.push_back(std::move(policy));
    }
    for (int r_index = 0; r_index <= maximum; ++r_index) {
      const long double r = static_cast<long double>(r_index)
                            / static_cast<long double>(denominator);
      for (int s_index = 0; s_index <= r_index; ++s_index) {
        const long double s = static_cast<long double>(s_index)
                              / static_cast<long double>(denominator);
        for (int t_index = 0; t_index <= s_index; ++t_index) {
          const long double t = static_cast<long double>(t_index)
                                / static_cast<long double>(denominator);
          long double demand = 0.0L;
          if (critical_demand(r, s, t, demand)) {
            update(
                bellman_margin,
                bellman[flat_index(r_index, s_index, t_index)] - demand,
                r_index, s_index, t_index, -1
            );
          }
        }
      }
    }
    // The per-iteration policies record only a strict improvement over the
    // preceding finite-horizon value.  Recompute the minimizing action from
    // the final horizon values before reading off the control structure.
    std::vector<int> final_policy(entries, -1);
    for (int r_index = 0; r_index <= maximum; ++r_index) {
      const long double r = static_cast<long double>(r_index)
                            / static_cast<long double>(denominator);
      const long double r_four = r * r * r * r;
      for (int s_index = 0; s_index <= r_index; ++s_index) {
        const long double s = static_cast<long double>(s_index)
                              / static_cast<long double>(denominator);
        for (int t_index = 0; t_index <= s_index; ++t_index) {
          const long double t = static_cast<long double>(t_index)
                                / static_cast<long double>(denominator);
          long double best = terminal_current(r, s, t);
          for (int q_index = 0; q_index <= t_index; ++q_index) {
            const long double q = static_cast<long double>(q_index)
                                  / static_cast<long double>(denominator);
            const long double candidate =
                renewal(r, s, t, q)
                + r_four * bellman[flat_index(s_index, t_index, q_index)];
            if (candidate < best) {
              best = candidate;
              final_policy[flat_index(r_index, s_index, t_index)] = q_index;
            }
          }
        }
      }
    }
    std::uint64_t policy_terminal = 0U;
    std::uint64_t policy_zero = 0U;
    std::uint64_t policy_upper = 0U;
    std::uint64_t policy_interior = 0U;
    std::uint64_t critical_policy_terminal = 0U;
    std::uint64_t critical_policy_zero = 0U;
    std::uint64_t critical_policy_upper = 0U;
    std::uint64_t critical_policy_interior = 0U;
    std::uint64_t critical_escape_failures = 0U;
    std::uint64_t nongeometric_critical_escape_failures = 0U;
    std::uint64_t critical_geometric_tails = 0U;
    std::uint64_t nongeometric_critical_tail_failures = 0U;
    int maximum_critical_escape_depth = 0;
    std::array<int, 3> maximum_critical_escape_state{0, 0, 0};
    std::array<int, 3> first_nongeometric_escape_state{0, 0, 0};
    std::array<int, 3> first_nongeometric_tail_failure{0, 0, 0};
    std::vector<std::string> critical_nonterminal_policies;
    if (!policies.empty()) {
      for (int r_index = 0; r_index <= maximum; ++r_index) {
        for (int s_index = 0; s_index <= r_index; ++s_index) {
          for (int t_index = 0; t_index <= s_index; ++t_index) {
            const int q_index = final_policy[
                flat_index(r_index, s_index, t_index)];
            const long double r = static_cast<long double>(r_index)
                                  / static_cast<long double>(denominator);
            const long double s = static_cast<long double>(s_index)
                                  / static_cast<long double>(denominator);
            const long double t = static_cast<long double>(t_index)
                                  / static_cast<long double>(denominator);
            long double ignored_demand = 0.0L;
            const bool critical = critical_demand(r, s, t, ignored_demand);
            if (q_index < 0) {
              ++policy_terminal;
              if (critical) {
                ++critical_policy_terminal;
              }
            } else if (q_index == 0) {
              ++policy_zero;
              if (critical) {
                ++critical_policy_zero;
              }
            } else if (q_index == t_index) {
              ++policy_upper;
              if (critical) {
                ++critical_policy_upper;
              }
            } else {
              ++policy_interior;
              if (critical) {
                ++critical_policy_interior;
              }
            }
            if (critical && q_index >= 0
                && critical_nonterminal_policies.size() < 12U) {
              critical_nonterminal_policies.push_back(
                  std::to_string(r_index) + '/' + std::to_string(denominator)
                  + ',' + std::to_string(s_index) + '/'
                  + std::to_string(denominator) + ','
                  + std::to_string(t_index) + '/'
                  + std::to_string(denominator) + "->"
                  + std::to_string(q_index) + '/'
                  + std::to_string(denominator));
            }
            if (critical) {
              int trace_r = r_index;
              int trace_s = s_index;
              int trace_t = t_index;
              int escape_depth = 0;
              bool escaped = false;
              bool reached_geometric_tail = false;
              for (int step = 0; step < depth; ++step) {
                if (trace_r == trace_s && trace_s == trace_t) {
                  reached_geometric_tail = true;
                  break;
                }
                long double trace_demand = 0.0L;
                if (!critical_demand(
                        static_cast<long double>(trace_r)
                            / static_cast<long double>(denominator),
                        static_cast<long double>(trace_s)
                            / static_cast<long double>(denominator),
                        static_cast<long double>(trace_t)
                            / static_cast<long double>(denominator),
                        trace_demand
                    )) {
                  escaped = true;
                  break;
                }
                const int trace_q = final_policy[
                    flat_index(trace_r, trace_s, trace_t)];
                if (trace_q < 0) {
                  escaped = true;
                  break;
                }
                ++escape_depth;
                trace_r = trace_s;
                trace_s = trace_t;
                trace_t = trace_q;
              }
              if (!escaped) {
                ++critical_escape_failures;
                if (r_index != s_index || s_index != t_index) {
                  if (nongeometric_critical_escape_failures == 0U) {
                    first_nongeometric_escape_state = {
                        r_index, s_index, t_index};
                  }
                  ++nongeometric_critical_escape_failures;
                }
              }
              if (reached_geometric_tail) {
                ++critical_geometric_tails;
              }
              if (!escaped && !reached_geometric_tail) {
                if (trace_r != trace_s || trace_s != trace_t) {
                  if (nongeometric_critical_tail_failures == 0U) {
                    first_nongeometric_tail_failure = {
                        trace_r, trace_s, trace_t};
                  }
                  ++nongeometric_critical_tail_failures;
                }
              }
              if (escape_depth > maximum_critical_escape_depth) {
                maximum_critical_escape_depth = escape_depth;
                maximum_critical_escape_state = {
                    r_index, s_index, t_index};
              }
            }
          }
        }
      }
    }
    const char* candidate_name = "1+t";
    if (candidate_one) {
      candidate_name = "1";
    } else if (candidate_clipped) {
      candidate_name = "min(T4,1+r)";
    } else if (candidate_saturated) {
      candidate_name = "saturated-2";
    } else if (candidate_saturated_deep) {
      candidate_name = "saturated-8";
    } else if (candidate_tangent) {
      candidate_name = "min(T4,B_tan)";
    } else if (candidate_shift) {
      candidate_name = "min(T4,B_shift)";
    } else if (candidate_downward) {
      candidate_name = "min(T4,1+2t-r)";
    } else if (candidate_last) {
      candidate_name = "min(T4,1+t)";
    } else if (candidate_steep) {
      candidate_name = "min(T4,B_steep)";
    } else if (candidate_quadratic) {
      candidate_name = "min(T4,B_quadratic)";
    } else if (candidate_r) {
      candidate_name = "1+r";
    }
    std::cout << "SU2_WALL_121_BELLMAN_BARRIER_PROBE"
              << " candidate=" << candidate_name
              << " grid_states=" << states
              << " bellman_depth=" << depth;
    print_minimum("terminal_residual", terminal, denominator);
    print_minimum("terminal_floor_residual", terminal_floor, denominator);
    print_minimum("renewal_residual", renewal_barrier, denominator);
    print_minimum("critical_renewal_residual", critical_renewal_barrier,
                  denominator);
    print_minimum("append_delta_q_ge_1", append_above_one, denominator);
    print_minimum("payment_terminal_residual", terminal_payment, denominator);
    print_minimum("payment_barrier_residual", barrier_payment, denominator);
    print_minimum("payment_renewal_residual", renewal_payment, denominator);
    print_minimum("critical_payment_renewal_residual",
                  critical_renewal_payment, denominator);
    print_minimum("critical_to_noncritical_payment_residual",
                  critical_to_noncritical_payment, denominator);
    print_minimum("critical_to_nongeometric_payment_residual",
                  critical_to_nongeometric_payment, denominator);
    print_minimum("critical_unit_floor_renewal_residual",
                  critical_unit_floor_payment, denominator);
    print_minimum("critical_unit_floor_q_le_one_residual",
                  critical_unit_floor_q_le_one_payment, denominator);
    print_minimum("critical_unit_floor_q_lt_one_residual",
                  critical_unit_floor_q_lt_one_payment, denominator);
    if (renewal_barrier.q >= 0) {
      const long double r = static_cast<long double>(renewal_barrier.r)
                            / static_cast<long double>(denominator);
      const long double s = static_cast<long double>(renewal_barrier.s)
                            / static_cast<long double>(denominator);
      const long double t = static_cast<long double>(renewal_barrier.t)
                            / static_cast<long double>(denominator);
      const long double q = static_cast<long double>(renewal_barrier.q)
                            / static_cast<long double>(denominator);
      long double current_demand = 0.0L;
      const bool current_critical = critical_demand(r, s, t, current_demand);
      long double next_demand = 0.0L;
      const bool next_critical = critical_demand(s, t, q, next_demand);
      std::cout << " renewal_components=(T4=" << terminal_current(r, s, t)
                << ",B=" << barrier(r, s, t)
                << ",R=" << renewal(r, s, t, q)
                << ",next_T4=" << terminal_current(s, t, q)
                << ",next_B=" << barrier(s, t, q)
                << ",critical=" << current_critical
                << ",demand=" << current_demand
                << ",next_critical=" << next_critical
                << ",next_demand=" << next_demand << ')';
    }
    print_minimum("critical_bellman_margin", bellman_margin, denominator);
    std::cout << " final_policy=(terminal=" << policy_terminal
              << ",q0=" << policy_zero
              << ",qt=" << policy_upper
              << ",interior=" << policy_interior << ')';
    std::cout << " critical_final_policy=(terminal="
              << critical_policy_terminal << ",q0="
              << critical_policy_zero << ",qt=" << critical_policy_upper
              << ",interior=" << critical_policy_interior << ')';
    std::cout << " finite_horizon_critical_nonexit=(max_steps="
              << maximum_critical_escape_depth << "@"
              << maximum_critical_escape_state[0] << '/' << denominator
              << ',' << maximum_critical_escape_state[1] << '/'
              << denominator << ',' << maximum_critical_escape_state[2]
              << '/' << denominator << ",unescaped="
              << critical_escape_failures << ",nongeometric_unescaped="
              << nongeometric_critical_escape_failures;
    if (nongeometric_critical_escape_failures != 0U) {
      std::cout << "@" << first_nongeometric_escape_state[0] << '/'
                << denominator << ',' << first_nongeometric_escape_state[1]
                << '/' << denominator << ','
                << first_nongeometric_escape_state[2] << '/' << denominator;
    }
    std::cout << ')';
    std::cout << " finite_horizon_critical_resolution=(geometric_tails="
              << critical_geometric_tails << ",unresolved="
              << nongeometric_critical_tail_failures;
    if (nongeometric_critical_tail_failures != 0U) {
      std::cout << "@" << first_nongeometric_tail_failure[0] << '/'
                << denominator << ',' << first_nongeometric_tail_failure[1]
                << '/' << denominator << ','
                << first_nongeometric_tail_failure[2] << '/'
                << denominator;
    }
    std::cout << ')';
    std::cout << " critical_nonterminal_sample={";
    for (std::size_t index = 0U;
         index < critical_nonterminal_policies.size(); ++index) {
      if (index != 0U) {
        std::cout << ';';
      }
      std::cout << critical_nonterminal_policies[index];
    }
    std::cout << '}';
    std::cout << " minimum_margin_final_horizon_policy={";
    int trace_r = bellman_margin.r;
    int trace_s = bellman_margin.s;
    int trace_t = bellman_margin.t;
    bool first_control = true;
    for (int iteration = depth - 1; iteration >= 0; --iteration) {
      const int q_index = final_policy[
          flat_index(trace_r, trace_s, trace_t)];
      if (q_index < 0) {
        break;
      }
      if (!first_control) {
        std::cout << ',';
      }
      first_control = false;
      std::cout << q_index << '/' << denominator;
      trace_r = trace_s;
      trace_s = trace_t;
      trace_t = q_index;
    }
    std::cout << '}';
    if (nongeometric_critical_escape_failures != 0U) {
      std::cout << " first_nongeometric_final_horizon_policy={";
      int escape_r = first_nongeometric_escape_state[0];
      int escape_s = first_nongeometric_escape_state[1];
      int escape_t = first_nongeometric_escape_state[2];
      bool first_state = true;
      for (int iteration = 0; iteration < depth; ++iteration) {
        const int q_index = final_policy[
            flat_index(escape_r, escape_s, escape_t)];
        if (!first_state) {
          std::cout << ';';
        }
        first_state = false;
        std::cout << escape_r << '/' << denominator << ',' << escape_s
                  << '/' << denominator << ',' << escape_t << '/'
                  << denominator;
        if (q_index < 0) {
          std::cout << "->terminal";
          break;
        }
        std::cout << "->" << q_index << '/' << denominator;
        escape_r = escape_s;
        escape_s = escape_t;
        escape_t = q_index;
      }
      std::cout << '}';
    }
    std::cout << " longest_final_horizon_critical_policy={";
    int longest_r = maximum_critical_escape_state[0];
    int longest_s = maximum_critical_escape_state[1];
    int longest_t = maximum_critical_escape_state[2];
    bool first_longest_control = true;
    for (int iteration = 0; iteration < depth; ++iteration) {
      if (longest_r == longest_s && longest_s == longest_t) {
        std::cout << "->geometric";
        break;
      }
      const int q_index = final_policy[
          flat_index(longest_r, longest_s, longest_t)];
      if (q_index < 0) {
        std::cout << "->terminal";
        break;
      }
      if (!first_longest_control) {
        std::cout << ',';
      }
      first_longest_control = false;
      std::cout << q_index << '/' << denominator;
      longest_r = longest_s;
      longest_s = longest_t;
      longest_t = q_index;
    }
    std::cout << '}';
    const bool candidate_passes_grid =
        terminal.value >= -1.0e-15L
        && renewal_barrier.value >= -1.0e-15L;
    const bool bellman_margin_passes_grid =
        bellman_margin.value >= -1.0e-15L;
    std::cout << " candidate_grid_result="
              << (candidate_passes_grid ? "NO_COUNTEREXAMPLE"
                                        : "COUNTEREXAMPLE")
              << " bellman_grid_result="
              << (bellman_margin_passes_grid
                      ? "NO_NEGATIVE_CRITICAL_MARGIN"
                      : "NEGATIVE_CRITICAL_MARGIN")
              << '\n';
    return candidate_passes_grid ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
