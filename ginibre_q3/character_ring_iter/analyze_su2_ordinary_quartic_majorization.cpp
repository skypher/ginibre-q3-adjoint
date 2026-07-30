#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Monomial = std::array<int, 4>;
using CoefficientMap = std::map<Monomial, std::int64_t>;
using Integer = boost::multiprecision::cpp_int;
using Exponent = std::vector<int>;
using Polynomial = std::map<Exponent, Integer>;

bool triangle(const int first, const int second, const int label) {
  return std::abs(first - second) <= label && label <= first + second;
}

Monomial monomial(const int first, const int second, const int third,
                    const int fourth) {
  Monomial result{first, second, third, fourth};
  std::sort(result.begin(), result.end(), std::greater<int>());
  return result;
}

CoefficientMap current_polynomial(const int left, const int right,
                                  const int radius, const int target) {
  CoefficientMap coefficients;
  const int output_right = right + std::max(radius, target);

  for (int repeated = left; repeated <= right; ++repeated) {
    for (int output = 0; output <= output_right; ++output) {
      for (int first = left; first <= right; ++first) {
        if (!triangle(first, output, radius)) {
          continue;
        }
        for (int second = left; second <= right; ++second) {
          if (triangle(second, output, target)) {
            ++coefficients[monomial(repeated, repeated, first, second)];
          }
        }
      }
    }
  }

  for (int first = left; first <= right; ++first) {
    for (int second = left; second <= right; ++second) {
      if (!triangle(first, second, radius)) {
        continue;
      }
      for (int third = left; third <= right; ++third) {
        for (int fourth = left; fourth <= right; ++fourth) {
          if (triangle(third, fourth, target)) {
            --coefficients[monomial(first, second, third, fourth)];
          }
        }
      }
    }
  }
  return coefficients;
}

bool majorizes(const Monomial& spread, const Monomial& balanced) {
  int spread_prefix = 0;
  int balanced_prefix = 0;
  for (std::size_t index = 0; index + 1U < spread.size(); ++index) {
    spread_prefix += spread[index];
    balanced_prefix += balanced[index];
    if (spread_prefix < balanced_prefix) {
      return false;
    }
  }
  return true;
}

struct Edge {
  int destination;
  int reverse;
  std::int64_t capacity;
};

class Dinic {
 public:
  explicit Dinic(const int vertices)
      : graph_(static_cast<std::size_t>(vertices)),
        level_(static_cast<std::size_t>(vertices)),
        next_(static_cast<std::size_t>(vertices)) {}

  void add_edge(const int source, const int destination,
                const std::int64_t capacity) {
    const int source_reverse =
        static_cast<int>(graph_[static_cast<std::size_t>(destination)].size());
    const int destination_reverse =
        static_cast<int>(graph_[static_cast<std::size_t>(source)].size());
    graph_[static_cast<std::size_t>(source)].push_back(
        Edge{destination, source_reverse, capacity});
    graph_[static_cast<std::size_t>(destination)].push_back(
        Edge{source, destination_reverse, 0});
  }

  std::int64_t maximum_flow(const int source, const int sink) {
    std::int64_t result = 0;
    while (build_levels(source, sink)) {
      std::fill(next_.begin(), next_.end(), 0);
      while (true) {
        const std::int64_t sent =
            send(source, sink, std::numeric_limits<std::int64_t>::max());
        if (sent == 0) {
          break;
        }
        result += sent;
      }
    }
    return result;
  }

 private:
  bool build_levels(const int source, const int sink) {
    std::fill(level_.begin(), level_.end(), -1);
    std::queue<int> queue;
    level_[static_cast<std::size_t>(source)] = 0;
    queue.push(source);
    while (!queue.empty()) {
      const int vertex = queue.front();
      queue.pop();
      for (const Edge& edge : graph_[static_cast<std::size_t>(vertex)]) {
        if (edge.capacity > 0 &&
            level_[static_cast<std::size_t>(edge.destination)] < 0) {
          level_[static_cast<std::size_t>(edge.destination)] =
              level_[static_cast<std::size_t>(vertex)] + 1;
          queue.push(edge.destination);
        }
      }
    }
    return level_[static_cast<std::size_t>(sink)] >= 0;
  }

  std::int64_t send(const int vertex, const int sink,
                    const std::int64_t available) {
    if (vertex == sink) {
      return available;
    }
    int& next = next_[static_cast<std::size_t>(vertex)];
    auto& edges = graph_[static_cast<std::size_t>(vertex)];
    while (next < static_cast<int>(edges.size())) {
      Edge& edge = edges[static_cast<std::size_t>(next)];
      if (edge.capacity > 0 &&
          level_[static_cast<std::size_t>(edge.destination)] ==
              level_[static_cast<std::size_t>(vertex)] + 1) {
        const std::int64_t sent =
            send(edge.destination, sink, std::min(available, edge.capacity));
        if (sent > 0) {
          edge.capacity -= sent;
          Edge& reverse =
              graph_[static_cast<std::size_t>(edge.destination)]
                    [static_cast<std::size_t>(edge.reverse)];
          reverse.capacity += sent;
          return sent;
        }
      }
      ++next;
    }
    return 0;
  }

  std::vector<std::vector<Edge>> graph_;
  std::vector<int> level_;
  std::vector<int> next_;
};

struct TransportResult {
  std::int64_t demand = 0;
  std::int64_t flow = 0;
  int totals = 0;
};

TransportResult check_transport(const CoefficientMap& coefficients) {
  std::map<int, std::vector<std::pair<Monomial, std::int64_t>>> negative;
  std::map<int, std::vector<std::pair<Monomial, std::int64_t>>> positive;
  for (const auto& [term, coefficient] : coefficients) {
    const int total = term[0] + term[1] + term[2] + term[3];
    if (coefficient < 0) {
      negative[total].push_back({term, -coefficient});
    } else if (coefficient > 0) {
      positive[total].push_back({term, coefficient});
    }
  }

  TransportResult result;
  for (const auto& [total, demands] : negative) {
    ++result.totals;
    const auto supply_iterator = positive.find(total);
    if (supply_iterator == positive.end()) {
      for (const auto& [term, capacity] : demands) {
        static_cast<void>(term);
        result.demand += capacity;
      }
      continue;
    }
    const auto& supplies = supply_iterator->second;
    const int source = 0;
    const int demand_begin = 1;
    const int supply_begin =
        demand_begin + static_cast<int>(demands.size());
    const int sink = supply_begin + static_cast<int>(supplies.size());
    Dinic flow(sink + 1);
    std::int64_t total_demand = 0;
    for (std::size_t index = 0; index < demands.size(); ++index) {
      const std::int64_t capacity = demands[index].second;
      total_demand += capacity;
      flow.add_edge(source, demand_begin + static_cast<int>(index), capacity);
    }
    for (std::size_t index = 0; index < supplies.size(); ++index) {
      flow.add_edge(supply_begin + static_cast<int>(index), sink,
                    supplies[index].second);
    }
    for (std::size_t demand = 0; demand < demands.size(); ++demand) {
      for (std::size_t supply = 0; supply < supplies.size(); ++supply) {
        if (majorizes(demands[demand].first, supplies[supply].first)) {
          flow.add_edge(demand_begin + static_cast<int>(demand),
                        supply_begin + static_cast<int>(supply),
                        total_demand);
        }
      }
    }
    result.demand += total_demand;
    result.flow += flow.maximum_flow(source, sink);
  }
  return result;
}

int parse_nonnegative(const char* text, const std::string& name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value, &consumed, 10);
  if (consumed != value.size() || parsed < 0) {
    throw std::invalid_argument(name + " must be a nonnegative integer");
  }
  return static_cast<int>(parsed);
}

std::string render(const Monomial& term) {
  return "[" + std::to_string(term[0]) + "," + std::to_string(term[1]) +
         "," + std::to_string(term[2]) + "," + std::to_string(term[3]) +
         "]";
}

Polynomial constant_polynomial(const int variables, const Integer& value) {
  Polynomial result;
  if (value != 0) {
    result.emplace(Exponent(static_cast<std::size_t>(variables), 0), value);
  }
  return result;
}

Polynomial ratio_polynomial(const int variables, const int ratio) {
  Polynomial result;
  for (int drop = ratio; drop < variables; ++drop) {
    Exponent exponent(static_cast<std::size_t>(variables), 0);
    exponent[static_cast<std::size_t>(drop)] = 1;
    result.emplace(std::move(exponent), 1);
  }
  return result;
}

Polynomial multiply(const Polynomial& left, const Polynomial& right) {
  Polynomial result;
  for (const auto& [left_exponent, left_coefficient] : left) {
    for (const auto& [right_exponent, right_coefficient] : right) {
      Exponent exponent = left_exponent;
      for (std::size_t index = 0; index < exponent.size(); ++index) {
        exponent[index] += right_exponent[index];
      }
      result[exponent] += left_coefficient * right_coefficient;
    }
  }
  return result;
}

Polynomial ratio_drop_current(const CoefficientMap& coefficients,
                              const int left, const int right) {
  const int variables = right - left;
  std::vector<Polynomial> profile;
  profile.reserve(static_cast<std::size_t>(variables + 1));
  Polynomial entry = constant_polynomial(variables, 1);
  profile.push_back(entry);
  for (int ratio = 0; ratio < variables; ++ratio) {
    entry = multiply(entry, ratio_polynomial(variables, ratio));
    profile.push_back(entry);
  }

  Polynomial result;
  for (const auto& [term, coefficient] : coefficients) {
    if (coefficient == 0) {
      continue;
    }
    Polynomial expanded = constant_polynomial(variables, coefficient);
    for (const int index : term) {
      expanded =
          multiply(expanded, profile[static_cast<std::size_t>(index - left)]);
    }
    for (const auto& [exponent, expanded_coefficient] : expanded) {
      result[exponent] += expanded_coefficient;
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

std::string render(const Exponent& exponent) {
  std::string result = "[";
  for (std::size_t index = 0; index < exponent.size(); ++index) {
    if (index != 0U) {
      result += ",";
    }
    result += std::to_string(exponent[index]);
  }
  return result + "]";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const int maximum_label =
        argc >= 2 ? parse_nonnegative(argv[1], "maximum_label") : 5;
    const int maximum_left =
        argc >= 3 ? parse_nonnegative(argv[2], "maximum_left") : 3;
    const int maximum_width =
        argc >= 4 ? parse_nonnegative(argv[3], "maximum_width") : 7;
    if (argc > 4 || maximum_label < 1 || maximum_width < 1) {
      throw std::invalid_argument(
          "usage: analyze_su2_ordinary_quartic_majorization "
          "[maximum_label>=1] [maximum_left>=0] [maximum_width>=1]");
    }

    std::uint64_t cases = 0;
    std::uint64_t failures = 0;
    std::int64_t total_demand = 0;
    std::int64_t total_flow = 0;
    std::uint64_t ratio_drop_coefficients = 0;
    std::uint64_t negative_ratio_drop_coefficients = 0;
    std::uint64_t negative_curvature_coefficients = 0;
    int first_left = -1;
    int first_right = -1;
    int first_radius = -1;
    int first_target = -1;
    TransportResult first_result;
    CoefficientMap first_coefficients;
    int first_ratio_left = -1;
    int first_ratio_right = -1;
    int first_ratio_radius = -1;
    int first_ratio_target = -1;
    Exponent first_ratio_exponent;
    Integer first_ratio_coefficient = 0;
    int first_curvature_left = -1;
    int first_curvature_right = -1;
    int first_curvature_radius = -1;
    int first_curvature_target = -1;
    Exponent first_curvature_exponent;
    Integer first_curvature_coefficient = 0;

    for (int left = 0; left <= maximum_left; ++left) {
      for (int width = 1; width <= maximum_width; ++width) {
        const int right = left + width - 1;
        for (int radius = 1; radius <= maximum_label; ++radius) {
          for (int target = radius + 1; target <= maximum_label; ++target) {
            const CoefficientMap coefficients =
                current_polynomial(left, right, radius, target);
            const TransportResult result = check_transport(coefficients);
            const Polynomial ratio_drop =
                ratio_drop_current(coefficients, left, right);
            ++cases;
            total_demand += result.demand;
            total_flow += result.flow;
            if (result.flow != result.demand) {
              ++failures;
              if (first_left < 0) {
                first_left = left;
                first_right = right;
                first_radius = radius;
                first_target = target;
                first_result = result;
                first_coefficients = coefficients;
              }
            }
            ratio_drop_coefficients += ratio_drop.size();
            for (const auto& [exponent, coefficient] : ratio_drop) {
              if (coefficient < 0) {
                ++negative_ratio_drop_coefficients;
                if (first_ratio_left < 0) {
                  first_ratio_left = left;
                  first_ratio_right = right;
                  first_ratio_radius = radius;
                  first_ratio_target = target;
                  first_ratio_exponent = exponent;
                  first_ratio_coefficient = coefficient;
                }
                const bool uses_curvature =
                    exponent.size() >= 2U &&
                    std::any_of(exponent.begin(), exponent.end() - 1,
                                [](const int power) { return power > 0; });
                if (uses_curvature) {
                  ++negative_curvature_coefficients;
                  if (first_curvature_left < 0) {
                    first_curvature_left = left;
                    first_curvature_right = right;
                    first_curvature_radius = radius;
                    first_curvature_target = target;
                    first_curvature_exponent = exponent;
                    first_curvature_coefficient = coefficient;
                  }
                }
              }
            }
          }
        }
      }
    }

    std::cout << "SU2_ORDINARY_QUARTIC_MAJORIZATION"
              << " maximum_label=" << maximum_label
              << " maximum_left=" << maximum_left
              << " maximum_width=" << maximum_width << " cases=" << cases
              << " demand=" << total_demand << " flow=" << total_flow
              << " failures=" << failures
              << " ratio_drop_coefficients=" << ratio_drop_coefficients
              << " negative_ratio_drop_coefficients="
              << negative_ratio_drop_coefficients
              << " negative_curvature_coefficients="
              << negative_curvature_coefficients << '\n';
    if (first_left >= 0) {
      std::cout << "first_left=" << first_left
                << " first_right=" << first_right
                << " first_radius=" << first_radius
                << " first_target=" << first_target
                << " first_demand=" << first_result.demand
                << " first_flow=" << first_result.flow << '\n';
      for (const auto& [term, coefficient] : first_coefficients) {
        if (coefficient != 0) {
          std::cout << "coefficient term=" << render(term)
                    << " value=" << coefficient << '\n';
        }
      }
    }
    if (first_ratio_left >= 0) {
      std::cout << "first_ratio_left=" << first_ratio_left
                << " first_ratio_right=" << first_ratio_right
                << " first_ratio_radius=" << first_ratio_radius
                << " first_ratio_target=" << first_ratio_target
                << " first_ratio_exponent=" << render(first_ratio_exponent)
                << " first_ratio_coefficient=" << first_ratio_coefficient
                << '\n';
    }
    if (first_curvature_left >= 0) {
      std::cout << "first_curvature_left=" << first_curvature_left
                << " first_curvature_right=" << first_curvature_right
                << " first_curvature_radius=" << first_curvature_radius
                << " first_curvature_target=" << first_curvature_target
                << " first_curvature_exponent="
                << render(first_curvature_exponent)
                << " first_curvature_coefficient="
                << first_curvature_coefficient << '\n';
    }
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
