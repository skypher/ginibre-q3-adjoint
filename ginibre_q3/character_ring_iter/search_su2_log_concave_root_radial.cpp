#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

namespace {

int parse_positive(const char* text, const std::string& name) {
  const std::string value{text};
  std::size_t consumed = 0U;
  const long parsed = std::stol(value, &consumed, 10);
  if (consumed != value.size() || parsed <= 0) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  return static_cast<int>(parsed);
}

std::uint64_t splitmix64(std::uint64_t& state) {
  state += 0x9e3779b97f4a7c15ULL;
  std::uint64_t value = state;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

long double uniform01(std::uint64_t& state) {
  const std::uint64_t value = splitmix64(state) >> 11U;
  return static_cast<long double>(value) /
         static_cast<long double>(std::uint64_t{1} << 53U);
}

long double profile_value(const std::vector<long double>& profile,
                          const int index) {
  return index >= 0 && index < static_cast<int>(profile.size())
             ? profile[static_cast<std::size_t>(index)]
             : 0.0L;
}

struct Best {
  long double margin = std::numeric_limits<long double>::infinity();
  int support = -1;
  int antidiagonal = -1;
  int depth = -1;
  std::vector<long double> root;
};

Best inspect(const std::vector<long double>& root) {
  const int support = static_cast<int>(root.size()) - 1;
  std::vector<long double> half(root.size());
  long double tail = 0.0L;
  for (int index = support; index >= 0; --index) {
    tail += root[static_cast<std::size_t>(index)];
    half[static_cast<std::size_t>(index)] = tail;
  }

  std::vector<long double> autocorrelation(
      static_cast<std::size_t>(2 * support + 2), 0.0L);
  for (int shift = 0; shift <= 2 * support; ++shift) {
    long double value = 0.0L;
    for (int index = -support; index + shift <= support; ++index) {
      value +=
          half[static_cast<std::size_t>(std::abs(index))] *
          half[static_cast<std::size_t>(std::abs(index + shift))];
    }
    autocorrelation[static_cast<std::size_t>(shift)] = value;
  }
  std::vector<long double> square(
      static_cast<std::size_t>(2 * support + 1));
  for (int index = 0; index <= 2 * support; ++index) {
    square[static_cast<std::size_t>(index)] =
        autocorrelation[static_cast<std::size_t>(index)] -
        autocorrelation[static_cast<std::size_t>(index + 1)];
  }
  const long double scale = square.front() * square.front();
  Best best;
  best.support = support;
  for (int antidiagonal = 1; antidiagonal <= 2 * support - 2;
       ++antidiagonal) {
    for (int depth = 0; 2 * depth < antidiagonal; ++depth) {
      const int complement = antidiagonal - depth;
      const long double radial =
          square.front() *
              (profile_value(square, antidiagonal + 1) +
               profile_value(square, antidiagonal + 2)) +
          profile_value(square, depth) *
              profile_value(square, complement) -
          profile_value(square, depth + 1) *
              profile_value(square, complement + 1);
      const long double normalized = radial / scale;
      if (normalized < best.margin) {
        best.margin = normalized;
        best.antidiagonal = antidiagonal;
        best.depth = depth;
      }
    }
  }
  return best;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const int samples =
        argc >= 2 ? parse_positive(argv[1], "samples") : 200000;
    const int maximum_support =
        argc >= 3 ? parse_positive(argv[2], "maximum_support") : 50;
    if (argc > 3 || maximum_support < 2) {
      throw std::invalid_argument(
          "usage: search_su2_log_concave_root_radial "
          "[samples>=1] [maximum_support>=2]");
    }

    Best global;
    std::atomic<std::uint64_t> checked{0U};
#pragma omp parallel
    {
      const int thread = omp_get_thread_num();
      std::uint64_t state =
          std::uint64_t{0x6a09e667f3bcc909ULL} ^
          (static_cast<std::uint64_t>(thread) *
           std::uint64_t{0x9e3779b97f4a7c15ULL});
      Best local;
#pragma omp for schedule(dynamic, 64)
      for (int sample = 0; sample < samples; ++sample) {
        const int support =
            2 + static_cast<int>(
                    uniform01(state) *
                    static_cast<long double>(maximum_support - 1));
        const long double log_range =
            std::exp(
                std::log(0.05L) +
                uniform01(state) *
                    (std::log(30.0L) - std::log(0.05L)));
        std::vector<long double> slopes(
            static_cast<std::size_t>(support));
        for (long double& slope : slopes) {
          slope =
              log_range * (2.0L * uniform01(state) - 1.0L);
        }
        std::sort(slopes.begin(), slopes.end(), std::greater<>());
        std::vector<long double> logarithms(
            static_cast<std::size_t>(support + 1), 0.0L);
        for (int index = 1; index <= support; ++index) {
          logarithms[static_cast<std::size_t>(index)] =
              logarithms[static_cast<std::size_t>(index - 1)] +
              slopes[static_cast<std::size_t>(index - 1)];
        }
        const long double maximum =
            *std::max_element(logarithms.begin(), logarithms.end());
        std::vector<long double> root(logarithms.size());
        for (std::size_t index = 0; index < logarithms.size(); ++index) {
          root[index] = std::exp(logarithms[index] - maximum);
        }
        Best candidate = inspect(root);
        if (candidate.margin < local.margin) {
          local = std::move(candidate);
          local.root = std::move(root);
        }
        ++checked;
      }
#pragma omp critical
      {
        if (local.margin < global.margin) {
          global = std::move(local);
        }
      }
    }

    std::cout
        << std::setprecision(20)
        << "SU2_LOG_CONCAVE_ROOT_RADIAL_SEARCH"
        << " samples=" << checked.load()
        << " maximum_support=" << maximum_support
        << " threads=" << omp_get_max_threads()
        << " minimum_normalized_margin=" << global.margin
        << " support=" << global.support
        << " A=" << global.antidiagonal
        << " L=" << global.depth
        << " result="
        << (global.margin < -1.0e-12L ? "NEGATIVE_FOUND" : "NO_NEGATIVE")
        << '\n';
    if (global.margin < -1.0e-12L) {
      std::cout << "root=[";
      for (std::size_t index = 0; index < global.root.size(); ++index) {
        if (index != 0U) {
          std::cout << ',';
        }
        std::cout << global.root[index];
      }
      std::cout << "]\n";
    }
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
