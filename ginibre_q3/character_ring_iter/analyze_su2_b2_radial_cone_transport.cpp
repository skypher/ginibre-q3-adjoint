#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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

bool triangle(const int first, const int second, const int output) {
  return std::abs(first - second) <= output &&
         output <= first + second;
}

int square_factor_entry(const int q, const int row, const int column) {
  int result = 0;
  for (int label = 0; label <= 2 * q; ++label) {
    if (triangle(label, column, row)) {
      ++result;
    }
  }
  return result;
}

long long compound_entry(const int q, const int out_first,
                         const int out_second, const int in_first,
                         const int in_second) {
  return static_cast<long long>(
             square_factor_entry(q, out_first, in_first)) *
             square_factor_entry(q, out_second, in_second) -
         static_cast<long long>(
             square_factor_entry(q, out_first, in_second)) *
             square_factor_entry(q, out_second, in_first);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const int maximum_q =
        argc >= 2 ? parse_positive(argv[1], "maximum_q") : 8;
    const int maximum_a =
        argc >= 3 ? parse_positive(argv[2], "maximum_a") : 32;
    const bool dump_first_output =
        argc >= 4 && std::string(argv[3]) == "--dump-first-output";
    if (argc >= 4 && !dump_first_output) {
      throw std::invalid_argument(
          "third argument must be --dump-first-output");
    }
    if (argc > 4) {
      throw std::invalid_argument(
          "usage: analyze_su2_b2_radial_cone_transport "
          "[maximum_q] [maximum_a] [--dump-first-output]");
    }

    std::size_t output_functionals = 0U;
    std::size_t input_antidiagonals = 0U;
    std::size_t dual_coefficients = 0U;
    std::size_t failures = 0U;
    int first_q = -1;
    int first_a = -1;
    int first_l = -1;
    int first_b = -1;
    int first_r = -1;
    long long first_value = 0;
    std::string first_kind;

    for (int q = 1; q <= maximum_q; ++q) {
      for (int a = 1; a <= maximum_a; ++a) {
        for (int l = 0; 2 * l < a; ++l) {
          ++output_functionals;
          const int maximum_input_index = a + 4 * q + 2;
          const int maximum_b = 2 * maximum_input_index - 1;
          std::vector<long long> central_coefficient(
              static_cast<std::size_t>(maximum_input_index), 0);
          for (int b = 0; b <= maximum_b; ++b) {
            const int maximum_r = b / 2;
            std::vector<long long> coefficient(
                static_cast<std::size_t>(maximum_r + 1), 0);
            bool nonzero = false;
            for (int r = 0; r <= maximum_r; ++r) {
              const int s = b - r + 1;
              if (s <= r || s > maximum_input_index) {
                continue;
              }
              long long value = 0;
              for (int i = 0; i <= l; ++i) {
                value += compound_entry(q, i, a - i + 1, r, s);
              }
              coefficient[static_cast<std::size_t>(r)] = value;
              nonzero = nonzero || value != 0;
            }
            if (!nonzero) {
              continue;
            }
            ++input_antidiagonals;
            if (dump_first_output && q == 1 && a == 1 && l == 0) {
              std::cout << "raw_antidiagonal b=" << b
                        << " coefficients=[";
              for (int r = 0; r <= maximum_r; ++r) {
                if (r != 0) {
                  std::cout << ',';
                }
                std::cout
                    << coefficient[static_cast<std::size_t>(r)];
              }
              std::cout << "]\n";
            }

            const bool even = b % 2 == 0;
            const int last_constrained =
                even ? maximum_r - 1 : maximum_r;
            if (even) {
              central_coefficient[static_cast<std::size_t>(maximum_r)] =
                  coefficient[static_cast<std::size_t>(maximum_r)];
            }
            for (int r = 0; r <= last_constrained; ++r) {
              const long long next =
                  r == maximum_r
                      ? 0
                      : coefficient[static_cast<std::size_t>(r + 1)];
              const long long prefix_coefficient =
                  coefficient[static_cast<std::size_t>(r)] - next;
              if (dump_first_output && q == 1 && a == 1 && l == 0 &&
                  prefix_coefficient != 0) {
                std::cout << "strict_prefix_dual b=" << b
                          << " r=" << r
                          << " value=" << prefix_coefficient << '\n';
              }
              ++dual_coefficients;
              if (prefix_coefficient < 0) {
                ++failures;
                if (first_q < 0) {
                  first_q = q;
                  first_a = a;
                  first_l = l;
                  first_b = b;
                  first_r = r;
                  first_value = prefix_coefficient;
                  first_kind = "negative_prefix_dual";
                }
              }
            }
          }
          for (int r = 0; r < maximum_input_index; ++r) {
            const long long next =
                r + 1 < maximum_input_index
                    ? central_coefficient[static_cast<std::size_t>(r + 1)]
                    : 0;
            const long long diagonal_dual =
                central_coefficient[static_cast<std::size_t>(r)] - next;
            if (dump_first_output && q == 1 && a == 1 && l == 0 &&
                diagonal_dual != 0) {
              std::cout << "diagonal_dual r=" << r
                        << " value=" << diagonal_dual << '\n';
            }
            ++dual_coefficients;
            if (diagonal_dual < 0) {
              ++failures;
              if (first_q < 0) {
                first_q = q;
                first_a = a;
                first_l = l;
                first_b = 2 * r;
                first_r = r;
                first_value = diagonal_dual;
                first_kind = "negative_diagonal_dual";
              }
            }
          }
        }
      }
    }

    std::cout << "SU2_B2_RADIAL_CONE_TRANSPORT"
              << " maximum_q=" << maximum_q
              << " maximum_a=" << maximum_a
              << " output_functionals=" << output_functionals
              << " input_antidiagonals=" << input_antidiagonals
              << " dual_coefficients=" << dual_coefficients
              << " failures=" << failures;
    if (first_q >= 0) {
      std::cout << " first_kind=" << first_kind
                << " first_q=" << first_q
                << " first_a=" << first_a
                << " first_l=" << first_l
                << " first_b=" << first_b
                << " first_r=" << first_r
                << " first_value=" << first_value;
    }
    std::cout << " result=" << (failures == 0U ? "PASS" : "FAIL")
              << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
