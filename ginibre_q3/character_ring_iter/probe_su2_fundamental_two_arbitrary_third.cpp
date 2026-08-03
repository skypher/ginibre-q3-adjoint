#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

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
  const int parsed = parse_nonnegative(text, name);
  if (parsed == 0) {
    throw std::invalid_argument(name + " must be positive");
  }
  return parsed;
}

Integer at(const std::vector<Integer>& profile, const int index) {
  return index >= 0 && index < static_cast<int>(profile.size())
             ? profile[static_cast<std::size_t>(index)]
             : Integer(0);
}

std::vector<Integer> multiply_by_square(
    const std::vector<Integer>& profile, const int label) {
  const int old_support = static_cast<int>(profile.size()) - 1;
  std::vector<Integer> result(
      static_cast<std::size_t>(old_support + label + 1));
  for (int source = 0; source <= old_support; ++source) {
    for (int shell = 0; shell <= label; ++shell) {
      for (int target = std::abs(source - shell);
           target <= source + shell; ++target) {
        result[static_cast<std::size_t>(target)] +=
            profile[static_cast<std::size_t>(source)];
      }
    }
  }
  return result;
}

std::vector<Integer> profile_for(const int fundamentals, const int first,
                                 const int second) {
  std::vector<Integer> result{Integer(1)};
  for (int index = 0; index < fundamentals; ++index) {
    result = multiply_by_square(result, 1);
  }
  result = multiply_by_square(result, first);
  return multiply_by_square(result, second);
}

Integer radial(const std::vector<Integer>& profile, const int contraction) {
  const Integer zero = at(profile, 0);
  if (contraction == 0) {
    return zero * (at(profile, 3) + at(profile, 4) + at(profile, 5)) -
           at(profile, 1) * at(profile, 4);
  }
  return zero * (at(profile, 4) + at(profile, 5)) +
         at(profile, 1) * at(profile, 2) -
         at(profile, 2) * at(profile, 3);
}

bool zero_threshold(const std::vector<Integer>& profile) {
  return 7 * at(profile, 1) >= 19 * at(profile, 0);
}

bool one_threshold(const std::vector<Integer>& profile) {
  const Integer zero = at(profile, 0);
  const Integer delta = 3 * at(profile, 1) - 2 * zero;
  return delta >= 0 && delta * delta >= 40 * zero * zero;
}

void print_word(const int fundamentals, const int first, const int second) {
  std::cout << "[";
  for (int index = 0; index < fundamentals; ++index) {
    if (index != 0) {
      std::cout << ',';
    }
    std::cout << 1;
  }
  if (fundamentals != 0) {
    std::cout << ',';
  }
  std::cout << first << ',' << second << "]";
}

void inspect_case(const int fundamentals, const int first, const int second) {
  const std::vector<Integer> profile =
      profile_for(fundamentals, first, second);
  const Integer zero = at(profile, 0);
  const Integer one = at(profile, 1);
  const Integer zero_margin = 7 * one - 19 * zero;
  const Integer delta = 3 * one - 2 * zero;
  const Integer one_margin = delta * delta - 40 * zero * zero;
  std::cout << "SU2_FUNDAMENTAL_TWO_ARBITRARY_THIRD_CASE word=";
  print_word(fundamentals, first, second);
  std::cout << " c0=" << zero
            << " c1=" << one
            << " zero_threshold_margin=" << zero_margin
            << " one_threshold_margin=" << one_margin
            << " K30=" << radial(profile, 0)
            << " K31=" << radial(profile, 1)
            << " result=PASS_EXACT\n";
}

using Polynomial = std::vector<Integer>;

void trim(Polynomial& polynomial) {
  while (!polynomial.empty() && polynomial.back() == 0) {
    polynomial.pop_back();
  }
}

Polynomial linear(const int slope, const int intercept) {
  Polynomial result{Integer(intercept), Integer(slope)};
  trim(result);
  return result;
}

Polynomial multiply(const Polynomial& first, const Polynomial& second) {
  if (first.empty() || second.empty()) {
    return {};
  }
  Polynomial result(first.size() + second.size() - 1U);
  for (std::size_t left = 0U; left < first.size(); ++left) {
    for (std::size_t right = 0U; right < second.size(); ++right) {
      result[left + right] += first[left] * second[right];
    }
  }
  trim(result);
  return result;
}

void add_scaled(Polynomial& target, const Polynomial& source,
                const int scalar) {
  if (target.size() < source.size()) {
    target.resize(source.size());
  }
  for (std::size_t index = 0U; index < source.size(); ++index) {
    target[index] += Integer(scalar) * source[index];
  }
  trim(target);
}

Integer binomial(const int size, const int selection) {
  if (selection < 0 || selection > size) {
    return 0;
  }
  const int reduced = std::min(selection, size - selection);
  Integer result = 1;
  for (int index = 1; index <= reduced; ++index) {
    result *= size - reduced + index;
    result /= index;
  }
  return result;
}

Integer power(const int base, const int exponent) {
  Integer result = 1;
  for (int index = 0; index < exponent; ++index) {
    result *= base;
  }
  return result;
}

Polynomial shifted_coefficients(const Polynomial& polynomial,
                                const int shift) {
  Polynomial result(polynomial.size());
  for (std::size_t degree = 0U; degree < polynomial.size(); ++degree) {
    for (std::size_t lower = 0U; lower <= degree; ++lower) {
      const int exponent = static_cast<int>(degree - lower);
      result[lower] += polynomial[degree] * binomial(
          static_cast<int>(degree), static_cast<int>(lower)) *
          power(shift, exponent);
    }
  }
  trim(result);
  return result;
}

Integer evaluate(const Polynomial& polynomial, const int value) {
  Integer result = 0;
  for (auto iterator = polynomial.rbegin(); iterator != polynomial.rend();
       ++iterator) {
    result *= value;
    result += *iterator;
  }
  return result;
}

bool nonnegative_coefficients(const Polynomial& polynomial) {
  return std::all_of(
      polynomial.begin(), polynomial.end(),
      [](const Integer& coefficient) { return coefficient >= 0; });
}

void print_polynomial(const Polynomial& polynomial) {
  std::cout << '[';
  for (std::size_t index = 0U; index < polynomial.size(); ++index) {
    if (index != 0U) {
      std::cout << ',';
    }
    std::cout << polynomial[index];
  }
  std::cout << ']';
}

void verify_b2c3_threshold_polynomials() {
  constexpr int denominator_top = 7;
  constexpr int zero_coefficients[] = {3, 7, 8, 6, 3, 1};
  constexpr int one_coefficients[] = {7, 18, 21, 17, 10, 4, 1};
  Polynomial scaled_zero;
  Polynomial scaled_one;
  for (int shell = 0; shell <= 6; ++shell) {
    Polynomial ratio{Integer(2 * shell + 1)};
    for (int index = 0; index < shell; ++index) {
      ratio = multiply(ratio, linear(1, -index));
    }
    for (int index = shell + 2; index <= denominator_top; ++index) {
      ratio = multiply(ratio, linear(1, index));
    }
    if (shell < 6) {
      add_scaled(scaled_zero, ratio, zero_coefficients[shell]);
    }
    add_scaled(scaled_one, ratio, one_coefficients[shell]);
  }
  Polynomial zero_margin = scaled_one;
  for (Integer& coefficient : zero_margin) {
    coefficient *= 7;
  }
  add_scaled(zero_margin, scaled_zero, -19);
  Polynomial delta = scaled_one;
  for (Integer& coefficient : delta) {
    coefficient *= 3;
  }
  add_scaled(delta, scaled_zero, -2);
  Polynomial one_margin = multiply(delta, delta);
  const Polynomial zero_square = multiply(scaled_zero, scaled_zero);
  add_scaled(one_margin, zero_square, -40);

  Polynomial denominator{Integer(1)};
  for (int index = 2; index <= denominator_top; ++index) {
    denominator = multiply(denominator, linear(1, index));
  }
  for (const int fundamentals : {6, 10, 11, 16, 17}) {
    const std::vector<Integer> profile = profile_for(fundamentals, 2, 3);
    const Integer zero = at(profile, 0);
    const Integer direct_zero_margin =
        7 * at(profile, 1) - 19 * zero;
    const Integer direct_delta = 3 * at(profile, 1) - 2 * zero;
    const Integer direct_one_margin =
        direct_delta * direct_delta - 40 * zero * zero;
    const Integer divisor = evaluate(denominator, fundamentals);
    const Integer ballot_zero =
        binomial(2 * fundamentals, fundamentals) / (fundamentals + 1);
    if (evaluate(zero_margin, fundamentals) * ballot_zero
            != direct_zero_margin * divisor
        || evaluate(one_margin, fundamentals) * ballot_zero * ballot_zero
               != direct_one_margin * divisor * divisor) {
      throw std::runtime_error("threshold-polynomial reconstruction failed");
    }
  }
  const Polynomial zero_shift = shifted_coefficients(zero_margin, 11);
  const Polynomial one_shift = shifted_coefficients(one_margin, 17);
  if (!nonnegative_coefficients(zero_shift)
      || !nonnegative_coefficients(one_shift)
      || evaluate(zero_margin, 10) >= 0
      || evaluate(one_margin, 16) >= 0) {
    throw std::runtime_error("threshold-polynomial sign certificate failed");
  }
  std::cout << "SU2_FUNDAMENTAL_TWO_ARBITRARY_B23_THRESHOLDS"
            << " zero_shift=11 coefficients=";
  print_polynomial(zero_shift);
  std::cout << " one_shift=17 coefficients=";
  print_polynomial(one_shift);
  std::cout << " result=PASS_EXACT\n";
}

Polynomial scaled_b23_shell(const int output, const int denominator_top) {
  constexpr int factor_coefficients[] = {3, 7, 8, 6, 3, 1};
  Polynomial result;
  for (int shell = 0; shell <= output + 5; ++shell) {
    int coefficient = 0;
    for (int factor_shell = 0; factor_shell <= 5; ++factor_shell) {
      if (std::abs(shell - factor_shell) <= output
          && output <= shell + factor_shell) {
        coefficient += factor_coefficients[factor_shell];
      }
    }
    if (coefficient == 0) {
      continue;
    }
    Polynomial ratio{Integer(2 * shell + 1)};
    for (int index = 0; index < shell; ++index) {
      ratio = multiply(ratio, linear(1, -index));
    }
    for (int index = shell + 2; index <= denominator_top; ++index) {
      ratio = multiply(ratio, linear(1, index));
    }
    add_scaled(result, ratio, coefficient);
  }
  return result;
}

void verify_b2c3_radial_polynomials() {
  constexpr int denominator_top = 11;
  std::vector<Polynomial> shells;
  for (int output = 0; output <= 5; ++output) {
    shells.push_back(scaled_b23_shell(output, denominator_top));
  }
  Polynomial k30 = multiply(shells[0U], shells[3U]);
  add_scaled(k30, multiply(shells[0U], shells[4U]), 1);
  add_scaled(k30, multiply(shells[0U], shells[5U]), 1);
  add_scaled(k30, multiply(shells[1U], shells[4U]), -1);
  Polynomial k31 = multiply(shells[0U], shells[4U]);
  add_scaled(k31, multiply(shells[0U], shells[5U]), 1);
  add_scaled(k31, multiply(shells[1U], shells[2U]), 1);
  add_scaled(k31, multiply(shells[2U], shells[3U]), -1);

  Polynomial denominator{Integer(1)};
  for (int index = 2; index <= denominator_top; ++index) {
    denominator = multiply(denominator, linear(1, index));
  }
  for (const int fundamentals : {0, 1, 6, 16}) {
    const std::vector<Integer> profile = profile_for(fundamentals, 2, 3);
    const Integer ballot_zero =
        binomial(2 * fundamentals, fundamentals) / (fundamentals + 1);
    const Integer divisor = evaluate(denominator, fundamentals);
    if (evaluate(k30, fundamentals) * ballot_zero * ballot_zero
            != radial(profile, 0) * divisor * divisor
        || evaluate(k31, fundamentals) * ballot_zero * ballot_zero
               != radial(profile, 1) * divisor * divisor) {
      throw std::runtime_error("radial-polynomial reconstruction failed");
    }
  }
  const Polynomial k30_shift = shifted_coefficients(k30, 0);
  const Polynomial k31_shift = shifted_coefficients(k31, 0);
  if (!nonnegative_coefficients(k30_shift)
      || !nonnegative_coefficients(k31_shift)) {
    throw std::runtime_error("radial-polynomial sign certificate failed");
  }
  std::cout << "SU2_FUNDAMENTAL_TWO_ARBITRARY_B23_RADIAL"
            << " K30 coefficients=";
  print_polynomial(k30_shift);
  std::cout << " K31 coefficients=";
  print_polynomial(k31_shift);
  std::cout << " result=PASS_EXACT\n";
}

struct RatioMinimum {
  bool initialized = false;
  Integer numerator = 0;
  Integer denominator = 1;
  int first = 0;
  int second = 0;
};

void observe_ratio(RatioMinimum& minimum, const Integer& numerator,
                   const Integer& denominator, const int first,
                   const int second) {
  if (!minimum.initialized
      || numerator * minimum.denominator
             < minimum.numerator * denominator) {
    minimum = RatioMinimum{true, numerator, denominator, first, second};
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2
        && std::string{argv[1]} == "--b23-threshold-polynomials") {
      verify_b2c3_threshold_polynomials();
      return EXIT_SUCCESS;
    }
    if (argc == 2 && std::string{argv[1]} == "--b23-radial-polynomials") {
      verify_b2c3_radial_polynomials();
      return EXIT_SUCCESS;
    }
    if (argc == 5 && std::string{argv[1]} == "--case") {
      inspect_case(
          parse_nonnegative(argv[2], "fundamental count"),
          parse_positive(argv[3], "first arbitrary label"),
          parse_positive(argv[4], "second arbitrary label"));
      return EXIT_SUCCESS;
    }
    if (argc != 3) {
      throw std::invalid_argument(
          "usage: probe_su2_fundamental_two_arbitrary_third "
          "MAXIMUM_FUNDAMENTALS MAXIMUM_LABEL | "
          "--case FUNDAMENTALS FIRST SECOND | "
          "--b23-threshold-polynomials | --b23-radial-polynomials");
    }
    const int maximum_fundamentals =
        parse_nonnegative(argv[1], "maximum fundamentals");
    const int maximum_label = parse_positive(argv[2], "maximum label");

    std::size_t words = 0U;
    std::size_t zero_below = 0U;
    std::size_t one_below = 0U;
    std::size_t both_below = 0U;
    std::size_t mixed_both_below = 0U;
    std::size_t third_failures = 0U;
    bool first_both = true;
    int first_m = 0;
    int first_b = 0;
    int first_c = 0;
    Integer first_k30 = 0;
    Integer first_k31 = 0;
    int largest_mixed_m = -1;
    int largest_mixed_b = 0;
    int largest_mixed_c = 0;
    int largest_zero_below_m = -1;
    int largest_zero_below_b = 0;
    int largest_zero_below_c = 0;
    int largest_one_below_m = -1;
    int largest_one_below_b = 0;
    int largest_one_below_c = 0;
    int largest_mixed_zero_below_m = -1;
    int largest_mixed_zero_below_b = 0;
    int largest_mixed_zero_below_c = 0;
    int largest_mixed_one_below_m = -1;
    int largest_mixed_one_below_b = 0;
    int largest_mixed_one_below_c = 0;
    std::vector<RatioMinimum> mixed_minima(
        static_cast<std::size_t>(maximum_fundamentals + 1));
    for (int fundamentals = 0; fundamentals <= maximum_fundamentals;
         ++fundamentals) {
      for (int first = 1; first <= maximum_label; ++first) {
        for (int second = first; second <= maximum_label; ++second) {
          const std::vector<Integer> profile =
              profile_for(fundamentals, first, second);
          const Integer k30 = radial(profile, 0);
          const Integer k31 = radial(profile, 1);
          ++words;
          if (k30 < 0 || k31 < 0) {
            ++third_failures;
          }
          const bool zero_below_threshold = !zero_threshold(profile);
          const bool one_below_threshold = !one_threshold(profile);
          if (first > 1 && second > first) {
            observe_ratio(
                mixed_minima[static_cast<std::size_t>(fundamentals)],
                at(profile, 1), at(profile, 0), first, second);
          }
          zero_below += zero_below_threshold ? 1U : 0U;
          one_below += one_below_threshold ? 1U : 0U;
          if (zero_below_threshold && fundamentals > largest_zero_below_m) {
            largest_zero_below_m = fundamentals;
            largest_zero_below_b = first;
            largest_zero_below_c = second;
          }
          if (one_below_threshold && fundamentals > largest_one_below_m) {
            largest_one_below_m = fundamentals;
            largest_one_below_b = first;
            largest_one_below_c = second;
          }
          if (first > 1 && second > first) {
            if (zero_below_threshold
                && fundamentals > largest_mixed_zero_below_m) {
              largest_mixed_zero_below_m = fundamentals;
              largest_mixed_zero_below_b = first;
              largest_mixed_zero_below_c = second;
            }
            if (one_below_threshold
                && fundamentals > largest_mixed_one_below_m) {
              largest_mixed_one_below_m = fundamentals;
              largest_mixed_one_below_b = first;
              largest_mixed_one_below_c = second;
            }
          }
          if (zero_below_threshold && one_below_threshold) {
            ++both_below;
            if (first > 1 && second > first) {
              ++mixed_both_below;
              if (fundamentals > largest_mixed_m) {
                largest_mixed_m = fundamentals;
                largest_mixed_b = first;
                largest_mixed_c = second;
              }
            }
            if (first_both) {
              first_both = false;
              first_m = fundamentals;
              first_b = first;
              first_c = second;
              first_k30 = k30;
              first_k31 = k31;
            }
          }
        }
      }
    }
    std::cout << "SU2_FUNDAMENTAL_TWO_ARBITRARY_THIRD"
              << " maximum_fundamentals=" << maximum_fundamentals
              << " maximum_label=" << maximum_label
              << " words=" << words
              << " zero_below=" << zero_below
              << " one_below=" << one_below
              << " both_below=" << both_below
              << " mixed_both_below=" << mixed_both_below
              << " third_failures=" << third_failures;
    if (!first_both) {
      std::cout << " first_both_below=";
      print_word(first_m, first_b, first_c);
      std::cout << " first_K30=" << first_k30
                << " first_K31=" << first_k31;
    }
    if (largest_mixed_m >= 0) {
      std::cout << " largest_mixed_m=" << largest_mixed_m
                << " largest_mixed_labels=[" << largest_mixed_b
                << ',' << largest_mixed_c << ']';
    }
    if (largest_zero_below_m >= 0) {
      std::cout << " largest_zero_below_m=" << largest_zero_below_m
                << " largest_zero_below_labels=[" << largest_zero_below_b
                << ',' << largest_zero_below_c << ']';
    }
    if (largest_one_below_m >= 0) {
      std::cout << " largest_one_below_m=" << largest_one_below_m
                << " largest_one_below_labels=[" << largest_one_below_b
                << ',' << largest_one_below_c << ']';
    }
    if (largest_mixed_zero_below_m >= 0) {
      std::cout << " largest_mixed_zero_below_m="
                << largest_mixed_zero_below_m
                << " largest_mixed_zero_below_labels=["
                << largest_mixed_zero_below_b << ','
                << largest_mixed_zero_below_c << ']';
    }
    if (largest_mixed_one_below_m >= 0) {
      std::cout << " largest_mixed_one_below_m="
                << largest_mixed_one_below_m
                << " largest_mixed_one_below_labels=["
                << largest_mixed_one_below_b << ','
                << largest_mixed_one_below_c << ']';
    }
    std::cout << " result="
              << (third_failures == 0U ? "PASS_EXACT_BOX" : "FAIL")
              << '\n';
    for (int fundamentals = 0; fundamentals <= maximum_fundamentals;
         ++fundamentals) {
      const RatioMinimum& minimum =
          mixed_minima[static_cast<std::size_t>(fundamentals)];
      if (!minimum.initialized) {
        continue;
      }
      std::cout << "SU2_FUNDAMENTAL_TWO_ARBITRARY_THIRD_MIXED_MINIMUM"
                << " fundamentals=" << fundamentals
                << " labels=[" << minimum.first << ',' << minimum.second
                << "] c1=" << minimum.numerator
                << " c0=" << minimum.denominator << '\n';
    }
    return third_failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
