#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <iostream>
#include <map>
#include <utility>

using boost::multiprecision::cpp_int;

namespace {

class Polynomial {
 public:
  Polynomial() = default;
  explicit Polynomial(int constant) {
    if (constant != 0) {
      coefficients_[{0, 0}] = constant;
    }
  }

  static Polynomial variable(int index) {
    Polynomial result;
    result.coefficients_[index == 0 ? std::pair<int, int>{1, 0}
                                    : std::pair<int, int>{0, 1}] = 1;
    return result;
  }

  Polynomial& operator+=(const Polynomial& other) {
    for (const auto& [degree, coefficient] : other.coefficients_) {
      coefficients_[degree] += coefficient;
      if (coefficients_[degree] == 0) {
        coefficients_.erase(degree);
      }
    }
    return *this;
  }

  Polynomial& operator-=(const Polynomial& other) {
    for (const auto& [degree, coefficient] : other.coefficients_) {
      coefficients_[degree] -= coefficient;
      if (coefficients_[degree] == 0) {
        coefficients_.erase(degree);
      }
    }
    return *this;
  }

  Polynomial& operator*=(const Polynomial& other) {
    std::map<std::pair<int, int>, cpp_int> product;
    for (const auto& [left_degree, left_coefficient] : coefficients_) {
      for (const auto& [right_degree, right_coefficient] :
           other.coefficients_) {
        const std::pair<int, int> degree{
            left_degree.first + right_degree.first,
            left_degree.second + right_degree.second};
        product[degree] += left_coefficient * right_coefficient;
      }
    }
    coefficients_ = std::move(product);
    return *this;
  }

  const auto& coefficients() const { return coefficients_; }

 private:
  std::map<std::pair<int, int>, cpp_int> coefficients_;
};

Polynomial operator+(Polynomial left, const Polynomial& right) {
  left += right;
  return left;
}

Polynomial operator-(Polynomial left, const Polynomial& right) {
  left -= right;
  return left;
}

Polynomial operator*(Polynomial left, const Polynomial& right) {
  left *= right;
  return left;
}

Polynomial operator*(int scalar, Polynomial value) {
  value *= Polynomial(scalar);
  return value;
}

Polynomial square(const Polynomial& value) { return value * value; }

}  // namespace

int main() {
  const Polynomial u = Polynomial::variable(0);
  const Polynomial v = Polynomial::variable(1);
  const Polynomial s = Polynomial(7) + 2 * u + v;
  const Polynomial d = Polynomial(2) + v;
  const Polynomial a = 2 * square(d) - square(s) + s;
  const Polynomial b = 3 * square(s) - 3 * s - Polynomial(2) -
                       4 * square(d);
  const Polynomial numerator =
      6 * (s - Polynomial(2)) * square(s + Polynomial(2)) * square(a) +
      10 * square(d) * (s + Polynomial(2)) * a * b +
      4 * (s + Polynomial(3)) * square(d) * square(b);

  bool nonnegative = true;
  std::cout << "CONTINUOUS_TP2_NUMERATOR terms="
            << numerator.coefficients().size() << '\n';
  for (const auto& [degree, coefficient] : numerator.coefficients()) {
    std::cout << "  u^" << degree.first << " v^" << degree.second
              << " coefficient=" << coefficient << '\n';
    if (coefficient < 0) {
      nonnegative = false;
    }
  }
  std::cout << "CONTINUOUS_TP2_NUMERATOR result="
            << (nonnegative ? "PASS" : "FAIL") << '\n';

  const auto numerator_for = [](const Polynomial& s_value,
                                const Polynomial& d_value) {
    const Polynomial a_value = 2 * square(d_value) - square(s_value)
                               + s_value;
    const Polynomial b_value = 3 * square(s_value) - 3 * s_value
                               - Polynomial(2) - 4 * square(d_value);
    return 6 * (s_value - Polynomial(2))
               * square(s_value + Polynomial(2)) * square(a_value)
           + 10 * square(d_value) * (s_value + Polynomial(2))
               * a_value * b_value
           + 4 * (s_value + Polynomial(3)) * square(d_value)
               * square(b_value);
  };
  const auto denominator_for = [](const Polynomial& s_value) {
    return square(s_value) * square(s_value + Polynomial(1))
           * square(s_value + Polynomial(2)) * (s_value + Polynomial(3))
           * (s_value - Polynomial(2));
  };
  const Polynomial n0 = numerator_for(s, d);
  const Polynomial d0 = denominator_for(s);
  const Polynomial n_minus = numerator_for(
      s + Polynomial(1), d + Polynomial(1)
  );
  const Polynomial d_minus = denominator_for(s + Polynomial(1));
  const Polynomial n_plus = numerator_for(
      s + Polynomial(1), d - Polynomial(1)
  );
  const Polynomial d_plus = denominator_for(s + Polynomial(1));
  const Polynomial cross_numerator =
      16 * n0 * d_minus * d_plus
      - n_minus * d0 * d_plus
      - n_plus * d0 * d_minus;
  bool cross_nonnegative = true;
  cpp_int minimum_cross_coefficient = 0;
  bool first_cross_coefficient = true;
  std::cout << "CONTINUOUS_TP2_PASCAL_CROSS terms="
            << cross_numerator.coefficients().size() << '\n';
  for (const auto& [degree, coefficient] :
       cross_numerator.coefficients()) {
    (void)degree;
    if (first_cross_coefficient
        || coefficient < minimum_cross_coefficient) {
      minimum_cross_coefficient = coefficient;
      first_cross_coefficient = false;
    }
    if (coefficient < 0) {
      cross_nonnegative = false;
    }
  }
  std::cout << "CONTINUOUS_TP2_PASCAL_CROSS minimum_coefficient="
            << minimum_cross_coefficient << '\n';
  std::cout << "CONTINUOUS_TP2_PASCAL_CROSS result="
            << (cross_nonnegative ? "PASS" : "FAIL") << '\n';

  const Polynomial symmetric_s = Polynomial(5) + u + v;
  const Polynomial symmetric_d = u - v;
  const Polynomial symmetric_n0 = numerator_for(symmetric_s, symmetric_d);
  const Polynomial symmetric_d0 = denominator_for(symmetric_s);
  const Polynomial symmetric_n_minus = numerator_for(
      symmetric_s + Polynomial(1), symmetric_d + Polynomial(1)
  );
  const Polynomial symmetric_d_minus = denominator_for(
      symmetric_s + Polynomial(1)
  );
  const Polynomial symmetric_n_plus = numerator_for(
      symmetric_s + Polynomial(1), symmetric_d - Polynomial(1)
  );
  const Polynomial symmetric_d_plus = denominator_for(
      symmetric_s + Polynomial(1)
  );
  const Polynomial symmetric_cross_numerator =
      16 * symmetric_n0 * symmetric_d_minus * symmetric_d_plus
      - symmetric_n_minus * symmetric_d0 * symmetric_d_plus
      - symmetric_n_plus * symmetric_d0 * symmetric_d_minus;
  bool symmetric_cross_nonnegative = true;
  cpp_int symmetric_minimum_cross_coefficient = 0;
  bool first_symmetric_cross_coefficient = true;
  std::cout << "CONTINUOUS_TP2_PASCAL_CROSS_SYMMETRIC terms="
            << symmetric_cross_numerator.coefficients().size() << '\n';
  for (const auto& [degree, coefficient] :
       symmetric_cross_numerator.coefficients()) {
    (void)degree;
    if (first_symmetric_cross_coefficient
        || coefficient < symmetric_minimum_cross_coefficient) {
      symmetric_minimum_cross_coefficient = coefficient;
      first_symmetric_cross_coefficient = false;
    }
    if (coefficient < 0) {
      symmetric_cross_nonnegative = false;
    }
  }
  std::cout << "CONTINUOUS_TP2_PASCAL_CROSS_SYMMETRIC "
               "minimum_coefficient="
            << symmetric_minimum_cross_coefficient << '\n';
  std::cout << "CONTINUOUS_TP2_PASCAL_CROSS_SYMMETRIC result="
            << (symmetric_cross_nonnegative ? "PASS" : "FAIL") << '\n';
  return nonnegative && cross_nonnegative && symmetric_cross_nonnegative
      ? 0 : 1;
}
