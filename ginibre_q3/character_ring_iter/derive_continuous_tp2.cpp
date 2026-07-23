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
  return nonnegative ? 0 : 1;
}
