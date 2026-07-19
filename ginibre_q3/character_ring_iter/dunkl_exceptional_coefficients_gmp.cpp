#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
using Key = std::uint32_t;
using Poly = std::unordered_map<Key, mpq_class>;
using Vector = std::vector<mpq_class>;
using Matrix = std::vector<Vector>;

struct GroupData {
  int rank;
  int positive_roots;
  int A;
  int dimension;
};

const std::map<std::string, GroupData> kData = {
    {"F4", {4, 24, 36, 52}}, {"E6", {6, 36, 72, 78}},
    {"E7", {7, 63, 126, 133}}, {"E8", {8, 120, 240, 248}},
};

[[noreturn]] void fail(const std::string& message) { throw std::runtime_error(message); }

int exponent(Key key, int variable) { return (key >> (4 * variable)) & 15U; }
Key with_exponent(Key key, int variable, int value) {
  const Key mask = Key(15) << (4 * variable);
  return (key & ~mask) | (Key(value) << (4 * variable));
}

void add_term(Poly& polynomial, Key key, const mpq_class& value) {
  if (value == 0) return;
  const auto found = polynomial.find(key);
  if (found == polynomial.end()) {
    polynomial.emplace(key, value);
  } else {
    found->second += value;
    if (found->second == 0) polynomial.erase(found);
  }
}

void add_scaled(Poly& target, const Poly& source, const mpq_class& scale) {
  if (scale == 0) return;
  for (const auto& [key, coefficient] : source) add_term(target, key, scale * coefficient);
}

Poly multiply(const Poly& first, const Poly& second) {
  Poly out;
  out.reserve(first.size() * std::min<size_t>(second.size(), 32));
  for (const auto& [a, ca] : first) {
    for (const auto& [b, cb] : second) {
      Key key = 0;
      for (int i = 0; i < 8; ++i) {
        const int power = exponent(a, i) + exponent(b, i);
        if (power > 15) fail("packed monomial exponent overflow");
        key = with_exponent(key, i, power);
      }
      add_term(out, key, ca * cb);
    }
  }
  return out;
}

Poly linear_power(const Vector& alpha, int power) {
  Poly linear;
  for (int i = 0; i < static_cast<int>(alpha.size()); ++i)
    if (alpha[i] != 0) linear.emplace(with_exponent(0, i, 1), alpha[i]);
  Poly out{{0, mpq_class(1)}};
  for (int i = 0; i < power; ++i) out = multiply(out, linear);
  return out;
}

Poly root_power_sum(const std::vector<Vector>& roots, int power) {
  Poly out;
  for (const auto& root : roots) add_scaled(out, linear_power(root, power), 1);
  return out;
}

Poly ordinary_laplacian(const Poly& polynomial, int variables) {
  Poly out;
  for (const auto& [key, coefficient] : polynomial) {
    for (int i = 0; i < variables; ++i) {
      const int e = exponent(key, i);
      if (e >= 2)
        add_term(out, with_exponent(key, i, e - 2), coefficient * e * (e - 1));
    }
  }
  return out;
}

Poly directional_derivative(const Poly& polynomial, const Vector& alpha) {
  Poly out;
  for (const auto& [key, coefficient] : polynomial) {
    for (int i = 0; i < static_cast<int>(alpha.size()); ++i) {
      const int e = exponent(key, i);
      if (e != 0 && alpha[i] != 0)
        add_term(out, with_exponent(key, i, e - 1), coefficient * e * alpha[i]);
    }
  }
  return out;
}

Poly scaled(Poly polynomial, const mpq_class& scale) {
  if (scale == 0) return {};
  for (auto& [key, coefficient] : polynomial) coefficient *= scale;
  return polynomial;
}

Poly subtract(Poly first, const Poly& second) {
  add_scaled(first, second, -1);
  return first;
}

Poly multiply_nonpivot_linear(const Poly& polynomial, const Vector& alpha, int pivot) {
  Poly out;
  for (const auto& [key, coefficient] : polynomial) {
    for (int i = 0; i < static_cast<int>(alpha.size()); ++i) {
      if (i == pivot || alpha[i] == 0) continue;
      const int e = exponent(key, i);
      add_term(out, with_exponent(key, i, e + 1), coefficient * alpha[i]);
    }
  }
  return out;
}

Poly divide_by_linear(const Poly& polynomial, const Vector& alpha) {
  if (polynomial.empty()) return {};
  int pivot = -1;
  for (int i = 0; i < static_cast<int>(alpha.size()); ++i)
    if (alpha[i] != 0) {
      pivot = i;
      break;
    }
  if (pivot < 0) fail("zero root in linear division");
  int maximum = 0;
  for (const auto& [key, coefficient] : polynomial)
    maximum = std::max(maximum, exponent(key, pivot));
  if (maximum == 0) fail("nonzero polynomial is not divisible by pivot form");

  std::vector<Poly> buckets(maximum + 1);
  for (const auto& [key, coefficient] : polynomial) {
    const int degree = exponent(key, pivot);
    add_term(buckets[degree], with_exponent(key, pivot, 0), coefficient);
  }
  std::vector<Poly> quotient_buckets(maximum);
  Poly next;
  for (int degree = maximum; degree >= 1; --degree) {
    quotient_buckets[degree - 1] = scaled(subtract(buckets[degree], next), 1 / alpha[pivot]);
    next = multiply_nonpivot_linear(quotient_buckets[degree - 1], alpha, pivot);
  }
  if (!subtract(buckets[0], next).empty()) fail("nonzero root-linear division remainder");
  Poly quotient;
  for (int degree = 0; degree < maximum; ++degree)
    for (const auto& [key, coefficient] : quotient_buckets[degree])
      add_term(quotient, with_exponent(key, pivot, degree), coefficient);
  return quotient;
}

Poly dunkl_laplacian(const Poly& polynomial, const std::vector<Vector>& roots,
                     int variables) {
  Poly result = ordinary_laplacian(polynomial, variables);
  const int thread_count = std::max(1, omp_get_max_threads());
  std::vector<Poly> local(thread_count);
#pragma omp parallel for schedule(dynamic, 1)
  for (int index = 0; index < static_cast<int>(roots.size()); ++index) {
    Poly derivative = directional_derivative(polynomial, roots[index]);
    Poly quotient = divide_by_linear(derivative, roots[index]);
    add_scaled(local[omp_get_thread_num()], quotient, 2);
  }
  for (const auto& part : local) add_scaled(result, part, 1);
  return result;
}

mpq_class expectation(Poly polynomial, const std::vector<Vector>& roots, int degree,
                      const std::string& label) {
  if (degree & 1) return 0;
  const int steps = degree / 2;
  for (int step = 1; step <= steps; ++step) {
    std::cout << label << " Dunkl step " << step << '/' << steps
              << "; terms=" << polynomial.size() << '\n';
    polynomial = dunkl_laplacian(polynomial, roots, roots[0].size());
  }
  const auto found = polynomial.find(0);
  mpq_class constant = found == polynomial.end() ? mpq_class(0) : found->second;
  mpz_class denominator;
  mpz_ui_pow_ui(denominator.get_mpz_t(), 4, steps);
  mpz_class fact;
  mpz_fac_ui(fact.get_mpz_t(), steps);
  return constant / (denominator * fact);
}

Poly quadratic_form(const Matrix& matrix) {
  Poly out;
  const int n = matrix.size();
  for (int i = 0; i < n; ++i) {
    if (matrix[i][i] != 0) add_term(out, with_exponent(0, i, 2), matrix[i][i]);
    for (int j = i + 1; j < n; ++j) {
      const mpq_class coefficient = matrix[i][j] + matrix[j][i];
      if (coefficient != 0)
        add_term(out, with_exponent(with_exponent(0, i, 1), j, 1), coefficient);
    }
  }
  return out;
}

Poly quadratic_power(const Matrix& matrix, int power) {
  const Poly quadratic = quadratic_form(matrix);
  Poly out{{0, mpq_class(1)}};
  for (int i = 0; i < power; ++i) out = multiply(out, quadratic);
  return out;
}

Matrix identity(int n) {
  Matrix out(n, Vector(n));
  for (int i = 0; i < n; ++i) out[i][i] = 1;
  return out;
}

mpq_class inner(const Vector& a, const Vector& b) {
  mpq_class out = 0;
  for (size_t i = 0; i < a.size(); ++i) out += a[i] * b[i];
  return out;
}

struct VectorLess {
  bool operator()(const Vector& a, const Vector& b) const {
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
  }
};

std::vector<Vector> generate_roots(const std::vector<Vector>& simple) {
  std::map<Vector, bool, VectorLess> seen;
  std::queue<Vector> frontier;
  for (const auto& root : simple) {
    seen.emplace(root, true);
    frontier.push(root);
  }
  while (!frontier.empty()) {
    const Vector beta = frontier.front();
    frontier.pop();
    for (const auto& alpha : simple) {
      const mpq_class coefficient = 2 * inner(beta, alpha) / inner(alpha, alpha);
      Vector reflected(beta.size());
      for (size_t i = 0; i < beta.size(); ++i) reflected[i] = beta[i] - coefficient * alpha[i];
      if (seen.emplace(reflected, true).second) frontier.push(std::move(reflected));
    }
  }
  std::vector<Vector> positive;
  for (const auto& [root, unused] : seen) {
    for (const auto& coordinate : root) {
      if (coordinate == 0) continue;
      if (coordinate > 0) positive.push_back(root);
      break;
    }
  }
  return positive;
}

std::vector<Vector> exceptional_simple(int rank) {
  const mpq_class h(1, 2), mh(-1, 2);
  const Vector a1{h, mh, mh, mh, mh, mh, mh, h};
  const Vector a2{1, 1, 0, 0, 0, 0, 0, 0};
  const Vector a3{-1, 1, 0, 0, 0, 0, 0, 0};
  const Vector a4{0, -1, 1, 0, 0, 0, 0, 0};
  const Vector a5{0, 0, -1, 1, 0, 0, 0, 0};
  const Vector a6{0, 0, 0, -1, 1, 0, 0, 0};
  const Vector a7{0, 0, 0, 0, -1, 1, 0, 0};
  return rank == 6 ? std::vector<Vector>{a1, a2, a3, a4, a5, a6}
                   : std::vector<Vector>{a1, a2, a3, a4, a5, a6, a7};
}

std::vector<Vector> f4_roots() {
  std::vector<Vector> roots;
  for (int i = 0; i < 4; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      for (int sign : {1, -1}) {
        Vector root(4);
        root[i] = 1;
        root[j] = sign;
        roots.push_back(std::move(root));
      }
    }
  }
  for (int i = 0; i < 4; ++i) {
    Vector root(4);
    root[i] = 1;
    roots.push_back(std::move(root));
  }
  for (int s1 : {1, -1})
    for (int s2 : {1, -1})
      for (int s3 : {1, -1}) roots.push_back(Vector{mpq_class(1, 2), mpq_class(s1, 2),
                                                     mpq_class(s2, 2), mpq_class(s3, 2)});
  return roots;
}

std::vector<Vector> e8_roots() {
  std::vector<Vector> full;
  for (int i = 0; i < 8; ++i)
    for (int j = i + 1; j < 8; ++j)
      for (int si : {1, -1})
        for (int sj : {1, -1}) {
          Vector root(8);
          root[i] = si;
          root[j] = sj;
          full.push_back(std::move(root));
        }
  for (int mask = 0; mask < 256; ++mask) {
    if (__builtin_popcount(static_cast<unsigned>(mask)) & 1) continue;
    Vector root(8);
    for (int i = 0; i < 8; ++i) root[i] = (mask & (1 << i)) ? mpq_class(-1, 2) : mpq_class(1, 2);
    full.push_back(std::move(root));
  }
  std::vector<Vector> positive;
  for (const auto& root : full) {
    for (const auto& coordinate : root) {
      if (coordinate == 0) continue;
      if (coordinate > 0) positive.push_back(root);
      break;
    }
  }
  return positive;
}

Matrix invert(Matrix matrix) {
  const int n = matrix.size();
  Matrix augmented(n, Vector(2 * n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) augmented[i][j] = matrix[i][j];
    augmented[i][n + i] = 1;
  }
  for (int column = 0; column < n; ++column) {
    int pivot = column;
    while (pivot < n && augmented[pivot][column] == 0) ++pivot;
    if (pivot == n) fail("singular Gram matrix");
    std::swap(augmented[pivot], augmented[column]);
    const mpq_class scale = augmented[column][column];
    for (auto& value : augmented[column]) value /= scale;
    for (int row = 0; row < n; ++row) {
      if (row == column || augmented[row][column] == 0) continue;
      const mpq_class factor = augmented[row][column];
      for (int j = 0; j < 2 * n; ++j) augmented[row][j] -= factor * augmented[column][j];
    }
  }
  Matrix out(n, Vector(n));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) out[i][j] = augmented[i][n + j];
  return out;
}

Matrix projector(const std::vector<Vector>& simple) {
  const int rank = simple.size(), n = simple[0].size();
  Matrix gram(rank, Vector(rank));
  for (int i = 0; i < rank; ++i)
    for (int j = 0; j < rank; ++j) gram[i][j] = inner(simple[i], simple[j]);
  const Matrix inverse = invert(std::move(gram));
  Matrix out(n, Vector(n));
  for (int a = 0; a < n; ++a)
    for (int b = 0; b < n; ++b)
      for (int i = 0; i < rank; ++i)
        for (int j = 0; j < rank; ++j)
          out[a][b] += simple[i][a] * inverse[i][j] * simple[j][b];
  return out;
}

mpq_class gamma_moment(const std::string& group, int k) {
  mpq_class out = 1;
  const mpq_class a(kData.at(group).dimension, 2);
  for (int i = 0; i < k; ++i) out *= a + i;
  return out;
}

mpq_class root_q_moment(const std::vector<Vector>& roots, const Matrix& projection,
                        int root_power, int q_power, const std::string& label) {
  Poly polynomial = root_power_sum(roots, root_power);
  if (q_power) polynomial = multiply(polynomial, quadratic_power(projection, q_power));
  return expectation(std::move(polynomial), roots, root_power + 2 * q_power, label);
}

mpq_class c2_from_r(const std::string& group, const mpq_class& r_moment) {
  const auto data = kData.at(group);
  const mpq_class d = data.dimension, cp(data.A, data.rank), a(data.dimension, 2);
  const mpq_class coef_r = -d * d / (360 * cp * cp * cp);
  const mpq_class coef_w6 = 5 * d / (12 * (d + 2)) - mpq_class(1, 3) +
                            (d + 12) * d / (144 * (d + 2));
  const mpq_class coef_w4 = d * d * (d + 1) / (288 * (d + 2));
  const mpq_class coef_w8 = (d + 12) * (d + 12) / (288 * (d + 2) * (d + 2));
  return coef_r * r_moment + coef_w6 * gamma_moment(group, 3) +
         coef_w4 * gamma_moment(group, 2) + coef_w8 * gamma_moment(group, 4);
}

mpq_class c3_from_moments(const std::string& group, const mpq_class& r0,
                          const mpq_class& r1, const mpq_class& r2,
                          const mpq_class& s0) {
  const mpq_class d = kData.at(group).dimension;
  auto q = [&](int k) -> mpq_class { return gamma_moment(group, k); };
  auto p_q = [&](int k) -> mpq_class {
    const mpq_class moment = q(k + 2);
    return mpq_class(mpq_class(5, kData.at(group).dimension + 2) * moment);
  };
  auto pj_q = [&](int j, int k) -> mpq_class {
    const mpq_class factor = [&] {
      mpq_class factor = 1;
      for (int i = 0; i < j; ++i) factor *= mpq_class(5, kData.at(group).dimension + 2);
      return factor;
    }();
    const mpq_class moment = q(2 * j + k);
    return mpq_class(moment * factor);
  };
  const mpq_class f(5, kData.at(group).dimension + 2);
  const mpq_class d2 = d * d, d3 = d2 * d;
  const mpq_class f2 = f * f, f3 = f2 * f;
  const mpq_class q3 = q(3), q4 = q(4), q5 = q(5), q6 = q(6);
  const mpq_class p1 = p_q(1), p2 = p_q(2), pp2 = pj_q(2, 0);
  mpq_class total = 0;
  total += d3 * s0 / 20160;
  total -= d2 * r1 / 360;
  total -= d2 * pp2 / 288;
  total += d * p2 / 12;
  total -= q4 / 4;
  total -= d3 * f * r2 / 4320;
  total += d2 * f2 * q5 / 144;
  total -= d * f * q5 / 36;
  total += d2 * r2 / 720;
  total -= d * f * q5 / 24;
  total += q5 / 6;
  total += d3 * f3 * q6 / (1728 * 6);
  total -= 3 * d2 * f2 * q6 / (288 * 6);
  total += 3 * d * f * q6 / (48 * 6);
  total -= q6 / (8 * 6);
  total += d3 * r1 / (360 * 12);
  total -= d2 * p2 / 144;
  total += d * q4 / 36;
  total -= d3 * f2 * q5 / (144 * 12 * 2);
  total += d2 * f * q5 / (144 * 2);
  total -= d * q5 / (48 * 2);
  total += d3 * p2 / (12 * 288);
  total -= d3 * f2 * q4 / (12 * 1440);
  total -= d2 * q4 / (2 * 288);
  total += d2 * p2 / (2 * 1440);
  total -= d3 * q3 / 10368;
  total += d3 * p1 / 17280;
  total -= d3 * r0 / 60480;
  return total;
}

void verify_group(const std::string& group, const std::vector<Vector>& roots,
                  const Matrix& projection, const mpq_class& expected_c2,
                  const mpq_class& expected_c3) {
  const auto start = Clock::now();
  if (static_cast<int>(roots.size()) != kData.at(group).positive_roots)
    fail(group + " positive-root count mismatch");
  mpq_class norm_sum = 0;
  for (const auto& root : roots) norm_sum += inner(root, root);
  if (norm_sum != kData.at(group).A) fail(group + " root norm sum mismatch");
  const mpq_class radial = expectation(quadratic_power(projection, 3), roots, 6,
                                       group + " radial6");
  if (radial != gamma_moment(group, 3)) fail(group + " radial normalization failed");
  const mpq_class r0w = root_q_moment(roots, projection, 6, 0, group + " R");
  const mpq_class r1w = root_q_moment(roots, projection, 6, 1, group + " RQ");
  const mpq_class r2w = root_q_moment(roots, projection, 6, 2, group + " RQ2");
  const mpq_class s0w = root_q_moment(roots, projection, 8, 0, group + " S");
  const mpq_class cp(kData.at(group).A, kData.at(group).rank);
  const mpq_class c2 = c2_from_r(group, r0w);
  const mpq_class c3 = c3_from_moments(group, r0w / (cp * cp * cp),
                                       r1w / (cp * cp * cp), r2w / (cp * cp * cp),
                                       s0w / (cp * cp * cp * cp));
  if (c2 != expected_c2) fail(group + " c2 mismatch: " + c2.get_str());
  if (c3 != expected_c3) fail(group + " c3 mismatch: " + c3.get_str());
  std::cout << group << " exact coefficients: c2=" << c2 << " c3=" << c3
            << "; elapsed=" << std::chrono::duration<double>(Clock::now() - start).count()
            << " s\n";
}

void verify_f4() {
  const auto roots = f4_roots();
  const auto projection = identity(4);
  const mpq_class s = root_q_moment(roots, projection, 8, 0, "F4 S normalization");
  if (s != 6265350) fail("F4 Dunkl/Wick S normalization mismatch");
  verify_group("F4", roots, projection, mpq_class(8401705, 486),
               mpq_class(mpz_class("-44918139383"), mpz_class(39366)));
}

void verify_e8() {
  const auto roots = e8_roots();
  if (roots.size() != 120) fail("E8 positive-root count mismatch");
  const auto projection = identity(8);
  const mpq_class radial = expectation(quadratic_power(projection, 4), roots, 8,
                                       "E8 radial8");
  if (radial != gamma_moment("E8", 4)) fail("E8 radial normalization failed");
  const mpq_class s0w = root_q_moment(roots, projection, 8, 0, "E8 S");
  const mpq_class cp(kData.at("E8").A, kData.at("E8").rank);
  const mpq_class r0 = 15 * gamma_moment("E8", 3) / (cp * cp * cp);
  const mpq_class r1 = 15 * gamma_moment("E8", 4) / (cp * cp * cp);
  const mpq_class r2 = 15 * gamma_moment("E8", 5) / (cp * cp * cp);
  const mpq_class c2 = c2_from_r("E8", r0 * cp * cp * cp);
  const mpq_class c3 = c3_from_moments("E8", r0, r1, r2,
                                       s0w / (cp * cp * cp * cp));
  if (c2 != mpq_class(mpz_class("346377539"), mpz_class(45))) fail("E8 c2 mismatch");
  if (c3 != mpq_class(mpz_class("-11478332145838"), mpz_class(1125)))
    fail("E8 c3 mismatch: " + c3.get_str());
  std::cout << "E8 exact coefficients: c2=" << c2 << " c3=" << c3 << '\n';
}
}  // namespace

int main(int argc, char** argv) {
  try {
    int threads = omp_get_max_threads();
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--threads" && i + 1 < argc)
        threads = std::stoi(argv[++i]);
      else
        fail("usage: dunkl_exceptional_coefficients_gmp [--threads N]");
    }
    omp_set_num_threads(std::max(1, threads));
    const auto start = Clock::now();
    verify_f4();
    const auto e6_simple = exceptional_simple(6);
    verify_group("E6", generate_roots(e6_simple), projector(e6_simple),
                 mpq_class(mpz_class("10524215"), mpz_class(128)),
                 mpq_class(mpz_class("-11872316417"), mpz_class(1024)));
    const auto e7_simple = exceptional_simple(7);
    verify_group("E7", generate_roots(e7_simple), projector(e7_simple),
                 mpq_class(mpz_class("82042395295"), mpz_class(124416)),
                 mpq_class(mpz_class("-41730554417849633"), mpz_class(161243136)));
    verify_e8();
    std::cout << "All exceptional Dunkl coefficient certificates verified with exact GMP "
              << "arithmetic; OpenMP threads=" << omp_get_max_threads() << ", elapsed="
              << std::chrono::duration<double>(Clock::now() - start).count() << " s.\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "Exceptional Dunkl GMP certificate failure: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
