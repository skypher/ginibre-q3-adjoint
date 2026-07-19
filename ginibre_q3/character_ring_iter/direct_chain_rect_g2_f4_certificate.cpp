#include <gmpxx.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Exponent = std::array<int, 4>;
using Polynomial = std::map<Exponent, mpq_class>;
using Vector = std::vector<mpq_class>;

struct Rectangle {
  mpq_class s0;
  mpq_class s1;
  mpq_class t0;
  mpq_class t1;
};

struct GroupData {
  std::string name;
  int rank;
  int positive_roots;
  int dimension;
  int negative_constant;
  int kappa;
  std::vector<unsigned long> degrees;
  mpq_class quartic_denominator;
  int tail_start;
  Rectangle rectangle;
  mpq_class expected_discriminant_scale;
  mpq_class expected_coroot_covolume_squared;
  mpz_class expected_center_index;
  mpq_class expected_weyl_lattice_factor;
  mpq_class expected_gap;
  mpq_class expected_h;
  mpz_class base_threshold;
  mpz_class step_threshold;
};

[[noreturn]] void fail(const std::string& message) {
  throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
  if (!condition) fail(message);
}

mpq_class qpow(mpq_class base, unsigned long exponent) {
  mpq_class out = 1;
  while (exponent != 0) {
    if (exponent & 1UL) out *= base;
    exponent >>= 1UL;
    if (exponent != 0) base *= base;
  }
  return out;
}

mpz_class factorial(unsigned long n) {
  mpz_class out;
  mpz_fac_ui(out.get_mpz_t(), n);
  return out;
}

unsigned long ceiling(const mpq_class& value) {
  require(value >= 0, "ceiling called on a negative rational");
  mpz_class quotient;
  mpz_cdiv_q(quotient.get_mpz_t(), value.get_num_mpz_t(), value.get_den_mpz_t());
  return quotient.get_ui();
}

mpq_class inner(const Vector& first, const Vector& second) {
  require(first.size() == second.size(), "vector-size mismatch");
  mpq_class out = 0;
  for (size_t i = 0; i < first.size(); ++i) out += first[i] * second[i];
  return out;
}

mpq_class determinant(std::vector<Vector> matrix) {
  require(!matrix.empty(), "determinant called on an empty matrix");
  const std::size_t size = matrix.size();
  for (const Vector& row : matrix)
    require(row.size() == size, "determinant matrix is not square");
  mpq_class out = 1;
  for (std::size_t column = 0; column < size; ++column) {
    std::size_t pivot = column;
    while (pivot < size && matrix[pivot][column] == 0) ++pivot;
    require(pivot < size, "singular determinant matrix");
    if (pivot != column) {
      std::swap(matrix[pivot], matrix[column]);
      out = -out;
    }
    const mpq_class pivot_value = matrix[column][column];
    out *= pivot_value;
    for (std::size_t row = column + 1; row < size; ++row) {
      const mpq_class multiplier = matrix[row][column] / pivot_value;
      for (std::size_t entry = column; entry < size; ++entry)
        matrix[row][entry] -= multiplier * matrix[column][entry];
    }
  }
  out.canonicalize();
  return out;
}

mpq_class gram_determinant(const std::vector<Vector>& vectors) {
  std::vector<Vector> gram(vectors.size(), Vector(vectors.size()));
  for (std::size_t i = 0; i < vectors.size(); ++i)
    for (std::size_t j = 0; j < vectors.size(); ++j)
      gram[i][j] = inner(vectors[i], vectors[j]);
  return determinant(std::move(gram));
}

std::vector<Vector> inverse_matrix(const std::vector<Vector>& input) {
  require(!input.empty(), "inverse called on an empty matrix");
  const std::size_t size = input.size();
  std::vector<Vector> augmented(size, Vector(2 * size));
  for (std::size_t i = 0; i < size; ++i) {
    require(input[i].size() == size, "inverse matrix is not square");
    for (std::size_t j = 0; j < size; ++j) augmented[i][j] = input[i][j];
    augmented[i][size + i] = 1;
  }
  for (std::size_t column = 0; column < size; ++column) {
    std::size_t pivot = column;
    while (pivot < size && augmented[pivot][column] == 0) ++pivot;
    require(pivot < size, "singular inverse matrix");
    if (pivot != column) std::swap(augmented[pivot], augmented[column]);
    const mpq_class pivot_value = augmented[column][column];
    for (mpq_class& entry : augmented[column]) entry /= pivot_value;
    for (std::size_t row = 0; row < size; ++row) {
      if (row == column) continue;
      const mpq_class multiplier = augmented[row][column];
      for (std::size_t entry = 0; entry < 2 * size; ++entry)
        augmented[row][entry] -= multiplier * augmented[column][entry];
    }
  }
  std::vector<Vector> inverse(size, Vector(size));
  for (std::size_t i = 0; i < size; ++i)
    for (std::size_t j = 0; j < size; ++j)
      inverse[i][j] = augmented[i][size + j];
  return inverse;
}

mpq_class minimum_dual_lattice_norm_squared(const std::vector<Vector>& simple_roots) {
  std::vector<Vector> root_gram(simple_roots.size(), Vector(simple_roots.size()));
  for (std::size_t i = 0; i < simple_roots.size(); ++i)
    for (std::size_t j = 0; j < simple_roots.size(); ++j)
      root_gram[i][j] = inner(simple_roots[i], simple_roots[j]);
  const std::vector<Vector> dual_gram = inverse_matrix(root_gram);
  unsigned long cases = 1;
  for (std::size_t i = 0; i < simple_roots.size(); ++i) cases *= 3;
  mpq_class minimum = 2;
  std::vector<int> coordinates(simple_roots.size());
  for (unsigned long code = 0; code < cases; ++code) {
    unsigned long cursor = code;
    bool nonzero = false;
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
      coordinates[i] = static_cast<int>(cursor % 3) - 1;
      cursor /= 3;
      nonzero = nonzero || coordinates[i] != 0;
    }
    if (!nonzero) continue;
    mpq_class value = 0;
    for (std::size_t i = 0; i < coordinates.size(); ++i)
      for (std::size_t j = 0; j < coordinates.size(); ++j)
        value += coordinates[i] * coordinates[j] * dual_gram[i][j];
    if (value < minimum) minimum = value;
  }
  // Outside the ternary box some simple-root pairing has absolute value
  // at least two.  Since every simple root has squared length at most two,
  // Cauchy's inequality gives norm squared at least two.
  return minimum;
}

mpz_class cartan_determinant(const std::vector<Vector>& simple_roots,
                             const std::vector<Vector>& simple_coroots) {
  require(simple_roots.size() == simple_coroots.size(),
          "simple-root/coroot size mismatch");
  std::vector<Vector> cartan(simple_roots.size(), Vector(simple_roots.size()));
  for (std::size_t i = 0; i < simple_roots.size(); ++i)
    for (std::size_t j = 0; j < simple_roots.size(); ++j)
      cartan[i][j] = inner(simple_roots[i], simple_coroots[j]);
  mpq_class value = determinant(std::move(cartan));
  if (value < 0) value = -value;
  require(value.get_den() == 1, "Cartan determinant is not integral");
  return value.get_num();
}

void add_term(Polynomial& polynomial, const Exponent& exponent,
              const mpq_class& coefficient) {
  if (coefficient == 0) return;
  polynomial[exponent] += coefficient;
  if (polynomial[exponent] == 0) polynomial.erase(exponent);
}

Polynomial linear_power_four(const Vector& linear, int variables) {
  Polynomial out;
  for (int i = 0; i < variables; ++i)
    for (int j = 0; j < variables; ++j)
      for (int k = 0; k < variables; ++k)
        for (int l = 0; l < variables; ++l) {
          Exponent exponent{};
          ++exponent[i];
          ++exponent[j];
          ++exponent[k];
          ++exponent[l];
          add_term(out, exponent, linear[i] * linear[j] * linear[k] * linear[l]);
        }
  return out;
}

Polynomial quadratic_root_sum(const std::vector<Vector>& roots, int variables) {
  Polynomial out;
  for (const auto& root : roots)
    for (int i = 0; i < variables; ++i)
      for (int j = 0; j < variables; ++j) {
        Exponent exponent{};
        ++exponent[i];
        ++exponent[j];
        add_term(out, exponent, root[i] * root[j]);
      }
  return out;
}

Polynomial multiply(const Polynomial& first, const Polynomial& second) {
  Polynomial out;
  for (const auto& [a, ca] : first)
    for (const auto& [b, cb] : second) {
      Exponent exponent{};
      for (int i = 0; i < 4; ++i) exponent[i] = a[i] + b[i];
      add_term(out, exponent, ca * cb);
    }
  return out;
}

void verify_quartic_identity(const GroupData& data,
                             const std::vector<Vector>& coordinate_roots,
                             int variables) {
  Polynomial quartic;
  for (const auto& root : coordinate_roots) {
    const Polynomial fourth = linear_power_four(root, variables);
    for (const auto& [exponent, coefficient] : fourth)
      add_term(quartic, exponent, data.quartic_denominator * coefficient);
  }
  const Polynomial quadratic = quadratic_root_sum(coordinate_roots, variables);
  const Polynomial square = multiply(quadratic, quadratic);
  require(quartic == square, data.name + " quartic root identity failed");
}

std::vector<Vector> g2_metric_roots() {
  return {
      {1, -1, 0}, {1, 0, -1}, {0, 1, -1},
      {mpq_class(2, 3), mpq_class(-1, 3), mpq_class(-1, 3)},
      {mpq_class(1, 3), mpq_class(-2, 3), mpq_class(1, 3)},
      {mpq_class(1, 3), mpq_class(1, 3), mpq_class(-2, 3)},
  };
}

std::vector<Vector> g2_plane_forms() {
  // Restrict v=(x,y,-x-y) to the G2 root plane.
  return {{1, -1}, {2, 1}, {1, 2}, {1, 0}, {0, -1}, {1, 1}};
}

std::vector<Vector> g2_simple_roots() {
  return {{1, -1, 0},
          {mpq_class(-2, 3), mpq_class(1, 3), mpq_class(1, 3)}};
}

std::vector<Vector> g2_simple_coroots() {
  return {{1, -1, 0}, {-2, 1, 1}};
}

std::vector<Vector> f4_roots() {
  std::vector<Vector> roots;
  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j)
      for (int sign : {1, -1}) {
        Vector root(4);
        root[i] = 1;
        root[j] = sign;
        roots.push_back(std::move(root));
      }
  for (int i = 0; i < 4; ++i) {
    Vector root(4);
    root[i] = 1;
    roots.push_back(std::move(root));
  }
  for (int s1 : {1, -1})
    for (int s2 : {1, -1})
      for (int s3 : {1, -1})
        roots.push_back(Vector{mpq_class(1, 2), mpq_class(s1, 2),
                               mpq_class(s2, 2), mpq_class(s3, 2)});
  return roots;
}

std::vector<Vector> f4_simple_roots() {
  return {{0, 1, -1, 0},
          {0, 0, 1, -1},
          {0, 0, 0, 1},
          {mpq_class(1, 2), mpq_class(-1, 2), mpq_class(-1, 2),
           mpq_class(-1, 2)}};
}

std::vector<Vector> f4_simple_coroots() {
  return {{0, 1, -1, 0}, {0, 0, 1, -1}, {0, 0, 0, 2}, {1, -1, -1, -1}};
}

mpq_class discriminant_scale(const std::vector<Vector>& metric_roots) {
  // The Macdonald--Mehta discriminant uses a length-squared-two normal for
  // each reflection hyperplane.  Replacing those normals by the actual roots
  // multiplies Delta^2 by product |alpha|^2/2.
  mpq_class out = 1;
  for (const auto& root : metric_roots) out *= inner(root, root) / 2;
  return out;
}

mpq_class a0_prefactor(const GroupData& data, const mpq_class& scale_squared,
                       const mpq_class& weyl_lattice_factor) {
  mpz_class degree_product = 1;
  for (unsigned long degree : data.degrees) degree_product *= factorial(degree - 1);
  const mpq_class a(data.dimension, 2);
  return 8 * a * data.dimension * data.dimension *
         qpow(mpq_class(data.dimension, data.kappa), data.dimension) *
         degree_product * degree_product * scale_squared * weyl_lattice_factor /
         qpow(mpq_class(2), data.rank);
}

mpq_class rectangle_mass_lower(const GroupData& data, const mpq_class& gap) {
  const Rectangle& box = data.rectangle;
  const unsigned long shape = static_cast<unsigned long>(data.dimension / 2);
  const mpq_class s_integral =
      (qpow(box.s1, shape) - qpow(box.s0, shape)) / shape;
  const mpq_class t_integral =
      (qpow(box.t1, shape) - qpow(box.t0, shape)) / shape;
  return gap * gap * qpow(mpq_class(1, 3), ceiling(box.s1 + box.t1)) *
         s_integral * t_integral /
         (factorial(shape - 1) * factorial(shape - 1) * shape);
}

void verify_group(const GroupData& data, const std::vector<Vector>& metric_roots,
                  const std::vector<Vector>& coordinate_roots, int variables,
                  const std::vector<Vector>& simple_roots,
                  const std::vector<Vector>& simple_coroots) {
  require(static_cast<int>(metric_roots.size()) == data.positive_roots,
          data.name + " positive-root count mismatch");
  require(metric_roots.size() == coordinate_roots.size(),
          data.name + " coordinate-root count mismatch");

  mpq_class norm_sum = 0;
  for (const auto& root : metric_roots) norm_sum += inner(root, root);
  require(norm_sum == data.rank * data.kappa,
          data.name + " quadratic-frame constant mismatch");
  verify_quartic_identity(data, coordinate_roots, variables);

  const mpq_class scale = discriminant_scale(metric_roots);
  require(scale == data.expected_discriminant_scale,
          data.name + " short-root discriminant scale mismatch");

  const mpq_class coroot_covolume_squared = gram_determinant(simple_coroots);
  require(coroot_covolume_squared == data.expected_coroot_covolume_squared,
          data.name + " coroot-lattice covolume mismatch");
  const mpz_class center_index = cartan_determinant(simple_roots, simple_coroots);
  require(center_index == data.expected_center_index,
          data.name + " center-index mismatch");
  const mpq_class weyl_lattice_factor =
      mpq_class(center_index * center_index) / coroot_covolume_squared;
  require(weyl_lattice_factor == data.expected_weyl_lattice_factor,
          data.name + " normalized Weyl lattice factor mismatch");
  const mpq_class minimum_coweight_norm =
      minimum_dual_lattice_norm_squared(simple_roots);
  require(minimum_coweight_norm == 2,
          data.name + " adjoint-torus injectivity minimum mismatch");

  const int n = data.tail_start;
  const int alpha = 2 + 2 * data.positive_roots + data.rank;
  const int shape = data.dimension / 2;
  require(data.dimension % 2 == 0 && alpha == data.dimension + 2,
          data.name + " radial exponent mismatch");
  require(shape > 0, data.name + " invalid radial shape");

  const mpq_class pi_lower(333, 106), pi_upper(355, 113);
  const mpq_class cap = pi_lower * pi_lower * data.kappa * n /
                        (4 * data.dimension);
  require(data.rectangle.s1 < cap && data.rectangle.t1 < cap,
          data.name + " rectangle leaves the root-angle cap");
  require(data.rectangle.s0 >= data.rectangle.t1,
          data.name + " rectangle overlaps its transpose");

  const mpq_class gap =
      data.rectangle.s0 - data.rectangle.t1 -
      data.dimension * data.rectangle.s1 * data.rectangle.s1 /
          (6 * data.quartic_denominator * n);
  require(gap == data.expected_gap && gap > 0,
          data.name + " character-gap lower bound mismatch");

  const mpq_class a_coefficient = data.rectangle.s1 + data.rectangle.t1;
  const mpq_class b_coefficient =
      mpq_class(1, 18) * mpq_class(2 * data.dimension) /
      data.quartic_denominator *
      (data.rectangle.s0 * data.rectangle.s0 +
       data.rectangle.t0 * data.rectangle.t0);
  const mpq_class h = 1 - a_coefficient / n + b_coefficient / (n * n);
  require(h == data.expected_h && h > 0,
          data.name + " character-average lower bound mismatch");
  require(n * a_coefficient > 2 * b_coefficient,
          data.name + " h(n) monotonicity check failed");

  const mpq_class direct = qpow(mpq_class(2 * data.dimension) * h, 2) - 4;
  require(direct > 0, data.name + " direct Chain factor is not positive");

  const mpq_class a0_lower =
      a0_prefactor(data, scale * scale, weyl_lattice_factor) /
      qpow(pi_upper, data.rank);
  const mpq_class mass_lower = rectangle_mass_lower(data, gap);
  const mpq_class sine_exponent =
      mpq_class(4 * data.dimension, 18 * n) * a_coefficient;
  const mpq_class sine_lower = qpow(mpq_class(1, 3), ceiling(sine_exponent));
  const mpq_class positive_lower =
      a0_lower * qpow(mpq_class(2 * data.dimension), n) /
      qpow(mpq_class(n), alpha) * mass_lower * qpow(h, n) * direct * sine_lower;
  const mpz_class negative_factor =
      2 * ((2 * data.negative_constant) * (2 * data.negative_constant) + 4);
  const mpq_class negative =
      negative_factor * qpow(mpq_class(2 * data.negative_constant), n);
  const mpq_class base_ratio = positive_lower / negative;
  require(base_ratio > data.base_threshold,
          data.name + " base rectangular domination failed");

  const mpq_class step_ratio =
      qpow(mpq_class(data.dimension, data.negative_constant) * h, 2) *
      qpow(mpq_class(n, n + 2), alpha);
  require(step_ratio > data.step_threshold,
          data.name + " odd-step propagation ratio failed");

  std::cout << "group=" << data.name << " n0=" << n
            << " roots=" << data.positive_roots << " kappa=" << data.kappa
            << " quartic_denominator=" << data.quartic_denominator << '\n';
  std::cout << "rectangle_s=[" << data.rectangle.s0 << ',' << data.rectangle.s1
            << "] rectangle_t=[" << data.rectangle.t0 << ',' << data.rectangle.t1
            << "] gap=" << gap << " h=" << h << '\n';
  std::cout << "discriminant_scale=" << scale
            << " coroot_covolume_squared=" << coroot_covolume_squared
            << " center_index=" << center_index
            << " weyl_lattice_factor=" << weyl_lattice_factor
            << " minimum_coweight_norm_squared=" << minimum_coweight_norm
            << " radial_exp_ceiling=" << ceiling(a_coefficient)
            << " sine_exp_ceiling=" << ceiling(sine_exponent) << '\n';
  std::cout << "base_ratio_gt_" << data.base_threshold << ": OK\n";
  std::cout << "odd_step_ratio_gt_" << data.step_threshold << ": OK\n";
  std::cout << "quartic_identity_exact: OK\n"
            << "root_angle_cap_exact: OK\n"
            << "h_monotone_from_n0: OK\n";
}

}  // namespace

int main() {
  try {
    const GroupData g2{
        "G2", 2, 6, 14, 2, 4, {2, 6}, mpq_class(16, 5), 21,
        {mpq_class(19, 5), mpq_class(21, 5), mpq_class(7, 5), mpq_class(9, 5)},
        mpq_class(1, 27), mpq_class(3), mpz_class(1), mpq_class(1, 3),
        mpq_class(111, 80), mpq_class(1661, 2268), 1, 6};
    const GroupData f4{
        "F4", 4, 24, 52, 4, 9, {2, 6, 8, 12}, mpq_class(54, 5), 65,
        {mpq_class(51, 4), mpq_class(53, 4), mpq_class(31, 4), mpq_class(33, 4)},
        mpq_class(1, 4096), mpq_class(4), mpz_class(1), mpq_class(1, 4),
        mpq_class(3023, 1296), mpq_class(44063, 63180),
        mpz_class("10000000000000"), 10};

    verify_group(g2, g2_metric_roots(), g2_plane_forms(), 2,
                 g2_simple_roots(), g2_simple_coroots());
    const auto roots_f4 = f4_roots();
    verify_group(f4, roots_f4, roots_f4, 4,
                 f4_simple_roots(), f4_simple_coroots());
    std::cout << "RESULT: G2/F4 RECTANGULAR TAILS PASS (exact GMP rationals)\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "G2/F4 rectangular certificate failure: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
