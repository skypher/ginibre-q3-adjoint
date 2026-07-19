#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Rational = mpq_class;
using Polynomial = std::vector<Rational>;

namespace {

Rational q(long long numerator, long long denominator = 1) {
    Rational out(mpz_class(std::to_string(numerator)), mpz_class(std::to_string(denominator)));
    out.canonicalize();
    return out;
}

Rational pow_q(Rational base, unsigned int exponent) {
    Rational result(1);
    while (exponent > 0) {
        if (exponent & 1U) result *= base;
        exponent >>= 1U;
        if (exponent) base *= base;
    }
    result.canonicalize();
    return result;
}

mpz_class pow_z(unsigned int base, unsigned int exponent) {
    mpz_class result;
    mpz_ui_pow_ui(result.get_mpz_t(), base, exponent);
    return result;
}

mpz_class factorial_z(unsigned int n) {
    mpz_class result = 1;
    for (unsigned int k = 2; k <= n; ++k) result *= k;
    return result;
}

mpz_class binomial_z(unsigned int n, unsigned int k) {
    mpz_class result;
    mpz_bin_uiui(result.get_mpz_t(), n, k);
    return result;
}

double log10_z(const mpz_class& input) {
    mpz_class value = input >= 0 ? input : -input;
    if (value == 0) return -INFINITY;
    std::string digits = value.get_str();
    const std::size_t take = std::min<std::size_t>(16, digits.size());
    const double leading = std::stod(digits.substr(0, take));
    return std::log10(leading) + static_cast<double>(digits.size() - take);
}

double log10_q(const Rational& value) {
    if (value <= 0) return -INFINITY;
    return log10_z(value.get_num()) - log10_z(value.get_den());
}

std::string fixed3(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(3) << value;
    return out.str();
}

void trim(Polynomial& poly) {
    while (!poly.empty() && poly.back() == 0) poly.pop_back();
}

Polynomial poly_add(const Polynomial& left, const Polynomial& right) {
    Polynomial out(std::max(left.size(), right.size()), q(0));
    for (std::size_t i = 0; i < left.size(); ++i) out[i] += left[i];
    for (std::size_t i = 0; i < right.size(); ++i) out[i] += right[i];
    trim(out);
    return out;
}

Polynomial poly_scale(Polynomial poly, const Rational& scalar) {
    for (Rational& coefficient : poly) coefficient *= scalar;
    trim(poly);
    return poly;
}

Polynomial poly_mul(const Polynomial& left, const Polynomial& right) {
    if (left.empty() || right.empty()) return {};
    Polynomial out(left.size() + right.size() - 1, q(0));
    for (std::size_t i = 0; i < left.size(); ++i)
        for (std::size_t j = 0; j < right.size(); ++j)
            out[i + j] += left[i] * right[j];
    trim(out);
    return out;
}

Polynomial poly_pow(Polynomial base, unsigned int exponent) {
    Polynomial result{q(1)};
    while (exponent > 0) {
        if (exponent & 1U) result = poly_mul(result, base);
        exponent >>= 1U;
        if (exponent) base = poly_mul(base, base);
    }
    return result;
}

Rational poly_eval(const Polynomial& poly, const Rational& value) {
    Rational result = 0;
    for (auto it = poly.rbegin(); it != poly.rend(); ++it) result = result * value + *it;
    return result;
}

Polynomial compose_linear(const Polynomial& poly, const Rational& scale, const Rational& shift) {
    Polynomial result;
    Polynomial power{q(1)};
    const Polynomial linear{shift, scale};
    for (const Rational& coefficient : poly) {
        result = poly_add(result, poly_scale(power, coefficient));
        power = poly_mul(power, linear);
    }
    return result;
}

Polynomial chebyshev_t(unsigned int degree) {
    if (degree == 0) return {q(1)};
    if (degree == 1) return {q(0), q(1)};
    Polynomial previous{q(1)};
    Polynomial current{q(0), q(1)};
    const Polynomial variable{q(0), q(1)};
    for (unsigned int k = 2; k <= degree; ++k) {
        Polynomial next = poly_add(poly_scale(poly_mul(variable, current), q(2)),
                                   poly_scale(previous, q(-1)));
        previous = std::move(current);
        current = std::move(next);
    }
    return current;
}

std::pair<Polynomial, Polynomial> poly_divmod(Polynomial numerator, Polynomial denominator) {
    trim(numerator);
    trim(denominator);
    if (denominator.empty()) throw std::runtime_error("polynomial division by zero");
    if (numerator.size() < denominator.size()) return {{}, numerator};
    Polynomial quotient(numerator.size() - denominator.size() + 1, q(0));
    while (!numerator.empty() && numerator.size() >= denominator.size()) {
        const std::size_t degree = numerator.size() - denominator.size();
        const Rational coefficient = numerator.back() / denominator.back();
        quotient[degree] += coefficient;
        for (std::size_t i = 0; i < denominator.size(); ++i)
            numerator[i + degree] -= coefficient * denominator[i];
        trim(numerator);
    }
    trim(quotient);
    return {quotient, numerator};
}

std::vector<Rational> bernstein_coefficients_on_unit_interval(const Polynomial& poly) {
    if (poly.empty()) return {};
    const unsigned int degree = static_cast<unsigned int>(poly.size() - 1);
    std::vector<Rational> out(degree + 1, q(0));
    for (unsigned int i = 0; i <= degree; ++i) {
        for (unsigned int k = 0; k <= i; ++k) {
            if (k < poly.size())
                out[i] += poly[k] * Rational(binomial_z(i, k), binomial_z(degree, k));
        }
    }
    return out;
}

bool bernstein_positive_on(const Polynomial& poly, const Rational& left, const Rational& right) {
    const Polynomial unit = compose_linear(poly, right - left, left);
    const std::vector<Rational> coefficients = bernstein_coefficients_on_unit_interval(unit);
    return !coefficients.empty()
        && std::all_of(coefficients.begin(), coefficients.end(), [](const Rational& value) {
               return value > 0;
           });
}

bool bernstein_positive_subdivided(const Polynomial& poly, const Rational& left,
                                   const Rational& right, int parts) {
    for (int index = 0; index < parts; ++index) {
        const Rational a = left + (right - left) * index / parts;
        const Rational b = left + (right - left) * (index + 1) / parts;
        if (!bernstein_positive_on(poly, a, b)) return false;
    }
    return true;
}

int floor_q(const Rational& value) {
    mpz_class quotient;
    mpz_fdiv_q(quotient.get_mpz_t(), value.get_num_mpz_t(), value.get_den_mpz_t());
    return quotient.get_si();
}

Rational exp_fractional_septic_lower(const Rational& value) {
    const Rational e_upper = q(1457, 536);
    const int whole = floor_q(value);
    const Rational frac = value - whole;
    Rational poly = q(1) - frac + frac * frac / 2 - pow_q(frac, 3) / 6
        + pow_q(frac, 4) / 24 - pow_q(frac, 5) / 120
        + pow_q(frac, 6) / 720 - pow_q(frac, 7) / 5040;
    Rational out = poly / pow_q(e_upper, static_cast<unsigned int>(whole));
    out.canonicalize();
    return out;
}

struct Config {
    std::string group;
    int d;
    int rank;
    int r_plus;
    int kappa;
    int quartic_denom;
    int n;
    int alpha;
    int c_negative;
    int weyl_lattice_factor;
    Rational shape;
    Rational grid;
    Rational s_start;
    Rational s_end;
    Rational t_start;
    Rational t_end;
    Rational step_floor;
    int expected_cells;
    std::vector<int> degrees;
};

const Rational DELTA = q(14, 3);
const Rational DELTA_SQUARED = DELTA * DELTA;
const Rational GAP_LINEAR = q(99119225, 100000000);
const Rational GAP_QUARTIC = q(-7802858, 100000000);
const Rational GAP_SEXTIC = q(169078, 100000000);
const Rational COS_QUARTIC = q(83, 2000);

Config e6_config() {
    return {"E6", 78, 6, 36, 12, 16, 75, 80, 6, 3, q(39), q(1, 10),
            q(10), q(30), q(5), q(20), q(1), 23088, {2, 5, 6, 8, 9, 12}};
}

Config e7_config() {
    return {"E7", 133, 7, 63, 18, 27, 65, 135, 7, 2, q(133, 2), q(1, 20),
            q(15), q(35), q(8), q(28), q(101, 100), 66925,
            {2, 6, 8, 10, 12, 14, 18}};
}

struct Cell {
    Rational s0;
    Rational s1;
    Rational t0;
    Rational t1;
};

Rational gap_p_coefficient(const Config& cfg, const Cell& cell) {
    return GAP_QUARTIC * q(2 * cfg.d, cfg.quartic_denom) * cell.s1 * cell.s1
        + COS_QUARTIC * q(2 * cfg.d, cfg.quartic_denom) * cell.t0 * cell.t0;
}

Rational gap_q_coefficient(const Config& cfg, const Cell& cell) {
    return GAP_SEXTIC * q((2 * cfg.d) * (2 * cfg.d),
                          cfg.quartic_denom * cfg.quartic_denom)
        * pow_q(cell.s0, 3);
}

Rational gap_lower(const Config& cfg, const Cell& cell) {
    return GAP_LINEAR * cell.s0 - cell.t1
        + gap_p_coefficient(cfg, cell) / cfg.n
        + gap_q_coefficient(cfg, cell) / (cfg.n * cfg.n);
}

Rational h_b_coefficient(const Config& cfg, const Cell& cell) {
    return COS_QUARTIC * q(2 * cfg.d, cfg.quartic_denom)
        * (cell.s0 * cell.s0 + cell.t0 * cell.t0);
}

Rational h_lower(const Config& cfg, const Cell& cell) {
    return q(1) - (cell.s1 + cell.t1) / cfg.n
        + h_b_coefficient(cfg, cell) / (cfg.n * cfg.n);
}

Rational direct_lower(const Config& cfg, const Rational& h) {
    return q((2 * cfg.d) * (2 * cfg.d)) * h * h - 4;
}

Rational odd_step_lower(const Config& cfg, const Rational& h) {
    return pow_q(q(cfg.d, cfg.c_negative) * h, 2)
        * pow_q(q(cfg.n, cfg.n + 2), cfg.alpha);
}

Rational sine_exponent(const Config& cfg, const Cell& cell) {
    const Rational quadratic = q(2, 24) * q(2 * cfg.d)
        * (cell.s1 + cell.t1) / cfg.n;
    const Rational quartic_coefficient = q(1, 2880) + DELTA_SQUARED / 103800;
    const Rational quartic = 2 * quartic_coefficient
        * q((2 * cfg.d) * (2 * cfg.d), cfg.quartic_denom)
        * (cell.s1 * cell.s1 + cell.t1 * cell.t1) / (cfg.n * cfg.n);
    return quadratic + quartic;
}

Rational gamma_inv_square_lower(const Config& cfg) {
    if (cfg.group == "E6") {
        const mpz_class value = factorial_z(38);
        return Rational(1, value * value);
    }
    const mpz_class numerator = pow_z(4, 132) * factorial_z(66) * factorial_z(66);
    const mpz_class denominator = factorial_z(132) * factorial_z(132);
    return Rational(numerator, denominator) / q(355, 113);
}

Rational a0_lower(const Config& cfg) {
    mpz_class product = 1;
    for (int degree : cfg.degrees)
        product *= factorial_z(static_cast<unsigned int>(degree - 1));
    const Rational a = q(cfg.d, 2);
    Rational prefactor = q(8) * a * cfg.d * cfg.d
        * pow_q(q(cfg.d, cfg.kappa), cfg.d) * Rational(product * product)
        * cfg.weyl_lattice_factor
        / pow_z(2, cfg.rank);
    return prefactor / pow_q(q(355, 113), cfg.rank);
}

Rational half_integer_radial_integral(const Rational& low, const Rational& high) {
    // On [low,high], s^(131/2) >= low^65*sqrt(low).  Bound the last square
    // root downward by an exact six-decimal rational; no floating arithmetic
    // enters the certificate.
    const mpz_class scale = 1000000;
    const mpz_class scaled_numerator = low.get_num() * scale * scale;
    mpz_class radicand;
    mpz_fdiv_q(radicand.get_mpz_t(),
               scaled_numerator.get_mpz_t(),
               low.get_den_mpz_t());
    mpz_class root;
    mpz_sqrt(root.get_mpz_t(), radicand.get_mpz_t());
    const Rational sqrt_lower(root, scale);
    const Rational width = high - low;
    return width * pow_q(low, 65) * sqrt_lower;
}

Rational rectangle_fraction(const Config& cfg, const Cell& cell, const Rational& gap) {
    Rational s_integral;
    Rational t_integral;
    if (cfg.group == "E6") {
        s_integral = (pow_q(cell.s1, 39) - pow_q(cell.s0, 39)) / 39;
        t_integral = (pow_q(cell.t1, 39) - pow_q(cell.t0, 39)) / 39;
    } else {
        s_integral = half_integer_radial_integral(cell.s0, cell.s1);
        t_integral = half_integer_radial_integral(cell.t0, cell.t1);
    }
    Rational result = gap * gap * exp_fractional_septic_lower(cell.s1 + cell.t1)
        * s_integral * t_integral * gamma_inv_square_lower(cfg) / cfg.shape;
    result.canonicalize();
    return result;
}

bool cell_is_retained(const Config& cfg, const Cell& cell) {
    const Rational cap = DELTA_SQUARED * cfg.kappa * cfg.n / (4 * cfg.d);
    if (!(cell.s1 < cap && cell.t1 < cap)) return false;
    const Rational gap = gap_lower(cfg, cell);
    if (gap <= 0) return false;
    const Rational p = gap_p_coefficient(cfg, cell);
    const Rational qq = gap_q_coefficient(cfg, cell);
    if (!(-p * cfg.n > 2 * qq)) return false;
    const Rational h = h_lower(cfg, cell);
    if (h <= 0 || direct_lower(cfg, h) <= 0) return false;
    const Rational a = cell.s1 + cell.t1;
    const Rational b = h_b_coefficient(cfg, cell);
    if (!(cfg.n * a > 2 * b)) return false;
    return odd_step_lower(cfg, h) > cfg.step_floor;
}

std::vector<Cell> retained_cells(const Config& cfg) {
    std::vector<Cell> cells;
    for (Rational s0 = cfg.s_start; s0 + cfg.grid <= cfg.s_end; s0 += cfg.grid) {
        for (Rational t0 = cfg.t_start; t0 + cfg.grid <= cfg.t_end; t0 += cfg.grid) {
            Cell cell{s0, s0 + cfg.grid, t0, t0 + cfg.grid};
            if (cell_is_retained(cfg, cell)) cells.push_back(std::move(cell));
        }
    }
    return cells;
}

struct Stats {
    Rational positive = 0;
    Rational rect_total = 0;
    Rational min_gap;
    Rational min_h;
    Rational min_direct;
    Rational min_step;
    Rational max_sine = 0;
    bool initialized = false;

    void add(const Rational& term, const Rational& rect, const Rational& gap,
             const Rational& h, const Rational& direct, const Rational& step,
             const Rational& sine) {
        positive += term;
        rect_total += rect;
        if (!initialized) {
            min_gap = gap;
            min_h = h;
            min_direct = direct;
            min_step = step;
            max_sine = sine;
            initialized = true;
        } else {
            if (gap < min_gap) min_gap = gap;
            if (h < min_h) min_h = h;
            if (direct < min_direct) min_direct = direct;
            if (step < min_step) min_step = step;
            if (sine > max_sine) max_sine = sine;
        }
    }

    void merge(const Stats& other) {
        if (!other.initialized) return;
        if (!initialized) {
            *this = other;
            return;
        }
        positive += other.positive;
        rect_total += other.rect_total;
        if (other.min_gap < min_gap) min_gap = other.min_gap;
        if (other.min_h < min_h) min_h = other.min_h;
        if (other.min_direct < min_direct) min_direct = other.min_direct;
        if (other.min_step < min_step) min_step = other.min_step;
        if (other.max_sine > max_sine) max_sine = other.max_sine;
    }
};

Stats compute_stats(const Config& cfg, const std::vector<Cell>& cells) {
    const Rational common = a0_lower(cfg) * pow_q(q(2 * cfg.d), cfg.n)
        / pow_q(q(cfg.n), cfg.alpha);
    const int threads = omp_get_max_threads();
    std::vector<Stats> partial(static_cast<std::size_t>(threads));
    const auto started = std::chrono::steady_clock::now();

#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
#pragma omp for schedule(dynamic, 32)
        for (std::size_t index = 0; index < cells.size(); ++index) {
            const Cell& cell = cells[index];
            const Rational gap = gap_lower(cfg, cell);
            const Rational h = h_lower(cfg, cell);
            const Rational direct = direct_lower(cfg, h);
            const Rational rect = rectangle_fraction(cfg, cell, gap);
            const Rational sine = sine_exponent(cfg, cell);
            Rational term = common * rect * pow_q(h, cfg.n) * direct
                * exp_fractional_septic_lower(sine);
            term.canonicalize();
            partial[static_cast<std::size_t>(tid)].add(
                term, rect, gap, h, direct, odd_step_lower(cfg, h), sine);

            if ((index + 1) % 20000 == 0 && tid == 0) {
                const double elapsed = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - started).count();
                std::cout << "progress group=" << cfg.group << " index~" << (index + 1)
                          << "/" << cells.size() << " elapsed=" << elapsed << "s\n";
            }
        }
    }
    Stats total;
    for (const Stats& item : partial) total.merge(item);
    return total;
}

std::vector<mpz_class> read_moments(const std::string& path, int maximum) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("could not open moment log: " + path);
    std::vector<mpz_class> moments(static_cast<std::size_t>(maximum + 1));
    std::vector<bool> seen(static_cast<std::size_t>(maximum + 1), false);
    const std::regex pattern(R"(^m_([0-9]+) = ([0-9]+))");
    std::string line;
    while (std::getline(input, line)) {
        std::smatch match;
        if (!std::regex_search(line, match, pattern)) continue;
        const int index = std::stoi(match[1].str());
        if (index <= maximum) {
            moments[static_cast<std::size_t>(index)] = mpz_class(match[2].str());
            seen[static_cast<std::size_t>(index)] = true;
        }
    }
    if (!std::all_of(seen.begin(), seen.end(), [](bool value) { return value; }))
        throw std::runtime_error("moment log is missing a required prefix entry");
    return moments;
}

mpz_class sum_moment(const std::vector<mpz_class>& moments, int power) {
    mpz_class total = 0;
    for (int k = 0; k <= power; ++k) {
        const mpz_class choose = binomial_z(static_cast<unsigned int>(power),
                                            static_cast<unsigned int>(k));
        total += choose * moments[k + 2] * moments[power - k];
        total += choose * moments[k] * moments[power - k + 2];
        total -= 2 * choose * moments[k + 1] * moments[power - k + 1];
    }
    return total;
}

Polynomial e7_chebyshev_factor() {
    const Polynomial raw = compose_linear(chebyshev_t(26), q(2, 259), q(-273, 259));
    const Rational at_minus_fourteen = poly_eval(raw, q(-14));
    return poly_scale(raw, q(1) / at_minus_fourteen);
}

Rational e7_negative_integral(const std::vector<mpz_class>& moments) {
    const Polynomial square = poly_mul(e7_chebyshev_factor(), e7_chebyshev_factor());
    Rational total = 0;
    for (std::size_t i = 0; i < square.size(); ++i) {
        total += square[i] * (sum_moment(moments, 16 + static_cast<int>(i))
                              + 4 * sum_moment(moments, 14 + static_cast<int>(i)));
    }
    total.canonicalize();
    return total;
}

bool e7_negative_interval_check() {
    constexpr int n = 65;
    constexpr int p = 14;
    const Polynomial factor = e7_chebyshev_factor();
    const Polynomial x{q(0), q(1)};
    const Polynomial x2 = poly_mul(x, x);
    const Polynomial xp = poly_pow(x, p);
    const Rational scale = q(24, 25) * pow_z(14, n - p);

    const Polynomial factor_negative = compose_linear(factor, q(-1), q(0));
    const Polynomial lhs_negative = poly_scale(
        poly_mul(poly_mul(xp, poly_add(x2, {q(4)})),
                 poly_mul(factor_negative, factor_negative)), scale);
    const Polynomial rhs_negative = poly_mul(poly_pow(x, n), poly_add(x2, {q(-4)}));
    const Polynomial residual_negative = poly_add(lhs_negative, poly_scale(rhs_negative, q(-1)));
    auto divided_x = poly_divmod(residual_negative, xp);
    if (!divided_x.second.empty()) return false;
    auto divided_endpoint = poly_divmod(divided_x.first, {q(14), q(-1)});
    if (!divided_endpoint.second.empty()) return false;

    const Polynomial lhs_positive = poly_scale(
        poly_mul(poly_mul(xp, poly_add(x2, {q(4)})), poly_mul(factor, factor)), scale);
    const Polynomial rhs_positive = poly_mul(
        poly_pow(x, n), poly_add({q(4)}, poly_scale(x2, q(-1))));
    const Polynomial residual_positive = poly_add(lhs_positive, poly_scale(rhs_positive, q(-1)));
    auto positive_divided = poly_divmod(residual_positive, xp);
    if (!positive_divided.second.empty()) return false;

    return bernstein_positive_on(divided_endpoint.first, q(2), q(14))
        && bernstein_positive_on(positive_divided.first, q(0), q(2));
}

bool scalar_polynomial_checks() {
    const Rational endpoint = DELTA_SQUARED;

    // The alternating Taylor lower bound for 2(1-cos(sqrt(w))) ends at -w^6/12!.
    Polynomial gap_residual{
        q(0), q(1) - GAP_LINEAR, q(-1, 12) - GAP_QUARTIC,
        q(1, 360) - GAP_SEXTIC, q(-1, 20160), q(1, 1814400),
        q(-1, 239500800)};
    gap_residual.erase(gap_residual.begin());

    // The alternating Taylor lower bound for 2cos(sqrt(w)) is taken through -2w^9/18!.
    Polynomial cosine_residual;
    for (int k = 2; k <= 9; ++k) {
        Rational coefficient(2 * ((k % 2 == 0) ? 1 : -1), factorial_z(2 * k));
        if (k == 2) coefficient -= COS_QUARTIC;
        cosine_residual.push_back(coefficient);
    }

    // sinc(sqrt(y))*exp(y/6+y^2/180+64y^3/103800)-1 >= 0.
    const Polynomial exponent{q(0), q(1, 6), q(1, 180), q(64, 103800)};
    Polynomial exp_lower;
    for (int k = 0; k <= 7; ++k)
        exp_lower = poly_add(exp_lower, poly_scale(poly_pow(exponent, k),
                                                   Rational(1, factorial_z(k))));
    Polynomial sinc_lower;
    for (int k = 0; k <= 7; ++k)
        sinc_lower.push_back(Rational((k % 2 == 0) ? 1 : -1, factorial_z(2 * k + 1)));
    Polynomial sine_residual = poly_add(poly_mul(sinc_lower, exp_lower), {q(-1)});
    while (!sine_residual.empty() && sine_residual.front() == 0)
        sine_residual.erase(sine_residual.begin());

    return bernstein_positive_subdivided(gap_residual, q(0), endpoint, 4)
        && bernstein_positive_subdivided(cosine_residual, q(0), endpoint, 4)
        && bernstein_positive_on(sine_residual, q(0), endpoint / 4);
}

using Root = std::vector<int>;
using IntegerPolynomial = std::map<std::vector<int>, mpz_class>;

std::vector<std::vector<int>> simply_laced_cartan(int rank) {
    std::vector<std::vector<int>> matrix(rank, std::vector<int>(rank, 0));
    for (int i = 0; i < rank; ++i) matrix[i][i] = 2;
    std::vector<std::pair<int, int>> edges;
    if (rank == 6) {
        edges = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {2, 5}};
    } else if (rank == 7) {
        edges = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {2, 6}};
    } else {
        throw std::runtime_error("unsupported exceptional rank");
    }
    for (const auto& edge : edges) {
        matrix[edge.first][edge.second] = -1;
        matrix[edge.second][edge.first] = -1;
    }
    return matrix;
}

mpz_class integer_determinant(std::vector<std::vector<mpz_class>> matrix) {
    if (matrix.empty()) throw std::runtime_error("empty determinant matrix");
    const std::size_t size = matrix.size();
    for (const auto& row : matrix)
        if (row.size() != size) throw std::runtime_error("nonsquare determinant matrix");
    mpz_class sign = 1;
    mpz_class previous = 1;
    for (std::size_t column = 0; column + 1 < size; ++column) {
        std::size_t pivot = column;
        while (pivot < size && matrix[pivot][column] == 0) ++pivot;
        if (pivot == size) return 0;
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            sign = -sign;
        }
        const mpz_class pivot_value = matrix[column][column];
        for (std::size_t row = column + 1; row < size; ++row)
            for (std::size_t entry = column + 1; entry < size; ++entry) {
                mpz_class numerator = matrix[row][entry] * pivot_value
                    - matrix[row][column] * matrix[column][entry];
                if (column != 0) {
                    if (numerator % previous != 0)
                        throw std::runtime_error("Bareiss division was not exact");
                    numerator /= previous;
                }
                matrix[row][entry] = numerator;
            }
        previous = pivot_value;
    }
    return sign * matrix.back().back();
}

bool weyl_lattice_factor_check(int rank, int expected_factor) {
    const auto cartan = simply_laced_cartan(rank);
    std::vector<std::vector<mpz_class>> exact(
        cartan.size(), std::vector<mpz_class>(cartan.size()));
    for (std::size_t i = 0; i < cartan.size(); ++i)
        for (std::size_t j = 0; j < cartan.size(); ++j)
            exact[i][j] = cartan[i][j];
    mpz_class determinant = integer_determinant(std::move(exact));
    if (determinant < 0) determinant = -determinant;
    // In the simply-laced adjoint form, |P^vee/Q^vee|=det(A) and
    // covol(Q^vee)^2=det(A), so the normalized pair factor is det(A).
    return determinant == expected_factor;
}

std::vector<std::vector<Rational>> inverse_matrix(
    const std::vector<std::vector<int>>& input) {
    const std::size_t size = input.size();
    std::vector<std::vector<Rational>> augmented(
        size, std::vector<Rational>(2 * size, q(0)));
    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = 0; j < size; ++j) augmented[i][j] = q(input[i][j]);
        augmented[i][size + i] = q(1);
    }
    for (std::size_t column = 0; column < size; ++column) {
        std::size_t pivot = column;
        while (pivot < size && augmented[pivot][column] == 0) ++pivot;
        if (pivot == size) throw std::runtime_error("singular rational matrix");
        if (pivot != column) std::swap(augmented[pivot], augmented[column]);
        const Rational pivot_value = augmented[column][column];
        for (Rational& entry : augmented[column]) entry /= pivot_value;
        for (std::size_t row = 0; row < size; ++row) {
            if (row == column) continue;
            const Rational multiplier = augmented[row][column];
            for (std::size_t entry = 0; entry < 2 * size; ++entry)
                augmented[row][entry] -= multiplier * augmented[column][entry];
        }
    }
    std::vector<std::vector<Rational>> inverse(
        size, std::vector<Rational>(size, q(0)));
    for (std::size_t i = 0; i < size; ++i)
        for (std::size_t j = 0; j < size; ++j)
            inverse[i][j] = augmented[i][size + j];
    return inverse;
}

Rational coweight_norm_squared(const std::vector<int>& coordinates,
                               const std::vector<std::vector<Rational>>& gram) {
    Rational result = q(0);
    for (std::size_t i = 0; i < coordinates.size(); ++i)
        for (std::size_t j = 0; j < coordinates.size(); ++j)
            result += coordinates[i] * coordinates[j] * gram[i][j];
    result.canonicalize();
    return result;
}

Rational minimum_coweight_norm_squared(int rank) {
    const auto gram = inverse_matrix(simply_laced_cartan(rank));
    std::vector<int> coordinates(static_cast<std::size_t>(rank), -1);
    Rational minimum = q(2);
    const unsigned long cases = pow_z(3, static_cast<unsigned int>(rank)).get_ui();
    for (unsigned long code = 0; code < cases; ++code) {
        unsigned long cursor = code;
        bool nonzero = false;
        for (int i = 0; i < rank; ++i) {
            coordinates[static_cast<std::size_t>(i)] = static_cast<int>(cursor % 3) - 1;
            cursor /= 3;
            nonzero = nonzero || coordinates[static_cast<std::size_t>(i)] != 0;
        }
        if (!nonzero) continue;
        const Rational value = coweight_norm_squared(coordinates, gram);
        if (value < minimum) minimum = value;
    }
    // If some simple-root pairing has absolute value at least two, Cauchy
    // and |alpha_i|^2=2 give |lambda|^2>=2.  Thus the ternary enumeration
    // is global whenever it finds a value below two.
    if (!(minimum < 2)) throw std::runtime_error("coweight minimum not isolated");
    return minimum;
}

std::vector<Root> positive_roots(int rank) {
    const auto cartan = simply_laced_cartan(rank);
    std::set<Root> roots;
    std::vector<Root> queue;
    for (int i = 0; i < rank; ++i) {
        Root root(rank, 0);
        root[i] = 1;
        roots.insert(root);
        queue.push_back(root);
    }
    for (std::size_t cursor = 0; cursor < queue.size(); ++cursor) {
        const Root root = queue[cursor];
        for (int i = 0; i < rank; ++i) {
            int pairing = 0;
            for (int j = 0; j < rank; ++j) pairing += root[j] * cartan[j][i];
            Root reflected = root;
            reflected[i] -= pairing;
            if (roots.insert(reflected).second) queue.push_back(std::move(reflected));
        }
    }
    std::vector<Root> positive;
    for (const Root& root : roots) {
        if (std::all_of(root.begin(), root.end(), [](int value) { return value >= 0; }))
            positive.push_back(root);
    }
    return positive;
}

IntegerPolynomial linear_power_sum(const std::vector<Root>& roots, int power) {
    if (power != 2 && power != 4) throw std::runtime_error("unsupported root power");
    const int rank = static_cast<int>(roots.front().size());
    IntegerPolynomial result;
    for (const Root& root : roots) {
        if (power == 2) {
            for (int i = 0; i < rank; ++i)
                for (int j = 0; j < rank; ++j) {
                    std::vector<int> exponent(rank, 0);
                    ++exponent[i];
                    ++exponent[j];
                    result[exponent] += root[i] * root[j];
                }
        } else {
            for (int i = 0; i < rank; ++i)
                for (int j = 0; j < rank; ++j)
                    for (int k = 0; k < rank; ++k)
                        for (int ell = 0; ell < rank; ++ell) {
                            std::vector<int> exponent(rank, 0);
                            ++exponent[i];
                            ++exponent[j];
                            ++exponent[k];
                            ++exponent[ell];
                            result[exponent] += root[i] * root[j] * root[k] * root[ell];
                        }
        }
    }
    return result;
}

IntegerPolynomial integer_poly_square(const IntegerPolynomial& input) {
    IntegerPolynomial result;
    for (const auto& left : input)
        for (const auto& right : input) {
            std::vector<int> exponent(left.first.size(), 0);
            for (std::size_t i = 0; i < exponent.size(); ++i)
                exponent[i] = left.first[i] + right.first[i];
            result[exponent] += left.second * right.second;
        }
    return result;
}

bool root_quartic_identity_check(int rank, int expected_positive, int quartic_denom) {
    const std::vector<Root> roots = positive_roots(rank);
    if (static_cast<int>(roots.size()) != expected_positive) return false;
    const IntegerPolynomial quadratic = linear_power_sum(roots, 2);
    const IntegerPolynomial quartic = linear_power_sum(roots, 4);
    IntegerPolynomial residual = integer_poly_square(quadratic);
    for (const auto& term : quartic) residual[term.first] -= quartic_denom * term.second;
    return std::all_of(residual.begin(), residual.end(), [](const auto& term) {
        return term.second == 0;
    });
}

Rational negative_bound(const Config& cfg, const std::vector<mpz_class>& e7_moments) {
    if (cfg.group == "E6") {
        return q(2 * ((2 * cfg.c_negative) * (2 * cfg.c_negative) + 4))
            * pow_z(2 * cfg.c_negative, cfg.n);
    }
    return q(24, 25) * pow_z(14, cfg.n - 14) * e7_negative_integral(e7_moments);
}

bool run_case(const Config& cfg, const std::vector<mpz_class>& e7_moments) {
    const std::vector<Cell> cells = retained_cells(cfg);
    const bool cell_count_ok = static_cast<int>(cells.size()) == cfg.expected_cells;
    const bool transpose_disjoint = std::all_of(
        cells.begin(), cells.end(), [](const Cell& cell) { return cell.s0 >= cell.t1; });
    const Stats stats = compute_stats(cfg, cells);
    const Rational negative = negative_bound(cfg, e7_moments);
    const Rational ratio = stats.positive / negative;
    const Rational cap = DELTA_SQUARED * cfg.kappa * cfg.n / (4 * cfg.d);
    const bool pass = cell_count_ok && transpose_disjoint
        && stats.initialized && stats.min_gap > 0
        && stats.min_h > 0 && stats.min_direct > 0 && stats.min_step > 1 && ratio > 1;

    std::cout << "group=" << cfg.group << " n0=" << cfg.n << " delta=" << DELTA
              << " grid=" << cfg.grid << " cap=" << cap << "\n";
    std::cout << "retained_cells=" << cells.size() << " expected=" << cfg.expected_cells << "\n";
    std::cout << "gap_min=" << stats.min_gap << "\n";
    std::cout << "h_min=" << stats.min_h << " direct_min=" << stats.min_direct << "\n";
    std::cout << "sine_exponent_max=" << stats.max_sine << "\n";
    std::cout << "log10_positive_lower=" << fixed3(log10_q(stats.positive)) << "\n";
    std::cout << "log10_negative_upper=" << fixed3(log10_q(negative)) << "\n";
    std::cout << "log10_base_ratio_lower=" << fixed3(log10_q(ratio)) << "\n";
    std::cout << "log10_odd_step_ratio_lower=" << fixed3(log10_q(stats.min_step)) << "\n";
    std::cout << "cell_count: " << (cell_count_ok ? "OK" : "FAIL") << "\n";
    std::cout << "transposed_cells_disjoint: "
              << (transpose_disjoint ? "OK" : "FAIL") << "\n";
    std::cout << "base_ratio_gt_1: " << (ratio > 1 ? "OK" : "FAIL") << "\n";
    std::cout << "odd_step_ratio_gt_1: " << (stats.min_step > 1 ? "OK" : "FAIL") << "\n";
    std::cout << "case_result: " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

}  // namespace

int main(int argc, char** argv) {
    std::string e7_moment_path = "ginibre_q3/character_ring_iter/logs/e7_70.log";
    if (argc == 3 && std::string(argv[1]) == "--e7-moments") {
        e7_moment_path = argv[2];
    } else if (argc != 1) {
        throw std::runtime_error("usage: direct_chain_rect_e6_e7_rank_certificate [--e7-moments PATH]");
    }

    std::cout << "COMMAND: direct_chain_rect_e6_e7_rank_certificate --e7-moments "
              << e7_moment_path << "\n";
    std::cout << "OpenMP max threads=" << omp_get_max_threads() << "\n";
    const bool scalar_ok = scalar_polynomial_checks();
    std::cout << "scalar_gap_cosine_sine_bernstein: " << (scalar_ok ? "OK" : "FAIL") << "\n";
    const bool root_ok = root_quartic_identity_check(6, 36, 16)
        && root_quartic_identity_check(7, 63, 27);
    std::cout << "e6_e7_root_counts_and_quartic_identities: "
              << (root_ok ? "OK" : "FAIL") << "\n";
    const bool lattice_ok = weyl_lattice_factor_check(6, 3)
        && weyl_lattice_factor_check(7, 2)
        && minimum_coweight_norm_squared(6) == q(4, 3)
        && minimum_coweight_norm_squared(7) == q(3, 2);
    std::cout << "e6_e7_normalized_weyl_lattice_factors_3_2: "
              << (lattice_ok ? "OK" : "FAIL") << "\n";
    std::cout << "e6_e7_adjoint_torus_injectivity_minima_4/3_3/2: "
              << (lattice_ok ? "OK" : "FAIL") << "\n";

    const std::vector<mpz_class> moments = read_moments(e7_moment_path, 70);
    const bool moment_endpoints_ok = moments[0] == 1 && moments[1] == 0
        && moments[70] == mpz_class(
            "129078404848467795345644841293487698699547974717178497889513558681954344108627813155932870113361484");
    const bool interval_ok = e7_negative_interval_check();
    const Rational integral = e7_negative_integral(moments);
    const bool integral_positive = integral > 0;
    std::cout << "e7_moment_prefix_0_70: " << (moment_endpoints_ok ? "OK" : "FAIL") << "\n";
    std::cout << "e7_chebyshev_negative_intervals_bernstein: " << (interval_ok ? "OK" : "FAIL") << "\n";
    std::cout << "e7_negative_integral_positive: " << (integral_positive ? "OK" : "FAIL") << "\n";
    std::cout << "log10_e7_negative_integral=" << fixed3(log10_q(integral)) << "\n";

    const bool e6_ok = run_case(e6_config(), moments);
    const bool e7_ok = run_case(e7_config(), moments);
    const bool all_ok = scalar_ok && root_ok && lattice_ok
        && moment_endpoints_ok && interval_ok
        && integral_positive && e6_ok && e7_ok;
    std::cout << "RESULT: " << (all_ok ? "ALL PASS" : "FAIL") << "\n";
    return all_ok ? 0 : 1;
}
