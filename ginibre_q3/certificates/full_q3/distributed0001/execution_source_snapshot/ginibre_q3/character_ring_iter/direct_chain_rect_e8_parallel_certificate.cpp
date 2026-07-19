#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Rational = mpq_class;

namespace {

constexpr int D = 248;
constexpr int C_NEG = 8;
constexpr int RANK = 8;
constexpr int ALPHA = 250;
constexpr int SHAPE = 124;
constexpr int QUARTIC_DENOM = 50;
constexpr int SEXTIC_DENOM = 1800;
constexpr int EIGHTH_LOWER_DENOM = 300000;

Rational q(long long numerator, long long denominator = 1) {
    Rational out(mpz_class(std::to_string(numerator)), mpz_class(std::to_string(denominator)));
    out.canonicalize();
    return out;
}

Rational parse_q(const std::string& text) {
    Rational out(text);
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

int floor_q(const Rational& value) {
    mpz_class quotient;
    mpz_fdiv_q(quotient.get_mpz_t(), value.get_num_mpz_t(), value.get_den_mpz_t());
    return quotient.get_si();
}

int ceil_q_nonnegative(const Rational& value) {
    if (value < 0) throw std::runtime_error("ceil_q_nonnegative received a negative value");
    mpz_class quotient;
    mpz_cdiv_q(quotient.get_mpz_t(), value.get_num_mpz_t(), value.get_den_mpz_t());
    return quotient.get_si();
}

double log10_z(const mpz_class& input) {
    mpz_class value = input >= 0 ? input : -input;
    if (value == 0) return -INFINITY;
    std::string digits = value.get_str();
    std::size_t take = std::min<std::size_t>(16, digits.size());
    double leading = std::stod(digits.substr(0, take));
    return std::log10(leading) + static_cast<double>(digits.size() - take);
}

double log10_q(const Rational& value) {
    if (value <= 0) return -INFINITY;
    return log10_z(value.get_num()) - log10_z(value.get_den());
}

std::string fixed2(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

bool e8_weyl_lattice_normalization_check() {
    std::vector<std::vector<Rational>> cartan(8, std::vector<Rational>(8, q(0)));
    for (int i = 0; i < 8; ++i) cartan[i][i] = q(2);
    for (const auto& edge : std::vector<std::pair<int, int>>{
             {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {2, 7}}) {
        cartan[edge.first][edge.second] = q(-1);
        cartan[edge.second][edge.first] = q(-1);
    }
    Rational determinant = q(1);
    for (int column = 0; column < 8; ++column) {
        int pivot = column;
        while (pivot < 8 && cartan[pivot][column] == 0) ++pivot;
        if (pivot == 8) return false;
        if (pivot != column) {
            std::swap(cartan[pivot], cartan[column]);
            determinant = -determinant;
        }
        const Rational pivot_value = cartan[column][column];
        determinant *= pivot_value;
        for (int row = column + 1; row < 8; ++row) {
            const Rational multiplier = cartan[row][column] / pivot_value;
            for (int entry = column; entry < 8; ++entry)
                cartan[row][entry] -= multiplier * cartan[column][entry];
        }
    }
    // det(A)=1 makes P^vee=Q^vee and the normalized pair factor one.
    // The original integral Cartan form is even with diagonal two, so every
    // nonzero coweight has squared norm at least two; this also validates the
    // injectivity of every delta<=14/3 cap used below.
    return determinant == 1;
}

Rational exp_fractional_septic_lower(const Rational& value, const Rational& e_upper) {
    int whole = floor_q(value);
    Rational frac = value - whole;
    Rational poly = q(1) - frac + frac * frac / 2 - pow_q(frac, 3) / 6
        + pow_q(frac, 4) / 24 - pow_q(frac, 5) / 120
        + pow_q(frac, 6) / 720 - pow_q(frac, 7) / 5040;
    Rational out = poly / pow_q(e_upper, static_cast<unsigned int>(whole));
    out.canonicalize();
    return out;
}

Rational exp_fractional_cubic_lower(const Rational& value, const Rational& e_upper) {
    int whole = floor_q(value);
    Rational frac = value - whole;
    Rational poly = q(1) - frac + frac * frac / 2 - pow_q(frac, 3) / 6;
    Rational out = poly / pow_q(e_upper, static_cast<unsigned int>(whole));
    out.canonicalize();
    return out;
}

Rational exp_ceiling_lower(const Rational& value, const Rational& e_upper) {
    return q(1) / pow_q(e_upper, static_cast<unsigned int>(ceil_q_nonnegative(value)));
}

Rational gamma_inv_square() {
    mpz_class f = factorial_z(SHAPE - 1);
    return Rational(1, f * f);
}

Rational a0_lower() {
    const int degrees[] = {2, 8, 12, 14, 18, 20, 24, 30};
    mpz_class prod = 1;
    for (int degree : degrees) prod *= factorial_z(static_cast<unsigned int>(degree - 1));

    Rational prefactor = q(8) * q(D, 2) * D * D * pow_q(q(D, 30), D)
        * Rational(prod * prod) / pow_z(2, RANK);
    return prefactor / pow_q(q(355, 113), RANK);
}

enum class Formula {
    Cubic,
    CubicCap14,
    CubicN91,
    CubicWide93,
    QuarticMatched119,
    QuarticEighth121,
    QuarticSexticCubicExp,
    QuarticSexticCeiling,
    Tail133
};

struct CaseConfig {
    std::string label;
    int n;
    int expected_cells;
    Formula formula;
    Rational delta = q(4);
    Rational grid = q(1, 20);
    Rational s_start = q(40);
    Rational s_end = q(-1);
    Rational t_start = q(30);
    Rational t_end = q(48);
    bool s_end_is_cap = false;
    bool t_end_is_cap = true;
    Rational gap6_lower = q(0);
    Rational sine_denom = q(1);
};

struct Cell {
    Rational s0;
    Rational s1;
    Rational t0;
    Rational t1;
};

struct Stats {
    Rational positive_total = q(0);
    Rational rect_total = q(0);
    Rational best_term = q(0);
    Rational min_gap;
    Rational min_h;
    Rational min_direct;
    Rational max_sine_exponent = q(0);
    bool initialized = false;

    void add(
        const Rational& term,
        const Rational& rect,
        const Rational& gap,
        const Rational& h,
        const Rational& direct,
        const Rational& sine_exponent
    ) {
        positive_total += term;
        rect_total += rect;
        if (term > best_term) best_term = term;
        if (!initialized) {
            min_gap = gap;
            min_h = h;
            min_direct = direct;
            max_sine_exponent = sine_exponent;
            initialized = true;
        } else {
            if (gap < min_gap) min_gap = gap;
            if (h < min_h) min_h = h;
            if (direct < min_direct) min_direct = direct;
            if (sine_exponent > max_sine_exponent) max_sine_exponent = sine_exponent;
        }
    }

    void merge(const Stats& other) {
        if (!other.initialized) return;
        if (!initialized) {
            *this = other;
            return;
        }
        positive_total += other.positive_total;
        rect_total += other.rect_total;
        if (other.best_term > best_term) best_term = other.best_term;
        if (other.min_gap < min_gap) min_gap = other.min_gap;
        if (other.min_h < min_h) min_h = other.min_h;
        if (other.min_direct < min_direct) min_direct = other.min_direct;
        if (other.max_sine_exponent > max_sine_exponent) {
            max_sine_exponent = other.max_sine_exponent;
        }
    }
};

Rational cap_value(const CaseConfig& cfg) {
    if (cfg.formula == Formula::Tail133) {
        return q(333, 106) * q(333, 106) * 30 * cfg.n / (4 * D);
    }
    return cfg.delta * cfg.delta * 30 * cfg.n / (4 * D);
}

Rational gap_lower(const CaseConfig& cfg, const Cell& cell) {
    const Rational s0 = cell.s0;
    const Rational s1 = cell.s1;
    const Rational t0 = cell.t0;
    const Rational t1 = cell.t1;

    if (cfg.formula == Formula::QuarticMatched119) {
        Rational rho = q(D, 6 * QUARTIC_DENOM * cfg.n);
        Rational sixth_scale = q((2 * D) * (2 * D), SEXTIC_DENOM);
        return s0 - t1 - rho * s1 * s1 + rho * t0 * t0
            + cfg.gap6_lower * sixth_scale * pow_q(s0, 3) / (cfg.n * cfg.n)
            - q(1, 360) * sixth_scale * pow_q(t1, 3) / (cfg.n * cfg.n);
    }

    if (cfg.formula == Formula::QuarticEighth121
        || cfg.formula == Formula::QuarticSexticCubicExp
        || cfg.formula == Formula::QuarticSexticCeiling) {
        Rational rho = q(D, 6 * QUARTIC_DENOM * cfg.n);
        Rational sixth_gain = cfg.gap6_lower * q((2 * D) * (2 * D), SEXTIC_DENOM);
        return s0 - t1 - rho * s1 * s1 + sixth_gain * pow_q(s0, 3) / (cfg.n * cfg.n);
    }

    if (cfg.formula == Formula::Tail133) {
        Rational rho = q(D, 6 * QUARTIC_DENOM * cfg.n);
        return s0 - t1 - rho * s1 * s1;
    }

    Rational c0 = q(99119225, 100000000);
    Rational c1 = q(-7802858, 100000000);
    Rational c2 = q(181097, 100000000);
    if (cfg.formula == Formula::CubicCap14) c2 = q(169078, 100000000);
    if (cfg.formula == Formula::CubicN91) c2 = q(1562, 890661);
    if (cfg.formula == Formula::CubicWide93) c2 = q(7, 4000);

    Rational upper0 = q(1);
    Rational upper1 = q(-8228444, 100000000);
    Rational upper2 = q(233397, 100000000);
    Rational lower = c0 * s0
        + c1 * q(2 * D, QUARTIC_DENOM) * s0 * s0 / cfg.n
        + c2 * q((2 * D) * (2 * D), SEXTIC_DENOM) * pow_q(s0, 3) / (cfg.n * cfg.n);
    Rational upper = upper0 * t1
        + upper1 * q(2 * D, QUARTIC_DENOM) * t1 * t1 / cfg.n
        + upper2 * q((2 * D) * (2 * D), SEXTIC_DENOM) * pow_q(t1, 3) / (cfg.n * cfg.n);
    return lower - upper;
}

Rational local_h_lower(const CaseConfig& cfg, const Cell& cell) {
    const Rational s0 = cell.s0;
    const Rational s1 = cell.s1;
    const Rational t0 = cell.t0;
    const Rational t1 = cell.t1;

    if (cfg.formula == Formula::QuarticMatched119 || cfg.formula == Formula::QuarticEighth121) {
        Rational quartic_gain = q(1, 12) * q(2 * D, QUARTIC_DENOM);
        Rational sixth_loss = q(1, 360) * q((2 * D) * (2 * D), SEXTIC_DENOM);
        Rational eighth_gain = q(1, 30000) * q((2 * D) * (2 * D) * (2 * D), EIGHTH_LOWER_DENOM);
        return q(1) - (s1 + t1) / cfg.n
            + quartic_gain * (s0 * s0 + t0 * t0) / (cfg.n * cfg.n)
            - sixth_loss * (pow_q(s1, 3) + pow_q(t1, 3)) / (cfg.n * cfg.n * cfg.n)
            + eighth_gain * (pow_q(s0, 4) + pow_q(t0, 4))
                * q(1, static_cast<long long>(cfg.n) * cfg.n * cfg.n * cfg.n);
    }

    if (cfg.formula == Formula::QuarticSexticCubicExp || cfg.formula == Formula::QuarticSexticCeiling) {
        Rational quartic_gain = q(1, 12) * q(2 * D, QUARTIC_DENOM);
        Rational sixth_loss = q(1, 360) * q((2 * D) * (2 * D), SEXTIC_DENOM);
        return q(1) - (s1 + t1) / cfg.n
            + quartic_gain * (s0 * s0 + t0 * t0) / (cfg.n * cfg.n)
            - sixth_loss * (pow_q(s1, 3) + pow_q(t1, 3)) / (cfg.n * cfg.n * cfg.n);
    }

    if (cfg.formula == Formula::Tail133) {
        Rational quartic_gain = q(1, 18) * q(2 * D, QUARTIC_DENOM);
        return q(1) - (s1 + t1) / cfg.n
            + quartic_gain * (s0 * s0 + t0 * t0) / (cfg.n * cfg.n);
    }

    const Rational cos0 = q(199909270, 100000000);
    const Rational cos1 = q(-99730681, 100000000);
    const Rational cos2 = q(8092566, 100000000);
    const Rational cos3 = q(-216849, 100000000);
    return (q(16) + q(240) * cos0) / (2 * D)
        + cos1 * (s1 + t1) / cfg.n
        + cos2 * q(2 * D, QUARTIC_DENOM) * (s0 * s0 + t0 * t0) / (cfg.n * cfg.n)
        + cos3 * q((2 * D) * (2 * D), SEXTIC_DENOM) * (pow_q(s1, 3) + pow_q(t1, 3))
            / (cfg.n * cfg.n * cfg.n);
}

Rational local_direct_factor_lower(const Rational& h) {
    return q((2 * D) * (2 * D)) * h * h - 4;
}

Rational radial_exp_lower(const CaseConfig& cfg, const Rational& value) {
    if (cfg.formula == Formula::QuarticSexticCeiling || cfg.formula == Formula::Tail133) {
        return exp_ceiling_lower(value, q(11, 4));
    }
    if (cfg.formula == Formula::QuarticSexticCubicExp) {
        return exp_fractional_cubic_lower(value, q(1457, 536));
    }
    return exp_fractional_septic_lower(value, q(1457, 536));
}

Rational sine_exponent(const CaseConfig& cfg, const Cell& cell) {
    const Rational s1 = cell.s1;
    const Rational t1 = cell.t1;

    if (cfg.formula == Formula::Cubic && cfg.sine_denom != q(1)) {
        return q(4 * D, cfg.n) * (s1 + t1) / cfg.sine_denom;
    }

    if (cfg.formula == Formula::CubicCap14) {
        return 2 * q(1, 24) * q(2 * D) * (s1 + t1) / cfg.n
            + 2 * q(1, 2880) * q((2 * D) * (2 * D), QUARTIC_DENOM)
                * (s1 * s1 + t1 * t1) / (cfg.n * cfg.n)
            + 2 * q(1, 103800) * q((2 * D) * (2 * D), SEXTIC_DENOM)
                * (pow_q(s1, 3) + pow_q(t1, 3)) / (cfg.n * cfg.n * cfg.n);
    }

    if (cfg.formula == Formula::CubicN91) {
        return 2 * q(1, 24) * q(2 * D) * (s1 + t1) / cfg.n
            + 2 * q(1, 2880) * q((2 * D) * (2 * D), QUARTIC_DENOM)
                * (s1 * s1 + t1 * t1) / (cfg.n * cfg.n)
            + 2 * q(1, 115000) * q((2 * D) * (2 * D), SEXTIC_DENOM)
                * (pow_q(s1, 3) + pow_q(t1, 3)) / (cfg.n * cfg.n * cfg.n);
    }

    if (cfg.formula == Formula::CubicWide93) {
        return 2 * q(1, 24) * q(2 * D) * (s1 + t1) / cfg.n
            + 2 * q(1, 2000) * q((2 * D) * (2 * D), QUARTIC_DENOM)
                * (s1 * s1 + t1 * t1) / (cfg.n * cfg.n);
    }

    if (cfg.formula == Formula::QuarticMatched119 || cfg.formula == Formula::QuarticEighth121) {
        return q(4 * D, cfg.n) * (s1 + t1) / q(81, 4);
    }

    if (cfg.formula == Formula::QuarticSexticCubicExp || cfg.formula == Formula::QuarticSexticCeiling
        || cfg.formula == Formula::Tail133) {
        return q(4 * D, 1) * (s1 + t1) / (cfg.sine_denom * cfg.n);
    }

    return 2 * q(1, 24) * q(2 * D) * (s1 + t1) / cfg.n
        + 2 * q(1, 2100) * q((2 * D) * (2 * D), QUARTIC_DENOM)
            * (s1 * s1 + t1 * t1) / (cfg.n * cfg.n);
}

Rational sine_density_lower(const CaseConfig& cfg, const Rational& exponent) {
    if (cfg.formula == Formula::QuarticSexticCeiling || cfg.formula == Formula::Tail133) {
        return exp_ceiling_lower(exponent, q(11, 4));
    }
    if (cfg.formula == Formula::QuarticSexticCubicExp) {
        return exp_fractional_cubic_lower(exponent, q(1457, 536));
    }
    return exp_fractional_septic_lower(exponent, q(1457, 536));
}

Rational rectangle_fraction_lower(const CaseConfig& cfg, const Cell& cell, const Rational& gap) {
    Rational s_integral = (pow_q(cell.s1, SHAPE) - pow_q(cell.s0, SHAPE)) / SHAPE;
    Rational t_integral = (pow_q(cell.t1, SHAPE) - pow_q(cell.t0, SHAPE)) / SHAPE;
    Rational out = gap * gap * radial_exp_lower(cfg, cell.s1 + cell.t1)
        * s_integral * t_integral * gamma_inv_square() / SHAPE;
    out.canonicalize();
    return out;
}

std::vector<Cell> candidate_cells(const CaseConfig& cfg) {
    std::vector<Cell> cells;
    Rational cap = cap_value(cfg);
    auto limit_with_cap = [&cap](const Rational& endpoint, bool use_cap) {
        if (endpoint < 0) return cap;
        if (use_cap && cap < endpoint) return cap;
        return endpoint;
    };
    Rational s_limit = limit_with_cap(cfg.s_end, cfg.s_end_is_cap);
    Rational t_limit = limit_with_cap(cfg.t_end, cfg.t_end_is_cap);

    for (Rational s0 = cfg.s_start; s0 + cfg.grid <= s_limit; s0 += cfg.grid) {
        for (Rational t0 = cfg.t_start; t0 + cfg.grid <= t_limit; t0 += cfg.grid) {
            Cell cell{s0, s0 + cfg.grid, t0, t0 + cfg.grid};
            Rational gap = gap_lower(cfg, cell);
            if (gap <= 0) continue;
            Rational h = local_h_lower(cfg, cell);
            if (h <= 0) continue;
            if (local_direct_factor_lower(h) <= 0) continue;
            cells.push_back(std::move(cell));
        }
    }
    return cells;
}

Stats compute_stats(const CaseConfig& cfg, const std::vector<Cell>& cells) {
    const Rational a0 = a0_lower();
    const int threads = omp_get_max_threads();
    std::vector<Stats> partial(static_cast<std::size_t>(threads));
    auto started = std::chrono::steady_clock::now();
    std::vector<double> last_print(static_cast<std::size_t>(threads), 0.0);

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
#pragma omp for schedule(dynamic, 16)
        for (std::size_t index = 0; index < cells.size(); ++index) {
            const Cell& cell = cells[index];
            Rational gap = gap_lower(cfg, cell);
            Rational h = local_h_lower(cfg, cell);
            Rational direct = local_direct_factor_lower(h);
            Rational rect = rectangle_fraction_lower(cfg, cell, gap);
            Rational exponent = sine_exponent(cfg, cell);
            Rational term = a0 * pow_q(q(2 * D), cfg.n) / pow_q(q(cfg.n), ALPHA)
                * rect * pow_q(h, cfg.n) * direct * sine_density_lower(cfg, exponent);
            term.canonicalize();
            partial[static_cast<std::size_t>(tid)].add(term, rect, gap, h, direct, exponent);

            if ((index + 1) % 5000 == 0) {
                auto now = std::chrono::steady_clock::now();
                double elapsed = std::chrono::duration<double>(now - started).count();
                if (elapsed - last_print[static_cast<std::size_t>(tid)] >= 15.0) {
                    last_print[static_cast<std::size_t>(tid)] = elapsed;
#pragma omp critical
                    {
                        std::cout << "progress n=" << cfg.n << " processed~" << (index + 1)
                                  << "/" << cells.size() << " elapsed=" << elapsed << "s"
                                  << std::endl;
                    }
                }
            }
        }
    }

    Stats total;
    for (const Stats& item : partial) total.merge(item);
    return total;
}

std::map<std::string, Rational> read_negative_bounds(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("could not open negative-bound table: " + path);
    std::map<std::string, Rational> out;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> parts;
        std::stringstream stream(line);
        std::string part;
        while (std::getline(stream, part, '|')) parts.push_back(part);
        if (parts.size() != 4) throw std::runtime_error("bad negative-bound row: " + line);
        out.emplace(parts[0], parse_q(parts[3]));
    }
    return out;
}

std::vector<CaseConfig> cases() {
    std::vector<CaseConfig> out;
    auto cubic = [](std::string label, int n, int cells, Rational s_end, Rational t_end) {
        CaseConfig cfg;
        cfg.label = std::move(label);
        cfg.n = n;
        cfg.expected_cells = cells;
        cfg.formula = Formula::Cubic;
        cfg.s_end = s_end;
        cfg.t_end = t_end;
        cfg.s_end_is_cap = true;
        cfg.t_end_is_cap = true;
        return cfg;
    };

    for (auto item : std::vector<std::tuple<int, int>>{
             {77, 57817}, {79, 68095}, {81, 79060}, {83, 91138},
             {85, 103495}, {87, 116526}, {89, 130740}}) {
        CaseConfig cfg;
        cfg.label = "n" + std::to_string(std::get<0>(item));
        cfg.n = std::get<0>(item);
        cfg.expected_cells = std::get<1>(item);
        cfg.formula = Formula::CubicCap14;
        cfg.delta = q(14, 3);
        cfg.s_end = q(-1);
        cfg.t_end = q(-1);
        cfg.s_end_is_cap = true;
        cfg.t_end_is_cap = true;
        out.push_back(cfg);
    }

    CaseConfig n91;
    n91.label = "n91";
    n91.n = 91;
    n91.expected_cells = 69502;
    n91.formula = Formula::CubicN91;
    n91.delta = q(13, 3);
    n91.s_end = q(-1);
    n91.t_end = q(-1);
    n91.s_end_is_cap = true;
    n91.t_end_is_cap = true;
    out.push_back(n91);

    CaseConfig n93;
    n93.label = "n93";
    n93.n = 93;
    n93.expected_cells = 53822;
    n93.formula = Formula::CubicWide93;
    n93.delta = q(21, 5);
    n93.s_end = q(248, 5);
    n93.t_end = q(969, 20);
    n93.t_end_is_cap = true;
    out.push_back(n93);

    out.push_back(cubic("n95", 95, 29784, q(919, 20), q(907, 20)));
    out.push_back(cubic("n97", 97, 35869, q(469, 10), q(469, 10)));
    out.push_back(cubic("n99", 99, 42662, q(479, 10), q(479, 10)));
    out.push_back(cubic("n101", 101, 49478, q(49), q(48)));
    out.push_back(cubic("n103", 103, 56387, q(50), q(48)));
    out.push_back(cubic("n105", 105, 63660, q(254, 5), q(48)));
    out.push_back(cubic("n107", 107, 70562, q(207, 4), q(48)));
    out.push_back(cubic("n109", 109, 77454, q(527, 10), q(48)));
    out.push_back(cubic("n111", 111, 84702, q(537, 10), q(48)));
    out.push_back(cubic("n113", 113, 91581, q(1093, 20), q(48)));
    out.push_back(cubic("n115", 115, 98455, q(278, 5), q(48)));

    CaseConfig n117 = cubic("n117", 117, 111427, q(283, 5), q(50));
    n117.formula = Formula::Cubic;
    n117.sine_denom = q(81, 4);
    out.push_back(n117);

    CaseConfig n119;
    n119.label = "n119";
    n119.n = 119;
    n119.expected_cells = 15000;
    n119.formula = Formula::QuarticMatched119;
    n119.grid = q(1, 10);
    n119.s_start = q(45);
    n119.s_end = q(60);
    n119.t_start = q(30);
    n119.t_end = q(42);
    n119.s_end_is_cap = true;
    n119.t_end_is_cap = false;
    n119.gap6_lower = q(1, 475);
    out.push_back(n119);

    CaseConfig n121 = n119;
    n121.label = "n121";
    n121.n = 121;
    n121.expected_cells = 74053;
    n121.formula = Formula::QuarticEighth121;
    n121.grid = q(1, 30);
    out.push_back(n121);

    auto sextic_cubic = [](std::string label, int n, int cells, Rational delta,
                           Rational grid, Rational s_start, Rational s_end,
                           Rational t_start, Rational t_end,
                           Rational gap6, Rational sine_denom) {
        CaseConfig cfg;
        cfg.label = std::move(label);
        cfg.n = n;
        cfg.expected_cells = cells;
        cfg.formula = Formula::QuarticSexticCubicExp;
        cfg.delta = delta;
        cfg.grid = grid;
        cfg.s_start = s_start;
        cfg.s_end = s_end;
        cfg.t_start = t_start;
        cfg.t_end = t_end;
        cfg.s_end_is_cap = true;
        cfg.t_end_is_cap = false;
        cfg.gap6_lower = gap6;
        cfg.sine_denom = sine_denom;
        return cfg;
    };

    out.push_back(sextic_cubic("n123", 123, 28121, q(4), q(1, 40), q(53), q(59), q(69, 2), q(75, 2), q(1, 480), q(20)));
    out.push_back(sextic_cubic("n125", 125, 8497, q(7, 2), q(1, 40), q(44), q(49), q(31), q(67, 2), q(1, 450), q(21)));
    out.push_back(sextic_cubic("n127", 127, 10528, q(10, 3), q(1, 40), q(40), q(-1), q(29), q(32), q(1, 450), q(21)));
    out.push_back(sextic_cubic("n129", 129, 12719, q(19, 6), q(1, 40), q(37), q(-1), q(26), q(31), q(1, 432), q(21)));

    CaseConfig n131 = sextic_cubic("n131", 131, 17797, q(19, 6), q(1, 40), q(37), q(-1), q(26), q(31), q(1, 432), q(21));
    n131.formula = Formula::QuarticSexticCeiling;
    out.push_back(n131);

    CaseConfig tail;
    tail.label = "n133";
    tail.n = 133;
    tail.expected_cells = 13558;
    tail.formula = Formula::Tail133;
    tail.grid = q(1, 40);
    tail.s_start = q(37);
    tail.s_end = q(-1);
    tail.t_start = q(26);
    tail.t_end = q(31);
    tail.t_end_is_cap = false;
    tail.sine_denom = q(21);
    out.push_back(tail);
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    std::string negative_path = "ginibre_q3/certificates/exceptional_rect/e8_rect_negative_bounds.tsv";
    std::string only;
    for (int index = 1; index < argc; ++index) {
        std::string arg = argv[index];
        if (arg == "--negative-bounds") {
            if (index + 1 >= argc) throw std::runtime_error("missing --negative-bounds value");
            negative_path = argv[++index];
        } else if (arg == "--only") {
            if (index + 1 >= argc) throw std::runtime_error("missing --only value");
            only = argv[++index];
        } else {
            throw std::runtime_error("unknown argument: " + arg);
        }
    }

    std::map<std::string, Rational> negative_bounds = read_negative_bounds(negative_path);
    bool all_ok = true;
    std::cout << "COMMAND: direct_chain_rect_e8_parallel_certificate --negative-bounds "
              << negative_path;
    if (!only.empty()) std::cout << " --only " << only;
    std::cout << std::endl;
    std::cout << "OpenMP max threads=" << omp_get_max_threads() << std::endl;
    const bool lattice_ok = e8_weyl_lattice_normalization_check();
    all_ok = all_ok && lattice_ok;
    std::cout << "e8_normalized_weyl_lattice_factor_1_and_minimum_2: "
              << (lattice_ok ? "OK" : "FAIL") << std::endl;

    for (const CaseConfig& cfg : cases()) {
        if (!only.empty() && cfg.label != only) continue;
        auto neg_it = negative_bounds.find(cfg.label);
        if (neg_it == negative_bounds.end()) {
            throw std::runtime_error("missing negative bound for " + cfg.label);
        }
        Rational negative = neg_it->second;
        std::cout << "case " << cfg.label << " n=" << cfg.n << " start" << std::endl;
        std::vector<Cell> cells = candidate_cells(cfg);
        bool cell_count_ok = static_cast<int>(cells.size()) == cfg.expected_cells;
        bool transpose_disjoint = std::all_of(
            cells.begin(), cells.end(), [](const Cell& cell) { return cell.s0 >= cell.t1; });
        Stats stats = compute_stats(cfg, cells);
        Rational ratio = stats.positive_total / negative;
        Rational best_ratio = stats.best_term / negative;
        bool ratio_ok = ratio > 1;
        bool gaps_ok = stats.min_gap > 0;
        bool h_ok = stats.min_h > 0;
        bool direct_ok = stats.min_direct > 0;
        bool cap_ok = true;
        Rational cap = cap_value(cfg);
        for (const Cell& cell : cells) {
            if (!(cell.s1 < cap && cell.t1 < cap)) {
                cap_ok = false;
                break;
            }
        }
        all_ok = all_ok && cell_count_ok && transpose_disjoint
            && ratio_ok && gaps_ok && h_ok && direct_ok && cap_ok;

        std::cout << "group=E8 n=" << cfg.n << std::endl;
        std::cout << "retained_cells=" << cells.size() << " expected=" << cfg.expected_cells << std::endl;
        std::cout << "cap=" << cap << std::endl;
        std::cout << "negative_bound=" << negative << std::endl;
        std::cout << "log10_union_fraction_lower=" << fixed2(log10_q(stats.rect_total)) << std::endl;
        std::cout << "gap" << cfg.n << "_min=" << stats.min_gap << std::endl;
        std::cout << "h" << cfg.n << "_min=" << stats.min_h
                  << " direct_factor" << cfg.n << "_min=" << stats.min_direct << std::endl;
        std::cout << "sine_exponent_max=" << stats.max_sine_exponent << std::endl;
        std::cout << "log10_positive" << cfg.n << "_lower=" << fixed2(log10_q(stats.positive_total)) << std::endl;
        std::cout << "log10_negative" << cfg.n << "_upper=" << fixed2(log10_q(negative)) << std::endl;
        std::cout << "log10_ratio" << cfg.n << "_lower=" << fixed2(log10_q(ratio)) << std::endl;
        std::cout << "log10_best_cell_ratio" << cfg.n << "=" << fixed2(log10_q(best_ratio)) << std::endl;
        std::cout << "cell_count_" << cfg.expected_cells << ": " << (cell_count_ok ? "OK" : "FAIL") << std::endl;
        std::cout << "transposed_cells_disjoint: "
                  << (transpose_disjoint ? "OK" : "FAIL") << std::endl;
        std::cout << "all_cells_in_cap: " << (cap_ok ? "OK" : "FAIL") << std::endl;
        std::cout << "positive_gaps: " << (gaps_ok ? "OK" : "FAIL") << std::endl;
        std::cout << "positive_h" << cfg.n << ": " << (h_ok ? "OK" : "FAIL") << std::endl;
        std::cout << "positive_direct" << cfg.n << ": " << (direct_ok ? "OK" : "FAIL") << std::endl;
        std::cout << "ratio" << cfg.n << "_gt_1: " << (ratio_ok ? "OK" : "FAIL") << std::endl;

        if (cfg.formula == Formula::Tail133) {
            Rational h0 = stats.min_h;
            Rational step = q(D, C_NEG) * h0 * q(D, C_NEG) * h0
                * pow_q(q(cfg.n, cfg.n + 2), ALPHA);
            bool step_ok = step > 1;
            all_ok = all_ok && step_ok;
            std::cout << "log10_odd_step_ratio_lower=" << fixed2(log10_q(step)) << std::endl;
            std::cout << "odd_step_ratio_gt_1: " << (step_ok ? "OK" : "FAIL") << std::endl;
        }
        std::cout << "case " << cfg.label << " result="
                  << ((cell_count_ok && transpose_disjoint && ratio_ok && gaps_ok
                       && h_ok && direct_ok && cap_ok) ? "OK" : "FAIL")
                  << std::endl;
    }

    std::cout << "RESULT: " << (all_ok ? "ALL PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}
