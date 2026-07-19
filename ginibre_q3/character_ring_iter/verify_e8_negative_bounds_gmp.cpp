#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using Rational = mpq_class;
using Polynomial = std::vector<Rational>;

namespace {

Rational q(const std::string& text) {
    Rational result(text);
    result.canonicalize();
    return result;
}

Rational q(long long numerator, long long denominator = 1) {
    return q(std::to_string(numerator) + "/" + std::to_string(denominator));
}

mpz_class pow_z(unsigned long base, unsigned long exponent) {
    mpz_class result;
    mpz_ui_pow_ui(result.get_mpz_t(), base, exponent);
    return result;
}

mpz_class binomial_z(unsigned long n, unsigned long k) {
    mpz_class result;
    mpz_bin_uiui(result.get_mpz_t(), n, k);
    return result;
}

std::vector<std::string> split(const std::string& text, char delimiter) {
    std::vector<std::string> result;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, delimiter)) result.push_back(item);
    return result;
}

void trim(Polynomial& polynomial) {
    while (!polynomial.empty() && polynomial.back() == 0) polynomial.pop_back();
}

Polynomial add(const Polynomial& left, const Polynomial& right) {
    Polynomial result(std::max(left.size(), right.size()), q(0));
    for (std::size_t i = 0; i < left.size(); ++i) result[i] += left[i];
    for (std::size_t i = 0; i < right.size(); ++i) result[i] += right[i];
    trim(result);
    return result;
}

Polynomial scale(const Polynomial& polynomial, const Rational& factor) {
    Polynomial result = polynomial;
    for (Rational& coefficient : result) coefficient *= factor;
    trim(result);
    return result;
}

Polynomial multiply(const Polynomial& left, const Polynomial& right) {
    if (left.empty() || right.empty()) return {};
    Polynomial result(left.size() + right.size() - 1, q(0));
    for (std::size_t i = 0; i < left.size(); ++i)
        for (std::size_t j = 0; j < right.size(); ++j)
            result[i + j] += left[i] * right[j];
    trim(result);
    return result;
}

Polynomial power(Polynomial base, unsigned int exponent) {
    Polynomial result{q(1)};
    while (exponent > 0) {
        if (exponent & 1U) result = multiply(result, base);
        exponent >>= 1U;
        if (exponent) base = multiply(base, base);
    }
    return result;
}

Rational evaluate(const Polynomial& polynomial, const Rational& point) {
    Rational result = q(0);
    for (auto it = polynomial.rbegin(); it != polynomial.rend(); ++it)
        result = result * point + *it;
    return result;
}

Polynomial compose_linear(
    const Polynomial& polynomial,
    const Rational& linear_scale,
    const Rational& shift
) {
    Polynomial result;
    Polynomial term{q(1)};
    const Polynomial linear{shift, linear_scale};
    for (const Rational& coefficient : polynomial) {
        result = add(result, scale(term, coefficient));
        term = multiply(term, linear);
    }
    return result;
}

std::pair<Polynomial, Polynomial> divide(
    Polynomial numerator,
    Polynomial denominator
) {
    trim(numerator);
    trim(denominator);
    if (denominator.empty()) throw std::runtime_error("polynomial division by zero");
    if (numerator.size() < denominator.size()) return {{}, numerator};
    Polynomial quotient(numerator.size() - denominator.size() + 1, q(0));
    while (!numerator.empty() && numerator.size() >= denominator.size()) {
        const std::size_t shift_by = numerator.size() - denominator.size();
        const Rational coefficient = numerator.back() / denominator.back();
        quotient[shift_by] = coefficient;
        for (std::size_t i = 0; i < denominator.size(); ++i)
            numerator[i + shift_by] -= coefficient * denominator[i];
        trim(numerator);
    }
    trim(quotient);
    return {quotient, numerator};
}

Polynomial derivative(const Polynomial& polynomial) {
    if (polynomial.size() < 2) return {};
    Polynomial result(polynomial.size() - 1, q(0));
    for (std::size_t i = 1; i < polynomial.size(); ++i)
        result[i - 1] = polynomial[i] * static_cast<unsigned long>(i);
    trim(result);
    return result;
}

int sign(const Rational& value) {
    return value > 0 ? 1 : (value < 0 ? -1 : 0);
}

int sign_variations(const std::vector<Rational>& values) {
    int previous = 0;
    int variations = 0;
    for (const Rational& value : values) {
        const int current = sign(value);
        if (current == 0) continue;
        if (previous != 0 && current != previous) ++variations;
        previous = current;
    }
    return variations;
}

std::vector<Polynomial> sturm_sequence(const Polynomial& polynomial) {
    Polynomial first = polynomial;
    trim(first);
    Polynomial second = derivative(first);
    std::vector<Polynomial> sequence{first};
    if (second.empty()) return sequence;
    sequence.push_back(second);
    while (!sequence.back().empty()) {
        auto [unused, remainder] = divide(
            sequence[sequence.size() - 2], sequence.back());
        (void)unused;
        if (remainder.empty()) break;
        sequence.push_back(scale(remainder, q(-1)));
    }
    return sequence;
}

bool sturm_positive(
    const Polynomial& polynomial,
    const Rational& left,
    const Rational& right
) {
    if (evaluate(polynomial, left) <= 0 || evaluate(polynomial, right) <= 0)
        return false;
    const std::vector<Polynomial> sequence = sturm_sequence(polynomial);
    std::vector<Rational> at_left;
    std::vector<Rational> at_right;
    at_left.reserve(sequence.size());
    at_right.reserve(sequence.size());
    for (const Polynomial& item : sequence) {
        at_left.push_back(evaluate(item, left));
        at_right.push_back(evaluate(item, right));
    }
    return sign_variations(at_left) - sign_variations(at_right) == 0;
}

bool bernstein_positive(
    const Polynomial& polynomial,
    const Rational& left,
    const Rational& right
) {
    if (polynomial.empty()) return false;
    const Polynomial unit = compose_linear(polynomial, right - left, left);
    const unsigned long degree = static_cast<unsigned long>(unit.size() - 1);
    for (unsigned long i = 0; i <= degree; ++i) {
        Rational coefficient = q(0);
        for (unsigned long k = 0; k <= i; ++k) {
            if (k >= unit.size()) break;
            coefficient += unit[k]
                * Rational(binomial_z(i, k), binomial_z(degree, k));
        }
        if (coefficient <= 0) return false;
    }
    return true;
}

bool bernstein_subdivide(
    const Polynomial& polynomial,
    const Rational& left,
    const Rational& right,
    int depth,
    int maximum_depth,
    int& deepest
) {
    deepest = std::max(deepest, depth);
    if (bernstein_positive(polynomial, left, right)) return true;
    if (depth == maximum_depth) return false;
    const Rational midpoint = (left + right) / 2;
    return bernstein_subdivide(
               polynomial, left, midpoint, depth + 1, maximum_depth, deepest)
        && bernstein_subdivide(
               polynomial, midpoint, right, depth + 1, maximum_depth, deepest);
}

struct PositivityResult {
    bool passed = false;
    std::string method;
    int bernstein_depth = 0;
};

PositivityResult certify_positive(
    const Polynomial& polynomial,
    const Rational& left,
    const Rational& right
) {
    PositivityResult result;
    if (evaluate(polynomial, left) <= 0 || evaluate(polynomial, right) <= 0)
        return result;
    if (bernstein_subdivide(polynomial, left, right, 0, 8, result.bernstein_depth)) {
        result.passed = true;
        result.method = "Bernstein";
        return result;
    }
    result.passed = sturm_positive(polynomial, left, right);
    result.method = "Sturm";
    return result;
}

struct TableRow {
    int cells = 0;
    Rational bound = q(0);
};

std::map<int, TableRow> read_bound_table(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("could not open E8 bound table: " + path);
    std::map<int, TableRow> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> parts = split(line, '|');
        if (parts.size() != 4) throw std::runtime_error("bad E8 bound row: " + line);
        const int n = std::stoi(parts[1]);
        if (parts[0] != "n" + std::to_string(n) || rows.count(n) != 0)
            throw std::runtime_error("bad or duplicate E8 bound label");
        rows.emplace(n, TableRow{std::stoi(parts[2]), q(parts[3])});
    }
    return rows;
}

struct FactorCase {
    int n = 0;
    int moment_weight = 0;
    Rational bound_scale = q(0);
    bool divide_negative_endpoint = false;
    Polynomial factor;
};

std::vector<FactorCase> read_factor_table(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("could not open E8 factor table: " + path);
    std::vector<FactorCase> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> parts = split(line, '|');
        if (parts.size() != 5) throw std::runtime_error("bad E8 factor row");
        FactorCase row;
        row.n = std::stoi(parts[0]);
        row.moment_weight = std::stoi(parts[1]);
        row.bound_scale = q(parts[2]);
        row.divide_negative_endpoint = std::stoi(parts[3]) != 0;
        for (const std::string& coefficient : split(parts[4], ','))
            row.factor.push_back(q(coefficient));
        trim(row.factor);
        if (row.factor.empty()) throw std::runtime_error("empty E8 factor polynomial");
        rows.push_back(std::move(row));
    }
    return rows;
}

std::map<int, mpz_class> read_m_equals(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("could not open moment source: " + path);
    const std::regex expression(R"(m_\s*(\d+)\s*=\s*(-?\d+))");
    std::map<int, mpz_class> moments;
    std::string line;
    while (std::getline(input, line)) {
        std::smatch match;
        if (std::regex_search(line, match, expression))
            moments[std::stoi(match[1].str())] = mpz_class(match[2].str());
    }
    return moments;
}

std::map<int, mpz_class> read_oeis(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("could not open OEIS moment source: " + path);
    const std::regex expression(R"(^\s*(\d+)\s+(-?\d+)\s*$)");
    std::map<int, mpz_class> moments;
    std::string line;
    while (std::getline(input, line)) {
        std::smatch match;
        if (std::regex_match(line, match, expression))
            moments[std::stoi(match[1].str())] = mpz_class(match[2].str());
    }
    return moments;
}

std::vector<mpz_class> audited_moments(
    const std::string& local_path,
    const std::string& ancillary_path,
    const std::string& oeis_path
) {
    const auto local = read_m_equals(local_path);
    const auto ancillary = read_m_equals(ancillary_path);
    const auto oeis = read_oeis(oeis_path);
    std::vector<mpz_class> moments(101);
    for (int index = 0; index <= 70; ++index) {
        const auto iterator = local.find(index);
        if (iterator == local.end())
            throw std::runtime_error("local E8 moment source stops before m_70");
        moments[static_cast<std::size_t>(index)] = iterator->second;
    }
    for (int index = 71; index <= 100; ++index) {
        const auto iterator = ancillary.find(index);
        if (iterator == ancillary.end())
            throw std::runtime_error("ancillary E8 moment source stops before m_100");
        moments[static_cast<std::size_t>(index)] = iterator->second;
    }
    for (int index = 0; index <= 30; ++index) {
        const auto iterator = oeis.find(index);
        if (iterator == oeis.end() || iterator->second != moments[static_cast<std::size_t>(index)])
            throw std::runtime_error("OEIS/local E8 moment prefix mismatch");
    }
    return moments;
}

mpz_class moment_integral(int power_value, const std::vector<mpz_class>& moments) {
    if (power_value < 0 || power_value + 2 >= static_cast<int>(moments.size()))
        throw std::runtime_error("E8 moment integral exceeds audited source range");
    mpz_class total = 0;
    for (int k = 0; k <= power_value; ++k) {
        const mpz_class coefficient = binomial_z(
            static_cast<unsigned long>(power_value), static_cast<unsigned long>(k));
        total += coefficient * moments[static_cast<std::size_t>(k + 2)]
            * moments[static_cast<std::size_t>(power_value - k)];
        total += coefficient * moments[static_cast<std::size_t>(k)]
            * moments[static_cast<std::size_t>(power_value - k + 2)];
        total -= 2 * coefficient * moments[static_cast<std::size_t>(k + 1)]
            * moments[static_cast<std::size_t>(power_value - k + 1)];
    }
    return total;
}

Rational factor_integral(
    const FactorCase& row,
    const std::vector<mpz_class>& moments
) {
    const Polynomial square = multiply(row.factor, row.factor);
    Rational total = q(0);
    for (std::size_t i = 0; i < square.size(); ++i) {
        const int shift = static_cast<int>(i);
        total += square[i] * Rational(
            moment_integral(row.moment_weight + 2 + shift, moments)
            + 4 * moment_integral(row.moment_weight + shift, moments));
    }
    total.canonicalize();
    return total;
}

struct FactorResult {
    int n = 0;
    bool pointwise_negative = false;
    bool pointwise_positive = false;
    bool bound_match = false;
    bool positive_integral = false;
    std::string negative_method;
    std::string positive_method;
    int negative_depth = 0;
    int positive_depth = 0;
};

FactorResult verify_factor_case(
    const FactorCase& row,
    const TableRow& table,
    const std::vector<mpz_class>& moments
) {
    const Polynomial variable{q(0), q(1)};
    const Polynomial square = multiply(variable, variable);
    const Polynomial plus_four = add(square, Polynomial{q(4)});
    const Polynomial minus_four = add(square, Polynomial{q(-4)});
    const Polynomial four_minus_square = add(Polynomial{q(4)}, scale(square, q(-1)));
    const Rational prefactor = row.bound_scale
        * Rational(pow_z(16, static_cast<unsigned long>(row.n - row.moment_weight)));

    const Polynomial reflected_factor = compose_linear(row.factor, q(-1), q(0));
    Polynomial negative = add(
        scale(multiply(plus_four, multiply(reflected_factor, reflected_factor)), prefactor),
        scale(multiply(
            power(variable, static_cast<unsigned int>(row.n - row.moment_weight)),
            minus_four), q(-1)));
    if (row.divide_negative_endpoint) {
        auto [quotient, remainder] = divide(negative, Polynomial{q(16), q(-1)});
        if (!remainder.empty()) {
            FactorResult failed;
            failed.n = row.n;
            return failed;
        }
        negative = std::move(quotient);
    }
    const PositivityResult negative_check = certify_positive(negative, q(2), q(16));

    const Polynomial positive = add(
        scale(multiply(plus_four, multiply(row.factor, row.factor)), prefactor),
        scale(multiply(
            power(variable, static_cast<unsigned int>(row.n - row.moment_weight)),
            four_minus_square), q(-1)));
    const PositivityResult positive_check = certify_positive(positive, q(0), q(2));

    const Rational integral = factor_integral(row, moments);
    const Rational computed_bound = prefactor * integral;
    FactorResult result;
    result.n = row.n;
    result.pointwise_negative = negative_check.passed;
    result.pointwise_positive = positive_check.passed;
    result.bound_match = computed_bound == table.bound;
    result.positive_integral = integral > 0;
    result.negative_method = negative_check.method;
    result.positive_method = positive_check.method;
    result.negative_depth = negative_check.bernstein_depth;
    result.positive_depth = positive_check.bernstein_depth;
    return result;
}

const std::map<int, int>& expected_cells() {
    static const std::map<int, int> values{
        {77, 57817}, {79, 68095}, {81, 79060}, {83, 91138},
        {85, 103495}, {87, 116526}, {89, 130740}, {91, 69502},
        {93, 53822}, {95, 29784}, {97, 35869}, {99, 42662},
        {101, 49478}, {103, 56387}, {105, 63660}, {107, 70562},
        {109, 77454}, {111, 84702}, {113, 91581}, {115, 98455},
        {117, 111427}, {119, 15000}, {121, 74053}, {123, 28121},
        {125, 8497}, {127, 10528}, {129, 12719}, {131, 17797},
        {133, 13558},
    };
    return values;
}

}  // namespace

int main(int argc, char** argv) {
    std::string table_path = "../certificates/exceptional_rect/e8_rect_negative_bounds.tsv";
    std::string factor_path = "../certificates/exceptional_rect/e8_negative_majorant_factors.tsv";
    std::string local_moments = "logs/e8_70.log";
    std::string ancillary_moments = "../references/arxiv_2412_21189_e8_m71_m100.txt";
    std::string oeis_moments = "../references/oeis_A179663.txt";
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        auto value = [&](const std::string& name) {
            if (index + 1 >= argc) throw std::runtime_error("missing value for " + name);
            return std::string(argv[++index]);
        };
        if (argument == "--table") table_path = value(argument);
        else if (argument == "--factors") factor_path = value(argument);
        else if (argument == "--local-moments") local_moments = value(argument);
        else if (argument == "--ancillary-moments") ancillary_moments = value(argument);
        else if (argument == "--oeis-moments") oeis_moments = value(argument);
        else throw std::runtime_error("unknown argument: " + argument);
    }

    const std::map<int, TableRow> table = read_bound_table(table_path);
    const std::vector<FactorCase> factors = read_factor_table(factor_path);
    const std::vector<mpz_class> moments = audited_moments(
        local_moments, ancillary_moments, oeis_moments);
    bool all_ok = table.size() == expected_cells().size() && factors.size() == 16;
    for (const auto& [n, cells] : expected_cells()) {
        const auto iterator = table.find(n);
        all_ok = all_ok && iterator != table.end() && iterator->second.cells == cells
            && iterator->second.bound > 0;
    }

    std::vector<FactorCase> active = factors;
    active.push_back(FactorCase{109, 12, q(1), false, Polynomial{q(1)}});
    active.push_back(FactorCase{111, 2, q(1), false, Polynomial{q(1)}});
    active.push_back(FactorCase{113, 2, q(1), false, Polynomial{q(1)}});
    std::sort(active.begin(), active.end(), [](const FactorCase& a, const FactorCase& b) {
        return a.n < b.n;
    });
    std::vector<FactorResult> results(active.size());

#pragma omp parallel for schedule(dynamic, 1)
    for (std::size_t index = 0; index < active.size(); ++index) {
        const auto iterator = table.find(active[index].n);
        if (iterator != table.end())
            results[index] = verify_factor_case(active[index], iterator->second, moments);
    }

    std::cout << "E8 exact negative-bound audit (C++/GMP)" << std::endl;
    std::cout << "moment_sources_m0_m100_and_oeis_overlap: OK" << std::endl;
    for (const FactorResult& result : results) {
        const bool passed = result.pointwise_negative && result.pointwise_positive
            && result.bound_match && result.positive_integral;
        all_ok = all_ok && passed;
        std::cout << "n=" << result.n
                  << " negative=" << (result.pointwise_negative ? "OK" : "FAIL")
                  << "(" << result.negative_method << ",depth=" << result.negative_depth << ")"
                  << " positive=" << (result.pointwise_positive ? "OK" : "FAIL")
                  << "(" << result.positive_method << ",depth=" << result.positive_depth << ")"
                  << " moment_integral=" << (result.positive_integral ? "OK" : "FAIL")
                  << " bound_match=" << (result.bound_match ? "OK" : "FAIL")
                  << std::endl;
    }

    for (int n = 115; n <= 133; n += 2) {
        const Rational endpoint_bound = Rational(2 * (16 * 16 + 4))
            * Rational(pow_z(16, static_cast<unsigned long>(n)));
        const auto iterator = table.find(n);
        const bool passed = iterator != table.end() && iterator->second.bound == endpoint_bound;
        all_ok = all_ok && passed;
        std::cout << "n=" << n << " endpoint_bound_match="
                  << (passed ? "OK" : "FAIL") << std::endl;
    }

    std::cout << "E8_NEGATIVE_BOUND_TABLE_ROWS=" << table.size() << std::endl;
    std::cout << "E8_NEGATIVE_BOUND_AUDIT_GMP: "
              << (all_ok ? "ALL PASS" : "FAIL") << std::endl;
    return all_ok ? 0 : 1;
}
