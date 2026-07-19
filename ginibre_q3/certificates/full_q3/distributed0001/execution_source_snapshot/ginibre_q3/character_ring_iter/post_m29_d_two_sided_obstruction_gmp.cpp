#include <gmpxx.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

constexpr int RANK = 98;
constexpr int CHAIN_M = 102;
constexpr double LOG10_2 = 0.301029995663981195213738894724493027;

double log10_abs(const mpz_class& value) {
    if (value == 0) return -std::numeric_limits<double>::infinity();
    mpz_class absolute = value >= 0 ? value : -value;
    long exponent = 0;
    const double mantissa = mpz_get_d_2exp(&exponent, absolute.get_mpz_t());
    return std::log10(mantissa) + static_cast<double>(exponent) * LOG10_2;
}

std::vector<mpz_class> stable_moments(int maximum) {
    std::vector<mpz_class> moments(static_cast<std::size_t>(maximum + 1));
    moments[0] = 1;
    if (maximum >= 1) moments[1] = 0;
    for (int index = 1; index < maximum; ++index) {
        mpz_class next = index * moments[static_cast<std::size_t>(index)]
            + index * moments[static_cast<std::size_t>(index - 1)];
        if (index >= 2) {
            next -= (mpz_class(index) * (index - 1) / 2)
                * moments[static_cast<std::size_t>(index - 2)];
        }
        if (next <= 0) {
            std::cerr << "stable recurrence lost positivity\n";
            std::exit(2);
        }
        moments[static_cast<std::size_t>(index + 1)] = std::move(next);
    }
    return moments;
}

std::vector<mpz_class> binomial_row(int n) {
    std::vector<mpz_class> row(static_cast<std::size_t>(n + 1));
    row[0] = 1;
    for (int k = 0; k < n; ++k) {
        row[static_cast<std::size_t>(k + 1)] =
            row[static_cast<std::size_t>(k)] * (n - k) / (k + 1);
    }
    return row;
}

void add_if_valid(
    mpz_class& coefficient,
    const std::vector<mpz_class>& row,
    int index,
    int factor
) {
    if (index < 0 || index >= static_cast<int>(row.size())) return;
    coefficient += factor * row[static_cast<std::size_t>(index)];
}

mpz_class pair_coefficient(
    const std::vector<mpz_class>& row,
    int first,
    int second,
    int scale
) {
    const int factor = 2 * scale;
    mpz_class coefficient = 0;
    add_if_valid(coefficient, row, first - 2, factor);
    if (second != first) add_if_valid(coefficient, row, second - 2, factor);
    add_if_valid(coefficient, row, first - 1, -factor);
    if (second != first) add_if_valid(coefficient, row, second - 1, -factor);
    return coefficient;
}

void accumulate_diagonal(
    int n,
    int scale,
    const std::vector<mpz_class>& stable,
    std::vector<mpz_class>& linear,
    mpz_class& stable_value,
    mpz_class& one_sided_quadratic,
    mpz_class& two_sided_quadratic
) {
    const std::vector<mpz_class> row = binomial_row(n);
    const int pair_sum = n + 2;
    for (int first = 0; first <= pair_sum / 2; ++first) {
        const int second = pair_sum - first;
        const mpz_class coefficient = pair_coefficient(row, first, second, scale);
        if (coefficient == 0) continue;

        const mpz_class stable_product =
            stable[static_cast<std::size_t>(first)]
            * stable[static_cast<std::size_t>(second)];
        stable_value += coefficient * stable_product;

        if (first >= RANK && second >= RANK) {
            if (coefficient < 0) {
                one_sided_quadratic += coefficient * stable_product;
            }
            if (first == second) {
                if (coefficient < 0) {
                    two_sided_quadratic += coefficient * stable_product;
                }
            } else {
                const mpz_class absolute =
                    coefficient >= 0 ? coefficient : -coefficient;
                two_sided_quadratic -= absolute * stable_product;
            }
        }

        if (first == second) {
            linear[static_cast<std::size_t>(first)] +=
                2 * coefficient * stable[static_cast<std::size_t>(first)];
        } else {
            linear[static_cast<std::size_t>(first)] +=
                coefficient * stable[static_cast<std::size_t>(second)];
            linear[static_cast<std::size_t>(second)] +=
                coefficient * stable[static_cast<std::size_t>(first)];
        }
    }
}

}  // namespace

int main() {
    const int maximum = 2 * CHAIN_M + 5;
    const std::vector<mpz_class> stable = stable_moments(maximum);
    std::vector<mpz_class> linear(static_cast<std::size_t>(maximum + 1));
    mpz_class stable_value = 0;
    mpz_class one_sided_quadratic = 0;
    mpz_class two_sided_quadratic = 0;

    accumulate_diagonal(
        2 * CHAIN_M + 3,
        1,
        stable,
        linear,
        stable_value,
        one_sided_quadratic,
        two_sided_quadratic
    );
    accumulate_diagonal(
        2 * CHAIN_M + 1,
        -4,
        stable,
        linear,
        stable_value,
        one_sided_quadratic,
        two_sided_quadratic
    );

    mpz_class one_sided_linear = 0;
    mpz_class two_sided_linear = 0;
    for (int index = RANK; index <= maximum; ++index) {
        const mpz_class& coefficient = linear[static_cast<std::size_t>(index)];
        if (coefficient < 0) {
            one_sided_linear += coefficient * stable[static_cast<std::size_t>(index)];
        }
        const mpz_class absolute = coefficient >= 0 ? coefficient : -coefficient;
        two_sided_linear -= absolute * stable[static_cast<std::size_t>(index)];
    }

    const mpz_class one_sided =
        stable_value + one_sided_linear + one_sided_quadratic;
    const mpz_class two_sided =
        stable_value + two_sided_linear + two_sided_quadratic;

    std::cout << "D two-sided stable-envelope obstruction GMP audit\n"
              << "rank=D_" << RANK
              << " chain_m=" << CHAIN_M
              << " odd_n=" << (2 * CHAIN_M + 3)
              << " maximum_moment=" << maximum << '\n'
              << std::setprecision(18)
              << "stable_sign=" << sgn(stable_value)
              << " stable_log10=" << log10_abs(stable_value) << '\n'
              << "one_sided_sign=" << sgn(one_sided)
              << " one_sided_log10=" << log10_abs(one_sided) << '\n'
              << "two_sided_sign=" << sgn(two_sided)
              << " two_sided_log10=" << log10_abs(two_sided) << '\n';

    if (one_sided <= 0 || two_sided >= 0) {
        std::cerr << "expected obstruction signs were not reproduced\n";
        return 1;
    }
    std::cout << "RESULT: one-sided envelope passes; sign-free two-sided envelope fails\n";
    return 0;
}
