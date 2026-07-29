#define main anchored_lobe_probe_main
#include "probe_su2_anchored_lobe_all_powers.cpp"
#undef main

#include <array>
#include <cstddef>
#include <iostream>
#include <map>
#include <vector>

namespace {

Integer binomial(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    k = std::min(k, n - k);
    Integer value = 1;
    for (int index = 1; index <= k; ++index) {
        value *= n - k + index;
        value /= index;
    }
    return value;
}

Rational coefficient_at(
    const RationalPolynomial& polynomial,
    std::size_t degree
) {
    if (degree >= polynomial.numerator.coefficients().size()) {
        return 0;
    }
    return polynomial.numerator.coefficients()[degree];
}

using Grid = std::vector<std::vector<RationalPolynomial>>;

Grid chamber_grid(bool high_chamber, int grid_degree, int power) {
    Grid grid(
        static_cast<std::size_t>(grid_degree + 1),
        std::vector<RationalPolynomial>(
            static_cast<std::size_t>(grid_degree + 1)
        )
    );
    for (int x = 0; x <= grid_degree; ++x) {
        for (int y = 0; y <= grid_degree; ++y) {
            // Low: S=x+1, Q=x+y+1.
            // High: B=S-Q-1=x, Q=x+y+1.
            const int q_half = x + y + 1;
            const int target_half =
                high_chamber ? 2 * x + y + 2 : x + 1;
            grid[static_cast<std::size_t>(x)]
                [static_cast<std::size_t>(y)] =
                    first_block_determinant(
                        q_half,
                        target_half,
                        power
                    );
        }
    }
    int denominator_power = 0;
    for (const auto& row : grid) {
        for (const RationalPolynomial& polynomial : row) {
            denominator_power = std::max(
                denominator_power,
                polynomial.denominator_power
            );
        }
    }
    for (auto& row : grid) {
        for (RationalPolynomial& polynomial : row) {
            polynomial = align_denominator(
                std::move(polynomial),
                denominator_power
            );
        }
    }
    return grid;
}

Rational forward_coefficient(
    const Grid& grid,
    int x_order,
    int y_order,
    std::size_t z_degree
) {
    Rational result = 0;
    for (int x = 0; x <= x_order; ++x) {
        for (int y = 0; y <= y_order; ++y) {
            Rational term{
                binomial(x_order, x) * binomial(y_order, y)
            };
            term *= coefficient_at(
                grid[static_cast<std::size_t>(x)]
                    [static_cast<std::size_t>(y)],
                z_degree
            );
            if ((x_order - x + y_order - y) % 2 != 0) {
                result -= term;
            } else {
                result += term;
            }
        }
    }
    return result;
}

Polynomial binomial_basis(int order) {
    Polynomial result{Rational{1}};
    for (int root = 0; root < order; ++root) {
        result =
            result
            * Polynomial{
                std::vector<Rational>{
                    Rational{-root},
                    Rational{1}
                }
            };
        result = Rational{1, root + 1} * result;
    }
    return result;
}

bool report_chamber(bool high_chamber, int grid_degree, int power) {
    const Grid grid = chamber_grid(
        high_chamber,
        grid_degree,
        power
    );
    std::size_t max_z_degree = 0;
    for (const auto& row : grid) {
        for (const RationalPolynomial& polynomial : row) {
            max_z_degree = std::max(
                max_z_degree,
                polynomial.numerator.coefficients().size()
            );
        }
    }

    std::size_t negative = 0;
    std::size_t printed_negative = 0;
    std::size_t nonzero = 0;
    std::size_t degree_bound_violations = 0;
    int max_x_order = 0;
    int max_y_order = 0;
    int max_total_order = 0;
    for (std::size_t z_degree = 0;
         z_degree < max_z_degree;
         ++z_degree) {
        for (int x_order = 0; x_order <= grid_degree; ++x_order) {
            for (int y_order = 0; y_order <= grid_degree; ++y_order) {
                const Rational coefficient = forward_coefficient(
                    grid,
                    x_order,
                    y_order,
                    z_degree
                );
                if (coefficient != 0) {
                    ++nonzero;
                    max_x_order = std::max(max_x_order, x_order);
                    max_y_order = std::max(max_y_order, y_order);
                    max_total_order = std::max(
                        max_total_order,
                        x_order + y_order
                    );
                    if (x_order + y_order > 4 * power + 2) {
                        ++degree_bound_violations;
                    }
                }
                if (coefficient < 0) {
                    ++negative;
                    if (printed_negative < 3) {
                        ++printed_negative;
                        std::cout
                            << "SU2_ANCHORED_LOBE_FIRST_BLOCK_POWERS"
                            << " chamber="
                            << (high_chamber ? "high" : "low")
                            << " negative_binomial"
                            << " z_degree=" << z_degree
                            << " x_order=" << x_order
                            << " y_order=" << y_order
                            << " coefficient=" << coefficient
                            << '\n';
                    }
                }
            }
        }
    }

    std::size_t negative_monomial = 0;
    std::size_t nonzero_monomial = 0;
    for (std::size_t z_degree = 0;
         z_degree < max_z_degree;
         ++z_degree) {
        std::map<std::array<int, 2>, Rational> monomial;
        for (int x_order = 0; x_order <= grid_degree; ++x_order) {
            const Polynomial x_basis = binomial_basis(x_order);
            for (int y_order = 0; y_order <= grid_degree; ++y_order) {
                const Rational coefficient = forward_coefficient(
                    grid,
                    x_order,
                    y_order,
                    z_degree
                );
                if (coefficient == 0) {
                    continue;
                }
                const Polynomial y_basis = binomial_basis(y_order);
                for (std::size_t x_degree = 0;
                     x_degree < x_basis.coefficients().size();
                     ++x_degree) {
                    for (std::size_t y_degree = 0;
                         y_degree < y_basis.coefficients().size();
                         ++y_degree) {
                        monomial[std::array<int, 2>{
                            static_cast<int>(x_degree),
                            static_cast<int>(y_degree)
                        }] +=
                            coefficient
                            * x_basis.coefficients()[x_degree]
                            * y_basis.coefficients()[y_degree];
                    }
                }
            }
        }
        for (const auto& [exponent, coefficient] : monomial) {
            static_cast<void>(exponent);
            if (coefficient != 0) {
                ++nonzero_monomial;
            }
            if (coefficient < 0) {
                ++negative_monomial;
            }
        }
    }
    std::cout
        << "SU2_ANCHORED_LOBE_FIRST_BLOCK_POWERS"
        << " chamber=" << (high_chamber ? "high" : "low")
        << " power=" << power
        << " denominator_power="
        << grid[0][0].denominator_power
        << " z_degrees=" << max_z_degree
        << " nonzero_binomial_coefficients=" << nonzero
        << " max_x_order=" << max_x_order
        << " max_y_order=" << max_y_order
        << " max_total_order=" << max_total_order
        << " degree_bound_violations=" << degree_bound_violations
        << " negative_binomial_coefficients=" << negative
        << " nonzero_monomial_coefficients=" << nonzero_monomial
        << " negative_monomial_coefficients=" << negative_monomial
        << '\n';
    const int expected_denominator_power = high_chamber ? 1 : 0;
    const std::size_t expected_z_degrees = static_cast<std::size_t>(
        high_chamber ? power + 2 : power + 1
    );
    return negative == 0
        && degree_bound_violations == 0
        && grid[0][0].denominator_power
            == expected_denominator_power
        && max_z_degree == expected_z_degrees;
}

}  // namespace

int main() {
    bool ok = true;
    for (int power = 2; power <= 4; ++power) {
        // The explicit two-spline recurrence has total parameter degree
        // at most 4p+2.  The extra grid layers both recover its Newton
        // coefficients and audit the degree cutoff.
        const int grid_degree = 4 * power + 8;
        ok =
            report_chamber(false, grid_degree, power)
            && ok;
        ok =
            report_chamber(true, grid_degree, power)
            && ok;
    }
    std::cout
        << "SU2_ANCHORED_LOBE_FIRST_BLOCK_POWERS"
        << " result="
        << (ok ? "PASS_EXACT_BINOMIAL_CERTIFICATE" : "FAIL")
        << '\n';
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
