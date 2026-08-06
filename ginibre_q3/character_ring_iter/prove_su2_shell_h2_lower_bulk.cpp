// Exact bulk certificate for the separated-lower H_2 current.
//
// In the chamber y>=5d+8 all terminal paths remain above the lower wall;
// the P_2,P_3,P_4 profiles are ordinary b-fold interval convolutions.  This
// source sums their residue-four moments by truncated-binomial identities
// and proves the resulting lower current coefficientwise on every residue
// cone d=8a+delta, y=8b+epsilon, z=4c+zeta.

#define main prove_su2_k4_intermediate_embedded_main
#include "prove_su2_k4_intermediate.cpp"
#undef main

namespace {

int floor_divide(int numerator, int denominator) {
    if (denominator <= 0) {
        throw std::runtime_error("nonpositive divisor");
    }
    if (numerator >= 0) {
        return numerator / denominator;
    }
    return -((-numerator + denominator - 1) / denominator);
}

int ceil_divide(int numerator, int denominator) {
    return -floor_divide(-numerator, denominator);
}

int positive_modulo(int value, int modulus) {
    return (value % modulus + modulus) % modulus;
}

Polynomial power_sum(int exponent, const Polynomial& upper) {
    if (exponent == 0) {
        return upper + constant(1);
    }
    if (exponent == 1) {
        return upper * (upper + constant(1)) * Polynomial(Rational(1, 2));
    }
    if (exponent == 2) {
        return upper * (upper + constant(1))
            * (scale(upper, 2) + constant(1)) * Polynomial(Rational(1, 6));
    }
    if (exponent == 3) {
        const Polynomial triangular = upper * (upper + constant(1))
            * Polynomial(Rational(1, 2));
        return triangular * triangular;
    }
    throw std::runtime_error("unexpected summation degree");
}

Polynomial sum_over_index(
    const Polynomial& value,
    const Polynomial& upper
) {
    const Polynomial first = Polynomial::variable(0);
    const Polynomial second = Polynomial::variable(1);
    Polynomial result;
    for (const auto& [exponent, coefficient] : value.terms()) {
        if (exponent[2] > 3) {
            throw std::runtime_error("bulk summand degree above three");
        }
        Polynomial term(coefficient);
        term *= power(first, exponent[0]);
        term *= power(second, exponent[1]);
        term *= power_sum(exponent[2], upper);
        result += term;
    }
    return result;
}

Polynomial lift_av(const Polynomial& value) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial v = Polynomial::variable(1);
    return substitute(value, std::array<Polynomial, 3>{
        a,
        v,
        constant(0)
    });
}

bool coefficientwise_nonnegative(const Polynomial& value) {
    return std::all_of(
        value.terms().begin(),
        value.terms().end(),
        [](const auto& term) { return term.second >= 0; }
    );
}

std::size_t negative_coefficients(const Polynomial& value) {
    return static_cast<std::size_t>(std::count_if(
        value.terms().begin(),
        value.terms().end(),
        [](const auto& term) { return term.second < 0; }
    ));
}

void print_polynomial(const Polynomial& value) {
    for (const auto& [exponent, coefficient] : value.terms()) {
        std::cout
            << " term=(" << exponent[0] << ',' << exponent[1] << ','
            << exponent[2] << ") coefficient=" << coefficient << '\n';
    }
}

// Sum W^moment (N_d^power)_(W,y) on the fixed W mod 8 rail.  In the
// bulk chamber W=y+power*d-2t stays positive, and the profile is the
// b-fold interval convolution indexed by 0<=t<=power*d.
Polynomial bulk_moment(
    int power_value,
    int d_residue,
    int y_residue,
    int depth_residue,
    int y_offset,
    int moment
) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial v = Polynomial::variable(1);
    const Polynomial index = Polynomial::variable(2);
    const Polynomial d = scale(a, 8) + constant(d_residue);
    const Polynomial y = scale(
        scale(a, 5) + constant(y_offset) + v,
        8
    ) + constant(y_residue);
    Polynomial result;
    const int upper_constant = floor_divide(
        power_value * d_residue - depth_residue,
        4
    );
    for (int image = 0; image < power_value; ++image) {
        const int lower_constant = ceil_divide(
            image * (d_residue + 1) - depth_residue,
            4
        );
        const Polynomial lower = scale(a, 2 * image)
            + constant(lower_constant);
        const Polynomial upper = scale(a, 2 * power_value)
            + constant(upper_constant);
        const Polynomial depth = scale(index + lower, 4)
            + constant(depth_residue);
        const Polynomial profile = binomial(
            depth - scale(d + constant(1), image)
                + constant(power_value - 1),
            power_value - 1
        );
        const Polynomial label = y + scale(d, power_value) - scale(depth, 2);
        const long coefficient = binomial_long(power_value, image)
            * (image % 2 == 0 ? 1 : -1);
        result += scale(sum_over_index(
            profile * power(label, moment),
            upper - lower
        ), coefficient);
    }
    return lift_av(result);
}

Polynomial bulk_current(
    int d_residue,
    int y_residue,
    int z_residue,
    int rail,
    int y_offset
) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial v = Polynomial::variable(1);
    const Polynomial c = Polynomial::variable(2);
    const Polynomial d = scale(a, 8) + constant(d_residue);
    const Polynomial y = scale(
        scale(a, 5) + constant(y_offset) + v,
        8
    ) + constant(y_residue);
    const Polynomial q = y + scale(d, 3) + scale(c, 4)
        + constant(z_residue + 1);
    const int q_residue = positive_modulo(
        y_residue + 3 * d_residue + z_residue + 1,
        4
    );
    const int p2_depth = positive_modulo(
        (y_residue + d_residue) / 2 - q_residue + rail,
        4
    );
    const int p4_depth = positive_modulo(
        (y_residue + 3 * d_residue) / 2 - q_residue + rail,
        4
    );
    const int p3_depth = positive_modulo(
        (y_residue + 3 * d_residue) / 2 - rail,
        4
    );
    const Polynomial p20 = bulk_moment(
        2, d_residue, y_residue, p2_depth, y_offset, 0
    );
    const Polynomial p21 = bulk_moment(
        2, d_residue, y_residue, p2_depth, y_offset, 1
    );
    const Polynomial p22 = bulk_moment(
        2, d_residue, y_residue, p2_depth, y_offset, 2
    );
    const Polynomial p40 = bulk_moment(
        4, d_residue, y_residue, p4_depth, y_offset, 0
    );
    const Polynomial p30 = bulk_moment(
        3, d_residue, y_residue, p3_depth, y_offset, 0
    );
    const Polynomial p31 = bulk_moment(
        3, d_residue, y_residue, p3_depth, y_offset, 1
    );

    // In this chamber W>d and every P_3 label is in the fixed adverse
    // suffix.  These are the moment forms of (P5A.102BUK36).
    Polynomial result = (scale(p21, 4) + scale(p20, 4)) * q * q;
    result += (
        (scale(d * d, 2) + scale(d, 4) + constant(4)) * p20
        - scale(p22, 2)
    ) * q;
    result -= p22;
    result -= scale(p21, 4) * (d + constant(1)) * (d + constant(1));
    result -= (scale(d * d, 3) + scale(d, 6) + constant(2)) * p20;
    result += scale(p40, 4) * q + scale(p40, 2);
    result -= scale(p30, 4) * q * q;
    result += (
        scale(p31, 4) - scale(p30, 8) * d - scale(p30, 8)
    ) * q;
    result += scale(p31, 2)
        + (scale(d * d, 4) + scale(d, 4)) * p30;
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 6 && std::string(argv[1]) == "--evaluate") {
            const int d_value = std::atoi(argv[2]);
            const int y_value = std::atoi(argv[3]);
            const int z_value = std::atoi(argv[4]);
            const int rail = std::atoi(argv[5]);
            if (
                d_value < 11
                || y_value < 5 * d_value + 8
                || z_value < 0
                || rail < 0
                || rail >= 4
                || (d_value - y_value) % 2 != 0
            ) {
                throw std::runtime_error("input is outside the bulk cone");
            }
            const int delta = positive_modulo(d_value, 8);
            const int epsilon = positive_modulo(y_value, 8);
            const int zeta = positive_modulo(z_value, 4);
            const int y_offset = ceil_divide(
                5 * delta + 8 - epsilon,
                8
            );
            const int a_value = (d_value - delta) / 8;
            const int b_value = (y_value - epsilon) / 8;
            const int v_value = b_value - 5 * a_value - y_offset;
            const int c_value = (z_value - zeta) / 4;
            if (v_value < 0) {
                throw std::runtime_error("bulk coordinate conversion failed");
            }
            const Rational value = evaluate(
                bulk_current(delta, epsilon, zeta, rail, y_offset),
                std::array<int, 3>{a_value, v_value, c_value}
            );
            if (value.denominator() != 1) {
                throw std::runtime_error("nonintegral bulk evaluation");
            }
            std::cout
                << "SU2_SHELL_H2_LOWER_BULK_EVALUATION"
                << " d=" << d_value
                << " y=" << y_value
                << " z=" << z_value
                << " rail=" << rail
                << " value=" << value.numerator() << '\n';
            return EXIT_SUCCESS;
        }
        const Polynomial a = Polynomial::variable(0);
        const Polynomial v = Polynomial::variable(1);
        const Polynomial c = Polynomial::variable(2);
        std::size_t cases = 0U;
        std::size_t certified = 0U;
        for (int delta = 0; delta < 8; ++delta) {
            const int minimum_a = delta <= 2 ? 2 : 1;
            for (int epsilon = delta % 2; epsilon < 8; epsilon += 2) {
                const int y_offset = ceil_divide(
                    5 * delta + 8 - epsilon,
                    8
                );
                for (int zeta = 0; zeta < 4; ++zeta) {
                    for (int rail = 0; rail < 4; ++rail) {
                        const Polynomial current = bulk_current(
                            delta, epsilon, zeta, rail, y_offset
                        );
                        const Polynomial shifted = substitute(
                            current,
                            std::array<Polynomial, 3>{
                                a + constant(minimum_a),
                                v,
                                c
                            }
                        );
                        ++cases;
                        if (!coefficientwise_nonnegative(shifted)) {
                            std::cout
                                << "SU2_SHELL_H2_LOWER_BULK_FAILURE"
                                << " delta=" << delta
                                << " epsilon=" << epsilon
                                << " zeta=" << zeta
                                << " rail=" << rail
                                << " negative_coefficients="
                                << negative_coefficients(shifted) << '\n';
                            print_polynomial(shifted);
                            return EXIT_FAILURE;
                        }
                        ++certified;
                    }
                }
            }
        }
        std::cout
            << "SU2_SHELL_H2_LOWER_BULK"
            << " cases=" << cases
            << " certified=" << certified
            << " result=PASS_EXACT_FINITE_CONE_CERTIFICATE\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
