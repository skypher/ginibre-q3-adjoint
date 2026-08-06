// Exact certificates for y=3d+t, t=0,2,...,10, d>=12.

#define main prove_su2_k4_intermediate_embedded_main
#include "prove_su2_k4_intermediate.cpp"
#undef main

namespace {

int floor_divide(int numerator, int denominator) {
    if (denominator <= 0) throw std::runtime_error("nonpositive divisor");
    if (numerator >= 0) return numerator / denominator;
    return -((-numerator + denominator - 1) / denominator);
}

int ceil_divide(int numerator, int denominator) {
    return -floor_divide(-numerator, denominator);
}

int positive_modulo(int value, int modulus) {
    return (value % modulus + modulus) % modulus;
}

Polynomial power_sum(int exponent, const Polynomial& upper) {
    if (exponent == 0) return upper + constant(1);
    if (exponent == 1) {
        return upper * (upper + constant(1)) * Polynomial(Rational(1, 2));
    }
    if (exponent == 2) {
        return upper * (upper + constant(1))
            * (scale(upper, 2) + constant(1)) * Polynomial(Rational(1, 6));
    }
    if (exponent == 3) {
        const Polynomial triangle = upper * (upper + constant(1))
            * Polynomial(Rational(1, 2));
        return triangle * triangle;
    }
    throw std::runtime_error("unexpected summation degree");
}

Polynomial sum_over_index(
    const Polynomial& value,
    const Polynomial& upper
) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial c = Polynomial::variable(1);
    Polynomial result;
    for (const auto& [exponent, coefficient] : value.terms()) {
        if (exponent[2] > 3) {
            throw std::runtime_error("mid-lip summand degree above three");
        }
        Polynomial term(coefficient);
        term *= power(a, exponent[0]);
        term *= power(c, exponent[1]);
        term *= power_sum(exponent[2], upper);
        result += term;
    }
    return result;
}

Polynomial lift_ac(const Polynomial& value) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial c = Polynomial::variable(1);
    return substitute(value, std::array<Polynomial, 3>{
        a, c, constant(0)
    });
}

bool coefficientwise_nonnegative(const Polynomial& value) {
    return std::all_of(
        value.terms().begin(), value.terms().end(),
        [](const auto& term) { return term.second >= 0; }
    );
}

std::size_t negative_coefficients(const Polynomial& value) {
    return static_cast<std::size_t>(std::count_if(
        value.terms().begin(), value.terms().end(),
        [](const auto& term) { return term.second < 0; }
    ));
}

void print_polynomial(const Polynomial& value) {
    for (const auto& [exponent, coefficient] : value.terms()) {
        std::cout << " term=(" << exponent[0] << ',' << exponent[1]
                  << ',' << exponent[2] << ") coefficient="
                  << coefficient << '\n';
    }
}

Polynomial raw_moment(
    int power_value,
    int delta,
    int t,
    int rail,
    int depth_residue,
    int moment,
    int mode
) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial index = Polynomial::variable(2);
    const Polynomial d = scale(a, 8) + constant(delta);
    const Polynomial y = scale(d, 3) + constant(t);
    Polynomial upper;
    if (mode == 0) {
        upper = scale(a, 2 * power_value) + constant(floor_divide(
            power_value * delta - depth_residue, 4
        ));
    } else if (mode == 1) {
        const int rho_shift = 1 + positive_modulo(rail - delta - 1, 4);
        upper = scale(a, 4) + constant(floor_divide(
            2 * delta + t / 2 - rho_shift - depth_residue, 4
        ));
    } else {
        upper = scale(a, 7) + constant(floor_divide(
            7 * delta + t - 2 * depth_residue, 8
        ));
    }
    Polynomial result;
    for (int image = 0; image < power_value; ++image) {
        const Polynomial lower = scale(a, 2 * image) + constant(
            ceil_divide(image * (delta + 1) - depth_residue, 4)
        );
        const Polynomial range = upper - lower;
        // For the P_3 j=2 image, an endpoint can be strictly before its
        // first congruent depth.  It is then identically absent on this ray.
        if (mode == 1 && image == 2) {
            const Rational endpoint = evaluate(range, {0, 0, 0});
            if (endpoint < 0) continue;
        }
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
            profile * power(label, moment), range
        ), coefficient);
    }
    return lift_ac(result);
}

Polynomial p4_lower_wall_mass(int delta, int t, int depth_residue) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial index = Polynomial::variable(2);
    const Polynomial d = scale(a, 8) + constant(delta);
    const Polynomial y = scale(d, 3) + constant(t);
    const Polynomial upper = scale(a, 7) + constant(floor_divide(
        7 * delta + t - 2 * depth_residue, 8
    ));
    const Polynomial lower = scale(a, 6) + constant(ceil_divide(
        3 * delta + t + 1 - depth_residue, 4
    ));
    const Polynomial depth = scale(index + lower, 4) + constant(depth_residue);
    return lift_ac(sum_over_index(
        binomial(depth - y + constant(2), 3), upper - lower
    ));
}

Polynomial mid_lip_current(int delta, int t, int zeta, int rail) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial c = Polynomial::variable(1);
    const Polynomial d = scale(a, 8) + constant(delta);
    const Polynomial y = scale(d, 3) + constant(t);
    const Polynomial q = y + scale(d, 3) + scale(c, 4)
        + constant(zeta + 1);
    const int y_residue = positive_modulo(3 * delta + t, 8);
    const int q_residue = positive_modulo(
        y_residue + 3 * delta + zeta + 1, 4
    );
    const int p2_depth = positive_modulo(
        (y_residue + delta) / 2 - q_residue + rail, 4
    );
    const int p4_depth = positive_modulo(
        (y_residue + 3 * delta) / 2 - q_residue + rail, 4
    );
    const int p3_depth = positive_modulo(
        (y_residue + 3 * delta) / 2 - rail, 4
    );
    const Polynomial p20 = raw_moment(2, delta, t, rail, p2_depth, 0, 0);
    const Polynomial p21 = raw_moment(2, delta, t, rail, p2_depth, 1, 0);
    const Polynomial p22 = raw_moment(2, delta, t, rail, p2_depth, 2, 0);
    const Polynomial p40 = raw_moment(4, delta, t, rail, p4_depth, 0, 2)
        - p4_lower_wall_mass(delta, t, p4_depth);
    const Polynomial p30 = raw_moment(3, delta, t, rail, p3_depth, 0, 1);
    const Polynomial p31 = raw_moment(3, delta, t, rail, p3_depth, 1, 1);
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
    result += (scale(p31, 4) - scale(p30, 8) * d - scale(p30, 8)) * q;
    result += scale(p31, 2) + (scale(d * d, 4) + scale(d, 4)) * p30;
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 6 && std::string(argv[1]) == "--evaluate") {
            const int d_value = std::atoi(argv[2]);
            const int t = std::atoi(argv[3]);
            const int z_value = std::atoi(argv[4]);
            const int rail = std::atoi(argv[5]);
            if (
                d_value < 12
                || t < 0
                || t > 10
                || t % 2 != 0
                || t >= d_value
                || z_value < 0
                || rail < 0
                || rail >= 4
            ) {
                throw std::runtime_error("input is outside the mid-lip cone");
            }
            const int delta = positive_modulo(d_value, 8);
            const int a_value = (d_value - delta) / 8;
            const int zeta = positive_modulo(z_value, 4);
            const int c_value = (z_value - zeta) / 4;
            const Rational value = evaluate(
                mid_lip_current(delta, t, zeta, rail),
                std::array<int, 3>{a_value, c_value, 0}
            );
            if (value.denominator() != 1) {
                throw std::runtime_error("nonintegral mid-lip evaluation");
            }
            std::cout
                << "SU2_SHELL_H2_LOWER_MID_LIPS_EVALUATION"
                << " d=" << d_value
                << " t=" << t
                << " z=" << z_value
                << " rail=" << rail
                << " value=" << value.numerator() << '\n';
            return EXIT_SUCCESS;
        }
        const Polynomial a = Polynomial::variable(0);
        const Polynomial c = Polynomial::variable(1);
        std::size_t cases = 0U;
        std::size_t certified = 0U;
        for (int delta = 0; delta < 8; ++delta) {
            const int minimum_a = delta <= 3 ? 2 : 1;
            for (int t : {0, 2, 4, 6, 8, 10}) {
                for (int zeta = 0; zeta < 4; ++zeta) {
                    for (int rail = 0; rail < 4; ++rail) {
                        const Polynomial current = mid_lip_current(
                            delta, t, zeta, rail
                        );
                        const Polynomial shifted = substitute(
                            current,
                            std::array<Polynomial, 3>{
                                a + constant(minimum_a), c, constant(0)
                            }
                        );
                        ++cases;
                        if (!coefficientwise_nonnegative(shifted)) {
                            std::cout
                                << "SU2_SHELL_H2_LOWER_MID_LIPS_FAILURE"
                                << " delta=" << delta
                                << " t=" << t
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
            << "SU2_SHELL_H2_LOWER_MID_LIPS"
            << " cases=" << cases
            << " certified=" << certified
            << " result=PASS_EXACT_FINITE_CONE_CERTIFICATE\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
