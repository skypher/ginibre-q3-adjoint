// Exact certificate for 3d+12<=y<4d in the separated-lower H_2 current.
// Here P_2 and P_3 are ordinary.  P_4 has one lower-wall image, which is
// summed explicitly together with the truncated ordinary range.

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
    const Polynomial v = Polynomial::variable(1);
    Polynomial result;
    for (const auto& [exponent, coefficient] : value.terms()) {
        if (exponent[2] > 3) {
            throw std::runtime_error("mid-high summand degree above three");
        }
        Polynomial term(coefficient);
        term *= power(a, exponent[0]);
        term *= power(v, exponent[1]);
        term *= power_sum(exponent[2], upper);
        result += term;
    }
    return result;
}

Polynomial lift_av(const Polynomial& value) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial v = Polynomial::variable(1);
    return substitute(value, std::array<Polynomial, 3>{
        a, v, constant(0)
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

struct Coordinates {
    int delta;
    int lambda;
    int eta;
    int carry;
    int b_offset;
    int rho_shift;
    Polynomial d;
    Polynomial y;
    Polynomial t;
};

Coordinates coordinates(
    int delta,
    int lambda,
    int b_offset,
    int rail
) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial v = Polynomial::variable(1);
    const int eta = positive_modulo(delta - lambda, 8);
    const int carry = (eta + lambda - delta) / 8;
    const Polynomial t = scale(v + constant(b_offset), 8)
        + constant(lambda);
    const Polynomial d = scale(a + v + constant(carry + b_offset), 8)
        + constant(delta);
    const Polynomial y = scale(d, 3) + t;
    const int rho_shift = 1 + positive_modulo(rail - delta - 1, 4);
    return Coordinates{
        delta, lambda, eta, carry, b_offset, rho_shift, d, y, t
    };
}

// Raw P_b moment on a fixed depth rail.  P_2 is full; P_3 ends at its
// adverse suffix endpoint; the ordinary contribution to P_4 ends at W=0.
Polynomial raw_moment(
    int power_value,
    int depth_residue,
    int moment,
    int mode,
    const Coordinates& data
) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial v = Polynomial::variable(1);
    const Polynomial index = Polynomial::variable(2);
    Polynomial upper;
    if (mode == 0) {
        upper = scale(
            a + v + constant(data.b_offset + data.carry),
            2 * power_value
        )
            + constant(floor_divide(
                power_value * data.delta - depth_residue, 4
            ));
    } else if (mode == 1) {
        // (y+3d-2rho)/2=2d+t/2-rho_shift.
        upper = scale(a, 4) + scale(v + constant(data.b_offset), 5)
            + constant(4 * data.carry)
            + constant(floor_divide(
                2 * data.delta + data.lambda / 2 - data.rho_shift
                    - depth_residue,
                4
            ));
    } else {
        // (y+4d)/2=(7d+t)/2, the W=0 endpoint for P_4.
        upper = scale(a, 7) + scale(v + constant(data.b_offset), 8)
            + constant(7 * data.carry)
            + constant(floor_divide(
                7 * data.delta + data.lambda - depth_residue * 2,
                8
            ));
    }
    Polynomial result;
    for (int image = 0; image < power_value; ++image) {
        const Polynomial lower = scale(
            a + v + constant(data.b_offset + data.carry), 2 * image
        ) + constant(ceil_divide(
            image * (data.delta + 1) - depth_residue, 4
        ));
        const Polynomial depth = scale(index + lower, 4)
            + constant(depth_residue);
        const Polynomial profile = binomial(
            depth - scale(data.d + constant(1), image)
                + constant(power_value - 1),
            power_value - 1
        );
        const Polynomial label = data.y + scale(data.d, power_value)
            - scale(depth, 2);
        const long coefficient = binomial_long(power_value, image)
            * (image % 2 == 0 ? 1 : -1);
        result += scale(sum_over_index(
            profile * power(label, moment), upper - lower
        ), coefficient);
    }
    return lift_av(result);
}

// The sole active lower-wall contribution to P_4 is
// -binom(t-y+2,3) for t>=y+1.
Polynomial p4_lower_wall_mass(
    int depth_residue,
    const Coordinates& data
) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial v = Polynomial::variable(1);
    const Polynomial index = Polynomial::variable(2);
    const Polynomial upper = scale(a, 7)
        + scale(v + constant(data.b_offset), 8)
        + constant(7 * data.carry)
        + constant(floor_divide(
            7 * data.delta + data.lambda - depth_residue * 2,
            8
        ));
    const Polynomial lower = scale(a, 6)
        + scale(v + constant(data.b_offset), 8)
        + constant(6 * data.carry)
        + constant(ceil_divide(
            3 * data.delta + data.lambda + 1 - depth_residue,
            4
        ));
    const Polynomial depth = scale(index + lower, 4)
        + constant(depth_residue);
    return lift_av(sum_over_index(
        binomial(depth - data.y + constant(2), 3), upper - lower
    ));
}

Polynomial mid_high_current(
    int delta,
    int lambda,
    int b_offset,
    int zeta,
    int rail
) {
    const Polynomial c = Polynomial::variable(2);
    const Coordinates data = coordinates(delta, lambda, b_offset, rail);
    const int y_residue = positive_modulo(3 * delta + lambda, 8);
    const Polynomial q = data.y + scale(data.d, 3) + scale(c, 4)
        + constant(zeta + 1);
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
    const Polynomial p20 = raw_moment(2, p2_depth, 0, 0, data);
    const Polynomial p21 = raw_moment(2, p2_depth, 1, 0, data);
    const Polynomial p22 = raw_moment(2, p2_depth, 2, 0, data);
    const Polynomial p40 = raw_moment(4, p4_depth, 0, 2, data)
        - p4_lower_wall_mass(p4_depth, data);
    const Polynomial p30 = raw_moment(3, p3_depth, 0, 1, data);
    const Polynomial p31 = raw_moment(3, p3_depth, 1, 1, data);
    Polynomial result = (scale(p21, 4) + scale(p20, 4)) * q * q;
    result += (
        (scale(data.d * data.d, 2) + scale(data.d, 4) + constant(4))
            * p20 - scale(p22, 2)
    ) * q;
    result -= p22;
    result -= scale(p21, 4)
        * (data.d + constant(1)) * (data.d + constant(1));
    result -= (scale(data.d * data.d, 3) + scale(data.d, 6) + constant(2))
        * p20;
    result += scale(p40, 4) * q + scale(p40, 2);
    result -= scale(p30, 4) * q * q;
    result += (scale(p31, 4) - scale(p30, 8) * data.d - scale(p30, 8)) * q;
    result += scale(p31, 2) + (scale(data.d * data.d, 4) + scale(data.d, 4))
        * p30;
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 6 && std::string(argv[1]) == "--evaluate") {
            const int d_value = std::atoi(argv[2]);
            const int t_value = std::atoi(argv[3]);
            const int z_value = std::atoi(argv[4]);
            const int rail = std::atoi(argv[5]);
            if (
                d_value < 12
                || t_value < 12
                || t_value >= d_value
                || t_value % 2 != 0
                || z_value < 0
                || rail < 0
                || rail >= 4
            ) {
                throw std::runtime_error("input is outside the mid-high cone");
            }
            const int delta = positive_modulo(d_value, 8);
            const int lambda = positive_modulo(t_value, 8);
            const int b_offset = ceil_divide(12 - lambda, 8);
            const int v_value = (t_value - lambda) / 8 - b_offset;
            const int eta = positive_modulo(delta - lambda, 8);
            const int carry = (eta + lambda - delta) / 8;
            const int a_value = (d_value - delta) / 8
                - v_value - b_offset - carry;
            const int zeta = positive_modulo(z_value, 4);
            const int c_value = (z_value - zeta) / 4;
            if (a_value < 0 || v_value < 0) {
                throw std::runtime_error("mid-high coordinate conversion failed");
            }
            const Rational value = evaluate(
                mid_high_current(delta, lambda, b_offset, zeta, rail),
                std::array<int, 3>{a_value, v_value, c_value}
            );
            if (value.denominator() != 1) {
                throw std::runtime_error("nonintegral mid-high evaluation");
            }
            std::cout
                << "SU2_SHELL_H2_LOWER_MID_HIGH_EVALUATION"
                << " d=" << d_value
                << " t=" << t_value
                << " z=" << z_value
                << " rail=" << rail
                << " value=" << value.numerator() << '\n';
            return EXIT_SUCCESS;
        }
        std::size_t cases = 0U;
        std::size_t certified = 0U;
        for (int delta = 0; delta < 8; ++delta) {
            for (int lambda : {0, 2, 4, 6}) {
                const int b_offset = ceil_divide(12 - lambda, 8);
                for (int zeta = 0; zeta < 4; ++zeta) {
                    for (int rail = 0; rail < 4; ++rail) {
                        const Polynomial current = mid_high_current(
                            delta, lambda, b_offset, zeta, rail
                        );
                        ++cases;
                        if (!coefficientwise_nonnegative(current)) {
                            std::cout
                                << "SU2_SHELL_H2_LOWER_MID_HIGH_FAILURE"
                                << " delta=" << delta
                                << " lambda=" << lambda
                                << " zeta=" << zeta
                                << " rail=" << rail
                                << " negative_coefficients="
                                << negative_coefficients(current) << '\n';
                            print_polynomial(current);
                            return EXIT_FAILURE;
                        }
                        ++certified;
                    }
                }
            }
        }
        std::cout
            << "SU2_SHELL_H2_LOWER_MID_HIGH"
            << " cases=" << cases
            << " certified=" << certified
            << " result=PASS_EXACT_FINITE_CONE_CERTIFICATE\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
