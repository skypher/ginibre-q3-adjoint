// Exact certificate for 2d<=y<=3d in the separated-lower H_2 current.

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
    const Polynomial b = Polynomial::variable(1);
    Polynomial result;
    for (const auto& [exponent, coefficient] : value.terms()) {
        if (exponent[2] > 3) {
            throw std::runtime_error("mid-slope summand degree above three");
        }
        Polynomial term(coefficient);
        term *= power(a, exponent[0]);
        term *= power(b, exponent[1]);
        term *= power_sum(exponent[2], upper);
        result += term;
    }
    return result;
}

Polynomial lift_ab(const Polynomial& value) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial b = Polynomial::variable(1);
    return substitute(value, std::array<Polynomial, 3>{
        a, b, constant(0)
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

bool nonnegative_sum_floor_certificate(
    const Polynomial& value,
    int minimum_sum
) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial b = Polynomial::variable(1);
    const Polynomial c = Polynomial::variable(2);
    if (minimum_sum <= 0) return coefficientwise_nonnegative(value);
    for (int fixed_b = 0; fixed_b < minimum_sum; ++fixed_b) {
        if (!coefficientwise_nonnegative(substitute(
                value,
                std::array<Polynomial, 3>{
                    a + constant(minimum_sum - fixed_b),
                    constant(fixed_b), c
                }
            ))) {
            return false;
        }
    }
    return coefficientwise_nonnegative(substitute(
        value,
        std::array<Polynomial, 3>{a, b + constant(minimum_sum), c}
    ));
}

struct Coordinates {
    int delta;
    int eta;
    int carry;
    int rho_shift;
    Polynomial d;
    Polynomial s;
    Polynomial y;
};

Coordinates coordinates(int delta, int eta, int rail) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial b = Polynomial::variable(1);
    const int lambda = positive_modulo(delta - eta, 8);
    const int carry = (eta + lambda - delta) / 8;
    const Polynomial d = scale(a + b + constant(carry), 8)
        + constant(delta);
    const Polynomial s = scale(a, 8) + constant(eta);
    const Polynomial y = scale(d, 3) - s;
    const int rho_shift = 1 + positive_modulo(rail - delta - 1, 4);
    return Coordinates{delta, eta, carry, rho_shift, d, s, y};
}

// mode 0: full raw P_2; mode 1: P_3 adverse suffix (j=0,1);
// mode 2: ordinary P_4 through W=0.
Polynomial raw_moment(
    int power_value,
    int depth_residue,
    int moment,
    int mode,
    const Coordinates& data
) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial b = Polynomial::variable(1);
    const Polynomial index = Polynomial::variable(2);
    Polynomial upper;
    if (mode == 0) {
        upper = scale(a + b + constant(data.carry), 2 * power_value)
            + constant(floor_divide(
                power_value * data.delta - depth_residue, 4
            ));
    } else if (mode == 1) {
        // (y+3d-2rho)/2=2d-s/2-rho_shift.
        upper = scale(a, 3) + scale(b, 4)
            + constant(4 * data.carry + floor_divide(
                2 * data.delta - data.eta / 2 - data.rho_shift
                    - depth_residue,
                4
            ));
    } else {
        // (y+4d)/2=(7d-s)/2.
        upper = scale(a, 6) + scale(b, 7)
            + constant(7 * data.carry + floor_divide(
                7 * data.delta - data.eta - 2 * depth_residue,
                8
            ));
    }
    const int images = mode == 1 ? 2 : power_value;
    Polynomial result;
    for (int image = 0; image < images; ++image) {
        const Polynomial lower = scale(
            a + b + constant(data.carry), 2 * image
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
    return lift_ab(result);
}

// Sum P_2(d-W) on the W<=d part of its rail; it corrects g=2W to g=d+W.
Polynomial p2_low_gap(
    int depth_residue,
    const Coordinates& data
) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial b = Polynomial::variable(1);
    const Polynomial index = Polynomial::variable(2);
    const Polynomial upper = scale(a + b + constant(data.carry), 4)
        + constant(floor_divide(2 * data.delta - depth_residue, 4));
    // r>=(4d-s)/2=2d-s/2.
    const Polynomial base = scale(a, 3) + scale(b, 4)
        + constant(4 * data.carry + ceil_divide(
            2 * data.delta - data.eta / 2 - depth_residue, 4
        ));
    Polynomial result;
    for (int image = 0; image < 2; ++image) {
        // On 2d<=y, the cutoff W<=d lies beyond both image thresholds.
        const Polynomial depth = scale(index + base, 4)
            + constant(depth_residue);
        const Polynomial profile = binomial(
            depth - scale(data.d + constant(1), image) + constant(1), 1
        );
        const Polynomial label = data.y + scale(data.d, 2) - scale(depth, 2);
        const long coefficient = image == 0 ? 1 : -2;
        result += scale(sum_over_index(
            profile * (data.d - label), upper - base
        ), coefficient);
    }
    return lift_ab(result);
}

Polynomial p4_lower_wall_mass(
    int depth_residue,
    const Coordinates& data
) {
    const Polynomial a = Polynomial::variable(0);
    const Polynomial b = Polynomial::variable(1);
    const Polynomial index = Polynomial::variable(2);
    const Polynomial upper = scale(a, 6) + scale(b, 7)
        + constant(7 * data.carry + floor_divide(
            7 * data.delta - data.eta - 2 * depth_residue, 8
        ));
    const Polynomial lower = scale(a, 4) + scale(b, 6)
        + constant(6 * data.carry + ceil_divide(
            3 * data.delta - data.eta + 1 - depth_residue, 4
        ));
    const Polynomial depth = scale(index + lower, 4) + constant(depth_residue);
    return lift_ab(sum_over_index(
        binomial(depth - data.y + constant(2), 3), upper - lower
    ));
}

Polynomial mid_current(int delta, int eta, int zeta, int rail) {
    const Polynomial c = Polynomial::variable(2);
    const Coordinates data = coordinates(delta, eta, rail);
    const int y_residue = positive_modulo(3 * delta - eta, 8);
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
    const Polynomial gap = p2_low_gap(p2_depth, data);
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
    // g=d+W rather than 2W on W<=d.
    result += gap * (
        scale(q * q, 2) + scale(q, 2) * (data.d + constant(2))
        - (scale(data.d * data.d, 2) + scale(data.d, 3))
    );
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
            const int s_value = std::atoi(argv[3]);
            const int z_value = std::atoi(argv[4]);
            const int rail = std::atoi(argv[5]);
            if (
                d_value < 11
                || s_value < 0
                || s_value > d_value
                || s_value % 2 != 0
                || z_value < 0
                || rail < 0
                || rail >= 4
            ) {
                throw std::runtime_error("input is outside the mid-slope cone");
            }
            const int delta = positive_modulo(d_value, 8);
            const int eta = positive_modulo(s_value, 8);
            const int lambda = positive_modulo(delta - eta, 8);
            const int carry = (eta + lambda - delta) / 8;
            const int a_value = (s_value - eta) / 8;
            const int b_value = (d_value - delta) / 8 - a_value - carry;
            const int zeta = positive_modulo(z_value, 4);
            const int c_value = (z_value - zeta) / 4;
            if (b_value < 0) {
                throw std::runtime_error("mid-slope coordinate conversion failed");
            }
            const Rational value = evaluate(
                mid_current(delta, eta, zeta, rail),
                std::array<int, 3>{a_value, b_value, c_value}
            );
            if (value.denominator() != 1) {
                throw std::runtime_error("nonintegral mid-slope evaluation");
            }
            std::cout
                << "SU2_SHELL_H2_LOWER_MID_EVALUATION"
                << " d=" << d_value
                << " s=" << s_value
                << " z=" << z_value
                << " rail=" << rail
                << " value=" << value.numerator() << '\n';
            return EXIT_SUCCESS;
        }
        std::size_t cases = 0U;
        std::size_t certified = 0U;
        for (int delta = 0; delta < 8; ++delta) {
            const int minimum_a_plus_b = delta <= 2 ? 2 : 1;
            for (int eta = 0; eta < 8; eta += 2) {
                const int lambda = positive_modulo(delta - eta, 8);
                const int carry = (eta + lambda - delta) / 8;
                const int minimum_sum = std::max(
                    0, minimum_a_plus_b - carry
                );
                for (int zeta = 0; zeta < 4; ++zeta) {
                    for (int rail = 0; rail < 4; ++rail) {
                        const Polynomial current = mid_current(
                            delta, eta, zeta, rail
                        );
                        ++cases;
                        if (!nonnegative_sum_floor_certificate(
                                current, minimum_sum
                            )) {
                            std::cout
                                << "SU2_SHELL_H2_LOWER_MID_FAILURE"
                                << " delta=" << delta
                                << " eta=" << eta
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
            << "SU2_SHELL_H2_LOWER_MID"
            << " cases=" << cases
            << " certified=" << certified
            << " result=PASS_EXACT_FINITE_CONE_CERTIFICATE\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
