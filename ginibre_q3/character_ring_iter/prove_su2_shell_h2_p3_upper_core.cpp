// Exact bivariate certificate for the separated-upper reflected P_3 core.
//
// In Lemma 5A8H28P3F put d=4a+delta, x=2b+epsilon and restrict to
// 0<=s<=d on one rail s=sigma (mod 4).  The ordinary P_3 profile is
//
//   binom(s+2,2)-binom(s-2x+1,2)_+.
//
// This source sums both quadratic pieces exactly, by Faulhaber identities,
// and proves the resulting polynomial current coefficientwise on the two
// natural cones where the delayed second piece is present or absent.

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

Polynomial power_sum(int exponent, const Polynomial& upper) {
    if (exponent == 0) {
        return upper + constant(1);
    }
    if (exponent == 1) {
        return upper * (upper + constant(1)) * Polynomial(Rational(1, 2));
    }
    if (exponent == 2) {
        return upper * (upper + constant(1))
            * (scale(upper, 2) + constant(1))
            * Polynomial(Rational(1, 6));
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
            throw std::runtime_error("core current has degree above three");
        }
        Polynomial term(coefficient);
        term *= power(first, exponent[0]);
        term *= power(second, exponent[1]);
        term *= power_sum(exponent[2], upper);
        result += term;
    }
    return result;
}

Polynomial rail_sum(
    int first_value,
    const Polynomial& upper,
    const Polynomial& slope,
    const Polynomial& intercept
) {
    const Polynomial index = Polynomial::variable(2);
    const Polynomial s = scale(index, 4) + constant(first_value);
    const Polynomial profile = binomial(s + constant(2), 2);
    return sum_over_index(profile * (slope * s - intercept), upper);
}

bool coefficientwise_nonnegative(const Polynomial& value) {
    return std::all_of(
        value.terms().begin(),
        value.terms().end(),
        [](const auto& term) { return term.second >= 0; }
    );
}

std::optional<int> nonnegative_integer_tail_certificate(
    const Polynomial& value
) {
    if (coefficientwise_nonnegative(value)) {
        return 0;
    }
    const Polynomial first = Polynomial::variable(0);
    const Polynomial second = Polynomial::variable(1);
    for (int threshold = 1; threshold <= 32; ++threshold) {
        bool finite_prefix = true;
        for (int value_of_second = 0;
             value_of_second < threshold;
             ++value_of_second) {
            const Polynomial slice = substitute(
                value,
                std::array<Polynomial, 3>{
                    first,
                    constant(value_of_second),
                    constant(0)
                }
            );
            if (!coefficientwise_nonnegative(slice)) {
                finite_prefix = false;
                break;
            }
        }
        if (!finite_prefix) {
            continue;
        }
        const Polynomial shifted = substitute(
            value,
            std::array<Polynomial, 3>{
                first,
                second + constant(threshold),
                constant(0)
            }
        );
        if (coefficientwise_nonnegative(shifted)) {
            return threshold;
        }
    }
    return std::nullopt;
}

bool nonnegative_sum_floor_certificate(
    const Polynomial& value,
    int minimum_sum
) {
    if (minimum_sum <= 0) {
        return coefficientwise_nonnegative(value);
    }
    const Polynomial first = Polynomial::variable(0);
    const Polynomial second = Polynomial::variable(1);
    for (int fixed_second = 0;
         fixed_second < minimum_sum;
         ++fixed_second) {
        const Polynomial slice = substitute(
            value,
            std::array<Polynomial, 3>{
                first + constant(minimum_sum - fixed_second),
                constant(fixed_second),
                constant(0)
            }
        );
        if (!coefficientwise_nonnegative(slice)) {
            return false;
        }
    }
    return coefficientwise_nonnegative(substitute(
        value,
        std::array<Polynomial, 3>{
            first,
            second + constant(minimum_sum),
            constant(0)
        }
    ));
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

Polynomial replace_ab(
    const Polynomial& value,
    const Polynomial& a,
    const Polynomial& b
) {
    return substitute(value, std::array<Polynomial, 3>{
        a,
        b,
        constant(0)
    });
}

}  // namespace

int main() {
    try {
        const Polynomial a = Polynomial::variable(0);
        const Polynomial b = Polynomial::variable(1);
        const Polynomial u = Polynomial::variable(0);
        const Polynomial v = Polynomial::variable(1);
        std::size_t cases = 0U;
        std::size_t certified = 0U;
        std::size_t tail_certified = 0U;

        for (int delta = 0; delta < 4; ++delta) {
            const int minimum_a = delta == 3 ? 2 : 3;
            for (int epsilon = 0; epsilon < 2; ++epsilon) {
                const Polynomial d = scale(a, 4) + constant(delta);
                const Polynomial x = scale(b, 2) + constant(epsilon);
                const Polynomial slope = scale(d, 6) + scale(x, 4)
                    + constant(3);
                const Polynomial intercept =
                    (scale(d, 2) + constant(1)) * (d + x);
                for (int sigma = 0; sigma < 4; ++sigma) {
                    const int tau =
                        (sigma - 2 * epsilon - 1 + 8) % 4;
                    const int first_count_shift = floor_divide(
                        delta - sigma,
                        4
                    );
                    const int second_count_shift = floor_divide(
                        delta - 2 * epsilon - 1 - tau,
                        4
                    );
                    const Polynomial first_upper =
                        a + constant(first_count_shift);
                    const Polynomial second_upper =
                        a - b + constant(second_count_shift);
                    const Polynomial total = rail_sum(
                        sigma,
                        first_upper,
                        slope,
                        intercept
                    ) - rail_sum(
                        tau,
                        second_upper,
                        slope,
                        intercept - slope * (
                            scale(x, 2) + constant(1)
                        )
                    );

                    // The delayed quadratic term is present exactly when
                    // a-b+second_count_shift>=0.  Put that count equal to
                    // v; this enlarges the physical cone but preserves the
                    // exact polynomial identity on it.
                    const Polynomial present = replace_ab(
                        total,
                        u + v - constant(second_count_shift),
                        u
                    );
                    ++cases;
                    const int present_minimum_sum = std::max(
                        0,
                        minimum_a + second_count_shift
                    );
                    if (!nonnegative_sum_floor_certificate(
                            present,
                            present_minimum_sum
                        )) {
                        std::cout
                            << "SU2_SHELL_H2_P3_UPPER_CORE_FAILURE"
                            << " branch=present"
                            << " delta=" << delta
                            << " epsilon=" << epsilon
                            << " sigma=" << sigma
                            << " negative_coefficients="
                            << negative_coefficients(present) << '\n';
                        print_polynomial(present);
                        return EXIT_FAILURE;
                    }
                    ++certified;

                    // When it is absent, write b=a+r+c+1, where c is the
                    // second count shift.  Here a=u+minimum_a and r>=0.
                    const Polynomial absent = replace_ab(
                        total,
                        u + constant(minimum_a),
                        u + v + constant(
                            minimum_a + second_count_shift + 1
                        )
                    );
                    ++cases;
                    const std::optional<int> absent_certificate =
                        nonnegative_integer_tail_certificate(absent);
                    if (!absent_certificate.has_value()) {
                        std::cout
                            << "SU2_SHELL_H2_P3_UPPER_CORE_FAILURE"
                            << " branch=absent"
                            << " delta=" << delta
                            << " epsilon=" << epsilon
                            << " sigma=" << sigma
                            << " negative_coefficients="
                            << negative_coefficients(absent) << '\n';
                        print_polynomial(absent);
                        return EXIT_FAILURE;
                    }
                    if (*absent_certificate > 0) {
                        ++tail_certified;
                    }
                    ++certified;
                }
            }
        }

        std::cout
            << "SU2_SHELL_H2_P3_UPPER_CORE"
            << " cases=" << cases
            << " certified=" << certified
            << " integer_tail_certificates=" << tail_certified
            << " result=PASS_EXACT_FINITE_CONE_CERTIFICATE\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
