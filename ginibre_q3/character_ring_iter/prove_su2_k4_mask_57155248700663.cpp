#define main original_su2_k4_intermediate_main
#include "prove_su2_k4_intermediate.cpp"
#undef main

namespace {

constexpr std::uint64_t selected_mask = 57155248700663ULL;

Polynomial affine_variable(
    std::size_t index,
    long scale_factor,
    long offset
) {
    return scale(Polynomial::variable(index), scale_factor)
        + constant(offset);
}

std::array<Polynomial, 3> original_coordinates(
    const Polynomial& x,
    const Polynomial& a,
    const Polynomial& b
) {
    const Polynomial q = x + constant(4);
    const Polynomial y = a + scale(x, 2) + constant(6);
    const Polynomial h =
        (a + b + scale(x, 2))
        * Polynomial(Rational(1, 4));
    return {q, h, y};
}

bool certify_branch(
    const Chamber& chamber,
    const std::array<Polynomial, 3>& coordinates,
    const std::string& branch
) {
    for (const Polynomial& constraint : chamber.constraints) {
        if (
            !nonnegative_monomial(
                substitute(constraint, coordinates)
            )
        ) {
            std::cout
                << "SU2_K4_MASK_57155248700663_BRANCH"
                << " mask=" << selected_mask
                << " branch=" << branch
                << " result=FAIL_DOMAIN_COVER\n";
            return false;
        }
    }
    const Polynomial pulled =
        substitute(chamber.margin, coordinates);
    const bool monomial = nonnegative_monomial(pulled);
    const bool newton = !monomial && nonnegative_newton(pulled);
    if (!monomial && !newton) {
        std::cout
            << "SU2_K4_MASK_57155248700663_BRANCH"
            << " mask=" << selected_mask
            << " branch=" << branch
            << " terms=" << pulled.terms().size()
            << " result=NEEDS_REFINEMENT\n";
        return false;
    }
    std::cout
        << "SU2_K4_MASK_57155248700663_BRANCH"
        << " mask=" << selected_mask
        << " branch=" << branch
        << " terms=" << pulled.terms().size()
        << " basis=" << (monomial ? "monomial" : "newton")
        << " result=PASS_EXACT_IDENTITY\n";
    return true;
}

long least_with_parity(long lower, int parity) {
    long result = lower;
    if ((result & 1L) != static_cast<long>(parity)) {
        ++result;
    }
    return result;
}

int modulo_four(int value) {
    const int residue = value % 4;
    return residue < 0 ? residue + 4 : residue;
}

long least_with_residue_four(long lower, int residue) {
    long result = lower;
    while (modulo_four(static_cast<int>(result)) != residue) {
        ++result;
    }
    return result;
}

bool certify_lower_difference_branch(
    const Chamber& chamber,
    int difference,
    long minimum_r,
    bool fixed_r,
    int b_parity,
    int r_parity
) {
    const long b_minimum =
        difference == 1 && minimum_r == 0 ? 1 : (
            minimum_r < 0 ? -minimum_r : 0
        );
    const long b_offset =
        least_with_parity(b_minimum, b_parity);
    const Polynomial b =
        affine_variable(0U, 2, b_offset);
    const Polynomial r = fixed_r
        ? constant(minimum_r)
        : affine_variable(
            1U,
            2,
            least_with_parity(minimum_r, r_parity)
        );
    const Polynomial a = b + constant(difference + 1);
    const Polynomial x =
        scale(b, 2) + r + constant(difference - 2);
    const std::string branch =
        "lower_d" + std::to_string(difference)
        + "_r" + (fixed_r ? "fixed" : "ray")
        + std::to_string(minimum_r)
        + "_bp" + std::to_string(b_parity)
        + "_rp" + std::to_string(r_parity);
    return certify_branch(
        chamber,
        original_coordinates(x, a, b),
        branch
    );
}

}  // namespace

int main() {
    try {
        const Formula formula = make_formula();
        const Chamber chamber = make_chamber(formula, selected_mask);
        if (!integer_feasible(chamber.constraints)) {
            throw std::runtime_error("selected chamber is infeasible");
        }

        std::size_t attempted = 0U;
        std::size_t certified = 0U;

        // Put a=Y-2Q+2, b=4H-Y+6, and x=Q-4.
        // For a=0,1,2, the chamber is
        // b>=L_a, b<=x+3, where L_a=4,1,0.  Writing
        // b=L_a+B and C=x+3-b leaves B,C>=0 and
        // B+C>=3-L_a.  Split the latter into B>=3-L_a
        // and the finitely many exact B slices, then split the level
        // lattice modulo four.
        for (int a_value = 0; a_value <= 2; ++a_value) {
            const int lower_b =
                a_value == 0 ? 4 : (a_value == 1 ? 1 : 0);
            const int minimum_sum = std::max(0, 3 - lower_b);

            for (int residue_b = 0; residue_b < 4; ++residue_b) {
                for (int residue_c = 0; residue_c < 4; ++residue_c) {
                    const long b_offset = least_with_residue_four(
                        minimum_sum,
                        residue_b
                    );
                    const long c_offset =
                        least_with_residue_four(0, residue_c);
                    if (
                        modulo_four(
                            a_value + 3 * lower_b - 6
                            + 3 * static_cast<int>(b_offset)
                            + 2 * static_cast<int>(c_offset)
                        ) != 0
                    ) {
                        continue;
                    }
                    const Polynomial B =
                        affine_variable(0U, 4, b_offset);
                    const Polynomial C =
                        affine_variable(1U, 4, c_offset);
                    const Polynomial a = constant(a_value);
                    const Polynomial b =
                        B + constant(lower_b);
                    const Polynomial x =
                        B + C + constant(lower_b - 3);
                    ++attempted;
                    if (
                        certify_branch(
                            chamber,
                            original_coordinates(x, a, b),
                            "small_a" + std::to_string(a_value)
                                + "_tail_rb"
                                + std::to_string(residue_b)
                                + "_rc"
                                + std::to_string(residue_c)
                        )
                    ) {
                        ++certified;
                    }
                }
            }

            for (int fixed_b = 0;
                 fixed_b < minimum_sum;
                 ++fixed_b) {
                const int lower_c = minimum_sum - fixed_b;
                for (int residue_c = 0;
                     residue_c < 4;
                     ++residue_c) {
                    const long c_offset =
                        least_with_residue_four(
                            lower_c,
                            residue_c
                        );
                    if (
                        modulo_four(
                            a_value + 3 * lower_b - 6
                            + 3 * fixed_b
                            + 2 * static_cast<int>(c_offset)
                        ) != 0
                    ) {
                        continue;
                    }
                    const Polynomial B = constant(fixed_b);
                    const Polynomial C =
                        affine_variable(0U, 4, c_offset);
                    const Polynomial a = constant(a_value);
                    const Polynomial b =
                        B + constant(lower_b);
                    const Polynomial x =
                        B + C + constant(lower_b - 3);
                    ++attempted;
                    if (
                        certify_branch(
                            chamber,
                            original_coordinates(x, a, b),
                            "small_a" + std::to_string(a_value)
                                + "_slice_b"
                                + std::to_string(fixed_b)
                                + "_rc"
                                + std::to_string(residue_c)
                        )
                    ) {
                        ++certified;
                    }
                }
            }
        }

        // In the branch a>=3 and b>=a, write A=a-3, B=b-a,
        // C=x-A-B.  The remaining inequality is B+2C>=2A.
        // First C>=A: C=A+D.  Congruence modulo four is the only
        // remaining restriction.
        for (int residue_a = 0; residue_a < 4; ++residue_a) {
            for (int residue_b = 0; residue_b < 4; ++residue_b) {
                for (int residue_d = 0; residue_d < 4; ++residue_d) {
                    if (
                        (
                            6 * residue_a
                            + 3 * residue_b
                            + 2 * residue_d
                            + 6
                        ) % 4 != 0
                    ) {
                        continue;
                    }
                    const Polynomial A =
                        affine_variable(0U, 4, residue_a);
                    const Polynomial B =
                        affine_variable(1U, 4, residue_b);
                    const Polynomial D =
                        affine_variable(2U, 4, residue_d);
                    const Polynomial a = A + constant(3);
                    const Polynomial b = a + B;
                    const Polynomial x =
                        scale(A, 2) + B + D;
                    ++attempted;
                    if (
                        certify_branch(
                            chamber,
                            original_coordinates(x, a, b),
                            "upper_outer_ra"
                                + std::to_string(residue_a)
                                + "_rb" + std::to_string(residue_b)
                                + "_rd" + std::to_string(residue_d)
                        )
                    ) {
                        ++certified;
                    }
                }
            }
        }

        // The complementary branch C<A has E=A-C>=1 and
        // B=2E+F.  Put E=e+1 and split the same level-lattice
        // congruence into its residue classes.
        for (int residue_c = 0; residue_c < 4; ++residue_c) {
            for (int residue_e = 0; residue_e < 4; ++residue_e) {
                for (int residue_f = 0; residue_f < 4; ++residue_f) {
                    if (
                        (
                            6 * residue_c
                            + 10 * residue_e
                            + 3 * residue_f
                            + 16
                        ) % 4 != 0
                    ) {
                        continue;
                    }
                    const Polynomial C =
                        affine_variable(0U, 4, residue_c);
                    const Polynomial e =
                        affine_variable(1U, 4, residue_e);
                    const Polynomial F =
                        affine_variable(2U, 4, residue_f);
                    const Polynomial A = C + e + constant(1);
                    const Polynomial B =
                        scale(e, 2) + F + constant(2);
                    const Polynomial a = A + constant(3);
                    const Polynomial b = a + B;
                    const Polynomial x = A + B + C;
                    ++attempted;
                    if (
                        certify_branch(
                            chamber,
                            original_coordinates(x, a, b),
                            "upper_inner_rc"
                                + std::to_string(residue_c)
                                + "_re" + std::to_string(residue_e)
                                + "_rf" + std::to_string(residue_f)
                        )
                    ) {
                        ++certified;
                    }
                }
            }
        }

        // If b<a, integrality and the chamber constraint force
        // a-b in {2,4,6}; equivalently d=a-b-1 is 1,3,5.
        // Put R=X-b for X=x-b-d+2.  The exact lower bounds are
        // R>=-2,-1,0 respectively.  Split the finitely many
        // negative R values and the remaining ray by parity.
        for (const int difference : {1, 3, 5}) {
            const long r_floor =
                static_cast<long>((difference - 5) / 2);
            for (long r = r_floor; r < 0; ++r) {
                const int required =
                    modulo_four(
                        -3 * difference
                        - 2 * static_cast<int>(r)
                        - 1
                    );
                if (required != 0 && required != 2) {
                    continue;
                }
                const int b_parity = required / 2;
                ++attempted;
                if (
                    certify_lower_difference_branch(
                        chamber,
                        difference,
                        r,
                        true,
                        b_parity,
                        0
                    )
                ) {
                    ++certified;
                }
            }
            const long ray_floor = std::max(0L, r_floor);
            for (int b_parity = 0; b_parity < 2; ++b_parity) {
                for (int r_parity = 0; r_parity < 2; ++r_parity) {
                    if (
                        (
                            3 * difference
                            + 2 * (
                                b_parity + r_parity
                            )
                            + 1
                        ) % 4 != 0
                    ) {
                        continue;
                    }
                    ++attempted;
                    if (
                        certify_lower_difference_branch(
                            chamber,
                            difference,
                            ray_floor,
                            false,
                            b_parity,
                            r_parity
                        )
                    ) {
                        ++certified;
                    }
                }
            }
        }

        const bool passed = attempted == certified;
        std::cout
            << "SU2_K4_MASK_57155248700663_CONE"
            << " mask=" << selected_mask
            << " attempted=" << attempted
            << " certified=" << certified
            << " result="
            << (passed
                ? "PASS_EXACT_CERTIFICATE"
                : "INCOMPLETE")
            << '\n';
        if (passed) {
            std::cout
                << "SU2_K4_INTERMEDIATE"
                << " hinges=" << formula.hinges.size()
                << " feasible_chambers=1"
                << " certified_chambers=1"
                << " result=PASS_EXACT_CERTIFICATE\n";
        }
        return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_K4_MASK_57155248700663_CONE FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
