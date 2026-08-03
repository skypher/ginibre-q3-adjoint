#include <set>

#define main prove_su2_k4_embedded_main
#include "prove_su2_k4_intermediate.cpp"
#undef main

namespace {

using WideMask = boost::multiprecision::uint128_t;

struct GroupFormula {
    std::vector<Hinge> hinges;
    std::map<int, std::vector<Term>> u;
    std::map<int, std::vector<Term>> f;
};

void append_multiplicity(
    GroupFormula& formula,
    std::vector<Term>& terms,
    const Polynomial& q,
    const Polynomial& period,
    const Polynomial& endpoint,
    int tensor_power
) {
    const int maximum_branch = (tensor_power - 1) / 2;
    for (int branch = 0; branch <= maximum_branch; ++branch) {
        const bool direct = branch % 2 == 0;
        const int translate =
            direct ? branch / 2 : (branch + 1) / 2;
        const Polynomial image = direct
            ? scale(period, translate) + endpoint
            : scale(period, translate) - constant(1) - endpoint;
        const long affine_sign = direct ? 1L : -1L;
        const int maximum_index = maximum_branch - branch;
        for (int index = 0; index <= maximum_index; ++index) {
            const Polynomial slack =
                scale(q, tensor_power - 2 * index)
                - image
                - constant(index);
            const std::size_t hinge = formula.hinges.size();
            formula.hinges.push_back(Hinge{slack});
            const long ordinary_sign =
                index % 2 == 0 ? 1L : -1L;
            terms.push_back(Term{
                hinge,
                affine_sign
                    * ordinary_sign
                    * binomial_long(tensor_power, index)
            });
        }
    }
}

GroupFormula make_group_formula(const std::string& target) {
    const Polynomial q = Polynomial::variable(0);
    const Polynomial h = Polynomial::variable(1);
    const Polynomial y = Polynomial::variable(2);
    const Polynomial period =
        scale(q, 4) + scale(h, 2) + constant(4);

    std::set<int> f_powers;
    std::set<int> u_powers;
    if (target == "g0") {
        f_powers = {4, 5, 6, 7, 8, 9};
        u_powers = {3, 4, 5, 6, 7, 8};
    } else if (target == "g1") {
        f_powers = {4, 5, 6, 7, 8, 9};
        u_powers = {5, 6, 7, 8, 9, 10};
    } else if (target == "g2") {
        f_powers = {6, 7, 8, 9};
        u_powers = {7, 8, 9, 10};
    } else if (target == "c5") {
        // The next anchored kernel is
        // f_10 u_11-f_11 u_10.  It is the exact finite C_5 target,
        // not a low-group relaxation.
        f_powers = {10, 11};
        u_powers = {10, 11};
    } else {
        throw std::runtime_error("target must be c5, g0, g1, or g2");
    }

    GroupFormula formula;
    for (const int power : u_powers) {
        append_multiplicity(
            formula,
            formula.u[power],
            q,
            period,
            y,
            power
        );
    }
    for (const int power : f_powers) {
        append_multiplicity(
            formula,
            formula.f[power],
            q,
            period,
            constant(0),
            power
        );
    }
    if (formula.hinges.size() > 128U) {
        throw std::runtime_error("group formula exceeds wide mask");
    }
    return formula;
}

bool mask_active(
    const WideMask& mask,
    std::size_t hinge
) {
    return (
        mask & (WideMask(1) << hinge)
    ) != 0;
}

Polynomial selected_multiplicity(
    const GroupFormula& formula,
    const WideMask& mask,
    const std::map<int, std::vector<Term>>& families,
    int power
) {
    const auto found = families.find(power);
    if (found == families.end()) {
        throw std::runtime_error("missing multiplicity family");
    }
    Polynomial result;
    for (const Term& term : found->second) {
        if (mask_active(mask, term.hinge)) {
            result += scale(
                binomial(
                    formula.hinges[term.hinge].slack
                        + constant(power - 2),
                    power - 2
                ),
                term.coefficient
            );
        }
    }
    return result;
}

Polynomial group_margin(
    const GroupFormula& formula,
    const WideMask& mask,
    const std::string& target
) {
    const auto f = [&](int power) {
        return selected_multiplicity(
            formula,
            mask,
            formula.f,
            power
        );
    };
    const auto u = [&](int power) {
        return selected_multiplicity(
            formula,
            mask,
            formula.u,
            power
        );
    };
    const auto block = [&](int scalar_power, int walk_power) {
        return
            f(scalar_power) * u(walk_power)
            - f(scalar_power + 1) * u(walk_power - 1);
    };
    if (target == "g0") {
        return
            scale(block(4, 8), 330)
            + scale(block(6, 6), 462)
            + scale(block(8, 4), 165);
    }
    if (target == "g1") {
        return
            scale(block(4, 10), 385)
            + scale(block(6, 8), 1254)
            + scale(block(8, 6), 1122);
    }
    if (target == "g2") {
        return
            scale(block(6, 10), 2035)
            + scale(block(8, 8), 4026);
    }
    if (target == "c5") {
        return block(10, 11);
    }
    throw std::runtime_error("unknown group target");
}

std::vector<WideMask> feasible_group_masks(
    const GroupFormula& formula
) {
    z3::context context;
    z3::solver solver(context);
    const z3::expr q = context.int_const("Q");
    const z3::expr h = context.int_const("H");
    const z3::expr y = context.int_const("Y");
    const std::array<z3::expr, 3> variables{q, h, y};
    solver.add(q >= 1);
    solver.add(h >= 0);
    solver.add(y >= 0);
    solver.add(y <= 2 * q + 1 + h);

    std::vector<WideMask> masks;
    while (solver.check() == z3::sat) {
        const z3::model model = solver.get_model();
        WideMask mask = 0;
        z3::expr block = context.bool_val(false);
        for (std::size_t index = 0;
             index < formula.hinges.size();
             ++index) {
            const z3::expr active = z3_affine(
                context,
                formula.hinges[index].slack,
                variables
            ) >= 0;
            const bool value =
                model.eval(active, true).bool_value() == Z3_L_TRUE;
            if (value) {
                mask |= WideMask(1) << index;
            }
            block = block || (value ? !active : active);
        }
        masks.push_back(mask);
        solver.add(block);
        if (masks.size() % 1000U == 0U) {
            std::cerr
                << "SU2_T4_GROUP_MASKS"
                << " progress=" << masks.size()
                << std::endl;
        }
    }
    std::sort(masks.begin(), masks.end());
    if (std::unique(masks.begin(), masks.end()) != masks.end()) {
        throw std::runtime_error("duplicate group activation mask");
    }
    return masks;
}

Chamber make_group_chamber(
    const GroupFormula& formula,
    const WideMask& mask,
    const std::string& target,
    std::uint64_t identifier
) {
    const Polynomial q = Polynomial::variable(0);
    const Polynomial h = Polynomial::variable(1);
    const Polynomial y = Polynomial::variable(2);
    Chamber chamber{
        identifier,
        {
            q - constant(1),
            h,
            y,
            scale(q, 2) + constant(1) + h - y
        },
        group_margin(formula, mask, target)
    };
    for (std::size_t index = 0;
         index < formula.hinges.size();
         ++index) {
        const Polynomial& slack = formula.hinges[index].slack;
        chamber.constraints.push_back(
            mask_active(mask, index)
                ? slack
                : constant(-1) - slack
        );
    }
    return chamber;
}

bool bounded_group_integer_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    const std::optional<Rational> q_bound =
        affine_upper_bound(constraints, 0U);
    const std::optional<Rational> h_bound =
        affine_upper_bound(constraints, 1U);
    if (
        !q_bound.has_value()
        || !h_bound.has_value()
        || *q_bound < 1
        || *h_bound < 0
    ) {
        return false;
    }
    const Integer maximum_q_integer =
        q_bound->numerator() / q_bound->denominator();
    const Integer maximum_h_integer =
        h_bound->numerator() / h_bound->denominator();
    if (
        maximum_q_integer > 1000
        || maximum_h_integer > 10000
    ) {
        return false;
    }
    const int maximum_q = maximum_q_integer.convert_to<int>();
    const int maximum_h = maximum_h_integer.convert_to<int>();
    std::uint64_t points = 0U;
    for (int q_value = 1; q_value <= maximum_q; ++q_value) {
        for (int h_value = 0; h_value <= maximum_h; ++h_value) {
            for (int y_value = 0;
                 y_value <= 2 * q_value + 1 + h_value;
                 ++y_value) {
                const std::array<int, 3> values{
                    q_value,
                    h_value,
                    y_value
                };
                const bool in_chamber = std::all_of(
                    chamber.constraints.begin(),
                    chamber.constraints.end(),
                    [&values](const Polynomial& constraint) {
                        return evaluate(constraint, values) >= 0;
                    }
                );
                if (!in_chamber) {
                    continue;
                }
                ++points;
                if (evaluate(chamber.margin, values) < 0) {
                    throw std::runtime_error(
                        "bounded group chamber contains a negative point"
                    );
                }
            }
        }
    }
    if (points == 0U) {
        throw std::runtime_error(
            "bounded feasible group chamber has no enumerated point"
        );
    }
    std::cout
        << "SU2_T4_GROUP_BOUNDED"
        << " position=" << chamber.mask
        << " Q_bound=" << *q_bound
        << " H_bound=" << *h_bound
        << " points=" << points
        << " result=PASS_EXACT_FINITE"
        << std::endl;
    return true;
}

Integer floor_rational(const Rational& value) {
    const Integer numerator = value.numerator();
    const Integer denominator = value.denominator();
    if (numerator >= 0) {
        return numerator / denominator;
    }
    return -((-numerator + denominator - 1) / denominator);
}

Integer ceil_rational(const Rational& value) {
    return -floor_rational(-value);
}

bool bounded_qy_h_ray_newton_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    const std::optional<Rational> q_bound =
        affine_upper_bound(constraints, 0U);
    const std::optional<Rational> y_bound =
        affine_upper_bound(constraints, 2U);
    if (
        !q_bound.has_value()
        || !y_bound.has_value()
        || *q_bound < 1
        || *y_bound < 0
    ) {
        return false;
    }
    const Integer maximum_q_integer = floor_rational(*q_bound);
    const Integer maximum_y_integer = floor_rational(*y_bound);
    if (
        maximum_q_integer > 1000
        || maximum_y_integer > 10000
    ) {
        return false;
    }
    const int maximum_q = maximum_q_integer.convert_to<int>();
    const int maximum_y = maximum_y_integer.convert_to<int>();
    std::uint64_t rays = 0U;
    for (int q_value = 1; q_value <= maximum_q; ++q_value) {
        for (int y_value = 0; y_value <= maximum_y; ++y_value) {
            Integer minimum_h = 0;
            std::optional<Integer> maximum_h;
            bool feasible = true;
            for (const Polynomial& constraint : constraints) {
                const Affine affine = affine_coefficients(constraint);
                const Rational constant_part =
                    affine[0]
                    + affine[1] * q_value
                    + affine[3] * y_value;
                const Rational h_coefficient = affine[2];
                if (h_coefficient == 0) {
                    if (constant_part < 0) {
                        feasible = false;
                        break;
                    }
                    continue;
                }
                if (h_coefficient > 0) {
                    minimum_h = std::max(
                        minimum_h,
                        ceil_rational(-constant_part / h_coefficient)
                    );
                } else {
                    const Integer bound = floor_rational(
                        constant_part / (-h_coefficient)
                    );
                    if (!maximum_h.has_value() || bound < *maximum_h) {
                        maximum_h = bound;
                    }
                }
            }
            if (
                !feasible
                || (maximum_h.has_value() && *maximum_h < minimum_h)
            ) {
                continue;
            }
            if (minimum_h > std::numeric_limits<long>::max()) {
                return false;
            }
            const std::array<Polynomial, 3> substitution{
                constant(q_value),
                Polynomial::variable(1)
                    + constant(minimum_h.convert_to<long>()),
                constant(y_value)
            };
            const Polynomial translated = substitute(
                chamber.margin,
                substitution
            );
            if (nonnegative_newton(translated)) {
                ++rays;
                continue;
            }
            if (!maximum_h.has_value() || *maximum_h > 10000) {
                return false;
            }
            for (Integer h_value = minimum_h;
                 h_value <= *maximum_h;
                 ++h_value) {
                if (h_value > std::numeric_limits<int>::max()) {
                    return false;
                }
                const std::array<int, 3> values{
                    q_value,
                    h_value.convert_to<int>(),
                    y_value
                };
                if (evaluate(chamber.margin, values) < 0) {
                    throw std::runtime_error(
                        "bounded Q,Y chamber contains a negative point"
                    );
                }
            }
            ++rays;
        }
    }
    if (rays == 0U) {
        return false;
    }
    std::cout
        << "SU2_T4_GROUP_H_RAY_NEWTON"
        << " position=" << chamber.mask
        << " Q_bound=" << *q_bound
        << " Y_bound=" << *y_bound
        << " rays=" << rays
        << " result=PASS_EXACT"
        << std::endl;
    return true;
}

bool certify_group_chamber(
    const Chamber& chamber,
    const GroupFormula& formula
) {
    static_cast<void>(formula);
    const std::vector<Polynomial> constraints =
        irredundant_constraints(chamber);
    const std::vector<bool> forced_zero =
        forced_zero_constraints(constraints);
    bool passed = bounded_group_integer_certificate(chamber, constraints);
    if (!passed) {
        passed = bounded_qy_h_ray_newton_certificate(
            chamber,
            constraints
        );
    }
    if (!passed) {
        passed = direct_facet_certificate(
            chamber,
            constraints,
            forced_zero
        );
    }
    if (!passed) {
        const std::vector<Polynomial> domain_constraints(
            chamber.constraints.begin(),
            chamber.constraints.begin() + 4
        );
        passed = direct_facet_certificate(chamber, domain_constraints);
    }
    if (!passed) {
        passed = sum_cone_certificate(chamber, constraints);
    }
    if (!passed) {
        passed = constant_three_sum_certificate(
            chamber,
            constraints
        );
    }
    if (!passed) {
        passed = equal_sum_square_cone_certificate(
            chamber,
            constraints
        );
    }
    if (!passed) {
        passed = double_sum_cone_certificate(
            chamber,
            constraints
        );
    }
    if (!passed) {
        passed = pair_cut_certificate(chamber, constraints);
    }
    if (!passed) {
        passed = pair_cut_tree_certificate(
            chamber,
            constraints
        );
    }
    return passed;
}

}  // namespace

#ifndef SU2_T4_GROUPS_EMBEDDED
int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 4) {
            throw std::runtime_error(
                "usage: TARGET [--masks-only|--position INDEX]"
            );
        }
        const std::string target = argv[1];
        bool masks_only = false;
        std::optional<std::size_t> selected_position;
        if (argc == 3 && std::string(argv[2]) == "--masks-only") {
            masks_only = true;
        } else if (
            argc == 4
            && std::string(argv[2]) == "--position"
        ) {
            selected_position =
                static_cast<std::size_t>(std::stoull(argv[3]));
        } else if (argc != 2) {
            throw std::runtime_error(
                "usage: TARGET [--masks-only|--position INDEX]"
            );
        }

        const GroupFormula formula = make_group_formula(target);
        const std::vector<WideMask> masks =
            feasible_group_masks(formula);
        std::cout
            << "SU2_T4_GROUP_MASKS"
            << " target=" << target
            << " hinges=" << formula.hinges.size()
            << " masks=" << masks.size()
            << " result=PASS_EXACT_CENSUS"
            << std::endl;
        if (masks_only) {
            return EXIT_SUCCESS;
        }
        if (
            selected_position.has_value()
            && *selected_position >= masks.size()
        ) {
            throw std::runtime_error(
                "selected chamber position is out of range"
            );
        }

        std::size_t attempted = 0U;
        std::size_t certified = 0U;
        for (std::size_t position = 0;
             position < masks.size();
             ++position) {
            if (
                selected_position.has_value()
                && position != *selected_position
            ) {
                continue;
            }
            ++attempted;
            const Chamber chamber = make_group_chamber(
                formula,
                masks[position],
                target,
                static_cast<std::uint64_t>(position)
            );
            if (certify_group_chamber(chamber, formula)) {
                ++certified;
            } else {
                std::cout
                    << "SU2_T4_GROUP_UNRESOLVED"
                    << " target=" << target
                    << " position=" << position
                    << " mask="
                    << masks[position].convert_to<std::string>()
                    << " result=NEEDS_STRONGER_CERTIFICATE"
                    << std::endl;
            }
            if ((position + 1U) % 25U == 0U) {
                std::cerr
                    << "SU2_T4_GROUP"
                    << " target=" << target
                    << " progress=" << position + 1U
                    << '/' << masks.size()
                    << " certified=" << certified
                    << std::endl;
            }
        }
        const bool complete =
            certified == attempted
            && (
                selected_position.has_value()
                || attempted == masks.size()
            );
        std::cout
            << "SU2_T4_GROUP"
            << " target=" << target
            << " attempted=" << attempted
            << " certified=" << certified
            << " result="
            << (complete ? "PASS_EXACT_CERTIFICATE" : "INCOMPLETE")
            << std::endl;
        return complete ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_T4_GROUP FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
#endif
