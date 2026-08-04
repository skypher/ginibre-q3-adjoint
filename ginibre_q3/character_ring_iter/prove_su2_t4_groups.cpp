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
    } else if (target == "c6") {
        // This is the next anchored member
        // f_12 u_13-f_13 u_12.  Its fan census is a feasibility
        // diagnostic for a uniform anchored-kernel mechanism.
        f_powers = {12, 13};
        u_powers = {12, 13};
    } else {
        throw std::runtime_error("target must be c5, c6, g0, g1, or g2");
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
    if (target == "c6") {
        return block(12, 13);
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

std::optional<int> slice_bound_on_forced_equality(
    const Polynomial& first,
    const Polynomial& second,
    const Polynomial& equality
) {
    const Affine sum = affine_coefficients(first + second);
    const Affine forced = affine_coefficients(equality);
    std::optional<Rational> multiple;
    for (std::size_t index = 1U; index < sum.size(); ++index) {
        if (forced[index] == 0) {
            if (sum[index] != 0) {
                return std::nullopt;
            }
            continue;
        }
        const Rational candidate = -sum[index] / forced[index];
        if (multiple.has_value() && candidate != *multiple) {
            return std::nullopt;
        }
        multiple = candidate;
    }
    if (!multiple.has_value()) {
        return std::nullopt;
    }
    const Rational bound = sum[0] + *multiple * forced[0];
    if (
        bound.denominator() != 1
        || bound.numerator() < 0
        || bound.numerator() > std::numeric_limits<int>::max()
    ) {
        return std::nullopt;
    }
    return bound.numerator().convert_to<int>();
}

Integer positive_gcd(Integer left, Integer right) {
    if (left < 0) {
        left = -left;
    }
    if (right < 0) {
        right = -right;
    }
    while (right != 0) {
        const Integer remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

Integer positive_lcm(const Integer& left, const Integer& right) {
    if (left == 0 || right == 0) {
        return 0;
    }
    return (left / positive_gcd(left, right)) * right;
}

std::optional<int> affine_lattice_period(
    const std::array<Polynomial, 3>& inverse
) {
    Integer period = 1;
    for (const Polynomial& coordinate : inverse) {
        const Affine affine = affine_coefficients(coordinate);
        period = positive_lcm(period, affine[3].denominator());
        if (period > std::numeric_limits<int>::max()) {
            return std::nullopt;
        }
    }
    return period.convert_to<int>();
}

bool integral_inverse_point(
    const std::array<Polynomial, 3>& inverse,
    int slice,
    int residue
) {
    const std::array<int, 3> values{slice, 0, residue};
    return std::all_of(
        inverse.begin(),
        inverse.end(),
        [&values](const Polynomial& coordinate) {
            return evaluate(coordinate, values).denominator() == 1;
        }
    );
}

// If an affine chamber has one forced equality and one bounded complementary
// facet pair modulo that equality, its integer points split into finitely many
// congruence rays.  The first coordinate is a facet a in 0<=a<=M, the second
// is an equality facet e=0, and the third is a nonnegative ray facet c.  The
// inverse facet map makes all admissible c values periodic modulo the exact
// lattice period.  A nonnegative Newton expansion on every resulting ray
// certifies the whole chamber without relaxing its integer lattice.
bool forced_equality_slice_ray_newton_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    for (std::size_t equality_first = 0U;
         equality_first < constraints.size();
         ++equality_first) {
        for (std::size_t equality_second = equality_first + 1U;
             equality_second < constraints.size();
             ++equality_second) {
            const Polynomial equality_sum =
                constraints[equality_first] + constraints[equality_second];
            if (!equality_sum.terms().empty()) {
                continue;
            }
            for (std::size_t slice_first = 0U;
                 slice_first < constraints.size();
                 ++slice_first) {
                if (
                    slice_first == equality_first
                    || slice_first == equality_second
                ) {
                    continue;
                }
                for (std::size_t slice_second = slice_first + 1U;
                     slice_second < constraints.size();
                     ++slice_second) {
                    if (
                        slice_second == equality_first
                        || slice_second == equality_second
                    ) {
                        continue;
                    }
                    const std::optional<int> slice_bound =
                        slice_bound_on_forced_equality(
                            constraints[slice_first],
                            constraints[slice_second],
                            constraints[equality_first]
                        );
                    if (!slice_bound.has_value() || *slice_bound > 64) {
                        continue;
                    }
                    for (std::size_t ray = 0U;
                         ray < constraints.size();
                         ++ray) {
                        if (
                            ray == equality_first
                            || ray == equality_second
                            || ray == slice_first
                            || ray == slice_second
                        ) {
                            continue;
                        }
                        std::array<Polynomial, 3> inverse;
                        try {
                            inverse = inverse_facet_map(
                                std::array<Polynomial, 3>{
                                    constraints[slice_first],
                                    constraints[equality_first],
                                    constraints[ray]
                                }
                            );
                        } catch (const std::runtime_error&) {
                            continue;
                        }
                        const std::optional<int> period =
                            affine_lattice_period(inverse);
                        if (!period.has_value()) {
                            continue;
                        }
                        const Polynomial pulled = substitute(
                            chamber.margin,
                            inverse
                        );
                        std::uint64_t rays = 0U;
                        bool passed = true;
                        for (int slice = 0;
                             slice <= *slice_bound && passed;
                             ++slice) {
                            for (int residue = 0;
                                 residue < *period;
                                 ++residue) {
                                if (!integral_inverse_point(
                                        inverse,
                                        slice,
                                        residue
                                    )) {
                                    continue;
                                }
                                const Polynomial parameter =
                                    Polynomial::variable(0);
                                const Polynomial ray_margin = substitute(
                                    pulled,
                                    std::array<Polynomial, 3>{
                                        constant(slice),
                                        constant(0),
                                        scale(parameter, *period)
                                            + constant(residue)
                                    }
                                );
                                if (!nonnegative_newton(ray_margin)) {
                                    passed = false;
                                    break;
                                }
                                ++rays;
                            }
                        }
                        if (!passed || rays == 0U) {
                            continue;
                        }
                        std::cout
                            << "SU2_T4_GROUP_FORCED_EQUALITY_RAYS"
                            << " position=" << chamber.mask
                            << " equality=("
                            << equality_first << ',' << equality_second << ')'
                            << " slice=("
                            << slice_first << ',' << slice_second << ')'
                            << " ray=" << ray
                            << " bound=" << *slice_bound
                            << " period=" << *period
                            << " rays=" << rays
                            << " result=PASS_EXACT_NEWTON"
                            << std::endl;
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

// A square-cone certificate may lose a congruence sublattice.  If a-b has
// fixed residue modulo g, the two integer branches a>=b and b>=a+1 are better
// parameterised by their residue-corrected multiple-of-g difference.  This
// routine checks that finite refinement for equal pair sums and for a
// pair-sum offset of one.  On each branch the remaining condition is
// x+g*y+z>=minimum, which is reduced to finitely many translated orthants.
std::vector<int> integral_affine_moduli(const Polynomial& polynomial) {
    const Affine affine = affine_coefficients(polynomial);
    for (const Rational& coefficient : affine) {
        if (coefficient.denominator() != 1) {
            return {};
        }
    }
    Integer common = 0;
    for (std::size_t index = 1U; index < affine.size(); ++index) {
        common = positive_gcd(common, affine[index].numerator());
    }
    if (
        common < 2
        || common > std::numeric_limits<int>::max()
    ) {
        return {};
    }
    const int maximum = common.convert_to<int>();
    std::set<int> moduli;
    for (int divisor = 2; divisor <= maximum / divisor; ++divisor) {
        if (maximum % divisor != 0) {
            continue;
        }
        moduli.insert(divisor);
        moduli.insert(maximum / divisor);
    }
    moduli.insert(maximum);
    return {moduli.begin(), moduli.end()};
}

int affine_residue(const Polynomial& polynomial, int modulus) {
    const Affine affine = affine_coefficients(polynomial);
    const int residue = (affine[0].numerator() % modulus).convert_to<int>();
    return residue < 0 ? residue + modulus : residue;
}

bool weighted_modular_orthant_certificate(
    const Polynomial& margin,
    const std::array<Polynomial, 3>& facets,
    int minimum,
    int modulus
) {
    if (modulus < 2) {
        throw std::invalid_argument("modulus must be at least two");
    }
    if (minimum < 0) {
        minimum = 0;
    }
    std::array<Polynomial, 3> inverse;
    try {
        inverse = inverse_facet_map(facets);
    } catch (const std::runtime_error&) {
        return false;
    }
    const Polynomial x = Polynomial::variable(0);
    const Polynomial y = Polynomial::variable(1);
    const Polynomial z = Polynomial::variable(2);
    const Polynomial pulled = substitute(
        substitute(margin, inverse),
        std::array<Polynomial, 3>{x, scale(y, modulus), z}
    );
    const int y_floor = (minimum + modulus - 1) / modulus;
    if (
        !nonnegative_basis(
            substitute(
                pulled,
                std::array<Polynomial, 3>{
                    x,
                    y + constant(y_floor),
                    z
                }
            )
        )
    ) {
        return false;
    }
    for (int y_value = 0; y_value < y_floor; ++y_value) {
        const int remainder = minimum - modulus * y_value;
        if (
            !nonnegative_basis(
                substitute(
                    pulled,
                    std::array<Polynomial, 3>{
                        x + constant(remainder),
                        constant(y_value),
                        z
                    }
                )
            )
        ) {
            return false;
        }
        for (int x_value = 0; x_value < remainder; ++x_value) {
            if (
                !nonnegative_basis(
                    substitute(
                        pulled,
                        std::array<Polynomial, 3>{
                            constant(x_value),
                            constant(y_value),
                            z + constant(remainder - x_value)
                        }
                    )
                )
            ) {
                return false;
            }
        }
    }
    return true;
}

bool parity_unit_offset_square_cone_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    struct Pair {
        std::size_t first;
        std::size_t second;
    };
    std::vector<Pair> pairs;
    for (std::size_t first = 0; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < constraints.size();
             ++second) {
            pairs.push_back(Pair{first, second});
        }
    }
    const auto pair_minimum = [&constraints](const Pair& pair) {
        int minimum = 0;
        for (int candidate = 1; candidate <= 20; ++candidate) {
            if (
                constraints_imply_sum(
                    constraints,
                    constraints[pair.first],
                    constraints[pair.second],
                    candidate
                )
            ) {
                minimum = candidate;
            } else {
                break;
            }
        }
        return minimum;
    };
    for (std::size_t left = 0; left < pairs.size(); ++left) {
        for (std::size_t right = left + 1U;
             right < pairs.size();
             ++right) {
            const Pair* left_pair = &pairs[left];
            const Pair* right_pair = &pairs[right];
            if (
                left_pair->first == right_pair->first
                || left_pair->first == right_pair->second
                || left_pair->second == right_pair->first
                || left_pair->second == right_pair->second
            ) {
                continue;
            }
            Polynomial offset_polynomial =
                constraints[right_pair->first]
                + constraints[right_pair->second]
                - constraints[left_pair->first]
                - constraints[left_pair->second];
            int offset = 0;
            if (offset_polynomial.terms().empty()) {
                offset = 0;
            } else if (
                offset_polynomial.terms().size() == 1U
                && offset_polynomial.terms().begin()->first
                    == Exponent{0, 0, 0}
                && offset_polynomial.terms().begin()->second.denominator()
                    == 1
            ) {
                const Integer numerator = offset_polynomial.terms().begin()
                    ->second.numerator();
                if (numerator == 1) {
                    offset = 1;
                } else if (numerator == -1) {
                    offset = 1;
                    std::swap(left_pair, right_pair);
                } else {
                    continue;
                }
            } else {
                continue;
            }
            const int left_minimum = pair_minimum(*left_pair);
            if (left_minimum == 0) {
                continue;
            }
            const int right_minimum = pair_minimum(*right_pair);
            if (right_minimum != left_minimum + offset) {
                continue;
            }
            const std::array<std::size_t, 2> left_indices{
                left_pair->first,
                left_pair->second
            };
            const std::array<std::size_t, 2> right_indices{
                right_pair->first,
                right_pair->second
            };
            for (std::size_t left_choice = 0;
                 left_choice < 2U;
                 ++left_choice) {
                for (std::size_t right_choice = 0;
                     right_choice < 2U;
                     ++right_choice) {
                    const Polynomial& a =
                        constraints[left_indices[left_choice]];
                    const Polynomial& c =
                        constraints[left_indices[1U - left_choice]];
                    const Polynomial& b =
                        constraints[right_indices[right_choice]];
                    const Polynomial& d =
                        constraints[right_indices[1U - right_choice]];
                    for (const int modulus : integral_affine_moduli(a - b)) {
                        const int residue = affine_residue(a - b, modulus);
                        const int first_shift = residue;
                        const int second_shift = modulus - residue;
                        if (
                            !weighted_modular_orthant_certificate(
                                chamber.margin,
                                std::array<Polynomial, 3>{
                                    b,
                                    a - b - constant(first_shift),
                                    c
                                },
                                left_minimum - first_shift,
                                modulus
                            )
                        ) {
                            continue;
                        }
                        if (
                            !weighted_modular_orthant_certificate(
                                chamber.margin,
                                std::array<Polynomial, 3>{
                                    a,
                                    b - a - constant(second_shift),
                                    d
                                },
                                left_minimum + offset - second_shift,
                                modulus
                            )
                        ) {
                            continue;
                        }
                        std::cout
                            << "SU2_T4_GROUP_MODULAR_PAIR_SUM_CONE"
                            << " position=" << chamber.mask
                            << " pairs=("
                            << left_pair->first << ','
                            << left_pair->second << ';'
                            << right_pair->first << ','
                            << right_pair->second << ')'
                            << " orientations=("
                            << left_choice << ','
                            << right_choice << ')'
                            << " minimum=" << left_minimum
                            << " offset=" << offset
                            << " modulus=" << modulus
                            << " residue=" << residue
                            << " result=PASS_EXACT_MODULAR_SQUARE"
                            << std::endl;
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

using IntegralRecessionRay = std::array<int, 3>;

std::optional<IntegralRecessionRay> primitive_recession_ray(
    const Polynomial& first,
    const Polynomial& second,
    const std::vector<Polynomial>& constraints
) {
    const Affine left = affine_coefficients(first);
    const Affine right = affine_coefficients(second);
    const std::array<Rational, 3> cross{
        left[2] * right[3] - left[3] * right[2],
        left[3] * right[1] - left[1] * right[3],
        left[1] * right[2] - left[2] * right[1]
    };
    if (cross[0] == 0 && cross[1] == 0 && cross[2] == 0) {
        return std::nullopt;
    }
    Integer denominator = 1;
    for (const Rational& coordinate : cross) {
        denominator = positive_lcm(denominator, coordinate.denominator());
    }
    std::array<Integer, 3> scaled{};
    Integer divisor = 0;
    for (std::size_t coordinate = 0U;
         coordinate < scaled.size();
         ++coordinate) {
        scaled[coordinate] = cross[coordinate].numerator()
            * (denominator / cross[coordinate].denominator());
        divisor = positive_gcd(divisor, scaled[coordinate]);
    }
    if (divisor == 0) {
        return std::nullopt;
    }
    for (Integer& coordinate : scaled) {
        coordinate /= divisor;
    }
    if (scaled[0] < 0) {
        for (Integer& coordinate : scaled) {
            coordinate = -coordinate;
        }
    }
    if (
        scaled[0] <= 0
        || scaled[0] > std::numeric_limits<int>::max()
        || scaled[1] < std::numeric_limits<int>::min()
        || scaled[1] > std::numeric_limits<int>::max()
        || scaled[2] < std::numeric_limits<int>::min()
        || scaled[2] > std::numeric_limits<int>::max()
    ) {
        return std::nullopt;
    }
    const IntegralRecessionRay ray{
        scaled[0].convert_to<int>(),
        scaled[1].convert_to<int>(),
        scaled[2].convert_to<int>()
    };
    for (const Polynomial& constraint : constraints) {
        const Affine affine = affine_coefficients(constraint);
        const Rational slope =
            affine[1] * ray[0]
            + affine[2] * ray[1]
            + affine[3] * ray[2];
        if (slope < 0) {
            return std::nullopt;
        }
    }
    return ray;
}

bool single_recession_ray_descends(
    const std::vector<Polynomial>& constraints,
    const IntegralRecessionRay& ray,
    int threshold
) {
    z3::context context;
    z3::solver solver(context);
    const std::array<z3::expr, 3> variables{
        context.int_const("Q"),
        context.int_const("H"),
        context.int_const("Y")
    };
    for (const Polynomial& constraint : constraints) {
        solver.add(z3_affine(context, constraint, variables) >= 0);
    }
    solver.add(variables[0] >= threshold);
    const std::array<Polynomial, 3> shift{
        Polynomial::variable(0) - constant(ray[0]),
        Polynomial::variable(1) - constant(ray[1]),
        Polynomial::variable(2) - constant(ray[2])
    };
    z3::expr leaves_chamber = context.bool_val(false);
    for (const Polynomial& constraint : constraints) {
        leaves_chamber = leaves_chamber
            || z3_affine(
                context,
                substitute(constraint, shift),
                variables
            ) < 0;
    }
    solver.add(leaves_chamber);
    return solver.check() == z3::unsat;
}

bool single_recession_ray_newton_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    if (constraints.size() < 3U || constraints.size() > 8U) {
        return false;
    }
    std::set<IntegralRecessionRay> rays;
    for (std::size_t first = 0U; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < constraints.size();
             ++second) {
            const std::optional<IntegralRecessionRay> ray =
                primitive_recession_ray(
                    constraints[first], constraints[second], constraints
                );
            if (ray.has_value()) {
                rays.insert(*ray);
            }
        }
    }
    if (rays.size() != 1U) {
        return false;
    }
    const IntegralRecessionRay ray = *rays.begin();
    constexpr int maximum_threshold = 96;
    if (!single_recession_ray_descends(
            constraints, ray, maximum_threshold
        )) {
        return false;
    }
    int lower = 1;
    int upper = maximum_threshold;
    while (lower < upper) {
        const int middle = lower + (upper - lower) / 2;
        if (single_recession_ray_descends(constraints, ray, middle)) {
            upper = middle;
        } else {
            lower = middle + 1;
        }
    }
    const int threshold = lower;
    std::vector<Polynomial> bounded_constraints = constraints;
    bounded_constraints.push_back(
        constant(threshold - 1) - Polynomial::variable(0)
    );
    const std::optional<Rational> h_bound =
        affine_upper_bound(bounded_constraints, 1U);
    const std::optional<Rational> y_bound =
        affine_upper_bound(bounded_constraints, 2U);
    if (
        !h_bound.has_value()
        || !y_bound.has_value()
        || *h_bound < 0
        || *y_bound < 0
    ) {
        return false;
    }
    const Integer maximum_h_integer = floor_rational(*h_bound);
    const Integer maximum_y_integer = floor_rational(*y_bound);
    if (maximum_h_integer > 256 || maximum_y_integer > 1024) {
        return false;
    }
    const int maximum_h = maximum_h_integer.convert_to<int>();
    const int maximum_y = maximum_y_integer.convert_to<int>();
    std::uint64_t base_points = 0U;
    for (int q_value = 1; q_value < threshold; ++q_value) {
        for (int h_value = 0; h_value <= maximum_h; ++h_value) {
            for (int y_value = 0; y_value <= maximum_y; ++y_value) {
                const std::array<int, 3> point{q_value, h_value, y_value};
                if (!std::all_of(
                        constraints.begin(), constraints.end(),
                        [&point](const Polynomial& constraint) {
                            return evaluate(constraint, point) >= 0;
                        }
                    )) {
                    continue;
                }
                const Polynomial parameter = Polynomial::variable(0);
                const std::array<Polynomial, 3> substitution{
                    constant(q_value) + scale(parameter, ray[0]),
                    constant(h_value) + scale(parameter, ray[1]),
                    constant(y_value) + scale(parameter, ray[2])
                };
                if (!nonnegative_basis(
                        substitute(chamber.margin, substitution)
                    )) {
                    return false;
                }
                ++base_points;
            }
        }
    }
    if (base_points == 0U) {
        return false;
    }
    std::cout
        << "SU2_T4_GROUP_SINGLE_RECESSION_RAY"
        << " position=" << chamber.mask
        << " ray=(" << ray[0] << ',' << ray[1] << ',' << ray[2] << ')'
        << " threshold=" << threshold
        << " base_points=" << base_points
        << " H_bound=" << *h_bound
        << " Y_bound=" << *y_bound
        << " result=PASS_EXACT_NEWTON"
        << std::endl;
    return true;
}

bool recession_ray_family_descends(
    const std::vector<Polynomial>& constraints,
    const std::vector<IntegralRecessionRay>& rays,
    int threshold
) {
    z3::context context;
    z3::solver solver(context);
    const std::array<z3::expr, 3> variables{
        context.int_const("Q"),
        context.int_const("H"),
        context.int_const("Y")
    };
    for (const Polynomial& constraint : constraints) {
        solver.add(z3_affine(context, constraint, variables) >= 0);
    }
    solver.add(variables[0] >= threshold);
    for (const IntegralRecessionRay& ray : rays) {
        const std::array<Polynomial, 3> shift{
            Polynomial::variable(0) - constant(ray[0]),
            Polynomial::variable(1) - constant(ray[1]),
            Polynomial::variable(2) - constant(ray[2])
        };
        z3::expr leaves_chamber = context.bool_val(false);
        for (const Polynomial& constraint : constraints) {
            leaves_chamber = leaves_chamber
                || z3_affine(
                    context,
                    substitute(constraint, shift),
                    variables
                ) < 0;
        }
        solver.add(leaves_chamber);
    }
    return solver.check() == z3::unsat;
}

bool recession_fan_newton_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints
) {
    if (constraints.size() < 3U || constraints.size() > 8U) {
        return false;
    }
    std::set<IntegralRecessionRay> ray_set;
    for (std::size_t first = 0U; first < constraints.size(); ++first) {
        for (std::size_t second = first + 1U;
             second < constraints.size();
             ++second) {
            const std::optional<IntegralRecessionRay> ray =
                primitive_recession_ray(
                    constraints[first], constraints[second], constraints
                );
            if (ray.has_value()) {
                ray_set.insert(*ray);
            }
        }
    }
    if (ray_set.size() < 2U || ray_set.size() > 3U) {
        return false;
    }
    const std::vector<IntegralRecessionRay> rays{
        ray_set.begin(), ray_set.end()
    };
    constexpr int maximum_threshold = 96;
    if (!recession_ray_family_descends(
            constraints, rays, maximum_threshold
        )) {
        return false;
    }
    int lower = 1;
    int upper = maximum_threshold;
    while (lower < upper) {
        const int middle = lower + (upper - lower) / 2;
        if (recession_ray_family_descends(constraints, rays, middle)) {
            upper = middle;
        } else {
            lower = middle + 1;
        }
    }
    const int threshold = lower;
    std::vector<Polynomial> bounded_constraints = constraints;
    bounded_constraints.push_back(
        constant(threshold - 1) - Polynomial::variable(0)
    );
    const std::optional<Rational> h_bound =
        affine_upper_bound(bounded_constraints, 1U);
    const std::optional<Rational> y_bound =
        affine_upper_bound(bounded_constraints, 2U);
    if (
        !h_bound.has_value()
        || !y_bound.has_value()
        || *h_bound < 0
        || *y_bound < 0
    ) {
        return false;
    }
    const Integer maximum_h_integer = floor_rational(*h_bound);
    const Integer maximum_y_integer = floor_rational(*y_bound);
    if (maximum_h_integer > 256 || maximum_y_integer > 1024) {
        return false;
    }
    const int maximum_h = maximum_h_integer.convert_to<int>();
    const int maximum_y = maximum_y_integer.convert_to<int>();
    constexpr std::size_t maximum_base_points = 64U;
    std::vector<std::array<int, 3>> bases;
    for (int q_value = 1; q_value < threshold; ++q_value) {
        for (int h_value = 0; h_value <= maximum_h; ++h_value) {
            for (int y_value = 0; y_value <= maximum_y; ++y_value) {
                const std::array<int, 3> point{q_value, h_value, y_value};
                if (!std::all_of(
                        constraints.begin(), constraints.end(),
                        [&point](const Polynomial& constraint) {
                            return evaluate(constraint, point) >= 0;
                        }
                    )) {
                    continue;
                }
                bool reducible = false;
                for (const IntegralRecessionRay& ray : rays) {
                    const std::array<int, 3> predecessor{
                        q_value - ray[0],
                        h_value - ray[1],
                        y_value - ray[2]
                    };
                    if (std::all_of(
                            constraints.begin(), constraints.end(),
                            [&predecessor](const Polynomial& constraint) {
                                return evaluate(constraint, predecessor) >= 0;
                            }
                        )) {
                        reducible = true;
                        break;
                    }
                }
                if (reducible) {
                    continue;
                }
                bases.push_back(point);
                if (bases.size() > maximum_base_points) {
                    return false;
                }
            }
        }
    }
    if (bases.empty()) {
        return false;
    }
    std::cerr
        << "SU2_T4_GROUP_RECESSION_FAN_PROGRESS"
        << " position=" << chamber.mask
        << " phase=base_newton"
        << " bases=" << bases.size()
        << std::endl;
    for (std::size_t base_index = 0U; base_index < bases.size(); ++base_index) {
        const std::array<int, 3>& point = bases[base_index];
        std::cerr
            << "SU2_T4_GROUP_RECESSION_FAN_PROGRESS"
            << " position=" << chamber.mask
            << " phase=base_newton"
            << " base=" << base_index + 1U
            << '/' << bases.size()
            << std::endl;
        std::array<Polynomial, 3> substitution{
            constant(point[0]), constant(point[1]), constant(point[2])
        };
        for (std::size_t index = 0U; index < rays.size(); ++index) {
            const Polynomial parameter = Polynomial::variable(index);
            for (std::size_t coordinate = 0U;
                 coordinate < substitution.size();
                 ++coordinate) {
                substitution[coordinate] += scale(
                    parameter, rays[index][coordinate]
                );
            }
        }
        if (!nonnegative_basis(substitute(chamber.margin, substitution))) {
            return false;
        }
    }
    std::cout
        << "SU2_T4_GROUP_RECESSION_FAN"
        << " position=" << chamber.mask
        << " rays=" << rays.size()
        << " ray_vectors=(";
    for (std::size_t ray_index = 0U; ray_index < rays.size(); ++ray_index) {
        if (ray_index != 0U) {
            std::cout << ',';
        }
        const IntegralRecessionRay& ray = rays[ray_index];
        std::cout << '(' << ray[0] << ',' << ray[1] << ',' << ray[2] << ')';
    }
    std::cout
        << ')'
        << " threshold=" << threshold
        << " base_points=" << bases.size()
        << " H_bound=" << *h_bound
        << " Y_bound=" << *y_bound
        << " result=PASS_EXACT_NEWTON"
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
        passed = forced_equality_slice_ray_newton_certificate(
            chamber,
            constraints
        );
    }
    if (!passed) {
        passed = parity_unit_offset_square_cone_certificate(
            chamber,
            constraints
        );
    }
    if (!passed) {
        passed = single_recession_ray_newton_certificate(
            chamber,
            constraints
        );
    }
    if (!passed) {
        passed = recession_fan_newton_certificate(chamber, constraints);
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
