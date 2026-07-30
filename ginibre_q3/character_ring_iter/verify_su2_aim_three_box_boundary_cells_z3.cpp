#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include <z3++.h>

namespace {

int parse_case(const char* text, const char* name, int maximum) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long parsed = std::stol(value, &consumed, 10);
    if (consumed != value.size() || parsed < 0 || parsed >= maximum) {
        throw std::invalid_argument(
            std::string(name) + " is outside its case range");
    }
    return static_cast<int>(parsed);
}

struct GuardedFormula {
    explicit GuardedFormula(z3::context& context_value)
        : context(context_value) {}

    z3::expr record(const z3::expr& guard) {
        const unsigned identifier = Z3_get_ast_id(context, guard);
        if (identifiers.insert(identifier).second) {
            guards.push_back(guard);
        }
        return guard;
    }

    z3::expr twice_quadratic_spline(const z3::expr& degree) {
        return z3::ite(
            record(degree >= 0),
            (degree + 1) * (degree + 2),
            context.int_val(0));
    }

    z3::expr twice_three_box_coefficient(
        const z3::expr& shortest,
        const z3::expr& middle,
        const z3::expr& longest,
        const z3::expr& support_end,
        const z3::expr& index) {
        const z3::expr reflected = z3::ite(
            record(2 * index <= support_end),
            index,
            support_end - index);
        const z3::expr degree = reflected - 1;
        const z3::expr half_coefficient
            = (degree + 1) * (degree + 2)
              - twice_quadratic_spline(degree - shortest)
              - twice_quadratic_spline(degree - middle)
              - twice_quadratic_spline(degree - longest)
              + twice_quadratic_spline(
                    degree - shortest - middle);
        return z3::ite(
            record(index < support_end),
            half_coefficient,
            context.int_val(0));
    }

    z3::context& context;
    std::vector<z3::expr> guards;
    std::unordered_set<unsigned> identifiers;
};

const char* render(const z3::check_result& result) {
    if (result == z3::unsat) {
        return "unsat";
    }
    if (result == z3::sat) {
        return "sat";
    }
    return "unknown";
}

bool contains_ite(const z3::expr& expression) {
    if (expression.is_app()
        && expression.decl().decl_kind() == Z3_OP_ITE) {
        return true;
    }
    for (unsigned index = 0;
         index < expression.num_args();
         ++index) {
        if (contains_ite(expression.arg(index))) {
            return true;
        }
    }
    return false;
}

void collect_nonnegative_affine(
    const z3::expr& input,
    bool truth,
    std::vector<z3::expr>& constraints,
    std::set<std::string>& seen) {
    const z3::expr formula = input.simplify();
    const Z3_decl_kind kind = formula.decl().decl_kind();
    if (kind == Z3_OP_TRUE) {
        if (!truth) {
            throw std::logic_error("infeasible false truth assignment");
        }
        return;
    }
    if (kind == Z3_OP_FALSE) {
        if (truth) {
            throw std::logic_error("infeasible true truth assignment");
        }
        return;
    }
    if (kind == Z3_OP_NOT) {
        collect_nonnegative_affine(
            formula.arg(0),
            !truth,
            constraints,
            seen);
        return;
    }
    if (kind == Z3_OP_AND && truth) {
        for (unsigned index = 0; index < formula.num_args(); ++index) {
            collect_nonnegative_affine(
                formula.arg(index),
                true,
                constraints,
                seen);
        }
        return;
    }

    z3::expr nonnegative = formula.ctx().int_val(0);
    if (kind == Z3_OP_GE) {
        nonnegative = truth
            ? formula.arg(0) - formula.arg(1)
            : formula.arg(1) - formula.arg(0) - 1;
    } else if (kind == Z3_OP_LE) {
        nonnegative = truth
            ? formula.arg(1) - formula.arg(0)
            : formula.arg(0) - formula.arg(1) - 1;
    } else if (kind == Z3_OP_GT) {
        nonnegative = truth
            ? formula.arg(0) - formula.arg(1) - 1
            : formula.arg(1) - formula.arg(0);
    } else if (kind == Z3_OP_LT) {
        nonnegative = truth
            ? formula.arg(1) - formula.arg(0) - 1
            : formula.arg(0) - formula.arg(1);
    } else if (kind == Z3_OP_EQ && truth) {
        collect_nonnegative_affine(
            formula.arg(0) >= formula.arg(1),
            true,
            constraints,
            seen);
        collect_nonnegative_affine(
            formula.arg(1) >= formula.arg(0),
            true,
            constraints,
            seen);
        return;
    } else {
        throw std::logic_error(
            "unsupported affine guard: " + formula.to_string());
    }
    nonnegative = nonnegative.simplify();
    if (contains_ite(nonnegative)) {
        throw std::logic_error(
            "residual ite in affine guard: "
            + nonnegative.to_string());
    }
    const std::string key = nonnegative.to_string();
    if (seen.insert(key).second) {
        constraints.push_back(nonnegative);
    }
}

bool handelman_certificate(
    z3::context& context,
    const z3::expr_vector& base,
    const std::vector<z3::expr>& guards,
    const std::vector<bool>& values,
    const z3::expr& margin,
    const z3::expr& parity,
    int parity_case,
    const std::vector<z3::expr>& variables,
    std::size_t& constraint_count,
    std::size_t& basis_count,
    bool& real_certificate,
    std::string& rendered_constraints,
    std::string& relaxation_witness,
    std::size_t& residue_queries,
    int& maximum_scale) {
    z3::expr_vector parity_from(context);
    z3::expr_vector parity_to(context);
    parity_from.push_back(parity);
    parity_to.push_back(context.int_val(parity_case));

    std::vector<z3::expr> constraints;
    std::set<std::string> seen;
    for (unsigned index = 0; index < base.size(); ++index) {
        z3::expr condition = base[static_cast<int>(index)];
        collect_nonnegative_affine(
            condition.substitute(parity_from, parity_to),
            true,
            constraints,
            seen);
    }
    const std::size_t protected_constraints = constraints.size();
    for (std::size_t index = 0; index < guards.size(); ++index) {
        z3::expr condition = guards[index];
        z3::expr_vector branch_from(context);
        z3::expr_vector branch_to(context);
        for (std::size_t other = 0;
             other < guards.size();
             ++other) {
            if (other != index) {
                branch_from.push_back(guards[other]);
                branch_to.push_back(
                    context.bool_val(values[other]));
            }
        }
        condition = condition
                        .substitute(branch_from, branch_to)
                        .simplify();
        collect_nonnegative_affine(
            condition.substitute(parity_from, parity_to),
            values[index],
            constraints,
            seen);
    }
    for (std::size_t candidate = protected_constraints;
         candidate < constraints.size();) {
        z3::solver redundancy(context, "QF_LIA");
        for (std::size_t index = 0;
             index < constraints.size();
             ++index) {
            if (index != candidate) {
                redundancy.add(constraints[index] >= 0);
            }
        }
        redundancy.add(constraints[candidate] < 0);
        if (redundancy.check() == z3::unsat) {
            constraints.erase(
                constraints.begin()
                + static_cast<std::ptrdiff_t>(candidate));
        } else {
            ++candidate;
        }
    }
    z3::expr fixed_margin = margin;
    fixed_margin
        = fixed_margin.substitute(parity_from, parity_to).simplify();

    std::vector<z3::expr> basis;
    basis.push_back(context.int_val(1));
    for (const z3::expr& constraint : constraints) {
        basis.push_back(constraint);
    }
    for (std::size_t first = 0;
         first < constraints.size();
         ++first) {
        for (std::size_t second = first;
             second < constraints.size();
             ++second) {
            basis.push_back(
                constraints[first] * constraints[second]);
        }
    }
    constraint_count = constraints.size();
    basis_count = basis.size();
    rendered_constraints.clear();
    for (const z3::expr& constraint : constraints) {
        rendered_constraints += "\n  ";
        rendered_constraints += constraint.to_string();
        rendered_constraints += " >= 0";
    }

    std::vector<z3::expr> multipliers;
    multipliers.reserve(basis.size());
    z3::solver solver(context, "QF_LRA");
    z3::params parameters(context);
    parameters.set("timeout", 60000U);
    solver.set(parameters);
    for (std::size_t index = 0; index < basis.size(); ++index) {
        const z3::expr multiplier = context.real_const(
            ("lambda_" + std::to_string(index)).c_str());
        multipliers.push_back(multiplier);
        solver.add(multiplier >= 0);
    }

    std::vector<std::vector<int>> samples;
    samples.push_back(std::vector<int>(variables.size(), 0));
    for (std::size_t index = 0; index < variables.size(); ++index) {
        std::vector<int> first(variables.size(), 0);
        first[index] = 1;
        samples.push_back(first);
        std::vector<int> second(variables.size(), 0);
        second[index] = 2;
        samples.push_back(second);
    }
    for (std::size_t first = 0;
         first < variables.size();
         ++first) {
        for (std::size_t second = first + 1;
             second < variables.size();
             ++second) {
            std::vector<int> sample(variables.size(), 0);
            sample[first] = 1;
            sample[second] = 1;
            samples.push_back(sample);
        }
    }

    z3::expr_vector variable_vector(context);
    for (const z3::expr& variable : variables) {
        variable_vector.push_back(variable);
    }
    for (const std::vector<int>& sample : samples) {
        z3::expr_vector values_vector(context);
        for (int value : sample) {
            values_vector.push_back(context.int_val(value));
        }
        z3::expr margin_value = fixed_margin;
        margin_value = margin_value
                           .substitute(
                               variable_vector,
                               values_vector)
                           .simplify();
        z3::expr certificate_value = context.real_val(0);
        for (std::size_t index = 0; index < basis.size(); ++index) {
            z3::expr basis_value = basis[index];
            basis_value = basis_value
                              .substitute(
                                  variable_vector,
                                  values_vector)
                              .simplify();
            certificate_value = certificate_value
                + multipliers[index]
                      * basis_value;
        }
        solver.add(certificate_value == margin_value);
    }
    if (solver.check() == z3::sat) {
        const z3::model certificate_model = solver.get_model();
        std::vector<z3::expr> coefficients;
        coefficients.reserve(multipliers.size());
        for (const z3::expr& multiplier : multipliers) {
            const z3::expr coefficient
                = certificate_model.eval(multiplier, true);
            if (!(coefficient >= 0).simplify().is_true()) {
                throw std::logic_error(
                    "negative Handelman coefficient");
            }
            coefficients.push_back(coefficient);
        }
        for (const std::vector<int>& sample : samples) {
            z3::expr_vector values_vector(context);
            for (int value : sample) {
                values_vector.push_back(context.int_val(value));
            }
            z3::expr margin_value = fixed_margin;
            margin_value = margin_value
                               .substitute(
                                   variable_vector,
                                   values_vector)
                               .simplify();
            z3::expr certificate_value = context.real_val(0);
            for (std::size_t index = 0;
                 index < basis.size();
                 ++index) {
                z3::expr basis_value = basis[index];
                basis_value = basis_value
                                  .substitute(
                                      variable_vector,
                                      values_vector)
                                  .simplify();
                certificate_value = certificate_value
                    + coefficients[index] * basis_value;
            }
            if (!(certificate_value == margin_value)
                     .simplify()
                     .is_true()) {
                throw std::logic_error(
                    "Handelman interpolation identity mismatch");
            }
        }
        real_certificate = false;
        return true;
    }

    std::vector<z3::expr> real_variables;
    real_variables.reserve(variables.size());
    for (std::size_t index = 0; index < variables.size(); ++index) {
        real_variables.push_back(context.real_const(
            ("real_" + std::to_string(index)).c_str()));
    }
    const auto lift_degree_two =
        [&context, &variables, &real_variables](
            const z3::expr& polynomial) {
            z3::expr_vector lifted_variables(context);
            for (const z3::expr& variable : variables) {
                lifted_variables.push_back(variable);
            }
            const auto evaluate =
                [&context, &lifted_variables](
                    const z3::expr& expression,
                    const std::vector<int>& point) {
                    z3::expr_vector values_vector(context);
                    for (int value : point) {
                        values_vector.push_back(
                            context.int_val(value));
                    }
                    z3::expr copy = expression;
                    return z3::to_real(
                        copy.substitute(
                                lifted_variables,
                                values_vector)
                            .simplify());
                };
            const std::vector<int> zero(variables.size(), 0);
            const z3::expr constant = evaluate(polynomial, zero);
            z3::expr lifted = constant;
            std::vector<z3::expr> linear;
            std::vector<z3::expr> diagonal;
            linear.reserve(variables.size());
            diagonal.reserve(variables.size());
            for (std::size_t index = 0;
                 index < variables.size();
                 ++index) {
                std::vector<int> first(variables.size(), 0);
                std::vector<int> second(variables.size(), 0);
                first[index] = 1;
                second[index] = 2;
                const z3::expr first_value
                    = evaluate(polynomial, first);
                const z3::expr second_value
                    = evaluate(polynomial, second);
                const z3::expr quadratic
                    = (second_value - 2 * first_value + constant)
                    / context.real_val(2);
                const z3::expr affine
                    = first_value - constant - quadratic;
                diagonal.push_back(quadratic);
                linear.push_back(affine);
                lifted = lifted
                    + affine * real_variables[index]
                    + quadratic * real_variables[index]
                          * real_variables[index];
            }
            for (std::size_t first = 0;
                 first < variables.size();
                 ++first) {
                for (std::size_t second = first + 1;
                     second < variables.size();
                     ++second) {
                    std::vector<int> point(variables.size(), 0);
                    point[first] = 1;
                    point[second] = 1;
                    const z3::expr cross
                        = evaluate(polynomial, point)
                          - constant
                          - linear[first]
                          - linear[second]
                          - diagonal[first]
                          - diagonal[second];
                    lifted = lifted
                        + cross * real_variables[first]
                              * real_variables[second];
                }
            }
            return lifted.simplify();
        };

    std::vector<z3::expr> lifted_constraints;
    lifted_constraints.reserve(constraints.size());
    for (const z3::expr& constraint : constraints) {
        lifted_constraints.push_back(
            lift_degree_two(constraint));
    }
    z3::expr lifted_margin = lift_degree_two(fixed_margin);
    std::vector<z3::expr> coordinate_variables;
    coordinate_variables.reserve(variables.size());
    for (std::size_t index = 0; index < variables.size(); ++index) {
        coordinate_variables.push_back(context.real_const(
            ("coordinate_" + std::to_string(index)).c_str()));
    }
    const z3::expr& coordinate_m = coordinate_variables[0];
    const z3::expr& coordinate_padding = coordinate_variables[1];
    const z3::expr& coordinate_load = coordinate_variables[2];
    const z3::expr& coordinate_width_slack
        = coordinate_variables[3];
    const z3::expr& coordinate_parity_slack
        = coordinate_variables[4];
    z3::expr_vector coordinate_from(context);
    for (const z3::expr& variable : real_variables) {
        coordinate_from.push_back(variable);
    }
    z3::expr_vector coordinate_to(context);
    coordinate_to.push_back(coordinate_m);
    coordinate_to.push_back(
        2 * coordinate_padding
        + coordinate_parity_slack
        + coordinate_width_slack);
    coordinate_to.push_back(
        2 * coordinate_padding + coordinate_parity_slack);
    coordinate_to.push_back(coordinate_padding);
    coordinate_to.push_back(
        4 * coordinate_m + 2 * parity_case + 1
        + coordinate_load - 2 * coordinate_padding
        - coordinate_parity_slack);
    for (z3::expr& constraint : lifted_constraints) {
        constraint = constraint
                         .substitute(coordinate_from, coordinate_to)
                         .simplify();
    }
    lifted_margin = lifted_margin
                        .substitute(coordinate_from, coordinate_to)
                        .simplify();
    const unsigned residue_cases
        = 1U << static_cast<unsigned>(variables.size());
    const int maximum_depth = 3;
    std::function<bool(int, const std::vector<int>&, int)>
        prove_residue;
    prove_residue =
        [&](int scale,
            const std::vector<int>& residues,
            int depth) {
        ++residue_queries;
        maximum_scale = std::max(maximum_scale, scale);
        z3::expr_vector residue_from(context);
        z3::expr_vector residue_to(context);
        for (std::size_t index = 0;
             index < coordinate_variables.size();
             ++index) {
            residue_from.push_back(coordinate_variables[index]);
            const z3::expr quotient = context.real_const(
                ("quotient_" + std::to_string(index)).c_str());
            residue_to.push_back(
                scale * quotient
                + context.int_val(residues[index]));
        }
        z3::solver real_solver(context, "QF_NRA");
        z3::params real_parameters(context);
        real_parameters.set("timeout", 10000U);
        real_solver.set(real_parameters);
        for (const z3::expr& constraint : lifted_constraints) {
            z3::expr copy = constraint;
            real_solver.add(
                copy.substitute(residue_from, residue_to) >= 0);
        }
        z3::expr margin_copy = lifted_margin;
        real_solver.add(
            margin_copy.substitute(residue_from, residue_to)
            <= context.real_val(-2));
        const z3::check_result real_result = real_solver.check();
        if (real_result == z3::unsat) {
            return true;
        }
        if (depth >= maximum_depth) {
            relaxation_witness
                = "scale=" + std::to_string(scale)
                  + " result=" + render(real_result)
                  + " residues=";
            for (int residue : residues) {
                relaxation_witness += std::to_string(residue) + ",";
            }
            if (real_result == z3::sat) {
                relaxation_witness += " model="
                    + real_solver.get_model().to_string();
            } else {
                relaxation_witness += " reason="
                    + real_solver.reason_unknown();
            }
            relaxation_witness += " polynomial="
                + lifted_margin.to_string();
            return false;
        }
        for (unsigned mask = 0; mask < residue_cases; ++mask) {
            std::vector<int> children = residues;
            for (std::size_t index = 0;
                 index < variables.size();
                 ++index) {
                if ((mask
                     & (1U << static_cast<unsigned>(index)))
                    != 0U) {
                    children[index] += scale;
                }
            }
            if (!prove_residue(2 * scale, children, depth + 1)) {
                return false;
            }
        }
        return true;
    };
    real_certificate = prove_residue(
        1,
        std::vector<int>(variables.size(), 0),
        0);
    return real_certificate;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 6) {
            throw std::invalid_argument(
                "usage: verify_su2_aim_three_box_boundary_cells_z3 "
                "support_band width_order shift_band pair_order parity");
        }
        const int support_band
            = parse_case(argv[1], "support_band", 2);
        const int width_order
            = parse_case(argv[2], "width_order", 3);
        const int shift_band
            = parse_case(argv[3], "shift_band", 2);
        const int pair_order
            = parse_case(argv[4], "pair_order", 2);
        const int parity_case
            = parse_case(argv[5], "parity", 2);

        z3::context context;
        const z3::expr half_quotient = context.int_const("m");
        const z3::expr parity = context.int_const("p");
        const z3::expr half_level
            = 2 * half_quotient + parity;
        const z3::expr period = 2 * half_level + 2;
        const z3::expr left = context.int_const("ell");
        const z3::expr width = context.int_const("w");
        const z3::expr padding = context.int_const("A");
        const z3::expr shift = left + 2 * padding;
        const z3::expr label = context.int_const("s");
        const z3::expr maximum_label = half_level - parity;
        const z3::expr support_end
            = left + width + shift - 1;
        z3::expr shortest = left;
        z3::expr middle = width;
        z3::expr longest = shift;

        z3::expr_vector base(context);
        base.push_back(half_quotient >= 1);
        base.push_back(parity == parity_case);
        base.push_back(half_level >= 2);
        base.push_back(left >= 1 && left < period);
        base.push_back(width >= 1 && width < period);
        base.push_back(padding >= 0 && padding <= half_level);
        base.push_back(label >= 1 && label <= maximum_label);
        base.push_back(support_end > 2 * period);
        base.push_back(
            support_band == 0
                ? support_end <= 3 * period
                : support_end > 3 * period);
        if (width_order == 0) {
            base.push_back(width <= left);
            shortest = width;
            middle = left;
            longest = shift;
        } else if (width_order == 1) {
            base.push_back(width > left && width <= shift);
            shortest = left;
            middle = width;
            longest = shift;
        } else {
            base.push_back(width > shift);
            shortest = left;
            middle = shift;
            longest = width;
        }
        base.push_back(
            shift_band == 0 ? shift < period : shift >= period);
        base.push_back(
            pair_order == 0
                ? longest <= shortest + middle
                : longest > shortest + middle);

        GuardedFormula formula(context);
        z3::expr twice_margin
            = formula.twice_three_box_coefficient(
                shortest,
                middle,
                longest,
                support_end,
                period - label - 1);
        for (int wall = 1; wall <= 3; ++wall) {
            twice_margin = twice_margin
                - formula.twice_three_box_coefficient(
                      shortest,
                      middle,
                      longest,
                      support_end,
                      wall * period + label)
                + formula.twice_three_box_coefficient(
                      shortest,
                      middle,
                      longest,
                      support_end,
                      (wall + 1) * period - label - 1);
        }

        z3::solver enumerator(context, "QF_LIA");
        enumerator.add(base);
        std::uint64_t cells = 0U;
        std::uint64_t handelman_cells = 0U;
        std::uint64_t real_cells = 0U;
        std::uint64_t nia_cells = 0U;
        std::uint64_t residue_queries = 0U;
        int maximum_scale = 1;
        while (enumerator.check() == z3::sat) {
            const z3::model model = enumerator.get_model();
            std::vector<bool> values;
            values.reserve(formula.guards.size());
            z3::expr_vector block(context);
            z3::expr_vector substitutions(context);
            z3::expr_vector truth_values(context);
            for (const z3::expr& guard : formula.guards) {
                const bool value
                    = model.eval(guard, true).is_true();
                values.push_back(value);
                block.push_back(value ? !guard : guard);
                substitutions.push_back(guard);
                truth_values.push_back(context.bool_val(value));
            }
            enumerator.add(z3::mk_or(block));
            const z3::expr cell_margin = twice_margin
                .substitute(substitutions, truth_values)
                .simplify();
            std::size_t constraint_count = 0U;
            std::size_t basis_count = 0U;
            bool real_certificate = false;
            std::string rendered_constraints;
            std::string relaxation_witness;
            std::size_t cell_residue_queries = 0U;
            int cell_maximum_scale = 1;
            const std::vector<z3::expr> variables{
                half_quotient,
                left,
                width,
                padding,
                label};
            if (handelman_certificate(
                    context,
                    base,
                    formula.guards,
                    values,
                    cell_margin,
                    parity,
                    parity_case,
                    variables,
                    constraint_count,
                    basis_count,
                    real_certificate,
                    rendered_constraints,
                    relaxation_witness,
                    cell_residue_queries,
                    cell_maximum_scale)) {
                ++cells;
                residue_queries += cell_residue_queries;
                maximum_scale = std::max(
                    maximum_scale,
                    cell_maximum_scale);
                if (real_certificate) {
                    ++real_cells;
                } else {
                    ++handelman_cells;
                }
                continue;
            }

            z3::solver cell_solver(context, "QF_NIA");
            z3::params parameters(context);
            parameters.set("timeout", 10000U);
            cell_solver.set(parameters);
            cell_solver.add(base);
            for (std::size_t index = 0;
                 index < formula.guards.size();
                 ++index) {
                cell_solver.add(
                    values[index]
                        ? formula.guards[index]
                        : !formula.guards[index]);
            }
            cell_solver.add(cell_margin < 0);
            const z3::check_result result = cell_solver.check();
            ++cells;
            residue_queries += cell_residue_queries;
            maximum_scale = std::max(
                maximum_scale,
                cell_maximum_scale);
            if (result != z3::unsat) {
                std::cout
                    << "SU2_AIM_THREE_BOX_BOUNDARY_CELLS_Z3"
                    << " logic=QF_LIA+QF_LRA+QF_NRA+QF_NIA"
                    << " support_band=" << support_band
                    << " width_order=" << width_order
                    << " shift_band=" << shift_band
                    << " pair_order=" << pair_order
                    << " parity=" << parity_case
                    << " guards=" << formula.guards.size()
                    << " cells=" << cells
                    << " handelman_cells=" << handelman_cells
                    << " real_cells=" << real_cells
                    << " nia_cells=" << nia_cells
                    << " residue_queries=" << residue_queries
                    << " maximum_scale=" << maximum_scale
                    << " constraints=" << constraint_count
                    << " basis=" << basis_count
                    << " result=" << render(result) << '\n';
                if (result == z3::sat) {
                    std::cout << "MODEL "
                              << cell_solver.get_model() << '\n';
                } else {
                    std::cout << "REASON "
                              << cell_solver.reason_unknown() << '\n';
                    std::cout << "MARGIN " << cell_margin << '\n';
                    z3::params normalization(context);
                    normalization.set("som", true);
                    std::cout << "EXPANDED_MARGIN "
                              << cell_margin.simplify(normalization)
                              << '\n';
                    std::cout << "CONSTRAINTS"
                              << rendered_constraints << '\n';
                    std::cout << "RELAXATION "
                              << relaxation_witness << '\n';
                }
                return EXIT_FAILURE;
            }
            ++nia_cells;
        }
        std::cout << "SU2_AIM_THREE_BOX_BOUNDARY_CELLS_Z3"
                  << " logic=QF_LIA+QF_LRA+QF_NRA+QF_NIA"
                  << " support_band=" << support_band
                  << " width_order=" << width_order
                  << " shift_band=" << shift_band
                  << " pair_order=" << pair_order
                  << " parity=" << parity_case
                  << " guards=" << formula.guards.size()
                  << " cells=" << cells
                  << " handelman_cells=" << handelman_cells
                  << " real_cells=" << real_cells
                  << " nia_cells=" << nia_cells
                  << " residue_queries=" << residue_queries
                  << " maximum_scale=" << maximum_scale
                  << " result=unsat\n";
        return EXIT_SUCCESS;
    } catch (const z3::exception& error) {
        std::cerr << "z3 error: " << error.msg() << '\n';
        return EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
