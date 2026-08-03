#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#define SU2_T4_GROUPS_EMBEDDED
#include "prove_su2_t4_groups.cpp"
#undef SU2_T4_GROUPS_EMBEDDED

namespace {

std::size_t parse_position(const char* text) {
    const std::string value{text};
    std::size_t consumed = 0U;
    const unsigned long long parsed = std::stoull(value, &consumed, 10);
    if (consumed != value.size()) {
        throw std::invalid_argument("position must be a nonnegative integer");
    }
    if (parsed > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument("position is out of range");
    }
    return static_cast<std::size_t>(parsed);
}

std::string render(const Rational& value) {
    if (value.denominator() == 1) {
        return value.numerator().convert_to<std::string>();
    }
    return value.numerator().convert_to<std::string>() + "/"
           + value.denominator().convert_to<std::string>();
}

void print_affine(const Polynomial& polynomial, const std::string& name) {
    const Affine coefficients = affine_coefficients(polynomial);
    std::cout << name
              << " constant=" << render(coefficients[0])
              << " Q=" << render(coefficients[1])
              << " H=" << render(coefficients[2])
              << " Y=" << render(coefficients[3]) << '\n';
}

void print_model(const std::vector<Polynomial>& constraints) {
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
    if (solver.check() != z3::sat) {
        throw std::runtime_error("selected chamber is infeasible");
    }
    const z3::model model = solver.get_model();
    std::cout << "SU2_T4_GROUP_CHAMBER_MODEL"
              << " Q=" << model.eval(variables[0], true)
              << " H=" << model.eval(variables[1], true)
              << " Y=" << model.eval(variables[2], true) << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        static_cast<void>(&certify_group_chamber);
        if (argc != 2 && argc != 3) {
            throw std::invalid_argument(
                "usage: analyze_su2_t4_group_chamber <c5-position>"
                "\n   or: analyze_su2_t4_group_chamber --scan-low-cost"
                "\n   or: analyze_su2_t4_group_chamber --scan-parity"
                "\n   or: analyze_su2_t4_group_chamber --scan-forced-equality"
                "\n   or: analyze_su2_t4_group_chamber --parity-only"
                " <c5-position>"
            );
        }
        const bool scan_low_cost = std::string(argv[1]) == "--scan-low-cost";
        const bool scan_parity = std::string(argv[1]) == "--scan-parity";
        const bool scan_forced_equality =
            std::string(argv[1]) == "--scan-forced-equality";
        const bool parity_only = std::string(argv[1]) == "--parity-only";
        if (parity_only && argc != 3) {
            throw std::invalid_argument(
                "--parity-only requires a c5 position"
            );
        }
        const GroupFormula formula = make_group_formula("c5");
        const std::vector<WideMask> masks = feasible_group_masks(formula);
        if (scan_low_cost) {
            std::size_t certified = 0U;
            for (std::size_t position = 0U; position < masks.size(); ++position) {
                const Chamber chamber = make_group_chamber(
                    formula,
                    masks[position],
                    "c5",
                    static_cast<std::uint64_t>(position)
                );
                const std::vector<Polynomial> reduced =
                    irredundant_constraints(chamber);
                bool passed = bounded_group_integer_certificate(
                    chamber,
                    reduced
                );
                if (!passed) {
                    passed = bounded_qy_h_ray_newton_certificate(
                        chamber,
                        reduced
                    );
                }
                if (passed) {
                    ++certified;
                }
                if ((position + 1U) % 25U == 0U) {
                    std::cerr
                        << "SU2_T4_GROUP_LOW_COST"
                        << " progress=" << position + 1U
                        << '/' << masks.size()
                        << " certified=" << certified << '\n';
                }
            }
            std::cout
                << "SU2_T4_GROUP_LOW_COST"
                << " masks=" << masks.size()
                << " certified=" << certified
                << " unresolved=" << masks.size() - certified
                << " result=PASS_EXACT_PARTITION"
                << '\n';
            return EXIT_SUCCESS;
        }
        if (scan_parity) {
            std::size_t certified = 0U;
            for (std::size_t position = 0U; position < masks.size(); ++position) {
                const Chamber chamber = make_group_chamber(
                    formula,
                    masks[position],
                    "c5",
                    static_cast<std::uint64_t>(position)
                );
                const std::vector<Polynomial> reduced =
                    irredundant_constraints(chamber);
                if (parity_unit_offset_square_cone_certificate(
                        chamber,
                        reduced
                    )) {
                    ++certified;
                }
                if ((position + 1U) % 25U == 0U) {
                    std::cerr
                        << "SU2_T4_GROUP_PARITY_SCAN"
                        << " progress=" << position + 1U
                        << '/' << masks.size()
                        << " certified=" << certified << '\n';
                }
            }
            std::cout
                << "SU2_T4_GROUP_PARITY_SCAN"
                << " masks=" << masks.size()
                << " certified=" << certified
                << " unresolved=" << masks.size() - certified
                << " result=PASS_EXACT_PARTITION"
                << '\n';
            return EXIT_SUCCESS;
        }
        if (scan_forced_equality) {
            std::size_t certified = 0U;
            for (std::size_t position = 0U; position < masks.size(); ++position) {
                const Chamber chamber = make_group_chamber(
                    formula,
                    masks[position],
                    "c5",
                    static_cast<std::uint64_t>(position)
                );
                const std::vector<Polynomial> reduced =
                    irredundant_constraints(chamber);
                if (forced_equality_slice_ray_newton_certificate(
                        chamber,
                        reduced
                    )) {
                    ++certified;
                }
                if ((position + 1U) % 25U == 0U) {
                    std::cerr
                        << "SU2_T4_GROUP_FORCED_EQUALITY_SCAN"
                        << " progress=" << position + 1U
                        << '/' << masks.size()
                        << " certified=" << certified << '\n';
                }
            }
            std::cout
                << "SU2_T4_GROUP_FORCED_EQUALITY_SCAN"
                << " masks=" << masks.size()
                << " certified=" << certified
                << " unresolved=" << masks.size() - certified
                << " result=PASS_EXACT_PARTITION"
                << '\n';
            return EXIT_SUCCESS;
        }
        const std::size_t position = parse_position(
            parity_only ? argv[2] : argv[1]
        );
        if (position >= masks.size()) {
            throw std::invalid_argument("position is out of range");
        }
        const Chamber chamber = make_group_chamber(
            formula, masks[position], "c5", static_cast<std::uint64_t>(position)
        );
        const std::vector<Polynomial> reduced = irredundant_constraints(chamber);
        if (parity_only) {
            const bool passed = parity_unit_offset_square_cone_certificate(
                chamber,
                reduced
            );
            std::cout
                << "SU2_T4_GROUP_PARITY_ONLY"
                << " position=" << position
                << " result="
                << (passed ? "PASS_EXACT" : "NO_CERTIFICATE")
                << '\n';
            return passed ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        std::cout << "SU2_T4_GROUP_CHAMBER"
                  << " target=c5"
                  << " position=" << position
                  << " mask=" << masks[position]
                  << " constraints=" << chamber.constraints.size()
                  << " irredundant=" << reduced.size()
                  << " margin_terms=" << chamber.margin.terms().size() << '\n';
        for (std::size_t index = 0U; index < reduced.size(); ++index) {
            print_affine(
                reduced[index], "SU2_T4_GROUP_CHAMBER_CONSTRAINT"
                    " index=" + std::to_string(index)
            );
        }
        print_model(reduced);
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
