#include <z3++.h>

#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Signature = std::array<int, 6>;

void generate_signatures(
    const std::size_t position,
    const int minimum_branch,
    Signature& signature,
    std::vector<Signature>& signatures) {
    if (position == signature.size()) {
        signatures.push_back(signature);
        return;
    }
    for (int branch = minimum_branch; branch <= 3; ++branch) {
        signature[position] = branch;
        generate_signatures(
            position + 1U, branch, signature, signatures);
    }
}

z3::expr twice_coefficient(
    const bool late_wall,
    const int branch,
    const z3::expr& a,
    const z3::expr& shift,
    const z3::expr& label) {
    const z3::expr width = 2 * a;
    const z3::expr size = width + 1;
    const z3::expr band
        = 2 * size * (2 * label + 1) - 2 * label * (label + 1);
    const z3::expr wall_offset = label - 2 * shift;
    if (late_wall) {
        if (branch == 0) {
            return band;
        }
        if (branch == 1) {
            return 2 * size * size;
        }
        if (branch == 2) {
            return 2 * size * size
                   - wall_offset * (wall_offset + 1);
        }
    } else {
        if (branch == 0) {
            return band;
        }
        if (branch == 1) {
            return band - wall_offset * (wall_offset + 1);
        }
        if (branch == 2) {
            return 2 * size * size
                   - wall_offset * (wall_offset + 1);
        }
    }
    const z3::expr distance = 2 * shift + 2 * width - label;
    return (distance + 1) * (distance + 2);
}

void add_branch_constraints(
    z3::solver& solver,
    const bool late_wall,
    const int branch,
    const z3::expr& a,
    const z3::expr& shift,
    const z3::expr& label) {
    const z3::expr width = 2 * a;
    if (late_wall) {
        if (branch == 0) {
            solver.add(label <= width);
        } else if (branch == 1) {
            solver.add(label >= width + 1);
            solver.add(label <= 2 * shift);
        } else if (branch == 2) {
            solver.add(label >= 2 * shift + 1);
            solver.add(label <= 2 * shift + width);
        } else {
            solver.add(label >= 2 * shift + width + 1);
            solver.add(label <= 2 * shift + 2 * width);
        }
    } else {
        if (branch == 0) {
            solver.add(label <= 2 * shift);
        } else if (branch == 1) {
            solver.add(label >= 2 * shift + 1);
            solver.add(label <= width);
        } else if (branch == 2) {
            solver.add(label >= width + 1);
            solver.add(label <= 2 * shift + width);
        } else {
            solver.add(label >= 2 * shift + width + 1);
            solver.add(label <= 2 * shift + 2 * width);
        }
    }
}

std::string render(const Signature& signature) {
    std::string result;
    for (const int branch : signature) {
        result += static_cast<char>('0' + branch);
    }
    return result;
}

}  // namespace

int main() {
    try {
        Signature partial{};
        std::vector<Signature> signatures;
        generate_signatures(0U, 0, partial, signatures);
        if (signatures.size() != 84U) {
            throw std::runtime_error("unexpected branch-signature count");
        }

        std::size_t attempted = 0U;
        std::size_t certified = 0U;
        for (const bool late_wall : {false, true}) {
            for (const Signature& signature : signatures) {
                ++attempted;
                z3::context context;
                z3::solver solver(context);
                z3::params parameters(context);
                parameters.set("timeout", 5000U);
                solver.set(parameters);

                const z3::expr a = context.real_const("a");
                const z3::expr shift = context.real_const("u");
                const z3::expr antidiagonal = context.real_const("A");
                const z3::expr contraction = context.real_const("L");
                const z3::expr width = 2 * a;
                const std::array<z3::expr, 6> labels{
                    contraction,
                    contraction + 1,
                    antidiagonal - contraction,
                    antidiagonal - contraction + 1,
                    antidiagonal + 1,
                    antidiagonal + 2,
                };

                solver.add(a >= 1);
                solver.add(shift >= 1);
                solver.add(shift <= width - 1);
                solver.add(antidiagonal >= 1);
                solver.add(contraction >= 0);
                solver.add(2 * contraction + 1 <= antidiagonal);
                solver.add(
                    antidiagonal <= 2 * (shift + width) - 2);
                if (late_wall) {
                    solver.add(shift >= a);
                } else {
                    solver.add(shift <= a - 1);
                }
                for (std::size_t index = 0U;
                     index < labels.size();
                     ++index) {
                    add_branch_constraints(
                        solver,
                        late_wall,
                        signature[index],
                        a,
                        shift,
                        labels[index]);
                }

                const z3::expr c_zero = 2 * (width + 1);
                std::array<z3::expr, 6> coefficients{
                    twice_coefficient(
                        late_wall,
                        signature[0],
                        a,
                        shift,
                        labels[0]),
                    twice_coefficient(
                        late_wall,
                        signature[1],
                        a,
                        shift,
                        labels[1]),
                    twice_coefficient(
                        late_wall,
                        signature[2],
                        a,
                        shift,
                        labels[2]),
                    twice_coefficient(
                        late_wall,
                        signature[3],
                        a,
                        shift,
                        labels[3]),
                    twice_coefficient(
                        late_wall,
                        signature[4],
                        a,
                        shift,
                        labels[4]),
                    twice_coefficient(
                        late_wall,
                        signature[5],
                        a,
                        shift,
                        labels[5]),
                };
                const z3::expr four_radial
                    = c_zero * (coefficients[4] + coefficients[5])
                      + coefficients[0] * coefficients[2]
                      - coefficients[1] * coefficients[3];
                solver.add(four_radial < 0);

                const z3::check_result result = solver.check();
                if (result == z3::unsat) {
                    ++certified;
                    continue;
                }
                std::cout
                    << "SU2_TWO_FACTOR_INTERMEDIATE_RADIAL_CELL"
                    << " regime="
                    << (late_wall ? "n_le_2u" : "2u_lt_n")
                    << " signature=" << render(signature)
                    << " result=";
                if (result == z3::sat) {
                    const z3::model model = solver.get_model();
                    std::cout
                        << "REAL_COUNTEREXAMPLE"
                        << " a=" << model.eval(a, true)
                        << " u=" << model.eval(shift, true)
                        << " A=" << model.eval(antidiagonal, true)
                        << " L=" << model.eval(contraction, true)
                        << " four_radial="
                        << model.eval(four_radial, true)
                        << '\n';
                    return EXIT_FAILURE;
                }
                std::cout
                    << "UNKNOWN reason=" << solver.reason_unknown()
                    << '\n';
                return 2;
            }
        }

        std::cout
            << "SU2_TWO_FACTOR_INTERMEDIATE_RADIAL_CELLS"
            << " signatures_per_regime=" << signatures.size()
            << " attempted=" << attempted
            << " certified=" << certified
            << " result="
            << (certified == attempted
                    ? "PASS_UNSAT_EXACT_REAL_CELLS"
                    : "FAIL")
            << '\n';
        return certified == attempted ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const z3::exception& error) {
        std::cerr << "z3 error: " << error.msg() << '\n';
        return 3;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 4;
    }
}
