#include <cstdlib>
#include <iostream>

#include <z3++.h>

int main() {
    z3::context context;
    z3::solver solver(context, "QF_NRA");
    z3::params parameters(context);
    parameters.set("timeout", 300000U);
    solver.set(parameters);
    const z3::expr m = context.real_const("m");
    const z3::expr padding = context.real_const("A");
    const z3::expr load = context.real_const("a");
    const z3::expr width_slack = context.real_const("u");
    const z3::expr parity_slack = context.real_const("v");
    const z3::expr width = 2 * padding + parity_slack;
    const z3::expr left = width + width_slack;
    const z3::expr label
        = 4 * m + 1 + load - width;

    const z3::expr margin
        = 6 + 36 * m + 48 * m * m
          - 2 * left - 7 * width + width * width
          + 10 * padding - 4 * padding * padding
          - 11 * label + label * label
          - 8 * m * left - 16 * m * width
          + 16 * m * padding - 24 * m * label
          + 2 * left * width - 4 * left * padding
          + 4 * left * label + 2 * width * label
          + 4 * padding * label;

    solver.add(4 * m + 1 - left >= 0);
    solver.add(
        12 * m + 7 - 2 * left - width - 2 * padding >= 0);
    solver.add(left - width >= 0);
    solver.add(width - 2 * padding >= 0);
    solver.add(-1 - 4 * m + label + width >= 0);
    solver.add(
        16 * m + 6 - 2 * label - 2 * left
        - width - 2 * padding >= 0);
    solver.add(4 - left - width + 8 * m - label >= 0);
    solver.add(
        -5 + 2 * left + width + 2 * padding
        - 8 * m - 2 * label >= 0);
    solver.add(
        -2 - 4 * m - label + left + 2 * padding >= 0);
    solver.add(
        -7 + 2 * left + width + 2 * padding
        - 12 * m + label >= 0);
    solver.add(
        5 - 2 * left - 2 * padding + 8 * m + label >= 0);
    solver.add(margin <= -2);

    z3::params normalization(context);
    normalization.set("som", true);
    std::cout << "MARGIN " << margin.simplify(normalization) << '\n';
    const z3::check_result result = solver.check();
    std::cout << "SU2_AIM_THREE_BOX_BOUNDARY_SLACKS_Z3 result="
              << (result == z3::unsat
                      ? "unsat"
                      : result == z3::sat ? "sat" : "unknown")
              << '\n';
    if (result == z3::sat) {
        std::cout << solver.get_model() << '\n';
    } else if (result == z3::unknown) {
        std::cout << solver.reason_unknown() << '\n';
    }
    return result == z3::unsat ? EXIT_SUCCESS : EXIT_FAILURE;
}
