#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <z3++.h>

namespace {

z3::expr maximum(const z3::expr& left, const z3::expr& right) {
    return z3::ite(left >= right, left, right);
}

z3::expr minimum(const z3::expr& left, const z3::expr& right) {
    return z3::ite(left <= right, left, right);
}

z3::expr triangular_count(
    const z3::expr& left_length,
    const z3::expr& right_length,
    const z3::expr& total) {
    const z3::expr begin = maximum(
        total - right_length + 1,
        total.ctx().int_val(0));
    const z3::expr end = minimum(left_length - 1, total);
    return maximum(
        end - begin + 1,
        total.ctx().int_val(0));
}

z3::expr adjacent_margin_difference(
    const z3::expr& period,
    const z3::expr& left,
    const z3::expr& width,
    const z3::expr& shift,
    const z3::expr& label) {
    z3::expr difference = triangular_count(
                              left,
                              width,
                              period - label - 2 - shift)
        - triangular_count(
                              left,
                              width,
                              period - label - 2);
    for (int wall = 1; wall <= 3; ++wall) {
        difference = difference
            + triangular_count(
                  left,
                  width,
                  wall * period + label - shift)
            - triangular_count(
                  left,
                  width,
                  wall * period + label)
            + triangular_count(
                  left,
                  width,
                  (wall + 1) * period - label - 2 - shift)
            - triangular_count(
                  left,
                  width,
                  (wall + 1) * period - label - 2);
    }
    return difference;
}

const char* render(const z3::check_result& result) {
    if (result == z3::unsat) {
        return "unsat";
    }
    if (result == z3::sat) {
        return "sat";
    }
    return "unknown";
}

}  // namespace

int main() {
    try {
        z3::context context;
        z3::solver solver(context, "QF_LIA");
        const z3::expr half_quotient = context.int_const("m");
        const z3::expr parity = context.int_const("p");
        const z3::expr half_level
            = 2 * half_quotient + parity;
        const z3::expr period = 2 * half_level + 2;
        const z3::expr left = context.int_const("ell");
        const z3::expr width = context.int_const("w");
        const z3::expr padding = context.int_const("A");
        const z3::expr shift = left + 2 * padding;
        const z3::expr first_label = context.int_const("s");
        const z3::expr second_label = context.int_const("t");
        const z3::expr maximum_label = half_level - parity;

        solver.add(half_quotient >= 1);
        solver.add(parity >= 0 && parity <= 1);
        solver.add(half_level >= 2);
        solver.add(left >= 1 && left < period);
        solver.add(width >= 1 && width < period);
        solver.add(padding >= 0 && padding <= half_level);
        solver.add(left + width + shift - 1 > 2 * period);
        solver.add(first_label >= 1);
        solver.add(first_label < second_label);
        solver.add(second_label < maximum_label);
        solver.add(adjacent_margin_difference(
                       period,
                       left,
                       width,
                       shift,
                       first_label)
                   < 0);
        solver.add(adjacent_margin_difference(
                       period,
                       left,
                       width,
                       shift,
                       second_label)
                   > 0);

        const z3::check_result result = solver.check();
        std::cout << "SU2_AIM_THREE_BOX_LABEL_ENDPOINTS_Z3"
                  << " logic=QF_LIA"
                  << " forbidden_recrossing=" << render(result)
                  << '\n';
        if (result == z3::sat) {
            std::cout << "MODEL " << solver.get_model() << '\n';
        }
        return result == z3::unsat ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const z3::exception& error) {
        std::cerr << "z3 error: " << error.msg() << '\n';
        return EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
