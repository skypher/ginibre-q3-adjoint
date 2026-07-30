#include <z3++.h>

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

z3::expr twice_coefficient(
    z3::context& context,
    const z3::expr& shift,
    const z3::expr& width,
    const z3::expr& label) {
    const z3::expr size = width + 1;
    const z3::expr twice_band = z3::ite(
        label <= width,
        2 * size * (2 * label + 1) - 2 * label * (label + 1),
        2 * size * size);
    const z3::expr triangle_index = label - 2 * shift - 1;
    const z3::expr twice_triangle = z3::ite(
        triangle_index < 0,
        context.int_val(0),
        z3::ite(
            triangle_index <= width,
            (triangle_index + 1) * (triangle_index + 2),
            z3::ite(
                triangle_index <= 2 * width,
                2 * size * size
                    - (2 * width - triangle_index)
                        * (2 * width - triangle_index + 1),
                2 * size * size)));
    return twice_band - twice_triangle;
}

}  // namespace

int main() {
    try {
        z3::context context;
        z3::solver solver(context);
        z3::params parameters(context);
        parameters.set("timeout", 30000U);
        solver.set(parameters);

        const z3::expr a = context.int_const("a");
        const z3::expr shift = context.int_const("u");
        const z3::expr antidiagonal = context.int_const("A");
        const z3::expr contraction = context.int_const("L");
        const z3::expr width = 2 * a;
        const z3::expr middle = antidiagonal - contraction;

        solver.add(a >= 1);
        solver.add(shift >= 1);
        solver.add(shift <= width - 1);
        solver.add(antidiagonal >= 1);
        solver.add(contraction >= 0);
        solver.add(2 * contraction < antidiagonal);
        solver.add(antidiagonal <= 2 * (shift + width) - 2);

        const z3::expr c_zero = twice_coefficient(
            context, shift, width, context.int_val(0));
        const z3::expr c_l = twice_coefficient(
            context, shift, width, contraction);
        const z3::expr c_l_next = twice_coefficient(
            context, shift, width, contraction + 1);
        const z3::expr c_middle = twice_coefficient(
            context, shift, width, middle);
        const z3::expr c_middle_next = twice_coefficient(
            context, shift, width, middle + 1);
        const z3::expr c_top = twice_coefficient(
            context, shift, width, antidiagonal + 1);
        const z3::expr c_top_next = twice_coefficient(
            context, shift, width, antidiagonal + 2);
        const z3::expr four_radial
            = c_zero * (c_top + c_top_next)
              + c_l * c_middle - c_l_next * c_middle_next;
        solver.add(four_radial < 0);

        const z3::check_result result = solver.check();
        std::cout
            << "SU2_TWO_FACTOR_INTERMEDIATE_RADIAL_Z3"
            << " domain=1<=u<=2a-1"
            << " query=four_radial<0"
            << " result=";
        if (result == z3::unsat) {
            std::cout << "PASS_UNSAT_EXACT\n";
            return EXIT_SUCCESS;
        }
        if (result == z3::sat) {
            const z3::model model = solver.get_model();
            std::cout
                << "COUNTEREXAMPLE"
                << " a=" << model.eval(a, true)
                << " u=" << model.eval(shift, true)
                << " A=" << model.eval(antidiagonal, true)
                << " L=" << model.eval(contraction, true)
                << " four_radial=" << model.eval(four_radial, true)
                << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "UNKNOWN reason=" << solver.reason_unknown() << '\n';
        return 2;
    } catch (const z3::exception& error) {
        std::cerr << "z3 error: " << error.msg() << '\n';
        return 3;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 4;
    }
}
