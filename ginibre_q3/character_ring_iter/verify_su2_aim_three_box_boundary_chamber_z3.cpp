#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

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

z3::expr twice_quadratic_spline(const z3::expr& degree) {
    return z3::ite(
        degree >= 0,
        (degree + 1) * (degree + 2),
        degree.ctx().int_val(0));
}

z3::expr twice_three_box_coefficient(
    const z3::expr& shortest,
    const z3::expr& middle,
    const z3::expr& longest,
    const z3::expr& support_end,
    const z3::expr& index) {
    const z3::expr reflected = z3::ite(
        2 * index <= support_end,
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
        index >= 1 && index < support_end,
        half_coefficient,
        index.ctx().int_val(0));
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

int main(int argc, char** argv) {
    try {
        if (argc != 6) {
            throw std::invalid_argument(
                "usage: verify_su2_aim_three_box_boundary_chamber_z3 "
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
        z3::solver solver(context, "QF_NIA");
        z3::params parameters(context);
        parameters.set("timeout", 600000U);
        solver.set(parameters);

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

        solver.add(half_quotient >= 1);
        solver.add(parity == parity_case);
        solver.add(half_level >= 2);
        solver.add(left >= 1 && left < period);
        solver.add(width >= 1 && width < period);
        solver.add(padding >= 0 && padding <= half_level);
        solver.add(label >= 1 && label <= maximum_label);
        solver.add(support_end > 2 * period);
        if (support_band == 0) {
            solver.add(support_end <= 3 * period);
        } else {
            solver.add(support_end > 3 * period);
        }
        if (width_order == 0) {
            solver.add(width <= left);
            shortest = width;
            middle = left;
            longest = shift;
        } else if (width_order == 1) {
            solver.add(width > left && width <= shift);
            shortest = left;
            middle = width;
            longest = shift;
        } else {
            solver.add(width > shift);
            shortest = left;
            middle = shift;
            longest = width;
        }
        if (shift_band == 0) {
            solver.add(shift < period);
        } else {
            solver.add(shift >= period);
        }
        if (pair_order == 0) {
            solver.add(longest <= shortest + middle);
        } else {
            solver.add(longest > shortest + middle);
        }

        z3::expr twice_margin = twice_three_box_coefficient(
            shortest,
            middle,
            longest,
            support_end,
            period - label - 1);
        for (int wall = 1; wall <= 3; ++wall) {
            twice_margin = twice_margin
                - twice_three_box_coefficient(
                      shortest,
                      middle,
                      longest,
                      support_end,
                      wall * period + label)
                + twice_three_box_coefficient(
                      shortest,
                      middle,
                      longest,
                      support_end,
                      (wall + 1) * period - label - 1);
        }
        solver.add(twice_margin < 0);

        const z3::check_result result = solver.check();
        std::cout << "SU2_AIM_THREE_BOX_BOUNDARY_CHAMBER_Z3"
                  << " logic=QF_NIA"
                  << " support_band=" << support_band
                  << " width_order=" << width_order
                  << " shift_band=" << shift_band
                  << " pair_order=" << pair_order
                  << " parity=" << parity_case
                  << " result=" << render(result)
                  << '\n';
        if (result == z3::sat) {
            std::cout << "MODEL " << solver.get_model() << '\n';
        }
        if (result == z3::unknown) {
            std::cout << "REASON " << solver.reason_unknown() << '\n';
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
