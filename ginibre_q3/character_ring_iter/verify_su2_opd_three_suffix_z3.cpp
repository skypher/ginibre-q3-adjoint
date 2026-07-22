#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include <z3++.h>

namespace {

z3::expr absolute(const z3::expr& value) {
    return z3::ite(value >= 0, value, -value);
}

z3::expr maximum(const z3::expr& first, const z3::expr& second) {
    return z3::ite(first >= second, first, second);
}

z3::expr indicator(const z3::expr& proposition) {
    return z3::ite(proposition, proposition.ctx().int_val(1),
                   proposition.ctx().int_val(0));
}

z3::expr fusion_contains(
    const z3::expr& left,
    const z3::expr& right,
    const z3::expr& output
) {
    return output >= absolute(left - right)
        && output <= left + right
        && z3::mod(left + right - output, 2) == 0;
}

void print_model_value(
    const z3::model& model,
    const char* name,
    const z3::expr& value
) {
    std::cout << ' ' << name << '=' << model.eval(value, true);
}

}  // namespace

int main() {
    try {
        z3::context context;
        z3::solver solver(context);
        const z3::expr p = context.int_const("p");
        const z3::expr q = context.int_const("q");
        const z3::expr r = context.int_const("r");
        const z3::expr s = context.int_const("s");
        const z3::expr a = context.int_const("a");
        solver.add(p >= 2 && p <= q && q <= r && r <= s && a >= 1);

        const z3::expr labels[3]{q, r, s};
        z3::expr negative = context.int_val(0);
        for (int omitted = 0; omitted < 3; ++omitted) {
            const int first_index = (omitted + 1) % 3;
            const int second_index = (omitted + 2) % 3;
            const z3::expr& first = labels[first_index];
            const z3::expr& second = labels[second_index];
            const z3::expr& remaining = labels[omitted];
            negative = negative + indicator(
                a == remaining && fusion_contains(first, second, p - 1)
            );
            negative = negative + indicator(
                fusion_contains(first, second, context.int_val(1))
                && fusion_contains(remaining, p, a)
            );
        }
        const z3::expr triple_to_one
            = indicator(fusion_contains(q, r, s - 1))
            + indicator(fusion_contains(q, r, s + 1));
        negative = negative + z3::ite(
            a == p, triple_to_one, context.int_val(0)
        );

        z3::expr candidates = context.int_val(0);
        const z3::expr first_low = q - p + 1;
        const z3::expr first_high = q + p - 1;
        for (int first_index = 0; first_index < 6; ++first_index) {
            const z3::expr first = first_low + 2 * first_index;
            const z3::expr parity
                = z3::mod(first + r - a - s, 2) == 0;
            const z3::expr second_low = maximum(
                absolute(first - r), absolute(a - s)
            );
            for (int second_index = 0; second_index < 6; ++second_index) {
                const z3::expr second
                    = second_low + 2 * second_index;
                const z3::expr valid = first <= first_high && parity
                    && second <= first + r && second <= a + s;
                candidates = candidates + indicator(valid);
            }
        }

        const z3::expr exceptional = q == p && r == p && s == p + 1
            && a == p && z3::mod(p, 2) == 0;
        solver.add(!exceptional);
        solver.add(candidates < negative);
        const z3::check_result result = solver.check();
        if (result == z3::sat) {
            const z3::model model = solver.get_model();
            std::cout << "SU2_OPD_THREE_SUFFIX_Z3 FAIL";
            print_model_value(model, "p", p);
            print_model_value(model, "q", q);
            print_model_value(model, "r", r);
            print_model_value(model, "s", s);
            print_model_value(model, "a", a);
            print_model_value(model, "candidates", candidates);
            print_model_value(model, "negative", negative);
            std::cout << '\n';
            return EXIT_FAILURE;
        }
        if (result == z3::unknown) {
            std::cout << "SU2_OPD_THREE_SUFFIX_Z3 UNKNOWN reason="
                      << solver.reason_unknown() << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "SU2_OPD_THREE_SUFFIX_Z3 PASS result=unsat "
                  << "candidate_channels=36\n";
        return EXIT_SUCCESS;
    } catch (const z3::exception& error) {
        std::cerr << error.msg() << '\n';
        return EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
