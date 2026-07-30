#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

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

struct QueryResult {
    z3::check_result result;
    std::string model;
};

QueryResult run_query(int shell, bool middle) {
    z3::context context;
    z3::solver solver(context, "QF_LIA");
    const z3::expr half_level = context.int_const("K");
    const z3::expr left = context.int_const("ell");
    const z3::expr right = context.int_const("w");
    const z3::expr padding = context.int_const("A");
    const z3::expr label = context.int_const("s");
    const z3::expr period = 2 * half_level + 2;
    const z3::expr shift = left + 2 * padding;

    solver.add(half_level >= 2);
    solver.add(left >= 1 && left < period);
    solver.add(right >= 1 && right < period);
    solver.add(padding >= 0 && padding <= half_level);
    if (middle) {
        solver.add(label == half_level);
    } else {
        solver.add(label >= 0 && label < half_level);
    }

    z3::expr margin = context.int_val(0);
    for (int wall = shell; wall <= 3; ++wall) {
        const z3::expr first = wall * period + label;
        margin = margin
            + triangular_count(left, right, first - shift)
            - triangular_count(left, right, first);
        if (!middle) {
            const z3::expr reflected
                = (wall + 1) * period - label - 2;
            margin = margin
                + triangular_count(
                    left,
                    right,
                    reflected - shift)
                - triangular_count(left, right, reflected);
        }
    }
    solver.add(margin < 0);
    const z3::check_result result = solver.check();
    return {
        result,
        result == z3::sat ? solver.get_model().to_string() : ""};
}

const char* render_result(const z3::check_result& result) {
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
        const QueryResult first_adjacent = run_query(1, false);
        const QueryResult second_adjacent = run_query(2, false);
        const QueryResult first_middle = run_query(1, true);
        const QueryResult second_middle = run_query(2, true);

        std::cout
            << "SU2_AIM_THREE_BOX_MONOTONICITY_Z3"
            << " logic=QF_LIA"
            << " first_adjacent="
            << render_result(first_adjacent.result)
            << " second_adjacent="
            << render_result(second_adjacent.result)
            << " first_middle="
            << render_result(first_middle.result)
            << " second_middle="
            << render_result(second_middle.result)
            << '\n';
        const auto print_model = [](const char* name,
                                    const QueryResult& query) {
            if (query.result == z3::sat) {
                std::cout << name << "_MODEL " << query.model << '\n';
            }
        };
        print_model("FIRST_ADJACENT", first_adjacent);
        print_model("SECOND_ADJACENT", second_adjacent);
        print_model("FIRST_MIDDLE", first_middle);
        print_model("SECOND_MIDDLE", second_middle);

        return first_adjacent.result == z3::unsat
                && second_adjacent.result == z3::unsat
                && first_middle.result == z3::unsat
                && second_middle.result == z3::unsat
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const z3::exception& error) {
        std::cerr << "z3 error: " << error.msg() << '\n';
        return EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
