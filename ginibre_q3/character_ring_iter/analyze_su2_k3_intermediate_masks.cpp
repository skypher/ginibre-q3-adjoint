#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

#include <z3++.h>

int main() {
    try {
        z3::context context;
        z3::solver solver(context);
        const z3::expr q = context.int_const("Q");
        const z3::expr h = context.int_const("h");
        const z3::expr y = context.int_const("Y");
        solver.add(q >= 1);
        solver.add(h >= 1);
        solver.add(h <= 5 * q - 2);
        solver.add(y >= 0);
        solver.add(y <= 2 * q + 1 + h);

        const std::vector<z3::expr> slacks{
            6 * q - y,
            4 * q - y - 1,
            2 * q - y - 2,
            2 * q - 2 * h + y - 3,
            y - 2 * h - 4,
            2 * q - 2 * h - y - 4,
            7 * q - y,
            5 * q - y - 1,
            3 * q - y - 2,
            q - y - 3,
            3 * q - 2 * h + y - 3,
            q - 2 * h + y - 4,
            y - q - 2 * h - 5,
            3 * q - 2 * h - y - 4,
            q - 2 * h - y - 5,
            y - q - 4 * h - 7,
            2 * q - 2 * h - 3,
            3 * q - 2 * h - 3,
            q - 2 * h - 4,
            3 * q - 2 * h - 4,
            q - 2 * h - 5,
            q - 3
        };
        if (slacks.size() > 63U) {
            throw std::runtime_error("activation mask exceeds uint64");
        }

        std::set<std::uint64_t> masks;
        std::uint64_t iterations = 0U;
        std::size_t minimum_facets = 1000U;
        std::size_t maximum_facets = 0U;
        while (solver.check() == z3::sat) {
            const z3::model model = solver.get_model();
            std::uint64_t mask = 0U;
            z3::expr block = context.bool_val(false);
            for (std::size_t index = 0U;
                 index < slacks.size();
                 ++index) {
                const z3::expr active_expression =
                    slacks[index] >= 0;
                const bool active = z3::eq(
                    model.eval(active_expression, true),
                    context.bool_val(true)
                );
                if (active) {
                    mask |= std::uint64_t{1} << index;
                    block = block || !active_expression;
                } else {
                    block = block || active_expression;
                }
            }
            if (!masks.insert(mask).second) {
                throw std::runtime_error(
                    "blocking clause repeated an activation mask"
                );
            }
            ++iterations;
            std::vector<z3::expr> constraints{
                q - 1,
                h - 1,
                5 * q - 2 - h,
                y,
                2 * q + 1 + h - y
            };
            for (std::size_t index = 0U;
                 index < slacks.size();
                 ++index) {
                constraints.push_back(
                    (mask & (std::uint64_t{1} << index)) != 0U
                        ? slacks[index]
                        : -slacks[index] - 1
                );
            }
            std::vector<std::size_t> retained(constraints.size());
            std::iota(retained.begin(), retained.end(), 0U);
            for (std::size_t position = 0U;
                 position < retained.size();) {
                z3::solver redundancy_solver(context);
                for (std::size_t other = 0U;
                     other < retained.size();
                     ++other) {
                    if (other != position) {
                        redundancy_solver.add(
                            constraints[retained[other]] >= 0
                        );
                    }
                }
                redundancy_solver.add(
                    constraints[retained[position]] < 0
                );
                if (redundancy_solver.check() == z3::unsat) {
                    retained.erase(
                        retained.begin()
                            + static_cast<std::ptrdiff_t>(position)
                    );
                } else {
                    ++position;
                }
            }
            minimum_facets = std::min(
                minimum_facets,
                retained.size()
            );
            maximum_facets = std::max(
                maximum_facets,
                retained.size()
            );
            std::cout
                << "SU2_K3_INTERMEDIATE_MASK"
                << " index=" << iterations
                << " mask=" << mask
                << " facets=" << retained.size()
                << " witness=("
                << model.eval(q, true) << ','
                << model.eval(h, true) << ','
                << model.eval(y, true) << ")\n";
            solver.add(block);
        }
        std::cout
            << "SU2_K3_INTERMEDIATE_MASKS"
            << " hinges=" << slacks.size()
            << " masks=" << masks.size()
            << " minimum_facets=" << minimum_facets
            << " maximum_facets=" << maximum_facets
            << " result=PASS_EXACT_CENSUS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_K3_INTERMEDIATE_MASKS FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
