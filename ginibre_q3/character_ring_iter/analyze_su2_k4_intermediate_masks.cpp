#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <z3++.h>

namespace {

void append_block(
    std::vector<z3::expr>& slacks,
    const z3::expr& q,
    const z3::expr& image,
    int power,
    int maximum_image
) {
    for (int index = 0; index <= maximum_image; ++index) {
        slacks.push_back(
            (power - 2 * index) * q - image - index
        );
    }
}

bool model_boolean(
    const z3::model& model,
    const z3::expr& expression,
    z3::context& context
) {
    return z3::eq(
        model.eval(expression, true),
        context.bool_val(true)
    );
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const bool list_masks =
            argc == 2 && std::string(argv[1]) == "--list";
        if (argc != 1 && !list_masks) {
            throw std::runtime_error("usage: [--list]");
        }
        z3::context context;
        z3::solver solver(context);
        const z3::expr q = context.int_const("Q");
        const z3::expr h = context.int_const("H");
        const z3::expr y = context.int_const("Y");
        const z3::expr zero = context.int_val(0);
        const z3::expr period = 4 * q + 4 + 2 * h;

        solver.add(q >= 1);
        solver.add(h >= 0);
        solver.add(y >= 0);
        solver.add(y <= 2 * q + 1 + h);

        std::vector<z3::expr> slacks;
        const auto append_eight = [&](const z3::expr& endpoint) {
            append_block(slacks, q, endpoint, 8, 3);
            append_block(
                slacks,
                q,
                period - 1 - endpoint,
                8,
                2
            );
            append_block(slacks, q, period + endpoint, 8, 1);
            append_block(
                slacks,
                q,
                2 * period - 1 - endpoint,
                8,
                0
            );
        };
        const auto append_nine = [&](const z3::expr& endpoint) {
            append_block(slacks, q, endpoint, 9, 4);
            append_block(
                slacks,
                q,
                period - 1 - endpoint,
                9,
                3
            );
            append_block(slacks, q, period + endpoint, 9, 2);
            append_block(
                slacks,
                q,
                2 * period - 1 - endpoint,
                9,
                1
            );
            append_block(slacks, q, 2 * period + endpoint, 9, 0);
        };

        append_eight(y);
        append_nine(y);
        append_eight(zero);
        append_nine(zero);
        if (slacks.size() != 50U) {
            throw std::runtime_error(
                "unexpected number of C4 activation predicates"
            );
        }

        std::set<std::string> masks;
        std::size_t minimum_active = slacks.size();
        std::size_t maximum_active = 0U;
        while (solver.check() == z3::sat) {
            const z3::model model = solver.get_model();
            std::string mask;
            mask.reserve(slacks.size());
            std::size_t active_count = 0U;
            z3::expr block = context.bool_val(false);
            for (const z3::expr& slack : slacks) {
                const z3::expr active_expression = slack >= 0;
                const bool active = model_boolean(
                    model,
                    active_expression,
                    context
                );
                mask.push_back(active ? '1' : '0');
                active_count += active ? 1U : 0U;
                block = block
                    || (active ? !active_expression : active_expression);
            }
            if (!masks.insert(mask).second) {
                throw std::runtime_error(
                    "blocking clause repeated a C4 activation mask"
                );
            }
            minimum_active = std::min(minimum_active, active_count);
            maximum_active = std::max(maximum_active, active_count);
            solver.add(block);
            if (masks.size() % 1000U == 0U) {
                std::cerr
                    << "SU2_K4_INTERMEDIATE_MASKS"
                    << " progress=" << masks.size()
                    << std::endl;
            }
        }

        std::cout
            << "SU2_K4_INTERMEDIATE_MASKS"
            << " hinges=" << slacks.size()
            << " masks=" << masks.size()
            << " minimum_active=" << minimum_active
            << " maximum_active=" << maximum_active
            << " result=PASS_EXACT_CENSUS\n";
        if (list_masks) {
            std::vector<std::uint64_t> numeric_masks;
            numeric_masks.reserve(masks.size());
            for (const std::string& mask : masks) {
                std::uint64_t numeric = 0U;
                for (std::size_t index = 0U;
                     index < mask.size();
                     ++index) {
                    if (mask[index] == '1') {
                        numeric |= std::uint64_t{1} << index;
                    }
                }
                numeric_masks.push_back(numeric);
            }
            std::sort(numeric_masks.begin(), numeric_masks.end());
            for (const std::uint64_t mask : numeric_masks) {
                std::cout
                    << "SU2_K4_INTERMEDIATE_MASK value="
                    << mask << '\n';
            }
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_K4_INTERMEDIATE_MASKS FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
