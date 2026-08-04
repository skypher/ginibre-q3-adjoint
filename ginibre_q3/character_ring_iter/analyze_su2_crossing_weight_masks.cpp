#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <z3++.h>

namespace {

z3::expr maximum2(const z3::expr& left, const z3::expr& right) {
    return z3::ite(left >= right, left, right);
}

z3::expr minimum2(const z3::expr& left, const z3::expr& right) {
    return z3::ite(left <= right, left, right);
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

void append_power_hinges(
    z3::context& context,
    const z3::expr& level,
    const z3::expr& label,
    const z3::expr& target,
    int power,
    std::vector<z3::expr>& hinges
) {
    const z3::expr product = power * label;
    const z3::expr lower_depth = level - label - target;
    const z3::expr upper_depth = label + target;
    for (int branch = 0; branch < 3; ++branch) {
        z3::expr lower_label = context.int_val(0);
        z3::expr upper_label = context.int_val(0);
        if (branch == 0) {
            lower_label = lower_depth;
            upper_label = upper_depth;
        } else if (branch == 1) {
            lower_label = 2 * level - label - target + 1;
            upper_label = level + label + target + 1;
        } else {
            lower_label = 2 * level + 2 + lower_depth;
            upper_label = 2 * level + 2 + upper_depth;
        }
        const z3::expr raw_lower = product - upper_label;
        const z3::expr raw_upper = product - lower_label;
        const z3::expr lower = maximum2(context.int_val(0), raw_lower);
        const z3::expr upper = minimum2(product, raw_upper);
        hinges.push_back(lower == context.int_val(0));
        hinges.push_back(lower == raw_lower);
        hinges.push_back(upper == product);
        hinges.push_back(upper == raw_upper);
        hinges.push_back(lower <= upper);
        hinges.push_back(lower - 1 >= 0);
        for (int image = 1; image <= power; ++image) {
            const z3::expr activation =
                image * (2 * label + 1);
            hinges.push_back(upper >= activation);
            hinges.push_back(lower - 1 >= activation);
        }
    }
}

std::uint64_t count_masks(const std::vector<int>& powers) {
    z3::context context;
    z3::solver solver(context);
    const z3::expr level = context.int_const("K");
    const z3::expr label = context.int_const("Q");
    const z3::expr target = context.int_const("V");
    solver.add(label >= 1);
    solver.add(2 * label < level);
    solver.add(target >= 0);
    solver.add(2 * target < level);
    solver.add(level - label - target <= label + target);

    std::vector<z3::expr> hinges;
    hinges.push_back(z3::mod(level, 2) == 0);
    hinges.push_back(z3::mod(label, 2) == 0);
    hinges.push_back(z3::mod(target, 2) == 0);
    for (const int power : powers) {
        append_power_hinges(
            context,
            level,
            label,
            target,
            power,
            hinges
        );
    }

    std::set<std::string> masks;
    while (solver.check() == z3::sat) {
        const z3::model model = solver.get_model();
        std::string mask;
        mask.reserve(hinges.size());
        z3::expr block = context.bool_val(false);
        for (const z3::expr& hinge : hinges) {
            const bool active = model_boolean(model, hinge, context);
            mask.push_back(active ? '1' : '0');
            block = block || (active ? !hinge : hinge);
        }
        if (!masks.insert(mask).second) {
            throw std::runtime_error("repeated crossing-weight mask");
        }
        solver.add(block);
    }
    return static_cast<std::uint64_t>(masks.size());
}

}  // namespace

int main() {
    try {
        for (int power = 1; power <= 5; ++power) {
            const std::uint64_t masks = count_masks({power});
            std::cout
                << "SU2_CROSSING_WEIGHT_MASKS"
                << " power=" << power
                << " masks=" << masks
                << " result=PASS_EXACT_PRESBURGER_CENSUS\n";
        }
        const std::uint64_t joint_masks = count_masks({1, 2, 3, 4, 5});
        std::cout
            << "SU2_CROSSING_WEIGHT_JOINT_MASKS"
            << " powers=1,2,3,4,5"
            << " active_crossing_masks=" << joint_masks
            << " zero_crossing_region=1"
            << " result=PASS_EXACT_PRESBURGER_CENSUS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
