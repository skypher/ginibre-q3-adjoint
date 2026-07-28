#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <z3++.h>

namespace {

z3::expr maximum3(
    const z3::expr& first,
    const z3::expr& second,
    const z3::expr& third
) {
    const z3::expr first_two = z3::ite(
        first >= second,
        first,
        second
    );
    return z3::ite(first_two >= third, first_two, third);
}

z3::expr minimum3(
    const z3::expr& first,
    const z3::expr& second,
    const z3::expr& third
) {
    const z3::expr first_two = z3::ite(
        first <= second,
        first,
        second
    );
    return z3::ite(first_two <= third, first_two, third);
}

z3::expr ceiling_half(const z3::expr& value) {
    return -((-value) / 2);
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

int main() {
    try {
        std::uint64_t total_masks = 0U;
        for (int power = 1; power <= 4; ++power) {
            z3::context context;
            z3::solver solver(context);
            const z3::expr level = context.int_const("ell");
            const z3::expr label = context.int_const("d");
            const z3::expr source = context.int_const("s");
            const z3::expr target = context.int_const("t");
            const z3::expr half_period = level + 2;
            const z3::expr total = power * label;

            solver.add(label >= 1);
            solver.add(level - label >= 14);
            solver.add(z3::mod(level - label, 2) == 0);
            solver.add(source >= 0 && source <= level);
            solver.add(target >= 0 && target <= level);
            solver.add(
                z3::mod(source + total - target, 2) == 0
            );

            const z3::expr source_above =
                source >= target;
            const z3::expr lower_label = z3::ite(
                source_above,
                source - target,
                target - source
            );
            const z3::expr ordinary_upper =
                source + target;
            const z3::expr affine_upper =
                2 * level - source - target;
            const z3::expr ordinary_upper_active =
                ordinary_upper <= affine_upper;
            const z3::expr upper_label = z3::ite(
                ordinary_upper_active,
                ordinary_upper,
                affine_upper
            );

            std::vector<z3::expr> hinges;
            hinges.push_back(z3::mod(level, 2) == 0);
            hinges.push_back(z3::mod(label, 2) == 0);
            hinges.push_back(z3::mod(source, 2) == 0);
            hinges.push_back(z3::mod(target, 2) == 0);
            hinges.push_back(source_above);
            hinges.push_back(ordinary_upper_active);

            for (int branch = 0; branch < 4; ++branch) {
                const z3::expr branch_lower_twice =
                    total - (branch + 1) * half_period + 2;
                const z3::expr branch_upper_twice =
                    total - branch * half_period;
                z3::expr fusion_lower_twice = context.int_val(0);
                z3::expr fusion_upper_twice = context.int_val(0);
                if (branch % 2 == 0) {
                    fusion_lower_twice =
                        total - branch * half_period - upper_label;
                    fusion_upper_twice =
                        total - branch * half_period - lower_label;
                } else {
                    fusion_lower_twice =
                        lower_label + total
                        - (branch + 1) * half_period + 2;
                    fusion_upper_twice =
                        upper_label + total
                        - (branch + 1) * half_period + 2;
                }

                const z3::expr lower_zero = context.int_val(0);
                const z3::expr lower_branch =
                    ceiling_half(branch_lower_twice);
                const z3::expr lower_fusion =
                    ceiling_half(fusion_lower_twice);
                const z3::expr upper_weight = total / 2;
                const z3::expr upper_branch =
                    branch_upper_twice / 2;
                const z3::expr upper_fusion =
                    fusion_upper_twice / 2;
                const z3::expr lower = maximum3(
                    lower_zero,
                    lower_branch,
                    lower_fusion
                );
                const z3::expr upper = minimum3(
                    upper_weight,
                    upper_branch,
                    upper_fusion
                );

                hinges.push_back(lower == lower_zero);
                hinges.push_back(lower == lower_branch);
                hinges.push_back(lower == lower_fusion);
                hinges.push_back(upper == upper_weight);
                hinges.push_back(upper == upper_branch);
                hinges.push_back(upper == upper_fusion);
                hinges.push_back(lower <= upper);
                hinges.push_back(lower - 1 >= 0);
                if (power >= 3) {
                    hinges.push_back(upper >= label + 1);
                    hinges.push_back(lower - 1 >= label + 1);
                }
            }
            if (hinges.size() > 63U) {
                throw std::runtime_error(
                    "endpoint mask exceeds uint64"
                );
            }

            std::set<std::uint64_t> masks;
            while (solver.check() == z3::sat) {
                const z3::model model = solver.get_model();
                std::uint64_t mask = 0U;
                z3::expr block = context.bool_val(false);
                for (std::size_t index = 0U;
                     index < hinges.size();
                     ++index) {
                    const bool active = model_boolean(
                        model,
                        hinges[index],
                        context
                    );
                    if (active) {
                        mask |= std::uint64_t{1} << index;
                        block = block || !hinges[index];
                    } else {
                        block = block || hinges[index];
                    }
                }
                if (!masks.insert(mask).second) {
                    throw std::runtime_error(
                        "blocking clause repeated endpoint mask"
                    );
                }
                solver.add(block);
            }

            total_masks += masks.size();
            std::cout
                << "SU2_FUSION_ENDPOINT_MASKS"
                << " power=" << power
                << " hinges=" << hinges.size()
                << " masks=" << masks.size()
                << " result=PASS_EXACT_CENSUS\n";
        }
        std::uint64_t total_joint_masks = 0U;
        const std::vector<std::vector<int>> packets{
            {1, 3},
            {2, 4}
        };
        for (const auto& powers : packets) {
            z3::context context;
            z3::solver solver(context);
            const z3::expr level = context.int_const("joint_ell");
            const z3::expr label = context.int_const("joint_d");
            const z3::expr source = context.int_const("joint_s");
            const z3::expr target = context.int_const("joint_t");
            const z3::expr half_period = level + 2;

            solver.add(label >= 1);
            solver.add(level - label >= 14);
            solver.add(z3::mod(level - label, 2) == 0);
            solver.add(source >= 0 && source <= level);
            solver.add(target >= 0 && target <= level);
            solver.add(
                z3::mod(
                    source + powers.front() * label - target,
                    2
                ) == 0
            );

            const z3::expr source_above = source >= target;
            const z3::expr lower_label = z3::ite(
                source_above,
                source - target,
                target - source
            );
            const z3::expr ordinary_upper = source + target;
            const z3::expr affine_upper =
                2 * level - source - target;
            const z3::expr ordinary_upper_active =
                ordinary_upper <= affine_upper;
            const z3::expr upper_label = z3::ite(
                ordinary_upper_active,
                ordinary_upper,
                affine_upper
            );

            std::vector<z3::expr> hinges;
            hinges.push_back(z3::mod(level, 2) == 0);
            hinges.push_back(z3::mod(label, 2) == 0);
            hinges.push_back(z3::mod(source, 2) == 0);
            hinges.push_back(z3::mod(target, 2) == 0);
            hinges.push_back(source_above);
            hinges.push_back(ordinary_upper_active);

            for (const int power : powers) {
                const z3::expr total = power * label;
                for (int branch = 0; branch < 4; ++branch) {
                    const z3::expr branch_lower_twice =
                        total - (branch + 1) * half_period + 2;
                    const z3::expr branch_upper_twice =
                        total - branch * half_period;
                    z3::expr fusion_lower_twice =
                        context.int_val(0);
                    z3::expr fusion_upper_twice =
                        context.int_val(0);
                    if (branch % 2 == 0) {
                        fusion_lower_twice =
                            total - branch * half_period
                            - upper_label;
                        fusion_upper_twice =
                            total - branch * half_period
                            - lower_label;
                    } else {
                        fusion_lower_twice =
                            lower_label + total
                            - (branch + 1) * half_period + 2;
                        fusion_upper_twice =
                            upper_label + total
                            - (branch + 1) * half_period + 2;
                    }

                    const z3::expr lower_zero =
                        context.int_val(0);
                    const z3::expr lower_branch =
                        ceiling_half(branch_lower_twice);
                    const z3::expr lower_fusion =
                        ceiling_half(fusion_lower_twice);
                    const z3::expr upper_weight = total / 2;
                    const z3::expr upper_branch =
                        branch_upper_twice / 2;
                    const z3::expr upper_fusion =
                        fusion_upper_twice / 2;
                    const z3::expr lower = maximum3(
                        lower_zero,
                        lower_branch,
                        lower_fusion
                    );
                    const z3::expr upper = minimum3(
                        upper_weight,
                        upper_branch,
                        upper_fusion
                    );

                    hinges.push_back(lower == lower_zero);
                    hinges.push_back(lower == lower_branch);
                    hinges.push_back(lower == lower_fusion);
                    hinges.push_back(upper == upper_weight);
                    hinges.push_back(upper == upper_branch);
                    hinges.push_back(upper == upper_fusion);
                    hinges.push_back(lower <= upper);
                    hinges.push_back(lower - 1 >= 0);
                    if (power >= 3) {
                        hinges.push_back(upper >= label + 1);
                        hinges.push_back(lower - 1 >= label + 1);
                    }
                }
            }

            std::set<std::string> masks;
            while (solver.check() == z3::sat) {
                const z3::model model = solver.get_model();
                std::string mask;
                mask.reserve(hinges.size());
                z3::expr block = context.bool_val(false);
                for (const auto& hinge : hinges) {
                    const bool active = model_boolean(
                        model,
                        hinge,
                        context
                    );
                    mask.push_back(active ? '1' : '0');
                    block = block || (active ? !hinge : hinge);
                }
                if (!masks.insert(mask).second) {
                    throw std::runtime_error(
                        "blocking clause repeated joint endpoint mask"
                    );
                }
                solver.add(block);
            }
            total_joint_masks += masks.size();
            std::cout
                << "SU2_FUSION_ENDPOINT_JOINT_MASKS"
                << " packet="
                << (powers.front() % 2 == 0 ? "even" : "odd")
                << " powers=" << powers.front() << ',' << powers.back()
                << " hinges=" << hinges.size()
                << " masks=" << masks.size()
                << " result=PASS_EXACT_CENSUS\n";
        }
        std::cout
            << "SU2_FUSION_ENDPOINT_MASKS_TOTAL"
            << " powers=4"
            << " masks=" << total_masks
            << " result=PASS_EXACT_CENSUS\n";
        std::cout
            << "SU2_FUSION_ENDPOINT_JOINT_MASKS_TOTAL"
            << " packets=2"
            << " masks=" << total_joint_masks
            << " result=PASS_EXACT_CENSUS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
