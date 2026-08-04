#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <z3++.h>

namespace {

z3::expr maximum2(const z3::expr& left, const z3::expr& right) {
    return z3::ite(left >= right, left, right);
}

z3::expr maximum3(
    const z3::expr& first,
    const z3::expr& second,
    const z3::expr& third
) {
    return z3::ite(
        first >= second,
        z3::ite(first >= third, first, third),
        z3::ite(second >= third, second, third)
    );
}

z3::expr minimum2(const z3::expr& left, const z3::expr& right) {
    return z3::ite(left <= right, left, right);
}

z3::expr minimum3(
    const z3::expr& first,
    const z3::expr& second,
    const z3::expr& third
) {
    return z3::ite(
        first <= second,
        z3::ite(first <= third, first, third),
        z3::ite(second <= third, second, third)
    );
}

z3::expr ceiling_half(const z3::expr& value) {
    return -((-value) / 2);
}

int parse_parity_class(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value < 0
        || value > 15
    ) {
        throw std::runtime_error("parity class must be an integer in 0..15");
    }
    return static_cast<int>(value);
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

void append_crossing_hinges(
    z3::context& context,
    const z3::expr& level,
    const z3::expr& label,
    const z3::expr& crossing,
    int power,
    bool active_only,
    std::vector<z3::expr>& hinges
) {
    const z3::expr product = power * label;
    const z3::expr lower_depth = level - label - crossing;
    const z3::expr upper_depth = label + crossing;
    for (int branch = 0; branch < 3; ++branch) {
        z3::expr lower_label = context.int_val(0);
        z3::expr upper_label = context.int_val(0);
        if (branch == 0) {
            lower_label = lower_depth;
            upper_label = upper_depth;
        } else if (branch == 1) {
            lower_label = 2 * level - label - crossing + 1;
            upper_label = level + label + crossing + 1;
        } else {
            lower_label = 2 * level + 2 + lower_depth;
            upper_label = 2 * level + 2 + upper_depth;
        }
        const z3::expr raw_lower = product - upper_label;
        const z3::expr raw_upper = product - lower_label;
        const z3::expr lower = maximum2(context.int_val(0), raw_lower);
        const z3::expr upper = minimum2(product, raw_upper);
        const z3::expr nonempty = lower <= upper;
        hinges.push_back(nonempty);
        // Empty endpoint runs contribute zero, so active-packet modes retain
        // only the selectors that can change an endpoint polynomial.
        const auto append_selector = [&hinges, &nonempty, active_only](
            const z3::expr& selector
        ) {
            hinges.push_back(
                active_only ? z3::implies(nonempty, selector) : selector
            );
        };
        append_selector(lower == context.int_val(0));
        append_selector(lower == raw_lower);
        append_selector(upper == product);
        append_selector(upper == raw_upper);
        append_selector(lower - 1 >= 0);
        for (int image = 1; image <= power; ++image) {
            const z3::expr activation = image * (2 * label + 1);
            append_selector(upper >= activation);
            append_selector(lower - 1 >= activation);
        }
    }
}

void append_terminal_hinges(
    z3::context& context,
    const z3::expr& level,
    const z3::expr& label,
    const z3::expr& crossing,
    const z3::expr& target,
    int power,
    bool active_only,
    std::vector<z3::expr>& hinges
) {
    const z3::expr ell = level - 1;
    const z3::expr terminal_label = 2 * label;
    const z3::expr source_label = 2 * crossing;
    const z3::expr target_label = 2 * target;
    const z3::expr half_period = ell + 2;
    const z3::expr total = power * terminal_label;
    // Lemma 5A8H28P1A cancels the simple-current shifts, so the terminal
    // packet is the unshifted lower-level power N_(2Q)^power at 2V,2x.
    const z3::expr source_above = source_label >= target_label;
    const z3::expr lower_label = z3::ite(
        source_above,
        source_label - target_label,
        target_label - source_label
    );
    const z3::expr ordinary_upper = source_label + target_label;
    const z3::expr affine_upper =
        2 * ell - source_label - target_label;
    const z3::expr ordinary_upper_active =
        ordinary_upper <= affine_upper;
    const z3::expr upper_label = z3::ite(
        ordinary_upper_active,
        ordinary_upper,
        affine_upper
    );
    hinges.push_back(source_above);
    hinges.push_back(ordinary_upper_active);

    for (int branch = 0; branch < 4; ++branch) {
        const z3::expr branch_lower_twice =
            total - (branch + 1) * half_period + 2;
        const z3::expr branch_upper_twice = total - branch * half_period;
        z3::expr fusion_lower_twice = context.int_val(0);
        z3::expr fusion_upper_twice = context.int_val(0);
        if (branch % 2 == 0) {
            fusion_lower_twice =
                total - branch * half_period - upper_label;
            fusion_upper_twice =
                total - branch * half_period - lower_label;
        } else {
            fusion_lower_twice = lower_label + total
                - (branch + 1) * half_period + 2;
            fusion_upper_twice = upper_label + total
                - (branch + 1) * half_period + 2;
        }
        const z3::expr lower_zero = context.int_val(0);
        const z3::expr lower_branch = ceiling_half(branch_lower_twice);
        const z3::expr lower_fusion = ceiling_half(fusion_lower_twice);
        const z3::expr upper_weight = total / 2;
        const z3::expr upper_branch = branch_upper_twice / 2;
        const z3::expr upper_fusion = fusion_upper_twice / 2;
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
        const z3::expr nonempty = lower <= upper;
        hinges.push_back(nonempty);
        const auto append_selector = [&hinges, &nonempty, active_only](
            const z3::expr& selector
        ) {
            hinges.push_back(
                active_only ? z3::implies(nonempty, selector) : selector
            );
        };
        append_selector(lower == lower_zero);
        append_selector(lower == lower_branch);
        append_selector(lower == lower_fusion);
        append_selector(upper == upper_weight);
        append_selector(upper == upper_branch);
        append_selector(upper == upper_fusion);
        append_selector(lower - 1 >= 0);
        if (power >= 3) {
            append_selector(upper >= terminal_label + 1);
            append_selector(lower - 1 >= terminal_label + 1);
        }
    }
}

std::uint64_t enumerate_masks(const std::string& mode, int parity_class) {
    z3::context context;
    z3::solver solver(context);
    const z3::expr level = context.int_const("K");
    const z3::expr label = context.int_const("Q");
    const z3::expr crossing = context.int_const("V");
    const z3::expr target = context.int_const("x");
    const z3::expr width = level - 1 - 2 * label;

    solver.add(label >= 7);
    solver.add(width >= 10);
    solver.add(crossing >= 0 && 2 * crossing < level);
    solver.add(target >= 0 && 2 * target < level);
    solver.add(level - label - crossing <= label + crossing);
    const bool diagonal_g0 = mode == "g0";
    if (diagonal_g0) {
        solver.add(target == crossing);
    }

    std::vector<z3::expr> hinges;
    if (parity_class < 0) {
        hinges.push_back(z3::mod(level, 2) == 0);
        hinges.push_back(z3::mod(label, 2) == 0);
        hinges.push_back(z3::mod(crossing, 2) == 0);
        if (!diagonal_g0) {
            hinges.push_back(z3::mod(target, 2) == 0);
        }
    } else {
        solver.add(
            z3::mod(level, 2)
            == ((parity_class >> 0) & 1)
        );
        solver.add(
            z3::mod(label, 2)
            == ((parity_class >> 1) & 1)
        );
        solver.add(
            z3::mod(crossing, 2)
            == ((parity_class >> 2) & 1)
        );
        solver.add(
            z3::mod(target, 2)
            == ((parity_class >> 3) & 1)
        );
    }
    if (!diagonal_g0) {
        hinges.push_back(crossing == target);
    }
    hinges.push_back(label - 2 * width - 3 >= 0);
    const bool active_only = mode == "even-active" || mode == "odd-active";

    const auto append_crossing_powers = [&context, &level, &label,
        &crossing, &hinges, active_only](const std::vector<int>& powers) {
        for (const int power : powers) {
            append_crossing_hinges(
                context,
                level,
                label,
                crossing,
                power,
                active_only,
                hinges
            );
        }
    };
    const auto append_terminal_powers = [&context, &level, &label,
        &crossing, &target, &hinges, active_only](
            const std::vector<int>& powers
        ) {
        for (const int power : powers) {
            append_terminal_hinges(
                context,
                level,
                label,
                crossing,
                target,
                power,
                active_only,
                hinges
            );
        }
    };
    if (mode == "crossing" || mode == "full") {
        append_crossing_powers({1, 2, 3, 4, 5});
    }
    if (mode == "terminal" || mode == "full") {
        append_terminal_powers({1, 2, 3, 4});
    }
    if (mode == "even") {
        append_crossing_powers({1, 2, 3, 4, 5});
        append_terminal_powers({2, 4});
    }
    if (mode == "odd") {
        append_crossing_powers({1, 2, 3, 4});
        append_terminal_powers({1, 3});
    }
    if (mode == "even-active") {
        append_crossing_powers({1, 2, 3, 4, 5});
        append_terminal_powers({2, 4});
    }
    if (mode == "odd-active") {
        append_crossing_powers({1, 2, 3, 4});
        append_terminal_powers({1, 3});
    }
    if (mode == "g0") {
        append_crossing_powers({4, 5});
    }
    if (hinges.size() > 1000U) {
        throw std::runtime_error("joint activation mask is unexpectedly wide");
    }

    std::uint64_t masks = 0U;
    while (solver.check() == z3::sat) {
        const z3::model model = solver.get_model();
        z3::expr block = context.bool_val(false);
        for (const z3::expr& hinge : hinges) {
            const bool active = model_boolean(model, hinge, context);
            block = block || (active ? !hinge : hinge);
        }
        solver.add(block);
        if (masks == std::numeric_limits<std::uint64_t>::max()) {
            throw std::runtime_error("activation-mask counter overflow");
        }
        ++masks;
        if (masks % 100U == 0U) {
            std::cerr
                << "SU2_SHELL_JOINT_ACTIVATION_MASKS_PROGRESS"
                << " mode=" << mode
                << " parity_class="
                << (
                    parity_class < 0
                        ? std::string("all")
                        : std::to_string(parity_class)
                )
                << " masks=" << masks << '\n';
        }
    }
    std::cout
        << "SU2_SHELL_JOINT_ACTIVATION_MASKS"
        << " mode=" << mode
        << " parity_class="
        << (
            parity_class < 0
                ? std::string("all")
                : std::to_string(parity_class)
        )
        << " hinges=" << hinges.size()
        << " masks=" << masks
        << " result=PASS_EXACT_PRESBURGER_CENSUS\n";
    return masks;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2 && argc != 4) {
            throw std::runtime_error(
                "usage: MODE [--parity 0..15]"
            );
        }
        const std::string mode(argv[1]);
        int parity_class = -1;
        if (argc == 4) {
            if (std::string(argv[2]) != "--parity") {
                throw std::runtime_error(
                    "the only optional selector is --parity 0..15"
                );
            }
            parity_class = parse_parity_class(argv[3]);
        }
        if (
            mode != "crossing"
            && mode != "terminal"
            && mode != "even"
            && mode != "odd"
            && mode != "even-active"
            && mode != "odd-active"
            && mode != "g0"
            && mode != "full"
        ) {
            throw std::runtime_error(
                "mode must be crossing, terminal, even, odd, even-active, "
                "odd-active, g0, or full"
            );
        }
        static_cast<void>(enumerate_masks(mode, parity_class));
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
