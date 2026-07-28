#include <limits>

#define main original_su2_k4_intermediate_main
#include "prove_su2_k4_intermediate.cpp"
#undef main

namespace {

std::size_t parse_index(const char* text) {
    std::size_t parsed = 0U;
    const std::string input(text);
    const unsigned long long value = std::stoull(input, &parsed);
    if (
        parsed != input.size()
        || value > std::numeric_limits<std::size_t>::max()
    ) {
        throw std::runtime_error("tree-pair ordinal is invalid");
    }
    return static_cast<std::size_t>(value);
}

bool pair_cut_tree_root_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints,
    std::size_t selected_pair
) {
    std::size_t ordinal = 0U;
    for (std::size_t first = 0U;
         first < constraints.size();
         ++first) {
        for (std::size_t second = first + 1U;
             second < constraints.size();
             ++second) {
            if (ordinal++ != selected_pair) {
                continue;
            }
            std::vector<Polynomial> first_branch = constraints;
            first_branch.push_back(
                constraints[first] - constraints[second]
            );
            if (
                !pair_cut_tree_constraints(
                    chamber,
                    first_branch,
                    1
                )
            ) {
                return false;
            }
            std::vector<Polynomial> second_branch = constraints;
            second_branch.push_back(
                constraints[second]
                    - constraints[first]
                    - constant(1)
            );
            if (
                !pair_cut_tree_constraints(
                    chamber,
                    second_branch,
                    1
                )
            ) {
                return false;
            }
            std::cout
                << "SU2_K4_INTERMEDIATE_PAIR_TREE_ROOT"
                << " mask=" << chamber.mask
                << " pair=(" << first << ',' << second << ')'
                << " ordinal=" << selected_pair
                << " depth=2"
                << " result=PASS_EXACT_INTEGER_TREE\n";
            return true;
        }
    }
    throw std::runtime_error("tree-pair ordinal is out of range");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (
            argc != 3
            && !(
                argc == 4
                && std::string(argv[2]) == "tree-pair"
            )
        ) {
            throw std::runtime_error(
                "usage: MASK "
                "constant|bounded|sum|equal|double|pair|tree "
                "| MASK tree-pair ORDINAL"
            );
        }
        std::size_t parsed = 0U;
        const std::string mask_text = argv[1];
        const unsigned long long raw_mask =
            std::stoull(mask_text, &parsed);
        if (
            parsed != mask_text.size()
            || raw_mask > std::numeric_limits<std::uint64_t>::max()
        ) {
            throw std::runtime_error("invalid mask");
        }
        const std::uint64_t mask = static_cast<std::uint64_t>(raw_mask);
        const std::string method = argv[2];
        const Formula formula = make_formula();
        const Chamber chamber = make_chamber(formula, mask);
        if (!integer_feasible(chamber.constraints)) {
            throw std::runtime_error("selected mask is infeasible");
        }
        const std::vector<Polynomial> constraints =
            irredundant_constraints(chamber);

        bool passed = false;
        if (method == "constant") {
            passed = constant_three_sum_certificate(
                chamber,
                constraints
            );
        } else if (method == "bounded") {
            passed = bounded_integer_certificate(
                chamber,
                formula,
                constraints
            );
        } else if (method == "sum") {
            passed = sum_cone_certificate(chamber, constraints);
        } else if (method == "equal") {
            passed = equal_sum_square_cone_certificate(
                chamber,
                constraints
            );
        } else if (method == "double") {
            passed = double_sum_cone_certificate(chamber, constraints);
        } else if (method == "pair") {
            passed = pair_cut_certificate(chamber, constraints);
        } else if (method == "tree") {
            passed = pair_cut_tree_certificate(chamber, constraints);
        } else if (method == "tree-pair") {
            passed = pair_cut_tree_root_certificate(
                chamber,
                constraints,
                parse_index(argv[3])
            );
        } else {
            throw std::runtime_error("unknown certificate method");
        }

        if (!passed) {
            std::cout
                << "SU2_K4_INTERMEDIATE_METHOD"
                << " mask=" << mask
                << " method=" << method
                << " result=NO_CERTIFICATE\n";
            return EXIT_FAILURE;
        }
        std::cout
            << "SU2_K4_INTERMEDIATE"
            << " hinges=" << formula.hinges.size()
            << " feasible_chambers=1"
            << " certified_chambers=1"
            << " result=PASS_EXACT_CERTIFICATE\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_K4_INTERMEDIATE_METHOD FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
