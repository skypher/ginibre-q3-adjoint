#include <limits>

#define main original_su2_k4_intermediate_main
#include "prove_su2_k4_intermediate.cpp"
#undef main

namespace {

std::size_t parse_size_argument(const char* text, const char* name) {
    std::size_t parsed = 0U;
    const std::string input(text);
    const unsigned long long value = std::stoull(input, &parsed);
    if (
        parsed != input.size()
        || value > std::numeric_limits<std::size_t>::max()
    ) {
        throw std::runtime_error(std::string(name) + " is invalid");
    }
    return static_cast<std::size_t>(value);
}

bool rooted_tree_certificate(
    const Chamber& chamber,
    const std::vector<Polynomial>& constraints,
    std::size_t selected_pair,
    int depth
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
                    depth - 1
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
                    depth - 1
                )
            ) {
                return false;
            }

            std::cout
                << "SU2_K4_INTERMEDIATE_PAIR_TREE_ROOT"
                << " mask=" << chamber.mask
                << " pair=(" << first << ',' << second << ')'
                << " ordinal=" << selected_pair
                << " depth=" << depth
                << " result=PASS_EXACT_INTEGER_TREE\n";
            return true;
        }
    }
    throw std::runtime_error("root-pair ordinal is out of range");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4) {
            throw std::runtime_error(
                "usage: MASK ROOT_PAIR_ORDINAL DEPTH"
            );
        }
        const std::size_t mask_value =
            parse_size_argument(argv[1], "mask");
        const std::size_t selected_pair =
            parse_size_argument(argv[2], "root-pair ordinal");
        const std::size_t depth_value =
            parse_size_argument(argv[3], "depth");
        if (
            mask_value > std::numeric_limits<std::uint64_t>::max()
            || depth_value == 0U
            || depth_value > static_cast<std::size_t>(
                std::numeric_limits<int>::max()
            )
        ) {
            throw std::runtime_error("mask or depth is out of range");
        }

        const std::uint64_t mask =
            static_cast<std::uint64_t>(mask_value);
        const Formula formula = make_formula();
        const Chamber chamber = make_chamber(formula, mask);
        if (!integer_feasible(chamber.constraints)) {
            throw std::runtime_error("selected mask is infeasible");
        }
        const std::vector<Polynomial> constraints =
            irredundant_constraints(chamber);
        const int depth = static_cast<int>(depth_value);

        if (
            !rooted_tree_certificate(
                chamber,
                constraints,
                selected_pair,
                depth
            )
        ) {
            std::cout
                << "SU2_K4_INTERMEDIATE_TREE_DEPTH"
                << " mask=" << mask
                << " ordinal=" << selected_pair
                << " depth=" << depth
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
            << "SU2_K4_INTERMEDIATE_TREE_DEPTH FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
