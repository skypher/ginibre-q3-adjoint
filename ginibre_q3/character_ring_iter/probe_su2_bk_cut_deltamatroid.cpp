#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int parse_nonnegative(const char* text, const char* name) {
    const std::string input(text);
    std::size_t used = 0U;
    const long long value = std::stoll(input, &used, 10);
    if (used != input.size() || value < 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
}

bool cut_feasible(
    const std::vector<int>& alpha,
    const std::vector<int>& beta,
    std::uint64_t cut
) {
    int heights[2]{0, 0};
    for (std::size_t position = 0U; position < alpha.size(); ++position) {
        const std::size_t colour = (cut & (std::uint64_t{1} << position))
            == 0U ? 0U : 1U;
        if (heights[colour] < beta[position]) {
            return false;
        }
        heights[colour] += alpha[position] - beta[position];
    }
    return heights[0] == 0 && heights[1] == 0;
}

std::vector<std::uint64_t> feasible_cuts(
    const std::vector<int>& alpha,
    const std::vector<int>& beta
) {
    if (alpha.size() >= std::numeric_limits<std::uint64_t>::digits) {
        throw std::runtime_error("too many vertices for cut masks");
    }
    const std::uint64_t limit = std::uint64_t{1} << alpha.size();
    std::vector<std::uint64_t> output;
    for (std::uint64_t cut = 0U; cut < limit; ++cut) {
        if (cut_feasible(alpha, beta, cut)) {
            output.push_back(cut);
        }
    }
    return output;
}

void print_vector(const std::vector<int>& values) {
    std::cout << '[';
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << ']';
}

[[noreturn]] void report_failure(
    const std::vector<int>& alpha,
    const std::vector<int>& beta,
    const std::vector<std::uint64_t>& cuts,
    std::uint64_t first,
    std::uint64_t second,
    std::size_t element
) {
    std::cout << "SU2_BK_CUT_DELTAMATROID result=FAIL alpha=";
    print_vector(alpha);
    std::cout << " beta=";
    print_vector(beta);
    std::cout << " cuts={";
    for (std::size_t index = 0U; index < cuts.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << cuts[index];
    }
    std::cout << "} first=" << first
              << " second=" << second
              << " element=" << element << '\n';
    std::exit(EXIT_FAILURE);
}

void test_symmetric_exchange(
    const std::vector<int>& alpha,
    const std::vector<int>& beta,
    std::uint64_t& states,
    std::uint64_t& fibres
) {
    ++states;
    const std::vector<std::uint64_t> cuts = feasible_cuts(alpha, beta);
    ++fibres;
    const auto contains = [&cuts](std::uint64_t cut) {
        return std::binary_search(cuts.begin(), cuts.end(), cut);
    };
    for (const std::uint64_t first : cuts) {
        for (const std::uint64_t second : cuts) {
            const std::uint64_t difference = first ^ second;
            for (std::size_t element = 0U; element < alpha.size(); ++element) {
                const std::uint64_t element_bit = std::uint64_t{1} << element;
                if ((difference & element_bit) == 0U) {
                    continue;
                }
                bool exchanged = false;
                for (std::size_t partner = 0U;
                     partner < alpha.size(); ++partner) {
                    const std::uint64_t partner_bit
                        = std::uint64_t{1} << partner;
                    if ((difference & partner_bit) == 0U) {
                        continue;
                    }
                    const std::uint64_t candidate = partner == element
                        ? first ^ element_bit
                        : first ^ element_bit ^ partner_bit;
                    if (contains(candidate)) {
                        exchanged = true;
                        break;
                    }
                }
                if (!exchanged) {
                    report_failure(
                        alpha, beta, cuts, first, second, element
                    );
                }
            }
        }
    }
}

void enumerate_paths(
    std::size_t position,
    int height,
    int maximum_part,
    std::vector<int>& alpha,
    std::vector<int>& beta,
    std::uint64_t& states,
    std::uint64_t& fibres
) {
    if (position == alpha.size()) {
        if (height == 0) {
            test_symmetric_exchange(alpha, beta, states, fibres);
        }
        return;
    }
    const std::size_t remaining = alpha.size() - position - 1U;
    const int maximum_future_drop
        = static_cast<int>(remaining) * maximum_part;
    for (int top = 0; top <= maximum_part; ++top) {
        for (int bottom = 0; bottom <= maximum_part - top; ++bottom) {
            if (top + bottom == 0) {
                continue;
            }
            const int next_height = height + top - bottom;
            if (next_height < 0 || next_height > maximum_future_drop) {
                continue;
            }
            alpha[position] = top;
            beta[position] = bottom;
            enumerate_paths(
                position + 1U, next_height, maximum_part,
                alpha, beta, states, fibres
            );
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: probe_su2_bk_cut_deltamatroid vertices maximum_part"
            );
        }
        const int vertices_input = parse_nonnegative(argv[1], "vertex count");
        const int maximum_part = parse_nonnegative(argv[2], "maximum part");
        if (vertices_input < 1 || maximum_part < 1 || vertices_input > 20) {
            throw std::runtime_error("search bounds are unsupported");
        }
        const std::size_t vertices
            = static_cast<std::size_t>(vertices_input);
        std::vector<int> alpha(vertices, 0);
        std::vector<int> beta(vertices, 0);
        std::uint64_t states = 0U;
        std::uint64_t fibres = 0U;
        enumerate_paths(
            0U, 0, maximum_part, alpha, beta, states, fibres
        );
        std::cout << "SU2_BK_CUT_DELTAMATROID vertices=" << vertices
                  << " maximum_part=" << maximum_part
                  << " queue_paths=" << states
                  << " fibres=" << fibres
                  << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
