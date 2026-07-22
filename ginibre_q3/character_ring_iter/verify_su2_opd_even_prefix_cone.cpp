#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Wedge = std::pair<int, int>;
using State = std::map<Wedge, cpp_int>;

void add_wedge(State& state, int first, int second, const cpp_int& value) {
    if (first == second || value == 0) {
        return;
    }
    if (first > second) {
        state[{first, second}] += value;
    } else {
        state[{second, first}] -= value;
    }
}

State update(const State& state, int label) {
    State result;
    for (const auto& [wedge, coefficient] : state) {
        const auto [first, second] = wedge;
        for (int output = std::abs(first - label);
             output <= first + label;
             output += 2) {
            add_wedge(result, output, second, coefficient);
        }
        for (int output = std::abs(second - label);
             output <= second + label;
             output += 2) {
            add_wedge(result, first, output, coefficient);
        }
    }
    return result;
}

cpp_int coefficient(const State& state, int first, int second) {
    if (first == second) {
        return 0;
    }
    const bool reversed = first < second;
    const Wedge key = reversed ? Wedge{second, first} : Wedge{first, second};
    const auto found = state.find(key);
    if (found == state.end()) {
        return 0;
    }
    return reversed ? -found->second : found->second;
}

void print_word(int p, const std::vector<int>& suffix, int odd) {
    std::cout << "p=" << p << " even_prefix=[";
    for (std::size_t index = 0; index < suffix.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << suffix[index];
    }
    std::cout << "] odd=" << odd;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4) {
            throw std::runtime_error(
                "usage: verify_su2_opd_even_prefix_cone MAXIMUM_EVEN_LABEL "
                "MAXIMUM_EVEN_FACTORS MAXIMUM_ODD_LABEL"
            );
        }
        const int maximum_even = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        const int maximum_odd = std::stoi(argv[3]);
        if (maximum_even < 2 || maximum_factors < 0 || maximum_odd < 3) {
            throw std::runtime_error("invalid bound");
        }
        long long words = 0;
        long long cone_entries = 0;
        long long odd_extensions = 0;
        long long boundary_targets = 0;
        for (int p = 2; p <= maximum_even; p += 2) {
            std::vector<int> suffix;
            State initial;
            add_wedge(initial, p - 1, 0, 1);
            add_wedge(initial, p, 1, -1);
            const auto visit = [&](const auto& self, int next_even,
                                   const State& state) -> void {
                ++words;
                for (const auto& [wedge, value] : state) {
                    const auto [larger, smaller] = wedge;
                    ++cone_entries;
                    if (((larger - smaller) & 1) == 0
                        || (((larger & 1) == 0) ? value > 0 : value < 0)) {
                        std::cout << "FAIL signed-cone wedge=(" << larger
                                  << ',' << smaller << ") value=" << value
                                  << ' ';
                        print_word(p, suffix, -1);
                        std::cout << '\n';
                        std::exit(EXIT_FAILURE);
                    }
                }
                for (int odd = p | 1; odd <= maximum_odd; odd += 2) {
                    ++odd_extensions;
                    const State extended = update(state, odd);
                    int support = p + 1 + odd;
                    for (const int label : suffix) {
                        support += label;
                    }
                    for (int target = 1; target <= support; ++target) {
                        ++boundary_targets;
                        const cpp_int boundary
                            = coefficient(extended, target, 0);
                        if (boundary < 0) {
                            std::cout << "FAIL odd-boundary target=" << target
                                      << " value=" << boundary << ' ';
                            print_word(p, suffix, odd);
                            std::cout << '\n';
                            std::exit(EXIT_FAILURE);
                        }
                        cpp_int reservoir = 0;
                        for (int output = std::abs(target - odd);
                             output <= target + odd;
                             output += 2) {
                            reservoir += coefficient(state, output, 0);
                        }
                        const cpp_int carrier
                            = coefficient(state, target, odd);
                        if (boundary != reservoir + carrier) {
                            std::cout << "FAIL transfer-identity target="
                                      << target << ' ';
                            print_word(p, suffix, odd);
                            std::cout << '\n';
                            std::exit(EXIT_FAILURE);
                        }
                    }
                }
                if (suffix.size()
                    == static_cast<std::size_t>(maximum_factors)) {
                    return;
                }
                for (int label = next_even;
                     label <= maximum_even;
                     label += 2) {
                    suffix.push_back(label);
                    const State next = update(state, label);
                    self(self, label, next);
                    suffix.pop_back();
                }
            };
            visit(visit, p, initial);
        }
        std::cout << "SU2_OPD_EVEN_PREFIX_CONE PASS words=" << words
                  << " cone_entries=" << cone_entries
                  << " odd_extensions=" << odd_extensions
                  << " boundary_targets=" << boundary_targets
                  << " maximum_even_label=" << maximum_even
                  << " maximum_even_factors=" << maximum_factors
                  << " maximum_odd_label=" << maximum_odd << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
