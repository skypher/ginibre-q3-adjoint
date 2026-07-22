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

template <class Function>
void for_each_output(int level, int left, int right, Function function) {
    const int maximum = std::min(left + right, 2 * level - left - right);
    for (int output = std::abs(left - right);
         output <= maximum;
         output += 2) {
        function(output);
    }
}

State update(int level, const State& state, int label) {
    State result;
    for (const auto& [wedge, coefficient] : state) {
        const auto [first, second] = wedge;
        for_each_output(level, first, label, [&](const int output) {
            add_wedge(result, output, second, coefficient);
        });
        for_each_output(level, second, label, [&](const int output) {
            add_wedge(result, first, output, coefficient);
        });
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

void print_word(int level, int p, const std::vector<int>& suffix, int odd) {
    std::cout << "level=" << level << " p=" << p << " even_prefix=[";
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
        if (argc != 3) {
            throw std::runtime_error(
                "usage: verify_su2_level_even_prefix_cone MAXIMUM_LEVEL "
                "MAXIMUM_EVEN_FACTORS"
            );
        }
        const int maximum_level = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_level < 2 || maximum_factors < 0) {
            throw std::runtime_error("invalid bound");
        }
        long long words = 0;
        long long cone_entries = 0;
        long long odd_extensions = 0;
        long long boundary_targets = 0;
        bool lower_packet_nonnegative = true;
        for (int level = 2; level <= maximum_level; ++level) {
            for (int p = 2; p <= level; p += 2) {
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
                            || (((larger & 1) == 0)
                                    ? value > 0 : value < 0)) {
                            std::cout << "FAIL signed-cone wedge=(" << larger
                                      << ',' << smaller << ") value=" << value
                                      << ' ';
                            print_word(level, p, suffix, -1);
                            std::cout << '\n';
                            std::exit(EXIT_FAILURE);
                        }
                    }
                    for (int odd = p | 1; odd <= level; odd += 2) {
                        ++odd_extensions;
                        const State extended = update(level, state, odd);
                        State outer_extended;
                        if (p < level) {
                            add_wedge(outer_extended, p + 1, 0, 1);
                        }
                        for (const int label : suffix) {
                            outer_extended = update(
                                level, outer_extended, label
                            );
                        }
                        outer_extended = update(level, outer_extended, odd);
                        for (int target = 1; target <= level; ++target) {
                            ++boundary_targets;
                            const cpp_int boundary
                                = coefficient(extended, target, 0);
                            const cpp_int outer_boundary
                                = coefficient(outer_extended, target, 0);
                            if (boundary < 0) {
                                if (lower_packet_nonnegative) {
                                    lower_packet_nonnegative = false;
                                    std::cout << "LOWER_PACKET_FAIL target="
                                              << target << " value=" << boundary
                                              << " outer=" << outer_boundary
                                              << " total="
                                              << boundary + outer_boundary
                                              << ' ';
                                    print_word(level, p, suffix, odd);
                                    std::cout << '\n';
                                }
                            }
                            if (boundary + outer_boundary < 0) {
                                std::cout << "FAIL combined-boundary target="
                                          << target << " value="
                                          << boundary + outer_boundary << ' ';
                                print_word(level, p, suffix, odd);
                                std::cout << '\n';
                                std::exit(EXIT_FAILURE);
                            }
                            cpp_int reservoir = 0;
                            for_each_output(
                                level, target, odd, [&](const int output) {
                                    reservoir += coefficient(
                                        state, output, 0
                                    );
                                }
                            );
                            const cpp_int carrier
                                = coefficient(state, target, odd);
                            if (boundary != reservoir + carrier) {
                                std::cout << "FAIL transfer-identity target="
                                          << target << ' ';
                                print_word(level, p, suffix, odd);
                                std::cout << '\n';
                                std::exit(EXIT_FAILURE);
                            }
                        }
                    }
                    if (suffix.size()
                        == static_cast<std::size_t>(maximum_factors)) {
                        return;
                    }
                    for (int label = next_even; label <= level; label += 2) {
                        suffix.push_back(label);
                        const State next = update(level, state, label);
                        self(self, label, next);
                        suffix.pop_back();
                    }
                };
                visit(visit, p, initial);
            }
        }
        std::cout << "SU2_LEVEL_EVEN_PREFIX_CONE PASS words=" << words
                  << " cone_entries=" << cone_entries
                  << " odd_extensions=" << odd_extensions
                  << " boundary_targets=" << boundary_targets
                  << " lower_packet_nonnegative="
                  << (lower_packet_nonnegative ? "true" : "false")
                  << " maximum_level=" << maximum_level
                  << " maximum_even_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
