#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
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

cpp_int boundary(const State& state, int target) {
    const auto found = state.find({target, 0});
    return found == state.end() ? cpp_int{0} : found->second;
}

void print_word(const std::vector<int>& word) {
    std::cout << '[';
    for (std::size_t index = 0; index < word.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << word[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4 && argc != 5) {
            throw std::runtime_error(
                "usage: search_su2_even_wedge_boundary MAXIMUM_EVEN_LABEL "
                "MAXIMUM_FACTORS MAXIMUM_INITIAL_INDEX_OR_ODD_LABEL "
                "[--ordered-atoms]"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        const int maximum_initial = std::stoi(argv[3]);
        const bool ordered_atoms = argc == 5
            && std::string(argv[4]) == "--ordered-atoms";
        if (argc == 5 && !ordered_atoms) {
            throw std::runtime_error("unknown mode");
        }
        if (maximum_label < 2 || maximum_factors < 0 || maximum_initial < 2) {
            throw std::runtime_error("invalid bound");
        }
        if (ordered_atoms) {
            long long words = 0;
            long long atoms = 0;
            long long targets = 0;
            for (int p = 2; p <= maximum_label; p += 2) {
                std::vector<int> word;
                const auto visit = [&](const auto& self, int next_label) -> void {
                    ++words;
                    int total = 0;
                    for (const int label : word) {
                        total += label;
                    }
                    for (int odd = p | 1; odd <= maximum_initial; odd += 2) {
                        std::vector<Wedge> initial_atoms;
                        for (int output = std::abs(odd - (p - 1));
                             output <= odd + p - 1;
                             output += 2) {
                            if (output > 0) {
                                initial_atoms.emplace_back(output, 0);
                            }
                        }
                        for (const int output : {odd - 1, odd + 1}) {
                            if (output > p) {
                                initial_atoms.emplace_back(output, p);
                            }
                        }
                        for (const auto& atom : initial_atoms) {
                            ++atoms;
                            State state;
                            add_wedge(state, atom.first, atom.second, 1);
                            for (const int label : word) {
                                state = update(state, label);
                            }
                            const int support = atom.first + atom.second + total;
                            for (int target = 1; target <= support; ++target) {
                                ++targets;
                                const cpp_int value = boundary(state, target);
                                if (value < 0) {
                                    std::cout << "FAIL ordered-atom p=" << p
                                              << " odd=" << odd << " atom=("
                                              << atom.first << ',' << atom.second
                                              << ") even_word=";
                                    print_word(word);
                                    std::cout << " target=" << target
                                              << " value=" << value << '\n';
                                    std::exit(EXIT_FAILURE);
                                }
                            }
                        }
                    }
                    if (word.size()
                        == static_cast<std::size_t>(maximum_factors)) {
                        return;
                    }
                    for (int label = next_label;
                         label <= maximum_label;
                         label += 2) {
                        word.push_back(label);
                        self(self, label);
                        word.pop_back();
                    }
                };
                visit(visit, p);
            }
            std::cout << "SU2_EVEN_WEDGE_ORDERED_ATOMS PASS words=" << words
                      << " atoms=" << atoms << " boundary_targets=" << targets
                      << " maximum_even_label=" << maximum_label
                      << " maximum_factors=" << maximum_factors
                      << " maximum_odd_label=" << maximum_initial << '\n';
            return EXIT_SUCCESS;
        }
        long long words = 0;
        long long initial_wedges = 0;
        long long boundary_targets = 0;
        std::vector<int> word;
        const auto visit = [&](const auto& self, int next_label) -> void {
            ++words;
            int total = 0;
            for (const int label : word) {
                total += label;
            }
            for (int larger = 2; larger <= maximum_initial; larger += 2) {
                for (int smaller = 0; smaller < larger; smaller += 2) {
                    ++initial_wedges;
                    State state;
                    add_wedge(state, larger, smaller, 1);
                    for (const int label : word) {
                        state = update(state, label);
                    }
                    const int support = larger + smaller + total;
                    for (int target = 1; target <= support; ++target) {
                        ++boundary_targets;
                        const cpp_int value = boundary(state, target);
                        if (value < 0) {
                            std::cout << "FAIL initial=(" << larger << ','
                                      << smaller << ") even_word=";
                            print_word(word);
                            std::cout << " target=" << target
                                      << " value=" << value << '\n';
                            std::exit(EXIT_FAILURE);
                        }
                    }
                }
            }
            if (word.size() == static_cast<std::size_t>(maximum_factors)) {
                return;
            }
            for (int label = next_label; label <= maximum_label; label += 2) {
                word.push_back(label);
                self(self, label);
                word.pop_back();
            }
        };
        visit(visit, 2);
        std::cout << "SU2_EVEN_WEDGE_BOUNDARY PASS words=" << words
                  << " initial_wedges=" << initial_wedges
                  << " boundary_targets=" << boundary_targets
                  << " maximum_even_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors
                  << " maximum_initial_index=" << maximum_initial << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
