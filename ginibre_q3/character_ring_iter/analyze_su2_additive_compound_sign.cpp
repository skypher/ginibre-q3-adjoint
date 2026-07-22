#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace {

struct DisjointSet {
    explicit DisjointSet(std::size_t size)
        : parent(size), rank(size, 0), parity(size, 0) {
        for (std::size_t index = 0; index < size; ++index) {
            parent[index] = index;
        }
    }

    std::pair<std::size_t, int> find(std::size_t vertex) {
        if (parent[vertex] == vertex) {
            return {vertex, 0};
        }
        const auto [root, parent_parity] = find(parent[vertex]);
        parity[vertex] ^= parent_parity;
        parent[vertex] = root;
        return {root, parity[vertex]};
    }

    bool impose(std::size_t first, std::size_t second, int difference) {
        auto [first_root, first_parity] = find(first);
        auto [second_root, second_parity] = find(second);
        if (first_root == second_root) {
            return (first_parity ^ second_parity) == difference;
        }
        if (rank[first_root] < rank[second_root]) {
            std::swap(first_root, second_root);
            std::swap(first_parity, second_parity);
        }
        parent[second_root] = first_root;
        parity[second_root]
            = first_parity ^ second_parity ^ difference;
        if (rank[first_root] == rank[second_root]) {
            ++rank[first_root];
        }
        return true;
    }

    std::vector<std::size_t> parent;
    std::vector<int> rank;
    std::vector<int> parity;
};

using State = std::pair<int, int>;

bool fuses(int input, int label, int output, int level) {
    if (output < std::abs(input - label)
        || output > input + label
        || ((input + label + output) & 1) != 0) {
        return false;
    }
    return level < 0 || output <= 2 * level - input - label;
}

void add_wedge(
    std::map<State, int>& column,
    int first,
    int second,
    int coefficient
) {
    if (first == second) {
        return;
    }
    if (first < second) {
        std::swap(first, second);
        coefficient = -coefficient;
    }
    column[{first, second}] += coefficient;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 3 || argc > 5) {
            throw std::runtime_error(
                "usage: analyze_su2_additive_compound_sign "
                "MAXIMUM_STATE MAXIMUM_LABEL [LEVEL [MINIMUM_LABEL]]"
            );
        }
        const int maximum_state = std::stoi(argv[1]);
        const int maximum_label = std::stoi(argv[2]);
        const int level = argc >= 4 ? std::stoi(argv[3]) : -1;
        const int minimum_label = argc == 5 ? std::stoi(argv[4]) : 2;
        if (maximum_state < 1 || maximum_label < 2
            || minimum_label < 0 || minimum_label > maximum_label
            || (level >= 0 && maximum_state != level)) {
            throw std::runtime_error("invalid bounds");
        }

        std::vector<State> states;
        std::map<State, std::size_t> index;
        for (int first = 1; first <= maximum_state; ++first) {
            for (int second = 0; second < first; ++second) {
                index[{first, second}] = states.size();
                states.emplace_back(first, second);
            }
        }
        DisjointSet signs(states.size());
        std::size_t constraints = 0U;
        for (int label = minimum_label;
             label <= maximum_label; ++label) {
            for (std::size_t input_index = 0U;
                 input_index < states.size(); ++input_index) {
                const auto [first, second] = states[input_index];
                std::map<State, int> column;
                for (int output = 0; output <= maximum_state; ++output) {
                    if (fuses(first, label, output, level)) {
                        add_wedge(column, output, second, 1);
                    }
                    if (fuses(second, label, output, level)) {
                        add_wedge(column, first, output, 1);
                    }
                }
                for (const auto& [output, coefficient] : column) {
                    if (coefficient == 0) {
                        continue;
                    }
                    const std::size_t output_index = index.at(output);
                    const int difference = coefficient < 0 ? 1 : 0;
                    ++constraints;
                    if (!signs.impose(
                            input_index, output_index, difference
                        )) {
                        std::cout
                            << "SU2_ADDITIVE_COMPOUND_COMMON_SIGN result=FAIL"
                            << " label=" << label
                            << " input=(" << first << ',' << second << ')'
                            << " output=(" << output.first << ','
                            << output.second << ')'
                            << " coefficient=" << coefficient
                            << " maximum_state=" << maximum_state
                            << " maximum_label=" << maximum_label
                            << " minimum_label=" << minimum_label
                            << " level=" << level << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
        }

        const std::size_t initial = index.at({1, 0});
        for (int target = 1; target <= maximum_state; ++target) {
            const std::size_t boundary = index.at({target, 0});
            auto [initial_root, initial_parity] = signs.find(initial);
            auto [boundary_root, boundary_parity] = signs.find(boundary);
            if (initial_root == boundary_root
                && initial_parity != boundary_parity) {
                std::cout
                    << "SU2_ADDITIVE_COMPOUND_COMMON_SIGN result=FAIL"
                    << " boundary_target=" << target
                    << " forced_boundary_sign=negative"
                    << " maximum_state=" << maximum_state
                    << " maximum_label=" << maximum_label
                    << " minimum_label=" << minimum_label
                    << " level=" << level << '\n';
                return EXIT_FAILURE;
            }
        }
        std::cout << "SU2_ADDITIVE_COMPOUND_COMMON_SIGN result=PASS"
                  << " states=" << states.size()
                  << " constraints=" << constraints
                  << " maximum_state=" << maximum_state
                  << " maximum_label=" << maximum_label
                  << " minimum_label=" << minimum_label
                  << " level=" << level << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
