#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

using boost::multiprecision::cpp_int;

namespace {

using Matrix = std::vector<std::vector<cpp_int>>;

Matrix zero_matrix(int dimension) {
    return Matrix(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(static_cast<std::size_t>(dimension))
    );
}

void for_each_output(int left, int right, const auto& function) {
    for (int output = std::abs(left - right); output <= left + right; output += 2) {
        function(output);
    }
}

std::vector<std::vector<int>> multisets(int maximum_label, int maximum_factors) {
    std::vector<std::vector<int>> result;
    std::vector<int> current;
    const auto visit = [&](const auto& self, int first, int remaining) -> void {
        result.push_back(current);
        if (remaining == 0) {
            return;
        }
        for (int label = first; label <= maximum_label; ++label) {
            current.push_back(label);
            self(self, label, remaining - 1);
            current.pop_back();
        }
    };
    visit(visit, 2, maximum_factors);
    return result;
}

Matrix additive_update(const Matrix& wedge, int label) {
    const int dimension = static_cast<int>(wedge.size());
    Matrix result = zero_matrix(dimension);
    for (int first = 0; first < dimension; ++first) {
        for (int second = 0; second < dimension; ++second) {
            const cpp_int& value = wedge[static_cast<std::size_t>(first)]
                [static_cast<std::size_t>(second)];
            if (value == 0) {
                continue;
            }
            for_each_output(first, label, [&](int output) {
                if (output < dimension) {
                    result[static_cast<std::size_t>(output)]
                        [static_cast<std::size_t>(second)] += value;
                }
            });
            for_each_output(second, label, [&](int output) {
                if (output < dimension) {
                    result[static_cast<std::size_t>(first)]
                        [static_cast<std::size_t>(output)] += value;
                }
            });
        }
    }
    return result;
}

void add_wedge(Matrix& matrix, int first, int second, const cpp_int& value) {
    matrix[static_cast<std::size_t>(first)]
        [static_cast<std::size_t>(second)] += value;
    matrix[static_cast<std::size_t>(second)]
        [static_cast<std::size_t>(first)] -= value;
}

void print_labels(const std::vector<int>& labels) {
    std::cout << '[';
    for (std::size_t index = 0; index < labels.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << labels[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 3 || argc > 5) {
            throw std::runtime_error(
                "usage: search_su2_edge_packet_domination MAXIMUM_LABEL "
                "MAXIMUM_FACTORS [--all-targets] [--largest-first]"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 2 || maximum_factors < 1) {
            throw std::runtime_error("invalid search bound");
        }
        bool all_targets = false;
        bool largest_first = false;
        for (int argument = 3; argument < argc; ++argument) {
            const std::string option = argv[argument];
            if (option == "--all-targets") {
                all_targets = true;
            } else if (option == "--largest-first") {
                largest_first = true;
            } else {
                throw std::runtime_error("unknown mode");
            }
        }
        const auto words = multisets(maximum_label, maximum_factors);
        std::size_t first_failure = words.size();
        int failure_target = -1;
        cpp_int failure_value = 0;
        bool failure_is_outer = false;

#pragma omp parallel
        {
            std::size_t local_failure = words.size();
            int local_target = -1;
            cpp_int local_value = 0;
            bool local_is_outer = false;

#pragma omp for schedule(dynamic)
            for (std::size_t index = 0; index < words.size(); ++index) {
                const auto& word = words[index];
                if (word.empty()) {
                    continue;
                }
                int total = 0;
                for (int label : word) {
                    total += label;
                }
                const int dimension = 2 * total + 5;
                const std::size_t first_position = largest_first
                    ? word.size() - 1U : 0U;
                const int first_label = word[first_position];
                Matrix difference = zero_matrix(dimension);
                add_wedge(difference, first_label - 1, 0, 1);
                add_wedge(difference, first_label, 1, -1);
                Matrix outer = zero_matrix(dimension);
                add_wedge(outer, first_label + 1, 0, 1);
                for (std::size_t position = 0; position < word.size(); ++position) {
                    if (position == first_position) {
                        continue;
                    }
                    difference = additive_update(difference, word[position]);
                    outer = additive_update(outer, word[position]);
                }
                for (int target = 2; target < dimension; ++target) {
                    if (!all_targets
                        && std::binary_search(word.begin(), word.end(), target)) {
                        continue;
                    }
                    const cpp_int& lower_value = difference
                        [static_cast<std::size_t>(target)][0];
                    const cpp_int& outer_value = outer
                        [static_cast<std::size_t>(target)][0];
                    if (lower_value < 0 || outer_value < 0) {
                        if (index < local_failure) {
                            local_failure = index;
                            local_target = target;
                            local_is_outer = outer_value < 0;
                            local_value = local_is_outer
                                ? outer_value : lower_value;
                        }
                        break;
                    }
                }
            }

#pragma omp critical
            {
                if (local_failure < first_failure) {
                    first_failure = local_failure;
                    failure_target = local_target;
                    failure_value = local_value;
                    failure_is_outer = local_is_outer;
                }
            }
        }

        std::cout << "SU2_EDGE_PACKET_DOMINATION tested_words=" << words.size()
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors << '\n';
        if (first_failure != words.size()) {
            std::cout << "FAIL plus=";
            print_labels(words[first_failure]);
            std::cout << " target=" << failure_target
                      << (failure_is_outer
                              ? " outer_row=" : " lower_minus_crossing=")
                      << failure_value << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
