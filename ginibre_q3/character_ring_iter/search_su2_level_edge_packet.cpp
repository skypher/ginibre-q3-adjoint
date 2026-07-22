#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
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

void for_each_output(int left, int right, int level, const auto& function) {
    const int last = std::min(left + right, 2 * level - left - right);
    for (int output = std::abs(left - right); output <= last; output += 2) {
        function(output);
    }
}

std::vector<std::vector<int>> multisets(
    int first_label,
    int maximum_label,
    int maximum_factors
) {
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
    visit(visit, first_label, maximum_factors);
    return result;
}

Matrix additive_update(const Matrix& wedge, int label, int level) {
    const int dimension = level + 1;
    Matrix result = zero_matrix(dimension);
    for (int first = 0; first < dimension; ++first) {
        for (int second = 0; second < dimension; ++second) {
            const cpp_int& value = wedge[static_cast<std::size_t>(first)]
                [static_cast<std::size_t>(second)];
            if (value == 0) {
                continue;
            }
            for_each_output(first, label, level, [&](int output) {
                result[static_cast<std::size_t>(output)]
                    [static_cast<std::size_t>(second)] += value;
            });
            for_each_output(second, label, level, [&](int output) {
                result[static_cast<std::size_t>(first)]
                    [static_cast<std::size_t>(output)] += value;
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
        if (argc != 3) {
            throw std::runtime_error(
                "usage: search_su2_level_edge_packet MAXIMUM_LEVEL "
                "MAXIMUM_FACTORS"
            );
        }
        const int maximum_level = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_level < 2 || maximum_factors < 1) {
            throw std::runtime_error("invalid search bound");
        }

        struct Failure {
            int level = -1;
            std::vector<int> word;
            int target = -1;
            cpp_int value = 0;
            cpp_int companion = 0;
            bool outer = false;
        };
        Failure failure;
        std::size_t tested_words = 0;

        for (int level = 2; level <= maximum_level && failure.level < 0; ++level) {
            const auto words = multisets(2, level, maximum_factors);
            tested_words += words.size();
            Failure level_failure;

#pragma omp parallel
            {
                Failure local_failure;

#pragma omp for schedule(dynamic)
                for (std::size_t index = 0; index < words.size(); ++index) {
                    const auto& word = words[index];
                    if (local_failure.level >= 0 || word.empty()) {
                        continue;
                    }
                    const int p = word.front();
                    const int dimension = level + 1;
                    Matrix lower = zero_matrix(dimension);
                    add_wedge(lower, p - 1, 0, 1);
                    add_wedge(lower, p, 1, -1);
                    Matrix outer = zero_matrix(dimension);
                    if (p < level) {
                        add_wedge(outer, p + 1, 0, 1);
                    }
                    for (std::size_t position = 1; position < word.size(); ++position) {
                        lower = additive_update(lower, word[position], level);
                        outer = additive_update(outer, word[position], level);
                    }
                    for (int target = 2; target <= level; ++target) {
                        const cpp_int& lower_value = lower
                            [static_cast<std::size_t>(target)][0];
                        const cpp_int& outer_value = outer
                            [static_cast<std::size_t>(target)][0];
                        if (lower_value < 0 || outer_value < 0) {
                            local_failure.level = level;
                            local_failure.word = word;
                            local_failure.target = target;
                            local_failure.outer = outer_value < 0;
                            local_failure.value = local_failure.outer
                                ? outer_value : lower_value;
                            local_failure.companion = local_failure.outer
                                ? lower_value : outer_value;
                            break;
                        }
                    }
                }

#pragma omp critical
                {
                    if (local_failure.level >= 0 && level_failure.level < 0) {
                        level_failure = std::move(local_failure);
                    }
                }
            }
            if (level_failure.level >= 0) {
                failure = std::move(level_failure);
            }
        }

        std::cout << "SU2_LEVEL_EDGE_PACKET levels=2.." << maximum_level
                  << " maximum_factors=" << maximum_factors
                  << " tested_words=" << tested_words << '\n';
        if (failure.level >= 0) {
            std::cout << "FAIL level=" << failure.level << " plus=";
            print_labels(failure.word);
            std::cout << " target=" << failure.target
                      << (failure.outer
                              ? " outer_row=" : " lower_minus_crossing=")
                      << failure.value
                      << " companion=" << failure.companion
                      << " total=" << failure.value + failure.companion
                      << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
