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

Matrix fusion_matrix(int label, int level) {
    const int dimension = level + 1;
    Matrix matrix = zero_matrix(dimension);
    for (int input = 0; input < dimension; ++input) {
        for_each_output(input, label, level, [&](int output) {
            matrix[static_cast<std::size_t>(output)]
                [static_cast<std::size_t>(input)] = 1;
        });
    }
    return matrix;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const int dimension = static_cast<int>(left.size());
    Matrix result = zero_matrix(dimension);
    for (int row = 0; row < dimension; ++row) {
        for (int middle = 0; middle < dimension; ++middle) {
            const cpp_int& first = left[static_cast<std::size_t>(row)]
                [static_cast<std::size_t>(middle)];
            if (first == 0) {
                continue;
            }
            for (int column = 0; column < dimension; ++column) {
                const cpp_int& second = right[static_cast<std::size_t>(middle)]
                    [static_cast<std::size_t>(column)];
                if (second != 0) {
                    result[static_cast<std::size_t>(row)]
                        [static_cast<std::size_t>(column)] += first * second;
                }
            }
        }
    }
    return result;
}

Matrix additive_update(const Matrix& wedge, const Matrix& fusion) {
    const int dimension = static_cast<int>(wedge.size());
    Matrix result = zero_matrix(dimension);
    for (int row = 0; row < dimension; ++row) {
        for (int column = 0; column < dimension; ++column) {
            for (int middle = 0; middle < dimension; ++middle) {
                result[static_cast<std::size_t>(row)]
                    [static_cast<std::size_t>(column)] +=
                    fusion[static_cast<std::size_t>(row)]
                        [static_cast<std::size_t>(middle)]
                        * wedge[static_cast<std::size_t>(middle)]
                            [static_cast<std::size_t>(column)]
                    + wedge[static_cast<std::size_t>(row)]
                        [static_cast<std::size_t>(middle)]
                        * fusion[static_cast<std::size_t>(middle)]
                            [static_cast<std::size_t>(column)];
            }
        }
    }
    return result;
}

Matrix suffix_operator(
    const std::vector<int>& word,
    std::size_t first,
    int level
) {
    const int dimension = level + 1;
    Matrix result = zero_matrix(dimension);
    for (int state = 0; state < dimension; ++state) {
        result[static_cast<std::size_t>(state)]
            [static_cast<std::size_t>(state)] = 1;
    }
    for (std::size_t index = first; index < word.size(); ++index) {
        result = multiply(fusion_matrix(word[index], level), result);
    }
    return result;
}

cpp_int cut_coefficient(const Matrix& wedge, const Matrix& suffix, int target) {
    const Matrix first = multiply(suffix, wedge);
    const Matrix second = multiply(wedge, suffix);
    return first[static_cast<std::size_t>(target)][0]
        + second[static_cast<std::size_t>(target)][0];
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
                "usage: search_su2_level_edge_exterior MAXIMUM_LEVEL "
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
            std::size_t cut = 0;
            int target = -1;
            cpp_int value = 0;
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
                for (std::size_t word_index = 0; word_index < words.size(); ++word_index) {
                    if (local_failure.level >= 0 || words[word_index].empty()) {
                        continue;
                    }
                    const auto& word = words[word_index];
                    const int dimension = level + 1;
                    Matrix wedge = zero_matrix(dimension);
                    wedge[1][0] = 1;
                    wedge[0][1] = -1;
                    for (std::size_t cut = 0; cut <= word.size(); ++cut) {
                        const Matrix suffix = suffix_operator(word, cut, level);
                        for (int target = 2; target <= level; ++target) {
                            if (std::binary_search(word.begin(), word.end(), target)) {
                                continue;
                            }
                            const cpp_int value = cut_coefficient(
                                wedge, suffix, target
                            );
                            if (value < 0) {
                                local_failure.level = level;
                                local_failure.word = word;
                                local_failure.cut = cut;
                                local_failure.target = target;
                                local_failure.value = value;
                                break;
                            }
                        }
                        if (local_failure.level >= 0 || cut == word.size()) {
                            break;
                        }
                        wedge = additive_update(
                            wedge, fusion_matrix(word[cut], level)
                        );
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

        std::cout << "SU2_LEVEL_EDGE_EXTERIOR levels=2.." << maximum_level
                  << " maximum_factors=" << maximum_factors
                  << " tested_words=" << tested_words << '\n';
        if (failure.level >= 0) {
            std::cout << "FAIL level=" << failure.level << " plus=";
            print_labels(failure.word);
            std::cout << " cut=" << failure.cut
                      << " target=" << failure.target
                      << " coefficient=" << failure.value << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
