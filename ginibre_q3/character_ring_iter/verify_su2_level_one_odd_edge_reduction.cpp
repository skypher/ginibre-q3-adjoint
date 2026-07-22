#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Matrix = std::vector<std::vector<cpp_int>>;

Matrix zero_matrix(int dimension) {
    return Matrix(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(static_cast<std::size_t>(dimension))
    );
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

Matrix plus_update(int level, const Matrix& current, int label) {
    const int dimension = level + 1;
    Matrix result = zero_matrix(dimension);
    for (int first = 0; first < dimension; ++first) {
        for (int second = 0; second < dimension; ++second) {
            const cpp_int value = current[static_cast<std::size_t>(first)]
                [static_cast<std::size_t>(second)];
            if (value == 0) {
                continue;
            }
            for_each_output(level, first, label, [&](const int output) {
                result[static_cast<std::size_t>(output)]
                    [static_cast<std::size_t>(second)] += value;
            });
            for_each_output(level, second, label, [&](const int output) {
                result[static_cast<std::size_t>(first)]
                    [static_cast<std::size_t>(output)] += value;
            });
        }
    }
    return result;
}

Matrix exterior_update(int level, const Matrix& current, int label) {
    return plus_update(level, current, label);
}

cpp_int defect(
    int level,
    const Matrix& coefficients,
    int first,
    int second
) {
    cpp_int boundary = 0;
    for_each_output(level, first, second, [&](const int output) {
        boundary += coefficients[static_cast<std::size_t>(output)][0];
    });
    return boundary
        - coefficients[static_cast<std::size_t>(first)]
            [static_cast<std::size_t>(second)];
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
        if (argc != 3) {
            throw std::runtime_error(
                "usage: verify_su2_level_one_odd_edge_reduction "
                "MAXIMUM_LEVEL MAXIMUM_EVEN_FACTORS"
            );
        }
        const int maximum_level = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_level < 2 || maximum_factors < 0) {
            throw std::runtime_error("invalid bound");
        }
        long long words = 0;
        long long identities = 0;
        for (int level = 2; level <= maximum_level; ++level) {
            std::vector<int> word;
            const auto visit = [&](const auto& self, int next_even) -> void {
                ++words;
                const int dimension = level + 1;
                Matrix coefficients = zero_matrix(dimension);
                coefficients[0][0] = 1;
                Matrix exterior = zero_matrix(dimension);
                exterior[1][0] = 1;
                exterior[0][1] = -1;
                for (const int label : word) {
                    coefficients = plus_update(level, coefficients, label);
                    exterior = exterior_update(level, exterior, label);
                }
                for (int odd = 1; odd <= level; odd += 2) {
                    const Matrix extended = exterior_update(
                        level, exterior, odd
                    );
                    for (int target = 1; target <= level; ++target) {
                        ++identities;
                        const cpp_int left = extended
                            [static_cast<std::size_t>(target)][0];
                        cpp_int right = 0;
                        for_each_output(level, odd, 1, [&](const int neighbor) {
                            right += defect(
                                level, coefficients, target, neighbor
                            );
                        });
                        if (left != right) {
                            std::cout << "FAIL level=" << level
                                      << " even_word=";
                            print_word(word);
                            std::cout << " odd=" << odd
                                      << " target=" << target
                                      << " left=" << left
                                      << " right=" << right << '\n';
                            std::exit(EXIT_FAILURE);
                        }
                    }
                }
                if (word.size()
                    == static_cast<std::size_t>(maximum_factors)) {
                    return;
                }
                for (int label = next_even; label <= level; label += 2) {
                    word.push_back(label);
                    self(self, label);
                    word.pop_back();
                }
            };
            visit(visit, 2);
        }
        std::cout << "SU2_LEVEL_ONE_ODD_EDGE_REDUCTION PASS words=" << words
                  << " identities=" << identities
                  << " maximum_level=" << maximum_level
                  << " maximum_even_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
