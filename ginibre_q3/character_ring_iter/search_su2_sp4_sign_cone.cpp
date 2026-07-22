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

Matrix plus_coefficients(const std::vector<int>& word) {
    int total = 0;
    for (int label : word) {
        total += label;
    }
    const int dimension = total + 4;
    Matrix current(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(static_cast<std::size_t>(dimension))
    );
    current[0][0] = 1;
    for (int label : word) {
        Matrix next(
            static_cast<std::size_t>(dimension),
            std::vector<cpp_int>(static_cast<std::size_t>(dimension))
        );
        for (int first = 0; first < dimension; ++first) {
            for (int second = 0; second < dimension; ++second) {
                const cpp_int& value = current[static_cast<std::size_t>(first)]
                    [static_cast<std::size_t>(second)];
                if (value == 0) {
                    continue;
                }
                for_each_output(first, label, [&](int output) {
                    if (output < dimension) {
                        next[static_cast<std::size_t>(output)]
                            [static_cast<std::size_t>(second)] += value;
                    }
                });
                for_each_output(second, label, [&](int output) {
                    if (output < dimension) {
                        next[static_cast<std::size_t>(first)]
                            [static_cast<std::size_t>(output)] += value;
                    }
                });
            }
        }
        current = std::move(next);
    }
    return current;
}

cpp_int at(const Matrix& matrix, int first, int second) {
    if (first < 0 || second < 0 || first >= static_cast<int>(matrix.size())
        || second >= static_cast<int>(matrix.size())) {
        return 0;
    }
    return matrix[static_cast<std::size_t>(first)]
        [static_cast<std::size_t>(second)];
}

cpp_int sp4_coefficient(const Matrix& coefficients, int first, int second) {
    const int row = first + 1;
    return at(coefficients, row - 1, second)
        + at(coefficients, row + 1, second)
        - at(coefficients, row, second - 1)
        - at(coefficients, row, second + 1);
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
                "usage: search_su2_sp4_sign_cone MAXIMUM_LABEL MAXIMUM_FACTORS"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        const auto words = multisets(maximum_label, maximum_factors);
        for (const auto& word : words) {
            const Matrix coefficients = plus_coefficients(word);
            for (int first = 0; first + 1 < static_cast<int>(coefficients.size()); ++first) {
                for (int second = 0; second <= first; ++second) {
                    const cpp_int value = sp4_coefficient(
                        coefficients, first, second
                    );
                    const cpp_int signed_value = (second & 1) == 0
                        ? value : -value;
                    if (signed_value < 0) {
                        std::cout << "SU2_SP4_SIGN_CONE FAIL plus=";
                        print_labels(word);
                        std::cout << " weight=(" << first << ',' << second
                                  << ") coefficient=" << value << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
        }
        std::cout << "SU2_SP4_SIGN_CONE PASS tested_words=" << words.size()
                  << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
