#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

using Matrix = std::vector<std::vector<cpp_int>>;

Matrix zero_matrix(int dimension) {
    return Matrix(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(static_cast<std::size_t>(dimension))
    );
}

Matrix identity_matrix(int dimension) {
    Matrix result = zero_matrix(dimension);
    for (int index = 0; index < dimension; ++index) {
        result[static_cast<std::size_t>(index)]
            [static_cast<std::size_t>(index)] = 1;
    }
    return result;
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

Matrix fusion_matrix(int level, int label) {
    const int dimension = level + 1;
    Matrix result = zero_matrix(dimension);
    for (int input = 0; input < dimension; ++input) {
        const int upper = std::min(input + label, 2 * level - input - label);
        for (int output = std::abs(input - label); output <= upper; output += 2) {
            result[static_cast<std::size_t>(output)]
                [static_cast<std::size_t>(input)] = 1;
        }
    }
    return result;
}

Matrix additive_update(const Matrix& wedge, const Matrix& operator_matrix) {
    const Matrix left = multiply(operator_matrix, wedge);
    const Matrix right = multiply(wedge, operator_matrix);
    const int dimension = static_cast<int>(wedge.size());
    Matrix result = zero_matrix(dimension);
    for (int row = 0; row < dimension; ++row) {
        for (int column = 0; column < dimension; ++column) {
            result[static_cast<std::size_t>(row)]
                [static_cast<std::size_t>(column)] =
                left[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
                + right[static_cast<std::size_t>(row)]
                    [static_cast<std::size_t>(column)];
        }
    }
    return result;
}

void visit_suffixes(
    int first_label,
    int maximum_label,
    int remaining,
    std::vector<int>& suffix,
    const auto& visit
) {
    visit(suffix);
    if (remaining == 0) {
        return;
    }
    for (int label = first_label; label <= maximum_label; ++label) {
        suffix.push_back(label);
        visit_suffixes(label, maximum_label, remaining - 1, suffix, visit);
        suffix.pop_back();
    }
}

cpp_int closed_boundary_coefficient(
    const Matrix& suffix, int level, int first_label, int target
) {
    cpp_int result = 0;
    for (int neighbor = std::abs(first_label - 1);
         neighbor <= std::min(first_label + 1, 2 * level - first_label - 1);
         neighbor += 2) {
        if (neighbor == 0) {
            continue;
        }
        result += suffix[static_cast<std::size_t>(target)]
            [static_cast<std::size_t>(neighbor)];
        if (target == neighbor) {
            result += suffix[0][0];
        }
    }
    if (first_label >= 2) {
        if (target == 1) {
            result += suffix[0][static_cast<std::size_t>(first_label)];
        }
        if (target == first_label) {
            result -= suffix[0][1];
        }
    }
    return result;
}

void print_labels(const std::vector<int>& labels) {
    std::cout << '[';
    for (std::size_t index = 0U; index < labels.size(); ++index) {
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
        if (argc == 3
            && std::string(argv[1]) == "--probe-two-first-mono") {
            const int maximum_level = std::stoi(argv[2]);
            if (maximum_level < 1) {
                throw std::runtime_error("invalid maximum level");
            }
            std::size_t checked_coefficients = 0U;
            for (int level = 1; level <= maximum_level; ++level) {
                const int dimension = level + 1;
                Matrix seed = zero_matrix(dimension);
                seed[1][0] = 1;
                seed[0][1] = -1;
                for (int first = 1; first <= level; ++first) {
                    const Matrix after_first = additive_update(
                        seed, fusion_matrix(level, first)
                    );
                    for (int second = first; second <= level; ++second) {
                        const Matrix after_second = additive_update(
                            after_first, fusion_matrix(level, second)
                        );
                        for (int suffix = 1; suffix <= level; ++suffix) {
                            const Matrix output = additive_update(
                                after_second, fusion_matrix(level, suffix)
                            );
                            for (int target = 1; target <= level; ++target) {
                                const cpp_int& value = output[
                                    static_cast<std::size_t>(target)][0];
                                if (value < 0) {
                                    std::cout
                                        << "SU2_TWO_FIRST_MONO result=NEGATIVE"
                                        << " level=" << level
                                        << " first=" << first
                                        << " second=" << second
                                        << " suffix=" << suffix
                                        << " target=" << target
                                        << " value=" << value << '\n';
                                    return EXIT_FAILURE;
                                }
                                ++checked_coefficients;
                            }
                        }
                    }
                }
            }
            std::cout << "SU2_TWO_FIRST_MONO result=PASS"
                      << " maximum_level=" << maximum_level
                      << " coefficients=" << checked_coefficients << '\n';
            return EXIT_SUCCESS;
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: verify_su2_first_prefix_hall_cut MAXIMUM_LEVEL "
                "MAXIMUM_SUFFIX_FACTORS | --probe-two-first-mono MAXIMUM_LEVEL"
            );
        }
        const int maximum_level = std::stoi(argv[1]);
        const int maximum_suffix_factors = std::stoi(argv[2]);
        if (maximum_level < 1 || maximum_suffix_factors < 0) {
            throw std::runtime_error("invalid bounds");
        }

        std::size_t checked_suffixes = 0U;
        std::size_t checked_coefficients = 0U;
        for (int level = 1; level <= maximum_level; ++level) {
            const int dimension = level + 1;
            for (int first_label = 1; first_label <= level; ++first_label) {
                std::vector<int> labels;
                visit_suffixes(1, level, maximum_suffix_factors, labels,
                    [&](const std::vector<int>& word) -> void {
                        Matrix suffix = identity_matrix(dimension);
                        for (const int label : word) {
                            suffix = multiply(fusion_matrix(level, label), suffix);
                        }
                        Matrix seed = zero_matrix(dimension);
                        seed[1][0] = 1;
                        seed[0][1] = -1;
                        const Matrix first_wedge = additive_update(
                            seed, fusion_matrix(level, first_label)
                        );
                        const Matrix direct = additive_update(first_wedge, suffix);
                        ++checked_suffixes;
                        for (int target = 1; target <= level; ++target) {
                            const cpp_int formula = closed_boundary_coefficient(
                                suffix, level, first_label, target
                            );
                            const cpp_int& direct_value = direct[
                                static_cast<std::size_t>(target)][0];
                            if (formula != direct_value || direct_value < 0) {
                                std::cout << "SU2_FIRST_PREFIX_HALL_CUT result=FAIL"
                                          << " level=" << level
                                          << " first=" << first_label
                                          << " suffix=";
                                print_labels(word);
                                std::cout << " target=" << target
                                          << " direct=" << direct_value
                                          << " formula=" << formula << '\n';
                                throw std::runtime_error(
                                    "first-prefix Hall formula or sign failure"
                                );
                            }
                            if (first_label >= 2 && target == first_label
                                && suffix[static_cast<std::size_t>(first_label)]
                                    [static_cast<std::size_t>(first_label - 1)]
                                    < suffix[0][1]) {
                                std::cout << "SU2_FIRST_PREFIX_HALL_CUT "
                                          << "result=FAIL_ANCHOR"
                                          << " level=" << level
                                          << " first=" << first_label
                                          << " suffix=";
                                print_labels(word);
                                std::cout << '\n';
                                throw std::runtime_error(
                                    "first-prefix Hall anchor failure"
                                );
                            }
                            ++checked_coefficients;
                        }
                    }
                );
            }
        }
        std::cout << "SU2_FIRST_PREFIX_HALL_CUT result=PASS"
                  << " maximum_level=" << maximum_level
                  << " maximum_suffix_factors=" << maximum_suffix_factors
                  << " suffixes=" << checked_suffixes
                  << " coefficients=" << checked_coefficients << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
