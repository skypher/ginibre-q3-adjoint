#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Character = std::vector<cpp_int>;
using Matrix = std::vector<Character>;

void for_each_output(int left, int right, const auto& function) {
    for (int output = std::abs(left - right);
         output <= left + right;
         output += 2) {
        function(output);
    }
}

Matrix zero_matrix(int dimension) {
    return Matrix(
        static_cast<std::size_t>(dimension),
        Character(static_cast<std::size_t>(dimension))
    );
}

Matrix add_plus(const Matrix& current, int label) {
    const int dimension = static_cast<int>(current.size());
    Matrix result = zero_matrix(dimension);
    for (int first = 0; first < dimension; ++first) {
        for (int second = 0; second < dimension; ++second) {
            const cpp_int value = current[static_cast<std::size_t>(first)]
                [static_cast<std::size_t>(second)];
            if (value == 0) {
                continue;
            }
            for_each_output(first, label, [&](const int output) {
                if (output < dimension) {
                    result[static_cast<std::size_t>(output)]
                        [static_cast<std::size_t>(second)] += value;
                }
            });
            for_each_output(second, label, [&](const int output) {
                if (output < dimension) {
                    result[static_cast<std::size_t>(first)]
                        [static_cast<std::size_t>(output)] += value;
                }
            });
        }
    }
    return result;
}

Character column(const Matrix& matrix, int second) {
    const int dimension = static_cast<int>(matrix.size());
    Character result(static_cast<std::size_t>(dimension));
    if (second < 0 || second >= dimension) {
        return result;
    }
    for (int first = 0; first < dimension; ++first) {
        result[static_cast<std::size_t>(first)]
            = matrix[static_cast<std::size_t>(first)]
                [static_cast<std::size_t>(second)];
    }
    return result;
}

Character multiply(const Character& input, int label) {
    const int dimension = static_cast<int>(input.size());
    Character result(static_cast<std::size_t>(dimension));
    for (int first = 0; first < dimension; ++first) {
        const cpp_int value = input[static_cast<std::size_t>(first)];
        if (value == 0) {
            continue;
        }
        for_each_output(first, label, [&](const int output) {
            if (output < dimension) {
                result[static_cast<std::size_t>(output)] += value;
            }
        });
    }
    return result;
}

void add(Character& target, const Character& source, int sign) {
    for (std::size_t index = 0; index < target.size(); ++index) {
        target[index] += sign * source[index];
    }
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
        if (argc != 4) {
            throw std::runtime_error(
                "usage: analyze_su2_eoh_decomposition MAXIMUM_EVEN_LABEL "
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
        long long cases = 0;
        bool first_packet_nonnegative = true;
        bool second_packet_nonnegative = true;
        bool defect_reservoir_nonnegative = true;
        for (int p = 2; p <= maximum_even; p += 2) {
            std::vector<int> word;
            const auto visit = [&](const auto& self, int next_even) -> void {
                ++words;
                int total = 0;
                for (const int label : word) {
                    total += label;
                }
                const int dimension = total + p + maximum_odd
                    + maximum_even + 8;
                Matrix coefficients = zero_matrix(dimension);
                coefficients[0][0] = 1;
                for (const int label : word) {
                    coefficients = add_plus(coefficients, label);
                }
                const Character c_zero = column(coefficients, 0);
                const Character c_p = column(coefficients, p);
                for (int q = p | 1; q <= maximum_odd; q += 2) {
                    ++cases;
                    Character first(static_cast<std::size_t>(dimension));
                    for_each_output(q, p - 1, [&](const int b) {
                        Character defect = multiply(c_zero, b);
                        add(defect, column(coefficients, b), -1);
                        add(first, defect, 1);
                    });
                    Character inner = multiply(c_p, q - 1);
                    add(inner, multiply(c_p, q + 1), 1);
                    const Character c_q_minus = column(coefficients, q - 1);
                    const Character c_q_plus = column(coefficients, q + 1);
                    add(inner, c_q_minus, -1);
                    add(inner, c_q_plus, -1);
                    const Character second = multiply(inner, p);
                    Character negative_second = c_q_minus;
                    add(negative_second, c_q_plus, 1);
                    negative_second = multiply(negative_second, p);
                    Character defect_reservoir = first;
                    add(defect_reservoir, negative_second, -1);
                    Character total_packet = first;
                    add(total_packet, second, 1);
                    for (int target = 0; target < dimension; ++target) {
                        const std::size_t index = static_cast<std::size_t>(target);
                        if (first[index] < 0 && first_packet_nonnegative) {
                            first_packet_nonnegative = false;
                            std::cout << "FIRST_PACKET_FAIL p=" << p
                                      << " q=" << q << " even_word=";
                            print_word(word);
                            std::cout << " target=" << target
                                      << " value=" << first[index] << '\n';
                        }
                        if (second[index] < 0 && second_packet_nonnegative) {
                            second_packet_nonnegative = false;
                            std::cout << "SECOND_PACKET_FAIL p=" << p
                                      << " q=" << q << " even_word=";
                            print_word(word);
                            std::cout << " target=" << target
                                      << " value=" << second[index]
                                      << " first_packet=" << first[index]
                                      << " total=" << total_packet[index]
                                      << '\n';
                        }
                        if (target > 0 && defect_reservoir[index] < 0
                            && defect_reservoir_nonnegative) {
                            defect_reservoir_nonnegative = false;
                            std::cout << "DEFECT_RESERVOIR_FAIL p=" << p
                                      << " q=" << q << " even_word=";
                            print_word(word);
                            std::cout << " target=" << target
                                      << " value=" << defect_reservoir[index]
                                      << '\n';
                        }
                        if (total_packet[index] < 0) {
                            std::cout << "TOTAL_FAIL p=" << p << " q=" << q
                                      << " even_word=";
                            print_word(word);
                            std::cout << " target=" << target
                                      << " value=" << total_packet[index]
                                      << '\n';
                            std::exit(EXIT_FAILURE);
                        }
                    }
                }
                if (word.size()
                    == static_cast<std::size_t>(maximum_factors)) {
                    return;
                }
                for (int label = next_even;
                     label <= maximum_even;
                     label += 2) {
                    word.push_back(label);
                    self(self, label);
                    word.pop_back();
                }
            };
            visit(visit, p);
        }
        std::cout << "SU2_EOH_DECOMPOSITION PASS words=" << words
                  << " cases=" << cases
                  << " first_packet_nonnegative="
                  << (first_packet_nonnegative ? "true" : "false")
                  << " second_packet_nonnegative="
                  << (second_packet_nonnegative ? "true" : "false")
                  << " defect_reservoir_nonnegative="
                  << (defect_reservoir_nonnegative ? "true" : "false")
                  << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
