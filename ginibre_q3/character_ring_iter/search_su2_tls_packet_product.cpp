#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

using Matrix = std::vector<std::vector<cpp_int>>;

void for_each_output(int left, int right, const auto& function) {
    for (int output = std::abs(left - right); output <= left + right; output += 2) {
        function(output);
    }
}

Matrix zero_matrix(int dimension) {
    return Matrix(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(static_cast<std::size_t>(dimension))
    );
}

Matrix row_character(int label, int dimension) {
    Matrix answer = zero_matrix(dimension);
    if (label <= 0) {
        return answer;
    }
    for (int first = 0; first < label; ++first) {
        answer[static_cast<std::size_t>(first)]
            [static_cast<std::size_t>(label - 1 - first)] = 1;
    }
    return answer;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const int dimension = static_cast<int>(left.size());
    Matrix answer = zero_matrix(dimension);
    for (int a = 0; a < dimension; ++a) {
        for (int b = 0; b < dimension; ++b) {
            const cpp_int& left_value = left[static_cast<std::size_t>(a)]
                [static_cast<std::size_t>(b)];
            if (left_value == 0) {
                continue;
            }
            for (int c = 0; c < dimension; ++c) {
                for (int d = 0; d < dimension; ++d) {
                    const cpp_int& right_value = right[static_cast<std::size_t>(c)]
                        [static_cast<std::size_t>(d)];
                    if (right_value == 0) {
                        continue;
                    }
                    for_each_output(a, c, [&](int first) {
                        for_each_output(b, d, [&](int second) {
                            if (first < dimension && second < dimension) {
                                answer[static_cast<std::size_t>(first)]
                                    [static_cast<std::size_t>(second)] +=
                                    left_value * right_value;
                            }
                        });
                    });
                }
            }
        }
    }
    return answer;
}

cpp_int at(const Matrix& matrix, int row, int column) {
    if (row < 0 || column < 0 || row >= static_cast<int>(matrix.size())
        || column >= static_cast<int>(matrix.size())) {
        return 0;
    }
    return matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)];
}

Matrix packet_product(int plus_label, int left, int right) {
    const int dimension = plus_label + left + right + 5;
    Matrix packet = multiply(
        row_character(left, dimension), row_character(right, dimension)
    );
    const Matrix lower = multiply(
        row_character(left - 2, dimension),
        row_character(right - 2, dimension)
    );
    for (int first = 0; first < dimension; ++first) {
        for (int second = 0; second < dimension; ++second) {
            packet[static_cast<std::size_t>(first)]
                [static_cast<std::size_t>(second)] -=
                lower[static_cast<std::size_t>(first)]
                    [static_cast<std::size_t>(second)];
        }
    }
    Matrix plus = zero_matrix(dimension);
    plus[static_cast<std::size_t>(plus_label)][0] = 1;
    plus[0][static_cast<std::size_t>(plus_label)] = 1;
    return multiply(packet, plus);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: search_su2_tls_packet_product MAXIMUM_LABEL"
            );
        }
        const int maximum = std::stoi(argv[1]);
        if (maximum < 2) {
            throw std::runtime_error("MAXIMUM_LABEL must be at least two");
        }
        std::uint64_t tested = 0;
        for (int plus = 1; plus <= maximum; ++plus) {
            for (int left = 1; left <= maximum; ++left) {
                for (int right = 2; right <= left; ++right) {
                    if (plus == left || plus == right) {
                        continue;
                    }
                    const Matrix coefficients = packet_product(plus, left, right);
                    for (int first = 0;
                         first < static_cast<int>(coefficients.size()) - 1;
                         ++first) {
                        for (int second = 0; second <= first; ++second) {
                            const int row = first + 1;
                            const cpp_int value = at(coefficients, row - 1, second)
                                + at(coefficients, row + 1, second)
                                - at(coefficients, row, second - 1)
                                - at(coefficients, row, second + 1);
                            ++tested;
                            if (value < 0) {
                                std::cout << "SU2_TLS_PACKET_PRODUCT FAIL plus="
                                          << plus << " target=(" << left << ','
                                          << right << ") character=(" << first
                                          << ',' << second << ") coefficient="
                                          << value << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                    }
                }
            }
        }
        std::cout << "SU2_TLS_PACKET_PRODUCT PASS tested=" << tested
                  << " maximum_label=" << maximum << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
