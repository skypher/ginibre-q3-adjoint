#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
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

Matrix row_character(int d_label, int dimension) {
    Matrix result = zero_matrix(dimension);
    for (int first = 0; first < d_label; ++first) {
        result[static_cast<std::size_t>(first)]
            [static_cast<std::size_t>(d_label - 1 - first)] = 1;
    }
    return result;
}

Matrix packet_y(int p, int dimension) {
    Matrix result = zero_matrix(dimension);
    for (int first = 1; first < p; ++first) {
        result[static_cast<std::size_t>(first)]
            [static_cast<std::size_t>(p - first)] = -1;
    }
    return result;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const int dimension = static_cast<int>(left.size());
    Matrix result = zero_matrix(dimension);
    for (int a = 0; a < dimension; ++a) {
        for (int b = 0; b < dimension; ++b) {
            const cpp_int& first_value = left[static_cast<std::size_t>(a)]
                [static_cast<std::size_t>(b)];
            if (first_value == 0) {
                continue;
            }
            for (int c = 0; c < dimension; ++c) {
                for (int d = 0; d < dimension; ++d) {
                    const cpp_int& second_value = right[static_cast<std::size_t>(c)]
                        [static_cast<std::size_t>(d)];
                    if (second_value == 0) {
                        continue;
                    }
                    for_each_output(a, c, [&](int first) {
                        for_each_output(b, d, [&](int second) {
                            result[static_cast<std::size_t>(first)]
                                [static_cast<std::size_t>(second)] +=
                                first_value * second_value;
                        });
                    });
                }
            }
        }
    }
    return result;
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

bool in_dual_packet(int p, int a, int first, int second) {
    if (second < 1 || second > first) {
        return false;
    }
    const bool first_family = a >= p + 1 && second <= p - 1
        && first == a + second - p - 1;
    const int twice_u = p + first + second - a - 1;
    if ((twice_u & 1) != 0) {
        return first_family;
    }
    const int u = twice_u / 2;
    const bool second_family = u >= 1 && u <= p - 1
        && second <= u && u <= first
        && (first > u || u == p - 1);
    return first_family || second_family;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4 && std::string(argv[1]) == "--stencil") {
            const int p = std::stoi(argv[2]);
            const int a = std::stoi(argv[3]);
            if (p < 2 || a < 1) {
                throw std::runtime_error("invalid labels");
            }
            const int dimension = p + a + 5;
            std::map<std::pair<int, int>, int> stencil;
            for (int first = 0; first < dimension; ++first) {
                for (int second = 1; second <= first; ++second) {
                    if (!in_dual_packet(p, a, first, second)) {
                        continue;
                    }
                    ++stencil[{first, second}];
                    ++stencil[{first + 2, second}];
                    --stencil[{first + 1, second - 1}];
                    --stencil[{first + 1, second + 1}];
                }
            }
            std::cout << "PACKET_DUAL_STENCIL p=" << p << " a=" << a
                      << '\n';
            for (const auto& [weight, coefficient] : stencil) {
                if (coefficient != 0) {
                    std::cout << '(' << weight.first << ',' << weight.second
                              << ") " << coefficient << '\n';
                }
            }
            return EXIT_SUCCESS;
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: inspect_su2_edge_packet_dual P A | --verify MAXIMUM "
                "| --stencil P A"
            );
        }
        if (std::string(argv[1]) == "--verify") {
            const int maximum = std::stoi(argv[2]);
            std::uint64_t tested = 0;
            for (int p = 2; p <= maximum; ++p) {
                for (int a = 1; a <= maximum; ++a) {
                    const int dimension = p + a + 5;
                    const Matrix product = multiply(
                        row_character(a, dimension), packet_y(p, dimension)
                    );
                    bool positive_seen = false;
                    for (int first = 0; first + 1 < dimension; ++first) {
                        for (int second = 0; second <= first; ++second) {
                            const cpp_int value = sp4_coefficient(
                                product, first, second
                            );
                            if (value == 0) {
                                continue;
                            }
                            const bool expected_positive = a < p
                                && first == p - a - 1 && second == 0;
                            const bool expected_negative = in_dual_packet(
                                p, a, first, second
                            );
                            if (expected_positive && value == 1) {
                                positive_seen = true;
                            } else if (!(expected_negative && value == -1)) {
                                std::cout << "PACKET_DUAL_VERIFY FAIL p=" << p
                                          << " a=" << a << " weight=("
                                          << first << ',' << second
                                          << ") coefficient=" << value
                                          << " expected_negative="
                                          << expected_negative << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                    }
                    if (positive_seen != (a < p)) {
                        std::cout << "PACKET_DUAL_VERIFY FAIL positive_term p="
                                  << p << " a=" << a << '\n';
                        return EXIT_FAILURE;
                    }
                    for (int first = 0; first + 1 < dimension; ++first) {
                        for (int second = 1; second <= first; ++second) {
                            if (!in_dual_packet(p, a, first, second)) {
                                continue;
                            }
                            if (sp4_coefficient(product, first, second) != -1) {
                                std::cout
                                    << "PACKET_DUAL_VERIFY FAIL missing weight p="
                                    << p << " a=" << a << " weight=("
                                    << first << ',' << second << ")\n";
                                return EXIT_FAILURE;
                            }
                        }
                    }
                    ++tested;
                }
            }
            std::cout << "PACKET_DUAL_VERIFY PASS tested=" << tested
                      << " maximum=" << maximum << '\n';
            return EXIT_SUCCESS;
        }
        const int p = std::stoi(argv[1]);
        const int a = std::stoi(argv[2]);
        if (p < 2 || a < 1) {
            throw std::runtime_error("invalid labels");
        }
        const int dimension = p + a + 5;
        const Matrix product = multiply(
            row_character(a, dimension), packet_y(p, dimension)
        );
        std::cout << "PACKET_DUAL p=" << p << " a=" << a << '\n';
        for (int first = 0; first + 1 < dimension; ++first) {
            for (int second = 0; second <= first; ++second) {
                const cpp_int value = sp4_coefficient(product, first, second);
                if (value != 0) {
                    std::cout << '(' << first << ',' << second << ") "
                              << value << '\n';
                }
            }
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
