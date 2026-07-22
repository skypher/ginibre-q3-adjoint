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

void for_each_output(int left, int right, const auto& function) {
    for (int output = std::abs(left - right); output <= left + right;
         output += 2) {
        function(output);
    }
}

bool fusion_contains(int left, int right, int output) {
    return output >= std::abs(left - right) && output <= left + right
        && ((output - left - right) & 1) == 0;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const int dimension = static_cast<int>(left.size());
    Matrix result = zero_matrix(dimension);
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
                                result[static_cast<std::size_t>(first)]
                                    [static_cast<std::size_t>(second)] +=
                                    left_value * right_value;
                            }
                        });
                    });
                }
            }
        }
    }
    return result;
}

Matrix row_restriction(int label, int dimension) {
    Matrix result = zero_matrix(dimension);
    for (int first = 0; first <= label; ++first) {
        result[static_cast<std::size_t>(first)]
            [static_cast<std::size_t>(label - first)] = 1;
    }
    return result;
}

Matrix mixed_part(int label, int dimension) {
    Matrix result = zero_matrix(dimension);
    for (int first = 1; first < label; ++first) {
        result[static_cast<std::size_t>(first)]
            [static_cast<std::size_t>(label - first)] = 1;
    }
    return result;
}

Matrix endpoint_update(const Matrix& input, int label) {
    const int dimension = static_cast<int>(input.size());
    Matrix result = zero_matrix(dimension);
    for (int first = 0; first < dimension; ++first) {
        for (int second = 0; second < dimension; ++second) {
            const cpp_int& value = input[static_cast<std::size_t>(first)]
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

cpp_int wedge_packet_value(
    int p, int a, const std::vector<int>& suffix, int dimension
) {
    Matrix packet = zero_matrix(dimension);
    add_wedge(packet, p - 1, 0, 1);
    add_wedge(packet, p, 1, -1);
    for (int label : suffix) {
        packet = endpoint_update(packet, label);
    }
    return packet[static_cast<std::size_t>(a)][0];
}

cpp_int short_suffix_formula(
    int p, int a, const std::vector<int>& suffix, int dimension
) {
    if (suffix.empty()) {
        return a == p - 1 ? cpp_int{1} : cpp_int{0};
    }
    const int q = suffix[0];
    if (suffix.size() == 1U) {
        cpp_int answer = fusion_contains(q, p - 1, a)
            ? cpp_int{1} : cpp_int{0};
        if (q == p && a == 1) {
            ++answer;
        }
        return answer;
    }
    if (suffix.size() != 2U) {
        throw std::runtime_error("short suffix formula requires at most two labels");
    }
    const int r = suffix[1];
    cpp_int answer = 0;
    for (int channel = 0; channel < dimension; ++channel) {
        if (fusion_contains(a, r, channel)
            && fusion_contains(q, p - 1, channel)) {
            ++answer;
        }
    }
    if (q == p && fusion_contains(a, r, 1)) {
        ++answer;
    }
    if (r == q && a == p - 1) {
        ++answer;
    }
    if (a == 1 && fusion_contains(q, p, r)) {
        ++answer;
    }
    if (a == p && fusion_contains(q, 1, r)) {
        --answer;
    }
    if (r == p && fusion_contains(q, 1, a)) {
        ++answer;
    }
    return answer;
}

cpp_int k_index_value(
    int p, int a, const std::vector<int>& suffix, int dimension,
    Matrix* output
) {
    Matrix module = multiply(
        row_restriction(a - 1, dimension), mixed_part(p, dimension)
    );
    for (int label : suffix) {
        module = endpoint_update(module, label);
    }
    const cpp_int index = module[1][1] - module[0][0] - module[2][0];
    if (output != nullptr) {
        *output = std::move(module);
    }
    return index;
}

std::vector<std::vector<int>> suffixes(int p, int maximum_label, int length) {
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
    visit(visit, p, length);
    return result;
}

void print_suffix(const std::vector<int>& suffix) {
    std::cout << '[';
    for (std::size_t index = 0; index < suffix.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << suffix[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc >= 2 && std::string(argv[1]) == "--verify") {
            if (argc != 4) {
                throw std::runtime_error(
                    "usage: inspect_su2_edge_packet_kindex --verify "
                    "MAXIMUM_LABEL MAXIMUM_SUFFIX_LENGTH"
                );
            }
            const int maximum_label = std::stoi(argv[2]);
            const int maximum_length = std::stoi(argv[3]);
            if (maximum_label < 2 || maximum_length < 0) {
                throw std::runtime_error("invalid verification bound");
            }
            std::uint64_t tested = 0;
            for (int p = 2; p <= maximum_label; ++p) {
                for (const auto& suffix : suffixes(
                         p, maximum_label, maximum_length
                     )) {
                    int suffix_sum = 0;
                    for (int label : suffix) {
                        suffix_sum += label;
                    }
                    const int maximum_target = p + suffix_sum + 1;
                    const int dimension = 2 * maximum_target + 6;
                    for (int a = 1; a <= maximum_target; ++a) {
                        const cpp_int wedge = wedge_packet_value(
                            p, a, suffix, dimension
                        );
                        const cpp_int index = k_index_value(
                            p, a, suffix, dimension, nullptr
                        );
                        if (wedge != index) {
                            std::cout << "EDGE_PACKET_KINDEX FAIL p=" << p
                                      << " a=" << a << " suffix=";
                            print_suffix(suffix);
                            std::cout << " wedge=" << wedge
                                      << " index=" << index << '\n';
                            return EXIT_FAILURE;
                        }
                        if (suffix.size() <= 2U) {
                            const cpp_int formula = short_suffix_formula(
                                p, a, suffix, dimension
                            );
                            if (wedge != formula) {
                                std::cout << "EDGE_PACKET_SHORT_FORMULA FAIL p="
                                          << p << " a=" << a << " suffix=";
                                print_suffix(suffix);
                                std::cout << " wedge=" << wedge
                                          << " formula=" << formula << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                        ++tested;
                    }
                }
            }
            std::cout << "EDGE_PACKET_KINDEX PASS tested=" << tested
                      << " maximum_label=" << maximum_label
                      << " maximum_suffix_length=" << maximum_length
                      << '\n';
            return EXIT_SUCCESS;
        }
        if (argc < 3) {
            throw std::runtime_error(
                "usage: inspect_su2_edge_packet_kindex P A [SUFFIX ...]"
            );
        }
        const int p = std::stoi(argv[1]);
        const int a = std::stoi(argv[2]);
        if (p < 2 || a < 1) {
            throw std::runtime_error("invalid packet or target label");
        }
        std::vector<int> suffix;
        int suffix_sum = 0;
        for (int argument = 3; argument < argc; ++argument) {
            const int label = std::stoi(argv[argument]);
            if (label < p) {
                throw std::runtime_error("suffix labels must be at least P");
            }
            suffix.push_back(label);
            suffix_sum += label;
        }
        const int dimension = 2 * (a + p + suffix_sum) + 6;
        Matrix module;
        const cpp_int index = k_index_value(
            p, a, suffix, dimension, &module
        );
        const cpp_int wedge = wedge_packet_value(p, a, suffix, dimension);
        std::cout << "EDGE_PACKET_KINDEX p=" << p << " a=" << a
                  << " suffix=";
        print_suffix(suffix);
        std::cout << '\n'
                  << "m00=" << module[0][0]
                  << " m20=" << module[2][0]
                  << " m11=" << module[1][1] << '\n'
                  << "index=" << index << " wedge=" << wedge << '\n';
        const int slack_limit = std::min(10, dimension - 2);
        for (int diagonal = 1; diagonal <= slack_limit; ++diagonal) {
            const cpp_int slack =
                module[static_cast<std::size_t>(diagonal)]
                    [static_cast<std::size_t>(diagonal)]
                - module[static_cast<std::size_t>(diagonal - 1)]
                    [static_cast<std::size_t>(diagonal - 1)]
                - module[static_cast<std::size_t>(diagonal + 1)]
                    [static_cast<std::size_t>(diagonal - 1)];
            std::cout << "diagonal_slack[" << diagonal << "]=" << slack
                      << '\n';
        }
        return index == wedge ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
