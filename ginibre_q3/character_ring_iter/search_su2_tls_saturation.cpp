#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

using boost::multiprecision::cpp_int;

namespace {

using Flags = std::array<cpp_int, 4>;
using Matrix = std::vector<std::vector<Flags>>;

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
    visit(visit, 1, maximum_factors);
    return result;
}

Matrix flagged_coefficients(const std::vector<int>& labels) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    const int dimension = total + 3;
    Matrix current(
        static_cast<std::size_t>(dimension),
        std::vector<Flags>(static_cast<std::size_t>(dimension))
    );
    // A set bit means that the most recent update of that coordinate was
    // saturated.  Both coordinates are initially treated as saturated.
    current[0][0][3] = 1;
    for (int label : labels) {
        Matrix next(
            static_cast<std::size_t>(dimension),
            std::vector<Flags>(static_cast<std::size_t>(dimension))
        );
        for (int left = 0; left < dimension; ++left) {
            for (int right = 0; right < dimension; ++right) {
                for (int flags = 0; flags < 4; ++flags) {
                    const cpp_int& value = current[static_cast<std::size_t>(left)]
                        [static_cast<std::size_t>(right)]
                        [static_cast<std::size_t>(flags)];
                    if (value == 0) {
                        continue;
                    }
                    for_each_output(left, label, [&](int output) {
                        const int saturated = output == std::abs(left - label) ? 1 : 0;
                        const int next_flags = (flags & 2) | saturated;
                        next[static_cast<std::size_t>(output)]
                            [static_cast<std::size_t>(right)]
                            [static_cast<std::size_t>(next_flags)] += value;
                    });
                    for_each_output(right, label, [&](int output) {
                        const int saturated = output == std::abs(right - label) ? 2 : 0;
                        const int next_flags = (flags & 1) | saturated;
                        next[static_cast<std::size_t>(left)]
                            [static_cast<std::size_t>(output)]
                            [static_cast<std::size_t>(next_flags)] += value;
                    });
                }
            }
        }
        current = std::move(next);
    }
    return current;
}

cpp_int total_at(const Matrix& matrix, int left, int right) {
    if (left < 0 || right < 0 || left >= static_cast<int>(matrix.size())
        || right >= static_cast<int>(matrix.size())) {
        return 0;
    }
    cpp_int answer = 0;
    for (const cpp_int& value : matrix[static_cast<std::size_t>(left)]
                                      [static_cast<std::size_t>(right)]) {
        answer += value;
    }
    return answer;
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
        if (argc != 3 && argc != 4) {
            throw std::runtime_error(
                "usage: search_su2_tls_saturation MAXIMUM_LABEL MAXIMUM_FACTORS "
                "[--reverse|--total]"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 1 || maximum_factors < 0) {
            throw std::runtime_error("invalid search bound");
        }
        const bool reverse = argc == 4 && std::string(argv[3]) == "--reverse";
        const bool total_mode = argc == 4 && std::string(argv[3]) == "--total";
        if (argc == 4 && !reverse && !total_mode) {
            throw std::runtime_error("unknown ordering option");
        }
        auto words = multisets(maximum_label, maximum_factors);
        if (reverse) {
            for (auto& word : words) {
                std::reverse(word.begin(), word.end());
            }
        }
        std::size_t first_failure = words.size();
        int failure_left = -1;
        int failure_right = -1;
        cpp_int failure_residual = 0;
        cpp_int failure_boundary = 0;

#pragma omp parallel
        {
            std::size_t local_first = words.size();
            int local_left = -1;
            int local_right = -1;
            cpp_int local_residual = 0;
            cpp_int local_boundary = 0;

#pragma omp for schedule(dynamic)
            for (std::size_t index = 0; index < words.size(); ++index) {
                const Matrix matrix = flagged_coefficients(words[index]);
                bool failed = false;
                for (int left = 2;
                     left < static_cast<int>(matrix.size()) && !failed;
                     ++left) {
                    if (std::find(
                            words[index].begin(), words[index].end(), left
                        ) != words[index].end()) {
                        continue;
                    }
                    for (int right = 2; right <= left; ++right) {
                        if (std::find(
                                words[index].begin(), words[index].end(), right
                            ) != words[index].end()) {
                            continue;
                        }
                        const cpp_int total = total_at(matrix, left, right);
                        const cpp_int both_unsaturated =
                            matrix[static_cast<std::size_t>(left)]
                                [static_cast<std::size_t>(right)][0];
                        const cpp_int residual = total_mode
                            ? total
                            : total - both_unsaturated;
                        const cpp_int boundary = total_at(matrix, left + right, 0)
                            + total_at(matrix, left + right - 2, 0);
                        if (residual > boundary) {
                            if (index < local_first) {
                                local_first = index;
                                local_left = left;
                                local_right = right;
                                local_residual = residual;
                                local_boundary = boundary;
                            }
                            failed = true;
                            break;
                        }
                    }
                }
            }

#pragma omp critical
            {
                if (local_first < first_failure) {
                    first_failure = local_first;
                    failure_left = local_left;
                    failure_right = local_right;
                    failure_residual = local_residual;
                    failure_boundary = local_boundary;
                }
            }
        }

        std::cout << "SU2_TLS_SATURATION tested=" << words.size()
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors
                  << " order=" << (reverse ? "descending" : "ascending")
                  << " source=" << (total_mode ? "total" : "saturated")
                  << '\n';
        if (first_failure != words.size()) {
            std::cout << "FAIL plus=";
            print_labels(words[first_failure]);
            std::cout << " target=(" << failure_left << ',' << failure_right
                      << ") source_count=" << failure_residual
                      << " boundary=" << failure_boundary << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
