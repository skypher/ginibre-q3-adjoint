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

Matrix plus_coefficients(const std::vector<int>& labels, int dimension) {
    Matrix current(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(static_cast<std::size_t>(dimension))
    );
    current[0][0] = 1;
    for (int label : labels) {
        Matrix next(
            static_cast<std::size_t>(dimension),
            std::vector<cpp_int>(static_cast<std::size_t>(dimension))
        );
        for (int left = 0; left < dimension; ++left) {
            for (int right = 0; right < dimension; ++right) {
                const cpp_int& coefficient = current[static_cast<std::size_t>(left)]
                    [static_cast<std::size_t>(right)];
                if (coefficient == 0) {
                    continue;
                }
                for_each_output(left, label, [&](int output) {
                    if (output < dimension) {
                        next[static_cast<std::size_t>(output)]
                            [static_cast<std::size_t>(right)] += coefficient;
                    }
                });
                for_each_output(right, label, [&](int output) {
                    if (output < dimension) {
                        next[static_cast<std::size_t>(left)]
                            [static_cast<std::size_t>(output)] += coefficient;
                    }
                });
            }
        }
        current = std::move(next);
    }
    return current;
}

std::vector<cpp_int> tensor_multiplicities(
    const std::vector<int>& labels,
    int dimension
) {
    std::vector<cpp_int> current(static_cast<std::size_t>(dimension));
    current[0] = 1;
    for (int label : labels) {
        std::vector<cpp_int> next(static_cast<std::size_t>(dimension));
        for (int input = 0; input < dimension; ++input) {
            if (current[static_cast<std::size_t>(input)] == 0) {
                continue;
            }
            for_each_output(input, label, [&](int output) {
                if (output < dimension) {
                    next[static_cast<std::size_t>(output)]
                        += current[static_cast<std::size_t>(input)];
                }
            });
        }
        current = std::move(next);
    }
    return current;
}

Matrix defect(const std::vector<int>& labels, int minimum_dimension = 0) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    // The extra diagonal range detects stabilization of every tail whose
    // difference can meet the boundary support 0,...,total.
    const int dimension = std::max(2 * total + 5, minimum_dimension);
    const Matrix coefficients = plus_coefficients(labels, dimension);
    std::vector<cpp_int> boundary(static_cast<std::size_t>(dimension));
    for (int label = 0; label < dimension; ++label) {
        boundary[static_cast<std::size_t>(label)] =
            coefficients[static_cast<std::size_t>(label)][0];
    }

    Matrix answer(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(static_cast<std::size_t>(dimension))
    );
    for (int left = 0; left < dimension; ++left) {
        for (int right = 0; right < dimension; ++right) {
            cpp_int hankel = 0;
            for_each_output(left, right, [&](int output) {
                if (output < dimension) {
                    hankel += boundary[static_cast<std::size_t>(output)];
                }
            });
            answer[static_cast<std::size_t>(left)][static_cast<std::size_t>(right)] =
                hankel - coefficients[static_cast<std::size_t>(left)]
                    [static_cast<std::size_t>(right)];
        }
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
        if (argc != 4) {
            throw std::runtime_error(
                "usage: search_su2_defect_cone "
                "(--step-one|--step-one-disjoint|--tls|--cross|--local"
                "|--edge-global|--edge-step) "
                "MAXIMUM_LABEL MAXIMUM_FACTORS"
            );
        }
        const std::string mode = argv[1];
        if (mode != "--step-one" && mode != "--step-one-disjoint"
            && mode != "--tls" && mode != "--cross") {
            if (mode != "--local" && mode != "--edge-global"
                && mode != "--edge-step") {
            throw std::runtime_error("unknown search mode");
            }
        }
        const bool disjoint_only = mode != "--step-one";
        const bool include_tls_edge = mode == "--tls";
        const bool cross_difference = mode == "--cross";
        const bool local_domination = mode == "--local";
        const bool edge_global = mode == "--edge-global";
        const bool edge_step = mode == "--edge-step";
        const int diagonal_step = include_tls_edge ? 2 : 1;
        const int maximum_label = std::stoi(argv[2]);
        const int maximum_factors = std::stoi(argv[3]);
        if (maximum_label < 1 || maximum_factors < 0) {
            throw std::runtime_error("invalid search bound");
        }
        const auto words = multisets(maximum_label, maximum_factors);
        std::size_t first_failure = words.size();
        int failure_left = -1;
        int failure_right = -1;
        cpp_int failure_value = 0;

#pragma omp parallel
        {
            std::size_t local_first = words.size();
            int local_left = -1;
            int local_right = -1;
            cpp_int local_value = 0;

#pragma omp for schedule(dynamic)
            for (std::size_t index = 0; index < words.size(); ++index) {
                if (edge_global) {
                    if (std::binary_search(
                            words[index].begin(), words[index].end(), 1
                        )) {
                        continue;
                    }
                    int total = 0;
                    for (int label : words[index]) {
                        total += label;
                    }
                    const int dimension = 2 * total + 5;
                    const Matrix coefficients = plus_coefficients(
                        words[index], dimension
                    );
                    const auto multiplicities = tensor_multiplicities(
                        words[index], dimension
                    );
                    bool failed = false;
                    for (int left = 2; left + 1 < dimension; ++left) {
                        if (std::binary_search(
                                words[index].begin(), words[index].end(), left
                            )) {
                            continue;
                        }
                        const cpp_int supply = multiplicities
                                [static_cast<std::size_t>(left - 1)]
                            + multiplicities[static_cast<std::size_t>(left + 1)];
                        const cpp_int demand = coefficients
                            [static_cast<std::size_t>(left)][1];
                        if (supply < demand) {
                            if (index < local_first) {
                                local_first = index;
                                local_left = left;
                                local_right = 1;
                                local_value = supply - demand;
                            }
                            failed = true;
                            break;
                        }
                    }
                    if (failed) {
                        continue;
                    }
                    continue;
                }
                const Matrix matrix = defect(
                    words[index], 4 * (maximum_label + maximum_factors
                        * maximum_label) + 10
                );
                bool failed = false;
                if (edge_step) {
                    if (std::binary_search(
                            words[index].begin(), words[index].end(), 1
                        )) {
                        continue;
                    }
                    const int first_p = words[index].empty()
                        ? 2 : std::max(2, words[index].back());
                    for (int p = first_p; p <= maximum_label && !failed; ++p) {
                        for (int left = 2;
                             left + p + 2 < static_cast<int>(matrix.size());
                             ++left) {
                            if (left == p || std::binary_search(
                                    words[index].begin(), words[index].end(), left
                                )) {
                                continue;
                            }
                            const cpp_int difference =
                                matrix[static_cast<std::size_t>(left)]
                                    [static_cast<std::size_t>(p - 1)]
                                + matrix[static_cast<std::size_t>(left)]
                                    [static_cast<std::size_t>(p + 1)]
                                - matrix[static_cast<std::size_t>(left - 1)]
                                    [static_cast<std::size_t>(p)]
                                - matrix[static_cast<std::size_t>(left + 1)]
                                    [static_cast<std::size_t>(p)];
                            if (difference < 0) {
                                if (index < local_first) {
                                    local_first = index;
                                    local_left = left;
                                    local_right = p;
                                    local_value = difference;
                                }
                                failed = true;
                                break;
                            }
                        }
                    }
                    continue;
                }
                if (local_domination) {
                    const int first_p = words[index].empty()
                        ? 1 : words[index].back();
                    for (int p = first_p; p <= maximum_label && !failed; ++p) {
                        for (int left = 2;
                             left + p < static_cast<int>(matrix.size()) && !failed;
                             ++left) {
                            if (std::binary_search(
                                    words[index].begin(), words[index].end(), left
                                ) || left == p) {
                                continue;
                            }
                            for (int right = 2; right <= left; ++right) {
                                if (left + right >= static_cast<int>(matrix.size())
                                    || std::binary_search(
                                        words[index].begin(), words[index].end(), right
                                    ) || right == p) {
                                    continue;
                                }
                                cpp_int first_supply = 0;
                                for_each_output(left, p, [&](int output) {
                                    first_supply += matrix
                                        [static_cast<std::size_t>(output)]
                                        [static_cast<std::size_t>(right)];
                                });
                                for_each_output(left - 2, p, [&](int output) {
                                    first_supply -= matrix
                                        [static_cast<std::size_t>(output)]
                                        [static_cast<std::size_t>(right - 2)];
                                });
                                cpp_int second_supply = 0;
                                for_each_output(right, p, [&](int output) {
                                    second_supply += matrix
                                        [static_cast<std::size_t>(left)]
                                        [static_cast<std::size_t>(output)];
                                });
                                for_each_output(right - 2, p, [&](int output) {
                                    second_supply -= matrix
                                        [static_cast<std::size_t>(left - 2)]
                                        [static_cast<std::size_t>(output)];
                                });
                                const int first_demand_label = left + right;
                                const int second_demand_label = left + right - 2;
                                const cpp_int first_demand = matrix
                                    [static_cast<std::size_t>(first_demand_label)]
                                    [static_cast<std::size_t>(p)];
                                const cpp_int second_demand = matrix
                                    [static_cast<std::size_t>(second_demand_label)]
                                    [static_cast<std::size_t>(p)];
                                if (first_supply < first_demand
                                    || second_supply < second_demand) {
                                    if (index < local_first) {
                                        local_first = index;
                                        local_left = left;
                                        local_right = right;
                                        local_value = std::min(
                                            first_supply - first_demand,
                                            second_supply - second_demand
                                        );
                                    }
                                    failed = true;
                                    break;
                                }
                            }
                        }
                    }
                    continue;
                }
                for (int left = 1;
                     left < static_cast<int>(matrix.size()) && !failed;
                     ++left) {
                    for (int right = 1;
                         right < static_cast<int>(matrix.size());
                         ++right) {
                        if (disjoint_only
                            && (std::binary_search(
                                    words[index].begin(), words[index].end(), left
                                )
                                || std::binary_search(
                                    words[index].begin(), words[index].end(), right
                                ))) {
                            continue;
                        }
                        cpp_int difference = 0;
                        if (cross_difference) {
                            if (left < 2 || right < 2 || left == right) {
                                continue;
                            }
                            difference = matrix[static_cast<std::size_t>(left)]
                                    [static_cast<std::size_t>(right)]
                                + matrix[static_cast<std::size_t>(left - 2)]
                                    [static_cast<std::size_t>(right - 2)]
                                - matrix[static_cast<std::size_t>(left - 2)]
                                    [static_cast<std::size_t>(right)]
                                - matrix[static_cast<std::size_t>(left)]
                                    [static_cast<std::size_t>(right - 2)];
                        } else if (include_tls_edge
                                   && (left == 1 || right == 1)) {
                            difference = matrix[static_cast<std::size_t>(left)]
                                [static_cast<std::size_t>(right)];
                        } else {
                            if (left < diagonal_step || right < diagonal_step) {
                                continue;
                            }
                            difference = matrix[static_cast<std::size_t>(left)]
                                    [static_cast<std::size_t>(right)]
                                - matrix[static_cast<std::size_t>(
                                    left - diagonal_step
                                )][static_cast<std::size_t>(
                                    right - diagonal_step
                                )];
                        }
                        const bool violates = cross_difference
                            ? difference > 0
                            : difference < 0;
                        if (violates) {
                            if (index < local_first) {
                                local_first = index;
                                local_left = left;
                                local_right = right;
                                local_value = difference;
                            }
                            failed = true;
                            break;
                        }
                    }
                }
                if (failed) {
                    continue;
                }
            }

#pragma omp critical
            {
                if (local_first < first_failure) {
                    first_failure = local_first;
                    failure_left = local_left;
                    failure_right = local_right;
                    failure_value = local_value;
                }
            }
        }

        std::cout << "SU2_DEFECT_CONE mode=" << mode
                  << " tested=" << words.size()
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors << '\n';
        if (first_failure != words.size()) {
            std::cout << "FAIL plus=";
            print_labels(words[first_failure]);
            std::cout << " at=(" << failure_left << ',' << failure_right
                      << ") difference=" << failure_value << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
