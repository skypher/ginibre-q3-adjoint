#include <algorithm>
#include <cstdint>
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

Matrix zero_matrix(int dimension) {
    return Matrix(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(static_cast<std::size_t>(dimension))
    );
}

const cpp_int& matrix_entry(const Matrix& matrix, int row, int column) {
    static const cpp_int zero = 0;
    if (row < 0 || column < 0
        || row >= static_cast<int>(matrix.size())
        || column >= static_cast<int>(matrix.size())) {
        return zero;
    }
    return matrix[static_cast<std::size_t>(row)]
        [static_cast<std::size_t>(column)];
}

Matrix add_plus(const Matrix& current, int label) {
    const int dimension = static_cast<int>(current.size());
    Matrix next = zero_matrix(dimension);
    for (int left = 0; left < dimension; ++left) {
        for (int right = 0; right < dimension; ++right) {
            const cpp_int& value = current[static_cast<std::size_t>(left)]
                [static_cast<std::size_t>(right)];
            if (value == 0) {
                continue;
            }
            for_each_output(left, label, [&](int output) {
                if (output < dimension) {
                    next[static_cast<std::size_t>(output)]
                        [static_cast<std::size_t>(right)] += value;
                }
            });
            for_each_output(right, label, [&](int output) {
                if (output < dimension) {
                    next[static_cast<std::size_t>(left)]
                        [static_cast<std::size_t>(output)] += value;
                }
            });
        }
    }
    return next;
}

Matrix forced_kernel_range(
    const std::vector<int>& labels,
    std::size_t first,
    std::size_t last,
    int dimension
) {
    Matrix current = zero_matrix(dimension);
    for (int state = 0; state < dimension; ++state) {
        current[static_cast<std::size_t>(state)]
            [static_cast<std::size_t>(state)] = 1;
    }
    for (std::size_t index = first; index < last; ++index) {
        Matrix next = zero_matrix(dimension);
        for (int start = 0; start < dimension; ++start) {
            for (int middle = 0; middle < dimension; ++middle) {
                const cpp_int& value = current[static_cast<std::size_t>(start)]
                    [static_cast<std::size_t>(middle)];
                if (value == 0) {
                    continue;
                }
                for_each_output(middle, labels[index], [&](int output) {
                    if (output < dimension) {
                        next[static_cast<std::size_t>(start)]
                            [static_cast<std::size_t>(output)] += value;
                    }
                });
            }
        }
        current = std::move(next);
    }
    return current;
}

bool transition_exists(int input, int label, int output) {
    return output >= std::abs(input - label) && output <= input + label
        && ((input + label - output) & 1) == 0;
}

cpp_int saturated_suffix(
    const Matrix& middle_kernel,
    int start,
    int output,
    int last_label
) {
    cpp_int answer = matrix_entry(
        middle_kernel, start, output + last_label
    );
    if (last_label > output) {
        answer += matrix_entry(
            middle_kernel, start, last_label - output
        );
    }
    return answer;
}

cpp_int oriented_residual(
    const Matrix& prefix,
    const Matrix& suffix,
    const Matrix& middle_suffix,
    int split_label,
    int last_label,
    int first_target,
    int second_target
) {
    const int dimension = static_cast<int>(prefix.size());
    cpp_int answer = 0;
    for (int first = 0; first < dimension; ++first) {
        for (int second = 0; second < dimension; ++second) {
            const cpp_int& multiplicity = prefix[static_cast<std::size_t>(first)]
                [static_cast<std::size_t>(second)];
            if (multiplicity == 0
                || !transition_exists(second, split_label, second_target)) {
                continue;
            }
            const bool second_saturated =
                second_target == std::abs(second - split_label);
            const cpp_int continuation = second_saturated
                ? suffix[static_cast<std::size_t>(first)]
                    [static_cast<std::size_t>(first_target)]
                : saturated_suffix(
                    middle_suffix, first, first_target, last_label
                );
            answer += multiplicity * continuation;
        }
    }
    return answer;
}

cpp_int boundary_first_orientation(
    const Matrix& prefix,
    const Matrix& suffix,
    int split_label,
    int target
) {
    const int dimension = static_cast<int>(prefix.size());
    cpp_int answer = 0;
    for (int first = 0; first < dimension; ++first) {
        answer += prefix[static_cast<std::size_t>(first)]
            [static_cast<std::size_t>(split_label)]
            * matrix_entry(suffix, first, target);
    }
    return answer;
}

cpp_int boundary_second_orientation(
    const Matrix& prefix,
    const Matrix& suffix,
    int split_label,
    int target
) {
    const int dimension = static_cast<int>(prefix.size());
    cpp_int answer = 0;
    for (int first = 0; first < dimension; ++first) {
        if (!transition_exists(first, split_label, target)) {
            continue;
        }
        for (int second = 0; second < dimension; ++second) {
            answer += prefix[static_cast<std::size_t>(first)]
                [static_cast<std::size_t>(second)]
                * suffix[static_cast<std::size_t>(second)][0];
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
        if (argc != 3 && argc != 4) {
            throw std::runtime_error(
                "usage: search_su2_tls_last_update MAXIMUM_LABEL "
                "MAXIMUM_FACTORS [--prefix|--prefix-edge|--exterior]"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 1 || maximum_factors < 2) {
            throw std::runtime_error("invalid search bound");
        }
        const bool prefix_mode = argc == 4;
        const bool exterior_mode = prefix_mode
            && std::string(argv[3]) == "--exterior";
        const bool edge_only = prefix_mode
            && std::string(argv[3]) == "--prefix-edge";
        if (prefix_mode && std::string(argv[3]) != "--prefix"
            && !edge_only && !exterior_mode) {
            throw std::runtime_error("unknown mode");
        }
        const auto words = multisets(maximum_label, maximum_factors);
        std::size_t first_failure = words.size();
        std::size_t failure_split = 0;
        int failure_left = -1;
        int failure_right = -1;
        cpp_int failure_residual = 0;
        cpp_int failure_boundary = 0;
        bool minimum_set = false;
        bool active_minimum_set = false;
        std::size_t minimum_word = words.size();
        std::size_t minimum_split = 0;
        int minimum_left = -1;
        int minimum_right = -1;
        cpp_int minimum_residual = 0;
        cpp_int minimum_boundary = 0;
        std::size_t active_minimum_word = words.size();
        std::size_t active_minimum_split = 0;
        int active_minimum_left = -1;
        int active_minimum_right = -1;
        cpp_int active_minimum_residual = 0;
        cpp_int active_minimum_boundary = 0;

#pragma omp parallel
        {
            std::size_t local_first = words.size();
            std::size_t local_split = 0;
            int local_left = -1;
            int local_right = -1;
            cpp_int local_residual = 0;
            cpp_int local_boundary = 0;
            bool local_minimum_set = false;
            bool local_active_minimum_set = false;
            std::size_t local_minimum_word = words.size();
            std::size_t local_minimum_split = 0;
            int local_minimum_left = -1;
            int local_minimum_right = -1;
            cpp_int local_minimum_residual = 0;
            cpp_int local_minimum_boundary = 0;
            std::size_t local_active_minimum_word = words.size();
            std::size_t local_active_minimum_split = 0;
            int local_active_minimum_left = -1;
            int local_active_minimum_right = -1;
            cpp_int local_active_minimum_residual = 0;
            cpp_int local_active_minimum_boundary = 0;

#pragma omp for schedule(dynamic)
            for (std::size_t word_index = 0;
                 word_index < words.size();
                 ++word_index) {
                const auto& word = words[word_index];
                if (word.size() < 2U) {
                    continue;
                }
                int total = 0;
                for (int label : word) {
                    total += label;
                }
                const int dimension = total + 5;
                const Matrix whole_word = forced_kernel_range(
                    word, 0U, word.size(), dimension
                );
                Matrix prefix = zero_matrix(dimension);
                prefix[0][0] = 1;
                Matrix cumulative_residual = zero_matrix(dimension);
                Matrix cumulative_boundary = zero_matrix(dimension);
                bool failed = false;
                for (std::size_t split = 0;
                     split + 1U < word.size() && !failed;
                     ++split) {
                    const Matrix suffix = forced_kernel_range(
                        word, split + 1U, word.size(), dimension
                    );
                    const Matrix middle = forced_kernel_range(
                        word, split + 1U, word.size() - 1U, dimension
                    );
                    for (int left = 2; left <= total + 2 && !failed; ++left) {
                        if (std::binary_search(word.begin(), word.end(), left)) {
                            continue;
                        }
                        const int right_end = edge_only ? 1 : left;
                        for (int right = 1; right <= right_end; ++right) {
                            if (std::binary_search(word.begin(), word.end(), right)) {
                                continue;
                            }
                            const cpp_int residual = oriented_residual(
                                prefix, suffix, middle, word[split], word.back(),
                                left, right
                            ) + oriented_residual(
                                prefix, suffix, middle, word[split], word.back(),
                                right, left
                            );
                            const int sum = left + right;
                            const cpp_int boundary =
                                boundary_first_orientation(
                                    prefix, suffix, word[split], sum
                                )
                                + boundary_first_orientation(
                                    prefix, suffix, word[split], sum - 2
                                )
                                + boundary_second_orientation(
                                    prefix, suffix, word[split], sum
                                )
                                + boundary_second_orientation(
                                    prefix, suffix, word[split], sum - 2
                                );
                            cpp_int& residual_sum = cumulative_residual
                                [static_cast<std::size_t>(left)]
                                [static_cast<std::size_t>(right)];
                            cpp_int& boundary_sum = cumulative_boundary
                                [static_cast<std::size_t>(left)]
                                [static_cast<std::size_t>(right)];
                            if (split == 0U) {
                                boundary_sum += matrix_entry(
                                    whole_word, 0, sum
                                ) + matrix_entry(whole_word, 0, sum - 2);
                            }
                            residual_sum += residual;
                            if (!exterior_mode) {
                                boundary_sum += boundary;
                            }
                            if (prefix_mode) {
                                const cpp_int slack = boundary_sum - residual_sum;
                                const cpp_int local_slack = local_minimum_boundary
                                    - local_minimum_residual;
                                if (!local_minimum_set || slack < local_slack) {
                                    local_minimum_set = true;
                                    local_minimum_word = word_index;
                                    local_minimum_split = split;
                                    local_minimum_left = left;
                                    local_minimum_right = right;
                                    local_minimum_residual = residual_sum;
                                    local_minimum_boundary = boundary_sum;
                                }
                                const cpp_int local_active_slack =
                                    local_active_minimum_boundary
                                    - local_active_minimum_residual;
                                if (residual_sum > 0
                                    && (!local_active_minimum_set
                                        || slack < local_active_slack)) {
                                    local_active_minimum_set = true;
                                    local_active_minimum_word = word_index;
                                    local_active_minimum_split = split;
                                    local_active_minimum_left = left;
                                    local_active_minimum_right = right;
                                    local_active_minimum_residual = residual_sum;
                                    local_active_minimum_boundary = boundary_sum;
                                }
                            }
                            const bool inequality_fails = prefix_mode
                                ? residual_sum > boundary_sum
                                : residual > boundary;
                            if (inequality_fails) {
                                if (word_index < local_first) {
                                    local_first = word_index;
                                    local_split = split;
                                    local_left = left;
                                    local_right = right;
                                    local_residual = prefix_mode
                                        ? residual_sum : residual;
                                    local_boundary = prefix_mode
                                        ? boundary_sum : boundary;
                                }
                                failed = true;
                                break;
                            }
                        }
                    }
                    prefix = add_plus(prefix, word[split]);
                }
            }

#pragma omp critical
            {
                if (local_first < first_failure) {
                    first_failure = local_first;
                    failure_split = local_split;
                    failure_left = local_left;
                    failure_right = local_right;
                    failure_residual = local_residual;
                    failure_boundary = local_boundary;
                }
                const cpp_int minimum_slack = minimum_boundary - minimum_residual;
                const cpp_int local_minimum_slack = local_minimum_boundary
                    - local_minimum_residual;
                if (local_minimum_set
                    && (!minimum_set || local_minimum_slack < minimum_slack)) {
                    minimum_set = true;
                    minimum_word = local_minimum_word;
                    minimum_split = local_minimum_split;
                    minimum_left = local_minimum_left;
                    minimum_right = local_minimum_right;
                    minimum_residual = local_minimum_residual;
                    minimum_boundary = local_minimum_boundary;
                }
                const cpp_int active_minimum_slack = active_minimum_boundary
                    - active_minimum_residual;
                const cpp_int local_active_minimum_slack =
                    local_active_minimum_boundary - local_active_minimum_residual;
                if (local_active_minimum_set
                    && (!active_minimum_set
                        || local_active_minimum_slack < active_minimum_slack)) {
                    active_minimum_set = true;
                    active_minimum_word = local_active_minimum_word;
                    active_minimum_split = local_active_minimum_split;
                    active_minimum_left = local_active_minimum_left;
                    active_minimum_right = local_active_minimum_right;
                    active_minimum_residual = local_active_minimum_residual;
                    active_minimum_boundary = local_active_minimum_boundary;
                }
            }
        }

        std::cout << "SU2_TLS_LAST_UPDATE tested_words=" << words.size()
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors
                  << " mode="
                  << (exterior_mode
                          ? "exterior"
                          : (edge_only
                                ? "prefix-edge"
                                : (prefix_mode ? "prefix" : "pointwise")))
                  << '\n';
        if (first_failure != words.size()) {
            std::cout << "FAIL plus=";
            print_labels(words[first_failure]);
            std::cout << " split_index=" << failure_split
                      << " target=(" << failure_left << ',' << failure_right
                      << ") residual=" << failure_residual
                      << " boundary=" << failure_boundary << '\n';

            const auto& word = words[first_failure];
            int total = 0;
            for (int label : word) {
                total += label;
            }
            const int dimension = total + 5;
            const Matrix whole_word = forced_kernel_range(
                word, 0U, word.size(), dimension
            );
            Matrix prefix = zero_matrix(dimension);
            prefix[0][0] = 1;
            cpp_int cumulative_residual = 0;
            const int failure_sum = failure_left + failure_right;
            const cpp_int exterior_boundary = matrix_entry(
                whole_word, 0, failure_sum
            ) + matrix_entry(whole_word, 0, failure_sum - 2);
            cpp_int cumulative_boundary = prefix_mode ? exterior_boundary : 0;
            if (prefix_mode) {
                std::cout << "  never_updated_boundary="
                          << exterior_boundary << '\n';
            }
            for (std::size_t split = 0; split + 1U < word.size(); ++split) {
                const Matrix suffix = forced_kernel_range(
                    word, split + 1U, word.size(), dimension
                );
                const Matrix middle = forced_kernel_range(
                    word, split + 1U, word.size() - 1U, dimension
                );
                const cpp_int residual = oriented_residual(
                    prefix, suffix, middle, word[split], word.back(),
                    failure_left, failure_right
                ) + oriented_residual(
                    prefix, suffix, middle, word[split], word.back(),
                    failure_right, failure_left
                );
                const int sum = failure_left + failure_right;
                const cpp_int boundary = boundary_first_orientation(
                    prefix, suffix, word[split], sum
                ) + boundary_first_orientation(
                    prefix, suffix, word[split], sum - 2
                ) + boundary_second_orientation(
                    prefix, suffix, word[split], sum
                ) + boundary_second_orientation(
                    prefix, suffix, word[split], sum - 2
                );
                cumulative_residual += residual;
                if (!exterior_mode) {
                    cumulative_boundary += boundary;
                }
                std::cout << "  split_index=" << split
                          << " residual=" << residual
                          << " boundary=" << boundary
                          << " prefix_residual=" << cumulative_residual
                          << " prefix_boundary=" << cumulative_boundary << '\n';
                prefix = add_plus(prefix, word[split]);
            }
            return EXIT_FAILURE;
        }
        if (prefix_mode && minimum_set) {
            std::cout << "minimum_prefix_slack="
                      << minimum_boundary - minimum_residual << " plus=";
            print_labels(words[minimum_word]);
            std::cout << " split_index=" << minimum_split
                      << " target=(" << minimum_left << ',' << minimum_right
                      << ") residual=" << minimum_residual
                      << " boundary=" << minimum_boundary << '\n';
        }
        if (prefix_mode && active_minimum_set) {
            std::cout << "active_minimum_prefix_slack="
                      << active_minimum_boundary - active_minimum_residual
                      << " plus=";
            print_labels(words[active_minimum_word]);
            std::cout << " split_index=" << active_minimum_split
                      << " target=(" << active_minimum_left << ','
                      << active_minimum_right << ") residual="
                      << active_minimum_residual
                      << " boundary=" << active_minimum_boundary << '\n';
        }
        std::cout << "PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
