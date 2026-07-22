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
using Vector = std::vector<cpp_int>;

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

Matrix plus_coefficients(
    const std::vector<int>& labels,
    std::size_t last,
    int dimension
) {
    Matrix current(
        static_cast<std::size_t>(dimension),
        Vector(static_cast<std::size_t>(dimension))
    );
    current[0][0] = 1;
    for (std::size_t index = 0; index < last; ++index) {
        const int label = labels[index];
        Matrix next(
            static_cast<std::size_t>(dimension),
            Vector(static_cast<std::size_t>(dimension))
        );
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
        current = std::move(next);
    }
    return current;
}

Vector suffix_multiplicities(
    const std::vector<int>& labels,
    std::size_t first,
    int dimension
) {
    Vector current(static_cast<std::size_t>(dimension));
    current[0] = 1;
    for (std::size_t index = first; index < labels.size(); ++index) {
        Vector next(static_cast<std::size_t>(dimension));
        for (int input = 0; input < dimension; ++input) {
            const cpp_int& value = current[static_cast<std::size_t>(input)];
            if (value == 0) {
                continue;
            }
            for_each_output(input, labels[index], [&](int output) {
                if (output < dimension) {
                    next[static_cast<std::size_t>(output)] += value;
                }
            });
        }
        current = std::move(next);
    }
    return current;
}

cpp_int commutator_at(const Matrix& coefficients, const Vector& suffix, int target) {
    const int dimension = static_cast<int>(suffix.size());
    cpp_int answer = 0;
    for (int state = 0; state < dimension; ++state) {
        const cpp_int left_coefficient =
            (target > 0
                 ? coefficients[static_cast<std::size_t>(target - 1)]
                       [static_cast<std::size_t>(state)]
                 : cpp_int(0))
            + (target + 1 < dimension
                   ? coefficients[static_cast<std::size_t>(target + 1)]
                         [static_cast<std::size_t>(state)]
                   : cpp_int(0));
        cpp_int shifted_suffix = 0;
        if (state > 0) {
            shifted_suffix += suffix[static_cast<std::size_t>(state - 1)];
        }
        if (state + 1 < dimension) {
            shifted_suffix += suffix[static_cast<std::size_t>(state + 1)];
        }
        answer += left_coefficient * suffix[static_cast<std::size_t>(state)]
            - coefficients[static_cast<std::size_t>(target)]
                  [static_cast<std::size_t>(state)]
                * shifted_suffix;
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
        if (argc != 3) {
            throw std::runtime_error(
                "usage: search_su2_edge_commutator MAXIMUM_LABEL MAXIMUM_FACTORS"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 1 || maximum_factors < 2) {
            throw std::runtime_error("invalid search bound");
        }
        const auto words = multisets(maximum_label, maximum_factors);
        std::size_t first_failure = words.size();
        std::size_t failure_cut = 0;
        int failure_target = -1;
        cpp_int failure_value = 0;

#pragma omp parallel
        {
            std::size_t local_failure = words.size();
            std::size_t local_cut = 0;
            int local_target = -1;
            cpp_int local_value = 0;

#pragma omp for schedule(dynamic)
            for (std::size_t word_index = 0; word_index < words.size(); ++word_index) {
                const auto& word = words[word_index];
                if (word.size() < 2U || std::binary_search(word.begin(), word.end(), 1)) {
                    continue;
                }
                int total = 0;
                for (int label : word) {
                    total += label;
                }
                const int dimension = 2 * total + 5;
                bool failed = false;
                for (std::size_t cut = 1; cut < word.size() && !failed; ++cut) {
                    const Matrix coefficients = plus_coefficients(word, cut, dimension);
                    const Vector suffix = suffix_multiplicities(word, cut, dimension);
                    for (int target = 2; target + 1 < dimension; ++target) {
                        if (std::binary_search(word.begin(), word.end(), target)) {
                            continue;
                        }
                        const cpp_int value = commutator_at(
                            coefficients, suffix, target
                        );
                        if (value < 0) {
                            if (word_index < local_failure) {
                                local_failure = word_index;
                                local_cut = cut;
                                local_target = target;
                                local_value = value;
                            }
                            failed = true;
                            break;
                        }
                    }
                }
            }

#pragma omp critical
            {
                if (local_failure < first_failure) {
                    first_failure = local_failure;
                    failure_cut = local_cut;
                    failure_target = local_target;
                    failure_value = local_value;
                }
            }
        }

        std::cout << "SU2_EDGE_COMMUTATOR tested_words=" << words.size()
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors << '\n';
        if (first_failure != words.size()) {
            std::cout << "FAIL plus=";
            print_labels(words[first_failure]);
            std::cout << " cut=" << failure_cut
                      << " target=" << failure_target
                      << " commutator=" << failure_value << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
