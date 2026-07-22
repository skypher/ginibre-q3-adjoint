#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

using boost::multiprecision::cpp_int;

namespace {

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
    visit(visit, 2, maximum_factors);
    return result;
}

Vector apply_label(const Vector& current, int label) {
    const int dimension = static_cast<int>(current.size());
    Vector next(static_cast<std::size_t>(dimension));
    for (int input = 0; input < dimension; ++input) {
        const cpp_int& value = current[static_cast<std::size_t>(input)];
        if (value == 0) {
            continue;
        }
        for_each_output(input, label, [&](int output) {
            if (output < dimension) {
                next[static_cast<std::size_t>(output)] += value;
            }
        });
    }
    return next;
}

Vector kernel_column(
    const std::vector<int>& word,
    unsigned long long mask,
    bool selected,
    int start,
    int dimension
) {
    Vector current(static_cast<std::size_t>(dimension));
    current[static_cast<std::size_t>(start)] = 1;
    for (std::size_t position = 1; position < word.size(); ++position) {
        const bool in_mask = ((mask >> (position - 1U)) & 1ULL) != 0ULL;
        if (in_mask == selected) {
            current = apply_label(current, word[position]);
        }
    }
    return current;
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
                "usage: search_su2_edge_packet_split MAXIMUM_LABEL "
                "MAXIMUM_FACTORS"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 2 || maximum_factors < 1 || maximum_factors > 63) {
            throw std::runtime_error("invalid search bound");
        }
        const auto words = multisets(maximum_label, maximum_factors);
        std::size_t first_failure = words.size();
        unsigned long long failure_mask = 0;
        int failure_target = -1;
        cpp_int failure_value = 0;

#pragma omp parallel
        {
            std::size_t local_failure = words.size();
            unsigned long long local_mask = 0;
            int local_target = -1;
            cpp_int local_value = 0;

#pragma omp for schedule(dynamic)
            for (std::size_t index = 0; index < words.size(); ++index) {
                const auto& word = words[index];
                if (word.empty()) {
                    continue;
                }
                int total = 0;
                for (int label : word) {
                    total += label;
                }
                const int dimension = 2 * total + 5;
                const int p = word.front();
                const std::size_t suffix_size = word.size() - 1U;
                const unsigned long long mask_end = 1ULL << suffix_size;
                bool failed = false;
                for (unsigned long long mask = 0; mask < mask_end && !failed; ++mask) {
                    const Vector k_pm1 = kernel_column(
                        word, mask, true, p - 1, dimension
                    );
                    const Vector k_1 = kernel_column(
                        word, mask, true, 1, dimension
                    );
                    const Vector k_0 = kernel_column(
                        word, mask, true, 0, dimension
                    );
                    const Vector k_p = kernel_column(
                        word, mask, true, p, dimension
                    );
                    const Vector ell_pm1 = kernel_column(
                        word, mask, false, p - 1, dimension
                    );
                    const Vector ell_1 = kernel_column(
                        word, mask, false, 1, dimension
                    );
                    const Vector ell_0 = kernel_column(
                        word, mask, false, 0, dimension
                    );
                    const Vector ell_p = kernel_column(
                        word, mask, false, p, dimension
                    );
                    for (int target = 2; target < dimension; ++target) {
                        const cpp_int value =
                            k_pm1[static_cast<std::size_t>(target)] * ell_0[0]
                            + k_1[static_cast<std::size_t>(target)]
                                * ell_0[static_cast<std::size_t>(p)]
                            - k_0[static_cast<std::size_t>(target)]
                                * ell_0[static_cast<std::size_t>(p - 1)]
                            - k_p[static_cast<std::size_t>(target)] * ell_0[1]
                            + ell_pm1[static_cast<std::size_t>(target)] * k_0[0]
                            + ell_1[static_cast<std::size_t>(target)]
                                * k_0[static_cast<std::size_t>(p)]
                            - ell_0[static_cast<std::size_t>(target)]
                                * k_0[static_cast<std::size_t>(p - 1)]
                            - ell_p[static_cast<std::size_t>(target)] * k_0[1];
                        if (value < 0) {
                            if (index < local_failure) {
                                local_failure = index;
                                local_mask = mask;
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
                    failure_mask = local_mask;
                    failure_target = local_target;
                    failure_value = local_value;
                }
            }
        }

        std::cout << "SU2_EDGE_PACKET_SPLIT tested_words=" << words.size()
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors << '\n';
        if (first_failure != words.size()) {
            std::cout << "FAIL plus=";
            print_labels(words[first_failure]);
            std::cout << " suffix_mask=" << failure_mask
                      << " target=" << failure_target
                      << " endpoint_slack=" << failure_value << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
