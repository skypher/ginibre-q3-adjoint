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

using Vector = std::vector<cpp_int>;

void for_each_output(int left, int right, const auto& function) {
    for (int output = std::abs(left - right); output <= left + right;
         output += 2) {
        function(output);
    }
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

Vector oriented_values(
    const std::vector<int>& word,
    unsigned long long mask,
    int dimension
) {
    const int p = word.front();
    const Vector k_pm1 = kernel_column(word, mask, true, p - 1, dimension);
    const Vector k_1 = kernel_column(word, mask, true, 1, dimension);
    const Vector k_0 = kernel_column(word, mask, true, 0, dimension);
    const Vector k_p = kernel_column(word, mask, true, p, dimension);
    const Vector ell_0 = kernel_column(word, mask, false, 0, dimension);
    Vector answer(static_cast<std::size_t>(dimension));
    for (int target = 1; target < dimension; ++target) {
        answer[static_cast<std::size_t>(target)] =
            k_pm1[static_cast<std::size_t>(target)] * ell_0[0]
            + k_1[static_cast<std::size_t>(target)]
                * ell_0[static_cast<std::size_t>(p)]
            - k_0[static_cast<std::size_t>(target)]
                * ell_0[static_cast<std::size_t>(p - 1)]
            - k_p[static_cast<std::size_t>(target)] * ell_0[1];
    }
    return answer;
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
        if (argc < 3 || argc > 4) {
            throw std::runtime_error(
                "usage: search_su2_edge_packet_toggle MAXIMUM_LABEL "
                "MAXIMUM_FACTORS [--even-summed|--p2|--p2-even-summed]"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 2 || maximum_factors < 2
            || maximum_factors > 63) {
            throw std::runtime_error("invalid search bound");
        }
        const std::string mode = argc == 4 ? argv[3] : "--local";
        const bool even_summed = mode == "--even-summed"
            || mode == "--p2-even-summed";
        const bool restrict_p2 = mode == "--p2"
            || mode == "--p2-even-summed";
        if (mode != "--local" && mode != "--even-summed"
            && mode != "--p2" && mode != "--p2-even-summed") {
            throw std::runtime_error("unknown mode");
        }
        const auto words = multisets(maximum_label, maximum_factors);
        std::size_t first_failure = words.size();
        unsigned long long failure_mask = 0;
        int failure_target = -1;
        int failure_position = -1;
        cpp_int failure_value = 0;

#pragma omp parallel
        {
            std::size_t local_failure = words.size();
            unsigned long long local_mask = 0;
            int local_target = -1;
            int local_position = -1;
            cpp_int local_value = 0;

#pragma omp for schedule(dynamic)
            for (std::size_t index = 0; index < words.size(); ++index) {
                const auto& word = words[index];
                if (word.size() < 2U || (restrict_p2 && word.front() != 2)) {
                    continue;
                }
                std::size_t odd_position = word.size();
                for (std::size_t position = 1; position < word.size(); ++position) {
                    if ((word[position] & 1) != 0) {
                        odd_position = position;
                        break;
                    }
                }
                if (odd_position == word.size()) {
                    continue;
                }
                int total = 0;
                for (int label : word) {
                    total += label;
                }
                const int dimension = 2 * total + 5;
                const std::size_t suffix_size = word.size() - 1U;
                const unsigned long long mask_end = 1ULL << suffix_size;
                const unsigned long long toggle =
                    1ULL << (odd_position - 1U);
                bool failed = false;
                if (even_summed) {
                    std::vector<std::size_t> other_odd_positions;
                    for (std::size_t position = 1; position < word.size();
                         ++position) {
                        if (position != odd_position
                            && (word[position] & 1) != 0) {
                            other_odd_positions.push_back(position);
                        }
                    }
                    const std::size_t group_count = std::size_t{1}
                        << other_odd_positions.size();
                    std::vector<Vector> groups(
                        group_count,
                        Vector(static_cast<std::size_t>(dimension))
                    );
                    for (unsigned long long mask = 0; mask < mask_end; ++mask) {
                        std::size_t group = 0U;
                        for (std::size_t bit = 0;
                             bit < other_odd_positions.size(); ++bit) {
                            const std::size_t suffix_bit =
                                other_odd_positions[bit] - 1U;
                            if (((mask >> suffix_bit) & 1ULL) != 0ULL) {
                                group |= std::size_t{1} << bit;
                            }
                        }
                        const Vector values = oriented_values(
                            word, mask, dimension
                        );
                        for (int target = 1; target < dimension; ++target) {
                            groups[group][static_cast<std::size_t>(target)] +=
                                values[static_cast<std::size_t>(target)];
                        }
                    }
                    for (std::size_t group = 0; group < group_count && !failed;
                         ++group) {
                        for (int target = 1; target < dimension; ++target) {
                            const cpp_int& value = groups[group]
                                [static_cast<std::size_t>(target)];
                            if (value < 0) {
                                local_failure = index;
                                local_mask = static_cast<unsigned long long>(group);
                                local_target = target;
                                local_position = static_cast<int>(odd_position);
                                local_value = value;
                                failed = true;
                                break;
                            }
                        }
                    }
                } else {
                    for (unsigned long long mask = 0;
                         mask < mask_end && !failed; ++mask) {
                        if ((mask & toggle) != 0ULL) {
                            continue;
                        }
                        const Vector first = oriented_values(
                            word, mask, dimension
                        );
                        const Vector second = oriented_values(
                            word, mask | toggle, dimension
                        );
                        for (int target = 1; target < dimension; ++target) {
                            const cpp_int value = first
                                    [static_cast<std::size_t>(target)]
                                + second[static_cast<std::size_t>(target)];
                            if (value < 0) {
                                local_failure = index;
                                local_mask = mask;
                                local_target = target;
                                local_position = static_cast<int>(odd_position);
                                local_value = value;
                                failed = true;
                                break;
                            }
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
                    failure_position = local_position;
                    failure_value = local_value;
                }
            }
        }

        std::cout << "SU2_EDGE_PACKET_TOGGLE tested_words=" << words.size()
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors
                  << " mode=" << mode
                  << '\n';
        if (first_failure != words.size()) {
            std::cout << "FAIL plus=";
            print_labels(words[first_failure]);
            std::cout << " odd_suffix_position=" << failure_position
                      << " base_mask=" << failure_mask
                      << " target=" << failure_target
                      << " toggle_slack=" << failure_value << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
