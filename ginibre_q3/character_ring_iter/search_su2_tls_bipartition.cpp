#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

using Multiplicities = std::vector<cpp_int>;

void for_each_output(int left, int right, const auto& function) {
    for (int output = std::abs(left - right); output <= left + right; output += 2) {
        function(output);
    }
}

Multiplicities tensor_multiplicities(const std::vector<int>& labels) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    Multiplicities current(static_cast<std::size_t>(total + 1));
    current[0] = 1;
    int current_total = 0;
    for (int label : labels) {
        Multiplicities next(static_cast<std::size_t>(total + 1));
        for (int input = 0; input <= current_total; ++input) {
            const cpp_int& coefficient = current[static_cast<std::size_t>(input)];
            if (coefficient == 0) {
                continue;
            }
            for_each_output(input, label, [&](int output) {
                next[static_cast<std::size_t>(output)] += coefficient;
            });
        }
        current = std::move(next);
        current_total += label;
    }
    return current;
}

cpp_int at(const Multiplicities& values, int label) {
    if (label < 0 || label >= static_cast<int>(values.size())) {
        return 0;
    }
    return values[static_cast<std::size_t>(label)];
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

void print_oriented_profile(
    const std::vector<int>& word,
    int left,
    int right
) {
    std::map<std::vector<int>, std::pair<std::uint64_t, cpp_int>> grouped;
    const std::uint64_t subset_count = std::uint64_t{1} << word.size();
    for (std::uint64_t mask = 0; mask < subset_count; ++mask) {
        std::vector<int> first;
        std::vector<int> second;
        for (std::size_t index = 0; index < word.size(); ++index) {
            if ((mask & (std::uint64_t{1} << index)) != 0U) {
                first.push_back(word[index]);
            } else {
                second.push_back(word[index]);
            }
        }
        const Multiplicities x = tensor_multiplicities(first);
        const Multiplicities y = tensor_multiplicities(second);
        const int sum = left + right;
        const cpp_int value = (at(x, sum) + at(x, sum - 2)) * at(y, 0)
            + at(x, left - 2) * at(y, right - 2)
            - at(x, left) * at(y, right);
        auto& entry = grouped[first];
        ++entry.first;
        entry.second += value;
    }
    cpp_int total = 0;
    std::cout << "oriented profile:\n";
    for (const auto& [labels, entry] : grouped) {
        if (entry.second == 0) {
            continue;
        }
        std::cout << "  first=";
        print_labels(labels);
        std::cout << " copies=" << entry.first << " total=" << entry.second << '\n';
        total += entry.second;
    }
    std::cout << "  aggregate=" << total << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: search_su2_tls_bipartition MAXIMUM_LABEL MAXIMUM_FACTORS"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 1 || maximum_factors < 0 || maximum_factors > 62) {
            throw std::runtime_error("invalid search bound");
        }

        const auto words = multisets(maximum_label, maximum_factors);
        std::uint64_t tested = 0;
        for (const auto& word : words) {
            int total = 0;
            for (int label : word) {
                total += label;
            }
            const std::uint64_t subset_count = std::uint64_t{1} << word.size();
            // Test each unordered bipartition only once.
            for (std::uint64_t mask = 0; mask < (subset_count + 1U) / 2U; ++mask) {
                std::vector<int> first;
                std::vector<int> second;
                for (std::size_t index = 0; index < word.size(); ++index) {
                    if ((mask & (std::uint64_t{1} << index)) != 0U) {
                        first.push_back(word[index]);
                    } else {
                        second.push_back(word[index]);
                    }
                }
                const Multiplicities first_multiplicity = tensor_multiplicities(first);
                const Multiplicities second_multiplicity = tensor_multiplicities(second);
                for (int left = 2; left <= total + 2; ++left) {
                    if (std::binary_search(word.begin(), word.end(), left)) {
                        continue;
                    }
                    for (int right = 2; right <= left; ++right) {
                        if (std::binary_search(word.begin(), word.end(), right)) {
                            continue;
                        }
                        const int sum = left + right;
                        const auto oriented = [&](const Multiplicities& x,
                                                  const Multiplicities& y) {
                            return (at(x, sum) + at(x, sum - 2)) * at(y, 0)
                                + at(x, left - 2) * at(y, right - 2)
                                - at(x, left) * at(y, right);
                        };
                        const cpp_int value = oriented(
                            first_multiplicity, second_multiplicity
                        ) + oriented(second_multiplicity, first_multiplicity);
                        ++tested;
                        if (value < 0) {
                            std::cout << "SU2_TLS_BIPARTITION FAIL plus=";
                            print_labels(word);
                            std::cout << " first=";
                            print_labels(first);
                            std::cout << " second=";
                            print_labels(second);
                            std::cout << " target=(" << left << ',' << right
                                      << ") contribution=" << value << '\n';
                            print_oriented_profile(word, left, right);
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }
        std::cout << "SU2_TLS_BIPARTITION PASS tested=" << tested
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
