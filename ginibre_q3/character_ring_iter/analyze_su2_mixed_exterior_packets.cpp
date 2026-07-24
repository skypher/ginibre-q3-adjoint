#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

using Matrix = std::vector<std::vector<cpp_int>>;

Matrix identity_matrix(int maximum_state) {
    Matrix result(
        static_cast<std::size_t>(maximum_state + 1),
        std::vector<cpp_int>(static_cast<std::size_t>(maximum_state + 1))
    );
    for (int state = 0; state <= maximum_state; ++state) {
        result[static_cast<std::size_t>(state)]
              [static_cast<std::size_t>(state)] = 1;
    }
    return result;
}

Matrix multiply_by_fusion_label(const Matrix& input, int label) {
    const int maximum_state = static_cast<int>(input.size()) - 1;
    Matrix result(
        input.size(), std::vector<cpp_int>(input.size())
    );
    for (int row = 0; row <= maximum_state; ++row) {
        for (int middle = 0; middle <= maximum_state; ++middle) {
            const cpp_int& coefficient =
                input[static_cast<std::size_t>(row)]
                     [static_cast<std::size_t>(middle)];
            if (coefficient == 0) {
                continue;
            }
            for (int column = std::abs(middle - label);
                 column <= middle + label && column <= maximum_state;
                 column += 2) {
                result[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(column)] += coefficient;
            }
        }
    }
    return result;
}

cpp_int mixed_minor(
    const Matrix& first,
    const Matrix& second,
    int target,
    int source
) {
    const auto index = [](int value) {
        return static_cast<std::size_t>(value);
    };
    return first[index(target)][index(source)] * second[0U][0U]
        + first[0U][0U] * second[index(target)][index(source)]
        - first[index(target)][0U] * second[index(source)][0U]
        - first[index(source)][0U] * second[index(target)][0U];
}

void print_word(const std::vector<int>& word) {
    std::cout << '[';
    for (std::size_t index = 0U; index < word.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << word[index];
    }
    std::cout << ']';
}

void print_subset(
    const std::vector<int>& word,
    std::size_t mask
) {
    std::cout << '[';
    bool first = true;
    for (std::size_t index = 0U; index < word.size(); ++index) {
        if ((mask & (std::size_t{1} << index)) == 0U) {
            continue;
        }
        if (!first) {
            std::cout << ',';
        }
        first = false;
        std::cout << word[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_mixed_exterior_packets "
                "MAXIMUM_LABEL MAXIMUM_FACTORS"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 1 || maximum_factors < 1
            || maximum_factors >= static_cast<int>(
                std::numeric_limits<std::size_t>::digits
            )) {
            throw std::runtime_error("invalid bound");
        }

        std::vector<int> word;
        std::size_t words = 0U;
        std::size_t packets = 0U;
        bool found = false;

        const auto visit = [&](const auto& self, int first_label) -> void {
            if (!word.empty() && !found) {
                ++words;
                int total = 0;
                for (int label : word) {
                    total += label;
                }
                const int maximum_state = total + maximum_label;
                const std::size_t subset_count =
                    std::size_t{1} << word.size();
                std::vector<Matrix> products(
                    subset_count, identity_matrix(maximum_state)
                );
                for (std::size_t mask = 1U; mask < subset_count; ++mask) {
                    const std::size_t bit = mask & (~mask + 1U);
                    const unsigned raw_index =
                        static_cast<unsigned>(__builtin_ctzll(bit));
                    const std::size_t index =
                        static_cast<std::size_t>(raw_index);
                    products[mask] = multiply_by_fusion_label(
                        products[mask ^ bit], word[index]
                    );
                }
                const std::size_t full_mask = subset_count - 1U;
                for (std::size_t mask = 0U;
                     mask < subset_count && !found; ++mask) {
                    const std::size_t complement = full_mask ^ mask;
                    if (mask > complement) {
                        continue;
                    }
                    for (int source = 1;
                         source <= maximum_label && !found; ++source) {
                        if (std::find(word.begin(), word.end(), source)
                            != word.end()) {
                            continue;
                        }
                        for (int target = 1;
                             target <= maximum_label; ++target) {
                            if (std::find(word.begin(), word.end(), target)
                                != word.end()) {
                                continue;
                            }
                            ++packets;
                            const cpp_int value = mixed_minor(
                                products[mask], products[complement],
                                target, source
                            );
                            if (value >= 0) {
                                continue;
                            }
                            cpp_int doubled_total = 0;
                            std::size_t negative_packets = 0U;
                            for (std::size_t other = 0U;
                                 other < subset_count; ++other) {
                                const cpp_int other_value = mixed_minor(
                                    products[other],
                                    products[full_mask ^ other],
                                    target,
                                    source
                                );
                                doubled_total += other_value;
                                if (other_value < 0) {
                                    ++negative_packets;
                                }
                            }
                            if ((doubled_total & 1) != 0) {
                                throw std::runtime_error(
                                    "odd doubled exterior total"
                                );
                            }
                            std::cout
                                << "SU2_MIXED_EXTERIOR_PACKET result=FAIL"
                                << " word=";
                            print_word(word);
                            std::cout << " first=";
                            print_subset(word, mask);
                            std::cout << " second=";
                            print_subset(word, complement);
                            std::cout << " source=" << source
                                      << " target=" << target
                                      << " packet=" << value
                                      << " total=" << doubled_total / 2
                                      << " negative_oriented_packets="
                                      << negative_packets
                                      << " words=" << words
                                      << " packets=" << packets << '\n';
                            found = true;
                            break;
                        }
                    }
                }
            }
            if (found
                || static_cast<int>(word.size()) == maximum_factors) {
                return;
            }
            for (int label = first_label;
                 label <= maximum_label && !found; ++label) {
                word.push_back(label);
                self(self, label);
                word.pop_back();
            }
        };
        visit(visit, 1);

        if (!found) {
            std::cout << "SU2_MIXED_EXTERIOR_PACKET result=PASS"
                      << " maximum_label=" << maximum_label
                      << " maximum_factors=" << maximum_factors
                      << " words=" << words
                      << " packets=" << packets << '\n';
        }
        return found ? EXIT_FAILURE : EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
