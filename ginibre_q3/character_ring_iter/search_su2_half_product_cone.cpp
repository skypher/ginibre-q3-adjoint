#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <omp.h>

namespace {

using Coefficient = std::int64_t;
using Expansion = std::unordered_map<std::uint64_t, Coefficient>;

std::uint64_t pair_key(int left, int right) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(left)) << 32U)
        | static_cast<std::uint32_t>(right);
}

std::pair<int, int> decode_pair(std::uint64_t key) {
    return {
        static_cast<int>(static_cast<std::uint32_t>(key >> 32U)),
        static_cast<int>(static_cast<std::uint32_t>(key))
    };
}

template <class Function>
void for_each_output(int left, int right, Function function) {
    for (int output = std::abs(left - right); output <= left + right; output += 2) {
        function(output);
    }
}

Expansion half_product(const std::vector<int>& labels) {
    Expansion states;
    states.emplace(pair_key(0, 0), 1);
    for (int label : labels) {
        Expansion next;
        for (const auto& [key, coefficient] : states) {
            const auto [left, right] = decode_pair(key);
            for_each_output(left, label, [&](int output) {
                next[pair_key(output, right)] += coefficient;
            });
            for_each_output(right, label, [&](int output) {
                next[pair_key(left, output)] -= coefficient;
            });
        }
        states = std::move(next);
    }
    return states;
}

Coefficient coefficient_at(
    const std::vector<int>& labels,
    int left,
    int right
) {
    const Expansion expansion = half_product(labels);
    const auto found = expansion.find(pair_key(left, right));
    return found == expansion.end() ? Coefficient{0} : found->second;
}

Coefficient inner_product(
    const std::vector<int>& left_labels,
    const std::vector<int>& right_labels
) {
    const Expansion left = half_product(left_labels);
    const Expansion right = half_product(right_labels);
    Coefficient answer = 0;
    for (const auto& [key, coefficient] : left) {
        const auto found = right.find(key);
        if (found != right.end()) {
            answer += coefficient * found->second;
        }
    }
    return answer;
}

bool checkerboard_wedge_positive(
    const std::vector<int>& labels,
    int& witness_left,
    int& witness_right,
    Coefficient& witness_coefficient
) {
    if ((labels.size() & 1U) == 0U) {
        throw std::runtime_error("wedge test requires an odd number of factors");
    }
    const Expansion expansion = half_product(labels);
    for (const auto& [key, coefficient] : expansion) {
        const auto [left, right] = decode_pair(key);
        if (left <= right) {
            continue;
        }
        const Coefficient adjusted = (right & 1) == 0 ? coefficient : -coefficient;
        if (adjusted < 0) {
            witness_left = left;
            witness_right = right;
            witness_coefficient = coefficient;
            return false;
        }
    }
    return true;
}

void combinations_rec(
    int maximum_label,
    int remaining,
    int first,
    std::vector<int>& current,
    std::vector<std::vector<int>>& output
) {
    if (remaining == 0) {
        output.push_back(current);
        return;
    }
    for (int label = first; label <= maximum_label; ++label) {
        current.push_back(label);
        combinations_rec(maximum_label, remaining - 1, label, current, output);
        current.pop_back();
    }
}

std::vector<std::vector<int>> combinations(int maximum_label, int size) {
    std::vector<std::vector<int>> output;
    std::vector<int> current;
    combinations_rec(maximum_label, size, 1, current, output);
    return output;
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

int parse_positive(const char* text, const char* name) {
    const long value = std::strtol(text, nullptr, 10);
    if (value <= 0 || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_label = argc > 1
            ? parse_positive(argv[1], "maximum label") : 12;
        const int factor_count = argc > 2
            ? parse_positive(argv[2], "factor count") : 3;
        if (argc > 3 || (factor_count & 1) == 0) {
            throw std::runtime_error(
                "usage: search_su2_half_product_cone [maximum_label odd_factor_count]"
            );
        }

        const auto tuples = combinations(maximum_label, factor_count);
        std::int64_t failures = 0;
        std::int64_t first_failure = -1;
        int witness_left = 0;
        int witness_right = 0;
        Coefficient witness_coefficient = 0;

        if (factor_count == 3 && maximum_label >= 2) {
            std::cout << "SU2_HALF_PRODUCT_CONE probes"
                      << " labels=[1,1,1] coefficient=(2,1)->"
                      << coefficient_at({1, 1, 1}, 2, 1)
                      << " labels=[1,2,2] coefficient=(2,1)->"
                      << coefficient_at({1, 2, 2}, 2, 1)
                      << " cross_inner_product="
                      << inner_product({1, 1, 1}, {1, 2, 2}) << '\n';
        }

        #pragma omp parallel for schedule(dynamic) reduction(+:failures)
        for (std::int64_t raw_index = 0;
             raw_index < static_cast<std::int64_t>(tuples.size());
             ++raw_index) {
            int local_left = 0;
            int local_right = 0;
            Coefficient local_coefficient = 0;
            if (!checkerboard_wedge_positive(
                    tuples[static_cast<std::size_t>(raw_index)],
                    local_left, local_right, local_coefficient)) {
                ++failures;
                #pragma omp critical
                {
                    if (first_failure < 0 || raw_index < first_failure) {
                        first_failure = raw_index;
                        witness_left = local_left;
                        witness_right = local_right;
                        witness_coefficient = local_coefficient;
                    }
                }
            }
        }

        std::cout << "SU2_HALF_PRODUCT_CONE"
                  << " maximum_label=" << maximum_label
                  << " factor_count=" << factor_count
                  << " cases=" << tuples.size()
                  << " failures=" << failures;
        if (first_failure >= 0) {
            std::cout << " first_labels=";
            print_labels(tuples[static_cast<std::size_t>(first_failure)]);
            std::cout << " coefficient=(" << witness_left << ',' << witness_right
                      << ")->" << witness_coefficient;
        }
        std::cout << (failures == 0 ? " PASS\n" : " FAIL\n");
        return failures == 0 ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
