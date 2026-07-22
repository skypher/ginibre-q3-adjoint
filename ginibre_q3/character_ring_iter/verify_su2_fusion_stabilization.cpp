#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

using boost::multiprecision::cpp_int;

namespace {

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
void for_each_output(int level, int left, int right, Function function) {
    int upper = left + right;
    if (level >= 0) {
        if (left > level || right > level) {
            throw std::runtime_error("fusion label exceeds the level");
        }
        upper = std::min(upper, 2 * level - left - right);
    }
    for (int output = std::abs(left - right); output <= upper; output += 2) {
        function(output);
    }
}

cpp_int q3_dynamic(
    int level,
    const std::vector<int>& minus,
    const std::vector<int>& plus
) {
    std::unordered_map<std::uint64_t, cpp_int> states;
    states.emplace(pair_key(0, 0), 1);

    auto apply = [&states, level](int label, int sign) {
        std::unordered_map<std::uint64_t, cpp_int> next;
        for (const auto& [key, coefficient] : states) {
            const auto [left, right] = decode_pair(key);
            for_each_output(level, left, label, [&](int output) {
                next[pair_key(output, right)] += coefficient;
            });
            for_each_output(level, right, label, [&](int output) {
                next[pair_key(left, output)] += sign * coefficient;
            });
        }
        states = std::move(next);
    };

    for (int label : minus) {
        apply(label, -1);
    }
    for (int label : plus) {
        apply(label, 1);
    }
    const auto found = states.find(pair_key(0, 0));
    return found == states.end() ? cpp_int(0) : found->second;
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

struct Case {
    std::vector<int> minus;
    std::vector<int> plus;
};

std::vector<Case> make_cases(int maximum_label, int maximum_factors) {
    std::vector<Case> cases;
    for (int factors = 2; factors <= maximum_factors; ++factors) {
        for (int minus_count = 2; minus_count <= factors; minus_count += 2) {
            const auto minus_lists = combinations(maximum_label, minus_count);
            const auto plus_lists = combinations(maximum_label, factors - minus_count);
            for (const auto& minus : minus_lists) {
                for (const auto& plus : plus_lists) {
                    cases.push_back(Case{minus, plus});
                }
            }
        }
    }
    return cases;
}

int stabilization_level(const Case& test_case) {
    int total = 0;
    int maximum = 0;
    for (int label : test_case.minus) {
        total += label;
        maximum = std::max(maximum, label);
    }
    for (int label : test_case.plus) {
        total += label;
        maximum = std::max(maximum, label);
    }
    return std::max(maximum, (total + 1) / 2);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_label = argc > 1 ? std::atoi(argv[1]) : 6;
        const int maximum_factors = argc > 2 ? std::atoi(argv[2]) : 8;
        if (maximum_label < 1 || maximum_factors < 2) {
            throw std::runtime_error("expected maximum_label >= 1 and maximum_factors >= 2");
        }

        const std::vector<Case> cases = make_cases(maximum_label, maximum_factors);
        std::int64_t failures = 0;
        std::int64_t checked = 0;

        #pragma omp parallel for schedule(dynamic) reduction(+:failures,checked)
        for (std::int64_t raw_index = 0;
             raw_index < static_cast<std::int64_t>(cases.size());
             ++raw_index) {
            const auto& test_case = cases[static_cast<std::size_t>(raw_index)];
            const int level = stabilization_level(test_case);
            const cpp_int ordinary = q3_dynamic(-1, test_case.minus, test_case.plus);
            const cpp_int fusion = q3_dynamic(level, test_case.minus, test_case.plus);
            const cpp_int next_level = q3_dynamic(level + 1, test_case.minus, test_case.plus);
            if (ordinary != fusion || ordinary != next_level) {
                ++failures;
            }
            ++checked;
        }

        std::cout << "SU2_FUSION_STABILIZATION"
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors
                  << " cases=" << checked
                  << " failures=" << failures
                  << (failures == 0 ? " PASS\n" : " FAIL\n");
        return failures == 0 ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
