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

using State = std::unordered_map<std::uint64_t, cpp_int>;

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

State signed_product(
    int level,
    const std::vector<int>& minus,
    const std::vector<int>& plus
) {
    State states;
    states.emplace(pair_key(0, 0), 1);

    auto apply = [&states, level](int label, int sign) {
        State next;
        for (const auto& [key, coefficient] : states) {
            const auto [left, right] = decode_pair(key);
            for_each_output(level, left, label, [&](int output) {
                next[pair_key(output, right)] += coefficient;
            });
            for_each_output(level, right, label, [&](int output) {
                next[pair_key(left, output)] += sign * coefficient;
            });
        }
        std::erase_if(next, [](const auto& entry) {
            return entry.second == 0;
        });
        states = std::move(next);
    };

    for (int label : minus) {
        apply(label, -1);
    }
    for (int label : plus) {
        apply(label, 1);
    }
    return states;
}

cpp_int corner(const State& state) {
    const auto found = state.find(pair_key(0, 0));
    return found == state.end() ? cpp_int(0) : found->second;
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
    for (int factors = 0; factors <= maximum_factors; ++factors) {
        for (int minus_count = 0; minus_count <= factors; ++minus_count) {
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

int total_label(const Case& test_case) {
    int total = 0;
    for (int label : test_case.minus) {
        total += label;
    }
    for (int label : test_case.plus) {
        total += label;
    }
    return total;
}

int scalar_stabilization_level(const Case& test_case) {
    int maximum = 0;
    for (int label : test_case.minus) {
        maximum = std::max(maximum, label);
    }
    for (int label : test_case.plus) {
        maximum = std::max(maximum, label);
    }
    return std::max(maximum, (total_label(test_case) + 1) / 2);
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
        std::int64_t full_state_failures = 0;
        std::int64_t corner_failures = 0;
        std::int64_t checked = 0;

        #pragma omp parallel for schedule(dynamic) \
            reduction(+:full_state_failures,corner_failures,checked)
        for (std::int64_t raw_index = 0;
             raw_index < static_cast<std::int64_t>(cases.size());
             ++raw_index) {
            const auto& test_case = cases[static_cast<std::size_t>(raw_index)];
            const int full_level = std::max(1, total_label(test_case));
            const int scalar_level =
                std::max(1, scalar_stabilization_level(test_case));
            const State ordinary =
                signed_product(-1, test_case.minus, test_case.plus);
            const State fusion_full =
                signed_product(full_level, test_case.minus, test_case.plus);
            const State fusion_full_next =
                signed_product(full_level + 1, test_case.minus, test_case.plus);
            if (ordinary != fusion_full || ordinary != fusion_full_next) {
                ++full_state_failures;
            }
            const State fusion_scalar =
                signed_product(scalar_level, test_case.minus, test_case.plus);
            const State fusion_scalar_next =
                signed_product(
                    scalar_level + 1, test_case.minus, test_case.plus);
            if (corner(ordinary) != corner(fusion_scalar) ||
                corner(ordinary) != corner(fusion_scalar_next)) {
                ++corner_failures;
            }
            ++checked;
        }

        std::cout << "SU2_FUSION_STABILIZATION"
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors
                  << " cases=" << checked
                  << " full_state_failures=" << full_state_failures
                  << " corner_failures=" << corner_failures
                  << (full_state_failures == 0 && corner_failures == 0
                          ? " PASS\n"
                          : " FAIL\n");
        return full_state_failures == 0 && corner_failures == 0 ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
