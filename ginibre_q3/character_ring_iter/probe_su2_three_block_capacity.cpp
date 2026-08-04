#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Count = unsigned long long;

int parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0U;
    const long parsed = std::stol(value, &consumed);
    if (consumed != value.size() || parsed <= 0L || parsed > 12L) {
        throw std::invalid_argument(name + " must be an integer in [1,12]");
    }
    return static_cast<int>(parsed);
}

Count three_step_count(
    int start,
    int endpoint,
    const std::vector<int>& labels,
    int state_limit
) {
    std::vector<Count> current(
        static_cast<std::size_t>(state_limit + 1), 0U
    );
    current[static_cast<std::size_t>(start)] = 1U;
    for (const int label : labels) {
        std::vector<Count> next(
            static_cast<std::size_t>(state_limit + 1), 0U
        );
        for (int source = 0; source <= state_limit; ++source) {
            const Count paths = current[static_cast<std::size_t>(source)];
            if (paths == 0U) {
                continue;
            }
            for (int target = std::abs(source - label);
                 target <= source + label;
                 ++target) {
                if (target <= state_limit) {
                    next[static_cast<std::size_t>(target)] += paths;
                }
            }
        }
        current = std::move(next);
    }
    return current[static_cast<std::size_t>(endpoint)];
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::invalid_argument(
                "usage: probe_su2_three_block_capacity MAX_LABEL MAX_STATE"
            );
        }
        const int maximum_label = parse_positive(argv[1], "MAX_LABEL");
        const int maximum_state = parse_positive(argv[2], "MAX_STATE");
        unsigned long long checked = 0U;
        for (int first = 1; first <= maximum_label; ++first) {
            for (int second = 1; second <= maximum_label; ++second) {
                for (int third = 1; third <= maximum_label; ++third) {
                    const std::vector<int> labels{first, second, third};
                    const int state_limit = maximum_state + first + second + third;
                    for (int upper_start = 1;
                         upper_start <= maximum_state;
                         ++upper_start) {
                        for (int lower_start = 0;
                             lower_start < upper_start;
                             ++lower_start) {
                            for (int lower_end = 0;
                                 lower_end <= maximum_state;
                                 ++lower_end) {
                                for (int upper_end = lower_end + 1;
                                     upper_end <= maximum_state;
                                     ++upper_end) {
                                    const Count straight =
                                        three_step_count(
                                            upper_start,
                                            lower_end,
                                            labels,
                                            state_limit
                                        )
                                        * three_step_count(
                                            lower_start,
                                            upper_end,
                                            labels,
                                            state_limit
                                        );
                                    if (straight == 0U) {
                                        continue;
                                    }
                                    const Count crossed =
                                        three_step_count(
                                            lower_start,
                                            lower_end,
                                            labels,
                                            state_limit
                                        )
                                        * three_step_count(
                                            upper_start,
                                            upper_end,
                                            labels,
                                            state_limit
                                        );
                                    ++checked;
                                    if (crossed < straight) {
                                        std::cout
                                            << "SU2_THREE_BLOCK_CAPACITY_DEFICIT"
                                            << " labels=(" << first << ','
                                            << second << ',' << third << ')'
                                            << " starts=(" << upper_start << ','
                                            << lower_start << ')'
                                            << " ends=(" << lower_end << ','
                                            << upper_end << ')'
                                            << " straight=" << straight
                                            << " crossed=" << crossed
                                            << " checked=" << checked
                                            << " result=PASS"
                                            << '\n';
                                        return EXIT_SUCCESS;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << "SU2_THREE_BLOCK_CAPACITY_SCAN"
                  << " maximum_label=" << maximum_label
                  << " maximum_state=" << maximum_state
                  << " checked=" << checked
                  << " result=NO_DEFICIT"
                  << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
