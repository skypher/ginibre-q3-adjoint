#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Multiplicities = std::vector<Integer>;

void for_each_output(
    int left,
    int right,
    int level,
    const auto& function
) {
    const int upper = level < 0
        ? left + right
        : std::min(left + right, 2 * level - left - right);
    for (int output = std::abs(left - right);
         output <= upper; output += 2) {
        function(output);
    }
}

Multiplicities tensor_multiplicities(
    const std::vector<int>& labels,
    std::uint64_t mask,
    bool selected,
    int level
) {
    int total = 0;
    for (std::size_t index = 0U; index < labels.size(); ++index) {
        if ((((mask >> index) & 1U) != 0U) == selected) {
            total += labels[index];
        }
    }
    const int maximum_state = level < 0 ? total : level;
    Multiplicities current(
        static_cast<std::size_t>(maximum_state + 1), 0
    );
    current[0] = 1;
    int current_total = 0;
    for (std::size_t index = 0U; index < labels.size(); ++index) {
        if ((((mask >> index) & 1U) != 0U) != selected) {
            continue;
        }
        Multiplicities next(current.size(), 0);
        const int maximum_input = std::min(current_total, maximum_state);
        for (int input = 0; input <= maximum_input; ++input) {
            const Integer& coefficient
                = current[static_cast<std::size_t>(input)];
            if (coefficient == 0) {
                continue;
            }
            for_each_output(input, labels[index], level, [&](int output) {
                next[static_cast<std::size_t>(output)] += coefficient;
            });
        }
        current = std::move(next);
        current_total += labels[index];
    }
    return current;
}

Integer at(const Multiplicities& values, int label) {
    if (label < 0 || label >= static_cast<int>(values.size())) {
        return 0;
    }
    return values[static_cast<std::size_t>(label)];
}

void print_labels(const std::vector<int>& labels) {
    std::cout << '[';
    for (std::size_t index = 0U; index < labels.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << labels[index];
    }
    std::cout << ']';
}

struct CheckResult {
    bool stronger = true;
    bool required = true;
    int failed_target = -1;
    Integer a_1p = 0;
    Integer a_0_previous = 0;
    Integer a_0_next = 0;
    Integer wall_slack = 0;
    Integer wall_demand = 0;
    Integer wall_supply = 0;
};

CheckResult inspect(const std::vector<int>& labels, int level = -1) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    const int maximum_state = level < 0 ? total : level;
    std::vector<Integer> a_0(
        static_cast<std::size_t>(maximum_state + 2), 0
    );
    std::vector<Integer> a_1(
        static_cast<std::size_t>(maximum_state + 1), 0
    );
    const std::uint64_t subset_count
        = std::uint64_t{1} << labels.size();
    for (std::uint64_t mask = 0; mask < subset_count; ++mask) {
        const Multiplicities first
            = tensor_multiplicities(labels, mask, true, level);
        const Multiplicities second
            = tensor_multiplicities(labels, mask, false, level);
        const Integer first_zero = at(first, 0);
        const Integer first_one = at(first, 1);
        for (int target = 0; target <= maximum_state; ++target) {
            a_0[static_cast<std::size_t>(target)]
                += first_zero * at(second, target);
            a_1[static_cast<std::size_t>(target)]
                += first_one * at(second, target);
        }
    }

    CheckResult result;
    if (level >= 1) {
        result.wall_demand = a_1[static_cast<std::size_t>(level)];
        result.wall_supply = a_0[static_cast<std::size_t>(level - 1)];
        result.wall_slack = result.wall_supply - result.wall_demand;
    }
    for (int target = 0; target <= maximum_state; ++target) {
        const Integer previous
            = target == 0 ? Integer{0}
                          : a_0[static_cast<std::size_t>(target - 1)];
        const Integer next = a_0[static_cast<std::size_t>(target + 1)];
        const Integer demand = a_1[static_cast<std::size_t>(target)];
        if (demand > next && result.stronger) {
            result.stronger = false;
            result.failed_target = target;
            result.a_1p = demand;
            result.a_0_previous = previous;
            result.a_0_next = next;
        }
        if (demand > previous + next) {
            result.required = false;
            result.failed_target = target;
            result.a_1p = demand;
            result.a_0_previous = previous;
            result.a_0_next = next;
            return result;
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4 && std::string(argv[1]) == "finite") {
            const int maximum_level = std::stoi(argv[2]);
            const int maximum_factors = std::stoi(argv[3]);
            if (maximum_level < 2 || maximum_factors < 0
                || maximum_factors >= 63) {
                throw std::runtime_error("invalid finite search bounds");
            }
            std::uint64_t cases = 0U;
            std::uint64_t active_wall_cases = 0U;
            std::uint64_t tight_wall_cases = 0U;
            std::uint64_t simple_current_cases = 0U;
            std::uint64_t printed_tight = 0U;
            std::uint64_t printed_positive = 0U;
            Integer minimum_wall_slack = 0;
            bool have_wall_slack = false;
            for (int level = 2; level <= maximum_level; ++level) {
                std::vector<int> labels;
                const auto visit = [&](const auto& self, int first) -> bool {
                    const CheckResult result = inspect(labels, level);
                    ++cases;
                    if (!result.required) {
                        std::cout
                            << "SU2_LEVEL_FUNDAMENTAL_TRANSFER required=FAIL"
                            << " level=" << level << " labels=";
                        print_labels(labels);
                        std::cout << " target=" << result.failed_target
                                  << " A_1p=" << result.a_1p
                                  << " A_0_previous="
                                  << result.a_0_previous
                                  << " A_0_next=" << result.a_0_next
                                  << '\n';
                        return false;
                    }
                    if (std::binary_search(
                            labels.begin(), labels.end(), level
                        )) {
                        ++simple_current_cases;
                        if (result.wall_slack != 0) {
                            std::cout
                                << "SU2_LEVEL_FUNDAMENTAL_TRANSFER "
                                << "simple_current_wall=FAIL level="
                                << level << " labels=";
                            print_labels(labels);
                            std::cout << " demand="
                                      << result.wall_demand
                                      << " supply="
                                      << result.wall_supply
                                      << " slack="
                                      << result.wall_slack << '\n';
                            return false;
                        }
                    }
                    if (!labels.empty() && result.wall_demand != 0) {
                        ++active_wall_cases;
                        if (!have_wall_slack
                            || result.wall_slack < minimum_wall_slack) {
                            minimum_wall_slack = result.wall_slack;
                            have_wall_slack = true;
                        }
                        if (result.wall_slack == 0) {
                            ++tight_wall_cases;
                            if (printed_tight < 12U) {
                                std::cout << "tight_wall level=" << level
                                          << " labels=";
                                print_labels(labels);
                                std::cout << " demand="
                                          << result.wall_demand
                                          << " supply="
                                          << result.wall_supply << '\n';
                                ++printed_tight;
                            }
                        } else if (printed_positive < 12U) {
                            std::cout << "positive_wall level=" << level
                                      << " labels=";
                            print_labels(labels);
                            std::cout << " demand="
                                      << result.wall_demand
                                      << " supply="
                                      << result.wall_supply
                                      << " slack="
                                      << result.wall_slack << '\n';
                            ++printed_positive;
                        }
                    }
                    if (labels.size()
                            == static_cast<std::size_t>(maximum_factors)) {
                        return true;
                    }
                    for (int label = first; label <= level; ++label) {
                        labels.push_back(label);
                        if (!self(self, label)) {
                            return false;
                        }
                        labels.pop_back();
                    }
                    return true;
                };
                if (!visit(visit, 2)) {
                    return EXIT_FAILURE;
                }
                std::cout << "progress level=" << level
                          << " cases=" << cases
                          << " active_wall_cases=" << active_wall_cases
                          << " tight_wall_cases=" << tight_wall_cases
                          << " simple_current_cases="
                          << simple_current_cases
                          << '\n';
            }
            std::cout << "SU2_LEVEL_FUNDAMENTAL_TRANSFER levels=2.."
                      << maximum_level
                      << " maximum_factors=" << maximum_factors
                      << " cases=" << cases
                      << " active_wall_cases=" << active_wall_cases
                      << " tight_wall_cases=" << tight_wall_cases
                      << " simple_current_cases=" << simple_current_cases
                      << " minimum_wall_slack="
                      << (have_wall_slack ? minimum_wall_slack : Integer{0})
                      << " required=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc != 4) {
            throw std::runtime_error(
                "usage: search_su2_fundamental_transfer "
                "MAXIMUM_LABEL MAXIMUM_FACTORS MAXIMUM_TOTAL\n"
                "   or: search_su2_fundamental_transfer finite "
                "MAXIMUM_LEVEL MAXIMUM_FACTORS"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        const int maximum_total = std::stoi(argv[3]);
        if (maximum_label < 2 || maximum_factors < 0
            || maximum_factors >= 63 || maximum_total < 0) {
            throw std::runtime_error("invalid search bounds");
        }

        std::uint64_t cases = 0U;
        std::uint64_t stronger_failures = 0U;
        std::vector<int> labels;
        const auto visit = [&](const auto& self, int first, int total) -> bool {
            const CheckResult result = inspect(labels);
            ++cases;
            if (!result.required) {
                std::cout << "SU2_FUNDAMENTAL_TRANSFER required=FAIL labels=";
                print_labels(labels);
                std::cout << " target=" << result.failed_target
                          << " A_1p=" << result.a_1p
                          << " A_0_previous=" << result.a_0_previous
                          << " A_0_next=" << result.a_0_next << '\n';
                return false;
            }
            if (!result.stronger) {
                ++stronger_failures;
                if (stronger_failures <= 12U) {
                    std::cout << "stronger_failure labels=";
                    print_labels(labels);
                    std::cout << " target=" << result.failed_target
                              << " A_1p=" << result.a_1p
                              << " A_0_previous=" << result.a_0_previous
                              << " A_0_next=" << result.a_0_next << '\n';
                }
            }
            if (labels.size()
                    == static_cast<std::size_t>(maximum_factors)) {
                return true;
            }
            for (int label = first; label <= maximum_label; ++label) {
                if (total + label > maximum_total) {
                    break;
                }
                labels.push_back(label);
                if (!self(self, label, total + label)) {
                    return false;
                }
                labels.pop_back();
            }
            return true;
        };

        if (!visit(visit, 2, 0)) {
            return EXIT_FAILURE;
        }
        std::cout << "SU2_FUNDAMENTAL_TRANSFER cases=" << cases
                  << " stronger_failures=" << stronger_failures
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors
                  << " maximum_total=" << maximum_total
                  << " required=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
