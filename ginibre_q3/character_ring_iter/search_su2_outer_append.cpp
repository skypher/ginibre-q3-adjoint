#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

cpp_int partial_coefficient(
    const std::vector<int>& labels,
    std::uint64_t minus_mask,
    int target
) {
    const int total = std::accumulate(labels.begin(), labels.end(), 0);
    const int width = total + 1;
    std::vector<cpp_int> current(
        static_cast<std::size_t>(width * width)
    );
    current[0] = 1;
    int support = 0;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        const int label = labels[index];
        const bool minus = ((minus_mask >> index) & 1U) != 0U;
        std::vector<cpp_int> next(current.size());
        for (int left = 0; left <= support; ++left) {
            for (int right = 0; right <= support; ++right) {
                const cpp_int& value = current[static_cast<std::size_t>(
                    left * width + right
                )];
                if (value == 0) {
                    continue;
                }
                for (int output = std::abs(left - label);
                     output <= left + label;
                     output += 2) {
                    next[static_cast<std::size_t>(
                        output * width + right
                    )] += value;
                }
                for (int output = std::abs(right - label);
                     output <= right + label;
                     output += 2) {
                    next[static_cast<std::size_t>(
                        left * width + output
                    )] += minus ? -value : value;
                }
            }
        }
        support += label;
        current = std::move(next);
    }
    return current[static_cast<std::size_t>(target * width)];
}

void print_word(const std::vector<int>& labels, std::uint64_t minus_mask) {
    std::cout << '[';
    for (std::size_t index = 0; index < labels.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << (((minus_mask >> index) & 1U) != 0U ? '-' : '+')
                  << labels[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4) {
            throw std::runtime_error(
                "usage: search_su2_outer_append MAXIMUM_LABEL "
                "MAXIMUM_FACTORS MAXIMUM_DEFICIT"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        const int maximum_deficit = std::stoi(argv[3]);
        if (maximum_label < 2 || maximum_factors < 1
            || maximum_factors >= 62 || maximum_deficit < 0) {
            throw std::runtime_error("invalid search bound");
        }

        std::uint64_t tested = 0;
        bool failed = false;
        std::vector<int> labels;
        const auto visit = [&](const auto& self, int first) -> void {
            if (failed) {
                return;
            }
            if (!labels.empty()) {
                const int sum = std::accumulate(labels.begin(), labels.end(), 0);
                const std::uint64_t masks
                    = std::uint64_t{1} << labels.size();
                for (int deficit = 0;
                     deficit <= maximum_deficit;
                     ++deficit) {
                    const int old_target = sum - 2 * deficit;
                    if (old_target < labels.back()) {
                        continue;
                    }
                    for (std::uint64_t mask = 0; mask < masks; ++mask) {
                        const cpp_int old_value = partial_coefficient(
                            labels, mask, old_target
                        );
                        for (int appended = 2;
                             appended <= maximum_label;
                             ++appended) {
                            const int new_target = old_target + appended;
                            if (new_target < std::max(labels.back(), appended)) {
                                continue;
                            }
                            std::vector<int> extended = labels;
                            extended.push_back(appended);
                            for (int appended_minus = 0;
                                 appended_minus <= 1;
                                 ++appended_minus) {
                                const std::uint64_t extended_mask = mask
                                    | (static_cast<std::uint64_t>(
                                        appended_minus
                                    ) << labels.size());
                                const cpp_int new_value = partial_coefficient(
                                    extended, extended_mask, new_target
                                );
                                ++tested;
                                if (new_value < old_value) {
                                    std::cout << "FAIL deficit=" << deficit
                                              << " old=";
                                    print_word(labels, mask);
                                    std::cout << " appended="
                                              << (appended_minus != 0 ? '-' : '+')
                                              << appended
                                              << " old_value=" << old_value
                                              << " new_value=" << new_value
                                              << '\n';
                                    failed = true;
                                    return;
                                }
                            }
                        }
                    }
                }
            }
            if (labels.size() == static_cast<std::size_t>(maximum_factors)) {
                return;
            }
            for (int label = first; label <= maximum_label; ++label) {
                labels.push_back(label);
                self(self, label);
                labels.pop_back();
                if (failed) {
                    return;
                }
            }
        };
        visit(visit, 1);
        if (failed) {
            std::cout << "tested=" << tested << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "SU2_OUTER_APPEND PASS tested=" << tested
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors
                  << " maximum_deficit=" << maximum_deficit << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
