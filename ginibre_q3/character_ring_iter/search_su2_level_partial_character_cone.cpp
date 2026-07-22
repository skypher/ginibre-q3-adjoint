#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

template <class Function>
void for_each_output(int level, int left, int right, Function function) {
    const int lower = std::abs(left - right);
    const int upper = std::min(left + right, 2 * level - left - right);
    for (int output = lower; output <= upper; output += 2) {
        function(output);
    }
}

std::vector<std::vector<int>> multisets(int level, int maximum_size) {
    std::vector<std::vector<int>> result;
    std::vector<int> current;
    const auto visit = [&](const auto& self, int first) -> void {
        if (!current.empty()) {
            result.push_back(current);
        }
        if (current.size() == static_cast<std::size_t>(maximum_size)) {
            return;
        }
        for (int label = first; label <= level; ++label) {
            current.push_back(label);
            self(self, label);
            current.pop_back();
        }
    };
    visit(visit, 1);
    return result;
}

bool inspect_word(
    int level,
    const std::vector<int>& labels,
    std::uint64_t minus_mask,
    int& failed_target,
    cpp_int& failed_value
) {
    const int width = level + 1;
    std::vector<cpp_int> current(
        static_cast<std::size_t>(width * width)
    );
    current[0] = 1;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        const int label = labels[index];
        const bool minus = ((minus_mask >> index) & 1U) != 0U;
        std::vector<cpp_int> next(current.size());
        for (int left = 0; left <= level; ++left) {
            for (int right = 0; right <= level; ++right) {
                const cpp_int& value = current[static_cast<std::size_t>(
                    left * width + right
                )];
                if (value == 0) {
                    continue;
                }
                for_each_output(level, left, label, [&](int output) {
                    next[static_cast<std::size_t>(
                        output * width + right
                    )] += value;
                });
                for_each_output(level, right, label, [&](int output) {
                    next[static_cast<std::size_t>(
                        left * width + output
                    )] += minus ? -value : value;
                });
            }
        }
        current = std::move(next);
    }
    for (int target = 0; target <= level; ++target) {
        const cpp_int& value = current[static_cast<std::size_t>(
            target * width
        )];
        if (value < 0) {
            failed_target = target;
            failed_value = value;
            return false;
        }
    }
    return true;
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
        if (argc != 3) {
            throw std::runtime_error(
                "usage: search_su2_level_partial_character_cone "
                "MAXIMUM_LEVEL MAXIMUM_FACTORS"
            );
        }
        const int maximum_level = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_level < 1 || maximum_factors < 1
            || maximum_factors >= 63) {
            throw std::runtime_error("invalid search bound");
        }

        std::uint64_t total_tested = 0;
        for (int level = 1; level <= maximum_level; ++level) {
            const auto words = multisets(level, maximum_factors);
            std::uint64_t level_tested = 0;
            for (const auto& labels : words) {
                const std::uint64_t masks
                    = std::uint64_t{1} << labels.size();
                for (std::uint64_t mask = 0; mask < masks; ++mask) {
                    int failed_target = -1;
                    cpp_int failed_value = 0;
                    ++level_tested;
                    if (!inspect_word(
                            level, labels, mask,
                            failed_target, failed_value
                        )) {
                        std::cout << "FAIL level=" << level << " factors=";
                        print_word(labels, mask);
                        std::cout << " target=" << failed_target
                                  << " coefficient=" << failed_value << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
            total_tested += level_tested;
            std::cout << "SU2_LEVEL_PARTIAL_CHARACTER_CONE level=" << level
                      << " tested=" << level_tested << " PASS\n";
        }
        std::cout << "SU2_LEVEL_PARTIAL_CHARACTER_CONE total_tested="
                  << total_tested << " maximum_level=" << maximum_level
                  << " maximum_factors=" << maximum_factors << " PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
