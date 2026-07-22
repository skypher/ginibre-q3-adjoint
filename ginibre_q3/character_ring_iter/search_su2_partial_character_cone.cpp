#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

bool inspect_word(
    const std::vector<int>& labels,
    std::uint64_t minus_mask,
    int& failed_target,
    cpp_int& failed_value
) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
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
    for (int target = 0; target <= total; ++target) {
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
                "usage: search_su2_partial_character_cone MAXIMUM_LABEL "
                "MAXIMUM_FACTORS"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 1 || maximum_factors < 1
            || maximum_factors >= 63) {
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
                const std::uint64_t masks
                    = std::uint64_t{1} << labels.size();
                for (std::uint64_t mask = 0; mask < masks; ++mask) {
                    int failed_target = -1;
                    cpp_int failed_value = 0;
                    ++tested;
                    if (!inspect_word(
                            labels, mask, failed_target, failed_value
                        )) {
                        std::cout << "FAIL factors=";
                        print_word(labels, mask);
                        std::cout << " target=" << failed_target
                                  << " coefficient=" << failed_value
                                  << " tested=" << tested << '\n';
                        failed = true;
                        return;
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
            return EXIT_FAILURE;
        }
        std::cout << "SU2_PARTIAL_CHARACTER_CONE PASS tested=" << tested
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
