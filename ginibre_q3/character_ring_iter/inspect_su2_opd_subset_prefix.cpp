#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

cpp_int multiplicity(const std::vector<int>& labels, int target) {
    int total = 0;
    for (const int label : labels) {
        total += label;
    }
    if (target < 0 || target > total) {
        return 0;
    }
    std::vector<cpp_int> current(static_cast<std::size_t>(total + 1));
    current[0] = 1;
    int support = 0;
    for (const int label : labels) {
        std::vector<cpp_int> next(current.size());
        for (int input = 0; input <= support; ++input) {
            const cpp_int coefficient = current[static_cast<std::size_t>(input)];
            if (coefficient == 0) {
                continue;
            }
            for (int output = std::abs(input - label);
                 output <= input + label;
                 output += 2) {
                next[static_cast<std::size_t>(output)] += coefficient;
            }
        }
        support += label;
        current = std::move(next);
    }
    return current[static_cast<std::size_t>(target)];
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 4) {
            throw std::runtime_error(
                "usage: inspect_su2_opd_subset_prefix P TARGET SUFFIX..."
            );
        }
        const int p = std::stoi(argv[1]);
        const int target = std::stoi(argv[2]);
        std::vector<int> suffix;
        for (int index = 3; index < argc; ++index) {
            suffix.push_back(std::stoi(argv[index]));
        }
        if (p < 2 || target < 1 || suffix.size() >= 63U
            || std::ranges::any_of(suffix, [p](const int q) { return q < p; })) {
            throw std::runtime_error("invalid ordered-packet data");
        }
        const std::size_t count = suffix.size();
        std::vector<cpp_int> layers(count + 1U);
        const std::uint64_t masks = std::uint64_t{1} << count;
        for (std::uint64_t mask = 0; mask < masks; ++mask) {
            std::vector<int> first;
            std::vector<int> second;
            for (std::size_t index = 0; index < count; ++index) {
                if (((mask >> index) & 1U) == 0U) {
                    first.push_back(suffix[index]);
                } else {
                    second.push_back(suffix[index]);
                }
            }
            auto first_p_minus = first;
            first_p_minus.push_back(p - 1);
            auto first_p = first;
            first_p.push_back(p);
            auto first_one = first;
            first_one.push_back(1);
            const cpp_int value
                = multiplicity(second, 0) * multiplicity(first_p_minus, target)
                - multiplicity(second, p - 1) * multiplicity(first, target)
                - multiplicity(second, 1) * multiplicity(first_p, target)
                + multiplicity(second, p) * multiplicity(first_one, target);
            const auto layer = static_cast<std::size_t>(
                __builtin_popcountll(mask)
            );
            layers[layer] += value;
        }
        cpp_int prefix = 0;
        for (std::size_t layer = 0; layer <= count; ++layer) {
            prefix += layers[layer];
            std::cout << "layer=" << layer << " coefficient=" << layers[layer]
                      << " prefix=" << prefix << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
