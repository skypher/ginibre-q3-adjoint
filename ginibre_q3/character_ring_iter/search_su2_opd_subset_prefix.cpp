#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using MultiplicityKey = std::pair<std::vector<int>, int>;

cpp_int multiplicity(std::vector<int> labels, int target) {
    static std::map<MultiplicityKey, cpp_int> cache;
    std::sort(labels.begin(), labels.end());
    const MultiplicityKey key{labels, target};
    const auto found = cache.find(key);
    if (found != cache.end()) {
        return found->second;
    }
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    if (target < 0 || target > total) {
        return 0;
    }
    std::vector<cpp_int> current(static_cast<std::size_t>(total + 1));
    current[0] = 1;
    int support = 0;
    for (int label : labels) {
        std::vector<cpp_int> next(current.size());
        for (int input = 0; input <= support; ++input) {
            const cpp_int& coefficient
                = current[static_cast<std::size_t>(input)];
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
    const cpp_int result = current[static_cast<std::size_t>(target)];
    cache.emplace(key, result);
    return result;
}

std::vector<cpp_int> layers(
    int p,
    const std::vector<int>& suffix,
    int target
) {
    const int count = static_cast<int>(suffix.size());
    std::vector<cpp_int> result(static_cast<std::size_t>(count + 1));
    const std::uint64_t masks = std::uint64_t{1} << count;
    for (std::uint64_t mask = 0; mask < masks; ++mask) {
        std::vector<int> first;
        std::vector<int> second;
        for (int index = 0; index < count; ++index) {
            if (((mask >> index) & 1U) != 0U) {
                second.push_back(suffix[static_cast<std::size_t>(index)]);
            } else {
                first.push_back(suffix[static_cast<std::size_t>(index)]);
            }
        }
        std::vector<int> first_p_minus = first;
        first_p_minus.push_back(p - 1);
        std::vector<int> first_p = first;
        first_p.push_back(p);
        std::vector<int> first_one = first;
        first_one.push_back(1);
        const cpp_int positive = multiplicity(second, 0)
                * multiplicity(first_p_minus, target)
            + multiplicity(second, p)
                * multiplicity(first_one, target);
        const cpp_int negative = multiplicity(second, p - 1)
                * multiplicity(first, target)
            + multiplicity(second, 1)
                * multiplicity(first_p, target);
        const int layer = __builtin_popcountll(mask);
        result[static_cast<std::size_t>(layer)] += positive - negative;
    }
    return result;
}

void print_word(int p, const std::vector<int>& suffix, int target) {
    std::cout << "p=" << p << " suffix=[";
    for (std::size_t index = 0; index < suffix.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << suffix[index];
    }
    std::cout << "] target=" << target;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 4) {
            throw std::runtime_error(
                "usage: search_su2_opd_subset_prefix MAXIMUM_LABEL "
                "MAXIMUM_SUFFIX_FACTORS [prefix|shell]"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        const std::string mode = argc == 4 ? argv[3] : "prefix";
        if (maximum_label < 2 || maximum_factors < 0
            || maximum_factors >= 63
            || (mode != "prefix" && mode != "shell")) {
            throw std::runtime_error("invalid bound");
        }
        long long words = 0;
        long long cases = 0;
        for (int p = 2; p <= maximum_label; ++p) {
            std::vector<int> suffix;
            const auto visit = [&](const auto& self, int first) -> void {
                ++words;
                int total = p;
                for (int label : suffix) {
                    total += label;
                }
                for (int target = 1; target <= total; ++target) {
                    ++cases;
                    const auto current = layers(p, suffix, target);
                    cpp_int cumulative = 0;
                    if (mode == "prefix") {
                        for (std::size_t layer = 0;
                             layer < current.size();
                             ++layer) {
                            cumulative += current[layer];
                            if (cumulative < 0) {
                                std::cout << "FAIL mode=prefix layer=" << layer
                                          << " value=" << cumulative << ' ';
                                print_word(p, suffix, target);
                                std::cout << '\n';
                                std::exit(EXIT_FAILURE);
                            }
                        }
                    } else {
                        const std::size_t count = suffix.size();
                        for (std::size_t depth = 0;
                             depth <= count / 2U;
                             ++depth) {
                            cumulative += current[depth];
                            const std::size_t complement = count - depth;
                            if (complement != depth) {
                                cumulative += current[complement];
                            }
                            if (cumulative < 0) {
                                std::cout << "FAIL mode=shell depth=" << depth
                                          << " value=" << cumulative << ' ';
                                print_word(p, suffix, target);
                                std::cout << '\n';
                                std::exit(EXIT_FAILURE);
                            }
                        }
                    }
                }
                if (suffix.size()
                    == static_cast<std::size_t>(maximum_factors)) {
                    return;
                }
                for (int label = first; label <= maximum_label; ++label) {
                    suffix.push_back(label);
                    self(self, label);
                    suffix.pop_back();
                }
            };
            visit(visit, p);
        }
        std::cout << "SU2_OPD_SUBSET_";
        if (mode == "prefix") {
            std::cout << "PREFIX";
        } else {
            std::cout << "SHELL";
        }
        std::cout << " PASS words=" << words
                  << " cases=" << cases
                  << " maximum_label=" << maximum_label
                  << " maximum_suffix_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
