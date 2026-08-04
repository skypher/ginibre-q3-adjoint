#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

constexpr int deficit = 12;
constexpr int label_count = 6;
using Counts = std::array<int, label_count>;
using MultiplicityKey = std::array<int, label_count + 1>;
using OuterKey = std::array<int, 2 * label_count>;

std::map<MultiplicityKey, cpp_int>& multiplicity_cache() {
    static std::map<MultiplicityKey, cpp_int> cache;
    return cache;
}

std::map<OuterKey, cpp_int>& outer_cache() {
    static std::map<OuterKey, cpp_int> cache;
    return cache;
}

cpp_int binomial(int top, int bottom) {
    if (bottom < 0 || bottom > top) {
        return 0;
    }
    bottom = std::min(bottom, top - bottom);
    cpp_int result = 1;
    for (int factor = 1; factor <= bottom; ++factor) {
        result *= top - bottom + factor;
        result /= factor;
    }
    return result;
}

int total_label(const Counts& counts) {
    int result = 0;
    for (int index = 0; index < label_count; ++index) {
        result += (index + 1) * counts[static_cast<std::size_t>(index)];
    }
    return result;
}

int total_count(const Counts& counts) {
    return std::accumulate(counts.begin(), counts.end(), 0);
}

int maximum_label(const Counts& counts) {
    for (int index = label_count - 1; index >= 0; --index) {
        if (counts[static_cast<std::size_t>(index)] != 0) {
            return index + 1;
        }
    }
    return 0;
}

bool admissible(const Counts& counts) {
    const int maximum = maximum_label(counts);
    return maximum != 0 && total_label(counts) - deficit >= maximum;
}

cpp_int multiplicity(const Counts& counts, int target) {
    MultiplicityKey key{};
    for (std::size_t index = 0U; index < counts.size(); ++index) {
        key[index] = counts[index];
    }
    key[counts.size()] = target;
    auto& cache = multiplicity_cache();
    const auto found = cache.find(key);
    if (found != cache.end()) {
        return found->second;
    }

    const int total = total_label(counts);
    if (target < 0 || target > total) {
        return 0;
    }
    std::vector<cpp_int> current(static_cast<std::size_t>(total + 1));
    current[0U] = 1;
    int support = 0;
    for (int label = 1; label <= label_count; ++label) {
        for (int copy = 0;
             copy < counts[static_cast<std::size_t>(label - 1)]; ++copy) {
            std::vector<cpp_int> next(current.size());
            for (int input = 0; input <= support; ++input) {
                if (current[static_cast<std::size_t>(input)] == 0) {
                    continue;
                }
                for (int output = std::abs(input - label);
                     output <= input + label; output += 2) {
                    next[static_cast<std::size_t>(output)]
                        += current[static_cast<std::size_t>(input)];
                }
            }
            support += label;
            current = std::move(next);
        }
    }
    const cpp_int result = current[static_cast<std::size_t>(target)];
    cache.emplace(key, result);
    return result;
}

cpp_int outer_value(const Counts& counts, const Counts& signs) {
    OuterKey key{};
    for (std::size_t index = 0U; index < counts.size(); ++index) {
        key[index] = counts[index];
        key[counts.size() + index] = signs[index];
    }
    auto& cache = outer_cache();
    const auto found = cache.find(key);
    if (found != cache.end()) {
        return found->second;
    }

    const int target = total_label(counts) - deficit;
    cpp_int result = 0;
    Counts selected{};
    const auto visit = [&](const auto& self, int coordinate,
                           int selected_label) -> void {
        if (coordinate == label_count) {
            const cpp_int invariant = multiplicity(selected, 0);
            if (invariant == 0) {
                return;
            }
            Counts complement{};
            cpp_int weight = invariant;
            int sign_parity = 0;
            for (int index = 0; index < label_count; ++index) {
                const std::size_t position = static_cast<std::size_t>(index);
                complement[position] = counts[position] - selected[position];
                weight *= binomial(counts[position], selected[position]);
                if (signs[position] < 0) {
                    sign_parity += selected[position];
                }
            }
            weight *= multiplicity(complement, target);
            result += (sign_parity & 1) == 0 ? weight : -weight;
            return;
        }

        const int label = coordinate + 1;
        const int maximum_selected = std::min(
            counts[static_cast<std::size_t>(coordinate)],
            (deficit - selected_label) / label
        );
        for (int amount = 0; amount <= maximum_selected; ++amount) {
            selected[static_cast<std::size_t>(coordinate)] = amount;
            self(self, coordinate + 1, selected_label + label * amount);
        }
        selected[static_cast<std::size_t>(coordinate)] = 0;
    };
    visit(visit, 0, 0);
    cache.emplace(key, result);
    return result;
}

Counts signs_for_pattern(int pattern) {
    Counts signs{};
    signs[0U] = 1;
    for (int index = 1; index < label_count; ++index) {
        const int bit = 1 << (label_count - 1 - index);
        signs[static_cast<std::size_t>(index)] = (pattern & bit) == 0 ? 1 : -1;
    }
    return signs;
}

void print_counts(const Counts& counts) {
    std::cout << '[';
    for (std::size_t index = 0U; index < counts.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << counts[index];
    }
    std::cout << ']';
}

cpp_int finite_difference(
    const Counts& base,
    const Counts& degree,
    const Counts& signs
) {
    cpp_int result = 0;
    Counts offset{};
    int total_degree = 0;
    for (int value : degree) {
        total_degree += value;
    }
    const auto visit = [&](const auto& self, int coordinate,
                           int used_degree, cpp_int coefficient) -> void {
        if (coordinate == label_count) {
            Counts point{};
            for (std::size_t index = 0U; index < point.size(); ++index) {
                point[index] = base[index] + offset[index];
            }
            const cpp_int term = coefficient * outer_value(point, signs);
            result += ((total_degree - used_degree) & 1) == 0
                ? term : -term;
            return;
        }
        const int maximum = degree[static_cast<std::size_t>(coordinate)];
        for (int value = 0; value <= maximum; ++value) {
            offset[static_cast<std::size_t>(coordinate)] = value;
            self(
                self,
                coordinate + 1,
                used_degree + value,
                coefficient * binomial(maximum, value)
            );
        }
    };
    visit(visit, 0, 0, 1);
    return result;
}

bool first_negative_newton_coefficient(
    const Counts& base,
    const Counts& signs,
    Counts& negative_degree,
    cpp_int& negative_value,
    std::uint64_t& tested
) {
    Counts degree{};
    for (int total_degree = 0; total_degree <= deficit; ++total_degree) {
        const auto visit = [&](const auto& self, int coordinate,
                               int remaining) -> bool {
            if (coordinate == label_count) {
                if (remaining != 0) {
                    return false;
                }
                const cpp_int coefficient = finite_difference(
                    base, degree, signs
                );
                ++tested;
                if (coefficient < 0) {
                    negative_degree = degree;
                    negative_value = coefficient;
                    return true;
                }
                return false;
            }
            for (int value = 0; value <= remaining; ++value) {
                degree[static_cast<std::size_t>(coordinate)] = value;
                if (self(self, coordinate + 1, remaining - value)) {
                    return true;
                }
            }
            return false;
        };
        if (visit(visit, 0, total_degree)) {
            return true;
        }
    }
    return false;
}

std::vector<Counts> minimal_bases() {
    std::vector<Counts> bases;
    Counts counts{};
    const auto visit = [&](const auto& self, int coordinate,
                           int used_label) -> void {
        if (coordinate == label_count) {
            if (!admissible(counts)) {
                return;
            }
            for (int index = 0; index < label_count; ++index) {
                const std::size_t position = static_cast<std::size_t>(index);
                if (counts[position] == 0) {
                    continue;
                }
                --counts[position];
                const bool predecessor_admissible = admissible(counts);
                ++counts[position];
                if (predecessor_admissible) {
                    return;
                }
            }
            bases.push_back(counts);
            return;
        }
        const int label = coordinate + 1;
        for (int amount = 0; used_label + label * amount <= 24; ++amount) {
            counts[static_cast<std::size_t>(coordinate)] = amount;
            self(self, coordinate + 1, used_label + label * amount);
        }
        counts[static_cast<std::size_t>(coordinate)] = 0;
    };
    visit(visit, 0, 0);
    return bases;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 5) {
            throw std::runtime_error(
                "usage: probe_su2_outer_twelve MAXIMUM_FACTORS "
                "[minimal-bases|newton PATTERN [BASE_INDEX]]"
            );
        }
        const int maximum_factors = std::stoi(argv[1]);
        if (maximum_factors < 1 || maximum_factors > 60) {
            throw std::runtime_error("invalid maximum factor count");
        }
        const std::string mode = argc >= 3 ? argv[2] : "search";
        if (mode != "search" && mode != "minimal-bases" && mode != "newton") {
            throw std::runtime_error("invalid probe mode");
        }
        if ((mode == "search" && argc != 2)
            || (mode == "minimal-bases" && argc != 3)
            || (mode == "newton" && argc != 4 && argc != 5)) {
            throw std::runtime_error("invalid probe arguments");
        }

        if (mode == "minimal-bases") {
            const std::vector<Counts> bases = minimal_bases();
            int largest_total_label = 0;
            for (const Counts& base : bases) {
                largest_total_label = std::max(
                    largest_total_label, total_label(base)
                );
            }
            std::cout << "SU2_OUTER_TWELVE_MINIMAL_BASES bases=" << bases.size()
                      << " maximum_total_label=" << largest_total_label
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }

        if (mode == "newton") {
            const int pattern = std::stoi(argv[3]);
            if (pattern < 0 || pattern >= (1 << (label_count - 1))) {
                throw std::runtime_error("invalid Newton sign pattern");
            }
            const std::vector<Counts> bases = minimal_bases();
            const int requested_base = argc == 5 ? std::stoi(argv[4]) : -1;
            if (requested_base < -1
                || requested_base >= static_cast<int>(bases.size())) {
                throw std::runtime_error("invalid Newton base index");
            }
            const Counts signs = signs_for_pattern(pattern);
            const int first_base = requested_base < 0 ? 0 : requested_base;
            const int last_base = requested_base < 0
                ? static_cast<int>(bases.size()) - 1 : requested_base;
            std::uint64_t tested = 0U;
            for (int index = first_base; index <= last_base; ++index) {
                Counts negative_degree{};
                cpp_int negative_value = 0;
                if (first_negative_newton_coefficient(
                        bases[static_cast<std::size_t>(index)],
                        signs,
                        negative_degree,
                        negative_value,
                        tested
                    )) {
                    std::cout << "SU2_OUTER_TWELVE_NEWTON_FAIL pattern="
                              << pattern << " base_index=" << index
                              << " base=";
                    print_counts(bases[static_cast<std::size_t>(index)]);
                    std::cout << " degree=";
                    print_counts(negative_degree);
                    std::cout << " value=" << negative_value
                              << " tested=" << tested << '\n';
                    return EXIT_FAILURE;
                }
            }
            std::cout << "SU2_OUTER_TWELVE_NEWTON_PASS pattern=" << pattern
                      << " first_base=" << first_base
                      << " last_base=" << last_base
                      << " bases=" << last_base - first_base + 1
                      << " coefficients=" << tested << '\n';
            return EXIT_SUCCESS;
        }

        std::uint64_t tested = 0U;
        cpp_int minimum = 0;
        Counts minimum_counts{};
        int minimum_pattern = 0;
        bool initialized = false;
        Counts counts{};
        const auto visit = [&](const auto& self, int coordinate) -> bool {
            if (coordinate == label_count) {
                if (!admissible(counts)) {
                    return true;
                }
                for (int pattern = 0; pattern < (1 << (label_count - 1));
                     ++pattern) {
                    const cpp_int value = outer_value(
                        counts, signs_for_pattern(pattern)
                    );
                    ++tested;
                    if (!initialized || value < minimum) {
                        initialized = true;
                        minimum = value;
                        minimum_counts = counts;
                        minimum_pattern = pattern;
                    }
                    if (value < 0) {
                        std::cout << "SU2_OUTER_TWELVE result=FAIL pattern="
                                  << pattern << " counts=";
                        print_counts(counts);
                        std::cout << " value=" << value << '\n';
                        return false;
                    }
                }
                return true;
            }

            const int used = total_count(counts);
            for (int amount = 0; amount + used <= maximum_factors; ++amount) {
                counts[static_cast<std::size_t>(coordinate)] = amount;
                if (!self(self, coordinate + 1)) {
                    return false;
                }
            }
            counts[static_cast<std::size_t>(coordinate)] = 0;
            return true;
        };

        if (!visit(visit, 0)) {
            return EXIT_FAILURE;
        }
        std::cout << "SU2_OUTER_TWELVE result=PASS tested=" << tested
                  << " maximum_factors=" << maximum_factors
                  << " minimum=" << minimum
                  << " minimum_pattern=" << minimum_pattern
                  << " minimum_counts=";
        print_counts(minimum_counts);
        std::cout << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
