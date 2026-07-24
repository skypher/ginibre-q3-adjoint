#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

struct FusionInterval {
    int lower;
    int upper;
    int parity;
};

FusionInterval fusion_interval(
    const int first,
    const int second,
    const int level
) {
    return {
        std::abs(first - second),
        std::min(first + second, 2 * level - first - second),
        (first + second) & 1
    };
}

int intersection_size(
    const FusionInterval& first,
    const FusionInterval& second
) {
    if (first.parity != second.parity) {
        return 0;
    }
    const int lower = std::max(first.lower, second.lower);
    const int upper = std::min(first.upper, second.upper);
    return lower <= upper ? (upper - lower) / 2 + 1 : 0;
}

int triple_multiplicity(
    const int first,
    const int second,
    const int third,
    const int target,
    const int level
) {
    return intersection_size(
        fusion_interval(first, second, level),
        fusion_interval(third, target, level)
    );
}

bool invariant_triple(
    const std::array<int, 3>& labels,
    const int level
) {
    return triple_multiplicity(
        labels[0], labels[1], labels[2], 0, level
    ) != 0;
}

cpp_int sixfold_invariant(
    const std::array<int, 6>& labels,
    const int level
) {
    cpp_int result = 0;
    for (int target = 0; target <= level; ++target) {
        const int left = triple_multiplicity(
            labels[0], labels[1], labels[2], target, level
        );
        const int right = triple_multiplicity(
            labels[3], labels[4], labels[5], target, level
        );
        result += cpp_int{left} * right;
    }
    return result;
}

int admissible_split_count(
    const std::array<int, 6>& labels,
    const int level
) {
    constexpr unsigned int full_mask = (1U << 6U) - 1U;
    int result = 0;
    for (unsigned int mask = 0; mask <= full_mask; ++mask) {
        if (std::popcount(mask) != 3) {
            continue;
        }
        const unsigned int complement = full_mask ^ mask;
        if (mask > complement) {
            continue;
        }
        std::array<int, 3> left{};
        std::array<int, 3> right{};
        std::size_t left_index = 0;
        std::size_t right_index = 0;
        for (std::size_t index = 0; index < labels.size(); ++index) {
            if ((mask & (1U << index)) != 0U) {
                left[left_index++] = labels[index];
            } else {
                right[right_index++] = labels[index];
            }
        }
        if (invariant_triple(left, level)
            && invariant_triple(right, level)) {
            ++result;
        }
    }
    return result;
}

cpp_int binomial(int n, int r) {
    if (r < 0 || r > n) {
        return 0;
    }
    r = std::min(r, n - r);
    cpp_int result = 1;
    for (int index = 1; index <= r; ++index) {
        result *= n - r + index;
        result /= index;
    }
    return result;
}

int parse_maximum_level(const int argc, char** argv) {
    if (argc > 2) {
        throw std::runtime_error(
            "usage: verify_su2_finite_distinct_six [MAX_K]"
        );
    }
    if (argc == 1) {
        return 30;
    }
    std::size_t parsed = 0;
    const std::string argument = argv[1];
    const int value = std::stoi(argument, &parsed);
    if (parsed != argument.size() || value < 1) {
        throw std::runtime_error("MAX_K must be a positive integer");
    }
    return value;
}

void print_labels(const std::array<int, 6>& labels) {
    std::cerr << '[';
    for (std::size_t index = 0; index < labels.size(); ++index) {
        if (index != 0) {
            std::cerr << ',';
        }
        std::cerr << labels[index];
    }
    std::cerr << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = parse_maximum_level(argc, argv);
        std::uint64_t cases = 0;
        std::uint64_t active_cases = 0;
        cpp_int minimum_margin = 0;
        bool have_active_case = false;

        for (int level = 1; level <= maximum_level; ++level) {
            for (int a = 1; a <= level; ++a) {
                for (int b = a + 1; b <= level; ++b) {
                    for (int c = b + 1; c <= level; ++c) {
                        for (int d = c + 1; d <= level; ++d) {
                            for (int e = d + 1; e <= level; ++e) {
                                for (int f = e + 1; f <= level; ++f) {
                                    const std::array<int, 6> labels{
                                        a, b, c, d, e, f
                                    };
                                    ++cases;
                                    const int split_count =
                                        admissible_split_count(labels, level);
                                    const cpp_int invariant =
                                        sixfold_invariant(labels, level);
                                    if (invariant < split_count) {
                                        std::cerr
                                            << "SU2_FINITE_DISTINCT_SIX FAIL"
                                            << " k=" << level
                                            << " labels=";
                                        print_labels(labels);
                                        std::cerr
                                            << " invariant=" << invariant
                                            << " splits=" << split_count
                                            << '\n';
                                        return EXIT_FAILURE;
                                    }
                                    if (split_count != 0) {
                                        ++active_cases;
                                        const cpp_int margin =
                                            invariant - split_count;
                                        if (!have_active_case
                                            || margin < minimum_margin) {
                                            minimum_margin = margin;
                                            have_active_case = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        const cpp_int expected_cases = binomial(maximum_level + 1, 7);
        if (cpp_int{cases} != expected_cases) {
            throw std::runtime_error(
                "distinct-sextuple enumeration count mismatch"
            );
        }
        std::cout << "SU2_FINITE_DISTINCT_SIX PASS"
                  << " max_k=" << maximum_level
                  << " cases=" << cases
                  << " active_cases=" << active_cases
                  << " minimum_margin="
                  << (have_active_case ? minimum_margin : cpp_int{0})
                  << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_FINITE_DISTINCT_SIX ERROR "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
