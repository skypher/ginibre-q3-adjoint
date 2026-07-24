#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

using Laurent = std::map<int, cpp_int>;

int parse_positive(const char* text, const char* name) {
    const int value = std::stoi(text);
    if (value <= 0) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return value;
}

Laurent multiply(const Laurent& left, const Laurent& right) {
    Laurent result;
    for (const auto& [left_exponent, left_coefficient] : left) {
        for (const auto& [right_exponent, right_coefficient] : right) {
            result[left_exponent + right_exponent]
                += left_coefficient * right_coefficient;
        }
    }
    return result;
}

Laurent power(Laurent base, int exponent) {
    Laurent result{{0, 1}};
    while (exponent > 0) {
        if ((exponent & 1) != 0) {
            result = multiply(result, base);
        }
        exponent /= 2;
        if (exponent > 0) {
            base = multiply(base, base);
        }
    }
    return result;
}

cpp_int moment(
    int order,
    int sector,
    int minus_pairs,
    int plus_pairs,
    int degree
) {
    const Laurent two_minus_c{{-1, -1}, {0, 2}, {1, -1}};
    const Laurent two_plus_c{{-1, 1}, {0, 2}, {1, 1}};
    const Laurent c{{-1, 1}, {1, 1}};
    const Laurent polynomial = multiply(
        multiply(
            power(two_minus_c, minus_pairs),
            power(two_plus_c, plus_pairs)
        ),
        power(c, degree)
    );

    cpp_int answer = 0;
    for (const auto& [exponent, coefficient] : polynomial) {
        if (exponent % order != 0) {
            continue;
        }
        cpp_int root_sum = order;
        if (sector < 0 && ((exponent / order) & 1) != 0) {
            root_sum = -root_sum;
        }
        answer += coefficient * root_sum;
    }
    return answer;
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

cpp_int sector_sum(
    int order,
    int sector,
    int minus_pairs,
    int plus_pairs,
    int label_two_plus_count
) {
    std::vector<cpp_int> moments(
        static_cast<std::size_t>(label_two_plus_count + 3)
    );
    for (int degree = 0; degree <= label_two_plus_count + 2; ++degree) {
        moments[static_cast<std::size_t>(degree)] = moment(
            order, sector, minus_pairs, plus_pairs, degree
        );
    }
    cpp_int answer = 0;
    cpp_int power_of_two = cpp_int(1) << label_two_plus_count;
    for (int degree = 0; degree <= label_two_plus_count; ++degree) {
        const cpp_int determinant
            = moments[static_cast<std::size_t>(degree)]
                * moments[static_cast<std::size_t>(degree + 2)]
            - moments[static_cast<std::size_t>(degree + 1)]
                * moments[static_cast<std::size_t>(degree + 1)];
        answer += 2 * binomial(label_two_plus_count, degree)
            * power_of_two * determinant;
        if (degree != label_two_plus_count) {
            power_of_two /= 2;
        }
    }
    return answer;
}

template <typename Function>
void for_each_fusion_output(
    int level,
    int left,
    int right,
    Function&& function
) {
    const int lower = std::abs(left - right);
    const int upper = std::min(left + right, 2 * level - left - right);
    for (int output = lower; output <= upper; output += 2) {
        function(output);
    }
}

using Matrix = std::vector<cpp_int>;

Matrix apply_factor(
    const Matrix& state,
    int level,
    int label,
    int sign
) {
    const int rank = level + 1;
    Matrix next(static_cast<std::size_t>(rank * rank));
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            const cpp_int& coefficient = state[
                static_cast<std::size_t>(left * rank + right)
            ];
            if (coefficient == 0) {
                continue;
            }
            for_each_fusion_output(
                level,
                left,
                label,
                [&](int output) {
                    next[static_cast<std::size_t>(output * rank + right)]
                        += coefficient;
                }
            );
            for_each_fusion_output(
                level,
                right,
                label,
                [&](int output) {
                    next[static_cast<std::size_t>(left * rank + output)]
                        += sign * coefficient;
                }
            );
        }
    }
    return next;
}

cpp_int direct_value(
    int level,
    int minus_pairs,
    int plus_pairs,
    int label_two_plus_count
) {
    const int rank = level + 1;
    Matrix state(static_cast<std::size_t>(rank * rank));
    state[0] = 1;
    for (int index = 0; index < 2 * minus_pairs; ++index) {
        state = apply_factor(state, level, 1, -1);
    }
    for (int index = 0; index < 2 * plus_pairs; ++index) {
        state = apply_factor(state, level, 1, 1);
    }
    for (int index = 0; index < label_two_plus_count; ++index) {
        state = apply_factor(state, level, 2, 1);
    }
    return state[0];
}

Matrix canonical_state(
    int level,
    int minus_pairs,
    int plus_pairs,
    int label_two_plus_count
) {
    const int rank = level + 1;
    Matrix state(static_cast<std::size_t>(rank * rank));
    state[0] = 1;
    for (int index = 0; index < 2 * minus_pairs; ++index) {
        state = apply_factor(state, level, 1, -1);
    }
    for (int index = 0; index < 2 * plus_pairs; ++index) {
        state = apply_factor(state, level, 1, 1);
    }
    for (int index = 0; index < label_two_plus_count; ++index) {
        state = apply_factor(state, level, 2, 1);
    }
    return state;
}

struct PrefixWitness {
    bool negative = false;
    int parity = 0;
    int maximum_left = 0;
    int maximum_right = 0;
    cpp_int value = 0;
};

PrefixWitness first_negative_parity_prefix_rectangle(
    const Matrix& state,
    int level
) {
    const int rank = level + 1;
    for (int parity = 0; parity <= 1; ++parity) {
        for (int maximum_left = parity;
             maximum_left <= level;
             maximum_left += 2) {
            for (int maximum_right = parity;
                 maximum_right <= level;
                 maximum_right += 2) {
                cpp_int sum = 0;
                for (int left = parity;
                     left <= maximum_left;
                     left += 2) {
                    for (int right = parity;
                         right <= maximum_right;
                         right += 2) {
                        sum += state[static_cast<std::size_t>(
                            left * rank + right
                        )];
                    }
                }
                if (sum < 0) {
                    return PrefixWitness{
                        true,
                        parity,
                        maximum_left,
                        maximum_right,
                        sum
                    };
                }
            }
        }
    }
    return PrefixWitness{};
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5) {
            throw std::runtime_error(
                "usage: analyze_su2_two_generator_cyclic "
                "MAX_LEVEL MAX_M MAX_N MAX_C"
            );
        }
        const int maximum_level = parse_positive(argv[1], "maximum level");
        const int maximum_m = parse_positive(argv[2], "maximum m");
        const int maximum_n = parse_positive(argv[3], "maximum n");
        const int maximum_c = parse_positive(argv[4], "maximum c");

        std::uint64_t cases = 0U;
        std::uint64_t negative_sectors = 0U;
        std::uint64_t prefix_cone_failures = 0U;
        bool printed_negative_sector = false;
        bool printed_prefix_failure = false;
        bool printed_residual_prefix_failure = false;
        for (int level = 2; level <= maximum_level; ++level) {
            const int order = level + 2;
            for (int minus_pairs = 1;
                 minus_pairs <= maximum_m;
                 ++minus_pairs) {
                for (int plus_pairs = 1;
                     plus_pairs <= maximum_n;
                     ++plus_pairs) {
                    for (int label_two_plus_count = 0;
                         label_two_plus_count <= maximum_c;
                         ++label_two_plus_count) {
                        const cpp_int plus_sector = sector_sum(
                            order,
                            1,
                            minus_pairs,
                            plus_pairs,
                            label_two_plus_count
                        );
                        const cpp_int minus_sector = sector_sum(
                            order,
                            -1,
                            minus_pairs,
                            plus_pairs,
                            label_two_plus_count
                        );
                        const cpp_int numerator
                            = plus_sector + minus_sector;
                        const cpp_int denominator = 8 * order * order;
                        const cpp_int direct = direct_value(
                            level,
                            minus_pairs,
                            plus_pairs,
                            label_two_plus_count
                        );
                        if (numerator % denominator != 0
                            || numerator / denominator != direct) {
                            std::cerr
                                << "FAIL identity level=" << level
                                << " m=" << minus_pairs
                                << " n=" << plus_pairs
                                << " c=" << label_two_plus_count
                                << " direct=" << direct
                                << " numerator=" << numerator
                                << " denominator=" << denominator << '\n';
                            return EXIT_FAILURE;
                        }
                        if (direct < 0) {
                            std::cerr
                                << "FAIL negative total level=" << level
                                << " m=" << minus_pairs
                                << " n=" << plus_pairs
                                << " c=" << label_two_plus_count
                                << " value=" << direct << '\n';
                            return EXIT_FAILURE;
                        }
                        const Matrix canonical = canonical_state(
                            level,
                            minus_pairs,
                            plus_pairs,
                            label_two_plus_count
                        );
                        const PrefixWitness prefix_witness
                            = first_negative_parity_prefix_rectangle(
                                canonical, level
                            );
                        if (prefix_witness.negative) {
                            ++prefix_cone_failures;
                            if (!printed_prefix_failure) {
                                std::cout
                                    << "FIRST_PREFIX_CONE_FAILURE"
                                    << " level=" << level
                                    << " m=" << minus_pairs
                                    << " n=" << plus_pairs
                                    << " c=" << label_two_plus_count
                                    << " parity=" << prefix_witness.parity
                                    << " max_left="
                                    << prefix_witness.maximum_left
                                    << " max_right="
                                    << prefix_witness.maximum_right
                                    << " value=" << prefix_witness.value
                                    << '\n';
                                printed_prefix_failure = true;
                            }
                            if (minus_pairs >= 2
                                && plus_pairs >= 2
                                && !printed_residual_prefix_failure) {
                                std::cout
                                    << "FIRST_RESIDUAL_PREFIX_CONE_FAILURE"
                                    << " level=" << level
                                    << " m=" << minus_pairs
                                    << " n=" << plus_pairs
                                    << " c=" << label_two_plus_count
                                    << " parity=" << prefix_witness.parity
                                    << " max_left="
                                    << prefix_witness.maximum_left
                                    << " max_right="
                                    << prefix_witness.maximum_right
                                    << " value=" << prefix_witness.value
                                    << '\n';
                                printed_residual_prefix_failure = true;
                            }
                        }
                        for (const auto& [name, value] : {
                                 std::pair<const char*, cpp_int>{
                                     "plus", plus_sector
                                 },
                                 std::pair<const char*, cpp_int>{
                                     "minus", minus_sector
                                 }
                             }) {
                            if (value < 0) {
                                ++negative_sectors;
                                if (!printed_negative_sector) {
                                    std::cout
                                        << "FIRST_NEGATIVE_SECTOR"
                                        << " level=" << level
                                        << " m=" << minus_pairs
                                        << " n=" << plus_pairs
                                        << " c=" << label_two_plus_count
                                        << " sector=" << name
                                        << " value=" << value << '\n';
                                    printed_negative_sector = true;
                                }
                            }
                        }
                        ++cases;
                    }
                }
            }
        }
        std::cout
            << "PASS max_level=" << maximum_level
            << " max_m=" << maximum_m
            << " max_n=" << maximum_n
            << " max_c=" << maximum_c
            << " cases=" << cases
            << " negative_sectors=" << negative_sectors
            << " prefix_cone_failures=" << prefix_cone_failures << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
