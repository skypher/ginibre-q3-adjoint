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
    for (const auto& [first_exponent, first_coefficient] : left) {
        for (const auto& [second_exponent, second_coefficient] : right) {
            result[first_exponent + second_exponent]
                += first_coefficient * second_coefficient;
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

cpp_int cyclic_moment(int order, int sector, int minus_pairs, int degree) {
    const Laurent four_minus_c_squared{
        {-2, -1},
        {0, 2},
        {2, -1}
    };
    const Laurent c{{-1, 1}, {1, 1}};
    const Laurent polynomial = multiply(
        power(four_minus_c_squared, minus_pairs),
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

cpp_int sector_sum(int order, int sector, int minus_pairs, int plus_count) {
    std::vector<cpp_int> moments(
        static_cast<std::size_t>(plus_count + 3)
    );
    for (int degree = 0; degree <= plus_count + 2; ++degree) {
        moments[static_cast<std::size_t>(degree)]
            = cyclic_moment(order, sector, minus_pairs, degree);
    }

    cpp_int answer = 0;
    cpp_int power_of_two = cpp_int(1) << plus_count;
    for (int degree = 0; degree <= plus_count; ++degree) {
        const cpp_int determinant
            = moments[static_cast<std::size_t>(degree)]
                * moments[static_cast<std::size_t>(degree + 2)]
            - moments[static_cast<std::size_t>(degree + 1)]
                * moments[static_cast<std::size_t>(degree + 1)];
        answer += 2 * binomial(plus_count, degree)
            * power_of_two * determinant;
        if (degree != plus_count) {
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

cpp_int direct_ray(
    int level,
    int minus_pairs,
    int plus_count
) {
    const int rank = level + 1;
    Matrix state(static_cast<std::size_t>(rank * rank));
    state[0] = 1;
    for (int index = 0; index < 2 * minus_pairs; ++index) {
        state = apply_factor(state, level, 2, -1);
    }
    for (int index = 0; index < plus_count; ++index) {
        state = apply_factor(state, level, 2, 1);
    }
    return state[0];
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4) {
            throw std::runtime_error(
                "usage: verify_su2_label2_ray_cyclic "
                "MAX_LEVEL MAX_MINUS_PAIRS MAX_PLUS_COUNT"
            );
        }
        const int maximum_level = parse_positive(argv[1], "maximum level");
        const int maximum_minus_pairs
            = parse_positive(argv[2], "maximum minus pairs");
        const int maximum_plus_count
            = parse_positive(argv[3], "maximum plus count");

        std::uint64_t cases = 0U;
        std::uint64_t negative_individual_minors = 0U;
        cpp_int minimum_value = 0;
        bool minimum_initialized = false;

        for (int level = 2; level <= maximum_level; ++level) {
            const int order = level + 2;
            for (int minus_pairs = 1;
                 minus_pairs <= maximum_minus_pairs;
                 ++minus_pairs) {
                for (int plus_count = 0;
                     plus_count <= maximum_plus_count;
                     ++plus_count) {
                    const cpp_int direct
                        = direct_ray(level, minus_pairs, plus_count);
                    cpp_int numerator = sector_sum(
                        order, 1, minus_pairs, plus_count
                    );
                    cpp_int denominator = 4 * order * order;
                    if ((order & 1) == 0) {
                        numerator += sector_sum(
                            order, -1, minus_pairs, plus_count
                        );
                        denominator *= 2;
                    }
                    if (numerator % denominator != 0
                        || numerator / denominator != direct) {
                        std::cerr
                            << "FAIL identity level=" << level
                            << " minus_pairs=" << minus_pairs
                            << " plus_count=" << plus_count
                            << " direct=" << direct
                            << " numerator=" << numerator
                            << " denominator=" << denominator << '\n';
                        return EXIT_FAILURE;
                    }
                    if (direct < 0) {
                        std::cerr
                            << "FAIL negative level=" << level
                            << " minus_pairs=" << minus_pairs
                            << " plus_count=" << plus_count
                            << " value=" << direct << '\n';
                        return EXIT_FAILURE;
                    }
                    if (!minimum_initialized || direct < minimum_value) {
                        minimum_initialized = true;
                        minimum_value = direct;
                    }
                    ++cases;
                }

                for (int sector : {1, -1}) {
                    if ((order & 1) != 0 && sector < 0) {
                        continue;
                    }
                    for (int degree = 0;
                         degree <= maximum_plus_count;
                         ++degree) {
                        const cpp_int first = cyclic_moment(
                            order, sector, minus_pairs, degree
                        );
                        const cpp_int middle = cyclic_moment(
                            order, sector, minus_pairs, degree + 1
                        );
                        const cpp_int last = cyclic_moment(
                            order, sector, minus_pairs, degree + 2
                        );
                        if (first * last - middle * middle < 0) {
                            ++negative_individual_minors;
                        }
                    }
                }
            }
        }

        std::cout
            << "PASS max_level=" << maximum_level
            << " max_minus_pairs=" << maximum_minus_pairs
            << " max_plus_count=" << maximum_plus_count
            << " cases=" << cases
            << " negative_individual_minors="
            << negative_individual_minors
            << " minimum_value=" << minimum_value << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
