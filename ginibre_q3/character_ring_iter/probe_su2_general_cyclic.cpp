#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

using Laurent = std::map<int, cpp_int>;

Laurent multiply(const Laurent& left, const Laurent& right) {
    Laurent answer;
    for (const auto& [i, a] : left) {
        for (const auto& [j, b] : right) {
            answer[i + j] += a * b;
        }
    }
    return answer;
}

Laurent power(Laurent base, int exponent) {
    Laurent answer{{0, 1}};
    while (exponent > 0) {
        if ((exponent & 1) != 0) {
            answer = multiply(answer, base);
        }
        exponent /= 2;
        if (exponent > 0) {
            base = multiply(base, base);
        }
    }
    return answer;
}

cpp_int binomial(int n, int r) {
    r = std::min(r, n - r);
    cpp_int answer = 1;
    for (int j = 1; j <= r; ++j) {
        answer *= n - r + j;
        answer /= j;
    }
    return answer;
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

cpp_int sector_sum(
    int order,
    int sector,
    int minus_pairs,
    int plus_pairs,
    int label_two_plus
) {
    std::vector<cpp_int> moments(
        static_cast<std::size_t>(label_two_plus + 3)
    );
    for (int degree = 0; degree <= label_two_plus + 2; ++degree) {
        moments[static_cast<std::size_t>(degree)] = moment(
            order, sector, minus_pairs, plus_pairs, degree
        );
    }
    cpp_int answer = 0;
    cpp_int power_of_two = cpp_int(1) << label_two_plus;
    for (int degree = 0; degree <= label_two_plus; ++degree) {
        const cpp_int determinant
            = moments[static_cast<std::size_t>(degree)]
                * moments[static_cast<std::size_t>(degree + 2)]
            - moments[static_cast<std::size_t>(degree + 1)]
                * moments[static_cast<std::size_t>(degree + 1)];
        answer += 2 * binomial(label_two_plus, degree)
            * power_of_two * determinant;
        if (degree != label_two_plus) {
            power_of_two /= 2;
        }
    }
    return answer;
}

long double spectral_transport_margin(
    int order,
    int sector,
    int minus_pairs,
    int plus_pairs
) {
    long double positive = 0.0L;
    long double negative = 0.0L;
    for (int i = 0; i < order; ++i) {
        const long double theta_i
            = (2.0L * static_cast<long double>(i)
                + (sector < 0 ? 1.0L : 0.0L))
            * std::numbers::pi_v<long double>
            / static_cast<long double>(order);
        const long double first = 2.0L * std::cos(theta_i);
        const long double first_weight = std::pow(
            2.0L - first, minus_pairs
        ) * std::pow(2.0L + first, plus_pairs);
        for (int j = 0; j < order; ++j) {
            const long double theta_j
                = (2.0L * static_cast<long double>(j)
                    + (sector < 0 ? 1.0L : 0.0L))
                * std::numbers::pi_v<long double>
                / static_cast<long double>(order);
            const long double second = 2.0L * std::cos(theta_j);
            const long double second_weight = std::pow(
                2.0L - second, minus_pairs
            ) * std::pow(2.0L + second, plus_pairs);
            const long double base
                = (first - second) * (first - second)
                * first_weight * second_weight;
            const long double eigenvalue = 2.0L + first * second;
            if (first * second >= 0.0L) {
                positive += base * eigenvalue;
            } else if (eigenvalue < 0.0L) {
                negative += base * (-eigenvalue);
            }
        }
    }
    return positive - negative;
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
    int label_two_plus
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
    for (int index = 0; index < label_two_plus; ++index) {
        state = apply_factor(state, level, 2, 1);
    }
    return state[0];
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5) {
            throw std::runtime_error(
                "usage: probe_su2_general_cyclic "
                "MAX_LEVEL MAX_M MAX_N MAX_C"
            );
        }
        const int maximum_level = std::stoi(argv[1]);
        const int maximum_m = std::stoi(argv[2]);
        const int maximum_n = std::stoi(argv[3]);
        const int maximum_c = std::stoi(argv[4]);
        std::uint64_t cases = 0;
        std::uint64_t negative_sectors = 0;
        std::uint64_t first_alias_cases = 0;
        std::uint64_t negative_first_alias_corrections = 0;
        bool reported_hard_first_alias_correction = false;
        std::uint64_t reported_transport_failures = 0;
        bool reported_variance_failure = false;
        bool reported_jacobi_sign_failure = false;
        long double minimum_transport_margin
            = std::numeric_limits<long double>::infinity();
        long double minimum_odd_rounded_ratio
            = std::numeric_limits<long double>::infinity();
        for (int level = 2; level <= maximum_level; ++level) {
            const int order = level + 2;
            if ((order & 1) != 0) {
                const long double delta
                    = std::numbers::pi_v<long double>
                    / (2.0L * static_cast<long double>(order));
                for (int first = 1; 2 * first < order; ++first) {
                    const long double alpha
                        = static_cast<long double>(first)
                        * std::numbers::pi_v<long double>
                        / static_cast<long double>(order);
                    for (int second = 1; 2 * second < order; ++second) {
                        const long double beta
                            = static_cast<long double>(second)
                            * std::numbers::pi_v<long double>
                            / static_cast<long double>(order);
                        if (2.0L * std::cos(alpha) * std::cos(beta)
                            <= 1.0L) {
                            continue;
                        }
                        const long double ratio
                            = std::cos(alpha - delta)
                            * std::cos(beta - delta)
                            * std::sin(
                                (alpha + beta) / 2.0L - delta
                            )
                            / (
                                std::sin(alpha) * std::sin(beta)
                                * std::cos((alpha + beta) / 2.0L)
                            );
                        minimum_odd_rounded_ratio = std::min(
                            minimum_odd_rounded_ratio, ratio
                        );
                    }
                }
            }
            for (int m = 0; m <= maximum_m; ++m) {
                for (int n = 0; n <= maximum_n; ++n) {
                    for (int sector : {1, -1}) {
                        const cpp_int a0 = moment(
                            order, sector, m, n, 0
                        );
                        const cpp_int a1 = moment(
                            order, sector, m, n, 1
                        );
                        const cpp_int a2 = moment(
                            order, sector, m, n, 2
                        );
                        const cpp_int a3 = moment(
                            order, sector, m, n, 3
                        );
                        const cpp_int determinant = a0 * a2 - a1 * a1;
                        if (determinant > 2 * a0 * a0
                            && m >= 1 && n >= 1
                            && !reported_variance_failure) {
                            std::cout
                                << "VARIANCE_FAILURE level=" << level
                                << " sector=" << sector
                                << " m=" << m << " n=" << n << '\n';
                            reported_variance_failure = true;
                        }
                        const cpp_int b1_numerator
                            = a0 * a0 * a3
                            - 2 * a0 * a1 * a2
                            + a1 * a1 * a1;
                        if (a1 * b1_numerator < 0
                            && m >= 2 && n >= 2
                            && level > 2
                            && !reported_jacobi_sign_failure) {
                            std::cout
                                << "JACOBI_SIGN_FAILURE level=" << level
                                << " sector=" << sector
                                << " m=" << m << " n=" << n << '\n';
                            reported_jacobi_sign_failure = true;
                        }
                        const long double margin = spectral_transport_margin(
                            order, sector, m, n
                        );
                        minimum_transport_margin = std::min(
                            minimum_transport_margin, margin
                        );
                        if (margin < -1.0e-7L
                            && m >= 2 && n >= 2
                            && reported_transport_failures < 500) {
                            std::cout
                                << "TRANSPORT_FAILURE level=" << level
                                << " sector=" << sector
                                << " m=" << m << " n=" << n
                                << " margin="
                                << static_cast<double>(margin) << '\n';
                            ++reported_transport_failures;
                        }
                    }
                    for (int c = 0; c <= maximum_c; ++c) {
                        const cpp_int plus = sector_sum(
                            order, 1, m, n, c
                        );
                        const cpp_int minus = sector_sum(
                            order, -1, m, n, c
                        );
                        if (plus < 0 || minus < 0) {
                            ++negative_sectors;
                            if (negative_sectors <= 20) {
                                std::cout
                                    << "NEGATIVE_SECTOR level=" << level
                                    << " m=" << m << " n=" << n
                                    << " c=" << c << " plus=" << plus
                                    << " minus=" << minus << '\n';
                            }
                        }
                        const cpp_int numerator = plus + minus;
                        const cpp_int denominator
                            = 8 * order * order;
                        const cpp_int direct = direct_value(
                            level, m, n, c
                        );
                        if (numerator % denominator != 0
                            || numerator / denominator != direct) {
                            std::cout
                                << "IDENTITY_FAILURE level=" << level
                                << " m=" << m << " n=" << n
                                << " c=" << c
                                << " numerator=" << numerator
                                << " denominator=" << denominator
                                << " direct=" << direct << '\n';
                            return EXIT_FAILURE;
                        }
                        if (numerator < 0) {
                            std::cout
                                << "NEGATIVE_TOTAL level=" << level
                                << " m=" << m << " n=" << n
                                << " c=" << c << " numerator="
                                << numerator << '\n';
                            return EXIT_FAILURE;
                        }
                        const int degree = m + n + c + 2;
                        if (order <= degree && degree < 2 * order) {
                            const int no_alias_order = degree + 1;
                            const cpp_int no_alias_sector = sector_sum(
                                no_alias_order, 1, m, n, c
                            );
                            if (no_alias_sector
                                % (no_alias_order * no_alias_order) != 0) {
                                throw std::runtime_error(
                                    "nonintegral constant coefficient"
                                );
                            }
                            const cpp_int constant_coefficient
                                = no_alias_sector
                                / (no_alias_order * no_alias_order);
                            const cpp_int correction = numerator
                                - 2 * order * order * constant_coefficient;
                            if (correction % (8 * order * order) != 0) {
                                throw std::runtime_error(
                                    "nonintegral first-alias coefficient"
                                );
                            }
                            if (correction < 0) {
                                ++negative_first_alias_corrections;
                                if (negative_first_alias_corrections <= 20) {
                                    std::cout
                                        << "NEGATIVE_FIRST_ALIAS_CORRECTION"
                                        << " level=" << level
                                        << " m=" << m << " n=" << n
                                        << " c=" << c
                                        << " correction=" << correction
                                        << '\n';
                                }
                                if (!reported_hard_first_alias_correction
                                    && m >= 2 && n >= 2
                                    && (c & 1) != 0) {
                                    std::cout
                                        << "HARD_NEGATIVE_FIRST_ALIAS_CORRECTION"
                                        << " level=" << level
                                        << " m=" << m << " n=" << n
                                        << " c=" << c
                                        << " total=" << numerator
                                        << " continuous_part="
                                        << 2 * order * order
                                            * constant_coefficient
                                        << " correction=" << correction
                                        << " corner_coefficient="
                                        << correction
                                            / (8 * order * order)
                                        << '\n';
                                    reported_hard_first_alias_correction = true;
                                }
                            }
                            ++first_alias_cases;
                        }
                        ++cases;
                    }
                }
            }
        }
        std::cout
            << "PASS cases=" << cases
            << " negative_sectors=" << negative_sectors
            << " first_alias_cases=" << first_alias_cases
            << " negative_first_alias_corrections="
            << negative_first_alias_corrections
            << " minimum_transport_margin="
            << static_cast<double>(minimum_transport_margin)
            << " minimum_odd_rounded_ratio="
            << static_cast<double>(minimum_odd_rounded_ratio) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
