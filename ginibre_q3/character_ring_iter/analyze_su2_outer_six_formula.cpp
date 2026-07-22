#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

long long elementary(const std::vector<int>& signs, int degree) {
    std::vector<long long> coefficients(
        static_cast<std::size_t>(degree + 1), 0
    );
    coefficients[0] = 1;
    int used = 0;
    for (int sign : signs) {
        const int upper = std::min(used + 1, degree);
        for (int index = upper; index >= 1; --index) {
            coefficients[static_cast<std::size_t>(index)]
                += static_cast<long long>(sign)
                    * coefficients[static_cast<std::size_t>(index - 1)];
        }
        ++used;
    }
    return coefficients[static_cast<std::size_t>(degree)];
}

std::vector<int> signs(int count, int negative) {
    std::vector<int> result(static_cast<std::size_t>(count), 1);
    for (int index = 0; index < negative; ++index) {
        result[static_cast<std::size_t>(index)] = -1;
    }
    return result;
}

long long value(
    int ones,
    int twos,
    int threes,
    int negative_ones,
    int negative_twos,
    int negative_threes
) {
    const auto first = signs(ones, negative_ones);
    const auto second = signs(twos, negative_twos);
    const auto third = signs(threes, negative_threes);
    const long long n = ones + twos + threes;
    const long long a1 = elementary(first, 1);
    const long long a2 = elementary(first, 2);
    const long long a3 = elementary(first, 3);
    const long long a4 = elementary(first, 4);
    const long long a6 = elementary(first, 6);
    const long long b1 = elementary(second, 1);
    const long long b2 = elementary(second, 2);
    const long long b3 = elementary(second, 3);
    const long long c1 = elementary(third, 1);
    const long long c2 = elementary(third, 2);
    const long long empty = (n + 1) * n * (n - 1) / 6
        - static_cast<long long>(ones) * (n - 1) - twos;
    const long long after_two
        = (n - 2) * (n - 5) / 2 + n - ones;
    return empty + after_two * a2
        + (n - 3) * b2 + (n - 4) * a2 * b1
        + 2 * (n - 5) * a4 + c2 + a1 * b1 * c1
        + b3 + a3 * c1 + 2 * a2 * b2
        + 3 * a4 * b1 + 5 * a6;
}

long long binomial(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    long long result = 1;
    for (int index = 1; index <= k; ++index) {
        result = result * (n - k + index) / index;
    }
    return result;
}

long long outer_constant(int ones) {
    const long long first_two = binomial(ones, 2);
    const long long first_three = binomial(ones, 3);
    const long long first_four = binomial(ones, 4);
    const long long first_five = binomial(ones, 5);
    const long long first_six = binomial(ones, 6);
    return -first_two - 2 * first_three + 4 * first_four
        + 10 * first_five + 5 * first_six;
}

long long hostile_base_polynomial(int ones, int twos) {
    const long long first_two = binomial(ones, 2);
    const long long first_four = binomial(ones, 4);
    return 3 * binomial(twos, 3)
        + (first_two + 2LL * ones) * binomial(twos, 2)
        + (2 * first_two - first_four - 1) * twos
        + outer_constant(ones);
}

long long friendly_base_polynomial(int ones, int twos) {
    const long long first_two = binomial(ones, 2);
    const long long first_three = binomial(ones, 3);
    const long long first_four = binomial(ones, 4);
    return 5 * binomial(twos, 3)
        + (2LL * ones + 5 * first_two) * binomial(twos, 2)
        + (-1 + 6 * first_three + 5 * first_four) * twos
        + outer_constant(ones);
}

long long hostile_discriminant_polynomial(int ones) {
    const long long first_two = binomial(ones, 2);
    const long long first_four = binomial(ones, 4);
    const long long a = first_two + 2LL * ones;
    const long long ell = 2 * first_two - first_four - 1;
    const long long constant = outer_constant(ones);
    const long long doubled_linear = 2 * ell - a;
    return 8 * a * constant - doubled_linear * doubled_linear;
}

long long univariate_difference(int base, int degree) {
    long long result = 0;
    for (int index = 0; index <= degree; ++index) {
        const long long term
            = binomial(degree, index)
                * hostile_discriminant_polynomial(base + index);
        result += ((degree - index) & 1) != 0 ? -term : term;
    }
    return result;
}

long long base_t_difference(
    int ones,
    int base_twos,
    int degree,
    bool friendly
) {
    long long result = 0;
    for (int index = 0; index <= degree; ++index) {
        const long long value_at_index = friendly
            ? friendly_base_polynomial(ones, base_twos + index)
            : hostile_base_polynomial(ones, base_twos + index);
        const long long term = binomial(degree, index) * value_at_index;
        result += ((degree - index) & 1) != 0 ? -term : term;
    }
    return result;
}

long long outer_constant_difference(int base, int degree) {
    long long result = 0;
    for (int index = 0; index <= degree; ++index) {
        const long long term
            = binomial(degree, index) * outer_constant(base + index);
        result += ((degree - index) & 1) != 0 ? -term : term;
    }
    return result;
}

bool admissible(int ones, int twos, int threes) {
    const int total_label = ones + 2 * twos + 3 * threes;
    const int target = total_label - 6;
    const int maximum_label = threes != 0 ? 3
        : twos != 0 ? 2 : ones != 0 ? 1 : 0;
    return maximum_label != 0 && target >= maximum_label;
}

long long uniform_value(
    int ones,
    int twos,
    int threes,
    int beta,
    int delta
) {
    const int first_sign = 1;
    const int second_sign = beta;
    const int third_sign = delta;
    return value(
        ones, twos, threes,
        first_sign < 0 ? ones : 0,
        second_sign < 0 ? twos : 0,
        third_sign < 0 ? threes : 0
    );
}

long long finite_difference(
    int base_ones,
    int base_twos,
    int base_threes,
    int degree_ones,
    int degree_twos,
    int degree_threes,
    int beta,
    int delta
) {
    long long result = 0;
    for (int first = 0; first <= degree_ones; ++first) {
        for (int second = 0; second <= degree_twos; ++second) {
            for (int third = 0; third <= degree_threes; ++third) {
                const int parity = degree_ones - first
                    + degree_twos - second + degree_threes - third;
                const long long coefficient
                    = binomial(degree_ones, first)
                    * binomial(degree_twos, second)
                    * binomial(degree_threes, third);
                const long long term = coefficient * uniform_value(
                    base_ones + first,
                    base_twos + second,
                    base_threes + third,
                    beta,
                    delta
                );
                result += (parity & 1) != 0 ? -term : term;
            }
        }
    }
    return result;
}

constexpr std::array<int, 3> polynomial_degrees{6, 3, 2};

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_outer_six_formula MAXIMUM_COUNT"
            );
        }
        const int maximum_count = std::stoi(argv[1]);
        if (maximum_count < 1 || maximum_count > 100) {
            throw std::runtime_error("invalid maximum count");
        }
        long long minimum = std::numeric_limits<long long>::max();
        long long tested = 0;
        int best[6]{};
        long long uniform_minimum[4]{
            std::numeric_limits<long long>::max(),
            std::numeric_limits<long long>::max(),
            std::numeric_limits<long long>::max(),
            std::numeric_limits<long long>::max()
        };
        int uniform_best[4][6]{};
        for (int ones = 0; ones <= maximum_count; ++ones) {
            for (int twos = 0; twos + ones <= maximum_count; ++twos) {
                for (int threes = 0;
                     threes + twos + ones <= maximum_count;
                     ++threes) {
                    const int total_label = ones + 2 * twos + 3 * threes;
                    const int target = total_label - 6;
                    const int maximum_label = threes != 0 ? 3
                        : twos != 0 ? 2 : ones != 0 ? 1 : 0;
                    if (target < maximum_label || target < 1) {
                        continue;
                    }
                    for (int first_negative = 0;
                         first_negative <= ones;
                         ++first_negative) {
                        for (int second_negative = 0;
                             second_negative <= twos;
                             ++second_negative) {
                            for (int third_negative = 0;
                                 third_negative <= threes;
                                 ++third_negative) {
                                const long long current = value(
                                    ones, twos, threes,
                                    first_negative, second_negative,
                                    third_negative
                                );
                                ++tested;
                                if (current < minimum) {
                                    minimum = current;
                                    best[0] = ones;
                                    best[1] = twos;
                                    best[2] = threes;
                                    best[3] = first_negative;
                                    best[4] = second_negative;
                                    best[5] = third_negative;
                                }
                                const bool uniform_first
                                    = first_negative == 0
                                    || first_negative == ones;
                                const bool uniform_second
                                    = second_negative == 0
                                    || second_negative == twos;
                                const bool uniform_third
                                    = third_negative == 0
                                    || third_negative == threes;
                                if (uniform_first && uniform_second
                                    && uniform_third) {
                                    const int first_sign
                                        = first_negative == 0 ? 1 : -1;
                                    const int second_sign
                                        = second_negative == 0 ? 1 : -1;
                                    const int third_sign
                                        = third_negative == 0 ? 1 : -1;
                                    const int beta_index
                                        = second_sign < 0 ? 2 : 0;
                                    const int delta_index
                                        = first_sign * third_sign < 0 ? 1 : 0;
                                    const int pattern
                                        = beta_index + delta_index;
                                    if (current < uniform_minimum[pattern]) {
                                        uniform_minimum[pattern] = current;
                                        uniform_best[pattern][0] = ones;
                                        uniform_best[pattern][1] = twos;
                                        uniform_best[pattern][2] = threes;
                                        uniform_best[pattern][3]
                                            = first_negative;
                                        uniform_best[pattern][4]
                                            = second_negative;
                                        uniform_best[pattern][5]
                                            = third_negative;
                                    }
                                }
                                if (current < 0) {
                                    std::cout << "FAIL value=" << current
                                              << " ones=" << ones
                                              << " twos=" << twos
                                              << " threes=" << threes
                                              << " negative_ones="
                                              << first_negative
                                              << " negative_twos="
                                              << second_negative
                                              << " negative_threes="
                                              << third_negative << '\n';
                                    return EXIT_FAILURE;
                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << "SU2_OUTER_SIX_FORMULA PASS tested=" << tested
                  << " maximum_count=" << maximum_count
                  << " minimum=" << minimum
                  << " ones=" << best[0]
                  << " twos=" << best[1]
                  << " threes=" << best[2]
                  << " negative_ones=" << best[3]
                  << " negative_twos=" << best[4]
                  << " negative_threes=" << best[5] << '\n';
        for (int pattern = 0; pattern < 4; ++pattern) {
            std::cout << "pattern=" << pattern
                      << " minimum=" << uniform_minimum[pattern]
                      << " ones=" << uniform_best[pattern][0]
                      << " twos=" << uniform_best[pattern][1]
                      << " threes=" << uniform_best[pattern][2]
                      << " negative_ones=" << uniform_best[pattern][3]
                      << " negative_twos=" << uniform_best[pattern][4]
                      << " negative_threes=" << uniform_best[pattern][5]
                      << '\n';
        }
        std::vector<std::array<int, 3>> minimal_bases;
        for (int ones = 0; ones <= 9; ++ones) {
            for (int twos = 0; twos <= 5; ++twos) {
                for (int threes = 0; threes <= 3; ++threes) {
                    if (!admissible(ones, twos, threes)) {
                        continue;
                    }
                    bool minimal = true;
                    if (ones != 0 && admissible(ones - 1, twos, threes)) {
                        minimal = false;
                    }
                    if (twos != 0 && admissible(ones, twos - 1, threes)) {
                        minimal = false;
                    }
                    if (threes != 0 && admissible(ones, twos, threes - 1)) {
                        minimal = false;
                    }
                    if (minimal) {
                        minimal_bases.push_back({ones, twos, threes});
                    }
                }
            }
        }
        std::cout << "minimal_bases=";
        for (const auto& base : minimal_bases) {
            std::cout << ' ' << base[0] << ',' << base[1] << ',' << base[2];
        }
        std::cout << '\n';
        for (int pattern = 0; pattern < 2; ++pattern) {
            const int beta = pattern >= 2 ? -1 : 1;
            const int delta = (pattern & 1) != 0 ? -1 : 1;
            bool certificate = true;
            for (const auto& base : minimal_bases) {
                for (int first = 0; first <= 6; ++first) {
                    for (int second = 0; second <= 3; ++second) {
                        for (int third = 0; third <= 2; ++third) {
                            const long long difference = finite_difference(
                                base[0], base[1], base[2],
                                first, second, third, beta, delta
                            );
                            if (difference < 0) {
                                std::cout << "BINOMIAL_CERTIFICATE_FAIL pattern="
                                          << pattern << " base=" << base[0]
                                          << ',' << base[1] << ',' << base[2]
                                          << " degree=" << first << ','
                                          << second << ',' << third
                                          << " value=" << difference << '\n';
                                certificate = false;
                                break;
                            }
                        }
                        if (!certificate) {
                            break;
                        }
                    }
                    if (!certificate) {
                        break;
                    }
                }
                if (!certificate) {
                    break;
                }
            }
            if (!certificate) {
                return EXIT_FAILURE;
            }
            std::cout << "BINOMIAL_CERTIFICATE_PASS pattern=" << pattern
                      << " bases=" << minimal_bases.size() << '\n';
        }
        const std::vector<std::array<int, 2>> friendly_exception_pairs{
            {0, 1}, {0, 2}, {1, 1}, {2, 0}, {2, 1}, {3, 0}, {4, 0}
        };
        std::cout << "FRIENDLY_EXCEPTION_RAYS";
        for (const auto& pair : friendly_exception_pairs) {
            if (friendly_base_polynomial(pair[0], pair[1]) >= 0) {
                std::cout << " CLASSIFICATION_FAIL ones=" << pair[0]
                          << " twos=" << pair[1] << '\n';
                return EXIT_FAILURE;
            }
            int minimum_threes = 0;
            while (!admissible(pair[0], pair[1], minimum_threes)) {
                ++minimum_threes;
            }
            std::cout << ' ' << pair[0] << ',' << pair[1] << ','
                      << minimum_threes << ':';
            for (int degree = 0; degree <= 2; ++degree) {
                const long long coefficient = finite_difference(
                    pair[0], pair[1], minimum_threes,
                    0, 0, degree, 1, -1
                );
                if (coefficient < 0) {
                    std::cout << " FAIL degree=" << degree
                              << " value=" << coefficient << '\n';
                    return EXIT_FAILURE;
                }
                std::cout << (degree == 0 ? '[' : ',') << coefficient;
            }
            std::cout << ']';
        }
        std::cout << '\n';
        for (int pattern : {2, 3}) {
            const int beta = -1;
            const int delta = pattern == 2 ? 1 : -1;
            for (int first = 0; first <= polynomial_degrees[0]; ++first) {
                for (int second = 0;
                     second <= polynomial_degrees[1];
                     ++second) {
                    for (int third = 1;
                         third <= polynomial_degrees[2];
                         ++third) {
                        const long long coefficient = finite_difference(
                            0, 0, 0, first, second, third, beta, delta
                        );
                        if (coefficient < 0) {
                            std::cout << "U_COEFFICIENT_FAIL pattern="
                                      << pattern << " degree=" << first << ','
                                      << second << ',' << third << " value="
                                      << coefficient << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
            for (int ones = 0;
                 ones <= polynomial_degrees[0];
                 ++ones) {
                for (int twos = 0;
                     twos <= polynomial_degrees[1];
                     ++twos) {
                    const long long formula
                        = uniform_value(ones, twos, 0, beta, delta);
                    const long long reduced
                        = hostile_base_polynomial(ones, twos);
                    if (formula != reduced) {
                        std::cout << "HOSTILE_REDUCTION_FAIL pattern="
                                  << pattern << " ones=" << ones
                                  << " twos=" << twos << " formula="
                                  << formula << " reduced=" << reduced
                                  << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
            std::cout << "HOSTILE_U_REDUCTION_PASS pattern=" << pattern
                      << '\n';
        }
        constexpr std::array<int, 8> minimum_twos{
            4, 4, 3, 3, 2, 2, 1, 0
        };
        std::cout << "HOSTILE_SMALL_R_NEWTON";
        std::vector<std::array<int, 2>> negative_residual_pairs;
        for (int ones = 0; ones < 8; ++ones) {
            std::cout << ' ' << ones << ":[";
            for (int degree = 0; degree <= 3; ++degree) {
                const long long coefficient = base_t_difference(
                    ones,
                    minimum_twos[static_cast<std::size_t>(ones)],
                    degree,
                    false
                );
                if (coefficient < 0) {
                    std::cout << " FAIL ones=" << ones
                              << " degree=" << degree
                              << " value=" << coefficient << '\n';
                    return EXIT_FAILURE;
                }
                if (degree != 0) {
                    std::cout << ',';
                }
                std::cout << coefficient;
            }
            std::cout << ']';
            const int boundary
                = minimum_twos[static_cast<std::size_t>(ones)];
            for (int twos = 0; twos < boundary; ++twos) {
                if (hostile_base_polynomial(ones, twos) < 0) {
                    negative_residual_pairs.push_back({ones, twos});
                }
            }
        }
        std::cout << '\n';
        const std::vector<std::array<int, 2>> expected_negative_pairs{
            {0, 1}, {0, 2}, {1, 1}, {2, 0}, {3, 0}, {4, 0}
        };
        if (negative_residual_pairs != expected_negative_pairs) {
            std::cout << "HOSTILE_NEGATIVE_PAIR_CLASSIFICATION_FAIL\n";
            return EXIT_FAILURE;
        }
        std::cout << "HOSTILE_EXCEPTION_RAYS";
        for (const auto& pair : negative_residual_pairs) {
            int minimum_threes = 0;
            while (!admissible(pair[0], pair[1], minimum_threes)) {
                ++minimum_threes;
            }
            std::cout << ' ' << pair[0] << ',' << pair[1] << ','
                      << minimum_threes << ':';
            for (int pattern : {2, 3}) {
                const int delta = pattern == 2 ? 1 : -1;
                if (pattern == 3) {
                    std::cout << '|';
                }
                std::cout << '[';
                for (int degree = 0; degree <= 2; ++degree) {
                    const long long coefficient = finite_difference(
                        pair[0], pair[1], minimum_threes,
                        0, 0, degree, -1, delta
                    );
                    if (coefficient < 0) {
                        std::cout << " FAIL pattern=" << pattern
                                  << " degree=" << degree
                                  << " value=" << coefficient << '\n';
                        return EXIT_FAILURE;
                    }
                    if (degree != 0) {
                        std::cout << ',';
                    }
                    std::cout << coefficient;
                }
                std::cout << ']';
            }
        }
        std::cout << '\n';
        std::cout << "HOSTILE_DISCRIMINANT_NEWTON base=8";
        for (int degree = 0; degree <= 8; ++degree) {
            const long long coefficient = univariate_difference(8, degree);
            if (coefficient < 0) {
                std::cout << " FAIL degree=" << degree
                          << " value=" << coefficient << '\n';
                return EXIT_FAILURE;
            }
            std::cout << ' ' << degree << ':' << coefficient;
        }
        std::cout << '\n';
        std::cout << "OUTER_CONSTANT_NEWTON base=5";
        for (int degree = 0; degree <= 6; ++degree) {
            const long long coefficient
                = outer_constant_difference(5, degree);
            if (coefficient < 0) {
                std::cout << " FAIL degree=" << degree
                          << " value=" << coefficient << '\n';
                return EXIT_FAILURE;
            }
            std::cout << ' ' << degree << ':' << coefficient;
        }
        std::cout << '\n';
        std::cout << "SU2_OUTER_SIX_UNIFORM_CERTIFICATE PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
