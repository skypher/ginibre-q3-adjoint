#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

Integer binomial(int top, int order) {
    if (top < order) {
        return 0;
    }
    Integer result = 1;
    for (int index = 0; index < order; ++index) {
        result *= top - index;
        result /= index + 1;
    }
    return result;
}

Integer ordinary_six(int q, int y) {
    return binomial(6 * q - y + 4, 4)
        - 6 * binomial(4 * q - y + 3, 4)
        + 15 * binomial(2 * q - y + 2, 4);
}

Integer ordinary_seven(int q, int y) {
    return binomial(7 * q - y + 5, 5)
        - 7 * binomial(5 * q - y + 4, 5)
        + 21 * binomial(3 * q - y + 3, 5)
        - 35 * binomial(q - y + 2, 5);
}

Integer ordinary(int power, int q, int y) {
    if (y < 0) {
        return 0;
    }
    if (power == 6) {
        return ordinary_six(q, y);
    }
    if (power == 7) {
        return ordinary_seven(q, y);
    }
    throw std::runtime_error("unsupported power");
}

Integer kac_walton(
    int power,
    int q,
    int h,
    int y
) {
    const int period = 4 * q + 4 + 2 * h;
    const int support = power * q;
    Integer result = 0;
    for (int translation = 0;; ++translation) {
        const int positive = translation * period + y;
        const int negative =
            (translation + 1) * period - 1 - y;
        if (positive > support && negative > support) {
            break;
        }
        result += ordinary(power, q, positive);
        result -= ordinary(power, q, negative);
    }
    return result;
}

Integer closed_six(int q, int h, int y) {
    return ordinary_six(q, y)
        - binomial(2 * q - 2 * h + y + 1, 4)
        + 6 * binomial(y - 2 * h, 4)
        + binomial(2 * q - 2 * h - y, 4);
}

Integer closed_seven(int q, int h, int y) {
    return ordinary_seven(q, y)
        - binomial(3 * q - 2 * h + y + 2, 5)
        + 7 * binomial(q - 2 * h + y + 1, 5)
        - 21 * binomial(y - q - 2 * h, 5)
        + binomial(3 * q - 2 * h - y + 1, 5)
        - 7 * binomial(q - 2 * h - y, 5)
        - binomial(y - q - 4 * h - 2, 5);
}

std::vector<Integer> multiply(
    int half_level,
    int half_label,
    const std::vector<Integer>& state
) {
    std::vector<Integer> next(
        static_cast<std::size_t>(half_level + 1),
        Integer(0)
    );
    for (int source = 0; source <= half_level; ++source) {
        const Integer& coefficient =
            state[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        const int lower = std::abs(source - half_label);
        const int upper = std::min(
            source + half_label,
            2 * half_level - source - half_label
        );
        for (int target = lower; target <= upper; ++target) {
            next[static_cast<std::size_t>(target)] += coefficient;
        }
    }
    return next;
}

}  // namespace

int main() {
    try {
        constexpr int maximum_formula_q = 100;
        std::uint64_t formula_comparisons = 0U;
        for (int q = 1; q <= maximum_formula_q; ++q) {
            for (int h = 1; h <= 5 * q - 2; ++h) {
                const int half_level = 2 * q + 1 + h;
                for (int y = 0; y <= half_level; ++y) {
                    const Integer generic_six =
                        kac_walton(6, q, h, y);
                    const Integer generic_seven =
                        kac_walton(7, q, h, y);
                    if (
                        generic_six != closed_six(q, h, y)
                        || generic_seven != closed_seven(q, h, y)
                    ) {
                        std::cout
                            << "SU2_K3_INTERMEDIATE_FORMULA_MISMATCH"
                            << " Q=" << q
                            << " h=" << h
                            << " Y=" << y
                            << " generic6=" << generic_six
                            << " closed6=" << closed_six(q, h, y)
                            << " generic7=" << generic_seven
                            << " closed7=" << closed_seven(q, h, y)
                            << '\n';
                        throw std::runtime_error(
                            "closed formula disagrees with "
                            "generic Kac-Walton sum"
                        );
                    }
                    formula_comparisons += 2U;
                }
            }
            if (q % 10 == 0) {
                std::cout
                    << "SU2_K3_INTERMEDIATE_FORMULA_PROGRESS"
                    << " Q=" << q
                    << " maximum_Q=" << maximum_formula_q
                    << " comparisons=" << formula_comparisons
                    << std::endl;
            }
        }
        std::cout
            << "SU2_K3_INTERMEDIATE_FORMULA"
            << " maximum_Q=" << maximum_formula_q
            << " comparisons=" << formula_comparisons
            << " result=PASS_EXACT_KAC_WALTON\n";

        constexpr int maximum_recurrence_q = 15;
        std::uint64_t recurrence_comparisons = 0U;
        for (int q = 1; q <= maximum_recurrence_q; ++q) {
            for (int h = 1; h <= 5 * q - 2; ++h) {
                const int half_level = 2 * q + 1 + h;
                std::vector<Integer> state(
                    static_cast<std::size_t>(half_level + 1),
                    Integer(0)
                );
                state[0] = 1;
                for (int power = 1; power <= 7; ++power) {
                    state = multiply(half_level, q, state);
                    if (power != 6 && power != 7) {
                        continue;
                    }
                    for (int y = 0; y <= half_level; ++y) {
                        const Integer expected = power == 6
                            ? closed_six(q, h, y)
                            : closed_seven(q, h, y);
                        if (
                            state[static_cast<std::size_t>(y)]
                            != expected
                        ) {
                            std::cout
                                << "SU2_K3_INTERMEDIATE_RECURRENCE_"
                                   "MISMATCH"
                                << " Q=" << q
                                << " h=" << h
                                << " power=" << power
                                << " Y=" << y
                                << " recurrence="
                                << state[static_cast<std::size_t>(y)]
                                << " closed=" << expected << '\n';
                            throw std::runtime_error(
                                "closed formula disagrees with "
                                "finite-fusion recurrence"
                            );
                        }
                        ++recurrence_comparisons;
                    }
                }
            }
        }
        std::cout
            << "SU2_K3_INTERMEDIATE_RECURRENCE"
            << " maximum_Q=" << maximum_recurrence_q
            << " comparisons=" << recurrence_comparisons
            << " result=PASS_EXACT_FUSION\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_K3_INTERMEDIATE_AUDIT FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
