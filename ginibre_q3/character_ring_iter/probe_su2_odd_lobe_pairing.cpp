#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value <= 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error(
            std::string(name) + " must be a positive integer"
        );
    }
    return static_cast<int>(value);
}

long double pairing_weight(int half_label, long double theta) {
    const long double order =
        static_cast<long double>(half_label);
    const long double sine = std::sin(theta);
    const long double product =
        std::sin(order * theta)
        * std::cos((order + 1.0L) * theta);
    const long double weight =
        2.0L * order - 2.0L * product / sine;
    return weight / (sine * sine * sine);
}

long double odd_margin_weight(
    int half_label,
    int half_power,
    long double theta
) {
    const long double order =
        static_cast<long double>(half_label);
    const long double sine = std::sin(theta);
    const long double product =
        std::sin(order * theta)
        * std::cos((order + 1.0L) * theta);
    const long double weight =
        2.0L * order - 2.0L * product / sine;
    return weight / std::pow(sine, 2 * half_power - 1);
}

int analyze_twisted_components(
    int maximum_q_half,
    int maximum_power,
    int grid
) {
    const long double pi = std::acos(-1.0L);
    std::size_t rows = 0;
    std::size_t negative = 0;
    long double minimum =
        std::numeric_limits<long double>::infinity();
    int witness_q_half = 0;
    int witness_power = 0;
    int witness_frequency = 0;
    int witness_grid = 0;

    for (int q_half = 1; q_half <= maximum_q_half; ++q_half) {
        const int lobes = 2 * q_half + 1;
        for (int power = 1; power <= maximum_power; ++power) {
            for (int frequency = 1;
                 frequency <= 2 * lobes - 1; frequency += 2) {
                for (int point = 1; point <= grid; ++point) {
                    const long double phase =
                        pi * static_cast<long double>(point)
                        / static_cast<long double>(grid + 1);
                    long double value = 0.0L;
                    for (int lobe = 0; lobe < lobes; ++lobe) {
                        const long double theta =
                            (
                                phase
                                + static_cast<long double>(lobe) * pi
                            ) / static_cast<long double>(lobes);
                        const long double term =
                            std::sin(
                                static_cast<long double>(frequency)
                                * theta
                            )
                            / std::pow(std::sin(theta), 2 * power);
                        value += (lobe & 1) == 0 ? term : -term;
                    }
                    ++rows;
                    if (value < minimum) {
                        minimum = value;
                        witness_q_half = q_half;
                        witness_power = power;
                        witness_frequency = frequency;
                        witness_grid = point;
                    }
                    if (value < -1.0e-10L) {
                        ++negative;
                    }
                }
            }
        }
    }

    std::cout
        << std::setprecision(18)
        << "SU2_ODD_TWISTED_COMPONENT_PROBE"
        << " maximum_q_half=" << maximum_q_half
        << " maximum_power=" << maximum_power
        << " grid=" << grid
        << " rows=" << rows
        << " negative_rows=" << negative
        << " minimum=" << minimum
        << " witness={q_half=" << witness_q_half
        << " power=" << witness_power
        << " frequency=" << witness_frequency
        << " grid_point=" << witness_grid << '}'
        << " result="
        << (
            negative == 0
                ? "PASS_BOUNDED_NUMERICAL_DIAGNOSTIC"
                : "COUNTEREXAMPLE"
        )
        << '\n';
    return negative == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

int analyze_global_lobe_payment(
    int maximum_q_half,
    int maximum_half_power,
    int grid
) {
    const long double pi = std::acos(-1.0L);
    std::size_t rows = 0;
    std::size_t negative = 0;
    std::size_t negative_increment_rows = 0;
    long double minimum_margin =
        std::numeric_limits<long double>::infinity();
    int witness_q_half = 0;
    int witness_half_label = 0;
    int witness_half_power = 0;
    int witness_grid = 0;

    for (int q_half = 1; q_half <= maximum_q_half; ++q_half) {
        const int q = 2 * q_half;
        const int lobes = q + 1;
        for (int half_label = 1;
             half_label <= q; ++half_label) {
            for (int half_power = 2;
                 half_power <= maximum_half_power; ++half_power) {
                for (int point = 1; point <= grid; ++point) {
                    const long double phase =
                        pi * static_cast<long double>(point)
                        / static_cast<long double>(grid + 1);
                    long double margin = 0.0L;
                    long double increment = 0.0L;
                    for (int lobe = 0; lobe < lobes; ++lobe) {
                        const long double theta =
                            (
                                phase
                                + static_cast<long double>(lobe) * pi
                            ) / static_cast<long double>(lobes);
                        const long double term = odd_margin_weight(
                            half_label,
                            half_power,
                            theta
                        );
                        margin += (lobe & 1) == 0 ? term : -term;
                        const long double sine = std::sin(theta);
                        const long double frequency_sine = std::sin(
                            static_cast<long double>(half_label) * theta
                        );
                        const long double increment_term =
                            4.0L * frequency_sine * frequency_sine
                            / std::pow(
                                sine,
                                2 * half_power - 1
                            );
                        increment += (lobe & 1) == 0
                            ? increment_term
                            : -increment_term;
                    }
                    ++rows;
                    if (margin < minimum_margin) {
                        minimum_margin = margin;
                        witness_q_half = q_half;
                        witness_half_label = half_label;
                        witness_half_power = half_power;
                        witness_grid = point;
                    }
                    if (margin < -1.0e-10L) {
                        ++negative;
                    }
                    if (increment < -1.0e-10L) {
                        ++negative_increment_rows;
                    }
                }
            }
        }
    }

    std::cout
        << std::setprecision(18)
        << "SU2_ODD_GLOBAL_LOBE_PAYMENT_PROBE"
        << " maximum_q_half=" << maximum_q_half
        << " maximum_half_power=" << maximum_half_power
        << " grid=" << grid
        << " rows=" << rows
        << " negative_rows=" << negative
        << " negative_increment_rows=" << negative_increment_rows
        << " minimum_margin=" << minimum_margin
        << " witness={q_half=" << witness_q_half
        << " half_label=" << witness_half_label
        << " half_power=" << witness_half_power
        << " grid_point=" << witness_grid << '}'
        << " result="
        << (
            negative == 0
                ? "PASS_BOUNDED_NUMERICAL_DIAGNOSTIC"
                : "COUNTEREXAMPLE"
        )
        << '\n';
    return negative == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

int analyze_shift_pairing(int maximum_q_half, int grid) {
    const long double pi = std::acos(-1.0L);
    std::size_t rows = 0;
    std::size_t negative = 0;
    long double minimum_margin =
        std::numeric_limits<long double>::infinity();
    int witness_q_half = 0;
    int witness_half_label = 0;
    int witness_lobe = 0;
    int witness_grid = 0;

    for (int q_half = 1; q_half <= maximum_q_half; ++q_half) {
        const int q = 2 * q_half;
        const long double shift =
            pi / static_cast<long double>(q + 1);
        for (int half_label = 1;
             half_label <= q; ++half_label) {
            for (int lobe = 1;
                 lobe <= q_half; lobe += 2) {
                const long double lower =
                    static_cast<long double>(lobe) * shift;
                const long double upper = std::min(
                    static_cast<long double>(lobe + 1) * shift,
                    pi / 2.0L
                );
                for (int point = 1; point <= grid; ++point) {
                    const long double theta =
                        lower
                        + (upper - lower)
                            * static_cast<long double>(point)
                            / static_cast<long double>(grid + 1);
                    long double paired =
                        pairing_weight(half_label, theta - shift);
                    if (
                        lobe == q_half
                        && (q_half & 1) != 0
                    ) {
                        paired += pairing_weight(
                            half_label,
                            2.0L
                                * static_cast<long double>(q_half)
                                * shift
                                - theta
                        );
                    }
                    const long double margin =
                        paired - pairing_weight(half_label, theta);
                    ++rows;
                    if (margin < minimum_margin) {
                        minimum_margin = margin;
                        witness_q_half = q_half;
                        witness_half_label = half_label;
                        witness_lobe = lobe;
                        witness_grid = point;
                    }
                    if (margin < -1.0e-12L) {
                        ++negative;
                    }
                }
            }
        }
    }

    std::cout
        << std::setprecision(18)
        << "SU2_ODD_LOBE_SHIFT_PAIRING_PROBE"
        << " maximum_q_half=" << maximum_q_half
        << " grid=" << grid
        << " rows=" << rows
        << " negative_rows=" << negative
        << " minimum_margin=" << minimum_margin
        << " witness={q_half=" << witness_q_half
        << " half_label=" << witness_half_label
        << " lobe=" << witness_lobe
        << " grid_point=" << witness_grid << '}'
        << " result="
        << (
            negative == 0
                ? "PASS_BOUNDED_NUMERICAL_DIAGNOSTIC"
                : "COUNTEREXAMPLE"
        )
        << '\n';
    return negative == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (
            argc == 5
            && std::string(argv[1]) == "--components"
        ) {
            return analyze_twisted_components(
                parse_positive(argv[2], "maximum q half"),
                parse_positive(argv[3], "maximum power"),
                parse_positive(argv[4], "grid")
            );
        }
        if (
            argc == 5
            && std::string(argv[1]) == "--global"
        ) {
            return analyze_global_lobe_payment(
                parse_positive(argv[2], "maximum q half"),
                parse_positive(argv[3], "maximum half-power"),
                parse_positive(argv[4], "grid")
            );
        }
        if (
            argc == 4
            && std::string(argv[1]) == "--shift"
        ) {
            return analyze_shift_pairing(
                parse_positive(argv[2], "maximum q half"),
                parse_positive(argv[3], "grid")
            );
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: probe_su2_odd_lobe_pairing "
                "MAXIMUM_HALF_LABEL GRID | "
                "--components MAXIMUM_Q_HALF MAXIMUM_POWER GRID | "
                "--global MAXIMUM_Q_HALF MAXIMUM_HALF_POWER GRID | "
                "--shift MAXIMUM_Q_HALF GRID"
            );
        }
        const int maximum_half_label =
            parse_positive(argv[1], "maximum half-label");
        const int grid = parse_positive(argv[2], "grid");
        const long double half_pi = std::acos(-1.0L) / 2.0L;
        long double minimum_margin =
            std::numeric_limits<long double>::infinity();
        int witness_half_label = 0;
        int witness_grid = 0;
        std::size_t negative = 0;

        for (int half_label = 1;
             half_label <= maximum_half_label; ++half_label) {
            for (int point = 1; point <= grid; ++point) {
                const long double theta =
                    half_pi * static_cast<long double>(point)
                    / static_cast<long double>(grid);
                const long double order =
                    static_cast<long double>(half_label);
                const long double sine = std::sin(theta);
                const long double cosine = std::cos(theta);
                const long double left = std::sin(order * theta);
                const long double right =
                    std::cos((order + 1.0L) * theta);
                const long double product = left * right;
                const long double product_derivative =
                    order * std::cos(order * theta) * right
                    - (order + 1.0L) * left
                        * std::sin((order + 1.0L) * theta);
                const long double weight =
                    2.0L * order - 2.0L * product / sine;
                const long double derivative =
                    -2.0L * (
                        product_derivative * sine
                        - product * cosine
                    ) / (sine * sine);
                const long double margin =
                    3.0L * cosine * weight - sine * derivative;
                if (margin < minimum_margin) {
                    minimum_margin = margin;
                    witness_half_label = half_label;
                    witness_grid = point;
                }
                if (margin < -1.0e-12L) {
                    ++negative;
                }
            }
        }

        std::cout
            << std::setprecision(18)
            << "SU2_ODD_LOBE_PAIRING_PROBE"
            << " maximum_half_label=" << maximum_half_label
            << " grid=" << grid
            << " negative_rows=" << negative
            << " minimum_margin=" << minimum_margin
            << " witness={half_label=" << witness_half_label
            << " grid_point=" << witness_grid << '}'
            << " result="
            << (
                negative == 0
                    ? "PASS_BOUNDED_NUMERICAL_DIAGNOSTIC"
                    : "COUNTEREXAMPLE"
            )
            << '\n';
        return negative == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_ODD_LOBE_PAIRING_PROBE FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
