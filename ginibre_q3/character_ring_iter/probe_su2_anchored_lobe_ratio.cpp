#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

namespace {

long double even_lobe_sum(int n, int power, int frequency, long double phi) {
    long double value = 0.0L;
    for (int ell = 0; ell < n; ++ell) {
        const long double theta =
            (phi + static_cast<long double>(ell) * std::acos(-1.0L)) / n;
        value += std::sin(frequency * theta)
            / std::pow(std::sin(theta), 2 * power - 1);
    }
    return value;
}

long double odd_lobe_sum(int n, int power, int frequency, long double phi) {
    long double value = 0.0L;
    for (int ell = 0; ell < n; ++ell) {
        const long double theta =
            (phi + static_cast<long double>(ell) * std::acos(-1.0L)) / n;
        const long double sign = (ell % 2 == 0) ? 1.0L : -1.0L;
        value += sign * std::sin(frequency * theta)
            / std::pow(std::sin(theta), 2 * power);
    }
    return value;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--pair-tails") {
        constexpr int samples = 400;
        std::size_t tails = 0;
        long double minimum_relative = 0.0L;
        for (int n = 3; n <= 31; n += 2) {
            for (int target = 0; target <= n - 1; ++target) {
                const int frequency = 2 * target + 1;
                for (int sample = 1; sample < samples; ++sample) {
                    const long double phi =
                        std::acos(-1.0L) * sample / samples;
                    std::vector<std::tuple<long double, long double>> atoms;
                    for (int ell = 0; ell < n; ++ell) {
                        const long double theta_ell =
                            (phi + ell * std::acos(-1.0L)) / n;
                        const long double h_ell =
                            (ell % 2 == 0 ? 1.0L : -1.0L)
                            / std::sin(theta_ell);
                        const long double character_ell =
                            std::sin(frequency * theta_ell)
                            / std::sin(theta_ell);
                        for (int other = ell + 1; other < n; ++other) {
                            const long double theta_other =
                                (phi + other * std::acos(-1.0L)) / n;
                            const long double h_other =
                                (other % 2 == 0 ? 1.0L : -1.0L)
                                / std::sin(theta_other);
                            const long double character_other =
                                std::sin(frequency * theta_other)
                                / std::sin(theta_other);
                            const long double base =
                                h_ell * h_ell * h_other * h_other;
                            const long double coefficient =
                                (h_ell - h_other)
                                * (character_ell - character_other);
                            atoms.emplace_back(base, coefficient);
                        }
                    }
                    std::sort(
                        atoms.begin(),
                        atoms.end(),
                        [](const auto& left, const auto& right) {
                            return std::get<0>(left) < std::get<0>(right);
                        }
                    );
                    long double suffix = 0.0L;
                    long double absolute_suffix = 0.0L;
                    for (std::size_t index = atoms.size(); index > 0;) {
                        --index;
                        suffix += std::get<1>(atoms[index]);
                        absolute_suffix += std::abs(std::get<1>(atoms[index]));
                        const long double relative =
                            suffix / (absolute_suffix + 1.0L);
                        minimum_relative = std::min(
                            minimum_relative,
                            relative
                        );
                        ++tails;
                        if (relative < -1.0e-10L) {
                            std::cout
                                << std::setprecision(20)
                                << "PAIR_TAIL_COUNTEREXAMPLE"
                                << " n=" << n
                                << " target=" << target
                                << " sample=" << sample
                                << " cutoff_base="
                                << std::get<0>(atoms[index])
                                << " suffix=" << suffix
                                << " relative=" << relative << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }
        std::cout
            << std::setprecision(20)
            << "PAIR_TAIL_PROBE"
            << " tails=" << tails
            << " minimum_relative=" << minimum_relative
            << " result=PASS_BOUNDED_NUMERICAL_DIAGNOSTIC\n";
        return EXIT_SUCCESS;
    }
    if (argc == 2 && std::string(argv[1]) == "--cumulative") {
        constexpr int intervals = 8000;
        const long double step = std::acos(-1.0L) / intervals;
        std::size_t tested = 0;
        long double minimum_relative = 0.0L;
        for (int n = 3; n <= 19; n += 2) {
            for (int power = 1; power <= 10; ++power) {
                const int maximum_frequency =
                    std::min(4 * n - 1, 2 * power * (n - 1) + 1);
                const int frequencies = (maximum_frequency + 1) / 2;
                std::vector<long double> even_integrals(frequencies, 0.0L);
                std::vector<long double> odd_integrals(frequencies, 0.0L);
                for (int index = 0; index < intervals; ++index) {
                    const long double phi = (index + 0.5L) * step;
                    const long double sine = std::sin(phi);
                    const long double even_weight =
                        std::pow(sine, 2 * power);
                    const long double odd_weight = even_weight * sine;
                    for (
                        int frequency = 1;
                        frequency <= maximum_frequency;
                        frequency += 2
                    ) {
                        const std::size_t slot =
                            static_cast<std::size_t>((frequency - 1) / 2);
                        even_integrals[slot] += even_weight
                            * even_lobe_sum(n, power, frequency, phi);
                        odd_integrals[slot] += odd_weight
                            * odd_lobe_sum(n, power, frequency, phi);
                    }
                    if (index < 4) {
                        continue;
                    }
                    for (
                        int frequency = 1;
                        frequency <= maximum_frequency;
                        frequency += 2
                    ) {
                        const std::size_t slot =
                            static_cast<std::size_t>((frequency - 1) / 2);
                        const long double positive =
                            even_integrals[0] * odd_integrals[slot];
                        const long double negative =
                            odd_integrals[0] * even_integrals[slot];
                        const long double determinant = positive - negative;
                        const long double scale =
                            std::abs(positive) + std::abs(negative) + 1.0L;
                        const long double relative = determinant / scale;
                        ++tested;
                        minimum_relative = std::min(minimum_relative, relative);
                        if (relative < -1.0e-10L) {
                            std::cout
                                << std::setprecision(20)
                                << "CUMULATIVE_COUNTEREXAMPLE"
                                << " n=" << n
                                << " power=" << power
                                << " frequency=" << frequency
                                << " cutoff=" << index + 1
                                << "*pi/" << intervals
                                << " determinant=" << determinant
                                << " relative=" << relative << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }
        std::cout
            << std::setprecision(20)
            << "CUMULATIVE_PROBE"
            << " tested=" << tested
            << " minimum_relative=" << minimum_relative
            << " result=PASS_BOUNDED_NUMERICAL_DIAGNOSTIC\n";
        return EXIT_SUCCESS;
    }
    if (argc == 2 && std::string(argv[1]) == "--cumulative-special") {
        constexpr int intervals = 200000;
        const long double step = std::acos(-1.0L) / intervals;
        long double e1_integral = 0.0L;
        long double o1_integral = 0.0L;
        long double ek_integral = 0.0L;
        long double ok_integral = 0.0L;
        long double minimum = 0.0L;
        int witness = 0;
        for (int index = 0; index < intervals; ++index) {
            const long double phi = (index + 0.5L) * step;
            const long double sine = std::sin(phi);
            const long double even_weight = std::pow(sine, 4);
            const long double odd_weight = even_weight * sine;
            e1_integral += even_weight * even_lobe_sum(3, 2, 1, phi);
            o1_integral += odd_weight * odd_lobe_sum(3, 2, 1, phi);
            ek_integral += even_weight * even_lobe_sum(3, 2, 9, phi);
            ok_integral += odd_weight * odd_lobe_sum(3, 2, 9, phi);
            const long double determinant =
                e1_integral * ok_integral - o1_integral * ek_integral;
            if (determinant < minimum) {
                minimum = determinant;
                witness = index + 1;
            }
        }
        std::cout
            << std::setprecision(20)
            << "CUMULATIVE_SPECIAL"
            << " minimum=" << minimum
            << " cutoff=" << witness << "*pi/" << intervals
            << " final="
            << e1_integral * ok_integral - o1_integral * ek_integral
            << '\n';
        return EXIT_SUCCESS;
    }
    if (argc == 2 && std::string(argv[1]) == "--special") {
        {
            const long double phi = std::acos(-1.0L) / 2.0L;
            const long double e1 = even_lobe_sum(5, 2, 1, phi);
            const long double o1 = odd_lobe_sum(5, 2, 1, phi);
            const long double ek = even_lobe_sum(5, 2, 9, phi);
            const long double ok = odd_lobe_sum(5, 2, 9, phi);
            std::cout
                << std::setprecision(20)
                << "FIRST_BLOCK_CHECK"
                << " determinant=" << e1 * ok - o1 * ek
                << " e1=" << e1
                << " o1=" << o1
                << " ek=" << ek
                << " ok=" << ok << '\n';
        }
        for (int denominator = 3; denominator <= 24; ++denominator) {
            for (int numerator = 1; numerator < denominator; ++numerator) {
                const long double phi =
                    std::acos(-1.0L) * numerator / denominator;
                const long double e1 = even_lobe_sum(3, 2, 1, phi);
                const long double o1 = odd_lobe_sum(3, 2, 1, phi);
                const long double ek = even_lobe_sum(3, 2, 9, phi);
                const long double ok = odd_lobe_sum(3, 2, 9, phi);
                const long double determinant = e1 * ok - o1 * ek;
                if (determinant < -1.0e-8L) {
                    std::cout
                        << std::setprecision(20)
                        << "SPECIAL_COUNTEREXAMPLE"
                        << " phi=" << numerator << "*pi/" << denominator
                        << " determinant=" << determinant
                        << " e1=" << e1
                        << " o1=" << o1
                        << " ek=" << ek
                        << " ok=" << ok << '\n';
                    return EXIT_FAILURE;
                }
            }
        }
        return EXIT_SUCCESS;
    }
    constexpr int samples = 160;
    long double minimum_relative = std::numeric_limits<long double>::infinity();
    long double witness_determinant = 0.0L;
    int witness_n = 0;
    int witness_power = 0;
    int witness_frequency = 0;
    int witness_sample = 0;
    std::size_t tested = 0;
    std::size_t negative = 0;

    for (int n = 3; n <= 19; n += 2) {
        for (int power = 1; power <= 8; ++power) {
            const int maximum_frequency =
                std::min(2 * n - 1, 2 * power * (n - 1) + 1);
            for (
                int frequency = 1;
                frequency <= maximum_frequency;
                frequency += 2
            ) {
                for (int sample = 1; sample < samples; ++sample) {
                    const long double phi =
                        std::acos(-1.0L) * static_cast<long double>(sample)
                        / samples;
                    const long double e1 =
                        even_lobe_sum(n, power, 1, phi);
                    const long double o1 =
                        odd_lobe_sum(n, power, 1, phi);
                    const long double ek =
                        even_lobe_sum(n, power, frequency, phi);
                    const long double ok =
                        odd_lobe_sum(n, power, frequency, phi);
                    const long double determinant = e1 * ok - o1 * ek;
                    ++tested;
                    const long double scale =
                        std::abs(e1 * ok) + std::abs(o1 * ek) + 1.0L;
                    const long double relative = determinant / scale;
                    if (relative < minimum_relative) {
                        minimum_relative = relative;
                        witness_determinant = determinant;
                        witness_n = n;
                        witness_power = power;
                        witness_frequency = frequency;
                        witness_sample = sample;
                    }
                    if (relative < -1.0e-12L) {
                        ++negative;
                        std::cout
                            << std::setprecision(20)
                            << "COUNTEREXAMPLE"
                            << " n=" << n
                            << " power=" << power
                            << " frequency=" << frequency
                            << " sample=" << sample
                            << " phi=" << phi
                            << " determinant=" << determinant
                            << " relative=" << relative
                            << " e1=" << e1
                            << " o1=" << o1
                            << " ek=" << ek
                            << " ok=" << ok << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
        }
    }

    std::cout
        << std::setprecision(20)
        << "SU2_ANCHORED_LOBE_RATIO_PROBE"
        << " tested=" << tested
        << " negative=" << negative
        << " minimum_relative=" << minimum_relative
        << " witness_determinant=" << witness_determinant
        << " witness={n=" << witness_n
        << " power=" << witness_power
        << " frequency=" << witness_frequency
        << " sample=" << witness_sample << "}"
        << " result=PASS_BOUNDED_NUMERICAL_DIAGNOSTIC\n";
    return EXIT_SUCCESS;
}
