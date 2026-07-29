#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(
            std::string(name) + " must be a positive integer"
        );
    }
    return static_cast<int>(value);
}

int folded_mode(int mode_number, int multiplier, int period) {
    int value = (mode_number * multiplier) % period;
    if (value < 0) {
        value += period;
    }
    const int half_period = period / 2;
    if (value > half_period) {
        value = period - value;
    }
    if (value <= 0 || value >= half_period) {
        throw std::runtime_error("Galois fold left the mode set");
    }
    return value - 1;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 4) {
            throw std::runtime_error(
                "usage: probe_su2_terminal_galois_orbits "
                "MAXIMUM_LEVEL MAXIMUM_PAIR_POWER [MINIMUM_LABEL]"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_power =
            parse_positive(argv[2], "maximum pair power");
        const int minimum_label =
            argc == 4 ? parse_positive(argv[3], "minimum label") : 2;
        if (maximum_level < 6 || maximum_power < 4) {
            throw std::runtime_error(
                "require maximum level>=6 and maximum pair power>=4"
            );
        }
        if ((minimum_label & 1) != 0) {
            throw std::runtime_error("minimum label must be even");
        }

        const long double pi = std::acos(-1.0L);
        std::size_t parameter_rows = 0;
        std::size_t orbit_rows = 0;
        std::size_t orbit_power_rows = 0;
        std::size_t negative_orbit_powers = 0;
        long double minimum_normalized_orbit =
            std::numeric_limits<long double>::infinity();
        for (int level = 6; level <= maximum_level; level += 2) {
            const int width = level + 1;
            const int period = 2 * (level + 2);
            std::vector<int> units;
            for (int multiplier = 1;
                 multiplier < period; ++multiplier) {
                if (std::gcd(multiplier, period) == 1) {
                    units.push_back(multiplier);
                }
            }
            const long double normalization =
                std::sqrt(2.0L / static_cast<long double>(level + 2));
            std::vector<std::vector<long double>> transform(
                static_cast<std::size_t>(width),
                std::vector<long double>(
                    static_cast<std::size_t>(width)
                )
            );
            for (int row = 0; row <= level; ++row) {
                for (int mode = 0; mode <= level; ++mode) {
                    transform[static_cast<std::size_t>(row)]
                             [static_cast<std::size_t>(mode)] =
                        normalization
                        * std::sin(
                            static_cast<long double>(
                                (row + 1) * (mode + 1)
                            ) * pi
                            / static_cast<long double>(level + 2)
                        );
                }
            }

            for (int label = minimum_label;
                 2 * label < level; label += 2) {
                ++parameter_rows;
                std::vector<long double> fusion_eigenvalue(
                    static_cast<std::size_t>(width)
                );
                for (int mode = 0; mode <= level; ++mode) {
                    fusion_eigenvalue[static_cast<std::size_t>(mode)] =
                        transform[static_cast<std::size_t>(label)]
                                 [static_cast<std::size_t>(mode)]
                        / transform[0][static_cast<std::size_t>(mode)];
                }

                std::set<std::pair<int, int>> unseen;
                for (int first = 0; first <= level; ++first) {
                    for (int second = first + 1;
                         second <= level; ++second) {
                        unseen.emplace(first, second);
                    }
                }
                while (!unseen.empty()) {
                    const auto seed = *unseen.begin();
                    std::set<std::pair<int, int>> orbit;
                    for (const int multiplier : units) {
                        int first = folded_mode(
                            seed.first + 1,
                            multiplier,
                            period
                        );
                        int second = folded_mode(
                            seed.second + 1,
                            multiplier,
                            period
                        );
                        if (first == second) {
                            throw std::runtime_error(
                                "Galois orbit collapsed a mode pair"
                            );
                        }
                        if (first > second) {
                            std::swap(first, second);
                        }
                        orbit.emplace(first, second);
                    }
                    for (const auto& pair : orbit) {
                        unseen.erase(pair);
                    }
                    ++orbit_rows;

                    std::vector<long double> eigen_sums;
                    std::vector<long double> coefficients;
                    for (const auto& pair : orbit) {
                        const int first = pair.first;
                        const int second = pair.second;
                        const long double sum =
                            fusion_eigenvalue[
                                static_cast<std::size_t>(first)
                            ]
                            + fusion_eigenvalue[
                                static_cast<std::size_t>(second)
                            ];
                        const long double seed_minor =
                            transform[
                                static_cast<std::size_t>(label)
                            ][static_cast<std::size_t>(first)]
                            * transform[0][
                                static_cast<std::size_t>(second)
                            ]
                            - transform[
                                static_cast<std::size_t>(label)
                            ][static_cast<std::size_t>(second)]
                            * transform[0][
                                static_cast<std::size_t>(first)
                            ];
                        const int first_sign =
                            (first & 1) == 0 ? 1 : -1;
                        const int second_sign =
                            (second & 1) == 0 ? 1 : -1;
                        const long double target_minor =
                            static_cast<long double>(
                                first_sign - second_sign
                            )
                            * transform[0][
                                static_cast<std::size_t>(first)
                            ]
                            * transform[0][
                                static_cast<std::size_t>(second)
                            ];
                        eigen_sums.push_back(sum);
                        coefficients.push_back(
                            seed_minor * target_minor
                        );
                    }

                    for (int power = 4;
                         power <= maximum_power; ++power) {
                        const int exponent = 2 * power + 1;
                        long double contribution = 0.0L;
                        long double absolute_sum = 0.0L;
                        for (std::size_t index = 0;
                             index < orbit.size(); ++index) {
                            const long double term =
                                coefficients[index]
                                * std::pow(
                                    eigen_sums[index],
                                    exponent
                                );
                            contribution += term;
                            absolute_sum += std::abs(term);
                        }
                        const long double normalized =
                            contribution / (1.0L + absolute_sum);
                        minimum_normalized_orbit = std::min(
                            minimum_normalized_orbit,
                            normalized
                        );
                        ++orbit_power_rows;
                        if (normalized < -1.0e-12L) {
                            ++negative_orbit_powers;
                            std::cout
                                << "SU2_TERMINAL_GALOIS_ORBITS"
                                << " result=FAIL_NUMERICAL"
                                << " level=" << level
                                << " label=" << label
                                << " pair_power=" << power
                                << " exponent=" << exponent
                                << " orbit_size=" << orbit.size()
                                << " seed_pair=(" << seed.first
                                << ',' << seed.second << ')'
                                << " contribution="
                                << static_cast<double>(contribution)
                                << " normalized="
                                << static_cast<double>(normalized)
                                << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_TERMINAL_GALOIS_ORBITS"
            << " maximum_level=" << maximum_level
            << " maximum_pair_power=" << maximum_power
            << " minimum_label=" << minimum_label
            << " parameter_rows=" << parameter_rows
            << " orbit_rows=" << orbit_rows
            << " orbit_power_rows=" << orbit_power_rows
            << " negative_orbit_powers=" << negative_orbit_powers
            << " minimum_normalized_orbit="
            << static_cast<double>(minimum_normalized_orbit)
            << " result=PASS_NUMERICAL_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_GALOIS_ORBITS FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
