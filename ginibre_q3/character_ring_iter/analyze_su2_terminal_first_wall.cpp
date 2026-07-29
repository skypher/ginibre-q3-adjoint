#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

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

bool fuses(int level, int first, int second, int output) {
    return std::abs(first - second) <= output
        && output <= std::min(
            first + second, 2 * level - first - second
        )
        && ((first + second + output) & 1) == 0;
}

std::vector<std::vector<Integer>> path_table(
    int level,
    int label,
    int maximum_power
) {
    std::vector<std::vector<Integer>> table(
        static_cast<std::size_t>(maximum_power + 1),
        std::vector<Integer>(static_cast<std::size_t>(level + 1))
    );
    table[0][0] = 1;
    for (int power = 1; power <= maximum_power; ++power) {
        for (int source = 0; source <= level; ++source) {
            const Integer& coefficient =
                table[static_cast<std::size_t>(power - 1)]
                     [static_cast<std::size_t>(source)];
            if (coefficient == 0) {
                continue;
            }
            for (int output = 0; output <= level; ++output) {
                if (fuses(level, label, source, output)) {
                    table[static_cast<std::size_t>(power)]
                         [static_cast<std::size_t>(output)]
                        += coefficient;
                }
            }
        }
    }
    return table;
}

std::vector<std::vector<Integer>> binomial_table(int maximum) {
    std::vector<std::vector<Integer>> choose(
        static_cast<std::size_t>(maximum + 1)
    );
    for (int n = 0; n <= maximum; ++n) {
        choose[static_cast<std::size_t>(n)].assign(
            static_cast<std::size_t>(n + 1), Integer{0}
        );
        choose[static_cast<std::size_t>(n)][0] = 1;
        choose[static_cast<std::size_t>(n)]
              [static_cast<std::size_t>(n)] = 1;
        for (int r = 1; r < n; ++r) {
            choose[static_cast<std::size_t>(n)]
                  [static_cast<std::size_t>(r)] =
                choose[static_cast<std::size_t>(n - 1)]
                      [static_cast<std::size_t>(r - 1)]
                + choose[static_cast<std::size_t>(n - 1)]
                        [static_cast<std::size_t>(r)];
        }
    }
    return choose;
}

Integer weight(
    const std::vector<std::vector<Integer>>& choose,
    int length,
    int complementary
) {
    return Integer(length - 2 * complementary)
        * choose[static_cast<std::size_t>(length)]
                [static_cast<std::size_t>(complementary)]
        / length;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_terminal_first_wall MAXIMUM_LEVEL"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        if (maximum_level < 6) {
            throw std::runtime_error("require maximum level>=6");
        }
        const auto choose = binomial_table(maximum_level + 2);

        std::size_t rows = 0;
        std::size_t profile_identity_rows = 0;
        std::size_t negative_determinants = 0;
        std::size_t minimal_profile_minor_rows = 0;
        std::size_t negative_minimal_profile_minors = 0;
        std::size_t adjacent_profile_rows = 0;
        std::size_t negative_adjacent_profile_minors = 0;
        std::size_t failed_low_dimension_bounds = 0;
        std::size_t failed_high_dimension_bounds = 0;
        std::size_t high_increment_rows = 0;
        std::size_t negative_high_increment_bounds = 0;
        std::size_t odd_low_increment_rows = 0;
        std::size_t failed_odd_low_increment_bounds = 0;
        std::size_t negative_unit_payments = 0;
        std::size_t negative_adjacent_payments = 0;
        std::size_t negative_currents = 0;
        int maximum_payment_depth = 0;
        int adjacent_witness_level = 0;
        int adjacent_witness_label = 0;
        int depth_witness_level = 0;
        int depth_witness_label = 0;
        int maximum_negative_distance = 0;
        int maximum_high_dimension_failure_distance = 0;
        int high_bound_witness_level = 0;
        int high_bound_witness_label = 0;
        int high_bound_witness_channel = 0;
        Integer high_bound_witness_left = 0;
        Integer high_bound_witness_right = 0;
        bool closest_profile_present = false;
        int closest_profile_level = 0;
        int closest_profile_label = 0;
        int closest_profile_distance = 0;
        int closest_profile_channel = 0;
        Integer closest_profile_numerator = 0;
        Integer closest_profile_denominator = 0;
        Integer adjacent_witness_value = 0;

        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                const int distance = (level + label - 1) / label;
                const int length = 2 * distance + 2;
                const auto paths = path_table(level, label, length);
                const auto f = [&](int power) -> const Integer& {
                    return paths[static_cast<std::size_t>(power)][0];
                };
                const auto g = [&](int power) -> const Integer& {
                    return paths[static_cast<std::size_t>(power)]
                                [static_cast<std::size_t>(level)];
                };

                const Integer determinant =
                    g(length - distance) * f(distance)
                    - g(distance) * f(length - distance);
                Integer profile_sum = 0;
                for (int channel = 0;
                     channel <= 2 * label; channel += 2) {
                    const Integer profile_minor =
                        f(distance)
                            * paths[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(
                                level - channel
                            )]
                        - g(distance)
                            * paths[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(channel)];
                    profile_sum += profile_minor;
                    if (distance >= 4) {
                        ++minimal_profile_minor_rows;
                        if (profile_minor < 0) {
                            ++negative_minimal_profile_minors;
                        }
                        if (
                            paths[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(channel)]
                            > Integer(channel + 1) * f(distance)
                        ) {
                            ++failed_low_dimension_bounds;
                        }
                        if (
                            paths[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(
                                level - channel
                            )]
                            < Integer(channel + 1) * g(distance)
                        ) {
                            ++failed_high_dimension_bounds;
                            maximum_high_dimension_failure_distance =
                                std::max(
                                    maximum_high_dimension_failure_distance,
                                    distance
                                );
                            if (high_bound_witness_level == 0) {
                                high_bound_witness_level = level;
                                high_bound_witness_label = label;
                                high_bound_witness_channel = channel;
                                high_bound_witness_left =
                                    paths[
                                        static_cast<std::size_t>(
                                            distance
                                        )
                                    ][static_cast<std::size_t>(
                                        level - channel
                                    )];
                                high_bound_witness_right =
                                    Integer(channel + 1) * g(distance);
                            }
                        }
                        const Integer profile_denominator =
                            g(distance)
                            * paths[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(channel)];
                        const Integer profile_numerator =
                            f(distance)
                            * paths[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(
                                level - channel
                            )];
                        if (
                            channel > 0
                            && profile_denominator > 0
                            && (
                                !closest_profile_present
                                || profile_numerator
                                        * closest_profile_denominator
                                    < closest_profile_numerator
                                        * profile_denominator
                            )
                        ) {
                            closest_profile_present = true;
                            closest_profile_level = level;
                            closest_profile_label = label;
                            closest_profile_distance = distance;
                            closest_profile_channel = channel;
                            closest_profile_numerator =
                                profile_numerator;
                            closest_profile_denominator =
                                profile_denominator;
                        }
                    }
                }
                if (profile_sum != determinant) {
                    throw std::runtime_error(
                        "first-wall profile identity failed"
                    );
                }
                ++profile_identity_rows;
                if (distance >= 4) {
                    for (int channel = 0;
                         channel < 2 * label; channel += 2) {
                        const Integer adjacent_profile_minor =
                            paths[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(
                                level - channel - 2
                            )]
                                * paths[
                                    static_cast<std::size_t>(distance)
                                ][static_cast<std::size_t>(channel)]
                            - paths[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(
                                level - channel
                            )]
                                * paths[
                                    static_cast<std::size_t>(distance)
                                ][static_cast<std::size_t>(
                                    channel + 2
                                )];
                        ++adjacent_profile_rows;
                        if (adjacent_profile_minor < 0) {
                            ++negative_adjacent_profile_minors;
                        }
                    }
                }
                if (distance >= 5) {
                    for (int channel = 2;
                         channel <= 2 * label; channel += 2) {
                        const Integer high_increment =
                            paths[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(
                                level - channel
                            )]
                            - paths[
                                static_cast<std::size_t>(distance)
                            ][static_cast<std::size_t>(
                                level - channel + 2
                            )];
                        ++high_increment_rows;
                        if (high_increment < 2 * g(distance)) {
                            ++negative_high_increment_bounds;
                        }
                        if ((distance & 1) != 0) {
                            const Integer low_increment =
                                paths[
                                    static_cast<std::size_t>(distance)
                                ][static_cast<std::size_t>(channel)]
                                - paths[
                                    static_cast<std::size_t>(distance)
                                ][static_cast<std::size_t>(
                                    channel - 2
                                )];
                            ++odd_low_increment_rows;
                            if (low_increment > 2 * f(distance)) {
                                ++failed_odd_low_increment_bounds;
                            }
                        }
                    }
                }
                if (determinant < 0) {
                    ++negative_determinants;
                    maximum_negative_distance =
                        std::max(maximum_negative_distance, distance);
                }
                const Integer unit_payment =
                    determinant
                    + g(length - distance + 1)
                        * f(distance - 1);
                if (unit_payment < 0) {
                    ++negative_unit_payments;
                }

                Integer suffix = weight(choose, length, distance)
                    * determinant;
                const Integer adjacent =
                    suffix
                    + weight(choose, length, distance - 1)
                        * g(length - distance + 1)
                        * f(distance - 1);
                if (adjacent < 0) {
                    ++negative_adjacent_payments;
                    if (adjacent_witness_level == 0) {
                        adjacent_witness_level = level;
                        adjacent_witness_label = label;
                        adjacent_witness_value = adjacent;
                    }
                }

                int payment_depth = 0;
                if (suffix < 0) {
                    payment_depth = -1;
                    for (int complementary = distance - 1;
                         complementary >= 0; --complementary) {
                        suffix += weight(
                                      choose, length, complementary
                                  )
                            * g(length - complementary)
                            * f(complementary);
                        if (suffix >= 0) {
                            payment_depth =
                                distance - complementary;
                            break;
                        }
                    }
                    if (payment_depth < 0) {
                        ++negative_currents;
                    }
                }
                if (payment_depth > maximum_payment_depth) {
                    maximum_payment_depth = payment_depth;
                    depth_witness_level = level;
                    depth_witness_label = label;
                }
                ++rows;
            }
        }

        std::cout
            << "SU2_TERMINAL_FIRST_WALL"
            << " maximum_level=" << maximum_level
            << " rows=" << rows
            << " profile_identity_rows=" << profile_identity_rows
            << " negative_determinants=" << negative_determinants
            << " minimal_profile_minor_rows="
            << minimal_profile_minor_rows
            << " negative_minimal_profile_minors="
            << negative_minimal_profile_minors
            << " adjacent_profile_rows=" << adjacent_profile_rows
            << " negative_adjacent_profile_minors="
            << negative_adjacent_profile_minors
            << " failed_low_dimension_bounds="
            << failed_low_dimension_bounds
            << " failed_high_dimension_bounds="
            << failed_high_dimension_bounds
            << " high_increment_rows=" << high_increment_rows
            << " negative_high_increment_bounds="
            << negative_high_increment_bounds
            << " odd_low_increment_rows=" << odd_low_increment_rows
            << " failed_odd_low_increment_bounds="
            << failed_odd_low_increment_bounds
            << " maximum_high_dimension_failure_distance="
            << maximum_high_dimension_failure_distance
            << " maximum_negative_distance="
            << maximum_negative_distance
            << " negative_unit_payments=" << negative_unit_payments
            << " negative_adjacent_payments="
            << negative_adjacent_payments
            << " negative_currents=" << negative_currents
            << " maximum_payment_depth=" << maximum_payment_depth;
        if (adjacent_witness_level != 0) {
            std::cout
                << " first_adjacent_failure={level="
                << adjacent_witness_level
                << " label=" << adjacent_witness_label
                << " value=" << adjacent_witness_value << '}';
        }
        std::cout
            << " depth_witness={level=" << depth_witness_level
            << " label=" << depth_witness_label << '}';
        if (high_bound_witness_level != 0) {
            std::cout
                << " first_high_dimension_failure={level="
                << high_bound_witness_level
                << " label=" << high_bound_witness_label
                << " channel=" << high_bound_witness_channel
                << " left=" << high_bound_witness_left
                << " right=" << high_bound_witness_right << '}';
        }
        if (closest_profile_present) {
            std::cout
                << " closest_profile_ratio={level="
                << closest_profile_level
                << " label=" << closest_profile_label
                << " distance=" << closest_profile_distance
                << " channel=" << closest_profile_channel
                << " numerator=" << closest_profile_numerator
                << " denominator=" << closest_profile_denominator
                << '}';
        }
        std::cout
            << " result="
            << (
                negative_currents == 0
                    && negative_high_increment_bounds == 0
                    ? "PASS_BOUNDED_DIAGNOSTIC"
                    : "FAIL"
            )
            << '\n';
        return negative_currents == 0
                && negative_high_increment_bounds == 0
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_FIRST_WALL FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
