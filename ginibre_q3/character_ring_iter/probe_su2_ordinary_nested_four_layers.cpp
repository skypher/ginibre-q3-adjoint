#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = std::int64_t;
using Profile = std::vector<Integer>;

int positive_argument(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long parsed = std::stol(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0) {
        throw std::invalid_argument(
            std::string(name) + " must be a positive integer");
    }
    return static_cast<int>(parsed);
}

Profile interval_profile(int left, int right, int maximum) {
    Profile result(static_cast<std::size_t>(maximum + 1), 0);
    for (int index = left; index <= right; ++index) {
        result[static_cast<std::size_t>(index)] = 1;
    }
    return result;
}

Profile product(const Profile& left, const Profile& right) {
    const int left_maximum = static_cast<int>(left.size()) - 1;
    const int right_maximum = static_cast<int>(right.size()) - 1;
    Profile result(
        static_cast<std::size_t>(left_maximum + right_maximum + 1),
        0);
    for (int first = 0; first <= left_maximum; ++first) {
        if (left[static_cast<std::size_t>(first)] == 0) {
            continue;
        }
        for (int second = 0; second <= right_maximum; ++second) {
            if (right[static_cast<std::size_t>(second)] == 0) {
                continue;
            }
            for (int target = std::abs(first - second);
                 target <= first + second;
                 ++target) {
                result[static_cast<std::size_t>(target)]
                    += left[static_cast<std::size_t>(first)]
                       * right[static_cast<std::size_t>(second)];
            }
        }
    }
    return result;
}

Integer value(const Profile& profile, int index) {
    if (index < 0
        || index >= static_cast<int>(profile.size())) {
        return 0;
    }
    return profile[static_cast<std::size_t>(index)];
}

Integer cross_current(
    const Profile& first,
    const Profile& second,
    int radius,
    int target) {
    Integer first_interval = 0;
    Integer second_interval = 0;
    for (int label = std::abs(radius - target);
         label <= radius + target;
         ++label) {
        first_interval += value(first, label);
        second_interval += value(second, label);
    }
    return value(first, 0) * second_interval
        + value(second, 0) * first_interval
        - value(first, radius) * value(second, target)
        - value(second, radius) * value(first, target);
}

std::string render(const std::array<int, 8>& endpoints) {
    std::string output = "[";
    for (std::size_t index = 0; index < endpoints.size(); ++index) {
        if (index != 0U) {
            output += ",";
        }
        output += std::to_string(endpoints[index]);
    }
    return output + "]";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum = argc >= 2
            ? positive_argument(argv[1], "maximum_endpoint")
            : 10;
        if (argc > 2) {
            throw std::invalid_argument(
                "usage: probe_su2_ordinary_nested_four_layers "
                "[maximum_endpoint]");
        }

        std::uint64_t chains = 0U;
        std::uint64_t current_checks = 0U;
        std::uint64_t current_failures = 0U;
        std::array<int, 8> first_endpoints{};
        int first_radius = -1;
        int first_target = -1;
        Integer first_value = 0;

        for (int left_1 = 0; left_1 <= maximum; ++left_1) {
            for (int left_2 = left_1; left_2 <= maximum; ++left_2) {
                for (int left_3 = left_2; left_3 <= maximum; ++left_3) {
                    for (int left_4 = left_3;
                         left_4 <= maximum;
                         ++left_4) {
                        for (int right_4 = left_4;
                             right_4 <= maximum;
                             ++right_4) {
                            for (int right_3 = right_4;
                                 right_3 <= maximum;
                                 ++right_3) {
                                for (int right_2 = right_3;
                                     right_2 <= maximum;
                                     ++right_2) {
                                    for (int right_1 = right_2;
                                         right_1 <= maximum;
                                         ++right_1) {
                                        const std::array<int, 8> endpoints{
                                            left_1,
                                            left_2,
                                            left_3,
                                            left_4,
                                            right_4,
                                            right_3,
                                            right_2,
                                            right_1};
                                        const std::array<Profile, 4> layers{
                                            interval_profile(
                                                left_1,
                                                right_1,
                                                maximum),
                                            interval_profile(
                                                left_2,
                                                right_2,
                                                maximum),
                                            interval_profile(
                                                left_3,
                                                right_3,
                                                maximum),
                                            interval_profile(
                                                left_4,
                                                right_4,
                                                maximum)};
                                        const Profile product_12
                                            = product(layers[0], layers[1]);
                                        const Profile product_13
                                            = product(layers[0], layers[2]);
                                        const Profile product_14
                                            = product(layers[0], layers[3]);
                                        const Profile product_23
                                            = product(layers[1], layers[2]);
                                        const Profile product_24
                                            = product(layers[1], layers[3]);
                                        const Profile product_34
                                            = product(layers[2], layers[3]);
                                        ++chains;
                                        for (int radius = 0;
                                             radius <= 2 * maximum;
                                             ++radius) {
                                            for (int target = radius;
                                                 target <= 2 * maximum;
                                                 ++target) {
                                                const Integer polarized
                                                    = cross_current(
                                                          product_12,
                                                          product_34,
                                                          radius,
                                                          target)
                                                      + cross_current(
                                                          product_13,
                                                          product_24,
                                                          radius,
                                                          target)
                                                      + cross_current(
                                                          product_14,
                                                          product_23,
                                                          radius,
                                                          target);
                                                ++current_checks;
                                                if (polarized < 0) {
                                                    ++current_failures;
                                                    if (first_radius < 0) {
                                                        first_endpoints
                                                            = endpoints;
                                                        first_radius = radius;
                                                        first_target = target;
                                                        first_value
                                                            = polarized;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_ORDINARY_NESTED_FOUR_LAYERS"
            << " maximum_endpoint=" << maximum
            << " chains=" << chains
            << " current_checks=" << current_checks
            << " current_failures=" << current_failures
            << " first_endpoints=" << render(first_endpoints)
            << " first_radius=" << first_radius
            << " first_target=" << first_target
            << " first_value=" << first_value
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
