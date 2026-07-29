#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;

int parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0;
    const long long parsed = std::stoll(value, &consumed);
    if (consumed != value.size() || parsed <= 0) {
        throw std::invalid_argument(name + " must be positive");
    }
    return static_cast<int>(parsed);
}

std::vector<Integer> multiply_character(
    const std::vector<Integer>& profile,
    int q
) {
    std::vector<Integer> next(profile.size() + static_cast<std::size_t>(q));
    for (int source = 0; source < static_cast<int>(profile.size()); ++source) {
        for (int target = std::abs(source - q);
             target <= source + q;
             ++target) {
            next[static_cast<std::size_t>(target)]
                += profile[static_cast<std::size_t>(source)];
        }
    }
    return next;
}

Integer at(const std::vector<Integer>& profile, int index) {
    return
        index >= 0 && index < static_cast<int>(profile.size())
        ? profile[static_cast<std::size_t>(index)]
        : Integer{0};
}

bool dimension_normalized_log_concave(
    const std::vector<Integer>& profile,
    int index
) {
    const Integer center = at(profile, index);
    const Integer left = at(profile, index - 1);
    const Integer right = at(profile, index + 1);
    const Integer center_dimension = 2 * index + 1;
    const Integer left_dimension = 2 * index - 1;
    const Integer right_dimension = 2 * index + 3;
    return
        center * center * left_dimension * right_dimension
        >= left * right * center_dimension * center_dimension;
}

bool dimension_normalized_decreasing(
    const std::vector<Integer>& profile,
    int index
) {
    return
        at(profile, index + 1) * (2 * index + 1)
        <= at(profile, index) * (2 * index + 3);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        int maximum_q = 15;
        int maximum_half_power = 30;
        if (argc >= 2) {
            maximum_q = parse_positive(argv[1], "maximum_q");
        }
        if (argc >= 3) {
            maximum_half_power =
                parse_positive(argv[2], "maximum_half_power");
        }
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_dimension_normalized_profile"
                " [maximum_q] [maximum_half_power]"
            );
        }

        unsigned long long profiles = 0;
        unsigned long long inequalities = 0;
        unsigned long long log_concavity_failures = 0;
        unsigned long long decrease_failures = 0;
        std::string log_concavity_witness;
        std::string decrease_witness;
        for (int q = 1; q <= maximum_q; ++q) {
            std::vector<Integer> profile{Integer{1}};
            for (int step = 1; step <= 2 * maximum_half_power; ++step) {
                profile = multiply_character(profile, q);
                if ((step & 1) != 0 || step < 4) {
                    continue;
                }
                ++profiles;
                for (int index = 1;
                     index + 1 < static_cast<int>(profile.size());
                     ++index) {
                    ++inequalities;
                    if (!dimension_normalized_log_concave(profile, index)) {
                        ++log_concavity_failures;
                        if (log_concavity_witness.empty()) {
                            const Integer center =
                                at(profile, index);
                            const Integer left =
                                at(profile, index - 1);
                            const Integer right =
                                at(profile, index + 1);
                            const Integer lhs =
                                center * center
                                * (2 * index - 1)
                                * (2 * index + 3);
                            const Integer rhs =
                                left * right
                                * (2 * index + 1)
                                * (2 * index + 1);
                            log_concavity_witness =
                                "{q=" + std::to_string(q)
                                + ",power=" + std::to_string(step)
                                + ",index=" + std::to_string(index) + "}";
                            log_concavity_witness +=
                                " values={" + left.convert_to<std::string>()
                                + "," + center.convert_to<std::string>()
                                + "," + right.convert_to<std::string>()
                                + "} lhs=" + lhs.convert_to<std::string>()
                                + " rhs=" + rhs.convert_to<std::string>();
                        }
                    }
                }
                for (int index = 0;
                     index + 1 < static_cast<int>(profile.size());
                     ++index) {
                    if (!dimension_normalized_decreasing(profile, index)) {
                        ++decrease_failures;
                        if (decrease_witness.empty()) {
                            decrease_witness =
                                "{q=" + std::to_string(q)
                                + ",power=" + std::to_string(step)
                                + ",index=" + std::to_string(index) + "}";
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_DIMENSION_NORMALIZED_PROFILE"
            << " maximum_q=" << maximum_q
            << " maximum_half_power=" << maximum_half_power
            << " profiles=" << profiles
            << " inequalities=" << inequalities
            << " log_concavity_failures=" << log_concavity_failures
            << " decrease_failures=" << decrease_failures
            << " log_concavity_witness="
            << (log_concavity_witness.empty() ? "{}" : log_concavity_witness)
            << " decrease_witness="
            << (decrease_witness.empty() ? "{}" : decrease_witness)
            << " result="
            << (
                log_concavity_failures == 0 && decrease_failures == 0
                ? "PASS_BOUNDED_EXACT_DIAGNOSTIC"
                : "COUNTEREXAMPLE"
            )
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
