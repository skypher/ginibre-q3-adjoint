#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Vector = std::vector<cpp_int>;

struct Failure {
    int level = -1;
    int left = -1;
    int right = -1;
    cpp_int value = 0;
    Vector profile;
    Vector square;
};

int parse_positive(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0LL
        || parsed > static_cast<long long>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    }
    return static_cast<int>(parsed);
}

Vector fusion_step(const Vector& current, int level, int factor) {
    Vector next(static_cast<std::size_t>(level + 1), 0);
    for (int source = 0; source <= level; ++source) {
        const cpp_int& coefficient = current[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        const int lower = std::abs(source - factor);
        const int upper = std::min(source + factor, 2 * level - source - factor);
        for (int target = lower; target <= upper; ++target) {
            next[static_cast<std::size_t>(target)] += coefficient;
        }
    }
    return next;
}

Vector fusion_product(const Vector& left, const Vector& right, int level) {
    Vector result(static_cast<std::size_t>(level + 1), 0);
    for (int factor = 0; factor <= level; ++factor) {
        const cpp_int& coefficient = left[static_cast<std::size_t>(factor)];
        if (coefficient == 0) {
            continue;
        }
        const Vector transformed = fusion_step(right, level, factor);
        for (int target = 0; target <= level; ++target) {
            result[static_cast<std::size_t>(target)]
                += coefficient * transformed[static_cast<std::size_t>(target)];
        }
    }
    return result;
}

bool log_concave(const Vector& values) {
    bool support_started = false;
    bool support_ended = false;
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (values[index] != 0) {
            if (support_ended) {
                return false;
            }
            support_started = true;
        } else if (support_started) {
            support_ended = true;
        }
        if (index > 0U && index + 1U < values.size()) {
            if (values[index] * values[index]
                < values[index - 1U] * values[index + 1U]) {
                return false;
            }
        }
    }
    return true;
}

void print_vector(const Vector& values) {
    std::cout << '[';
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index > 0U) {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << ']';
}

void print_failure(const char* label, const Failure& failure) {
    std::cout
        << ' ' << label << "_level=" << 2 * failure.level
        << ' ' << label << "_left=" << 2 * failure.left
        << ' ' << label << "_right=" << 2 * failure.right
        << ' ' << label << "_value=" << failure.value
        << ' ' << label << "_profile=";
    print_vector(failure.profile);
    std::cout << ' ' << label << "_square=";
    print_vector(failure.square);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? parse_positive(argv[1], "maximum_half_level")
            : 6;
        const int maximum_coefficient = argc >= 3
            ? parse_positive(argv[2], "maximum_coefficient")
            : 3;
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_finite_square_profile_cone "
                "[maximum_half_level [maximum_coefficient]]");
        }

        std::uint64_t profiles = 0U;
        std::uint64_t log_concave_squares = 0U;
        std::uint64_t current_checks = 0U;
        std::uint64_t negative_currents = 0U;
        std::uint64_t open_checks = 0U;
        std::uint64_t negative_open_currents = 0U;
        std::uint64_t log_concave_current_checks = 0U;
        std::uint64_t negative_log_concave_currents = 0U;
        std::uint64_t open_log_concave_checks = 0U;
        std::uint64_t negative_open_log_concave_currents = 0U;
        std::uint64_t boundary_admissible_squares = 0U;
        std::uint64_t boundary_admissible_checks = 0U;
        std::uint64_t negative_boundary_admissible_currents = 0U;
        std::uint64_t doubly_log_concave_profiles = 0U;
        std::uint64_t doubly_log_concave_checks = 0U;
        std::uint64_t negative_doubly_log_concave_currents = 0U;
        std::uint64_t profile_log_concave_boundary_checks = 0U;
        std::uint64_t negative_profile_log_concave_boundary_currents = 0U;
        std::uint64_t profile_log_concave_reflection_checks = 0U;
        std::uint64_t negative_profile_log_concave_reflection_margins = 0U;
        std::uint64_t profile_log_concave_boundary_sequences = 0U;
        std::uint64_t profile_log_concave_boundary_recrossings = 0U;
        Failure first_failure;
        Failure first_open_failure;
        Failure first_log_concave_failure;
        Failure first_open_log_concave_failure;
        Failure first_boundary_admissible_failure;
        Failure first_doubly_log_concave_failure;
        Failure first_profile_log_concave_boundary_failure;
        Failure first_profile_log_concave_reflection_failure;
        Failure first_profile_log_concave_boundary_recrossing;

        for (int level = 1; level <= maximum_level; ++level) {
            std::uint64_t profile_count = 1U;
            for (int index = 0; index <= level; ++index) {
                profile_count *= static_cast<std::uint64_t>(
                    maximum_coefficient + 1);
            }
            for (std::uint64_t code = 1U; code < profile_count; ++code) {
                std::uint64_t remaining = code;
                Vector profile(static_cast<std::size_t>(level + 1), 0);
                for (int index = 0; index <= level; ++index) {
                    profile[static_cast<std::size_t>(index)]
                        = remaining
                            % static_cast<std::uint64_t>(
                                maximum_coefficient + 1);
                    remaining /= static_cast<std::uint64_t>(
                        maximum_coefficient + 1);
                }
                const Vector square
                    = fusion_product(profile, profile, level);
                const bool profile_is_log_concave = log_concave(profile);
                const bool square_is_log_concave = log_concave(square);
                const bool doubly_log_concave
                    = profile_is_log_concave && square_is_log_concave;
                ++profiles;
                if (square_is_log_concave) {
                    ++log_concave_squares;
                }

                std::vector<Vector> transforms;
                transforms.reserve(static_cast<std::size_t>(level + 1));
                for (int radius = 0; radius <= level; ++radius) {
                    transforms.push_back(
                        fusion_step(square, level, radius));
                }

                bool boundary_admissible = square_is_log_concave;
                for (int left = 1; left <= level; ++left) {
                    const cpp_int boundary_current
                        = square[0]
                            * transforms[static_cast<std::size_t>(left)][
                                static_cast<std::size_t>(level)]
                        - square[static_cast<std::size_t>(left)]
                            * square[static_cast<std::size_t>(level)];
                    if (boundary_current < 0) {
                        boundary_admissible = false;
                        break;
                    }
                }
                if (boundary_admissible) {
                    ++boundary_admissible_squares;
                }
                if (doubly_log_concave) {
                    ++doubly_log_concave_profiles;
                }
                if (profile_is_log_concave) {
                    ++profile_log_concave_boundary_sequences;
                    cpp_int previous_boundary = 0;
                    bool negative_increment_seen = false;
                    for (int index = 1; index <= level; ++index) {
                        ++profile_log_concave_reflection_checks;
                        const cpp_int margin
                            = square[static_cast<std::size_t>(index - 1)]
                                * square[static_cast<std::size_t>(
                                    level - index)]
                            - square[static_cast<std::size_t>(index)]
                                * square[static_cast<std::size_t>(
                                    level - index + 1)];
                        if (margin < 0) {
                            ++negative_profile_log_concave_reflection_margins;
                            if (
                                first_profile_log_concave_reflection_failure
                                    .level
                                < 0) {
                                first_profile_log_concave_reflection_failure = {
                                    level,
                                    index,
                                    level - index + 1,
                                    margin,
                                    profile,
                                    square};
                            }
                        }
                        const cpp_int boundary
                            = square[0]
                                * square[static_cast<std::size_t>(
                                    level - index)]
                            - square[static_cast<std::size_t>(index)]
                                * square[static_cast<std::size_t>(level)];
                        const cpp_int increment
                            = boundary - previous_boundary;
                        if (increment < 0) {
                            negative_increment_seen = true;
                        } else if (
                            increment > 0 && negative_increment_seen) {
                            ++profile_log_concave_boundary_recrossings;
                            if (
                                first_profile_log_concave_boundary_recrossing
                                    .level
                                < 0) {
                                first_profile_log_concave_boundary_recrossing = {
                                    level,
                                    index,
                                    level,
                                    increment,
                                    profile,
                                    square};
                            }
                            negative_increment_seen = false;
                        }
                        previous_boundary = boundary;
                    }
                }

                for (int left = 1; left <= level; ++left) {
                    for (int right = 1; right <= level; ++right) {
                        const cpp_int current
                            = square[0]
                                * transforms[
                                    static_cast<std::size_t>(left)][
                                    static_cast<std::size_t>(right)]
                            - square[static_cast<std::size_t>(left)]
                                * square[static_cast<std::size_t>(right)];
                        ++current_checks;
                        if (left + right <= level) {
                            ++open_checks;
                        }
                        if (current < 0) {
                            ++negative_currents;
                            if (first_failure.level < 0) {
                                first_failure = {
                                    level,
                                    left,
                                    right,
                                    current,
                                    profile,
                                    square};
                            }
                            if (left + right <= level) {
                                ++negative_open_currents;
                                if (first_open_failure.level < 0) {
                                    first_open_failure = {
                                        level,
                                        left,
                                        right,
                                        current,
                                        profile,
                                        square};
                                }
                            }
                        }
                        if (square_is_log_concave) {
                            ++log_concave_current_checks;
                            if (left + right <= level) {
                                ++open_log_concave_checks;
                            }
                            if (current < 0) {
                                ++negative_log_concave_currents;
                                if (first_log_concave_failure.level < 0) {
                                    first_log_concave_failure = {
                                        level,
                                        left,
                                        right,
                                        current,
                                        profile,
                                        square};
                                }
                                if (left + right <= level) {
                                    ++negative_open_log_concave_currents;
                                    if (first_open_log_concave_failure.level < 0) {
                                        first_open_log_concave_failure = {
                                            level,
                                            left,
                                            right,
                                            current,
                                            profile,
                                            square};
                                    }
                                }
                            }
                        }
                        if (boundary_admissible) {
                            ++boundary_admissible_checks;
                            if (current < 0) {
                                ++negative_boundary_admissible_currents;
                                if (first_boundary_admissible_failure.level < 0) {
                                    first_boundary_admissible_failure = {
                                        level,
                                        left,
                                        right,
                                        current,
                                        profile,
                                        square};
                                }
                            }
                        }
                        if (doubly_log_concave) {
                            ++doubly_log_concave_checks;
                            if (current < 0) {
                                ++negative_doubly_log_concave_currents;
                                if (first_doubly_log_concave_failure.level < 0) {
                                    first_doubly_log_concave_failure = {
                                        level,
                                        left,
                                        right,
                                        current,
                                        profile,
                                        square};
                                }
                            }
                        }
                        if (profile_is_log_concave && right == level) {
                            ++profile_log_concave_boundary_checks;
                            if (current < 0) {
                                ++negative_profile_log_concave_boundary_currents;
                                if (
                                    first_profile_log_concave_boundary_failure
                                        .level
                                    < 0) {
                                    first_profile_log_concave_boundary_failure = {
                                        level,
                                        left,
                                        right,
                                        current,
                                        profile,
                                        square};
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_FINITE_SQUARE_PROFILE_CONE"
            << " maximum_level=" << 2 * maximum_level
            << " maximum_coefficient=" << maximum_coefficient
            << " profiles=" << profiles
            << " log_concave_squares=" << log_concave_squares
            << " current_checks=" << current_checks
            << " negative_currents=" << negative_currents
            << " open_checks=" << open_checks
            << " negative_open_currents=" << negative_open_currents
            << " log_concave_current_checks="
            << log_concave_current_checks
            << " negative_log_concave_currents="
            << negative_log_concave_currents
            << " open_log_concave_checks=" << open_log_concave_checks
            << " negative_open_log_concave_currents="
            << negative_open_log_concave_currents
            << " boundary_admissible_squares="
            << boundary_admissible_squares
            << " boundary_admissible_checks="
            << boundary_admissible_checks
            << " negative_boundary_admissible_currents="
            << negative_boundary_admissible_currents
            << " doubly_log_concave_profiles="
            << doubly_log_concave_profiles
            << " doubly_log_concave_checks="
            << doubly_log_concave_checks
            << " negative_doubly_log_concave_currents="
            << negative_doubly_log_concave_currents
            << " profile_log_concave_boundary_checks="
            << profile_log_concave_boundary_checks
            << " negative_profile_log_concave_boundary_currents="
            << negative_profile_log_concave_boundary_currents
            << " profile_log_concave_reflection_checks="
            << profile_log_concave_reflection_checks
            << " negative_profile_log_concave_reflection_margins="
            << negative_profile_log_concave_reflection_margins
            << " profile_log_concave_boundary_sequences="
            << profile_log_concave_boundary_sequences
            << " profile_log_concave_boundary_recrossings="
            << profile_log_concave_boundary_recrossings;
        if (first_failure.level >= 0) {
            print_failure("first_failure", first_failure);
        }
        if (first_open_failure.level >= 0) {
            print_failure("first_open_failure", first_open_failure);
        }
        if (first_log_concave_failure.level >= 0) {
            print_failure(
                "first_log_concave_failure",
                first_log_concave_failure);
        }
        if (first_boundary_admissible_failure.level >= 0) {
            print_failure(
                "first_boundary_admissible_failure",
                first_boundary_admissible_failure);
        }
        if (first_open_log_concave_failure.level >= 0) {
            print_failure(
                "first_open_log_concave_failure",
                first_open_log_concave_failure);
        }
        if (first_doubly_log_concave_failure.level >= 0) {
            print_failure(
                "first_doubly_log_concave_failure",
                first_doubly_log_concave_failure);
        }
        if (first_profile_log_concave_boundary_failure.level >= 0) {
            print_failure(
                "first_profile_log_concave_boundary_failure",
                first_profile_log_concave_boundary_failure);
        }
        if (first_profile_log_concave_reflection_failure.level >= 0) {
            print_failure(
                "first_profile_log_concave_reflection_failure",
                first_profile_log_concave_reflection_failure);
        }
        if (first_profile_log_concave_boundary_recrossing.level >= 0) {
            print_failure(
                "first_profile_log_concave_boundary_recrossing",
                first_profile_log_concave_boundary_recrossing);
        }
        std::cout
            << " result="
            << (negative_currents == 0U
                    ? "NO_NEGATIVE_SQUARE_PROFILE_CURRENT"
                    : "NEGATIVE_SQUARE_PROFILE_CURRENT")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
