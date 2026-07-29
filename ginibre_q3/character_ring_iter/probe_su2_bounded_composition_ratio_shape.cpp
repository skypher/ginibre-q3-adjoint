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

class Profile {
public:
    Profile(int q, int power)
        : radius_(q * power),
          values_(static_cast<std::size_t>(2 * radius_ + 1)) {
        std::vector<Integer> current{Integer{1}};
        int current_radius = 0;
        for (int step = 0; step < power; ++step) {
            const int next_radius = current_radius + q;
            std::vector<Integer> next(
                static_cast<std::size_t>(2 * next_radius + 1)
            );
            for (int index = -current_radius;
                 index <= current_radius;
                 ++index) {
                for (int increment = -q; increment <= q; ++increment) {
                    next[static_cast<std::size_t>(
                        index + increment + next_radius
                    )] += current[static_cast<std::size_t>(
                        index + current_radius
                    )];
                }
            }
            current = std::move(next);
            current_radius = next_radius;
        }
        values_ = std::move(current);
    }

    Integer at(int index) const {
        if (index < -radius_ || index > radius_) {
            return 0;
        }
        return values_[static_cast<std::size_t>(index + radius_)];
    }

    int radius() const {
        return radius_;
    }

private:
    int radius_;
    std::vector<Integer> values_;
};

Integer psi(const Profile& profile, int row, int index) {
    return profile.at(index - row) - profile.at(index + row + 1);
}

Integer adjacent_wronskian(
    const Profile& profile,
    int row,
    int index
) {
    return
        psi(profile, 0, index) * psi(profile, row, index + 1)
        - psi(profile, 0, index + 1) * psi(profile, row, index);
}

int sign(const Integer& value) {
    return value > 0 ? 1 : (value < 0 ? -1 : 0);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        int maximum_q = 8;
        int maximum_power = 12;
        if (argc >= 2) {
            maximum_q = parse_positive(argv[1], "maximum_q");
        }
        if (argc >= 3) {
            maximum_power = parse_positive(argv[2], "maximum_power");
        }
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_bounded_composition_ratio_shape"
                " [maximum_q] [maximum_power]"
            );
        }

        unsigned long long rows = 0;
        unsigned long long rows_with_multiple_turns = 0;
        unsigned long long negative_to_positive_turns = 0;
        unsigned long long positive_to_negative_turns = 0;
        int maximum_turns = 0;
        std::string witness;

        for (int q = 1; q <= maximum_q; ++q) {
            for (int power = 1; power <= maximum_power; ++power) {
                const Profile profile{q, power};
                const int maximum_row = 2 * profile.radius();
                const int maximum_index = profile.radius() + maximum_row;
                for (int row = 1; row <= maximum_row; ++row) {
                    ++rows;
                    int previous_sign = 0;
                    int turns = 0;
                    for (int index = 0; index < maximum_index; ++index) {
                        const int current_sign =
                            sign(adjacent_wronskian(profile, row, index));
                        if (
                            current_sign != 0
                            && previous_sign != 0
                            && current_sign != previous_sign
                        ) {
                            ++turns;
                            if (previous_sign < current_sign) {
                                ++negative_to_positive_turns;
                            } else {
                                ++positive_to_negative_turns;
                            }
                        }
                        if (current_sign != 0) {
                            previous_sign = current_sign;
                        }
                    }
                    maximum_turns = std::max(maximum_turns, turns);
                    if (turns > 1) {
                        ++rows_with_multiple_turns;
                        if (witness.empty()) {
                            witness =
                                "{q=" + std::to_string(q)
                                + ",power=" + std::to_string(power)
                                + ",row=" + std::to_string(row) + "}";
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_BOUNDED_COMPOSITION_RATIO_SHAPE"
            << " maximum_q=" << maximum_q
            << " maximum_power=" << maximum_power
            << " rows=" << rows
            << " rows_with_multiple_turns=" << rows_with_multiple_turns
            << " negative_to_positive_turns=" << negative_to_positive_turns
            << " positive_to_negative_turns=" << positive_to_negative_turns
            << " maximum_turns=" << maximum_turns
            << " witness=" << (witness.empty() ? "{}" : witness)
            << " result="
            << (
                rows_with_multiple_turns == 0
                    ? "NO_MULTIPLE_ADJACENT_RATIO_TURNS"
                    : "COUNTEREXAMPLE"
            )
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
