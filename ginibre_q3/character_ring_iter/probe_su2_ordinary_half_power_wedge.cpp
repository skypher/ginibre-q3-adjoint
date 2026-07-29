#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

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

class SymmetricProfile {
public:
    SymmetricProfile(int q_half, int half_power)
        : radius_(q_half * half_power),
          values_(static_cast<std::size_t>(2 * radius_ + 1)) {
        std::vector<Integer> current{Integer{1}};
        int current_radius = 0;
        for (int step = 0; step < half_power; ++step) {
            const int next_radius = current_radius + q_half;
            std::vector<Integer> next(
                static_cast<std::size_t>(2 * next_radius + 1)
            );
            for (int exponent = -current_radius;
                 exponent <= current_radius;
                 ++exponent) {
                const Integer& coefficient =
                    current[static_cast<std::size_t>(
                        exponent + current_radius
                    )];
                for (int increment = -q_half;
                     increment <= q_half;
                     ++increment) {
                    next[static_cast<std::size_t>(
                        exponent + increment + next_radius
                    )] += coefficient;
                }
            }
            current = std::move(next);
            current_radius = next_radius;
        }
        values_ = std::move(current);
    }

    Integer at(int exponent) const {
        if (exponent < -radius_ || exponent > radius_) {
            return 0;
        }
        return values_[static_cast<std::size_t>(exponent + radius_)];
    }

    int radius() const {
        return radius_;
    }

private:
    int radius_;
    std::vector<Integer> values_;
};

Integer psi(const SymmetricProfile& profile, int row, int index) {
    return profile.at(index - row) - profile.at(index + row + 1);
}

Integer wronskian(
    const SymmetricProfile& profile,
    int row,
    int left,
    int right
) {
    return
        psi(profile, 0, left) * psi(profile, row, right)
        - psi(profile, 0, right) * psi(profile, row, left);
}

bool scan_case(
    int q_half,
    int half_power,
    int target,
    std::size_t& cases,
    std::size_t& tails,
    std::size_t& prefixes,
    bool& prefix_counterexample
) {
    const SymmetricProfile profile{q_half, half_power};
    const int max_index =
        profile.radius() + std::max(q_half, target);
    std::vector<Integer> by_min(
        static_cast<std::size_t>(max_index + 1)
    );
    std::vector<Integer> by_max(
        static_cast<std::size_t>(max_index + 1)
    );
    for (int left = 0; left <= max_index; ++left) {
        for (int right = left + 1; right <= max_index; ++right) {
            const Integer contribution =
                wronskian(profile, q_half, left, right)
                * wronskian(profile, target, left, right);
            by_min[static_cast<std::size_t>(left)] += contribution;
            by_max[static_cast<std::size_t>(right)] += contribution;
        }
    }

    Integer tail = 0;
    for (int cutoff = max_index; cutoff >= 0; --cutoff) {
        tail += by_min[static_cast<std::size_t>(cutoff)];
        ++tails;
        if (tail < 0) {
            std::cout
                << "SU2_ORDINARY_HALF_POWER_WEDGE"
                << " negative_tail"
                << " q_half=" << q_half
                << " half_power=" << half_power
                << " target=" << target
                << " cutoff=" << cutoff
                << " value=" << tail
                << '\n';
            return false;
        }
    }

    Integer prefix = 0;
    for (int cutoff = 0; cutoff <= max_index; ++cutoff) {
        prefix += by_max[static_cast<std::size_t>(cutoff)];
        ++prefixes;
        if (prefix < 0) {
            if (!prefix_counterexample) {
                std::cout
                    << "SU2_ORDINARY_HALF_POWER_WEDGE"
                    << " negative_prefix"
                    << " q_half=" << q_half
                    << " half_power=" << half_power
                    << " target=" << target
                    << " cutoff=" << cutoff
                    << " value=" << prefix
                    << '\n';
            }
            prefix_counterexample = true;
            break;
        }
    }

    ++cases;
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    int max_q_half = 12;
    int max_half_power = 14;
    try {
        if (argc >= 2) {
            max_q_half = parse_positive(argv[1], "max_q_half");
        }
        if (argc >= 3) {
            max_half_power =
                parse_positive(argv[2], "max_half_power");
        }
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_ordinary_half_power_wedge"
                " [max_q_half] [max_half_power]"
            );
        }
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }

    std::size_t cases = 0;
    std::size_t tails = 0;
    std::size_t prefixes = 0;
    bool prefix_counterexample = false;
    for (int q_half = 1; q_half <= max_q_half; ++q_half) {
        for (int half_power = 1;
             half_power <= max_half_power;
             ++half_power) {
            const int max_target = 2 * q_half * half_power;
            for (int target = 1; target <= max_target; ++target) {
                if (!scan_case(
                        q_half,
                        half_power,
                        target,
                        cases,
                        tails,
                        prefixes,
                        prefix_counterexample
                    )) {
                    return EXIT_SUCCESS;
                }
            }
        }
    }

    std::cout
        << "SU2_ORDINARY_HALF_POWER_WEDGE"
        << " cases=" << cases
        << " tails=" << tails
        << " prefixes=" << prefixes
        << " max_q_half=" << max_q_half
        << " max_half_power=" << max_half_power
        << " prefix_counterexample="
        << (prefix_counterexample ? 1 : 0)
        << " result=NO_NEGATIVE_WEDGE_TAIL"
        << '\n';
    return EXIT_SUCCESS;
}
