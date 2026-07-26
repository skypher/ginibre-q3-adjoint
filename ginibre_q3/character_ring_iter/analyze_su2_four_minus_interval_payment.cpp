#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Interval {
    int lower = 0;
    int upper = -2;

    auto operator<=>(const Interval&) const = default;
};

struct Witness {
    bool initialized = false;
    std::int64_t margin = 0;
    int level = 0;
    Interval first;
    Interval second;
    std::array<std::int64_t, 4> profile{};
};

int parse_positive(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error("maximum level must be positive");
    }
    return static_cast<int>(value);
}

bool fuses(int level, int first, int second, int output) {
    return std::abs(first - second) <= output
        && output <= std::min(
            first + second, 2 * level - first - second
        )
        && ((first + second + output) % 2) == 0;
}

std::int64_t fourfold_output(
    int level,
    const Interval& first,
    const Interval& second,
    int output
) {
    std::int64_t result = 0;
    for (int left = first.lower; left <= first.upper; left += 2) {
        for (int right = second.lower; right <= second.upper; right += 2) {
            result += fuses(level, left, right, output) ? 1 : 0;
        }
    }
    return result;
}

void print_interval(const Interval& interval) {
    std::cout << '[' << interval.lower << ',' << interval.upper << "]_2";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_four_minus_interval_payment maximum_level"
            );
        }
        const int maximum_level = parse_positive(argv[1]);
        std::array<Witness, 6> witnesses{};
        std::uint64_t interval_pairs = 0;

        for (int level = 4; level <= maximum_level; ++level) {
            std::set<Interval> interval_set;
            for (int first = 2; first <= level; first += 2) {
                if (std::min(first, level - first) < 2) {
                    continue;
                }
                for (int second = first; second <= level; second += 2) {
                    if (std::min(second, level - second) < 2) {
                        continue;
                    }
                    interval_set.insert(Interval{
                        std::abs(first - second),
                        std::min(
                            first + second, 2 * level - first - second
                        )
                    });
                }
            }
            const std::vector<Interval> intervals(
                interval_set.begin(), interval_set.end()
            );
            for (const Interval& first : intervals) {
                for (const Interval& second : intervals) {
                    const int lower = std::max(first.lower, second.lower);
                    const int upper = std::min(first.upper, second.upper);
                    if (lower <= 0 || lower > upper) {
                        continue;
                    }
                    const int d = (upper - lower) / 2 + 1;
                    if (d < 1 || d > 5) {
                        continue;
                    }
                    ++interval_pairs;
                    const std::array<std::int64_t, 4> profile{
                        fourfold_output(level, first, second, 0),
                        fourfold_output(level, first, second, 2),
                        fourfold_output(level, first, second, 4),
                        fourfold_output(level, first, second, 6)
                    };
                    if (profile[0] != d) {
                        throw std::runtime_error(
                            "zero-output profile disagrees with intersection"
                        );
                    }
                    const std::int64_t payment =
                        profile[0] + 2 * profile[1]
                        + 2 * profile[2] + profile[3];
                    const std::int64_t margin = payment - 16 * d;
                    Witness& witness =
                        witnesses[static_cast<std::size_t>(d)];
                    if (!witness.initialized || margin < witness.margin) {
                        witness = Witness{
                            true, margin, level, first, second, profile
                        };
                    }
                }
            }
        }

        std::cout << "SU2_FOUR_MINUS_INTERVAL_PAYMENT maximum_level="
                  << maximum_level << " interval_pairs=" << interval_pairs
                  << '\n';
        for (int d = 1; d <= 5; ++d) {
            const Witness& witness =
                witnesses[static_cast<std::size_t>(d)];
            std::cout << "d=" << d;
            if (!witness.initialized) {
                std::cout << " no-cases\n";
                continue;
            }
            std::cout << " minimum_margin=" << witness.margin
                      << " level=" << witness.level << " intervals=";
            print_interval(witness.first);
            std::cout << ',';
            print_interval(witness.second);
            std::cout << " profile=("
                      << witness.profile[0] << ','
                      << witness.profile[1] << ','
                      << witness.profile[2] << ','
                      << witness.profile[3] << ")\n";
        }
        std::cout << "SU2_FOUR_MINUS_INTERVAL_PAYMENT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_FOUR_MINUS_INTERVAL_PAYMENT FAILURE: "
                  << error.what() << '\n';
        return 1;
    }
}
