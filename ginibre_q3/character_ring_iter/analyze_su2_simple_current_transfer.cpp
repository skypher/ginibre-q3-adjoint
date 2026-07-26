#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
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
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(value);
}

bool fuses(int level, int first, int second, int output) {
    return std::abs(first - second) <= output
        && output <= std::min(
            first + second,
            2 * level - first - second
        )
        && ((first + second + output) & 1) == 0;
}

void multiply(int level, int label, std::vector<Integer>& state) {
    std::vector<Integer> next(static_cast<std::size_t>(level + 1));
    for (int source = 0; source <= level; ++source) {
        if (state[static_cast<std::size_t>(source)] == 0) {
            continue;
        }
        for (int output = 0; output <= level; ++output) {
            if (fuses(level, label, source, output)) {
                next[static_cast<std::size_t>(output)]
                    += state[static_cast<std::size_t>(source)];
            }
        }
    }
    state.swap(next);
}

struct Moments {
    std::vector<Integer> closed;
    std::vector<Integer> wall;
};

Moments moments(int level, int label, int maximum_power) {
    Moments result{
        std::vector<Integer>(
            static_cast<std::size_t>(maximum_power + 1)
        ),
        std::vector<Integer>(
            static_cast<std::size_t>(maximum_power + 1)
        )
    };
    std::vector<Integer> state(static_cast<std::size_t>(level + 1));
    state[0] = 1;
    for (int power = 0; power <= maximum_power; ++power) {
        result.closed[static_cast<std::size_t>(power)] = state[0];
        result.wall[static_cast<std::size_t>(power)] =
            state[static_cast<std::size_t>(level)];
        if (power != maximum_power) {
            multiply(level, label, state);
        }
    }
    return result;
}

long double endpoint_value(
    const Moments& data,
    const std::vector<long double>& elementary,
    int factors,
    long double& absolute_sum
) {
    long double value = 0.0L;
    absolute_sum = 0.0L;
    for (int count = 0; count <= factors; ++count) {
        const int left_even = 2 * factors - 2 * count + 2;
        const int left_odd = left_even - 1;
        const int right_even = 2 * count;
        const int right_odd = right_even + 1;
        const long double term = elementary[
            static_cast<std::size_t>(count)
        ] * (
            data.wall[static_cast<std::size_t>(left_even)]
                .convert_to<long double>()
            * data.closed[static_cast<std::size_t>(right_even)]
                .convert_to<long double>()
            - data.wall[static_cast<std::size_t>(left_odd)]
                .convert_to<long double>()
            * data.closed[static_cast<std::size_t>(right_odd)]
                .convert_to<long double>()
        );
        value += term;
        absolute_sum += std::abs(term);
    }
    return value;
}

bool check_order(
    int level,
    int label,
    int m,
    const Moments& data,
    const std::vector<long double>& slopes,
    const char* order,
    long long& prefixes
) {
    std::vector<long double> elementary(slopes.size() + 1U);
    elementary[0] = 1.0L;
    for (std::size_t index = 0; index < slopes.size(); ++index) {
        const long double slope = slopes[index];
        for (std::size_t degree = index + 1U;
             degree > 0U; --degree) {
            elementary[degree] +=
                slope * elementary[degree - 1U];
        }
        long double absolute_sum = 0.0L;
        const long double value = endpoint_value(
            data,
            elementary,
            static_cast<int>(index + 1U),
            absolute_sum
        );
        ++prefixes;
        const long double tolerance =
            1.0e-15L * std::max(1.0L, absolute_sum);
        if (value < -tolerance) {
            std::cout << std::setprecision(20)
                << "SIMPLE_CURRENT_TRANSFER_NEGATIVE"
                << " level=" << level
                << " label=" << label
                << " m=" << m
                << " order=" << order
                << " prefix=" << index + 1U
                << " value=" << value
                << " absolute_sum=" << absolute_sum << '\n';
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_simple_current_transfer "
                "MAXIMUM_LEVEL MAXIMUM_M"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_m = parse_positive(argv[2], "maximum m");
        const long double pi = std::acos(-1.0L);
        long long cases = 0;
        long long prefixes = 0;
        bool descending = true;
        bool ascending = true;
        for (int level = 4; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                const Moments data = moments(
                    level, label, 2 * maximum_m + 2
                );
                for (int m = 1; m <= maximum_m; ++m) {
                    std::vector<long double> slopes;
                    slopes.reserve(static_cast<std::size_t>(m));
                    for (int index = 1; index <= m; ++index) {
                        const long double angle =
                            static_cast<long double>(2 * index - 1) * pi
                            / static_cast<long double>(2 * (2 * m + 1));
                        const long double tangent = std::tan(angle);
                        slopes.push_back(
                            1.0L / (tangent * tangent)
                        );
                    }
                    ++cases;
                    if (descending
                        && !check_order(
                            level, label, m, data, slopes,
                            "descending", prefixes
                        )) {
                        descending = false;
                    }
                    std::reverse(slopes.begin(), slopes.end());
                    if (ascending
                        && !check_order(
                            level, label, m, data, slopes,
                            "ascending", prefixes
                        )) {
                        ascending = false;
                    }
                }
            }
        }
        std::cout
            << "SU2_SIMPLE_CURRENT_TRANSFER"
            << " maximum_level=" << maximum_level
            << " maximum_m=" << maximum_m
            << " cases=" << cases
            << " prefixes=" << prefixes
            << " descending_prefix_nonnegative="
            << (descending ? "true" : "false")
            << " ascending_prefix_nonnegative="
            << (ascending ? "true" : "false")
            << " result=PASS_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SIMPLE_CURRENT_TRANSFER FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
