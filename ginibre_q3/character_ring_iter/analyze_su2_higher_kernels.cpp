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
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(value);
}

bool fuses(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= std::min(
            source + label,
            2 * level - source - label
        )
        && ((source + label + target) & 1) == 0;
}

std::vector<Integer> multiply(
    int level,
    int label,
    const std::vector<Integer>& state
) {
    std::vector<Integer> next(static_cast<std::size_t>(level + 1));
    for (int source = 0; source <= level; ++source) {
        const Integer& coefficient =
            state[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        for (int target = 0; target <= level; ++target) {
            if (fuses(level, label, source, target)) {
                next[static_cast<std::size_t>(target)] += coefficient;
            }
        }
    }
    return next;
}

struct Result {
    bool initialized = false;
    bool failed = false;
    Integer minimum = 0;
    int level = 0;
    int label = 0;
    int vertex = 0;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_higher_kernels "
                "MAXIMUM_LEVEL MAXIMUM_S"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_s = parse_positive(argv[2], "maximum s");
        if (maximum_s < 2) {
            throw std::runtime_error("maximum s must be at least two");
        }
        std::vector<Result> results(
            static_cast<std::size_t>(maximum_s + 1)
        );
        std::vector<Result> dual_results(
            static_cast<std::size_t>(maximum_s + 1)
        );
        std::uint64_t parameters = 0U;
        std::uint64_t coordinates = 0U;
        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameters;
                std::vector<std::vector<Integer>> powers(
                    static_cast<std::size_t>(2 * maximum_s + 2)
                );
                powers[0].assign(
                    static_cast<std::size_t>(level + 1),
                    Integer(0)
                );
                powers[0][0] = 1;
                for (int power = 1;
                     power <= 2 * maximum_s + 1;
                     ++power) {
                    powers[static_cast<std::size_t>(power)] =
                        multiply(
                            level,
                            label,
                            powers[static_cast<std::size_t>(power - 1)]
                        );
                }
                for (int s = 2; s <= maximum_s; ++s) {
                    const int even_power = 2 * s;
                    const int odd_power = even_power + 1;
                    const Integer& f_even =
                        powers[static_cast<std::size_t>(even_power)][0];
                    const Integer& f_odd =
                        powers[static_cast<std::size_t>(odd_power)][0];
                    const Integer& g_even =
                        powers[static_cast<std::size_t>(even_power)][
                            static_cast<std::size_t>(level)
                        ];
                    const Integer& g_previous =
                        powers[static_cast<std::size_t>(even_power - 1)][
                            static_cast<std::size_t>(level)
                        ];
                    Result& result =
                        results[static_cast<std::size_t>(s)];
                    Result& dual_result =
                        dual_results[static_cast<std::size_t>(s)];
                    for (int vertex = 0; vertex <= level; vertex += 2) {
                        coordinates += 2U;
                        const Integer margin =
                            f_even
                                * powers[
                                    static_cast<std::size_t>(odd_power)
                                ][static_cast<std::size_t>(vertex)]
                            - f_odd
                                * powers[
                                    static_cast<std::size_t>(even_power)
                                ][static_cast<std::size_t>(vertex)];
                        if (!result.initialized
                            || margin < result.minimum) {
                            result.initialized = true;
                            result.minimum = margin;
                            result.level = level;
                            result.label = label;
                            result.vertex = vertex;
                        }
                        if (margin < 0) {
                            result.failed = true;
                        }
                        const Integer dual_margin =
                            g_even
                                * powers[
                                    static_cast<std::size_t>(even_power)
                                ][static_cast<std::size_t>(vertex)]
                            - g_previous
                                * powers[
                                    static_cast<std::size_t>(odd_power)
                                ][static_cast<std::size_t>(vertex)];
                        if (!dual_result.initialized
                            || dual_margin < dual_result.minimum) {
                            dual_result.initialized = true;
                            dual_result.minimum = dual_margin;
                            dual_result.level = level;
                            dual_result.label = label;
                            dual_result.vertex = vertex;
                        }
                        if (dual_margin < 0) {
                            dual_result.failed = true;
                        }
                    }
                }
            }
        }
        std::cout
            << "SU2_HIGHER_KERNELS"
            << " maximum_level=" << maximum_level
            << " maximum_s=" << maximum_s
            << " parameters=" << parameters
            << " coordinates=" << coordinates
            << " result=PASS_DISCOVERY\n";
        for (int s = 2; s <= maximum_s; ++s) {
            const Result& result = results[static_cast<std::size_t>(s)];
            std::cout
                << "KERNEL s=" << s
                << " minimum=" << result.minimum
                << " witness=(" << result.level
                << ',' << result.label
                << ',' << result.vertex << ')'
                << " result="
                << (result.failed ? "FAIL_POINTWISE" : "PASS_POINTWISE")
                << '\n';
            const Result& dual =
                dual_results[static_cast<std::size_t>(s)];
            std::cout
                << "DUAL_KERNEL p=" << s
                << " minimum=" << dual.minimum
                << " witness=(" << dual.level
                << ',' << dual.label
                << ',' << dual.vertex << ')'
                << " result="
                << (dual.failed ? "FAIL_POINTWISE" : "PASS_POINTWISE")
                << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_HIGHER_KERNELS FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
