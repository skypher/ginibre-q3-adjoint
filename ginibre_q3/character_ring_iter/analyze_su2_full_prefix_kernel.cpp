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

Integer binomial(int n, int r) {
    if (r < 0 || r > n) {
        return 0;
    }
    r = std::min(r, n - r);
    Integer result = 1;
    for (int index = 1; index <= r; ++index) {
        result *= n - r + index;
        result /= index;
    }
    return result;
}

struct Minimum {
    bool initialized = false;
    Integer value = 0;
    int level = 0;
    int label = 0;
    int prefix = 0;
    int truncation = 0;
    int vertex = 0;

    void observe(
        const Integer& candidate,
        int candidate_level,
        int candidate_label,
        int candidate_prefix,
        int candidate_truncation,
        int candidate_vertex
    ) {
        if (!initialized || candidate < value) {
            initialized = true;
            value = candidate;
            level = candidate_level;
            label = candidate_label;
            prefix = candidate_prefix;
            truncation = candidate_truncation;
            vertex = candidate_vertex;
        }
    }
};

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 6 && std::string(argv[1]) == "--case") {
            const int level = parse_positive(argv[2], "level");
            const int label = parse_positive(argv[3], "label");
            const int prefix = parse_positive(argv[4], "prefix");
            const int truncation =
                parse_positive(argv[5], "truncation");
            if ((level & 1) != 0 || (label & 1) != 0
                || 2 * label >= level || prefix < 4
                || truncation < 2 || truncation >= prefix) {
                throw std::runtime_error(
                    "case parameters are outside the FBPC range"
                );
            }
            const int total_length = 2 * prefix + 2;
            std::vector<std::vector<Integer>> powers(
                static_cast<std::size_t>(total_length + 1)
            );
            powers[0].assign(
                static_cast<std::size_t>(level + 1),
                Integer(0)
            );
            powers[0][0] = 1;
            for (int power = 1;
                 power <= total_length;
                 ++power) {
                powers[static_cast<std::size_t>(power)] =
                    multiply(
                        level,
                        label,
                        powers[static_cast<std::size_t>(power - 1)]
                    );
            }
            std::vector<Integer> kernel(
                static_cast<std::size_t>(level + 1)
            );
            for (int index = 0; index <= truncation; ++index) {
                const int even_power = 2 * index;
                const int odd_power = even_power + 1;
                const Integer weight =
                    binomial(total_length - 1, even_power);
                const Integer& even_return =
                    powers[static_cast<std::size_t>(even_power)][0];
                const Integer& odd_return =
                    powers[static_cast<std::size_t>(odd_power)][0];
                for (int vertex = 0; vertex <= level; vertex += 2) {
                    kernel[static_cast<std::size_t>(vertex)] +=
                        weight
                        * (
                            even_return
                                * powers[
                                    static_cast<std::size_t>(
                                        total_length - even_power
                                    )
                                ][static_cast<std::size_t>(vertex)]
                            - odd_return
                                * powers[
                                    static_cast<std::size_t>(
                                        total_length - odd_power
                                    )
                                ][static_cast<std::size_t>(vertex)]
                        );
                }
            }
            std::cout
                << "SU2_FULL_PREFIX_KERNEL_CASE"
                << " level=" << level
                << " label=" << label
                << " prefix=" << prefix
                << " truncation=" << truncation
                << " values=";
            for (int vertex = 0; vertex <= level; vertex += 2) {
                std::cout
                    << (vertex == 0 ? "" : ",")
                    << vertex << ':'
                    << kernel[static_cast<std::size_t>(vertex)];
            }
            std::cout << '\n';
            return EXIT_SUCCESS;
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_full_prefix_kernel "
                "MAXIMUM_LEVEL MAXIMUM_PREFIX | "
                "--case LEVEL LABEL PREFIX TRUNCATION"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_prefix =
            parse_positive(argv[2], "maximum prefix");
        if (maximum_prefix < 4) {
            throw std::runtime_error(
                "maximum prefix must be at least four"
            );
        }

        Minimum global_minimum;
        Minimum endpoint_minimum;
        std::uint64_t parameter_rows = 0U;
        std::uint64_t kernel_rows = 0U;
        std::uint64_t coordinates = 0U;
        std::uint64_t negative_coordinates = 0U;
        std::uint64_t negative_endpoint_rows = 0U;
        bool printed_first_negative = false;

        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameter_rows;
                std::vector<std::vector<Integer>> powers(
                    static_cast<std::size_t>(2 * maximum_prefix + 3)
                );
                powers[0].assign(
                    static_cast<std::size_t>(level + 1),
                    Integer(0)
                );
                powers[0][0] = 1;
                for (int power = 1;
                     power <= 2 * maximum_prefix + 2;
                     ++power) {
                    powers[static_cast<std::size_t>(power)] =
                        multiply(
                            level,
                            label,
                            powers[static_cast<std::size_t>(power - 1)]
                        );
                }

                for (int prefix = 4;
                     prefix <= maximum_prefix;
                     ++prefix) {
                    const int total_length = 2 * prefix + 2;
                    std::vector<Integer> kernel(
                        static_cast<std::size_t>(level + 1)
                    );
                    for (int truncation = 0;
                         truncation <= prefix - 1;
                         ++truncation) {
                        const int even_power = 2 * truncation;
                        const int odd_power = even_power + 1;
                        const Integer weight =
                            binomial(total_length - 1, even_power);
                        const Integer& even_return =
                            powers[
                                static_cast<std::size_t>(even_power)
                            ][0];
                        const Integer& odd_return =
                            powers[
                                static_cast<std::size_t>(odd_power)
                            ][0];
                        for (int vertex = 0;
                             vertex <= level;
                             vertex += 2) {
                            kernel[static_cast<std::size_t>(vertex)] +=
                                weight
                                * (
                                    even_return
                                        * powers[
                                            static_cast<std::size_t>(
                                                total_length - even_power
                                            )
                                        ][static_cast<std::size_t>(vertex)]
                                    - odd_return
                                        * powers[
                                            static_cast<std::size_t>(
                                                total_length - odd_power
                                            )
                                        ][static_cast<std::size_t>(vertex)]
                                );
                        }
                        if (truncation < 2) {
                            continue;
                        }
                        ++kernel_rows;
                        const Integer& endpoint =
                            kernel[static_cast<std::size_t>(level)];
                        endpoint_minimum.observe(
                            endpoint,
                            level,
                            label,
                            prefix,
                            truncation,
                            level
                        );
                        if (endpoint < 0) {
                            ++negative_endpoint_rows;
                        }
                        for (int vertex = 0;
                             vertex <= level;
                             vertex += 2) {
                            ++coordinates;
                            const Integer& value =
                                kernel[static_cast<std::size_t>(vertex)];
                            global_minimum.observe(
                                value,
                                level,
                                label,
                                prefix,
                                truncation,
                                vertex
                            );
                            if (value >= 0) {
                                continue;
                            }
                            ++negative_coordinates;
                            if (!printed_first_negative) {
                                std::cout
                                    << "SU2_FULL_PREFIX_KERNEL"
                                    << " first_negative"
                                    << " level=" << level
                                    << " label=" << label
                                    << " prefix=" << prefix
                                    << " truncation=" << truncation
                                    << " vertex=" << vertex
                                    << " value=" << value << '\n';
                                printed_first_negative = true;
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_FULL_PREFIX_KERNEL"
            << " maximum_level=" << maximum_level
            << " maximum_prefix=" << maximum_prefix
            << " parameter_rows=" << parameter_rows
            << " kernel_rows=" << kernel_rows
            << " coordinates=" << coordinates
            << " negative_coordinates=" << negative_coordinates
            << " negative_endpoint_rows=" << negative_endpoint_rows
            << " minimum=" << global_minimum.value
            << " witness=("
            << global_minimum.level << ','
            << global_minimum.label << ','
            << global_minimum.prefix << ','
            << global_minimum.truncation << ','
            << global_minimum.vertex << ')'
            << " endpoint_minimum=" << endpoint_minimum.value
            << " endpoint_witness=("
            << endpoint_minimum.level << ','
            << endpoint_minimum.label << ','
            << endpoint_minimum.prefix << ','
            << endpoint_minimum.truncation << ')'
            << " result="
            << (
                negative_coordinates == 0U
                    ? "PASS_POINTWISE_DISCOVERY"
                    : "FAIL_POINTWISE"
            )
            << '\n';
        return negative_endpoint_rows == 0U
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_FULL_PREFIX_KERNEL FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
