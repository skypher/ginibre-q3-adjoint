#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Matrix = std::vector<std::vector<int>>;
using IntegerMatrix = std::vector<std::vector<Integer>>;

int parse_positive(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value <= 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error("bound must be a positive integer");
    }
    return static_cast<int>(value);
}

bool fuses_half(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= source + label
        && source + target + label <= 2 * level;
}

bool fuses_standard(int level, int left, int right, int target) {
    return (left + right + target) % 2 == 0
        && std::abs(left - right) <= target
        && target <= left + right
        && left + right + target <= 2 * level;
}

IntegerMatrix multiply(const IntegerMatrix& left, const Matrix& right) {
    IntegerMatrix result(
        left.size(),
        std::vector<Integer>(right.size())
    );
    for (std::size_t source = 0U; source < left.size(); ++source) {
        for (std::size_t middle = 0U; middle < right.size(); ++middle) {
            if (left[source][middle] == 0) {
                continue;
            }
            for (std::size_t target = 0U; target < right.size(); ++target) {
                if (right[middle][target] != 0) {
                    result[source][target] +=
                        left[source][middle] * right[middle][target];
                }
            }
        }
    }
    return result;
}

IntegerMatrix identity(int size) {
    IntegerMatrix result(
        static_cast<std::size_t>(size),
        std::vector<Integer>(static_cast<std::size_t>(size))
    );
    for (int index = 0; index < size; ++index) {
        result[static_cast<std::size_t>(index)]
              [static_cast<std::size_t>(index)] = 1;
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error("usage: MAXIMUM_K");
        }
        const int maximum_k = parse_positive(argv[1]);
        std::uint64_t parameters = 0U;
        std::uint64_t block_entries = 0U;
        std::uint64_t power_entries = 0U;

        for (int k = 3; k <= maximum_k; ++k) {
            const int lower_level = k - 1;
            const int paired = (k + 1) / 2;
            for (int q = 1; 2 * q < k; ++q) {
                ++parameters;
                const int width = lower_level - 2 * q;
                if (width < 0) {
                    throw std::runtime_error("negative terminal width");
                }
                Matrix odd(
                    static_cast<std::size_t>(paired),
                    std::vector<int>(static_cast<std::size_t>(paired))
                );
                for (int source = 0; source < paired; ++source) {
                    for (int target = 0; target < paired; ++target) {
                        const int same = fuses_half(k, q, source, target)
                            ? 1 : 0;
                        const int crossed = fuses_half(
                            k,
                            q,
                            source,
                            k - target
                        ) ? 1 : 0;
                        odd[static_cast<std::size_t>(source)]
                           [static_cast<std::size_t>(target)] =
                            same - crossed;
                        const int expected = fuses_standard(
                            lower_level,
                            2 * q,
                            2 * source,
                            2 * target
                        ) ? 1 : 0;
                        ++block_entries;
                        if (
                            odd[static_cast<std::size_t>(source)]
                               [static_cast<std::size_t>(target)]
                            != expected
                        ) {
                            throw std::runtime_error(
                                "odd quotient coordinate mismatch"
                            );
                        }
                    }
                }

                Matrix terminal(
                    static_cast<std::size_t>(lower_level + 1),
                    std::vector<int>(
                        static_cast<std::size_t>(lower_level + 1)
                    )
                );
                for (int source = 0; source <= lower_level; ++source) {
                    for (int target = 0; target <= lower_level; ++target) {
                        terminal[static_cast<std::size_t>(source)]
                                [static_cast<std::size_t>(target)] =
                            fuses_standard(
                                lower_level,
                                width,
                                source,
                                target
                            ) ? 1 : 0;
                    }
                }
                IntegerMatrix odd_power = identity(paired);
                IntegerMatrix terminal_power = identity(lower_level + 1);
                for (int power = 1; power <= 4; ++power) {
                    odd_power = multiply(odd_power, odd);
                    terminal_power = multiply(terminal_power, terminal);
                    for (int source = 0; source < paired; ++source) {
                        for (int target = 0; target < paired; ++target) {
                            const int reflected_target = power % 2 == 0
                                ? 2 * target
                                : lower_level - 2 * target;
                            const Integer& expected = terminal_power[
                                static_cast<std::size_t>(2 * source)
                            ][static_cast<std::size_t>(reflected_target)];
                            ++power_entries;
                            if (
                                odd_power[
                                    static_cast<std::size_t>(source)
                                ][static_cast<std::size_t>(target)]
                                != expected
                            ) {
                                throw std::runtime_error(
                                    "terminal-coordinate power mismatch"
                                );
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_SHELL_TERMINAL_COORDINATES"
            << " maximum_K=" << maximum_k
            << " parameters=" << parameters
            << " block_entries=" << block_entries
            << " power_entries=" << power_entries
            << " powers=1..4"
            << " result=PASS_EXACT_IDENTITIES\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
