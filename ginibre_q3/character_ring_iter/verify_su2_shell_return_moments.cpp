#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

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

bool fuses(int level, int label, int source, int target) {
    return (source + label + target) % 2 == 0
        && std::abs(source - label) <= target
        && target <= source + label
        && source + label + target <= 2 * level;
}

std::vector<Integer> multiply_row(
    const std::vector<Integer>& row,
    int level,
    int label
) {
    std::vector<Integer> result(static_cast<std::size_t>(level + 1));
    for (int source = 0; source <= level; ++source) {
        if (row[static_cast<std::size_t>(source)] == 0) {
            continue;
        }
        for (int target = 0; target <= level; ++target) {
            if (fuses(level, label, source, target)) {
                result[static_cast<std::size_t>(target)] +=
                    row[static_cast<std::size_t>(source)];
            }
        }
    }
    return result;
}

Integer binomial_integer(int top, int bottom) {
    if (bottom < 0 || top < bottom) {
        return 0;
    }
    Integer result = 1;
    for (int index = 1; index <= bottom; ++index) {
        result *= top - bottom + index;
        result /= index;
    }
    return result;
}

Integer tensor_weight_five(int label_half, int depth) {
    Integer result = 0;
    for (int image = 0; image <= 5; ++image) {
        const Integer term = binomial_integer(5, image)
            * binomial_integer(
                depth - image * (2 * label_half + 1) + 4,
                4
            );
        if (image % 2 == 0) {
            result += term;
        } else {
            result -= term;
        }
    }
    return result;
}

Integer classical_five_multiplicity(int label_half, int depth) {
    return tensor_weight_five(label_half, depth)
        - tensor_weight_five(label_half, depth - 1);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error("usage: MAXIMUM_K");
        }
        const int maximum_k = parse_positive(argv[1]);
        std::uint64_t parameters = 0U;
        std::uint64_t direct_entries = 0U;
        std::uint64_t fold_entries = 0U;
        std::uint64_t active_corrections = 0U;

        for (int k = 3; k <= maximum_k; ++k) {
            const int level = 2 * k;
            for (int q = 1; 2 * q < k; ++q) {
                ++parameters;
                const int label = 2 * q;
                const int width = k - 1 - 2 * q;
                std::vector<Integer> row(
                    static_cast<std::size_t>(level + 1)
                );
                row[0] = 1;
                Integer f4_direct = 0;
                Integer f5_direct = 0;
                for (int power = 1; power <= 5; ++power) {
                    row = multiply_row(row, level, label);
                    if (power == 4) {
                        f4_direct = row[0];
                    } else if (power == 5) {
                        f5_direct = row[0];
                    }
                }
                direct_entries += 2U;

                const Integer f4_expected = 2 * q + 1;
                const Integer f5_stable = (
                    5 * Integer(q) * q + 5 * q + 2
                ) / 2;
                const Integer correction = binomial_integer(
                    q - 2 * width - 1,
                    2
                );
                const Integer f5_expected = f5_stable - correction;
                if (correction != 0) {
                    ++active_corrections;
                }
                if (f4_direct != f4_expected || f5_direct != f5_expected) {
                    throw std::runtime_error(
                        "return-moment mismatch K=" + std::to_string(k)
                        + " Q=" + std::to_string(q)
                    );
                }

                const int negative_depth = q - 2 * width - 3;
                const int positive_depth = negative_depth - 1;
                const Integer folded_five = classical_five_multiplicity(
                    q,
                    5 * q
                ) - classical_five_multiplicity(q, negative_depth)
                    + classical_five_multiplicity(q, positive_depth);
                ++fold_entries;
                if (folded_five != f5_direct) {
                    throw std::runtime_error(
                        "Kac-Walton return-fold mismatch K="
                        + std::to_string(k)
                        + " Q=" + std::to_string(q)
                    );
                }
            }
        }

        std::cout
            << "SU2_SHELL_RETURN_MOMENTS"
            << " maximum_K=" << maximum_k
            << " parameters=" << parameters
            << " direct_entries=" << direct_entries
            << " Kac_Walton_entries=" << fold_entries
            << " active_f5_corrections=" << active_corrections
            << " powers=4,5"
            << " result=PASS_EXACT_IDENTITIES\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
