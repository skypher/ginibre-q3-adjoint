#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Vector = std::vector<cpp_int>;
using Matrix = std::vector<Vector>;

int parse_positive(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0LL) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

Matrix fusion_matrix(int level, int factor) {
    Matrix result(
        static_cast<std::size_t>(level + 1),
        Vector(static_cast<std::size_t>(level + 1), 0)
    );
    for (int source = 0; source <= level; ++source) {
        const int lower = std::abs(source - factor);
        const int upper
            = std::min(source + factor, 2 * level - source - factor);
        for (int target = lower; target <= upper; ++target) {
            result[static_cast<std::size_t>(target)]
                  [static_cast<std::size_t>(source)] = 1;
        }
    }
    return result;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const std::size_t size = left.size();
    Matrix result(size, Vector(size, 0));
    for (std::size_t row = 0U; row < size; ++row) {
        for (std::size_t middle = 0U; middle < size; ++middle) {
            if (left[row][middle] == 0) {
                continue;
            }
            for (std::size_t column = 0U; column < size; ++column) {
                result[row][column]
                    += left[row][middle] * right[middle][column];
            }
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_half_level = argc == 2
            ? parse_positive(argv[1], "maximum half level")
            : 24;
        if (argc > 2) {
            throw std::invalid_argument(
                "usage: probe_su2_parity_block_tn [maximum_half_level]"
            );
        }

        std::uint64_t checks = 0U;
        for (int level = 1; level <= maximum_half_level; ++level) {
            for (int factor = 1; factor <= level; ++factor) {
                const Matrix fusion = fusion_matrix(level, factor);
                const Matrix square = multiply(fusion, fusion);
                for (int parity = 0; parity < 2; ++parity) {
                    std::vector<int> labels;
                    for (int label = parity; label <= level; label += 2) {
                        labels.push_back(label);
                    }
                    for (std::size_t i = 0U; i < labels.size(); ++i) {
                        for (std::size_t j = i + 1U;
                             j < labels.size(); ++j) {
                            for (std::size_t a = 0U; a < labels.size(); ++a) {
                                for (std::size_t b = a + 1U;
                                     b < labels.size(); ++b) {
                                    const cpp_int minor
                                        = square[static_cast<std::size_t>(
                                                     labels[i]
                                                 )][static_cast<std::size_t>(
                                                     labels[a]
                                                 )]
                                            * square[static_cast<std::size_t>(
                                                         labels[j]
                                                     )][static_cast<std::size_t>(
                                                         labels[b]
                                                     )]
                                            - square[static_cast<std::size_t>(
                                                         labels[i]
                                                     )][static_cast<std::size_t>(
                                                         labels[b]
                                                     )]
                                                * square[static_cast<std::size_t>(
                                                             labels[j]
                                                         )][static_cast<std::size_t>(
                                                             labels[a]
                                                         )];
                                    ++checks;
                                    if (minor < 0) {
                                        std::cout
                                            << "SU2_PARITY_BLOCK_TN"
                                            << " result=COUNTEREXAMPLE"
                                            << " level=" << 2 * level
                                            << " factor=" << 2 * factor
                                            << " parity=" << parity
                                            << " rows=[" << 2 * labels[i]
                                            << ',' << 2 * labels[j] << ']'
                                            << " columns=[" << 2 * labels[a]
                                            << ',' << 2 * labels[b] << ']'
                                            << " minor=" << minor
                                            << " checks=" << checks << '\n';
                                        return EXIT_SUCCESS;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << "SU2_PARITY_BLOCK_TN"
                  << " maximum_level=" << 2 * maximum_half_level
                  << " checks=" << checks
                  << " result=PASS_2X2\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
