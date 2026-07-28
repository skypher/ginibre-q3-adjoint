#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using Matrix = std::vector<std::vector<int>>;

int parse_positive(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value <= 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error("maximum half-level must be positive");
    }
    return static_cast<int>(value);
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const int size = static_cast<int>(left.size());
    Matrix product(
        static_cast<std::size_t>(size),
        std::vector<int>(static_cast<std::size_t>(size))
    );
    for (int row = 0; row < size; ++row) {
        for (int middle = 0; middle < size; ++middle) {
            if (left[static_cast<std::size_t>(row)]
                    [static_cast<std::size_t>(middle)] == 0) {
                continue;
            }
            for (int column = 0; column < size; ++column) {
                product[static_cast<std::size_t>(row)]
                       [static_cast<std::size_t>(column)] +=
                    left[static_cast<std::size_t>(row)]
                        [static_cast<std::size_t>(middle)]
                    * right[static_cast<std::size_t>(middle)]
                           [static_cast<std::size_t>(column)];
            }
        }
    }
    return product;
}

bool fuses_half(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= source + label
        && source + target + label <= 2 * level;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2 && argc != 3) {
            throw std::runtime_error(
                "usage: MAXIMUM_HALF_LEVEL [OPERATOR_POWER]"
            );
        }
        const int maximum_half_level = parse_positive(argv[1]);
        const int operator_power =
            argc == 3 ? parse_positive(argv[2]) : 1;
        std::uint64_t parameters = 0U;
        std::uint64_t suffix_rows = 0U;
        std::uint64_t differences = 0U;
        int minimum_difference = 0;
        int maximum_difference = 0;

        for (int level = 3; level <= maximum_half_level; ++level) {
            for (int label = 1; 2 * label < level; ++label) {
                ++parameters;
                const int paired = (level + 1) / 2;
                Matrix odd(
                    static_cast<std::size_t>(paired),
                    std::vector<int>(static_cast<std::size_t>(paired))
                );
                for (int source = 0; source < paired; ++source) {
                    for (int target = 0; target < paired; ++target) {
                        const int same = fuses_half(
                            level,
                            label,
                            source,
                            target
                        ) ? 1 : 0;
                        const int crossed = fuses_half(
                            level,
                            label,
                            source,
                            level - target
                        ) ? 1 : 0;
                        if (crossed > same) {
                            throw std::runtime_error(
                                "crossed edge lacks same-side edge"
                            );
                        }
                        odd[static_cast<std::size_t>(source)]
                           [static_cast<std::size_t>(target)] =
                            same - crossed;
                    }
                }
                Matrix powered(
                    static_cast<std::size_t>(paired),
                    std::vector<int>(static_cast<std::size_t>(paired))
                );
                for (int vertex = 0; vertex < paired; ++vertex) {
                    powered[static_cast<std::size_t>(vertex)]
                           [static_cast<std::size_t>(vertex)] = 1;
                }
                for (int exponent = 0;
                     exponent < operator_power;
                     ++exponent) {
                    powered = multiply(powered, odd);
                }

                std::vector<int> column_counts(
                    static_cast<std::size_t>(paired)
                );
                for (int rho = paired - 1; rho >= 0; --rho) {
                    ++suffix_rows;
                    for (int target = 0; target < paired; ++target) {
                        column_counts[static_cast<std::size_t>(target)]
                            += powered[static_cast<std::size_t>(rho)]
                                      [static_cast<std::size_t>(
                                          (
                                              level % 2 == 1
                                              && (operator_power & 1) != 0
                                          )
                                              ? paired - 1 - target
                                              : target
                                      )];
                    }
                    for (int target = 1; target < paired; ++target) {
                        const int difference =
                            column_counts[static_cast<std::size_t>(target)]
                            - column_counts[
                                static_cast<std::size_t>(target - 1)
                            ];
                        ++differences;
                        minimum_difference = std::min(
                            minimum_difference,
                            difference
                        );
                        maximum_difference = std::max(
                            maximum_difference,
                            difference
                        );
                        if (difference < 0) {
                            std::cout
                                << "SU2_SUFFIX_CONE_INVARIANCE"
                                << " counterexample"
                                << " level=" << level
                                << " label=" << label
                                << " operator_power=" << operator_power
                                << " rho=" << rho
                                << " target=" << target
                                << " previous="
                                << column_counts[
                                    static_cast<std::size_t>(target - 1)
                                ]
                                << " current="
                                << column_counts[
                                    static_cast<std::size_t>(target)
                                ]
                                << " difference=" << difference
                                << " result=FAIL_INVARIANCE\n";
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_SUFFIX_CONE_INVARIANCE"
            << " maximum_half_level=" << maximum_half_level
            << " operator_power=" << operator_power
            << " parameters=" << parameters
            << " suffix_rows=" << suffix_rows
            << " differences=" << differences
            << " minimum_difference=" << minimum_difference
            << " maximum_difference=" << maximum_difference
            << " result=PASS_EXACT_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_SUFFIX_CONE_INVARIANCE FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
