#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using IntegerVector = std::vector<Integer>;
using IntegerMatrix = std::vector<IntegerVector>;
using RationalVector = std::vector<Rational>;
using RationalMatrix = std::vector<RationalVector>;

int parse_positive(const char* text, const char* name) {
    const std::string value(text);
    std::size_t used = 0U;
    const long long parsed = std::stoll(value, &used);
    if (used != value.size() || parsed <= 0LL) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

struct FoldedLabel {
    int label = -1;
    int sign = 0;
};

FoldedLabel fold_character(int degree, int level) {
    if (degree < 0) {
        return {};
    }
    const int period = 2 * (level + 2);
    const int residue = degree % period;
    if (residue <= level) {
        return {residue, 1};
    }
    if (residue == level + 1) {
        return {};
    }
    const int reflected = period - residue - 2;
    return reflected < 0 ? FoldedLabel{} : FoldedLabel{reflected, -1};
}

IntegerMatrix fusion_square_matrix(int level, int factor) {
    const int dimension = level + 1;
    IntegerMatrix first(
        static_cast<std::size_t>(dimension),
        IntegerVector(static_cast<std::size_t>(dimension), 0)
    );
    for (int source = 0; source <= level; ++source) {
        const int lower = std::abs(source - factor);
        const int upper = std::min(
            source + factor,
            2 * level - source - factor
        );
        for (int target = lower; target <= upper; target += 2) {
            first[static_cast<std::size_t>(target)]
                 [static_cast<std::size_t>(source)] = 1;
        }
    }
    IntegerMatrix square(
        static_cast<std::size_t>(dimension),
        IntegerVector(static_cast<std::size_t>(dimension), 0)
    );
    for (int row = 0; row < dimension; ++row) {
        for (int middle = 0; middle < dimension; ++middle) {
            if (first[static_cast<std::size_t>(row)]
                     [static_cast<std::size_t>(middle)] == 0) {
                continue;
            }
            for (int column = 0; column < dimension; ++column) {
                square[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(column)] +=
                    first[static_cast<std::size_t>(middle)]
                         [static_cast<std::size_t>(column)];
            }
        }
    }
    return square;
}

RationalMatrix inverse(const IntegerMatrix& input) {
    const int dimension = static_cast<int>(input.size());
    RationalMatrix augmented(
        static_cast<std::size_t>(dimension),
        RationalVector(static_cast<std::size_t>(2 * dimension), 0)
    );
    for (int row = 0; row < dimension; ++row) {
        for (int column = 0; column < dimension; ++column) {
            augmented[static_cast<std::size_t>(row)]
                     [static_cast<std::size_t>(column)] =
                input[static_cast<std::size_t>(row)]
                     [static_cast<std::size_t>(column)];
        }
        augmented[static_cast<std::size_t>(row)]
                 [static_cast<std::size_t>(dimension + row)] = 1;
    }
    for (int column = 0; column < dimension; ++column) {
        int pivot = column;
        while (pivot < dimension
               && augmented[static_cast<std::size_t>(pivot)]
                           [static_cast<std::size_t>(column)] == 0) {
            ++pivot;
        }
        if (pivot == dimension) {
            throw std::runtime_error("Pieri vectors are not a basis");
        }
        std::swap(
            augmented[static_cast<std::size_t>(column)],
            augmented[static_cast<std::size_t>(pivot)]
        );
        const Rational pivot_value = augmented[static_cast<std::size_t>(column)]
            [static_cast<std::size_t>(column)];
        for (int entry = 0; entry < 2 * dimension; ++entry) {
            augmented[static_cast<std::size_t>(column)]
                     [static_cast<std::size_t>(entry)] /= pivot_value;
        }
        for (int row = 0; row < dimension; ++row) {
            if (row == column) {
                continue;
            }
            const Rational scale = augmented[static_cast<std::size_t>(row)]
                [static_cast<std::size_t>(column)];
            if (scale == 0) {
                continue;
            }
            for (int entry = 0; entry < 2 * dimension; ++entry) {
                augmented[static_cast<std::size_t>(row)]
                         [static_cast<std::size_t>(entry)] -= scale
                    * augmented[static_cast<std::size_t>(column)]
                               [static_cast<std::size_t>(entry)];
            }
        }
    }
    RationalMatrix result(
        static_cast<std::size_t>(dimension),
        RationalVector(static_cast<std::size_t>(dimension), 0)
    );
    for (int row = 0; row < dimension; ++row) {
        for (int column = 0; column < dimension; ++column) {
            result[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(column)] =
                augmented[static_cast<std::size_t>(row)]
                         [static_cast<std::size_t>(dimension + column)];
        }
    }
    return result;
}

IntegerVector multiply(const IntegerMatrix& matrix, const IntegerVector& input) {
    IntegerVector result(matrix.size(), 0);
    for (std::size_t row = 0U; row < matrix.size(); ++row) {
        for (std::size_t column = 0U; column < input.size(); ++column) {
            result[row] += matrix[row][column] * input[column];
        }
    }
    return result;
}

RationalVector multiply(const RationalMatrix& matrix, const IntegerVector& input) {
    RationalVector result(matrix.size(), 0);
    for (std::size_t row = 0U; row < matrix.size(); ++row) {
        for (std::size_t column = 0U; column < input.size(); ++column) {
            result[row] += matrix[row][column] * input[column];
        }
    }
    return result;
}

std::string render_pair(std::pair<int, int> pair) {
    return "(" + std::to_string(pair.first) + ','
        + std::to_string(pair.second) + ')';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: probe_su2_finite_c2_pieri_cone MAXIMUM_LEVEL"
            );
        }
        const int maximum_level = parse_positive(argv[1], "maximum level");
        std::uint64_t generators = 0U;
        std::uint64_t coefficient_checks = 0U;
        for (int level = 1; level <= maximum_level; ++level) {
            const int dimension = level * (level + 1) / 2;
            std::vector<std::pair<int, int>> wedge_pairs;
            std::vector<std::vector<int>> wedge_index(
                static_cast<std::size_t>(level + 1),
                std::vector<int>(static_cast<std::size_t>(level + 1), -1)
            );
            for (int high = 1; high <= level; ++high) {
                for (int low = 0; low < high; ++low) {
                    wedge_index[static_cast<std::size_t>(high)]
                               [static_cast<std::size_t>(low)] =
                        static_cast<int>(wedge_pairs.size());
                    wedge_pairs.emplace_back(high, low);
                }
            }
            std::vector<std::pair<int, int>> pieri_pairs;
            IntegerMatrix pieri(
                static_cast<std::size_t>(dimension),
                IntegerVector(static_cast<std::size_t>(dimension), 0)
            );
            for (int left = 1; left <= level; ++left) {
                for (int right = left; right <= level; ++right) {
                    const int column = static_cast<int>(pieri_pairs.size());
                    pieri_pairs.emplace_back(left, right);
                    for (int index = 0; index < left; ++index) {
                        for (int contraction = 0;
                             contraction <= left - 1 - index;
                             ++contraction) {
                            const int first = left + right - 1 - index
                                - 2 * contraction;
                            const int second = index;
                            const FoldedLabel folded_first = fold_character(
                                first, level
                            );
                            const FoldedLabel folded_second = fold_character(
                                second, level
                            );
                            if (folded_first.sign == 0
                                || folded_second.sign == 0
                                || folded_first.label == folded_second.label) {
                                continue;
                            }
                            int high = folded_first.label;
                            int low = folded_second.label;
                            int sign = folded_first.sign * folded_second.sign;
                            if (high < low) {
                                std::swap(high, low);
                                sign = -sign;
                            }
                            const int row = wedge_index[
                                static_cast<std::size_t>(high)
                            ][static_cast<std::size_t>(low)];
                            if (row < 0) {
                                throw std::runtime_error("invalid wedge row");
                            }
                            pieri[static_cast<std::size_t>(row)]
                                 [static_cast<std::size_t>(column)] += sign;
                        }
                    }
                }
            }
            if (static_cast<int>(pieri_pairs.size()) != dimension) {
                throw std::runtime_error("Pieri dimension mismatch");
            }
            const RationalMatrix inverse_pieri = inverse(pieri);
            for (int factor = 1; factor <= level; ++factor) {
                const IntegerMatrix square = fusion_square_matrix(level, factor);
                IntegerMatrix compound(
                    static_cast<std::size_t>(dimension),
                    IntegerVector(static_cast<std::size_t>(dimension), 0)
                );
                for (int source = 0; source < dimension; ++source) {
                    const auto [source_high, source_low] = wedge_pairs[
                        static_cast<std::size_t>(source)
                    ];
                    for (int target = 0; target < dimension; ++target) {
                        const auto [target_high, target_low] = wedge_pairs[
                            static_cast<std::size_t>(target)
                        ];
                        compound[static_cast<std::size_t>(target)]
                                [static_cast<std::size_t>(source)] =
                            square[static_cast<std::size_t>(target_high)]
                                  [static_cast<std::size_t>(source_high)]
                              * square[static_cast<std::size_t>(target_low)]
                                      [static_cast<std::size_t>(source_low)]
                            - square[static_cast<std::size_t>(target_high)]
                                  [static_cast<std::size_t>(source_low)]
                              * square[static_cast<std::size_t>(target_low)]
                                      [static_cast<std::size_t>(source_high)];
                    }
                }
                ++generators;
                for (int source = 0; source < dimension; ++source) {
                    IntegerVector column(static_cast<std::size_t>(dimension), 0);
                    for (int row = 0; row < dimension; ++row) {
                        column[static_cast<std::size_t>(row)] = pieri[
                            static_cast<std::size_t>(row)
                        ][static_cast<std::size_t>(source)];
                    }
                    const IntegerVector image = multiply(compound, column);
                    const RationalVector coefficients = multiply(
                        inverse_pieri, image
                    );
                    for (int target = 0; target < dimension; ++target) {
                        ++coefficient_checks;
                        const Rational& coefficient = coefficients[
                            static_cast<std::size_t>(target)
                        ];
                        if (coefficient < 0) {
                            std::cout
                                << "FINITE_C2_PIERI_CONE_FAIL"
                                << " level=" << level
                                << " factor=" << factor
                                << " source=" << render_pair(pieri_pairs[
                                    static_cast<std::size_t>(source)
                                ])
                                << " target=" << render_pair(pieri_pairs[
                                    static_cast<std::size_t>(target)
                                ])
                                << " coefficient=" << coefficient.numerator();
                            if (coefficient.denominator() != 1) {
                                std::cout << '/' << coefficient.denominator();
                            }
                            std::cout << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }
        std::cout
            << "SU2_FINITE_C2_PIERI_CONE"
            << " maximum_level=" << maximum_level
            << " generators=" << generators
            << " coefficient_checks=" << coefficient_checks
            << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
