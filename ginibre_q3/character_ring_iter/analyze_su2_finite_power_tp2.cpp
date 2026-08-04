#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Matrix = std::vector<std::vector<Integer>>;

int positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<int>(value);
}

bool fuses(int level, int factor, int left, int right) {
    // Labels are compressed by two: this is fusion by the even full label
    // 2*factor on the even full-label sector, so parity is already fixed.
    return std::abs(left - factor) <= right
        && right <= std::min(left + factor, 2 * level - left - factor);
}

Matrix identity(int dimension) {
    Matrix result(
        static_cast<std::size_t>(dimension),
        std::vector<Integer>(static_cast<std::size_t>(dimension))
    );
    for (int index = 0; index < dimension; ++index) {
        result[static_cast<std::size_t>(index)]
              [static_cast<std::size_t>(index)] = 1;
    }
    return result;
}

Matrix fusion_matrix(int level, int factor) {
    Matrix result(
        static_cast<std::size_t>(level + 1),
        std::vector<Integer>(static_cast<std::size_t>(level + 1))
    );
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            if (fuses(level, factor, left, right)) {
                result[static_cast<std::size_t>(left)]
                      [static_cast<std::size_t>(right)] = 1;
            }
        }
    }
    return result;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const int dimension = static_cast<int>(left.size());
    Matrix result(
        static_cast<std::size_t>(dimension),
        std::vector<Integer>(static_cast<std::size_t>(dimension))
    );
    for (int row = 0; row < dimension; ++row) {
        for (int inner = 0; inner < dimension; ++inner) {
            const Integer& coefficient =
                left[static_cast<std::size_t>(row)]
                    [static_cast<std::size_t>(inner)];
            if (coefficient == 0) {
                continue;
            }
            for (int column = 0; column < dimension; ++column) {
                result[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(column)]
                    += coefficient
                       * right[static_cast<std::size_t>(inner)]
                              [static_cast<std::size_t>(column)];
            }
        }
    }
    return result;
}

struct Witness {
    bool found = false;
    int level = 0;
    int factor = 0;
    int power = 0;
    int row = 0;
    int column = 0;
    Integer determinant = 0;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::invalid_argument(
                "usage: analyze_su2_finite_power_tp2 "
                "MAXIMUM_HALF_LEVEL MAXIMUM_POWER"
            );
        }
        const int maximum_level = positive(argv[1], "maximum half level");
        const int maximum_power = positive(argv[2], "maximum power");
        if (maximum_level < 2 || maximum_power < 2) {
            throw std::invalid_argument(
                "maximum half level and maximum power must both be at least two"
            );
        }

        std::vector<std::uint64_t> checks(
            static_cast<std::size_t>(maximum_power + 1)
        );
        std::vector<std::uint64_t> negatives(
            static_cast<std::size_t>(maximum_power + 1)
        );
        std::vector<std::uint64_t> anchored_checks(
            static_cast<std::size_t>(maximum_power + 1)
        );
        std::vector<std::uint64_t> anchored_negatives(
            static_cast<std::size_t>(maximum_power + 1)
        );
        std::vector<Witness> first(
            static_cast<std::size_t>(maximum_power + 1)
        );
        std::vector<Witness> first_anchored(
            static_cast<std::size_t>(maximum_power + 1)
        );

        for (int level = 2; level <= maximum_level; ++level) {
            for (int factor = 1; 2 * factor < level; ++factor) {
                const Matrix fusion = fusion_matrix(level, factor);
                Matrix power = identity(level + 1);
                for (int exponent = 1;
                     exponent <= maximum_power;
                     ++exponent) {
                    power = multiply(power, fusion);
                    if (exponent < 2) {
                        continue;
                    }
                    for (int row = 0; row < level; ++row) {
                        for (int column = 0; column < level; ++column) {
                            const Integer determinant =
                                power[static_cast<std::size_t>(row)]
                                     [static_cast<std::size_t>(column)]
                                * power[static_cast<std::size_t>(row + 1)]
                                       [static_cast<std::size_t>(column + 1)]
                                - power[static_cast<std::size_t>(row)]
                                       [static_cast<std::size_t>(column + 1)]
                                * power[static_cast<std::size_t>(row + 1)]
                                       [static_cast<std::size_t>(column)];
                            ++checks[static_cast<std::size_t>(exponent)];
                            if (determinant < 0) {
                                ++negatives[static_cast<std::size_t>(exponent)];
                                Witness& witness = first[
                                    static_cast<std::size_t>(exponent)
                                ];
                                if (!witness.found) {
                                    witness = Witness{
                                        true,
                                        level,
                                        factor,
                                        exponent,
                                        row,
                                        column,
                                        determinant
                                    };
                                }
                            }
                        }
                    }
                    if ((exponent & 1) == 0) {
                        for (int target = 0; target <= level; ++target) {
                            const Integer determinant =
                                power[0U][0U]
                                * power[static_cast<std::size_t>(factor)]
                                       [static_cast<std::size_t>(target)]
                                - power[0U][static_cast<std::size_t>(target)]
                                * power[static_cast<std::size_t>(factor)][0U];
                            ++anchored_checks[
                                static_cast<std::size_t>(exponent)
                            ];
                            if (determinant < 0) {
                                ++anchored_negatives[
                                    static_cast<std::size_t>(exponent)
                                ];
                                Witness& witness = first_anchored[
                                    static_cast<std::size_t>(exponent)
                                ];
                                if (!witness.found) {
                                    witness = Witness{
                                        true,
                                        level,
                                        factor,
                                        exponent,
                                        0,
                                        target,
                                        determinant
                                    };
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_FINITE_POWER_TP2"
            << " maximum_half_level=" << maximum_level
            << " maximum_power=" << maximum_power
            << '\n';
        for (int exponent = 2; exponent <= maximum_power; ++exponent) {
            const Witness& witness = first[static_cast<std::size_t>(exponent)];
            const Witness& anchored = first_anchored[
                static_cast<std::size_t>(exponent)
            ];
            std::cout
                << "POWER m=" << exponent
                << " adjacent_checks="
                << checks[static_cast<std::size_t>(exponent)]
                << " negative_adjacent_minors="
                << negatives[static_cast<std::size_t>(exponent)]
                << " first_negative=";
            if (witness.found) {
                std::cout
                    << "{K=" << witness.level
                    << " Q=" << witness.factor
                    << " rows=[" << witness.row << ',' << witness.row + 1
                    << "] columns=[" << witness.column << ','
                    << witness.column + 1 << "] value="
                    << witness.determinant << '}';
            } else {
                std::cout << "{}";
            }
            if ((exponent & 1) == 0) {
                std::cout
                    << " anchored_checks="
                    << anchored_checks[static_cast<std::size_t>(exponent)]
                    << " negative_anchored_minors="
                    << anchored_negatives[
                        static_cast<std::size_t>(exponent)
                    ]
                    << " first_negative_anchored=";
                if (anchored.found) {
                    std::cout
                        << "{K=" << anchored.level
                        << " Q=" << anchored.factor
                        << " columns=[0," << anchored.column
                        << "] value=" << anchored.determinant << '}';
                } else {
                    std::cout << "{}";
                }
            }
            std::cout << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_FINITE_POWER_TP2 FAILURE: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
