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

int parse_nonnegative(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value < 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string(name) + " must be nonnegative");
    }
    return static_cast<int>(value);
}

Matrix identity(int size) {
    Matrix result(
        static_cast<std::size_t>(size),
        std::vector<Integer>(static_cast<std::size_t>(size))
    );
    for (int index = 0; index < size; ++index) {
        result[static_cast<std::size_t>(index)]
              [static_cast<std::size_t>(index)] = 1;
    }
    return result;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const int size = static_cast<int>(left.size());
    Matrix result(
        static_cast<std::size_t>(size),
        std::vector<Integer>(static_cast<std::size_t>(size))
    );
    for (int row = 0; row < size; ++row) {
        for (int middle = 0; middle < size; ++middle) {
            const Integer& factor = left[static_cast<std::size_t>(row)]
                                            [static_cast<std::size_t>(middle)];
            if (factor == 0) {
                continue;
            }
            for (int column = 0; column < size; ++column) {
                result[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(column)] +=
                    factor * right[static_cast<std::size_t>(middle)]
                                  [static_cast<std::size_t>(column)];
            }
        }
    }
    return result;
}

Matrix fusion_matrix(int level, int label) {
    const int size = level + 1;
    Matrix result(
        static_cast<std::size_t>(size),
        std::vector<Integer>(static_cast<std::size_t>(size))
    );
    for (int source = 0; source <= level; ++source) {
        const int upper = std::min(
            source + label, 2 * level - source - label
        );
        for (int target = std::abs(source - label);
             target <= upper; target += 2) {
            result[static_cast<std::size_t>(source)]
                  [static_cast<std::size_t>(target)] = 1;
        }
    }
    return result;
}

Matrix fundamental_square(int level) {
    const Matrix fundamental = fusion_matrix(level, 1);
    return multiply(fundamental, fundamental);
}

Matrix label_two_square(int level) {
    const Matrix label_two = fusion_matrix(level, 2);
    return multiply(label_two, label_two);
}

void verify_case(
    const Matrix& even_power,
    int level,
    int factors
) {
    const Integer d0 = even_power[0U][0U];
    for (int row = 0; row <= level; ++row) {
        for (int column = 0; column <= level; ++column) {
            const Integer current = d0
                * even_power[static_cast<std::size_t>(row)]
                            [static_cast<std::size_t>(column)]
                - even_power[static_cast<std::size_t>(row)][0U]
                    * even_power[static_cast<std::size_t>(column)][0U];
            if (current < 0) {
                throw std::runtime_error(
                    "negative anchored current at level="
                    + std::to_string(level)
                    + " factors=" + std::to_string(factors)
                    + " row=" + std::to_string(row)
                    + " column=" + std::to_string(column)
                );
            }
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4 && std::string(argv[1]) == "--one-two") {
            const int maximum_level = parse_nonnegative(
                argv[2], "maximum level"
            );
            const int maximum_factors = parse_nonnegative(
                argv[3], "maximum factors"
            );
            std::uint64_t cases = 0U;
            std::uint64_t currents = 0U;
            for (int level = 2; level <= maximum_level; ++level) {
                Matrix fundamental_power = identity(level + 1);
                const Matrix fundamental_step = fundamental_square(level);
                const Matrix label_two_step = label_two_square(level);
                for (int fundamental_factors = 0;
                     fundamental_factors <= maximum_factors;
                     ++fundamental_factors) {
                    Matrix label_two_power = identity(level + 1);
                    for (int label_two_factors = 0;
                         label_two_factors <= maximum_factors;
                         ++label_two_factors) {
                        const Matrix word_square = multiply(
                            fundamental_power, label_two_power
                        );
                        verify_case(
                            word_square, level,
                            fundamental_factors + label_two_factors
                        );
                        ++cases;
                        const std::uint64_t width = static_cast<std::uint64_t>(
                            level + 1
                        );
                        currents += width * width;
                        label_two_power = multiply(
                            label_two_step, label_two_power
                        );
                    }
                    fundamental_power = multiply(
                        fundamental_step, fundamental_power
                    );
                }
            }
            std::cout << "SU2_ONE_TWO_ANCHORED"
                      << " levels=" << maximum_level
                      << " factors=" << maximum_factors
                      << " cases=" << cases
                      << " currents=" << currents
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc == 3 && std::string(argv[1]) == "--one-factor") {
            const int maximum_level = parse_nonnegative(
                argv[2], "maximum level"
            );
            std::uint64_t cases = 0U;
            std::uint64_t currents = 0U;
            for (int level = 0; level <= maximum_level; ++level) {
                for (int label = 0; label <= level; ++label) {
                    const Matrix fusion = fusion_matrix(level, label);
                    const Matrix square = multiply(fusion, fusion);
                    verify_case(square, level, label);
                    ++cases;
                    const std::uint64_t width = static_cast<std::uint64_t>(
                        level + 1
                    );
                    currents += width * width;
                }
            }
            std::cout << "SU2_ONE_FACTOR_ANCHORED"
                      << " levels=" << maximum_level
                      << " cases=" << cases
                      << " currents=" << currents
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: verify_su2_fundamental_anchored_tn "
                "MAXIMUM_LEVEL MAXIMUM_FACTORS"
                " | --one-factor MAXIMUM_LEVEL"
                " | --one-two MAXIMUM_LEVEL MAXIMUM_FACTORS"
            );
        }
        const int maximum_level = parse_nonnegative(
            argv[1], "maximum level"
        );
        const int maximum_factors = parse_nonnegative(
            argv[2], "maximum factors"
        );
        std::uint64_t cases = 0U;
        std::uint64_t currents = 0U;
        for (int level = 1; level <= maximum_level; ++level) {
            Matrix even_power = identity(level + 1);
            const Matrix step = fundamental_square(level);
            for (int factors = 0; factors <= maximum_factors; ++factors) {
                verify_case(even_power, level, factors);
                ++cases;
                const std::uint64_t width = static_cast<std::uint64_t>(
                    level + 1
                );
                currents += width * width;
                even_power = multiply(step, even_power);
            }
        }
        std::cout << "SU2_FUNDAMENTAL_ANCHORED_TN"
                  << " levels=" << maximum_level
                  << " factors=" << maximum_factors
                  << " cases=" << cases
                  << " currents=" << currents
                  << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_FUNDAMENTAL_ANCHORED_TN FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
