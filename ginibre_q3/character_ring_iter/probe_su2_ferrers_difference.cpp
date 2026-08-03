#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Matrix = std::vector<std::vector<int>>;
using Integer = boost::multiprecision::cpp_int;
using BigMatrix = std::vector<std::vector<Integer>>;

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(value);
}

bool fuses_half(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= source + label
        && source + target + label <= 2 * level;
}

Matrix zero_matrix(int size) {
    return Matrix(
        static_cast<std::size_t>(size),
        std::vector<int>(static_cast<std::size_t>(size))
    );
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const int size = static_cast<int>(left.size());
    Matrix result = zero_matrix(size);
    for (int row = 0; row < size; ++row) {
        for (int middle = 0; middle < size; ++middle) {
            const int factor = left[static_cast<std::size_t>(row)]
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

BigMatrix promote(const Matrix& matrix) {
    const int size = static_cast<int>(matrix.size());
    BigMatrix result(
        static_cast<std::size_t>(size),
        std::vector<Integer>(static_cast<std::size_t>(size))
    );
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            result[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(column)] =
                matrix[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(column)];
        }
    }
    return result;
}

BigMatrix multiply(const BigMatrix& left, const BigMatrix& right) {
    const int size = static_cast<int>(left.size());
    BigMatrix result(
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

BigMatrix identity(int size) {
    BigMatrix result(
        static_cast<std::size_t>(size),
        std::vector<Integer>(static_cast<std::size_t>(size))
    );
    for (int index = 0; index < size; ++index) {
        result[static_cast<std::size_t>(index)]
              [static_cast<std::size_t>(index)] = 1;
    }
    return result;
}

Matrix first_difference(int size) {
    Matrix result = zero_matrix(size);
    for (int row = 0; row < size; ++row) {
        result[static_cast<std::size_t>(row)]
              [static_cast<std::size_t>(row)] = 1;
        if (row != 0) {
            result[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(row - 1)] = -1;
        }
    }
    return result;
}

Matrix cumulative(int size) {
    Matrix result = zero_matrix(size);
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column <= row; ++column) {
            result[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(column)] = 1;
        }
    }
    return result;
}

Matrix reversal(int size) {
    Matrix result = zero_matrix(size);
    for (int row = 0; row < size; ++row) {
        result[static_cast<std::size_t>(row)]
              [static_cast<std::size_t>(size - 1 - row)] = 1;
    }
    return result;
}

int negative_entries(const Matrix& matrix) {
    int negatives = 0;
    for (const auto& row : matrix) {
        for (const int value : row) {
            if (value < 0) {
                ++negatives;
            }
        }
    }
    return negatives;
}

void print_matrix(const char* name, const Matrix& matrix) {
    std::cout << name << '=';
    for (const auto& row : matrix) {
        std::cout << '[';
        for (std::size_t column = 0U; column < row.size(); ++column) {
            std::cout << row[column];
            if (column + 1U != row.size()) {
                std::cout << ',';
            }
        }
        std::cout << ']';
    }
    std::cout << '\n';
}

struct Result {
    int matrices = 0;
    int crossing_failures = 0;
    int plus_cone_failures = 0;
    int minus_cone_failures = 0;
    int anchored_entries = 0;
    int anchored_failures = 0;
    bool has_witness = false;
    int witness_level = 0;
    int witness_label = 0;
    int witness_row = 0;
    int witness_column = 0;
    int witness_value = 0;
    bool has_anchor_witness = false;
    int anchor_level = 0;
    int anchor_label = 0;
    int anchor_power = 0;
    int anchor_row = 0;
    int anchor_column = 0;
    Integer anchor_value = 0;
};

void inspect(
    int level,
    int label,
    int maximum_anchor_power,
    bool print,
    Result& result
) {
    if (2 * label >= level) {
        return;
    }
    const int size = (level + 1) / 2;
    Matrix plus = zero_matrix(size);
    Matrix minus = zero_matrix(size);
    for (int source = 0; source < size; ++source) {
        for (int target = 0; target < size; ++target) {
            const int same = fuses_half(level, label, source, target) ? 1 : 0;
            const int crossed = fuses_half(
                level,
                label,
                source,
                level - target
            ) ? 1 : 0;
            if (crossed > same) {
                throw std::runtime_error("crossed edge lacks same-side edge");
            }
            plus[static_cast<std::size_t>(source)]
                [static_cast<std::size_t>(target)] = same + crossed;
            minus[static_cast<std::size_t>(source)]
                 [static_cast<std::size_t>(target)] = same - crossed;
        }
    }
    const Matrix difference = first_difference(size);
    const Matrix sum = cumulative(size);
    const Matrix reflect = reversal(size);
    Matrix raw_crossing = zero_matrix(size);
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            raw_crossing[static_cast<std::size_t>(row)]
                        [static_cast<std::size_t>(column)] =
                plus[static_cast<std::size_t>(row)]
                    [static_cast<std::size_t>(column)]
                - minus[static_cast<std::size_t>(row)]
                    [static_cast<std::size_t>(column)];
        }
    }
    const Matrix shifted_crossing = multiply(
        difference,
        multiply(raw_crossing, reflect)
    );
    const Matrix plus_cone = multiply(
        difference,
        multiply(plus, sum)
    );
    const Matrix minus_cone = multiply(
        difference,
        multiply(minus, sum)
    );
    ++result.matrices;
    result.crossing_failures += negative_entries(shifted_crossing);
    result.plus_cone_failures += negative_entries(plus_cone);
    result.minus_cone_failures += negative_entries(minus_cone);
    int crossing_shift = size;
    bool shifted_ferrers = true;
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            const int value = shifted_crossing[static_cast<std::size_t>(row)]
                                               [static_cast<std::size_t>(column)];
            if (value == 0) {
                continue;
            }
            crossing_shift = std::min(crossing_shift, row - column);
            if (value != 2 || row < column) {
                shifted_ferrers = false;
            }
        }
    }
    if (crossing_shift != size) {
        for (int row = 0; row < size; ++row) {
            for (int column = 0; column < size; ++column) {
                const int expected = row == column + crossing_shift ? 2 : 0;
                if (shifted_crossing[static_cast<std::size_t>(row)]
                                    [static_cast<std::size_t>(column)] != expected) {
                    shifted_ferrers = false;
                }
            }
        }
    }
    if (!shifted_ferrers) {
        throw std::runtime_error("crossing is not a shifted Ferrers kernel");
    }
    if (maximum_anchor_power >= 0 && crossing_shift != size) {
        const BigMatrix difference_big = promote(difference);
        const BigMatrix cumulative_big = promote(sum);
        const BigMatrix plus_big = promote(plus);
        BigMatrix power = identity(size);
        for (int exponent = 0; exponent <= maximum_anchor_power; ++exponent) {
            const BigMatrix transformed = multiply(
                difference_big,
                multiply(power, cumulative_big)
            );
            for (int row = 0; row < size; ++row) {
                for (int column = crossing_shift;
                     column < size;
                     ++column) {
                    ++result.anchored_entries;
                    const Integer& value = transformed[
                        static_cast<std::size_t>(row)
                    ][static_cast<std::size_t>(column)];
                    if (value < 0) {
                        ++result.anchored_failures;
                        if (!result.has_anchor_witness) {
                            result.has_anchor_witness = true;
                            result.anchor_level = level;
                            result.anchor_label = label;
                            result.anchor_power = exponent;
                            result.anchor_row = row;
                            result.anchor_column = column;
                            result.anchor_value = value;
                        }
                    }
                }
            }
            power = multiply(power, plus_big);
        }
    }
    const Matrix* candidates[] = {&shifted_crossing, &plus_cone, &minus_cone};
    for (const Matrix* candidate : candidates) {
        for (int row = 0; row < size; ++row) {
            for (int column = 0; column < size; ++column) {
                const int value = (*candidate)[static_cast<std::size_t>(row)]
                                              [static_cast<std::size_t>(column)];
                if (value < 0 && !result.has_witness) {
                    result.has_witness = true;
                    result.witness_level = level;
                    result.witness_label = label;
                    result.witness_row = row;
                    result.witness_column = column;
                    result.witness_value = value;
                }
            }
        }
    }
    if (print) {
        print_matrix("A_PLUS", plus);
        print_matrix("A_MINUS", minus);
        print_matrix("D_DELTA_J", shifted_crossing);
        print_matrix("D_A_PLUS_L", plus_cone);
        print_matrix("D_A_MINUS_L", minus_cone);
        if (crossing_shift != size) {
            std::cout << "CROSSING_SHIFT=" << crossing_shift << '\n';
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4 && std::string(argv[1]) == "--case") {
            const int level = parse_positive(argv[2], "half level");
            const int label = parse_positive(argv[3], "half label");
            Result result;
            inspect(level, label, -1, true, result);
            std::cout << "SU2_FERRERS_DIFFERENCE_CASE"
                      << " half_level=" << level
                      << " half_label=" << label
                      << " crossing_negative=" << result.crossing_failures
                      << " plus_cone_negative=" << result.plus_cone_failures
                      << " minus_cone_negative=" << result.minus_cone_failures
                      << " result="
                      << (
                          result.crossing_failures == 0
                          && result.plus_cone_failures == 0
                          && result.minus_cone_failures == 0
                              ? "PASS"
                              : "COUNTEREXAMPLE"
                      )
                      << '\n';
            return result.crossing_failures == 0
                && result.plus_cone_failures == 0
                && result.minus_cone_failures == 0
                ? EXIT_SUCCESS
                : EXIT_FAILURE;
        }
        if (argc == 4 && std::string(argv[1]) == "--anchor-scan") {
            const int maximum_level = parse_positive(
                argv[2],
                "maximum half level"
            );
            const int maximum_power = parse_positive(
                argv[3],
                "maximum anchored power"
            );
            Result result;
            for (int level = 3; level <= maximum_level; ++level) {
                for (int label = 1; 2 * label < level; ++label) {
                    inspect(level, label, maximum_power, false, result);
                }
            }
            std::cout << "SU2_FERRERS_ANCHORED_SCAN"
                      << " maximum_half_level=" << maximum_level
                      << " maximum_power=" << maximum_power
                      << " entries=" << result.anchored_entries
                      << " negative=" << result.anchored_failures;
            if (result.has_anchor_witness) {
                std::cout << " first_witness=("
                          << result.anchor_level << ','
                          << result.anchor_label << ','
                          << result.anchor_power << ','
                          << result.anchor_row << ','
                          << result.anchor_column << ','
                          << result.anchor_value << ')';
            }
            std::cout << " result="
                      << (
                          result.anchored_failures == 0
                              ? "PASS"
                              : "COUNTEREXAMPLE"
                      )
                      << '\n';
            return result.anchored_failures == 0
                ? EXIT_SUCCESS
                : EXIT_FAILURE;
        }
        if (argc != 3 || std::string(argv[1]) != "--scan") {
            throw std::runtime_error(
                "usage: probe_su2_ferrers_difference "
                "--case HALF_LEVEL HALF_LABEL | --scan MAXIMUM_HALF_LEVEL "
                "| --anchor-scan MAXIMUM_HALF_LEVEL MAXIMUM_POWER"
            );
        }
        const int maximum_level = parse_positive(argv[2], "maximum half level");
        Result result;
        for (int level = 3; level <= maximum_level; ++level) {
            for (int label = 1; 2 * label < level; ++label) {
                inspect(level, label, -1, false, result);
            }
        }
        std::cout << "SU2_FERRERS_DIFFERENCE_SCAN"
                  << " maximum_half_level=" << maximum_level
                  << " matrices=" << result.matrices
                  << " crossing_negative=" << result.crossing_failures
                  << " plus_cone_negative=" << result.plus_cone_failures
                  << " minus_cone_negative=" << result.minus_cone_failures;
        if (result.has_witness) {
            std::cout << " first_witness=("
                      << result.witness_level << ','
                      << result.witness_label << ','
                      << result.witness_row << ','
                      << result.witness_column << ','
                      << result.witness_value << ')';
        }
        std::cout << " result="
                  << (
                      result.crossing_failures == 0
                      && result.plus_cone_failures == 0
                      && result.minus_cone_failures == 0
                          ? "PASS"
                          : "COUNTEREXAMPLE"
                  )
                  << '\n';
        return result.crossing_failures == 0
            && result.plus_cone_failures == 0
            && result.minus_cone_failures == 0
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
