#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Matrix = std::vector<std::vector<int>>;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(const bool condition, const std::string& message) {
    if (!condition) {
        fail(message);
    }
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    require(left.size() == right.size(), "matrix-size mismatch");
    const std::size_t size = left.size();
    Matrix result(size, std::vector<int>(size));
    for (std::size_t row = 0U; row < size; ++row) {
        require(left[row].size() == size, "left matrix is not square");
        for (std::size_t middle = 0U; middle < size; ++middle) {
            require(
                right[middle].size() == size,
                "right matrix is not square"
            );
            if (left[row][middle] == 0) {
                continue;
            }
            for (std::size_t column = 0U; column < size; ++column) {
                result[row][column] += left[row][middle] * right[middle][column];
            }
        }
    }
    return result;
}

Matrix transpose(const Matrix& matrix) {
    const std::size_t size = matrix.size();
    Matrix result(size, std::vector<int>(size));
    for (std::size_t row = 0U; row < size; ++row) {
        require(matrix[row].size() == size, "matrix is not square");
        for (std::size_t column = 0U; column < size; ++column) {
            result[column][row] = matrix[row][column];
        }
    }
    return result;
}

bool fuses_half(
    const int half_level,
    const int half_label,
    const int source,
    const int target
) {
    return std::abs(source - half_label) <= target
        && target <= source + half_label
        && source + target + half_label <= 2 * half_level;
}

Matrix core_fusion_matrix(
    const int half_level,
    const int half_label,
    const int block_count
) {
    const int size = block_count * half_label;
    Matrix result(
        static_cast<std::size_t>(size),
        std::vector<int>(static_cast<std::size_t>(size))
    );
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            result[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(column)] = fuses_half(
                half_level,
                half_label,
                half_label + row,
                half_label + column
            ) ? 1 : 0;
        }
    }
    return result;
}

Matrix expected_core_matrix(const int half_label, const int block_count) {
    const int size = block_count * half_label;
    Matrix result(
        static_cast<std::size_t>(size),
        std::vector<int>(static_cast<std::size_t>(size))
    );
    for (int row = 0; row < size; ++row) {
        const int row_block = row / half_label;
        const int row_local = row % half_label;
        for (int column = 0; column < size; ++column) {
            const int column_block = column / half_label;
            const int column_local = column % half_label;
            const bool same_block = row_block == column_block;
            const bool upper_block = column_block == row_block + 1
                && column_local <= row_local;
            const bool lower_block = row_block == column_block + 1
                && row_local <= column_local;
            result[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(column)] =
                same_block || upper_block || lower_block ? 1 : 0;
        }
    }
    return result;
}

Matrix block_difference(const int half_label, const int block_count) {
    const int size = block_count * half_label;
    Matrix result(
        static_cast<std::size_t>(size),
        std::vector<int>(static_cast<std::size_t>(size))
    );
    for (int block = 0; block < block_count; ++block) {
        for (int local = 0; local < half_label; ++local) {
            const int row = block * half_label + local;
            result[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(row)] = 1;
            if (local > 0) {
                result[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(row - 1)] = -1;
            }
        }
    }
    return result;
}

Matrix expected_difference_form(const int half_label, const int block_count) {
    const int size = block_count * half_label;
    Matrix result(
        static_cast<std::size_t>(size),
        std::vector<int>(static_cast<std::size_t>(size))
    );
    for (int block = 0; block < block_count; ++block) {
        const int anchor = block * half_label;
        result[static_cast<std::size_t>(anchor)]
              [static_cast<std::size_t>(anchor)] = 1;
    }
    for (int block = 0; block + 1 < block_count; ++block) {
        for (int local = 0; local < half_label; ++local) {
            const int left = block * half_label + local;
            const int right = (block + 1) * half_label + local;
            result[static_cast<std::size_t>(left)]
                  [static_cast<std::size_t>(right)] = 1;
            result[static_cast<std::size_t>(right)]
                  [static_cast<std::size_t>(left)] = 1;
            if (local + 1 < half_label) {
                result[static_cast<std::size_t>(left)]
                      [static_cast<std::size_t>(right + 1)] = -1;
                result[static_cast<std::size_t>(right + 1)]
                      [static_cast<std::size_t>(left)] = -1;
            }
        }
    }
    return result;
}

}  // namespace

int main() {
    try {
        int core_checks = 0;
        int difference_checks = 0;
        for (int half_label = 1; half_label <= 12; ++half_label) {
            for (int half_level = 2 * half_label + 1;
                 half_level <= 80;
                 ++half_level) {
                const int available = half_level - 2 * half_label + 1;
                const int block_count = available / half_label;
                if (block_count < 2) {
                    continue;
                }
                const Matrix core = core_fusion_matrix(
                    half_level, half_label, block_count
                );
                require(
                    core == expected_core_matrix(half_label, block_count),
                    "block-Jacobi core mismatch at K="
                        + std::to_string(half_level)
                        + " Q=" + std::to_string(half_label)
                );
                ++core_checks;

                const Matrix difference = block_difference(
                    half_label, block_count
                );
                const Matrix transformed = multiply(
                    multiply(difference, core), transpose(difference)
                );
                require(
                    transformed == expected_difference_form(
                        half_label, block_count
                    ),
                    "block difference form mismatch at K="
                        + std::to_string(half_level)
                        + " Q=" + std::to_string(half_label)
                );
                ++difference_checks;
            }
        }
        std::cout
            << "SU2_BLOCK_JACOBI_CORE"
            << " half_levels<=80"
            << " half_labels<=12"
            << " core_checks=" << core_checks
            << " difference_checks=" << difference_checks
            << " result=PASS_LOCAL_IDENTITIES"
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
