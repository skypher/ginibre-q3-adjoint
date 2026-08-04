#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Matrix = std::vector<std::vector<Integer>>;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(const bool condition, const std::string& message) {
    if (!condition) {
        fail(message);
    }
}

Matrix fusion_matrix(const int level, const int label) {
    const std::size_t size = static_cast<std::size_t>(level + 1);
    Matrix result(size, std::vector<Integer>(size));
    for (int source = 0; source <= level; ++source) {
        const int lower = std::abs(source - label);
        const int upper = std::min(
            source + label,
            2 * level - source - label
        );
        for (int target = lower; target <= upper; target += 2) {
            result[static_cast<std::size_t>(target)]
                  [static_cast<std::size_t>(source)] = 1;
        }
    }
    return result;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    require(left.size() == right.size(), "matrix-size mismatch");
    const std::size_t size = left.size();
    Matrix result(size, std::vector<Integer>(size));
    for (std::size_t row = 0U; row < size; ++row) {
        require(left[row].size() == size, "left matrix is not square");
        for (std::size_t middle = 0U; middle < size; ++middle) {
            const Integer& coefficient = left[row][middle];
            if (coefficient == 0) {
                continue;
            }
            require(
                right[middle].size() == size,
                "right matrix is not square"
            );
            for (std::size_t column = 0U; column < size; ++column) {
                result[row][column] += coefficient * right[middle][column];
            }
        }
    }
    return result;
}

Matrix add(const Matrix& left, const Matrix& right) {
    require(left.size() == right.size(), "matrix-size mismatch");
    Matrix result = left;
    for (std::size_t row = 0U; row < result.size(); ++row) {
        require(right[row].size() == result[row].size(), "row-size mismatch");
        for (std::size_t column = 0U; column < result[row].size(); ++column) {
            result[row][column] += right[row][column];
        }
    }
    return result;
}

Matrix identity(const int level) {
    const std::size_t size = static_cast<std::size_t>(level + 1);
    Matrix result(size, std::vector<Integer>(size));
    for (std::size_t index = 0U; index < size; ++index) {
        result[index][index] = 1;
    }
    return result;
}

void verify_recurrence(const int level, const int label) {
    const Matrix fundamental = fusion_matrix(level, 1);
    const Matrix predecessor = fusion_matrix(level, label - 1);
    const Matrix expected = add(
        fusion_matrix(level, label),
        fusion_matrix(level, label - 2)
    );
    require(
        multiply(fundamental, predecessor) == expected,
        "fundamental fusion recurrence failure at level="
            + std::to_string(level)
            + " label=" + std::to_string(label)
    );
}

void verify_mixed_incidence_obstruction(const int level, const int label) {
    const Matrix sum = add(
        fusion_matrix(level, label),
        fusion_matrix(level, label - 2)
    );
    const std::size_t vacuum = 0U;
    const std::size_t neighbour = static_cast<std::size_t>(label - 2);
    require(sum[vacuum][vacuum] == 0, "nonzero vacuum diagonal");
    require(sum[vacuum][neighbour] == 1, "wrong vacuum off-diagonal");
    const Integer determinant = sum[vacuum][vacuum] * sum[neighbour][neighbour]
        - sum[vacuum][neighbour] * sum[neighbour][vacuum];
    require(determinant == -1, "wrong obstructing principal minor");

    const Matrix mixed_incidence = fusion_matrix(level, label - 1);
    const std::size_t lower_row = static_cast<std::size_t>(label - 3);
    const std::size_t upper_row = static_cast<std::size_t>(label - 1);
    const std::size_t first_column = 0U;
    const std::size_t second_column = 2U;
    const Integer mixed_determinant = mixed_incidence[lower_row][first_column]
            * mixed_incidence[upper_row][second_column]
        - mixed_incidence[lower_row][second_column]
            * mixed_incidence[upper_row][first_column];
    require(mixed_determinant == -1, "wrong mixed-incidence minor");
}

void verify_minimal_gram_identity(const int level) {
    const Matrix fundamental = fusion_matrix(level, 1);
    const Matrix minimal_sum = add(fusion_matrix(level, 2), identity(level));
    require(
        multiply(fundamental, fundamental) == minimal_sum,
        "minimal-label Gram identity failure at level=" + std::to_string(level)
    );
}

}  // namespace

int main() {
    try {
        int recurrence_checks = 0;
        int obstruction_checks = 0;
        int minimal_gram_checks = 0;
        for (int level = 4; level <= 32; ++level) {
            verify_minimal_gram_identity(level);
            ++minimal_gram_checks;
            for (int label = 2; label <= level; ++label) {
                verify_recurrence(level, label);
                ++recurrence_checks;
                if (label >= 4) {
                    verify_mixed_incidence_obstruction(level, label);
                    ++obstruction_checks;
                }
            }
        }
        std::cout
            << "SU2_MIXED_INCIDENCE_OBSTRUCTION"
            << " levels=4..32"
            << " recurrence_checks=" << recurrence_checks
            << " minimal_gram_checks=" << minimal_gram_checks
            << " obstruction_checks=" << obstruction_checks
            << " result=PASS_LOCAL_IDENTITIES"
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
