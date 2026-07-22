#include <array>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

constexpr int kRank = 3;
constexpr int kTensorDimension = kRank * kRank;
constexpr int kSymmetricDimension = 6;

using Matrix = std::array<long long, kTensorDimension * kTensorDimension>;
using SymmetricMatrix
    = std::array<long long, kSymmetricDimension * kSymmetricDimension>;
using Vector = std::array<long long, kTensorDimension>;

constexpr std::array<std::array<int, 2>, kSymmetricDimension> kPairs{{
    {{0, 0}}, {{1, 0}}, {{2, 0}}, {{1, 1}}, {{2, 1}}, {{2, 2}}
}};

long long& entry(Matrix& matrix, int row, int column) {
    return matrix[static_cast<std::size_t>(
        row * kTensorDimension + column
    )];
}

long long entry(const Matrix& matrix, int row, int column) {
    return matrix[static_cast<std::size_t>(
        row * kTensorDimension + column
    )];
}

Matrix identity() {
    Matrix result{};
    for (int index = 0; index < kTensorDimension; ++index) {
        entry(result, index, index) = 1;
    }
    return result;
}

Matrix add(const Matrix& left, const Matrix& right, long long right_scale) {
    Matrix result{};
    for (std::size_t index = 0; index < result.size(); ++index) {
        result[index] = left[index] + right_scale * right[index];
    }
    return result;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    Matrix result{};
    for (int row = 0; row < kTensorDimension; ++row) {
        for (int inner = 0; inner < kTensorDimension; ++inner) {
            const long long coefficient = entry(left, row, inner);
            if (coefficient == 0) {
                continue;
            }
            for (int column = 0; column < kTensorDimension; ++column) {
                entry(result, row, column)
                    += coefficient * entry(right, inner, column);
            }
        }
    }
    return result;
}

Matrix tensor_left(const std::array<long long, 9>& factor) {
    Matrix result{};
    for (int left_output = 0; left_output < kRank; ++left_output) {
        for (int left_input = 0; left_input < kRank; ++left_input) {
            const long long coefficient = factor[static_cast<std::size_t>(
                left_output * kRank + left_input
            )];
            for (int right = 0; right < kRank; ++right) {
                entry(
                    result,
                    left_output * kRank + right,
                    left_input * kRank + right
                ) = coefficient;
            }
        }
    }
    return result;
}

Matrix tensor_right(const std::array<long long, 9>& factor) {
    Matrix result{};
    for (int right_output = 0; right_output < kRank; ++right_output) {
        for (int right_input = 0; right_input < kRank; ++right_input) {
            const long long coefficient = factor[static_cast<std::size_t>(
                right_output * kRank + right_input
            )];
            for (int left = 0; left < kRank; ++left) {
                entry(
                    result,
                    left * kRank + right_output,
                    left * kRank + right_input
                ) = coefficient;
            }
        }
    }
    return result;
}

Vector symmetric_basis_vector(int index) {
    Vector result{};
    const auto [first, second] = kPairs[static_cast<std::size_t>(index)];
    result[static_cast<std::size_t>(first * kRank + second)] = 1;
    if (first != second) {
        result[static_cast<std::size_t>(second * kRank + first)] = 1;
    }
    return result;
}

Vector apply_matrix(const Matrix& matrix, const Vector& vector) {
    Vector result{};
    for (int row = 0; row < kTensorDimension; ++row) {
        for (int column = 0; column < kTensorDimension; ++column) {
            result[static_cast<std::size_t>(row)]
                += entry(matrix, row, column)
                    * vector[static_cast<std::size_t>(column)];
        }
    }
    return result;
}

SymmetricMatrix restrict_symmetric(const Matrix& matrix) {
    SymmetricMatrix result{};
    for (int column = 0; column < kSymmetricDimension; ++column) {
        const Vector image = apply_matrix(
            matrix, symmetric_basis_vector(column)
        );
        for (int row = 0; row < kSymmetricDimension; ++row) {
            const auto [first, second]
                = kPairs[static_cast<std::size_t>(row)];
            const long long first_value = image[static_cast<std::size_t>(
                first * kRank + second
            )];
            const long long second_value = image[static_cast<std::size_t>(
                second * kRank + first
            )];
            if (first_value != second_value) {
                std::cerr << "operator failed to preserve symmetry\n";
                std::exit(EXIT_FAILURE);
            }
            result[static_cast<std::size_t>(
                row * kSymmetricDimension + column
            )] = first_value;
        }
    }
    return result;
}

void print_matrix(const std::string& name, const Matrix& matrix) {
    const SymmetricMatrix restricted = restrict_symmetric(matrix);
    std::cout << name << '\n';
    for (int row = 0; row < kSymmetricDimension; ++row) {
        for (int column = 0; column < kSymmetricDimension; ++column) {
            if (column != 0) {
                std::cout << ' ';
            }
            std::cout << restricted[static_cast<std::size_t>(
                row * kSymmetricDimension + column
            )];
        }
        std::cout << '\n';
    }
}

}  // namespace

int main() {
    const std::array<long long, 9> first{{
        0, 1, 0,
        1, 1, 1,
        0, 1, 1
    }};
    const std::array<long long, 9> second{{
        0, 0, 1,
        0, 1, 1,
        1, 1, 0
    }};

    const Matrix first_left = tensor_left(first);
    const Matrix first_right = tensor_right(first);
    const Matrix second_left = tensor_left(second);
    const Matrix second_right = tensor_right(second);
    const Matrix plus_first = add(first_left, first_right, 1);
    const Matrix minus_first = add(first_left, first_right, -1);
    const Matrix plus_second = add(second_left, second_right, 1);
    const Matrix minus_second = add(second_left, second_right, -1);
    const Matrix first_minus_pair = multiply(minus_first, minus_first);
    const Matrix second_minus_pair = multiply(minus_second, minus_second);
    const Matrix mixed_minus_pair = multiply(minus_first, minus_second);

    print_matrix("U=T1_plus", plus_first);
    print_matrix("Q=T2_plus", plus_second);
    print_matrix("P=T1_minus_squared", first_minus_pair);
    print_matrix("R=T2_minus_squared", second_minus_pair);
    print_matrix("T=T1_minus_T2_minus", mixed_minus_pair);

    const Matrix unit = identity();
    const Matrix relation = add(
        second_left,
        add(multiply(first_left, first_left), first_left, -1),
        -1
    );
    const Matrix relation_with_unit = add(relation, unit, 1);
    for (long long coefficient : relation_with_unit) {
        if (coefficient != 0) {
            std::cerr << "rank-three fusion relation failed\n";
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
