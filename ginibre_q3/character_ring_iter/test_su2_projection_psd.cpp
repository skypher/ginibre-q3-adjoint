#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

extern "C" void dsyev_(
    const char* jobz,
    const char* uplo,
    const int* n,
    double* a,
    const int* lda,
    double* w,
    double* work,
    const int* lwork,
    int* info
);

namespace {

using Matrix = std::vector<double>;

double& entry(Matrix& matrix, int dimension, int row, int column) {
    return matrix[static_cast<std::size_t>(row + dimension * column)];
}

double entry(const Matrix& matrix, int dimension, int row, int column) {
    return matrix[static_cast<std::size_t>(row + dimension * column)];
}

void diagonalize(Matrix& matrix, int dimension, std::vector<double>& eigenvalues) {
    eigenvalues.resize(static_cast<std::size_t>(dimension));
    const char jobz = 'V';
    const char uplo = 'U';
    const int lda = dimension;
    int lwork = -1;
    double workspace_query = 0.0;
    int info = 0;
    dsyev_(
        &jobz, &uplo, &dimension, matrix.data(), &lda, eigenvalues.data(),
        &workspace_query, &lwork, &info
    );
    if (info != 0) {
        throw std::runtime_error("LAPACK workspace query failed");
    }
    lwork = std::max(1, static_cast<int>(workspace_query));
    std::vector<double> workspace(static_cast<std::size_t>(lwork));
    dsyev_(
        &jobz, &uplo, &dimension, matrix.data(), &lda, eigenvalues.data(),
        workspace.data(), &lwork, &info
    );
    if (info != 0) {
        throw std::runtime_error("LAPACK eigensolver failed");
    }
}

Matrix singlet_projector(int factors) {
    const int dimension = 1 << factors;
    Matrix spin_squared(
        static_cast<std::size_t>(dimension * dimension), 0.0
    );
    for (int state = 0; state < dimension; ++state) {
        const int up = std::popcount(static_cast<unsigned int>(state));
        const double magnetic = static_cast<double>(2 * up - factors) / 2.0;
        entry(spin_squared, dimension, state, state) =
            magnetic * magnetic + static_cast<double>(factors) / 2.0;
        for (int first = 0; first < factors; ++first) {
            for (int second = first + 1; second < factors; ++second) {
                const int first_bit = (state >> first) & 1;
                const int second_bit = (state >> second) & 1;
                if (first_bit == second_bit) {
                    continue;
                }
                const int flipped = state ^ (1 << first) ^ (1 << second);
                entry(spin_squared, dimension, flipped, state) += 1.0;
            }
        }
    }

    std::vector<double> eigenvalues;
    diagonalize(spin_squared, dimension, eigenvalues);
    Matrix projector(static_cast<std::size_t>(dimension * dimension), 0.0);
    for (int vector = 0; vector < dimension; ++vector) {
        if (std::abs(eigenvalues[static_cast<std::size_t>(vector)]) > 1e-8) {
            continue;
        }
        for (int row = 0; row < dimension; ++row) {
            for (int column = 0; column < dimension; ++column) {
                entry(projector, dimension, row, column) +=
                    entry(spin_squared, dimension, row, vector)
                    * entry(spin_squared, dimension, column, vector);
            }
        }
    }
    return projector;
}

int compressed_bits(int state, int mask, int factors) {
    int result = 0;
    int output_bit = 0;
    for (int bit = 0; bit < factors; ++bit) {
        if (((mask >> bit) & 1) == 0) {
            continue;
        }
        result |= ((state >> bit) & 1) << output_bit;
        ++output_bit;
    }
    return result;
}

Matrix signed_projection_sum(
    int factors,
    int minus_factors,
    const std::vector<Matrix>& projectors
) {
    const int dimension = 1 << factors;
    const int full_mask = dimension - 1;
    const int minus_mask = (1 << minus_factors) - 1;
    Matrix answer(static_cast<std::size_t>(dimension * dimension), 0.0);
    for (int subset = 0; subset < dimension; ++subset) {
        const int subset_size = std::popcount(static_cast<unsigned int>(subset));
        if ((subset_size & 1) != 0 || ((factors - subset_size) & 1) != 0) {
            continue;
        }
        const int complement = full_mask ^ subset;
        const int sign = (std::popcount(
            static_cast<unsigned int>(subset & minus_mask)
        ) & 1) == 0 ? 1 : -1;
        const int subset_dimension = 1 << subset_size;
        const int complement_dimension = 1 << (factors - subset_size);
        const Matrix& subset_projector =
            projectors[static_cast<std::size_t>(subset_size)];
        const Matrix& complement_projector =
            projectors[static_cast<std::size_t>(factors - subset_size)];
        for (int row = 0; row < dimension; ++row) {
            const int subset_row = compressed_bits(row, subset, factors);
            const int complement_row = compressed_bits(row, complement, factors);
            for (int column = 0; column < dimension; ++column) {
                const int subset_column = compressed_bits(column, subset, factors);
                const int complement_column = compressed_bits(
                    column, complement, factors
                );
                entry(answer, dimension, row, column) += static_cast<double>(sign)
                    * entry(
                        subset_projector,
                        subset_dimension,
                        subset_row,
                        subset_column
                    )
                    * entry(
                        complement_projector,
                        complement_dimension,
                        complement_row,
                        complement_column
                    );
            }
        }
    }
    return answer;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error("usage: test_su2_projection_psd MAXIMUM_FACTORS");
        }
        const int maximum_factors = std::stoi(argv[1]);
        if (maximum_factors < 2 || maximum_factors > 10) {
            throw std::runtime_error("MAXIMUM_FACTORS must lie between 2 and 10");
        }
        std::vector<Matrix> projectors(
            static_cast<std::size_t>(maximum_factors + 1)
        );
        for (int factors = 0; factors <= maximum_factors; factors += 2) {
            projectors[static_cast<std::size_t>(factors)] =
                singlet_projector(factors);
        }
        std::cout << std::setprecision(17);
        for (int factors = 2; factors <= maximum_factors; factors += 2) {
            const int dimension = 1 << factors;
            for (int minus_factors = 2;
                 minus_factors <= factors;
                 minus_factors += 2) {
                Matrix sum = signed_projection_sum(
                    factors, minus_factors, projectors
                );
                std::vector<double> eigenvalues;
                diagonalize(sum, dimension, eigenvalues);
                std::cout << "factors=" << factors
                          << " minus=" << minus_factors
                          << " minimum_eigenvalue=" << eigenvalues.front()
                          << " maximum_eigenvalue=" << eigenvalues.back() << '\n';
                if (eigenvalues.front() < -1e-8) {
                    std::cout << "FAIL\n";
                    return EXIT_FAILURE;
                }
            }
        }
        std::cout << "PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
