#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Matrix = std::vector<std::vector<std::int64_t>>;

int parse_nonnegative(const char* text, const char* name) {
    const std::string input(text);
    std::size_t used = 0U;
    const long long value = std::stoll(input, &used, 10);
    if (used != input.size() || value < 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
}

std::int64_t mod_normalize(std::int64_t value) {
    constexpr std::int64_t modulus = 1'000'000'007;
    value %= modulus;
    if (value < 0) {
        value += modulus;
    }
    return value;
}

std::int64_t mod_power(std::int64_t base, std::int64_t exponent) {
    constexpr std::int64_t modulus = 1'000'000'007;
    std::int64_t result = 1;
    base = mod_normalize(base);
    while (exponent > 0) {
        if ((exponent & 1) != 0) {
            result = result * base % modulus;
        }
        base = base * base % modulus;
        exponent >>= 1;
    }
    return result;
}

std::size_t modular_rank(Matrix matrix) {
    constexpr std::int64_t modulus = 1'000'000'007;
    if (matrix.empty()) {
        return 0U;
    }
    const std::size_t rows = matrix.size();
    const std::size_t columns = matrix.front().size();
    std::size_t pivot_row = 0U;
    for (std::size_t column = 0U;
         column < columns && pivot_row < rows; ++column) {
        std::size_t pivot = pivot_row;
        while (pivot < rows
               && mod_normalize(matrix[pivot][column]) == 0) {
            ++pivot;
        }
        if (pivot == rows) {
            continue;
        }
        std::swap(matrix[pivot_row], matrix[pivot]);
        const std::int64_t inverse = mod_power(
            matrix[pivot_row][column], modulus - 2
        );
        for (std::size_t entry = column; entry < columns; ++entry) {
            matrix[pivot_row][entry]
                = mod_normalize(
                    matrix[pivot_row][entry] * inverse
                );
        }
        for (std::size_t row = 0U; row < rows; ++row) {
            if (row == pivot_row) {
                continue;
            }
            const std::int64_t factor
                = mod_normalize(matrix[row][column]);
            if (factor == 0) {
                continue;
            }
            for (std::size_t entry = column;
                 entry < columns; ++entry) {
                matrix[row][entry] = mod_normalize(
                    matrix[row][entry]
                    - factor * matrix[pivot_row][entry]
                );
            }
        }
        ++pivot_row;
    }
    return pivot_row;
}

Matrix two_job_transfer(
    int incoming_total,
    int alpha_a,
    int beta_a,
    int sign_a,
    int alpha_b,
    int beta_b,
    int sign_b
) {
    const int delta_a = alpha_a - beta_a;
    const int delta_b = alpha_b - beta_b;
    const int outgoing_total
        = incoming_total + delta_a + delta_b;
    if (incoming_total < 0 || outgoing_total < 0) {
        throw std::runtime_error("negative local total height");
    }
    Matrix transfer(
        static_cast<std::size_t>(incoming_total + 1),
        std::vector<std::int64_t>(
            static_cast<std::size_t>(outgoing_total + 1), 0
        )
    );
    for (int incoming_cut = 0;
         incoming_cut <= incoming_total; ++incoming_cut) {
        std::map<int, std::int64_t> after_a;
        if (incoming_total - incoming_cut >= beta_a) {
            ++after_a[incoming_cut];
        }
        if (incoming_cut >= beta_a) {
            after_a[incoming_cut + delta_a] += sign_a;
        }
        const int middle_total = incoming_total + delta_a;
        for (const auto& [middle_cut, weight] : after_a) {
            if (middle_total - middle_cut >= beta_b) {
                transfer[static_cast<std::size_t>(incoming_cut)]
                        [static_cast<std::size_t>(middle_cut)]
                    += weight;
            }
            if (middle_cut >= beta_b) {
                const int outgoing_cut = middle_cut + delta_b;
                transfer[static_cast<std::size_t>(incoming_cut)]
                        [static_cast<std::size_t>(outgoing_cut)]
                    += static_cast<std::int64_t>(sign_b) * weight;
            }
        }
    }
    return transfer;
}

Matrix difference(const Matrix& first, const Matrix& second) {
    if (first.size() != second.size()
        || (!first.empty()
            && first.front().size() != second.front().size())) {
        throw std::runtime_error("local matrices have different shapes");
    }
    Matrix result = first;
    for (std::size_t row = 0U; row < result.size(); ++row) {
        for (std::size_t column = 0U;
             column < result[row].size(); ++column) {
            result[row][column] -= second[row][column];
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_bk_local_operator"
                " MAXIMUM_QUEUE_PART MAXIMUM_INCOMING_HEIGHT"
            );
        }
        const int part_bound = parse_nonnegative(
            argv[1], "maximum queue part"
        );
        const int height_bound = parse_nonnegative(
            argv[2], "maximum incoming height"
        );
        std::uint64_t cases = 0U;
        std::map<std::size_t, std::uint64_t> signed_rank_histogram;
        std::map<std::size_t, std::uint64_t> unsigned_rank_histogram;
        std::size_t maximum_signed_rank = 0U;
        std::size_t maximum_unsigned_rank = 0U;
        std::uint64_t unsigned_nonnegative_defects = 0U;
        std::uint64_t unsigned_nonpositive_defects = 0U;
        std::uint64_t unsigned_zero_defects = 0U;
        std::uint64_t unsigned_mixed_defects = 0U;
        std::uint64_t unsigned_delta_orientation_failures = 0U;
        std::uint64_t unsigned_determinant_orientation_failures = 0U;
        std::uint64_t unsigned_threshold_orientation_failures = 0U;
        std::uint64_t scheduling_tie_classification_failures = 0U;
        std::vector<int> unsigned_delta_orientation_witness;
        std::vector<int> unsigned_determinant_orientation_witness;
        std::vector<int> unsigned_threshold_orientation_witness;
        std::vector<int> signed_witness;
        std::vector<int> unsigned_witness;
        std::vector<int> unsigned_mixed_witness;
        for (int height = 0; height <= height_bound; ++height) {
            for (int alpha_a = 0; alpha_a <= part_bound; ++alpha_a) {
                for (int beta_a = 0; beta_a <= part_bound; ++beta_a) {
                    if (height < beta_a) {
                        continue;
                    }
                    const int middle
                        = height + alpha_a - beta_a;
                    for (int alpha_b = 0;
                         alpha_b <= part_bound; ++alpha_b) {
                        for (int beta_b = 0;
                             beta_b <= part_bound; ++beta_b) {
                            if (middle < beta_b) {
                                continue;
                            }
                            const int z = std::max(
                                0, beta_a + beta_b - height
                            );
                            const int switched_alpha_b
                                = alpha_b + z;
                            const int switched_beta_b = beta_b - z;
                            const int switched_alpha_a
                                = alpha_a - z;
                            const int switched_beta_a = beta_a + z;
                            if (switched_alpha_a < 0
                                || switched_beta_b < 0) {
                                throw std::runtime_error(
                                    "invalid switched local job"
                                );
                            }
                            for (const int sign_a : {-1, 1}) {
                                for (const int sign_b : {-1, 1}) {
                                    const Matrix original
                                        = two_job_transfer(
                                            height,
                                            alpha_a, beta_a, sign_a,
                                            alpha_b, beta_b, sign_b
                                        );
                                    const Matrix switched
                                        = two_job_transfer(
                                            height,
                                            switched_alpha_b,
                                            switched_beta_b,
                                            sign_b,
                                            switched_alpha_a,
                                            switched_beta_a,
                                            sign_a
                                        );
                                    const Matrix defect = difference(
                                        original, switched
                                    );
                                    const std::size_t signed_rank
                                        = modular_rank(defect);
                                    ++signed_rank_histogram[signed_rank];
                                    if (signed_rank
                                        > maximum_signed_rank) {
                                        maximum_signed_rank = signed_rank;
                                        signed_witness = {
                                            height, alpha_a, beta_a,
                                            alpha_b, beta_b,
                                            sign_a, sign_b, z
                                        };
                                    }
                                    if (sign_a == 1 && sign_b == 1) {
                                        const int original_threshold
                                            = std::max(
                                                beta_a,
                                                beta_b
                                                  - (alpha_a - beta_a)
                                            );
                                        const int switched_threshold
                                            = std::max(
                                                switched_beta_b,
                                                switched_beta_a
                                                  - (switched_alpha_b
                                                     - switched_beta_b)
                                            );
                                        if (z == 0) {
                                            const bool productive_a
                                                = alpha_a >= beta_a;
                                            const bool productive_b
                                                = alpha_b >= beta_b;
                                            const bool scheduled
                                                = productive_a
                                                  != productive_b
                                                ? productive_a
                                                : productive_a
                                                  ? beta_a <= beta_b
                                                  : alpha_a >= alpha_b;
                                            if (scheduled) {
                                                const bool predicted_tie
                                                    = productive_a
                                                      && productive_b
                                                    ? beta_a == beta_b
                                                      || alpha_a == beta_a
                                                    : !productive_a
                                                      && !productive_b
                                                    ? alpha_a == alpha_b
                                                    : alpha_a == beta_a
                                                      && beta_a <= alpha_b;
                                                const bool actual_tie
                                                    = original_threshold
                                                      == switched_threshold;
                                                if (predicted_tie
                                                    != actual_tie) {
                                                    ++scheduling_tie_classification_failures;
                                                }
                                            }
                                        }
                                        bool has_positive = false;
                                        bool has_negative = false;
                                        for (const auto& row : defect) {
                                            for (const std::int64_t entry
                                                 : row) {
                                                has_positive
                                                    = has_positive
                                                      || entry > 0;
                                                has_negative
                                                    = has_negative
                                                      || entry < 0;
                                            }
                                        }
                                        if (!has_positive
                                            && !has_negative) {
                                            ++unsigned_zero_defects;
                                        } else if (has_positive
                                                   && has_negative) {
                                            ++unsigned_mixed_defects;
                                            if (unsigned_mixed_witness
                                                .empty()) {
                                                unsigned_mixed_witness = {
                                                    height, alpha_a,
                                                    beta_a, alpha_b,
                                                    beta_b, z
                                                };
                                            }
                                        } else if (has_positive) {
                                            ++unsigned_nonnegative_defects;
                                            if (original_threshold
                                                >= switched_threshold) {
                                                ++unsigned_threshold_orientation_failures;
                                                if (unsigned_threshold_orientation_witness
                                                    .empty()) {
                                                    unsigned_threshold_orientation_witness = {
                                                        height,
                                                        original_threshold,
                                                        switched_threshold,
                                                        1
                                                    };
                                                }
                                            }
                                            const int determinant
                                                = alpha_a * beta_b
                                                  - alpha_b * beta_a;
                                            if (z == 0
                                                && determinant < 0) {
                                                ++unsigned_determinant_orientation_failures;
                                                if (unsigned_determinant_orientation_witness
                                                    .empty()) {
                                                    unsigned_determinant_orientation_witness = {
                                                        height, alpha_a,
                                                        beta_a, alpha_b,
                                                        beta_b, z, 1
                                                    };
                                                }
                                            }
                                            if (alpha_a - beta_a
                                                <= alpha_b - beta_b) {
                                                ++unsigned_delta_orientation_failures;
                                                if (unsigned_delta_orientation_witness
                                                    .empty()) {
                                                    unsigned_delta_orientation_witness = {
                                                        height, alpha_a,
                                                        beta_a, alpha_b,
                                                        beta_b, z, 1
                                                    };
                                                }
                                            }
                                        } else {
                                            ++unsigned_nonpositive_defects;
                                            if (original_threshold
                                                <= switched_threshold) {
                                                ++unsigned_threshold_orientation_failures;
                                                if (unsigned_threshold_orientation_witness
                                                    .empty()) {
                                                    unsigned_threshold_orientation_witness = {
                                                        height,
                                                        original_threshold,
                                                        switched_threshold,
                                                        -1
                                                    };
                                                }
                                            }
                                            const int determinant
                                                = alpha_a * beta_b
                                                  - alpha_b * beta_a;
                                            if (z == 0
                                                && determinant > 0) {
                                                ++unsigned_determinant_orientation_failures;
                                                if (unsigned_determinant_orientation_witness
                                                    .empty()) {
                                                    unsigned_determinant_orientation_witness = {
                                                        height, alpha_a,
                                                        beta_a, alpha_b,
                                                        beta_b, z, -1
                                                    };
                                                }
                                            }
                                            if (alpha_a - beta_a
                                                >= alpha_b - beta_b) {
                                                ++unsigned_delta_orientation_failures;
                                                if (unsigned_delta_orientation_witness
                                                    .empty()) {
                                                    unsigned_delta_orientation_witness = {
                                                        height, alpha_a,
                                                        beta_a, alpha_b,
                                                        beta_b, z, -1
                                                    };
                                                }
                                            }
                                        }
                                        ++unsigned_rank_histogram[
                                            signed_rank
                                        ];
                                        if (signed_rank
                                            > maximum_unsigned_rank) {
                                            maximum_unsigned_rank
                                                = signed_rank;
                                            unsigned_witness = {
                                                height, alpha_a, beta_a,
                                                alpha_b, beta_b, z
                                            };
                                        }
                                    }
                                    ++cases;
                                }
                            }
                        }
                    }
                }
            }
        }
        const auto print_histogram = [](
            const std::map<std::size_t, std::uint64_t>& histogram
        ) {
            std::cout << '{';
            bool first = true;
            for (const auto& [rank, multiplicity] : histogram) {
                if (!first) {
                    std::cout << ',';
                }
                first = false;
                std::cout << rank << ':' << multiplicity;
            }
            std::cout << '}';
        };
        std::cout
            << "SU2_BK_LOCAL_OPERATOR cases=" << cases
            << " signed_rank_histogram=";
        print_histogram(signed_rank_histogram);
        std::cout << " unsigned_rank_histogram=";
        print_histogram(unsigned_rank_histogram);
        std::cout
            << " maximum_signed_rank=" << maximum_signed_rank
            << " signed_witness=[";
        for (std::size_t index = 0U;
             index < signed_witness.size(); ++index) {
            if (index != 0U) {
                std::cout << ',';
            }
            std::cout << signed_witness[index];
        }
        std::cout
            << "] maximum_unsigned_rank=" << maximum_unsigned_rank
            << " unsigned_witness=[";
        for (std::size_t index = 0U;
             index < unsigned_witness.size(); ++index) {
            if (index != 0U) {
                std::cout << ',';
            }
            std::cout << unsigned_witness[index];
        }
        std::cout
            << "] unsigned_nonnegative_defects="
            << unsigned_nonnegative_defects
            << " unsigned_nonpositive_defects="
            << unsigned_nonpositive_defects
            << " unsigned_zero_defects=" << unsigned_zero_defects
            << " unsigned_mixed_defects=" << unsigned_mixed_defects
            << " unsigned_delta_orientation_failures="
            << unsigned_delta_orientation_failures
            << " unsigned_determinant_orientation_failures="
            << unsigned_determinant_orientation_failures
            << " unsigned_threshold_orientation_failures="
            << unsigned_threshold_orientation_failures
            << " scheduling_tie_classification_failures="
            << scheduling_tie_classification_failures
            << " unsigned_mixed_witness=[";
        for (std::size_t index = 0U;
             index < unsigned_mixed_witness.size(); ++index) {
            if (index != 0U) {
                std::cout << ',';
            }
            std::cout << unsigned_mixed_witness[index];
        }
        std::cout << "] unsigned_delta_orientation_witness=[";
        for (std::size_t index = 0U;
             index < unsigned_delta_orientation_witness.size();
             ++index) {
            if (index != 0U) {
                std::cout << ',';
            }
            std::cout << unsigned_delta_orientation_witness[index];
        }
        std::cout << "] unsigned_determinant_orientation_witness=[";
        for (std::size_t index = 0U;
             index < unsigned_determinant_orientation_witness.size();
             ++index) {
            if (index != 0U) {
                std::cout << ',';
            }
            std::cout << unsigned_determinant_orientation_witness[
                index
            ];
        }
        std::cout << "] unsigned_threshold_orientation_witness=[";
        for (std::size_t index = 0U;
             index < unsigned_threshold_orientation_witness.size();
             ++index) {
            if (index != 0U) {
                std::cout << ',';
            }
            std::cout << unsigned_threshold_orientation_witness[
                index
            ];
        }
        std::cout << "] result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
