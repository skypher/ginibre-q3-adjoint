#include "full_q3_bcd_remaining_data.hpp"

#include <gmpxx.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using BigInt = mpz_class;

constexpr int kCertifiedMaximum =
    full_q3_bcd_remaining::required_maximum_moment;
constexpr int kMaximumSeriesDegree = 2 * kCertifiedMaximum;
constexpr int kMaximumNodes = kCertifiedMaximum + 1;

// All arithmetic in the bounded-Littlewood determinant stage is performed in
// exact prime fields.  The residues are reconstructed with GMP only after the
// product of the primes exceeds the elementary character bound dim(G)^j.
class MontgomeryField {
public:
    explicit MontgomeryField(std::uint32_t modulus)
        : modulus_(modulus), negative_inverse_(negative_inverse(modulus)) {
        if (modulus_ < 3U || modulus_ % 2U == 0U) {
            throw std::runtime_error("invalid Montgomery modulus");
        }
        const std::uint64_t r_mod =
            (std::uint64_t{1} << 32U) % modulus_;
        r_squared_ = static_cast<std::uint32_t>((r_mod * r_mod) % modulus_);
        one_ = to_montgomery(1U);
    }

    [[nodiscard]] std::uint32_t modulus() const { return modulus_; }
    [[nodiscard]] std::uint32_t zero() const { return 0U; }
    [[nodiscard]] std::uint32_t one() const { return one_; }

    [[nodiscard]] std::uint32_t add(
        std::uint32_t left,
        std::uint32_t right
    ) const {
        const std::uint64_t sum =
            static_cast<std::uint64_t>(left) + right;
        return static_cast<std::uint32_t>(
            sum >= modulus_ ? sum - modulus_ : sum
        );
    }

    [[nodiscard]] std::uint32_t subtract(
        std::uint32_t left,
        std::uint32_t right
    ) const {
        return left >= right
            ? left - right
            : static_cast<std::uint32_t>(
                  static_cast<std::uint64_t>(left) + modulus_ - right
              );
    }

    [[nodiscard]] std::uint32_t negate(std::uint32_t value) const {
        return value == 0U ? 0U : modulus_ - value;
    }

    [[nodiscard]] std::uint32_t multiply(
        std::uint32_t left,
        std::uint32_t right
    ) const {
        return reduce(static_cast<std::uint64_t>(left) * right);
    }

    [[nodiscard]] std::uint32_t power(
        std::uint32_t base,
        std::uint64_t exponent
    ) const {
        std::uint32_t answer = one_;
        while (exponent > 0U) {
            if ((exponent & 1U) != 0U) answer = multiply(answer, base);
            exponent >>= 1U;
            if (exponent > 0U) base = multiply(base, base);
        }
        return answer;
    }

    [[nodiscard]] std::uint32_t inverse(std::uint32_t value) const {
        if (value == 0U) throw std::runtime_error("division by zero in F_p");
        return power(value, static_cast<std::uint64_t>(modulus_) - 2U);
    }

    [[nodiscard]] std::uint32_t from_unsigned(std::uint64_t value) const {
        const std::uint32_t reduced =
            static_cast<std::uint32_t>(value % modulus_);
        return multiply(reduced, r_squared_);
    }

    [[nodiscard]] std::uint32_t from_signed(std::int64_t value) const {
        if (value >= 0) {
            return from_unsigned(static_cast<std::uint64_t>(value));
        }
        const std::uint64_t magnitude =
            static_cast<std::uint64_t>(-(value + 1)) + 1U;
        return negate(from_unsigned(magnitude));
    }

    [[nodiscard]] std::uint32_t to_unsigned(std::uint32_t value) const {
        return reduce(value);
    }

private:
    static std::uint32_t negative_inverse(std::uint32_t odd) {
        // Newton iteration in Z/2^32 Z.  Six iterations cover all 32 bits.
        std::uint32_t inverse = odd;
        for (int iteration = 0; iteration < 6; ++iteration) {
            inverse *= 2U - odd * inverse;
        }
        return 0U - inverse;
    }

    [[nodiscard]] std::uint32_t reduce(std::uint64_t value) const {
        const std::uint32_t multiplier =
            static_cast<std::uint32_t>(value) * negative_inverse_;
        const __uint128_t sum = static_cast<__uint128_t>(value)
            + static_cast<__uint128_t>(multiplier) * modulus_;
        std::uint64_t quotient = static_cast<std::uint64_t>(sum >> 32U);
        if (quotient >= modulus_) quotient -= modulus_;
        return static_cast<std::uint32_t>(quotient);
    }

    [[nodiscard]] std::uint32_t to_montgomery(std::uint32_t value) const {
        return multiply(value, r_squared_);
    }

    std::uint32_t modulus_;
    std::uint32_t negative_inverse_;
    std::uint32_t r_squared_ = 0U;
    std::uint32_t one_ = 0U;
};

struct Series {
    std::array<std::uint32_t, kMaximumSeriesDegree + 1> coefficient{};
    int lower = kMaximumSeriesDegree + 1;
    int upper = -1;

    [[nodiscard]] bool is_zero() const { return upper < lower; }
};

void normalize(Series& series, int degree_limit) {
    series.lower = 0;
    while (series.lower <= degree_limit
           && series.coefficient[static_cast<std::size_t>(series.lower)] == 0U) {
        ++series.lower;
    }
    if (series.lower > degree_limit) {
        series.lower = kMaximumSeriesDegree + 1;
        series.upper = -1;
        return;
    }
    series.upper = degree_limit;
    while (series.upper >= series.lower
           && series.coefficient[static_cast<std::size_t>(series.upper)] == 0U) {
        --series.upper;
    }
}

Series constant_series(
    std::uint32_t value,
    int degree_limit
) {
    Series answer;
    answer.coefficient[0] = value;
    normalize(answer, degree_limit);
    return answer;
}

Series add_series(
    const Series& left,
    const Series& right,
    const MontgomeryField& field,
    int degree_limit,
    bool subtract_right = false
) {
    Series answer;
    for (int degree = 0; degree <= degree_limit; ++degree) {
        const std::size_t index = static_cast<std::size_t>(degree);
        answer.coefficient[index] = subtract_right
            ? field.subtract(left.coefficient[index], right.coefficient[index])
            : field.add(left.coefficient[index], right.coefficient[index]);
    }
    normalize(answer, degree_limit);
    return answer;
}

void subtract_series_in_place(
    Series& left,
    const Series& right,
    const MontgomeryField& field,
    int degree_limit
) {
    if (right.is_zero()) return;
    for (int degree = right.lower; degree <= right.upper; ++degree) {
        const std::size_t index = static_cast<std::size_t>(degree);
        left.coefficient[index] = field.subtract(
            left.coefficient[index], right.coefficient[index]
        );
    }
    normalize(left, degree_limit);
}

Series scale_series(
    const Series& input,
    std::uint32_t scalar,
    const MontgomeryField& field,
    int degree_limit
) {
    if (input.is_zero() || scalar == 0U) return Series{};
    Series answer;
    for (int degree = input.lower; degree <= input.upper; ++degree) {
        const std::size_t index = static_cast<std::size_t>(degree);
        answer.coefficient[index] = field.multiply(
            input.coefficient[index], scalar
        );
    }
    normalize(answer, degree_limit);
    return answer;
}

Series multiply_series(
    const Series& left,
    const Series& right,
    const MontgomeryField& field,
    int degree_limit
) {
    if (left.is_zero() || right.is_zero()
        || left.lower + right.lower > degree_limit) {
        return Series{};
    }
    Series answer;
    const int left_upper = std::min(left.upper, degree_limit - right.lower);
    for (int left_degree = left.lower;
         left_degree <= left_upper;
         ++left_degree) {
        const std::uint32_t left_value =
            left.coefficient[static_cast<std::size_t>(left_degree)];
        if (left_value == 0U) continue;
        const int right_upper = std::min(
            right.upper, degree_limit - left_degree
        );
        for (int right_degree = right.lower;
             right_degree <= right_upper;
             ++right_degree) {
            const std::uint32_t right_value =
                right.coefficient[static_cast<std::size_t>(right_degree)];
            if (right_value == 0U) continue;
            const std::size_t target = static_cast<std::size_t>(
                left_degree + right_degree
            );
            answer.coefficient[target] = field.add(
                answer.coefficient[target],
                field.multiply(left_value, right_value)
            );
        }
    }
    normalize(answer, degree_limit);
    return answer;
}

Series inverse_series(
    const Series& input,
    const MontgomeryField& field,
    int degree_limit
) {
    if (input.coefficient[0] == 0U) {
        throw std::runtime_error("nonunit pivot in bounded-Littlewood determinant");
    }
    Series answer;
    const std::uint32_t inverse_constant = field.inverse(input.coefficient[0]);
    answer.coefficient[0] = inverse_constant;
    for (int degree = 1; degree <= degree_limit; ++degree) {
        std::uint32_t convolution = 0U;
        const int upper = std::min(degree, input.upper);
        for (int index = 1; index <= upper; ++index) {
            convolution = field.add(
                convolution,
                field.multiply(
                    input.coefficient[static_cast<std::size_t>(index)],
                    answer.coefficient[static_cast<std::size_t>(degree - index)]
                )
            );
        }
        answer.coefficient[static_cast<std::size_t>(degree)] = field.negate(
            field.multiply(inverse_constant, convolution)
        );
    }
    normalize(answer, degree_limit);
    return answer;
}

std::vector<Series> leading_principal_determinants(
    std::vector<Series> matrix,
    int size,
    const MontgomeryField& field,
    int degree_limit
) {
    if (size < 1
        || matrix.size() != static_cast<std::size_t>(size * size)) {
        throw std::runtime_error("invalid determinant matrix");
    }
    std::vector<Series> determinants(static_cast<std::size_t>(size + 1));
    determinants[0] = constant_series(field.one(), degree_limit);
    Series cumulative = determinants[0];
    const auto at = [size, &matrix](int row, int column) -> Series& {
        return matrix[static_cast<std::size_t>(row * size + column)];
    };

    for (int pivot_index = 0; pivot_index < size; ++pivot_index) {
        const Series pivot = at(pivot_index, pivot_index);
        if (pivot.coefficient[0] == 0U) {
            throw std::runtime_error(
                "zero constant pivot in bounded-Littlewood determinant"
            );
        }
        cumulative = multiply_series(
            cumulative, pivot, field, degree_limit
        );
        determinants[static_cast<std::size_t>(pivot_index + 1)] = cumulative;
        const Series pivot_inverse = inverse_series(pivot, field, degree_limit);

        for (int row = pivot_index + 1; row < size; ++row) {
            if (at(row, pivot_index).is_zero()) continue;
            const Series factor = multiply_series(
                at(row, pivot_index), pivot_inverse, field, degree_limit
            );
            at(row, pivot_index) = Series{};
            for (int column = pivot_index + 1; column < size; ++column) {
                const Series product = multiply_series(
                    factor,
                    at(pivot_index, column),
                    field,
                    degree_limit
                );
                subtract_series_in_place(
                    at(row, column), product, field, degree_limit
                );
            }
        }
    }
    return determinants;
}

struct NodeData {
    std::vector<std::uint32_t> elementary;
    std::vector<std::uint32_t> complete;
    std::vector<Series> elementary_pair;
    std::vector<Series> complete_pair;
    Series elementary_sum;
    Series alternating_elementary_sum;
};

NodeData build_node_data(
    std::uint32_t node,
    int largest_f_index,
    const MontgomeryField& field,
    int degree_limit
) {
    NodeData data;
    data.elementary.assign(static_cast<std::size_t>(degree_limit + 1), 0U);
    data.complete.assign(static_cast<std::size_t>(degree_limit + 1), 0U);
    data.elementary[0] = field.one();
    data.complete[0] = field.one();
    for (int degree = 1; degree <= degree_limit; ++degree) {
        std::uint32_t elementary_numerator =
            data.elementary[static_cast<std::size_t>(degree - 1)];
        std::uint32_t complete_numerator =
            data.complete[static_cast<std::size_t>(degree - 1)];
        if (degree >= 2) {
            elementary_numerator = field.subtract(
                elementary_numerator,
                field.multiply(
                    node,
                    data.elementary[static_cast<std::size_t>(degree - 2)]
                )
            );
            complete_numerator = field.add(
                complete_numerator,
                field.multiply(
                    node,
                    data.complete[static_cast<std::size_t>(degree - 2)]
                )
            );
        }
        const std::uint32_t degree_inverse = field.inverse(
            field.from_unsigned(static_cast<std::uint64_t>(degree))
        );
        data.elementary[static_cast<std::size_t>(degree)] = field.multiply(
            elementary_numerator, degree_inverse
        );
        data.complete[static_cast<std::size_t>(degree)] = field.multiply(
            complete_numerator, degree_inverse
        );
    }

    data.elementary_sum = Series{};
    data.alternating_elementary_sum = Series{};
    for (int degree = 0; degree <= degree_limit; ++degree) {
        const std::size_t index = static_cast<std::size_t>(degree);
        data.elementary_sum.coefficient[index] = data.elementary[index];
        data.alternating_elementary_sum.coefficient[index] = degree % 2 == 0
            ? data.elementary[index]
            : field.negate(data.elementary[index]);
    }
    normalize(data.elementary_sum, degree_limit);
    normalize(data.alternating_elementary_sum, degree_limit);

    const auto pair_series = [degree_limit, largest_f_index, &field](
        const std::vector<std::uint32_t>& sequence
    ) {
        std::vector<Series> answer(
            static_cast<std::size_t>(largest_f_index + 1)
        );
        for (int shift = 0; shift <= largest_f_index; ++shift) {
            Series value;
            for (int lower = 0; 2 * lower + shift <= degree_limit; ++lower) {
                const int degree = 2 * lower + shift;
                value.coefficient[static_cast<std::size_t>(degree)] = field.add(
                    value.coefficient[static_cast<std::size_t>(degree)],
                    field.multiply(
                        sequence[static_cast<std::size_t>(lower)],
                        sequence[static_cast<std::size_t>(lower + shift)]
                    )
                );
            }
            normalize(value, degree_limit);
            answer[static_cast<std::size_t>(shift)] = value;
        }
        return answer;
    };
    data.elementary_pair = pair_series(data.elementary);
    data.complete_pair = pair_series(data.complete);
    return data;
}

Series pair_at(const std::vector<Series>& pairs, int index) {
    if (index < 0) index = -index;
    if (index >= static_cast<int>(pairs.size())) {
        throw std::runtime_error("bounded-Littlewood f-index overflow");
    }
    return pairs[static_cast<std::size_t>(index)];
}

std::vector<Series> build_b_matrix(
    const std::vector<Series>& pairs,
    int size,
    bool plus,
    const MontgomeryField& field,
    int degree_limit
) {
    std::vector<Series> matrix(static_cast<std::size_t>(size * size));
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            matrix[static_cast<std::size_t>(row * size + column)] = add_series(
                pair_at(pairs, row - column),
                pair_at(pairs, row + column + 1),
                field,
                degree_limit,
                !plus
            );
        }
    }
    return matrix;
}

std::vector<Series> build_c_matrix(
    const std::vector<Series>& pairs,
    int size,
    const MontgomeryField& field,
    int degree_limit
) {
    std::vector<Series> matrix(static_cast<std::size_t>(size * size));
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            matrix[static_cast<std::size_t>(row * size + column)] = add_series(
                pair_at(pairs, row - column),
                pair_at(pairs, row + column + 2),
                field,
                degree_limit,
                true
            );
        }
    }
    return matrix;
}

std::vector<Series> build_d_matrix(
    const std::vector<Series>& pairs,
    int size,
    const MontgomeryField& field,
    int degree_limit
) {
    std::vector<Series> matrix(static_cast<std::size_t>(size * size));
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            matrix[static_cast<std::size_t>(row * size + column)] = add_series(
                pair_at(pairs, row - column),
                pair_at(pairs, row + column),
                field,
                degree_limit,
                false
            );
        }
    }
    return matrix;
}

std::vector<std::vector<std::uint32_t>> gaussian_functional_weights(
    int maximum_moment,
    const MontgomeryField& field
) {
    const int node_count = maximum_moment + 1;
    std::vector<std::vector<std::uint32_t>> lagrange(
        static_cast<std::size_t>(node_count),
        std::vector<std::uint32_t>(static_cast<std::size_t>(node_count), 0U)
    );

    for (int node = 0; node < node_count; ++node) {
        std::vector<std::uint32_t> polynomial(1U, field.one());
        std::uint32_t denominator = field.one();
        for (int other = 0; other < node_count; ++other) {
            if (other == node) continue;
            std::vector<std::uint32_t> next(polynomial.size() + 1U, 0U);
            const std::uint32_t minus_other = field.negate(
                field.from_unsigned(static_cast<std::uint64_t>(other))
            );
            for (std::size_t degree = 0; degree < polynomial.size(); ++degree) {
                next[degree] = field.add(
                    next[degree],
                    field.multiply(polynomial[degree], minus_other)
                );
                next[degree + 1U] = field.add(
                    next[degree + 1U], polynomial[degree]
                );
            }
            polynomial = std::move(next);
            denominator = field.multiply(
                denominator,
                field.from_signed(static_cast<std::int64_t>(node - other))
            );
        }
        const std::uint32_t denominator_inverse = field.inverse(denominator);
        for (int degree = 0; degree < node_count; ++degree) {
            lagrange[static_cast<std::size_t>(node)]
                     [static_cast<std::size_t>(degree)] = field.multiply(
                polynomial[static_cast<std::size_t>(degree)],
                denominator_inverse
            );
        }
    }

    std::vector<std::uint32_t> factorial(
        static_cast<std::size_t>(maximum_moment + 1), field.one()
    );
    std::vector<std::uint32_t> odd_double_factorial(
        static_cast<std::size_t>(maximum_moment + 1), field.one()
    );
    for (int degree = 1; degree <= maximum_moment; ++degree) {
        factorial[static_cast<std::size_t>(degree)] = field.multiply(
            factorial[static_cast<std::size_t>(degree - 1)],
            field.from_unsigned(static_cast<std::uint64_t>(degree))
        );
        odd_double_factorial[static_cast<std::size_t>(degree)] = field.multiply(
            odd_double_factorial[static_cast<std::size_t>(degree - 1)],
            field.from_unsigned(static_cast<std::uint64_t>(2 * degree - 1))
        );
    }

    std::vector<std::vector<std::uint32_t>> weights(
        static_cast<std::size_t>(maximum_moment + 1),
        std::vector<std::uint32_t>(static_cast<std::size_t>(node_count), 0U)
    );
    for (int moment = 0; moment <= maximum_moment; ++moment) {
        for (int node = 0; node < node_count; ++node) {
            std::uint32_t value = 0U;
            for (int q_power = 0; q_power <= moment; ++q_power) {
                const int gaussian_half_degree = moment - q_power;
                const std::uint32_t target = field.multiply(
                    factorial[static_cast<std::size_t>(moment)],
                    odd_double_factorial[
                        static_cast<std::size_t>(gaussian_half_degree)
                    ]
                );
                value = field.add(
                    value,
                    field.multiply(
                        target,
                        lagrange[static_cast<std::size_t>(node)]
                                 [static_cast<std::size_t>(q_power)]
                    )
                );
            }
            weights[static_cast<std::size_t>(moment)]
                   [static_cast<std::size_t>(node)] = value;
        }
    }
    return weights;
}

struct MomentRow {
    char family;
    int rank;
    int moment_through;
    std::vector<BigInt> moments;
};

void validate_row_ledger() {
    std::size_t polynomial_count = 0U;
    std::size_t rational_cap_count = 0U;
    std::size_t directed_interval_count = 0U;
    int expected_b_rank = 2;
    int expected_c_rank = 2;
    int expected_d_rank = 4;
    int largest_moment = 0;

    for (const auto& cutoff : full_q3_bcd_remaining::row_cutoffs) {
        if ((cutoff.family == 'B' && cutoff.rank != expected_b_rank++)
            || (cutoff.family == 'C' && cutoff.rank != expected_c_rank++)
            || (cutoff.family == 'D' && cutoff.rank != expected_d_rank++)
            || (cutoff.family != 'B' && cutoff.family != 'C'
                && cutoff.family != 'D')) {
            throw std::runtime_error("bounded row-rank sequence mismatch");
        }
        const int stable_through = cutoff.family == 'D'
            ? cutoff.rank - 1
            : 2 * cutoff.rank + 1;
        if (cutoff.tail_onset % 2 == 0
            || cutoff.tail_onset <= stable_through
            || cutoff.moment_through < cutoff.tail_onset - 2) {
            throw std::runtime_error("invalid bounded row coverage ledger");
        }
        largest_moment = std::max(largest_moment, cutoff.moment_through);
        switch (cutoff.tail_method) {
            case full_q3_bcd_remaining::TailMethod::polynomial:
                ++polynomial_count;
                break;
            case full_q3_bcd_remaining::TailMethod::rational_cap:
                ++rational_cap_count;
                break;
            case full_q3_bcd_remaining::TailMethod::directed_interval:
                ++directed_interval_count;
                break;
        }
    }
    if (expected_b_rank != 22 || expected_c_rank != 29
        || expected_d_rank != 71
        || polynomial_count != full_q3_bcd_remaining::polynomial_rows
        || rational_cap_count != full_q3_bcd_remaining::rational_cap_rows
        || directed_interval_count
            != full_q3_bcd_remaining::directed_interval_rows
        || largest_moment != kCertifiedMaximum) {
        throw std::runtime_error("bounded row ledger scope mismatch");
    }
}

int group_dimension(char family, int rank) {
    if (family == 'B' || family == 'C') return rank * (2 * rank + 1);
    if (family == 'D') return rank * (2 * rank - 1);
    throw std::runtime_error("unknown classical family");
}

bool is_prime_32(std::uint32_t value) {
    if (value < 2U) return false;
    if (value % 2U == 0U) return value == 2U;
    for (std::uint32_t divisor = 3U;
         static_cast<std::uint64_t>(divisor) * divisor <= value;
         divisor += 2U) {
        if (value % divisor == 0U) return false;
    }
    return true;
}

std::vector<std::uint32_t> descending_primes(std::size_t count) {
    std::vector<std::uint32_t> primes;
    std::uint32_t candidate = 2147483647U;
    while (primes.size() < count) {
        if (is_prime_32(candidate)) primes.push_back(candidate);
        if (candidate <= 3U) throw std::runtime_error("prime search exhausted");
        candidate -= 2U;
    }
    return primes;
}

BigInt integer_power(int base, int exponent) {
    BigInt answer;
    mpz_ui_pow_ui(
        answer.get_mpz_t(),
        static_cast<unsigned long>(base),
        static_cast<unsigned long>(exponent)
    );
    return answer;
}

std::uint32_t normal_mod_inverse(std::uint32_t value, std::uint32_t prime) {
    std::uint64_t base = value;
    std::uint64_t exponent = static_cast<std::uint64_t>(prime) - 2U;
    std::uint64_t answer = 1U;
    while (exponent > 0U) {
        if ((exponent & 1U) != 0U) answer = (answer * base) % prime;
        exponent >>= 1U;
        if (exponent > 0U) base = (base * base) % prime;
    }
    return static_cast<std::uint32_t>(answer);
}

std::vector<std::vector<std::uint32_t>> moment_residues_for_prime(
    const std::vector<MomentRow>& rows,
    int maximum_moment,
    std::uint32_t prime,
    bool progress
) {
    const int degree_limit = 2 * maximum_moment;
    const MontgomeryField field(prime);
    const auto weights = gaussian_functional_weights(maximum_moment, field);
    std::vector<std::vector<std::uint32_t>> residues(
        rows.size(),
        std::vector<std::uint32_t>(
            static_cast<std::size_t>(maximum_moment + 1), 0U
        )
    );
    const std::uint32_t inverse_two = field.inverse(field.from_unsigned(2U));

    int b_maximum = 0;
    int c_maximum = 0;
    int d_maximum = 0;
    for (const MomentRow& row : rows) {
        if (row.family == 'B') b_maximum = std::max(b_maximum, row.rank);
        if (row.family == 'C') c_maximum = std::max(c_maximum, row.rank);
        if (row.family == 'D') d_maximum = std::max(d_maximum, row.rank);
    }
    const int largest_f_index = std::max({
        2 * b_maximum - 1,
        2 * c_maximum,
        2 * d_maximum - 2
    });

    for (int node_index = 0; node_index <= maximum_moment; ++node_index) {
        const std::uint32_t node = field.from_unsigned(
            static_cast<std::uint64_t>(node_index)
        );
        const NodeData data = build_node_data(
            node, largest_f_index, field, degree_limit
        );

        const auto b_minus = leading_principal_determinants(
            build_b_matrix(
                data.elementary_pair,
                b_maximum,
                false,
                field,
                degree_limit
            ),
            b_maximum,
            field,
            degree_limit
        );
        const auto b_plus = leading_principal_determinants(
            build_b_matrix(
                data.elementary_pair,
                b_maximum,
                true,
                field,
                degree_limit
            ),
            b_maximum,
            field,
            degree_limit
        );
        const auto c_determinants = leading_principal_determinants(
            build_c_matrix(
                data.complete_pair,
                c_maximum,
                field,
                degree_limit
            ),
            c_maximum,
            field,
            degree_limit
        );
        const auto d_determinants = leading_principal_determinants(
            build_d_matrix(
                data.elementary_pair,
                d_maximum,
                field,
                degree_limit
            ),
            d_maximum,
            field,
            degree_limit
        );

        for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
            const MomentRow& row = rows[row_index];
            Series generating;
            if (row.family == 'B') {
                const Series first = multiply_series(
                    data.elementary_sum,
                    b_minus[static_cast<std::size_t>(row.rank)],
                    field,
                    degree_limit
                );
                const Series second = multiply_series(
                    data.alternating_elementary_sum,
                    b_plus[static_cast<std::size_t>(row.rank)],
                    field,
                    degree_limit
                );
                generating = scale_series(
                    add_series(first, second, field, degree_limit),
                    inverse_two,
                    field,
                    degree_limit
                );
            } else if (row.family == 'C') {
                generating = c_determinants[
                    static_cast<std::size_t>(row.rank)
                ];
            } else if (row.family == 'D') {
                generating = scale_series(
                    d_determinants[static_cast<std::size_t>(row.rank)],
                    inverse_two,
                    field,
                    degree_limit
                );
            } else {
                throw std::runtime_error("unknown bounded-Littlewood family");
            }

            const int through = std::min(row.moment_through, maximum_moment);
            for (int moment = 0; moment <= through; ++moment) {
                const std::uint32_t coefficient = generating.coefficient[
                    static_cast<std::size_t>(2 * moment)
                ];
                residues[row_index][static_cast<std::size_t>(moment)] =
                    field.add(
                        residues[row_index][static_cast<std::size_t>(moment)],
                        field.multiply(
                            weights[static_cast<std::size_t>(moment)]
                                   [static_cast<std::size_t>(node_index)],
                            coefficient
                        )
                    );
            }
        }
        if (progress) {
            std::cout << "FULL_Q3_BOUNDED_LITTLEWOOD prime=" << prime
                      << " completed_node=" << node_index << '/'
                      << maximum_moment << std::endl;
        }
    }

    for (auto& row : residues) {
        for (std::uint32_t& value : row) value = field.to_unsigned(value);
    }
    return residues;
}

std::vector<MomentRow> reconstruct_bounded_littlewood_moments(
    int maximum_moment,
    bool progress
) {
    std::vector<MomentRow> rows;
    BigInt largest_bound = 0;
    for (const auto& cutoff : full_q3_bcd_remaining::row_cutoffs) {
        const int through = std::min(cutoff.moment_through, maximum_moment);
        rows.push_back(MomentRow{
            cutoff.family,
            cutoff.rank,
            through,
            std::vector<BigInt>(static_cast<std::size_t>(through + 1), 0)
        });
        largest_bound = std::max(
            largest_bound,
            integer_power(group_dimension(cutoff.family, cutoff.rank), through)
        );
    }

    std::vector<std::uint32_t> primes;
    BigInt planned_modulus_product = 1;
    for (std::uint32_t prime : descending_primes(32U)) {
        primes.push_back(prime);
        planned_modulus_product *= prime;
        if (planned_modulus_product > largest_bound) break;
    }
    if (planned_modulus_product <= largest_bound) {
        throw std::runtime_error("CRT prime table does not dominate character bound");
    }

    using ResidueTable = std::vector<std::vector<std::uint32_t>>;
    std::vector<ResidueTable> residues_by_prime(primes.size());
    std::vector<std::string> errors(primes.size());
    const int prime_count = static_cast<int>(primes.size());

    // Prime-field evaluations are independent.  Parallelizing at this level
    // keeps every field computation deterministic and leaves the GMP CRT merge
    // in the fixed descending-prime order below.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int prime_index = 0; prime_index < prime_count; ++prime_index) {
        const std::size_t index = static_cast<std::size_t>(prime_index);
        try {
            residues_by_prime[index] = moment_residues_for_prime(
                rows, maximum_moment, primes[index], false
            );
            if (progress) {
#ifdef _OPENMP
#pragma omp critical(full_q3_bounded_littlewood_progress)
#endif
                std::cout << "FULL_Q3_BOUNDED_LITTLEWOOD evaluated_prime="
                          << primes[index] << " prime_index=" << prime_index
                          << '/' << prime_count << std::endl;
            }
        } catch (const std::exception& error) {
            errors[index] = error.what();
        } catch (...) {
            errors[index] = "unknown prime-field evaluation failure";
        }
    }
    for (std::size_t index = 0; index < errors.size(); ++index) {
        if (!errors[index].empty()) {
            throw std::runtime_error(
                std::string("prime ") + std::to_string(primes[index])
                + " failed: " + errors[index]
            );
        }
    }

    BigInt modulus_product = 1;
    std::size_t primes_consumed = 0U;
    for (std::size_t prime_index = 0;
         prime_index < primes.size();
         ++prime_index) {
        const std::uint32_t prime = primes[prime_index];
        const ResidueTable& residues = residues_by_prime[prime_index];
        const std::uint32_t modulus_residue = static_cast<std::uint32_t>(
            mpz_fdiv_ui(modulus_product.get_mpz_t(), prime)
        );
        const std::uint32_t inverse_modulus = normal_mod_inverse(
            modulus_residue, prime
        );
        for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
            MomentRow& row = rows[row_index];
            for (int moment = 0; moment <= row.moment_through; ++moment) {
                BigInt& value = row.moments[static_cast<std::size_t>(moment)];
                const std::uint32_t current = static_cast<std::uint32_t>(
                    mpz_fdiv_ui(value.get_mpz_t(), prime)
                );
                const std::uint32_t residue =
                    residues[row_index][static_cast<std::size_t>(moment)];
                const std::uint32_t difference = residue >= current
                    ? residue - current
                    : static_cast<std::uint32_t>(
                          static_cast<std::uint64_t>(residue) + prime - current
                      );
                const std::uint32_t multiplier = static_cast<std::uint32_t>(
                    (static_cast<std::uint64_t>(difference) * inverse_modulus)
                    % prime
                );
                value += modulus_product * multiplier;
            }
        }
        modulus_product *= prime;
        ++primes_consumed;
        std::cout << "FULL_Q3_BOUNDED_LITTLEWOOD completed_prime=" << prime
                  << " primes_consumed=" << primes_consumed
                  << " modulus_bits=" << mpz_sizeinbase(
                         modulus_product.get_mpz_t(), 2
                     ) << std::endl;
    }
    if (modulus_product <= largest_bound) {
        throw std::runtime_error("CRT modulus does not dominate character bound");
    }

    for (const MomentRow& row : rows) {
        const int dimension = group_dimension(row.family, row.rank);
        for (int moment = 0; moment <= row.moment_through; ++moment) {
            const BigInt& value = row.moments[static_cast<std::size_t>(moment)];
            if (value < 0 || value > integer_power(dimension, moment)) {
                throw std::runtime_error(
                    "bounded-Littlewood moment left its character bound"
                );
            }
        }
    }
    return rows;
}

mpq_class rational(unsigned long numerator, unsigned long denominator) {
    if (denominator == 0UL) {
        throw std::runtime_error("zero exact-rational denominator");
    }
    mpq_class answer(numerator);
    answer /= denominator;
    answer.canonicalize();
    return answer;
}

mpq_class rational_power(mpq_class base, int exponent) {
    if (exponent < 0) {
        throw std::runtime_error("negative exact-rational exponent");
    }
    mpq_class answer(1);
    while (exponent > 0) {
        if (exponent % 2 != 0) answer *= base;
        exponent /= 2;
        if (exponent > 0) base *= base;
    }
    return answer;
}

BigInt factorial_exact(int value) {
    if (value < 0) throw std::runtime_error("negative factorial argument");
    BigInt answer = 1;
    for (int factor = 2; factor <= value; ++factor) answer *= factor;
    return answer;
}

BigInt power_of_two_exact(int exponent) {
    if (exponent < 0) {
        throw std::runtime_error("negative power-of-two exponent");
    }
    BigInt answer;
    mpz_ui_pow_ui(
        answer.get_mpz_t(),
        2UL,
        static_cast<unsigned long>(exponent)
    );
    return answer;
}

BigInt degree_factor_exact(std::initializer_list<int> degrees) {
    BigInt answer = 1;
    for (int degree : degrees) answer *= factorial_exact(degree);
    return answer;
}

int negative_endpoint(char family, int rank) {
    if (family == 'B' || family == 'C') return rank;
    if (family == 'D') return rank % 2 == 0 ? rank : rank - 2;
    throw std::runtime_error("unknown classical family in negative endpoint");
}

mpq_class exact_tail_margin(
    const mpq_class& probability,
    const mpq_class& push,
    int endpoint,
    int onset
) {
    if (probability <= 0 || probability > 1) {
        throw std::runtime_error("invalid exact tail probability");
    }
    if (push <= 2 * endpoint) {
        throw std::runtime_error("exact tail push does not exceed negative cap");
    }
    return probability * rational_power(push, onset)
        - rational_power(mpq_class(2 * endpoint), onset);
}

struct PolynomialTailSpec {
    char family;
    int rank;
    int power;
    unsigned long threshold_numerator;
    unsigned long threshold_denominator;
    unsigned long central_numerator;
    unsigned long central_denominator;
    int onset;
};

constexpr std::array<PolynomialTailSpec, 6> polynomial_tail_specs{{
    {'B', 3, 8, 227, 20, 23, 20, 37},
    {'B', 4, 10, 331, 20, 5, 4, 41},
    {'B', 5, 10, 75, 4, 6, 5, 55},
    {'C', 3, 8, 57, 5, 6, 5, 37},
    {'D', 4, 9, 57, 4, 23, 20, 47},
    {'D', 5, 7, 131, 10, 6, 5, 33},
}};

const PolynomialTailSpec* find_polynomial_tail(char family, int rank) {
    for (const PolynomialTailSpec& spec : polynomial_tail_specs) {
        if (spec.family == family && spec.rank == rank) return &spec;
    }
    return nullptr;
}

void verify_polynomial_tail(
    const MomentRow& row,
    const full_q3_bcd_remaining::RowCutoff& cutoff
) {
    const PolynomialTailSpec* spec = find_polynomial_tail(row.family, row.rank);
    if (spec == nullptr || spec->onset != cutoff.tail_onset) {
        throw std::runtime_error("polynomial-tail ledger mismatch");
    }
    const int k = spec->power;
    const int required_moment = 4 * k + 2;
    if (required_moment > row.moment_through
        || required_moment >= static_cast<int>(row.moments.size())) {
        throw std::runtime_error("polynomial tail lacks required exact moments");
    }
    const mpq_class threshold = rational(
        spec->threshold_numerator, spec->threshold_denominator
    );
    const mpq_class central = rational(
        spec->central_numerator, spec->central_denominator
    );
    const mpq_class a_value =
        mpq_class(row.moments[static_cast<std::size_t>(2 * k + 1)])
        - threshold
            * mpq_class(row.moments[static_cast<std::size_t>(2 * k)]);
    const mpq_class b_value =
        mpq_class(row.moments[static_cast<std::size_t>(4 * k + 2)])
        - 2 * threshold
            * mpq_class(row.moments[static_cast<std::size_t>(4 * k + 1)])
        + threshold * threshold
            * mpq_class(row.moments[static_cast<std::size_t>(4 * k)]);
    if (a_value <= 0 || b_value <= 0 || central <= 1) {
        throw std::runtime_error("invalid polynomial one-sided-tail inputs");
    }
    const mpq_class probability =
        (1 - 1 / (central * central)) * a_value * a_value / b_value;
    const mpq_class push = threshold - central;
    const mpq_class margin = exact_tail_margin(
        probability,
        push,
        negative_endpoint(row.family, row.rank),
        cutoff.tail_onset
    );
    if (margin <= 0) {
        throw std::runtime_error("polynomial-tail comparison is not positive");
    }
    std::cout << "FULL_Q3_BOUNDED_EXACT_TAIL row=" << row.family << '_'
              << row.rank << " method=polynomial onset=" << cutoff.tail_onset
              << " power=" << k
              << " margin_numerator=" << margin.get_num()
              << " margin_denominator=" << margin.get_den() << '\n';
}

void verify_rational_cap_tail(
    const MomentRow& row,
    const full_q3_bcd_remaining::RowCutoff& cutoff
) {
    const char family = row.family;
    const int rank = row.rank;
    int dimension = 0;
    int kappa = 0;
    int endpoint = 0;
    mpq_class radius;
    mpq_class central;
    mpq_class push;
    BigInt degree_factor;
    BigInt weyl_order;
    mpq_class cap_probability;
    int expected_onset = 0;

    if ((family == 'B' || family == 'C') && rank == 2) {
        expected_onset = 47;
        dimension = 10;
        endpoint = 2;
        radius = rational(1, 1);
        central = rational(6, 5);
        push = rational(39, 5);
        const mpq_class width = rational(1, 4);
        const mpq_class cap_integral =
            rational_power(width, 10)
            * (rational(1, 21) - rational(2, 25) + rational(1, 21));
        cap_probability =
            rational(1, 512) * rational_power(width, 4) * cap_integral;
    } else if (family == 'B' && rank == 6) {
        expected_onset = 51;
        dimension = 78;
        kappa = 11;
        endpoint = 6;
        radius = rational(1053, 50);
        central = rational(3, 2);
        push = rational(1386, 25);
        degree_factor = degree_factor_exact({2, 4, 6, 8, 10, 12});
        weyl_order = power_of_two_exact(6) * factorial_exact(6);
        const int shape = dimension / 2;
        const mpq_class density = 1 / (
            rational_power(rational(44, 7), rank / 2)
            * factorial_exact(shape)
        );
        cap_probability =
            rational(1, 1) / weyl_order
            * rational_power(1 - radius / 24, 2)
            * degree_factor * density
            * rational_power(radius / (2 * kappa), shape)
            / power_of_two_exact(rank);
        if (radius >= 24) {
            throw std::runtime_error("B6 product-sine cap left its domain");
        }
    } else if (family == 'D' && rank == 6) {
        expected_onset = 49;
        dimension = 66;
        kappa = 10;
        endpoint = 6;
        radius = rational(371, 20);
        central = rational(7, 5);
        push = rational(921, 20);
        degree_factor = degree_factor_exact({2, 4, 6, 8, 10, 6});
        weyl_order = power_of_two_exact(5) * factorial_exact(6);
        const int shape = dimension / 2;
        const mpq_class density = 1 / (
            rational_power(rational(44, 7), rank / 2)
            * factorial_exact(shape)
        );
        cap_probability =
            rational(1, 1) / weyl_order
            * rational_power(1 - radius / 24, 2)
            * degree_factor * density
            * rational_power(radius / (2 * kappa), shape);
        if (radius >= 24) {
            throw std::runtime_error("D6 product-sine cap left its domain");
        }
    } else if (family == 'D' && rank == 7) {
        expected_onset = 53;
        dimension = 91;
        kappa = 12;
        endpoint = 5;
        radius = rational(18, 1);
        central = rational(8, 5);
        push = rational(357, 5);
        degree_factor = degree_factor_exact({2, 4, 6, 8, 10, 12, 7});
        weyl_order = power_of_two_exact(6) * factorial_exact(7);
        const int half_index = (dimension + 1) / 2;
        const mpq_class density =
            rational_power(rational(4, 1), half_index)
            * factorial_exact(half_index)
            / (
                12 * rational_power(rational(22, 7), 4)
                * factorial_exact(2 * half_index)
            );
        cap_probability =
            rational(1, 1) / weyl_order
            * rational_power(1 - radius / 24, 2)
            * degree_factor * density
            * rational_power(radius / (2 * kappa), (dimension - 1) / 2)
            * rational(433, 500);
        if (radius >= 24
            || rational_power(rational(433, 500), 2)
                >= radius / (2 * kappa)
            || rational_power(rational(3, 2), 2) <= 2) {
            throw std::runtime_error("D7 rational cap bound is invalid");
        }
    } else if (family == 'D' && rank == 9) {
        expected_onset = 61;
        dimension = 153;
        kappa = 16;
        endpoint = 7;
        radius = rational(763, 10);
        central = rational(8, 5);
        push = rational(751, 10);
        degree_factor = degree_factor_exact({2, 4, 6, 8, 10, 12, 14, 16, 9});
        weyl_order = power_of_two_exact(8) * factorial_exact(9);
        const int half_index = (dimension + 1) / 2;
        const mpq_class density =
            rational(7, 5)
            * rational_power(rational(4, 1), half_index)
            * factorial_exact(half_index)
            / (
                rational_power(rational(44, 7), 5)
                * factorial_exact(2 * half_index)
            );
        cap_probability =
            rational(1, 1) / weyl_order
            * rational_power(1 - 2 * radius / (24 * kappa), kappa)
            * degree_factor * density
            * rational_power(radius / (2 * kappa), (dimension - 1) / 2)
            * rational(193, 125);
        if (radius >= 12 * kappa
            || rational_power(rational(193, 125), 2)
                >= radius / (2 * kappa)
            || rational_power(rational(10, 7), 2) <= 2) {
            throw std::runtime_error("D9 rational cap bound is invalid");
        }
    } else {
        throw std::runtime_error("unknown exact rational-cap row");
    }

    if (push != dimension - radius - central
        || endpoint != negative_endpoint(family, rank)
        || cutoff.tail_onset != expected_onset) {
        throw std::runtime_error("rational-cap geometry ledger mismatch");
    }
    if (kappa > 0 && radius >= 9 * kappa) {
        throw std::runtime_error("rational cap leaves the eigenangle chart");
    }
    const mpq_class central_probability = 1 - 1 / (central * central);
    const mpq_class probability = cap_probability * central_probability;
    const mpq_class margin = exact_tail_margin(
        probability, push, endpoint, cutoff.tail_onset
    );
    if (margin <= 0) {
        throw std::runtime_error("rational-cap tail comparison is not positive");
    }
    std::cout << "FULL_Q3_BOUNDED_EXACT_TAIL row=" << family << '_' << rank
              << " method=rational_cap onset=" << cutoff.tail_onset
              << " margin_numerator=" << margin.get_num()
              << " margin_denominator=" << margin.get_den() << '\n';
}

void verify_exact_tail_ledgers(const std::vector<MomentRow>& rows) {
    std::size_t polynomial_checked = 0U;
    std::size_t rational_cap_checked = 0U;
    for (const MomentRow& row : rows) {
        const auto* cutoff = full_q3_bcd_remaining::find_row(
            row.family, row.rank
        );
        if (cutoff == nullptr) {
            throw std::runtime_error("missing exact-tail cutoff row");
        }
        if (cutoff->tail_method
            == full_q3_bcd_remaining::TailMethod::polynomial) {
            verify_polynomial_tail(row, *cutoff);
            ++polynomial_checked;
        } else if (cutoff->tail_method
                   == full_q3_bcd_remaining::TailMethod::rational_cap) {
            verify_rational_cap_tail(row, *cutoff);
            ++rational_cap_checked;
        }
    }
    if (polynomial_checked != full_q3_bcd_remaining::polynomial_rows
        || rational_cap_checked
            != full_q3_bcd_remaining::rational_cap_rows) {
        throw std::runtime_error("bounded exact-tail scope mismatch");
    }
    std::cout << "FULL_Q3_BOUNDED_EXACT_TAILS polynomial_rows="
              << polynomial_checked << " rational_cap_rows="
              << rational_cap_checked << '\n';
}

std::vector<BigInt> stable_moments(int maximum) {
    std::vector<BigInt> moments(static_cast<std::size_t>(maximum + 1), 0);
    moments[0] = 1;
    if (maximum >= 1) moments[1] = 0;
    for (int n = 1; n < maximum; ++n) {
        moments[static_cast<std::size_t>(n + 1)] =
            n * moments[static_cast<std::size_t>(n)]
            + n * moments[static_cast<std::size_t>(n - 1)];
        if (n >= 2) {
            moments[static_cast<std::size_t>(n + 1)] -=
                (n * (n - 1) / 2)
                * moments[static_cast<std::size_t>(n - 2)];
        }
    }
    return moments;
}

BigInt binomial(int top, int bottom) {
    if (top < 0 || bottom < 0 || bottom > top) return 0;
    BigInt answer;
    mpz_bin_uiui(
        answer.get_mpz_t(),
        static_cast<unsigned long>(top),
        static_cast<unsigned long>(bottom)
    );
    return answer;
}

BigInt hierarchy_value(
    const std::vector<BigInt>& moments,
    int a,
    int n
) {
    BigInt answer = 0;
    for (int j = 0; j <= 2 * a; ++j) {
        for (int k = 0; k <= n; ++k) {
            const BigInt term = binomial(2 * a, j) * binomial(n, k)
                * moments[static_cast<std::size_t>(2 * a - j + k)]
                * moments[static_cast<std::size_t>(j + n - k)];
            if (j % 2 == 0) answer += term;
            else answer -= term;
        }
    }
    return answer;
}

int first_odd_above(int value) {
    int answer = value + 1;
    if (answer % 2 == 0) ++answer;
    return answer;
}

void verify_moments_and_prefix(
    const std::vector<MomentRow>& rows,
    int maximum_moment
) {
    const auto stable = stable_moments(maximum_moment);
    const MomentRow* b2 = nullptr;
    const MomentRow* c2 = nullptr;
    std::size_t checked_pairs = 0U;
    bool have_minimum = false;
    BigInt minimum = 0;
    char minimum_family = 0;
    int minimum_rank = 0;
    int minimum_a = 0;
    int minimum_n = 0;
    std::size_t bd_checked_pairs = 0U;
    bool have_bd_minimum = false;
    BigInt bd_minimum = 0;
    char bd_minimum_family = 0;
    int bd_minimum_rank = 0;
    int bd_minimum_a = 0;
    int bd_minimum_n = 0;

    for (const MomentRow& row : rows) {
        const int stable_through = row.family == 'D'
            ? row.rank - 1
            : 2 * row.rank + 1;
        for (int degree = 0; degree <= row.moment_through; ++degree) {
            const BigInt& value = row.moments[static_cast<std::size_t>(degree)];
            if (degree <= stable_through
                && value != stable[static_cast<std::size_t>(degree)]) {
                throw std::runtime_error("bounded-Littlewood stable-prefix mismatch");
            }
            std::cout << "FULL_Q3_BOUNDED_MOMENT row=" << row.family << '_'
                      << row.rank << " degree=" << degree
                      << " value=" << value << '\n';
        }
        if (row.family == 'B' && row.rank == 2) b2 = &row;
        if (row.family == 'C' && row.rank == 2) c2 = &row;

        const auto* cutoff = full_q3_bcd_remaining::find_row(
            row.family, row.rank
        );
        if (cutoff == nullptr) throw std::runtime_error("missing cutoff row");
        std::size_t row_checked = 0U;
        for (int total = first_odd_above(stable_through);
             total <= std::min(maximum_moment, cutoff->tail_onset - 2);
             total += 2) {
            for (int a = 2; 2 * a < total; ++a) {
                const int n = total - 2 * a;
                const BigInt value = hierarchy_value(row.moments, a, n);
                if (value <= 0) {
                    throw std::runtime_error(
                        std::string("nonpositive bounded-Littlewood hierarchy at ")
                        + row.family + '_' + std::to_string(row.rank)
                        + " a=" + std::to_string(a)
                        + " n=" + std::to_string(n)
                    );
                }
                if (!have_minimum || value < minimum) {
                    minimum = value;
                    minimum_family = row.family;
                    minimum_rank = row.rank;
                    minimum_a = a;
                    minimum_n = n;
                    have_minimum = true;
                }
                ++checked_pairs;
                ++row_checked;
                const bool in_bd_overlap =
                    (row.family == 'B' && 18 <= row.rank && row.rank <= 21)
                    || (row.family == 'D' && 31 <= row.rank && row.rank <= 70);
                if (in_bd_overlap) {
                    ++bd_checked_pairs;
                    if (!have_bd_minimum || value < bd_minimum) {
                        bd_minimum = value;
                        bd_minimum_family = row.family;
                        bd_minimum_rank = row.rank;
                        bd_minimum_a = a;
                        bd_minimum_n = n;
                        have_bd_minimum = true;
                    }
                }
            }
        }
        std::cout << "FULL_Q3_BOUNDED_ROW row=" << row.family << '_'
                  << row.rank << " stable_through=" << stable_through
                  << " checked_through="
                  << std::min(maximum_moment, cutoff->tail_onset - 2)
                  << " moment_through=" << row.moment_through
                  << " tail_onset=" << cutoff->tail_onset
                  << " pairs=" << row_checked << '\n';
    }
    if (b2 == nullptr || c2 == nullptr) {
        throw std::runtime_error("missing B2/C2 moment rows");
    }
    const int common = std::min(b2->moment_through, c2->moment_through);
    for (int degree = 0; degree <= common; ++degree) {
        if (b2->moments[static_cast<std::size_t>(degree)]
            != c2->moments[static_cast<std::size_t>(degree)]) {
            throw std::runtime_error("B2/C2 bounded-Littlewood mismatch");
        }
    }
    std::cout << "FULL_Q3_BOUNDED_LITTLEWOOD pairs=" << checked_pairs;
    if (have_minimum) {
        std::cout << " minimum=" << minimum_family << '_' << minimum_rank
                  << " a=" << minimum_a << " n=" << minimum_n
                  << " value=" << minimum;
    }
    std::cout << '\n';
    if (maximum_moment == kCertifiedMaximum) {
        const BigInt expected_bd_minimum(
            "1272988116891284109857142182380802", 10
        );
        if (bd_checked_pairs != 4869U || !have_bd_minimum
            || bd_minimum != expected_bd_minimum
            || bd_minimum_family != 'D' || bd_minimum_rank != 31
            || bd_minimum_a != 8 || bd_minimum_n != 15) {
            throw std::runtime_error("bounded-Littlewood B/D overlap mismatch");
        }
        std::cout << "FULL_Q3_BOUNDED_BD_OVERLAP pairs=" << bd_checked_pairs
                  << " minimum=" << bd_minimum_family << '_'
                  << bd_minimum_rank << " a=" << bd_minimum_a
                  << " n=" << bd_minimum_n << " value=" << bd_minimum
                  << '\n';
    }
    if (maximum_moment == kCertifiedMaximum
        && checked_pairs != full_q3_bcd_remaining::full_residual_pairs) {
        throw std::runtime_error("bounded-Littlewood hierarchy scope mismatch");
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        int maximum_moment = kCertifiedMaximum;
        bool progress = false;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--maximum" && index + 1 < argc) {
                maximum_moment = std::stoi(argv[++index]);
            } else if (argument == "--progress") {
                progress = true;
            } else {
                std::cerr << "usage: " << argv[0]
                          << " [--maximum N] [--progress]\n";
                return EXIT_FAILURE;
            }
        }
        if (maximum_moment < 2 || maximum_moment > kCertifiedMaximum) {
            throw std::runtime_error("invalid bounded-Littlewood maximum");
        }
        int parallel_threads = 1;
#ifdef _OPENMP
        parallel_threads = omp_get_max_threads();
#endif
        std::cout << "FULL_Q3_BOUNDED_LITTLEWOOD maximum=" << maximum_moment
                  << " parallel_threads=" << parallel_threads << '\n';
        validate_row_ledger();
        const auto rows = reconstruct_bounded_littlewood_moments(
            maximum_moment, progress
        );
        if (maximum_moment == kCertifiedMaximum) {
            verify_exact_tail_ledgers(rows);
        }
        verify_moments_and_prefix(rows, maximum_moment);
        if (maximum_moment == kCertifiedMaximum) {
            std::cout << "FULL_Q3_BCD_BOUNDED_LITTLEWOOD VERIFICATION: ALL PASS\n";
        } else {
            std::cout << "FULL_Q3_BCD_BOUNDED_LITTLEWOOD DIAGNOSTIC PREFIX: PASS\n";
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "FULL_Q3_BCD_BOUNDED_LITTLEWOOD VERIFICATION: FAIL: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
