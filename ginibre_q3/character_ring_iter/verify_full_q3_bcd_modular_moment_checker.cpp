#include "full_q3_bcd_remaining_data.hpp"

#include <gmpxx.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using BigInt = mpz_class;

// This independent control intentionally freezes the original 70-row
// degree-59 subledger.  The promotion verifier now also covers the separate
// B18--B21/D31--D70 box, whose reverse-Pieri transcript is checked by its own
// optional control.
constexpr int kMaximumMoment = 59;
constexpr int kSeriesDegree = 2 * kMaximumMoment;

class PrimeField {
public:
    explicit PrimeField(std::uint32_t modulus) : modulus_(modulus) {
        if (modulus_ < 127U || modulus_ >= (std::uint32_t{1} << 31U)) {
            throw std::runtime_error("independent checker received invalid prime");
        }
    }

    [[nodiscard]] std::uint32_t modulus() const { return modulus_; }
    [[nodiscard]] std::uint32_t one() const { return 1U; }

    [[nodiscard]] std::uint32_t add(
        std::uint32_t left,
        std::uint32_t right
    ) const {
        const std::uint64_t sum = static_cast<std::uint64_t>(left) + right;
        return static_cast<std::uint32_t>(sum % modulus_);
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
        return static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(left) * right) % modulus_
        );
    }

    [[nodiscard]] std::uint32_t from_signed(std::int64_t value) const {
        if (value >= 0) {
            return static_cast<std::uint32_t>(
                static_cast<std::uint64_t>(value) % modulus_
            );
        }
        const std::uint64_t magnitude =
            static_cast<std::uint64_t>(-(value + 1)) + 1U;
        return negate(static_cast<std::uint32_t>(magnitude % modulus_));
    }

    [[nodiscard]] std::uint32_t power(
        std::uint32_t base,
        std::uint64_t exponent
    ) const {
        std::uint32_t answer = one();
        while (exponent > 0U) {
            if ((exponent & 1U) != 0U) answer = multiply(answer, base);
            exponent >>= 1U;
            if (exponent > 0U) base = multiply(base, base);
        }
        return answer;
    }

    [[nodiscard]] std::uint32_t inverse(std::uint32_t value) const {
        if (value == 0U) throw std::runtime_error("division by zero modulo prime");
        return power(value, static_cast<std::uint64_t>(modulus_) - 2U);
    }

private:
    std::uint32_t modulus_;
};

struct Polynomial {
    std::array<std::uint32_t, kSeriesDegree + 1> coefficient{};
    int lower = kSeriesDegree + 1;
    int upper = -1;

    [[nodiscard]] bool zero() const { return upper < lower; }
};

void trim(Polynomial& value, int limit) {
    value.lower = 0;
    while (value.lower <= limit
           && value.coefficient[static_cast<std::size_t>(value.lower)] == 0U) {
        ++value.lower;
    }
    if (value.lower > limit) {
        value.lower = kSeriesDegree + 1;
        value.upper = -1;
        return;
    }
    value.upper = limit;
    while (value.upper >= value.lower
           && value.coefficient[static_cast<std::size_t>(value.upper)] == 0U) {
        --value.upper;
    }
}

Polynomial constant_polynomial(std::uint32_t value, int limit) {
    Polynomial answer;
    answer.coefficient[0] = value;
    trim(answer, limit);
    return answer;
}

Polynomial sum_polynomial(
    const Polynomial& left,
    const Polynomial& right,
    const PrimeField& field,
    int limit,
    bool subtract_right = false
) {
    Polynomial answer;
    for (int degree = 0; degree <= limit; ++degree) {
        const std::size_t index = static_cast<std::size_t>(degree);
        answer.coefficient[index] = subtract_right
            ? field.subtract(left.coefficient[index], right.coefficient[index])
            : field.add(left.coefficient[index], right.coefficient[index]);
    }
    trim(answer, limit);
    return answer;
}

void subtract_polynomial_in_place(
    Polynomial& left,
    const Polynomial& right,
    const PrimeField& field,
    int limit
) {
    if (right.zero()) return;
    for (int degree = right.lower; degree <= right.upper; ++degree) {
        const std::size_t index = static_cast<std::size_t>(degree);
        left.coefficient[index] = field.subtract(
            left.coefficient[index], right.coefficient[index]
        );
    }
    trim(left, limit);
}

Polynomial scale_polynomial(
    const Polynomial& input,
    std::uint32_t scalar,
    const PrimeField& field,
    int limit
) {
    Polynomial answer;
    if (input.zero() || scalar == 0U) return answer;
    for (int degree = input.lower; degree <= input.upper; ++degree) {
        answer.coefficient[static_cast<std::size_t>(degree)] = field.multiply(
            input.coefficient[static_cast<std::size_t>(degree)], scalar
        );
    }
    trim(answer, limit);
    return answer;
}

Polynomial product_polynomial(
    const Polynomial& left,
    const Polynomial& right,
    const PrimeField& field,
    int limit
) {
    Polynomial answer;
    if (left.zero() || right.zero() || left.lower + right.lower > limit) {
        return answer;
    }
    const int left_upper = std::min(left.upper, limit - right.lower);
    for (int left_degree = left.lower;
         left_degree <= left_upper;
         ++left_degree) {
        const std::uint32_t left_value =
            left.coefficient[static_cast<std::size_t>(left_degree)];
        if (left_value == 0U) continue;
        const int right_upper = std::min(right.upper, limit - left_degree);
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
    trim(answer, limit);
    return answer;
}

Polynomial reciprocal_polynomial(
    const Polynomial& input,
    const PrimeField& field,
    int limit
) {
    if (input.coefficient[0] == 0U) {
        throw std::runtime_error("nonunit formal-series pivot");
    }
    Polynomial answer;
    const std::uint32_t inverse_constant = field.inverse(input.coefficient[0]);
    answer.coefficient[0] = inverse_constant;
    for (int degree = 1; degree <= limit; ++degree) {
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
    trim(answer, limit);
    return answer;
}

// This implementation performs an ordinary LU factorization over
// F_p[[t]]/(t^(limit+1)).  It is deliberately separate from the production
// Montgomery implementation and uses a dynamic polynomial representation.
std::vector<Polynomial> principal_determinants(
    std::vector<Polynomial> matrix,
    int size,
    const PrimeField& field,
    int limit
) {
    if (size < 1
        || matrix.size() != static_cast<std::size_t>(size * size)) {
        throw std::runtime_error("invalid independent determinant matrix");
    }
    const auto at = [size, &matrix](int row, int column) -> Polynomial& {
        return matrix[static_cast<std::size_t>(row * size + column)];
    };
    std::vector<Polynomial> output(static_cast<std::size_t>(size + 1));
    output[0] = constant_polynomial(field.one(), limit);
    Polynomial determinant = output[0];
    for (int pivot_index = 0; pivot_index < size; ++pivot_index) {
        const Polynomial pivot = at(pivot_index, pivot_index);
        if (pivot.coefficient[0] == 0U) {
            throw std::runtime_error("zero constant pivot in independent LU");
        }
        determinant = product_polynomial(determinant, pivot, field, limit);
        output[static_cast<std::size_t>(pivot_index + 1)] = determinant;
        const Polynomial reciprocal = reciprocal_polynomial(pivot, field, limit);
        for (int row = pivot_index + 1; row < size; ++row) {
            if (at(row, pivot_index).zero()) continue;
            const Polynomial factor = product_polynomial(
                at(row, pivot_index), reciprocal, field, limit
            );
            at(row, pivot_index) = Polynomial{};
            for (int column = pivot_index + 1; column < size; ++column) {
                subtract_polynomial_in_place(
                    at(row, column),
                    product_polynomial(
                        factor, at(pivot_index, column), field, limit
                    ),
                    field,
                    limit
                );
            }
        }
    }
    return output;
}

struct NodeData {
    std::vector<std::uint32_t> elementary;
    std::vector<std::uint32_t> complete;
    std::vector<Polynomial> elementary_pairs;
    std::vector<Polynomial> complete_pairs;
    Polynomial elementary_sum;
    Polynomial alternating_elementary_sum;
};

NodeData make_node_data(
    std::uint32_t node,
    int largest_pair_index,
    const PrimeField& field,
    int limit
) {
    NodeData data;
    data.elementary.assign(static_cast<std::size_t>(limit + 1), 0U);
    data.complete.assign(static_cast<std::size_t>(limit + 1), 0U);
    data.elementary[0] = field.one();
    data.complete[0] = field.one();
    for (int degree = 1; degree <= limit; ++degree) {
        std::uint32_t elementary =
            data.elementary[static_cast<std::size_t>(degree - 1)];
        std::uint32_t complete =
            data.complete[static_cast<std::size_t>(degree - 1)];
        if (degree >= 2) {
            elementary = field.subtract(
                elementary,
                field.multiply(
                    node,
                    data.elementary[static_cast<std::size_t>(degree - 2)]
                )
            );
            complete = field.add(
                complete,
                field.multiply(
                    node,
                    data.complete[static_cast<std::size_t>(degree - 2)]
                )
            );
        }
        const std::uint32_t inverse_degree = field.inverse(
            static_cast<std::uint32_t>(degree)
        );
        data.elementary[static_cast<std::size_t>(degree)] =
            field.multiply(elementary, inverse_degree);
        data.complete[static_cast<std::size_t>(degree)] =
            field.multiply(complete, inverse_degree);
    }

    for (int degree = 0; degree <= limit; ++degree) {
        const std::size_t index = static_cast<std::size_t>(degree);
        data.elementary_sum.coefficient[index] = data.elementary[index];
        data.alternating_elementary_sum.coefficient[index] = degree % 2 == 0
            ? data.elementary[index]
            : field.negate(data.elementary[index]);
    }
    trim(data.elementary_sum, limit);
    trim(data.alternating_elementary_sum, limit);

    const auto build_pairs = [largest_pair_index, limit, &field](
        const std::vector<std::uint32_t>& sequence
    ) {
        std::vector<Polynomial> pairs(
            static_cast<std::size_t>(largest_pair_index + 1)
        );
        for (int shift = 0; shift <= largest_pair_index; ++shift) {
            Polynomial value;
            for (int lower = 0; 2 * lower + shift <= limit; ++lower) {
                const int degree = 2 * lower + shift;
                const std::size_t index = static_cast<std::size_t>(degree);
                value.coefficient[index] = field.add(
                    value.coefficient[index],
                    field.multiply(
                        sequence[static_cast<std::size_t>(lower)],
                        sequence[static_cast<std::size_t>(lower + shift)]
                    )
                );
            }
            trim(value, limit);
            pairs[static_cast<std::size_t>(shift)] = value;
        }
        return pairs;
    };
    data.elementary_pairs = build_pairs(data.elementary);
    data.complete_pairs = build_pairs(data.complete);
    return data;
}

Polynomial symmetric_pair(
    const std::vector<Polynomial>& pairs,
    int index
) {
    if (index < 0) index = -index;
    if (index >= static_cast<int>(pairs.size())) {
        throw std::runtime_error("independent pair index overflow");
    }
    return pairs[static_cast<std::size_t>(index)];
}

std::vector<Polynomial> type_b_matrix(
    const std::vector<Polynomial>& pairs,
    int size,
    bool plus,
    const PrimeField& field,
    int limit
) {
    std::vector<Polynomial> matrix(static_cast<std::size_t>(size * size));
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            matrix[static_cast<std::size_t>(row * size + column)] =
                sum_polynomial(
                    symmetric_pair(pairs, row - column),
                    symmetric_pair(pairs, row + column + 1),
                    field,
                    limit,
                    !plus
                );
        }
    }
    return matrix;
}

std::vector<Polynomial> type_c_matrix(
    const std::vector<Polynomial>& pairs,
    int size,
    const PrimeField& field,
    int limit
) {
    std::vector<Polynomial> matrix(static_cast<std::size_t>(size * size));
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            matrix[static_cast<std::size_t>(row * size + column)] =
                sum_polynomial(
                    symmetric_pair(pairs, row - column),
                    symmetric_pair(pairs, row + column + 2),
                    field,
                    limit,
                    true
                );
        }
    }
    return matrix;
}

std::vector<Polynomial> type_d_matrix(
    const std::vector<Polynomial>& pairs,
    int size,
    const PrimeField& field,
    int limit
) {
    std::vector<Polynomial> matrix(static_cast<std::size_t>(size * size));
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            matrix[static_cast<std::size_t>(row * size + column)] =
                sum_polynomial(
                    symmetric_pair(pairs, row - column),
                    symmetric_pair(pairs, row + column),
                    field,
                    limit,
                    false
                );
        }
    }
    return matrix;
}

struct CandidateRow {
    char family;
    int rank;
    int tail_onset;
    int moment_through;
    std::vector<BigInt> moments;
};

int dimension(char family, int rank) {
    if (family == 'B' || family == 'C') return rank * (2 * rank + 1);
    if (family == 'D') return rank * (2 * rank - 1);
    throw std::runtime_error("unknown family in dimension");
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

std::vector<CandidateRow> empty_rows() {
    std::vector<CandidateRow> rows;
    for (const auto& cutoff : full_q3_bcd_remaining::row_cutoffs) {
        const bool original_subledger =
            (cutoff.family == 'B' && cutoff.rank <= 17)
            || cutoff.family == 'C'
            || (cutoff.family == 'D' && cutoff.rank <= 30);
        if (!original_subledger) continue;
        rows.push_back(CandidateRow{
            cutoff.family,
            cutoff.rank,
            cutoff.tail_onset,
            cutoff.moment_through,
            std::vector<BigInt>(
                static_cast<std::size_t>(cutoff.moment_through + 1), -1
            )
        });
    }
    return rows;
}

std::size_t row_index(
    const std::vector<CandidateRow>& rows,
    char family,
    int rank
) {
    for (std::size_t index = 0; index < rows.size(); ++index) {
        if (rows[index].family == family && rows[index].rank == rank) {
            return index;
        }
    }
    throw std::runtime_error("moment log contains an unknown row");
}

bool parse_moment_line(
    const std::string& line,
    const std::string& prefix,
    char& family,
    int& rank,
    int& degree,
    BigInt& value
) {
    if (line.rfind(prefix, 0U) != 0U) return false;
    std::istringstream input(line.substr(prefix.size()));
    std::string row_token;
    std::string degree_token;
    std::string value_token;
    if (!(input >> row_token >> degree_token >> value_token)) {
        throw std::runtime_error("malformed moment line");
    }
    if (row_token.rfind("row=", 0U) != 0U
        || degree_token.rfind("degree=", 0U) != 0U
        || value_token.rfind("value=", 0U) != 0U) {
        throw std::runtime_error("malformed moment fields");
    }
    const std::string row = row_token.substr(4U);
    if (row.size() < 3U || row[1] != '_') {
        throw std::runtime_error("malformed moment row");
    }
    family = row[0];
    rank = std::stoi(row.substr(2U));
    degree = std::stoi(degree_token.substr(7U));
    value = BigInt(value_token.substr(6U));
    return true;
}

std::vector<CandidateRow> read_candidates(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open bounded moment log");
    std::vector<CandidateRow> rows = empty_rows();
    std::string line;
    std::size_t count = 0U;
    while (std::getline(input, line)) {
        char family = 0;
        int rank = 0;
        int degree = 0;
        BigInt value;
        if (!parse_moment_line(
                line,
                "FULL_Q3_BOUNDED_MOMENT ",
                family,
                rank,
                degree,
                value
            )) {
            continue;
        }
        CandidateRow& row = rows[row_index(rows, family, rank)];
        if (degree < 0 || degree > row.moment_through) {
            throw std::runtime_error("bounded moment degree outside row ledger");
        }
        BigInt& target = row.moments[static_cast<std::size_t>(degree)];
        if (target >= 0) throw std::runtime_error("duplicate bounded moment");
        target = value;
        ++count;
    }
    std::size_t expected = 0U;
    for (const CandidateRow& row : rows) {
        expected += static_cast<std::size_t>(row.moment_through + 1);
        const int bound_base = dimension(row.family, row.rank);
        for (int degree = 0; degree <= row.moment_through; ++degree) {
            const BigInt& value = row.moments[static_cast<std::size_t>(degree)];
            if (value < 0 || value > integer_power(bound_base, degree)) {
                throw std::runtime_error("candidate moment missing or outside character bound");
            }
        }
    }
    if (count != expected || expected != 3383U) {
        throw std::runtime_error("bounded moment log scope mismatch");
    }
    return rows;
}

std::size_t verify_reverse_overlap(
    const std::vector<CandidateRow>& candidates,
    const std::string& path
) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open reverse-Pieri prefix log");
    std::vector<std::vector<bool>> seen(candidates.size());
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        seen[index].assign(candidates[index].moments.size(), false);
    }
    std::string line;
    std::size_t checked = 0U;
    while (std::getline(input, line)) {
        char family = 0;
        int rank = 0;
        int degree = 0;
        BigInt value;
        if (!parse_moment_line(
                line,
                "FULL_Q3_MOMENT ",
                family,
                rank,
                degree,
                value
            )) {
            continue;
        }
        if (degree > 40) continue;
        const std::size_t index = row_index(candidates, family, rank);
        const CandidateRow& row = candidates[index];
        if (degree < 0 || degree > row.moment_through) {
            throw std::runtime_error("reverse-Pieri degree outside row ledger");
        }
        if (seen[index][static_cast<std::size_t>(degree)]) {
            throw std::runtime_error("duplicate reverse-Pieri moment");
        }
        if (value != row.moments[static_cast<std::size_t>(degree)]) {
            throw std::runtime_error("reverse-Pieri overlap mismatch");
        }
        seen[index][static_cast<std::size_t>(degree)] = true;
        ++checked;
    }
    std::size_t expected = 0U;
    for (const CandidateRow& row : candidates) {
        expected += static_cast<std::size_t>(std::min(40, row.moment_through) + 1);
    }
    if (checked != expected || expected != 2845U) {
        throw std::runtime_error("reverse-Pieri overlap scope mismatch");
    }
    return checked;
}

bool is_prime(std::uint32_t value) {
    if (value < 2U) return false;
    if (value % 2U == 0U) return value == 2U;
    for (std::uint32_t divisor = 3U;
         static_cast<std::uint64_t>(divisor) * divisor <= value;
         divisor += 2U) {
        if (value % divisor == 0U) return false;
    }
    return true;
}

std::vector<std::uint32_t> independent_primes(const BigInt& required_product) {
    std::vector<std::uint32_t> primes;
    BigInt product = 1;
    std::uint32_t candidate = 1000000007U;
    while (product <= required_product) {
        if (is_prime(candidate)) {
            primes.push_back(candidate);
            product *= candidate;
        }
        if (candidate <= 1000003U) {
            throw std::runtime_error("independent prime search exhausted");
        }
        candidate -= 2U;
    }
    return primes;
}

// Convert values P(0),...,P(j) into the functional
// j! sum_a (2(j-a)-1)!! [z^a]P(z).  The implementation uses Newton forward
// differences P(z)=sum_k Delta^k P(0) binom(z,k), unlike the production
// Lagrange-node weight table.
std::uint32_t gaussian_functional_newton(
    std::vector<std::uint32_t> values,
    int moment,
    const PrimeField& field
) {
    if (values.size() != static_cast<std::size_t>(moment + 1)) {
        throw std::runtime_error("invalid Newton value count");
    }
    std::vector<std::uint32_t> differences(
        static_cast<std::size_t>(moment + 1), 0U
    );
    differences[0] = values[0];
    for (int order = 1; order <= moment; ++order) {
        for (int index = 0; index <= moment - order; ++index) {
            values[static_cast<std::size_t>(index)] = field.subtract(
                values[static_cast<std::size_t>(index + 1)],
                values[static_cast<std::size_t>(index)]
            );
        }
        differences[static_cast<std::size_t>(order)] = values[0];
    }

    std::vector<std::uint32_t> factorial(
        static_cast<std::size_t>(moment + 1), field.one()
    );
    std::vector<std::uint32_t> odd_double_factorial(
        static_cast<std::size_t>(moment + 1), field.one()
    );
    for (int degree = 1; degree <= moment; ++degree) {
        factorial[static_cast<std::size_t>(degree)] = field.multiply(
            factorial[static_cast<std::size_t>(degree - 1)],
            static_cast<std::uint32_t>(degree)
        );
        odd_double_factorial[static_cast<std::size_t>(degree)] = field.multiply(
            odd_double_factorial[static_cast<std::size_t>(degree - 1)],
            static_cast<std::uint32_t>(2 * degree - 1)
        );
    }

    std::vector<std::uint32_t> binomial_basis(1U, field.one());
    std::uint32_t answer = 0U;
    for (int order = 0; order <= moment; ++order) {
        std::uint32_t basis_functional = 0U;
        for (int power = 0;
             power < static_cast<int>(binomial_basis.size());
             ++power) {
            basis_functional = field.add(
                basis_functional,
                field.multiply(
                    binomial_basis[static_cast<std::size_t>(power)],
                    odd_double_factorial[
                        static_cast<std::size_t>(moment - power)
                    ]
                )
            );
        }
        basis_functional = field.multiply(
            basis_functional,
            factorial[static_cast<std::size_t>(moment)]
        );
        answer = field.add(
            answer,
            field.multiply(
                differences[static_cast<std::size_t>(order)],
                basis_functional
            )
        );

        if (order == moment) break;
        std::vector<std::uint32_t> next(binomial_basis.size() + 1U, 0U);
        const std::uint32_t negative_order = field.negate(
            static_cast<std::uint32_t>(order)
        );
        for (std::size_t degree = 0; degree < binomial_basis.size(); ++degree) {
            next[degree] = field.add(
                next[degree],
                field.multiply(binomial_basis[degree], negative_order)
            );
            next[degree + 1U] = field.add(
                next[degree + 1U], binomial_basis[degree]
            );
        }
        const std::uint32_t inverse_order = field.inverse(
            static_cast<std::uint32_t>(order + 1)
        );
        for (std::uint32_t& coefficient : next) {
            coefficient = field.multiply(coefficient, inverse_order);
        }
        binomial_basis = std::move(next);
    }
    return answer;
}

using NodeCoefficients = std::vector<std::vector<std::vector<std::uint32_t>>>;

NodeCoefficients determinant_coefficients(
    const std::vector<CandidateRow>& rows,
    const PrimeField& field
) {
    NodeCoefficients values(rows.size());
    for (std::size_t row = 0; row < rows.size(); ++row) {
        values[row].resize(static_cast<std::size_t>(rows[row].moment_through + 1));
        for (int moment = 0; moment <= rows[row].moment_through; ++moment) {
            values[row][static_cast<std::size_t>(moment)].assign(
                static_cast<std::size_t>(moment + 1), 0U
            );
        }
    }

    constexpr int b_maximum = 17;
    constexpr int c_maximum = 28;
    constexpr int d_maximum = 30;
    constexpr int largest_pair_index = 2 * d_maximum - 2;
    const std::uint32_t inverse_two = field.inverse(2U);

    for (int node_index = 0; node_index <= kMaximumMoment; ++node_index) {
        const NodeData node = make_node_data(
            static_cast<std::uint32_t>(node_index),
            largest_pair_index,
            field,
            kSeriesDegree
        );
        const auto b_minus = principal_determinants(
            type_b_matrix(
                node.elementary_pairs,
                b_maximum,
                false,
                field,
                kSeriesDegree
            ),
            b_maximum,
            field,
            kSeriesDegree
        );
        const auto b_plus = principal_determinants(
            type_b_matrix(
                node.elementary_pairs,
                b_maximum,
                true,
                field,
                kSeriesDegree
            ),
            b_maximum,
            field,
            kSeriesDegree
        );
        const auto c_values = principal_determinants(
            type_c_matrix(
                node.complete_pairs,
                c_maximum,
                field,
                kSeriesDegree
            ),
            c_maximum,
            field,
            kSeriesDegree
        );
        const auto d_values = principal_determinants(
            type_d_matrix(
                node.elementary_pairs,
                d_maximum,
                field,
                kSeriesDegree
            ),
            d_maximum,
            field,
            kSeriesDegree
        );

        for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
            const CandidateRow& row = rows[row_index];
            Polynomial generating;
            if (row.family == 'B') {
                const Polynomial minus = product_polynomial(
                    node.elementary_sum,
                    b_minus[static_cast<std::size_t>(row.rank)],
                    field,
                    kSeriesDegree
                );
                const Polynomial plus = product_polynomial(
                    node.alternating_elementary_sum,
                    b_plus[static_cast<std::size_t>(row.rank)],
                    field,
                    kSeriesDegree
                );
                generating = scale_polynomial(
                    sum_polynomial(
                        minus, plus, field, kSeriesDegree
                    ),
                    inverse_two,
                    field,
                    kSeriesDegree
                );
            } else if (row.family == 'C') {
                generating = c_values[static_cast<std::size_t>(row.rank)];
            } else if (row.family == 'D') {
                generating = scale_polynomial(
                    d_values[static_cast<std::size_t>(row.rank)],
                    inverse_two,
                    field,
                    kSeriesDegree
                );
            } else {
                throw std::runtime_error("unknown family in determinant evaluation");
            }
            for (int moment = node_index; moment <= row.moment_through; ++moment) {
                values[row_index][static_cast<std::size_t>(moment)]
                      [static_cast<std::size_t>(node_index)] =
                    generating.coefficient[
                        static_cast<std::size_t>(2 * moment)
                    ];
            }
        }
    }
    return values;
}

std::size_t check_prime(
    const std::vector<CandidateRow>& rows,
    std::uint32_t prime
) {
    const PrimeField field(prime);
    const NodeCoefficients values = determinant_coefficients(rows, field);
    std::size_t checked = 0U;
    std::size_t high_checked = 0U;
    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        const CandidateRow& row = rows[row_index];
        for (int moment = 0; moment <= row.moment_through; ++moment) {
            const std::uint32_t expected = gaussian_functional_newton(
                values[row_index][static_cast<std::size_t>(moment)],
                moment,
                field
            );
            const std::uint32_t candidate = static_cast<std::uint32_t>(
                mpz_fdiv_ui(
                    row.moments[static_cast<std::size_t>(moment)].get_mpz_t(),
                    prime
                )
            );
            if (expected != candidate) {
                throw std::runtime_error(
                    std::string("modular moment mismatch at ")
                    + row.family + '_' + std::to_string(row.rank)
                    + " degree=" + std::to_string(moment)
                    + " prime=" + std::to_string(prime)
                );
            }
            ++checked;
            if (moment > 40) ++high_checked;
        }
    }
    if (checked != 3383U || high_checked != 538U) {
        throw std::runtime_error("modular moment scope mismatch");
    }
    return high_checked;
}

BigInt binomial(int top, int bottom) {
    if (bottom < 0 || bottom > top) return 0;
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
    const int total = 2 * a + n;
    if (a < 0 || n < 0 || total >= static_cast<int>(moments.size())) {
        throw std::runtime_error("invalid compressed hierarchy index");
    }
    BigInt answer = 0;
    for (int j = 0; j <= 2 * a; ++j) {
        for (int k = 0; k <= n; ++k) {
            const BigInt term =
                binomial(2 * a, j) * binomial(n, k)
                * moments[static_cast<std::size_t>(2 * a - j + k)]
                * moments[static_cast<std::size_t>(j + n - k)];
            if (j % 2 == 0) {
                answer += term;
            } else {
                answer -= term;
            }
        }
    }
    return answer;
}

int first_odd_above(int value) {
    int answer = value + 1;
    if (answer % 2 == 0) ++answer;
    return answer;
}

std::pair<std::size_t, std::size_t> verify_hierarchy(
    const std::vector<CandidateRow>& rows
) {
    std::size_t prefix = 0U;
    std::size_t high = 0U;
    for (const CandidateRow& row : rows) {
        const int stable_through = row.family == 'D'
            ? row.rank - 1
            : 2 * row.rank + 1;
        for (int total = first_odd_above(stable_through);
             total <= row.tail_onset - 2;
             total += 2) {
            for (int a = 2; 2 * a < total; ++a) {
                const int n = total - 2 * a;
                if (hierarchy_value(row.moments, a, n) <= 0) {
                    throw std::runtime_error(
                        std::string("nonpositive compressed hierarchy value at ")
                        + row.family + '_' + std::to_string(row.rank)
                        + " a=" + std::to_string(a)
                        + " n=" + std::to_string(n)
                    );
                }
                if (total <= 40) {
                    ++prefix;
                } else {
                    ++high;
                }
            }
        }
    }
    if (prefix != 7484U || high != 5509U) {
        throw std::runtime_error("compressed hierarchy scope mismatch");
    }
    return {prefix, high};
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::string moments_log;
        std::string reverse_log;
        bool progress = false;
        int prime_limit = -1;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--moments-log" && index + 1 < argc) {
                moments_log = argv[++index];
            } else if (argument == "--reverse-log" && index + 1 < argc) {
                reverse_log = argv[++index];
            } else if (argument == "--prime-limit" && index + 1 < argc) {
                prime_limit = std::stoi(argv[++index]);
            } else if (argument == "--progress") {
                progress = true;
            } else {
                std::cerr << "usage: " << argv[0]
                          << " --moments-log PATH --reverse-log PATH"
                          << " [--prime-limit N] [--progress]\n";
                return EXIT_FAILURE;
            }
        }
        if (moments_log.empty() || reverse_log.empty()) {
            throw std::runtime_error("both certificate logs are required");
        }

        const std::vector<CandidateRow> rows = read_candidates(moments_log);
        const std::size_t overlap = verify_reverse_overlap(rows, reverse_log);
        const auto [prefix_pairs, high_pairs] = verify_hierarchy(rows);

        BigInt largest_bound = 0;
        for (const CandidateRow& row : rows) {
            for (int moment = 0; moment <= row.moment_through; ++moment) {
                largest_bound = std::max(
                    largest_bound,
                    integer_power(dimension(row.family, row.rank), moment)
                );
            }
        }
        std::vector<std::uint32_t> primes = independent_primes(largest_bound);
        if (prime_limit >= 0) {
            if (prime_limit < 1 || prime_limit > static_cast<int>(primes.size())) {
                throw std::runtime_error("prime limit outside required table");
            }
            primes.resize(static_cast<std::size_t>(prime_limit));
        }

        std::vector<std::string> errors(primes.size());
        const int prime_count = static_cast<int>(primes.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int prime_index = 0; prime_index < prime_count; ++prime_index) {
            const std::size_t index = static_cast<std::size_t>(prime_index);
            try {
                check_prime(rows, primes[index]);
                if (progress) {
#ifdef _OPENMP
#pragma omp critical(full_q3_modular_checker_progress)
#endif
                    std::cout << "FULL_Q3_MODULAR_CHECKER completed_prime="
                              << primes[index]
                              << " index=" << (prime_index + 1)
                              << '/' << prime_count << std::endl;
                }
            } catch (const std::exception& error) {
                errors[index] = error.what();
            } catch (...) {
                errors[index] = "unknown modular checker failure";
            }
        }
        for (const std::string& error : errors) {
            if (!error.empty()) throw std::runtime_error(error);
        }

        BigInt modulus_product = 1;
        for (std::uint32_t prime : primes) modulus_product *= prime;
        if (prime_limit < 0 && modulus_product <= largest_bound) {
            throw std::runtime_error("independent modulus does not dominate bound");
        }
        const std::size_t modulus_bits = mpz_sizeinbase(
            modulus_product.get_mpz_t(), 2
        );
        const std::size_t bound_bits = mpz_sizeinbase(
            largest_bound.get_mpz_t(), 2
        );

        std::cout << "FULL_Q3_MODULAR_CHECKER reverse_overlap_moments="
                  << overlap
                  << " modular_moments=3383 modular_high_moments=538"
                  << " prefix_pairs=" << prefix_pairs
                  << " high_pairs=" << high_pairs
                  << " primes=" << primes.size()
                  << " modulus_bits=" << modulus_bits
                  << " required_bound_bits=" << bound_bits
                  << " newton_interpolation=1 plain_modular_arithmetic=1\n";
        if (prime_limit >= 0) {
            std::cout << "FULL_Q3_MODULAR_CHECKER DIAGNOSTIC PRIME PREFIX: PASS\n";
        } else {
            std::cout << "FULL_Q3_MODULAR_CHECKER VERIFICATION: ALL PASS\n";
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "FULL_Q3_MODULAR_CHECKER VERIFICATION: FAIL: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
