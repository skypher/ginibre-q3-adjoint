#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gmpxx.h>

namespace cyclo {

using Poly = std::vector<mpz_class>;
using RationalPoly = std::vector<mpq_class>;

void trim(Poly& polynomial) {
    while (!polynomial.empty() && polynomial.back() == 0) {
        polynomial.pop_back();
    }
}

void trim(RationalPoly& polynomial) {
    while (!polynomial.empty() && polynomial.back() == 0) {
        polynomial.pop_back();
    }
}

Poly multiply_polynomials(const Poly& left, const Poly& right) {
    if (left.empty() || right.empty()) {
        return {};
    }
    Poly result(
        left.size() + right.size() - 1U,
        mpz_class{0}
    );
    for (std::size_t first = 0U; first < left.size(); ++first) {
        for (std::size_t second = 0U;
             second < right.size();
             ++second) {
            result[first + second] +=
                left[first] * right[second];
        }
    }
    trim(result);
    return result;
}

Poly divide_exact(const Poly& numerator, const Poly& denominator) {
    if (denominator.empty()) {
        throw std::runtime_error(
            "division by the zero polynomial"
        );
    }
    Poly remainder = numerator;
    Poly quotient(
        numerator.size() >= denominator.size()
            ? numerator.size() - denominator.size() + 1U
            : 0U,
        mpz_class{0}
    );
    while (
        !remainder.empty()
        && remainder.size() >= denominator.size()
    ) {
        if (remainder.back() % denominator.back() != 0) {
            throw std::runtime_error(
                "nonexact cyclotomic polynomial division"
            );
        }
        const mpz_class coefficient =
            remainder.back() / denominator.back();
        const std::size_t shift =
            remainder.size() - denominator.size();
        quotient[shift] = coefficient;
        for (std::size_t index = 0U;
             index < denominator.size();
             ++index) {
            remainder[shift + index] -=
                coefficient * denominator[index];
        }
        trim(remainder);
    }
    if (!remainder.empty()) {
        throw std::runtime_error(
            "cyclotomic polynomial division left a remainder"
        );
    }
    trim(quotient);
    return quotient;
}

Poly cyclotomic_polynomial(int order) {
    Poly result(
        static_cast<std::size_t>(order + 1),
        mpz_class{0}
    );
    result[0U] = -1;
    result[static_cast<std::size_t>(order)] = 1;
    for (int divisor = 1; divisor < order; ++divisor) {
        if (order % divisor == 0) {
            result = divide_exact(
                result,
                cyclotomic_polynomial(divisor)
            );
        }
    }
    return result;
}

RationalPoly derivative(const RationalPoly& polynomial) {
    if (polynomial.size() <= 1U) {
        return {};
    }
    RationalPoly result(polynomial.size() - 1U, mpq_class{0});
    for (std::size_t degree = 1U;
         degree < polynomial.size();
         ++degree) {
        result[degree - 1U] =
            polynomial[degree] * mpq_class{degree};
    }
    trim(result);
    return result;
}

RationalPoly remainder(
    RationalPoly numerator,
    const RationalPoly& denominator
) {
    while (
        !numerator.empty()
        && numerator.size() >= denominator.size()
    ) {
        const mpq_class coefficient =
            numerator.back() / denominator.back();
        const std::size_t shift =
            numerator.size() - denominator.size();
        for (std::size_t index = 0U;
             index < denominator.size();
             ++index) {
            numerator[shift + index] -=
                coefficient * denominator[index];
        }
        trim(numerator);
    }
    return numerator;
}

mpq_class evaluate(
    const RationalPoly& polynomial,
    const mpq_class& point
) {
    mpq_class result{0};
    for (std::size_t index = polynomial.size();
         index > 0U;
         --index) {
        result *= point;
        result += polynomial[index - 1U];
    }
    return result;
}

std::vector<RationalPoly> sturm_sequence(const Poly& polynomial) {
    RationalPoly first;
    first.reserve(polynomial.size());
    for (const mpz_class& coefficient : polynomial) {
        first.emplace_back(coefficient);
    }
    std::vector<RationalPoly> result;
    result.push_back(first);
    result.push_back(derivative(first));
    while (result.back().size() > 1U) {
        RationalPoly next = remainder(
            result[result.size() - 2U],
            result.back()
        );
        for (mpq_class& coefficient : next) {
            coefficient = -coefficient;
        }
        trim(next);
        result.push_back(std::move(next));
    }
    return result;
}

int sign_of(const mpq_class& value) {
    return mpq_sgn(value.get_mpq_t());
}

int sturm_variations(
    const std::vector<RationalPoly>& sequence,
    const mpq_class& point
) {
    int previous = 0;
    int variations = 0;
    for (const RationalPoly& polynomial : sequence) {
        const int current = sign_of(evaluate(polynomial, point));
        if (current == 0) {
            continue;
        }
        if (previous != 0 && current != previous) {
            ++variations;
        }
        previous = current;
    }
    return variations;
}

int roots_between(
    const std::vector<RationalPoly>& sequence,
    const mpq_class& lower,
    const mpq_class& upper
) {
    return
        sturm_variations(sequence, lower)
        - sturm_variations(sequence, upper);
}

std::pair<mpq_class, mpq_class> evaluate_interval(
    const Poly& polynomial,
    const mpq_class& lower,
    const mpq_class& upper
) {
    mpq_class value_lower{0};
    mpq_class value_upper{0};
    for (std::size_t index = polynomial.size();
         index > 0U;
         --index) {
        const std::vector<mpq_class> products{
            value_lower * lower,
            value_lower * upper,
            value_upper * lower,
            value_upper * upper
        };
        const auto bounds = std::minmax_element(
            products.begin(),
            products.end()
        );
        value_lower =
            *bounds.first + mpq_class{polynomial[index - 1U]};
        value_upper =
            *bounds.second + mpq_class{polynomial[index - 1U]};
    }
    return {value_lower, value_upper};
}

struct Field {
    using Elt = Poly;

    int degree = 0;
    Poly minimal_polynomial;
    std::vector<RationalPoly> sturm;

    void init(int period) {
        Poly cyclotomic = cyclotomic_polynomial(2 * period);
        degree = (static_cast<int>(cyclotomic.size()) - 1) / 2;
        Poly residual = cyclotomic;
        minimal_polynomial.assign(
            static_cast<std::size_t>(degree + 1),
            mpz_class{0}
        );
        const Poly z_squared_plus_one{
            mpz_class{1},
            mpz_class{0},
            mpz_class{1}
        };
        for (int index = degree; index >= 0; --index) {
            const std::size_t top =
                static_cast<std::size_t>(degree + index);
            const mpz_class coefficient =
                top < residual.size()
                    ? residual[top]
                    : mpz_class{0};
            minimal_polynomial[
                static_cast<std::size_t>(index)
            ] = coefficient;
            if (coefficient == 0) {
                continue;
            }
            Poly term{mpz_class{1}};
            for (int exponent = 0;
                 exponent < index;
                 ++exponent) {
                term = multiply_polynomials(
                    term,
                    z_squared_plus_one
                );
            }
            Poly shifted(
                static_cast<std::size_t>(degree - index),
                mpz_class{0}
            );
            shifted.insert(shifted.end(), term.begin(), term.end());
            if (shifted.size() > residual.size()) {
                residual.resize(shifted.size(), mpz_class{0});
            }
            for (std::size_t position = 0U;
                 position < shifted.size();
                 ++position) {
                residual[position] -=
                    coefficient * shifted[position];
            }
            trim(residual);
        }
        if (
            !residual.empty()
            || minimal_polynomial.back() != 1
        ) {
            throw std::runtime_error(
                "real cyclotomic polynomial extraction failed"
            );
        }
        sturm = sturm_sequence(minimal_polynomial);
        const mpq_class lower{199, 100};
        const mpq_class upper{2};
        if (roots_between(sturm, lower, upper) != 1) {
            throw std::runtime_error(
                "the rational interval does not isolate zeta"
            );
        }
    }

    Elt fromInt(long value) const {
        return
            value == 0
                ? Elt{}
                : Elt{mpz_class{value}};
    }

    Elt reduce(Poly polynomial) const {
        while (
            static_cast<int>(polynomial.size()) > degree
        ) {
            const mpz_class coefficient = polynomial.back();
            const std::size_t shift =
                polynomial.size() - minimal_polynomial.size();
            for (std::size_t index = 0U;
                 index < minimal_polynomial.size();
                 ++index) {
                polynomial[shift + index] -=
                    coefficient * minimal_polynomial[index];
            }
            trim(polynomial);
        }
        return polynomial;
    }

    Elt add(const Elt& left, const Elt& right) const {
        Elt result = left;
        if (right.size() > result.size()) {
            result.resize(right.size(), mpz_class{0});
        }
        for (std::size_t index = 0U;
             index < right.size();
             ++index) {
            result[index] += right[index];
        }
        trim(result);
        return result;
    }

    Elt sub(const Elt& left, const Elt& right) const {
        Elt result = left;
        if (right.size() > result.size()) {
            result.resize(right.size(), mpz_class{0});
        }
        for (std::size_t index = 0U;
             index < right.size();
             ++index) {
            result[index] -= right[index];
        }
        trim(result);
        return result;
    }

    Elt mulE(const Elt& left, const Elt& right) const {
        return reduce(multiply_polynomials(left, right));
    }

    Elt powE(Elt base, unsigned long exponent) const {
        Elt result = fromInt(1);
        while (exponent != 0UL) {
            if ((exponent & 1UL) != 0UL) {
                result = mulE(result, base);
            }
            base = mulE(base, base);
            exponent >>= 1U;
        }
        return result;
    }

    static bool isZero(const Elt& element) {
        return element.empty();
    }

    Elt dickson(int index) const {
        Elt previous = fromInt(2);
        Elt current{
            mpz_class{0},
            mpz_class{1}
        };
        if (index == 0) {
            return previous;
        }
        const Elt zeta{
            mpz_class{0},
            mpz_class{1}
        };
        for (int step = 1; step < index; ++step) {
            Elt next = sub(mulE(current, zeta), previous);
            previous = std::move(current);
            current = std::move(next);
        }
        return current;
    }

    Elt chebyshevS(int index, const Elt& point) const {
        Elt previous = fromInt(1);
        Elt current = point;
        if (index == 0) {
            return previous;
        }
        for (int step = 1; step < index; ++step) {
            Elt next = sub(mulE(point, current), previous);
            previous = std::move(current);
            current = std::move(next);
        }
        return current;
    }

    int sign(const Elt& element) const {
        if (isZero(element)) {
            return 0;
        }
        mpq_class lower{199, 100};
        mpq_class upper{2};
        for (int iteration = 0; iteration < 4096; ++iteration) {
            const auto interval =
                evaluate_interval(element, lower, upper);
            if (sign_of(interval.first) > 0) {
                return 1;
            }
            if (sign_of(interval.second) < 0) {
                return -1;
            }
            const mpq_class midpoint = (lower + upper) / 2;
            const int left_roots =
                roots_between(sturm, lower, midpoint);
            if (left_roots == 1) {
                upper = midpoint;
            } else {
                lower = midpoint;
            }
        }
        throw std::runtime_error(
            "exact Sturm sign isolation did not terminate"
        );
    }
};

}  // namespace cyclo

namespace {

using ExactMatrix = std::vector<std::vector<mpz_class>>;

struct SpectralPair {
    int left = 0;
    int right = 0;
    long double approximate_base = 0.0L;
    cyclo::Field::Elt base;
    cyclo::Field::Elt residue;
};

struct ExactCallResult {
    std::size_t negative = 0U;
    std::size_t zero = 0U;
    bool has_negative = false;
    std::size_t negative_index = 0U;
    cyclo::Field::Elt negative_call;
};

ExactMatrix multiply(
    const ExactMatrix& left,
    const ExactMatrix& right
) {
    const std::size_t dimension = left.size();
    ExactMatrix result(
        dimension,
        std::vector<mpz_class>(dimension, mpz_class{0})
    );
    for (std::size_t row = 0U; row < dimension; ++row) {
        for (std::size_t inner = 0U; inner < dimension; ++inner) {
            if (left[row][inner] == 0) {
                continue;
            }
            for (std::size_t column = 0U;
                 column < dimension;
                 ++column) {
                result[row][column] +=
                    left[row][inner] * right[inner][column];
            }
        }
    }
    return result;
}

ExactMatrix fusion_matrix(int half_level, int factor) {
    const int dimension = half_level + 1;
    ExactMatrix result(
        static_cast<std::size_t>(dimension),
        std::vector<mpz_class>(
            static_cast<std::size_t>(dimension),
            mpz_class{0}
        )
    );
    for (int row = 0; row < dimension; ++row) {
        const int lower = std::abs(row - factor);
        const int upper = std::min(
            row + factor,
            2 * half_level - row - factor
        );
        for (int column = lower; column <= upper; ++column) {
            result[static_cast<std::size_t>(row)][
                static_cast<std::size_t>(column)
            ] = 1;
        }
    }
    return result;
}

cyclo::Field::Elt sine_pair_scaled(
    const cyclo::Field& field,
    int mode,
    int first_frequency,
    int second_frequency
) {
    const int difference =
        std::abs(first_frequency - second_frequency) * mode;
    const int sum =
        (first_frequency + second_frequency) * mode;
    return field.sub(
        field.dickson(difference),
        field.dickson(sum)
    );
}

cyclo::Field::Elt scaled_wedge_residue(
    const cyclo::Field& field,
    int half_level,
    int factor,
    int target,
    int left_mode,
    int right_mode
) {
    const int left = left_mode + 1;
    const int right = right_mode + 1;
    const int factor_frequency = 2 * factor + 1;
    const int target_frequency = 2 * target + 1;

    const cyclo::Field::Elt first = field.mulE(
        sine_pair_scaled(field, left, 1, 1),
        sine_pair_scaled(
            field,
            right,
            factor_frequency,
            target_frequency
        )
    );
    const cyclo::Field::Elt second = field.mulE(
        sine_pair_scaled(
            field,
            left,
            1,
            target_frequency
        ),
        sine_pair_scaled(
            field,
            right,
            factor_frequency,
            1
        )
    );
    const cyclo::Field::Elt third = field.mulE(
        sine_pair_scaled(
            field,
            left,
            factor_frequency,
            1
        ),
        sine_pair_scaled(
            field,
            right,
            1,
            target_frequency
        )
    );
    const cyclo::Field::Elt fourth = field.mulE(
        sine_pair_scaled(
            field,
            left,
            factor_frequency,
            target_frequency
        ),
        sine_pair_scaled(field, right, 1, 1)
    );

    cyclo::Field::Elt result = field.add(
        field.sub(field.sub(first, second), third),
        fourth
    );
    const long left_pairing =
        left_mode == half_level ? 1L : 2L;
    const long right_pairing =
        right_mode == half_level ? 1L : 2L;
    result = field.mulE(
        result,
        field.fromInt(left_pairing * right_pairing)
    );
    return result;
}

cyclo::Field::Elt squared_pair_eigenvalue(
    const cyclo::Field& field,
    int factor,
    int left_mode,
    int right_mode
) {
    const cyclo::Field::Elt left_node =
        field.dickson(left_mode + 1);
    const cyclo::Field::Elt right_node =
        field.dickson(right_mode + 1);
    const cyclo::Field::Elt left_eigenvalue =
        field.chebyshevS(2 * factor, left_node);
    const cyclo::Field::Elt right_eigenvalue =
        field.chebyshevS(2 * factor, right_node);
    const cyclo::Field::Elt product =
        field.mulE(left_eigenvalue, right_eigenvalue);
    return field.mulE(product, product);
}

long double approximate_squared_pair_eigenvalue(
    int period,
    int factor,
    int left_mode,
    int right_mode
) {
    const long double pi = std::acos(-1.0L);
    const auto eigenvalue =
        [period, factor, pi](int mode) {
            const long double theta =
                static_cast<long double>(mode + 1)
                * pi
                / static_cast<long double>(period);
            return
                std::sin(
                    static_cast<long double>(2 * factor + 1)
                    * theta
                )
                / std::sin(theta);
        };
    const long double product =
        eigenvalue(left_mode) * eigenvalue(right_mode);
    return product * product;
}

void print_element(
    const std::string& name,
    const cyclo::Field::Elt& element
) {
    std::cout << ' ' << name << "=[";
    for (std::size_t index = 0U;
         index < element.size();
         ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << element[index];
    }
    std::cout << ']';
}

std::vector<SpectralPair> build_exact_groups(
    cyclo::Field& field,
    int level,
    int factor,
    int target
) {
    const int half_level = level / 2;
    const int period = level + 2;
    std::vector<SpectralPair> pairs;
    for (int left = 0; left <= half_level; ++left) {
        for (int right = left + 1;
             right <= half_level;
             ++right) {
            pairs.push_back(
                SpectralPair{
                    left,
                    right,
                    approximate_squared_pair_eigenvalue(
                        period,
                        factor,
                        left,
                        right
                    ),
                    squared_pair_eigenvalue(
                        field,
                        factor,
                        left,
                        right
                    ),
                    scaled_wedge_residue(
                        field,
                        half_level,
                        factor,
                        target,
                        left,
                        right
                    )
                }
            );
        }
    }
    std::sort(
        pairs.begin(),
        pairs.end(),
        [](const SpectralPair& left, const SpectralPair& right) {
            return left.approximate_base < right.approximate_base;
        }
    );
    std::vector<SpectralPair> groups;
    for (const SpectralPair& pair : pairs) {
        if (
            groups.empty()
            || !cyclo::Field::isZero(
                field.sub(groups.back().base, pair.base)
            )
        ) {
            groups.push_back(pair);
        } else {
            groups.back().residue = field.add(
                groups.back().residue,
                pair.residue
            );
        }
    }
    for (std::size_t index = 1U;
         index < groups.size();
         ++index) {
        if (
            field.sign(
                field.sub(
                    groups[index].base,
                    groups[index - 1U].base
                )
            ) <= 0
        ) {
            throw std::runtime_error(
                "exact spectral-base ordering failed"
            );
        }
    }
    return groups;
}

cyclo::Field::Elt spectral_total(
    const cyclo::Field& field,
    const std::vector<SpectralPair>& groups,
    int shift
) {
    cyclo::Field::Elt total;
    for (const SpectralPair& group : groups) {
        total = field.add(
            total,
            field.mulE(
                group.residue,
                field.powE(
                    group.base,
                    static_cast<unsigned long>(shift)
                )
            )
        );
    }
    return total;
}

ExactCallResult power_convex_call_signs(
    const cyclo::Field& field,
    const std::vector<SpectralPair>& groups,
    int shift
) {
    ExactCallResult result;
    for (std::size_t cutoff = 0U;
         cutoff < groups.size();
         ++cutoff) {
        cyclo::Field::Elt call;
        for (std::size_t index = cutoff;
             index < groups.size();
             ++index) {
            call = field.add(
                call,
                field.mulE(
                    field.mulE(
                        groups[index].residue,
                        field.powE(
                            groups[index].base,
                            static_cast<unsigned long>(shift)
                        )
                    ),
                    field.sub(
                        groups[index].base,
                        groups[cutoff].base
                    )
                )
            );
        }
        const int sign = field.sign(call);
        if (sign < 0) {
            ++result.negative;
            if (!result.has_negative) {
                result.has_negative = true;
                result.negative_index = cutoff;
                result.negative_call = call;
            }
        } else if (sign == 0) {
            ++result.zero;
        }
    }
    return result;
}

mpz_class exact_anchored_determinant(
    int half_level,
    int factor,
    int target,
    int shift
) {
    const ExactMatrix fusion = fusion_matrix(half_level, factor);
    const ExactMatrix square = multiply(fusion, fusion);
    ExactMatrix power(
        static_cast<std::size_t>(half_level + 1),
        std::vector<mpz_class>(
            static_cast<std::size_t>(half_level + 1),
            mpz_class{0}
        )
    );
    for (int index = 0; index <= half_level; ++index) {
        power[static_cast<std::size_t>(index)][
            static_cast<std::size_t>(index)
        ] = 1;
    }
    for (int exponent = 0; exponent < shift; ++exponent) {
        power = multiply(power, square);
    }
    return
        power[0U][0U]
            * power[static_cast<std::size_t>(factor)][
                static_cast<std::size_t>(target)
            ]
        - power[0U][static_cast<std::size_t>(target)]
            * power[static_cast<std::size_t>(factor)][0U];
}

}  // namespace

int main() {
    try {
        constexpr int level = 38;
        constexpr int half_level = level / 2;
        constexpr int period = level + 2;
        constexpr int factor = 2;
        constexpr int target = 17;
        constexpr int shift = 5;

        cyclo::Field field;
        field.init(period);

        std::vector<SpectralPair> pairs;
        for (int left = 0; left <= half_level; ++left) {
            for (int right = left + 1;
                 right <= half_level;
                 ++right) {
                pairs.push_back(
                    SpectralPair{
                        left,
                        right,
                        approximate_squared_pair_eigenvalue(
                            period,
                            factor,
                            left,
                            right
                        ),
                        squared_pair_eigenvalue(
                            field,
                            factor,
                            left,
                            right
                        ),
                        scaled_wedge_residue(
                            field,
                            half_level,
                            factor,
                            target,
                            left,
                            right
                        )
                    }
                );
            }
        }
        std::sort(
            pairs.begin(),
            pairs.end(),
            [](const SpectralPair& left, const SpectralPair& right) {
                return
                    left.approximate_base < right.approximate_base;
            }
        );

        std::vector<SpectralPair> groups;
        for (const SpectralPair& pair : pairs) {
            if (
                groups.empty()
                || !cyclo::Field::isZero(
                    field.sub(groups.back().base, pair.base)
                )
            ) {
                groups.push_back(pair);
            } else {
                groups.back().residue = field.add(
                    groups.back().residue,
                    pair.residue
                );
            }
        }
        for (std::size_t index = 1U;
             index < groups.size();
             ++index) {
            if (
                field.sign(
                    field.sub(
                        groups[index].base,
                        groups[index - 1U].base
                    )
                ) <= 0
            ) {
                throw std::runtime_error(
                    "exact spectral-base ordering failed"
                );
            }
        }

        cyclo::Field::Elt total;
        cyclo::Field::Elt tail;
        bool found_negative = false;
        std::size_t negative_index = 0U;
        cyclo::Field::Elt negative_tail;
        for (std::size_t reverse = groups.size();
             reverse > 0U;
             --reverse) {
            const SpectralPair& group = groups[reverse - 1U];
            const cyclo::Field::Elt weighted = field.mulE(
                group.residue,
                field.powE(
                    group.base,
                    static_cast<unsigned long>(shift)
                )
            );
            tail = field.add(tail, weighted);
            if (!found_negative && field.sign(tail) < 0) {
                found_negative = true;
                negative_index = reverse - 1U;
                negative_tail = tail;
            }
        }
        for (const SpectralPair& group : groups) {
            total = field.add(
                total,
                field.mulE(
                    group.residue,
                    field.powE(
                        group.base,
                        static_cast<unsigned long>(shift)
                    )
                )
            );
        }

        ExactMatrix fusion = fusion_matrix(half_level, factor);
        const ExactMatrix square = multiply(fusion, fusion);
        ExactMatrix power(
            static_cast<std::size_t>(half_level + 1),
            std::vector<mpz_class>(
                static_cast<std::size_t>(half_level + 1),
                mpz_class{0}
            )
        );
        for (int index = 0; index <= half_level; ++index) {
            power[static_cast<std::size_t>(index)][
                static_cast<std::size_t>(index)
            ] = 1;
        }
        for (int exponent = 0; exponent < shift; ++exponent) {
            power = multiply(power, square);
        }
        const mpz_class determinant =
            power[0U][0U]
                * power[static_cast<std::size_t>(factor)][
                    static_cast<std::size_t>(target)
                ]
            - power[0U][static_cast<std::size_t>(target)]
                * power[static_cast<std::size_t>(factor)][0U];
        const mpz_class scale =
            16L * (half_level + 1) * (half_level + 1);
        const cyclo::Field::Elt expected{
            scale * determinant
        };
        if (!cyclo::Field::isZero(field.sub(total, expected))) {
            throw std::runtime_error(
                "exact spectral total disagrees with fusion determinant"
            );
        }
        if (!found_negative) {
            throw std::runtime_error(
                "expected a negative shifted upper tail"
            );
        }

        const SpectralPair& cutoff = groups[negative_index];
        std::cout
            << std::setprecision(20)
            << "SU2_FINITE_ANCHORED_SHIFTED_TAIL_COUNTEREXAMPLE"
            << " level=" << level
            << " factor=" << factor
            << " target=" << target
            << " shift=" << shift
            << " field_degree=" << field.degree
            << " pairs=" << pairs.size()
            << " exact_groups=" << groups.size()
            << " cutoff_pair=("
            << cutoff.left << ',' << cutoff.right << ')'
            << " cutoff_approx="
            << cutoff.approximate_base
            << " determinant=" << determinant
            << " total_scale=" << scale
            << " tail_sign=" << field.sign(negative_tail)
            << " result=PASS_EXACT_COUNTEREXAMPLE";
        print_element("cutoff", cutoff.base);
        print_element("negative_tail", negative_tail);
        std::cout << '\n';

        const auto interior_call_signs =
            power_convex_call_signs(field, groups, 4);
        if (interior_call_signs.negative != 0U) {
            throw std::runtime_error(
                "the exact interior power-convex case failed"
            );
        }
        std::cout
            << "SU2_FINITE_ANCHORED_POWER_CONVEX_CASE"
            << " level=38 factor=2 target=17 shift=4"
            << " exact_groups=" << groups.size()
            << " negative_calls=" << interior_call_signs.negative
            << " zero_calls=" << interior_call_signs.zero
            << " result=PASS_EXACT_CASE"
            << '\n';

        constexpr int wall_factor = 9;
        constexpr int wall_target = 19;
        constexpr int wall_shift = 4;
        const std::vector<SpectralPair> wall_groups =
            build_exact_groups(
                field,
                level,
                wall_factor,
                wall_target
            );
        const mpz_class wall_determinant =
            exact_anchored_determinant(
                half_level,
                wall_factor,
                wall_target,
                wall_shift
            );
        const cyclo::Field::Elt wall_expected{
            scale * wall_determinant
        };
        if (
            !cyclo::Field::isZero(
                field.sub(
                    spectral_total(
                        field,
                        wall_groups,
                        wall_shift
                    ),
                    wall_expected
                )
            )
        ) {
            throw std::runtime_error(
                "the exact wall spectral total is misnormalized"
            );
        }
        const auto wall_call_signs =
            power_convex_call_signs(
                field,
                wall_groups,
                wall_shift
            );
        if (wall_call_signs.negative != 0U) {
            throw std::runtime_error(
                "the exact wall power-convex case failed"
            );
        }
        std::cout
            << "SU2_FINITE_ANCHORED_POWER_CONVEX_CASE"
            << " level=38 factor=" << wall_factor
            << " target=" << wall_target
            << " shift=" << wall_shift
            << " exact_groups=" << wall_groups.size()
            << " determinant=" << wall_determinant
            << " negative_calls=" << wall_call_signs.negative
            << " zero_calls=" << wall_call_signs.zero
            << " result=PASS_EXACT_CASE"
            << '\n';

        constexpr int closure_level = 40;
        constexpr int closure_half_level = closure_level / 2;
        constexpr int closure_factor = 1;
        constexpr int closure_target = 18;
        constexpr int closure_shift = 10;
        cyclo::Field closure_field;
        closure_field.init(closure_level + 2);
        const std::vector<SpectralPair> closure_groups =
            build_exact_groups(
                closure_field,
                closure_level,
                closure_factor,
                closure_target
            );
        const ExactCallResult closure_call_signs =
            power_convex_call_signs(
                closure_field,
                closure_groups,
                closure_shift
            );
        if (!closure_call_signs.has_negative) {
            throw std::runtime_error(
                "expected an exact closure-anchor call counterexample"
            );
        }
        const mpz_class closure_determinant =
            exact_anchored_determinant(
                closure_half_level,
                closure_factor,
                closure_target,
                closure_shift
            );
        const mpz_class closure_scale =
            16L
            * (closure_half_level + 1)
            * (closure_half_level + 1);
        const cyclo::Field::Elt closure_expected{
            closure_scale * closure_determinant
        };
        if (
            !cyclo::Field::isZero(
                closure_field.sub(
                    spectral_total(
                        closure_field,
                        closure_groups,
                        closure_shift
                    ),
                    closure_expected
                )
            )
        ) {
            throw std::runtime_error(
                "closure-anchor spectral total is misnormalized"
            );
        }
        const SpectralPair& closure_cutoff =
            closure_groups[closure_call_signs.negative_index];
        std::cout
            << "SU2_FINITE_ANCHORED_POWER_CONVEX_COUNTEREXAMPLE"
            << " level=" << closure_level
            << " factor=" << closure_factor
            << " target=" << closure_target
            << " shift=" << closure_shift
            << " field_degree=" << closure_field.degree
            << " exact_groups=" << closure_groups.size()
            << " cutoff_pair=("
            << closure_cutoff.left << ','
            << closure_cutoff.right << ')'
            << " cutoff_approx="
            << closure_cutoff.approximate_base
            << " determinant=" << closure_determinant
            << " total_scale=" << closure_scale
            << " negative_calls=" << closure_call_signs.negative
            << " zero_calls=" << closure_call_signs.zero
            << " call_sign="
            << closure_field.sign(
                closure_call_signs.negative_call
            )
            << " result=PASS_EXACT_COUNTEREXAMPLE";
        print_element("cutoff", closure_cutoff.base);
        print_element(
            "negative_call",
            closure_call_signs.negative_call
        );
        std::cout << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr
            << "SU2_FINITE_ANCHORED_SHIFTED_TAIL_COUNTEREXAMPLE"
            << " result=FAILURE error=" << exception.what()
            << '\n';
        return EXIT_FAILURE;
    }
}
