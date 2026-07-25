#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_rational;
using Rat = cpp_rational;

static Rat quotient(const cpp_int& numerator, const cpp_int& denominator) {
    Rat result = numerator;
    result /= denominator;
    return result;
}

static Rat quotient(long numerator, long denominator) {
    return quotient(cpp_int(numerator), cpp_int(denominator));
}

struct Interval {
    Rat low;
    Rat high;
};

static Interval point(const Rat& value) {
    return {value, value};
}

static Interval add(const Interval& left, const Interval& right) {
    return {left.low + right.low, left.high + right.high};
}

static Interval subtract(const Interval& left, const Interval& right) {
    return {left.low - right.high, left.high - right.low};
}

static Interval negate(const Interval& value) {
    return {-value.high, -value.low};
}

static Interval multiply(const Interval& left, const Interval& right) {
    const std::array<Rat, 4> values{
        left.low * right.low,
        left.low * right.high,
        left.high * right.low,
        left.high * right.high,
    };
    return {
        *std::min_element(values.begin(), values.end()),
        *std::max_element(values.begin(), values.end()),
    };
}

static Interval scale(const Interval& value, const Rat& scalar) {
    if (scalar >= 0) {
        return {value.low * scalar, value.high * scalar};
    }
    return {value.high * scalar, value.low * scalar};
}

static Interval square(const Interval& value) {
    if (value.low >= 0) {
        return {value.low * value.low, value.high * value.high};
    }
    if (value.high <= 0) {
        return {value.high * value.high, value.low * value.low};
    }
    return {0, std::max(value.low * value.low, value.high * value.high)};
}

using Polynomial = std::vector<Rat>; // ascending coefficients

static void trim(Polynomial& polynomial) {
    while (polynomial.size() > 1U && polynomial.back() == 0) {
        polynomial.pop_back();
    }
}

static Polynomial derivative(const Polynomial& polynomial) {
    if (polynomial.size() <= 1U) {
        return {0};
    }
    Polynomial result(polynomial.size() - 1U);
    for (std::size_t index = 1; index < polynomial.size(); ++index) {
        result[index - 1U] = Rat(static_cast<unsigned long>(index)) * polynomial[index];
    }
    trim(result);
    return result;
}

static std::pair<Polynomial, Polynomial> divide_with_remainder(
    Polynomial dividend,
    const Polynomial& divisor
) {
    trim(dividend);
    Polynomial quotient_polynomial(
        dividend.size() >= divisor.size()
            ? dividend.size() - divisor.size() + 1U
            : 1U,
        0
    );
    while (
        !(dividend.size() == 1U && dividend[0] == 0) &&
        dividend.size() >= divisor.size()
    ) {
        const std::size_t offset = dividend.size() - divisor.size();
        const Rat coefficient = dividend.back() / divisor.back();
        quotient_polynomial[offset] = coefficient;
        for (std::size_t index = 0; index < divisor.size(); ++index) {
            dividend[offset + index] -= coefficient * divisor[index];
        }
        trim(dividend);
    }
    trim(quotient_polynomial);
    trim(dividend);
    return {quotient_polynomial, dividend};
}

static std::vector<Polynomial> sturm_sequence(const Polynomial& polynomial) {
    std::vector<Polynomial> sequence{polynomial, derivative(polynomial)};
    while (!(sequence.back().size() == 1U && sequence.back()[0] == 0)) {
        auto [quotient_polynomial, remainder] = divide_with_remainder(
            sequence[sequence.size() - 2U],
            sequence.back()
        );
        static_cast<void>(quotient_polynomial);
        for (Rat& coefficient : remainder) {
            coefficient = -coefficient;
        }
        trim(remainder);
        if (remainder.size() == 1U && remainder[0] == 0) {
            break;
        }
        sequence.push_back(remainder);
    }
    return sequence;
}

static Rat evaluate(const Polynomial& polynomial, const Rat& value) {
    Rat result = 0;
    for (auto iterator = polynomial.rbegin(); iterator != polynomial.rend(); ++iterator) {
        result = result * value + *iterator;
    }
    return result;
}

static int variations(const std::vector<Polynomial>& sequence, const Rat& value) {
    int previous = 0;
    int count = 0;
    for (const Polynomial& polynomial : sequence) {
        const Rat evaluated = evaluate(polynomial, value);
        const int sign = evaluated > 0 ? 1 : (evaluated < 0 ? -1 : 0);
        if (sign == 0) {
            continue;
        }
        if (previous != 0 && sign != previous) {
            ++count;
        }
        previous = sign;
    }
    return count;
}

static Rat decimal(const char* text) {
    std::string value(text);
    bool negative = false;
    if (!value.empty() && value.front() == '-') {
        negative = true;
        value.erase(value.begin());
    }
    const std::size_t point_position = value.find('.');
    const std::string integral = point_position == std::string::npos
        ? value
        : value.substr(0, point_position);
    const std::string fractional = point_position == std::string::npos
        ? std::string()
        : value.substr(point_position + 1U);
    cpp_int numerator = 0;
    for (char character : integral + fractional) {
        numerator = numerator * 10 + (character - '0');
    }
    cpp_int denominator = 1;
    for (std::size_t index = 0; index < fractional.size(); ++index) {
        denominator *= 10;
    }
    Rat result = quotient(numerator, denominator);
    return negative ? -result : result;
}

static Interval isolate_root(
    const std::vector<Polynomial>& sequence,
    Rat low,
    Rat high
) {
    if (variations(sequence, low) - variations(sequence, high) != 1) {
        throw std::runtime_error("invalid root bracket");
    }
    for (int iteration = 0; iteration < 180; ++iteration) {
        const Rat middle = (low + high) / 2;
        if (variations(sequence, low) - variations(sequence, middle) == 1) {
            high = middle;
        } else {
            low = middle;
        }
    }
    return {low, high};
}

static Interval evaluate_interval(const Polynomial& polynomial, const Interval& value) {
    Interval result = point(0);
    for (auto iterator = polynomial.rbegin(); iterator != polynomial.rend(); ++iterator) {
        result = add(multiply(result, value), point(*iterator));
    }
    return result;
}

struct Term {
    bool positive;
    Interval amount;
    Interval top_difference_squared;
    Interval first_sum_squared;
};

int main() {
    // f(x)=x^6-5x^5+5x^4+6x^3-7x^2-2x+1.
    const Polynomial minimal{1, -2, -7, 6, 5, -5, 1};
    const std::vector<Polynomial> sequence = sturm_sequence(minimal);
    const std::array<std::pair<const char*, const char*>, 6> brackets{{
        {"2.77", "2.78"},
        {"2.13", "2.14"},
        {"1.24", "1.25"},
        {"0.29", "0.30"},
        {"-0.50", "-0.49"},
        {"-0.95", "-0.94"},
    }};

    std::array<Interval, 6> first{};
    for (std::size_t index = 0; index < first.size(); ++index) {
        first[index] = isolate_root(
            sequence,
            decimal(brackets[index].first),
            decimal(brackets[index].second)
        );
    }

    // B_5=x^5-4x^4+2x^3+5x^2-2x-1 in the B_1 coordinate.
    const Polynomial top{-1, -2, 5, 2, -4, 1};
    std::array<Interval, 6> top_value{};
    std::array<Interval, 6> weight{};
    for (std::size_t index = 0; index < first.size(); ++index) {
        top_value[index] = evaluate_interval(top, first[index]);
        weight[index] = scale(
            subtract(point(3), first[index]),
            quotient(1, 13)
        );
    }

    std::vector<Term> terms;
    for (std::size_t left = 0; left < first.size(); ++left) {
        for (std::size_t right = left + 1U; right < first.size(); ++right) {
            const Interval difference = subtract(top_value[left], top_value[right]);
            const Interval sum = add(first[left], first[right]);
            if (!(sum.low > 0 || sum.high < 0)) {
                throw std::runtime_error("undetermined spectral sign");
            }
            const Interval absolute_sum = sum.low > 0 ? sum : negate(sum);
            const Interval amount = multiply(
                scale(multiply(weight[left], weight[right]), 2),
                multiply(square(difference), absolute_sum)
            );
            terms.push_back({
                sum.low > 0,
                amount,
                square(difference),
                square(sum),
            });
        }
    }

    std::vector<int> positives;
    std::vector<int> negatives;
    for (std::size_t index = 0; index < terms.size(); ++index) {
        (terms[index].positive ? positives : negatives).push_back(
            static_cast<int>(index)
        );
    }
    if (negatives.size() != 3U) {
        throw std::runtime_error("unexpected negative-pair count");
    }

    std::vector<std::vector<int>> edges(negatives.size());
    for (std::size_t negative_index = 0; negative_index < negatives.size(); ++negative_index) {
        const Term& negative = terms[static_cast<std::size_t>(negatives[negative_index])];
        for (std::size_t positive_index = 0; positive_index < positives.size(); ++positive_index) {
            const Term& positive = terms[static_cast<std::size_t>(positives[positive_index])];
            if (
                positive.top_difference_squared.low >= negative.top_difference_squared.high &&
                positive.first_sum_squared.low >= negative.first_sum_squared.high
            ) {
                edges[negative_index].push_back(static_cast<int>(positive_index));
            }
        }
        if (edges[negative_index].empty()) {
            throw std::runtime_error("negative pair has no dominant neighbour");
        }
    }

    for (unsigned mask = 1U; mask < (1U << negatives.size()); ++mask) {
        Interval demand = point(0);
        std::vector<bool> neighbour(positives.size(), false);
        for (std::size_t negative_index = 0; negative_index < negatives.size(); ++negative_index) {
            if ((mask & (1U << negative_index)) == 0U) {
                continue;
            }
            demand = add(
                demand,
                terms[static_cast<std::size_t>(negatives[negative_index])].amount
            );
            for (int positive_index : edges[negative_index]) {
                neighbour[static_cast<std::size_t>(positive_index)] = true;
            }
        }
        Interval capacity = point(0);
        for (std::size_t positive_index = 0; positive_index < positives.size(); ++positive_index) {
            if (neighbour[positive_index]) {
                capacity = add(
                    capacity,
                    terms[static_cast<std::size_t>(positives[positive_index])].amount
                );
            }
        }
        if (capacity.low < demand.high) {
            throw std::runtime_error("capacitated Hall inequality failed");
        }
    }

    std::cout
        << "SU2_O11_B1PLUS_TOPMINUS_TRANSPORT PASS"
        << " roots=6 pairs=15 negative_pairs=" << negatives.size()
        << " hall_subsets=7 residual_dimensions=2\n";
}
