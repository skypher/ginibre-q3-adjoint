#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

using boost::multiprecision::cpp_int;

namespace {

constexpr int deficit = 10;
constexpr int label_count = 5;
using Counts = std::array<int, label_count>;
using Directions = std::array<Counts, label_count>;
using MultiplicityKey = std::array<int, label_count + 1>;
using OuterKey = std::array<int, 2 * label_count>;
using RecursiveKey = std::array<int, label_count + 1>;
using BivariateMonomial
    = std::array<std::array<cpp_int, deficit + 1>, deficit + 1>;
using MultivariateMonomial = std::map<Counts, cpp_int>;
using Rational = boost::rational<cpp_int>;
using RationalPolynomial = std::vector<Rational>;

std::map<MultiplicityKey, cpp_int>& multiplicity_cache() {
    static std::map<MultiplicityKey, cpp_int> cache;
    return cache;
}

std::map<OuterKey, cpp_int>& outer_cache() {
    static std::map<OuterKey, cpp_int> cache;
    return cache;
}

long long binomial(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    long long result = 1;
    for (int index = 1; index <= k; ++index) {
        result = result * (n - k + index) / index;
    }
    return result;
}

int total_label(const Counts& counts) {
    int result = 0;
    for (int index = 0; index < label_count; ++index) {
        result += (index + 1) * counts[static_cast<std::size_t>(index)];
    }
    return result;
}

int total_count(const Counts& counts) {
    int result = 0;
    for (int count : counts) {
        result += count;
    }
    return result;
}

int maximum_label(const Counts& counts) {
    for (int index = label_count - 1; index >= 0; --index) {
        if (counts[static_cast<std::size_t>(index)] != 0) {
            return index + 1;
        }
    }
    return 0;
}

bool admissible(const Counts& counts) {
    const int maximum = maximum_label(counts);
    return maximum != 0 && total_label(counts) - deficit >= maximum;
}

cpp_int multiplicity(const Counts& counts, int target) {
    MultiplicityKey key{};
    for (std::size_t index = 0U; index < counts.size(); ++index) {
        key[index] = counts[index];
    }
    key[counts.size()] = target;
    auto& cache = multiplicity_cache();
    const auto found = cache.find(key);
    if (found != cache.end()) {
        return found->second;
    }
    const int total = total_label(counts);
    if (target < 0 || target > total) {
        return 0;
    }
    std::vector<cpp_int> current(static_cast<std::size_t>(total + 1));
    current[0U] = 1;
    int support = 0;
    for (int label = 1; label <= label_count; ++label) {
        for (int copy = 0;
             copy < counts[static_cast<std::size_t>(label - 1)]; ++copy) {
            std::vector<cpp_int> next(current.size());
            for (int input = 0; input <= support; ++input) {
                if (current[static_cast<std::size_t>(input)] == 0) {
                    continue;
                }
                for (int output = std::abs(input - label);
                     output <= input + label; output += 2) {
                    next[static_cast<std::size_t>(output)]
                        += current[static_cast<std::size_t>(input)];
                }
            }
            support += label;
            current = std::move(next);
        }
    }
    const cpp_int result = current[static_cast<std::size_t>(target)];
    cache.emplace(key, result);
    return result;
}

cpp_int outer_value(const Counts& counts, const Counts& signs) {
    OuterKey key{};
    for (std::size_t index = 0U; index < counts.size(); ++index) {
        key[index] = counts[index];
        key[counts.size() + index] = signs[index];
    }
    auto& cache = outer_cache();
    const auto found = cache.find(key);
    if (found != cache.end()) {
        return found->second;
    }
    const int target = total_label(counts) - deficit;
    cpp_int result = 0;
    Counts selected{};
    const auto visit = [&](const auto& self, int coordinate,
                           int selected_label) -> void {
        if (coordinate == label_count) {
            const cpp_int invariant = multiplicity(selected, 0);
            if (invariant == 0) {
                return;
            }
            Counts complement{};
            cpp_int weight = invariant;
            int parity = 0;
            for (int index = 0; index < label_count; ++index) {
                const std::size_t position = static_cast<std::size_t>(index);
                complement[position] = counts[position] - selected[position];
                weight *= binomial(counts[position], selected[position]);
                if (signs[position] < 0) {
                    parity += selected[position];
                }
            }
            weight *= multiplicity(complement, target);
            result += (parity & 1) != 0 ? -weight : weight;
            return;
        }
        const int label = coordinate + 1;
        const int maximum_selected = std::min(
            counts[static_cast<std::size_t>(coordinate)],
            (deficit - selected_label) / label
        );
        for (int count = 0; count <= maximum_selected; ++count) {
            selected[static_cast<std::size_t>(coordinate)] = count;
            self(self, coordinate + 1, selected_label + label * count);
        }
        selected[static_cast<std::size_t>(coordinate)] = 0;
    };
    visit(visit, 0, 0);
    cache.emplace(key, result);
    return result;
}

cpp_int finite_difference(
    const Counts& base,
    const Counts& degree,
    const Counts& signs
) {
    cpp_int result = 0;
    Counts offset{};
    int total_degree = 0;
    for (int value : degree) {
        total_degree += value;
    }
    const auto visit = [&](const auto& self, int coordinate,
                           int used_degree, long long coefficient) -> void {
        if (coordinate == label_count) {
            Counts point{};
            for (std::size_t index = 0U; index < point.size(); ++index) {
                point[index] = base[index] + offset[index];
            }
            const cpp_int term = coefficient * outer_value(point, signs);
            result += ((total_degree - used_degree) & 1) != 0
                ? -term : term;
            return;
        }
        const int maximum = degree[static_cast<std::size_t>(coordinate)];
        for (int value = 0; value <= maximum; ++value) {
            offset[static_cast<std::size_t>(coordinate)] = value;
            self(
                self, coordinate + 1, used_degree + value,
                coefficient * binomial(maximum, value)
            );
        }
    };
    visit(visit, 0, 0, 1);
    return result;
}

cpp_int cone_finite_difference(
    const Counts& base,
    const Directions& directions,
    const Counts& degree,
    const Counts& signs
) {
    cpp_int result = 0;
    Counts sample{};
    int total_degree = 0;
    for (int value : degree) {
        total_degree += value;
    }
    const auto visit = [&](const auto& self, int axis,
                           int used_degree, long long coefficient) -> void {
        if (axis == label_count) {
            Counts point = base;
            for (std::size_t direction = 0U;
                 direction < directions.size(); ++direction) {
                for (std::size_t coordinate = 0U;
                     coordinate < point.size(); ++coordinate) {
                    point[coordinate] += sample[direction]
                        * directions[direction][coordinate];
                }
            }
            const cpp_int term = coefficient * outer_value(point, signs);
            result += ((total_degree - used_degree) & 1) != 0
                ? -term : term;
            return;
        }
        const int maximum = degree[static_cast<std::size_t>(axis)];
        for (int value = 0; value <= maximum; ++value) {
            sample[static_cast<std::size_t>(axis)] = value;
            self(
                self, axis + 1, used_degree + value,
                coefficient * binomial(maximum, value)
            );
        }
    };
    visit(visit, 0, 0, 1);
    return result;
}

bool first_negative_cone_coefficient(
    const Counts& base,
    const Directions& directions,
    const Counts& signs,
    Counts& negative_degree,
    cpp_int& negative_value,
    long long& tested,
    unsigned int free_mask = (1U << label_count) - 1U
) {
    Counts degree{};
    for (int total_degree = 0; total_degree <= deficit; ++total_degree) {
        const auto visit = [&](const auto& self, int coordinate,
                               int remaining) -> bool {
            if (coordinate == label_count) {
                if (remaining != 0) {
                    return false;
                }
                const cpp_int coefficient = cone_finite_difference(
                    base, directions, degree, signs
                );
                ++tested;
                if (coefficient < 0) {
                    negative_degree = degree;
                    negative_value = coefficient;
                    return true;
                }
                return false;
            }
            const unsigned int bit = 1U << coordinate;
            if ((free_mask & bit) == 0U) {
                degree[static_cast<std::size_t>(coordinate)] = 0;
                return self(self, coordinate + 1, remaining);
            }
            for (int value = 0; value <= remaining; ++value) {
                degree[static_cast<std::size_t>(coordinate)] = value;
                if (self(self, coordinate + 1, remaining - value)) {
                    return true;
                }
            }
            return false;
        };
        if (visit(visit, 0, total_degree)) {
            return true;
        }
    }
    return false;
}

Counts signs_for_pattern(int pattern) {
    Counts signs{};
    signs[0U] = 1;
    for (int index = 1; index < label_count; ++index) {
        const int bit = 1 << (label_count - 1 - index);
        signs[static_cast<std::size_t>(index)]
            = (pattern & bit) != 0 ? -1 : 1;
    }
    return signs;
}

void print_counts(const Counts& counts) {
    for (std::size_t index = 0U; index < counts.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << counts[index];
    }
}

BivariateMonomial upper_bivariate_monomial(
    const Counts& base,
    const Counts& signs
) {
    Directions directions{};
    directions[0U] = Counts{1, 2, 0, 0, 0};
    directions[1U] = Counts{0, 1, 0, 0, 0};
    for (std::size_t coordinate = 2U; coordinate < directions.size();
         ++coordinate) {
        directions[coordinate][coordinate] = 1;
    }
    std::array<cpp_int, deficit + 1> factorial{};
    factorial[0U] = 1;
    for (int degree = 1; degree <= deficit; ++degree) {
        factorial[static_cast<std::size_t>(degree)]
            = degree * factorial[static_cast<std::size_t>(degree - 1)];
    }
    std::array<std::array<cpp_int, deficit + 1>, deficit + 1> stirling{};
    stirling[0U][0U] = 1;
    for (int degree = 1; degree <= deficit; ++degree) {
        for (int power = 1; power <= degree; ++power) {
            stirling[static_cast<std::size_t>(degree)]
                     [static_cast<std::size_t>(power)]
                = stirling[static_cast<std::size_t>(degree - 1)]
                          [static_cast<std::size_t>(power - 1)]
                - (degree - 1)
                    * stirling[static_cast<std::size_t>(degree - 1)]
                               [static_cast<std::size_t>(power)];
        }
    }
    BivariateMonomial monomial{};
    for (int first_degree = 0; first_degree <= deficit; ++first_degree) {
        for (int second_degree = 0;
             first_degree + second_degree <= deficit; ++second_degree) {
            const Counts degree{first_degree, second_degree, 0, 0, 0};
            const cpp_int newton = cone_finite_difference(
                base, directions, degree, signs
            );
            const cpp_int scale
                = factorial[static_cast<std::size_t>(deficit)]
                / (factorial[static_cast<std::size_t>(first_degree)]
                   * factorial[static_cast<std::size_t>(second_degree)]);
            for (int first_power = 0; first_power <= first_degree;
                 ++first_power) {
                for (int second_power = 0; second_power <= second_degree;
                     ++second_power) {
                    monomial[static_cast<std::size_t>(first_power)]
                            [static_cast<std::size_t>(second_power)]
                        += newton * scale
                        * stirling[static_cast<std::size_t>(first_degree)]
                                   [static_cast<std::size_t>(first_power)]
                        * stirling[static_cast<std::size_t>(second_degree)]
                                   [static_cast<std::size_t>(second_power)];
                }
            }
        }
    }
    return monomial;
}

void trim_polynomial(RationalPolynomial& polynomial) {
    while (!polynomial.empty() && polynomial.back() == 0) {
        polynomial.pop_back();
    }
}

RationalPolynomial derivative(const RationalPolynomial& polynomial) {
    RationalPolynomial result;
    if (polynomial.size() <= 1U) {
        return result;
    }
    result.resize(polynomial.size() - 1U);
    for (std::size_t degree = 1U; degree < polynomial.size(); ++degree) {
        result[degree - 1U] = static_cast<long long>(degree)
            * polynomial[degree];
    }
    trim_polynomial(result);
    return result;
}

RationalPolynomial polynomial_remainder(
    RationalPolynomial dividend,
    const RationalPolynomial& divisor
) {
    if (divisor.empty()) {
        throw std::runtime_error("zero polynomial divisor");
    }
    trim_polynomial(dividend);
    while (!dividend.empty() && dividend.size() >= divisor.size()) {
        const std::size_t shift = dividend.size() - divisor.size();
        const Rational factor = dividend.back() / divisor.back();
        for (std::size_t degree = 0U; degree < divisor.size(); ++degree) {
            dividend[degree + shift] -= factor * divisor[degree];
        }
        trim_polynomial(dividend);
    }
    return dividend;
}

int rational_sign(const Rational& value) {
    if (value > 0) {
        return 1;
    }
    if (value < 0) {
        return -1;
    }
    return 0;
}

int sturm_variations(
    const std::vector<RationalPolynomial>& sequence,
    bool at_positive_infinity
) {
    int previous = 0;
    int variations = 0;
    for (const RationalPolynomial& polynomial : sequence) {
        int sign = 0;
        if (at_positive_infinity) {
            sign = rational_sign(polynomial.back());
        } else {
            for (const Rational& coefficient : polynomial) {
                sign = rational_sign(coefficient);
                if (sign != 0) {
                    break;
                }
            }
        }
        if (sign != 0 && previous != 0 && sign != previous) {
            ++variations;
        }
        if (sign != 0) {
            previous = sign;
        }
    }
    return variations;
}

int positive_root_count(const std::vector<cpp_int>& coefficients) {
    RationalPolynomial polynomial;
    polynomial.reserve(coefficients.size());
    for (const cpp_int& coefficient : coefficients) {
        polynomial.emplace_back(coefficient);
    }
    trim_polynomial(polynomial);
    if (polynomial.empty() || polynomial[0U] == 0) {
        throw std::runtime_error("Sturm endpoint is a root");
    }
    if (polynomial.size() == 1U) {
        return 0;
    }
    std::vector<RationalPolynomial> sequence;
    sequence.push_back(polynomial);
    sequence.push_back(derivative(polynomial));
    while (!sequence.back().empty()) {
        RationalPolynomial remainder = polynomial_remainder(
            sequence[sequence.size() - 2U], sequence.back()
        );
        for (Rational& coefficient : remainder) {
            coefficient = -coefficient;
        }
        trim_polynomial(remainder);
        if (remainder.empty()) {
            break;
        }
        sequence.push_back(std::move(remainder));
    }
    return sturm_variations(sequence, false)
        - sturm_variations(sequence, true);
}

template <std::size_t Size>
int positive_root_count(const std::array<cpp_int, Size>& coefficients) {
    return positive_root_count(
        std::vector<cpp_int>(coefficients.begin(), coefficients.end())
    );
}

struct ParabolicSliceStats {
    int cutoff = 0;
    int maximum_positive_roots = 0;
    std::size_t negative_monomials = 0U;
};

cpp_int integer_power(int base, int exponent) {
    cpp_int result = 1;
    for (int index = 0; index < exponent; ++index) {
        result *= base;
    }
    return result;
}

cpp_int expected_parabolic_leading_coefficient(int first, int second) {
    if (first == 10 && second == 0) {
        return 42;
    }
    if (first == 8 && second == 1) {
        return -1260;
    }
    if (first == 6 && second == 2) {
        return 15120;
    }
    if (first == 4 && second == 3) {
        return -50400;
    }
    if (first == 2 && second == 4) {
        return 151200;
    }
    if (first == 0 && second == 5) {
        return 302400;
    }
    return 0;
}

bool certify_upper_bivariate_slice(
    const Counts& base,
    const Counts& signs,
    ParabolicSliceStats& stats
) {
    const BivariateMonomial monomial
        = upper_bivariate_monomial(base, signs);
    constexpr std::array<long long, 6> h_coefficients{
        1, -30, 360, -1200, 3600, 7200
    };
    for (int power = 3; power <= 5; ++power) {
        std::array<cpp_int, 6> difference{};
        for (std::size_t degree = 0U; degree < difference.size(); ++degree) {
            difference[degree] = h_coefficients[degree];
        }
        --difference[static_cast<std::size_t>(power)];
        if (positive_root_count(difference) != 0) {
            return false;
        }
    }
    struct NegativeTerm {
        cpp_int magnitude;
        int second_power = 0;
        int weight_gap = 0;
    };
    std::vector<NegativeTerm> negative_terms;
    for (int first = 0; first <= deficit; ++first) {
        for (int second = 0; second <= deficit; ++second) {
            const cpp_int& coefficient
                = monomial[static_cast<std::size_t>(first)]
                          [static_cast<std::size_t>(second)];
            const int weight = first + 2 * second;
            if (weight > deficit && coefficient != 0) {
                return false;
            }
            if (weight == deficit) {
                if (coefficient
                    != expected_parabolic_leading_coefficient(
                        first, second
                    )) {
                    return false;
                }
                continue;
            }
            if (coefficient < 0) {
                if (second > 5) {
                    return false;
                }
                negative_terms.push_back(
                    NegativeTerm{-coefficient, second, deficit - weight}
                );
            }
        }
    }
    stats.negative_monomials = negative_terms.size() + 2U;
    // Since K>=3, H>=1/6, H>=2z, and H>=45z^2.
    // Exact Sturm checks above additionally give H>=z^b for b=3,4,5.
    constexpr std::array<long long, 6> leading_lower_bounds{
        7, 84, 1890, 42, 42, 42
    };
    constexpr int maximum_cutoff = 131072;
    int cutoff = 1;
    for (;; cutoff *= 2) {
        Rational ratio_sum{0};
        for (const NegativeTerm& term : negative_terms) {
            const cpp_int denominator
                = leading_lower_bounds[
                    static_cast<std::size_t>(term.second_power)
                ] * integer_power(cutoff, term.weight_gap);
            ratio_sum += Rational{term.magnitude, denominator};
        }
        if (ratio_sum <= 1) {
            break;
        }
        if (cutoff >= maximum_cutoff) {
            return false;
        }
    }
    stats.cutoff = cutoff;
    for (int first = 0; first < cutoff; ++first) {
        std::array<cpp_int, 6> coefficients{};
        cpp_int power = 1;
        for (int first_power = 0; first_power <= deficit; ++first_power) {
            for (int second_power = 0; second_power <= 5; ++second_power) {
                coefficients[static_cast<std::size_t>(second_power)]
                    += monomial[static_cast<std::size_t>(first_power)]
                               [static_cast<std::size_t>(second_power)]
                    * power;
            }
            power *= first;
        }
        if (coefficients[0U] <= 0) {
            return false;
        }
        const int roots = positive_root_count(coefficients);
        stats.maximum_positive_roots = std::max(
            stats.maximum_positive_roots, roots
        );
        if (roots != 0) {
            return false;
        }
    }
    return true;
}

MultivariateMonomial transformed_monomial_polynomial(
    const Counts& base,
    const Directions& directions,
    const Counts& signs
) {
    std::array<cpp_int, deficit + 1> factorial{};
    factorial[0U] = 1;
    for (int degree = 1; degree <= deficit; ++degree) {
        factorial[static_cast<std::size_t>(degree)]
            = degree * factorial[static_cast<std::size_t>(degree - 1)];
    }
    std::array<std::array<cpp_int, deficit + 1>, deficit + 1> stirling{};
    stirling[0U][0U] = 1;
    for (int degree = 1; degree <= deficit; ++degree) {
        for (int power = 1; power <= degree; ++power) {
            stirling[static_cast<std::size_t>(degree)]
                     [static_cast<std::size_t>(power)]
                = stirling[static_cast<std::size_t>(degree - 1)]
                          [static_cast<std::size_t>(power - 1)]
                - (degree - 1)
                    * stirling[static_cast<std::size_t>(degree - 1)]
                               [static_cast<std::size_t>(power)];
        }
    }
    MultivariateMonomial monomial;
    Counts degree{};
    for (int total_degree = 0; total_degree <= deficit; ++total_degree) {
        const auto visit_degree = [&](const auto& self, int coordinate,
                                      int remaining) -> void {
            if (coordinate == label_count) {
                if (remaining != 0) {
                    return;
                }
                const cpp_int newton = cone_finite_difference(
                    base, directions, degree, signs
                );
                cpp_int scale = 1;
                for (int integer = 2; integer <= deficit; ++integer) {
                    scale *= integer;
                }
                for (int value : degree) {
                    for (int integer = 2; integer <= value; ++integer) {
                        scale /= integer;
                    }
                }
                Counts power{};
                const auto visit_power = [&](const auto& power_self,
                                             int axis) -> void {
                    if (axis == label_count) {
                        cpp_int term = newton * scale;
                        for (std::size_t index = 0U;
                             index < power.size(); ++index) {
                            term *= stirling[
                                static_cast<std::size_t>(degree[index])
                            ][static_cast<std::size_t>(power[index])];
                        }
                        monomial[power] += term;
                        return;
                    }
                    for (int value = 0;
                         value <= degree[static_cast<std::size_t>(axis)];
                         ++value) {
                        power[static_cast<std::size_t>(axis)] = value;
                        power_self(power_self, axis + 1);
                    }
                };
                visit_power(visit_power, 0);
                return;
            }
            for (int value = 0; value <= remaining; ++value) {
                degree[static_cast<std::size_t>(coordinate)] = value;
                self(self, coordinate + 1, remaining - value);
            }
        };
        visit_degree(visit_degree, 0, total_degree);
    }
    return monomial;
}

struct CylindricalStats {
    std::size_t polynomials = 0U;
    std::size_t nonzero_polynomials = 0U;
    int maximum_degree = 0;
    int maximum_positive_roots = 0;
    Counts failure_degree{};
};

bool certify_cylindrical_orthant(
    const Counts& base,
    const Directions& directions,
    const Counts& signs,
    std::size_t parameter_axis,
    CylindricalStats& stats
) {
    const MultivariateMonomial monomial
        = transformed_monomial_polynomial(base, directions, signs);
    std::array<std::array<cpp_int, deficit + 1>, deficit + 1> second{};
    second[0U][0U] = 1;
    for (int power = 1; power <= deficit; ++power) {
        for (int degree = 1; degree <= power; ++degree) {
            second[static_cast<std::size_t>(power)]
                  [static_cast<std::size_t>(degree)]
                = second[static_cast<std::size_t>(power - 1)]
                        [static_cast<std::size_t>(degree - 1)]
                + degree
                    * second[static_cast<std::size_t>(power - 1)]
                            [static_cast<std::size_t>(degree)];
        }
    }
    std::array<cpp_int, deficit + 1> factorial{};
    factorial[0U] = 1;
    for (int degree = 1; degree <= deficit; ++degree) {
        factorial[static_cast<std::size_t>(degree)]
            = degree * factorial[static_cast<std::size_t>(degree - 1)];
    }
    using ParameterPolynomial = std::array<cpp_int, deficit + 1>;
    std::map<Counts, ParameterPolynomial> coefficient_polynomials;
    for (const auto& [power, coefficient] : monomial) {
        if (coefficient == 0) {
            continue;
        }
        Counts newton_degree{};
        const auto visit = [&](const auto& self, int axis,
                               cpp_int factor) -> void {
            if (axis == label_count) {
                coefficient_polynomials[newton_degree]
                    [static_cast<std::size_t>(power[parameter_axis])]
                    += coefficient * factor;
                return;
            }
            const std::size_t position = static_cast<std::size_t>(axis);
            if (position == parameter_axis) {
                newton_degree[position] = 0;
                self(self, axis + 1, factor);
                return;
            }
            for (int degree = 0; degree <= power[position]; ++degree) {
                newton_degree[position] = degree;
                self(
                    self, axis + 1,
                    factor
                        * factorial[static_cast<std::size_t>(degree)]
                        * second[static_cast<std::size_t>(power[position])]
                                [static_cast<std::size_t>(degree)]
                );
            }
        };
        visit(visit, 0, cpp_int{1});
    }
    stats.polynomials += coefficient_polynomials.size();
    for (const auto& [degree, coefficients] : coefficient_polynomials) {
        int last = deficit;
        while (last >= 0
               && coefficients[static_cast<std::size_t>(last)] == 0) {
            --last;
        }
        if (last < 0) {
            continue;
        }
        ++stats.nonzero_polynomials;
        int first = 0;
        while (first <= last
               && coefficients[static_cast<std::size_t>(first)] == 0) {
            ++first;
        }
        if (coefficients[static_cast<std::size_t>(first)] < 0
            || coefficients[static_cast<std::size_t>(last)] < 0) {
            stats.failure_degree = degree;
            return false;
        }
        std::vector<cpp_int> reduced;
        for (int power = first; power <= last; ++power) {
            reduced.push_back(
                coefficients[static_cast<std::size_t>(power)]
            );
        }
        const int roots = positive_root_count(reduced);
        stats.maximum_positive_roots = std::max(
            stats.maximum_positive_roots, roots
        );
        stats.maximum_degree = std::max(
            stats.maximum_degree, last - first
        );
        if (roots != 0) {
            stats.failure_degree = degree;
            return false;
        }
    }
    return true;
}

bool strictly_positive_on_positive_axis(std::vector<cpp_int> coefficients) {
    while (!coefficients.empty() && coefficients.back() == 0) {
        coefficients.pop_back();
    }
    if (coefficients.empty()) {
        return false;
    }
    std::size_t first = 0U;
    while (first < coefficients.size() && coefficients[first] == 0) {
        ++first;
    }
    if (first == coefficients.size() || coefficients[first] < 0
        || coefficients.back() < 0) {
        return false;
    }
    coefficients.erase(
        coefficients.begin(),
        coefficients.begin() + static_cast<std::ptrdiff_t>(first)
    );
    return positive_root_count(coefficients) == 0;
}

struct WeightedBivariateStats {
    std::size_t polynomials = 0U;
    std::size_t direct_polynomials = 0U;
    std::size_t sturm_polynomials = 0U;
    std::size_t sturm_slices = 0U;
    std::size_t finite_newton_slices = 0U;
    int maximum_cutoff = 0;
    int maximum_weight = 0;
    Counts failure_degree{};
    std::string failure_reason;
};

bool certify_weighted_bivariate_polynomial(
    const BivariateMonomial& polynomial,
    int second_weight,
    WeightedBivariateStats& stats
) {
    ++stats.polynomials;
    bool has_negative = false;
    bool nonzero = false;
    int maximum_weight = -1;
    for (int first = 0; first <= deficit; ++first) {
        for (int second = 0; second <= deficit; ++second) {
            const cpp_int& coefficient
                = polynomial[static_cast<std::size_t>(first)]
                            [static_cast<std::size_t>(second)];
            if (coefficient == 0) {
                continue;
            }
            nonzero = true;
            has_negative = has_negative || coefficient < 0;
            maximum_weight = std::max(
                maximum_weight, first + second_weight * second
            );
        }
    }
    if (!nonzero || !has_negative) {
        ++stats.direct_polynomials;
        return true;
    }
    stats.maximum_weight = std::max(stats.maximum_weight, maximum_weight);
    std::vector<cpp_int> leading(static_cast<std::size_t>(deficit + 1));
    for (int first = 0; first <= deficit; ++first) {
        for (int second = 0; second <= deficit; ++second) {
            if (first + second_weight * second == maximum_weight) {
                leading[static_cast<std::size_t>(second)]
                    += polynomial[static_cast<std::size_t>(first)]
                                 [static_cast<std::size_t>(second)];
            }
        }
    }
    if (!strictly_positive_on_positive_axis(leading)) {
        stats.failure_reason = "leading";
        return false;
    }
    struct NegativeTerm {
        cpp_int magnitude;
        int gap = 0;
        int bound_power = 0;
    };
    std::map<int, int> bound_powers;
    std::vector<NegativeTerm> negative_terms;
    for (int first = 0; first <= deficit; ++first) {
        for (int second = 0; second <= deficit; ++second) {
            const cpp_int& coefficient
                = polynomial[static_cast<std::size_t>(first)]
                            [static_cast<std::size_t>(second)];
            const int weight = first + second_weight * second;
            if (coefficient >= 0 || weight == maximum_weight) {
                continue;
            }
            auto found = bound_powers.find(second);
            if (found == bound_powers.end()) {
                int bound_power = 30;
                bool bound_found = false;
                for (; bound_power >= -60; --bound_power) {
                    std::vector<cpp_int> difference = leading;
                    if (bound_power < 0) {
                        const cpp_int scale = cpp_int{1} << -bound_power;
                        for (cpp_int& value : difference) {
                            value *= scale;
                        }
                    }
                    if (difference.size()
                        <= static_cast<std::size_t>(second)) {
                        difference.resize(
                            static_cast<std::size_t>(second + 1)
                        );
                    }
                    difference[static_cast<std::size_t>(second)]
                        -= bound_power >= 0
                        ? cpp_int{1} << bound_power : cpp_int{1};
                    if (strictly_positive_on_positive_axis(difference)) {
                        bound_found = true;
                        break;
                    }
                }
                if (!bound_found) {
                    stats.failure_reason = "bound_second_"
                        + std::to_string(second);
                    return false;
                }
                found = bound_powers.emplace(second, bound_power).first;
            }
            negative_terms.push_back(
                NegativeTerm{
                    -coefficient, maximum_weight - weight, found->second
                }
            );
        }
    }
    constexpr int maximum_cutoff = 262144;
    int cutoff = 1;
    for (;; cutoff *= 2) {
        Rational ratio_sum{0};
        for (const NegativeTerm& term : negative_terms) {
            cpp_int numerator = term.magnitude;
            cpp_int denominator = integer_power(cutoff, term.gap);
            if (term.bound_power >= 0) {
                denominator *= cpp_int{1} << term.bound_power;
            } else {
                numerator *= cpp_int{1} << -term.bound_power;
            }
            ratio_sum += Rational{numerator, denominator};
        }
        if (ratio_sum <= 1) {
            break;
        }
        if (cutoff >= maximum_cutoff) {
            stats.failure_reason = "cutoff";
            return false;
        }
    }
    stats.maximum_cutoff = std::max(stats.maximum_cutoff, cutoff);
    ++stats.sturm_polynomials;
    std::array<std::array<cpp_int, deficit + 1>, deficit + 1> second{};
    second[0U][0U] = 1;
    for (int power_index = 1; power_index <= deficit; ++power_index) {
        for (int degree = 1; degree <= power_index; ++degree) {
            second[static_cast<std::size_t>(power_index)]
                  [static_cast<std::size_t>(degree)]
                = second[static_cast<std::size_t>(power_index - 1)]
                        [static_cast<std::size_t>(degree - 1)]
                + degree
                    * second[static_cast<std::size_t>(power_index - 1)]
                            [static_cast<std::size_t>(degree)];
        }
    }
    std::array<cpp_int, deficit + 1> factorial{};
    factorial[0U] = 1;
    for (int degree = 1; degree <= deficit; ++degree) {
        factorial[static_cast<std::size_t>(degree)]
            = degree * factorial[static_cast<std::size_t>(degree - 1)];
    }
    for (int first = 0; first < cutoff; ++first) {
        std::vector<cpp_int> univariate(
            static_cast<std::size_t>(deficit + 1)
        );
        cpp_int power = 1;
        for (int first_power = 0; first_power <= deficit; ++first_power) {
            for (int second_power = 0; second_power <= deficit;
                 ++second_power) {
                univariate[static_cast<std::size_t>(second_power)]
                    += polynomial[static_cast<std::size_t>(first_power)]
                                 [static_cast<std::size_t>(second_power)]
                    * power;
            }
            power *= first;
        }
        bool zero = true;
        for (const cpp_int& coefficient : univariate) {
            zero = zero && coefficient == 0;
        }
        if (!zero) {
            bool newton_nonnegative = true;
            for (int degree = 0; degree <= deficit; ++degree) {
                cpp_int newton = 0;
                for (int power_index = degree; power_index <= deficit;
                     ++power_index) {
                    newton += univariate[
                        static_cast<std::size_t>(power_index)
                    ] * factorial[static_cast<std::size_t>(degree)]
                    * second[static_cast<std::size_t>(power_index)]
                            [static_cast<std::size_t>(degree)];
                }
                if (newton < 0) {
                    newton_nonnegative = false;
                    break;
                }
            }
            if (newton_nonnegative) {
                ++stats.finite_newton_slices;
                continue;
            }
            if (!strictly_positive_on_positive_axis(univariate)) {
                stats.failure_reason = "finite_first_"
                    + std::to_string(first);
                return false;
            }
        }
        ++stats.sturm_slices;
    }
    return true;
}

bool certify_bivariate_newton_orthant(
    const Counts& base,
    const Directions& directions,
    const Counts& signs,
    int second_weight,
    WeightedBivariateStats& stats
) {
    const MultivariateMonomial monomial
        = transformed_monomial_polynomial(base, directions, signs);
    std::array<std::array<cpp_int, deficit + 1>, deficit + 1> second{};
    second[0U][0U] = 1;
    for (int power = 1; power <= deficit; ++power) {
        for (int degree = 1; degree <= power; ++degree) {
            second[static_cast<std::size_t>(power)]
                  [static_cast<std::size_t>(degree)]
                = second[static_cast<std::size_t>(power - 1)]
                        [static_cast<std::size_t>(degree - 1)]
                + degree
                    * second[static_cast<std::size_t>(power - 1)]
                            [static_cast<std::size_t>(degree)];
        }
    }
    std::array<cpp_int, deficit + 1> factorial{};
    factorial[0U] = 1;
    for (int degree = 1; degree <= deficit; ++degree) {
        factorial[static_cast<std::size_t>(degree)]
            = degree * factorial[static_cast<std::size_t>(degree - 1)];
    }
    std::map<Counts, BivariateMonomial> polynomials;
    for (const auto& [power, coefficient] : monomial) {
        if (coefficient == 0) {
            continue;
        }
        Counts newton_degree{};
        const auto visit = [&](const auto& self, int axis,
                               cpp_int factor) -> void {
            if (axis == label_count) {
                polynomials[newton_degree]
                    [static_cast<std::size_t>(power[0U])]
                    [static_cast<std::size_t>(power[1U])]
                    += coefficient * factor;
                return;
            }
            if (axis < 2) {
                self(self, axis + 1, factor);
                return;
            }
            const std::size_t position = static_cast<std::size_t>(axis);
            for (int degree = 0; degree <= power[position]; ++degree) {
                newton_degree[position] = degree;
                self(
                    self, axis + 1,
                    factor
                        * factorial[static_cast<std::size_t>(degree)]
                        * second[static_cast<std::size_t>(power[position])]
                                [static_cast<std::size_t>(degree)]
                );
            }
        };
        visit(visit, 0, cpp_int{1});
    }
    for (const auto& [degree, polynomial] : polynomials) {
        if (!certify_weighted_bivariate_polynomial(
                polynomial, second_weight, stats
            )) {
            stats.failure_degree = degree;
            if (stats.failure_reason.empty()) {
                stats.failure_reason = "bivariate_polynomial";
            }
            return false;
        }
    }
    return true;
}

std::vector<Counts> minimal_bases() {
    std::vector<Counts> result;
    // If an admissible point had total label at least 21, removal of any
    // one present factor would leave total label at least 16, still at least
    // deficit plus the largest possible remaining label five.  Thus every
    // product-order minimal point has total label at most 20.
    for (int first = 0; first <= 20; ++first) {
        for (int second = 0; second <= 10; ++second) {
            for (int third = 0; third <= 6; ++third) {
                for (int fourth = 0; fourth <= 5; ++fourth) {
                    for (int fifth = 0; fifth <= 4; ++fifth) {
                        const Counts counts{
                            first, second, third, fourth, fifth
                        };
                        if (total_label(counts) > 20 || !admissible(counts)) {
                            continue;
                        }
                        bool minimal = true;
                        for (std::size_t coordinate = 0U;
                             coordinate < counts.size(); ++coordinate) {
                            if (counts[coordinate] == 0) {
                                continue;
                            }
                            Counts predecessor = counts;
                            --predecessor[coordinate];
                            if (admissible(predecessor)) {
                                minimal = false;
                                break;
                            }
                        }
                        if (minimal) {
                            result.push_back(counts);
                        }
                    }
                }
            }
        }
    }
    return result;
}

Counts apply_directions(
    const Counts& transformed,
    const Directions& directions
) {
    Counts counts{};
    for (std::size_t direction = 0U; direction < directions.size();
         ++direction) {
        for (std::size_t coordinate = 0U; coordinate < counts.size();
             ++coordinate) {
            counts[coordinate] += transformed[direction]
                * directions[direction][coordinate];
        }
    }
    return counts;
}

std::vector<Counts> transformed_minimal_bases(
    const Directions& directions
) {
    Counts weights{};
    int maximum_weight = 0;
    for (std::size_t direction = 0U; direction < directions.size();
         ++direction) {
        weights[direction] = total_label(directions[direction]);
        if (weights[direction] <= 0) {
            throw std::runtime_error("nonpositive chamber direction");
        }
        maximum_weight = std::max(maximum_weight, weights[direction]);
    }
    // If S exceeds deficit+maximum_label+maximum_direction_weight, removal
    // of any present direction remains admissible.  The inclusive bound is
    // deliberately one unit conservative.
    const int weight_bound = 15 + maximum_weight;
    std::vector<Counts> result;
    Counts transformed{};
    const auto visit = [&](const auto& self, int coordinate,
                           int used_weight) -> void {
        if (coordinate == label_count) {
            const Counts counts = apply_directions(transformed, directions);
            if (!admissible(counts)) {
                return;
            }
            bool minimal = true;
            for (std::size_t axis = 0U; axis < transformed.size(); ++axis) {
                if (transformed[axis] == 0) {
                    continue;
                }
                Counts predecessor = transformed;
                --predecessor[axis];
                if (admissible(apply_directions(
                        predecessor, directions
                    ))) {
                    minimal = false;
                    break;
                }
            }
            if (minimal) {
                result.push_back(transformed);
            }
            return;
        }
        const std::size_t position = static_cast<std::size_t>(coordinate);
        const int maximum = (weight_bound - used_weight) / weights[position];
        for (int value = 0; value <= maximum; ++value) {
            transformed[position] = value;
            self(
                self, coordinate + 1,
                used_weight + value * weights[position]
            );
        }
        transformed[position] = 0;
    };
    visit(visit, 0, 0);
    return result;
}

bool first_negative_newton_coefficient(
    const Counts& base,
    const Counts& signs,
    Counts& negative_degree,
    cpp_int& negative_value,
    long long& tested,
    unsigned int free_mask = (1U << label_count) - 1U
) {
    Counts degree{};
    for (int total_degree = 0; total_degree <= deficit; ++total_degree) {
        const auto visit = [&](const auto& self, int coordinate,
                               int remaining) -> bool {
            if (coordinate == label_count) {
                if (remaining != 0) {
                    return false;
                }
                const cpp_int coefficient = finite_difference(
                    base, degree, signs
                );
                ++tested;
                if (coefficient < 0) {
                    negative_degree = degree;
                    negative_value = coefficient;
                    return true;
                }
                return false;
            }
            const unsigned int bit = 1U << coordinate;
            if ((free_mask & bit) == 0U) {
                degree[static_cast<std::size_t>(coordinate)] = 0;
                return self(self, coordinate + 1, remaining);
            }
            for (int value = 0; value <= remaining; ++value) {
                degree[static_cast<std::size_t>(coordinate)] = value;
                if (self(self, coordinate + 1, remaining - value)) {
                    return true;
                }
            }
            return false;
        };
        if (visit(visit, 0, total_degree)) {
            return true;
        }
    }
    return false;
}

struct RecursiveStats {
    std::size_t nodes = 0U;
    std::size_t leaves = 0U;
    std::size_t splits = 0U;
    std::size_t maximum_depth = 0U;
    long long coefficients = 0;
    Counts maximum_base{};
    Counts failure_degree{};
};

struct RecursiveConeStats {
    std::size_t nodes = 0U;
    std::size_t leaves = 0U;
    std::size_t splits = 0U;
    std::size_t maximum_depth = 0U;
    long long coefficients = 0;
    Counts maximum_transformed_base{};
    Counts maximum_base{};
    Counts failure_degree{};
};

bool recursive_cone_certificate(
    const Counts& transformed_base,
    unsigned int free_mask,
    const Directions& directions,
    const Counts& signs,
    std::size_t depth,
    RecursiveConeStats& stats,
    std::map<RecursiveKey, bool>& memo
) {
    RecursiveKey key{};
    for (std::size_t coordinate = 0U;
         coordinate < transformed_base.size(); ++coordinate) {
        key[coordinate] = transformed_base[coordinate];
    }
    key[transformed_base.size()] = static_cast<int>(free_mask);
    const auto found = memo.find(key);
    if (found != memo.end()) {
        return found->second;
    }
    ++stats.nodes;
    stats.maximum_depth = std::max(stats.maximum_depth, depth);
    const Counts base = apply_directions(transformed_base, directions);
    for (std::size_t coordinate = 0U; coordinate < base.size();
         ++coordinate) {
        stats.maximum_transformed_base[coordinate] = std::max(
            stats.maximum_transformed_base[coordinate],
            transformed_base[coordinate]
        );
        stats.maximum_base[coordinate] = std::max(
            stats.maximum_base[coordinate], base[coordinate]
        );
    }
    Counts negative_degree{};
    cpp_int negative_value = 0;
    if (!first_negative_cone_coefficient(
            base, directions, signs, negative_degree, negative_value,
            stats.coefficients, free_mask
        )) {
        ++stats.leaves;
        memo.emplace(key, true);
        return true;
    }
    if (free_mask == 0U || depth >= 50U || stats.nodes >= 10000U) {
        stats.failure_degree = negative_degree;
        memo.emplace(key, false);
        return false;
    }
    std::array<std::size_t, label_count> candidates{};
    for (std::size_t index = 0U; index < candidates.size(); ++index) {
        candidates[index] = index;
    }
    std::stable_sort(
        candidates.begin(), candidates.end(),
        [&negative_degree](std::size_t first, std::size_t second) {
            return negative_degree[first] > negative_degree[second];
        }
    );
    for (std::size_t axis : candidates) {
        const unsigned int bit = 1U << axis;
        if ((free_mask & bit) == 0U) {
            continue;
        }
        if (!recursive_cone_certificate(
                transformed_base, free_mask & ~bit, directions, signs,
                depth + 1U, stats, memo
            )) {
            continue;
        }
        Counts shifted = transformed_base;
        ++shifted[axis];
        if (recursive_cone_certificate(
                shifted, free_mask, directions, signs, depth + 1U, stats,
                memo
            )) {
            ++stats.splits;
            memo[key] = true;
            return true;
        }
    }
    stats.failure_degree = negative_degree;
    memo[key] = false;
    return false;
}

struct RecursiveHybridStats {
    std::size_t nodes = 0U;
    std::size_t leaves = 0U;
    std::size_t splits = 0U;
    std::size_t parabolic_checks = 0U;
    std::size_t maximum_depth = 0U;
    long long coefficients = 0;
    int maximum_parabolic_cutoff = 0;
    Counts maximum_transformed_base{};
    Counts maximum_base{};
    Counts failure_degree{};
    std::string failure_reason;
};

bool first_negative_high_cone_coefficient(
    const Counts& base,
    const Directions& directions,
    const Counts& signs,
    Counts& negative_degree,
    cpp_int& negative_value,
    long long& tested,
    unsigned int free_mask
) {
    Counts degree{};
    for (int total_degree = 0; total_degree <= deficit; ++total_degree) {
        const auto visit = [&](const auto& self, int coordinate,
                               int remaining) -> bool {
            if (coordinate == label_count) {
                if (remaining != 0
                    || degree[2U] + degree[3U] + degree[4U] == 0) {
                    return false;
                }
                const cpp_int coefficient = cone_finite_difference(
                    base, directions, degree, signs
                );
                ++tested;
                if (coefficient < 0) {
                    negative_degree = degree;
                    negative_value = coefficient;
                    return true;
                }
                return false;
            }
            const unsigned int bit = 1U << coordinate;
            if ((free_mask & bit) == 0U) {
                degree[static_cast<std::size_t>(coordinate)] = 0;
                return self(self, coordinate + 1, remaining);
            }
            for (int value = 0; value <= remaining; ++value) {
                degree[static_cast<std::size_t>(coordinate)] = value;
                if (self(self, coordinate + 1, remaining - value)) {
                    return true;
                }
            }
            return false;
        };
        if (visit(visit, 0, total_degree)) {
            return true;
        }
    }
    return false;
}

bool recursive_upper_hybrid_certificate(
    const Counts& transformed_base,
    unsigned int free_mask,
    const Directions& directions,
    const Counts& signs,
    std::size_t depth,
    RecursiveHybridStats& stats,
    std::map<RecursiveKey, bool>& memo
) {
    RecursiveKey key{};
    for (std::size_t coordinate = 0U;
         coordinate < transformed_base.size(); ++coordinate) {
        key[coordinate] = transformed_base[coordinate];
    }
    key[transformed_base.size()] = static_cast<int>(free_mask);
    const auto found = memo.find(key);
    if (found != memo.end()) {
        return found->second;
    }
    ++stats.nodes;
    stats.maximum_depth = std::max(stats.maximum_depth, depth);
    const Counts base = apply_directions(transformed_base, directions);
    for (std::size_t coordinate = 0U; coordinate < base.size();
         ++coordinate) {
        stats.maximum_transformed_base[coordinate] = std::max(
            stats.maximum_transformed_base[coordinate],
            transformed_base[coordinate]
        );
        stats.maximum_base[coordinate] = std::max(
            stats.maximum_base[coordinate], base[coordinate]
        );
    }
    Counts negative_degree{};
    cpp_int negative_value = 0;
    if (!first_negative_high_cone_coefficient(
            base, directions, signs, negative_degree, negative_value,
            stats.coefficients, free_mask
        )) {
        ++stats.parabolic_checks;
        ParabolicSliceStats slice_stats;
        if (certify_upper_bivariate_slice(base, signs, slice_stats)) {
            ++stats.leaves;
            stats.maximum_parabolic_cutoff = std::max(
                stats.maximum_parabolic_cutoff, slice_stats.cutoff
            );
            memo.emplace(key, true);
            return true;
        }
        stats.failure_reason = "parabolic_slice";
        memo.emplace(key, false);
        return false;
    }
    if (depth >= 50U || stats.nodes >= 10000U) {
        stats.failure_degree = negative_degree;
        stats.failure_reason = depth >= 50U ? "depth_cap" : "node_cap";
        memo.emplace(key, false);
        return false;
    }
    std::array<std::size_t, 3> candidates{2U, 3U, 4U};
    std::stable_sort(
        candidates.begin(), candidates.end(),
        [&negative_degree](std::size_t first, std::size_t second) {
            return negative_degree[first] > negative_degree[second];
        }
    );
    for (std::size_t axis : candidates) {
        const unsigned int bit = 1U << axis;
        if ((free_mask & bit) == 0U) {
            continue;
        }
        if (!recursive_upper_hybrid_certificate(
                transformed_base, free_mask & ~bit, directions, signs,
                depth + 1U, stats, memo
            )) {
            continue;
        }
        Counts shifted = transformed_base;
        ++shifted[axis];
        if (recursive_upper_hybrid_certificate(
                shifted, free_mask, directions, signs, depth + 1U, stats,
                memo
            )) {
            ++stats.splits;
            memo[key] = true;
            return true;
        }
    }
    stats.failure_degree = negative_degree;
    if (stats.failure_reason.empty()) {
        stats.failure_reason = "no_high_axis_split";
    }
    memo[key] = false;
    return false;
}

bool recursive_newton_certificate(
    const Counts& base,
    unsigned int free_mask,
    const Counts& signs,
    std::size_t depth,
    RecursiveStats& stats,
    std::map<RecursiveKey, bool>& memo
) {
    RecursiveKey key{};
    for (std::size_t coordinate = 0U; coordinate < base.size(); ++coordinate) {
        key[coordinate] = base[coordinate];
    }
    key[base.size()] = static_cast<int>(free_mask);
    const auto found = memo.find(key);
    if (found != memo.end()) {
        return found->second;
    }
    ++stats.nodes;
    stats.maximum_depth = std::max(stats.maximum_depth, depth);
    for (std::size_t coordinate = 0U; coordinate < base.size(); ++coordinate) {
        stats.maximum_base[coordinate] = std::max(
            stats.maximum_base[coordinate], base[coordinate]
        );
    }
    Counts negative_degree{};
    cpp_int negative_value = 0;
    if (!first_negative_newton_coefficient(
            base, signs, negative_degree, negative_value,
            stats.coefficients, free_mask
        )) {
        ++stats.leaves;
        memo.emplace(key, true);
        return true;
    }
    if (free_mask == 0U || depth >= 50U || stats.nodes >= 10000U) {
        stats.failure_degree = negative_degree;
        memo.emplace(key, false);
        return false;
    }
    std::array<std::size_t, label_count> candidates{};
    for (std::size_t index = 0U; index < candidates.size(); ++index) {
        candidates[index] = index;
    }
    std::stable_sort(
        candidates.begin(), candidates.end(),
        [&negative_degree](std::size_t first, std::size_t second) {
            return negative_degree[first] > negative_degree[second];
        }
    );
    for (std::size_t coordinate : candidates) {
        const unsigned int bit = 1U << coordinate;
        if ((free_mask & bit) == 0U) {
            continue;
        }
        if (!recursive_newton_certificate(
                base, free_mask & ~bit, signs, depth + 1U, stats, memo
            )) {
            continue;
        }
        Counts shifted = base;
        ++shifted[coordinate];
        if (recursive_newton_certificate(
                shifted, free_mask, signs, depth + 1U, stats, memo
            )) {
            ++stats.splits;
            memo[key] = true;
            return true;
        }
    }
    stats.failure_degree = negative_degree;
    memo[key] = false;
    return false;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 4) {
            throw std::runtime_error(
                "usage: analyze_su2_outer_ten MAXIMUM_TOTAL_COUNT "
                "[search|minimal-bases|certificate|recursive|probe|"
                "fan-probe|slope-probe|slope-global|chamber-negatives|"
                "chamber-slices|chamber-recursive|parabolic-tail|"
                "upper-monomial|parabolic-slice-certificate|"
                "upper-slice-all|upper-hybrid|parabolic-leading-full|"
                "upper-fixed-slices|cylindrical-sturm|bivariate-sturm] "
                "[PATTERN]"
            );
        }
        const int maximum_total_count = std::stoi(argv[1]);
        if (maximum_total_count < 1 || maximum_total_count > 40) {
            throw std::runtime_error("invalid maximum total count");
        }
        const std::string mode = argc >= 3 ? argv[2] : "search";
        if (mode != "search" && mode != "minimal-bases"
            && mode != "certificate" && mode != "recursive"
            && mode != "probe" && mode != "fan-probe"
            && mode != "slope-probe" && mode != "slope-global"
            && mode != "chamber-negatives" && mode != "chamber-slices"
            && mode != "chamber-recursive" && mode != "parabolic-tail"
            && mode != "upper-monomial"
            && mode != "parabolic-slice-certificate"
            && mode != "upper-slice-all" && mode != "upper-hybrid"
            && mode != "parabolic-leading-full"
            && mode != "upper-fixed-slices"
            && mode != "cylindrical-sturm" && mode != "bivariate-sturm") {
            throw std::runtime_error("invalid mode");
        }
        const int requested_pattern = argc == 4 ? std::stoi(argv[3]) : -1;
        if (argc == 4
            && ((mode != "certificate" && mode != "recursive"
                    && mode != "probe" && mode != "fan-probe"
                    && mode != "slope-probe" && mode != "slope-global"
                    && mode != "chamber-negatives"
                    && mode != "chamber-slices"
                    && mode != "chamber-recursive"
                    && mode != "parabolic-tail"
                    && mode != "upper-monomial"
                    && mode != "parabolic-slice-certificate"
                    && mode != "upper-slice-all"
                    && mode != "upper-hybrid"
                    && mode != "parabolic-leading-full"
                    && mode != "upper-fixed-slices"
                    && mode != "cylindrical-sturm"
                    && mode != "bivariate-sturm")
                || requested_pattern < 0 || requested_pattern >= 16)) {
            throw std::runtime_error("invalid certificate pattern");
        }
        const std::vector<Counts> bases = minimal_bases();
        std::cout << "SU2_OUTER_TEN_MINIMAL_BASES count=" << bases.size()
                  << '\n';
        if (mode == "minimal-bases") {
            for (const Counts& base : bases) {
                std::cout << "base=";
                print_counts(base);
                std::cout << " total_label=" << total_label(base) << '\n';
            }
            return EXIT_SUCCESS;
        }

        if (mode == "probe") {
            if (requested_pattern < 0) {
                throw std::runtime_error("probe mode requires a pattern");
            }
            const Counts base = bases.front();
            const Counts signs = signs_for_pattern(requested_pattern);
            Counts degree{};
            std::size_t negatives = 0U;
            long long tested = 0;
            for (int total_degree = 0; total_degree <= deficit;
                 ++total_degree) {
                const auto visit = [&](const auto& self, int coordinate,
                                       int remaining) -> void {
                    if (coordinate == label_count) {
                        if (remaining != 0) {
                            return;
                        }
                        const cpp_int coefficient = finite_difference(
                            base, degree, signs
                        );
                        ++tested;
                        if (coefficient < 0) {
                            ++negatives;
                            std::cout
                                << "SU2_OUTER_TEN_NEGATIVE pattern="
                                << requested_pattern << " base=";
                            print_counts(base);
                            std::cout << " degree=";
                            print_counts(degree);
                            std::cout << " value=" << coefficient << '\n';
                        }
                        return;
                    }
                    for (int value = 0; value <= remaining; ++value) {
                        degree[static_cast<std::size_t>(coordinate)] = value;
                        self(self, coordinate + 1, remaining - value);
                    }
                };
                visit(visit, 0, total_degree);
            }
            std::cout << "SU2_OUTER_TEN_NEGATIVE_SUMMARY pattern="
                      << requested_pattern << " base=";
            print_counts(base);
            std::cout << " negatives=" << negatives
                      << " coefficients=" << tested << '\n';
            return EXIT_SUCCESS;
        }

        if (mode == "fan-probe") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "fan-probe mode requires a pattern"
                );
            }
            const Counts base = bases.front();
            const Counts signs = signs_for_pattern(requested_pattern);
            for (std::size_t first = 0U; first < base.size(); ++first) {
                for (std::size_t second = first + 1U;
                     second < base.size(); ++second) {
                    Directions lower{};
                    Directions upper{};
                    for (std::size_t coordinate = 0U;
                         coordinate < base.size(); ++coordinate) {
                        lower[coordinate][coordinate] = 1;
                        upper[coordinate][coordinate] = 1;
                    }
                    ++lower[second][first];
                    ++upper[first][second];
                    for (int cone = 0; cone < 2; ++cone) {
                        const Directions& directions
                            = cone == 0 ? lower : upper;
                        Counts failure{};
                        cpp_int value = 0;
                        long long tested = 0;
                        const bool fails = first_negative_cone_coefficient(
                            base, directions, signs, failure, value, tested
                        );
                        std::cout << "SU2_OUTER_TEN_FAN_"
                                  << (fails ? "FAIL" : "PASS")
                                  << " pattern=" << requested_pattern
                                  << " pair=" << first << ',' << second
                                  << " cone=" << (cone == 0 ? "lower" : "upper")
                                  << " tested=" << tested;
                        if (fails) {
                            std::cout << " degree=";
                            print_counts(failure);
                            std::cout << " value=" << value;
                        }
                        std::cout << '\n';
                    }
                }
            }
            return EXIT_SUCCESS;
        }

        if (mode == "slope-probe") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "slope-probe mode requires a pattern"
                );
            }
            constexpr int maximum_slope = 16;
            std::vector<std::array<int, 2>> rays;
            rays.push_back({1, 0});
            for (int slope = 1; slope <= maximum_slope; ++slope) {
                rays.push_back({1, slope});
            }
            rays.push_back({0, 1});
            const Counts base = bases.front();
            const Counts signs = signs_for_pattern(requested_pattern);
            for (std::size_t ray = 0U; ray + 1U < rays.size(); ++ray) {
                Directions directions{};
                for (std::size_t coordinate = 2U;
                     coordinate < base.size(); ++coordinate) {
                    directions[coordinate][coordinate] = 1;
                }
                directions[0U][0U] = rays[ray][0U];
                directions[0U][1U] = rays[ray][1U];
                directions[1U][0U] = rays[ray + 1U][0U];
                directions[1U][1U] = rays[ray + 1U][1U];
                Counts failure{};
                cpp_int value = 0;
                long long tested = 0;
                const bool fails = first_negative_cone_coefficient(
                    base, directions, signs, failure, value, tested
                );
                std::cout << "SU2_OUTER_TEN_SLOPE_"
                          << (fails ? "FAIL" : "PASS")
                          << " pattern=" << requested_pattern << " rays="
                          << rays[ray][0U] << ',' << rays[ray][1U] << ':'
                          << rays[ray + 1U][0U] << ','
                          << rays[ray + 1U][1U] << " tested=" << tested;
                if (fails) {
                    std::cout << " degree=";
                    print_counts(failure);
                    std::cout << " value=" << value;
                }
                std::cout << '\n';
            }
            return EXIT_SUCCESS;
        }

        if (mode == "slope-global") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "slope-global mode requires a pattern"
                );
            }
            constexpr int maximum_slope = 16;
            std::vector<std::array<int, 2>> rays;
            rays.push_back({1, 0});
            for (int slope = 1; slope <= maximum_slope; ++slope) {
                rays.push_back({1, slope});
            }
            rays.push_back({0, 1});
            const Counts signs = signs_for_pattern(requested_pattern);
            std::size_t total_bases = 0U;
            long long total_coefficients = 0;
            for (std::size_t ray = 0U; ray + 1U < rays.size(); ++ray) {
                Directions directions{};
                for (std::size_t coordinate = 2U;
                     coordinate < directions.size(); ++coordinate) {
                    directions[coordinate][coordinate] = 1;
                }
                directions[0U][0U] = rays[ray][0U];
                directions[0U][1U] = rays[ray][1U];
                directions[1U][0U] = rays[ray + 1U][0U];
                directions[1U][1U] = rays[ray + 1U][1U];
                const std::vector<Counts> transformed_bases
                    = transformed_minimal_bases(directions);
                total_bases += transformed_bases.size();
                bool cone_passes = true;
                for (const Counts& transformed_base : transformed_bases) {
                    const Counts base = apply_directions(
                        transformed_base, directions
                    );
                    Counts failure{};
                    cpp_int value = 0;
                    long long tested = 0;
                    if (first_negative_cone_coefficient(
                            base, directions, signs, failure, value, tested
                        )) {
                        std::cout << "SU2_OUTER_TEN_SLOPE_GLOBAL_FAIL pattern="
                                  << requested_pattern << " rays="
                                  << rays[ray][0U] << ',' << rays[ray][1U]
                                  << ':' << rays[ray + 1U][0U] << ','
                                  << rays[ray + 1U][1U]
                                  << " transformed_base=";
                        print_counts(transformed_base);
                        std::cout << " base=";
                        print_counts(base);
                        std::cout << " degree=";
                        print_counts(failure);
                        std::cout << " value=" << value << '\n';
                        cone_passes = false;
                        break;
                    }
                    total_coefficients += tested;
                }
                if (!cone_passes) {
                    return EXIT_SUCCESS;
                }
                std::cout << "SU2_OUTER_TEN_SLOPE_GLOBAL_CONE_PASS pattern="
                          << requested_pattern << " rays=" << rays[ray][0U]
                          << ',' << rays[ray][1U] << ':'
                          << rays[ray + 1U][0U] << ','
                          << rays[ray + 1U][1U] << " bases="
                          << transformed_bases.size() << '\n';
            }
            std::cout << "SU2_OUTER_TEN_SLOPE_GLOBAL_PASS pattern="
                      << requested_pattern << " cones=" << rays.size() - 1U
                      << " bases=" << total_bases << " coefficients="
                      << total_coefficients << '\n';
            return EXIT_SUCCESS;
        }

        if (mode == "chamber-negatives") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "chamber-negatives mode requires a pattern"
                );
            }
            struct NegativeFamily {
                std::size_t occurrences = 0U;
                cpp_int minimum = 0;
                Counts transformed_base{};
                Counts base{};
            };
            constexpr std::array<std::array<int, 2>, 4> rays{
                std::array<int, 2>{1, 0},
                std::array<int, 2>{1, 1},
                std::array<int, 2>{1, 2},
                std::array<int, 2>{0, 1}
            };
            constexpr std::array<const char*, 3> names{
                "lower", "middle", "upper"
            };
            const Counts signs = signs_for_pattern(requested_pattern);
            for (std::size_t chamber = 0U; chamber < names.size();
                 ++chamber) {
                Directions directions{};
                for (std::size_t coordinate = 2U;
                     coordinate < directions.size(); ++coordinate) {
                    directions[coordinate][coordinate] = 1;
                }
                directions[0U][0U] = rays[chamber][0U];
                directions[0U][1U] = rays[chamber][1U];
                directions[1U][0U] = rays[chamber + 1U][0U];
                directions[1U][1U] = rays[chamber + 1U][1U];
                const std::vector<Counts> transformed_bases
                    = transformed_minimal_bases(directions);
                std::map<Counts, NegativeFamily> families;
                std::size_t negative_coefficients = 0U;
                long long tested = 0;
                for (const Counts& transformed_base : transformed_bases) {
                    const Counts base = apply_directions(
                        transformed_base, directions
                    );
                    Counts degree{};
                    for (int total_degree = 0; total_degree <= deficit;
                         ++total_degree) {
                        const auto visit = [&](const auto& self,
                                               int coordinate,
                                               int remaining) -> void {
                            if (coordinate == label_count) {
                                if (remaining != 0) {
                                    return;
                                }
                                const cpp_int coefficient
                                    = cone_finite_difference(
                                        base, directions, degree, signs
                                    );
                                ++tested;
                                if (coefficient >= 0) {
                                    return;
                                }
                                ++negative_coefficients;
                                NegativeFamily& family = families[degree];
                                ++family.occurrences;
                                if (family.occurrences == 1U
                                    || coefficient < family.minimum) {
                                    family.minimum = coefficient;
                                    family.transformed_base
                                        = transformed_base;
                                    family.base = base;
                                }
                                return;
                            }
                            for (int value = 0; value <= remaining;
                                 ++value) {
                                degree[static_cast<std::size_t>(coordinate)]
                                    = value;
                                self(
                                    self, coordinate + 1,
                                    remaining - value
                                );
                            }
                        };
                        visit(visit, 0, total_degree);
                    }
                }
                std::cout << "SU2_OUTER_TEN_CHAMBER_NEGATIVE_SUMMARY pattern="
                          << requested_pattern << " chamber=" << names[chamber]
                          << " bases=" << transformed_bases.size()
                          << " tested=" << tested
                          << " negatives=" << negative_coefficients
                          << " families=" << families.size() << '\n';
                for (const auto& [degree, family] : families) {
                    std::cout << "SU2_OUTER_TEN_CHAMBER_NEGATIVE_FAMILY pattern="
                              << requested_pattern << " chamber="
                              << names[chamber] << " degree=";
                    print_counts(degree);
                    std::cout << " occurrences=" << family.occurrences
                              << " minimum=" << family.minimum
                              << " transformed_base=";
                    print_counts(family.transformed_base);
                    std::cout << " base=";
                    print_counts(family.base);
                    std::cout << '\n';
                }
            }
            return EXIT_SUCCESS;
        }

        if (mode == "chamber-slices") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "chamber-slices mode requires a pattern"
                );
            }
            constexpr std::array<std::array<int, 2>, 4> rays{
                std::array<int, 2>{1, 0},
                std::array<int, 2>{1, 1},
                std::array<int, 2>{1, 2},
                std::array<int, 2>{0, 1}
            };
            constexpr std::array<const char*, 3> names{
                "lower", "middle", "upper"
            };
            const Counts signs = signs_for_pattern(requested_pattern);
            const Counts transformed_base{0, 0, 0, 0, 3};
            for (std::size_t chamber = 0U; chamber < names.size();
                 ++chamber) {
                Directions directions{};
                for (std::size_t coordinate = 2U;
                     coordinate < directions.size(); ++coordinate) {
                    directions[coordinate][coordinate] = 1;
                }
                directions[0U][0U] = rays[chamber][0U];
                directions[0U][1U] = rays[chamber][1U];
                directions[1U][0U] = rays[chamber + 1U][0U];
                directions[1U][1U] = rays[chamber + 1U][1U];
                const Counts base = apply_directions(
                    transformed_base, directions
                );
                std::cout << "SU2_OUTER_TEN_CHAMBER_SLICE pattern="
                          << requested_pattern << " chamber=" << names[chamber]
                          << " transformed_base=";
                print_counts(transformed_base);
                std::cout << " base=";
                print_counts(base);
                std::cout << '\n';
                for (int first = 0; first <= deficit; ++first) {
                    for (int second = 0; first + second <= deficit;
                         ++second) {
                        const Counts degree{first, second, 0, 0, 0};
                        const cpp_int coefficient = cone_finite_difference(
                            base, directions, degree, signs
                        );
                        if (coefficient == 0) {
                            continue;
                        }
                        std::cout << "degree=";
                        print_counts(degree);
                        std::cout << " value=" << coefficient << '\n';
                    }
                }
            }
            return EXIT_SUCCESS;
        }

        if (mode == "chamber-recursive") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "chamber-recursive mode requires a pattern"
                );
            }
            constexpr std::array<std::array<int, 2>, 4> rays{
                std::array<int, 2>{1, 0},
                std::array<int, 2>{1, 1},
                std::array<int, 2>{1, 2},
                std::array<int, 2>{0, 1}
            };
            constexpr std::array<const char*, 3> names{
                "lower", "middle", "upper"
            };
            const Counts signs = signs_for_pattern(requested_pattern);
            for (std::size_t chamber = 0U; chamber < names.size();
                 ++chamber) {
                Directions directions{};
                for (std::size_t coordinate = 2U;
                     coordinate < directions.size(); ++coordinate) {
                    directions[coordinate][coordinate] = 1;
                }
                directions[0U][0U] = rays[chamber][0U];
                directions[0U][1U] = rays[chamber][1U];
                directions[1U][0U] = rays[chamber + 1U][0U];
                directions[1U][1U] = rays[chamber + 1U][1U];
                const std::vector<Counts> transformed_bases
                    = transformed_minimal_bases(directions);
                RecursiveConeStats stats;
                bool passes = true;
                Counts failed_base{};
                for (const Counts& transformed_base : transformed_bases) {
                    std::map<RecursiveKey, bool> memo;
                    if (!recursive_cone_certificate(
                            transformed_base,
                            (1U << label_count) - 1U,
                            directions, signs, 0U, stats, memo
                        )) {
                        passes = false;
                        failed_base = transformed_base;
                        break;
                    }
                }
                std::cout << "SU2_OUTER_TEN_CHAMBER_RECURSIVE_"
                          << (passes ? "PASS" : "FAIL")
                          << " pattern=" << requested_pattern
                          << " chamber=" << names[chamber]
                          << " bases=" << transformed_bases.size()
                          << " nodes=" << stats.nodes
                          << " leaves=" << stats.leaves
                          << " splits=" << stats.splits
                          << " maximum_depth=" << stats.maximum_depth
                          << " coefficients=" << stats.coefficients
                          << " maximum_transformed_base=";
                print_counts(stats.maximum_transformed_base);
                std::cout << " maximum_base=";
                print_counts(stats.maximum_base);
                if (!passes) {
                    std::cout << " failed_transformed_base=";
                    print_counts(failed_base);
                    std::cout << " failure_degree=";
                    print_counts(stats.failure_degree);
                }
                std::cout << '\n';
                outer_cache().clear();
            }
            return EXIT_SUCCESS;
        }

        if (mode == "parabolic-tail") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "parabolic-tail mode requires a pattern"
                );
            }
            Directions directions{};
            directions[0U] = Counts{1, 2, 0, 0, 0};
            directions[1U] = Counts{0, 1, 0, 0, 0};
            for (std::size_t coordinate = 2U;
                 coordinate < directions.size(); ++coordinate) {
                directions[coordinate][coordinate] = 1;
            }
            const Counts base{0, 0, 0, 0, 3};
            const Counts signs = signs_for_pattern(requested_pattern);
            constexpr std::array<int, 6> expected{
                42, -14, 6, -2, 2, 10
            };
            constexpr std::array<long long, 6> scaled{
                1, -30, 360, -1200, 3600, 7200
            };
            bool passes = true;
            for (int second = 0; second <= 5; ++second) {
                const int first = 10 - 2 * second;
                const Counts degree{first, second, 0, 0, 0};
                const cpp_int coefficient = cone_finite_difference(
                    base, directions, degree, signs
                );
                if (coefficient != expected[static_cast<std::size_t>(second)]) {
                    passes = false;
                }
                std::cout << "SU2_OUTER_TEN_PARABOLIC_COEFFICIENT pattern="
                          << requested_pattern << " degree=";
                print_counts(degree);
                std::cout << " value=" << coefficient
                          << " expected="
                          << expected[static_cast<std::size_t>(second)]
                          << '\n';
            }
            constexpr std::array<long long, 6> decomposition{
                1, -30, 225 + 15 * 9, 15 * -80,
                15 * 240, 15 * 480
            };
            if (scaled != decomposition) {
                passes = false;
            }
            constexpr long long radical_margin
                = 241LL * 241LL - 3LL * 120LL * 120LL;
            if (radical_margin <= 0) {
                passes = false;
            }
            std::cout << "SU2_OUTER_TEN_PARABOLIC_TAIL_"
                      << (passes ? "PASS" : "FAIL")
                      << " pattern=" << requested_pattern
                      << " scaled_polynomial=1,-30,360,-1200,3600,7200"
                      << " radical_margin=" << radical_margin << '\n';
            return passes ? EXIT_SUCCESS : EXIT_FAILURE;
        }

        if (mode == "upper-monomial") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "upper-monomial mode requires a pattern"
                );
            }
            Directions directions{};
            directions[0U] = Counts{1, 2, 0, 0, 0};
            directions[1U] = Counts{0, 1, 0, 0, 0};
            for (std::size_t coordinate = 2U;
                 coordinate < directions.size(); ++coordinate) {
                directions[coordinate][coordinate] = 1;
            }
            const Counts base{0, 0, 0, 0, 3};
            const Counts signs = signs_for_pattern(requested_pattern);
            std::array<cpp_int, deficit + 1> factorial{};
            factorial[0U] = 1;
            for (int degree = 1; degree <= deficit; ++degree) {
                factorial[static_cast<std::size_t>(degree)]
                    = degree * factorial[static_cast<std::size_t>(degree - 1)];
            }
            std::array<std::array<cpp_int, deficit + 1>, deficit + 1>
                stirling{};
            stirling[0U][0U] = 1;
            for (int degree = 1; degree <= deficit; ++degree) {
                for (int power = 1; power <= degree; ++power) {
                    stirling[static_cast<std::size_t>(degree)]
                             [static_cast<std::size_t>(power)]
                        = stirling[static_cast<std::size_t>(degree - 1)]
                                  [static_cast<std::size_t>(power - 1)]
                        - (degree - 1)
                            * stirling[static_cast<std::size_t>(degree - 1)]
                                       [static_cast<std::size_t>(power)];
                }
            }
            std::array<std::array<cpp_int, deficit + 1>, deficit + 1>
                monomial{};
            for (int first_degree = 0; first_degree <= deficit;
                 ++first_degree) {
                for (int second_degree = 0;
                     first_degree + second_degree <= deficit;
                     ++second_degree) {
                    const Counts degree{
                        first_degree, second_degree, 0, 0, 0
                    };
                    const cpp_int newton = cone_finite_difference(
                        base, directions, degree, signs
                    );
                    const cpp_int scale
                        = factorial[static_cast<std::size_t>(deficit)]
                        / (factorial[static_cast<std::size_t>(first_degree)]
                           * factorial[static_cast<std::size_t>(second_degree)]);
                    for (int first_power = 0;
                         first_power <= first_degree; ++first_power) {
                        for (int second_power = 0;
                             second_power <= second_degree; ++second_power) {
                            monomial[static_cast<std::size_t>(first_power)]
                                    [static_cast<std::size_t>(second_power)]
                                += newton * scale
                                * stirling[static_cast<std::size_t>(first_degree)]
                                           [static_cast<std::size_t>(first_power)]
                                * stirling[static_cast<std::size_t>(second_degree)]
                                           [static_cast<std::size_t>(second_power)];
                        }
                    }
                }
            }
            std::cout << "SU2_OUTER_TEN_UPPER_MONOMIAL pattern="
                      << requested_pattern << " scale=10!\n";
            for (int first_power = 0; first_power <= deficit;
                 ++first_power) {
                for (int second_power = 0; second_power <= deficit;
                     ++second_power) {
                    const cpp_int& coefficient
                        = monomial[static_cast<std::size_t>(first_power)]
                                  [static_cast<std::size_t>(second_power)];
                    if (coefficient == 0) {
                        continue;
                    }
                    std::cout << "power=" << first_power << ','
                              << second_power << " coefficient="
                              << coefficient << '\n';
                }
            }
            return EXIT_SUCCESS;
        }

        if (mode == "parabolic-slice-certificate") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "parabolic-slice-certificate mode requires a pattern"
                );
            }
            const Counts signs = signs_for_pattern(requested_pattern);
            const Counts base{0, 0, 0, 0, 3};
            const BivariateMonomial monomial
                = upper_bivariate_monomial(base, signs);
            const std::map<std::array<int, 2>, cpp_int> expected_negatives{
                {{4, 3}, -50400}, {{5, 2}, -378000},
                {{6, 1}, -995400}, {{7, 0}, -913500},
                {{8, 1}, -1260}, {{9, 0}, -3150}
            };
            bool passes = true;
            std::size_t negative_count = 0U;
            for (int first = 0; first <= deficit; ++first) {
                for (int second = 0; second <= deficit; ++second) {
                    const cpp_int& coefficient
                        = monomial[static_cast<std::size_t>(first)]
                                  [static_cast<std::size_t>(second)];
                    if (coefficient >= 0) {
                        continue;
                    }
                    ++negative_count;
                    const auto found = expected_negatives.find(
                        std::array<int, 2>{first, second}
                    );
                    if (found == expected_negatives.end()
                        || coefficient != found->second) {
                        passes = false;
                    }
                }
            }
            if (negative_count != expected_negatives.size()
                || monomial[10U][0U] != 42
                || monomial[8U][1U] != -1260
                || monomial[6U][2U] != 15120
                || monomial[4U][3U] != -50400
                || monomial[2U][4U] != 151200
                || monomial[0U][5U] != 302400) {
                passes = false;
            }
            // K(z)>=3 follows from its exact minimum because
            // 214-120*sqrt(3)>0.  Hence H(z)>=1/6, H(z)>=2z,
            // and H(z)>=45z^2.  These bounds dominate all four
            // lower-weight negative monomials once u>=700.
            constexpr long long cubic_margin_above_three
                = 214LL * 214LL - 3LL * 120LL * 120LL;
            constexpr long long tail_denominator = 700LL * 700LL * 700LL;
            constexpr long long tail_numerator
                = 650LL * 700LL * 700LL
                + 11850LL * 700LL + 130500LL;
            if (cubic_margin_above_three <= 0
                || tail_numerator > tail_denominator) {
                passes = false;
            }
            constexpr int cutoff = 700;
            int maximum_positive_roots = 0;
            int failed_u = -1;
            for (int first = 0; first < cutoff && passes; ++first) {
                std::array<cpp_int, 6> coefficients{};
                cpp_int power = 1;
                for (int first_power = 0; first_power <= deficit;
                     ++first_power) {
                    for (int second_power = 0; second_power <= 5;
                         ++second_power) {
                        coefficients[static_cast<std::size_t>(second_power)]
                            += monomial[static_cast<std::size_t>(first_power)]
                                       [static_cast<std::size_t>(second_power)]
                            * power;
                    }
                    power *= first;
                }
                if (coefficients[0U] <= 0) {
                    passes = false;
                    failed_u = first;
                    break;
                }
                const int roots = positive_root_count(coefficients);
                maximum_positive_roots = std::max(
                    maximum_positive_roots, roots
                );
                if (roots != 0) {
                    passes = false;
                    failed_u = first;
                    break;
                }
            }
            std::cout << "SU2_OUTER_TEN_PARABOLIC_SLICE_"
                      << (passes ? "PASS" : "FAIL")
                      << " pattern=" << requested_pattern
                      << " integer_u_checked=" << cutoff
                      << " tail_u=" << cutoff
                      << " negative_monomials=" << negative_count
                      << " maximum_positive_roots="
                      << maximum_positive_roots
                      << " cubic_margin_above_three="
                      << cubic_margin_above_three
                      << " tail_margin="
                      << tail_denominator - tail_numerator;
            if (!passes) {
                std::cout << " failed_u=" << failed_u;
            }
            std::cout << '\n';
            return passes ? EXIT_SUCCESS : EXIT_FAILURE;
        }

        if (mode == "upper-slice-all") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "upper-slice-all mode requires a pattern"
                );
            }
            Directions directions{};
            directions[0U] = Counts{1, 2, 0, 0, 0};
            directions[1U] = Counts{0, 1, 0, 0, 0};
            for (std::size_t coordinate = 2U;
                 coordinate < directions.size(); ++coordinate) {
                directions[coordinate][coordinate] = 1;
            }
            const std::vector<Counts> transformed_bases
                = transformed_minimal_bases(directions);
            const Counts signs = signs_for_pattern(requested_pattern);
            int maximum_cutoff = 0;
            int maximum_positive_roots = 0;
            std::size_t maximum_negative_monomials = 0U;
            std::size_t passed = 0U;
            for (const Counts& transformed_base : transformed_bases) {
                const Counts base = apply_directions(
                    transformed_base, directions
                );
                ParabolicSliceStats stats;
                if (!certify_upper_bivariate_slice(base, signs, stats)) {
                    std::cout << "SU2_OUTER_TEN_UPPER_SLICE_FAIL pattern="
                              << requested_pattern << " transformed_base=";
                    print_counts(transformed_base);
                    std::cout << " base=";
                    print_counts(base);
                    std::cout << '\n';
                    return EXIT_SUCCESS;
                }
                ++passed;
                maximum_cutoff = std::max(maximum_cutoff, stats.cutoff);
                maximum_positive_roots = std::max(
                    maximum_positive_roots,
                    stats.maximum_positive_roots
                );
                maximum_negative_monomials = std::max(
                    maximum_negative_monomials,
                    stats.negative_monomials
                );
            }
            std::cout << "SU2_OUTER_TEN_UPPER_SLICE_PASS pattern="
                      << requested_pattern << " bases=" << passed
                      << " maximum_cutoff=" << maximum_cutoff
                      << " maximum_positive_roots="
                      << maximum_positive_roots
                      << " maximum_negative_monomials="
                      << maximum_negative_monomials << '\n';
            return EXIT_SUCCESS;
        }

        if (mode == "upper-hybrid") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "upper-hybrid mode requires a pattern"
                );
            }
            Directions directions{};
            directions[0U] = Counts{1, 2, 0, 0, 0};
            directions[1U] = Counts{0, 1, 0, 0, 0};
            for (std::size_t coordinate = 2U;
                 coordinate < directions.size(); ++coordinate) {
                directions[coordinate][coordinate] = 1;
            }
            const std::vector<Counts> transformed_bases
                = transformed_minimal_bases(directions);
            const Counts signs = signs_for_pattern(requested_pattern);
            RecursiveHybridStats stats;
            bool passes = true;
            Counts failed_base{};
            for (const Counts& transformed_base : transformed_bases) {
                std::map<RecursiveKey, bool> memo;
                if (!recursive_upper_hybrid_certificate(
                        transformed_base,
                        (1U << label_count) - 1U,
                        directions, signs, 0U, stats, memo
                    )) {
                    passes = false;
                    failed_base = transformed_base;
                    break;
                }
            }
            std::cout << "SU2_OUTER_TEN_UPPER_HYBRID_"
                      << (passes ? "PASS" : "FAIL")
                      << " pattern=" << requested_pattern
                      << " bases=" << transformed_bases.size()
                      << " nodes=" << stats.nodes
                      << " leaves=" << stats.leaves
                      << " splits=" << stats.splits
                      << " parabolic_checks=" << stats.parabolic_checks
                      << " maximum_depth=" << stats.maximum_depth
                      << " coefficients=" << stats.coefficients
                      << " maximum_parabolic_cutoff="
                      << stats.maximum_parabolic_cutoff
                      << " maximum_transformed_base=";
            print_counts(stats.maximum_transformed_base);
            std::cout << " maximum_base=";
            print_counts(stats.maximum_base);
            if (!passes) {
                std::cout << " failed_transformed_base=";
                print_counts(failed_base);
                std::cout << " failure_degree=";
                print_counts(stats.failure_degree);
                std::cout << " failure_reason=" << stats.failure_reason;
            }
            std::cout << '\n';
            return passes ? EXIT_SUCCESS : EXIT_FAILURE;
        }

        if (mode == "parabolic-leading-full") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "parabolic-leading-full mode requires a pattern"
                );
            }
            Directions directions{};
            directions[0U] = Counts{1, 2, 0, 0, 0};
            directions[1U] = Counts{0, 1, 0, 0, 0};
            for (std::size_t coordinate = 2U;
                 coordinate < directions.size(); ++coordinate) {
                directions[coordinate][coordinate] = 1;
            }
            const Counts base{0, 0, 0, 0, 3};
            const Counts signs = signs_for_pattern(requested_pattern);
            std::array<cpp_int, deficit + 1> factorial{};
            factorial[0U] = 1;
            for (int degree = 1; degree <= deficit; ++degree) {
                factorial[static_cast<std::size_t>(degree)]
                    = degree
                    * factorial[static_cast<std::size_t>(degree - 1)];
            }
            std::array<std::array<cpp_int, deficit + 1>, deficit + 1>
                stirling{};
            stirling[0U][0U] = 1;
            for (int degree = 1; degree <= deficit; ++degree) {
                for (int power = 1; power <= degree; ++power) {
                    stirling[static_cast<std::size_t>(degree)]
                             [static_cast<std::size_t>(power)]
                        = stirling[static_cast<std::size_t>(degree - 1)]
                                  [static_cast<std::size_t>(power - 1)]
                        - (degree - 1)
                            * stirling[static_cast<std::size_t>(degree - 1)]
                                       [static_cast<std::size_t>(power)];
                }
            }
            std::map<Counts, cpp_int> monomial;
            std::map<Counts, cpp_int> leading;
            Counts degree{};
            for (int total_degree = 0; total_degree <= deficit;
                 ++total_degree) {
                const auto visit_degree = [&](const auto& self,
                                              int coordinate,
                                              int remaining) -> void {
                    if (coordinate == label_count) {
                        if (remaining != 0) {
                            return;
                        }
                        const cpp_int newton = cone_finite_difference(
                            base, directions, degree, signs
                        );
                        cpp_int scale = factorial[deficit];
                        for (int value : degree) {
                            scale /= factorial[
                                static_cast<std::size_t>(value)
                            ];
                        }
                        Counts power{};
                        const auto visit_power = [&](const auto& power_self,
                                                     int axis) -> void {
                            if (axis == label_count) {
                                cpp_int term = newton * scale;
                                for (std::size_t index = 0U;
                                     index < power.size(); ++index) {
                                    term *= stirling[
                                        static_cast<std::size_t>(degree[index])
                                    ][static_cast<std::size_t>(power[index])];
                                }
                                monomial[power] += term;
                                int weight = power[0U];
                                for (std::size_t index = 1U;
                                     index < power.size(); ++index) {
                                    weight += 2 * power[index];
                                }
                                if (weight == deficit) {
                                    leading[power] += term;
                                }
                                return;
                            }
                            for (int value = 0;
                                 value <= degree[
                                     static_cast<std::size_t>(axis)
                                 ]; ++value) {
                                power[static_cast<std::size_t>(axis)] = value;
                                power_self(power_self, axis + 1);
                            }
                        };
                        visit_power(visit_power, 0);
                        return;
                    }
                    for (int value = 0; value <= remaining; ++value) {
                        degree[static_cast<std::size_t>(coordinate)] = value;
                        self(self, coordinate + 1, remaining - value);
                    }
                };
                visit_degree(visit_degree, 0, total_degree);
            }
            std::size_t nonzero = 0U;
            std::size_t negative = 0U;
            for (const auto& [power, coefficient] : leading) {
                if (coefficient == 0) {
                    continue;
                }
                ++nonzero;
                if (coefficient < 0) {
                    ++negative;
                }
                std::cout << "SU2_OUTER_TEN_PARABOLIC_LEADING pattern="
                          << requested_pattern << " power=";
                print_counts(power);
                std::cout << " coefficient=" << coefficient << '\n';
            }
            const auto scalar_coefficient = [](int first, int second,
                                               int high_degree) -> long long {
                if (first == 10 && second == 0 && high_degree == 0) {
                    return 42;
                }
                if (first == 8) {
                    if (second == 0 && high_degree == 1) {
                        return 1260;
                    }
                    if (second == 1 && high_degree == 0) {
                        return -1260;
                    }
                }
                if (first == 6) {
                    if (second == 0 && high_degree == 2) {
                        return 12600;
                    }
                    if (second == 1 && high_degree == 1) {
                        return -20160;
                    }
                    if (second == 2 && high_degree == 0) {
                        return 15120;
                    }
                }
                if (first == 4) {
                    if (second == 0 && high_degree == 3) {
                        return 50400;
                    }
                    if (second == 1 && high_degree == 2) {
                        return -75600;
                    }
                    if (second == 2 && high_degree == 1) {
                        return 151200;
                    }
                    if (second == 3 && high_degree == 0) {
                        return -50400;
                    }
                }
                if (first == 2) {
                    if (second == 0 && high_degree == 4) {
                        return 75600;
                    }
                    if (second == 2 && high_degree == 2) {
                        return 453600;
                    }
                    if (second == 4 && high_degree == 0) {
                        return 151200;
                    }
                }
                if (first == 0) {
                    constexpr std::array<long long, 6> coefficients{
                        30240, 151200, 604800, 907200, 907200, 302400
                    };
                    if (second + high_degree == 5) {
                        return coefficients[
                            static_cast<std::size_t>(second)
                        ];
                    }
                }
                return 0;
            };
            bool identity_passes = true;
            Counts power{};
            const auto audit_power = [&](const auto& self, int coordinate,
                                         int remaining_weight) -> void {
                if (coordinate == label_count) {
                    if (remaining_weight != 0) {
                        return;
                    }
                    const int high_degree = power[2U] + power[3U] + power[4U];
                    const long long scalar = scalar_coefficient(
                        power[0U], power[1U], high_degree
                    );
                    cpp_int expected = scalar;
                    expected *= factorial[
                        static_cast<std::size_t>(high_degree)
                    ];
                    for (std::size_t index = 2U; index < power.size();
                         ++index) {
                        expected /= factorial[
                            static_cast<std::size_t>(power[index])
                        ];
                    }
                    if (leading[power] != expected) {
                        identity_passes = false;
                    }
                    return;
                }
                const int weight = coordinate == 0 ? 1 : 2;
                for (int value = 0; value * weight <= remaining_weight;
                     ++value) {
                    power[static_cast<std::size_t>(coordinate)] = value;
                    self(
                        self, coordinate + 1,
                        remaining_weight - value * weight
                    );
                }
            };
            audit_power(audit_power, 0, deficit);
            const auto lower_bound = [](int second,
                                        int high_degree) -> long long {
                if (high_degree == 0) {
                    constexpr std::array<long long, 5> bounds{
                        7, 84, 1890, 42, 42
                    };
                    return second <= 4
                        ? bounds[static_cast<std::size_t>(second)] : 0;
                }
                if (high_degree == 1) {
                    constexpr std::array<long long, 4> bounds{
                        588, 1260, 1260, 1260
                    };
                    return second <= 3
                        ? bounds[static_cast<std::size_t>(second)] : 0;
                }
                if (high_degree == 2) {
                    constexpr std::array<long long, 3> bounds{
                        9450, 12600, 12600
                    };
                    return second <= 2
                        ? bounds[static_cast<std::size_t>(second)] : 0;
                }
                if (high_degree == 3) {
                    constexpr std::array<long long, 2> bounds{50400, 50400};
                    return second <= 1
                        ? bounds[static_cast<std::size_t>(second)] : 0;
                }
                if (high_degree == 4 && second == 0) {
                    return 75600;
                }
                return 0;
            };
            constexpr std::array<long long, 6> a_coefficients{
                1, -6, 36, 72, 0, 0
            };
            constexpr std::array<long long, 6> b_coefficients{
                1, -16, 120, 0, 720, 0
            };
            for (int power_index = 1; power_index <= 2; ++power_index) {
                std::array<cpp_int, 6> difference{};
                for (std::size_t degree_index = 0U;
                     degree_index < difference.size(); ++degree_index) {
                    difference[degree_index] = a_coefficients[degree_index];
                }
                --difference[static_cast<std::size_t>(power_index)];
                if (positive_root_count(difference) != 0) {
                    identity_passes = false;
                }
            }
            for (int power_index = 1; power_index <= 3; ++power_index) {
                std::array<cpp_int, 6> difference{};
                for (std::size_t degree_index = 0U;
                     degree_index < difference.size(); ++degree_index) {
                    difference[degree_index] = b_coefficients[degree_index];
                }
                --difference[static_cast<std::size_t>(power_index)];
                if (positive_root_count(difference) != 0) {
                    identity_passes = false;
                }
            }
            struct LowerNegativeTerm {
                cpp_int magnitude;
                long long bound = 0;
                int gap = 0;
            };
            std::vector<LowerNegativeTerm> lower_negatives;
            for (const auto& [monomial_power, coefficient] : monomial) {
                if (coefficient == 0) {
                    continue;
                }
                int weight = monomial_power[0U];
                int high_degree = 0;
                for (std::size_t index = 1U;
                     index < monomial_power.size(); ++index) {
                    weight += 2 * monomial_power[index];
                    if (index >= 2U) {
                        high_degree += monomial_power[index];
                    }
                }
                if (weight > deficit) {
                    identity_passes = false;
                }
                if (coefficient < 0 && weight < deficit) {
                    const long long bound = lower_bound(
                        monomial_power[1U], high_degree
                    );
                    if (bound <= 0) {
                        identity_passes = false;
                    } else {
                        lower_negatives.push_back(
                            LowerNegativeTerm{
                                -coefficient, bound, deficit - weight
                            }
                        );
                        std::cout
                            << "SU2_OUTER_TEN_PARABOLIC_LOWER_NEGATIVE pattern="
                            << requested_pattern << " power=";
                        print_counts(monomial_power);
                        std::cout << " coefficient=" << coefficient
                                  << " weight=" << weight << '\n';
                    }
                }
            }
            int full_tail_cutoff = 1;
            constexpr int maximum_full_tail_cutoff = 1048576;
            while (identity_passes) {
                Rational ratio_sum{0};
                for (const LowerNegativeTerm& term : lower_negatives) {
                    ratio_sum += Rational{
                        term.magnitude,
                        term.bound
                            * integer_power(full_tail_cutoff, term.gap)
                    };
                }
                if (ratio_sum <= 1) {
                    break;
                }
                if (full_tail_cutoff >= maximum_full_tail_cutoff) {
                    identity_passes = false;
                    break;
                }
                full_tail_cutoff *= 2;
            }
            std::cout << "SU2_OUTER_TEN_PARABOLIC_LEADING_SUMMARY pattern="
                      << requested_pattern << " nonzero=" << nonzero
                      << " negative=" << negative << " scale=10!"
                      << " two_variable_identity="
                      << (identity_passes ? "PASS" : "FAIL")
                      << " positivity_decomposition="
                      << (identity_passes ? "PASS" : "FAIL")
                      << " lower_negative_monomials="
                      << lower_negatives.size()
                      << " full_tail_u=" << full_tail_cutoff << '\n';
            return identity_passes ? EXIT_SUCCESS : EXIT_FAILURE;
        }

        if (mode == "upper-fixed-slices") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "upper-fixed-slices mode requires a pattern"
                );
            }
            Directions directions{};
            directions[0U] = Counts{1, 2, 0, 0, 0};
            directions[1U] = Counts{0, 1, 0, 0, 0};
            for (std::size_t coordinate = 2U;
                 coordinate < directions.size(); ++coordinate) {
                directions[coordinate][coordinate] = 1;
            }
            const Counts signs = signs_for_pattern(requested_pattern);
            RecursiveConeStats totals;
            int passed = 0;
            for (int first = 0; first < maximum_total_count; ++first) {
                const Counts transformed_base{first, 0, 0, 0, 3};
                RecursiveConeStats stats;
                std::map<RecursiveKey, bool> memo;
                constexpr unsigned int free_mask
                    = ((1U << label_count) - 1U) & ~1U;
                if (!recursive_cone_certificate(
                        transformed_base, free_mask, directions, signs, 0U,
                        stats, memo
                    )) {
                    std::cout << "SU2_OUTER_TEN_UPPER_FIXED_SLICES_FAIL pattern="
                              << requested_pattern << " first=" << first
                              << " nodes=" << stats.nodes
                              << " maximum_depth=" << stats.maximum_depth
                              << " failure_degree=";
                    print_counts(stats.failure_degree);
                    std::cout << '\n';
                    return EXIT_SUCCESS;
                }
                ++passed;
                totals.nodes += stats.nodes;
                totals.leaves += stats.leaves;
                totals.splits += stats.splits;
                totals.coefficients += stats.coefficients;
                totals.maximum_depth = std::max(
                    totals.maximum_depth, stats.maximum_depth
                );
                for (std::size_t coordinate = 0U;
                     coordinate < totals.maximum_base.size(); ++coordinate) {
                    totals.maximum_transformed_base[coordinate] = std::max(
                        totals.maximum_transformed_base[coordinate],
                        stats.maximum_transformed_base[coordinate]
                    );
                    totals.maximum_base[coordinate] = std::max(
                        totals.maximum_base[coordinate],
                        stats.maximum_base[coordinate]
                    );
                }
            }
            std::cout << "SU2_OUTER_TEN_UPPER_FIXED_SLICES_PASS pattern="
                      << requested_pattern << " slices=" << passed
                      << " nodes=" << totals.nodes
                      << " leaves=" << totals.leaves
                      << " splits=" << totals.splits
                      << " maximum_depth=" << totals.maximum_depth
                      << " coefficients=" << totals.coefficients
                      << " maximum_transformed_base=";
            print_counts(totals.maximum_transformed_base);
            std::cout << " maximum_base=";
            print_counts(totals.maximum_base);
            std::cout << '\n';
            return EXIT_SUCCESS;
        }

        if (mode == "cylindrical-sturm") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "cylindrical-sturm mode requires a pattern"
                );
            }
            constexpr std::array<std::array<int, 2>, 4> rays{
                std::array<int, 2>{1, 0},
                std::array<int, 2>{1, 1},
                std::array<int, 2>{1, 2},
                std::array<int, 2>{0, 1}
            };
            constexpr std::array<const char*, 3> names{
                "lower", "middle", "upper"
            };
            const Counts signs = signs_for_pattern(requested_pattern);
            bool global_passes = true;
            for (std::size_t chamber = 0U; chamber < names.size();
                 ++chamber) {
                Directions directions{};
                for (std::size_t coordinate = 2U;
                     coordinate < directions.size(); ++coordinate) {
                    directions[coordinate][coordinate] = 1;
                }
                directions[0U][0U] = rays[chamber][0U];
                directions[0U][1U] = rays[chamber][1U];
                directions[1U][0U] = rays[chamber + 1U][0U];
                directions[1U][1U] = rays[chamber + 1U][1U];
                const std::vector<Counts> transformed_bases
                    = transformed_minimal_bases(directions);
                std::size_t direct_passes = 0U;
                std::size_t cylindrical_passes = 0U;
                CylindricalStats stats;
                bool chamber_passes = true;
                for (const Counts& transformed_base : transformed_bases) {
                    const Counts base = apply_directions(
                        transformed_base, directions
                    );
                    Counts failure{};
                    cpp_int value = 0;
                    long long tested = 0;
                    if (!first_negative_cone_coefficient(
                            base, directions, signs, failure, value, tested
                        )) {
                        ++direct_passes;
                        continue;
                    }
                    const std::array<std::size_t, 2> parameter_axes
                        = chamber == 1U
                        ? std::array<std::size_t, 2>{1U, 0U}
                        : std::array<std::size_t, 2>{0U, 1U};
                    bool base_passes = false;
                    for (std::size_t parameter_axis : parameter_axes) {
                        CylindricalStats attempt;
                        if (certify_cylindrical_orthant(
                                base, directions, signs, parameter_axis,
                                attempt
                            )) {
                            stats.polynomials += attempt.polynomials;
                            stats.nonzero_polynomials
                                += attempt.nonzero_polynomials;
                            stats.maximum_degree = std::max(
                                stats.maximum_degree,
                                attempt.maximum_degree
                            );
                            stats.maximum_positive_roots = std::max(
                                stats.maximum_positive_roots,
                                attempt.maximum_positive_roots
                            );
                            ++cylindrical_passes;
                            base_passes = true;
                            break;
                        }
                        std::cout
                            << "SU2_OUTER_TEN_CYLINDRICAL_AXIS_FAIL pattern="
                            << requested_pattern << " chamber="
                            << names[chamber] << " parameter_axis="
                            << parameter_axis << " transformed_base=";
                        print_counts(transformed_base);
                        std::cout << " failure_degree=";
                        print_counts(attempt.failure_degree);
                        std::cout << " maximum_positive_roots="
                                  << attempt.maximum_positive_roots << '\n';
                    }
                    if (!base_passes) {
                        chamber_passes = false;
                        global_passes = false;
                        std::cout
                            << "SU2_OUTER_TEN_CYLINDRICAL_STURM_FAIL pattern="
                            << requested_pattern << " chamber="
                            << names[chamber] << " transformed_base=";
                        print_counts(transformed_base);
                        std::cout << " base=";
                        print_counts(base);
                        std::cout << '\n';
                        break;
                    }
                }
                std::cout << "SU2_OUTER_TEN_CYLINDRICAL_STURM_"
                          << (chamber_passes ? "PASS" : "FAIL")
                          << " pattern=" << requested_pattern
                          << " chamber=" << names[chamber]
                          << " bases=" << transformed_bases.size()
                          << " direct_passes=" << direct_passes
                          << " cylindrical_passes=" << cylindrical_passes
                          << " polynomials=" << stats.polynomials
                          << " nonzero_polynomials="
                          << stats.nonzero_polynomials
                          << " maximum_degree=" << stats.maximum_degree
                          << " maximum_positive_roots="
                          << stats.maximum_positive_roots << '\n';
                if (!chamber_passes) {
                    outer_cache().clear();
                    continue;
                }
                outer_cache().clear();
            }
            std::cout << "SU2_OUTER_TEN_CYLINDRICAL_STURM_GLOBAL_"
                      << (global_passes ? "PASS" : "FAIL")
                      << " pattern=" << requested_pattern << '\n';
            return EXIT_SUCCESS;
        }

        if (mode == "bivariate-sturm") {
            if (requested_pattern < 0) {
                throw std::runtime_error(
                    "bivariate-sturm mode requires a pattern"
                );
            }
            constexpr std::array<std::array<int, 2>, 4> rays{
                std::array<int, 2>{1, 0},
                std::array<int, 2>{1, 1},
                std::array<int, 2>{1, 2},
                std::array<int, 2>{0, 1}
            };
            constexpr std::array<const char*, 3> names{
                "lower", "middle", "upper"
            };
            const Counts signs = signs_for_pattern(requested_pattern);
            bool global_passes = true;
            for (std::size_t chamber = 0U; chamber < names.size();
                 ++chamber) {
                Directions directions{};
                for (std::size_t coordinate = 2U;
                     coordinate < directions.size(); ++coordinate) {
                    directions[coordinate][coordinate] = 1;
                }
                directions[0U][0U] = rays[chamber][0U];
                directions[0U][1U] = rays[chamber][1U];
                directions[1U][0U] = rays[chamber + 1U][0U];
                directions[1U][1U] = rays[chamber + 1U][1U];
                const std::vector<Counts> transformed_bases
                    = transformed_minimal_bases(directions);
                std::size_t direct_passes = 0U;
                std::size_t bivariate_passes = 0U;
                WeightedBivariateStats stats;
                bool chamber_passes = true;
                for (const Counts& transformed_base : transformed_bases) {
                    const Counts base = apply_directions(
                        transformed_base, directions
                    );
                    Counts failure{};
                    cpp_int value = 0;
                    long long tested = 0;
                    if (!first_negative_cone_coefficient(
                            base, directions, signs, failure, value, tested
                        )) {
                        ++direct_passes;
                        continue;
                    }
                    WeightedBivariateStats attempt;
                    const int second_weight = chamber == 2U ? 2 : 1;
                    if (!certify_bivariate_newton_orthant(
                            base, directions, signs, second_weight, attempt
                        )) {
                        chamber_passes = false;
                        global_passes = false;
                        std::cout
                            << "SU2_OUTER_TEN_BIVARIATE_STURM_FAIL pattern="
                            << requested_pattern << " chamber="
                            << names[chamber] << " transformed_base=";
                        print_counts(transformed_base);
                        std::cout << " failure_degree=";
                        print_counts(attempt.failure_degree);
                        std::cout << " polynomials=" << attempt.polynomials
                                  << " maximum_cutoff="
                                  << attempt.maximum_cutoff
                                  << " failure_reason="
                                  << attempt.failure_reason << '\n';
                        break;
                    }
                    ++bivariate_passes;
                    stats.polynomials += attempt.polynomials;
                    stats.direct_polynomials += attempt.direct_polynomials;
                    stats.sturm_polynomials += attempt.sturm_polynomials;
                    stats.sturm_slices += attempt.sturm_slices;
                    stats.finite_newton_slices
                        += attempt.finite_newton_slices;
                    stats.maximum_cutoff = std::max(
                        stats.maximum_cutoff, attempt.maximum_cutoff
                    );
                    stats.maximum_weight = std::max(
                        stats.maximum_weight, attempt.maximum_weight
                    );
                }
                std::cout << "SU2_OUTER_TEN_BIVARIATE_STURM_"
                          << (chamber_passes ? "PASS" : "FAIL")
                          << " pattern=" << requested_pattern
                          << " chamber=" << names[chamber]
                          << " bases=" << transformed_bases.size()
                          << " direct_passes=" << direct_passes
                          << " bivariate_passes=" << bivariate_passes
                          << " polynomials=" << stats.polynomials
                          << " direct_polynomials="
                          << stats.direct_polynomials
                          << " sturm_polynomials=" << stats.sturm_polynomials
                          << " sturm_slices=" << stats.sturm_slices
                          << " finite_newton_slices="
                          << stats.finite_newton_slices
                          << " maximum_cutoff=" << stats.maximum_cutoff
                          << " maximum_weight=" << stats.maximum_weight
                          << '\n';
                outer_cache().clear();
            }
            std::cout << "SU2_OUTER_TEN_BIVARIATE_STURM_GLOBAL_"
                      << (global_passes ? "PASS" : "FAIL")
                      << " pattern=" << requested_pattern << '\n';
            return EXIT_SUCCESS;
        }

        if (mode == "certificate" || mode == "recursive") {
            const int first_pattern = requested_pattern < 0
                ? 0 : requested_pattern;
            const int last_pattern = requested_pattern < 0
                ? 15 : requested_pattern;
            for (int pattern = first_pattern; pattern <= last_pattern;
                 ++pattern) {
                const Counts signs = signs_for_pattern(pattern);
                if (mode == "recursive") {
                    RecursiveStats stats;
                    bool passes = true;
                    for (const Counts& base : bases) {
                        std::map<RecursiveKey, bool> memo;
                        if (!recursive_newton_certificate(
                                base, (1U << label_count) - 1U, signs, 0U,
                                stats, memo
                            )) {
                            std::cout
                                << "SU2_OUTER_TEN_RECURSIVE_FAIL pattern="
                                << pattern << " base=";
                            print_counts(base);
                            std::cout << " nodes=" << stats.nodes
                                      << " maximum_depth="
                                      << stats.maximum_depth
                                      << " failure_degree=";
                            print_counts(stats.failure_degree);
                            std::cout << '\n';
                            passes = false;
                            break;
                        }
                    }
                    if (passes) {
                        std::cout
                            << "SU2_OUTER_TEN_RECURSIVE_PASS pattern="
                            << pattern << " bases=" << bases.size()
                            << " nodes=" << stats.nodes
                            << " leaves=" << stats.leaves
                            << " splits=" << stats.splits
                            << " maximum_depth=" << stats.maximum_depth
                            << " coefficients=" << stats.coefficients
                            << " maximum_base=";
                        print_counts(stats.maximum_base);
                        std::cout << '\n';
                    }
                    outer_cache().clear();
                    continue;
                }
                bool passes = true;
                long long tested = 0;
                for (const Counts& base : bases) {
                    Counts negative_degree{};
                    cpp_int negative_value = 0;
                    if (first_negative_newton_coefficient(
                            base, signs, negative_degree, negative_value,
                            tested
                        )) {
                        std::cout << "SU2_OUTER_TEN_NEWTON_FAIL pattern="
                                  << pattern << " base=";
                        print_counts(base);
                        std::cout << " degree=";
                        print_counts(negative_degree);
                        std::cout << " value=" << negative_value
                                  << " tested=" << tested << '\n';
                        passes = false;
                        break;
                    }
                }
                if (passes) {
                    std::cout << "SU2_OUTER_TEN_NEWTON_PASS pattern="
                              << pattern << " bases=" << bases.size()
                              << " coefficients=" << tested << '\n';
                }
                outer_cache().clear();
            }
            return EXIT_SUCCESS;
        }

        cpp_int minimum = 0;
        Counts minimum_counts{};
        int minimum_pattern = -1;
        bool initialized = false;
        long long tested = 0;
        Counts counts{};
        const auto visit = [&](const auto& self, int coordinate) -> bool {
            if (coordinate == label_count) {
                if (total_count(counts) > maximum_total_count
                    || !admissible(counts)) {
                    return true;
                }
                for (int pattern = 0; pattern < 16; ++pattern) {
                    const cpp_int value = outer_value(
                        counts, signs_for_pattern(pattern)
                    );
                    ++tested;
                    if (!initialized || value < minimum) {
                        initialized = true;
                        minimum = value;
                        minimum_counts = counts;
                        minimum_pattern = pattern;
                    }
                    if (value < 0) {
                        std::cout << "SU2_OUTER_TEN_FAIL pattern=" << pattern
                                  << " counts=";
                        print_counts(counts);
                        std::cout << " value=" << value << '\n';
                        return false;
                    }
                }
                return true;
            }
            const int used = total_count(counts);
            for (int count = 0; count + used <= maximum_total_count; ++count) {
                counts[static_cast<std::size_t>(coordinate)] = count;
                if (!self(self, coordinate + 1)) {
                    return false;
                }
            }
            counts[static_cast<std::size_t>(coordinate)] = 0;
            return true;
        };
        if (!visit(visit, 0)) {
            return EXIT_FAILURE;
        }
        std::cout << "SU2_OUTER_TEN_SEARCH PASS tested=" << tested
                  << " maximum_total_count=" << maximum_total_count
                  << " minimum=" << minimum << " pattern="
                  << minimum_pattern << " counts=";
        print_counts(minimum_counts);
        std::cout << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
