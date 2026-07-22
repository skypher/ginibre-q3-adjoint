#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

using Counts = std::array<int, 4>;
using MultiplicityKey = std::array<int, 5>;
using OuterKey = std::array<int, 8>;
using RecursiveKey = std::array<int, 5>;
using Directions = std::array<Counts, 4>;
using AxisPolynomial = std::vector<std::pair<Counts, cpp_int>>;

std::map<MultiplicityKey, cpp_int>& multiplicity_cache() {
    static std::map<MultiplicityKey, cpp_int> cache;
    return cache;
}

std::map<OuterKey, cpp_int>& outer_cache() {
    static std::map<OuterKey, cpp_int> cache;
    return cache;
}

std::map<OuterKey, AxisPolynomial>& axis_polynomial_cache() {
    static std::map<OuterKey, AxisPolynomial> cache;
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

cpp_int binomial_big(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    cpp_int result = 1;
    for (int index = 1; index <= k; ++index) {
        result *= n - k + index;
        result /= index;
    }
    return result;
}

int total_label(const Counts& counts) {
    int result = 0;
    for (int index = 0; index < 4; ++index) {
        result += (index + 1) * counts[static_cast<std::size_t>(index)];
    }
    return result;
}

int maximum_label(const Counts& counts) {
    for (int index = 3; index >= 0; --index) {
        if (counts[static_cast<std::size_t>(index)] != 0) {
            return index + 1;
        }
    }
    return 0;
}

bool admissible(const Counts& counts) {
    const int maximum = maximum_label(counts);
    return maximum != 0 && total_label(counts) - 8 >= maximum;
}

cpp_int multiplicity(const Counts& counts, int target) {
    auto& cache = multiplicity_cache();
    const MultiplicityKey key{
        counts[0], counts[1], counts[2], counts[3], target
    };
    const auto found = cache.find(key);
    if (found != cache.end()) {
        return found->second;
    }
    const int total = total_label(counts);
    if (target < 0 || target > total) {
        return 0;
    }
    std::vector<cpp_int> current(static_cast<std::size_t>(total + 1));
    current[0] = 1;
    int support = 0;
    for (int label = 1; label <= 4; ++label) {
        for (int copy = 0;
             copy < counts[static_cast<std::size_t>(label - 1)];
             ++copy) {
            std::vector<cpp_int> next(current.size());
            for (int input = 0; input <= support; ++input) {
                if (current[static_cast<std::size_t>(input)] == 0) {
                    continue;
                }
                for (int output = std::abs(input - label);
                     output <= input + label;
                     output += 2) {
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
    auto& cache = outer_cache();
    const OuterKey key{
        counts[0], counts[1], counts[2], counts[3],
        signs[0], signs[1], signs[2], signs[3]
    };
    const auto found = cache.find(key);
    if (found != cache.end()) {
        return found->second;
    }
    const int target = total_label(counts) - 8;
    cpp_int result = 0;
    Counts selected{};
    const auto visit = [&](const auto& self, int coordinate,
                           int selected_label) -> void {
        if (coordinate == 4) {
            const cpp_int invariant = multiplicity(selected, 0);
            if (invariant == 0) {
                return;
            }
            Counts complement{};
            cpp_int weight = invariant;
            int parity = 0;
            for (int index = 0; index < 4; ++index) {
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
            (8 - selected_label) / label
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
    const int total_degree
        = degree[0] + degree[1] + degree[2] + degree[3];
    const auto visit = [&](const auto& self, int coordinate,
                           int used_degree, long long coefficient) -> void {
        if (coordinate == 4) {
            Counts point{};
            for (int index = 0; index < 4; ++index) {
                point[static_cast<std::size_t>(index)]
                    = base[static_cast<std::size_t>(index)]
                    + offset[static_cast<std::size_t>(index)];
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
                self,
                coordinate + 1,
                used_degree + value,
                coefficient * binomial(maximum, value)
            );
        }
    };
    visit(visit, 0, 0, 1);
    return result;
}

const AxisPolynomial& axis_polynomial(
    const Counts& base,
    const Counts& signs
) {
    const OuterKey key{
        base[0], base[1], base[2], base[3],
        signs[0], signs[1], signs[2], signs[3]
    };
    auto& cache = axis_polynomial_cache();
    const auto found = cache.find(key);
    if (found != cache.end()) {
        return found->second;
    }
    constexpr Counts degrees{8, 4, 4, 4};
    AxisPolynomial polynomial;
    for (int first = 0; first <= degrees[0]; ++first) {
        for (int second = 0; second <= degrees[1]; ++second) {
            for (int third = 0; third <= degrees[2]; ++third) {
                for (int fourth = 0; fourth <= degrees[3]; ++fourth) {
                    if (first + second + third + fourth > 8) {
                        continue;
                    }
                    const Counts degree{first, second, third, fourth};
                    cpp_int coefficient = finite_difference(
                        base, degree, signs
                    );
                    if (coefficient != 0) {
                        polynomial.emplace_back(
                            degree, std::move(coefficient)
                        );
                    }
                }
            }
        }
    }
    return cache.emplace(key, std::move(polynomial)).first->second;
}

cpp_int evaluate_axis_polynomial(
    const Counts& base,
    const Counts& point,
    const AxisPolynomial& polynomial
) {
    Counts offset{};
    for (std::size_t coordinate = 0U; coordinate < offset.size();
         ++coordinate) {
        offset[coordinate] = point[coordinate] - base[coordinate];
        if (offset[coordinate] < 0) {
            throw std::runtime_error("axis polynomial evaluated below base");
        }
    }
    cpp_int value = 0;
    for (const auto& [degree, coefficient] : polynomial) {
        cpp_int term = coefficient;
        for (std::size_t coordinate = 0U; coordinate < degree.size();
             ++coordinate) {
            term *= binomial_big(offset[coordinate], degree[coordinate]);
        }
        value += term;
    }
    return value;
}

void audit_axis_polynomial(
    const Counts& base,
    const Counts& signs,
    const AxisPolynomial& polynomial
) {
    const OuterKey key{
        base[0], base[1], base[2], base[3],
        signs[0], signs[1], signs[2], signs[3]
    };
    static std::map<OuterKey, bool> audited;
    if (audited.find(key) != audited.end()) {
        return;
    }
    constexpr std::array<Counts, 3> offsets{
        Counts{9, 0, 0, 0},
        Counts{2, 5, 1, 0},
        Counts{3, 2, 2, 5}
    };
    for (const Counts& offset : offsets) {
        Counts point{};
        for (std::size_t coordinate = 0U; coordinate < point.size();
             ++coordinate) {
            point[coordinate] = base[coordinate] + offset[coordinate];
        }
        if (evaluate_axis_polynomial(base, point, polynomial)
            != outer_value(point, signs)) {
            throw std::runtime_error("axis polynomial audit failed");
        }
    }
    audited.emplace(key, true);
}

cpp_int directional_difference_from_polynomial(
    const Counts& base,
    const Counts& direction,
    int shift,
    int degree,
    const Counts& signs
) {
    const AxisPolynomial& polynomial = axis_polynomial(base, signs);
    cpp_int result = 0;
    for (int sample = 0; sample <= degree; ++sample) {
        Counts point = base;
        for (std::size_t coordinate = 0U; coordinate < point.size();
             ++coordinate) {
            point[coordinate] += (shift + sample) * direction[coordinate];
        }
        const cpp_int term = binomial(degree, sample)
            * evaluate_axis_polynomial(base, point, polynomial);
        result += ((degree - sample) & 1) != 0 ? -term : term;
    }
    return result;
}

cpp_int cone_difference_from_polynomial(
    const Counts& base,
    const Directions& directions,
    const Counts& degree,
    const Counts& signs
) {
    const AxisPolynomial& polynomial = axis_polynomial(base, signs);
    cpp_int result = 0;
    Counts sample{};
    const int total_degree
        = degree[0] + degree[1] + degree[2] + degree[3];
    const auto visit = [&](const auto& self, std::size_t axis,
                           int used_degree,
                           const cpp_int& coefficient) -> void {
        if (axis == sample.size()) {
            Counts point = base;
            for (std::size_t direction = 0U;
                 direction < directions.size(); ++direction) {
                for (std::size_t coordinate = 0U;
                     coordinate < point.size(); ++coordinate) {
                    point[coordinate] += sample[direction]
                        * directions[direction][coordinate];
                }
            }
            const cpp_int term = coefficient
                * evaluate_axis_polynomial(base, point, polynomial);
            result += ((total_degree - used_degree) & 1) != 0
                ? -term : term;
            return;
        }
        for (int value = 0; value <= degree[axis]; ++value) {
            sample[axis] = value;
            self(
                self, axis + 1U, used_degree + value,
                coefficient * binomial(degree[axis], value)
            );
        }
    };
    visit(visit, 0U, 0, cpp_int{1});
    return result;
}

std::vector<Counts> minimal_bases() {
    std::vector<Counts> result;
    for (int first = 0; first <= 12; ++first) {
        for (int second = 0; second <= 6; ++second) {
            for (int third = 0; third <= 4; ++third) {
                for (int fourth = 0; fourth <= 3; ++fourth) {
                    const Counts counts{first, second, third, fourth};
                    if (!admissible(counts)) {
                        continue;
                    }
                    bool minimal = true;
                    for (int coordinate = 0; coordinate < 4; ++coordinate) {
                        Counts predecessor = counts;
                        const std::size_t position
                            = static_cast<std::size_t>(coordinate);
                        if (predecessor[position] != 0) {
                            --predecessor[position];
                            if (admissible(predecessor)) {
                                minimal = false;
                            }
                        }
                    }
                    if (minimal) {
                        result.push_back(counts);
                    }
                }
            }
        }
    }
    return result;
}

Counts apply_directions(
    const Counts& coordinates,
    const Directions& directions
) {
    Counts counts{};
    for (std::size_t direction = 0U; direction < directions.size();
         ++direction) {
        for (std::size_t coordinate = 0U; coordinate < counts.size();
             ++coordinate) {
            counts[coordinate] += coordinates[direction]
                * directions[direction][coordinate];
        }
    }
    return counts;
}

std::vector<Counts> transformed_minimal_bases(
    const Directions& directions
) {
    std::vector<Counts> result;
    // Every direction has positive total label at most five.  If an
    // admissible transformed point had total label at least 17, removing
    // any present direction would leave total label at least 12 and hence
    // would still be admissible (the largest possible label is four).
    // Thus every product-order minimal point has total label at most 16,
    // and the coordinate bound 16 is exhaustive.
    constexpr int bound = 16;
    for (int first = 0; first <= bound; ++first) {
        for (int second = 0; second <= bound; ++second) {
            for (int third = 0; third <= bound; ++third) {
                for (int fourth = 0; fourth <= bound; ++fourth) {
                    const Counts transformed{first, second, third, fourth};
                    if (!admissible(apply_directions(
                            transformed, directions
                        ))) {
                        continue;
                    }
                    bool minimal = true;
                    for (std::size_t coordinate = 0U;
                         coordinate < transformed.size(); ++coordinate) {
                        if (transformed[coordinate] == 0) {
                            continue;
                        }
                        Counts predecessor = transformed;
                        --predecessor[coordinate];
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
                }
            }
        }
    }
    return result;
}

cpp_int shifted_middle_template_coefficient(
    const Counts& transformed_base,
    const Counts& degree,
    bool label_three_negative = false
) {
    if (degree[3] != 0) {
        return 0;
    }
    cpp_int result = 0;
    if (degree[2] == 0 && transformed_base[0] == 0
        && transformed_base[1] == 1) {
        const int sum = degree[0] + degree[1];
        if (sum == 6) {
            result = degree[0] * degree[0] - 8 * degree[0] + 26;
        } else if (sum == 7) {
            result = 4 * degree[0] - 7;
        } else if (sum == 8) {
            result = 14;
        }
    } else if (degree[2] == 0) {
        for (int first = degree[0]; first <= 8; ++first) {
            for (int second = degree[1]; second <= 8 - first; ++second) {
                const int sum = first + second;
                cpp_int coefficient = 0;
                if (sum == 6) {
                    coefficient = first * first - 12 * first + 52;
                } else if (sum == 7) {
                    coefficient = 4 * first - 21;
                } else if (sum == 8) {
                    coefficient = 14;
                }
                result += coefficient
                    * binomial_big(
                        transformed_base[0], first - degree[0]
                    )
                    * binomial_big(
                        transformed_base[1], second - degree[1]
                    );
            }
        }
    }
    if (label_three_negative && transformed_base[0] == 0
        && transformed_base[1] == 0) {
        const int ab_degree = degree[0] + degree[1];
        if (degree[2] == 0 && ab_degree == 4) {
            result += 1;
        } else if (degree[2] == 1 && ab_degree == 5) {
            result += degree[0] - 4;
        } else if (degree[2] == 2 && ab_degree == 4) {
            result += 2;
        }
    }
    return result;
}

using UnivariatePolynomial = std::array<cpp_int, 5>;

cpp_int evaluate_newton_polynomial(
    const UnivariatePolynomial& polynomial,
    int argument
) {
    cpp_int result = 0;
    for (std::size_t degree = 0U; degree < polynomial.size(); ++degree) {
        result += polynomial[degree]
            * binomial_big(argument, static_cast<int>(degree));
    }
    return result;
}

cpp_int shifted_newton_coefficient(
    const UnivariatePolynomial& polynomial,
    int shift,
    std::size_t degree
) {
    cpp_int result = 0;
    for (std::size_t source = degree; source < polynomial.size(); ++source) {
        result += polynomial[source]
            * binomial_big(
                shift, static_cast<int>(source - degree)
            );
    }
    return result;
}

UnivariatePolynomial newton_polynomial_from_values(
    std::array<cpp_int, 5> values
) {
    UnivariatePolynomial result{};
    for (std::size_t degree = 0U; degree < result.size(); ++degree) {
        result[degree] = values[0U];
        for (std::size_t index = 0U;
             index + degree + 1U < values.size(); ++index) {
            values[index] = values[index + 1U] - values[index];
        }
    }
    return result;
}

bool certify_nonnegative_integer_tail(
    const UnivariatePolynomial& polynomial,
    int first_argument,
    int& tail_shift
) {
    constexpr int maximum_shift = 512;
    for (int shift = first_argument; shift <= maximum_shift; ++shift) {
        if (evaluate_newton_polynomial(polynomial, shift) < 0) {
            return false;
        }
        bool coefficientwise = true;
        for (std::size_t degree = 0U; degree < polynomial.size(); ++degree) {
            if (shifted_newton_coefficient(
                    polynomial, shift, degree
                ) < 0) {
                coefficientwise = false;
                break;
            }
        }
        if (coefficientwise) {
            tail_shift = shift;
            return true;
        }
    }
    return false;
}

using MonomialPolynomial = std::array<cpp_int, 5>;

MonomialPolynomial multiply_monomial_polynomials(
    const MonomialPolynomial& first,
    const MonomialPolynomial& second
) {
    MonomialPolynomial result{};
    for (std::size_t left = 0U; left < first.size(); ++left) {
        for (std::size_t right = 0U;
             right < second.size() && left + right < result.size(); ++right) {
            result[left + right] += first[left] * second[right];
        }
    }
    return result;
}

cpp_int evaluate_monomial_polynomial(
    const MonomialPolynomial& polynomial,
    int argument
) {
    cpp_int result = 0;
    cpp_int power = 1;
    for (const cpp_int& coefficient : polynomial) {
        result += coefficient * power;
        power *= argument;
    }
    return result;
}

bool certify_middle_discriminants() {
    constexpr std::array<std::array<long long, 5>, 2> b_data{
        std::array<long long, 5>{240, -376, 16, 0, 0},
        std::array<long long, 5>{144, -280, 16, 0, 0}
    };
    constexpr std::array<std::array<long long, 5>, 2> c_data{
        std::array<long long, 5>{0, -322, 347, -26, 1},
        std::array<long long, 5>{0, -170, 187, -18, 1}
    };
    constexpr std::array<std::array<long long, 5>, 2> k_data{
        std::array<long long, 5>{-3600, 1620, 1094, -28, 14},
        std::array<long long, 5>{-1296, -60, 422, 20, 14}
    };
    for (std::size_t family = 0U; family < b_data.size(); ++family) {
        MonomialPolynomial b{};
        MonomialPolynomial c{};
        MonomialPolynomial k{};
        for (std::size_t degree = 0U; degree < b.size(); ++degree) {
            b[degree] = b_data[family][degree];
            c[degree] = c_data[family][degree];
            k[degree] = k_data[family][degree];
        }
        MonomialPolynomial discriminant
            = multiply_monomial_polynomials(b, b);
        for (std::size_t degree = 0U; degree < discriminant.size();
             ++degree) {
            discriminant[degree] -= 480 * c[degree];
            if (discriminant[degree] != -16 * k[degree]) {
                return false;
            }
        }
        std::array<cpp_int, 5> values{};
        for (int offset = 0; offset <= 4; ++offset) {
            values[static_cast<std::size_t>(offset)]
                = evaluate_monomial_polynomial(k, 6 + offset);
        }
        const UnivariatePolynomial shifted
            = newton_polynomial_from_values(values);
        if (shifted[0U] <= 0) {
            return false;
        }
        for (const cpp_int& coefficient : shifted) {
            if (coefficient < 0) {
                return false;
            }
        }
    }
    return true;
}

bool certify_label_three_middle_block() {
    const auto discriminant_margin = [](int first, int second) -> cpp_int {
        const cpp_int n = first + second;
        const cpp_int p = 120 * first * first
            + 8 * (2 * n * n - 47 * n + 30) * first
            + n * (n - 1) * (n * n - 25 * n + 322);
        const cpp_int q = 5 * first * (n - 4) - 4 * n * n + 11 * n;
        return 5 * n * (n - 4) * (n - 5) * p
            - 6 * (n - 1) * q * q + 600 * n * n * (n - 1);
    };
    const auto difference = [&discriminant_margin](
        int first_base, int second_base, int first_degree, int second_degree
    ) {
        cpp_int result = 0;
        for (int first = 0; first <= first_degree; ++first) {
            for (int second = 0; second <= second_degree; ++second) {
                cpp_int term = binomial_big(first_degree, first)
                    * binomial_big(second_degree, second)
                    * discriminant_margin(
                        first_base + first, second_base + second
                    );
                if (((first_degree - first) + (second_degree - second))
                    % 2 != 0) {
                    term = -term;
                }
                result += term;
            }
        }
        return result;
    };
    for (int split = 0; split <= 6; ++split) {
        for (int first_degree = 0; first_degree <= 7; ++first_degree) {
            for (int second_degree = 0;
                 first_degree + second_degree <= 7; ++second_degree) {
                const cpp_int coefficient = difference(
                    split, 6 - split, first_degree, second_degree
                );
                if (coefficient < 0) {
                    std::cout
                        << "LABEL_THREE_MIDDLE_DISCRIMINANT_FAIL base="
                        << split << ',' << 6 - split << " degree="
                        << first_degree << ',' << second_degree
                        << " coefficient=" << coefficient << '\n';
                    return false;
                }
            }
        }
    }
    return true;
}

struct UpperTemplateCertificate {
    std::map<Counts, cpp_int> coefficients;
    int a_tail = -1;
    int negative_linear_threshold = -1;
    int discriminant_tail = -1;
};

struct UpperQuadraticCertificate {
    std::map<Counts, cpp_int> coefficients;
    int nonnegative_label_three_tail = -1;
    int discriminant_tail = -1;
    std::string failure_reason;
};

bool build_upper_template_certificate(
    const Counts& base,
    const Directions& directions,
    const Counts& signs,
    UpperTemplateCertificate& certificate
) {
    const auto coefficient = [&base, &directions, &signs](
        const Counts& degree
    ) {
        return cone_difference_from_polynomial(
            base, directions, degree, signs
        );
    };
    constexpr std::array<int, 5> ratio_scales{210, 42, 14, 6, 3};
    UnivariatePolynomial a_scaled{};
    for (int offset = 0; offset <= 4; ++offset) {
        const Counts degree{4 + offset, 0, 0, 0};
        cpp_int value = coefficient(degree);
        certificate.coefficients[degree] = value;
        if (offset == 0) {
            if (value < 3) {
                return false;
            }
            value -= 3;
        }
        a_scaled[static_cast<std::size_t>(offset)]
            = ratio_scales[static_cast<std::size_t>(offset)] * value;
    }
    UnivariatePolynomial linear_scaled{};
    for (int offset = 0; offset <= 2; ++offset) {
        const Counts degree{4 + offset, 1, 0, 0};
        const cpp_int value = coefficient(degree);
        certificate.coefficients[degree] = value;
        constexpr std::array<int, 3> linear_scales{210, 42, 14};
        linear_scaled[static_cast<std::size_t>(offset)]
            = linear_scales[static_cast<std::size_t>(offset)] * value;
    }
    linear_scaled[0U] -= 140;

    constexpr std::array<Counts, 6> mixed_degrees{
        Counts{4, 2, 0, 0}, Counts{4, 1, 0, 1},
        Counts{4, 1, 1, 0}, Counts{4, 0, 0, 2},
        Counts{4, 0, 1, 1}, Counts{4, 0, 2, 0}
    };
    constexpr std::array<int, 6> mixed_coefficients{2, -1, -1, 2, 2, 2};
    for (std::size_t index = 0U; index < mixed_degrees.size(); ++index) {
        certificate.coefficients[mixed_degrees[index]]
            = mixed_coefficients[index];
        if (coefficient(mixed_degrees[index]) < mixed_coefficients[index]) {
            return false;
        }
    }

    if (!certify_nonnegative_integer_tail(
            a_scaled, 0, certificate.a_tail
        )) {
        return false;
    }
    constexpr int maximum_threshold = 512;
    for (int argument = 0; argument <= maximum_threshold; ++argument) {
        if (evaluate_newton_polynomial(linear_scaled, argument) < 0) {
            certificate.negative_linear_threshold = argument;
            break;
        }
    }
    if (certificate.negative_linear_threshold < 0) {
        return false;
    }
    const int threshold = certificate.negative_linear_threshold;
    if (evaluate_newton_polynomial(linear_scaled, threshold + 1)
        - evaluate_newton_polynomial(linear_scaled, threshold) > 0
        || linear_scaled[2U] > 0) {
        return false;
    }
    std::array<cpp_int, 5> discriminant_values{};
    for (int argument = 0; argument <= 4; ++argument) {
        const cpp_int a_value = evaluate_newton_polynomial(
            a_scaled, argument
        );
        const cpp_int linear_value = evaluate_newton_polynomial(
            linear_scaled, argument
        );
        discriminant_values[static_cast<std::size_t>(argument)]
            = 560 * a_value - linear_value * linear_value;
    }
    const UnivariatePolynomial negative_discriminant
        = newton_polynomial_from_values(discriminant_values);
    if (!certify_nonnegative_integer_tail(
            negative_discriminant, threshold,
            certificate.discriminant_tail
        )) {
        return false;
    }
    return true;
}

bool build_upper_quadratic_certificate(
    const Counts& base,
    const Directions& directions,
    const Counts& signs,
    UpperQuadraticCertificate& certificate
) {
    const auto coefficient = [&base, &directions, &signs](
        const Counts& degree
    ) {
        return cone_difference_from_polynomial(
            base, directions, degree, signs
        );
    };
    constexpr std::array<int, 5> ratio_scales{210, 42, 14, 6, 3};
    UnivariatePolynomial a_scaled{};
    for (int offset = 0; offset <= 4; ++offset) {
        const Counts degree{4 + offset, 0, 0, 0};
        const cpp_int value = coefficient(degree);
        certificate.coefficients[degree] = value;
        a_scaled[static_cast<std::size_t>(offset)]
            = ratio_scales[static_cast<std::size_t>(offset)] * value;
    }
    constexpr std::array<int, 3> linear_scales{210, 42, 14};
    UnivariatePolynomial m_linear_scaled{};
    UnivariatePolynomial u_linear_scaled{};
    for (int offset = 0; offset <= 2; ++offset) {
        const Counts m_degree{4 + offset, 1, 0, 0};
        const Counts u_degree{4 + offset, 0, 1, 0};
        const cpp_int m_value = coefficient(m_degree);
        const cpp_int u_value = coefficient(u_degree);
        certificate.coefficients[m_degree] = m_value;
        certificate.coefficients[u_degree] = u_value;
        const std::size_t position = static_cast<std::size_t>(offset);
        m_linear_scaled[position] = linear_scales[position] * m_value;
        u_linear_scaled[position] = linear_scales[position] * u_value;
    }
    m_linear_scaled[0U] -= 210;

    constexpr std::array<Counts, 6> quadratic_degrees{
        Counts{4, 2, 0, 0}, Counts{4, 0, 2, 0},
        Counts{4, 0, 0, 2}, Counts{4, 0, 1, 1},
        Counts{4, 1, 1, 0}, Counts{4, 1, 0, 1}
    };
    constexpr std::array<int, 6> quadratic_coefficients{2, 2, 2, 2, -1, -1};
    for (std::size_t index = 0U; index < quadratic_degrees.size(); ++index) {
        certificate.coefficients[quadratic_degrees[index]]
            = quadratic_coefficients[index];
        if (coefficient(quadratic_degrees[index])
            < quadratic_coefficients[index]) {
            certificate.failure_reason = "quadratic_resource";
            return false;
        }
    }

    constexpr int maximum_shift = 512;
    for (int argument = 0; argument <= maximum_shift; ++argument) {
        const cpp_int value = evaluate_newton_polynomial(
            u_linear_scaled, argument
        );
        const cpp_int difference
            = evaluate_newton_polynomial(u_linear_scaled, argument + 1)
            - value;
        if (value >= 0 && difference >= 0
            && u_linear_scaled[2U] >= 0) {
            certificate.nonnegative_label_three_tail = argument;
            break;
        }
    }
    if (certificate.nonnegative_label_three_tail < 0) {
        certificate.failure_reason = "label_three_tail";
        return false;
    }

    const int u_tail = certificate.nonnegative_label_three_tail;
    std::array<cpp_int, 5> discriminant_values{};
    for (int argument = 0; argument <= 4; ++argument) {
        const cpp_int a_value = evaluate_newton_polynomial(
            a_scaled, argument
        );
        const cpp_int p_value = evaluate_newton_polynomial(
            m_linear_scaled, argument
        );
        constexpr int q_value = -210;
        discriminant_values[static_cast<std::size_t>(argument)]
            = 630 * a_value - p_value * p_value
            - p_value * q_value - q_value * q_value;
    }
    const UnivariatePolynomial discriminant_margin
        = newton_polynomial_from_values(discriminant_values);
    for (int shift = u_tail; shift <= maximum_shift; ++shift) {
        bool coefficientwise = true;
        for (std::size_t degree = 0U;
             degree < discriminant_margin.size(); ++degree) {
            if (shifted_newton_coefficient(
                    discriminant_margin, shift, degree
                ) < 0) {
                coefficientwise = false;
                break;
            }
        }
        if (coefficientwise) {
            certificate.discriminant_tail = shift;
            break;
        }
    }
    if (certificate.discriminant_tail < 0) {
        certificate.failure_reason = "discriminant_tail";
        return false;
    }
    for (int argument = 0; argument < certificate.discriminant_tail;
         ++argument) {
        const cpp_int a_value = evaluate_newton_polynomial(
            a_scaled, argument
        );
        const cpp_int p_value = evaluate_newton_polynomial(
            m_linear_scaled, argument
        );
        const cpp_int g_value = evaluate_newton_polynomial(
            u_linear_scaled, argument
        );
        const cpp_int q_value = -210 + std::min(g_value, cpp_int{0});
        const cpp_int absolute_p = p_value < 0 ? -p_value : p_value;
        const cpp_int absolute_q = q_value < 0 ? -q_value : q_value;
        const cpp_int bound_big = (absolute_p + absolute_q) / 105 + 2;
        if (bound_big > 5000) {
            certificate.failure_reason = "integer_minimum_bound";
            return false;
        }
        const int bound = bound_big.convert_to<int>();
        cpp_int minimum = a_value;
        for (int m = 0; m <= bound; ++m) {
            for (int s = 0; s <= bound; ++s) {
                const cpp_int value = a_value
                    + cpp_int{210} * (m * m + s * s - m * s)
                    + p_value * m + q_value * s;
                minimum = std::min(minimum, value);
            }
        }
        if (minimum < 0) {
            certificate.failure_reason = "integer_minimum_at_"
                + std::to_string(argument);
            return false;
        }
    }
    return true;
}

void print_counts(const Counts& counts) {
    std::cout << counts[0] << ',' << counts[1] << ',' << counts[2] << ','
              << counts[3];
}

struct RecursiveStats {
    std::size_t nodes = 0U;
    std::size_t positive_leaves = 0U;
    std::size_t split_nodes = 0U;
    std::size_t maximum_depth = 0U;
    Counts maximum_base{};
};

bool newton_node_is_nonnegative(
    const Counts& base,
    unsigned int free_mask,
    const Counts& signs,
    Counts& first_negative_degree
) {
    constexpr Counts degrees{8, 4, 4, 4};
    for (int first = 0;
         first <= (((free_mask & 1U) != 0U) ? degrees[0] : 0);
         ++first) {
        for (int second = 0;
             second <= (((free_mask & 2U) != 0U) ? degrees[1] : 0);
             ++second) {
            for (int third = 0;
                 third <= (((free_mask & 4U) != 0U) ? degrees[2] : 0);
                 ++third) {
                for (int fourth = 0;
                     fourth <= (((free_mask & 8U) != 0U) ? degrees[3] : 0);
                     ++fourth) {
                    const Counts degree{first, second, third, fourth};
                    if (finite_difference(base, degree, signs) < 0) {
                        first_negative_degree = degree;
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

bool recursive_newton_certificate(
    const Counts& base,
    unsigned int free_mask,
    const Counts& signs,
    std::size_t depth,
    RecursiveStats& stats,
    std::map<RecursiveKey, bool>& memo
) {
    const RecursiveKey key{
        base[0], base[1], base[2], base[3],
        static_cast<int>(free_mask)
    };
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
    if (newton_node_is_nonnegative(
            base, free_mask, signs, negative_degree
        )) {
        ++stats.positive_leaves;
        memo.emplace(key, true);
        return true;
    }
    if (free_mask == 0U || depth >= 80U || stats.nodes >= 20000U) {
        memo.emplace(key, false);
        return false;
    }
    std::array<std::size_t, 4> candidates{0U, 1U, 2U, 3U};
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
        const unsigned int slice_mask = free_mask & ~bit;
        if (!recursive_newton_certificate(
                base, slice_mask, signs, depth + 1U, stats, memo
            )) {
            continue;
        }
        Counts shifted = base;
        ++shifted[coordinate];
        if (recursive_newton_certificate(
                shifted, free_mask, signs, depth + 1U, stats, memo
            )) {
            ++stats.split_nodes;
            memo[key] = true;
            return true;
        }
    }
    memo[key] = false;
    return false;
}

std::size_t grid_index(
    const Counts& coordinate,
    const Counts& bounds
) {
    std::size_t index = 0U;
    std::size_t stride = 1U;
    for (std::size_t axis = 0U; axis < coordinate.size(); ++axis) {
        index += static_cast<std::size_t>(coordinate[axis]) * stride;
        stride *= static_cast<std::size_t>(bounds[axis] + 1);
    }
    return index;
}

bool cone_newton_is_nonnegative(
    const Counts& base,
    const Directions& directions,
    const Counts& degree_bounds,
    const Counts& signs,
    Counts& first_negative_degree,
    const Counts* expansion_base = nullptr
) {
    std::size_t grid_size = 1U;
    for (int bound : degree_bounds) {
        grid_size *= static_cast<std::size_t>(bound + 1);
    }
    std::vector<cpp_int> coefficients(grid_size);
    const Counts& polynomial_base = expansion_base == nullptr
        ? base : *expansion_base;
    const AxisPolynomial& polynomial = axis_polynomial(
        polynomial_base, signs
    );
    audit_axis_polynomial(polynomial_base, signs, polynomial);
    Counts coordinate{};
    for (std::size_t linear = 0U; linear < grid_size; ++linear) {
        std::size_t residual = linear;
        Counts point = base;
        for (std::size_t axis = 0U; axis < coordinate.size(); ++axis) {
            const std::size_t width = static_cast<std::size_t>(
                degree_bounds[axis] + 1
            );
            coordinate[axis] = static_cast<int>(residual % width);
            residual /= width;
            for (std::size_t original = 0U; original < point.size();
                 ++original) {
                point[original] += coordinate[axis]
                    * directions[axis][original];
            }
        }
        coefficients[linear] = evaluate_axis_polynomial(
            polynomial_base, point, polynomial
        );
    }
    for (std::size_t axis = 0U; axis < coordinate.size(); ++axis) {
        std::vector<cpp_int> transformed(grid_size);
        for (std::size_t linear = 0U; linear < grid_size; ++linear) {
            std::size_t residual = linear;
            for (std::size_t decode = 0U; decode < coordinate.size();
                 ++decode) {
                const std::size_t width = static_cast<std::size_t>(
                    degree_bounds[decode] + 1
                );
                coordinate[decode] = static_cast<int>(residual % width);
                residual /= width;
            }
            cpp_int value = 0;
            const int order = coordinate[axis];
            for (int sample = 0; sample <= order; ++sample) {
                Counts source = coordinate;
                source[axis] = sample;
                const cpp_int term = binomial(order, sample)
                    * coefficients[grid_index(source, degree_bounds)];
                value += ((order - sample) & 1) != 0 ? -term : term;
            }
            transformed[linear] = value;
        }
        coefficients = std::move(transformed);
    }
    for (std::size_t linear = 0U; linear < grid_size; ++linear) {
        if (coefficients[linear] >= 0) {
            continue;
        }
        std::size_t residual = linear;
        for (std::size_t axis = 0U; axis < coordinate.size(); ++axis) {
            const std::size_t width = static_cast<std::size_t>(
                degree_bounds[axis] + 1
            );
            first_negative_degree[axis] = static_cast<int>(residual % width);
            residual /= width;
        }
        return false;
    }
    return true;
}

bool two_coordinate_fan_certificate(
    const Counts& base,
    std::size_t first,
    std::size_t second,
    const Counts& signs,
    Counts& first_failure,
    Counts& second_failure
) {
    constexpr Counts original_degrees{8, 4, 4, 4};
    std::array<std::size_t, 2> remaining{};
    std::size_t retained = 0U;
    for (std::size_t coordinate = 0U; coordinate < 4U; ++coordinate) {
        if (coordinate != first && coordinate != second) {
            remaining[retained] = coordinate;
            ++retained;
        }
    }
    Directions first_directions{};
    first_directions[0U][first] = 1;
    first_directions[0U][second] = 1;
    first_directions[1U][first] = 1;
    first_directions[2U][remaining[0U]] = 1;
    first_directions[3U][remaining[1U]] = 1;
    const Counts first_bounds{
        original_degrees[first] + original_degrees[second],
        original_degrees[first],
        original_degrees[remaining[0U]],
        original_degrees[remaining[1U]]
    };
    if (!cone_newton_is_nonnegative(
            base, first_directions, first_bounds, signs, first_failure
        )) {
        return false;
    }
    Counts shifted = base;
    ++shifted[second];
    Directions second_directions = first_directions;
    second_directions[1U][first] = 0;
    second_directions[1U][second] = 1;
    const Counts second_bounds{
        original_degrees[first] + original_degrees[second],
        original_degrees[second],
        original_degrees[remaining[0U]],
        original_degrees[remaining[1U]]
    };
    return cone_newton_is_nonnegative(
        shifted, second_directions, second_bounds, signs, second_failure
    );
}

Counts direction_degree_bounds(const Directions& directions) {
    constexpr Counts original_degrees{8, 4, 4, 4};
    constexpr int total_degree = 8;
    Counts bounds{};
    for (std::size_t direction = 0U; direction < directions.size();
         ++direction) {
        for (std::size_t coordinate = 0U;
             coordinate < directions[direction].size(); ++coordinate) {
            if (directions[direction][coordinate] != 0) {
                bounds[direction] += original_degrees[coordinate];
            }
        }
        bounds[direction] = std::min(bounds[direction], total_degree);
    }
    return bounds;
}

struct FanStats {
    std::size_t nodes = 0U;
    std::size_t leaves = 0U;
    std::size_t splits = 0U;
    std::size_t maximum_depth = 0U;
    std::array<int, 2> failure_left{};
    std::array<int, 2> failure_right{};
    Counts failure_degree{};
};

bool adaptive_slope_certificate(
    const Counts& base,
    const std::array<int, 2>& left,
    const std::array<int, 2>& right,
    const Counts& signs,
    std::size_t depth,
    FanStats& stats
) {
    ++stats.nodes;
    stats.maximum_depth = std::max(stats.maximum_depth, depth);
    Directions directions{};
    directions[0U][0U] = left[0U];
    directions[0U][1U] = left[1U];
    directions[1U][0U] = right[0U];
    directions[1U][1U] = right[1U];
    directions[2U][2U] = 1;
    directions[3U][3U] = 1;
    Counts failure{};
    if (cone_newton_is_nonnegative(
            base, directions, direction_degree_bounds(directions), signs,
            failure
        )) {
        ++stats.leaves;
        return true;
    }
    if (depth >= 24U || stats.nodes >= 2000U) {
        stats.failure_left = left;
        stats.failure_right = right;
        stats.failure_degree = failure;
        return false;
    }
    const std::array<int, 2> middle{
        left[0U] + right[0U], left[1U] + right[1U]
    };
    if (!adaptive_slope_certificate(
            base, left, middle, signs, depth + 1U, stats
        )) {
        return false;
    }
    if (!adaptive_slope_certificate(
            base, middle, right, signs, depth + 1U, stats
        )) {
        return false;
    }
    ++stats.splits;
    return true;
}

bool recursive_fan_certificate(
    const Counts& base,
    const Directions& directions,
    const Counts& signs,
    std::size_t depth,
    FanStats& stats
) {
    ++stats.nodes;
    stats.maximum_depth = std::max(stats.maximum_depth, depth);
    Counts negative_degree{};
    if (cone_newton_is_nonnegative(
            base, directions, direction_degree_bounds(directions), signs,
            negative_degree
        )) {
        ++stats.leaves;
        return true;
    }
    if (depth >= 10U || stats.nodes >= 2000U) {
        return false;
    }
    std::array<std::pair<std::size_t, std::size_t>, 6> pairs{};
    std::size_t pair_count = 0U;
    for (std::size_t first = 0U; first < 4U; ++first) {
        for (std::size_t second = first + 1U; second < 4U; ++second) {
            pairs[pair_count] = {first, second};
            ++pair_count;
        }
    }
    std::stable_sort(
        pairs.begin(), pairs.end(),
        [&negative_degree](const auto& first, const auto& second) {
            const int first_score = negative_degree[first.first]
                * negative_degree[first.second];
            const int second_score = negative_degree[second.first]
                * negative_degree[second.second];
            return first_score > second_score;
        }
    );
    for (const auto& [first, second] : pairs) {
        Directions first_child{};
        for (std::size_t coordinate = 0U; coordinate < 4U; ++coordinate) {
            first_child[0U][coordinate] = directions[first][coordinate]
                + directions[second][coordinate];
            first_child[1U][coordinate] = directions[first][coordinate];
        }
        std::size_t retained = 2U;
        for (std::size_t direction = 0U; direction < 4U; ++direction) {
            if (direction != first && direction != second) {
                first_child[retained] = directions[direction];
                ++retained;
            }
        }
        if (!recursive_fan_certificate(
                base, first_child, signs, depth + 1U, stats
            )) {
            continue;
        }
        Counts second_base = base;
        for (std::size_t coordinate = 0U; coordinate < 4U; ++coordinate) {
            second_base[coordinate] += directions[second][coordinate];
        }
        Directions second_child = first_child;
        second_child[1U] = directions[second];
        if (recursive_fan_certificate(
                second_base, second_child, signs, depth + 1U, stats
            )) {
            ++stats.splits;
            return true;
        }
    }
    return false;
}

bool recursive_cone_shift_certificate(
    const Counts& base,
    const Directions& directions,
    unsigned int free_mask,
    const Counts& signs,
    std::size_t depth,
    FanStats& stats
) {
    ++stats.nodes;
    stats.maximum_depth = std::max(stats.maximum_depth, depth);
    Directions active = directions;
    for (std::size_t direction = 0U; direction < active.size(); ++direction) {
        if ((free_mask & (1U << direction)) == 0U) {
            active[direction] = Counts{};
        }
    }
    Counts negative_degree{};
    if (cone_newton_is_nonnegative(
            base, active, direction_degree_bounds(active), signs,
            negative_degree
        )) {
        ++stats.leaves;
        return true;
    }
    if (free_mask == 0U || depth >= 40U || stats.nodes >= 5000U) {
        stats.failure_degree = negative_degree;
        return false;
    }
    std::array<std::size_t, 4> candidates{0U, 1U, 2U, 3U};
    std::stable_sort(
        candidates.begin(), candidates.end(),
        [&negative_degree](std::size_t first, std::size_t second) {
            return negative_degree[first] > negative_degree[second];
        }
    );
    for (std::size_t direction : candidates) {
        const unsigned int bit = 1U << direction;
        if ((free_mask & bit) == 0U) {
            continue;
        }
        if (!recursive_cone_shift_certificate(
                base, directions, free_mask & ~bit, signs,
                depth + 1U, stats
            )) {
            continue;
        }
        Counts shifted = base;
        for (std::size_t coordinate = 0U; coordinate < shifted.size();
             ++coordinate) {
            shifted[coordinate] += directions[direction][coordinate];
        }
        if (recursive_cone_shift_certificate(
                shifted, directions, free_mask, signs, depth + 1U, stats
            )) {
            ++stats.splits;
            return true;
        }
    }
    stats.failure_degree = negative_degree;
    return false;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 4) {
            throw std::runtime_error(
                "usage: analyze_su2_outer_eight MAXIMUM_TOTAL_COUNT "
                "[both|search|certificate|recursive|fan|recursive-fan|"
                "slope-fan|adaptive-slope|shifted-slope|ray-probe|"
                "pattern4-chambers|deficit8-chambers] "
                "[PATTERN]"
            );
        }
        const int maximum_count = std::stoi(argv[1]);
        if (maximum_count < 1 || maximum_count > 40) {
            throw std::runtime_error("invalid maximum count");
        }
        const std::string mode = argc >= 3 ? argv[2] : "both";
        const int requested_pattern = argc == 4 ? std::stoi(argv[3]) : -1;
        if (mode != "both" && mode != "search" && mode != "certificate"
            && mode != "recursive" && mode != "fan") {
            if (mode != "recursive-fan" && mode != "slope-fan"
                && mode != "adaptive-slope" && mode != "shifted-slope") {
                if (mode != "ray-probe" && mode != "pattern4-chambers"
                    && mode != "deficit8-chambers") {
                    throw std::runtime_error("invalid mode");
                }
            }
        }
        const bool pattern_mode = mode == "recursive" || mode == "fan"
            || mode == "recursive-fan" || mode == "slope-fan"
            || mode == "adaptive-slope" || mode == "shifted-slope";
        const bool selected_pattern_mode = pattern_mode || mode == "ray-probe"
            || mode == "pattern4-chambers" || mode == "deficit8-chambers";
        if (argc == 4
            && (!selected_pattern_mode || requested_pattern < 0
                || requested_pattern >= 8)) {
            throw std::runtime_error("invalid recursive pattern");
        }
        const auto bases = minimal_bases();
        std::cout << "minimal_bases=" << bases.size() << '\n';
        if (mode != "certificate" && !selected_pattern_mode) {
                cpp_int minimum = 0;
                Counts minimum_counts{};
                int minimum_pattern = -1;
                long long tested = 0;
                bool initialized = false;
                for (int first = 0; first <= maximum_count; ++first) {
                    for (int second = 0;
                         first + second <= maximum_count; ++second) {
                        for (int third = 0;
                             first + second + third <= maximum_count;
                             ++third) {
                            for (int fourth = 0;
                                 first + second + third + fourth
                                     <= maximum_count;
                                 ++fourth) {
                                const Counts counts{
                                    first, second, third, fourth
                                };
                                if (!admissible(counts)) {
                                    continue;
                                }
                                for (int pattern = 0; pattern < 8;
                                     ++pattern) {
                                    const Counts signs{
                                        1,
                                        (pattern & 4) != 0 ? -1 : 1,
                                        (pattern & 2) != 0 ? -1 : 1,
                                        (pattern & 1) != 0 ? -1 : 1
                                    };
                                    const cpp_int value = outer_value(
                                        counts, signs
                                    );
                                    ++tested;
                                    if (!initialized || value < minimum) {
                                        initialized = true;
                                        minimum = value;
                                        minimum_counts = counts;
                                        minimum_pattern = pattern;
                                    }
                                    if (value < 0) {
                                        std::cout << "FAIL value=" << value
                                                  << " pattern=" << pattern
                                                  << " counts=";
                                        print_counts(counts);
                                        std::cout << '\n';
                                        return EXIT_FAILURE;
                                    }
                                }
                            }
                        }
                    }
                }
                std::cout << "SU2_OUTER_EIGHT_SEARCH PASS tested=" << tested
                          << " maximum_count=" << maximum_count
                          << " minimum=" << minimum
                          << " pattern=" << minimum_pattern << " counts=";
                print_counts(minimum_counts);
                std::cout << '\n';
        }
        if (mode == "search") {
            return EXIT_SUCCESS;
        }
        outer_cache().clear();
        multiplicity_cache().clear();

        if (mode == "pattern4-chambers" || mode == "deficit8-chambers") {
            const int chamber_pattern = mode == "pattern4-chambers"
                ? 4 : requested_pattern;
            const Counts signs{
                1,
                (chamber_pattern & 4) != 0 ? -1 : 1,
                (chamber_pattern & 2) != 0 ? -1 : 1,
                (chamber_pattern & 1) != 0 ? -1 : 1
            };
            const std::string tag = "PATTERN"
                + std::to_string(chamber_pattern);
            const bool label_three_negative = (chamber_pattern & 2) != 0;
            if (!certify_middle_discriminants()) {
                throw std::runtime_error(
                    "pattern-four middle discriminant audit failed"
                );
            }
            if (label_three_negative
                && !certify_label_three_middle_block()) {
                throw std::runtime_error(
                    "negative-label-three middle block audit failed"
                );
            }
            std::cout << tag << "_TEMPLATE_SELF_AUDIT middle_boundary="
                      << shifted_middle_template_coefficient(
                             Counts{0, 1, 0, 2}, Counts{0, 6, 0, 0}
                         )
                      << " middle_discriminants=PASS"
                      << (label_three_negative
                              ? " label_three_middle_block=PASS" : "")
                      << '\n';
            std::array<Directions, 3> chambers{};
            chambers[0U][0U] = Counts{1, 0, 0, 0};
            chambers[0U][1U] = Counts{1, 1, 0, 0};
            chambers[1U][0U] = Counts{1, 1, 0, 0};
            chambers[1U][1U] = Counts{1, 2, 0, 0};
            chambers[2U][0U] = Counts{1, 2, 0, 0};
            chambers[2U][1U] = Counts{0, 1, 0, 0};
            constexpr std::array<const char*, 3> names{
                "lower", "middle", "upper"
            };
            for (Directions& chamber : chambers) {
                chamber[2U][2U] = 1;
                chamber[3U][3U] = 1;
            }
            std::array<std::size_t, 3> passes{};
            std::array<std::size_t, 3> failures{};
            std::array<std::size_t, 3> analytic_exceptions{};
            bool certificate_passes = true;
            for (std::size_t chamber = 0U; chamber < chambers.size();
                 ++chamber) {
                const std::vector<Counts> transformed_bases
                    = transformed_minimal_bases(chambers[chamber]);
                std::cout << tag << "_CHAMBER_BASES chamber="
                          << names[chamber] << " count="
                          << transformed_bases.size() << '\n';
                for (const Counts& transformed_base : transformed_bases) {
                    const Counts base = apply_directions(
                        transformed_base, chambers[chamber]
                    );
                    Counts failure{};
                    const bool pass = cone_newton_is_nonnegative(
                        base, chambers[chamber],
                        direction_degree_bounds(chambers[chamber]), signs,
                        failure
                    );
                    if (pass) {
                        ++passes[chamber];
                    } else {
                        ++failures[chamber];
                        std::cout
                            << tag << "_CHAMBER_EXCEPTION transformed_base=";
                        print_counts(transformed_base);
                        std::cout << " base=";
                        print_counts(base);
                        std::cout << " chamber=" << names[chamber]
                                  << " degree=";
                        print_counts(failure);
                        std::cout << '\n';
                        bool dominates_template = chamber != 0U;
                        UpperTemplateCertificate upper_certificate;
                        UpperQuadraticCertificate upper_quadratic_certificate;
                        bool upper_analytic_pass = true;
                        if (chamber == 2U) {
                            upper_analytic_pass = label_three_negative
                                ? build_upper_quadratic_certificate(
                                    base, chambers[chamber], signs,
                                    upper_quadratic_certificate
                                )
                                : build_upper_template_certificate(
                                    base, chambers[chamber], signs,
                                    upper_certificate
                                );
                        }
                        if (chamber == 2U && !upper_analytic_pass) {
                            dominates_template = false;
                            std::cout << tag << "_UPPER_ANALYTIC_FAIL "
                                      << "transformed_base=";
                            print_counts(transformed_base);
                            if (label_three_negative) {
                                std::cout << " reason="
                                          << upper_quadratic_certificate
                                                 .failure_reason;
                            }
                            std::cout << '\n';
                        }
                        const Counts bounds
                            = direction_degree_bounds(chambers[chamber]);
                        for (int first = 0; first <= bounds[0]; ++first) {
                            for (int second = 0; second <= bounds[1];
                                 ++second) {
                                for (int third = 0; third <= bounds[2];
                                     ++third) {
                                    for (int fourth = 0; fourth <= bounds[3];
                                         ++fourth) {
                                        const int degree_sum = first + second
                                            + third + fourth;
                                        if (degree_sum > 8) {
                                            continue;
                                        }
                                        const Counts degree{
                                            first, second, third, fourth
                                        };
                                        const cpp_int coefficient
                                            = cone_difference_from_polynomial(
                                                base, chambers[chamber],
                                                degree, signs
                                            );
                                        cpp_int template_coefficient = 0;
                                        if (chamber == 1U) {
                                            template_coefficient
                                                = shifted_middle_template_coefficient(
                                                    transformed_base, degree,
                                                    label_three_negative
                                                );
                                        } else if (chamber == 2U) {
                                            if (label_three_negative) {
                                                const auto found
                                                    = upper_quadratic_certificate
                                                          .coefficients.find(degree);
                                                if (found
                                                    != upper_quadratic_certificate
                                                           .coefficients.end()) {
                                                    template_coefficient
                                                        = found->second;
                                                }
                                            } else {
                                                const auto found
                                                    = upper_certificate
                                                          .coefficients.find(degree);
                                                if (found
                                                    != upper_certificate
                                                           .coefficients.end()) {
                                                    template_coefficient
                                                        = found->second;
                                                }
                                            }
                                        }
                                        if (coefficient < template_coefficient) {
                                            dominates_template = false;
                                            std::cout
                                                << tag << "_TEMPLATE_FAIL "
                                                << "transformed_base=";
                                            print_counts(transformed_base);
                                            std::cout << " chamber="
                                                      << names[chamber]
                                                      << " degree=";
                                            print_counts(degree);
                                            std::cout << " coefficient="
                                                      << coefficient
                                                      << " template="
                                                      << template_coefficient
                                                      << '\n';
                                        }
                                        if (std::getenv(
                                                "OUTER_CHAMBER_VERBOSE"
                                            ) != nullptr
                                            && coefficient < 0) {
                                            std::cout
                                                << tag << "_CHAMBER_NEGATIVE "
                                                << "transformed_base=";
                                            print_counts(transformed_base);
                                            std::cout << " chamber="
                                                      << names[chamber]
                                                      << " degree=";
                                            print_counts(degree);
                                            std::cout << " value="
                                                      << coefficient << '\n';
                                        }
                                    }
                                }
                            }
                        }
                        if (dominates_template) {
                            ++analytic_exceptions[chamber];
                            if (chamber == 2U) {
                                std::cout
                                    << tag << "_UPPER_ANALYTIC_PASS "
                                    << "transformed_base=";
                                print_counts(transformed_base);
                                if (label_three_negative) {
                                    std::cout
                                        << " label_three_tail="
                                        << upper_quadratic_certificate
                                               .nonnegative_label_three_tail
                                        << " discriminant_tail="
                                        << upper_quadratic_certificate
                                               .discriminant_tail;
                                } else {
                                    std::cout
                                        << " a_tail="
                                        << upper_certificate.a_tail
                                        << " negative_linear_threshold="
                                        << upper_certificate
                                               .negative_linear_threshold
                                        << " discriminant_tail="
                                        << upper_certificate
                                               .discriminant_tail;
                                }
                                std::cout << '\n';
                            }
                        } else {
                            certificate_passes = false;
                        }
                    }
                }
                std::cout << tag << "_CHAMBER_SUMMARY chamber="
                          << names[chamber] << " passes=" << passes[chamber]
                          << " exceptions=" << failures[chamber]
                          << " analytic_exceptions="
                          << analytic_exceptions[chamber] << '\n';
            }
            std::cout << tag << "_GLOBAL_CHAMBER_CERTIFICATE "
                      << (certificate_passes ? "PASS" : "FAIL") << '\n';
            return certificate_passes ? EXIT_SUCCESS : EXIT_FAILURE;
        }

        if (mode == "ray-probe") {
            const Counts base{0, 0, 0, 3};
            const Counts direction{1, 2, 0, 0};
            for (int pattern = 0; pattern < 8; ++pattern) {
                if (requested_pattern >= 0 && pattern != requested_pattern) {
                    continue;
                }
                const Counts signs{
                    1,
                    (pattern & 4) != 0 ? -1 : 1,
                    (pattern & 2) != 0 ? -1 : 1,
                    (pattern & 1) != 0 ? -1 : 1
                };
                std::cout << "RAY_PROBE pattern=" << pattern << " base=";
                print_counts(base);
                std::cout << " direction=";
                print_counts(direction);
                for (int degree = 0; degree <= 8; ++degree) {
                    std::cout << " delta" << degree << '='
                              << directional_difference_from_polynomial(
                                     base, direction, 0, degree, signs
                                 );
                }
                const cpp_int delta7_at_zero
                    = directional_difference_from_polynomial(
                        base, direction, 0, 7, signs
                    );
                const cpp_int delta8
                    = directional_difference_from_polynomial(
                        base, direction, 0, 8, signs
                    );
                std::cout << " delta7_shift1="
                          << directional_difference_from_polynomial(
                                 base, direction, 1, 7, signs
                             );
                if (delta7_at_zero < 0 && delta8 > 0) {
                    const cpp_int threshold = (-delta7_at_zero + delta8 - 1)
                        / delta8;
                    std::cout << " affine_nonnegative_threshold="
                              << threshold;
                }
                Directions boundary_cone{};
                boundary_cone[0U] = direction;
                boundary_cone[1U][1U] = 1;
                boundary_cone[2U][2U] = 1;
                boundary_cone[3U][3U] = 1;
                for (int shift = 0; shift <= 2; ++shift) {
                    Counts shifted = base;
                    for (std::size_t coordinate = 0U;
                         coordinate < shifted.size(); ++coordinate) {
                        shifted[coordinate] += shift * direction[coordinate];
                    }
                    Directions tested = boundary_cone;
                    const char* kind = "tail";
                    if (shift < 2) {
                        tested[0U] = Counts{};
                        kind = "slice";
                    }
                    Counts failure{};
                    const bool pass = cone_newton_is_nonnegative(
                        shifted, tested, direction_degree_bounds(tested),
                        signs, failure, &base
                    );
                    std::cout << ' ' << kind << shift << '='
                              << (pass ? "PASS" : "FAIL");
                    if (!pass) {
                        std::cout << ':';
                        print_counts(failure);
                    }
                }
                int tail_threshold = -1;
                bool slices_pass = true;
                Counts tail_failure{};
                for (int shift = 0; shift <= 64; ++shift) {
                    Counts shifted = base;
                    for (std::size_t coordinate = 0U;
                         coordinate < shifted.size(); ++coordinate) {
                        shifted[coordinate] += shift * direction[coordinate];
                    }
                    if (tail_threshold < 0
                        && cone_newton_is_nonnegative(
                            shifted, boundary_cone,
                            direction_degree_bounds(boundary_cone), signs,
                            tail_failure, &base
                        )) {
                        tail_threshold = shift;
                    }
                    if (tail_threshold < 0) {
                        Directions slice = boundary_cone;
                        slice[0U] = Counts{};
                        Counts slice_failure{};
                        if (!cone_newton_is_nonnegative(
                                shifted, slice,
                                direction_degree_bounds(slice), signs,
                                slice_failure, &base
                            )) {
                            slices_pass = false;
                        }
                    }
                    if (tail_threshold >= 0) {
                        break;
                    }
                }
                std::cout << " full_tail_threshold=" << tail_threshold
                          << " preceding_slices="
                          << (slices_pass ? "PASS" : "FAIL");
                if (tail_threshold < 0) {
                    std::cout << " last_tail_degree=";
                    print_counts(tail_failure);
                }
                std::cout << '\n';
                std::size_t negative_coefficients = 0U;
                for (int first = 0; first <= 8; ++first) {
                    for (int second = 0; second <= 4; ++second) {
                        if (first + second > 8) {
                            continue;
                        }
                        const Counts degree{first, second, 0, 0};
                        const cpp_int coefficient
                            = cone_difference_from_polynomial(
                                base, boundary_cone, degree, signs
                            );
                        if (coefficient != 0) {
                            std::cout << "RAY_BIVARIATE_COEFFICIENT pattern="
                                      << pattern << " degree=" << first
                                      << ',' << second << " value="
                                      << coefficient << '\n';
                        }
                    }
                }
                for (int first = 0; first <= 8; ++first) {
                    for (int second = 0; second <= 4; ++second) {
                        for (int third = 0; third <= 4; ++third) {
                            for (int fourth = 0; fourth <= 4; ++fourth) {
                                if (first + second + third + fourth > 8) {
                                    continue;
                                }
                                const Counts degree{
                                    first, second, third, fourth
                                };
                                const cpp_int coefficient
                                    = cone_difference_from_polynomial(
                                        base, boundary_cone, degree, signs
                                    );
                                if (coefficient != 0
                                    && first + second + third + fourth >= 6) {
                                    std::cout
                                        << "RAY_HIGH_COEFFICIENT pattern="
                                        << pattern << " degree=";
                                    print_counts(degree);
                                    std::cout << " value=" << coefficient
                                              << '\n';
                                }
                                if (coefficient < 0) {
                                    ++negative_coefficients;
                                    std::cout
                                        << "RAY_NEGATIVE_COEFFICIENT pattern="
                                        << pattern << " degree=";
                                    print_counts(degree);
                                    std::cout << " value=" << coefficient
                                              << '\n';
                                }
                            }
                        }
                    }
                }
                std::cout << "RAY_NEGATIVE_SUMMARY pattern=" << pattern
                          << " count=" << negative_coefficients << '\n';
                const auto print_chamber_negatives = [
                    &base, &signs, pattern
                ](const char* name, const Directions& chamber) {
                    std::size_t count = 0U;
                    const bool verbose
                        = std::getenv("OUTER_CHAMBER_VERBOSE") != nullptr;
                    for (int first = 0; first <= 8; ++first) {
                        for (int second = 0; second <= 8; ++second) {
                            for (int third = 0; third <= 4; ++third) {
                                for (int fourth = 0; fourth <= 4; ++fourth) {
                                    if (first + second + third + fourth > 8) {
                                        continue;
                                    }
                                    const Counts degree{
                                        first, second, third, fourth
                                    };
                                    const cpp_int coefficient
                                        = cone_difference_from_polynomial(
                                            base, chamber, degree, signs
                                        );
                                    if (verbose && coefficient != 0
                                        && first + second + third + fourth
                                            >= 6) {
                                        std::cout
                                            << "CHAMBER_HIGH pattern="
                                            << pattern << " chamber=" << name
                                            << " degree=";
                                        print_counts(degree);
                                        std::cout << " value=" << coefficient
                                                  << '\n';
                                    }
                                    if (coefficient < 0) {
                                        ++count;
                                        std::cout
                                            << "CHAMBER_NEGATIVE pattern="
                                            << pattern << " chamber=" << name
                                            << " degree=";
                                        print_counts(degree);
                                        std::cout << " value=" << coefficient
                                                  << '\n';
                                    }
                                }
                            }
                        }
                    }
                    std::cout << "CHAMBER_SUMMARY pattern=" << pattern
                              << " chamber=" << name << " negatives="
                              << count << '\n';
                };
                Directions lower_chamber{};
                lower_chamber[0U] = Counts{1, 0, 0, 0};
                lower_chamber[1U] = Counts{1, 1, 0, 0};
                lower_chamber[2U][2U] = 1;
                lower_chamber[3U][3U] = 1;
                print_chamber_negatives("lower", lower_chamber);
                Directions middle_chamber{};
                middle_chamber[0U] = Counts{1, 1, 0, 0};
                middle_chamber[1U] = Counts{1, 2, 0, 0};
                middle_chamber[2U][2U] = 1;
                middle_chamber[3U][3U] = 1;
                print_chamber_negatives("middle", middle_chamber);
            }
            return EXIT_SUCCESS;
        }

        if (mode == "shifted-slope") {
            constexpr std::array<std::array<int, 2>, 4> rays{
                std::array<int, 2>{1, 0},
                std::array<int, 2>{1, 1},
                std::array<int, 2>{1, 2},
                std::array<int, 2>{0, 1}
            };
            for (int pattern = 0; pattern < 8; ++pattern) {
                if (requested_pattern >= 0 && pattern != requested_pattern) {
                    continue;
                }
                const Counts signs{
                    1,
                    (pattern & 4) != 0 ? -1 : 1,
                    (pattern & 2) != 0 ? -1 : 1,
                    (pattern & 1) != 0 ? -1 : 1
                };
                FanStats stats;
                bool pass = true;
                for (const Counts& base : bases) {
                    for (std::size_t ray = 0U; ray + 1U < rays.size();
                         ++ray) {
                        Directions directions{};
                        directions[0U][0U] = rays[ray][0U];
                        directions[0U][1U] = rays[ray][1U];
                        directions[1U][0U] = rays[ray + 1U][0U];
                        directions[1U][1U] = rays[ray + 1U][1U];
                        directions[2U][2U] = 1;
                        directions[3U][3U] = 1;
                        if (!recursive_cone_shift_certificate(
                                base, directions, 15U, signs, 0U, stats
                            )) {
                            std::cout << "SHIFTED_SLOPE_FAIL pattern="
                                      << pattern << " base=";
                            print_counts(base);
                            std::cout << " rays=" << rays[ray][0U] << ','
                                      << rays[ray][1U] << ':'
                                      << rays[ray + 1U][0U] << ','
                                      << rays[ray + 1U][1U]
                                      << " nodes=" << stats.nodes
                                      << " maximum_depth="
                                      << stats.maximum_depth
                                      << " degree=";
                            print_counts(stats.failure_degree);
                            std::cout << '\n';
                            pass = false;
                            break;
                        }
                    }
                    if (!pass) {
                        break;
                    }
                }
                if (pass) {
                    std::cout << "SHIFTED_SLOPE_PASS pattern=" << pattern
                              << " bases=" << bases.size()
                              << " nodes=" << stats.nodes
                              << " leaves=" << stats.leaves
                              << " splits=" << stats.splits
                              << " maximum_depth=" << stats.maximum_depth
                              << '\n';
                }
            }
            return EXIT_SUCCESS;
        }

        if (mode == "adaptive-slope") {
            for (int pattern = 0; pattern < 8; ++pattern) {
                if (requested_pattern >= 0 && pattern != requested_pattern) {
                    continue;
                }
                const Counts signs{
                    1,
                    (pattern & 4) != 0 ? -1 : 1,
                    (pattern & 2) != 0 ? -1 : 1,
                    (pattern & 1) != 0 ? -1 : 1
                };
                FanStats stats;
                bool pass = true;
                for (const Counts& base : bases) {
                    if (!adaptive_slope_certificate(
                            base, std::array<int, 2>{1, 0},
                            std::array<int, 2>{0, 1}, signs, 0U, stats
                        )) {
                        std::cout << "ADAPTIVE_SLOPE_FAIL pattern=" << pattern
                                  << " base=";
                        print_counts(base);
                        std::cout << " nodes=" << stats.nodes
                                  << " maximum_depth="
                                  << stats.maximum_depth
                                  << " rays=" << stats.failure_left[0U]
                                  << ',' << stats.failure_left[1U] << ':'
                                  << stats.failure_right[0U] << ','
                                  << stats.failure_right[1U]
                                  << " degree=";
                        print_counts(stats.failure_degree);
                        std::cout << '\n';
                        pass = false;
                        break;
                    }
                }
                if (pass) {
                    std::cout << "ADAPTIVE_SLOPE_PASS pattern=" << pattern
                              << " bases=" << bases.size()
                              << " nodes=" << stats.nodes
                              << " leaves=" << stats.leaves
                              << " splits=" << stats.splits
                              << " maximum_depth=" << stats.maximum_depth
                              << '\n';
                }
            }
            return EXIT_SUCCESS;
        }

        if (mode == "slope-fan") {
            constexpr int maximum_slope = 16;
            for (int pattern = 0; pattern < 8; ++pattern) {
                if (requested_pattern >= 0 && pattern != requested_pattern) {
                    continue;
                }
                const Counts signs{
                    1,
                    (pattern & 4) != 0 ? -1 : 1,
                    (pattern & 2) != 0 ? -1 : 1,
                    (pattern & 1) != 0 ? -1 : 1
                };
                bool pass = true;
                std::size_t cones = 0U;
                for (const Counts& base : bases) {
                    std::vector<std::array<int, 2>> rays;
                    rays.push_back({1, 0});
                    for (int slope = 1; slope <= maximum_slope; ++slope) {
                        rays.push_back({1, slope});
                    }
                    rays.push_back({0, 1});
                    for (std::size_t ray = 0U; ray + 1U < rays.size();
                         ++ray) {
                        Directions directions{};
                        directions[0U][0U] = rays[ray][0U];
                        directions[0U][1U] = rays[ray][1U];
                        directions[1U][0U] = rays[ray + 1U][0U];
                        directions[1U][1U] = rays[ray + 1U][1U];
                        directions[2U][2U] = 1;
                        directions[3U][3U] = 1;
                        Counts failure{};
                        ++cones;
                        if (!cone_newton_is_nonnegative(
                                base, directions,
                                direction_degree_bounds(directions), signs,
                                failure
                            )) {
                            std::cout << "SLOPE_FAN_FAIL pattern=" << pattern
                                      << " base=";
                            print_counts(base);
                            std::cout << " rays=" << rays[ray][0U] << ','
                                      << rays[ray][1U] << ':'
                                      << rays[ray + 1U][0U] << ','
                                      << rays[ray + 1U][1U]
                                      << " degree=";
                            print_counts(failure);
                            std::cout << " cones=" << cones << '\n';
                            pass = false;
                            break;
                        }
                    }
                    if (!pass) {
                        break;
                    }
                }
                if (pass) {
                    std::cout << "SLOPE_FAN_PASS pattern=" << pattern
                              << " bases=" << bases.size()
                              << " maximum_slope=" << maximum_slope
                              << " cones=" << cones << '\n';
                }
            }
            return EXIT_SUCCESS;
        }

        if (mode == "recursive-fan") {
            Directions identity{};
            for (std::size_t coordinate = 0U; coordinate < 4U; ++coordinate) {
                identity[coordinate][coordinate] = 1;
            }
            for (int pattern = 0; pattern < 8; ++pattern) {
                if (requested_pattern >= 0 && pattern != requested_pattern) {
                    continue;
                }
                const Counts signs{
                    1,
                    (pattern & 4) != 0 ? -1 : 1,
                    (pattern & 2) != 0 ? -1 : 1,
                    (pattern & 1) != 0 ? -1 : 1
                };
                FanStats stats;
                bool pass = true;
                for (const Counts& base : bases) {
                    if (!recursive_fan_certificate(
                            base, identity, signs, 0U, stats
                        )) {
                        std::cout << "RECURSIVE_FAN_FAIL pattern=" << pattern
                                  << " base=";
                        print_counts(base);
                        std::cout << " nodes=" << stats.nodes
                                  << " maximum_depth="
                                  << stats.maximum_depth << '\n';
                        pass = false;
                        break;
                    }
                }
                if (pass) {
                    std::cout << "RECURSIVE_FAN_PASS pattern=" << pattern
                              << " bases=" << bases.size()
                              << " nodes=" << stats.nodes
                              << " leaves=" << stats.leaves
                              << " splits=" << stats.splits
                              << " maximum_depth=" << stats.maximum_depth
                              << '\n';
                }
            }
            return EXIT_SUCCESS;
        }

        if (mode == "fan") {
            for (int pattern = 0; pattern < 8; ++pattern) {
                if (requested_pattern >= 0 && pattern != requested_pattern) {
                    continue;
                }
                const Counts signs{
                    1,
                    (pattern & 4) != 0 ? -1 : 1,
                    (pattern & 2) != 0 ? -1 : 1,
                    (pattern & 1) != 0 ? -1 : 1
                };
                bool pass = true;
                std::array<std::size_t, 6> pair_counts{};
                for (const Counts& base : bases) {
                    bool base_pass = false;
                    std::size_t pair_index = 0U;
                    Counts first_failure{};
                    Counts second_failure{};
                    for (std::size_t first = 0U; first < 4U && !base_pass;
                         ++first) {
                        for (std::size_t second = first + 1U;
                             second < 4U; ++second) {
                            if (two_coordinate_fan_certificate(
                                    base, first, second, signs,
                                    first_failure, second_failure
                                )) {
                                ++pair_counts[pair_index];
                                base_pass = true;
                                break;
                            }
                            if (std::getenv("OUTER_FAN_VERBOSE") != nullptr) {
                                std::cout << "FAN_PAIR_FAIL pattern="
                                          << pattern << " base=";
                                print_counts(base);
                                std::cout << " base_value="
                                          << outer_value(base, signs);
                                std::cout << " pair=" << first << ','
                                          << second << " first_degree=";
                                print_counts(first_failure);
                                std::cout << " second_degree=";
                                print_counts(second_failure);
                                std::cout << '\n';
                            }
                            ++pair_index;
                        }
                    }
                    if (!base_pass) {
                        std::cout << "FAN_NEWTON_FAIL pattern=" << pattern
                                  << " base=";
                        print_counts(base);
                        std::cout << " last_first_degree=";
                        print_counts(first_failure);
                        std::cout << " last_second_degree=";
                        print_counts(second_failure);
                        std::cout << '\n';
                        pass = false;
                        break;
                    }
                }
                if (pass) {
                    std::cout << "FAN_NEWTON_PASS pattern=" << pattern
                              << " bases=" << bases.size()
                              << " pair_counts=";
                    for (std::size_t index = 0U; index < pair_counts.size();
                         ++index) {
                        if (index != 0U) {
                            std::cout << ',';
                        }
                        std::cout << pair_counts[index];
                    }
                    std::cout << '\n';
                }
            }
            return EXIT_SUCCESS;
        }

        if (mode == "recursive") {
            for (int pattern = 0; pattern < 8; ++pattern) {
                if (requested_pattern >= 0 && pattern != requested_pattern) {
                    continue;
                }
                const Counts signs{
                    1,
                    (pattern & 4) != 0 ? -1 : 1,
                    (pattern & 2) != 0 ? -1 : 1,
                    (pattern & 1) != 0 ? -1 : 1
                };
                RecursiveStats stats;
                bool pass = true;
                for (const Counts& base : bases) {
                    std::map<RecursiveKey, bool> memo;
                    if (!recursive_newton_certificate(
                            base, 15U, signs, 0U, stats, memo
                        )) {
                        std::cout << "RECURSIVE_NEWTON_FAIL pattern="
                                  << pattern << " base=";
                        print_counts(base);
                        std::cout << " nodes=" << stats.nodes
                                  << " maximum_depth="
                                  << stats.maximum_depth
                                  << " maximum_base=";
                        print_counts(stats.maximum_base);
                        std::cout << '\n';
                        pass = false;
                        break;
                    }
                }
                if (pass) {
                    std::cout << "RECURSIVE_NEWTON_PASS pattern=" << pattern
                              << " bases=" << bases.size()
                              << " nodes=" << stats.nodes
                              << " leaves=" << stats.positive_leaves
                              << " splits=" << stats.split_nodes
                              << " maximum_depth=" << stats.maximum_depth
                              << " maximum_base=";
                    print_counts(stats.maximum_base);
                    std::cout << '\n';
                }
            }
            return EXIT_SUCCESS;
        }

        constexpr Counts degrees{8, 4, 4, 4};
        for (int pattern = 0; pattern < 8; ++pattern) {
            const Counts signs{
                1,
                (pattern & 4) != 0 ? -1 : 1,
                (pattern & 2) != 0 ? -1 : 1,
                (pattern & 1) != 0 ? -1 : 1
            };
            bool passes = true;
            for (const Counts& base : bases) {
                for (int first = 0; first <= degrees[0] && passes; ++first) {
                    for (int second = 0;
                         second <= degrees[1] && passes;
                         ++second) {
                        for (int third = 0;
                             third <= degrees[2] && passes;
                             ++third) {
                            for (int fourth = 0;
                                 fourth <= degrees[3];
                                 ++fourth) {
                                const Counts degree{
                                    first, second, third, fourth
                                };
                                const cpp_int coefficient = finite_difference(
                                    base, degree, signs
                                );
                                if (coefficient < 0) {
                                    std::cout << "NEWTON_FAIL pattern="
                                              << pattern << " base=";
                                    print_counts(base);
                                    std::cout << " degree=";
                                    print_counts(degree);
                                    std::cout << " value=" << coefficient
                                              << '\n';
                                    passes = false;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (!passes) {
                    break;
                }
            }
            if (passes) {
                std::cout << "NEWTON_PASS pattern=" << pattern
                          << " bases=" << bases.size() << '\n';
            }
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
