#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace exact {

using Integer = boost::multiprecision::cpp_int;
using Poly = std::vector<Integer>;

void trim(Poly& polynomial) {
    while (!polynomial.empty() && polynomial.back() == 0) {
        polynomial.pop_back();
    }
}

Poly multiply(const Poly& first, const Poly& second) {
    if (first.empty() || second.empty()) {
        return {};
    }
    Poly product(first.size() + second.size() - 1);
    for (std::size_t first_index = 0;
         first_index < first.size(); ++first_index) {
        for (std::size_t second_index = 0;
             second_index < second.size(); ++second_index) {
            product[first_index + second_index] +=
                first[first_index] * second[second_index];
        }
    }
    trim(product);
    return product;
}

Poly divide_exact(Poly numerator, const Poly& denominator) {
    if (denominator.empty()) {
        throw std::runtime_error("division by zero polynomial");
    }
    Poly quotient(
        numerator.size() >= denominator.size()
            ? numerator.size() - denominator.size() + 1
            : 0
    );
    while (
        !numerator.empty()
        && numerator.size() >= denominator.size()
    ) {
        if (numerator.back() % denominator.back() != 0) {
            throw std::runtime_error("nonexact polynomial division");
        }
        const Integer coefficient =
            numerator.back() / denominator.back();
        const std::size_t shift =
            numerator.size() - denominator.size();
        quotient[shift] = coefficient;
        for (std::size_t index = 0;
             index < denominator.size(); ++index) {
            numerator[shift + index] -=
                coefficient * denominator[index];
        }
        trim(numerator);
    }
    if (!numerator.empty()) {
        throw std::runtime_error(
            "polynomial division left a remainder"
        );
    }
    trim(quotient);
    return quotient;
}

Poly cyclotomic(int order) {
    Poly polynomial(static_cast<std::size_t>(order + 1));
    polynomial[0] = -1;
    polynomial[static_cast<std::size_t>(order)] = 1;
    for (int divisor = 1; divisor < order; ++divisor) {
        if (order % divisor == 0) {
            polynomial = divide_exact(
                polynomial,
                cyclotomic(divisor)
            );
        }
    }
    return polynomial;
}

struct Field {
    int degree = 0;
    Poly minimal_polynomial;
    using Element = Poly;

    void initialize(int half_period) {
        const Poly phi = cyclotomic(2 * half_period);
        degree = static_cast<int>(phi.size()) - 1;
        if ((degree & 1) != 0) {
            throw std::runtime_error(
                "real cyclotomic degree is not even"
            );
        }
        degree /= 2;
        Poly residual = phi;
        minimal_polynomial.assign(
            static_cast<std::size_t>(degree + 1),
            Integer{0}
        );
        const Poly z_squared_plus_one = {1, 0, 1};
        for (int index = degree; index >= 0; --index) {
            const std::size_t top =
                static_cast<std::size_t>(degree + index);
            const Integer coefficient =
                top < residual.size() ? residual[top] : Integer{0};
            minimal_polynomial[
                static_cast<std::size_t>(index)
            ] = coefficient;
            if (coefficient == 0) {
                continue;
            }
            Poly term = {1};
            for (int exponent = 0;
                 exponent < index; ++exponent) {
                term = multiply(term, z_squared_plus_one);
            }
            Poly shifted(
                static_cast<std::size_t>(degree - index)
            );
            shifted.insert(
                shifted.end(),
                term.begin(),
                term.end()
            );
            if (shifted.size() > residual.size()) {
                residual.resize(shifted.size());
            }
            for (std::size_t position = 0;
                 position < shifted.size(); ++position) {
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
                "minimal-polynomial extraction failed"
            );
        }
    }

    Element from_integer(long long value) const {
        return value == 0 ? Element{} : Element{Integer{value}};
    }

    Element reduce(Poly polynomial) const {
        while (static_cast<int>(polynomial.size()) > degree) {
            const Integer coefficient = polynomial.back();
            const std::size_t shift =
                polynomial.size() - minimal_polynomial.size();
            for (std::size_t index = 0;
                 index < minimal_polynomial.size(); ++index) {
                polynomial[shift + index] -=
                    coefficient * minimal_polynomial[index];
            }
            trim(polynomial);
        }
        return polynomial;
    }

    Element add(const Element& first, const Element& second) const {
        Element sum = first;
        if (second.size() > sum.size()) {
            sum.resize(second.size());
        }
        for (std::size_t index = 0;
             index < second.size(); ++index) {
            sum[index] += second[index];
        }
        trim(sum);
        return sum;
    }

    Element subtract(
        const Element& first,
        const Element& second
    ) const {
        Element difference = first;
        if (second.size() > difference.size()) {
            difference.resize(second.size());
        }
        for (std::size_t index = 0;
             index < second.size(); ++index) {
            difference[index] -= second[index];
        }
        trim(difference);
        return difference;
    }

    Element multiply_elements(
        const Element& first,
        const Element& second
    ) const {
        return reduce(multiply(first, second));
    }

    Element power(Element base, unsigned long exponent) const {
        Element result = from_integer(1);
        while (exponent != 0UL) {
            if ((exponent & 1UL) != 0UL) {
                result = multiply_elements(result, base);
            }
            base = multiply_elements(base, base);
            exponent >>= 1U;
        }
        return result;
    }

    Element dickson(int index) const {
        const Element generator = {0, 1};
        Element previous = from_integer(2);
        Element current = generator;
        if (index == 0) {
            return previous;
        }
        for (int step = 1; step < index; ++step) {
            const Element next = subtract(
                multiply_elements(current, generator),
                previous
            );
            previous = current;
            current = next;
        }
        return current;
    }

    Element chebyshev(int index, const Element& value) const {
        Element previous = from_integer(1);
        Element current = value;
        if (index == 0) {
            return previous;
        }
        for (int step = 1; step < index; ++step) {
            const Element next = subtract(
                multiply_elements(value, current),
                previous
            );
            previous = current;
            current = next;
        }
        return current;
    }
};

}  // namespace exact

namespace {

int folded_mode(int mode_number, int multiplier, int period) {
    int value = (mode_number * multiplier) % period;
    if (value < 0) {
        value += period;
    }
    const int half_period = period / 2;
    if (value > half_period) {
        value = period - value;
    }
    if (value <= 0 || value >= half_period) {
        throw std::runtime_error("Galois fold left the mode set");
    }
    return value - 1;
}

std::set<std::pair<int, int>> orbit_of(
    int half_period,
    int first,
    int second
) {
    const int period = 2 * half_period;
    std::set<std::pair<int, int>> orbit;
    for (int multiplier = 1;
         multiplier < period; ++multiplier) {
        if (std::gcd(multiplier, period) != 1) {
            continue;
        }
        int image_first = folded_mode(
            first + 1,
            multiplier,
            period
        );
        int image_second = folded_mode(
            second + 1,
            multiplier,
            period
        );
        if (image_first == image_second) {
            throw std::runtime_error(
                "Galois orbit collapsed a mode pair"
            );
        }
        if (image_first > image_second) {
            std::swap(image_first, image_second);
        }
        orbit.emplace(image_first, image_second);
    }
    return orbit;
}

exact::Field::Element scaled_orbit_contribution(
    const exact::Field& field,
    int label,
    int exponent,
    const std::set<std::pair<int, int>>& orbit
) {
    exact::Field::Element total;
    for (const auto& pair : orbit) {
        const int first_number = pair.first + 1;
        const int second_number = pair.second + 1;
        const auto first_node = field.dickson(first_number);
        const auto second_node = field.dickson(second_number);
        const auto first_lambda =
            field.chebyshev(label, first_node);
        const auto second_lambda =
            field.chebyshev(label, second_node);
        const auto first_weight = field.subtract(
            field.from_integer(2),
            field.dickson(2 * first_number)
        );
        const auto second_weight = field.subtract(
            field.from_integer(2),
            field.dickson(2 * second_number)
        );
        const int first_sign =
            (pair.first & 1) == 0 ? 1 : -1;
        const int second_sign =
            (pair.second & 1) == 0 ? 1 : -1;
        auto term = field.multiply_elements(
            first_weight,
            second_weight
        );
        term = field.multiply_elements(
            term,
            field.from_integer(first_sign - second_sign)
        );
        term = field.multiply_elements(
            term,
            field.subtract(first_lambda, second_lambda)
        );
        term = field.multiply_elements(
            term,
            field.power(
                field.add(first_lambda, second_lambda),
                static_cast<unsigned long>(exponent)
            )
        );
        total = field.add(total, term);
    }
    return total;
}

}  // namespace

int main() {
    try {
        exact::Field field;
        field.initialize(14);
        const auto orbit = orbit_of(14, 1, 6);
        const std::set<std::pair<int, int>> expected_orbit = {
            {1, 6},
            {5, 6},
            {6, 9}
        };
        if (orbit != expected_orbit) {
            throw std::runtime_error("unexpected Galois orbit");
        }
        const auto scaled =
            scaled_orbit_contribution(field, 4, 9, orbit);
        if (scaled != field.from_integer(-50960)) {
            throw std::runtime_error(
                "exact scaled contribution mismatch"
            );
        }
        std::cout
            << "SU2_TERMINAL_GALOIS_WITNESS"
            << " level=12 label=4 pair_power=4 exponent=9"
            << " seed_pair=(1,6) orbit_size=3"
            << " normalization_denominator=784"
            << " scaled_contribution=-50960"
            << " contribution=-65"
            << " result=PASS_EXACT_NEGATIVE_CONTROL\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_GALOIS_WITNESS FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
