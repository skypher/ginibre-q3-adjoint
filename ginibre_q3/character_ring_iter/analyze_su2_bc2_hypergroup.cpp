#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

using boost::multiprecision::cpp_int;

namespace {

using Rational = boost::rational<cpp_int>;
using Exponent = std::pair<int, int>;
using Polynomial = std::map<Exponent, Rational>;

struct Weight {
    int first = 0;
    int second = 0;
};

void add_term(Polynomial& polynomial, int x_degree, int y_degree,
              const Rational& coefficient) {
    if (coefficient == 0) {
        return;
    }
    const Exponent exponent{x_degree, y_degree};
    polynomial[exponent] += coefficient;
    if (polynomial[exponent] == 0) {
        polynomial.erase(exponent);
    }
}

Polynomial add(const Polynomial& left, const Polynomial& right,
               const Rational& right_scale = Rational{1}) {
    Polynomial result = left;
    for (const auto& [exponent, coefficient] : right) {
        add_term(result, exponent.first, exponent.second,
                 right_scale * coefficient);
    }
    return result;
}

Polynomial multiply(const Polynomial& left, const Polynomial& right) {
    Polynomial result;
    for (const auto& [left_exponent, left_coefficient] : left) {
        for (const auto& [right_exponent, right_coefficient] : right) {
            add_term(result, left_exponent.first + right_exponent.first,
                     left_exponent.second + right_exponent.second,
                     left_coefficient * right_coefficient);
        }
    }
    return result;
}

Polynomial monomial(int x_degree, int y_degree,
                    const Rational& coefficient = Rational{1}) {
    Polynomial result;
    add_term(result, x_degree, y_degree, coefficient);
    return result;
}

Polynomial character(int label, bool x_variable) {
    if (label < 0) {
        return {};
    }
    Polynomial previous = monomial(0, 0);
    if (label == 0) {
        return previous;
    }
    Polynomial current = x_variable ? monomial(1, 0) : monomial(0, 1);
    for (int degree = 1; degree < label; ++degree) {
        const Polynomial variable = x_variable ? monomial(1, 0)
                                               : monomial(0, 1);
        Polynomial next = add(multiply(variable, current), previous,
                              Rational{-1});
        previous = std::move(current);
        current = std::move(next);
    }
    return current;
}

int maximum_x_degree(const Polynomial& polynomial) {
    int maximum = -1;
    for (const auto& [exponent, coefficient] : polynomial) {
        if (coefficient != 0) {
            maximum = std::max(maximum, exponent.first);
        }
    }
    return maximum;
}

Polynomial divide_by_x_minus_y(const Polynomial& numerator) {
    Polynomial remainder = numerator;
    Polynomial quotient;
    while (!remainder.empty()) {
        const int degree = maximum_x_degree(remainder);
        if (degree <= 0) {
            throw std::runtime_error("nonzero remainder in x-y division");
        }
        std::vector<std::pair<int, Rational>> leading;
        for (const auto& [exponent, coefficient] : remainder) {
            if (exponent.first == degree) {
                leading.emplace_back(exponent.second, coefficient);
            }
        }
        for (const auto& [y_degree, coefficient] : leading) {
            add_term(quotient, degree - 1, y_degree, coefficient);
            add_term(remainder, degree, y_degree, -coefficient);
            add_term(remainder, degree - 1, y_degree + 1, coefficient);
        }
    }
    return quotient;
}

Polynomial sp4_character(int first, int second) {
    if (first < second || second < 0) {
        return {};
    }
    const Polynomial forward = multiply(character(first + 1, true),
                                        character(second, false));
    const Polynomial reverse = multiply(character(second, true),
                                        character(first + 1, false));
    return divide_by_x_minus_y(add(forward, reverse, Rational{-1}));
}

cpp_int catalan(int index) {
    if (index < 0) {
        return 0;
    }
    cpp_int binomial = 1;
    for (int j = 1; j <= index; ++j) {
        binomial *= index + j;
        binomial /= j;
    }
    return binomial / (index + 1);
}

Rational semicircle_moment(int degree) {
    if (degree < 0 || (degree & 1) != 0) {
        return 0;
    }
    return Rational{catalan(degree / 2)};
}

Polynomial delta_power(int exponent) {
    Polynomial result = monomial(0, 0);
    const Polynomial delta = add(monomial(1, 0), monomial(0, 1),
                                 Rational{-1});
    for (int index = 0; index < exponent; ++index) {
        result = multiply(result, delta);
    }
    return result;
}

Rational integral(const Polynomial& polynomial) {
    Rational result{0};
    for (const auto& [exponent, coefficient] : polynomial) {
        result += coefficient * semicircle_moment(exponent.first)
                  * semicircle_moment(exponent.second);
    }
    return result;
}

Rational inner_product(const Polynomial& left, const Polynomial& right,
                       const Polynomial& density) {
    return integral(multiply(multiply(left, right), density));
}

std::string rational_string(const Rational& value) {
    if (value.denominator() == 1) {
        return value.numerator().convert_to<std::string>();
    }
    return value.numerator().convert_to<std::string>() + "/"
           + value.denominator().convert_to<std::string>();
}

std::string weight_string(const Weight& weight) {
    return "(" + std::to_string(weight.first) + ","
           + std::to_string(weight.second) + ")";
}

std::vector<Weight> weights_through(int maximum_degree) {
    std::vector<Weight> result;
    for (int total = 0; total <= maximum_degree; ++total) {
        for (int second = total / 2; second >= 0; --second) {
            const int first = total - second;
            if (first >= second) {
                result.push_back(Weight{first, second});
            }
        }
    }
    return result;
}

struct OrthogonalBasis {
    std::vector<Weight> weights;
    std::vector<Polynomial> polynomials;
    std::vector<Rational> norms;
};

OrthogonalBasis make_basis(int maximum_degree, const Polynomial& density) {
    OrthogonalBasis basis;
    basis.weights = weights_through(maximum_degree);
    for (const Weight& weight : basis.weights) {
        Polynomial next = sp4_character(weight.first, weight.second);
        for (std::size_t index = 0; index < basis.polynomials.size(); ++index) {
            const Rational projection = inner_product(
                next, basis.polynomials[index], density
            ) / basis.norms[index];
            next = add(next, basis.polynomials[index], -projection);
        }
        const Rational norm = inner_product(next, next, density);
        if (norm <= 0) {
            throw std::runtime_error("nonpositive Gram-Schmidt norm");
        }
        basis.polynomials.push_back(std::move(next));
        basis.norms.push_back(norm);
    }
    return basis;
}

std::vector<Rational> expand(const Polynomial& polynomial,
                             const OrthogonalBasis& basis,
                             const Polynomial& density) {
    std::vector<Rational> result;
    result.reserve(basis.polynomials.size());
    for (std::size_t index = 0; index < basis.polynomials.size(); ++index) {
        result.push_back(inner_product(polynomial, basis.polynomials[index],
                                       density) / basis.norms[index]);
    }
    return result;
}

bool report_first_negative(const std::string& name,
                           const std::vector<Rational>& coefficients,
                           const OrthogonalBasis& basis) {
    for (std::size_t index = 0; index < coefficients.size(); ++index) {
        if (coefficients[index] < 0) {
            std::cout << name << " negative_weight="
                      << weight_string(basis.weights[index])
                      << " coefficient="
                      << rational_string(coefficients[index]) << '\n';
            return true;
        }
    }
    return false;
}

void row_products(const OrthogonalBasis& basis, const Polynomial& density,
                  int half_minus_count, int maximum_label, int depth,
                  int first_label, Polynomial product,
                  std::vector<int>& labels, bool& failure) {
    if (depth == half_minus_count) {
        const auto coefficients = expand(product, basis, density);
        std::string name = "row_product labels=[";
        for (std::size_t index = 0; index < labels.size(); ++index) {
            if (index != 0U) {
                name += ',';
            }
            name += std::to_string(labels[index]);
        }
        name += ']';
        failure = report_first_negative(name, coefficients, basis) || failure;
        return;
    }
    for (int label = first_label; label <= maximum_label; ++label) {
        labels.push_back(label);
        row_products(basis, density, half_minus_count, maximum_label,
                     depth + 1, label,
                     multiply(product, sp4_character(label - 1, 0)),
                     labels, failure);
        labels.pop_back();
        if (failure) {
            return;
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_bc2_hypergroup HALF_MINUS_COUNT "
                "MAXIMUM_LABEL"
            );
        }
        const int half_minus_count = std::stoi(argv[1]);
        const int maximum_label = std::stoi(argv[2]);
        if (half_minus_count < 1 || maximum_label < 1) {
            throw std::runtime_error("arguments must be positive");
        }
        const int maximum_degree = half_minus_count * (maximum_label - 1)
                                   + maximum_label;
        const Polynomial density = delta_power(2 * half_minus_count);
        const OrthogonalBasis basis = make_basis(maximum_degree, density);

        bool failure = false;
        std::vector<int> labels;
        row_products(basis, density, half_minus_count, maximum_label, 0, 1,
                     monomial(0, 0), labels, failure);

        if (!failure) {
            for (int row_label = 1; row_label <= maximum_label; ++row_label) {
                Polynomial row = sp4_character(row_label - 1, 0);
                for (int plus_label = 1; plus_label <= maximum_label;
                     ++plus_label) {
                    const Polynomial plus = add(character(plus_label, true),
                                                character(plus_label, false));
                    const auto coefficients = expand(
                        multiply(row, plus), basis, density
                    );
                    const std::string name
                        = "plus_action row=" + std::to_string(row_label)
                          + " plus=" + std::to_string(plus_label);
                    if (report_first_negative(name, coefficients, basis)) {
                        failure = true;
                        break;
                    }
                }
                if (failure) {
                    break;
                }
            }
        }

        std::cout << "SU2_BC2_HYPERGROUP half_minus_count="
                  << half_minus_count << " maximum_label=" << maximum_label
                  << " result=" << (failure ? "FAIL" : "PASS") << '\n';
        return failure ? EXIT_FAILURE : EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
