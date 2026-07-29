#include <boost/multiprecision/cpp_int.hpp>

#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Polynomial = std::vector<Integer>;

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value <= 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error(
            std::string(name) + " must be a positive integer"
        );
    }
    return static_cast<int>(value);
}

Polynomial raise_power(
    const Polynomial& coefficients,
    int exponent,
    int lobes
) {
    const int degree = static_cast<int>(coefficients.size()) - 1;
    Polynomial result(coefficients.size() + 1, Integer(0));
    const Integer n_squared = Integer(lobes) * Integer(lobes);
    const Integer exponent_squared =
        Integer(exponent) * Integer(exponent);

    for (int j = 0; j <= degree + 1; ++j) {
        const std::size_t index = static_cast<std::size_t>(j);
        Integer l_value = 0;
        if (j + 1 <= degree) {
            l_value +=
                Integer(2 * (j + 1) * (2 * j + 1))
                * coefficients[static_cast<std::size_t>(j + 1)];
        }
        if (j <= degree) {
            l_value +=
                Integer(8 * j * j + 4 * j + 1)
                * coefficients[index];
            result[index] += exponent_squared * coefficients[index];
        }
        if (j >= 1 && j - 1 <= degree) {
            l_value +=
                Integer(2 * j * (2 * j - 1))
                * coefficients[static_cast<std::size_t>(j - 1)];
        }
        result[index] += n_squared * l_value;
    }
    return result;
}

Integer coefficient_at(const Polynomial& coefficients, int index) {
    if (
        index < 0
        || index >= static_cast<int>(coefficients.size())
    ) {
        return 0;
    }
    return coefficients[static_cast<std::size_t>(index)];
}

Integer shape_defect(
    const Polynomial& coefficients,
    int degree,
    int index
) {
    if (index < 0 || index > degree) {
        return 0;
    }
    return
        Integer(index + 1) * coefficient_at(coefficients, index + 1)
        - Integer(degree - index)
            * coefficient_at(coefficients, index);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: prove_su2_odd_twisted_endpoint "
                "maximum_odd_lobes maximum_power"
            );
        }
        const int maximum_lobes =
            parse_positive(argv[1], "maximum_odd_lobes");
        const int maximum_power =
            parse_positive(argv[2], "maximum_power");

        std::size_t polynomial_rows = 0;
        std::size_t shape_rows = 0;
        std::size_t induction_rows = 0;
        std::size_t endpoint_rows = 0;
        std::size_t failures = 0;

        for (int lobes = 3; lobes <= maximum_lobes; lobes += 2) {
            Polynomial coefficients{Integer(1)};
            for (int power = 1; power <= maximum_power; ++power) {
                const int degree = power - 1;
                const int exponent = 2 * power - 1;
                ++polynomial_rows;

                for (const Integer& coefficient : coefficients) {
                    if (coefficient <= 0) {
                        ++failures;
                    }
                }

                for (int j = 0; j < degree; ++j) {
                    const std::size_t index =
                        static_cast<std::size_t>(j);
                    const Integer shape =
                        Integer(j + 1)
                            * coefficients[
                                static_cast<std::size_t>(j + 1)
                            ]
                        - Integer(degree - j) * coefficients[index];
                    ++shape_rows;
                    if (shape < 0) {
                        ++failures;
                        if (failures <= 8) {
                            std::cout
                                << "SHAPE_FAILURE"
                                << " lobes=" << lobes
                                << " power=" << power
                                << " coefficient=" << j
                                << " value=" << shape << '\n';
                        }
                    }
                }

                for (int j = 0; j <= degree; ++j) {
                    const std::size_t index =
                        static_cast<std::size_t>(j);
                    Integer endpoint =
                        Integer(exponent + 4 * lobes * j)
                        * coefficients[index];
                    if (j >= 1) {
                        endpoint +=
                            Integer(
                                2 * lobes * (2 * j - 1) - exponent
                            ) * coefficients[
                                static_cast<std::size_t>(j - 1)
                            ];
                    }
                    ++endpoint_rows;
                    if (endpoint <= 0) {
                        ++failures;
                        if (failures <= 8) {
                            std::cout
                                << "ENDPOINT_FAILURE"
                                << " lobes=" << lobes
                                << " power=" << power
                                << " coefficient=" << j
                                << " value=" << endpoint << '\n';
                        }
                    }
                }

                const Polynomial raised =
                    raise_power(coefficients, exponent, lobes);
                const Integer n_squared =
                    Integer(lobes) * Integer(lobes);
                const Integer exponent_squared =
                    Integer(exponent) * Integer(exponent);
                for (int j = 0; j <= degree; ++j) {
                    const Integer left =
                        Integer(j + 1)
                            * coefficient_at(raised, j + 1)
                        - Integer(degree + 1 - j)
                            * coefficient_at(raised, j);
                    const Integer right =
                        Integer(2 * (j + 1) * (2 * j + 3))
                            * n_squared
                            * shape_defect(
                                coefficients,
                                degree,
                                j + 1
                            )
                        + (
                            n_squared
                                * Integer(
                                    8 * j * j
                                    + 8 * j
                                    + 4 * degree
                                    + 5
                                )
                            + exponent_squared
                        ) * shape_defect(
                            coefficients,
                            degree,
                            j
                        )
                        + Integer(2 * j * (2 * j - 1))
                            * n_squared
                            * shape_defect(
                                coefficients,
                                degree,
                                j - 1
                            )
                        + exponent_squared * (n_squared - 1)
                            * coefficient_at(coefficients, j);
                    ++induction_rows;
                    if (left != right) {
                        ++failures;
                        if (failures <= 8) {
                            std::cout
                                << "INDUCTION_IDENTITY_FAILURE"
                                << " lobes=" << lobes
                                << " power=" << power
                                << " coefficient=" << j
                                << " left=" << left
                                << " right=" << right << '\n';
                        }
                    }
                }
                coefficients = raised;
            }
        }

        std::cout
            << "SU2_ODD_TWISTED_ENDPOINT_EXACT"
            << " maximum_odd_lobes=" << maximum_lobes
            << " maximum_power=" << maximum_power
            << " polynomial_rows=" << polynomial_rows
            << " shape_rows=" << shape_rows
            << " induction_rows=" << induction_rows
            << " endpoint_rows=" << endpoint_rows
            << " failures=" << failures
            << " result="
            << (
                failures == 0
                    ? "PASS_EXACT_INTEGER_DIAGNOSTIC"
                    : "COUNTEREXAMPLE"
            )
            << '\n';
        return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
