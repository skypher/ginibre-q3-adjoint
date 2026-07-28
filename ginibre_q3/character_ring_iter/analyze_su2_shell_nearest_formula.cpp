#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using Matrix = std::vector<std::vector<int>>;
using Exponent = std::array<int, 2>;

std::vector<bool>* activation_recorder = nullptr;

int parse_positive(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error("bound must be a positive integer");
    }
    return static_cast<int>(value);
}

Integer binomial_nonnegative(int top, int order) {
    if (activation_recorder != nullptr) {
        activation_recorder->push_back(top >= order);
    }
    if (top < order || order < 0) {
        return 0;
    }
    Integer result = 1;
    for (int index = 1; index <= order; ++index) {
        result *= top - order + index;
        result /= index;
    }
    return result;
}

Integer binomial_small(int top, int order) {
    if (order < 0 || top < order) {
        return 0;
    }
    Integer result = 1;
    for (int index = 1; index <= order; ++index) {
        result *= top - order + index;
        result /= index;
    }
    return result;
}

Integer ordinary_root_multiplicity(
    int power,
    int label,
    int target
) {
    if (target < 0) {
        return 0;
    }
    if (power == 0) {
        return target == 0 ? Integer(1) : Integer(0);
    }
    if (power == 1) {
        return target == label ? Integer(1) : Integer(0);
    }
    const int degree = power - 2;
    const int coefficient = power * label - target;
    Integer result = 0;
    for (int selected = 0; selected <= power; ++selected) {
        const int top =
            coefficient - selected * (2 * label + 1) + degree;
        Integer term =
            binomial_small(power, selected)
            * binomial_nonnegative(top, degree);
        if ((selected % 2) == 0) {
            result += term;
        } else {
            result -= term;
        }
    }
    return result;
}

Integer finite_root_multiplicity_nearest(
    int power,
    int label,
    int target
) {
    const int half_level = 2 * label + 1;
    return ordinary_root_multiplicity(power, label, target)
        - ordinary_root_multiplicity(
            power,
            label,
            2 * half_level + 1 - target
        )
        + ordinary_root_multiplicity(
            power,
            label,
            2 * half_level + 2 + target
        );
}

Integer plus_root_multiplicity_nearest(
    int power,
    int label,
    int target
) {
    const int half_level = 2 * label + 1;
    return finite_root_multiplicity_nearest(power, label, target)
        + finite_root_multiplicity_nearest(
            power,
            label,
            half_level - target
        );
}

Integer crossing_column_prefix(
    int power,
    int label,
    int target
) {
    return plus_root_multiplicity_nearest(
        power + 1,
        label,
        target
    ) - plus_root_multiplicity_nearest(
        power,
        label,
        label - target
    );
}

std::pair<Integer, Integer> endpoint_margins(int label, int target) {
    std::vector<Integer> prefix(6U);
    for (int power = 0; power <= 5; ++power) {
        prefix[static_cast<std::size_t>(power)] =
            crossing_column_prefix(power, label, target);
    }
    const Integer f4 = 2 * label + 1;
    const Integer f5 = 2 * Integer(label) * (label + 2);
    return {
        f4 * (prefix[1] + prefix[3] + prefix[5])
            - f5 * (prefix[2] + prefix[4]),
        f4 * (prefix[2] + prefix[4])
            - f5 * (prefix[1] + prefix[3])
    };
}

Integer integer_power(int base, int exponent) {
    Integer result = 1;
    for (int index = 0; index < exponent; ++index) {
        result *= base;
    }
    return result;
}

std::vector<Exponent> degree_five_monomials() {
    std::vector<Exponent> result;
    for (int total = 0; total <= 5; ++total) {
        for (int q_degree = 0; q_degree <= total; ++q_degree) {
            result.push_back(
                Exponent{q_degree, total - q_degree}
            );
        }
    }
    return result;
}

std::vector<Rational> interpolate_margin(bool reflected) {
    const std::vector<Exponent> monomials = degree_five_monomials();
    std::vector<std::vector<Rational>> equations;
    for (int label = 1; label <= 12; ++label) {
        for (int target = 1; target <= label; ++target) {
            std::vector<Rational> row;
            row.reserve(monomials.size() + 1U);
            for (const Exponent& exponent : monomials) {
                row.emplace_back(
                    integer_power(label, exponent[0])
                    * integer_power(target, exponent[1])
                );
            }
            const auto margins = endpoint_margins(label, target);
            row.emplace_back(reflected ? margins.second : margins.first);
            equations.push_back(std::move(row));
        }
    }

    std::size_t pivot_row = 0U;
    for (std::size_t column = 0U;
         column < monomials.size();
         ++column) {
        std::size_t selected = pivot_row;
        while (
            selected < equations.size()
            && equations[selected][column] == 0
        ) {
            ++selected;
        }
        if (selected == equations.size()) {
            throw std::runtime_error("polynomial interpolation lost rank");
        }
        std::swap(equations[pivot_row], equations[selected]);
        const Rational pivot = equations[pivot_row][column];
        for (Rational& value : equations[pivot_row]) {
            value /= pivot;
        }
        for (std::size_t row = 0U; row < equations.size(); ++row) {
            if (row == pivot_row || equations[row][column] == 0) {
                continue;
            }
            const Rational factor = equations[row][column];
            for (std::size_t entry = column;
                 entry <= monomials.size();
                 ++entry) {
                equations[row][entry] -=
                    factor * equations[pivot_row][entry];
            }
        }
        ++pivot_row;
    }
    for (std::size_t row = pivot_row; row < equations.size(); ++row) {
        bool zero_left = true;
        for (std::size_t column = 0U;
             column < monomials.size();
             ++column) {
            zero_left = zero_left && equations[row][column] == 0;
        }
        if (zero_left && equations[row].back() != 0) {
            throw std::runtime_error(
                "nearest margin is not one degree-five polynomial"
            );
        }
    }
    std::vector<Rational> coefficients(monomials.size());
    for (std::size_t row = 0U; row < monomials.size(); ++row) {
        coefficients[row] = equations[row].back();
    }
    return coefficients;
}

Integer factorial(int value) {
    Integer result = 1;
    for (int index = 2; index <= value; ++index) {
        result *= index;
    }
    return result;
}

std::map<Exponent, Rational> pull_back_to_slacks(
    const std::vector<Rational>& coefficients
) {
    const std::vector<Exponent> monomials = degree_five_monomials();
    std::map<Exponent, Rational> result;
    for (std::size_t index = 0U; index < monomials.size(); ++index) {
        const int q_degree = monomials[index][0];
        const int target_degree = monomials[index][1];
        for (int x_from_q = 0; x_from_q <= q_degree; ++x_from_q) {
            for (int y_from_q = 0;
                 y_from_q <= q_degree - x_from_q;
                 ++y_from_q) {
                const int constant_from_q =
                    q_degree - x_from_q - y_from_q;
                const Rational q_coefficient(
                    factorial(q_degree),
                    factorial(x_from_q)
                        * factorial(y_from_q)
                        * factorial(constant_from_q)
                );
                for (int x_from_target = 0;
                     x_from_target <= target_degree;
                     ++x_from_target) {
                    const Rational target_coefficient =
                        binomial_small(
                            target_degree,
                            x_from_target
                        );
                    result[Exponent{
                        x_from_q + x_from_target,
                        y_from_q
                    }] += coefficients[index]
                        * q_coefficient * target_coefficient;
                }
            }
        }
    }
    return result;
}

Rational evaluate_polynomial(
    const std::vector<Rational>& coefficients,
    int label,
    int target
) {
    const std::vector<Exponent> monomials = degree_five_monomials();
    Rational result = 0;
    for (std::size_t index = 0U; index < monomials.size(); ++index) {
        result += coefficients[index]
            * integer_power(label, monomials[index][0])
            * integer_power(target, monomials[index][1]);
    }
    return result;
}

void print_certificate(
    const std::string& name,
    const std::map<Exponent, Rational>& coefficients
) {
    std::size_t nonzero = 0U;
    for (const auto& [exponent, coefficient] : coefficients) {
        if (coefficient == 0) {
            continue;
        }
        ++nonzero;
        if (coefficient < 0) {
            throw std::runtime_error(
                name + " has a negative slack coefficient"
            );
        }
        std::cout
            << "NEAREST_CERTIFICATE_TERM"
            << " margin=" << name
            << " x_degree=" << exponent[0]
            << " y_degree=" << exponent[1]
            << " coefficient=" << coefficient << '\n';
    }
    std::cout
        << "NEAREST_CERTIFICATE"
        << " margin=" << name
        << " nonzero_terms=" << nonzero
        << " result=PASS_NONNEGATIVE_SLACK_EXPANSION\n";
}

bool fuses_half(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= source + label
        && source + target + label <= 2 * level;
}

std::vector<Integer> multiply_row(
    const std::vector<Integer>& state,
    const Matrix& matrix
) {
    std::vector<Integer> next(matrix.size());
    for (std::size_t source = 0; source < matrix.size(); ++source) {
        for (std::size_t target = 0; target < matrix.size(); ++target) {
            next[target] += state[source] * matrix[source][target];
        }
    }
    return next;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2 && argc != 3) {
            throw std::runtime_error(
                "usage: MAXIMUM_LABEL [--table|--certificate]"
            );
        }
        const int maximum_label = parse_positive(argv[1]);
        const bool print_table =
            argc == 3 && std::string(argv[2]) == "--table";
        const bool print_slack_certificate =
            argc == 3 && std::string(argv[2]) == "--certificate";
        if (argc == 3 && !print_table && !print_slack_certificate) {
            throw std::runtime_error(
                "optional flag must be --table or --certificate"
            );
        }
        std::uint64_t formula_rows = 0U;
        std::uint64_t recurrence_checks = 0U;
        std::uint64_t negative_same_endpoint = 0U;
        std::uint64_t negative_reflected_endpoint = 0U;
        std::uint64_t polynomial_identity_checks = 0U;
        std::set<std::vector<bool>> activation_masks;
        Integer minimum_same = 0;
        Integer minimum_reflected = 0;
        bool initialized_minimum = false;
        std::vector<Rational> same_polynomial;
        std::vector<Rational> reflected_polynomial;
        if (print_slack_certificate) {
            same_polynomial = interpolate_margin(false);
            reflected_polynomial = interpolate_margin(true);
        }

        for (int label = 1; label <= maximum_label; ++label) {
            const int level = 2 * label + 1;
            const int size = label + 1;
            Matrix plus(
                static_cast<std::size_t>(size),
                std::vector<int>(static_cast<std::size_t>(size))
            );
            Matrix minus = plus;
            for (int source = 0; source < size; ++source) {
                for (int target = 0; target < size; ++target) {
                    const int same = fuses_half(
                        level,
                        label,
                        source,
                        target
                    ) ? 1 : 0;
                    const int crossed = fuses_half(
                        level,
                        label,
                        source,
                        level - target
                    ) ? 1 : 0;
                    plus[static_cast<std::size_t>(source)]
                        [static_cast<std::size_t>(target)] =
                            same + crossed;
                    minus[static_cast<std::size_t>(source)]
                         [static_cast<std::size_t>(target)] =
                            same - crossed;
                }
            }
            std::vector<std::vector<Integer>> plus_rows(
                7U,
                std::vector<Integer>(static_cast<std::size_t>(size))
            );
            plus_rows[0][0] = 1;
            for (int power = 1; power <= 6; ++power) {
                plus_rows[static_cast<std::size_t>(power)] =
                    multiply_row(
                        plus_rows[static_cast<std::size_t>(power - 1)],
                        plus
                    );
            }
            for (int source = 0; source < size; ++source) {
                for (int target = 0; target < size; ++target) {
                    const int expected =
                        target == label - source ? 1 : 0;
                    if (
                        minus[static_cast<std::size_t>(source)]
                             [static_cast<std::size_t>(target)]
                        != expected
                    ) {
                        throw std::runtime_error(
                            "nearest odd quotient is not the simple current"
                        );
                    }
                }
            }

            const Integer f4 = 2 * label + 1;
            const Integer f5 = 2 * Integer(label) * (label + 2);
            for (int target = 0; target <= label; ++target) {
                ++formula_rows;
                std::vector<bool> activation_mask;
                if (print_slack_certificate && target >= 1) {
                    activation_recorder = &activation_mask;
                }
                std::vector<Integer> prefix(6U);
                for (int power = 0; power <= 5; ++power) {
                    prefix[static_cast<std::size_t>(power)] =
                        crossing_column_prefix(power, label, target);
                    Integer recurrence = 0;
                    for (int source = 0; source < size; ++source) {
                        const int crossing =
                            plus[static_cast<std::size_t>(source)]
                                [static_cast<std::size_t>(target)]
                            - minus[static_cast<std::size_t>(source)]
                                   [static_cast<std::size_t>(target)];
                        recurrence += plus_rows[
                            static_cast<std::size_t>(power)
                        ][static_cast<std::size_t>(source)] * crossing;
                    }
                    ++recurrence_checks;
                    if (
                        prefix[static_cast<std::size_t>(power)]
                        != recurrence
                    ) {
                        throw std::runtime_error(
                            "closed formula disagrees with quotient recurrence"
                        );
                    }
                }
                activation_recorder = nullptr;
                if (print_slack_certificate && target >= 1) {
                    activation_masks.insert(activation_mask);
                }
                const Integer same =
                    f4 * (prefix[1] + prefix[3] + prefix[5])
                    - f5 * (prefix[2] + prefix[4]);
                const Integer reflected =
                    f4 * (prefix[2] + prefix[4])
                    - f5 * (prefix[1] + prefix[3]);
                if (print_table && label == maximum_label) {
                    std::cout
                        << "NEAREST_COLUMN"
                        << " Q=" << label
                        << " V=" << target;
                    for (int power = 1; power <= 5; ++power) {
                        std::cout
                            << " d" << power << '='
                            << prefix[static_cast<std::size_t>(power)];
                    }
                    std::cout
                        << " same=" << same
                        << " reflected=" << reflected << '\n';
                }
                if (!initialized_minimum) {
                    minimum_same = same;
                    minimum_reflected = reflected;
                    initialized_minimum = true;
                } else {
                    minimum_same = std::min(minimum_same, same);
                    minimum_reflected =
                        std::min(minimum_reflected, reflected);
                }
                if (same < 0) {
                    ++negative_same_endpoint;
                }
                if (reflected < 0) {
                    ++negative_reflected_endpoint;
                }
                if (print_slack_certificate && target >= 1) {
                    ++polynomial_identity_checks;
                    if (
                        evaluate_polynomial(
                            same_polynomial,
                            label,
                            target
                        ) != Rational(same)
                        || evaluate_polynomial(
                            reflected_polynomial,
                            label,
                            target
                        ) != Rational(reflected)
                    ) {
                        std::cout
                            << "FIRST_POLYNOMIAL_IDENTITY_MISMATCH"
                            << " Q=" << label
                            << " V=" << target
                            << " same=" << same
                            << " same_polynomial="
                            << evaluate_polynomial(
                                same_polynomial,
                                label,
                                target
                            )
                            << " reflected=" << reflected
                            << " reflected_polynomial="
                            << evaluate_polynomial(
                                reflected_polynomial,
                                label,
                                target
                            )
                            << '\n';
                        throw std::runtime_error(
                            "slack polynomial identity mismatch"
                        );
                    }
                }
            }
        }

        if (print_slack_certificate) {
            print_certificate(
                "same",
                pull_back_to_slacks(same_polynomial)
            );
            print_certificate(
                "reflected",
                pull_back_to_slacks(reflected_polynomial)
            );
        }

        std::cout
            << "SU2_SHELL_NEAREST_FORMULA"
            << " maximum_label=" << maximum_label
            << " formula_rows=" << formula_rows
            << " recurrence_checks=" << recurrence_checks
            << " polynomial_identity_checks="
            << polynomial_identity_checks
            << " activation_masks=" << activation_masks.size()
            << " negative_same_endpoint=" << negative_same_endpoint
            << " negative_reflected_endpoint="
            << negative_reflected_endpoint
            << " minimum_same=" << minimum_same
            << " minimum_reflected=" << minimum_reflected
            << " result="
            << (
                negative_same_endpoint == 0U
                    && negative_reflected_endpoint == 0U
                    ? "PASS_NEAREST_COLUMN_DISCOVERY"
                    : "FAIL_NEAREST_COLUMN_CANDIDATE"
            )
            << '\n';
        return negative_same_endpoint == 0U
                && negative_reflected_endpoint == 0U
            ? 0
            : 1;
    } catch (const std::exception& exception) {
        std::cerr << "error: " << exception.what() << '\n';
        return 2;
    }
}
