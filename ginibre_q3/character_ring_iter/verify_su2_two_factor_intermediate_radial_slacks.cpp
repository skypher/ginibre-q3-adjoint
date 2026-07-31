#include <z3++.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using Vector = std::array<Rational, 4>;
using Exponent = std::array<int, 4>;
using Signature = std::array<int, 6>;

class Polynomial {
public:
    Polynomial() = default;

    explicit Polynomial(const Rational& constant) {
        if (constant != 0) {
            terms_[Exponent{0, 0, 0, 0}] = constant;
        }
    }

    static Polynomial variable(const std::size_t index) {
        if (index >= 4U) {
            throw std::runtime_error("polynomial variable out of range");
        }
        Polynomial result;
        Exponent exponent{0, 0, 0, 0};
        exponent[index] = 1;
        result.terms_[exponent] = 1;
        return result;
    }

    Polynomial& operator+=(const Polynomial& other) {
        for (const auto& [exponent, coefficient] : other.terms_) {
            terms_[exponent] += coefficient;
            if (terms_[exponent] == 0) {
                terms_.erase(exponent);
            }
        }
        return *this;
    }

    Polynomial& operator-=(const Polynomial& other) {
        for (const auto& [exponent, coefficient] : other.terms_) {
            terms_[exponent] -= coefficient;
            if (terms_[exponent] == 0) {
                terms_.erase(exponent);
            }
        }
        return *this;
    }

    Polynomial& operator*=(const Polynomial& other) {
        std::map<Exponent, Rational> product;
        for (const auto& [left_exponent, left_coefficient] : terms_) {
            for (const auto& [right_exponent, right_coefficient]
                 : other.terms_) {
                Exponent exponent{};
                for (std::size_t index = 0U;
                     index < exponent.size();
                     ++index) {
                    exponent[index]
                        = left_exponent[index] + right_exponent[index];
                }
                product[exponent]
                    += left_coefficient * right_coefficient;
            }
        }
        terms_ = std::move(product);
        return *this;
    }

    [[nodiscard]] bool nonnegative_coefficients() const {
        return std::all_of(
            terms_.begin(),
            terms_.end(),
            [](const auto& term) { return term.second >= 0; });
    }

    [[nodiscard]] std::size_t term_count() const {
        return terms_.size();
    }

    [[nodiscard]] Rational minimum_coefficient() const {
        if (terms_.empty()) {
            return Rational(0);
        }
        Rational result = terms_.begin()->second;
        for (const auto& [exponent, coefficient] : terms_) {
            static_cast<void>(exponent);
            result = std::min(result, coefficient);
        }
        return result;
    }

    [[nodiscard]] const std::map<Exponent, Rational>& terms() const {
        return terms_;
    }

private:
    std::map<Exponent, Rational> terms_;
};

Polynomial operator+(Polynomial left, const Polynomial& right) {
    left += right;
    return left;
}

Polynomial operator-(Polynomial left, const Polynomial& right) {
    left -= right;
    return left;
}

Polynomial operator*(Polynomial left, const Polynomial& right) {
    left *= right;
    return left;
}

Polynomial operator*(const Rational& scalar, Polynomial value) {
    value *= Polynomial(scalar);
    return value;
}

struct Affine {
    Rational constant{0};
    Vector coefficients{Rational(0), Rational(0), Rational(0), Rational(0)};
};

struct Inequality {
    std::string name;
    Affine slack;
};

Affine add(const Affine& left, const Affine& right) {
    Affine result;
    result.constant = left.constant + right.constant;
    for (std::size_t index = 0U; index < 4U; ++index) {
        result.coefficients[index]
            = left.coefficients[index] + right.coefficients[index];
    }
    return result;
}

Affine subtract(const Affine& left, const Affine& right) {
    Affine result;
    result.constant = left.constant - right.constant;
    for (std::size_t index = 0U; index < 4U; ++index) {
        result.coefficients[index]
            = left.coefficients[index] - right.coefficients[index];
    }
    return result;
}

Affine scale(const Rational& scalar, const Affine& value) {
    Affine result;
    result.constant = scalar * value.constant;
    for (std::size_t index = 0U; index < 4U; ++index) {
        result.coefficients[index]
            = scalar * value.coefficients[index];
    }
    return result;
}

Affine constant(const long long value) {
    Affine result;
    result.constant = Rational(value);
    return result;
}

Affine variable(const std::size_t index) {
    Affine result;
    result.coefficients[index] = 1;
    return result;
}

void append_inequality(
    std::vector<Inequality>& inequalities,
    const std::string& name,
    const Affine& slack) {
    inequalities.push_back(Inequality{name, slack});
}

std::vector<Inequality> build_inequalities(
    const bool late_wall,
    const Signature& signature) {
    // Variables are (n,d,L,g), where n=2a, d=2u, and
    // A=2L+g+1.  Passing to this real relaxation is stronger than
    // the original even-integer problem.
    const Affine n = variable(0U);
    const Affine d = variable(1U);
    const Affine contraction = variable(2U);
    const Affine gap = variable(3U);
    const Affine antidiagonal = add(
        add(scale(Rational(2), contraction), gap),
        constant(1));
    const std::array<Affine, 6> labels{
        contraction,
        add(contraction, constant(1)),
        add(add(contraction, gap), constant(1)),
        add(add(contraction, gap), constant(2)),
        add(antidiagonal, constant(1)),
        add(antidiagonal, constant(2)),
    };

    std::vector<Inequality> inequalities;
    append_inequality(
        inequalities, "n_ge_2", subtract(n, constant(2)));
    append_inequality(
        inequalities, "d_ge_2", subtract(d, constant(2)));
    append_inequality(
        inequalities,
        "d_le_2n_minus_2",
        subtract(
            subtract(scale(Rational(2), n), d),
            constant(2)));
    append_inequality(inequalities, "L_ge_0", contraction);
    append_inequality(inequalities, "g_ge_0", gap);
    append_inequality(
        inequalities,
        "A_le_d_plus_2n_minus_2",
        subtract(
            add(d, scale(Rational(2), n)),
            add(antidiagonal, constant(2))));
    if (late_wall) {
        append_inequality(
            inequalities, "d_ge_n", subtract(d, n));
    } else {
        append_inequality(
            inequalities,
            "d_le_n_minus_2",
            subtract(subtract(n, d), constant(2)));
    }

    for (std::size_t index = 0U; index < labels.size(); ++index) {
        const Affine& label = labels[index];
        const int branch = signature[index];
        const std::string prefix
            = "t" + std::to_string(index) + "_b"
              + std::to_string(branch);
        if (late_wall) {
            if (branch == 0) {
                append_inequality(
                    inequalities,
                    prefix + "_upper",
                    subtract(n, label));
            } else if (branch == 1) {
                append_inequality(
                    inequalities,
                    prefix + "_lower",
                    subtract(subtract(label, n), constant(1)));
                append_inequality(
                    inequalities,
                    prefix + "_upper",
                    subtract(d, label));
            } else if (branch == 2) {
                append_inequality(
                    inequalities,
                    prefix + "_lower",
                    subtract(subtract(label, d), constant(1)));
                append_inequality(
                    inequalities,
                    prefix + "_upper",
                    subtract(add(d, n), label));
            } else {
                append_inequality(
                    inequalities,
                    prefix + "_lower",
                    subtract(
                        subtract(label, add(d, n)),
                        constant(1)));
                append_inequality(
                    inequalities,
                    prefix + "_upper",
                    subtract(
                        add(d, scale(Rational(2), n)),
                        label));
            }
        } else {
            if (branch == 0) {
                append_inequality(
                    inequalities,
                    prefix + "_upper",
                    subtract(d, label));
            } else if (branch == 1) {
                append_inequality(
                    inequalities,
                    prefix + "_lower",
                    subtract(subtract(label, d), constant(1)));
                append_inequality(
                    inequalities,
                    prefix + "_upper",
                    subtract(n, label));
            } else if (branch == 2) {
                append_inequality(
                    inequalities,
                    prefix + "_lower",
                    subtract(subtract(label, n), constant(1)));
                append_inequality(
                    inequalities,
                    prefix + "_upper",
                    subtract(add(d, n), label));
            } else {
                append_inequality(
                    inequalities,
                    prefix + "_lower",
                    subtract(
                        subtract(label, add(d, n)),
                        constant(1)));
                append_inequality(
                    inequalities,
                    prefix + "_upper",
                    subtract(
                        add(d, scale(Rational(2), n)),
                        label));
            }
        }
    }
    return inequalities;
}

z3::expr rational_expr(z3::context& context, const Rational& value) {
    const std::string numerator = value.numerator().convert_to<std::string>();
    const std::string denominator
        = value.denominator().convert_to<std::string>();
    if (value.denominator() == 1) {
        return context.real_val(numerator.c_str());
    }
    return context.real_val(numerator.c_str())
           / context.real_val(denominator.c_str());
}

bool feasible(
    const std::vector<Inequality>& inequalities,
    z3::context& context) {
    const std::array<z3::expr, 4> variables{
        context.real_const("n"),
        context.real_const("d"),
        context.real_const("L"),
        context.real_const("g"),
    };
    z3::solver solver(context);
    for (const Inequality& inequality : inequalities) {
        z3::expr expression
            = rational_expr(context, inequality.slack.constant);
        for (std::size_t index = 0U; index < 4U; ++index) {
            expression
                = expression
                  + rational_expr(
                        context,
                        inequality.slack.coefficients[index])
                        * variables[index];
        }
        solver.add(expression >= 0);
    }
    const z3::check_result result = solver.check();
    if (result == z3::unknown) {
        throw std::runtime_error(
            "exact linear feasibility returned unknown");
    }
    return result == z3::sat;
}

bool can_be_positive(
    const std::vector<Inequality>& inequalities,
    const std::size_t selected,
    z3::context& context) {
    const std::array<z3::expr, 4> variables{
        context.real_const("n"),
        context.real_const("d"),
        context.real_const("L"),
        context.real_const("g"),
    };
    z3::solver solver(context);
    std::vector<z3::expr> expressions;
    expressions.reserve(inequalities.size());
    for (const Inequality& inequality : inequalities) {
        z3::expr expression
            = rational_expr(context, inequality.slack.constant);
        for (std::size_t index = 0U; index < 4U; ++index) {
            expression
                = expression
                  + rational_expr(
                        context,
                        inequality.slack.coefficients[index])
                        * variables[index];
        }
        expressions.push_back(expression);
        solver.add(expression >= 0);
    }
    solver.add(expressions[selected] > 0);
    const z3::check_result result = solver.check();
    if (result == z3::unknown) {
        throw std::runtime_error(
            "exact affine-hull query returned unknown");
    }
    return result == z3::sat;
}

using Matrix = std::array<std::array<Rational, 8>, 4>;

bool inverse(
    const std::array<Affine, 4>& rows,
    std::array<std::array<Rational, 4>, 4>& result) {
    Matrix matrix{};
    for (std::size_t row = 0U; row < 4U; ++row) {
        for (std::size_t column = 0U; column < 4U; ++column) {
            matrix[row][column] = rows[row].coefficients[column];
        }
        matrix[row][row + 4U] = 1;
    }
    for (std::size_t column = 0U; column < 4U; ++column) {
        std::size_t pivot = column;
        while (pivot < 4U && matrix[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == 4U) {
            return false;
        }
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
        }
        const Rational divisor = matrix[column][column];
        for (std::size_t entry = 0U; entry < 8U; ++entry) {
            matrix[column][entry] /= divisor;
        }
        for (std::size_t row = 0U; row < 4U; ++row) {
            if (row == column) {
                continue;
            }
            const Rational multiplier = matrix[row][column];
            for (std::size_t entry = 0U; entry < 8U; ++entry) {
                matrix[row][entry]
                    -= multiplier * matrix[column][entry];
            }
        }
    }
    for (std::size_t row = 0U; row < 4U; ++row) {
        for (std::size_t column = 0U; column < 4U; ++column) {
            result[row][column] = matrix[row][column + 4U];
        }
    }
    return true;
}

std::size_t row_rank(const std::vector<Affine>& rows) {
    std::vector<Vector> matrix;
    matrix.reserve(rows.size());
    for (const Affine& row : rows) {
        matrix.push_back(row.coefficients);
    }
    std::size_t rank = 0U;
    for (std::size_t column = 0U;
         column < 4U && rank < matrix.size();
         ++column) {
        std::size_t pivot = rank;
        while (pivot < matrix.size() && matrix[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == matrix.size()) {
            continue;
        }
        std::swap(matrix[pivot], matrix[rank]);
        const Rational divisor = matrix[rank][column];
        for (std::size_t entry = column; entry < 4U; ++entry) {
            matrix[rank][entry] /= divisor;
        }
        for (std::size_t row = 0U; row < matrix.size(); ++row) {
            if (row == rank) {
                continue;
            }
            const Rational multiplier = matrix[row][column];
            for (std::size_t entry = column; entry < 4U; ++entry) {
                matrix[row][entry]
                    -= multiplier * matrix[rank][entry];
            }
        }
        ++rank;
    }
    return rank;
}

bool orthant_map(
    const std::vector<std::size_t>& equality_basis,
    const std::vector<std::size_t>& slack_basis,
    const std::vector<Inequality>& inequalities,
    std::array<Polynomial, 4>& variables) {
    if (equality_basis.size() + slack_basis.size() != 4U) {
        throw std::runtime_error("orthant basis has wrong size");
    }
    std::array<Affine, 4> rows{};
    std::size_t row = 0U;
    for (const std::size_t index : equality_basis) {
        rows[row] = inequalities[index].slack;
        ++row;
    }
    for (const std::size_t index : slack_basis) {
        rows[row] = inequalities[index].slack;
        ++row;
    }
    std::array<std::array<Rational, 4>, 4> matrix_inverse{};
    if (!inverse(rows, matrix_inverse)) {
        return false;
    }
    for (std::size_t variable_index = 0U;
         variable_index < 4U;
         ++variable_index) {
        Polynomial value(Rational(0));
        for (std::size_t row_index = 0U;
             row_index < 4U;
             ++row_index) {
            Polynomial right_hand_side(-rows[row_index].constant);
            if (row_index >= equality_basis.size()) {
                right_hand_side += Polynomial::variable(
                    row_index - equality_basis.size());
            }
            value += matrix_inverse[variable_index][row_index]
                     * right_hand_side;
        }
        variables[variable_index] = value;
    }
    for (const Inequality& inequality : inequalities) {
        Polynomial slack(inequality.slack.constant);
        for (std::size_t index = 0U; index < 4U; ++index) {
            slack += inequality.slack.coefficients[index]
                     * variables[index];
        }
        if (!slack.nonnegative_coefficients()) {
            return false;
        }
    }
    return true;
}

Polynomial twice_coefficient(
    const bool late_wall,
    const int branch,
    const Polynomial& n,
    const Polynomial& d,
    const Polynomial& label) {
    const Polynomial size = n + Polynomial(Rational(1));
    const Polynomial band
        = Rational(2) * size
              * (Rational(2) * label + Polynomial(Rational(1)))
          - Rational(2) * label
                * (label + Polynomial(Rational(1)));
    const Polynomial wall_offset = label - d;
    if (late_wall) {
        if (branch == 0) {
            return band;
        }
        if (branch == 1) {
            return Rational(2) * size * size;
        }
        if (branch == 2) {
            return Rational(2) * size * size
                   - wall_offset
                         * (
                             wall_offset
                             + Polynomial(Rational(1)));
        }
    } else {
        if (branch == 0) {
            return band;
        }
        if (branch == 1) {
            return band
                   - wall_offset
                         * (
                             wall_offset
                             + Polynomial(Rational(1)));
        }
        if (branch == 2) {
            return Rational(2) * size * size
                   - wall_offset
                         * (
                             wall_offset
                             + Polynomial(Rational(1)));
        }
    }
    const Polynomial distance
        = d + Rational(2) * n - label;
    return (distance + Polynomial(Rational(1)))
           * (distance + Polynomial(Rational(2)));
}

Polynomial radial(
    const bool late_wall,
    const Signature& signature,
    const std::array<Polynomial, 4>& variables) {
    const Polynomial& n = variables[0];
    const Polynomial& d = variables[1];
    const Polynomial& contraction = variables[2];
    const Polynomial& gap = variables[3];
    const Polynomial antidiagonal
        = Rational(2) * contraction + gap
          + Polynomial(Rational(1));
    const std::array<Polynomial, 6> labels{
        contraction,
        contraction + Polynomial(Rational(1)),
        contraction + gap + Polynomial(Rational(1)),
        contraction + gap + Polynomial(Rational(2)),
        antidiagonal + Polynomial(Rational(1)),
        antidiagonal + Polynomial(Rational(2)),
    };
    std::array<Polynomial, 6> coefficients{};
    for (std::size_t index = 0U; index < labels.size(); ++index) {
        coefficients[index] = twice_coefficient(
            late_wall,
            signature[index],
            n,
            d,
            labels[index]);
    }
    const Polynomial c_zero
        = Rational(2) * (n + Polynomial(Rational(1)));
    return c_zero * (coefficients[4] + coefficients[5])
           + coefficients[0] * coefficients[2]
           - coefficients[1] * coefficients[3];
}

std::string render(const Signature& signature) {
    std::string result;
    for (const int branch : signature) {
        result += static_cast<char>('0' + branch);
    }
    return result;
}

void generate_signatures(
    const std::size_t position,
    const int minimum_branch,
    Signature& signature,
    std::vector<Signature>& signatures) {
    if (position == signature.size()) {
        signatures.push_back(signature);
        return;
    }
    for (int branch = minimum_branch; branch <= 3; ++branch) {
        signature[position] = branch;
        generate_signatures(
            position + 1U, branch, signature, signatures);
    }
}

std::vector<std::array<Polynomial, 4>> late_mixed_maps() {
    std::vector<std::array<Polynomial, 4>> result;
    const Polynomial alpha = Polynomial::variable(0U);
    const Polynomial beta = Polynomial::variable(1U);
    const Polynomial gamma = Polynomial::variable(2U);
    const Polynomial z = Polynomial::variable(3U);
    for (int base_x = 0; base_x <= 2; ++base_x) {
        const int base_w = 2 - base_x;
        for (int cone = 0; cone < 2; ++cone) {
            Polynomial x;
            Polynomial y;
            Polynomial w;
            if (cone == 0) {
                x = Polynomial(Rational(base_x)) + alpha + beta;
                y = alpha;
                w = Polynomial(Rational(base_w)) + gamma;
            } else {
                x = Polynomial(Rational(base_x)) + alpha;
                y = alpha + beta;
                w = Polynomial(Rational(base_w)) + beta + gamma;
            }
            result.push_back(std::array<Polynomial, 4>{
                x + z + w + Polynomial(Rational(2)),
                x + y + w + Rational(2) * z
                    + Polynomial(Rational(4)),
                z + w + Polynomial(Rational(1)),
                x + y + Polynomial(Rational(1)),
            });
        }
    }
    return result;
}

std::vector<std::array<Polynomial, 4>> late_transport_maps() {
    std::vector<std::array<Polynomial, 4>> result;
    const Polynomial alpha = Polynomial::variable(0U);
    const Polynomial beta = Polynomial::variable(1U);
    const Polynomial gamma = Polynomial::variable(2U);
    const Polynomial delta = Polynomial::variable(3U);
    for (int cone = 0; cone < 3; ++cone) {
        Polynomial x;
        Polynomial w;
        Polynomial z;
        Polynomial delay;
        if (cone == 0) {
            delay = alpha;
            x = alpha + beta;
            w = gamma;
            z = delta;
        } else if (cone == 1) {
            x = alpha;
            delay = alpha + beta;
            w = beta + gamma;
            z = delta;
        } else {
            x = alpha;
            w = beta;
            delay = alpha + beta + gamma;
            z = gamma + delta;
        }
        const Polynomial n
            = x + w + z + Polynomial(Rational(2));
        const Polynomial contraction
            = z + w + Polynomial(Rational(1));
        const Polynomial d = n + delay;
        const Polynomial gap
            = d + x - z - Polynomial(Rational(1));
        result.push_back(
            std::array<Polynomial, 4>{
                n, d, contraction, gap});
    }
    return result;
}

std::vector<std::array<Polynomial, 4>> late_plateau_decline_maps() {
    std::vector<std::array<Polynomial, 4>> result;
    const Polynomial alpha = Polynomial::variable(0U);
    const Polynomial beta = Polynomial::variable(1U);
    const Polynomial gamma = Polynomial::variable(2U);
    const Polynomial delta = Polynomial::variable(3U);
    for (int base_delay = 0; base_delay <= 1; ++base_delay) {
        for (int cone = 0; cone < 3; ++cone) {
            Polynomial y;
            Polynomial contraction;
            Polynomial tail_slack;
            Polynomial delay;
            if (cone == 0) {
                delay = Polynomial(Rational(base_delay)) + alpha;
                y = alpha + beta;
                contraction = gamma;
                tail_slack = delta;
            } else if (cone == 1) {
                y = alpha;
                delay = Polynomial(Rational(base_delay))
                        + alpha + beta;
                contraction = beta + gamma;
                tail_slack = delta;
            } else {
                y = alpha;
                contraction = beta;
                delay = Polynomial(Rational(base_delay))
                        + alpha + beta + gamma;
                tail_slack = gamma + delta;
            }
            const Polynomial z
                = contraction + tail_slack
                  + Polynomial(Rational(1));
            const Polynomial n
                = y + contraction + tail_slack
                  + Polynomial(Rational(3));
            const Polynomial d = n + delay;
            const Polynomial gap = d - contraction + y;
            result.push_back(
                std::array<Polynomial, 4>{
                    n, d, contraction, gap});
        }
    }
    return result;
}

std::vector<std::array<Polynomial, 4>> late_boundary_tail_maps() {
    std::vector<std::array<Polynomial, 4>> result;
    const Polynomial alpha = Polynomial::variable(0U);
    const Polynomial beta = Polynomial::variable(1U);
    const Polynomial gamma = Polynomial::variable(2U);
    for (int cone = 0; cone < 2; ++cone) {
        Polynomial contraction;
        Polynomial support_slack;
        Polynomial delay;
        if (cone == 0) {
            delay = alpha;
            contraction = alpha + beta;
            support_slack = gamma;
        } else {
            contraction = alpha;
            delay = alpha + beta;
            support_slack = beta + gamma;
        }
        const Polynomial n
            = contraction + support_slack
              + Polynomial(Rational(2));
        const Polynomial d = n + delay;
        const Polynomial gap
            = d + n - contraction - Polynomial(Rational(1));
        result.push_back(
            std::array<Polynomial, 4>{
                n, d, contraction, gap});
    }
    return result;
}

std::vector<std::array<Polynomial, 4>> late_lower_boundary_maps() {
    std::vector<std::array<Polynomial, 4>> result;
    const Polynomial alpha = Polynomial::variable(0U);
    const Polynomial beta = Polynomial::variable(1U);
    const Polynomial gamma = Polynomial::variable(2U);
    for (int cone = 0; cone < 2; ++cone) {
        Polynomial support;
        Polynomial delay;
        Polynomial excess;
        if (cone == 0) {
            excess = alpha;
            delay = alpha + beta + Polynomial(Rational(1));
            support = alpha + beta + gamma
                      + Polynomial(Rational(3));
        } else {
            delay = alpha + Polynomial(Rational(1));
            excess = alpha + beta + Polynomial(Rational(1));
            support = alpha + beta + gamma
                      + Polynomial(Rational(4));
        }
        const Polynomial d = support + delay;
        const Polynomial gap = delay + excess;
        result.push_back(
            std::array<Polynomial, 4>{
                support, d, support, gap});
    }
    return result;
}

std::vector<std::array<Polynomial, 4>> late_inner_transport_maps() {
    std::vector<std::array<Polynomial, 4>> result;
    const Polynomial alpha = Polynomial::variable(0U);
    const Polynomial beta = Polynomial::variable(1U);
    const Polynomial gamma = Polynomial::variable(2U);
    const Polynomial x = Polynomial::variable(3U);
    for (int cone = 0; cone < 2; ++cone) {
        Polynomial y;
        Polynomial w;
        Polynomial excess_delay;
        if (cone == 0) {
            excess_delay = alpha;
            y = alpha + beta;
            w = gamma;
        } else {
            y = alpha;
            excess_delay = alpha + beta;
            w = beta + gamma;
        }
        const Polynomial n
            = x + y + w + Polynomial(Rational(4));
        const Polynomial contraction
            = n + x + Polynomial(Rational(1));
        const Polynomial delay
            = x + excess_delay + Polynomial(Rational(2));
        const Polynomial d = n + delay;
        const Polynomial gap
            = d - contraction + y;
        result.push_back(
            std::array<Polynomial, 4>{
                n, d, contraction, gap});
    }
    return result;
}

struct OrthantCertificate {
    std::vector<std::size_t> slack_basis;
    std::size_t terms = 0U;
    Rational minimum_coefficient{0};
};

struct LeafCertificate {
    std::string path;
    std::size_t dimension = 0U;
    std::vector<std::string> equalities;
    std::vector<std::string> basis;
    std::size_t terms = 0U;
    Rational minimum_coefficient{0};
};

struct CellStatistics {
    std::size_t nodes = 0U;
    std::size_t empty_nodes = 0U;
    std::size_t bases_examined = 0U;
    std::size_t maximum_depth = 0U;
};

bool certify_late_mixed_cell(
    const std::vector<Inequality>& inequalities,
    const Signature& signature,
    const std::string& path,
    std::vector<LeafCertificate>& leaves) {
    if (render(signature) != "001122") {
        return false;
    }
    const auto maps = late_mixed_maps();
    if (maps.size() != 6U) {
        throw std::runtime_error("unexpected mixed-cell map count");
    }
    for (std::size_t map_index = 0U;
         map_index < maps.size();
         ++map_index) {
        const auto& variables = maps[map_index];
        for (const Inequality& inequality : inequalities) {
            Polynomial slack(inequality.slack.constant);
            for (std::size_t index = 0U; index < 4U; ++index) {
                slack += inequality.slack.coefficients[index]
                         * variables[index];
            }
            if (!slack.nonnegative_coefficients()) {
                throw std::runtime_error(
                    "mixed-cell map leaves branch polyhedron");
            }
        }
        const Polynomial margin = radial(
            true, signature, variables);
        if (!margin.nonnegative_coefficients()) {
            throw std::runtime_error(
                "mixed-cell margin has a negative coefficient");
        }
        LeafCertificate leaf;
        leaf.path = path + "/semigroup_map_"
                    + std::to_string(map_index);
        leaf.dimension = 4U;
        leaf.basis = {
            "alpha", "beta", "gamma", "z"};
        leaf.terms = margin.term_count();
        leaf.minimum_coefficient = margin.minimum_coefficient();
        leaves.push_back(std::move(leaf));
    }
    return true;
}

bool certify_late_transport_cell(
    const std::vector<Inequality>& inequalities,
    const Signature& signature,
    const std::string& path,
    std::vector<LeafCertificate>& leaves) {
    if (render(signature) != "002233") {
        return false;
    }
    const auto maps = late_transport_maps();
    if (maps.size() != 3U) {
        throw std::runtime_error("unexpected transport-cell map count");
    }
    for (std::size_t map_index = 0U;
         map_index < maps.size();
         ++map_index) {
        const auto& variables = maps[map_index];
        for (const Inequality& inequality : inequalities) {
            Polynomial slack(inequality.slack.constant);
            for (std::size_t index = 0U; index < 4U; ++index) {
                slack += inequality.slack.coefficients[index]
                         * variables[index];
            }
            if (!slack.nonnegative_coefficients()) {
                throw std::runtime_error(
                    "transport-cell map leaves branch polyhedron");
            }
        }
        const Polynomial margin = radial(
            true, signature, variables);
        if (!margin.nonnegative_coefficients()) {
            throw std::runtime_error(
                "transport-cell margin has a negative coefficient");
        }
        LeafCertificate leaf;
        leaf.path = path + "/transport_map_"
                    + std::to_string(map_index);
        leaf.dimension = 4U;
        leaf.basis = {
            "alpha", "beta", "gamma", "delta"};
        leaf.terms = margin.term_count();
        leaf.minimum_coefficient = margin.minimum_coefficient();
        leaves.push_back(std::move(leaf));
    }
    return true;
}

bool certify_late_plateau_decline_cell(
    const std::vector<Inequality>& inequalities,
    const Signature& signature,
    const std::string& path,
    std::vector<LeafCertificate>& leaves) {
    if (render(signature) != "002222") {
        return false;
    }
    const auto maps = late_plateau_decline_maps();
    if (maps.size() != 6U) {
        throw std::runtime_error(
            "unexpected plateau-decline map count");
    }
    for (std::size_t map_index = 0U;
         map_index < maps.size();
         ++map_index) {
        const auto& variables = maps[map_index];
        for (const Inequality& inequality : inequalities) {
            Polynomial slack(inequality.slack.constant);
            for (std::size_t index = 0U; index < 4U; ++index) {
                slack += inequality.slack.coefficients[index]
                         * variables[index];
            }
            if (!slack.nonnegative_coefficients()) {
                throw std::runtime_error(
                    "plateau-decline map leaves branch polyhedron");
            }
        }
        const Polynomial margin = radial(
            true, signature, variables);
        if (!margin.nonnegative_coefficients()) {
            throw std::runtime_error(
                "plateau-decline margin has a negative coefficient");
        }
        LeafCertificate leaf;
        leaf.path = path + "/plateau_decline_map_"
                    + std::to_string(map_index);
        leaf.dimension = 4U;
        leaf.basis = {
            "alpha", "beta", "gamma", "delta"};
        leaf.terms = margin.term_count();
        leaf.minimum_coefficient = margin.minimum_coefficient();
        leaves.push_back(std::move(leaf));
    }
    return true;
}

bool certify_late_boundary_tail_cell(
    const std::vector<Inequality>& inequalities,
    const Signature& signature,
    const std::string& path,
    std::vector<LeafCertificate>& leaves) {
    if (render(signature) != "002333") {
        return false;
    }
    const auto maps = late_boundary_tail_maps();
    if (maps.size() != 2U) {
        throw std::runtime_error(
            "unexpected boundary-tail map count");
    }
    for (std::size_t map_index = 0U;
         map_index < maps.size();
         ++map_index) {
        const auto& variables = maps[map_index];
        for (const Inequality& inequality : inequalities) {
            Polynomial slack(inequality.slack.constant);
            for (std::size_t index = 0U; index < 4U; ++index) {
                slack += inequality.slack.coefficients[index]
                         * variables[index];
            }
            if (!slack.nonnegative_coefficients()) {
                throw std::runtime_error(
                    "boundary-tail map leaves branch polyhedron");
            }
        }
        const Polynomial margin = radial(
            true, signature, variables);
        if (!margin.nonnegative_coefficients()) {
            throw std::runtime_error(
                "boundary-tail margin has a negative coefficient");
        }
        LeafCertificate leaf;
        leaf.path = path + "/boundary_tail_map_"
                    + std::to_string(map_index);
        leaf.dimension = 3U;
        leaf.basis = {"alpha", "beta", "gamma"};
        leaf.terms = margin.term_count();
        leaf.minimum_coefficient = margin.minimum_coefficient();
        leaves.push_back(std::move(leaf));
    }
    return true;
}

bool certify_late_lower_boundary_cell(
    const std::vector<Inequality>& inequalities,
    const Signature& signature,
    const std::string& path,
    std::vector<LeafCertificate>& leaves) {
    if (render(signature) != "012233") {
        return false;
    }
    const auto maps = late_lower_boundary_maps();
    if (maps.size() != 2U) {
        throw std::runtime_error(
            "unexpected lower-boundary map count");
    }
    for (std::size_t map_index = 0U;
         map_index < maps.size();
         ++map_index) {
        const auto& variables = maps[map_index];
        for (const Inequality& inequality : inequalities) {
            Polynomial slack(inequality.slack.constant);
            for (std::size_t index = 0U; index < 4U; ++index) {
                slack += inequality.slack.coefficients[index]
                         * variables[index];
            }
            if (!slack.nonnegative_coefficients()) {
                throw std::runtime_error(
                    "lower-boundary map leaves branch polyhedron");
            }
        }
        const Polynomial margin = radial(
            true, signature, variables);
        if (!margin.nonnegative_coefficients()) {
            throw std::runtime_error(
                "lower-boundary margin has a negative coefficient");
        }
        LeafCertificate leaf;
        leaf.path = path + "/lower_boundary_map_"
                    + std::to_string(map_index);
        leaf.dimension = 3U;
        leaf.basis = {"alpha", "beta", "gamma"};
        leaf.terms = margin.term_count();
        leaf.minimum_coefficient = margin.minimum_coefficient();
        leaves.push_back(std::move(leaf));
    }
    return true;
}

bool certify_late_inner_transport_cell(
    const std::vector<Inequality>& inequalities,
    const Signature& signature,
    const std::string& path,
    std::vector<LeafCertificate>& leaves) {
    if (render(signature) != "112233") {
        return false;
    }
    const auto maps = late_inner_transport_maps();
    if (maps.size() != 2U) {
        throw std::runtime_error(
            "unexpected inner-transport map count");
    }
    for (std::size_t map_index = 0U;
         map_index < maps.size();
         ++map_index) {
        const auto& variables = maps[map_index];
        for (const Inequality& inequality : inequalities) {
            Polynomial slack(inequality.slack.constant);
            for (std::size_t index = 0U; index < 4U; ++index) {
                slack += inequality.slack.coefficients[index]
                         * variables[index];
            }
            if (!slack.nonnegative_coefficients()) {
                throw std::runtime_error(
                    "inner-transport map leaves branch polyhedron");
            }
        }
        const Polynomial margin = radial(
            true, signature, variables);
        if (!margin.nonnegative_coefficients()) {
            throw std::runtime_error(
                "inner-transport margin has a negative coefficient");
        }
        LeafCertificate leaf;
        leaf.path = path + "/inner_transport_map_"
                    + std::to_string(map_index);
        leaf.dimension = 4U;
        leaf.basis = {
            "alpha", "beta", "gamma", "x"};
        leaf.terms = margin.term_count();
        leaf.minimum_coefficient = margin.minimum_coefficient();
        leaves.push_back(std::move(leaf));
    }
    return true;
}

bool search_orthant_certificate(
    const std::size_t start,
    const std::size_t dimension,
    const std::vector<std::size_t>& candidates,
    const std::vector<std::size_t>& equality_basis,
    const std::vector<Inequality>& inequalities,
    const bool late_wall,
    const Signature& signature,
    std::vector<std::size_t>& selected,
    std::size_t& bases_examined,
    OrthantCertificate& certificate) {
    if (selected.size() == dimension) {
        ++bases_examined;
        std::array<Polynomial, 4> variables{};
        if (!orthant_map(
                equality_basis,
                selected,
                inequalities,
                variables)) {
            return false;
        }
        const Polynomial margin = radial(
            late_wall, signature, variables);
        if (!margin.nonnegative_coefficients()) {
            return false;
        }
        certificate.slack_basis = selected;
        certificate.terms = margin.term_count();
        certificate.minimum_coefficient
            = margin.minimum_coefficient();
        return true;
    }
    const std::size_t remaining = dimension - selected.size();
    for (std::size_t index = start;
         index + remaining <= candidates.size();
         ++index) {
        selected.push_back(candidates[index]);
        if (search_orthant_certificate(
                index + 1U,
                dimension,
                candidates,
                equality_basis,
                inequalities,
                late_wall,
                signature,
                selected,
                bases_examined,
                certificate)) {
            return true;
        }
        selected.pop_back();
    }
    return false;
}

bool direct_certificate(
    const std::vector<Inequality>& inequalities,
    const bool late_wall,
    const Signature& signature,
    const std::string& path,
    LeafCertificate& leaf,
    CellStatistics& statistics) {
    z3::context context;
    std::vector<std::size_t> equality_basis;
    std::vector<Affine> equality_rows;
    std::vector<std::size_t> candidates;
    for (std::size_t index = 0U;
         index < inequalities.size();
         ++index) {
        if (can_be_positive(inequalities, index, context)) {
            candidates.push_back(index);
            continue;
        }
        std::vector<Affine> enlarged = equality_rows;
        enlarged.push_back(inequalities[index].slack);
        if (row_rank(enlarged) > equality_rows.size()) {
            equality_basis.push_back(index);
            equality_rows.push_back(inequalities[index].slack);
        }
    }
    const std::size_t dimension = 4U - equality_basis.size();
    std::vector<std::size_t> selected;
    OrthantCertificate certificate;
    std::size_t bases_examined = 0U;
    const bool found = search_orthant_certificate(
        0U,
        dimension,
        candidates,
        equality_basis,
        inequalities,
        late_wall,
        signature,
        selected,
        bases_examined,
        certificate);
    statistics.bases_examined += bases_examined;
    if (!found) {
        return false;
    }
    leaf.path = path;
    leaf.dimension = dimension;
    for (const std::size_t index : equality_basis) {
        leaf.equalities.push_back(inequalities[index].name);
    }
    for (const std::size_t index : certificate.slack_basis) {
        leaf.basis.push_back(inequalities[index].name);
    }
    leaf.terms = certificate.terms;
    leaf.minimum_coefficient = certificate.minimum_coefficient;
    return true;
}

bool projected_certificate(
    const std::vector<Inequality>& inequalities,
    const bool late_wall,
    const Signature& signature,
    const std::string& path,
    LeafCertificate& leaf,
    CellStatistics& statistics) {
    std::array<Polynomial, 4> identity{
        Polynomial::variable(0U),
        Polynomial::variable(1U),
        Polynomial::variable(2U),
        Polynomial::variable(3U),
    };
    const Polynomial margin = radial(
        late_wall, signature, identity);
    std::array<bool, 4> active{false, false, false, false};
    for (const auto& [exponent, coefficient] : margin.terms()) {
        static_cast<void>(coefficient);
        for (std::size_t index = 0U; index < 4U; ++index) {
            active[index] = active[index] || exponent[index] != 0;
        }
    }
    if (std::all_of(active.begin(), active.end(), [](const bool value) {
            return value;
        })) {
        return false;
    }

    std::vector<Inequality> projected;
    for (const Inequality& inequality : inequalities) {
        bool depends_on_inactive = false;
        for (std::size_t index = 0U; index < 4U; ++index) {
            depends_on_inactive
                = depends_on_inactive
                  || (
                      !active[index]
                      && inequality.slack.coefficients[index] != 0);
        }
        if (!depends_on_inactive) {
            projected.push_back(inequality);
        }
    }
    for (std::size_t index = 0U; index < 4U; ++index) {
        if (active[index]) {
            continue;
        }
        const Affine coordinate = variable(index);
        projected.push_back(Inequality{
            "project_v" + std::to_string(index) + "_ge_0",
            coordinate});
        projected.push_back(Inequality{
            "project_v" + std::to_string(index) + "_le_0",
            scale(Rational(-1), coordinate)});
    }
    return direct_certificate(
        projected,
        late_wall,
        signature,
        path + "/inactive_projection",
        leaf,
        statistics);
}

bool integer_affine(const Affine& value) {
    if (value.constant.denominator() != 1) {
        return false;
    }
    return std::all_of(
        value.coefficients.begin(),
        value.coefficients.end(),
        [](const Rational& coefficient) {
            return coefficient.denominator() == 1;
        });
}

std::vector<std::size_t> split_order(
    const std::vector<Inequality>& inequalities) {
    std::vector<std::size_t> result;
    const std::array<std::string, 4> preferred{
        "L_ge_0", "g_ge_0", "d_ge_2", "n_ge_2"};
    for (const std::string& name : preferred) {
        for (std::size_t index = 0U;
             index < inequalities.size();
             ++index) {
            if (inequalities[index].name == name) {
                result.push_back(index);
            }
        }
    }
    for (std::size_t index = 0U;
         index < inequalities.size();
         ++index) {
        if (std::find(result.begin(), result.end(), index)
            == result.end()) {
            result.push_back(index);
        }
    }
    return result;
}

bool certify_integer_cell(
    const std::vector<Inequality>& inequalities,
    const bool late_wall,
    const Signature& signature,
    const std::size_t depth,
    const std::size_t maximum_depth,
    const std::string& path,
    std::vector<LeafCertificate>& leaves,
    CellStatistics& statistics) {
    ++statistics.nodes;
    statistics.maximum_depth
        = std::max(statistics.maximum_depth, depth);
    {
        z3::context context;
        if (!feasible(inequalities, context)) {
            ++statistics.empty_nodes;
            return true;
        }
    }

    if (
        late_wall && depth == 0U
        && certify_late_mixed_cell(
            inequalities, signature, path, leaves)) {
        return true;
    }
    if (
        late_wall && depth == 0U
        && certify_late_transport_cell(
            inequalities, signature, path, leaves)) {
        return true;
    }
    if (
        late_wall && depth == 0U
        && certify_late_plateau_decline_cell(
            inequalities, signature, path, leaves)) {
        return true;
    }
    if (
        late_wall && depth == 0U
        && certify_late_boundary_tail_cell(
            inequalities, signature, path, leaves)) {
        return true;
    }
    if (
        late_wall && depth == 0U
        && certify_late_lower_boundary_cell(
            inequalities, signature, path, leaves)) {
        return true;
    }
    if (
        late_wall && depth == 0U
        && certify_late_inner_transport_cell(
            inequalities, signature, path, leaves)) {
        return true;
    }

    LeafCertificate leaf;
    if (direct_certificate(
            inequalities,
            late_wall,
            signature,
            path,
            leaf,
            statistics)) {
        leaves.push_back(std::move(leaf));
        return true;
    }
    LeafCertificate projected_leaf;
    if (projected_certificate(
            inequalities,
            late_wall,
            signature,
            path,
            projected_leaf,
            statistics)) {
        leaves.push_back(std::move(projected_leaf));
        return true;
    }
    if (depth == maximum_depth) {
        return false;
    }

    const auto try_partition = [&](
                                   const Affine& left_slack,
                                   const std::string& left_name,
                                   const Affine& right_slack,
                                   const std::string& right_name) {
        std::vector<Inequality> left_child = inequalities;
        left_child.push_back(Inequality{left_name, left_slack});
        std::vector<Inequality> right_child = inequalities;
        right_child.push_back(Inequality{right_name, right_slack});
        {
            z3::context left_context;
            z3::context right_context;
            if (!feasible(left_child, left_context)
                || !feasible(right_child, right_context)) {
                return false;
            }
        }

        std::vector<LeafCertificate> left_leaves;
        std::vector<LeafCertificate> right_leaves;
        CellStatistics trial_statistics = statistics;
        if (!certify_integer_cell(
                left_child,
                late_wall,
                signature,
                depth + 1U,
                maximum_depth,
                path + "/" + left_name,
                left_leaves,
                trial_statistics)) {
            return false;
        }
        if (!certify_integer_cell(
                right_child,
                late_wall,
                signature,
                depth + 1U,
                maximum_depth,
                path + "/" + right_name,
                right_leaves,
                trial_statistics)) {
            return false;
        }
        statistics = trial_statistics;
        leaves.insert(
            leaves.end(), left_leaves.begin(), left_leaves.end());
        leaves.insert(
            leaves.end(),
            right_leaves.begin(),
            right_leaves.end());
        return true;
    };

    const std::vector<std::size_t> ordered
        = split_order(inequalities);
    for (const std::size_t split_index : ordered) {
        const Inequality& split = inequalities[split_index];
        if (
            (split.name != "L_ge_0" && split.name != "g_ge_0")
            || !integer_affine(split.slack)) {
            continue;
        }
        if (try_partition(
                scale(Rational(-1), split.slack),
                "split_" + split.name + "_eq_0",
                subtract(split.slack, constant(1)),
                "split_" + split.name + "_ge_1")) {
            return true;
        }
    }

    for (std::size_t left_position = 0U;
         left_position < ordered.size();
         ++left_position) {
        for (std::size_t right_position = left_position + 1U;
             right_position < ordered.size();
             ++right_position) {
            const Inequality& left
                = inequalities[ordered[left_position]];
            const Inequality& right
                = inequalities[ordered[right_position]];
            const Affine difference
                = subtract(left.slack, right.slack);
            if (!integer_affine(difference)) {
                continue;
            }
            const std::string prefix
                = "split_" + left.name + "_minus_" + right.name;
            if (try_partition(
                    difference,
                    prefix + "_ge_0",
                    subtract(
                        scale(Rational(-1), difference),
                        constant(1)),
                    prefix + "_le_minus_1")) {
                return true;
            }
        }
    }

    for (const std::size_t split_index : ordered) {
        const Inequality& split = inequalities[split_index];
        if (!integer_affine(split.slack)) {
            continue;
        }
        if (try_partition(
                scale(Rational(-1), split.slack),
                "split_" + split.name + "_eq_0",
                subtract(split.slack, constant(1)),
                "split_" + split.name + "_ge_1")) {
            return true;
        }
    }
    return false;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const bool filtered = argc == 4
                              && std::string(argv[1]) == "--cell";
        if (argc != 1 && !filtered) {
            throw std::runtime_error(
                "usage: verify_su2_two_factor_intermediate_radial_slacks "
                "[--cell early|late SIGNATURE]");
        }
        const std::string filtered_regime
            = filtered ? std::string(argv[2]) : "";
        const std::string filtered_signature
            = filtered ? std::string(argv[3]) : "";
        if (
            filtered
            && (
                (filtered_regime != "early"
                 && filtered_regime != "late")
                || filtered_signature.size() != 6U)) {
            throw std::runtime_error("invalid cell filter");
        }
        Signature partial{};
        std::vector<Signature> signatures;
        generate_signatures(0U, 0, partial, signatures);
        if (signatures.size() != 84U) {
            throw std::runtime_error("unexpected signature count");
        }

        std::size_t attempted = 0U;
        std::size_t infeasible = 0U;
        std::size_t certified = 0U;
        std::size_t total_terms = 0U;
        std::size_t total_leaves = 0U;
        std::size_t total_nodes = 0U;
        std::size_t total_bases_examined = 0U;
        std::size_t maximum_depth = 0U;
        for (const bool late_wall : {false, true}) {
            for (const Signature& signature : signatures) {
                if (
                    filtered
                    && (
                        (late_wall ? "late" : "early")
                            != filtered_regime
                        || render(signature) != filtered_signature)) {
                    continue;
                }
                ++attempted;
                const std::vector<Inequality> inequalities
                    = build_inequalities(late_wall, signature);
                z3::context context;
                if (!feasible(inequalities, context)) {
                    ++infeasible;
                    continue;
                }

                std::vector<LeafCertificate> leaves;
                CellStatistics statistics;
                const bool found = certify_integer_cell(
                    inequalities,
                    late_wall,
                    signature,
                    0U,
                    2U,
                    "root",
                    leaves,
                    statistics);
                if (!found) {
                    std::cout
                        << "UNCERTIFIED regime="
                        << (late_wall ? "n_le_d" : "d_le_n_minus_2")
                        << " signature=" << render(signature)
                        << " inequalities=" << inequalities.size()
                        << " nodes=" << statistics.nodes
                        << " bases_examined="
                        << statistics.bases_examined
                        << '\n';
                    return EXIT_FAILURE;
                }
                ++certified;
                std::size_t cell_terms = 0U;
                Rational cell_minimum(0);
                bool first_leaf = true;
                for (const LeafCertificate& leaf : leaves) {
                    cell_terms += leaf.terms;
                    if (first_leaf
                        || leaf.minimum_coefficient < cell_minimum) {
                        cell_minimum = leaf.minimum_coefficient;
                    }
                    first_leaf = false;
                }
                total_terms += cell_terms;
                total_leaves += leaves.size();
                total_nodes += statistics.nodes;
                total_bases_examined += statistics.bases_examined;
                maximum_depth = std::max(
                    maximum_depth, statistics.maximum_depth);
                std::cout
                    << "regime="
                    << (
                        late_wall
                            ? "n_le_d"
                            : "d_le_n_minus_2")
                    << " signature=" << render(signature)
                    << " leaves=" << leaves.size()
                    << " nodes=" << statistics.nodes
                    << " maximum_depth="
                    << statistics.maximum_depth
                    << " terms=" << cell_terms
                    << " minimum_coefficient="
                    << cell_minimum
                    << '\n';
            }
        }
        std::cout
            << "SU2_TWO_FACTOR_INTERMEDIATE_RADIAL_SLACKS"
            << " attempted=" << attempted
            << " infeasible=" << infeasible
            << " feasible=" << attempted - infeasible
            << " certified=" << certified
            << " leaves=" << total_leaves
            << " nodes=" << total_nodes
            << " maximum_depth=" << maximum_depth
            << " bases_examined=" << total_bases_examined
            << " terms=" << total_terms
            << " result="
            << (
                certified + infeasible == attempted
                    ? "PASS_EXACT_INTEGER_ORTHANTS"
                    : "FAIL")
            << '\n';
        return certified + infeasible == attempted
                   ? EXIT_SUCCESS
                   : EXIT_FAILURE;
    } catch (const z3::exception& error) {
        std::cerr << "z3 error: " << error.msg() << '\n';
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 3;
    }
}
