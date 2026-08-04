#include <boost/multiprecision/cpp_int.hpp>
#include <z3++.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

constexpr int kLevel = 5;
constexpr std::size_t kVariables = 6U;
constexpr std::size_t kQuadratics = 21U;
constexpr std::size_t kCubics = 56U;
constexpr std::size_t kQuartics = 126U;

using Linear = std::array<cpp_int, kVariables>;
using Quadratic = std::array<cpp_int, kQuadratics>;
using Cubic = std::array<cpp_int, kCubics>;
using Quartic = std::array<cpp_int, kQuartics>;
using LinearMatrix = std::array<std::array<Linear, kVariables>, kVariables>;

struct Constraint {
    std::string name;
    Quadratic polynomial{};
};

struct CubicConstraint {
    std::string name;
    Cubic polynomial{};
};

struct QuarticConstraint {
    std::string name;
    Quartic polynomial{};
};

std::size_t monomial(const int left, const int right) {
    const int low = std::min(left, right);
    const int high = std::max(left, right);
    std::size_t result = 0U;
    for (int row = 0; row < low; ++row) {
        result += kVariables - static_cast<std::size_t>(row);
    }
    return result + static_cast<std::size_t>(high - low);
}

std::size_t cubic_monomial(int first, int second, int third) {
    std::array<int, 3U> indices{first, second, third};
    std::sort(indices.begin(), indices.end());
    std::size_t result = 0U;
    for (int left = 0; left <= kLevel; ++left) {
        for (int middle = left; middle <= kLevel; ++middle) {
            for (int right = middle; right <= kLevel; ++right) {
                if (indices == std::array<int, 3U>{left, middle, right}) {
                    return result;
                }
                ++result;
            }
        }
    }
    throw std::logic_error("invalid cubic monomial");
}

std::size_t quartic_monomial(int first, int second, int third, int fourth) {
    std::array<int, 4U> indices{first, second, third, fourth};
    std::sort(indices.begin(), indices.end());
    std::size_t result = 0U;
    for (int left = 0; left <= kLevel; ++left) {
        for (int middle_left = left; middle_left <= kLevel; ++middle_left) {
            for (int middle_right = middle_left;
                 middle_right <= kLevel;
                 ++middle_right) {
                for (int right = middle_right; right <= kLevel; ++right) {
                    if (indices == std::array<int, 4U>{
                            left, middle_left, middle_right, right}) {
                        return result;
                    }
                    ++result;
                }
            }
        }
    }
    throw std::logic_error("invalid quartic monomial");
}

Linear unit(const int index) {
    Linear result{};
    result[static_cast<std::size_t>(index)] = 1;
    return result;
}

Linear add(const Linear& left, const Linear& right) {
    Linear result{};
    for (std::size_t index = 0U; index < kVariables; ++index) {
        result[index] = left[index] + right[index];
    }
    return result;
}

Linear subtract(const Linear& left, const Linear& right) {
    Linear result{};
    for (std::size_t index = 0U; index < kVariables; ++index) {
        result[index] = left[index] - right[index];
    }
    return result;
}

Quadratic product(const Linear& left, const Linear& right) {
    Quadratic result{};
    for (int row = 0; row <= kLevel; ++row) {
        for (int column = 0; column <= kLevel; ++column) {
            result[monomial(row, column)]
                += left[static_cast<std::size_t>(row)]
                   * right[static_cast<std::size_t>(column)];
        }
    }
    return result;
}

Quadratic subtract(const Quadratic& left, const Quadratic& right) {
    Quadratic result{};
    for (std::size_t index = 0U; index < kQuadratics; ++index) {
        result[index] = left[index] - right[index];
    }
    return result;
}

Cubic product(const Quadratic& quadratic, const Linear& linear) {
    Cubic result{};
    for (int left = 0; left <= kLevel; ++left) {
        for (int right = left; right <= kLevel; ++right) {
            const cpp_int coefficient = quadratic[monomial(left, right)];
            if (coefficient == 0) {
                continue;
            }
            for (int factor = 0; factor <= kLevel; ++factor) {
                result[cubic_monomial(left, right, factor)] += coefficient
                    * linear[static_cast<std::size_t>(factor)];
            }
        }
    }
    return result;
}

Cubic product(
    const Linear& first,
    const Linear& second,
    const Linear& third
) {
    Cubic result{};
    for (int first_index = 0; first_index <= kLevel; ++first_index) {
        for (int second_index = 0; second_index <= kLevel; ++second_index) {
            for (int third_index = 0; third_index <= kLevel; ++third_index) {
                result[cubic_monomial(first_index, second_index, third_index)]
                    += first[static_cast<std::size_t>(first_index)]
                     * second[static_cast<std::size_t>(second_index)]
                     * third[static_cast<std::size_t>(third_index)];
            }
        }
    }
    return result;
}

Quartic product(const Quadratic& first, const Quadratic& second) {
    Quartic result{};
    for (int first_left = 0; first_left <= kLevel; ++first_left) {
        for (int first_right = first_left;
             first_right <= kLevel;
             ++first_right) {
            const cpp_int first_coefficient = first[monomial(
                first_left, first_right
            )];
            if (first_coefficient == 0) {
                continue;
            }
            for (int second_left = 0; second_left <= kLevel; ++second_left) {
                for (int second_right = second_left;
                     second_right <= kLevel;
                     ++second_right) {
                    result[quartic_monomial(
                        first_left, first_right, second_left, second_right
                    )] += first_coefficient * second[monomial(
                        second_left, second_right
                    )];
                }
            }
        }
    }
    return result;
}

Quartic product(const Cubic& cubic, const Linear& linear) {
    Quartic result{};
    for (int first = 0; first <= kLevel; ++first) {
        for (int second = first; second <= kLevel; ++second) {
            for (int third = second; third <= kLevel; ++third) {
                const cpp_int coefficient = cubic[cubic_monomial(
                    first, second, third
                )];
                if (coefficient == 0) {
                    continue;
                }
                for (int fourth = 0; fourth <= kLevel; ++fourth) {
                    result[quartic_monomial(first, second, third, fourth)]
                        += coefficient
                         * linear[static_cast<std::size_t>(fourth)];
                }
            }
        }
    }
    return result;
}

void add_scaled(Cubic& total, const Cubic& value, const int scale) {
    for (std::size_t index = 0U; index < kCubics; ++index) {
        total[index] += scale * value[index];
    }
}

LinearMatrix multiplication_matrix() {
    LinearMatrix matrix{};
    for (int left = 0; left <= kLevel; ++left) {
        for (int right = 0; right <= kLevel; ++right) {
            const int lower = std::abs(left - right);
            const int upper = std::min(
                left + right,
                2 * kLevel - left - right
            );
            for (int label = lower; label <= upper; ++label) {
                matrix[static_cast<std::size_t>(left)]
                      [static_cast<std::size_t>(right)]
                      [static_cast<std::size_t>(label)] += 1;
            }
        }
    }
    return matrix;
}

Linear simple_current_block_entry(
    const LinearMatrix& matrix,
    const int sign,
    const int row,
    const int column
) {
    const Linear direct = matrix[static_cast<std::size_t>(row)]
                                [static_cast<std::size_t>(column)];
    const Linear reflected = matrix[static_cast<std::size_t>(row)]
                                   [static_cast<std::size_t>(kLevel - column)];
    return sign > 0 ? add(direct, reflected) : subtract(direct, reflected);
}

Linear bilinear_form(
    const LinearMatrix& matrix,
    const Linear& left,
    const Linear& right
) {
    Linear result{};
    for (int row = 0; row <= kLevel; ++row) {
        for (int column = 0; column <= kLevel; ++column) {
            const cpp_int weight
                = left[static_cast<std::size_t>(row)]
                  * right[static_cast<std::size_t>(column)];
            if (weight == 0) {
                continue;
            }
            for (std::size_t variable = 0U;
                 variable < kVariables;
                 ++variable) {
                result[variable]
                    += weight
                       * matrix[static_cast<std::size_t>(row)]
                               [static_cast<std::size_t>(column)][variable];
            }
        }
    }
    return result;
}

std::string render_vector(const Linear& vector) {
    std::string result;
    for (std::size_t index = 0U; index < kVariables; ++index) {
        if (index != 0U) {
            result += '_';
        }
        result += vector[index].convert_to<std::string>();
    }
    return result;
}

std::vector<Linear> sparse_vectors() {
    std::vector<Linear> result;
    for (int index = 0; index <= kLevel; ++index) {
        result.push_back(unit(index));
    }
    for (int left = 0; left <= kLevel; ++left) {
        for (int right = left + 1; right <= kLevel; ++right) {
            Linear plus = unit(left);
            plus[static_cast<std::size_t>(right)] = 1;
            result.push_back(plus);
            Linear minus = unit(left);
            minus[static_cast<std::size_t>(right)] = -1;
            result.push_back(minus);
        }
    }
    return result;
}

bool nonzero(const Quadratic& polynomial) {
    return std::any_of(
        polynomial.begin(), polynomial.end(),
        [](const cpp_int& coefficient) { return coefficient != 0; }
    );
}

std::vector<Constraint> constraints(const LinearMatrix& matrix) {
    std::vector<Constraint> result;
    for (int left = 0; left <= kLevel; ++left) {
        for (int right = left; right <= kLevel; ++right) {
            Quadratic polynomial{};
            polynomial[monomial(left, right)] = 1;
            result.push_back({
                "monomial_" + std::to_string(left)
                    + "_" + std::to_string(right),
                polynomial
            });
        }
    }
    for (int middle = 1; middle < kLevel; ++middle) {
        result.push_back({
            "log_concavity_" + std::to_string(middle),
            subtract(
                product(unit(middle), unit(middle)),
                product(unit(middle - 1), unit(middle + 1))
            )
        });
    }
    for (int radius = 1; radius <= kLevel; ++radius) {
        result.push_back({
            "boundary_" + std::to_string(radius),
            subtract(
                product(unit(0), unit(kLevel - radius)),
                product(unit(radius), unit(kLevel))
            )
        });
    }
    for (int left = 0; left <= kLevel; ++left) {
        for (int right = left + 1; right <= kLevel; ++right) {
            result.push_back({
                "psd_minor_" + std::to_string(left)
                    + "_" + std::to_string(right),
                subtract(
                    product(
                        matrix[static_cast<std::size_t>(left)]
                              [static_cast<std::size_t>(left)],
                        matrix[static_cast<std::size_t>(right)]
                              [static_cast<std::size_t>(right)]
                    ),
                    product(
                        matrix[static_cast<std::size_t>(left)]
                              [static_cast<std::size_t>(right)],
                        matrix[static_cast<std::size_t>(left)]
                              [static_cast<std::size_t>(right)]
                    )
                )
            });
        }
    }
    for (const std::pair<int, std::string>& block
         : std::array<std::pair<int, std::string>, 2U>{{
               {1, "even"}, {-1, "odd"}}}) {
        for (int left = 0; left < 3; ++left) {
            for (int right = left + 1; right < 3; ++right) {
                const Linear left_left = simple_current_block_entry(
                    matrix, block.first, left, left
                );
                const Linear right_right = simple_current_block_entry(
                    matrix, block.first, right, right
                );
                const Linear off_diagonal = simple_current_block_entry(
                    matrix, block.first, left, right
                );
                result.push_back({
                    "psd_" + block.second + "_minor_"
                        + std::to_string(left) + "_" + std::to_string(right),
                    subtract(
                        product(left_left, right_right),
                        product(off_diagonal, off_diagonal)
                    )
                });
            }
        }
    }
    const std::vector<Linear> sparse = sparse_vectors();
    for (std::size_t left = 0U; left < sparse.size(); ++left) {
        const Linear left_norm = bilinear_form(
            matrix, sparse[left], sparse[left]
        );
        for (int label = 0; label <= kLevel; ++label) {
            result.push_back({
                "psd_norm_" + render_vector(sparse[left])
                    + "_times_" + std::to_string(label),
                product(left_norm, unit(label))
            });
        }
        for (std::size_t right = left + 1U;
             right < sparse.size();
             ++right) {
            const Linear right_norm = bilinear_form(
                matrix, sparse[right], sparse[right]
            );
            const Linear inner = bilinear_form(
                matrix, sparse[left], sparse[right]
            );
            result.push_back({
                "psd_gram_" + render_vector(sparse[left])
                    + "_" + render_vector(sparse[right]),
                subtract(
                    product(left_norm, right_norm),
                    product(inner, inner)
                )
            });
        }
    }
    result.erase(
        std::remove_if(
            result.begin(), result.end(),
            [](const Constraint& constraint) {
                return !nonzero(constraint.polynomial);
            }
        ),
        result.end()
    );
    return result;
}

bool sparse_gram_constraint(const Constraint& constraint) {
    return constraint.name.rfind("psd_norm_", 0U) == 0U
        || constraint.name.rfind("psd_gram_", 0U) == 0U;
}

Cubic determinant(
    const std::array<std::array<Linear, 3U>, 3U>& block
) {
    Cubic result{};
    add_scaled(result, product(block[0U][0U], block[1U][1U], block[2U][2U]), 1);
    add_scaled(result, product(block[0U][1U], block[0U][2U], block[1U][2U]), 2);
    add_scaled(result, product(block[0U][0U], block[1U][2U], block[1U][2U]), -1);
    add_scaled(result, product(block[1U][1U], block[0U][2U], block[0U][2U]), -1);
    add_scaled(result, product(block[2U][2U], block[0U][1U], block[0U][1U]), -1);
    return result;
}

std::vector<CubicConstraint> cubic_constraints(
    const LinearMatrix& matrix,
    bool include_sparse_gram,
    bool include_block_determinants
) {
    std::vector<CubicConstraint> result;
    for (const Constraint& constraint : constraints(matrix)) {
        if (!include_sparse_gram && sparse_gram_constraint(constraint)) {
            continue;
        }
        for (int factor = 0; factor <= kLevel; ++factor) {
            result.push_back({
                "factor_" + std::to_string(factor)
                    + "_times_" + constraint.name,
                product(constraint.polynomial, unit(factor))
            });
        }
    }
    if (include_block_determinants) {
        for (const std::pair<int, std::string>& block
             : std::array<std::pair<int, std::string>, 2U>{
                   {{1, "even"}, {-1, "odd"}}}) {
            std::array<std::array<Linear, 3U>, 3U> entries{};
            for (int row = 0; row < 3; ++row) {
                for (int column = 0; column < 3; ++column) {
                    entries[static_cast<std::size_t>(row)]
                           [static_cast<std::size_t>(column)] =
                        simple_current_block_entry(
                            matrix, block.first, row, column
                        );
                }
            }
            result.push_back({
                "psd_" + block.second + "_determinant",
                determinant(entries)
            });
        }
    }
    return result;
}

std::vector<QuarticConstraint> quartic_constraints(
    const LinearMatrix& matrix
) {
    std::vector<QuarticConstraint> result;
    for (const Constraint& constraint : constraints(matrix)) {
        if (sparse_gram_constraint(constraint)) {
            continue;
        }
        for (int first = 0; first <= kLevel; ++first) {
            for (int second = first; second <= kLevel; ++second) {
                result.push_back({
                    "monomial_" + std::to_string(first)
                        + "_" + std::to_string(second)
                        + "_times_" + constraint.name,
                    product(
                        product(unit(first), unit(second)),
                        constraint.polynomial
                    )
                });
            }
        }
    }
    for (const std::pair<int, std::string>& block
         : std::array<std::pair<int, std::string>, 2U>{
               {{1, "even"}, {-1, "odd"}}}) {
        std::array<std::array<Linear, 3U>, 3U> entries{};
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column) {
                entries[static_cast<std::size_t>(row)]
                       [static_cast<std::size_t>(column)] =
                    simple_current_block_entry(
                        matrix, block.first, row, column
                    );
            }
        }
        const Cubic block_determinant = determinant(entries);
        for (int factor = 0; factor <= kLevel; ++factor) {
            result.push_back({
                "factor_" + std::to_string(factor)
                    + "_times_psd_" + block.second + "_determinant",
                product(block_determinant, unit(factor))
            });
        }
    }
    return result;
}

Quadratic current(
    const LinearMatrix& matrix,
    const int row,
    const int column
) {
    return subtract(
        product(
            unit(0),
            matrix[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(column)]
        ),
        product(unit(row), unit(column))
    );
}

std::pair<int, int> target(const std::string& name) {
    if (name.size() != 3U || name[0U] != 'a'
        || name[1U] < '0' || name[1U] > '5'
        || name[2U] < '0' || name[2U] > '5') {
        throw std::invalid_argument("target must have form aRS with 0<=R,S<=5");
    }
    return {name[1U] - '0', name[2U] - '0'};
}

z3::expr rational(z3::context& context, const cpp_int& value) {
    return context.real_val(value.convert_to<std::string>().c_str());
}

bool zero_numeral(const z3::expr& value) {
    if (!value.is_numeral()) {
        return false;
    }
    const std::string rendered = value.to_string();
    return rendered == "0" || rendered == "0.0";
}

cpp_int pairing(const Cubic& left, const Cubic& right) {
    cpp_int result = 0;
    for (std::size_t index = 0U; index < kCubics; ++index) {
        result += left[index] * right[index];
    }
    return result;
}

int cubic_farkas(
    const std::string& target_name,
    bool include_sparse_gram,
    bool include_block_determinants
) {
    const std::pair<int, int> coordinates = target(target_name);
    const LinearMatrix matrix = multiplication_matrix();
    const std::vector<CubicConstraint> cone = cubic_constraints(
        matrix, include_sparse_gram, include_block_determinants
    );
    const Cubic desired = product(
        current(matrix, coordinates.first, coordinates.second), unit(0)
    );

    z3::context context;
    z3::solver solver(context);
    std::vector<z3::expr> weights;
    weights.reserve(cone.size());
    for (std::size_t index = 0U; index < cone.size(); ++index) {
        weights.push_back(context.real_const(
            ("lambda3_" + std::to_string(index)).c_str()
        ));
        solver.add(weights.back() >= 0);
    }
    for (std::size_t monomial_index = 0U;
         monomial_index < kCubics;
         ++monomial_index) {
        z3::expr equation = context.real_val(0);
        for (std::size_t index = 0U; index < cone.size(); ++index) {
            equation = equation + weights[index] * rational(
                context, cone[index].polynomial[monomial_index]
            );
        }
        solver.add(equation == rational(context, desired[monomial_index]));
    }
    const z3::check_result result = solver.check();
    std::cout
        << "SU2_K5_CURRENT_CUBIC_FARKAS"
        << " target=" << target_name
        << " sparse_gram=" << (include_sparse_gram ? "on" : "off")
        << " block_psd=" << (include_block_determinants ? "on" : "off")
        << " constraints=" << cone.size()
        << " result=" << (result == z3::sat ? "SAT" : "UNSAT")
        << '\n';
    if (result != z3::sat) {
        return EXIT_SUCCESS;
    }
    const z3::model model = solver.get_model();
    for (std::size_t index = 0U; index < cone.size(); ++index) {
        const z3::expr value = model.eval(weights[index], true);
        if (zero_numeral(value)) {
            continue;
        }
        std::cout
            << "SU2_K5_CURRENT_CUBIC_FARKAS_TERM"
            << " weight=" << value
            << " constraint=" << cone[index].name
            << '\n';
    }
    return EXIT_SUCCESS;
}

int cubic_farkas_dual(
    const std::string& target_name,
    bool include_sparse_gram,
    bool include_block_determinants
) {
    const std::pair<int, int> coordinates = target(target_name);
    const LinearMatrix matrix = multiplication_matrix();
    const std::vector<CubicConstraint> cone = cubic_constraints(
        matrix, include_sparse_gram, include_block_determinants
    );
    const Cubic desired = product(
        current(matrix, coordinates.first, coordinates.second), unit(0)
    );

    z3::context context;
    z3::solver solver(context);
    std::vector<z3::expr> functional;
    functional.reserve(kCubics);
    for (std::size_t monomial_index = 0U;
         monomial_index < kCubics;
         ++monomial_index) {
        functional.push_back(context.real_const(
            ("dual_" + std::to_string(monomial_index)).c_str()
        ));
    }
    for (const CubicConstraint& constraint : cone) {
        z3::expr value = context.real_val(0);
        for (std::size_t monomial_index = 0U;
             monomial_index < kCubics;
             ++monomial_index) {
            value = value + functional[monomial_index] * rational(
                context, constraint.polynomial[monomial_index]
            );
        }
        solver.add(value >= 0);
    }
    z3::expr target_value = context.real_val(0);
    for (std::size_t monomial_index = 0U;
         monomial_index < kCubics;
         ++monomial_index) {
        target_value = target_value + functional[monomial_index] * rational(
            context, desired[monomial_index]
        );
    }
    solver.add(target_value <= context.real_val(-1));
    const z3::check_result result = solver.check();
    std::cout
        << "SU2_K5_CURRENT_CUBIC_FARKAS_DUAL"
        << " target=" << target_name
        << " sparse_gram=" << (include_sparse_gram ? "on" : "off")
        << " block_psd=" << (include_block_determinants ? "on" : "off")
        << " constraints=" << cone.size()
        << " result=" << (result == z3::sat ? "SAT" : "UNSAT")
        << '\n';
    if (result != z3::sat) {
        return EXIT_SUCCESS;
    }
    const z3::model model = solver.get_model();
    for (std::size_t monomial_index = 0U;
         monomial_index < kCubics;
         ++monomial_index) {
        const z3::expr value = model.eval(
            functional[monomial_index], true
        );
        if (zero_numeral(value)) {
            continue;
        }
        std::cout
            << "SU2_K5_CURRENT_CUBIC_FARKAS_DUAL_TERM"
            << " monomial=" << monomial_index
            << " weight=" << value
            << '\n';
    }
    return EXIT_SUCCESS;
}

int quartic_farkas(const std::string& target_name) {
    const std::pair<int, int> coordinates = target(target_name);
    const LinearMatrix matrix = multiplication_matrix();
    const std::vector<QuarticConstraint> cone = quartic_constraints(matrix);
    const Quadratic a_squared = product(unit(0), unit(0));
    const Quartic desired = product(
        a_squared, current(matrix, coordinates.first, coordinates.second)
    );

    z3::context context;
    z3::solver solver(context);
    std::vector<z3::expr> weights;
    weights.reserve(cone.size());
    for (std::size_t index = 0U; index < cone.size(); ++index) {
        weights.push_back(context.real_const(
            ("lambda4_" + std::to_string(index)).c_str()
        ));
        solver.add(weights.back() >= 0);
    }
    for (std::size_t monomial_index = 0U;
         monomial_index < kQuartics;
         ++monomial_index) {
        z3::expr equation = context.real_val(0);
        for (std::size_t index = 0U; index < cone.size(); ++index) {
            equation = equation + weights[index] * rational(
                context, cone[index].polynomial[monomial_index]
            );
        }
        solver.add(equation == rational(context, desired[monomial_index]));
    }
    const z3::check_result result = solver.check();
    std::cout
        << "SU2_K5_CURRENT_QUARTIC_FARKAS"
        << " target=" << target_name
        << " constraints=" << cone.size()
        << " result=" << (result == z3::sat ? "SAT" : "UNSAT")
        << '\n';
    if (result != z3::sat) {
        return EXIT_SUCCESS;
    }
    const z3::model model = solver.get_model();
    for (std::size_t index = 0U; index < cone.size(); ++index) {
        const z3::expr value = model.eval(weights[index], true);
        if (zero_numeral(value)) {
            continue;
        }
        std::cout
            << "SU2_K5_CURRENT_QUARTIC_FARKAS_TERM"
            << " weight=" << value
            << " constraint=" << cone[index].name
            << '\n';
    }
    return EXIT_SUCCESS;
}

int verify_a12_cubic_dual() {
    const LinearMatrix matrix = multiplication_matrix();
    const std::vector<CubicConstraint> cone = cubic_constraints(
        matrix, false, true
    );
    Cubic functional{};
    functional[cubic_monomial(0, 0, 0)] = 1;
    functional[cubic_monomial(0, 0, 5)] = 1;
    functional[cubic_monomial(0, 1, 2)] = 4;
    functional[cubic_monomial(0, 5, 5)] = 1;
    functional[cubic_monomial(1, 1, 1)] = 4;
    functional[cubic_monomial(5, 5, 5)] = 1;

    const Cubic desired = product(current(matrix, 1, 2), unit(0));
    const cpp_int target_value = pairing(functional, desired);
    std::size_t negative_generators = 0U;
    for (const CubicConstraint& constraint : cone) {
        if (pairing(functional, constraint.polynomial) < 0) {
            ++negative_generators;
        }
    }
    std::cout
        << "SU2_K5_CURRENT_CUBIC_FARKAS_DUAL_VERIFY"
        << " target=a12"
        << " constraints=" << cone.size()
        << " target_value=" << target_value
        << " negative_generators=" << negative_generators
        << " result="
        << (target_value < 0 && negative_generators == 0U
            ? "PASS_EXACT_DUAL" : "FAIL")
        << '\n';
    return target_value < 0 && negative_generators == 0U
        ? EXIT_SUCCESS : EXIT_FAILURE;
}

std::vector<std::vector<z3::expr>> z3_multiplication_matrix(
    z3::context& context,
    const std::vector<z3::expr>& profile
) {
    std::vector<std::vector<z3::expr>> matrix(
        kVariables,
        std::vector<z3::expr>(kVariables, context.real_val(0))
    );
    for (int left = 0; left <= kLevel; ++left) {
        for (int right = 0; right <= kLevel; ++right) {
            const int lower = std::abs(left - right);
            const int upper = std::min(
                left + right,
                2 * kLevel - left - right
            );
            for (int label = lower; label <= upper; ++label) {
                matrix[static_cast<std::size_t>(left)]
                      [static_cast<std::size_t>(right)]
                    = matrix[static_cast<std::size_t>(left)]
                            [static_cast<std::size_t>(right)]
                      + profile[static_cast<std::size_t>(label)];
            }
        }
    }
    return matrix;
}

void add_psd_three_by_three(
    z3::solver& solver,
    const std::array<std::array<z3::expr, 3U>, 3U>& block
) {
    for (std::size_t row = 0U; row < 3U; ++row) {
        solver.add(block[row][row] >= 0);
    }
    for (std::size_t left = 0U; left < 3U; ++left) {
        for (std::size_t right = left + 1U; right < 3U; ++right) {
            solver.add(
                block[left][left] * block[right][right]
                - block[left][right] * block[left][right] >= 0
            );
        }
    }
    const z3::expr determinant
        = block[0U][0U] * block[1U][1U] * block[2U][2U]
          + 2 * block[0U][1U] * block[0U][2U] * block[1U][2U]
          - block[0U][0U] * block[1U][2U] * block[1U][2U]
          - block[1U][1U] * block[0U][2U] * block[0U][2U]
          - block[2U][2U] * block[0U][1U] * block[0U][1U];
    solver.add(determinant >= 0);
}

int nonlinear_counterexample(const std::string& target_name) {
    const std::pair<int, int> coordinates = target(target_name);
    z3::context context;
    z3::solver solver(context);
    std::vector<z3::expr> profile;
    profile.reserve(kVariables);
    for (int label = 0; label <= kLevel; ++label) {
        profile.push_back(context.real_const(
            ("d_" + std::to_string(label)).c_str()
        ));
        solver.add(profile.back() >= 0);
    }
    // A negative anchored current forces d_0>0: otherwise positive
    // semidefiniteness makes the zeroth row vanish.  The constraints and
    // the current are homogeneous, so this removes the irrelevant scale.
    solver.add(profile[0U] == 1);
    for (int label = 1; label < kLevel; ++label) {
        solver.add(
            profile[static_cast<std::size_t>(label)]
                * profile[static_cast<std::size_t>(label)]
            >= profile[static_cast<std::size_t>(label - 1)]
                * profile[static_cast<std::size_t>(label + 1)]
        );
    }
    for (int radius = 1; radius <= kLevel; ++radius) {
        solver.add(
            profile[0U] * profile[static_cast<std::size_t>(kLevel - radius)]
            >= profile[static_cast<std::size_t>(radius)] * profile.back()
        );
    }
    const std::vector<std::vector<z3::expr>> matrix
        = z3_multiplication_matrix(context, profile);
    for (const int sign : std::array<int, 2U>{1, -1}) {
        std::array<std::array<z3::expr, 3U>, 3U> block{
            std::array<z3::expr, 3U>{
                context.real_val(0), context.real_val(0), context.real_val(0)},
            std::array<z3::expr, 3U>{
                context.real_val(0), context.real_val(0), context.real_val(0)},
            std::array<z3::expr, 3U>{
                context.real_val(0), context.real_val(0), context.real_val(0)}};
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column) {
                const z3::expr direct = matrix[static_cast<std::size_t>(row)]
                                                  [static_cast<std::size_t>(column)];
                const z3::expr reflected
                    = matrix[static_cast<std::size_t>(row)]
                            [static_cast<std::size_t>(kLevel - column)];
                block[static_cast<std::size_t>(row)]
                     [static_cast<std::size_t>(column)]
                    = sign > 0 ? direct + reflected : direct - reflected;
            }
        }
        add_psd_three_by_three(solver, block);
    }
    solver.add(
        profile[0U]
            * matrix[static_cast<std::size_t>(coordinates.first)]
                    [static_cast<std::size_t>(coordinates.second)]
        < profile[static_cast<std::size_t>(coordinates.first)]
            * profile[static_cast<std::size_t>(coordinates.second)]
    );
    const z3::check_result result = solver.check();
    std::cout
        << "SU2_K5_CURRENT_COUNTEREXAMPLE"
        << " target=" << target_name
        << " result="
        << (result == z3::sat
                ? "SAT"
                : result == z3::unsat ? "UNSAT" : "UNKNOWN")
        << '\n';
    if (result == z3::sat) {
        const z3::model model = solver.get_model();
        std::cout << "SU2_K5_CURRENT_COUNTEREXAMPLE_PROFILE";
        for (const z3::expr& coordinate : profile) {
            std::cout << ' ' << model.eval(coordinate, true);
        }
        std::cout << '\n';
    }
    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--counterexample") {
            return nonlinear_counterexample(argv[2]);
        }
        if (argc == 3 && std::string(argv[1]) == "--degree3") {
            return cubic_farkas(argv[2], true, true);
        }
        if (argc == 3 && std::string(argv[1]) == "--degree3-basic") {
            return cubic_farkas(argv[2], false, false);
        }
        if (argc == 3 && std::string(argv[1]) == "--degree3-block-psd") {
            return cubic_farkas(argv[2], false, true);
        }
        if (argc == 3 && std::string(argv[1])
            == "--degree3-block-psd-dual") {
            return cubic_farkas_dual(argv[2], false, true);
        }
        if (argc == 3 && std::string(argv[1])
            == "--degree4-block-psd") {
            return quartic_farkas(argv[2]);
        }
        if (argc == 2 && std::string(argv[1])
            == "--verify-degree3-block-psd-a12-dual") {
            return verify_a12_cubic_dual();
        }
        const std::string target_name = argc == 2 ? argv[1] : "a12";
        if (argc > 2) {
            throw std::invalid_argument(
                "usage: search_su2_k5_current_farkas [aRS] | "
                "--counterexample aRS | --degree3 aRS | "
                "--degree3-basic aRS | --degree3-block-psd aRS | "
                "--degree3-block-psd-dual aRS | "
                "--degree4-block-psd aRS | "
                "--verify-degree3-block-psd-a12-dual"
            );
        }
        const std::pair<int, int> coordinates = target(target_name);
        const LinearMatrix matrix = multiplication_matrix();
        const std::vector<Constraint> cone = constraints(matrix);
        const Quadratic desired = current(
            matrix, coordinates.first, coordinates.second
        );

        z3::context context;
        z3::solver solver(context);
        std::vector<z3::expr> weights;
        weights.reserve(cone.size());
        for (std::size_t index = 0U; index < cone.size(); ++index) {
            weights.push_back(context.real_const(
                ("lambda_" + std::to_string(index)).c_str()
            ));
            solver.add(weights.back() >= 0);
        }
        for (std::size_t monomial_index = 0U;
             monomial_index < kQuadratics;
             ++monomial_index) {
            z3::expr equation = context.real_val(0);
            for (std::size_t index = 0U; index < cone.size(); ++index) {
                equation = equation
                    + weights[index] * rational(
                        context, cone[index].polynomial[monomial_index]
                    );
            }
            solver.add(equation == rational(
                context, desired[monomial_index]
            ));
        }
        const z3::check_result result = solver.check();
        std::cout
            << "SU2_K5_CURRENT_FARKAS"
            << " target=" << target_name
            << " constraints=" << cone.size()
            << " result=" << (result == z3::sat ? "SAT" : "UNSAT")
            << '\n';
        if (result != z3::sat) {
            return EXIT_SUCCESS;
        }
        const z3::model model = solver.get_model();
        for (std::size_t index = 0U; index < cone.size(); ++index) {
            const z3::expr value = model.eval(weights[index], true);
            if (value.is_numeral() && value.to_string() == "0") {
                continue;
            }
            std::cout
                << "SU2_K5_CURRENT_FARKAS_TERM"
                << " weight=" << value
                << " constraint=" << cone[index].name
                << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
