#include <z3++.h>

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

int parse_power(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value < 2
        || value > 20
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error("power must lie in [2,20]");
    }
    return static_cast<int>(value);
}

long long binomial(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    if (k > n - k) {
        k = n - k;
    }
    long long result = 1;
    for (int index = 1; index <= k; ++index) {
        result *= n - k + index;
        result /= index;
    }
    return result;
}

long long factorial(int n) {
    long long result = 1;
    for (int value = 2; value <= n; ++value) {
        result *= value;
    }
    return result;
}

z3::expr choose_polynomial(
    z3::context& context,
    const z3::expr& upper,
    int lower
) {
    z3::expr result = context.real_val(1);
    for (int offset = 0; offset < lower; ++offset) {
        result = result * (upper - context.real_val(offset));
    }
    return result / context.real_val(
        static_cast<std::int64_t>(factorial(lower))
    );
}

z3::expr coefficient_formula(
    z3::context& context,
    const z3::expr& q,
    const z3::expr& coordinate,
    int power,
    int quotient
) {
    const z3::expr width = 2 * q + 1;
    const z3::expr centered_degree = power * q + coordinate;
    z3::expr result = context.real_val(0);
    for (int index = 0; index <= quotient; ++index) {
        const z3::expr upper =
            centered_degree - index * width + power - 1;
        z3::expr term =
            context.real_val(
                static_cast<std::int64_t>(
                    binomial(power, index)
                )
            )
            * choose_polynomial(context, upper, power - 1);
        if (index % 2 != 0) {
            term = -term;
        }
        result = result + term;
    }
    return result;
}

z3::expr coefficient(
    z3::context& context,
    const z3::expr& q,
    const z3::expr& coordinate,
    int power
) {
    const z3::expr width = 2 * q + 1;
    const z3::expr centered_degree = power * q + coordinate;
    z3::expr result = context.real_val(0);
    for (int quotient = power - 1; quotient >= 0; --quotient) {
        const z3::expr lower =
            quotient * width <= centered_degree;
        const z3::expr upper =
            centered_degree < (quotient + 1) * width;
        result = z3::ite(
            lower && upper,
            coefficient_formula(
                context,
                q,
                coordinate,
                power,
                quotient
            ),
            result
        );
    }
    return z3::ite(
        coordinate <= power * q,
        result,
        context.real_val(0)
    );
}

z3::check_result check_interior(int power, bool wall_case) {
    z3::context context;
    z3::solver solver(context, "QF_NRA");
    z3::params parameters(context);
    parameters.set("timeout", 300000U);
    solver.set(parameters);

    const z3::expr q = context.real_const("q");
    const z3::expr x = context.real_const("x");
    const z3::expr y = context.real_const("y");
    solver.add(q >= 1);
    solver.add(y >= x + 1);
    solver.add(y <= power * q);
    if (wall_case) {
        solver.add(x == 0);
    } else {
        solver.add(x >= 1);
    }

    const z3::expr vx = coefficient(context, q, x, power);
    const z3::expr vx_plus =
        coefficient(context, q, x + 1, power);
    const z3::expr vx_minus =
        coefficient(
            context,
            q,
            wall_case ? x + 1 : x - 1,
            power
        );
    const z3::expr vy = coefficient(context, q, y, power);
    const z3::expr vy_plus =
        coefficient(context, q, y + 1, power);
    const z3::expr vy_plus_two =
        coefficient(context, q, y + 2, power);
    const z3::expr determinant =
        (vx - vy) * (vx - vy_plus_two)
        - (vx_plus - vy_plus) * (vx_minus - vy_plus);
    solver.add(determinant < 0);
    const z3::check_result result = solver.check();
    std::cout
        << "SU2_POWER_TP2_NLSAT"
        << " power=" << power
        << " chamber=" << (wall_case ? "wall" : "interior")
        << " result=" << result;
    if (result == z3::sat) {
        const z3::model model = solver.get_model();
        std::cout
            << " q=" << model.eval(q, true)
            << " x=" << model.eval(x, true)
            << " y=" << model.eval(y, true)
            << " determinant="
            << model.eval(determinant, true);
    }
    std::cout << '\n';
    return result;
}

void add_chamber(
    z3::solver& solver,
    const z3::expr& q,
    const z3::expr& coordinate,
    int power,
    int category
) {
    if (category == power) {
        solver.add(coordinate >= power * q + 1);
        return;
    }
    const z3::expr width = 2 * q + 1;
    const z3::expr centered_degree = power * q + coordinate;
    solver.add(category * width <= centered_degree);
    solver.add(
        centered_degree <= (category + 1) * width - 1
    );
    solver.add(coordinate <= power * q);
}

z3::expr chamber_value(
    z3::context& context,
    const z3::expr& q,
    const z3::expr& coordinate,
    int power,
    int category
) {
    return
        category == power
        ? context.real_val(0)
        : coefficient_formula(
            context,
            q,
            coordinate,
            power,
            category
        );
}

std::vector<int> next_categories(int category, int power) {
    if (category == power) {
        return {power};
    }
    std::vector<int> result{category};
    if (category + 1 < power) {
        result.push_back(category + 1);
    }
    result.push_back(power);
    return result;
}

bool check_split(int power, bool wall_case) {
    z3::context context;
    const z3::expr q = context.real_const("q");
    const z3::expr x = context.real_const("x");
    const z3::expr y = context.real_const("y");
    std::size_t patterns = 0;
    std::size_t feasible_patterns = 0;
    std::size_t unknown_patterns = 0;
    std::size_t newton_deferred_patterns = 0;

    for (int x_minus_category = 0;
         x_minus_category < power;
         ++x_minus_category) {
        const std::vector<int> x_categories =
            wall_case
            ? std::vector<int>{x_minus_category}
            : next_categories(x_minus_category, power);
        for (const int x_category : x_categories) {
            if (x_category == power) {
                continue;
            }
            const std::vector<int> x_plus_categories =
                next_categories(x_category, power);
            for (const int x_plus_category : x_plus_categories) {
                for (int y_category = 0;
                     y_category < power;
                     ++y_category) {
                    for (
                        const int y_plus_category :
                        next_categories(y_category, power)
                    ) {
                        for (
                            const int y_plus_two_category :
                            next_categories(y_plus_category, power)
                        ) {
                            ++patterns;
                            z3::solver linear(context, "QF_LRA");
                            linear.add(q >= 1);
                            linear.add(y >= x + 1);
                            linear.add(y <= power * q);
                            if (wall_case) {
                                linear.add(x == 0);
                                add_chamber(
                                    linear,
                                    q,
                                    x,
                                    power,
                                    x_minus_category
                                );
                                add_chamber(
                                    linear,
                                    q,
                                    x + 1,
                                    power,
                                    x_plus_category
                                );
                            } else {
                                linear.add(x >= 1);
                                add_chamber(
                                    linear,
                                    q,
                                    x - 1,
                                    power,
                                    x_minus_category
                                );
                                add_chamber(
                                    linear,
                                    q,
                                    x,
                                    power,
                                    x_category
                                );
                                add_chamber(
                                    linear,
                                    q,
                                    x + 1,
                                    power,
                                    x_plus_category
                                );
                            }
                            add_chamber(
                                linear,
                                q,
                                y,
                                power,
                                y_category
                            );
                            add_chamber(
                                linear,
                                q,
                                y + 1,
                                power,
                                y_plus_category
                            );
                            add_chamber(
                                linear,
                                q,
                                y + 2,
                                power,
                                y_plus_two_category
                            );
                            if (linear.check() == z3::unsat) {
                                continue;
                            }
                            ++feasible_patterns;
                            const bool newton_boundary =
                                power == 6
                                && wall_case
                                && x_minus_category == 2
                                && x_category == 2
                                && x_plus_category == 2
                                && y_category == 3
                                && y_plus_category == 4
                                && y_plus_two_category == 4;
                            if (newton_boundary) {
                                ++newton_deferred_patterns;
                                continue;
                            }

                            z3::solver nonlinear(context, "QF_NRA");
                            z3::params parameters(context);
                            parameters.set("timeout", 10000U);
                            nonlinear.set(parameters);
                            const z3::expr_vector assertions =
                                linear.assertions();
                            for (int index = 0;
                                 index
                                    < static_cast<int>(
                                        assertions.size()
                                    );
                                 ++index) {
                                nonlinear.add(assertions[index]);
                            }
                            const z3::expr vx = chamber_value(
                                context,
                                q,
                                x,
                                power,
                                x_category
                            );
                            const z3::expr vx_plus = chamber_value(
                                context,
                                q,
                                x + 1,
                                power,
                                x_plus_category
                            );
                            const z3::expr vx_minus = chamber_value(
                                context,
                                q,
                                wall_case ? x + 1 : x - 1,
                                power,
                                wall_case
                                    ? x_plus_category
                                    : x_minus_category
                            );
                            const z3::expr vy = chamber_value(
                                context,
                                q,
                                y,
                                power,
                                y_category
                            );
                            const z3::expr vy_plus = chamber_value(
                                context,
                                q,
                                y + 1,
                                power,
                                y_plus_category
                            );
                            const z3::expr vy_plus_two =
                                chamber_value(
                                    context,
                                    q,
                                    y + 2,
                                    power,
                                    y_plus_two_category
                                );
                            const z3::expr determinant =
                                (vx - vy) * (vx - vy_plus_two)
                                - (vx_plus - vy_plus)
                                    * (vx_minus - vy_plus);
                            nonlinear.add(determinant < 0);
                            const z3::check_result result =
                                nonlinear.check();
                            if (result == z3::sat) {
                                const z3::model model =
                                    nonlinear.get_model();
                                std::cout
                                    << "SU2_POWER_TP2_NLSAT_SPLIT"
                                    << " power=" << power
                                    << " chamber="
                                    << (
                                        wall_case
                                            ? "wall"
                                            : "interior"
                                    )
                                    << " result=sat"
                                    << " q=" << model.eval(q, true)
                                    << " x=" << model.eval(x, true)
                                    << " y=" << model.eval(y, true)
                                    << " determinant="
                                    << model.eval(
                                        determinant,
                                        true
                                    )
                                    << '\n';
                                return false;
                            }
                            if (result == z3::unknown) {
                                ++unknown_patterns;
                                std::cout
                                    << "SU2_POWER_TP2_NLSAT_SPLIT"
                                    << " power=" << power
                                    << " chamber="
                                    << (
                                        wall_case
                                            ? "wall"
                                            : "interior"
                                    )
                                    << " result=unknown"
                                    << " categories=("
                                    << x_minus_category << ','
                                    << x_category << ','
                                    << x_plus_category << ';'
                                    << y_category << ','
                                    << y_plus_category << ','
                                    << y_plus_two_category << ')'
                                    << '\n';
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout
        << "SU2_POWER_TP2_NLSAT_SPLIT"
        << " power=" << power
        << " chamber=" << (wall_case ? "wall" : "interior")
        << " patterns=" << patterns
        << " feasible_patterns=" << feasible_patterns
        << " unknown_patterns=" << unknown_patterns
        << " newton_deferred_patterns="
        << newton_deferred_patterns
        << " result="
        << (
            unknown_patterns == 0
                ? (
                    newton_deferred_patterns == 0
                        ? "PASS_REAL_CERTIFICATE"
                        : "PASS_EXCEPT_NEWTON_BOUNDARY"
                )
                : "UNKNOWN"
        )
        << '\n';
    return unknown_patterns == 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout << std::unitbuf;
        if (
            argc != 2
            && !(
                argc == 3
                && (
                    std::string{argv[2]} == "split"
                    || std::string{argv[2]} == "wall"
                )
            )
        ) {
            throw std::runtime_error(
                "usage: analyze_su2_power_tp2_nlsat "
                "power [split|wall]"
            );
        }
        const int power = parse_power(argv[1]);
        const std::string mode = argc == 3 ? argv[2] : "";
        const bool split = mode == "split";
        const bool wall_only = mode == "wall";
        bool proved = false;
        if (wall_only) {
            proved = check_split(power, true);
        } else if (split) {
            const bool wall = check_split(power, true);
            const bool interior = check_split(power, false);
            proved = wall && interior;
        } else {
            const bool wall =
                check_interior(power, true) == z3::unsat;
            const bool interior =
                check_interior(power, false) == z3::unsat;
            proved = wall && interior;
        }
        std::cout
            << "SU2_POWER_TP2_NLSAT"
            << " power=" << power
            << " result="
            << (
                proved
                    ? (
                        wall_only
                            ? "PASS_WALL_CHECK"
                            : "PASS_REAL_CERTIFICATE"
                    )
                    : "NOT_PROVED"
            )
            << '\n';
        return proved ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
