#include <z3++.h>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#define main analyze_su2_aim_square_root_diagonals_embedded_main
#include "analyze_su2_aim_square_root_diagonals.cpp"
#undef main

namespace {

struct LpFailure {
    int level = -1;
    int shell = -1;
    int factor = -1;
    int radius = -1;
};

std::string render_lp_failure(const LpFailure& failure) {
    return "level=" + std::to_string(
        failure.level < 0 ? -1 : 2 * failure.level)
        + " shell=" + std::to_string(failure.shell)
        + " factor=" + std::to_string(failure.factor)
        + " radius=" + std::to_string(failure.radius);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 10;
        const int maximum_shell = argc >= 3
            ? positive_argument(argv[2], "maximum_shell")
            : 2;
        const int maximum_periods = argc >= 4
            ? positive_argument(argv[3], "maximum_periods")
            : 12;
        if (argc > 4 || maximum_level < 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim_square_root_diagonal_lp "
                "[maximum_half_level] [maximum_shell] "
                "[maximum_periods]");
        }

        std::uint64_t payments = 0U;
        std::uint64_t constraints = 0U;
        std::uint64_t satisfiable = 0U;
        std::uint64_t unsatisfiable = 0U;
        LpFailure first_unsatisfiable;
        std::string first_model;

        for (int level = 2; level <= maximum_level; ++level) {
            const int period = 2 * level + 2;
            const int maximum_total = maximum_periods * period;
            for (int factor = 1;
                 factor <= level / 2;
                 ++factor) {
                const Matrix transform
                    = reserve_transform(level, factor);
                for (int radius = 0; radius <= level; ++radius) {
                    for (int shell = 1;
                         shell <= maximum_shell;
                         ++shell) {
                        ++payments;
                        z3::context context;
                        z3::solver solver(context);
                        z3::expr_vector lambda(context);
                        for (int source = 0;
                             source <= level;
                             ++source) {
                            lambda.push_back(
                                context.real_const(
                                    ("lambda_"
                                     + std::to_string(source))
                                        .c_str()));
                            solver.add(
                                lambda[source] >= 0);
                        }

                        Vector payment_prefix(
                            static_cast<std::size_t>(
                                maximum_total + 1),
                            0);
                        std::vector<Vector> reserve_prefix(
                            static_cast<std::size_t>(level + 1),
                            Vector(
                                static_cast<std::size_t>(
                                    maximum_total + 1),
                                0));
                        long long payment_sum = 0;
                        Vector reserve_sum(
                            static_cast<std::size_t>(level + 1),
                            0);
                        for (int index = 0;
                             index <= maximum_total;
                             ++index) {
                            payment_sum += payment_coefficient(
                                index,
                                level,
                                shell,
                                factor,
                                radius,
                                transform[
                                    static_cast<std::size_t>(
                                        radius)]);
                            payment_prefix[
                                static_cast<std::size_t>(index)]
                                = payment_sum;
                            for (int source = 0;
                                 source <= level;
                                 ++source) {
                                reserve_sum[
                                    static_cast<std::size_t>(source)]
                                    += reserve_coefficient(
                                        index,
                                        level,
                                        shell,
                                        source);
                                reserve_prefix[
                                    static_cast<std::size_t>(source)]
                                    [static_cast<std::size_t>(index)]
                                    = reserve_sum[
                                        static_cast<std::size_t>(
                                            source)];
                            }
                        }

                        for (int total = 0;
                             total <= maximum_total;
                             ++total) {
                            const int middle = total / 2;
                            long long payment_suffix = 0;
                            Vector reserve_suffix(
                                static_cast<std::size_t>(level + 1),
                                0);
                            for (int left = middle;
                                 left >= 0;
                                 --left) {
                                const int right = total - left;
                                const int gap = right - left;
                                const int multiplicity
                                    = left == right ? 1 : 2;
                                payment_suffix += multiplicity
                                    * (payment_prefix[
                                           static_cast<std::size_t>(
                                               total)]
                                       - (gap == 0
                                              ? 0
                                              : payment_prefix[
                                                    static_cast<
                                                        std::size_t>(
                                                        gap - 1)]));
                                for (int source = 0;
                                     source <= level;
                                     ++source) {
                                    reserve_suffix[
                                        static_cast<std::size_t>(
                                            source)]
                                        += multiplicity
                                           * (reserve_prefix[
                                                  static_cast<
                                                      std::size_t>(
                                                      source)]
                                                  [static_cast<
                                                      std::size_t>(
                                                      total)]
                                              - (gap == 0
                                                     ? 0
                                                     : reserve_prefix[
                                                           static_cast<
                                                               std::size_t>(
                                                               source)]
                                                           [static_cast<
                                                               std::size_t>(
                                                               gap - 1)]));
                                }
                                z3::expr inequality
                                    = context.int_val(
                                        static_cast<std::int64_t>(
                                            payment_suffix));
                                for (int source = 0;
                                     source <= level;
                                     ++source) {
                                    inequality
                                        = inequality
                                          - lambda[source]
                                            * context.int_val(
                                                static_cast<
                                                    std::int64_t>(
                                                    reserve_suffix[
                                                        static_cast<
                                                            std::size_t>(
                                                            source)]));
                                }
                                solver.add(inequality >= 0);
                                ++constraints;
                            }
                        }

                        if (solver.check() == z3::sat) {
                            ++satisfiable;
                            if (first_model.empty()) {
                                const z3::model model
                                    = solver.get_model();
                                first_model = "level="
                                    + std::to_string(2 * level)
                                    + " shell="
                                    + std::to_string(shell)
                                    + " factor="
                                    + std::to_string(factor)
                                    + " radius="
                                    + std::to_string(radius)
                                    + " lambda=[";
                                for (int source = 0;
                                     source <= level;
                                     ++source) {
                                    if (source != 0) {
                                        first_model += ',';
                                    }
                                    first_model += model.eval(
                                        lambda[source],
                                        true).to_string();
                                }
                                first_model += ']';
                            }
                        } else {
                            ++unsatisfiable;
                            if (first_unsatisfiable.level < 0) {
                                first_unsatisfiable = {
                                    level,
                                    shell,
                                    factor,
                                    radius};
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_AIM_SQUARE_ROOT_DIAGONAL_LP"
            << " maximum_level=" << 2 * maximum_level
            << " maximum_shell=" << maximum_shell
            << " maximum_periods=" << maximum_periods
            << " payments=" << payments
            << " constraints=" << constraints
            << " satisfiable=" << satisfiable
            << " unsatisfiable=" << unsatisfiable
            << '\n'
            << "FIRST_UNSATISFIABLE "
            << render_lp_failure(first_unsatisfiable)
            << '\n'
            << "FIRST_MODEL " << first_model
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
