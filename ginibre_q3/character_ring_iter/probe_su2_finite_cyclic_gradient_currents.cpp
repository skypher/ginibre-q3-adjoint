#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Vector = std::vector<cpp_int>;

struct Witness {
    int level = -1;
    int factor = -1;
    int power = -1;
    int target = -1;
    int radius = -1;
    cpp_int value = 0;
};

int parse_nonnegative(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed < 0LL
        || parsed > static_cast<long long>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument(std::string(name) + " must be a nonnegative integer");
    }
    return static_cast<int>(parsed);
}

int residue(int value, int modulus) {
    const int reduced = value % modulus;
    return reduced < 0 ? reduced + modulus : reduced;
}

Vector finite_step(const Vector& current, int level, int factor) {
    const std::size_t dimension = static_cast<std::size_t>(level + 1);
    Vector next(dimension, 0);
    for (int source = 0; source <= level; ++source) {
        const cpp_int& coefficient = current[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        const int lower = std::abs(source - factor);
        const int upper = std::min(source + factor, 2 * level - source - factor);
        for (int target = lower; target <= upper; ++target) {
            next[static_cast<std::size_t>(target)] += coefficient;
        }
    }
    return next;
}

Vector cyclic_interval_step(const Vector& current, int factor) {
    const int order = static_cast<int>(current.size());
    Vector next(current.size(), 0);
    for (int position = 0; position < order; ++position) {
        const cpp_int& coefficient = current[static_cast<std::size_t>(position)];
        if (coefficient == 0) {
            continue;
        }
        for (int shift = -factor; shift <= factor; ++shift) {
            const int target = residue(position + shift, order);
            next[static_cast<std::size_t>(target)] += coefficient;
        }
    }
    return next;
}

cpp_int cyclic_gradient_at(const Vector& cyclic, int position) {
    const int order = static_cast<int>(cyclic.size());
    const int left = residue(position, order);
    const int right = residue(position + 1, order);
    return cyclic[static_cast<std::size_t>(left)]
        - cyclic[static_cast<std::size_t>(right)];
}

Vector finite_gradient(const Vector& cyclic, int level) {
    Vector result(static_cast<std::size_t>(level + 1), 0);
    for (int position = 0; position <= level; ++position) {
        result[static_cast<std::size_t>(position)]
            = cyclic_gradient_at(cyclic, position);
    }
    return result;
}

void print_witness(const char* label, const Witness& witness) {
    std::cout << ' ' << label
              << "_level=" << 2 * witness.level
              << ' ' << label << "_factor=" << 2 * witness.factor
              << ' ' << label << "_half_power=" << witness.power
              << ' ' << label << "_target=" << 2 * witness.target
              << ' ' << label << "_radius=" << 2 * witness.radius
              << ' ' << label << "_value=" << witness.value;
}

void print_shell_and_prefix_sequences(
    const char* label,
    const Witness& witness
) {
    Vector finite(static_cast<std::size_t>(witness.level + 1), 0);
    finite[0] = 1;
    for (int step = 0; step < 2 * witness.power; ++step) {
        finite = finite_step(finite, witness.level, witness.factor);
    }

    cpp_int previous = 0;
    std::cout << ' ' << label << "_shells=[";
    Vector prefixes;
    prefixes.reserve(static_cast<std::size_t>(witness.level));
    for (int radius = 1; radius <= witness.level; ++radius) {
        const Vector image = finite_step(finite, witness.level, radius);
        const cpp_int current
            = finite[0] * image[static_cast<std::size_t>(witness.target)]
            - finite[static_cast<std::size_t>(radius)]
                * finite[static_cast<std::size_t>(witness.target)];
        if (radius > 1) {
            std::cout << ',';
        }
        std::cout << current - previous;
        prefixes.push_back(current);
        previous = current;
    }
    std::cout << "] " << label << "_prefixes=[";
    for (std::size_t index = 0U; index < prefixes.size(); ++index) {
        if (index > 0U) {
            std::cout << ',';
        }
        std::cout << prefixes[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? parse_nonnegative(argv[1], "maximum_level")
            : 50;
        const int maximum_half_power = argc >= 3
            ? parse_nonnegative(argv[2], "maximum_half_power")
            : 30;
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_finite_cyclic_gradient_currents "
                "[maximum_half_level [maximum_half_power]]");
        }

        std::uint64_t fusion_product_checks = 0U;
        std::uint64_t gradient_checks = 0U;
        std::uint64_t antisymmetry_checks = 0U;
        std::uint64_t anchor_checks = 0U;
        std::uint64_t determinant_checks = 0U;
        std::uint64_t shell_identity_checks = 0U;
        std::uint64_t current_checks = 0U;
        std::uint64_t negative_currents = 0U;
        std::uint64_t zero_currents = 0U;
        std::uint64_t all_radius_checks = 0U;
        std::uint64_t negative_all_radius = 0U;
        std::uint64_t zero_all_radius = 0U;
        std::uint64_t adjacent_tp2_checks = 0U;
        std::uint64_t negative_adjacent_tp2 = 0U;
        std::uint64_t negative_adjacent_tp2_after_seeds = 0U;
        std::uint64_t log_concavity_checks = 0U;
        std::uint64_t negative_log_concavity = 0U;
        std::uint64_t shell_sign_sequences = 0U;
        std::uint64_t shell_sign_changes = 0U;
        std::uint64_t shell_negative_to_positive_changes = 0U;
        std::uint64_t shell_multiple_changes = 0U;
        std::uint64_t reduced_shell_negative_to_positive_changes = 0U;
        std::uint64_t reduced_shell_multiple_changes = 0U;
        std::uint64_t chord_checks = 0U;
        std::uint64_t negative_chord_margins = 0U;
        std::uint64_t reduced_chord_checks = 0U;
        std::uint64_t negative_reduced_chord_margins = 0U;
        int maximum_shell_sign_changes = 0;
        std::uint64_t negative_shells = 0U;
        std::uint64_t negative_suffixes = 0U;
        Witness first_negative_current;
        Witness first_negative_all_radius;
        Witness first_negative_adjacent_tp2;
        Witness first_negative_adjacent_tp2_after_seeds;
        Witness first_shell_negative_to_positive;
        Witness first_reduced_shell_negative_to_positive;
        Witness first_negative_chord_margin;
        Witness first_negative_reduced_chord_margin;
        Witness first_negative_shell;
        Witness first_negative_suffix;

        for (int level = 0; level <= maximum_level; ++level) {
            const int order = 2 * level + 2;

            for (int factor = 0; factor <= level; ++factor) {
                for (int source = 0; source <= level; ++source) {
                    Vector cyclic_basis(static_cast<std::size_t>(order), 0);
                    cyclic_basis[static_cast<std::size_t>(source)] = 1;
                    cyclic_basis[static_cast<std::size_t>(
                        residue(-source - 1, order))] = -1;
                    const Vector cyclic_product
                        = cyclic_interval_step(cyclic_basis, factor);

                    Vector finite_basis(
                        static_cast<std::size_t>(level + 1), 0);
                    finite_basis[static_cast<std::size_t>(source)] = 1;
                    const Vector finite_product
                        = finite_step(finite_basis, level, factor);

                    Vector reconstructed(static_cast<std::size_t>(order), 0);
                    for (int target = 0; target <= level; ++target) {
                        const cpp_int& coefficient
                            = finite_product[static_cast<std::size_t>(target)];
                        reconstructed[static_cast<std::size_t>(target)]
                            += coefficient;
                        reconstructed[static_cast<std::size_t>(
                            residue(-target - 1, order))] -= coefficient;
                    }
                    ++fusion_product_checks;
                    if (cyclic_product != reconstructed) {
                        std::cerr
                            << "CYCLIC_BASIS_MODEL_MISMATCH"
                            << " half_level=" << level
                            << " half_factor=" << factor
                            << " half_source=" << source << '\n';
                        return EXIT_FAILURE;
                    }
                }

                Vector finite(static_cast<std::size_t>(level + 1), 0);
                finite[0] = 1;
                Vector cyclic(static_cast<std::size_t>(order), 0);
                cyclic[0] = 1;

                const int maximum_power = 2 * maximum_half_power + 1;
                for (int power = 0; power <= maximum_power; ++power) {
                    const Vector gradient = finite_gradient(cyclic, level);
                    ++gradient_checks;
                    if (gradient != finite) {
                        std::cerr
                            << "CYCLIC_GRADIENT_POWER_MISMATCH"
                            << " half_level=" << level
                            << " half_factor=" << factor
                            << " power=" << power << '\n';
                        return EXIT_FAILURE;
                    }

                    for (int position = 0; position < order; ++position) {
                        ++antisymmetry_checks;
                        const cpp_int reflected = cyclic_gradient_at(
                            cyclic, -position - 1);
                        if (reflected
                            != -cyclic_gradient_at(cyclic, position)) {
                            std::cerr
                                << "CYCLIC_GRADIENT_ANTISYMMETRY_MISMATCH"
                                << " half_level=" << level
                                << " half_factor=" << factor
                                << " power=" << power
                                << " position=" << position << '\n';
                            return EXIT_FAILURE;
                        }
                    }

                    if (power % 2 == 0
                        && power >= 2
                        && factor >= 1) {
                        const int half_power = power / 2;
                        const Vector next_finite
                            = finite_step(finite, level, factor);
                        std::vector<Vector> radius_images;
                        radius_images.reserve(
                            static_cast<std::size_t>(level + 1));
                        for (int radius = 0; radius <= level; ++radius) {
                            radius_images.push_back(
                                finite_step(finite, level, radius));
                        }
                        for (int position = 1; position < level; ++position) {
                            const cpp_int log_concavity_minor
                                = finite[static_cast<std::size_t>(position)]
                                    * finite[
                                        static_cast<std::size_t>(position)]
                                - finite[
                                    static_cast<std::size_t>(position - 1)]
                                    * finite[
                                        static_cast<std::size_t>(
                                            position + 1)];
                            ++log_concavity_checks;
                            if (log_concavity_minor < 0) {
                                ++negative_log_concavity;
                            }
                        }
                        for (int row = 0; row < level; ++row) {
                            for (int column = 0;
                                 column < level;
                                 ++column) {
                                const cpp_int adjacent_minor
                                    = radius_images[
                                        static_cast<std::size_t>(column)][
                                        static_cast<std::size_t>(row)]
                                        * radius_images[
                                            static_cast<std::size_t>(
                                                column + 1)][
                                            static_cast<std::size_t>(
                                                row + 1)]
                                    - radius_images[
                                        static_cast<std::size_t>(
                                            column + 1)][
                                        static_cast<std::size_t>(row)]
                                        * radius_images[
                                            static_cast<std::size_t>(
                                                column)][
                                            static_cast<std::size_t>(
                                                row + 1)];
                                ++adjacent_tp2_checks;
                                if (adjacent_minor < 0) {
                                    ++negative_adjacent_tp2;
                                    if (first_negative_adjacent_tp2.level < 0) {
                                        first_negative_adjacent_tp2 = {
                                            level,
                                            factor,
                                            half_power,
                                            row,
                                            column,
                                            adjacent_minor};
                                    }
                                    if (half_power >= 3) {
                                        ++negative_adjacent_tp2_after_seeds;
                                        if (
                                            first_negative_adjacent_tp2_after_seeds
                                                .level
                                            < 0) {
                                            first_negative_adjacent_tp2_after_seeds
                                                = {
                                                    level,
                                                    factor,
                                                    half_power,
                                                    row,
                                                    column,
                                                    adjacent_minor};
                                        }
                                    }
                                }
                            }
                        }
                        const cpp_int anchor_from_gradient
                            = cyclic_gradient_at(cyclic, factor);
                        ++anchor_checks;
                        if (next_finite[0] != anchor_from_gradient) {
                            std::cerr
                                << "CYCLIC_GRADIENT_ANCHOR_MISMATCH"
                                << " half_level=" << level
                                << " half_factor=" << factor
                                << " half_power=" << half_power << '\n';
                            return EXIT_FAILURE;
                        }

                        for (int target = 0; target <= level; ++target) {
                            const cpp_int endpoint_current
                                = finite[0]
                                    * radius_images[
                                        static_cast<std::size_t>(level)][
                                        static_cast<std::size_t>(target)]
                                - finite[
                                    static_cast<std::size_t>(level)]
                                    * finite[
                                        static_cast<std::size_t>(target)];
                            int previous_shell_sign = 0;
                            int sequence_sign_changes = 0;
                            ++shell_sign_sequences;
                            for (int radius = 1;
                                 radius <= level;
                                 ++radius) {
                                const cpp_int current_radius
                                    = finite[0]
                                        * radius_images[
                                            static_cast<std::size_t>(
                                                radius)][
                                            static_cast<std::size_t>(
                                                target)]
                                    - finite[
                                        static_cast<std::size_t>(radius)]
                                        * finite[
                                            static_cast<std::size_t>(
                                                target)];
                                const cpp_int previous_radius
                                    = radius == 1
                                    ? cpp_int(0)
                                    : finite[0]
                                        * radius_images[
                                            static_cast<std::size_t>(
                                                radius - 1)][
                                            static_cast<std::size_t>(
                                                target)]
                                        - finite[
                                            static_cast<std::size_t>(
                                                radius - 1)]
                                            * finite[
                                                static_cast<std::size_t>(
                                                    target)];
                                const cpp_int shell
                                    = current_radius - previous_radius;
                                const int shell_sign
                                    = shell < 0 ? -1 : (shell > 0 ? 1 : 0);
                                if (shell_sign != 0) {
                                    if (previous_shell_sign != 0
                                        && shell_sign
                                            != previous_shell_sign) {
                                        ++shell_sign_changes;
                                        ++sequence_sign_changes;
                                        if (previous_shell_sign < 0
                                            && shell_sign > 0) {
                                            ++shell_negative_to_positive_changes;
                                            if (
                                                first_shell_negative_to_positive
                                                    .level
                                                < 0) {
                                                first_shell_negative_to_positive
                                                    = {
                                                        level,
                                                        factor,
                                                        half_power,
                                                        target,
                                                        radius,
                                                        shell};
                                            }
                                            if (2 * factor <= level) {
                                                ++reduced_shell_negative_to_positive_changes;
                                                if (
                                                    first_reduced_shell_negative_to_positive
                                                        .level
                                                    < 0) {
                                                    first_reduced_shell_negative_to_positive
                                                        = {
                                                            level,
                                                            factor,
                                                            half_power,
                                                            target,
                                                            radius,
                                                            shell};
                                                }
                                            }
                                        }
                                    }
                                    previous_shell_sign = shell_sign;
                                }
                            }
                            maximum_shell_sign_changes = std::max(
                                maximum_shell_sign_changes,
                                sequence_sign_changes);
                            if (sequence_sign_changes > 1) {
                                ++shell_multiple_changes;
                                if (2 * factor <= level) {
                                    ++reduced_shell_multiple_changes;
                                }
                            }

                            for (int radius = 0;
                                 radius <= level;
                                 ++radius) {
                                const cpp_int general_current
                                    = finite[0]
                                        * radius_images[
                                            static_cast<std::size_t>(
                                                radius)][
                                            static_cast<std::size_t>(
                                                target)]
                                    - finite[
                                        static_cast<std::size_t>(radius)]
                                        * finite[
                                            static_cast<std::size_t>(
                                                target)];
                                ++all_radius_checks;
                                if (general_current < 0) {
                                    ++negative_all_radius;
                                    if (first_negative_all_radius.level < 0) {
                                        first_negative_all_radius = {
                                            level,
                                            factor,
                                            half_power,
                                            target,
                                            radius,
                                            general_current};
                                    }
                                } else if (general_current == 0) {
                                    ++zero_all_radius;
                                }
                                const cpp_int chord_margin
                                    = level * general_current
                                    - radius * endpoint_current;
                                ++chord_checks;
                                if (chord_margin < 0) {
                                    ++negative_chord_margins;
                                    if (
                                        first_negative_chord_margin.level
                                        < 0) {
                                        first_negative_chord_margin = {
                                            level,
                                            factor,
                                            half_power,
                                            target,
                                            radius,
                                            chord_margin};
                                    }
                                }
                                if (2 * factor <= level) {
                                    ++reduced_chord_checks;
                                    if (chord_margin < 0) {
                                        ++negative_reduced_chord_margins;
                                        if (
                                            first_negative_reduced_chord_margin
                                                .level
                                            < 0) {
                                            first_negative_reduced_chord_margin
                                                = {
                                                    level,
                                                    factor,
                                                    half_power,
                                                    target,
                                                    radius,
                                                    chord_margin};
                                        }
                                    }
                                }
                            }

                            const cpp_int determinant
                                = finite[0]
                                    * next_finite[
                                        static_cast<std::size_t>(target)]
                                - next_finite[0]
                                    * finite[
                                        static_cast<std::size_t>(target)];
                            cpp_int shell_sum = 0;
                            cpp_int previous_prefix = 0;
                            for (int radius = 1;
                                 radius <= factor;
                                 ++radius) {
                                const cpp_int shell
                                    = finite[0]
                                        * (
                                            cyclic_gradient_at(
                                                cyclic, target - radius)
                                            + cyclic_gradient_at(
                                                cyclic, target + radius))
                                    + (
                                        finite[
                                            static_cast<std::size_t>(
                                                radius - 1)]
                                        - finite[
                                            static_cast<std::size_t>(
                                                radius)])
                                        * finite[
                                            static_cast<std::size_t>(
                                                target)];
                                shell_sum += shell;
                                const cpp_int direct_prefix
                                    = finite[0]
                                        * radius_images[
                                            static_cast<std::size_t>(
                                                radius)][
                                                static_cast<std::size_t>(
                                                    target)]
                                    - finite[
                                        static_cast<std::size_t>(radius)]
                                        * finite[
                                            static_cast<std::size_t>(
                                                target)];
                                ++shell_identity_checks;
                                if (shell_sum != direct_prefix
                                    || shell_sum - previous_prefix != shell) {
                                    std::cerr
                                        << "CUMULATIVE_SHELL_IDENTITY_MISMATCH"
                                        << " half_level=" << level
                                        << " half_factor=" << factor
                                        << " half_power=" << half_power
                                        << " half_target=" << target
                                        << " radius=" << radius << '\n';
                                    return EXIT_FAILURE;
                                }
                                previous_prefix = shell_sum;

                                ++current_checks;
                                if (direct_prefix < 0) {
                                    ++negative_currents;
                                    if (first_negative_current.level < 0) {
                                        first_negative_current = {
                                            level,
                                            factor,
                                            half_power,
                                            target,
                                            radius,
                                            direct_prefix};
                                    }
                                } else if (direct_prefix == 0) {
                                    ++zero_currents;
                                }
                                if (shell < 0) {
                                    ++negative_shells;
                                    if (first_negative_shell.level < 0) {
                                        first_negative_shell = {
                                            level,
                                            factor,
                                            half_power,
                                            target,
                                            radius,
                                            shell};
                                    }
                                }
                            }
                            ++determinant_checks;
                            if (shell_sum != determinant) {
                                std::cerr
                                    << "ANCHORED_DETERMINANT_SHELL_MISMATCH"
                                    << " half_level=" << level
                                    << " half_factor=" << factor
                                    << " half_power=" << half_power
                                    << " half_target=" << target << '\n';
                                return EXIT_FAILURE;
                            }

                            cpp_int prefix = 0;
                            for (int radius = 1;
                                 radius <= factor;
                                 ++radius) {
                                const cpp_int shell
                                    = finite[0]
                                        * (
                                            cyclic_gradient_at(
                                                cyclic, target - radius)
                                            + cyclic_gradient_at(
                                                cyclic, target + radius))
                                    + (
                                        finite[
                                            static_cast<std::size_t>(
                                                radius - 1)]
                                        - finite[
                                            static_cast<std::size_t>(
                                                radius)])
                                        * finite[
                                            static_cast<std::size_t>(
                                                target)];
                                prefix += shell;
                                const cpp_int suffix = determinant - prefix;
                                if (radius < factor && suffix < 0) {
                                    ++negative_suffixes;
                                    if (first_negative_suffix.level < 0) {
                                        first_negative_suffix = {
                                            level,
                                            factor,
                                            half_power,
                                            target,
                                            radius + 1,
                                            suffix};
                                    }
                                }
                            }
                        }
                    }

                    if (power < maximum_power) {
                        finite = finite_step(finite, level, factor);
                        cyclic = cyclic_interval_step(cyclic, factor);
                    }
                }
            }
        }

        std::cout
            << "SU2_FINITE_CYCLIC_GRADIENT_CURRENTS"
            << " maximum_level=" << 2 * maximum_level
            << " maximum_half_power=" << maximum_half_power
            << " fusion_product_checks=" << fusion_product_checks
            << " gradient_checks=" << gradient_checks
            << " antisymmetry_checks=" << antisymmetry_checks
            << " anchor_checks=" << anchor_checks
            << " determinant_checks=" << determinant_checks
            << " shell_identity_checks=" << shell_identity_checks
            << " current_checks=" << current_checks
            << " negative_currents=" << negative_currents
            << " zero_currents=" << zero_currents
            << " all_radius_checks=" << all_radius_checks
            << " negative_all_radius=" << negative_all_radius
            << " zero_all_radius=" << zero_all_radius
            << " adjacent_tp2_checks=" << adjacent_tp2_checks
            << " negative_adjacent_tp2=" << negative_adjacent_tp2
            << " negative_adjacent_tp2_after_seeds="
            << negative_adjacent_tp2_after_seeds
            << " log_concavity_checks=" << log_concavity_checks
            << " negative_log_concavity=" << negative_log_concavity
            << " shell_sign_sequences=" << shell_sign_sequences
            << " shell_sign_changes=" << shell_sign_changes
            << " shell_negative_to_positive_changes="
            << shell_negative_to_positive_changes
            << " shell_multiple_changes=" << shell_multiple_changes
            << " reduced_shell_negative_to_positive_changes="
            << reduced_shell_negative_to_positive_changes
            << " reduced_shell_multiple_changes="
            << reduced_shell_multiple_changes
            << " chord_checks=" << chord_checks
            << " negative_chord_margins="
            << negative_chord_margins
            << " reduced_chord_checks=" << reduced_chord_checks
            << " negative_reduced_chord_margins="
            << negative_reduced_chord_margins
            << " maximum_shell_sign_changes="
            << maximum_shell_sign_changes
            << " negative_shells=" << negative_shells
            << " negative_suffixes=" << negative_suffixes;
        if (first_negative_current.level >= 0) {
            print_witness("first_negative_current", first_negative_current);
        }
        if (first_negative_all_radius.level >= 0) {
            print_witness(
                "first_negative_all_radius",
                first_negative_all_radius);
        }
        if (first_negative_adjacent_tp2.level >= 0) {
            print_witness(
                "first_negative_adjacent_tp2",
                first_negative_adjacent_tp2);
        }
        if (first_negative_adjacent_tp2_after_seeds.level >= 0) {
            print_witness(
                "first_negative_adjacent_tp2_after_seeds",
                first_negative_adjacent_tp2_after_seeds);
        }
        if (first_shell_negative_to_positive.level >= 0) {
            print_witness(
                "first_shell_negative_to_positive",
                first_shell_negative_to_positive);
        }
        if (first_reduced_shell_negative_to_positive.level >= 0) {
            print_witness(
                "first_reduced_shell_negative_to_positive",
                first_reduced_shell_negative_to_positive);
            print_shell_and_prefix_sequences(
                "first_reduced_shell_negative_to_positive",
                first_reduced_shell_negative_to_positive);
        }
        if (first_negative_chord_margin.level >= 0) {
            print_witness(
                "first_negative_chord_margin",
                first_negative_chord_margin);
        }
        if (first_negative_reduced_chord_margin.level >= 0) {
            print_witness(
                "first_negative_reduced_chord_margin",
                first_negative_reduced_chord_margin);
        }
        if (first_negative_shell.level >= 0) {
            print_witness("first_negative_shell", first_negative_shell);
        }
        if (first_negative_suffix.level >= 0) {
            print_witness("first_negative_suffix", first_negative_suffix);
        }
        std::cout
            << " model_result=PASS_EXACT_IDENTITIES"
            << " current_result="
            << (negative_currents == 0U
                    ? "NO_NEGATIVE_CUMULATIVE_CURRENT"
                    : "NEGATIVE_CUMULATIVE_CURRENT")
            << " all_radius_result="
            << (negative_all_radius == 0U
                    ? "NO_NEGATIVE_ANCHORED_MINOR"
                    : "NEGATIVE_ANCHORED_MINOR")
            << " tp2_result="
            << (negative_adjacent_tp2 == 0U
                    ? "NO_NEGATIVE_ADJACENT_MINOR"
                    : "NEGATIVE_ADJACENT_MINOR")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
