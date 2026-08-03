#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Vector = std::vector<Integer>;

int parse_positive(const char* text, const char* name) {
    const std::string input(text);
    std::size_t consumed = 0U;
    const long parsed = std::stol(input, &consumed, 10);
    if (consumed != input.size() || parsed <= 0L) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

Vector fuse(const Vector& input, const int level, const int factor) {
    Vector output(input.size(), 0);
    for (int source = 0; source <= level; ++source) {
        const Integer& coefficient = input[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        const int lower = std::abs(source - factor);
        const int upper = std::min(source + factor,
                                   2 * level - source - factor);
        for (int target = lower; target <= upper; target += 2) {
            output[static_cast<std::size_t>(target)] += coefficient;
        }
    }
    return output;
}

Vector square_profile(const Vector& root, const int level) {
    Vector square(root.size(), 0);
    for (int factor = 0; factor <= level; ++factor) {
        const Integer& coefficient = root[static_cast<std::size_t>(factor)];
        if (coefficient == 0) {
            continue;
        }
        const Vector image = fuse(root, level, factor);
        for (std::size_t index = 0U; index < square.size(); ++index) {
            square[index] += coefficient * image[index];
        }
    }
    return square;
}

void print_vector(const Vector& values) {
    std::cout << '[';
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << ']';
}

struct Counters {
    std::uint64_t roots = 0U;
    std::uint64_t boundary_admissible = 0U;
    std::uint64_t residual_checks = 0U;
    std::uint64_t negative_residuals = 0U;
    std::uint64_t negative_total_updates = 0U;
    std::uint64_t update_identity_mismatches = 0U;
    bool have_first_negative_residual = false;
    int first_level = -1;
    int first_parity = -1;
    int first_radius = -1;
    Vector first_root;
    Vector first_square;
    Integer first_reserve = 0;
    Integer first_residual = 0;
    Integer first_total = 0;
};

void inspect(const Vector& root, const int level, const int parity,
             Counters& counters) {
    const int half_level = level / 2;
    const Vector square = square_profile(root, level);
    const Vector inserted_root = fuse(root, level, 1);
    const Vector inserted_square = square_profile(inserted_root, level);
    Vector e(static_cast<std::size_t>(half_level + 1), 0);
    Vector inserted_e(static_cast<std::size_t>(half_level + 1), 0);
    for (int index = 0; index <= half_level; ++index) {
        e[static_cast<std::size_t>(index)] =
            square[static_cast<std::size_t>(2 * index)];
        inserted_e[static_cast<std::size_t>(index)] =
            inserted_square[static_cast<std::size_t>(2 * index)];
    }
    Vector boundary(static_cast<std::size_t>(half_level + 1), 0);
    bool admissible = true;
    for (int radius = 0; radius <= half_level; ++radius) {
        boundary[static_cast<std::size_t>(radius)] =
            e[0] * e[static_cast<std::size_t>(half_level - radius)]
            - e[static_cast<std::size_t>(radius)]
                * e[static_cast<std::size_t>(half_level)];
        if (boundary[static_cast<std::size_t>(radius)] < 0) {
            admissible = false;
        }
    }
    ++counters.roots;
    if (!admissible) {
        return;
    }
    ++counters.boundary_admissible;
    for (int radius = 1; radius < half_level; ++radius) {
        const Integer c_radius =
            e[static_cast<std::size_t>(radius - 1)]
            + 2 * e[static_cast<std::size_t>(radius)]
            + e[static_cast<std::size_t>(radius + 1)];
        const int reflected = half_level - radius;
        const Integer c_reflected =
            e[static_cast<std::size_t>(reflected - 1)]
            + 2 * e[static_cast<std::size_t>(reflected)]
            + e[static_cast<std::size_t>(reflected + 1)];
        const Integer reserve =
            boundary[static_cast<std::size_t>(radius - 1)]
            + 2 * boundary[static_cast<std::size_t>(radius)]
            + boundary[static_cast<std::size_t>(radius + 1)];
        const Integer residual = e[1] * c_reflected
            - e[static_cast<std::size_t>(half_level - 1)] * c_radius;
        const Integer total = reserve + residual;
        const Integer direct = inserted_e[0]
            * inserted_e[static_cast<std::size_t>(half_level - radius)]
            - inserted_e[static_cast<std::size_t>(radius)]
                * inserted_e[static_cast<std::size_t>(half_level)];
        ++counters.residual_checks;
        if (direct != total) {
            ++counters.update_identity_mismatches;
        }
        if (residual < 0) {
            ++counters.negative_residuals;
            if (!counters.have_first_negative_residual) {
                counters.have_first_negative_residual = true;
                counters.first_level = level;
                counters.first_parity = parity;
                counters.first_radius = radius;
                counters.first_root = root;
                counters.first_square = square;
                counters.first_reserve = reserve;
                counters.first_residual = residual;
                counters.first_total = total;
            }
        }
        if (total < 0) {
            ++counters.negative_total_updates;
        }
    }
}

void enumerate_coefficients(const int position, const int length,
                            const int maximum_coordinate, const bool falling,
                            const int previous, std::vector<int>& values,
                            const int level, const int parity,
                            const int first, Counters& counters) {
    if (position == length) {
        Vector root(static_cast<std::size_t>(level + 1), 0);
        for (int index = 0; index < length; ++index) {
            root[static_cast<std::size_t>(first + 2 * index)] =
                values[static_cast<std::size_t>(index)];
        }
        inspect(root, level, parity, counters);
        return;
    }
    for (int coefficient = 1; coefficient <= maximum_coordinate;
         ++coefficient) {
        if (position > 0 && falling && coefficient > previous) {
            continue;
        }
        const bool next_falling = falling
            || (position > 0 && coefficient < previous);
        values[static_cast<std::size_t>(position)] = coefficient;
        enumerate_coefficients(position + 1, length, maximum_coordinate,
                               next_falling, coefficient, values, level,
                               parity, first, counters);
    }
}

void enumerate_roots(const int level, const int maximum_coordinate,
                     Counters& counters) {
    for (int parity = 0; parity <= 1; ++parity) {
        const int first_label = parity;
        const int last_label = level - ((level - parity) & 1);
        for (int first = first_label; first <= last_label; first += 2) {
            for (int last = first; last <= last_label; last += 2) {
                const int length = (last - first) / 2 + 1;
                std::vector<int> values(static_cast<std::size_t>(length), 0);
                enumerate_coefficients(0, length, maximum_coordinate, false,
                                       0, values, level, parity, first,
                                       counters);
            }
        }
    }
}

void replay_structured_residual_obstruction() {
    constexpr int level = 6;
    const Vector root{1, 0, 1, 0, 1, 0, 0};
    const Vector expected_square{3, 0, 6, 0, 6, 0, 2};
    const Vector square = square_profile(root, level);
    if (square != expected_square) {
        throw std::runtime_error("structured square replay mismatch");
    }
    const Vector e{3, 6, 6, 2};
    const Vector boundary{0, 6, 6, 0};
    const Integer c_one = e[0] + 2 * e[1] + e[2];
    const Integer c_two = e[1] + 2 * e[2] + e[3];
    const Integer reserve = boundary[0] + 2 * boundary[1] + boundary[2];
    const Integer residual = e[1] * c_two - e[2] * c_one;
    const Integer total = reserve + residual;
    if (c_one != 21 || c_two != 20 || reserve != 18 || residual != -6
        || total != 12) {
        throw std::runtime_error("structured residual replay mismatch");
    }
    std::cout << "SU2_USM_STRUCTURED_RESIDUAL_OBSTRUCTION"
              << " level=6 root=";
    print_vector(root);
    std::cout << " square=";
    print_vector(square);
    std::cout << " C1=" << c_one
              << " C2=" << c_two
              << " reserve=" << reserve
              << " residual=" << residual
              << " total=" << total
              << " result=PASS_EXACT\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2
            && std::string(argv[1]) == "--replay-structured-residual") {
            replay_structured_residual_obstruction();
            return EXIT_SUCCESS;
        }
        const int maximum_level = argc >= 2
            ? parse_positive(argv[1], "maximum_half_level") : 12;
        const int maximum_coordinate = argc >= 3
            ? parse_positive(argv[2], "maximum_coordinate") : 4;
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_usm_reserve "
                "[maximum_half_level] [maximum_coordinate] | "
                "--replay-structured-residual");
        }
        Counters counters;
        for (int level = 2; level <= maximum_level; level += 2) {
            enumerate_roots(level, maximum_coordinate, counters);
        }
        std::cout << "SU2_USM_FUNDAMENTAL_RESERVE"
                  << " maximum_half_level=" << maximum_level
                  << " maximum_coordinate=" << maximum_coordinate
                  << " roots=" << counters.roots
                  << " boundary_admissible=" << counters.boundary_admissible
                  << " residual_checks=" << counters.residual_checks
                  << " negative_residuals=" << counters.negative_residuals
                  << " negative_total_updates="
                  << counters.negative_total_updates
                  << " update_identity_mismatches="
                  << counters.update_identity_mismatches;
        if (counters.have_first_negative_residual) {
            std::cout << " first_level=" << counters.first_level
                      << " first_parity=" << counters.first_parity
                      << " first_radius=" << counters.first_radius
                      << " first_root=";
            print_vector(counters.first_root);
            std::cout << " first_square=";
            print_vector(counters.first_square);
            std::cout << " first_reserve=" << counters.first_reserve
                      << " first_residual=" << counters.first_residual
                      << " first_total=" << counters.first_total;
        }
        std::cout << " result="
                  << (counters.negative_total_updates == 0U
                          ? "TOTAL_UPDATES_NONNEGATIVE"
                          : "TOTAL_UPDATE_COUNTEREXAMPLE")
                  << '\n';
        return counters.negative_total_updates == 0U
                    && counters.update_identity_mismatches == 0U
            ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
