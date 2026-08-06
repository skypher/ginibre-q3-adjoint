// Exact common-source evaluator for the separated-lower H_2 rail current.
//
// The crossing weights have the closed form recorded in Lemma
// 5A8H28P3J.  Consequently this current is quadratic in Q; this source
// evaluates its three exact coefficients directly from ordinary N_d paths.

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

int positive_modulo(int value, int modulus) {
    return (value % modulus + modulus) % modulus;
}

Integer binomial_integer(int top, int bottom) {
    if (bottom < 0 || top < bottom) {
        return 0;
    }
    Integer result = 1;
    for (int index = 1; index <= bottom; ++index) {
        result *= top - bottom + index;
        result /= index;
    }
    return result;
}

Integer ordinary_profile_formula(int power, int d, int y, int label) {
    if (
        power < 1
        || label < 0
        || (label - y - power * d) % 2 != 0
    ) {
        return 0;
    }
    const int r = (y + power * d - label) / 2;
    if (r < 0) {
        return 0;
    }
    Integer result = 0;
    for (int j = 0; j <= power; ++j) {
        for (int k = 0; k <= 1; ++k) {
            const int top = r - j * (d + 1) - k * (y + 1)
                + power - 1;
            const Integer term = binomial_integer(top, power - 1);
            result += ((j + k) % 2 == 0 ? 1 : -1)
                * binomial_integer(power, j) * term;
        }
    }
    return result;
}

std::vector<Integer> multiply_by_character(
    const std::vector<Integer>& state,
    int d
) {
    std::vector<Integer> result(state.size());
    for (std::size_t source = 0U; source < state.size(); ++source) {
        if (state[source] == 0) {
            continue;
        }
        const int source_label = static_cast<int>(source);
        for (int target = std::abs(source_label - d);
             target <= source_label + d;
             target += 2) {
            if (static_cast<std::size_t>(target) >= result.size()) {
                throw std::runtime_error("ordinary support overflow");
            }
            result[static_cast<std::size_t>(target)] += state[source];
        }
    }
    return result;
}

struct Quadratic {
    Integer square = 0;
    Integer linear = 0;
    Integer constant = 0;

    Integer evaluate(int q) const {
        return square * q * q + linear * q + constant;
    }
};

struct ProfileMoments {
    Integer p20 = 0;
    Integer p21 = 0;
    Integer p22 = 0;
    Integer p30 = 0;
    Integer p31 = 0;
    Integer p40 = 0;
};

ProfileMoments profile_moments(int q, int d, int y, int rail) {
    const int maximum_label = y + 4 * d;
    std::vector<Integer> profile(
        static_cast<std::size_t>(maximum_label + 1)
    );
    profile[static_cast<std::size_t>(y)] = 1;
    profile = multiply_by_character(profile, d);
    profile = multiply_by_character(profile, d);
    const std::vector<Integer> p2 = profile;
    profile = multiply_by_character(profile, d);
    const std::vector<Integer> p3 = profile;
    profile = multiply_by_character(profile, d);
    const std::vector<Integer> p4 = profile;
    const int rho = d + 1 + positive_modulo(rail - (d + 1), 4);
    const int p2_residue = positive_modulo(2 * q + d - 2 * rail, 8);
    const int p3_residue = positive_modulo(2 * rail, 8);
    ProfileMoments result;
    for (int label = 0; label <= maximum_label; ++label) {
        if (label % 8 == p2_residue) {
            const Integer& p2_value = p2[static_cast<std::size_t>(label)];
            result.p20 += p2_value;
            result.p21 += label * p2_value;
            result.p22 += label * label * p2_value;
            result.p40 += p4[static_cast<std::size_t>(label)];
        }
        if (label % 8 == p3_residue && label >= 2 * rho) {
            const Integer& p3_value = p3[static_cast<std::size_t>(label)];
            result.p30 += p3_value;
            result.p31 += label * p3_value;
        }
    }
    return result;
}

Quadratic lower_current(int q, int d, int y, int rail) {
    if (d < 0 || y < 0 || rail < 0 || rail >= 4) {
        throw std::runtime_error("invalid lower-current parameters");
    }
    if ((d - y) % 2 != 0) {
        throw std::runtime_error("source has incompatible parity");
    }
    const int maximum_label = y + 4 * d;
    std::vector<Integer> profile(
        static_cast<std::size_t>(maximum_label + 1)
    );
    profile[static_cast<std::size_t>(y)] = 1;
    profile = multiply_by_character(profile, d);
    profile = multiply_by_character(profile, d);
    const std::vector<Integer> p2 = profile;
    profile = multiply_by_character(profile, d);
    const std::vector<Integer> p3 = profile;
    profile = multiply_by_character(profile, d);
    const std::vector<Integer> p4 = profile;

    const int rho = d + 1
        + positive_modulo(rail - (d + 1), 4);
    Quadratic result;
    for (int label = 0; label <= maximum_label; ++label) {
        if (label % 8 == positive_modulo(2 * rail, 8) && label >= 2 * rho) {
            // M_1(L/2)=-4Q^2+4(L-2d-2)Q+2L+4d^2+4d.
            result.square -= 4 * p3[static_cast<std::size_t>(label)];
            result.linear += 4 * (label - 2 * d - 2)
                * p3[static_cast<std::size_t>(label)];
            result.constant += (2 * label + 4 * d * d + 4 * d)
                * p3[static_cast<std::size_t>(label)];
        }
        if (label % 8 == positive_modulo(2 * q + d - 2 * rail, 8)) {
            // d_2=4Q-g, d_3=4Q^2+(4d+8)Q-h.
            const int g = label <= d ? d + label : 2 * label;
            const Integer h = Integer(label) * label
                + Integer(2 * d + 4) * label
                + Integer(3 * d * d + 6 * d + 2);
            const Integer p2_value = p2[static_cast<std::size_t>(label)];
            result.square += (2 * g + 4) * p2_value;
            result.linear += (
                8 * d * d + 2 * d * g + 4 * g + 16 * d + 8
                - 2 * h
            ) * p2_value;
            result.constant -= (
                h + Integer(2 * d * d + 3 * d) * g
            ) * p2_value;
            // f_4d_1=4Q+2.
            result.linear += 4 * p4[static_cast<std::size_t>(label)];
            result.constant += 2 * p4[static_cast<std::size_t>(label)];
        }
    }
    return result;
}

void bounded_scan(int maximum_d, int y_multiple) {
    if (maximum_d < 11 || y_multiple < 0) {
        throw std::runtime_error("invalid scan bounds");
    }
    std::uint64_t cases = 0U;
    std::uint64_t negative_values = 0U;
    std::uint64_t negative_four_step_increments = 0U;
    bool has_minimum_increment = false;
    Integer minimum_increment = 0;
    int witness_d = -1;
    int witness_y = -1;
    int witness_z = -1;
    int witness_rail = -1;
    for (int d = 11; d <= maximum_d; ++d) {
        for (int y = d % 2; y <= y_multiple * d; y += 2) {
            for (int z = 0; z < 4; ++z) {
                const int q = y + 3 * d + 1 + z;
                for (int rail = 0; rail < 4; ++rail) {
                    const Quadratic current = lower_current(q, d, y, rail);
                    const Integer value = current.evaluate(q);
                    const Integer increment = current.evaluate(q + 4) - value;
                    ++cases;
                    if (value < 0) {
                        ++negative_values;
                    }
                    if (increment < 0) {
                        ++negative_four_step_increments;
                    }
                    if (!has_minimum_increment || increment < minimum_increment) {
                        has_minimum_increment = true;
                        minimum_increment = increment;
                        witness_d = d;
                        witness_y = y;
                        witness_z = z;
                        witness_rail = rail;
                    }
                }
            }
        }
    }
    std::cout
        << "SU2_SHELL_H2_LOWER_CURRENT_SCAN"
        << " maximum_d=" << maximum_d
        << " y_multiple=" << y_multiple
        << " cases=" << cases
        << " negative_values=" << negative_values
        << " negative_four_step_increments="
        << negative_four_step_increments
        << " minimum_four_step_increment=" << minimum_increment
        << " minimum_four_step_increment_witness=d:" << witness_d
        << ",y:" << witness_y
        << ",z:" << witness_z
        << ",rail:" << witness_rail
        << " result="
        << (
            negative_values == 0U && negative_four_step_increments == 0U
                ? "PASS_BOUNDED_DISCOVERY"
                : "FAIL_BOUNDED_DISCOVERY"
        ) << '\n';
}

void profile_audit(int maximum_d, int y_multiple) {
    if (maximum_d < 1 || y_multiple < 0) {
        throw std::runtime_error("invalid profile-audit bounds");
    }
    std::uint64_t entries = 0U;
    for (int d = 1; d <= maximum_d; ++d) {
        for (int y = d % 2; y <= y_multiple * d; y += 2) {
            std::vector<Integer> profile(
                static_cast<std::size_t>(y + 4 * d + 1)
            );
            profile[static_cast<std::size_t>(y)] = 1;
            for (int power = 1; power <= 4; ++power) {
                profile = multiply_by_character(profile, d);
                for (int label = 0; label <= y + power * d; ++label) {
                    if (
                        profile[static_cast<std::size_t>(label)]
                        != ordinary_profile_formula(power, d, y, label)
                    ) {
                        throw std::runtime_error(
                            "ordinary terminal-profile formula mismatch"
                        );
                    }
                    ++entries;
                }
            }
        }
    }
    std::cout
        << "SU2_SHELL_H2_LOWER_PROFILE_AUDIT"
        << " maximum_d=" << maximum_d
        << " y_multiple=" << y_multiple
        << " entries=" << entries
        << " result=PASS_FORMULA_AUDIT\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4 && std::string(argv[1]) == "--scan") {
            bounded_scan(std::atoi(argv[2]), std::atoi(argv[3]));
            return EXIT_SUCCESS;
        }
        if (argc == 4 && std::string(argv[1]) == "--profile-audit") {
            profile_audit(std::atoi(argv[2]), std::atoi(argv[3]));
            return EXIT_SUCCESS;
        }
        if (argc == 6 && std::string(argv[1]) == "--moments") {
            const int q = std::atoi(argv[2]);
            const int d = std::atoi(argv[3]);
            const int y = std::atoi(argv[4]);
            const int rail = std::atoi(argv[5]);
            const ProfileMoments moments = profile_moments(q, d, y, rail);
            std::cout
                << "SU2_SHELL_H2_LOWER_PROFILE_MOMENTS"
                << " Q=" << q
                << " d=" << d
                << " y=" << y
                << " rail=" << rail
                << " p20=" << moments.p20
                << " p21=" << moments.p21
                << " p22=" << moments.p22
                << " p30=" << moments.p30
                << " p31=" << moments.p31
                << " p40=" << moments.p40 << '\n';
            return EXIT_SUCCESS;
        }
        if (argc != 5) {
            throw std::runtime_error(
                "usage: Q D Y RAIL | --scan MAX_D Y_MULTIPLE | "
                "--profile-audit MAX_D Y_MULTIPLE"
            );
        }
        const int q = std::atoi(argv[1]);
        const int d = std::atoi(argv[2]);
        const int y = std::atoi(argv[3]);
        const int rail = std::atoi(argv[4]);
        if (q - y <= 3 * d || q <= 0) {
            throw std::runtime_error("input is not strictly separated-lower");
        }
        const Quadratic current = lower_current(q, d, y, rail);
        std::cout
            << "SU2_SHELL_H2_LOWER_CURRENT"
            << " Q=" << q
            << " d=" << d
            << " y=" << y
            << " rail=" << rail
            << " square=" << current.square
            << " linear=" << current.linear
            << " constant=" << current.constant
            << " value=" << current.evaluate(q)
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
