#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Atom {
    long double eigenvalue = 0.0L;
    long double residue = 0.0L;
};

int parse_bounded(
    const char* text,
    int minimum,
    int maximum,
    const std::string& name
) {
    const std::string value{text};
    std::size_t consumed = 0;
    const long parsed = std::stol(value, &consumed);
    if (
        consumed != value.size()
        || parsed < minimum
        || parsed > maximum
    ) {
        throw std::runtime_error(
            name + " must lie in ["
            + std::to_string(minimum) + ","
            + std::to_string(maximum) + "]"
        );
    }
    return static_cast<int>(parsed);
}

bool same_eigenvalue(long double left, long double right) {
    constexpr long double tolerance = 2.0e-15L;
    return std::abs(left - right)
        <= tolerance
            * std::max(
                1.0L,
                std::max(std::abs(left), std::abs(right))
            );
}

std::vector<long double> mode_vector(int level, int mode) {
    const long double pi = std::acos(-1.0L);
    const long double normalization =
        std::sqrt(2.0L / static_cast<long double>(level + 2));
    std::vector<long double> result(
        static_cast<std::size_t>(level + 1)
    );
    for (int label = 0; label <= level; ++label) {
        result[static_cast<std::size_t>(label)] =
            normalization
            * std::sin(
                static_cast<long double>((label + 1) * (mode + 1))
                * pi
                / static_cast<long double>(level + 2)
            );
    }
    return result;
}

long double fusion_eigenvalue(int level, int factor, int mode) {
    const long double pi = std::acos(-1.0L);
    const long double theta =
        static_cast<long double>(mode + 1)
        * pi
        / static_cast<long double>(level + 2);
    return
        std::sin(static_cast<long double>(factor + 1) * theta)
        / std::sin(theta);
}

long double wedge_coordinate(
    const std::vector<long double>& left_mode,
    const std::vector<long double>& right_mode,
    int first,
    int second
) {
    return
        left_mode[static_cast<std::size_t>(first)]
            * right_mode[static_cast<std::size_t>(second)]
        - right_mode[static_cast<std::size_t>(first)]
            * left_mode[static_cast<std::size_t>(second)];
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level =
            argc >= 2
                ? parse_bounded(argv[1], 4, 400, "maximum_level")
                : 80;
        if (argc > 2) {
            throw std::runtime_error(
                "usage: probe_su2_finite_anchored_spectral_tails "
                "[maximum_level]"
            );
        }

        std::size_t parameter_rows = 0;
        std::size_t target_rows = 0;
        std::size_t grouped_atoms = 0;
        long double minimum_tail = 0.0L;
        bool minimum_initialized = false;

        for (int level = 4; level <= maximum_level; level += 2) {
            std::vector<std::vector<long double>> modes;
            modes.reserve(static_cast<std::size_t>(level + 1));
            for (int mode = 0; mode <= level; ++mode) {
                modes.push_back(mode_vector(level, mode));
            }

            for (
                int factor = 2;
                2 * factor < level;
                factor += 2
            ) {
                ++parameter_rows;
                std::vector<long double> eigenvalues(
                    static_cast<std::size_t>(level + 1)
                );
                for (int mode = 0; mode <= level; ++mode) {
                    eigenvalues[static_cast<std::size_t>(mode)] =
                        fusion_eigenvalue(level, factor, mode);
                }

                for (int target = 2; target <= level; target += 2) {
                    ++target_rows;
                    std::vector<Atom> atoms;
                    for (int left = 0; left <= level; ++left) {
                        for (
                            int right = left + 1;
                            right <= level;
                            ++right
                        ) {
                            const long double product =
                                eigenvalues[
                                    static_cast<std::size_t>(left)
                                ]
                                * eigenvalues[
                                    static_cast<std::size_t>(right)
                                ];
                            const long double source =
                                wedge_coordinate(
                                    modes[
                                        static_cast<std::size_t>(left)
                                    ],
                                    modes[
                                        static_cast<std::size_t>(right)
                                    ],
                                    0,
                                    factor
                                );
                            const long double destination =
                                wedge_coordinate(
                                    modes[
                                        static_cast<std::size_t>(left)
                                    ],
                                    modes[
                                        static_cast<std::size_t>(right)
                                    ],
                                    0,
                                    target
                                );
                            atoms.push_back(
                                Atom{
                                    product * product,
                                    source * destination
                                }
                            );
                        }
                    }
                    std::sort(
                        atoms.begin(),
                        atoms.end(),
                        [](const Atom& left, const Atom& right) {
                            return left.eigenvalue < right.eigenvalue;
                        }
                    );

                    std::vector<Atom> groups;
                    for (const Atom& atom : atoms) {
                        if (
                            groups.empty()
                            || !same_eigenvalue(
                                groups.back().eigenvalue,
                                atom.eigenvalue
                            )
                        ) {
                            groups.push_back(atom);
                        } else {
                            groups.back().residue += atom.residue;
                        }
                    }
                    grouped_atoms += groups.size();

                    long double tail = 0.0L;
                    long double scale = 0.0L;
                    for (
                        std::size_t index = groups.size();
                        index > 0U;
                        --index
                    ) {
                        tail += groups[index - 1U].residue;
                        scale += std::abs(groups[index - 1U].residue);
                        if (
                            !minimum_initialized
                            || tail < minimum_tail
                        ) {
                            minimum_tail = tail;
                            minimum_initialized = true;
                        }
                        const long double tolerance =
                            2.0e-12L * std::max(1.0L, scale);
                        if (tail < -tolerance) {
                            std::cout
                                << std::setprecision(20)
                                << "SU2_FINITE_ANCHORED_SPECTRAL_TAILS"
                                << " counterexample"
                                << " level=" << level
                                << " factor=" << factor
                                << " target=" << target
                                << " cutoff_eigenvalue="
                                << groups[index - 1U].eigenvalue
                                << " tail=" << tail
                                << " scale=" << scale
                                << " result=NEGATIVE_UPPER_TAIL"
                                << '\n';
                            return EXIT_SUCCESS;
                        }
                    }
                }
            }
        }

        std::cout
            << std::setprecision(20)
            << "SU2_FINITE_ANCHORED_SPECTRAL_TAILS"
            << " maximum_level=" << maximum_level
            << " parameter_rows=" << parameter_rows
            << " target_rows=" << target_rows
            << " grouped_atoms=" << grouped_atoms
            << " minimum_tail=" << minimum_tail
            << " result=NO_NEGATIVE_UPPER_TAIL"
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
