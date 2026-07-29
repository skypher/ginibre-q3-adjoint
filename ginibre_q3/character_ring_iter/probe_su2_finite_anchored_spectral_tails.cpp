#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using boost::multiprecision::cpp_int;
using IntegerMatrix = std::vector<std::vector<cpp_int>>;

struct Atom {
    long double eigenvalue = 0.0L;
    long double residue = 0.0L;
};

struct ShiftResult {
    std::size_t tested_tails = 0U;
    std::size_t negative_tails = 0U;
    long double minimum_tail = 0.0L;
    bool has_witness = false;
    long double witness_tail = 0.0L;
    int level = 0;
    int factor = 0;
    int target = 0;
    long double cutoff_eigenvalue = 0.0L;
    long double scale = 0.0L;
};

struct ConvexResult {
    std::size_t tested_calls = 0U;
    std::size_t negative_calls = 0U;
    long double minimum_call = 0.0L;
    bool has_witness = false;
    long double witness_call = 0.0L;
    int level = 0;
    int factor = 0;
    int target = 0;
    long double cutoff = 0.0L;
    long double scale = 0.0L;
};

void record_convex_call(
    ConvexResult& result,
    long double call,
    long double scale,
    int level,
    int factor,
    int target,
    long double cutoff
) {
    ++result.tested_calls;
    result.minimum_call = std::min(result.minimum_call, call);
    const long double tolerance =
        2.0e-10L * std::max(1.0L, std::abs(scale));
    if (call < -tolerance) {
        ++result.negative_calls;
        if (!result.has_witness || call < result.witness_call) {
            result.has_witness = true;
            result.witness_call = call;
            result.level = level;
            result.factor = factor;
            result.target = target;
            result.cutoff = cutoff;
            result.scale = scale;
        }
    }
}

int parse_bounded(
    const char* text,
    int minimum,
    int maximum,
    const std::string& name
) {
    const std::string value{text};
    std::size_t consumed = 0U;
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

IntegerMatrix identity_matrix(int dimension) {
    IntegerMatrix result(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(
            static_cast<std::size_t>(dimension),
            cpp_int{0}
        )
    );
    for (int index = 0; index < dimension; ++index) {
        result[static_cast<std::size_t>(index)][
            static_cast<std::size_t>(index)
        ] = 1;
    }
    return result;
}

IntegerMatrix multiply(
    const IntegerMatrix& left,
    const IntegerMatrix& right
) {
    const std::size_t dimension = left.size();
    IntegerMatrix result(
        dimension,
        std::vector<cpp_int>(dimension, cpp_int{0})
    );
    for (std::size_t row = 0U; row < dimension; ++row) {
        for (std::size_t inner = 0U; inner < dimension; ++inner) {
            if (left[row][inner] == 0) {
                continue;
            }
            for (std::size_t column = 0U;
                 column < dimension;
                 ++column) {
                result[row][column] +=
                    left[row][inner] * right[inner][column];
            }
        }
    }
    return result;
}

IntegerMatrix fusion_matrix(int half_level, int factor) {
    const int dimension = half_level + 1;
    IntegerMatrix result(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(
            static_cast<std::size_t>(dimension),
            cpp_int{0}
        )
    );
    for (int row = 0; row < dimension; ++row) {
        const int lower = std::abs(row - factor);
        const int upper = std::min(
            row + factor,
            2 * half_level - row - factor
        );
        for (int column = lower; column <= upper; ++column) {
            result[static_cast<std::size_t>(row)][
                static_cast<std::size_t>(column)
            ] = 1;
        }
    }
    return result;
}

cpp_int anchored_minor(
    const IntegerMatrix& matrix,
    int factor,
    int target
) {
    return
        matrix[0U][0U]
            * matrix[static_cast<std::size_t>(factor)][
                static_cast<std::size_t>(target)
            ]
        - matrix[0U][static_cast<std::size_t>(target)]
            * matrix[static_cast<std::size_t>(factor)][0U];
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

std::vector<long double> mode_vector(
    int level,
    int half_level,
    int mode
) {
    const long double pi = std::acos(-1.0L);
    const long double normalization =
        std::sqrt(2.0L / static_cast<long double>(level + 2));
    const long double pairing_factor =
        mode == half_level ? 1.0L : std::sqrt(2.0L);
    std::vector<long double> result(
        static_cast<std::size_t>(half_level + 1)
    );
    for (int vertex = 0; vertex <= half_level; ++vertex) {
        result[static_cast<std::size_t>(vertex)] =
            pairing_factor
            * normalization
            * std::sin(
                static_cast<long double>(
                    (2 * vertex + 1) * (mode + 1)
                )
                * pi
                / static_cast<long double>(level + 2)
            );
    }
    return result;
}

long double fusion_eigenvalue(
    int level,
    int factor,
    int mode
) {
    const long double pi = std::acos(-1.0L);
    const long double theta =
        static_cast<long double>(mode + 1)
        * pi
        / static_cast<long double>(level + 2);
    return
        std::sin(
            static_cast<long double>(2 * factor + 1) * theta
        )
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

void crosscheck_moment(
    const std::vector<Atom>& groups,
    int shift,
    const cpp_int& exact,
    int level,
    int factor,
    int target
) {
    long double moment = 0.0L;
    long double moment_scale = 0.0L;
    for (const Atom& group : groups) {
        const long double summand =
            group.residue * std::pow(group.eigenvalue, shift);
        moment += summand;
        moment_scale += std::abs(summand);
    }
    const long double expected = exact.convert_to<long double>();
    const long double tolerance =
        2.0e-9L
        * std::max(
            1.0L,
            std::max(std::abs(expected), moment_scale)
        );
    if (std::abs(moment - expected) > tolerance) {
        throw std::runtime_error(
            "spectral moment disagrees with exact fusion minor"
            " at level=" + std::to_string(level)
            + " factor=" + std::to_string(factor)
            + " target=" + std::to_string(target)
            + " shift=" + std::to_string(shift)
            + " moment="
            + std::to_string(static_cast<double>(moment))
            + " exact=" + exact.convert_to<std::string>()
        );
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level =
            argc >= 2
                ? parse_bounded(argv[1], 4, 400, "maximum_level")
                : 80;
        const int maximum_shift =
            argc >= 3
                ? parse_bounded(argv[2], 0, 40, "maximum_shift")
                : 12;
        const int minimum_level =
            argc >= 4
                ? parse_bounded(
                    argv[3],
                    4,
                    maximum_level,
                    "minimum_level"
                )
                : 4;
        if (argc > 4) {
            throw std::runtime_error(
                "usage: probe_su2_finite_anchored_spectral_tails "
                "[maximum_level [maximum_shift [minimum_level]]]"
            );
        }
        const int first_level =
            (minimum_level % 2 == 0)
                ? minimum_level
                : minimum_level + 1;

        std::vector<ShiftResult> results(
            static_cast<std::size_t>(maximum_shift + 1)
        );
        std::vector<ShiftResult> affine_results(
            static_cast<std::size_t>(maximum_shift + 1)
        );
        std::vector<ConvexResult> convex_results(
            static_cast<std::size_t>(maximum_shift + 1)
        );
        std::vector<ConvexResult> affine_convex_results(
            static_cast<std::size_t>(maximum_shift + 1)
        );
        std::vector<ConvexResult> power_convex_results(
            static_cast<std::size_t>(maximum_shift + 1)
        );
        std::vector<ConvexResult> affine_power_convex_results(
            static_cast<std::size_t>(maximum_shift + 1)
        );
        ConvexResult closure_anchor_interior_log_convex;
        ConvexResult closure_anchor_interior_power_convex;
        std::size_t parameter_rows = 0U;
        std::size_t target_rows = 0U;
        std::size_t grouped_atoms = 0U;
        std::size_t exact_moment_crosschecks = 0U;

        for (
            int level = first_level;
            level <= maximum_level;
            level += 2
        ) {
            const int half_level = level / 2;
            std::vector<std::vector<long double>> modes;
            modes.reserve(static_cast<std::size_t>(half_level + 1));
            for (int mode = 0; mode <= half_level; ++mode) {
                modes.push_back(
                    mode_vector(level, half_level, mode)
                );
            }

            for (
                int factor = 1;
                2 * factor < half_level;
                ++factor
            ) {
                ++parameter_rows;
                std::vector<long double> eigenvalues(
                    static_cast<std::size_t>(half_level + 1)
                );
                for (int mode = 0; mode <= half_level; ++mode) {
                    eigenvalues[static_cast<std::size_t>(mode)] =
                        fusion_eigenvalue(level, factor, mode);
                }

                const IntegerMatrix fusion =
                    fusion_matrix(half_level, factor);
                const IntegerMatrix square =
                    multiply(fusion, fusion);
                std::vector<IntegerMatrix> powers;
                powers.reserve(
                    static_cast<std::size_t>(maximum_shift + 1)
                );
                powers.push_back(identity_matrix(half_level + 1));
                for (int shift = 1;
                     shift <= maximum_shift;
                     ++shift) {
                    powers.push_back(
                        multiply(powers.back(), square)
                    );
                }

                for (int target = 1;
                     target <= half_level;
                     ++target) {
                    ++target_rows;
                    std::vector<Atom> atoms;
                    for (int left = 0;
                         left <= half_level;
                         ++left) {
                        for (int right = left + 1;
                             right <= half_level;
                             ++right) {
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
                            return
                                left.eigenvalue < right.eigenvalue;
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

                    for (int shift = 0;
                         shift <= maximum_shift;
                         ++shift) {
                        const cpp_int exact = anchored_minor(
                            powers[static_cast<std::size_t>(shift)],
                            factor,
                            target
                        );
                        crosscheck_moment(
                            groups,
                            shift,
                            exact,
                            level,
                            factor,
                            target
                        );
                        ++exact_moment_crosschecks;

                        long double tail = 0.0L;
                        long double scale = 0.0L;
                        ShiftResult& result =
                            results[static_cast<std::size_t>(shift)];
                        ShiftResult& affine_result =
                            affine_results[
                                static_cast<std::size_t>(shift)
                            ];
                        for (
                            std::size_t index = groups.size();
                            index > 0U;
                            --index
                        ) {
                            const long double weighted_residue =
                                groups[index - 1U].residue
                                * std::pow(
                                    groups[index - 1U].eigenvalue,
                                    shift
                                );
                            tail += weighted_residue;
                            scale += std::abs(weighted_residue);
                            ++result.tested_tails;
                            result.minimum_tail =
                                std::min(result.minimum_tail, tail);
                            const bool affine_relevant =
                                2 * shift * factor > half_level;
                            if (affine_relevant) {
                                ++affine_result.tested_tails;
                                affine_result.minimum_tail = std::min(
                                    affine_result.minimum_tail,
                                    tail
                                );
                            }
                            const long double tolerance =
                                2.0e-10L
                                * std::max(1.0L, scale);
                            if (tail < -tolerance) {
                                ++result.negative_tails;
                                if (
                                    !result.has_witness
                                    || tail < result.witness_tail
                                ) {
                                    result.has_witness = true;
                                    result.witness_tail = tail;
                                    result.level = level;
                                    result.factor = factor;
                                    result.target = target;
                                    result.cutoff_eigenvalue =
                                        groups[index - 1U].eigenvalue;
                                    result.scale = scale;
                                }
                                if (affine_relevant) {
                                    ++affine_result.negative_tails;
                                    if (
                                        !affine_result.has_witness
                                        || tail
                                            < affine_result.witness_tail
                                    ) {
                                        affine_result.has_witness = true;
                                        affine_result.witness_tail = tail;
                                        affine_result.level = level;
                                        affine_result.factor = factor;
                                        affine_result.target = target;
                                        affine_result.cutoff_eigenvalue =
                                            groups[
                                                index - 1U
                                            ].eigenvalue;
                                        affine_result.scale = scale;
                                    }
                                }
                            }
                        }

                        if (shift > 0) {
                            const std::size_t count = groups.size();
                            std::vector<long double> suffix_weight(
                                count + 1U,
                                0.0L
                            );
                            std::vector<long double> suffix_weight_log(
                                count + 1U,
                                0.0L
                            );
                            std::vector<long double> suffix_absolute(
                                count + 1U,
                                0.0L
                            );
                            std::vector<long double>
                                suffix_absolute_log(
                                    count + 1U,
                                    0.0L
                                );
                            std::vector<long double> suffix_weight_base(
                                count + 1U,
                                0.0L
                            );
                            std::vector<long double>
                                suffix_absolute_base(
                                    count + 1U,
                                    0.0L
                                );
                            for (std::size_t index = count;
                                 index > 0U;
                                 --index) {
                                suffix_weight[index - 1U] =
                                    suffix_weight[index];
                                suffix_weight_log[index - 1U] =
                                    suffix_weight_log[index];
                                suffix_absolute[index - 1U] =
                                    suffix_absolute[index];
                                suffix_absolute_log[index - 1U] =
                                    suffix_absolute_log[index];
                                suffix_weight_base[index - 1U] =
                                    suffix_weight_base[index];
                                suffix_absolute_base[index - 1U] =
                                    suffix_absolute_base[index];
                                const long double eigenvalue =
                                    groups[index - 1U].eigenvalue;
                                if (eigenvalue <= 1.0e-24L) {
                                    continue;
                                }
                                const long double logarithm =
                                    std::log(eigenvalue);
                                const long double weight =
                                    groups[index - 1U].residue
                                    * std::pow(eigenvalue, shift);
                                suffix_weight[index - 1U] += weight;
                                suffix_weight_log[index - 1U] +=
                                    weight * logarithm;
                                suffix_absolute[index - 1U] +=
                                    std::abs(weight);
                                suffix_absolute_log[index - 1U] +=
                                    std::abs(weight) * logarithm;
                                suffix_weight_base[index - 1U] +=
                                    weight * eigenvalue;
                                suffix_absolute_base[index - 1U] +=
                                    std::abs(weight) * eigenvalue;
                            }

                            ConvexResult& convex =
                                convex_results[
                                    static_cast<std::size_t>(shift)
                                ];
                            ConvexResult& affine_convex =
                                affine_convex_results[
                                    static_cast<std::size_t>(shift)
                                ];
                            const bool affine_relevant =
                                2 * shift * factor > half_level
                                && target < half_level;
                            for (std::size_t index = 0U;
                                 index < count;
                                 ++index) {
                                const long double eigenvalue =
                                    groups[index].eigenvalue;
                                if (eigenvalue <= 1.0e-24L) {
                                    continue;
                                }
                                const long double cutoff =
                                    std::log(eigenvalue);
                                const long double call =
                                    suffix_weight_log[index]
                                    - cutoff * suffix_weight[index];
                                const long double call_scale =
                                    suffix_absolute_log[index]
                                    - cutoff * suffix_absolute[index];
                                ++convex.tested_calls;
                                convex.minimum_call =
                                    std::min(
                                        convex.minimum_call,
                                        call
                                    );
                                if (affine_relevant) {
                                    ++affine_convex.tested_calls;
                                    affine_convex.minimum_call =
                                        std::min(
                                            affine_convex.minimum_call,
                                            call
                                        );
                                }
                                const long double tolerance =
                                    2.0e-10L
                                    * std::max(
                                        1.0L,
                                        std::abs(call_scale)
                                    );
                                if (call < -tolerance) {
                                    ++convex.negative_calls;
                                    if (
                                        !convex.has_witness
                                        || call < convex.witness_call
                                    ) {
                                        convex.has_witness = true;
                                        convex.witness_call = call;
                                        convex.level = level;
                                        convex.factor = factor;
                                        convex.target = target;
                                        convex.cutoff = cutoff;
                                        convex.scale = call_scale;
                                    }
                                    if (affine_relevant) {
                                        ++affine_convex.negative_calls;
                                        if (
                                            !affine_convex.has_witness
                                            || call
                                                < affine_convex
                                                    .witness_call
                                        ) {
                                            affine_convex.has_witness =
                                                true;
                                            affine_convex.witness_call =
                                                call;
                                            affine_convex.level = level;
                                            affine_convex.factor =
                                                factor;
                                            affine_convex.target =
                                                target;
                                            affine_convex.cutoff =
                                                cutoff;
                                            affine_convex.scale =
                                                call_scale;
                                        }
                                    }
                                }
                                const int closure_anchor = std::max(
                                    4,
                                    half_level / (2 * factor)
                                );
                                if (
                                    shift == closure_anchor
                                    && target < half_level
                                ) {
                                    record_convex_call(
                                        closure_anchor_interior_log_convex,
                                        call,
                                        call_scale,
                                        level,
                                        factor,
                                        target,
                                        cutoff
                                    );
                                }

                                const long double power_call =
                                    suffix_weight_base[index]
                                    - eigenvalue
                                        * suffix_weight[index];
                                const long double power_scale =
                                    suffix_absolute_base[index]
                                    - eigenvalue
                                        * suffix_absolute[index];
                                ConvexResult& power_convex =
                                    power_convex_results[
                                        static_cast<std::size_t>(shift)
                                    ];
                                record_convex_call(
                                    power_convex,
                                    power_call,
                                    power_scale,
                                    level,
                                    factor,
                                    target,
                                    eigenvalue
                                );
                                if (affine_relevant) {
                                    ConvexResult&
                                        affine_power_convex =
                                            affine_power_convex_results[
                                                static_cast<std::size_t>(
                                                    shift
                                                )
                                            ];
                                    record_convex_call(
                                        affine_power_convex,
                                        power_call,
                                        power_scale,
                                        level,
                                        factor,
                                        target,
                                        eigenvalue
                                    );
                                }
                                if (
                                    shift == closure_anchor
                                    && target < half_level
                                ) {
                                    record_convex_call(
                                        closure_anchor_interior_power_convex,
                                        power_call,
                                        power_scale,
                                        level,
                                        factor,
                                        target,
                                        eigenvalue
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << std::setprecision(20)
            << "SU2_FINITE_ANCHORED_SPECTRAL_TAILS"
            << " minimum_level=" << first_level
            << " maximum_level=" << maximum_level
            << " maximum_shift=" << maximum_shift
            << " parameter_rows=" << parameter_rows
            << " target_rows=" << target_rows
            << " grouped_atoms=" << grouped_atoms
            << " exact_moment_crosschecks="
            << exact_moment_crosschecks
            << '\n';
        for (int shift = 0; shift <= maximum_shift; ++shift) {
            const ShiftResult& result =
                results[static_cast<std::size_t>(shift)];
            std::cout
                << std::setprecision(20)
                << "SU2_FINITE_ANCHORED_SPECTRAL_TAIL_SHIFT"
                << " shift=" << shift
                << " tested_tails=" << result.tested_tails
                << " negative_tails=" << result.negative_tails
                << " minimum_tail=" << result.minimum_tail;
            if (result.has_witness) {
                std::cout
                    << " worst_level=" << result.level
                    << " worst_factor=" << result.factor
                    << " worst_target=" << result.target
                    << " worst_tail=" << result.witness_tail
                    << " worst_cutoff_eigenvalue="
                    << result.cutoff_eigenvalue
                    << " worst_scale=" << result.scale;
            }
            std::cout
                << " result="
                << (
                    result.negative_tails == 0U
                        ? "NO_NEGATIVE_UPPER_TAIL"
                        : "NEGATIVE_UPPER_TAIL"
                )
                << '\n';
            const ShiftResult& affine_result =
                affine_results[static_cast<std::size_t>(shift)];
            std::cout
                << std::setprecision(20)
                << "SU2_FINITE_ANCHORED_AFFINE_TAIL_SHIFT"
                << " shift=" << shift
                << " tested_tails=" << affine_result.tested_tails
                << " negative_tails="
                << affine_result.negative_tails
                << " minimum_tail="
                << affine_result.minimum_tail;
            if (affine_result.has_witness) {
                std::cout
                    << " worst_level=" << affine_result.level
                    << " worst_factor=" << affine_result.factor
                    << " worst_target=" << affine_result.target
                    << " worst_tail="
                    << affine_result.witness_tail
                    << " worst_cutoff_eigenvalue="
                    << affine_result.cutoff_eigenvalue
                    << " worst_scale=" << affine_result.scale;
            }
            std::cout
                << " result="
                << (
                    affine_result.negative_tails == 0U
                        ? "NO_NEGATIVE_UPPER_TAIL"
                        : "NEGATIVE_UPPER_TAIL"
                )
                << '\n';
            const ConvexResult& convex =
                convex_results[static_cast<std::size_t>(shift)];
            std::cout
                << std::setprecision(20)
                << "SU2_FINITE_ANCHORED_LOG_CONVEX_SHIFT"
                << " shift=" << shift
                << " tested_calls=" << convex.tested_calls
                << " negative_calls=" << convex.negative_calls
                << " minimum_call=" << convex.minimum_call;
            if (convex.has_witness) {
                std::cout
                    << " worst_level=" << convex.level
                    << " worst_factor=" << convex.factor
                    << " worst_target=" << convex.target
                    << " worst_call=" << convex.witness_call
                    << " worst_cutoff_log_eigenvalue="
                    << convex.cutoff
                    << " worst_scale=" << convex.scale;
            }
            std::cout
                << " result="
                << (
                    convex.negative_calls == 0U
                        ? "NO_NEGATIVE_CALL"
                        : "NEGATIVE_CALL"
                )
                << '\n';
            const ConvexResult& affine_convex =
                affine_convex_results[
                    static_cast<std::size_t>(shift)
                ];
            std::cout
                << std::setprecision(20)
                << "SU2_FINITE_ANCHORED_AFFINE_INTERIOR_LOG_CONVEX_SHIFT"
                << " shift=" << shift
                << " tested_calls="
                << affine_convex.tested_calls
                << " negative_calls="
                << affine_convex.negative_calls
                << " minimum_call="
                << affine_convex.minimum_call;
            if (affine_convex.has_witness) {
                std::cout
                    << " worst_level=" << affine_convex.level
                    << " worst_factor=" << affine_convex.factor
                    << " worst_target=" << affine_convex.target
                    << " worst_call="
                    << affine_convex.witness_call
                    << " worst_cutoff_log_eigenvalue="
                    << affine_convex.cutoff
                    << " worst_scale=" << affine_convex.scale;
            }
            std::cout
                << " result="
                << (
                    affine_convex.negative_calls == 0U
                        ? "NO_NEGATIVE_CALL"
                        : "NEGATIVE_CALL"
                )
                << '\n';
            const ConvexResult& power_convex =
                power_convex_results[
                    static_cast<std::size_t>(shift)
                ];
            std::cout
                << std::setprecision(20)
                << "SU2_FINITE_ANCHORED_POWER_CONVEX_SHIFT"
                << " shift=" << shift
                << " tested_calls="
                << power_convex.tested_calls
                << " negative_calls="
                << power_convex.negative_calls
                << " minimum_call="
                << power_convex.minimum_call;
            if (power_convex.has_witness) {
                std::cout
                    << " worst_level=" << power_convex.level
                    << " worst_factor=" << power_convex.factor
                    << " worst_target=" << power_convex.target
                    << " worst_call="
                    << power_convex.witness_call
                    << " worst_cutoff_eigenvalue="
                    << power_convex.cutoff
                    << " worst_scale=" << power_convex.scale;
            }
            std::cout
                << " result="
                << (
                    power_convex.negative_calls == 0U
                        ? "NO_NEGATIVE_CALL"
                        : "NEGATIVE_CALL"
                )
                << '\n';
            const ConvexResult& affine_power_convex =
                affine_power_convex_results[
                    static_cast<std::size_t>(shift)
                ];
            std::cout
                << std::setprecision(20)
                << "SU2_FINITE_ANCHORED_AFFINE_INTERIOR_POWER_CONVEX_SHIFT"
                << " shift=" << shift
                << " tested_calls="
                << affine_power_convex.tested_calls
                << " negative_calls="
                << affine_power_convex.negative_calls
                << " minimum_call="
                << affine_power_convex.minimum_call;
            if (affine_power_convex.has_witness) {
                std::cout
                    << " worst_level="
                    << affine_power_convex.level
                    << " worst_factor="
                    << affine_power_convex.factor
                    << " worst_target="
                    << affine_power_convex.target
                    << " worst_call="
                    << affine_power_convex.witness_call
                    << " worst_cutoff_eigenvalue="
                    << affine_power_convex.cutoff
                    << " worst_scale="
                    << affine_power_convex.scale;
            }
            std::cout
                << " result="
                << (
                    affine_power_convex.negative_calls == 0U
                        ? "NO_NEGATIVE_CALL"
                        : "NEGATIVE_CALL"
                )
                << '\n';
        }
        std::cout
            << std::setprecision(20)
            << "SU2_FINITE_ANCHORED_LOG_CLOSURE_ANCHOR"
            << " tested_calls="
            << closure_anchor_interior_log_convex.tested_calls
            << " negative_calls="
            << closure_anchor_interior_log_convex.negative_calls
            << " minimum_call="
            << closure_anchor_interior_log_convex.minimum_call;
        if (closure_anchor_interior_log_convex.has_witness) {
            std::cout
                << " worst_level="
                << closure_anchor_interior_log_convex.level
                << " worst_factor="
                << closure_anchor_interior_log_convex.factor
                << " worst_target="
                << closure_anchor_interior_log_convex.target
                << " worst_call="
                << closure_anchor_interior_log_convex.witness_call
                << " worst_cutoff_log_eigenvalue="
                << closure_anchor_interior_log_convex.cutoff
                << " worst_scale="
                << closure_anchor_interior_log_convex.scale;
        }
        std::cout
            << " result="
            << (
                closure_anchor_interior_log_convex.negative_calls == 0U
                    ? "NO_NEGATIVE_CALL"
                    : "NEGATIVE_CALL"
            )
            << '\n';
        std::cout
            << std::setprecision(20)
            << "SU2_FINITE_ANCHORED_CLOSURE_ANCHOR"
            << " tested_calls="
            << closure_anchor_interior_power_convex.tested_calls
            << " negative_calls="
            << closure_anchor_interior_power_convex.negative_calls
            << " minimum_call="
            << closure_anchor_interior_power_convex.minimum_call;
        if (closure_anchor_interior_power_convex.has_witness) {
            std::cout
                << " worst_level="
                << closure_anchor_interior_power_convex.level
                << " worst_factor="
                << closure_anchor_interior_power_convex.factor
                << " worst_target="
                << closure_anchor_interior_power_convex.target
                << " worst_call="
                << closure_anchor_interior_power_convex.witness_call
                << " worst_cutoff_eigenvalue="
                << closure_anchor_interior_power_convex.cutoff
                << " worst_scale="
                << closure_anchor_interior_power_convex.scale;
        }
        std::cout
            << " result="
            << (
                closure_anchor_interior_power_convex.negative_calls == 0U
                    ? "NO_NEGATIVE_CALL"
                    : "NEGATIVE_CALL"
            )
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
