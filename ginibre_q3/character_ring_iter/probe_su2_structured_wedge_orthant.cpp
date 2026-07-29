#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value <= 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error(
            std::string(name) + " must be a positive integer"
        );
    }
    return static_cast<int>(value);
}

class SymmetricProfile {
public:
    SymmetricProfile(int q, int power)
        : radius_(q * power),
          values_(static_cast<std::size_t>(2 * radius_ + 1)) {
        std::vector<Integer> current{Integer(1)};
        int current_radius = 0;
        for (int step = 0; step < power; ++step) {
            const int next_radius = current_radius + q;
            std::vector<Integer> next(
                static_cast<std::size_t>(2 * next_radius + 1)
            );
            for (int exponent = -current_radius;
                 exponent <= current_radius;
                 ++exponent) {
                const Integer& coefficient = current[
                    static_cast<std::size_t>(
                        exponent + current_radius
                    )
                ];
                for (int increment = -q; increment <= q; ++increment) {
                    next[static_cast<std::size_t>(
                        exponent + increment + next_radius
                    )] += coefficient;
                }
            }
            current = std::move(next);
            current_radius = next_radius;
        }
        values_ = std::move(current);
    }

    Integer at(int exponent) const {
        if (exponent < -radius_ || exponent > radius_) {
            return 0;
        }
        return values_[static_cast<std::size_t>(exponent + radius_)];
    }

    int radius() const {
        return radius_;
    }

private:
    int radius_;
    std::vector<Integer> values_;
};

Integer psi(const SymmetricProfile& profile, int row, int index) {
    return profile.at(index - row) - profile.at(index + row + 1);
}

Integer unfolded_defect(
    const SymmetricProfile& profile,
    int row,
    int left,
    int gap
) {
    const Integer left_difference =
        profile.at(left) - profile.at(left + 1);
    const Integer right_difference =
        profile.at(left + gap) - profile.at(left + gap + 1);
    return
        left_difference * profile.at(left + gap - row)
        - right_difference * profile.at(left - row);
}

Integer current(
    const std::vector<Integer>& base,
    const std::vector<Integer>& row,
    int left,
    int right
) {
    return
        base[static_cast<std::size_t>(left)]
            * row[static_cast<std::size_t>(right)]
        - base[static_cast<std::size_t>(right)]
            * row[static_cast<std::size_t>(left)];
}

Integer absolute_value(const Integer& value) {
    return value < 0 ? -value : value;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = argc >= 4 ? argv[3] : "";
        const bool single_case = mode == "single" || mode == "tp2";
        const bool tp2_range = mode == "tp2-range" && argc == 5;
        const bool tp2_only =
            mode == "tp2" || mode == "tp2-all" || tp2_range;
        if (
            argc != 3
            && !(argc == 4 && (single_case || tp2_only))
            && !tp2_range
        ) {
            throw std::runtime_error(
                "usage: probe_su2_structured_wedge_orthant "
                "maximum_half_label maximum_half_power "
                "[single|tp2|tp2-all|tp2-range minimum_half_power]"
            );
        }
        const int maximum_q =
            parse_positive(argv[1], "maximum_half_label");
        const int maximum_power =
            parse_positive(argv[2], "maximum_half_power");
        const int minimum_power =
            tp2_range
            ? parse_positive(argv[4], "minimum_half_power")
            : 1;
        if (minimum_power > maximum_power) {
            throw std::runtime_error(
                "minimum_half_power exceeds maximum_half_power"
            );
        }

        std::size_t cases = 0;
        std::size_t kernel_profiles = 0;
        std::size_t orthants = 0;
        std::size_t negative_fixed_gap_suffixes = 0;
        std::size_t negative_fixed_base_suffixes = 0;
        std::size_t negative_direct_suffixes = 0;
        std::size_t positive_cross_suffixes = 0;
        std::size_t target_label_decreases = 0;
        std::size_t boundary_reserve_failures = 0;
        int maximum_product_negative_blocks = 0;
        int maximum_current_negative_blocks = 0;
        std::size_t negative_endpoint_mismatches = 0;
        std::size_t nonnested_negative_intervals = 0;
        std::size_t layer_payment_failures = 0;
        std::size_t negative_orientation_failures = 0;
        std::size_t loss_magnitude_decreases = 0;
        std::size_t post_block_payment_failures = 0;
        std::size_t first_post_payment_failures = 0;
        std::size_t negative_full_adjacent_minors = 0;
        int maximum_post_terms_needed = 0;
        bool have_minimum_active_adjacent_minor = false;
        Integer minimum_active_adjacent_minor = 0;
        bool have_first_fixed_gap_negative = false;
        bool have_first_fixed_base_negative = false;
        bool have_first_direct_negative = false;
        bool have_first_cross_positive = false;
        bool have_first_target_label_decrease = false;
        bool have_first_boundary_reserve_failure = false;
        std::string first_fixed_gap_negative;
        std::string first_fixed_base_negative;
        std::string first_direct_negative;
        std::string first_cross_positive;
        std::string first_target_label_decrease;
        std::string first_boundary_reserve_failure;
        std::string first_layer_payment_failure;
        std::string first_post_block_payment_failure;
        std::string first_first_post_payment_failure;
        std::string maximum_post_terms_case;
        std::string first_negative_full_adjacent_minor;
        std::string minimum_active_adjacent_minor_case;
        for (int q = single_case ? maximum_q : 1;
             q <= maximum_q;
             ++q) {
            for (
                int power =
                    single_case ? maximum_power : minimum_power;
                 power <= maximum_power;
                 ++power) {
                const SymmetricProfile profile(q, power);
                const int maximum_target = 2 * q * power;
                const int maximum_index =
                    profile.radius() + maximum_target;
                const int width = maximum_index + 2;
                std::vector<Integer> base(
                    static_cast<std::size_t>(width)
                );
                std::vector<Integer> distinguished(
                    static_cast<std::size_t>(width)
                );
                for (int index = 0; index <= maximum_index; ++index) {
                    base[static_cast<std::size_t>(index)] =
                        psi(profile, 0, index);
                    distinguished[static_cast<std::size_t>(index)] =
                        psi(profile, q, index);
                }
                for (int row = 0; row < maximum_target; ++row) {
                    for (int column = 0;
                         column < maximum_index;
                         ++column) {
                        const Integer minor =
                            psi(profile, row, column)
                                * psi(
                                    profile,
                                    row + 1,
                                    column + 1
                                )
                            - psi(profile, row, column + 1)
                                * psi(profile, row + 1, column);
                        const bool active =
                            psi(profile, row, column) > 0
                            && psi(profile, row, column + 1) > 0
                            && psi(profile, row + 1, column) > 0
                            && psi(
                                profile,
                                row + 1,
                                column + 1
                            ) > 0;
                        if (
                            active
                            && (
                                !have_minimum_active_adjacent_minor
                                || minor
                                    < minimum_active_adjacent_minor
                            )
                        ) {
                            have_minimum_active_adjacent_minor = true;
                            minimum_active_adjacent_minor = minor;
                            minimum_active_adjacent_minor_case =
                                "{q=" + std::to_string(q)
                                + " half_power="
                                + std::to_string(power)
                                + " row="
                                + std::to_string(row)
                                + " column="
                                + std::to_string(column)
                                + " value=" + minor.str()
                                + "}";
                        }
                        if (minor < 0) {
                            ++negative_full_adjacent_minors;
                            if (
                                first_negative_full_adjacent_minor
                                    .empty()
                            ) {
                                first_negative_full_adjacent_minor =
                                    "{q=" + std::to_string(q)
                                    + " half_power="
                                    + std::to_string(power)
                                    + " row="
                                    + std::to_string(row)
                                    + " column="
                                    + std::to_string(column)
                                    + " value=" + minor.str()
                                    + "}";
                            }
                        }
                    }
                }
                ++kernel_profiles;
                if (tp2_only) {
                    continue;
                }
                std::vector<std::vector<Integer>>
                    previous_fixed_gap_suffix(
                        static_cast<std::size_t>(width),
                        std::vector<Integer>(
                            static_cast<std::size_t>(width)
                        )
                    );

                for (int target = 1;
                     target <= maximum_target;
                     ++target) {
                    std::vector<Integer> target_row(
                        static_cast<std::size_t>(width)
                    );
                    for (int index = 0;
                         index <= maximum_index;
                         ++index) {
                        target_row[static_cast<std::size_t>(index)] =
                            psi(profile, target, index);
                    }
                    std::vector<std::vector<Integer>> suffix(
                        static_cast<std::size_t>(width),
                        std::vector<Integer>(
                            static_cast<std::size_t>(width)
                        )
                    );
                    for (int gap = 1;
                         gap <= maximum_index;
                         ++gap) {
                        Integer direct_suffix = 0;
                        Integer cross_suffix = 0;
                        Integer coupled_suffix = 0;
                        Integer negative_mass = 0;
                        Integer boundary_reserve = 0;
                        int product_negative_blocks = 0;
                        bool inside_product_negative_block = false;
                        int distinguished_negative_blocks = 0;
                        int target_negative_blocks = 0;
                        bool inside_distinguished_negative = false;
                        bool inside_target_negative = false;
                        int distinguished_outer_negative = -1;
                        int target_outer_negative = -1;
                        int distinguished_inner_negative =
                            maximum_index + 1;
                        int target_inner_negative = maximum_index + 1;
                        std::vector<Integer> distinguished_currents(
                            static_cast<std::size_t>(width)
                        );
                        std::vector<Integer> target_currents(
                            static_cast<std::size_t>(width)
                        );
                        std::vector<Integer> products(
                            static_cast<std::size_t>(width)
                        );
                        for (int left = maximum_index - gap;
                             left >= 0;
                             --left) {
                            const int reflected_left =
                                -left - gap - 1;
                            const Integer distinguished_left =
                                unfolded_defect(
                                    profile,
                                    q,
                                    left,
                                    gap
                                );
                            const Integer distinguished_reflected =
                                unfolded_defect(
                                    profile,
                                    q,
                                    reflected_left,
                                    gap
                                );
                            const Integer target_left =
                                unfolded_defect(
                                    profile,
                                    target,
                                    left,
                                    gap
                                );
                            const Integer target_reflected =
                                unfolded_defect(
                                    profile,
                                    target,
                                    reflected_left,
                                    gap
                                );
                            direct_suffix +=
                                distinguished_left * target_left
                                + distinguished_reflected
                                    * target_reflected;
                            cross_suffix +=
                                distinguished_left * target_reflected
                                + distinguished_reflected * target_left;
                            coupled_suffix +=
                                current(
                                    base,
                                    distinguished,
                                    left,
                                    left + gap
                                )
                                * current(
                                    base,
                                    target_row,
                                    left,
                                    left + gap
                                );
                            const Integer product =
                                current(
                                    base,
                                    distinguished,
                                    left,
                                    left + gap
                                )
                                * current(
                                    base,
                                    target_row,
                                    left,
                                    left + gap
                                );
                            const Integer distinguished_current =
                                current(
                                    base,
                                    distinguished,
                                    left,
                                    left + gap
                                );
                            const Integer target_current =
                                current(
                                    base,
                                    target_row,
                                    left,
                                    left + gap
                                );
                            distinguished_currents[
                                static_cast<std::size_t>(left)
                            ] = distinguished_current;
                            target_currents[
                                static_cast<std::size_t>(left)
                            ] = target_current;
                            products[static_cast<std::size_t>(left)] =
                                product;
                            if (distinguished_current < 0) {
                                distinguished_inner_negative = left;
                                if (distinguished_outer_negative < 0) {
                                    distinguished_outer_negative = left;
                                }
                                if (!inside_distinguished_negative) {
                                    ++distinguished_negative_blocks;
                                    inside_distinguished_negative = true;
                                }
                            } else if (distinguished_current > 0) {
                                inside_distinguished_negative = false;
                            }
                            if (target_current < 0) {
                                target_inner_negative = left;
                                if (target_outer_negative < 0) {
                                    target_outer_negative = left;
                                }
                                if (!inside_target_negative) {
                                    ++target_negative_blocks;
                                    inside_target_negative = true;
                                }
                            } else if (target_current > 0) {
                                inside_target_negative = false;
                            }
                            if (product < 0) {
                                negative_mass -= product;
                                if (!inside_product_negative_block) {
                                    ++product_negative_blocks;
                                    inside_product_negative_block = true;
                                }
                            } else if (product > 0) {
                                inside_product_negative_block = false;
                            }
                            if (
                                base[
                                    static_cast<std::size_t>(
                                        left + gap
                                    )
                                ] == 0
                            ) {
                                boundary_reserve += product;
                            }
                            if (
                                coupled_suffix
                                != direct_suffix - cross_suffix
                            ) {
                                throw std::runtime_error(
                                    "Toeplitz-image decomposition mismatch"
                                );
                            }
                            if (direct_suffix < 0) {
                                ++negative_direct_suffixes;
                                if (!have_first_direct_negative) {
                                    have_first_direct_negative = true;
                                    first_direct_negative =
                                        "{q=" + std::to_string(q)
                                        + " half_power="
                                        + std::to_string(power)
                                        + " target="
                                        + std::to_string(target)
                                        + " cutoff="
                                        + std::to_string(left)
                                        + " gap="
                                        + std::to_string(gap)
                                        + " direct="
                                        + direct_suffix.str()
                                        + " cross="
                                        + cross_suffix.str()
                                        + " coupled="
                                        + coupled_suffix.str()
                                        + "}";
                                }
                            }
                            if (cross_suffix > 0) {
                                ++positive_cross_suffixes;
                                if (!have_first_cross_positive) {
                                    have_first_cross_positive = true;
                                    first_cross_positive =
                                        "{q=" + std::to_string(q)
                                        + " half_power="
                                        + std::to_string(power)
                                        + " target="
                                        + std::to_string(target)
                                        + " cutoff="
                                        + std::to_string(left)
                                        + " gap="
                                        + std::to_string(gap)
                                        + " direct="
                                        + direct_suffix.str()
                                        + " cross="
                                        + cross_suffix.str()
                                        + " coupled="
                                        + coupled_suffix.str()
                                    + "}";
                                }
                            }
                        }
                        maximum_product_negative_blocks = std::max(
                            maximum_product_negative_blocks,
                            product_negative_blocks
                        );
                        maximum_current_negative_blocks = std::max(
                            maximum_current_negative_blocks,
                            std::max(
                                distinguished_negative_blocks,
                                target_negative_blocks
                            )
                        );
                        if (
                            distinguished_outer_negative >= 0
                            && target_outer_negative >= 0
                            && distinguished_outer_negative
                                != target_outer_negative
                        ) {
                            ++negative_endpoint_mismatches;
                        }
                        if (
                            distinguished_outer_negative >= 0
                            && target_outer_negative >= 0
                        ) {
                            const bool distinguished_contains_target =
                                distinguished_inner_negative
                                    <= target_inner_negative
                                && distinguished_outer_negative
                                    >= target_outer_negative;
                            const bool target_contains_distinguished =
                                target_inner_negative
                                    <= distinguished_inner_negative
                                && target_outer_negative
                                    >= distinguished_outer_negative;
                            if (
                                !distinguished_contains_target
                                && !target_contains_distinguished
                            ) {
                                ++nonnested_negative_intervals;
                            }
                        }
                        int negative_inner = -1;
                        int negative_outer = -1;
                        int losing_row = 0;
                        for (int left = 0;
                             left + gap <= maximum_index;
                             ++left) {
                            if (
                                products[
                                    static_cast<std::size_t>(left)
                                ] >= 0
                            ) {
                                continue;
                            }
                            if (negative_inner < 0) {
                                negative_inner = left;
                            }
                            negative_outer = left;
                            const int current_losing_row =
                                distinguished_currents[
                                    static_cast<std::size_t>(left)
                                ] < 0
                                ? 1
                                : 2;
                            if (losing_row == 0) {
                                losing_row = current_losing_row;
                            } else if (
                                losing_row != current_losing_row
                            ) {
                                ++negative_orientation_failures;
                            }
                        }
                        if (negative_inner >= 0) {
                            Integer post_block_reserve = 0;
                            Integer first_post_reserve = 0;
                            Integer accumulated_post_reserve = 0;
                            int post_terms_needed = 0;
                            bool reached_negative_mass = false;
                            for (int left = negative_outer + 1;
                                 left + gap <= maximum_index;
                                 ++left) {
                                const Integer& product =
                                    products[
                                        static_cast<std::size_t>(left)
                                    ];
                                if (product <= 0) {
                                    continue;
                                }
                                post_block_reserve += product;
                                if (first_post_reserve == 0) {
                                    first_post_reserve = product;
                                }
                                if (!reached_negative_mass) {
                                    accumulated_post_reserve += product;
                                    ++post_terms_needed;
                                    if (
                                        accumulated_post_reserve
                                            >= negative_mass
                                    ) {
                                        reached_negative_mass = true;
                                    }
                                }
                            }
                            if (
                                post_terms_needed
                                    > maximum_post_terms_needed
                            ) {
                                maximum_post_terms_needed =
                                    post_terms_needed;
                                maximum_post_terms_case =
                                    "{q=" + std::to_string(q)
                                    + " half_power="
                                    + std::to_string(power)
                                    + " target="
                                    + std::to_string(target)
                                    + " gap="
                                    + std::to_string(gap)
                                    + " negative_mass="
                                    + negative_mass.str()
                                    + "}";
                            }
                            if (post_block_reserve < negative_mass) {
                                ++post_block_payment_failures;
                                if (
                                    first_post_block_payment_failure
                                        .empty()
                                ) {
                                    first_post_block_payment_failure =
                                        "{q=" + std::to_string(q)
                                        + " half_power="
                                        + std::to_string(power)
                                        + " target="
                                        + std::to_string(target)
                                        + " gap="
                                        + std::to_string(gap)
                                        + " negative_mass="
                                        + negative_mass.str()
                                        + " reserve="
                                        + post_block_reserve.str()
                                        + "}";
                                }
                            }
                            if (first_post_reserve < negative_mass) {
                                ++first_post_payment_failures;
                                if (
                                    first_first_post_payment_failure
                                        .empty()
                                ) {
                                    first_first_post_payment_failure =
                                        "{q=" + std::to_string(q)
                                        + " half_power="
                                        + std::to_string(power)
                                        + " target="
                                        + std::to_string(target)
                                        + " gap="
                                        + std::to_string(gap)
                                        + " negative_mass="
                                        + negative_mass.str()
                                        + " first_reserve="
                                        + first_post_reserve.str()
                                        + " terms_needed="
                                        + std::to_string(
                                            post_terms_needed
                                        )
                                        + "}";
                                }
                            }
                            Integer previous_loss = 0;
                            for (int left = negative_inner;
                                 left <= negative_outer;
                                 ++left) {
                                if (
                                    products[
                                        static_cast<std::size_t>(left)
                                    ] >= 0
                                ) {
                                    continue;
                                }
                                const Integer loss =
                                    losing_row == 1
                                    ? absolute_value(
                                        distinguished_currents[
                                            static_cast<std::size_t>(
                                                left
                                            )
                                        ]
                                    )
                                    : absolute_value(
                                        target_currents[
                                            static_cast<std::size_t>(
                                                left
                                            )
                                        ]
                                    );
                                if (loss < previous_loss) {
                                    ++loss_magnitude_decreases;
                                }
                                previous_loss = loss;
                            }
                            for (int cutoff = negative_inner;
                                 cutoff <= negative_outer;
                                 ++cutoff) {
                                for (
                                    int threshold_index = cutoff;
                                    threshold_index <= negative_outer;
                                    ++threshold_index
                                ) {
                                    if (
                                        products[
                                            static_cast<std::size_t>(
                                                threshold_index
                                            )
                                        ] >= 0
                                    ) {
                                        continue;
                                    }
                                    const Integer threshold =
                                        losing_row == 1
                                        ? absolute_value(
                                            distinguished_currents[
                                                static_cast<std::size_t>(
                                                    threshold_index
                                                )
                                            ]
                                        )
                                        : absolute_value(
                                            target_currents[
                                                static_cast<std::size_t>(
                                                    threshold_index
                                                )
                                            ]
                                        );
                                    Integer demand = 0;
                                    for (int left = cutoff;
                                         left <= negative_outer;
                                         ++left) {
                                        if (
                                            products[
                                                static_cast<std::size_t>(
                                                    left
                                                )
                                            ] >= 0
                                        ) {
                                            continue;
                                        }
                                        const Integer loss =
                                            losing_row == 1
                                            ? absolute_value(
                                                distinguished_currents[
                                                    static_cast<
                                                        std::size_t
                                                    >(left)
                                                ]
                                            )
                                            : absolute_value(
                                                target_currents[
                                                    static_cast<
                                                        std::size_t
                                                    >(left)
                                                ]
                                            );
                                        if (loss >= threshold) {
                                            demand +=
                                                losing_row == 1
                                                ? absolute_value(
                                                    target_currents[
                                                        static_cast<
                                                            std::size_t
                                                        >(left)
                                                    ]
                                                )
                                                : absolute_value(
                                                    distinguished_currents[
                                                        static_cast<
                                                            std::size_t
                                                        >(left)
                                                    ]
                                                );
                                        }
                                    }
                                    Integer capacity = 0;
                                    for (
                                        int left = negative_outer + 1;
                                        left + gap <= maximum_index;
                                        ++left
                                    ) {
                                        if (
                                            products[
                                                static_cast<std::size_t>(
                                                    left
                                                )
                                            ] <= 0
                                        ) {
                                            continue;
                                        }
                                        const Integer resource_loss =
                                            losing_row == 1
                                            ? absolute_value(
                                                distinguished_currents[
                                                    static_cast<
                                                        std::size_t
                                                    >(left)
                                                ]
                                            )
                                            : absolute_value(
                                                target_currents[
                                                    static_cast<
                                                        std::size_t
                                                    >(left)
                                                ]
                                            );
                                        if (resource_loss >= threshold) {
                                            capacity +=
                                                losing_row == 1
                                                ? absolute_value(
                                                    target_currents[
                                                        static_cast<
                                                            std::size_t
                                                        >(left)
                                                    ]
                                                )
                                                : absolute_value(
                                                    distinguished_currents[
                                                        static_cast<
                                                            std::size_t
                                                        >(left)
                                                    ]
                                                );
                                        }
                                    }
                                    if (capacity < demand) {
                                        ++layer_payment_failures;
                                        Integer fixed_gap_suffix = 0;
                                        Integer post_block_positive = 0;
                                        for (
                                            int left = cutoff;
                                            left + gap <= maximum_index;
                                            ++left
                                        ) {
                                            fixed_gap_suffix +=
                                                products[
                                                    static_cast<std::size_t>(
                                                        left
                                                    )
                                                ];
                                            if (
                                                left > negative_outer
                                                && products[
                                                    static_cast<std::size_t>(
                                                        left
                                                    )
                                                ] > 0
                                            ) {
                                                post_block_positive +=
                                                    products[
                                                        static_cast<
                                                            std::size_t
                                                        >(left)
                                                    ];
                                            }
                                        }
                                        std::cout
                                            << "SU2_LAYER_PAYMENT_FAILURE"
                                            << " q=" << q
                                            << " half_power=" << power
                                            << " target=" << target
                                            << " gap=" << gap
                                            << " cutoff=" << cutoff
                                            << " negative_inner="
                                            << negative_inner
                                            << " negative_outer="
                                            << negative_outer
                                            << " losing_row=" << losing_row
                                            << " threshold=" << threshold
                                            << " demand=" << demand
                                            << " capacity=" << capacity
                                            << " negative_mass="
                                            << negative_mass
                                            << " post_block_positive="
                                            << post_block_positive
                                            << " boundary_reserve="
                                            << boundary_reserve
                                            << " coupled_suffix="
                                            << fixed_gap_suffix
                                            << '\n';
                                        if (
                                            first_layer_payment_failure
                                                .empty()
                                        ) {
                                            for (
                                                int left = 0;
                                                left + gap <= maximum_index;
                                                ++left
                                            ) {
                                                if (
                                                    distinguished_currents[
                                                        static_cast<
                                                            std::size_t
                                                        >(left)
                                                    ] == 0
                                                    && target_currents[
                                                        static_cast<
                                                            std::size_t
                                                        >(left)
                                                    ] == 0
                                                ) {
                                                    continue;
                                                }
                                                std::cout
                                                    << "SU2_LAYER_PROFILE"
                                                    << " left=" << left
                                                    << " distinguished="
                                                    << distinguished_currents[
                                                        static_cast<
                                                            std::size_t
                                                        >(left)
                                                    ]
                                                    << " target="
                                                    << target_currents[
                                                        static_cast<
                                                            std::size_t
                                                        >(left)
                                                    ]
                                                    << " product="
                                                    << products[
                                                        static_cast<
                                                            std::size_t
                                                        >(left)
                                                    ]
                                                    << '\n';
                                            }
                                            first_layer_payment_failure =
                                                "{q=" + std::to_string(q)
                                                + " half_power="
                                                + std::to_string(power)
                                                + " target="
                                                + std::to_string(target)
                                                + " gap="
                                                + std::to_string(gap)
                                                + " cutoff="
                                                + std::to_string(cutoff)
                                                + " losing_row="
                                                + std::to_string(
                                                    losing_row
                                                )
                                                + " threshold="
                                                + threshold.str()
                                                + " demand="
                                                + demand.str()
                                                + " capacity="
                                                + capacity.str()
                                                + " coupled_suffix="
                                                + fixed_gap_suffix.str()
                                                + "}";
                                        }
                                    }
                                }
                            }
                        }
                        if (boundary_reserve < negative_mass) {
                            ++boundary_reserve_failures;
                            if (
                                !have_first_boundary_reserve_failure
                            ) {
                                have_first_boundary_reserve_failure =
                                    true;
                                first_boundary_reserve_failure =
                                    "{q=" + std::to_string(q)
                                    + " half_power="
                                    + std::to_string(power)
                                    + " target="
                                    + std::to_string(target)
                                    + " gap="
                                    + std::to_string(gap)
                                    + " negative_mass="
                                    + negative_mass.str()
                                    + " boundary_reserve="
                                    + boundary_reserve.str()
                                    + "}";
                            }
                        }
                    }
                    for (int left = maximum_index;
                         left >= 0;
                         --left) {
                        for (int gap = maximum_index - left;
                             gap >= 1;
                             --gap) {
                            const int right = left + gap;
                            const Integer contribution =
                                current(
                                    base,
                                    distinguished,
                                    left,
                                    right
                                )
                                * current(
                                    base,
                                    target_row,
                                    left,
                                    right
                                );
                            suffix[static_cast<std::size_t>(left)]
                                  [static_cast<std::size_t>(gap)] =
                                contribution
                                + suffix[
                                    static_cast<std::size_t>(left + 1)
                                  ][static_cast<std::size_t>(gap)]
                                + suffix[
                                    static_cast<std::size_t>(left)
                                  ][static_cast<std::size_t>(gap + 1)]
                                - suffix[
                                    static_cast<std::size_t>(left + 1)
                                  ][static_cast<std::size_t>(gap + 1)];
                            const Integer fixed_gap_suffix =
                                suffix[static_cast<std::size_t>(left)]
                                      [static_cast<std::size_t>(gap)]
                                - suffix[
                                    static_cast<std::size_t>(left)
                                  ][static_cast<std::size_t>(gap + 1)];
                            if (fixed_gap_suffix < 0) {
                                ++negative_fixed_gap_suffixes;
                                if (!have_first_fixed_gap_negative) {
                                    have_first_fixed_gap_negative = true;
                                    first_fixed_gap_negative =
                                        "{q=" + std::to_string(q)
                                        + " half_power="
                                        + std::to_string(power)
                                        + " target="
                                        + std::to_string(target)
                                        + " cutoff="
                                        + std::to_string(left)
                                        + " gap="
                                        + std::to_string(gap)
                                        + " value="
                                        + fixed_gap_suffix.str()
                                        + "}";
                                }
                            }
                            const Integer& previous_target =
                                previous_fixed_gap_suffix[
                                    static_cast<std::size_t>(left)
                                ][static_cast<std::size_t>(gap)];
                            if (fixed_gap_suffix < previous_target) {
                                ++target_label_decreases;
                                if (
                                    !have_first_target_label_decrease
                                ) {
                                    have_first_target_label_decrease =
                                        true;
                                    first_target_label_decrease =
                                        "{q=" + std::to_string(q)
                                        + " half_power="
                                        + std::to_string(power)
                                        + " target="
                                        + std::to_string(target)
                                        + " cutoff="
                                        + std::to_string(left)
                                        + " gap="
                                        + std::to_string(gap)
                                        + " previous="
                                        + previous_target.str()
                                        + " value="
                                        + fixed_gap_suffix.str()
                                        + "}";
                                }
                            }
                            previous_fixed_gap_suffix[
                                static_cast<std::size_t>(left)
                            ][static_cast<std::size_t>(gap)] =
                                fixed_gap_suffix;
                            const Integer fixed_base_suffix =
                                suffix[static_cast<std::size_t>(left)]
                                      [static_cast<std::size_t>(gap)]
                                - suffix[
                                    static_cast<std::size_t>(left + 1)
                                  ][static_cast<std::size_t>(gap)];
                            if (fixed_base_suffix < 0) {
                                ++negative_fixed_base_suffixes;
                                if (!have_first_fixed_base_negative) {
                                    have_first_fixed_base_negative = true;
                                    first_fixed_base_negative =
                                        "{q=" + std::to_string(q)
                                        + " half_power="
                                        + std::to_string(power)
                                        + " target="
                                        + std::to_string(target)
                                        + " base="
                                        + std::to_string(left)
                                        + " minimum_gap="
                                        + std::to_string(gap)
                                        + " value="
                                        + fixed_base_suffix.str()
                                        + "}";
                                }
                            }
                            ++orthants;
                            if (
                                suffix[static_cast<std::size_t>(left)]
                                      [static_cast<std::size_t>(gap)] < 0
                            ) {
                                std::cout
                                    << "SU2_STRUCTURED_WEDGE_ORTHANT"
                                    << " negative"
                                    << " q=" << q
                                    << " half_power=" << power
                                    << " target=" << target
                                    << " cutoff=" << left
                                    << " minimum_gap=" << gap
                                    << " value="
                                    << suffix[
                                        static_cast<std::size_t>(left)
                                      ][static_cast<std::size_t>(gap)]
                                    << '\n';
                                return EXIT_SUCCESS;
                            }
                        }
                    }
                    ++cases;
                }
            }
        }

        std::cout
            << "SU2_STRUCTURED_WEDGE_ORTHANT"
            << " maximum_half_label=" << maximum_q
            << " maximum_half_power=" << maximum_power
            << " minimum_half_power=" << minimum_power
            << " single_case=" << (single_case ? 1 : 0)
            << " tp2_only=" << (tp2_only ? 1 : 0)
            << " kernel_profiles=" << kernel_profiles
            << " cases=" << cases
            << " orthants=" << orthants
            << " negative_fixed_gap_suffixes="
            << negative_fixed_gap_suffixes
            << " negative_fixed_base_suffixes="
            << negative_fixed_base_suffixes
            << " negative_direct_suffixes="
            << negative_direct_suffixes
            << " positive_cross_suffixes="
            << positive_cross_suffixes
            << " target_label_decreases="
            << target_label_decreases
            << " maximum_product_negative_blocks="
            << maximum_product_negative_blocks
            << " maximum_current_negative_blocks="
            << maximum_current_negative_blocks
            << " negative_endpoint_mismatches="
            << negative_endpoint_mismatches
            << " nonnested_negative_intervals="
            << nonnested_negative_intervals
            << " negative_orientation_failures="
            << negative_orientation_failures
            << " loss_magnitude_decreases="
            << loss_magnitude_decreases
            << " layer_payment_failures="
            << layer_payment_failures
            << " post_block_payment_failures="
            << post_block_payment_failures
            << " first_post_payment_failures="
            << first_post_payment_failures
            << " negative_full_adjacent_minors="
            << negative_full_adjacent_minors
            << " minimum_active_adjacent_minor="
            << (
                have_minimum_active_adjacent_minor
                    ? minimum_active_adjacent_minor.str()
                    : "{}"
            )
            << " minimum_active_adjacent_minor_case="
            << (
                minimum_active_adjacent_minor_case.empty()
                    ? "{}"
                    : minimum_active_adjacent_minor_case
            )
            << " maximum_post_terms_needed="
            << maximum_post_terms_needed
            << " maximum_post_terms_case="
            << (
                maximum_post_terms_case.empty()
                    ? "{}"
                    : maximum_post_terms_case
            )
            << " boundary_reserve_failures="
            << boundary_reserve_failures
            << " first_fixed_gap_negative="
            << (
                have_first_fixed_gap_negative
                    ? first_fixed_gap_negative
                    : "{}"
            )
            << " first_fixed_base_negative="
            << (
                have_first_fixed_base_negative
                    ? first_fixed_base_negative
                    : "{}"
            )
            << " first_direct_negative="
            << (
                have_first_direct_negative
                    ? first_direct_negative
                    : "{}"
            )
            << " first_cross_positive="
            << (
                have_first_cross_positive
                    ? first_cross_positive
                    : "{}"
            )
            << " first_target_label_decrease="
            << (
                have_first_target_label_decrease
                    ? first_target_label_decrease
                    : "{}"
            )
            << " first_boundary_reserve_failure="
            << (
                have_first_boundary_reserve_failure
                    ? first_boundary_reserve_failure
                    : "{}"
            )
            << " first_layer_payment_failure="
            << (
                first_layer_payment_failure.empty()
                    ? "{}"
                    : first_layer_payment_failure
            )
            << " first_post_block_payment_failure="
            << (
                first_post_block_payment_failure.empty()
                    ? "{}"
                    : first_post_block_payment_failure
            )
            << " first_first_post_payment_failure="
            << (
                first_first_post_payment_failure.empty()
                    ? "{}"
                    : first_first_post_payment_failure
            )
            << " first_negative_full_adjacent_minor="
            << (
                first_negative_full_adjacent_minor.empty()
                    ? "{}"
                    : first_negative_full_adjacent_minor
            )
            << " result=NO_NEGATIVE_ORTHANT"
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
