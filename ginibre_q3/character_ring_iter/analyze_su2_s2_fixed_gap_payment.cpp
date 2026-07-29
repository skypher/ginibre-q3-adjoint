#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

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

long long fusion_interval_length(int q, int row, int index) {
    const int upper = std::min(2 * q, index + row);
    const int lower = std::abs(index - row);
    return upper >= lower
        ? static_cast<long long>(upper - lower + 1)
        : 0LL;
}

long long current(int q, int row, int gap, int left) {
    const int right = left + gap;
    const long long profile_left =
        left <= 2 * q ? 1LL : 0LL;
    const long long profile_right =
        right <= 2 * q ? 1LL : 0LL;
    return
        profile_left * fusion_interval_length(q, row, right)
        - profile_right * fusion_interval_length(q, row, left);
}

int sign(long long value) {
    return value > 0 ? 1 : (value < 0 ? -1 : 0);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_s2_fixed_gap_payment maximum_q"
            );
        }
        const int maximum_q = parse_positive(argv[1], "maximum_q");

        std::size_t cases = 0;
        std::size_t suffixes = 0;
        std::size_t negative_entries = 0;
        std::size_t reflection_failures = 0;
        std::size_t boundary_reserve_failures = 0;
        std::size_t pointwise_boundary_failures = 0;
        std::size_t ferrers_layer_failures = 0;
        int maximum_negative_blocks = 0;
        int maximum_sign_changes = 0;
        std::string first_negative_entry;
        std::string first_reflection_failure;
        std::string first_boundary_reserve_failure;
        std::string first_pointwise_boundary_failure;
        std::string first_ferrers_layer_failure;
        long long minimum_positive_case_margin =
            std::numeric_limits<long long>::max();
        std::string minimum_margin_case;
        std::string first_reverse_orientation;

        for (int q = 1; q <= maximum_q; ++q) {
            for (int target = 0; target <= 5 * q; ++target) {
                for (int gap = 1; gap <= 3 * q; ++gap) {
                    std::vector<long long> products(
                        static_cast<std::size_t>(2 * q + 1)
                    );
                    int negative_blocks = 0;
                    int sign_changes = 0;
                    int previous_sign = 0;
                    bool inside_negative = false;
                    for (int left = 0; left <= 2 * q; ++left) {
                        products[static_cast<std::size_t>(left)] =
                            current(q, q, gap, left)
                            * current(q, target, gap, left);
                    }
                    for (int left = 0; left <= 2 * q; ++left) {
                        const long long product =
                            products[static_cast<std::size_t>(left)];
                        const int product_sign = sign(product);
                        if (product_sign < 0) {
                            ++negative_entries;
                            if (first_negative_entry.empty()) {
                                first_negative_entry =
                                    "{q=" + std::to_string(q)
                                    + " target="
                                    + std::to_string(target)
                                    + " gap="
                                    + std::to_string(gap)
                                    + " left="
                                    + std::to_string(left)
                                    + " value="
                                    + std::to_string(product)
                                    + "}";
                            }
                            if (!inside_negative) {
                                ++negative_blocks;
                                inside_negative = true;
                            }
                            const int boundary_depth =
                                2 * q - gap - left + 1;
                            const long long distinguished_current =
                                current(q, q, gap, left);
                            const long long target_current =
                                current(q, target, gap, left);
                            const long long distinguished_boundary =
                                fusion_interval_length(
                                    q,
                                    q,
                                    2 * q + boundary_depth
                                );
                            const long long target_boundary =
                                fusion_interval_length(
                                    q,
                                    target,
                                    2 * q + boundary_depth
                                );
                            if (
                                boundary_depth < 1
                                || boundary_depth > gap
                                || std::abs(distinguished_current)
                                    > distinguished_boundary
                                || std::abs(target_current)
                                    > target_boundary
                            ) {
                                ++pointwise_boundary_failures;
                                if (
                                    first_pointwise_boundary_failure
                                        .empty()
                                ) {
                                    first_pointwise_boundary_failure =
                                        "{q=" + std::to_string(q)
                                        + " target="
                                        + std::to_string(target)
                                        + " gap="
                                        + std::to_string(gap)
                                        + " left="
                                        + std::to_string(left)
                                        + " depth="
                                        + std::to_string(
                                            boundary_depth
                                        )
                                        + " currents=("
                                        + std::to_string(
                                            distinguished_current
                                        )
                                        + ","
                                        + std::to_string(target_current)
                                        + ") boundaries=("
                                        + std::to_string(
                                            distinguished_boundary
                                        )
                                        + ","
                                        + std::to_string(
                                            target_boundary
                                        )
                                        + ")}";
                                }
                            }
                            const int reflected = 3 * q - left;
                            if (
                                reflected < 0
                                || reflected > 2 * q
                                || product
                                    + products[
                                        static_cast<std::size_t>(
                                            reflected
                                        )
                                      ] < 0
                            ) {
                                ++reflection_failures;
                                if (first_reflection_failure.empty()) {
                                    first_reflection_failure =
                                        "{q=" + std::to_string(q)
                                        + " target="
                                        + std::to_string(target)
                                        + " gap="
                                        + std::to_string(gap)
                                        + " left="
                                        + std::to_string(left)
                                        + " reflected="
                                        + std::to_string(reflected)
                                        + " value="
                                        + std::to_string(product)
                                        + "}";
                                }
                            }
                        } else if (product_sign > 0) {
                            inside_negative = false;
                        }
                        if (
                            product_sign != 0
                            && previous_sign != 0
                            && product_sign != previous_sign
                        ) {
                            ++sign_changes;
                        }
                        if (product_sign != 0) {
                            previous_sign = product_sign;
                        }
                    }
                    maximum_negative_blocks = std::max(
                        maximum_negative_blocks,
                        negative_blocks
                    );
                    maximum_sign_changes = std::max(
                        maximum_sign_changes,
                        sign_changes
                    );

                    long long negative_mass = 0;
                    long long boundary_reserve = 0;
                    std::vector<long long> distinguished_loss_layers(
                        static_cast<std::size_t>(gap + 1)
                    );
                    std::vector<long long>
                        maximum_target_gain_by_layer(
                            static_cast<std::size_t>(gap + 1)
                        );
                    std::vector<long long> target_loss_layers(
                        static_cast<std::size_t>(gap + 1)
                    );
                    std::vector<long long>
                        maximum_distinguished_gain_by_layer(
                        static_cast<std::size_t>(gap + 1)
                    );
                    for (int left = 0; left <= 2 * q; ++left) {
                        const long long product =
                            products[static_cast<std::size_t>(left)];
                        if (product < 0) {
                            negative_mass -= product;
                            const long long distinguished_loss =
                                -current(q, q, gap, left);
                            const long long target_gain =
                                current(q, target, gap, left);
                            if (
                                distinguished_loss > 0
                                && target_gain > 0
                            ) {
                                for (long long layer = 1;
                                     layer <= distinguished_loss;
                                     ++layer) {
                                    ++distinguished_loss_layers[
                                        static_cast<std::size_t>(layer)
                                    ];
                                    maximum_target_gain_by_layer[
                                        static_cast<std::size_t>(layer)
                                    ] = std::max(
                                        maximum_target_gain_by_layer[
                                            static_cast<std::size_t>(
                                                layer
                                            )
                                        ],
                                        target_gain
                                    );
                                }
                            } else if (
                                distinguished_loss < 0
                                && target_gain < 0
                            ) {
                                if (first_reverse_orientation.empty()) {
                                    first_reverse_orientation =
                                        "{q=" + std::to_string(q)
                                        + " target="
                                        + std::to_string(target)
                                        + " gap="
                                        + std::to_string(gap)
                                        + " left="
                                        + std::to_string(left)
                                        + " distinguished="
                                        + std::to_string(
                                            -distinguished_loss
                                        )
                                        + " target="
                                        + std::to_string(target_gain)
                                        + "}";
                                }
                                const long long distinguished_gain =
                                    -distinguished_loss;
                                const long long target_loss =
                                    -target_gain;
                                for (long long layer = 1;
                                     layer <= target_loss;
                                     ++layer) {
                                    ++target_loss_layers[
                                        static_cast<std::size_t>(layer)
                                    ];
                                    maximum_distinguished_gain_by_layer[
                                        static_cast<std::size_t>(layer)
                                    ] = std::max(
                                        maximum_distinguished_gain_by_layer[
                                            static_cast<std::size_t>(
                                                layer
                                            )
                                        ],
                                        distinguished_gain
                                    );
                                }
                            } else {
                                throw std::runtime_error(
                                    "invalid negative-product orientation"
                                );
                            }
                        }
                        if (left + gap > 2 * q) {
                            boundary_reserve += product;
                        }
                    }
                    for (int layer = 1; layer <= gap; ++layer) {
                        const long long distinguished_capacity =
                            fusion_interval_length(
                                q,
                                q,
                                2 * q + layer
                            );
                        const long long target_capacity =
                            fusion_interval_length(
                                q,
                                target,
                                2 * q + layer
                            );
                        if (
                            distinguished_loss_layers[
                                static_cast<std::size_t>(layer)
                            ] > distinguished_capacity
                            || maximum_target_gain_by_layer[
                                static_cast<std::size_t>(layer)
                            ] > target_capacity
                            || target_loss_layers[
                                static_cast<std::size_t>(layer)
                            ] > target_capacity
                            || maximum_distinguished_gain_by_layer[
                                static_cast<std::size_t>(layer)
                            ] > distinguished_capacity
                        ) {
                            ++ferrers_layer_failures;
                            if (first_ferrers_layer_failure.empty()) {
                                first_ferrers_layer_failure =
                                    "{q=" + std::to_string(q)
                                    + " target="
                                    + std::to_string(target)
                                    + " gap="
                                    + std::to_string(gap)
                                    + " layer="
                                    + std::to_string(layer)
                                    + " count="
                                    + std::to_string(
                                        distinguished_loss_layers[
                                            static_cast<std::size_t>(
                                                layer
                                            )
                                        ]
                                    )
                                    + " max_gain="
                                    + std::to_string(
                                        maximum_target_gain_by_layer[
                                            static_cast<std::size_t>(
                                                layer
                                            )
                                        ]
                                    )
                                    + " reverse_count="
                                    + std::to_string(
                                        target_loss_layers[
                                            static_cast<std::size_t>(
                                                layer
                                            )
                                        ]
                                    )
                                    + " reverse_max_gain="
                                    + std::to_string(
                                        maximum_distinguished_gain_by_layer[
                                            static_cast<std::size_t>(
                                                layer
                                            )
                                        ]
                                    )
                                    + " capacities=("
                                    + std::to_string(
                                        distinguished_capacity
                                    )
                                    + ","
                                    + std::to_string(target_capacity)
                                    + ")}";
                            }
                        }
                    }
                    if (boundary_reserve < negative_mass) {
                        ++boundary_reserve_failures;
                        if (first_boundary_reserve_failure.empty()) {
                            first_boundary_reserve_failure =
                                "{q=" + std::to_string(q)
                                + " target="
                                + std::to_string(target)
                                + " gap="
                                + std::to_string(gap)
                                + " negative_mass="
                                + std::to_string(negative_mass)
                                + " boundary_reserve="
                                + std::to_string(boundary_reserve)
                                + "}";
                        }
                    }
                    if (
                        negative_mass > 0
                        && boundary_reserve - negative_mass
                            < minimum_positive_case_margin
                    ) {
                        minimum_positive_case_margin =
                            boundary_reserve - negative_mass;
                        minimum_margin_case =
                            "{q=" + std::to_string(q)
                            + " target="
                            + std::to_string(target)
                            + " gap="
                            + std::to_string(gap)
                            + " negative_mass="
                            + std::to_string(negative_mass)
                            + " boundary_reserve="
                            + std::to_string(boundary_reserve)
                            + " margin="
                            + std::to_string(
                                boundary_reserve - negative_mass
                            )
                            + "}";
                    }

                    long long suffix = 0;
                    for (int left = 2 * q; left >= 0; --left) {
                        suffix += products[
                            static_cast<std::size_t>(left)
                        ];
                        ++suffixes;
                        if (suffix < 0) {
                            std::cout
                                << "SU2_S2_FIXED_GAP_PAYMENT"
                                << " negative_suffix"
                                << " q=" << q
                                << " target=" << target
                                << " gap=" << gap
                                << " cutoff=" << left
                                << " value=" << suffix
                                << '\n';
                            return EXIT_SUCCESS;
                        }
                    }
                    ++cases;
                }
            }
        }

        std::cout
            << "SU2_S2_FIXED_GAP_PAYMENT"
            << " maximum_q=" << maximum_q
            << " cases=" << cases
            << " suffixes=" << suffixes
            << " negative_entries=" << negative_entries
            << " maximum_negative_blocks=" << maximum_negative_blocks
            << " maximum_sign_changes=" << maximum_sign_changes
            << " reflection_failures=" << reflection_failures
            << " boundary_reserve_failures="
            << boundary_reserve_failures
            << " pointwise_boundary_failures="
            << pointwise_boundary_failures
            << " ferrers_layer_failures="
            << ferrers_layer_failures
            << " first_negative_entry="
            << (
                first_negative_entry.empty()
                    ? "{}"
                    : first_negative_entry
            )
            << " first_reflection_failure="
            << (
                first_reflection_failure.empty()
                    ? "{}"
                    : first_reflection_failure
            )
            << " first_boundary_reserve_failure="
            << (
                first_boundary_reserve_failure.empty()
                    ? "{}"
                    : first_boundary_reserve_failure
            )
            << " first_pointwise_boundary_failure="
            << (
                first_pointwise_boundary_failure.empty()
                    ? "{}"
                    : first_pointwise_boundary_failure
            )
            << " first_ferrers_layer_failure="
            << (
                first_ferrers_layer_failure.empty()
                    ? "{}"
                    : first_ferrers_layer_failure
            )
            << " minimum_margin_case="
            << (
                minimum_margin_case.empty()
                    ? "{}"
                    : minimum_margin_case
            )
            << " first_reverse_orientation="
            << (
                first_reverse_orientation.empty()
                    ? "{}"
                    : first_reverse_orientation
            )
            << " result=NO_NEGATIVE_SUFFIX"
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
