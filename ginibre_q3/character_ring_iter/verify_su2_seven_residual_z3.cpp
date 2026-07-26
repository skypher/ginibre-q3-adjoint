#include <map>

#define main verify_su2_four_minus_shallow_main
#include "verify_su2_four_minus_shallow_z3.cpp"
#undef main

namespace {

struct ResidualOrbit {
    int minus_count;
    int odd_count;
    int odd_minus_count;
};

constexpr std::array<ResidualOrbit, 8> residual_orbits{{
    {4, 0, 0},
    {4, 2, 0},
    {4, 2, 1},
    {4, 4, 1},
    {4, 4, 2},
    {4, 6, 3},
    {6, 0, 0},
    {6, 2, 1},
}};

int residual_bit_count(unsigned int value) {
    int result = 0;
    while (value != 0U) {
        result += static_cast<int>(value & 1U);
        value >>= 1U;
    }
    return result;
}

int orbit_label_parity(
    const ResidualOrbit& orbit,
    int index
) {
    if (index < orbit.minus_count) {
        return index < orbit.odd_minus_count ? 1 : 0;
    }
    const int plus_index = index - orbit.minus_count;
    const int odd_plus =
        orbit.odd_count - orbit.odd_minus_count;
    return plus_index < odd_plus ? 1 : 0;
}

int orbit_label_class(
    const ResidualOrbit& orbit,
    int index
) {
    if (index < orbit.minus_count) {
        return orbit_label_parity(orbit, index) == 1 ? 0 : 1;
    }
    return orbit_label_parity(orbit, index) == 1 ? 2 : 3;
}

std::array<int, 3> mask_triple(unsigned int mask) {
    std::array<int, 3> result{};
    int next = 0;
    for (int index = 0; index < 7; ++index) {
        if (((mask >> index) & 1U) != 0U) {
            result[static_cast<std::size_t>(next)] = index;
            ++next;
        }
    }
    return result;
}

std::array<int, 4> mask_complement(unsigned int mask) {
    std::array<int, 4> result{};
    int next = 0;
    for (int index = 0; index < 7; ++index) {
        if (((mask >> index) & 1U) == 0U) {
            result[static_cast<std::size_t>(next)] = index;
            ++next;
        }
    }
    return result;
}

std::vector<unsigned int> signed_cut_masks(
    const ResidualOrbit& orbit,
    bool negative
) {
    std::vector<unsigned int> result;
    const unsigned int minus_mask =
        (1U << orbit.minus_count) - 1U;
    for (unsigned int mask = 0U; mask < (1U << 7U); ++mask) {
        if (residual_bit_count(mask) != 3) {
            continue;
        }
        const bool is_negative =
            (residual_bit_count(mask & minus_mask) & 1)
            != 0;
        if (is_negative == negative) {
            result.push_back(mask);
        }
    }
    return result;
}

std::vector<unsigned int> selected_orbit_masks(
    const ResidualOrbit& orbit
) {
    std::map<std::array<int, 4>, unsigned int> representatives;
    for (const unsigned int mask :
         signed_cut_masks(orbit, true)) {
        std::array<int, 4> counts{};
        for (int index = 0; index < 7; ++index) {
            if (((mask >> index) & 1U) != 0U) {
                ++counts[static_cast<std::size_t>(
                    orbit_label_class(orbit, index)
                )];
            }
        }
        representatives.emplace(counts, mask);
    }
    std::vector<unsigned int> result;
    result.reserve(representatives.size());
    for (const auto& [counts, mask] : representatives) {
        static_cast<void>(counts);
        result.push_back(mask);
    }
    return result;
}

z3::expr residual_rank(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 7>& labels,
    unsigned int mask
) {
    const std::array<int, 3> cut = mask_triple(mask);
    const std::array<int, 4> rest = mask_complement(mask);
    const z3::expr active = fusion(
        k,
        labels[static_cast<std::size_t>(cut[0])],
        labels[static_cast<std::size_t>(cut[1])],
        labels[static_cast<std::size_t>(cut[2])]
    );
    return z3::ite(
        active,
        interval_rank(
            ctx,
            k,
            labels[static_cast<std::size_t>(rest[0])],
            labels[static_cast<std::size_t>(rest[1])],
            labels[static_cast<std::size_t>(rest[2])],
            labels[static_cast<std::size_t>(rest[3])]
        ),
        ctx.int_val(0)
    );
}

z3::expr residual_selected_rank_chamber(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 7>& labels,
    const ResidualOrbit& orbit,
    unsigned int mask,
    int selected_wall_mask,
    int selected_interval_mask
) {
    const std::array<int, 3> cut = mask_triple(mask);
    const std::array<int, 4> rest = mask_complement(mask);
    const int first_parity =
        orbit_label_parity(orbit, rest[0])
        ^ orbit_label_parity(orbit, rest[1]);
    const int second_parity =
        orbit_label_parity(orbit, rest[2])
        ^ orbit_label_parity(orbit, rest[3]);
    if (first_parity != second_parity) {
        return ctx.int_val(0);
    }

    const z3::expr& first =
        labels[static_cast<std::size_t>(rest[0])];
    const z3::expr& second =
        labels[static_cast<std::size_t>(rest[1])];
    const z3::expr& third =
        labels[static_cast<std::size_t>(rest[2])];
    const z3::expr& fourth =
        labels[static_cast<std::size_t>(rest[3])];
    const z3::expr first_lower =
        ((selected_interval_mask >> 0) & 1) != 0
            ? first - second
            : second - first;
    const z3::expr second_lower =
        ((selected_interval_mask >> 1) & 1) != 0
            ? third - fourth
            : fourth - third;
    const z3::expr first_upper =
        ((selected_wall_mask >> 1) & 1) != 0
            ? first + second
            : 2 * k - first - second;
    const z3::expr second_upper =
        ((selected_wall_mask >> 2) & 1) != 0
            ? third + fourth
            : 2 * k - third - fourth;
    const z3::expr lower =
        ((selected_interval_mask >> 2) & 1) != 0
            ? first_lower
            : second_lower;
    const z3::expr upper =
        ((selected_interval_mask >> 3) & 1) != 0
            ? first_upper
            : second_upper;
    const z3::expr active = fusion(
        k,
        labels[static_cast<std::size_t>(cut[0])],
        labels[static_cast<std::size_t>(cut[1])],
        labels[static_cast<std::size_t>(cut[2])]
    );
    return z3::ite(
        active && lower <= upper,
        (upper - lower) / 2 + 1,
        ctx.int_val(0)
    );
}

z3::expr residual_triple_output(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 7>& labels,
    const std::array<int, 3>& cut,
    int output,
    bool reflect
) {
    const z3::expr first = reflect
        ? k - labels[static_cast<std::size_t>(cut[0])]
        : labels[static_cast<std::size_t>(cut[0])];
    return interval_rank(
        ctx,
        k,
        first,
        labels[static_cast<std::size_t>(cut[1])],
        labels[static_cast<std::size_t>(cut[2])],
        ctx.int_val(output)
    );
}

z3::expr residual_fourfold_output(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 7>& labels,
    const std::array<int, 4>& rest,
    int output,
    bool reflect
) {
    const z3::expr first = reflect
        ? k - labels[static_cast<std::size_t>(rest[0])]
        : labels[static_cast<std::size_t>(rest[0])];
    return four_tensor_output(
        ctx,
        k,
        first,
        labels[static_cast<std::size_t>(rest[1])],
        labels[static_cast<std::size_t>(rest[2])],
        labels[static_cast<std::size_t>(rest[3])],
        output
    );
}

z3::expr residual_selected_local(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 7>& labels,
    const ResidualOrbit& orbit,
    unsigned int mask,
    int level_parity
) {
    const std::array<int, 3> cut = mask_triple(mask);
    const std::array<int, 4> rest = mask_complement(mask);
    int common_parity = 0;
    for (const int index : cut) {
        common_parity ^= orbit_label_parity(orbit, index);
    }
    const int reflected_parity =
        (level_parity + common_parity) & 1;
    const int low_maximum = common_parity + 8;
    z3::expr result = ctx.int_val(0);
    for (int j = 0; j <= 4; ++j) {
        const int low = common_parity + 2 * j;
        result = result + bounded_product(
            ctx,
            residual_triple_output(
                ctx, k, labels, cut, low, false
            ),
            residual_fourfold_output(
                ctx, k, labels, rest, low, false
            ),
            low + 1
        );

        const int reflected_low = reflected_parity + 2 * j;
        const z3::expr high = k - reflected_low;
        const z3::expr reflected_term = bounded_product(
            ctx,
            residual_triple_output(
                ctx, k, labels, cut, reflected_low, true
            ),
            residual_fourfold_output(
                ctx, k, labels, rest, reflected_low, true
            ),
            reflected_low + 1
        );
        result = result + z3::ite(
            high > low_maximum,
            reflected_term,
            ctx.int_val(0)
        );
    }
    return result;
}

z3::expr residual_pair_reservoir(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 7>& labels,
    const ResidualOrbit& orbit,
    int level_parity,
    int equality_mask = -1
) {
    z3::expr result = ctx.int_val(0);
    int equality_index = 0;
    for (int first = 0; first < 7; ++first) {
        for (int second = first + 1; second < 7; ++second) {
            const bool same_sign =
                (first < orbit.minus_count)
                == (second < orbit.minus_count);
            if (!same_sign) {
                continue;
            }
            const bool equality_active =
                equality_mask >= 0
                && ((equality_mask >> equality_index) & 1) != 0;
            std::array<int, 5> rest{};
            int next = 0;
            for (int index = 0; index < 7; ++index) {
                if (index != first && index != second) {
                    rest[static_cast<std::size_t>(next)] = index;
                    ++next;
                }
            }
            const int common_parity =
                orbit_label_parity(orbit, rest[0])
                ^ orbit_label_parity(orbit, rest[1]);
            const int reflected_parity =
                (level_parity + common_parity) & 1;
            const int low_maximum = common_parity + 8;
            z3::expr contribution = ctx.int_val(0);
            for (int j = 0; j <= 4; ++j) {
                const int low = common_parity + 2 * j;
                const z3::expr low_pair = fusion(
                    k,
                    labels[static_cast<std::size_t>(rest[0])],
                    labels[static_cast<std::size_t>(rest[1])],
                    ctx.int_val(low)
                );
                const z3::expr low_triple = interval_rank(
                    ctx,
                    k,
                    labels[static_cast<std::size_t>(rest[2])],
                    labels[static_cast<std::size_t>(rest[3])],
                    labels[static_cast<std::size_t>(rest[4])],
                    ctx.int_val(low)
                );
                contribution = contribution + z3::ite(
                    low_pair, low_triple, ctx.int_val(0)
                );

                const int reflected_low =
                    reflected_parity + 2 * j;
                const z3::expr high = k - reflected_low;
                const z3::expr high_pair = fusion(
                    k,
                    k - labels[
                        static_cast<std::size_t>(rest[0])
                    ],
                    labels[static_cast<std::size_t>(rest[1])],
                    ctx.int_val(reflected_low)
                );
                const z3::expr high_triple = interval_rank(
                    ctx,
                    k,
                    k - labels[
                        static_cast<std::size_t>(rest[2])
                    ],
                    labels[static_cast<std::size_t>(rest[3])],
                    labels[static_cast<std::size_t>(rest[4])],
                    ctx.int_val(reflected_low)
                );
                contribution = contribution + z3::ite(
                    high > low_maximum && high_pair,
                    high_triple,
                    ctx.int_val(0)
                );
            }
            if (equality_mask < 0) {
                result = result + z3::ite(
                    labels[static_cast<std::size_t>(first)]
                        == labels[static_cast<std::size_t>(second)],
                    contribution,
                    ctx.int_val(0)
                );
            } else if (equality_active) {
                result = result + contribution;
            }
            ++equality_index;
        }
    }
    return result;
}

z3::expr residual_positive_reservoir(
    z3::context& ctx,
    const z3::expr& k,
    const std::array<z3::expr, 7>& labels,
    const ResidualOrbit& orbit
) {
    z3::expr result = ctx.int_val(0);
    for (const unsigned int mask :
         signed_cut_masks(orbit, false)) {
        result = result + residual_rank(
            ctx, k, labels, mask
        );
    }
    return result;
}

QueryResult verify_residual_query(
    int orbit_index,
    int level_parity,
    int selected_orbit,
    int cap_mask = -1,
    int selected_wall_mask = -1,
    int selected_interval_mask = -1,
    int equality_mask = -1,
    int target_mode = 0,
    int selected_rank_value = 0
) {
    const ResidualOrbit& orbit =
        residual_orbits[static_cast<std::size_t>(orbit_index)];
    const std::vector<unsigned int> representatives =
        selected_orbit_masks(orbit);
    const unsigned int selected_mask =
        representatives[static_cast<std::size_t>(selected_orbit)];

    z3::context ctx;
    z3::solver solver(ctx);
    const z3::expr k =
        2 * ctx.int_const("k_half") + level_parity;
    std::array<z3::expr, 7> labels{
        ctx.int_const("l0"),
        ctx.int_const("l1"),
        ctx.int_const("l2"),
        ctx.int_const("l3"),
        ctx.int_const("l4"),
        ctx.int_const("l5"),
        ctx.int_const("l6")
    };
    for (int index = 0; index < 7; ++index) {
        const std::string name =
            "l" + std::to_string(index) + "_half";
        labels[static_cast<std::size_t>(index)] =
            2 * ctx.int_const(name.c_str())
            + orbit_label_parity(orbit, index);
    }

    solver.add(k >= 1);
    for (const z3::expr& label : labels) {
        solver.add(label >= 1 && label <= k);
        if (target_mode != 0) {
            solver.add(label >= 2 && label <= k - 2);
        }
    }
    if (cap_mask >= 0) {
        for (int index = 0; index < 7; ++index) {
            const z3::expr& label =
                labels[static_cast<std::size_t>(index)];
            solver.add(
                ((cap_mask >> index) & 1) != 0
                    ? 2 * label <= k
                    : 2 * label > k
            );
        }
    }
    if (selected_wall_mask >= 0) {
        const std::array<int, 3> cut =
            mask_triple(selected_mask);
        const std::array<int, 4> rest =
            mask_complement(selected_mask);
        const std::array<std::array<int, 2>, 3> pairs{{
            {cut[0], cut[1]},
            {rest[0], rest[1]},
            {rest[2], rest[3]},
        }};
        for (int index = 0; index < 3; ++index) {
            const auto& pair =
                pairs[static_cast<std::size_t>(index)];
            const z3::expr pair_sum =
                labels[static_cast<std::size_t>(pair[0])]
                + labels[static_cast<std::size_t>(pair[1])];
            solver.add(
                ((selected_wall_mask >> index) & 1) != 0
                    ? pair_sum <= k
                    : pair_sum > k
            );
        }
    }
    if (selected_interval_mask >= 0) {
        const std::array<int, 4> rest =
            mask_complement(selected_mask);
        const z3::expr& first =
            labels[static_cast<std::size_t>(rest[0])];
        const z3::expr& second =
            labels[static_cast<std::size_t>(rest[1])];
        const z3::expr& third =
            labels[static_cast<std::size_t>(rest[2])];
        const z3::expr& fourth =
            labels[static_cast<std::size_t>(rest[3])];
        const z3::expr first_lower = zabs(first - second);
        const z3::expr second_lower = zabs(third - fourth);
        const z3::expr first_upper =
            zmin(first + second, 2 * k - first - second);
        const z3::expr second_upper =
            zmin(third + fourth, 2 * k - third - fourth);
        solver.add(
            ((selected_interval_mask >> 0) & 1) != 0
                ? first >= second
                : first < second
        );
        solver.add(
            ((selected_interval_mask >> 1) & 1) != 0
                ? third >= fourth
                : third < fourth
        );
        solver.add(
            ((selected_interval_mask >> 2) & 1) != 0
                ? first_lower >= second_lower
                : first_lower < second_lower
        );
        solver.add(
            ((selected_interval_mask >> 3) & 1) != 0
                ? first_upper <= second_upper
                : first_upper > second_upper
        );
    }
    for (int minus = 0; minus < orbit.minus_count; ++minus) {
        for (int plus = orbit.minus_count; plus < 7; ++plus) {
            solver.add(
                labels[static_cast<std::size_t>(minus)]
                != labels[static_cast<std::size_t>(plus)]
            );
        }
    }
    if (equality_mask >= 0) {
        int equality_index = 0;
        for (int first = 0; first < 7; ++first) {
            for (int second = first + 1; second < 7; ++second) {
                const bool same_sign =
                    (first < orbit.minus_count)
                    == (second < orbit.minus_count);
                if (!same_sign) {
                    continue;
                }
                const z3::expr& first_label =
                    labels[static_cast<std::size_t>(first)];
                const z3::expr& second_label =
                    labels[static_cast<std::size_t>(second)];
                solver.add(
                    ((equality_mask >> equality_index) & 1) != 0
                        ? first_label == second_label
                        : first_label != second_label
                );
                ++equality_index;
            }
        }
    }

    for (int label_class = 0; label_class < 4; ++label_class) {
        for (int selected = 0; selected < 2; ++selected) {
            std::vector<int> group;
            for (int index = 0; index < 7; ++index) {
                const bool in_selected =
                    ((selected_mask >> index) & 1U) != 0U;
                if (orbit_label_class(orbit, index) == label_class
                    && in_selected == (selected != 0)) {
                    group.push_back(index);
                }
            }
            for (std::size_t i = 1U; i < group.size(); ++i) {
                solver.add(
                    labels[static_cast<std::size_t>(group[i - 1U])]
                    <= labels[static_cast<std::size_t>(group[i])]
                );
            }
        }
    }

    const std::vector<unsigned int> negative_masks =
        signed_cut_masks(orbit, true);
    const auto rank_for_mask = [&](unsigned int mask) {
        return mask == selected_mask
                && selected_wall_mask >= 0
                && selected_interval_mask >= 0
            ? residual_selected_rank_chamber(
                ctx,
                k,
                labels,
                orbit,
                mask,
                selected_wall_mask,
                selected_interval_mask
            )
            : residual_rank(ctx, k, labels, mask);
    };
    z3::expr selected_rank =
        rank_for_mask(selected_mask).simplify();
    z3::expr demand = ctx.int_val(0);
    std::vector<z3::expr> ranks;
    if (target_mode == 0) {
        ranks.reserve(negative_masks.size());
        for (const unsigned int mask : negative_masks) {
            const z3::expr rank =
                mask == selected_mask
                ? selected_rank
                : rank_for_mask(mask).simplify();
            ranks.push_back(rank);
            demand = demand + rank;
        }
    }
    solver.add(selected_rank > 0);
    if (selected_rank_value > 0) {
        solver.add(selected_rank == selected_rank_value);
    }
    if (target_mode == 0) {
        for (const z3::expr& rank : ranks) {
            solver.add(rank <= selected_rank);
        }
    }

    const z3::expr local = residual_selected_local(
        ctx,
        k,
        labels,
        orbit,
        selected_mask,
        level_parity
    ).simplify();
    z3::expr supply = local;
    if (target_mode == 0 || target_mode == 1) {
        supply = (
            supply + residual_pair_reservoir(
            ctx, k, labels, orbit, level_parity, equality_mask
            )
        ).simplify();
    }
    if (target_mode == 0) {
        supply = (
            supply + residual_positive_reservoir(
            ctx, k, labels, orbit
            )
        ).simplify();
        solver.add(supply < demand);
    } else {
        solver.add(
            supply
            < static_cast<int>(negative_masks.size()) * selected_rank
        );
    }

    const char* dump_path = std::getenv("SU2_SEVEN_DUMP_SMT2");
    if (dump_path != nullptr && dump_path[0] != '\0') {
        std::ofstream dump(dump_path);
        if (!dump) {
            return {
                false,
                std::string("cannot open SMT2 dump path: ") + dump_path
            };
        }
        dump << solver.to_smt2();
    }

    const z3::check_result result = solver.check();
    if (result == z3::unsat) {
        return {true, {}};
    }
    std::ostringstream diagnostic;
    diagnostic << "orbit=" << orbit_index
               << " m=" << orbit.minus_count
               << " o=" << orbit.odd_count
               << " r=" << orbit.odd_minus_count
               << " level_parity=" << level_parity
               << " selected_orbit=" << selected_orbit
               << " selected_mask=" << selected_mask
               << " result=" << result;
    if (result == z3::sat) {
        diagnostic << " model=" << solver.get_model();
    } else {
        diagnostic << " reason=" << solver.reason_unknown();
    }
    return {false, diagnostic.str()};
}

struct ResidualTask {
    int orbit_index;
    int level_parity;
    int selected_orbit;
};

std::vector<ResidualTask> residual_tasks() {
    std::vector<ResidualTask> result;
    for (int orbit_index = 0;
         orbit_index < static_cast<int>(residual_orbits.size());
         ++orbit_index) {
        const int selected_count = static_cast<int>(
            selected_orbit_masks(
                residual_orbits[
                    static_cast<std::size_t>(orbit_index)
                ]
            ).size()
        );
        for (int level_parity = 0;
             level_parity < 2; ++level_parity) {
            for (int selected_orbit = 0;
                 selected_orbit < selected_count;
                 ++selected_orbit) {
                result.push_back(ResidualTask{
                    orbit_index, level_parity, selected_orbit
                });
            }
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    const std::vector<ResidualTask> tasks = residual_tasks();
    if (argc >= 4 && argc <= 10) {
        const int orbit_index = std::atoi(argv[1]);
        const int level_parity = std::atoi(argv[2]);
        const int selected_orbit = std::atoi(argv[3]);
        const int cap_mask =
            argc >= 5 ? std::atoi(argv[4]) : -1;
        const int selected_wall_mask =
            argc >= 6 ? std::atoi(argv[5]) : -1;
        const int selected_interval_mask =
            argc >= 7 ? std::atoi(argv[6]) : -1;
        const int equality_mask =
            argc >= 8 ? std::atoi(argv[7]) : -1;
        const int target_mode =
            argc >= 9 ? std::atoi(argv[8]) : 0;
        const int selected_rank_value =
            argc == 10 ? std::atoi(argv[9]) : 0;
        if (orbit_index < 0
            || orbit_index >= static_cast<int>(
                residual_orbits.size()
            )
            || level_parity < 0 || level_parity >= 2
            || selected_orbit < 0
            || selected_orbit >= static_cast<int>(
                selected_orbit_masks(
                    residual_orbits[
                        static_cast<std::size_t>(orbit_index)
                    ]
                ).size()
            )
            || cap_mask < -1 || cap_mask >= (1 << 7)
            || selected_wall_mask < -1
            || selected_wall_mask >= (1 << 3)
            || selected_interval_mask < -1
            || selected_interval_mask >= (1 << 4)
            || equality_mask < -1
            || equality_mask >= (1 << 15)
            || target_mode < 0 || target_mode > 2
            || selected_rank_value < 0) {
            std::cerr
                << "usage: verify_su2_seven_residual_z3 "
                << "[ORBIT(0..7) LEVEL_PARITY(0..1) "
                << "SELECTED_ORBIT [CAP_MASK(0..127) "
                << "SELECTED_WALL_MASK(0..7) "
                << "SELECTED_INTERVAL_MASK(0..15) "
                << "EQUALITY_MASK(0..32767) "
                << "TARGET_MODE(0=direct,1=deep-pair-ceiling,"
                << "2=deep-local-ceiling) "
                << "SELECTED_RANK(0=unbounded)]]]]]]\n";
            return EXIT_FAILURE;
        }
        const QueryResult result = verify_residual_query(
            orbit_index,
            level_parity,
            selected_orbit,
            cap_mask,
            selected_wall_mask,
            selected_interval_mask,
            equality_mask,
            target_mode,
            selected_rank_value
        );
        if (!result.passed) {
            std::cerr << "SU2_SEVEN_RESIDUAL_Z3 FAIL "
                      << result.diagnostic << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "SU2_SEVEN_RESIDUAL_Z3"
                  << " orbit=" << orbit_index
                  << " level_parity=" << level_parity
                  << " selected_orbit=" << selected_orbit
                  << " cap_mask=" << cap_mask
                  << " selected_wall_mask=" << selected_wall_mask
                  << " selected_interval_mask="
                  << selected_interval_mask
                  << " equality_mask=" << equality_mask
                  << " target_mode=" << target_mode
                  << " selected_rank_value="
                  << selected_rank_value
                  << " counterexample=UNSAT result=PASS\n";
        return EXIT_SUCCESS;
    }
    if (argc != 1) {
        std::cerr
            << "usage: verify_su2_seven_residual_z3 "
            << "[ORBIT(0..7) LEVEL_PARITY(0..1) "
            << "SELECTED_ORBIT [CAP_MASK(0..127) "
            << "SELECTED_WALL_MASK(0..7) "
            << "SELECTED_INTERVAL_MASK(0..15) "
            << "EQUALITY_MASK(0..32767) "
            << "TARGET_MODE(0=direct,1=deep-pair-ceiling,"
            << "2=deep-local-ceiling) "
            << "SELECTED_RANK(0=unbounded)]]]]]]\n";
        return EXIT_FAILURE;
    }

    const unsigned workers = std::min<unsigned>(
        worker_limit(),
        static_cast<unsigned>(tasks.size())
    );
    std::atomic<std::size_t> next_task{0U};
    std::atomic<std::size_t> completed{0U};
    std::atomic<bool> failed{false};
    std::mutex diagnostic_mutex;
    std::string diagnostic;
    std::vector<std::jthread> pool;
    pool.reserve(workers);
    for (unsigned worker = 0; worker < workers; ++worker) {
        pool.emplace_back([&]() {
            while (!failed.load()) {
                const std::size_t task_index =
                    next_task.fetch_add(1U);
                if (task_index >= tasks.size()) {
                    return;
                }
                const ResidualTask& task = tasks[task_index];
                const QueryResult result = verify_residual_query(
                    task.orbit_index,
                    task.level_parity,
                    task.selected_orbit
                );
                if (!result.passed) {
                    {
                        std::lock_guard<std::mutex> lock(
                            diagnostic_mutex
                        );
                        if (diagnostic.empty()) {
                            diagnostic = result.diagnostic;
                        }
                    }
                    failed.store(true);
                    return;
                }
                completed.fetch_add(1U);
            }
        });
    }
    pool.clear();

    if (failed.load() || completed.load() != tasks.size()) {
        std::cerr << "SU2_SEVEN_RESIDUAL_Z3 FAIL"
                  << " completed=" << completed.load()
                  << " queries=" << tasks.size()
                  << " diagnostic=" << diagnostic << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "SU2_SEVEN_RESIDUAL_Z3"
              << " orbits=8 queries=" << tasks.size()
              << " workers=" << workers
              << " counterexamples=UNSAT result=PASS\n";
    return EXIT_SUCCESS;
}
