#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#define SU2_SEVEN_RESIDUAL_NO_MAIN
#include "verify_su2_seven_residual_z3.cpp"
#undef SU2_SEVEN_RESIDUAL_NO_MAIN

namespace {

struct ShallowTask {
    int orbit_index;
    int level_parity;
    int selected_orbit;
    int shallow_position;
    int shallow_kind;
    int selected_rank;
    bool selected_complement_deep;
    int endpoint_pattern;
    int selected_wall_mask;
    int selected_interval_mask;
    int selected_local_value = -1;
};

constexpr int choose_small(int n, int r) {
    if (r < 0 || r > n) {
        return 0;
    }
    if (r == 0 || r == n) {
        return 1;
    }
    int result = 1;
    for (int index = 1; index <= r; ++index) {
        result = result * (n - r + index) / index;
    }
    return result;
}

constexpr int shallow_cut_count(
    int minus_count,
    int shallow_count,
    int shallow_minus_count
) {
    int result = 0;
    for (int additional_minus = 0;
         additional_minus <= 3 - shallow_count;
         ++additional_minus) {
        if (((shallow_minus_count + additional_minus) & 1) == 0) {
            continue;
        }
        result += choose_small(
            minus_count - shallow_minus_count,
            additional_minus
        ) * choose_small(
            7 - minus_count - shallow_count + shallow_minus_count,
            3 - shallow_count - additional_minus
        );
    }
    return result;
}

static_assert(shallow_cut_count(4, 1, 0) == 8);
static_assert(shallow_cut_count(4, 1, 1) == 6);
static_assert(shallow_cut_count(4, 2, 0) == 4);
static_assert(shallow_cut_count(4, 2, 1) == 2);
static_assert(shallow_cut_count(4, 2, 2) == 2);
static_assert(shallow_cut_count(4, 3, 1) == 1);
static_assert(shallow_cut_count(4, 3, 3) == 1);
static_assert(shallow_cut_count(6, 1, 1) == 10);
static_assert(shallow_cut_count(6, 2, 2) == 4);
static_assert(shallow_cut_count(6, 3, 3) == 1);

bool selected_cut_can_be_active_shallow(
    const ResidualOrbit& orbit,
    unsigned int mask
) {
    int parity = 0;
    for (int index = 0; index < 7; ++index) {
        if (((mask >> index) & 1U) != 0U) {
            parity ^= orbit_label_parity(orbit, index);
        }
    }
    return parity == 0;
}

int shallow_kind_parity(int kind, int level_parity) {
    if (kind == 0) {
        return 1;
    }
    if (kind == 1) {
        return 1 - level_parity;
    }
    return level_parity;
}

int selected_pattern_threshold(
    const ResidualOrbit& orbit,
    int shallow_count,
    int shallow_minus_count,
    int sole_kind
) {
    if (orbit.minus_count == 4) {
        if (shallow_count == 1) {
            if (sole_kind != 2) {
                return shallow_minus_count == 0 ? 5 : 4;
            }
            return shallow_minus_count == 0 ? 28 : 11;
        }
        if (shallow_count == 2) {
            return 5;
        }
        return 11;
    }
    if (shallow_count == 1) {
        return sole_kind == 2 ? 9 : 7;
    }
    if (shallow_count == 2) {
        return 6;
    }
    return 14;
}

constexpr int small_popcount(unsigned int value) {
    int result = 0;
    while (value != 0U) {
        result += static_cast<int>(value & 1U);
        value >>= 1U;
    }
    return result;
}

constexpr int parity_active_negative_count(
    const ResidualOrbit& orbit,
    unsigned int required_mask
) {
    int result = 0;
    const unsigned int minus_mask =
        (1U << orbit.minus_count) - 1U;
    for (unsigned int mask = 0U; mask < (1U << 7U); ++mask) {
        if (small_popcount(mask) != 3
            || (mask & required_mask) != required_mask
            || (small_popcount(mask & minus_mask) & 1) == 0) {
            continue;
        }
        int parity = 0;
        for (int position = 0; position < 7; ++position) {
            if (((mask >> position) & 1U) != 0U) {
                parity ^= orbit_label_parity(orbit, position);
            }
        }
        if (parity == 0) {
            ++result;
        }
    }
    return result;
}

constexpr int parity_low_cut_constant(
    const ResidualOrbit& orbit,
    unsigned int shallow_mask,
    unsigned int cap_one_mask
) {
    int result = 0;
    const unsigned int minus_mask =
        (1U << orbit.minus_count) - 1U;
    for (unsigned int mask = 0U; mask < (1U << 7U); ++mask) {
        if (small_popcount(mask) != 3
            || (small_popcount(mask & minus_mask) & 1) == 0) {
            continue;
        }
        int parity = 0;
        for (int position = 0; position < 7; ++position) {
            if (((mask >> position) & 1U) != 0U) {
                parity ^= orbit_label_parity(orbit, position);
            }
        }
        if (parity != 0 || (mask & shallow_mask) == shallow_mask) {
            continue;
        }
        const unsigned int omitted = shallow_mask & ~mask;
        result += (omitted & cap_one_mask) != 0U ? 1 : 2;
    }
    return result;
}

int parity_refined_threshold(
    const ResidualOrbit& orbit,
    int shallow_count,
    int shallow_minus_count,
    int sole_kind,
    unsigned int shallow_mask,
    unsigned int cap_one_mask,
    bool top_depth_at_least_four
) {
    const int coarse = selected_pattern_threshold(
        orbit,
        shallow_count,
        shallow_minus_count,
        sole_kind
    );
    int slope = 0;
    int intercept = 0;
    if (shallow_count == 1) {
        if (sole_kind == 2) {
            if (top_depth_at_least_four) {
                slope = 16;
                intercept = 31;
            } else {
                slope = 9;
                intercept = 12;
            }
        } else {
            slope = 17;
            intercept = 24;
        }
    } else if (shallow_count == 2) {
        slope = 12;
        intercept = 15;
    } else {
        slope = 4;
        intercept = 3;
    }
    const int containing =
        parity_active_negative_count(orbit, shallow_mask);
    if (slope <= containing) {
        return coarse;
    }
    const int numerator = intercept + parity_low_cut_constant(
        orbit, shallow_mask, cap_one_mask
    );
    const int denominator = slope - containing;
    const int refined =
        (numerator + denominator - 1) / denominator;
    return std::min(coarse, refined);
}

bool closes_by_sharp_small_rank_profile(
    const ResidualOrbit& orbit,
    int shallow_count,
    int sole_kind,
    unsigned int shallow_mask,
    unsigned int cap_one_mask,
    bool top_depth_at_least_four,
    int rank
) {
    int local = -1;
    if (rank == 3) {
        if (shallow_count == 1) {
            if (sole_kind != 2) {
                local = 31;
            } else {
                local = top_depth_at_least_four ? 20 : 17;
            }
        } else if (shallow_count == 2) {
            local = 24;
        } else {
            local = 10;
        }
    } else if (rank == 4) {
        if (shallow_count == 1) {
            if (sole_kind != 2) {
                local = 48;
            } else {
                local = top_depth_at_least_four ? 36 : 26;
            }
        } else if (shallow_count == 2) {
            local = 36;
        } else {
            local = 14;
        }
    }
    if (local < 0) {
        return false;
    }
    const int containing =
        parity_active_negative_count(orbit, shallow_mask);
    const int ceiling = containing * rank
        + parity_low_cut_constant(
            orbit, shallow_mask, cap_one_mask
        );
    return local >= ceiling;
}

constexpr bool closes_by_three_endpoint_graph(
    const ResidualOrbit& orbit,
    int level_parity,
    int shallow_count,
    int pattern
) {
    if (orbit.minus_count != 4 || shallow_count != 3) {
        return false;
    }
    int one_position = -1;
    int top_neighbor_position = -1;
    int top_position = -1;
    for (int position = 0; position < 7; ++position) {
        const int code = (pattern >> (2 * position)) & 3;
        if (code == 1) {
            one_position = position;
        } else if (code == 2) {
            top_neighbor_position = position;
        } else if (code == 3) {
            top_position = position;
        }
    }
    if (orbit.odd_count == 2 && orbit.odd_minus_count == 0) {
        if (level_parity == 0) {
            return one_position >= orbit.minus_count
                && top_neighbor_position >= orbit.minus_count
                && top_position >= 0
                && top_position < orbit.minus_count;
        }
        return one_position >= orbit.minus_count
            && top_position >= orbit.minus_count
            && top_neighbor_position >= 0
            && top_neighbor_position < orbit.minus_count;
    }
    if (orbit.odd_count == 2
        && orbit.odd_minus_count == 1
        && level_parity == 0) {
        const bool one_is_minus =
            one_position >= 0
            && one_position < orbit.minus_count;
        const bool neighbor_is_minus =
            top_neighbor_position >= 0
            && top_neighbor_position < orbit.minus_count;
        return one_is_minus != neighbor_is_minus
            && top_position >= orbit.minus_count;
    }
    if (orbit.odd_count == 4
        && orbit.odd_minus_count == 1
        && level_parity == 0) {
        return one_position >= orbit.minus_count
            && top_neighbor_position >= orbit.minus_count
            && top_position >= 0
            && top_position < orbit.minus_count;
    }
    if (orbit.odd_count == 4
        && orbit.odd_minus_count == 1
        && level_parity == 1) {
        return one_position >= orbit.minus_count
            && top_position >= orbit.minus_count
            && top_neighbor_position >= 0
            && top_neighbor_position < orbit.minus_count;
    }
    if (orbit.odd_count == 4
        && orbit.odd_minus_count == 2) {
        const bool all_minus =
            one_position >= 0
            && one_position < orbit.minus_count
            && top_neighbor_position >= 0
            && top_neighbor_position < orbit.minus_count
            && top_position >= 0
            && top_position < orbit.minus_count;
        const bool one_minus_two_plus =
            one_position >= orbit.minus_count
            && top_neighbor_position >= orbit.minus_count
            && top_position >= 0
            && top_position < orbit.minus_count;
        const bool mixed_one_minus_two_plus =
            top_position >= orbit.minus_count
            && (
                (one_position < orbit.minus_count)
                != (
                    top_neighbor_position
                        < orbit.minus_count
                )
            );
        return all_minus
            || one_minus_two_plus
            || mixed_one_minus_two_plus;
    }
    if (orbit.odd_count == 6
        && orbit.odd_minus_count == 3
        && level_parity == 0) {
        const bool one_is_minus =
            one_position >= 0
            && one_position < orbit.minus_count;
        const bool neighbor_is_minus =
            top_neighbor_position >= 0
            && top_neighbor_position < orbit.minus_count;
        return one_is_minus == neighbor_is_minus
            && top_position >= 0
            && top_position < orbit.minus_count;
    }
    return false;
}

constexpr bool closes_by_single_top_graph(
    int orbit_index,
    int level_parity,
    int selected_orbit,
    int shallow_count,
    int sole_kind
) {
    if (orbit_index != 3
        && orbit_index != 6
        && orbit_index != 7) {
        return false;
    }
    if (shallow_count != 1
        || sole_kind != 2) {
        return false;
    }
    if (orbit_index == 3) {
        return (level_parity == 1 && selected_orbit == 0)
            || (level_parity == 0 && selected_orbit == 1);
    }
    if (orbit_index == 6) {
        return level_parity == 0 && selected_orbit == 0;
    }
    return level_parity == 0 && selected_orbit == 0;
}

constexpr bool closes_by_single_neighbor_graph(
    int orbit_index,
    int level_parity,
    int selected_orbit,
    int shallow_count,
    int sole_kind,
    int rank
) {
    return orbit_index == 6
        && level_parity == 1
        && selected_orbit == 0
        && (
            shallow_count == 1
            || shallow_count == 2
        )
        && sole_kind == 1
        && rank >= 3;
}

static_assert(closes_by_three_endpoint_graph(
    residual_orbits[2], 0, 3, 3330
));
static_assert(closes_by_three_endpoint_graph(
    residual_orbits[2], 0, 3, 3585
));

void append_selected_pattern_tasks(
    std::vector<ShallowTask>& result,
    int orbit_index,
    int level_parity,
    int selected_orbit,
    int maximum_rank,
    const ResidualOrbit& orbit,
    unsigned int selected_mask
) {
    const std::array<int, 3> cut = mask_triple(selected_mask);
    for (int first_code = 0; first_code < 4; ++first_code) {
        for (int second_code = 0; second_code < 4; ++second_code) {
            for (int third_code = 0; third_code < 4; ++third_code) {
                const std::array<int, 3> codes{
                    first_code, second_code, third_code
                };
                int shallow_count = 0;
                int shallow_minus_count = 0;
                int sole_kind = -1;
                int shallow_kind_mask = 0;
                int pattern = 0;
                unsigned int shallow_mask = 0U;
                unsigned int cap_one_mask = 0U;
                bool compatible = true;
                for (std::size_t offset = 0U;
                     offset < cut.size(); ++offset) {
                    const int position = cut[offset];
                    const int code = codes[offset];
                    for (std::size_t earlier = 0U;
                         earlier < offset; ++earlier) {
                        if (orbit_label_class(
                                orbit, cut[earlier]
                            ) == orbit_label_class(orbit, position)
                            && codes[earlier] > code) {
                            compatible = false;
                        }
                    }
                    if (!compatible) {
                        break;
                    }
                    if (code == 0) {
                        continue;
                    }
                    const int kind = code - 1;
                    if (shallow_kind_parity(
                            kind, level_parity
                        ) != orbit_label_parity(orbit, position)) {
                        compatible = false;
                        break;
                    }
                    ++shallow_count;
                    shallow_kind_mask |= 1 << kind;
                    if (position < orbit.minus_count) {
                        ++shallow_minus_count;
                    }
                    sole_kind = kind;
                    pattern |= code << (2 * position);
                    shallow_mask |= 1U << position;
                    if (kind == 2) {
                        cap_one_mask |= 1U << position;
                    }
                }
                if (!compatible || shallow_count == 0) {
                    continue;
                }
                if (shallow_count == 3
                    && shallow_kind_mask != 7) {
                    continue;
                }
                if (shallow_count == 2
                    && cap_one_mask != 0U) {
                    continue;
                }
                bool top_depth_at_least_four = false;
                if (shallow_count == 1
                    && sole_kind == 2
                    && level_parity == 0) {
                    top_depth_at_least_four = true;
                    for (std::size_t offset = 0U;
                         offset < cut.size(); ++offset) {
                        const int position = cut[offset];
                        if (((shallow_mask >> position) & 1U) != 0U) {
                            continue;
                        }
                        top_depth_at_least_four =
                            top_depth_at_least_four
                            && orbit_label_parity(
                                orbit, position
                            ) == 1;
                        }
                }
                if (closes_by_three_endpoint_graph(
                        orbit,
                        level_parity,
                        shallow_count,
                        pattern
                    )) {
                    continue;
                }
                if (closes_by_single_top_graph(
                        orbit_index,
                        level_parity,
                        selected_orbit,
                        shallow_count,
                        sole_kind
                    )) {
                    continue;
                }
                const int threshold = parity_refined_threshold(
                    orbit,
                    shallow_count,
                    shallow_minus_count,
                    sole_kind,
                    shallow_mask,
                    cap_one_mask,
                    top_depth_at_least_four
                );
                const int last_rank =
                    std::min(maximum_rank, threshold - 1);
                for (int rank = 3; rank <= last_rank; ++rank) {
                    if (closes_by_single_neighbor_graph(
                            orbit_index,
                            level_parity,
                            selected_orbit,
                            shallow_count,
                            sole_kind,
                            rank
                        )) {
                        continue;
                    }
                    if (closes_by_sharp_small_rank_profile(
                            orbit,
                            shallow_count,
                            sole_kind,
                            shallow_mask,
                            cap_one_mask,
                            top_depth_at_least_four,
                            rank
                        )) {
                        continue;
                    }
                    result.push_back(ShallowTask{
                        orbit_index,
                        level_parity,
                        selected_orbit,
                        -1,
                        -1,
                        rank,
                        true,
                        pattern,
                        -1,
                        -1
                    });
                }
            }
        }
    }
}

std::vector<ShallowTask> shallow_tasks(
    int maximum_rank,
    bool patterns_only
) {
    std::vector<ShallowTask> result;
    for (int orbit_index = 1;
         orbit_index < static_cast<int>(residual_orbits.size());
         ++orbit_index) {
        const ResidualOrbit& orbit =
            residual_orbits[static_cast<std::size_t>(orbit_index)];
        const std::vector<unsigned int> representatives =
            selected_orbit_masks(orbit);
        for (int selected_orbit = 0;
             selected_orbit
                 < static_cast<int>(representatives.size());
             ++selected_orbit) {
            if (!selected_cut_can_be_active_shallow(
                    orbit,
                    representatives[
                        static_cast<std::size_t>(selected_orbit)
                    ])) {
                continue;
            }
            for (int level_parity = 0;
                 level_parity < 2; ++level_parity) {
                if (!patterns_only) {
                    for (int position = 0; position < 7; ++position) {
                        bool canonical_position = true;
                        const unsigned int selected_mask =
                            representatives[
                                static_cast<std::size_t>(
                                    selected_orbit
                                )
                            ];
                        const bool position_is_selected =
                            ((selected_mask >> position) & 1U) != 0U;
                        for (int earlier = 0;
                             earlier < position; ++earlier) {
                            const bool earlier_is_selected =
                                ((selected_mask >> earlier) & 1U) != 0U;
                            if (earlier_is_selected
                                    == position_is_selected
                                && orbit_label_class(orbit, earlier)
                                    == orbit_label_class(
                                        orbit, position
                                    )) {
                                canonical_position = false;
                                break;
                            }
                        }
                        if (!canonical_position) {
                            continue;
                        }
                        const int required_parity =
                            orbit_label_parity(orbit, position);
                        for (int kind = 0; kind < 3; ++kind) {
                            if (shallow_kind_parity(
                                    kind, level_parity
                                ) != required_parity) {
                                continue;
                            }
                            std::vector<int> ranks;
                            if (maximum_rank == 0) {
                                ranks.push_back(0);
                            } else {
                                ranks.push_back(1);
                                ranks.push_back(2);
                            }
                            for (const int rank : ranks) {
                                if (
                                    rank == 1
                                    && orbit_index == 1
                                    && selected_orbit == 1
                                    && (
                                        (
                                            level_parity == 0
                                            && (
                                                (
                                                    position == 0
                                                    && kind == 2
                                                )
                                                || (
                                                    position == 1
                                                    && kind == 2
                                                )
                                                || (
                                                    position == 4
                                                    && (
                                                        kind == 0
                                                        || kind == 1
                                                    )
                                                )
                                                || (
                                                    position == 6
                                                    && kind == 2
                                                )
                                            )
                                        )
                                        || (
                                            level_parity == 1
                                            && (
                                                (
                                                    position == 0
                                                    && kind == 1
                                                )
                                                || (
                                                    position == 1
                                                    && kind == 1
                                                )
                                                || (
                                                    position == 4
                                                    && (
                                                        kind == 0
                                                        || kind == 2
                                                    )
                                                )
                                                || (
                                                    position == 6
                                                    && kind == 1
                                                )
                                            )
                                        )
                                    )
                                ) {
                                    continue;
                                }
                                if (
                                    rank == 1
                                    && orbit_index == 1
                                    && level_parity == 0
                                    && selected_orbit == 2
                                    && (
                                        position == 0
                                        || position == 3
                                    )
                                    && kind == 2
                                ) {
                                    continue;
                                }
                                if (rank == 2
                                    && orbit_index == 6
                                    && level_parity == 0
                                    && kind == 2
                                    && position_is_selected) {
                                    continue;
                                }
                                if (rank == 2
                                    && orbit_index == 6
                                    && level_parity == 1
                                    && kind == 1
                                    && position_is_selected) {
                                    continue;
                                }
                                if (rank == 2
                                    && orbit_index == 6
                                    && level_parity == 1
                                    && kind == 1
                                    && !position_is_selected
                                    && position < orbit.minus_count) {
                                    continue;
                                }
                                if (rank == 2
                                    && orbit_index == 6
                                    && level_parity == 1
                                    && kind == 1
                                    && !position_is_selected
                                    && position >= orbit.minus_count) {
                                    continue;
                                }
                                if (rank == 2
                                    && kind == 2
                                    && !position_is_selected) {
                                    continue;
                                }
                                result.push_back(ShallowTask{
                                    orbit_index,
                                    level_parity,
                                    selected_orbit,
                                    position,
                                    kind,
                                    rank,
                                    false,
                                    -1,
                                    -1,
                                    -1
                                });
                            }
                        }
                    }
                }
                if (maximum_rank > 2) {
                    append_selected_pattern_tasks(
                        result,
                        orbit_index,
                        level_parity,
                        selected_orbit,
                        maximum_rank,
                        orbit,
                        representatives[
                            static_cast<std::size_t>(selected_orbit)
                        ]
                    );
                }
            }
        }
    }
    return result;
}

std::string task_name(const ShallowTask& task) {
    return "orbit=" + std::to_string(task.orbit_index)
        + " parity=" + std::to_string(task.level_parity)
        + " selected=" + std::to_string(task.selected_orbit)
        + " position=" + std::to_string(task.shallow_position)
        + " kind=" + std::to_string(task.shallow_kind)
        + " rank=" + std::to_string(task.selected_rank)
        + " pattern=" + std::to_string(task.endpoint_pattern)
        + " wall=" + std::to_string(task.selected_wall_mask)
        + " interval=" + std::to_string(
            task.selected_interval_mask
        )
        + " local=" + std::to_string(task.selected_local_value);
}

std::vector<ShallowTask> rank_one_switch_cover(
    const std::vector<ShallowTask>& base_tasks
) {
    std::vector<ShallowTask> result;
    for (const ShallowTask& base : base_tasks) {
        if (base.selected_rank != 1) {
            continue;
        }
        const ResidualOrbit& orbit =
            residual_orbits[
                static_cast<std::size_t>(base.orbit_index)
            ];
        const int last_local =
            (orbit.minus_count == 4 ? 16 : 20) - 1;
        for (int local = 1; local <= last_local; ++local) {
            if (local >= 5) {
                ShallowTask task = base;
                task.selected_local_value = local;
                result.push_back(task);
                continue;
            }
            for (int wall = 0; wall < 8; ++wall) {
                ShallowTask task = base;
                task.selected_wall_mask = wall;
                task.selected_local_value = local;
                result.push_back(task);
            }
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    const bool show_pattern =
        argc == 3 && std::string(argv[1]) == "--show-pattern";
    const bool pattern_chambers =
        argc == 3 && std::string(argv[1]) == "--pattern-chambers";
    const bool small_task =
        argc == 3 && std::string(argv[1]) == "--small-task";
    const bool small_supply_task =
        argc == 3
        && std::string(argv[1]) == "--small-supply-task";
    const bool small_supply_chambers =
        argc == 3
        && std::string(argv[1]) == "--small-supply-chambers";
    const bool small_switch =
        argc == 2 && std::string(argv[1]) == "--small-switch";
    const bool small_switch_task =
        argc == 3
        && std::string(argv[1]) == "--small-switch-task";
    const bool small_switch_bad_bound =
        argc == 4
        && std::string(argv[1]) == "--small-switch-bad-bound";
    const bool small_switch_bad_bound_descriptor =
        argc == 8
        && std::string(argv[1])
            == "--small-switch-bad-bound-descriptor";
    const bool small_switch_cap =
        argc == 4
        && std::string(argv[1]) == "--small-switch-cap";
    const bool small_switch_cell =
        argc == 6
        && std::string(argv[1]) == "--small-switch-cell";
    const bool small_switch_local =
        argc == 7
        && std::string(argv[1]) == "--small-switch-local";
    const bool small_switch_cover =
        argc == 2
        && std::string(argv[1]) == "--small-switch-cover";
    const bool small_switch_cover_task =
        argc == 3
        && std::string(argv[1]) == "--small-switch-cover-task";
    const bool small_switch_cover_range =
        argc == 4
        && std::string(argv[1]) == "--small-switch-cover-range";
    const bool list_switch_cover =
        argc == 2
        && std::string(argv[1]) == "--list-small-switch-cover";
    const bool small_rank =
        argc == 2 && std::string(argv[1]) == "--small";
    const bool bounded_rank =
        argc == 2 && std::string(argv[1]) == "--bounded";
    const bool patterns_only =
        argc == 2 && std::string(argv[1]) == "--patterns";
    const bool count =
        argc == 2 && (
            std::string(argv[1]) == "--count"
            || std::string(argv[1]) == "--count-small"
            || std::string(argv[1]) == "--count-bounded"
            || std::string(argv[1]) == "--count-patterns"
            || std::string(argv[1])
                == "--count-small-switch-cover"
        );
    const bool count_small =
        argc == 2 && std::string(argv[1]) == "--count-small";
    const bool count_bounded =
        argc == 2 && std::string(argv[1]) == "--count-bounded";
    const bool count_patterns =
        argc == 2 && std::string(argv[1]) == "--count-patterns";
    const bool count_switch_cover =
        argc == 2
        && std::string(argv[1]) == "--count-small-switch-cover";
    const bool list_small =
        argc == 2 && std::string(argv[1]) == "--list-small";
    if (argc != 1
        && !small_rank
        && !bounded_rank
        && !patterns_only
        && !count
        && !list_small
        && !small_task
        && !small_supply_task
        && !small_supply_chambers
        && !small_switch
        && !small_switch_task
        && !small_switch_bad_bound
        && !small_switch_bad_bound_descriptor
        && !small_switch_cap
        && !small_switch_cell
        && !small_switch_local
        && !small_switch_cover
        && !small_switch_cover_task
        && !small_switch_cover_range
        && !list_switch_cover
        && !show_pattern
        && !pattern_chambers) {
        std::cerr
            << "usage: verify_su2_seven_shallow_z3 "
            << "[--count|--count-small|--count-bounded"
            << "|--count-patterns|--small|--bounded|--patterns"
            << "|--list-small"
            << "|--small-task INDEX"
            << "|--small-supply-task INDEX"
            << "|--small-supply-chambers INDEX"
            << "|--small-switch|--small-switch-task INDEX"
            << "|--small-switch-bad-bound INDEX THRESHOLD"
            << "|--small-switch-bad-bound-descriptor"
            << " ORBIT PARITY SELECTED POSITION KIND THRESHOLD"
            << "|--small-switch-cap INDEX CAP_MASK"
            << "|--small-switch-cell INDEX CAP_MASK WALL INTERVAL"
            << "|--small-switch-local INDEX CAP_MASK WALL INTERVAL LOCAL"
            << "|--small-switch-cover"
            << "|--count-small-switch-cover"
            << "|--small-switch-cover-task INDEX"
            << "|--small-switch-cover-range BEGIN END"
            << "|--list-small-switch-cover"
            << "|--show-pattern INDEX|--pattern-chambers INDEX]\n";
        return EXIT_FAILURE;
    }
    const int maximum_rank =
        bounded_rank || count_bounded || patterns_only || count_patterns
            || show_pattern || pattern_chambers
            ? 27
            : (
                small_rank || count_small || list_small || small_task
                    || small_supply_task
                    || small_supply_chambers
                    || small_switch || small_switch_task
                    || small_switch_bad_bound
                    || small_switch_bad_bound_descriptor
                    || small_switch_cap || small_switch_cell
                    || small_switch_local
                    || small_switch_cover
                    || small_switch_cover_task
                    || small_switch_cover_range
                    || list_switch_cover
                    || count_switch_cover
                    ? 2
                    : 0
            );
    const bool use_patterns_only =
        patterns_only || count_patterns || show_pattern
        || pattern_chambers;
    std::vector<ShallowTask> tasks =
        shallow_tasks(maximum_rank, use_patterns_only);
    if (
        small_switch_cover || small_switch_cover_task
        || small_switch_cover_range || list_switch_cover
        || count_switch_cover
    ) {
        tasks = rank_one_switch_cover(tasks);
    }
    int switch_cap_mask = -1;
    int switch_local_value = -1;
    int switch_bad_threshold = -1;
    if (small_switch_bad_bound_descriptor) {
        std::array<long, 6> values{};
        for (int index = 0; index < 6; ++index) {
            char* end = nullptr;
            values[static_cast<std::size_t>(index)] =
                std::strtol(argv[index + 2], &end, 10);
            if (
                end == argv[index + 2]
                || *end != '\0'
            ) {
                std::cerr << "invalid switch descriptor\n";
                return EXIT_FAILURE;
            }
        }
        const int orbit_index = static_cast<int>(values[0]);
        const int level_parity = static_cast<int>(values[1]);
        const int selected_orbit = static_cast<int>(values[2]);
        const int shallow_position = static_cast<int>(values[3]);
        const int shallow_kind = static_cast<int>(values[4]);
        const int threshold = static_cast<int>(values[5]);
        if (
            orbit_index < 0
            || orbit_index
                >= static_cast<int>(residual_orbits.size())
            || level_parity < 0 || level_parity > 1
            || shallow_position < 0 || shallow_position >= 7
            || shallow_kind < 0 || shallow_kind > 2
            || threshold < 1 || threshold > 20
        ) {
            std::cerr << "invalid switch descriptor\n";
            return EXIT_FAILURE;
        }
        const ResidualOrbit& orbit =
            residual_orbits[static_cast<std::size_t>(orbit_index)];
        const std::vector<unsigned int> representatives =
            selected_orbit_masks(orbit);
        if (
            selected_orbit < 0
            || selected_orbit
                >= static_cast<int>(representatives.size())
            || !selected_cut_can_be_active_shallow(
                orbit,
                representatives[
                    static_cast<std::size_t>(selected_orbit)
                ]
            )
            || shallow_kind_parity(shallow_kind, level_parity)
                != orbit_label_parity(orbit, shallow_position)
        ) {
            std::cerr << "invalid switch descriptor\n";
            return EXIT_FAILURE;
        }
        tasks.assign(1U, ShallowTask{
            orbit_index,
            level_parity,
            selected_orbit,
            shallow_position,
            shallow_kind,
            1,
            false,
            -1,
            -1,
            -1
        });
        switch_bad_threshold = threshold;
    }
    if (list_switch_cover) {
        for (std::size_t index = 0U; index < tasks.size(); ++index) {
            std::cout << index << ' ' << task_name(tasks[index]) << '\n';
        }
        return EXIT_SUCCESS;
    }
    if (small_switch_cover_range) {
        char* begin_end = nullptr;
        char* finish_end = nullptr;
        const long begin = std::strtol(argv[2], &begin_end, 10);
        const long finish = std::strtol(argv[3], &finish_end, 10);
        if (
            begin_end == argv[2] || *begin_end != '\0'
            || finish_end == argv[3] || *finish_end != '\0'
            || begin < 0 || finish <= begin
            || static_cast<std::size_t>(finish) > tasks.size()
        ) {
            std::cerr << "invalid switch-cover range\n";
            return EXIT_FAILURE;
        }
        tasks = std::vector<ShallowTask>(
            tasks.begin() + begin,
            tasks.begin() + finish
        );
    }
    if (list_small) {
        for (std::size_t index = 0U; index < tasks.size(); ++index) {
            std::cout << index << ' ' << task_name(tasks[index]) << '\n';
        }
        return EXIT_SUCCESS;
    }
    if (
        small_task || small_supply_task || small_switch_task
        || small_switch_bad_bound
        || small_switch_cap || small_switch_cell
        || small_switch_local
        || small_switch_cover_task
        || small_supply_chambers
        || show_pattern || pattern_chambers
    ) {
        char* end = nullptr;
        const long parsed = std::strtol(argv[2], &end, 10);
        if (end == argv[2]
            || *end != '\0'
            || parsed < 0
            || static_cast<std::size_t>(parsed) >= tasks.size()) {
            std::cerr << "invalid pattern task index\n";
            return EXIT_FAILURE;
        }
        ShallowTask base =
            tasks[static_cast<std::size_t>(parsed)];
        if (
            (
                small_switch_task || small_switch_bad_bound
                || small_switch_cap
                || small_switch_cell || small_switch_local
                || small_switch_cover_task
            )
            && base.selected_rank != 1
        ) {
            std::cerr << "switch task must have selected rank one\n";
            return EXIT_FAILURE;
        }
        if (
            small_switch_cap || small_switch_cell
            || small_switch_local
        ) {
            char* cap_end = nullptr;
            const long cap_parsed =
                std::strtol(argv[3], &cap_end, 10);
            if (
                cap_end == argv[3]
                || *cap_end != '\0'
                || (
                    small_switch_cap
                    && cap_parsed < 0
                )
                || cap_parsed < -1
                || cap_parsed >= 128
            ) {
                std::cerr << "invalid switch cap mask\n";
                return EXIT_FAILURE;
            }
            switch_cap_mask = static_cast<int>(cap_parsed);
        }
        if (small_switch_bad_bound) {
            char* threshold_end = nullptr;
            const long threshold_parsed =
                std::strtol(argv[3], &threshold_end, 10);
            if (
                threshold_end == argv[3]
                || *threshold_end != '\0'
                || threshold_parsed < 1
                || threshold_parsed > 20
            ) {
                std::cerr << "invalid bad-switch threshold\n";
                return EXIT_FAILURE;
            }
            switch_bad_threshold =
                static_cast<int>(threshold_parsed);
        }
        if (small_switch_cell || small_switch_local) {
            char* wall_end = nullptr;
            char* interval_end = nullptr;
            const long wall_parsed =
                std::strtol(argv[4], &wall_end, 10);
            const long interval_parsed =
                std::strtol(argv[5], &interval_end, 10);
            if (
                wall_end == argv[4] || *wall_end != '\0'
                || (
                    small_switch_cell
                    && wall_parsed < 0
                )
                || wall_parsed < -1 || wall_parsed >= 8
                || interval_end == argv[5]
                || *interval_end != '\0'
                || (
                    small_switch_cell
                    && interval_parsed < 0
                )
                || interval_parsed < -1
                || interval_parsed >= 16
            ) {
                std::cerr << "invalid switch wall/interval cell\n";
                return EXIT_FAILURE;
            }
            base.selected_wall_mask =
                static_cast<int>(wall_parsed);
            base.selected_interval_mask =
                static_cast<int>(interval_parsed);
        }
        if (small_switch_local) {
            char* local_end = nullptr;
            const long local_parsed =
                std::strtol(argv[6], &local_end, 10);
            if (
                local_end == argv[6] || *local_end != '\0'
                || local_parsed < 1 || local_parsed >= 20
            ) {
                std::cerr << "invalid switch local value\n";
                return EXIT_FAILURE;
            }
            switch_local_value = static_cast<int>(local_parsed);
        }
        if (
            small_task || small_supply_task || small_switch_task
            || small_switch_bad_bound
            || small_switch_cap || small_switch_cell
            || small_switch_local
            || small_switch_cover_task
        ) {
            tasks.assign(1U, base);
        }
        if (show_pattern) {
            std::cout << task_name(base) << '\n';
            return EXIT_SUCCESS;
        }
        if (pattern_chambers) {
            tasks.clear();
            tasks.reserve(8U * 16U);
            for (int wall = 0; wall < 8; ++wall) {
                for (int interval = 0; interval < 16; ++interval) {
                    ShallowTask chamber = base;
                    chamber.selected_wall_mask = wall;
                    chamber.selected_interval_mask = interval;
                    tasks.push_back(chamber);
                }
            }
        }
        if (small_supply_chambers) {
            tasks.clear();
            tasks.reserve(8U * 16U);
            for (int wall = 0; wall < 8; ++wall) {
                for (int interval = 0; interval < 16; ++interval) {
                    ShallowTask chamber = base;
                    chamber.selected_wall_mask = wall;
                    chamber.selected_interval_mask = interval;
                    tasks.push_back(chamber);
                }
            }
        }
    }
    if (small_switch) {
        tasks.erase(
            std::remove_if(
                tasks.begin(),
                tasks.end(),
                [](const ShallowTask& task) {
                    return task.selected_rank != 1;
                }
            ),
            tasks.end()
        );
    }
    const std::string tag = small_switch_cover
        || small_switch_cover_task || small_switch_cover_range
        || count_switch_cover
        ? "SU2_SEVEN_SHALLOW_SMALL_SWITCH_COVER_Z3"
        : (small_switch
        ? "SU2_SEVEN_SHALLOW_SMALL_SWITCH_Z3"
        : (
            small_switch_task || small_switch_bad_bound
                || small_switch_bad_bound_descriptor
                || small_switch_cap
                || small_switch_cell || small_switch_local
        ? "SU2_SEVEN_SHALLOW_SMALL_SWITCH_TASK_Z3"
        : (small_supply_chambers
        ? "SU2_SEVEN_SHALLOW_SMALL_SUPPLY_CHAMBERS_Z3"
        : (small_supply_task
        ? "SU2_SEVEN_SHALLOW_SMALL_SUPPLY_TASK_Z3"
        : (small_task
        ? "SU2_SEVEN_SHALLOW_SMALL_TASK_Z3"
        : (pattern_chambers
        ? "SU2_SEVEN_SHALLOW_PATTERN_CHAMBERS_Z3"
        : (use_patterns_only
        ? "SU2_SEVEN_SHALLOW_PATTERNS_Z3"
        : (maximum_rank == 27
            ? "SU2_SEVEN_SHALLOW_BOUNDED_Z3"
        : (maximum_rank == 2
            ? "SU2_SEVEN_SHALLOW_SMALL_Z3"
            : "SU2_SEVEN_SHALLOW_Z3")))))))));
    if (count) {
        std::cout << tag
                  << " tasks=" << tasks.size() << '\n';
        return EXIT_SUCCESS;
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

    std::mutex progress_mutex;
    std::condition_variable progress_condition;
    bool finished = false;
    std::jthread reporter([&]() {
        std::unique_lock<std::mutex> lock(progress_mutex);
        while (!progress_condition.wait_for(
            lock,
            std::chrono::seconds(30),
            [&]() { return finished; }
        )) {
            std::cout
                << tag << " progress="
                << completed.load() << '/' << tasks.size()
                << " assigned=" << std::min(
                    next_task.load(), tasks.size()
                )
                << " failed=" << (failed.load() ? 1 : 0)
                << '\n' << std::flush;
        }
    });

    std::cout << tag
              << " tasks=" << tasks.size()
              << " workers=" << workers
              << " start=1\n" << std::flush;

    std::vector<std::jthread> pool;
    pool.reserve(workers);
    for (unsigned worker = 0; worker < workers; ++worker) {
        pool.emplace_back([&]() {
            while (!failed.load()) {
                const std::size_t index = next_task.fetch_add(1U);
                if (index >= tasks.size()) {
                    return;
                }
                const ShallowTask& task = tasks[index];
                const QueryResult result = verify_residual_query(
                    task.orbit_index,
                    task.level_parity,
                    task.selected_orbit,
                    switch_cap_mask,
                    task.selected_wall_mask,
                    task.selected_interval_mask,
                    -1,
                    (
                        small_supply_task || small_supply_chambers
                            ? 3
                            : (
                                small_switch || small_switch_task
                                    || small_switch_bad_bound
                                    || small_switch_bad_bound_descriptor
                                    || small_switch_cap
                                    || small_switch_cell
                                    || small_switch_local
                                    || small_switch_cover
                                    || small_switch_cover_task
                                    || small_switch_cover_range
                                    ? 4
                                    : 0
                            )
                    ),
                    task.selected_rank,
                    task.shallow_position,
                    task.shallow_kind,
                    task.selected_complement_deep,
                    task.endpoint_pattern,
                    task.selected_local_value >= 0
                        ? task.selected_local_value
                        : switch_local_value,
                    switch_bad_threshold
                );
                if (!result.passed) {
                    {
                        std::lock_guard<std::mutex> lock(
                            diagnostic_mutex
                        );
                        if (diagnostic.empty()) {
                            diagnostic =
                                task_name(task) + " "
                                + result.diagnostic;
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

    {
        std::lock_guard<std::mutex> lock(progress_mutex);
        finished = true;
    }
    progress_condition.notify_all();
    reporter.join();

    if (failed.load() || completed.load() != tasks.size()) {
        std::cerr << tag << " FAIL"
                  << " completed=" << completed.load()
                  << " tasks=" << tasks.size()
                  << " diagnostic=" << diagnostic << '\n';
        return EXIT_FAILURE;
    }
    std::cout << tag
              << " tasks=" << tasks.size()
              << " counterexamples=UNSAT result=PASS\n";
    return EXIT_SUCCESS;
}
