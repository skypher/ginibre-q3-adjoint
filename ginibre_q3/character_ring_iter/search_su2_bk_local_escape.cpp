#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

std::size_t maximum_source_escape
    = std::numeric_limits<std::size_t>::max();
bool require_almost_square_free = false;

struct State {
    std::vector<std::size_t> order;
    std::vector<int> alpha;
};

struct Potential {
    std::int64_t delta = 0;
    std::uint64_t total = 0U;
};

struct SourceTotals {
    std::uint64_t empty = 0U;
    std::uint64_t odd = 0U;
    std::uint64_t positive = 0U;
};

struct HardSourceShrinkObstruction {
    bool present = false;
    std::uint64_t queue_path_index = 0U;
    std::size_t target = 0U;
    std::uint64_t sign_mask = 0U;
    std::vector<int> alpha;
    std::vector<int> beta;
    std::vector<int> signs;
    Potential transfer;
    Potential direct;
    SourceTotals totals;
};

struct SearchCounts {
    std::uint64_t queue_paths = 0U;
    std::uint64_t selected_queue_paths = 0U;
    std::uint64_t signed_states = 0U;
    std::uint64_t negative_states = 0U;
    std::uint64_t direct_escapes = 0U;
    std::uint64_t two_step_escapes = 0U;
    std::uint64_t three_step_braid_escapes = 0U;
    std::uint64_t longer_plateau_escapes = 0U;
    std::uint64_t source_deficit_beyond_two_steps = 0U;
    std::uint64_t source_forest_word_failures = 0U;
    std::uint64_t source_prefix_forest_word_failures = 0U;
    std::uint64_t source_nonzero_transfer_prefix_states = 0U;
    std::uint64_t source_nonzero_transfer_prefix_steps = 0U;
    std::uint64_t source_multiple_transfer_prefix_states = 0U;
    std::uint64_t source_prefix_threshold_tie_steps = 0U;
    std::uint64_t source_prefix_positive_scheduling_steps = 0U;
    std::uint64_t source_prefix_negative_scheduling_steps = 0U;
    std::uint64_t source_prefix_mixed_scheduling_states = 0U;
    std::uint64_t source_prefix_multiple_strict_scheduling_states = 0U;
    std::size_t source_maximum_strict_scheduling_steps = 0U;
    std::uint64_t source_extremal_shell_states = 0U;
    std::uint64_t source_unidirectional_shell_escapes = 0U;
    std::uint64_t source_extremal_shell_escapes = 0U;
    std::uint64_t source_extremal_sink_escapes = 0U;
    std::uint64_t source_extremal_source_escapes = 0U;
    std::uint64_t source_extremal_source_shrink_escapes = 0U;
    std::uint64_t source_extremal_source_enlarge_escapes = 0U;
    std::uint64_t source_extremal_source_positive_escapes = 0U;
    std::uint64_t source_extremal_shell_obstructions = 0U;
    std::uint64_t source_extremal_positive_plateau_states = 0U;
    std::uint64_t source_normal_zero_transfer_states = 0U;
    std::uint64_t source_normal_zero_transfer_escapes = 0U;
    std::uint64_t source_normal_zero_transfer_failures = 0U;
    std::size_t source_maximum_zero_transfer_component = 0U;
    std::size_t source_maximum_threshold_tie_quotient = 0U;
    std::uint64_t hard_extremal_shell_states = 0U;
    std::uint64_t hard_extremal_source_shrink_escapes = 0U;
    std::uint64_t hard_extremal_source_shrink_obstructions = 0U;
    std::map<std::int64_t, std::uint64_t>
        hard_source_shrink_obstruction_margin_histogram;
    std::uint64_t source_nonprimitive_beyond_two_steps = 0U;
    std::uint64_t nonprimitive_beyond_two_steps = 0U;
    std::uint64_t source_nonprimitive_negative_states = 0U;
    std::uint64_t source_return_shell_direct_escapes = 0U;
    std::uint64_t source_return_shell_two_step_escapes = 0U;
    std::uint64_t source_return_shell_failures = 0U;
    std::uint64_t source_composite_direct_escapes = 0U;
    std::uint64_t source_composite_two_step_escapes = 0U;
    std::uint64_t source_composite_escape_failures = 0U;
    std::uint64_t source_composite_two_step_return_shell = 0U;
    std::uint64_t
        source_composite_no_normal_negative_excursion = 0U;
    std::uint64_t trivial_sector_beyond_two_steps = 0U;
    std::uint64_t closed_odd_minus_plateaus = 0U;
    std::uint64_t closed_even_minus_trivial_plateaus = 0U;
    std::size_t maximum_escape_length = 0U;
    std::map<std::vector<std::size_t>, std::uint64_t>
        longer_escape_word_histogram;
    std::map<std::vector<std::size_t>, std::uint64_t>
        source_deficit_escape_word_histogram;
    std::map<std::pair<std::int64_t, std::uint64_t>, std::uint64_t>
        source_composite_two_step_potential_histogram;
    std::map<
        std::tuple<
            std::int64_t, std::uint64_t,
            std::int64_t, std::uint64_t
        >,
        std::uint64_t
    > source_composite_two_step_factor_histogram;
    std::map<
        std::tuple<std::vector<int>, std::vector<int>, std::size_t>,
        SourceTotals
    > source_totals_cache;
    HardSourceShrinkObstruction first_hard_source_shrink_obstruction;
};

bool almost_square_free(
    const std::vector<std::size_t>& word
) {
    std::map<std::size_t, std::size_t> multiplicities;
    for (const std::size_t generator : word) {
        ++multiplicities[generator];
    }
    std::size_t repeated = 0U;
    for (const auto& [generator, multiplicity] : multiplicities) {
        (void)generator;
        if (multiplicity > 2U) {
            return false;
        }
        if (multiplicity == 2U) {
            ++repeated;
        }
    }
    return repeated <= 1U
        && word.size() <= multiplicities.size() + 1U;
}

bool square_free_proper_prefix(
    const std::vector<std::size_t>& word
) {
    std::set<std::size_t> used;
    for (std::size_t index = 0U;
         index + 1U < word.size(); ++index) {
        if (!used.insert(word[index]).second) {
            return false;
        }
    }
    return true;
}

int parse_nonnegative(const char* text, const char* name) {
    const std::string input(text);
    std::size_t used = 0U;
    const long long value = std::stoll(input, &used, 10);
    if (used != input.size() || value < 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
}

Potential transfer_potential(
    const State& state,
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target
) {
    std::map<int, std::pair<std::uint64_t, std::int64_t>>
        current{{0, {1U, 1}}};
    int total_height = 0;
    for (const std::size_t vertex : state.order) {
        const int alpha = state.alpha[vertex];
        const int beta = degrees[vertex] - alpha;
        std::map<int, std::pair<std::uint64_t, std::int64_t>> next;
        const auto add = [](
            std::pair<std::uint64_t, std::int64_t>& destination,
            std::uint64_t unsigned_count,
            std::int64_t signed_count
        ) {
            if (destination.first
                > std::numeric_limits<std::uint64_t>::max()
                    - unsigned_count) {
                throw std::runtime_error("unsigned transfer overflow");
            }
            if ((signed_count > 0
                 && destination.second
                    > std::numeric_limits<std::int64_t>::max()
                        - signed_count)
                || (signed_count < 0
                    && destination.second
                        < std::numeric_limits<std::int64_t>::min()
                            - signed_count)) {
                throw std::runtime_error("signed transfer overflow");
            }
            destination.first += unsigned_count;
            destination.second += signed_count;
        };
        for (const auto& [cut_height, counts] : current) {
            if (total_height - cut_height >= beta) {
                add(next[cut_height], counts.first, counts.second);
            }
            if (vertex != target && cut_height >= beta) {
                add(
                    next[cut_height + 2 * alpha - degrees[vertex]],
                    counts.first,
                    static_cast<std::int64_t>(signs[vertex])
                        * counts.second
                );
            }
        }
        total_height += 2 * alpha - degrees[vertex];
        current = std::move(next);
    }
    if (total_height != 0) {
        throw std::runtime_error("nonzero terminal queue height");
    }
    const auto found = current.find(0);
    if (found == current.end()) {
        return {};
    }
    return Potential{found->second.second, found->second.first};
}

Potential direct_potential(
    const State& state,
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target
) {
    if (state.order.size()
        >= std::numeric_limits<std::uint64_t>::digits) {
        throw std::runtime_error("too many vertices for direct check");
    }
    const std::uint64_t subset_limit
        = std::uint64_t{1} << state.order.size();
    Potential result;
    for (std::uint64_t subset = 0U;
         subset < subset_limit; ++subset) {
        if ((subset & (std::uint64_t{1} << target)) != 0U) {
            continue;
        }
        int cut_height = 0;
        int complement_height = 0;
        int subset_sign = 1;
        bool feasible = true;
        for (const std::size_t vertex : state.order) {
            const int beta = degrees[vertex] - state.alpha[vertex];
            int& height = (subset
                & (std::uint64_t{1} << vertex)) != 0U
                ? cut_height : complement_height;
            if (height < beta) {
                feasible = false;
                break;
            }
            height += 2 * state.alpha[vertex] - degrees[vertex];
            if ((subset & (std::uint64_t{1} << vertex)) != 0U) {
                subset_sign *= signs[vertex];
            }
        }
        if (!feasible || cut_height != 0 || complement_height != 0) {
            continue;
        }
        ++result.total;
        result.delta += subset_sign;
    }
    return result;
}

std::uint64_t singlet_multiplicity(
    const std::vector<int>& labels
) {
    std::map<int, std::uint64_t> current{{0, 1U}};
    for (const int label : labels) {
        std::map<int, std::uint64_t> next;
        for (const auto& [incoming, multiplicity] : current) {
            for (int outgoing = std::abs(incoming - label);
                 outgoing <= incoming + label; outgoing += 2) {
                if (next[outgoing]
                    > std::numeric_limits<std::uint64_t>::max()
                        - multiplicity) {
                    throw std::runtime_error(
                        "singlet multiplicity overflow"
                    );
                }
                next[outgoing] += multiplicity;
            }
        }
        current = std::move(next);
    }
    const auto found = current.find(0);
    return found == current.end() ? 0U : found->second;
}

SourceTotals source_totals(
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target
) {
    std::vector<std::size_t> factors;
    for (std::size_t vertex = 0U; vertex < degrees.size(); ++vertex) {
        if (vertex != target) {
            factors.push_back(vertex);
        }
    }
    if (factors.size()
        >= std::numeric_limits<std::uint64_t>::digits) {
        throw std::runtime_error(
            "too many factors for source-total enumeration"
        );
    }
    SourceTotals totals;
    const std::uint64_t subset_limit
        = std::uint64_t{1} << factors.size();
    for (std::uint64_t subset = 0U;
         subset < subset_limit; ++subset) {
        std::vector<int> left;
        std::vector<int> right{degrees[target]};
        std::size_t minus_count = 0U;
        for (std::size_t index = 0U; index < factors.size(); ++index) {
            const std::size_t vertex = factors[index];
            if ((subset & (std::uint64_t{1} << index)) != 0U) {
                left.push_back(degrees[vertex]);
                if (signs[vertex] < 0) {
                    ++minus_count;
                }
            } else {
                right.push_back(degrees[vertex]);
            }
        }
        const std::uint64_t left_multiplicity
            = singlet_multiplicity(left);
        const std::uint64_t right_multiplicity
            = singlet_multiplicity(right);
        if (left_multiplicity != 0U
            && right_multiplicity
                > std::numeric_limits<std::uint64_t>::max()
                    / left_multiplicity) {
            throw std::runtime_error("source product overflow");
        }
        const std::uint64_t product
            = left_multiplicity * right_multiplicity;
        std::uint64_t& destination = subset == 0U
            ? totals.empty
            : ((minus_count & 1U) != 0U
                ? totals.odd : totals.positive);
        if (destination
            > std::numeric_limits<std::uint64_t>::max() - product) {
            throw std::runtime_error("source total overflow");
        }
        destination += product;
    }
    return totals;
}

const SourceTotals& cached_source_totals(
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target,
    SearchCounts& counts
) {
    std::vector<std::pair<int, int>> signed_degrees;
    signed_degrees.reserve(degrees.size() - 1U);
    for (std::size_t vertex = 0U; vertex < degrees.size(); ++vertex) {
        if (vertex != target) {
            signed_degrees.emplace_back(
                degrees[vertex], signs[vertex]
            );
        }
    }
    std::sort(signed_degrees.begin(), signed_degrees.end());
    std::vector<int> canonical_degrees;
    std::vector<int> canonical_signs;
    canonical_degrees.reserve(signed_degrees.size());
    canonical_signs.reserve(signed_degrees.size());
    for (const auto& [degree, sign] : signed_degrees) {
        canonical_degrees.push_back(degree);
        canonical_signs.push_back(sign);
    }
    const auto key = std::make_tuple(
        std::move(canonical_degrees),
        std::move(canonical_signs),
        static_cast<std::size_t>(degrees[target])
    );
    const auto found = counts.source_totals_cache.find(key);
    if (found != counts.source_totals_cache.end()) {
        return found->second;
    }
    return counts.source_totals_cache.emplace(
        key, source_totals(degrees, signs, target)
    ).first->second;
}

int bk_transfer_amount(
    const State& state,
    const std::vector<int>& degrees,
    std::size_t first
) {
    if (first + 1U >= state.order.size()) {
        throw std::runtime_error("BK generator out of range");
    }
    int height = 0;
    for (std::size_t position = 0U; position < first; ++position) {
        const std::size_t vertex = state.order[position];
        height += 2 * state.alpha[vertex] - degrees[vertex];
    }
    const std::size_t a = state.order[first];
    const std::size_t b = state.order[first + 1U];
    const int beta_a = degrees[a] - state.alpha[a];
    const int beta_b = degrees[b] - state.alpha[b];
    return std::max(0, beta_a + beta_b - height);
}

bool scheduling_precedes(
    int alpha_c,
    int beta_c,
    int alpha_d,
    int beta_d
) {
    const bool productive_c = alpha_c >= beta_c;
    const bool productive_d = alpha_d >= beta_d;
    if (productive_c != productive_d) {
        return productive_c;
    }
    if (productive_c) {
        return beta_c <= beta_d;
    }
    return alpha_c >= alpha_d;
}

int scheduling_move_direction(
    const State& state,
    const std::vector<int>& degrees,
    std::size_t first
) {
    const std::size_t left = state.order[first];
    const std::size_t right = state.order[first + 1U];
    const int alpha_left = state.alpha[left];
    const int alpha_right = state.alpha[right];
    const int beta_left = degrees[left] - alpha_left;
    const int beta_right = degrees[right] - alpha_right;
    const int left_right_threshold = beta_left + std::max(
        0, beta_right - alpha_left
    );
    const int right_left_threshold = beta_right + std::max(
        0, beta_left - alpha_right
    );
    if (left_right_threshold == right_left_threshold) {
        return 0;
    }
    const bool left_before_right = scheduling_precedes(
        alpha_left, beta_left, alpha_right, beta_right
    );
    const bool right_before_left = scheduling_precedes(
        alpha_right, beta_right, alpha_left, beta_left
    );
    if (left_before_right && right_before_left) {
        return 0;
    }
    if (right_before_left) {
        return 1;
    }
    if (left_before_right) {
        return -1;
    }
    throw std::runtime_error("scheduling preorder is not total");
}

State bk_neighbor(
    const State& state,
    const std::vector<int>& degrees,
    std::size_t first
) {
    const int z = bk_transfer_amount(state, degrees, first);
    const std::size_t a = state.order[first];
    const std::size_t b = state.order[first + 1U];
    State neighbor = state;
    neighbor.alpha[a] -= z;
    neighbor.alpha[b] += z;
    std::swap(neighbor.order[first], neighbor.order[first + 1U]);
    return neighbor;
}

void record_source_prefix_transfers(
    State state,
    const std::vector<int>& degrees,
    const std::vector<std::size_t>& word,
    SearchCounts& counts
) {
    std::size_t nonzero_transfers = 0U;
    bool has_positive_direction = false;
    bool has_negative_direction = false;
    std::size_t strict_scheduling_steps = 0U;
    for (std::size_t index = 0U;
         index + 1U < word.size(); ++index) {
        const std::size_t generator = word[index];
        const int transfer = bk_transfer_amount(
            state, degrees, generator
        );
        if (transfer > 0) {
            ++nonzero_transfers;
            ++counts.source_nonzero_transfer_prefix_steps;
        } else {
            const int direction = scheduling_move_direction(
                state, degrees, generator
            );
            if (direction > 0) {
                has_positive_direction = true;
                ++strict_scheduling_steps;
                ++counts.source_prefix_positive_scheduling_steps;
            } else if (direction < 0) {
                has_negative_direction = true;
                ++strict_scheduling_steps;
                ++counts.source_prefix_negative_scheduling_steps;
            } else {
                ++counts.source_prefix_threshold_tie_steps;
            }
        }
        state = bk_neighbor(state, degrees, generator);
    }
    if (nonzero_transfers != 0U) {
        ++counts.source_nonzero_transfer_prefix_states;
    }
    if (nonzero_transfers > 1U) {
        ++counts.source_multiple_transfer_prefix_states;
    }
    if (has_positive_direction && has_negative_direction) {
        ++counts.source_prefix_mixed_scheduling_states;
    }
    if (strict_scheduling_steps > 1U) {
        ++counts.source_prefix_multiple_strict_scheduling_states;
    }
    counts.source_maximum_strict_scheduling_steps = std::max(
        counts.source_maximum_strict_scheduling_steps,
        strict_scheduling_steps
    );
}

bool primitive_queue_path(
    const State& state,
    const std::vector<int>& degrees
) {
    int height = 0;
    for (std::size_t position = 0U;
         position < state.order.size(); ++position) {
        const std::size_t vertex = state.order[position];
        height += 2 * state.alpha[vertex] - degrees[vertex];
        if (position + 1U < state.order.size() && height == 0) {
            return false;
        }
    }
    return height == 0;
}

std::vector<std::size_t> return_boundaries(
    const State& state,
    const std::vector<int>& degrees
) {
    std::vector<std::size_t> result;
    int height = 0;
    for (std::size_t position = 0U;
         position + 1U < state.order.size(); ++position) {
        const std::size_t vertex = state.order[position];
        height += 2 * state.alpha[vertex] - degrees[vertex];
        if (height == 0) {
            result.push_back(position);
        }
    }
    return result;
}

std::tuple<
    std::int64_t, std::uint64_t,
    std::int64_t, std::uint64_t
> two_excursion_factor_signature(
    const State& state,
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target
) {
    const std::vector<std::size_t> boundaries
        = return_boundaries(state, degrees);
    if (boundaries.size() != 1U) {
        throw std::runtime_error(
            "two-excursion signature needs exactly one return"
        );
    }
    const std::size_t split = boundaries.front() + 1U;
    State left;
    State right;
    left.alpha = state.alpha;
    right.alpha = state.alpha;
    left.order.assign(
        state.order.begin(),
        state.order.begin() + static_cast<std::ptrdiff_t>(split)
    );
    right.order.assign(
        state.order.begin() + static_cast<std::ptrdiff_t>(split),
        state.order.end()
    );
    const bool target_on_left = std::find(
        left.order.begin(), left.order.end(), target
    ) != left.order.end();
    const State& target_excursion = target_on_left ? left : right;
    const State& other_excursion = target_on_left ? right : left;
    const Potential target_factor = transfer_potential(
        target_excursion, degrees, signs, target
    );
    const Potential other_target_free = transfer_potential(
        other_excursion, degrees, signs, degrees.size()
    );
    if ((other_target_free.delta & 1) != 0
        || (other_target_free.total & 1U) != 0U) {
        throw std::runtime_error(
            "target-free excursion factor is not even"
        );
    }
    return {
        target_factor.delta,
        target_factor.total,
        other_target_free.delta / 2,
        other_target_free.total / 2U
    };
}

bool lower(const Potential& candidate, const Potential& original);
bool plateau(const Potential& candidate, const Potential& original);
std::vector<std::size_t> plateau_escape_word(
    const State& initial,
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target,
    const Potential& original
);

bool has_normal_escaping_negative_excursion(
    const State& state,
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target
) {
    std::vector<std::size_t> endpoints
        = return_boundaries(state, degrees);
    endpoints.push_back(state.order.size() - 1U);
    std::size_t begin = 0U;
    bool found_negative = false;
    for (const std::size_t endpoint : endpoints) {
        State excursion;
        excursion.alpha = state.alpha;
        excursion.order.assign(
            state.order.begin()
                + static_cast<std::ptrdiff_t>(begin),
            state.order.begin()
                + static_cast<std::ptrdiff_t>(endpoint + 1U)
        );
        const auto found_target = std::find(
            excursion.order.begin(), excursion.order.end(), target
        );
        const std::size_t local_target
            = found_target != excursion.order.end()
                ? target : excursion.order.front();
        const Potential local = transfer_potential(
            excursion, degrees, signs, local_target
        );
        if (local.delta < 0) {
            found_negative = true;
            const std::vector<std::size_t> word
                = plateau_escape_word(
                    excursion, degrees, signs, local_target, local
                );
            State prefix_state = excursion;
            std::size_t positive_transfers = 0U;
            for (std::size_t index = 0U;
                 index + 1U < word.size(); ++index) {
                if (bk_transfer_amount(
                        prefix_state, degrees, word[index]
                    ) > 0) {
                    ++positive_transfers;
                }
                prefix_state = bk_neighbor(
                    prefix_state, degrees, word[index]
                );
            }
            if (!word.empty()
                && word.size() <= excursion.order.size()
                && square_free_proper_prefix(word)
                && positive_transfers <= 1U) {
                return true;
            }
        }
        begin = endpoint + 1U;
    }
    if (!found_negative) {
        throw std::runtime_error(
            "negative global state has no negative excursion"
        );
    }
    return false;
}

int return_shell_escape_depth(
    const State& state,
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target,
    const Potential& original
) {
    for (const std::size_t boundary :
         return_boundaries(state, degrees)) {
        const State neighbor = bk_neighbor(
            state, degrees, boundary
        );
        if (lower(
                transfer_potential(
                    neighbor, degrees, signs, target
                ),
                original
            )) {
            return 1;
        }
    }
    for (std::size_t first = 0U;
         first + 1U < state.order.size(); ++first) {
        const State middle = bk_neighbor(
            state, degrees, first
        );
        if (!plateau(
                transfer_potential(
                    middle, degrees, signs, target
                ),
                original
            )) {
            continue;
        }
        for (const std::size_t boundary :
             return_boundaries(middle, degrees)) {
            const State neighbor = bk_neighbor(
                middle, degrees, boundary
            );
            if (lower(
                    transfer_potential(
                        neighbor, degrees, signs, target
                    ),
                    original
                )) {
                return 2;
            }
        }
    }
    return 0;
}

int generic_two_step_escape_depth(
    const State& state,
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target,
    const Potential& original
) {
    std::vector<State> neighbors;
    std::vector<Potential> values;
    for (std::size_t first = 0U;
         first + 1U < state.order.size(); ++first) {
        neighbors.push_back(bk_neighbor(
            state, degrees, first
        ));
        values.push_back(transfer_potential(
            neighbors.back(), degrees, signs, target
        ));
        if (lower(values.back(), original)) {
            return 1;
        }
    }
    for (std::size_t first = 0U;
         first + 1U < state.order.size(); ++first) {
        if (!plateau(values[first], original)) {
            continue;
        }
        for (const int direction : {-1, 1}) {
            const int second_signed
                = static_cast<int>(first) + direction;
            if (second_signed < 0
                || second_signed + 1
                    >= static_cast<int>(state.order.size())) {
                continue;
            }
            const std::size_t second
                = static_cast<std::size_t>(second_signed);
            const State neighbor = bk_neighbor(
                neighbors[first], degrees, second
            );
            if (lower(
                    transfer_potential(
                        neighbor, degrees, signs, target
                    ),
                    original
                )) {
                return 2;
            }
        }
    }
    return 0;
}

bool lower(const Potential& candidate, const Potential& original) {
    return candidate.delta > original.delta
        || (candidate.delta == original.delta
            && candidate.total < original.total);
}

bool plateau(const Potential& candidate, const Potential& original) {
    return candidate.delta == original.delta
        && candidate.total == original.total;
}

struct ExtremalShellAnalysis {
    std::size_t states = 0U;
    std::size_t tie_components = 0U;
    bool has_unidirectional_escape = false;
    bool has_extremal_escape = false;
    bool has_sink_escape = false;
    bool has_source_escape = false;
    bool has_source_shrink_escape = false;
    bool has_source_enlarge_escape = false;
    bool has_source_positive_escape = false;
    bool has_reachable_extremal_positive_plateau = false;
};

class DisjointSet {
public:
    std::size_t add() {
        const std::size_t index = parent_.size();
        parent_.push_back(index);
        rank_.push_back(0U);
        return index;
    }

    std::size_t find(std::size_t value) {
        if (parent_[value] != value) {
            parent_[value] = find(parent_[value]);
        }
        return parent_[value];
    }

    void unite(std::size_t first, std::size_t second) {
        first = find(first);
        second = find(second);
        if (first == second) {
            return;
        }
        if (rank_[first] < rank_[second]) {
            std::swap(first, second);
        }
        parent_[second] = first;
        if (rank_[first] == rank_[second]) {
            ++rank_[first];
        }
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<std::size_t> rank_;
};

ExtremalShellAnalysis analyze_zero_transfer_extrema(
    const State& initial,
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target,
    const Potential& original
) {
    using Key = std::pair<std::vector<std::size_t>, std::vector<int>>;
    struct Edge {
        std::size_t first = 0U;
        std::size_t second = 0U;
        int direction = 0;
    };

    std::vector<State> states{initial};
    std::map<Key, std::size_t> indices;
    indices.emplace(
        Key{initial.order, initial.alpha}, std::size_t{0}
    );
    std::queue<std::size_t> frontier;
    frontier.push(0U);
    std::vector<Edge> edges;
    DisjointSet tie_sets;
    tie_sets.add();

    while (!frontier.empty()) {
        const std::size_t current_index = frontier.front();
        frontier.pop();
        const State current = states[current_index];
        for (std::size_t generator = 0U;
             generator + 1U < current.order.size(); ++generator) {
            if (bk_transfer_amount(current, degrees, generator) != 0) {
                continue;
            }
            State neighbor = bk_neighbor(
                current, degrees, generator
            );
            if (!plateau(
                    transfer_potential(
                        neighbor, degrees, signs, target
                    ),
                    original
                )) {
                continue;
            }
            const Key key{neighbor.order, neighbor.alpha};
            const auto [found, inserted] = indices.emplace(
                key, states.size()
            );
            if (inserted) {
                states.push_back(std::move(neighbor));
                tie_sets.add();
                frontier.push(found->second);
            }
            if (current_index < found->second) {
                const int direction = scheduling_move_direction(
                    current, degrees, generator
                );
                edges.push_back(
                    Edge{current_index, found->second, direction}
                );
                if (direction == 0) {
                    tie_sets.unite(current_index, found->second);
                }
            }
        }
    }

    std::map<std::size_t, std::size_t> quotient_indices;
    std::vector<std::size_t> state_quotient(states.size(), 0U);
    for (std::size_t index = 0U; index < states.size(); ++index) {
        const std::size_t root = tie_sets.find(index);
        const auto [found, inserted] = quotient_indices.emplace(
            root, quotient_indices.size()
        );
        (void)inserted;
        state_quotient[index] = found->second;
    }
    const std::size_t quotient_size = quotient_indices.size();
    std::vector<std::set<std::size_t>> outgoing(quotient_size);
    std::vector<std::set<std::size_t>> incoming(quotient_size);
    for (const Edge& edge : edges) {
        if (edge.direction == 0) {
            continue;
        }
        std::size_t source = state_quotient[edge.first];
        std::size_t destination = state_quotient[edge.second];
        if (edge.direction < 0) {
            std::swap(source, destination);
        }
        if (source == destination) {
            throw std::runtime_error(
                "strict scheduling edge collapsed by threshold ties"
            );
        }
        outgoing[source].insert(destination);
        incoming[destination].insert(source);
    }

    const std::size_t initial_quotient = state_quotient[0U];
    const auto reachable = [](
        const std::vector<std::set<std::size_t>>& adjacency,
        std::size_t start
    ) {
        std::vector<bool> seen(adjacency.size(), false);
        std::queue<std::size_t> queue;
        seen[start] = true;
        queue.push(start);
        while (!queue.empty()) {
            const std::size_t current = queue.front();
            queue.pop();
            for (const std::size_t next : adjacency[current]) {
                if (!seen[next]) {
                    seen[next] = true;
                    queue.push(next);
                }
            }
        }
        return seen;
    };
    const std::vector<bool> forward = reachable(
        outgoing, initial_quotient
    );
    const std::vector<bool> backward = reachable(
        incoming, initial_quotient
    );
    std::vector<bool> sink(quotient_size, false);
    std::vector<bool> source(quotient_size, false);
    std::vector<bool> extremal(quotient_size, false);
    for (std::size_t component = 0U;
         component < quotient_size; ++component) {
        sink[component]
            = forward[component] && outgoing[component].empty();
        source[component]
            = backward[component] && incoming[component].empty();
        extremal[component] = sink[component] || source[component];
    }

    ExtremalShellAnalysis result;
    result.states = states.size();
    result.tie_components = quotient_size;
    for (std::size_t index = 0U; index < states.size(); ++index) {
        const std::size_t component = state_quotient[index];
        const bool unidirectional
            = forward[component] || backward[component];
        if (!unidirectional) {
            continue;
        }
        const State& state = states[index];
        for (std::size_t generator = 0U;
             generator + 1U < state.order.size(); ++generator) {
            const int transfer = bk_transfer_amount(
                state, degrees, generator
            );
            State neighbor = bk_neighbor(
                state, degrees, generator
            );
            const Potential candidate = transfer_potential(
                neighbor, degrees, signs, target
            );
            if (lower(candidate, original)) {
                result.has_unidirectional_escape = true;
                if (extremal[component]) {
                    result.has_extremal_escape = true;
                }
                if (sink[component]) {
                    result.has_sink_escape = true;
                }
                if (source[component]) {
                    result.has_source_escape = true;
                    if (transfer > 0) {
                        result.has_source_positive_escape = true;
                    } else {
                        const int direction
                            = scheduling_move_direction(
                                state, degrees, generator
                            );
                        if (direction < 0) {
                            result.has_source_shrink_escape = true;
                        } else if (direction > 0) {
                            result.has_source_enlarge_escape = true;
                        } else {
                            throw std::runtime_error(
                                "threshold-tie edge is not a plateau"
                            );
                        }
                    }
                }
            }
            if (extremal[component]
                && transfer > 0 && plateau(candidate, original)) {
                result.has_reachable_extremal_positive_plateau = true;
            }
        }
    }
    return result;
}

bool has_normal_zero_transfer_escape(
    const State& initial,
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target,
    const Potential& original
) {
    if (initial.order.size()
        > std::numeric_limits<std::uint64_t>::digits) {
        throw std::runtime_error(
            "too many generators for square-free mask"
        );
    }
    struct Node {
        State state;
        std::uint64_t used = 0U;
        std::size_t strict = 0U;
    };
    using Key = std::tuple<
        std::vector<std::size_t>,
        std::vector<int>,
        std::uint64_t,
        std::size_t
    >;
    for (const int required_direction : {-1, 1}) {
        std::queue<Node> frontier;
        std::set<Key> visited;
        frontier.push(Node{initial, 0U, 0U});
        visited.emplace(
            initial.order, initial.alpha, 0U, std::size_t{0}
        );
        while (!frontier.empty()) {
            Node current = std::move(frontier.front());
            frontier.pop();
            for (std::size_t generator = 0U;
                 generator + 1U < current.state.order.size();
                 ++generator) {
                State neighbor = bk_neighbor(
                    current.state, degrees, generator
                );
                const Potential candidate = transfer_potential(
                    neighbor, degrees, signs, target
                );
                if (lower(candidate, original)) {
                    return true;
                }
                const std::uint64_t bit
                    = std::uint64_t{1} << generator;
                if (!plateau(candidate, original)
                    || bk_transfer_amount(
                        current.state, degrees, generator
                    ) != 0
                    || (current.used & bit) != 0U) {
                    continue;
                }
                const int direction = scheduling_move_direction(
                    current.state, degrees, generator
                );
                if (direction != 0
                    && direction != required_direction) {
                    continue;
                }
                const std::size_t strict
                    = current.strict
                      + (direction == 0 ? 0U : 1U);
                if (strict > 2U) {
                    continue;
                }
                const std::uint64_t used = current.used | bit;
                if (!visited.emplace(
                        neighbor.order, neighbor.alpha, used, strict
                    ).second) {
                    continue;
                }
                frontier.push(Node{
                    std::move(neighbor), used, strict
                });
            }
        }
    }
    return false;
}

std::vector<std::size_t> plateau_escape_word(
    const State& initial,
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target,
    const Potential& original
) {
    using Key = std::pair<std::vector<std::size_t>, std::vector<int>>;
    struct SearchNode {
        State state;
        std::vector<std::size_t> word;
    };
    std::queue<SearchNode> frontier;
    std::set<Key> visited;
    frontier.push(SearchNode{initial, {}});
    visited.emplace(initial.order, initial.alpha);
    while (!frontier.empty()) {
        SearchNode current = std::move(frontier.front());
        frontier.pop();
        for (std::size_t generator = 0U;
             generator + 1U < current.state.order.size();
             ++generator) {
            State neighbor = bk_neighbor(
                current.state, degrees, generator
            );
            const Potential value = transfer_potential(
                neighbor, degrees, signs, target
            );
            std::vector<std::size_t> word = current.word;
            word.push_back(generator);
            if (lower(value, original)) {
                return word;
            }
            if (!plateau(value, original)
                || !visited.emplace(
                    neighbor.order, neighbor.alpha
                ).second) {
                continue;
            }
            frontier.push(SearchNode{
                std::move(neighbor), std::move(word)
            });
        }
    }
    return {};
}

void print_vector(const std::vector<int>& values) {
    std::cout << '[';
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << ']';
}

[[noreturn]] void report_failure(
    const State& state,
    const std::vector<int>& degrees,
    const std::vector<int>& signs,
    std::size_t target,
    const Potential& original
) {
    std::vector<int> ordered_alpha;
    std::vector<int> ordered_beta;
    std::vector<int> ordered_signs;
    ordered_alpha.reserve(state.order.size());
    ordered_beta.reserve(state.order.size());
    ordered_signs.reserve(state.order.size());
    for (const std::size_t vertex : state.order) {
        ordered_alpha.push_back(state.alpha[vertex]);
        ordered_beta.push_back(degrees[vertex] - state.alpha[vertex]);
        ordered_signs.push_back(
            vertex == target ? 0 : signs[vertex]
        );
    }
    std::cout << "SU2_BK_LOCAL_ESCAPE result=FAIL target_position=";
    const auto found = std::find(
        state.order.begin(), state.order.end(), target
    );
    std::cout << std::distance(state.order.begin(), found)
              << " delta=" << original.delta
              << " total=" << original.total << " alpha=";
    print_vector(ordered_alpha);
    std::cout << " beta=";
    print_vector(ordered_beta);
    std::cout << " signs=";
    print_vector(ordered_signs);
    const Potential direct = direct_potential(
        state, degrees, signs, target
    );
    std::cout << " direct_delta=" << direct.delta
              << " direct_total=" << direct.total;
    if (direct.delta != original.delta
        || direct.total != original.total) {
        throw std::runtime_error(
            "direct and transfer potentials disagree"
        );
    }
    const std::vector<std::size_t> escape_word
        = plateau_escape_word(
            state, degrees, signs, target, original
        );
    std::cout << " plateau_escape_word=[";
    for (std::size_t index = 0U; index < escape_word.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << escape_word[index];
    }
    std::cout << ']';
    const SourceTotals totals = source_totals(
        degrees, signs, target
    );
    std::cout << " empty_sources=" << totals.empty
              << " odd_sources=" << totals.odd
              << " positive_sources=" << totals.positive;
    if (totals.empty > std::numeric_limits<std::uint64_t>::max()
            - totals.positive) {
        throw std::runtime_error("signed source total overflow");
    }
    const std::uint64_t nonnegative_sources
        = totals.empty + totals.positive;
    if (nonnegative_sources >= totals.odd) {
        std::cout << " signed_source_integer="
                  << nonnegative_sources - totals.odd;
    } else {
        std::cout << " signed_source_integer=-"
                  << totals.odd - nonnegative_sources;
    }
    std::cout << " neighbors={";
    for (std::size_t first = 0U;
         first + 1U < state.order.size(); ++first) {
        if (first != 0U) {
            std::cout << ',';
        }
        const State neighbor = bk_neighbor(state, degrees, first);
        const Potential value = transfer_potential(
            neighbor, degrees, signs, target
        );
        std::cout << first << ":(" << value.delta
                  << ',' << value.total << ')';
    }
    std::cout << "}\n";
    std::exit(EXIT_FAILURE);
}

void test_signed_state(
    const std::vector<int>& alpha,
    const std::vector<int>& beta,
    std::size_t target,
    const std::vector<int>& signs,
    std::uint64_t queue_path_index,
    std::uint64_t sign_mask,
    SearchCounts& counts,
    bool source_filtered
) {
    const std::size_t vertices = alpha.size();
    std::vector<int> degrees(vertices, 0);
    for (std::size_t index = 0U; index < vertices; ++index) {
        degrees[index] = alpha[index] + beta[index];
    }
    State state;
    state.alpha = alpha;
    state.order.resize(vertices);
    for (std::size_t index = 0U; index < vertices; ++index) {
        state.order[index] = index;
    }
    ++counts.signed_states;
    const Potential original = transfer_potential(
        state, degrees, signs, target
    );
    if (original.delta >= 0) {
        return;
    }
    ++counts.negative_states;
    if (source_filtered
        && !primitive_queue_path(state, degrees)) {
        ++counts.source_nonprimitive_negative_states;
        const int depth = return_shell_escape_depth(
            state, degrees, signs, target, original
        );
        if (depth == 1) {
            ++counts.source_return_shell_direct_escapes;
        } else if (depth == 2) {
            ++counts.source_return_shell_two_step_escapes;
        } else {
            ++counts.source_return_shell_failures;
        }
        const int generic_depth = generic_two_step_escape_depth(
            state, degrees, signs, target, original
        );
        if (generic_depth == 1) {
            ++counts.source_composite_direct_escapes;
        } else if (generic_depth == 2) {
            ++counts.source_composite_two_step_escapes;
            ++counts.source_composite_two_step_potential_histogram[
                {original.delta, original.total}
            ];
            if (original.delta == -2) {
                ++counts.source_composite_two_step_factor_histogram[
                    two_excursion_factor_signature(
                        state, degrees, signs, target
                    )
                ];
            }
            if (depth == 2) {
                ++counts.source_composite_two_step_return_shell;
            }
        } else {
            ++counts.source_composite_escape_failures;
        }
        if (!has_normal_escaping_negative_excursion(
                state, degrees, signs, target
            )) {
            ++counts
                .source_composite_no_normal_negative_excursion;
        }
    }
    std::size_t minus_factors = 0U;
    for (std::size_t vertex = 0U; vertex < vertices; ++vertex) {
        if (vertex != target && signs[vertex] < 0) {
            ++minus_factors;
        }
    }
    const Potential target_free = transfer_potential(
        state, degrees, signs, vertices
    );
    if (original.total
            > std::numeric_limits<std::uint64_t>::max() / 2U
        || original.delta
            > std::numeric_limits<std::int64_t>::max() / 2
        || original.delta
            < std::numeric_limits<std::int64_t>::min() / 2) {
        throw std::runtime_error(
            "target-free complement comparison overflow"
        );
    }
    if (target_free.total != 2U * original.total
        || ((minus_factors & 1U) == 0U
            && target_free.delta != 2 * original.delta)
        || ((minus_factors & 1U) != 0U
            && target_free.delta != 0)) {
        throw std::runtime_error(
            "target-free complement identity failed"
        );
    }
    std::vector<State> neighbors;
    std::vector<Potential> values;
    neighbors.reserve(vertices - 1U);
    values.reserve(vertices - 1U);
    for (std::size_t first = 0U; first + 1U < vertices; ++first) {
        neighbors.push_back(bk_neighbor(state, degrees, first));
        values.push_back(transfer_potential(
            neighbors.back(), degrees, signs, target
        ));
        if (lower(values.back(), original)) {
            ++counts.direct_escapes;
            counts.maximum_escape_length = std::max(
                counts.maximum_escape_length, std::size_t{1}
            );
            return;
        }
    }
    bool normal_zero_transfer_checked = false;
    const auto check_normal_zero_transfer = [&]() {
        if (normal_zero_transfer_checked) {
            return;
        }
        normal_zero_transfer_checked = true;
        ++counts.source_normal_zero_transfer_states;
        if (has_normal_zero_transfer_escape(
                state, degrees, signs, target, original
            )) {
            ++counts.source_normal_zero_transfer_escapes;
        } else {
            ++counts.source_normal_zero_transfer_failures;
            report_failure(
                state, degrees, signs, target, original
            );
        }
    };
    if (source_filtered) {
        check_normal_zero_transfer();
    }
    for (std::size_t first = 0U; first + 1U < vertices; ++first) {
        if (!plateau(values[first], original)) {
            continue;
        }
        for (const int direction : {-1, 1}) {
            const int second_signed
                = static_cast<int>(first) + direction;
            if (second_signed < 0
                || second_signed + 1 >= static_cast<int>(vertices)) {
                continue;
            }
            const std::size_t second
                = static_cast<std::size_t>(second_signed);
            const State second_neighbor = bk_neighbor(
                neighbors[first], degrees, second
            );
            const Potential second_value = transfer_potential(
                second_neighbor, degrees, signs, target
            );
            if (lower(second_value, original)) {
                ++counts.two_step_escapes;
                counts.maximum_escape_length = std::max(
                    counts.maximum_escape_length, std::size_t{2}
                );
                return;
            }
        }
    }
    const SourceTotals& totals = cached_source_totals(
        degrees, signs, target, counts
    );
    if (!primitive_queue_path(state, degrees)) {
        ++counts.nonprimitive_beyond_two_steps;
    }
    const bool source_deficit = totals.odd > totals.positive;
    const ExtremalShellAnalysis extremal
        = analyze_zero_transfer_extrema(
            state, degrees, signs, target, original
        );
    ++counts.hard_extremal_shell_states;
    if (extremal.has_source_shrink_escape) {
        ++counts.hard_extremal_source_shrink_escapes;
    } else {
        ++counts.hard_extremal_source_shrink_obstructions;
        if (!counts.first_hard_source_shrink_obstruction.present) {
            const Potential direct = direct_potential(
                state, degrees, signs, target
            );
            if (direct.delta != original.delta
                || direct.total != original.total) {
                throw std::runtime_error(
                    "hard-obstruction direct and transfer potentials disagree"
                );
            }
            HardSourceShrinkObstruction& obstruction
                = counts.first_hard_source_shrink_obstruction;
            obstruction.present = true;
            obstruction.queue_path_index = queue_path_index;
            obstruction.target = target;
            obstruction.sign_mask = sign_mask;
            obstruction.alpha = alpha;
            obstruction.beta = beta;
            obstruction.signs = signs;
            obstruction.signs[target] = 0;
            obstruction.transfer = original;
            obstruction.direct = direct;
            obstruction.totals = totals;
        }
        if (totals.odd
                > static_cast<std::uint64_t>(
                    std::numeric_limits<std::int64_t>::max())
            || totals.positive
                > static_cast<std::uint64_t>(
                    std::numeric_limits<std::int64_t>::max())) {
            throw std::runtime_error(
                "source-shell obstruction margin overflow"
            );
        }
        const std::int64_t margin
            = static_cast<std::int64_t>(totals.positive)
              - static_cast<std::int64_t>(totals.odd);
        ++counts.hard_source_shrink_obstruction_margin_histogram[
            margin
        ];
    }
    if (source_deficit) {
        ++counts.source_deficit_beyond_two_steps;
        if (!primitive_queue_path(state, degrees)) {
            ++counts.source_nonprimitive_beyond_two_steps;
        }
        ++counts.source_extremal_shell_states;
        if (extremal.has_unidirectional_escape) {
            ++counts.source_unidirectional_shell_escapes;
        }
        if (extremal.has_extremal_escape) {
            ++counts.source_extremal_shell_escapes;
        } else {
            ++counts.source_extremal_shell_obstructions;
        }
        if (extremal.has_sink_escape) {
            ++counts.source_extremal_sink_escapes;
        }
        if (extremal.has_source_escape) {
            ++counts.source_extremal_source_escapes;
        }
        if (extremal.has_source_shrink_escape) {
            ++counts.source_extremal_source_shrink_escapes;
        }
        if (extremal.has_source_enlarge_escape) {
            ++counts.source_extremal_source_enlarge_escapes;
        }
        if (extremal.has_source_positive_escape) {
            ++counts.source_extremal_source_positive_escapes;
        }
        if (extremal.has_reachable_extremal_positive_plateau) {
            ++counts.source_extremal_positive_plateau_states;
        }
        check_normal_zero_transfer();
        counts.source_maximum_zero_transfer_component = std::max(
            counts.source_maximum_zero_transfer_component,
            extremal.states
        );
        counts.source_maximum_threshold_tie_quotient = std::max(
            counts.source_maximum_threshold_tie_quotient,
            extremal.tie_components
        );
    } else {
        ++counts.trivial_sector_beyond_two_steps;
    }
    for (std::size_t first = 0U; first + 1U < vertices; ++first) {
        if (!plateau(values[first], original)) {
            continue;
        }
        for (const int direction : {-1, 1}) {
            const int second_signed
                = static_cast<int>(first) + direction;
            if (second_signed < 0
                || second_signed + 1 >= static_cast<int>(vertices)) {
                continue;
            }
            const std::size_t second
                = static_cast<std::size_t>(second_signed);
            const State second_neighbor = bk_neighbor(
                neighbors[first], degrees, second
            );
            const Potential second_value = transfer_potential(
                second_neighbor, degrees, signs, target
            );
            if (!plateau(second_value, original)) {
                continue;
            }
            const State third_neighbor = bk_neighbor(
                second_neighbor, degrees, first
            );
            const Potential third_value = transfer_potential(
                third_neighbor, degrees, signs, target
            );
            if (lower(third_value, original)) {
                ++counts.three_step_braid_escapes;
                if (source_deficit) {
                    const std::vector<std::size_t> word{
                        first, second, first
                    };
                    ++counts.source_deficit_escape_word_histogram[word];
                    if (!almost_square_free(word)) {
                        ++counts.source_forest_word_failures;
                    }
                    if (!square_free_proper_prefix(word)) {
                        ++counts.source_prefix_forest_word_failures;
                    }
                    record_source_prefix_transfers(
                        state, degrees, word, counts
                    );
                }
                counts.maximum_escape_length = std::max(
                    counts.maximum_escape_length, std::size_t{3}
                );
                return;
            }
        }
    }
    const std::vector<std::size_t> escape_word
        = plateau_escape_word(
            state, degrees, signs, target, original
        );
    if (!escape_word.empty()) {
        if (source_deficit
            && (escape_word.size() > vertices
                || escape_word.size() > maximum_source_escape)) {
            report_failure(state, degrees, signs, target, original);
        }
        ++counts.longer_plateau_escapes;
        ++counts.longer_escape_word_histogram[escape_word];
        if (source_deficit) {
            ++counts.source_deficit_escape_word_histogram[
                escape_word
            ];
            if (!almost_square_free(escape_word)) {
                ++counts.source_forest_word_failures;
                if (require_almost_square_free) {
                    report_failure(
                        state, degrees, signs, target, original
                    );
                }
            }
            if (!square_free_proper_prefix(escape_word)) {
                ++counts.source_prefix_forest_word_failures;
                if (require_almost_square_free) {
                    report_failure(
                        state, degrees, signs, target, original
                    );
                }
            }
            record_source_prefix_transfers(
                state, degrees, escape_word, counts
            );
        }
        counts.maximum_escape_length = std::max(
            counts.maximum_escape_length, escape_word.size()
        );
        return;
    }
    if (source_deficit) {
        report_failure(state, degrees, signs, target, original);
    }
    if ((minus_factors & 1U) != 0U) {
        ++counts.closed_odd_minus_plateaus;
        return;
    }
    ++counts.closed_even_minus_trivial_plateaus;
}

void test_queue_path(
    const std::vector<int>& alpha,
    const std::vector<int>& beta,
    SearchCounts& counts,
    bool source_only,
    std::uint64_t shard_count,
    std::uint64_t shard_index
) {
    const std::uint64_t path_index = counts.queue_paths;
    ++counts.queue_paths;
    if (path_index % shard_count != shard_index) {
        return;
    }
    ++counts.selected_queue_paths;
    const std::size_t vertices = alpha.size();
    std::vector<int> degrees(vertices, 0);
    for (std::size_t index = 0U; index < vertices; ++index) {
        degrees[index] = alpha[index] + beta[index];
    }
    for (std::size_t target = 0U; target < vertices; ++target) {
        bool valid = true;
        std::vector<int> distinct_degrees;
        for (std::size_t index = 0U; index < vertices; ++index) {
            if (index == target) {
                continue;
            }
            if (degrees[index] == 0) {
                valid = false;
                break;
            }
            if (std::find(
                    distinct_degrees.begin(), distinct_degrees.end(),
                    degrees[index]
                ) == distinct_degrees.end()) {
                distinct_degrees.push_back(degrees[index]);
            }
        }
        if (!valid
            || distinct_degrees.size()
                >= std::numeric_limits<std::uint64_t>::digits) {
            continue;
        }
        const std::uint64_t sign_limit
            = std::uint64_t{1} << distinct_degrees.size();
        for (std::uint64_t sign_mask = 1U;
             sign_mask < sign_limit; ++sign_mask) {
            std::vector<int> signs(vertices, 1);
            for (std::size_t index = 0U; index < vertices; ++index) {
                if (index == target) {
                    continue;
                }
                const auto found = std::find(
                    distinct_degrees.begin(), distinct_degrees.end(),
                    degrees[index]
                );
                const std::size_t sign_index = static_cast<std::size_t>(
                    std::distance(distinct_degrees.begin(), found)
                );
                if ((sign_mask & (std::uint64_t{1} << sign_index))
                    != 0U) {
                    signs[index] = -1;
                }
            }
            if (source_only) {
                std::size_t minus_factors = 0U;
                for (std::size_t vertex = 0U;
                     vertex < vertices; ++vertex) {
                    if (vertex != target && signs[vertex] < 0) {
                        ++minus_factors;
                    }
                }
                if ((minus_factors & 1U) != 0U) {
                    continue;
                }
                const SourceTotals& totals = cached_source_totals(
                    degrees, signs, target, counts
                );
                if (totals.odd <= totals.positive) {
                    continue;
                }
            }
            test_signed_state(
                alpha, beta, target, signs, path_index, sign_mask,
                counts, source_only
            );
        }
    }
}

void enumerate_queue_paths(
    std::size_t position,
    int height,
    int part_bound,
    std::vector<int>& alpha,
    std::vector<int>& beta,
    SearchCounts& counts,
    bool source_only,
    std::uint64_t shard_count,
    std::uint64_t shard_index
) {
    if (position == alpha.size()) {
        if (height == 0) {
            test_queue_path(
                alpha, beta, counts, source_only,
                shard_count, shard_index
            );
        }
        return;
    }
    for (int bottom = 0; bottom <= std::min(part_bound, height);
         ++bottom) {
        for (int top = 0; top <= part_bound; ++top) {
            const int next_height = height + top - bottom;
            const std::size_t remaining
                = alpha.size() - position - 1U;
            if (next_height
                > static_cast<int>(remaining)
                    * part_bound) {
                continue;
            }
            alpha[position] = top;
            beta[position] = bottom;
            enumerate_queue_paths(
                position + 1U, next_height, part_bound,
                alpha, beta, counts, source_only,
                shard_count, shard_index
            );
        }
    }
}

void test_long_escape_family(int middle_jobs, int sign_mask) {
    if (middle_jobs < 2 || (middle_jobs & 1) != 0) {
        throw std::runtime_error(
            "family parameter must be a positive even integer"
        );
    }
    const int sink_jobs = (middle_jobs + 4) / 2;
    const std::size_t vertices = static_cast<std::size_t>(
        middle_jobs + sink_jobs + 3
    );
    std::vector<int> alpha(vertices, 0);
    std::vector<int> beta(vertices, 0);
    if (sign_mask < 0 || sign_mask >= 8) {
        throw std::runtime_error("family sign mask must be below eight");
    }
    std::vector<int> signs(vertices, 1);
    alpha[0U] = 1;
    alpha[1U] = 2;
    for (int index = 0; index < middle_jobs; ++index) {
        alpha[static_cast<std::size_t>(index + 2)] = 1;
    }
    const std::size_t positive
        = static_cast<std::size_t>(middle_jobs + 2);
    alpha[positive] = 2;
    beta[positive] = 1;
    for (int index = 0; index < sink_jobs; ++index) {
        beta[positive + 1U + static_cast<std::size_t>(index)] = 2;
    }
    std::vector<int> degrees(vertices, 0);
    for (std::size_t index = 0U; index < vertices; ++index) {
        degrees[index] = alpha[index] + beta[index];
        if ((sign_mask & (1 << (degrees[index] - 1))) != 0) {
            signs[index] = -1;
        }
    }
    State state;
    state.alpha = alpha;
    state.order.resize(vertices);
    for (std::size_t index = 0U; index < vertices; ++index) {
        state.order[index] = index;
    }
    const std::size_t target = 1U;
    const Potential original = transfer_potential(
        state, degrees, signs, target
    );
    const Potential direct = direct_potential(
        state, degrees, signs, target
    );
    if (direct.delta != original.delta
        || direct.total != original.total) {
        throw std::runtime_error(
            "family direct and transfer potentials disagree"
        );
    }
    const std::vector<std::size_t> word = plateau_escape_word(
        state, degrees, signs, target, original
    );
    const SourceTotals totals = source_totals(
        degrees, signs, target
    );
    std::size_t minus_factors = 0U;
    for (std::size_t vertex = 0U; vertex < vertices; ++vertex) {
        if (vertex != target && signs[vertex] < 0) {
            ++minus_factors;
        }
    }
    std::cout
        << "SU2_BK_LONG_ESCAPE_FAMILY middle_jobs=" << middle_jobs
        << " vertices=" << vertices
        << " sign_mask=" << sign_mask
        << " minus_factors=" << minus_factors
        << " delta=" << original.delta
        << " total=" << original.total
        << " empty=" << totals.empty
        << " O=" << totals.odd
        << " P=" << totals.positive
        << " escape_length=" << word.size()
        << " escape_word=[";
    for (std::size_t index = 0U; index < word.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << word[index];
    }
    std::cout << "] neighbors={";
    for (std::size_t generator = 0U;
         generator + 1U < vertices; ++generator) {
        if (generator != 0U) {
            std::cout << ',';
        }
        const State neighbor = bk_neighbor(
            state, degrees, generator
        );
        const Potential transfer = transfer_potential(
            neighbor, degrees, signs, target
        );
        const Potential direct_neighbor = direct_potential(
            neighbor, degrees, signs, target
        );
        if (direct_neighbor.delta != transfer.delta
            || direct_neighbor.total != transfer.total) {
            throw std::runtime_error(
                "family neighbor direct and transfer disagree"
            );
        }
        std::cout << generator << ":(" << transfer.delta
                  << ',' << transfer.total << ')';
    }
    std::cout << "} result=PASS\n";
}

void test_random_composite_states(
    int trials_input,
    int vertices_input,
    int part_bound,
    std::uint64_t seed,
    bool require_source_deficit,
    bool force_composite,
    bool minimal_deficit_only
) {
    if (trials_input < 1 || vertices_input < 3 || part_bound < 1) {
        throw std::runtime_error("random search bounds are too small");
    }
    const std::size_t vertices
        = static_cast<std::size_t>(vertices_input);
    const std::uint64_t trials
        = static_cast<std::uint64_t>(trials_input);
    std::mt19937_64 generator(seed);
    std::uint64_t valid_paths = 0U;
    std::uint64_t tested_negative_states = 0U;
    std::uint64_t direct_escapes = 0U;
    std::uint64_t two_step_escapes = 0U;
    std::uint64_t longer_escapes = 0U;
    std::uint64_t prefix_forest_failures = 0U;
    std::uint64_t multiple_transfer_failures = 0U;
    std::uint64_t nonprimitive_negative_states = 0U;
    std::uint64_t no_normal_negative_excursion = 0U;
    std::uint64_t mixed_scheduling_direction_failures = 0U;
    std::uint64_t over_two_strict_scheduling_failures = 0U;
    std::uint64_t normal_zero_transfer_states = 0U;
    std::uint64_t normal_zero_transfer_failures = 0U;
    std::size_t maximum_strict_scheduling_steps = 0U;
    std::size_t maximum_escape_length = 0U;
    std::map<
        std::tuple<
            std::int64_t, std::uint64_t,
            std::int64_t, std::uint64_t
        >,
        std::uint64_t
    > two_step_factor_histogram;

    for (std::uint64_t trial = 0U; trial < trials; ++trial) {
        std::vector<int> alpha(vertices, 0);
        std::vector<int> beta(vertices, 0);
        int height = 0;
        bool valid_path = true;
        std::uniform_int_distribution<std::size_t> choose_split(
            1U, vertices - 2U
        );
        const std::size_t split = force_composite
            ? choose_split(generator) : vertices;
        for (std::size_t position = 0U;
             position < vertices; ++position) {
            const std::size_t segment_end
                = position < split ? split : vertices;
            const std::size_t remaining
                = segment_end - position - 1U;
            std::vector<std::pair<int, int>> choices;
            for (int bottom = 0;
                 bottom <= std::min(part_bound, height); ++bottom) {
                for (int top = 0; top <= part_bound; ++top) {
                    if (top == 0 && bottom == 0) {
                        continue;
                    }
                    const int next_height = height + top - bottom;
                    if (next_height
                        <= static_cast<int>(remaining) * part_bound) {
                        choices.emplace_back(top, bottom);
                    }
                }
            }
            if (choices.empty()) {
                valid_path = false;
                break;
            }
            std::uniform_int_distribution<std::size_t> choose(
                0U, choices.size() - 1U
            );
            const auto [top, bottom] = choices[choose(generator)];
            alpha[position] = top;
            beta[position] = bottom;
            height += top - bottom;
            if (force_composite
                && position + 1U == split && height != 0) {
                valid_path = false;
                break;
            }
        }
        if (!valid_path || height != 0) {
            continue;
        }
        ++valid_paths;

        std::vector<int> degrees(vertices, 0);
        for (std::size_t index = 0U; index < vertices; ++index) {
            degrees[index] = alpha[index] + beta[index];
        }
        State state;
        state.alpha = alpha;
        state.order.resize(vertices);
        for (std::size_t index = 0U; index < vertices; ++index) {
            state.order[index] = index;
        }
        if (force_composite
            && primitive_queue_path(state, degrees)) {
            throw std::runtime_error(
                "forced-composite random path is primitive"
            );
        }
        std::uniform_int_distribution<std::size_t> choose_target(
            0U, vertices - 1U
        );
        constexpr std::size_t assignment_attempts = 32U;
        for (std::size_t assignment = 0U;
             assignment < assignment_attempts; ++assignment) {
            const std::size_t target = choose_target(generator);
            std::vector<int> distinct_degrees;
            for (std::size_t index = 0U; index < vertices; ++index) {
                if (index == target) {
                    continue;
                }
                if (std::find(
                        distinct_degrees.begin(),
                        distinct_degrees.end(), degrees[index]
                    ) == distinct_degrees.end()) {
                    distinct_degrees.push_back(degrees[index]);
                }
            }
            if (distinct_degrees.empty()
                || distinct_degrees.size()
                    >= std::numeric_limits<std::uint64_t>::digits) {
                continue;
            }
            std::uniform_int_distribution<std::uint64_t> choose_signs(
                1U,
                (std::uint64_t{1} << distinct_degrees.size()) - 1U
            );
            const std::uint64_t sign_mask = choose_signs(generator);
            std::vector<int> signs(vertices, 1);
            std::size_t minus_factors = 0U;
            for (std::size_t index = 0U; index < vertices; ++index) {
                if (index == target) {
                    continue;
                }
                const auto found = std::find(
                    distinct_degrees.begin(),
                    distinct_degrees.end(), degrees[index]
                );
                const std::size_t sign_index
                    = static_cast<std::size_t>(std::distance(
                        distinct_degrees.begin(), found
                    ));
                if ((sign_mask
                     & (std::uint64_t{1} << sign_index)) != 0U) {
                    signs[index] = -1;
                    ++minus_factors;
                }
            }
            if ((minus_factors & 1U) != 0U) {
                continue;
            }
            const Potential original = transfer_potential(
                state, degrees, signs, target
            );
            if (original.delta >= 0) {
                continue;
            }
            if (minimal_deficit_only && original.delta != -2) {
                continue;
            }
            if (require_source_deficit) {
                const SourceTotals totals = source_totals(
                    degrees, signs, target
                );
                if (totals.odd <= totals.positive) {
                    continue;
                }
            }
            ++tested_negative_states;
            if (require_source_deficit
                && !primitive_queue_path(state, degrees)) {
                ++nonprimitive_negative_states;
                if (!has_normal_escaping_negative_excursion(
                        state, degrees, signs, target
                    )) {
                    ++no_normal_negative_excursion;
                    report_failure(
                        state, degrees, signs, target, original
                    );
                }
            }
            const int depth = generic_two_step_escape_depth(
                state, degrees, signs, target, original
            );
            if (depth != 1) {
                ++normal_zero_transfer_states;
                if (!has_normal_zero_transfer_escape(
                        state, degrees, signs, target, original
                    )) {
                    ++normal_zero_transfer_failures;
                    report_failure(
                        state, degrees, signs, target, original
                    );
                }
            }
            if (depth == 1) {
                ++direct_escapes;
                maximum_escape_length = std::max(
                    maximum_escape_length, std::size_t{1}
                );
            } else if (depth == 2) {
                ++two_step_escapes;
                maximum_escape_length = std::max(
                    maximum_escape_length, std::size_t{2}
                );
                if (original.delta == -2
                    && return_boundaries(state, degrees).size() == 1U) {
                    ++two_step_factor_histogram[
                        two_excursion_factor_signature(
                            state, degrees, signs, target
                        )
                    ];
                }
                if (two_step_escapes <= 5U) {
                    const std::vector<std::size_t> word
                        = plateau_escape_word(
                            state, degrees, signs, target, original
                        );
                    std::vector<int> displayed_signs = signs;
                    displayed_signs[target] = 0;
                    std::cout
                        << "SU2_BK_RANDOM_COMPOSITE_TWO_STEP alpha=";
                    print_vector(alpha);
                    std::cout << " beta=";
                    print_vector(beta);
                    std::cout << " signs=";
                    print_vector(displayed_signs);
                    std::cout << " target=" << target
                              << " delta=" << original.delta
                              << " total=" << original.total
                              << " word=[";
                    for (std::size_t index = 0U;
                         index < word.size(); ++index) {
                        if (index != 0U) {
                            std::cout << ',';
                        }
                        std::cout << word[index];
                    }
                    std::cout << "]\n";
                }
            } else {
                const std::vector<std::size_t> word
                    = plateau_escape_word(
                        state, degrees, signs, target, original
                    );
                std::size_t positive_transfers = 0U;
                std::size_t strict_scheduling_steps = 0U;
                bool has_positive_direction = false;
                bool has_negative_direction = false;
                State prefix_state = state;
                for (std::size_t index = 0U;
                     index + 1U < word.size(); ++index) {
                    const int transfer = bk_transfer_amount(
                        prefix_state, degrees, word[index]
                    );
                    if (transfer > 0) {
                        ++positive_transfers;
                    } else {
                        const int direction = scheduling_move_direction(
                            prefix_state, degrees, word[index]
                        );
                        if (direction > 0) {
                            has_positive_direction = true;
                            ++strict_scheduling_steps;
                        } else if (direction < 0) {
                            has_negative_direction = true;
                            ++strict_scheduling_steps;
                        }
                    }
                    prefix_state = bk_neighbor(
                        prefix_state, degrees, word[index]
                    );
                }
                const bool prefix_forest
                    = square_free_proper_prefix(word);
                if (!prefix_forest) {
                    ++prefix_forest_failures;
                }
                if (positive_transfers > 1U) {
                    ++multiple_transfer_failures;
                }
                if (has_positive_direction
                    && has_negative_direction) {
                    ++mixed_scheduling_direction_failures;
                }
                if (strict_scheduling_steps > 2U) {
                    ++over_two_strict_scheduling_failures;
                }
                maximum_strict_scheduling_steps = std::max(
                    maximum_strict_scheduling_steps,
                    strict_scheduling_steps
                );
                if (word.empty() || word.size() > vertices
                    || !prefix_forest
                    || positive_transfers > 1U
                    || (has_positive_direction
                        && has_negative_direction)
                    || strict_scheduling_steps > 2U) {
                    report_failure(
                        state, degrees, signs, target, original
                    );
                }
                ++longer_escapes;
                maximum_escape_length = std::max(
                    maximum_escape_length, word.size()
                );
            }
        }
    }

    std::cout
        << "SU2_BK_RANDOM_COMPOSITE source_only="
        << (require_source_deficit ? 1 : 0)
        << " forced_composite=" << (force_composite ? 1 : 0)
        << " minimal_deficit_only="
        << (minimal_deficit_only ? 1 : 0)
        << " trials=" << trials
        << " vertices=" << vertices
        << " maximum_queue_part=" << part_bound
        << " seed=" << seed
        << " valid_paths=" << valid_paths
        << " tested_negative_states="
        << tested_negative_states
        << " direct_escapes=" << direct_escapes
        << " two_step_escapes=" << two_step_escapes
        << " longer_escapes=" << longer_escapes
        << " prefix_forest_failures="
        << prefix_forest_failures
        << " multiple_transfer_failures="
        << multiple_transfer_failures
        << " nonprimitive_negative_states="
        << nonprimitive_negative_states
        << " no_normal_negative_excursion="
        << no_normal_negative_excursion
        << " mixed_scheduling_direction_failures="
        << mixed_scheduling_direction_failures
        << " over_two_strict_scheduling_failures="
        << over_two_strict_scheduling_failures
        << " normal_zero_transfer_states="
        << normal_zero_transfer_states
        << " normal_zero_transfer_failures="
        << normal_zero_transfer_failures
        << " maximum_strict_scheduling_steps="
        << maximum_strict_scheduling_steps
        << " maximum_escape_length="
        << maximum_escape_length
        << " two_step_factor_histogram={";
    bool first_factor = true;
    for (const auto& [factor, multiplicity] :
         two_step_factor_histogram) {
        if (!first_factor) {
            std::cout << ',';
        }
        first_factor = false;
        std::cout
            << '(' << std::get<0>(factor)
            << ',' << std::get<1>(factor)
            << ';' << std::get<2>(factor)
            << ',' << std::get<3>(factor)
            << "):" << multiplicity;
    }
    std::cout << "} result=PASS\n";
}

void analyze_closed_core_order_slice() {
    const std::vector<int> degrees{1,1,1,1,1,3,2,2,2};
    const std::vector<int> signs{1,1,1,1,1,-1,-1,-1,-1};
    const std::size_t target = 0U;
    State state;
    state.alpha.resize(degrees.size(), 0);
    state.order.resize(degrees.size(), 0U);
    for (std::size_t index = 0U; index < degrees.size(); ++index) {
        state.order[index] = index;
    }
    std::uint64_t states = 0U;
    std::uint64_t negative_states = 0U;
    std::int64_t delta_sum = 0;
    std::map<std::int64_t, std::uint64_t> histogram;
    const std::function<void(std::size_t, int)> enumerate
        = [&](std::size_t position, int height) {
            if (position == degrees.size()) {
                if (height != 0) {
                    return;
                }
                const Potential value = transfer_potential(
                    state, degrees, signs, target
                );
                ++states;
                negative_states += value.delta < 0 ? 1U : 0U;
                delta_sum += value.delta;
                ++histogram[value.delta];
                return;
            }
            for (int alpha = 0;
                 alpha <= degrees[position]; ++alpha) {
                const int beta = degrees[position] - alpha;
                if (height < beta) {
                    continue;
                }
                state.alpha[position] = alpha;
                enumerate(
                    position + 1U,
                    height + alpha - beta
                );
            }
        };
    enumerate(0U, 0);
    const SourceTotals totals = source_totals(
        degrees, signs, target
    );
    std::cout
        << "SU2_BK_CLOSED_CORE_SLICE states=" << states
        << " negative_states=" << negative_states
        << " delta_sum=" << delta_sum
        << " empty=" << totals.empty
        << " O=" << totals.odd
        << " P=" << totals.positive
        << " histogram={";
    bool first = true;
    for (const auto& [delta, multiplicity] : histogram) {
        if (!first) {
            std::cout << ',';
        }
        first = false;
        std::cout << delta << ':' << multiplicity;
    }
    const std::int64_t expected
        = static_cast<std::int64_t>(totals.empty)
          + static_cast<std::int64_t>(totals.positive)
          - static_cast<std::int64_t>(totals.odd);
    std::cout << "} result="
              << (states == totals.empty && delta_sum == expected
                  ? "PASS" : "FAIL")
              << '\n';
}

void test_rigid_excursion_embedding() {
    const std::vector<int> alpha{
        2, 3, 3, 2, 0, 0, 0, 0,
        1, 0
    };
    const std::vector<int> beta{
        0, 1, 0, 0, 1, 2, 3, 3,
        0, 1
    };
    const std::vector<int> signs{
        1, 1, 1, 1, 1, 1, -1, -1,
        1, 1
    };
    const std::size_t target = 2U;
    std::vector<int> degrees(alpha.size(), 0);
    for (std::size_t index = 0U; index < alpha.size(); ++index) {
        degrees[index] = alpha[index] + beta[index];
    }
    State state;
    state.alpha = alpha;
    state.order.resize(alpha.size());
    for (std::size_t index = 0U; index < alpha.size(); ++index) {
        state.order[index] = index;
    }
    const Potential original = transfer_potential(
        state, degrees, signs, target
    );
    const SourceTotals totals = source_totals(
        degrees, signs, target
    );
    const std::vector<std::size_t> word = plateau_escape_word(
        state, degrees, signs, target, original
    );
    const auto factor = two_excursion_factor_signature(
        state, degrees, signs, target
    );
    std::cout
        << "SU2_BK_RIGID_EXCURSION_EMBEDDING delta="
        << original.delta
        << " total=" << original.total
        << " O=" << totals.odd
        << " P=" << totals.positive
        << " factors=("
        << std::get<0>(factor) << ',' << std::get<1>(factor)
        << ';' << std::get<2>(factor) << ',' << std::get<3>(factor)
        << ") escape_word=[";
    for (std::size_t index = 0U; index < word.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << word[index];
    }
    std::cout << "] result=PASS\n";
}

void scan_rigid_excursion_embeddings(int maximum_label) {
    if (maximum_label < 1) {
        throw std::runtime_error(
            "embedding scan maximum label is too small"
        );
    }
    const std::vector<int> base_alpha{
        2, 3, 3, 2, 0, 0, 0, 0
    };
    const std::vector<int> base_beta{
        0, 1, 0, 0, 1, 2, 3, 3
    };
    const std::vector<int> base_signs{
        1, 1, 1, 1, 1, 1, -1, -1
    };
    const std::size_t target = 2U;
    std::uint64_t tested = 0U;
    std::uint64_t source_deficit = 0U;
    std::uint64_t source_deficit_beyond_two = 0U;
    for (int label = 1; label <= maximum_label; ++label) {
        for (const int appended_sign : {-1, 1}) {
            bool sign_compatible = true;
            for (std::size_t index = 0U;
                 index < base_alpha.size(); ++index) {
                if (index == target) {
                    continue;
                }
                const int degree
                    = base_alpha[index] + base_beta[index];
                if (degree == label
                    && base_signs[index] != appended_sign) {
                    sign_compatible = false;
                    break;
                }
            }
            if (!sign_compatible) {
                continue;
            }
            std::vector<int> alpha = base_alpha;
            std::vector<int> beta = base_beta;
            std::vector<int> signs = base_signs;
            alpha.push_back(label);
            alpha.push_back(0);
            beta.push_back(0);
            beta.push_back(label);
            signs.push_back(appended_sign);
            signs.push_back(appended_sign);
            std::vector<int> degrees(alpha.size(), 0);
            for (std::size_t index = 0U;
                 index < alpha.size(); ++index) {
                degrees[index] = alpha[index] + beta[index];
            }
            State state;
            state.alpha = alpha;
            state.order.resize(alpha.size());
            for (std::size_t index = 0U;
                 index < alpha.size(); ++index) {
                state.order[index] = index;
            }
            const Potential original = transfer_potential(
                state, degrees, signs, target
            );
            const SourceTotals totals = source_totals(
                degrees, signs, target
            );
            const int depth = generic_two_step_escape_depth(
                state, degrees, signs, target, original
            );
            const std::vector<std::size_t> word
                = plateau_escape_word(
                    state, degrees, signs, target, original
                );
            ++tested;
            const bool is_source_deficit
                = totals.odd > totals.positive;
            if (is_source_deficit) {
                ++source_deficit;
                if (depth == 0) {
                    ++source_deficit_beyond_two;
                }
            }
            std::cout
                << "SU2_BK_RIGID_EMBEDDING label=" << label
                << " sign=" << appended_sign
                << " O=" << totals.odd
                << " P=" << totals.positive
                << " source_deficit="
                << (is_source_deficit ? 1 : 0)
                << " two_step_depth=" << depth
                << " escape_length=" << word.size()
                << "\n";
        }
    }
    std::cout
        << "SU2_BK_RIGID_EMBEDDING_SCAN maximum_label="
        << maximum_label
        << " tested=" << tested
        << " source_deficit=" << source_deficit
        << " source_deficit_beyond_two="
        << source_deficit_beyond_two
        << " result="
        << (source_deficit_beyond_two == 0U ? "PASS" : "FAIL")
        << "\n";
    if (source_deficit_beyond_two != 0U) {
        std::exit(EXIT_FAILURE);
    }
}

void test_composite_two_step_counterexample() {
    const std::vector<int> alpha{
        2, 3, 3, 2, 0, 0, 0, 0, 3, 0
    };
    const std::vector<int> beta{
        0, 1, 0, 0, 1, 2, 3, 3, 0, 3
    };
    const std::vector<int> signs{
        1, 1, 1, 1, 1, 1, -1, -1, -1, -1
    };
    const std::size_t target = 2U;
    std::vector<int> degrees(alpha.size(), 0);
    for (std::size_t index = 0U; index < alpha.size(); ++index) {
        degrees[index] = alpha[index] + beta[index];
    }
    State state;
    state.alpha = alpha;
    state.order.resize(alpha.size());
    for (std::size_t index = 0U; index < alpha.size(); ++index) {
        state.order[index] = index;
    }
    const Potential original = transfer_potential(
        state, degrees, signs, target
    );
    const Potential direct = direct_potential(
        state, degrees, signs, target
    );
    const SourceTotals totals = source_totals(
        degrees, signs, target
    );
    const int two_step_depth = generic_two_step_escape_depth(
        state, degrees, signs, target, original
    );
    const std::vector<std::size_t> word = plateau_escape_word(
        state, degrees, signs, target, original
    );
    const auto factor = two_excursion_factor_signature(
        state, degrees, signs, target
    );
    bool prefix_is_plateau = true;
    State current = state;
    std::size_t positive_transfers = 0U;
    std::vector<int> scheduling_directions;
    std::vector<Potential> word_potentials;
    for (std::size_t index = 0U; index < word.size(); ++index) {
        if (index + 1U < word.size()) {
            const int transfer = bk_transfer_amount(
                current, degrees, word[index]
            );
            if (transfer > 0) {
                ++positive_transfers;
            } else {
                scheduling_directions.push_back(
                    scheduling_move_direction(
                        current, degrees, word[index]
                    )
                );
            }
        }
        current = bk_neighbor(current, degrees, word[index]);
        word_potentials.push_back(transfer_potential(
            current, degrees, signs, target
        ));
        if (index + 1U < word.size()
            && !plateau(word_potentials.back(), original)) {
            prefix_is_plateau = false;
        }
    }
    if (direct.delta != original.delta
        || direct.total != original.total
        || original.delta != -2
        || totals.odd <= totals.positive
        || primitive_queue_path(state, degrees)
        || two_step_depth != 0
        || word.size() != 4U
        || !prefix_is_plateau
        || word_potentials.empty()
        || !lower(word_potentials.back(), original)) {
        throw std::runtime_error(
            "composite two-step counterexample audit failed"
        );
    }
    std::cout
        << "SU2_BK_COMPOSITE_TWO_STEP_COUNTEREXAMPLE alpha=";
    print_vector(alpha);
    std::cout << " beta=";
    print_vector(beta);
    std::vector<int> displayed_signs = signs;
    displayed_signs[target] = 0;
    std::cout << " signs=";
    print_vector(displayed_signs);
    std::cout
        << " target=" << target
        << " delta=" << original.delta
        << " total=" << original.total
        << " O=" << totals.odd
        << " P=" << totals.positive
        << " factors=("
        << std::get<0>(factor) << ',' << std::get<1>(factor)
        << ';' << std::get<2>(factor) << ',' << std::get<3>(factor)
        << ") two_step_depth=" << two_step_depth
        << " positive_prefix_transfers=" << positive_transfers
        << " scheduling_directions=[";
    for (std::size_t index = 0U;
         index < scheduling_directions.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << scheduling_directions[index];
    }
    std::cout << "] word=[";
    for (std::size_t index = 0U; index < word.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << word[index];
    }
    std::cout << "] word_potentials=[";
    for (std::size_t index = 0U;
         index < word_potentials.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << '(' << word_potentials[index].delta
                  << ',' << word_potentials[index].total << ')';
    }
    std::cout << "] result=PASS\n";
}

void scan_closed_rigid_source_totals(
    int maximum_label,
    int maximum_pairs
) {
    if (maximum_label < 1 || maximum_pairs < 1) {
        throw std::runtime_error(
            "closed embedding scan maximum label is too small"
        );
    }
    const std::vector<int> base_alpha{
        1, 1, 1, 1, 1, 2, 0, 0, 0
    };
    const std::vector<int> base_beta{
        0, 0, 0, 0, 0, 1, 2, 2, 2
    };
    const std::vector<int> base_signs{
        1, 1, 1, 1, 1, -1, -1, -1, -1
    };
    const std::size_t target = 0U;
    std::uint64_t tested = 0U;
    std::uint64_t source_deficit = 0U;
    std::uint64_t source_without_two_step = 0U;
    for (int pair_count = 1;
         pair_count <= maximum_pairs; ++pair_count) {
        for (int label = 1; label <= maximum_label; ++label) {
            for (const int appended_sign : {-1, 1}) {
            bool sign_compatible = true;
            for (std::size_t index = 0U;
                 index < base_alpha.size(); ++index) {
                if (index == target) {
                    continue;
                }
                const int degree
                    = base_alpha[index] + base_beta[index];
                if (degree == label
                    && base_signs[index] != appended_sign) {
                    sign_compatible = false;
                    break;
                }
            }
            if (!sign_compatible) {
                continue;
            }
            std::vector<int> alpha = base_alpha;
            std::vector<int> beta = base_beta;
            std::vector<int> signs = base_signs;
            for (int pair = 0; pair < pair_count; ++pair) {
                alpha.push_back(label);
                alpha.push_back(0);
                beta.push_back(0);
                beta.push_back(label);
                signs.push_back(appended_sign);
                signs.push_back(appended_sign);
            }
            std::vector<int> degrees(alpha.size(), 0);
            for (std::size_t index = 0U;
                 index < alpha.size(); ++index) {
                degrees[index] = alpha[index] + beta[index];
            }
            State state;
            state.alpha = alpha;
            state.order.resize(alpha.size());
            for (std::size_t index = 0U;
                 index < alpha.size(); ++index) {
                state.order[index] = index;
            }
            const Potential original = transfer_potential(
                state, degrees, signs, target
            );
            const SourceTotals totals = source_totals(
                degrees, signs, target
            );
            const int depth = generic_two_step_escape_depth(
                state, degrees, signs, target, original
            );
            const bool is_source_deficit
                = totals.odd > totals.positive;
            ++tested;
            if (is_source_deficit) {
                ++source_deficit;
                if (depth == 0) {
                    ++source_without_two_step;
                }
            }
            std::cout
                << "SU2_BK_CLOSED_RIGID_SOURCE label=" << label
                << " sign=" << appended_sign
                << " pair_count=" << pair_count
                << " delta=" << original.delta
                << " total=" << original.total
                << " O=" << totals.odd
                << " P=" << totals.positive
                << " source_deficit="
                << (is_source_deficit ? 1 : 0)
                << " two_step_depth=" << depth
                << "\n";
            }
        }
    }
    std::cout
        << "SU2_BK_CLOSED_RIGID_SOURCE_SCAN maximum_label="
        << maximum_label
        << " maximum_pairs=" << maximum_pairs
        << " tested=" << tested
        << " source_deficit=" << source_deficit
        << " source_without_two_step="
        << source_without_two_step
        << " result=PASS\n";
}

void scan_closed_mixed_rigid_source_totals(
    int maximum_label,
    int pair_count
) {
    if (maximum_label < 1 || pair_count < 1) {
        throw std::runtime_error(
            "mixed padding scan bounds are too small"
        );
    }
    const std::vector<int> base_alpha{
        1, 1, 1, 1, 1, 2, 0, 0, 0
    };
    const std::vector<int> base_beta{
        0, 0, 0, 0, 0, 1, 2, 2, 2
    };
    const std::vector<int> base_signs{
        1, 1, 1, 1, 1, -1, -1, -1, -1
    };
    const std::size_t target = 0U;
    std::vector<std::pair<int, int>> options;
    for (int label = 1; label <= maximum_label; ++label) {
        bool has_fixed_sign = false;
        int fixed_sign = 1;
        for (std::size_t index = 0U;
             index < base_alpha.size(); ++index) {
            if (index == target) {
                continue;
            }
            const int degree
                = base_alpha[index] + base_beta[index];
            if (degree == label) {
                has_fixed_sign = true;
                fixed_sign = base_signs[index];
                break;
            }
        }
        if (has_fixed_sign) {
            options.emplace_back(label, fixed_sign);
        } else {
            options.emplace_back(label, -1);
            options.emplace_back(label, 1);
        }
    }

    std::vector<std::pair<int, int>> selected;
    std::uint64_t tested = 0U;
    std::uint64_t source_deficit = 0U;
    const std::function<void(std::size_t)> enumerate
        = [&](std::size_t first_option) {
            if (selected.size()
                == static_cast<std::size_t>(pair_count)) {
                std::vector<int> alpha = base_alpha;
                std::vector<int> beta = base_beta;
                std::vector<int> signs = base_signs;
                for (const auto& [label, sign] : selected) {
                    alpha.push_back(label);
                    alpha.push_back(0);
                    beta.push_back(0);
                    beta.push_back(label);
                    signs.push_back(sign);
                    signs.push_back(sign);
                }
                std::vector<int> degrees(alpha.size(), 0);
                for (std::size_t index = 0U;
                     index < alpha.size(); ++index) {
                    degrees[index] = alpha[index] + beta[index];
                }
                const SourceTotals totals = source_totals(
                    degrees, signs, target
                );
                ++tested;
                if (totals.odd > totals.positive) {
                    ++source_deficit;
                    std::cout
                        << "SU2_BK_CLOSED_MIXED_SOURCE choices=[";
                    for (std::size_t index = 0U;
                         index < selected.size(); ++index) {
                        if (index != 0U) {
                            std::cout << ',';
                        }
                        std::cout << selected[index].first
                                  << '^'
                                  << selected[index].second;
                    }
                    std::cout
                        << "] O=" << totals.odd
                        << " P=" << totals.positive
                        << "\n";
                }
                return;
            }
            for (std::size_t option = first_option;
                 option < options.size(); ++option) {
                const auto [label, sign] = options[option];
                bool compatible = true;
                for (const auto& [old_label, old_sign] : selected) {
                    if (old_label == label && old_sign != sign) {
                        compatible = false;
                        break;
                    }
                }
                if (!compatible) {
                    continue;
                }
                selected.emplace_back(label, sign);
                enumerate(option);
                selected.pop_back();
            }
        };
    enumerate(0U);
    std::cout
        << "SU2_BK_CLOSED_MIXED_SOURCE_SCAN maximum_label="
        << maximum_label
        << " pair_count=" << pair_count
        << " tested=" << tested
        << " source_deficit=" << source_deficit
        << " result=" << (source_deficit == 0U ? "PASS" : "FAIL")
        << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            throw std::runtime_error(
                "missing search mode or vertex count"
            );
        }
        if (argc == 2 && std::string(argv[1]) == "embedding") {
            test_rigid_excursion_embedding();
            return EXIT_SUCCESS;
        }
        if (argc == 3 && std::string(argv[1]) == "embedding-scan") {
            scan_rigid_excursion_embeddings(
                parse_nonnegative(argv[2], "maximum label")
            );
            return EXIT_SUCCESS;
        }
        if (argc == 2
            && std::string(argv[1]) == "composite-counterexample") {
            test_composite_two_step_counterexample();
            return EXIT_SUCCESS;
        }
        if (argc == 2
            && std::string(argv[1]) == "closed-core-slice") {
            analyze_closed_core_order_slice();
            return EXIT_SUCCESS;
        }
        if (argc == 3
            && std::string(argv[1]) == "closed-embedding-scan") {
            scan_closed_rigid_source_totals(
                parse_nonnegative(argv[2], "maximum label"),
                1
            );
            return EXIT_SUCCESS;
        }
        if (argc == 4
            && std::string(argv[1]) == "closed-padding-scan") {
            scan_closed_rigid_source_totals(
                parse_nonnegative(argv[2], "maximum label"),
                parse_nonnegative(argv[3], "maximum pair count")
            );
            return EXIT_SUCCESS;
        }
        if (argc == 4
            && std::string(argv[1]) == "closed-mixed-padding-scan") {
            scan_closed_mixed_rigid_source_totals(
                parse_nonnegative(argv[2], "maximum label"),
                parse_nonnegative(argv[3], "pair count")
            );
            return EXIT_SUCCESS;
        }
        if ((argc == 3 || argc == 4)
            && std::string(argv[1]) == "family") {
            test_long_escape_family(
                parse_nonnegative(argv[2], "middle-job count"),
                argc == 4
                    ? parse_nonnegative(argv[3], "sign mask")
                    : 3
            );
            return EXIT_SUCCESS;
        }
        if (argc == 6
            && (std::string(argv[1]) == "random-composite"
                || std::string(argv[1]) == "random-composite-any"
                || std::string(argv[1]) == "random-source")) {
            const std::string random_mode(argv[1]);
            test_random_composite_states(
                parse_nonnegative(argv[2], "trial count"),
                parse_nonnegative(argv[3], "vertex count"),
                parse_nonnegative(argv[4], "maximum queue part"),
                static_cast<std::uint64_t>(
                    parse_nonnegative(argv[5], "seed")
                ),
                random_mode != "random-composite-any",
                random_mode != "random-source",
                random_mode == "random-composite-any"
            );
            return EXIT_SUCCESS;
        }
        const std::string mode(argv[1]);
        const bool unrestricted_sharded
            = argc == 6 && mode == "shard";
        const bool strict_three
            = mode == "source3" || mode == "source3-shard";
        const bool strict_forest
            = mode == "source-forest"
              || mode == "source-forest-shard";
        const bool source_only
            = (argc == 4
               && (mode == "source" || mode == "source3"
                   || mode == "source-forest"))
              || (argc == 6
                  && (mode == "source-shard"
                      || mode == "source3-shard"
                      || mode == "source-forest-shard"));
        const bool sharded
            = unrestricted_sharded
              || (argc == 6
                  && (mode == "source-shard"
                      || mode == "source3-shard"
                      || mode == "source-forest-shard"));
        if ((!source_only && !unrestricted_sharded && argc != 3)
            || (source_only && !sharded && argc != 4)
            || (sharded && argc != 6)) {
            throw std::runtime_error(
                "usage: search_su2_bk_local_escape"
                " vertices maximum_queue_part |"
                " source vertices maximum_queue_part |"
                " source3 vertices maximum_queue_part |"
                " source-forest vertices maximum_queue_part |"
                " shard vertices maximum_queue_part"
                " shard_count shard_index |"
                " source-shard vertices maximum_queue_part"
                " shard_count shard_index |"
                " source3-shard vertices maximum_queue_part"
                " shard_count shard_index |"
                " source-forest-shard vertices maximum_queue_part"
                " shard_count shard_index |"
                " family even_middle_jobs [sign_mask] |"
                " embedding | embedding-scan maximum_label |"
                " composite-counterexample |"
                " closed-core-slice |"
                " closed-embedding-scan maximum_label |"
                " closed-padding-scan maximum_label maximum_pairs |"
                " closed-mixed-padding-scan maximum_label pair_count |"
                " random-composite trials vertices"
                " maximum_queue_part seed |"
                " random-composite-any trials vertices"
                " maximum_queue_part seed |"
                " random-source trials vertices"
                " maximum_queue_part seed"
            );
        }
        const int argument_offset
            = (source_only || unrestricted_sharded) ? 2 : 1;
        const int vertices_input = parse_nonnegative(
            argv[argument_offset], "vertex count"
        );
        const int part_bound = parse_nonnegative(
            argv[argument_offset + 1], "maximum queue part"
        );
        if (vertices_input < 2 || part_bound < 1) {
            throw std::runtime_error("search bounds are too small");
        }
        const int shard_count_input = sharded
            ? parse_nonnegative(argv[4], "shard count") : 1;
        const int shard_index_input = sharded
            ? parse_nonnegative(argv[5], "shard index") : 0;
        if (shard_count_input < 1
            || shard_index_input >= shard_count_input) {
            throw std::runtime_error("invalid shard specification");
        }
        const std::uint64_t shard_count
            = static_cast<std::uint64_t>(shard_count_input);
        const std::uint64_t shard_index
            = static_cast<std::uint64_t>(shard_index_input);
        if (strict_three) {
            maximum_source_escape = 3U;
        }
        require_almost_square_free = strict_forest;
        const std::size_t vertices
            = static_cast<std::size_t>(vertices_input);
        std::vector<int> alpha(vertices, 0);
        std::vector<int> beta(vertices, 0);
        SearchCounts counts;
        enumerate_queue_paths(
            0U, 0, part_bound, alpha, beta, counts, source_only,
            shard_count, shard_index
        );
        std::cout
            << "SU2_BK_LOCAL_ESCAPE source_only="
            << (source_only ? 1 : 0)
            << " vertices=" << vertices
            << " maximum_queue_part=" << part_bound
            << " queue_paths=" << counts.queue_paths
            << " selected_queue_paths="
            << counts.selected_queue_paths
            << " shard_count=" << shard_count
            << " shard_index=" << shard_index
            << " signed_states=" << counts.signed_states
            << " negative_states=" << counts.negative_states
            << " direct_escapes=" << counts.direct_escapes
            << " two_step_escapes=" << counts.two_step_escapes
            << " three_step_braid_escapes="
            << counts.three_step_braid_escapes
            << " longer_plateau_escapes="
            << counts.longer_plateau_escapes
            << " source_deficit_beyond_two_steps="
            << counts.source_deficit_beyond_two_steps
            << " source_forest_word_failures="
            << counts.source_forest_word_failures
            << " source_prefix_forest_word_failures="
            << counts.source_prefix_forest_word_failures
            << " source_nonzero_transfer_prefix_states="
            << counts.source_nonzero_transfer_prefix_states
            << " source_nonzero_transfer_prefix_steps="
            << counts.source_nonzero_transfer_prefix_steps
            << " source_multiple_transfer_prefix_states="
            << counts.source_multiple_transfer_prefix_states
            << " source_prefix_threshold_tie_steps="
            << counts.source_prefix_threshold_tie_steps
            << " source_prefix_positive_scheduling_steps="
            << counts.source_prefix_positive_scheduling_steps
            << " source_prefix_negative_scheduling_steps="
            << counts.source_prefix_negative_scheduling_steps
            << " source_prefix_mixed_scheduling_states="
            << counts.source_prefix_mixed_scheduling_states
            << " source_prefix_multiple_strict_scheduling_states="
            << counts.source_prefix_multiple_strict_scheduling_states
            << " source_maximum_strict_scheduling_steps="
            << counts.source_maximum_strict_scheduling_steps
            << " source_extremal_shell_states="
            << counts.source_extremal_shell_states
            << " source_unidirectional_shell_escapes="
            << counts.source_unidirectional_shell_escapes
            << " source_extremal_shell_escapes="
            << counts.source_extremal_shell_escapes
            << " source_extremal_sink_escapes="
            << counts.source_extremal_sink_escapes
            << " source_extremal_source_escapes="
            << counts.source_extremal_source_escapes
            << " source_extremal_source_shrink_escapes="
            << counts.source_extremal_source_shrink_escapes
            << " source_extremal_source_enlarge_escapes="
            << counts.source_extremal_source_enlarge_escapes
            << " source_extremal_source_positive_escapes="
            << counts.source_extremal_source_positive_escapes
            << " source_extremal_shell_obstructions="
            << counts.source_extremal_shell_obstructions
            << " source_extremal_positive_plateau_states="
            << counts.source_extremal_positive_plateau_states
            << " source_normal_zero_transfer_states="
            << counts.source_normal_zero_transfer_states
            << " source_normal_zero_transfer_escapes="
            << counts.source_normal_zero_transfer_escapes
            << " source_normal_zero_transfer_failures="
            << counts.source_normal_zero_transfer_failures
            << " source_maximum_zero_transfer_component="
            << counts.source_maximum_zero_transfer_component
            << " source_maximum_threshold_tie_quotient="
            << counts.source_maximum_threshold_tie_quotient
            << " hard_extremal_shell_states="
            << counts.hard_extremal_shell_states
            << " hard_extremal_source_shrink_escapes="
            << counts.hard_extremal_source_shrink_escapes
            << " hard_extremal_source_shrink_obstructions="
            << counts.hard_extremal_source_shrink_obstructions
            << " source_nonprimitive_beyond_two_steps="
            << counts.source_nonprimitive_beyond_two_steps
            << " nonprimitive_beyond_two_steps="
            << counts.nonprimitive_beyond_two_steps
            << " source_nonprimitive_negative_states="
            << counts.source_nonprimitive_negative_states
            << " source_return_shell_direct_escapes="
            << counts.source_return_shell_direct_escapes
            << " source_return_shell_two_step_escapes="
            << counts.source_return_shell_two_step_escapes
            << " source_return_shell_failures="
            << counts.source_return_shell_failures
            << " source_composite_direct_escapes="
            << counts.source_composite_direct_escapes
            << " source_composite_two_step_escapes="
            << counts.source_composite_two_step_escapes
            << " source_composite_escape_failures="
            << counts.source_composite_escape_failures
            << " source_composite_two_step_return_shell="
            << counts.source_composite_two_step_return_shell
            << " source_composite_no_normal_negative_excursion="
            << counts
                .source_composite_no_normal_negative_excursion
            << " trivial_sector_beyond_two_steps="
            << counts.trivial_sector_beyond_two_steps
            << " closed_odd_minus_plateaus="
            << counts.closed_odd_minus_plateaus
            << " closed_even_minus_trivial_plateaus="
            << counts.closed_even_minus_trivial_plateaus
            << " maximum_escape_length="
            << counts.maximum_escape_length
            << " longer_escape_word_histogram={";
        bool first_word = true;
        for (const auto& [word, multiplicity]
             : counts.longer_escape_word_histogram) {
            if (!first_word) {
                std::cout << ',';
            }
            first_word = false;
            std::cout << '[';
            for (std::size_t index = 0U; index < word.size(); ++index) {
                if (index != 0U) {
                    std::cout << '.';
                }
                std::cout << word[index];
            }
            std::cout << "]:" << multiplicity;
        }
        std::cout << "} source_deficit_escape_word_histogram={";
        first_word = true;
        for (const auto& [word, multiplicity]
             : counts.source_deficit_escape_word_histogram) {
            if (!first_word) {
                std::cout << ',';
            }
            first_word = false;
            std::cout << '[';
            for (std::size_t index = 0U; index < word.size(); ++index) {
                if (index != 0U) {
                    std::cout << '.';
                }
                std::cout << word[index];
            }
            std::cout << "]:" << multiplicity;
        }
        std::cout
            << "} hard_source_shrink_obstruction_margin_histogram={";
        bool first_margin = true;
        for (const auto& [margin, multiplicity] :
             counts.hard_source_shrink_obstruction_margin_histogram) {
            if (!first_margin) {
                std::cout << ',';
            }
            first_margin = false;
            std::cout << margin << ':' << multiplicity;
        }
        std::cout << "} source_composite_two_step_potential_histogram={";
        bool first_potential = true;
        for (const auto& [potential, multiplicity] :
             counts.source_composite_two_step_potential_histogram) {
            if (!first_potential) {
                std::cout << ',';
            }
            first_potential = false;
            std::cout << '(' << potential.first << ','
                      << potential.second << "):" << multiplicity;
        }
        std::cout << "} source_composite_two_step_factor_histogram={";
        bool first_factor = true;
        for (const auto& [factor, multiplicity] :
             counts.source_composite_two_step_factor_histogram) {
            if (!first_factor) {
                std::cout << ',';
            }
            first_factor = false;
            std::cout
                << '(' << std::get<0>(factor)
                << ',' << std::get<1>(factor)
                << ';' << std::get<2>(factor)
                << ',' << std::get<3>(factor)
                << "):" << multiplicity;
        }
        std::cout << "} first_hard_source_shrink_obstruction=";
        const HardSourceShrinkObstruction& obstruction
            = counts.first_hard_source_shrink_obstruction;
        if (!obstruction.present) {
            std::cout << "none";
        } else {
            std::cout << "{queue_path_index="
                      << obstruction.queue_path_index
                      << ",target_position=" << obstruction.target
                      << ",sign_mask=" << obstruction.sign_mask
                      << ",alpha=";
            print_vector(obstruction.alpha);
            std::cout << ",beta=";
            print_vector(obstruction.beta);
            std::cout << ",signs=";
            print_vector(obstruction.signs);
            std::cout << ",transfer_delta="
                      << obstruction.transfer.delta
                      << ",transfer_total="
                      << obstruction.transfer.total
                      << ",direct_delta="
                      << obstruction.direct.delta
                      << ",direct_total="
                      << obstruction.direct.total
                      << ",empty_sources="
                      << obstruction.totals.empty
                      << ",O=" << obstruction.totals.odd
                      << ",P=" << obstruction.totals.positive
                      << '}';
        }
        std::cout << "} result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
