#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct State {
    int first;
    int second;
};

struct Preimage {
    std::vector<State> states;
    std::vector<char> moves;
    int cut_first;
    int cut_second;
};

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(value);
}

bool fuses(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= std::min(
            source + label,
            2 * level - source - label
        )
        && ((source + label + target) & 1) == 0;
}

std::vector<int> neighbors(int level, int label, int source) {
    std::vector<int> result;
    for (int target = 0; target <= level; ++target) {
        if (fuses(level, label, source, target)) {
            result.push_back(target);
        }
    }
    return result;
}

std::vector<int> flatten(
    int model,
    const std::vector<State>& states,
    const std::vector<char>& moves
) {
    std::vector<int> result;
    result.reserve(2U + 3U * moves.size());
    result.push_back(model);
    result.push_back(static_cast<int>(moves.size()));
    for (std::size_t index = 0; index < moves.size(); ++index) {
        result.push_back(states[index].first);
        result.push_back(states[index].second);
        result.push_back(moves[index] == 'Y' ? 1 : 0);
    }
    result.push_back(states.back().first);
    result.push_back(states.back().second);
    return result;
}

bool valid_path(
    int level,
    int label,
    int model,
    const std::vector<State>& states,
    const std::vector<char>& moves,
    int minimum_y,
    int maximum_y
) {
    if (states.size() != moves.size() + 1U) {
        return false;
    }
    const State expected_start{0, model == 0 ? 0 : level};
    const State expected_end{level, model == 0 ? 0 : level};
    if (states.front().first != expected_start.first
        || states.front().second != expected_start.second
        || states.back().first != expected_end.first
        || states.back().second != expected_end.second
        || moves.empty() || moves.back() != 'X') {
        return false;
    }
    int y_moves = 0;
    for (std::size_t index = 0; index < moves.size(); ++index) {
        const State& source = states[index];
        const State& target = states[index + 1U];
        if (moves[index] == 'X') {
            if (source.second != target.second
                || !fuses(level, label, source.first, target.first)) {
                return false;
            }
        } else {
            ++y_moves;
            if (source.first != target.first
                || !fuses(level, label, source.second, target.second)) {
                return false;
            }
        }
    }
    return (y_moves & 1) == 0
        && minimum_y <= y_moves && y_moves <= maximum_y;
}

std::uint64_t binomial(int n, int r) {
    if (r < 0 || r > n) {
        return 0U;
    }
    r = std::min(r, n - r);
    std::uint64_t result = 1U;
    for (int index = 1; index <= r; ++index) {
        std::uint64_t numerator =
            static_cast<std::uint64_t>(n - r + index);
        std::uint64_t denominator = static_cast<std::uint64_t>(index);
        const std::uint64_t first_gcd =
            std::gcd(result, denominator);
        result /= first_gcd;
        denominator /= first_gcd;
        const std::uint64_t second_gcd =
            std::gcd(numerator, denominator);
        numerator /= second_gcd;
        denominator /= second_gcd;
        if (denominator != 1U) {
            throw std::runtime_error("binomial reduction failed");
        }
        if (numerator != 0U
            && result
                > std::numeric_limits<std::uint64_t>::max() / numerator) {
            throw std::runtime_error("binomial coefficient exceeds uint64");
        }
        result *= numerator;
    }
    return result;
}

struct FlowEdge {
    int target;
    int reverse;
    std::uint64_t capacity;
};

class Dinic {
public:
    explicit Dinic(int vertex_count)
        : graph_(static_cast<std::size_t>(vertex_count)),
          level_(static_cast<std::size_t>(vertex_count), -1),
          next_(static_cast<std::size_t>(vertex_count), 0) {}

    void add_edge(int source, int target, std::uint64_t capacity) {
        const int source_reverse = static_cast<int>(
            graph_[static_cast<std::size_t>(target)].size()
        );
        const int target_reverse = static_cast<int>(
            graph_[static_cast<std::size_t>(source)].size()
        );
        graph_[static_cast<std::size_t>(source)].push_back(
            FlowEdge{target, source_reverse, capacity}
        );
        graph_[static_cast<std::size_t>(target)].push_back(
            FlowEdge{source, target_reverse, 0U}
        );
    }

    std::uint64_t max_flow(int source, int sink) {
        std::uint64_t flow = 0U;
        while (breadth_first(source, sink)) {
            std::fill(next_.begin(), next_.end(), 0);
            while (true) {
                const std::uint64_t augmentation = depth_first(
                    source,
                    sink,
                    std::numeric_limits<std::uint64_t>::max()
                );
                if (augmentation == 0U) {
                    break;
                }
                if (flow
                    > std::numeric_limits<std::uint64_t>::max()
                        - augmentation) {
                    throw std::runtime_error("maximum flow exceeds uint64");
                }
                flow += augmentation;
            }
        }
        return flow;
    }

    std::vector<bool> residual_reachable(int source) const {
        std::vector<bool> reached(graph_.size(), false);
        std::queue<int> queue;
        reached[static_cast<std::size_t>(source)] = true;
        queue.push(source);
        while (!queue.empty()) {
            const int vertex = queue.front();
            queue.pop();
            for (const FlowEdge& edge :
                 graph_[static_cast<std::size_t>(vertex)]) {
                if (edge.capacity > 0U
                    && !reached[static_cast<std::size_t>(edge.target)]) {
                    reached[static_cast<std::size_t>(edge.target)] = true;
                    queue.push(edge.target);
                }
            }
        }
        return reached;
    }

    std::uint64_t flow_sent(int vertex, int edge_index) const {
        const FlowEdge& edge =
            graph_[static_cast<std::size_t>(vertex)][
                static_cast<std::size_t>(edge_index)
            ];
        return graph_[static_cast<std::size_t>(edge.target)][
            static_cast<std::size_t>(edge.reverse)
        ].capacity;
    }

private:
    bool breadth_first(int source, int sink) {
        std::fill(level_.begin(), level_.end(), -1);
        std::queue<int> queue;
        level_[static_cast<std::size_t>(source)] = 0;
        queue.push(source);
        while (!queue.empty()) {
            const int vertex = queue.front();
            queue.pop();
            for (const FlowEdge& edge :
                 graph_[static_cast<std::size_t>(vertex)]) {
                if (edge.capacity > 0U
                    && level_[static_cast<std::size_t>(edge.target)] < 0) {
                    level_[static_cast<std::size_t>(edge.target)] =
                        level_[static_cast<std::size_t>(vertex)] + 1;
                    queue.push(edge.target);
                }
            }
        }
        return level_[static_cast<std::size_t>(sink)] >= 0;
    }

    std::uint64_t depth_first(
        int vertex,
        int sink,
        std::uint64_t available
    ) {
        if (vertex == sink) {
            return available;
        }
        int& edge_index = next_[static_cast<std::size_t>(vertex)];
        while (
            edge_index
            < static_cast<int>(
                graph_[static_cast<std::size_t>(vertex)].size()
            )
        ) {
            FlowEdge& edge =
                graph_[static_cast<std::size_t>(vertex)][
                    static_cast<std::size_t>(edge_index)
                ];
            if (edge.capacity > 0U
                && level_[static_cast<std::size_t>(edge.target)]
                    == level_[static_cast<std::size_t>(vertex)] + 1) {
                const std::uint64_t pushed = depth_first(
                    edge.target,
                    sink,
                    std::min(available, edge.capacity)
                );
                if (pushed > 0U) {
                    edge.capacity -= pushed;
                    FlowEdge& reverse =
                        graph_[static_cast<std::size_t>(edge.target)][
                            static_cast<std::size_t>(edge.reverse)
                        ];
                    reverse.capacity += pushed;
                    return pushed;
                }
            }
            ++edge_index;
        }
        return 0U;
    }

    std::vector<std::vector<FlowEdge>> graph_;
    std::vector<int> level_;
    std::vector<int> next_;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5 && argc != 6) {
            throw std::runtime_error(
                "usage: analyze_su2_bpc_suffix_reflection "
                "LEVEL LABEL PREFIX TRUNCATION "
                "[--crossed|--aligned|--matching|--matching-local|"
                "--matching-core|--matching-capacity|"
                "--matching-capacity-local|--matching-capacity-down|"
                "--matching-capacity-down-core|"
                "--matching-capacity-down-segment|"
                "--matching-capacity-down-word|"
                "--matching-capacity-down-concat|"
                "--matching-capacity-down-rotation|"
                "--matching-capacity-down-shell|"
                "--matching-capacity-down-basic|"
                "--matching-capacity-adjacent]"
            );
        }
        const bool crossed_splice =
            argc == 6 && std::string(argv[5]) == "--crossed";
        const bool aligned_splice =
            argc == 6 && std::string(argv[5]) == "--aligned";
        const bool capacity_down_core_matching =
            argc == 6
            && std::string(argv[5])
                == "--matching-capacity-down-core";
        const bool capacity_down_segment_matching =
            argc == 6
            && std::string(argv[5])
                == "--matching-capacity-down-segment";
        const bool capacity_down_word_matching =
            argc == 6
            && std::string(argv[5])
                == "--matching-capacity-down-word";
        const bool capacity_down_concat_matching =
            argc == 6
            && std::string(argv[5])
                == "--matching-capacity-down-concat";
        const bool capacity_down_rotation_matching =
            argc == 6
            && std::string(argv[5])
                == "--matching-capacity-down-rotation";
        const bool capacity_down_shell_matching =
            argc == 6
            && std::string(argv[5])
                == "--matching-capacity-down-shell";
        const bool capacity_down_basic_matching =
            argc == 6
            && std::string(argv[5])
                == "--matching-capacity-down-basic";
        const bool minimal_matching =
            (
                argc == 6 && std::string(argv[5]) == "--matching-core"
            ) || capacity_down_core_matching
                || capacity_down_basic_matching;
        const bool capacity_local_matching =
            argc == 6
            && std::string(argv[5]) == "--matching-capacity-local";
        const bool capacity_down_matching =
            (
                argc == 6
                && std::string(argv[5]) == "--matching-capacity-down"
            ) || capacity_down_core_matching
                || capacity_down_segment_matching
                || capacity_down_word_matching
                || capacity_down_concat_matching
                || capacity_down_rotation_matching
                || capacity_down_shell_matching
                || capacity_down_basic_matching;
        const bool capacity_adjacent_matching =
            argc == 6
            && std::string(argv[5]) == "--matching-capacity-adjacent";
        const bool local_matching =
            (
                argc == 6 && std::string(argv[5]) == "--matching-local"
            ) || capacity_local_matching;
        const bool capacity_matching =
            (
                argc == 6 && std::string(argv[5]) == "--matching-capacity"
            ) || capacity_local_matching || capacity_down_matching
                || capacity_adjacent_matching;
        const bool cumulative_matching =
            (
                argc == 6 && std::string(argv[5]) == "--matching"
            ) || minimal_matching || local_matching || capacity_matching;
        if (argc == 6 && !crossed_splice && !aligned_splice
            && !cumulative_matching) {
            throw std::runtime_error(
                "the optional argument is --crossed, --aligned, "
                "--matching, --matching-local, --matching-core, or "
                "--matching-capacity"
                "[-local|-down|-down-core|-down-segment|-down-word|"
                "-down-concat|-down-rotation|-down-shell|-down-basic|"
                "-adjacent]"
            );
        }
        const std::string tag = capacity_matching
            ? (
                capacity_local_matching
                    ? "SU2_BPC_LOCAL_CAPACITY_MATCHING"
                    : (
                        capacity_down_matching
                            ? (
                                capacity_down_core_matching
                                    ? "SU2_BPC_DOWN_CORE_CAPACITY_MATCHING"
                                    : (
                                        capacity_down_segment_matching
                                            ? "SU2_BPC_DOWN_SEGMENT_CAPACITY_MATCHING"
                                            : (
                                                capacity_down_word_matching
                                                    ? "SU2_BPC_DOWN_WORD_CAPACITY_MATCHING"
                                                    : (
                                                        capacity_down_concat_matching
                                                            ? "SU2_BPC_DOWN_CONCAT_CAPACITY_MATCHING"
                                                            : (
                                                                capacity_down_rotation_matching
                                                                    ? "SU2_BPC_DOWN_ROTATION_CAPACITY_MATCHING"
                                                                    : (
                                                                        capacity_down_shell_matching
                                                                            ? "SU2_BPC_DOWN_SHELL_CAPACITY_MATCHING"
                                                                            : (
                                                                                capacity_down_basic_matching
                                                                                    ? "SU2_BPC_DOWN_BASIC_CAPACITY_MATCHING"
                                                                                    : "SU2_BPC_DOWN_CAPACITY_MATCHING"
                                                                            )
                                                                    )
                                                            )
                                                    )
                                            )
                                    )
                            )
                            : (
                                capacity_adjacent_matching
                                    ? "SU2_BPC_ADJACENT_CAPACITY_MATCHING"
                                    : "SU2_BPC_CAPACITY_MATCHING"
                            )
                    )
            )
            : (
                minimal_matching
                    ? "SU2_BPC_CORE_MATCHING"
                    : (
                        local_matching
                            ? "SU2_BPC_LOCAL_MATCHING"
                            : (
                                cumulative_matching
                                    ? "SU2_BPC_CUMULATIVE_MATCHING"
                                    : (
                                        aligned_splice
                                            ? "SU2_BPC_ALIGNED_SPLICE"
                                            : (
                                                crossed_splice
                                                    ? "SU2_BPC_CROSSED_SPLICE"
                                                    : "SU2_BPC_SUFFIX_REFLECTION"
                                            )
                                    )
                            )
                    )
            );
        const int level = parse_positive(argv[1], "level");
        const int label = parse_positive(argv[2], "label");
        const int prefix = parse_positive(argv[3], "prefix");
        const int truncation = parse_positive(argv[4], "truncation");
        if ((level & 1) != 0 || (label & 1) != 0
            || 2 * label >= level || prefix < 4
            || truncation < 2 || truncation >= prefix) {
            throw std::runtime_error("parameters are outside the BPC range");
        }

        const int length = 2 * prefix + 2;
        const int minimum_negative_y = 1;
        const int maximum_negative_y = 2 * truncation + 1;
        const int minimum_positive_y = 0;
        const int maximum_positive_y = 2 * truncation;

        std::vector<std::vector<int>> adjacency(
            static_cast<std::size_t>(level + 1)
        );
        for (int source = 0; source <= level; ++source) {
            adjacency[static_cast<std::size_t>(source)] =
                neighbors(level, label, source);
        }

        std::uint64_t negative_paths = 0;
        std::uint64_t covered_paths = 0;
        std::uint64_t uncovered_paths = 0;
        std::uint64_t invalid_images = 0;
        std::uint64_t collisions = 0;
        std::map<std::vector<int>, Preimage> images;
        std::map<std::vector<int>, int> matching_positive_ids;
        std::map<std::vector<int>, std::vector<int>>
            matching_candidate_cache;
        std::map<std::vector<int>, int> capacity_left_ids;
        std::vector<std::uint64_t> capacity_left_demand;
        std::vector<std::uint64_t> capacity_right_supply;
        std::vector<std::vector<int>> matching_edges;
        std::vector<int> matching_negative_y;
        std::vector<int> matching_positive_y;
        std::vector<Preimage> matching_negative_objects;
        std::vector<std::vector<int>> matching_positive_keys;
        std::uint64_t matching_edge_count = 0;
        std::uint64_t matching_invalid_path_pairs = 0;
        bool printed_uncovered = false;
        bool printed_collision = false;

        for (int model = 0; model <= 1; ++model) {
            std::vector<State> states;
            std::vector<char> moves;
            states.push_back({0, model == 0 ? 0 : level});

            std::function<void(int, int)> visit =
                [&](int step, int y_count) {
                    if (y_count > maximum_negative_y) {
                        return;
                    }
                    if (step == length) {
                        const State expected{
                            level, model == 0 ? 0 : level
                        };
                        if (states.back().first != expected.first
                            || states.back().second != expected.second
                            || moves.back() != 'Y'
                            || y_count < minimum_negative_y
                            || (y_count & 1) == 0) {
                            return;
                        }
                        ++negative_paths;

                        std::vector<int> first_path;
                        std::vector<int> second_path;
                        first_path.push_back(states.front().first);
                        second_path.push_back(
                            model == 0
                                ? level - states.front().second
                                : states.front().second
                        );
                        for (std::size_t index = 0U;
                             index < moves.size();
                             ++index) {
                            if (moves[index] == 'X') {
                                first_path.push_back(
                                    states[index + 1U].first
                                );
                            } else {
                                second_path.push_back(
                                    model == 0
                                        ? level
                                            - states[index + 1U].second
                                        : states[index + 1U].second
                                );
                            }
                        }

                        if (cumulative_matching) {
                            std::vector<int> source_pair_key;
                            source_pair_key.push_back(model);
                            source_pair_key.push_back(
                                static_cast<int>(first_path.size())
                            );
                            source_pair_key.insert(
                                source_pair_key.end(),
                                first_path.begin(),
                                first_path.end()
                            );
                            source_pair_key.push_back(
                                static_cast<int>(second_path.size())
                            );
                            source_pair_key.insert(
                                source_pair_key.end(),
                                second_path.begin(),
                                second_path.end()
                            );
                            if (capacity_matching) {
                                const auto [position, inserted] =
                                    capacity_left_ids.emplace(
                                        source_pair_key,
                                        static_cast<int>(
                                            capacity_left_ids.size()
                                        )
                                    );
                                if (!inserted) {
                                    ++capacity_left_demand[
                                        static_cast<std::size_t>(
                                            position->second
                                        )
                                    ];
                                    return;
                                }
                                capacity_left_demand.push_back(1U);
                                matching_negative_y.push_back(y_count);
                            } else {
                                matching_negative_objects.push_back(
                                    Preimage{states, moves, -1, -1}
                                );
                            }
                            const auto cached =
                                matching_candidate_cache.find(
                                    source_pair_key
                                );
                            if (!capacity_matching
                                && cached
                                != matching_candidate_cache.end()) {
                                matching_edges.push_back(cached->second);
                                matching_negative_y.push_back(y_count);
                                matching_edge_count +=
                                    cached->second.size();
                                return;
                            }
                            std::set<int> candidate_ids;
                            std::set<std::vector<int>> candidate_path_pairs;
                            const auto add_path_pair =
                                [&](const std::vector<int>& image_first,
                                    const std::vector<int>& image_second) {
                                    if (image_first.empty()
                                        || image_second.empty()
                                        || image_first.front() != 0
                                        || image_first.back() != level
                                        || image_second.front() != level
                                        || image_second.back() != level) {
                                        ++matching_invalid_path_pairs;
                                        return;
                                    }
                                    for (std::size_t edge = 0;
                                         edge + 1U < image_first.size();
                                         ++edge) {
                                        if (!fuses(
                                                level,
                                                label,
                                                image_first[edge],
                                                image_first[edge + 1U]
                                            )) {
                                            ++matching_invalid_path_pairs;
                                            return;
                                        }
                                    }
                                    for (std::size_t edge = 0;
                                         edge + 1U < image_second.size();
                                         ++edge) {
                                        if (!fuses(
                                                level,
                                                label,
                                                image_second[edge],
                                                image_second[edge + 1U]
                                            )) {
                                            ++matching_invalid_path_pairs;
                                            return;
                                        }
                                    }
                                    const int image_first_length =
                                        static_cast<int>(
                                            image_first.size()
                                        ) - 1;
                                    const int image_second_length =
                                        static_cast<int>(
                                            image_second.size()
                                        ) - 1;
                                    if (image_first_length
                                            + image_second_length
                                            != length
                                        || (image_second_length & 1) != 0
                                        || image_second_length
                                            < minimum_positive_y
                                        || image_second_length
                                            > maximum_positive_y
                                        || (
                                            capacity_down_matching
                                            && image_second_length
                                                > y_count - 1
                                        ) || (
                                            capacity_adjacent_matching
                                            && image_second_length
                                                != y_count - 1
                                        )) {
                                        return;
                                    }
                                    std::vector<int> path_pair_key;
                                    path_pair_key.push_back(
                                        static_cast<int>(
                                            image_first.size()
                                        )
                                    );
                                    path_pair_key.insert(
                                        path_pair_key.end(),
                                        image_first.begin(),
                                        image_first.end()
                                    );
                                    path_pair_key.push_back(
                                        static_cast<int>(
                                            image_second.size()
                                        )
                                    );
                                    path_pair_key.insert(
                                        path_pair_key.end(),
                                        image_second.begin(),
                                        image_second.end()
                                    );
                                    if (!candidate_path_pairs
                                            .insert(path_pair_key)
                                            .second) {
                                        return;
                                    }
                                    if (capacity_matching) {
                                        path_pair_key.insert(
                                            path_pair_key.begin(),
                                            model
                                        );
                                        const auto [position, inserted] =
                                            matching_positive_ids.emplace(
                                                path_pair_key,
                                                static_cast<int>(
                                                    matching_positive_ids
                                                        .size()
                                                )
                                            );
                                        if (inserted) {
                                            capacity_right_supply.push_back(
                                                binomial(
                                                    length - 1,
                                                    image_second_length
                                                )
                                            );
                                            matching_positive_y.push_back(
                                                image_second_length
                                            );
                                        }
                                        candidate_ids.insert(
                                            position->second
                                        );
                                        return;
                                    }
                                    std::vector<char> image_moves;
                                    image_moves.reserve(
                                        static_cast<std::size_t>(length)
                                    );
                                    std::function<void(int, int, int)>
                                        enumerate_shuffles =
                                            [&](int shuffle_step,
                                                int used_first,
                                                int used_second) {
                                                if (
                                                    shuffle_step
                                                        == length - 1
                                                ) {
                                                    if (
                                                        used_first
                                                            != image_first_length
                                                                - 1
                                                        || used_second
                                                            != image_second_length
                                                    ) {
                                                        return;
                                                    }
                                                    image_moves.push_back('X');
                                                    std::vector<State>
                                                        image_states;
                                                    image_states.reserve(
                                                        image_moves.size()
                                                            + 1U
                                                    );
                                                    int first_index = 0;
                                                    int second_index = 0;
                                                    image_states.push_back({
                                                        image_first.front(),
                                                        model == 0
                                                            ? level
                                                                - image_second
                                                                    .front()
                                                            : image_second
                                                                .front()
                                                    });
                                                    for (const char move :
                                                         image_moves) {
                                                        if (move == 'X') {
                                                            ++first_index;
                                                        } else {
                                                            ++second_index;
                                                        }
                                                        image_states.push_back({
                                                            image_first[
                                                                static_cast<
                                                                    std::size_t
                                                                >(first_index)
                                                            ],
                                                            model == 0
                                                                ? level
                                                                    - image_second[
                                                                        static_cast<
                                                                            std::size_t
                                                                        >(
                                                                            second_index
                                                                        )
                                                                    ]
                                                                : image_second[
                                                                    static_cast<
                                                                        std::size_t
                                                                    >(
                                                                        second_index
                                                                    )
                                                                ]
                                                        });
                                                    }
                                                    const std::vector<int>
                                                        key = flatten(
                                                            model,
                                                            image_states,
                                                            image_moves
                                                        );
                                                    const auto [
                                                        position,
                                                        inserted
                                                    ] =
                                                        matching_positive_ids
                                                            .emplace(
                                                                key,
                                                                static_cast<
                                                                    int
                                                                >(
                                                                    matching_positive_ids
                                                                        .size()
                                                                )
                                                            );
                                                    static_cast<void>(
                                                        inserted
                                                    );
                                                    if (inserted) {
                                                        matching_positive_y
                                                            .push_back(
                                                                image_second_length
                                                            );
                                                        matching_positive_keys
                                                            .push_back(key);
                                                    }
                                                    candidate_ids.insert(
                                                        position->second
                                                    );
                                                    image_moves.pop_back();
                                                    return;
                                                }
                                                if (
                                                    used_first
                                                        < image_first_length
                                                            - 1
                                                ) {
                                                    image_moves.push_back('X');
                                                    enumerate_shuffles(
                                                        shuffle_step + 1,
                                                        used_first + 1,
                                                        used_second
                                                    );
                                                    image_moves.pop_back();
                                                }
                                                if (
                                                    used_second
                                                        < image_second_length
                                                ) {
                                                    image_moves.push_back('Y');
                                                    enumerate_shuffles(
                                                        shuffle_step + 1,
                                                        used_first,
                                                        used_second + 1
                                                    );
                                                    image_moves.pop_back();
                                                }
                                            };
                                    enumerate_shuffles(0, 0, 0);
                                };
                            const auto add_splice_pair =
                                [&](int cut_first, int cut_second) {
                                    std::vector<int> image_first(
                                        first_path.begin(),
                                        first_path.begin() + cut_first + 1
                                    );
                                    image_first.insert(
                                        image_first.end(),
                                        second_path.begin()
                                            + cut_second + 1,
                                        second_path.end()
                                    );
                                    std::vector<int> image_second(
                                        second_path.begin(),
                                        second_path.begin()
                                            + cut_second + 1
                                    );
                                    image_second.insert(
                                        image_second.end(),
                                        first_path.begin() + cut_first + 1,
                                        first_path.end()
                                    );
                                    add_path_pair(
                                        image_first,
                                        image_second
                                    );
                                };
                            const auto add_segment_pair =
                                [&](int first_begin,
                                    int first_end,
                                    int second_begin,
                                    int second_end) {
                                    std::vector<int> image_first(
                                        first_path.begin(),
                                        first_path.begin()
                                            + first_begin + 1
                                    );
                                    image_first.insert(
                                        image_first.end(),
                                        second_path.begin()
                                            + second_begin + 1,
                                        second_path.begin()
                                            + second_end + 1
                                    );
                                    image_first.insert(
                                        image_first.end(),
                                        first_path.begin()
                                            + first_end + 1,
                                        first_path.end()
                                    );
                                    std::vector<int> image_second(
                                        second_path.begin(),
                                        second_path.begin()
                                            + second_begin + 1
                                    );
                                    image_second.insert(
                                        image_second.end(),
                                        first_path.begin()
                                            + first_begin + 1,
                                        first_path.begin()
                                            + first_end + 1
                                    );
                                    image_second.insert(
                                        image_second.end(),
                                        second_path.begin()
                                            + second_end + 1,
                                        second_path.end()
                                    );
                                    add_path_pair(
                                        image_first,
                                        image_second
                                    );
                                };

                            const int first_length =
                                static_cast<int>(first_path.size()) - 1;
                            const int second_length =
                                static_cast<int>(second_path.size()) - 1;
                            const int first_offset =
                                first_length - second_length + 1;
                            const int first_second_index =
                                std::max(0, -first_offset);
                            for (int second_index = first_second_index;
                                 second_index <= second_length - 2;
                                 ++second_index) {
                                const int first_index =
                                    first_offset + second_index;
                                if (fuses(
                                        level,
                                        label,
                                        first_path[
                                            static_cast<std::size_t>(
                                                first_index
                                            )
                                        ],
                                        second_path[
                                            static_cast<std::size_t>(
                                                second_index + 1
                                            )
                                        ]
                                    )
                                    && fuses(
                                        level,
                                        label,
                                        second_path[
                                            static_cast<std::size_t>(
                                                second_index
                                            )
                                        ],
                                        first_path[
                                            static_cast<std::size_t>(
                                                first_index + 1
                                            )
                                        ]
                                    )) {
                                    add_splice_pair(
                                        first_index,
                                        second_index
                                    );
                                }
                            }

                            for (int first_index = 0;
                                 first_index < first_length;
                                 ++first_index) {
                                for (int second_index = 0;
                                     second_index < second_length;
                                     ++second_index) {
                                    const int image_second_length =
                                        second_index
                                        + first_length
                                        - first_index;
                                    if (
                                        (image_second_length & 1) != 0
                                        || image_second_length
                                            < minimum_positive_y
                                        || image_second_length
                                            > maximum_positive_y
                                    ) {
                                        continue;
                                    }
                                    if (fuses(
                                            level,
                                            label,
                                            first_path[
                                                static_cast<std::size_t>(
                                                    first_index
                                                )
                                            ],
                                            second_path[
                                                static_cast<std::size_t>(
                                                    second_index + 1
                                                )
                                            ]
                                        )
                                        && fuses(
                                            level,
                                            label,
                                            second_path[
                                                static_cast<std::size_t>(
                                                    second_index
                                                )
                                            ],
                                            first_path[
                                                static_cast<std::size_t>(
                                                    first_index + 1
                                                )
                                            ]
                                        )) {
                                        add_splice_pair(
                                            first_index,
                                            second_index
                                        );
                                    }
                                }
                            }
                            for (int first_index = 0;
                                 first_index <= first_length;
                                 ++first_index) {
                                for (int second_index = 0;
                                     second_index <= second_length;
                                     ++second_index) {
                                    const int image_second_length =
                                        second_index
                                        + first_length
                                        - first_index;
                                    if (
                                        first_path[
                                            static_cast<std::size_t>(
                                                first_index
                                            )
                                        ] == second_path[
                                            static_cast<std::size_t>(
                                                second_index
                                            )
                                        ]
                                        && (image_second_length & 1) == 0
                                        && minimum_positive_y
                                            <= image_second_length
                                        && image_second_length
                                            <= maximum_positive_y
                                    ) {
                                        add_splice_pair(
                                            first_index,
                                            second_index
                                        );
                                    }
                                }
                            }
                            if (!local_matching) {
                                for (int second_index = 0;
                                     second_index < second_length;
                                     ++second_index) {
                                    if (second_path[
                                            static_cast<std::size_t>(
                                                second_index
                                            )
                                        ] != second_path[
                                            static_cast<std::size_t>(
                                                second_index + 1
                                            )
                                        ]) {
                                        continue;
                                    }
                                    for (int first_index = 0;
                                         first_index <= first_length;
                                         ++first_index) {
                                        const int vertex = first_path[
                                            static_cast<std::size_t>(
                                                first_index
                                            )
                                        ];
                                        if (!fuses(
                                                level,
                                                label,
                                                vertex,
                                                vertex
                                            )) {
                                            continue;
                                        }
                                        std::vector<int> image_first =
                                            first_path;
                                        image_first.insert(
                                            image_first.begin()
                                                + first_index + 1,
                                            vertex
                                        );
                                        std::vector<int> image_second =
                                            second_path;
                                        image_second.erase(
                                            image_second.begin()
                                                + second_index + 1
                                        );
                                        add_path_pair(
                                            image_first,
                                            image_second
                                        );
                                    }
                                }
                            }
                            if (capacity_down_basic_matching) {
                                const auto add_basic_shells =
                                    [&](const std::vector<int>& loop) {
                                        if (
                                            loop.size() < 3U
                                            || first_path.size() < 2U
                                        ) {
                                            return;
                                        }
                                        std::vector<int> upper_first(
                                            first_path.begin(),
                                            first_path.end() - 1
                                        );
                                        upper_first.insert(
                                            upper_first.end(),
                                            loop.begin() + 2,
                                            loop.end() - 1
                                        );
                                        upper_first.push_back(level);
                                        add_path_pair(
                                            upper_first,
                                            std::vector<int>{
                                                level,
                                                loop[1],
                                                level
                                            }
                                        );

                                        std::vector<int> reflected_loop;
                                        reflected_loop.reserve(loop.size());
                                        for (const int vertex : loop) {
                                            reflected_loop.push_back(
                                                level - vertex
                                            );
                                        }
                                        std::vector<int> lower_first{
                                            first_path[0],
                                            first_path[1]
                                        };
                                        lower_first.insert(
                                            lower_first.end(),
                                            reflected_loop.begin() + 2,
                                            reflected_loop.end() - 1
                                        );
                                        lower_first.insert(
                                            lower_first.end(),
                                            first_path.begin() + 2,
                                            first_path.end()
                                        );
                                        add_path_pair(
                                            lower_first,
                                            std::vector<int>{
                                                level,
                                                level - reflected_loop[1],
                                                level
                                            }
                                        );
                                    };
                                add_basic_shells(second_path);
                                add_basic_shells(
                                    std::vector<int>(
                                        second_path.rbegin(),
                                        second_path.rend()
                                    )
                                );
                            }
                            if (!minimal_matching) {
                                for (int first_begin = 0;
                                   first_begin < first_length;
                                   ++first_begin) {
                                for (int first_end = first_begin + 1;
                                     first_end < first_length;
                                     ++first_end) {
                                    for (int second_begin = 0;
                                         second_begin < second_length;
                                         ++second_begin) {
                                        for (
                                            int second_end =
                                                second_begin + 1;
                                            second_end < second_length;
                                            ++second_end
                                        ) {
                                            const int image_second_length =
                                                second_length
                                                + first_end
                                                - first_begin
                                                - second_end
                                                + second_begin;
                                            if (
                                                (image_second_length & 1)
                                                    != 0
                                                || image_second_length
                                                    < minimum_positive_y
                                                || image_second_length
                                                    > maximum_positive_y
                                            ) {
                                                continue;
                                            }
                                            if (
                                                fuses(
                                                    level,
                                                    label,
                                                    first_path[
                                                        static_cast<
                                                            std::size_t
                                                        >(first_begin)
                                                    ],
                                                    second_path[
                                                        static_cast<
                                                            std::size_t
                                                        >(second_begin + 1)
                                                    ]
                                                )
                                                && fuses(
                                                    level,
                                                    label,
                                                    second_path[
                                                        static_cast<
                                                            std::size_t
                                                        >(second_begin)
                                                    ],
                                                    first_path[
                                                        static_cast<
                                                            std::size_t
                                                        >(first_begin + 1)
                                                    ]
                                                )
                                                && fuses(
                                                    level,
                                                    label,
                                                    second_path[
                                                        static_cast<
                                                            std::size_t
                                                        >(second_end)
                                                    ],
                                                    first_path[
                                                        static_cast<
                                                            std::size_t
                                                        >(first_end + 1)
                                                    ]
                                                )
                                                && fuses(
                                                    level,
                                                    label,
                                                    first_path[
                                                        static_cast<
                                                            std::size_t
                                                        >(first_end)
                                                    ],
                                                    second_path[
                                                        static_cast<
                                                            std::size_t
                                                        >(second_end + 1)
                                                    ]
                                                )
                                            ) {
                                                add_segment_pair(
                                                    first_begin,
                                                    first_end,
                                                    second_begin,
                                                    second_end
                                                );
                                            }
                                        }
                                    }
                                }
                                }

                            if (!capacity_down_segment_matching) {
                            int prefix_x = 0;
                            int prefix_y = 0;
                            for (int index = 0; index < length; ++index) {
                                const int suffix_x =
                                    first_length - prefix_x;
                                const int suffix_y =
                                    second_length - prefix_y;
                                const int new_y = prefix_y + suffix_x;
                                if (suffix_y > suffix_x
                                    && (new_y & 1) == 0
                                    && minimum_positive_y <= new_y
                                    && new_y <= maximum_positive_y
                                    && suffix_y > 0) {
                                    const bool first_connector =
                                        fuses(
                                            level,
                                            label,
                                            first_path[
                                                static_cast<std::size_t>(
                                                    prefix_x
                                                )
                                            ],
                                            second_path[
                                                static_cast<std::size_t>(
                                                    prefix_y + 1
                                                )
                                            ]
                                        );
                                    const bool second_connector =
                                        suffix_x == 0
                                            ? first_path[
                                                static_cast<std::size_t>(
                                                    prefix_x
                                                )
                                            ] == second_path[
                                                static_cast<std::size_t>(
                                                    prefix_y
                                                )
                                            ]
                                            : fuses(
                                                level,
                                                label,
                                                second_path[
                                                    static_cast<std::size_t>(
                                                        prefix_y
                                                    )
                                                ],
                                                first_path[
                                                    static_cast<std::size_t>(
                                                        prefix_x + 1
                                                    )
                                                ]
                                            );
                                    if (first_connector
                                        && second_connector) {
                                        add_splice_pair(
                                            prefix_x,
                                            prefix_y
                                        );
                                    }
                                }
                                if (moves[
                                        static_cast<std::size_t>(index)
                                    ] == 'X') {
                                    ++prefix_x;
                                } else {
                                    ++prefix_y;
                                }
                            }
                            if (!capacity_down_word_matching) {
                            std::vector<int> concatenated = first_path;
                            concatenated.insert(
                                concatenated.end(),
                                second_path.begin() + 1,
                                second_path.end()
                            );
                            add_path_pair(
                                concatenated,
                                std::vector<int>{level}
                            );
                            std::vector<int> reverse_concatenated =
                                first_path;
                            reverse_concatenated.insert(
                                reverse_concatenated.end(),
                                second_path.rbegin() + 1,
                                second_path.rend()
                            );
                            add_path_pair(
                                reverse_concatenated,
                                std::vector<int>{level}
                            );
                            const auto add_loop_extractions =
                                [&](const std::vector<int>& full_path) {
                                    for (
                                        int loop_start = 0;
                                        loop_start <= length;
                                        ++loop_start
                                    ) {
                                        if (full_path[
                                                static_cast<std::size_t>(
                                                    loop_start
                                                )
                                            ] != level) {
                                            continue;
                                        }
                                        for (
                                            int loop_end = loop_start;
                                            loop_end <= length;
                                            ++loop_end
                                        ) {
                                            const int loop_length =
                                                loop_end - loop_start;
                                            if (
                                                (loop_length & 1) != 0
                                                || loop_length
                                                    > maximum_positive_y
                                                || full_path[
                                                    static_cast<std::size_t>(
                                                        loop_end
                                                    )
                                                ] != level
                                            ) {
                                                continue;
                                            }
                                            std::vector<int> image_first(
                                                full_path.begin(),
                                                full_path.begin()
                                                    + loop_start + 1
                                            );
                                            image_first.insert(
                                                image_first.end(),
                                                full_path.begin()
                                                    + loop_end + 1,
                                                full_path.end()
                                            );
                                            std::vector<int> image_second(
                                                full_path.begin()
                                                    + loop_start,
                                                full_path.begin()
                                                    + loop_end + 1
                                            );
                                            add_path_pair(
                                                image_first,
                                                image_second
                                            );
                                        }
                                    }
                                };
                            const auto
                                add_reflected_lower_loop_extractions =
                                    [&](const std::vector<int>& full_path) {
                                        for (
                                            int loop_start = 0;
                                            loop_start <= length;
                                            ++loop_start
                                        ) {
                                            if (full_path[
                                                    static_cast<std::size_t>(
                                                        loop_start
                                                    )
                                                ] != 0) {
                                                continue;
                                            }
                                            for (
                                                int loop_end = loop_start;
                                                loop_end <= length;
                                                ++loop_end
                                            ) {
                                                const int loop_length =
                                                    loop_end - loop_start;
                                                if (
                                                    (loop_length & 1) != 0
                                                    || loop_length
                                                        > maximum_positive_y
                                                    || full_path[
                                                        static_cast<
                                                            std::size_t
                                                        >(loop_end)
                                                    ] != 0
                                                ) {
                                                    continue;
                                                }
                                                std::vector<int> image_first(
                                                    full_path.begin(),
                                                    full_path.begin()
                                                        + loop_start + 1
                                                );
                                                image_first.insert(
                                                    image_first.end(),
                                                    full_path.begin()
                                                        + loop_end + 1,
                                                    full_path.end()
                                                );
                                                std::vector<int>
                                                    image_second;
                                                image_second.reserve(
                                                    static_cast<std::size_t>(
                                                        loop_length + 1
                                                    )
                                                );
                                                for (
                                                    int index = loop_start;
                                                    index <= loop_end;
                                                    ++index
                                                ) {
                                                    image_second.push_back(
                                                        level - full_path[
                                                            static_cast<
                                                                std::size_t
                                                            >(index)
                                                        ]
                                                    );
                                                }
                                                add_path_pair(
                                                    image_first,
                                                    image_second
                                                );
                                            }
                                        }
                                    };
                            add_loop_extractions(concatenated);
                            add_loop_extractions(reverse_concatenated);
                            if (!capacity_down_concat_matching) {
                            for (int rotation_start = 0;
                                 rotation_start < second_length;
                                 ++rotation_start) {
                                if (second_path[
                                        static_cast<std::size_t>(
                                            rotation_start
                                        )
                                    ] != level) {
                                    continue;
                                }
                                std::vector<int> rotated_loop(
                                    second_path.begin() + rotation_start,
                                    second_path.end()
                                );
                                rotated_loop.insert(
                                    rotated_loop.end(),
                                    second_path.begin() + 1,
                                    second_path.begin()
                                        + rotation_start + 1
                                );
                                for (int orientation = 0;
                                     orientation < 2;
                                     ++orientation) {
                                    const std::vector<int> oriented_loop =
                                        orientation == 0
                                            ? rotated_loop
                                            : std::vector<int>(
                                                rotated_loop.rbegin(),
                                                rotated_loop.rend()
                                            );
                                    std::vector<int> full_path =
                                        first_path;
                                    full_path.insert(
                                        full_path.end(),
                                        oriented_loop.begin() + 1,
                                        oriented_loop.end()
                                    );
                                    add_path_pair(
                                        full_path,
                                        std::vector<int>{level}
                                    );
                                    add_loop_extractions(full_path);

                                    if (
                                        !local_matching
                                        && !capacity_down_rotation_matching
                                    ) {
                                      if (oriented_loop.size() >= 3U
                                          && first_path.size() >= 2U) {
                                        std::vector<int> shell_first(
                                            first_path.begin(),
                                            first_path.end() - 1
                                        );
                                        shell_first.insert(
                                            shell_first.end(),
                                            oriented_loop.begin() + 2,
                                            oriented_loop.end() - 1
                                        );
                                        shell_first.push_back(level);
                                        const std::vector<int> shell_second{
                                            level,
                                            oriented_loop[1],
                                            level
                                        };
                                        add_path_pair(
                                            shell_first,
                                            shell_second
                                        );
                                      }

                                      std::vector<int> lower_loop;
                                      lower_loop.reserve(
                                          oriented_loop.size()
                                      );
                                      for (const int vertex :
                                           oriented_loop) {
                                          lower_loop.push_back(
                                              level - vertex
                                          );
                                      }
                                      if (capacity_down_shell_matching) {
                                          if (
                                              lower_loop.size() >= 3U
                                              && first_path.size() >= 2U
                                          ) {
                                              std::vector<int>
                                                  lower_shell_first{
                                                      first_path[0],
                                                      first_path[1]
                                                  };
                                              lower_shell_first.insert(
                                                  lower_shell_first.end(),
                                                  lower_loop.begin() + 2,
                                                  lower_loop.end() - 1
                                              );
                                              lower_shell_first.insert(
                                                  lower_shell_first.end(),
                                                  first_path.begin() + 2,
                                                  first_path.end()
                                              );
                                              const std::vector<int>
                                                  lower_shell_second{
                                                      level,
                                                      level - lower_loop[1],
                                                      level
                                                  };
                                              add_path_pair(
                                                  lower_shell_first,
                                                  lower_shell_second
                                              );
                                          }
                                          continue;
                                      }
                                      if (lower_loop.size() >= 3U
                                          && first_path.size() >= 2U) {
                                          std::vector<int> lower_shell_first{
                                              first_path[0],
                                              first_path[1]
                                          };
                                          lower_shell_first.insert(
                                              lower_shell_first.end(),
                                              lower_loop.begin() + 2,
                                              lower_loop.end() - 1
                                          );
                                          lower_shell_first.insert(
                                              lower_shell_first.end(),
                                              first_path.begin() + 2,
                                              first_path.end()
                                          );
                                          const std::vector<int>
                                              lower_shell_second{
                                                  level,
                                                  level - lower_loop[1],
                                                  level
                                              };
                                          add_path_pair(
                                              lower_shell_first,
                                              lower_shell_second
                                          );

                                          std::vector<int> interior_loop(
                                              lower_loop.begin() + 1,
                                              lower_loop.end() - 1
                                          );
                                          std::vector<int> upper_tail(
                                              first_path.begin() + 1,
                                              first_path.end()
                                          );
                                          const int interior_length =
                                              static_cast<int>(
                                                  interior_loop.size()
                                              ) - 1;
                                          const int tail_length =
                                              static_cast<int>(
                                                  upper_tail.size()
                                              ) - 1;
                                          const auto add_interior_recoupling =
                                              [&](int tail_cut,
                                                  int interior_cut) {
                                                  std::vector<int> new_loop(
                                                      upper_tail.begin(),
                                                      upper_tail.begin()
                                                          + tail_cut + 1
                                                  );
                                                  new_loop.insert(
                                                      new_loop.end(),
                                                      interior_loop.begin()
                                                          + interior_cut + 1,
                                                      interior_loop.end()
                                                  );
                                                  std::vector<int> new_tail(
                                                      interior_loop.begin(),
                                                      interior_loop.begin()
                                                          + interior_cut + 1
                                                  );
                                                  new_tail.insert(
                                                      new_tail.end(),
                                                      upper_tail.begin()
                                                          + tail_cut + 1,
                                                      upper_tail.end()
                                                  );
                                                  std::vector<int>
                                                      image_first{
                                                          first_path.front()
                                                      };
                                                  image_first.insert(
                                                      image_first.end(),
                                                      new_loop.begin(),
                                                      new_loop.end()
                                                  );
                                                  image_first.insert(
                                                      image_first.end(),
                                                      new_tail.begin() + 1,
                                                      new_tail.end()
                                                  );
                                                  add_path_pair(
                                                      image_first,
                                                      lower_shell_second
                                                  );
                                              };
                                          for (int tail_cut = 0;
                                               tail_cut < tail_length;
                                               ++tail_cut) {
                                              for (
                                                  int interior_cut = 0;
                                                  interior_cut
                                                      < interior_length;
                                                  ++interior_cut
                                              ) {
                                                  if (
                                                      fuses(
                                                          level,
                                                          label,
                                                          upper_tail[
                                                              static_cast<
                                                                  std::size_t
                                                              >(tail_cut)
                                                          ],
                                                          interior_loop[
                                                              static_cast<
                                                                  std::size_t
                                                              >(
                                                                  interior_cut
                                                                  + 1
                                                              )
                                                          ]
                                                      )
                                                      && fuses(
                                                          level,
                                                          label,
                                                          interior_loop[
                                                              static_cast<
                                                                  std::size_t
                                                              >(interior_cut)
                                                          ],
                                                          upper_tail[
                                                              static_cast<
                                                                  std::size_t
                                                              >(tail_cut + 1)
                                                          ]
                                                      )
                                                  ) {
                                                      add_interior_recoupling(
                                                          tail_cut,
                                                          interior_cut
                                                      );
                                                  }
                                              }
                                          }
                                          for (int tail_cut = 0;
                                               tail_cut <= tail_length;
                                               ++tail_cut) {
                                              for (
                                                  int interior_cut = 0;
                                                  interior_cut
                                                      <= interior_length;
                                                  ++interior_cut
                                              ) {
                                                  if (
                                                      upper_tail[
                                                          static_cast<
                                                              std::size_t
                                                          >(tail_cut)
                                                      ] == interior_loop[
                                                          static_cast<
                                                              std::size_t
                                                          >(interior_cut)
                                                      ]
                                                  ) {
                                                      add_interior_recoupling(
                                                          tail_cut,
                                                          interior_cut
                                                      );
                                                  }
                                              }
                                          }
                                      }
                                      std::vector<int> relocated_path =
                                          lower_loop;
                                      relocated_path.insert(
                                          relocated_path.end(),
                                          first_path.begin() + 1,
                                          first_path.end()
                                      );
                                      add_path_pair(
                                          relocated_path,
                                          std::vector<int>{level}
                                      );
                                      add_loop_extractions(relocated_path);
                                      add_reflected_lower_loop_extractions(
                                          relocated_path
                                      );

                                      const int lower_length =
                                          static_cast<int>(
                                              lower_loop.size()
                                          ) - 1;
                                      for (int first_cut = 0;
                                           first_cut < first_length;
                                           ++first_cut) {
                                          for (int lower_cut = 0;
                                               lower_cut < lower_length;
                                               ++lower_cut) {
                                              if (
                                                  !fuses(
                                                      level,
                                                      label,
                                                      first_path[
                                                          static_cast<
                                                              std::size_t
                                                          >(first_cut)
                                                      ],
                                                      lower_loop[
                                                          static_cast<
                                                              std::size_t
                                                          >(lower_cut + 1)
                                                      ]
                                                  )
                                                  || !fuses(
                                                      level,
                                                      label,
                                                      lower_loop[
                                                          static_cast<
                                                              std::size_t
                                                          >(lower_cut)
                                                      ],
                                                      first_path[
                                                          static_cast<
                                                              std::size_t
                                                          >(first_cut + 1)
                                                      ]
                                                  )
                                              ) {
                                                  continue;
                                              }
                                              std::vector<int> new_lower(
                                                  first_path.begin(),
                                                  first_path.begin()
                                                      + first_cut + 1
                                              );
                                              new_lower.insert(
                                                  new_lower.end(),
                                                  lower_loop.begin()
                                                      + lower_cut + 1,
                                                  lower_loop.end()
                                              );
                                              std::vector<int> new_path(
                                                  lower_loop.begin(),
                                                  lower_loop.begin()
                                                      + lower_cut + 1
                                              );
                                              new_path.insert(
                                                  new_path.end(),
                                                  first_path.begin()
                                                      + first_cut + 1,
                                                  first_path.end()
                                              );
                                              std::vector<int> recoupled =
                                                  new_lower;
                                              recoupled.insert(
                                                  recoupled.end(),
                                                  new_path.begin() + 1,
                                                  new_path.end()
                                              );
                                              add_path_pair(
                                                  recoupled,
                                                  std::vector<int>{level}
                                              );
                                              add_loop_extractions(
                                                  recoupled
                                              );
                                              add_reflected_lower_loop_extractions(
                                                  recoupled
                                              );
                                          }
                                      }
                                      for (int first_cut = 0;
                                           first_cut <= first_length;
                                           ++first_cut) {
                                          for (int lower_cut = 0;
                                               lower_cut <= lower_length;
                                               ++lower_cut) {
                                              if (
                                                  first_path[
                                                      static_cast<std::size_t>(
                                                          first_cut
                                                      )
                                                  ] != lower_loop[
                                                      static_cast<std::size_t>(
                                                          lower_cut
                                                      )
                                                  ]
                                              ) {
                                                  continue;
                                              }
                                              std::vector<int> new_lower(
                                                  first_path.begin(),
                                                  first_path.begin()
                                                      + first_cut + 1
                                              );
                                              new_lower.insert(
                                                  new_lower.end(),
                                                  lower_loop.begin()
                                                      + lower_cut + 1,
                                                  lower_loop.end()
                                              );
                                              std::vector<int> new_path(
                                                  lower_loop.begin(),
                                                  lower_loop.begin()
                                                      + lower_cut + 1
                                              );
                                              new_path.insert(
                                                  new_path.end(),
                                                  first_path.begin()
                                                      + first_cut + 1,
                                                  first_path.end()
                                              );
                                              std::vector<int> recoupled =
                                                  new_lower;
                                              recoupled.insert(
                                                  recoupled.end(),
                                                  new_path.begin() + 1,
                                                  new_path.end()
                                              );
                                              add_path_pair(
                                                  recoupled,
                                                  std::vector<int>{level}
                                              );
                                              add_loop_extractions(
                                                  recoupled
                                              );
                                              add_reflected_lower_loop_extractions(
                                                  recoupled
                                              );
                                          }
                                      }
                                    }
                                }
                            }
                            }
                            }
                            }
                            }
                            matching_edges.emplace_back(
                                candidate_ids.begin(),
                                candidate_ids.end()
                            );
                            if (!capacity_matching) {
                                matching_candidate_cache.emplace(
                                    std::move(source_pair_key),
                                    matching_edges.back()
                                );
                                matching_negative_y.push_back(y_count);
                            }
                            matching_edge_count += candidate_ids.size();
                            if (candidate_ids.empty()
                                && !printed_uncovered) {
                                std::cout
                                    << tag
                                    << " first_isolated"
                                    << " model=" << model
                                    << " y_moves=" << y_count
                                    << " moves=";
                                for (const char move : moves) {
                                    std::cout << move;
                                }
                                std::cout << " states=";
                                for (const State& state : states) {
                                    std::cout
                                        << '(' << state.first
                                        << ',' << state.second << ')';
                                }
                                std::cout << '\n';
                                printed_uncovered = true;
                            }
                            return;
                        }

                        int cut = -1;
                        int cut_first = -1;
                        int cut_second = -1;
                        int suffix_x = 0;
                        int suffix_y = 0;
                        if (aligned_splice) {
                            const int first_length =
                                static_cast<int>(first_path.size()) - 1;
                            const int second_length =
                                static_cast<int>(second_path.size()) - 1;
                            const int first_offset =
                                first_length - second_length + 1;
                            const int first_second_index =
                                std::max(0, -first_offset);
                            for (int second_index = first_second_index;
                                 second_index <= second_length - 2;
                                 ++second_index) {
                                const int first_index =
                                    first_offset + second_index;
                                const int a = first_path[
                                    static_cast<std::size_t>(first_index)
                                ];
                                const int b = second_path[
                                    static_cast<std::size_t>(second_index)
                                ];
                                const int d = first_path[
                                    static_cast<std::size_t>(
                                        first_index + 1
                                    )
                                ];
                                const int c = second_path[
                                    static_cast<std::size_t>(
                                        second_index + 1
                                    )
                                ];
                                if (a <= b && c <= d) {
                                    if (fuses(level, label, a, c)
                                        && fuses(level, label, b, d)) {
                                        cut = length - 1;
                                        cut_first = first_index;
                                        cut_second = second_index;
                                    }
                                    break;
                                }
                            }
                        } else {
                            for (int index = length - 1;
                                 index >= 0;
                                 --index) {
                                if (moves[static_cast<std::size_t>(index)]
                                    == 'X') {
                                    ++suffix_x;
                                } else {
                                    ++suffix_y;
                                }
                                const State& state =
                                    states[
                                        static_cast<std::size_t>(index)
                                    ];
                                const int new_y =
                                    y_count - suffix_y + suffix_x;
                                const int first_index =
                                    static_cast<int>(first_path.size()) - 1
                                    - suffix_x;
                                const int second_index =
                                    static_cast<int>(second_path.size()) - 1
                                    - suffix_y;
                                bool splice_is_valid = model == 0
                                    ? state.first + state.second == level
                                    : state.first == state.second;
                                if (crossed_splice
                                    && suffix_y > 0
                                    && first_index >= 0
                                    && second_index >= 0) {
                                    const bool first_connector =
                                        fuses(
                                            level,
                                            label,
                                            first_path[
                                                static_cast<std::size_t>(
                                                    first_index
                                                )
                                            ],
                                            second_path[
                                                static_cast<std::size_t>(
                                                    second_index + 1
                                                )
                                            ]
                                        );
                                    const bool second_connector =
                                        suffix_x == 0
                                            ? first_path[
                                                static_cast<std::size_t>(
                                                    first_index
                                                )
                                            ] == second_path[
                                                static_cast<std::size_t>(
                                                    second_index
                                                )
                                            ]
                                            : fuses(
                                                level,
                                                label,
                                                second_path[
                                                    static_cast<std::size_t>(
                                                        second_index
                                                    )
                                                ],
                                                first_path[
                                                    static_cast<std::size_t>(
                                                        first_index + 1
                                                    )
                                                ]
                                            );
                                    splice_is_valid =
                                        first_connector && second_connector;
                                }
                                if (splice_is_valid
                                    && suffix_y > suffix_x
                                    && (new_y & 1) == 0
                                    && minimum_positive_y <= new_y
                                    && new_y <= maximum_positive_y) {
                                    cut = index;
                                    cut_first = first_index;
                                    cut_second = second_index;
                                    break;
                                }
                            }
                        }
                        if (cut < 0) {
                            ++uncovered_paths;
                            if (!printed_uncovered) {
                                std::cout
                                    << tag
                                    << " first_uncovered"
                                    << " model=" << model
                                    << " y_moves=" << y_count
                                    << " moves=";
                                for (const char move : moves) {
                                    std::cout << move;
                                }
                                std::cout << " states=";
                                for (const State& state : states) {
                                    std::cout << '(' << state.first
                                              << ',' << state.second
                                              << ')';
                                }
                                std::cout << '\n';
                                printed_uncovered = true;
                            }
                            return;
                        }

                        std::vector<int> image_first(
                            first_path.begin(),
                            first_path.begin() + cut_first + 1
                        );
                        image_first.insert(
                            image_first.end(),
                            second_path.begin() + cut_second + 1,
                            second_path.end()
                        );
                        std::vector<int> image_second(
                            second_path.begin(),
                            second_path.begin() + cut_second + 1
                        );
                        image_second.insert(
                            image_second.end(),
                            first_path.begin() + cut_first + 1,
                            first_path.end()
                        );

                        std::vector<char> image_moves = moves;
                        for (int index = cut; index < length; ++index) {
                            image_moves[static_cast<std::size_t>(index)] =
                                moves[static_cast<std::size_t>(index)]
                                    == 'X' ? 'Y' : 'X';
                        }
                        std::vector<State> image_states;
                        image_states.reserve(
                            image_moves.size() + 1U
                        );
                        std::size_t first_index = 0U;
                        std::size_t second_index = 0U;
                        image_states.push_back({
                            image_first.front(),
                            model == 0
                                ? level - image_second.front()
                                : image_second.front()
                        });
                        for (const char move : image_moves) {
                            if (move == 'X') {
                                ++first_index;
                            } else {
                                ++second_index;
                            }
                            image_states.push_back({
                                image_first[first_index],
                                model == 0
                                    ? level - image_second[second_index]
                                    : image_second[second_index]
                            });
                        }
                        if (!valid_path(
                                level,
                                label,
                                model,
                                image_states,
                                image_moves,
                                minimum_positive_y,
                                maximum_positive_y
                            )) {
                            ++invalid_images;
                            return;
                        }
                        ++covered_paths;
                        const std::vector<int> key = flatten(
                            model, image_states, image_moves
                        );
                        const auto [image_position, inserted] =
                            images.emplace(
                                key,
                                Preimage{
                                    states,
                                    moves,
                                    cut_first,
                                    cut_second
                                }
                            );
                        if (!inserted) {
                            ++collisions;
                            if (!printed_collision) {
                                const auto print_preimage =
                                    [&](const char* name,
                                        const Preimage& preimage) {
                                        std::cout
                                            << ' ' << name
                                            << "_cut=("
                                            << preimage.cut_first << ','
                                            << preimage.cut_second << ')'
                                            << ' ' << name << "_moves=";
                                        for (const char move :
                                             preimage.moves) {
                                            std::cout << move;
                                        }
                                        std::cout
                                            << ' ' << name << "_states=";
                                        for (const State& state :
                                             preimage.states) {
                                            std::cout
                                                << '(' << state.first
                                                << ',' << state.second
                                                << ')';
                                        }
                                    };
                                std::cout
                                    << tag
                                    << " first_collision"
                                    << " model=" << model
                                    << " y_moves=" << y_count;
                                print_preimage(
                                    "prior", image_position->second
                                );
                                print_preimage(
                                    "current",
                                    Preimage{
                                        states,
                                        moves,
                                        cut_first,
                                        cut_second
                                    }
                                );
                                std::cout << " image_moves=";
                                for (const char move : image_moves) {
                                    std::cout << move;
                                }
                                std::cout << " image_states=";
                                for (const State& state : image_states) {
                                    std::cout
                                        << '(' << state.first
                                        << ',' << state.second << ')';
                                }
                                std::cout << '\n';
                                printed_collision = true;
                            }
                        }
                        return;
                    }

                    const State source = states.back();
                    const bool force_y = step == length - 1;
                    if (!force_y) {
                        for (const int target :
                             adjacency[
                                 static_cast<std::size_t>(source.first)
                             ]) {
                            states.push_back({target, source.second});
                            moves.push_back('X');
                            visit(step + 1, y_count);
                            moves.pop_back();
                            states.pop_back();
                        }
                    }
                    for (const int target :
                         adjacency[
                             static_cast<std::size_t>(source.second)
                         ]) {
                        states.push_back({source.first, target});
                        moves.push_back('Y');
                        visit(step + 1, y_count + 1);
                        moves.pop_back();
                        states.pop_back();
                    }
                };
            visit(0, 0);
        }

        if (capacity_matching) {
            const int left_size =
                static_cast<int>(matching_edges.size());
            const int right_size =
                static_cast<int>(capacity_right_supply.size());
            if (
                capacity_left_demand.size()
                    != static_cast<std::size_t>(left_size)
                || matching_negative_y.size()
                    != static_cast<std::size_t>(left_size)
                || matching_positive_y.size()
                    != static_cast<std::size_t>(right_size)
                || matching_positive_ids.size()
                    != static_cast<std::size_t>(right_size)
            ) {
                throw std::runtime_error(
                    "projected capacity arrays have inconsistent sizes"
                );
            }
            std::uint64_t total_demand = 0U;
            std::uint64_t isolated_demand = 0U;
            int isolated_types = 0;
            for (int left = 0; left < left_size; ++left) {
                const std::uint64_t demand =
                    capacity_left_demand[static_cast<std::size_t>(left)];
                const std::uint64_t expected = binomial(
                    length - 1,
                    matching_negative_y[static_cast<std::size_t>(left)] - 1
                );
                if (demand != expected) {
                    throw std::runtime_error(
                        "enumerated demand disagrees with binomial demand"
                    );
                }
                if (total_demand
                    > std::numeric_limits<std::uint64_t>::max() - demand) {
                    throw std::runtime_error("total demand exceeds uint64");
                }
                total_demand += demand;
                if (matching_edges[static_cast<std::size_t>(left)].empty()) {
                    ++isolated_types;
                    isolated_demand += demand;
                }
            }
            if (total_demand != negative_paths) {
                throw std::runtime_error(
                    "projected demands do not sum to negative path count"
                );
            }
            std::uint64_t total_supply = 0U;
            for (const std::uint64_t supply : capacity_right_supply) {
                if (total_supply
                    > std::numeric_limits<std::uint64_t>::max() - supply) {
                    throw std::runtime_error("total supply exceeds uint64");
                }
                total_supply += supply;
            }
            if (total_demand
                == std::numeric_limits<std::uint64_t>::max()) {
                throw std::runtime_error(
                    "projected infinite capacity exceeds uint64"
                );
            }
            const std::uint64_t infinite_capacity = total_demand + 1U;
            const int source = 0;
            const int left_base = 1;
            const int right_base = left_base + left_size;
            const int sink = right_base + right_size;
            Dinic flow_graph(sink + 1);
            for (int left = 0; left < left_size; ++left) {
                flow_graph.add_edge(
                    source,
                    left_base + left,
                    capacity_left_demand[static_cast<std::size_t>(left)]
                );
                for (const int right :
                     matching_edges[static_cast<std::size_t>(left)]) {
                    flow_graph.add_edge(
                        left_base + left,
                        right_base + right,
                        infinite_capacity
                    );
                }
            }
            for (int right = 0; right < right_size; ++right) {
                flow_graph.add_edge(
                    right_base + right,
                    sink,
                    capacity_right_supply[
                        static_cast<std::size_t>(right)
                    ]
                );
            }
            const std::uint64_t matched =
                flow_graph.max_flow(source, sink);
            const std::vector<bool> reachable =
                flow_graph.residual_reachable(source);
            std::vector<std::vector<std::uint64_t>> flow_by_y(
                static_cast<std::size_t>(maximum_negative_y + 1),
                std::vector<std::uint64_t>(
                    static_cast<std::size_t>(maximum_positive_y + 1),
                    0U
                )
            );
            for (int left = 0; left < left_size; ++left) {
                const std::vector<int>& targets =
                    matching_edges[static_cast<std::size_t>(left)];
                for (std::size_t edge = 0; edge < targets.size(); ++edge) {
                    const std::uint64_t sent = flow_graph.flow_sent(
                        left_base + left,
                        static_cast<int>(edge) + 1
                    );
                    const int target_y = matching_positive_y[
                        static_cast<std::size_t>(targets[edge])
                    ];
                    flow_by_y[
                        static_cast<std::size_t>(
                            matching_negative_y[
                                static_cast<std::size_t>(left)
                            ]
                        )
                    ][static_cast<std::size_t>(target_y)] += sent;
                }
            }
            std::uint64_t hall_demand = 0U;
            std::uint64_t hall_supply = 0U;
            int hall_left_types = 0;
            int hall_right_types = 0;
            std::vector<std::uint64_t> hall_left_by_y(
                static_cast<std::size_t>(maximum_negative_y + 1),
                0U
            );
            std::vector<std::uint64_t> hall_right_by_y(
                static_cast<std::size_t>(maximum_positive_y + 1),
                0U
            );
            for (int left = 0; left < left_size; ++left) {
                if (!reachable[
                        static_cast<std::size_t>(left_base + left)
                    ]) {
                    continue;
                }
                ++hall_left_types;
                const std::uint64_t demand =
                    capacity_left_demand[static_cast<std::size_t>(left)];
                hall_demand += demand;
                hall_left_by_y[
                    static_cast<std::size_t>(
                        matching_negative_y[
                            static_cast<std::size_t>(left)
                        ]
                    )
                ] += demand;
            }
            for (int right = 0; right < right_size; ++right) {
                if (!reachable[
                        static_cast<std::size_t>(right_base + right)
                    ]) {
                    continue;
                }
                ++hall_right_types;
                const std::uint64_t supply =
                    capacity_right_supply[static_cast<std::size_t>(right)];
                hall_supply += supply;
                hall_right_by_y[
                    static_cast<std::size_t>(
                        matching_positive_y[
                            static_cast<std::size_t>(right)
                        ]
                    )
                ] += supply;
            }
            std::cout
                << tag
                << " level=" << level
                << " label=" << label
                << " prefix=" << prefix
                << " truncation=" << truncation
                << " negative_paths=" << negative_paths
                << " negative_types=" << left_size
                << " positive_types=" << right_size
                << " projected_edges=" << matching_edge_count
                << " reached_supply=" << total_supply
                << " rejected_invalid_path_pairs="
                    << matching_invalid_path_pairs
                << " isolated_types=" << isolated_types
                << " isolated_demand=" << isolated_demand
                << " matched=" << matched
                << " unmatched=" << total_demand - matched
                << " hall_left_types=" << hall_left_types
                << " hall_right_types=" << hall_right_types
                << " hall_demand=" << hall_demand
                << " hall_supply=" << hall_supply
                << " hall_deficit=" << hall_demand - hall_supply
                << " result="
                << (
                    matched == total_demand
                        ? "PASS_CAPACITY"
                        : "FAIL_CAPACITY"
                )
                << '\n';
            std::cout << tag << " hall_demand_by_y=";
            for (int y = 1; y <= maximum_negative_y; y += 2) {
                std::cout
                    << (y == 1 ? "" : ",")
                    << y << ':'
                    << hall_left_by_y[static_cast<std::size_t>(y)];
            }
            std::cout << " hall_supply_by_y=";
            for (int y = 0; y <= maximum_positive_y; y += 2) {
                std::cout
                    << (y == 0 ? "" : ",")
                    << y << ':'
                    << hall_right_by_y[static_cast<std::size_t>(y)];
            }
            std::cout << '\n';
            std::cout << tag << " flow_by_y=";
            bool first_flow = true;
            for (int negative_y = 1;
                 negative_y <= maximum_negative_y;
                 negative_y += 2) {
                for (int positive_y = 0;
                     positive_y <= maximum_positive_y;
                     positive_y += 2) {
                    const std::uint64_t sent = flow_by_y[
                        static_cast<std::size_t>(negative_y)
                    ][static_cast<std::size_t>(positive_y)];
                    if (sent == 0U) {
                        continue;
                    }
                    std::cout
                        << (first_flow ? "" : ",")
                        << negative_y << "->" << positive_y
                        << ':' << sent;
                    first_flow = false;
                }
            }
            std::cout << '\n';
            return matched == total_demand ? EXIT_SUCCESS : EXIT_FAILURE;
        }

        if (cumulative_matching) {
            std::set<std::vector<int>> positive_universe;
            std::vector<std::uint64_t> positive_universe_by_y(
                static_cast<std::size_t>(maximum_positive_y + 1),
                0U
            );
            for (int model = 0; model <= 1; ++model) {
                std::vector<State> positive_states;
                std::vector<char> positive_moves;
                positive_states.push_back({
                    0, model == 0 ? 0 : level
                });
                std::function<void(int, int)> enumerate_positive =
                    [&](int path_step, int y_count) {
                        if (y_count > maximum_positive_y) {
                            return;
                        }
                        if (path_step == length) {
                            const State expected{
                                level, model == 0 ? 0 : level
                            };
                            if (
                                positive_states.back().first
                                    == expected.first
                                && positive_states.back().second
                                    == expected.second
                                && positive_moves.back() == 'X'
                                && (y_count & 1) == 0
                            ) {
                                const std::vector<int> key = flatten(
                                    model,
                                    positive_states,
                                    positive_moves
                                );
                                if (positive_universe.insert(key).second) {
                                    ++positive_universe_by_y[
                                        static_cast<std::size_t>(y_count)
                                    ];
                                }
                            }
                            return;
                        }
                        const State source = positive_states.back();
                        for (const int target :
                             adjacency[
                                 static_cast<std::size_t>(source.first)
                             ]) {
                            positive_states.push_back({
                                target, source.second
                            });
                            positive_moves.push_back('X');
                            enumerate_positive(path_step + 1, y_count);
                            positive_moves.pop_back();
                            positive_states.pop_back();
                        }
                        if (path_step != length - 1) {
                            for (const int target :
                                 adjacency[
                                     static_cast<std::size_t>(
                                         source.second
                                     )
                                 ]) {
                                positive_states.push_back({
                                    source.first, target
                                });
                                positive_moves.push_back('Y');
                                enumerate_positive(
                                    path_step + 1,
                                    y_count + 1
                                );
                                positive_moves.pop_back();
                                positive_states.pop_back();
                            }
                        }
                    };
                enumerate_positive(0, 0);
            }
            for (const std::vector<int>& key : matching_positive_keys) {
                if (!positive_universe.contains(key)) {
                    throw std::runtime_error(
                        "generated matching node is not a positive path"
                    );
                }
            }

            const int left_size =
                static_cast<int>(matching_edges.size());
            const int right_size =
                static_cast<int>(matching_positive_ids.size());
            std::vector<int> left_match(
                static_cast<std::size_t>(left_size),
                -1
            );
            std::vector<int> right_match(
                static_cast<std::size_t>(right_size),
                -1
            );
            std::vector<int> distance(
                static_cast<std::size_t>(left_size),
                -1
            );
            const auto breadth_first = [&]() {
                std::queue<int> queue;
                for (int left = 0; left < left_size; ++left) {
                    if (left_match[static_cast<std::size_t>(left)] < 0) {
                        distance[static_cast<std::size_t>(left)] = 0;
                        queue.push(left);
                    } else {
                        distance[static_cast<std::size_t>(left)] = -1;
                    }
                }
                bool found = false;
                while (!queue.empty()) {
                    const int left = queue.front();
                    queue.pop();
                    for (const int right :
                         matching_edges[static_cast<std::size_t>(left)]) {
                        const int next_left =
                            right_match[static_cast<std::size_t>(right)];
                        if (next_left < 0) {
                            found = true;
                        } else if (
                            distance[
                                static_cast<std::size_t>(next_left)
                            ] < 0
                        ) {
                            distance[
                                static_cast<std::size_t>(next_left)
                            ] = distance[
                                static_cast<std::size_t>(left)
                            ] + 1;
                            queue.push(next_left);
                        }
                    }
                }
                return found;
            };
            std::function<bool(int)> augment = [&](int left) {
                for (const int right :
                     matching_edges[static_cast<std::size_t>(left)]) {
                    const int next_left =
                        right_match[static_cast<std::size_t>(right)];
                    if (next_left < 0
                        || (
                            distance[
                                static_cast<std::size_t>(next_left)
                            ] == distance[
                                static_cast<std::size_t>(left)
                            ] + 1
                            && augment(next_left)
                        )) {
                        left_match[static_cast<std::size_t>(left)] = right;
                        right_match[static_cast<std::size_t>(right)] = left;
                        return true;
                    }
                }
                distance[static_cast<std::size_t>(left)] = -1;
                return false;
            };

            int matched = 0;
            while (breadth_first()) {
                for (int left = 0; left < left_size; ++left) {
                    if (left_match[static_cast<std::size_t>(left)] < 0
                        && augment(left)) {
                        ++matched;
                    }
                }
            }
            std::size_t isolated = 0U;
            for (const std::vector<int>& edge_list : matching_edges) {
                if (edge_list.empty()) {
                    ++isolated;
                }
            }
            std::vector<bool> reachable_left(
                static_cast<std::size_t>(left_size),
                false
            );
            std::vector<bool> reachable_right(
                static_cast<std::size_t>(right_size),
                false
            );
            std::queue<int> alternating_queue;
            for (int left = 0; left < left_size; ++left) {
                if (left_match[static_cast<std::size_t>(left)] < 0) {
                    reachable_left[static_cast<std::size_t>(left)] = true;
                    alternating_queue.push(left);
                }
            }
            while (!alternating_queue.empty()) {
                const int left = alternating_queue.front();
                alternating_queue.pop();
                for (const int right :
                     matching_edges[static_cast<std::size_t>(left)]) {
                    if (left_match[static_cast<std::size_t>(left)]
                            == right
                        || reachable_right[
                            static_cast<std::size_t>(right)
                        ]) {
                        continue;
                    }
                    reachable_right[static_cast<std::size_t>(right)] = true;
                    const int next_left =
                        right_match[static_cast<std::size_t>(right)];
                    if (next_left >= 0
                        && !reachable_left[
                            static_cast<std::size_t>(next_left)
                        ]) {
                        reachable_left[
                            static_cast<std::size_t>(next_left)
                        ] = true;
                        alternating_queue.push(next_left);
                    }
                }
            }
            std::vector<std::uint64_t> hall_left_by_y(
                static_cast<std::size_t>(maximum_negative_y + 1),
                0U
            );
            std::vector<std::uint64_t> hall_right_by_y(
                static_cast<std::size_t>(maximum_positive_y + 1),
                0U
            );
            std::size_t hall_left = 0U;
            std::size_t hall_right = 0U;
            for (int left = 0; left < left_size; ++left) {
                if (reachable_left[static_cast<std::size_t>(left)]) {
                    ++hall_left;
                    ++hall_left_by_y[
                        static_cast<std::size_t>(
                            matching_negative_y[
                                static_cast<std::size_t>(left)
                            ]
                        )
                    ];
                }
            }
            for (int right = 0; right < right_size; ++right) {
                if (reachable_right[static_cast<std::size_t>(right)]) {
                    ++hall_right;
                    ++hall_right_by_y[
                        static_cast<std::size_t>(
                            matching_positive_y[
                                static_cast<std::size_t>(right)
                            ]
                        )
                    ];
                }
            }
            std::vector<int> first_unreachable_positive;
            for (const std::vector<int>& key : positive_universe) {
                if (matching_positive_ids.find(key)
                    == matching_positive_ids.end()) {
                    int key_y = 0;
                    for (int index = 0; index < key[1]; ++index) {
                        key_y += key[
                            2U + 3U * static_cast<std::size_t>(index) + 2U
                        ];
                    }
                    if (key_y == 2) {
                        first_unreachable_positive = key;
                        break;
                    }
                }
            }
            std::cout
                << tag
                << " level=" << level
                << " label=" << label
                << " prefix=" << prefix
                << " truncation=" << truncation
                << " negative_paths=" << negative_paths
                << " positive_paths=" << positive_universe.size()
                << " positive_nodes=" << right_size
                << " unreachable_positive="
                    << positive_universe.size()
                        - matching_positive_ids.size()
                << " splice_edges=" << matching_edge_count
                << " rejected_invalid_path_pairs="
                    << matching_invalid_path_pairs
                << " isolated_negative=" << isolated
                << " matched=" << matched
                << " unmatched=" << left_size - matched
                << " hall_left=" << hall_left
                << " hall_right=" << hall_right
                << " hall_deficit=" << hall_left - hall_right
                << " result="
                << (
                    matched == left_size
                        ? "PASS_MATCHING"
                        : "FAIL_MATCHING"
                )
                << '\n';
            std::cout << tag << " hall_left_by_y=";
            for (int y = 1; y <= maximum_negative_y; y += 2) {
                std::cout
                    << (y == 1 ? "" : ",")
                    << y << ':'
                    << hall_left_by_y[static_cast<std::size_t>(y)];
            }
            std::cout << " hall_right_by_y=";
            for (int y = 0; y <= maximum_positive_y; y += 2) {
                std::cout
                    << (y == 0 ? "" : ",")
                    << y << ':'
                    << hall_right_by_y[static_cast<std::size_t>(y)];
            }
            std::cout << '\n';
            const auto print_preimage =
                [&](const char* name, const Preimage& preimage) {
                    std::cout
                        << tag << ' ' << name << "_moves=";
                    for (const char move : preimage.moves) {
                        std::cout << move;
                    }
                    std::cout << ' ' << name << "_states=";
                    for (const State& state : preimage.states) {
                        std::cout
                            << '(' << state.first
                            << ',' << state.second << ')';
                    }
                    std::cout << '\n';
                };
            for (int left = 0; left < left_size; ++left) {
                if (left_match[static_cast<std::size_t>(left)] < 0) {
                    print_preimage(
                        "first_unmatched",
                        matching_negative_objects[
                            static_cast<std::size_t>(left)
                        ]
                    );
                    break;
                }
            }
            if (!first_unreachable_positive.empty()) {
                const int object_length =
                    first_unreachable_positive[1];
                Preimage unreachable{
                    {},
                    {},
                    -1,
                    -1
                };
                unreachable.states.reserve(
                    static_cast<std::size_t>(object_length + 1)
                );
                unreachable.moves.reserve(
                    static_cast<std::size_t>(object_length)
                );
                for (int index = 0; index < object_length; ++index) {
                    const std::size_t offset =
                        2U + 3U * static_cast<std::size_t>(index);
                    unreachable.states.push_back({
                        first_unreachable_positive[offset],
                        first_unreachable_positive[offset + 1U]
                    });
                    unreachable.moves.push_back(
                        first_unreachable_positive[offset + 2U] == 0
                            ? 'X'
                            : 'Y'
                    );
                }
                const std::size_t final_offset =
                    2U
                    + 3U * static_cast<std::size_t>(object_length);
                unreachable.states.push_back({
                    first_unreachable_positive[final_offset],
                    first_unreachable_positive[final_offset + 1U]
                });
                print_preimage(
                    "first_unreachable_positive",
                    unreachable
                );
            }
            std::cout << tag << " positive_universe_by_y=";
            for (int y = 0; y <= maximum_positive_y; y += 2) {
                std::cout
                    << (y == 0 ? "" : ",")
                    << y << ':'
                    << positive_universe_by_y[
                        static_cast<std::size_t>(y)
                    ];
            }
            std::cout << '\n';
            return matched == left_size ? EXIT_SUCCESS : EXIT_FAILURE;
        }

        std::cout
            << tag
            << " level=" << level
            << " label=" << label
            << " prefix=" << prefix
            << " truncation=" << truncation
            << " negative_paths=" << negative_paths
            << " covered_paths=" << covered_paths
            << " uncovered_paths=" << uncovered_paths
            << " invalid_images=" << invalid_images
            << " collisions=" << collisions
            << " result="
            << (
                uncovered_paths == 0
                    && invalid_images == 0
                    && collisions == 0
                    ? "PASS_INJECTION"
                    : "FAIL_INJECTION"
            )
            << '\n';
        return uncovered_paths == 0
                && invalid_images == 0
                && collisions == 0
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_BPC_SUFFIX_REFLECTION FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
