#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

struct Task {
    int q = 0;
    int length = 0;
};

struct Counts {
    unsigned long long pairs = 0;
    unsigned long long good = 0;
    unsigned long long bad = 0;
    unsigned long long collisions = 0;
    unsigned long long all_switch_candidates = 0;
    unsigned long long all_switch_no_candidate = 0;
    unsigned long long all_switch_matching = 0;
    unsigned long long all_switch_unmatched = 0;
    unsigned long long bridge_candidates = 0;
    unsigned long long bridge_no_candidate = 0;
    unsigned long long bridge_matching = 0;
    unsigned long long bridge_unmatched = 0;
    unsigned long long bridge_components = 0;
    unsigned long long bridge_deficient_components = 0;
    unsigned long long bridge_complete_components = 0;
    unsigned long long bridge_edge_degree_reversals = 0;
    unsigned long long bridge_left_reciprocal_failures = 0;
    unsigned long long bridge_right_reciprocal_failures = 0;
    unsigned long long bridge_labeled_edge_degree_reversals = 0;
    unsigned long long bridge_labeled_left_reciprocal_failures = 0;
    unsigned long long bridge_labeled_right_reciprocal_failures = 0;
    unsigned long long capacity_no_window = 0;
    std::string witness;
    std::string no_candidate_witness;
    std::string bridge_no_candidate_witness;
    std::string capacity_no_window_witness;
};

class UnionFind {
public:
    explicit UnionFind(std::size_t size)
        : parent_(size), rank_(size, 0) {
        for (std::size_t index = 0; index < size; ++index) {
            parent_[index] = static_cast<int>(index);
        }
    }

    int find(int value) {
        int root = value;
        while (parent_[static_cast<std::size_t>(root)] != root) {
            root = parent_[static_cast<std::size_t>(root)];
        }
        while (value != root) {
            const int next = parent_[static_cast<std::size_t>(value)];
            parent_[static_cast<std::size_t>(value)] = root;
            value = next;
        }
        return root;
    }

    void unite(int left, int right) {
        int left_root = find(left);
        int right_root = find(right);
        if (left_root == right_root) {
            return;
        }
        if (
            rank_[static_cast<std::size_t>(left_root)]
            < rank_[static_cast<std::size_t>(right_root)]
        ) {
            std::swap(left_root, right_root);
        }
        parent_[static_cast<std::size_t>(right_root)] = left_root;
        if (
            rank_[static_cast<std::size_t>(left_root)]
            == rank_[static_cast<std::size_t>(right_root)]
        ) {
            ++rank_[static_cast<std::size_t>(left_root)];
        }
    }

private:
    std::vector<int> parent_;
    std::vector<int> rank_;
};

struct ComponentCounts {
    unsigned long long left = 0;
    unsigned long long right = 0;
    unsigned long long edges = 0;
};

int parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0;
    const long long parsed = std::stoll(value, &consumed);
    if (consumed != value.size() || parsed <= 0) {
        throw std::invalid_argument(name + " must be positive");
    }
    return static_cast<int>(parsed);
}

unsigned int adaptive_threads() {
    const unsigned int hardware = std::max(1U, std::thread::hardware_concurrency());
    double load[1] = {0.0};
    const int read = getloadavg(load, 1);
    const unsigned int occupied =
        read == 1
        ? static_cast<unsigned int>(std::ceil(std::max(0.0, load[0])))
        : 0U;
    return std::max(
        1U,
        hardware > occupied + 2U ? hardware - occupied - 2U : 1U
    );
}

bool edge(int x, int y, int q) {
    return std::abs(x - q) <= y && y <= x + q;
}

int two_step_count(int x, int y, int q) {
    const int lower = std::max(std::abs(x - q), std::abs(y - q));
    const int upper = std::min(x + q, y + q);
    return std::max(0, upper - lower + 1);
}

void enumerate_from(
    int q,
    int steps,
    std::vector<int>& path,
    std::vector<std::vector<int>>& paths
) {
    if (static_cast<int>(path.size()) == steps + 1) {
        paths.push_back(path);
        return;
    }
    const int x = path.back();
    for (int y = std::abs(x - q); y <= x + q; ++y) {
        path.push_back(y);
        enumerate_from(q, steps, path, paths);
        path.pop_back();
    }
}

std::vector<std::vector<int>> enumerate_paths(int q, int steps) {
    std::vector<std::vector<int>> paths;
    std::vector<int> path{0};
    enumerate_from(q, steps, path, paths);
    return paths;
}

std::vector<int> refine(const std::vector<int>& path, int q) {
    std::vector<int> heights;
    heights.reserve(
        1U
        + static_cast<std::size_t>(2 * q)
            * (path.size() - 1U)
    );
    int height = 2 * path.front();
    heights.push_back(height);
    for (std::size_t t = 0; t + 1U < path.size(); ++t) {
        const int downs = path[t] + q - path[t + 1U];
        for (int step = 0; step < downs; ++step) {
            --height;
            heights.push_back(height);
        }
        for (int step = downs; step < 2 * q; ++step) {
            ++height;
            heights.push_back(height);
        }
    }
    return heights;
}

std::string encode(
    const std::vector<int>& short_path,
    const std::vector<int>& long_path
) {
    std::string key;
    for (int state : short_path) {
        key += std::to_string(state);
        key.push_back(',');
    }
    key.push_back('|');
    for (int state : long_path) {
        key += std::to_string(state);
        key.push_back(',');
    }
    return key;
}

std::string show_pair(
    const std::vector<int>& left,
    const std::vector<int>& right
) {
    std::string text = "{";
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (i != 0U) {
            text.push_back(',');
        }
        text += std::to_string(left[i]);
    }
    text += "};{";
    for (std::size_t i = 0; i < right.size(); ++i) {
        if (i != 0U) {
            text.push_back(',');
        }
        text += std::to_string(right[i]);
    }
    text.push_back('}');
    return text;
}

Counts check_task(const Task& task) {
    const int q = task.q;
    const int n = task.length;
    const int block = 2 * q;
    const auto short_paths = enumerate_paths(q, n);
    const auto all_long_paths = enumerate_paths(q, n + 1);
    std::vector<std::vector<int>> returning_long;
    for (const auto& path : all_long_paths) {
        if (path.back() == 0) {
            returning_long.push_back(path);
        }
    }

    Counts counts;
    std::map<std::string, std::string> images;
    std::map<std::string, int> all_image_ids;
    std::vector<std::vector<int>> all_adjacency;
    std::map<std::string, int> bridge_image_ids;
    std::vector<std::vector<int>> bridge_adjacency;
    std::vector<std::vector<int>> bridge_labeled_adjacency;
    for (const auto& a : returning_long) {
        const auto refined_a = refine(a, q);
        std::vector<int> aligned_a(
            refined_a.begin() + block,
            refined_a.end()
        );
        for (const auto& b : short_paths) {
            if (b.back() == 0) {
                continue;
            }
            ++counts.pairs;
            const auto refined_b = refine(b, q);
            std::vector<std::size_t> meetings;
            for (
                std::size_t time = 1U;
                time < aligned_a.size();
                ++time
            ) {
                if (aligned_a[time] == refined_b[time]) {
                    meetings.push_back(time);
                }
            }
            if (meetings.empty()) {
                counts.witness = "missing_meeting";
                ++counts.bad;
                all_adjacency.emplace_back();
                bridge_adjacency.emplace_back();
                bridge_labeled_adjacency.emplace_back();
                continue;
            }

            auto switched_image = [&](
                std::size_t meeting,
                std::vector<int>& short_output,
                std::vector<int>& long_output
            ) {
                short_output.assign(
                    static_cast<std::size_t>(n + 1),
                    0
                );
                long_output.assign(
                    static_cast<std::size_t>(n + 2),
                    0
                );
                for (int t = 0; t <= n; ++t) {
                    const std::size_t time =
                        static_cast<std::size_t>(t * block);
                    const int short_height =
                        time <= meeting
                        ? refined_b[time]
                        : aligned_a[time];
                    const int long_height =
                        time <= meeting
                        ? aligned_a[time]
                        : refined_b[time];
                    short_output[static_cast<std::size_t>(t)] =
                        short_height / 2;
                    long_output[static_cast<std::size_t>(t + 1)] =
                        long_height / 2;
                }

                bool valid = true;
                for (int t = 0; t < n; ++t) {
                    valid =
                        valid
                        && edge(
                            short_output[static_cast<std::size_t>(t)],
                            short_output[
                                static_cast<std::size_t>(t + 1)
                            ],
                            q
                        );
                }
                for (int t = 0; t <= n; ++t) {
                    valid =
                        valid
                        && edge(
                            long_output[static_cast<std::size_t>(t)],
                            long_output[
                                static_cast<std::size_t>(t + 1)
                            ],
                            q
                        );
                }
                return valid;
            };

            std::vector<int> short_output;
            std::vector<int> long_output;
            const bool first_valid = switched_image(
                meetings.front(),
                short_output,
                long_output
            );
            if (!first_valid) {
                ++counts.bad;
            } else {
                ++counts.good;
                const std::string key =
                    encode(short_output, long_output);
                const std::string input = show_pair(a, b);
                const auto [iterator, inserted] =
                    images.emplace(key, input);
                if (!inserted) {
                    ++counts.collisions;
                    if (counts.witness.empty()) {
                        counts.witness =
                            "first=" + iterator->second
                            + " second=" + input
                            + " image=" + key;
                    }
                }
            }

            std::vector<int> adjacency;
            for (const std::size_t meeting : meetings) {
                if (!switched_image(
                        meeting,
                        short_output,
                        long_output
                    )) {
                    continue;
                }
                const std::string key =
                    encode(short_output, long_output);
                const auto [iterator, inserted] =
                    all_image_ids.emplace(
                        key,
                        static_cast<int>(all_image_ids.size())
                    );
                adjacency.push_back(iterator->second);
            }
            std::sort(adjacency.begin(), adjacency.end());
            adjacency.erase(
                std::unique(adjacency.begin(), adjacency.end()),
                adjacency.end()
            );
            if (adjacency.empty()) {
                ++counts.all_switch_no_candidate;
                if (counts.no_candidate_witness.empty()) {
                    counts.no_candidate_witness = show_pair(a, b);
                }
            }
            counts.all_switch_candidates += adjacency.size();
            all_adjacency.push_back(std::move(adjacency));

            std::vector<int> bridge_edges;
            auto add_bridge = [&](
                int start,
                int length,
                int short_middle,
                int long_middle
            ) {
                std::vector<int> short_bridge(
                    static_cast<std::size_t>(n + 1)
                );
                std::vector<int> long_bridge(
                    static_cast<std::size_t>(n + 2)
                );
                long_bridge[0] = 0;
                for (int time = 0; time <= n; ++time) {
                    if (time <= start) {
                        short_bridge[static_cast<std::size_t>(time)] =
                            b[static_cast<std::size_t>(time)];
                        long_bridge[
                            static_cast<std::size_t>(time + 1)
                        ] = a[static_cast<std::size_t>(time + 1)];
                    } else if (time >= start + length) {
                        short_bridge[static_cast<std::size_t>(time)] =
                            a[static_cast<std::size_t>(time + 1)];
                        long_bridge[
                            static_cast<std::size_t>(time + 1)
                        ] = b[static_cast<std::size_t>(time)];
                    } else {
                        short_bridge[static_cast<std::size_t>(time)] =
                            short_middle;
                        long_bridge[
                            static_cast<std::size_t>(time + 1)
                        ] = long_middle;
                    }
                }
                const std::string key =
                    encode(short_bridge, long_bridge);
                const auto insertion = bridge_image_ids.emplace(
                    key,
                    static_cast<int>(bridge_image_ids.size())
                );
                bridge_edges.push_back(insertion.first->second);
            };

            for (int start = 0; start < n; ++start) {
                if (
                    edge(
                        b[static_cast<std::size_t>(start)],
                        a[static_cast<std::size_t>(start + 2)],
                        q
                    )
                    && edge(
                        a[static_cast<std::size_t>(start + 1)],
                        b[static_cast<std::size_t>(start + 1)],
                        q
                    )
                ) {
                    add_bridge(start, 1, 0, 0);
                }
            }
            const int max_state = (n + 1) * q;
            bool has_capacity_window = false;
            for (int start = 0; start + 2 <= n; ++start) {
                const int short_start =
                    b[static_cast<std::size_t>(start)];
                const int short_end =
                    a[static_cast<std::size_t>(start + 3)];
                const int long_start =
                    a[static_cast<std::size_t>(start + 1)];
                const int long_end =
                    b[static_cast<std::size_t>(start + 2)];
                const int straight_capacity =
                    two_step_count(short_start, long_end, q)
                    * two_step_count(long_start, short_end, q);
                const int crossed_capacity =
                    two_step_count(short_start, short_end, q)
                    * two_step_count(long_start, long_end, q);
                if (crossed_capacity >= straight_capacity) {
                    has_capacity_window = true;
                }
                for (
                    int short_middle = 0;
                    short_middle <= max_state;
                    ++short_middle
                ) {
                    if (
                        !edge(short_start, short_middle, q)
                        || !edge(short_middle, short_end, q)
                    ) {
                        continue;
                    }
                    for (
                        int long_middle = 0;
                        long_middle <= max_state;
                        ++long_middle
                    ) {
                        if (
                            edge(long_start, long_middle, q)
                            && edge(long_middle, long_end, q)
                        ) {
                            add_bridge(
                                start,
                                2,
                                short_middle,
                                long_middle
                            );
                        }
                    }
                }
            }
            if (!has_capacity_window) {
                ++counts.capacity_no_window;
                if (counts.capacity_no_window_witness.empty()) {
                    counts.capacity_no_window_witness =
                        show_pair(a, b);
                }
            }
            const std::vector<int> labeled_bridge_edges = bridge_edges;
            std::sort(bridge_edges.begin(), bridge_edges.end());
            bridge_edges.erase(
                std::unique(bridge_edges.begin(), bridge_edges.end()),
                bridge_edges.end()
            );
            if (bridge_edges.empty()) {
                ++counts.bridge_no_candidate;
                if (counts.bridge_no_candidate_witness.empty()) {
                    counts.bridge_no_candidate_witness =
                        show_pair(a, b);
                }
            }
            counts.bridge_candidates += bridge_edges.size();
            bridge_labeled_adjacency.push_back(labeled_bridge_edges);
            bridge_adjacency.push_back(std::move(bridge_edges));
        }
    }

    auto matching_size = [](
        const std::vector<std::vector<int>>& graph,
        std::size_t image_count
    ) {
        std::vector<int> matched_image(graph.size(), -1);
        std::vector<int> matched_domain(image_count, -1);
        std::vector<int> distance(graph.size(), -1);
        auto augment = [&](auto&& self, int domain) -> bool {
            for (
                const int image :
                graph[static_cast<std::size_t>(domain)]
            ) {
                const int previous =
                    matched_domain[static_cast<std::size_t>(image)];
                if (
                    previous < 0
                    || (
                        distance[static_cast<std::size_t>(previous)]
                            == distance[static_cast<std::size_t>(domain)] + 1
                        && self(self, previous)
                    )
                ) {
                    matched_image[static_cast<std::size_t>(domain)] =
                        image;
                    matched_domain[static_cast<std::size_t>(image)] =
                        domain;
                    return true;
                }
            }
            distance[static_cast<std::size_t>(domain)] = -1;
            return false;
        };
        unsigned long long matching = 0;
        while (true) {
            std::queue<int> queue;
            for (
                int domain = 0;
                domain < static_cast<int>(graph.size());
                ++domain
            ) {
                if (
                    matched_image[static_cast<std::size_t>(domain)] < 0
                ) {
                    distance[static_cast<std::size_t>(domain)] = 0;
                    queue.push(domain);
                } else {
                    distance[static_cast<std::size_t>(domain)] = -1;
                }
            }
            bool reaches_free_image = false;
            while (!queue.empty()) {
                const int domain = queue.front();
                queue.pop();
                for (
                    const int image :
                    graph[static_cast<std::size_t>(domain)]
                ) {
                    const int next =
                        matched_domain[static_cast<std::size_t>(image)];
                    if (next < 0) {
                        reaches_free_image = true;
                    } else if (
                        distance[static_cast<std::size_t>(next)] < 0
                    ) {
                        distance[static_cast<std::size_t>(next)] =
                            distance[static_cast<std::size_t>(domain)] + 1;
                        queue.push(next);
                    }
                }
            }
            if (!reaches_free_image) {
                break;
            }
            unsigned long long augmented = 0;
            for (
                int domain = 0;
                domain < static_cast<int>(graph.size());
                ++domain
            ) {
                if (
                    matched_image[static_cast<std::size_t>(domain)] < 0
                    && augment(augment, domain)
                ) {
                    ++augmented;
                }
            }
            if (augmented == 0) {
                break;
            }
            matching += augmented;
        }
        return matching;
    };
    counts.all_switch_matching =
        matching_size(all_adjacency, all_image_ids.size());
    counts.all_switch_unmatched =
        counts.pairs - counts.all_switch_matching;
    counts.bridge_matching =
        matching_size(bridge_adjacency, bridge_image_ids.size());
    counts.bridge_unmatched =
        counts.pairs - counts.bridge_matching;

    const int left_size = static_cast<int>(bridge_adjacency.size());
    UnionFind components(
        bridge_adjacency.size() + bridge_image_ids.size()
    );
    for (int left = 0; left < left_size; ++left) {
        for (
            const int right :
            bridge_adjacency[static_cast<std::size_t>(left)]
        ) {
            components.unite(left, left_size + right);
        }
    }
    std::map<int, ComponentCounts> component_counts;
    for (int left = 0; left < left_size; ++left) {
        ++component_counts[components.find(left)].left;
    }
    for (
        int right = 0;
        right < static_cast<int>(bridge_image_ids.size());
        ++right
    ) {
        ++component_counts[components.find(left_size + right)].right;
    }
    for (int left = 0; left < left_size; ++left) {
        for (
            const int right :
            bridge_adjacency[static_cast<std::size_t>(left)]
        ) {
            ++component_counts[
                components.find(left_size + right)
            ].edges;
        }
    }
    counts.bridge_components = component_counts.size();
    for (const auto& [root, component] : component_counts) {
        static_cast<void>(root);
        if (component.right < component.left) {
            ++counts.bridge_deficient_components;
        }
        if (
            component.edges
            == component.left * component.right
        ) {
            ++counts.bridge_complete_components;
        }
    }
    std::vector<unsigned long long> right_degrees(
        bridge_image_ids.size(),
        0
    );
    for (const auto& neighbors : bridge_adjacency) {
        for (const int right : neighbors) {
            ++right_degrees[static_cast<std::size_t>(right)];
        }
    }
    std::vector<long double> right_reciprocal_load(
        bridge_image_ids.size(),
        0.0L
    );
    for (const auto& neighbors : bridge_adjacency) {
        long double left_load = 0.0L;
        const auto left_degree =
            static_cast<unsigned long long>(neighbors.size());
        for (const int right : neighbors) {
            const auto right_degree =
                right_degrees[static_cast<std::size_t>(right)];
            if (left_degree < right_degree) {
                ++counts.bridge_edge_degree_reversals;
            }
            left_load += 1.0L / static_cast<long double>(right_degree);
            right_reciprocal_load[static_cast<std::size_t>(right)]
                += 1.0L / static_cast<long double>(left_degree);
        }
        if (left_load + 1.0e-18L < 1.0L) {
            ++counts.bridge_left_reciprocal_failures;
        }
    }
    for (const long double load : right_reciprocal_load) {
        if (load > 1.0L + 1.0e-18L) {
            ++counts.bridge_right_reciprocal_failures;
        }
    }
    std::fill(right_degrees.begin(), right_degrees.end(), 0);
    for (const auto& neighbors : bridge_labeled_adjacency) {
        for (const int right : neighbors) {
            ++right_degrees[static_cast<std::size_t>(right)];
        }
    }
    std::fill(
        right_reciprocal_load.begin(),
        right_reciprocal_load.end(),
        0.0L
    );
    for (const auto& neighbors : bridge_labeled_adjacency) {
        long double left_load = 0.0L;
        const auto left_degree =
            static_cast<unsigned long long>(neighbors.size());
        for (const int right : neighbors) {
            const auto right_degree =
                right_degrees[static_cast<std::size_t>(right)];
            if (left_degree < right_degree) {
                ++counts.bridge_labeled_edge_degree_reversals;
            }
            left_load += 1.0L / static_cast<long double>(right_degree);
            right_reciprocal_load[static_cast<std::size_t>(right)]
                += 1.0L / static_cast<long double>(left_degree);
        }
        if (left_load + 1.0e-18L < 1.0L) {
            ++counts.bridge_labeled_left_reciprocal_failures;
        }
    }
    for (const long double load : right_reciprocal_load) {
        if (load > 1.0L + 1.0e-18L) {
            ++counts.bridge_labeled_right_reciprocal_failures;
        }
    }
    return counts;
}

}  // namespace

int main(int argc, char** argv) {
    int max_q = 3;
    int max_length = 4;
    try {
        if (argc >= 2) {
            max_q = parse_positive(argv[1], "max_q");
        }
        if (argc >= 3) {
            max_length = parse_positive(argv[2], "max_length");
        }
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_kostka_crossing_switch [max_q] [max_length]"
            );
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }

    std::vector<Task> tasks;
    for (int q = 1; q <= max_q; ++q) {
        for (int length = 2; length <= max_length; length += 2) {
            tasks.push_back({q, length});
        }
    }
    std::vector<Counts> results(tasks.size());
    std::atomic<std::size_t> next{0U};
    const unsigned int thread_count = std::min(
        adaptive_threads(),
        static_cast<unsigned int>(std::max<std::size_t>(1U, tasks.size()))
    );
    auto worker = [&]() {
        while (true) {
            const std::size_t index = next.fetch_add(1U);
            if (index >= tasks.size()) {
                return;
            }
            results[index] = check_task(tasks[index]);
        }
    };
    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    for (unsigned int thread = 0; thread < thread_count; ++thread) {
        workers.emplace_back(worker);
    }
    for (std::thread& thread : workers) {
        thread.join();
    }

    unsigned long long pairs = 0;
    unsigned long long good = 0;
    unsigned long long bad = 0;
    unsigned long long collisions = 0;
    unsigned long long all_switch_candidates = 0;
    unsigned long long all_switch_no_candidate = 0;
    unsigned long long all_switch_matching = 0;
    unsigned long long all_switch_unmatched = 0;
    unsigned long long bridge_candidates = 0;
    unsigned long long bridge_no_candidate = 0;
    unsigned long long bridge_matching = 0;
    unsigned long long bridge_unmatched = 0;
    unsigned long long bridge_components = 0;
    unsigned long long bridge_deficient_components = 0;
    unsigned long long bridge_complete_components = 0;
    unsigned long long bridge_edge_degree_reversals = 0;
    unsigned long long bridge_left_reciprocal_failures = 0;
    unsigned long long bridge_right_reciprocal_failures = 0;
    unsigned long long bridge_labeled_edge_degree_reversals = 0;
    unsigned long long bridge_labeled_left_reciprocal_failures = 0;
    unsigned long long bridge_labeled_right_reciprocal_failures = 0;
    unsigned long long capacity_no_window = 0;
    for (std::size_t i = 0; i < tasks.size(); ++i) {
        pairs += results[i].pairs;
        good += results[i].good;
        bad += results[i].bad;
        collisions += results[i].collisions;
        all_switch_candidates += results[i].all_switch_candidates;
        all_switch_no_candidate += results[i].all_switch_no_candidate;
        all_switch_matching += results[i].all_switch_matching;
        all_switch_unmatched += results[i].all_switch_unmatched;
        bridge_candidates += results[i].bridge_candidates;
        bridge_no_candidate += results[i].bridge_no_candidate;
        bridge_matching += results[i].bridge_matching;
        bridge_unmatched += results[i].bridge_unmatched;
        bridge_components += results[i].bridge_components;
        bridge_deficient_components +=
            results[i].bridge_deficient_components;
        bridge_complete_components +=
            results[i].bridge_complete_components;
        bridge_edge_degree_reversals +=
            results[i].bridge_edge_degree_reversals;
        bridge_left_reciprocal_failures +=
            results[i].bridge_left_reciprocal_failures;
        bridge_right_reciprocal_failures +=
            results[i].bridge_right_reciprocal_failures;
        bridge_labeled_edge_degree_reversals +=
            results[i].bridge_labeled_edge_degree_reversals;
        bridge_labeled_left_reciprocal_failures +=
            results[i].bridge_labeled_left_reciprocal_failures;
        bridge_labeled_right_reciprocal_failures +=
            results[i].bridge_labeled_right_reciprocal_failures;
        capacity_no_window += results[i].capacity_no_window;
        if (!results[i].witness.empty()) {
            std::cout
                << "SU2_KOSTKA_CROSSING_SWITCH witness"
                << " q=" << tasks[i].q
                << " length=" << tasks[i].length
                << " " << results[i].witness
                << '\n';
        }
        if (!results[i].no_candidate_witness.empty()) {
            std::cout
                << "SU2_KOSTKA_CROSSING_SWITCH no_candidate"
                << " q=" << tasks[i].q
                << " length=" << tasks[i].length
                << " input=" << results[i].no_candidate_witness
                << '\n';
        }
        if (!results[i].bridge_no_candidate_witness.empty()) {
            std::cout
                << "SU2_KOSTKA_CROSSING_SWITCH bridge_no_candidate"
                << " q=" << tasks[i].q
                << " length=" << tasks[i].length
                << " input="
                << results[i].bridge_no_candidate_witness
                << '\n';
        }
        if (!results[i].capacity_no_window_witness.empty()) {
            std::cout
                << "SU2_KOSTKA_CROSSING_SWITCH capacity_no_window"
                << " q=" << tasks[i].q
                << " length=" << tasks[i].length
                << " input="
                << results[i].capacity_no_window_witness
                << '\n';
        }
    }
    std::cout
        << "SU2_KOSTKA_CROSSING_SWITCH"
        << " tasks=" << tasks.size()
        << " pairs=" << pairs
        << " good=" << good
        << " bad=" << bad
        << " collisions=" << collisions
        << " all_switch_candidates=" << all_switch_candidates
        << " all_switch_no_candidate=" << all_switch_no_candidate
        << " all_switch_matching=" << all_switch_matching
        << " all_switch_unmatched=" << all_switch_unmatched
        << " bridge_candidates=" << bridge_candidates
        << " bridge_no_candidate=" << bridge_no_candidate
        << " bridge_matching=" << bridge_matching
        << " bridge_unmatched=" << bridge_unmatched
        << " bridge_components=" << bridge_components
        << " bridge_deficient_components="
        << bridge_deficient_components
        << " bridge_complete_components="
        << bridge_complete_components
        << " bridge_edge_degree_reversals="
        << bridge_edge_degree_reversals
        << " bridge_left_reciprocal_failures="
        << bridge_left_reciprocal_failures
        << " bridge_right_reciprocal_failures="
        << bridge_right_reciprocal_failures
        << " bridge_labeled_edge_degree_reversals="
        << bridge_labeled_edge_degree_reversals
        << " bridge_labeled_left_reciprocal_failures="
        << bridge_labeled_left_reciprocal_failures
        << " bridge_labeled_right_reciprocal_failures="
        << bridge_labeled_right_reciprocal_failures
        << " capacity_no_window=" << capacity_no_window
        << " max_q=" << max_q
        << " max_length=" << max_length
        << " threads=" << thread_count
        << " result="
        << (
            bridge_unmatched == 0
            ? "TWO_BLOCK_MATCHING"
            : "NO_TWO_BLOCK_MATCHING"
        )
        << '\n';
    return EXIT_SUCCESS;
}
