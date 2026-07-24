#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using Key = std::vector<int>;

enum class OrderMode { cyclic, all_permutations };

enum class CarrierMode { first_descent, all_h2_prefixes };

struct EdgeMode {
    OrderMode order;
    CarrierMode carrier;
    bool include_evaluation;
};

struct LeftVertex {
    unsigned active_mask;
    std::vector<int> scalar_grades;
    std::vector<int> canonical_heights;
};

struct Graph {
    std::vector<LeftVertex> left;
    std::vector<std::vector<int>> adjacency;
    std::map<Key, int> right_ids;
    std::vector<Key> right_keys;
};

struct Deficit {
    int scalar_index_grade;
    int matching_rank;
    std::vector<int> hall_left;
    std::vector<int> hall_right;
    Graph graph;
};

struct Configuration {
    int q;
    std::vector<int> plus_labels;
};

struct SearchResult {
    Configuration configuration;
    Deficit deficit;
};

struct ThreadChoice {
    unsigned hardware_threads;
    unsigned ram_threads;
    unsigned load_threads;
    unsigned chosen_threads;
    double load_average;
    unsigned long long available_bytes;
};

void enumerate_paths(const std::vector<int>& plus_labels,
                     const std::vector<int>& order,
                     int position,
                     int height,
                     std::vector<int>& heights,
                     std::vector<std::vector<int>>& paths) {
    if (position == static_cast<int>(order.size())) {
        if (height == 0) {
            paths.push_back(heights);
        }
        return;
    }

    const int factor = order[position];
    const int label = plus_labels[factor];
    for (int next = std::abs(height - label); next <= height + label;
         next += 2) {
        heights[factor] = next;
        enumerate_paths(plus_labels, order, position + 1, next, heights,
                        paths);
    }
}

std::vector<std::vector<int>> fusion_paths(
    int initial_height,
    const std::vector<int>& plus_labels,
    const std::vector<int>& order) {
    std::vector<std::vector<int>> paths;
    std::vector<int> heights(plus_labels.size(), -1);
    enumerate_paths(plus_labels, order, 0, initial_height, heights, paths);
    return paths;
}

std::vector<int> active_indices(unsigned mask, int factor_count) {
    std::vector<int> indices;
    for (int i = 0; i < factor_count; ++i) {
        if ((mask & (1U << i)) != 0U) {
            indices.push_back(i);
        }
    }
    return indices;
}

std::vector<std::vector<int>> allowed_orders(std::vector<int> indices,
                                             OrderMode mode) {
    std::vector<std::vector<int>> orders;
    if (mode == OrderMode::cyclic) {
        for (int shift = 0; shift < static_cast<int>(indices.size());
             ++shift) {
            std::vector<int> order;
            order.insert(order.end(), indices.begin() + shift, indices.end());
            order.insert(order.end(), indices.begin(), indices.begin() + shift);
            orders.push_back(std::move(order));
        }
        return orders;
    }

    do {
        orders.push_back(indices);
    } while (std::next_permutation(indices.begin(), indices.end()));
    return orders;
}

std::set<Key> saturated_targets(int q,
                                const std::vector<int>& plus_labels,
                                unsigned active_mask,
                                const std::vector<int>& scalar_grades,
                                EdgeMode mode) {
    const std::vector<int> active =
        active_indices(active_mask, static_cast<int>(plus_labels.size()));
    const std::vector<std::vector<int>> orders =
        allowed_orders(active, mode.order);
    std::set<Key> targets;

    for (const std::vector<int>& order : orders) {
        for (const std::vector<int>& heights :
             fusion_paths(q, plus_labels, order)) {
            std::vector<int> stopping_positions;
            for (int position = 0;
                 position < static_cast<int>(order.size()); ++position) {
                const int height = heights[order[position]];
                if (height < q) {
                    stopping_positions.push_back(position);
                    if (mode.carrier == CarrierMode::first_descent) {
                        break;
                    }
                }
            }

            for (const int stopping_position : stopping_positions) {
                const int carrier_label = heights[order[stopping_position]];
                std::vector<int> output_status = scalar_grades;
                int previous_height = q;
                for (int position = 0; position <= stopping_position;
                     ++position) {
                    const int factor = order[position];
                    const int height = heights[factor];
                    output_status[factor] =
                        (plus_labels[factor] + height - previous_height) / 2;
                    previous_height = height;
                }

                unsigned suffix_mask = 0U;
                std::vector<int> suffix_order;
                for (int position = stopping_position + 1;
                     position < static_cast<int>(order.size()); ++position) {
                    const int factor = order[position];
                    suffix_mask |= 1U << factor;
                    output_status[factor] = -1;
                    suffix_order.push_back(factor);
                }
                std::sort(suffix_order.begin(), suffix_order.end());

                for (const std::vector<int>& suffix_heights :
                     fusion_paths(carrier_label, plus_labels, suffix_order)) {
                    Key key{1, carrier_label,
                            static_cast<int>(suffix_mask)};
                    key.insert(key.end(), output_status.begin(),
                               output_status.end());
                    key.insert(key.end(), suffix_heights.begin(),
                               suffix_heights.end());
                    targets.insert(std::move(key));
                }
            }
        }
    }

    if (mode.include_evaluation) {
        int active_label_sum = q;
        for (int factor = 0; factor < static_cast<int>(plus_labels.size());
             ++factor) {
            if ((active_mask & (1U << factor)) != 0U) {
                active_label_sum += plus_labels[factor];
            }
        }
        const int required_sum = active_label_sum / 2;
        for (int minus_grade = 0; minus_grade <= q; ++minus_grade) {
            std::vector<int> allocation = scalar_grades;
            const auto recurse = [&](const auto& self,
                                     int factor,
                                     int remaining) -> void {
                if (factor == static_cast<int>(plus_labels.size())) {
                    if (remaining == 0) {
                        Key key{0, minus_grade};
                        key.insert(key.end(), allocation.begin(),
                                   allocation.end());
                        targets.insert(std::move(key));
                    }
                    return;
                }
                if ((active_mask & (1U << factor)) == 0U) {
                    self(self, factor + 1, remaining);
                    return;
                }
                for (int grade = 0;
                     grade <= plus_labels[factor] && grade <= remaining;
                     ++grade) {
                    allocation[factor] = grade;
                    self(self, factor + 1, remaining - grade);
                }
            };
            recurse(recurse, 0, required_sum - minus_grade);
        }
    }
    return targets;
}

void add_left_block(Graph& graph,
                    int q,
                    const std::vector<int>& plus_labels,
                    unsigned active_mask,
                    const std::vector<int>& scalar_grades,
                    const std::vector<std::vector<int>>& source_paths,
                    EdgeMode mode) {
    const std::set<Key> targets =
        saturated_targets(q, plus_labels, active_mask, scalar_grades, mode);
    std::vector<int> neighbors;
    for (const Key& key : targets) {
        const auto [iterator, inserted] =
            graph.right_ids.emplace(key, graph.right_keys.size());
        if (inserted) {
            graph.right_keys.push_back(key);
        }
        neighbors.push_back(iterator->second);
    }

    for (const std::vector<int>& path : source_paths) {
        graph.left.push_back(LeftVertex{active_mask, scalar_grades, path});
        graph.adjacency.push_back(neighbors);
    }
}

void enumerate_scalar_blocks(
    int q,
    const std::vector<int>& plus_labels,
    unsigned active_mask,
    const std::vector<std::vector<int>>& source_paths,
    int factor,
    std::vector<int>& scalar_grades,
    int inactive_grade_sum,
    int active_label_sum,
    EdgeMode mode,
    std::map<int, Graph>& graphs) {
    if (factor == static_cast<int>(plus_labels.size())) {
        if (active_label_sum % 2 != 0) {
            return;
        }
        const int scalar_index_grade =
            active_label_sum / 2 + inactive_grade_sum;
        add_left_block(graphs[scalar_index_grade], q, plus_labels, active_mask,
                       scalar_grades, source_paths, mode);
        return;
    }

    if ((active_mask & (1U << factor)) != 0U) {
        scalar_grades[factor] = -1;
        enumerate_scalar_blocks(
            q, plus_labels, active_mask, source_paths, factor + 1,
            scalar_grades, inactive_grade_sum,
            active_label_sum + plus_labels[factor], mode, graphs);
        return;
    }

    for (int grade = 0; grade <= plus_labels[factor]; ++grade) {
        scalar_grades[factor] = grade;
        enumerate_scalar_blocks(
            q, plus_labels, active_mask, source_paths, factor + 1,
            scalar_grades, inactive_grade_sum + grade, active_label_sum, mode,
            graphs);
    }
}

std::map<int, Graph> build_graphs(int q,
                                  const std::vector<int>& plus_labels,
                                  EdgeMode mode) {
    std::map<int, Graph> graphs;
    const unsigned mask_limit = 1U << plus_labels.size();
    for (unsigned active_mask = 1U; active_mask < mask_limit; ++active_mask) {
        const std::vector<int> canonical_order =
            active_indices(active_mask, static_cast<int>(plus_labels.size()));
        const std::vector<std::vector<int>> source_paths =
            fusion_paths(q, plus_labels, canonical_order);
        if (source_paths.empty()) {
            continue;
        }

        std::vector<int> scalar_grades(plus_labels.size(), -1);
        enumerate_scalar_blocks(q, plus_labels, active_mask, source_paths, 0,
                                scalar_grades, 0, q, mode, graphs);
    }
    return graphs;
}

bool augment(int left,
             const Graph& graph,
             std::vector<int>& seen,
             int stamp,
             std::vector<int>& matched_right) {
    for (const int right : graph.adjacency[left]) {
        if (seen[right] == stamp) {
            continue;
        }
        seen[right] = stamp;
        if (matched_right[right] < 0 ||
            augment(matched_right[right], graph, seen, stamp,
                    matched_right)) {
            matched_right[right] = left;
            return true;
        }
    }
    return false;
}

std::optional<Deficit> hall_deficit(int scalar_index_grade, Graph graph) {
    std::vector<int> matched_right(graph.right_keys.size(), -1);
    std::vector<int> seen(graph.right_keys.size(), 0);
    int rank = 0;
    for (int left = 0; left < static_cast<int>(graph.left.size()); ++left) {
        if (augment(left, graph, seen, left + 1, matched_right)) {
            ++rank;
        }
    }
    if (rank == static_cast<int>(graph.left.size())) {
        return std::nullopt;
    }

    std::vector<int> matched_left(graph.left.size(), -1);
    for (int right = 0; right < static_cast<int>(matched_right.size());
         ++right) {
        if (matched_right[right] >= 0) {
            matched_left[matched_right[right]] = right;
        }
    }

    std::vector<bool> reachable_left(graph.left.size(), false);
    std::vector<bool> reachable_right(graph.right_keys.size(), false);
    std::queue<int> queue;
    for (int left = 0; left < static_cast<int>(graph.left.size()); ++left) {
        if (matched_left[left] < 0) {
            reachable_left[left] = true;
            queue.push(left);
        }
    }

    while (!queue.empty()) {
        const int left = queue.front();
        queue.pop();
        for (const int right : graph.adjacency[left]) {
            if (matched_left[left] == right || reachable_right[right]) {
                continue;
            }
            reachable_right[right] = true;
            const int next_left = matched_right[right];
            if (next_left >= 0 && !reachable_left[next_left]) {
                reachable_left[next_left] = true;
                queue.push(next_left);
            }
        }
    }

    Deficit result;
    result.scalar_index_grade = scalar_index_grade;
    result.matching_rank = rank;
    for (int left = 0; left < static_cast<int>(reachable_left.size());
         ++left) {
        if (reachable_left[left]) {
            result.hall_left.push_back(left);
        }
    }
    for (int right = 0; right < static_cast<int>(reachable_right.size());
         ++right) {
        if (reachable_right[right]) {
            result.hall_right.push_back(right);
        }
    }
    result.graph = std::move(graph);
    return result;
}

std::optional<Deficit> check_configuration(const Configuration& configuration,
                                           EdgeMode mode) {
    std::map<int, Graph> graphs =
        build_graphs(configuration.q, configuration.plus_labels, mode);
    for (auto& [grade, graph] : graphs) {
        std::optional<Deficit> deficit =
            hall_deficit(grade, std::move(graph));
        if (deficit.has_value()) {
            return deficit;
        }
    }
    return std::nullopt;
}

unsigned long long available_memory_bytes() {
    std::ifstream input("/proc/meminfo");
    std::string key;
    unsigned long long kibibytes = 0;
    std::string unit;
    while (input >> key >> kibibytes >> unit) {
        if (key == "MemAvailable:") {
            return kibibytes * 1024ULL;
        }
    }
    return 512ULL * 1024ULL * 1024ULL;
}

ThreadChoice choose_threads() {
    constexpr unsigned long long bytes_per_worker =
        512ULL * 1024ULL * 1024ULL;
    const unsigned hardware = std::max(1U, std::thread::hardware_concurrency());
    const unsigned long long available = available_memory_bytes();
    const unsigned ram_threads = std::max(
        1U, static_cast<unsigned>(available / bytes_per_worker));
    double load = 0.0;
    if (getloadavg(&load, 1) != 1) {
        load = 0.0;
    }
    const unsigned occupied =
        static_cast<unsigned>(std::min<double>(hardware - 1, std::ceil(load)));
    const unsigned load_threads = std::max(1U, hardware - occupied);
    const unsigned chosen =
        std::max(1U, std::min({hardware, ram_threads, load_threads}));
    return ThreadChoice{hardware, ram_threads, load_threads, chosen, load,
                        available};
}

std::vector<Configuration> configurations(int factor_count,
                                          int maximum_q,
                                          int maximum_label) {
    std::vector<Configuration> result;
    std::vector<int> labels(factor_count, 1);
    while (true) {
        for (int q = 1; q <= maximum_q; ++q) {
            if (std::find(labels.begin(), labels.end(), q) == labels.end()) {
                result.push_back(Configuration{q, labels});
            }
        }

        int position = factor_count - 1;
        while (position >= 0 && labels[position] == maximum_label) {
            --position;
        }
        if (position < 0) {
            break;
        }
        ++labels[position];
        for (int i = position + 1; i < factor_count; ++i) {
            labels[i] = 1;
        }
    }
    return result;
}

std::optional<SearchResult> search_mode(int maximum_q,
                                        int maximum_label,
                                        int maximum_factors,
                                        EdgeMode mode,
                                        unsigned thread_count,
                                        std::size_t& checked) {
    checked = 0;
    for (int factor_count = 1; factor_count <= maximum_factors;
         ++factor_count) {
        const std::vector<Configuration> batch =
            configurations(factor_count, maximum_q, maximum_label);
        std::vector<std::optional<Deficit>> results(batch.size());
        std::atomic<std::size_t> next{0};

        const unsigned active_threads = std::max(
            1U, std::min<unsigned>(thread_count, batch.size()));
        std::vector<std::thread> workers;
        workers.reserve(active_threads);
        for (unsigned worker = 0; worker < active_threads; ++worker) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t index = next.fetch_add(1);
                    if (index >= batch.size()) {
                        return;
                    }
                    results[index] = check_configuration(batch[index], mode);
                }
            });
        }
        for (std::thread& worker : workers) {
            worker.join();
        }

        for (std::size_t index = 0; index < batch.size(); ++index) {
            ++checked;
            if (results[index].has_value()) {
                return SearchResult{batch[index],
                                    std::move(*results[index])};
            }
        }
    }
    return std::nullopt;
}

void print_values(const std::vector<int>& values) {
    for (const int value : values) {
        std::cout << value << ',';
    }
}

std::string mode_name(EdgeMode mode) {
    if (mode.order == OrderMode::cyclic &&
        mode.carrier == CarrierMode::first_descent &&
        !mode.include_evaluation) {
        return "CYCLIC_FIRST_DESCENT_SATURATED";
    }
    if (mode.order == OrderMode::all_permutations &&
        mode.carrier == CarrierMode::first_descent &&
        !mode.include_evaluation) {
        return "ALL_PERMUTATIONS_FIRST_DESCENT_SATURATED";
    }
    if (mode.include_evaluation) {
        return "OPTIMISTIC_ESSENTIAL_H1_0_HYBRID_UPPER_GRAPH";
    }
    return "ALL_PERMUTATIONS_ALL_H2_PREFIXES_SATURATED";
}

void print_result(EdgeMode mode,
                  const std::optional<SearchResult>& result,
                  std::size_t checked,
                  int maximum_q,
                  int maximum_label,
                  int maximum_factors) {
    const std::string name = mode_name(mode);
    if (!result.has_value()) {
        std::cout << name << "_BOUNDED_PASS"
                  << " q<=" << maximum_q << " label<=" << maximum_label
                  << " factors<=" << maximum_factors
                  << " configurations=" << checked << '\n';
        return;
    }

    const SearchResult& found = *result;
    const Deficit& deficit = found.deficit;
    std::cout << name << "_HALL_DEFICIT\n";
    std::cout << "q=" << found.configuration.q << " plus_labels=";
    print_values(found.configuration.plus_labels);
    std::cout << " scalar_index_grade=" << deficit.scalar_index_grade << '\n';
    std::cout << "left_total=" << deficit.graph.left.size()
              << " right_total=" << deficit.graph.right_keys.size()
              << " matching_rank=" << deficit.matching_rank << '\n';
    std::cout << "hall_left_size=" << deficit.hall_left.size()
              << " hall_neighborhood_size=" << deficit.hall_right.size()
              << " hall_deficit="
              << deficit.hall_left.size() - deficit.hall_right.size()
              << '\n';
    std::cout << "first_collision_factor_count="
              << found.configuration.plus_labels.size() << '\n';
    std::cout << "minimality_order=factor_count,ordered_plus_word,q,grade\n";
    std::cout << "configurations_checked=" << checked << '\n';

    const std::size_t sample_count =
        std::min<std::size_t>(8, deficit.hall_left.size());
    for (std::size_t i = 0; i < sample_count; ++i) {
        const LeftVertex& left =
            deficit.graph.left[deficit.hall_left[i]];
        std::cout << " hall_left_sample mask=" << left.active_mask
                  << " scalar_grades=";
        print_values(left.scalar_grades);
        std::cout << " canonical_heights=";
        print_values(left.canonical_heights);
        std::cout << '\n';
    }
    const std::size_t right_sample_count =
        std::min<std::size_t>(8, deficit.hall_right.size());
    for (std::size_t i = 0; i < right_sample_count; ++i) {
        std::cout << " hall_right_sample key=";
        print_values(deficit.graph.right_keys[deficit.hall_right[i]]);
        std::cout << '\n';
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) {
        std::cerr << "usage: search_su2_fdc_h2_recoupling_hall "
                     "MAXIMUM_Q MAXIMUM_LABEL MAXIMUM_FACTORS "
                     "[HYBRID_ONLY]\n";
        return 1;
    }

    const int maximum_q = std::stoi(argv[1]);
    const int maximum_label = std::stoi(argv[2]);
    const int maximum_factors = std::stoi(argv[3]);
    if (maximum_q < 1 || maximum_label < 1 || maximum_factors < 1 ||
        maximum_factors >= 8 * static_cast<int>(sizeof(unsigned)) ||
        (argc == 5 && std::string(argv[4]) != "HYBRID_ONLY")) {
        std::cerr << "invalid bound\n";
        return 1;
    }
    const bool hybrid_only = argc == 5;

    const ThreadChoice threads = choose_threads();
    std::cout << "thread_policy hardware=" << threads.hardware_threads
              << " ram_limit=" << threads.ram_threads
              << " load_limit=" << threads.load_threads
              << " chosen=" << threads.chosen_threads
              << " load1=" << threads.load_average
              << " mem_available_bytes=" << threads.available_bytes << '\n';
    std::cout
        << "graph=recoupling-saturated upper graph on essential H1(0) "
           "sources\n";
    std::cout
        << "scope=paired H1(j+1) sources displaced by H2(j) use are excluded; "
           "bounded pass is not a full-differential rank proof\n";

    std::vector<EdgeMode> modes;
    if (!hybrid_only) {
        modes.push_back(
            EdgeMode{OrderMode::cyclic, CarrierMode::first_descent, false});
        modes.push_back(EdgeMode{OrderMode::all_permutations,
                                 CarrierMode::first_descent, false});
        modes.push_back(EdgeMode{OrderMode::all_permutations,
                                 CarrierMode::all_h2_prefixes, false});
    }
    modes.push_back(EdgeMode{OrderMode::all_permutations,
                             CarrierMode::all_h2_prefixes, true});
    for (const EdgeMode mode : modes) {
        std::size_t checked = 0;
        const std::optional<SearchResult> result =
            search_mode(maximum_q, maximum_label, maximum_factors, mode,
                        threads.chosen_threads, checked);
        print_result(mode, result, checked, maximum_q, maximum_label,
                     maximum_factors);
    }
    return 0;
}
