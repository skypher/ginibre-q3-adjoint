#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Path = std::vector<int>;

int parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0U;
    const long parsed = std::stol(value, &consumed);
    if (consumed != value.size() || parsed <= 0L || parsed > 8L) {
        throw std::invalid_argument(name + " must be an integer in [1,8]");
    }
    return static_cast<int>(parsed);
}

int parse_positive_index(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0U;
    const long parsed = std::stol(value, &consumed);
    if (consumed != value.size() || parsed <= 0L
        || parsed > static_cast<long>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument(name + " must be a positive int");
    }
    return static_cast<int>(parsed);
}

bool edge(int left, int right, int label) {
    return std::abs(left - label) <= right && right <= left + label;
}

void extend_paths(
    const std::vector<int>& labels,
    std::size_t position,
    Path& path,
    std::vector<Path>& paths
) {
    if (position == labels.size()) {
        paths.push_back(path);
        return;
    }
    const int source = path.back();
    const int label = labels[position];
    for (int target = std::abs(source - label);
         target <= source + label;
         ++target) {
        path.push_back(target);
        extend_paths(labels, position + 1U, path, paths);
        path.pop_back();
    }
}

std::vector<Path> paths_from(int start, const std::vector<int>& labels) {
    std::vector<Path> paths;
    Path path;
    path.reserve(labels.size() + 1U);
    path.push_back(start);
    extend_paths(labels, 0U, path, paths);
    return paths;
}

bool two_step_support(int left, int right, int first, int second) {
    return std::abs(left - right) <= first + second
        && left + right >= std::abs(first - second);
}

bool valid_path(const Path& path, const std::vector<int>& labels) {
    if (path.size() != labels.size() + 1U) {
        return false;
    }
    for (std::size_t time = 0U; time < labels.size(); ++time) {
        if (!edge(path[time], path[time + 1U], labels[time])) {
            return false;
        }
    }
    return true;
}

std::string show(const std::vector<int>& values) {
    std::string result{"("};
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            result.push_back(',');
        }
        result += std::to_string(values[index]);
    }
    result.push_back(')');
    return result;
}

std::string encode(const Path& returning, const Path& outgoing) {
    return show(returning) + "|" + show(outgoing);
}

int image_id(
    const Path& returning,
    const Path& outgoing,
    std::map<std::string, int>& ids
) {
    const std::string key = encode(returning, outgoing);
    const auto insertion = ids.emplace(key, static_cast<int>(ids.size()));
    return insertion.first->second;
}

void add_one_block_images(
    const Path& long_path,
    const Path& short_path,
    const std::vector<int>& labels,
    std::vector<int>& images,
    std::map<std::string, int>& ids
) {
    for (std::size_t time = 0U; time < labels.size(); ++time) {
        const int label = labels[time];
        if (!edge(short_path[time], long_path[time + 1U], label)
            || !edge(long_path[time], short_path[time + 1U], label)) {
            continue;
        }
        Path returning(labels.size() + 1U);
        Path outgoing(labels.size() + 1U);
        for (std::size_t index = 0U; index <= labels.size(); ++index) {
            returning[index] = index <= time
                ? short_path[index]
                : long_path[index];
            outgoing[index] = index <= time
                ? long_path[index]
                : short_path[index];
        }
        if (!valid_path(returning, labels) || !valid_path(outgoing, labels)) {
            throw std::runtime_error("invalid one-block image");
        }
        images.push_back(image_id(returning, outgoing, ids));
    }
}

void add_two_block_images(
    const Path& long_path,
    const Path& short_path,
    const std::vector<int>& labels,
    std::vector<int>& images,
    std::map<std::string, int>& ids
) {
    for (std::size_t time = 0U; time + 1U < labels.size(); ++time) {
        const int first = labels[time];
        const int second = labels[time + 1U];
        if (!two_step_support(
                short_path[time], long_path[time + 2U], first, second
            )
            || !two_step_support(
                long_path[time], short_path[time + 2U], first, second
            )) {
            continue;
        }
        for (int returning_middle = std::abs(short_path[time] - first);
             returning_middle <= short_path[time] + first;
             ++returning_middle) {
            if (!edge(returning_middle, long_path[time + 2U], second)) {
                continue;
            }
            for (int outgoing_middle = std::abs(long_path[time] - first);
                 outgoing_middle <= long_path[time] + first;
                 ++outgoing_middle) {
                if (!edge(outgoing_middle, short_path[time + 2U], second)) {
                    continue;
                }
                Path returning(labels.size() + 1U);
                Path outgoing(labels.size() + 1U);
                for (std::size_t index = 0U;
                     index <= labels.size();
                     ++index) {
                    if (index <= time) {
                        returning[index] = short_path[index];
                        outgoing[index] = long_path[index];
                    } else if (index == time + 1U) {
                        returning[index] = returning_middle;
                        outgoing[index] = outgoing_middle;
                    } else {
                        returning[index] = long_path[index];
                        outgoing[index] = short_path[index];
                    }
                }
                if (!valid_path(returning, labels)
                    || !valid_path(outgoing, labels)) {
                    throw std::runtime_error("invalid two-block image");
                }
                images.push_back(image_id(returning, outgoing, ids));
            }
        }
    }
}

void extend_segment_paths(
    const std::vector<int>& labels,
    std::size_t start,
    std::size_t length,
    std::size_t step,
    int endpoint,
    Path& segment,
    std::vector<Path>& segments
) {
    if (step == length) {
        if (segment.back() == endpoint) {
            segments.emplace_back();
            Path& output = segments.back();
            output.reserve(segment.size());
            for (const int state : segment) {
                output.push_back(state);
            }
        }
        return;
    }
    const int source = segment.back();
    const int label = labels[start + step];
    for (int target = std::abs(source - label);
         target <= source + label;
         ++target) {
        segment.push_back(target);
        extend_segment_paths(
            labels,
            start,
            length,
            step + 1U,
            endpoint,
            segment,
            segments
        );
        segment.pop_back();
    }
}

std::vector<Path> segment_paths(
    int start_state,
    int endpoint,
    const std::vector<int>& labels,
    std::size_t start,
    std::size_t length
) {
    std::vector<Path> segments;
    Path segment;
    segment.reserve(length + 1U);
    segment.push_back(start_state);
    extend_segment_paths(
        labels,
        start,
        length,
        0U,
        endpoint,
        segment,
        segments
    );
    return segments;
}

void add_wide_block_images(
    const Path& long_path,
    const Path& short_path,
    const std::vector<int>& labels,
    int maximum_blocks,
    std::vector<int>& images,
    std::map<std::string, int>& ids
) {
    for (std::size_t length = 3U;
         length <= static_cast<std::size_t>(maximum_blocks)
             && length <= labels.size();
         ++length) {
        for (std::size_t time = 0U;
             time + length <= labels.size();
             ++time) {
            const std::vector<Path> returning_segments = segment_paths(
                short_path[time],
                long_path[time + length],
                labels,
                time,
                length
            );
            const std::vector<Path> outgoing_segments = segment_paths(
                long_path[time],
                short_path[time + length],
                labels,
                time,
                length
            );
            for (const Path& returning_segment : returning_segments) {
                for (const Path& outgoing_segment : outgoing_segments) {
                    Path returning(labels.size() + 1U);
                    Path outgoing(labels.size() + 1U);
                    for (std::size_t index = 0U;
                         index <= labels.size();
                         ++index) {
                        if (index <= time) {
                            returning[index] = short_path[index];
                            outgoing[index] = long_path[index];
                        } else if (index < time + length) {
                            returning[index] =
                                returning_segment[index - time];
                            outgoing[index] =
                                outgoing_segment[index - time];
                        } else {
                            returning[index] = long_path[index];
                            outgoing[index] = short_path[index];
                        }
                    }
                    if (!valid_path(returning, labels)
                        || !valid_path(outgoing, labels)) {
                        throw std::runtime_error("invalid wide-block image");
                    }
                    images.push_back(image_id(returning, outgoing, ids));
                }
            }
        }
    }
}

struct Matching {
    std::size_t size = 0U;
    std::size_t hall_left = 0U;
    std::size_t hall_right = 0U;
    std::vector<bool> hall_left_membership;
};

Matching maximum_matching(const std::vector<std::vector<int>>& graph,
                          std::size_t image_count) {
    if (graph.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || image_count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("matching graph exceeds int indexing");
    }
    std::vector<int> matched_left(graph.size(), -1);
    std::vector<int> matched_right(image_count, -1);
    std::vector<int> distance(graph.size(), -1);
    const auto augment = [&](auto&& self, int left) -> bool {
        for (const int right : graph[static_cast<std::size_t>(left)]) {
            const int previous = matched_right[static_cast<std::size_t>(right)];
            if (previous < 0
                || (distance[static_cast<std::size_t>(previous)]
                        == distance[static_cast<std::size_t>(left)] + 1
                    && self(self, previous))) {
                matched_left[static_cast<std::size_t>(left)] = right;
                matched_right[static_cast<std::size_t>(right)] = left;
                return true;
            }
        }
        distance[static_cast<std::size_t>(left)] = -1;
        return false;
    };

    Matching result;
    while (true) {
        std::queue<int> queue;
        for (int left = 0; left < static_cast<int>(graph.size()); ++left) {
            if (matched_left[static_cast<std::size_t>(left)] < 0) {
                distance[static_cast<std::size_t>(left)] = 0;
                queue.push(left);
            } else {
                distance[static_cast<std::size_t>(left)] = -1;
            }
        }
        bool reaches_free_right = false;
        while (!queue.empty()) {
            const int left = queue.front();
            queue.pop();
            for (const int right : graph[static_cast<std::size_t>(left)]) {
                const int next = matched_right[static_cast<std::size_t>(right)];
                if (next < 0) {
                    reaches_free_right = true;
                } else if (distance[static_cast<std::size_t>(next)] < 0) {
                    distance[static_cast<std::size_t>(next)] =
                        distance[static_cast<std::size_t>(left)] + 1;
                    queue.push(next);
                }
            }
        }
        if (!reaches_free_right) {
            break;
        }
        std::size_t augmented = 0U;
        for (int left = 0; left < static_cast<int>(graph.size()); ++left) {
            if (matched_left[static_cast<std::size_t>(left)] < 0
                && augment(augment, left)) {
                ++augmented;
            }
        }
        if (augmented == 0U) {
            break;
        }
        result.size += augmented;
    }
    std::vector<bool> reachable_left(graph.size(), false);
    std::vector<bool> reachable_right(image_count, false);
    std::queue<int> alternating_queue;
    for (int left = 0; left < static_cast<int>(graph.size()); ++left) {
        if (matched_left[static_cast<std::size_t>(left)] < 0) {
            reachable_left[static_cast<std::size_t>(left)] = true;
            alternating_queue.push(left);
        }
    }
    while (!alternating_queue.empty()) {
        const int left = alternating_queue.front();
        alternating_queue.pop();
        for (const int right : graph[static_cast<std::size_t>(left)]) {
            if (matched_left[static_cast<std::size_t>(left)] == right
                || reachable_right[static_cast<std::size_t>(right)]) {
                continue;
            }
            reachable_right[static_cast<std::size_t>(right)] = true;
            const int next = matched_right[static_cast<std::size_t>(right)];
            if (next >= 0 && !reachable_left[static_cast<std::size_t>(next)]) {
                reachable_left[static_cast<std::size_t>(next)] = true;
                alternating_queue.push(next);
            }
        }
    }
    for (const bool reached : reachable_left) {
        if (reached) {
            ++result.hall_left;
        }
    }
    for (const bool reached : reachable_right) {
        if (reached) {
            ++result.hall_right;
        }
    }
    result.hall_left_membership = std::move(reachable_left);
    return result;
}

struct Statistics {
    unsigned long long instances = 0U;
    unsigned long long domain_vertices = 0U;
    unsigned long long image_vertices = 0U;
    unsigned long long incidences = 0U;
    unsigned long long matching_failures = 0U;
    std::size_t first_hall_left = 0U;
    std::size_t first_hall_right = 0U;
    std::string witness;
    std::vector<std::string> failures;
    std::vector<std::string> first_hall_vertices;
    std::vector<std::pair<Path, Path>> first_hall_pairs;
    std::string first_three_block_repair;
};

void inspect_instance(
    int boundary,
    int target,
    const std::vector<Path>& long_paths,
    const std::vector<Path>& short_paths,
    const std::vector<int>& labels,
    int maximum_blocks,
    bool print_hall,
    bool summarize_hall,
    bool repair_hall_three,
    Statistics& statistics
) {
    std::map<std::string, int> ids;
    std::vector<std::vector<int>> graph;
    std::vector<std::string> domain_keys;
    std::vector<std::pair<Path, Path>> domain_pairs;
    for (const Path& long_path : long_paths) {
        for (const Path& short_path : short_paths) {
            std::vector<int> images;
            add_one_block_images(long_path, short_path, labels, images, ids);
            add_two_block_images(long_path, short_path, labels, images, ids);
            if (maximum_blocks >= 3) {
                add_wide_block_images(
                    long_path,
                    short_path,
                    labels,
                    maximum_blocks,
                    images,
                    ids
                );
            }
            std::sort(images.begin(), images.end());
            images.erase(std::unique(images.begin(), images.end()), images.end());
            statistics.incidences +=
                static_cast<unsigned long long>(images.size());
            graph.push_back(std::move(images));
            domain_keys.push_back(encode(long_path, short_path));
            domain_pairs.emplace_back(long_path, short_path);
            if (repair_hall_three && graph.size() % 100000U == 0U) {
                std::cerr << "SU2_HALL_REPAIR_BUILD domain=" << graph.size()
                          << '\n' << std::flush;
            }
        }
    }
    if (graph.empty()) {
        return;
    }
    ++statistics.instances;
    statistics.domain_vertices += static_cast<unsigned long long>(graph.size());
    statistics.image_vertices += static_cast<unsigned long long>(ids.size());
    const Matching matching = maximum_matching(graph, ids.size());
    if (matching.size != graph.size()) {
        ++statistics.matching_failures;
        const std::string failure = "R=" + std::to_string(boundary)
            + " S=" + std::to_string(target)
            + " domain=" + std::to_string(graph.size())
            + " image=" + std::to_string(ids.size())
            + " matching=" + std::to_string(matching.size)
            + " hall_left=" + std::to_string(matching.hall_left)
            + " hall_right=" + std::to_string(matching.hall_right);
        if (statistics.failures.size() < 16U) {
            statistics.failures.push_back(failure);
        }
        if (statistics.witness.empty()) {
            statistics.first_hall_left = matching.hall_left;
            statistics.first_hall_right = matching.hall_right;
            statistics.witness = failure;
            if (print_hall || summarize_hall) {
                for (std::size_t index = 0U;
                     index < matching.hall_left_membership.size();
                     ++index) {
                    if (matching.hall_left_membership[index]) {
                        if (print_hall) {
                            statistics.first_hall_vertices.push_back(
                                domain_keys[index]
                            );
                        }
                        if (summarize_hall) {
                            statistics.first_hall_pairs.push_back(
                                domain_pairs[index]
                            );
                        }
                    }
                }
            }
            if (repair_hall_three) {
                std::map<std::string, int> repair_ids;
                std::vector<std::vector<int>> repair_graph;
                repair_graph.reserve(matching.hall_left);
                for (std::size_t index = 0U;
                     index < matching.hall_left_membership.size();
                     ++index) {
                    if (!matching.hall_left_membership[index]) {
                        continue;
                    }
                    const Path& returning = domain_pairs[index].first;
                    const Path& outgoing = domain_pairs[index].second;
                    std::vector<int> images;
                    add_one_block_images(
                        returning,
                        outgoing,
                        labels,
                        images,
                        repair_ids
                    );
                    add_two_block_images(
                        returning,
                        outgoing,
                        labels,
                        images,
                        repair_ids
                    );
                    add_wide_block_images(
                        returning,
                        outgoing,
                        labels,
                        3,
                        images,
                        repair_ids
                    );
                    std::sort(images.begin(), images.end());
                    images.erase(
                        std::unique(images.begin(), images.end()),
                        images.end()
                    );
                    repair_graph.push_back(std::move(images));
                    if (repair_graph.size() % 10000U == 0U) {
                        std::cerr << "SU2_HALL_REPAIR_BUILD selected="
                                  << repair_graph.size()
                                  << " images=" << repair_ids.size()
                                  << '\n' << std::flush;
                    }
                }
                const Matching repair_matching = maximum_matching(
                    repair_graph,
                    repair_ids.size()
                );
                std::size_t repair_width_two_targets = 0U;
                for (const auto& entry : repair_ids) {
                    if (ids.contains(entry.first)) {
                        ++repair_width_two_targets;
                    }
                }
                statistics.first_three_block_repair =
                    "domain=" + std::to_string(repair_graph.size())
                    + " image=" + std::to_string(repair_ids.size())
                    + " width_two_targets="
                    + std::to_string(repair_width_two_targets)
                    + " three_block_only_targets="
                    + std::to_string(
                        repair_ids.size() - repair_width_two_targets
                    )
                    + " matching=" + std::to_string(repair_matching.size)
                    + " result=" + (repair_matching.size == repair_graph.size()
                        ? "PASS"
                        : "FAIL");
            }
        }
    }
}

void print_hall_summary(const Statistics& statistics) {
    if (statistics.first_hall_pairs.empty()) {
        return;
    }
    std::set<Path> returning_paths;
    std::set<Path> outgoing_paths;
    std::set<std::string> hall_pair_keys;
    for (const auto& paths : statistics.first_hall_pairs) {
        returning_paths.insert(paths.first);
        outgoing_paths.insert(paths.second);
        hall_pair_keys.insert(encode(paths.first, paths.second));
    }
    const std::size_t product_size =
        returning_paths.size() * outgoing_paths.size();
    std::cout << " hall_returning_paths=" << returning_paths.size()
              << " hall_outgoing_paths=" << outgoing_paths.size()
              << " hall_cartesian_product_size=" << product_size
              << " hall_is_cartesian="
              << (product_size == statistics.first_hall_pairs.size()
                  ? "YES"
                  : "NO");
    if (product_size != statistics.first_hall_pairs.size()) {
        for (const Path& returning : returning_paths) {
            bool found = false;
            for (const Path& outgoing : outgoing_paths) {
                if (hall_pair_keys.contains(encode(returning, outgoing))) {
                    continue;
                }
                const auto returning_witness = std::find_if(
                    statistics.first_hall_pairs.begin(),
                    statistics.first_hall_pairs.end(),
                    [&returning](const std::pair<Path, Path>& paths) {
                        return paths.first == returning;
                    }
                );
                const auto outgoing_witness = std::find_if(
                    statistics.first_hall_pairs.begin(),
                    statistics.first_hall_pairs.end(),
                    [&outgoing](const std::pair<Path, Path>& paths) {
                        return paths.second == outgoing;
                    }
                );
                if (returning_witness == statistics.first_hall_pairs.end()
                    || outgoing_witness == statistics.first_hall_pairs.end()) {
                    throw std::runtime_error("missing Hall projection witness");
                }
                std::cout << " hall_missing_cross="
                          << '[' << encode(returning, outgoing) << ']'
                          << " hall_returning_projection_witness="
                          << '[' << encode(
                              returning_witness->first,
                              returning_witness->second
                          ) << ']'
                          << " hall_outgoing_projection_witness="
                          << '[' << encode(
                              outgoing_witness->first,
                              outgoing_witness->second
                          ) << ']';
                found = true;
                break;
            }
            if (found) {
                break;
            }
        }
    }

    const std::size_t path_size =
        statistics.first_hall_pairs.front().first.size();
    std::cout << " hall_returning_coordinate_ranges=";
    for (std::size_t index = 0U; index < path_size; ++index) {
        int minimum = statistics.first_hall_pairs.front().first[index];
        int maximum = minimum;
        std::set<int> values;
        for (const auto& paths : statistics.first_hall_pairs) {
            const int value = paths.first[index];
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
            values.insert(value);
        }
        std::cout << '[' << minimum << ',' << maximum << ';';
        for (const int value : values) {
            std::cout << value << ',';
        }
        std::cout << ']';
    }
    std::cout << " hall_outgoing_coordinate_ranges=";
    for (std::size_t index = 0U; index < path_size; ++index) {
        int minimum = statistics.first_hall_pairs.front().second[index];
        int maximum = minimum;
        std::set<int> values;
        for (const auto& paths : statistics.first_hall_pairs) {
            const int value = paths.second[index];
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
            values.insert(value);
        }
        std::cout << '[' << minimum << ',' << maximum << ';';
        for (const int value : values) {
            std::cout << value << ',';
        }
        std::cout << ']';
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        int argument_end = argc;
        bool print_hall = false;
        bool summarize_hall = false;
        bool repair_hall_three = false;
        while (argument_end >= 2) {
            const std::string flag{argv[argument_end - 1]};
            if (flag == "--print-hall") {
                print_hall = true;
                --argument_end;
            } else if (flag == "--summarize-hall") {
                summarize_hall = true;
                --argument_end;
            } else if (flag == "--repair-hall-three") {
                repair_hall_three = true;
                --argument_end;
            } else {
                break;
            }
        }
        bool fixed_instance = false;
        int fixed_boundary = 0;
        int fixed_target = 0;
        int maximum_blocks = 2;
        int root_start = 2;
        if (argument_end >= 7 && std::string(argv[1]) == "--instance") {
            if (std::string(argv[4]) != "--maximum-blocks"
                || std::string(argv[6]) != "--paired-root") {
                throw std::invalid_argument(
                    "instance mode requires --instance R S --maximum-blocks "
                    "L --paired-root LABEL... [--print-hall] "
                    "[--summarize-hall] [--repair-hall-three]"
                );
            }
            fixed_instance = true;
            fixed_boundary = parse_positive_index(argv[2], "R");
            fixed_target = parse_positive_index(argv[3], "S");
            maximum_blocks = parse_positive_index(argv[5], "L");
            root_start = 7;
        } else if (argument_end < 3 || std::string(argv[1]) != "--paired-root") {
            throw std::invalid_argument(
                "usage: probe_su2_heterogeneous_bridge_matching "
                "--paired-root LABEL... or --instance R S --maximum-blocks "
                "L --paired-root LABEL... [--print-hall] "
                "[--summarize-hall] [--repair-hall-three]"
            );
        }
        std::vector<int> root;
        root.reserve(static_cast<std::size_t>(argument_end - root_start));
        for (int index = root_start; index < argument_end; ++index) {
            root.push_back(parse_positive(argv[index], "root label"));
        }
        if (root.size() > 4U) {
            throw std::invalid_argument(
                "the exact matching diagnostic accepts at most four root labels"
            );
        }
        std::vector<int> labels;
        labels.reserve(2U * root.size());
        for (const int label : root) {
            labels.push_back(label);
            labels.push_back(label);
        }
        int total_label = 0;
        for (const int label : labels) {
            total_label += label;
        }
        if (maximum_blocks > static_cast<int>(labels.size())) {
            throw std::invalid_argument("L exceeds the paired word length");
        }
        if (repair_hall_three && maximum_blocks != 2) {
            throw std::invalid_argument(
                "--repair-hall-three requires a width-two source graph"
            );
        }
        if (repair_hall_three && !fixed_instance) {
            throw std::invalid_argument(
                "--repair-hall-three requires one fixed R,S instance"
            );
        }
        const std::vector<Path> all_short_paths = paths_from(0, labels);
        Statistics statistics;
        const int first_boundary = fixed_instance ? fixed_boundary : 1;
        const int last_boundary = fixed_instance ? fixed_boundary : total_label;
        for (int boundary = first_boundary; boundary <= last_boundary; ++boundary) {
            const std::vector<Path> all_long_paths =
                paths_from(boundary, labels);
            std::vector<Path> long_paths;
            for (const Path& path : all_long_paths) {
                if (path.back() == 0) {
                    long_paths.push_back(path);
                }
            }
            if (long_paths.empty()) {
                continue;
            }
            const int first_target = fixed_instance ? fixed_target : 1;
            const int last_target = fixed_instance ? fixed_target : total_label;
            for (int target = first_target; target <= last_target; ++target) {
                std::vector<Path> short_paths;
                for (const Path& path : all_short_paths) {
                    if (path.back() == target) {
                        short_paths.push_back(path);
                    }
                }
                if (!short_paths.empty()) {
                    inspect_instance(
                        boundary,
                        target,
                        long_paths,
                        short_paths,
                        labels,
                        maximum_blocks,
                        print_hall,
                        summarize_hall,
                        repair_hall_three,
                        statistics
                    );
                }
            }
        }
        std::cout << "SU2_HETEROGENEOUS_PAIRED_BRIDGE_MATCHING"
                  << " root=" << show(root)
                  << " word=" << show(labels)
                  << " maximum_blocks=" << maximum_blocks
                  << " instances=" << statistics.instances
                  << " domain_vertices=" << statistics.domain_vertices
                  << " image_vertices=" << statistics.image_vertices
                  << " incidences=" << statistics.incidences
                  << " matching_failures=" << statistics.matching_failures
                  << " witness=" << statistics.witness
                  << " failures=";
        for (const std::string& failure : statistics.failures) {
            std::cout << '[' << failure << ']';
        }
        std::cout << " hall_vertices=";
        for (const std::string& vertex : statistics.first_hall_vertices) {
            std::cout << '[' << vertex << ']';
        }
        if (summarize_hall) {
            print_hall_summary(statistics);
        }
        if (repair_hall_three) {
            std::cout << " three_block_hall_repair="
                      << statistics.first_three_block_repair;
        }
        std::cout
                  << " result=PASS"
                  << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
