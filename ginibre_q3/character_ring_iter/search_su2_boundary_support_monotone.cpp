#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <queue>
#include <string>
#include <utility>
#include <vector>

using Count = std::uint64_t;

struct Edge {
    int to;
    int reverse;
    Count capacity;
};

class Dinic {
  public:
    explicit Dinic(int vertices)
        : graph_(static_cast<std::size_t>(vertices)),
          level_(static_cast<std::size_t>(vertices)),
          next_(static_cast<std::size_t>(vertices)) {}

    void add_edge(int from, int to, Count capacity) {
        Edge forward{to, static_cast<int>(graph_[static_cast<std::size_t>(to)].size()),
                     capacity};
        Edge reverse{from,
                     static_cast<int>(graph_[static_cast<std::size_t>(from)].size()),
                     0};
        graph_[static_cast<std::size_t>(from)].push_back(forward);
        graph_[static_cast<std::size_t>(to)].push_back(reverse);
    }

    Count flow(int source, int sink) {
        Count answer = 0;
        while (bfs(source, sink)) {
            std::fill(next_.begin(), next_.end(), 0);
            while (Count pushed = dfs(source, sink, std::numeric_limits<Count>::max())) {
                answer += pushed;
            }
        }
        return answer;
    }

    std::vector<bool> reachable(int source) const {
        std::vector<bool> seen(graph_.size(), false);
        std::queue<int> pending;
        seen[static_cast<std::size_t>(source)] = true;
        pending.push(source);
        while (!pending.empty()) {
            const int from = pending.front();
            pending.pop();
            for (const Edge& edge : graph_[static_cast<std::size_t>(from)]) {
                if (edge.capacity == 0 || seen[static_cast<std::size_t>(edge.to)]) {
                    continue;
                }
                seen[static_cast<std::size_t>(edge.to)] = true;
                pending.push(edge.to);
            }
        }
        return seen;
    }

  private:
    bool bfs(int source, int sink) {
        std::fill(level_.begin(), level_.end(), -1);
        std::queue<int> pending;
        level_[static_cast<std::size_t>(source)] = 0;
        pending.push(source);
        while (!pending.empty()) {
            const int from = pending.front();
            pending.pop();
            for (const Edge& edge : graph_[static_cast<std::size_t>(from)]) {
                if (edge.capacity == 0 || level_[static_cast<std::size_t>(edge.to)] >= 0) {
                    continue;
                }
                level_[static_cast<std::size_t>(edge.to)] =
                    level_[static_cast<std::size_t>(from)] + 1;
                pending.push(edge.to);
            }
        }
        return level_[static_cast<std::size_t>(sink)] >= 0;
    }

    Count dfs(int from, int sink, Count limit) {
        if (from == sink) {
            return limit;
        }
        int& index = next_[static_cast<std::size_t>(from)];
        while (index < static_cast<int>(graph_[static_cast<std::size_t>(from)].size())) {
            Edge& edge = graph_[static_cast<std::size_t>(from)]
                               [static_cast<std::size_t>(index)];
            if (edge.capacity != 0
                && level_[static_cast<std::size_t>(edge.to)]
                       == level_[static_cast<std::size_t>(from)] + 1) {
                const Count pushed =
                    dfs(edge.to, sink, std::min(limit, edge.capacity));
                if (pushed != 0) {
                    edge.capacity -= pushed;
                    graph_[static_cast<std::size_t>(edge.to)]
                          [static_cast<std::size_t>(edge.reverse)]
                              .capacity += pushed;
                    return pushed;
                }
            }
            ++index;
        }
        return 0;
    }

    std::vector<std::vector<Edge>> graph_;
    std::vector<int> level_;
    std::vector<int> next_;
};

std::vector<Count> multiply(
    const std::vector<Count>& input,
    int label,
    int level
) {
    std::vector<Count> output(static_cast<std::size_t>(level + 1), 0);
    for (int left = 0; left <= level; ++left) {
        const Count value = input[static_cast<std::size_t>(left)];
        if (value == 0) {
            continue;
        }
        const int lower = std::abs(left - label);
        const int upper = std::min(left + label, 2 * level - left - label);
        for (int right = lower; right <= upper; right += 2) {
            output[static_cast<std::size_t>(right)] += value;
        }
    }
    return output;
}

std::string mask_string(int mask, int factors) {
    std::string answer = "{";
    bool first = true;
    for (int index = 0; index < factors; ++index) {
        if ((mask & (1 << index)) == 0) {
            continue;
        }
        if (!first) {
            answer += ",";
        }
        first = false;
        answer += std::to_string(index);
    }
    answer += "}";
    return answer;
}

enum class Routing {
    subsets,
    one_toggle,
    one_toggle_q1,
    one_import,
    two_import,
    three_import
};

bool routing_allowed(int source, int target, Routing routing) {
    if (routing == Routing::subsets) {
        return (target & ~source) == 0;
    }
    if (routing == Routing::one_import) {
        return __builtin_popcount(
                   static_cast<unsigned>(target & ~source)
               ) <= 1;
    }
    if (routing == Routing::two_import) {
        return __builtin_popcount(
                   static_cast<unsigned>(target & ~source)
               ) <= 2;
    }
    if (routing == Routing::three_import) {
        return __builtin_popcount(
                   static_cast<unsigned>(target & ~source)
               ) <= 3;
    }
    return __builtin_popcount(
               static_cast<unsigned>(source ^ target)
           ) <= 1;
}

bool inspect_word(
    const std::vector<int>& labels,
    int level,
    bool print_pass,
    Routing routing
) {
    const int factors = static_cast<int>(labels.size());
    const int masks = 1 << factors;
    std::vector<std::vector<Count>> product(
        static_cast<std::size_t>(masks),
        std::vector<Count>(static_cast<std::size_t>(level + 1), 0)
    );
    product[0][0] = 1;
    for (int mask = 1; mask < masks; ++mask) {
        const int bit = __builtin_ctz(static_cast<unsigned>(mask));
        product[static_cast<std::size_t>(mask)] = multiply(
            product[static_cast<std::size_t>(mask ^ (1 << bit))],
            labels[static_cast<std::size_t>(bit)],
            level
        );
    }

    const int maximum_q = routing == Routing::one_toggle_q1 ? 1 : level;
    for (int q = 1; q <= maximum_q; ++q) {
        for (int a = 0; a <= level; ++a) {
            std::vector<Count> demand(static_cast<std::size_t>(masks), 0);
            std::vector<Count> capacity(static_cast<std::size_t>(masks), 0);
            Count total_demand = 0;
            Count total_capacity = 0;
            for (int mask = 0; mask < masks; ++mask) {
                const int complement = (masks - 1) ^ mask;
                demand[static_cast<std::size_t>(mask)] =
                    product[static_cast<std::size_t>(mask)]
                           [static_cast<std::size_t>(q)]
                    * product[static_cast<std::size_t>(complement)]
                             [static_cast<std::size_t>(a)];
                std::vector<Count> fused = multiply(
                    product[static_cast<std::size_t>(complement)], q, level
                );
                capacity[static_cast<std::size_t>(mask)] =
                    product[static_cast<std::size_t>(mask)][0]
                    * fused[static_cast<std::size_t>(a)];
                total_demand += demand[static_cast<std::size_t>(mask)];
                total_capacity += capacity[static_cast<std::size_t>(mask)];
            }
            if (total_demand == 0 || total_demand > total_capacity) {
                continue;
            }

            const int source = 2 * masks;
            const int sink = source + 1;
            Dinic network(sink + 1);
            const Count infinity = total_demand;
            for (int left = 0; left < masks; ++left) {
                if (demand[static_cast<std::size_t>(left)] == 0) {
                    continue;
                }
                network.add_edge(source, left, demand[static_cast<std::size_t>(left)]);
                for (int right = 0; right < masks; ++right) {
                    if (routing_allowed(left, right, routing)
                        && capacity[static_cast<std::size_t>(right)] != 0) {
                        network.add_edge(left, masks + right, infinity);
                    }
                }
            }
            for (int right = 0; right < masks; ++right) {
                if (capacity[static_cast<std::size_t>(right)] != 0) {
                    network.add_edge(
                        masks + right, sink,
                        capacity[static_cast<std::size_t>(right)]
                    );
                }
            }
            const Count rank = network.flow(source, sink);
            if (rank == total_demand) {
                continue;
            }
            const auto reachable = network.reachable(source);
            Count witness_demand = 0;
            Count witness_capacity = 0;
            std::cout << "FAIL level=" << level << " labels=";
            for (int label : labels) {
                std::cout << label << ',';
            }
            std::cout << " q=" << q << " a=" << a
                      << " total_demand=" << total_demand
                      << " total_capacity=" << total_capacity
                      << " monotone_rank=" << rank << "\nleft=";
            for (int mask = 0; mask < masks; ++mask) {
                if (reachable[static_cast<std::size_t>(mask)]) {
                    witness_demand += demand[static_cast<std::size_t>(mask)];
                    std::cout << mask_string(mask, factors) << ':'
                              << demand[static_cast<std::size_t>(mask)] << ' ';
                }
            }
            std::cout << "\nneighbor=";
            for (int mask = 0; mask < masks; ++mask) {
                if (reachable[static_cast<std::size_t>(masks + mask)]) {
                    witness_capacity += capacity[static_cast<std::size_t>(mask)];
                    std::cout << mask_string(mask, factors) << ':'
                              << capacity[static_cast<std::size_t>(mask)] << ' ';
                }
            }
            std::cout << "\nwitness_demand=" << witness_demand
                      << " witness_capacity=" << witness_capacity << '\n';
            return false;
        }
    }
    if (print_pass) {
        std::cout << "PASS level=" << level << " labels=";
        for (int label : labels) {
            std::cout << label << ',';
        }
        std::cout << '\n';
    }
    return true;
}

bool enumerate_words(
    std::vector<int>& labels,
    int position,
    int maximum_label,
    int level,
    Routing routing
) {
    if (position == static_cast<int>(labels.size())) {
        return inspect_word(labels, level, false, routing);
    }
    const int lower = position == 0 ? 1 : labels[static_cast<std::size_t>(position - 1)];
    for (int label = lower; label <= maximum_label; ++label) {
        labels[static_cast<std::size_t>(position)] = label;
        if (!enumerate_words(
                labels, position + 1, maximum_label, level, routing
            )) {
            return false;
        }
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc >= 2
        && (std::string(argv[1]) == "--word"
            || std::string(argv[1]) == "--one-import-word"
            || std::string(argv[1]) == "--two-import-word"
            || std::string(argv[1]) == "--three-import-word")) {
        if (argc < 5) {
            std::cerr << "usage: --word LEVEL LABEL...\n"
                         "       --one-import-word LEVEL LABEL...\n";
            return 2;
        }
        Routing routing = Routing::subsets;
        if (std::string(argv[1]) == "--one-import-word") {
            routing = Routing::one_import;
        } else if (std::string(argv[1]) == "--two-import-word") {
            routing = Routing::two_import;
        } else if (std::string(argv[1]) == "--three-import-word") {
            routing = Routing::three_import;
        }
        const int level = std::stoi(argv[2]);
        std::vector<int> labels;
        for (int index = 3; index < argc; ++index) {
            labels.push_back(std::stoi(argv[index]));
        }
        return inspect_word(labels, level, true, routing) ? 0 : 1;
    }
    Routing routing = Routing::subsets;
    int argument = 1;
    if (argc > 1 && std::string(argv[1]) == "--one-toggle") {
        routing = Routing::one_toggle;
        ++argument;
    } else if (argc > 1 && std::string(argv[1]) == "--one-toggle-q1") {
        routing = Routing::one_toggle_q1;
        ++argument;
    } else if (argc > 1 && std::string(argv[1]) == "--one-import") {
        routing = Routing::one_import;
        ++argument;
    } else if (argc > 1 && std::string(argv[1]) == "--two-import") {
        routing = Routing::two_import;
        ++argument;
    } else if (argc > 1 && std::string(argv[1]) == "--three-import") {
        routing = Routing::three_import;
        ++argument;
    }
    const int maximum_factors =
        argc > argument ? std::stoi(argv[argument]) : 7;
    const int maximum_label =
        argc > argument + 1 ? std::stoi(argv[argument + 1]) : 5;
    const int level =
        argc > argument + 2 ? std::stoi(argv[argument + 2]) : 20;
    for (int factors = 1; factors <= maximum_factors; ++factors) {
        std::vector<int> labels(static_cast<std::size_t>(factors), 1);
        if (!enumerate_words(labels, 0, maximum_label, level, routing)) {
            return 1;
        }
        std::cout << "PASS factors=" << factors << '\n';
    }
    return 0;
}
