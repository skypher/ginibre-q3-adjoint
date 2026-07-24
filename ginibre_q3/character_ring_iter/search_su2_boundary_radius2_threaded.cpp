#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mutex>
#include <queue>
#include <random>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

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
        Edge forward{
            to,
            static_cast<int>(graph_[static_cast<std::size_t>(to)].size()),
            capacity
        };
        Edge reverse{
            from,
            static_cast<int>(graph_[static_cast<std::size_t>(from)].size()),
            0
        };
        graph_[static_cast<std::size_t>(from)].push_back(forward);
        graph_[static_cast<std::size_t>(to)].push_back(reverse);
    }

    Count flow(int source, int sink) {
        Count answer = 0;
        while (bfs(source, sink)) {
            std::fill(next_.begin(), next_.end(), 0);
            while (const Count pushed =
                       dfs(source, sink, std::numeric_limits<Count>::max())) {
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

struct Failure {
    int level = 0;
    int q = 0;
    int a = 0;
    Count demand = 0;
    Count capacity = 0;
    Count rank = 0;
    std::vector<int> labels;
    std::vector<std::pair<int, Count>> left;
    std::vector<std::pair<int, Count>> neighbor;
};

bool inspect(
    const std::vector<int>& labels,
    int level,
    int q,
    int a,
    Failure& failure
) {
    const int factors = static_cast<int>(labels.size());
    const int masks = 1 << factors;
    const int full = masks - 1;
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

    std::vector<Count> demand(static_cast<std::size_t>(masks), 0);
    std::vector<Count> capacity(static_cast<std::size_t>(masks), 0);
    Count total_demand = 0;
    Count total_capacity = 0;
    for (int mask = 0; mask < masks; ++mask) {
        const int complement = full ^ mask;
        demand[static_cast<std::size_t>(mask)] =
            product[static_cast<std::size_t>(mask)][static_cast<std::size_t>(q)]
            * product[static_cast<std::size_t>(complement)]
                     [static_cast<std::size_t>(a)];
        const auto fused = multiply(
            product[static_cast<std::size_t>(complement)], q, level
        );
        capacity[static_cast<std::size_t>(mask)] =
            product[static_cast<std::size_t>(mask)][0]
            * fused[static_cast<std::size_t>(a)];
        total_demand += demand[static_cast<std::size_t>(mask)];
        total_capacity += capacity[static_cast<std::size_t>(mask)];
    }
    if (total_demand == 0) {
        return true;
    }

    // Nodes: sources, deletion layer, one-import layer, two-import layer,
    // targets, global source, global sink.
    const int source_offset = 0;
    const int deletion_offset = masks;
    const int one_offset = 2 * masks;
    const int two_offset = 3 * masks;
    const int target_offset = 4 * masks;
    const int source = 5 * masks;
    const int sink = source + 1;
    Dinic network(sink + 1);
    const Count infinity = total_demand;

    for (int mask = 0; mask < masks; ++mask) {
        if (demand[static_cast<std::size_t>(mask)] != 0) {
            network.add_edge(
                source, source_offset + mask,
                demand[static_cast<std::size_t>(mask)]
            );
            network.add_edge(source_offset + mask, deletion_offset + mask, infinity);
        }
        if (capacity[static_cast<std::size_t>(mask)] != 0) {
            network.add_edge(
                target_offset + mask, sink,
                capacity[static_cast<std::size_t>(mask)]
            );
        }
        network.add_edge(deletion_offset + mask, target_offset + mask, infinity);
        network.add_edge(one_offset + mask, target_offset + mask, infinity);
        network.add_edge(two_offset + mask, target_offset + mask, infinity);
        for (int bit = 0; bit < factors; ++bit) {
            const int flag = 1 << bit;
            if ((mask & flag) != 0) {
                network.add_edge(
                    deletion_offset + mask,
                    deletion_offset + (mask ^ flag),
                    infinity
                );
            } else {
                network.add_edge(
                    deletion_offset + mask, one_offset + (mask | flag), infinity
                );
                network.add_edge(
                    one_offset + mask, two_offset + (mask | flag), infinity
                );
            }
        }
    }

    const Count rank = network.flow(source, sink);
    if (rank == total_demand) {
        return true;
    }
    failure.level = level;
    failure.q = q;
    failure.a = a;
    failure.demand = total_demand;
    failure.capacity = total_capacity;
    failure.rank = rank;
    failure.labels = labels;
    const auto reachable = network.reachable(source);
    for (int mask = 0; mask < masks; ++mask) {
        if (reachable[static_cast<std::size_t>(source_offset + mask)]) {
            failure.left.emplace_back(
                mask, demand[static_cast<std::size_t>(mask)]
            );
        }
        if (reachable[static_cast<std::size_t>(target_offset + mask)]) {
            failure.neighbor.emplace_back(
                mask, capacity[static_cast<std::size_t>(mask)]
            );
        }
    }
    return false;
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

void print_failure(const Failure& failure) {
    std::cout << "FAIL level=" << failure.level << " labels=";
    for (int label : failure.labels) {
        std::cout << label << ',';
    }
    std::cout << " q=" << failure.q << " a=" << failure.a
              << " total_demand=" << failure.demand
              << " total_capacity=" << failure.capacity
              << " radius2_rank=" << failure.rank << "\nleft=";
    for (const auto& [mask, value] : failure.left) {
        std::cout << mask_string(mask, static_cast<int>(failure.labels.size()))
                  << ':' << value << ' ';
    }
    std::cout << "\nneighbor=";
    for (const auto& [mask, value] : failure.neighbor) {
        std::cout << mask_string(mask, static_cast<int>(failure.labels.size()))
                  << ':' << value << ' ';
    }
    std::cout << '\n';
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr
            << "usage: LEVEL MAXIMUM_FACTORS TRIALS SEED [TOP_WINDOW]\n";
        return 2;
    }
    const int level = std::stoi(argv[1]);
    const int maximum_factors = std::stoi(argv[2]);
    const std::uint64_t trials = std::stoull(argv[3]);
    const std::uint64_t seed = std::stoull(argv[4]);
    const int top_window = argc > 5 ? std::stoi(argv[5]) : 5;
    if (level < 2 || maximum_factors < 2 || maximum_factors > 15
        || top_window < 1) {
        std::cerr << "invalid bounds\n";
        return 2;
    }

    std::atomic<bool> failed{false};
    std::atomic<std::uint64_t> completed{0};
    std::mutex failure_mutex;
    Failure first_failure;

#pragma omp parallel
    {
#ifdef _OPENMP
        const std::uint64_t thread =
            static_cast<std::uint64_t>(omp_get_thread_num());
#else
        const std::uint64_t thread = 0;
#endif
        std::mt19937_64 generator(
            seed + UINT64_C(0x9e3779b97f4a7c15) * (thread + 1)
        );
        std::uniform_int_distribution<int> factor_distribution(2, maximum_factors);
        std::uniform_int_distribution<int> label_distribution(1, level);
        std::uniform_int_distribution<int> coin(0, 1);

#pragma omp for schedule(dynamic)
        for (std::uint64_t trial = 0; trial < trials; ++trial) {
            if (failed.load(std::memory_order_relaxed)) {
                continue;
            }
            const int factors = factor_distribution(generator);
            std::vector<int> labels(static_cast<std::size_t>(factors));
            for (int& label : labels) {
                if (coin(generator) == 0) {
                    label = std::min(level, label_distribution(generator));
                } else {
                    const int width = std::min(top_window, level);
                    std::uniform_int_distribution<int> shell(0, width - 1);
                    label = coin(generator) == 0
                                ? 1 + shell(generator)
                                : level - shell(generator);
                }
            }
            std::sort(labels.begin(), labels.end());

            const int lower = std::max(1, level - top_window);
            for (int q = lower; q <= level && !failed.load(); ++q) {
                for (int a = lower; a <= level && !failed.load(); ++a) {
                    Failure witness;
                    if (!inspect(labels, level, q, a, witness)) {
                        bool expected = false;
                        if (failed.compare_exchange_strong(expected, true)) {
                            std::lock_guard<std::mutex> lock(failure_mutex);
                            first_failure = std::move(witness);
                        }
                    }
                }
            }
            ++completed;
        }
    }

    if (failed.load()) {
        print_failure(first_failure);
        std::cout << "completed_before_failure=" << completed.load() << '\n';
        return 1;
    }
#ifdef _OPENMP
    const int threads = omp_get_max_threads();
#else
    const int threads = 1;
#endif
    std::cout << "PASS level=" << level << " maximum_factors="
              << maximum_factors << " trials=" << completed.load()
              << " threads=" << threads << " top_window=" << top_window
              << '\n';
    return 0;
}
