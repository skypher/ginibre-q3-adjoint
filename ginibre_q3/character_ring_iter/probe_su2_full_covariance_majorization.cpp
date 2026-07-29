#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Pair = std::array<int, 2>;
using Quartet = std::array<int, 4>;
using Polynomial2 = std::map<Pair, long long>;
using Polynomial4 = std::map<Quartet, long long>;

int parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0;
    const long long parsed = std::stoll(value, &consumed);
    if (consumed != value.size() || parsed <= 0) {
        throw std::invalid_argument(name + " must be positive");
    }
    return static_cast<int>(parsed);
}

Polynomial2 coefficient(int support, int label) {
    Polynomial2 result;
    for (int left = 0; left <= support; ++left) {
        for (int right = 0; right <= support; ++right) {
            if (
                std::abs(left - right) <= label
                && label <= left + right
            ) {
                Pair indices{left, right};
                std::sort(indices.begin(), indices.end());
                ++result[indices];
            }
        }
    }
    return result;
}

Polynomial2 fusion_sum(int support, int q, int a) {
    Polynomial2 result;
    for (int label = std::abs(q - a); label <= q + a; ++label) {
        for (const auto& [indices, multiplicity] :
             coefficient(support, label)) {
            result[indices] += multiplicity;
        }
    }
    return result;
}

void add_product(
    Polynomial4& result,
    const Polynomial2& left,
    const Polynomial2& right,
    long long sign
) {
    for (const auto& [first, first_coefficient] : left) {
        for (const auto& [second, second_coefficient] : right) {
            Quartet indices{
                first[0], first[1], second[0], second[1]
            };
            std::sort(indices.begin(), indices.end());
            result[indices] +=
                sign * first_coefficient * second_coefficient;
        }
    }
}

bool majorizes(const Quartet& outer, const Quartet& inner) {
    int outer_total = 0;
    int inner_total = 0;
    for (int index = 0; index < 4; ++index) {
        outer_total += outer[index];
        inner_total += inner[index];
    }
    if (outer_total != inner_total) {
        return false;
    }
    int outer_suffix = 0;
    int inner_suffix = 0;
    for (int index = 3; index >= 1; --index) {
        outer_suffix += outer[index];
        inner_suffix += inner[index];
        if (outer_suffix < inner_suffix) {
            return false;
        }
    }
    return true;
}

struct Edge {
    int target;
    long long capacity;
    int reverse;
};

class Dinic {
public:
    explicit Dinic(int vertices)
        : graph_(static_cast<std::size_t>(vertices)),
          level_(static_cast<std::size_t>(vertices)),
          next_(static_cast<std::size_t>(vertices)) {}

    void add_edge(int source, int target, long long capacity) {
        const int source_reverse =
            static_cast<int>(graph_[static_cast<std::size_t>(target)].size());
        const int target_reverse =
            static_cast<int>(graph_[static_cast<std::size_t>(source)].size());
        graph_[static_cast<std::size_t>(source)].push_back(
            {target, capacity, source_reverse}
        );
        graph_[static_cast<std::size_t>(target)].push_back(
            {source, 0, target_reverse}
        );
    }

    long long maximum_flow(int source, int sink) {
        long long result = 0;
        while (build_levels(source, sink)) {
            std::fill(next_.begin(), next_.end(), 0);
            while (true) {
                const long long pushed = push(
                    source,
                    sink,
                    std::numeric_limits<long long>::max()
                );
                if (pushed == 0) {
                    break;
                }
                result += pushed;
            }
        }
        return result;
    }

private:
    bool build_levels(int source, int sink) {
        std::fill(level_.begin(), level_.end(), -1);
        std::queue<int> queue;
        level_[static_cast<std::size_t>(source)] = 0;
        queue.push(source);
        while (!queue.empty()) {
            const int vertex = queue.front();
            queue.pop();
            for (const Edge& edge :
                 graph_[static_cast<std::size_t>(vertex)]) {
                if (
                    edge.capacity > 0
                    && level_[static_cast<std::size_t>(edge.target)] < 0
                ) {
                    level_[static_cast<std::size_t>(edge.target)] =
                        level_[static_cast<std::size_t>(vertex)] + 1;
                    queue.push(edge.target);
                }
            }
        }
        return level_[static_cast<std::size_t>(sink)] >= 0;
    }

    long long push(int vertex, int sink, long long limit) {
        if (vertex == sink) {
            return limit;
        }
        int& begin = next_[static_cast<std::size_t>(vertex)];
        auto& edges = graph_[static_cast<std::size_t>(vertex)];
        for (; begin < static_cast<int>(edges.size()); ++begin) {
            Edge& edge = edges[static_cast<std::size_t>(begin)];
            if (
                edge.capacity == 0
                || level_[static_cast<std::size_t>(edge.target)]
                    != level_[static_cast<std::size_t>(vertex)] + 1
            ) {
                continue;
            }
            const long long pushed =
                push(edge.target, sink, std::min(limit, edge.capacity));
            if (pushed == 0) {
                continue;
            }
            edge.capacity -= pushed;
            graph_[static_cast<std::size_t>(edge.target)]
                [static_cast<std::size_t>(edge.reverse)].capacity += pushed;
            return pushed;
        }
        return 0;
    }

    std::vector<std::vector<Edge>> graph_;
    std::vector<int> level_;
    std::vector<int> next_;
};

std::string show(const Quartet& indices) {
    return
        "{" + std::to_string(indices[0])
        + "," + std::to_string(indices[1])
        + "," + std::to_string(indices[2])
        + "," + std::to_string(indices[3]) + "}";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        int maximum_support = 10;
        int maximum_label = 10;
        if (argc >= 2) {
            maximum_support = parse_positive(argv[1], "maximum_support");
        }
        if (argc >= 3) {
            maximum_label = parse_positive(argv[2], "maximum_label");
        }
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_full_covariance_majorization"
                " [maximum_support] [maximum_label]"
            );
        }

        unsigned long long cases = 0;
        unsigned long long negative_units = 0;
        for (int support = 1; support <= maximum_support; ++support) {
            const Polynomial2 norm = coefficient(support, 0);
            for (int q = 1; q <= maximum_label; ++q) {
                const Polynomial2 q_coefficient =
                    coefficient(support, q);
                for (int a = q + 1; a <= maximum_label; ++a) {
                    const Polynomial2 a_coefficient =
                        coefficient(support, a);
                    Polynomial4 determinant;
                    add_product(
                        determinant,
                        norm,
                        fusion_sum(support, q, a),
                        1
                    );
                    add_product(
                        determinant,
                        q_coefficient,
                        a_coefficient,
                        -1
                    );

                    std::vector<std::pair<Quartet, long long>> negative;
                    std::vector<std::pair<Quartet, long long>> positive;
                    long long demand = 0;
                    for (const auto& [indices, coefficient_value] :
                         determinant) {
                        if (coefficient_value < 0) {
                            negative.push_back(
                                {indices, -coefficient_value}
                            );
                            demand -= coefficient_value;
                        } else if (coefficient_value > 0) {
                            positive.push_back(
                                {indices, coefficient_value}
                            );
                        }
                    }

                    const int source = 0;
                    const int negative_start = 1;
                    const int positive_start =
                        negative_start + static_cast<int>(negative.size());
                    const int sink =
                        positive_start + static_cast<int>(positive.size());
                    Dinic flow(sink + 1);
                    for (int index = 0;
                         index < static_cast<int>(negative.size());
                         ++index) {
                        flow.add_edge(
                            source,
                            negative_start + index,
                            negative[static_cast<std::size_t>(index)].second
                        );
                    }
                    for (int left = 0;
                         left < static_cast<int>(negative.size());
                         ++left) {
                        for (int right = 0;
                             right < static_cast<int>(positive.size());
                             ++right) {
                            if (
                                majorizes(
                                    negative[static_cast<std::size_t>(left)]
                                        .first,
                                    positive[static_cast<std::size_t>(right)]
                                        .first
                                )
                            ) {
                                flow.add_edge(
                                    negative_start + left,
                                    positive_start + right,
                                    demand
                                );
                            }
                        }
                    }
                    for (int index = 0;
                         index < static_cast<int>(positive.size());
                         ++index) {
                        flow.add_edge(
                            positive_start + index,
                            sink,
                            positive[static_cast<std::size_t>(index)].second
                        );
                    }
                    const long long paid = flow.maximum_flow(source, sink);
                    ++cases;
                    negative_units += static_cast<unsigned long long>(demand);
                    if (paid != demand) {
                        std::cout
                            << "SU2_FULL_COVARIANCE_MAJORIZATION"
                            << " support=" << support
                            << " q=" << q
                            << " a=" << a
                            << " negative_monomials=" << negative.size()
                            << " positive_monomials=" << positive.size()
                            << " demand=" << demand
                            << " paid=" << paid;
                        for (const auto& [indices, coefficient_value] :
                             negative) {
                            bool has_target = false;
                            for (const auto& [target, capacity] : positive) {
                                if (
                                    capacity > 0
                                    && majorizes(indices, target)
                                ) {
                                    has_target = true;
                                    break;
                                }
                            }
                            if (!has_target) {
                                std::cout
                                    << " isolated_negative="
                                    << show(indices)
                                    << ":" << coefficient_value;
                                break;
                            }
                        }
                        std::cout << " polynomial={";
                        bool first_term = true;
                        for (const auto& [indices, coefficient_value] :
                             determinant) {
                            if (coefficient_value == 0) {
                                continue;
                            }
                            if (!first_term) {
                                std::cout << ',';
                            }
                            first_term = false;
                            std::cout
                                << show(indices)
                                << ':' << coefficient_value;
                        }
                        std::cout << '}';
                        std::cout << " result=COUNTEREXAMPLE_TO_ROUTE\n";
                        return EXIT_SUCCESS;
                    }
                }
            }
        }
        std::cout
            << "SU2_FULL_COVARIANCE_MAJORIZATION"
            << " maximum_support=" << maximum_support
            << " maximum_label=" << maximum_label
            << " cases=" << cases
            << " negative_units=" << negative_units
            << " result=ALL_MATCHED_BOUNDED\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
