#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

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

struct CaseResult {
    bool balanced = true;
    bool boundary_compatible = true;
    int pair_count = 0;
    std::uint64_t nonzero_compound_entries = 0U;
    int conflict_first = -1;
    int conflict_second = -1;
    int bad_boundary_vertex = -1;
};

CaseResult analyze_case(int level, int label, bool print_detail) {
    std::vector<int> vertices;
    for (int vertex = 0; vertex <= level; vertex += 2) {
        vertices.push_back(vertex);
    }
    const int rank = static_cast<int>(vertices.size());
    std::vector<std::vector<long long>> fusion(
        static_cast<std::size_t>(rank),
        std::vector<long long>(static_cast<std::size_t>(rank), 0)
    );
    for (int row = 0; row < rank; ++row) {
        for (int column = 0; column < rank; ++column) {
            fusion[static_cast<std::size_t>(row)][
                static_cast<std::size_t>(column)
            ] = fuses(
                level,
                label,
                vertices[static_cast<std::size_t>(row)],
                vertices[static_cast<std::size_t>(column)]
            ) ? 1 : 0;
        }
    }
    std::vector<std::vector<long long>> square(
        static_cast<std::size_t>(rank),
        std::vector<long long>(static_cast<std::size_t>(rank), 0)
    );
    for (int row = 0; row < rank; ++row) {
        for (int middle = 0; middle < rank; ++middle) {
            if (fusion[static_cast<std::size_t>(row)][
                    static_cast<std::size_t>(middle)] == 0) {
                continue;
            }
            for (int column = 0; column < rank; ++column) {
                square[static_cast<std::size_t>(row)][
                    static_cast<std::size_t>(column)
                ] += fusion[static_cast<std::size_t>(middle)][
                    static_cast<std::size_t>(column)
                ];
            }
        }
    }
    std::vector<std::pair<int, int>> pairs;
    std::vector<std::vector<int>> pair_id(
        static_cast<std::size_t>(rank),
        std::vector<int>(static_cast<std::size_t>(rank), -1)
    );
    for (int first = 0; first < rank; ++first) {
        for (int second = first + 1; second < rank; ++second) {
            pair_id[static_cast<std::size_t>(first)][
                static_cast<std::size_t>(second)
            ] = static_cast<int>(pairs.size());
            pairs.emplace_back(first, second);
        }
    }
    const int pair_count = static_cast<int>(pairs.size());
    std::vector<std::vector<std::pair<int, int>>> signed_adjacency(
        static_cast<std::size_t>(pair_count)
    );
    CaseResult result;
    result.pair_count = pair_count;
    for (int left = 0; left < pair_count; ++left) {
        const auto [a, b] = pairs[static_cast<std::size_t>(left)];
        for (int right = 0; right < pair_count; ++right) {
            const auto [c, d] = pairs[static_cast<std::size_t>(right)];
            const long long minor =
                square[static_cast<std::size_t>(a)][
                    static_cast<std::size_t>(c)
                ] * square[static_cast<std::size_t>(b)][
                    static_cast<std::size_t>(d)
                ]
                - square[static_cast<std::size_t>(a)][
                    static_cast<std::size_t>(d)
                ] * square[static_cast<std::size_t>(b)][
                    static_cast<std::size_t>(c)
                ];
            if (minor == 0) {
                continue;
            }
            ++result.nonzero_compound_entries;
            signed_adjacency[static_cast<std::size_t>(left)].push_back(
                {right, minor > 0 ? 1 : -1}
            );
        }
    }

    const int zero_index = 0;
    const int q_index = label / 2;
    const int start = pair_id[static_cast<std::size_t>(zero_index)][
        static_cast<std::size_t>(q_index)
    ];
    std::vector<int> sign(static_cast<std::size_t>(pair_count), 0);
    sign[static_cast<std::size_t>(start)] = 1;
    std::queue<int> queue;
    queue.push(start);
    while (!queue.empty() && result.balanced) {
        const int current = queue.front();
        queue.pop();
        for (const auto& [target, edge_sign] :
             signed_adjacency[static_cast<std::size_t>(current)]) {
            const int required =
                sign[static_cast<std::size_t>(current)] * edge_sign;
            int& target_sign = sign[static_cast<std::size_t>(target)];
            if (target_sign == 0) {
                target_sign = required;
                queue.push(target);
            } else if (target_sign != required) {
                result.balanced = false;
                result.conflict_first = current;
                result.conflict_second = target;
                break;
            }
        }
    }
    if (result.balanced) {
        for (int vertex = 1; vertex < rank; ++vertex) {
            const int boundary = pair_id[0][
                static_cast<std::size_t>(vertex)
            ];
            if (sign[static_cast<std::size_t>(boundary)] != 0
                && sign[static_cast<std::size_t>(boundary)] != 1) {
                result.boundary_compatible = false;
                result.bad_boundary_vertex =
                    vertices[static_cast<std::size_t>(vertex)];
                break;
            }
        }
    }
    if (print_detail) {
        std::cout
            << "SU2_M2_COMPOUND_SIGN_CASE"
            << " level=" << level
            << " label=" << label
            << " rank=" << rank
            << " pair_count=" << pair_count
            << " nonzero_entries=" << result.nonzero_compound_entries
            << " balanced=" << (result.balanced ? 1 : 0)
            << " boundary_compatible="
                << (result.boundary_compatible ? 1 : 0)
            << " conflict=(" << result.conflict_first
            << ',' << result.conflict_second << ')'
            << " bad_boundary_vertex=" << result.bad_boundary_vertex
            << '\n';
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4 && std::string(argv[1]) == "--case") {
            const int level = parse_positive(argv[2], "level");
            const int label = parse_positive(argv[3], "label");
            if ((level & 1) != 0 || (label & 1) != 0
                || 2 * label >= level) {
                throw std::runtime_error(
                    "case requires even level/label and 2*label<level"
                );
            }
            analyze_case(level, label, true);
            return EXIT_SUCCESS;
        }
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_m2_compound_sign MAXIMUM_LEVEL "
                "| --case LEVEL LABEL"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        std::uint64_t cases = 0U;
        std::uint64_t balanced = 0U;
        std::uint64_t compatible = 0U;
        bool printed_failure = false;
        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                ++cases;
                const CaseResult result =
                    analyze_case(level, label, false);
                if (result.balanced) {
                    ++balanced;
                }
                if (result.balanced && result.boundary_compatible) {
                    ++compatible;
                }
                if ((!result.balanced || !result.boundary_compatible)
                    && !printed_failure) {
                    analyze_case(level, label, true);
                    printed_failure = true;
                }
            }
        }
        std::cout
            << "SU2_M2_COMPOUND_SIGN"
            << " maximum_level=" << maximum_level
            << " cases=" << cases
            << " balanced=" << balanced
            << " compatible=" << compatible
            << " result="
            << (compatible == cases ? "PASS_DISCOVERY" : "FAIL_SIGNING")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_M2_COMPOUND_SIGN FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
