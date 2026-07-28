#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
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

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_fusion_chordality MAXIMUM_LEVEL"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        std::uint64_t parameters = 0U;
        std::uint64_t eliminated_vertices = 0U;
        std::uint64_t clique_tests = 0U;

        for (int level = 6; level <= maximum_level; level += 2) {
            const int vertex_count = level / 2 + 1;
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameters;
                std::vector<int> order;
                order.reserve(static_cast<std::size_t>(vertex_count));
                int lower = 0;
                int upper = level;
                while (lower <= upper) {
                    order.push_back(lower);
                    if (lower != upper) {
                        order.push_back(upper);
                    }
                    lower += 2;
                    upper -= 2;
                }
                std::vector<bool> removed(
                    static_cast<std::size_t>(level + 1),
                    false
                );
                for (const int vertex : order) {
                    std::vector<int> remaining_neighbors;
                    for (int candidate = 0;
                         candidate <= level;
                         candidate += 2) {
                        if (candidate != vertex
                            && !removed[
                                static_cast<std::size_t>(candidate)
                            ]
                            && fuses(
                                level,
                                label,
                                vertex,
                                candidate
                            )) {
                            remaining_neighbors.push_back(candidate);
                        }
                    }
                    for (std::size_t first = 0U;
                         first < remaining_neighbors.size();
                         ++first) {
                        for (std::size_t second = first + 1U;
                             second < remaining_neighbors.size();
                             ++second) {
                            ++clique_tests;
                            if (!fuses(
                                    level,
                                    label,
                                    remaining_neighbors[first],
                                    remaining_neighbors[second]
                                )) {
                                std::cout
                                    << "SU2_FUSION_CHORDALITY"
                                    << " first_failure"
                                    << " level=" << level
                                    << " label=" << label
                                    << " eliminated=" << vertex
                                    << " neighbors=("
                                    << remaining_neighbors[first] << ','
                                    << remaining_neighbors[second] << ')'
                                    << " result=FAIL_OUTSIDE_IN_PEO\n";
                                return EXIT_FAILURE;
                            }
                        }
                    }
                    removed[static_cast<std::size_t>(vertex)] = true;
                    ++eliminated_vertices;
                }
            }
        }

        std::cout
            << "SU2_FUSION_CHORDALITY"
            << " maximum_level=" << maximum_level
            << " parameters=" << parameters
            << " eliminated_vertices=" << eliminated_vertices
            << " clique_tests=" << clique_tests
            << " result=PASS_OUTSIDE_IN_PEO_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_FUSION_CHORDALITY FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
