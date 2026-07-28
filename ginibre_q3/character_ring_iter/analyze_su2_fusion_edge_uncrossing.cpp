#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

int parse_positive(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error("maximum level must be positive");
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
                "usage: analyze_su2_fusion_edge_uncrossing MAXIMUM_LEVEL"
            );
        }
        const int maximum_level = parse_positive(argv[1]);
        std::uint64_t crossed_pairs = 0;
        std::uint64_t failures = 0;
        std::uint64_t endpoint_failures = 0;
        std::uint64_t interior_failures = 0;
        std::uint64_t lower_wall_failures = 0;
        std::uint64_t upper_wall_failures = 0;
        std::uint64_t two_wall_failures = 0;
        std::uint64_t unclassified_failures = 0;
        std::uint64_t reflection_failures = 0;
        bool printed_first = false;
        bool printed_first_interior = false;

        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                for (int lower_source = 0;
                     lower_source < level;
                     lower_source += 2) {
                    for (int upper_source = lower_source + 2;
                         upper_source <= level;
                         upper_source += 2) {
                        for (int lower_target = 0;
                             lower_target < level;
                             lower_target += 2) {
                            for (int upper_target = lower_target + 2;
                                 upper_target <= level;
                                 upper_target += 2) {
                                if (!fuses(
                                        level,
                                        label,
                                        lower_source,
                                        upper_target
                                    )
                                    || !fuses(
                                        level,
                                        label,
                                        upper_source,
                                        lower_target
                                    )) {
                                    continue;
                                }
                                ++crossed_pairs;
                                const bool lower_edge = fuses(
                                    level,
                                    label,
                                    lower_source,
                                    lower_target
                                );
                                const bool upper_edge = fuses(
                                    level,
                                    label,
                                    upper_source,
                                    upper_target
                                );
                                if (lower_edge && upper_edge) {
                                    continue;
                                }
                                ++failures;
                                const bool touches_endpoint =
                                    lower_source == 0
                                    || upper_source == level
                                    || lower_target == 0
                                    || upper_target == level;
                                if (touches_endpoint) {
                                    ++endpoint_failures;
                                } else {
                                    ++interior_failures;
                                }
                                const bool lower_wall =
                                    lower_target
                                    < std::abs(lower_source - label);
                                const bool upper_wall =
                                    upper_target
                                    > 2 * level - upper_source - label;
                                if (lower_wall) {
                                    ++lower_wall_failures;
                                }
                                if (upper_wall) {
                                    ++upper_wall_failures;
                                }
                                if (lower_wall && upper_wall) {
                                    ++two_wall_failures;
                                }
                                if (!lower_wall && !upper_wall) {
                                    ++unclassified_failures;
                                }
                                const int reflected_lower_source =
                                    level - upper_target;
                                const int reflected_upper_source =
                                    level - lower_target;
                                const int reflected_lower_target =
                                    level - upper_source;
                                const int reflected_upper_target =
                                    level - lower_source;
                                const bool reflected_crossed =
                                    fuses(
                                        level,
                                        label,
                                        reflected_lower_source,
                                        reflected_upper_target
                                    )
                                    && fuses(
                                        level,
                                        label,
                                        reflected_upper_source,
                                        reflected_lower_target
                                    );
                                const bool reflected_lower_wall =
                                    reflected_lower_target
                                    < std::abs(
                                        reflected_lower_source - label
                                    );
                                const bool reflected_upper_wall =
                                    reflected_upper_target
                                    > 2 * level
                                        - reflected_upper_source
                                        - label;
                                if (!reflected_crossed
                                    || reflected_lower_wall != upper_wall
                                    || reflected_upper_wall != lower_wall) {
                                    ++reflection_failures;
                                }
                                if (!printed_first) {
                                    std::cout
                                        << "SU2_FUSION_EDGE_UNCROSSING"
                                        << " first_failure"
                                        << " level=" << level
                                        << " label=" << label
                                        << " sources=(" << lower_source
                                        << ',' << upper_source << ')'
                                        << " targets=(" << lower_target
                                        << ',' << upper_target << ')'
                                        << " lower_edge="
                                        << (lower_edge ? 1 : 0)
                                        << " upper_edge="
                                        << (upper_edge ? 1 : 0)
                                        << '\n';
                                    printed_first = true;
                                }
                                if (!touches_endpoint
                                    && !printed_first_interior) {
                                    std::cout
                                        << "SU2_FUSION_EDGE_UNCROSSING"
                                        << " first_interior_failure"
                                        << " level=" << level
                                        << " label=" << label
                                        << " sources=(" << lower_source
                                        << ',' << upper_source << ')'
                                        << " targets=(" << lower_target
                                        << ',' << upper_target << ')'
                                        << " lower_edge="
                                        << (lower_edge ? 1 : 0)
                                        << " upper_edge="
                                        << (upper_edge ? 1 : 0)
                                        << '\n';
                                    printed_first_interior = true;
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_FUSION_EDGE_UNCROSSING"
            << " maximum_level=" << maximum_level
            << " crossed_pairs=" << crossed_pairs
            << " failures=" << failures
            << " endpoint_failures=" << endpoint_failures
            << " interior_failures=" << interior_failures
            << " lower_wall_failures=" << lower_wall_failures
            << " upper_wall_failures=" << upper_wall_failures
            << " two_wall_failures=" << two_wall_failures
            << " unclassified_failures=" << unclassified_failures
            << " reflection_failures=" << reflection_failures
            << " result="
            << (
                unclassified_failures == 0
                    && reflection_failures == 0
                    ? "PASS_CLASSIFICATION"
                    : "FAIL_CLASSIFICATION"
            )
            << '\n';
        return unclassified_failures == 0
                && reflection_failures == 0
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_FUSION_EDGE_UNCROSSING FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
