#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
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
        if (argc != 2 && argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_fusion_markov_monotonicity "
                "MAXIMUM_LEVEL [--raw]"
            );
        }
        const bool lazy = argc == 2;
        if (argc == 3 && std::string(argv[2]) != "--raw") {
            throw std::runtime_error(
                "the only optional argument is --raw"
            );
        }
        const std::string tag = lazy
            ? "SU2_FUSION_LAZY_MARKOV_MONOTONICITY"
            : "SU2_FUSION_MARKOV_MONOTONICITY";
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        long double minimum_margin =
            std::numeric_limits<long double>::infinity();
        int witness_level = 0;
        int witness_label = 0;
        int witness_source = 0;
        int witness_threshold = 0;
        std::uint64_t comparisons = 0U;

        for (int level = 6; level <= maximum_level; level += 2) {
            const long double angle =
                std::acos(-1.0L)
                / static_cast<long double>(level + 2);
            std::vector<long double> dimension(
                static_cast<std::size_t>(level + 1)
            );
            for (int vertex = 0; vertex <= level; ++vertex) {
                dimension[static_cast<std::size_t>(vertex)] =
                    std::sin(
                        static_cast<long double>(vertex + 1) * angle
                    ) / std::sin(angle);
            }
            for (int label = 2;
                 2 * label < level;
                 label += 2) {
                for (int source = 0;
                     source + 2 <= level;
                     source += 2) {
                    std::vector<long double> lower_row(
                        static_cast<std::size_t>(level + 1),
                        0.0L
                    );
                    std::vector<long double> upper_row(
                        static_cast<std::size_t>(level + 1),
                        0.0L
                    );
                    long double lower_total = 0.0L;
                    long double upper_total = 0.0L;
                    for (int target = 0;
                         target <= level;
                         target += 2) {
                        if (fuses(level, label, source, target)) {
                            lower_row[
                                static_cast<std::size_t>(target)
                            ] = dimension[
                                static_cast<std::size_t>(target)
                            ];
                            lower_total += dimension[
                                static_cast<std::size_t>(target)
                            ];
                        }
                        if (fuses(
                                level,
                                label,
                                source + 2,
                                target
                            )) {
                            upper_row[
                                static_cast<std::size_t>(target)
                            ] = dimension[
                                static_cast<std::size_t>(target)
                            ];
                            upper_total += dimension[
                                static_cast<std::size_t>(target)
                            ];
                        }
                    }
                    long double lower_cdf = 0.0L;
                    long double upper_cdf = 0.0L;
                    for (int threshold = 0;
                         threshold <= level;
                         threshold += 2) {
                        const long double move_weight =
                            lazy ? 0.5L : 1.0L;
                        lower_cdf += move_weight * lower_row[
                            static_cast<std::size_t>(threshold)
                        ] / lower_total;
                        upper_cdf += move_weight * upper_row[
                            static_cast<std::size_t>(threshold)
                        ] / upper_total;
                        if (lazy && threshold == source) {
                            lower_cdf += 0.5L;
                        }
                        if (lazy && threshold == source + 2) {
                            upper_cdf += 0.5L;
                        }
                        const long double margin =
                            lower_cdf - upper_cdf;
                        ++comparisons;
                        if (margin < minimum_margin) {
                            minimum_margin = margin;
                            witness_level = level;
                            witness_label = label;
                            witness_source = source;
                            witness_threshold = threshold;
                        }
                        if (margin < -1.0e-15L) {
                            std::cout
                                << tag
                                << " counterexample"
                                << " level=" << level
                                << " label=" << label
                                << " source=" << source
                                << " threshold=" << threshold
                                << " margin="
                                << std::setprecision(20)
                                << margin
                                << " result=FAIL\n";
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }
        std::cout
            << tag
            << " maximum_level=" << maximum_level
            << " comparisons=" << comparisons
            << " minimum_margin=" << std::setprecision(20)
            << minimum_margin
            << " witness=("
            << witness_level << ','
            << witness_label << ','
            << witness_source << ','
            << witness_threshold << ')'
            << " result=PASS_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_FUSION_MARKOV_MONOTONICITY FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
