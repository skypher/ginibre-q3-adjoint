#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
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

long double character_value(int label, long double theta) {
    return std::sin(
        static_cast<long double>(label + 1) * theta
    ) / std::sin(theta);
}

struct SpectralPair {
    long double base;
    long double coefficient;
};

struct CaseResult {
    bool failed = false;
    bool tail_failed = false;
    bool shifted_tail_failed = false;
    long double minimum_grouped_residue = 0.0L;
    long double minimum_tail = 0.0L;
    long double minimum_shifted_tail = 0.0L;
    int vertex = 0;
    int tail_vertex = 0;
    int shifted_tail_vertex = 0;
    long double base = 0.0L;
    long double tail_base = 0.0L;
    long double shifted_tail_base = 0.0L;
    int negative_groups = 0;
    int negative_tails = 0;
    int negative_shifted_tails = 0;
    int groups = 0;
};

CaseResult analyze_case(int level, int label, bool print_detail) {
    const int half_level = level / 2;
    const int mode_count = half_level + 1;
    const long double normalization =
        std::sqrt(2.0L / static_cast<long double>(level + 2));
    std::vector<long double> eigenvalues(
        static_cast<std::size_t>(mode_count)
    );
    std::vector<std::vector<long double>> eigenvectors(
        static_cast<std::size_t>(half_level + 1),
        std::vector<long double>(static_cast<std::size_t>(mode_count))
    );
    for (int mode = 0; mode < mode_count; ++mode) {
        const long double theta =
            static_cast<long double>(mode + 1)
            * std::acos(-1.0L)
            / static_cast<long double>(level + 2);
        eigenvalues[static_cast<std::size_t>(mode)] =
            character_value(label, theta);
        const long double pairing_factor =
            mode == half_level ? 1.0L : std::sqrt(2.0L);
        for (int half_vertex = 0;
             half_vertex <= half_level;
             ++half_vertex) {
            const int vertex = 2 * half_vertex;
            eigenvectors[static_cast<std::size_t>(half_vertex)][
                static_cast<std::size_t>(mode)
            ] = pairing_factor * normalization * std::sin(
                static_cast<long double>(vertex + 1) * theta
            );
        }
    }

    CaseResult result;
    bool initialized = false;
    const int q_index = label / 2;
    for (int half_vertex = 1;
         half_vertex <= half_level;
         ++half_vertex) {
        std::vector<SpectralPair> pairs;
        for (int first = 0; first < mode_count; ++first) {
            for (int second = first + 1;
                 second < mode_count;
                 ++second) {
                const long double lambda_first =
                    eigenvalues[static_cast<std::size_t>(first)];
                const long double lambda_second =
                    eigenvalues[static_cast<std::size_t>(second)];
                const long double row_minor =
                    eigenvectors[0][static_cast<std::size_t>(first)]
                        * eigenvectors[
                            static_cast<std::size_t>(q_index)
                        ][static_cast<std::size_t>(second)]
                    - eigenvectors[0][static_cast<std::size_t>(second)]
                        * eigenvectors[
                            static_cast<std::size_t>(q_index)
                        ][static_cast<std::size_t>(first)];
                const long double column_minor =
                    eigenvectors[0][static_cast<std::size_t>(first)]
                        * eigenvectors[
                            static_cast<std::size_t>(half_vertex)
                        ][static_cast<std::size_t>(second)]
                    - eigenvectors[0][static_cast<std::size_t>(second)]
                        * eigenvectors[
                            static_cast<std::size_t>(half_vertex)
                        ][static_cast<std::size_t>(first)];
                pairs.push_back(SpectralPair{
                    lambda_first * lambda_first
                        * lambda_second * lambda_second,
                    row_minor * column_minor
                });
            }
        }
        std::sort(
            pairs.begin(),
            pairs.end(),
            [](const SpectralPair& left, const SpectralPair& right) {
                return left.base < right.base;
            }
        );
        std::size_t begin = 0U;
        std::vector<SpectralPair> grouped;
        while (begin < pairs.size()) {
            std::size_t end = begin + 1U;
            long double coefficient = pairs[begin].coefficient;
            while (end < pairs.size()) {
                const long double scale = std::max(
                    1.0L,
                    std::max(
                        std::abs(pairs[begin].base),
                        std::abs(pairs[end].base)
                    )
                );
                if (std::abs(pairs[end].base - pairs[begin].base)
                    > 1.0e-14L * scale) {
                    break;
                }
                coefficient += pairs[end].coefficient;
                ++end;
            }
            ++result.groups;
            grouped.push_back(SpectralPair{
                pairs[begin].base,
                coefficient
            });
            if (!initialized
                || coefficient < result.minimum_grouped_residue) {
                initialized = true;
                result.minimum_grouped_residue = coefficient;
                result.vertex = 2 * half_vertex;
                result.base = pairs[begin].base;
            }
            if (coefficient < -1.0e-14L) {
                result.failed = true;
                ++result.negative_groups;
            }
            begin = end;
        }
        long double tail = 0.0L;
        long double shifted_tail = 0.0L;
        for (std::size_t reverse = grouped.size();
             reverse > 0U;
             --reverse) {
            const SpectralPair& group = grouped[reverse - 1U];
            tail += group.coefficient;
            shifted_tail += group.coefficient
                * group.base * group.base;
            if (tail < result.minimum_tail) {
                result.minimum_tail = tail;
                result.tail_vertex = 2 * half_vertex;
                result.tail_base = group.base;
            }
            if (tail < -1.0e-14L) {
                result.tail_failed = true;
                ++result.negative_tails;
            }
            if (shifted_tail < result.minimum_shifted_tail) {
                result.minimum_shifted_tail = shifted_tail;
                result.shifted_tail_vertex = 2 * half_vertex;
                result.shifted_tail_base = group.base;
            }
            if (shifted_tail < -1.0e-12L) {
                result.shifted_tail_failed = true;
                ++result.negative_shifted_tails;
            }
        }
    }
    if (print_detail) {
        std::cout
            << std::setprecision(20)
            << "SU2_KERNEL_SPECTRAL_RESIDUES_CASE"
            << " level=" << level
            << " label=" << label
            << " modes=" << mode_count
            << " groups=" << result.groups
            << " negative_groups=" << result.negative_groups
            << " minimum=" << result.minimum_grouped_residue
            << " witness_vertex=" << result.vertex
            << " witness_base=" << result.base
            << " negative_tails=" << result.negative_tails
            << " minimum_tail=" << result.minimum_tail
            << " tail_vertex=" << result.tail_vertex
            << " tail_base=" << result.tail_base
            << " negative_shifted_tails="
                << result.negative_shifted_tails
            << " minimum_shifted_tail="
                << result.minimum_shifted_tail
            << " shifted_tail_vertex="
                << result.shifted_tail_vertex
            << " shifted_tail_base=" << result.shifted_tail_base
            << " result="
            << (result.shifted_tail_failed
                ? "FAIL_SHIFTED_TAIL_POSITIVITY"
                : "PASS_NUMERICAL_SHIFTED_TAILS")
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
                "usage: analyze_su2_kernel_spectral_residues "
                "MAXIMUM_LEVEL | --case LEVEL LABEL"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        std::uint64_t cases = 0U;
        std::uint64_t failures = 0U;
        std::uint64_t tail_failures = 0U;
        std::uint64_t shifted_tail_failures = 0U;
        bool printed_failure = false;
        bool printed_shifted_failure = false;
        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                ++cases;
                const CaseResult result =
                    analyze_case(level, label, false);
                if (result.failed) {
                    ++failures;
                }
                if (result.tail_failed) {
                    ++tail_failures;
                }
                if (result.shifted_tail_failed) {
                    ++shifted_tail_failures;
                }
                if ((result.failed || result.tail_failed
                        || result.shifted_tail_failed)
                    && !printed_failure) {
                    analyze_case(level, label, true);
                    printed_failure = true;
                }
                if (result.shifted_tail_failed
                    && !printed_shifted_failure) {
                    analyze_case(level, label, true);
                    printed_shifted_failure = true;
                }
            }
        }
        std::cout
            << "SU2_KERNEL_SPECTRAL_RESIDUES"
            << " maximum_level=" << maximum_level
            << " cases=" << cases
            << " failures=" << failures
            << " tail_failures=" << tail_failures
            << " shifted_tail_failures=" << shifted_tail_failures
            << " result="
            << (shifted_tail_failures == 0U
                ? "PASS_NUMERICAL_SHIFTED_TAILS"
                : "FAIL_SHIFTED_TAIL_POSITIVITY")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_KERNEL_SPECTRAL_RESIDUES FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
