#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int positive_argument(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0
        || parsed > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

long long diagonal_count(int left, int right, int total) {
    const int begin = std::max(0, total - right + 1);
    const int end = std::min(left - 1, total);
    return std::max(0, end - begin + 1);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_half_level = argc >= 2
            ? positive_argument(argv[1], "maximum_half_level")
            : 100;
        if (argc > 2 || maximum_half_level < 2) {
            throw std::invalid_argument(
                "usage: analyze_su2_aim_three_box_cyclic_order "
                "[maximum_half_level]");
        }

        std::uint64_t comparisons = 0U;
        std::uint64_t failures = 0U;
        std::uint64_t prefix_comparisons = 0U;
        std::uint64_t prefix_failures = 0U;
        std::uint64_t distance_failures = 0U;
        for (int half_level = 2;
             half_level <= maximum_half_level;
             ++half_level) {
            const int period = 2 * half_level + 2;
            for (int left = 1; left < period; ++left) {
                for (int padding = 0;
                     padding <= half_level;
                     ++padding) {
                    const int right = left + 2 * padding;
                    const int span = left + right;
                    const int maximum_complement
                        = std::min(period - 1, span - period - 2);
                    if (maximum_complement < 1) {
                        continue;
                    }
                    std::vector<long long> cyclic(
                        static_cast<std::size_t>(period),
                        0);
                    for (int residue = 0;
                         residue < period;
                         ++residue) {
                        for (int lift = -1; lift <= 3; ++lift) {
                            cyclic[static_cast<std::size_t>(residue)]
                                += diagonal_count(
                                    left,
                                    right,
                                    residue + lift * period - 1);
                        }
                    }
                    for (int label = 1;
                         label <= 2 * (half_level / 2);
                         ++label) {
                        for (int complement = 1;
                             complement <= maximum_complement;
                             ++complement) {
                            long long prefix = 0;
                            long long boundary = 0;
                            for (int index = 1;
                                 index <= complement;
                                 ++index) {
                                ++comparisons;
                                const int first
                                    = (label + index) % period;
                                const int second
                                    = ((index - label - 1) % period
                                       + period)
                                    % period;
                                if (cyclic[
                                        static_cast<std::size_t>(first)]
                                    < cyclic[
                                        static_cast<std::size_t>(second)]) {
                                    ++failures;
                                }
                                prefix += cyclic[
                                              static_cast<std::size_t>(
                                                  first)]
                                          - cyclic[
                                              static_cast<std::size_t>(
                                                  second)];
                            }
                            const int width = period - complement;
                            for (int layer = 0;
                                 layer < std::min(width, label);
                                 ++layer) {
                                boundary += diagonal_count(
                                    left,
                                    right,
                                    label - layer - 1);
                            }
                            ++prefix_comparisons;
                            if (prefix + boundary < 0) {
                                ++prefix_failures;
                                std::cout
                                    << "PREFIX_FAIL half_level="
                                    << half_level
                                    << " label=" << label
                                    << " left=" << left
                                    << " padding=" << padding
                                    << " complement=" << complement
                                    << " cyclic_value=" << prefix
                                    << " boundary=" << boundary << '\n';
                                return EXIT_FAILURE;
                            }
                            const int twice_center
                                = span + complement - 1;
                            const int first_point
                                = complement + label;
                            const int second_point
                                = complement - label - 1;
                            const auto cyclic_distance =
                                [period](int twice_point,
                                         int center) {
                                    int residue
                                        = (twice_point - center)
                                        % (2 * period);
                                    if (residue < 0) {
                                        residue += 2 * period;
                                    }
                                    return std::min(
                                        residue,
                                        2 * period - residue);
                                };
                            if (cyclic_distance(
                                    2 * first_point,
                                    twice_center)
                                > cyclic_distance(
                                    2 * second_point,
                                    twice_center)) {
                                ++distance_failures;
                            }
                        }
                    }
                }
            }
        }
        std::cout << "SU2_AIM_THREE_BOX_CYCLIC_ORDER"
                  << " maximum_half_level=" << maximum_half_level
                  << " comparisons=" << comparisons
                  << " failures=" << failures
                  << " prefix_comparisons=" << prefix_comparisons
                  << " prefix_failures=" << prefix_failures
                  << " distance_failures=" << distance_failures
                  << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
