#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

using Count = std::uint64_t;
using Multiplicities = std::map<int, Count>;

Multiplicities fuse(const Multiplicities& input, const int label) {
    Multiplicities output;
    for (const auto& [source, multiplicity] : input) {
        for (int target = std::abs(source - label);
             target <= source + label;
             target += 2) {
            output[target] += multiplicity;
        }
    }
    return output;
}

Count invariant_multiplicity(const int number_of_twos) {
    Multiplicities multiplicities{{3, 1}};
    multiplicities = fuse(multiplicities, 1);
    for (int index = 0; index < number_of_twos; ++index) {
        multiplicities = fuse(multiplicities, 2);
    }
    return multiplicities[0];
}

Count two_power_multiplicity(
    const int target,
    const int number_of_twos
) {
    Multiplicities multiplicities{{0, 1}};
    for (int index = 0; index < number_of_twos; ++index) {
        multiplicities = fuse(multiplicities, 2);
    }
    return multiplicities[target];
}

Count binomial(const int n, const int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    Count value = 1;
    for (int index = 1; index <= k; ++index) {
        value = value * static_cast<Count>(n - k + index)
            / static_cast<Count>(index);
    }
    return value;
}

Count bounded_compositions(
    const int number_of_variables,
    const int maximum_part,
    const int total
) {
    std::vector<Count> coefficients(
        static_cast<std::size_t>(total + 1), 0
    );
    coefficients[0] = 1;
    for (int variable = 0; variable < number_of_variables; ++variable) {
        std::vector<Count> next(
            static_cast<std::size_t>(total + 1), 0
        );
        for (int degree = 0; degree <= total; ++degree) {
            for (int part = 0;
                 part <= maximum_part && degree + part <= total;
                 ++part) {
                next[static_cast<std::size_t>(degree + part)]
                    += coefficients[static_cast<std::size_t>(degree)];
            }
        }
        coefficients = std::move(next);
    }
    return coefficients[static_cast<std::size_t>(total)];
}

bool has_evaluation_neighbor(const std::array<int, 5>& grades) {
    for (unsigned int mask = 1U; mask < 31U; ++mask) {
        const int active_size = static_cast<int>(std::popcount(mask));
        if (active_size > 4) {
            continue;
        }
        int inactive_sum = 0;
        for (int index = 0; index < 5; ++index) {
            if ((mask & (1U << index)) == 0U) {
                inactive_sum += grades[static_cast<std::size_t>(index)];
            }
        }
        if (inactive_sum == 4 - active_size) {
            return true;
        }
    }
    return false;
}

struct ResidualLeft {
    std::vector<std::size_t> neighbors;
};

void generate_off_allocations(
    const std::vector<int>& inactive,
    const std::size_t position,
    const int remaining,
    std::array<int, 5>& allocation,
    std::vector<std::array<int, 5>>& output
) {
    if (position == inactive.size()) {
        if (remaining == 0) {
            output.push_back(allocation);
        }
        return;
    }
    const int coordinate = inactive[position];
    for (int value = 0; value <= 2 && value <= remaining; ++value) {
        allocation[static_cast<std::size_t>(coordinate)] = value;
        generate_off_allocations(
            inactive, position + 1, remaining - value, allocation, output
        );
    }
    allocation[static_cast<std::size_t>(coordinate)] = -1;
}

bool augment(
    const std::size_t left,
    const std::vector<ResidualLeft>& graph,
    std::vector<int>& matched_right,
    std::vector<bool>& seen
) {
    for (const std::size_t right : graph[left].neighbors) {
        if (seen[right]) {
            continue;
        }
        seen[right] = true;
        if (matched_right[right] < 0
            || augment(
                static_cast<std::size_t>(matched_right[right]),
                graph,
                matched_right,
                seen
            )) {
            matched_right[right] = static_cast<int>(left);
            return true;
        }
    }
    return false;
}

}  // namespace

int main() {
    try {
        const std::array<Count, 4> expected_multiplicities{1, 2, 5, 12};
        const std::array<Count, 4> expected_off_allocations{16, 6, 2, 1};
        const std::array<Count, 4> expected_left_counts{80, 120, 100, 60};

        std::array<Count, 4> multiplicities{};
        std::array<Count, 4> off_allocations{};
        std::array<Count, 4> left_counts{};
        Count left_total = 0;
        for (int active_twos = 1; active_twos <= 4; ++active_twos) {
            const std::size_t index
                = static_cast<std::size_t>(active_twos - 1);
            multiplicities[index] = invariant_multiplicity(active_twos);
            off_allocations[index] = bounded_compositions(
                5 - active_twos, 2, 4 - active_twos
            );
            left_counts[index] = binomial(5, active_twos)
                * multiplicities[index] * off_allocations[index];
            left_total += left_counts[index];
        }

        Count ambient_scalar = 0;
        Count neighborhood = 0;
        Count zero_characterization = 0;
        std::vector<std::array<int, 7>> scalar_allocations;
        for (int minus_grade = 0; minus_grade <= 3; ++minus_grade) {
            for (int one_grade = 0; one_grade <= 1; ++one_grade) {
                for (int first = 0; first <= 2; ++first) {
                    for (int second = 0; second <= 2; ++second) {
                        for (int third = 0; third <= 2; ++third) {
                            for (int fourth = 0; fourth <= 2; ++fourth) {
                                for (int fifth = 0; fifth <= 2; ++fifth) {
                                    const std::array<int, 5> grades{
                                        first, second, third, fourth, fifth
                                    };
                                    int total = minus_grade + one_grade;
                                    for (const int grade : grades) {
                                        total += grade;
                                    }
                                    if (total != 6) {
                                        continue;
                                    }
                                    ++ambient_scalar;
                                    scalar_allocations.push_back({
                                        minus_grade,
                                        one_grade,
                                        first,
                                        second,
                                        third,
                                        fourth,
                                        fifth
                                    });
                                    if (has_evaluation_neighbor(grades)) {
                                        ++neighborhood;
                                    }
                                    bool has_zero = false;
                                    for (const int grade : grades) {
                                        has_zero = has_zero || grade == 0;
                                    }
                                    if (has_zero) {
                                        ++zero_characterization;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        const std::array<Count, 4> expected_carrier_multiplicities{
            1, 1, 3, 6
        };
        const std::array<Count, 4> expected_residual_multiplicities{
            0, 1, 2, 6
        };
        const std::array<Count, 4> expected_carrier_counts{80, 60, 60, 30};
        std::array<Count, 4> carrier_multiplicities{};
        std::array<Count, 4> residual_multiplicities{};
        std::array<Count, 4> carrier_counts{};
        Count carrier_capacity = 0;
        Count residual_total = 0;
        for (int active_twos = 1; active_twos <= 4; ++active_twos) {
            const std::size_t index
                = static_cast<std::size_t>(active_twos - 1);
            carrier_multiplicities[index]
                = two_power_multiplicity(2, active_twos);
            residual_multiplicities[index]
                = two_power_multiplicity(4, active_twos);
            carrier_counts[index] = binomial(5, active_twos)
                * off_allocations[index] * carrier_multiplicities[index];
            carrier_capacity += carrier_counts[index];
            residual_total += binomial(5, active_twos)
                * off_allocations[index] * residual_multiplicities[index];
        }

        std::vector<ResidualLeft> residual_graph;
        for (unsigned int mask = 1U; mask < 31U; ++mask) {
            const int active_twos = static_cast<int>(std::popcount(mask));
            if (active_twos > 4) {
                continue;
            }
            const Count copies = residual_multiplicities[
                static_cast<std::size_t>(active_twos - 1)
            ];
            if (copies == 0) {
                continue;
            }
            std::vector<int> inactive;
            for (int coordinate = 0; coordinate < 5; ++coordinate) {
                if ((mask & (1U << coordinate)) == 0U) {
                    inactive.push_back(coordinate);
                }
            }
            std::array<int, 5> allocation{-1, -1, -1, -1, -1};
            std::vector<std::array<int, 5>> off_allocations_for_mask;
            generate_off_allocations(
                inactive,
                0,
                4 - active_twos,
                allocation,
                off_allocations_for_mask
            );
            for (const auto& fixed : off_allocations_for_mask) {
                ResidualLeft left;
                for (std::size_t right = 0;
                     right < scalar_allocations.size();
                     ++right) {
                    bool compatible = true;
                    for (int coordinate = 0;
                         coordinate < 5;
                         ++coordinate) {
                        if ((mask & (1U << coordinate)) == 0U
                            && scalar_allocations[right][
                                static_cast<std::size_t>(coordinate + 2)
                            ] != fixed[static_cast<std::size_t>(coordinate)]) {
                            compatible = false;
                        }
                    }
                    if (compatible) {
                        left.neighbors.push_back(right);
                    }
                }
                for (Count copy = 0; copy < copies; ++copy) {
                    residual_graph.push_back(left);
                }
            }
        }

        std::vector<int> matched_right(scalar_allocations.size(), -1);
        Count residual_matching_rank = 0;
        for (std::size_t left = 0; left < residual_graph.size(); ++left) {
            std::vector<bool> seen(scalar_allocations.size(), false);
            if (augment(left, residual_graph, matched_right, seen)) {
                ++residual_matching_rank;
            }
        }
        const Count reduced_essential_rank
            = carrier_capacity + residual_matching_rank;

        if (multiplicities != expected_multiplicities
            || off_allocations != expected_off_allocations
            || left_counts != expected_left_counts
            || left_total != 360
            || ambient_scalar != 312
            || neighborhood != 305
            || zero_characterization != neighborhood
            || carrier_multiplicities != expected_carrier_multiplicities
            || residual_multiplicities != expected_residual_multiplicities
            || carrier_counts != expected_carrier_counts
            || carrier_capacity != 230
            || residual_total != 130
            || residual_graph.size() != residual_total
            || residual_matching_rank != residual_total
            || reduced_essential_rank != left_total) {
            throw std::runtime_error("False Target 23A6A witness mismatch");
        }

        std::cout << "FALSE_TARGET_23A6A PASS\n"
                  << "fusion_multiplicities=1,2,5,12\n"
                  << "off_support_allocations=16,6,2,1\n"
                  << "left_by_support_size=80,120,100,60\n"
                  << "left=360\n"
                  << "ambient_scalar=312\n"
                  << "excluded_scalar=7\n"
                  << "neighborhood=305\n"
                  << "hall_deficit=55\n"
                  << "carrier_multiplicities=1,1,3,6\n"
                  << "residual_multiplicities=0,1,2,6\n"
                  << "carrier_by_support_size=80,60,60,30\n"
                  << "carrier_capacity=230\n"
                  << "residual_left=130\n"
                  << "residual_matching_rank=130\n"
                  << "reduced_essential_rank=360\n"
                  << "paired_h1_columns_not_included=230\n"
                  << "full_schur_complement_checked=NO\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
