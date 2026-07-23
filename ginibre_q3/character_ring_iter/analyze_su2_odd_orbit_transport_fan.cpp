#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Term {
    int sign = 1;
    long double coefficient = 0.0L;
    std::vector<long double> lambda;
};

struct Chamber {
    int rank = 0;
    int support_code = 0;
    int parity_code = 0;
    std::vector<int> signs;
    std::vector<int> minimum_power;
    std::vector<Term> terms;
    std::vector<int> negative;
    std::vector<int> positive;
};

struct Statistics {
    std::uint64_t chambers = 0U;
    std::uint64_t pointwise_chambers = 0U;
    std::uint64_t regimes = 0U;
    std::uint64_t direct_nodes = 0U;
    std::uint64_t braid_nodes = 0U;
    std::uint64_t rational_fan_nodes = 0U;
    std::uint64_t shifted_nodes = 0U;
    std::uint64_t exact_leaves = 0U;
    std::uint64_t separated_pole_regimes = 0U;
    int maximum_shift = 0;
    int maximum_depth = 0;
    int maximum_fan_bound = 0;
};

bool is_separated_pole_regime(const Chamber& chamber) {
    const int last_label = chamber.rank - 2;
    for (int label = 0; label < chamber.rank - 1; ++label) {
        const bool active = label == 0 || label == last_label;
        const std::size_t index = static_cast<std::size_t>(label);
        if (active) {
            if (chamber.signs[index] != -1 ||
                chamber.minimum_power[index] != 1) {
                return false;
            }
        } else if (chamber.signs[index] != 0) {
            return false;
        }
    }
    return true;
}

long double magnitude_at(
    const Term& term,
    const std::vector<int>& floors
) {
    long double result = term.coefficient;
    for (std::size_t label = 0; label < floors.size(); ++label) {
        result *= std::pow(
            term.lambda[label],
            static_cast<long double>(floors[label])
        );
    }
    return result;
}

long double transformed_base(
    const Term& term,
    const std::vector<int>& coordinate
) {
    long double result = 1.0L;
    for (int label : coordinate) {
        result *= term.lambda[static_cast<std::size_t>(label)];
    }
    return result;
}

bool approximately_at_least(long double left, long double right) {
    const long double tolerance
        = 1.0e-12L * (1.0L + std::abs(left) + std::abs(right));
    return left + tolerance >= right;
}

bool hall_transport(
    const Chamber& chamber,
    const std::vector<int>& floors,
    const std::vector<std::vector<int>>& coordinates
) {
    std::vector<long double> magnitudes;
    magnitudes.reserve(chamber.terms.size());
    for (const Term& term : chamber.terms) {
        magnitudes.push_back(magnitude_at(term, floors));
    }

    std::vector<std::vector<int>> edges(chamber.negative.size());
    for (std::size_t ni = 0; ni < chamber.negative.size(); ++ni) {
        const Term& negative = chamber.terms[static_cast<std::size_t>(
            chamber.negative[ni]
        )];
        for (std::size_t pi_index = 0;
             pi_index < chamber.positive.size();
             ++pi_index) {
            const Term& positive = chamber.terms[static_cast<std::size_t>(
                chamber.positive[pi_index]
            )];
            bool dominates = true;
            for (const std::vector<int>& coordinate : coordinates) {
                if (!approximately_at_least(
                        transformed_base(positive, coordinate),
                        transformed_base(negative, coordinate)
                    )) {
                    dominates = false;
                    break;
                }
            }
            if (dominates) {
                edges[ni].push_back(static_cast<int>(pi_index));
            }
        }
    }

    if (chamber.negative.size() >= 63U) {
        throw std::runtime_error("too many negative spectral pairs");
    }
    const std::uint64_t subsets
        = std::uint64_t{1} << chamber.negative.size();
    for (std::uint64_t subset = 1U; subset < subsets; ++subset) {
        long double demand = 0.0L;
        std::vector<bool> neighbors(chamber.positive.size(), false);
        for (std::size_t ni = 0; ni < chamber.negative.size(); ++ni) {
            if (((subset >> ni) & 1U) == 0U) {
                continue;
            }
            demand += magnitudes[static_cast<std::size_t>(
                chamber.negative[ni]
            )];
            for (int pi_index : edges[ni]) {
                neighbors[static_cast<std::size_t>(pi_index)] = true;
            }
        }
        long double capacity = 0.0L;
        for (std::size_t pi_index = 0;
             pi_index < chamber.positive.size();
             ++pi_index) {
            if (neighbors[pi_index]) {
                capacity += magnitudes[static_cast<std::size_t>(
                    chamber.positive[pi_index]
                )];
            }
        }
        if (!approximately_at_least(capacity, demand)) {
            return false;
        }
    }
    return true;
}

std::vector<std::vector<int>> cartesian_coordinates(
    const std::vector<int>& variables
) {
    std::vector<std::vector<int>> result;
    result.reserve(variables.size());
    for (int variable : variables) {
        result.push_back(std::vector<int>{variable});
    }
    return result;
}

bool braid_fan_transport(
    const Chamber& chamber,
    const std::vector<int>& floors,
    std::vector<int> variables
) {
    std::sort(variables.begin(), variables.end());
    do {
        std::vector<std::vector<int>> coordinates;
        std::vector<int> cumulative;
        for (int variable : variables) {
            cumulative.push_back(variable);
            coordinates.push_back(cumulative);
        }
        if (!hall_transport(chamber, floors, coordinates)) {
            return false;
        }
    } while (std::next_permutation(variables.begin(), variables.end()));
    return true;
}

std::vector<int> ray_coordinate(
    const std::vector<int>& variables,
    int first_multiplicity,
    int second_multiplicity
) {
    std::vector<int> coordinate;
    coordinate.insert(
        coordinate.end(),
        static_cast<std::size_t>(first_multiplicity),
        variables[0]
    );
    coordinate.insert(
        coordinate.end(),
        static_cast<std::size_t>(second_multiplicity),
        variables[1]
    );
    return coordinate;
}

bool rational_fan_transport(
    const Chamber& chamber,
    const std::vector<int>& floors,
    const std::vector<int>& variables,
    int bound
) {
    if (variables.size() != 2U || bound < 1) {
        return false;
    }
    std::vector<std::pair<int, int>> rays;
    rays.emplace_back(1, 0);
    rays.emplace_back(0, 1);
    for (int first = 1; first <= bound; ++first) {
        for (int second = 1; second <= bound; ++second) {
            if (std::gcd(first, second) == 1) {
                rays.emplace_back(first, second);
            }
        }
    }
    std::sort(
        rays.begin(),
        rays.end(),
        [](const auto& left, const auto& right) {
            return static_cast<long long>(left.second) * right.first
                < static_cast<long long>(right.second) * left.first;
        }
    );
    rays.erase(std::unique(rays.begin(), rays.end()), rays.end());
    for (std::size_t ray = 0; ray + 1U < rays.size(); ++ray) {
        const auto [first_left, first_right] = rays[ray];
        const auto [second_left, second_right] = rays[ray + 1U];
        const std::vector<std::vector<int>> coordinates{
            ray_coordinate(
                variables, first_left, first_right
            ),
            ray_coordinate(
                variables, second_left, second_right
            )
        };
        if (!hall_transport(chamber, floors, coordinates)) {
            return false;
        }
    }
    return true;
}

bool exact_leaf(const Chamber& chamber, const std::vector<int>& floors) {
    long double value = 0.0L;
    long double absolute_sum = 0.0L;
    for (const Term& term : chamber.terms) {
        const long double magnitude = magnitude_at(term, floors);
        value += static_cast<long double>(term.sign) * magnitude;
        absolute_sum += magnitude;
    }
    return value >= -1.0e-11L * (1.0L + absolute_sum);
}

bool certify_region(
    const Chamber& chamber,
    const std::vector<int>& floors,
    const std::vector<int>& variables,
    int maximum_shift,
    int depth,
    Statistics& statistics
) {
    statistics.maximum_depth = std::max(statistics.maximum_depth, depth);
    if (variables.empty()) {
        ++statistics.exact_leaves;
        return exact_leaf(chamber, floors);
    }
    if (hall_transport(
            chamber, floors, cartesian_coordinates(variables)
        )) {
        ++statistics.direct_nodes;
        return true;
    }
    if (braid_fan_transport(chamber, floors, variables)) {
        ++statistics.braid_nodes;
        return true;
    }
    if (variables.size() == 2U) {
        constexpr int fan_bound = 100;
        if (rational_fan_transport(
                chamber, floors, variables, fan_bound
            )) {
            ++statistics.rational_fan_nodes;
            statistics.maximum_fan_bound = std::max(
                statistics.maximum_fan_bound, fan_bound
            );
            return true;
        }
    }

    int successful_shift = -1;
    int successful_fan_bound = 0;
    for (int shift = 1; shift <= maximum_shift; ++shift) {
        std::vector<int> shifted = floors;
        for (int variable : variables) {
            shifted[static_cast<std::size_t>(variable)] += shift;
        }
        if (braid_fan_transport(chamber, shifted, variables)) {
            successful_shift = shift;
            break;
        }
        if (variables.size() == 2U) {
            constexpr int fan_bound = 100;
            if (rational_fan_transport(
                    chamber, shifted, variables, fan_bound
                )) {
                successful_shift = shift;
                successful_fan_bound = fan_bound;
                break;
            }
        }
    }
    if (successful_shift < 0) {
        return false;
    }
    ++statistics.shifted_nodes;
    statistics.maximum_shift = std::max(
        statistics.maximum_shift, successful_shift
    );
    statistics.maximum_fan_bound = std::max(
        statistics.maximum_fan_bound, successful_fan_bound
    );

    for (std::size_t variable_index = 0;
         variable_index < variables.size();
         ++variable_index) {
        const int variable = variables[variable_index];
        std::vector<int> face_variables = variables;
        face_variables.erase(
            face_variables.begin()
                + static_cast<std::ptrdiff_t>(variable_index)
        );
        for (int offset = 0; offset < successful_shift; ++offset) {
            std::vector<int> face_floors = floors;
            face_floors[static_cast<std::size_t>(variable)] += offset;
            if (!certify_region(
                    chamber,
                    face_floors,
                    face_variables,
                    maximum_shift,
                    depth + 1,
                    statistics
                )) {
                return false;
            }
        }
    }
    return true;
}

Chamber make_chamber(
    int rank,
    int support_code,
    int parity_code,
    const std::vector<int>& signs,
    const std::vector<int>& minimum_power,
    const std::vector<std::vector<long double>>& values,
    const std::vector<long double>& weights
) {
    Chamber chamber;
    chamber.rank = rank;
    chamber.support_code = support_code;
    chamber.parity_code = parity_code;
    chamber.signs = signs;
    chamber.minimum_power = minimum_power;
    for (int left = 0; left < rank; ++left) {
        for (int right = left + 1; right < rank; ++right) {
            Term term;
            term.coefficient = 2.0L
                * weights[static_cast<std::size_t>(left)]
                * weights[static_cast<std::size_t>(right)];
            term.lambda.resize(static_cast<std::size_t>(rank - 1), 1.0L);
            bool zero = false;
            for (int label = 0; label < rank - 1; ++label) {
                if (signs[static_cast<std::size_t>(label)] == 0) {
                    continue;
                }
                const long double base
                    = values[static_cast<std::size_t>(left)]
                            [static_cast<std::size_t>(label)]
                    + static_cast<long double>(
                          signs[static_cast<std::size_t>(label)]
                      )
                        * values[static_cast<std::size_t>(right)]
                              [static_cast<std::size_t>(label)];
                if (std::abs(base) < 1.0e-16L) {
                    zero = true;
                    break;
                }
                const int power
                    = minimum_power[static_cast<std::size_t>(label)];
                if (base < 0.0L && (power & 1) != 0) {
                    term.sign = -term.sign;
                }
                term.coefficient *= std::pow(
                    std::abs(base), static_cast<long double>(power)
                );
                term.lambda[static_cast<std::size_t>(label)] = base * base;
            }
            if (!zero) {
                chamber.terms.push_back(std::move(term));
            }
        }
    }
    for (std::size_t term = 0; term < chamber.terms.size(); ++term) {
        if (chamber.terms[term].sign < 0) {
            chamber.negative.push_back(static_cast<int>(term));
        } else {
            chamber.positive.push_back(static_cast<int>(term));
        }
    }
    return chamber;
}

bool analyze_rank(int rank, int maximum_shift, Statistics& statistics) {
    const int order = 2 * rank + 1;
    const int label_count = rank - 1;
    const long double pi = std::acos(-1.0L);
    std::vector<long double> weights(static_cast<std::size_t>(rank), 0.0L);
    std::vector<std::vector<long double>> values(
        static_cast<std::size_t>(rank),
        std::vector<long double>(static_cast<std::size_t>(label_count), 0.0L)
    );
    for (int node = 0; node < rank; ++node) {
        const long double angle
            = pi * static_cast<long double>(node + 1)
            / static_cast<long double>(order);
        weights[static_cast<std::size_t>(node)]
            = 4.0L * std::sin(angle) * std::sin(angle)
            / static_cast<long double>(order);
        for (int label = 0; label < label_count; ++label) {
            values[static_cast<std::size_t>(node)]
                  [static_cast<std::size_t>(label)]
                = std::sin(
                      static_cast<long double>(2 * label + 3) * angle
                  ) / std::sin(angle);
        }
    }

    int support_patterns = 1;
    for (int label = 0; label < label_count; ++label) {
        support_patterns *= 3;
    }
    for (int support_code = 0;
         support_code < support_patterns;
         ++support_code) {
        int code = support_code;
        std::vector<int> signs(static_cast<std::size_t>(label_count), 0);
        int active_count = 0;
        int minus_count = 0;
        for (int label = 0; label < label_count; ++label) {
            const int digit = code % 3;
            code /= 3;
            if (digit != 0) {
                signs[static_cast<std::size_t>(label)]
                    = digit == 1 ? 1 : -1;
                ++active_count;
                if (digit == 2) {
                    ++minus_count;
                }
            }
        }
        if (active_count == 0 || minus_count == 0) {
            continue;
        }
        const int parity_patterns = 1 << active_count;
        for (int parity_code = 0;
             parity_code < parity_patterns;
             ++parity_code) {
            std::vector<int> minimum_power(
                static_cast<std::size_t>(label_count), 0
            );
            int active_index = 0;
            int minus_parity = 0;
            for (int label = 0; label < label_count; ++label) {
                if (signs[static_cast<std::size_t>(label)] == 0) {
                    continue;
                }
                const bool odd
                    = ((parity_code >> active_index) & 1) != 0;
                minimum_power[static_cast<std::size_t>(label)]
                    = odd ? 1 : 2;
                if (signs[static_cast<std::size_t>(label)] == -1 && odd) {
                    minus_parity ^= 1;
                }
                ++active_index;
            }
            if (minus_parity != 0) {
                continue;
            }
            const Chamber chamber = make_chamber(
                rank,
                support_code,
                parity_code,
                signs,
                minimum_power,
                values,
                weights
            );
            if (chamber.negative.empty()) {
                ++statistics.pointwise_chambers;
                continue;
            }
            ++statistics.chambers;
            int residual_patterns = 1;
            for (int label = 0; label < label_count; ++label) {
                if (signs[static_cast<std::size_t>(label)] != 0) {
                    residual_patterns *= 2;
                }
            }
            for (int residual_code = 0;
                 residual_code < residual_patterns;
                 ++residual_code) {
                std::vector<int> floors(
                    static_cast<std::size_t>(label_count), 0
                );
                std::vector<int> variables;
                int residual_index = 0;
                for (int label = 0; label < label_count; ++label) {
                    if (signs[static_cast<std::size_t>(label)] == 0) {
                        continue;
                    }
                    if (((residual_code >> residual_index) & 1) != 0) {
                        floors[static_cast<std::size_t>(label)] = 1;
                        variables.push_back(label);
                    }
                    ++residual_index;
                }
                ++statistics.regimes;
                if (is_separated_pole_regime(chamber)) {
                    ++statistics.separated_pole_regimes;
                    continue;
                }
                if (!certify_region(
                        chamber,
                        floors,
                        variables,
                        maximum_shift,
                        0,
                        statistics
                    )) {
                    std::cout << "SU2_ODD_ORBIT_TRANSPORT_FAN result=FAIL"
                              << " rank=" << rank
                              << " support_code=" << support_code
                              << " parity_code=" << parity_code
                              << " residual_code=" << residual_code
                              << " variables=";
                    for (int variable : variables) {
                        std::cout << variable + 1 << ',';
                    }
                    std::cout << '\n';
                    for (std::size_t term_index = 0;
                         term_index < chamber.terms.size();
                         ++term_index) {
                        const Term& term = chamber.terms[term_index];
                        std::cout << "  term=" << term_index
                                  << " sign=" << term.sign
                                  << " coefficient=" << term.coefficient
                                  << " lambda=";
                        for (int variable : variables) {
                            std::cout << term.lambda[
                                static_cast<std::size_t>(variable)
                            ] << ',';
                        }
                        std::cout << '\n';
                    }
                    return false;
                }
            }
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout << std::setprecision(15);
        if (argc != 3) {
            std::cerr << "usage: " << argv[0]
                      << " RANK MAXIMUM_SHIFT\n";
            return EXIT_FAILURE;
        }
        const int rank = std::stoi(argv[1]);
        const int maximum_shift = std::stoi(argv[2]);
        if (rank < 3 || rank > 8 || maximum_shift < 1) {
            throw std::runtime_error("invalid transport-fan bounds");
        }
        Statistics statistics;
        if (!analyze_rank(rank, maximum_shift, statistics)) {
            return EXIT_FAILURE;
        }
        std::cout << "SU2_ODD_ORBIT_TRANSPORT_FAN"
                  << " rank=" << rank
                  << " level=" << 2 * rank - 1
                  << " chambers=" << statistics.chambers
                  << " pointwise_chambers="
                  << statistics.pointwise_chambers
                  << " regimes=" << statistics.regimes
                  << " direct_nodes=" << statistics.direct_nodes
                  << " braid_nodes=" << statistics.braid_nodes
                  << " rational_fan_nodes="
                  << statistics.rational_fan_nodes
                  << " shifted_nodes=" << statistics.shifted_nodes
                  << " exact_leaves=" << statistics.exact_leaves
                  << " separated_pole_regimes="
                  << statistics.separated_pole_regimes
                  << " maximum_shift=" << statistics.maximum_shift
                  << " maximum_depth=" << statistics.maximum_depth
                  << " maximum_fan_bound="
                  << statistics.maximum_fan_bound
                  << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
