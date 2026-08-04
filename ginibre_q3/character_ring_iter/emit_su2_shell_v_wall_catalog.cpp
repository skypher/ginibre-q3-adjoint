#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <z3++.h>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Matrix = std::vector<std::vector<int>>;
using IntegerMatrix = std::vector<std::vector<Integer>>;

enum Coordinate : std::size_t {
    K = 0U,
    Q = 1U,
    V = 2U,
    X = 3U,
    Constant = 4U
};

struct Affine {
    std::array<int, 5U> coefficient{};
};

int parse_positive(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value <= 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error("bound must be a positive integer");
    }
    return static_cast<int>(value);
}

int parse_nonnegative(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value < 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error("bound must be a nonnegative integer");
    }
    return static_cast<int>(value);
}

bool fuses_half(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= source + label
        && source + target + label <= 2 * level;
}

std::vector<Integer> multiply_row(
    const std::vector<Integer>& row,
    const Matrix& matrix
) {
    std::vector<Integer> result(matrix.size());
    for (std::size_t source = 0U; source < matrix.size(); ++source) {
        if (row[source] == 0) {
            continue;
        }
        for (std::size_t target = 0U; target < matrix.size(); ++target) {
            const int entry = matrix[source][target];
            if (entry != 0) {
                result[target] += row[source] * entry;
            }
        }
    }
    return result;
}

IntegerMatrix multiply_right(const IntegerMatrix& left, const Matrix& right) {
    IntegerMatrix result(
        left.size(),
        std::vector<Integer>(right.size())
    );
    for (std::size_t source = 0U; source < left.size(); ++source) {
        for (std::size_t middle = 0U; middle < right.size(); ++middle) {
            if (left[source][middle] == 0) {
                continue;
            }
            for (std::size_t target = 0U; target < right.size(); ++target) {
                const int entry = right[middle][target];
                if (entry != 0) {
                    result[source][target] += left[source][middle] * entry;
                }
            }
        }
    }
    return result;
}

IntegerMatrix identity_matrix(int size) {
    IntegerMatrix result(
        static_cast<std::size_t>(size),
        std::vector<Integer>(static_cast<std::size_t>(size))
    );
    for (int index = 0; index < size; ++index) {
        result[static_cast<std::size_t>(index)]
              [static_cast<std::size_t>(index)] = 1;
    }
    return result;
}

Affine coordinate(Coordinate index, int coefficient = 1) {
    Affine result;
    result.coefficient[index] = coefficient;
    return result;
}

Affine constant(int value) {
    return coordinate(Constant, value);
}

Affine add(const Affine& left, const Affine& right) {
    Affine result;
    for (std::size_t index = 0U; index < result.coefficient.size(); ++index) {
        result.coefficient[index] =
            left.coefficient[index] + right.coefficient[index];
    }
    return result;
}

Affine subtract(const Affine& left, const Affine& right) {
    Affine result;
    for (std::size_t index = 0U; index < result.coefficient.size(); ++index) {
        result.coefficient[index] =
            left.coefficient[index] - right.coefficient[index];
    }
    return result;
}

Affine scale(int scalar, const Affine& value) {
    Affine result;
    for (std::size_t index = 0U; index < result.coefficient.size(); ++index) {
        result.coefficient[index] = scalar * value.coefficient[index];
    }
    return result;
}

bool normalize(Affine& wall) {
    int divisor = 0;
    for (const int coefficient : wall.coefficient) {
        divisor = std::gcd(divisor, std::abs(coefficient));
    }
    if (divisor == 0) {
        return false;
    }
    for (int& coefficient : wall.coefficient) {
        coefficient /= divisor;
    }
    for (const int coefficient : wall.coefficient) {
        if (coefficient == 0) {
            continue;
        }
        if (coefficient < 0) {
            for (int& entry : wall.coefficient) {
                entry = -entry;
            }
        }
        break;
    }
    return true;
}

struct AffineLess {
    bool operator()(const Affine& left, const Affine& right) const {
        return left.coefficient < right.coefficient;
    }
};

void insert_wall(std::set<Affine, AffineLess>& walls, Affine wall) {
    if (normalize(wall)) {
        walls.insert(wall);
    }
}

void insert_comparisons(
    std::set<Affine, AffineLess>& walls,
    const std::vector<Affine>& first,
    const std::vector<Affine>& second
) {
    for (const Affine& left : first) {
        for (const Affine& right : second) {
            insert_wall(walls, subtract(left, right));
        }
    }
}

void insert_internal_comparisons(
    std::set<Affine, AffineLess>& walls,
    const std::vector<Affine>& candidates
) {
    for (std::size_t left = 0U; left < candidates.size(); ++left) {
        for (std::size_t right = left + 1U;
             right < candidates.size();
             ++right) {
            insert_wall(walls, subtract(candidates[left], candidates[right]));
        }
    }
}

void insert_rounded_comparisons(
    std::set<Affine, AffineLess>& walls,
    const std::vector<Affine>& first,
    const std::vector<Affine>& second
) {
    for (const Affine& left : first) {
        for (const Affine& right : second) {
            for (int shift = -1; shift <= 1; ++shift) {
                insert_wall(
                    walls,
                    add(subtract(left, right), constant(shift))
                );
            }
        }
    }
}

void insert_internal_rounded_comparisons(
    std::set<Affine, AffineLess>& walls,
    const std::vector<Affine>& candidates
) {
    for (std::size_t left = 0U; left < candidates.size(); ++left) {
        for (std::size_t right = left + 1U;
             right < candidates.size();
             ++right) {
            insert_rounded_comparisons(
                walls,
                {candidates[left]},
                {candidates[right]}
            );
        }
    }
}

void append_crossing_walls(std::set<Affine, AffineLess>& walls) {
    const Affine k = coordinate(K);
    const Affine q = coordinate(Q);
    const Affine v = coordinate(V);
    for (int power = 1; power <= 5; ++power) {
        const Affine product = scale(power, q);
        for (int branch = 0; branch < 3; ++branch) {
            Affine raw_lower;
            Affine raw_upper;
            if (branch == 0) {
                raw_lower = subtract(
                    scale(power - 1, q),
                    v
                );
                raw_upper = add(
                    add(scale(power + 1, q), scale(-1, k)),
                    v
                );
            } else if (branch == 1) {
                raw_lower = add(
                    add(scale(power - 1, q), scale(-1, k)),
                    add(scale(-1, v), constant(-1))
                );
                raw_upper = add(
                    add(scale(power + 1, q), scale(-2, k)),
                    add(v, constant(-1))
                );
            } else {
                raw_lower = add(
                    add(scale(power - 1, q), scale(-2, k)),
                    add(scale(-1, v), constant(-2))
                );
                raw_upper = add(
                    add(scale(power + 1, q), scale(-3, k)),
                    add(v, constant(-2))
                );
            }
            const std::vector<Affine> lower{constant(0), raw_lower};
            const std::vector<Affine> upper{product, raw_upper};
            insert_internal_comparisons(walls, lower);
            insert_internal_comparisons(walls, upper);
            insert_comparisons(walls, lower, upper);
            for (int image = 1; image <= power; ++image) {
                const Affine activation = add(
                    scale(2 * image, q),
                    constant(image)
                );
                insert_comparisons(walls, lower, {activation});
                insert_comparisons(walls, upper, {activation});
            }
        }
    }
}

void append_terminal_walls(std::set<Affine, AffineLess>& walls) {
    const Affine k = coordinate(K);
    const Affine q = coordinate(Q);
    const Affine v = coordinate(V);
    const Affine x = coordinate(X);
    const Affine lower_level = add(k, constant(-1));
    const Affine width = add(lower_level, scale(-2, q));
    const Affine half_period = add(k, constant(1));
    for (int power = 1; power <= 4; ++power) {
        const Affine target = power % 2 == 0
            ? scale(2, x)
            : subtract(lower_level, scale(2, x));
        const Affine source = scale(2, v);
        const std::vector<Affine> lower_labels{
            subtract(source, target),
            subtract(target, source)
        };
        const std::vector<Affine> upper_labels{
            add(source, target),
            subtract(scale(2, lower_level), add(source, target))
        };
        insert_internal_rounded_comparisons(walls, lower_labels);
        insert_internal_rounded_comparisons(walls, upper_labels);
        const Affine total = scale(power, width);
        for (int branch = 0; branch < 4; ++branch) {
            const Affine branch_lower = add(
                subtract(total, scale(branch + 1, half_period)),
                constant(2)
            );
            const Affine branch_upper = subtract(
                total,
                scale(branch, half_period)
            );
            for (const Affine& lower_label : lower_labels) {
                for (const Affine& upper_label : upper_labels) {
                    Affine fusion_lower;
                    Affine fusion_upper;
                    if (branch % 2 == 0) {
                        fusion_lower = subtract(
                            subtract(total, scale(branch, half_period)),
                            upper_label
                        );
                        fusion_upper = subtract(
                            subtract(total, scale(branch, half_period)),
                            lower_label
                        );
                    } else {
                        fusion_lower = add(
                            add(lower_label, total),
                            add(scale(-(branch + 1), half_period), constant(2))
                        );
                        fusion_upper = add(
                            add(upper_label, total),
                            add(scale(-(branch + 1), half_period), constant(2))
                        );
                    }
                    const std::vector<Affine> lower{
                        constant(0), branch_lower, fusion_lower
                    };
                    const std::vector<Affine> upper{
                        total, branch_upper, fusion_upper
                    };
                    insert_internal_rounded_comparisons(walls, lower);
                    insert_internal_rounded_comparisons(walls, upper);
                    insert_rounded_comparisons(walls, lower, upper);
                    if (power >= 3) {
                        const Affine upper_activation = scale(
                            2,
                            add(width, constant(1))
                        );
                        const Affine lower_activation = add(
                            scale(2, width),
                            constant(3)
                        );
                        insert_comparisons(
                            walls,
                            lower,
                            {lower_activation}
                        );
                        insert_comparisons(
                            walls,
                            upper,
                            {upper_activation}
                        );
                    }
                }
            }
        }
    }
}

void append_return_walls(std::set<Affine, AffineLess>& walls) {
    const Affine k = coordinate(K);
    const Affine q = coordinate(Q);
    insert_wall(
        walls,
        add(add(scale(5, q), scale(-2, k)), constant(-1))
    );
}

void print_wall(const Affine& wall) {
    static constexpr std::array<const char*, 5U> names{
        "K", "Q", "V", "x", "1"
    };
    bool first = true;
    for (std::size_t index = 0U; index < wall.coefficient.size(); ++index) {
        const int coefficient = wall.coefficient[index];
        if (coefficient == 0) {
            continue;
        }
        if (!first && coefficient > 0) {
            std::cout << '+';
        }
        if (coefficient == -1) {
            std::cout << '-';
        } else if (coefficient != 1) {
            std::cout << coefficient << '*';
        }
        std::cout << names[index];
        first = false;
    }
    std::cout << "=0\n";
}

std::size_t count_v_walls(const std::set<Affine, AffineLess>& walls) {
    return static_cast<std::size_t>(std::count_if(
        walls.begin(),
        walls.end(),
        [](const Affine& wall) {
            return wall.coefficient[V] != 0;
        }
    ));
}

std::string v_coefficient_set(const std::set<Affine, AffineLess>& walls) {
    std::set<int> coefficients;
    for (const Affine& wall : walls) {
        const int coefficient = std::abs(wall.coefficient[V]);
        if (coefficient != 0) {
            coefficients.insert(coefficient);
        }
    }
    std::string result;
    for (const int coefficient : coefficients) {
        if (!result.empty()) {
            result += ',';
        }
        result += std::to_string(coefficient);
    }
    return result;
}

z3::expr z3_affine(
    z3::context& context,
    const Affine& wall,
    const z3::expr& level,
    const z3::expr& label,
    const z3::expr& crossing,
    const z3::expr& target
) {
    z3::expr result = context.int_val(wall.coefficient[Constant]);
    const std::array<z3::expr, 4U> variables{
        level, label, crossing, target
    };
    for (std::size_t index = 0U; index < variables.size(); ++index) {
        const int coefficient = wall.coefficient[index];
        if (coefficient != 0) {
            result = result + coefficient * variables[index];
        }
    }
    return result;
}

struct FeasibilitySummary {
    std::size_t moving = 0U;
    std::size_t feasible_moving = 0U;
    std::size_t parameter = 0U;
    std::size_t feasible_parameter = 0U;
    std::vector<Affine> moving_walls;
    std::vector<Affine> parameter_walls;
};

FeasibilitySummary assess_feasibility(
    const std::set<Affine, AffineLess>& walls
) {
    FeasibilitySummary summary;
    for (const Affine& wall : walls) {
        z3::context context;
        z3::solver solver(context);
        const z3::expr level = context.int_const("K");
        const z3::expr label = context.int_const("Q");
        const z3::expr crossing = context.int_const("V");
        const z3::expr target = context.int_const("x");
        solver.add(label >= 7);
        solver.add(level - 1 - 2 * label >= 10);
        solver.add(crossing >= 0 && 2 * crossing < level);
        solver.add(target >= 0 && 2 * target < level);
        solver.add(level - label - crossing <= label + crossing);
        const z3::expr expression = z3_affine(
            context,
            wall,
            level,
            label,
            crossing,
            target
        );
        const int crossing_coefficient = wall.coefficient[V];
        if (crossing_coefficient == 0) {
            ++summary.parameter;
            solver.add(expression == 0);
            if (solver.check() == z3::sat) {
                ++summary.feasible_parameter;
                summary.parameter_walls.push_back(wall);
            }
        } else {
            ++summary.moving;
            if (crossing_coefficient > 0) {
                solver.add(
                    expression <= 0
                    && expression >= -crossing_coefficient
                );
            } else {
                solver.add(
                    expression >= 0
                    && expression <= -crossing_coefficient
                );
            }
            if (solver.check() == z3::sat) {
                ++summary.feasible_moving;
                summary.moving_walls.push_back(wall);
            }
        }
    }
    return summary;
}

long long evaluate_wall(
    const Affine& wall,
    int level,
    int label,
    int crossing,
    int target
) {
    return static_cast<long long>(wall.coefficient[K]) * level
        + static_cast<long long>(wall.coefficient[Q]) * label
        + static_cast<long long>(wall.coefficient[V]) * crossing
        + static_cast<long long>(wall.coefficient[X]) * target
        + wall.coefficient[Constant];
}

bool has_moving_event(
    const std::set<Affine, AffineLess>& walls,
    int level,
    int label,
    int target,
    int lower,
    int upper
) {
    for (const Affine& wall : walls) {
        if (wall.coefficient[V] == 0) {
            continue;
        }
        long long previous = evaluate_wall(
            wall,
            level,
            label,
            lower,
            target
        );
        if (previous == 0) {
            return true;
        }
        for (int crossing = lower + 1; crossing <= upper; ++crossing) {
            const long long current = evaluate_wall(
                wall,
                level,
                label,
                crossing,
                target
            );
            if (
                current == 0
                || (previous < 0 && current > 0)
                || (previous > 0 && current < 0)
            ) {
                return true;
            }
            previous = current;
        }
    }
    return false;
}

bool hits_moving_wall(
    const std::set<Affine, AffineLess>& walls,
    int level,
    int label,
    int crossing,
    int target
) {
    for (const Affine& wall : walls) {
        if (
            wall.coefficient[V] != 0
            && evaluate_wall(wall, level, label, crossing, target) == 0
        ) {
            return true;
        }
    }
    return false;
}

bool crosses_moving_wall(
    const std::set<Affine, AffineLess>& walls,
    int level,
    int label,
    int left,
    int target
) {
    for (const Affine& wall : walls) {
        if (wall.coefficient[V] == 0) {
            continue;
        }
        const long long first = evaluate_wall(
            wall,
            level,
            label,
            left,
            target
        );
        const long long second = evaluate_wall(
            wall,
            level,
            label,
            left + 1,
            target
        );
        if (
            first == 0
            || second == 0
            || (first < 0 && second > 0)
            || (first > 0 && second < 0)
        ) {
            return true;
        }
    }
    return false;
}

void emit_segments(
    const std::set<Affine, AffineLess>& walls,
    int level,
    int label,
    int target
) {
    if (2 * label >= level) {
        throw std::runtime_error("segments require 2Q<K");
    }
    const int paired = (level + 1) / 2;
    if (target >= paired) {
        throw std::runtime_error("target lies outside the lower orbit");
    }
    std::vector<bool> singleton(static_cast<std::size_t>(paired));
    std::vector<bool> edge_cut(
        static_cast<std::size_t>(std::max(0, paired - 1))
    );
    for (int crossing = 0; crossing < paired; ++crossing) {
        singleton[static_cast<std::size_t>(crossing)] = hits_moving_wall(
            walls,
            level,
            label,
            crossing,
            target
        );
        if (crossing + 1 < paired) {
            edge_cut[static_cast<std::size_t>(crossing)] =
                crosses_moving_wall(
                    walls,
                    level,
                    label,
                    crossing,
                    target
                );
        }
    }
    std::uint64_t segment_count = 0U;
    std::uint64_t wall_vertex_count = 0U;
    std::uint64_t rail_count = 0U;
    for (int lower = 0; lower < paired;) {
        int upper = lower;
        if (singleton[static_cast<std::size_t>(lower)]) {
            ++wall_vertex_count;
        } else {
            while (
                upper + 1 < paired
                && !singleton[static_cast<std::size_t>(upper + 1)]
                && !edge_cut[static_cast<std::size_t>(upper)]
            ) {
                ++upper;
            }
        }
        ++segment_count;
        std::cout
            << "SEGMENT"
            << " V_lower=" << lower
            << " V_upper=" << upper
            << " singleton=" << (lower == upper ? 1 : 0)
            << '\n';
        for (int residue = 0; residue < 4; ++residue) {
            const int offset = (residue - lower % 4 + 4) % 4;
            const int rail_lower = lower + offset;
            if (rail_lower > upper) {
                continue;
            }
            const int rail_size = 1 + (upper - rail_lower) / 4;
            ++rail_count;
            std::cout
                << "RAIL"
                << " V_residue=" << residue
                << " V_lower=" << rail_lower
                << " V_upper=" << rail_lower + 4 * (rail_size - 1)
                << " terms=" << rail_size
                << '\n';
        }
        lower = upper + 1;
    }
    std::cout
        << "SU2_SHELL_V_SEGMENTS"
        << " K=" << level
        << " Q=" << label
        << " x=" << target
        << " lower_orbit_size=" << paired
        << " safe_segments=" << segment_count
        << " wall_vertices=" << wall_vertex_count
        << " rails=" << rail_count
        << " result=PASS_SAFE_WALL_PARTITION\n";
}

struct SegmentAuditSummary {
    std::uint64_t parameters = 0U;
    std::uint64_t targets = 0U;
    std::uint64_t candidate_windows = 0U;
    std::uint64_t wall_censored_windows = 0U;
    std::uint64_t cubic_windows = 0U;
    std::uint64_t direct_cubic_rails = 0U;
    std::uint64_t direct_cubic_pieces = 0U;
    std::uint64_t maximum_direct_piece_terms = 0U;
    std::uint64_t safe_rails = 0U;
    std::uint64_t negative_safe_rails = 0U;
    std::uint64_t free_safe_rails = 0U;
    std::uint64_t negative_free_safe_rails = 0U;
    std::uint64_t full_tails = 0U;
    std::uint64_t negative_full_tails = 0U;
    std::uint64_t free_full_tails = 0U;
    std::uint64_t negative_free_full_tails = 0U;
    std::uint64_t residue_tails = 0U;
    std::uint64_t negative_residue_tails = 0U;
    std::uint64_t free_residue_tails = 0U;
    std::uint64_t negative_free_residue_tails = 0U;
    std::uint64_t free_residue_negative_atoms = 0U;
    int maximum_free_residue_payment_span = 0;
    bool has_maximum_free_residue_payment_span = false;
    int maximum_free_residue_payment_level = -1;
    int maximum_free_residue_payment_label = -1;
    int maximum_free_residue_payment_target = -1;
    int maximum_free_residue_payment_residue = -1;
    int maximum_free_residue_payment_debit = -1;
    int maximum_free_residue_payment_credit = -1;
    bool has_negative_safe_rail = false;
    bool has_negative_free_safe_rail = false;
    int first_negative_level = -1;
    int first_negative_label = -1;
    int first_negative_target = -1;
    int first_negative_residue = -1;
    int first_negative_lower = -1;
    int first_negative_upper = -1;
    Integer first_negative_sum = 0;
    int first_negative_free_level = -1;
    int first_negative_free_label = -1;
    int first_negative_free_target = -1;
    int first_negative_free_residue = -1;
    int first_negative_free_lower = -1;
    int first_negative_free_upper = -1;
    Integer first_negative_free_sum = 0;
    bool has_negative_free_full_tail = false;
    int first_negative_free_full_level = -1;
    int first_negative_free_full_label = -1;
    int first_negative_free_full_target = -1;
    int first_negative_free_full_rho = -1;
    Integer first_negative_free_full_sum = 0;
    bool has_negative_free_residue_tail = false;
    int first_negative_free_residue_level = -1;
    int first_negative_free_residue_label = -1;
    int first_negative_free_residue_target = -1;
    int first_negative_free_residue_class = -1;
    int first_negative_free_residue_rho = -1;
    Integer first_negative_free_residue_sum = 0;
};

Integer third_difference(
    const std::vector<Integer>& values,
    std::size_t index
) {
    return values[index + 3U]
        - 3 * values[index + 2U]
        + 3 * values[index + 1U]
        - values[index];
}

void record_direct_cubic_partition(
    const std::vector<Integer>& summand,
    SegmentAuditSummary& summary
) {
    for (int residue = 0; residue < 4; ++residue) {
        std::vector<Integer> rail;
        for (std::size_t index = static_cast<std::size_t>(residue);
             index < summand.size();
             index += 4U) {
            rail.push_back(summand[index]);
        }
        if (rail.empty()) {
            continue;
        }
        ++summary.direct_cubic_rails;
        std::size_t lower = 0U;
        while (lower < rail.size()) {
            std::size_t upper = std::min(lower + 3U, rail.size() - 1U);
            if (lower + 3U < rail.size()) {
                const Integer base = third_difference(rail, lower);
                std::size_t next = lower + 1U;
                while (
                    next + 3U < rail.size()
                    && third_difference(rail, next) == base
                ) {
                    upper = next + 3U;
                    ++next;
                }
            }
            ++summary.direct_cubic_pieces;
            const std::uint64_t terms = static_cast<std::uint64_t>(
                upper - lower + 1U
            );
            summary.maximum_direct_piece_terms = std::max(
                summary.maximum_direct_piece_terms,
                terms
            );
            lower = upper + 1U;
        }
    }
}

void record_safe_rail_signs(
    const std::set<Affine, AffineLess>& walls,
    int level,
    int label,
    int target,
    const std::vector<Integer>& summand,
    SegmentAuditSummary& summary
) {
    const int paired = static_cast<int>(summand.size());
    std::vector<bool> singleton(static_cast<std::size_t>(paired));
    std::vector<bool> edge_cut(
        static_cast<std::size_t>(std::max(0, paired - 1))
    );
    for (int crossing = 0; crossing < paired; ++crossing) {
        singleton[static_cast<std::size_t>(crossing)] = hits_moving_wall(
            walls,
            level,
            label,
            crossing,
            target
        );
        if (crossing + 1 < paired) {
            edge_cut[static_cast<std::size_t>(crossing)] =
                crosses_moving_wall(
                    walls,
                    level,
                    label,
                    crossing,
                    target
                );
        }
    }
    for (int lower = 0; lower < paired;) {
        int upper = lower;
        if (!singleton[static_cast<std::size_t>(lower)]) {
            while (
                upper + 1 < paired
                && !singleton[static_cast<std::size_t>(upper + 1)]
                && !edge_cut[static_cast<std::size_t>(upper)]
            ) {
                ++upper;
            }
        }
        for (int residue = 0; residue < 4; ++residue) {
            const int offset = (residue - lower % 4 + 4) % 4;
            const int rail_lower = lower + offset;
            if (rail_lower > upper) {
                continue;
            }
            Integer rail_sum = 0;
            for (int crossing = rail_lower;
                 crossing <= upper;
                 crossing += 4) {
                rail_sum += summand[static_cast<std::size_t>(crossing)];
            }
            ++summary.safe_rails;
            if (rail_sum < 0) {
                ++summary.negative_safe_rails;
                if (!summary.has_negative_safe_rail) {
                    summary.has_negative_safe_rail = true;
                    summary.first_negative_level = level;
                    summary.first_negative_label = label;
                    summary.first_negative_target = target;
                    summary.first_negative_residue = residue;
                    summary.first_negative_lower = rail_lower;
                    summary.first_negative_upper = upper;
                    summary.first_negative_sum = rail_sum;
                }
            }
            if (label >= 7 && level - 1 - 2 * label >= 10) {
                ++summary.free_safe_rails;
                if (rail_sum < 0) {
                    ++summary.negative_free_safe_rails;
                    if (!summary.has_negative_free_safe_rail) {
                        summary.has_negative_free_safe_rail = true;
                        summary.first_negative_free_level = level;
                        summary.first_negative_free_label = label;
                        summary.first_negative_free_target = target;
                        summary.first_negative_free_residue = residue;
                        summary.first_negative_free_lower = rail_lower;
                        summary.first_negative_free_upper = upper;
                        summary.first_negative_free_sum = rail_sum;
                    }
                }
            }
        }
        lower = upper + 1;
    }
}

void record_tail_signs(
    int level,
    int label,
    int target,
    const std::vector<Integer>& summand,
    SegmentAuditSummary& summary
) {
    const int paired = static_cast<int>(summand.size());
    const bool free_band = label >= 7 && level - 1 - 2 * label >= 10;
    Integer full_tail = 0;
    for (int crossing = paired - 1; crossing >= 0; --crossing) {
        full_tail += summand[static_cast<std::size_t>(crossing)];
        ++summary.full_tails;
        if (full_tail < 0) {
            ++summary.negative_full_tails;
        }
        if (free_band) {
            ++summary.free_full_tails;
            if (full_tail < 0) {
                ++summary.negative_free_full_tails;
                if (!summary.has_negative_free_full_tail) {
                    summary.has_negative_free_full_tail = true;
                    summary.first_negative_free_full_level = level;
                    summary.first_negative_free_full_label = label;
                    summary.first_negative_free_full_target = target;
                    summary.first_negative_free_full_rho = crossing;
                    summary.first_negative_free_full_sum = full_tail;
                }
            }
        }
    }
    for (int residue = 0; residue < 4; ++residue) {
        int crossing = paired - 1;
        while (crossing >= 0 && crossing % 4 != residue) {
            --crossing;
        }
        Integer residue_tail = 0;
        for (; crossing >= 0; crossing -= 4) {
            residue_tail += summand[static_cast<std::size_t>(crossing)];
            ++summary.residue_tails;
            if (residue_tail < 0) {
                ++summary.negative_residue_tails;
            }
            if (free_band) {
                ++summary.free_residue_tails;
                if (residue_tail < 0) {
                    ++summary.negative_free_residue_tails;
                    if (!summary.has_negative_free_residue_tail) {
                        summary.has_negative_free_residue_tail = true;
                        summary.first_negative_free_residue_level = level;
                        summary.first_negative_free_residue_label = label;
                        summary.first_negative_free_residue_target = target;
                        summary.first_negative_free_residue_class = residue;
                        summary.first_negative_free_residue_rho = crossing;
                        summary.first_negative_free_residue_sum = residue_tail;
                    }
                }
            }
        }
        if (free_band) {
            struct Reserve {
                int crossing = 0;
                Integer amount = 0;
            };
            std::deque<Reserve> reserve;
            for (int payment_crossing = paired - 1;
                 payment_crossing >= 0;
                 --payment_crossing) {
                if (payment_crossing % 4 != residue) {
                    continue;
                }
                const Integer& value =
                    summand[static_cast<std::size_t>(payment_crossing)];
                if (value >= 0) {
                    if (value != 0) {
                        reserve.push_back(Reserve{payment_crossing, value});
                    }
                    continue;
                }
                ++summary.free_residue_negative_atoms;
                Integer debt = -value;
                while (debt != 0) {
                    if (reserve.empty()) {
                        throw std::runtime_error(
                            "negative residue tail has no later reserve"
                        );
                    }
                    Reserve& nearest = reserve.back();
                    const Integer payment = std::min(nearest.amount, debt);
                    nearest.amount -= payment;
                    debt -= payment;
                    const int span = nearest.crossing - payment_crossing;
                    if (span > summary.maximum_free_residue_payment_span) {
                        summary.maximum_free_residue_payment_span = span;
                        summary.has_maximum_free_residue_payment_span = true;
                        summary.maximum_free_residue_payment_level = level;
                        summary.maximum_free_residue_payment_label = label;
                        summary.maximum_free_residue_payment_target = target;
                        summary.maximum_free_residue_payment_residue = residue;
                        summary.maximum_free_residue_payment_debit =
                            payment_crossing;
                        summary.maximum_free_residue_payment_credit =
                            nearest.crossing;
                    }
                    if (nearest.amount == 0) {
                        reserve.pop_back();
                    }
                }
            }
        }
    }
}

SegmentAuditSummary audit_cubic_segments(
    const std::set<Affine, AffineLess>& walls,
    int maximum_level
) {
    SegmentAuditSummary summary;
    for (int level = 3; level <= maximum_level; ++level) {
        const int paired = (level + 1) / 2;
        const bool has_center = level % 2 == 0;
        const int quotient_size = paired + (has_center ? 1 : 0);
        const int center = has_center ? paired : -1;
        for (int label = 1; 2 * label < level; ++label) {
            ++summary.parameters;
            Matrix plus(
                static_cast<std::size_t>(quotient_size),
                std::vector<int>(static_cast<std::size_t>(quotient_size))
            );
            Matrix minus(
                static_cast<std::size_t>(paired),
                std::vector<int>(static_cast<std::size_t>(paired))
            );
            Matrix delta(
                static_cast<std::size_t>(quotient_size),
                std::vector<int>(static_cast<std::size_t>(paired))
            );
            for (int source = 0; source < paired; ++source) {
                for (int target = 0; target < paired; ++target) {
                    const int same = fuses_half(
                        level,
                        label,
                        source,
                        target
                    ) ? 1 : 0;
                    const int crossed = fuses_half(
                        level,
                        label,
                        source,
                        level - target
                    ) ? 1 : 0;
                    plus[static_cast<std::size_t>(source)]
                        [static_cast<std::size_t>(target)] = same + crossed;
                    minus[static_cast<std::size_t>(source)]
                         [static_cast<std::size_t>(target)] = same - crossed;
                    delta[static_cast<std::size_t>(source)]
                         [static_cast<std::size_t>(target)] = 2 * crossed;
                }
                if (has_center) {
                    const int joins = fuses_half(
                        level,
                        label,
                        source,
                        level / 2
                    ) ? 1 : 0;
                    plus[static_cast<std::size_t>(source)]
                        [static_cast<std::size_t>(center)] = joins;
                    plus[static_cast<std::size_t>(center)]
                        [static_cast<std::size_t>(source)] = 2 * joins;
                    delta[static_cast<std::size_t>(center)]
                         [static_cast<std::size_t>(source)] = 2 * joins;
                }
            }
            if (has_center) {
                plus[static_cast<std::size_t>(center)]
                    [static_cast<std::size_t>(center)] = 1;
            }

            std::vector<std::vector<Integer>> plus_rows(
                6U,
                std::vector<Integer>(static_cast<std::size_t>(quotient_size))
            );
            std::vector<std::vector<Integer>> minus_rows(
                6U,
                std::vector<Integer>(static_cast<std::size_t>(paired))
            );
            plus_rows[0][0] = 1;
            minus_rows[0][0] = 1;
            for (int power = 1; power <= 5; ++power) {
                plus_rows[static_cast<std::size_t>(power)] = multiply_row(
                    plus_rows[static_cast<std::size_t>(power - 1)],
                    plus
                );
                minus_rows[static_cast<std::size_t>(power)] = multiply_row(
                    minus_rows[static_cast<std::size_t>(power - 1)],
                    minus
                );
            }
            const Integer f4 = (
                plus_rows[4][0] + minus_rows[4][0]
            ) / 2;
            const Integer f5 = (
                plus_rows[5][0] + minus_rows[5][0]
            ) / 2;

            std::vector<std::vector<Integer>> crossing_weights(
                6U,
                std::vector<Integer>(static_cast<std::size_t>(paired))
            );
            for (int power = 1; power <= 5; ++power) {
                for (int crossing = 0; crossing < paired; ++crossing) {
                    Integer value = 0;
                    for (int source = 0; source < quotient_size; ++source) {
                        value += plus_rows[static_cast<std::size_t>(power)]
                            [static_cast<std::size_t>(source)]
                            * delta[static_cast<std::size_t>(source)]
                                   [static_cast<std::size_t>(crossing)];
                    }
                    crossing_weights[static_cast<std::size_t>(power)]
                                    [static_cast<std::size_t>(crossing)] = value;
                }
            }

            std::vector<IntegerMatrix> terminal_powers(
                5U,
                IntegerMatrix(
                    static_cast<std::size_t>(paired),
                    std::vector<Integer>(static_cast<std::size_t>(paired))
                )
            );
            terminal_powers[0] = identity_matrix(paired);
            for (int power = 1; power <= 4; ++power) {
                terminal_powers[static_cast<std::size_t>(power)] =
                    multiply_right(
                        terminal_powers[static_cast<std::size_t>(power - 1)],
                        minus
                    );
            }

            for (int target = 0; target < paired; ++target) {
                ++summary.targets;
                std::vector<Integer> summand(
                    static_cast<std::size_t>(paired)
                );
                for (int crossing = 0; crossing < paired; ++crossing) {
                    const std::size_t coordinate =
                        static_cast<std::size_t>(crossing);
                    const std::size_t target_coordinate =
                        static_cast<std::size_t>(target);
                    const Integer even = f4 * (
                        crossing_weights[1][coordinate]
                            * terminal_powers[4][coordinate][target_coordinate]
                        + crossing_weights[2][coordinate]
                            * terminal_powers[3][coordinate][target_coordinate]
                        + crossing_weights[3][coordinate]
                            * terminal_powers[2][coordinate][target_coordinate]
                        + crossing_weights[4][coordinate]
                            * terminal_powers[1][coordinate][target_coordinate]
                        + crossing_weights[5][coordinate]
                            * terminal_powers[0][coordinate][target_coordinate]
                    );
                    const Integer odd = f5 * (
                        crossing_weights[1][coordinate]
                            * terminal_powers[3][coordinate][target_coordinate]
                        + crossing_weights[2][coordinate]
                            * terminal_powers[2][coordinate][target_coordinate]
                        + crossing_weights[3][coordinate]
                            * terminal_powers[1][coordinate][target_coordinate]
                        + crossing_weights[4][coordinate]
                            * terminal_powers[0][coordinate][target_coordinate]
                    );
                    summand[coordinate] = even - odd;
                }
                record_direct_cubic_partition(summand, summary);
                record_safe_rail_signs(
                    walls,
                    level,
                    label,
                    target,
                    summand,
                    summary
                );
                record_tail_signs(
                    level,
                    label,
                    target,
                    summand,
                    summary
                );

                for (int lower = 0; lower + 16 < paired; ++lower) {
                    ++summary.candidate_windows;
                    if (has_moving_event(
                            walls,
                            level,
                            label,
                            target,
                            lower,
                            lower + 16
                        )) {
                        ++summary.wall_censored_windows;
                        continue;
                    }
                    const Integer fourth_difference =
                        summand[static_cast<std::size_t>(lower + 16)]
                        - 4 * summand[static_cast<std::size_t>(lower + 12)]
                        + 6 * summand[static_cast<std::size_t>(lower + 8)]
                        - 4 * summand[static_cast<std::size_t>(lower + 4)]
                        + summand[static_cast<std::size_t>(lower)];
                    if (fourth_difference != 0) {
                        throw std::runtime_error(
                            "noncubic segment witness K="
                            + std::to_string(level)
                            + " Q=" + std::to_string(label)
                            + " x=" + std::to_string(target)
                            + " V=" + std::to_string(lower)
                        );
                    }
                    ++summary.cubic_windows;
                }
            }
        }
    }
    return summary;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const bool list = argc == 2 && std::string(argv[1]) == "--list";
        const bool feasible =
            argc == 2 && std::string(argv[1]) == "--feasible";
        const bool feasible_list =
            argc == 2 && std::string(argv[1]) == "--feasible-list";
        const bool audit = argc == 3 && std::string(argv[1]) == "--audit";
        const bool segments =
            argc == 5 && std::string(argv[1]) == "--segments";
        if (
            argc != 1
            && !list
            && !feasible
            && !feasible_list
            && !audit
            && !segments
        ) {
            throw std::runtime_error(
                "usage: [--list|--feasible|--feasible-list|--audit MAXIMUM_K|"
                "--segments K Q x]"
            );
        }
        std::set<Affine, AffineLess> crossing_walls;
        std::set<Affine, AffineLess> terminal_walls;
        std::set<Affine, AffineLess> return_walls;
        append_crossing_walls(crossing_walls);
        append_terminal_walls(terminal_walls);
        append_return_walls(return_walls);
        std::set<Affine, AffineLess> joint_walls = crossing_walls;
        joint_walls.insert(terminal_walls.begin(), terminal_walls.end());
        joint_walls.insert(return_walls.begin(), return_walls.end());
        std::cout
            << "SU2_SHELL_V_WALL_CATALOG"
            << " crossing_walls=" << crossing_walls.size()
            << " crossing_v_walls=" << count_v_walls(crossing_walls)
            << " terminal_walls=" << terminal_walls.size()
            << " terminal_v_walls=" << count_v_walls(terminal_walls)
            << " return_walls=" << return_walls.size()
            << " joint_walls=" << joint_walls.size()
            << " joint_v_walls=" << count_v_walls(joint_walls)
            << " v_coefficients=" << v_coefficient_set(joint_walls)
            << " result=PASS_SYMBOLIC_CATALOG\n";
        if (feasible || feasible_list) {
            const FeasibilitySummary summary = assess_feasibility(joint_walls);
            std::cout
                << "SU2_SHELL_V_WALL_FEASIBILITY"
                << " moving=" << summary.moving
                << " feasible_moving=" << summary.feasible_moving
                << " parameter=" << summary.parameter
                << " feasible_parameter=" << summary.feasible_parameter
                << " result=PASS_EXACT_PRESBURGER_FILTER\n";
            if (feasible_list) {
                for (const Affine& wall : summary.moving_walls) {
                    std::cout << "MOVING ";
                    print_wall(wall);
                }
                for (const Affine& wall : summary.parameter_walls) {
                    std::cout << "PARAMETER ";
                    print_wall(wall);
                }
            }
        }
        if (list) {
            for (const Affine& wall : joint_walls) {
                print_wall(wall);
            }
        }
        if (audit) {
            const SegmentAuditSummary summary = audit_cubic_segments(
                joint_walls,
                parse_positive(argv[2])
            );
            std::cout
                << "SU2_SHELL_V_SEGMENT_AUDIT"
                << " maximum_K=" << argv[2]
                << " parameters=" << summary.parameters
                << " targets=" << summary.targets
                << " candidate_windows=" << summary.candidate_windows
                << " wall_censored_windows=" << summary.wall_censored_windows
                << " cubic_windows=" << summary.cubic_windows
                << " direct_cubic_rails=" << summary.direct_cubic_rails
                << " direct_cubic_pieces=" << summary.direct_cubic_pieces
                << " maximum_direct_piece_terms="
                    << summary.maximum_direct_piece_terms
                << " safe_rails=" << summary.safe_rails
                << " negative_safe_rails=" << summary.negative_safe_rails
                << " free_safe_rails=" << summary.free_safe_rails
                << " negative_free_safe_rails="
                    << summary.negative_free_safe_rails
                << " full_tails=" << summary.full_tails
                << " negative_full_tails=" << summary.negative_full_tails
                << " free_full_tails=" << summary.free_full_tails
                << " negative_free_full_tails="
                    << summary.negative_free_full_tails
                << " residue_tails=" << summary.residue_tails
                << " negative_residue_tails="
                    << summary.negative_residue_tails
                << " free_residue_tails=" << summary.free_residue_tails
                << " negative_free_residue_tails="
                    << summary.negative_free_residue_tails
                << " free_residue_negative_atoms="
                    << summary.free_residue_negative_atoms
                << " maximum_free_residue_payment_span="
                    << summary.maximum_free_residue_payment_span
                << " maximum_free_residue_payment_witness=";
            if (summary.has_maximum_free_residue_payment_span) {
                std::cout
                    << "K:" << summary.maximum_free_residue_payment_level
                    << ",Q:" << summary.maximum_free_residue_payment_label
                    << ",x:" << summary.maximum_free_residue_payment_target
                    << ",residue:"
                    << summary.maximum_free_residue_payment_residue
                    << ",debit:" << summary.maximum_free_residue_payment_debit
                    << ",credit:"
                    << summary.maximum_free_residue_payment_credit;
            } else {
                std::cout << "none";
            }
            std::cout
                << " first_negative_safe_rail=";
            if (summary.has_negative_safe_rail) {
                std::cout
                    << "K:" << summary.first_negative_level
                    << ",Q:" << summary.first_negative_label
                    << ",x:" << summary.first_negative_target
                    << ",residue:" << summary.first_negative_residue
                    << ",lower:" << summary.first_negative_lower
                    << ",upper:" << summary.first_negative_upper
                    << ",sum:" << summary.first_negative_sum;
            } else {
                std::cout << "none";
            }
            std::cout << " first_negative_free_safe_rail=";
            if (summary.has_negative_free_safe_rail) {
                std::cout
                    << "K:" << summary.first_negative_free_level
                    << ",Q:" << summary.first_negative_free_label
                    << ",x:" << summary.first_negative_free_target
                    << ",residue:" << summary.first_negative_free_residue
                    << ",lower:" << summary.first_negative_free_lower
                    << ",upper:" << summary.first_negative_free_upper
                    << ",sum:" << summary.first_negative_free_sum;
            } else {
                std::cout << "none";
            }
            std::cout << " first_negative_free_full_tail=";
            if (summary.has_negative_free_full_tail) {
                std::cout
                    << "K:" << summary.first_negative_free_full_level
                    << ",Q:" << summary.first_negative_free_full_label
                    << ",x:" << summary.first_negative_free_full_target
                    << ",rho:" << summary.first_negative_free_full_rho
                    << ",sum:" << summary.first_negative_free_full_sum;
            } else {
                std::cout << "none";
            }
            std::cout << " first_negative_free_residue_tail=";
            if (summary.has_negative_free_residue_tail) {
                std::cout
                    << "K:" << summary.first_negative_free_residue_level
                    << ",Q:" << summary.first_negative_free_residue_label
                    << ",x:" << summary.first_negative_free_residue_target
                    << ",residue:"
                    << summary.first_negative_free_residue_class
                    << ",rho:" << summary.first_negative_free_residue_rho
                    << ",sum:" << summary.first_negative_free_residue_sum;
            } else {
                std::cout << "none";
            }
            std::cout
                << " residue_step=4"
                << " result=PASS_CUBIC_SEGMENT_AUDIT\n";
        }
        if (segments) {
            emit_segments(
                joint_walls,
                parse_positive(argv[2]),
                parse_positive(argv[3]),
                parse_nonnegative(argv[4])
            );
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
