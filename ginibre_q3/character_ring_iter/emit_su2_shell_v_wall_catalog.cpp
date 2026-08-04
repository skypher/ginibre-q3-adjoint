#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <z3++.h>

namespace {

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

}  // namespace

int main(int argc, char** argv) {
    try {
        const bool list = argc == 2 && std::string(argv[1]) == "--list";
        const bool feasible =
            argc == 2 && std::string(argv[1]) == "--feasible";
        const bool feasible_list =
            argc == 2 && std::string(argv[1]) == "--feasible-list";
        if (argc != 1 && !list && !feasible && !feasible_list) {
            throw std::runtime_error("usage: [--list|--feasible|--feasible-list]");
        }
        std::set<Affine, AffineLess> crossing_walls;
        std::set<Affine, AffineLess> terminal_walls;
        append_crossing_walls(crossing_walls);
        append_terminal_walls(terminal_walls);
        std::set<Affine, AffineLess> joint_walls = crossing_walls;
        joint_walls.insert(terminal_walls.begin(), terminal_walls.end());
        std::cout
            << "SU2_SHELL_V_WALL_CATALOG"
            << " crossing_walls=" << crossing_walls.size()
            << " crossing_v_walls=" << count_v_walls(crossing_walls)
            << " terminal_walls=" << terminal_walls.size()
            << " terminal_v_walls=" << count_v_walls(terminal_walls)
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
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
