#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using Wedge = std::pair<int, int>;
using State = std::map<Wedge, int>;

void add_wedge(State& state, int first, int second, int coefficient) {
    if (first == second || coefficient == 0) {
        return;
    }
    if (first > second) {
        state[{first, second}] += coefficient;
    } else {
        state[{second, first}] -= coefficient;
    }
}

std::vector<int> ordinary_fusion(int first, int second) {
    std::vector<int> outputs;
    for (int label = std::abs(first - second); label <= first + second;
         label += 2) {
        outputs.push_back(label);
    }
    return outputs;
}

std::vector<int> level_fusion(int level, int first, int second) {
    std::vector<int> outputs;
    const int upper =
        std::min(first + second, 2 * level - first - second);
    for (int label = std::abs(first - second); label <= upper; label += 2) {
        outputs.push_back(label);
    }
    return outputs;
}

State ordinary_update(int label, int larger, int smaller) {
    State result;
    for (const int output : ordinary_fusion(label, larger)) {
        add_wedge(result, output, smaller, 1);
    }
    for (const int output : ordinary_fusion(label, smaller)) {
        add_wedge(result, larger, output, 1);
    }
    return result;
}

State ordinary_update(int label, const State& state) {
    State result;
    for (const auto& [wedge, coefficient] : state) {
        for (const int output : ordinary_fusion(label, wedge.first)) {
            add_wedge(result, output, wedge.second, coefficient);
        }
        for (const int output : ordinary_fusion(label, wedge.second)) {
            add_wedge(result, wedge.first, output, coefficient);
        }
    }
    return result;
}

State level_update(int level, int label, int larger, int smaller) {
    State result;
    for (const int output : level_fusion(level, label, larger)) {
        add_wedge(result, output, smaller, 1);
    }
    for (const int output : level_fusion(level, label, smaller)) {
        add_wedge(result, larger, output, 1);
    }
    return result;
}

void require_nonnegative(const State& state) {
    for (const auto& [wedge, coefficient] : state) {
        (void)wedge;
        if (coefficient < 0) {
            throw std::runtime_error("negative oriented-wedge coefficient");
        }
    }
}

void require_same_parity(const State& state) {
    require_nonnegative(state);
    for (const auto& [wedge, coefficient] : state) {
        if (coefficient != 0 && (wedge.first - wedge.second) % 2 != 0) {
            throw std::runtime_error("parity bridge left fixed-parity cone");
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: verify_su2_label2_wedge_cone "
                "MAXIMUM_ORDINARY_LABEL MAXIMUM_LEVEL");
        }
        const int maximum_ordinary_label = std::stoi(argv[1]);
        const int maximum_level = std::stoi(argv[2]);
        if (maximum_ordinary_label < 1 || maximum_level < 2) {
            throw std::runtime_error("invalid bound");
        }

        long long ordinary_columns = 0;
        long long fundamental_columns = 0;
        for (int larger = 1; larger <= maximum_ordinary_label; ++larger) {
            for (int smaller = 0; smaller < larger; ++smaller) {
                require_nonnegative(ordinary_update(1, larger, smaller));
                ++fundamental_columns;
                if ((larger - smaller) % 2 != 0) {
                    continue;
                }
                require_nonnegative(ordinary_update(2, larger, smaller));
                ++ordinary_columns;
            }
        }

        long long level_columns = 0;
        long long level_fundamental_columns = 0;
        long long parity_bridges = 0;
        for (int level = 2; level <= maximum_level; ++level) {
            for (int larger = 1; larger <= level; ++larger) {
                for (int smaller = 0; smaller < larger; ++smaller) {
                    require_nonnegative(
                        level_update(level, 1, larger, smaller));
                    ++level_fundamental_columns;
                    if ((larger - smaller) % 2 != 0) {
                        continue;
                    }
                    require_nonnegative(
                        level_update(level, 2, larger, smaller));
                    ++level_columns;
                }
            }
            for (int odd = 1; odd <= level; odd += 2) {
                require_same_parity(level_update(level, 1, odd, 0));
                ++parity_bridges;
            }
        }

        for (int label = 4; label <= maximum_ordinary_label; label += 2) {
            const State witness = ordinary_update(label, label - 2, 0);
            const auto found = witness.find({label, label - 2});
            if (found == witness.end() || found->second != -1) {
                throw std::runtime_error("missing higher-label obstruction");
            }
        }

        int grouped_power_obstructions = 0;
        for (int label = 4;
             label <= std::min(maximum_ordinary_label, 12); label += 2) {
            for (int power = 1; power <= 6; ++power) {
                bool found_negative = false;
                for (int larger = 2;
                     larger <= 2 * label && !found_negative; larger += 2) {
                    for (int smaller = 0;
                         smaller < larger && !found_negative; smaller += 2) {
                        State state;
                        add_wedge(state, larger, smaller, 1);
                        for (int step = 0; step < power; ++step) {
                            state = ordinary_update(label, state);
                        }
                        for (const auto& [wedge, coefficient] : state) {
                            (void)wedge;
                            if (coefficient < 0) {
                                found_negative = true;
                                ++grouped_power_obstructions;
                                break;
                            }
                        }
                    }
                }
                if (!found_negative) {
                    throw std::runtime_error(
                        "grouped power unexpectedly preserves tested cone");
                }
            }
        }

        std::cout << "SU2_LABEL2_WEDGE_CONE PASS\n";
        std::cout << "ordinary_columns=" << ordinary_columns
                  << " maximum_ordinary_label=" << maximum_ordinary_label
                  << '\n';
        std::cout << "ordinary_fundamental_columns=" << fundamental_columns
                  << '\n';
        std::cout << "finite_level_columns=" << level_columns
                  << " maximum_level=" << maximum_level << '\n';
        std::cout << "finite_level_fundamental_columns="
                  << level_fundamental_columns
                  << " odd_seed_parity_bridges=" << parity_bridges << '\n';
        std::cout << "higher_even_label_obstruction="
                     "A_p(e_(p-2) wedge e_0) has coefficient -1 at "
                     "e_p wedge e_(p-2), p>=4\n";
        std::cout << "bounded_grouped_power_obstructions="
                  << grouped_power_obstructions
                  << " labels=4,6,...,12 powers=1,...,6\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
