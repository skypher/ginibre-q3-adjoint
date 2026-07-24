#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

struct LeftState {
    unsigned active_mask;
    std::vector<int> scalar_grades;
    std::vector<int> fusion_heights;
};

struct SearchStats {
    std::size_t configurations = 0;
    std::size_t left_states = 0;
};

void enumerate_fusion_paths(const std::vector<int>& plus_labels,
                            unsigned active_mask,
                            int factor,
                            int height,
                            std::vector<int>& fusion_heights,
                            std::vector<std::vector<int>>& paths) {
    if (factor == static_cast<int>(plus_labels.size())) {
        if (height == 0) {
            paths.push_back(fusion_heights);
        }
        return;
    }

    if ((active_mask & (1U << factor)) == 0U) {
        fusion_heights[factor] = -1;
        enumerate_fusion_paths(plus_labels, active_mask, factor + 1, height,
                               fusion_heights, paths);
        return;
    }

    const int label = plus_labels[factor];
    for (int next = std::abs(height - label); next <= height + label;
         next += 2) {
        fusion_heights[factor] = next;
        enumerate_fusion_paths(plus_labels, active_mask, factor + 1, next,
                               fusion_heights, paths);
    }
}

void print_vector(const std::vector<int>& values) {
    for (const int value : values) {
        std::cout << value << ',';
    }
}

void print_left_state(const LeftState& state) {
    std::cout << " mask=" << state.active_mask << " scalar_grades=";
    print_vector(state.scalar_grades);
    std::cout << " fusion_heights=";
    print_vector(state.fusion_heights);
    std::cout << '\n';
}

bool enumerate_scalar_grades(
    const std::vector<int>& plus_labels,
    unsigned active_mask,
    int factor,
    std::vector<int>& scalar_grades,
    const std::vector<int>& fusion_heights,
    int q,
    std::map<std::vector<int>, LeftState>& seen_targets,
    SearchStats& stats) {
    if (factor < static_cast<int>(plus_labels.size())) {
        if ((active_mask & (1U << factor)) != 0U) {
            scalar_grades[factor] = -1;
            return enumerate_scalar_grades(
                plus_labels, active_mask, factor + 1, scalar_grades,
                fusion_heights, q, seen_targets, stats);
        }

        for (int grade = 0; grade <= plus_labels[factor]; ++grade) {
            scalar_grades[factor] = grade;
            if (enumerate_scalar_grades(
                    plus_labels, active_mask, factor + 1, scalar_grades,
                    fusion_heights, q, seen_targets, stats)) {
                return true;
            }
        }
        return false;
    }

    int active_label_sum = q;
    int inactive_grade_sum = 0;
    for (int i = 0; i < static_cast<int>(plus_labels.size()); ++i) {
        if ((active_mask & (1U << i)) != 0U) {
            active_label_sum += plus_labels[i];
        } else {
            inactive_grade_sum += scalar_grades[i];
        }
    }
    if (active_label_sum % 2 != 0) {
        return false;
    }
    const int scalar_index_grade = active_label_sum / 2 + inactive_grade_sum;

    int previous_height = q;
    int first_descent = -1;
    int carrier_label = -1;
    std::vector<int> output_status(plus_labels.size(), -3);
    std::vector<int> suffix_heights(plus_labels.size(), -3);

    for (int i = 0; i < static_cast<int>(plus_labels.size()); ++i) {
        if ((active_mask & (1U << i)) == 0U) {
            output_status[i] = scalar_grades[i];
            continue;
        }

        const int height = fusion_heights[i];
        if (first_descent < 0 && height < q) {
            first_descent = i;
            carrier_label = height;
        }

        if (first_descent < 0 || i == first_descent) {
            output_status[i] =
                (plus_labels[i] + height - previous_height) / 2;
        } else {
            output_status[i] = -1;
            suffix_heights[i] = height;
        }
        previous_height = height;
    }

    if (first_descent < 0) {
        return false;
    }

    std::vector<int> target_key{scalar_index_grade, carrier_label};
    target_key.insert(target_key.end(), output_status.begin(),
                      output_status.end());
    target_key.insert(target_key.end(), suffix_heights.begin(),
                      suffix_heights.end());

    ++stats.left_states;
    const LeftState state{active_mask, scalar_grades, fusion_heights};
    const auto [iterator, inserted] =
        seen_targets.emplace(target_key, state);
    if (inserted) {
        return false;
    }

    std::cout << "H2_FIRST_DESCENT_COLLISION\n";
    std::cout << "q=" << q << " plus_labels=";
    print_vector(plus_labels);
    std::cout << " scalar_index_grade=" << scalar_index_grade
              << " carrier_label=" << carrier_label << '\n';
    print_left_state(iterator->second);
    print_left_state(state);
    std::cout << " target_status=";
    print_vector(output_status);
    std::cout << " suffix_heights=";
    print_vector(suffix_heights);
    std::cout << '\n';
    std::cout << "neighborhood_size=1 left_subset_size=2 hall_deficit=1\n";
    std::cout << "first_collision_factor_count=" << plus_labels.size()
              << '\n';
    std::cout << "minimality_order=factor_count,ordered_plus_word,q\n";
    std::cout << "configurations_checked=" << stats.configurations
              << " left_states_checked=" << stats.left_states << '\n';
    return true;
}

bool has_collision(int q,
                   const std::vector<int>& plus_labels,
                   SearchStats& stats) {
    ++stats.configurations;
    std::map<std::vector<int>, LeftState> seen_targets;
    const unsigned mask_limit = 1U << plus_labels.size();

    for (unsigned active_mask = 1; active_mask < mask_limit; ++active_mask) {
        std::vector<int> fusion_heights(plus_labels.size(), -1);
        std::vector<std::vector<int>> paths;
        enumerate_fusion_paths(plus_labels, active_mask, 0, q, fusion_heights,
                               paths);

        for (const std::vector<int>& path : paths) {
            std::vector<int> scalar_grades(plus_labels.size(), -1);
            if (enumerate_scalar_grades(
                    plus_labels, active_mask, 0, scalar_grades, path, q,
                    seen_targets, stats)) {
                return true;
            }
        }
    }
    return false;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: search_su2_fdc_h2_first_descent "
                     "MAXIMUM_Q MAXIMUM_LABEL MAXIMUM_FACTORS\n";
        return 1;
    }

    const int maximum_q = std::stoi(argv[1]);
    const int maximum_label = std::stoi(argv[2]);
    const int maximum_factors = std::stoi(argv[3]);
    if (maximum_q < 1 || maximum_label < 1 || maximum_factors < 1 ||
        maximum_factors >= 8 * static_cast<int>(sizeof(unsigned))) {
        std::cerr << "invalid bound\n";
        return 1;
    }

    SearchStats stats;
    for (int factor_count = 1; factor_count <= maximum_factors;
         ++factor_count) {
        std::vector<int> plus_labels(factor_count, 1);
        while (true) {
            for (int q = 1; q <= maximum_q; ++q) {
                const bool support_disjoint =
                    std::find(plus_labels.begin(), plus_labels.end(), q) ==
                    plus_labels.end();
                if (support_disjoint &&
                    has_collision(q, plus_labels, stats)) {
                    return 0;
                }
            }

            int position = factor_count - 1;
            while (position >= 0 &&
                   plus_labels[position] == maximum_label) {
                --position;
            }
            if (position < 0) {
                break;
            }
            ++plus_labels[position];
            for (int i = position + 1; i < factor_count; ++i) {
                plus_labels[i] = 1;
            }
        }
    }

    std::cout << "H2_FIRST_DESCENT_SEARCH_PASS no_collision\n";
    std::cout << "configurations_checked=" << stats.configurations
              << " left_states_checked=" << stats.left_states << '\n';
    return 0;
}
