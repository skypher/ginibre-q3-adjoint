#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <z3++.h>

namespace {

using Exponents = std::array<z3::expr, 4U>;

struct Factor {
    z3::expr active;
    Exponents exponents;
};

struct Source {
    z3::expr active;
    Exponents q_exponents;
    Exponents r_exponents;
};

Exponents zero_exponents(z3::context& context) {
    return {
        context.int_val(0), context.int_val(0),
        context.int_val(0), context.int_val(0)
    };
}

Factor selected_factor(
    z3::context& context,
    const z3::expr& terminal,
    const z3::expr& label,
    std::size_t root
) {
    const z3::expr lower = (label == terminal - 1).simplify();
    const z3::expr upper = (label == terminal + 1).simplify();
    Exponents exponents = zero_exponents(context);
    exponents[0U] = z3::ite(
        lower, context.int_val(1), context.int_val(0)
    ).simplify();
    exponents[root] = z3::ite(lower, terminal - 1, terminal).simplify();
    return {(lower || upper).simplify(), std::move(exponents)};
}

Factor complement_factor(
    z3::context& context,
    const z3::expr& terminal,
    const z3::expr& lower_label,
    const z3::expr& upper_label,
    std::size_t lower_root,
    std::size_t upper_root,
    bool upper_channel
) {
    Exponents exponents = zero_exponents(context);
    const z3::expr gap = upper_label - lower_label;
    if (!upper_channel) {
        const z3::expr active
            = gap <= terminal - 1
            && lower_label + upper_label >= terminal - 1;
        exponents[0U] = context.int_val(1);
        exponents[lower_root]
            = ((terminal - 1 + lower_label - upper_label) / 2).simplify();
        exponents[upper_root]
            = ((terminal - 1 - lower_label + upper_label) / 2).simplify();
        return {active.simplify(), std::move(exponents)};
    }
    const z3::expr interior
        = gap <= terminal - 1
        && lower_label + upper_label >= terminal + 1;
    const z3::expr wall = gap == terminal + 1;
    exponents[lower_root] = z3::ite(
        wall, context.int_val(0),
        ((terminal + lower_label - upper_label - 1) / 2).simplify()
    );
    exponents[upper_root] = z3::ite(
        wall, terminal,
        ((terminal - lower_label + upper_label + 1) / 2).simplify()
    );
    return {(interior || wall).simplify(), std::move(exponents)};
}

std::vector<Source> build_sources(
    z3::context& context,
    const z3::expr& q,
    const z3::expr& r,
    const std::array<z3::expr, 3U>& even
) {
    std::vector<Source> sources;
    sources.reserve(12U);
    for (std::size_t selected_terminal = 0U;
         selected_terminal < 2U; ++selected_terminal) {
        const z3::expr& terminal = selected_terminal == 0U ? q : r;
        const z3::expr& complement = selected_terminal == 0U ? r : q;
        for (std::size_t selected = 0U; selected < even.size(); ++selected) {
            std::array<std::size_t, 2U> remaining{};
            std::size_t position = 0U;
            for (std::size_t vertex = 0U; vertex < even.size(); ++vertex) {
                if (vertex != selected) {
                    remaining[position] = vertex;
                    ++position;
                }
            }
            const Factor selected_part = selected_factor(
                context, terminal, even[selected], selected + 1U
            );
            for (int channel = 0; channel < 2; ++channel) {
                const Factor complement_part = complement_factor(
                    context, complement,
                    even[remaining[0U]], even[remaining[1U]],
                    remaining[0U] + 1U, remaining[1U] + 1U,
                    channel == 1
                );
                Source source{
                    (selected_part.active
                        && complement_part.active).simplify(),
                    zero_exponents(context), zero_exponents(context)
                };
                if (selected_terminal == 0U) {
                    source.q_exponents = selected_part.exponents;
                    source.r_exponents = complement_part.exponents;
                } else {
                    source.q_exponents = complement_part.exponents;
                    source.r_exponents = selected_part.exponents;
                }
                sources.push_back(std::move(source));
            }
        }
    }
    return sources;
}

z3::expr range_at_most_one(
    z3::context& context,
    const std::vector<Source>& sources,
    const std::vector<z3::expr>& members,
    bool remaining_is_q,
    std::size_t first_coordinate,
    std::size_t second_coordinate = 4U
) {
    z3::expr answer = context.bool_val(true);
    for (std::size_t first = 0U; first < sources.size(); ++first) {
        const Exponents& first_exponents = remaining_is_q
            ? sources[first].q_exponents
            : sources[first].r_exponents;
        const z3::expr first_value = second_coordinate == 4U
            ? first_exponents[first_coordinate]
            : first_exponents[first_coordinate]
                + first_exponents[second_coordinate];
        for (std::size_t second = first + 1U;
             second < sources.size(); ++second) {
            const Exponents& second_exponents = remaining_is_q
                ? sources[second].q_exponents
                : sources[second].r_exponents;
            const z3::expr second_value = second_coordinate == 4U
                ? second_exponents[first_coordinate]
                : second_exponents[first_coordinate]
                    + second_exponents[second_coordinate];
            answer = answer && z3::implies(
                members[first] && members[second],
                first_value - second_value <= 1
                    && second_value - first_value <= 1
            );
        }
    }
    return answer;
}

z3::expr v2_shape_certificate(
    z3::context& context,
    const std::vector<Source>& sources,
    const std::vector<z3::expr>& members,
    bool remaining_is_q
) {
    std::array<z3::expr, 4U> coordinate_binary{
        range_at_most_one(context, sources, members, remaining_is_q, 0U),
        range_at_most_one(context, sources, members, remaining_is_q, 1U),
        range_at_most_one(context, sources, members, remaining_is_q, 2U),
        range_at_most_one(context, sources, members, remaining_is_q, 3U)
    };
    z3::expr direct = context.bool_val(false);
    for (std::size_t first = 0U; first < 4U; ++first) {
        for (std::size_t second = first + 1U; second < 4U; ++second) {
            direct = direct
                || (coordinate_binary[first] && coordinate_binary[second]);
        }
    }
    z3::expr after_merger = context.bool_val(false);
    for (std::size_t merged_first = 0U; merged_first < 4U; ++merged_first) {
        for (std::size_t merged_second = merged_first + 1U;
             merged_second < 4U; ++merged_second) {
            z3::expr one_binary = range_at_most_one(
                context, sources, members, remaining_is_q,
                merged_first, merged_second
            );
            for (std::size_t coordinate = 0U; coordinate < 4U;
                 ++coordinate) {
                if (coordinate != merged_first
                    && coordinate != merged_second) {
                    one_binary = one_binary || coordinate_binary[coordinate];
                }
            }
            after_merger = after_merger || one_binary;
        }
    }
    return direct || after_merger;
}

z3::expr filtration_has_v2_shape(
    z3::context& context,
    const std::vector<Source>& sources,
    bool jet_is_q,
    std::size_t jet_root
) {
    z3::expr success = context.bool_val(true);
    for (std::size_t seed = 0U; seed < sources.size(); ++seed) {
        const Exponents& seed_jet = jet_is_q
            ? sources[seed].q_exponents
            : sources[seed].r_exponents;
        std::vector<z3::expr> members;
        members.reserve(sources.size());
        for (const Source& source : sources) {
            const Exponents& source_jet = jet_is_q
                ? source.q_exponents : source.r_exponents;
            members.push_back(
                source.active
                && source_jet[jet_root] == seed_jet[jet_root]
            );
        }
        success = success && z3::implies(
            sources[seed].active,
            v2_shape_certificate(
                context, sources, members, !jet_is_q
            )
        );
    }
    return success;
}

struct RootTransform {
    bool merged;
    std::size_t first;
    std::size_t second;
};

std::vector<RootTransform> root_transforms() {
    std::vector<RootTransform> transforms{{false, 0U, 0U}};
    for (std::size_t first = 0U; first < 4U; ++first) {
        for (std::size_t second = first + 1U; second < 4U; ++second) {
            transforms.push_back({true, first, second});
        }
    }
    return transforms;
}

std::vector<z3::expr> transformed_exponents(
    const Source& source,
    bool remaining_is_q,
    const RootTransform& transform
) {
    const Exponents& original = remaining_is_q
        ? source.q_exponents : source.r_exponents;
    if (!transform.merged) {
        return {
            original[0U], original[1U], original[2U], original[3U]
        };
    }
    std::vector<z3::expr> answer;
    answer.reserve(3U);
    for (std::size_t coordinate = 0U; coordinate < 4U; ++coordinate) {
        if (coordinate == transform.second) {
            continue;
        }
        answer.push_back(
            coordinate == transform.first
                ? original[transform.first] + original[transform.second]
                : original[coordinate]
        );
    }
    return answer;
}

z3::expr transformed_coordinate_binary(
    z3::context& context,
    const std::vector<Source>& sources,
    const std::vector<z3::expr>& members,
    bool remaining_is_q,
    const RootTransform& transform,
    std::size_t coordinate
) {
    z3::expr answer = context.bool_val(true);
    for (std::size_t first = 0U; first < sources.size(); ++first) {
        const std::vector<z3::expr> first_exponents
            = transformed_exponents(
                sources[first], remaining_is_q, transform
            );
        for (std::size_t second = first + 1U;
             second < sources.size(); ++second) {
            const std::vector<z3::expr> second_exponents
                = transformed_exponents(
                    sources[second], remaining_is_q, transform
                );
            answer = answer && z3::implies(
                members[first] && members[second],
                first_exponents[coordinate]
                        - second_exponents[coordinate] <= 1
                    && second_exponents[coordinate]
                        - first_exponents[coordinate] <= 1
            );
        }
    }
    return answer;
}

z3::expr transformed_v2_shape(
    z3::context& context,
    const std::vector<Source>& sources,
    const std::vector<z3::expr>& members,
    bool remaining_is_q,
    const RootTransform& transform
) {
    const std::size_t coordinates = transform.merged ? 3U : 4U;
    std::vector<z3::expr> binary;
    binary.reserve(coordinates);
    for (std::size_t coordinate = 0U; coordinate < coordinates;
         ++coordinate) {
        binary.push_back(transformed_coordinate_binary(
            context, sources, members, remaining_is_q,
            transform, coordinate
        ));
    }
    if (transform.merged) {
        z3::expr answer = context.bool_val(false);
        for (const z3::expr& coordinate_binary : binary) {
            answer = answer || coordinate_binary;
        }
        return answer;
    }
    z3::expr answer = context.bool_val(false);
    for (std::size_t first = 0U; first < binary.size(); ++first) {
        for (std::size_t second = first + 1U; second < binary.size();
             ++second) {
            answer = answer || (binary[first] && binary[second]);
        }
    }
    return answer;
}

z3::expr hall_failure_witness(
    z3::context& context,
    const std::vector<Source>& sources,
    const std::vector<z3::expr>& members,
    bool remaining_is_q,
    const RootTransform& transform,
    const z3::expr& degree,
    const std::string& prefix
) {
    const std::size_t coordinates = transform.merged ? 3U : 4U;

    // A chamber contains at most six potentially active source rows.  In
    // that regime enumerate the nonempty subsets literally.  This removes
    // all auxiliary Boolean selectors and minimum variables from the SMT
    // query and makes the encoded Hall failure visibly identical to
    //
    //   |I| + |coordinatewise_min_{i in I} w_i| > degree + 1.
    //
    // Keep the selector encoding below for the unfiltered twelve-source
    // modes, where literal enumeration would create 4095 disjuncts.
    if (sources.size() <= 8U) {
        z3::expr failure = context.bool_val(false);
        const std::uint64_t subset_count
            = std::uint64_t{1} << sources.size();
        for (std::uint64_t mask = 1U; mask < subset_count; ++mask) {
            z3::expr selected_are_members = context.bool_val(true);
            z3::expr gcd_degree = context.int_val(0);
            bool first_selected = true;
            std::vector<z3::expr> minima;
            minima.reserve(coordinates);
            for (std::size_t coordinate = 0U; coordinate < coordinates;
                 ++coordinate) {
                minima.push_back(context.int_val(0));
            }
            for (std::size_t row = 0U; row < sources.size(); ++row) {
                if ((mask & (std::uint64_t{1} << row)) == 0U) {
                    continue;
                }
                selected_are_members = selected_are_members && members[row];
                const std::vector<z3::expr> exponents
                    = transformed_exponents(
                        sources[row], remaining_is_q, transform
                    );
                for (std::size_t coordinate = 0U;
                     coordinate < coordinates; ++coordinate) {
                    if (first_selected) {
                        minima[coordinate] = exponents[coordinate];
                    } else {
                        minima[coordinate] = z3::ite(
                            minima[coordinate] <= exponents[coordinate],
                            minima[coordinate], exponents[coordinate]
                        ).simplify();
                    }
                }
                first_selected = false;
            }
            for (const z3::expr& minimum : minima) {
                gcd_degree = gcd_degree + minimum;
            }
            const int cardinality = std::popcount(mask);
            failure = failure
                || (selected_are_members
                    && context.int_val(cardinality) + gcd_degree
                        > degree + 1);
        }
        return failure;
    }

    std::vector<z3::expr> selected;
    selected.reserve(sources.size());
    z3::expr nonempty = context.bool_val(false);
    z3::expr cardinality = context.int_val(0);
    z3::expr constraints = context.bool_val(true);
    for (std::size_t row = 0U; row < sources.size(); ++row) {
        const z3::expr indicator = context.bool_const(
            (prefix + "_row_" + std::to_string(row)).c_str()
        );
        selected.push_back(indicator);
        nonempty = nonempty || indicator;
        cardinality = cardinality + z3::ite(
            indicator, context.int_val(1), context.int_val(0)
        );
        constraints = constraints && z3::implies(
            indicator, members[row]
        );
    }
    z3::expr gcd_degree = context.int_val(0);
    for (std::size_t coordinate = 0U; coordinate < coordinates;
         ++coordinate) {
        const z3::expr minimum = context.int_const(
            (prefix + "_min_" + std::to_string(coordinate)).c_str()
        );
        z3::expr attained = context.bool_val(false);
        for (std::size_t row = 0U; row < sources.size(); ++row) {
            const std::vector<z3::expr> exponents
                = transformed_exponents(
                    sources[row], remaining_is_q, transform
                );
            constraints = constraints && z3::implies(
                selected[row], minimum <= exponents[coordinate]
            );
            attained = attained
                || (selected[row] && minimum == exponents[coordinate]);
        }
        constraints = constraints && attained;
        gcd_degree = gcd_degree + minimum;
    }
    return nonempty && constraints
        && cardinality + gcd_degree > degree + 1;
}

z3::expr filtration_fails_full_v2(
    z3::context& context,
    const std::vector<Source>& sources,
    bool jet_is_q,
    std::size_t jet_root,
    const z3::expr& q,
    const z3::expr& r,
    const std::string& prefix
) {
    const std::vector<RootTransform> transforms = root_transforms();
    z3::expr failure = context.bool_val(false);
    for (std::size_t seed = 0U; seed < sources.size(); ++seed) {
        const Exponents& seed_jet = jet_is_q
            ? sources[seed].q_exponents
            : sources[seed].r_exponents;
        std::vector<z3::expr> members;
        members.reserve(sources.size());
        for (const Source& source : sources) {
            const Exponents& source_jet = jet_is_q
                ? source.q_exponents : source.r_exponents;
            members.push_back(
                source.active
                && source_jet[jet_root] == seed_jet[jet_root]
            );
        }
        z3::expr group_bad = context.bool_val(true);
        for (std::size_t transform_index = 0U;
             transform_index < transforms.size(); ++transform_index) {
            const RootTransform& transform = transforms[transform_index];
            const z3::expr shape = transformed_v2_shape(
                context, sources, members, !jet_is_q, transform
            );
            const z3::expr hall_failure = hall_failure_witness(
                context, sources, members, !jet_is_q, transform,
                jet_is_q ? r : q,
                prefix + "_seed_" + std::to_string(seed)
                    + "_transform_" + std::to_string(transform_index)
            );
            group_bad = group_bad && (!shape || hall_failure);
        }
        failure = failure || (sources[seed].active && group_bad);
    }
    return failure;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = argc >= 2 ? argv[1] : "--shape";
        const bool single_q3_large_chamber
            = mode == "--full-q3-large-chamber";
        const bool equal_strict_min_classification
            = mode == "--full-equal-strict-min-classification";
        const bool equal_classification
            = equal_strict_min_classification;
        const bool full_query = mode == "--full"
            || mode == "--full-q3"
            || mode == "--full-q3-r5"
            || mode == "--full-q3-r5-chambers"
            || mode == "--full-q3-large"
            || mode == "--full-q3-large-chambers"
            || single_q3_large_chamber
            || mode == "--full-adjacent"
            || mode == "--full-adjacent-chambers"
            || mode == "--full-separated"
            || mode == "--full-separated-chambers"
            || equal_classification;
        const bool valid_arguments
            = (single_q3_large_chamber && argc == 5)
            || (!single_q3_large_chamber && argc <= 2 && full_query)
            || (argc == 1);
        if (!valid_arguments) {
            std::cerr
                << "usage: verify_su2_v2_presburger"
                << " [--full|--full-q3|--full-q3-r5"
                << "|--full-q3-r5-chambers|--full-q3-large"
                << "|--full-q3-large-chambers"
                << "|--full-q3-large-chamber A B C"
                << "|--full-adjacent|--full-adjacent-chambers"
                << "|--full-separated|--full-separated-chambers"
                << "|--full-equal-strict-min-classification]\n";
            return EXIT_FAILURE;
        }
        std::array<int, 3U> requested_category{0, 0, 0};
        if (single_q3_large_chamber) {
            for (std::size_t index = 0U; index < requested_category.size();
                 ++index) {
                requested_category[index] = std::stoi(argv[index + 2U]);
            }
            if (requested_category[0U] < 1
                || requested_category[0U] > requested_category[1U]
                || requested_category[1U] > requested_category[2U]
                || requested_category[2U] > 6) {
                std::cerr << "invalid ordered q3-large chamber\n";
                return EXIT_FAILURE;
            }
        }
        std::cout << std::unitbuf;
        z3::context context;
        const z3::expr q_half = context.int_const("q_half");
        const z3::expr r_half = context.int_const("r_half");
        const bool q_is_three = mode == "--full-q3"
            || mode == "--full-q3-r5"
            || mode == "--full-q3-r5-chambers"
            || mode == "--full-q3-large"
            || mode == "--full-q3-large-chambers"
            || single_q3_large_chamber;
        const z3::expr q = q_is_three
            ? context.int_val(3) : 2 * q_half + 1;
        const bool r_is_five = mode == "--full-q3-r5"
            || mode == "--full-q3-r5-chambers";
        const z3::expr r = equal_classification
            ? q : (r_is_five ? context.int_val(5) : 2 * r_half + 1);
        const std::array<z3::expr, 3U> even_half{
            context.int_const("e1_half"),
            context.int_const("e2_half"),
            context.int_const("e3_half")
        };
        const std::array<z3::expr, 3U> even{
            2 * even_half[0U],
            2 * even_half[1U],
            2 * even_half[2U]
        };
        const std::vector<Source> sources
            = build_sources(context, q, r, even);

        z3::tactic qflia(context, "qflia");
        z3::solver solver = qflia.mk_solver();
        if (equal_classification) {
            solver.add(q_half >= 1);
        } else if (r_is_five) {
            // Both odd labels are already fixed syntactically.
        } else if (q_is_three) {
            solver.add(r_half >= 2);
        } else {
            solver.add(q_half >= 1 && q_half < r_half);
        }
        solver.add(even_half[0U] >= 1
            && even_half[0U] <= even_half[1U]
            && even_half[1U] <= even_half[2U]);
        z3::expr nonempty = context.bool_val(false);
        for (const Source& source : sources) {
            nonempty = nonempty || source.active;
        }
        solver.add(nonempty);

        const z3::expr circuit_family
            = even[0U] == 2 && even[1U] == q + 1
            && even[2U] == r + 1;
        const z3::expr adjacent_neighbors
            = z3::ite(
                  even[0U] == q + 1,
                  context.int_val(1), context.int_val(0)
              )
            + z3::ite(
                  even[1U] == q + 1,
                  context.int_val(1), context.int_val(0)
              )
            + z3::ite(
                  even[2U] == q + 1,
                  context.int_val(1), context.int_val(0)
              );
        const z3::expr adjacent_two_neighbor_family
            = r == q + 2 && adjacent_neighbors >= 2;
        const z3::expr cubic_repeated_neighbor_family
            = q == 3 && even[0U] == even[1U]
            && even[1U] == even[2U]
            && (even[0U] == r - 1 || even[0U] == r + 1);
        const z3::expr cubic_special_family
            = q == 3 && r == 5 && even[0U] == 4
            && even[1U] == 4 && even[2U] == 8;
        if (!equal_classification) {
            solver.add(!(circuit_family || adjacent_two_neighbor_family
                || cubic_repeated_neighbor_family || cubic_special_family));
        }

        const auto all_filtrations_fail_for
            = [&](const std::vector<Source>& candidate_sources,
                  const std::string& prefix) {
                z3::expr all_fail = context.bool_val(true);
                if (q_is_three) {
                    for (std::size_t root = 0U; root < 3U; ++root) {
                        all_fail = all_fail
                            && filtration_fails_full_v2(
                                   context, candidate_sources, true,
                                   root, q, r,
                                   prefix + "_q_" + std::to_string(root)
                               );
                    }
                    all_fail = all_fail
                        && filtration_fails_full_v2(
                               context, candidate_sources, false,
                               0U, q, r, prefix + "_r_0"
                           );
                } else {
                    for (std::size_t root = 0U; root < 4U; ++root) {
                        all_fail = all_fail
                            && filtration_fails_full_v2(
                                   context, candidate_sources, true,
                                   root, q, r,
                                   prefix + "_q_" + std::to_string(root)
                               )
                            && filtration_fails_full_v2(
                                   context, candidate_sources, false,
                                   root, q, r,
                                   prefix + "_r_" + std::to_string(root)
                               );
                    }
                }
                return all_fail;
            };
        const auto uniform_filtration_choice
            = [&](const std::vector<Source>& candidate_sources,
                  const z3::expr_vector& assumptions,
                  const std::string& prefix) {
                for (std::size_t choice = 0U; choice < 8U; ++choice) {
                    const bool jet_is_q = choice < 4U;
                    const std::size_t root = choice % 4U;
                    solver.push();
                    solver.add(assumptions);
                    solver.add(filtration_fails_full_v2(
                        context, candidate_sources, jet_is_q, root,
                        q, r, prefix + "_choice_" + std::to_string(choice)
                    ));
                    const z3::check_result result = solver.check();
                    solver.pop();
                    if (result == z3::unsat) {
                        return static_cast<int>(choice);
                    }
                }
                return -1;
            };

        const z3::expr equal_low_three_residual = q == 3 && (
            (even[0U] == 2 && even[1U] == 2
                && (even[2U] == 2 || even[2U] == 4))
            || (even[0U] == 2 && even[1U] == 4
                && (even[2U] == 4 || even[2U] == 6))
            || (even[0U] == 4 && even[1U] == 4
                && (even[2U] == 4 || even[2U] == 6
                    || even[2U] == 8))
        );
        const z3::expr equal_low_five_residual = q == 5 && (
            (even[0U] == 2 && even[1U] == 4 && even[2U] == 6)
            || (even[0U] == 2 && even[1U] == 6 && even[2U] == 6)
            || (even[0U] == 4 && even[1U] == 4
                && (even[2U] == 4 || even[2U] == 6
                    || even[2U] == 8))
            || (even[0U] == 4 && even[1U] == 6
                && (even[2U] == 6 || even[2U] == 8
                    || even[2U] == 10))
            || (even[0U] == 6 && even[1U] == 6
                && (even[2U] == 6 || even[2U] == 8
                    || even[2U] == 10 || even[2U] == 12))
        );
        const z3::expr equal_high_residual = q >= 7 && (
            (even[0U] == 2 && even[1U] == q - 1
                && even[2U] == q + 1)
            || (even[0U] == 2 && even[1U] == q + 1
                && even[2U] == q + 1)
            || (even[0U] == q - 1 && even[1U] == q - 1
                && even[2U] == q - 1)
            || (even[0U] == q - 1 && even[1U] == q - 1
                && even[2U] == q + 1)
            || (even[0U] == q - 1 && even[1U] == q + 1
                && even[2U] == q + 1)
            || (even[0U] == q + 1 && even[1U] == q + 1
                && even[2U] == q + 1)
            || (even[0U] == q - 1 && even[1U] == q - 1
                && even[2U] == 2 * q - 4)
            || (even[0U] == q - 1 && even[1U] == q - 1
                && even[2U] == 2 * q - 2)
            || (even[0U] == q - 1 && even[1U] == q + 1
                && even[2U] == 2 * q - 2)
            || (even[0U] == q - 1 && even[1U] == q + 1
                && even[2U] == 2 * q)
            || (even[0U] == q + 1 && even[1U] == q + 1
                && even[2U] == 2 * q - 2)
            || (even[0U] == q + 1 && even[1U] == q + 1
                && even[2U] == 2 * q)
            || (even[0U] == q + 1 && even[1U] == q + 1
                && even[2U] == 2 * q + 2)
        );
        const z3::expr equal_predicted_residual
            = equal_low_three_residual || equal_low_five_residual
                || equal_high_residual;

        if (equal_strict_min_classification) {
            // At q=r, sources 0,...,5 pair under terminal exchange with
            // sources 6,...,11.  Filter an orbit sum by the minimum order at
            // one root; its leading terminal form is the factor with the
            // larger exponent.  This is exactly the strict-minimum V2 test
            // used by analyze_equal_terminal_factors.
            z3::expr orbit_nonempty = context.bool_val(false);
            z3::expr all_roots_fail = context.bool_val(true);
            std::vector<z3::expr> root_ties;
            std::vector<z3::expr> root_filtration_failures;
            for (std::size_t root = 0U; root < 4U; ++root) {
                z3::expr tied_pair = context.bool_val(false);
                std::vector<Source> root_sources;
                root_sources.reserve(6U);
                for (std::size_t source = 0U; source < 6U; ++source) {
                    const Source& first = sources[source];
                    const Source& second = sources[source + 6U];
                    const z3::expr active
                        = (first.active && second.active).simplify();
                    const z3::expr first_jet
                        = first.q_exponents[root];
                    const z3::expr second_jet
                        = first.r_exponents[root];
                    Source orbit{
                        active, zero_exponents(context),
                        zero_exponents(context)
                    };
                    orbit.q_exponents[root] = z3::ite(
                        first_jet <= second_jet, first_jet, second_jet
                    ).simplify();
                    for (std::size_t coordinate = 0U; coordinate < 4U;
                         ++coordinate) {
                        // Choose one whole leading factor.  A coordinatewise
                        // maximum would splice two different binary forms.
                        orbit.r_exponents[coordinate] = z3::ite(
                            first_jet < second_jet,
                            first.r_exponents[coordinate],
                            first.q_exponents[coordinate]
                        ).simplify();
                    }
                    root_sources.push_back(std::move(orbit));
                    tied_pair = tied_pair
                        || (active && first_jet == second_jet);
                    orbit_nonempty = orbit_nonempty
                        || active;
                }
                const z3::expr filtration_failure
                    = filtration_fails_full_v2(
                          context, root_sources, true, root, q, q,
                          "equal_strict_min_" + std::to_string(root)
                      );
                root_ties.push_back(tied_pair);
                root_filtration_failures.push_back(filtration_failure);
                all_roots_fail = all_roots_fail
                    && (tied_pair || filtration_failure);
            }
            const z3::expr predicted_residual
                = equal_predicted_residual;
            solver.add(orbit_nonempty && all_roots_fail
                && !predicted_residual);
            const z3::check_result result = solver.check();
            std::cout
                << "SU2_EQUAL_STRICT_MIN_CLASSIFICATION q>=3 result="
                << result << " low_cases=19 residual_families=13\n";
            if (result == z3::sat) {
                const z3::model model = solver.get_model();
                std::cout << model << '\n';
                for (std::size_t root = 0U; root < 4U; ++root) {
                    std::cout << "root=" << root << " tied="
                              << model.eval(root_ties[root], true)
                              << " filtration_failure="
                              << model.eval(
                                     root_filtration_failures[root], true
                                 ) << '\n';
                }
            }
            return result == z3::unsat ? EXIT_SUCCESS : EXIT_FAILURE;
        }

        if (full_query) {
            if (mode == "--full-q3") {
                solver.add(q_half == 1);
            } else if (r_is_five) {
                solver.add(q_half == 1 && r_half == 2);
            } else if (mode == "--full-q3-large"
                || mode == "--full-q3-large-chambers"
                || single_q3_large_chamber) {
                solver.add(q_half == 1 && r_half >= 3);
            } else if (mode == "--full-adjacent"
                || mode == "--full-adjacent-chambers") {
                solver.add(q_half >= 2 && r_half == q_half + 1);
            } else if (mode == "--full-separated"
                || mode == "--full-separated-chambers") {
                solver.add(q_half >= 2 && r_half >= q_half + 2);
            }
            const bool chamber_query
                = mode == "--full-q3-r5-chambers"
                || mode == "--full-q3-large-chambers"
                || single_q3_large_chamber
                || mode == "--full-adjacent-chambers"
                || mode == "--full-separated-chambers";
            if (!chamber_query) {
                solver.add(all_filtrations_fail_for(sources, "global"));
            }
        } else {
            z3::expr some_filtration = context.bool_val(false);
            for (std::size_t root = 0U; root < 4U; ++root) {
                some_filtration = some_filtration
                    || filtration_has_v2_shape(
                           context, sources, true, root
                       )
                    || filtration_has_v2_shape(
                           context, sources, false, root
                       );
            }
            solver.add(!some_filtration);
        }

        if (mode == "--full-q3-r5-chambers") {
            std::size_t chamber_count = 0U;
            for (int first = 1; first <= 4; ++first) {
                for (int second = first; second <= 4; ++second) {
                    for (int third = second; third <= 4; ++third) {
                        const std::array<int, 3U> category{
                            first, second, third
                        };
                        z3::expr_vector assumptions(context);
                        for (std::size_t index = 0U;
                             index < even_half.size(); ++index) {
                            assumptions.push_back(
                                category[index] < 4
                                    ? even_half[index] == category[index]
                                    : even_half[index] >= 4
                            );
                        }
                        std::vector<Source> chamber_sources;
                        for (std::size_t index = 0U;
                             index < category.size(); ++index) {
                            if (category[index] == 1
                                || category[index] == 2) {
                                chamber_sources.push_back(
                                    sources[2U * index]
                                );
                                chamber_sources.push_back(
                                    sources[2U * index + 1U]
                                );
                            }
                            if (category[index] == 2
                                || category[index] == 3) {
                                chamber_sources.push_back(
                                    sources[6U + 2U * index]
                                );
                                chamber_sources.push_back(
                                    sources[6U + 2U * index + 1U]
                                );
                            }
                        }
                        ++chamber_count;
                        if (chamber_sources.empty()) {
                            std::cout
                                << "SU2_V2_PRESBURGER_CHAMBER categories="
                                << first << ',' << second << ',' << third
                                << " result=empty\n";
                            continue;
                        }
                        int uniform_choice = -1;
                        constexpr std::array<int, 4U> q3_choices{
                            0, 1, 2, 4
                        };
                        for (int choice : q3_choices) {
                            const bool jet_is_q = choice < 4;
                            const std::size_t root
                                = static_cast<std::size_t>(choice % 4);
                            solver.push();
                            solver.add(assumptions);
                            solver.add(filtration_fails_full_v2(
                                context, chamber_sources, jet_is_q, root,
                                q, r, "r5_uniform_"
                                    + std::to_string(first)
                                    + std::to_string(second)
                                    + std::to_string(third)
                                    + "_choice_" + std::to_string(choice)
                            ));
                            const z3::check_result result = solver.check();
                            solver.pop();
                            if (result == z3::unsat) {
                                uniform_choice = choice;
                                break;
                            }
                        }
                        if (uniform_choice >= 0) {
                            std::cout
                                << "SU2_V2_PRESBURGER_CHAMBER categories="
                                << first << ',' << second << ',' << third
                                << " result=unsat uniform_choice="
                                << (uniform_choice < 4 ? 'q' : 'r')
                                << (uniform_choice % 4) << '\n';
                            continue;
                        }
                        solver.push();
                        solver.add(assumptions);
                        solver.add(all_filtrations_fail_for(
                            chamber_sources,
                            "r5_" + std::to_string(first)
                                + std::to_string(second)
                                + std::to_string(third)
                        ));
                        const z3::check_result chamber_result = solver.check();
                        std::cout << "SU2_V2_PRESBURGER_CHAMBER categories="
                                  << first << ',' << second << ',' << third
                                  << " result=" << chamber_result << '\n';
                        if (chamber_result != z3::unsat) {
                            if (chamber_result == z3::sat) {
                                std::cout << solver.get_model() << '\n';
                            }
                            solver.pop();
                            return EXIT_FAILURE;
                        }
                        solver.pop();
                    }
                }
            }
            std::cout << "SU2_V2_PRESBURGER_FULL mode=" << mode
                      << " chambers=" << chamber_count
                      << " result=unsat\n";
            return EXIT_SUCCESS;
        }
        if (mode == "--full-q3-large-chambers"
            || single_q3_large_chamber) {
            std::size_t chamber_count = 0U;
            for (int first = 1; first <= 6; ++first) {
                for (int second = first; second <= 6; ++second) {
                    for (int third = second; third <= 6; ++third) {
                        const std::array<int, 3U> category{
                            first, second, third
                        };
                        if (single_q3_large_chamber
                            && category != requested_category) {
                            continue;
                        }
                        z3::expr_vector assumptions(context);
                        for (std::size_t index = 0U;
                             index < even_half.size(); ++index) {
                            z3::expr condition = context.bool_val(false);
                            switch (category[index]) {
                                case 1:
                                    condition = even_half[index] == 1;
                                    break;
                                case 2:
                                    condition = even_half[index] == 2;
                                    break;
                                case 3:
                                    condition = even_half[index] >= 3
                                        && even_half[index] <= r_half - 1;
                                    break;
                                case 4:
                                    condition = even_half[index] == r_half;
                                    break;
                                case 5:
                                    condition = even_half[index]
                                        == r_half + 1;
                                    break;
                                case 6:
                                    condition = even_half[index]
                                        >= r_half + 2;
                                    break;
                                default:
                                    throw std::runtime_error(
                                        "invalid q3-large chamber category"
                                    );
                            }
                            assumptions.push_back(condition);
                        }
                        std::vector<Source> chamber_sources;
                        for (std::size_t index = 0U;
                             index < category.size(); ++index) {
                            if (category[index] == 1
                                || category[index] == 2) {
                                chamber_sources.push_back(
                                    sources[2U * index]
                                );
                                chamber_sources.push_back(
                                    sources[2U * index + 1U]
                                );
                            }
                            if (category[index] == 4
                                || category[index] == 5) {
                                chamber_sources.push_back(
                                    sources[6U + 2U * index]
                                );
                                chamber_sources.push_back(
                                    sources[6U + 2U * index + 1U]
                                );
                            }
                        }
                        solver.push();
                        solver.add(assumptions);
                        solver.add(all_filtrations_fail_for(
                            chamber_sources,
                            "q3_" + std::to_string(first)
                                + std::to_string(second)
                                + std::to_string(third)
                        ));
                        const z3::check_result chamber_result = solver.check();
                        ++chamber_count;
                        std::cout << "SU2_V2_PRESBURGER_CHAMBER categories="
                                  << first << ',' << second << ',' << third
                                  << " result=" << chamber_result << '\n';
                        if (chamber_result != z3::unsat) {
                            if (chamber_result == z3::sat) {
                                std::cout << solver.get_model() << '\n';
                            }
                            solver.pop();
                            return EXIT_FAILURE;
                        }
                        solver.pop();
                    }
                }
            }
            std::cout << "SU2_V2_PRESBURGER_FULL mode=" << mode
                      << " chambers=" << chamber_count
                      << " result=unsat\n";
            return EXIT_SUCCESS;
        }
        if (mode == "--full-adjacent-chambers") {
            std::size_t chamber_count = 0U;
            for (int first = 1; first <= 5; ++first) {
                for (int second = first; second <= 5; ++second) {
                    for (int third = second; third <= 5; ++third) {
                        const std::array<int, 3U> category{
                            first, second, third
                        };
                        z3::expr_vector assumptions(context);
                        for (std::size_t index = 0U;
                             index < even_half.size(); ++index) {
                            z3::expr condition = context.bool_val(false);
                            switch (category[index]) {
                                case 1:
                                    condition = even_half[index]
                                        <= q_half - 1;
                                    break;
                                case 2:
                                    condition = even_half[index] == q_half;
                                    break;
                                case 3:
                                    condition = even_half[index]
                                        == q_half + 1;
                                    break;
                                case 4:
                                    condition = even_half[index]
                                        == q_half + 2;
                                    break;
                                case 5:
                                    condition = even_half[index]
                                        >= q_half + 3;
                                    break;
                                default:
                                    throw std::runtime_error(
                                        "invalid adjacent chamber category"
                                    );
                            }
                            assumptions.push_back(condition);
                        }
                        std::vector<Source> chamber_sources;
                        std::size_t shared_neighbor_count = 0U;
                        for (std::size_t index = 0U;
                             index < category.size(); ++index) {
                            if (category[index] == 3) {
                                ++shared_neighbor_count;
                            }
                            if (category[index] == 2
                                || category[index] == 3) {
                                chamber_sources.push_back(
                                    sources[2U * index]
                                );
                                chamber_sources.push_back(
                                    sources[2U * index + 1U]
                                );
                            }
                            if (category[index] == 3
                                || category[index] == 4) {
                                chamber_sources.push_back(
                                    sources[6U + 2U * index]
                                );
                                chamber_sources.push_back(
                                    sources[6U + 2U * index + 1U]
                                );
                            }
                        }
                        ++chamber_count;
                        if (shared_neighbor_count >= 2U) {
                            std::cout
                                << "SU2_V2_PRESBURGER_CHAMBER categories="
                                << first << ',' << second << ',' << third
                                << " result=excluded-adjacent-two-neighbor\n";
                            continue;
                        }
                        if (chamber_sources.empty()) {
                            std::cout
                                << "SU2_V2_PRESBURGER_CHAMBER categories="
                                << first << ',' << second << ',' << third
                                << " result=empty\n";
                            continue;
                        }
                        const int uniform_choice = uniform_filtration_choice(
                            chamber_sources, assumptions,
                            "adj_uniform_" + std::to_string(first)
                                + std::to_string(second)
                                + std::to_string(third)
                        );
                        if (uniform_choice >= 0) {
                            std::cout
                                << "SU2_V2_PRESBURGER_CHAMBER categories="
                                << first << ',' << second << ',' << third
                                << " result=unsat uniform_choice="
                                << (uniform_choice < 4 ? 'q' : 'r')
                                << (uniform_choice % 4) << '\n';
                            continue;
                        }
                        solver.push();
                        solver.add(assumptions);
                        solver.add(all_filtrations_fail_for(
                            chamber_sources,
                            "adj_" + std::to_string(first)
                                + std::to_string(second)
                                + std::to_string(third)
                        ));
                        const z3::check_result chamber_result = solver.check();
                        std::cout << "SU2_V2_PRESBURGER_CHAMBER categories="
                                  << first << ',' << second << ',' << third
                                  << " result=" << chamber_result << '\n';
                        if (chamber_result != z3::unsat) {
                            if (chamber_result == z3::sat) {
                                std::cout << solver.get_model() << '\n';
                            }
                            solver.pop();
                            return EXIT_FAILURE;
                        }
                        solver.pop();
                    }
                }
            }
            std::cout << "SU2_V2_PRESBURGER_FULL mode=" << mode
                      << " chambers=" << chamber_count
                      << " result=unsat\n";
            return EXIT_SUCCESS;
        }
        if (mode == "--full-separated-chambers") {
            std::size_t chamber_count = 0U;
            for (int first = 1; first <= 7; ++first) {
                for (int second = first; second <= 7; ++second) {
                    for (int third = second; third <= 7; ++third) {
                        const std::array<int, 3U> category{
                            first, second, third
                        };
                        z3::expr_vector assumptions(context);
                        for (std::size_t index = 0U;
                             index < even_half.size(); ++index) {
                            z3::expr condition = context.bool_val(false);
                            switch (category[index]) {
                                case 1:
                                    condition = even_half[index]
                                        <= q_half - 1;
                                    break;
                                case 2:
                                    condition = even_half[index] == q_half;
                                    break;
                                case 3:
                                    condition = even_half[index]
                                        == q_half + 1;
                                    break;
                                case 4:
                                    condition = even_half[index]
                                            >= q_half + 2
                                        && even_half[index] <= r_half - 1;
                                    break;
                                case 5:
                                    condition = even_half[index] == r_half;
                                    break;
                                case 6:
                                    condition = even_half[index]
                                        == r_half + 1;
                                    break;
                                case 7:
                                    condition = even_half[index]
                                        >= r_half + 2;
                                    break;
                                default:
                                    throw std::runtime_error(
                                        "invalid separated chamber category"
                                    );
                            }
                            assumptions.push_back(condition);
                        }
                        std::vector<Source> chamber_sources;
                        for (std::size_t index = 0U;
                             index < category.size(); ++index) {
                            if (category[index] == 2
                                || category[index] == 3) {
                                chamber_sources.push_back(
                                    sources[2U * index]
                                );
                                chamber_sources.push_back(
                                    sources[2U * index + 1U]
                                );
                            }
                            if (category[index] == 5
                                || category[index] == 6) {
                                chamber_sources.push_back(
                                    sources[6U + 2U * index]
                                );
                                chamber_sources.push_back(
                                    sources[6U + 2U * index + 1U]
                                );
                            }
                        }
                        ++chamber_count;
                        if (chamber_sources.empty()) {
                            std::cout
                                << "SU2_V2_PRESBURGER_CHAMBER categories="
                                << first << ',' << second << ',' << third
                                << " result=empty\n";
                            continue;
                        }
                        const int uniform_choice = uniform_filtration_choice(
                            chamber_sources, assumptions,
                            "sep_uniform_" + std::to_string(first)
                                + std::to_string(second)
                                + std::to_string(third)
                        );
                        if (uniform_choice >= 0) {
                            std::cout
                                << "SU2_V2_PRESBURGER_CHAMBER categories="
                                << first << ',' << second << ',' << third
                                << " result=unsat uniform_choice="
                                << (uniform_choice < 4 ? 'q' : 'r')
                                << (uniform_choice % 4) << '\n';
                            continue;
                        }
                        solver.push();
                        solver.add(assumptions);
                        solver.add(all_filtrations_fail_for(
                            chamber_sources,
                            "sep_" + std::to_string(first)
                                + std::to_string(second)
                                + std::to_string(third)
                        ));
                        const z3::check_result chamber_result = solver.check();
                        std::cout << "SU2_V2_PRESBURGER_CHAMBER categories="
                                  << first << ',' << second << ',' << third
                                  << " result=" << chamber_result << '\n';
                        if (chamber_result != z3::unsat) {
                            if (chamber_result == z3::sat) {
                                std::cout << solver.get_model() << '\n';
                            }
                            solver.pop();
                            return EXIT_FAILURE;
                        }
                        solver.pop();
                    }
                }
            }
            std::cout << "SU2_V2_PRESBURGER_FULL mode=" << mode
                      << " chambers=" << chamber_count
                      << " result=unsat\n";
            return EXIT_SUCCESS;
        }
        const z3::check_result result = solver.check();
        std::cout << (full_query
            ? "SU2_V2_PRESBURGER_FULL mode=" + mode + " result="
            : "SU2_V2_PRESBURGER_SHAPE result=")
                  << result << '\n';
        if (result == z3::sat) {
            std::cout << solver.get_model() << '\n';
            return EXIT_FAILURE;
        }
        return result == z3::unsat ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const z3::exception& error) {
        std::cerr << "z3 error: " << error.msg() << '\n';
        return EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
