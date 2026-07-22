#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace {

int multiplicity(const std::vector<int>& labels, int target) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    if (target < 0 || target > total) {
        return 0;
    }
    std::vector<int> current(static_cast<std::size_t>(total + 1));
    current[0] = 1;
    int support = 0;
    for (int label : labels) {
        std::vector<int> next(current.size());
        for (int input = 0; input <= support; ++input) {
            const int coefficient = current[static_cast<std::size_t>(input)];
            if (coefficient == 0) {
                continue;
            }
            for (int output = std::abs(input - label);
                 output <= input + label;
                 output += 2) {
                next[static_cast<std::size_t>(output)] += coefficient;
            }
        }
        support += label;
        current = std::move(next);
    }
    return current[static_cast<std::size_t>(target)];
}

struct Decomposition {
    std::array<int, 5> layers{};
    std::array<int, 5> positive_layers{};
    std::array<int, 5> negative_layers{};
    std::array<int, 5> first_negative_layers{};
    std::array<int, 5> second_negative_layers{};
    int positive = 0;
    int negative = 0;
};

Decomposition decomposition(
    int p,
    const std::array<int, 4>& suffix,
    int target
) {
    Decomposition result;
    for (std::uint32_t mask = 0; mask < 16U; ++mask) {
        std::vector<int> first;
        std::vector<int> second;
        for (int index = 0; index < 4; ++index) {
            const int label = suffix[static_cast<std::size_t>(index)];
            if (((mask >> index) & 1U) != 0U) {
                second.push_back(label);
            } else {
                first.push_back(label);
            }
        }
        std::vector<int> first_p_minus = first;
        first_p_minus.push_back(p - 1);
        std::vector<int> first_p = first;
        first_p.push_back(p);
        std::vector<int> first_one = first;
        first_one.push_back(1);
        const int positive = multiplicity(second, 0)
                * multiplicity(first_p_minus, target)
            + multiplicity(second, p)
                * multiplicity(first_one, target);
        const int first_negative = multiplicity(second, p - 1)
            * multiplicity(first, target);
        const int second_negative = multiplicity(second, 1)
            * multiplicity(first_p, target);
        const int negative = first_negative + second_negative;
        const int layer = std::popcount(mask);
        result.layers[static_cast<std::size_t>(layer)]
            += positive - negative;
        result.positive_layers[static_cast<std::size_t>(layer)] += positive;
        result.negative_layers[static_cast<std::size_t>(layer)] += negative;
        result.first_negative_layers[static_cast<std::size_t>(layer)]
            += first_negative;
        result.second_negative_layers[static_cast<std::size_t>(layer)]
            += second_negative;
        result.positive += positive;
        result.negative += negative;
    }
    return result;
}

int one_terminal_candidate_count(
    int p,
    int q,
    int r,
    int s,
    int t,
    int target,
    int width
) {
    const int total = s + t + target;
    const int largest = std::max({s, t, target});
    const int defect = 2 * largest - total;
    const int support_low = defect > 0 ? defect : total % 2;
    int result = 0;
    for (int first_index = 0; first_index < p; ++first_index) {
        const int first = q - p + 1 + 2 * first_index;
        if ((first + r - total) % 2 != 0) {
            continue;
        }
        const int low = std::max(std::abs(first - r), support_low);
        const int high = std::min(first + r, total);
        if (low <= high) {
            const int second_count
                = std::min(width, (high - low) / 2 + 1);
            for (int second_index = 0;
                 second_index < second_count;
                 ++second_index) {
                const int second = low + 2 * second_index;
                if ((second + s - target - t) % 2 != 0) {
                    continue;
                }
                const int third_low = std::max(
                    std::abs(second - s), std::abs(target - t)
                );
                const int third_high = std::min(
                    second + s, target + t
                );
                if (third_low <= third_high) {
                    const int terminal_width = p >= 3 ? 3 : 4;
                    result += std::min(
                        terminal_width,
                        (third_high - third_low) / 2 + 1
                    );
                }
            }
        }
    }
    return result;
}

void print_case(const std::tuple<int, int, int, int, int, int>& value) {
    const auto [p, q, r, s, t, target] = value;
    std::cout << p << ',' << q << ',' << r << ',' << s << ',' << t << ','
              << target;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2 && argc != 7) {
            throw std::runtime_error(
                "usage: analyze_su2_opd_four_suffix MAXIMUM_LABEL\n"
                "   or: analyze_su2_opd_four_suffix p q r s t target"
            );
        }
        if (argc == 7) {
            const int p = std::stoi(argv[1]);
            const int q = std::stoi(argv[2]);
            const int r = std::stoi(argv[3]);
            const int s = std::stoi(argv[4]);
            const int t = std::stoi(argv[5]);
            const int target = std::stoi(argv[6]);
            if (p < 2 || p > q || q > r || r > s || s > t
                || target < 1) {
                throw std::runtime_error("invalid ordered case");
            }
            const Decomposition current = decomposition(
                p, {q, r, s, t}, target
            );
            int prefix = 0;
            std::cout << "SU2_OPD_FOUR_SUFFIX_CASE p=" << p
                      << " q=" << q << " r=" << r << " s=" << s
                      << " t=" << t << " target=" << target
                      << " positive=" << current.positive
                      << " negative=" << current.negative << '\n';
            for (int layer = 0; layer <= 4; ++layer) {
                prefix += current.layers[static_cast<std::size_t>(layer)];
                std::cout << "layer=" << layer << " contribution="
                          << current.layers[static_cast<std::size_t>(layer)]
                          << " positive="
                          << current.positive_layers[
                                 static_cast<std::size_t>(layer)]
                          << " negative="
                          << current.negative_layers[
                                 static_cast<std::size_t>(layer)]
                          << " first_negative="
                          << current.first_negative_layers[
                                 static_cast<std::size_t>(layer)]
                          << " second_negative="
                          << current.second_negative_layers[
                                 static_cast<std::size_t>(layer)]
                          << " prefix=" << prefix << '\n';
            }
            return prefix < 0 ? EXIT_FAILURE : EXIT_SUCCESS;
        }
        const int maximum_label = std::stoi(argv[1]);
        if (maximum_label < 2 || maximum_label > 30) {
            throw std::runtime_error("invalid maximum label");
        }
        long long tested = 0;
        long long disjoint_tested = 0;
        int minimum_total = std::numeric_limits<int>::max();
        int minimum_main_slack = std::numeric_limits<int>::max();
        int minimum_disjoint_total = std::numeric_limits<int>::max();
        int minimum_disjoint_main_slack = std::numeric_limits<int>::max();
        int maximum_disjoint_negative = 0;
        int minimum_disjoint_main_singleton_slack
            = std::numeric_limits<int>::max();
        int minimum_disjoint_endpoint_slack
            = std::numeric_limits<int>::max();
        int minimum_disjoint_boundary_slack
            = std::numeric_limits<int>::max();
        int minimum_active_disjoint_boundary_slack
            = std::numeric_limits<int>::max();
        int minimum_amplification_slack = std::numeric_limits<int>::max();
        int minimum_disjoint_interior_main_singleton_slack
            = std::numeric_limits<int>::max();
        int minimum_disjoint_strict_interior_main_slack
            = std::numeric_limits<int>::max();
        int minimum_active_disjoint_interior_main_singleton_slack
            = std::numeric_limits<int>::max();
        int minimum_active_disjoint_strict_interior_main_slack
            = std::numeric_limits<int>::max();
        std::array<int, 3> minimum_shell{
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max()
        };
        std::array<int, 2> minimum_parity{
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max()
        };
        std::array<int, 5> minimum_prefix{
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max()
        };
        std::array<int, 31> minimum_active_main_by_p{};
        std::array<int, 31> maximum_active_negative_by_p{};
        std::array<int, 31> minimum_active_main_slack_by_p{};
        std::array<int, 31> minimum_active_boundary_slack_by_p{};
        std::array<long long, 31> residue_cases_by_p{};
        std::array<long long, 31> residue_equal_minimum_by_p{};
        std::array<long long, 31> residue_endpoint_by_p{};
        std::array<int, 31> minimum_residue_coarse_slack_by_p{};
        std::array<int, 31> minimum_residue_coarse_boundary_slack_by_p{};
        minimum_residue_coarse_slack_by_p.fill(
            std::numeric_limits<int>::max()
        );
        minimum_residue_coarse_boundary_slack_by_p.fill(
            std::numeric_limits<int>::max()
        );
        std::array<std::tuple<int, int, int, int, int, int>, 31>
            minimum_residue_coarse_case_by_p{};
        std::array<std::tuple<int, int, int, int, int, int>, 31>
            minimum_active_main_case_by_p{};
        std::array<std::array<int, 31>, 3> minimum_category_main_by_p{};
        std::array<
            std::array<std::tuple<int, int, int, int, int, int>, 31>,
            3
        > minimum_category_main_case_by_p{};
        for (auto& category : minimum_category_main_by_p) {
            category.fill(std::numeric_limits<int>::max());
        }
        minimum_active_main_by_p.fill(std::numeric_limits<int>::max());
        minimum_active_main_slack_by_p.fill(
            std::numeric_limits<int>::max()
        );
        minimum_active_boundary_slack_by_p.fill(
            std::numeric_limits<int>::max()
        );
        std::tuple<int, int, int, int, int, int> total_case{};
        std::tuple<int, int, int, int, int, int> main_case{};
        std::tuple<int, int, int, int, int, int> disjoint_total_case{};
        std::tuple<int, int, int, int, int, int> disjoint_main_case{};
        std::tuple<int, int, int, int, int, int> disjoint_negative_case{};
        std::tuple<int, int, int, int, int, int>
            disjoint_main_singleton_case{};
        std::tuple<int, int, int, int, int, int> disjoint_endpoint_case{};
        std::tuple<int, int, int, int, int, int> disjoint_boundary_case{};
        std::tuple<int, int, int, int, int, int>
            active_disjoint_boundary_case{};
        std::tuple<int, int, int, int, int, int> amplification_case{};
        std::tuple<int, int, int, int, int, int>
            disjoint_interior_main_singleton_case{};
        std::tuple<int, int, int, int, int, int>
            disjoint_strict_interior_main_case{};
        std::tuple<int, int, int, int, int, int>
            active_disjoint_interior_main_singleton_case{};
        std::tuple<int, int, int, int, int, int>
            active_disjoint_strict_interior_main_case{};
        std::array<std::tuple<int, int, int, int, int, int>, 3> shell_case{};
        std::array<std::tuple<int, int, int, int, int, int>, 2> parity_case{};
        std::array<std::tuple<int, int, int, int, int, int>, 5> prefix_case{};
        for (int p = 2; p <= maximum_label; ++p) {
            for (int q = p; q <= maximum_label; ++q) {
                for (int r = q; r <= maximum_label; ++r) {
                    for (int s = r; s <= maximum_label; ++s) {
                        for (int t = s; t <= maximum_label; ++t) {
                            const std::array<int, 4> suffix{q, r, s, t};
                            const int maximum_target = p + q + r + s + t;
                            for (int target = 1;
                                 target <= maximum_target;
                                 ++target) {
                                const Decomposition current = decomposition(
                                    p, suffix, target
                                );
                                int total = 0;
                                for (int layer : current.layers) {
                                    total += layer;
                                }
                                const int main = multiplicity(
                                    {p - 1, q, r, s, t}, target
                                );
                                const int main_slack
                                    = main - current.negative;
                                const auto current_case
                                    = std::tuple{p, q, r, s, t, target};
                                const bool disjoint = target != p
                                    && target != q && target != r
                                    && target != s && target != t;
                                const std::size_t p_index
                                    = static_cast<std::size_t>(p);
                                const std::array<int, 3> shells{
                                    current.layers[0U] + current.layers[4U],
                                    current.layers[1U] + current.layers[3U],
                                    current.layers[2U]
                                };
                                const std::array<int, 2> parity_blocks{
                                    current.layers[0U] + current.layers[2U]
                                        + current.layers[4U],
                                    current.layers[1U] + current.layers[3U]
                                };
                                ++tested;
                                if (total < minimum_total) {
                                    minimum_total = total;
                                    total_case = current_case;
                                }
                                if (main_slack < minimum_main_slack) {
                                    minimum_main_slack = main_slack;
                                    main_case = current_case;
                                }
                                if (disjoint) {
                                    ++disjoint_tested;
                                    if (total < minimum_disjoint_total) {
                                        minimum_disjoint_total = total;
                                        disjoint_total_case = current_case;
                                    }
                                    if (main_slack
                                        < minimum_disjoint_main_slack) {
                                        minimum_disjoint_main_slack
                                            = main_slack;
                                        disjoint_main_case = current_case;
                                    }
                                    if (current.negative
                                        > maximum_disjoint_negative) {
                                        maximum_disjoint_negative
                                            = current.negative;
                                        disjoint_negative_case = current_case;
                                    }
                                    const int main_singleton_slack = main
                                        + current.layers[1U]
                                        - current.negative;
                                    if (main_singleton_slack
                                        < minimum_disjoint_main_singleton_slack) {
                                        minimum_disjoint_main_singleton_slack
                                            = main_singleton_slack;
                                        disjoint_main_singleton_case
                                            = current_case;
                                    }
                                    const int boundary_slack = main
                                        + current.layers[1U]
                                        + current.layers[4U]
                                        - current.negative;
                                    const int endpoint_slack = main
                                        + current.layers[4U]
                                        - current.negative;
                                    if (endpoint_slack
                                        < minimum_disjoint_endpoint_slack) {
                                        minimum_disjoint_endpoint_slack
                                            = endpoint_slack;
                                        disjoint_endpoint_case = current_case;
                                    }
                                    if (boundary_slack
                                        < minimum_disjoint_boundary_slack) {
                                        minimum_disjoint_boundary_slack
                                            = boundary_slack;
                                        disjoint_boundary_case = current_case;
                                    }
                                    if (current.negative > 0
                                        && boundary_slack
                                            < minimum_active_disjoint_boundary_slack) {
                                        minimum_active_disjoint_boundary_slack
                                            = boundary_slack;
                                        active_disjoint_boundary_case
                                            = current_case;
                                    }
                                    if (current.negative > 0) {
                                        minimum_active_main_by_p[p_index]
                                            = std::min(
                                                minimum_active_main_by_p[
                                                    p_index],
                                                main
                                            );
                                        if (main
                                            == minimum_active_main_by_p[
                                                p_index]) {
                                            minimum_active_main_case_by_p[
                                                p_index] = current_case;
                                        }
                                        maximum_active_negative_by_p[p_index]
                                            = std::max(
                                                maximum_active_negative_by_p[
                                                    p_index],
                                                current.negative
                                            );
                                        minimum_active_main_slack_by_p[p_index]
                                            = std::min(
                                                minimum_active_main_slack_by_p[
                                                    p_index],
                                                main_slack
                                            );
                                        minimum_active_boundary_slack_by_p[
                                            p_index] = std::min(
                                                minimum_active_boundary_slack_by_p[
                                                    p_index],
                                                boundary_slack
                                            );
                                    }
                                    const std::array<int, 3> categories{
                                        current.first_negative_layers[2U],
                                        current.second_negative_layers[2U],
                                        current.second_negative_layers[3U]
                                    };
                                    for (std::size_t category = 0;
                                         category < categories.size();
                                         ++category) {
                                        if (categories[category] > 0) {
                                            if (main
                                                < minimum_category_main_by_p[
                                                    category][p_index]) {
                                                minimum_category_main_by_p[
                                                    category][p_index] = main;
                                                minimum_category_main_case_by_p[
                                                    category][p_index]
                                                    = current_case;
                                            }
                                        }
                                    }
                                    if (categories[1U] > 0
                                        || categories[2U] > 0) {
                                        const int amplification_slack = main
                                            - p * (p + 2);
                                        if (amplification_slack
                                            < minimum_amplification_slack) {
                                            minimum_amplification_slack
                                                = amplification_slack;
                                            amplification_case = current_case;
                                        }
                                        if (amplification_slack < 0) {
                                            std::cout
                                                << "FAIL amplification_slack="
                                                << amplification_slack
                                                << " case=";
                                            print_case(current_case);
                                            std::cout << '\n';
                                            return EXIT_FAILURE;
                                        }
                                    }
                                    if (current.negative > p * (p + 2)) {
                                        ++residue_cases_by_p[p_index];
                                        if (q == p) {
                                            ++residue_equal_minimum_by_p[
                                                p_index];
                                        }
                                        if (target == 1
                                            || target == p - 1) {
                                            ++residue_endpoint_by_p[p_index];
                                        }
                                        const int width
                                            = (4 * p + 12 + p - 1) / p;
                                        const int coarse_slack
                                            = one_terminal_candidate_count(
                                                p, q, r, s, t, target, width
                                            ) - current.negative;
                                        if (coarse_slack
                                            < minimum_residue_coarse_slack_by_p[
                                                p_index]) {
                                            minimum_residue_coarse_slack_by_p[
                                                p_index] = coarse_slack;
                                            minimum_residue_coarse_case_by_p[
                                                p_index] = current_case;
                                        }
                                        const int coarse_boundary_slack
                                            = coarse_slack
                                            + current.positive_layers[1U]
                                            + current.positive_layers[4U];
                                        minimum_residue_coarse_boundary_slack_by_p[
                                            p_index] = std::min(
                                                minimum_residue_coarse_boundary_slack_by_p[
                                                    p_index],
                                                coarse_boundary_slack
                                            );
                                    }
                                    if (target != 1 && target != p - 1
                                        && main_singleton_slack
                                            < minimum_disjoint_interior_main_singleton_slack) {
                                        minimum_disjoint_interior_main_singleton_slack
                                            = main_singleton_slack;
                                        disjoint_interior_main_singleton_case
                                            = current_case;
                                    }
                                    if (current.negative > 0 && target != 1
                                        && target != p - 1
                                        && main_singleton_slack
                                            < minimum_active_disjoint_interior_main_singleton_slack) {
                                        minimum_active_disjoint_interior_main_singleton_slack
                                            = main_singleton_slack;
                                        active_disjoint_interior_main_singleton_case
                                            = current_case;
                                    }
                                    if (q > p && target != 1
                                        && target != p - 1
                                        && main_slack
                                            < minimum_disjoint_strict_interior_main_slack) {
                                        minimum_disjoint_strict_interior_main_slack
                                            = main_slack;
                                        disjoint_strict_interior_main_case
                                            = current_case;
                                    }
                                    if (current.negative > 0 && q > p
                                        && target != 1
                                        && target != p - 1
                                        && main_slack
                                            < minimum_active_disjoint_strict_interior_main_slack) {
                                        minimum_active_disjoint_strict_interior_main_slack
                                            = main_slack;
                                        active_disjoint_strict_interior_main_case
                                            = current_case;
                                    }
                                }
                                for (int shell = 0; shell < 3; ++shell) {
                                    if (shells[static_cast<std::size_t>(shell)]
                                        < minimum_shell[
                                            static_cast<std::size_t>(shell)]) {
                                        minimum_shell[
                                            static_cast<std::size_t>(shell)]
                                            = shells[static_cast<std::size_t>(
                                                shell
                                            )];
                                        shell_case[
                                            static_cast<std::size_t>(shell)]
                                            = current_case;
                                    }
                                }
                                for (int parity = 0; parity < 2; ++parity) {
                                    const int value = parity_blocks[
                                        static_cast<std::size_t>(parity)];
                                    if (value < minimum_parity[
                                            static_cast<std::size_t>(parity)]) {
                                        minimum_parity[
                                            static_cast<std::size_t>(parity)]
                                            = value;
                                        parity_case[
                                            static_cast<std::size_t>(parity)]
                                            = current_case;
                                    }
                                }
                                int prefix = 0;
                                for (int layer = 0; layer <= 4; ++layer) {
                                    prefix += current.layers[
                                        static_cast<std::size_t>(layer)
                                    ];
                                    if (prefix < minimum_prefix[
                                            static_cast<std::size_t>(layer)]) {
                                        minimum_prefix[
                                            static_cast<std::size_t>(layer)]
                                            = prefix;
                                        prefix_case[
                                            static_cast<std::size_t>(layer)]
                                            = current_case;
                                    }
                                }
                                if (total < 0) {
                                    std::cout << "FAIL total=" << total
                                              << " case=";
                                    print_case(current_case);
                                    std::cout << '\n';
                                    return EXIT_FAILURE;
                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << "SU2_OPD_FOUR_SUFFIX PASS tested=" << tested
                  << " maximum_label=" << maximum_label
                  << " minimum_total=" << minimum_total << " total_case=";
        print_case(total_case);
        std::cout << " minimum_main_slack=" << minimum_main_slack
                  << " main_case=";
        print_case(main_case);
        std::cout << '\n';
        std::cout << "disjoint_tested=" << disjoint_tested
                  << " minimum_disjoint_total=" << minimum_disjoint_total
                  << " disjoint_total_case=";
        print_case(disjoint_total_case);
        std::cout << " minimum_disjoint_main_slack="
                  << minimum_disjoint_main_slack
                  << " disjoint_main_case=";
        print_case(disjoint_main_case);
        std::cout << " minimum_disjoint_main_singleton_slack="
                  << minimum_disjoint_main_singleton_slack
                  << " disjoint_main_singleton_case=";
        print_case(disjoint_main_singleton_case);
        std::cout << " minimum_disjoint_endpoint_slack="
                  << minimum_disjoint_endpoint_slack
                  << " disjoint_endpoint_case=";
        print_case(disjoint_endpoint_case);
        std::cout << " minimum_disjoint_boundary_slack="
                  << minimum_disjoint_boundary_slack
                  << " disjoint_boundary_case=";
        print_case(disjoint_boundary_case);
        std::cout << " minimum_active_disjoint_boundary_slack="
                  << minimum_active_disjoint_boundary_slack
                  << " active_disjoint_boundary_case=";
        print_case(active_disjoint_boundary_case);
        std::cout << " minimum_amplification_slack="
                  << minimum_amplification_slack
                  << " amplification_case=";
        print_case(amplification_case);
        std::cout << " minimum_disjoint_interior_main_singleton_slack="
                  << minimum_disjoint_interior_main_singleton_slack
                  << " disjoint_interior_main_singleton_case=";
        print_case(disjoint_interior_main_singleton_case);
        std::cout << " minimum_disjoint_strict_interior_main_slack="
                  << minimum_disjoint_strict_interior_main_slack
                  << " disjoint_strict_interior_main_case=";
        print_case(disjoint_strict_interior_main_case);
        std::cout
            << " minimum_active_disjoint_interior_main_singleton_slack="
            << minimum_active_disjoint_interior_main_singleton_slack
            << " active_disjoint_interior_main_singleton_case=";
        print_case(active_disjoint_interior_main_singleton_case);
        std::cout << " minimum_active_disjoint_strict_interior_main_slack="
                  << minimum_active_disjoint_strict_interior_main_slack
                  << " active_disjoint_strict_interior_main_case=";
        print_case(active_disjoint_strict_interior_main_case);
        std::cout << " maximum_disjoint_negative="
                  << maximum_disjoint_negative
                  << " disjoint_negative_case=";
        print_case(disjoint_negative_case);
        std::cout << '\n';
        for (int layer = 0; layer <= 4; ++layer) {
            std::cout << "prefix_through_layer=" << layer
                      << " minimum="
                      << minimum_prefix[static_cast<std::size_t>(layer)]
                      << " case=";
            print_case(prefix_case[static_cast<std::size_t>(layer)]);
            std::cout << '\n';
        }
        for (int shell = 0; shell < 3; ++shell) {
            std::cout << "shell=" << shell << " minimum="
                      << minimum_shell[static_cast<std::size_t>(shell)]
                      << " case=";
            print_case(shell_case[static_cast<std::size_t>(shell)]);
            std::cout << '\n';
        }
        for (int parity = 0; parity < 2; ++parity) {
            std::cout << "parity_block=" << parity << " minimum="
                      << minimum_parity[static_cast<std::size_t>(parity)]
                      << " case=";
            print_case(parity_case[static_cast<std::size_t>(parity)]);
            std::cout << '\n';
        }
        for (int p = 2; p <= maximum_label; ++p) {
            const std::size_t p_index = static_cast<std::size_t>(p);
            if (minimum_active_main_by_p[p_index]
                == std::numeric_limits<int>::max()) {
                continue;
            }
            std::cout << "p=" << p
                      << " minimum_active_main="
                      << minimum_active_main_by_p[p_index]
                      << " minimum_active_main_case=";
            print_case(minimum_active_main_case_by_p[p_index]);
            std::cout
                      << " maximum_active_negative="
                      << maximum_active_negative_by_p[p_index]
                      << " minimum_active_main_slack="
                      << minimum_active_main_slack_by_p[p_index]
                      << " minimum_active_boundary_slack="
                      << minimum_active_boundary_slack_by_p[p_index]
                      << " residue_cases="
                      << residue_cases_by_p[p_index]
                      << " residue_equal_minimum="
                      << residue_equal_minimum_by_p[p_index]
                      << " residue_endpoint="
                      << residue_endpoint_by_p[p_index]
                      << " minimum_residue_coarse_slack="
                      << minimum_residue_coarse_slack_by_p[p_index]
                      << " minimum_residue_coarse_case=";
            print_case(minimum_residue_coarse_case_by_p[p_index]);
            std::cout
                      << " minimum_residue_coarse_boundary_slack="
                      << minimum_residue_coarse_boundary_slack_by_p[p_index]
                      << " minimum_category_main=";
            for (std::size_t category = 0; category < 3; ++category) {
                if (category != 0U) {
                    std::cout << ',';
                }
                const int value
                    = minimum_category_main_by_p[category][p_index];
                if (value == std::numeric_limits<int>::max()) {
                    std::cout << "none";
                } else {
                    std::cout << value << '@';
                    print_case(minimum_category_main_case_by_p[
                        category][p_index]);
                }
            }
            std::cout << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
