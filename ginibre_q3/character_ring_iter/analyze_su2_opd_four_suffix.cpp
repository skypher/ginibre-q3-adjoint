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
        const int negative = multiplicity(second, p - 1)
                * multiplicity(first, target)
            + multiplicity(second, 1)
                * multiplicity(first_p, target);
        const int layer = std::popcount(mask);
        result.layers[static_cast<std::size_t>(layer)]
            += positive - negative;
        result.positive += positive;
        result.negative += negative;
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
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_opd_four_suffix MAXIMUM_LABEL"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        if (maximum_label < 2 || maximum_label > 30) {
            throw std::runtime_error("invalid maximum label");
        }
        long long tested = 0;
        int minimum_total = std::numeric_limits<int>::max();
        int minimum_main_slack = std::numeric_limits<int>::max();
        std::array<int, 5> minimum_prefix{
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max(),
            std::numeric_limits<int>::max()
        };
        std::tuple<int, int, int, int, int, int> total_case{};
        std::tuple<int, int, int, int, int, int> main_case{};
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
                                ++tested;
                                if (total < minimum_total) {
                                    minimum_total = total;
                                    total_case = current_case;
                                }
                                if (main_slack < minimum_main_slack) {
                                    minimum_main_slack = main_slack;
                                    main_case = current_case;
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
        for (int layer = 0; layer <= 4; ++layer) {
            std::cout << "prefix_through_layer=" << layer
                      << " minimum="
                      << minimum_prefix[static_cast<std::size_t>(layer)]
                      << " case=";
            print_case(prefix_case[static_cast<std::size_t>(layer)]);
            std::cout << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
