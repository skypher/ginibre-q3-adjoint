#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Wedge = std::pair<int, int>;
using State = std::map<Wedge, cpp_int>;

void add_wedge(State& state, int first, int second, const cpp_int& value) {
    if (first == second || value == 0) {
        return;
    }
    if (first > second) {
        state[{first, second}] += value;
    } else {
        state[{second, first}] -= value;
    }
}

State update(const State& state, int label) {
    State result;
    for (const auto& [wedge, coefficient] : state) {
        const auto [first, second] = wedge;
        for (int output = std::abs(first - label);
             output <= first + label;
             output += 2) {
            add_wedge(result, output, second, coefficient);
        }
        for (int output = std::abs(second - label);
             output <= second + label;
             output += 2) {
            add_wedge(result, first, output, coefficient);
        }
    }
    return result;
}

cpp_int coefficient(const State& state, int first, int second) {
    if (first == second) {
        return 0;
    }
    const bool reversed = first < second;
    const Wedge key = reversed ? Wedge{second, first} : Wedge{first, second};
    const auto found = state.find(key);
    if (found == state.end()) {
        return 0;
    }
    return reversed ? -found->second : found->second;
}

void print_word(int p, const std::vector<int>& suffix, int q) {
    std::cout << "p=" << p << " prefix=[";
    for (std::size_t index = 0; index < suffix.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << suffix[index];
    }
    std::cout << "] q=" << q;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: verify_su2_opd_dominant_suffix MAXIMUM_PREFIX_LABEL "
                "MAXIMUM_PREFIX_FACTORS"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 2 || maximum_factors < 0) {
            throw std::runtime_error("invalid bound");
        }
        long long words = 0;
        long long extensions = 0;
        long long boundary_cases = 0;
        for (int p = 2; p <= maximum_label; ++p) {
            std::vector<int> suffix;
            State initial;
            add_wedge(initial, p - 1, 0, 1);
            add_wedge(initial, p, 1, -1);
            const auto visit = [&](const auto& self, int next_label,
                                   int total, const State& state) -> void {
                ++words;
                const int maximum_support = p + 1 + total;
                for (int target = 1; target <= maximum_support; ++target) {
                    if (coefficient(state, target, 0) < 0) {
                        std::cout << "FAIL prefix-boundary target=" << target
                                  << ' ';
                        print_word(p, suffix, next_label);
                        std::cout << '\n';
                        std::exit(EXIT_FAILURE);
                    }
                }
                const int first_q = std::max(next_label, p + total - 3);
                const int last_q = std::max(first_q, p + total + 2);
                for (int q = first_q; q <= last_q; ++q) {
                    ++extensions;
                    const int count_twos = static_cast<int>(std::ranges::count(
                        suffix, 2
                    ));
                    const int count_threes = static_cast<int>(std::ranges::count(
                        suffix, 3
                    ));
                    for (int target = 1;
                         target <= maximum_support + 1;
                         ++target) {
                        cpp_int expected = 0;
                        if (q == p + total && target == 1) {
                            expected = 1;
                        } else if (q == p + total - 1 && p == 2
                                   && target == 2) {
                            expected = -1;
                        } else if (q == p + total - 2 && target == 1) {
                            expected = static_cast<int>(suffix.size());
                            if (p == 2) {
                                expected += 1 + count_twos;
                            }
                        } else if (q == p + total - 2 && target == 3) {
                            if (p == 2) {
                                expected = count_twos;
                            } else if (p == 3) {
                                expected = -1;
                            }
                        } else if (q == p + total - 3 && target == 2) {
                            if (p == 2) {
                                expected = count_threes - 2 * count_twos
                                    - static_cast<int>(suffix.size());
                            } else if (p == 3) {
                                expected = 1 + count_threes;
                            }
                        } else if (q == p + total - 3 && target == 4) {
                            if (p == 2) {
                                expected = count_threes - count_twos;
                            } else if (p == 3) {
                                expected = count_threes;
                            } else if (p == 4) {
                                expected = -1;
                            }
                        }
                        const cpp_int actual = coefficient(state, target, q);
                        if (actual != expected) {
                            std::cout << "FAIL carrier target=" << target
                                      << " actual=" << actual
                                      << " expected=" << expected << ' ';
                            print_word(p, suffix, q);
                            std::cout << '\n';
                            std::exit(EXIT_FAILURE);
                        }
                    }
                    if (q == p + total - 1 && p == 2
                        && coefficient(state, q, 0) != 1) {
                        std::cout << "FAIL reservoir value="
                                  << coefficient(state, q, 0) << ' ';
                        print_word(p, suffix, q);
                        std::cout << '\n';
                        std::exit(EXIT_FAILURE);
                    }
                    if (q == p + total - 2 && p == 3
                        && coefficient(state, q + 1, 0) != 1) {
                        std::cout << "FAIL second-reservoir value="
                                  << coefficient(state, q + 1, 0) << ' ';
                        print_word(p, suffix, q);
                        std::cout << '\n';
                        std::exit(EXIT_FAILURE);
                    }
                    if (q == p + total - 3 && p == 2) {
                        const cpp_int factor_count
                            = static_cast<int>(suffix.size());
                        const cpp_int n_two = count_twos;
                        const cpp_int n_three = count_threes;
                        const cpp_int expected_lower
                            = factor_count * (factor_count + 1) / 2 - 1
                            + n_two * (n_two - 1)
                            + n_two * (factor_count - 1)
                            - n_two * n_three;
                        if (coefficient(state, q + 2, 0) != 1
                            || coefficient(state, q, 0)
                                != factor_count + n_two
                            || coefficient(state, q - 2, 0)
                                != expected_lower) {
                            std::cout << "FAIL third-reservoir values="
                                      << coefficient(state, q + 2, 0) << ','
                                      << coefficient(state, q, 0) << ','
                                      << coefficient(state, q - 2, 0) << ' ';
                            print_word(p, suffix, q);
                            std::cout << '\n';
                            std::exit(EXIT_FAILURE);
                        }
                    }
                    if (q == p + total - 3 && p == 4
                        && coefficient(state, q + 2, 0) != 1) {
                        std::cout << "FAIL fourth-reservoir value="
                                  << coefficient(state, q + 2, 0) << ' ';
                        print_word(p, suffix, q);
                        std::cout << '\n';
                        std::exit(EXIT_FAILURE);
                    }
                    const State extended = update(state, q);
                    const int extended_support = maximum_support + q;
                    for (int target = 1;
                         target <= extended_support;
                         ++target) {
                        ++boundary_cases;
                        if (coefficient(extended, target, 0) < 0) {
                            std::cout << "FAIL extended-boundary target="
                                      << target << " value="
                                      << coefficient(extended, target, 0)
                                      << ' ';
                            print_word(p, suffix, q);
                            std::cout << '\n';
                            std::exit(EXIT_FAILURE);
                        }
                    }
                }
                if (suffix.size()
                    == static_cast<std::size_t>(maximum_factors)) {
                    return;
                }
                for (int label = next_label; label <= maximum_label; ++label) {
                    suffix.push_back(label);
                    const State next = update(state, label);
                    self(self, label, total + label, next);
                    suffix.pop_back();
                }
            };
            visit(visit, p, 0, initial);
        }
        std::cout << "SU2_OPD_DOMINANT_SUFFIX PASS words=" << words
                  << " extensions=" << extensions
                  << " boundary_cases=" << boundary_cases
                  << " maximum_prefix_label=" << maximum_label
                  << " maximum_prefix_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
