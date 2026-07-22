#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Wedge = std::pair<int, int>;
using State = std::map<Wedge, cpp_int>;

void add_wedge(
    State& state,
    const int first,
    const int second,
    const cpp_int& value
) {
    if (first == second || value == 0) {
        return;
    }
    if (first > second) {
        state[{first, second}] += value;
        if (state[{first, second}] == 0) {
            state.erase({first, second});
        }
    } else {
        state[{second, first}] -= value;
        if (state[{second, first}] == 0) {
            state.erase({second, first});
        }
    }
}

State update(const State& state, const int label) {
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

State source_curl_seed(const int q, const int r) {
    State answer;
    State fundamental;
    add_wedge(fundamental, 1, 0, 1);
    for (int output = std::abs(q - r);
         output <= q + r;
         output += 2) {
        const State reservoir = update(fundamental, output);
        for (const auto& [wedge, coefficient] : reservoir) {
            add_wedge(
                answer, wedge.first, wedge.second, coefficient
            );
        }
    }
    for (int output = std::abs(q - 1);
         output <= q + 1;
         output += 2) {
        add_wedge(answer, output, r, 1);
    }
    for (int output = std::abs(r - 1);
         output <= r + 1;
         output += 2) {
        add_wedge(answer, output, q, 1);
    }
    return answer;
}

cpp_int twisted_coefficient(const Wedge& wedge, const cpp_int& coefficient) {
    // For opposite parity, central translation sends
    //   e_u wedge e_v  -> (-1)^v (e_u tensor e_v + e_v tensor e_u).
    return (wedge.second & 1) == 0 ? coefficient : -coefficient;
}

cpp_int twisted_at(const State& state, const int odd, const int even) {
    const Wedge wedge{
        std::max(odd, even), std::min(odd, even)
    };
    const auto found = state.find(wedge);
    if (found == state.end()) {
        return 0;
    }
    return twisted_coefficient(wedge, found->second);
}

void print_word(const std::vector<int>& word) {
    std::cout << '[';
    for (std::size_t index = 0; index < word.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << word[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4 && argc != 5) {
            throw std::runtime_error(
                "usage: search_su2_two_odd_twisted_cone "
                "MAXIMUM_ODD_LABEL MAXIMUM_EVEN_LABEL MAXIMUM_EVEN_FACTORS "
                "[--dominance]"
            );
        }
        const int maximum_odd = std::stoi(argv[1]);
        const int maximum_even = std::stoi(argv[2]);
        const int maximum_factors = std::stoi(argv[3]);
        const bool dominance_mode = argc == 5
            && std::string(argv[4]) == "--dominance";
        if (argc == 5 && !dominance_mode) {
            throw std::runtime_error("unknown mode");
        }
        if (maximum_odd < 1 || maximum_even < 2 || maximum_factors < 0) {
            throw std::runtime_error("invalid search bound");
        }
        long long states = 0;
        long long source_identities = 0;
        for (int q = 3; q <= maximum_odd; q += 2) {
            for (int r = q; r <= maximum_odd; r += 2) {
                State direct;
                add_wedge(direct, 1, 0, 1);
                direct = update(direct, q);
                direct = update(direct, r);
                if (direct != source_curl_seed(q, r)) {
                    throw std::runtime_error(
                        "source-curl seed identity failed"
                    );
                }
                ++source_identities;
            }
        }
        for (int q = 3; q <= maximum_odd; q += 2) {
            for (int r = q; r <= maximum_odd; r += 2) {
                State seed;
                add_wedge(seed, 1, 0, 1);
                seed = update(seed, q);
                seed = update(seed, r);
                std::vector<int> word;
                const auto visit = [&](const auto& self, const int next_even) -> void {
                    ++states;
                    State current = seed;
                    for (const int label : word) {
                        current = update(current, label);
                    }
                    if (dominance_mode) {
                        int support = q + r + 1;
                        for (const int label : word) {
                            support += label;
                        }
                        for (int a = 1; a <= support; a += 2) {
                            const cpp_int boundary = twisted_at(current, a, 0);
                            if (boundary < 0) {
                                std::cout << "BOUNDARY_FAIL q=" << q
                                          << " r=" << r << " even_word=";
                                print_word(word);
                                std::cout << " a=" << a << " coefficient="
                                          << boundary << '\n';
                                std::exit(EXIT_SUCCESS);
                            }
                            for (int b = 2; b <= support; b += 2) {
                                cpp_int value = twisted_at(current, a, b);
                                for (int output = std::abs(a - b);
                                     output <= a + b;
                                     output += 2) {
                                    value += twisted_at(current, output, 0);
                                }
                                if (value < 0) {
                                    std::cout << "INTERVAL_DOMINANCE_FAIL q="
                                              << q << " r=" << r
                                              << " even_word=";
                                    print_word(word);
                                    std::cout << " a=" << a << " b=" << b
                                              << " value=" << value << '\n';
                                    std::exit(EXIT_SUCCESS);
                                }
                            }
                        }
                    } else if (!word.empty()) {
                        for (const auto& [wedge, coefficient] : current) {
                            const cpp_int value = twisted_coefficient(
                                wedge, coefficient
                            );
                            if (value < 0) {
                                std::cout << "TWISTED_CONE_FAIL q=" << q
                                          << " r=" << r << " even_word=";
                                print_word(word);
                                std::cout << " wedge=(" << wedge.first << ','
                                          << wedge.second << ") coefficient="
                                          << value << " source_identities="
                                          << source_identities << '\n';
                                std::exit(EXIT_SUCCESS);
                            }
                        }
                    }
                    if (word.size()
                        == static_cast<std::size_t>(maximum_factors)) {
                        return;
                    }
                    for (int label = next_even;
                         label <= maximum_even;
                         label += 2) {
                        word.push_back(label);
                        self(self, label);
                        word.pop_back();
                    }
                };
                visit(visit, 2);
            }
        }
        std::cout << (dominance_mode
                          ? "SU2_TWO_ODD_INTERVAL_DOMINANCE PASS states="
                          : "SU2_TWO_ODD_TWISTED_CONE PASS states=")
                  << states
                  << " source_identities=" << source_identities
                  << " maximum_odd_label=" << maximum_odd
                  << " maximum_even_label=" << maximum_even
                  << " maximum_even_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
