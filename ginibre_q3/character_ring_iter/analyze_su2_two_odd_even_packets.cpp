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
    const Wedge wedge{
        std::max(first, second), std::min(first, second)
    };
    state[wedge] += first > second ? value : -value;
    if (state[wedge] == 0) {
        state.erase(wedge);
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

void print_state(const State& state) {
    bool first = true;
    std::cout << '{';
    for (const auto& [wedge, coefficient] : state) {
        if (!first) {
            std::cout << ',';
        }
        first = false;
        std::cout << wedge.first << '^' << wedge.second << ':' << coefficient;
    }
    std::cout << '}';
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

void enumerate_words(
    const int remaining,
    const int next,
    std::vector<int>& word,
    std::vector<std::vector<int>>& output
) {
    if (remaining == 0) {
        output.push_back(word);
        return;
    }
    for (int label = next; label <= remaining; label += 2) {
        word.push_back(label);
        enumerate_words(remaining - label, label, word, output);
        word.pop_back();
    }
}

State packet(const bool outer, const int p, const std::vector<int>& word) {
    State state;
    if (outer) {
        add_wedge(state, p + 1, 0, 1);
    } else {
        add_wedge(state, p - 1, 0, 1);
        add_wedge(state, p, 1, -1);
    }
    for (const int label : word) {
        state = update(state, label);
    }
    return state;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_two_odd_even_packets Q R"
            );
        }
        const int q = std::stoi(argv[1]);
        const int r = std::stoi(argv[2]);
        if (q < 1 || r < 1 || (q & 1) == 0 || (r & 1) == 0) {
            throw std::runtime_error("Q and R must be positive odd labels");
        }
        State target;
        add_wedge(target, 1, 0, 1);
        target = update(update(target, q), r);
        std::cout << "TARGET q=" << q << " r=" << r << ' ';
        print_state(target);
        std::cout << '\n';

        const int total = q + r;
        for (int p = 2; p <= total; p += 2) {
            std::vector<std::vector<int>> words;
            std::vector<int> word;
            enumerate_words(total - p, 2, word, words);
            for (const auto& suffix : words) {
                const State outer = packet(true, p, suffix);
                std::cout << "X p=" << p << " Q=";
                print_word(suffix);
                std::cout << ' ';
                print_state(outer);
                std::cout << '\n';
                if (suffix.empty() || suffix.front() >= p) {
                    const State ordered = packet(false, p, suffix);
                    std::cout << "Y p=" << p << " Q=";
                    print_word(suffix);
                    std::cout << ' ';
                    print_state(ordered);
                    std::cout << '\n';
                }
            }
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
