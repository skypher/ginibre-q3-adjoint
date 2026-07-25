#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <atomic>
#include <iostream>
#include <map>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;
using Key = std::pair<int, int>;
using State = std::map<Key, cpp_int>;

static std::vector<int> outputs(int left, int label, int level) {
    const int low = std::abs(left - label);
    const int high = level < 0
        ? left + label
        : std::min(left + label, 2 * level - left - label);
    std::vector<int> result;
    for (int value = low; value <= high; value += 2) {
        result.push_back(value);
    }
    return result;
}

static State step(const State& input, int label, int sign, int level) {
    State output;
    for (const auto& [key, coefficient] : input) {
        const auto [left, right] = key;
        for (int value : outputs(left, label, level)) {
            output[{value, right}] += coefficient;
        }
        for (int value : outputs(right, label, level)) {
            output[{left, value}] += sign * coefficient;
        }
    }
    for (auto iterator = output.begin(); iterator != output.end();) {
        if (iterator->second == 0) {
            iterator = output.erase(iterator);
        } else {
            ++iterator;
        }
    }
    return output;
}

static State word(int adjoint_minus, int fundamental_minus, int fundamental_plus, int level) {
    State state{{{0, 0}, 1}};
    for (int index = 0; index < adjoint_minus; ++index) {
        state = step(state, 2, -1, level);
    }
    for (int index = 0; index < fundamental_minus; ++index) {
        state = step(state, 1, -1, level);
    }
    for (int index = 0; index < fundamental_plus; ++index) {
        state = step(state, 1, +1, level);
    }
    return state;
}

static cpp_int coefficient(const State& state, int left, int right) {
    const auto iterator = state.find({left, right});
    return iterator == state.end() ? cpp_int(0) : iterator->second;
}

int main() {
    std::vector<std::tuple<int, int, int>> tasks;
    for (int adjoint_minus = 0; adjoint_minus <= 12; ++adjoint_minus) {
        for (int fundamental_minus = 0; fundamental_minus <= 12; ++fundamental_minus) {
            for (int fundamental_plus = 0; fundamental_plus <= 12; ++fundamental_plus) {
                if ((adjoint_minus + fundamental_minus) % 2 == 0) {
                    tasks.emplace_back(adjoint_minus, fundamental_minus, fundamental_plus);
                }
            }
        }
    }
    if (tasks.size() != 1105U) {
        throw std::runtime_error("unexpected task count");
    }

    std::atomic<std::size_t> next{0};
    std::atomic<bool> failed{false};
    auto worker = [&]() {
        while (true) {
            const std::size_t index = next.fetch_add(1);
            if (index >= tasks.size()) {
                break;
            }
            const auto [adjoint_minus, fundamental_minus, fundamental_plus] = tasks[index];
            const int total = 2 * adjoint_minus + fundamental_minus + fundamental_plus;
            int level = std::max(3, (total + 2) / 2);
            if (level % 2 == 0) {
                ++level;
            }
            const State ordinary = word(
                adjoint_minus,
                fundamental_minus,
                fundamental_plus,
                -1
            );
            const State finite = word(
                adjoint_minus,
                fundamental_minus,
                fundamental_plus,
                level
            );
            const cpp_int corner = coefficient(ordinary, 0, 0);
            const cpp_int partial = coefficient(ordinary, 1, 0);
            if (
                corner < 0 ||
                partial < 0 ||
                corner != coefficient(finite, 0, 0) ||
                partial != coefficient(finite, 1, 0)
            ) {
                failed.store(true);
                break;
            }
        }
    };

    std::vector<std::thread> workers;
    for (int index = 0; index < 4; ++index) {
        workers.emplace_back(worker);
    }
    for (auto& thread : workers) {
        thread.join();
    }
    if (failed.load()) {
        throw std::runtime_error("regression failure");
    }

    std::cout
        << "SU2_V1_V2_MINUS_SECTOR PASS"
        << " cases=" << tasks.size()
        << " max_exponent=12 threads=4"
        << " scalar=1 partial_v1=1 stable_transfer=1\n";
}
