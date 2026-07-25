#include <boost/multiprecision/cpp_int.hpp>

#include <atomic>
#include <iostream>
#include <map>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;
using Key = std::pair<int, int>;
using State = std::map<Key, cpp_int>;

static State step(const State& input, int label, int sign) {
    State output;
    for (const auto& [key, coefficient] : input) {
        const auto [left, right] = key;
        for (int value = std::abs(left - label); value <= left + label; value += 2) {
            output[{value, right}] += coefficient;
        }
        for (int value = std::abs(right - label); value <= right + label; value += 2) {
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

static State word(int half_fundamental_minus, int adjoint_plus, bool transformed) {
    State state{{{0, 0}, 1}};
    const int fundamental_sign = transformed ? 1 : -1;
    for (int index = 0; index < 2 * half_fundamental_minus; ++index) {
        state = step(state, 1, fundamental_sign);
    }
    for (int index = 0; index < adjoint_plus; ++index) {
        state = step(state, 2, 1);
    }
    return state;
}

int main() {
    std::vector<std::pair<int, int>> tasks;
    for (int half_fundamental_minus = 0; half_fundamental_minus <= 12; ++half_fundamental_minus) {
        for (int adjoint_plus = 0; adjoint_plus <= 12; ++adjoint_plus) {
            tasks.emplace_back(half_fundamental_minus, adjoint_plus);
        }
    }

    std::atomic<std::size_t> next{0};
    std::atomic<bool> failed{false};
    auto worker = [&]() {
        while (true) {
            const std::size_t index = next.fetch_add(1);
            if (index >= tasks.size()) {
                break;
            }
            const auto [half_fundamental_minus, adjoint_plus] = tasks[index];
            const State original = word(half_fundamental_minus, adjoint_plus, false);
            const State all_plus = word(half_fundamental_minus, adjoint_plus, true);
            std::map<Key, cpp_int> keys = original;
            for (const auto& item : all_plus) {
                keys.emplace(item.first, 0);
            }
            for (const auto& item : keys) {
                const auto original_iterator = original.find(item.first);
                const auto positive_iterator = all_plus.find(item.first);
                cpp_int left = original_iterator == original.end()
                    ? cpp_int(0)
                    : original_iterator->second;
                cpp_int right = positive_iterator == all_plus.end()
                    ? cpp_int(0)
                    : positive_iterator->second;
                if ((item.first.second & 1) != 0) {
                    right = -right;
                }
                if (left != right) {
                    failed.store(true);
                    return;
                }
            }
        }
    };

    std::vector<std::thread> threads;
    for (int index = 0; index < 4; ++index) {
        threads.emplace_back(worker);
    }
    for (auto& thread : threads) {
        thread.join();
    }
    if (failed.load()) {
        throw std::runtime_error("parity automorphism regression failed");
    }

    std::cout
        << "SU2_V1_V2_PARITY_AUTOMORPHISM PASS"
        << " cases=" << tasks.size()
        << " max_half_minus=12 max_adjoint_plus=12"
        << " threads=4 full_tensor_identity=1\n";
}
