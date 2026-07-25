#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>

using boost::multiprecision::cpp_int;
using Key = std::pair<int, int>;
using State = std::map<Key, cpp_int>;

static State step(const State& input, int label, int sign, int rank) {
    State output;
    const int level = 2 * rank - 1;
    for (const auto& [key, coefficient] : input) {
        const auto [left, right] = key;
        int low = std::abs(left - label);
        int high = std::min(left + label, level - left - label);
        for (int result = low; result <= high; ++result) {
            output[{result, right}] += coefficient;
        }
        low = std::abs(right - label);
        high = std::min(right + label, level - right - label);
        for (int result = low; result <= high; ++result) {
            output[{left, result}] += sign * coefficient;
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

static State word(
    int rank,
    int first_minus,
    int top_minus,
    int top_plus,
    bool reduced
) {
    State state{{{0, 0}, 1}};
    const int top = rank - 1;
    if (reduced) {
        for (int index = 0; index < first_minus + top_minus; ++index) {
            state = step(state, top, -1, rank);
        }
        for (int index = 0; index < first_minus + top_plus; ++index) {
            state = step(state, top, +1, rank);
        }
    } else {
        for (int index = 0; index < first_minus; ++index) {
            state = step(state, 1, -1, rank);
        }
        for (int index = 0; index < top_minus; ++index) {
            state = step(state, top, -1, rank);
        }
        for (int index = 0; index < top_plus; ++index) {
            state = step(state, top, +1, rank);
        }
    }
    return state;
}

int main() {
    std::size_t cases = 0;
    for (int rank = 2; rank <= 12; ++rank) {
        for (int first_minus = 0; first_minus <= 6; ++first_minus) {
            for (int top_minus = 0; top_minus <= 6; ++top_minus) {
                for (int top_plus = 0; top_plus <= 6; ++top_plus) {
                    if ((first_minus + top_minus) % 2 != 0) {
                        continue;
                    }
                    const State original = word(
                        rank,
                        first_minus,
                        top_minus,
                        top_plus,
                        false
                    );
                    const State reduced = word(
                        rank,
                        first_minus,
                        top_minus,
                        top_plus,
                        true
                    );
                    if (original != reduced) {
                        throw std::runtime_error("top reduction mismatch");
                    }
                    ++cases;
                }
            }
        }
    }
    if (cases != 1925U) {
        throw std::runtime_error("unexpected case count");
    }
    std::cout
        << "SU2_ODD_ORBIT_B1_TOP_REDUCTION PASS"
        << " levels=11 cases=" << cases
        << " full_tensor_identity=1\n";
}
