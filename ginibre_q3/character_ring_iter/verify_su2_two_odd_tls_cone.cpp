#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>

namespace {

using Pair = std::pair<int, int>;

Pair ordered_pair(int first, int second) {
    return first >= second ? Pair{first, second} : Pair{second, first};
}

template <class Function>
void for_each_output(int left, int right, Function function) {
    for (int output = std::abs(left - right);
         output <= left + right;
         output += 2) {
        function(output);
    }
}

void add_weight(std::map<Pair, int>& weights, int first, int second, int value) {
    weights[ordered_pair(first, second)] += value;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: verify_su2_two_odd_tls_cone MAXIMUM_ODD_LABEL"
            );
        }
        const int maximum_odd = std::stoi(argv[1]);
        if (maximum_odd < 1) {
            throw std::runtime_error("invalid bound");
        }
        long long triples = 0;
        long long rays = 0;
        for (int a = 1; a <= maximum_odd; a += 2) {
            for (int q = 1; q <= maximum_odd; q += 2) {
                for (int r = q; r <= maximum_odd; r += 2) {
                    if (a + 1 < q + r) {
                        continue;
                    }
                    ++triples;
                    std::map<Pair, int> weights;
                    for_each_output(r, 1, [&](const int neighbor) {
                        for_each_output(a, q, [&](const int first) {
                            add_weight(weights, first, neighbor, 1);
                        });
                    });
                    for_each_output(q, 1, [&](const int neighbor) {
                        for_each_output(a, r, [&](const int first) {
                            add_weight(weights, first, neighbor, 1);
                        });
                    });
                    for_each_output(q, r, [&](const int second) {
                        for_each_output(a, 1, [&](const int first) {
                            add_weight(weights, first, second, -1);
                        });
                    });
                    const int maximum_index = a + q + r + 1;
                    for (int first = 0; first <= maximum_index; first += 2) {
                        for (int second = 2; second <= first; second += 2) {
                            ++rays;
                            int coefficient = 0;
                            for (int shift = 0;
                                 first + shift <= maximum_index;
                                 shift += 2) {
                                const auto found = weights.find(
                                    {first + shift, second + shift}
                                );
                                if (found != weights.end()) {
                                    coefficient += found->second;
                                }
                            }
                            if (coefficient < 0) {
                                std::cout << "FAIL a=" << a << " q=" << q
                                          << " r=" << r << " ray=("
                                          << first << ',' << second
                                          << ") coefficient=" << coefficient
                                          << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                    }
                }
            }
        }
        std::cout << "SU2_TWO_ODD_TLS_CONE PASS triples=" << triples
                  << " rays=" << rays
                  << " maximum_odd_label=" << maximum_odd << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
