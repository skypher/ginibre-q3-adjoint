#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>

namespace {

using Weight = std::pair<int, int>;
using Character = std::map<Weight, int>;

void add_scaled(Character& target, const Character& source, int scale) {
    for (const auto& [weight, coefficient] : source) {
        target[weight] += scale * coefficient;
    }
}

Character row_product(int first_d, int second_d) {
    Character answer;
    if (first_d <= 0 || second_d <= 0) {
        return answer;
    }
    int m = first_d - 1;
    int n = second_d - 1;
    if (m < n) {
        std::swap(m, n);
    }
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= i; ++j) {
            ++answer[{m + n - i - j, i - j}];
        }
    }
    return answer;
}

Character local_packet(int a, int b, int p) {
    Character answer;
    for (int output = std::abs(a - p); output <= a + p; output += 2) {
        add_scaled(answer, row_product(output, b), 1);
        add_scaled(answer, row_product(output - 2, b - 2), -1);
    }
    add_scaled(answer, row_product(a + b, p), -1);
    return answer;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: search_su2_sp4_local_packet MAXIMUM_LABEL"
            );
        }
        const int maximum = std::stoi(argv[1]);
        if (maximum < 1) {
            throw std::runtime_error("invalid maximum label");
        }
        for (int p = 1; p <= maximum; ++p) {
            for (int a = 2; a <= maximum; ++a) {
                if (a == p) {
                    continue;
                }
                for (int b = 2; b <= a; ++b) {
                    if (b == p) {
                        continue;
                    }
                    const Character packet = local_packet(a, b, p);
                    for (const auto& [weight, coefficient] : packet) {
                        if (coefficient < 0) {
                            std::cout << "SP4_LOCAL_PACKET FALSE a=" << a
                                      << " b=" << b << " p=" << p
                                      << " weight=(" << weight.first << ','
                                      << weight.second << ") coefficient="
                                      << coefficient << '\n';
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
        }
        std::cout << "SP4_LOCAL_PACKET maximum_label=" << maximum
                  << " CHARACTER_POSITIVE\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SP4_LOCAL_PACKET FAILURE: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
