#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

template <class Function>
void outputs(const int k, const int a, const int b, Function function) {
    const int upper = std::min(a + b, 2 * k - a - b);
    for (int c = std::abs(a - b); c <= upper; c += 2) {
        function(c);
    }
}

int multiplicity(
    const int k,
    const std::array<int, 3>& labels,
    const int target
) {
    std::vector<int> current(static_cast<std::size_t>(k + 1), 0);
    std::vector<int> next(static_cast<std::size_t>(k + 1), 0);
    current[0] = 1;
    for (const int label : labels) {
        std::fill(next.begin(), next.end(), 0);
        for (int source = 0; source <= k; ++source) {
            const int value = current[static_cast<std::size_t>(source)];
            if (value == 0) {
                continue;
            }
            outputs(k, source, label, [&](const int output) {
                next[static_cast<std::size_t>(output)] += value;
            });
        }
        current.swap(next);
    }
    return current[static_cast<std::size_t>(target)];
}

bool fuses_to(const int k, const int a, const int b, const int target) {
    bool found = false;
    outputs(k, a, b, [&](const int output) {
        found = found || output == target;
    });
    return found;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error("usage: audit_su2_h13_offsets MAX_LEVEL");
        }
        const int maximum_level = std::stoi(argv[1]);
        for (int k = 12; k <= maximum_level; k += 2) {
            const int a = k - 4;
            int rows = 0;
            for (int r0 = 2; r0 <= k; r0 += 2) {
                for (int r1 = r0; r1 <= k; r1 += 2) {
                    for (int r2 = r1; r2 <= k; r2 += 2) {
                        const std::array<int, 3> r{r0, r1, r2};
                        if (r0 == a || r1 == a || r2 == a
                            || r0 == k || r1 == k || r2 == k) {
                            continue;
                        }
                        if (multiplicity(k, r, a) != 3
                            || multiplicity(k, r, k) == 0) {
                            continue;
                        }
                        int c = 0;
                        int low = 0;
                        int high = 0;
                        for (int i = 0; i < 3; ++i) {
                            const int j = (i + 1) % 3;
                            const int ell = (i + 2) % 3;
                            c += (r[static_cast<std::size_t>(i)] == 2
                                  || r[static_cast<std::size_t>(i)] == 4)
                                && fuses_to(
                                    k,
                                    r[static_cast<std::size_t>(j)],
                                    r[static_cast<std::size_t>(ell)],
                                    a
                                );
                            low += std::abs(
                                       r[static_cast<std::size_t>(i)] - a
                                   ) == 2
                                && std::abs(
                                       r[static_cast<std::size_t>(j)]
                                       - r[static_cast<std::size_t>(ell)]
                                   ) <= 2;
                            high += std::abs(
                                        r[static_cast<std::size_t>(i)] - 4
                                    ) <= 2
                                && std::abs(
                                       r[static_cast<std::size_t>(j)]
                                       + r[static_cast<std::size_t>(ell)] - k
                                   ) <= 2;
                        }
                        if (c + low + high < 4) {
                            continue;
                        }
                        ++rows;
                        std::cout << "k=" << k
                                  << " R=[" << r0 << ',' << r1 << ',' << r2
                                  << "] C=" << c << " L=" << low
                                  << " H=" << high << '\n';
                    }
                }
            }
            std::cerr << "finished k=" << k << " rows=" << rows << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
