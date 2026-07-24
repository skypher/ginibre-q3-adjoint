#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Count = std::int64_t;

Count invariant_multiplicity(const std::vector<int>& labels) {
    int total = 0;
    for (const int label : labels) {
        total += label;
    }
    std::vector<Count> current(static_cast<std::size_t>(total + 1), 0);
    std::vector<Count> next(static_cast<std::size_t>(total + 1), 0);
    current[0] = 1;
    int maximum = 0;
    for (const int label : labels) {
        std::fill(next.begin(), next.end(), 0);
        for (int source = 0; source <= maximum; ++source) {
            const Count multiplicity =
                current[static_cast<std::size_t>(source)];
            if (multiplicity == 0) {
                continue;
            }
            for (int target = std::abs(source - label);
                 target <= source + label; target += 2) {
                next[static_cast<std::size_t>(target)] += multiplicity;
            }
        }
        maximum += label;
        current.swap(next);
    }
    return current[0];
}

bool disjoint_support(
    const std::array<int, 2>& minus,
    const std::array<int, 5>& plus
) {
    for (const int negative : minus) {
        for (const int positive : plus) {
            if (negative == positive) {
                return false;
            }
        }
    }
    return true;
}

struct Statistics {
    Count invariant = 0;
    Count negative = 0;
    Count equal_reservoir = 0;
    int active = 0;
    Count maximum_cut = 0;
};

Statistics statistics(
    const std::array<int, 2>& minus,
    const std::array<int, 5>& plus
) {
    Statistics result;
    result.invariant = invariant_multiplicity({
        minus[0], minus[1],
        plus[0], plus[1], plus[2], plus[3], plus[4]
    });

    for (int orientation = 0; orientation < 2; ++orientation) {
        for (int first = 0; first < 5; ++first) {
            for (int second = first + 1; second < 5; ++second) {
                const Count triple = invariant_multiplicity({
                    minus[static_cast<std::size_t>(orientation)],
                    plus[static_cast<std::size_t>(first)],
                    plus[static_cast<std::size_t>(second)]
                });
                if (triple == 0) {
                    continue;
                }
                std::vector<int> complement{
                    minus[static_cast<std::size_t>(1 - orientation)]
                };
                for (int index = 0; index < 5; ++index) {
                    if (index != first && index != second) {
                        complement.push_back(
                            plus[static_cast<std::size_t>(index)]
                        );
                    }
                }
                const Count cut = triple
                    * invariant_multiplicity(complement);
                if (cut != 0) {
                    ++result.active;
                    result.negative += cut;
                    result.maximum_cut =
                        std::max(result.maximum_cut, cut);
                }
            }
        }
    }

    if (minus[0] == minus[1]) {
        result.equal_reservoir += invariant_multiplicity({
            plus[0], plus[1], plus[2], plus[3], plus[4]
        });
    }
    for (int first = 0; first < 5; ++first) {
        for (int second = first + 1; second < 5; ++second) {
            if (plus[static_cast<std::size_t>(first)]
                != plus[static_cast<std::size_t>(second)]) {
                continue;
            }
            std::vector<int> remaining{minus[0], minus[1]};
            for (int index = 0; index < 5; ++index) {
                if (index != first && index != second) {
                    remaining.push_back(
                        plus[static_cast<std::size_t>(index)]
                    );
                }
            }
            result.equal_reservoir += invariant_multiplicity(remaining);
        }
    }
    return result;
}

void print_case(
    const std::array<int, 2>& minus,
    const std::array<int, 5>& plus,
    const Statistics& value
) {
    std::cout << " minus=[" << minus[0] << ',' << minus[1] << "] plus=[";
    for (std::size_t index = 0; index < plus.size(); ++index) {
        if (index != 0) {
            std::cout << ',';
        }
        std::cout << plus[index];
    }
    std::cout << "] N=" << value.invariant
              << " T=" << value.negative
              << " P_eq=" << value.equal_reservoir
              << " raw_margin=" << value.invariant - value.negative
              << " paid_margin="
              << value.invariant + value.equal_reservoir - value.negative
              << " active=" << value.active
              << " max_d=" << value.maximum_cut;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: search_su2_all_even_two_minus_dense "
                "MAXIMUM_EVEN_LABEL"
            );
        }
        const long parsed = std::strtol(argv[1], nullptr, 10);
        if (parsed < 2 || parsed > std::numeric_limits<int>::max()
            || parsed % 2 != 0) {
            throw std::runtime_error(
                "MAXIMUM_EVEN_LABEL must be a positive even integer"
            );
        }
        const int maximum = static_cast<int>(parsed);
        std::atomic<std::uint64_t> cases{0};
        std::atomic<std::uint64_t> dense_cases{0};
        std::mutex result_mutex;
        bool found = false;
        std::array<int, 2> found_minus{};
        std::array<int, 5> found_plus{};
        Statistics found_statistics;
        Count minimum_raw_margin = std::numeric_limits<Count>::max();
        std::array<int, 2> minimum_minus{};
        std::array<int, 5> minimum_plus{};
        Statistics minimum_statistics;
        std::array<Count, 4> minimum_twenty_cut_margin{
            std::numeric_limits<Count>::max(),
            std::numeric_limits<Count>::max(),
            std::numeric_limits<Count>::max(),
            std::numeric_limits<Count>::max()
        };
        std::array<std::array<int, 2>, 4> endpoint_minus{};
        std::array<std::array<int, 5>, 4> endpoint_plus{};
        std::array<Statistics, 4> endpoint_statistics{};

#pragma omp parallel for schedule(dynamic)
        for (int first_minus = 2;
             first_minus <= maximum; first_minus += 2) {
            Count local_minimum = std::numeric_limits<Count>::max();
            std::array<int, 2> local_minus{};
            std::array<int, 5> local_plus{};
            Statistics local_statistics;
            std::array<Count, 4> local_endpoint_margin{
                std::numeric_limits<Count>::max(),
                std::numeric_limits<Count>::max(),
                std::numeric_limits<Count>::max(),
                std::numeric_limits<Count>::max()
            };
            std::array<std::array<int, 2>, 4> local_endpoint_minus{};
            std::array<std::array<int, 5>, 4> local_endpoint_plus{};
            std::array<Statistics, 4> local_endpoint_statistics{};
            for (int second_minus = first_minus;
                 second_minus <= maximum; second_minus += 2) {
                const std::array<int, 2> minus{
                    first_minus, second_minus
                };
                for (int a = 2; a <= maximum; a += 2) {
                    for (int b = a; b <= maximum; b += 2) {
                        for (int c = b; c <= maximum; c += 2) {
                            for (int d = c; d <= maximum; d += 2) {
                                for (int e = d; e <= maximum; e += 2) {
                                    const std::array<int, 5> plus{
                                        a, b, c, d, e
                                    };
                                    if (!disjoint_support(minus, plus)) {
                                        continue;
                                    }
                                    const Statistics value =
                                        statistics(minus, plus);
                                    if (value.invariant == 0) {
                                        continue;
                                    }
                                    ++cases;
                                    if (value.active < 13) {
                                        continue;
                                    }
                                    ++dense_cases;
                                    const Count raw_margin =
                                        value.invariant - value.negative;
                                    if (value.maximum_cut <= 3) {
                                        const std::size_t endpoint =
                                            static_cast<std::size_t>(
                                                value.maximum_cut
                                            );
                                        const Count twenty_cut_margin =
                                            value.invariant
                                            - 20 * value.maximum_cut;
                                        if (twenty_cut_margin
                                            < local_endpoint_margin[
                                                endpoint]) {
                                            local_endpoint_margin[endpoint] =
                                                twenty_cut_margin;
                                            local_endpoint_minus[endpoint] =
                                                minus;
                                            local_endpoint_plus[endpoint] =
                                                plus;
                                            local_endpoint_statistics[
                                                endpoint] = value;
                                        }
                                    }
                                    if (raw_margin < local_minimum) {
                                        local_minimum = raw_margin;
                                        local_minus = minus;
                                        local_plus = plus;
                                        local_statistics = value;
                                    }
                                    if (raw_margin < 0) {
                                        std::lock_guard<std::mutex> lock(
                                            result_mutex
                                        );
                                        if (!found) {
                                            found = true;
                                            found_minus = minus;
                                            found_plus = plus;
                                            found_statistics = value;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (local_minimum != std::numeric_limits<Count>::max()) {
                std::lock_guard<std::mutex> lock(result_mutex);
                if (local_minimum < minimum_raw_margin) {
                    minimum_raw_margin = local_minimum;
                    minimum_minus = local_minus;
                    minimum_plus = local_plus;
                    minimum_statistics = local_statistics;
                }
                for (std::size_t endpoint = 1;
                     endpoint <= 3; ++endpoint) {
                    if (local_endpoint_margin[endpoint]
                        < minimum_twenty_cut_margin[endpoint]) {
                        minimum_twenty_cut_margin[endpoint] =
                            local_endpoint_margin[endpoint];
                        endpoint_minus[endpoint] =
                            local_endpoint_minus[endpoint];
                        endpoint_plus[endpoint] =
                            local_endpoint_plus[endpoint];
                        endpoint_statistics[endpoint] =
                            local_endpoint_statistics[endpoint];
                    }
                }
            }
        }

        std::cout << "SU2_ALL_EVEN_TWO_MINUS_DENSE maximum="
                  << maximum
                  << " cases=" << cases.load()
                  << " dense_cases=" << dense_cases.load();
#ifdef _OPENMP
        std::cout << " threads=" << omp_get_max_threads();
#else
        std::cout << " threads=1";
#endif
        if (found) {
            std::cout << " result=RAW_FAIL first";
            print_case(found_minus, found_plus, found_statistics);
        } else {
            std::cout << " result=RAW_PASS minimum";
            if (minimum_raw_margin == std::numeric_limits<Count>::max()) {
                std::cout << " none";
            } else {
                print_case(
                    minimum_minus, minimum_plus, minimum_statistics
                );
            }
        }
        std::cout << '\n';
        for (std::size_t endpoint = 1; endpoint <= 3; ++endpoint) {
            std::cout << "MAX_D_ENDPOINT d=" << endpoint << " minimum";
            if (minimum_twenty_cut_margin[endpoint]
                == std::numeric_limits<Count>::max()) {
                std::cout << " none";
            } else {
                print_case(
                    endpoint_minus[endpoint],
                    endpoint_plus[endpoint],
                    endpoint_statistics[endpoint]
                );
                std::cout << " N_minus_20d="
                          << minimum_twenty_cut_margin[endpoint];
            }
            std::cout << '\n';
        }
        return found ? EXIT_FAILURE : EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
