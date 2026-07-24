#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

namespace {

using Labels = std::array<int, 6>;

struct Case {
    int level;
    Labels labels;
};

struct Failure {
    bool present = false;
    std::string reason;
    int level = 0;
    Labels labels{};
    unsigned int minus_mask = 0;
};

struct Statistics {
    std::uint64_t sextuples = 0;
    std::uint64_t active_sextuples = 0;
    std::uint64_t signed_cases = 0;
    std::int64_t minimum_margin =
        std::numeric_limits<std::int64_t>::max();
    std::int64_t minimum_active_margin =
        std::numeric_limits<std::int64_t>::max();
    std::int64_t minimum_signed_half =
        std::numeric_limits<std::int64_t>::max();
    std::int64_t minimum_active_signed_half =
        std::numeric_limits<std::int64_t>::max();
};

template <class Function>
void fusion_outputs(
    const int level,
    const int first,
    const int second,
    Function function
) {
    const int upper = std::min(
        first + second,
        2 * level - first - second
    );
    for (int output = std::abs(first - second);
         output <= upper;
         output += 2) {
        function(output);
    }
}

std::int64_t invariant(
    const int level,
    const Labels& labels,
    const unsigned int mask
) {
    std::vector<std::int64_t> current(
        static_cast<std::size_t>(level + 1),
        0
    );
    current[0] = 1;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        if ((mask & (1U << index)) == 0U) {
            continue;
        }
        std::vector<std::int64_t> next(
            static_cast<std::size_t>(level + 1),
            0
        );
        for (int source = 0; source <= level; ++source) {
            const std::int64_t multiplicity =
                current[static_cast<std::size_t>(source)];
            if (multiplicity == 0) {
                continue;
            }
            fusion_outputs(
                level,
                source,
                labels[index],
                [&](const int target) {
                    next[static_cast<std::size_t>(target)]
                        += multiplicity;
                }
            );
        }
        current = std::move(next);
    }
    return current[0];
}

bool support_disjoint(
    const Labels& labels,
    const unsigned int minus_mask
) {
    for (std::size_t first = 0; first < labels.size(); ++first) {
        for (std::size_t second = first + 1;
             second < labels.size();
             ++second) {
            if (labels[first] == labels[second]
                && ((minus_mask >> first) & 1U)
                    != ((minus_mask >> second) & 1U)) {
                return false;
            }
        }
    }
    return true;
}

std::uint64_t binomial(std::uint64_t n, std::uint64_t r) {
    if (r > n) {
        return 0;
    }
    r = std::min(r, n - r);
    std::uint64_t result = 1;
    for (std::uint64_t index = 1; index <= r; ++index) {
        result = result * (n - r + index) / index;
    }
    return result;
}

int parse_maximum_level(const int argc, char** argv) {
    if (argc > 2) {
        throw std::runtime_error(
            "usage: verify_su2_finite_six_pair_reservoir [MAX_K]"
        );
    }
    if (argc == 1) {
        return 14;
    }
    std::size_t parsed = 0;
    const std::string argument = argv[1];
    const int value = std::stoi(argument, &parsed);
    if (parsed != argument.size() || value < 1) {
        throw std::runtime_error("MAX_K must be a positive integer");
    }
    return value;
}

std::vector<Case> enumerate_cases(const int maximum_level) {
    std::vector<Case> cases;
    for (int level = 1; level <= maximum_level; ++level) {
        Labels labels{};
        const auto visit = [&](const auto& self,
                               const int index,
                               const int lower) -> void {
            if (index == 6) {
                cases.push_back({level, labels});
                return;
            }
            for (int label = lower; label <= level; ++label) {
                labels[static_cast<std::size_t>(index)] = label;
                self(self, index + 1, label);
            }
        };
        visit(visit, 0, 1);
    }
    return cases;
}

void merge_statistics(Statistics& target, const Statistics& source) {
    target.sextuples += source.sextuples;
    target.active_sextuples += source.active_sextuples;
    target.signed_cases += source.signed_cases;
    target.minimum_margin =
        std::min(target.minimum_margin, source.minimum_margin);
    target.minimum_active_margin = std::min(
        target.minimum_active_margin,
        source.minimum_active_margin
    );
    target.minimum_signed_half = std::min(
        target.minimum_signed_half,
        source.minimum_signed_half
    );
    target.minimum_active_signed_half = std::min(
        target.minimum_active_signed_half,
        source.minimum_active_signed_half
    );
}

void record_failure(
    Failure& failure,
    const std::string& reason,
    const Case& current_case,
    const unsigned int minus_mask = 0
) {
    if (!failure.present) {
        failure.present = true;
        failure.reason = reason;
        failure.level = current_case.level;
        failure.labels = current_case.labels;
        failure.minus_mask = minus_mask;
    }
}

void print_labels(const Labels& labels) {
    std::cerr << '[';
    for (std::size_t index = 0; index < labels.size(); ++index) {
        if (index != 0) {
            std::cerr << ',';
        }
        std::cerr << labels[index];
    }
    std::cerr << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = parse_maximum_level(argc, argv);
        const std::vector<Case> cases = enumerate_cases(maximum_level);
        Statistics statistics;
        Failure failure;
        constexpr unsigned int full_mask = (1U << 6U) - 1U;

#pragma omp parallel
        {
            Statistics local;
            Failure local_failure;

#pragma omp for schedule(dynamic, 32)
            for (std::size_t case_index = 0;
                 case_index < cases.size();
                 ++case_index) {
                const Case& current_case = cases[case_index];
                std::array<std::int64_t, 64> invariants{};
                for (unsigned int mask = 0; mask <= full_mask; ++mask) {
                    invariants[mask] = invariant(
                        current_case.level,
                        current_case.labels,
                        mask
                    );
                }

                const std::int64_t sixfold = invariants[full_mask];
                std::int64_t pair_reservoir = 0;
                for (std::size_t first = 0;
                     first < current_case.labels.size();
                     ++first) {
                    for (std::size_t second = first + 1;
                         second < current_case.labels.size();
                         ++second) {
                        if (current_case.labels[first]
                            != current_case.labels[second]) {
                            continue;
                        }
                        const unsigned int pair_mask =
                            (1U << first) | (1U << second);
                        pair_reservoir +=
                            invariants[full_mask ^ pair_mask];
                    }
                }

                int split_count = 0;
                for (unsigned int mask = 0; mask <= full_mask; ++mask) {
                    if (std::popcount(mask) != 3
                        || mask > (full_mask ^ mask)) {
                        continue;
                    }
                    if (invariants[mask] != 0
                        && invariants[full_mask ^ mask] != 0) {
                        ++split_count;
                    }
                }

                const std::int64_t margin =
                    sixfold + pair_reservoir - split_count;
                ++local.sextuples;
                local.minimum_margin =
                    std::min(local.minimum_margin, margin);
                if (split_count != 0) {
                    ++local.active_sextuples;
                    local.minimum_active_margin = std::min(
                        local.minimum_active_margin,
                        margin
                    );
                }
                if (margin < 0) {
                    record_failure(
                        local_failure,
                        "negative pair-reservoir margin",
                        current_case
                    );
                }

                for (unsigned int minus_mask = 0;
                     minus_mask <= full_mask;
                     ++minus_mask) {
                    if ((std::popcount(minus_mask) & 1) != 0
                        || !support_disjoint(
                            current_case.labels,
                            minus_mask
                        )) {
                        continue;
                    }
                    ++local.signed_cases;
                    std::int64_t signed_contraction = 0;
                    for (unsigned int mask = 0;
                         mask <= full_mask;
                         ++mask) {
                        const std::int64_t term =
                            invariants[mask]
                            * invariants[full_mask ^ mask];
                        signed_contraction +=
                            (std::popcount(mask & minus_mask) & 1) == 0
                            ? term : -term;
                    }

                    int positive_middle = 0;
                    int negative_middle = 0;
                    for (unsigned int mask = 0;
                         mask <= full_mask;
                         ++mask) {
                        if (std::popcount(mask) != 3
                            || mask > (full_mask ^ mask)
                            || invariants[mask] == 0
                            || invariants[full_mask ^ mask] == 0) {
                            continue;
                        }
                        if ((std::popcount(mask & minus_mask) & 1) == 0) {
                            ++positive_middle;
                        } else {
                            ++negative_middle;
                        }
                    }
                    const std::int64_t signed_half =
                        sixfold + pair_reservoir
                        + positive_middle - negative_middle;
                    local.minimum_signed_half = std::min(
                        local.minimum_signed_half,
                        signed_half
                    );
                    if (negative_middle != 0) {
                        local.minimum_active_signed_half = std::min(
                            local.minimum_active_signed_half,
                            signed_half
                        );
                    }
                    if (signed_half < 0
                        || signed_contraction != 2 * signed_half) {
                        record_failure(
                            local_failure,
                            signed_half < 0
                                ? "negative signed half-contraction"
                                : "subset identity mismatch",
                            current_case,
                            minus_mask
                        );
                    }
                }
            }

#pragma omp critical
            {
                merge_statistics(statistics, local);
                if (local_failure.present && !failure.present) {
                    failure = local_failure;
                }
            }
        }

        const std::uint64_t expected_cases = binomial(
            static_cast<std::uint64_t>(maximum_level + 6),
            7
        );
        if (statistics.sextuples != expected_cases) {
            throw std::runtime_error(
                "sorted-sextuple enumeration count mismatch"
            );
        }
        if (failure.present) {
            std::cerr << "SU2_FINITE_SIX_PAIR_RESERVOIR FAIL"
                      << " reason=" << failure.reason
                      << " k=" << failure.level
                      << " labels=";
            print_labels(failure.labels);
            std::cerr << " minus_mask=" << failure.minus_mask << '\n';
            return EXIT_FAILURE;
        }

        const auto display_minimum = [](
            const std::int64_t value
        ) -> std::int64_t {
            return value == std::numeric_limits<std::int64_t>::max()
                ? 0 : value;
        };
        std::cout
            << "SU2_FINITE_SIX_PAIR_RESERVOIR PASS"
            << " max_k=" << maximum_level
            << " threads=" << omp_get_max_threads()
            << " sextuples=" << statistics.sextuples
            << " active_sextuples=" << statistics.active_sextuples
            << " signed_support_disjoint_cases="
            << statistics.signed_cases
            << " minimum_margin="
            << display_minimum(statistics.minimum_margin)
            << " minimum_active_margin="
            << display_minimum(statistics.minimum_active_margin)
            << " minimum_signed_half="
            << display_minimum(statistics.minimum_signed_half)
            << " minimum_active_signed_half="
            << display_minimum(statistics.minimum_active_signed_half)
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_FINITE_SIX_PAIR_RESERVOIR ERROR "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
