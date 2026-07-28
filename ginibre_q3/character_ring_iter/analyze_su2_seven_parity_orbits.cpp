#include <atomic>
#include <cmath>
#include <fstream>
#include <mutex>
#include <thread>

#define main analyze_su2_four_minus_shallow_main
#include "analyze_su2_four_minus_shallow.cpp"
#undef main

namespace {

struct Orbit {
    int minus_count;
    int odd_count;
    int odd_minus_count;
};

struct OrbitStratum : Stratum {
    Witness local_only;
    Witness local_positive;
    Witness maximal_ceiling;
    Witness all_cut_local_ceiling;
    Witness all_cut_pair_ceiling;
    Witness all_cut_positive_ceiling;
    Witness all_cut_ceiling;
    Witness zero_pair_ceiling;
    Witness low_band_bound;
    Witness direct;
};

constexpr std::array<Orbit, 8> orbits{{
    {4, 0, 0},
    {4, 2, 0},
    {4, 2, 1},
    {4, 4, 1},
    {4, 4, 2},
    {4, 6, 3},
    {6, 0, 0},
    {6, 2, 1},
}};

Labels merge_labels(const Labels& first, const Labels& second) {
    Labels result;
    result.reserve(first.size() + second.size());
    std::merge(
        first.begin(), first.end(),
        second.begin(), second.end(),
        std::back_inserter(result)
    );
    return result;
}

std::vector<int> probe_outputs(
    const Labels& first,
    int level
) {
    int common_parity = 0;
    for (const int label : first) {
        common_parity ^= label & 1;
    }
    std::vector<int> result;
    for (int j = 0; j <= 4; ++j) {
        const int low = common_parity + 2 * j;
        if (low <= level) {
            result.push_back(low);
        }
    }
    const int reflected_parity =
        (level + common_parity) & 1;
    for (int j = 0; j <= 4; ++j) {
        const int reflected_low = reflected_parity + 2 * j;
        const int high = level - reflected_low;
        if (high < 0
            || std::find(result.begin(), result.end(), high)
                != result.end()) {
            continue;
        }
        result.push_back(high);
    }
    return result;
}

std::int64_t local_channels(
    const Labels& first,
    const Labels& second,
    int level
) {
    std::int64_t result = 0;
    for (const int output : probe_outputs(first, level)) {
        result += multiplicity(first, output, level)
            * multiplicity(second, output, level);
    }
    return result;
}

std::int64_t low_band_channels(
    const Labels& first,
    const Labels& second,
    int level
) {
    std::int64_t result = 0;
    constexpr std::array<int, 4> outputs{0, 2, 4, 6};
    for (const int output : outputs) {
        result += multiplicity(first, output, level)
            * multiplicity(second, output, level);
    }
    return result;
}

void merge_witness(Witness& target, const Witness& source) {
    if (source.initialized
        && (!target.initialized || source.value < target.value)) {
        target = source;
    }
}

[[maybe_unused]] void merge_stratum(
    OrbitStratum& target,
    const OrbitStratum& source
) {
    target.cases += source.cases;
    target.active_cases += source.active_cases;
    merge_witness(target.raw, source.raw);
    merge_witness(target.pair, source.pair);
    merge_witness(target.active_pair, source.active_pair);
    merge_witness(target.local_pair, source.local_pair);
    merge_witness(target.local_pair_low, source.local_pair_low);
    merge_witness(target.local_only, source.local_only);
    merge_witness(target.local_positive, source.local_positive);
    merge_witness(target.maximal_ceiling, source.maximal_ceiling);
    merge_witness(
        target.all_cut_local_ceiling,
        source.all_cut_local_ceiling
    );
    merge_witness(
        target.all_cut_pair_ceiling,
        source.all_cut_pair_ceiling
    );
    merge_witness(
        target.all_cut_positive_ceiling,
        source.all_cut_positive_ceiling
    );
    merge_witness(target.all_cut_ceiling, source.all_cut_ceiling);
    merge_witness(target.zero_pair_ceiling, source.zero_pair_ceiling);
    merge_witness(target.low_band_bound, source.low_band_bound);
    merge_witness(target.direct, source.direct);
    merge_witness(target.exact, source.exact);
    target.maximum_ratio = std::max(
        target.maximum_ratio, source.maximum_ratio
    );
}

[[maybe_unused]] std::array<OrbitStratum, 2> scan_level(
    const Orbit& orbit,
    int level
) {
    Labels even;
    Labels odd;
    for (int label = 1; label <= level; ++label) {
        ((label & 1) == 0 ? even : odd).push_back(label);
    }

    const int even_minus =
        orbit.minus_count - orbit.odd_minus_count;
    const int odd_plus =
        orbit.odd_count - orbit.odd_minus_count;
    const int plus_count = 7 - orbit.minus_count;
    const int even_plus = plus_count - odd_plus;
    const std::vector<Labels> minus_even =
        sorted_lists(even, even_minus);
    const std::vector<Labels> minus_odd =
        sorted_lists(odd, orbit.odd_minus_count);
    const std::vector<Labels> plus_even =
        sorted_lists(even, even_plus);
    const std::vector<Labels> plus_odd =
        sorted_lists(odd, odd_plus);

    std::vector<Labels> minus_lists;
    std::vector<Labels> plus_lists;
    minus_lists.reserve(minus_even.size() * minus_odd.size());
    plus_lists.reserve(plus_even.size() * plus_odd.size());
    for (const Labels& first : minus_even) {
        for (const Labels& second : minus_odd) {
            minus_lists.push_back(merge_labels(first, second));
        }
    }
    for (const Labels& first : plus_even) {
        for (const Labels& second : plus_odd) {
            plus_lists.push_back(merge_labels(first, second));
        }
    }

    std::array<OrbitStratum, 2> strata{};
    constexpr unsigned int full_mask = (1U << 7U) - 1U;
    for (const Labels& minus : minus_lists) {
        for (const Labels& plus : plus_lists) {
            if (!disjoint_support(minus, plus)) {
                continue;
            }
            Labels labels = minus;
            labels.insert(labels.end(), plus.begin(), plus.end());

            const std::int64_t sevenfold =
                invariant(labels, level);
            std::int64_t pair = 0;
            std::int64_t pair_low = 0;
            for (std::size_t i = 0U; i < labels.size(); ++i) {
                for (std::size_t j = i + 1U;
                     j < labels.size(); ++j) {
                    if (labels[i] != labels[j]) {
                        continue;
                    }
                    Labels rest;
                    for (std::size_t t = 0U;
                         t < labels.size(); ++t) {
                        if (t != i && t != j) {
                            rest.push_back(labels[t]);
                        }
                    }
                    pair += invariant(rest, level);
                    const Labels first{rest[0], rest[1]};
                    const Labels second{
                        rest[2], rest[3], rest[4]
                    };
                    pair_low += local_channels(
                        first, second, level
                    );
                }
            }

            std::int64_t positive = 0;
            std::int64_t negative = 0;
            std::int64_t maximum_cut = 0;
            int negative_cut_count = 0;
            std::vector<unsigned int> maximum_masks;
            std::vector<std::pair<unsigned int, std::int64_t>>
                active_negative_cuts;
            for (unsigned int mask = 0U;
                 mask <= full_mask; ++mask) {
                if (popcount(mask) != 3) {
                    continue;
                }
                const std::int64_t term =
                    invariant(subset(labels, mask, true), level)
                    * invariant(
                        subset(labels, mask, false), level
                    );
                int minus_in_cut = 0;
                for (int i = 0; i < orbit.minus_count; ++i) {
                    minus_in_cut += static_cast<int>(
                        (mask >> i) & 1U
                    );
                }
                if ((minus_in_cut & 1) == 0) {
                    positive += term;
                } else {
                    ++negative_cut_count;
                    negative += term;
                    if (term > 0) {
                        active_negative_cuts.emplace_back(mask, term);
                    }
                    if (term > maximum_cut) {
                        maximum_cut = term;
                        maximum_masks.clear();
                        maximum_masks.push_back(mask);
                    } else if (
                        term != 0 && term == maximum_cut
                    ) {
                        maximum_masks.push_back(mask);
                    }
                }
            }

            std::int64_t best_local = 0;
            for (const unsigned int mask : maximum_masks) {
                best_local = std::max(
                    best_local,
                    local_channels(
                        subset(labels, mask, true),
                        subset(labels, mask, false),
                        level
                    )
                );
            }
            std::int64_t all_cut_ceiling =
                std::numeric_limits<std::int64_t>::max();
            std::int64_t all_cut_local_ceiling = all_cut_ceiling;
            std::int64_t all_cut_pair_ceiling = all_cut_ceiling;
            std::int64_t all_cut_positive_ceiling = all_cut_ceiling;
            std::int64_t low_band_bound = all_cut_ceiling;
            for (const auto& [mask, term] : active_negative_cuts) {
                const Labels cut_labels =
                    subset(labels, mask, true);
                const Labels rest_labels =
                    subset(labels, mask, false);
                const std::int64_t cut_local =
                    local_channels(cut_labels, rest_labels, level);
                const std::int64_t cut_low =
                    low_band_channels(
                        cut_labels, rest_labels, level
                    );
                const std::int64_t ceiling =
                    negative_cut_count * term;
                low_band_bound = std::min(
                    low_band_bound,
                    cut_low - (24 * term - 43)
                );
                all_cut_local_ceiling = std::min(
                    all_cut_local_ceiling,
                    cut_local - ceiling
                );
                all_cut_pair_ceiling = std::min(
                    all_cut_pair_ceiling,
                    cut_local + pair_low - ceiling
                );
                all_cut_positive_ceiling = std::min(
                    all_cut_positive_ceiling,
                    cut_local + positive - ceiling
                );
                all_cut_ceiling = std::min(
                    all_cut_ceiling,
                    cut_local + pair_low + positive - ceiling
                );
            }

            bool shallow = false;
            for (const int label : labels) {
                shallow = shallow
                    || std::min(label, level - label) <= 1;
            }
            OrbitStratum& stratum =
                strata[static_cast<std::size_t>(shallow)];
            ++stratum.cases;
            consider(
                stratum.raw, sevenfold - negative, level,
                minus, plus, sevenfold, pair, positive,
                negative, maximum_cut
            );
            consider(
                stratum.pair, sevenfold + pair - negative, level,
                minus, plus, sevenfold, pair, positive,
                negative, maximum_cut
            );
            if (negative > 0) {
                ++stratum.active_cases;
                consider(
                    stratum.active_pair,
                    sevenfold + pair - negative, level,
                    minus, plus, sevenfold, pair, positive,
                    negative, maximum_cut
                );
                consider(
                    stratum.local_pair,
                    best_local + pair - negative, level,
                    minus, plus, sevenfold, pair, positive,
                    negative, maximum_cut
                );
                consider(
                    stratum.local_pair_low,
                    best_local + pair_low - negative, level,
                    minus, plus, sevenfold, pair_low, positive,
                    negative, maximum_cut
                );
                consider(
                    stratum.local_only,
                    best_local - negative, level,
                    minus, plus, sevenfold, pair_low, positive,
                    negative, maximum_cut
                );
                consider(
                    stratum.local_positive,
                    best_local + positive - negative, level,
                    minus, plus, sevenfold, pair_low, positive,
                    negative, maximum_cut
                );
                consider(
                    stratum.maximal_ceiling,
                    best_local + pair_low + positive
                        - negative_cut_count * maximum_cut,
                    level, minus, plus, sevenfold, pair_low,
                    positive, negative, maximum_cut
                );
                consider(
                    stratum.all_cut_ceiling,
                    all_cut_ceiling, level,
                    minus, plus, sevenfold, pair_low,
                    positive, negative, maximum_cut
                );
                consider(
                    stratum.all_cut_local_ceiling,
                    all_cut_local_ceiling, level,
                    minus, plus, sevenfold, pair_low,
                    positive, negative, maximum_cut
                );
                consider(
                    stratum.all_cut_pair_ceiling,
                    all_cut_pair_ceiling, level,
                    minus, plus, sevenfold, pair_low,
                    positive, negative, maximum_cut
                );
                consider(
                    stratum.all_cut_positive_ceiling,
                    all_cut_positive_ceiling, level,
                    minus, plus, sevenfold, pair_low,
                    positive, negative, maximum_cut
                );
                if (pair_low == 0) {
                    consider(
                        stratum.zero_pair_ceiling,
                        all_cut_local_ceiling, level,
                        minus, plus, sevenfold, pair_low,
                        positive, negative, maximum_cut
                    );
                }
                consider(
                    stratum.low_band_bound,
                    low_band_bound, level,
                    minus, plus, sevenfold, pair_low,
                    positive, negative, maximum_cut
                );
                consider(
                    stratum.direct,
                    best_local + pair_low + positive - negative,
                    level, minus, plus, sevenfold, pair_low,
                    positive, negative, maximum_cut
                );
                const int ratio = static_cast<int>(
                    (negative + maximum_cut - 1) / maximum_cut
                );
                stratum.maximum_ratio = std::max(
                    stratum.maximum_ratio, ratio
                );
            }
            consider(
                stratum.exact,
                sevenfold + pair + positive - negative, level,
                minus, plus, sevenfold, pair, positive,
                negative, maximum_cut
            );
        }
    }
    return strata;
}

[[maybe_unused]] unsigned worker_limit_orbits() {
    const unsigned hardware = std::max(
        1U, std::thread::hardware_concurrency()
    );
    double load_average = 0.0;
    const int load_status = getloadavg(&load_average, 1);
    const unsigned loaded = load_status == 1
        ? static_cast<unsigned>(std::max(0.0, std::ceil(load_average)))
        : 0U;
    const unsigned cpu_limit =
        hardware > loaded ? hardware - loaded : 1U;

    std::ifstream memory("/proc/meminfo");
    std::string key;
    std::uint64_t value = 0;
    std::string unit;
    std::uint64_t available_kib = 0;
    while (memory >> key >> value >> unit) {
        if (key == "MemAvailable:") {
            available_kib = value;
            break;
        }
    }
    constexpr std::uint64_t kib_per_worker = 256U * 1024U;
    const unsigned memory_limit = available_kib == 0
        ? hardware
        : static_cast<unsigned>(std::max<std::uint64_t>(
            1U, available_kib / kib_per_worker
        ));
    return std::max(1U, std::min(cpu_limit, memory_limit));
}

}  // namespace

#ifndef SU2_SEVEN_PARITY_ORBITS_NO_MAIN
int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_seven_parity_orbits maximum_level"
            );
        }
        const int maximum_level = parse_positive(argv[1]);
        if (maximum_level < 4) {
            throw std::runtime_error(
                "maximum level must be at least four"
            );
        }

        constexpr int first_level = 4;
        const int levels = maximum_level - first_level + 1;
        const int task_count =
            static_cast<int>(orbits.size()) * levels;
        std::array<
            std::array<OrbitStratum, 2>,
            orbits.size()
        > totals{};
        std::array<std::mutex, orbits.size()> locks;
        std::atomic<int> next_task{0};
        const unsigned workers = std::min<unsigned>(
            worker_limit_orbits(),
            static_cast<unsigned>(task_count)
        );
        std::vector<std::jthread> pool;
        pool.reserve(workers);
        for (unsigned worker = 0; worker < workers; ++worker) {
            pool.emplace_back([&]() {
                while (true) {
                    const int task = next_task.fetch_add(1);
                    if (task >= task_count) {
                        return;
                    }
                    const std::size_t orbit_index =
                        static_cast<std::size_t>(task / levels);
                    const int level = first_level + task % levels;
                    const auto local = scan_level(
                        orbits[orbit_index], level
                    );
                    std::lock_guard<std::mutex> lock(
                        locks[orbit_index]
                    );
                    for (std::size_t depth = 0U;
                         depth < 2U; ++depth) {
                        merge_stratum(
                            totals[orbit_index][depth],
                            local[depth]
                        );
                    }
                }
            });
        }
        pool.clear();

        std::cout << "SU2_SEVEN_PARITY_ORBITS maximum_level="
                  << maximum_level
                  << " workers=" << workers << '\n';
        OrbitStratum aggregate{};
        std::array<OrbitStratum, 2> depth_aggregates{};
        for (std::size_t index = 0U;
             index < orbits.size(); ++index) {
            const Orbit& orbit = orbits[index];
            for (std::size_t depth = 0U; depth < 2U; ++depth) {
                const OrbitStratum& stratum =
                    totals[index][depth];
                if (stratum.cases == 0U) {
                    continue;
                }
                merge_stratum(aggregate, stratum);
                merge_stratum(depth_aggregates[depth], stratum);
                std::cout
                    << "m=" << orbit.minus_count
                    << " o=" << orbit.odd_count
                    << " r=" << orbit.odd_minus_count
                    << " depth=" << (depth == 0U ? "deep" : "shallow")
                    << " cases=" << stratum.cases
                    << " active=" << stratum.active_cases;
                print_witness("raw", stratum.raw);
                print_witness("pair", stratum.pair);
                print_witness("active_pair", stratum.active_pair);
                print_witness("local_pair", stratum.local_pair);
                print_witness(
                    "local_pair_truncated",
                    stratum.local_pair_low
                );
                print_witness("local_only", stratum.local_only);
                print_witness(
                    "local_positive", stratum.local_positive
                );
                print_witness(
                    "maximal_ceiling", stratum.maximal_ceiling
                );
                print_witness(
                    "all_cut_ceiling", stratum.all_cut_ceiling
                );
                print_witness(
                    "all_cut_local_ceiling",
                    stratum.all_cut_local_ceiling
                );
                print_witness(
                    "all_cut_pair_ceiling",
                    stratum.all_cut_pair_ceiling
                );
                print_witness(
                    "all_cut_positive_ceiling",
                    stratum.all_cut_positive_ceiling
                );
                print_witness(
                    "zero_pair_ceiling",
                    stratum.zero_pair_ceiling
                );
                print_witness(
                    "low_band_bound", stratum.low_band_bound
                );
                print_witness("direct", stratum.direct);
                print_witness("exact", stratum.exact);
                std::cout << " max_ceil_T_over_d="
                          << stratum.maximum_ratio << '\n';
            }
        }
        std::cout
            << "ABLATION_SUMMARY cases=" << aggregate.cases
            << " active=" << aggregate.active_cases;
        print_witness("local_only", aggregate.local_only);
        print_witness(
            "local_pair_truncated", aggregate.local_pair_low
        );
        print_witness(
            "local_positive", aggregate.local_positive
        );
        print_witness(
            "maximal_ceiling", aggregate.maximal_ceiling
        );
        print_witness(
            "all_cut_ceiling", aggregate.all_cut_ceiling
        );
        print_witness(
            "all_cut_local_ceiling",
            aggregate.all_cut_local_ceiling
        );
        print_witness(
            "all_cut_pair_ceiling",
            aggregate.all_cut_pair_ceiling
        );
        print_witness(
            "all_cut_positive_ceiling",
            aggregate.all_cut_positive_ceiling
        );
        print_witness(
            "zero_pair_ceiling", aggregate.zero_pair_ceiling
        );
        print_witness(
            "low_band_bound", aggregate.low_band_bound
        );
        print_witness("direct", aggregate.direct);
        std::cout << '\n';
        for (std::size_t depth = 0U; depth < 2U; ++depth) {
            const OrbitStratum& depth_aggregate =
                depth_aggregates[depth];
            std::cout
                << "ABLATION_DEPTH depth="
                << (depth == 0U ? "deep" : "shallow")
                << " cases=" << depth_aggregate.cases
                << " active=" << depth_aggregate.active_cases;
            print_witness(
                "local_only", depth_aggregate.local_only
            );
            print_witness(
                "local_pair_truncated",
                depth_aggregate.local_pair_low
            );
            print_witness(
                "local_positive",
                depth_aggregate.local_positive
            );
            print_witness(
                "maximal_ceiling",
                depth_aggregate.maximal_ceiling
            );
            print_witness(
                "all_cut_ceiling",
                depth_aggregate.all_cut_ceiling
            );
            print_witness(
                "all_cut_local_ceiling",
                depth_aggregate.all_cut_local_ceiling
            );
            print_witness(
                "all_cut_pair_ceiling",
                depth_aggregate.all_cut_pair_ceiling
            );
            print_witness(
                "all_cut_positive_ceiling",
                depth_aggregate.all_cut_positive_ceiling
            );
            print_witness(
                "zero_pair_ceiling",
                depth_aggregate.zero_pair_ceiling
            );
            print_witness(
                "low_band_bound",
                depth_aggregate.low_band_bound
            );
            print_witness("direct", depth_aggregate.direct);
            std::cout << '\n';
        }
        std::cout << "SU2_SEVEN_PARITY_ORBITS PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SEVEN_PARITY_ORBITS FAILURE: "
                  << error.what() << '\n';
        return 1;
    }
}
#endif
