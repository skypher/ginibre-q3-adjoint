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

void merge_witness(Witness& target, const Witness& source) {
    if (source.initialized
        && (!target.initialized || source.value < target.value)) {
        target = source;
    }
}

void merge_stratum(
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
    merge_witness(target.direct, source.direct);
    merge_witness(target.exact, source.exact);
    target.maximum_ratio = std::max(
        target.maximum_ratio, source.maximum_ratio
    );
}

std::array<OrbitStratum, 2> scan_level(
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
            std::vector<unsigned int> maximum_masks;
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
                    negative += term;
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

unsigned worker_limit_orbits() {
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
        for (std::size_t index = 0U;
             index < orbits.size(); ++index) {
            const Orbit& orbit = orbits[index];
            for (std::size_t depth = 0U; depth < 2U; ++depth) {
                const OrbitStratum& stratum =
                    totals[index][depth];
                if (stratum.cases == 0U) {
                    continue;
                }
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
                print_witness("direct", stratum.direct);
                print_witness("exact", stratum.exact);
                std::cout << " max_ceil_T_over_d="
                          << stratum.maximum_ratio << '\n';
            }
        }
        std::cout << "SU2_SEVEN_PARITY_ORBITS PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_SEVEN_PARITY_ORBITS FAILURE: "
                  << error.what() << '\n';
        return 1;
    }
}
