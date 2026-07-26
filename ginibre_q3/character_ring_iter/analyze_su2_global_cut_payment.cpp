#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using Integer = __int128_t;

struct Evaluation {
    bool active = false;
    bool exact_failure = false;
    Integer exact_margin = 0;
    Integer direct_margin = 0;
    Integer invariant = 0;
    Integer positive = 0;
    Integer negative = 0;
    Integer local = 0;
    std::uint64_t selected_mask = 0U;
};

struct Witness {
    bool initialized = false;
    Integer margin = 0;
    int level = 0;
    std::uint64_t minus_mask = 0U;
    std::vector<int> labels;
    Evaluation evaluation;
};

std::uint64_t parse_u64(const char* text, const char* name) {
    char* end = nullptr;
    const unsigned long long value = std::strtoull(text, &end, 10);
    if (end == text || *end != '\0') {
        throw std::runtime_error(
            std::string("invalid ") + name
        );
    }
    return static_cast<std::uint64_t>(value);
}

int parse_positive(const char* text, const char* name) {
    const std::uint64_t value = parse_u64(text, name);
    if (value == 0U
        || value > static_cast<std::uint64_t>(
            std::numeric_limits<int>::max()
        )) {
        throw std::runtime_error(
            std::string(name) + " must be positive"
        );
    }
    return static_cast<int>(value);
}

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U))
        * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U))
        * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

std::string integer_string(Integer value) {
    if (value == 0) {
        return "0";
    }
    const bool negative = value < 0;
    if (negative) {
        value = -value;
    }
    std::string result;
    while (value != 0) {
        const int digit = static_cast<int>(value % 10);
        result.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    if (negative) {
        result.push_back('-');
    }
    std::reverse(result.begin(), result.end());
    return result;
}

bool disjoint_support(
    const std::vector<int>& labels,
    std::uint64_t minus_mask
) {
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        for (std::size_t j = i + 1U; j < labels.size(); ++j) {
            if (labels[i] == labels[j]
                && (((minus_mask >> i) & 1U)
                    != ((minus_mask >> j) & 1U))) {
                return false;
            }
        }
    }
    return true;
}

std::vector<int> probe_outputs(
    int level,
    int common_parity
) {
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
        const int high =
            level - reflected_parity - 2 * j;
        if (high >= 0
            && std::find(result.begin(), result.end(), high)
                == result.end()) {
            result.push_back(high);
        }
    }
    return result;
}

Evaluation evaluate(
    int level,
    const std::vector<int>& labels,
    std::uint64_t minus_mask
) {
    const int factors = static_cast<int>(labels.size());
    const std::uint64_t state_count = 1ULL << factors;
    const std::uint64_t full_mask = state_count - 1U;
    const std::size_t width =
        static_cast<std::size_t>(level + 1);
    std::vector<std::int64_t> states(
        static_cast<std::size_t>(state_count) * width,
        0
    );
    const auto at = [
        &states, width
    ](std::uint64_t mask, int output) -> std::int64_t& {
        return states[
            static_cast<std::size_t>(mask) * width
            + static_cast<std::size_t>(output)
        ];
    };
    at(0U, 0) = 1;
    for (std::uint64_t mask = 1U; mask < state_count; ++mask) {
        const int index = std::countr_zero(mask);
        const std::uint64_t previous =
            mask & (mask - 1U);
        const int label = labels[static_cast<std::size_t>(index)];
        for (int source = 0; source <= level; ++source) {
            const std::int64_t multiplicity =
                at(previous, source);
            if (multiplicity == 0) {
                continue;
            }
            const int upper = std::min(
                source + label,
                2 * level - source - label
            );
            for (int output = std::abs(source - label);
                 output <= upper; output += 2) {
                std::int64_t updated = 0;
                if (__builtin_add_overflow(
                        at(mask, output),
                        multiplicity,
                        &updated
                    )) {
                    throw std::overflow_error(
                        "fusion multiplicity overflow"
                    );
                }
                at(mask, output) = updated;
            }
        }
    }

    Evaluation result;
    result.invariant = at(full_mask, 0);
    Integer maximum_negative = 0;
    for (std::uint64_t mask = 1U; mask < full_mask; ++mask) {
        const std::uint64_t complement = full_mask ^ mask;
        if (mask > complement) {
            continue;
        }
        const Integer weight =
            static_cast<Integer>(at(mask, 0))
            * static_cast<Integer>(at(complement, 0));
        const bool negative =
            (std::popcount(mask & minus_mask) & 1) != 0;
        if (negative) {
            result.negative += weight;
            if (weight > maximum_negative) {
                maximum_negative = weight;
                result.selected_mask = mask;
            }
        } else {
            result.positive += weight;
        }
    }
    result.exact_margin =
        result.invariant + result.positive - result.negative;
    result.exact_failure = result.exact_margin < 0;
    if (maximum_negative == 0) {
        result.direct_margin =
            result.positive - result.negative;
        return result;
    }

    result.active = true;
    const std::uint64_t complement =
        full_mask ^ result.selected_mask;
    int common_parity = 0;
    for (int index = 0; index < factors; ++index) {
        if (((result.selected_mask >> index) & 1U) != 0U) {
            common_parity ^= (
                labels[static_cast<std::size_t>(index)] & 1
            );
        }
    }
    for (const int output :
         probe_outputs(level, common_parity)) {
        result.local +=
            static_cast<Integer>(
                at(result.selected_mask, output)
            ) * static_cast<Integer>(
                at(complement, output)
            );
    }
    result.direct_margin =
        result.local + result.positive - result.negative;
    return result;
}

void print_labels(const std::vector<int>& labels) {
    std::cout << '[';
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        if (i != 0U) {
            std::cout << ',';
        }
        std::cout << labels[i];
    }
    std::cout << ']';
}

void print_witness(const char* name, const Witness& witness) {
    std::cout << name << '=';
    if (!witness.initialized) {
        std::cout << "none\n";
        return;
    }
    std::cout << integer_string(witness.margin)
              << " level=" << witness.level
              << " minus_mask=" << witness.minus_mask
              << " labels=";
    print_labels(witness.labels);
    const Evaluation& e = witness.evaluation;
    std::cout << " selected_mask=" << e.selected_mask
              << " (N,U,T,L)=("
              << integer_string(e.invariant) << ','
              << integer_string(e.positive) << ','
              << integer_string(e.negative) << ','
              << integer_string(e.local) << ")\n";
}

unsigned worker_limit(
    int maximum_level,
    int maximum_factors
) {
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
    const std::uint64_t state_bytes =
        (std::uint64_t{1} << maximum_factors)
        * static_cast<std::uint64_t>(maximum_level + 1)
        * sizeof(std::int64_t);
    const std::uint64_t kib_per_worker = std::max<std::uint64_t>(
        64U * 1024U,
        (state_bytes + 1023U) / 1024U
    );
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
        if (argc != 5) {
            throw std::runtime_error(
                "usage: analyze_su2_global_cut_payment "
                "seed trials maximum_level maximum_factors"
            );
        }
        const std::uint64_t seed = parse_u64(argv[1], "seed");
        const std::uint64_t trials = parse_u64(argv[2], "trials");
        const int maximum_level =
            parse_positive(argv[3], "maximum level");
        const int maximum_factors =
            parse_positive(argv[4], "maximum factors");
        if (trials == 0U || maximum_level < 2
            || maximum_factors < 2 || maximum_factors > 20) {
            throw std::runtime_error(
                "require trials>0, level>=2, and 2<=factors<=20"
            );
        }

        const unsigned workers = worker_limit(
            maximum_level, maximum_factors
        );
        std::atomic<std::uint64_t> next_trial{0U};
        std::atomic<std::uint64_t> evaluated{0U};
        std::atomic<std::uint64_t> active{0U};
        std::atomic<std::uint64_t> interior_active{0U};
        std::atomic<bool> failed{false};
        std::mutex witness_mutex;
        Witness minimum_exact;
        Witness minimum_direct;
        Witness minimum_interior_direct;
        std::vector<Witness> minimum_direct_by_factors(
            static_cast<std::size_t>(maximum_factors + 1)
        );
        Witness failure;
        std::vector<std::jthread> pool;
        pool.reserve(workers);
        for (unsigned worker = 0; worker < workers; ++worker) {
            pool.emplace_back([&]() {
                while (!failed.load()) {
                    const std::uint64_t trial =
                        next_trial.fetch_add(1U);
                    if (trial >= trials) {
                        return;
                    }
                    std::uint64_t random = splitmix64(
                        seed ^ (trial * 0xd1342543de82ef95ULL)
                    );
                    const int level = 2 + static_cast<int>(
                        random % static_cast<std::uint64_t>(
                            maximum_level - 1
                        )
                    );
                    random = splitmix64(random);
                    const int factors = 2 + static_cast<int>(
                        random % static_cast<std::uint64_t>(
                            maximum_factors - 1
                        )
                    );
                    std::vector<int> labels(
                        static_cast<std::size_t>(factors)
                    );
                    std::uint64_t minus_mask = 0U;
                    int total_parity = 0;
                    for (int index = 0; index < factors; ++index) {
                        random = splitmix64(random);
                        labels[static_cast<std::size_t>(index)] =
                            1 + static_cast<int>(
                                random
                                % static_cast<std::uint64_t>(level)
                            );
                        total_parity ^=
                            labels[static_cast<std::size_t>(index)] & 1;
                        random = splitmix64(random);
                        if ((random & 1U) != 0U) {
                            minus_mask |= 1ULL << index;
                        }
                    }
                    if ((std::popcount(minus_mask) & 1) != 0) {
                        minus_mask ^= 1ULL << (factors - 1);
                    }
                    if (total_parity != 0
                        || !disjoint_support(labels, minus_mask)) {
                        continue;
                    }

                    const Evaluation result = evaluate(
                        level, labels, minus_mask
                    );
                    evaluated.fetch_add(1U);
                    if (result.active) {
                        active.fetch_add(1U);
                        if (result.local < result.invariant) {
                            interior_active.fetch_add(1U);
                        }
                    }
                    std::lock_guard<std::mutex> lock(witness_mutex);
                    const auto update = [&](
                        Witness& target, Integer margin
                    ) {
                        if (!target.initialized
                            || margin < target.margin) {
                            target = Witness{
                                true,
                                margin,
                                level,
                                minus_mask,
                                labels,
                                result
                            };
                        }
                    };
                    update(minimum_exact, result.exact_margin);
                    if (result.active) {
                        update(
                            minimum_direct,
                            result.direct_margin
                        );
                        update(
                            minimum_direct_by_factors[
                                static_cast<std::size_t>(factors)
                            ],
                            result.direct_margin
                        );
                        if (result.local < result.invariant) {
                            update(
                                minimum_interior_direct,
                                result.direct_margin
                            );
                        }
                    }
                    if (result.exact_failure) {
                        failure = Witness{
                            true,
                            result.exact_margin,
                            level,
                            minus_mask,
                            labels,
                            result
                        };
                        failed.store(true);
                    }
                }
            });
        }
        pool.clear();

        std::cout << "SU2_GLOBAL_CUT_PAYMENT"
                  << " seed=" << seed
                  << " trials=" << trials
                  << " evaluated=" << evaluated.load()
                  << " active=" << active.load()
                  << " interior_active=" << interior_active.load()
                  << " maximum_level=" << maximum_level
                  << " maximum_factors=" << maximum_factors
                  << " workers=" << workers << '\n';
        print_witness("minimum_exact", minimum_exact);
        print_witness("minimum_direct", minimum_direct);
        print_witness(
            "minimum_interior_direct",
            minimum_interior_direct
        );
        for (int factors = 2;
             factors <= maximum_factors; ++factors) {
            std::ostringstream name;
            name << "minimum_direct_n" << factors;
            print_witness(
                name.str().c_str(),
                minimum_direct_by_factors[
                    static_cast<std::size_t>(factors)
                ]
            );
        }
        if (failure.initialized) {
            print_witness("exact_counterexample", failure);
            std::cout << "SU2_GLOBAL_CUT_PAYMENT FAIL\n";
            return EXIT_FAILURE;
        }
        std::cout << "SU2_GLOBAL_CUT_PAYMENT"
                  << " exact_counterexamples=0 result=PASS_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_GLOBAL_CUT_PAYMENT FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
