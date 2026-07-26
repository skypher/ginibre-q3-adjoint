#include <algorithm>
#include <atomic>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

struct Evaluation {
    bool active = false;
    bool gks_failure = false;
    bool direct_failure = false;
    Integer twice_exact_margin = 0;
    Integer twice_anchor_margin = 0;
    Integer twice_direct_margin = 0;
    Integer invariant = 0;
    Integer positive_oriented = 0;
    Integer negative_oriented = 0;
    Integer twice_local = 0;
    Integer maximum_negative = 0;
    std::size_t selected_state = 0U;
};

struct Case {
    int level = 0;
    std::vector<int> labels;
    std::vector<int> counts;
    std::vector<bool> minus;
};

struct Witness {
    bool initialized = false;
    Integer margin = 0;
    Case input;
    Evaluation evaluation;
};

std::uint64_t parse_u64(const char* text, const char* name) {
    char* end = nullptr;
    const unsigned long long value = std::strtoull(text, &end, 10);
    if (end == text || *end != '\0') {
        throw std::runtime_error(std::string("invalid ") + name);
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

std::vector<int> probe_outputs(int level, int common_parity) {
    std::vector<int> result;
    for (int j = 0; j <= 4; ++j) {
        const int low = common_parity + 2 * j;
        if (low <= level) {
            result.push_back(low);
        }
    }
    const int reflected_parity = (level + common_parity) & 1;
    for (int j = 0; j <= 4; ++j) {
        const int high = level - reflected_parity - 2 * j;
        if (high >= 0
            && std::find(result.begin(), result.end(), high)
                == result.end()) {
            result.push_back(high);
        }
    }
    return result;
}

std::vector<std::vector<Integer>> binomial_table(int maximum) {
    std::vector<std::vector<Integer>> choose(
        static_cast<std::size_t>(maximum + 1)
    );
    for (int n = 0; n <= maximum; ++n) {
        choose[static_cast<std::size_t>(n)].assign(
            static_cast<std::size_t>(n + 1), Integer{0}
        );
        choose[static_cast<std::size_t>(n)][0] = 1;
        choose[static_cast<std::size_t>(n)][
            static_cast<std::size_t>(n)
        ] = 1;
        for (int r = 1; r < n; ++r) {
            choose[static_cast<std::size_t>(n)][
                static_cast<std::size_t>(r)
            ] =
                choose[static_cast<std::size_t>(n - 1)][
                    static_cast<std::size_t>(r - 1)
                ]
                + choose[static_cast<std::size_t>(n - 1)][
                    static_cast<std::size_t>(r)
                ];
        }
    }
    return choose;
}

Evaluation evaluate(
    const Case& input,
    const std::vector<std::vector<Integer>>& choose
) {
    const int level = input.level;
    const std::size_t types = input.labels.size();
    const std::size_t width = static_cast<std::size_t>(level + 1);
    std::vector<std::size_t> strides(types, 1U);
    std::size_t state_count = 1U;
    for (std::size_t type = 0U; type < types; ++type) {
        strides[type] = state_count;
        const std::size_t radix = static_cast<std::size_t>(
            input.counts[type] + 1
        );
        if (state_count
            > std::numeric_limits<std::size_t>::max() / radix) {
            throw std::overflow_error("count-state space overflow");
        }
        state_count *= radix;
    }
    if (state_count
        > std::numeric_limits<std::size_t>::max() / width) {
        throw std::overflow_error("fusion-state space overflow");
    }

    std::vector<Integer> states(state_count * width);
    const auto at = [
        &states, width
    ](std::size_t state, int output) -> Integer& {
        return states[
            state * width + static_cast<std::size_t>(output)
        ];
    };
    at(0U, 0) = 1;
    for (std::size_t state = 1U; state < state_count; ++state) {
        std::size_t selected = types;
        for (std::size_t type = 0U; type < types; ++type) {
            const std::size_t radix = static_cast<std::size_t>(
                input.counts[type] + 1
            );
            if (((state / strides[type]) % radix) != 0U) {
                selected = type;
                break;
            }
        }
        if (selected == types) {
            throw std::logic_error("nonzero state has no digit");
        }
        const std::size_t previous = state - strides[selected];
        const int label = input.labels[selected];
        for (int source = 0; source <= level; ++source) {
            const Integer& multiplicity = at(previous, source);
            if (multiplicity == 0) {
                continue;
            }
            const int upper = std::min(
                source + label,
                2 * level - source - label
            );
            for (int output = std::abs(source - label);
                 output <= upper; output += 2) {
                at(state, output) += multiplicity;
            }
        }
    }

    Evaluation result;
    const std::size_t full_state = state_count - 1U;
    result.invariant = at(full_state, 0);
    Integer maximum_negative = 0;
    for (std::size_t state = 1U; state < full_state; ++state) {
        std::size_t work = state;
        bool negative = false;
        Integer orbit_multiplicity = 1;
        for (std::size_t type = types; type-- > 0U;) {
            const std::size_t radix = static_cast<std::size_t>(
                input.counts[type] + 1
            );
            const int selected = static_cast<int>(
                (work / strides[type]) % radix
            );
            if (input.minus[type] && ((selected & 1) != 0)) {
                negative = !negative;
            }
            orbit_multiplicity *=
                choose[static_cast<std::size_t>(input.counts[type])][
                    static_cast<std::size_t>(selected)
                ];
        }
        const std::size_t complement = full_state - state;
        const Integer weight =
            at(state, 0) * at(complement, 0);
        const Integer indexed_weight = orbit_multiplicity * weight;
        if (negative) {
            result.negative_oriented += indexed_weight;
            if (weight > maximum_negative) {
                maximum_negative = weight;
                result.selected_state = state;
            }
        } else {
            result.positive_oriented += indexed_weight;
        }
    }

    result.twice_exact_margin =
        2 * result.invariant
        + result.positive_oriented
        - result.negative_oriented;
    result.maximum_negative = maximum_negative;
    result.twice_anchor_margin =
        result.positive_oriented
        + 2 * maximum_negative
        - result.negative_oriented;
    result.gks_failure = result.twice_exact_margin < 0;
    if (maximum_negative == 0) {
        result.twice_direct_margin =
            result.positive_oriented - result.negative_oriented;
        result.direct_failure = result.twice_direct_margin < 0;
        return result;
    }

    result.active = true;
    const std::size_t complement =
        full_state - result.selected_state;
    int common_parity = 0;
    for (std::size_t type = 0U; type < types; ++type) {
        const std::size_t radix = static_cast<std::size_t>(
            input.counts[type] + 1
        );
        const int selected = static_cast<int>(
            (result.selected_state / strides[type]) % radix
        );
        if ((selected & 1) != 0) {
            common_parity ^= input.labels[type] & 1;
        }
    }
    Integer local = 0;
    for (const int output : probe_outputs(level, common_parity)) {
        local += at(result.selected_state, output)
            * at(complement, output);
    }
    result.twice_local = 2 * local;
    result.twice_direct_margin =
        result.positive_oriented
        + result.twice_local
        - result.negative_oriented;
    result.direct_failure = result.twice_direct_margin < 0;
    return result;
}

std::vector<int> state_digits(
    const Case& input,
    std::size_t state
) {
    std::vector<int> result(input.labels.size(), 0);
    for (std::size_t type = 0U; type < input.labels.size(); ++type) {
        const std::size_t radix = static_cast<std::size_t>(
            input.counts[type] + 1
        );
        result[type] = static_cast<int>(state % radix);
        state /= radix;
    }
    return result;
}

void print_vector(const std::vector<int>& values) {
    std::cout << '[';
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << ']';
}

void print_signs(const std::vector<bool>& values) {
    std::cout << '[';
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << (values[index] ? '-' : '+');
    }
    std::cout << ']';
}

void print_witness(const char* name, const Witness& witness) {
    std::cout << name << '=';
    if (!witness.initialized) {
        std::cout << "none\n";
        return;
    }
    const Case& input = witness.input;
    const Evaluation& result = witness.evaluation;
    std::cout << witness.margin
              << " level=" << input.level
              << " labels=";
    print_vector(input.labels);
    std::cout << " counts=";
    print_vector(input.counts);
    std::cout << " signs=";
    print_signs(input.minus);
    std::cout << " selected_counts=";
    print_vector(state_digits(input, result.selected_state));
    std::cout << " (2N,Uo,To,2w,2L)=("
              << 2 * result.invariant << ','
              << result.positive_oriented << ','
              << result.negative_oriented << ','
              << 2 * result.maximum_negative << ','
              << result.twice_local << ")\n";
}

unsigned worker_limit(
    int maximum_level,
    int maximum_factors,
    int maximum_types
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
    std::uint64_t value = 0U;
    std::string unit;
    std::uint64_t available_kib = 0U;
    while (memory >> key >> value >> unit) {
        if (key == "MemAvailable:") {
            available_kib = value;
            break;
        }
    }
    const std::uint64_t balanced_radix =
        static_cast<std::uint64_t>(
            maximum_factors / std::max(1, maximum_types) + 2
        );
    std::uint64_t states = 1U;
    for (int type = 0; type < maximum_types; ++type) {
        if (states > std::numeric_limits<std::uint64_t>::max()
                / balanced_radix) {
            states = std::numeric_limits<std::uint64_t>::max();
            break;
        }
        states *= balanced_radix;
    }
    const std::uint64_t estimated_bytes =
        states * static_cast<std::uint64_t>(maximum_level + 1)
        * 48U;
    const std::uint64_t kib_per_worker = std::max<std::uint64_t>(
        128U * 1024U,
        (estimated_bytes + 1023U) / 1024U
    );
    const unsigned memory_limit = available_kib == 0U
        ? hardware
        : static_cast<unsigned>(std::max<std::uint64_t>(
            1U, available_kib / kib_per_worker
        ));
    return std::max(1U, std::min(cpu_limit, memory_limit));
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2
            && std::string(argv[1]) == "--anchor-counterexample") {
            const Case input{
                4,
                {1, 2, 3, 4},
                {10, 18, 18, 8},
                {false, false, true, false}
            };
            const Evaluation result = evaluate(
                input, binomial_table(54)
            );
            print_witness(
                "anchor_counterexample",
                Witness{
                    true,
                    result.twice_anchor_margin,
                    input,
                    result
                }
            );
            const bool passed =
                result.invariant == Integer{"208971104256"}
                && result.positive_oriented
                    == Integer{"156853697876420027319578880"}
                && result.negative_oriented
                    == Integer{"156853697876420445125738496"}
                && result.maximum_negative
                    == Integer{"69657034752"}
                && result.twice_local
                    == Integer{"417942208512"}
                && result.twice_anchor_margin
                    == Integer{"-278492090112"}
                && result.twice_direct_margin
                    == Integer{"136048896"}
                && result.twice_exact_margin
                    == Integer{"136048896"}
                && state_digits(input, result.selected_state)
                    == std::vector<int>({1, 0, 1, 1});
            std::cout
                << "SU2_GCP_ANCHOR_COUNTEREXAMPLE "
                << "anchor_only=FAIL gcp=PASS gks=PASS result="
                << (passed ? "PASS\n" : "FAIL\n");
            return passed ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        if (argc != 6) {
            throw std::runtime_error(
                "usage: analyze_su2_global_cut_payment_counts "
                "seed trials maximum_level maximum_factors maximum_types"
                " | --anchor-counterexample"
            );
        }
        const std::uint64_t seed = parse_u64(argv[1], "seed");
        const std::uint64_t trials = parse_u64(argv[2], "trials");
        const int maximum_level =
            parse_positive(argv[3], "maximum level");
        const int maximum_factors =
            parse_positive(argv[4], "maximum factors");
        const int maximum_types =
            parse_positive(argv[5], "maximum types");
        if (trials == 0U || maximum_level < 2
            || maximum_factors < 2 || maximum_factors > 200
            || maximum_types > 6) {
            throw std::runtime_error(
                "require trials>0, level>=2, 2<=factors<=200, "
                "and 1<=types<=6"
            );
        }

        const auto choose = binomial_table(maximum_factors);
        const unsigned workers = worker_limit(
            maximum_level, maximum_factors, maximum_types
        );
        std::atomic<std::uint64_t> next_trial{0U};
        std::atomic<std::uint64_t> evaluated{0U};
        std::atomic<std::uint64_t> active{0U};
        std::atomic<bool> stop{false};
        std::mutex witness_mutex;
        Witness minimum_exact;
        Witness minimum_anchor;
        Witness minimum_direct;
        Witness gks_failure;
        Witness direct_failure;
        std::vector<std::jthread> pool;
        pool.reserve(workers);
        for (unsigned worker = 0U; worker < workers; ++worker) {
            pool.emplace_back([&]() {
                while (!stop.load()) {
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
                    const int type_limit = std::min({
                        maximum_types, level, factors
                    });
                    random = splitmix64(random);
                    const int types = 1 + static_cast<int>(
                        random % static_cast<std::uint64_t>(type_limit)
                    );

                    Case input;
                    input.level = level;
                    std::vector<int> label_pool(
                        static_cast<std::size_t>(level)
                    );
                    std::iota(label_pool.begin(), label_pool.end(), 1);
                    for (int index = 0; index < types; ++index) {
                        random = splitmix64(random);
                        const int other = index + static_cast<int>(
                            random % static_cast<std::uint64_t>(
                                level - index
                            )
                        );
                        std::swap(
                            label_pool[static_cast<std::size_t>(index)],
                            label_pool[static_cast<std::size_t>(other)]
                        );
                    }
                    input.labels.assign(
                        label_pool.begin(),
                        label_pool.begin() + types
                    );
                    std::sort(input.labels.begin(), input.labels.end());
                    input.counts.assign(
                        static_cast<std::size_t>(types), 1
                    );
                    for (int remaining = factors - types;
                         remaining > 0; --remaining) {
                        random = splitmix64(random);
                        const std::size_t type =
                            static_cast<std::size_t>(
                                random
                                % static_cast<std::uint64_t>(types)
                            );
                        ++input.counts[type];
                    }
                    input.minus.assign(
                        static_cast<std::size_t>(types), false
                    );
                    int minus_count = 0;
                    int total_label_parity = 0;
                    for (int type = 0; type < types; ++type) {
                        random = splitmix64(random);
                        input.minus[static_cast<std::size_t>(type)] =
                            (random & 1U) != 0U;
                        if (input.minus[static_cast<std::size_t>(type)]) {
                            minus_count += input.counts[
                                static_cast<std::size_t>(type)
                            ];
                        }
                        total_label_parity ^=
                            (input.labels[static_cast<std::size_t>(type)]
                             * input.counts[static_cast<std::size_t>(type)])
                            & 1;
                    }
                    if (minus_count == 0 || (minus_count & 1) != 0
                        || total_label_parity != 0) {
                        continue;
                    }

                    const Evaluation result = evaluate(input, choose);
                    evaluated.fetch_add(1U);
                    if (result.active) {
                        active.fetch_add(1U);
                    }
                    std::lock_guard<std::mutex> lock(witness_mutex);
                    const auto update = [&](
                        Witness& target, const Integer& margin
                    ) {
                        if (!target.initialized
                            || margin < target.margin) {
                            target = Witness{
                                true, margin, input, result
                            };
                        }
                    };
                    update(
                        minimum_exact, result.twice_exact_margin
                    );
                    if (result.active) {
                        update(
                            minimum_anchor,
                            result.twice_anchor_margin
                        );
                        update(
                            minimum_direct,
                            result.twice_direct_margin
                        );
                    }
                    if (result.gks_failure) {
                        update(
                            gks_failure, result.twice_exact_margin
                        );
                        stop.store(true);
                    } else if (result.direct_failure) {
                        update(
                            direct_failure,
                            result.twice_direct_margin
                        );
                        stop.store(true);
                    }
                }
            });
        }
        pool.clear();

        std::cout << "SU2_GLOBAL_CUT_PAYMENT_COUNTS"
                  << " seed=" << seed
                  << " trials=" << trials
                  << " evaluated=" << evaluated.load()
                  << " active=" << active.load()
                  << " maximum_level=" << maximum_level
                  << " maximum_factors=" << maximum_factors
                  << " maximum_types=" << maximum_types
                  << " workers=" << workers << '\n';
        print_witness("minimum_twice_exact", minimum_exact);
        print_witness("minimum_twice_anchor", minimum_anchor);
        print_witness("minimum_twice_direct", minimum_direct);
        print_witness("gks_counterexample", gks_failure);
        print_witness("gcp_counterexample", direct_failure);
        if (gks_failure.initialized) {
            std::cout
                << "SU2_GLOBAL_CUT_PAYMENT_COUNTS GKS_FAIL\n";
            return EXIT_FAILURE;
        }
        if (direct_failure.initialized) {
            std::cout
                << "SU2_GLOBAL_CUT_PAYMENT_COUNTS "
                << "GCP_FAIL_GKS_NONNEGATIVE\n";
            return EXIT_SUCCESS;
        }
        std::cout
            << "SU2_GLOBAL_CUT_PAYMENT_COUNTS "
            << "counterexamples=0 result=PASS_DISCOVERY\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_GLOBAL_CUT_PAYMENT_COUNTS FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
