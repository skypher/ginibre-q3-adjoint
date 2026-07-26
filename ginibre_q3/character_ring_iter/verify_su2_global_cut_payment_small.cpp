#include <functional>

#define main analyze_su2_global_cut_payment_main
#include "analyze_su2_global_cut_payment.cpp"
#undef main

namespace {

struct ExhaustiveState {
    std::atomic<std::uint64_t> rows{0U};
    std::atomic<std::uint64_t> active{0U};
    std::atomic<bool> failed{false};
    std::mutex mutex;
    Witness minimum_direct;
    Witness failure;
};

void check_multiset(
    int level,
    int factors,
    const std::vector<int>& counts,
    ExhaustiveState& state
) {
    int total_parity = 0;
    int used_labels = 0;
    for (int label = 1; label <= level; ++label) {
        const int count = counts[static_cast<std::size_t>(label)];
        total_parity ^= (label & 1) & (count & 1);
        if (count != 0) {
            ++used_labels;
        }
    }
    if (total_parity != 0) {
        return;
    }

    std::vector<int> labels;
    labels.reserve(static_cast<std::size_t>(factors));
    std::vector<int> used;
    for (int label = 1; label <= level; ++label) {
        const int count = counts[static_cast<std::size_t>(label)];
        if (count != 0) {
            used.push_back(label);
        }
        for (int copy = 0; copy < count; ++copy) {
            labels.push_back(label);
        }
    }

    const std::uint64_t assignments =
        std::uint64_t{1} << used_labels;
    for (std::uint64_t assignment = 0U;
         assignment < assignments && !state.failed.load();
         ++assignment) {
        std::uint64_t minus_mask = 0U;
        int position = 0;
        int minus_count = 0;
        int used_index = 0;
        for (int label = 1; label <= level; ++label) {
            const int count =
                counts[static_cast<std::size_t>(label)];
            if (count == 0) {
                continue;
            }
            const bool minus =
                ((assignment >> used_index) & 1U) != 0U;
            if (minus) {
                minus_count += count;
                for (int copy = 0; copy < count; ++copy) {
                    minus_mask |= std::uint64_t{1}
                        << (position + copy);
                }
            }
            position += count;
            ++used_index;
        }
        if ((minus_count & 1) != 0) {
            continue;
        }

        const Evaluation result = evaluate(
            level, labels, minus_mask
        );
        state.rows.fetch_add(1U);
        if (!result.active) {
            continue;
        }
        state.active.fetch_add(1U);
        std::lock_guard<std::mutex> lock(state.mutex);
        if (!state.minimum_direct.initialized
            || result.direct_margin
                < state.minimum_direct.margin) {
            state.minimum_direct = Witness{
                true,
                result.direct_margin,
                level,
                minus_mask,
                labels,
                result
            };
        }
        if (result.direct_margin < 0) {
            state.failure = Witness{
                true,
                result.direct_margin,
                level,
                minus_mask,
                labels,
                result
            };
            state.failed.store(true);
        }
    }
}

void enumerate_task(
    int level,
    int factors,
    ExhaustiveState& state
) {
    std::vector<int> counts(
        static_cast<std::size_t>(level + 1), 0
    );
    const std::function<void(int, int)> recurse =
        [&](int label, int remaining) {
            if (state.failed.load()) {
                return;
            }
            if (label == level) {
                counts[static_cast<std::size_t>(label)] =
                    remaining;
                check_multiset(
                    level, factors, counts, state
                );
                return;
            }
            for (int count = 0; count <= remaining; ++count) {
                counts[static_cast<std::size_t>(label)] = count;
                recurse(label + 1, remaining - count);
            }
        };
    recurse(1, factors);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: verify_su2_global_cut_payment_small "
                "maximum_level maximum_factors"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_factors =
            parse_positive(argv[2], "maximum factors");
        if (maximum_level < 2 || maximum_factors < 2
            || maximum_factors > 20) {
            throw std::runtime_error(
                "require level>=2 and 2<=factors<=20"
            );
        }

        struct Task {
            int level;
            int factors;
        };
        std::vector<Task> tasks;
        for (int level = 2; level <= maximum_level; ++level) {
            for (int factors = 2;
                 factors <= maximum_factors; ++factors) {
                tasks.push_back(Task{level, factors});
            }
        }
        const unsigned workers = std::min<unsigned>(
            worker_limit(maximum_level, maximum_factors),
            static_cast<unsigned>(tasks.size())
        );
        ExhaustiveState state;
        std::atomic<std::size_t> next_task{0U};
        std::vector<std::jthread> pool;
        pool.reserve(workers);
        for (unsigned worker = 0; worker < workers; ++worker) {
            pool.emplace_back([&]() {
                while (!state.failed.load()) {
                    const std::size_t index =
                        next_task.fetch_add(1U);
                    if (index >= tasks.size()) {
                        return;
                    }
                    enumerate_task(
                        tasks[index].level,
                        tasks[index].factors,
                        state
                    );
                }
            });
        }
        pool.clear();

        std::cout << "SU2_GLOBAL_CUT_PAYMENT_SMALL"
                  << " maximum_level=" << maximum_level
                  << " maximum_factors=" << maximum_factors
                  << " rows=" << state.rows.load()
                  << " active=" << state.active.load()
                  << " workers=" << workers << '\n';
        print_witness("minimum_direct", state.minimum_direct);
        if (state.failure.initialized) {
            print_witness("counterexample", state.failure);
            std::cout
                << "SU2_GLOBAL_CUT_PAYMENT_SMALL FAIL\n";
            return EXIT_FAILURE;
        }
        std::cout << "SU2_GLOBAL_CUT_PAYMENT_SMALL"
                  << " counterexamples=0 result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_GLOBAL_CUT_PAYMENT_SMALL FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
