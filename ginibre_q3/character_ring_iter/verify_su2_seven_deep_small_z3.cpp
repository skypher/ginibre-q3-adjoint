#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#define SU2_SEVEN_RESIDUAL_NO_MAIN
#include "verify_su2_seven_residual_z3.cpp"
#undef SU2_SEVEN_RESIDUAL_NO_MAIN

namespace {

struct DeepSmallTask {
    int orbit_index;
    int level_parity;
    int selected_orbit;
    int equality_mask;
    int selected_rank;
};

void append_partition_masks(
    const std::vector<int>& indices,
    const std::array<std::array<int, 7>, 7>& pair_bits,
    std::size_t position,
    int maximum_block,
    std::vector<int>& blocks,
    std::vector<int>& masks
) {
    if (position == indices.size()) {
        int mask = 0;
        for (std::size_t first = 0U;
             first < indices.size(); ++first) {
            for (std::size_t second = first + 1U;
                 second < indices.size(); ++second) {
                if (blocks[first] != blocks[second]) {
                    continue;
                }
                const int bit = pair_bits[
                    static_cast<std::size_t>(indices[first])
                ][static_cast<std::size_t>(indices[second])];
                if (bit < 0) {
                    std::abort();
                }
                mask |= 1 << bit;
            }
        }
        masks.push_back(mask);
        return;
    }

    for (int block = 0; block <= maximum_block + 1; ++block) {
        blocks[position] = block;
        append_partition_masks(
            indices,
            pair_bits,
            position + 1U,
            std::max(maximum_block, block),
            blocks,
            masks
        );
    }
}

std::vector<int> group_partition_masks(
    const std::vector<int>& indices,
    const std::array<std::array<int, 7>, 7>& pair_bits
) {
    if (indices.size() <= 1U) {
        return {0};
    }
    std::vector<int> blocks(indices.size(), 0);
    std::vector<int> result;
    blocks[0] = 0;
    append_partition_masks(
        indices, pair_bits, 1U, 0, blocks, result
    );
    std::sort(result.begin(), result.end());
    result.erase(
        std::unique(result.begin(), result.end()),
        result.end()
    );
    return result;
}

std::vector<int> equality_partition_masks(
    const ResidualOrbit& orbit
) {
    std::array<std::array<int, 7>, 7> pair_bits{};
    for (auto& row : pair_bits) {
        row.fill(-1);
    }
    int equality_index = 0;
    for (int first = 0; first < 7; ++first) {
        for (int second = first + 1; second < 7; ++second) {
            const bool same_sign =
                (first < orbit.minus_count)
                == (second < orbit.minus_count);
            if (!same_sign) {
                continue;
            }
            pair_bits[static_cast<std::size_t>(first)]
                     [static_cast<std::size_t>(second)] =
                equality_index;
            pair_bits[static_cast<std::size_t>(second)]
                     [static_cast<std::size_t>(first)] =
                equality_index;
            ++equality_index;
        }
    }

    std::vector<int> masks{0};
    for (int sign = 0; sign < 2; ++sign) {
        for (int parity = 0; parity < 2; ++parity) {
            std::vector<int> group;
            for (int index = 0; index < 7; ++index) {
                const bool is_minus = index < orbit.minus_count;
                if (is_minus == (sign == 0)
                    && orbit_label_parity(orbit, index) == parity) {
                    group.push_back(index);
                }
            }
            const std::vector<int> group_masks =
                group_partition_masks(group, pair_bits);
            std::vector<int> combined;
            combined.reserve(masks.size() * group_masks.size());
            for (const int prefix : masks) {
                for (const int suffix : group_masks) {
                    combined.push_back(prefix | suffix);
                }
            }
            masks = std::move(combined);
        }
    }
    std::sort(masks.begin(), masks.end());
    masks.erase(
        std::unique(masks.begin(), masks.end()),
        masks.end()
    );
    return masks;
}

bool selected_cut_can_be_active(
    const ResidualOrbit& orbit,
    unsigned int mask
) {
    int parity = 0;
    for (int index = 0; index < 7; ++index) {
        if (((mask >> index) & 1U) != 0U) {
            parity ^= orbit_label_parity(orbit, index);
        }
    }
    return parity == 0;
}

std::vector<DeepSmallTask> deep_small_tasks() {
    std::vector<DeepSmallTask> result;
    for (int orbit_index = 1;
         orbit_index < static_cast<int>(residual_orbits.size());
         ++orbit_index) {
        const ResidualOrbit& orbit =
            residual_orbits[static_cast<std::size_t>(orbit_index)];
        const std::vector<unsigned int> representatives =
            selected_orbit_masks(orbit);
        const std::vector<int> equality_masks =
            equality_partition_masks(orbit);
        const int maximum_rank =
            orbit.minus_count == 4 ? 5 : 10;
        for (int selected_orbit = 0;
             selected_orbit
                 < static_cast<int>(representatives.size());
             ++selected_orbit) {
            if (!selected_cut_can_be_active(
                    orbit,
                    representatives[
                        static_cast<std::size_t>(selected_orbit)
                    ])) {
                continue;
            }
            for (int level_parity = 0;
                 level_parity < 2; ++level_parity) {
                for (int selected_rank = 1;
                     selected_rank <= maximum_rank;
                     ++selected_rank) {
                    for (const int equality_mask : equality_masks) {
                        result.push_back(DeepSmallTask{
                            orbit_index,
                            level_parity,
                            selected_orbit,
                            equality_mask,
                            selected_rank
                        });
                    }
                }
            }
        }
    }
    return result;
}

std::string task_name(const DeepSmallTask& task) {
    return "orbit=" + std::to_string(task.orbit_index)
        + " parity=" + std::to_string(task.level_parity)
        + " selected=" + std::to_string(task.selected_orbit)
        + " equality=" + std::to_string(task.equality_mask)
        + " rank=" + std::to_string(task.selected_rank);
}

}  // namespace

int main(int argc, char** argv) {
    const std::vector<DeepSmallTask> tasks = deep_small_tasks();
    if (argc == 2 && std::string(argv[1]) == "--count") {
        std::cout << "SU2_SEVEN_DEEP_SMALL_Z3"
                  << " tasks=" << tasks.size() << '\n';
        for (int orbit_index = 1;
             orbit_index < static_cast<int>(residual_orbits.size());
             ++orbit_index) {
            const auto masks = equality_partition_masks(
                residual_orbits[
                    static_cast<std::size_t>(orbit_index)
                ]
            );
            std::cout << "orbit=" << orbit_index
                      << " equality_partitions=" << masks.size()
                      << '\n';
        }
        return EXIT_SUCCESS;
    }
    if (argc != 1) {
        std::cerr
            << "usage: verify_su2_seven_deep_small_z3 [--count]\n";
        return EXIT_FAILURE;
    }

    const unsigned workers = std::min<unsigned>(
        worker_limit(),
        static_cast<unsigned>(tasks.size())
    );
    std::atomic<std::size_t> next_task{0U};
    std::atomic<std::size_t> completed{0U};
    std::atomic<bool> failed{false};
    std::mutex diagnostic_mutex;
    std::string diagnostic;

    std::mutex progress_mutex;
    std::condition_variable progress_condition;
    bool finished = false;
    std::jthread reporter([&]() {
        std::unique_lock<std::mutex> lock(progress_mutex);
        while (!progress_condition.wait_for(
            lock,
            std::chrono::seconds(30),
            [&]() { return finished; }
        )) {
            std::cout
                << "SU2_SEVEN_DEEP_SMALL_Z3 progress="
                << completed.load() << '/' << tasks.size()
                << " assigned=" << std::min(
                    next_task.load(), tasks.size()
                )
                << " failed=" << (failed.load() ? 1 : 0)
                << '\n' << std::flush;
        }
    });

    std::cout << "SU2_SEVEN_DEEP_SMALL_Z3"
              << " tasks=" << tasks.size()
              << " workers=" << workers
              << " start=1\n" << std::flush;

    std::vector<std::jthread> pool;
    pool.reserve(workers);
    for (unsigned worker = 0; worker < workers; ++worker) {
        pool.emplace_back([&]() {
            while (!failed.load()) {
                const std::size_t index = next_task.fetch_add(1U);
                if (index >= tasks.size()) {
                    return;
                }
                const DeepSmallTask& task = tasks[index];
                const QueryResult result = verify_residual_query(
                    task.orbit_index,
                    task.level_parity,
                    task.selected_orbit,
                    -1,
                    -1,
                    -1,
                    task.equality_mask,
                    1,
                    task.selected_rank
                );
                if (!result.passed) {
                    {
                        std::lock_guard<std::mutex> lock(
                            diagnostic_mutex
                        );
                        if (diagnostic.empty()) {
                            diagnostic =
                                task_name(task) + " " + result.diagnostic;
                        }
                    }
                    failed.store(true);
                    return;
                }
                completed.fetch_add(1U);
            }
        });
    }
    pool.clear();

    {
        std::lock_guard<std::mutex> lock(progress_mutex);
        finished = true;
    }
    progress_condition.notify_all();
    reporter.join();

    if (failed.load() || completed.load() != tasks.size()) {
        std::cerr << "SU2_SEVEN_DEEP_SMALL_Z3 FAIL"
                  << " completed=" << completed.load()
                  << " tasks=" << tasks.size()
                  << " diagnostic=" << diagnostic << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "SU2_SEVEN_DEEP_SMALL_Z3"
              << " tasks=" << tasks.size()
              << " counterexamples=UNSAT result=PASS\n";
    return EXIT_SUCCESS;
}
