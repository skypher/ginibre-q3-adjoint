#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#define SU2_SEVEN_RESIDUAL_NO_MAIN
#include "verify_su2_seven_residual_z3.cpp"
#undef SU2_SEVEN_RESIDUAL_NO_MAIN

namespace {

std::size_t parse_size(const char* text, const char* name) {
    std::size_t parsed = 0U;
    const std::string input(text);
    const unsigned long long value = std::stoull(input, &parsed);
    if (
        parsed != input.size()
        || value > std::numeric_limits<std::size_t>::max()
    ) {
        throw std::runtime_error(std::string(name) + " is invalid");
    }
    return static_cast<std::size_t>(value);
}

unsigned adaptive_workers(std::size_t task_count) {
    const unsigned hardware =
        std::max(1U, std::thread::hardware_concurrency());
    double load_average = 0.0;
    const int load_status = getloadavg(&load_average, 1);
    const unsigned loaded = load_status == 1
        ? static_cast<unsigned>(
            std::max(0.0, std::ceil(load_average))
        )
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
    constexpr std::uint64_t reserve_kib = 8U * 1024U * 1024U;
    constexpr std::uint64_t kib_per_worker = 1536U * 1024U;
    const std::uint64_t usable_kib =
        available_kib > reserve_kib
            ? available_kib - reserve_kib
            : 0U;
    const unsigned memory_limit = static_cast<unsigned>(
        std::max<std::uint64_t>(1U, usable_kib / kib_per_worker)
    );
    return std::max(
        1U,
        std::min({
            cpu_limit,
            memory_limit,
            static_cast<unsigned>(task_count)
        })
    );
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 9) {
            throw std::runtime_error(
                "usage: TASK_INDEX ORBIT PARITY SELECTED "
                "POSITION KIND BEGIN_CELL END_CELL"
            );
        }
        const std::size_t task_index =
            parse_size(argv[1], "task index");
        const std::size_t orbit =
            parse_size(argv[2], "orbit");
        const std::size_t parity =
            parse_size(argv[3], "parity");
        const std::size_t selected =
            parse_size(argv[4], "selected orbit");
        const std::size_t position =
            parse_size(argv[5], "position");
        const std::size_t kind =
            parse_size(argv[6], "kind");
        const std::size_t begin =
            parse_size(argv[7], "begin cell");
        const std::size_t end =
            parse_size(argv[8], "end cell");
        if (
            task_index >= 308U
            || orbit >= residual_orbits.size()
            || parity > 1U
            || position >= 7U
            || kind > 2U
            || begin >= end
            || end > 16384U
        ) {
            throw std::runtime_error(
                "invalid task descriptor or chamber range"
            );
        }
        const std::vector<unsigned int> representatives =
            selected_orbit_masks(residual_orbits[orbit]);
        if (selected >= representatives.size()) {
            throw std::runtime_error(
                "selected orbit is outside the residual descriptor"
            );
        }
        const std::size_t count = end - begin;
        const unsigned workers = adaptive_workers(count);
        std::atomic<std::size_t> next(begin);
        std::atomic<std::size_t> complete(0U);
        std::atomic<bool> failed(false);
        std::mutex diagnostic_mutex;
        std::string diagnostic;

        std::cout
            << "SU2_SEVEN_SHALLOW_RANK_TWO_CELLS_Z3"
            << " task=" << task_index
            << " descriptor={orbit=" << orbit
            << " parity=" << parity
            << " selected=" << selected
            << " position=" << position
            << " kind=" << kind
            << " rank=2}"
            << " begin=" << begin
            << " end=" << end
            << " workers=" << workers
            << " start=1\n"
            << std::flush;

        std::vector<std::jthread> pool;
        pool.reserve(workers);
        for (unsigned worker = 0U; worker < workers; ++worker) {
            pool.emplace_back([&]() {
                while (!failed.load()) {
                    const std::size_t index = next.fetch_add(1U);
                    if (index >= end) {
                        return;
                    }
                    const int cap = static_cast<int>(index / 128U);
                    const int wall =
                        static_cast<int>((index / 16U) % 8U);
                    const int interval =
                        static_cast<int>(index % 16U);
                    const QueryResult result = verify_residual_query(
                        static_cast<int>(orbit),
                        static_cast<int>(parity),
                        static_cast<int>(selected),
                        cap,
                        wall,
                        interval,
                        -1,
                        0,
                        2,
                        static_cast<int>(position),
                        static_cast<int>(kind),
                        false,
                        -1,
                        -1,
                        -1
                    );
                    if (!result.passed) {
                        std::lock_guard<std::mutex> lock(
                            diagnostic_mutex
                        );
                        if (diagnostic.empty()) {
                            diagnostic =
                                "cell=" + std::to_string(index)
                                + " cap=" + std::to_string(cap)
                                + " wall=" + std::to_string(wall)
                                + " interval="
                                + std::to_string(interval)
                                + ' ' + result.diagnostic;
                        }
                        failed.store(true);
                        return;
                    }
                    const std::size_t done =
                        complete.fetch_add(1U) + 1U;
                    if (done % 1024U == 0U || done == count) {
                        std::cout
                            << "SU2_SEVEN_SHALLOW_RANK_TWO_CELLS_Z3"
                            << " task=" << task_index
                            << " progress=" << done << '/' << count
                            << " failed=" << (failed.load() ? 1 : 0)
                            << '\n'
                            << std::flush;
                    }
                }
            });
        }
        pool.clear();

        if (failed.load() || complete.load() != count) {
            std::cerr
                << "SU2_SEVEN_SHALLOW_RANK_TWO_CELLS_Z3 FAIL"
                << " task=" << task_index
                << " complete=" << complete.load()
                << " count=" << count
                << " diagnostic=" << diagnostic
                << '\n';
            return EXIT_FAILURE;
        }
        std::cout
            << "SU2_SEVEN_SHALLOW_RANK_TWO_CELLS_Z3"
            << " task=" << task_index
            << " begin=" << begin
            << " end=" << end
            << " counterexamples=UNSAT"
            << " result=PASS_EXACT_CERTIFICATE_SET\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_SEVEN_SHALLOW_RANK_TWO_CELLS_Z3 ERROR "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
