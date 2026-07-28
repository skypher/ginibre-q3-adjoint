#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#define SU2_SEVEN_RESIDUAL_NO_MAIN
#include "verify_su2_seven_residual_z3.cpp"
#undef SU2_SEVEN_RESIDUAL_NO_MAIN

namespace {

int parse_index(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value < 0 || value > 16384) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
}

unsigned worker_count(int task_count) {
    const char* configured = std::getenv("Q3_MAX_THREADS");
    if (configured == nullptr || configured[0] == '\0') {
        return 1U;
    }
    char* end = nullptr;
    const long value = std::strtol(configured, &end, 10);
    if (end == configured || *end != '\0' || value <= 0) {
        throw std::runtime_error("invalid Q3_MAX_THREADS");
    }
    return std::min<unsigned>(
        static_cast<unsigned>(value), static_cast<unsigned>(task_count)
    );
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error("usage: BEGIN END");
        }
        const int begin = parse_index(argv[1], "begin");
        const int end = parse_index(argv[2], "end");
        if (end <= begin) {
            throw std::runtime_error("END must exceed BEGIN");
        }
        const int count = end - begin;
        const unsigned workers = worker_count(count);
        std::atomic<int> next{begin};
        std::atomic<int> complete{0};
        std::atomic<bool> failed{false};
        std::mutex diagnostic_mutex;
        std::string diagnostic;
        std::cout << "SU2_RANK_TWO_CHAMBERS_Z3 begin=" << begin
                  << " end=" << end << " workers=" << workers
                  << " start=1\n" << std::flush;
        std::vector<std::jthread> pool;
        pool.reserve(workers);
        for (unsigned worker = 0; worker < workers; ++worker) {
            pool.emplace_back([&]() {
                while (!failed.load()) {
                    const int index = next.fetch_add(1);
                    if (index >= end) {
                        return;
                    }
                    const int cap = index / (8 * 16);
                    const int wall = (index / 16) % 8;
                    const int interval = index % 16;
                    const QueryResult result = verify_residual_query(
                        1, 1, 1, cap, wall, interval,
                        -1, 0, 2, 1, 1, false, -1, -1, -1
                    );
                    if (!result.passed) {
                        std::lock_guard<std::mutex> lock(diagnostic_mutex);
                        if (diagnostic.empty()) {
                            diagnostic = "cap=" + std::to_string(cap)
                                + " wall=" + std::to_string(wall)
                                + " interval=" + std::to_string(interval)
                                + " " + result.diagnostic;
                        }
                        failed.store(true);
                        return;
                    }
                    complete.fetch_add(1);
                }
            });
        }
        pool.clear();
        if (failed.load() || complete.load() != count) {
            std::cerr << "SU2_RANK_TWO_CHAMBERS_Z3 FAIL complete="
                      << complete.load() << " count=" << count
                      << " diagnostic=" << diagnostic << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "SU2_RANK_TWO_CHAMBERS_Z3 begin=" << begin
                  << " end=" << end
                  << " counterexamples=UNSAT result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_RANK_TWO_CHAMBERS_Z3 ERROR "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
