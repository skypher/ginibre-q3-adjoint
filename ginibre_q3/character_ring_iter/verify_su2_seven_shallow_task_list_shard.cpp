#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

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

bool run_task(
    std::size_t position,
    const std::filesystem::path& log_directory
) {
    const std::filesystem::path log_path =
        log_directory
        / ("small_task_" + std::to_string(position) + ".log");
    const pid_t child = fork();
    if (child < 0) {
        throw std::runtime_error("fork failed");
    }
    if (child == 0) {
        const int descriptor = open(
            log_path.c_str(),
            O_WRONLY | O_CREAT | O_TRUNC,
            static_cast<mode_t>(0644)
        );
        if (descriptor < 0) {
            _exit(125);
        }
        if (
            dup2(descriptor, STDOUT_FILENO) < 0
            || dup2(descriptor, STDERR_FILENO) < 0
        ) {
            _exit(125);
        }
        close(descriptor);
        const std::string position_text = std::to_string(position);
        execl(
            "./verify_su2_seven_shallow_z3",
            "verify_su2_seven_shallow_z3",
            "--small-task",
            position_text.c_str(),
            static_cast<char*>(nullptr)
        );
        _exit(errno == ENOENT ? 127 : 126);
    }

    int status = 0;
    while (waitpid(child, &status, 0) < 0) {
        if (errno != EINTR) {
            throw std::runtime_error("waitpid failed");
        }
    }
    return WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

std::size_t adaptive_thread_limit() {
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
    constexpr std::uint64_t kib_per_worker = 512U * 1024U;
    const std::uint64_t memory_limit = available_kib == 0U
        ? static_cast<std::uint64_t>(hardware)
        : std::max<std::uint64_t>(
            1U, available_kib / kib_per_worker
        );
    return std::max<std::size_t>(
        1U,
        std::min<std::uint64_t>(cpu_limit, memory_limit)
    );
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 4) {
            throw std::runtime_error(
                "usage: THREADS|auto LOG_DIRECTORY TASK_INDEX..."
            );
        }
        const std::size_t requested_threads =
            std::string(argv[1]) == "auto"
                ? adaptive_thread_limit()
                : parse_size(argv[1], "thread count");
        if (requested_threads == 0U) {
            throw std::runtime_error("thread count must be positive");
        }
        const std::filesystem::path log_directory(argv[2]);
        std::filesystem::create_directories(log_directory);

        std::vector<std::size_t> tasks;
        std::unordered_set<std::size_t> seen;
        tasks.reserve(static_cast<std::size_t>(argc - 3));
        for (int argument = 3; argument < argc; ++argument) {
            const std::size_t task =
                parse_size(argv[argument], "task index");
            if (task >= 308U || !seen.insert(task).second) {
                throw std::runtime_error(
                    "task indices must be distinct and below 308"
                );
            }
            tasks.push_back(task);
        }

        const std::size_t thread_count =
            std::min(requested_threads, tasks.size());
        std::atomic<std::size_t> next(0U);
        std::atomic<std::size_t> completed(0U);
        std::atomic<std::size_t> failed(0U);
        std::vector<std::thread> workers;
        workers.reserve(thread_count);
        for (std::size_t worker = 0U; worker < thread_count; ++worker) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t ordinal = next.fetch_add(1U);
                    if (ordinal >= tasks.size()) {
                        return;
                    }
                    try {
                        if (!run_task(tasks[ordinal], log_directory)) {
                            ++failed;
                        }
                    } catch (...) {
                        ++failed;
                    }
                    const std::size_t done = ++completed;
                    if (done % 25U == 0U || done == tasks.size()) {
                        std::cerr
                            << "SU2_SEVEN_SHALLOW_TASK_LIST_SHARD"
                            << " completed=" << done
                            << '/' << tasks.size()
                            << " failed=" << failed.load()
                            << std::endl;
                    }
                }
            });
        }
        for (std::thread& worker : workers) {
            worker.join();
        }

        const bool passed = failed.load() == 0U;
        std::cout
            << "SU2_SEVEN_SHALLOW_TASK_LIST_SHARD"
            << " threads=" << thread_count
            << " attempted=" << completed.load()
            << " failed=" << failed.load()
            << " result="
            << (passed ? "PASS_EXACT_CERTIFICATE_SET" : "INCOMPLETE")
            << '\n';
        return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_SEVEN_SHALLOW_TASK_LIST_SHARD FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
