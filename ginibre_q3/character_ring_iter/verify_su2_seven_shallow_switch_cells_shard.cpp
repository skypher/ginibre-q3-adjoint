#include <atomic>
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
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

bool run_cell(
    std::size_t task,
    std::size_t cap,
    std::size_t cell,
    const std::filesystem::path& log_directory
) {
    const std::size_t wall = cell / 16U;
    const std::size_t interval = cell % 16U;
    const std::filesystem::path log_path =
        log_directory / ("cell_" + std::to_string(cell) + ".log");
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
        const std::string task_text = std::to_string(task);
        const std::string cap_text = std::to_string(cap);
        const std::string wall_text = std::to_string(wall);
        const std::string interval_text = std::to_string(interval);
        execl(
            "./verify_su2_seven_shallow_z3",
            "verify_su2_seven_shallow_z3",
            "--small-switch-cell",
            task_text.c_str(),
            cap_text.c_str(),
            wall_text.c_str(),
            interval_text.c_str(),
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

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 7) {
            throw std::runtime_error(
                "usage: TASK_INDEX CAP_MASK BEGIN_CELL END_CELL "
                "THREADS LOG_DIRECTORY"
            );
        }
        const std::size_t task = parse_size(argv[1], "task index");
        const std::size_t cap = parse_size(argv[2], "cap mask");
        const std::size_t begin = parse_size(argv[3], "begin cell");
        const std::size_t end = parse_size(argv[4], "end cell");
        const std::size_t requested_threads =
            parse_size(argv[5], "thread count");
        if (
            task >= 308U
            || cap >= 128U
            || begin >= end
            || end > 128U
            || requested_threads == 0U
        ) {
            throw std::runtime_error(
                "require task<308, cap<128, "
                "0<=begin<end<=128, and positive threads"
            );
        }
        const std::filesystem::path log_directory(argv[6]);
        std::filesystem::create_directories(log_directory);

        const std::size_t count = end - begin;
        const std::size_t thread_count =
            std::min(requested_threads, count);
        std::atomic<std::size_t> next(begin);
        std::atomic<std::size_t> completed(0U);
        std::atomic<std::size_t> failed(0U);
        std::vector<std::thread> workers;
        workers.reserve(thread_count);
        for (std::size_t worker = 0U; worker < thread_count; ++worker) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t cell = next.fetch_add(1U);
                    if (cell >= end) {
                        return;
                    }
                    try {
                        if (!run_cell(task, cap, cell, log_directory)) {
                            ++failed;
                        }
                    } catch (...) {
                        ++failed;
                    }
                    const std::size_t done = ++completed;
                    if (done % 16U == 0U || done == count) {
                        std::cerr
                            << "SU2_SEVEN_SHALLOW_SWITCH_CELL_SHARD"
                            << " task=" << task
                            << " cap=" << cap
                            << " completed=" << done
                            << '/' << count
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
            << "SU2_SEVEN_SHALLOW_SWITCH_CELL_SHARD"
            << " task=" << task
            << " cap=" << cap
            << " begin=" << begin
            << " end=" << end
            << " threads=" << thread_count
            << " attempted=" << completed.load()
            << " failed=" << failed.load()
            << " result="
            << (passed ? "PASS_EXACT_CERTIFICATE_SET" : "INCOMPLETE")
            << '\n';
        return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_SEVEN_SHALLOW_SWITCH_CELL_SHARD FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
