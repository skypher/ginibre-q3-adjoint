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

}  // namespace

int main(int argc, char** argv) {
    try {
        if (
            argc != 5
            && !(
                argc == 6
                && std::string(argv[5]) == "--reverse"
            )
        ) {
            throw std::runtime_error(
                "usage: BEGIN_POSITION END_POSITION THREADS "
                "LOG_DIRECTORY [--reverse]"
            );
        }
        const std::size_t begin =
            parse_size(argv[1], "begin position");
        const std::size_t end =
            parse_size(argv[2], "end position");
        const std::size_t requested_threads =
            parse_size(argv[3], "thread count");
        if (
            begin >= end
            || end > 308U
            || requested_threads == 0U
        ) {
            throw std::runtime_error(
                "require 0<=begin<end<=308 and positive threads"
            );
        }
        const std::size_t thread_count =
            std::min(requested_threads, end - begin);
        const std::filesystem::path log_directory(argv[4]);
        std::filesystem::create_directories(log_directory);
        const bool reverse = argc == 6;

        std::atomic<std::size_t> next(begin);
        std::atomic<std::size_t> completed(0U);
        std::atomic<std::size_t> failed(0U);
        std::vector<std::thread> workers;
        workers.reserve(thread_count);
        for (std::size_t worker = 0U;
             worker < thread_count;
             ++worker) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t ordinal = next.fetch_add(1U);
                    if (ordinal >= end) {
                        return;
                    }
                    const std::size_t position = reverse
                        ? end - 1U - (ordinal - begin)
                        : ordinal;
                    try {
                        if (!run_task(position, log_directory)) {
                            ++failed;
                        }
                    } catch (...) {
                        ++failed;
                    }
                    const std::size_t count = ++completed;
                    if (count % 25U == 0U || count == end - begin) {
                        std::cerr
                            << "SU2_SEVEN_SHALLOW_TASK_SHARD"
                            << " completed=" << count
                            << '/' << end - begin
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
            << "SU2_SEVEN_SHALLOW_TASK_SHARD"
            << " begin=" << begin
            << " end=" << end
            << " threads=" << thread_count
            << " attempted=" << completed.load()
            << " failed=" << failed.load()
            << " result="
            << (passed ? "PASS_EXACT_CERTIFICATE_SET" : "INCOMPLETE")
            << std::endl;
        return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_SEVEN_SHALLOW_TASK_SHARD FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
