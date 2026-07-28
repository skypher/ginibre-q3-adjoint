#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
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

std::vector<std::uint64_t> read_masks(
    const std::filesystem::path& path,
    bool subset
) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open mask list");
    }
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("mask list has no header");
    }
    std::size_t expected_size = 302U;
    if (subset) {
        const std::string subset_prefix =
            "SU2_K4_INTERMEDIATE_MASK_SUBSET count=";
        if (!line.starts_with(subset_prefix)) {
            throw std::runtime_error("mask subset has no subset header");
        }
        expected_size = parse_size(
            line.substr(subset_prefix.size()).c_str(),
            "mask subset count"
        );
        if (expected_size == 0U) {
            throw std::runtime_error("mask subset must be nonempty");
        }
    } else if (
        line
        != "SU2_K4_INTERMEDIATE_MASKS hinges=50 masks=302 "
           "minimum_active=7 maximum_active=44 "
           "result=PASS_EXACT_CENSUS"
    ) {
        throw std::runtime_error("mask list has no exact census header");
    }
    const std::string prefix = "SU2_K4_INTERMEDIATE_MASK value=";
    std::vector<std::uint64_t> masks;
    while (std::getline(input, line)) {
        if (!line.starts_with(prefix)) {
            throw std::runtime_error("malformed mask-list row");
        }
        std::size_t parsed = 0U;
        const std::string value_text = line.substr(prefix.size());
        const unsigned long long value =
            std::stoull(value_text, &parsed);
        if (
            parsed != value_text.size()
            || value > std::numeric_limits<std::uint64_t>::max()
        ) {
            throw std::runtime_error("invalid mask-list value");
        }
        masks.push_back(static_cast<std::uint64_t>(value));
    }
    const std::set<std::uint64_t> unique(masks.begin(), masks.end());
    if (
        masks.size() != expected_size
        || unique.size() != masks.size()
        || !std::is_sorted(masks.begin(), masks.end())
    ) {
        throw std::runtime_error(
            "mask list count, uniqueness, or order is invalid"
        );
    }
    return masks;
}

bool run_mask(
    std::uint64_t mask,
    const std::filesystem::path& log_directory
) {
    const std::filesystem::path log_path =
        log_directory / ("k4_mask_" + std::to_string(mask) + ".log");
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
        const std::string mask_text = std::to_string(mask);
        execl(
            "./prove_su2_k4_intermediate",
            "prove_su2_k4_intermediate",
            "--mask",
            mask_text.c_str(),
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
            argc != 4
            && !(
                argc == 5
                && (
                    std::string(argv[4]) == "--reverse"
                    || std::string(argv[4]) == "--subset"
                )
            )
        ) {
            throw std::runtime_error(
                "usage: MASK_LIST THREADS LOG_DIRECTORY "
                "[--reverse|--subset]"
            );
        }
        const bool reverse =
            argc == 5 && std::string(argv[4]) == "--reverse";
        const bool subset =
            argc == 5 && std::string(argv[4]) == "--subset";
        std::vector<std::uint64_t> masks = read_masks(argv[1], subset);
        if (reverse) {
            std::reverse(masks.begin(), masks.end());
        }
        const std::size_t requested_threads =
            parse_size(argv[2], "thread count");
        if (requested_threads == 0U) {
            throw std::runtime_error("thread count must be positive");
        }
        const std::size_t thread_count =
            std::min(requested_threads, masks.size());
        const std::filesystem::path log_directory(argv[3]);
        std::filesystem::create_directories(log_directory);

        std::atomic<std::size_t> next(0U);
        std::atomic<std::size_t> completed(0U);
        std::atomic<std::size_t> failed(0U);
        std::vector<std::thread> workers;
        workers.reserve(thread_count);
        for (std::size_t worker = 0U;
             worker < thread_count;
             ++worker) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t position = next.fetch_add(1U);
                    if (position >= masks.size()) {
                        return;
                    }
                    try {
                        if (!run_mask(
                                masks[position],
                                log_directory
                            )) {
                            ++failed;
                        }
                    } catch (...) {
                        ++failed;
                    }
                    const std::size_t count = ++completed;
                    if (
                        count % 25U == 0U
                        || count == masks.size()
                    ) {
                        std::cerr
                            << "SU2_K4_INTERMEDIATE_SHARD"
                            << " completed=" << count
                            << '/' << masks.size()
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
            << "SU2_K4_INTERMEDIATE_SHARD"
            << " masks=" << masks.size()
            << " threads=" << thread_count
            << " attempted=" << completed.load()
            << " failed=" << failed.load()
            << " result="
            << (passed ? "PASS_EXACT_CERTIFICATE_SET" : "INCOMPLETE")
            << std::endl;
        return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_K4_INTERMEDIATE_SHARD FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
