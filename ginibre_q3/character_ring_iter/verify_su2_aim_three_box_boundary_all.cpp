#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>

extern char** environ;

namespace {

struct Case {
    std::array<int, 5> coordinates{};
};

struct Result {
    int exit_code = EXIT_FAILURE;
    std::string output;
};

std::size_t available_memory_threads() {
    std::ifstream input("/proc/meminfo");
    std::string name;
    std::uint64_t value = 0U;
    std::string unit;
    while (input >> name >> value >> unit) {
        if (name == "MemAvailable:") {
            constexpr std::uint64_t kibibytes_per_worker
                = 512U * 1024U;
            return static_cast<std::size_t>(
                std::max<std::uint64_t>(
                    1U,
                    value / kibibytes_per_worker));
        }
    }
    return 1U;
}

double current_load() {
    std::ifstream input("/proc/loadavg");
    double load = 0.0;
    input >> load;
    return load;
}

std::size_t worker_count(std::size_t cases, double load) {
    const unsigned detected = std::thread::hardware_concurrency();
    const std::size_t processors
        = detected == 0U ? 1U : static_cast<std::size_t>(detected);
    const std::size_t occupied = static_cast<std::size_t>(
        std::max(0.0, std::ceil(load)));
    const std::size_t load_limit
        = processors > occupied ? processors - occupied : 1U;
    std::size_t workers = std::min(
        {cases, load_limit, available_memory_threads()});
    if (const char* cap_text = std::getenv("Q3_MAX_THREADS");
        cap_text != nullptr) {
        const std::string cap_value(cap_text);
        std::size_t consumed = 0U;
        const unsigned long long parsed
            = std::stoull(cap_value, &consumed, 10);
        if (consumed != cap_value.size() || parsed == 0U) {
            throw std::invalid_argument(
                "Q3_MAX_THREADS must be a positive integer");
        }
        const auto maximum
            = static_cast<unsigned long long>(
                std::numeric_limits<std::size_t>::max());
        workers = std::min(
            workers,
            static_cast<std::size_t>(std::min(parsed, maximum)));
    }
    return std::max<std::size_t>(1U, workers);
}

Result run_case(const std::string& verifier, const Case& test_case) {
    std::array<std::string, 5> arguments;
    std::array<char*, 7> argv{};
    argv[0] = const_cast<char*>(verifier.c_str());
    for (std::size_t index = 0; index < arguments.size(); ++index) {
        arguments[index] = std::to_string(
            test_case.coordinates[index]);
        argv[index + 1U] = arguments[index].data();
    }
    argv[6] = nullptr;

    int descriptors[2]{-1, -1};
    if (pipe(descriptors) != 0) {
        throw std::runtime_error("pipe failed");
    }
    posix_spawn_file_actions_t actions;
    if (posix_spawn_file_actions_init(&actions) != 0) {
        close(descriptors[0]);
        close(descriptors[1]);
        throw std::runtime_error("spawn actions initialization failed");
    }
    posix_spawn_file_actions_adddup2(
        &actions,
        descriptors[1],
        STDOUT_FILENO);
    posix_spawn_file_actions_adddup2(
        &actions,
        descriptors[1],
        STDERR_FILENO);
    posix_spawn_file_actions_addclose(&actions, descriptors[0]);
    posix_spawn_file_actions_addclose(&actions, descriptors[1]);

    pid_t process = -1;
    const int spawn_result = posix_spawn(
        &process,
        verifier.c_str(),
        &actions,
        nullptr,
        argv.data(),
        environ);
    posix_spawn_file_actions_destroy(&actions);
    close(descriptors[1]);
    if (spawn_result != 0) {
        close(descriptors[0]);
        throw std::runtime_error("posix_spawn failed");
    }

    std::string output;
    std::array<char, 4096> buffer{};
    while (true) {
        const ssize_t count = read(
            descriptors[0],
            buffer.data(),
            buffer.size());
        if (count == 0) {
            break;
        }
        if (count < 0) {
            close(descriptors[0]);
            throw std::runtime_error("pipe read failed");
        }
        output.append(
            buffer.data(),
            static_cast<std::size_t>(count));
    }
    close(descriptors[0]);

    int status = 0;
    if (waitpid(process, &status, 0) < 0) {
        throw std::runtime_error("waitpid failed");
    }
    const int exit_code = WIFEXITED(status)
        ? WEXITSTATUS(status)
        : EXIT_FAILURE;
    return {exit_code, output};
}

std::uint64_t field(
    const std::string& output,
    const std::string& name) {
    const std::size_t begin = output.find(name);
    if (begin == std::string::npos) {
        throw std::runtime_error("missing output field " + name);
    }
    const std::size_t value_begin = begin + name.size();
    std::size_t value_end = value_begin;
    while (value_end < output.size()
           && output[value_end] >= '0'
           && output[value_end] <= '9') {
        ++value_end;
    }
    return std::stoull(
        output.substr(value_begin, value_end - value_begin));
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: verify_su2_aim_three_box_boundary_all "
                "[case_verifier] [certificate]");
        }
        const std::string verifier = argc >= 2
            ? std::string(argv[1])
            : "./verify_su2_aim_three_box_boundary_cells_z3";
        std::vector<Case> cases;
        for (int support_band = 0; support_band < 2; ++support_band) {
            for (int width_order = 0; width_order < 3; ++width_order) {
                for (int shift_band = 0; shift_band < 2; ++shift_band) {
                    for (int pair_order = 0;
                         pair_order < 2;
                         ++pair_order) {
                        for (int parity = 0; parity < 2; ++parity) {
                            cases.push_back(
                                {{support_band,
                                  width_order,
                                  shift_band,
                                  pair_order,
                                  parity}});
                        }
                    }
                }
            }
        }

        const double load = current_load();
        const std::size_t workers = worker_count(cases.size(), load);
        std::vector<Result> results(cases.size());
        std::atomic<std::size_t> next{0U};
        std::vector<std::thread> threads;
        threads.reserve(workers);
        for (std::size_t worker = 0; worker < workers; ++worker) {
            threads.emplace_back([&]() {
                while (true) {
                    const std::size_t index = next.fetch_add(1U);
                    if (index >= cases.size()) {
                        return;
                    }
                    results[index] = run_case(verifier, cases[index]);
                }
            });
        }
        for (std::thread& thread : threads) {
            thread.join();
        }

        std::uint64_t total_cells = 0U;
        std::uint64_t total_handelman = 0U;
        std::uint64_t total_real = 0U;
        std::uint64_t total_nia = 0U;
        std::uint64_t total_residue_queries = 0U;
        std::uint64_t maximum_scale = 1U;
        bool pass = true;
        std::ostringstream report;
        report << "SU2_AIM_THREE_BOX_BOUNDARY_ALL"
               << " cases=" << cases.size()
               << " workers=" << workers
               << " initial_load=" << load << '\n';
        for (std::size_t index = 0; index < cases.size(); ++index) {
            report << results[index].output;
            const bool success = results[index].exit_code == EXIT_SUCCESS
                && results[index].output.find("result=unsat")
                    != std::string::npos;
            pass = pass && success;
            if (success) {
                total_cells += field(results[index].output, "cells=");
                total_handelman += field(
                    results[index].output,
                    "handelman_cells=");
                total_real += field(
                    results[index].output,
                    "real_cells=");
                total_nia += field(
                    results[index].output,
                    "nia_cells=");
                total_residue_queries += field(
                    results[index].output,
                    "residue_queries=");
                maximum_scale = std::max(
                    maximum_scale,
                    field(results[index].output, "maximum_scale="));
            }
        }
        report << "SU2_AIM_THREE_BOX_BOUNDARY_ALL_SUMMARY"
               << " cells=" << total_cells
               << " handelman_cells=" << total_handelman
               << " real_cells=" << total_real
               << " nia_cells=" << total_nia
               << " residue_queries=" << total_residue_queries
               << " maximum_scale=" << maximum_scale
               << " result=" << (pass ? "PASS" : "FAIL") << '\n';
        const std::string report_text = report.str();
        std::cout << report_text;
        if (argc == 3) {
            std::ofstream certificate(argv[2]);
            certificate << report_text;
            if (!certificate) {
                throw std::runtime_error(
                    "could not write certificate");
            }
        }
        return pass ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
