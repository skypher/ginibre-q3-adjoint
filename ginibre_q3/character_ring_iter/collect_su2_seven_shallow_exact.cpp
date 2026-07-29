#include <array>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/wait.h>

namespace {

namespace fs = std::filesystem;

struct Task {
    int index = -1;
    int orbit = -1;
    int parity = -1;
    int selected = -1;
    int position = -1;
    int kind = -1;
    int rank = -1;
    int pattern = -1;
    int wall = -1;
    int interval = -1;
    int local = -1;
};

std::string read_file(const fs::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(
            "cannot open receipt: " + path.string()
        );
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    if (input.bad()) {
        throw std::runtime_error(
            "cannot read receipt: " + path.string()
        );
    }
    return buffer.str();
}

bool has_bad_result(const std::string& text) {
    static const std::vector<std::string> bad = {
        "counterexamples=SAT",
        "result=FAIL",
        "result=UNKNOWN",
        " FAILURE:",
        " ERROR:"
    };
    for (const auto& token : bad) {
        if (text.find(token) != std::string::npos) {
            return true;
        }
    }
    return false;
}

std::vector<Task> read_tasks(const fs::path& verifier) {
    if (
        !fs::exists(verifier)
        || !fs::is_regular_file(verifier)
        || !std::regex_match(
            verifier.string(),
            std::regex(R"(^[A-Za-z0-9_./-]+$)")
        )
    ) {
        throw std::runtime_error(
            "invalid verifier path: " + verifier.string()
        );
    }
    const std::string command =
        verifier.string() + " --list-small";
    FILE* pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) {
        throw std::runtime_error("cannot execute verifier task list");
    }
    std::string output;
    std::array<char, 4096> buffer{};
    while (std::fgets(
               buffer.data(),
               static_cast<int>(buffer.size()),
               pipe
           ) != nullptr) {
        output += buffer.data();
    }
    const int status = pclose(pipe);
    if (
        status == -1
        || !WIFEXITED(status)
        || WEXITSTATUS(status) != 0
    ) {
        throw std::runtime_error("verifier task-list command failed");
    }
    std::istringstream input(output);
    const std::regex row(
        R"(^([0-9]+) orbit=([0-9]+) parity=([0-9]+) selected=([0-9]+) position=([0-9]+) kind=([0-9]+) rank=([0-9]+) pattern=(-?[0-9]+) wall=(-?[0-9]+) interval=(-?[0-9]+) local=(-?[0-9]+)$)"
    );
    std::vector<Task> tasks;
    std::string line;
    std::smatch match;
    while (std::getline(input, line)) {
        if (!std::regex_match(line, match, row)) {
            throw std::runtime_error(
                "malformed task-list row: " + line
            );
        }
        Task task;
        task.index = std::stoi(match[1].str());
        task.orbit = std::stoi(match[2].str());
        task.parity = std::stoi(match[3].str());
        task.selected = std::stoi(match[4].str());
        task.position = std::stoi(match[5].str());
        task.kind = std::stoi(match[6].str());
        task.rank = std::stoi(match[7].str());
        task.pattern = std::stoi(match[8].str());
        task.wall = std::stoi(match[9].str());
        task.interval = std::stoi(match[10].str());
        task.local = std::stoi(match[11].str());
        tasks.push_back(task);
    }
    if (tasks.size() != 308U) {
        throw std::runtime_error(
            "task list must contain exactly 308 rows"
        );
    }
    int rank_one = 0;
    int rank_two = 0;
    for (std::size_t index = 0; index < tasks.size(); ++index) {
        if (tasks[index].index != static_cast<int>(index)) {
            throw std::runtime_error(
                "task indices are not exactly 0,...,307"
            );
        }
        if (tasks[index].rank == 1) {
            ++rank_one;
        } else if (tasks[index].rank == 2) {
            ++rank_two;
        } else {
            throw std::runtime_error("unexpected task rank");
        }
    }
    if (rank_one != 168 || rank_two != 140) {
        throw std::runtime_error(
            "task-rank split must be exactly 168+140"
        );
    }
    return tasks;
}

bool direct_pass(const std::string& text) {
    return !has_bad_result(text)
        && text.find(
            "SU2_SEVEN_SHALLOW_SMALL_TASK_Z3 tasks=1 "
            "counterexamples=UNSAT result=PASS"
        ) != std::string::npos;
}

bool chamber_pass(
    const std::string& text,
    const Task& task
) {
    if (task.rank != 2 || has_bad_result(text)) {
        return false;
    }
    std::ostringstream descriptor;
    descriptor
        << "SU2_SEVEN_SHALLOW_RANK_TWO_CELLS_Z3 task="
        << task.index
        << " descriptor={orbit=" << task.orbit
        << " parity=" << task.parity
        << " selected=" << task.selected
        << " position=" << task.position
        << " kind=" << task.kind
        << " rank=2} begin=0 end=16384";
    std::ostringstream terminal;
    terminal
        << "SU2_SEVEN_SHALLOW_RANK_TWO_CELLS_Z3 task="
        << task.index
        << " begin=0 end=16384 counterexamples=UNSAT "
        << "result=PASS_EXACT_CERTIFICATE_SET";
    return text.find(descriptor.str()) != std::string::npos
        && text.find(terminal.str()) != std::string::npos;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 4) {
            throw std::runtime_error(
                "usage: collect_su2_seven_shallow_exact "
                "VERIFIER RANK_ONE_LOG "
                "[--export NEW_DIRECTORY] RECEIPT_ROOT..."
            );
        }
        int root_begin = 3;
        fs::path export_root;
        if (std::string(argv[3]) == "--export") {
            if (argc < 6) {
                throw std::runtime_error(
                    "--export requires a new directory and "
                    "at least one receipt root"
                );
            }
            export_root = argv[4];
            root_begin = 5;
            if (fs::exists(export_root)) {
                throw std::runtime_error(
                    "export path already exists: "
                    + export_root.string()
                );
            }
        }
        const std::vector<Task> tasks = read_tasks(argv[1]);
        const std::string rank_one_text = read_file(argv[2]);
        if (
            has_bad_result(rank_one_text)
            || rank_one_text.find(
                "SU2_SEVEN_SHALLOW_SMALL_SWITCH_Z3 tasks=168 "
                "counterexamples=UNSAT result=PASS"
            ) == std::string::npos
        ) {
            throw std::runtime_error(
                "rank-one aggregate receipt is not terminal PASS"
            );
        }

        std::vector<bool> covered(tasks.size(), false);
        std::set<int> direct_tasks;
        std::set<int> chamber_tasks;
        std::vector<fs::path> direct_paths(tasks.size());
        std::vector<fs::path> chamber_paths(tasks.size());
        for (const Task& task : tasks) {
            if (task.rank == 1) {
                covered[static_cast<std::size_t>(task.index)] = true;
            }
        }

        const std::regex direct_name(
            R"(^small_task_([0-9]+)\.log$)"
        );
        const std::regex chamber_name(
            R"(^rank_two_task([0-9]+)_cells\.log$)"
        );
        for (int argument = root_begin; argument < argc; ++argument) {
            const fs::path root = argv[argument];
            if (!fs::exists(root) || !fs::is_directory(root)) {
                throw std::runtime_error(
                    "receipt root is not a directory: "
                    + root.string()
                );
            }
            for (const auto& entry :
                 fs::recursive_directory_iterator(root)) {
                if (!entry.is_regular_file()) {
                    continue;
                }
                const std::string name =
                    entry.path().filename().string();
                std::smatch match;
                if (std::regex_match(name, match, direct_name)) {
                    const int index = std::stoi(match[1].str());
                    if (index < 0
                        || index >= static_cast<int>(tasks.size())) {
                        throw std::runtime_error(
                            "out-of-range direct receipt: "
                            + entry.path().string()
                        );
                    }
                    const std::string text = read_file(entry.path());
                    if (direct_pass(text)) {
                        direct_tasks.insert(index);
                        covered[static_cast<std::size_t>(index)] = true;
                        fs::path& selected =
                            direct_paths[
                                static_cast<std::size_t>(index)
                            ];
                        if (
                            selected.empty()
                            || entry.path().string()
                                < selected.string()
                        ) {
                            selected = entry.path();
                        }
                    }
                    continue;
                }
                if (std::regex_match(name, match, chamber_name)) {
                    const int index = std::stoi(match[1].str());
                    if (index < 0
                        || index >= static_cast<int>(tasks.size())) {
                        throw std::runtime_error(
                            "out-of-range chamber receipt: "
                            + entry.path().string()
                        );
                    }
                    const std::string text = read_file(entry.path());
                    if (chamber_pass(
                            text,
                            tasks[static_cast<std::size_t>(index)]
                        )) {
                        chamber_tasks.insert(index);
                        covered[static_cast<std::size_t>(index)] = true;
                        fs::path& selected =
                            chamber_paths[
                                static_cast<std::size_t>(index)
                            ];
                        if (
                            selected.empty()
                            || entry.path().string()
                                < selected.string()
                        ) {
                            selected = entry.path();
                        }
                    }
                }
            }
        }

        std::vector<int> missing;
        int covered_rank_one = 0;
        int covered_rank_two = 0;
        for (const Task& task : tasks) {
            if (covered[static_cast<std::size_t>(task.index)]) {
                if (task.rank == 1) {
                    ++covered_rank_one;
                } else {
                    ++covered_rank_two;
                }
            } else {
                missing.push_back(task.index);
            }
        }
        std::size_t overlap = 0;
        for (const int index : chamber_tasks) {
            if (direct_tasks.contains(index)) {
                ++overlap;
            }
        }

        std::cout
            << "SU2_SEVEN_SHALLOW_EXACT_COLLECTOR"
            << " tasks=308"
            << " rank_one=" << covered_rank_one << "/168"
            << " rank_two=" << covered_rank_two << "/140"
            << " direct_receipts=" << direct_tasks.size()
            << " chamber_receipts=" << chamber_tasks.size()
            << " method_overlap=" << overlap
            << " missing=" << missing.size();
        if (!missing.empty()) {
            std::cout << " missing_indices={";
            for (std::size_t index = 0;
                 index < missing.size(); ++index) {
                if (index != 0U) {
                    std::cout << ',';
                }
                std::cout << missing[index];
            }
            std::cout << '}';
        }
        std::cout
            << " result="
            << (missing.empty() ? "PASS_EXACT_UNION" : "INCOMPLETE")
            << '\n';
        if (!export_root.empty()) {
            if (!missing.empty()) {
                throw std::runtime_error(
                    "cannot export an incomplete receipt set"
                );
            }
            fs::create_directories(export_root);
            fs::copy_file(
                argv[2],
                export_root / "rank_one_aggregate.log"
            );
            std::ofstream index(export_root / "receipt_index.txt");
            if (!index) {
                throw std::runtime_error(
                    "cannot create export receipt index"
                );
            }
            index
                << "SU2_SEVEN_SHALLOW_EXACT_RECEIPT_INDEX"
                << " tasks=308 rank_one=168 rank_two=140\n";
            for (const Task& task : tasks) {
                index
                    << task.index
                    << " orbit=" << task.orbit
                    << " parity=" << task.parity
                    << " selected=" << task.selected
                    << " position=" << task.position
                    << " kind=" << task.kind
                    << " rank=" << task.rank
                    << " pattern=" << task.pattern
                    << " wall=" << task.wall
                    << " interval=" << task.interval
                    << " local=" << task.local;
                if (task.rank == 1) {
                    index
                        << " method=rank_one_aggregate"
                        << " receipt=rank_one_aggregate.log\n";
                    continue;
                }
                const std::size_t task_index =
                    static_cast<std::size_t>(task.index);
                const bool use_chamber =
                    !chamber_paths[task_index].empty();
                const fs::path& source = use_chamber
                    ? chamber_paths[task_index]
                    : direct_paths[task_index];
                if (source.empty()) {
                    throw std::runtime_error(
                        "internal error: covered task lacks receipt"
                    );
                }
                const std::string filename = use_chamber
                    ? (
                        "rank_two_task"
                        + std::to_string(task.index)
                        + "_cells.log"
                    )
                    : (
                        "small_task_"
                        + std::to_string(task.index)
                        + ".log"
                    );
                fs::copy_file(source, export_root / filename);
                index
                    << " method="
                    << (use_chamber ? "chamber" : "direct")
                    << " receipt=" << filename << '\n';
            }
            index.flush();
            if (!index) {
                throw std::runtime_error(
                    "cannot finish export receipt index"
                );
            }
            std::cout
                << "SU2_SEVEN_SHALLOW_EXACT_EXPORT"
                << " directory=" << export_root.string()
                << " files=142 result=PASS\n";
        }
        return missing.empty() ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_SEVEN_SHALLOW_EXACT_COLLECTOR FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
