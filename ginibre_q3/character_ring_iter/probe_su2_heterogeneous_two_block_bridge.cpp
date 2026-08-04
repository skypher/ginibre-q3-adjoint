#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

int parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0U;
    const long parsed = std::stol(value, &consumed);
    if (consumed != value.size() || parsed <= 0L || parsed > 8L) {
        throw std::invalid_argument(name + " must be an integer in [1,8]");
    }
    return static_cast<int>(parsed);
}

bool edge(int left, int right, int label) {
    return std::abs(left - label) <= right && right <= left + label;
}

void extend_paths(
    const std::vector<int>& labels,
    std::size_t position,
    std::vector<int>& path,
    std::vector<std::vector<int>>& paths
) {
    if (position == labels.size()) {
        paths.push_back(path);
        return;
    }
    const int source = path.back();
    const int label = labels[position];
    for (int target = std::abs(source - label);
         target <= source + label;
         ++target) {
        path.push_back(target);
        extend_paths(labels, position + 1U, path, paths);
        path.pop_back();
    }
}

std::vector<std::vector<int>> paths_from(
    int start,
    const std::vector<int>& labels
) {
    std::vector<std::vector<int>> paths;
    std::vector<int> path;
    path.reserve(labels.size() + 1U);
    path.push_back(start);
    extend_paths(labels, 0U, path, paths);
    return paths;
}

bool two_step_support(
    int left,
    int right,
    int first_label,
    int second_label
) {
    for (int middle = std::abs(left - first_label);
         middle <= left + first_label;
         ++middle) {
        if (edge(middle, right, second_label)) {
            return true;
        }
    }
    return false;
}

bool has_one_block_bridge(
    const std::vector<int>& long_path,
    const std::vector<int>& short_path,
    const std::vector<int>& labels
) {
    for (std::size_t time = 0U; time < labels.size(); ++time) {
        if (edge(
                short_path[time],
                long_path[time + 1U],
                labels[time]
            )
            && edge(
                long_path[time],
                short_path[time + 1U],
                labels[time]
            )) {
            return true;
        }
    }
    return false;
}

bool has_two_block_bridge(
    const std::vector<int>& long_path,
    const std::vector<int>& short_path,
    const std::vector<int>& labels
) {
    for (std::size_t time = 0U; time + 1U < labels.size(); ++time) {
        if (two_step_support(
                short_path[time],
                long_path[time + 2U],
                labels[time],
                labels[time + 1U]
            )
            && two_step_support(
                long_path[time],
                short_path[time + 2U],
                labels[time],
                labels[time + 1U]
            )) {
            return true;
        }
    }
    return false;
}

std::string show(const std::vector<int>& values) {
    std::string result{"("};
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            result.push_back(',');
        }
        result += std::to_string(values[index]);
    }
    result.push_back(')');
    return result;
}

struct Scan {
    unsigned long long domain_pairs = 0U;
    unsigned long long candidate_free_pairs = 0U;
    bool has_obstruction = false;
    std::string witness;
};

class CandidateFreeSearch {
public:
    explicit CandidateFreeSearch(const std::vector<int>& labels)
        : labels_(labels) {
        int total_label = 0;
        for (const int label : labels_) {
            total_label += label;
        }
        if (labels_.size() > 255U || total_label > 1023) {
            throw std::invalid_argument(
                "the exact dynamic search requires at most 255 blocks "
                "and total label at most 1023"
            );
        }
    }

    bool find(std::vector<int>& long_path, std::vector<int>& short_path) {
        int total_label = 0;
        for (const int label : labels_) {
            total_label += label;
        }
        for (int initial = 1; initial <= total_label; ++initial) {
            long_path.assign(1U, initial);
            short_path.assign(1U, 0);
            if (extend(long_path, short_path)) {
                return true;
            }
        }
        long_path.clear();
        short_path.clear();
        return false;
    }

    std::uint64_t explored_states() const {
        return explored_states_;
    }

private:
    std::uint64_t state_key(
        std::size_t steps,
        const std::vector<int>& long_path,
        const std::vector<int>& short_path
    ) const {
        const auto pack = [](int value) {
            return static_cast<std::uint64_t>(value);
        };
        const std::size_t current = steps;
        const std::size_t previous = steps - 1U;
        return (static_cast<std::uint64_t>(steps) << 48U)
            | (pack(long_path[previous]) << 36U)
            | (pack(long_path[current]) << 24U)
            | (pack(short_path[previous]) << 12U)
            | pack(short_path[current]);
    }

    bool extend(
        std::vector<int>& long_path,
        std::vector<int>& short_path
    ) {
        const std::size_t steps = long_path.size() - 1U;
        if (steps == labels_.size()) {
            return long_path.back() == 0 && short_path.back() > 0;
        }
        const int label = labels_[steps];
        const int long_source = long_path.back();
        const int short_source = short_path.back();
        for (int long_target = std::abs(long_source - label);
             long_target <= long_source + label;
             ++long_target) {
            for (int short_target = std::abs(short_source - label);
                 short_target <= short_source + label;
                 ++short_target) {
                if (edge(short_source, long_target, label)
                    && edge(long_source, short_target, label)) {
                    continue;
                }
                if (steps >= 1U
                    && two_step_support(
                        short_path[steps - 1U],
                        long_target,
                        labels_[steps - 1U],
                        label
                    )
                    && two_step_support(
                        long_path[steps - 1U],
                        short_target,
                        labels_[steps - 1U],
                        label
                    )) {
                    continue;
                }
                long_path.push_back(long_target);
                short_path.push_back(short_target);
                const std::size_t next_steps = steps + 1U;
                bool known_dead = false;
                std::uint64_t key = 0U;
                if (next_steps >= 1U && next_steps < labels_.size()) {
                    key = state_key(next_steps, long_path, short_path);
                    known_dead = dead_states_.contains(key);
                }
                if (!known_dead) {
                    ++explored_states_;
                    if (extend(long_path, short_path)) {
                        return true;
                    }
                    if (next_steps >= 1U && next_steps < labels_.size()) {
                        dead_states_.insert(key);
                    }
                }
                long_path.pop_back();
                short_path.pop_back();
            }
        }
        return false;
    }

    const std::vector<int>& labels_;
    std::unordered_set<std::uint64_t> dead_states_;
    std::uint64_t explored_states_ = 0U;
};

struct MountainScan {
    std::uint64_t words = 0U;
    std::uint64_t explored_states = 0U;
    bool found_obstruction = false;
    std::vector<int> word;
    std::vector<int> long_path;
    std::vector<int> short_path;
};

void scan_sorted_roots(
    int maximum_label,
    int root_length,
    int minimum_label,
    std::vector<int>& root,
    MountainScan& scan
) {
    if (scan.found_obstruction) {
        return;
    }
    if (static_cast<int>(root.size()) == root_length) {
        std::vector<int> word = root;
        word.insert(word.end(), root.rbegin(), root.rend());
        CandidateFreeSearch search(word);
        std::vector<int> long_path;
        std::vector<int> short_path;
        const bool found = search.find(long_path, short_path);
        ++scan.words;
        scan.explored_states += search.explored_states();
        if (found) {
            scan.found_obstruction = true;
            scan.word = std::move(word);
            scan.long_path = std::move(long_path);
            scan.short_path = std::move(short_path);
        }
        return;
    }
    for (int label = minimum_label; label <= maximum_label; ++label) {
        root.push_back(label);
        scan_sorted_roots(
            maximum_label,
            root_length,
            label,
            root,
            scan
        );
        root.pop_back();
        if (scan.found_obstruction) {
            return;
        }
    }
}

Scan scan_order(const std::vector<int>& labels) {
    Scan result;
    int total_label = 0;
    for (const int label : labels) {
        total_label += label;
    }
    const auto short_paths = paths_from(0, labels);
    for (int initial = 1; initial <= total_label; ++initial) {
        const auto long_paths = paths_from(initial, labels);
        for (const auto& long_path : long_paths) {
            if (long_path.back() != 0) {
                continue;
            }
            for (const auto& short_path : short_paths) {
                if (short_path.back() == 0) {
                    continue;
                }
                ++result.domain_pairs;
                if (has_one_block_bridge(long_path, short_path, labels)
                    || has_two_block_bridge(long_path, short_path, labels)) {
                    continue;
                }
                ++result.candidate_free_pairs;
                if (!result.has_obstruction) {
                    result.has_obstruction = true;
                    result.witness = "word=" + show(labels)
                        + " long=" + show(long_path)
                        + " short=" + show(short_path);
                }
            }
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4 && std::string(argv[1]) == "--sorted-root-scan") {
            const int maximum_label = parse_positive(argv[2], "MAX_LABEL");
            const int maximum_root_length =
                parse_positive(argv[3], "MAX_ROOT_LENGTH");
            MountainScan scan;
            for (int root_length = 1;
                 root_length <= maximum_root_length;
                 ++root_length) {
                std::vector<int> root;
                root.reserve(static_cast<std::size_t>(root_length));
                scan_sorted_roots(
                    maximum_label,
                    root_length,
                    1,
                    root,
                    scan
                );
                if (scan.found_obstruction) {
                    break;
                }
            }
            std::cout << "SU2_HETEROGENEOUS_TWO_BLOCK_BRIDGE_MOUNTAIN_SCAN"
                      << " maximum_label=" << maximum_label
                      << " maximum_root_length=" << maximum_root_length
                      << " words=" << scan.words
                      << " explored_states=" << scan.explored_states
                      << " candidate_free="
                      << (scan.found_obstruction ? 1 : 0);
            if (scan.found_obstruction) {
                std::cout << " word=" << show(scan.word)
                          << " long=" << show(scan.long_path)
                          << " short=" << show(scan.short_path);
            }
            std::cout << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc >= 3 && std::string(argv[1]) == "--find-obstruction") {
            std::vector<int> word;
            word.reserve(static_cast<std::size_t>(argc - 2));
            for (int index = 2; index < argc; ++index) {
                word.push_back(parse_positive(argv[index], "path label"));
            }
            CandidateFreeSearch search(word);
            std::vector<int> long_path;
            std::vector<int> short_path;
            const bool found = search.find(long_path, short_path);
            std::cout << "SU2_HETEROGENEOUS_TWO_BLOCK_BRIDGE_DYNAMIC_SCAN"
                      << " word=" << show(word)
                      << " candidate_free=" << (found ? 1 : 0)
                      << " explored_states=" << search.explored_states();
            if (found) {
                std::cout << " long=" << show(long_path)
                          << " short=" << show(short_path);
            }
            std::cout << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc >= 3
            && (std::string(argv[1]) == "--root"
                || std::string(argv[1]) == "--paired-root"
                || std::string(argv[1]) == "--palindromic-root"
                || std::string(argv[1]) == "--word")) {
            std::vector<int> labels;
            labels.reserve(static_cast<std::size_t>(argc - 2));
            for (int index = 2; index < argc; ++index) {
                labels.push_back(parse_positive(argv[index], "path label"));
            }
            const bool root_mode = std::string(argv[1]) == "--root";
            const bool paired_root_mode =
                std::string(argv[1]) == "--paired-root";
            const bool palindromic_root_mode =
                std::string(argv[1]) == "--palindromic-root";
            std::vector<int> word = labels;
            if (root_mode) {
                word.insert(word.end(), labels.begin(), labels.end());
            } else if (paired_root_mode) {
                word.clear();
                word.reserve(2U * labels.size());
                for (const int label : labels) {
                    word.push_back(label);
                    word.push_back(label);
                }
            } else if (palindromic_root_mode) {
                word.insert(word.end(), labels.rbegin(), labels.rend());
            }
            const Scan scan = scan_order(word);
            std::cout << "SU2_HETEROGENEOUS_TWO_BLOCK_BRIDGE_"
                      << (root_mode
                              ? "ROOT_SCAN"
                              : (paired_root_mode
                                      ? "PAIRED_ROOT_SCAN"
                                      : (palindromic_root_mode
                                              ? "PALINDROMIC_ROOT_SCAN"
                                              : "WORD_SCAN")))
                      << " input=" << show(labels)
                      << " word=" << show(word)
                      << " domain_pairs=" << scan.domain_pairs
                      << " candidate_free_pairs=" << scan.candidate_free_pairs
                      << " witness=" << scan.witness
                      << " result=PASS"
                      << '\n';
            return EXIT_SUCCESS;
        }
        if (argc != 2) {
            throw std::invalid_argument(
                "usage: probe_su2_heterogeneous_two_block_bridge MAX_LABEL "
                "or --root LABEL... or --paired-root LABEL... "
                "or --palindromic-root LABEL... "
                "or --word LABEL... or --find-obstruction LABEL... "
                "or --sorted-root-scan MAX_LABEL MAX_ROOT_LENGTH"
            );
        }
        const int maximum_label = parse_positive(argv[1], "MAX_LABEL");
        unsigned long long label_pairs = 0U;
        unsigned long long word_orders = 0U;
        unsigned long long domain_pairs = 0U;
        unsigned long long candidate_free_pairs = 0U;
        unsigned long long obstructed_orders = 0U;
        unsigned long long covered_orders = 0U;
        std::string witness;

        for (int lower = 1; lower <= maximum_label; ++lower) {
            for (int upper = lower + 1; upper <= maximum_label; ++upper) {
                ++label_pairs;
                std::vector<int> word{lower, lower, upper, upper};
                do {
                    ++word_orders;
                    const Scan scan = scan_order(word);
                    domain_pairs += scan.domain_pairs;
                    candidate_free_pairs += scan.candidate_free_pairs;
                    if (!scan.has_obstruction) {
                        ++covered_orders;
                        continue;
                    }
                    ++obstructed_orders;
                    if (witness.empty()) {
                        witness = scan.witness;
                    }
                } while (std::next_permutation(word.begin(), word.end()));
            }
        }

        std::cout << "SU2_HETEROGENEOUS_TWO_BLOCK_BRIDGE_SCAN"
                  << " maximum_label=" << maximum_label
                  << " label_pairs=" << label_pairs
                  << " word_orders=" << word_orders
                  << " obstructed_orders=" << obstructed_orders
                  << " covered_orders=" << covered_orders
                  << " domain_pairs=" << domain_pairs
                  << " candidate_free_pairs=" << candidate_free_pairs
                  << " witness=" << witness
                  << " result=PASS"
                  << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
