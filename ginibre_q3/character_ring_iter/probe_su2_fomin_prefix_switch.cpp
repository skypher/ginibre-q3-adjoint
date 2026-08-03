#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Result {
    std::size_t sources = 0U;
    std::size_t prefix_failures = 0U;
    std::size_t collisions = 0U;
    std::vector<int> first_source;
    std::vector<int> first_image;
    std::vector<int> first_preimage;
};

int parse_nonnegative(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value < 0L
        || value > static_cast<long>(1 << 20)) {
        throw std::runtime_error(std::string(name) + " must be nonnegative");
    }
    return static_cast<int>(value);
}

bool fuses(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= std::min(source + label, 2 * level - source - label)
        && ((source + label + target) & 1) == 0;
}

std::vector<int> key(const std::vector<int>& first, const std::vector<int>& second) {
    std::vector<int> result;
    result.reserve(first.size() + second.size() + 2U);
    result.push_back(static_cast<int>(first.size()));
    result.insert(result.end(), first.begin(), first.end());
    result.push_back(static_cast<int>(second.size()));
    result.insert(result.end(), second.begin(), second.end());
    return result;
}

std::pair<std::vector<int>, std::vector<int>> fomin_switch(
    const std::vector<int>& first,
    const std::vector<int>& second,
    bool last_intersection
) {
    std::vector<int> loop_erasure;
    std::vector<int> exit_index;
    std::map<int, std::size_t> position;
    for (std::size_t index = 0U; index < first.size(); ++index) {
        const int vertex = first[index];
        const auto found = position.find(vertex);
        if (found == position.end()) {
            position.emplace(vertex, loop_erasure.size());
            loop_erasure.push_back(vertex);
            exit_index.push_back(static_cast<int>(index));
            continue;
        }
        const std::size_t keep = found->second;
        for (std::size_t erased = keep + 1U; erased < loop_erasure.size(); ++erased) {
            position.erase(loop_erasure[erased]);
        }
        loop_erasure.resize(keep + 1U);
        exit_index.resize(keep + 1U);
        exit_index[keep] = static_cast<int>(index);
    }
    int chosen_second_index = -1;
    std::size_t chosen_loop_position = 0U;
    for (std::size_t second_index = 0U;
         second_index < second.size();
         ++second_index) {
        const auto found = position.find(second[second_index]);
        if (found == position.end()) {
            continue;
        }
        chosen_second_index = static_cast<int>(second_index);
        chosen_loop_position = found->second;
        if (!last_intersection) {
            break;
        }
    }
    if (chosen_second_index >= 0) {
        const int first_index = exit_index[chosen_loop_position];
        const std::size_t second_index =
            static_cast<std::size_t>(chosen_second_index);
        std::vector<int> switched_first(
            first.begin(), first.begin() + first_index + 1
        );
        switched_first.insert(
            switched_first.end(), std::next(
                second.begin(), static_cast<std::ptrdiff_t>(second_index) + 1
            ),
            second.end()
        );
        std::vector<int> switched_second(
            second.begin(), std::next(
                second.begin(), static_cast<std::ptrdiff_t>(second_index) + 1
            )
        );
        switched_second.insert(
            switched_second.end(), first.begin() + first_index + 1,
            first.end()
        );
        return {switched_first, switched_second};
    }
    throw std::runtime_error("the two paths have no loop-erased intersection");
}

bool valid_path(int level, int label, const std::vector<int>& path) {
    for (std::size_t index = 0U; index + 1U < path.size(); ++index) {
        if (!fuses(level, label, path[index], path[index + 1U])) {
            return false;
        }
    }
    return true;
}

Result probe(
    int level,
    int label,
    int length,
    int truncation,
    bool reverse_orientation,
    bool last_intersection
) {
    if (level <= 0 || label <= 0 || label >= level || (level & 1) != 0
        || (label & 1) != 0 || 2 * label >= level || (length & 1) == 0
        || truncation < 0 || 2 * truncation >= length) {
        throw std::runtime_error("parameters are outside the reflected-prefix range");
    }
    Result result;
    std::map<std::vector<int>, std::vector<int>> images;
    for (int second_length = 0; second_length <= 2 * truncation;
         second_length += 2) {
        const int first_length = length - second_length;
        std::vector<int> first_path{level};
        std::vector<int> second_path{0};
        const auto enumerate_second = [&](const auto& self, int remaining) -> void {
            if (remaining == 0) {
                if (second_path.back() != label) {
                    return;
                }
                const auto switched = reverse_orientation
                    ? fomin_switch(second_path, first_path, last_intersection)
                    : fomin_switch(first_path, second_path, last_intersection);
                const std::vector<int>& image_first = reverse_orientation
                    ? switched.second
                    : switched.first;
                const std::vector<int>& image_second = reverse_orientation
                    ? switched.first
                    : switched.second;
                ++result.sources;
                const bool valid = image_first.front() == level
                    && image_first.back() == label
                    && image_second.front() == 0
                    && image_second.back() == 0
                    && valid_path(level, label, image_first)
                    && valid_path(level, label, image_second)
                    && ((static_cast<int>(image_second.size()) - 1) & 1) == 0
                    && static_cast<int>(image_second.size()) - 1
                        <= 2 * truncation;
                const std::vector<int> source_key = key(first_path, second_path);
                const std::vector<int> image_key = key(image_first, image_second);
                if (!valid) {
                    ++result.prefix_failures;
                    if (result.first_source.empty()) {
                        result.first_source = source_key;
                        result.first_image = image_key;
                    }
                    return;
                }
                const auto [position, inserted] = images.emplace(image_key, source_key);
                if (!inserted && position->second != source_key) {
                    ++result.collisions;
                    if (result.first_source.empty()) {
                        result.first_source = source_key;
                        result.first_image = image_key;
                        result.first_preimage = position->second;
                    }
                }
                return;
            }
            for (int target = 0; target <= level; ++target) {
                if (!fuses(level, label, second_path.back(), target)) {
                    continue;
                }
                second_path.push_back(target);
                self(self, remaining - 1);
                second_path.pop_back();
            }
        };
        const auto enumerate_first = [&](const auto& self, int remaining) -> void {
            if (remaining == 0) {
                if (first_path.back() != 0) {
                    return;
                }
                enumerate_second(enumerate_second, second_length);
                return;
            }
            for (int target = 0; target <= level; ++target) {
                if (!fuses(level, label, first_path.back(), target)) {
                    continue;
                }
                first_path.push_back(target);
                self(self, remaining - 1);
                first_path.pop_back();
            }
        };
        enumerate_first(enumerate_first, first_length);
    }
    return result;
}

void replay_first_prefix_collision() {
    const Result result = probe(6, 2, 9, 3, false, false);
    if (result.sources != 106U || result.prefix_failures != 0U
        || result.collisions != 14U) {
        throw std::runtime_error("Fomin-prefix collision replay mismatch");
    }
    std::cout << "SU2_FOMIN_PREFIX_COLLISION"
              << " sources=106"
              << " prefix_failures=0"
              << " collisions=14"
              << " result=PASS_EXACT\n";
}

void replay_reverse_prefix_obstruction() {
    const Result result = probe(6, 2, 9, 3, true, false);
    if (result.sources != 106U || result.prefix_failures != 91U
        || result.collisions != 3U) {
        throw std::runtime_error("reverse Fomin-prefix replay mismatch");
    }
    std::cout << "SU2_FOMIN_REVERSE_PREFIX_OBSTRUCTION"
              << " sources=106"
              << " prefix_failures=91"
              << " collisions=3"
              << " result=PASS_EXACT\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2
            && std::string(argv[1]) == "--replay-first-prefix-collision") {
            replay_first_prefix_collision();
            return 0;
        }
        if (argc == 2
            && std::string(argv[1]) == "--replay-reverse-prefix-obstruction") {
            replay_reverse_prefix_obstruction();
            return 0;
        }
        const bool switched_mode = argc == 6
            && (std::string(argv[1]) == "--reverse"
                || std::string(argv[1]) == "--last"
                || std::string(argv[1]) == "--reverse-last");
        if (argc != 5 && !switched_mode) {
            throw std::runtime_error(
                "usage: probe_su2_fomin_prefix_switch LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --reverse LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --last LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --reverse-last LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --replay-first-prefix-collision"
                " | --replay-reverse-prefix-obstruction"
            );
        }
        const std::string mode = switched_mode ? std::string(argv[1]) : "";
        const bool reverse_orientation = mode == "--reverse"
            || mode == "--reverse-last";
        const bool last_intersection = mode == "--last"
            || mode == "--reverse-last";
        const int offset = switched_mode ? 1 : 0;
        const Result result = probe(
            parse_nonnegative(argv[1 + offset], "level"),
            parse_nonnegative(argv[2 + offset], "label"),
            parse_nonnegative(argv[3 + offset], "odd length"),
            parse_nonnegative(argv[4 + offset], "truncation"),
            reverse_orientation,
            last_intersection
        );
        std::cout << "SU2_FOMIN_PREFIX_SWITCH"
                  << " orientation="
                  << (reverse_orientation ? "reverse" : "direct")
                  << " contact=" << (last_intersection ? "last" : "first")
                  << " sources=" << result.sources
                  << " prefix_failures=" << result.prefix_failures
                  << " collisions=" << result.collisions;
        if (!result.first_source.empty()) {
            std::cout << " first_source=";
            for (const int value : result.first_source) {
                std::cout << value << ',';
            }
            std::cout << " first_image=";
            for (const int value : result.first_image) {
                std::cout << value << ',';
            }
            if (!result.first_preimage.empty()) {
                std::cout << " first_preimage=";
                for (const int value : result.first_preimage) {
                    std::cout << value << ',';
                }
            }
        }
        std::cout << '\n';
        return (result.prefix_failures == 0U && result.collisions == 0U) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
