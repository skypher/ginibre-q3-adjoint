#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Result {
    std::size_t sources = 0U;
    std::size_t prefix_failures = 0U;
    std::size_t collisions = 0U;
    std::size_t missing_parity_contacts = 0U;
    std::size_t transfer_failures = 0U;
    int witness_level = -1;
    int witness_label = -1;
    int witness_length = -1;
    int witness_truncation = -1;
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
    int contact_mode
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
    std::optional<std::size_t> chosen_second_index;
    std::size_t chosen_loop_position = 0U;
    for (std::size_t second_index = 0U;
         second_index < second.size();
         ++second_index) {
        const auto found = position.find(second[second_index]);
        if (found == position.end()) {
            continue;
        }
        const int first_index = exit_index[found->second];
        const bool parity_compatible =
            ((first_index + static_cast<int>(second_index)) & 1) != 0;
        if (contact_mode == 2 && !parity_compatible) {
            continue;
        }
        chosen_second_index = second_index;
        chosen_loop_position = found->second;
        if (contact_mode == 0) {
            break;
        }
    }
    if (chosen_second_index.has_value()) {
        const int first_index = exit_index[chosen_loop_position];
        const std::size_t second_index = *chosen_second_index;
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
    throw std::runtime_error(
        contact_mode == 2
            ? "the two paths have no parity-compatible loop-erased intersection"
            : "the two paths have no loop-erased intersection"
    );
}

bool valid_path(int level, int label, const std::vector<int>& path) {
    for (std::size_t index = 0U; index + 1U < path.size(); ++index) {
        if (!fuses(level, label, path[index], path[index + 1U])) {
            return false;
        }
    }
    return true;
}

bool transfer_two_step_returns(
    std::vector<int>& first,
    std::vector<int>& second,
    int label,
    int maximum_second_length,
    bool terminal_only
) {
    while (static_cast<int>(second.size()) - 1 > maximum_second_length) {
        if (first.empty() || first.back() != label || second.size() < 3U) {
            return false;
        }
        if (terminal_only) {
            if (second[second.size() - 3U] != 0
                || second[second.size() - 2U] != label
                || second.back() != 0) {
                return false;
            }
            second.erase(second.end() - 2, second.end());
            first.push_back(0);
            first.push_back(label);
            continue;
        }
        std::optional<std::size_t> last_return;
        for (std::size_t index = 0U; index + 2U < second.size(); ++index) {
            if (second[index] == second[index + 2U]) {
                last_return = index;
            }
        }
        if (!last_return.has_value()) {
            return false;
        }
        const auto erase_begin = second.begin()
            + static_cast<std::ptrdiff_t>(*last_return + 1U);
        second.erase(erase_begin, erase_begin + 2);
        first.push_back(label);
        first.push_back(label);
    }
    return true;
}

Result probe(
    int level,
    int label,
    int length,
    int truncation,
    bool reverse_orientation,
    int contact_mode,
    int transfer_mode
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
                ++result.sources;
                std::pair<std::vector<int>, std::vector<int>> switched;
                try {
                    switched = reverse_orientation
                        ? fomin_switch(second_path, first_path, contact_mode)
                        : fomin_switch(first_path, second_path, contact_mode);
                } catch (const std::runtime_error&) {
                    ++result.missing_parity_contacts;
                    ++result.prefix_failures;
                    if (result.first_source.empty()) {
                        result.first_source = key(first_path, second_path);
                    }
                    return;
                }
                std::vector<int> image_first = reverse_orientation
                    ? switched.second
                    : switched.first;
                std::vector<int> image_second = reverse_orientation
                    ? switched.first
                    : switched.second;
                if (transfer_mode != 0
                    && !transfer_two_step_returns(
                        image_first,
                        image_second,
                        label,
                        2 * truncation,
                        transfer_mode == 1
                    )) {
                    ++result.transfer_failures;
                    ++result.prefix_failures;
                    if (result.first_source.empty()) {
                        result.first_source = key(first_path, second_path);
                        result.first_image = key(image_first, image_second);
                    }
                    return;
                }
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

void accumulate(Result& total, const Result& part) {
    total.sources += part.sources;
    total.prefix_failures += part.prefix_failures;
    total.collisions += part.collisions;
    total.missing_parity_contacts += part.missing_parity_contacts;
    total.transfer_failures += part.transfer_failures;
    if (total.first_source.empty() && !part.first_source.empty()) {
        total.first_source = part.first_source;
        total.first_image = part.first_image;
        total.first_preimage = part.first_preimage;
    }
}

Result scan_direct_parity_last(
    int maximum_level,
    int maximum_length,
    int transfer_mode
) {
    Result result;
    for (int level = 4; level <= maximum_level; level += 2) {
        for (int label = 2; 2 * label < level; label += 2) {
            for (int length = 3; length <= maximum_length; length += 2) {
                for (int truncation = 0;
                     2 * truncation <= length - 3;
                     ++truncation) {
                    const bool have_witness = !result.first_source.empty();
                    const Result part = probe(
                        level,
                        label,
                        length,
                        truncation,
                        false,
                        2,
                        transfer_mode
                    );
                    accumulate(
                        result,
                        part
                    );
                    if (!have_witness && !part.first_source.empty()) {
                        result.witness_level = level;
                        result.witness_label = label;
                        result.witness_length = length;
                        result.witness_truncation = truncation;
                    }
                }
            }
        }
    }
    return result;
}

void replay_first_prefix_collision() {
    const Result result = probe(6, 2, 9, 3, false, 0, false);
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
    const Result result = probe(6, 2, 9, 3, true, 0, false);
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
        if (argc == 4 && (std::string(argv[1]) == "--scan"
            || std::string(argv[1]) == "--transfer-scan"
            || std::string(argv[1]) == "--excursion-transfer-scan")) {
            const int maximum_level =
                parse_nonnegative(argv[2], "maximum level");
            const int maximum_length =
                parse_nonnegative(argv[3], "maximum odd length");
            if (maximum_level < 4 || maximum_length < 3
                || (maximum_length & 1) == 0) {
                throw std::runtime_error(
                    "scan bounds require even level at least four and odd length at least three"
                );
            }
            const std::string scan_mode(argv[1]);
            const int transfer_mode = scan_mode == "--transfer-scan"
                ? 1
                : (scan_mode == "--excursion-transfer-scan" ? 2 : 0);
            const Result result = scan_direct_parity_last(
                maximum_level,
                maximum_length,
                transfer_mode
            );
            std::cout
                << "SU2_FOMIN_PARITY_LAST_SCAN"
                << " transfer_mode=" << transfer_mode
                << " maximum_level=" << maximum_level
                << " maximum_odd_length=" << maximum_length
                << " sources=" << result.sources
                << " prefix_failures=" << result.prefix_failures
                << " collisions=" << result.collisions
                << " missing_parity_contacts="
                << result.missing_parity_contacts
                << " transfer_failures=" << result.transfer_failures;
            if (!result.first_source.empty()) {
                std::cout
                    << " witness_level=" << result.witness_level
                    << " witness_label=" << result.witness_label
                    << " witness_odd_length=" << result.witness_length
                    << " witness_truncation=" << result.witness_truncation;
                std::cout << " first_source=";
                for (const int value : result.first_source) {
                    std::cout << value << ',';
                }
                std::cout << " first_image=";
                for (const int value : result.first_image) {
                    std::cout << value << ',';
                }
            }
            std::cout
                << " result="
                << (result.prefix_failures == 0U && result.collisions == 0U
                    ? "PASS_EXACT_BOX" : "COUNTEREXAMPLE")
                << '\n';
            return result.prefix_failures == 0U && result.collisions == 0U
                ? 0 : 1;
        }
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
                || std::string(argv[1]) == "--reverse-last"
                || std::string(argv[1]) == "--parity-last"
                || std::string(argv[1]) == "--reverse-parity-last"
                || std::string(argv[1]) == "--parity-last-transfer"
                || std::string(argv[1]) == "--parity-last-excursion-transfer");
        if (argc != 5 && !switched_mode) {
            throw std::runtime_error(
                "usage: probe_su2_fomin_prefix_switch LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --reverse LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --last LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --reverse-last LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --parity-last LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --reverse-parity-last LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --parity-last-transfer LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --parity-last-excursion-transfer LEVEL LABEL ODD_LENGTH TRUNCATION"
                " | --scan MAXIMUM_EVEN_LEVEL MAXIMUM_ODD_LENGTH"
                " | --transfer-scan MAXIMUM_EVEN_LEVEL MAXIMUM_ODD_LENGTH"
                " | --excursion-transfer-scan MAXIMUM_EVEN_LEVEL MAXIMUM_ODD_LENGTH"
                " | --replay-first-prefix-collision"
                " | --replay-reverse-prefix-obstruction"
            );
        }
        const std::string mode = switched_mode ? std::string(argv[1]) : "";
        const bool reverse_orientation = mode == "--reverse"
            || mode == "--reverse-last"
            || mode == "--reverse-parity-last";
        const int contact_mode = (mode == "--last" || mode == "--reverse-last")
            ? 1
            : ((mode == "--parity-last" || mode == "--reverse-parity-last"
                || mode == "--parity-last-transfer"
                || mode == "--parity-last-excursion-transfer")
                ? 2
                : 0);
        const int transfer_mode = mode == "--parity-last-transfer"
            ? 1
            : (mode == "--parity-last-excursion-transfer" ? 2 : 0);
        const int offset = switched_mode ? 1 : 0;
        const Result result = probe(
            parse_nonnegative(argv[1 + offset], "level"),
            parse_nonnegative(argv[2 + offset], "label"),
            parse_nonnegative(argv[3 + offset], "odd length"),
            parse_nonnegative(argv[4 + offset], "truncation"),
            reverse_orientation,
            contact_mode,
            transfer_mode
        );
        std::cout << "SU2_FOMIN_PREFIX_SWITCH"
                  << " orientation="
                  << (reverse_orientation ? "reverse" : "direct")
                  << " contact="
                  << (contact_mode == 0 ? "first"
                      : (contact_mode == 1 ? "last" : "parity-last"))
                  << " sources=" << result.sources
                  << " prefix_failures=" << result.prefix_failures
                  << " collisions=" << result.collisions
                  << " missing_parity_contacts="
                  << result.missing_parity_contacts
                  << " transfer_mode=" << transfer_mode
                  << " transfer_failures=" << result.transfer_failures;
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
