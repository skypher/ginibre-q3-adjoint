#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Lines = std::vector<std::string>;

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

Lines read_lines(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open " + path.string());
    }
    Lines lines;
    std::string line;
    while (std::getline(input, line)) {
        lines.push_back(line);
    }
    return lines;
}

bool has_bad_status(const Lines& lines) {
    for (const std::string& line : lines) {
        if (
            line.find("FAILURE") != std::string::npos
            || line.find("result=FAIL") != std::string::npos
            || line.find("UNRESOLVED") != std::string::npos
            || line.find("INCOMPLETE") != std::string::npos
            || line.find("ERROR") != std::string::npos
        ) {
            return true;
        }
    }
    return false;
}

struct MaskList {
    std::string header;
    std::vector<std::uint64_t> masks;
};

MaskList read_masks(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open mask list");
    }
    MaskList result;
    if (
        !std::getline(input, result.header)
        || result.header
            != "SU2_K4_INTERMEDIATE_MASKS hinges=50 masks=302 "
               "minimum_active=7 maximum_active=44 "
               "result=PASS_EXACT_CENSUS"
    ) {
        throw std::runtime_error("mask list has no exact census header");
    }
    const std::string prefix = "SU2_K4_INTERMEDIATE_MASK value=";
    std::string line;
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
        result.masks.push_back(static_cast<std::uint64_t>(value));
    }
    const std::set<std::uint64_t> unique(
        result.masks.begin(),
        result.masks.end()
    );
    if (
        result.masks.size() != 302U
        || unique.size() != result.masks.size()
        || !std::is_sorted(result.masks.begin(), result.masks.end())
    ) {
        throw std::runtime_error(
            "mask list does not contain 302 unique sorted masks"
        );
    }
    return result;
}

bool valid_k4_log(const Lines& lines, std::uint64_t mask) {
    if (
        lines.empty()
        || has_bad_status(lines)
        || lines.back()
            != "SU2_K4_INTERMEDIATE hinges=50 feasible_chambers=1 "
               "certified_chambers=1 result=PASS_EXACT_CERTIFICATE"
    ) {
        return false;
    }
    const std::string field =
        " mask=" + std::to_string(mask);
    return std::any_of(
        lines.begin(),
        lines.end(),
        [&field](const std::string& line) {
            const std::size_t position = line.find(field);
            if (position == std::string::npos) {
                return false;
            }
            const std::size_t after = position + field.size();
            return (after == line.size() || line[after] == ' ')
                && line.find("result=PASS_") != std::string::npos;
        }
    );
}

std::string group_census(const std::string& target) {
    if (target == "g0") {
        return "SU2_T4_GROUP_MASKS target=g0 hinges=88 masks=541 "
               "result=PASS_EXACT_CENSUS";
    }
    if (target == "g1") {
        return "SU2_T4_GROUP_MASKS target=g1 hinges=112 masks=1065 "
               "result=PASS_EXACT_CENSUS";
    }
    if (target == "g2") {
        return "SU2_T4_GROUP_MASKS target=g2 hinges=91 masks=848 "
               "result=PASS_EXACT_CENSUS";
    }
    throw std::runtime_error("target must be g0, g1, or g2");
}

bool valid_group_log(
    const Lines& lines,
    const std::string& target
) {
    const std::string terminal =
        "SU2_T4_GROUP target=" + target
        + " attempted=1 certified=1 result=PASS_EXACT_CERTIFICATE";
    const std::string census = group_census(target);
    return !lines.empty()
        && !has_bad_status(lines)
        && lines.back() == terminal
        && std::find(lines.begin(), lines.end(), census) != lines.end();
}

template <typename Validator>
Lines select_log(
    const std::vector<std::filesystem::path>& directories,
    const std::string& filename,
    Validator validator
) {
    for (const std::filesystem::path& directory : directories) {
        const std::filesystem::path path = directory / filename;
        if (!std::filesystem::is_regular_file(path)) {
            continue;
        }
        const Lines lines = read_lines(path);
        if (validator(lines)) {
            return lines;
        }
    }
    throw std::runtime_error("no complete exact log for " + filename);
}

void collect_k4(
    const std::filesystem::path& mask_path,
    const std::vector<std::filesystem::path>& directories
) {
    const MaskList list = read_masks(mask_path);
    std::cout << list.header << '\n';
    for (const std::uint64_t mask : list.masks) {
        const Lines lines = select_log(
            directories,
            "k4_mask_" + std::to_string(mask) + ".log",
            [mask](const Lines& candidate) {
                return valid_k4_log(candidate, mask);
            }
        );
        for (std::size_t index = 0U;
             index + 1U < lines.size();
             ++index) {
            std::cout << lines[index] << '\n';
        }
    }
    std::cout
        << "SU2_K4_INTERMEDIATE hinges=50 feasible_chambers="
        << list.masks.size()
        << " certified_chambers=" << list.masks.size()
        << " result=PASS_EXACT_CERTIFICATE\n";
}

void list_missing_k4(
    const std::filesystem::path& mask_path,
    const std::vector<std::filesystem::path>& directories
) {
    const MaskList list = read_masks(mask_path);
    std::vector<std::uint64_t> missing;
    for (const std::uint64_t mask : list.masks) {
        bool found = false;
        const std::string filename =
            "k4_mask_" + std::to_string(mask) + ".log";
        for (const std::filesystem::path& directory : directories) {
            const std::filesystem::path path = directory / filename;
            if (!std::filesystem::is_regular_file(path)) {
                continue;
            }
            if (valid_k4_log(read_lines(path), mask)) {
                found = true;
                break;
            }
        }
        if (!found) {
            missing.push_back(mask);
        }
    }
    std::cout
        << "SU2_K4_INTERMEDIATE_MASK_SUBSET count="
        << missing.size() << '\n';
    for (const std::uint64_t mask : missing) {
        std::cout
            << "SU2_K4_INTERMEDIATE_MASK value="
            << mask << '\n';
    }
}

void collect_group(
    const std::string& target,
    std::size_t count,
    const std::vector<std::filesystem::path>& directories
) {
    const std::string census = group_census(target);
    const std::size_t expected =
        target == "g0" ? 541U : (target == "g1" ? 1065U : 848U);
    if (count != expected) {
        throw std::runtime_error("group count disagrees with census");
    }
    std::cout << census << '\n';
    for (std::size_t position = 0U; position < count; ++position) {
        const Lines lines = select_log(
            directories,
            target + "_position_" + std::to_string(position) + ".log",
            [&target](const Lines& candidate) {
                return valid_group_log(candidate, target);
            }
        );
        for (const std::string& line : lines) {
            if (
                line.starts_with("SU2_T4_GROUP_MASKS")
                || line.starts_with("SU2_T4_GROUP target=")
            ) {
                continue;
            }
            std::cout << line << '\n';
        }
    }
    std::cout
        << "SU2_T4_GROUP target=" << target
        << " attempted=" << count
        << " certified=" << count
        << " result=PASS_EXACT_CERTIFICATE\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 5) {
            throw std::runtime_error(
                "usage: k4|missing-k4 MASK_LIST LOG_DIRECTORY... "
                "| group TARGET COUNT LOG_DIRECTORY..."
            );
        }
        const std::string mode = argv[1];
        if (mode == "k4" || mode == "missing-k4") {
            std::vector<std::filesystem::path> directories;
            for (int index = 3; index < argc; ++index) {
                directories.emplace_back(argv[index]);
            }
            if (mode == "k4") {
                collect_k4(argv[2], directories);
            } else {
                list_missing_k4(argv[2], directories);
            }
        } else if (mode == "group") {
            if (argc < 6) {
                throw std::runtime_error(
                    "group mode needs target, count, and directories"
                );
            }
            std::vector<std::filesystem::path> directories;
            for (int index = 4; index < argc; ++index) {
                directories.emplace_back(argv[index]);
            }
            collect_group(
                argv[2],
                parse_size(argv[3], "group count"),
                directories
            );
        } else {
            throw std::runtime_error(
                "mode must be k4, missing-k4, or group"
            );
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_EXACT_SHARD_COLLECTOR FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
