#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kChambers = 1181U;

struct RangeLog {
    std::size_t begin = 0U;
    std::size_t end = 0U;
    std::filesystem::path path;
};

std::string read_file(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open " + path.string());
    }
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()
    );
}

bool contains(const std::string& text, const std::string& needle) {
    return text.find(needle) != std::string::npos;
}

std::size_t parse_size(const std::ssub_match& match) {
    const std::string value = match.str();
    std::size_t used = 0U;
    const unsigned long long parsed = std::stoull(value, &used);
    if (used != value.size()) {
        throw std::runtime_error("invalid range endpoint");
    }
    return static_cast<std::size_t>(parsed);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: collect_su2_t4_c6_ranges LOG_DIRECTORY"
            );
        }
        const std::filesystem::path directory{argv[1]};
        if (!std::filesystem::is_directory(directory)) {
            throw std::runtime_error("log path is not a directory");
        }

        const std::regex name_pattern{
            R"(^shard_[0-9]+_([0-9]+)_([0-9]+)\.log$)"
        };
        const std::string census =
            "SU2_T4_GROUP_MASKS target=c6 hinges=98 masks=1181 "
            "result=PASS_EXACT_CENSUS";
        std::vector<RangeLog> logs;
        for (const std::filesystem::directory_entry& entry
             : std::filesystem::directory_iterator(directory)) {
            if (!entry.is_regular_file()) {
                continue;
            }
            std::smatch match;
            const std::string name = entry.path().filename().string();
            if (!std::regex_match(name, match, name_pattern)) {
                continue;
            }
            const std::size_t begin = parse_size(match[1]);
            const std::size_t end = parse_size(match[2]);
            if (begin >= end || end > kChambers) {
                throw std::runtime_error("invalid shard interval " + name);
            }
            const std::string transcript = read_file(entry.path());
            if (
                contains(transcript, "UNRESOLVED")
                || contains(transcript, "INCOMPLETE")
                || contains(transcript, "FAILURE")
            ) {
                throw std::runtime_error(
                    "nonterminal failure marker in " + name
                );
            }
            const std::string terminal =
                "SU2_T4_GROUP target=c6 attempted="
                + std::to_string(end - begin)
                + " certified=" + std::to_string(end - begin)
                + " scope=range result=PASS_EXACT_RANGE";
            if (!contains(transcript, census) || !contains(transcript, terminal)) {
                throw std::runtime_error(
                    "missing exact terminal record in " + name
                );
            }
            logs.push_back(RangeLog{begin, end, entry.path()});
        }

        std::sort(
            logs.begin(), logs.end(),
            [](const RangeLog& left, const RangeLog& right) {
                return left.begin < right.begin;
            }
        );
        std::size_t expected_begin = 0U;
        for (const RangeLog& log : logs) {
            if (log.begin != expected_begin) {
                throw std::runtime_error(
                    "missing or overlapping range before " + log.path.string()
                );
            }
            expected_begin = log.end;
        }
        if (expected_begin != kChambers) {
            throw std::runtime_error("range union does not cover all chambers");
        }

        std::cout
            << "SU2_C6_RANGE_AGGREGATE"
            << " shards=" << logs.size()
            << " chambers=" << kChambers
            << " exact_pass=" << kChambers
            << " bad=0 missing=0"
            << " result=PASS_EXACT_CERTIFICATE"
            << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_C6_RANGE_AGGREGATE FAILURE: " << error.what() << '\n';
        return 1;
    }
}
