#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Pair {
    std::vector<int> short_path;
    std::vector<int> long_path;

    friend bool operator<(const Pair& left, const Pair& right) {
        if (left.short_path != right.short_path) {
            return left.short_path < right.short_path;
        }
        return left.long_path < right.long_path;
    }
};

int parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0U;
    const long parsed = std::stol(value, &consumed);
    if (consumed != value.size() || parsed <= 0L) {
        throw std::runtime_error(name + " must be positive");
    }
    return static_cast<int>(parsed);
}

bool edge(int left, int right, int q) {
    return std::abs(left - q) <= right && right <= left + q;
}

void enumerate_from(
    int q,
    int steps,
    std::vector<int>& path,
    std::vector<std::vector<int>>& paths
) {
    if (static_cast<int>(path.size()) == steps + 1) {
        paths.push_back(path);
        return;
    }
    const int source = path.back();
    for (int target = std::abs(source - q);
         target <= source + q;
         ++target) {
        path.push_back(target);
        enumerate_from(q, steps, path, paths);
        path.pop_back();
    }
}

std::vector<std::vector<int>> paths_from_zero(int q, int steps) {
    std::vector<std::vector<int>> paths;
    std::vector<int> path{0};
    enumerate_from(q, steps, path, paths);
    return paths;
}

std::string show(const Pair& pair) {
    std::string result{"{"};
    for (std::size_t index = 0U; index < pair.short_path.size(); ++index) {
        if (index != 0U) {
            result.push_back(',');
        }
        result += std::to_string(pair.short_path[index]);
    }
    result += "};{";
    for (std::size_t index = 0U; index < pair.long_path.size(); ++index) {
        if (index != 0U) {
            result.push_back(',');
        }
        result += std::to_string(pair.long_path[index]);
    }
    result.push_back('}');
    return result;
}

bool consecutive(const std::set<int>& indices) {
    if (indices.empty()) {
        return true;
    }
    return static_cast<std::size_t>(*indices.rbegin() - *indices.begin() + 1)
        == indices.size();
}

int area(const Pair& pair) {
    int result = 0;
    for (const int value : pair.short_path) {
        result += value;
    }
    for (const int value : pair.long_path) {
        result += value;
    }
    return result;
}

bool area_then_lex(const Pair& left, const Pair& right) {
    if (area(left) != area(right)) {
        return area(left) < area(right);
    }
    return left < right;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 4) {
            throw std::runtime_error(
                "usage: probe_su2_kostka_bridge_biconvex MAX_Q MAX_LENGTH "
                "[--lex|--area]"
            );
        }
        const int maximum_q = parse_positive(argv[1], "MAX_Q");
        const int maximum_length = parse_positive(argv[2], "MAX_LENGTH");
        const bool area_order = argc == 4 && std::string(argv[3]) == "--area";
        const bool lex_order = argc == 3 || (argc == 4 && std::string(argv[3]) == "--lex");
        if (!area_order && !lex_order) {
            throw std::runtime_error("the optional order is --lex or --area");
        }
        if ((maximum_length & 1) != 0 || maximum_length < 2) {
            throw std::runtime_error("MAX_LENGTH must be an even integer at least two");
        }

        std::size_t tasks = 0U;
        std::size_t left_rows = 0U;
        std::size_t right_rows = 0U;
        for (int q = 1; q <= maximum_q; ++q) {
            for (int length = 2; length <= maximum_length; length += 2) {
                ++tasks;
                const auto short_paths = paths_from_zero(q, length);
                const auto all_long_paths = paths_from_zero(q, length + 1);
                std::vector<std::vector<int>> returning;
                for (const auto& path : all_long_paths) {
                    if (path.back() == 0) {
                        returning.push_back(path);
                    }
                }

                std::vector<Pair> left;
                std::map<Pair, std::set<Pair>> relation;
                for (const auto& long_path : returning) {
                    for (const auto& short_path : short_paths) {
                        if (short_path.back() == 0) {
                            continue;
                        }
                        const Pair source{short_path, long_path};
                        left.push_back(source);
                        std::set<Pair>& targets = relation[source];
                        for (int start = 0; start < length; ++start) {
                            if (edge(short_path[static_cast<std::size_t>(start)],
                                     long_path[static_cast<std::size_t>(start + 2)], q)
                                && edge(long_path[static_cast<std::size_t>(start + 1)],
                                        short_path[static_cast<std::size_t>(start + 1)], q)) {
                                Pair image;
                                image.short_path.resize(
                                    static_cast<std::size_t>(length + 1), 0
                                );
                                image.long_path.resize(
                                    static_cast<std::size_t>(length + 2), 0
                                );
                                image.long_path[0] = 0;
                                for (int time = 0; time <= length; ++time) {
                                    if (time <= start) {
                                        image.short_path[static_cast<std::size_t>(time)] =
                                            short_path[static_cast<std::size_t>(time)];
                                        image.long_path[static_cast<std::size_t>(time + 1)] =
                                            long_path[static_cast<std::size_t>(time + 1)];
                                    } else {
                                        image.short_path[static_cast<std::size_t>(time)] =
                                            long_path[static_cast<std::size_t>(time + 1)];
                                        image.long_path[static_cast<std::size_t>(time + 1)] =
                                            short_path[static_cast<std::size_t>(time)];
                                    }
                                }
                                targets.insert(std::move(image));
                            }
                        }
                        for (int start = 0; start + 2 <= length; ++start) {
                            const int short_start =
                                short_path[static_cast<std::size_t>(start)];
                            const int short_end =
                                long_path[static_cast<std::size_t>(start + 3)];
                            const int long_start =
                                long_path[static_cast<std::size_t>(start + 1)];
                            const int long_end =
                                short_path[static_cast<std::size_t>(start + 2)];
                            for (int short_middle = 0;
                                 short_middle <= (length + 1) * q;
                                 ++short_middle) {
                                if (!edge(short_start, short_middle, q)
                                    || !edge(short_middle, short_end, q)) {
                                    continue;
                                }
                                for (int long_middle = 0;
                                     long_middle <= (length + 1) * q;
                                     ++long_middle) {
                                    if (!edge(long_start, long_middle, q)
                                        || !edge(long_middle, long_end, q)) {
                                        continue;
                                    }
                                    Pair image;
                                    image.short_path.resize(
                                        static_cast<std::size_t>(length + 1), 0
                                    );
                                    image.long_path.resize(
                                        static_cast<std::size_t>(length + 2), 0
                                    );
                                    image.long_path[0] = 0;
                                    for (int time = 0; time <= length; ++time) {
                                        if (time <= start) {
                                            image.short_path[static_cast<std::size_t>(time)] =
                                                short_path[static_cast<std::size_t>(time)];
                                            image.long_path[static_cast<std::size_t>(time + 1)] =
                                                long_path[static_cast<std::size_t>(time + 1)];
                                        } else if (time >= start + 2) {
                                            image.short_path[static_cast<std::size_t>(time)] =
                                                long_path[static_cast<std::size_t>(time + 1)];
                                            image.long_path[static_cast<std::size_t>(time + 1)] =
                                                short_path[static_cast<std::size_t>(time)];
                                        } else {
                                            image.short_path[static_cast<std::size_t>(time)] =
                                                short_middle;
                                            image.long_path[static_cast<std::size_t>(time + 1)] =
                                                long_middle;
                                        }
                                    }
                                    targets.insert(std::move(image));
                                }
                            }
                        }
                    }
                }

                std::set<Pair> right_set;
                for (const auto& [source, targets] : relation) {
                    static_cast<void>(source);
                    right_set.insert(targets.begin(), targets.end());
                }
                if (area_order) {
                    std::sort(left.begin(), left.end(), area_then_lex);
                }
                std::vector<Pair> right(right_set.begin(), right_set.end());
                if (area_order) {
                    std::sort(right.begin(), right.end(), area_then_lex);
                }
                std::map<Pair, int> right_index;
                for (std::size_t index = 0U; index < right.size(); ++index) {
                    right_index.emplace(right[index], static_cast<int>(index));
                }
                std::map<Pair, int> left_index;
                for (std::size_t index = 0U; index < left.size(); ++index) {
                    left_index.emplace(left[index], static_cast<int>(index));
                }
                std::vector<std::set<int>> reverse(right.size());
                for (const Pair& source : left) {
                    std::set<int> neighborhood;
                    for (const Pair& image : relation.at(source)) {
                        const int image_index = right_index.at(image);
                        neighborhood.insert(image_index);
                        reverse[static_cast<std::size_t>(image_index)].insert(
                            left_index.at(source)
                        );
                    }
                    ++left_rows;
                    if (!consecutive(neighborhood)) {
                        std::cout
                            << "SU2_KOSTKA_BRIDGE_BICONVEX left_failure"
                            << " q=" << q << " length=" << length
                            << " order=" << (area_order ? "area" : "lex")
                            << " source=" << show(source)
                            << " degree=" << neighborhood.size()
                            << " first=" << *neighborhood.begin()
                            << " last=" << *neighborhood.rbegin()
                            << " neighborhood={";
                        bool first = true;
                        for (const int index : neighborhood) {
                            if (!first) {
                                std::cout << ',';
                            }
                            std::cout << index;
                            first = false;
                        }
                        std::cout << "} result=NOT_BICONVEX\n";
                        return EXIT_FAILURE;
                    }
                }
                for (std::size_t index = 0U; index < right.size(); ++index) {
                    ++right_rows;
                    if (!consecutive(reverse[index])) {
                        std::cout
                            << "SU2_KOSTKA_BRIDGE_BICONVEX right_failure"
                            << " q=" << q << " length=" << length
                            << " order=" << (area_order ? "area" : "lex")
                            << " image=" << show(right[index])
                            << " degree=" << reverse[index].size()
                            << " first=" << *reverse[index].begin()
                            << " last=" << *reverse[index].rbegin()
                            << " result=NOT_BICONVEX\n";
                        return EXIT_FAILURE;
                    }
                }
            }
        }
        std::cout << "SU2_KOSTKA_BRIDGE_BICONVEX"
                  << " tasks=" << tasks
                  << " left_rows=" << left_rows
                  << " right_rows=" << right_rows
                  << " max_q=" << maximum_q
                  << " max_length=" << maximum_length
                  << " order=" << (area_order ? "area" : "lex")
                  << " result=PASS_NATURAL_LEX_BICONVEX\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_KOSTKA_BRIDGE_BICONVEX failure: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
