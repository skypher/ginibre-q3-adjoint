#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct CrystalStats {
    int epsilon = 0;
    int phi = 0;
};

CrystalStats crystal_stats(
    const std::vector<int>& labels,
    const std::vector<int>& state,
    std::size_t end
) {
    CrystalStats answer;
    for (std::size_t i = 0U; i < end; ++i) {
        const int factor_epsilon = state[i];
        const int factor_phi = labels[i] - state[i];
        const int next_epsilon = answer.epsilon
            + std::max(0, factor_epsilon - answer.phi);
        const int next_phi = factor_phi
            + std::max(0, answer.phi - factor_epsilon);
        answer = {next_epsilon, next_phi};
    }
    return answer;
}

int crystal_weight(
    const std::vector<int>& labels,
    const std::vector<int>& state
) {
    int answer = 0;
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        answer += labels[i] - 2 * state[i];
    }
    return answer;
}

std::size_t lowering_index_rec(
    const std::vector<int>& labels,
    const std::vector<int>& state,
    std::size_t end
) {
    if (end == 0U) {
        throw std::runtime_error("lowering an empty crystal word");
    }
    if (end == 1U) {
        if (state[0] >= labels[0]) {
            throw std::runtime_error("invalid crystal lowering");
        }
        return 0U;
    }
    const CrystalStats prefix = crystal_stats(labels, state, end - 1U);
    const int last_epsilon = state[end - 1U];
    if (prefix.phi > last_epsilon) {
        return lowering_index_rec(labels, state, end - 1U);
    }
    if (state[end - 1U] >= labels[end - 1U]) {
        throw std::runtime_error("invalid terminal crystal lowering");
    }
    return end - 1U;
}

void lower_word(
    const std::vector<int>& labels,
    std::vector<int>& state,
    int steps
) {
    for (int step = 0; step < steps; ++step) {
        const std::size_t index = lowering_index_rec(
            labels, state, labels.size()
        );
        ++state[index];
    }
}

void enumerate_highest_rec(
    const std::vector<int>& labels,
    int target,
    std::size_t index,
    std::vector<int>& state,
    std::vector<std::vector<int>>& output
) {
    if (index == labels.size()) {
        if (crystal_weight(labels, state) == target
            && crystal_stats(labels, state, labels.size()).epsilon == 0) {
            output.push_back(state);
        }
        return;
    }
    for (int down = 0; down <= labels[index]; ++down) {
        state[index] = down;
        enumerate_highest_rec(labels, target, index + 1U, state, output);
    }
}

std::vector<std::vector<int>> highest_states(
    const std::vector<int>& labels,
    int target
) {
    std::vector<std::vector<int>> output;
    std::vector<int> state(labels.size(), 0);
    enumerate_highest_rec(labels, target, 0U, state, output);
    return output;
}

struct Token {
    std::size_t original = 0U;
    int label = 0;
    int down = 0;
};

void swap_by_crystal_r(Token& left, Token& right) {
    const int channel = std::min(right.down, left.label - left.down);
    const int distance = left.down + right.down - channel;
    const int new_left_down = std::min(
        distance, right.label - channel
    );
    const int new_right_down = channel
        + std::max(0, distance - (right.label - channel));
    Token new_left{right.original, right.label, new_left_down};
    Token new_right{left.original, left.label, new_right_down};
    left = new_left;
    right = new_right;
}

std::vector<int> crystal_unshuffle(
    const std::vector<int>& labels,
    unsigned int red_mask,
    const std::vector<int>& red_state,
    const std::vector<int>& blue_lowest_state
) {
    std::vector<Token> tokens;
    std::size_t red_index = 0U;
    std::size_t blue_index = 0U;
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        if (((red_mask >> i) & 1U) != 0U) {
            tokens.push_back({i, labels[i], red_state[red_index]});
            ++red_index;
        }
    }
    for (std::size_t i = 0U; i < labels.size(); ++i) {
        if (((red_mask >> i) & 1U) == 0U) {
            tokens.push_back({i, labels[i], blue_lowest_state[blue_index]});
            ++blue_index;
        }
    }
    bool changed = true;
    while (changed) {
        changed = false;
        for (std::size_t i = 0U; i + 1U < tokens.size(); ++i) {
            if (tokens[i].original > tokens[i + 1U].original) {
                swap_by_crystal_r(tokens[i], tokens[i + 1U]);
                changed = true;
            }
        }
    }
    std::vector<int> output(labels.size(), 0);
    for (const Token& token : tokens) {
        output[token.original] = token.down;
    }
    return output;
}

struct CaseResult {
    std::size_t sources = 0U;
    std::size_t distinct_images = 0U;
    bool channels_valid = true;
};

std::vector<int> merge_coloured_state(
    std::size_t factor_count,
    unsigned int first_mask,
    const std::vector<int>& first_state,
    const std::vector<int>& second_state
) {
    std::vector<int> result(factor_count, 0);
    std::size_t first_index = 0U;
    std::size_t second_index = 0U;
    for (std::size_t i = 0U; i < factor_count; ++i) {
        if (((first_mask >> i) & 1U) != 0U) {
            result[i] = first_state[first_index];
            ++first_index;
        } else {
            result[i] = second_state[second_index];
            ++second_index;
        }
    }
    return result;
}

std::map<std::vector<int>, std::size_t> coloured_highest_profile(
    const std::vector<int>& labels,
    int first_target,
    int second_target
) {
    const unsigned int subset_count = 1U << labels.size();
    std::map<std::vector<int>, std::size_t> profile;
    for (unsigned int first_mask = 0U; first_mask < subset_count; ++first_mask) {
        std::vector<int> first_labels;
        std::vector<int> second_labels;
        for (std::size_t i = 0U; i < labels.size(); ++i) {
            if (((first_mask >> i) & 1U) != 0U) {
                first_labels.push_back(labels[i]);
            } else {
                second_labels.push_back(labels[i]);
            }
        }
        const auto first_states = highest_states(first_labels, first_target);
        const auto second_states = highest_states(second_labels, second_target);
        for (const auto& first_state : first_states) {
            for (const auto& second_state : second_states) {
                ++profile[merge_coloured_state(
                    labels.size(), first_mask, first_state, second_state
                )];
            }
        }
    }
    return profile;
}

bool analyze_recolouring(
    const std::vector<int>& labels,
    int target_a,
    int target_b,
    bool verbose
) {
    const auto source = coloured_highest_profile(labels, target_a, target_b);
    std::map<std::vector<int>, std::size_t> boundary;
    for (int channel = std::abs(target_a - target_b);
         channel <= target_a + target_b; channel += 2) {
        const auto profile = coloured_highest_profile(labels, channel, 0);
        for (const auto& [state, count] : profile) {
            boundary[state] += count;
        }
    }
    std::size_t source_total = 0U;
    std::size_t boundary_total = 0U;
    for (const auto& [state, count] : source) {
        source_total += count;
        const auto found = boundary.find(state);
        const std::size_t supply = found == boundary.end() ? 0U : found->second;
        if (count > supply) {
            if (verbose) {
                std::cout << "recolouring_deficit state=";
                for (int value : state) {
                    std::cout << value << ',';
                }
                std::cout << " demand=" << count << " supply=" << supply
                          << '\n';
            }
            return false;
        }
    }
    for (const auto& [state, count] : boundary) {
        (void)state;
        boundary_total += count;
    }
    if (verbose) {
        std::cout << "recolouring_source=" << source_total
                  << " recolouring_boundary=" << boundary_total << '\n';
    }
    return true;
}

CaseResult analyze_case(
    const std::vector<int>& labels,
    int target_a,
    int target_b,
    bool verbose
) {
    if (labels.size() >= std::numeric_limits<unsigned int>::digits) {
        throw std::runtime_error("too many factors for subset mask");
    }
    const unsigned int subset_count = 1U << labels.size();
    std::set<std::vector<int>> images;
    std::map<std::vector<int>, std::vector<unsigned int>> preimages;
    std::size_t sources = 0U;
    bool channels_valid = true;
    for (unsigned int red_mask = 0U; red_mask < subset_count; ++red_mask) {
        std::vector<int> red_labels;
        std::vector<int> blue_labels;
        for (std::size_t i = 0U; i < labels.size(); ++i) {
            if (((red_mask >> i) & 1U) != 0U) {
                red_labels.push_back(labels[i]);
            } else {
                blue_labels.push_back(labels[i]);
            }
        }
        const auto red_states = highest_states(red_labels, target_a);
        const auto blue_states = highest_states(blue_labels, target_b);
        for (const auto& red_state : red_states) {
            for (const auto& blue_state : blue_states) {
                std::vector<int> blue_lowest = blue_state;
                lower_word(blue_labels, blue_lowest, target_b);
                const std::vector<int> image = crystal_unshuffle(
                    labels, red_mask, red_state, blue_lowest
                );
                const CrystalStats stats = crystal_stats(
                    labels, image, labels.size()
                );
                const int top_weight = crystal_weight(labels, image)
                    + 2 * stats.epsilon;
                if (top_weight < std::abs(target_a - target_b)
                    || top_weight > target_a + target_b
                    || ((top_weight - target_a - target_b) & 1) != 0) {
                    channels_valid = false;
                }
                images.insert(image);
                preimages[image].push_back(red_mask);
                ++sources;
            }
        }
    }
    if (verbose) {
        for (const auto& [image, masks] : preimages) {
            if (masks.size() <= 1U) {
                continue;
            }
            std::cout << "collision image=";
            for (int value : image) {
                std::cout << value << ',';
            }
            std::cout << " masks=";
            for (unsigned int mask : masks) {
                std::cout << mask << ',';
            }
            std::cout << '\n';
            break;
        }
    }
    return {sources, images.size(), channels_valid};
}

void combinations_rec(
    int maximum_label,
    int remaining,
    int first,
    std::vector<int>& current,
    std::vector<std::vector<int>>& output
) {
    if (remaining == 0) {
        output.push_back(current);
        return;
    }
    for (int label = first; label <= maximum_label; ++label) {
        current.push_back(label);
        combinations_rec(
            maximum_label, remaining - 1, label, current, output
        );
        current.pop_back();
    }
}

std::vector<std::vector<int>> combinations(int maximum_label, int size) {
    std::vector<std::vector<int>> output;
    std::vector<int> current;
    combinations_rec(maximum_label, size, 1, current, output);
    return output;
}

void print_vector(const std::vector<int>& values) {
    std::cout << '[';
    for (std::size_t i = 0U; i < values.size(); ++i) {
        if (i != 0U) {
            std::cout << ',';
        }
        std::cout << values[i];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc >= 5 && std::string(argv[1]) == "--recolour") {
            const int target_a = static_cast<int>(std::strtol(argv[2], nullptr, 10));
            const int target_b = static_cast<int>(std::strtol(argv[3], nullptr, 10));
            std::vector<int> labels;
            for (int i = 4; i < argc; ++i) {
                labels.push_back(
                    static_cast<int>(std::strtol(argv[i], nullptr, 10))
                );
            }
            return analyze_recolouring(
                labels, target_a, target_b, true
            ) ? 0 : 1;
        }
        if (argc >= 5 && std::string(argv[1]) == "--case") {
            const int target_a = static_cast<int>(std::strtol(argv[2], nullptr, 10));
            const int target_b = static_cast<int>(std::strtol(argv[3], nullptr, 10));
            std::vector<int> labels;
            for (int i = 4; i < argc; ++i) {
                labels.push_back(
                    static_cast<int>(std::strtol(argv[i], nullptr, 10))
                );
            }
            const CaseResult result = analyze_case(
                labels, target_a, target_b, true
            );
            std::cout << "sources=" << result.sources
                      << " distinct_images=" << result.distinct_images
                      << " channels_valid=" << result.channels_valid << '\n';
            return result.sources == result.distinct_images
                    && result.channels_valid
                ? 0 : 1;
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: search_su2_crystal_shuffle MAX_LABEL MAX_FACTORS"
                " | --case A B LABEL... | --recolour A B LABEL..."
            );
        }
        const int maximum_label = static_cast<int>(
            std::strtol(argv[1], nullptr, 10)
        );
        const int maximum_factors = static_cast<int>(
            std::strtol(argv[2], nullptr, 10)
        );
        if (maximum_label <= 0 || maximum_factors < 0
            || maximum_factors >= std::numeric_limits<unsigned int>::digits) {
            throw std::runtime_error("invalid search bound");
        }
        std::size_t tested = 0U;
        for (int factor_count = 0; factor_count <= maximum_factors;
             ++factor_count) {
            for (const auto& labels : combinations(maximum_label, factor_count)) {
                int total = 0;
                for (int label : labels) {
                    total += label;
                }
                for (int a = 1; a <= total; ++a) {
                    if (std::binary_search(labels.begin(), labels.end(), a)) {
                        continue;
                    }
                    for (int b = 1; b <= a; ++b) {
                        if (std::binary_search(labels.begin(), labels.end(), b)) {
                            continue;
                        }
                        const CaseResult result = analyze_case(
                            labels, a, b, false
                        );
                        ++tested;
                        if (!result.channels_valid
                            || result.sources != result.distinct_images) {
                            std::cout << "CRYSTAL_SHUFFLE_COUNTEREXAMPLE labels=";
                            print_vector(labels);
                            std::cout << " target=(" << a << ',' << b << ')'
                                      << " sources=" << result.sources
                                      << " distinct_images="
                                      << result.distinct_images
                                      << " channels_valid="
                                      << result.channels_valid << '\n';
                            (void)analyze_case(labels, a, b, true);
                            return 1;
                        }
                    }
                }
            }
        }
        std::cout << "SU2_CRYSTAL_SHUFFLE tested=" << tested
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors << " PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_CRYSTAL_SHUFFLE FAILURE: " << error.what() << '\n';
        return 1;
    }
}
