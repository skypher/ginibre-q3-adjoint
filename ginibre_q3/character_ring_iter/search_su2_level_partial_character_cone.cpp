#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

template <class Function>
void for_each_output(int level, int left, int right, Function function) {
    const int lower = std::abs(left - right);
    const int upper = std::min(left + right, 2 * level - left - right);
    for (int output = lower; output <= upper; output += 2) {
        function(output);
    }
}

std::vector<std::vector<int>> multisets(int level, int maximum_size) {
    std::vector<std::vector<int>> result;
    std::vector<int> current;
    const auto visit = [&](const auto& self, int first) -> void {
        if (!current.empty()) {
            result.push_back(current);
        }
        if (current.size() == static_cast<std::size_t>(maximum_size)) {
            return;
        }
        for (int label = first; label <= level; ++label) {
            current.push_back(label);
            self(self, label);
            current.pop_back();
        }
    };
    visit(visit, 1);
    return result;
}

std::vector<cpp_int> signed_product(
    int level,
    const std::vector<int>& labels,
    std::uint64_t minus_mask
) {
    const int width = level + 1;
    std::vector<cpp_int> current(
        static_cast<std::size_t>(width * width)
    );
    current[0] = 1;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        const int label = labels[index];
        const bool minus = ((minus_mask >> index) & 1U) != 0U;
        std::vector<cpp_int> next(current.size());
        for (int left = 0; left <= level; ++left) {
            for (int right = 0; right <= level; ++right) {
                const cpp_int& value = current[static_cast<std::size_t>(
                    left * width + right
                )];
                if (value == 0) {
                    continue;
                }
                for_each_output(level, left, label, [&](int output) {
                    next[static_cast<std::size_t>(
                        output * width + right
                    )] += value;
                });
                for_each_output(level, right, label, [&](int output) {
                    next[static_cast<std::size_t>(
                        left * width + output
                    )] += minus ? -value : value;
                });
            }
        }
        current = std::move(next);
    }
    return current;
}

std::vector<cpp_int> fusion_hankel_defect(
    int level,
    const std::vector<cpp_int>& coefficient
) {
    const int width = level + 1;
    if (coefficient.size()
        != static_cast<std::size_t>(width * width)) {
        throw std::runtime_error("coefficient matrix has wrong size");
    }
    std::vector<cpp_int> defect(coefficient.size(), 0);
    for (int left = 0; left <= level; ++left) {
        for (int right = 0; right <= level; ++right) {
            cpp_int hankel = 0;
            for_each_output(level, left, right, [&](int channel) {
                hankel += coefficient[static_cast<std::size_t>(
                    channel * width
                )];
            });
            const std::size_t entry = static_cast<std::size_t>(
                left * width + right
            );
            defect[entry] = hankel - coefficient[entry];
        }
    }
    return defect;
}

bool inspect_word(
    int level,
    const std::vector<int>& labels,
    std::uint64_t minus_mask,
    int& failed_target,
    cpp_int& failed_value
) {
    const int width = level + 1;
    const std::vector<cpp_int> current
        = signed_product(level, labels, minus_mask);
    for (int target = 0; target <= level; ++target) {
        const cpp_int& value = current[static_cast<std::size_t>(
            target * width
        )];
        if (value < 0) {
            failed_target = target;
            failed_value = value;
            return false;
        }
    }
    return true;
}

int popcount(std::uint64_t value) {
    int answer = 0;
    while (value != 0U) {
        answer += static_cast<int>(value & 1U);
        value >>= 1U;
    }
    return answer;
}

void print_word(const std::vector<int>& labels, std::uint64_t minus_mask) {
    std::cout << '[';
    for (std::size_t index = 0; index < labels.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << (((minus_mask >> index) & 1U) != 0U ? '-' : '+')
                  << labels[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4
            && std::string(argv[1]) == "simple-current-boundary") {
            const int maximum_level = std::stoi(argv[2]);
            const int maximum_factors = std::stoi(argv[3]);
            if (maximum_level < 2 || maximum_factors < 0
                || maximum_factors >= 63) {
                throw std::runtime_error(
                    "invalid simple-current boundary bounds"
                );
            }
            std::uint64_t words_tested = 0U;
            std::uint64_t entries_tested = 0U;
            std::uint64_t boundary_only_failures = 0U;
            std::uint64_t lower_only_failures = 0U;
            std::uint64_t strict_boundary_failures = 0U;
            std::uint64_t strict_joint_failures = 0U;
            std::uint64_t fully_disjoint_joint_failures = 0U;
            bool printed_boundary_failure = false;
            bool printed_lower_failure = false;
            bool printed_strict_boundary_failure = false;
            bool printed_strict_joint_failure = false;
            bool printed_fully_disjoint_joint_failure = false;
            cpp_int minimum_combined = 0;
            bool have_combined = false;
            for (int level = 2; level <= maximum_level; ++level) {
                auto words = multisets(level - 1, maximum_factors);
                words.insert(words.begin(), std::vector<int>{});
                for (const auto& labels : words) {
                    const int width = level + 1;
                    const std::vector<cpp_int> coefficient
                        = signed_product(level, labels, 0U);
                    const std::vector<cpp_int> defect
                        = fusion_hankel_defect(level, coefficient);
                    ++words_tested;
                    const int first_next
                        = labels.empty() ? 1 : labels.back();
                    for (int next_label = first_next;
                         next_label < level; ++next_label) {
                        const bool strict_next = labels.empty()
                            || next_label > labels.back();
                        std::vector<int> extended = labels;
                        extended.push_back(next_label);
                        const std::vector<cpp_int> next_defect
                            = fusion_hankel_defect(
                                level,
                                signed_product(level, extended, 0U)
                        );
                        for (int target = 0; target <= level; ++target) {
                            if (std::binary_search(
                                    extended.begin(), extended.end(), target
                                )) {
                                continue;
                            }
                            cpp_int boundary = 0;
                            for_each_output(
                                level, target, next_label,
                                [&](int output) {
                                    boundary += defect[
                                        static_cast<std::size_t>(
                                            output * width + level
                                        )
                                    ];
                                }
                            );
                            const cpp_int& lower = defect[
                                static_cast<std::size_t>(
                                    target * width + level - next_label
                                )
                            ];
                            const cpp_int& sink = defect[
                                static_cast<std::size_t>(
                                    (level - target) * width + next_label
                                )
                            ];
                            const cpp_int combined = boundary + lower - sink;
                            const cpp_int& direct = next_defect[
                                static_cast<std::size_t>(
                                    target * width + level
                                )
                            ];
                            ++entries_tested;
                            if (combined != direct || direct < 0) {
                                std::cout
                                    << "SU2_SIMPLE_CURRENT_BOUNDARY "
                                    << "result=FAIL level=" << level
                                    << " prefix=";
                                print_word(labels, 0U);
                                std::cout << " next=" << next_label
                                          << " target=" << target
                                          << " boundary=" << boundary
                                          << " lower=" << lower
                                          << " sink=" << sink
                                          << " combined=" << combined
                                          << " direct=" << direct << '\n';
                                return EXIT_FAILURE;
                            }
                            if (!have_combined || combined < minimum_combined) {
                                minimum_combined = combined;
                                have_combined = true;
                            }
                            if (boundary < sink) {
                                ++boundary_only_failures;
                                if (!printed_boundary_failure) {
                                    std::cout
                                        << "boundary_only_failure level="
                                        << level << " prefix=";
                                    print_word(labels, 0U);
                                    std::cout << " next=" << next_label
                                              << " target=" << target
                                              << " boundary=" << boundary
                                              << " lower=" << lower
                                              << " sink=" << sink << '\n';
                                    printed_boundary_failure = true;
                                }
                                if (strict_next) {
                                    ++strict_boundary_failures;
                                    if (!printed_strict_boundary_failure) {
                                        std::cout
                                            << "strict_boundary_only_failure"
                                            << " level=" << level
                                            << " prefix=";
                                        print_word(labels, 0U);
                                        std::cout << " next=" << next_label
                                                  << " target=" << target
                                                  << " boundary=" << boundary
                                                  << " lower=" << lower
                                                  << " sink=" << sink << '\n';
                                        printed_strict_boundary_failure = true;
                                    }
                                }
                            }
                            if (lower < sink) {
                                ++lower_only_failures;
                                if (!printed_lower_failure) {
                                    std::cout
                                        << "lower_only_failure level="
                                        << level << " prefix=";
                                    print_word(labels, 0U);
                                    std::cout << " next=" << next_label
                                              << " target=" << target
                                              << " boundary=" << boundary
                                              << " lower=" << lower
                                              << " sink=" << sink << '\n';
                                    printed_lower_failure = true;
                                }
                            }
                            if (strict_next && boundary < sink
                                && lower < sink) {
                                ++strict_joint_failures;
                                if (!printed_strict_joint_failure) {
                                    std::cout << "strict_joint_failure"
                                              << " level=" << level
                                              << " prefix=";
                                    print_word(labels, 0U);
                                    std::cout << " next=" << next_label
                                              << " target=" << target
                                              << " boundary=" << boundary
                                              << " lower=" << lower
                                              << " sink=" << sink << '\n';
                                    printed_strict_joint_failure = true;
                                }
                                if (!std::binary_search(
                                        labels.begin(), labels.end(),
                                        level - target
                                    )) {
                                    ++fully_disjoint_joint_failures;
                                    if (!printed_fully_disjoint_joint_failure) {
                                        std::cout
                                            << "fully_disjoint_joint_failure"
                                            << " level=" << level
                                            << " prefix=";
                                        print_word(labels, 0U);
                                        std::cout << " next=" << next_label
                                                  << " target=" << target
                                                  << " boundary=" << boundary
                                                  << " lower=" << lower
                                                  << " sink=" << sink << '\n';
                                        printed_fully_disjoint_joint_failure
                                            = true;
                                    }
                                }
                            }
                        }
                    }
                }
                std::cout << "progress level=" << level
                          << " words=" << words_tested
                          << " entries=" << entries_tested
                          << " boundary_only_failures="
                          << boundary_only_failures
                          << " lower_only_failures="
                          << lower_only_failures
                          << " strict_boundary_only_failures="
                          << strict_boundary_failures
                          << " strict_joint_failures="
                          << strict_joint_failures
                          << " fully_disjoint_joint_failures="
                          << fully_disjoint_joint_failures << '\n';
            }
            std::cout << "SU2_SIMPLE_CURRENT_BOUNDARY levels=2.."
                      << maximum_level
                      << " maximum_prefix_factors=" << maximum_factors
                      << " words=" << words_tested
                      << " entries=" << entries_tested
                      << " boundary_only_failures="
                      << boundary_only_failures
                      << " lower_only_failures=" << lower_only_failures
                      << " strict_boundary_only_failures="
                      << strict_boundary_failures
                      << " strict_joint_failures="
                      << strict_joint_failures
                      << " fully_disjoint_joint_failures="
                      << fully_disjoint_joint_failures
                      << " minimum_combined="
                      << (have_combined ? minimum_combined : cpp_int{0})
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc == 4
            && std::string(argv[1]) == "simple-current-wall") {
            const int maximum_level = std::stoi(argv[2]);
            const int maximum_factors = std::stoi(argv[3]);
            if (maximum_level < 2 || maximum_factors < 0
                || maximum_factors >= 63) {
                throw std::runtime_error(
                    "invalid simple-current wall bounds"
                );
            }
            std::uint64_t total_tested = 0U;
            std::uint64_t total_active = 0U;
            std::uint64_t total_tight = 0U;
            std::uint64_t printed_tight = 0U;
            cpp_int minimum_slack = 0;
            bool have_slack = false;
            const std::size_t sector_count
                = static_cast<std::size_t>(maximum_factors + 1);
            std::vector<std::uint64_t> sector_tested(sector_count, 0U);
            std::vector<std::uint64_t> sector_active(sector_count, 0U);
            std::vector<std::uint64_t> sector_tight(sector_count, 0U);
            std::vector<cpp_int> sector_minimum(sector_count, 0);
            std::vector<bool> sector_has_minimum(sector_count, false);
            for (int level = 2; level <= maximum_level; ++level) {
                auto words = multisets(level - 1, maximum_factors);
                words.insert(words.begin(), std::vector<int>{});
                for (const auto& labels : words) {
                    const std::uint64_t masks
                        = std::uint64_t{1} << labels.size();
                    for (std::uint64_t mask = 0U; mask < masks; ++mask) {
                        if ((popcount(mask) & 1) != 0) {
                            continue;
                        }
                        const int width = level + 1;
                        const std::vector<cpp_int> product
                            = signed_product(level, labels, mask);
                        const cpp_int& h_00 = product[0];
                        const cpp_int& h_kk = product[
                            static_cast<std::size_t>(
                                level * width + level
                            )
                        ];
                        const cpp_int slack = h_00 - h_kk;
                        const std::size_t minus_count
                            = static_cast<std::size_t>(popcount(mask));
                        ++total_tested;
                        ++sector_tested[minus_count];
                        if (h_00 < 0 || slack < 0) {
                            std::cout
                                << "SU2_SIMPLE_CURRENT_WALL result=FAIL"
                                << " level=" << level << " remainder=";
                            print_word(labels, mask);
                            std::cout << " h_00=" << h_00
                                      << " h_kk=" << h_kk
                                      << " slack=" << slack << '\n';
                            return EXIT_FAILURE;
                        }
                        if (h_kk == 0) {
                            continue;
                        }
                        ++total_active;
                        ++sector_active[minus_count];
                        if (!have_slack || slack < minimum_slack) {
                            minimum_slack = slack;
                            have_slack = true;
                        }
                        if (!sector_has_minimum[minus_count]
                            || slack < sector_minimum[minus_count]) {
                            sector_minimum[minus_count] = slack;
                            sector_has_minimum[minus_count] = true;
                        }
                        if (slack == 0) {
                            ++total_tight;
                            ++sector_tight[minus_count];
                            if (printed_tight < 12U) {
                                std::cout << "tight_simple_current_wall"
                                          << " level=" << level
                                          << " remainder=";
                                print_word(labels, mask);
                                std::cout << " h_00=" << h_00
                                          << " h_kk=" << h_kk << '\n';
                                ++printed_tight;
                            }
                        }
                    }
                }
                std::cout << "progress level=" << level
                          << " tested=" << total_tested
                          << " active=" << total_active
                          << " tight=" << total_tight << '\n';
            }
            for (std::size_t minus_count = 0U;
                 minus_count < sector_count; ++minus_count) {
                if (sector_tested[minus_count] == 0U) {
                    continue;
                }
                std::cout << "minus_sector count=" << minus_count
                          << " tested=" << sector_tested[minus_count]
                          << " active=" << sector_active[minus_count]
                          << " tight=" << sector_tight[minus_count]
                          << " minimum_active_slack="
                          << (sector_has_minimum[minus_count]
                                  ? sector_minimum[minus_count]
                                  : cpp_int{0})
                          << '\n';
            }
            std::cout << "SU2_SIMPLE_CURRENT_WALL levels=2.."
                      << maximum_level
                      << " maximum_remainder_factors=" << maximum_factors
                      << " tested=" << total_tested
                      << " active=" << total_active
                      << " tight=" << total_tight
                      << " minimum_active_slack="
                      << (have_slack ? minimum_slack : cpp_int{0})
                      << " result=PASS\n";
            return EXIT_SUCCESS;
        }
        if (argc != 3) {
            throw std::runtime_error(
                "usage: search_su2_level_partial_character_cone "
                "MAXIMUM_LEVEL MAXIMUM_FACTORS\n"
                "   or: search_su2_level_partial_character_cone "
                "simple-current-wall MAXIMUM_LEVEL MAXIMUM_FACTORS\n"
                "   or: search_su2_level_partial_character_cone "
                "simple-current-boundary MAXIMUM_LEVEL MAXIMUM_FACTORS"
            );
        }
        const int maximum_level = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_level < 1 || maximum_factors < 1
            || maximum_factors >= 63) {
            throw std::runtime_error("invalid search bound");
        }

        std::uint64_t total_tested = 0;
        for (int level = 1; level <= maximum_level; ++level) {
            const auto words = multisets(level, maximum_factors);
            std::uint64_t level_tested = 0;
            for (const auto& labels : words) {
                const std::uint64_t masks
                    = std::uint64_t{1} << labels.size();
                for (std::uint64_t mask = 0; mask < masks; ++mask) {
                    int failed_target = -1;
                    cpp_int failed_value = 0;
                    ++level_tested;
                    if (!inspect_word(
                            level, labels, mask,
                            failed_target, failed_value
                        )) {
                        std::cout << "FAIL level=" << level << " factors=";
                        print_word(labels, mask);
                        std::cout << " target=" << failed_target
                                  << " coefficient=" << failed_value << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
            total_tested += level_tested;
            std::cout << "SU2_LEVEL_PARTIAL_CHARACTER_CONE level=" << level
                      << " tested=" << level_tested << " PASS\n";
        }
        std::cout << "SU2_LEVEL_PARTIAL_CHARACTER_CONE total_tested="
                  << total_tested << " maximum_level=" << maximum_level
                  << " maximum_factors=" << maximum_factors << " PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
