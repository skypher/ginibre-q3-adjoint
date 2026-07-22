#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace {

int wedge_coefficient(int left, int right, int first, int second) {
    return (left == first && right == second ? 1 : 0)
        - (left == second && right == first ? 1 : 0);
}

int one_suffix_entry(int p, int q, int first, int second) {
    int result = 0;
    for (int output = std::abs(q - (p - 1));
         output <= q + p - 1;
         output += 2) {
        result += wedge_coefficient(output, 0, first, second);
    }
    result += wedge_coefficient(p - 1, q, first, second);
    for (int output = std::abs(q - p);
         output <= q + p;
         output += 2) {
        result -= wedge_coefficient(output, 1, first, second);
    }
    result -= wedge_coefficient(p, q - 1, first, second);
    result -= wedge_coefficient(p, q + 1, first, second);
    return result;
}

int two_suffix_entry(
    int p,
    int q,
    int r,
    int first,
    int second
) {
    int result = 0;
    for (int output = std::abs(first - r);
         output <= first + r;
         output += 2) {
        result += one_suffix_entry(p, q, output, second);
    }
    for (int output = std::abs(second - r);
         output <= second + r;
         output += 2) {
        result += one_suffix_entry(p, q, first, output);
    }
    return result;
}

int three_suffix_boundary(
    int p,
    int q,
    int r,
    int s,
    int target
) {
    int result = two_suffix_entry(p, q, r, target, s);
    for (int output = std::abs(target - s);
         output <= target + s;
         output += 2) {
        result += two_suffix_entry(p, q, r, output, 0);
    }
    return result;
}

int multiplicity(std::vector<int> labels, int target) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    std::vector<int> current(static_cast<std::size_t>(total + 1));
    current[0] = 1;
    int support = 0;
    for (int label : labels) {
        std::vector<int> next(current.size());
        for (int input = 0; input <= support; ++input) {
            for (int output = std::abs(input - label);
                 output <= input + label;
                 output += 2) {
                next[static_cast<std::size_t>(output)]
                    += current[static_cast<std::size_t>(input)];
            }
        }
        support += label;
        current = std::move(next);
    }
    if (target < 0 || target > total) {
        return 0;
    }
    return current[static_cast<std::size_t>(target)];
}

struct NegativeProfile {
    int main = 0;
    int positive = 0;
    int negative = 0;
    int first_mask = 0;
    int second_mask = 0;
    int triple = 0;
};

NegativeProfile negative_profile(
    int p,
    int q,
    int r,
    int s,
    int target
) {
    const std::array<int, 3> suffix{q, r, s};
    NegativeProfile result;
    result.main = multiplicity({p - 1, q, r, s}, target);
    for (int omitted = 0; omitted < 3; ++omitted) {
        std::vector<int> pair;
        for (int index = 0; index < 3; ++index) {
            if (index != omitted) {
                pair.push_back(suffix[static_cast<std::size_t>(index)]);
            }
        }
        const int remaining
            = suffix[static_cast<std::size_t>(omitted)];
        if (remaining == p) {
            std::vector<int> labels{1};
            labels.insert(labels.end(), pair.begin(), pair.end());
            result.positive += multiplicity(labels, target);
        }
        result.positive += multiplicity(pair, 0)
            * multiplicity({remaining, p - 1}, target);
        result.positive += multiplicity(pair, p)
            * multiplicity({remaining, 1}, target);
        if (target == remaining) {
            const int contribution = multiplicity(pair, p - 1);
            result.negative += contribution;
            if (contribution != 0) {
                result.first_mask |= 1 << omitted;
            }
        }
        const int second = multiplicity(pair, 1) * multiplicity(
            {remaining, p}, target
        );
        result.negative += second;
        if (second != 0) {
            result.second_mask |= 1 << omitted;
        }
    }
    if (target == p) {
        result.triple = multiplicity({q, r, s}, 1);
        result.negative += result.triple;
    }
    if (target == p - 1) {
        result.positive += multiplicity({q, r, s}, 0);
    }
    if (target == 1) {
        result.positive += multiplicity({q, r, s}, p);
    }
    return result;
}

int low_corner_candidate_count(
    int p,
    int q,
    int r,
    int s,
    int target
) {
    int result = 0;
    const int first_low = q - p + 1;
    const int first_high = q + p - 1;
    for (int first_index = 0; first_index < 6; ++first_index) {
        const int first = first_low + 2 * first_index;
        if (first > first_high) {
            continue;
        }
        if (((first + r - target - s) & 1) != 0) {
            continue;
        }
        const int second_low = std::max(
            std::abs(first - r), std::abs(target - s)
        );
        const int second_high = std::min(first + r, target + s);
        for (int second_index = 0; second_index < 6; ++second_index) {
            const int second = second_low + 2 * second_index;
            if (second <= second_high) {
                ++result;
            }
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_opd_three_suffix MAXIMUM_LABEL"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        if (maximum_label < 2 || maximum_label > 300) {
            throw std::runtime_error("invalid maximum label");
        }
        long long tested = 0;
        long long negative_interior = 0;
        int minimum = std::numeric_limits<int>::max();
        int minimum_main_slack = std::numeric_limits<int>::max();
        std::tuple<int, int, int, int, int> minimum_case{};
        std::tuple<int, int, int, int, int> minimum_main_case{};
        std::map<int, std::pair<int, std::tuple<int, int, int, int, int>>>
            boundary_by_negative;
        std::map<int, std::pair<int, std::tuple<int, int, int, int, int>>>
            main_by_negative_terms;
        using Pattern = std::tuple<int, int, int>;
        std::map<Pattern,
                 std::pair<int, std::tuple<int, int, int, int, int>>>
            main_by_pattern;
        for (int p = 2; p <= maximum_label; ++p) {
            for (int q = p; q <= maximum_label; ++q) {
                for (int r = q; r <= maximum_label; ++r) {
                    for (int s = r; s <= maximum_label; ++s) {
                        const int maximum_target = p + q + r + s + 2;
                        for (int target = 1;
                             target <= maximum_target;
                             ++target) {
                            const int interior = two_suffix_entry(
                                p, q, r, target, s
                            );
                            int boundary = 0;
                            for (int output = std::abs(target - s);
                                 output <= target + s;
                                 output += 2) {
                                boundary += two_suffix_entry(
                                    p, q, r, output, 0
                                );
                            }
                            const int total = boundary + interior;
                            const NegativeProfile profile = negative_profile(
                                p, q, r, s, target
                            );
                            const int main = profile.main;
                            const int negative = profile.negative;
                            if (main + profile.positive - negative != total) {
                                std::cout << "SUBSET_DECOMPOSITION_FAIL p="
                                          << p << " q=" << q << " r=" << r
                                          << " s=" << s
                                          << " target=" << target
                                          << " main=" << main
                                          << " positive=" << profile.positive
                                          << " negative=" << negative
                                          << " direct=" << total << '\n';
                                return EXIT_FAILURE;
                            }
                            const Pattern profile_key{
                                profile.first_mask,
                                profile.second_mask,
                                profile.triple
                            };
                            auto profile_found
                                = main_by_pattern.find(profile_key);
                            if (profile_found == main_by_pattern.end()
                                || main < profile_found->second.first) {
                                main_by_pattern[profile_key] = {
                                    main,
                                    {p, q, r, s, target}
                                };
                            }
                            const int main_slack = main - negative;
                            const int candidate_count
                                = low_corner_candidate_count(
                                    p, q, r, s, target
                                );
                            const bool exceptional = q == p && r == p
                                && s == p + 1 && target == p
                                && (p & 1) == 0;
                            if (candidate_count < negative
                                && !exceptional) {
                                std::cout << "CANDIDATE_POOL_FAIL p=" << p
                                          << " q=" << q << " r=" << r
                                          << " s=" << s
                                          << " target=" << target
                                          << " candidates="
                                          << candidate_count
                                          << " negative=" << negative
                                          << '\n';
                                return EXIT_FAILURE;
                            }
                            auto main_found = main_by_negative_terms.find(
                                negative
                            );
                            if (main_found == main_by_negative_terms.end()
                                || main < main_found->second.first) {
                                main_by_negative_terms[negative] = {
                                    main,
                                    {p, q, r, s, target}
                                };
                            }
                            if (main_slack < minimum_main_slack) {
                                minimum_main_slack = main_slack;
                                minimum_main_case = {p, q, r, s, target};
                            }
                            if (main_slack < 0
                                && !(q == p && r == p && s == p + 1
                                     && target == p)) {
                                std::cout << "MAIN_SLACK_SUPPORT_FAIL p="
                                          << p << " q=" << q << " r=" << r
                                          << " s=" << s
                                          << " target=" << target
                                          << " slack=" << main_slack << '\n';
                                return EXIT_FAILURE;
                            }
                            ++tested;
                            if (total < minimum) {
                                minimum = total;
                                minimum_case = {p, q, r, s, target};
                            }
                            if (total < 0) {
                                std::cout << "FAIL p=" << p << " q=" << q
                                          << " r=" << r << " s=" << s
                                          << " target=" << target
                                          << " boundary=" << boundary
                                          << " interior=" << interior
                                          << " total=" << total << '\n';
                                return EXIT_FAILURE;
                            }
                            if (interior < 0) {
                                ++negative_interior;
                                auto found = boundary_by_negative.find(
                                    interior
                                );
                                if (found == boundary_by_negative.end()
                                    || boundary < found->second.first) {
                                    boundary_by_negative[interior] = {
                                        boundary,
                                        {p, q, r, s, target}
                                    };
                                }
                            }
                            const int direct = three_suffix_boundary(
                                p, q, r, s, target
                            );
                            if (direct != total) {
                                throw std::runtime_error(
                                    "internal update mismatch"
                                );
                            }
                        }
                    }
                }
            }
        }
        const auto [p, q, r, s, target] = minimum_case;
        const auto [mp, mq, mr, ms, mt] = minimum_main_case;
        std::cout << "SU2_OPD_THREE_SUFFIX PASS tested=" << tested
                  << " maximum_label=" << maximum_label
                  << " minimum=" << minimum << " minimum_case="
                  << p << ',' << q << ',' << r << ',' << s << ','
                  << target << " negative_interior=" << negative_interior
                  << '\n';
        std::cout << "minimum_main_slack=" << minimum_main_slack
                  << " witness=" << mp << ',' << mq << ',' << mr << ','
                  << ms << ',' << mt << '\n';
        for (const auto& [interior, witness] : boundary_by_negative) {
            const auto [wp, wq, wr, ws, wt] = witness.second;
            std::cout << "interior=" << interior
                      << " minimum_boundary=" << witness.first
                      << " witness=" << wp << ',' << wq << ',' << wr << ','
                      << ws << ',' << wt << '\n';
        }
        for (const auto& [negative, witness] : main_by_negative_terms) {
            const auto [wp, wq, wr, ws, wt] = witness.second;
            std::cout << "negative_terms=" << negative
                      << " minimum_main=" << witness.first
                      << " witness=" << wp << ',' << wq << ',' << wr << ','
                      << ws << ',' << wt << '\n';
        }
        for (const auto& [pattern, witness] : main_by_pattern) {
            const auto [first_mask, second_mask, triple] = pattern;
            const auto [wp, wq, wr, ws, wt] = witness.second;
            std::cout << "profile=" << first_mask << ',' << second_mask
                      << ',' << triple << " minimum_main=" << witness.first
                      << " witness=" << wp << ',' << wq << ',' << wr << ','
                      << ws << ',' << wt << '\n';
        }
        std::cout << "exceptional_family";
        for (int even_p = 2; even_p <= maximum_label; even_p += 2) {
            if (even_p + 1 > maximum_label) {
                break;
            }
            const NegativeProfile profile = negative_profile(
                even_p, even_p, even_p, even_p + 1, even_p
            );
            const int interior = two_suffix_entry(
                even_p, even_p, even_p, even_p, even_p + 1
            );
            int boundary = 0;
            for (int output = 1;
                 output <= 2 * even_p + 1;
                 output += 2) {
                boundary += two_suffix_entry(
                    even_p, even_p, even_p, output, 0
                );
            }
            std::cout << ' ' << even_p << ":[" << profile.main << ','
                      << profile.negative << ',' << boundary << ','
                      << interior << ',' << boundary + interior << ']';
        }
        std::cout << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
