#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Labels5 = std::array<int, 5>;

template <class Function>
void fusion_outputs(const int k, const int a, const int b, Function function) {
    const int upper = std::min(a + b, 2 * k - a - b);
    for (int c = std::abs(a - b); c <= upper; c += 2) {
        function(c);
    }
}

std::vector<std::int64_t> product(
    const int k,
    const std::vector<int>& labels
) {
    std::vector<std::int64_t> current(static_cast<std::size_t>(k + 1), 0);
    std::vector<std::int64_t> next(static_cast<std::size_t>(k + 1), 0);
    current[0] = 1;
    for (const int label : labels) {
        std::fill(next.begin(), next.end(), 0);
        for (int source = 0; source <= k; ++source) {
            const auto value = current[static_cast<std::size_t>(source)];
            if (value == 0) {
                continue;
            }
            fusion_outputs(k, source, label, [&](const int output) {
                next[static_cast<std::size_t>(output)] += value;
            });
        }
        current.swap(next);
    }
    return current;
}

std::int64_t multiplicity(
    const int k,
    const std::vector<int>& labels,
    const int output = 0
) {
    return product(k, labels)[static_cast<std::size_t>(output)];
}

struct Profile {
    std::array<int, 6> active{};
    std::array<int, 6> rank{};
    std::array<std::int64_t, 5> supply{};
    std::int64_t global_m0 = 0;
    std::int64_t global_m2 = 0;
    std::int64_t demand = 0;

    auto operator<=>(const Profile&) const = default;
};

struct StructuralProfile {
    std::array<int, 6> active{};
    std::array<int, 6> rank{};

    auto operator<=>(const StructuralProfile&) const = default;
};

StructuralProfile canonicalize(const Profile& profile) {
    bool initialized = false;
    StructuralProfile best;
    for (int swap_rows = 0; swap_rows < 2; ++swap_rows) {
        std::array<int, 3> columns{0, 1, 2};
        do {
            StructuralProfile candidate;
            for (int i = 0; i < 2; ++i) {
                const int source_i = swap_rows == 0 ? i : 1 - i;
                for (int j = 0; j < 3; ++j) {
                    const int source_j = columns[static_cast<std::size_t>(j)];
                    const auto target =
                        static_cast<std::size_t>(3 * i + j);
                    const auto source =
                        static_cast<std::size_t>(3 * source_i + source_j);
                    candidate.active[target] = profile.active[source];
                    candidate.rank[target] = profile.rank[source];
                }
            }
            if (!initialized || candidate < best) {
                best = candidate;
                initialized = true;
            }
        } while (std::next_permutation(columns.begin(), columns.end()));
    }
    return best;
}

struct Representative {
    int k = 0;
    std::array<int, 2> odd{};
    std::array<int, 3> even{};
    std::int64_t margin = 0;
    Profile full{};
};

std::int64_t pi_pair(
    const int k,
    const std::array<int, 2>& pair,
    const std::array<int, 3>& triple
) {
    const auto pair_product = product(k, {pair[0], pair[1]});
    const auto triple_product = product(k, {triple[0], triple[1], triple[2]});
    return pair_product[0]
            * (2 * triple_product[0] + triple_product[2])
        + pair_product[2] * triple_product[0];
}

Profile make_profile(
    const int k,
    const std::array<int, 2>& odd,
    const std::array<int, 3>& even
) {
    Profile result;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            const auto index = static_cast<std::size_t>(3 * i + j);
            result.active[index] =
                std::abs(odd[static_cast<std::size_t>(i)]
                         - even[static_cast<std::size_t>(j)]) == 1
                ? 1
                : 0;
            std::vector<int> complement{
                odd[static_cast<std::size_t>(1 - i)]
            };
            for (int h = 0; h < 3; ++h) {
                if (h != j) {
                    complement.push_back(even[static_cast<std::size_t>(h)]);
                }
            }
            result.rank[index] = static_cast<int>(
                multiplicity(k, complement, 1)
            );
            if (result.active[index] == 0) {
                result.rank[index] = 0;
            } else {
                result.demand += 2 * result.rank[index];
            }
        }
    }

    const std::vector<int> all{
        odd[0], odd[1], even[0], even[1], even[2]
    };
    result.global_m0 = multiplicity(k, all, 0);
    result.global_m2 = multiplicity(k, all, 2);
    result.supply[0] =
        2 * result.global_m0 + result.global_m2;
    for (int j = 0; j < 3; ++j) {
        if (even[static_cast<std::size_t>(j)] != 2) {
            continue;
        }
        result.supply[0] += multiplicity(
            k,
            {
                odd[0], odd[1],
                even[static_cast<std::size_t>((j + 1) % 3)],
                even[static_cast<std::size_t>((j + 2) % 3)]
            }
        );
    }
    result.supply[1] = pi_pair(
        k,
        {odd[0], odd[1]},
        {even[0], even[1], even[2]}
    );
    for (int omitted = 0; omitted < 3; ++omitted) {
        std::array<int, 2> pair{};
        int at = 0;
        for (int j = 0; j < 3; ++j) {
            if (j != omitted) {
                pair[static_cast<std::size_t>(at++)] =
                    even[static_cast<std::size_t>(j)];
            }
        }
        result.supply[static_cast<std::size_t>(2 + omitted)] = pi_pair(
            k,
            pair,
            {odd[0], odd[1], even[static_cast<std::size_t>(omitted)]}
        );
    }
    return result;
}

void print_array(const auto& values) {
    std::cout << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        std::cout << (i == 0U ? "" : ",") << values[i];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 6 && std::string(argv[1]) == "--w-table") {
            const int k = std::stoi(argv[2]);
            const int f = std::stoi(argv[3]);
            const int u = std::stoi(argv[4]);
            const int v = std::stoi(argv[5]);
            const auto coefficient = [&](const int target) {
                std::int64_t value = 0;
                for (int x = 1; x <= u; x += 2) {
                    for (int y = 1; y <= v; y += 2) {
                        fusion_outputs(k, x, y, [&](const int output) {
                            if (output == target) {
                                ++value;
                            }
                        });
                    }
                }
                return value;
            };
            const auto c_minus = f >= 2 ? coefficient(f - 2) : 0;
            const auto c_zero = coefficient(f);
            const auto c_plus = f + 2 <= k ? coefficient(f + 2) : 0;
            std::cout << "C=(" << c_minus << ',' << c_zero << ','
                      << c_plus << ") S="
                      << c_minus + 3 * c_zero + c_plus << '\n';
            return EXIT_SUCCESS;
        }
        if (argc == 3 && std::string(argv[1]) == "--cap-grid") {
            const int k = std::stoi(argv[2]);
            if (k < 3) {
                throw std::runtime_error("invalid level");
            }
            std::vector<std::array<int, 2>> odd_multisets;
            std::vector<std::array<int, 3>> even_multisets;
            for (int o0 = 3; o0 <= k; o0 += 2) {
                for (int o1 = o0; o1 <= k; o1 += 2) {
                    const int maximum_cap = std::max(
                        std::min(o0, k - o0) + 1,
                        std::min(o1, k - o1) + 1
                    );
                    if (maximum_cap <= 4) {
                        odd_multisets.push_back({o0, o1});
                    }
                }
            }
            for (int e0 = 2; e0 <= k; e0 += 2) {
                for (int e1 = e0; e1 <= k; e1 += 2) {
                    for (int e2 = e1; e2 <= k; e2 += 2) {
                        const int maximum_cap = std::max({
                            std::min(e0, k - e0) + 1,
                            std::min(e1, k - e1) + 1,
                            std::min(e2, k - e2) + 1
                        });
                        if (maximum_cap <= 4) {
                            even_multisets.push_back({e0, e1, e2});
                        }
                    }
                }
            }
            std::cout << "k=" << k << " columns";
            for (const auto& even : even_multisets) {
                std::cout << ' ';
                print_array(even);
            }
            std::cout << '\n';
            for (const auto& odd : odd_multisets) {
                std::cout << "row=";
                print_array(odd);
                for (const auto& even : even_multisets) {
                    const auto profile = make_profile(k, odd, even);
                    const int sigma = odd[0] + odd[1]
                        + even[0] + even[1] + even[2];
                    std::cout << ' '
                              << (profile.supply[0] - profile.demand)
                              << (sigma >= 2 * k ? "w" : "s");
                }
                std::cout << '\n';
            }
            return EXIT_SUCCESS;
        }
        if (argc == 3 && std::string(argv[1]) == "--cap-enum") {
            const int k = std::stoi(argv[2]);
            if (k < 3) {
                throw std::runtime_error("invalid level");
            }
            std::uint64_t cap_cases = 0;
            std::uint64_t wall_cases = 0;
            std::uint64_t active_cases = 0;
            std::uint64_t raw_deficits = 0;
            for (int o0 = 3; o0 <= k; o0 += 2) {
                for (int o1 = o0; o1 <= k; o1 += 2) {
                    for (int e0 = 2; e0 <= k; e0 += 2) {
                        for (int e1 = e0; e1 <= k; e1 += 2) {
                            for (int e2 = e1; e2 <= k; e2 += 2) {
                                const int maximum_cap = std::max({
                                    std::min(o0, k - o0) + 1,
                                    std::min(o1, k - o1) + 1,
                                    std::min(e0, k - e0) + 1,
                                    std::min(e1, k - e1) + 1,
                                    std::min(e2, k - e2) + 1
                                });
                                if (maximum_cap > 4) {
                                    continue;
                                }
                                ++cap_cases;
                                if (o0 + o1 + e0 + e1 + e2 < 2 * k) {
                                    continue;
                                }
                                ++wall_cases;
                                const auto profile = make_profile(
                                    k, {o0, o1}, {e0, e1, e2}
                                );
                                if (profile.demand == 0) {
                                    continue;
                                }
                                ++active_cases;
                                if (profile.supply[0] >= profile.demand) {
                                    continue;
                                }
                                ++raw_deficits;
                                std::cout << "deficit odd=";
                                print_array(std::array<int, 2>{o0, o1});
                                std::cout << " even=";
                                print_array(std::array<int, 3>{e0, e1, e2});
                                std::cout << " supply0=" << profile.supply[0]
                                          << " T=" << profile.demand
                                          << '\n';
                            }
                        }
                    }
                }
            }
            std::cout << "k=" << k << " cap_cases=" << cap_cases
                      << " wall_cases=" << wall_cases
                      << " active_cases=" << active_cases
                      << " raw_deficits=" << raw_deficits << '\n';
            return EXIT_SUCCESS;
        }
        if (argc == 7) {
            const int k = std::stoi(argv[1]);
            const std::array<int, 2> odd{
                std::stoi(argv[2]), std::stoi(argv[3])
            };
            const std::array<int, 3> even{
                std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6])
            };
            const auto profile = make_profile(k, odd, even);
            std::cout << "case k=" << k << " odd=";
            print_array(odd);
            std::cout << " even=";
            print_array(even);
            std::cout << " A=";
            print_array(profile.active);
            std::cout << " D=";
            print_array(profile.rank);
            std::cout << " supply=";
            print_array(profile.supply);
            std::cout << " global=(" << profile.global_m0 << ','
                      << profile.global_m2 << ')';
            std::cout << " T=" << profile.demand
                      << " margin="
                      << std::accumulate(
                             profile.supply.begin(),
                             profile.supply.end(),
                             -profile.demand
                         )
                      << '\n';
            return EXIT_SUCCESS;
        }
        const bool print_raw_deficits =
            argc == 4 && std::string(argv[1]) == "--s0-deficit";
        const bool print_cap_table =
            argc == 4 && std::string(argv[1]) == "--cap-table";
        if (argc != 3 && !print_raw_deficits && !print_cap_table) {
            throw std::runtime_error(
                "usage: analyze_su2_q1_two_odd_profiles MIN_LEVEL MAX_LEVEL\n"
                "   or: analyze_su2_q1_two_odd_profiles"
                " --s0-deficit MIN_LEVEL MAX_LEVEL\n"
                "   or: analyze_su2_q1_two_odd_profiles"
                " --cap-table MIN_LEVEL MAX_LEVEL\n"
                "   or: analyze_su2_q1_two_odd_profiles --cap-enum K\n"
                "   or: analyze_su2_q1_two_odd_profiles --cap-grid K\n"
                "   or: analyze_su2_q1_two_odd_profiles"
                " --w-table K F U V\n"
                "   or: analyze_su2_q1_two_odd_profiles"
                " K ODD1 ODD2 EVEN1 EVEN2 EVEN3"
            );
        }
        const int argument_offset =
            (print_raw_deficits || print_cap_table) ? 1 : 0;
        const int minimum_level = std::stoi(argv[1 + argument_offset]);
        const int maximum_level = std::stoi(argv[2 + argument_offset]);
        if (minimum_level < 3 || maximum_level < minimum_level) {
            throw std::runtime_error("invalid level range");
        }

        std::map<StructuralProfile, Representative> profiles;
        std::vector<Representative> raw_deficits;
        std::map<std::array<int, 6>, Representative> cap_table;
        std::uint64_t cases = 0;
        for (int k = minimum_level; k <= maximum_level; ++k) {
            std::map<StructuralProfile, Representative> level_profiles;
            std::vector<Representative> level_raw_deficits;
            std::map<std::array<int, 6>, Representative> level_cap_table;
#pragma omp parallel
            {
                std::map<StructuralProfile, Representative> local_profiles;
                std::vector<Representative> local_raw_deficits;
                std::map<std::array<int, 6>, Representative> local_cap_table;
#pragma omp for schedule(dynamic)
                for (int o0 = 3; o0 <= k; o0 += 2) {
                    for (int o1 = o0; o1 <= k; o1 += 2) {
                        for (int e0 = 2; e0 <= k; e0 += 2) {
                            for (int e1 = e0; e1 <= k; e1 += 2) {
                                for (int e2 = e1; e2 <= k; e2 += 2) {
                                    const std::array<int, 2> odd{o0, o1};
                                    const std::array<int, 3> even{e0, e1, e2};
                                    if (o0 + o1 + e0 + e1 + e2 < 2 * k) {
                                        continue;
                                    }
                                    const auto profile =
                                        make_profile(k, odd, even);
                                    if (profile.demand == 0) {
                                        continue;
                                    }
                                    const auto supply = std::accumulate(
                                        profile.supply.begin(),
                                        profile.supply.end(),
                                        std::int64_t{0}
                                    );
                                    const Representative representative{
                                        k, odd, even,
                                        supply - profile.demand, profile
                                    };
                                    if (print_raw_deficits
                                        && profile.supply[0]
                                            < profile.demand) {
                                        local_raw_deficits.push_back(
                                            representative
                                        );
                                    }
                                    const auto structural =
                                        canonicalize(profile);
                                    const int maximum_cap = std::max({
                                        std::min(o0, k - o0) + 1,
                                        std::min(o1, k - o1) + 1,
                                        std::min(e0, k - e0) + 1,
                                        std::min(e1, k - e1) + 1,
                                        std::min(e2, k - e2) + 1
                                    });
                                    if (print_cap_table
                                        && maximum_cap >= 5) {
                                        const std::int64_t raw_margin =
                                            profile.supply[0]
                                            - profile.demand;
                                        Representative raw_representative{
                                            k, odd, even, raw_margin, profile
                                        };
                                        auto [raw_position, raw_inserted] =
                                            local_cap_table.emplace(
                                                structural.active,
                                                raw_representative
                                            );
                                        if (!raw_inserted
                                            && raw_margin
                                                < raw_position
                                                      ->second.margin) {
                                            raw_position->second =
                                                raw_representative;
                                        }
                                    }
                                    auto [position, inserted] =
                                        local_profiles.emplace(
                                            structural, representative
                                        );
                                    if (!inserted
                                        && representative.margin
                                            < position->second.margin) {
                                        position->second = representative;
                                    }
                                }
                            }
                        }
                    }
                }
#pragma omp critical
                {
                    for (const auto& [structural, representative]
                         : local_profiles) {
                        auto [position, inserted] = level_profiles.emplace(
                            structural, representative
                        );
                        if (!inserted
                            && representative.margin
                                < position->second.margin) {
                            position->second = representative;
                        }
                    }
                    level_raw_deficits.insert(
                        level_raw_deficits.end(),
                        local_raw_deficits.begin(),
                        local_raw_deficits.end()
                    );
                    for (const auto& [active, representative]
                         : local_cap_table) {
                        auto [position, inserted] =
                            level_cap_table.emplace(active, representative);
                        if (!inserted
                            && representative.margin
                                < position->second.margin) {
                            position->second = representative;
                        }
                    }
                }
            }
            raw_deficits.insert(
                raw_deficits.end(),
                level_raw_deficits.begin(),
                level_raw_deficits.end()
            );
            for (const auto& [active, representative] : level_cap_table) {
                auto [position, inserted] =
                    cap_table.emplace(active, representative);
                if (!inserted
                    && representative.margin < position->second.margin) {
                    position->second = representative;
                }
            }
            for (const auto& [structural, representative] : level_profiles) {
                auto [position, inserted] =
                    profiles.emplace(structural, representative);
                if (!inserted
                    && representative.margin < position->second.margin) {
                    position->second = representative;
                }
            }
            cases += static_cast<std::uint64_t>(level_profiles.size());
            std::cerr << "finished k=" << k
                      << " level_profiles=" << level_profiles.size()
                      << " total_profiles=" << profiles.size() << '\n';
        }

        if (print_raw_deficits) {
            std::sort(
                raw_deficits.begin(),
                raw_deficits.end(),
                [](const Representative& left, const Representative& right) {
                    return std::tie(left.k, left.odd, left.even)
                        < std::tie(right.k, right.odd, right.even);
                }
            );
            std::cout << "raw_deficits=" << raw_deficits.size() << '\n';
            for (const auto& representative : raw_deficits) {
                const auto& profile = representative.full;
                std::cout << "case k=" << representative.k << " odd=";
                print_array(representative.odd);
                std::cout << " even=";
                print_array(representative.even);
                std::cout << " A=";
                print_array(profile.active);
                std::cout << " D=";
                print_array(profile.rank);
                std::cout << " supply=";
                print_array(profile.supply);
                std::cout << " global=(" << profile.global_m0 << ','
                          << profile.global_m2 << ')'
                          << " T=" << profile.demand
                          << " margin=" << representative.margin << '\n';
            }
            return EXIT_SUCCESS;
        }
        if (print_cap_table) {
            std::cout << "cap_table=" << cap_table.size() << '\n';
            for (const auto& [active, representative] : cap_table) {
                const auto& profile = representative.full;
                std::cout << "A=";
                print_array(active);
                std::cout << " min_raw_margin=" << representative.margin
                          << " witness k=" << representative.k << " odd=";
                print_array(representative.odd);
                std::cout << " even=";
                print_array(representative.even);
                std::cout << " D=";
                print_array(canonicalize(profile).rank);
                std::cout << " supply0=" << profile.supply[0]
                          << " T=" << profile.demand << '\n';
            }
            return EXIT_SUCCESS;
        }

        std::cout << "profiles=" << profiles.size()
                  << " level_profile_sum=" << cases << '\n';
        int index = 0;
        for (const auto& [structural, representative] : profiles) {
            const auto& profile = representative.full;
            std::cout << "profile=" << ++index
                      << " k=" << representative.k << " odd=";
            print_array(representative.odd);
            std::cout << " even=";
            print_array(representative.even);
            std::cout << " A=";
            print_array(structural.active);
            std::cout << " D=";
            print_array(structural.rank);
            std::cout << " supply=";
            print_array(profile.supply);
            std::cout << " global=(" << profile.global_m0 << ','
                      << profile.global_m2 << ')';
            std::cout << " T=" << profile.demand
                      << " margin=" << representative.margin << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
