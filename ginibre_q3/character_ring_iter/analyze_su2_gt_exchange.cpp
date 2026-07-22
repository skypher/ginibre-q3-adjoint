#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Vector = std::vector<int>;

int parse_nonnegative(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value < 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
}

bool valid_top(const Vector& degrees, const Vector& top) {
    int top_total = 0;
    int degree_total = 0;
    int top_before = 0;
    int bottom_through = 0;
    for (std::size_t i = 0U; i < degrees.size(); ++i) {
        top_total += top[i];
        degree_total += degrees[i];
        bottom_through += degrees[i] - top[i];
        if (bottom_through > top_before) {
            return false;
        }
        top_before += top[i];
    }
    return 2 * top_total == degree_total;
}

void enumerate_tops(
    const Vector& degrees,
    std::size_t vertex,
    Vector& top,
    std::vector<Vector>& output
) {
    if (vertex == degrees.size()) {
        if (valid_top(degrees, top)) {
            output.push_back(top);
        }
        return;
    }
    for (int value = 0; value <= degrees[vertex]; ++value) {
        top[vertex] = value;
        enumerate_tops(degrees, vertex + 1U, top, output);
    }
}

std::vector<Vector> tops(const Vector& degrees) {
    Vector top(degrees.size(), 0);
    std::vector<Vector> output;
    enumerate_tops(degrees, 0U, top, output);
    return output;
}

bool contains(const std::vector<Vector>& bases, const Vector& value) {
    return std::find(bases.begin(), bases.end(), value) != bases.end();
}

bool ordinary_exchange(
    const std::vector<Vector>& bases,
    Vector& first_witness,
    Vector& second_witness,
    std::size_t& coordinate
) {
    for (const Vector& first : bases) {
        for (const Vector& second : bases) {
            for (std::size_t i = 0U; i < first.size(); ++i) {
                if (first[i] <= second[i]) {
                    continue;
                }
                bool found = false;
                for (std::size_t j = 0U; j < first.size(); ++j) {
                    if (first[j] <= second[j]) {
                        Vector exchanged = first;
                        --exchanged[i];
                        ++exchanged[j];
                        if (contains(bases, exchanged)) {
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    first_witness = first;
                    second_witness = second;
                    coordinate = i;
                    return false;
                }
            }
        }
    }
    return true;
}

bool increment(Vector& values, int maximum) {
    for (std::size_t reverse = 0U; reverse < values.size(); ++reverse) {
        const std::size_t i = values.size() - 1U - reverse;
        if (values[i] < maximum) {
            ++values[i];
            std::fill(values.begin() + static_cast<std::ptrdiff_t>(i + 1U),
                      values.end(), 0);
            return true;
        }
    }
    return false;
}

void print_vector(const Vector& values) {
    std::cout << '[';
    for (std::size_t i = 0U; i < values.size(); ++i) {
        if (i != 0U) {
            std::cout << ',';
        }
        std::cout << values[i];
    }
    std::cout << ']';
}

bool admissible_cut(
    const Vector& degrees,
    const Vector& top,
    unsigned int cut
) {
    Vector left(degrees.size(), 0);
    Vector right(degrees.size(), 0);
    for (std::size_t i = 0U; i + 1U < degrees.size(); ++i) {
        if (((cut >> i) & 1U) != 0U) {
            left[i] = degrees[i];
        } else {
            right[i] = degrees[i];
        }
    }
    right.back() = degrees.back();
    Vector left_top(degrees.size(), 0);
    Vector right_top(degrees.size(), 0);
    for (std::size_t i = 0U; i < degrees.size(); ++i) {
        if (left[i] != 0) {
            left_top[i] = top[i];
        }
        if (right[i] != 0) {
            right_top[i] = top[i];
        }
    }
    return valid_top(left, left_top) && valid_top(right, right_top);
}

bool pairwise_xor_capacity(
    const Vector& degrees,
    const std::vector<Vector>& bases,
    bool require_support_disjoint_odd,
    unsigned int& first_cut,
    unsigned int& second_cut,
    std::uint64_t& common_fibres,
    std::uint64_t& xor_capacity
) {
    if (degrees.size() >= std::numeric_limits<unsigned int>::digits) {
        throw std::runtime_error("too many vertices for cut masks");
    }
    const unsigned int cut_limit = 1U << (degrees.size() - 1U);
    std::vector<std::uint64_t> capacities(cut_limit, 0U);
    std::map<std::pair<unsigned int, unsigned int>, std::uint64_t> common;
    for (const Vector& top : bases) {
        std::vector<unsigned int> cuts;
        for (unsigned int cut = 0U; cut < cut_limit; ++cut) {
            if (admissible_cut(degrees, top, cut)) {
                cuts.push_back(cut);
                if (cut != 0U) {
                    ++capacities[cut];
                }
            }
        }
        for (std::size_t i = 0U; i < cuts.size(); ++i) {
            for (std::size_t j = i + 1U; j < cuts.size(); ++j) {
                ++common[{cuts[i], cuts[j]}];
            }
        }
    }
    for (const auto& [cuts, count] : common) {
        if (require_support_disjoint_odd) {
            std::map<int, std::array<unsigned int, 2U>> group_parities;
            for (std::size_t i = 0U; i + 1U < degrees.size(); ++i) {
                if (degrees[i] == 0) {
                    continue;
                }
                group_parities[degrees[i]][0U]
                    ^= (cuts.first >> i) & 1U;
                group_parities[degrees[i]][1U]
                    ^= (cuts.second >> i) & 1U;
            }
            bool first_nonzero = false;
            bool second_nonzero = false;
            for (const auto& [degree, parities] : group_parities) {
                (void)degree;
                first_nonzero = first_nonzero || parities[0U] != 0U;
                second_nonzero = second_nonzero || parities[1U] != 0U;
            }
            if (!first_nonzero || !second_nonzero) {
                continue;
            }
        }
        const unsigned int xor_cut = cuts.first ^ cuts.second;
        if (count > capacities[xor_cut]) {
            first_cut = cuts.first;
            second_cut = cuts.second;
            common_fibres = count;
            xor_capacity = capacities[xor_cut];
            return false;
        }
    }
    return true;
}

bool source_balance(
    const Vector& degrees,
    const std::vector<Vector>& bases,
    unsigned int& minus_groups,
    std::uint64_t& odd_sources,
    std::uint64_t& positive_sources,
    std::uint64_t& empty_sources,
    std::uint64_t& largest_odd_channel
) {
    const unsigned int cut_limit = 1U << (degrees.size() - 1U);
    std::vector<std::uint64_t> capacities(cut_limit, 0U);
    for (const Vector& top : bases) {
        for (unsigned int cut = 0U; cut < cut_limit; ++cut) {
            if (admissible_cut(degrees, top, cut)) {
                ++capacities[cut];
            }
        }
    }
    std::map<int, unsigned int> group_index;
    for (std::size_t i = 0U; i + 1U < degrees.size(); ++i) {
        if (degrees[i] > 0 && !group_index.contains(degrees[i])) {
            const std::size_t next = group_index.size();
            if (next >= std::numeric_limits<unsigned int>::digits) {
                throw std::runtime_error("too many label groups");
            }
            group_index[degrees[i]] = static_cast<unsigned int>(next);
        }
    }
    const unsigned int sign_limit = 1U << group_index.size();
    for (unsigned int signs = 1U; signs < sign_limit; ++signs) {
        unsigned int minus_factor_parity = 0U;
        for (std::size_t i = 0U; i + 1U < degrees.size(); ++i) {
            if (degrees[i] > 0) {
                const unsigned int group = group_index.at(degrees[i]);
                minus_factor_parity ^= (signs >> group) & 1U;
            }
        }
        if (degrees.back() == 0 && minus_factor_parity != 0U) {
            continue;
        }
        std::uint64_t odd = 0U;
        std::uint64_t positive = 0U;
        std::uint64_t largest_odd = 0U;
        for (unsigned int cut = 1U; cut < cut_limit; ++cut) {
            unsigned int parity = 0U;
            for (std::size_t i = 0U; i + 1U < degrees.size(); ++i) {
                if (((cut >> i) & 1U) != 0U && degrees[i] > 0) {
                    const unsigned int group = group_index.at(degrees[i]);
                    parity ^= (signs >> group) & 1U;
                }
            }
            if (parity != 0U) {
                odd += capacities[cut];
                largest_odd = std::max(largest_odd, capacities[cut]);
            } else {
                positive += capacities[cut];
            }
        }
        if (odd > positive && odd - positive > largest_odd) {
            minus_groups = signs;
            odd_sources = odd;
            positive_sources = positive;
            empty_sources = capacities[0U];
            largest_odd_channel = largest_odd;
            return false;
        }
    }
    return true;
}

Vector bender_knuth_swap_top(
    const Vector& degrees,
    const Vector& top,
    std::size_t first
) {
    if (first + 1U >= degrees.size()) {
        throw std::runtime_error("Bender--Knuth swap is out of range");
    }
    int height = 0;
    for (std::size_t i = 0U; i < first; ++i) {
        height += 2 * top[i] - degrees[i];
    }
    const int first_bottom = degrees[first] - top[first];
    const int second_bottom = degrees[first + 1U] - top[first + 1U];
    const int paired = std::max(
        0, first_bottom + second_bottom - height
    );
    Vector swapped = top;
    swapped[first] = top[first + 1U] + paired;
    swapped[first + 1U] = top[first] - paired;
    return swapped;
}

bool bender_knuth_bijection(
    const Vector& degrees,
    const std::vector<Vector>& bases,
    std::size_t& coordinate,
    Vector& witness,
    Vector& image
) {
    for (std::size_t first = 0U; first + 1U < degrees.size(); ++first) {
        Vector swapped_degrees = degrees;
        std::swap(swapped_degrees[first], swapped_degrees[first + 1U]);
        const std::vector<Vector> swapped_bases = tops(swapped_degrees);
        std::set<Vector> images;
        for (const Vector& top : bases) {
            const Vector swapped = bender_knuth_swap_top(
                degrees, top, first
            );
            if (!valid_top(swapped_degrees, swapped)
                || bender_knuth_swap_top(
                    swapped_degrees, swapped, first
                ) != top) {
                coordinate = first;
                witness = top;
                image = swapped;
                return false;
            }
            images.insert(swapped);
        }
        if (images != std::set<Vector>(
                swapped_bases.begin(), swapped_bases.end()
            )) {
            coordinate = first;
            return false;
        }
    }
    return true;
}

void apply_bender_knuth_swap(
    Vector& degrees,
    Vector& top,
    std::size_t first
) {
    top = bender_knuth_swap_top(degrees, top, first);
    std::swap(degrees[first], degrees[first + 1U]);
}

bool bender_knuth_braid(
    const Vector& degrees,
    const std::vector<Vector>& bases,
    std::size_t& coordinate,
    Vector& witness,
    Vector& first_image,
    Vector& second_image
) {
    for (std::size_t first = 0U; first + 2U < degrees.size(); ++first) {
        for (const Vector& top : bases) {
            Vector first_degrees = degrees;
            first_image = top;
            apply_bender_knuth_swap(
                first_degrees, first_image, first
            );
            apply_bender_knuth_swap(
                first_degrees, first_image, first + 1U
            );
            apply_bender_knuth_swap(
                first_degrees, first_image, first
            );

            Vector second_degrees = degrees;
            second_image = top;
            apply_bender_knuth_swap(
                second_degrees, second_image, first + 1U
            );
            apply_bender_knuth_swap(
                second_degrees, second_image, first
            );
            apply_bender_knuth_swap(
                second_degrees, second_image, first + 1U
            );
            if (first_degrees != second_degrees
                || first_image != second_image) {
                coordinate = first;
                witness = top;
                return false;
            }
        }
    }
    return true;
}

using State = std::pair<Vector, Vector>;

Vector bender_knuth_raw_neighbor(
    const Vector& physical_degrees,
    const State& state,
    std::size_t first
) {
    const Vector& order = state.first;
    const Vector& raw_top = state.second;
    int height = 0;
    for (std::size_t position = 0U; position < first; ++position) {
        const std::size_t vertex = static_cast<std::size_t>(order[position]);
        height += 2 * raw_top[vertex] - physical_degrees[vertex];
    }
    const std::size_t first_vertex
        = static_cast<std::size_t>(order[first]);
    const std::size_t second_vertex
        = static_cast<std::size_t>(order[first + 1U]);
    const int first_bottom
        = physical_degrees[first_vertex] - raw_top[first_vertex];
    const int second_bottom
        = physical_degrees[second_vertex] - raw_top[second_vertex];
    const int transferred = std::max(
        0, first_bottom + second_bottom - height
    );
    Vector neighbor = raw_top;
    neighbor[first_vertex] -= transferred;
    neighbor[second_vertex] += transferred;
    return neighbor;
}

bool bender_knuth_connected(
    const Vector& physical_degrees,
    bool classify_by_gcd,
    std::uint64_t& state_count,
    std::uint64_t& component_count,
    State& witness
) {
    Vector order(physical_degrees.size(), 0);
    for (std::size_t i = 0U; i < order.size(); ++i) {
        order[i] = static_cast<int>(i);
    }
    std::set<State> states;
    do {
        Vector ordered_degrees(order.size(), 0);
        for (std::size_t position = 0U;
             position < order.size(); ++position) {
            ordered_degrees[position] = physical_degrees[
                static_cast<std::size_t>(order[position])
            ];
        }
        for (const Vector& ordered_top : tops(ordered_degrees)) {
            Vector raw_top(order.size(), 0);
            for (std::size_t position = 0U;
                 position < order.size(); ++position) {
                raw_top[static_cast<std::size_t>(order[position])]
                    = ordered_top[position];
            }
            states.emplace(order, std::move(raw_top));
        }
    } while (std::next_permutation(order.begin(), order.end()));

    state_count = static_cast<std::uint64_t>(states.size());
    std::set<State> visited;
    std::map<int, std::uint64_t> gcd_components;
    component_count = 0U;
    for (const State& initial : states) {
        if (visited.contains(initial)) {
            continue;
        }
        ++component_count;
        if (component_count == 2U) {
            witness = initial;
        }
        int state_gcd = 0;
        for (int degree : physical_degrees) {
            state_gcd = std::gcd(state_gcd, degree);
        }
        for (int value : initial.second) {
            state_gcd = std::gcd(state_gcd, value);
        }
        ++gcd_components[state_gcd];
        if (gcd_components[state_gcd] == 2U) {
            witness = initial;
        }
        std::queue<State> frontier;
        visited.insert(initial);
        frontier.push(initial);
        while (!frontier.empty()) {
            const State current = frontier.front();
            frontier.pop();
            for (std::size_t i = 0U; i + 1U < current.first.size(); ++i) {
                State neighbor{
                    current.first,
                    bender_knuth_raw_neighbor(physical_degrees, current, i)
                };
                std::swap(neighbor.first[i], neighbor.first[i + 1U]);
                if (!states.contains(neighbor)) {
                    throw std::runtime_error(
                        "Bender--Knuth connectivity neighbor is missing"
                    );
                }
                if (visited.insert(neighbor).second) {
                    frontier.push(std::move(neighbor));
                }
            }
        }
    }
    if (!classify_by_gcd) {
        return component_count == 1U;
    }
    for (const auto& [state_gcd, count] : gcd_components) {
        (void)state_gcd;
        if (count != 1U) {
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5) {
            throw std::runtime_error(
                "usage: analyze_su2_gt_exchange MODE MAXIMUM_VERTICES "
                "MAXIMUM_DEGREE MAXIMUM_TOTAL_DEGREE"
            );
        }
        const std::string mode = argv[1];
        if (mode != "bases" && mode != "pairs" && mode != "bk"
            && mode != "bk-braid"
            && mode != "bk-connected"
            && mode != "bk-gcd-connected"
            && mode != "support-pairs" && mode != "balance"
            && mode != "balance-min2") {
            throw std::runtime_error(
                "mode must be bases, bk, bk-braid, bk-connected, "
                "bk-gcd-connected, pairs, support-pairs, or balance"
            );
        }
        const int maximum_vertices = parse_nonnegative(argv[2], "vertices");
        const int maximum_degree = parse_nonnegative(argv[3], "degree");
        const int maximum_total = parse_nonnegative(argv[4], "total degree");
        std::uint64_t cases = 0U;
        for (int vertices = 1; vertices <= maximum_vertices; ++vertices) {
            Vector degrees(static_cast<std::size_t>(vertices), 0);
            bool more = true;
            while (more) {
                int total = 0;
                for (int degree : degrees) {
                    total += degree;
                }
                if (total <= maximum_total && (total & 1) == 0) {
                    if (mode == "balance-min2"
                        && std::any_of(
                            degrees.begin(), degrees.end(),
                            [](int degree) { return degree < 2; }
                        )) {
                        more = increment(degrees, maximum_degree);
                        continue;
                    }
                    const std::vector<Vector> bases = tops(degrees);
                    if (!bases.empty()) {
                        ++cases;
                        if (mode == "balance" || mode == "balance-min2") {
                            unsigned int minus_groups = 0U;
                            std::uint64_t odd_sources = 0U;
                            std::uint64_t positive_sources = 0U;
                            std::uint64_t empty_sources = 0U;
                            std::uint64_t largest_odd_channel = 0U;
                            if (!source_balance(
                                    degrees, bases, minus_groups, odd_sources,
                                    positive_sources, empty_sources,
                                    largest_odd_channel
                                )) {
                                std::cout
                                    << "SU2_GT_SOURCE_BALANCE result=DEFICIT "
                                    << "degrees=";
                                print_vector(degrees);
                                std::cout << " minus_groups=" << minus_groups
                                          << " odd_sources=" << odd_sources
                                          << " positive_sources="
                                          << positive_sources
                                          << " empty_sources=" << empty_sources
                                          << " largest_odd_channel="
                                          << largest_odd_channel
                                          << " deficit="
                                          << (odd_sources - positive_sources)
                                          << " coefficient="
                                          << (static_cast<std::int64_t>(
                                                  empty_sources
                                              ) + static_cast<std::int64_t>(
                                                  positive_sources
                                              ) - static_cast<std::int64_t>(
                                                  odd_sources
                                              ))
                                          << '\n';
                                return 1;
                            }
                        } else if (mode == "bases") {
                            Vector first;
                            Vector second;
                            std::size_t coordinate = 0U;
                            if (!ordinary_exchange(
                                    bases, first, second, coordinate
                                )) {
                                std::cout
                                    << "SU2_GT_EXCHANGE result=FAIL degrees=";
                                print_vector(degrees);
                                std::cout << " first=";
                                print_vector(first);
                                std::cout << " second=";
                                print_vector(second);
                                std::cout << " coordinate=" << coordinate
                                          << " bases=" << bases.size() << '\n';
                                return 1;
                            }
                        } else if (mode == "bk") {
                            std::size_t coordinate = 0U;
                            Vector witness;
                            Vector image;
                            if (!bender_knuth_bijection(
                                    degrees, bases, coordinate,
                                    witness, image
                                )) {
                                std::cout
                                    << "SU2_GT_BK_BIJECTION result=FAIL "
                                    << "degrees=";
                                print_vector(degrees);
                                std::cout << " coordinate=" << coordinate
                                          << " witness=";
                                print_vector(witness);
                                std::cout << " image=";
                                print_vector(image);
                                std::cout << '\n';
                                return 1;
                            }
                        } else if (mode == "bk-braid") {
                            std::size_t coordinate = 0U;
                            Vector witness;
                            Vector first_image;
                            Vector second_image;
                            if (!bender_knuth_braid(
                                    degrees, bases, coordinate, witness,
                                    first_image, second_image
                                )) {
                                std::cout
                                    << "SU2_GT_BK_BRAID result=FAIL degrees=";
                                print_vector(degrees);
                                std::cout << " coordinate=" << coordinate
                                          << " witness=";
                                print_vector(witness);
                                std::cout << " first_image=";
                                print_vector(first_image);
                                std::cout << " second_image=";
                                print_vector(second_image);
                                std::cout << '\n';
                                return 1;
                            }
                        } else if (mode == "bk-connected"
                                   || mode == "bk-gcd-connected") {
                            std::uint64_t state_count = 0U;
                            std::uint64_t component_count = 0U;
                            State witness;
                            if (!bender_knuth_connected(
                                    degrees, mode == "bk-gcd-connected",
                                    state_count, component_count, witness
                                )) {
                                std::cout
                                    << (mode == "bk-gcd-connected"
                                        ? "SU2_GT_BK_GCD_CONNECTED"
                                        : "SU2_GT_BK_CONNECTED")
                                    << " result=FAIL "
                                    << "degrees=";
                                print_vector(degrees);
                                std::cout << " states=" << state_count
                                          << " components=" << component_count
                                          << " witness_order=";
                                print_vector(witness.first);
                                std::cout << " witness_top=";
                                print_vector(witness.second);
                                std::cout << '\n';
                                return 1;
                            }
                        } else {
                            unsigned int first_cut = 0U;
                            unsigned int second_cut = 0U;
                            std::uint64_t common_fibres = 0U;
                            std::uint64_t xor_capacity = 0U;
                            if (!pairwise_xor_capacity(
                                    degrees, bases, mode == "support-pairs",
                                    first_cut, second_cut,
                                    common_fibres, xor_capacity
                                )) {
                                std::cout
                                    << (mode == "support-pairs"
                                            ? "SU2_GT_SUPPORT_PAIR_XOR"
                                            : "SU2_GT_PAIR_XOR")
                                    << " result=FAIL degrees=";
                                print_vector(degrees);
                                std::cout << " first_cut=" << first_cut
                                          << " second_cut=" << second_cut
                                          << " xor_cut="
                                          << (first_cut ^ second_cut)
                                          << " common_fibres=" << common_fibres
                                          << " xor_capacity=" << xor_capacity
                                          << " bases=" << bases.size() << '\n';
                                for (const Vector& top : bases) {
                                    if (!admissible_cut(
                                            degrees, top, first_cut
                                        ) || !admissible_cut(
                                            degrees, top, second_cut
                                        )) {
                                        continue;
                                    }
                                    std::cout << "witness_top=";
                                    print_vector(top);
                                    std::cout << " admissible_cuts=[";
                                    const unsigned int cut_limit
                                        = 1U << (degrees.size() - 1U);
                                    bool first_printed = false;
                                    for (unsigned int cut = 0U;
                                         cut < cut_limit; ++cut) {
                                        if (admissible_cut(degrees, top, cut)) {
                                            if (first_printed) {
                                                std::cout << ',';
                                            }
                                            std::cout << cut;
                                            first_printed = true;
                                        }
                                    }
                                    std::cout << "]\n";
                                    break;
                                }
                                return 1;
                            }
                        }
                    }
                }
                more = increment(degrees, maximum_degree);
            }
            std::cout << "progress vertices=" << vertices
                      << " cases=" << cases << '\n' << std::flush;
        }
        std::cout << (mode == "bases" ? "SU2_GT_EXCHANGE"
                      : mode == "bk" ? "SU2_GT_BK_BIJECTION"
                      : mode == "bk-braid" ? "SU2_GT_BK_BRAID"
                      : mode == "bk-connected" ? "SU2_GT_BK_CONNECTED"
                      : mode == "bk-gcd-connected"
                          ? "SU2_GT_BK_GCD_CONNECTED"
                      : mode == "balance" || mode == "balance-min2"
                          ? "SU2_GT_SOURCE_BALANCE"
                      : mode == "support-pairs"
                          ? "SU2_GT_SUPPORT_PAIR_XOR"
                          : "SU2_GT_PAIR_XOR")
                  << " cases=" << cases
                  << " result=PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
