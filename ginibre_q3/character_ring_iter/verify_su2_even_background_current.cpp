#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;

struct Fusion {
    explicit Fusion(int fusion_level)
        : level(fusion_level), channels(static_cast<std::size_t>(fusion_level + 1),
              std::vector<std::vector<int>>(
                  static_cast<std::size_t>(fusion_level + 1)
              )) {
        for (int left = 0; left <= fusion_level; ++left) {
            for (int right = 0; right <= fusion_level; ++right) {
                const int upper = std::min(
                    left + right, 2 * fusion_level - left - right
                );
                for (int output = std::abs(left - right);
                     output <= upper; output += 2) {
                    channels[static_cast<std::size_t>(left)]
                            [static_cast<std::size_t>(right)]
                        .push_back(output);
                }
            }
        }
    }

    int level;
    std::vector<std::vector<std::vector<int>>> channels;
};

struct Table {
    explicit Table(int table_dimension)
        : dimension(table_dimension), values(
              static_cast<std::size_t>(table_dimension)
                  * static_cast<std::size_t>(table_dimension)
          ) {}

    Integer& at(int row, int column) {
        return values[static_cast<std::size_t>(row)
                      * static_cast<std::size_t>(dimension)
                      + static_cast<std::size_t>(column)];
    }

    const Integer& at(int row, int column) const {
        return values[static_cast<std::size_t>(row)
                      * static_cast<std::size_t>(dimension)
                      + static_cast<std::size_t>(column)];
    }

    int dimension;
    std::vector<Integer> values;
};

bool next_word(std::vector<int>& word, int maximum_label) {
    for (std::size_t reverse = word.size(); reverse > 0U; --reverse) {
        int& entry = word[reverse - 1U];
        if (entry < maximum_label) {
            ++entry;
            return true;
        }
        entry = 1;
    }
    return false;
}

Table vacuum(const Fusion& fusion) {
    Table result(fusion.level + 1);
    result.at(0, 0) = 1;
    return result;
}

Table update(const Table& source, const Fusion& fusion, int label, int sign) {
    Table result(source.dimension);
    for (int row = 0; row < source.dimension; ++row) {
        for (int column = 0; column < source.dimension; ++column) {
            const Integer& value = source.at(row, column);
            if (value == 0) {
                continue;
            }
            for (const int output : fusion.channels[static_cast<std::size_t>(label)]
                                               [static_cast<std::size_t>(row)]) {
                result.at(output, column) += value;
            }
            for (const int output : fusion.channels[static_cast<std::size_t>(label)]
                                               [static_cast<std::size_t>(column)]) {
                result.at(row, output) += sign * value;
            }
        }
    }
    return result;
}

Table background_table(
    const Fusion& fusion,
    const std::vector<int>& labels,
    unsigned int minus_mask
) {
    Table result = vacuum(fusion);
    for (std::size_t index = 0U; index < labels.size(); ++index) {
        const unsigned int bit = 1U << index;
        const int sign = (minus_mask & bit) == 0U ? 1 : -1;
        result = update(result, fusion, labels[index], sign);
    }
    return result;
}

Table plus_leaves(
    Table result,
    const Fusion& fusion,
    const std::vector<int>& leaves,
    std::size_t first_skip,
    std::size_t second_skip
) {
    for (std::size_t index = 0U; index < leaves.size(); ++index) {
        if (index == first_skip || index == second_skip) {
            continue;
        }
        result = update(result, fusion, leaves[index], 1);
    }
    return result;
}

Integer boundary(const Table& table, const Fusion& fusion, int q, int a) {
    Integer result = -table.at(q, a);
    for (const int output : fusion.channels[static_cast<std::size_t>(q)]
                                             [static_cast<std::size_t>(a)]) {
        result += table.at(0, output);
    }
    return result;
}

Integer contraction(
    const Fusion& fusion,
    const std::vector<int>& labels,
    std::uint64_t minus_mask
) {
    Table result = vacuum(fusion);
    for (std::size_t index = 0U; index < labels.size(); ++index) {
        const std::uint64_t bit = std::uint64_t{1} << index;
        const int sign = (minus_mask & bit) == 0U ? 1 : -1;
        result = update(result, fusion, labels[index], sign);
    }
    return result.at(0, 0);
}

bool even_parity(std::uint64_t mask) {
    int parity = 0;
    while (mask != 0U) {
        parity ^= static_cast<int>(mask & std::uint64_t{1});
        mask >>= 1U;
    }
    return parity == 0;
}

bool verify_two_bit_hypercube(
    const Fusion& fusion,
    const std::vector<int>& labels,
    std::uint64_t& checks
) {
    if (labels.size() < 2U || labels.size() >= 63U) {
        return true;
    }
    const std::uint64_t masks = std::uint64_t{1} << labels.size();
    for (std::uint64_t minus_mask = 0U; minus_mask < masks; ++minus_mask) {
        if (!even_parity(minus_mask)) {
            continue;
        }
        for (std::size_t first = 0U; first < labels.size(); ++first) {
            for (std::size_t second = first + 1U;
                 second < labels.size(); ++second) {
                const std::uint64_t pair =
                    (std::uint64_t{1} << first)
                    | (std::uint64_t{1} << second);
                const Integer lhs = contraction(fusion, labels, minus_mask)
                    + contraction(fusion, labels, minus_mask ^ pair);
                const std::uint64_t pair_signs = minus_mask & pair;
                const bool fused_minus =
                    pair_signs != 0U && pair_signs != pair;
                Integer rhs = 0;
                for (const int fused
                     : fusion.channels[static_cast<std::size_t>(labels[first])]
                                      [static_cast<std::size_t>(labels[second])]) {
                    std::vector<int> child_labels;
                    child_labels.reserve(labels.size() - 1U);
                    std::uint64_t child_mask = 0U;
                    for (std::size_t index = 0U; index < labels.size(); ++index) {
                        if (index == first || index == second) {
                            continue;
                        }
                        const std::size_t child = child_labels.size();
                        child_labels.push_back(labels[index]);
                        if ((minus_mask & (std::uint64_t{1} << index)) != 0U) {
                            child_mask |= std::uint64_t{1} << child;
                        }
                    }
                    const std::size_t star = child_labels.size();
                    child_labels.push_back(fused);
                    if (fused_minus) {
                        child_mask |= std::uint64_t{1} << star;
                    }
                    rhs += 2 * contraction(fusion, child_labels, child_mask);
                }
                if (lhs != rhs) {
                    std::cerr << "TWO_BIT_HYPERCUBE failure"
                              << " level=" << fusion.level
                              << " first=" << first
                              << " second=" << second << '\n';
                    return false;
                }
                ++checks;
            }
        }
    }
    return true;
}

bool symmetric(const Table& table) {
    for (int row = 0; row < table.dimension; ++row) {
        for (int column = row + 1; column < table.dimension; ++column) {
            if (table.at(row, column) != table.at(column, row)) {
                return false;
            }
        }
    }
    return true;
}

int parse_positive(const char* text, const char* name) {
    const int value = std::atoi(text);
    if (value <= 0) {
        std::cerr << "invalid " << name << '\n';
        std::exit(2);
    }
    return value;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: verify_su2_even_background_current "
                  << "MAXIMUM_LEVEL MAXIMUM_LABEL "
                  << "MAXIMUM_BACKGROUND_FACTORS MAXIMUM_LEAVES\n";
        return 2;
    }
    const int maximum_level = parse_positive(argv[1], "maximum level");
    const int maximum_label = parse_positive(argv[2], "maximum label");
    const int maximum_background = parse_positive(
        argv[3], "maximum background factors"
    );
    const int maximum_leaves = parse_positive(argv[4], "maximum leaves");
    if (maximum_background >= 31 || maximum_leaves < 2
        || maximum_background + maximum_leaves >= 63) {
        std::cerr << "invalid factor bound\n";
        return 2;
    }

    std::uint64_t backgrounds = 0U;
    std::uint64_t currents = 0U;
    std::uint64_t contraction_checks = 0U;
    std::uint64_t hypercube_checks = 0U;
    for (int level = 2; level <= maximum_level; ++level) {
        const int label_limit = std::min(level, maximum_label);
        const Fusion fusion(level);
        for (int background_size = 0;
             background_size <= maximum_background; ++background_size) {
            std::vector<int> background(
                static_cast<std::size_t>(background_size), 1
            );
            bool more_backgrounds = true;
            while (more_backgrounds) {
                const unsigned int masks = 1U << background.size();
                for (unsigned int mask = 0U; mask < masks; ++mask) {
                    int minus_count = 0;
                    for (std::size_t index = 0U;
                         index < background.size(); ++index) {
                        minus_count += static_cast<int>(
                            (mask & (1U << index)) != 0U
                        );
                    }
                    if ((minus_count & 1) != 0) {
                        continue;
                    }
                    const Table base = background_table(
                        fusion, background, mask
                    );
                    if (!symmetric(base)) {
                        std::cerr << "EVEN_BACKGROUND symmetry failure"
                                  << " level=" << level << '\n';
                        return 1;
                    }
                    ++backgrounds;
                    for (int leaf_count = 2;
                         leaf_count <= maximum_leaves; ++leaf_count) {
                        std::vector<int> leaves(
                            static_cast<std::size_t>(leaf_count), 1
                        );
                        bool more_leaves = true;
                        while (more_leaves) {
                            if (mask == 0U) {
                                std::vector<int> full_labels = background;
                                full_labels.insert(
                                    full_labels.end(), leaves.begin(), leaves.end()
                                );
                                if (!verify_two_bit_hypercube(
                                        fusion, full_labels, hypercube_checks
                                    )) {
                                    return 1;
                                }
                            }
                            for (int q = 1; q <= label_limit; ++q) {
                                std::vector<Integer> beta(
                                    leaves.size()
                                );
                                Integer lhs = 0;
                                for (std::size_t index = 0U;
                                     index < leaves.size(); ++index) {
                                    const Table table = plus_leaves(
                                        base, fusion, leaves,
                                        index, leaves.size()
                                    );
                                    beta[index] = boundary(
                                        table, fusion, q, leaves[index]
                                    );
                                    lhs += beta[index];
                                    Table direct = update(
                                        table, fusion, q, -1
                                    );
                                    direct = update(
                                        direct, fusion, leaves[index], -1
                                    );
                                    if (direct.at(0, 0) != 2 * beta[index]) {
                                        std::cerr
                                            << "EVEN_BACKGROUND contraction failure"
                                            << " level=" << level
                                            << " q=" << q
                                            << " leaf=" << leaves[index]
                                            << '\n';
                                        return 1;
                                    }
                                    ++contraction_checks;
                                }
                                lhs *= static_cast<int>(leaves.size() - 1U);
                                Integer rhs = 0;
                                for (std::size_t first = 0U;
                                     first < leaves.size(); ++first) {
                                    for (std::size_t second = first + 1U;
                                         second < leaves.size(); ++second) {
                                        const Table table = plus_leaves(
                                            base, fusion, leaves,
                                            first, second
                                        );
                                        for (const int fused
                                             : fusion.channels[static_cast<std::size_t>(leaves[first])]
                                                              [static_cast<std::size_t>(leaves[second])]) {
                                            rhs += 2 * boundary(
                                                table, fusion, q, fused
                                            );
                                        }
                                    }
                                }
                                if (lhs != rhs) {
                                    std::cerr << "EVEN_BACKGROUND current failure"
                                              << " level=" << level
                                              << " q=" << q << '\n';
                                    return 1;
                                }
                                ++currents;
                            }
                            more_leaves = next_word(leaves, label_limit);
                        }
                    }
                }
                if (background.empty()) {
                    more_backgrounds = false;
                } else {
                    more_backgrounds = next_word(background, label_limit);
                }
            }
        }
    }
    std::cout << "SU2_EVEN_BACKGROUND_CURRENT"
              << " levels=2.." << maximum_level
              << " labels<=" << maximum_label
              << " background_factors<=" << maximum_background
              << " leaves<=" << maximum_leaves
              << " backgrounds=" << backgrounds
              << " currents=" << currents
              << " contractions=" << contraction_checks
              << " hypercube_edges=" << hypercube_checks
              << " result=PASS\n";
    return 0;
}
