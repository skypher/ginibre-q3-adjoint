#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <z3++.h>

namespace {

using Exponents = std::array<int, 4U>;
using Orbit = std::pair<Exponents, Exponents>;

struct Column {
    int first = 0;
    int second = 0;
};

struct DegreeTemplate {
    int kind = 0;  // 0: constant, 1: q-offset, 2: (q-1)/2+offset
    int offset = 0;

    int value(int q) const {
        if (kind == 0) {
            return offset;
        }
        if (kind == 1) {
            return q - offset;
        }
        return (q - 1) / 2 + offset;
    }
};

struct ColumnTemplate {
    DegreeTemplate first;
    DegreeTemplate second;

    Column value(int q) const {
        return {first.value(q), second.value(q)};
    }
};

struct AssignmentResult {
    bool unique = false;
    std::int64_t weight = 0;
    std::vector<std::size_t> row_for_column;
};

struct SymbolicOrbit {
    z3::expr active;
    std::vector<z3::expr> first;
    std::vector<z3::expr> second;
};

struct SymbolicEntry {
    z3::expr available;
    z3::expr weight;
};

std::vector<Orbit> build_orbits(
    int q, const std::array<int, 3U>& even
) {
    std::vector<Orbit> orbits;
    for (std::size_t selected = 0U; selected < even.size(); ++selected) {
        Exponents first{};
        if (even[selected] == q - 1) {
            first[0U] = 1;
            first[selected + 1U] = q - 1;
        } else if (even[selected] == q + 1) {
            first[selected + 1U] = q;
        } else {
            continue;
        }
        std::array<std::size_t, 2U> remaining{};
        std::size_t position = 0U;
        for (std::size_t root = 0U; root < even.size(); ++root) {
            if (root != selected) {
                remaining[position++] = root;
            }
        }
        const int x = even[remaining[0U]];
        const int y = even[remaining[1U]];
        for (int channel = 0; channel < 2; ++channel) {
            Exponents second{};
            bool active = false;
            if (channel == 0) {
                active = y - x <= q - 1 && x + y >= q - 1;
                if (active) {
                    second[0U] = 1;
                    second[remaining[0U] + 1U] = (q - 1 + x - y) / 2;
                    second[remaining[1U] + 1U] = (q - 1 - x + y) / 2;
                }
            } else if (y - x == q + 1) {
                active = true;
                second[remaining[1U] + 1U] = q;
            } else {
                active = y - x <= q - 1 && x + y >= q + 1;
                if (active) {
                    second[remaining[0U] + 1U] = (q + x - y - 1) / 2;
                    second[remaining[1U] + 1U] = (q - x + y + 1) / 2;
                }
            }
            if (active) {
                orbits.push_back({first, second});
            }
        }
    }
    return orbits;
}

std::optional<int> form_weight(const Exponents& exponents, int degree) {
    // Specialize (V,A,B,C)=(X,Y,X+Y,X+tY).  The coefficient of
    // X^(q-degree)Y^degree has highest t-power min(c,degree-a).
    const int a = exponents[1U];
    const int b = exponents[2U];
    const int c = exponents[3U];
    if (degree < a || degree > a + b + c) {
        return std::nullopt;
    }
    return std::min(c, degree - a);
}

std::optional<int> tensor_weight(
    const Orbit& orbit, const Column& column, bool antisymmetric
) {
    const std::optional<int> ff = form_weight(orbit.first, column.first);
    const std::optional<int> fs = form_weight(orbit.first, column.second);
    const std::optional<int> hf = form_weight(orbit.second, column.first);
    const std::optional<int> hs = form_weight(orbit.second, column.second);
    std::optional<int> forward;
    std::optional<int> reverse;
    if (ff.has_value() && hs.has_value()) {
        forward = *ff + *hs;
    }
    if (hf.has_value() && fs.has_value()) {
        reverse = *hf + *fs;
    }
    if (!forward.has_value()) {
        return reverse;
    }
    if (!reverse.has_value()) {
        return forward;
    }
    if (antisymmetric && *forward == *reverse) {
        // The two leading coefficients could cancel.  Reject this entry;
        // the search only emits certificates immune to such cancellation.
        return std::nullopt;
    }
    return std::max(*forward, *reverse);
}

AssignmentResult unique_max_assignment(
    const std::vector<Orbit>& rows,
    const std::vector<Column>& columns,
    bool antisymmetric
) {
    const std::size_t size = rows.size();
    if (columns.size() != size || size >= 20U) {
        throw std::runtime_error("invalid tropical assignment dimensions");
    }
    constexpr std::int64_t negative_infinity
        = std::numeric_limits<std::int64_t>::min() / 4;
    const std::size_t states = std::size_t{1} << size;
    std::vector<std::int64_t> best(states, negative_infinity);
    std::vector<unsigned int> count(states, 0U);
    std::vector<std::vector<std::size_t>> parents(
        size + 1U, std::vector<std::size_t>(states, size)
    );
    best[0U] = 0;
    count[0U] = 1U;
    for (std::size_t column = 0U; column < size; ++column) {
        std::vector<std::int64_t> next(states, negative_infinity);
        std::vector<unsigned int> next_count(states, 0U);
        for (std::size_t mask = 0U; mask < states; ++mask) {
            if (count[mask] == 0U
                || static_cast<std::size_t>(__builtin_popcountll(mask))
                    != column) {
                continue;
            }
            for (std::size_t row = 0U; row < size; ++row) {
                if ((mask & (std::size_t{1} << row)) != 0U) {
                    continue;
                }
                const std::optional<int> weight
                    = tensor_weight(rows[row], columns[column], antisymmetric);
                if (!weight.has_value()) {
                    continue;
                }
                const std::size_t next_mask
                    = mask | (std::size_t{1} << row);
                const std::int64_t candidate = best[mask] + *weight;
                if (candidate > next[next_mask]) {
                    next[next_mask] = candidate;
                    next_count[next_mask] = count[mask];
                    parents[column + 1U][next_mask] = row;
                } else if (candidate == next[next_mask]) {
                    next_count[next_mask] = std::min(
                        2U, next_count[next_mask] + count[mask]
                    );
                    parents[column + 1U][next_mask] = size;
                }
            }
        }
        best.swap(next);
        count.swap(next_count);
    }
    const std::size_t full = states - 1U;
    AssignmentResult result;
    result.unique = count[full] == 1U;
    result.weight = best[full];
    if (!result.unique) {
        return result;
    }
    result.row_for_column.resize(size, size);
    std::size_t mask = full;
    for (std::size_t column = size; column != 0U; --column) {
        const std::size_t row = parents[column][mask];
        if (row == size) {
            throw std::runtime_error("lost unique assignment parent");
        }
        result.row_for_column[column - 1U] = row;
        mask ^= std::size_t{1} << row;
    }
    return result;
}

std::optional<std::pair<std::vector<Column>, AssignmentResult>> search_minor(
    const std::vector<Orbit>& rows, int q, bool antisymmetric
) {
    std::set<int> degree_set;
    for (int offset = 0; offset <= 4; ++offset) {
        degree_set.insert(offset);
        degree_set.insert(q - offset);
    }
    const int half = q / 2;
    for (int offset = -2; offset <= 3; ++offset) {
        degree_set.insert(half + offset);
    }
    std::vector<int> degrees;
    for (int degree : degree_set) {
        if (degree >= 0 && degree <= q) {
            degrees.push_back(degree);
        }
    }
    std::vector<Column> candidates;
    for (int first : degrees) {
        for (int second : degrees) {
            if (second < first || (antisymmetric && second == first)) {
                continue;
            }
            bool usable = false;
            for (const Orbit& row : rows) {
                usable = usable
                    || tensor_weight(row, {first, second}, antisymmetric)
                        .has_value();
            }
            if (usable) {
                candidates.push_back({first, second});
            }
        }
    }
    if (candidates.size() < rows.size()) {
        return std::nullopt;
    }
    using Seed = std::mt19937_64::result_type;
    const Seed seed
        = static_cast<Seed>(UINT64_C(0x9e3779b97f4a7c15))
            ^ static_cast<Seed>(q)
            ^ (antisymmetric
                ? static_cast<Seed>(UINT64_C(0x51ed270b)) : Seed{0});
    std::mt19937_64 generator(seed);
    std::vector<std::size_t> indices(candidates.size());
    for (std::size_t index = 0U; index < indices.size(); ++index) {
        indices[index] = index;
    }
    constexpr std::size_t attempts = 250'000U;
    for (std::size_t attempt = 0U; attempt < attempts; ++attempt) {
        for (std::size_t chosen = 0U; chosen < rows.size(); ++chosen) {
            std::uniform_int_distribution<std::size_t> distribution(
                chosen, indices.size() - 1U
            );
            const std::size_t replacement = distribution(generator);
            std::swap(indices[chosen], indices[replacement]);
        }
        std::vector<Column> columns;
        columns.reserve(rows.size());
        for (std::size_t chosen = 0U; chosen < rows.size(); ++chosen) {
            columns.push_back(candidates[indices[chosen]]);
        }
        const AssignmentResult result = unique_max_assignment(
            rows, columns, antisymmetric
        );
        if (result.unique) {
            return std::make_pair(std::move(columns), result);
        }
    }
    return std::nullopt;
}

std::array<int, 3U> family_even(int q, int family) {
    const int qm = q - 1;
    const int qp = q + 1;
    const std::array<std::array<int, 3U>, 13U> families{
        std::array<int, 3U>{2, qm, qp},
        std::array<int, 3U>{2, qp, qp},
        std::array<int, 3U>{qm, qm, qm},
        std::array<int, 3U>{qm, qm, qp},
        std::array<int, 3U>{qm, qp, qp},
        std::array<int, 3U>{qp, qp, qp},
        std::array<int, 3U>{qm, qm, 2 * q - 4},
        std::array<int, 3U>{qm, qm, 2 * q - 2},
        std::array<int, 3U>{qm, qp, 2 * q - 2},
        std::array<int, 3U>{qm, qp, 2 * q},
        std::array<int, 3U>{qp, qp, 2 * q - 2},
        std::array<int, 3U>{qp, qp, 2 * q},
        std::array<int, 3U>{qp, qp, 2 * q + 2}
    };
    if (family < 0 || family >= static_cast<int>(families.size())) {
        throw std::runtime_error("invalid residual-family index");
    }
    return families[static_cast<std::size_t>(family)];
}

std::string degree_template_name(const DegreeTemplate& degree) {
    if (degree.kind == 0) {
        return std::to_string(degree.offset);
    }
    if (degree.kind == 1) {
        return degree.offset == 0 ? "q"
            : "q-" + std::to_string(degree.offset);
    }
    if (degree.offset == 0) {
        return "h";
    }
    return degree.offset > 0
        ? "h+" + std::to_string(degree.offset)
        : "h" + std::to_string(degree.offset);
}

std::vector<std::size_t> active_source_indices(
    int q, const std::array<int, 3U>& even
) {
    std::vector<std::size_t> indices;
    for (std::size_t selected = 0U; selected < even.size(); ++selected) {
        const bool selected_active
            = even[selected] == q - 1 || even[selected] == q + 1;
        std::array<std::size_t, 2U> remaining{};
        std::size_t position = 0U;
        for (std::size_t root = 0U; root < even.size(); ++root) {
            if (root != selected) {
                remaining[position++] = root;
            }
        }
        const int x = even[remaining[0U]];
        const int y = even[remaining[1U]];
        const bool lower_active
            = y - x <= q - 1 && x + y >= q - 1;
        const bool upper_active = y - x == q + 1
            || (y - x <= q - 1 && x + y >= q + 1);
        if (selected_active && lower_active) {
            indices.push_back(2U * selected);
        }
        if (selected_active && upper_active) {
            indices.push_back(2U * selected + 1U);
        }
    }
    return indices;
}

std::vector<z3::expr> symbolic_even_labels(
    z3::context& context, const z3::expr& q, int family
) {
    const z3::expr two = context.int_val(2);
    switch (family) {
        case 0: return {two, q - 1, q + 1};
        case 1: return {two, q + 1, q + 1};
        case 2: return {q - 1, q - 1, q - 1};
        case 3: return {q - 1, q - 1, q + 1};
        case 4: return {q - 1, q + 1, q + 1};
        case 5: return {q + 1, q + 1, q + 1};
        case 6: return {q - 1, q - 1, 2 * q - 4};
        case 7: return {q - 1, q - 1, 2 * q - 2};
        case 8: return {q - 1, q + 1, 2 * q - 2};
        case 9: return {q - 1, q + 1, 2 * q};
        case 10: return {q + 1, q + 1, 2 * q - 2};
        case 11: return {q + 1, q + 1, 2 * q};
        case 12: return {q + 1, q + 1, 2 * q + 2};
        default: throw std::runtime_error("invalid symbolic family");
    }
}

std::vector<SymbolicOrbit> build_symbolic_potential_orbits(
    z3::context& context, const z3::expr& q,
    const std::vector<z3::expr>& even
) {
    std::vector<SymbolicOrbit> sources;
    sources.reserve(6U);
    for (std::size_t selected = 0U; selected < 3U; ++selected) {
        const z3::expr lower = even[selected] == q - 1;
        const z3::expr upper = even[selected] == q + 1;
        std::vector<z3::expr> first(4U, context.int_val(0));
        first[0U] = z3::ite(lower, context.int_val(1), context.int_val(0));
        first[selected + 1U] = z3::ite(lower, q - 1, q);
        std::array<std::size_t, 2U> remaining{};
        std::size_t position = 0U;
        for (std::size_t root = 0U; root < 3U; ++root) {
            if (root != selected) {
                remaining[position++] = root;
            }
        }
        const z3::expr x = even[remaining[0U]];
        const z3::expr y = even[remaining[1U]];
        for (int channel = 0; channel < 2; ++channel) {
            std::vector<z3::expr> second(4U, context.int_val(0));
            z3::expr complement_active = context.bool_val(false);
            if (channel == 0) {
                complement_active = y - x <= q - 1 && x + y >= q - 1;
                second[0U] = context.int_val(1);
                second[remaining[0U] + 1U] = (q - 1 + x - y) / 2;
                second[remaining[1U] + 1U] = (q - 1 - x + y) / 2;
            } else {
                const z3::expr wall = y - x == q + 1;
                const z3::expr interior
                    = y - x <= q - 1 && x + y >= q + 1;
                complement_active = wall || interior;
                second[remaining[0U] + 1U] = z3::ite(
                    wall, context.int_val(0), (q + x - y - 1) / 2
                );
                second[remaining[1U] + 1U] = z3::ite(
                    wall, q, (q - x + y + 1) / 2
                );
            }
            sources.push_back({
                (lower || upper) && complement_active,
                first, std::move(second)
            });
        }
    }
    return sources;
}

z3::expr symbolic_degree(
    z3::context& context, const z3::expr& q,
    const z3::expr& half, const DegreeTemplate& degree
) {
    if (degree.kind == 0) {
        return context.int_val(degree.offset);
    }
    if (degree.kind == 1) {
        return q - degree.offset;
    }
    return half + degree.offset;
}

SymbolicEntry symbolic_form_entry(
    const std::vector<z3::expr>& exponents, const z3::expr& degree
) {
    const z3::expr residual = degree - exponents[1U];
    return {
        degree >= exponents[1U]
            && degree <= exponents[1U] + exponents[2U] + exponents[3U],
        z3::ite(exponents[3U] <= residual, exponents[3U], residual)
    };
}

SymbolicEntry symbolic_tensor_entry(
    const SymbolicOrbit& orbit, const z3::expr& first_degree,
    const z3::expr& second_degree, bool antisymmetric
) {
    const SymbolicEntry ff = symbolic_form_entry(orbit.first, first_degree);
    const SymbolicEntry fs = symbolic_form_entry(orbit.first, second_degree);
    const SymbolicEntry hf = symbolic_form_entry(orbit.second, first_degree);
    const SymbolicEntry hs = symbolic_form_entry(orbit.second, second_degree);
    const z3::expr forward_available = ff.available && hs.available;
    const z3::expr reverse_available = hf.available && fs.available;
    const z3::expr forward_weight = ff.weight + hs.weight;
    const z3::expr reverse_weight = hf.weight + fs.weight;
    const z3::expr choose_forward = forward_available
        && (!reverse_available || forward_weight >= reverse_weight);
    const z3::expr available = antisymmetric
        ? ((forward_available && !reverse_available)
            || (!forward_available && reverse_available)
            || (forward_available && reverse_available
                && forward_weight != reverse_weight))
        : (forward_available || reverse_available);
    return {
        available,
        z3::ite(choose_forward, forward_weight, reverse_weight)
    };
}

bool verify_uniform_template_exact(
    int family, bool antisymmetric,
    const std::array<std::size_t, 4U>& permutation,
    const std::vector<ColumnTemplate>& columns,
    const std::vector<std::size_t>& winning_rows
) {
    z3::context context;
    const z3::expr half = context.int_const("h");
    const z3::expr q = 2 * half + 1;
    const std::vector<z3::expr> even
        = symbolic_even_labels(context, q, family);
    const std::vector<SymbolicOrbit> potential
        = build_symbolic_potential_orbits(context, q, even);
    const std::vector<std::size_t> active
        = active_source_indices(7, family_even(7, family));
    std::vector<std::size_t> selected = active;
    if (antisymmetric) {
        const std::vector<Orbit> concrete
            = build_orbits(7, family_even(7, family));
        selected.clear();
        for (std::size_t position = 0U; position < concrete.size();
             ++position) {
            if (concrete[position].first != concrete[position].second) {
                selected.push_back(active[position]);
            }
        }
    }
    std::vector<SymbolicOrbit> rows;
    z3::expr pattern_mismatch = context.bool_val(false);
    for (std::size_t source = 0U; source < potential.size(); ++source) {
        const bool expected = std::find(active.begin(), active.end(), source)
            != active.end();
        pattern_mismatch = pattern_mismatch
            || (expected ? !potential[source].active
                         : potential[source].active);
        const bool included
            = std::find(selected.begin(), selected.end(), source)
                != selected.end();
        if (!included) {
            continue;
        }
        std::vector<z3::expr> first;
        std::vector<z3::expr> second;
        first.reserve(4U);
        second.reserve(4U);
        for (std::size_t coordinate = 0U; coordinate < 4U; ++coordinate) {
            first.push_back(
                potential[source].first[permutation[coordinate]]
            );
            second.push_back(
                potential[source].second[permutation[coordinate]]
            );
        }
        if (antisymmetric) {
            z3::expr unequal = context.bool_val(false);
            for (std::size_t coordinate = 0U; coordinate < 4U;
                 ++coordinate) {
                unequal = unequal || first[coordinate] != second[coordinate];
            }
            // Every retained antisymmetric source in these twelve families
            // is nonfixed.  Include failure of that fact in the query.
            pattern_mismatch = pattern_mismatch || !unequal;
        }
        rows.push_back({context.bool_val(true), first, second});
    }
    if (rows.size() != columns.size()
        || winning_rows.size() != columns.size()) {
        throw std::runtime_error("symbolic template dimension mismatch");
    }
    const std::size_t size = rows.size();
    std::vector<std::vector<SymbolicEntry>> matrix;
    matrix.reserve(size);
    for (std::size_t column = 0U; column < size; ++column) {
        const z3::expr first_degree = symbolic_degree(
            context, q, half, columns[column].first
        );
        const z3::expr second_degree = symbolic_degree(
            context, q, half, columns[column].second
        );
        std::vector<SymbolicEntry> entries;
        entries.reserve(size);
        for (const SymbolicOrbit& row : rows) {
            entries.push_back(symbolic_tensor_entry(
                row, first_degree, second_degree, antisymmetric
            ));
        }
        matrix.push_back(std::move(entries));
    }
    z3::expr winner_available = context.bool_val(true);
    z3::expr winner_weight = context.int_val(0);
    for (std::size_t column = 0U; column < size; ++column) {
        const std::size_t row = winning_rows[column];
        winner_available = winner_available && matrix[column][row].available;
        winner_weight = winner_weight + matrix[column][row].weight;
    }
    z3::expr competitor_wins = context.bool_val(false);
    std::vector<std::size_t> assignment(size);
    for (std::size_t row = 0U; row < size; ++row) {
        assignment[row] = row;
    }
    do {
        if (assignment == winning_rows) {
            continue;
        }
        z3::expr available = context.bool_val(true);
        z3::expr weight = context.int_val(0);
        for (std::size_t column = 0U; column < size; ++column) {
            available = available
                && matrix[column][assignment[column]].available;
            weight = weight + matrix[column][assignment[column]].weight;
        }
        competitor_wins = competitor_wins
            || (available && weight >= winner_weight);
    } while (std::next_permutation(assignment.begin(), assignment.end()));
    z3::tactic qflia(context, "qflia");
    z3::solver solver = qflia.mk_solver();
    solver.add(half >= 3);
    solver.add(pattern_mismatch || !winner_available || competitor_wins);
    return solver.check() == z3::unsat;
}

std::optional<std::pair<std::vector<ColumnTemplate>, AssignmentResult>>
search_uniform_family(
    int family, bool antisymmetric,
    const std::array<std::size_t, 4U>& permutation
) {
    const std::array<int, 8U> samples{7, 9, 11, 15, 21, 31, 41, 61};
    std::vector<std::vector<Orbit>> rows_by_sample;
    for (int q : samples) {
        std::vector<Orbit> rows = build_orbits(q, family_even(q, family));
        for (Orbit& row : rows) {
            Exponents first{};
            Exponents second{};
            for (std::size_t coordinate = 0U; coordinate < 4U;
                 ++coordinate) {
                first[coordinate] = row.first[permutation[coordinate]];
                second[coordinate] = row.second[permutation[coordinate]];
            }
            row = {first, second};
        }
        if (antisymmetric) {
            rows.erase(
                std::remove_if(
                    rows.begin(), rows.end(),
                    [](const Orbit& row) { return row.first == row.second; }
                ),
                rows.end()
            );
        }
        rows_by_sample.push_back(std::move(rows));
    }
    const std::size_t size = rows_by_sample.front().size();
    for (const std::vector<Orbit>& rows : rows_by_sample) {
        if (rows.size() != size) {
            throw std::runtime_error("family row count is not stable");
        }
    }
    std::vector<DegreeTemplate> degrees;
    for (int value = 0; value <= 6; ++value) {
        degrees.push_back({0, value});
        degrees.push_back({1, value});
    }
    for (int offset = -3; offset <= 4; ++offset) {
        degrees.push_back({2, offset});
    }
    std::vector<ColumnTemplate> candidates;
    for (const DegreeTemplate& first : degrees) {
        for (const DegreeTemplate& second : degrees) {
            bool valid = true;
            for (int q : samples) {
                const int first_value = first.value(q);
                const int second_value = second.value(q);
                valid = valid && first_value >= 0 && second_value <= q
                    && first_value <= second_value
                    && (!antisymmetric || first_value < second_value);
            }
            if (valid) {
                candidates.push_back({first, second});
            }
        }
    }
    std::vector<std::size_t> indices(candidates.size());
    for (std::size_t index = 0U; index < indices.size(); ++index) {
        indices[index] = index;
    }
    using Seed = std::mt19937_64::result_type;
    const Seed seed
        = static_cast<Seed>(UINT64_C(0xd1b54a32d192ed03))
            ^ static_cast<Seed>(family)
            ^ (static_cast<Seed>(permutation[0U]) << 8U)
            ^ (static_cast<Seed>(permutation[1U]) << 16U)
            ^ (static_cast<Seed>(permutation[2U]) << 24U)
            ^ (antisymmetric
                ? static_cast<Seed>(UINT64_C(0x94d049bb133111eb))
                : Seed{0});
    std::mt19937_64 generator(seed);
    constexpr std::size_t attempts = 500'000U;
    for (std::size_t attempt = 0U; attempt < attempts; ++attempt) {
        for (std::size_t chosen = 0U; chosen < size; ++chosen) {
            std::uniform_int_distribution<std::size_t> distribution(
                chosen, indices.size() - 1U
            );
            const std::size_t replacement = distribution(generator);
            std::swap(indices[chosen], indices[replacement]);
        }
        std::vector<ColumnTemplate> selected;
        selected.reserve(size);
        for (std::size_t chosen = 0U; chosen < size; ++chosen) {
            selected.push_back(candidates[indices[chosen]]);
        }
        std::optional<AssignmentResult> common;
        bool passes = true;
        for (std::size_t sample = 0U; sample < samples.size(); ++sample) {
            std::vector<Column> columns;
            columns.reserve(size);
            for (const ColumnTemplate& column : selected) {
                columns.push_back(column.value(samples[sample]));
            }
            std::sort(
                columns.begin(), columns.end(),
                [](const Column& first, const Column& second) {
                    return std::tie(first.first, first.second)
                        < std::tie(second.first, second.second);
                }
            );
            if (std::adjacent_find(
                    columns.begin(), columns.end(),
                    [](const Column& first, const Column& second) {
                        return first.first == second.first
                            && first.second == second.second;
                    }
                ) != columns.end()) {
                passes = false;
                break;
            }
            // Preserve template order for the common matching after the
            // duplicate check.
            columns.clear();
            for (const ColumnTemplate& column : selected) {
                columns.push_back(column.value(samples[sample]));
            }
            const AssignmentResult result = unique_max_assignment(
                rows_by_sample[sample], columns, antisymmetric
            );
            if (!result.unique) {
                passes = false;
                break;
            }
            if (!common.has_value()) {
                common = result;
            } else if (common->row_for_column != result.row_for_column) {
                passes = false;
                break;
            }
        }
        if (passes && common.has_value()
            && verify_uniform_template_exact(
                family, antisymmetric, permutation, selected,
                common->row_for_column
            )) {
            return std::make_pair(std::move(selected), *common);
        }
    }
    return std::nullopt;
}

void print_uniform_certificate(int family, bool antisymmetric) {
    std::cout << "SU2_EQUAL_TROPICAL_UNIFORM family=" << family
              << " sector=" << (antisymmetric ? "minus" : "plus");
    std::array<std::size_t, 4U> permutation{0U, 1U, 2U, 3U};
    std::optional<std::pair<std::vector<ColumnTemplate>, AssignmentResult>>
        certificate;
    std::array<std::size_t, 4U> winning_permutation = permutation;
    do {
        certificate = search_uniform_family(
            family, antisymmetric, permutation
        );
        if (certificate.has_value()) {
            winning_permutation = permutation;
            break;
        }
    } while (std::next_permutation(
        permutation.begin(), permutation.end()
    ));
    if (!certificate.has_value()) {
        std::cout << " result=not-found\n";
        return;
    }
    std::cout << " result=exact-uniform permutation=";
    for (std::size_t coordinate : winning_permutation) {
        std::cout << coordinate;
    }
    std::cout << " columns=";
    for (std::size_t index = 0U; index < certificate->first.size(); ++index) {
        if (index != 0U) {
            std::cout << ';';
        }
        std::cout << degree_template_name(certificate->first[index].first)
                  << ','
                  << degree_template_name(certificate->first[index].second)
                  << "->" << certificate->second.row_for_column[index];
    }
    std::cout << '\n';
}

void print_certificate(
    int q, const std::array<int, 3U>& even,
    const std::vector<Orbit>& rows, bool antisymmetric
) {
    const auto certificate = search_minor(rows, q, antisymmetric);
    std::cout << "SU2_EQUAL_TROPICAL q=" << q << " even="
              << even[0U] << ',' << even[1U] << ',' << even[2U]
              << " sector=" << (antisymmetric ? "minus" : "plus")
              << " rows=" << rows.size();
    if (!certificate.has_value()) {
        std::cout << " result=not-found\n";
        return;
    }
    std::cout << " result=unique weight=" << certificate->second.weight
              << " columns=";
    for (std::size_t index = 0U; index < certificate->first.size(); ++index) {
        if (index != 0U) {
            std::cout << ';';
        }
        std::cout << certificate->first[index].first << ','
                  << certificate->first[index].second << "->"
                  << certificate->second.row_for_column[index];
    }
    std::cout << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--uniform-family") {
            const int family = std::stoi(argv[2]);
            print_uniform_certificate(family, false);
            print_uniform_certificate(family, true);
            return EXIT_SUCCESS;
        }
        if (argc != 5) {
            std::cerr << "usage: search_su2_equal_terminal_tropical_minor"
                      << " Q E1 E2 E3\n"
                      << "   or: search_su2_equal_terminal_tropical_minor"
                      << " --uniform-family INDEX\n";
            return EXIT_FAILURE;
        }
        const int q = std::stoi(argv[1]);
        const std::array<int, 3U> even{
            std::stoi(argv[2]), std::stoi(argv[3]), std::stoi(argv[4])
        };
        if (q < 3 || (q & 1) == 0 || even[0U] < 2
            || even[0U] > even[1U] || even[1U] > even[2U]) {
            throw std::runtime_error("invalid labels");
        }
        const std::vector<Orbit> rows = build_orbits(q, even);
        std::vector<Orbit> minus_rows;
        for (const Orbit& row : rows) {
            if (row.first != row.second) {
                minus_rows.push_back(row);
            }
        }
        print_certificate(q, even, rows, false);
        print_certificate(q, even, minus_rows, true);
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
