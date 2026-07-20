#include "ankerl/unordered_dense.h"
#include "lie_data.h"

#include <gmpxx.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

template <int Rank>
struct Weight {
    std::array<std::int32_t, Rank> coordinate{};

    bool operator==(const Weight& other) const {
        return coordinate == other.coordinate;
    }
};

template <int Rank>
struct WeightHash {
    using is_avalanching = void;

    std::size_t operator()(const Weight<Rank>& weight) const noexcept {
        std::size_t hash = 1469598103934665603ULL;
        for (std::int32_t coordinate : weight.coordinate) {
            const std::uint32_t value = static_cast<std::uint32_t>(coordinate);
            for (int byte = 0; byte < 4; ++byte) {
                hash ^= static_cast<std::size_t>((value >> (8 * byte)) & 0xffU);
                hash *= 1099511628211ULL;
            }
        }
        return hash;
    }
};

template <int Rank>
using MultiplicityMap =
    ankerl::unordered_dense::map<Weight<Rank>, mpz_class, WeightHash<Rank>>;

template <int Rank>
using WeightSet = ankerl::unordered_dense::set<Weight<Rank>, WeightHash<Rank>>;

template <int Rank>
struct RootDatum {
    std::string name;
    std::array<std::array<int, Rank>, Rank> cartan{};
    std::vector<Weight<Rank>> adjoint_weights;
    std::vector<int> adjoint_multiplicities;
    int expected_root_count = 0;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

std::map<int, mpz_class> read_expected_moments(
    const std::vector<std::string>& paths) {
    const std::regex moment_pattern(R"(\bm_\s*(\d+)\s*=\s*(-?\d+))");
    std::map<int, mpz_class> moments;
    for (const std::string& path : paths) {
        std::ifstream input(path);
        if (!input) fail("cannot open expected-moment source: " + path);
        std::string line;
        while (std::getline(input, line)) {
            for (std::sregex_iterator it(line.begin(), line.end(), moment_pattern), end;
                 it != end; ++it) {
                const int degree = std::stoi((*it)[1].str());
                const mpz_class value((*it)[2].str());
                const auto [found, inserted] = moments.emplace(degree, value);
                if (!inserted && found->second != value) {
                    fail("conflicting source values for m_" + std::to_string(degree));
                }
            }
        }
    }
    return moments;
}

struct ClassicalClaims {
    std::map<int, mpz_class> moments;
    std::map<int, mpz_class> deltas;
    std::size_t ignored_determinant_components = 0;

    std::size_t size() const { return moments.size() + deltas.size(); }
};

void insert_consistent(
    std::map<int, mpz_class>& values,
    int degree,
    const mpz_class& value,
    const std::string& description) {
    const auto [found, inserted] = values.emplace(degree, value);
    if (!inserted && found->second != value) {
        fail("conflicting " + description + " values in degree " +
             std::to_string(degree));
    }
}

ClassicalClaims read_classical_claims(
    const std::vector<std::string>& paths,
    char family,
    int rank) {
    const std::regex delta_pattern(R"(([BCD])_(\d+).*?Delta_(\d+)\s*=\s*(-?\d+))");
    const std::regex moment_pattern(
        R"(D_(\d+)\s+m_(\d+).*?moment_\d+\s*=\s*(-?\d+))");
    const std::regex direct_moment_pattern(
        R"(D_(\d+).*?m_(\d+)\s*=\s*(-?\d+))");
    const std::regex legacy_determinant_pattern(
        R"(^\s*D_(\d+)\s+Delta_(\d+)\s*=\s*(-?\d+)\s*$)");
    ClassicalClaims claims;
    for (const std::string& path : paths) {
        std::ifstream input(path);
        if (!input) fail("cannot open classical source: " + path);
        std::string line;
        while (std::getline(input, line)) {
            std::smatch legacy_match;
            if (family == 'D'
                && std::regex_match(line, legacy_match,
                                    legacy_determinant_pattern)
                && std::stoi(legacy_match[1].str()) == rank
                && std::stoi(legacy_match[2].str()) > 2 * rank) {
                // A bare out-of-window D-row "Delta" in the historical
                // lower-bound supplier is the nonnegative determinant
                // component B_j, not the full correction B_j+A_j.  The
                // downstream interval checker makes the same distinction
                // and ignores it as an exact moment.  Full/exact-Weyl lines
                // remain checked through their moment fields.
                ++claims.ignored_determinant_components;
                continue;
            }
            for (std::sregex_iterator it(line.begin(), line.end(), delta_pattern), end;
                 it != end; ++it) {
                if ((*it)[1].str()[0] != family || std::stoi((*it)[2].str()) != rank) {
                    continue;
                }
                const int degree = std::stoi((*it)[3].str());
                insert_consistent(claims.deltas, degree,
                                  mpz_class((*it)[4].str()), "Delta");
            }
            if (family != 'D') continue;
            for (std::sregex_iterator it(line.begin(), line.end(), moment_pattern), end;
                 it != end; ++it) {
                if (std::stoi((*it)[1].str()) != rank) continue;
                insert_consistent(claims.moments, std::stoi((*it)[2].str()),
                                  mpz_class((*it)[3].str()), "moment");
            }
            for (std::sregex_iterator it(line.begin(), line.end(), direct_moment_pattern), end;
                 it != end; ++it) {
                if (std::stoi((*it)[1].str()) != rank) continue;
                insert_consistent(claims.moments, std::stoi((*it)[2].str()),
                                  mpz_class((*it)[3].str()), "moment");
            }
        }
    }
    if (claims.size() == 0) {
        fail("no classical source claims found for " + std::string(1, family) + "_" +
             std::to_string(rank));
    }
    return claims;
}

std::map<int, mpz_class> read_stable_moments(const std::string& path) {
    std::ifstream input(path);
    if (!input) fail("cannot open stable-moment comparison table: " + path);
    std::map<int, mpz_class> comparison;
    int degree = 0;
    std::string value;
    while (input >> degree >> value) {
        if (degree < 0) fail("stable-moment table contains a negative degree");
        const auto [_, inserted] = comparison.emplace(degree, mpz_class(value));
        if (!inserted) {
            fail("stable-moment table repeats degree " + std::to_string(degree));
        }
    }
    if (!input.eof()) fail("malformed stable-moment comparison table: " + path);
    if (comparison.empty()) fail("stable-moment comparison table is empty");

    const int maximum_degree = comparison.rbegin()->first;
    for (int j = 0; j <= maximum_degree; ++j) {
        if (!comparison.contains(j)) {
            fail("stable-moment comparison table omits degree " +
                 std::to_string(j));
        }
    }

    // Lemma bc-stable-moment-recurrence in the paper proves
    //   s_{n+1}=n s_n+n s_{n-1}-binom(n,2)s_{n-2},
    // with s_0=1 and s_1=0.  These regenerated values, rather than the
    // external table, are the values used in every correction comparison.
    std::map<int, mpz_class> regenerated;
    regenerated.emplace(0, mpz_class(1));
    if (maximum_degree >= 1) regenerated.emplace(1, mpz_class(0));
    for (int n = 1; n < maximum_degree; ++n) {
        mpz_class next = mpz_class(n) * regenerated.at(n)
                       + mpz_class(n) * regenerated.at(n - 1);
        if (n >= 2) {
            const mpz_class binomial = mpz_class(n) * mpz_class(n - 1) / 2;
            next -= binomial * regenerated.at(n - 2);
        }
        regenerated.emplace(n + 1, std::move(next));
    }
    for (int j = 0; j <= maximum_degree; ++j) {
        if (regenerated.at(j) != comparison.at(j)) {
            fail("stable-moment recurrence/table mismatch in degree " +
                 std::to_string(j));
        }
    }
    return regenerated;
}

template <int Rank>
Weight<Rank> reflect(
    const Weight<Rank>& weight,
    int simple_root,
    const std::array<std::array<int, Rank>, Rank>& cartan) {
    Weight<Rank> result = weight;
    const std::int64_t coefficient = weight.coordinate[simple_root];
    for (int column = 0; column < Rank; ++column) {
        const std::int64_t value =
            static_cast<std::int64_t>(weight.coordinate[column]) -
            coefficient * cartan[simple_root][column];
        if (value < std::numeric_limits<std::int32_t>::min() ||
            value > std::numeric_limits<std::int32_t>::max()) {
            fail("weight-coordinate overflow in simple reflection");
        }
        result.coordinate[column] = static_cast<std::int32_t>(value);
    }
    return result;
}

template <int Rank>
WeightSet<Rank> root_closure(
    const std::array<std::array<int, Rank>, Rank>& cartan) {
    WeightSet<Rank> roots;
    std::vector<Weight<Rank>> frontier;
    for (int row = 0; row < Rank; ++row) {
        Weight<Rank> root;
        for (int column = 0; column < Rank; ++column) {
            root.coordinate[column] = cartan[row][column];
        }
        if (roots.insert(root).second) frontier.push_back(root);
    }
    for (std::size_t index = 0; index < frontier.size(); ++index) {
        for (int simple_root = 0; simple_root < Rank; ++simple_root) {
            Weight<Rank> image = reflect<Rank>(frontier[index], simple_root, cartan);
            if (roots.insert(image).second) frontier.push_back(image);
        }
    }
    return roots;
}

template <int Rank>
void validate_root_datum(const RootDatum<Rank>& datum) {
    for (int row = 0; row < Rank; ++row) {
        if (datum.cartan[row][row] != 2) fail("Cartan diagonal is not two");
        for (int column = 0; column < Rank; ++column) {
            if (row != column && datum.cartan[row][column] > 0) {
                fail("Cartan off-diagonal entry is positive");
            }
            if ((datum.cartan[row][column] == 0) !=
                (datum.cartan[column][row] == 0)) {
                fail("Cartan zero pattern is not symmetric");
            }
        }
    }

    if (datum.adjoint_weights.size() != datum.adjoint_multiplicities.size()) {
        fail("adjoint weight/multiplicity arrays have different sizes");
    }
    WeightSet<Rank> listed_roots;
    int zero_entries = 0;
    for (std::size_t index = 0; index < datum.adjoint_weights.size(); ++index) {
        const Weight<Rank>& weight = datum.adjoint_weights[index];
        const bool zero = std::all_of(
            weight.coordinate.begin(), weight.coordinate.end(),
            [](std::int32_t value) { return value == 0; });
        if (zero) {
            ++zero_entries;
            if (datum.adjoint_multiplicities[index] != Rank) {
                fail("zero adjoint weight does not have rank multiplicity");
            }
        } else {
            if (datum.adjoint_multiplicities[index] != 1) {
                fail("nonzero adjoint root does not have multiplicity one");
            }
            if (!listed_roots.insert(weight).second) fail("duplicate adjoint root");
        }
    }
    if (zero_entries != 1) fail("adjoint datum does not have one zero-weight entry");

    const WeightSet<Rank> generated_roots = root_closure<Rank>(datum.cartan);
    if (static_cast<int>(generated_roots.size()) != datum.expected_root_count) {
        fail("unexpected Weyl root-closure size for " + datum.name);
    }
    if (generated_roots.size() != listed_roots.size()) {
        fail("listed and generated root sets have different sizes");
    }
    for (const Weight<Rank>& root : generated_roots) {
        if (!listed_roots.contains(root)) {
            fail("listed adjoint weights omit a generated root");
        }
    }
}

template <int Rank, std::size_t WeightCount>
RootDatum<Rank> exceptional_datum(
    const std::string& name,
    const int (&cartan)[Rank][Rank],
    const int (&weights)[WeightCount][Rank],
    const int (&multiplicities)[WeightCount],
    int expected_root_count) {
    RootDatum<Rank> datum;
    datum.name = name;
    datum.expected_root_count = expected_root_count;
    for (int row = 0; row < Rank; ++row) {
        for (int column = 0; column < Rank; ++column) {
            datum.cartan[row][column] = cartan[row][column];
        }
    }
    datum.adjoint_weights.reserve(WeightCount);
    datum.adjoint_multiplicities.reserve(WeightCount);
    for (std::size_t index = 0; index < WeightCount; ++index) {
        Weight<Rank> weight;
        for (int coordinate = 0; coordinate < Rank; ++coordinate) {
            weight.coordinate[coordinate] = weights[index][coordinate];
        }
        datum.adjoint_weights.push_back(weight);
        datum.adjoint_multiplicities.push_back(multiplicities[index]);
    }
    return datum;
}

template <int Rank>
RootDatum<Rank> classical_datum(char family) {
    RootDatum<Rank> datum;
    datum.name = std::string(1, family) + "_" + std::to_string(Rank);
    for (int index = 0; index < Rank; ++index) datum.cartan[index][index] = 2;

    if (family == 'B' || family == 'C') {
        if constexpr (Rank >= 2) {
            for (int index = 0; index + 1 < Rank; ++index) {
                datum.cartan[index][index + 1] = -1;
                datum.cartan[index + 1][index] = -1;
            }
            if (family == 'B') {
                datum.cartan[Rank - 2][Rank - 1] = -2;
                datum.expected_root_count = 2 * Rank * Rank;
            } else {
                datum.cartan[Rank - 1][Rank - 2] = -2;
                datum.expected_root_count = 2 * Rank * Rank;
            }
        } else {
            fail("B/C rank must be at least two");
        }
    } else if (family == 'D') {
        if constexpr (Rank >= 4) {
            for (int index = 0; index + 1 < Rank - 2; ++index) {
                datum.cartan[index][index + 1] = -1;
                datum.cartan[index + 1][index] = -1;
            }
            datum.cartan[Rank - 3][Rank - 2] = -1;
            datum.cartan[Rank - 2][Rank - 3] = -1;
            datum.cartan[Rank - 3][Rank - 1] = -1;
            datum.cartan[Rank - 1][Rank - 3] = -1;
            datum.expected_root_count = 2 * Rank * (Rank - 1);
        } else {
            fail("D rank must be at least four");
        }
    } else {
        fail("unsupported classical family");
    }

    const WeightSet<Rank> roots = root_closure<Rank>(datum.cartan);
    datum.adjoint_weights.reserve(roots.size() + 1);
    datum.adjoint_multiplicities.reserve(roots.size() + 1);
    for (const Weight<Rank>& root : roots) {
        datum.adjoint_weights.push_back(root);
        datum.adjoint_multiplicities.push_back(1);
    }
    datum.adjoint_weights.emplace_back();
    datum.adjoint_multiplicities.push_back(Rank);
    return datum;
}

template <int Rank>
int dominant_reflect(
    Weight<Rank>& shifted_weight,
    const std::array<std::array<int, Rank>, Rank>& cartan) {
    int sign = 1;
    int reflection_count = 0;
    while (true) {
        int bad_coordinate = -1;
        for (int index = 0; index < Rank; ++index) {
            if (shifted_weight.coordinate[index] < 0) {
                bad_coordinate = index;
                break;
            }
            if (shifted_weight.coordinate[index] == 0) return 0;
        }
        if (bad_coordinate < 0) return sign;
        shifted_weight = reflect<Rank>(shifted_weight, bad_coordinate, cartan);
        sign = -sign;
        if (++reflection_count > 10000) fail("dominant reflection did not terminate");
    }
}

template <int Rank>
MultiplicityMap<Rank> tensor_with_adjoint(
    const MultiplicityMap<Rank>& source,
    const RootDatum<Rank>& datum,
    std::int32_t& maximum_coordinate) {
    MultiplicityMap<Rank> target;
    target.reserve(std::max<std::size_t>(64, 2 * source.size()));

    for (const auto& [highest_weight, coefficient] : source) {
        if (coefficient <= 0) fail("source decomposition has a nonpositive coefficient");
        for (std::size_t index = 0; index < datum.adjoint_weights.size(); ++index) {
            Weight<Rank> shifted;
            for (int coordinate = 0; coordinate < Rank; ++coordinate) {
                const std::int64_t value =
                    static_cast<std::int64_t>(highest_weight.coordinate[coordinate]) +
                    datum.adjoint_weights[index].coordinate[coordinate] + 1;
                if (value < std::numeric_limits<std::int32_t>::min() ||
                    value > std::numeric_limits<std::int32_t>::max()) {
                    fail("weight-coordinate overflow before dominant reflection");
                }
                shifted.coordinate[coordinate] = static_cast<std::int32_t>(value);
            }
            const int sign = dominant_reflect<Rank>(shifted, datum.cartan);
            if (sign == 0) continue;

            Weight<Rank> dominant;
            for (int coordinate = 0; coordinate < Rank; ++coordinate) {
                const std::int64_t value =
                    static_cast<std::int64_t>(shifted.coordinate[coordinate]) - 1;
                if (value < 0 || value > std::numeric_limits<std::int32_t>::max()) {
                    fail("dominant highest weight is outside the checked coordinate range");
                }
                dominant.coordinate[coordinate] = static_cast<std::int32_t>(value);
                maximum_coordinate = std::max(maximum_coordinate, dominant.coordinate[coordinate]);
            }

            mpz_class contribution = coefficient * datum.adjoint_multiplicities[index];
            if (sign < 0) contribution = -contribution;
            const auto [found, inserted] = target.try_emplace(dominant, contribution);
            if (!inserted) found->second += contribution;
        }
    }

    MultiplicityMap<Rank> positive;
    positive.reserve(target.size());
    for (auto& [weight, coefficient] : target) {
        if (coefficient < 0) fail("Racah--Speiser result has a negative multiplicity");
        if (coefficient != 0) positive.emplace(weight, std::move(coefficient));
    }
    return positive;
}

template <int Rank>
mpz_class character_inner_product(
    const MultiplicityMap<Rank>& left,
    const MultiplicityMap<Rank>& right) {
    const MultiplicityMap<Rank>* smaller = &left;
    const MultiplicityMap<Rank>* larger = &right;
    if (smaller->size() > larger->size()) std::swap(smaller, larger);
    mpz_class result = 0;
    for (const auto& [weight, coefficient] : *smaller) {
        const auto found = larger->find(weight);
        if (found != larger->end()) {
            mpz_addmul(result.get_mpz_t(), coefficient.get_mpz_t(),
                       found->second.get_mpz_t());
        }
    }
    return result;
}

template <int Rank>
int verify_group(
    RootDatum<Rank> datum,
    int maximum_moment,
    const std::vector<std::string>& source_paths,
    const std::string& stable_moment_path = "",
    char classical_family = '\0',
    int classical_rank = 0) {
    validate_root_datum(datum);
    std::map<int, mpz_class> expected;
    ClassicalClaims classical_claims;
    std::map<int, mpz_class> stable_moments;
    if (classical_family == '\0') {
        expected = read_expected_moments(source_paths);
        for (int degree = 0; degree <= maximum_moment; ++degree) {
            if (!expected.contains(degree)) {
                fail("expected sources omit m_" + std::to_string(degree));
            }
        }
        // A source ledger may retain a longer historical overlap than the
        // theorem consumes.  The command-line maximum is the fail-closed
        // proof obligation; discard surplus rows so the summary records the
        // exact prefix that was regenerated and compared.
        std::erase_if(expected, [maximum_moment](const auto& entry) {
            return entry.first > maximum_moment;
        });
    } else {
        if (stable_moment_path.empty()) fail("classical verification requires stable moments");
        classical_claims =
            read_classical_claims(source_paths, classical_family, classical_rank);
        stable_moments = read_stable_moments(stable_moment_path);
        for (const auto& [degree, _] : classical_claims.moments) {
            if (degree > maximum_moment) fail("classical moment claim exceeds MAX_MOMENT");
        }
        for (const auto& [degree, _] : classical_claims.deltas) {
            if (degree > maximum_moment) fail("classical Delta claim exceeds MAX_MOMENT");
            if (!stable_moments.contains(degree)) fail("stable source omits claimed degree");
        }
    }

    MultiplicityMap<Rank> current;
    Weight<Rank> zero;
    current.emplace(zero, mpz_class(1));
    std::vector<mpz_class> computed(maximum_moment + 1);
    std::int32_t maximum_coordinate = 0;
    const int maximum_power = (maximum_moment + 1) / 2;

    for (int power = 0; power <= maximum_power; ++power) {
        const auto started = Clock::now();
        const int even_degree = 2 * power;
        if (even_degree <= maximum_moment) {
            computed[even_degree] = character_inner_product(current, current);
        }
        if (power == maximum_power) break;

        MultiplicityMap<Rank> next = tensor_with_adjoint(current, datum, maximum_coordinate);
        const int odd_degree = 2 * power + 1;
        if (odd_degree <= maximum_moment) {
            computed[odd_degree] = character_inner_product(current, next);
        }
        const double seconds =
            std::chrono::duration<double>(Clock::now() - started).count();
        std::cout << "power=" << power + 1 << " support=" << next.size()
                  << " elapsed_seconds=" << seconds << '\n';
        current = std::move(next);
    }

    int failures = 0;
    if (classical_family == '\0') {
        for (int degree = 0; degree <= maximum_moment; ++degree) {
            if (computed[degree] != expected.at(degree)) {
                ++failures;
                std::cerr << "mismatch group=" << datum.name << " degree=" << degree
                          << " computed=" << computed[degree]
                          << " expected=" << expected.at(degree) << '\n';
            }
        }
    } else {
        for (const auto& [degree, value] : classical_claims.moments) {
            if (computed[degree] != value) {
                ++failures;
                std::cerr << "moment mismatch group=" << datum.name
                          << " degree=" << degree << " computed=" << computed[degree]
                          << " expected=" << value << '\n';
            }
        }
        for (const auto& [degree, value] : classical_claims.deltas) {
            const mpz_class delta = computed[degree] - stable_moments.at(degree);
            if (delta != value) {
                ++failures;
                std::cerr << "Delta mismatch group=" << datum.name
                          << " degree=" << degree << " computed=" << delta
                          << " expected=" << value << '\n';
            }
        }
    }
    std::cout << "root_datum=" << datum.name
              << " rank=" << Rank
              << " roots=" << datum.expected_root_count
              << " maximum_coordinate=" << maximum_coordinate << '\n';
    std::cout << "SUMMARY group=" << datum.name
              << " moments=" << maximum_moment + 1
              << " source_values="
              << (classical_family == '\0' ? expected.size() : classical_claims.size())
              << " failures=" << failures
              << " ignored_determinant_components="
              << (classical_family == '\0'
                      ? 0
                      : classical_claims.ignored_determinant_components)
              << " stable_recurrence_values_checked="
              << (classical_family == '\0' ? 0 : stable_moments.size())
              << '\n';
    return failures == 0 ? 0 : 1;
}

int dispatch(
    const std::string& group,
    int maximum_moment,
    const std::vector<std::string>& source_paths,
    const std::string& stable_moment_path) {
    if (group == "G2") {
        return verify_group(
            exceptional_datum("G2", lie_G2::cartan, lie_G2::adj_weight_coord,
                              lie_G2::adj_weight_mult, 12),
            maximum_moment, source_paths);
    }
    if (group == "F4") {
        return verify_group(
            exceptional_datum("F4", lie_F4::cartan, lie_F4::adj_weight_coord,
                              lie_F4::adj_weight_mult, 48),
            maximum_moment, source_paths);
    }
    if (group == "E6") {
        return verify_group(
            exceptional_datum("E6", lie_E6::cartan, lie_E6::adj_weight_coord,
                              lie_E6::adj_weight_mult, 72),
            maximum_moment, source_paths);
    }
    if (group == "E7") {
        return verify_group(
            exceptional_datum("E7", lie_E7::cartan, lie_E7::adj_weight_coord,
                              lie_E7::adj_weight_mult, 126),
            maximum_moment, source_paths);
    }
    if (group == "E8") {
        return verify_group(
            exceptional_datum("E8", lie_E8::cartan, lie_E8::adj_weight_coord,
                              lie_E8::adj_weight_mult, 240),
            maximum_moment, source_paths);
    }
    if (group.size() >= 2 && (group[0] == 'B' || group[0] == 'C' || group[0] == 'D')) {
        const char family = group[0];
        const int rank = std::stoi(group.substr(1));
#define CLASSICAL_DISPATCH(RANK)                                                \
        if (rank == RANK) {                                                     \
            return verify_group(classical_datum<RANK>(family), maximum_moment, \
                                source_paths, stable_moment_path, family, rank); \
        }
        CLASSICAL_DISPATCH(2)
        CLASSICAL_DISPATCH(3)
        CLASSICAL_DISPATCH(4)
        CLASSICAL_DISPATCH(5)
        CLASSICAL_DISPATCH(6)
        CLASSICAL_DISPATCH(7)
        CLASSICAL_DISPATCH(8)
        CLASSICAL_DISPATCH(9)
        CLASSICAL_DISPATCH(10)
        CLASSICAL_DISPATCH(11)
        CLASSICAL_DISPATCH(12)
        CLASSICAL_DISPATCH(13)
        CLASSICAL_DISPATCH(14)
        CLASSICAL_DISPATCH(15)
        CLASSICAL_DISPATCH(16)
        CLASSICAL_DISPATCH(17)
        CLASSICAL_DISPATCH(18)
        CLASSICAL_DISPATCH(19)
        CLASSICAL_DISPATCH(20)
        CLASSICAL_DISPATCH(21)
        CLASSICAL_DISPATCH(22)
        CLASSICAL_DISPATCH(23)
        CLASSICAL_DISPATCH(24)
#undef CLASSICAL_DISPATCH
        fail("unsupported classical rank: " + std::to_string(rank));
    }
    fail("unsupported group: " + group);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 4) {
            std::cerr << "usage: " << argv[0]
                      << " GROUP MAX_MOMENT SOURCE_LOG [SOURCE_LOG ...]\n";
            return 2;
        }
        const std::string group = argv[1];
        const int maximum_moment = std::stoi(argv[2]);
        if (maximum_moment < 0) fail("MAX_MOMENT must be nonnegative");
        std::vector<std::string> source_paths;
        std::string stable_moment_path;
        for (int index = 3; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--stable-moments") {
                if (++index >= argc) fail("--stable-moments requires a path");
                stable_moment_path = argv[index];
            } else {
                source_paths.push_back(argument);
            }
        }
        if (source_paths.empty()) fail("at least one source log is required");
        return dispatch(group, maximum_moment, source_paths, stable_moment_path);
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
