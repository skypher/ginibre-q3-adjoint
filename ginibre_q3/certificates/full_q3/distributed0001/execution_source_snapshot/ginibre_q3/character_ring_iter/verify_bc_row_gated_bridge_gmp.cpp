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
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <omp.h>

#include "ankerl/unordered_dense.h"

using BigInt = mpz_class;

namespace {

constexpr int kMaximumMoment = 63;
constexpr int kMaximumPrimeCount = 5;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

BigInt binomial(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    BigInt value = 1;
    for (int index = 1; index <= k; ++index) {
        value *= n - k + index;
        value /= index;
    }
    return value;
}

std::vector<BigInt> stable_moments() {
    std::vector<BigInt> stable(kMaximumMoment + 1);
    stable[0] = 1;
    stable[1] = 0;
    for (int n = 1; n < kMaximumMoment; ++n) {
        stable[n + 1] = BigInt(n) * stable[n] + BigInt(n) * stable[n - 1];
        if (n >= 2) {
            stable[n + 1] -= BigInt(n) * BigInt(n - 1) / 2 * stable[n - 2];
        }
    }
    return stable;
}

void compare_stable_table(
    const std::string& path,
    const std::vector<BigInt>& stable
) {
    std::ifstream input(path);
    if (!input) fail("could not open stable comparison table: " + path);
    std::vector<bool> seen(stable.size(), false);
    int degree = -1;
    std::string value;
    while (input >> degree >> value) {
        if (degree < 0 || degree >= static_cast<int>(stable.size())) continue;
        if (seen[degree]) fail("duplicate stable-table degree " + std::to_string(degree));
        seen[degree] = true;
        if (BigInt(value) != stable[degree]) {
            fail("stable recurrence/table mismatch in degree " + std::to_string(degree));
        }
    }
    if (!input.eof()) fail("malformed stable comparison table: " + path);
    for (int degree_index = 0; degree_index <= kMaximumMoment; ++degree_index) {
        if (!seen[degree_index]) {
            fail("stable comparison table omits degree " + std::to_string(degree_index));
        }
    }
}

std::array<std::uint64_t, kMaximumPrimeCount> make_primes() {
    std::array<std::uint64_t, kMaximumPrimeCount> primes{};
    BigInt start = BigInt(1) << 61;
    for (int index = 0; index < kMaximumPrimeCount; ++index) {
        BigInt candidate = start + BigInt(index) * 1000000 + 1009;
        BigInt prime;
        mpz_nextprime(prime.get_mpz_t(), candidate.get_mpz_t());
        if (mpz_probab_prime_p(prime.get_mpz_t(), 40) == 0) {
            fail("GMP rejected a generated CRT prime");
        }
        if (prime >= (BigInt(1) << 62)) fail("CRT prime exceeds the safe addition range");
        primes[index] = prime.get_ui();
        for (int previous = 0; previous < index; ++previous) {
            if (primes[previous] == primes[index]) fail("duplicate CRT prime");
        }
    }
    return primes;
}

template <int PrimeCount>
struct Residues {
    std::array<std::uint64_t, PrimeCount> value{};
};

template <int PrimeCount>
Residues<PrimeCount> residue_one() {
    Residues<PrimeCount> result;
    result.value.fill(1);
    return result;
}

template <int PrimeCount>
void add_residues(
    Residues<PrimeCount>& destination,
    const Residues<PrimeCount>& source,
    const std::array<std::uint64_t, kMaximumPrimeCount>& primes
) {
    for (int index = 0; index < PrimeCount; ++index) {
        // Every prime is below 2^62, so this sum is strictly below 2^63.
        std::uint64_t sum = destination.value[index] + source.value[index];
        if (sum >= primes[index]) sum -= primes[index];
        destination.value[index] = sum;
    }
}

template <int PrimeCount>
BigInt crt_reconstruct(
    const Residues<PrimeCount>& residues,
    const std::array<std::uint64_t, kMaximumPrimeCount>& primes
) {
    BigInt result = residues.value[0];
    BigInt modulus = primes[0];
    for (int index = 1; index < PrimeCount; ++index) {
        const std::uint64_t prime = primes[index];
        const std::uint64_t result_mod = mpz_fdiv_ui(result.get_mpz_t(), prime);
        const std::uint64_t modulus_mod = mpz_fdiv_ui(modulus.get_mpz_t(), prime);
        BigInt modulus_mod_big = modulus_mod;
        BigInt prime_big = prime;
        BigInt inverse;
        if (mpz_invert(
                inverse.get_mpz_t(),
                modulus_mod_big.get_mpz_t(),
                prime_big.get_mpz_t()) == 0) {
            fail("CRT moduli are not coprime");
        }
        BigInt difference = BigInt(residues.value[index]) - result_mod;
        difference %= prime_big;
        if (difference < 0) difference += prime_big;
        BigInt multiplier = difference * inverse % prime_big;
        result += modulus * multiplier;
        modulus *= prime_big;
    }
    for (int index = 0; index < PrimeCount; ++index) {
        if (mpz_fdiv_ui(result.get_mpz_t(), primes[index]) != residues.value[index]) {
            fail("internal CRT residue mismatch");
        }
    }
    return result;
}

template <int PrimeCount>
BigInt crt_modulus(
    const std::array<std::uint64_t, kMaximumPrimeCount>& primes
) {
    BigInt modulus = 1;
    for (int index = 0; index < PrimeCount; ++index) modulus *= primes[index];
    return modulus;
}

// A state stores a partition lambda in the common e_2^j Pieri expansion.
// If transposed=0, parts is lambda and length(lambda) is within the current
// B-height cap.  If transposed=1, parts is lambda' and lambda is outside that
// height cap but remains within the current C-width cap.  This canonical union
// represents exactly {length(lambda)<=2b_B+1 or lambda_1<=2b_C} without a
// long 2j-entry partition key.
template <int Capacity>
struct PackedState {
    std::array<std::uint8_t, Capacity> parts{};
    std::uint8_t length = 0;
    std::uint8_t transposed = 0;

    bool operator==(const PackedState& other) const {
        if (length != other.length || transposed != other.transposed) return false;
        for (int index = 0; index < length; ++index) {
            if (parts[index] != other.parts[index]) return false;
        }
        return true;
    }
};

template <int Capacity>
struct PackedStateHash {
    std::size_t operator()(const PackedState<Capacity>& state) const {
        std::uint64_t hash = 1469598103934665603ull;
        hash ^= static_cast<std::uint64_t>(state.length)
              | (static_cast<std::uint64_t>(state.transposed) << 8);
        hash *= 1099511628211ull;
        for (int index = 0; index < state.length; ++index) {
            hash ^= static_cast<std::uint64_t>(state.parts[index])
                  + 0x9e3779b97f4a7c15ull + (hash << 6) + (hash >> 2);
            hash *= 1099511628211ull;
        }
        return static_cast<std::size_t>(hash);
    }
};

template <int Capacity>
bool all_even_rows(const PackedState<Capacity>& state) {
    if (state.transposed) fail("all_even_rows called on a transposed state");
    for (int index = 0; index < state.length; ++index) {
        if (state.parts[index] % 2 != 0) return false;
    }
    return true;
}

template <int Capacity>
bool paired_columns(const PackedState<Capacity>& state) {
    if (!state.transposed) fail("paired_columns called on an ordinary state");
    if (state.length % 2 != 0) return false;
    for (int index = 0; index < state.length; index += 2) {
        if (state.parts[index] != state.parts[index + 1]) return false;
    }
    return true;
}

template <int Capacity>
PackedState<Capacity> transpose_ordinary(const PackedState<Capacity>& ordinary) {
    if (ordinary.transposed) fail("transpose_ordinary received a transposed state");
    PackedState<Capacity> result;
    result.transposed = 1;
    if (ordinary.length == 0) return result;
    const int width = ordinary.parts[0];
    if (width > Capacity) fail("transposed state exceeds packed capacity");
    result.length = static_cast<std::uint8_t>(width);
    for (int column = 1; column <= width; ++column) {
        int height = 0;
        while (height < ordinary.length && ordinary.parts[height] >= column) ++height;
        result.parts[column - 1] = static_cast<std::uint8_t>(height);
    }
    return result;
}

template <int Capacity>
bool canonicalize_ordinary(
    const PackedState<Capacity>& ordinary,
    int b_height_cap,
    int c_width_cap,
    PackedState<Capacity>& output
) {
    if (ordinary.transposed) fail("canonicalize_ordinary received a transposed state");
    if (ordinary.length <= b_height_cap) {
        output = ordinary;
        return true;
    }
    const int width = ordinary.length == 0 ? 0 : ordinary.parts[0];
    if (width > c_width_cap) return false;
    output = transpose_ordinary(ordinary);
    return true;
}

template <int Capacity>
bool row_can_receive_one(const PackedState<Capacity>& state, int row) {
    return row == 0 || state.parts[row - 1] > state.parts[row];
}

template <int Capacity, int PrimeCount, class Map>
void insert_or_add(
    Map& destination,
    PackedState<Capacity>&& state,
    const Residues<PrimeCount>& coefficient,
    const std::array<std::uint64_t, kMaximumPrimeCount>& primes
) {
    auto [iterator, inserted] = destination.try_emplace(std::move(state), coefficient);
    if (!inserted) add_residues(iterator->second, coefficient, primes);
}

template <int Capacity, int PrimeCount>
using StateMap = ankerl::unordered_dense::map<
    PackedState<Capacity>,
    Residues<PrimeCount>,
    PackedStateHash<Capacity>>;

template <int Capacity, int PrimeCount>
StateMap<Capacity, PrimeCount> pieri_step(
    const StateMap<Capacity, PrimeCount>& current,
    int b_height_cap,
    int c_width_cap,
    const std::array<std::uint64_t, kMaximumPrimeCount>& primes
) {
    StateMap<Capacity, PrimeCount> next;
    const std::size_t reserve_size = current.size() < 100000
        ? current.size() * 2 + 32
        : current.size() + current.size() / 2 + 32;
    next.reserve(reserve_size);

    for (const auto& [state, coefficient] : current) {
        if (!state.transposed) {
            const int width = state.length == 0 ? 0 : state.parts[0];
            if (state.length > b_height_cap && width > c_width_cap) continue;
            const int slots = static_cast<int>(state.length) + 2;
            if (slots > Capacity) fail("ordinary Pieri successor exceeds packed capacity");
            for (int first = 0; first < slots; ++first) {
                if (!row_can_receive_one(state, first)) continue;
                for (int second = first + 1; second < slots; ++second) {
                    if (second != first + 1 && !row_can_receive_one(state, second)) continue;
                    PackedState<Capacity> candidate = state;
                    candidate.parts[first] += 1;
                    candidate.parts[second] += 1;
                    candidate.length = static_cast<std::uint8_t>(
                        std::max<int>(candidate.length, second + 1));
                    PackedState<Capacity> canonical;
                    if (canonicalize_ordinary(
                            candidate, b_height_cap, c_width_cap, canonical)) {
                        insert_or_add<Capacity, PrimeCount>(
                            next, std::move(canonical), coefficient, primes);
                    }
                }
            }
        } else {
            // Conjugation sends a vertical two-strip to a horizontal two-strip.
            // Generate the latter directly in the short, width-capped key.
            if (state.length > c_width_cap) continue;
            const int old_length = state.length;

            // Two boxes in one row (and therefore in two distinct columns).
            for (int row = 0; row <= old_length; ++row) {
                const int old_part = row < old_length ? state.parts[row] : 0;
                if (row > 0 && state.parts[row - 1] < old_part + 2) continue;
                PackedState<Capacity> candidate = state;
                candidate.parts[row] += 2;
                candidate.length = static_cast<std::uint8_t>(
                    std::max<int>(candidate.length, row + 1));
                if (candidate.length > c_width_cap) continue;
                insert_or_add<Capacity, PrimeCount>(
                    next, std::move(candidate), coefficient, primes);
            }

            // One box in each of two rows.  The old row lengths must differ,
            // exactly expressing that the two new boxes occupy distinct columns.
            for (int first = 0; first < old_length; ++first) {
                if (first > 0 && state.parts[first - 1] == state.parts[first]) continue;
                for (int second = first + 1; second <= old_length; ++second) {
                    const int second_part = second < old_length ? state.parts[second] : 0;
                    if (state.parts[first] == second_part) continue;
                    if (second > first + 1
                        && state.parts[second - 1] == second_part) continue;
                    PackedState<Capacity> candidate = state;
                    candidate.parts[first] += 1;
                    candidate.parts[second] += 1;
                    candidate.length = static_cast<std::uint8_t>(
                        std::max<int>(candidate.length, second + 1));
                    if (candidate.length > c_width_cap) continue;
                    insert_or_add<Capacity, PrimeCount>(
                        next, std::move(candidate), coefficient, primes);
                }
            }
        }
    }
    return next;
}

struct RankWindow {
    int rank;
    int exact_through;
};

struct FamilyMomentTable {
    std::vector<BigInt> value;
    std::vector<bool> seen;
};

struct GroupResult {
    std::string name;
    std::vector<RankWindow> windows;
    std::map<int, FamilyMomentTable> b_moments;
    std::map<int, FamilyMomentTable> c_moments;
    std::size_t peak_states = 0;
    double elapsed_seconds = 0.0;
    int prime_count = 0;
    BigInt modulus;
    std::string progress_log;
};

int active_max_rank(const std::vector<RankWindow>& windows, int degree) {
    int result = -1;
    for (const RankWindow& window : windows) {
        if (window.exact_through >= degree) result = std::max(result, window.rank);
    }
    return result;
}

template <int Capacity, int PrimeCount>
GroupResult run_group(
    const std::string& name,
    const std::vector<RankWindow>& windows,
    const std::vector<BigInt>& stable,
    const std::array<std::uint64_t, kMaximumPrimeCount>& primes
) {
    GroupResult result;
    result.name = name;
    result.windows = windows;
    result.prime_count = PrimeCount;
    result.modulus = crt_modulus<PrimeCount>(primes);
    const int maximum_degree = std::max_element(
        windows.begin(), windows.end(),
        [](const RankWindow& left, const RankWindow& right) {
            return left.exact_through < right.exact_through;
        })->exact_through;
    const int minimum_rank = std::min_element(
        windows.begin(), windows.end(),
        [](const RankWindow& left, const RankWindow& right) {
            return left.rank < right.rank;
        })->rank;
    if (result.modulus <= stable[maximum_degree]) {
        fail(name + " CRT modulus does not dominate the stable moment bound");
    }

    for (const RankWindow& window : windows) {
        result.b_moments[window.rank] = FamilyMomentTable{
            std::vector<BigInt>(window.exact_through + 1),
            std::vector<bool>(window.exact_through + 1, false)};
        result.c_moments[window.rank] = FamilyMomentTable{
            std::vector<BigInt>(window.exact_through + 1),
            std::vector<bool>(window.exact_through + 1, false)};
    }

    StateMap<Capacity, PrimeCount> current;
    PackedState<Capacity> empty;
    current.emplace(empty, residue_one<PrimeCount>());
    result.peak_states = 1;
    const auto started = std::chrono::steady_clock::now();
    std::ostringstream progress;
    progress << "group=" << name
             << " prime_count=" << PrimeCount
             << " modulus_bits=" << mpz_sizeinbase(result.modulus.get_mpz_t(), 2)
             << " max_degree=" << maximum_degree << '\n';

    for (int degree = 1; degree <= maximum_degree; ++degree) {
        const int maximum_rank = active_max_rank(windows, degree);
        if (maximum_rank < 0) fail(name + " has no active rank at a required degree");
        const int b_height_cap = 2 * maximum_rank + 1;
        const int c_width_cap = 2 * maximum_rank;
        current = pieri_step<Capacity, PrimeCount>(
            current, b_height_cap, c_width_cap, primes);
        result.peak_states = std::max(result.peak_states, current.size());

        if (degree >= 2 * minimum_rank + 1) {
            std::vector<Residues<PrimeCount>> b_bins(maximum_rank + 1);
            std::vector<Residues<PrimeCount>> c_bins(maximum_rank + 1);
            Residues<PrimeCount> all_even_in_region;

            for (const auto& [state, coefficient] : current) {
                bool is_even_target = false;
                int original_height = 0;
                int original_width = 0;
                if (!state.transposed) {
                    if (all_even_rows(state)) {
                        is_even_target = true;
                        original_height = state.length;
                        original_width = state.length == 0 ? 0 : state.parts[0];
                    }
                } else if (paired_columns(state)) {
                    is_even_target = true;
                    original_height = state.length == 0 ? 0 : state.parts[0];
                    original_width = state.length;
                }
                if (!is_even_target) continue;
                add_residues(all_even_in_region, coefficient, primes);

                const int b_required = std::max(0, (original_height - 1 + 1) / 2);
                const int c_required = (original_width + 1) / 2;
                if (b_required <= maximum_rank) {
                    add_residues(b_bins[b_required], coefficient, primes);
                }
                if (c_required <= maximum_rank) {
                    add_residues(c_bins[c_required], coefficient, primes);
                }
            }

            for (int rank = 1; rank <= maximum_rank; ++rank) {
                add_residues(b_bins[rank], b_bins[rank - 1], primes);
                add_residues(c_bins[rank], c_bins[rank - 1], primes);
            }

            const bool region_is_all_partitions =
                b_height_cap + c_width_cap + 1 > 2 * degree;
            if (region_is_all_partitions) {
                BigInt total = crt_reconstruct(all_even_in_region, primes);
                if (total != stable[degree]) {
                    fail(name + " full-region stable target sum mismatch in degree "
                         + std::to_string(degree)
                         + ": states=" + std::to_string(current.size())
                         + ", computed=" + total.get_str()
                         + ", expected=" + stable[degree].get_str());
                }
            }

            for (const RankWindow& window : windows) {
                if (degree > window.exact_through) continue;
                BigInt b_value = crt_reconstruct(b_bins[window.rank], primes);
                BigInt c_value = crt_reconstruct(c_bins[window.rank], primes);
                if (b_value > stable[degree] || c_value > stable[degree]) {
                    fail(name + " reconstructed finite moment exceeds its stable bound");
                }
                const int boundary = 2 * window.rank + 2;
                if (degree < boundary
                    && (b_value != stable[degree] || c_value != stable[degree])) {
                    fail(name + " pre-boundary stabilization mismatch for rank "
                         + std::to_string(window.rank) + " in degree "
                         + std::to_string(degree));
                }
                auto& b_table = result.b_moments.at(window.rank);
                auto& c_table = result.c_moments.at(window.rank);
                b_table.value[degree] = std::move(b_value);
                c_table.value[degree] = std::move(c_value);
                b_table.seen[degree] = true;
                c_table.seen[degree] = true;
            }
        }

        progress << "group=" << name
                 << " degree=" << degree << '/' << maximum_degree
                 << " states=" << current.size()
                 << " b_cap=" << b_height_cap
                 << " c_cap=" << c_width_cap << '\n';
    }

    for (const RankWindow& window : windows) {
        const int boundary = 2 * window.rank + 2;
        for (int degree = boundary; degree <= window.exact_through; ++degree) {
            if (!result.b_moments.at(window.rank).seen[degree]
                || !result.c_moments.at(window.rank).seen[degree]) {
                fail(name + " omitted a required finite moment");
            }
        }
    }
    result.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    result.progress_log = progress.str();
    return result;
}

BigInt q3(const std::vector<BigInt>& moments, int n) {
    BigInt total = 0;
    for (int k = 0; k <= n; ++k) {
        total += binomial(n, k)
               * (moments[k + 2] * moments[n - k]
                  - moments[k + 1] * moments[n - k + 1]);
    }
    return 2 * total;
}

BigInt chain_difference(const std::vector<BigInt>& moments, int chain_m) {
    return q3(moments, 2 * chain_m + 3) - 4 * q3(moments, 2 * chain_m + 1);
}

std::vector<BigInt> chain_linear_coefficients(
    const std::vector<BigInt>& moments,
    int chain_m
) {
    const int maximum_index = 2 * chain_m + 5;
    std::vector<BigInt> coefficients(maximum_index + 1);
    for (const auto& [n, scale] :
         std::array<std::pair<int, int>, 2>{{
             {2 * chain_m + 3, 1},
             {2 * chain_m + 1, -4}}}) {
        for (int k = 0; k <= n; ++k) {
            BigInt coefficient = 2 * scale * binomial(n, k);
            int first = k + 2;
            int second = n - k;
            coefficients[first] += coefficient * moments[second];
            coefficients[second] += coefficient * moments[first];
            first = k + 1;
            second = n - k + 1;
            coefficients[first] -= coefficient * moments[second];
            coefficients[second] -= coefficient * moments[first];
        }
    }
    return coefficients;
}

std::map<std::pair<int, int>, BigInt> chain_quadratic_coefficients(int chain_m) {
    std::map<std::pair<int, int>, BigInt> coefficients;
    for (const auto& [n, scale] :
         std::array<std::pair<int, int>, 2>{{
             {2 * chain_m + 3, 1},
             {2 * chain_m + 1, -4}}}) {
        for (int k = 0; k <= n; ++k) {
            BigInt coefficient = 2 * scale * binomial(n, k);
            for (const auto& [first_raw, second_raw, sign] :
                 std::array<std::tuple<int, int, int>, 2>{{
                     {k + 2, n - k, 1},
                     {k + 1, n - k + 1, -1}}}) {
                int first = std::min(first_raw, second_raw);
                int second = std::max(first_raw, second_raw);
                coefficients[{first, second}] += sign * coefficient;
            }
        }
    }
    return coefficients;
}

struct IntegerInterval {
    BigInt lower;
    BigInt upper;
};

BigInt product_term_lower(
    const IntegerInterval& first,
    const IntegerInterval& second,
    const BigInt& coefficient
) {
    std::array<BigInt, 4> products{
        first.lower * second.lower,
        first.lower * second.upper,
        first.upper * second.lower,
        first.upper * second.upper};
    if (coefficient >= 0) {
        return coefficient * *std::min_element(products.begin(), products.end());
    }
    return coefficient * *std::max_element(products.begin(), products.end());
}

struct BridgeSummary {
    int rows = 0;
    int steps = 0;
    int exact_steps = 0;
    int bounded_steps = 0;
    int failures = 0;
    bool have_minimum = false;
    BigInt minimum_margin;
    std::string minimum_label;
    BigInt ledger_checksum = 0;
};

const GroupResult& group_for_rank(
    const std::vector<GroupResult>& groups,
    int rank
) {
    for (const GroupResult& group : groups) {
        for (const RankWindow& window : group.windows) {
            if (window.rank == rank) return group;
        }
    }
    fail("no generated group for rank " + std::to_string(rank));
}

int exact_through_for_rank(const std::vector<GroupResult>& groups, int rank) {
    if (rank >= 20) return 2 * rank + 1;
    const GroupResult& group = group_for_rank(groups, rank);
    for (const RankWindow& window : group.windows) {
        if (window.rank == rank) return window.exact_through;
    }
    fail("missing exact window for rank " + std::to_string(rank));
}

std::optional<BigInt> exact_delta(
    const std::vector<GroupResult>& groups,
    const std::vector<BigInt>& stable,
    char family,
    int rank,
    int degree
) {
    const int boundary = 2 * rank + 2;
    if (degree < boundary) return BigInt(0);
    if (rank >= 20) return std::nullopt;
    const GroupResult& group = group_for_rank(groups, rank);
    const int exact_through = exact_through_for_rank(groups, rank);
    if (degree > exact_through) return std::nullopt;
    const auto& tables = family == 'B' ? group.b_moments : group.c_moments;
    const FamilyMomentTable& table = tables.at(rank);
    if (degree >= static_cast<int>(table.seen.size()) || !table.seen[degree]) {
        fail("required generated moment is absent");
    }
    return table.value[degree] - stable[degree];
}

BridgeSummary verify_bridge(
    const std::vector<GroupResult>& groups,
    const std::vector<BigInt>& stable,
    bool print_steps
) {
    BridgeSummary summary;
    for (char family : std::array<char, 2>{'B', 'C'}) {
        for (int rank = 2; rank <= 30; ++rank) {
            ++summary.rows;
            const int exact_through = exact_through_for_rank(groups, rank);
            int row_exact = 0;
            int row_bounded = 0;
            BigInt row_minimum;
            int row_minimum_m = -1;
            bool row_has_minimum = false;

            for (int chain_m = rank - 1; chain_m <= 29; ++chain_m) {
                ++summary.steps;
                const int maximum_index = 2 * chain_m + 5;
                const BigInt stable_difference = chain_difference(stable, chain_m);
                if (stable_difference <= 0) fail("stable Chain difference is not positive");
                const std::vector<BigInt> linear =
                    chain_linear_coefficients(stable, chain_m);
                const auto quadratic = chain_quadratic_coefficients(chain_m);

                std::vector<IntegerInterval> delta_interval(maximum_index + 1);
                bool all_exact = true;
                for (int degree = 0; degree <= maximum_index; ++degree) {
                    std::optional<BigInt> delta =
                        exact_delta(groups, stable, family, rank, degree);
                    if (delta.has_value()) {
                        if (*delta > 0 || *delta < -stable[degree]) {
                            fail("generated correction violates -s_j <= delta_j <= 0");
                        }
                        delta_interval[degree] = IntegerInterval{*delta, *delta};
                        summary.ledger_checksum +=
                            BigInt((family == 'B' ? 1 : 2) * 1000000
                                   + rank * 1000 + degree)
                            * (*delta);
                    } else {
                        all_exact = false;
                        delta_interval[degree] = IntegerInterval{-stable[degree], 0};
                    }
                }

                BigInt margin = stable_difference;
                for (int degree = 0; degree <= maximum_index; ++degree) {
                    const IntegerInterval& interval = delta_interval[degree];
                    margin += linear[degree]
                        * (linear[degree] >= 0 ? interval.lower : interval.upper);
                }
                for (const auto& [indices, coefficient] : quadratic) {
                    margin += product_term_lower(
                        delta_interval[indices.first],
                        delta_interval[indices.second],
                        coefficient);
                }

                if (all_exact) {
                    std::vector<BigInt> row(stable.begin(), stable.begin() + maximum_index + 1);
                    for (int degree = 0; degree <= maximum_index; ++degree) {
                        row[degree] += delta_interval[degree].lower;
                    }
                    const BigInt direct = chain_difference(row, chain_m);
                    if (direct != margin) fail("exact expansion/direct Chain mismatch");
                    ++summary.exact_steps;
                    ++row_exact;
                } else {
                    ++summary.bounded_steps;
                    ++row_bounded;
                }

                if (margin <= 0) {
                    ++summary.failures;
                    fail(std::string(1, family) + "_" + std::to_string(rank)
                         + " m=" + std::to_string(chain_m)
                         + " has nonpositive lower margin " + margin.get_str());
                }
                const std::string label = std::string(1, family) + "_"
                    + std::to_string(rank) + " m=" + std::to_string(chain_m);
                if (!summary.have_minimum || margin < summary.minimum_margin) {
                    summary.have_minimum = true;
                    summary.minimum_margin = margin;
                    summary.minimum_label = label;
                }
                if (!row_has_minimum || margin < row_minimum) {
                    row_has_minimum = true;
                    row_minimum = margin;
                    row_minimum_m = chain_m;
                }
                if (print_steps) {
                    std::cout << "step family=" << family
                              << " rank=" << rank
                              << " m=" << chain_m
                              << " mode=" << (all_exact ? "exact" : "interval")
                              << " margin=" << margin << '\n';
                }
            }
            std::cout << "row family=" << family
                      << " rank=" << rank
                      << " m_range=" << rank - 1 << "..29"
                      << " exact_through=" << exact_through
                      << " exact_steps=" << row_exact
                      << " interval_steps=" << row_bounded
                      << " minimum_m=" << row_minimum_m
                      << " minimum_margin=" << row_minimum << '\n';
        }
    }
    return summary;
}

struct Options {
    int threads = 3;
    int diagnostic_max_degree = -1;
    bool print_steps = true;
    std::string stable_path = "../references/oeis_A002137_stable.txt";
    std::vector<std::string> selected_groups;
};

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--threads") {
            if (++index >= argc) fail("--threads requires an integer");
            options.threads = std::stoi(argv[index]);
            if (options.threads <= 0) fail("--threads must be positive");
        } else if (argument == "--stable") {
            if (++index >= argc) fail("--stable requires a path");
            options.stable_path = argv[index];
        } else if (argument == "--no-step-lines") {
            options.print_steps = false;
        } else if (argument == "--diagnostic-max-degree") {
            if (++index >= argc) fail("--diagnostic-max-degree requires an integer");
            options.diagnostic_max_degree = std::stoi(argv[index]);
            if (options.diagnostic_max_degree < 5
                || options.diagnostic_max_degree > kMaximumMoment) {
                fail("--diagnostic-max-degree must lie in 5..63");
            }
        } else if (argument == "--group") {
            if (++index >= argc) fail("--group requires low, mid, or high");
            options.selected_groups.emplace_back(argv[index]);
        } else if (argument == "--help") {
            std::cout
                << "usage: verify_bc_row_gated_bridge_gmp [options]\n"
                << "  --threads N        OpenMP group workers (default 3)\n"
                << "  --stable PATH      optional stable-table comparison\n"
                << "  --no-step-lines    suppress the 870 per-step ledger lines\n"
                << "  --group NAME       diagnostic group selection: low|mid|high\n"
                << "  --diagnostic-max-degree N  truncate selected group traversals\n";
            std::exit(0);
        } else {
            fail("unknown argument: " + argument);
        }
    }
    return options;
}

bool group_selected(const Options& options, const std::string& name) {
    if (options.selected_groups.empty()) return true;
    return std::find(options.selected_groups.begin(), options.selected_groups.end(), name)
        != options.selected_groups.end();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        const std::vector<BigInt> stable = stable_moments();
        compare_stable_table(options.stable_path, stable);
        const auto primes = make_primes();

        std::cout << "B/C row-gated bridge source verifier\n";
        std::cout << "arithmetic=exact modular Pieri plus GMP CRT"
                  << " stable_recurrence_values_checked=64\n";
        std::cout << "crt_primes=";
        for (int index = 0; index < kMaximumPrimeCount; ++index) {
            if (index) std::cout << ',';
            std::cout << primes[index];
        }
        std::cout << '\n';

        struct Task {
            std::string name;
            std::vector<RankWindow> windows;
        };
        std::vector<Task> tasks;
        if (group_selected(options, "low")) {
            tasks.push_back(Task{
                "low",
                {{2, 61}, {3, 61}, {4, 61}}});
        }
        if (group_selected(options, "mid")) {
            tasks.push_back(Task{
                "mid",
                {{5, 47}, {6, 45}, {7, 43}, {8, 43}}});
        }
        if (group_selected(options, "high")) {
            tasks.push_back(Task{
                "high",
                {{9, 41}, {10, 41}, {11, 41}, {12, 41}, {13, 41},
                 {14, 41}, {15, 41}, {16, 41}, {17, 41}, {18, 41}, {19, 41}}});
        }
        if (tasks.empty()) fail("no valid group selected");
        if (options.diagnostic_max_degree >= 0) {
            for (Task& task : tasks) {
                for (RankWindow& window : task.windows) {
                    window.exact_through = std::min(
                        window.exact_through, options.diagnostic_max_degree);
                }
            }
        }

        std::vector<GroupResult> groups(tasks.size());
        std::vector<std::string> task_errors(tasks.size());
        omp_set_dynamic(0);
        omp_set_num_threads(std::min<int>(options.threads, tasks.size()));
        #pragma omp parallel for schedule(dynamic, 1)
        for (int task_index = 0; task_index < static_cast<int>(tasks.size()); ++task_index) {
            try {
                const Task& task = tasks[task_index];
                if (task.name == "low") {
                    groups[task_index] = run_group<11, 5>(
                        task.name, task.windows, stable, primes);
                } else if (task.name == "mid") {
                    groups[task_index] = run_group<19, 4>(
                        task.name, task.windows, stable, primes);
                } else if (task.name == "high") {
                    groups[task_index] = run_group<41, 3>(
                        task.name, task.windows, stable, primes);
                } else {
                    fail("internal unknown group: " + task.name);
                }
            } catch (const std::exception& error) {
                task_errors[task_index] = error.what();
            }
        }
        for (int task_index = 0; task_index < static_cast<int>(tasks.size()); ++task_index) {
            if (!task_errors[task_index].empty()) {
                fail("group " + tasks[task_index].name + " failed: "
                     + task_errors[task_index]);
            }
        }

        std::sort(groups.begin(), groups.end(), [](const GroupResult& left, const GroupResult& right) {
            return left.windows.front().rank < right.windows.front().rank;
        });
        for (const GroupResult& group : groups) {
            std::cout << group.progress_log;
            for (const RankWindow& window : group.windows) {
                const int boundary = 2 * window.rank + 2;
                if (boundary <= window.exact_through) {
                    const auto& b_table = group.b_moments.at(window.rank);
                    const auto& c_table = group.c_moments.at(window.rank);
                    std::cout << "correction_window rank=" << window.rank
                              << " degrees=" << boundary << ".." << window.exact_through
                              << " B_first=" << b_table.value[boundary] - stable[boundary]
                              << " C_first=" << c_table.value[boundary] - stable[boundary]
                              << " B_last=" << b_table.value[window.exact_through]
                                                - stable[window.exact_through]
                              << " C_last=" << c_table.value[window.exact_through]
                                                - stable[window.exact_through]
                              << '\n';
                }
            }
            std::cout << "group_summary name=" << group.name
                      << " ranks=" << group.windows.front().rank << ".."
                      << group.windows.back().rank
                      << " peak_states=" << group.peak_states
                      << " status=PASS\n";
        }

        if (!options.selected_groups.empty() || options.diagnostic_max_degree >= 0) {
            std::cout << "diagnostic_groups=" << groups.size()
                      << " full_bridge_verification=SKIPPED status=PASS\n";
            return 0;
        }

        int correction_entries = 0;
        BigInt correction_checksum = 0;
        for (const GroupResult& group : groups) {
            for (const RankWindow& window : group.windows) {
                const int boundary = 2 * window.rank + 2;
                for (int degree = boundary; degree <= window.exact_through; ++degree) {
                    const BigInt b_delta =
                        group.b_moments.at(window.rank).value[degree] - stable[degree];
                    const BigInt c_delta =
                        group.c_moments.at(window.rank).value[degree] - stable[degree];
                    correction_checksum += BigInt(1000000 + window.rank * 1000 + degree)
                        * b_delta;
                    correction_checksum += BigInt(2000000 + window.rank * 1000 + degree)
                        * c_delta;
                    correction_entries += 2;
                }
            }
        }
        if (correction_entries != 832) fail("generated correction-entry count mismatch");
        std::cout << "generated_correction_entries=" << correction_entries
                  << " correction_checksum=" << correction_checksum << '\n';

        const BridgeSummary summary = verify_bridge(groups, stable, options.print_steps);
        if (summary.rows != 58 || summary.steps != 870
            || summary.exact_steps + summary.bounded_steps != 870
            || summary.failures != 0) {
            fail("B/C bridge coverage count mismatch");
        }
        std::cout << "minimum_positive_margin_label=" << summary.minimum_label
                  << " minimum_positive_margin=" << summary.minimum_margin << '\n';
        std::cout << "ledger_checksum=" << summary.ledger_checksum << '\n';
        std::cout << "families=2 ranks=58 row_gated_steps=870"
                  << " exact_steps=" << summary.exact_steps
                  << " interval_steps=" << summary.bounded_steps
                  << " failures=0 status=PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "verify_bc_row_gated_bridge_gmp: ERROR: "
                  << error.what() << '\n';
        return 1;
    }
}
