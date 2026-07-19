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
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "ankerl/unordered_dense.h"

using BigInt = mpz_class;
using Partition = std::vector<int>;

constexpr int kPackedPartitionMaxRows = 32;

struct PackedPartition {
    std::array<std::uint8_t, kPackedPartitionMaxRows> parts{};
    std::uint8_t length = 0;

    bool operator==(const PackedPartition& other) const {
        if (length != other.length) return false;
        for (int index = 0; index < length; ++index) {
            if (parts[index] != other.parts[index]) return false;
        }
        return true;
    }
};

struct PartitionHash {
    std::size_t operator()(const Partition& partition) const {
        std::uint64_t h = 1469598103934665603ull;
        for (int part : partition) {
            std::uint64_t x = static_cast<std::uint64_t>(part);
            h ^= x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
            h *= 1099511628211ull;
        }
        return static_cast<std::size_t>(h);
    }
};

struct PackedPartitionHash {
    std::size_t operator()(const PackedPartition& partition) const {
        std::uint64_t h = 1469598103934665603ull ^ partition.length;
        for (int index = 0; index < partition.length; ++index) {
            std::uint64_t x = static_cast<std::uint64_t>(partition.parts[index]);
            h ^= x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
            h *= 1099511628211ull;
        }
        return static_cast<std::size_t>(h);
    }
};

struct CoeffKey {
    int power;
    Partition partition;

    bool operator==(const CoeffKey& other) const {
        return power == other.power && partition == other.partition;
    }
};

struct CoeffKeyHash {
    std::size_t operator()(const CoeffKey& key) const {
        std::uint64_t h = 7809847782465536322ull ^ static_cast<std::uint64_t>(key.power);
        PartitionHash partition_hash;
        h ^= static_cast<std::uint64_t>(partition_hash(key.partition)) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
        return static_cast<std::size_t>(h);
    }
};

ankerl::unordered_dense::map<CoeffKey, BigInt, CoeffKeyHash> coeff_cache;
ankerl::unordered_dense::map<Partition, std::vector<Partition>, PartitionHash> predecessor_cache;
std::unordered_map<int, std::vector<Partition>> partition_cache;
bool use_predecessor_cache = true;
bool clear_coefficient_cache_after_moment = false;
bool use_forward_coefficients = false;
bool use_reverse_target_sums = false;
bool coefficient_progress = false;
std::uint64_t coefficient_cache_misses = 0;
bool coefficient_progress_initialized = false;
std::chrono::steady_clock::time_point last_coefficient_progress;

void clear_predecessor_cache_after_moment() {
    if (use_predecessor_cache) predecessor_cache.clear();
    if (clear_coefficient_cache_after_moment) coeff_cache.clear();
}

void partition_rec(int total, int max_part, Partition& current, std::vector<Partition>& out) {
    if (total == 0) {
        out.push_back(current);
        return;
    }
    int first_max = std::min(total, max_part);
    for (int first = first_max; first >= 1; --first) {
        current.push_back(first);
        partition_rec(total - first, std::min(first, total - first), current, out);
        current.pop_back();
    }
}

std::vector<Partition> integer_partitions(int total) {
    std::vector<Partition> out;
    Partition current;
    partition_rec(total, total, current, out);
    return out;
}

const std::vector<Partition>& integer_partitions_cached(int total) {
    auto found = partition_cache.find(total);
    if (found != partition_cache.end()) return found->second;
    auto inserted = partition_cache.emplace(total, integer_partitions(total));
    return inserted.first->second;
}

bool should_print_progress(std::size_t index, std::size_t total) {
    return index + 1 == total || (index + 1) % 1000 == 0;
}

bool should_print_forward_inner_progress(std::size_t processed, std::size_t total) {
    return total >= 500000 && processed % 500000 == 0;
}

bool partition_all_even(const Partition& partition) {
    for (int part : partition) {
        if (part % 2 != 0) return false;
    }
    return true;
}

bool partition_all_odd(const Partition& partition) {
    for (int part : partition) {
        if (part % 2 == 0) return false;
    }
    return true;
}

bool partition_is_nonincreasing(const Partition& partition) {
    for (std::size_t index = 1; index < partition.size(); ++index) {
        if (partition[index - 1] < partition[index]) return false;
    }
    return true;
}

Partition conjugate_partition(const Partition& partition) {
    Partition out;
    if (partition.empty()) return out;
    for (int column = 1; column <= partition.front(); ++column) {
        int count = 0;
        for (int row : partition) {
            if (row >= column) ++count;
        }
        out.push_back(count);
    }
    return out;
}

int partition_sum(const Partition& partition) {
    int total = 0;
    for (int row : partition) total += row;
    return total;
}

bool is_vertical_two_strip(const Partition& lower, const Partition& upper) {
    if (partition_sum(upper) - partition_sum(lower) != 2) return false;
    const std::size_t length = std::max(lower.size(), upper.size());
    for (std::size_t i = 0; i < length; ++i) {
        int lo = i < lower.size() ? lower[i] : 0;
        int up = i < upper.size() ? upper[i] : 0;
        if (up < lo) return false;
        if (up - lo != 0 && up - lo != 1) return false;
    }
    return true;
}

Partition normalize_partition(Partition partition) {
    Partition out;
    for (int row : partition) {
        if (row > 0) out.push_back(row);
    }
    std::sort(out.begin(), out.end(), std::greater<int>());
    return out;
}

const std::vector<Partition>& pieri_e2_predecessors(const Partition& partition) {
    auto cached = predecessor_cache.find(partition);
    if (cached != predecessor_cache.end()) return cached->second;

    ankerl::unordered_dense::set<Partition, PartitionHash> found;
    for (std::size_t first = 0; first < partition.size(); ++first) {
        for (std::size_t second = first + 1; second < partition.size(); ++second) {
            Partition predecessor = partition;
            predecessor[first] -= 1;
            predecessor[second] -= 1;
            Partition normalized = normalize_partition(std::move(predecessor));
            if (is_vertical_two_strip(normalized, partition)) {
                found.insert(std::move(normalized));
            }
        }
    }

    std::vector<Partition> out(found.begin(), found.end());
    std::sort(out.begin(), out.end(), [](const Partition& a, const Partition& b) {
        return a > b;
    });
    auto inserted = predecessor_cache.emplace(partition, std::move(out));
    return inserted.first->second;
}

std::vector<Partition> pieri_e2_predecessors_uncached(const Partition& partition) {
    ankerl::unordered_dense::set<Partition, PartitionHash> found;
    for (std::size_t first = 0; first < partition.size(); ++first) {
        for (std::size_t second = first + 1; second < partition.size(); ++second) {
            Partition predecessor = partition;
            predecessor[first] -= 1;
            predecessor[second] -= 1;
            Partition normalized = normalize_partition(std::move(predecessor));
            if (is_vertical_two_strip(normalized, partition)) {
                found.insert(std::move(normalized));
            }
        }
    }

    std::vector<Partition> out(found.begin(), found.end());
    std::sort(out.begin(), out.end(), [](const Partition& a, const Partition& b) {
        return a > b;
    });
    return out;
}

BigInt pieri_e2_coefficient(int power, const Partition& partition) {
    if (power == 0) return partition.empty() ? BigInt(1) : BigInt(0);
    if (partition_sum(partition) != 2 * power) return 0;
    CoeffKey key{power, partition};
    auto cached = coeff_cache.find(key);
    if (cached != coeff_cache.end()) return cached->second;
    if (coefficient_progress) {
        ++coefficient_cache_misses;
        auto now = std::chrono::steady_clock::now();
        if (!coefficient_progress_initialized) {
            coefficient_progress_initialized = true;
            last_coefficient_progress = now;
        } else {
            double seconds = std::chrono::duration<double>(now - last_coefficient_progress).count();
            if (seconds >= 10.0) {
                std::cout << "    coeff heartbeat misses=" << coefficient_cache_misses
                          << " cache=" << coeff_cache.size()
                          << " power=" << power
                          << " partition-size=" << partition.size()
                          << std::endl;
                last_coefficient_progress = now;
            }
        }
    }

    BigInt total = 0;
    std::vector<Partition> uncached_predecessors;
    const std::vector<Partition>* predecessors = nullptr;
    if (use_predecessor_cache) {
        predecessors = &pieri_e2_predecessors(partition);
    } else {
        uncached_predecessors = pieri_e2_predecessors_uncached(partition);
        predecessors = &uncached_predecessors;
    }
    for (const Partition& predecessor : *predecessors) {
        total += pieri_e2_coefficient(power - 1, predecessor);
    }
    auto inserted = coeff_cache.emplace(std::move(key), total);
    return inserted.first->second;
}

BigInt pieri_e2_coefficient_sum_reverse(
    int power,
    const std::vector<Partition>& targets,
    bool progress,
    const std::string& label
) {
    if (power == 0) {
        BigInt total = 0;
        for (const Partition& target : targets) {
            if (target.empty()) ++total;
        }
        return total;
    }

    ankerl::unordered_dense::map<Partition, BigInt, PartitionHash> current;
    for (const Partition& target : targets) current[target] += 1;
    auto start = std::chrono::steady_clock::now();
    for (int step = power; step >= 1; --step) {
        ankerl::unordered_dense::map<Partition, BigInt, PartitionHash> next;
        for (const auto& [partition, weight] : current) {
            for (const Partition& predecessor : pieri_e2_predecessors_uncached(partition)) {
                next[predecessor] += weight;
            }
        }
        current = std::move(next);
        if (progress) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    " << label
                      << " reverse step " << (power - step + 1) << "/" << power
                      << " frontier=" << current.size()
                      << " elapsed=" << seconds << "s" << std::endl;
        }
    }
    auto found = current.find(Partition{});
    if (found == current.end()) return 0;
    return found->second;
}

std::vector<BigInt> pieri_e2_coefficient_sums_reverse(
    int power,
    const std::vector<std::vector<Partition>>& target_sets,
    bool progress,
    const std::string& label
) {
    std::size_t color_count = target_sets.size();
    std::vector<BigInt> zero(color_count);
    if (power == 0) {
        std::vector<BigInt> totals(color_count);
        for (std::size_t color = 0; color < color_count; ++color) {
            for (const Partition& target : target_sets[color]) {
                if (target.empty()) totals[color] += 1;
            }
        }
        return totals;
    }

    ankerl::unordered_dense::map<Partition, std::vector<BigInt>, PartitionHash> current;
    for (std::size_t color = 0; color < color_count; ++color) {
        for (const Partition& target : target_sets[color]) {
            auto& weights = current[target];
            if (weights.empty()) weights = zero;
            weights[color] += 1;
        }
    }

    auto start = std::chrono::steady_clock::now();
    for (int step = power; step >= 1; --step) {
        ankerl::unordered_dense::map<Partition, std::vector<BigInt>, PartitionHash> next;
        for (const auto& [partition, weights] : current) {
            for (const Partition& predecessor : pieri_e2_predecessors_uncached(partition)) {
                auto& next_weights = next[predecessor];
                if (next_weights.empty()) next_weights = zero;
                for (std::size_t color = 0; color < color_count; ++color) {
                    if (weights[color] != 0) next_weights[color] += weights[color];
                }
            }
        }
        current = std::move(next);
        if (progress) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    " << label
                      << " shared reverse step " << (power - step + 1) << "/" << power
                      << " frontier=" << current.size()
                      << " colors=" << color_count
                      << " elapsed=" << seconds << "s" << std::endl;
        }
    }

    auto found = current.find(Partition{});
    if (found == current.end()) return zero;
    return found->second;
}

BigInt pieri_e2_coefficient_sum(
    int power,
    const std::vector<Partition>& targets,
    bool progress,
    const std::string& label
) {
    if (use_reverse_target_sums) {
        return pieri_e2_coefficient_sum_reverse(power, targets, progress, label);
    }
    BigInt total = 0;
    auto start = std::chrono::steady_clock::now();
    for (std::size_t index = 0; index < targets.size(); ++index) {
        total += pieri_e2_coefficient(power, targets[index]);
        if (progress && should_print_progress(index, targets.size())) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    " << label << " target "
                      << (index + 1) << "/" << targets.size()
                      << " elapsed=" << seconds << "s" << std::endl;
        }
    }
    return total;
}

std::vector<Partition> pieri_e2_successors(const Partition& lower) {
    ankerl::unordered_dense::set<Partition, PartitionHash> found;
    int slots = static_cast<int>(lower.size()) + 2;
    for (int first = 0; first < slots; ++first) {
        for (int second = first + 1; second < slots; ++second) {
            Partition candidate = lower;
            candidate.resize(slots, 0);
            candidate[first] += 1;
            candidate[second] += 1;
            candidate = normalize_partition(std::move(candidate));
            if (is_vertical_two_strip(lower, candidate)) {
                found.insert(std::move(candidate));
            }
        }
    }

    std::vector<Partition> out(found.begin(), found.end());
    std::sort(out.begin(), out.end(), [](const Partition& a, const Partition& b) {
        return a > b;
    });
    return out;
}

std::vector<Partition> pieri_e2_successors_row_capped_fast(
    const Partition& lower,
    int row_cap
) {
    std::vector<Partition> out;
    int slots = std::min(row_cap, static_cast<int>(lower.size()) + 2);
    out.reserve(static_cast<std::size_t>(slots * std::max(0, slots - 1) / 2));
    for (int first = 0; first < slots; ++first) {
        for (int second = first + 1; second < slots; ++second) {
            Partition candidate = lower;
            candidate.resize(std::max<std::size_t>(
                lower.size(),
                static_cast<std::size_t>(second + 1)
            ), 0);
            candidate[first] += 1;
            candidate[second] += 1;
            if (partition_is_nonincreasing(candidate)) {
                out.push_back(std::move(candidate));
            }
        }
    }
    return out;
}

ankerl::unordered_dense::map<Partition, BigInt, PartitionHash>
pieri_e2_forward_step(
    const ankerl::unordered_dense::map<Partition, BigInt, PartitionHash>& current,
    bool progress = false,
    const std::string& label = "forward",
    int step = 0,
    int total_steps = 0,
    std::chrono::steady_clock::time_point global_start = std::chrono::steady_clock::now()
) {
    ankerl::unordered_dense::map<Partition, BigInt, PartitionHash> next;
    next.reserve(current.size() * 2 + 16);
    std::size_t processed = 0;
    auto step_start = std::chrono::steady_clock::now();
    for (const auto& [partition, coefficient] : current) {
        for (Partition successor : pieri_e2_successors(partition)) {
            next[std::move(successor)] += coefficient;
        }
        ++processed;
        if (progress && should_print_forward_inner_progress(processed, current.size())) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - global_start).count();
            double step_elapsed = std::chrono::duration<double>(now - step_start).count();
            std::cout << "      " << label << " forward step "
                      << step << "/" << total_steps
                      << " processed=" << processed << "/" << current.size()
                      << " next_states=" << next.size()
                      << " step_elapsed=" << step_elapsed << "s"
                      << " elapsed=" << elapsed << "s" << std::endl;
        }
    }
    return next;
}

ankerl::unordered_dense::map<Partition, BigInt, PartitionHash>
pieri_e2_forward_step_row_capped(
    const ankerl::unordered_dense::map<Partition, BigInt, PartitionHash>& current,
    int row_cap,
    bool progress = false,
    const std::string& label = "forward capped",
    int step = 0,
    int total_steps = 0,
    std::chrono::steady_clock::time_point global_start = std::chrono::steady_clock::now()
) {
    ankerl::unordered_dense::map<Partition, BigInt, PartitionHash> next;
    next.reserve(current.size() * 2 + 16);
    std::size_t processed = 0;
    auto step_start = std::chrono::steady_clock::now();
    for (const auto& [partition, coefficient] : current) {
        for (Partition successor : pieri_e2_successors_row_capped_fast(partition, row_cap)) {
            next[std::move(successor)] += coefficient;
        }
        ++processed;
        if (progress && should_print_forward_inner_progress(processed, current.size())) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - global_start).count();
            double step_elapsed = std::chrono::duration<double>(now - step_start).count();
            std::cout << "      " << label << " step "
                      << step << "/" << total_steps
                      << " processed=" << processed << "/" << current.size()
                      << " next_states=" << next.size()
                      << " step_elapsed=" << step_elapsed << "s"
                      << " elapsed=" << elapsed << "s" << std::endl;
        }
    }
    return next;
}

bool packed_row_can_receive_with_left_neighbor(const PackedPartition& partition, int row) {
    return row == 0 || partition.parts[row - 1] > partition.parts[row];
}

bool packed_partition_all_even(const PackedPartition& partition) {
    for (int index = 0; index < partition.length; ++index) {
        if (partition.parts[index] % 2 != 0) return false;
    }
    return true;
}

bool packed_partition_full_all_odd(const PackedPartition& partition, int row_cap) {
    if (partition.length != row_cap) return false;
    for (int index = 0; index < partition.length; ++index) {
        if (partition.parts[index] % 2 == 0) return false;
    }
    return true;
}

ankerl::unordered_dense::map<PackedPartition, BigInt, PackedPartitionHash>
pieri_e2_forward_step_row_capped_packed(
    const ankerl::unordered_dense::map<PackedPartition, BigInt, PackedPartitionHash>& current,
    int row_cap,
    bool progress = false,
    const std::string& label = "forward capped",
    int step = 0,
    int total_steps = 0,
    std::chrono::steady_clock::time_point global_start = std::chrono::steady_clock::now()
) {
    if (row_cap > kPackedPartitionMaxRows) {
        throw std::runtime_error("packed row cap exceeds static capacity");
    }
    if (total_steps > std::numeric_limits<std::uint8_t>::max()) {
        throw std::runtime_error("packed partition part exceeds 8-bit capacity");
    }
    ankerl::unordered_dense::map<PackedPartition, BigInt, PackedPartitionHash> next;
    next.reserve(current.size() + current.size() / 3 + 16);
    std::size_t processed = 0;
    auto step_start = std::chrono::steady_clock::now();
    for (const auto& [partition, coefficient] : current) {
        int slots = std::min(row_cap, static_cast<int>(partition.length) + 2);
        for (int first = 0; first < slots; ++first) {
            if (!packed_row_can_receive_with_left_neighbor(partition, first)) continue;
            for (int second = first + 1; second < slots; ++second) {
                if (second != first + 1 &&
                    !packed_row_can_receive_with_left_neighbor(partition, second)) {
                    continue;
                }
                PackedPartition candidate = partition;
                int candidate_length = std::max<int>(candidate.length, second + 1);
                candidate.parts[first] += 1;
                candidate.parts[second] += 1;
                candidate.length = static_cast<std::uint8_t>(candidate_length);
                next[std::move(candidate)] += coefficient;
            }
        }
        ++processed;
        if (progress && should_print_forward_inner_progress(processed, current.size())) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - global_start).count();
            double step_elapsed = std::chrono::duration<double>(now - step_start).count();
            std::cout << "      " << label << " packed step "
                      << step << "/" << total_steps
                      << " processed=" << processed << "/" << current.size()
                      << " next_states=" << next.size()
                      << " step_elapsed=" << step_elapsed << "s"
                      << " elapsed=" << elapsed << "s" << std::endl;
        }
    }
    return next;
}

std::uint64_t add_mod_u63(std::uint64_t left, std::uint64_t right, std::uint64_t modulus) {
    std::uint64_t sum = left + right;
    if (sum >= modulus || sum < left) sum -= modulus;
    return sum;
}

ankerl::unordered_dense::map<PackedPartition, std::uint64_t, PackedPartitionHash>
pieri_e2_forward_step_row_capped_packed_mod(
    const ankerl::unordered_dense::map<PackedPartition, std::uint64_t, PackedPartitionHash>& current,
    int row_cap,
    std::uint64_t modulus,
    bool progress = false,
    const std::string& label = "forward capped mod",
    int step = 0,
    int total_steps = 0,
    std::chrono::steady_clock::time_point global_start = std::chrono::steady_clock::now()
) {
    if (row_cap > kPackedPartitionMaxRows) {
        throw std::runtime_error("packed row cap exceeds static capacity");
    }
    if (total_steps > std::numeric_limits<std::uint8_t>::max()) {
        throw std::runtime_error("packed partition part exceeds 8-bit capacity");
    }
    if (modulus >= (std::uint64_t{1} << 63)) {
        throw std::runtime_error("modular replay requires modulus below 2^63");
    }
    ankerl::unordered_dense::map<PackedPartition, std::uint64_t, PackedPartitionHash> next;
    next.reserve(current.size() + current.size() / 3 + 16);
    std::size_t processed = 0;
    auto step_start = std::chrono::steady_clock::now();
    for (const auto& [partition, coefficient] : current) {
        int slots = std::min(row_cap, static_cast<int>(partition.length) + 2);
        for (int first = 0; first < slots; ++first) {
            if (!packed_row_can_receive_with_left_neighbor(partition, first)) continue;
            for (int second = first + 1; second < slots; ++second) {
                if (second != first + 1 &&
                    !packed_row_can_receive_with_left_neighbor(partition, second)) {
                    continue;
                }
                PackedPartition candidate = partition;
                int candidate_length = std::max<int>(candidate.length, second + 1);
                candidate.parts[first] += 1;
                candidate.parts[second] += 1;
                candidate.length = static_cast<std::uint8_t>(candidate_length);
                auto [iterator, inserted] = next.emplace(std::move(candidate), coefficient);
                if (!inserted) iterator->second = add_mod_u63(iterator->second, coefficient, modulus);
            }
        }
        ++processed;
        if (progress && should_print_forward_inner_progress(processed, current.size())) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - global_start).count();
            double step_elapsed = std::chrono::duration<double>(now - step_start).count();
            std::cout << "      " << label << " packed step "
                      << step << "/" << total_steps
                      << " processed=" << processed << "/" << current.size()
                      << " next_states=" << next.size()
                      << " step_elapsed=" << step_elapsed << "s"
                      << " elapsed=" << elapsed << "s" << std::endl;
        }
    }
    return next;
}

ankerl::unordered_dense::map<Partition, BigInt, PartitionHash>
pieri_e2_forward_coefficients(int power, bool progress) {
    ankerl::unordered_dense::map<Partition, BigInt, PartitionHash> current;
    current.emplace(Partition{}, BigInt(1));
    auto start = std::chrono::steady_clock::now();
    for (int step = 0; step < power; ++step) {
        current = pieri_e2_forward_step(current, progress, "forward", step + 1, power, start);
        if (progress) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    forward step " << (step + 1) << "/" << power
                      << " states=" << current.size()
                      << " elapsed=" << seconds << "s" << std::endl;
        }
    }
    return current;
}

std::vector<Partition> b_boundary_targets(int rank, int moment_index) {
    int matrix_size = 2 * rank + 1;
    std::vector<Partition> targets;
    for (const Partition& half_partition : integer_partitions_cached(moment_index)) {
        if (static_cast<int>(half_partition.size()) <= matrix_size) continue;
        Partition target;
        target.reserve(half_partition.size());
        for (int part : half_partition) target.push_back(2 * part);
        targets.push_back(std::move(target));
    }
    return targets;
}

std::vector<Partition> o_even_boundary_targets(int rank, int moment_index) {
    int matrix_size = 2 * rank;
    std::vector<Partition> targets;
    for (const Partition& half_partition : integer_partitions_cached(moment_index)) {
        if (static_cast<int>(half_partition.size()) <= matrix_size) continue;
        Partition target;
        target.reserve(half_partition.size());
        for (int part : half_partition) target.push_back(2 * part);
        targets.push_back(std::move(target));
    }
    return targets;
}

std::vector<Partition> o_even_kept_targets(int rank, int moment_index) {
    int matrix_size = 2 * rank;
    std::vector<Partition> targets;
    for (const Partition& half_partition : integer_partitions_cached(moment_index)) {
        if (static_cast<int>(half_partition.size()) > matrix_size) continue;
        Partition target;
        target.reserve(half_partition.size());
        for (int part : half_partition) target.push_back(2 * part);
        targets.push_back(std::move(target));
    }
    return targets;
}

std::vector<Partition> c_boundary_targets(int rank, int moment_index) {
    std::vector<Partition> targets;
    for (const Partition& paired_partition : integer_partitions_cached(moment_index)) {
        if (static_cast<int>(paired_partition.size()) <= rank) continue;
        Partition target;
        target.reserve(2 * paired_partition.size());
        for (int part : paired_partition) {
            target.push_back(part);
            target.push_back(part);
        }
        targets.push_back(conjugate_partition(target));
    }
    return targets;
}

std::vector<Partition> bc_boundary_targets(char family, int rank, int moment_index) {
    if (family == 'B') return b_boundary_targets(rank, moment_index);
    if (family == 'C') return c_boundary_targets(rank, moment_index);
    throw std::runtime_error("unknown B/C family");
}

std::vector<BigInt> read_stable_moments(const std::string& path, int max_index) {
    std::string resolved_path = path;
    std::ifstream input(resolved_path);
    if (!input && path == "ginibre_q3/references/oeis_A002137_stable.txt") {
        resolved_path = "references/oeis_A002137_stable.txt";
        input.open(resolved_path);
    }
    if (!input) throw std::runtime_error("could not open stable moments file: " + path);
    std::vector<BigInt> moments(max_index + 1);
    std::vector<bool> seen(max_index + 1, false);
    int index;
    std::string value;
    while (input >> index >> value) {
        if (0 <= index && index <= max_index) {
            moments[index] = BigInt(value);
            seen[index] = true;
        }
    }
    for (int i = 0; i <= max_index; ++i) {
        if (!seen[i]) throw std::runtime_error("missing stable moment m_" + std::to_string(i));
    }
    return moments;
}

BigInt binom_int(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    BigInt ans = 1;
    for (int i = 1; i <= k; ++i) {
        ans *= n - k + i;
        ans /= i;
    }
    return ans;
}

BigInt q3(const std::vector<BigInt>& moments, int n) {
    BigInt total = 0;
    for (int k = 0; k <= n; ++k) {
        total += binom_int(n, k) * (moments[k + 2] * moments[n - k] - moments[k + 1] * moments[n - k + 1]);
    }
    return 2 * total;
}

BigInt chain_diff(const std::vector<BigInt>& moments, int chain_m) {
    return q3(moments, 2 * chain_m + 3) - 4 * q3(moments, 2 * chain_m + 1);
}

BigInt interval_product_lower(
    const std::vector<BigInt>& lower,
    const std::vector<BigInt>& upper,
    int first,
    int second,
    const BigInt& coefficient
) {
    std::array<BigInt, 4> products{
        lower[first] * lower[second],
        lower[first] * upper[second],
        upper[first] * lower[second],
        upper[first] * upper[second],
    };
    if (coefficient >= 0) return coefficient * *std::min_element(products.begin(), products.end());
    return coefficient * *std::max_element(products.begin(), products.end());
}

BigInt chain_diff_interval_lower(
    const std::vector<BigInt>& lower,
    const std::vector<BigInt>& upper,
    int chain_m
) {
    BigInt total = 0;
    for (const auto& [n_value, scale_int] :
         std::vector<std::pair<int, int>>{{2 * chain_m + 3, 1}, {2 * chain_m + 1, -4}}) {
        for (int k = 0; k <= n_value; ++k) {
            BigInt c = 2 * scale_int * binom_int(n_value, k);
            total += interval_product_lower(lower, upper, k + 2, n_value - k, c);
            total += interval_product_lower(lower, upper, k + 1, n_value - k + 1, -c);
        }
    }
    return total;
}

std::vector<BigInt> chain_diff_linear_coefficients(
    const std::vector<BigInt>& moments,
    int chain_m
) {
    int max_index = 2 * chain_m + 5;
    std::vector<BigInt> coeffs(max_index + 1);
    for (const auto& [n_value, scale_int] :
         std::vector<std::pair<int, int>>{{2 * chain_m + 3, 1}, {2 * chain_m + 1, -4}}) {
        for (int k = 0; k <= n_value; ++k) {
            BigInt c = 2 * scale_int * binom_int(n_value, k);
            int first = k + 2;
            int second = n_value - k;
            coeffs[first] += c * moments[second];
            coeffs[second] += c * moments[first];
            first = k + 1;
            second = n_value - k + 1;
            coeffs[first] -= c * moments[second];
            coeffs[second] -= c * moments[first];
        }
    }
    return coeffs;
}

std::map<std::pair<int, int>, BigInt> chain_diff_quadratic_coefficients(int chain_m) {
    std::map<std::pair<int, int>, BigInt> coeffs;
    for (const auto& [n_value, scale_int] :
         std::vector<std::pair<int, int>>{{2 * chain_m + 3, 1}, {2 * chain_m + 1, -4}}) {
        for (int k = 0; k <= n_value; ++k) {
            BigInt c = 2 * scale_int * binom_int(n_value, k);
            for (const auto& [first_raw, second_raw, sign_int] :
                 std::vector<std::tuple<int, int, int>>{
                     {k + 2, n_value - k, 1},
                     {k + 1, n_value - k + 1, -1},
                 }) {
                int first = std::min(first_raw, second_raw);
                int second = std::max(first_raw, second_raw);
                coeffs[{first, second}] += sign_int * c;
            }
        }
    }
    return coeffs;
}

BigInt b_boundary_correction(int rank, int moment_index, bool progress) {
    int matrix_size = 2 * rank + 1;
    int first_boundary = matrix_size + 1;
    if (moment_index < first_boundary) return 0;

    std::vector<Partition> targets = b_boundary_targets(rank, moment_index);

    BigInt total = 0;
    auto start = std::chrono::steady_clock::now();
    if (use_forward_coefficients) {
        auto coefficients = pieri_e2_forward_coefficients(moment_index, progress);
        for (std::size_t index = 0; index < targets.size(); ++index) {
            auto found = coefficients.find(targets[index]);
            if (found != coefficients.end()) total += found->second;
            if (progress && should_print_progress(index, targets.size())) {
                auto now = std::chrono::steady_clock::now();
                double seconds = std::chrono::duration<double>(now - start).count();
                std::cout << "    m_" << moment_index << ": B forward target "
                          << (index + 1) << "/" << targets.size()
                          << " elapsed=" << seconds << "s" << std::endl;
            }
        }
    } else {
        for (std::size_t index = 0; index < targets.size(); ++index) {
            total += pieri_e2_coefficient(moment_index, targets[index]);
            if (progress && should_print_progress(index, targets.size())) {
                auto now = std::chrono::steady_clock::now();
                double seconds = std::chrono::duration<double>(now - start).count();
                std::cout << "    m_" << moment_index << ": B target " << (index + 1)
                          << "/" << targets.size() << " elapsed=" << seconds << "s" << std::endl;
            }
        }
    }
    return -total;
}

BigInt o_even_boundary_correction(int rank, int moment_index, bool progress) {
    int matrix_size = 2 * rank;
    int first_boundary = matrix_size + 1;
    if (moment_index < first_boundary) return 0;

    std::vector<Partition> targets = o_even_boundary_targets(rank, moment_index);

    BigInt total = 0;
    auto start = std::chrono::steady_clock::now();
    if (use_forward_coefficients) {
        auto coefficients = pieri_e2_forward_coefficients(moment_index, progress);
        for (std::size_t index = 0; index < targets.size(); ++index) {
            auto found = coefficients.find(targets[index]);
            if (found != coefficients.end()) total += found->second;
            if (progress && should_print_progress(index, targets.size())) {
                auto now = std::chrono::steady_clock::now();
                double seconds = std::chrono::duration<double>(now - start).count();
                std::cout << "    m_" << moment_index << ": O-even forward target "
                          << (index + 1) << "/" << targets.size()
                          << " elapsed=" << seconds << "s" << std::endl;
            }
        }
    } else {
        for (std::size_t index = 0; index < targets.size(); ++index) {
            total += pieri_e2_coefficient(moment_index, targets[index]);
            if (progress && should_print_progress(index, targets.size())) {
                auto now = std::chrono::steady_clock::now();
                double seconds = std::chrono::duration<double>(now - start).count();
                std::cout << "    m_" << moment_index << ": O-even target "
                          << (index + 1) << "/" << targets.size()
                          << " elapsed=" << seconds << "s" << std::endl;
            }
        }
    }
    return -total;
}

BigInt o_even_boundary_correction_from_stable(
    int rank,
    int moment_index,
    const BigInt& stable_moment,
    bool progress
) {
    int matrix_size = 2 * rank;
    int first_boundary = matrix_size + 1;
    if (moment_index < first_boundary) return 0;

    std::vector<Partition> overflow_targets = o_even_boundary_targets(rank, moment_index);
    std::vector<Partition> kept_targets = o_even_kept_targets(rank, moment_index);
    if (kept_targets.size() < overflow_targets.size()) {
        BigInt total = 0;
        auto start = std::chrono::steady_clock::now();
        if (use_forward_coefficients) {
            auto coefficients = pieri_e2_forward_coefficients(moment_index, progress);
            for (std::size_t index = 0; index < kept_targets.size(); ++index) {
                auto found = coefficients.find(kept_targets[index]);
                if (found != coefficients.end()) total += found->second;
                if (progress && should_print_progress(index, kept_targets.size())) {
                    auto now = std::chrono::steady_clock::now();
                    double seconds = std::chrono::duration<double>(now - start).count();
                    std::cout << "    m_" << moment_index << ": O-even kept forward target "
                              << (index + 1) << "/" << kept_targets.size()
                              << " elapsed=" << seconds << "s" << std::endl;
                }
            }
        } else {
            total = pieri_e2_coefficient_sum(
                moment_index,
                kept_targets,
                progress,
                "m_" + std::to_string(moment_index) + ": O-even kept"
            );
        }
        return total - stable_moment;
    }

    BigInt total = 0;
    auto start = std::chrono::steady_clock::now();
    if (use_forward_coefficients) {
        auto coefficients = pieri_e2_forward_coefficients(moment_index, progress);
        for (std::size_t index = 0; index < overflow_targets.size(); ++index) {
            auto found = coefficients.find(overflow_targets[index]);
            if (found != coefficients.end()) total += found->second;
            if (progress && should_print_progress(index, overflow_targets.size())) {
                auto now = std::chrono::steady_clock::now();
                double seconds = std::chrono::duration<double>(now - start).count();
                std::cout << "    m_" << moment_index << ": O-even overflow forward target "
                          << (index + 1) << "/" << overflow_targets.size()
                          << " elapsed=" << seconds << "s" << std::endl;
            }
        }
    } else {
        total = pieri_e2_coefficient_sum(
            moment_index,
            overflow_targets,
            progress,
            "m_" + std::to_string(moment_index) + ": O-even overflow"
        );
    }
    return -total;
}

BigInt c_boundary_correction(int rank, int moment_index, bool progress) {
    int first_boundary = 2 * rank + 2;
    if (moment_index < first_boundary) return 0;

    std::vector<Partition> targets = c_boundary_targets(rank, moment_index);

    BigInt total = 0;
    auto start = std::chrono::steady_clock::now();
    if (use_forward_coefficients) {
        auto coefficients = pieri_e2_forward_coefficients(moment_index, progress);
        for (std::size_t index = 0; index < targets.size(); ++index) {
            auto found = coefficients.find(targets[index]);
            if (found != coefficients.end()) total += found->second;
            if (progress && should_print_progress(index, targets.size())) {
                auto now = std::chrono::steady_clock::now();
                double seconds = std::chrono::duration<double>(now - start).count();
                std::cout << "    m_" << moment_index << ": C forward target "
                          << (index + 1) << "/" << targets.size()
                          << " elapsed=" << seconds << "s" << std::endl;
            }
        }
    } else {
        for (std::size_t index = 0; index < targets.size(); ++index) {
            total += pieri_e2_coefficient(moment_index, targets[index]);
            if (progress && should_print_progress(index, targets.size())) {
                auto now = std::chrono::steady_clock::now();
                double seconds = std::chrono::duration<double>(now - start).count();
                std::cout << "    m_" << moment_index << ": C target " << (index + 1)
                          << "/" << targets.size() << " elapsed=" << seconds << "s" << std::endl;
            }
        }
    }
    return -total;
}

std::vector<Partition> d_determinant_targets(int rank, int depth) {
    if (depth < 0) return {};
    int matrix_size = 2 * rank;
    const std::vector<Partition>& parts = integer_partitions_cached(depth);
    std::vector<Partition> targets;
    targets.reserve(parts.size());
    for (const Partition& excess_partition : parts) {
        if (static_cast<int>(excess_partition.size()) > matrix_size) continue;
        Partition target;
        target.reserve(matrix_size);
        for (int row : excess_partition) target.push_back(1 + 2 * row);
        for (int i = static_cast<int>(excess_partition.size()); i < matrix_size; ++i) {
            target.push_back(1);
        }
        targets.push_back(std::move(target));
    }
    return targets;
}

BigInt d_determinant_correction_pieri(int rank, int depth, bool progress) {
    if (depth < 0) return 0;
    int power = rank + depth;
    std::vector<Partition> targets = d_determinant_targets(rank, depth);
    if (use_reverse_target_sums) {
        return pieri_e2_coefficient_sum_reverse(
            power,
            targets,
            progress,
            "depth " + std::to_string(depth) + ": determinant"
        );
    }
    BigInt total = 0;
    auto start = std::chrono::steady_clock::now();
    for (std::size_t index = 0; index < targets.size(); ++index) {
        total += pieri_e2_coefficient(power, targets[index]);
        if (progress && should_print_progress(index, targets.size())) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    depth " << depth << ": partition " << (index + 1)
                      << "/" << targets.size() << " elapsed=" << seconds << "s" << std::endl;
        }
    }
    return total;
}

BigInt d_determinant_correction(int rank, int moment_index, bool progress) {
    int depth = moment_index - rank;
    if (depth < 0) return 0;
    if (depth == 0) return 1;
    if (depth == 1) return binom_int(rank, 2);
    if (depth == 2) {
        BigInt r = rank;
        return r * (r + 1) * (r * r + r + 2) / 8;
    }
    if (depth == 3) {
        BigInt r = rank;
        return (r - 1) * (r + 2)
            * (r * r * r * r + 8 * r * r * r + 25 * r * r + 18 * r - 24)
            / 48;
    }
    if (depth == 4) {
        BigInt r = rank;
        return (r + 3)
            * (r * r * r * r * r * r * r
               + 17 * r * r * r * r * r * r
               + 103 * r * r * r * r * r
               + 203 * r * r * r * r
               - 104 * r * r * r
               - 460 * r * r
               - 48 * r
               + 384)
            / 384;
    }
    return d_determinant_correction_pieri(rank, depth, progress);
}

void run_d_delta_shared_range_certificate(
    int rank_low,
    int rank_high,
    int moment_low,
    int moment_high,
    bool progress
) {
    if (rank_low > rank_high) throw std::runtime_error("empty D shared delta rank range");
    if (moment_low > moment_high) throw std::runtime_error("empty D shared delta moment range");
    for (int moment = moment_low; moment <= moment_high; ++moment) {
        std::vector<BigInt> corrections(rank_high - rank_low + 1);
        std::vector<std::vector<Partition>> target_sets(rank_high - rank_low + 1);
        bool need_reverse = false;
        for (int rank = rank_low; rank <= rank_high; ++rank) {
            int offset = rank - rank_low;
            int depth = moment - rank;
            if (depth < 0) {
                corrections[offset] = 0;
            } else if (depth <= 4) {
                corrections[offset] = d_determinant_correction(rank, moment, false);
            } else {
                target_sets[offset] = d_determinant_targets(rank, depth);
                need_reverse = true;
            }
        }
        if (need_reverse) {
            std::vector<BigInt> reverse_sums = pieri_e2_coefficient_sums_reverse(
                moment,
                target_sets,
                progress,
                "m_" + std::to_string(moment) + ": determinant"
            );
            for (int rank = rank_low; rank <= rank_high; ++rank) {
                int offset = rank - rank_low;
                if (moment - rank > 4) corrections[offset] = reverse_sums[offset];
            }
        }
        for (int rank = rank_low; rank <= rank_high; ++rank) {
            int offset = rank - rank_low;
            if (corrections[offset] < 0) throw std::runtime_error("D determinant correction is negative");
            std::cout << "D_" << rank << " Delta_" << moment
                      << " = " << corrections[offset] << std::endl;
        }
    }
}

BigInt known_chain_correction(char family, int rank, int chain_m) {
    if (chain_m != 15) return 0;
    if (family == 'B' && rank == 15) return BigInt("-339264390238328549868");
    if (family == 'B' && rank == 14) return BigInt("-737698907593136969475080");
    if (family == 'B' && rank == 13) return BigInt("-635328139851004313067237604");
    if (family == 'B' && rank == 12) return BigInt("-256707126830279829538508902680");
    if (family == 'B' && rank == 11) return BigInt("-52619004517878560530326968247436");
    if (family == 'B' && rank == 10) return BigInt("-5647409803148711830444230484267320");
    if (family == 'B' && rank == 9) return BigInt("-317557978055891990572818575684268260");
    if (family == 'B' && rank == 8) return BigInt("-9121041006433172256172146313596666964");
    if (family == 'B' && rank == 7) return BigInt("-127249345089157838812431532950376551000");
    if (family == 'B' && rank == 6) return BigInt("-804447184340788017806785311495257336096");
    if (family == 'B' && rank == 5) return BigInt("-2205880266500709219447769304943864381112");
    if (family == 'B' && rank == 4) return BigInt("-3029283134234602422641171284939031479176");
    if (family == 'B' && rank == 3) return BigInt("-3129776772761581844184843476221342890344");
    if (family == 'B' && rank == 2) return BigInt("-3131022199950767133420123332518997494192");
    if (family == 'B' && rank == 16) return BigInt("-40215585765779526");
    if (family == 'C' && rank == 15) return BigInt("-10337678835881613524400902");
    if (family == 'C' && rank == 14) return BigInt("-4225520088202992397777700472");
    if (family == 'C' && rank == 13) return BigInt("-754561571025515760464808930956");
    if (family == 'C' && rank == 12) return BigInt("-70721890698322742456043741347784");
    if (family == 'C' && rank == 11) return BigInt("-3820811508088050424306952330577532");
    if (family == 'C' && rank == 10) return BigInt("-125013333708779322051787748009650640");
    if (family == 'C' && rank == 9) return BigInt("-2530540051838536162939935842522210948");
    if (family == 'C' && rank == 8) return BigInt("-31637675464076401675130110253262660880");
    if (family == 'C' && rank == 7) return BigInt("-237961480808888649646563447638142418680");
    if (family == 'C' && rank == 6) return BigInt("-1018016221430356705833863174993672657336");
    if (family == 'C' && rank == 5) return BigInt("-2323348962227461750000694619447848635912");
    if (family == 'C' && rank == 4) return BigInt("-3040952058298754574362610465587221330608");
    if (family == 'C' && rank == 3) return BigInt("-3129864620820508469906220324474542092136");
    if (family == 'C' && rank == 2) return BigInt("-3131022199950767133420123332518997494192");
    if (family == 'C' && rank == 16) return BigInt("-7092579055254392700002");
    return 0;
}

void run_boundary(
    char family,
    int rank,
    int chain_m,
    const std::string& stable_path,
    bool progress
) {
    int moment_max = 2 * chain_m + 5;
    std::vector<BigInt> row = read_stable_moments(stable_path, moment_max);
    BigInt stable_diff = chain_diff(row, chain_m);
    std::cout << family << " boundary m=" << chain_m << " rank-" << rank << " C++ certificate" << std::endl;
    std::cout << "moment range m_0..m_" << moment_max << std::endl;
    std::cout << "stable Chain diff D(" << chain_m << ") = " << stable_diff << std::endl;
    if (family == 'B') {
        std::cout << "begin B_" << rank << "=SO(" << (2 * rank + 1) << ")" << std::endl;
    } else if (family == 'C') {
        std::cout << "begin C_" << rank << "=Sp(" << rank << ")" << std::endl;
    } else {
        throw std::runtime_error("unknown family");
    }
    int first_boundary = 2 * rank + 2;
    for (int moment_index = first_boundary; moment_index <= moment_max; ++moment_index) {
        BigInt correction =
            family == 'B'
                ? b_boundary_correction(rank, moment_index, progress)
                : c_boundary_correction(rank, moment_index, progress);
        row[moment_index] += correction;
        clear_predecessor_cache_after_moment();
        std::cout << "  m_" << moment_index << " += " << correction << std::endl;
    }
    BigInt row_diff = chain_diff(row, chain_m);
    BigInt chain_correction = row_diff - stable_diff;
    if (family == 'B') {
        std::cout << "B_" << rank << "=SO(" << (2 * rank + 1) << ")" << std::endl;
    } else {
        std::cout << "C_" << rank << "=Sp(" << rank << ")" << std::endl;
    }
    std::cout << "  boundary correction in Chain diff = " << chain_correction << std::endl;
    std::cout << "  " << family << "_" << rank << " Chain diff D(" << chain_m << ") = " << row_diff << std::endl;
    BigInt expected = known_chain_correction(family, rank, chain_m);
    if (expected != 0 && chain_correction != expected) {
        throw std::runtime_error("unexpected known Chain correction");
    }
    if (row_diff <= 0) throw std::runtime_error("boundary Chain diff is not positive");
}

void run_d_boundary(
    int rank,
    int chain_m,
    const std::string& stable_path,
    bool progress
) {
    int moment_max = 2 * chain_m + 5;
    std::vector<BigInt> row = read_stable_moments(stable_path, moment_max);
    BigInt stable_diff = chain_diff(row, chain_m);
    std::cout << "D boundary m=" << chain_m << " rank-" << rank << " C++ certificate" << std::endl;
    std::cout << "moment range m_0..m_" << moment_max << std::endl;
    std::cout << "stable Chain diff D(" << chain_m << ") = " << stable_diff << std::endl;
    std::cout << "begin D_" << rank << "=SO(" << (2 * rank) << ")" << std::endl;
    BigInt chain_correction = 0;
    for (int moment_index = rank; moment_index <= moment_max; ++moment_index) {
        BigInt correction = d_determinant_correction(rank, moment_index, progress);
        row[moment_index] += correction;
        clear_predecessor_cache_after_moment();
        std::cout << "  m_" << moment_index << " += " << correction << std::endl;
    }
    BigInt row_diff = chain_diff(row, chain_m);
    chain_correction = row_diff - stable_diff;
    std::cout << "D_" << rank << "=SO(" << (2 * rank) << ")" << std::endl;
    std::cout << "  boundary correction in Chain diff = " << chain_correction << std::endl;
    std::cout << "  D_" << rank << " Chain diff D(" << chain_m << ") = " << row_diff << std::endl;
    if (row_diff <= 0) throw std::runtime_error("D boundary Chain diff is not positive");
}

void run_d_exact_boundary(
    int rank,
    int chain_m,
    const std::string& stable_path,
    bool progress
) {
    int moment_max = 2 * chain_m + 5;
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_max);
    std::vector<BigInt> row = stable;
    BigInt stable_diff = chain_diff(row, chain_m);
    std::cout << "D exact boundary m=" << chain_m << " rank-" << rank << " C++ certificate" << std::endl;
    std::cout << "moment range m_0..m_" << moment_max << std::endl;
    std::cout << "stable Chain diff D(" << chain_m << ") = " << stable_diff << std::endl;
    std::cout << "begin D_" << rank << "=SO(" << (2 * rank) << ")" << std::endl;
    for (int moment_index = rank; moment_index <= moment_max; ++moment_index) {
        BigInt even_overflow =
            o_even_boundary_correction_from_stable(rank, moment_index, stable[moment_index], progress);
        BigInt determinant = d_determinant_correction(rank, moment_index, progress);
        row[moment_index] += even_overflow + determinant;
        clear_predecessor_cache_after_moment();
        std::cout << "  m_" << moment_index
                  << " += O_even " << even_overflow
                  << " + det " << determinant
                  << "; Delta_" << moment_index
                  << " = " << (row[moment_index] - stable[moment_index])
                  << std::endl;
    }
    BigInt row_diff = chain_diff(row, chain_m);
    BigInt chain_correction = row_diff - stable_diff;
    std::cout << "D_" << rank << "=SO(" << (2 * rank) << ") exact" << std::endl;
    std::cout << "  exact boundary correction in Chain diff = " << chain_correction << std::endl;
    std::cout << "  D_" << rank << " exact Chain diff D(" << chain_m << ") = " << row_diff << std::endl;
    if (row_diff <= 0) throw std::runtime_error("D exact boundary Chain diff is not positive");
}

struct DExactSharedRow {
    int rank;
    std::vector<BigInt> moments;
};

void run_d_exact_boundary_shared_range(
    int low_rank,
    int high_rank,
    int chain_m,
    const std::string& stable_path,
    bool progress
) {
    if (low_rank > high_rank) throw std::runtime_error("empty D exact rank range");
    int moment_max = 2 * chain_m + 5;
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_max);
    BigInt stable_diff = chain_diff(stable, chain_m);
    std::vector<DExactSharedRow> rows;
    for (int rank = low_rank; rank <= high_rank; ++rank) {
        rows.push_back(DExactSharedRow{rank, stable});
    }

    std::cout << "D exact shared boundary m=" << chain_m
              << " ranks " << low_rank << ".." << high_rank
              << " C++ certificate" << std::endl;
    std::cout << "moment range m_0..m_" << moment_max << std::endl;
    std::cout << "stable Chain diff D(" << chain_m << ") = " << stable_diff << std::endl;
    for (const DExactSharedRow& row : rows) {
        std::cout << "begin D_" << row.rank << "=SO(" << (2 * row.rank) << ")" << std::endl;
    }

    for (int moment_index = low_rank; moment_index <= moment_max; ++moment_index) {
        if (progress) {
            std::cout << "shared D exact moment m_" << moment_index
                      << " ranks " << low_rank << ".." << high_rank
                      << std::endl;
        }
        for (DExactSharedRow& row : rows) {
            if (moment_index < row.rank) continue;
            BigInt even_overflow =
                o_even_boundary_correction_from_stable(
                    row.rank,
                    moment_index,
                    stable[moment_index],
                    progress
                );
            BigInt determinant = d_determinant_correction(row.rank, moment_index, progress);
            row.moments[moment_index] += even_overflow + determinant;
            std::cout << "  D_" << row.rank
                      << " m_" << moment_index
                      << " += O_even " << even_overflow
                      << " + det " << determinant
                      << "; Delta_" << moment_index
                      << " = " << (row.moments[moment_index] - stable[moment_index])
                      << std::endl;
        }
        clear_predecessor_cache_after_moment();
    }

    for (const DExactSharedRow& row : rows) {
        BigInt row_diff = chain_diff(row.moments, chain_m);
        BigInt chain_correction = row_diff - stable_diff;
        std::cout << "D_" << row.rank << "=SO(" << (2 * row.rank) << ") exact" << std::endl;
        std::cout << "  exact boundary correction in Chain diff = " << chain_correction << std::endl;
        std::cout << "  D_" << row.rank << " exact Chain diff D(" << chain_m << ") = " << row_diff << std::endl;
        if (row_diff <= 0) throw std::runtime_error("D exact boundary Chain diff is not positive");
    }
}

void run_d_exact_moment_shared_range(
    int low_rank,
    int high_rank,
    int moment_low,
    int moment_high,
    const std::string& stable_path,
    bool progress
) {
    if (low_rank > high_rank) throw std::runtime_error("empty D exact rank range");
    if (moment_low > moment_high) throw std::runtime_error("empty D exact moment range");
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_high);

    std::cout << "D exact shared moment range ranks "
              << low_rank << ".." << high_rank
              << " moments m_" << moment_low << "..m_" << moment_high
              << " C++ certificate" << std::endl;
    for (int rank = low_rank; rank <= high_rank; ++rank) {
        std::cout << "begin D_" << rank << "=SO(" << (2 * rank) << ")" << std::endl;
    }

    for (int moment_index = moment_low; moment_index <= moment_high; ++moment_index) {
        if (progress) {
            std::cout << "shared D exact moment m_" << moment_index
                      << " ranks " << low_rank << ".." << high_rank
                      << std::endl;
        }
        std::vector<BigInt> even_sums(high_rank - low_rank + 1);
        std::vector<std::vector<Partition>> even_target_sets(high_rank - low_rank + 1);
        std::vector<bool> even_uses_kept(high_rank - low_rank + 1, false);
        bool need_even_reverse = false;
        for (int rank = low_rank; rank <= high_rank; ++rank) {
            int offset = rank - low_rank;
            int matrix_size = 2 * rank;
            int first_boundary = matrix_size + 1;
            if (moment_index < first_boundary) continue;
            std::vector<Partition> overflow_targets = o_even_boundary_targets(rank, moment_index);
            std::vector<Partition> kept_targets = o_even_kept_targets(rank, moment_index);
            if (kept_targets.size() < overflow_targets.size()) {
                even_target_sets[offset] = std::move(kept_targets);
                even_uses_kept[offset] = true;
            } else {
                even_target_sets[offset] = std::move(overflow_targets);
            }
            need_even_reverse = true;
        }
        if (need_even_reverse) {
            if (use_reverse_target_sums) {
                even_sums = pieri_e2_coefficient_sums_reverse(
                    moment_index,
                    even_target_sets,
                    progress,
                    "m_" + std::to_string(moment_index) + ": O-even"
                );
            } else {
                for (int rank = low_rank; rank <= high_rank; ++rank) {
                    int offset = rank - low_rank;
                    even_sums[offset] = pieri_e2_coefficient_sum(
                        moment_index,
                        even_target_sets[offset],
                        progress,
                        "m_" + std::to_string(moment_index) + ": O-even"
                    );
                }
            }
        }

        std::vector<BigInt> determinant_sums(high_rank - low_rank + 1);
        std::vector<std::vector<Partition>> determinant_target_sets(high_rank - low_rank + 1);
        bool need_determinant_reverse = false;
        for (int rank = low_rank; rank <= high_rank; ++rank) {
            int offset = rank - low_rank;
            int depth = moment_index - rank;
            if (depth < 0) {
                determinant_sums[offset] = 0;
            } else if (depth <= 4) {
                determinant_sums[offset] = d_determinant_correction(rank, moment_index, false);
            } else {
                determinant_target_sets[offset] = d_determinant_targets(rank, depth);
                need_determinant_reverse = true;
            }
        }
        if (need_determinant_reverse) {
            if (use_reverse_target_sums) {
                std::vector<BigInt> reverse_sums = pieri_e2_coefficient_sums_reverse(
                    moment_index,
                    determinant_target_sets,
                    progress,
                    "m_" + std::to_string(moment_index) + ": determinant"
                );
                for (int rank = low_rank; rank <= high_rank; ++rank) {
                    int offset = rank - low_rank;
                    if (moment_index - rank > 4) determinant_sums[offset] = reverse_sums[offset];
                }
            } else {
                for (int rank = low_rank; rank <= high_rank; ++rank) {
                    int offset = rank - low_rank;
                    if (moment_index - rank <= 4) continue;
                    determinant_sums[offset] = d_determinant_correction_pieri(
                        rank,
                        moment_index - rank,
                        progress
                    );
                }
            }
        }
        for (int rank = low_rank; rank <= high_rank; ++rank) {
            if (moment_index < rank) continue;
            int offset = rank - low_rank;
            BigInt even_overflow = 0;
            int matrix_size = 2 * rank;
            int first_boundary = matrix_size + 1;
            if (moment_index >= first_boundary) {
                if (even_uses_kept[offset]) {
                    even_overflow = even_sums[offset] - stable[moment_index];
                } else {
                    even_overflow = -even_sums[offset];
                }
            }
            BigInt determinant = determinant_sums[offset];
            BigInt delta = even_overflow + determinant;
            std::cout << "  D_" << rank
                      << " m_" << moment_index
                      << " += O_even " << even_overflow
                      << " + det " << determinant
                      << "; Delta_" << moment_index
                      << " = " << delta
                      << "; moment_" << moment_index
                      << " = " << (stable[moment_index] + delta)
                      << std::endl;
        }
        clear_predecessor_cache_after_moment();
    }
}

void run_d_o_even_moment_shared_range(
    int low_rank,
    int high_rank,
    int moment_low,
    int moment_high,
    const std::string& stable_path,
    bool progress
) {
    if (low_rank > high_rank) throw std::runtime_error("empty D O-even rank range");
    if (moment_low > moment_high) throw std::runtime_error("empty D O-even moment range");
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_high);

    std::cout << "D O-even shared moment range ranks "
              << low_rank << ".." << high_rank
              << " moments m_" << moment_low << "..m_" << moment_high
              << " C++ certificate" << std::endl;
    for (int rank = low_rank; rank <= high_rank; ++rank) {
        std::cout << "begin D_" << rank << "=SO(" << (2 * rank) << ")" << std::endl;
    }

    for (int moment_index = moment_low; moment_index <= moment_high; ++moment_index) {
        if (progress) {
            std::cout << "shared D O-even moment m_" << moment_index
                      << " ranks " << low_rank << ".." << high_rank
                      << std::endl;
        }
        std::vector<BigInt> even_sums(high_rank - low_rank + 1);
        std::vector<std::vector<Partition>> even_target_sets(high_rank - low_rank + 1);
        std::vector<bool> even_uses_kept(high_rank - low_rank + 1, false);
        bool need_even_reverse = false;
        for (int rank = low_rank; rank <= high_rank; ++rank) {
            int offset = rank - low_rank;
            int matrix_size = 2 * rank;
            int first_boundary = matrix_size + 1;
            if (moment_index < first_boundary) continue;
            std::vector<Partition> overflow_targets = o_even_boundary_targets(rank, moment_index);
            std::vector<Partition> kept_targets = o_even_kept_targets(rank, moment_index);
            if (kept_targets.size() < overflow_targets.size()) {
                even_target_sets[offset] = std::move(kept_targets);
                even_uses_kept[offset] = true;
            } else {
                even_target_sets[offset] = std::move(overflow_targets);
            }
            need_even_reverse = true;
        }
        if (need_even_reverse) {
            if (use_reverse_target_sums) {
                even_sums = pieri_e2_coefficient_sums_reverse(
                    moment_index,
                    even_target_sets,
                    progress,
                    "m_" + std::to_string(moment_index) + ": O-even"
                );
            } else {
                for (int rank = low_rank; rank <= high_rank; ++rank) {
                    int offset = rank - low_rank;
                    even_sums[offset] = pieri_e2_coefficient_sum(
                        moment_index,
                        even_target_sets[offset],
                        progress,
                        "m_" + std::to_string(moment_index) + ": O-even"
                    );
                }
            }
        }
        for (int rank = low_rank; rank <= high_rank; ++rank) {
            if (moment_index < rank) continue;
            int offset = rank - low_rank;
            BigInt even_overflow = 0;
            int matrix_size = 2 * rank;
            int first_boundary = matrix_size + 1;
            if (moment_index >= first_boundary) {
                if (even_uses_kept[offset]) {
                    even_overflow = even_sums[offset] - stable[moment_index];
                } else {
                    even_overflow = -even_sums[offset];
                }
            }
            std::cout << "D_" << rank
                      << " O_even_" << moment_index
                      << " = " << even_overflow
                      << std::endl;
        }
        clear_predecessor_cache_after_moment();
    }
}

void run_d_capped_suffix_moment_range(
    int rank,
    int moment_low,
    int moment_high,
    const std::string& stable_path,
    bool progress
) {
    if (moment_low > moment_high) throw std::runtime_error("empty D capped suffix moment range");
    int row_cap = 2 * rank;
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_high);
    ankerl::unordered_dense::map<PackedPartition, BigInt, PackedPartitionHash> current;
    current.emplace(PackedPartition{}, BigInt(1));

    std::cout << "D capped suffix moment range rank "
              << rank << " moments m_" << moment_low << "..m_" << moment_high
              << " row cap " << row_cap << " C++ certificate" << std::endl;
    std::cout << "begin D_" << rank << "=SO(" << row_cap << ")" << std::endl;

    auto start = std::chrono::steady_clock::now();
    for (int moment_index = 1; moment_index <= moment_high; ++moment_index) {
        current = pieri_e2_forward_step_row_capped_packed(
            current,
            row_cap,
            progress,
            "D_" + std::to_string(rank) + " capped forward",
            moment_index,
            moment_high,
            start
        );
        if (progress) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    D_" << rank << " capped forward step "
                      << moment_index << "/" << moment_high
                      << " states=" << current.size()
                      << " elapsed=" << seconds << "s" << std::endl;
        }
        if (moment_index < moment_low) continue;

        BigInt kept_even = 0;
        BigInt determinant = 0;
        std::size_t scanned = 0;
        auto scan_start = std::chrono::steady_clock::now();
        for (const auto& [partition, coefficient] : current) {
            if (packed_partition_all_even(partition)) kept_even += coefficient;
            if (moment_index >= rank && packed_partition_full_all_odd(partition, row_cap)) {
                determinant += coefficient;
            }
            ++scanned;
            if (progress && should_print_forward_inner_progress(scanned, current.size())) {
                auto now = std::chrono::steady_clock::now();
                double seconds = std::chrono::duration<double>(now - scan_start).count();
                std::cout << "      D_" << rank << " capped scan m_"
                          << moment_index
                          << " processed=" << scanned << "/" << current.size()
                          << " elapsed=" << seconds << "s" << std::endl;
            }
        }

        BigInt o_even = kept_even - stable[moment_index];
        BigInt delta = o_even + determinant;
        std::cout << "D_" << rank
                  << " O_even_" << moment_index
                  << " = " << o_even
                  << std::endl;
        std::cout << "D_" << rank
                  << " Delta_" << moment_index
                  << " = " << determinant
                  << std::endl;
        std::cout << "  D_" << rank
                  << " m_" << moment_index
                  << " += O_even " << o_even
                  << " + det " << determinant
                  << "; Delta_" << moment_index
                  << " = " << delta
                  << "; moment_" << moment_index
                  << " = " << (stable[moment_index] + delta)
                  << std::endl;
    }
}

void run_d_capped_suffix_moment_range_mod_prime(
    int rank,
    int moment_low,
    int moment_high,
    std::uint64_t modulus,
    bool progress
) {
    if (moment_low > moment_high) throw std::runtime_error("empty D capped suffix moment range");
    if (modulus < 3 || modulus % 2 == 0) throw std::runtime_error("modulus must be an odd integer");
    int row_cap = 2 * rank;
    ankerl::unordered_dense::map<PackedPartition, std::uint64_t, PackedPartitionHash> current;
    current.emplace(PackedPartition{}, 1 % modulus);

    std::cout << "D capped suffix modular moment range rank "
              << rank << " moments m_" << moment_low << "..m_" << moment_high
              << " row cap " << row_cap << " modulus " << modulus
              << " C++ certificate" << std::endl;
    std::cout << "begin D_" << rank << "=SO(" << row_cap << ") mod " << modulus << std::endl;

    auto start = std::chrono::steady_clock::now();
    for (int moment_index = 1; moment_index <= moment_high; ++moment_index) {
        current = pieri_e2_forward_step_row_capped_packed_mod(
            current,
            row_cap,
            modulus,
            progress,
            "D_" + std::to_string(rank) + " capped mod " + std::to_string(modulus),
            moment_index,
            moment_high,
            start
        );
        if (progress) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    D_" << rank << " capped mod " << modulus
                      << " forward step " << moment_index << "/" << moment_high
                      << " states=" << current.size()
                      << " elapsed=" << seconds << "s" << std::endl;
        }
        if (moment_index < moment_low) continue;

        std::uint64_t kept_even = 0;
        std::uint64_t determinant = 0;
        std::size_t scanned = 0;
        auto scan_start = std::chrono::steady_clock::now();
        for (const auto& [partition, coefficient] : current) {
            if (packed_partition_all_even(partition)) {
                kept_even = add_mod_u63(kept_even, coefficient, modulus);
            }
            if (moment_index >= rank && packed_partition_full_all_odd(partition, row_cap)) {
                determinant = add_mod_u63(determinant, coefficient, modulus);
            }
            ++scanned;
            if (progress && should_print_forward_inner_progress(scanned, current.size())) {
                auto now = std::chrono::steady_clock::now();
                double seconds = std::chrono::duration<double>(now - scan_start).count();
                std::cout << "      D_" << rank << " capped mod scan m_"
                          << moment_index
                          << " processed=" << scanned << "/" << current.size()
                          << " elapsed=" << seconds << "s" << std::endl;
            }
        }

        std::cout << "D_" << rank
                  << " MOD " << modulus
                  << " m_" << moment_index
                  << " kept_even=" << kept_even
                  << " det=" << determinant
                  << std::endl;
    }
}

void run_d_interval_boundary(
    int rank,
    int chain_m,
    const std::string& stable_path,
    bool progress
) {
    int moment_max = 2 * chain_m + 5;
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_max);
    std::vector<BigInt> determinant(moment_max + 1);
    std::cout << "D interval boundary m=" << chain_m << " rank-" << rank << " C++ certificate" << std::endl;
    std::cout << "moment range m_0..m_" << moment_max << std::endl;
    std::cout << "begin D_" << rank << "=SO(" << (2 * rank) << ")" << std::endl;
    for (int moment_index = rank; moment_index <= moment_max; ++moment_index) {
        determinant[moment_index] = d_determinant_correction(rank, moment_index, progress);
        clear_predecessor_cache_after_moment();
        std::cout << "  det_" << moment_index << " = " << determinant[moment_index] << std::endl;
    }

    std::vector<BigInt> lower(moment_max + 1);
    std::vector<BigInt> upper(moment_max + 1);
    for (int index = 0; index <= moment_max; ++index) {
        lower[index] = determinant[index];
        upper[index] = stable[index] + determinant[index];
    }
    lower[0] = stable[0];
    upper[0] = stable[0];

    BigInt interval_lower = chain_diff_interval_lower(lower, upper, chain_m);
    std::cout << "D_" << rank << "=SO(" << (2 * rank) << ") interval" << std::endl;
    std::cout << "  moment interval: det_s <= M_s <= stable_s + det_s" << std::endl;
    std::cout << "  D_" << rank << " interval lower Chain diff D(" << chain_m
              << ") = " << interval_lower << std::endl;
    if (interval_lower <= 0) throw std::runtime_error("D interval lower Chain diff is not positive");
}

void run_m22_remaining_lower_certificate(
    const std::string& stable_path,
    bool progress,
    bool run_bc_rows,
    bool run_d_rows,
    int bc_low_rank = 2,
    int bc_high_rank = 5,
    bool run_b_family = true,
    bool run_c_family = true,
    int bc_exact_through_override = -1,
    bool include_positive_quadratic_exact = false,
    int chain_m = 22,
    int d_low_rank = 4,
    int d_high_rank = 11
) {
    int moment_max = 2 * chain_m + 5;
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_max);
    BigInt stable_diff = chain_diff(stable, chain_m);
    std::vector<BigInt> linear_coeffs = chain_diff_linear_coefficients(stable, chain_m);
    auto quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m);

    std::vector<std::pair<int, BigInt>> negative_linear;
    for (int index = 0; index <= moment_max; ++index) {
        if (linear_coeffs[index] < 0) negative_linear.emplace_back(index, linear_coeffs[index]);
    }
    std::vector<std::tuple<int, int, BigInt>> negative_quadratic;
    for (const auto& [indices, value] : quadratic_coeffs) {
        if (value < 0) negative_quadratic.emplace_back(indices.first, indices.second, value);
    }

    std::cout << "m=" << chain_m << " remaining low-row lower certificate" << std::endl;
    std::cout << "stable Chain diff D(" << chain_m << ") = " << stable_diff << std::endl;
    std::cout << "negative linear coefficient indices = [";
    for (std::size_t i = 0; i < negative_linear.size(); ++i) {
        if (i) std::cout << ", ";
        std::cout << negative_linear[i].first;
    }
    std::cout << "]" << std::endl;
    std::cout << "negative quadratic pairs = [";
    for (std::size_t i = 0; i < negative_quadratic.size(); ++i) {
        if (i) std::cout << ", ";
        const auto& [first, second, value] = negative_quadratic[i];
        std::cout << "(" << first << "," << second << "," << value << ")";
    }
    std::cout << "]" << std::endl;
    if (stable_diff <= 0) throw std::runtime_error("stable Chain diff is not positive");

    BigInt minimum_margin;
    std::string minimum_label;
    bool have_margin = false;
    auto record_margin = [&](const std::string& label, const BigInt& margin) {
        if (!have_margin || margin < minimum_margin) {
            minimum_margin = margin;
            minimum_label = label;
            have_margin = true;
        }
    };

    auto run_bc_row = [&](char family, int rank) {
        std::string label = std::string(1, family) + "_" + std::to_string(rank);
        int first_boundary = 2 * rank + 2;
        std::vector<std::tuple<int, int, BigInt>> active_negative_quadratic;
        int exact_through = 31;
        for (const auto& [first, second, value] : negative_quadratic) {
            if (first >= first_boundary && second >= first_boundary) {
                active_negative_quadratic.emplace_back(first, second, value);
                exact_through = std::max(exact_through, std::max(first, second));
            }
        }
        if (rank <= 4) exact_through = std::max(exact_through, 41);
        if (bc_exact_through_override >= 0) {
            if (bc_exact_through_override < first_boundary) {
                throw std::runtime_error("exact-through override is below first boundary");
            }
            if (bc_exact_through_override > moment_max) {
                throw std::runtime_error("exact-through override exceeds moment range");
            }
            exact_through = bc_exact_through_override;
        }
        for (const auto& [first, second, value] : active_negative_quadratic) {
            if (first > exact_through || second > exact_through) {
                throw std::runtime_error(label + " active negative quadratic term outside exact window");
            }
        }

        std::cout << "begin " << label
                  << ": first boundary = " << first_boundary
                  << "; exact through m_" << exact_through << std::endl;
        std::map<int, BigInt> deltas;
        if (use_reverse_target_sums) {
            for (int moment_index = first_boundary; moment_index <= exact_through; ++moment_index) {
                std::vector<std::vector<Partition>> target_sets;
                target_sets.push_back(bc_boundary_targets(family, rank, moment_index));
                std::vector<BigInt> totals = pieri_e2_coefficient_sums_reverse(
                    moment_index,
                    target_sets,
                    progress,
                    label + " Delta_" + std::to_string(moment_index)
                );
                BigInt correction = -totals[0];
                if (correction > 0) throw std::runtime_error(label + " boundary correction is positive");
                deltas[moment_index] = correction;
                std::cout << "  " << label << " Delta_" << moment_index
                          << " = " << correction << std::endl;
            }
        } else if (use_forward_coefficients) {
            ankerl::unordered_dense::map<Partition, BigInt, PartitionHash> current;
            current.emplace(Partition{}, BigInt(1));
            auto start = std::chrono::steady_clock::now();
            for (int moment_index = 1; moment_index <= exact_through; ++moment_index) {
                current = pieri_e2_forward_step(
                    current, progress, label, moment_index, exact_through, start);
                if (progress) {
                    auto now = std::chrono::steady_clock::now();
                    double seconds = std::chrono::duration<double>(now - start).count();
                    std::cout << "    " << label << " forward step "
                              << moment_index << "/" << exact_through
                              << " states=" << current.size()
                              << " elapsed=" << seconds << "s" << std::endl;
                }
                if (moment_index < first_boundary) continue;
                BigInt total = 0;
                std::vector<Partition> targets =
                    bc_boundary_targets(family, rank, moment_index);
                auto target_start = std::chrono::steady_clock::now();
                for (std::size_t index = 0; index < targets.size(); ++index) {
                    auto found = current.find(targets[index]);
                    if (found != current.end()) total += found->second;
                    if (progress && should_print_progress(index, targets.size())) {
                        auto now = std::chrono::steady_clock::now();
                        double seconds = std::chrono::duration<double>(now - target_start).count();
                        std::cout << "    m_" << moment_index << ": "
                                  << family << " incremental target "
                                  << (index + 1) << "/" << targets.size()
                                  << " elapsed=" << seconds << "s" << std::endl;
                    }
                }
                BigInt correction = -total;
                if (correction > 0) throw std::runtime_error(label + " boundary correction is positive");
                deltas[moment_index] = correction;
                std::cout << "  " << label << " Delta_" << moment_index
                          << " = " << correction << std::endl;
            }
        } else {
            for (int moment_index = first_boundary; moment_index <= exact_through; ++moment_index) {
                BigInt correction =
                    family == 'B'
                        ? b_boundary_correction(rank, moment_index, progress)
                        : c_boundary_correction(rank, moment_index, progress);
                clear_predecessor_cache_after_moment();
                if (correction > 0) throw std::runtime_error(label + " boundary correction is positive");
                deltas[moment_index] = correction;
                std::cout << "  " << label << " Delta_" << moment_index
                          << " = " << correction << std::endl;
            }
        }

        int tail_start = std::max(exact_through + 1, first_boundary);
        BigInt negative_linear_bound = 0;
        for (int index = tail_start; index <= moment_max; ++index) {
            if (linear_coeffs[index] > 0) negative_linear_bound -= linear_coeffs[index] * stable[index];
        }
        BigInt linear_exact = 0;
        for (const auto& [index, delta] : deltas) {
            linear_exact += linear_coeffs[index] * delta;
        }
        BigInt quadratic_exact = 0;
        for (const auto& [first, second, value] : active_negative_quadratic) {
            quadratic_exact += value * deltas[first] * deltas[second];
        }
        BigInt positive_quadratic_exact = 0;
        if (include_positive_quadratic_exact) {
            for (const auto& [indices, value] : quadratic_coeffs) {
                if (value > 0
                    && deltas.find(indices.first) != deltas.end()
                    && deltas.find(indices.second) != deltas.end()) {
                    positive_quadratic_exact +=
                        value * deltas[indices.first] * deltas[indices.second];
                }
            }
        }
        BigInt lower_margin =
            stable_diff
            + negative_linear_bound
            + linear_exact
            + quadratic_exact
            + positive_quadratic_exact;
        std::cout << label
                  << " tail = " << negative_linear_bound
                  << "; linear exact = " << linear_exact
                  << "; negative quadratic exact = " << quadratic_exact;
        if (include_positive_quadratic_exact) {
            std::cout << "; positive quadratic exact = " << positive_quadratic_exact;
        }
        std::cout
                  << "; lower margin = " << lower_margin << std::endl;
        if (lower_margin <= 0) throw std::runtime_error(label + " lower margin is not positive");
        record_margin(label, lower_margin);
    };

    if (run_bc_rows) {
        std::vector<char> families;
        if (run_b_family) families.push_back('B');
        if (run_c_family) families.push_back('C');
        for (char family : families) {
            for (int rank = bc_high_rank; rank >= bc_low_rank; --rank) {
                run_bc_row(family, rank);
            }
        }
    }

    if (!run_d_rows) {
        std::cout << "minimum lower margin = " << minimum_label
                  << ": " << minimum_margin << std::endl;
        return;
    }

    for (int rank = d_high_rank; rank >= d_low_rank; --rank) {
        std::string label = "D_" + std::to_string(rank);
        std::cout << "begin " << label << " lower row" << std::endl;
        std::map<int, BigInt> deltas;
        std::map<int, bool> active_indices;
        for (const auto& [index, value] : negative_linear) {
            if (index >= rank) active_indices[index] = true;
        }
        for (const auto& [first, second, value] : negative_quadratic) {
            if (first >= rank && second >= rank) {
                active_indices[first] = true;
                active_indices[second] = true;
            }
        }
        for (const auto& [moment_index, _unused] : active_indices) {
            BigInt correction = d_determinant_correction(rank, moment_index, progress);
            clear_predecessor_cache_after_moment();
            if (correction < 0) throw std::runtime_error(label + " determinant correction is negative");
            deltas[moment_index] = correction;
            std::cout << "  " << label << " Delta_" << moment_index
                      << " = " << correction << std::endl;
        }

        BigInt negative_linear_sum = 0;
        for (const auto& [index, value] : negative_linear) {
            if (index >= rank) negative_linear_sum += value * deltas[index];
        }
        BigInt negative_quadratic_sum = 0;
        for (const auto& [first, second, value] : negative_quadratic) {
            if (first >= rank && second >= rank) {
                negative_quadratic_sum += value * deltas[first] * deltas[second];
            }
        }
        BigInt lower_margin = stable_diff + negative_linear_sum + negative_quadratic_sum;
        std::cout << label
                  << " negative linear = " << negative_linear_sum
                  << "; negative quadratic = " << negative_quadratic_sum
                  << "; lower margin = " << lower_margin << std::endl;
        if (lower_margin <= 0) throw std::runtime_error(label + " lower margin is not positive");
        record_margin(label, lower_margin);
    }

    std::cout << "minimum lower margin = " << minimum_label
              << ": " << minimum_margin << std::endl;
}

void run_b11_post_m29_bridge_range_certificate(
    const std::string& stable_path,
    bool progress,
    bool include_positive_quadratic_exact
) {
    constexpr int rank = 11;
    constexpr int first_boundary = 2 * rank + 2;
    constexpr int exact_through = 51;
    constexpr int low_chain_m = 32;
    constexpr int high_chain_m = 36;
    constexpr int max_moment = 2 * high_chain_m + 5;

    std::vector<BigInt> stable = read_stable_moments(stable_path, max_moment);
    std::map<int, BigInt> deltas;

    std::cout << "B_11 post-m=29 bridge range certificate" << std::endl;
    std::cout << "chain m range = " << low_chain_m << ".." << high_chain_m
              << "; exact deltas m_" << first_boundary << "..m_"
              << exact_through << std::endl;

    for (int moment_index = first_boundary; moment_index <= exact_through; ++moment_index) {
        std::vector<std::vector<Partition>> target_sets;
        target_sets.push_back(bc_boundary_targets('B', rank, moment_index));
        std::vector<BigInt> totals = pieri_e2_coefficient_sums_reverse(
            moment_index,
            target_sets,
            progress,
            "B_11 bridge Delta_" + std::to_string(moment_index)
        );
        BigInt correction = -totals[0];
        if (correction > 0) throw std::runtime_error("B_11 boundary correction is positive");
        deltas[moment_index] = correction;
        std::cout << "  B_11 Delta_" << moment_index
                  << " = " << correction << std::endl;
    }

    BigInt minimum_margin;
    int minimum_chain_m = 0;
    bool have_margin = false;

    for (int chain_m = low_chain_m; chain_m <= high_chain_m; ++chain_m) {
        int moment_max = 2 * chain_m + 5;
        BigInt stable_diff = chain_diff(stable, chain_m);
        std::vector<BigInt> linear_coeffs = chain_diff_linear_coefficients(stable, chain_m);
        auto quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m);
        if (stable_diff <= 0) throw std::runtime_error("stable Chain diff is not positive");

        std::vector<std::tuple<int, int, BigInt>> active_negative_quadratic;
        for (const auto& [indices, value] : quadratic_coeffs) {
            if (value < 0 && indices.first >= first_boundary && indices.second >= first_boundary) {
                if (indices.first > exact_through || indices.second > exact_through) {
                    throw std::runtime_error("B_11 active negative quadratic term outside exact window");
                }
                active_negative_quadratic.emplace_back(indices.first, indices.second, value);
            }
        }

        BigInt negative_linear_bound = 0;
        for (int index = exact_through + 1; index <= moment_max; ++index) {
            if (linear_coeffs[index] > 0) negative_linear_bound -= linear_coeffs[index] * stable[index];
        }

        BigInt linear_exact = 0;
        for (const auto& [index, delta] : deltas) {
            if (index <= moment_max) linear_exact += linear_coeffs[index] * delta;
        }

        BigInt negative_quadratic_exact = 0;
        for (const auto& [first, second, value] : active_negative_quadratic) {
            negative_quadratic_exact += value * deltas.at(first) * deltas.at(second);
        }

        BigInt positive_quadratic_exact = 0;
        if (include_positive_quadratic_exact) {
            for (const auto& [indices, value] : quadratic_coeffs) {
                if (value > 0
                    && deltas.find(indices.first) != deltas.end()
                    && deltas.find(indices.second) != deltas.end()) {
                    positive_quadratic_exact +=
                        value * deltas[indices.first] * deltas[indices.second];
                }
            }
        }

        BigInt lower_margin =
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + positive_quadratic_exact;

        int target_n = 2 * chain_m + 3;
        std::cout << "B_11 chain m=" << chain_m
                  << " target n=" << target_n
                  << " stable diff = " << stable_diff
                  << "; tail = " << negative_linear_bound
                  << "; linear exact = " << linear_exact
                  << "; negative quadratic exact = " << negative_quadratic_exact;
        if (include_positive_quadratic_exact) {
            std::cout << "; positive quadratic exact = " << positive_quadratic_exact;
        }
        std::cout << "; lower margin = " << lower_margin << std::endl;

        if (lower_margin <= 0) throw std::runtime_error("B_11 bridge lower margin is not positive");
        if (!have_margin || lower_margin < minimum_margin) {
            minimum_margin = lower_margin;
            minimum_chain_m = chain_m;
            have_margin = true;
        }
    }

    std::cout << "minimum lower margin = B_11 m=" << minimum_chain_m
              << ": " << minimum_margin << std::endl;
}

void run_m22_bc_rank_pair_positive_certificate(
    const std::string& stable_path,
    bool progress,
    int rank,
    int exact_through_override = -1,
    int chain_m = 22
) {
    int moment_max = 2 * chain_m + 5;
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_max);
    BigInt stable_diff = chain_diff(stable, chain_m);
    std::vector<BigInt> linear_coeffs = chain_diff_linear_coefficients(stable, chain_m);
    auto quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m);

    std::vector<std::tuple<int, int, BigInt>> negative_quadratic;
    for (const auto& [indices, value] : quadratic_coeffs) {
        if (value < 0) negative_quadratic.emplace_back(indices.first, indices.second, value);
    }

    struct RowState {
        char family;
        int rank;
        int first_boundary;
        int exact_through;
        std::string label;
        std::map<int, BigInt> deltas;
    };

    std::vector<RowState> rows;
    for (char family : std::vector<char>{'B', 'C'}) {
        RowState row{family, rank, 2 * rank + 2, 31, "", {}};
        row.label = std::string(1, family) + "_" + std::to_string(rank);
        for (const auto& [first, second, value] : negative_quadratic) {
            if (first >= row.first_boundary && second >= row.first_boundary) {
                row.exact_through = std::max(row.exact_through, std::max(first, second));
            }
        }
        if (exact_through_override >= 0) {
            if (exact_through_override < row.first_boundary) {
                throw std::runtime_error("exact-through override is below first boundary");
            }
            if (exact_through_override > moment_max) {
                throw std::runtime_error("exact-through override exceeds moment range");
            }
            row.exact_through = exact_through_override;
        }
        rows.push_back(std::move(row));
    }

    int max_exact_through = 0;
    for (const RowState& row : rows) {
        max_exact_through = std::max(max_exact_through, row.exact_through);
    }

    std::cout << "m=" << chain_m << " B/C rank-pair positive-quadratic certificate" << std::endl;
    std::cout << "rank = " << rank << std::endl;
    std::cout << "stable Chain diff D(" << chain_m << ") = " << stable_diff << std::endl;
    for (const RowState& row : rows) {
        std::cout << "begin " << row.label
                  << ": first boundary = " << row.first_boundary
                  << "; exact through m_" << row.exact_through << std::endl;
    }
    if (stable_diff <= 0) throw std::runtime_error("stable Chain diff is not positive");

    ankerl::unordered_dense::map<Partition, BigInt, PartitionHash> current;
    current.emplace(Partition{}, BigInt(1));
    auto start = std::chrono::steady_clock::now();
    for (int moment_index = 1; moment_index <= max_exact_through; ++moment_index) {
        current = pieri_e2_forward_step(
            current, progress, "shared", moment_index, max_exact_through, start);
        if (progress) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    shared forward step "
                      << moment_index << "/" << max_exact_through
                      << " states=" << current.size()
                      << " elapsed=" << seconds << "s" << std::endl;
        }
        for (RowState& row : rows) {
            if (moment_index < row.first_boundary || moment_index > row.exact_through) continue;
            BigInt total = 0;
            std::vector<Partition> targets =
                bc_boundary_targets(row.family, row.rank, moment_index);
            auto target_start = std::chrono::steady_clock::now();
            for (std::size_t index = 0; index < targets.size(); ++index) {
                auto found = current.find(targets[index]);
                if (found != current.end()) total += found->second;
                if (progress && should_print_progress(index, targets.size())) {
                    auto now = std::chrono::steady_clock::now();
                    double seconds = std::chrono::duration<double>(now - target_start).count();
                    std::cout << "    m_" << moment_index << ": "
                              << row.label << " shared target "
                              << (index + 1) << "/" << targets.size()
                              << " elapsed=" << seconds << "s" << std::endl;
                }
            }
            BigInt correction = -total;
            if (correction > 0) throw std::runtime_error(row.label + " boundary correction is positive");
            row.deltas[moment_index] = correction;
            std::cout << "  " << row.label << " Delta_" << moment_index
                      << " = " << correction << std::endl;
        }
    }

    BigInt minimum_margin;
    std::string minimum_label;
    bool have_margin = false;
    auto record_margin = [&](const std::string& label, const BigInt& margin) {
        if (!have_margin || margin < minimum_margin) {
            minimum_margin = margin;
            minimum_label = label;
            have_margin = true;
        }
    };

    for (const RowState& row : rows) {
        for (const auto& [first, second, value] : negative_quadratic) {
            if (first >= row.first_boundary && second >= row.first_boundary
                && (first > row.exact_through || second > row.exact_through)) {
                throw std::runtime_error(row.label + " active negative quadratic term outside exact window");
            }
        }

        int tail_start = row.exact_through + 1;
        BigInt negative_linear_bound = 0;
        for (int index = tail_start; index <= moment_max; ++index) {
            if (linear_coeffs[index] > 0) negative_linear_bound -= linear_coeffs[index] * stable[index];
        }
        BigInt linear_exact = 0;
        for (const auto& [index, delta] : row.deltas) {
            linear_exact += linear_coeffs[index] * delta;
        }
        BigInt negative_quadratic_exact = 0;
        for (const auto& [first, second, value] : negative_quadratic) {
            if (first >= row.first_boundary && second >= row.first_boundary) {
                negative_quadratic_exact += value * row.deltas.at(first) * row.deltas.at(second);
            }
        }
        BigInt positive_quadratic_exact = 0;
        for (const auto& [indices, value] : quadratic_coeffs) {
            if (value > 0
                && row.deltas.find(indices.first) != row.deltas.end()
                && row.deltas.find(indices.second) != row.deltas.end()) {
                positive_quadratic_exact +=
                    value * row.deltas.at(indices.first) * row.deltas.at(indices.second);
            }
        }
        BigInt lower_margin =
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + positive_quadratic_exact;
        std::cout << row.label
                  << " tail = " << negative_linear_bound
                  << "; linear exact = " << linear_exact
                  << "; negative quadratic exact = " << negative_quadratic_exact
                  << "; positive quadratic exact = " << positive_quadratic_exact
                  << "; lower margin = " << lower_margin << std::endl;
        if (lower_margin <= 0) throw std::runtime_error(row.label + " lower margin is not positive");
        record_margin(row.label, lower_margin);
    }

    std::cout << "minimum lower margin = " << minimum_label
              << ": " << minimum_margin << std::endl;
}

void run_bc_rank_pair_delta_certificate(
    bool progress,
    int rank,
    int moment_index
) {
    if (rank < 2) {
        throw std::runtime_error("B/C rank-pair delta rank must be at least 2");
    }
    int first_boundary = 2 * rank + 2;
    if (moment_index < first_boundary) {
        throw std::runtime_error("B/C rank-pair delta moment is below first boundary");
    }

    std::vector<std::vector<Partition>> target_sets;
    std::vector<std::string> labels;
    for (char family : std::vector<char>{'B', 'C'}) {
        target_sets.push_back(bc_boundary_targets(family, rank, moment_index));
        labels.push_back(std::string(1, family) + "_" + std::to_string(rank));
    }

    std::cout << "B/C rank-pair single-delta reverse certificate" << std::endl;
    std::cout << "rank = " << rank << std::endl;
    std::cout << "moment = " << moment_index << std::endl;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        std::cout << labels[index] << " target count = "
                  << target_sets[index].size() << std::endl;
    }

    std::vector<BigInt> totals = pieri_e2_coefficient_sums_reverse(
        moment_index,
        target_sets,
        progress,
        "B/C_" + std::to_string(rank) + " Delta_" + std::to_string(moment_index)
    );
    for (std::size_t index = 0; index < labels.size(); ++index) {
        BigInt correction = -totals[index];
        if (correction > 0) {
            throw std::runtime_error(labels[index] + " boundary correction is positive");
        }
        std::cout << "  " << labels[index] << " Delta_" << moment_index
                  << " = " << correction << std::endl;
    }
}

void run_bc_family_delta_certificate(
    bool progress,
    char family,
    int rank,
    int moment_index
) {
    if (family != 'B' && family != 'C') {
        throw std::runtime_error("B/C family delta family must be B or C");
    }
    if (rank < 2) {
        throw std::runtime_error("B/C family delta rank must be at least 2");
    }
    int first_boundary = 2 * rank + 2;
    if (moment_index < first_boundary) {
        throw std::runtime_error("B/C family delta moment is below first boundary");
    }

    std::vector<std::vector<Partition>> target_sets;
    target_sets.push_back(bc_boundary_targets(family, rank, moment_index));
    std::string label = std::string(1, family) + "_" + std::to_string(rank);

    std::cout << "B/C family single-delta reverse certificate" << std::endl;
    std::cout << "family = " << family << std::endl;
    std::cout << "rank = " << rank << std::endl;
    std::cout << "moment = " << moment_index << std::endl;
    std::cout << label << " target count = "
              << target_sets[0].size() << std::endl;

    std::vector<BigInt> totals = pieri_e2_coefficient_sums_reverse(
        moment_index,
        target_sets,
        progress,
        label + " Delta_" + std::to_string(moment_index)
    );
    BigInt correction = -totals[0];
    if (correction > 0) {
        throw std::runtime_error(label + " boundary correction is positive");
    }
    std::cout << "  " << label << " Delta_" << moment_index
              << " = " << correction << std::endl;
}

void run_bc_rank_pair_range_positive_certificate(
    const std::string& stable_path,
    bool progress,
    int low_rank,
    int high_rank,
    int exact_through_override = -1,
    int chain_m = 22
) {
    if (low_rank > high_rank) throw std::runtime_error("empty B/C rank-pair range");
    int moment_max = 2 * chain_m + 5;
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_max);
    BigInt stable_diff = chain_diff(stable, chain_m);
    std::vector<BigInt> linear_coeffs = chain_diff_linear_coefficients(stable, chain_m);
    auto quadratic_coeffs = chain_diff_quadratic_coefficients(chain_m);

    std::vector<std::tuple<int, int, BigInt>> negative_quadratic;
    for (const auto& [indices, value] : quadratic_coeffs) {
        if (value < 0) negative_quadratic.emplace_back(indices.first, indices.second, value);
    }

    struct RowState {
        char family;
        int rank;
        int first_boundary;
        int exact_through;
        std::string label;
        std::map<int, BigInt> deltas;
    };

    std::vector<RowState> rows;
    for (int rank = high_rank; rank >= low_rank; --rank) {
        for (char family : std::vector<char>{'B', 'C'}) {
            RowState row{family, rank, 2 * rank + 2, 31, "", {}};
            row.label = std::string(1, family) + "_" + std::to_string(rank);
            for (const auto& [first, second, value] : negative_quadratic) {
                if (first >= row.first_boundary && second >= row.first_boundary) {
                    row.exact_through = std::max(row.exact_through, std::max(first, second));
                }
            }
            if (rank <= 4) row.exact_through = std::max(row.exact_through, 41);
            if (exact_through_override >= 0) {
                if (exact_through_override < row.first_boundary) {
                    throw std::runtime_error("exact-through override is below first boundary");
                }
                if (exact_through_override > moment_max) {
                    throw std::runtime_error("exact-through override exceeds moment range");
                }
                row.exact_through = exact_through_override;
            }
            rows.push_back(std::move(row));
        }
    }

    int max_exact_through = 0;
    for (const RowState& row : rows) {
        max_exact_through = std::max(max_exact_through, row.exact_through);
    }

    std::cout << "m=" << chain_m << " B/C rank-pair range positive-quadratic certificate" << std::endl;
    std::cout << "ranks = " << low_rank << ".." << high_rank << std::endl;
    std::cout << "stable Chain diff D(" << chain_m << ") = " << stable_diff << std::endl;
    for (const RowState& row : rows) {
        std::cout << "begin " << row.label
                  << ": first boundary = " << row.first_boundary
                  << "; exact through m_" << row.exact_through << std::endl;
    }
    if (stable_diff <= 0) throw std::runtime_error("stable Chain diff is not positive");

    if (use_reverse_target_sums) {
        for (int moment_index = 1; moment_index <= max_exact_through; ++moment_index) {
            std::vector<std::vector<Partition>> target_sets;
            std::vector<std::size_t> row_indices;
            for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
                RowState& row = rows[row_index];
                if (moment_index < row.first_boundary || moment_index > row.exact_through) continue;
                target_sets.push_back(bc_boundary_targets(row.family, row.rank, moment_index));
                row_indices.push_back(row_index);
            }
            if (target_sets.empty()) continue;
            std::vector<BigInt> totals = pieri_e2_coefficient_sums_reverse(
                moment_index,
                target_sets,
                progress,
                "B/C range Delta_" + std::to_string(moment_index)
            );
            for (std::size_t index = 0; index < row_indices.size(); ++index) {
                RowState& row = rows[row_indices[index]];
                BigInt correction = -totals[index];
                if (correction > 0) throw std::runtime_error(row.label + " boundary correction is positive");
                row.deltas[moment_index] = correction;
                std::cout << "  " << row.label << " Delta_" << moment_index
                          << " = " << correction << std::endl;
            }
        }
    } else {
        ankerl::unordered_dense::map<Partition, BigInt, PartitionHash> current;
        current.emplace(Partition{}, BigInt(1));
        auto start = std::chrono::steady_clock::now();
        for (int moment_index = 1; moment_index <= max_exact_through; ++moment_index) {
            current = pieri_e2_forward_step(
                current, progress, "shared-range", moment_index, max_exact_through, start);
            if (progress) {
                auto now = std::chrono::steady_clock::now();
                double seconds = std::chrono::duration<double>(now - start).count();
                std::cout << "    shared-range forward step "
                          << moment_index << "/" << max_exact_through
                          << " states=" << current.size()
                          << " elapsed=" << seconds << "s" << std::endl;
            }
            for (RowState& row : rows) {
                if (moment_index < row.first_boundary || moment_index > row.exact_through) continue;
                BigInt total = 0;
                std::vector<Partition> targets =
                    bc_boundary_targets(row.family, row.rank, moment_index);
                auto target_start = std::chrono::steady_clock::now();
                for (std::size_t index = 0; index < targets.size(); ++index) {
                    auto found = current.find(targets[index]);
                    if (found != current.end()) total += found->second;
                    if (progress && should_print_progress(index, targets.size())) {
                        auto now = std::chrono::steady_clock::now();
                        double seconds = std::chrono::duration<double>(now - target_start).count();
                        std::cout << "    m_" << moment_index << ": "
                                  << row.label << " shared-range target "
                                  << (index + 1) << "/" << targets.size()
                                  << " elapsed=" << seconds << "s" << std::endl;
                    }
                }
                BigInt correction = -total;
                if (correction > 0) throw std::runtime_error(row.label + " boundary correction is positive");
                row.deltas[moment_index] = correction;
                std::cout << "  " << row.label << " Delta_" << moment_index
                          << " = " << correction << std::endl;
            }
        }
    }

    BigInt minimum_margin;
    std::string minimum_label;
    bool have_margin = false;
    auto record_margin = [&](const std::string& label, const BigInt& margin) {
        if (!have_margin || margin < minimum_margin) {
            minimum_margin = margin;
            minimum_label = label;
            have_margin = true;
        }
    };

    for (const RowState& row : rows) {
        for (const auto& [first, second, value] : negative_quadratic) {
            if (first >= row.first_boundary && second >= row.first_boundary
                && (first > row.exact_through || second > row.exact_through)) {
                throw std::runtime_error(row.label + " active negative quadratic term outside exact window");
            }
        }

        int tail_start = row.exact_through + 1;
        BigInt negative_linear_bound = 0;
        for (int index = tail_start; index <= moment_max; ++index) {
            if (linear_coeffs[index] > 0) negative_linear_bound -= linear_coeffs[index] * stable[index];
        }
        BigInt linear_exact = 0;
        for (const auto& [index, delta] : row.deltas) {
            linear_exact += linear_coeffs[index] * delta;
        }
        BigInt negative_quadratic_exact = 0;
        for (const auto& [first, second, value] : negative_quadratic) {
            if (first >= row.first_boundary && second >= row.first_boundary) {
                negative_quadratic_exact += value * row.deltas.at(first) * row.deltas.at(second);
            }
        }
        BigInt positive_quadratic_exact = 0;
        for (const auto& [indices, value] : quadratic_coeffs) {
            if (value > 0
                && row.deltas.find(indices.first) != row.deltas.end()
                && row.deltas.find(indices.second) != row.deltas.end()) {
                positive_quadratic_exact +=
                    value * row.deltas.at(indices.first) * row.deltas.at(indices.second);
            }
        }
        BigInt lower_margin =
            stable_diff
            + negative_linear_bound
            + linear_exact
            + negative_quadratic_exact
            + positive_quadratic_exact;
        std::cout << row.label
                  << " tail = " << negative_linear_bound
                  << "; linear exact = " << linear_exact
                  << "; negative quadratic exact = " << negative_quadratic_exact
                  << "; positive quadratic exact = " << positive_quadratic_exact
                  << "; lower margin = " << lower_margin << std::endl;
        if (lower_margin <= 0) throw std::runtime_error(row.label + " lower margin is not positive");
        record_margin(row.label, lower_margin);
    }

    std::cout << "minimum lower margin = " << minimum_label
              << ": " << minimum_margin << std::endl;
}

using Exp2 = std::pair<int, int>;
using Laurent2 = std::map<Exp2, BigInt>;

void add_laurent2(Laurent2& out, int first, int second, const BigInt& value) {
    if (value == 0) return;
    Exp2 key{first, second};
    BigInt& slot = out[key];
    slot += value;
    if (slot == 0) out.erase(key);
}

Laurent2 multiply_laurent2(const Laurent2& left, const Laurent2& right) {
    Laurent2 out;
    for (const auto& [left_exp, left_coeff] : left) {
        for (const auto& [right_exp, right_coeff] : right) {
            add_laurent2(
                out,
                left_exp.first + right_exp.first,
                left_exp.second + right_exp.second,
                left_coeff * right_coeff
            );
        }
    }
    return out;
}

Laurent2 b2_weyl_density() {
    Laurent2 density;
    density[{0, 0}] = 1;
    for (const Exp2& root : std::vector<Exp2>{{1, 0}, {0, 1}, {1, 1}, {1, -1}}) {
        Laurent2 factor;
        factor[{0, 0}] = 2;
        factor[root] = -1;
        factor[{-root.first, -root.second}] = -1;
        density = multiply_laurent2(density, factor);
    }
    return density;
}

Laurent2 multiply_by_b2_adjoint_character(const Laurent2& poly) {
    static const std::vector<std::tuple<int, int, int>> character_terms{
        {0, 0, 2},
        {1, 0, 1}, {-1, 0, 1}, {0, 1, 1}, {0, -1, 1},
        {1, 1, 1}, {-1, -1, 1}, {1, -1, 1}, {-1, 1, 1},
    };
    Laurent2 out;
    for (const auto& [exp, coeff] : poly) {
        for (const auto& [first, second, mult] : character_terms) {
            add_laurent2(out, exp.first + first, exp.second + second, coeff * mult);
        }
    }
    return out;
}

BigInt constant_term_against_density(const Laurent2& poly, const Laurent2& density) {
    BigInt total = 0;
    for (const auto& [exp, coeff] : poly) {
        auto found = density.find({-exp.first, -exp.second});
        if (found != density.end()) total += coeff * found->second;
    }
    return total;
}

std::vector<BigInt> b2_weyl_moments(int moment_max, bool progress) {
    Laurent2 density = b2_weyl_density();
    BigInt density_ct = constant_term_against_density(Laurent2{{{0, 0}, 1}}, density);
    if (density_ct != 8) throw std::runtime_error("B2 Weyl density constant term is not |W|=8");

    std::vector<BigInt> moments(moment_max + 1);
    Laurent2 power;
    power[{0, 0}] = 1;
    auto start = std::chrono::steady_clock::now();
    for (int index = 0; index <= moment_max; ++index) {
        BigInt raw = constant_term_against_density(power, density);
        if (raw % 8 != 0) throw std::runtime_error("B2 Weyl moment numerator is not divisible by 8");
        moments[index] = raw / 8;
        if (progress) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    B2 Weyl moment " << index << "/" << moment_max
                      << " support=" << power.size()
                      << " elapsed=" << seconds << "s" << std::endl;
        }
        if (index < moment_max) power = multiply_by_b2_adjoint_character(power);
    }
    return moments;
}

void run_b2_c2_weyl_certificate(const std::string& stable_path, bool progress, int chain_m) {
    int moment_max = 2 * chain_m + 5;
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_max);
    std::vector<BigInt> moments = b2_weyl_moments(moment_max, progress);
    BigInt row_diff = chain_diff(moments, chain_m);

    std::cout << "B2/C2 exact Weyl m=" << chain_m << " certificate" << std::endl;
    std::cout << "Weyl density constant term = 8" << std::endl;
    std::cout << "moment range m_0..m_" << moment_max << std::endl;
    int first_report = std::max(0, moment_max - 12);
    for (int index = first_report; index <= moment_max; ++index) {
        std::cout << "m_" << index << " = " << moments[index]
                  << "; Delta_" << index << " = " << (moments[index] - stable[index])
                  << std::endl;
    }
    std::cout << "B_2/C_2 Chain diff D(" << chain_m << ") = " << row_diff << std::endl;
    if (row_diff <= 0) throw std::runtime_error("B2/C2 Chain diff is not positive");
}

struct Exp3 {
    int first;
    int second;
    int third;

    bool operator==(const Exp3& other) const {
        return first == other.first && second == other.second && third == other.third;
    }
};

struct Exp3Hash {
    std::size_t operator()(const Exp3& exp) const {
        std::uint64_t h = 1469598103934665603ull;
        for (int value : std::array<int, 3>{exp.first, exp.second, exp.third}) {
            std::uint64_t x = static_cast<std::uint64_t>(
                static_cast<std::int64_t>(value) + 0x100000000ll
            );
            h ^= x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
            h *= 1099511628211ull;
        }
        return static_cast<std::size_t>(h);
    }
};

using Laurent3 = ankerl::unordered_dense::map<Exp3, BigInt, Exp3Hash>;

Exp3 add_exp3(const Exp3& left, const Exp3& right) {
    return {left.first + right.first, left.second + right.second, left.third + right.third};
}

Exp3 neg_exp3(const Exp3& exp) {
    return {-exp.first, -exp.second, -exp.third};
}

void add_laurent3(Laurent3& out, const Exp3& exp, const BigInt& value) {
    if (value == 0) return;
    BigInt& slot = out[exp];
    slot += value;
    if (slot == 0) out.erase(exp);
}

std::vector<Exp3> rank3_positive_roots(char family) {
    std::vector<Exp3> roots;
    if (family == 'B') {
        roots.push_back({1, 0, 0});
        roots.push_back({0, 1, 0});
        roots.push_back({0, 0, 1});
    } else if (family == 'C') {
        roots.push_back({2, 0, 0});
        roots.push_back({0, 2, 0});
        roots.push_back({0, 0, 2});
    } else {
        throw std::runtime_error("rank-3 Weyl family must be B or C");
    }
    roots.push_back({1, 1, 0});
    roots.push_back({1, -1, 0});
    roots.push_back({1, 0, 1});
    roots.push_back({1, 0, -1});
    roots.push_back({0, 1, 1});
    roots.push_back({0, 1, -1});
    return roots;
}

Laurent3 multiply_by_weyl_factor3(const Laurent3& poly, const Exp3& root) {
    Laurent3 out;
    out.reserve(poly.size() * 3 + 16);
    for (const auto& [exp, coeff] : poly) {
        add_laurent3(out, exp, 2 * coeff);
        add_laurent3(out, add_exp3(exp, root), -coeff);
        add_laurent3(out, add_exp3(exp, neg_exp3(root)), -coeff);
    }
    return out;
}

Laurent3 rank3_weyl_density(char family) {
    Laurent3 density;
    density[{0, 0, 0}] = 1;
    for (const Exp3& root : rank3_positive_roots(family)) {
        density = multiply_by_weyl_factor3(density, root);
    }
    return density;
}

Laurent3 multiply_by_rank3_adjoint_character(const Laurent3& poly, const std::vector<Exp3>& positive_roots) {
    Laurent3 out;
    out.reserve(poly.size() * 2 + 1024);
    for (const auto& [exp, coeff] : poly) {
        add_laurent3(out, exp, 3 * coeff);
        for (const Exp3& root : positive_roots) {
            add_laurent3(out, add_exp3(exp, root), coeff);
            add_laurent3(out, add_exp3(exp, neg_exp3(root)), coeff);
        }
    }
    return out;
}

BigInt constant_term_against_density3(const Laurent3& poly, const Laurent3& density) {
    BigInt total = 0;
    for (const auto& [exp, coeff] : poly) {
        auto found = density.find(neg_exp3(exp));
        if (found != density.end()) total += coeff * found->second;
    }
    return total;
}

std::vector<BigInt> rank3_weyl_moments(char family, int moment_max, bool progress) {
    std::vector<Exp3> positive_roots = rank3_positive_roots(family);
    if (positive_roots.size() != 9) throw std::runtime_error("rank-3 positive root count is not 9");
    Laurent3 density = rank3_weyl_density(family);
    BigInt density_ct = constant_term_against_density3(Laurent3{{{0, 0, 0}, 1}}, density);
    if (density_ct != 48) throw std::runtime_error("rank-3 Weyl density constant term is not |W|=48");

    std::vector<BigInt> moments(moment_max + 1);
    Laurent3 power;
    power[{0, 0, 0}] = 1;
    auto start = std::chrono::steady_clock::now();
    for (int index = 0; index <= moment_max; ++index) {
        BigInt raw = constant_term_against_density3(power, density);
        if (raw % 48 != 0) throw std::runtime_error("rank-3 Weyl moment numerator is not divisible by 48");
        moments[index] = raw / 48;
        if (progress) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    " << family << "3 Weyl moment " << index << "/" << moment_max
                      << " support=" << power.size()
                      << " elapsed=" << seconds << "s" << std::endl;
        }
        if (index < moment_max) power = multiply_by_rank3_adjoint_character(power, positive_roots);
    }
    return moments;
}

void run_b3_c3_weyl_certificate(const std::string& stable_path, bool progress, int chain_m) {
    int moment_max = 2 * chain_m + 5;
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_max);

    for (char family : std::vector<char>{'B', 'C'}) {
        std::vector<BigInt> moments = rank3_weyl_moments(family, moment_max, progress);
        BigInt row_diff = chain_diff(moments, chain_m);

        std::cout << family << "_3 exact Weyl m=" << chain_m << " certificate" << std::endl;
        std::cout << "Weyl density constant term = 48" << std::endl;
        std::cout << "moment range m_0..m_" << moment_max << std::endl;
        int first_report = std::max(0, moment_max - 12);
        for (int index = first_report; index <= moment_max; ++index) {
            std::cout << "m_" << index << " = " << moments[index]
                      << "; Delta_" << index << " = " << (moments[index] - stable[index])
                      << std::endl;
        }
        std::cout << family << "_3 Chain diff D(" << chain_m << ") = " << row_diff << std::endl;
        if (row_diff <= 0) throw std::runtime_error(std::string(1, family) + "3 Chain diff is not positive");
    }
}

struct Exp4 {
    int first;
    int second;
    int third;
    int fourth;

    bool operator==(const Exp4& other) const {
        return first == other.first
            && second == other.second
            && third == other.third
            && fourth == other.fourth;
    }
};

struct Exp4Hash {
    std::size_t operator()(const Exp4& exp) const {
        std::uint64_t h = 1469598103934665603ull;
        for (int value : std::array<int, 4>{exp.first, exp.second, exp.third, exp.fourth}) {
            std::uint64_t x = static_cast<std::uint64_t>(
                static_cast<std::int64_t>(value) + 0x100000000ll
            );
            h ^= x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
            h *= 1099511628211ull;
        }
        return static_cast<std::size_t>(h);
    }
};

using Laurent4 = ankerl::unordered_dense::map<Exp4, BigInt, Exp4Hash>;
using Laurent4Mod = ankerl::unordered_dense::map<Exp4, std::uint64_t, Exp4Hash>;

void set_exp4_component(Exp4& exp, int index, int value) {
    if (index == 0) exp.first = value;
    else if (index == 1) exp.second = value;
    else if (index == 2) exp.third = value;
    else if (index == 3) exp.fourth = value;
    else throw std::runtime_error("Exp4 component index is out of range");
}

Exp4 add_exp4(const Exp4& left, const Exp4& right) {
    return {
        left.first + right.first,
        left.second + right.second,
        left.third + right.third,
        left.fourth + right.fourth
    };
}

Exp4 neg_exp4(const Exp4& exp) {
    return {-exp.first, -exp.second, -exp.third, -exp.fourth};
}

void add_laurent4(Laurent4& out, const Exp4& exp, const BigInt& value) {
    if (value == 0) return;
    BigInt& slot = out[exp];
    slot += value;
    if (slot == 0) out.erase(exp);
}

std::uint64_t mul_mod_u63(std::uint64_t left, std::uint64_t right, std::uint64_t modulus) {
    return static_cast<std::uint64_t>(
        (static_cast<unsigned __int128>(left) * static_cast<unsigned __int128>(right)) % modulus
    );
}

std::uint64_t pow_mod_u63(std::uint64_t base, std::uint64_t exponent, std::uint64_t modulus) {
    std::uint64_t result = 1 % modulus;
    while (exponent > 0) {
        if (exponent & 1ull) result = mul_mod_u63(result, base, modulus);
        base = mul_mod_u63(base, base, modulus);
        exponent >>= 1;
    }
    return result;
}

std::uint64_t inverse_mod_prime_u63(std::uint64_t value, std::uint64_t modulus) {
    if (value == 0 || value >= modulus) {
        throw std::runtime_error("invalid modular inverse input");
    }
    return pow_mod_u63(value, modulus - 2, modulus);
}

void add_laurent4_mod(
    Laurent4Mod& out,
    const Exp4& exp,
    std::uint64_t value,
    std::uint64_t modulus
) {
    value %= modulus;
    if (value == 0) return;
    std::uint64_t& slot = out[exp];
    slot = add_mod_u63(slot, value, modulus);
    if (slot == 0) out.erase(exp);
}

void add_laurent4_mod_scaled(
    Laurent4Mod& out,
    const Exp4& exp,
    std::uint64_t coeff,
    int scale,
    std::uint64_t modulus
) {
    if (coeff == 0 || scale == 0) return;
    std::uint64_t scaled =
        mul_mod_u63(coeff, static_cast<std::uint64_t>(std::abs(scale)), modulus);
    if (scale < 0 && scaled != 0) scaled = modulus - scaled;
    add_laurent4_mod(out, exp, scaled, modulus);
}

std::vector<Exp4> rank4_positive_roots(char family) {
    std::vector<Exp4> roots;
    if (family != 'B' && family != 'C' && family != 'D') {
        throw std::runtime_error("rank-4 Weyl family must be B, C, or D");
    }
    if (family == 'B' || family == 'C') {
        for (int index = 0; index < 4; ++index) {
            Exp4 root{0, 0, 0, 0};
            set_exp4_component(root, index, family == 'B' ? 1 : 2);
            roots.push_back(root);
        }
    }
    for (int first = 0; first < 4; ++first) {
        for (int second = first + 1; second < 4; ++second) {
            Exp4 sum_root{0, 0, 0, 0};
            set_exp4_component(sum_root, first, 1);
            set_exp4_component(sum_root, second, 1);
            roots.push_back(sum_root);

            Exp4 difference_root{0, 0, 0, 0};
            set_exp4_component(difference_root, first, 1);
            set_exp4_component(difference_root, second, -1);
            roots.push_back(difference_root);
        }
    }
    return roots;
}

int rank4_weyl_order(char family) {
    if (family == 'B' || family == 'C') return 384;
    if (family == 'D') return 192;
    throw std::runtime_error("rank-4 Weyl family must be B, C, or D");
}

Laurent4 multiply_by_weyl_factor4(const Laurent4& poly, const Exp4& root) {
    Laurent4 out;
    out.reserve(poly.size() * 3 + 16);
    for (const auto& [exp, coeff] : poly) {
        add_laurent4(out, exp, 2 * coeff);
        add_laurent4(out, add_exp4(exp, root), -coeff);
        add_laurent4(out, add_exp4(exp, neg_exp4(root)), -coeff);
    }
    return out;
}

Laurent4 rank4_weyl_density(char family) {
    Laurent4 density;
    density[{0, 0, 0, 0}] = 1;
    for (const Exp4& root : rank4_positive_roots(family)) {
        density = multiply_by_weyl_factor4(density, root);
    }
    return density;
}

Laurent4 multiply_by_rank4_adjoint_character(const Laurent4& poly, const std::vector<Exp4>& positive_roots) {
    Laurent4 out;
    out.reserve(poly.size() * 2 + 4096);
    for (const auto& [exp, coeff] : poly) {
        add_laurent4(out, exp, 4 * coeff);
        for (const Exp4& root : positive_roots) {
            add_laurent4(out, add_exp4(exp, root), coeff);
            add_laurent4(out, add_exp4(exp, neg_exp4(root)), coeff);
        }
    }
    return out;
}

Laurent4Mod multiply_by_weyl_factor4_mod(
    const Laurent4Mod& poly,
    const Exp4& root,
    std::uint64_t modulus
) {
    Laurent4Mod out;
    out.reserve(poly.size() * 3 + 16);
    for (const auto& [exp, coeff] : poly) {
        add_laurent4_mod_scaled(out, exp, coeff, 2, modulus);
        add_laurent4_mod_scaled(out, add_exp4(exp, root), coeff, -1, modulus);
        add_laurent4_mod_scaled(out, add_exp4(exp, neg_exp4(root)), coeff, -1, modulus);
    }
    return out;
}

Laurent4Mod rank4_weyl_density_mod(char family, std::uint64_t modulus) {
    Laurent4Mod density;
    density[{0, 0, 0, 0}] = 1 % modulus;
    for (const Exp4& root : rank4_positive_roots(family)) {
        density = multiply_by_weyl_factor4_mod(density, root, modulus);
    }
    return density;
}

Laurent4Mod multiply_by_rank4_adjoint_character_mod(
    const Laurent4Mod& poly,
    const std::vector<Exp4>& positive_roots,
    std::uint64_t modulus
) {
    Laurent4Mod out;
    out.reserve(poly.size() * 2 + 4096);
    for (const auto& [exp, coeff] : poly) {
        add_laurent4_mod_scaled(out, exp, coeff, 4, modulus);
        for (const Exp4& root : positive_roots) {
            add_laurent4_mod(out, add_exp4(exp, root), coeff, modulus);
            add_laurent4_mod(out, add_exp4(exp, neg_exp4(root)), coeff, modulus);
        }
    }
    return out;
}

BigInt constant_term_against_density4(const Laurent4& poly, const Laurent4& density) {
    BigInt total = 0;
    for (const auto& [exp, coeff] : poly) {
        auto found = density.find(neg_exp4(exp));
        if (found != density.end()) total += coeff * found->second;
    }
    return total;
}

std::uint64_t constant_term_against_density4_mod(
    const Laurent4Mod& poly,
    const Laurent4Mod& density,
    std::uint64_t modulus
) {
    std::uint64_t total = 0;
    for (const auto& [exp, coeff] : poly) {
        auto found = density.find(neg_exp4(exp));
        if (found != density.end()) {
            total = add_mod_u63(total, mul_mod_u63(coeff, found->second, modulus), modulus);
        }
    }
    return total;
}

std::vector<BigInt> rank4_weyl_moments(char family, int moment_max, bool progress) {
    std::vector<Exp4> positive_roots = rank4_positive_roots(family);
    std::size_t expected_root_count = family == 'D' ? 12 : 16;
    if (positive_roots.size() != expected_root_count) throw std::runtime_error("rank-4 positive root count is incorrect");
    Laurent4 density = rank4_weyl_density(family);
    BigInt density_ct = constant_term_against_density4(Laurent4{{{0, 0, 0, 0}, 1}}, density);
    int weyl_order = rank4_weyl_order(family);
    if (density_ct != weyl_order) throw std::runtime_error("rank-4 Weyl density constant term is not |W|");

    std::vector<BigInt> moments(moment_max + 1);
    Laurent4 power;
    power[{0, 0, 0, 0}] = 1;
    auto start = std::chrono::steady_clock::now();
    for (int index = 0; index <= moment_max; ++index) {
        BigInt raw = constant_term_against_density4(power, density);
        if (raw % weyl_order != 0) throw std::runtime_error("rank-4 Weyl moment numerator is not divisible by |W|");
        moments[index] = raw / weyl_order;
        if (progress) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    " << family << "4 Weyl moment " << index << "/" << moment_max
                      << " support=" << power.size()
                      << " elapsed=" << seconds << "s" << std::endl;
        }
        if (index < moment_max) power = multiply_by_rank4_adjoint_character(power, positive_roots);
    }
    return moments;
}

std::vector<std::uint64_t> rank4_weyl_moments_mod(
    char family,
    int moment_max,
    std::uint64_t modulus,
    bool progress
) {
    std::vector<Exp4> positive_roots = rank4_positive_roots(family);
    std::size_t expected_root_count = family == 'D' ? 12 : 16;
    if (positive_roots.size() != expected_root_count) {
        throw std::runtime_error("rank-4 positive root count is incorrect");
    }
    Laurent4Mod density = rank4_weyl_density_mod(family, modulus);
    Laurent4Mod one;
    one[{0, 0, 0, 0}] = 1 % modulus;
    std::uint64_t density_ct = constant_term_against_density4_mod(one, density, modulus);
    int weyl_order = rank4_weyl_order(family);
    if (density_ct != static_cast<std::uint64_t>(weyl_order)) {
        throw std::runtime_error("rank-4 modular Weyl density constant term is not |W|");
    }
    std::uint64_t weyl_order_inv =
        inverse_mod_prime_u63(static_cast<std::uint64_t>(weyl_order), modulus);

    std::vector<std::uint64_t> moments(moment_max + 1);
    Laurent4Mod power;
    power[{0, 0, 0, 0}] = 1 % modulus;
    auto start = std::chrono::steady_clock::now();
    for (int index = 0; index <= moment_max; ++index) {
        std::uint64_t raw = constant_term_against_density4_mod(power, density, modulus);
        moments[index] = mul_mod_u63(raw, weyl_order_inv, modulus);
        if (progress) {
            auto now = std::chrono::steady_clock::now();
            double seconds = std::chrono::duration<double>(now - start).count();
            std::cout << "    " << family << "4 modular Weyl " << modulus
                      << " moment " << index << "/" << moment_max
                      << " support=" << power.size()
                      << " elapsed=" << seconds << "s" << std::endl;
        }
        if (index < moment_max) {
            power = multiply_by_rank4_adjoint_character_mod(power, positive_roots, modulus);
        }
    }
    return moments;
}

void run_rank4_weyl_family_certificate(char family, const std::string& stable_path, bool progress, int chain_m) {
    int moment_max = 2 * chain_m + 5;
    std::vector<BigInt> stable = read_stable_moments(stable_path, moment_max);
    std::vector<BigInt> moments = rank4_weyl_moments(family, moment_max, progress);
    BigInt row_diff = chain_diff(moments, chain_m);

    std::cout << family << "_4 exact Weyl m=" << chain_m << " certificate" << std::endl;
    std::cout << "Weyl density constant term = " << rank4_weyl_order(family) << std::endl;
    std::cout << "moment range m_0..m_" << moment_max << std::endl;
    int first_report = std::max(0, moment_max - 12);
    for (int index = first_report; index <= moment_max; ++index) {
        std::cout << "m_" << index << " = " << moments[index]
                  << "; Delta_" << index << " = " << (moments[index] - stable[index])
                  << std::endl;
    }
    std::cout << family << "_4 Chain diff D(" << chain_m << ") = " << row_diff << std::endl;
    if (row_diff <= 0) throw std::runtime_error(std::string(1, family) + "4 Chain diff is not positive");
}

void run_rank4_weyl_mod_prime_certificate(
    char family,
    int chain_m,
    std::uint64_t modulus,
    bool progress
) {
    if (family != 'B' && family != 'C' && family != 'D') {
        throw std::runtime_error("rank-4 modular family must be B, C, or D");
    }
    if (modulus < 3 || modulus % 2 == 0) {
        throw std::runtime_error("rank-4 modular modulus must be an odd prime");
    }
    int moment_max = 2 * chain_m + 5;
    std::vector<std::uint64_t> moments =
        rank4_weyl_moments_mod(family, moment_max, modulus, progress);
    std::cout << family << "_4 MOD " << modulus
              << " Weyl density constant term = " << rank4_weyl_order(family)
              << std::endl;
    for (int index = 0; index <= moment_max; ++index) {
        std::cout << family << "_4 MOD " << modulus
                  << " m_" << index
                  << " moment=" << moments[index]
                  << std::endl;
    }
}

void run_b4_c4_weyl_certificate(const std::string& stable_path, bool progress, int chain_m) {
    for (char family : std::vector<char>{'B', 'C'}) {
        run_rank4_weyl_family_certificate(family, stable_path, progress, chain_m);
    }
}

int main(int argc, char** argv) {
    bool run_b = false;
    bool run_c = false;
    bool run_d = false;
    bool run_d_exact = false;
    bool run_d_exact_range = false;
    bool run_d_exact_moment_range = false;
    bool run_d_o_even_moment_range = false;
    bool run_d_capped_suffix = false;
    bool run_d_capped_suffix_mod_prime = false;
    bool run_d_interval = false;
    bool run_d_delta_shared_range = false;
    bool run_m22_remaining = false;
    bool run_m22_remaining_bc = false;
    bool run_m22_remaining_d = false;
    bool run_m22_b5_lower = false;
    bool run_m22_bc_single_lower = false;
    bool run_m22_bc_rank_pair_positive = false;
    bool run_b11_post_m29_bridge_range = false;
    bool run_bc_rank_pair_delta = false;
    bool run_bc_family_delta = false;
    bool run_bc_rank_pair_range_positive = false;
    bool run_b2_c2_weyl = false;
    bool run_b3_c3_weyl = false;
    bool run_b4_c4_weyl = false;
    bool run_b4_weyl = false;
    bool run_c4_weyl = false;
    bool run_d4_weyl = false;
    bool run_rank4_weyl_mod = false;
    char rank4_weyl_mod_family = 0;
    std::uint64_t rank4_weyl_mod_prime = 0;
    bool run_generic_lower = false;
    bool generic_lower_bc = false;
    bool generic_lower_d = false;
    bool generic_lower_b_family = false;
    bool generic_lower_c_family = false;
    bool progress = false;
    bool no_predecessor_cache = false;
    bool include_positive_quadratic_exact = false;
    int b_rank = 16;
    int c_rank = 16;
    int d_rank = 35;
    int d_exact_rank = 0;
    int d_exact_low_rank = 0;
    int d_exact_high_rank = 0;
    int d_exact_moment_low_rank = 0;
    int d_exact_moment_high_rank = 0;
    int d_exact_moment_low = 0;
    int d_exact_moment_high = 0;
    int d_o_even_moment_low_rank = 0;
    int d_o_even_moment_high_rank = 0;
    int d_o_even_moment_low = 0;
    int d_o_even_moment_high = 0;
    int d_capped_suffix_rank = 0;
    int d_capped_suffix_moment_low = 0;
    int d_capped_suffix_moment_high = 0;
    std::uint64_t d_capped_suffix_mod_prime = 0;
    int d_interval_rank = 0;
    int d_delta_rank = 0;
    int d_delta_moment = 0;
    int d_delta_range_rank = 0;
    int d_delta_range_low = 0;
    int d_delta_range_high = 0;
    int d_delta_shared_low_rank = 0;
    int d_delta_shared_high_rank = 0;
    int d_delta_shared_low = 0;
    int d_delta_shared_high = 0;
    int chain_m = 15;
    char m22_bc_single_family = 'B';
    int m22_bc_single_rank = 0;
    int m22_bc_exact_through_override = -1;
    int m22_bc_rank_pair = 0;
    char bc_family_delta_family = 'B';
    int bc_family_delta_rank = 0;
    int bc_family_delta_moment = 0;
    int bc_rank_pair_delta_rank = 0;
    int bc_rank_pair_delta_moment = 0;
    int bc_rank_pair_low = 0;
    int bc_rank_pair_high = 0;
    int generic_bc_low_rank = 0;
    int generic_bc_high_rank = 0;
    int generic_d_low_rank = 0;
    int generic_d_high_rank = 0;
    std::string stable_path = "ginibre_q3/references/oeis_A002137_stable.txt";
    std::vector<std::pair<char, int>> rows;

    auto add_range = [&rows](char family, int low, int high) {
        if (low > high) throw std::runtime_error("empty rank range");
        for (int rank = low; rank <= high; ++rank) rows.emplace_back(family, rank);
    };

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--b-boundary-m15-rank16-certificate") {
            run_b = true;
        } else if (arg == "--c-boundary-m15-rank16-certificate") {
            run_c = true;
        } else if (arg == "--b-boundary-m15-rank" || arg == "--b-boundary-rank") {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + arg);
            run_b = true;
            b_rank = std::stoi(argv[++i]);
        } else if (arg == "--c-boundary-m15-rank" || arg == "--c-boundary-rank") {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + arg);
            run_c = true;
            c_rank = std::stoi(argv[++i]);
        } else if (arg == "--d-boundary-rank") {
            if (i + 1 >= argc) throw std::runtime_error("missing value for --d-boundary-rank");
            run_d = true;
            d_rank = std::stoi(argv[++i]);
        } else if (arg == "--d-exact-boundary-rank") {
            if (i + 1 >= argc) throw std::runtime_error("missing value for --d-exact-boundary-rank");
            run_d_exact = true;
            d_exact_rank = std::stoi(argv[++i]);
        } else if (arg == "--d-exact-boundary-range") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --d-exact-boundary-range");
            run_d_exact_range = true;
            d_exact_low_rank = std::stoi(argv[++i]);
            d_exact_high_rank = std::stoi(argv[++i]);
        } else if (arg == "--d-exact-moment-range") {
            if (i + 4 >= argc) throw std::runtime_error("missing values for --d-exact-moment-range");
            run_d_exact_moment_range = true;
            d_exact_moment_low_rank = std::stoi(argv[++i]);
            d_exact_moment_high_rank = std::stoi(argv[++i]);
            d_exact_moment_low = std::stoi(argv[++i]);
            d_exact_moment_high = std::stoi(argv[++i]);
        } else if (arg == "--d-o-even-moment-range") {
            if (i + 4 >= argc) throw std::runtime_error("missing values for --d-o-even-moment-range");
            run_d_o_even_moment_range = true;
            d_o_even_moment_low_rank = std::stoi(argv[++i]);
            d_o_even_moment_high_rank = std::stoi(argv[++i]);
            d_o_even_moment_low = std::stoi(argv[++i]);
            d_o_even_moment_high = std::stoi(argv[++i]);
        } else if (arg == "--d-capped-suffix-moment-range") {
            if (i + 3 >= argc) throw std::runtime_error("missing values for --d-capped-suffix-moment-range");
            run_d_capped_suffix = true;
            d_capped_suffix_rank = std::stoi(argv[++i]);
            d_capped_suffix_moment_low = std::stoi(argv[++i]);
            d_capped_suffix_moment_high = std::stoi(argv[++i]);
        } else if (arg == "--d-capped-suffix-moment-range-mod-prime") {
            if (i + 4 >= argc) throw std::runtime_error("missing values for --d-capped-suffix-moment-range-mod-prime");
            run_d_capped_suffix_mod_prime = true;
            d_capped_suffix_rank = std::stoi(argv[++i]);
            d_capped_suffix_moment_low = std::stoi(argv[++i]);
            d_capped_suffix_moment_high = std::stoi(argv[++i]);
            d_capped_suffix_mod_prime = static_cast<std::uint64_t>(std::stoull(argv[++i]));
        } else if (arg == "--d-interval-boundary-rank") {
            if (i + 1 >= argc) throw std::runtime_error("missing value for --d-interval-boundary-rank");
            run_d_interval = true;
            d_interval_rank = std::stoi(argv[++i]);
        } else if (arg == "--d-delta-certificate") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --d-delta-certificate");
            d_delta_rank = std::stoi(argv[++i]);
            d_delta_moment = std::stoi(argv[++i]);
        } else if (arg == "--d-delta-range-certificate") {
            if (i + 3 >= argc) throw std::runtime_error("missing values for --d-delta-range-certificate");
            d_delta_range_rank = std::stoi(argv[++i]);
            d_delta_range_low = std::stoi(argv[++i]);
            d_delta_range_high = std::stoi(argv[++i]);
        } else if (arg == "--d-delta-shared-range-certificate") {
            if (i + 4 >= argc) throw std::runtime_error("missing values for --d-delta-shared-range-certificate");
            run_d_delta_shared_range = true;
            d_delta_shared_low_rank = std::stoi(argv[++i]);
            d_delta_shared_high_rank = std::stoi(argv[++i]);
            d_delta_shared_low = std::stoi(argv[++i]);
            d_delta_shared_high = std::stoi(argv[++i]);
        } else if (arg == "--b-boundary-range") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --b-boundary-range");
            int low = std::stoi(argv[++i]);
            int high = std::stoi(argv[++i]);
            add_range('B', low, high);
        } else if (arg == "--c-boundary-range") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --c-boundary-range");
            int low = std::stoi(argv[++i]);
            int high = std::stoi(argv[++i]);
            add_range('C', low, high);
        } else if (arg == "--d-boundary-range") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --d-boundary-range");
            int low = std::stoi(argv[++i]);
            int high = std::stoi(argv[++i]);
            add_range('D', low, high);
        } else if (arg == "--chain-m") {
            if (i + 1 >= argc) throw std::runtime_error("missing value for --chain-m");
            chain_m = std::stoi(argv[++i]);
        } else if (arg == "--progress") {
            progress = true;
        } else if (arg == "--m22-remaining-lower-certificate") {
            run_m22_remaining = true;
        } else if (arg == "--m22-bc-remaining-lower-certificate") {
            run_m22_remaining_bc = true;
        } else if (arg == "--m22-d-remaining-lower-certificate") {
            run_m22_remaining_d = true;
        } else if (arg == "--m22-b5-lower-certificate") {
            run_m22_b5_lower = true;
        } else if (arg == "--m22-bc-row-lower-certificate") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --m22-bc-row-lower-certificate");
            std::string family = argv[++i];
            if (family != "B" && family != "C") {
                throw std::runtime_error("--m22-bc-row-lower-certificate family must be B or C");
            }
            m22_bc_single_family = family[0];
            m22_bc_single_rank = std::stoi(argv[++i]);
            run_m22_bc_single_lower = true;
        } else if (arg == "--bc-lower-range") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --bc-lower-range");
            generic_bc_low_rank = std::stoi(argv[++i]);
            generic_bc_high_rank = std::stoi(argv[++i]);
            run_generic_lower = true;
            generic_lower_bc = true;
            generic_lower_b_family = true;
            generic_lower_c_family = true;
        } else if (arg == "--b-lower-range") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --b-lower-range");
            generic_bc_low_rank = std::stoi(argv[++i]);
            generic_bc_high_rank = std::stoi(argv[++i]);
            run_generic_lower = true;
            generic_lower_bc = true;
            generic_lower_b_family = true;
        } else if (arg == "--c-lower-range") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --c-lower-range");
            generic_bc_low_rank = std::stoi(argv[++i]);
            generic_bc_high_rank = std::stoi(argv[++i]);
            run_generic_lower = true;
            generic_lower_bc = true;
            generic_lower_c_family = true;
        } else if (arg == "--d-lower-range") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --d-lower-range");
            generic_d_low_rank = std::stoi(argv[++i]);
            generic_d_high_rank = std::stoi(argv[++i]);
            run_generic_lower = true;
            generic_lower_d = true;
        } else if (arg == "--m22-bc-rank-pair-positive-certificate") {
            if (i + 1 >= argc) throw std::runtime_error("missing value for --m22-bc-rank-pair-positive-certificate");
            m22_bc_rank_pair = std::stoi(argv[++i]);
            run_m22_bc_rank_pair_positive = true;
        } else if (arg == "--bc-rank-pair-positive-certificate") {
            if (i + 1 >= argc) throw std::runtime_error("missing value for --bc-rank-pair-positive-certificate");
            m22_bc_rank_pair = std::stoi(argv[++i]);
            run_m22_bc_rank_pair_positive = true;
        } else if (arg == "--b11-post-m29-bridge-range-certificate") {
            run_b11_post_m29_bridge_range = true;
        } else if (arg == "--bc-family-delta") {
            if (i + 3 >= argc) throw std::runtime_error("missing values for --bc-family-delta");
            std::string family = argv[++i];
            if (family != "B" && family != "C") {
                throw std::runtime_error("--bc-family-delta family must be B or C");
            }
            bc_family_delta_family = family[0];
            bc_family_delta_rank = std::stoi(argv[++i]);
            bc_family_delta_moment = std::stoi(argv[++i]);
            run_bc_family_delta = true;
        } else if (arg == "--bc-rank-pair-delta") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --bc-rank-pair-delta");
            bc_rank_pair_delta_rank = std::stoi(argv[++i]);
            bc_rank_pair_delta_moment = std::stoi(argv[++i]);
            run_bc_rank_pair_delta = true;
        } else if (arg == "--bc-rank-pair-range-positive-certificate") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --bc-rank-pair-range-positive-certificate");
            bc_rank_pair_low = std::stoi(argv[++i]);
            bc_rank_pair_high = std::stoi(argv[++i]);
            run_bc_rank_pair_range_positive = true;
        } else if (arg == "--b2-c2-weyl-m22-certificate") {
            run_b2_c2_weyl = true;
            chain_m = 22;
        } else if (arg == "--b2-c2-weyl-certificate") {
            run_b2_c2_weyl = true;
        } else if (arg == "--b3-c3-weyl-certificate") {
            run_b3_c3_weyl = true;
        } else if (arg == "--b4-c4-weyl-certificate") {
            run_b4_c4_weyl = true;
        } else if (arg == "--b4-weyl-certificate") {
            run_b4_weyl = true;
        } else if (arg == "--c4-weyl-certificate") {
            run_c4_weyl = true;
        } else if (arg == "--d4-weyl-certificate") {
            run_d4_weyl = true;
        } else if (arg == "--rank4-weyl-mod-prime") {
            if (i + 2 >= argc) throw std::runtime_error("missing values for --rank4-weyl-mod-prime");
            std::string family = argv[++i];
            if (family != "B" && family != "C" && family != "D") {
                throw std::runtime_error("--rank4-weyl-mod-prime family must be B, C, or D");
            }
            rank4_weyl_mod_family = family[0];
            rank4_weyl_mod_prime = static_cast<std::uint64_t>(std::stoull(argv[++i]));
            run_rank4_weyl_mod = true;
        } else if (arg == "--m22-bc-exact-through") {
            if (i + 1 >= argc) throw std::runtime_error("missing value for --m22-bc-exact-through");
            m22_bc_exact_through_override = std::stoi(argv[++i]);
        } else if (arg == "--clear-coeff-cache-after-moment") {
            clear_coefficient_cache_after_moment = true;
        } else if (arg == "--forward-coefficients") {
            use_forward_coefficients = true;
        } else if (arg == "--reverse-target-sums") {
            use_reverse_target_sums = true;
        } else if (arg == "--include-positive-quadratic-exact") {
            include_positive_quadratic_exact = true;
        } else if (arg == "--no-predecessor-cache") {
            no_predecessor_cache = true;
        } else if (arg == "--stable-moments") {
            if (i + 1 >= argc) throw std::runtime_error("missing value for --stable-moments");
            stable_path = argv[++i];
        } else {
            throw std::runtime_error("unknown argument " + arg);
        }
    }

    if (
        !run_b && !run_c && !run_d && !run_d_exact && !run_d_exact_range
        && !run_d_exact_moment_range
        && !run_d_o_even_moment_range
        && !run_d_capped_suffix && !run_d_capped_suffix_mod_prime
        && !run_d_interval && !run_m22_remaining
        && !run_m22_remaining_bc && !run_m22_remaining_d && !run_m22_b5_lower
        && !run_m22_bc_single_lower
        && !run_m22_bc_rank_pair_positive
        && !run_b11_post_m29_bridge_range
        && !run_bc_family_delta
        && !run_bc_rank_pair_delta
        && !run_bc_rank_pair_range_positive
        && !run_b2_c2_weyl
        && !run_b3_c3_weyl
        && !run_b4_c4_weyl
        && !run_b4_weyl
        && !run_c4_weyl
        && !run_d4_weyl
        && !run_rank4_weyl_mod
        && d_delta_rank == 0
        && d_delta_range_rank == 0
        && !run_d_delta_shared_range
        && !run_generic_lower
        && rows.empty()
    ) {
        std::cerr << "usage: classical_boundary_certificate "
                  << "[--b-boundary-m15-rank16-certificate] "
                  << "[--c-boundary-m15-rank16-certificate] "
                  << "[--b-boundary-rank R] [--c-boundary-rank R] "
                  << "[--b-boundary-range LO HI] [--c-boundary-range LO HI] "
                  << "[--d-boundary-range LO HI] "
                  << "[--d-boundary-rank R] [--chain-m M] "
                  << "[--d-exact-boundary-rank R] "
                  << "[--d-exact-boundary-range LO HI] "
                  << "[--d-exact-moment-range RLO RHI MLO MHI] "
                  << "[--d-o-even-moment-range RLO RHI MLO MHI] "
                  << "[--d-capped-suffix-moment-range R MLO MHI] "
                  << "[--d-capped-suffix-moment-range-mod-prime R MLO MHI P] "
                  << "[--d-interval-boundary-rank R] "
                  << "[--d-delta-certificate R M] "
                  << "[--d-delta-range-certificate R LO HI] "
                  << "[--d-delta-shared-range-certificate RLO RHI MLO MHI] "
                  << "[--m22-remaining-lower-certificate] "
                  << "[--m22-bc-remaining-lower-certificate] "
                  << "[--m22-d-remaining-lower-certificate] "
                  << "[--m22-b5-lower-certificate] "
                  << "[--m22-bc-row-lower-certificate B|C R] "
                  << "[--bc-lower-range LO HI] "
                  << "[--b-lower-range LO HI] "
                  << "[--c-lower-range LO HI] "
                  << "[--d-lower-range LO HI] "
                  << "[--bc-rank-pair-positive-certificate R] "
                  << "[--b11-post-m29-bridge-range-certificate] "
                  << "[--bc-family-delta B|C R M] "
                  << "[--bc-rank-pair-delta R M] "
                  << "[--bc-rank-pair-range-positive-certificate LO HI] "
                  << "[--m22-bc-rank-pair-positive-certificate R] "
                  << "[--b2-c2-weyl-certificate] "
                  << "[--b2-c2-weyl-m22-certificate] "
                  << "[--b3-c3-weyl-certificate] "
                  << "[--b4-c4-weyl-certificate] "
                  << "[--b4-weyl-certificate] "
                  << "[--c4-weyl-certificate] "
                  << "[--d4-weyl-certificate] "
                  << "[--rank4-weyl-mod-prime B|C|D P] "
                  << "[--m22-bc-exact-through N] "
                  << "[--progress] [--no-predecessor-cache] "
                  << "[--clear-coeff-cache-after-moment] "
                  << "[--forward-coefficients] "
                  << "[--reverse-target-sums] "
                  << "[--include-positive-quadratic-exact]"
                  << std::endl;
        return 2;
    }

    try {
        use_predecessor_cache = !no_predecessor_cache;
        coefficient_progress = progress;
        if (run_m22_remaining) run_m22_remaining_lower_certificate(stable_path, progress, true, true);
        if (run_m22_remaining_bc) run_m22_remaining_lower_certificate(stable_path, progress, true, false);
        if (run_m22_remaining_d) run_m22_remaining_lower_certificate(stable_path, progress, false, true);
        if (run_m22_b5_lower) {
            run_m22_remaining_lower_certificate(
                stable_path,
                progress,
                true,
                false,
                5,
                5,
                true,
                false,
                m22_bc_exact_through_override,
                include_positive_quadratic_exact
            );
        }
        if (run_m22_bc_single_lower) {
            run_m22_remaining_lower_certificate(
                stable_path,
                progress,
                true,
                false,
                m22_bc_single_rank,
                m22_bc_single_rank,
                m22_bc_single_family == 'B',
                m22_bc_single_family == 'C',
                m22_bc_exact_through_override,
                include_positive_quadratic_exact
            );
        }
        if (run_m22_bc_rank_pair_positive) {
            run_m22_bc_rank_pair_positive_certificate(
                stable_path,
                progress,
                m22_bc_rank_pair,
                m22_bc_exact_through_override,
                chain_m
            );
        }
        if (run_b11_post_m29_bridge_range) {
            run_b11_post_m29_bridge_range_certificate(
                stable_path,
                progress,
                include_positive_quadratic_exact
            );
        }
        if (run_bc_family_delta) {
            run_bc_family_delta_certificate(
                progress,
                bc_family_delta_family,
                bc_family_delta_rank,
                bc_family_delta_moment
            );
        }
        if (run_bc_rank_pair_delta) {
            run_bc_rank_pair_delta_certificate(
                progress,
                bc_rank_pair_delta_rank,
                bc_rank_pair_delta_moment
            );
        }
        if (run_bc_rank_pair_range_positive) {
            run_bc_rank_pair_range_positive_certificate(
                stable_path,
                progress,
                bc_rank_pair_low,
                bc_rank_pair_high,
                m22_bc_exact_through_override,
                chain_m
            );
        }
        if (run_b2_c2_weyl) run_b2_c2_weyl_certificate(stable_path, progress, chain_m);
        if (run_b3_c3_weyl) run_b3_c3_weyl_certificate(stable_path, progress, chain_m);
        if (run_b4_c4_weyl) run_b4_c4_weyl_certificate(stable_path, progress, chain_m);
        if (run_b4_weyl) run_rank4_weyl_family_certificate('B', stable_path, progress, chain_m);
        if (run_c4_weyl) run_rank4_weyl_family_certificate('C', stable_path, progress, chain_m);
        if (run_d4_weyl) run_rank4_weyl_family_certificate('D', stable_path, progress, chain_m);
        if (run_rank4_weyl_mod) {
            run_rank4_weyl_mod_prime_certificate(
                rank4_weyl_mod_family,
                chain_m,
                rank4_weyl_mod_prime,
                progress
            );
        }
        if (d_delta_rank != 0) {
            if (d_delta_moment < d_delta_rank) {
                throw std::runtime_error("D delta moment is below the first boundary");
            }
            BigInt correction = d_determinant_correction(d_delta_rank, d_delta_moment, progress);
            if (correction < 0) throw std::runtime_error("D determinant correction is negative");
            std::cout << "D_" << d_delta_rank << " Delta_" << d_delta_moment
                      << " = " << correction << std::endl;
        }
        if (d_delta_range_rank != 0) {
            if (d_delta_range_low > d_delta_range_high) {
                throw std::runtime_error("empty D delta moment range");
            }
            if (d_delta_range_low < d_delta_range_rank) {
                throw std::runtime_error("D delta range begins below the first boundary");
            }
            for (int moment = d_delta_range_low; moment <= d_delta_range_high; ++moment) {
                BigInt correction = d_determinant_correction(d_delta_range_rank, moment, progress);
                if (correction < 0) throw std::runtime_error("D determinant correction is negative");
                std::cout << "D_" << d_delta_range_rank << " Delta_" << moment
                          << " = " << correction << std::endl;
            }
        }
        if (run_d_delta_shared_range) {
            run_d_delta_shared_range_certificate(
                d_delta_shared_low_rank,
                d_delta_shared_high_rank,
                d_delta_shared_low,
                d_delta_shared_high,
                progress
            );
        }
        if (run_generic_lower) {
            run_m22_remaining_lower_certificate(
                stable_path,
                progress,
                generic_lower_bc,
                generic_lower_d,
                generic_bc_low_rank,
                generic_bc_high_rank,
                generic_lower_b_family,
                generic_lower_c_family,
                m22_bc_exact_through_override,
                include_positive_quadratic_exact,
                chain_m,
                generic_d_low_rank,
                generic_d_high_rank
            );
        }
        if (run_b) run_boundary('B', b_rank, chain_m, stable_path, progress);
        if (run_c) run_boundary('C', c_rank, chain_m, stable_path, progress);
        if (run_d) run_d_boundary(d_rank, chain_m, stable_path, progress);
        if (run_d_exact) run_d_exact_boundary(d_exact_rank, chain_m, stable_path, progress);
        if (run_d_interval) run_d_interval_boundary(d_interval_rank, chain_m, stable_path, progress);
        if (run_d_exact_range) {
            if (d_exact_low_rank > d_exact_high_rank) {
                throw std::runtime_error("empty D exact rank range");
            }
            run_d_exact_boundary_shared_range(
                d_exact_low_rank,
                d_exact_high_rank,
                chain_m,
                stable_path,
                progress
            );
        }
        if (run_d_exact_moment_range) {
            run_d_exact_moment_shared_range(
                d_exact_moment_low_rank,
                d_exact_moment_high_rank,
                d_exact_moment_low,
                d_exact_moment_high,
                stable_path,
                progress
            );
        }
        if (run_d_o_even_moment_range) {
            run_d_o_even_moment_shared_range(
                d_o_even_moment_low_rank,
                d_o_even_moment_high_rank,
                d_o_even_moment_low,
                d_o_even_moment_high,
                stable_path,
                progress
            );
        }
        if (run_d_capped_suffix) {
            run_d_capped_suffix_moment_range(
                d_capped_suffix_rank,
                d_capped_suffix_moment_low,
                d_capped_suffix_moment_high,
                stable_path,
                progress
            );
        }
        if (run_d_capped_suffix_mod_prime) {
            run_d_capped_suffix_moment_range_mod_prime(
                d_capped_suffix_rank,
                d_capped_suffix_moment_low,
                d_capped_suffix_moment_high,
                d_capped_suffix_mod_prime,
                progress
            );
        }
        std::size_t row_index = 0;
        for (const auto& [family, rank] : rows) {
            ++row_index;
            std::cout << "batch row " << row_index << "/" << rows.size()
                      << " " << family << "_" << rank << std::endl;
            if (family == 'B' || family == 'C') {
                run_boundary(family, rank, chain_m, stable_path, progress);
            } else if (family == 'D') {
                run_d_boundary(rank, chain_m, stable_path, progress);
            } else {
                throw std::runtime_error("unknown batch family");
            }
        }
        std::cout.flush();
        std::cerr.flush();
        std::_Exit(0);
    } catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << std::endl;
        return 1;
    }
}
