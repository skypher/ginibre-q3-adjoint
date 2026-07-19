#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "ankerl/unordered_dense.h"

namespace {

constexpr int kMaxRows = 64;
constexpr long double kNegInf = -std::numeric_limits<long double>::infinity();
constexpr long double kLog10Half = -0.301029995663981195213738894724493027L;

struct PackedPartition {
    std::array<std::uint16_t, kMaxRows> parts{};
    std::uint16_t length = 0;

    bool operator==(const PackedPartition& other) const {
        if (length != other.length) return false;
        for (int i = 0; i < length; ++i) {
            if (parts[static_cast<std::size_t>(i)] != other.parts[static_cast<std::size_t>(i)]) {
                return false;
            }
        }
        return true;
    }
};

struct PackedPartitionHash {
    std::size_t operator()(const PackedPartition& partition) const {
        std::uint64_t h = 1469598103934665603ull ^ partition.length;
        for (int i = 0; i < partition.length; ++i) {
            std::uint64_t x = partition.parts[static_cast<std::size_t>(i)];
            h ^= x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
            h *= 1099511628211ull;
        }
        return static_cast<std::size_t>(h);
    }
};

struct SignedLog {
    int sign = 0;
    long double log_abs = kNegInf;
};

void add_log(long double& acc, long double value) {
    if (!std::isfinite(value)) return;
    if (!std::isfinite(acc)) {
        acc = value;
        return;
    }
    const long double hi = std::max(acc, value);
    const long double lo = std::min(acc, value);
    acc = hi - lo > 80.0L ? hi : hi + std::log10(1.0L + std::pow(10.0L, lo - hi));
}

void add_signed_log(SignedLog& acc, int sign, long double log_abs) {
    if (sign == 0 || !std::isfinite(log_abs)) return;
    const int term_sign = sign > 0 ? 1 : -1;
    if (acc.sign == 0) {
        acc.sign = term_sign;
        acc.log_abs = log_abs;
        return;
    }
    const long double hi = std::max(acc.log_abs, log_abs);
    const long double lo = std::min(acc.log_abs, log_abs);
    if (acc.sign == term_sign) {
        acc.log_abs = hi - lo > 80.0L
            ? hi
            : hi + std::log10(1.0L + std::pow(10.0L, lo - hi));
        return;
    }
    if (std::fabs(acc.log_abs - log_abs) < 1.0e-18L) {
        acc.sign = 0;
        acc.log_abs = kNegInf;
        return;
    }
    const int hi_sign = acc.log_abs > log_abs ? acc.sign : term_sign;
    if (hi - lo > 80.0L) {
        acc.sign = hi_sign;
        acc.log_abs = hi;
        return;
    }
    const long double diff = 1.0L - std::pow(10.0L, lo - hi);
    if (!(diff > 0.0L)) {
        acc.sign = 0;
        acc.log_abs = kNegInf;
        return;
    }
    acc.sign = hi_sign;
    acc.log_abs = hi + std::log10(diff);
}

std::string display(long double value) {
    if (std::isinf(value)) return value > 0 ? "+inf" : "-inf";
    std::ostringstream out;
    out << std::setprecision(18) << static_cast<double>(value);
    return out.str();
}

std::vector<long double> stable_moment_logs(int max_index) {
    std::vector<SignedLog> moments(static_cast<std::size_t>(max_index + 1));
    moments[0] = SignedLog{1, 0.0L};
    if (max_index >= 1) moments[1] = SignedLog{0, kNegInf};
    for (int index = 1; index < max_index; ++index) {
        SignedLog next;
        if (moments[static_cast<std::size_t>(index)].sign != 0) {
            add_signed_log(
                next,
                moments[static_cast<std::size_t>(index)].sign,
                std::log10(static_cast<long double>(index))
                    + moments[static_cast<std::size_t>(index)].log_abs
            );
        }
        if (moments[static_cast<std::size_t>(index - 1)].sign != 0) {
            add_signed_log(
                next,
                moments[static_cast<std::size_t>(index - 1)].sign,
                std::log10(static_cast<long double>(index))
                    + moments[static_cast<std::size_t>(index - 1)].log_abs
            );
        }
        if (index >= 2 && moments[static_cast<std::size_t>(index - 2)].sign != 0) {
            add_signed_log(
                next,
                -moments[static_cast<std::size_t>(index - 2)].sign,
                std::log10(static_cast<long double>(index))
                    + std::log10(static_cast<long double>(index - 1))
                    - std::log10(2.0L)
                    + moments[static_cast<std::size_t>(index - 2)].log_abs
            );
        }
        if (next.sign <= 0) {
            std::cerr << "stable moment recurrence lost positivity at index "
                      << (index + 1) << "\n";
            std::exit(1);
        }
        moments[static_cast<std::size_t>(index + 1)] = next;
    }
    std::vector<long double> logs(static_cast<std::size_t>(max_index + 1), kNegInf);
    for (int index = 0; index <= max_index; ++index) {
        logs[static_cast<std::size_t>(index)] = moments[static_cast<std::size_t>(index)].log_abs;
    }
    return logs;
}

bool row_can_receive(const PackedPartition& partition, int row) {
    return row == 0 || partition.parts[static_cast<std::size_t>(row - 1)]
        > partition.parts[static_cast<std::size_t>(row)];
}

bool all_rows_even(const PackedPartition& partition) {
    for (int i = 0; i < partition.length; ++i) {
        if (partition.parts[static_cast<std::size_t>(i)] % 2 != 0) return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    int rank = 14;
    int max_moment = 925;
    int progress_step = 25;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--rank" && i + 1 < argc) {
            rank = std::atoi(argv[++i]);
        } else if (arg == "--max-moment" && i + 1 < argc) {
            max_moment = std::atoi(argv[++i]);
        } else if (arg == "--progress-step" && i + 1 < argc) {
            progress_step = std::atoi(argv[++i]);
        } else {
            std::cerr << "usage: " << argv[0]
                      << " [--rank R] [--max-moment M] [--progress-step K]\n";
            return 2;
        }
    }
    const int row_cap = 2 * rank + 1;
    if (rank < 1 || row_cap > kMaxRows || max_moment < 0 || progress_step < 1) {
        std::cerr << "invalid B row-cap scan arguments\n";
        return 2;
    }

    const std::vector<long double> stable_logs = stable_moment_logs(max_moment);
    ankerl::unordered_dense::map<PackedPartition, long double, PackedPartitionHash> current;
    current.emplace(PackedPartition{}, 0.0L);

    long double worst_kept_minus_stable = 0.0L;
    int worst_moment = 0;
    std::size_t worst_states = 1;
    std::cout << "B row-capped Pieri kept-mass diagnostic\n"
              << "DIAGNOSTIC_ONLY exact recurrence in log domain; not proof evidence\n"
              << "rank=" << rank
              << " row_cap=" << row_cap
              << " max_moment=" << max_moment
              << " progress_step=" << progress_step
              << std::endl;

    for (int moment = 1; moment <= max_moment; ++moment) {
        ankerl::unordered_dense::map<PackedPartition, long double, PackedPartitionHash> next;
        next.reserve(current.size() + current.size() / 2 + 16);
        for (const auto& [partition, log_weight] : current) {
            const int slots = std::min(row_cap, static_cast<int>(partition.length) + 2);
            for (int first = 0; first < slots; ++first) {
                if (!row_can_receive(partition, first)) continue;
                for (int second = first + 1; second < slots; ++second) {
                    if (second != first + 1 && !row_can_receive(partition, second)) continue;
                    PackedPartition candidate = partition;
                    const int candidate_length = std::max<int>(candidate.length, second + 1);
                    candidate.parts[static_cast<std::size_t>(first)] += 1;
                    candidate.parts[static_cast<std::size_t>(second)] += 1;
                    candidate.length = static_cast<std::uint16_t>(candidate_length);
                    auto [iterator, inserted] = next.emplace(candidate, log_weight);
                    if (!inserted) add_log(iterator->second, log_weight);
                }
            }
        }
        current = std::move(next);

        long double kept_even_log = kNegInf;
        for (const auto& [partition, log_weight] : current) {
            if (all_rows_even(partition)) add_log(kept_even_log, log_weight);
        }
        const long double kept_minus_stable =
            kept_even_log - stable_logs[static_cast<std::size_t>(moment)];
        if (moment == 1 || kept_minus_stable < worst_kept_minus_stable) {
            worst_kept_minus_stable = kept_minus_stable;
            worst_moment = moment;
            worst_states = current.size();
        }
        if (moment == max_moment || moment % progress_step == 0 || moment == 2 * rank + 2) {
            long double omitted_ratio = 0.0L;
            if (kept_minus_stable < -1.0e-10L) {
                omitted_ratio = 1.0L - std::pow(10.0L, kept_minus_stable);
            }
            std::cout << "progress moment=" << moment
                      << "/" << max_moment
                      << " states=" << current.size()
                      << " kept_minus_stable_log10=" << display(kept_minus_stable)
                      << " omitted_ratio_est=" << display(omitted_ratio)
                      << " half_ok=" << (kept_minus_stable >= kLog10Half)
                      << std::endl;
            std::fflush(stdout);
        }
    }
    std::cout << "SUMMARY rank=" << rank
              << " row_cap=" << row_cap
              << " max_moment=" << max_moment
              << " worst_moment=" << worst_moment
              << " worst_states=" << worst_states
              << " worst_kept_minus_stable_log10=" << display(worst_kept_minus_stable)
              << " half_ok=" << (worst_kept_minus_stable >= kLog10Half)
              << std::endl;
    return 0;
}
