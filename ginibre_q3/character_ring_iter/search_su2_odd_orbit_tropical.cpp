#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct GeneratorData {
    int sign = 0;
    int minimum_power = 0;
    int direction = 0;
};

struct PairData {
    int left = 0;
    int right = 0;
    int term_sign = 1;
    long double score = 0.0L;
};

class XorShift64 {
public:
    explicit XorShift64(std::uint64_t seed) : state_(seed) {}

    std::uint64_t next() {
        state_ ^= state_ << 13U;
        state_ ^= state_ >> 7U;
        state_ ^= state_ << 17U;
        return state_;
    }

    int bounded(int upper_exclusive) {
        if (upper_exclusive <= 0) {
            throw std::runtime_error("invalid random bound");
        }
        return static_cast<int>(
            next() % static_cast<std::uint64_t>(upper_exclusive)
        );
    }

private:
    std::uint64_t state_;
};

void print_candidate(
    int rank,
    const std::vector<GeneratorData>& generators,
    const PairData& winner,
    long double margin,
    std::uint64_t trial
) {
    std::cout << "SU2_ODD_ORBIT_TROPICAL result=FAIL"
              << " rank=" << rank
              << " level=" << 2 * rank - 1
              << " trial=" << trial
              << " winner_pair=" << winner.left + 1
              << ',' << winner.right + 1
              << " margin=" << margin
              << " word=";
    for (std::size_t label = 0; label < generators.size(); ++label) {
        const GeneratorData& generator = generators[label];
        if (generator.sign == 0) {
            continue;
        }
        std::cout << '[' << (generator.sign == 1 ? '+' : '-')
                  << label + 1 << '^' << generator.minimum_power
                  << "+2N*" << generator.direction << ']';
    }
    std::cout << '\n';
}

bool search_rank(
    int rank,
    std::uint64_t trials,
    int maximum_direction,
    XorShift64& random
) {
    const int order = 2 * rank + 1;
    const int label_count = rank - 1;
    const int pair_count = rank * (rank - 1) / 2;
    const long double pi = std::acos(-1.0L);
    std::vector<std::vector<long double>> values(
        static_cast<std::size_t>(rank),
        std::vector<long double>(static_cast<std::size_t>(label_count), 0.0L)
    );
    for (int node = 0; node < rank; ++node) {
        const long double angle
            = pi * static_cast<long double>(node + 1)
            / static_cast<long double>(order);
        for (int label = 0; label < label_count; ++label) {
            values[static_cast<std::size_t>(node)]
                  [static_cast<std::size_t>(label)]
                = std::sin(
                      static_cast<long double>(2 * label + 3) * angle
                  ) / std::sin(angle);
        }
    }

    std::vector<std::vector<std::array<long double, 2>>> logarithms(
        static_cast<std::size_t>(pair_count),
        std::vector<std::array<long double, 2>>(
            static_cast<std::size_t>(label_count)
        )
    );
    std::vector<std::vector<std::array<int, 2>>> base_signs(
        static_cast<std::size_t>(pair_count),
        std::vector<std::array<int, 2>>(
            static_cast<std::size_t>(label_count)
        )
    );
    int pair = 0;
    for (int left = 0; left < rank; ++left) {
        for (int right = left + 1; right < rank; ++right) {
            for (int label = 0; label < label_count; ++label) {
                for (int sign_index = 0; sign_index < 2; ++sign_index) {
                    const int sign = sign_index == 0 ? 1 : -1;
                    const long double base
                        = values[static_cast<std::size_t>(left)]
                                [static_cast<std::size_t>(label)]
                        + static_cast<long double>(sign)
                            * values[static_cast<std::size_t>(right)]
                                  [static_cast<std::size_t>(label)];
                    base_signs[static_cast<std::size_t>(pair)]
                              [static_cast<std::size_t>(label)]
                              [static_cast<std::size_t>(sign_index)]
                        = base < 0.0L ? -1 : 1;
                    logarithms[static_cast<std::size_t>(pair)]
                              [static_cast<std::size_t>(label)]
                              [static_cast<std::size_t>(sign_index)]
                        = std::abs(base) < 1.0e-18L
                        ? -std::numeric_limits<long double>::infinity()
                        : 2.0L * std::log(std::abs(base));
                }
            }
            ++pair;
        }
    }

    for (std::uint64_t trial = 0U; trial < trials; ++trial) {
        std::vector<GeneratorData> generators(
            static_cast<std::size_t>(label_count)
        );
        bool has_minus = false;
        bool has_direction = false;
        int total_minus_parity = 0;
        int first_minus = -1;
        for (int label = 0; label < label_count; ++label) {
            GeneratorData& generator
                = generators[static_cast<std::size_t>(label)];
            const int support = random.bounded(3);
            if (support == 0) {
                continue;
            }
            generator.sign = support == 1 ? 1 : -1;
            generator.minimum_power = random.bounded(2) == 0 ? 1 : 2;
            generator.direction = random.bounded(maximum_direction + 1);
            has_direction = has_direction || generator.direction != 0;
            if (generator.sign == -1) {
                has_minus = true;
                if (first_minus < 0) {
                    first_minus = label;
                }
                if ((generator.minimum_power & 1) != 0) {
                    total_minus_parity ^= 1;
                }
            }
        }
        if (!has_minus || !has_direction) {
            continue;
        }
        if (total_minus_parity != 0) {
            GeneratorData& generator = generators[
                static_cast<std::size_t>(first_minus)
            ];
            generator.minimum_power
                = generator.minimum_power == 1 ? 2 : 1;
        }

        std::vector<PairData> pair_data;
        pair_data.reserve(static_cast<std::size_t>(pair_count));
        pair = 0;
        for (int left = 0; left < rank; ++left) {
            for (int right = left + 1; right < rank; ++right) {
                PairData item;
                item.left = left;
                item.right = right;
                for (int label = 0; label < label_count; ++label) {
                    const GeneratorData& generator
                        = generators[static_cast<std::size_t>(label)];
                    if (generator.sign == 0) {
                        continue;
                    }
                    const int sign_index = generator.sign == 1 ? 0 : 1;
                    const int base_sign
                        = base_signs[static_cast<std::size_t>(pair)]
                                    [static_cast<std::size_t>(label)]
                                    [static_cast<std::size_t>(sign_index)];
                    if (base_sign == -1
                        && (generator.minimum_power & 1) != 0) {
                        item.term_sign = -item.term_sign;
                    }
                    const long double logarithm
                        = logarithms[static_cast<std::size_t>(pair)]
                                    [static_cast<std::size_t>(label)]
                                    [static_cast<std::size_t>(sign_index)];
                    if (!std::isfinite(logarithm)) {
                        item.score
                            = -std::numeric_limits<long double>::infinity();
                        break;
                    }
                    if (generator.direction != 0) {
                        item.score += static_cast<long double>(
                            generator.direction
                        ) * logarithm;
                    }
                }
                pair_data.push_back(item);
                ++pair;
            }
        }
        std::sort(
            pair_data.begin(),
            pair_data.end(),
            [](const PairData& left, const PairData& right) {
                return left.score > right.score;
            }
        );
        if (pair_data.size() >= 2U
            && pair_data[0].term_sign < 0
            && std::isfinite(pair_data[0].score)) {
            const long double margin
                = pair_data[0].score - pair_data[1].score;
            if (margin > 1.0e-10L) {
                print_candidate(
                    rank, generators, pair_data[0], margin, trial
                );
                return false;
            }
        }
    }
    std::cout << "progress rank=" << rank
              << " level=" << 2 * rank - 1
              << " trials=" << trials
              << " result=PASS\n";
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 4) {
            std::cerr << "usage: " << argv[0]
                      << " MAXIMUM_RANK TRIALS MAXIMUM_DIRECTION\n";
            return EXIT_FAILURE;
        }
        const int maximum_rank = std::stoi(argv[1]);
        const std::uint64_t trials = std::stoull(argv[2]);
        const int maximum_direction = std::stoi(argv[3]);
        if (maximum_rank < 5 || trials == 0U || maximum_direction < 1) {
            throw std::runtime_error("invalid tropical search bounds");
        }
        XorShift64 random(0x9e3779b97f4a7c15ULL);
        for (int rank = 5; rank <= maximum_rank; ++rank) {
            if (!search_rank(rank, trials, maximum_direction, random)) {
                return EXIT_FAILURE;
            }
        }
        std::cout << "SU2_ODD_ORBIT_TROPICAL ranks=5.."
                  << maximum_rank
                  << " trials_per_rank=" << trials
                  << " maximum_direction=" << maximum_direction
                  << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
