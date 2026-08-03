#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = __int128_t;

struct Evaluation {
    bool active = false;
    bool stronger_failure = false;
    Integer positive_minor_sum = 0;
    Integer nonzero_anchor_payment = 0;
    Integer anchored_sum = 0;
    Integer anchor_weight = 0;
    std::uint64_t anchor = 0U;
};

std::uint64_t parse_u64(const char* text, const char* name) {
    char* end = nullptr;
    const unsigned long long value = std::strtoull(text, &end, 10);
    if (end == text || *end != '\0') {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<std::uint64_t>(value);
}

int parse_positive(const char* text, const char* name) {
    const std::uint64_t value = parse_u64(text, name);
    if (value == 0U
        || value > static_cast<std::uint64_t>(
            std::numeric_limits<int>::max()
        )) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(value);
}

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

std::string integer_string(Integer value) {
    if (value == 0) {
        return "0";
    }
    const bool negative = value < 0;
    if (negative) {
        value = -value;
    }
    std::string result;
    while (value != 0) {
        const int digit = static_cast<int>(value % 10);
        result.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    if (negative) {
        result.push_back('-');
    }
    std::reverse(result.begin(), result.end());
    return result;
}

std::vector<int> retained_outputs(int level, int parity) {
    std::vector<int> result;
    for (int index = 0; index <= 4; ++index) {
        const int output = parity + 2 * index;
        if (output <= level) {
            result.push_back(output);
        }
    }
    const int reflected_parity = (level + parity) & 1;
    for (int index = 0; index <= 4; ++index) {
        const int output = level - reflected_parity - 2 * index;
        if (output >= 0
            && std::find(result.begin(), result.end(), output)
                == result.end()) {
            result.push_back(output);
        }
    }
    return result;
}

Evaluation evaluate(
    int level,
    const std::vector<int>& labels,
    std::uint64_t minus_mask
) {
    const int factors = static_cast<int>(labels.size());
    const std::uint64_t state_count = 1ULL << factors;
    const std::uint64_t full = state_count - 1U;
    const std::size_t width = static_cast<std::size_t>(level + 1);
    std::vector<std::int64_t> states(
        static_cast<std::size_t>(state_count) * width, 0
    );
    const auto at = [&states, width](std::uint64_t mask, int output)
        -> std::int64_t& {
        return states[static_cast<std::size_t>(mask) * width
            + static_cast<std::size_t>(output)];
    };
    at(0U, 0) = 1;
    for (std::uint64_t mask = 1U; mask < state_count; ++mask) {
        const int index = static_cast<int>(std::countr_zero(mask));
        const std::uint64_t previous = mask & (mask - 1U);
        const int label = labels[static_cast<std::size_t>(index)];
        for (int source = 0; source <= level; ++source) {
            const std::int64_t multiplicity = at(previous, source);
            if (multiplicity == 0) {
                continue;
            }
            const int upper = std::min(
                source + label, 2 * level - source - label
            );
            for (int output = std::abs(source - label);
                 output <= upper; output += 2) {
                std::int64_t updated = 0;
                if (__builtin_add_overflow(
                        at(mask, output), multiplicity, &updated
                    )) {
                    throw std::overflow_error("fusion multiplicity overflow");
                }
                at(mask, output) = updated;
            }
        }
    }

    Evaluation result;
    Integer maximum_weight = 0;
    for (std::uint64_t mask = 1U; mask < full; ++mask) {
        const std::uint64_t complement = full ^ mask;
        if (mask > complement) {
            continue;
        }
        const Integer weight = static_cast<Integer>(at(mask, 0))
            * static_cast<Integer>(at(complement, 0));
        const bool negative =
            (std::popcount(mask & minus_mask) & 1) != 0;
        if (negative && weight > maximum_weight) {
            maximum_weight = weight;
            result.anchor = mask;
        }
    }
    if (maximum_weight == 0) {
        return result;
    }
    result.active = true;
    result.anchor_weight = maximum_weight;

    const std::uint64_t complement = full ^ result.anchor;
    int parity = 0;
    for (int index = 0; index < factors; ++index) {
        if (((result.anchor >> index) & 1U) != 0U) {
            parity ^= labels[static_cast<std::size_t>(index)] & 1;
        }
    }
    Integer retained = 0;
    for (const int output : retained_outputs(level, parity)) {
        retained += static_cast<Integer>(at(result.anchor, output))
            * static_cast<Integer>(at(complement, output));
    }
    result.nonzero_anchor_payment = retained - maximum_weight;

    for (std::uint64_t mask = 1U; mask < full; ++mask) {
        const std::uint64_t other_complement = full ^ mask;
        if (mask > other_complement || mask == result.anchor) {
            continue;
        }
        const bool negative =
            (std::popcount(mask & minus_mask) & 1) != 0;
        if (!negative) {
            continue;
        }
        const std::uint64_t switched = result.anchor ^ mask;
        const std::uint64_t switched_complement = full ^ switched;
        const Integer original = static_cast<Integer>(at(mask, 0))
            * static_cast<Integer>(at(other_complement, 0));
        const Integer replacement = static_cast<Integer>(at(switched, 0))
            * static_cast<Integer>(at(switched_complement, 0));
        const Integer defect = original - replacement;
        result.anchored_sum += defect;
        if (defect > 0) {
            result.positive_minor_sum += defect;
        }
    }
    result.stronger_failure =
        result.positive_minor_sum > result.nonzero_anchor_payment;
    return result;
}

void print_case(
    int level,
    const std::vector<int>& labels,
    std::uint64_t minus_mask,
    const Evaluation& result
) {
    std::cout << " level=" << level << " labels=[";
    for (std::size_t index = 0U; index < labels.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << labels[index];
    }
    std::cout << "] minus_mask=" << minus_mask
              << " anchor=" << result.anchor
              << " anchor_weight=" << integer_string(result.anchor_weight)
              << " positive_minor_sum="
              << integer_string(result.positive_minor_sum)
              << " anchored_sum=" << integer_string(result.anchored_sum)
              << " nonzero_anchor_payment="
              << integer_string(result.nonzero_anchor_payment)
              << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--replay") {
            const int level = 11;
            const std::vector<int> labels{
                1, 5, 11, 9, 3, 4, 8, 10, 7, 6, 2
            };
            constexpr std::uint64_t minus_mask = 670U;
            const Evaluation result = evaluate(level, labels, minus_mask);
            const bool passed = result.active
                && result.anchor == 292U
                && result.anchor_weight == Integer{105}
                && result.positive_minor_sum == Integer{1585}
                && result.anchored_sum == Integer{-449}
                && result.nonzero_anchor_payment == Integer{1482}
                && result.stronger_failure;
            std::cout << "SU2_POSITIVE_MINOR_PAYMENT_REPLAY result="
                      << (passed ? "PASS" : "FAIL");
            print_case(level, labels, minus_mask, result);
            return passed ? EXIT_SUCCESS : EXIT_FAILURE;
        }
        if (argc != 5) {
            throw std::runtime_error(
                "usage: probe_su2_positive_minor_payment "
                "seed trials maximum_level maximum_factors"
                " | --replay"
            );
        }
        const std::uint64_t seed = parse_u64(argv[1], "seed");
        const std::uint64_t trials = parse_u64(argv[2], "trials");
        const int maximum_level = parse_positive(argv[3], "maximum level");
        const int maximum_factors = parse_positive(argv[4], "maximum factors");
        if (trials == 0U || maximum_level < 2 || maximum_factors < 2
            || maximum_factors > 20 || maximum_factors > maximum_level) {
            throw std::runtime_error(
                "require trials>0, level>=2, and "
                "2<=factors<=min(20,level)"
            );
        }

        std::uint64_t active = 0U;
        for (std::uint64_t trial = 0U; trial < trials; ++trial) {
            std::uint64_t random = splitmix64(
                seed ^ (trial * 0xd1342543de82ef95ULL)
            );
            const int level = 2 + static_cast<int>(
                random % static_cast<std::uint64_t>(maximum_level - 1)
            );
            const int factor_limit = std::min(maximum_factors, level);
            random = splitmix64(random);
            const int factors = 2 + static_cast<int>(
                random % static_cast<std::uint64_t>(factor_limit - 1)
            );
            std::vector<int> pool(static_cast<std::size_t>(level));
            for (int label = 1; label <= level; ++label) {
                pool[static_cast<std::size_t>(label - 1)] = label;
            }
            for (int index = 0; index < factors; ++index) {
                random = splitmix64(random);
                const int other = index + static_cast<int>(
                    random % static_cast<std::uint64_t>(level - index)
                );
                std::swap(
                    pool[static_cast<std::size_t>(index)],
                    pool[static_cast<std::size_t>(other)]
                );
            }
            std::vector<int> labels(
                pool.begin(), pool.begin() + factors
            );
            int parity = 0;
            for (const int label : labels) {
                parity ^= label & 1;
            }
            if (parity != 0) {
                continue;
            }
            random = splitmix64(random);
            std::uint64_t minus_mask = random
                & ((1ULL << factors) - 1U);
            if ((std::popcount(minus_mask) & 1) != 0) {
                minus_mask ^= 1U;
            }
            if (minus_mask == 0U) {
                minus_mask = 3U;
            }
            const Evaluation result = evaluate(level, labels, minus_mask);
            if (!result.active) {
                continue;
            }
            ++active;
            if (result.stronger_failure) {
                std::cout << "SU2_POSITIVE_MINOR_PAYMENT "
                          << "separate_positive_minors=FAIL";
                print_case(level, labels, minus_mask, result);
                return EXIT_SUCCESS;
            }
        }
        std::cout << "SU2_POSITIVE_MINOR_PAYMENT "
                  << "separate_positive_minors=NO_COUNTEREXAMPLE"
                  << " trials=" << trials << " active=" << active << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_POSITIVE_MINOR_PAYMENT FAILURE: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
