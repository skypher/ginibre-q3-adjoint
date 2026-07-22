#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

using boost::multiprecision::cpp_int;

namespace {

struct Witness {
    bool initialized = false;
    cpp_int value = 0;
    std::vector<int> minus;
    std::vector<int> plus;
};

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

std::uint64_t pair_key(int left, int right) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(left)) << 32U)
        | static_cast<std::uint32_t>(right);
}

std::pair<int, int> decode_pair(std::uint64_t key) {
    return {
        static_cast<int>(static_cast<std::uint32_t>(key >> 32U)),
        static_cast<int>(static_cast<std::uint32_t>(key))
    };
}

template <class Function>
void for_each_fusion_output(int level, int left, int right, Function function) {
    const int upper = std::min(left + right, 2 * level - left - right);
    for (int output = std::abs(left - right); output <= upper; output += 2) {
        function(output);
    }
}

cpp_int q3_dynamic(
    int level,
    const std::vector<int>& minus,
    const std::vector<int>& plus
) {
    std::unordered_map<std::uint64_t, cpp_int> states;
    states.emplace(pair_key(0, 0), 1);

    auto apply = [&states, level](int label, int sign) {
        std::unordered_map<std::uint64_t, cpp_int> next;
        for (const auto& [key, coefficient] : states) {
            const auto [left, right] = decode_pair(key);
            for_each_fusion_output(level, left, label, [&](int output) {
                next[pair_key(output, right)] += coefficient;
            });
            for_each_fusion_output(level, right, label, [&](int output) {
                next[pair_key(left, output)] += sign * coefficient;
            });
        }
        states = std::move(next);
    };

    for (int label : minus) {
        apply(label, -1);
    }
    for (int label : plus) {
        apply(label, 1);
    }
    const auto found = states.find(pair_key(0, 0));
    return found == states.end() ? cpp_int(0) : found->second;
}

cpp_int q3_two_minus_boundary(
    int level,
    int first_minus,
    int second_minus,
    const std::vector<int>& plus
) {
    std::unordered_map<std::uint64_t, cpp_int> coefficients;
    coefficients.emplace(pair_key(0, 0), 1);
    for (int label : plus) {
        std::unordered_map<std::uint64_t, cpp_int> next;
        for (const auto& [key, coefficient] : coefficients) {
            const auto [left, right] = decode_pair(key);
            for_each_fusion_output(level, left, label, [&](int output) {
                next[pair_key(output, right)] += coefficient;
            });
            for_each_fusion_output(level, right, label, [&](int output) {
                next[pair_key(left, output)] += coefficient;
            });
        }
        coefficients = std::move(next);
    }

    const auto coefficient_at = [&coefficients](int left, int right) {
        const auto found = coefficients.find(pair_key(left, right));
        return found == coefficients.end() ? cpp_int(0) : found->second;
    };

    cpp_int answer = -coefficient_at(first_minus, second_minus)
        - coefficient_at(second_minus, first_minus);
    for_each_fusion_output(level, first_minus, second_minus, [&](int output) {
        answer += coefficient_at(output, 0) + coefficient_at(0, output);
    });
    return answer;
}

cpp_int invariant_multiplicity(int level, const std::vector<int>& labels) {
    std::unordered_map<int, cpp_int> states;
    states.emplace(0, 1);
    for (int label : labels) {
        std::unordered_map<int, cpp_int> next;
        for (const auto& [source, coefficient] : states) {
            for_each_fusion_output(level, source, label, [&](int output) {
                next[output] += coefficient;
            });
        }
        states = std::move(next);
    }
    const auto found = states.find(0);
    return found == states.end() ? cpp_int(0) : found->second;
}

cpp_int q3_subset(
    int level,
    const std::vector<int>& minus,
    const std::vector<int>& plus
) {
    std::vector<int> labels = minus;
    labels.insert(labels.end(), plus.begin(), plus.end());
    if (labels.size() >= 63U) {
        throw std::runtime_error("subset cross-check requires fewer than 63 factors");
    }

    cpp_int answer = 0;
    const std::uint64_t subset_count = std::uint64_t{1} << labels.size();
    for (std::uint64_t mask = 0; mask < subset_count; ++mask) {
        std::vector<int> left;
        std::vector<int> right;
        int parity = 0;
        for (std::size_t index = 0; index < labels.size(); ++index) {
            if (((mask >> index) & 1U) != 0U) {
                left.push_back(labels[index]);
                if (index < minus.size()) {
                    parity ^= 1;
                }
            } else {
                right.push_back(labels[index]);
            }
        }
        const cpp_int term = invariant_multiplicity(level, left)
            * invariant_multiplicity(level, right);
        answer += parity == 0 ? term : -term;
    }
    return answer;
}

void combinations_rec(
    const std::vector<int>& labels,
    int remaining,
    std::size_t first,
    std::vector<int>& current,
    std::vector<std::vector<int>>& output
) {
    if (remaining == 0) {
        output.push_back(current);
        return;
    }
    for (std::size_t index = first; index < labels.size(); ++index) {
        current.push_back(labels[index]);
        combinations_rec(labels, remaining - 1, index, current, output);
        current.pop_back();
    }
}

std::vector<std::vector<int>> combinations(
    const std::vector<int>& labels,
    int size
) {
    std::vector<std::vector<int>> output;
    std::vector<int> current;
    combinations_rec(labels, size, 0U, current, output);
    return output;
}

bool earlier_witness(
    const std::vector<int>& minus,
    const std::vector<int>& plus,
    const Witness& witness
) {
    if (minus != witness.minus) {
        return std::lexicographical_compare(
            minus.begin(), minus.end(), witness.minus.begin(), witness.minus.end()
        );
    }
    return std::lexicographical_compare(
        plus.begin(), plus.end(), witness.plus.begin(), witness.plus.end()
    );
}

void consider(
    Witness& witness,
    const cpp_int& value,
    const std::vector<int>& minus,
    const std::vector<int>& plus
) {
    if (!witness.initialized || value < witness.value
        || (value == witness.value && earlier_witness(minus, plus, witness))) {
        witness.initialized = true;
        witness.value = value;
        witness.minus = minus;
        witness.plus = plus;
    }
}

void print_vector(const std::vector<int>& values) {
    std::cout << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << ']';
}

void cross_check(int maximum_level) {
    for (std::uint64_t sample = 0; sample < 128U; ++sample) {
        std::uint64_t state = splitmix64(sample ^ 20260721U);
        const int level = 2 + static_cast<int>(
            state % static_cast<std::uint64_t>(maximum_level - 1)
        );
        state = splitmix64(state);
        const int factors = 2 + static_cast<int>(state % 6U);
        state = splitmix64(state);
        const int minus_count = 2 * (1 + static_cast<int>(
            state % static_cast<std::uint64_t>(factors / 2)
        ));
        std::vector<int> minus(static_cast<std::size_t>(minus_count));
        std::vector<int> plus(static_cast<std::size_t>(factors - minus_count));
        for (int& label : minus) {
            state = splitmix64(state);
            label = 1 + static_cast<int>(
                state % static_cast<std::uint64_t>(level)
            );
        }
        for (int& label : plus) {
            state = splitmix64(state);
            label = 1 + static_cast<int>(
                state % static_cast<std::uint64_t>(level)
            );
        }
        std::sort(minus.begin(), minus.end());
        std::sort(plus.begin(), plus.end());
        if (q3_dynamic(level, minus, plus) != q3_subset(level, minus, plus)) {
            throw std::runtime_error("dynamic and subset calculations disagree");
        }
    }
    for (std::uint64_t sample = 0; sample < 128U; ++sample) {
        std::uint64_t state = splitmix64(sample ^ 0x6a09e667f3bcc909ULL);
        const int level = 2 + static_cast<int>(
            state % static_cast<std::uint64_t>(maximum_level - 1)
        );
        state = splitmix64(state);
        const int first_minus = 1 + static_cast<int>(
            state % static_cast<std::uint64_t>(level)
        );
        state = splitmix64(state);
        const int second_minus = 1 + static_cast<int>(
            state % static_cast<std::uint64_t>(level)
        );
        state = splitmix64(state);
        const int plus_count = static_cast<int>(state % 7U);
        std::vector<int> plus(static_cast<std::size_t>(plus_count));
        for (int& label : plus) {
            state = splitmix64(state);
            label = 1 + static_cast<int>(
                state % static_cast<std::uint64_t>(level)
            );
        }
        std::sort(plus.begin(), plus.end());
        const std::vector<int> minus{first_minus, second_minus};
        if (q3_dynamic(level, minus, plus)
            != q3_two_minus_boundary(
                level, first_minus, second_minus, plus)) {
            throw std::runtime_error(
                "two-minus boundary formula disagrees with dynamic calculation"
            );
        }
    }
    std::cout << "SU2_LEVEL_CENTRAL_Q3 independent_cross_checks=128"
              << " boundary_cross_checks=128 PASS\n";
}

Witness exhaustive_level(int level, int maximum_factors, std::uint64_t& cases) {
    std::vector<int> labels;
    for (int label = 1; label <= level; ++label) {
        labels.push_back(label);
    }

    Witness overall;
    cases = 0;
    for (int factors = 2; factors <= maximum_factors; ++factors) {
        for (int minus_count = 2; minus_count <= factors; minus_count += 2) {
            const auto minus_lists = combinations(labels, minus_count);
            const auto plus_lists = combinations(labels, factors - minus_count);
            cases += static_cast<std::uint64_t>(minus_lists.size())
                * static_cast<std::uint64_t>(plus_lists.size());

            #pragma omp parallel
            {
                Witness local;
                #pragma omp for schedule(dynamic)
                for (std::int64_t raw_index = 0;
                     raw_index < static_cast<std::int64_t>(minus_lists.size());
                     ++raw_index) {
                    const auto& minus = minus_lists[static_cast<std::size_t>(raw_index)];
                    for (const auto& plus : plus_lists) {
                        consider(local, q3_dynamic(level, minus, plus), minus, plus);
                    }
                }
                #pragma omp critical
                {
                    if (local.initialized) {
                        consider(
                            overall, local.value, local.minus, local.plus
                        );
                    }
                }
            }
        }
    }
    return overall;
}

Witness random_level(
    int level,
    int maximum_factors,
    std::uint64_t samples,
    std::uint64_t seed
) {
    Witness overall;
    #pragma omp parallel
    {
        Witness local;
        #pragma omp for schedule(dynamic, 16)
        for (std::int64_t raw_sample = 0;
             raw_sample < static_cast<std::int64_t>(samples);
             ++raw_sample) {
            const std::uint64_t sample = static_cast<std::uint64_t>(raw_sample);
            std::mt19937_64 generator(splitmix64(seed ^ sample));
            const int factors = 4 + static_cast<int>(
                generator() % static_cast<std::uint64_t>(maximum_factors - 3)
            );
            const int minus_count = 2 * (1 + static_cast<int>(
                generator() % static_cast<std::uint64_t>(factors / 2)
            ));
            std::vector<int> minus(static_cast<std::size_t>(minus_count));
            std::vector<int> plus(static_cast<std::size_t>(factors - minus_count));
            for (int& label : minus) {
                label = 1 + static_cast<int>(
                    generator() % static_cast<std::uint64_t>(level)
                );
            }
            for (int& label : plus) {
                label = 1 + static_cast<int>(
                    generator() % static_cast<std::uint64_t>(level)
                );
            }
            std::sort(minus.begin(), minus.end());
            std::sort(plus.begin(), plus.end());
            consider(local, q3_dynamic(level, minus, plus), minus, plus);
        }
        #pragma omp critical
        {
            if (local.initialized) {
                consider(overall, local.value, local.minus, local.plus);
            }
        }
    }
    return overall;
}

int parse_positive(const char* text, const char* name) {
    const long value = std::strtol(text, nullptr, 10);
    if (value <= 0 || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
}

std::uint64_t parse_u64(const char* text, const char* name) {
    const unsigned long long value = std::strtoull(text, nullptr, 10);
    if (value == 0U) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<std::uint64_t>(value);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout << std::unitbuf;
        if (argc == 5 && std::string(argv[1]) == "--random") {
            const int level = parse_positive(argv[2], "level");
            const int maximum_factors = parse_positive(
                argv[3], "random maximum factor count"
            );
            const std::uint64_t samples = parse_u64(argv[4], "random sample count");
            if (level < 2 || maximum_factors < 4) {
                throw std::runtime_error(
                    "random mode requires level >= 2 and maximum_factors >= 4"
                );
            }
            cross_check(level);
            const Witness witness = random_level(
                level, maximum_factors, samples, 20260721U
            );
            std::cout << "SU2_LEVEL_CENTRAL_Q3 random level=" << level
                      << " samples=" << samples
                      << " maximum_factors=" << maximum_factors
                      << " minimum=" << witness.value << " minus=";
            print_vector(witness.minus);
            std::cout << " plus=";
            print_vector(witness.plus);
            std::cout << '\n';
            const bool negative = witness.value < 0;
            std::cout << "SU2_LEVEL_CENTRAL_Q3 RANDOM SEARCH: "
                      << (negative
                          ? "NEGATIVE WITNESS FOUND"
                          : "NO NEGATIVE IN TESTED DOMAIN")
                      << '\n';
            return negative ? 1 : 0;
        }

        int maximum_level = 6;
        int maximum_factors = 8;
        if (argc != 1 && argc != 3) {
            throw std::runtime_error(
                "usage: search_su2_level_central_q3 "
                "[maximum_level maximum_factors] or "
                "--random level maximum_factors samples"
            );
        }
        if (argc == 3) {
            maximum_level = parse_positive(argv[1], "maximum level");
            maximum_factors = parse_positive(argv[2], "maximum factor count");
        }
        if (maximum_level < 2 || maximum_factors < 2) {
            throw std::runtime_error(
                "the maximum level and factor count must both be at least 2"
            );
        }

        cross_check(maximum_level);
        bool negative = false;
        for (int level = 2; level <= maximum_level; ++level) {
            std::uint64_t cases = 0;
            const Witness witness = exhaustive_level(level, maximum_factors, cases);
            std::cout << "SU2_LEVEL_CENTRAL_Q3 level=" << level
                      << " cases=" << cases
                      << " minimum=" << witness.value << " minus=";
            print_vector(witness.minus);
            std::cout << " plus=";
            print_vector(witness.plus);
            std::cout << '\n';
            if (witness.value < 0) {
                negative = true;
                break;
            }
        }
        std::cout << "SU2_LEVEL_CENTRAL_Q3 SEARCH: "
                  << (negative
                      ? "NEGATIVE WITNESS FOUND"
                      : "NO NEGATIVE IN TESTED DOMAIN")
                  << '\n';
        return negative ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_LEVEL_CENTRAL_Q3 FAILURE: " << error.what() << '\n';
        return 1;
    }
}
