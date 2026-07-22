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

// Return the exact coefficient of 1 tensor 1 in
// product_i (chi_i tensor 1 + sign_i 1 tensor chi_i).
// Clebsch--Gordan is applied independently in the two tensor factors.
cpp_int q3_dynamic(
    const std::vector<int>& minus,
    const std::vector<int>& plus
) {
    std::unordered_map<std::uint64_t, cpp_int> states;
    states.emplace(pair_key(0, 0), 1);

    auto apply = [&states](int weight, int sign) {
        std::unordered_map<std::uint64_t, cpp_int> next;
        next.reserve(states.size() * static_cast<std::size_t>(weight + 2));
        for (const auto& [key, coefficient] : states) {
            const auto [left, right] = decode_pair(key);
            for (int out = std::abs(left - weight); out <= left + weight; out += 2) {
                next[pair_key(out, right)] += coefficient;
            }
            for (int out = std::abs(right - weight); out <= right + weight; out += 2) {
                next[pair_key(left, out)] += sign * coefficient;
            }
        }
        states = std::move(next);
    };

    for (int weight : minus) {
        apply(weight, -1);
    }
    for (int weight : plus) {
        apply(weight, 1);
    }
    const auto found = states.find(pair_key(0, 0));
    return found == states.end() ? cpp_int(0) : found->second;
}

cpp_int invariant_multiplicity(const std::vector<int>& weights) {
    std::unordered_map<int, cpp_int> current;
    current.emplace(0, 1);
    for (int weight : weights) {
        std::unordered_map<int, cpp_int> next;
        for (const auto& [source, multiplicity] : current) {
            for (int out = std::abs(source - weight); out <= source + weight; out += 2) {
                next[out] += multiplicity;
            }
        }
        current = std::move(next);
    }
    const auto found = current.find(0);
    return found == current.end() ? cpp_int(0) : found->second;
}

// Deliberately separate subset expansion used to check the dynamic program.
cpp_int q3_subset(
    const std::vector<int>& minus,
    const std::vector<int>& plus
) {
    std::vector<int> weights = minus;
    weights.insert(weights.end(), plus.begin(), plus.end());
    if (weights.size() >= 63U) {
        throw std::runtime_error("subset cross-check requires fewer than 63 factors");
    }
    cpp_int answer = 0;
    const std::uint64_t subsets = std::uint64_t{1} << weights.size();
    for (std::uint64_t mask = 0; mask < subsets; ++mask) {
        std::vector<int> left;
        std::vector<int> right;
        int minus_parity = 0;
        for (std::size_t index = 0; index < weights.size(); ++index) {
            if ((mask >> index) & 1U) {
                right.push_back(weights[index]);
                if (index < minus.size()) {
                    minus_parity ^= 1;
                }
            } else {
                left.push_back(weights[index]);
            }
        }
        cpp_int term = invariant_multiplicity(left) * invariant_multiplicity(right);
        answer += minus_parity == 0 ? term : -term;
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
    combinations_rec(labels, size, 0, current, output);
    return output;
}

void consider(
    Witness& best,
    const cpp_int& value,
    const std::vector<int>& minus,
    const std::vector<int>& plus
) {
    if (!best.initialized || value < best.value) {
        best.initialized = true;
        best.value = value;
        best.minus = minus;
        best.plus = plus;
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

void cross_check() {
    for (std::uint64_t sample = 0; sample < 256; ++sample) {
        std::mt19937_64 generator(splitmix64(sample));
        const int factors = 2 + static_cast<int>(generator() % 7U);
        int minus_count = 2 + 2 * static_cast<int>(
            generator() % static_cast<std::uint64_t>(factors / 2)
        );
        if (minus_count > factors) {
            minus_count -= 2;
        }
        std::vector<int> minus(static_cast<std::size_t>(minus_count));
        std::vector<int> plus(static_cast<std::size_t>(factors - minus_count));
        for (int& weight : minus) {
            weight = 1 + static_cast<int>(generator() % 8U);
        }
        for (int& weight : plus) {
            weight = 1 + static_cast<int>(generator() % 8U);
        }
        std::sort(minus.begin(), minus.end());
        std::sort(plus.begin(), plus.end());
        const cpp_int dynamic = q3_dynamic(minus, plus);
        const cpp_int subset = q3_subset(minus, plus);
        if (dynamic != subset) {
            throw std::runtime_error("dynamic and subset calculations disagree");
        }
    }
    std::cout << "SU2_CENTRAL_Q3 independent_cross_checks=256 PASS\n";
}

Witness exhaustive_search(
    const std::string& name,
    const std::vector<int>& labels,
    int maximum_factors
) {
    Witness overall;
    for (int factors = 2; factors <= maximum_factors; ++factors) {
        Witness degree_best;
        std::uint64_t degree_cases = 0;
        for (int minus_count = 2; minus_count <= factors; minus_count += 2) {
            const auto minus_lists = combinations(labels, minus_count);
            const auto plus_lists = combinations(labels, factors - minus_count);
            const std::uint64_t block_cases =
                static_cast<std::uint64_t>(minus_lists.size())
                * static_cast<std::uint64_t>(plus_lists.size());
            degree_cases += block_cases;

            #pragma omp parallel
            {
                Witness local;
                #pragma omp for schedule(dynamic)
                for (std::int64_t i = 0;
                     i < static_cast<std::int64_t>(minus_lists.size()); ++i) {
                    for (const auto& plus : plus_lists) {
                        const cpp_int value = q3_dynamic(
                            minus_lists[static_cast<std::size_t>(i)], plus
                        );
                        consider(
                            local, value,
                            minus_lists[static_cast<std::size_t>(i)], plus
                        );
                    }
                }
                #pragma omp critical
                {
                    if (local.initialized) {
                        consider(
                            degree_best, local.value, local.minus, local.plus
                        );
                    }
                }
            }
        }
        consider(
            overall, degree_best.value, degree_best.minus, degree_best.plus
        );
        std::cout << "SU2_CENTRAL_Q3 exhaustive group=" << name
                  << " factors=" << factors
                  << " cases=" << degree_cases
                  << " minimum=" << degree_best.value << " minus=";
        print_vector(degree_best.minus);
        std::cout << " plus=";
        print_vector(degree_best.plus);
        std::cout << '\n';
        if (degree_best.value < 0) {
            break;
        }
    }
    return overall;
}

Witness random_search(
    const std::string& name,
    int label_step,
    int maximum_label_index,
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
             raw_sample < static_cast<std::int64_t>(samples); ++raw_sample) {
            const std::uint64_t sample = static_cast<std::uint64_t>(raw_sample);
            std::mt19937_64 generator(splitmix64(seed ^ sample));
            const int factors = 4 + static_cast<int>(
                generator() % static_cast<std::uint64_t>(maximum_factors - 3)
            );
            const int even_choices = factors / 2;
            const int minus_count = 2 * (
                1 + static_cast<int>(generator() % static_cast<std::uint64_t>(even_choices))
            );
            std::vector<int> minus(static_cast<std::size_t>(minus_count));
            std::vector<int> plus(static_cast<std::size_t>(factors - minus_count));
            // A product of SU(2) irreducibles can contain the trivial module
            // only if its total highest weight is even and its largest weight
            // is no greater than the sum of all the others.  Rejection here
            // removes identically zero Q3 cases without imposing a sign bias.
            while (true) {
                int total_weight = 0;
                int largest_weight = 0;
                for (int& weight : minus) {
                    weight = label_step * (
                        1 + static_cast<int>(generator()
                            % static_cast<std::uint64_t>(maximum_label_index))
                    );
                    total_weight += weight;
                    largest_weight = std::max(largest_weight, weight);
                }
                for (int& weight : plus) {
                    weight = label_step * (
                        1 + static_cast<int>(generator()
                            % static_cast<std::uint64_t>(maximum_label_index))
                    );
                    total_weight += weight;
                    largest_weight = std::max(largest_weight, weight);
                }
                if ((total_weight & 1) == 0 && 2 * largest_weight <= total_weight) {
                    break;
                }
            }
            std::sort(minus.begin(), minus.end());
            std::sort(plus.begin(), plus.end());
            consider(local, q3_dynamic(minus, plus), minus, plus);
        }
        #pragma omp critical
        {
            if (local.initialized) {
                consider(overall, local.value, local.minus, local.plus);
            }
        }
    }
    std::cout << "SU2_CENTRAL_Q3 random group=" << name
              << " samples=" << samples << " maximum_factors=" << maximum_factors
              << " maximum_weight=" << label_step * maximum_label_index
              << " seed=" << seed << " minimum=" << overall.value << " minus=";
    print_vector(overall.minus);
    std::cout << " plus=";
    print_vector(overall.plus);
    std::cout << '\n';
    return overall;
}

int parse_int(const char* text, const char* name) {
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
        int exhaustive_weight_count = 6;
        int exhaustive_factors = 8;
        int random_weight_count = 16;
        int random_factors = 14;
        std::uint64_t random_samples = 20000;
        std::uint64_t seed = 20260721;
        if (argc != 1 && argc != 6) {
            throw std::runtime_error(
                "usage: search_su2_central_q3 "
                "[exhaustive_weight_count exhaustive_factors "
                "random_weight_count random_factors random_samples]"
            );
        }
        if (argc == 6) {
            exhaustive_weight_count = parse_int(argv[1], "exhaustive weight count");
            exhaustive_factors = parse_int(argv[2], "exhaustive factor count");
            random_weight_count = parse_int(argv[3], "random weight count");
            random_factors = parse_int(argv[4], "random factor count");
            random_samples = parse_u64(argv[5], "random sample count");
        }
        if (exhaustive_factors < 2) {
            throw std::runtime_error("exhaustive factor count must be at least 2");
        }
        if (random_factors < 4) {
            throw std::runtime_error("random factor count must be at least 4");
        }

        cross_check();
        std::vector<int> su2_labels;
        std::vector<int> so3_labels;
        for (int index = 1; index <= exhaustive_weight_count; ++index) {
            su2_labels.push_back(index);
            so3_labels.push_back(2 * index);
        }

        Witness su2 = exhaustive_search("SU2", su2_labels, exhaustive_factors);
        Witness so3 = exhaustive_search("SO3", so3_labels, exhaustive_factors);
        Witness su2_random = random_search(
            "SU2", 1, random_weight_count, random_factors,
            random_samples, seed
        );
        Witness so3_random = random_search(
            "SO3", 2, random_weight_count, random_factors,
            random_samples, seed ^ 0x5a5a5a5a5a5a5a5aULL
        );

        const bool negative = su2.value < 0 || so3.value < 0
            || su2_random.value < 0 || so3_random.value < 0;
        std::cout << "SU2_CENTRAL_Q3 SEARCH: "
                  << (negative ? "NEGATIVE WITNESS FOUND" : "NO NEGATIVE IN TESTED DOMAIN")
                  << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SU2_CENTRAL_Q3 FAILURE: " << error.what() << '\n';
        return 1;
    }
}
