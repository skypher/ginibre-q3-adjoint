#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

using boost::multiprecision::cpp_int;
using Rational = boost::rational<cpp_int>;

namespace {

struct Interval {
    Rational lower;
    Rational upper;
};

struct PairData {
    int first = 0;
    int second = 0;
    Interval difference_squared;
    Interval endpoint_minus;
    Interval endpoint_plus;
    Interval endpoint_product;
    Interval kernel;
};

int parse_positive(const char* text, const char* name) {
    const int value = std::stoi(text);
    if (value <= 0) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return value;
}

Rational minimum(
    const Rational& first,
    const Rational& second
) {
    return first < second ? first : second;
}

Rational maximum(
    const Rational& first,
    const Rational& second
) {
    return first > second ? first : second;
}

Interval add(const Interval& left, const Interval& right) {
    return Interval{
        left.lower + right.lower,
        left.upper + right.upper
    };
}

Interval negate(const Interval& value) {
    return Interval{-value.upper, -value.lower};
}

Interval subtract(const Interval& left, const Interval& right) {
    return add(left, negate(right));
}

Interval multiply(const Interval& left, const Interval& right) {
    const Rational products[4]{
        left.lower * right.lower,
        left.lower * right.upper,
        left.upper * right.lower,
        left.upper * right.upper
    };
    Rational lower = products[0];
    Rational upper = products[0];
    for (int index = 1; index < 4; ++index) {
        lower = minimum(lower, products[index]);
        upper = maximum(upper, products[index]);
    }
    return Interval{lower, upper};
}

Interval square_nonnegative(const Interval& value) {
    if (value.lower < 0) {
        throw std::runtime_error(
            "square_nonnegative received a negative interval"
        );
    }
    return Interval{
        value.lower * value.lower,
        value.upper * value.upper
    };
}

Interval exact(long long value) {
    return Interval{Rational(value), Rational(value)};
}

Interval cos_at_rational(const Rational& argument) {
    const Rational argument_squared = argument * argument;
    Rational term(1);
    Rational sum(1);
    constexpr int maximum_term = 16;
    for (int index = 1; index <= maximum_term; ++index) {
        term *= argument_squared;
        term /= Rational(
            static_cast<long long>(2 * index - 1)
            * static_cast<long long>(2 * index)
        );
        if ((index & 1) == 0) {
            sum += term;
        } else {
            sum -= term;
        }
    }
    Rational next = term * argument_squared;
    next /= Rational(
        static_cast<long long>(2 * maximum_term + 1)
        * static_cast<long long>(2 * maximum_term + 2)
    );
    if (((maximum_term + 1) & 1) != 0) {
        next = -next;
    }
    return Interval{minimum(sum, sum + next), maximum(sum, sum + next)};
}

std::pair<Rational, Rational> pi_bounds() {
    const cpp_int denominator(
        "100000000000000000000"
    );
    const cpp_int lower_numerator(
        "314159265358979323846"
    );
    return {
        Rational(lower_numerator, denominator),
        Rational(lower_numerator + 1, denominator)
    };
}

Interval twice_cosine(int index, int order) {
    if (index < 0 || index > order) {
        throw std::runtime_error("cosine index outside [0,N]");
    }
    if (index == 0) {
        return exact(2);
    }
    if (index == order) {
        return exact(-2);
    }
    if (2 * index == order) {
        return exact(0);
    }

    bool negative = false;
    int reduced_index = index;
    if (2 * index > order) {
        negative = true;
        reduced_index = order - index;
    }
    const auto [pi_lower, pi_upper] = pi_bounds();
    const Rational lower_argument
        = pi_lower * Rational(reduced_index, order);
    const Rational upper_argument
        = pi_upper * Rational(reduced_index, order);
    const Interval lower_cosine = cos_at_rational(upper_argument);
    const Interval upper_cosine = cos_at_rational(lower_argument);
    Interval result{
        2 * lower_cosine.lower,
        2 * upper_cosine.upper
    };
    if (negative) {
        result = negate(result);
    }
    return result;
}

PairData make_pair_data(
    int first,
    int second,
    const std::vector<Interval>& values
) {
    const Interval two = exact(2);
    const Interval difference = subtract(
        values[static_cast<std::size_t>(first)],
        values[static_cast<std::size_t>(second)]
    );
    if (difference.lower <= 0) {
        throw std::runtime_error("spectral values are not strictly ordered");
    }
    const Interval endpoint_minus = multiply(
        subtract(two, values[static_cast<std::size_t>(first)]),
        subtract(two, values[static_cast<std::size_t>(second)])
    );
    const Interval endpoint_plus = multiply(
        add(two, values[static_cast<std::size_t>(first)]),
        add(two, values[static_cast<std::size_t>(second)])
    );
    return PairData{
        first,
        second,
        square_nonnegative(difference),
        endpoint_minus,
        endpoint_plus,
        multiply(endpoint_minus, endpoint_plus),
        add(
            two,
            multiply(
                values[static_cast<std::size_t>(first)],
                values[static_cast<std::size_t>(second)]
            )
        )
    };
}

bool at_least(const Interval& left, const Interval& right) {
    return left.lower >= right.upper;
}

bool certified_edge(
    const PairData& load,
    const PairData& credit,
    bool minus_dominant
) {
    const Interval load_magnitude = negate(load.kernel);
    if (!at_least(credit.kernel, load_magnitude)) {
        return false;
    }
    if (!at_least(credit.endpoint_product, load.endpoint_product)) {
        return false;
    }
    if (minus_dominant) {
        if (!at_least(credit.endpoint_minus, load.endpoint_minus)) {
            return false;
        }
    } else if (!at_least(credit.endpoint_plus, load.endpoint_plus)) {
        return false;
    }

    const Interval credit_side = multiply(
        multiply(
            credit.difference_squared,
            credit.endpoint_product
        ),
        credit.kernel
    );
    const Interval load_side = multiply(
        multiply(
            load.difference_squared,
            load.endpoint_product
        ),
        load_magnitude
    );
    return at_least(credit_side, load_side);
}

bool augment(
    int load,
    const std::vector<std::vector<int>>& adjacency,
    std::vector<int>& credit_match,
    std::vector<bool>& seen
) {
    for (int credit : adjacency[static_cast<std::size_t>(load)]) {
        if (seen[static_cast<std::size_t>(credit)]) {
            continue;
        }
        seen[static_cast<std::size_t>(credit)] = true;
        if (credit_match[static_cast<std::size_t>(credit)] < 0
            || augment(
                credit_match[static_cast<std::size_t>(credit)],
                adjacency,
                credit_match,
                seen
            )) {
            credit_match[static_cast<std::size_t>(credit)] = load;
            return true;
        }
    }
    return false;
}

int matching_size(
    const std::vector<std::vector<int>>& adjacency,
    int credit_count
) {
    std::vector<int> credit_match(
        static_cast<std::size_t>(credit_count),
        -1
    );
    int result = 0;
    for (int load = 0;
         load < static_cast<int>(adjacency.size());
         ++load) {
        std::vector<bool> seen(
            static_cast<std::size_t>(credit_count),
            false
        );
        if (augment(load, adjacency, credit_match, seen)) {
            ++result;
        }
    }
    return result;
}

bool neighbor_sets_are_nested(
    const std::vector<std::vector<int>>& adjacency
) {
    for (std::size_t first = 0; first < adjacency.size(); ++first) {
        for (std::size_t second = first + 1;
             second < adjacency.size();
             ++second) {
            const auto subset = [&](std::size_t left, std::size_t right) {
                for (int value : adjacency[left]) {
                    if (!std::binary_search(
                            adjacency[right].begin(),
                            adjacency[right].end(),
                            value
                        )) {
                        return false;
                    }
                }
                return true;
            };
            if (!subset(first, second) && !subset(second, first)) {
                return false;
            }
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error(
                "usage: analyze_su2_cyclic_jacobi_matching MAX_LEVEL"
            );
        }
        const int maximum_level = parse_positive(
            argv[1], "maximum level"
        );
        std::uint64_t graphs = 0U;
        std::uint64_t loads = 0U;
        std::uint64_t edges = 0U;
        std::uint64_t non_ferrers_graphs = 0U;
        int first_non_ferrers_level = -1;
        int first_matching_failure_level = -1;
        int maximum_loads = 0;

        for (int level = 2; level <= maximum_level; ++level) {
            const int order = level + 2;
            for (int sector_parity = 0;
                 sector_parity <= 1;
                 ++sector_parity) {
                std::vector<Interval> values;
                for (int index = sector_parity;
                     index <= order;
                     index += 2) {
                    if (index == 0 || index == order) {
                        continue;
                    }
                    values.push_back(twice_cosine(index, order));
                }

                std::vector<PairData> negative_pairs;
                std::vector<PairData> positive_pairs;
                for (int first = 0;
                     first < static_cast<int>(values.size());
                     ++first) {
                    for (int second = first + 1;
                         second < static_cast<int>(values.size());
                         ++second) {
                        const PairData pair = make_pair_data(
                            first, second, values
                        );
                        if (pair.kernel.upper < 0) {
                            negative_pairs.push_back(pair);
                        } else if (pair.kernel.lower > 0) {
                            positive_pairs.push_back(pair);
                        }
                    }
                }
                maximum_loads = std::max(
                    maximum_loads,
                    static_cast<int>(negative_pairs.size())
                );

                for (bool minus_dominant : {true, false}) {
                    std::vector<std::vector<int>> adjacency(
                        negative_pairs.size()
                    );
                    for (std::size_t load = 0;
                         load < negative_pairs.size();
                         ++load) {
                        for (std::size_t credit = 0;
                             credit < positive_pairs.size();
                             ++credit) {
                            if (certified_edge(
                                    negative_pairs[load],
                                    positive_pairs[credit],
                                    minus_dominant
                                )) {
                                adjacency[load].push_back(
                                    static_cast<int>(credit)
                                );
                                ++edges;
                            }
                        }
                    }
                    const int matched = matching_size(
                        adjacency,
                        static_cast<int>(positive_pairs.size())
                    );
                    if (matched
                        != static_cast<int>(negative_pairs.size())
                        && first_matching_failure_level < 0) {
                        first_matching_failure_level = level;
                        std::cout
                            << "FIRST_MATCHING_FAILURE"
                            << " level=" << level
                            << " order=" << order
                            << " sector="
                            << (sector_parity == 0 ? "plus" : "minus")
                            << " slope="
                            << (minus_dominant ? "m>=n" : "n>=m")
                            << " loads=" << negative_pairs.size()
                            << " matched=" << matched << '\n';
                    }
                    if (!neighbor_sets_are_nested(adjacency)) {
                        ++non_ferrers_graphs;
                        if (first_non_ferrers_level < 0) {
                            first_non_ferrers_level = level;
                        }
                    }
                    ++graphs;
                    loads += negative_pairs.size();
                }
            }
        }

        std::cout
            << "CYCLIC_JACOBI_MATCHING"
            << " max_level=" << maximum_level
            << " graphs=" << graphs
            << " loads=" << loads
            << " edges=" << edges
            << " maximum_loads=" << maximum_loads
            << " first_matching_failure_level="
            << first_matching_failure_level
            << " non_ferrers_graphs=" << non_ferrers_graphs
            << " first_non_ferrers_level=" << first_non_ferrers_level
            << '\n';
        return first_matching_failure_level < 0
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
