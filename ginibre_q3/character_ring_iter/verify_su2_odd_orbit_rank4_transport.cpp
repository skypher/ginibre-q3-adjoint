#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

using boost::multiprecision::cpp_int;
using Rational = boost::rational<cpp_int>;

namespace {

constexpr int rank = 4;
constexpr int labels = rank - 1;
constexpr int pair_count = rank * (rank - 1) / 2;

struct Interval {
    Rational lower;
    Rational upper;
};

struct SpectralPair {
    int left = 0;
    int right = 0;
    int term_sign = 1;
    Interval initial_magnitude;
    std::array<Interval, labels> lambda{};
};

struct Chamber {
    int support_code = 0;
    int parity_code = 0;
    std::array<int, labels> signs{};
    std::array<bool, labels> active{};
    std::array<int, labels> minimum_power{};
    std::vector<SpectralPair> pairs;
    std::vector<int> negatives;
    std::vector<int> positives;
};

Rational rational(long long numerator, long long denominator = 1) {
    return Rational(cpp_int(numerator), cpp_int(denominator));
}

Interval point(long long value) {
    const Rational exact = rational(value);
    return Interval{exact, exact};
}

Interval interval(long long lower, long long upper, long long denominator) {
    return Interval{
        rational(lower, denominator),
        rational(upper, denominator)
    };
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
    const std::array<Rational, 4> products{
        left.lower * right.lower,
        left.lower * right.upper,
        left.upper * right.lower,
        left.upper * right.upper
    };
    const auto bounds = std::minmax_element(products.begin(), products.end());
    return Interval{*bounds.first, *bounds.second};
}

Interval scale(const Interval& value, long long factor) {
    return multiply(value, point(factor));
}

Interval power(Interval base, int exponent) {
    if (exponent < 0) {
        throw std::runtime_error("negative interval exponent");
    }
    Interval result = point(1);
    while (exponent != 0) {
        if ((exponent & 1) != 0) {
            result = multiply(result, base);
        }
        exponent /= 2;
        if (exponent != 0) {
            base = multiply(base, base);
        }
    }
    return result;
}

int strict_sign(const Interval& value) {
    if (value.lower > rational(0)) {
        return 1;
    }
    if (value.upper < rational(0)) {
        return -1;
    }
    throw std::runtime_error("interval does not determine a strict sign");
}

Interval absolute_value(const Interval& value) {
    const int sign = strict_sign(value);
    return sign == 1 ? value : negate(value);
}

bool proves_at_least(const Interval& left, const Interval& right) {
    return left.lower >= right.upper;
}

Rational cubic_value(const Rational& value) {
    return value * value * value - rational(3) * value + rational(1);
}

void verify_root_isolation(const std::array<Interval, rank>& nodes) {
    for (int node : {0, 1, 3}) {
        const Rational lower_image = cubic_value(
            nodes[static_cast<std::size_t>(node)].lower
        );
        const Rational upper_image = cubic_value(
            nodes[static_cast<std::size_t>(node)].upper
        );
        if (!((lower_image <= rational(0) && upper_image >= rational(0))
              || (lower_image >= rational(0)
                  && upper_image <= rational(0)))) {
            throw std::runtime_error("cubic root interval does not bracket");
        }
        const Interval derivative = subtract(
            scale(power(nodes[static_cast<std::size_t>(node)], 2), 3),
            point(3)
        );
        static_cast<void>(strict_sign(derivative));
    }
    if (!(nodes[0].lower > nodes[1].upper
          && nodes[1].lower > nodes[2].upper
          && nodes[2].lower > nodes[3].upper)) {
        throw std::runtime_error("spectral root order is not isolated");
    }
}

std::array<Interval, labels> orbit_values(const Interval& node) {
    const Interval square = power(node, 2);
    const Interval cube = multiply(square, node);
    return std::array<Interval, labels>{
        add(point(1), node),
        subtract(add(square, node), point(1)),
        subtract(add(add(cube, square), scale(node, -2)), point(1))
    };
}

std::vector<SpectralPair> make_pairs(
    const std::array<int, labels>& signs,
    const std::array<bool, labels>& active,
    const std::array<int, labels>& minimum_power,
    const std::array<Interval, rank>& weights,
    const std::array<std::array<Interval, labels>, rank>& values
) {
    std::vector<SpectralPair> result;
    for (int left = 0; left < rank; ++left) {
        for (int right = left + 1; right < rank; ++right) {
            SpectralPair pair;
            pair.left = left;
            pair.right = right;
            pair.initial_magnitude = scale(
                multiply(
                    weights[static_cast<std::size_t>(left)],
                    weights[static_cast<std::size_t>(right)]
                ),
                2
            );
            for (int label = 0; label < labels; ++label) {
                if (!active[static_cast<std::size_t>(label)]) {
                    pair.lambda[static_cast<std::size_t>(label)] = point(1);
                    continue;
                }
                const Interval signed_right
                    = signs[static_cast<std::size_t>(label)] == 1
                    ? values[static_cast<std::size_t>(right)]
                            [static_cast<std::size_t>(label)]
                    : negate(values[static_cast<std::size_t>(right)]
                                   [static_cast<std::size_t>(label)]);
                const Interval base = add(
                    values[static_cast<std::size_t>(left)]
                          [static_cast<std::size_t>(label)],
                    signed_right
                );
                const int base_sign = strict_sign(base);
                const int exponent
                    = minimum_power[static_cast<std::size_t>(label)];
                if (base_sign == -1 && (exponent & 1) != 0) {
                    pair.term_sign = -pair.term_sign;
                }
                const Interval magnitude = absolute_value(base);
                pair.initial_magnitude = multiply(
                    pair.initial_magnitude,
                    power(magnitude, exponent)
                );
                pair.lambda[static_cast<std::size_t>(label)]
                    = power(magnitude, 2);
            }
            result.push_back(pair);
        }
    }
    return result;
}

std::vector<cpp_int> apply_orbit_factor(
    const std::vector<cpp_int>& state,
    int radius,
    int sign
) {
    std::vector<cpp_int> next(state.size(), 0);
    for (int left = 0; left < rank; ++left) {
        for (int right = 0; right < rank; ++right) {
            const cpp_int& value = state[static_cast<std::size_t>(
                left * rank + right
            )];
            if (value == 0) {
                continue;
            }
            const auto apply_one = [&](int input, auto destination) {
                const int lower = std::abs(input - radius);
                const int upper = std::min(
                    input + radius,
                    2 * rank - 1 - input - radius
                );
                for (int output = lower; output <= upper; ++output) {
                    destination(output);
                }
            };
            apply_one(left, [&](int output) {
                next[static_cast<std::size_t>(output * rank + right)]
                    += value;
            });
            apply_one(right, [&](int output) {
                cpp_int& destination = next[static_cast<std::size_t>(
                    left * rank + output
                )];
                destination += sign == 1 ? value : -value;
            });
        }
    }
    return next;
}

cpp_int exact_corner(
    const std::array<int, labels>& signs,
    const std::array<int, labels>& powers
) {
    std::vector<cpp_int> state(
        static_cast<std::size_t>(rank * rank), 0
    );
    state[0] = 1;
    for (int label = 0; label < labels; ++label) {
        for (int copy = 0;
             copy < powers[static_cast<std::size_t>(label)];
             ++copy) {
            state = apply_orbit_factor(
                state,
                label + 1,
                signs[static_cast<std::size_t>(label)]
            );
        }
    }
    return state[0];
}

Interval floor_magnitude(
    const SpectralPair& pair,
    const std::array<int, labels>& floors
) {
    Interval result = pair.initial_magnitude;
    for (int label = 0; label < labels; ++label) {
        result = multiply(
            result,
            power(
                pair.lambda[static_cast<std::size_t>(label)],
                floors[static_cast<std::size_t>(label)]
            )
        );
    }
    return result;
}

Interval transformed_base(
    const SpectralPair& pair,
    const std::vector<int>& labels_to_multiply
) {
    Interval result = point(1);
    for (int label : labels_to_multiply) {
        result = multiply(
            result,
            pair.lambda[static_cast<std::size_t>(label)]
        );
    }
    return result;
}

bool hall_transport(
    const Chamber& chamber,
    const std::array<int, labels>& floors,
    const std::vector<std::vector<int>>& transformed_coordinates
) {
    std::vector<Interval> magnitudes;
    magnitudes.reserve(chamber.pairs.size());
    for (const SpectralPair& pair : chamber.pairs) {
        magnitudes.push_back(floor_magnitude(pair, floors));
    }

    std::vector<std::vector<int>> edges(chamber.negatives.size());
    for (std::size_t ni = 0; ni < chamber.negatives.size(); ++ni) {
        const SpectralPair& negative = chamber.pairs[
            static_cast<std::size_t>(chamber.negatives[ni])
        ];
        for (std::size_t pi_index = 0;
             pi_index < chamber.positives.size();
             ++pi_index) {
            const SpectralPair& positive = chamber.pairs[
                static_cast<std::size_t>(chamber.positives[pi_index])
            ];
            bool dominates = true;
            for (const std::vector<int>& coordinate
                 : transformed_coordinates) {
                if (!proves_at_least(
                        transformed_base(positive, coordinate),
                        transformed_base(negative, coordinate)
                    )) {
                    dominates = false;
                    break;
                }
            }
            if (dominates) {
                edges[ni].push_back(static_cast<int>(pi_index));
            }
        }
    }

    const int subsets = 1 << chamber.negatives.size();
    for (int subset = 1; subset < subsets; ++subset) {
        Interval demand = point(0);
        std::vector<bool> neighbors(chamber.positives.size(), false);
        for (std::size_t ni = 0; ni < chamber.negatives.size(); ++ni) {
            if (((subset >> ni) & 1) == 0) {
                continue;
            }
            demand = add(
                demand,
                magnitudes[static_cast<std::size_t>(
                    chamber.negatives[ni]
                )]
            );
            for (int pi_index : edges[ni]) {
                neighbors[static_cast<std::size_t>(pi_index)] = true;
            }
        }
        Interval capacity = point(0);
        for (std::size_t pi_index = 0;
             pi_index < chamber.positives.size();
             ++pi_index) {
            if (neighbors[pi_index]) {
                capacity = add(
                    capacity,
                    magnitudes[static_cast<std::size_t>(
                        chamber.positives[pi_index]
                    )]
                );
            }
        }
        if (!proves_at_least(capacity, demand)) {
            return false;
        }
    }
    return true;
}

int tail_shift(int support_code, int parity_code, int extra_mask) {
    if (extra_mask != 4) {
        return -1;
    }
    if (support_code == 20 && parity_code == 3) {
        return 2;
    }
    if ((support_code == 21 && parity_code == 1)
        || (support_code == 23
            && (parity_code == 2 || parity_code == 5))) {
        return 1;
    }
    return -1;
}

bool is_two_cone_case(int support_code, int parity_code, int extra_mask) {
    static_cast<void>(support_code);
    static_cast<void>(parity_code);
    return (extra_mask & 5) == 5;
}

int expected_common_shift(int support_code, int parity_code, int extra_mask) {
    if (extra_mask != 5) {
        return -1;
    }
    if (support_code == 20 && parity_code == 3) {
        return 4;
    }
    if (support_code == 23 && parity_code == 2) {
        return 3;
    }
    if (support_code == 23 && parity_code == 5) {
        return 2;
    }
    return -1;
}

}  // namespace

int main() {
    try {
        constexpr long long denominator = 1'000'000'000'000LL;
        const std::array<Interval, rank> nodes{
            interval(1'532'088'886'237LL, 1'532'088'886'238LL, denominator),
            interval(347'296'355'333LL, 347'296'355'334LL, denominator),
            point(-1),
            interval(-1'879'385'241'572LL, -1'879'385'241'571LL,
                     denominator)
        };
        verify_root_isolation(nodes);

        std::array<Interval, rank> weights{};
        std::array<std::array<Interval, labels>, rank> values{};
        for (int node = 0; node < rank; ++node) {
            weights[static_cast<std::size_t>(node)] = multiply(
                subtract(point(2), nodes[static_cast<std::size_t>(node)]),
                Interval{rational(1, 9), rational(1, 9)}
            );
            values[static_cast<std::size_t>(node)]
                = orbit_values(nodes[static_cast<std::size_t>(node)]);
        }

        std::uint64_t chambers = 0U;
        std::uint64_t pointwise = 0U;
        std::uint64_t base_points = 0U;
        std::uint64_t direct_regimes = 0U;
        std::uint64_t peeled_points = 0U;
        std::uint64_t tail_regimes = 0U;
        std::uint64_t two_cone_regimes = 0U;
        std::uint64_t shifted_two_cone_regimes = 0U;
        std::uint64_t strip_regimes = 0U;
        std::uint64_t strip_peeled_points = 0U;
        std::uint64_t unresolved_regimes = 0U;
        for (int support_code = 0; support_code < 27; ++support_code) {
            int code = support_code;
            std::array<int, labels> signs{};
            std::array<bool, labels> active{};
            int active_count = 0;
            int minus_labels = 0;
            for (int label = 0; label < labels; ++label) {
                const int digit = code % 3;
                code /= 3;
                if (digit != 0) {
                    active[static_cast<std::size_t>(label)] = true;
                    signs[static_cast<std::size_t>(label)]
                        = digit == 1 ? 1 : -1;
                    ++active_count;
                    if (digit == 2) {
                        ++minus_labels;
                    }
                }
            }
            if (active_count == 0 || minus_labels == 0) {
                continue;
            }
            for (int parity_code = 0;
                 parity_code < (1 << active_count);
                 ++parity_code) {
                std::array<int, labels> minimum_power{};
                int active_index = 0;
                int total_minus_parity = 0;
                for (int label = 0; label < labels; ++label) {
                    if (!active[static_cast<std::size_t>(label)]) {
                        continue;
                    }
                    const bool odd
                        = ((parity_code >> active_index) & 1) != 0;
                    minimum_power[static_cast<std::size_t>(label)]
                        = odd ? 1 : 2;
                    if (signs[static_cast<std::size_t>(label)] == -1
                        && odd) {
                        total_minus_parity ^= 1;
                    }
                    ++active_index;
                }
                if (total_minus_parity != 0) {
                    continue;
                }

                Chamber chamber;
                chamber.support_code = support_code;
                chamber.parity_code = parity_code;
                chamber.signs = signs;
                chamber.active = active;
                chamber.minimum_power = minimum_power;
                chamber.pairs = make_pairs(
                    signs, active, minimum_power, weights, values
                );
                for (int pair = 0; pair < pair_count; ++pair) {
                    if (chamber.pairs[static_cast<std::size_t>(pair)]
                            .term_sign < 0) {
                        chamber.negatives.push_back(pair);
                    } else {
                        chamber.positives.push_back(pair);
                    }
                }
                if (chamber.negatives.empty()) {
                    ++pointwise;
                    continue;
                }
                ++chambers;
                if (exact_corner(signs, minimum_power) < 0) {
                    throw std::runtime_error("negative chamber base point");
                }
                ++base_points;

                for (int extra_mask = 1;
                     extra_mask < (1 << labels);
                     ++extra_mask) {
                    bool compatible = true;
                    std::array<int, labels> floors{};
                    std::vector<std::vector<int>> coordinates;
                    for (int label = 0; label < labels; ++label) {
                        if (((extra_mask >> label) & 1) == 0) {
                            continue;
                        }
                        if (!active[static_cast<std::size_t>(label)]) {
                            compatible = false;
                        }
                        floors[static_cast<std::size_t>(label)] = 1;
                        coordinates.push_back(std::vector<int>{label});
                    }
                    if (!compatible) {
                        continue;
                    }
                    if (hall_transport(chamber, floors, coordinates)) {
                        ++direct_regimes;
                        continue;
                    }

                    const int shift = tail_shift(
                        support_code, parity_code, extra_mask
                    );
                    if (shift >= 0) {
                        for (int residual = 1;
                             residual <= shift;
                             ++residual) {
                            std::array<int, labels> powers = minimum_power;
                            powers[2] += 2 * residual;
                            if (exact_corner(signs, powers) < 0) {
                                throw std::runtime_error(
                                    "negative peeled tail point"
                                );
                            }
                            ++peeled_points;
                        }
                        floors[2] = 1 + shift;
                        if (!hall_transport(chamber, floors, coordinates)) {
                            throw std::runtime_error(
                                "exact tail transport failed"
                            );
                        }
                        ++tail_regimes;
                        continue;
                    }

                    if (is_two_cone_case(
                            support_code, parity_code, extra_mask
                        )) {
                        const std::vector<int> product_coordinate{0, 2};
                        bool cone_passed = true;
                        for (int preferred : {0, 2}) {
                            std::vector<std::vector<int>> cone{
                                product_coordinate,
                                std::vector<int>{preferred}
                            };
                            if ((extra_mask & 2) != 0) {
                                cone.push_back(std::vector<int>{1});
                            }
                            if (!hall_transport(chamber, floors, cone)) {
                                cone_passed = false;
                            }
                        }
                        if (cone_passed) {
                            ++two_cone_regimes;
                            continue;
                        }
                        int successful_common_shift = -1;
                        for (int common_shift = 1;
                             common_shift <= 6;
                             ++common_shift) {
                            std::array<int, labels> shifted_floors = floors;
                            shifted_floors[0] += common_shift;
                            shifted_floors[2] += common_shift;
                            bool shifted_passed = true;
                            for (int preferred : {0, 2}) {
                                std::vector<std::vector<int>> cone{
                                    product_coordinate,
                                    std::vector<int>{preferred}
                                };
                                if ((extra_mask & 2) != 0) {
                                    cone.push_back(std::vector<int>{1});
                                }
                                if (!hall_transport(
                                        chamber, shifted_floors, cone
                                    )) {
                                    shifted_passed = false;
                                }
                            }
                            if (shifted_passed) {
                                successful_common_shift = common_shift;
                                break;
                            }
                        }
                        if (successful_common_shift >= 0) {
                            if (successful_common_shift != expected_common_shift(
                                    support_code,
                                    parity_code,
                                    extra_mask
                                )) {
                                throw std::runtime_error(
                                    "unexpected shifted-cone classification"
                                );
                            }
                            bool strips_passed = true;
                            for (int fixed_residual = 1;
                                 fixed_residual <= successful_common_shift
                                    && strips_passed;
                                 ++fixed_residual) {
                                for (int fixed_label : {0, 2}) {
                                    const int variable_label
                                        = fixed_label == 0 ? 2 : 0;
                                    int variable_shift = -1;
                                    for (int candidate_shift = 0;
                                         candidate_shift <= 12;
                                         ++candidate_shift) {
                                        std::array<int, labels> strip_floors{};
                                        strip_floors[
                                            static_cast<std::size_t>(
                                                fixed_label
                                            )
                                        ] = fixed_residual;
                                        strip_floors[
                                            static_cast<std::size_t>(
                                                variable_label
                                            )
                                        ] = 1 + candidate_shift;
                                        if (hall_transport(
                                                chamber,
                                                strip_floors,
                                                std::vector<
                                                    std::vector<int>
                                                >{
                                                    std::vector<int>{
                                                        variable_label
                                                    }
                                                }
                                            )) {
                                            variable_shift = candidate_shift;
                                            break;
                                        }
                                    }
                                    if (variable_shift < 0) {
                                        strips_passed = false;
                                        break;
                                    }
                                    ++strip_regimes;
                                    for (int variable_residual = 1;
                                         variable_residual <= variable_shift;
                                         ++variable_residual) {
                                        std::array<int, labels> powers
                                            = minimum_power;
                                        powers[static_cast<std::size_t>(
                                            fixed_label
                                        )] += 2 * fixed_residual;
                                        powers[static_cast<std::size_t>(
                                            variable_label
                                        )] += 2 * variable_residual;
                                        if (exact_corner(signs, powers) < 0) {
                                            throw std::runtime_error(
                                                "negative strip point"
                                            );
                                        }
                                        ++strip_peeled_points;
                                    }
                                }
                            }
                            if (strips_passed) {
                                ++shifted_two_cone_regimes;
                                continue;
                            }
                        }
                    }
                    std::cerr << "unresolved exact regime support_code="
                              << support_code
                              << " parity_code=" << parity_code
                              << " extra_mask=" << extra_mask << '\n';
                    ++unresolved_regimes;
                }
            }
        }

        std::cout << "SU2_ODD_ORBIT_RANK4_TRANSPORT"
                  << " chambers=" << chambers
                  << " pointwise_chambers=" << pointwise
                  << " exact_base_points=" << base_points
                  << " direct_regimes=" << direct_regimes
                  << " exact_peeled_points=" << peeled_points
                  << " tail_regimes=" << tail_regimes
                  << " two_cone_regimes=" << two_cone_regimes
                  << " shifted_two_cone_regimes="
                  << shifted_two_cone_regimes
                  << " strip_regimes=" << strip_regimes
                  << " strip_peeled_points=" << strip_peeled_points
                  << " unresolved_regimes=" << unresolved_regimes
                  << " result="
                  << (unresolved_regimes == 0U ? "PASS" : "FAIL")
                  << '\n';
        return unresolved_regimes == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
