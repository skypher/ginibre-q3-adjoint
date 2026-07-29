#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;
using Exponent = std::array<int, 2>;

class Polynomial {
public:
    Polynomial() = default;

    explicit Polynomial(const Rational& constant) {
        if (constant != 0) {
            coefficients_[Exponent{0, 0}] = constant;
        }
    }

    static Polynomial variable(int index) {
        if (index < 0 || index >= 2) {
            throw std::runtime_error("invalid variable index");
        }
        Polynomial result;
        Exponent exponent{0, 0};
        exponent[static_cast<std::size_t>(index)] = 1;
        result.coefficients_[exponent] = 1;
        return result;
    }

    Polynomial& operator+=(const Polynomial& other) {
        for (const auto& [exponent, coefficient] : other.coefficients_) {
            coefficients_[exponent] += coefficient;
            if (coefficients_[exponent] == 0) {
                coefficients_.erase(exponent);
            }
        }
        return *this;
    }

    Polynomial& operator-=(const Polynomial& other) {
        return *this += Rational{-1} * other;
    }

    const std::map<Exponent, Rational>& coefficients() const {
        return coefficients_;
    }

private:
    std::map<Exponent, Rational> coefficients_;

    friend Polynomial operator*(
        const Polynomial& left,
        const Polynomial& right
    );
    friend Polynomial operator*(const Rational& scalar, Polynomial value);
};

Polynomial operator+(Polynomial left, const Polynomial& right) {
    left += right;
    return left;
}

Polynomial operator-(Polynomial left, const Polynomial& right) {
    left -= right;
    return left;
}

Polynomial operator*(
    const Polynomial& left,
    const Polynomial& right
) {
    Polynomial result;
    for (const auto& [left_exponent, left_coefficient] :
         left.coefficients_) {
        for (const auto& [right_exponent, right_coefficient] :
             right.coefficients_) {
            const Exponent exponent{
                left_exponent[0] + right_exponent[0],
                left_exponent[1] + right_exponent[1]
            };
            result.coefficients_[exponent] +=
                left_coefficient * right_coefficient;
            if (result.coefficients_[exponent] == 0) {
                result.coefficients_.erase(exponent);
            }
        }
    }
    return result;
}

Polynomial operator*(const Rational& scalar, Polynomial value) {
    if (scalar == 0) {
        return Polynomial{};
    }
    for (auto& [exponent, coefficient] : value.coefficients_) {
        static_cast<void>(exponent);
        coefficient *= scalar;
    }
    return value;
}

Polynomial operator*(Polynomial value, const Rational& scalar) {
    return scalar * value;
}

Polynomial choose_three(const Polynomial& value) {
    return (
        (value + Polynomial{Rational{1}})
        * (value + Polynomial{Rational{2}})
        * (value + Polynomial{Rational{3}})
    ) * Rational{1, 6};
}

Polynomial multiplicity(
    const Polynomial& q_half,
    const Polynomial& distance,
    bool second_cap
) {
    const Polynomial first =
        choose_three(Rational{5} * q_half - distance);
    const Polynomial second =
        choose_three(Rational{3} * q_half - distance
                     - Polynomial{Rational{1}});
    Polynomial result = first - Rational{5} * second;
    if (second_cap) {
        result += Rational{10}
            * choose_three(
                q_half - distance - Polynomial{Rational{2}}
            );
    }
    return result;
}

Polynomial margin(
    const Polynomial& q_half,
    const Polynomial& distance,
    bool second_cap
) {
    const Polynomial central = multiplicity(
        q_half,
        Polynomial{Rational{0}},
        true
    );
    return (
        Rational{2} * distance + Polynomial{Rational{1}}
    ) * central
        - multiplicity(q_half, distance, second_cap);
}

bool report(const std::string& name, const Polynomial& polynomial) {
    std::size_t negative = 0;
    for (const auto& [exponent, coefficient] :
         polynomial.coefficients()) {
        if (coefficient < 0) {
            ++negative;
        }
    }
    std::cout
        << "SU2_TERMINAL_ODD_LOW case=" << name
        << " terms=" << polynomial.coefficients().size()
        << " negative_coefficients=" << negative
        << " coefficients={";
    bool first = true;
    for (const auto& [exponent, coefficient] :
         polynomial.coefficients()) {
        if (!first) {
            std::cout << ',';
        }
        first = false;
        std::cout
            << '(' << exponent[0] << ',' << exponent[1]
            << "):" << coefficient;
    }
    std::cout << "}\n";
    return negative == 0;
}

Polynomial choose_five(const Polynomial& value) {
    Polynomial result{Rational{1}};
    for (int shift = 1; shift <= 5; ++shift) {
        result = result
            * (value + Polynomial{Rational{shift}});
    }
    return result * Rational{1, 120};
}

Polynomial multiplicity_d7(
    const Polynomial& q_half,
    const Polynomial& distance,
    bool final_cap
) {
    Polynomial result =
        choose_five(Rational{7} * q_half - distance)
        - Rational{7}
            * choose_five(
                Rational{5} * q_half - distance
                - Polynomial{Rational{1}}
            )
        + Rational{21}
            * choose_five(
                Rational{3} * q_half - distance
                - Polynomial{Rational{2}}
            );
    if (final_cap) {
        result -= Rational{35}
            * choose_five(
                q_half - distance - Polynomial{Rational{3}}
            );
    }
    return result;
}

Polynomial margin_d7(
    const Polynomial& q_half,
    const Polynomial& distance,
    bool final_cap
) {
    const Polynomial central = multiplicity_d7(
        q_half,
        Polynomial{Rational{0}},
        true
    );
    return (
        Rational{2} * distance + Polynomial{Rational{1}}
    ) * central
        - multiplicity_d7(q_half, distance, final_cap);
}

Integer binomial_integer(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    k = std::min(k, n - k);
    Integer result = 1;
    for (int index = 1; index <= k; ++index) {
        result *= n - k + index;
        result /= index;
    }
    return result;
}

Polynomial choose_degree(const Polynomial& value, int degree) {
    Polynomial result{Rational{1}};
    Integer factorial = 1;
    for (int shift = 1; shift <= degree; ++shift) {
        result = result
            * (value + Polynomial{Rational{shift}});
        factorial *= shift;
    }
    return result * Rational{1, factorial};
}

Polynomial multiplicity_odd(
    int half_distance,
    const Polynomial& q_half,
    const Polynomial& distance,
    bool final_cap
) {
    const int power = 2 * half_distance + 1;
    const int degree = power - 2;
    Polynomial result;
    const int last = final_cap
        ? half_distance
        : half_distance - 1;
    for (int image = 0; image <= last; ++image) {
        Rational coefficient{
            binomial_integer(power, image)
        };
        if ((image & 1) != 0) {
            coefficient *= -1;
        }
        result += coefficient * choose_degree(
            Rational{power - 2 * image} * q_half
                - distance - Polynomial{Rational{image}},
            degree
        );
    }
    return result;
}

Polynomial margin_odd(
    int half_distance,
    const Polynomial& q_half,
    const Polynomial& distance,
    bool final_cap
) {
    const Polynomial central = multiplicity_odd(
        half_distance,
        q_half,
        Polynomial{Rational{0}},
        true
    );
    return (
        Rational{2} * distance + Polynomial{Rational{1}}
    ) * central
        - multiplicity_odd(
            half_distance,
            q_half,
            distance,
            final_cap
        );
}

int parse_maximum(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value < 2
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error(
            "maximum half-distance must be an integer at least two"
        );
    }
    return static_cast<int>(value);
}

std::vector<Integer> stable_profile(int q_half, int power) {
    std::vector<Integer> profile(1, Integer{1});
    for (int step = 0; step < power; ++step) {
        std::vector<Integer> next(
            static_cast<std::size_t>((step + 1) * q_half + 1),
            Integer{0}
        );
        for (int source = 0; source <= step * q_half; ++source) {
            const Integer& coefficient =
                profile[static_cast<std::size_t>(source)];
            if (coefficient == 0) {
                continue;
            }
            for (int output = std::abs(source - q_half);
                 output <= source + q_half; ++output) {
                next[static_cast<std::size_t>(output)] += coefficient;
            }
        }
        profile = std::move(next);
    }
    return profile;
}

std::vector<Integer> multiply_profile(
    const std::vector<Integer>& profile,
    int label_half
) {
    const int maximum_source =
        static_cast<int>(profile.size()) - 1;
    std::vector<Integer> result(
        static_cast<std::size_t>(maximum_source + label_half + 1),
        Integer{0}
    );
    for (int source = 0; source <= maximum_source; ++source) {
        const Integer& coefficient =
            profile[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        for (int output = std::abs(source - label_half);
             output <= source + label_half; ++output) {
            result[static_cast<std::size_t>(output)] += coefficient;
        }
    }
    return result;
}

int analyze_distance_five_bootstrap(int maximum_q_half) {
    std::size_t cases = 0;
    std::size_t coefficient_rows = 0;
    std::size_t negative_coefficients = 0;
    int first_q_half = 0;
    int first_distance_half = 0;
    int first_output_half = 0;
    Integer first_value = 0;

    for (int q_half = 1;
         q_half <= maximum_q_half; ++q_half) {
        const std::vector<Integer> fifth =
            stable_profile(q_half, 5);
        for (int distance_half = 1;
             distance_half <= 2 * q_half; ++distance_half) {
            const std::vector<Integer> upper =
                multiply_profile(fifth, distance_half);
            const std::vector<Integer> lower =
                multiply_profile(fifth, distance_half - 1);
            ++cases;
            for (std::size_t output = 0;
                 output < upper.size(); ++output) {
                Integer coefficient =
                    Integer{2}
                        * (
                            output < fifth.size()
                                ? fifth[output]
                                : Integer{0}
                        )
                    - upper[output]
                    + (
                        output < lower.size()
                            ? lower[output]
                            : Integer{0}
                    );
                ++coefficient_rows;
                if (coefficient < 0) {
                    ++negative_coefficients;
                    if (first_q_half == 0) {
                        first_q_half = q_half;
                        first_distance_half = distance_half;
                        first_output_half =
                            static_cast<int>(output);
                        first_value = coefficient;
                    }
                }
            }
        }
    }

    std::cout
        << "SU2_TERMINAL_DISTANCE_FIVE_BOOTSTRAP"
        << " maximum_q_half=" << maximum_q_half
        << " cases=" << cases
        << " coefficient_rows=" << coefficient_rows
        << " negative_coefficients=" << negative_coefficients;
    if (first_q_half != 0) {
        std::cout
            << " first_negative={q_half=" << first_q_half
            << " distance_half=" << first_distance_half
            << " output_half=" << first_output_half
            << " coefficient=" << first_value << '}';
    }
    std::cout
        << " result="
        << (
            negative_coefficients == 0
                ? "PASS_BOUNDED_DIAGNOSTIC"
                : "COUNTEREXAMPLE"
        )
        << '\n';
    return negative_coefficients == 0
        ? EXIT_SUCCESS
        : EXIT_FAILURE;
}

int analyze_profile_concavity(int maximum_half_power, int maximum_q_half) {
    std::size_t profiles = 0;
    std::size_t concavity_rows = 0;
    std::size_t negative_concavity_rows = 0;
    std::size_t low_concavity_rows = 0;
    std::size_t negative_low_concavity_rows = 0;
    std::size_t paired_slope_rows = 0;
    std::size_t negative_paired_slope_rows = 0;
    std::size_t transfer_identity_rows = 0;
    std::size_t failed_transfer_identity_rows = 0;
    std::size_t negative_transfer_tail_rows = 0;
    std::size_t block_identity_rows = 0;
    std::size_t failed_block_identity_rows = 0;
    std::size_t negative_block_payment_rows = 0;
    std::size_t dimension_rows = 0;
    std::size_t failed_dimension_rows = 0;
    std::size_t odlp_rows = 0;
    std::size_t failed_odlp_rows = 0;
    int first_concavity_half_power = 0;
    int first_concavity_q_half = 0;
    int first_concavity_label_half = 0;
    Integer first_concavity_value = 0;

    for (int half_power = 2;
         half_power <= maximum_half_power; ++half_power) {
        for (int q_half = 1;
             q_half <= maximum_q_half; ++q_half) {
            const std::vector<Integer> profile =
                stable_profile(q_half, 2 * half_power);
            const std::vector<Integer> previous =
                stable_profile(q_half, 2 * half_power - 2);
            const auto previous_at = [&](int index) -> Integer {
                if (
                    index < 0
                    || static_cast<std::size_t>(index)
                        >= previous.size()
                ) {
                    return 0;
                }
                return previous[static_cast<std::size_t>(index)];
            };
            const int transfer_label = 2 * q_half;
            ++profiles;
            for (int label_half = 1;
                 label_half < 3 * q_half; ++label_half) {
                const Integer curvature =
                    Integer{2}
                        * profile[static_cast<std::size_t>(label_half)]
                    - profile[
                        static_cast<std::size_t>(label_half - 1)
                    ]
                    - profile[
                        static_cast<std::size_t>(label_half + 1)
                    ];
                ++concavity_rows;
                if (curvature < 0) {
                    ++negative_concavity_rows;
                    if (first_concavity_half_power == 0) {
                        first_concavity_half_power = half_power;
                        first_concavity_q_half = q_half;
                        first_concavity_label_half = label_half;
                        first_concavity_value = curvature;
                    }
                }
                if (label_half <= 2 * q_half) {
                    ++low_concavity_rows;
                    if (curvature < 0) {
                        ++negative_low_concavity_rows;
                    }
                    const Integer transfer_tail =
                        Integer{2} * previous_at(label_half)
                        + previous_at(transfer_label - label_half)
                        - previous_at(
                            transfer_label + label_half + 1
                        );
                    ++transfer_identity_rows;
                    if (curvature != transfer_tail) {
                        ++failed_transfer_identity_rows;
                    }
                    if (transfer_tail < 0) {
                        ++negative_transfer_tail_rows;
                    }
                }
            }

            for (int distance_half = q_half + 1;
                 distance_half < 2 * q_half; ++distance_half) {
                const int low = distance_half - q_half - 1;
                const int high = q_half + distance_half;
                const Integer slope_gap =
                    profile[static_cast<std::size_t>(low + 1)]
                    - profile[static_cast<std::size_t>(low)]
                    - profile[static_cast<std::size_t>(high + 1)]
                    + profile[static_cast<std::size_t>(high)];
                ++paired_slope_rows;
                if (slope_gap < 0) {
                    ++negative_paired_slope_rows;
                }
                Integer block_payment = 0;
                for (int index = low + 1;
                     index <= transfer_label - low - 1; ++index) {
                    block_payment += Integer{3}
                        * previous_at(index);
                }
                for (int index = transfer_label - low;
                     index <= transfer_label + low + 1; ++index) {
                    block_payment += Integer{2}
                        * previous_at(index);
                }
                for (int index = transfer_label + low + 2;
                     index <= 2 * transfer_label + low + 2; ++index) {
                    block_payment -= previous_at(index);
                }
                ++block_identity_rows;
                if (slope_gap != block_payment) {
                    ++failed_block_identity_rows;
                }
                if (block_payment < 0) {
                    ++negative_block_payment_rows;
                }
            }

            ++dimension_rows;
            if (
                profile[static_cast<std::size_t>(q_half)]
                > Integer{2 * q_half + 1} * profile[0]
            ) {
                ++failed_dimension_rows;
            }

            for (int distance_half = 1;
                 distance_half <= 2 * q_half; ++distance_half) {
                Integer left;
                if (distance_half <= q_half) {
                    left =
                        profile[static_cast<std::size_t>(
                            q_half - distance_half
                        )]
                        + profile[static_cast<std::size_t>(
                            q_half + distance_half
                        )];
                } else {
                    left =
                        profile[static_cast<std::size_t>(
                            q_half + distance_half
                        )]
                        - profile[static_cast<std::size_t>(
                            distance_half - q_half - 1
                        )];
                }
                ++odlp_rows;
                if (
                    left
                    > Integer{2}
                        * profile[static_cast<std::size_t>(q_half)]
                ) {
                    ++failed_odlp_rows;
                }
            }
        }
    }

    std::cout
        << "SU2_TERMINAL_EVEN_PROFILE_CONCAVITY"
        << " maximum_half_power=" << maximum_half_power
        << " maximum_q_half=" << maximum_q_half
        << " profiles=" << profiles
        << " concavity_rows=" << concavity_rows
        << " negative_concavity_rows=" << negative_concavity_rows
        << " low_concavity_rows=" << low_concavity_rows
        << " negative_low_concavity_rows="
        << negative_low_concavity_rows
        << " paired_slope_rows=" << paired_slope_rows
        << " negative_paired_slope_rows="
        << negative_paired_slope_rows
        << " transfer_identity_rows=" << transfer_identity_rows
        << " failed_transfer_identity_rows="
        << failed_transfer_identity_rows
        << " negative_transfer_tail_rows="
        << negative_transfer_tail_rows
        << " block_identity_rows=" << block_identity_rows
        << " failed_block_identity_rows="
        << failed_block_identity_rows
        << " negative_block_payment_rows="
        << negative_block_payment_rows
        << " dimension_rows=" << dimension_rows
        << " failed_dimension_rows=" << failed_dimension_rows
        << " odlp_rows=" << odlp_rows
        << " failed_odlp_rows=" << failed_odlp_rows;
    if (first_concavity_half_power != 0) {
        std::cout
            << " first_concavity_failure={half_power="
            << first_concavity_half_power
            << " q_half=" << first_concavity_q_half
            << " label_half=" << first_concavity_label_half
            << " curvature=" << first_concavity_value << '}';
    }
    std::cout
        << " result="
        << (
            negative_low_concavity_rows == 0
                && negative_paired_slope_rows == 0
                && failed_transfer_identity_rows == 0
                && negative_transfer_tail_rows == 0
                && failed_block_identity_rows == 0
                && negative_block_payment_rows == 0
                && failed_dimension_rows == 0
                && failed_odlp_rows == 0
                ? "PASS_BOUNDED_DIAGNOSTIC"
                : "COUNTEREXAMPLE"
        )
        << '\n';
    return negative_low_concavity_rows == 0
            && negative_paired_slope_rows == 0
            && failed_transfer_identity_rows == 0
            && negative_transfer_tail_rows == 0
            && failed_block_identity_rows == 0
            && negative_block_payment_rows == 0
            && failed_dimension_rows == 0
            && failed_odlp_rows == 0
        ? EXIT_SUCCESS
        : EXIT_FAILURE;
}

std::size_t negative_coefficients(const Polynomial& polynomial) {
    std::size_t negative = 0;
    for (const auto& [exponent, coefficient] :
         polynomial.coefficients()) {
        static_cast<void>(exponent);
        if (coefficient < 0) {
            ++negative;
        }
    }
    return negative;
}

int analyze_cones(int maximum) {
    const Polynomial x = Polynomial::variable(0);
    const Polynomial y = Polynomial::variable(1);
    std::size_t cases = 0;
    std::size_t coefficients = 0;
    std::size_t negative = 0;
    for (int half_distance = 2;
         half_distance <= maximum; ++half_distance) {
        const auto inspect = [&](const Polynomial& polynomial) {
            ++cases;
            coefficients += polynomial.coefficients().size();
            negative += negative_coefficients(polynomial);
        };
        inspect(margin_odd(
            half_distance,
            x + y + Polynomial{Rational{half_distance}},
            x,
            true
        ));
        for (int offset = -half_distance + 1;
             offset <= half_distance - 1; ++offset) {
            const int base = std::max(1, std::abs(offset));
            inspect(margin_odd(
                half_distance,
                x + Polynomial{Rational{base}},
                x + Polynomial{Rational{base + offset}},
                false
            ));
        }
        inspect(margin_odd(
            half_distance,
            x + y + Polynomial{Rational{half_distance}},
            Rational{2} * x + y
                + Polynomial{Rational{2 * half_distance}},
            false
        ));
    }
    std::cout
        << "SU2_TERMINAL_ODD_LOW_CONES"
        << " maximum_half_distance=" << maximum
        << " maximum_distance=" << 2 * maximum + 1
        << " cases=" << cases
        << " coefficients=" << coefficients
        << " negative_coefficients=" << negative
        << " result="
        << (
            negative == 0
                ? "PASS_BOUNDED_SYMBOLIC_DIAGNOSTIC"
                : "COUNTEREXAMPLE_TO_COEFFICIENT_CONE"
        )
        << '\n';
    return negative == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--bootstrap") {
            return analyze_distance_five_bootstrap(
                parse_maximum(argv[2])
            );
        }
        if (argc == 4 && std::string(argv[1]) == "--profile") {
            return analyze_profile_concavity(
                parse_maximum(argv[2]),
                parse_maximum(argv[3])
            );
        }
        if (argc == 2) {
            return analyze_cones(parse_maximum(argv[1]));
        }
        if (argc != 1) {
            throw std::runtime_error(
                "usage: prove_su2_terminal_odd_low_fixed "
                "[MAXIMUM_HALF_DISTANCE] | "
                "--bootstrap MAXIMUM_Q_HALF | "
                "--profile MAXIMUM_HALF_POWER MAXIMUM_Q_HALF"
            );
        }
        const Polynomial x = Polynomial::variable(0);
        const Polynomial y = Polynomial::variable(1);
        bool passed = true;

        passed = report(
            "d5_second_cap",
            margin(x + y + Polynomial{Rational{2}}, x, true)
        ) && passed;
        passed = report(
            "d5_first_cap_lower_0",
            margin(x + Polynomial{Rational{1}}, x, false)
        ) && passed;
        passed = report(
            "d5_first_cap_lower_1",
            margin(
                x + Polynomial{Rational{1}},
                x + Polynomial{Rational{1}},
                false
            )
        ) && passed;
        passed = report(
            "d5_first_cap_tail",
            margin(
                x + y + Polynomial{Rational{1}},
                Rational{2} * x + y + Polynomial{Rational{2}},
                false
            )
        ) && passed;
        passed = report(
            "d7_final_cap",
            margin_d7(x + y + Polynomial{Rational{3}}, x, true)
        ) && passed;
        passed = report(
            "d7_boundary_minus_2",
            margin_d7(x + Polynomial{Rational{2}}, x, false)
        ) && passed;
        passed = report(
            "d7_boundary_minus_1",
            margin_d7(x + Polynomial{Rational{1}}, x, false)
        ) && passed;
        passed = report(
            "d7_boundary_0",
            margin_d7(
                x + Polynomial{Rational{1}},
                x + Polynomial{Rational{1}},
                false
            )
        ) && passed;
        passed = report(
            "d7_boundary_plus_1",
            margin_d7(
                x + Polynomial{Rational{1}},
                x + Polynomial{Rational{2}},
                false
            )
        ) && passed;
        passed = report(
            "d7_boundary_plus_2",
            margin_d7(
                x + Polynomial{Rational{2}},
                x + Polynomial{Rational{4}},
                false
            )
        ) && passed;
        passed = report(
            "d7_tail",
            margin_d7(
                x + y + Polynomial{Rational{3}},
                Rational{2} * x + y + Polynomial{Rational{6}},
                false
            )
        ) && passed;

        std::cout
            << "SU2_TERMINAL_ODD_LOW result="
            << (passed ? "PASS_EXACT_CERTIFICATE" : "UNRESOLVED")
            << '\n';
        return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_TERMINAL_ODD_LOW FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
