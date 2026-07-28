#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(value);
}

bool fuses(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= std::min(
            source + label,
            2 * level - source - label
        )
        && ((source + label + target) & 1) == 0;
}

std::vector<Integer> multiply(
    int level,
    int label,
    const std::vector<Integer>& state
) {
    std::vector<Integer> next(static_cast<std::size_t>(level + 1));
    for (int source = 0; source <= level; ++source) {
        const Integer& coefficient =
            state[static_cast<std::size_t>(source)];
        if (coefficient == 0) {
            continue;
        }
        for (int target = 0; target <= level; ++target) {
            if (fuses(level, label, source, target)) {
                next[static_cast<std::size_t>(target)] += coefficient;
            }
        }
    }
    return next;
}

std::vector<Integer> multiply_series(
    const std::vector<Integer>& left,
    const std::vector<Integer>& right
) {
    const std::size_t size = left.size();
    std::vector<Integer> product(size);
    for (std::size_t i = 0U; i < size; ++i) {
        for (std::size_t j = 0U; i + j < size; ++j) {
            product[i + j] += left[i] * right[j];
        }
    }
    return product;
}

Integer binomial(int n, int r) {
    if (r < 0 || r > n) {
        return 0;
    }
    r = std::min(r, n - r);
    Integer result = 1;
    for (int index = 1; index <= r; ++index) {
        result *= n - r + index;
        result /= index;
    }
    return result;
}

struct Witness {
    int level = 0;
    int label = 0;
    int prefix = 0;
    int truncation = 0;
    int head = -1;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 7) {
            throw std::runtime_error(
                "usage: analyze_su2_reduced_current_payment "
                "MAXIMUM_LEVEL MAXIMUM_PREFIX "
                "[DUMP_LEVEL DUMP_LABEL DUMP_PREFIX DUMP_TRUNCATION]"
            );
        }
        const int maximum_level =
            parse_positive(argv[1], "maximum level");
        const int maximum_prefix =
            parse_positive(argv[2], "maximum prefix");
        if (maximum_prefix < 4) {
            throw std::runtime_error(
                "maximum prefix must be at least four"
            );
        }
        const int dump_level =
            argc == 7 ? parse_positive(argv[3], "dump level") : -1;
        const int dump_label =
            argc == 7 ? parse_positive(argv[4], "dump label") : -1;
        const int dump_prefix =
            argc == 7 ? parse_positive(argv[5], "dump prefix") : -1;
        const int dump_truncation =
            argc == 7 ? parse_positive(argv[6], "dump truncation") : -1;

        std::uint64_t parameters = 0U;
        std::uint64_t rows = 0U;
        std::uint64_t identities = 0U;
        std::uint64_t negative_targets = 0U;
        std::uint64_t unresolved_rows = 0U;
        std::uint64_t reduced_lag_two_coefficients = 0U;
        std::uint64_t reduced_lag_two_negatives = 0U;
        std::uint64_t reduced_second_lag_two_negatives = 0U;
        std::uint64_t reduced_second_lag_two_nonminimal_negatives = 0U;
        std::uint64_t milestone_lag_two_coefficients = 0U;
        std::uint64_t milestone_lag_two_negatives = 0U;
        std::uint64_t doubled_abel_coefficients = 0U;
        std::uint64_t doubled_abel_negatives = 0U;
        std::uint64_t doubled_abel_nonminimal_negatives = 0U;
        std::uint64_t fibonacci_residual_coefficients = 0U;
        std::uint64_t fibonacci_residual_negatives = 0U;
        std::uint64_t fibonacci_current_coefficients = 0U;
        std::uint64_t fibonacci_current_negatives = 0U;
        std::uint64_t safe_fibonacci_residual_coefficients = 0U;
        std::uint64_t safe_fibonacci_residual_negatives = 0U;
        std::uint64_t safe_fibonacci_current_coefficients = 0U;
        std::uint64_t safe_fibonacci_current_negatives = 0U;
        std::uint64_t qge2_fibonacci_current_negatives = 0U;
        std::uint64_t core_suffix_coordinates = 0U;
        std::uint64_t core_suffix_negatives = 0U;
        int maximum_required_head = -1;
        Witness maximum_witness;
        int maximum_safe_required_head = -1;
        Witness maximum_safe_witness;
        int maximum_qge2_safe_required_head = -1;
        Witness maximum_qge2_safe_witness;

        for (int level = 6; level <= maximum_level; level += 2) {
            for (int label = 2; 2 * label < level; label += 2) {
                ++parameters;
                const int maximum_degree = 2 * maximum_prefix + 2;
                std::vector<Integer> returns(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                std::vector<Integer> endpoints(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                std::vector<Integer> state(
                    static_cast<std::size_t>(level + 1)
                );
                state[0] = 1;
                returns[0] = 1;
                for (int power = 1;
                     power <= maximum_degree;
                     ++power) {
                    state = multiply(level, label, state);
                    returns[static_cast<std::size_t>(power)] = state[0];
                    endpoints[static_cast<std::size_t>(power)] =
                        state[static_cast<std::size_t>(level)];
                }

                std::vector<Integer> first_passage(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                for (int degree = 0;
                     degree <= maximum_degree;
                     ++degree) {
                    Integer value =
                        endpoints[static_cast<std::size_t>(degree)];
                    for (int index = 1; index <= degree; ++index) {
                        value -= returns[static_cast<std::size_t>(index)]
                            * first_passage[
                                static_cast<std::size_t>(degree - index)
                            ];
                    }
                    first_passage[static_cast<std::size_t>(degree)] =
                        value;
                }

                const int half_level = level / 2;
                const int half_label = label / 2;
                const int exponent =
                    (half_level - 1) / half_label;
                const int safe_fibonacci_exponent = exponent;
                const int minimum_first_passage_degree = exponent + 1;
                std::vector<Integer> reduced(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                for (int degree = 0;
                     degree <= maximum_degree;
                     ++degree) {
                    for (int difference = 0;
                         difference <= exponent
                            && difference <= degree;
                         ++difference) {
                        const Integer term =
                            binomial(exponent, difference)
                            * first_passage[
                                static_cast<std::size_t>(
                                    degree - difference
                                )];
                        reduced[static_cast<std::size_t>(degree)] +=
                            (difference & 1) == 0 ? term : -term;
                    }
                    if (reduced[static_cast<std::size_t>(degree)] < 0) {
                        throw std::runtime_error(
                            "negative reduced first-passage coefficient"
                        );
                    }
                }
                if (reduced[static_cast<std::size_t>(
                        minimum_first_passage_degree
                    )] == 0) {
                    throw std::runtime_error(
                        "incorrect minimum first-passage degree"
                    );
                }
                std::vector<Integer> denominator(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                denominator[0] = 1;
                std::vector<Integer> factor(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                factor[0] = 1;
                factor[1] = -1;
                denominator = multiply_series(denominator, factor);
                std::fill(factor.begin(), factor.end(), Integer(0));
                factor[0] = 1;
                factor[2] = -1;
                denominator = multiply_series(denominator, factor);
                std::fill(factor.begin(), factor.end(), Integer(0));
                factor[0] = 1;
                factor[1] = -1;
                factor[2] = -1;
                for (int copy = 1; copy < exponent; ++copy) {
                    denominator = multiply_series(denominator, factor);
                }
                const std::vector<Integer> fibonacci_residual =
                    multiply_series(first_passage, denominator);
                for (int degree = 0;
                     degree <= maximum_degree;
                     ++degree) {
                    ++fibonacci_residual_coefficients;
                    if (fibonacci_residual[
                            static_cast<std::size_t>(degree)
                        ] < 0) {
                        ++fibonacci_residual_negatives;
                        if (fibonacci_residual_negatives == 1U) {
                            std::cout
                                << "FIRST_NEGATIVE_FIBONACCI_RESIDUAL"
                                << " level=" << level
                                << " label=" << label
                                << " exponent=" << exponent
                                << " degree=" << degree
                                << " value="
                                    << fibonacci_residual[
                                        static_cast<std::size_t>(degree)
                                    ]
                                << '\n';
                        }
                    }
                }
                std::vector<Integer> safe_denominator(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                safe_denominator[0] = 1;
                std::fill(factor.begin(), factor.end(), Integer(0));
                factor[0] = 1;
                factor[1] = -1;
                factor[2] = -1;
                for (int copy = 0;
                     copy < safe_fibonacci_exponent;
                     ++copy) {
                    safe_denominator =
                        multiply_series(safe_denominator, factor);
                }
                const std::vector<Integer> safe_fibonacci_residual =
                    multiply_series(first_passage, safe_denominator);
                for (int degree = 0;
                     degree <= maximum_degree;
                     ++degree) {
                    ++safe_fibonacci_residual_coefficients;
                    if (safe_fibonacci_residual[
                            static_cast<std::size_t>(degree)
                        ] < 0) {
                        ++safe_fibonacci_residual_negatives;
                        if (safe_fibonacci_residual_negatives == 1U) {
                            std::cout
                                << "FIRST_NEGATIVE_SAFE_FIBONACCI_RESIDUAL"
                                << " level=" << level
                                << " label=" << label
                                << " exponent="
                                    << safe_fibonacci_exponent
                                << " degree=" << degree
                                << " value="
                                    << safe_fibonacci_residual[
                                        static_cast<std::size_t>(degree)
                                    ]
                                << '\n';
                        }
                    }
                }
                for (int degree = 0;
                     degree <= maximum_degree;
                     ++degree) {
                    const Integer lag_two_residual =
                        reduced[static_cast<std::size_t>(degree)]
                        - (
                            degree >= 2
                                ? reduced[
                                    static_cast<std::size_t>(degree - 2)
                                ]
                                : Integer(0)
                        );
                    ++reduced_lag_two_coefficients;
                    if (lag_two_residual < 0) {
                        ++reduced_lag_two_negatives;
                        if (reduced_lag_two_negatives == 1U) {
                            std::cout
                                << "FIRST_NEGATIVE_REDUCED_LAG_TWO"
                                << " level=" << level
                                << " label=" << label
                                << " degree=" << degree
                                << " value=" << lag_two_residual
                                << '\n';
                        }
                    }
                    Integer second_lag_two_residual =
                        reduced[static_cast<std::size_t>(degree)];
                    if (degree >= 2) {
                        second_lag_two_residual -=
                            2 * reduced[
                                static_cast<std::size_t>(degree - 2)
                            ];
                    }
                    if (degree >= 4) {
                        second_lag_two_residual +=
                            reduced[
                                static_cast<std::size_t>(degree - 4)
                            ];
                    }
                    if (second_lag_two_residual < 0) {
                        ++reduced_second_lag_two_negatives;
                        if (!(level == 6 && label == 2)) {
                            ++reduced_second_lag_two_nonminimal_negatives;
                        }
                        if (reduced_second_lag_two_negatives == 1U) {
                            std::cout
                                << "FIRST_NEGATIVE_REDUCED_SECOND_LAG_TWO"
                                << " level=" << level
                                << " label=" << label
                                << " degree=" << degree
                                << " value=" << second_lag_two_residual
                                << '\n';
                        }
                    }
                    Integer milestone_lag_two_residual = 0;
                    for (int difference = 0;
                         difference <= exponent
                            && 2 * difference <= degree;
                         ++difference) {
                        const Integer term =
                            binomial(exponent, difference)
                            * reduced[
                                static_cast<std::size_t>(
                                    degree - 2 * difference
                                )];
                        milestone_lag_two_residual +=
                            (difference & 1) == 0 ? term : -term;
                    }
                    ++milestone_lag_two_coefficients;
                    if (milestone_lag_two_residual < 0) {
                        ++milestone_lag_two_negatives;
                        if (milestone_lag_two_negatives == 1U) {
                            std::cout
                                << "FIRST_NEGATIVE_MILESTONE_LAG_TWO"
                                << " level=" << level
                                << " label=" << label
                                << " exponent=" << exponent
                                << " degree=" << degree
                                << " value=" << milestone_lag_two_residual
                                << '\n';
                        }
                    }
                }

                for (int prefix = 4;
                     prefix <= maximum_prefix;
                     ++prefix) {
                    const int n = 2 * prefix + 2;
                    std::vector<Integer> current(
                        static_cast<std::size_t>(n + 1)
                    );
                    for (int truncation = 0;
                         truncation <= prefix - 1;
                         ++truncation) {
                        const int even = 2 * truncation;
                        const int odd = even + 1;
                        const Integer weight = binomial(n - 1, even);
                        current[static_cast<std::size_t>(even)] +=
                            weight
                            * returns[static_cast<std::size_t>(even)];
                        current[static_cast<std::size_t>(odd)] -=
                            weight
                            * returns[static_cast<std::size_t>(odd)];
                        // Prefixes through t=3 are already proved
                        // uniformly; this diagnostic targets R6.
                        if (truncation < 4) {
                            continue;
                        }
                        ++rows;

                        std::vector<Integer> h(
                            static_cast<std::size_t>(n + 1)
                        );
                        for (int degree = 0; degree <= n; ++degree) {
                            for (int index = 0; index <= degree; ++index) {
                                h[static_cast<std::size_t>(degree)] +=
                                    current[
                                        static_cast<std::size_t>(index)
                                    ] * returns[
                                        static_cast<std::size_t>(
                                            degree - index
                                        )];
                            }
                        }
                        std::vector<Integer> abel(
                            static_cast<std::size_t>(n + 1)
                        );
                        for (int degree = 0; degree <= n; ++degree) {
                            for (int index = 0; index <= degree; ++index) {
                                abel[static_cast<std::size_t>(degree)] +=
                                    h[static_cast<std::size_t>(index)]
                                    * binomial(
                                        degree - index + exponent - 1,
                                        exponent - 1
                                );
                            }
                        }
                        std::vector<Integer> doubled_abel(
                            static_cast<std::size_t>(n + 1)
                        );
                        for (int degree = 0; degree <= n; ++degree) {
                            for (int pair_degree = 0;
                                 2 * pair_degree <= degree;
                                 ++pair_degree) {
                                doubled_abel[
                                    static_cast<std::size_t>(degree)
                                ] += abel[
                                    static_cast<std::size_t>(
                                        degree - 2 * pair_degree
                                    )] * binomial(
                                        pair_degree + exponent - 1,
                                        exponent - 1
                                    );
                            }
                            ++doubled_abel_coefficients;
                            if (doubled_abel[
                                    static_cast<std::size_t>(degree)
                                ] < 0) {
                                ++doubled_abel_negatives;
                                if (!(level == 6 && label == 2)) {
                                    ++doubled_abel_nonminimal_negatives;
                                }
                                if (doubled_abel_negatives == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_DOUBLED_ABEL"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " exponent=" << exponent
                                        << " degree=" << degree
                                        << " value="
                                            << doubled_abel[
                                                static_cast<std::size_t>(
                                                    degree
                                                )]
                                        << '\n';
                                }
                            }
                        }
                        std::vector<Integer> safe_fibonacci_current = h;
                        for (int copy = 0;
                             copy < safe_fibonacci_exponent;
                             ++copy) {
                            std::vector<Integer> quotient(
                                static_cast<std::size_t>(n + 1)
                            );
                            for (int degree = 0;
                                 degree <= n;
                                 ++degree) {
                                quotient[
                                    static_cast<std::size_t>(degree)
                                ] = safe_fibonacci_current[
                                    static_cast<std::size_t>(degree)
                                ];
                                if (degree >= 1) {
                                    quotient[
                                        static_cast<std::size_t>(degree)
                                    ] += quotient[
                                        static_cast<std::size_t>(
                                            degree - 1
                                        )];
                                }
                                if (degree >= 2) {
                                    quotient[
                                        static_cast<std::size_t>(degree)
                                    ] += quotient[
                                        static_cast<std::size_t>(
                                            degree - 2
                                        )];
                                }
                            }
                            safe_fibonacci_current = std::move(quotient);
                        }
                        for (int degree = 0; degree <= n; ++degree) {
                            ++safe_fibonacci_current_coefficients;
                            if (safe_fibonacci_current[
                                    static_cast<std::size_t>(degree)
                                ] < 0) {
                                ++safe_fibonacci_current_negatives;
                                if (half_label >= 2) {
                                    ++qge2_fibonacci_current_negatives;
                                    if (qge2_fibonacci_current_negatives
                                        == 1U) {
                                        std::cout
                                            << "FIRST_NEGATIVE_QGE2_FIBONACCI_CURRENT"
                                            << " level=" << level
                                            << " label=" << label
                                            << " prefix=" << prefix
                                            << " truncation="
                                                << truncation
                                            << " exponent="
                                                << safe_fibonacci_exponent
                                            << " degree=" << degree
                                            << " value="
                                                << safe_fibonacci_current[
                                                    static_cast<
                                                        std::size_t
                                                    >(degree)
                                                ]
                                            << '\n';
                                    }
                                }
                                if (safe_fibonacci_current_negatives
                                    == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_SAFE_FIBONACCI_CURRENT"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " exponent="
                                            << safe_fibonacci_exponent
                                        << " degree=" << degree
                                        << " value="
                                            << safe_fibonacci_current[
                                                static_cast<std::size_t>(
                                                    degree
                                                )]
                                        << '\n';
                                }
                            }
                        }
                        Integer core_suffix = 0;
                        for (int residual_degree = n;
                             residual_degree
                                >= minimum_first_passage_degree;
                             --residual_degree) {
                            core_suffix += safe_fibonacci_residual[
                                static_cast<std::size_t>(residual_degree)
                            ] * safe_fibonacci_current[
                                static_cast<std::size_t>(
                                    n - residual_degree
                                )];
                            ++core_suffix_coordinates;
                            if (core_suffix < 0) {
                                ++core_suffix_negatives;
                                if (core_suffix_negatives == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_CORE_SUFFIX"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " residual_degree="
                                            << residual_degree
                                        << " value=" << core_suffix
                                        << '\n';
                                }
                            }
                        }
                        std::vector<Integer> fibonacci_current = h;
                        for (int degree = 1; degree <= n; ++degree) {
                            fibonacci_current[
                                static_cast<std::size_t>(degree)
                            ] += fibonacci_current[
                                static_cast<std::size_t>(degree - 1)
                            ];
                        }
                        for (int degree = 2; degree <= n; ++degree) {
                            fibonacci_current[
                                static_cast<std::size_t>(degree)
                            ] += fibonacci_current[
                                static_cast<std::size_t>(degree - 2)
                            ];
                        }
                        for (int copy = 1; copy < exponent; ++copy) {
                            std::vector<Integer> quotient(
                                static_cast<std::size_t>(n + 1)
                            );
                            for (int degree = 0;
                                 degree <= n;
                                 ++degree) {
                                quotient[
                                    static_cast<std::size_t>(degree)
                                ] = fibonacci_current[
                                    static_cast<std::size_t>(degree)
                                ];
                                if (degree >= 1) {
                                    quotient[
                                        static_cast<std::size_t>(degree)
                                    ] += quotient[
                                        static_cast<std::size_t>(
                                            degree - 1
                                        )];
                                }
                                if (degree >= 2) {
                                    quotient[
                                        static_cast<std::size_t>(degree)
                                    ] += quotient[
                                        static_cast<std::size_t>(
                                            degree - 2
                                        )];
                                }
                            }
                            fibonacci_current = std::move(quotient);
                        }
                        for (int degree = 0; degree <= n; ++degree) {
                            ++fibonacci_current_coefficients;
                            if (fibonacci_current[
                                    static_cast<std::size_t>(degree)
                                ] < 0) {
                                ++fibonacci_current_negatives;
                                if (fibonacci_current_negatives == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_FIBONACCI_CURRENT"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " exponent=" << exponent
                                        << " degree=" << degree
                                        << " value="
                                            << fibonacci_current[
                                                static_cast<std::size_t>(
                                                    degree
                                                )]
                                        << '\n';
                                }
                            }
                        }

                        Integer factored_target = 0;
                        for (int degree = 0; degree <= n; ++degree) {
                            factored_target +=
                                reduced[static_cast<std::size_t>(degree)]
                                * abel[
                                    static_cast<std::size_t>(n - degree)
                                ];
                        }
                        Integer direct_target = 0;
                        for (int degree = 0; degree <= n; ++degree) {
                            direct_target +=
                                current[static_cast<std::size_t>(degree)]
                                * endpoints[
                                    static_cast<std::size_t>(n - degree)
                                ];
                        }
                        if (factored_target != direct_target) {
                            throw std::runtime_error(
                                "factored-current identity failed"
                            );
                        }
                        ++identities;
                        if (direct_target < 0) {
                            ++negative_targets;
                        }

                        const int available =
                            n - minimum_first_passage_degree;
                        int required_head = available < 0 ? 0 : -1;
                        Integer head_sum = 0;
                        for (int head = 0; head <= available; ++head) {
                            const int reduced_degree =
                                minimum_first_passage_degree + head;
                            head_sum += reduced[
                                static_cast<std::size_t>(reduced_degree)
                            ] * abel[
                                static_cast<std::size_t>(
                                    n - reduced_degree
                                )];
                            bool tail_safe = true;
                            const int maximum_tail_abel_degree =
                                available - head - 1;
                            for (int degree = 0;
                                 degree <= maximum_tail_abel_degree;
                                 ++degree) {
                                if (abel[
                                        static_cast<std::size_t>(degree)
                                    ] < 0) {
                                    tail_safe = false;
                                    break;
                                }
                            }
                            if (head_sum >= 0 && tail_safe) {
                                required_head = head;
                                break;
                            }
                        }
                        if (required_head < 0) {
                            ++unresolved_rows;
                        } else if (required_head > maximum_required_head) {
                            maximum_required_head = required_head;
                            maximum_witness = {
                                level,
                                label,
                                prefix,
                                truncation,
                                required_head
                            };
                        }
                        int safe_required_head =
                            available < 0 ? 0 : -1;
                        Integer safe_head_sum = 0;
                        for (int head = 0;
                             head <= available;
                             ++head) {
                            const int residual_degree =
                                minimum_first_passage_degree + head;
                            safe_head_sum += safe_fibonacci_residual[
                                static_cast<std::size_t>(residual_degree)
                            ] * safe_fibonacci_current[
                                static_cast<std::size_t>(
                                    n - residual_degree
                                )];
                            bool safe_tail = true;
                            const int maximum_tail_degree =
                                available - head - 1;
                            for (int degree = 0;
                                 degree <= maximum_tail_degree;
                                 ++degree) {
                                if (safe_fibonacci_current[
                                        static_cast<std::size_t>(degree)
                                    ] < 0) {
                                    safe_tail = false;
                                    break;
                                }
                            }
                            if (safe_head_sum >= 0 && safe_tail) {
                                safe_required_head = head;
                                break;
                            }
                        }
                        if (safe_required_head
                            > maximum_safe_required_head) {
                            maximum_safe_required_head =
                                safe_required_head;
                            maximum_safe_witness = {
                                level,
                                label,
                                prefix,
                                truncation,
                                safe_required_head
                            };
                        }
                        if (half_label >= 2
                            && safe_required_head
                                > maximum_qge2_safe_required_head) {
                            maximum_qge2_safe_required_head =
                                safe_required_head;
                            maximum_qge2_safe_witness = {
                                level,
                                label,
                                prefix,
                                truncation,
                                safe_required_head
                            };
                        }
                        if (level == dump_level
                            && label == dump_label
                            && prefix == dump_prefix
                            && truncation == dump_truncation) {
                            std::cout
                                << "REDUCED_PAYMENT_DUMP"
                                << " level=" << level
                                << " label=" << label
                                << " prefix=" << prefix
                                << " truncation=" << truncation
                                << " exponent=" << exponent
                                << " minimum_degree="
                                    << minimum_first_passage_degree
                                << " required_head=" << required_head
                                << " target=" << direct_target << '\n';
                            for (int head = 0;
                                 head <= available;
                                 ++head) {
                                const int reduced_degree =
                                    minimum_first_passage_degree + head;
                                const int abel_degree =
                                    n - reduced_degree;
                                std::cout
                                    << "head=" << head
                                    << " reduced="
                                        << reduced[
                                            static_cast<std::size_t>(
                                                reduced_degree
                                            )]
                                    << " safe_residual="
                                        << safe_fibonacci_residual[
                                            static_cast<std::size_t>(
                                                reduced_degree
                                            )]
                                    << " abel_degree=" << abel_degree
                                    << " abel="
                                        << abel[
                                            static_cast<std::size_t>(
                                                abel_degree
                                            )]
                                    << " safe_current="
                                        << safe_fibonacci_current[
                                            static_cast<std::size_t>(
                                                abel_degree
                                            )]
                                    << '\n';
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "SU2_REDUCED_CURRENT_PAYMENT"
            << " maximum_level=" << maximum_level
            << " maximum_prefix=" << maximum_prefix
            << " parameters=" << parameters
            << " rows=" << rows
            << " identities=" << identities
            << " negative_targets=" << negative_targets
            << " unresolved_rows=" << unresolved_rows
            << " reduced_lag_two_coefficients="
                << reduced_lag_two_coefficients
            << " reduced_lag_two_negatives="
                << reduced_lag_two_negatives
            << " reduced_second_lag_two_negatives="
                << reduced_second_lag_two_negatives
            << " reduced_second_lag_two_nonminimal_negatives="
                << reduced_second_lag_two_nonminimal_negatives
            << " milestone_lag_two_coefficients="
                << milestone_lag_two_coefficients
            << " milestone_lag_two_negatives="
                << milestone_lag_two_negatives
            << " doubled_abel_coefficients="
                << doubled_abel_coefficients
            << " doubled_abel_negatives="
                << doubled_abel_negatives
            << " doubled_abel_nonminimal_negatives="
                << doubled_abel_nonminimal_negatives
            << " fibonacci_residual_coefficients="
                << fibonacci_residual_coefficients
            << " fibonacci_residual_negatives="
                << fibonacci_residual_negatives
            << " fibonacci_current_coefficients="
                << fibonacci_current_coefficients
            << " fibonacci_current_negatives="
                << fibonacci_current_negatives
            << " safe_fibonacci_residual_coefficients="
                << safe_fibonacci_residual_coefficients
            << " safe_fibonacci_residual_negatives="
                << safe_fibonacci_residual_negatives
            << " safe_fibonacci_current_coefficients="
                << safe_fibonacci_current_coefficients
            << " safe_fibonacci_current_negatives="
                << safe_fibonacci_current_negatives
            << " qge2_fibonacci_current_negatives="
                << qge2_fibonacci_current_negatives
            << " core_suffix_coordinates="
                << core_suffix_coordinates
            << " core_suffix_negatives="
                << core_suffix_negatives
            << " maximum_required_head=" << maximum_required_head
            << " witness=("
            << maximum_witness.level << ','
            << maximum_witness.label << ','
            << maximum_witness.prefix << ','
            << maximum_witness.truncation << ','
            << maximum_witness.head << ')'
            << " maximum_safe_required_head="
                << maximum_safe_required_head
            << " safe_witness=("
            << maximum_safe_witness.level << ','
            << maximum_safe_witness.label << ','
            << maximum_safe_witness.prefix << ','
            << maximum_safe_witness.truncation << ','
            << maximum_safe_witness.head << ')'
            << " maximum_qge2_safe_required_head="
                << maximum_qge2_safe_required_head
            << " qge2_safe_witness=("
            << maximum_qge2_safe_witness.level << ','
            << maximum_qge2_safe_witness.label << ','
            << maximum_qge2_safe_witness.prefix << ','
            << maximum_qge2_safe_witness.truncation << ','
            << maximum_qge2_safe_witness.head << ')'
            << " result="
            << (
                negative_targets == 0U
                    && unresolved_rows == 0U
                    && reduced_lag_two_negatives == 0U
                    && milestone_lag_two_negatives == 0U
                    && fibonacci_residual_negatives == 0U
                    && safe_fibonacci_residual_negatives == 0U
                    ? "PASS_REDUCED_PAYMENT_DISCOVERY"
                    : "FAIL_REDUCED_PAYMENT_CANDIDATE"
            )
            << '\n';
        return negative_targets == 0U
                && unresolved_rows == 0U
                && reduced_lag_two_negatives == 0U
                && milestone_lag_two_negatives == 0U
                && fibonacci_residual_negatives == 0U
                && safe_fibonacci_residual_negatives == 0U
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_REDUCED_CURRENT_PAYMENT FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
