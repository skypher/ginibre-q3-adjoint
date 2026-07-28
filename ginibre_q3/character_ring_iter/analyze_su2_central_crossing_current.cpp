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
using Matrix = std::vector<std::vector<int>>;

int parse_positive(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0
        || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string(name) + " must be positive");
    }
    return static_cast<int>(value);
}

bool fuses_half(
    int half_level,
    int half_label,
    int source,
    int target
) {
    return std::abs(source - half_label) <= target
        && target <= source + half_label
        && source + target + half_label <= 2 * half_level;
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

std::vector<Integer> multiply_row(
    const std::vector<Integer>& state,
    const Matrix& matrix
) {
    const int size = static_cast<int>(matrix.size());
    std::vector<Integer> next(static_cast<std::size_t>(size));
    for (int source = 0; source < size; ++source) {
        if (state[static_cast<std::size_t>(source)] == 0) {
            continue;
        }
        for (int target = 0; target < size; ++target) {
            const int entry =
                matrix[static_cast<std::size_t>(source)]
                      [static_cast<std::size_t>(target)];
            if (entry != 0) {
                next[static_cast<std::size_t>(target)] +=
                    state[static_cast<std::size_t>(source)] * entry;
            }
        }
    }
    return next;
}

std::vector<Integer> multiply_column(
    const Matrix& matrix,
    const std::vector<Integer>& state
) {
    const int size = static_cast<int>(matrix.size());
    std::vector<Integer> next(static_cast<std::size_t>(size));
    for (int target = 0; target < size; ++target) {
        for (int source = 0; source < size; ++source) {
            const int entry =
                matrix[static_cast<std::size_t>(target)]
                      [static_cast<std::size_t>(source)];
            if (entry != 0
                && state[static_cast<std::size_t>(source)] != 0) {
                next[static_cast<std::size_t>(target)] +=
                    entry * state[static_cast<std::size_t>(source)];
            }
        }
    }
    return next;
}

struct Witness {
    bool initialized = false;
    Integer value = 0;
    int level = 0;
    int label = 0;
    int prefix = 0;
    int truncation = 0;
    int crossing_source = 0;
    int crossing_target = 0;

    void observe(
        const Integer& candidate,
        int candidate_level,
        int candidate_label,
        int candidate_prefix,
        int candidate_truncation,
        int candidate_source,
        int candidate_target
    ) {
        if (!initialized || candidate < value) {
            initialized = true;
            value = candidate;
            level = candidate_level;
            label = candidate_label;
            prefix = candidate_prefix;
            truncation = candidate_truncation;
            crossing_source = candidate_source;
            crossing_target = candidate_target;
        }
    }
};

}  // namespace

int main(int argc, char** argv) {
    try {
        const bool case_mode =
            argc == 6 && std::string(argv[1]) == "--case";
        const bool suffix_only =
            argc == 4 && std::string(argv[3]) == "--suffix-only";
        if (argc != 3 && !case_mode && !suffix_only) {
            throw std::runtime_error(
                "usage: analyze_su2_central_crossing_current "
                "MAXIMUM_LEVEL MAXIMUM_PREFIX [--suffix-only] | "
                "--case LEVEL LABEL PREFIX TRUNCATION"
            );
        }
        const int case_level =
            case_mode ? parse_positive(argv[2], "case level") : 0;
        const int case_label =
            case_mode ? parse_positive(argv[3], "case label") : 0;
        const int case_prefix =
            case_mode ? parse_positive(argv[4], "case prefix") : 0;
        const int case_truncation =
            case_mode ? parse_positive(argv[5], "case truncation") : 0;
        const int maximum_level =
            case_mode ? case_level
                      : parse_positive(argv[1], "maximum level");
        const int maximum_prefix =
            case_mode ? case_prefix
                      : parse_positive(argv[2], "maximum prefix");
        if (maximum_prefix < 4) {
            throw std::runtime_error(
                "maximum prefix must be at least four"
            );
        }

        std::uint64_t parameters = 0U;
        std::uint64_t quotient_entries = 0U;
        std::uint64_t crossing_entries = 0U;
        std::uint64_t component_rows = 0U;
        std::uint64_t component_negatives = 0U;
        std::uint64_t target_rows = 0U;
        std::uint64_t target_negatives = 0U;
        std::uint64_t target_sign_recrossings = 0U;
        std::uint64_t target_sign_reversals = 0U;
        int maximum_target_sign_changes = 0;
        bool printed_first_target_recrossing = false;
        std::uint64_t target_suffix_rows = 0U;
        std::uint64_t target_suffix_negatives = 0U;
        std::uint64_t adjacent_target_rows = 0U;
        std::uint64_t adjacent_target_negatives = 0U;
        std::uint64_t inner_reserve_rows = 0U;
        std::uint64_t inner_reserve_negatives = 0U;
        std::uint64_t endpoint_identity_rows = 0U;
        Witness minimum;
        Witness target_minimum;
        Witness target_suffix_minimum;
        Witness adjacent_target_minimum;
        Witness inner_reserve_minimum;
        Integer case_total = 0;

        for (int level = 6; level <= maximum_level; level += 2) {
            if (case_mode && level != case_level) {
                continue;
            }
            const int half_level = level / 2;
            for (int label = 2; 2 * label < level; label += 2) {
                if (case_mode && label != case_label) {
                    continue;
                }
                const int half_label = label / 2;
                ++parameters;
                const int paired = (half_level + 1) / 2;
                const bool has_center = (half_level & 1) == 0;
                const int plus_size = paired + (has_center ? 1 : 0);
                const int center = has_center ? paired : -1;

                Matrix even_quotient(
                    static_cast<std::size_t>(plus_size),
                    std::vector<int>(
                        static_cast<std::size_t>(plus_size),
                        0
                    )
                );
                Matrix odd_quotient(
                    static_cast<std::size_t>(plus_size),
                    std::vector<int>(
                        static_cast<std::size_t>(plus_size),
                        0
                    )
                );
                for (int source = 0; source < paired; ++source) {
                    for (int target = 0; target < paired; ++target) {
                        const int same = fuses_half(
                            half_level,
                            half_label,
                            source,
                            target
                        ) ? 1 : 0;
                        const int crossed = fuses_half(
                            half_level,
                            half_label,
                            source,
                            half_level - target
                        ) ? 1 : 0;
                        if (crossed > same) {
                            throw std::runtime_error(
                                "cross-midline edge lacks same-side edge"
                            );
                        }
                        even_quotient[
                            static_cast<std::size_t>(source)
                        ][static_cast<std::size_t>(target)] =
                            same + crossed;
                        odd_quotient[
                            static_cast<std::size_t>(source)
                        ][static_cast<std::size_t>(target)] =
                            same - crossed;
                        ++quotient_entries;
                    }
                    if (has_center) {
                        const int joins_center = fuses_half(
                            half_level,
                            half_label,
                            source,
                            half_level / 2
                        ) ? 1 : 0;
                        even_quotient[
                            static_cast<std::size_t>(source)
                        ][static_cast<std::size_t>(center)] =
                            joins_center;
                        even_quotient[
                            static_cast<std::size_t>(center)
                        ][static_cast<std::size_t>(source)] =
                            2 * joins_center;
                    }
                }
                if (has_center) {
                    even_quotient[
                        static_cast<std::size_t>(center)
                    ][static_cast<std::size_t>(center)] =
                        fuses_half(
                            half_level,
                            half_label,
                            half_level / 2,
                            half_level / 2
                        ) ? 1 : 0;
                }

                Matrix crossing = even_quotient;
                for (int source = 0;
                     source < plus_size;
                     ++source) {
                    for (int target = 0;
                         target < plus_size;
                         ++target) {
                        crossing[
                            static_cast<std::size_t>(source)
                        ][static_cast<std::size_t>(target)] -=
                            odd_quotient[
                                static_cast<std::size_t>(source)
                            ][static_cast<std::size_t>(target)];
                        if (crossing[
                                static_cast<std::size_t>(source)
                            ][static_cast<std::size_t>(target)] < 0) {
                            throw std::runtime_error(
                                "quotient difference is negative"
                            );
                        }
                    }
                }

                const int maximum_degree = 2 * maximum_prefix + 2;
                std::vector<std::vector<Integer>> even_rows(
                    static_cast<std::size_t>(maximum_degree + 1),
                    std::vector<Integer>(
                        static_cast<std::size_t>(plus_size)
                    )
                );
                std::vector<std::vector<Integer>> odd_columns(
                    static_cast<std::size_t>(maximum_degree + 1),
                    std::vector<Integer>(
                        static_cast<std::size_t>(plus_size)
                    )
                );
                even_rows[0][0] = 1;
                odd_columns[0][0] = 1;
                for (int degree = 1;
                     degree <= maximum_degree;
                     ++degree) {
                    even_rows[static_cast<std::size_t>(degree)] =
                        multiply_row(
                            even_rows[
                                static_cast<std::size_t>(degree - 1)
                            ],
                            even_quotient
                        );
                    odd_columns[static_cast<std::size_t>(degree)] =
                        multiply_column(
                            odd_quotient,
                            odd_columns[
                                static_cast<std::size_t>(degree - 1)
                            ]
                        );
                }

                std::vector<Integer> returns(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                std::vector<Integer> endpoints(
                    static_cast<std::size_t>(maximum_degree + 1)
                );
                for (int degree = 0;
                     degree <= maximum_degree;
                     ++degree) {
                    const Integer sum =
                        even_rows[
                            static_cast<std::size_t>(degree)
                        ][0]
                        + odd_columns[
                            static_cast<std::size_t>(degree)
                        ][0];
                    if ((sum & 1) != 0) {
                        throw std::runtime_error(
                            "nonintegral folded return coefficient"
                        );
                    }
                    returns[static_cast<std::size_t>(degree)] = sum / 2;
                    const Integer difference =
                        even_rows[
                            static_cast<std::size_t>(degree)
                        ][0]
                        - odd_columns[
                            static_cast<std::size_t>(degree)
                        ][0];
                    if ((difference & 1) != 0) {
                        throw std::runtime_error(
                            "nonintegral folded endpoint coefficient"
                        );
                    }
                    endpoints[static_cast<std::size_t>(degree)] =
                        difference / 2;
                }

                for (int crossing_source = 0;
                     !suffix_only && crossing_source < plus_size;
                     ++crossing_source) {
                    for (int crossing_target = 0;
                         crossing_target < paired;
                         ++crossing_target) {
                        const int crossing_weight = crossing[
                            static_cast<std::size_t>(crossing_source)
                        ][static_cast<std::size_t>(crossing_target)];
                        if (crossing_weight == 0) {
                            continue;
                        }
                        ++crossing_entries;
                        std::vector<Integer> crossing_kernel(
                            static_cast<std::size_t>(maximum_degree + 1)
                        );
                        for (int degree = 0;
                             degree <= maximum_degree;
                             ++degree) {
                            for (int left = 0;
                                 left <= degree;
                                 ++left) {
                                crossing_kernel[
                                    static_cast<std::size_t>(degree)
                                ] += even_rows[
                                    static_cast<std::size_t>(left)
                                ][static_cast<std::size_t>(
                                    crossing_source
                                )] * odd_columns[
                                    static_cast<std::size_t>(
                                        degree - left
                                    )
                                ][static_cast<std::size_t>(
                                    crossing_target
                                )] * crossing_weight;
                            }
                        }

                        for (int prefix = 4;
                             prefix <= maximum_prefix;
                             ++prefix) {
                            const int n = 2 * prefix + 2;
                            std::vector<Integer> polynomial(
                                static_cast<std::size_t>(n + 1)
                            );
                            for (int truncation = 0;
                                 truncation <= prefix - 1;
                                 ++truncation) {
                                const int even = 2 * truncation;
                                const int odd = even + 1;
                                const Integer weight =
                                    binomial(n - 1, even);
                                polynomial[
                                    static_cast<std::size_t>(even)
                                ] += weight * returns[
                                    static_cast<std::size_t>(even)
                                ];
                                polynomial[
                                    static_cast<std::size_t>(odd)
                                ] -= weight * returns[
                                    static_cast<std::size_t>(odd)
                                ];
                                if (truncation < 2) {
                                    continue;
                                }
                                Integer value = 0;
                                for (int degree = 0;
                                     degree <= n - 1;
                                     ++degree) {
                                    value += polynomial[
                                        static_cast<std::size_t>(degree)
                                    ] * crossing_kernel[
                                        static_cast<std::size_t>(
                                            n - 1 - degree
                                        )
                                    ];
                                }
                                ++component_rows;
                                minimum.observe(
                                    value,
                                    level,
                                    label,
                                    prefix,
                                    truncation,
                                    crossing_source,
                                    crossing_target
                                );
                                if (
                                    case_mode
                                    && prefix == case_prefix
                                    && truncation == case_truncation
                                ) {
                                    std::cout
                                        << "CROSSING_COMPONENT"
                                        << " source=" << crossing_source
                                        << " target=" << crossing_target
                                        << " weight=" << crossing_weight
                                        << " value=" << value << '\n';
                                    case_total += value;
                                }
                                if (value < 0) {
                                    ++component_negatives;
                                    if (case_mode) {
                                        continue;
                                    }
                                    if (component_negatives == 1U) {
                                        std::cout
                                            << "FIRST_NEGATIVE_COMPONENT"
                                            << " level=" << level
                                            << " label=" << label
                                            << " prefix=" << prefix
                                            << " truncation=" << truncation
                                            << " crossing=("
                                            << crossing_source << ','
                                            << crossing_target << ')'
                                            << " weight=" << crossing_weight
                                            << " value=" << value << '\n';
                                    }
                                }
                            }
                        }
                    }
                }

                std::vector<std::vector<std::vector<Integer>>>
                    target_values(
                        static_cast<std::size_t>(paired),
                        std::vector<std::vector<Integer>>(
                            static_cast<std::size_t>(
                                maximum_prefix + 1
                            ),
                            std::vector<Integer>(
                                static_cast<std::size_t>(
                                    maximum_prefix
                                )
                            )
                        )
                    );
                for (int crossing_target = 0;
                     crossing_target < paired;
                     ++crossing_target) {
                    std::vector<Integer> target_kernel(
                        static_cast<std::size_t>(maximum_degree + 1)
                    );
                    for (int degree = 0;
                         degree <= maximum_degree;
                         ++degree) {
                        for (int left = 0; left <= degree; ++left) {
                            Integer left_weight = 0;
                            for (int crossing_source = 0;
                                 crossing_source < plus_size;
                                 ++crossing_source) {
                                left_weight += crossing[
                                    static_cast<std::size_t>(
                                        crossing_source
                                    )
                                ][static_cast<std::size_t>(
                                    crossing_target
                                )] * even_rows[
                                    static_cast<std::size_t>(left)
                                ][static_cast<std::size_t>(
                                    crossing_source
                                )];
                            }
                            target_kernel[
                                static_cast<std::size_t>(degree)
                            ] += left_weight * odd_columns[
                                static_cast<std::size_t>(degree - left)
                            ][static_cast<std::size_t>(crossing_target)];
                        }
                    }
                    for (int prefix = 4;
                         prefix <= maximum_prefix;
                         ++prefix) {
                        const int n = 2 * prefix + 2;
                        std::vector<Integer> polynomial(
                            static_cast<std::size_t>(n + 1)
                        );
                        for (int truncation = 0;
                             truncation <= prefix - 1;
                             ++truncation) {
                            const int even = 2 * truncation;
                            const int odd = even + 1;
                            const Integer weight =
                                binomial(n - 1, even);
                            polynomial[
                                static_cast<std::size_t>(even)
                            ] += weight * returns[
                                static_cast<std::size_t>(even)
                            ];
                            polynomial[
                                static_cast<std::size_t>(odd)
                            ] -= weight * returns[
                                static_cast<std::size_t>(odd)
                            ];
                            if (truncation < 2) {
                                continue;
                            }
                            Integer value = 0;
                            for (int degree = 0;
                                 degree <= n - 1;
                                 ++degree) {
                                value += polynomial[
                                    static_cast<std::size_t>(degree)
                                ] * target_kernel[
                                    static_cast<std::size_t>(
                                        n - 1 - degree
                                    )
                                ];
                            }
                            ++target_rows;
                            target_values[
                                static_cast<std::size_t>(crossing_target)
                            ][static_cast<std::size_t>(prefix)]
                             [static_cast<std::size_t>(truncation)] =
                                value;
                            target_minimum.observe(
                                value,
                                level,
                                label,
                                prefix,
                                truncation,
                                crossing_target,
                                crossing_target
                            );
                            if (value < 0) {
                                ++target_negatives;
                                if (target_negatives == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_TARGET"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " target=" << crossing_target
                                        << " value=" << value << '\n';
                                }
                            }
                        }
                    }
                }
                for (int prefix = 4;
                     prefix <= maximum_prefix;
                     ++prefix) {
                    for (int truncation = 2;
                         truncation <= prefix - 1;
                         ++truncation) {
                        int last_target_sign = 0;
                        int target_sign_changes = 0;
                        bool seen_positive_target = false;
                        for (int target = 0;
                             target < paired;
                             ++target) {
                            const Integer& value = target_values[
                                static_cast<std::size_t>(target)
                            ][static_cast<std::size_t>(prefix)]
                             [static_cast<std::size_t>(truncation)];
                            const int sign =
                                value < 0 ? -1 : (value > 0 ? 1 : 0);
                            if (sign == 0) {
                                continue;
                            }
                            if (
                                last_target_sign != 0
                                && sign != last_target_sign
                            ) {
                                ++target_sign_changes;
                            }
                            if (seen_positive_target && sign < 0) {
                                ++target_sign_reversals;
                            }
                            if (sign > 0) {
                                seen_positive_target = true;
                            }
                            last_target_sign = sign;
                        }
                        maximum_target_sign_changes = std::max(
                            maximum_target_sign_changes,
                            target_sign_changes
                        );
                        if (target_sign_changes > 1) {
                            ++target_sign_recrossings;
                            if (!printed_first_target_recrossing) {
                                printed_first_target_recrossing = true;
                                std::cout
                                    << "FIRST_TARGET_SIGN_RECROSSING"
                                    << " level=" << level
                                    << " label=" << label
                                    << " prefix=" << prefix
                                    << " truncation=" << truncation
                                    << " values=";
                                for (int target = 0;
                                     target < paired;
                                     ++target) {
                                    if (target != 0) {
                                        std::cout << ',';
                                    }
                                    std::cout << target_values[
                                        static_cast<std::size_t>(target)
                                    ][static_cast<std::size_t>(prefix)]
                                     [static_cast<std::size_t>(
                                         truncation
                                     )];
                                }
                                std::cout << '\n';
                            }
                        }
                        Integer inner_reserve = target_values[
                            static_cast<std::size_t>(paired - 1)
                        ][static_cast<std::size_t>(prefix)]
                         [static_cast<std::size_t>(truncation)];
                        for (int target = 0;
                             target < paired - 1;
                             ++target) {
                            const Integer& increment = target_values[
                                static_cast<std::size_t>(target)
                            ][static_cast<std::size_t>(prefix)]
                             [static_cast<std::size_t>(truncation)];
                            if (increment < 0) {
                                inner_reserve += increment;
                            }
                        }
                        ++inner_reserve_rows;
                        inner_reserve_minimum.observe(
                            inner_reserve,
                            level,
                            label,
                            prefix,
                            truncation,
                            paired - 1,
                            paired - 1
                        );
                        if (inner_reserve < 0) {
                            ++inner_reserve_negatives;
                            if (inner_reserve_negatives == 1U) {
                                std::cout
                                    << "FIRST_NEGATIVE_INNER_RESERVE"
                                    << " level=" << level
                                    << " label=" << label
                                    << " prefix=" << prefix
                                    << " truncation=" << truncation
                                    << " value=" << inner_reserve << '\n';
                            }
                        }
                        for (int target = paired - 2;
                             target >= 0;
                             --target) {
                            const Integer adjacent =
                                target_values[
                                    static_cast<std::size_t>(target)
                                ][static_cast<std::size_t>(prefix)]
                                 [static_cast<std::size_t>(truncation)]
                                + target_values[
                                    static_cast<std::size_t>(target + 1)
                                ][static_cast<std::size_t>(prefix)]
                                 [static_cast<std::size_t>(truncation)];
                            ++adjacent_target_rows;
                            adjacent_target_minimum.observe(
                                adjacent,
                                level,
                                label,
                                prefix,
                                truncation,
                                target,
                                target + 1
                            );
                            if (adjacent < 0) {
                                ++adjacent_target_negatives;
                                if (adjacent_target_negatives == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_ADJACENT_TARGET"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " target=" << target
                                        << " value=" << adjacent << '\n';
                                }
                            }
                        }
                        Integer suffix = 0;
                        for (int target = paired - 1;
                             target >= 0;
                             --target) {
                            const Integer& target_value = target_values[
                                static_cast<std::size_t>(target)
                            ][static_cast<std::size_t>(prefix)]
                             [static_cast<std::size_t>(truncation)];
                            suffix += target_value;
                            if (
                                case_mode
                                && prefix == case_prefix
                                && truncation == case_truncation
                            ) {
                                std::cout
                                    << "TARGET_SHELL"
                                    << " target=" << target
                                    << " increment=" << target_value
                                    << " suffix=" << suffix << '\n';
                            }
                            ++target_suffix_rows;
                            target_suffix_minimum.observe(
                                suffix,
                                level,
                                label,
                                prefix,
                                truncation,
                                target,
                                target
                            );
                            if (suffix < 0) {
                                ++target_suffix_negatives;
                                if (target_suffix_negatives == 1U) {
                                    std::cout
                                        << "FIRST_NEGATIVE_TARGET_SUFFIX"
                                        << " level=" << level
                                        << " label=" << label
                                        << " prefix=" << prefix
                                        << " truncation=" << truncation
                                        << " target=" << target
                                        << " value=" << suffix << '\n';
                                }
                            }
                        }
                        Integer direct_current = 0;
                        for (int slice = 0;
                             slice <= truncation;
                             ++slice) {
                            const int even = 2 * slice;
                            const int odd = even + 1;
                            const Integer weight =
                                binomial(2 * prefix + 1, even);
                            direct_current += weight * (
                                returns[
                                    static_cast<std::size_t>(even)
                                ] * endpoints[
                                    static_cast<std::size_t>(
                                        2 * prefix + 2 - even
                                    )
                                ]
                                - returns[
                                    static_cast<std::size_t>(odd)
                                ] * endpoints[
                                    static_cast<std::size_t>(
                                        2 * prefix + 2 - odd
                                    )
                                ]
                            );
                        }
                        ++endpoint_identity_rows;
                        if (suffix != 2 * direct_current) {
                            throw std::runtime_error(
                                "central-crossing sum disagrees "
                                "with endpoint current"
                            );
                        }
                    }
                }
            }
        }

        if (case_mode) {
            std::cout
                << "SU2_CENTRAL_CROSSING_CURRENT_CASE"
                << " level=" << case_level
                << " label=" << case_label
                << " prefix=" << case_prefix
                << " truncation=" << case_truncation
                << " doubled_total=" << case_total
                << " result=COMPONENT_PROFILE\n";
            return EXIT_SUCCESS;
        }

        std::cout
            << "SU2_CENTRAL_CROSSING_CURRENT"
            << " maximum_level=" << maximum_level
            << " maximum_prefix=" << maximum_prefix
            << " parameters=" << parameters
            << " quotient_entries=" << quotient_entries
            << " crossing_entries=" << crossing_entries
            << " component_rows=" << component_rows
            << " component_negatives=" << component_negatives
            << " minimum=" << minimum.value
            << " witness=("
            << minimum.level << ','
            << minimum.label << ','
            << minimum.prefix << ','
            << minimum.truncation << ','
            << minimum.crossing_source << ','
            << minimum.crossing_target << ')'
            << " target_rows=" << target_rows
            << " target_negatives=" << target_negatives
            << " target_sign_recrossings="
                << target_sign_recrossings
            << " target_sign_reversals="
                << target_sign_reversals
            << " maximum_target_sign_changes="
                << maximum_target_sign_changes
            << " target_minimum=" << target_minimum.value
            << " target_witness=("
            << target_minimum.level << ','
            << target_minimum.label << ','
            << target_minimum.prefix << ','
            << target_minimum.truncation << ','
            << target_minimum.crossing_target << ')'
            << " target_suffix_rows=" << target_suffix_rows
            << " target_suffix_negatives=" << target_suffix_negatives
            << " target_suffix_minimum=" << target_suffix_minimum.value
            << " target_suffix_witness=("
            << target_suffix_minimum.level << ','
            << target_suffix_minimum.label << ','
            << target_suffix_minimum.prefix << ','
            << target_suffix_minimum.truncation << ','
            << target_suffix_minimum.crossing_target << ')'
            << " adjacent_target_rows=" << adjacent_target_rows
            << " adjacent_target_negatives="
                << adjacent_target_negatives
            << " adjacent_target_minimum="
                << adjacent_target_minimum.value
            << " adjacent_target_witness=("
            << adjacent_target_minimum.level << ','
            << adjacent_target_minimum.label << ','
            << adjacent_target_minimum.prefix << ','
            << adjacent_target_minimum.truncation << ','
            << adjacent_target_minimum.crossing_source << ','
            << adjacent_target_minimum.crossing_target << ')'
            << " inner_reserve_rows=" << inner_reserve_rows
            << " inner_reserve_negatives=" << inner_reserve_negatives
            << " inner_reserve_minimum=" << inner_reserve_minimum.value
            << " inner_reserve_witness=("
            << inner_reserve_minimum.level << ','
            << inner_reserve_minimum.label << ','
            << inner_reserve_minimum.prefix << ','
            << inner_reserve_minimum.truncation << ','
            << inner_reserve_minimum.crossing_target << ')'
            << " endpoint_identity_rows=" << endpoint_identity_rows
            << " result="
            << (
                target_suffix_negatives == 0U
                ? "PASS_TARGET_SUFFIX_DISCOVERY"
                : "FAIL_TARGET_SUFFIX"
            )
            << '\n';
        return target_suffix_negatives == 0U
            ? EXIT_SUCCESS
            : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr
            << "SU2_CENTRAL_CROSSING_CURRENT FAILURE: "
            << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
