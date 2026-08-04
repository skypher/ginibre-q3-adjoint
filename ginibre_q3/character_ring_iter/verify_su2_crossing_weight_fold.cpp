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

int parse_positive(const char* text) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (
        end == text
        || *end != '\0'
        || value <= 0
        || value > std::numeric_limits<int>::max()
    ) {
        throw std::runtime_error("bound must be a positive integer");
    }
    return static_cast<int>(value);
}

bool fuses(int level, int label, int source, int target) {
    return std::abs(source - label) <= target
        && target <= source + label
        && source + target + label <= 2 * level;
}

std::vector<Integer> multiply_row(
    const std::vector<Integer>& row,
    const Matrix& matrix
) {
    std::vector<Integer> result(matrix.size());
    for (std::size_t source = 0U; source < matrix.size(); ++source) {
        if (row[source] == 0) {
            continue;
        }
        for (std::size_t target = 0U; target < matrix.size(); ++target) {
            const int entry = matrix[source][target];
            if (entry != 0) {
                result[target] += row[source] * entry;
            }
        }
    }
    return result;
}

Integer binomial_integer(int top, int bottom) {
    if (bottom < 0 || top < bottom) {
        return 0;
    }
    Integer result = 1;
    for (int index = 1; index <= bottom; ++index) {
        result *= top - bottom + index;
        result /= index;
    }
    return result;
}

Integer tensor_weight(int power, int label, int index) {
    Integer result = 0;
    for (int summand = 0; summand <= power; ++summand) {
        const Integer term = binomial_integer(power, summand)
            * binomial_integer(
                index - summand * (label + 1) + power - 1,
                power - 1
            );
        if (summand % 2 == 0) {
            result += term;
        } else {
            result -= term;
        }
    }
    return result;
}

Integer folded_tail_endpoint(
    int level,
    int label,
    int power,
    int target
) {
    if (power == 0) {
        return 0;
    }
    const int product_label = power * label;
    const int lower_depth = level - label - target;
    const int upper_depth = label + target;
    if (lower_depth > upper_depth) {
        return 0;
    }
    Integer result = 0;
    for (int branch = 0; branch < 3; ++branch) {
        int lower_label = 0;
        int upper_label = 0;
        int sign = 1;
        if (branch == 0) {
            lower_label = lower_depth;
            upper_label = upper_depth;
        } else if (branch == 1) {
            lower_label = 2 * level - label - target + 1;
            upper_label = level + label + target + 1;
            sign = -1;
        } else {
            lower_label = 3 * level - label - target + 2;
            upper_label = 2 * level + label + target + 2;
        }
        const int lower_index = std::max(
            0,
            product_label - upper_label
        );
        const int upper_index = std::min(
            product_label,
            product_label - lower_label
        );
        if (lower_index > upper_index) {
            continue;
        }
        const Integer contribution = tensor_weight(
            power,
            2 * label,
            upper_index
        ) - tensor_weight(power, 2 * label, lower_index - 1);
        result += sign * contribution;
    }
    return 2 * result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::runtime_error("usage: MAXIMUM_LEVEL");
        }
        const int maximum_level = parse_positive(argv[1]);
        std::uint64_t parameters = 0U;
        std::uint64_t rows = 0U;
        std::uint64_t endpoint_rows = 0U;
        std::uint64_t entries = 0U;
        std::uint64_t central_entries = 0U;

        for (int level = 3; level <= maximum_level; ++level) {
            for (int label = 1; 2 * label < level; ++label) {
                ++parameters;
                const int paired = (level + 1) / 2;
                const bool has_center = level % 2 == 0;
                const int quotient_size = paired + (has_center ? 1 : 0);
                const int center = has_center ? paired : -1;

                Matrix full(
                    static_cast<std::size_t>(level + 1),
                    std::vector<int>(static_cast<std::size_t>(level + 1))
                );
                for (int source = 0; source <= level; ++source) {
                    for (int target = 0; target <= level; ++target) {
                        full[static_cast<std::size_t>(source)]
                            [static_cast<std::size_t>(target)] =
                            fuses(level, label, source, target) ? 1 : 0;
                    }
                }

                Matrix plus(
                    static_cast<std::size_t>(quotient_size),
                    std::vector<int>(
                        static_cast<std::size_t>(quotient_size)
                    )
                );
                Matrix delta(
                    static_cast<std::size_t>(quotient_size),
                    std::vector<int>(static_cast<std::size_t>(paired))
                );
                for (int source = 0; source < paired; ++source) {
                    for (int target = 0; target < paired; ++target) {
                        const int same = fuses(
                            level,
                            label,
                            source,
                            target
                        ) ? 1 : 0;
                        const int crossed = fuses(
                            level,
                            label,
                            source,
                            level - target
                        ) ? 1 : 0;
                        plus[static_cast<std::size_t>(source)]
                            [static_cast<std::size_t>(target)] =
                            same + crossed;
                        delta[static_cast<std::size_t>(source)]
                             [static_cast<std::size_t>(target)] =
                            2 * crossed;
                    }
                    if (has_center) {
                        const int joins = fuses(
                            level,
                            label,
                            source,
                            level / 2
                        ) ? 1 : 0;
                        plus[static_cast<std::size_t>(source)]
                            [static_cast<std::size_t>(center)] = joins;
                        plus[static_cast<std::size_t>(center)]
                            [static_cast<std::size_t>(source)] = 2 * joins;
                        delta[static_cast<std::size_t>(center)]
                             [static_cast<std::size_t>(source)] =
                            2 * joins;
                    }
                }
                if (has_center) {
                    plus[static_cast<std::size_t>(center)]
                        [static_cast<std::size_t>(center)] =
                        fuses(
                            level,
                            label,
                            level / 2,
                            level / 2
                        ) ? 1 : 0;
                }

                std::vector<Integer> quotient_row(
                    static_cast<std::size_t>(quotient_size)
                );
                std::vector<Integer> full_row(
                    static_cast<std::size_t>(level + 1)
                );
                quotient_row[0] = 1;
                full_row[0] = 1;
                for (int power = 0; power <= 5; ++power) {
                    for (int source = 0; source < quotient_size; ++source) {
                        Integer folded = 0;
                        if (source < paired) {
                            folded = full_row[
                                static_cast<std::size_t>(source)
                            ] + full_row[
                                static_cast<std::size_t>(level - source)
                            ];
                        } else {
                            folded = full_row[
                                static_cast<std::size_t>(level / 2)
                            ];
                        }
                        ++entries;
                        if (folded != quotient_row[
                                           static_cast<std::size_t>(source)
                                       ]) {
                            throw std::runtime_error(
                                "reflection-even fold mismatch"
                            );
                        }
                    }

                    for (int target = 0; target < paired; ++target) {
                        Integer direct = 0;
                        for (int source = 0;
                             source < quotient_size;
                             ++source) {
                            direct += quotient_row[
                                static_cast<std::size_t>(source)
                            ] * delta[static_cast<std::size_t>(source)]
                                     [static_cast<std::size_t>(target)];
                        }
                        Integer folded_tail = 0;
                        const int lower = level - label - target;
                        for (int source = 0; source < paired; ++source) {
                            if (source >= lower) {
                                folded_tail += 2 * (
                                    full_row[
                                        static_cast<std::size_t>(source)
                                    ] + full_row[
                                        static_cast<std::size_t>(
                                            level - source
                                        )]
                                );
                            }
                        }
                        if (has_center && center >= lower) {
                            folded_tail += 2 * full_row[
                                static_cast<std::size_t>(level / 2)
                            ];
                            ++central_entries;
                        }
                        ++rows;
                        if (direct != folded_tail) {
                            throw std::runtime_error(
                                "crossing-weight folded-tail mismatch"
                            );
                        }
                        if (power != 0) {
                            const Integer endpoint = folded_tail_endpoint(
                                level,
                                label,
                                power,
                                target
                            );
                            ++endpoint_rows;
                            if (direct != endpoint) {
                                throw std::runtime_error(
                                    "crossing-weight endpoint mismatch"
                                    " level=" + std::to_string(level)
                                    + " label=" + std::to_string(label)
                                    + " power=" + std::to_string(power)
                                    + " target=" + std::to_string(target)
                                    + " direct=" + direct.str()
                                    + " endpoint=" + endpoint.str()
                                );
                            }
                        }
                    }
                    if (power != 5) {
                        quotient_row = multiply_row(quotient_row, plus);
                        full_row = multiply_row(full_row, full);
                    }
                }
            }
        }

        std::cout
            << "SU2_CROSSING_WEIGHT_FOLD"
            << " maximum_level=" << maximum_level
            << " parameters=" << parameters
            << " fold_entries=" << entries
            << " tail_rows=" << rows
            << " endpoint_rows=" << endpoint_rows
            << " central_tail_rows=" << central_entries
            << " powers=0..5"
            << " result=PASS_EXACT_IDENTITIES\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
