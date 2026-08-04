#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Vector = std::vector<Integer>;
using Matrix = std::vector<Vector>;

int positive(const char* text, const char* name) {
    const std::string value{text};
    std::size_t used = 0U;
    const long parsed = std::stol(value, &used);
    if (used != value.size() || parsed <= 0) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

Integer binomial(int upper, int lower) {
    if (lower < 0 || lower > upper) {
        return 0;
    }
    lower = std::min(lower, upper - lower);
    Integer result = 1;
    for (int index = 1; index <= lower; ++index) {
        result *= upper - lower + index;
        result /= index;
    }
    return result;
}

Matrix fusion_matrix(int half_level, int half_factor) {
    const int size = half_level + 1;
    Matrix result(static_cast<std::size_t>(size), Vector(
        static_cast<std::size_t>(size)
    ));
    for (int source = 0; source <= half_level; ++source) {
        const int lower = std::abs(source - half_factor);
        const int upper = std::min(
            source + half_factor,
            2 * half_level - source - half_factor
        );
        for (int target = lower; target <= upper; ++target) {
            result[static_cast<std::size_t>(target)]
                  [static_cast<std::size_t>(source)] = 1;
        }
    }
    return result;
}

Vector multiply(const Matrix& matrix, const Vector& vector) {
    const int size = static_cast<int>(matrix.size());
    Vector result(static_cast<std::size_t>(size));
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            const Integer& entry = matrix[static_cast<std::size_t>(row)]
                                         [static_cast<std::size_t>(column)];
            if (entry != 0 && vector[static_cast<std::size_t>(column)] != 0) {
                result[static_cast<std::size_t>(row)]
                    += entry * vector[static_cast<std::size_t>(column)];
            }
        }
    }
    return result;
}

Vector add(const Vector& left, const Vector& right) {
    Vector result = left;
    for (std::size_t index = 0U; index < result.size(); ++index) {
        result[index] += right[index];
    }
    return result;
}

Vector scale(const Vector& vector, const Integer& scalar) {
    Vector result = vector;
    for (Integer& entry : result) {
        entry *= scalar;
    }
    return result;
}

Vector reflected(const Vector& vector) {
    Vector result = vector;
    std::reverse(result.begin(), result.end());
    return result;
}

bool nonnegative(const Vector& vector) {
    return std::all_of(vector.begin(), vector.end(), [](const Integer& value) {
        return value >= 0;
    });
}

Integer coefficient(int depth, int order, int slice) {
    Integer result = 0;
    for (int increment = 0; increment <= order; ++increment) {
        const Integer term = binomial(order, increment)
            * binomial(2 * depth + 3 + 2 * increment, 2 * slice);
        result += ((order - increment) & 1) == 0 ? term : -term;
    }
    return result;
}

struct Powers {
    std::vector<Vector> values;
    std::vector<Integer> returns;
};

Powers powers(const Matrix& matrix, int maximum) {
    const int size = static_cast<int>(matrix.size());
    Powers result;
    result.values.reserve(static_cast<std::size_t>(maximum + 1));
    result.returns.reserve(static_cast<std::size_t>(maximum + 1));
    Vector current(static_cast<std::size_t>(size));
    current[0] = 1;
    for (int exponent = 0; exponent <= maximum; ++exponent) {
        result.returns.push_back(current[0]);
        result.values.push_back(current);
        if (exponent != maximum) {
            current = multiply(matrix, current);
        }
    }
    return result;
}

Vector low_group(const Powers& powers, int depth, int order) {
    const int size = static_cast<int>(powers.values.front().size());
    Vector result(static_cast<std::size_t>(size));
    for (int slice = 0; slice <= depth; ++slice) {
        const Integer weight = coefficient(depth, order, slice);
        if (weight == 0) {
            continue;
        }
        const int even = 2 * slice;
        const int first_power = 2 * depth + 4 + 2 * order - 2 * slice;
        const int second_power = first_power - 1;
        const Vector first = scale(
            powers.values[static_cast<std::size_t>(first_power)],
            powers.returns[static_cast<std::size_t>(even)] * weight
        );
        const Vector second = scale(
            powers.values[static_cast<std::size_t>(second_power)],
            -powers.returns[static_cast<std::size_t>(even + 1)] * weight
        );
        result = add(result, add(first, second));
    }
    return result;
}

Vector terminal_defect(const Powers& powers, int depth, int order) {
    const Integer weight = coefficient(depth, order, depth);
    const int first_power = 2 * order + 4;
    const int second_power = first_power - 1;
    return add(
        scale(
            powers.values[static_cast<std::size_t>(first_power)],
            powers.returns[static_cast<std::size_t>(2 * depth)] * weight
        ),
        scale(
            powers.values[static_cast<std::size_t>(second_power)],
            -powers.returns[static_cast<std::size_t>(2 * depth + 1)] * weight
        )
    );
}

Vector direct_base_value(const Powers& powers, int depth, int index) {
    const int size = static_cast<int>(powers.values.front().size());
    Vector result(static_cast<std::size_t>(size));
    for (int slice = 0; slice <= depth; ++slice) {
        const Integer weight = binomial(2 * index + 1, 2 * slice);
        if (weight == 0) {
            continue;
        }
        const int first_power = 2 * index + 2 - 2 * slice;
        const int second_power = first_power - 1;
        result = add(
            result,
            add(
                scale(
                    reflected(powers.values[static_cast<std::size_t>(first_power)]),
                    powers.returns[static_cast<std::size_t>(2 * slice)] * weight
                ),
                scale(
                    reflected(powers.values[static_cast<std::size_t>(second_power)]),
                    -powers.returns[static_cast<std::size_t>(2 * slice + 1)] * weight
                )
            )
        );
    }
    return result;
}

Vector direct_newton_kernel(const Matrix& matrix, const Powers& powers,
                             int depth, int order) {
    const int size = static_cast<int>(powers.values.front().size());
    Vector result(static_cast<std::size_t>(size));
    for (int increment = 0; increment <= order; ++increment) {
        Vector term = direct_base_value(powers, depth, depth + 1 + increment);
        for (int step = 0; step < order - increment; ++step) {
            term = multiply(matrix, multiply(matrix, term));
        }
        const Integer weight = binomial(order, increment);
        result = add(
            result,
            scale(term, ((order - increment) & 1) == 0 ? weight : -weight)
        );
    }
    return reflected(result);
}

std::string render(const Vector& vector) {
    std::string result{"["};
    for (std::size_t index = 0U; index < vector.size(); ++index) {
        if (index != 0U) {
            result += ',';
        }
        result += vector[index].convert_to<std::string>();
    }
    return result + ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        int argument = 1;
        bool scalar_diagonal_only = false;
        if (argc >= 2 && std::string(argv[1]) == "--scalar-diagonal") {
            scalar_diagonal_only = true;
            ++argument;
        }
        const int maximum_half_level = argc > argument
            ? positive(argv[argument], "maximum half level")
            : 20;
        const int maximum_depth = argc > argument + 1
            ? positive(argv[argument + 1], "maximum depth")
            : 10;
        if (argc > argument + 2) {
            throw std::invalid_argument(
                "usage: probe_su2_newton_depth_recursion "
                "[--scalar-diagonal] [maximum_half_level] [maximum_depth]"
            );
        }

        std::uint64_t recurrence_checks = 0U;
        std::uint64_t paired_checks = 0U;
        bool paired_counterexample = false;
        int paired_half_level = 0;
        int paired_half_factor = 0;
        int paired_depth = 0;
        int paired_order = 0;
        Vector paired_witness;
        Vector paired_current;
        for (int half_level = 3;
             half_level <= maximum_half_level;
             ++half_level) {
            for (int half_factor = 1;
                 2 * half_factor < half_level;
                 ++half_factor) {
                const Matrix matrix = fusion_matrix(half_level, half_factor);
                const int maximum_power = 4 * maximum_depth + 4;
                const Powers data = powers(matrix, maximum_power);
                for (int depth = 2; depth <= maximum_depth; ++depth) {
                    for (int order = 0; order <= depth - 2; ++order) {
                        const Vector current = low_group(data, depth, order);
                        const Vector direct = direct_newton_kernel(
                            matrix, data, depth, order
                        );
                        if (current != direct) {
                            throw std::runtime_error(
                                "closed low group disagrees with direct Newton kernel"
                            );
                        }
                        const Vector previous = low_group(
                            data, depth - 1, order
                        );
                        const Vector shifted = low_group(
                            data, depth - 1, order + 1
                        );
                        const Vector defect = terminal_defect(
                            data, depth, order
                        );
                        const Vector reconstructed = add(
                            add(multiply(matrix, multiply(matrix, previous)), shifted),
                            defect
                        );
                        ++recurrence_checks;
                        if (current != reconstructed) {
                            throw std::runtime_error(
                                "Newton depth recurrence mismatch"
                            );
                        }
                        Vector smoothed = current;
                        for (int smoothing = 0;
                             smoothing <= depth - 2 - order;
                             ++smoothing) {
                            if (smoothed[0] < 0) {
                                std::cout
                                    << "SU2_NEWTON_DEPTH_RECURSION"
                                    << " result=SCALAR_DIAGONAL_COUNTEREXAMPLE"
                                    << " half_level=" << half_level
                                    << " half_factor=" << half_factor
                                    << " depth=" << depth
                                    << " order=" << order
                                    << " smoothing=" << smoothing
                                    << " value=" << smoothed[0]
                                    << " current=" << render(current)
                                    << '\n';
                                return EXIT_SUCCESS;
                            }
                            if (smoothing != depth - 2 - order) {
                                smoothed = multiply(
                                    matrix, multiply(matrix, smoothed)
                                );
                            }
                        }
                        if (!scalar_diagonal_only && !nonnegative(current)) {
                            std::cout
                                << "SU2_NEWTON_DEPTH_RECURSION"
                                << " result=LOW_GROUP_COUNTEREXAMPLE"
                                << " half_level=" << half_level
                                << " half_factor=" << half_factor
                                << " depth=" << depth
                                << " order=" << order
                                << " current=" << render(current)
                                << '\n';
                            return EXIT_SUCCESS;
                        }
                        const Vector paired = add(shifted, defect);
                        ++paired_checks;
                        if (!paired_counterexample && !nonnegative(paired)) {
                            paired_counterexample = true;
                            paired_half_level = half_level;
                            paired_half_factor = half_factor;
                            paired_depth = depth;
                            paired_order = order;
                            paired_witness = paired;
                            paired_current = current;
                        }
                    }
                }
            }
        }
        std::cout << "SU2_NEWTON_DEPTH_RECURSION";
        if (scalar_diagonal_only) {
            std::cout << " result=NO_SCALAR_DIAGONAL_COUNTEREXAMPLE";
        } else if (paired_counterexample) {
            std::cout
                << " result=PAIRED_COUNTEREXAMPLE_LOW_GROUPS_NONNEGATIVE"
                << " half_level=" << paired_half_level
                << " half_factor=" << paired_half_factor
                << " depth=" << paired_depth
                << " order=" << paired_order
                << " paired=" << render(paired_witness)
                << " current=" << render(paired_current);
        } else {
            std::cout << " result=NO_PAIRED_COUNTEREXAMPLE";
        }
        std::cout
            << " maximum_half_level=" << maximum_half_level
            << " maximum_depth=" << maximum_depth
            << " recurrence_checks=" << recurrence_checks
            << " paired_checks=" << paired_checks << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_NEWTON_DEPTH_RECURSION error="
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
