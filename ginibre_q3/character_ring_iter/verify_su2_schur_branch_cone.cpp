#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Matrix = std::vector<std::vector<Integer>>;
using Vector = std::vector<Integer>;

struct Compound {
    Matrix matrix;
    std::vector<std::pair<int, int>> states;
    std::vector<std::vector<int>> index;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) {
        fail(message);
    }
}

int path_entry(int first, int second) {
    return std::abs(first - second) <= 1 ? 1 : 0;
}

Compound second_compound(int level) {
    Compound result;
    result.index.assign(
        static_cast<std::size_t>(level),
        std::vector<int>(static_cast<std::size_t>(level), -1)
    );
    for (int first = 0; first < level; ++first) {
        for (int second = first + 1; second < level; ++second) {
            result.index[static_cast<std::size_t>(first)]
                        [static_cast<std::size_t>(second)] =
                static_cast<int>(result.states.size());
            result.states.emplace_back(first, second);
        }
    }
    const std::size_t size = result.states.size();
    result.matrix.assign(size, Vector(size));
    for (std::size_t row = 0; row < size; ++row) {
        const auto [first, second] = result.states[row];
        for (std::size_t column = 0; column < size; ++column) {
            const auto [third, fourth] = result.states[column];
            result.matrix[row][column] =
                path_entry(first, third) * path_entry(second, fourth)
                - path_entry(first, fourth) * path_entry(second, third);
        }
    }
    return result;
}

Vector multiply(const Matrix& matrix, const Vector& vector) {
    require(matrix.size() == vector.size(), "matrix/vector size mismatch");
    Vector result(vector.size());
    for (std::size_t row = 0; row < matrix.size(); ++row) {
        require(
            matrix[row].size() == vector.size(),
            "matrix is not square"
        );
        for (std::size_t column = 0; column < vector.size(); ++column) {
            result[row] += matrix[row][column] * vector[column];
        }
    }
    return result;
}

Vector iterate(const Matrix& matrix, Vector vector, int exponent) {
    for (int step = 0; step < exponent; ++step) {
        vector = multiply(matrix, vector);
    }
    return vector;
}

Vector boundary_seed(const Compound& compound) {
    Vector seed(compound.states.size());
    const int adjacent_zero =
        compound.index[0][1];
    const int adjacent_one =
        compound.index[1][2];
    seed[static_cast<std::size_t>(adjacent_zero)] = -2;
    seed[static_cast<std::size_t>(adjacent_one)] = 1;
    return seed;
}

Integer boundary_value(int level, int exponent) {
    const Compound compound = second_compound(level);
    const Vector value =
        iterate(compound.matrix, boundary_seed(compound), exponent);
    const int target =
        compound.index[static_cast<std::size_t>(level - 2)]
                      [static_cast<std::size_t>(level - 1)];
    return value[static_cast<std::size_t>(target)];
}

void verify_table(int level, const std::vector<long long>& expected) {
    const Compound compound = second_compound(level);
    const Vector value = iterate(compound.matrix, boundary_seed(compound), 5);
    require(
        expected.size() == value.size(),
        "hard-coded five-step table has the wrong size"
    );
    for (std::size_t coordinate = 0;
         coordinate < value.size();
         ++coordinate) {
        require(
            value[coordinate] == expected[coordinate],
            "five-step table mismatch at level " + std::to_string(level)
        );
        require(
            value[coordinate] >= 0,
            "negative five-step coefficient at level "
                + std::to_string(level)
        );
    }
}

Integer fibonacci(int index) {
    Integer previous = 0;
    Integer current = 1;
    for (int step = 0; step < index; ++step) {
        const Integer next = previous + current;
        previous = current;
        current = next;
    }
    return previous;
}

Integer lucas(int index) {
    if (index == 0) {
        return 2;
    }
    Integer previous = 2;
    Integer current = 1;
    for (int step = 1; step < index; ++step) {
        const Integer next = previous + current;
        previous = current;
        current = next;
    }
    return current;
}

std::vector<Integer> level_four_endpoints(int maximum_excess) {
    const int maximum_degree = 3 + maximum_excess;
    std::vector<Integer> endpoints(
        static_cast<std::size_t>(maximum_degree + 1)
    );
    Vector state(4);
    state[0] = 1;
    for (int degree = 0; degree <= maximum_degree; ++degree) {
        endpoints[static_cast<std::size_t>(degree)] = state[3];
        Vector next(4);
        for (int label = 0; label < 4; ++label) {
            next[static_cast<std::size_t>(label)] +=
                state[static_cast<std::size_t>(label)];
            if (label > 0) {
                next[static_cast<std::size_t>(label - 1)] +=
                    state[static_cast<std::size_t>(label)];
            }
            if (label + 1 < 4) {
                next[static_cast<std::size_t>(label + 1)] +=
                    state[static_cast<std::size_t>(label)];
            }
        }
        state = std::move(next);
    }
    return endpoints;
}

}  // namespace

int main() {
    try {
        for (int level = 4; level <= 9; ++level) {
            const Compound compound = second_compound(level);
            for (const auto& row : compound.matrix) {
                for (const Integer& entry : row) {
                    require(
                        entry >= 0,
                        "second compound is not entrywise nonnegative"
                    );
                }
            }
        }

        verify_table(
            5,
            {
                38, 82, 150, 112,
                70, 192, 154,
                138, 126,
                42
            }
        );
        verify_table(
            6,
            {
                38, 82, 180, 195, 145,
                70, 238, 280, 216,
                182, 245, 199,
                105, 115,
                40
            }
        );
        verify_table(
            7,
            {
                38, 82, 180, 195, 175, 83,
                70, 238, 280, 262, 126,
                182, 245, 243, 119,
                105, 145, 77,
                55, 35,
                8
            }
        );
        verify_table(
            8,
            {
                38, 82, 180, 195, 175, 83, 30,
                70, 238, 280, 262, 126, 46,
                182, 245, 243, 119, 44,
                105, 145, 77, 30,
                55, 35, 15,
                8, 5,
                1
            }
        );

        const Compound level_nine = second_compound(9);
        const Vector level_nine_value =
            iterate(level_nine.matrix, boundary_seed(level_nine), 5);
        const Compound level_eight = second_compound(8);
        const Vector level_eight_value =
            iterate(level_eight.matrix, boundary_seed(level_eight), 5);
        for (std::size_t coordinate = 0;
             coordinate < level_nine.states.size();
             ++coordinate) {
            const auto [first, second] = level_nine.states[coordinate];
            if (second < 8) {
                const int embedded =
                    level_eight.index[static_cast<std::size_t>(first)]
                                     [static_cast<std::size_t>(second)];
                require(
                    level_nine_value[coordinate]
                        == level_eight_value[
                            static_cast<std::size_t>(embedded)
                        ],
                    "stable five-step table does not embed"
                );
            } else {
                require(
                    level_nine_value[coordinate] == 0,
                    "five-step support escaped the stable table"
                );
            }
        }

        require(boundary_value(5, 3) == 1, "wrong K=5,N=3 value");
        require(boundary_value(5, 4) == 14, "wrong K=5,N=4 value");
        require(boundary_value(6, 4) == 4, "wrong K=6,N=4 value");

        const Compound level_four = second_compound(4);
        const Vector seed = boundary_seed(level_four);
        const Vector first = iterate(level_four.matrix, seed, 1);
        const Vector third = iterate(level_four.matrix, seed, 3);
        const Vector fifth = iterate(level_four.matrix, seed, 5);
        Vector residual(seed.size());
        for (std::size_t coordinate = 0;
             coordinate < seed.size();
             ++coordinate) {
            residual[coordinate] =
                fifth[coordinate] - 3 * third[coordinate] + first[coordinate];
        }
        Vector expected_residual(seed.size());
        expected_residual[
            static_cast<std::size_t>(level_four.index[0][2])
        ] = 1;
        expected_residual[
            static_cast<std::size_t>(level_four.index[1][3])
        ] = -1;
        require(
            residual == expected_residual,
            "wrong K=4 odd-time recurrence residual"
        );
        require(
            iterate(level_four.matrix, residual, 2) == residual,
            "K=4 recurrence residual is not fixed by H^2"
        );
        require(boundary_value(4, 1) == 1, "wrong K=4,N=1 value");
        require(boundary_value(4, 3) == 3, "wrong K=4,N=3 value");

        constexpr int formula_checks = 65;
        const std::vector<Integer> endpoints =
            level_four_endpoints(2 * formula_checks + 4);
        for (int index = 0; index < formula_checks; ++index) {
            const int branch = 2 * index + 1;
            const int offset = 3;
            const Integer shifted =
                endpoints[
                    static_cast<std::size_t>(offset + branch + 1)
                ] * (
                    endpoints[
                        static_cast<std::size_t>(offset + branch)
                    ]
                    - endpoints[
                        static_cast<std::size_t>(offset + branch - 1)
                    ]
                )
                - endpoints[
                    static_cast<std::size_t>(offset + branch + 2)
                ] * (
                    endpoints[
                        static_cast<std::size_t>(offset + branch - 1)
                    ]
                    - endpoints[
                        static_cast<std::size_t>(offset + branch - 2)
                    ]
                );
            const Integer closed_numerator =
                38 * lucas(6 * index)
                + 85 * fibonacci(6 * index)
                - 48 * lucas(2 * index)
                - 100 * fibonacci(2 * index)
                + 40;
            require(
                20 * shifted == closed_numerator,
                "K=4 shifted Fibonacci-Lucas formula mismatch"
            );
            require(shifted > 0, "K=4 shifted branch is not positive");
        }

        std::cout
            << "SU2_SCHUR_BRANCH_CONE"
            << " local_levels=5,6,7,stable_8"
            << " k4_formula_checks=" << formula_checks
            << " result=PASS_UNBOUNDED_LOCAL_IDENTITIES"
            << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
