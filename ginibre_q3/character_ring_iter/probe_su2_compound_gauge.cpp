#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Vector = std::vector<cpp_int>;
using Matrix = std::vector<Vector>;

int parse_positive(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0LL
        || parsed > static_cast<long long>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    }
    return static_cast<int>(parsed);
}

Matrix fusion_matrix(int level, int factor) {
    Matrix result(
        static_cast<std::size_t>(level + 1),
        Vector(static_cast<std::size_t>(level + 1), 0));
    for (int source = 0; source <= level; ++source) {
        const int lower = std::abs(source - factor);
        const int upper = std::min(source + factor, 2 * level - source - factor);
        for (int target = lower; target <= upper; ++target) {
            result[static_cast<std::size_t>(target)][
                static_cast<std::size_t>(source)] = 1;
        }
    }
    return result;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const std::size_t size = left.size();
    Matrix result(size, Vector(size, 0));
    for (std::size_t row = 0U; row < size; ++row) {
        for (std::size_t middle = 0U; middle < size; ++middle) {
            if (left[row][middle] == 0) {
                continue;
            }
            for (std::size_t column = 0U; column < size; ++column) {
                result[row][column]
                    += left[row][middle] * right[middle][column];
            }
        }
    }
    return result;
}

std::vector<std::pair<int, int>> pairs(int level) {
    std::vector<std::pair<int, int>> result;
    for (int left = 0; left <= level; ++left) {
        for (int right = left + 1; right <= level; ++right) {
            result.emplace_back(left, right);
        }
    }
    return result;
}

Matrix compound(const Matrix& matrix, const std::vector<std::pair<int, int>>& basis) {
    Matrix result(
        basis.size(),
        Vector(basis.size(), 0));
    for (std::size_t row = 0U; row < basis.size(); ++row) {
        const auto [row_left, row_right] = basis[row];
        for (std::size_t column = 0U; column < basis.size(); ++column) {
            const auto [column_left, column_right] = basis[column];
            result[row][column]
                = matrix[static_cast<std::size_t>(row_left)][
                    static_cast<std::size_t>(column_left)]
                    * matrix[static_cast<std::size_t>(row_right)][
                    static_cast<std::size_t>(column_right)]
                - matrix[static_cast<std::size_t>(row_left)][
                    static_cast<std::size_t>(column_right)]
                    * matrix[static_cast<std::size_t>(row_right)][
                    static_cast<std::size_t>(column_left)];
        }
    }
    return result;
}

int sign(const cpp_int& value) {
    return value < 0 ? -1 : (value > 0 ? 1 : 0);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_level = argc >= 2
            ? parse_positive(argv[1], "maximum_half_level")
            : 50;
        if (argc > 2) {
            throw std::invalid_argument(
                "usage: probe_su2_compound_gauge [maximum_half_level]");
        }

        std::uint64_t cases = 0U;
        std::uint64_t balanced_cases = 0U;
        std::uint64_t boundary_compatible_cases = 0U;
        int first_unbalanced_level = -1;
        int first_unbalanced_factor = -1;
        int first_boundary_conflict_level = -1;
        int first_boundary_conflict_factor = -1;
        std::pair<int, int> first_conflict_edge{-1, -1};
        std::string first_unbalanced_cycle = "none";

        for (int level = 1; level <= maximum_level; ++level) {
            const auto basis = pairs(level);
            for (int factor = 1; factor <= level; ++factor) {
                ++cases;
                const Matrix fusion = fusion_matrix(level, factor);
                const Matrix square = multiply(fusion, fusion);
                const Matrix lifted = compound(square, basis);

                std::vector<int> gauge(basis.size(), 0);
                std::vector<int> parent(basis.size(), -1);
                bool balanced = true;
                for (std::size_t start = 0U;
                     start < basis.size() && balanced;
                     ++start) {
                    if (gauge[start] != 0) {
                        continue;
                    }
                    gauge[start] = 1;
                    std::queue<std::size_t> pending;
                    pending.push(start);
                    while (!pending.empty() && balanced) {
                        const std::size_t left = pending.front();
                        pending.pop();
                        for (std::size_t right = 0U;
                             right < basis.size();
                             ++right) {
                            if (left == right || lifted[left][right] == 0) {
                                continue;
                            }
                            const int required
                                = gauge[left] * sign(lifted[left][right]);
                            if (gauge[right] == 0) {
                                gauge[right] = required;
                                parent[right] = static_cast<int>(left);
                                pending.push(right);
                            } else if (gauge[right] != required) {
                                balanced = false;
                                if (first_conflict_edge.first < 0) {
                                    first_conflict_edge = {
                                        static_cast<int>(left),
                                        static_cast<int>(right)};
                                    std::vector<std::size_t> left_path;
                                    std::vector<std::size_t> right_path;
                                    for (int vertex = static_cast<int>(left);
                                         vertex >= 0;
                                         vertex = parent[static_cast<std::size_t>(vertex)]) {
                                        left_path.push_back(
                                            static_cast<std::size_t>(vertex));
                                    }
                                    for (int vertex = static_cast<int>(right);
                                         vertex >= 0;
                                         vertex = parent[static_cast<std::size_t>(vertex)]) {
                                        right_path.push_back(
                                            static_cast<std::size_t>(vertex));
                                    }
                                    std::size_t left_tail = left_path.size();
                                    std::size_t right_tail = right_path.size();
                                    while (left_tail > 0U && right_tail > 0U
                                           && left_path[left_tail - 1U]
                                               == right_path[right_tail - 1U]) {
                                        --left_tail;
                                        --right_tail;
                                    }
                                    std::vector<std::size_t> cycle;
                                    for (std::size_t position = 0U;
                                         position <= left_tail;
                                         ++position) {
                                        cycle.push_back(left_path[position]);
                                    }
                                    for (std::size_t position = right_tail;
                                         position > 0U;
                                         --position) {
                                        cycle.push_back(right_path[position - 1U]);
                                    }
                                    cycle.push_back(left);
                                    std::ostringstream description;
                                    int cycle_sign = 1;
                                    {
                                        const auto [a, b] = basis[cycle.front()];
                                        description << '(' << a << ',' << b << ')';
                                    }
                                    for (std::size_t position = 0U;
                                         position + 1U < cycle.size();
                                         ++position) {
                                        const cpp_int& edge_value
                                            = lifted[cycle[position]][
                                                cycle[position + 1U]];
                                        const int edge_sign = sign(
                                            edge_value);
                                        cycle_sign *= edge_sign;
                                        const auto [a, b]
                                            = basis[cycle[position + 1U]];
                                        description << "-[" << edge_value << "]-"
                                                    << '(' << a << ',' << b << ')';
                                    }
                                    description << ";sign_product=" << cycle_sign;
                                    first_unbalanced_cycle = description.str();
                                }
                                break;
                            }
                        }
                    }
                }
                if (!balanced) {
                    if (first_unbalanced_level < 0) {
                        first_unbalanced_level = level;
                        first_unbalanced_factor = factor;
                    }
                    continue;
                }
                ++balanced_cases;

                bool boundary_compatible = true;
                int boundary_sign = 0;
                for (std::size_t index = 0U;
                     index < basis.size();
                     ++index) {
                    if (basis[index].first != 0) {
                        continue;
                    }
                    if (boundary_sign == 0) {
                        boundary_sign = gauge[index];
                    } else if (gauge[index] != boundary_sign) {
                        boundary_compatible = false;
                        break;
                    }
                }
                if (boundary_compatible) {
                    ++boundary_compatible_cases;
                } else if (first_boundary_conflict_level < 0) {
                    first_boundary_conflict_level = level;
                    first_boundary_conflict_factor = factor;
                }
            }
        }

        std::cout
            << "SU2_COMPOUND_GAUGE"
            << " maximum_level=" << 2 * maximum_level
            << " cases=" << cases
            << " balanced_cases=" << balanced_cases
            << " boundary_compatible_cases="
            << boundary_compatible_cases
            << " first_unbalanced_level="
            << (first_unbalanced_level < 0
                    ? -1
                    : 2 * first_unbalanced_level)
            << " first_unbalanced_factor="
            << (first_unbalanced_factor < 0
                    ? -1
                    : 2 * first_unbalanced_factor)
            << " first_boundary_conflict_level="
            << (first_boundary_conflict_level < 0
                    ? -1
                    : 2 * first_boundary_conflict_level)
            << " first_boundary_conflict_factor="
            << (first_boundary_conflict_factor < 0
                    ? -1
                    : 2 * first_boundary_conflict_factor)
            << " first_conflict_edge=("
            << first_conflict_edge.first << ','
            << first_conflict_edge.second << ')'
            << " first_unbalanced_cycle="
            << first_unbalanced_cycle
            << " result="
            << (boundary_compatible_cases == cases
                    ? "PASS_NONNEGATIVE_BOUNDARY_GAUGE"
                    : "FAIL_NONNEGATIVE_BOUNDARY_GAUGE")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
