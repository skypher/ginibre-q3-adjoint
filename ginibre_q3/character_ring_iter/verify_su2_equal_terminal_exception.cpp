#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Cubic = std::array<Integer, 4U>;
using Tensor = std::array<Integer, 16U>;

std::vector<Integer> multiply(
    const std::vector<Integer>& first,
    const std::vector<Integer>& second
) {
    std::vector<Integer> product(first.size() + second.size() - 1U, 0);
    for (std::size_t i = 0U; i < first.size(); ++i) {
        for (std::size_t j = 0U; j < second.size(); ++j) {
            product[i + j] += first[i] * second[j];
        }
    }
    return product;
}

Cubic cubic_product(const std::array<int, 3U>& slopes) {
    std::vector<Integer> polynomial{1};
    for (int slope : slopes) {
        polynomial = multiply(polynomial, {1, -slope});
    }
    if (polynomial.size() != 4U) {
        throw std::runtime_error("cubic product has the wrong degree");
    }
    return {polynomial[0U], polynomial[1U],
            polynomial[2U], polynomial[3U]};
}

Tensor symmetric_tensor(const Cubic& first, const Cubic& second) {
    Tensor tensor{};
    for (std::size_t i = 0U; i < first.size(); ++i) {
        for (std::size_t j = 0U; j < second.size(); ++j) {
            tensor[i * second.size() + j]
                = first[i] * second[j] + second[i] * first[j];
        }
    }
    return tensor;
}

Integer determinant_bareiss(std::vector<std::vector<Integer>> matrix) {
    const std::size_t size = matrix.size();
    if (size == 0U) {
        return 1;
    }
    Integer sign = 1;
    Integer previous = 1;
    for (std::size_t pivot = 0U; pivot + 1U < size; ++pivot) {
        std::size_t row = pivot;
        while (row < size && matrix[row][pivot] == 0) {
            ++row;
        }
        if (row == size) {
            return 0;
        }
        if (row != pivot) {
            std::swap(matrix[row], matrix[pivot]);
            sign = -sign;
        }
        const Integer pivot_value = matrix[pivot][pivot];
        for (std::size_t i = pivot + 1U; i < size; ++i) {
            for (std::size_t j = pivot + 1U; j < size; ++j) {
                const Integer numerator
                    = matrix[i][j] * pivot_value
                    - matrix[i][pivot] * matrix[pivot][j];
                if (numerator % previous != 0) {
                    throw std::runtime_error("nonexact Bareiss division");
                }
                matrix[i][j] = numerator / previous;
            }
        }
        previous = pivot_value;
    }
    return sign * matrix[size - 1U][size - 1U];
}

}  // namespace

int main() {
    try {
        // Specialize the four generic roots (v,a,b,c) to the distinct
        // affine slopes (1,2,3,4).  The five distinct symmetric terminal
        // tensors for (q;e1,e2,e3)=(3;2,4,4) are
        //
        //   (v a^2) odot (v b c), (v a^2) odot (b c^2),
        //   b^3 odot (v c^2), b^3 odot c^3,
        //   c^3 odot (v b^2).
        constexpr int v = 1;
        constexpr int a = 2;
        constexpr int b = 3;
        constexpr int c = 4;
        const Cubic vaa = cubic_product({v, a, a});
        const Cubic vbc = cubic_product({v, b, c});
        const Cubic bcc = cubic_product({b, c, c});
        const Cubic bbb = cubic_product({b, b, b});
        const Cubic vcc = cubic_product({v, c, c});
        const Cubic ccc = cubic_product({c, c, c});
        const Cubic vbb = cubic_product({v, b, b});
        const std::array<Tensor, 5U> rows{
            symmetric_tensor(vaa, vbc),
            symmetric_tensor(vaa, bcc),
            symmetric_tensor(bbb, vcc),
            symmetric_tensor(bbb, ccc),
            symmetric_tensor(ccc, vbb)
        };

        std::array<std::size_t, 5U> columns{};
        bool found = false;
        Integer determinant = 0;
        for (std::size_t c0 = 0U; c0 < 12U && !found; ++c0) {
            for (std::size_t c1 = c0 + 1U; c1 < 13U && !found; ++c1) {
                for (std::size_t c2 = c1 + 1U; c2 < 14U && !found; ++c2) {
                    for (std::size_t c3 = c2 + 1U; c3 < 15U && !found; ++c3) {
                        for (std::size_t c4 = c3 + 1U; c4 < 16U; ++c4) {
                            const std::array<std::size_t, 5U> candidate{
                                c0, c1, c2, c3, c4
                            };
                            std::vector<std::vector<Integer>> minor(
                                rows.size(),
                                std::vector<Integer>(candidate.size(), 0)
                            );
                            for (std::size_t row = 0U; row < rows.size(); ++row) {
                                for (std::size_t column = 0U;
                                     column < candidate.size(); ++column) {
                                    minor[row][column]
                                        = rows[row][candidate[column]];
                                }
                            }
                            determinant = determinant_bareiss(std::move(minor));
                            if (determinant != 0) {
                                columns = candidate;
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
        if (!found) {
            std::cerr << "SU2_EQUAL_TERMINAL_EXCEPTION result=rank_defect\n";
            return EXIT_FAILURE;
        }
        std::cout << "SU2_EQUAL_TERMINAL_EXCEPTION slopes=1,2,3,4"
                  << " columns=";
        for (std::size_t index = 0U; index < columns.size(); ++index) {
            if (index != 0U) {
                std::cout << ',';
            }
            std::cout << columns[index];
        }
        std::cout << " determinant=" << determinant
                  << " result=independent\n";

        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
