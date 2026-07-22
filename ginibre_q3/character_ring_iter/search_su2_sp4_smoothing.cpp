#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

using Matrix = std::vector<std::vector<cpp_int>>;

void for_each_output(int left, int right, const auto& function) {
    for (int output = std::abs(left - right); output <= left + right; output += 2) {
        function(output);
    }
}

Matrix add_plus(const Matrix& current, int label) {
    const int dimension = static_cast<int>(current.size());
    Matrix next(
        static_cast<std::size_t>(dimension),
        std::vector<cpp_int>(static_cast<std::size_t>(dimension))
    );
    for (int left = 0; left < dimension; ++left) {
        for (int right = 0; right < dimension; ++right) {
            const cpp_int& value = current[static_cast<std::size_t>(left)]
                [static_cast<std::size_t>(right)];
            if (value == 0) {
                continue;
            }
            for_each_output(left, label, [&](int output) {
                if (output < dimension) {
                    next[static_cast<std::size_t>(output)]
                        [static_cast<std::size_t>(right)] += value;
                }
            });
            for_each_output(right, label, [&](int output) {
                if (output < dimension) {
                    next[static_cast<std::size_t>(left)]
                        [static_cast<std::size_t>(output)] += value;
                }
            });
        }
    }
    return next;
}

cpp_int at(const Matrix& matrix, int row, int column) {
    if (row < 0 || column < 0 || row >= static_cast<int>(matrix.size())
        || column >= static_cast<int>(matrix.size())) {
        return 0;
    }
    return matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)];
}

bool character_positive(
    const Matrix& coefficients,
    int& failure_first,
    int& failure_second,
    cpp_int& failure_value
) {
    for (int first = 0; first < static_cast<int>(coefficients.size()) - 1; ++first) {
        for (int second = 0; second <= first; ++second) {
            const int row = first + 1;
            const cpp_int value = at(coefficients, row - 1, second)
                + at(coefficients, row + 1, second)
                - at(coefficients, row, second - 1)
                - at(coefficients, row, second + 1);
            if (value < 0) {
                failure_first = first;
                failure_second = second;
                failure_value = value;
                return false;
            }
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: search_su2_sp4_smoothing MAXIMUM_P MAXIMUM_ONES"
            );
        }
        const int maximum_p = std::stoi(argv[1]);
        const int maximum_ones = std::stoi(argv[2]);
        if (maximum_p < 1 || maximum_ones < 0) {
            throw std::runtime_error("invalid search bound");
        }
        for (int p = 1; p <= maximum_p; ++p) {
            int first_positive = -1;
            bool positivity_persisted = true;
            int last_failure_first = -1;
            int last_failure_second = -1;
            cpp_int last_failure_value = 0;
            const int dimension = p + maximum_ones + 3;
            Matrix coefficients(
                static_cast<std::size_t>(dimension),
                std::vector<cpp_int>(static_cast<std::size_t>(dimension))
            );
            coefficients[0][0] = 1;
            coefficients = add_plus(coefficients, p);
            for (int ones = 0; ones <= maximum_ones; ++ones) {
                int failure_first = -1;
                int failure_second = -1;
                cpp_int failure_value = 0;
                if (character_positive(
                        coefficients, failure_first, failure_second, failure_value
                    )) {
                    if (first_positive < 0) {
                        first_positive = ones;
                    }
                } else {
                    if (first_positive >= 0) {
                        positivity_persisted = false;
                    }
                    last_failure_first = failure_first;
                    last_failure_second = failure_second;
                    last_failure_value = failure_value;
                }
                if (ones != maximum_ones) {
                    coefficients = add_plus(coefficients, 1);
                }
            }
            std::cout << "p=" << p << " first_character_positive_ones="
                      << first_positive
                      << " persisted=" << (positivity_persisted ? "yes" : "no");
            if (last_failure_first >= 0) {
                std::cout << " last_failure=(" << last_failure_first << ','
                          << last_failure_second << ") coefficient="
                          << last_failure_value;
            }
            std::cout << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
