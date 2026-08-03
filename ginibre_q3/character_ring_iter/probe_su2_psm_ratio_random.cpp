#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Vector = std::vector<Integer>;
using Matrix = std::vector<Vector>;

std::uint64_t parse_positive_u64(const char* text, const char* name) {
    const std::string input(text);
    std::size_t used = 0U;
    const unsigned long long value = std::stoull(input, &used, 10);
    if (used != input.size() || value == 0ULL) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<std::uint64_t>(value);
}

int parse_positive_int(const char* text, const char* name) {
    const std::uint64_t value = parse_positive_u64(text, name);
    if (value > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
}

std::uint64_t splitmix64(std::uint64_t& state) {
    state += UINT64_C(0x9e3779b97f4a7c15);
    std::uint64_t value = state;
    value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

Matrix fusion_matrix(int level, int factor) {
    const int dimension = level + 1;
    Matrix result(
        static_cast<std::size_t>(dimension),
        Vector(static_cast<std::size_t>(dimension), 0)
    );
    for (int source = 0; source <= level; ++source) {
        const int lower = std::abs(source - factor);
        const int upper = std::min(
            source + factor, 2 * level - source - factor
        );
        for (int target = lower; target <= upper; target += 2) {
            result[static_cast<std::size_t>(target)]
                  [static_cast<std::size_t>(source)] = 1;
        }
    }
    return result;
}

Vector apply_matrix(const Matrix& matrix, const Vector& input) {
    Vector result(matrix.size(), 0);
    for (std::size_t row = 0U; row < matrix.size(); ++row) {
        for (std::size_t column = 0U; column < input.size(); ++column) {
            result[row] += matrix[row][column] * input[column];
        }
    }
    return result;
}

Vector fusion_square(const Vector& root, const std::vector<Matrix>& fusion) {
    Vector result(root.size(), 0);
    for (std::size_t label = 0U; label < root.size(); ++label) {
        if (root[label] == 0) {
            continue;
        }
        const Vector image = apply_matrix(fusion[label], root);
        for (std::size_t target = 0U; target < result.size(); ++target) {
            result[target] += root[label] * image[target];
        }
    }
    return result;
}

bool compressed_log_concave(const Vector& values) {
    int parity = -1;
    int first = -1;
    int last = -1;
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (values[index] == 0) {
            continue;
        }
        const int current = static_cast<int>(index & 1U);
        if (parity < 0) {
            parity = current;
        } else if (parity != current) {
            return false;
        }
        if (first < 0) {
            first = static_cast<int>(index);
        }
        last = static_cast<int>(index);
    }
    if (first < 0) {
        return false;
    }
    for (int index = first; index <= last; index += 2) {
        if (values[static_cast<std::size_t>(index)] == 0) {
            return false;
        }
    }
    for (int index = first + 2; index + 2 <= last; index += 2) {
        const Integer& left = values[static_cast<std::size_t>(index - 2)];
        const Integer& middle = values[static_cast<std::size_t>(index)];
        const Integer& right = values[static_cast<std::size_t>(index + 2)];
        if (middle * middle < left * right) {
            return false;
        }
    }
    return true;
}

bool boundary_admissible(const Vector& square) {
    const int level = static_cast<int>(square.size()) - 1;
    for (int radius = 0; radius <= level; ++radius) {
        if (square[0U] * square[static_cast<std::size_t>(level - radius)]
            < square[static_cast<std::size_t>(radius)] * square.back()) {
            return false;
        }
    }
    return true;
}

bool full_current_admissible(
    const Vector& square,
    const std::vector<Matrix>& fusion,
    int& left,
    int& right
) {
    const Integer& anchor = square[0U];
    for (std::size_t row = 0U; row < fusion.size(); ++row) {
        const Vector translated = apply_matrix(fusion[row], square);
        for (std::size_t column = 0U; column < fusion.size(); ++column) {
            if (anchor * translated[column]
                < square[row] * square[column]) {
                left = static_cast<int>(row);
                right = static_cast<int>(column);
                return false;
            }
        }
    }
    return true;
}

Integer power(Integer base, int exponent) {
    Integer result = 1;
    while (exponent > 0) {
        if ((exponent & 1) != 0) {
            result *= base;
        }
        base *= base;
        exponent /= 2;
    }
    return result;
}

Vector random_parity_log_concave_root(
    std::uint64_t& state, int level, int base, int maximum_step
) {
    int parity = static_cast<int>(splitmix64(state) & UINT64_C(1));
    if (parity > level) {
        parity = 0;
    }
    const int count = (level - parity) / 2 + 1;
    const int first = static_cast<int>(
        splitmix64(state) % static_cast<std::uint64_t>(count)
    );
    const int last = first + static_cast<int>(
        splitmix64(state)
        % static_cast<std::uint64_t>(count - first)
    );
    std::vector<int> steps;
    steps.reserve(static_cast<std::size_t>(last - first));
    for (int index = first; index < last; ++index) {
        const int span = 2 * maximum_step + 1;
        steps.push_back(
            static_cast<int>(splitmix64(state)
                % static_cast<std::uint64_t>(span)) - maximum_step
        );
    }
    std::sort(steps.begin(), steps.end(), std::greater<int>());
    std::vector<int> exponent(static_cast<std::size_t>(last - first + 1), 0);
    for (std::size_t index = 1U; index < exponent.size(); ++index) {
        exponent[index] = exponent[index - 1U] + steps[index - 1U];
    }
    const int minimum = *std::min_element(exponent.begin(), exponent.end());
    Vector root(static_cast<std::size_t>(level + 1), 0);
    const Integer integer_base(base);
    for (std::size_t index = 0U; index < exponent.size(); ++index) {
        const int label = parity + 2 * (first + static_cast<int>(index));
        root[static_cast<std::size_t>(label)] = power(
            integer_base, exponent[index] - minimum
        );
    }
    if (!compressed_log_concave(root)) {
        throw std::runtime_error("generated root is not log concave");
    }
    return root;
}

void print_vector(const Vector& values) {
    std::cout << '[';
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5) {
            throw std::runtime_error(
                "usage: probe_su2_psm_ratio_random SAMPLES MAXIMUM_LEVEL "
                "BASE MAXIMUM_STEP"
            );
        }
        const std::uint64_t samples = parse_positive_u64(argv[1], "samples");
        const int maximum_level = parse_positive_int(argv[2], "maximum level");
        const int base = parse_positive_int(argv[3], "base");
        const int maximum_step = parse_positive_int(argv[4], "maximum step");
        if (base <= 1) {
            throw std::runtime_error("base must exceed one");
        }
        std::vector<std::vector<Matrix>> fusion_by_level(
            static_cast<std::size_t>(maximum_level + 1)
        );
        for (int level = 1; level <= maximum_level; ++level) {
            auto& fusion = fusion_by_level[static_cast<std::size_t>(level)];
            fusion.reserve(static_cast<std::size_t>(level + 1));
            for (int label = 0; label <= level; ++label) {
                fusion.push_back(fusion_matrix(level, label));
            }
        }
        std::uint64_t state = UINT64_C(0x7446c6b151bf75d1);
        std::uint64_t square_shape = 0U;
        std::uint64_t boundary_profiles = 0U;
        std::uint64_t insertion_tests = 0U;
        for (std::uint64_t sample = 0U; sample < samples; ++sample) {
            const int level = 1 + static_cast<int>(
                splitmix64(state) % static_cast<std::uint64_t>(maximum_level)
            );
            const auto& fusion
                = fusion_by_level[static_cast<std::size_t>(level)];
            const Vector root = random_parity_log_concave_root(
                state, level, base, maximum_step
            );
            const Vector square = fusion_square(root, fusion);
            if (!compressed_log_concave(square)) {
                continue;
            }
            ++square_shape;
            if (!boundary_admissible(square)) {
                continue;
            }
            ++boundary_profiles;
            int left = -1;
            int right = -1;
            if (!full_current_admissible(square, fusion, left, right)) {
                std::cout << "SU2_PSM_RATIO_RANDOM result=STAR_COUNTEREXAMPLE"
                          << " sample=" << sample << " level=" << level
                          << " left=" << left << " right=" << right
                          << " root=";
                print_vector(root);
                std::cout << " square=";
                print_vector(square);
                std::cout << '\n';
                return EXIT_FAILURE;
            }
            for (int factor = 1; factor <= level; ++factor) {
                ++insertion_tests;
                const Vector image = apply_matrix(
                    fusion[static_cast<std::size_t>(factor)], root
                );
                const Vector image_square = fusion_square(image, fusion);
                if (!boundary_admissible(image_square)) {
                    std::cout << "SU2_PSM_RATIO_RANDOM result=INSERTION_COUNTEREXAMPLE"
                              << " sample=" << sample << " level=" << level
                              << " factor=" << factor << " root=";
                    print_vector(root);
                    std::cout << " image=";
                    print_vector(image);
                    std::cout << " square=";
                    print_vector(square);
                    std::cout << " image_square=";
                    print_vector(image_square);
                    std::cout << '\n';
                    return EXIT_FAILURE;
                }
            }
        }
        std::cout << "SU2_PSM_RATIO_RANDOM"
                  << " samples=" << samples
                  << " maximum_level=" << maximum_level
                  << " base=" << base
                  << " maximum_step=" << maximum_step
                  << " square_shape=" << square_shape
                  << " boundary_profiles=" << boundary_profiles
                  << " insertion_tests=" << insertion_tests
                  << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
