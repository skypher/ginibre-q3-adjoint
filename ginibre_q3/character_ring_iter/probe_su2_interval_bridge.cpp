#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Vector = std::vector<Integer>;
using Matrix = std::vector<Vector>;

int parse_positive(const char* text, const char* name) {
    const std::string input(text);
    std::size_t used = 0U;
    const long long value = std::stoll(input, &used, 10);
    if (used != input.size() || value <= 0LL) {
        throw std::runtime_error(std::string("invalid ") + name);
    }
    return static_cast<int>(value);
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
    Vector output(matrix.size(), 0);
    for (std::size_t row = 0U; row < matrix.size(); ++row) {
        for (std::size_t column = 0U; column < input.size(); ++column) {
            output[row] += matrix[row][column] * input[column];
        }
    }
    return output;
}

Vector square(const Vector& root, const std::vector<Matrix>& fusion) {
    Vector output(root.size(), 0);
    for (std::size_t label = 0U; label < root.size(); ++label) {
        if (root[label] == 0) {
            continue;
        }
        const Vector image = apply_matrix(fusion[label], root);
        for (std::size_t target = 0U; target < output.size(); ++target) {
            output[target] += root[label] * image[target];
        }
    }
    return output;
}

int boundary_failure(const Vector& profile) {
    const int level = static_cast<int>(profile.size()) - 1;
    for (int radius = 0; radius <= level; ++radius) {
        if (profile[0U] * profile[static_cast<std::size_t>(level - radius)]
            < profile[static_cast<std::size_t>(radius)]
                * profile[static_cast<std::size_t>(level)]) {
            return radius;
        }
    }
    return -1;
}

Integer boundary_margin(const Vector& profile, int radius) {
    const int level = static_cast<int>(profile.size()) - 1;
    return profile[0U] * profile[static_cast<std::size_t>(level - radius)]
        - profile[static_cast<std::size_t>(radius)] * profile.back();
}

std::pair<int, int> current_failure(
    const Vector& profile,
    const std::vector<Matrix>& fusion
) {
    const Integer anchor = profile[0U];
    for (std::size_t left = 0U; left < fusion.size(); ++left) {
        const Vector translated = apply_matrix(fusion[left], profile);
        for (std::size_t right = 0U; right < fusion.size(); ++right) {
            if (anchor * translated[right] < profile[left] * profile[right]) {
                return {
                    static_cast<int>(left), static_cast<int>(right)
                };
            }
        }
    }
    return {-1, -1};
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

Vector bridged_root(int scale, const std::string& family) {
    constexpr int level = 10;
    Vector root(static_cast<std::size_t>(level + 1), 0);
    if (family == "valley") {
        root[2U] = scale;
        root[4U] = 1;
        root[6U] = 1;
        root[8U] = 1;
        root[10U] = 2 * scale;
    } else if (family == "down") {
        root[2U] = scale;
        root[4U] = 1;
        root[6U] = 1;
        root[8U] = 1;
        root[10U] = 1;
    } else if (family == "up") {
        root[2U] = 1;
        root[4U] = 1;
        root[6U] = 1;
        root[8U] = 1;
        root[10U] = scale;
    } else if (family == "plateau") {
        root[2U] = 1;
        root[4U] = scale;
        root[6U] = scale;
        root[8U] = scale;
        root[10U] = 1;
    } else {
        throw std::runtime_error("invalid family");
    }
    return root;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2 && argc != 3) {
            throw std::runtime_error(
                "usage: probe_su2_interval_bridge MAXIMUM_SCALE "
                "[valley|down|up|plateau]"
            );
        }
        const int maximum_scale = parse_positive(argv[1], "maximum scale");
        const std::string family = argc == 3 ? argv[2] : "valley";
        constexpr int level = 10;
        std::vector<Matrix> fusion;
        fusion.reserve(static_cast<std::size_t>(level + 1));
        for (int label = 0; label <= level; ++label) {
            fusion.push_back(fusion_matrix(level, label));
        }
        int boundary_admissible = 0;
        for (int scale = 1; scale <= maximum_scale; ++scale) {
            const Vector root = bridged_root(scale, family);
            const Vector profile = square(root, fusion);
            if (boundary_failure(profile) >= 0) {
                continue;
            }
            ++boundary_admissible;
            const auto [left, right] = current_failure(profile, fusion);
            if (left >= 0) {
                std::cout << "SU2_INTERVAL_BRIDGE result=FULL_CURRENT_COUNTEREXAMPLE"
                          << " scale=" << scale
                          << " family=" << family
                          << " left=" << left
                          << " right=" << right
                          << " root=";
                print_vector(root);
                std::cout << " square=";
                print_vector(profile);
                std::cout << '\n';
                return EXIT_FAILURE;
            }
            const Vector image = apply_matrix(fusion[1U], root);
            const Vector image_square = square(image, fusion);
            const int failure = boundary_failure(image_square);
            if (failure >= 0) {
                std::cout << "SU2_INTERVAL_BRIDGE result=INSERTION_COUNTEREXAMPLE"
                          << " scale=" << scale
                          << " family=" << family
                          << " boundary=" << failure
                          << " root=";
                print_vector(root);
                std::cout << " square=";
                print_vector(profile);
                std::cout << " image=";
                print_vector(image);
                std::cout << " image_square=";
                print_vector(image_square);
                std::cout << " image_margin="
                          << boundary_margin(image_square, failure);
                std::cout << '\n';
                return EXIT_FAILURE;
            }
        }
        std::cout << "SU2_INTERVAL_BRIDGE"
                  << " maximum_scale=" << maximum_scale
                  << " family=" << family
                  << " boundary_admissible=" << boundary_admissible
                  << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
