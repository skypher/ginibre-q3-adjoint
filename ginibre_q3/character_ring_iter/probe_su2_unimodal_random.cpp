#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
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

std::uint64_t bounded_nonzero(std::uint64_t& state, std::uint64_t upper) {
    return UINT64_C(1) + splitmix64(state) % upper;
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

Vector fusion_square(const Vector& root, const std::vector<Matrix>& fusion) {
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
            < profile[static_cast<std::size_t>(radius)] * profile.back()) {
            return radius;
        }
    }
    return -1;
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

Vector random_unimodal_profile(
    std::uint64_t& state,
    int level,
    std::uint64_t coefficient_cap
) {
    int parity = static_cast<int>(splitmix64(state) & UINT64_C(1));
    if (parity > level) {
        parity = 0;
    }
    const int count = (level - parity) / 2 + 1;
    const int left = static_cast<int>(
        splitmix64(state) % static_cast<std::uint64_t>(count)
    );
    const int right = left + static_cast<int>(
        splitmix64(state)
        % static_cast<std::uint64_t>(count - left)
    );
    const int peak = left + static_cast<int>(
        splitmix64(state)
        % static_cast<std::uint64_t>(right - left + 1)
    );
    std::vector<std::uint64_t> compressed(
        static_cast<std::size_t>(count), 0U
    );
    compressed[static_cast<std::size_t>(peak)]
        = bounded_nonzero(state, coefficient_cap);
    for (int index = peak - 1; index >= left; --index) {
        compressed[static_cast<std::size_t>(index)] = bounded_nonzero(
            state, compressed[static_cast<std::size_t>(index + 1)]
        );
    }
    for (int index = peak + 1; index <= right; ++index) {
        compressed[static_cast<std::size_t>(index)] = bounded_nonzero(
            state, compressed[static_cast<std::size_t>(index - 1)]
        );
    }
    Vector output(static_cast<std::size_t>(level + 1), 0);
    for (int index = left; index <= right; ++index) {
        output[static_cast<std::size_t>(parity + 2 * index)]
            = compressed[static_cast<std::size_t>(index)];
    }
    return output;
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
        if (argc != 4 && argc != 5) {
            throw std::runtime_error(
                "usage: probe_su2_unimodal_random SAMPLES MAXIMUM_LEVEL "
                "MAXIMUM_COEFFICIENT [all-factors]"
            );
        }
        const std::uint64_t samples = parse_positive_u64(argv[1], "samples");
        const int maximum_level = parse_positive_int(argv[2], "maximum level");
        const std::uint64_t coefficient_cap = parse_positive_u64(
            argv[3], "maximum coefficient"
        );
        const bool all_factors = argc == 5;
        if (all_factors && std::string(argv[4]) != "all-factors") {
            throw std::runtime_error("invalid factor mode");
        }
        std::vector<std::vector<Matrix>> fusion_by_level(
            static_cast<std::size_t>(maximum_level + 1)
        );
        for (int level = 1; level <= maximum_level; ++level) {
            auto& fusion = fusion_by_level[static_cast<std::size_t>(level)];
            fusion.reserve(static_cast<std::size_t>(level + 1));
            for (int factor = 0; factor <= level; ++factor) {
                fusion.push_back(fusion_matrix(level, factor));
            }
        }
        std::uint64_t state = UINT64_C(0x4b1d2e3f5a697887);
        std::uint64_t boundary_admissible = 0U;
        std::uint64_t factor_images = 0U;
        for (std::uint64_t sample = 0U; sample < samples; ++sample) {
            const int level = 1 + static_cast<int>(
                splitmix64(state)
                % static_cast<std::uint64_t>(maximum_level)
            );
            const auto& fusion
                = fusion_by_level[static_cast<std::size_t>(level)];
            const Vector root = random_unimodal_profile(
                state, level, coefficient_cap
            );
            const Vector profile = fusion_square(root, fusion);
            if (boundary_failure(profile) >= 0) {
                continue;
            }
            ++boundary_admissible;
            const auto [left, right] = current_failure(profile, fusion);
            if (left >= 0) {
                std::cout << "SU2_UNIMODAL_RANDOM result=FULL_CURRENT_COUNTEREXAMPLE"
                          << " sample=" << sample
                          << " level=" << level
                          << " left=" << left
                          << " right=" << right
                          << " root=";
                print_vector(root);
                std::cout << " square=";
                print_vector(profile);
                std::cout << '\n';
                return EXIT_FAILURE;
            }
            const int first_factor = all_factors ? 1 : 1 + static_cast<int>(
                splitmix64(state) % static_cast<std::uint64_t>(level)
            );
            const int last_factor = all_factors ? level : first_factor;
            for (int factor = first_factor; factor <= last_factor; ++factor) {
                ++factor_images;
                const Vector image = apply_matrix(
                    fusion[static_cast<std::size_t>(factor)], root
                );
                const Vector image_square = fusion_square(image, fusion);
                const int boundary = boundary_failure(image_square);
                if (boundary >= 0) {
                    std::cout << "SU2_UNIMODAL_RANDOM result=INSERTION_COUNTEREXAMPLE"
                              << " sample=" << sample
                              << " level=" << level
                              << " factor=" << factor
                              << " boundary=" << boundary
                              << " root=";
                    print_vector(root);
                    std::cout << " image=";
                    print_vector(image);
                    std::cout << " square=";
                    print_vector(profile);
                    std::cout << " image_square=";
                    print_vector(image_square);
                    std::cout << '\n';
                    return EXIT_FAILURE;
                }
            }
        }
        std::cout << "SU2_UNIMODAL_RANDOM"
                  << " samples=" << samples
                  << " maximum_level=" << maximum_level
                  << " maximum_coefficient=" << coefficient_cap
                  << " all_factors=" << (all_factors ? "true" : "false")
                  << " boundary_admissible=" << boundary_admissible
                  << " factor_images=" << factor_images
                  << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
