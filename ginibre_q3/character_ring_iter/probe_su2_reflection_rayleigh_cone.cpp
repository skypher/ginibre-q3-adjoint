#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdint>
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
    Vector result(matrix.size(), 0);
    for (std::size_t row = 0U; row < matrix.size(); ++row) {
        for (std::size_t column = 0U; column < input.size(); ++column) {
            result[row] += matrix[row][column] * input[column];
        }
    }
    return result;
}

Integer inner(const Vector& left, const Vector& right) {
    Integer result = 0;
    for (std::size_t index = 0U; index < left.size(); ++index) {
        result += left[index] * right[index];
    }
    return result;
}

Vector reflect(const Vector& input) {
    Vector result(input.size(), 0);
    for (std::size_t index = 0U; index < input.size(); ++index) {
        result[index] = input[input.size() - 1U - index];
    }
    return result;
}

int first_failed_boundary(
    const Vector& vector,
    const std::vector<Matrix>& fusion
) {
    const Integer norm = inner(vector, vector);
    const Integer overlap = inner(vector, reflect(vector));
    for (std::size_t label = 0U; label < fusion.size(); ++label) {
        const Vector image = apply_matrix(fusion[label], vector);
        const Integer direct = inner(vector, image);
        const Integer reflected = inner(vector, reflect(image));
        if (norm * reflected < overlap * direct) {
            return static_cast<int>(label);
        }
    }
    return -1;
}

Vector boundary_numerators(
    const Vector& vector,
    const std::vector<Matrix>& fusion
) {
    const Integer norm = inner(vector, vector);
    const Integer overlap = inner(vector, reflect(vector));
    Vector result(fusion.size(), 0);
    for (std::size_t label = 0U; label < fusion.size(); ++label) {
        const Vector image = apply_matrix(fusion[label], vector);
        result[label] = norm * inner(vector, reflect(image))
            - overlap * inner(vector, image);
    }
    return result;
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

bool increment(Vector& values, int maximum_coordinate) {
    for (std::size_t reverse = 0U; reverse < values.size(); ++reverse) {
        const std::size_t index = values.size() - 1U - reverse;
        if (values[index] < maximum_coordinate) {
            ++values[index];
            std::fill(values.begin() + static_cast<std::ptrdiff_t>(index + 1U),
                      values.end(), 0);
            return true;
        }
    }
    return false;
}

bool supported_on_one_label_parity(const Vector& values) {
    int parity = -1;
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
    }
    return parity >= 0;
}

bool parity_interval(const Vector& values) {
    if (!supported_on_one_label_parity(values)) {
        return false;
    }
    std::size_t first = values.size();
    std::size_t last = 0U;
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (values[index] != 0) {
            first = std::min(first, index);
            last = std::max(last, index);
        }
    }
    if (first == values.size()) {
        return false;
    }
    for (std::size_t index = first; index <= last; index += 2U) {
        if (values[index] == 0) {
            return false;
        }
    }
    return true;
}

bool unimodal_parity_interval(const Vector& values) {
    if (!parity_interval(values)) {
        return false;
    }
    std::size_t first = values.size();
    std::size_t last = 0U;
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (values[index] != 0) {
            first = std::min(first, index);
            last = std::max(last, index);
        }
    }
    bool decreasing = false;
    for (std::size_t index = first + 2U; index <= last; index += 2U) {
        const Integer& previous = values[index - 2U];
        const Integer& current = values[index];
        if (current < previous) {
            decreasing = true;
        } else if (decreasing && current > previous) {
            return false;
        }
    }
    return true;
}

bool log_concave_parity_interval(const Vector& values) {
    if (!parity_interval(values)) {
        return false;
    }
    std::size_t first = values.size();
    std::size_t last = 0U;
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (values[index] != 0) {
            first = std::min(first, index);
            last = std::max(last, index);
        }
    }
    for (std::size_t index = first + 2U; index + 2U <= last; index += 2U) {
        if (values[index] * values[index]
            < values[index - 2U] * values[index + 2U]) {
            return false;
        }
    }
    return true;
}

Vector square_profile(
    const Vector& vector,
    const std::vector<Matrix>& fusion
) {
    Vector result(vector.size(), 0);
    for (std::size_t label = 0U; label < vector.size(); ++label) {
        if (vector[label] == 0) {
            continue;
        }
        const Vector product = apply_matrix(fusion[label], vector);
        for (std::size_t index = 0U; index < result.size(); ++index) {
            result[index] += vector[label] * product[index];
        }
    }
    return result;
}

std::pair<int, int> first_failed_anchored_current(
    const Vector& root,
    const std::vector<Matrix>& fusion
) {
    const Vector square = square_profile(root, fusion);
    const Integer d_zero = square[0U];
    for (std::size_t left = 0U; left < fusion.size(); ++left) {
        const Vector row = apply_matrix(fusion[left], square);
        for (std::size_t right = 0U; right < fusion.size(); ++right) {
            if (d_zero * row[right] < square[left] * square[right]) {
                return {
                    static_cast<int>(left), static_cast<int>(right)
                };
            }
        }
    }
    return {-1, -1};
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 4) {
            throw std::runtime_error(
                "usage: probe_su2_reflection_rayleigh_cone "
                "MAXIMUM_LEVEL MAXIMUM_COORDINATE "
                "[all|parity|parity-interval|parity-unimodal|parity-lc|"
                "parity-lc-square-lc|boundary-to-full|"
                "parity-interval-boundary-to-full|"
                "parity-unimodal-boundary-to-full]"
            );
        }
        const int maximum_level = parse_positive(argv[1], "maximum level");
        const int maximum_coordinate = parse_positive(
            argv[2], "maximum coordinate"
        );
        const std::string mode = argc == 4 ? argv[3] : "all";
        if (mode != "all" && mode != "parity"
            && mode != "parity-interval" && mode != "parity-lc"
            && mode != "parity-unimodal"
            && mode != "parity-lc-square-lc"
            && mode != "boundary-to-full"
            && mode != "parity-interval-boundary-to-full"
            && mode != "parity-unimodal-boundary-to-full") {
            throw std::runtime_error("invalid cone mode");
        }
        std::uint64_t tested_vectors = 0U;
        std::uint64_t cone_vectors = 0U;
        std::uint64_t tested_factor_images = 0U;
        for (int level = 1; level <= maximum_level; ++level) {
            std::vector<Matrix> fusion;
            fusion.reserve(static_cast<std::size_t>(level + 1));
            for (int label = 0; label <= level; ++label) {
                fusion.push_back(fusion_matrix(level, label));
            }
            Vector vector(static_cast<std::size_t>(level + 1), 0);
            bool more = true;
            while (more) {
                bool nonzero = false;
                for (const Integer& entry : vector) {
                    nonzero = nonzero || entry != 0;
                }
                const bool root_log_concave
                    = log_concave_parity_interval(vector);
                const bool square_log_concave = root_log_concave
                    && log_concave_parity_interval(square_profile(
                        vector, fusion
                    ));
                if (nonzero
                    && (mode == "all"
                        || (mode == "parity"
                            && supported_on_one_label_parity(vector))
                        || (mode == "parity-interval"
                            && parity_interval(vector))
                        || (mode == "parity-unimodal"
                            && unimodal_parity_interval(vector))
                        || (mode == "parity-lc" && root_log_concave)
                        || (mode == "parity-lc-square-lc"
                            && square_log_concave)
                        || (mode == "boundary-to-full"
                            && square_log_concave)
                        || (mode == "parity-interval-boundary-to-full"
                            && parity_interval(vector))
                        || (mode == "parity-unimodal-boundary-to-full"
                            && unimodal_parity_interval(vector)))) {
                    ++tested_vectors;
                    if (first_failed_boundary(vector, fusion) < 0) {
                        ++cone_vectors;
                        if (mode == "boundary-to-full"
                            || mode == "parity-interval-boundary-to-full"
                            || mode == "parity-unimodal-boundary-to-full") {
                            const auto [left, right]
                                = first_failed_anchored_current(
                                    vector, fusion
                                );
                            if (left >= 0) {
                                std::cout
                                    << "SU2_REFLECTION_RAYLEIGH_CONE"
                                    << " result=FULL_CURRENT_COUNTEREXAMPLE"
                                    << " level=" << level
                                    << " left=" << left
                                    << " right=" << right
                                    << " source=";
                                print_vector(vector);
                                std::cout << " square=";
                                print_vector(square_profile(vector, fusion));
                                std::cout << '\n';
                                return EXIT_FAILURE;
                            }
                        }
                        if (mode != "boundary-to-full"
                            && mode != "parity-interval-boundary-to-full"
                            && mode != "parity-unimodal-boundary-to-full") {
                            for (int factor = 0; factor <= level; ++factor) {
                                ++tested_factor_images;
                                const Vector image = apply_matrix(
                                    fusion[static_cast<std::size_t>(factor)],
                                    vector
                                );
                                const int failed_boundary
                                    = first_failed_boundary(image, fusion);
                                if (failed_boundary >= 0) {
                                    std::cout
                                        << "SU2_REFLECTION_RAYLEIGH_CONE"
                                        << " result=COUNTEREXAMPLE"
                                        << " level=" << level
                                        << " factor=" << factor
                                        << " boundary=" << failed_boundary
                                        << " source=";
                                    print_vector(vector);
                                    std::cout << " image=";
                                    print_vector(image);
                                    std::cout << " source_margins=";
                                    print_vector(boundary_numerators(
                                        vector, fusion
                                    ));
                                    std::cout << " image_margins=";
                                    print_vector(boundary_numerators(
                                        image, fusion
                                    ));
                                    std::cout << '\n';
                                    return EXIT_FAILURE;
                                }
                            }
                        }
                    }
                }
                more = increment(vector, maximum_coordinate);
            }
        }
        std::cout
            << "SU2_REFLECTION_RAYLEIGH_CONE"
            << " maximum_level=" << maximum_level
            << " maximum_coordinate=" << maximum_coordinate
            << " mode=" << mode
            << " tested_vectors=" << tested_vectors
            << " cone_vectors=" << cone_vectors
            << " tested_factor_images=" << tested_factor_images
            << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
