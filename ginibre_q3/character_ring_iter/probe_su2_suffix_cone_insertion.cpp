#include <boost/multiprecision/cpp_int.hpp>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Pair = std::pair<int, int>;
using Vector = std::vector<cpp_int>;
using Matrix = std::vector<Vector>;

struct ExtremeRay {
    int cutoff = -1;
    int gap = -1;
    std::size_t positive = 0U;
    bool has_negative = false;
    std::size_t negative = 0U;
};

int positive_argument(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0LL
        || parsed > static_cast<long long>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

std::vector<Pair> wedge_basis(int maximum_label) {
    std::vector<Pair> basis;
    for (int left = 0; left <= maximum_label; ++left) {
        for (int right = left + 1; right <= maximum_label; ++right) {
            basis.emplace_back(left, right);
        }
    }
    return basis;
}

std::vector<std::vector<std::size_t>> wedge_indices(
    int maximum_label,
    const std::vector<Pair>& basis) {
    const std::size_t size = static_cast<std::size_t>(maximum_label + 1);
    std::vector<std::vector<std::size_t>> indices(
        size,
        std::vector<std::size_t>(size, basis.size()));
    for (std::size_t index = 0U; index < basis.size(); ++index) {
        const auto [left, right] = basis[index];
        indices[static_cast<std::size_t>(left)]
               [static_cast<std::size_t>(right)] = index;
    }
    return indices;
}

std::vector<std::vector<int>> ordinary_fusion(
    int input_maximum,
    int factor) {
    const int output_maximum = input_maximum + factor;
    std::vector<std::vector<int>> matrix(
        static_cast<std::size_t>(output_maximum + 1),
        std::vector<int>(static_cast<std::size_t>(input_maximum + 1), 0));
    for (int source = 0; source <= input_maximum; ++source) {
        const int lower = std::abs(source - factor);
        const int upper = source + factor;
        for (int target = lower; target <= upper; ++target) {
            matrix[static_cast<std::size_t>(target)]
                  [static_cast<std::size_t>(source)] = 1;
        }
    }
    return matrix;
}

Matrix second_compound(
    const std::vector<std::vector<int>>& matrix,
    const std::vector<Pair>& output_basis,
    const std::vector<Pair>& input_basis) {
    Matrix result(output_basis.size(), Vector(input_basis.size(), 0));
    for (std::size_t row = 0U; row < output_basis.size(); ++row) {
        const auto [row_left, row_right] = output_basis[row];
        for (std::size_t column = 0U; column < input_basis.size(); ++column) {
            const auto [column_left, column_right] = input_basis[column];
            result[row][column]
                = matrix[static_cast<std::size_t>(row_left)]
                        [static_cast<std::size_t>(column_left)]
                    * matrix[static_cast<std::size_t>(row_right)]
                            [static_cast<std::size_t>(column_right)]
                - matrix[static_cast<std::size_t>(row_left)]
                        [static_cast<std::size_t>(column_right)]
                    * matrix[static_cast<std::size_t>(row_right)]
                            [static_cast<std::size_t>(column_left)];
        }
    }
    return result;
}

std::vector<ExtremeRay> suffix_extreme_rays(
    int maximum_label,
    const std::vector<std::vector<std::size_t>>& indices) {
    std::vector<ExtremeRay> rays;
    for (int gap = 1; gap <= maximum_label; ++gap) {
        for (int cutoff = 0; cutoff + gap <= maximum_label; ++cutoff) {
            ExtremeRay ray;
            ray.cutoff = cutoff;
            ray.gap = gap;
            ray.positive
                = indices[static_cast<std::size_t>(cutoff)]
                         [static_cast<std::size_t>(cutoff + gap)];
            if (cutoff > 0) {
                ray.has_negative = true;
                ray.negative
                    = indices[static_cast<std::size_t>(cutoff - 1)]
                             [static_cast<std::size_t>(cutoff + gap - 1)];
            }
            rays.push_back(ray);
        }
    }
    return rays;
}

Vector ray_vector(std::size_t size, const ExtremeRay& ray) {
    Vector result(size, 0);
    result[ray.positive] = 1;
    if (ray.has_negative) {
        result[ray.negative] = -1;
    }
    return result;
}

cpp_int minimum_suffix(
    const Vector& vector,
    int maximum_label,
    const std::vector<std::vector<std::size_t>>& indices) {
    cpp_int minimum = 0;
    for (int gap = 1; gap <= maximum_label; ++gap) {
        cpp_int suffix = 0;
        for (int cutoff = maximum_label - gap; cutoff >= 0; --cutoff) {
            suffix
                += vector[
                    indices[static_cast<std::size_t>(cutoff)]
                           [static_cast<std::size_t>(cutoff + gap)]];
            if (suffix < minimum) {
                minimum = suffix;
            }
        }
    }
    return minimum;
}

Vector apply_compound_to_ray(
    const Matrix& compound,
    const ExtremeRay& ray) {
    Vector result(compound.size(), 0);
    for (std::size_t row = 0U; row < compound.size(); ++row) {
        result[row] = compound[row][ray.positive];
        if (ray.has_negative) {
            result[row] -= compound[row][ray.negative];
        }
    }
    return result;
}

cpp_int dot(const Vector& left, const Vector& right) {
    cpp_int result = 0;
    for (std::size_t index = 0U; index < left.size(); ++index) {
        result += left[index] * right[index];
    }
    return result;
}

std::string ray_description(
    const ExtremeRay& ray,
    const std::vector<Pair>& basis) {
    const auto [positive_left, positive_right] = basis[ray.positive];
    std::string result
        = "+e_(" + std::to_string(positive_left) + ','
        + std::to_string(positive_right) + ')';
    if (ray.has_negative) {
        const auto [negative_left, negative_right] = basis[ray.negative];
        result
            += "-e_(" + std::to_string(negative_left) + ','
            + std::to_string(negative_right) + ')';
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_input_label = argc >= 2
            ? positive_argument(argv[1], "maximum_input_label")
            : 20;
        const int maximum_factor = argc >= 3
            ? positive_argument(argv[2], "maximum_factor")
            : 10;
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: probe_su2_suffix_cone_insertion "
                "[maximum_input_label] [maximum_factor]");
        }

        std::uint64_t cases = 0U;
        std::uint64_t bilinear_entries = 0U;
        for (int input_maximum = 1;
             input_maximum <= maximum_input_label;
             ++input_maximum) {
            const std::vector<Pair> input_basis
                = wedge_basis(input_maximum);
            const auto input_indices
                = wedge_indices(input_maximum, input_basis);
            const std::vector<ExtremeRay> rays
                = suffix_extreme_rays(input_maximum, input_indices);
            for (int factor = 1; factor <= maximum_factor; ++factor) {
                ++cases;
                const int output_maximum = input_maximum + factor;
                const std::vector<Pair> output_basis
                    = wedge_basis(output_maximum);
                const auto fusion
                    = ordinary_fusion(input_maximum, factor);
                const Matrix compound
                    = second_compound(fusion, output_basis, input_basis);

                std::vector<Vector> images;
                images.reserve(rays.size());
                for (const ExtremeRay& ray : rays) {
                    images.push_back(apply_compound_to_ray(compound, ray));
                }

                for (std::size_t left = 0U; left < rays.size(); ++left) {
                    for (std::size_t right = left;
                         right < rays.size();
                         ++right) {
                        ++bilinear_entries;
                        const cpp_int value = dot(images[left], images[right]);
                        if (value >= 0) {
                            continue;
                        }
                        const Vector left_ray
                            = ray_vector(input_basis.size(), rays[left]);
                        const Vector right_ray
                            = ray_vector(input_basis.size(), rays[right]);
                        const cpp_int left_minimum = minimum_suffix(
                            left_ray, input_maximum, input_indices);
                        const cpp_int right_minimum = minimum_suffix(
                            right_ray, input_maximum, input_indices);
                        if (left_minimum < 0 || right_minimum < 0) {
                            throw std::logic_error(
                                "constructed ray is outside suffix cone");
                        }
                        std::cout
                            << "SU2_SUFFIX_CONE_INSERTION"
                            << " input_maximum=" << input_maximum
                            << " factor=" << factor
                            << " cases=" << cases
                            << " bilinear_entries=" << bilinear_entries
                            << " left_cutoff=" << rays[left].cutoff
                            << " left_gap=" << rays[left].gap
                            << " left_ray="
                            << ray_description(rays[left], input_basis)
                            << " left_minimum_suffix=" << left_minimum
                            << " right_cutoff=" << rays[right].cutoff
                            << " right_gap=" << rays[right].gap
                            << " right_ray="
                            << ray_description(rays[right], input_basis)
                            << " right_minimum_suffix=" << right_minimum
                            << " bilinear_value=" << value
                            << " result=COUNTEREXAMPLE_EXACT"
                            << '\n';
                        return EXIT_SUCCESS;
                    }
                }
            }
        }

        std::cout
            << "SU2_SUFFIX_CONE_INSERTION"
            << " maximum_input_label=" << maximum_input_label
            << " maximum_factor=" << maximum_factor
            << " cases=" << cases
            << " bilinear_entries=" << bilinear_entries
            << " result=NO_COUNTEREXAMPLE_BOUNDED"
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
