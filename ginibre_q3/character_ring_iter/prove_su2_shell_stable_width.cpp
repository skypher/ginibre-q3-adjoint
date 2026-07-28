#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Pair = std::pair<int, int>;

struct DatasetKey {
    std::string family;
    int first_order = 0;
    int second_order = 0;

    bool operator<(const DatasetKey& other) const {
        if (family != other.family) {
            return family < other.family;
        }
        if (first_order != other.first_order) {
            return first_order < other.first_order;
        }
        return second_order < other.second_order;
    }
};

Integer binomial_integer(int top, int bottom) {
    if (bottom < 0 || top < bottom) {
        return 0;
    }
    Integer result = 1;
    for (int index = 1; index <= bottom; ++index) {
        result *= top - bottom + index;
        result /= index;
    }
    return result;
}

int parse_field(const std::string& token, const std::string& prefix) {
    if (token.rfind(prefix, 0U) != 0U) {
        throw std::runtime_error("missing field " + prefix);
    }
    return std::stoi(token.substr(prefix.size()));
}

std::string parse_string_field(
    const std::string& token,
    const std::string& prefix
) {
    if (token.rfind(prefix, 0U) != 0U) {
        throw std::runtime_error("missing field " + prefix);
    }
    return token.substr(prefix.size());
}

void insert_value(
    std::map<DatasetKey, std::map<Pair, Integer>>& cones,
    std::map<DatasetKey, std::map<int, Integer>>& rays,
    const DatasetKey& key,
    bool cone,
    int first,
    int second,
    const Integer& value
) {
    if (cone) {
        const auto inserted = cones[key].emplace(
            Pair{first, second},
            value
        );
        if (!inserted.second && inserted.first->second != value) {
            throw std::runtime_error("inconsistent cone value");
        }
    } else {
        const auto inserted = rays[key].emplace(first, value);
        if (!inserted.second && inserted.first->second != value) {
            throw std::runtime_error("inconsistent ray value");
        }
    }
}

}  // namespace

int main() {
    try {
        std::map<DatasetKey, std::map<Pair, Integer>> cones;
        std::map<DatasetKey, std::map<int, Integer>> rays;
        std::uint64_t input_coefficients = 0U;
        std::string line;
        while (std::getline(std::cin, line)) {
            if (
                line.rfind(
                    "FIXED_OFFSET_CONE_COEFFICIENT ",
                    0U
                ) != 0U
            ) {
                continue;
            }
            std::istringstream stream(line);
            std::string marker;
            std::string offset_token;
            std::string key_token;
            std::string first_token;
            std::string second_token;
            std::string value_token;
            stream
                >> marker
                >> offset_token
                >> key_token
                >> first_token
                >> second_token
                >> value_token;
            if (!stream || marker != "FIXED_OFFSET_CONE_COEFFICIENT") {
                throw std::runtime_error("malformed coefficient record");
            }
            const int offset = parse_field(offset_token, "offset=");
            const int d = offset - 1;
            const std::string key =
                parse_string_field(key_token, "key=");
            const int first_order =
                parse_field(first_token, "first_order=");
            const int second_order =
                parse_field(second_token, "second_order=");
            const Integer value(
                parse_string_field(value_token, "value=")
            );
            if (d < 2 || d > 26) {
                throw std::runtime_error("unexpected width");
            }
            const DatasetKey base{
                "",
                first_order,
                second_order
            };

            if (key.rfind("E_a", 0U) == 0U) {
                const int a = std::stoi(key.substr(3U));
                const int alpha = a + 2 * d;
                const int beta = 2 * d - a;
                if (
                    alpha < 0
                    || beta < 0
                    || alpha + beta != 4 * d
                ) {
                    throw std::runtime_error("invalid even profile");
                }
                const int residue = alpha % 4;
                const int beta_residue = (4 - residue) % 4;
                const int p = (alpha - residue) / 4;
                const int q = (beta - beta_residue) / 4;
                DatasetKey dataset = base;
                bool cone = false;
                int first = 0;
                int second = 0;
                if (residue == 0) {
                    if (p >= 2) {
                        dataset.family = "E_r0_cone";
                        cone = true;
                        first = p - 2;
                        second = q;
                    } else if (p == 1) {
                        dataset.family = "E_r0_p1_ray";
                        first = q - 1;
                    } else {
                        dataset.family = "E_r0_p0_ray";
                        first = q - 2;
                    }
                } else if (p >= 1) {
                    dataset.family =
                        "E_r" + std::to_string(residue) + "_cone";
                    cone = true;
                    first = p - 1;
                    second = q;
                } else {
                    dataset.family =
                        "E_r" + std::to_string(residue) + "_p0_ray";
                    first = q - 1;
                }
                if (first < 0 || second < 0) {
                    throw std::runtime_error(
                        "invalid even residue chart"
                    );
                }
                insert_value(
                    cones,
                    rays,
                    dataset,
                    cone,
                    first,
                    second,
                    value
                );
            } else if (key.rfind("O_b", 0U) == 0U) {
                DatasetKey dataset = base;
                dataset.family = "O_full_ray";
                insert_value(
                    cones,
                    rays,
                    dataset,
                    false,
                    d - 2,
                    0,
                    value
                );
            } else if (key.rfind("T_b", 0U) == 0U) {
                const std::size_t separator = key.find("_a", 3U);
                if (separator == std::string::npos) {
                    throw std::runtime_error("invalid total key");
                }
                const int b = std::stoi(
                    key.substr(3U, separator - 3U)
                );
                const std::string type = key.substr(separator + 2U);
                if (type != "full" && type != "empty") {
                    throw std::runtime_error(
                        "unexpected total cone type"
                    );
                }
                const int lower =
                    d % 2 == 0
                        ? -3 * d / 2
                        : (1 - 3 * d) / 2;
                const int normalized = b - lower - 1;
                const int complement = 3 * d - 1 - normalized;
                if (
                    normalized < 0
                    || complement < 0
                    || normalized + complement != 3 * d - 1
                ) {
                    throw std::runtime_error(
                        "invalid odd profile"
                    );
                }
                const int residue = normalized % 3;
                const int complement_residue = (2 - residue + 3) % 3;
                const int p = (normalized - residue) / 3;
                const int q =
                    (complement - complement_residue) / 3;
                DatasetKey dataset = base;
                dataset.family =
                    "T_" + type + "_r"
                    + std::to_string(residue);
                bool cone = false;
                int first = 0;
                int second = 0;
                if (p >= 1) {
                    dataset.family += "_cone";
                    cone = true;
                    first = p - 1;
                    second = q;
                } else {
                    dataset.family += "_p0_ray";
                    first = q - 1;
                }
                if (first < 0 || second < 0) {
                    throw std::runtime_error(
                        "invalid total residue chart"
                    );
                }
                insert_value(
                    cones,
                    rays,
                    dataset,
                    cone,
                    first,
                    second,
                    value
                );
            } else {
                throw std::runtime_error("unknown cone family");
            }
            ++input_coefficients;
        }

        std::uint64_t finite_width_coefficients = 0U;
        constexpr int boundary_depth = 3;
        std::map<DatasetKey, std::map<Pair, Integer>> parity_cones;
        std::map<DatasetKey, std::map<int, Integer>>
            cone_boundary_rays;
        for (const auto& [key, values] : cones) {
            for (const auto& [index, value] : values) {
                DatasetKey split = key;
                split.family +=
                    "_parity" + std::to_string(index.first % 2)
                    + std::to_string(index.second % 2);
                const int first = index.first / 2;
                const int second = index.second / 2;
                if (
                    first < boundary_depth
                    && second < boundary_depth
                ) {
                    if (value < 0) {
                        throw std::runtime_error(
                            "negative finite-width cone coefficient"
                        );
                    }
                    ++finite_width_coefficients;
                } else if (first < boundary_depth) {
                    split.family +=
                        "_second_axis_first"
                        + std::to_string(first);
                    cone_boundary_rays[split][
                        second - boundary_depth
                    ] = value;
                } else if (second < boundary_depth) {
                    split.family +=
                        "_first_axis_second"
                        + std::to_string(second);
                    cone_boundary_rays[split][
                        first - boundary_depth
                    ] = value;
                } else {
                    split.family += "_interior";
                    parity_cones[split][{
                        first - boundary_depth,
                        second - boundary_depth
                    }] = value;
                }
            }
        }
        cones = std::move(parity_cones);
        std::map<DatasetKey, std::map<int, Integer>> parity_rays =
            std::move(cone_boundary_rays);
        for (const auto& [key, values] : rays) {
            for (const auto& [index, value] : values) {
                DatasetKey split = key;
                split.family +=
                    "_parity" + std::to_string(index % 2);
                const int coordinate = index / 2;
                if (coordinate < boundary_depth) {
                    if (value < 0) {
                        throw std::runtime_error(
                            "negative finite-width coefficient"
                        );
                    }
                    ++finite_width_coefficients;
                } else {
                    parity_rays[split][
                        coordinate - boundary_depth
                    ] = value;
                }
            }
        }
        rays = std::move(parity_rays);

        std::uint64_t proved_coefficients = 0U;
        Integer minimum_coefficient = 0;
        bool initialized_minimum = false;
        for (const auto& [key, values] : rays) {
            const int degree = std::max(
                0,
                5 - key.first_order - key.second_order
            );
            std::vector<Integer> coefficients;
            bool invalid_ray = false;
            for (int order = 0; order <= degree + 1; ++order) {
                Integer coefficient = 0;
                for (int index = 0; index <= order; ++index) {
                    const auto found = values.find(index);
                    if (found == values.end()) {
                        throw std::runtime_error(
                            "incomplete uniform ray: " + key.family
                        );
                    }
                    const Integer term =
                        binomial_integer(order, index)
                        * found->second;
                    if ((order - index) % 2 == 0) {
                        coefficient += term;
                    } else {
                        coefficient -= term;
                    }
                }
                if (
                    (order <= degree && coefficient < 0)
                    || (order == degree + 1 && coefficient != 0)
                ) {
                    std::cerr
                        << "FAILED_UNIFORM_STABLE_RAY"
                        << " family=" << key.family
                        << " crossing_order=" << key.first_order
                        << ',' << key.second_order
                        << " order=" << order
                        << " residual_degree=" << degree
                        << " coefficient=" << coefficient;
                    for (int sample = 0; sample <= degree + 1; ++sample) {
                        const auto found = values.find(sample);
                        if (found != values.end()) {
                            std::cerr
                                << " v" << sample << '='
                                << found->second;
                        }
                    }
                    std::cerr << '\n';
                    invalid_ray = true;
                }
                if (order <= degree) {
                    coefficients.push_back(coefficient);
                    ++proved_coefficients;
                    if (
                        !initialized_minimum
                        || coefficient < minimum_coefficient
                    ) {
                        minimum_coefficient = coefficient;
                        initialized_minimum = true;
                    }
                }
            }
            if (invalid_ray) {
                throw std::runtime_error(
                    "uniform stable ray certificate failed"
                );
            }
            for (const auto& [index, value] : values) {
                Integer reconstructed = 0;
                for (
                    std::size_t order = 0U;
                    order < coefficients.size();
                    ++order
                ) {
                    reconstructed +=
                        coefficients[order]
                        * binomial_integer(
                            index,
                            static_cast<int>(order)
                        );
                }
                if (reconstructed != value) {
                    throw std::runtime_error(
                        "uniform ray reconstruction mismatch"
                    );
                }
            }
        }

        for (const auto& [key, values] : cones) {
            const int degree = std::max(
                0,
                5 - key.first_order - key.second_order
            );
            std::map<Pair, Integer> coefficients;
            for (int first_order = 0;
                 first_order <= degree + 1;
                 ++first_order) {
                for (int second_order = 0;
                     second_order <= degree + 1 - first_order;
                     ++second_order) {
                    Integer coefficient = 0;
                    for (int first = 0;
                         first <= first_order;
                         ++first) {
                        for (int second = 0;
                             second <= second_order;
                             ++second) {
                            const auto found = values.find({
                                first,
                                second
                            });
                            if (found == values.end()) {
                                throw std::runtime_error(
                                    "incomplete uniform cone: "
                                    + key.family
                                );
                            }
                            Integer term =
                                binomial_integer(first_order, first)
                                * binomial_integer(second_order, second)
                                * found->second;
                            if (
                                (
                                    first_order - first
                                    + second_order - second
                                ) % 2 == 0
                            ) {
                                coefficient += term;
                            } else {
                                coefficient -= term;
                            }
                        }
                    }
                    const int order = first_order + second_order;
                    if (
                        (order <= degree && coefficient < 0)
                        || (
                            order == degree + 1
                            && coefficient != 0
                        )
                    ) {
                        std::cerr
                            << "FAILED_UNIFORM_STABLE_CONE"
                            << " family=" << key.family
                            << " crossing_order=" << key.first_order
                            << ',' << key.second_order
                            << " order=" << first_order
                            << ',' << second_order
                            << " residual_degree=" << degree
                            << " coefficient=" << coefficient << '\n';
                        throw std::runtime_error(
                            "uniform stable cone certificate failed"
                        );
                    }
                    if (order <= degree) {
                        coefficients[{
                            first_order,
                            second_order
                        }] = coefficient;
                        ++proved_coefficients;
                        if (
                            !initialized_minimum
                            || coefficient < minimum_coefficient
                        ) {
                            minimum_coefficient = coefficient;
                            initialized_minimum = true;
                        }
                    }
                }
            }
            for (const auto& [index, value] : values) {
                Integer reconstructed = 0;
                for (const auto& [order, coefficient] : coefficients) {
                    reconstructed +=
                        coefficient
                        * binomial_integer(index.first, order.first)
                        * binomial_integer(index.second, order.second);
                }
                if (reconstructed != value) {
                    throw std::runtime_error(
                        "uniform cone reconstruction mismatch"
                    );
                }
            }
        }

        std::cout
            << "SU2_SHELL_STABLE_WIDTH_CERTIFICATE"
            << " input_coefficients=" << input_coefficients
            << " ray_families=" << rays.size()
            << " cone_families=" << cones.size()
            << " finite_width_coefficients="
            << finite_width_coefficients
            << " proved_newton_coefficients=" << proved_coefficients
            << " minimum_coefficient=" << minimum_coefficient
            << " maximum_total_degree=10"
            << " result=PASS_UNIFORM_NONNEGATIVE_NEWTON_EXPANSIONS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
