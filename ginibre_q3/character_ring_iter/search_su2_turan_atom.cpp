#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

cpp_int coefficient(
    const std::vector<cpp_int>& polynomial,
    int offset,
    int charge
) {
    const int index = offset + charge;
    if (index < 0 || index >= static_cast<int>(polynomial.size())) {
        return 0;
    }
    return polynomial[static_cast<std::size_t>(index)];
}

cpp_int turan(
    const std::vector<cpp_int>& polynomial,
    int offset,
    int charge
) {
    const cpp_int center = coefficient(polynomial, offset, charge);
    return center * center
        - coefficient(polynomial, offset, charge - 2)
            * coefficient(polynomial, offset, charge + 2);
}

cpp_int central_difference_after_factor(
    const std::vector<cpp_int>& polynomial,
    int offset,
    int label,
    bool selected_minus
) {
    cpp_int result = 0;
    if (!selected_minus && (label & 1) == 0) {
        const cpp_int c0 = coefficient(polynomial, offset, 0);
        const cpp_int c2 = coefficient(polynomial, offset, 2);
        result += 2 * (c0 * c0 - c2 * c2);
    }
    for (int feature = label; feature > 0; feature -= 2) {
        const cpp_int a0 = coefficient(polynomial, offset, -feature)
            + (selected_minus ? -1 : 1)
                * coefficient(polynomial, offset, feature);
        const cpp_int a2 = coefficient(polynomial, offset, 2 - feature)
            + (selected_minus ? -1 : 1)
                * coefficient(polynomial, offset, 2 + feature);
        result += a0 * a0 - a2 * a2;
    }
    return result;
}

void print_case(
    int outer_label,
    const std::vector<int>& magnitudes,
    std::uint64_t minus_mask,
    const cpp_int& lower,
    const cpp_int& upper
) {
    std::cout << "outer_label=" << outer_label << " features=[";
    for (std::size_t index = 0; index < magnitudes.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << (((minus_mask >> index) & 1U) != 0U ? 'S' : 'C')
                  << magnitudes[index];
    }
    std::cout << "] H_p=" << lower << " H_(p+2)=" << upper << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 3 || argc > 4) {
            throw std::runtime_error(
                "usage: search_su2_turan_atom MAXIMUM_LABEL MAXIMUM_FACTORS "
                "[--even|--central]"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 1 || maximum_factors < 0
            || maximum_factors >= 63) {
            throw std::runtime_error("invalid search bound");
        }
        const bool even_only = argc == 4
            && std::string(argv[3]) == "--even";
        const bool central_mode = argc == 4
            && std::string(argv[3]) == "--central";
        if (argc == 4 && !even_only && !central_mode) {
            throw std::runtime_error("unknown mode");
        }

        std::uint64_t tested = 0;
        bool failed = false;
        std::vector<int> magnitudes;
        const auto visit = [&](const auto& self, int first) -> void {
            if (failed) {
                return;
            }
            if (!magnitudes.empty()) {
                const std::uint64_t masks
                    = std::uint64_t{1} << magnitudes.size();
                int total = 0;
                for (int magnitude : magnitudes) {
                    total += magnitude;
                }
                for (int outer = magnitudes.back();
                     outer <= maximum_label;
                     outer += even_only ? 2 : 1) {
                    if (((outer + total) & 1) != 0) {
                        continue;
                    }
                    for (std::uint64_t mask = 0; mask < masks; ++mask) {
                        std::vector<cpp_int> polynomial(
                            static_cast<std::size_t>(2 * total + 1)
                        );
                        polynomial[static_cast<std::size_t>(total)] = 1;
                        int support = 0;
                        for (std::size_t index = 0;
                             index < magnitudes.size();
                             ++index) {
                            const int magnitude = magnitudes[index];
                            const bool minus
                                = ((mask >> index) & 1U) != 0U;
                            std::vector<cpp_int> next(polynomial.size());
                            for (int charge = -support;
                                 charge <= support;
                                 ++charge) {
                                const cpp_int& value = polynomial[
                                    static_cast<std::size_t>(total + charge)
                                ];
                                if (value == 0) {
                                    continue;
                                }
                                next[static_cast<std::size_t>(
                                    total + charge + magnitude
                                )] += value;
                                next[static_cast<std::size_t>(
                                    total + charge - magnitude
                                )] += minus ? -value : value;
                            }
                            support += magnitude;
                            polynomial = std::move(next);
                        }
                        const cpp_int lower = central_mode
                            ? central_difference_after_factor(
                                polynomial, total, outer,
                                (std::popcount(mask) & 1) != 0
                            )
                            : turan(polynomial, total, outer);
                        const cpp_int upper = central_mode ? cpp_int{0}
                            : turan(polynomial, total, outer + 2);
                        ++tested;
                        if (lower < upper) {
                            print_case(outer, magnitudes, mask, lower, upper);
                            failed = true;
                            return;
                        }
                    }
                }
            }
            if (magnitudes.size()
                == static_cast<std::size_t>(maximum_factors)) {
                return;
            }
            const int first_magnitude = even_only
                ? (first + (first & 1)) : first;
            for (int magnitude = first_magnitude;
                 magnitude <= maximum_label;
                 magnitude += even_only ? 2 : 1) {
                magnitudes.push_back(magnitude);
                self(self, magnitude);
                magnitudes.pop_back();
                if (failed) {
                    return;
                }
            }
        };
        visit(visit, 1);
        if (failed) {
            std::cout << "FAIL tested=" << tested << '\n';
            return EXIT_FAILURE;
        }
        std::cout << "SU2_TURAN_ATOM PASS tested=" << tested
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
