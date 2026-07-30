#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using boost::multiprecision::cpp_int;

namespace {

using Vector = std::vector<cpp_int>;

int positive_argument(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0U;
    const long long parsed = std::stoll(value, &consumed, 10);
    if (consumed != value.size() || parsed <= 0
        || parsed > std::numeric_limits<int>::max()) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<int>(parsed);
}

cpp_int at(const Vector& values, int index) {
    if (index < 0 || index >= static_cast<int>(values.size())) {
        return 0;
    }
    return values[static_cast<std::size_t>(index)];
}

Vector multiply_by_square(const Vector& input, int factor) {
    const int maximum_source = static_cast<int>(input.size()) - 1;
    Vector output(
        static_cast<std::size_t>(maximum_source + 2 * factor + 1),
        0);
    for (int source = 0; source <= maximum_source; ++source) {
        for (int character = 0; character <= 2 * factor; ++character) {
            for (int target = std::abs(source - character);
                 target <= source + character;
                 ++target) {
                output[static_cast<std::size_t>(target)]
                    += input[static_cast<std::size_t>(source)];
            }
        }
    }
    return output;
}

std::string render(const Vector& values) {
    std::string output = "[";
    for (std::size_t index = 0U; index < values.size(); ++index) {
        if (index != 0U) {
            output += ',';
        }
        output += values[index].convert_to<std::string>();
    }
    return output + ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const int maximum_factor = argc >= 2
            ? positive_argument(argv[1], "maximum_factor")
            : 50;
        const int maximum_column = argc >= 3
            ? positive_argument(argv[2], "maximum_column")
            : 200;
        const int display_factor = argc >= 4
            ? positive_argument(argv[3], "display_factor")
            : 2;
        if (argc > 4 || display_factor > maximum_factor) {
            throw std::invalid_argument(
                "usage: analyze_su2_normalized_dimension_invariance "
                "[maximum_factor] [maximum_column] [display_factor]");
        }

        std::uint64_t columns = 0U;
        std::uint64_t coefficient_checks = 0U;
        std::uint64_t coefficient_failures = 0U;
        int first_factor = -1;
        int first_column = -1;
        int first_row = -1;
        cpp_int first_value = 0;

        for (int factor = 1; factor <= maximum_factor; ++factor) {
            for (int column = 0;
                 column <= maximum_column;
                 ++column) {
                Vector input(
                    static_cast<std::size_t>(column + 1),
                    0);
                for (int label = 0; label <= column; ++label) {
                    input[static_cast<std::size_t>(label)]
                        = 2 * label + 1;
                }
                const Vector product
                    = multiply_by_square(input, factor);
                Vector output_margin(product.size(), 0);
                for (int row = 0;
                     row < static_cast<int>(product.size());
                     ++row) {
                    output_margin[static_cast<std::size_t>(row)]
                        = (2 * row + 3) * at(product, row)
                          - (2 * row + 1) * at(product, row + 1);
                    ++coefficient_checks;
                    if (output_margin[
                            static_cast<std::size_t>(row)] < 0) {
                        ++coefficient_failures;
                        if (first_factor < 0) {
                            first_factor = factor;
                            first_column = column;
                            first_row = row;
                            first_value = output_margin[
                                static_cast<std::size_t>(row)];
                        }
                    }
                }
                if (factor == display_factor
                    && column <= 2 * display_factor + 2) {
                    std::cout
                        << "DETAIL"
                        << " factor=" << factor
                        << " column=" << column
                        << " numerator=" << render(output_margin)
                        << " denominator="
                        << (2 * column + 1) * (2 * column + 3)
                        << '\n';
                }
                ++columns;
            }
        }

        std::cout
            << "SU2_NORMALIZED_DIMENSION_INVARIANCE"
            << " maximum_factor=" << maximum_factor
            << " maximum_column=" << maximum_column
            << " columns=" << columns
            << " coefficient_checks=" << coefficient_checks
            << " coefficient_failures=" << coefficient_failures
            << '\n'
            << "FIRST_FAILURE"
            << " factor=" << first_factor
            << " column=" << first_column
            << " row=" << first_row
            << " value=" << first_value
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
