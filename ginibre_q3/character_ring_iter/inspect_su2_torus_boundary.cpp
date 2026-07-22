#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

using Linear = std::map<int, long long>;
using Quadratic = std::map<std::pair<int, int>, long long>;

void add_linear(Linear& target, int charge, long long coefficient, int parity) {
    if (charge == 0 && parity < 0) {
        return;
    }
    if (charge < 0) {
        charge = -charge;
        coefficient *= static_cast<long long>(parity);
    }
    target[charge] += coefficient;
}

void add_product(
    Quadratic& target,
    const Linear& left,
    const Linear& right,
    long long scale
) {
    for (const auto& [left_charge, left_coefficient] : left) {
        for (const auto& [right_charge, right_coefficient] : right) {
            const auto key = std::minmax(left_charge, right_charge);
            target[key] += scale * left_coefficient * right_coefficient;
        }
    }
}

Linear multiplied_coefficient(int output, int feature, int sign) {
    Linear result;
    add_linear(result, output - feature, 1, sign);
    add_linear(result, output + feature, sign, sign);
    return result;
}

void add_q(
    Quadratic& target,
    const Linear& coefficient_zero,
    const Linear& coefficient_two,
    const Linear& coefficient_four,
    long long scale
) {
    add_product(target, coefficient_zero, coefficient_zero, scale);
    add_product(target, coefficient_two, coefficient_two, -2 * scale);
    add_product(target, coefficient_zero, coefficient_four, scale);
}

Quadratic boundary_form(int label, int sign) {
    Quadratic result;
    if (sign > 0 && (label & 1) == 0) {
        Linear zero{{0, 1}};
        Linear two{{2, 1}};
        Linear four{{4, 1}};
        add_q(result, zero, two, four, 2);
    }
    for (int feature = label; feature > 0; feature -= 2) {
        add_q(
            result,
            multiplied_coefficient(0, feature, sign),
            multiplied_coefficient(2, feature, sign),
            multiplied_coefficient(4, feature, sign),
            1
        );
    }
    for (auto iterator = result.begin(); iterator != result.end();) {
        if (iterator->second == 0) {
            iterator = result.erase(iterator);
        } else {
            ++iterator;
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: inspect_su2_torus_boundary LABEL plus|minus"
            );
        }
        const int label = std::stoi(argv[1]);
        if (label < 1) {
            throw std::runtime_error("LABEL must be positive");
        }
        const std::string mode = argv[2];
        const int sign = mode == "plus" ? 1 : mode == "minus" ? -1 : 0;
        if (sign == 0) {
            throw std::runtime_error("expected plus or minus");
        }
        const Quadratic form = boundary_form(label, sign);
        std::cout << "label=" << label << " mode=" << mode << '\n';
        for (const auto& [charges, coefficient] : form) {
            std::cout << coefficient << " b[" << charges.first << "]b["
                      << charges.second << "]\n";
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
