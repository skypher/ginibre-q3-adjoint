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

Linear multiplied_coefficient(
    int output,
    int feature,
    int local_sign,
    int base_parity
) {
    Linear result;
    add_linear(result, output - feature, 1, base_parity);
    add_linear(result, output + feature, local_sign, base_parity);
    return result;
}

void add_turan(
    Quadratic& target,
    int charge,
    int feature,
    int local_sign,
    int base_parity,
    long long scale
) {
    const Linear center = multiplied_coefficient(
        charge, feature, local_sign, base_parity
    );
    const Linear lower = multiplied_coefficient(
        charge - 2, feature, local_sign, base_parity
    );
    const Linear upper = multiplied_coefficient(
        charge + 2, feature, local_sign, base_parity
    );
    add_product(target, center, center, scale);
    add_product(target, lower, upper, -scale);
}

void add_identity_turan(
    Quadratic& target,
    int charge,
    int base_parity,
    long long scale
) {
    Linear center;
    Linear lower;
    Linear upper;
    add_linear(center, charge, 1, base_parity);
    add_linear(lower, charge - 2, 1, base_parity);
    add_linear(upper, charge + 2, 1, base_parity);
    add_product(target, center, center, scale);
    add_product(target, lower, upper, -scale);
}

Quadratic update_form(
    int outer_label,
    int local_label,
    int local_sign,
    int base_parity
) {
    Quadratic result;
    if (local_sign > 0 && (local_label & 1) == 0) {
        add_identity_turan(result, outer_label, base_parity, 2);
        add_identity_turan(result, outer_label + 2, base_parity, -2);
    }
    for (int feature = local_label; feature > 0; feature -= 2) {
        add_turan(
            result, outer_label, feature, local_sign, base_parity, 1
        );
        add_turan(
            result, outer_label + 2, feature, local_sign, base_parity, -1
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
        if (argc != 5) {
            throw std::runtime_error(
                "usage: inspect_su2_turan_update OUTER LOCAL "
                "plus|minus symmetric|antisymmetric"
            );
        }
        const int outer_label = std::stoi(argv[1]);
        const int local_label = std::stoi(argv[2]);
        const std::string local_mode = argv[3];
        const std::string base_mode = argv[4];
        const int local_sign = local_mode == "plus" ? 1
            : local_mode == "minus" ? -1 : 0;
        const int base_parity = base_mode == "symmetric" ? 1
            : base_mode == "antisymmetric" ? -1 : 0;
        if (outer_label < local_label || local_label < 1
            || local_sign == 0 || base_parity == 0) {
            throw std::runtime_error("invalid arguments");
        }
        const Quadratic form = update_form(
            outer_label, local_label, local_sign, base_parity
        );
        std::cout << "outer=" << outer_label << " local=" << local_label
                  << " local_mode=" << local_mode
                  << " base_mode=" << base_mode << '\n';
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
