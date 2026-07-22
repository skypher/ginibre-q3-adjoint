#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

using Charge = std::pair<int, int>;
using Tensor = std::map<Charge, cpp_int>;

Tensor multiply_fundamental_difference(const Tensor& input) {
    Tensor output;
    for (const auto& [charge, coefficient] : input) {
        const auto [left, right] = charge;
        if (left == 0) {
            output[{1, right}] += coefficient;
        } else {
            output[{left - 1, right}] += coefficient;
            output[{left + 1, right}] += coefficient;
        }
        if (right == 0) {
            output[{left, 1}] -= coefficient;
        } else {
            output[{left, right - 1}] -= coefficient;
            output[{left, right + 1}] -= coefficient;
        }
    }
    return output;
}

void print_contraction(const Tensor& tensor) {
    bool first = true;
    for (const auto& [charge, coefficient] : tensor) {
        if (charge.second != 0 || coefficient == 0) {
            continue;
        }
        std::cout << (first ? " [" : ",") << charge.first << ':'
                  << coefficient;
        first = false;
    }
    std::cout << (first ? " []" : "]");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: search_su2_fundamental_difference_cone "
                "MAXIMUM_CHARGE MAXIMUM_POWER"
            );
        }
        const int maximum_charge = std::stoi(argv[1]);
        const int maximum_power = std::stoi(argv[2]);
        if (maximum_charge < 0 || maximum_power < 0) {
            throw std::runtime_error("invalid bound");
        }
        long long tested = 0;
        for (int left = 0; left <= maximum_charge; ++left) {
            for (int right = left; right <= maximum_charge; ++right) {
                Tensor current;
                current[{left, right}] += 1;
                if (left != right) {
                    current[{right, left}] += 1;
                }
                for (int power = 0; power <= maximum_power; ++power) {
                    ++tested;
                    bool nonnegative = true;
                    for (const auto& [charge, coefficient] : current) {
                        if (charge.second == 0 && coefficient < 0) {
                            nonnegative = false;
                            break;
                        }
                    }
                    if (!nonnegative) {
                        std::cout << "FAIL left=" << left
                                  << " right=" << right
                                  << " power=" << power;
                        print_contraction(current);
                        std::cout << '\n';
                        return EXIT_FAILURE;
                    }
                    current = multiply_fundamental_difference(current);
                }
            }
        }
        std::cout << "SU2_FUNDAMENTAL_DIFFERENCE_ATOM_CONE PASS tested="
                  << tested << " maximum_charge=" << maximum_charge
                  << " maximum_power=" << maximum_power << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
