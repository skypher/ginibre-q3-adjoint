#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

using Laurent = std::map<int, long long>;

Laurent multiply(const Laurent& left, const Laurent& right) {
    Laurent product;
    for (const auto& [left_degree, left_coefficient] : left) {
        for (const auto& [right_degree, right_coefficient] : right) {
            product[left_degree + right_degree]
                += left_coefficient * right_coefficient;
        }
    }
    return product;
}

Laurent c_feature(const int degree) {
    if (degree <= 0) {
        throw std::runtime_error("C feature degree must be positive");
    }
    return {{-degree, 1}, {degree, 1}};
}

Laurent s_feature(const int degree) {
    if (degree <= 0) {
        throw std::runtime_error("S feature degree must be positive");
    }
    return {{-degree, -1}, {degree, 1}};
}

long long coefficient(const Laurent& polynomial, const int degree) {
    const auto iterator = polynomial.find(degree);
    return iterator == polynomial.end() ? 0 : iterator->second;
}

long long q_form(const Laurent& polynomial) {
    const long long a0 = coefficient(polynomial, 0);
    const long long a2 = coefficient(polynomial, 2);
    const long long a4 = coefficient(polynomial, 4);
    return a0 * a0 - 2 * a2 * a2 + a0 * a4;
}

Laurent terminal_feature(
    const int first_minus,
    const int second_minus,
    const int target_feature
) {
    return multiply(
        multiply(s_feature(first_minus), s_feature(second_minus)),
        multiply(c_feature(1), c_feature(target_feature))
    );
}

long long frozen_top_sum(const int target) {
    long long total = 0;
    for (const int first_minus : {1, 3}) {
        for (const int second_minus : {1, 3}) {
            for (int target_feature = 1;
                 target_feature <= target;
                 target_feature += 2) {
                total += q_form(multiply(
                    terminal_feature(
                        first_minus, second_minus, target_feature
                    ),
                    c_feature(target + 3)
                ));
            }
        }
    }
    return total;
}

long long background_feature_sum(
    const int target,
    const int background_feature
) {
    long long total = 0;
    for (const int first_minus : {1, 3}) {
        for (const int second_minus : {1, 3}) {
            for (int target_feature = 1;
                 target_feature <= target;
                 target_feature += 2) {
                const Laurent terminal = terminal_feature(
                    first_minus, second_minus, target_feature
                );
                total += background_feature == 0
                    ? 2 * q_form(terminal)
                    : q_form(multiply(
                        terminal, c_feature(background_feature)
                    ));
            }
        }
    }
    return total;
}

}  // namespace

int main() {
    try {
        for (int target = 5; target <= 501; target += 2) {
            const long long value = frozen_top_sum(target);
            if (value != -6) {
                std::cout << "SU2_TERMINAL_FEATURE_FREEZE FAIL target="
                          << target << " value=" << value << '\n';
                return EXIT_FAILURE;
            }
        }

        const std::vector<long long> expected{8, 2, 6, 0, -6};
        std::vector<long long> observed;
        observed.push_back(background_feature_sum(5, 0));
        for (int feature = 2; feature <= 8; feature += 2) {
            observed.push_back(background_feature_sum(5, feature));
        }
        if (observed != expected) {
            std::cout << "SU2_TERMINAL_FEATURE_FREEZE FAIL breakdown";
            for (const long long value : observed) {
                std::cout << ' ' << value;
            }
            std::cout << '\n';
            return EXIT_FAILURE;
        }

        long long complete_sum = 0;
        for (const long long value : observed) {
            complete_sum += value;
        }
        std::cout << "SU2_TERMINAL_FEATURE_FREEZE PASS"
                  << " odd_targets=249"
                  << " frozen_top=-6"
                  << " target5_breakdown=8,2,6,0,-6"
                  << " complete_sum=" << complete_sum << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << "error: " << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
