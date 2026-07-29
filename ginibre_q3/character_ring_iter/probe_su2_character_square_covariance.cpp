#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Integer = boost::multiprecision::cpp_int;

int parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0;
    const long long parsed = std::stoll(value, &consumed);
    if (consumed != value.size() || parsed <= 0) {
        throw std::invalid_argument(name + " must be positive");
    }
    return static_cast<int>(parsed);
}

std::vector<Integer> multiply(
    const std::vector<Integer>& left,
    const std::vector<Integer>& right
) {
    std::vector<Integer> product(left.size() + right.size() - 1U);
    for (int a = 0; a < static_cast<int>(left.size()); ++a) {
        for (int b = 0; b < static_cast<int>(right.size()); ++b) {
            for (int c = std::abs(a - b); c <= a + b; ++c) {
                product[static_cast<std::size_t>(c)]
                    += left[static_cast<std::size_t>(a)]
                    * right[static_cast<std::size_t>(b)];
            }
        }
    }
    return product;
}

Integer at(const std::vector<Integer>& profile, int index) {
    return
        index >= 0 && index < static_cast<int>(profile.size())
        ? profile[static_cast<std::size_t>(index)]
        : Integer{0};
}

std::string show(const std::vector<Integer>& profile) {
    std::string result = "{";
    for (std::size_t index = 0; index < profile.size(); ++index) {
        if (index != 0U) {
            result.push_back(',');
        }
        result += profile[index].convert_to<std::string>();
    }
    result.push_back('}');
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        int maximum_support = 7;
        int maximum_coefficient = 3;
        if (argc >= 2) {
            maximum_support =
                parse_positive(argv[1], "maximum_support");
        }
        if (argc >= 3) {
            maximum_coefficient =
                parse_positive(argv[2], "maximum_coefficient");
        }
        if (argc > 3 || maximum_coefficient > 9) {
            throw std::invalid_argument(
                "usage: probe_su2_character_square_covariance"
                " [maximum_support] [maximum_coefficient_at_most_9]"
            );
        }

        unsigned long long profiles = 0;
        unsigned long long determinants = 0;
        Integer minimum = 0;
        std::string witness;
        unsigned long long profile_count = 1;
        for (int index = 0; index <= maximum_support; ++index) {
            profile_count *=
                static_cast<unsigned long long>(maximum_coefficient + 1);
        }
        for (unsigned long long code = 1; code < profile_count; ++code) {
            unsigned long long remaining = code;
            std::vector<Integer> root(
                static_cast<std::size_t>(maximum_support + 1)
            );
            for (Integer& coefficient : root) {
                coefficient =
                    remaining
                    % static_cast<unsigned long long>(
                        maximum_coefficient + 1
                    );
                remaining /=
                    static_cast<unsigned long long>(
                        maximum_coefficient + 1
                    );
            }
            while (!root.empty() && root.back() == 0) {
                root.pop_back();
            }
            if (root.empty()) {
                continue;
            }
            bool log_concave = true;
            for (std::size_t index = 1U;
                 index + 1U < root.size();
                 ++index) {
                if (
                    root[index] * root[index]
                    < root[index - 1U] * root[index + 1U]
                ) {
                    log_concave = false;
                    break;
                }
            }
            if (
                !log_concave
                || std::find(root.begin(), root.end(), Integer{0})
                    != root.end()
            ) {
                continue;
            }
            ++profiles;
            const std::vector<Integer> square = multiply(root, root);
            for (int q = 1; q < static_cast<int>(square.size()); ++q) {
                std::vector<Integer> character(
                    static_cast<std::size_t>(q + 1)
                );
                character[static_cast<std::size_t>(q)] = 1;
                const std::vector<Integer> updated =
                    multiply(square, character);
                for (int target = 1;
                     target < static_cast<int>(square.size());
                     ++target) {
                    const Integer determinant =
                        at(square, 0) * at(updated, target)
                        - at(updated, 0) * at(square, target);
                    ++determinants;
                    if (determinant < minimum) {
                        minimum = determinant;
                        witness =
                            "root=" + show(root)
                            + " q=" + std::to_string(q)
                            + " target=" + std::to_string(target)
                            + " square=" + show(square);
                    }
                }
            }
        }

        std::cout
            << "SU2_CHARACTER_SQUARE_COVARIANCE"
            << " maximum_support=" << maximum_support
            << " maximum_coefficient=" << maximum_coefficient
            << " profiles=" << profiles
            << " determinants=" << determinants
            << " minimum=" << minimum
            << " witness=" << (witness.empty() ? "{}" : witness)
            << " result="
            << (
                minimum >= 0
                ? "PASS_BOUNDED_EXACT_DIAGNOSTIC"
                : "COUNTEREXAMPLE"
            )
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
