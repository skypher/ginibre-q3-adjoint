#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

std::vector<std::vector<int>> multisets(int maximum_label, int maximum_size) {
    std::vector<std::vector<int>> result;
    std::vector<int> current;
    const auto visit = [&](const auto& self, int first) -> void {
        if (!current.empty()) {
            result.push_back(current);
        }
        if (current.size() == static_cast<std::size_t>(maximum_size)) {
            return;
        }
        for (int label = first; label <= maximum_label; ++label) {
            current.push_back(label);
            self(self, label);
            current.pop_back();
        }
    };
    visit(visit, 1);
    return result;
}

cpp_int partial_coefficient(
    const std::vector<int>& labels,
    std::uint64_t minus_mask,
    int target
) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    const int width = total + 1;
    std::vector<cpp_int> current(
        static_cast<std::size_t>(width * width)
    );
    current[0] = 1;
    int support = 0;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        const int label = labels[index];
        const bool minus = ((minus_mask >> index) & 1U) != 0U;
        std::vector<cpp_int> next(current.size());
        for (int left = 0; left <= support; ++left) {
            for (int right = 0; right <= support; ++right) {
                const cpp_int& value = current[static_cast<std::size_t>(
                    left * width + right
                )];
                if (value == 0) {
                    continue;
                }
                for (int output = std::abs(left - label);
                     output <= left + label;
                     output += 2) {
                    next[static_cast<std::size_t>(
                        output * width + right
                    )] += value;
                }
                for (int output = std::abs(right - label);
                     output <= right + label;
                     output += 2) {
                    cpp_int& destination = next[static_cast<std::size_t>(
                        left * width + output
                    )];
                    destination += minus ? -value : value;
                }
            }
        }
        support += label;
        current = std::move(next);
    }
    return current[static_cast<std::size_t>(target * width)];
}

long long elementary(const std::vector<int>& signs, int degree) {
    std::vector<long long> coefficients(
        static_cast<std::size_t>(degree + 1), 0
    );
    coefficients[0] = 1;
    int used = 0;
    for (int sign : signs) {
        const int upper = std::min(used + 1, degree);
        for (int index = upper; index >= 1; --index) {
            coefficients[static_cast<std::size_t>(index)]
                += static_cast<long long>(sign)
                    * coefficients[static_cast<std::size_t>(index - 1)];
        }
        ++used;
    }
    return coefficients[static_cast<std::size_t>(degree)];
}

cpp_int closed_four_defect(
    const std::vector<int>& labels,
    std::uint64_t minus_mask
) {
    std::vector<int> signs_one;
    std::vector<int> signs_two;
    int at_least_three = 0;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        const int sign = ((minus_mask >> index) & 1U) != 0U ? -1 : 1;
        if (labels[index] == 1) {
            signs_one.push_back(sign);
        } else if (labels[index] == 2) {
            signs_two.push_back(sign);
        } else {
            ++at_least_three;
        }
    }
    const long long n = static_cast<long long>(labels.size());
    const long long e2 = elementary(signs_one, 2);
    const long long e4 = elementary(signs_one, 4);
    const long long y = std::accumulate(
        signs_two.begin(), signs_two.end(), 0LL
    );
    const long long f2 = elementary(signs_two, 2);
    const long long label_at_least_two
        = static_cast<long long>(signs_two.size()) + at_least_three;
    return n * (n - 3) / 2 + label_at_least_two
        + (n - 3 + y) * e2 + f2 + 2 * e4;
}

cpp_int closed_six_defect(
    const std::vector<int>& labels,
    std::uint64_t minus_mask
) {
    std::vector<int> signs_one;
    std::vector<int> signs_two;
    std::vector<int> signs_three;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        const int sign = ((minus_mask >> index) & 1U) != 0U ? -1 : 1;
        if (labels[index] == 1) {
            signs_one.push_back(sign);
        } else if (labels[index] == 2) {
            signs_two.push_back(sign);
        } else if (labels[index] == 3) {
            signs_three.push_back(sign);
        }
    }

    const long long n = static_cast<long long>(labels.size());
    const long long r = static_cast<long long>(signs_one.size());
    const long long t = static_cast<long long>(signs_two.size());
    const long long a1 = elementary(signs_one, 1);
    const long long a2 = elementary(signs_one, 2);
    const long long a3 = elementary(signs_one, 3);
    const long long a4 = elementary(signs_one, 4);
    const long long a6 = elementary(signs_one, 6);
    const long long b1 = elementary(signs_two, 1);
    const long long b2 = elementary(signs_two, 2);
    const long long b3 = elementary(signs_two, 3);
    const long long c1 = elementary(signs_three, 1);
    const long long c2 = elementary(signs_three, 2);

    const long long empty = (n + 1) * n * (n - 1) / 6
        - r * (n - 1) - t;
    const long long after_two
        = (n - 2) * (n - 5) / 2 + n - r;
    return empty + after_two * a2
        + (n - 3) * b2 + (n - 4) * a2 * b1
        + 2 * (n - 5) * a4 + c2 + a1 * b1 * c1
        + b3 + a3 * c1 + 2 * a2 * b2
        + 3 * a4 * b1 + 5 * a6;
}

void print_case(
    const std::vector<int>& labels,
    std::uint64_t minus_mask,
    int target
) {
    std::cout << "target=" << target << " factors=[";
    for (std::size_t index = 0; index < labels.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << (((minus_mask >> index) & 1U) != 0U ? '-' : '+')
                  << labels[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: verify_su2_outer_layers MAXIMUM_LABEL "
                "MAXIMUM_FACTORS"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 1 || maximum_factors < 1
            || maximum_factors >= 63) {
            throw std::runtime_error("invalid search bound");
        }
        const auto words = multisets(maximum_label, maximum_factors);
        std::uint64_t tested_four = 0;
        std::uint64_t tested_six = 0;
        for (const auto& labels : words) {
            const int sum = std::accumulate(labels.begin(), labels.end(), 0);
            const std::uint64_t masks
                = std::uint64_t{1} << labels.size();
            for (int defect : {4, 6}) {
                const int target = sum - defect;
                if (target < 1 || target < labels.back()) {
                    continue;
                }
                for (std::uint64_t mask = 0; mask < masks; ++mask) {
                    const cpp_int exact = partial_coefficient(
                        labels, mask, target
                    );
                    const cpp_int closed = defect == 4
                        ? closed_four_defect(labels, mask)
                        : closed_six_defect(labels, mask);
                    if (defect == 4) {
                        ++tested_four;
                    } else {
                        ++tested_six;
                    }
                    if (exact != closed || exact < 0) {
                        std::cout << "FAIL defect=" << defect << ' ';
                        print_case(labels, mask, target);
                        std::cout << " exact=" << exact
                                  << " closed=" << closed << '\n';
                        return EXIT_FAILURE;
                    }
                }
            }
        }
        std::cout << "SU2_OUTER_LAYERS PASS tested_four=" << tested_four
                  << " tested_six=" << tested_six
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
