#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Exponent = std::vector<int>;

struct Polynomial {
    std::map<Exponent, Integer> terms;
};

Polynomial constant(int variables, const Integer& value) {
    Polynomial result;
    if (value != 0) {
        result.terms.emplace(
            Exponent(static_cast<std::size_t>(variables), 0),
            value
        );
    }
    return result;
}

Polynomial variable(int variables, int index) {
    Polynomial result = constant(variables, 1);
    auto node = result.terms.extract(result.terms.begin());
    node.key()[static_cast<std::size_t>(index)] = 1;
    result.terms.insert(std::move(node));
    return result;
}

Polynomial& operator+=(Polynomial& left, const Polynomial& right) {
    for (const auto& [exponent, coefficient] : right.terms) {
        Integer& target = left.terms[exponent];
        target += coefficient;
        if (target == 0) {
            left.terms.erase(exponent);
        }
    }
    return left;
}

Polynomial operator-(Polynomial left, const Polynomial& right) {
    for (const auto& [exponent, coefficient] : right.terms) {
        Integer& target = left.terms[exponent];
        target -= coefficient;
        if (target == 0) {
            left.terms.erase(exponent);
        }
    }
    return left;
}

Polynomial operator*(const Polynomial& left, const Polynomial& right) {
    Polynomial result;
    for (const auto& [left_exponent, left_coefficient] : left.terms) {
        for (
            const auto& [right_exponent, right_coefficient]
            : right.terms
        ) {
            Exponent exponent = left_exponent;
            for (std::size_t index = 0; index < exponent.size(); ++index) {
                exponent[index] += right_exponent[index];
            }
            result.terms[exponent] +=
                left_coefficient * right_coefficient;
        }
    }
    for (auto iterator = result.terms.begin();
         iterator != result.terms.end();) {
        if (iterator->second == 0) {
            iterator = result.terms.erase(iterator);
        } else {
            ++iterator;
        }
    }
    return result;
}

Polynomial at(
    const std::vector<Polynomial>& profile,
    int index,
    int variables
) {
    if (index < 0 || index >= static_cast<int>(profile.size())) {
        return constant(variables, 0);
    }
    return profile[static_cast<std::size_t>(index)];
}

std::vector<Polynomial> transform(
    const std::vector<Polynomial>& profile,
    int label,
    int variables
) {
    std::vector<Polynomial> result(
        profile.size() + static_cast<std::size_t>(label),
        constant(variables, 0)
    );
    for (int target = 0; target < static_cast<int>(result.size());
         ++target) {
        for (int source = std::abs(target - label);
             source <= target + label;
             ++source) {
            result[static_cast<std::size_t>(target)] +=
                at(profile, source, variables);
        }
    }
    return result;
}

Polynomial inner(
    const std::vector<Polynomial>& left,
    const std::vector<Polynomial>& right,
    int cutoff,
    int variables
) {
    Polynomial result = constant(variables, 0);
    const int size = std::max(
        static_cast<int>(left.size()),
        static_cast<int>(right.size())
    );
    for (int index = cutoff; index < size; ++index) {
        result += at(left, index, variables)
            * at(right, index, variables);
    }
    return result;
}

std::vector<Polynomial> profile_from_ratio_drops(
    int core_length,
    int shift
) {
    const int variables = core_length - 1;
    std::vector<Polynomial> ratios;
    ratios.reserve(static_cast<std::size_t>(variables));
    for (int ratio = 0; ratio < variables; ++ratio) {
        Polynomial value = constant(variables, 0);
        for (int drop = ratio; drop < variables; ++drop) {
            value += variable(variables, drop);
        }
        ratios.push_back(std::move(value));
    }
    std::vector<Polynomial> result(
        static_cast<std::size_t>(shift),
        constant(variables, 0)
    );
    Polynomial entry = constant(variables, 1);
    result.push_back(entry);
    for (const Polynomial& ratio : ratios) {
        entry = entry * ratio;
        result.push_back(entry);
    }
    return result;
}

void print_exponent(const Exponent& exponent) {
    std::cout << '{';
    for (std::size_t index = 0; index < exponent.size(); ++index) {
        if (index != 0) {
            std::cout << ',';
        }
        std::cout << exponent[index];
    }
    std::cout << '}';
}

int parse_positive(const char* text, const std::string& name) {
    const std::string value{text};
    std::size_t consumed = 0;
    const long long parsed = std::stoll(value, &consumed);
    if (consumed != value.size() || parsed <= 0) {
        throw std::invalid_argument(name + " must be positive");
    }
    return static_cast<int>(parsed);
}

}  // namespace

int main(int argc, char** argv) {
    int max_length = 5;
    int max_label = 3;
    int max_shift = 3;
    try {
        if (argc >= 2) {
            max_length = parse_positive(argv[1], "max_length");
        }
        if (argc >= 3) {
            max_label = parse_positive(argv[2], "max_label");
        }
        if (argc >= 4) {
            max_shift = parse_positive(argv[3], "max_shift");
        }
        if (argc > 4 || max_length < 2) {
            throw std::invalid_argument(
                "usage: probe_su2_log_concave_ratio_polynomial"
                " [max_length>=2] [max_label] [max_shift]"
            );
        }
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }

    std::size_t determinants = 0;
    std::size_t coefficients = 0;
    for (int length = 2; length <= max_length; ++length) {
        const int variables = length - 1;
        for (int shift = 0; shift <= max_shift; ++shift) {
            const std::vector<Polynomial> profile =
                profile_from_ratio_drops(length, shift);
            std::vector<std::vector<Polynomial>> images;
            images.push_back(profile);
            for (int label = 1; label <= max_label; ++label) {
                images.push_back(transform(profile, label, variables));
            }
            for (int q = 1; q <= max_label; ++q) {
                for (int a = 1; a <= max_label; ++a) {
                    const int size = std::max(
                        static_cast<int>(
                            images[static_cast<std::size_t>(q)].size()
                        ),
                        static_cast<int>(
                            images[static_cast<std::size_t>(a)].size()
                        )
                    );
                    for (int cutoff = 0; cutoff <= size; ++cutoff) {
                        const Polynomial determinant =
                            inner(
                                profile,
                                profile,
                                cutoff,
                                variables
                            )
                                * inner(
                                    images[static_cast<std::size_t>(q)],
                                    images[static_cast<std::size_t>(a)],
                                    cutoff,
                                    variables
                                )
                            - inner(
                                  profile,
                                  images[static_cast<std::size_t>(q)],
                                  cutoff,
                                  variables
                              )
                                * inner(
                                    profile,
                                    images[static_cast<std::size_t>(a)],
                                    cutoff,
                                    variables
                                );
                        ++determinants;
                        coefficients += determinant.terms.size();
                        for (
                            const auto& [exponent, coefficient]
                            : determinant.terms
                        ) {
                            if (coefficient < 0) {
                                std::cout
                                    << "SU2_LOG_CONCAVE_RATIO_POLYNOMIAL"
                                    << " negative_coefficient"
                                    << " length=" << length
                                    << " shift=" << shift
                                    << " q=" << q
                                    << " a=" << a
                                    << " cutoff=" << cutoff
                                    << " coefficient=" << coefficient
                                    << " exponent=";
                                print_exponent(exponent);
                                std::cout << '\n';
                                return EXIT_SUCCESS;
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout
        << "SU2_LOG_CONCAVE_RATIO_POLYNOMIAL"
        << " determinants=" << determinants
        << " coefficients=" << coefficients
        << " max_length=" << max_length
        << " max_label=" << max_label
        << " max_shift=" << max_shift
        << " result=NO_NEGATIVE_COEFFICIENT"
        << '\n';
    return EXIT_SUCCESS;
}
