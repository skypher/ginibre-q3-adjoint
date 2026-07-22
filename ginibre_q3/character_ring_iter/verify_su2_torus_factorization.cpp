#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

using boost::multiprecision::cpp_int;

namespace {

struct Factor {
    int label;
    bool minus;
};

struct TorusStats {
    static constexpr std::size_t profile_size = 7U;

    cpp_int b00 = 0;
    cpp_int b22 = 0;
    cpp_int b04 = 0;
    std::array<cpp_int, profile_size> diagonal{};
    std::array<cpp_int, profile_size> offset_four{};

    cpp_int value() const {
        return b00 - 2 * b22 + b04;
    }
};

std::vector<std::vector<int>> multisets(int maximum_label, int maximum_size) {
    std::vector<std::vector<int>> result;
    std::vector<int> current;
    const auto visit = [&](const auto& self, int first, int remaining) -> void {
        result.push_back(current);
        if (remaining == 0) {
            return;
        }
        for (int label = first; label <= maximum_label; ++label) {
            current.push_back(label);
            self(self, label, remaining - 1);
            current.pop_back();
        }
    };
    visit(visit, 1, maximum_size);
    return result;
}

cpp_int invariant_multiplicity(const std::vector<int>& labels) {
    int total = 0;
    for (int label : labels) {
        total += label;
    }
    std::vector<cpp_int> current(static_cast<std::size_t>(total + 1));
    current[0] = 1;
    int support = 0;
    for (int label : labels) {
        std::vector<cpp_int> next(static_cast<std::size_t>(total + 1));
        for (int input = 0; input <= support; ++input) {
            const cpp_int& value = current[static_cast<std::size_t>(input)];
            if (value == 0) {
                continue;
            }
            for (int output = std::abs(input - label);
                 output <= input + label;
                 output += 2) {
                next[static_cast<std::size_t>(output)] += value;
            }
        }
        support += label;
        current = std::move(next);
    }
    return current[0];
}

cpp_int subset_value(const std::vector<Factor>& factors) {
    if (factors.size() >= 63U) {
        throw std::runtime_error("too many factors for subset cross-check");
    }
    const std::uint64_t count = std::uint64_t{1} << factors.size();
    cpp_int answer = 0;
    std::vector<int> left;
    std::vector<int> right;
    for (std::uint64_t mask = 0; mask < count; ++mask) {
        left.clear();
        right.clear();
        bool negative = false;
        for (std::size_t index = 0; index < factors.size(); ++index) {
            if ((mask & (std::uint64_t{1} << index)) != 0U) {
                left.push_back(factors[index].label);
                negative = negative != factors[index].minus;
            } else {
                right.push_back(factors[index].label);
            }
        }
        const cpp_int term = invariant_multiplicity(left)
            * invariant_multiplicity(right);
        answer += negative ? -term : term;
    }
    return answer;
}

cpp_int coefficient(
    const std::vector<cpp_int>& polynomial,
    int offset,
    int exponent
) {
    const int index = offset + exponent;
    if (index < 0 || index >= static_cast<int>(polynomial.size())) {
        return 0;
    }
    return polynomial[static_cast<std::size_t>(index)];
}

void accumulate_features(
    const std::vector<Factor>& factors,
    std::size_t index,
    int offset,
    const std::vector<cpp_int>& polynomial,
    const cpp_int& multiplicity,
    TorusStats& stats
) {
    if (index == factors.size()) {
        const cpp_int c0 = coefficient(polynomial, offset, 0);
        const cpp_int c2 = coefficient(polynomial, offset, 2);
        const cpp_int c4 = coefficient(polynomial, offset, 4);
        stats.b00 += multiplicity * c0 * c0;
        stats.b22 += multiplicity * c2 * c2;
        stats.b04 += multiplicity * c0 * c4;
        for (std::size_t profile_index = 0U;
             profile_index < TorusStats::profile_size;
             ++profile_index) {
            const int charge = 2 * static_cast<int>(profile_index);
            const cpp_int lower = coefficient(polynomial, offset, charge);
            const cpp_int upper = coefficient(
                polynomial, offset, charge + 4
            );
            stats.diagonal[profile_index]
                += multiplicity * lower * lower;
            stats.offset_four[profile_index]
                += multiplicity * lower * upper;
        }
        return;
    }

    const Factor factor = factors[index];
    if (!factor.minus && (factor.label & 1) == 0) {
        accumulate_features(
            factors, index + 1U, offset, polynomial,
            2 * multiplicity, stats
        );
    }
    for (int weight = factor.label; weight > 0; weight -= 2) {
        std::vector<cpp_int> next(polynomial.size());
        for (int exponent = -offset; exponent <= offset; ++exponent) {
            const cpp_int& value = polynomial[static_cast<std::size_t>(
                offset + exponent
            )];
            if (value == 0) {
                continue;
            }
            next[static_cast<std::size_t>(offset + exponent + weight)] += value;
            const cpp_int lower_value = factor.minus ? -value : value;
            next[static_cast<std::size_t>(offset + exponent - weight)]
                += lower_value;
        }
        accumulate_features(
            factors, index + 1U, offset, next, multiplicity, stats
        );
    }
}

TorusStats torus_stats(const std::vector<Factor>& factors) {
    int total = 0;
    for (const Factor& factor : factors) {
        total += factor.label;
    }
    std::vector<cpp_int> polynomial(static_cast<std::size_t>(2 * total + 1));
    polynomial[static_cast<std::size_t>(total)] = 1;
    TorusStats stats;
    accumulate_features(factors, 0U, total, polynomial, 1, stats);
    return stats;
}

void print_case(const std::vector<Factor>& factors) {
    std::cout << "minus=[";
    bool first = true;
    for (const Factor& factor : factors) {
        if (!factor.minus) {
            continue;
        }
        if (!first) {
            std::cout << ',';
        }
        first = false;
        std::cout << factor.label;
    }
    std::cout << "] plus=[";
    first = true;
    for (const Factor& factor : factors) {
        if (factor.minus) {
            continue;
        }
        if (!first) {
            std::cout << ',';
        }
        first = false;
        std::cout << factor.label;
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 3 || argc > 4) {
            throw std::runtime_error(
                "usage: verify_su2_torus_factorization "
                "MAXIMUM_LABEL MAXIMUM_FACTORS [--analyze]"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 1 || maximum_factors < 0) {
            throw std::runtime_error("invalid search bound");
        }
        const bool analyze = argc == 4
            && std::string(argv[3]) == "--analyze";
        if (argc == 4 && !analyze) {
            throw std::runtime_error("unknown mode");
        }
        const auto words = multisets(maximum_label, maximum_factors);
        struct SignedWord {
            std::size_t minus_index;
            std::size_t plus_index;
        };
        std::vector<SignedWord> cases;
        for (std::size_t minus_index = 0;
             minus_index < words.size();
             ++minus_index) {
            if ((words[minus_index].size() & 1U) != 0U) {
                continue;
            }
            for (std::size_t plus_index = 0;
                 plus_index < words.size();
                 ++plus_index) {
                if (words[minus_index].size() + words[plus_index].size()
                    <= static_cast<std::size_t>(maximum_factors)) {
                    cases.push_back({minus_index, plus_index});
                }
            }
        }

        std::size_t first_failure = cases.size();
        cpp_int failed_subset = 0;
        cpp_int failed_torus = 0;
        std::size_t first_b00_minus_b22 = cases.size();
        std::size_t first_b00_minus_2b22 = cases.size();
        std::size_t first_b04 = cases.size();
        std::size_t first_negative_d_with_minus = cases.size();
        std::size_t first_negative_b04_with_b22 = cases.size();
        std::size_t first_both_negative = cases.size();
        std::size_t first_radial_failure = cases.size();
        std::size_t first_shifted_q_failure = cases.size();
        std::size_t first_turan_negative = cases.size();
        std::size_t failed_profile_index = TorusStats::profile_size;
        std::size_t failed_q_profile_index = TorusStats::profile_size;
        std::size_t failed_turan_profile_index = TorusStats::profile_size;

#pragma omp parallel
        {
            std::size_t local_failure = cases.size();
            cpp_int local_subset = 0;
            cpp_int local_torus = 0;
            std::size_t local_b00_minus_b22 = cases.size();
            std::size_t local_b00_minus_2b22 = cases.size();
            std::size_t local_b04 = cases.size();
            std::size_t local_negative_d_with_minus = cases.size();
            std::size_t local_negative_b04_with_b22 = cases.size();
            std::size_t local_both_negative = cases.size();
            std::size_t local_radial_failure = cases.size();
            std::size_t local_shifted_q_failure = cases.size();
            std::size_t local_turan_negative = cases.size();
            std::size_t local_profile_index = TorusStats::profile_size;
            std::size_t local_q_profile_index = TorusStats::profile_size;
            std::size_t local_turan_profile_index
                = TorusStats::profile_size;

#pragma omp for schedule(dynamic)
            for (std::size_t case_index = 0;
                 case_index < cases.size();
                 ++case_index) {
                const SignedWord signed_word = cases[case_index];
                std::vector<Factor> factors;
                for (int label : words[signed_word.minus_index]) {
                    factors.push_back({label, true});
                }
                for (int label : words[signed_word.plus_index]) {
                    factors.push_back({label, false});
                }
                const cpp_int subset = subset_value(factors);
                const TorusStats stats = torus_stats(factors);
                const cpp_int torus = stats.value();
                if (subset != torus && case_index < local_failure) {
                    local_failure = case_index;
                    local_subset = subset;
                    local_torus = torus;
                }
                if (analyze) {
                    if (stats.b00 < stats.b22
                        && case_index < local_b00_minus_b22) {
                        local_b00_minus_b22 = case_index;
                    }
                    if (stats.b00 < 2 * stats.b22
                        && case_index < local_b00_minus_2b22) {
                        local_b00_minus_2b22 = case_index;
                    }
                    if (stats.b04 < 0 && case_index < local_b04) {
                        local_b04 = case_index;
                    }
                    const bool negative_d = stats.b00 < 2 * stats.b22;
                    if (negative_d
                        && !words[signed_word.minus_index].empty()
                        && case_index < local_negative_d_with_minus) {
                        local_negative_d_with_minus = case_index;
                    }
                    if (stats.b04 < 0 && stats.b22 > 0
                        && case_index < local_negative_b04_with_b22) {
                        local_negative_b04_with_b22 = case_index;
                    }
                    if (negative_d && stats.b04 < 0
                        && case_index < local_both_negative) {
                        local_both_negative = case_index;
                    }
                    for (std::size_t profile_index = 0U;
                         profile_index + 1U < TorusStats::profile_size;
                         ++profile_index) {
                        if (stats.diagonal[profile_index]
                                < stats.diagonal[profile_index + 1U]
                            && case_index < local_radial_failure) {
                            local_radial_failure = case_index;
                            local_profile_index = profile_index;
                        }
                        const cpp_int shifted_q
                            = stats.diagonal[profile_index]
                            - 2 * stats.diagonal[profile_index + 1U]
                            + stats.offset_four[profile_index];
                        if (shifted_q < 0
                            && case_index < local_shifted_q_failure) {
                            local_shifted_q_failure = case_index;
                            local_q_profile_index = profile_index;
                        }
                        const cpp_int turan_value
                            = stats.diagonal[profile_index + 1U]
                            - stats.offset_four[profile_index];
                        if (turan_value < 0
                            && case_index < local_turan_negative) {
                            local_turan_negative = case_index;
                            local_turan_profile_index = profile_index;
                        }
                    }
                }
            }

#pragma omp critical
            {
                if (local_failure < first_failure) {
                    first_failure = local_failure;
                    failed_subset = local_subset;
                    failed_torus = local_torus;
                }
                first_b00_minus_b22 = std::min(
                    first_b00_minus_b22, local_b00_minus_b22
                );
                first_b00_minus_2b22 = std::min(
                    first_b00_minus_2b22, local_b00_minus_2b22
                );
                first_b04 = std::min(
                    first_b04, local_b04
                );
                first_negative_d_with_minus = std::min(
                    first_negative_d_with_minus, local_negative_d_with_minus
                );
                first_negative_b04_with_b22 = std::min(
                    first_negative_b04_with_b22,
                    local_negative_b04_with_b22
                );
                first_both_negative = std::min(
                    first_both_negative, local_both_negative
                );
                if (local_radial_failure < first_radial_failure) {
                    first_radial_failure = local_radial_failure;
                    failed_profile_index = local_profile_index;
                }
                if (local_shifted_q_failure < first_shifted_q_failure) {
                    first_shifted_q_failure = local_shifted_q_failure;
                    failed_q_profile_index = local_q_profile_index;
                }
                if (local_turan_negative < first_turan_negative) {
                    first_turan_negative = local_turan_negative;
                    failed_turan_profile_index
                        = local_turan_profile_index;
                }
            }
        }

        std::cout << "SU2_TORUS_FACTORIZATION tested=" << cases.size()
                  << " maximum_label=" << maximum_label
                  << " maximum_factors=" << maximum_factors << '\n';
        if (first_failure != cases.size()) {
            const SignedWord signed_word = cases[first_failure];
            std::vector<Factor> factors;
            for (int label : words[signed_word.minus_index]) {
                factors.push_back({label, true});
            }
            for (int label : words[signed_word.plus_index]) {
                factors.push_back({label, false});
            }
            std::cout << "FAIL ";
            print_case(factors);
            std::cout << " subset=" << failed_subset
                      << " torus=" << failed_torus << '\n';
            return EXIT_FAILURE;
        }
        if (analyze) {
            const auto report = [&](const char* name, std::size_t index) {
                std::cout << name;
                if (index == cases.size()) {
                    std::cout << " NO_NEGATIVE\n";
                    return;
                }
                const SignedWord signed_word = cases[index];
                std::vector<Factor> factors;
                for (int label : words[signed_word.minus_index]) {
                    factors.push_back({label, true});
                }
                for (int label : words[signed_word.plus_index]) {
                    factors.push_back({label, false});
                }
                const TorusStats stats = torus_stats(factors);
                std::cout << " NEGATIVE ";
                print_case(factors);
                std::cout << " b00=" << stats.b00
                          << " b22=" << stats.b22
                          << " b04=" << stats.b04 << '\n';
            };
            report("B00_MINUS_B22", first_b00_minus_b22);
            report("B00_MINUS_2B22", first_b00_minus_2b22);
            report("B04", first_b04);
            report("NEGATIVE_D_WITH_MINUS", first_negative_d_with_minus);
            report(
                "NEGATIVE_B04_WITH_POSITIVE_B22",
                first_negative_b04_with_b22
            );
            report("BOTH_NEGATIVE", first_both_negative);
            const auto report_profile = [&](const char* name,
                                            std::size_t index,
                                            std::size_t profile_index,
                                            bool shifted_q) {
                std::cout << name;
                if (index == cases.size()) {
                    std::cout << " NO_NEGATIVE\n";
                    return;
                }
                const SignedWord signed_word = cases[index];
                std::vector<Factor> factors;
                for (int label : words[signed_word.minus_index]) {
                    factors.push_back({label, true});
                }
                for (int label : words[signed_word.plus_index]) {
                    factors.push_back({label, false});
                }
                const TorusStats stats = torus_stats(factors);
                const int charge = 2 * static_cast<int>(profile_index);
                std::cout << " NEGATIVE ";
                print_case(factors);
                std::cout << " charge=" << charge
                          << " diagonal=" << stats.diagonal[profile_index]
                          << " next_diagonal="
                          << stats.diagonal[profile_index + 1U];
                if (shifted_q) {
                    std::cout << " offset_four="
                              << stats.offset_four[profile_index];
                }
                std::cout << '\n';
            };
            report_profile(
                "RADIAL_DIAGONAL", first_radial_failure,
                failed_profile_index, false
            );
            report_profile(
                "SHIFTED_Q", first_shifted_q_failure,
                failed_q_profile_index, true
            );
            report_profile(
                "TURAN_VALUE", first_turan_negative,
                failed_turan_profile_index, true
            );
        }
        std::cout << "PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
