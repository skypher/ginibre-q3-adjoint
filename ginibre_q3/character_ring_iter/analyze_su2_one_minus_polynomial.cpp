#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

using boost::multiprecision::cpp_int;

namespace {

using Polynomial = std::vector<cpp_int>;

struct Failure {
    bool found = false;
    int minus_label = 0;
    std::vector<int> plus_labels;
    int index = -1;
    cpp_int left = 0;
    cpp_int right = 0;
};

long long binomial(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    long long result = 1;
    for (int j = 1; j <= k; ++j) {
        result = result * (n - k + j) / j;
    }
    return result;
}

Polynomial one_minus_polynomial(
    int minus_label,
    const std::vector<int>& plus_labels
) {
    int plus_sum = 0;
    for (int label : plus_labels) {
        plus_sum += label;
    }
    const int total = minus_label + plus_sum;
    std::vector<Polynomial> state(
        static_cast<std::size_t>(plus_sum + 1),
        Polynomial(static_cast<std::size_t>(2 * plus_sum + 1))
    );
    state[0U][0U] = 1;
    int support = 0;
    int degree = 0;
    for (int label : plus_labels) {
        std::vector<Polynomial> next(
            static_cast<std::size_t>(plus_sum + 1),
            Polynomial(static_cast<std::size_t>(2 * plus_sum + 1))
        );
        for (int character = 0; character <= support; ++character) {
            for (int exponent = 0; exponent <= degree; ++exponent) {
                const cpp_int& coefficient = state[static_cast<std::size_t>(
                    character
                )][static_cast<std::size_t>(exponent)];
                if (coefficient == 0) {
                    continue;
                }
                for (int j = 0; j <= label; ++j) {
                    next[static_cast<std::size_t>(character)]
                        [static_cast<std::size_t>(exponent + 2 * j)]
                        += coefficient;
                }
                for (int output = std::abs(character - label);
                     output <= character + label; output += 2) {
                    next[static_cast<std::size_t>(output)]
                        [static_cast<std::size_t>(exponent + label)]
                        += coefficient;
                }
            }
        }
        support += label;
        degree += 2 * label;
        state = std::move(next);
    }

    Polynomial result(static_cast<std::size_t>(total + 1));
    const Polynomial& invariant = state[0U];
    for (int exponent = 0; exponent <= 2 * plus_sum; exponent += 2) {
        const cpp_int& coefficient
            = invariant[static_cast<std::size_t>(exponent)];
        for (int j = 0; j <= minus_label; ++j) {
            result[static_cast<std::size_t>(exponent / 2 + j)]
                += coefficient;
        }
    }
    if (minus_label <= plus_sum) {
        const Polynomial& target
            = state[static_cast<std::size_t>(minus_label)];
        for (int exponent = 0; exponent <= 2 * plus_sum; exponent += 2) {
            const int shifted = exponent + minus_label;
            if ((shifted & 1) == 0) {
                result[static_cast<std::size_t>(shifted / 2)]
                    -= target[static_cast<std::size_t>(exponent)];
            }
        }
    }
    return result;
}

std::vector<std::vector<int>> plus_multisets(
    int maximum_label,
    int maximum_factors
) {
    std::vector<std::vector<int>> result;
    std::vector<int> current;
    const auto visit = [&](const auto& self, int first) -> void {
        result.push_back(current);
        if (static_cast<int>(current.size()) == maximum_factors) {
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

void record_failure(
    Failure& failure,
    int minus_label,
    const std::vector<int>& plus_labels,
    int index,
    const cpp_int& left,
    const cpp_int& right
) {
    if (failure.found) {
        return;
    }
    failure.found = true;
    failure.minus_label = minus_label;
    failure.plus_labels = plus_labels;
    failure.index = index;
    failure.left = left;
    failure.right = right;
}

void print_word(const Failure& failure) {
    std::cout << "minus=" << failure.minus_label << " plus=[";
    for (std::size_t index = 0; index < failure.plus_labels.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << failure.plus_labels[index];
    }
    std::cout << "] index=" << failure.index
              << " left=" << failure.left << " right=" << failure.right;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::runtime_error(
                "usage: analyze_su2_one_minus_polynomial MAXIMUM_LABEL "
                "MAXIMUM_PLUS_FACTORS"
            );
        }
        const int maximum_label = std::stoi(argv[1]);
        const int maximum_factors = std::stoi(argv[2]);
        if (maximum_label < 1 || maximum_factors < 0) {
            throw std::runtime_error("invalid bound");
        }
        const auto words = plus_multisets(maximum_label, maximum_factors);
        Failure coefficient_failure;
        Failure palindromy_failure;
        Failure unimodality_failure;
        Failure log_concavity_failure;
        Failure gamma_failure;
        std::uint64_t tested = 0;
        std::mutex failure_mutex;

#pragma omp parallel for schedule(dynamic) reduction(+:tested)
        for (std::size_t word_index = 0; word_index < words.size();
             ++word_index) {
            for (int minus_label = 1; minus_label <= maximum_label;
                 ++minus_label) {
                ++tested;
                const Polynomial polynomial = one_minus_polynomial(
                    minus_label, words[word_index]
                );
                const int degree = static_cast<int>(polynomial.size()) - 1;
                Failure local_coefficient;
                Failure local_palindromy;
                Failure local_unimodality;
                Failure local_log_concavity;
                Failure local_gamma;
                for (int index = 0; index <= degree; ++index) {
                    if (polynomial[static_cast<std::size_t>(index)] < 0) {
                        record_failure(
                            local_coefficient, minus_label, words[word_index],
                            index, polynomial[static_cast<std::size_t>(index)],
                            0
                        );
                        break;
                    }
                    if (polynomial[static_cast<std::size_t>(index)]
                        != polynomial[static_cast<std::size_t>(degree-index)]) {
                        record_failure(
                            local_palindromy, minus_label, words[word_index],
                            index, polynomial[static_cast<std::size_t>(index)],
                            polynomial[static_cast<std::size_t>(degree-index)]
                        );
                        break;
                    }
                }
                for (int index = 1; index <= degree / 2; ++index) {
                    if (polynomial[static_cast<std::size_t>(index)]
                        < polynomial[static_cast<std::size_t>(index-1)]) {
                        record_failure(
                            local_unimodality, minus_label, words[word_index],
                            index, polynomial[static_cast<std::size_t>(index)],
                            polynomial[static_cast<std::size_t>(index-1)]
                        );
                        break;
                    }
                }
                for (int index = 1; index < degree; ++index) {
                    const cpp_int square
                        = polynomial[static_cast<std::size_t>(index)]
                        * polynomial[static_cast<std::size_t>(index)];
                    const cpp_int adjacent
                        = polynomial[static_cast<std::size_t>(index-1)]
                        * polynomial[static_cast<std::size_t>(index+1)];
                    if (square < adjacent) {
                        record_failure(
                            local_log_concavity, minus_label,
                            words[word_index], index, square, adjacent
                        );
                        break;
                    }
                }
                Polynomial residual = polynomial;
                for (int index = 0; index <= degree / 2; ++index) {
                    const cpp_int gamma
                        = residual[static_cast<std::size_t>(index)];
                    if (gamma < 0) {
                        record_failure(
                            local_gamma, minus_label, words[word_index],
                            index, gamma, 0
                        );
                        break;
                    }
                    const int power = degree - 2 * index;
                    for (int j = 0; j <= power; ++j) {
                        residual[static_cast<std::size_t>(index+j)]
                            -= gamma * binomial(power, j);
                    }
                }
                if (local_coefficient.found || local_palindromy.found
                    || local_unimodality.found
                    || local_log_concavity.found || local_gamma.found) {
                    std::lock_guard<std::mutex> lock(failure_mutex);
                    if (local_coefficient.found && !coefficient_failure.found) {
                        coefficient_failure = local_coefficient;
                    }
                    if (local_palindromy.found && !palindromy_failure.found) {
                        palindromy_failure = local_palindromy;
                    }
                    if (local_unimodality.found && !unimodality_failure.found) {
                        unimodality_failure = local_unimodality;
                    }
                    if (local_log_concavity.found
                        && !log_concavity_failure.found) {
                        log_concavity_failure = local_log_concavity;
                    }
                    if (local_gamma.found && !gamma_failure.found) {
                        gamma_failure = local_gamma;
                    }
                }
            }
        }

        std::cout << "SU2_ONE_MINUS_POLYNOMIAL tested=" << tested
                  << " maximum_label=" << maximum_label
                  << " maximum_plus_factors=" << maximum_factors << '\n';
        const auto report = [](const std::string& property,
                               const Failure& failure) {
            std::cout << property << '=';
            if (!failure.found) {
                std::cout << "PASS\n";
                return;
            }
            std::cout << "FAIL ";
            print_word(failure);
            std::cout << '\n';
        };
        report("coefficient_nonnegative", coefficient_failure);
        report("palindromic", palindromy_failure);
        report("initial_unimodal", unimodality_failure);
        report("log_concave", log_concavity_failure);
        report("gamma_nonnegative", gamma_failure);
        return coefficient_failure.found || palindromy_failure.found
                || unimodality_failure.found
            ? EXIT_FAILURE : EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
