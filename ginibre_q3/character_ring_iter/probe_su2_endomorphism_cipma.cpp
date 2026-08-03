#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

using boost::multiprecision::cpp_int;

namespace {

enum class Mode { ordinary, finite };

enum class ProfileRestriction {
    none,
    parity_interval,
    parity_log_concave,
    parity_log_concave_real_rooted
};

using Rational = boost::rational<cpp_int>;
using Polynomial = std::vector<Rational>;

[[nodiscard]] int fusion_upper(const int left, const int right, const int level,
                               const Mode mode) {
    if (mode == Mode::ordinary) {
        return left + right;
    }
    return std::min(left + right, 2 * level - left - right);
}

template <class Function>
void for_each_fusion_output(const int left, const int right, const int level,
                            const Mode mode, Function&& function) {
    const int lower = std::abs(left - right);
    const int upper = fusion_upper(left, right, level, mode);
    for (int output = lower; output <= upper; output += 2) {
        function(output);
    }
}

[[nodiscard]] std::vector<cpp_int> endomorphism_character(
    const std::vector<unsigned int>& multiplicity, const int level, const Mode mode
) {
    const int largest_label = static_cast<int>(multiplicity.size()) - 1;
    const int output_bound = mode == Mode::ordinary ? 2 * largest_label : level;
    std::vector<cpp_int> character(static_cast<std::size_t>(output_bound + 1));
    for (int left = 0; left <= largest_label; ++left) {
        if (multiplicity[static_cast<std::size_t>(left)] == 0U) {
            continue;
        }
        for (int right = 0; right <= largest_label; ++right) {
            if (multiplicity[static_cast<std::size_t>(right)] == 0U) {
                continue;
            }
            const cpp_int weight = cpp_int(multiplicity[static_cast<std::size_t>(left)])
                * cpp_int(multiplicity[static_cast<std::size_t>(right)]);
            for_each_fusion_output(left, right, level, mode, [&](const int output) {
                character[static_cast<std::size_t>(output)] += weight;
            });
        }
    }
    return character;
}

[[nodiscard]] cpp_int fusion_by_character_entry(
    const std::vector<cpp_int>& character, const int left, const int right,
    const int level, const Mode mode
) {
    cpp_int result = 0;
    for_each_fusion_output(left, right, level, mode, [&](const int output) {
        if (output < static_cast<int>(character.size())) {
            result += character[static_cast<std::size_t>(output)];
        }
    });
    return result;
}

void print_vector(const std::vector<unsigned int>& value) {
    std::cout << '[';
    for (std::size_t index = 0; index < value.size(); ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << value[index];
    }
    std::cout << ']';
}

[[nodiscard]] bool has_cipma_counterexample(
    const std::vector<unsigned int>& multiplicity, const int level, const Mode mode,
    int& counter_left, int& counter_right, cpp_int& counter_value
) {
    const std::vector<cpp_int> character
        = endomorphism_character(multiplicity, level, mode);
    const cpp_int dimension_zero = character.front();
    const int test_bound = static_cast<int>(character.size()) - 1;
    for (int left = 0; left <= test_bound; ++left) {
        for (int right = 0; right <= test_bound; ++right) {
            const cpp_int value = dimension_zero
                    * fusion_by_character_entry(character, left, right, level, mode)
                - character[static_cast<std::size_t>(left)]
                    * character[static_cast<std::size_t>(right)];
            if (value < 0) {
                counter_left = left;
                counter_right = right;
                counter_value = value;
                return true;
            }
        }
    }
    return false;
}

[[nodiscard]] Mode parse_mode(const std::string& text) {
    if (text == "--ordinary") {
        return Mode::ordinary;
    }
    if (text == "--finite") {
        return Mode::finite;
    }
    throw std::runtime_error("mode must be --ordinary or --finite");
}

[[nodiscard]] ProfileRestriction parse_restriction(const std::string& text) {
    if (text == "--all") {
        return ProfileRestriction::none;
    }
    if (text == "--parity-interval") {
        return ProfileRestriction::parity_interval;
    }
    if (text == "--parity-logconcave") {
        return ProfileRestriction::parity_log_concave;
    }
    if (text == "--parity-logconcave-realrooted") {
        return ProfileRestriction::parity_log_concave_real_rooted;
    }
    throw std::runtime_error(
        "restriction must be --all, --parity-interval, --parity-logconcave, "
        "or --parity-logconcave-realrooted"
    );
}

void trim(Polynomial& polynomial) {
    while (!polynomial.empty() && polynomial.back() == Rational(0)) {
        polynomial.pop_back();
    }
}

[[nodiscard]] int degree(const Polynomial& polynomial) {
    return static_cast<int>(polynomial.size()) - 1;
}

[[nodiscard]] Polynomial add(const Polynomial& left, const Polynomial& right,
                             const Rational right_scale = Rational(1)) {
    Polynomial result(std::max(left.size(), right.size()), Rational(0));
    for (std::size_t index = 0; index < left.size(); ++index) {
        result[index] += left[index];
    }
    for (std::size_t index = 0; index < right.size(); ++index) {
        result[index] += right_scale * right[index];
    }
    trim(result);
    return result;
}

[[nodiscard]] Polynomial multiply_x_minus_one(const Polynomial& polynomial) {
    Polynomial result(polynomial.size() + 1U, Rational(0));
    for (std::size_t index = 0; index < polynomial.size(); ++index) {
        result[index] -= polynomial[index];
        result[index + 1U] += polynomial[index];
    }
    trim(result);
    return result;
}

[[nodiscard]] Polynomial multiply_by_linear(const Polynomial& polynomial,
                                            const Rational constant,
                                            const Rational linear) {
    Polynomial result(polynomial.size() + 1U, Rational(0));
    for (std::size_t index = 0; index < polynomial.size(); ++index) {
        result[index] += constant * polynomial[index];
        result[index + 1U] += linear * polynomial[index];
    }
    trim(result);
    return result;
}

[[nodiscard]] std::vector<Polynomial> character_basis(const int maximum_index,
                                                       const bool odd_parity) {
    std::vector<Polynomial> result(
        static_cast<std::size_t>(maximum_index + 1), Polynomial{});
    result[0] = Polynomial{Rational(1)};
    if (maximum_index == 0) {
        return result;
    }
    result[1] = odd_parity ? Polynomial{Rational(-1), Rational(1)}
                           : Polynomial{Rational(0), Rational(1)};
    for (int index = 1; index < maximum_index; ++index) {
        result[static_cast<std::size_t>(index + 1)]
            = add(multiply_x_minus_one(result[static_cast<std::size_t>(index)]),
                  result[static_cast<std::size_t>(index - 1)], Rational(-1));
    }
    return result;
}

[[nodiscard]] Polynomial factor_profile_polynomial(
    const std::vector<unsigned int>& multiplicity, const int parity
) {
    const int maximum_index = (static_cast<int>(multiplicity.size()) - 1 - parity) / 2;
    const std::vector<Polynomial> basis = character_basis(maximum_index, parity != 0);
    Polynomial result;
    for (int index = 0; index <= maximum_index; ++index) {
        const unsigned int coefficient
            = multiplicity[static_cast<std::size_t>(parity + 2 * index)];
        if (coefficient == 0U) {
            continue;
        }
        Polynomial term = basis[static_cast<std::size_t>(index)];
        for (Rational& value : term) {
            value *= Rational(cpp_int(coefficient));
        }
        result = add(result, term);
    }
    return result;
}

[[nodiscard]] Rational evaluate(const Polynomial& polynomial, const int argument) {
    Rational result(0);
    for (auto iterator = polynomial.rbegin(); iterator != polynomial.rend(); ++iterator) {
        result = result * Rational(argument) + *iterator;
    }
    return result;
}

[[nodiscard]] Polynomial derivative(const Polynomial& polynomial) {
    if (polynomial.size() <= 1U) {
        return {};
    }
    Polynomial result(polynomial.size() - 1U, Rational(0));
    for (std::size_t index = 1; index < polynomial.size(); ++index) {
        result[index - 1U] = Rational(cpp_int(index)) * polynomial[index];
    }
    trim(result);
    return result;
}

[[nodiscard]] std::pair<Polynomial, Polynomial> divide_with_remainder(
    Polynomial numerator, const Polynomial& denominator
) {
    if (denominator.empty()) {
        throw std::runtime_error("division by the zero polynomial");
    }
    Polynomial quotient(static_cast<std::size_t>(std::max(
        0, degree(numerator) - degree(denominator) + 1)), Rational(0));
    while (!numerator.empty() && degree(numerator) >= degree(denominator)) {
        const int shift = degree(numerator) - degree(denominator);
        const Rational factor = numerator.back() / denominator.back();
        quotient[static_cast<std::size_t>(shift)] += factor;
        for (int index = 0; index <= degree(denominator); ++index) {
            numerator[static_cast<std::size_t>(index + shift)]
                -= factor * denominator[static_cast<std::size_t>(index)];
        }
        trim(numerator);
    }
    trim(quotient);
    return {std::move(quotient), std::move(numerator)};
}

[[nodiscard]] Polynomial monic(Polynomial polynomial) {
    if (polynomial.empty()) {
        return polynomial;
    }
    const Rational leading = polynomial.back();
    for (Rational& value : polynomial) {
        value /= leading;
    }
    return polynomial;
}

[[nodiscard]] Polynomial polynomial_gcd(Polynomial left, Polynomial right) {
    while (!right.empty()) {
        Polynomial remainder = divide_with_remainder(left, right).second;
        left = std::move(right);
        right = std::move(remainder);
    }
    return monic(std::move(left));
}

[[nodiscard]] int sign(const Rational& value) {
    return value > Rational(0) ? 1 : (value < Rational(0) ? -1 : 0);
}

[[nodiscard]] int sturm_variations(const std::vector<Polynomial>& sequence,
                                    const int argument) {
    int previous = 0;
    int variations = 0;
    for (const Polynomial& polynomial : sequence) {
        const int current = sign(evaluate(polynomial, argument));
        if (current == 0) {
            continue;
        }
        if (previous != 0 && previous != current) {
            ++variations;
        }
        previous = current;
    }
    return variations;
}

[[nodiscard]] bool has_interval_roots(Polynomial polynomial) {
    if (polynomial.empty()) {
        return false;
    }
    for (const int endpoint : {-1, 3}) {
        const Polynomial divisor{Rational(-endpoint), Rational(1)};
        while (evaluate(polynomial, endpoint) == Rational(0)) {
            const auto [quotient, remainder] = divide_with_remainder(polynomial, divisor);
            if (!remainder.empty()) {
                throw std::runtime_error("nonexact endpoint-root division");
            }
            polynomial = quotient;
        }
    }
    if (degree(polynomial) <= 0) {
        return true;
    }
    const Polynomial gcd = polynomial_gcd(polynomial, derivative(polynomial));
    const auto [square_free, remainder] = divide_with_remainder(polynomial, gcd);
    if (!remainder.empty()) {
        throw std::runtime_error("nonexact square-free quotient");
    }
    std::vector<Polynomial> sturm{square_free, derivative(square_free)};
    while (!sturm.back().empty()) {
        Polynomial remainder_value
            = divide_with_remainder(sturm[sturm.size() - 2U], sturm.back()).second;
        for (Rational& value : remainder_value) {
            value = -value;
        }
        trim(remainder_value);
        if (remainder_value.empty()) {
            break;
        }
        sturm.push_back(std::move(remainder_value));
    }
    const int roots = sturm_variations(sturm, -1) - sturm_variations(sturm, 3);
    return roots == degree(square_free);
}

int replay_finite_counterexample_lift_scan() {
    constexpr int level = 22;
    constexpr int parity = 1;
    std::vector<unsigned int> multiplicity(static_cast<std::size_t>(level + 1));
    multiplicity[5] = 1U;
    multiplicity[7] = 2U;
    multiplicity[9] = 4U;
    multiplicity[11] = 8U;
    multiplicity[13] = 2U;
    const Polynomial canonical = factor_profile_polynomial(multiplicity, parity);
    if (has_interval_roots(canonical)) {
        throw std::runtime_error("counterexample canonical representative became real rooted");
    }
    const Polynomial modulus
        = character_basis(level / 2, true).back();
    if (!has_interval_roots(modulus)) {
        throw std::runtime_error("finite odd-character quotient lost interval roots");
    }
    for (int linear = -16; linear <= 16; ++linear) {
        for (int constant = -16; constant <= 16; ++constant) {
            if (linear == 0 && constant == 0) {
                continue;
            }
            const Polynomial correction = multiply_by_linear(
                modulus, Rational(constant), Rational(linear));
            const Polynomial lift = add(canonical, correction);
            if (has_interval_roots(lift)) {
                std::cout
                    << "SU2_ENDOMORPHISM_CIPMA_LIFT result=COUNTEREXAMPLE"
                    << " level=" << level
                    << " parity=" << parity
                    << " correction_linear=" << linear
                    << " correction_constant=" << constant
                    << '\n';
                return EXIT_FAILURE;
            }
        }
    }
    std::cout
        << "SU2_ENDOMORPHISM_CIPMA_LIFT result=NO_SMALL_LIFT"
        << " level=" << level
        << " parity=" << parity
        << " correction_degree<=1"
        << " coefficient_box=[-16,16]"
        << '\n';
    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--replay-finite-counterexample-lift") {
            return replay_finite_counterexample_lift_scan();
        }
        if (argc != 5) {
            throw std::runtime_error(
                "usage: probe_su2_endomorphism_cipma "
                "(--ordinary|--finite) "
                "(--all|--parity-interval|--parity-logconcave|"
                "--parity-logconcave-realrooted) "
                "MAXIMUM_LABEL MAXIMUM_MULTIPLICITY"
            );
        }
        const Mode mode = parse_mode(argv[1]);
        const ProfileRestriction restriction = parse_restriction(argv[2]);
        const int maximum_label = std::stoi(argv[3]);
        const unsigned long parsed_multiplicity = std::stoul(argv[4]);
        if (maximum_label < 0 || parsed_multiplicity == 0UL
            || parsed_multiplicity > static_cast<unsigned long>(UINT_MAX)) {
            throw std::runtime_error("invalid search bound");
        }
        const unsigned int maximum_multiplicity
            = static_cast<unsigned int>(parsed_multiplicity);
        const int level = maximum_label;
        const int vector_size = maximum_label + 1;
        std::vector<unsigned int> multiplicity(static_cast<std::size_t>(vector_size));
        cpp_int tested = 0;
        bool found = false;
        int counter_left = -1;
        int counter_right = -1;
        cpp_int counter_value = 0;
        std::vector<unsigned int> counterexample;

        const auto test_current_profile = [&]() {
            if (restriction == ProfileRestriction::parity_log_concave_real_rooted) {
                int parity = -1;
                for (int label = 0; label <= maximum_label; ++label) {
                    if (multiplicity[static_cast<std::size_t>(label)] != 0U) {
                        parity = label % 2;
                        break;
                    }
                }
                if (parity < 0
                    || !has_interval_roots(factor_profile_polynomial(multiplicity, parity))) {
                    return;
                }
            }
            ++tested;
            if (has_cipma_counterexample(multiplicity, level, mode, counter_left,
                                          counter_right, counter_value)) {
                found = true;
                counterexample = multiplicity;
            }
        };

        if (restriction == ProfileRestriction::none) {
            const auto enumerate_all = [&](const auto& self, const int index,
                                           const bool nonzero) -> void {
                if (found) {
                    return;
                }
                if (index == vector_size) {
                    if (nonzero) {
                        test_current_profile();
                    }
                    return;
                }
                for (unsigned int value = 0; value <= maximum_multiplicity; ++value) {
                    multiplicity[static_cast<std::size_t>(index)] = value;
                    self(self, index + 1, nonzero || value != 0U);
                    if (found) {
                        return;
                    }
                }
            };
            enumerate_all(enumerate_all, 0, false);
        } else {
            for (int parity = 0; parity <= 1 && !found; ++parity) {
                for (int first = parity; first <= maximum_label && !found; first += 2) {
                    for (int last = first; last <= maximum_label && !found; last += 2) {
                        std::fill(multiplicity.begin(), multiplicity.end(), 0U);
                        const auto enumerate_interval = [&] (const auto& self,
                                                             const int position) -> void {
                            if (found) {
                                return;
                            }
                            if (position > last) {
                                test_current_profile();
                                return;
                            }
                            for (unsigned int value = 1U;
                                 value <= maximum_multiplicity; ++value) {
                                multiplicity[static_cast<std::size_t>(position)] = value;
                                const bool violates_log_concavity
                                    = restriction
                                            == ProfileRestriction::parity_log_concave
                                        || restriction
                                            == ProfileRestriction::parity_log_concave_real_rooted;
                                const bool fails_local_log_concavity
                                    = violates_log_concavity
                                        && position >= first + 4
                                        && cpp_int(multiplicity[static_cast<std::size_t>(
                                                       position - 2)])
                                                * cpp_int(multiplicity[static_cast<std::size_t>(
                                                    position - 2)])
                                            < cpp_int(multiplicity[static_cast<std::size_t>(
                                                       position - 4)])
                                                * cpp_int(value);
                                if (!fails_local_log_concavity) {
                                    self(self, position + 2);
                                }
                                if (found) {
                                    return;
                                }
                            }
                        };
                        enumerate_interval(enumerate_interval, first);
                    }
                }
            }
        }

        std::cout << "SU2_ENDOMORPHISM_CIPMA mode="
                  << (mode == Mode::ordinary ? "ordinary" : "finite")
                  << " restriction="
                  << (restriction == ProfileRestriction::none
                          ? "all"
                          : restriction == ProfileRestriction::parity_interval
                              ? "parity_interval"
                              : restriction == ProfileRestriction::parity_log_concave
                                  ? "parity_logconcave"
                                  : "parity_logconcave_realrooted")
                  << " maximum_label=" << maximum_label
                  << " maximum_multiplicity=" << maximum_multiplicity
                  << " tested=" << tested;
        if (found) {
            std::cout << " result=COUNTEREXAMPLE multiplicity=";
            print_vector(counterexample);
            std::cout << " left=" << counter_left << " right=" << counter_right
                      << " value=" << counter_value << '\n';
            return EXIT_FAILURE;
        }
        std::cout << " result=PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_ENDOMORPHISM_CIPMA error=" << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
