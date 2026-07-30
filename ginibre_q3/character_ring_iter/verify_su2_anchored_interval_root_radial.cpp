#include <array>
#include <cstddef>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;
using Exponent = std::array<int, 3>;

class Polynomial {
public:
    Polynomial() = default;

    explicit Polynomial(const Integer& constant) {
        if (constant != 0) {
            terms_[Exponent{0, 0, 0}] = constant;
        }
    }

    static Polynomial variable(const std::size_t index) {
        if (index >= 3U) {
            throw std::runtime_error("polynomial variable index out of range");
        }
        Polynomial result;
        Exponent exponent{0, 0, 0};
        exponent[index] = 1;
        result.terms_[exponent] = 1;
        return result;
    }

    Polynomial& operator+=(const Polynomial& other) {
        for (const auto& [exponent, coefficient] : other.terms_) {
            terms_[exponent] += coefficient;
            if (terms_[exponent] == 0) {
                terms_.erase(exponent);
            }
        }
        return *this;
    }

    Polynomial& operator-=(const Polynomial& other) {
        for (const auto& [exponent, coefficient] : other.terms_) {
            terms_[exponent] -= coefficient;
            if (terms_[exponent] == 0) {
                terms_.erase(exponent);
            }
        }
        return *this;
    }

    Polynomial& operator*=(const Polynomial& other) {
        std::map<Exponent, Integer> product;
        for (const auto& [left_exponent, left_coefficient] : terms_) {
            for (const auto& [right_exponent, right_coefficient]
                 : other.terms_) {
                Exponent exponent{};
                for (std::size_t index = 0; index < exponent.size();
                     ++index) {
                    exponent[index]
                        = left_exponent[index] + right_exponent[index];
                }
                product[exponent] += left_coefficient * right_coefficient;
            }
        }
        terms_ = std::move(product);
        return *this;
    }

    [[nodiscard]] bool nonnegative_coefficients() const {
        for (const auto& [exponent, coefficient] : terms_) {
            static_cast<void>(exponent);
            if (coefficient < 0) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] std::size_t term_count() const {
        return terms_.size();
    }

    [[nodiscard]] Integer minimum_coefficient() const {
        if (terms_.empty()) {
            return 0;
        }
        auto iterator = terms_.begin();
        Integer minimum = iterator->second;
        ++iterator;
        for (; iterator != terms_.end(); ++iterator) {
            if (iterator->second < minimum) {
                minimum = iterator->second;
            }
        }
        return minimum;
    }

    [[nodiscard]] const std::map<Exponent, Integer>& terms() const {
        return terms_;
    }

    [[nodiscard]] Integer evaluate(
        const std::array<Integer, 3>& values) const {
        Integer result = 0;
        for (const auto& [exponent, coefficient] : terms_) {
            Integer monomial = coefficient;
            for (std::size_t index = 0; index < exponent.size(); ++index) {
                for (int power = 0; power < exponent[index]; ++power) {
                    monomial *= values[index];
                }
            }
            result += monomial;
        }
        return result;
    }

private:
    std::map<Exponent, Integer> terms_;
};

Polynomial operator+(Polynomial left, const Polynomial& right) {
    left += right;
    return left;
}

Polynomial operator-(Polynomial left, const Polynomial& right) {
    left -= right;
    return left;
}

Polynomial operator*(Polynomial left, const Polynomial& right) {
    left *= right;
    return left;
}

Polynomial operator*(const Integer& scalar, Polynomial value) {
    value *= Polynomial(scalar);
    return value;
}

Polynomial low_twice(
    const Polynomial& support,
    const Polynomial& label) {
    return Integer(2) * (support + Polynomial(1))
               * (Integer(2) * label + Polynomial(1))
           - Integer(3) * label * (label + Polynomial(1));
}

Polynomial high_twice(
    const Polynomial& support,
    const Polynomial& label) {
    const Polynomial distance = Integer(2) * support - label;
    return (distance + Polynomial(1)) * (distance + Polynomial(2));
}

Polynomial separated_low_twice(
    const Polynomial& support,
    const Polynomial& label) {
    return Integer(2) * (support + Polynomial(1))
               * (Integer(2) * label + Polynomial(1))
           - Integer(2) * label * (label + Polynomial(1));
}

Polynomial separated_high_twice(const Polynomial& support) {
    return Integer(2) * (support + Polynomial(1))
           * (support + Polynomial(1));
}

struct Chamber {
    std::string name;
    Polynomial support;
    Polynomial antidiagonal;
    Polynomial contraction;
    bool middle_high = false;
    bool middle_next_high = false;
    bool top_one_high = false;
    bool top_two_high = false;
};

Polynomial twice_coefficient(
    const Polynomial& support,
    const Polynomial& label,
    const bool high) {
    return high ? high_twice(support, label)
                : low_twice(support, label);
}

Polynomial scaled_radial(const Chamber& chamber) {
    const Polynomial middle
        = chamber.antidiagonal - chamber.contraction;
    const Polynomial c_zero
        = low_twice(chamber.support, Polynomial(0));
    const Polynomial c_l = low_twice(
        chamber.support, chamber.contraction);
    const Polynomial c_l_next = low_twice(
        chamber.support, chamber.contraction + Polynomial(1));
    const Polynomial c_middle = twice_coefficient(
        chamber.support, middle, chamber.middle_high);
    const Polynomial c_middle_next = twice_coefficient(
        chamber.support,
        middle + Polynomial(1),
        chamber.middle_next_high);
    const Polynomial c_top_one = twice_coefficient(
        chamber.support,
        chamber.antidiagonal + Polynomial(1),
        chamber.top_one_high);
    const Polynomial c_top_two = twice_coefficient(
        chamber.support,
        chamber.antidiagonal + Polynomial(2),
        chamber.top_two_high);
    return c_zero * (c_top_one + c_top_two)
           + c_l * c_middle - c_l_next * c_middle_next;
}

Polynomial scaled_separated_radial(const Chamber& chamber) {
    const Polynomial middle
        = chamber.antidiagonal - chamber.contraction;
    const auto coefficient = [&chamber](
                                 const Polynomial& label,
                                 const bool high) {
        return high
            ? separated_high_twice(chamber.support)
            : separated_low_twice(chamber.support, label);
    };
    const Polynomial c_zero
        = separated_low_twice(chamber.support, Polynomial(0));
    const Polynomial c_l = separated_low_twice(
        chamber.support, chamber.contraction);
    const Polynomial c_l_next = separated_low_twice(
        chamber.support, chamber.contraction + Polynomial(1));
    const Polynomial c_middle = coefficient(
        middle, chamber.middle_high);
    const Polynomial c_middle_next = coefficient(
        middle + Polynomial(1), chamber.middle_next_high);
    const Polynomial c_top_one = coefficient(
        chamber.antidiagonal + Polynomial(1),
        chamber.top_one_high);
    const Polynomial c_top_two = coefficient(
        chamber.antidiagonal + Polynomial(2),
        chamber.top_two_high);
    return c_zero * (c_top_one + c_top_two)
           + c_l * c_middle - c_l_next * c_middle_next;
}

void print_negative_terms(const Polynomial& polynomial) {
    for (const auto& [exponent, coefficient] : polynomial.terms()) {
        if (coefficient < 0) {
            std::cout << " negative=(" << exponent[0] << ','
                      << exponent[1] << ',' << exponent[2]
                      << "):" << coefficient;
        }
    }
}

Integer coefficient(const int support, const int label) {
    if (label < 0 || label > 2 * support) {
        return 0;
    }
    if (label <= support) {
        return (support + 1) * (2 * label + 1)
               - 3 * label * (label + 1) / 2;
    }
    const int distance = 2 * support - label;
    return (distance + 1) * (distance + 2) / 2;
}

Integer direct_coefficient(const int support, const int label) {
    Integer result = 0;
    for (int left = 0; left <= support; ++left) {
        for (int right = 0; right <= support; ++right) {
            if (std::abs(left - right) <= label
                && label <= left + right) {
                ++result;
            }
        }
    }
    return result;
}

Integer separated_coefficient(const int support, const int label) {
    if (label < 0 || label > 2 * support) {
        return 0;
    }
    if (label <= support) {
        return (support + 1) * (2 * label + 1)
               - label * (label + 1);
    }
    return Integer(support + 1) * (support + 1);
}

enum class TransitionBranch {
    Plateau,
    FirstDecline,
    Tail,
};

struct TransitionChamber {
    std::string name;
    Polynomial support;
    Polynomial offset;
    Polynomial contraction;
    TransitionBranch middle = TransitionBranch::Plateau;
    TransitionBranch middle_next = TransitionBranch::Plateau;
    TransitionBranch top = TransitionBranch::FirstDecline;
    TransitionBranch top_next = TransitionBranch::FirstDecline;
};

Polynomial transition_twice(
    const Polynomial& support,
    const Polynomial& offset,
    const TransitionBranch branch) {
    const Polynomial size = support + Polynomial(1);
    if (branch == TransitionBranch::Plateau) {
        return Integer(2) * size * size;
    }
    if (branch == TransitionBranch::FirstDecline) {
        return Integer(2) * size * size
               - offset * (offset + Polynomial(1));
    }
    const Polynomial distance = Integer(2) * support - offset;
    return (distance + Polynomial(1)) * (distance + Polynomial(2));
}

Polynomial scaled_transition_radial(const TransitionChamber& chamber) {
    const Polynomial difference
        = chamber.offset - chamber.contraction;
    const Polynomial c_zero
        = Integer(2) * (chamber.support + Polynomial(1));
    const Polynomial c_l = separated_low_twice(
        chamber.support, chamber.contraction);
    const Polynomial c_l_next = separated_low_twice(
        chamber.support, chamber.contraction + Polynomial(1));
    const Polynomial c_middle = transition_twice(
        chamber.support,
        difference - Polynomial(1),
        chamber.middle);
    const Polynomial c_middle_next = transition_twice(
        chamber.support, difference, chamber.middle_next);
    const Polynomial c_top = transition_twice(
        chamber.support, chamber.offset, chamber.top);
    const Polynomial c_top_next = transition_twice(
        chamber.support,
        chamber.offset + Polynomial(1),
        chamber.top_next);
    return c_zero * (c_top + c_top_next)
           + c_l * c_middle - c_l_next * c_middle_next;
}

Integer transition_coefficient(const int support, const int offset) {
    const Integer size = support + 1;
    if (offset <= 0) {
        return size * size;
    }
    if (offset <= support) {
        return size * size - offset * (offset + 1) / 2;
    }
    if (offset <= 2 * support) {
        const int distance = 2 * support - offset;
        return Integer(distance + 1) * (distance + 2) / 2;
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const bool separated_mode
            = argc == 2
              && std::string(argv[1]) == "--separated-lower-band";
        const bool transition_mode
            = argc == 2
              && std::string(argv[1]) == "--separated-transition-band";
        if (argc > 2
            || (argc == 2 && !separated_mode && !transition_mode)) {
            throw std::runtime_error(
                "usage: verify_su2_anchored_interval_root_radial "
                "[--separated-lower-band|--separated-transition-band]");
        }
        const Polynomial x = Polynomial::variable(0);
        const Polynomial y = Polynomial::variable(1);
        const Polynomial z = Polynomial::variable(2);

        if (transition_mode) {
            const std::array<TransitionChamber, 7> transition_chambers{{
                {
                    "B_le_L",
                    x + y + z + Polynomial(1),
                    x,
                    x + y,
                    TransitionBranch::Plateau,
                    TransitionBranch::Plateau,
                    TransitionBranch::FirstDecline,
                    TransitionBranch::FirstDecline,
                },
                {
                    "B_eq_L_plus_1_below_n",
                    x + y + Polynomial(2),
                    x + Polynomial(1),
                    x,
                    TransitionBranch::Plateau,
                    TransitionBranch::FirstDecline,
                    TransitionBranch::FirstDecline,
                    TransitionBranch::FirstDecline,
                },
                {
                    "B_eq_L_plus_1_eq_n",
                    x + Polynomial(1),
                    x + Polynomial(1),
                    x,
                    TransitionBranch::Plateau,
                    TransitionBranch::FirstDecline,
                    TransitionBranch::Tail,
                    TransitionBranch::Tail,
                },
                {
                    "two_le_gap_B_below_n",
                    x + y + z + Polynomial(3),
                    x + y + Polynomial(2),
                    x,
                    TransitionBranch::FirstDecline,
                    TransitionBranch::FirstDecline,
                    TransitionBranch::FirstDecline,
                    TransitionBranch::FirstDecline,
                },
                {
                    "two_le_gap_B_above_n_positive_excess",
                    x + y + z + Polynomial(2),
                    Integer(2) * x + y + z + Polynomial(3),
                    x + y + Polynomial(1),
                    TransitionBranch::FirstDecline,
                    TransitionBranch::FirstDecline,
                    TransitionBranch::Tail,
                    TransitionBranch::Tail,
                },
                {
                    "two_le_gap_B_eq_n",
                    x + y + Polynomial(2),
                    x + y + Polynomial(2),
                    y,
                    TransitionBranch::FirstDecline,
                    TransitionBranch::FirstDecline,
                    TransitionBranch::Tail,
                    TransitionBranch::Tail,
                },
                {
                    "gap_above_n",
                    x + y + z + Polynomial(2),
                    Integer(2) * x + Integer(2) * y + z
                        + Polynomial(3),
                    x,
                    TransitionBranch::Tail,
                    TransitionBranch::Tail,
                    TransitionBranch::Tail,
                    TransitionBranch::Tail,
                },
            }};

            bool passed = true;
            std::size_t total_terms = 0U;
            std::array<Polynomial, 7> transition_polynomials{};
            for (std::size_t index = 0U;
                 index < transition_chambers.size();
                 ++index) {
                const Polynomial radial = scaled_transition_radial(
                    transition_chambers[index]);
                transition_polynomials[index] = radial;
                const bool nonnegative = radial.nonnegative_coefficients();
                passed = passed && nonnegative;
                total_terms += radial.term_count();
                std::cout
                    << "transition_chamber="
                    << transition_chambers[index].name
                    << " terms=" << radial.term_count()
                    << " minimum_coefficient="
                    << radial.minimum_coefficient()
                    << " result=" << (nonnegative ? "PASS" : "FAIL");
                if (!nonnegative) {
                    print_negative_terms(radial);
                }
                std::cout << '\n';
            }

            std::size_t coverage_checks = 0U;
            for (int support = 1; support <= 64; ++support) {
                for (int offset = 0;
                     offset <= 2 * support - 1;
                     ++offset) {
                    for (int contraction = 0;
                         contraction <= support - 1;
                         ++contraction) {
                        const int difference = offset - contraction;
                        std::size_t selected = 0U;
                        std::array<Integer, 3> slacks{0, 0, 0};
                        if (difference <= 0) {
                            selected = 0U;
                            slacks = {
                                offset,
                                contraction - offset,
                                support - contraction - 1};
                        } else if (
                            difference == 1
                            && contraction <= support - 2) {
                            selected = 1U;
                            slacks = {
                                contraction,
                                support - contraction - 2,
                                0};
                        } else if (difference == 1) {
                            selected = 2U;
                            slacks = {contraction, 0, 0};
                        } else if (offset <= support - 1) {
                            selected = 3U;
                            slacks = {
                                contraction,
                                difference - 2,
                                support - offset - 1};
                        } else if (difference <= support) {
                            if (offset > support) {
                                selected = 4U;
                                slacks = {
                                    offset - support - 1,
                                    support - difference,
                                    difference - offset + support - 1};
                            } else {
                                selected = 5U;
                                slacks = {
                                    difference - 2,
                                    support - difference,
                                    0};
                            }
                        } else {
                            selected = 6U;
                            slacks = {
                                contraction,
                                difference - support - 1,
                                2 * support - contraction
                                    - difference - 1};
                        }
                        for (const Integer& slack : slacks) {
                            if (slack < 0) {
                                throw std::runtime_error(
                                    "negative transition chamber slack");
                            }
                        }
                        const Integer c_zero = support + 1;
                        const Integer c_l = separated_coefficient(
                            support, contraction);
                        const Integer c_l_next = separated_coefficient(
                            support, contraction + 1);
                        const Integer c_middle = transition_coefficient(
                            support, difference - 1);
                        const Integer c_middle_next
                            = transition_coefficient(
                                support, difference);
                        const Integer c_top = transition_coefficient(
                            support, offset);
                        const Integer c_top_next
                            = transition_coefficient(
                                support, offset + 1);
                        const Integer exact = Integer(4) * (
                            c_zero * (c_top + c_top_next)
                            + c_l * c_middle
                            - c_l_next * c_middle_next);
                        if (transition_polynomials[selected].evaluate(
                                slacks)
                            != exact) {
                            throw std::runtime_error(
                                "transition chamber coverage mismatch");
                        }
                        ++coverage_checks;
                    }
                }
            }
            std::cout
                << "SU2_SEPARATED_INTERVAL_TRANSITION_RADIAL"
                << " chambers=" << transition_chambers.size()
                << " terms=" << total_terms
                << " chamber_coverage_checks=" << coverage_checks
                << " result=" << (passed ? "PASS_EXACT" : "FAIL")
                << '\n';
            return passed ? 0 : 1;
        }

        // x, y, z are nonnegative integer slacks.  The initial chambers
        // separate the locations where the explicit coefficient formula
        // for (1_[0,n])^2 changes from its rising quadratic to its
        // reflected triangular tail.
        const std::array<Chamber, 7> chambers{{
            {
                "A_plus_2_le_n",
                Integer(2) * x + y + z + Polynomial(3),
                Integer(2) * x + y + Polynomial(1),
                x,
                false,
                false,
                false,
                false,
            },
            {
                "A_plus_1_eq_n",
                Integer(2) * x + y + Polynomial(2),
                Integer(2) * x + y + Polynomial(1),
                x,
                false,
                false,
                false,
                true,
            },
            {
                "A_eq_n_L_zero",
                y + Polynomial(1),
                y + Polynomial(1),
                Polynomial(0),
                false,
                true,
                true,
                true,
            },
            {
                "A_eq_n_L_positive",
                Integer(2) * x + y + Polynomial(3),
                Integer(2) * x + y + Polynomial(3),
                x + Polynomial(1),
                false,
                false,
                true,
                true,
            },
            {
                "A_gt_n_middle_le_n",
                Integer(2) * x + y + z + Polynomial(4),
                Integer(2) * x + y + Integer(2) * z
                    + Polynomial(5),
                x + z + Polynomial(2),
                false,
                false,
                true,
                true,
            },
            {
                "A_gt_n_middle_eq_n",
                x + y + Polynomial(3),
                Integer(2) * x + y + Polynomial(4),
                x + Polynomial(1),
                false,
                true,
                true,
                true,
            },
            {
                "A_gt_n_middle_gt_n",
                x + y + z + Polynomial(3),
                Integer(2) * x + y + Integer(2) * z
                    + Polynomial(4),
                x,
                true,
                true,
                true,
                true,
            },
        }};

        bool passed = true;
        std::size_t total_terms = 0U;
        std::array<Polynomial, 7> radial_polynomials{};
        std::size_t chamber_index = 0U;
        for (const Chamber& chamber : chambers) {
            const Polynomial radial = separated_mode
                ? scaled_separated_radial(chamber)
                : scaled_radial(chamber);
            radial_polynomials[chamber_index] = radial;
            ++chamber_index;
            const bool nonnegative = radial.nonnegative_coefficients();
            passed = passed && nonnegative;
            total_terms += radial.term_count();
            std::cout
                      << (separated_mode ? "separated_" : "")
                      << "chamber=" << chamber.name
                      << " terms=" << radial.term_count()
                      << " minimum_coefficient="
                      << radial.minimum_coefficient()
                      << " result="
                      << (nonnegative ? "PASS" : "FAIL");
            if (!nonnegative) {
                print_negative_terms(radial);
            }
            std::cout << '\n';
        }

        std::size_t coefficient_checks = 0U;
        std::size_t coverage_checks = 0U;
        for (int support = 0; support <= 64; ++support) {
            for (int label = 0; label <= 2 * support; ++label) {
                ++coefficient_checks;
                if (coefficient(support, label)
                    != direct_coefficient(support, label)) {
                    throw std::runtime_error(
                        "coefficient formula mismatch");
                }
            }
            for (int antidiagonal = 1;
                 antidiagonal <= 2 * support - 2;
                 ++antidiagonal) {
                for (int contraction = 0;
                     2 * contraction < antidiagonal;
                     ++contraction) {
                    const int middle = antidiagonal - contraction;
                    std::size_t selected = 0U;
                    std::array<Integer, 3> slacks{0, 0, 0};
                    if (antidiagonal + 2 <= support) {
                        selected = 0U;
                        slacks = {
                            contraction,
                            antidiagonal - 2 * contraction - 1,
                            support - antidiagonal - 2};
                    } else if (antidiagonal + 1 == support) {
                        selected = 1U;
                        slacks = {
                            contraction,
                            antidiagonal - 2 * contraction - 1,
                            0};
                    } else if (
                        antidiagonal == support && contraction == 0) {
                        selected = 2U;
                        slacks = {0, support - 1, 0};
                    } else if (antidiagonal == support) {
                        selected = 3U;
                        slacks = {
                            contraction - 1,
                            support - 2 * contraction - 1,
                            0};
                    } else if (middle < support) {
                        const int difference
                            = antidiagonal - support;
                        selected = 4U;
                        slacks = {
                            contraction - difference - 1,
                            support - difference
                                - 2 * (
                                    contraction - difference - 1)
                                - 3,
                            difference - 1};
                    } else if (middle == support) {
                        selected = 5U;
                        slacks = {
                            contraction - 1,
                            support - contraction - 2,
                            0};
                    } else {
                        selected = 6U;
                        slacks = {
                            contraction,
                            support - contraction
                                - (middle - support - 1) - 3,
                            middle - support - 1};
                    }
                    for (const Integer& slack : slacks) {
                        if (slack < 0) {
                            throw std::runtime_error(
                                "negative chamber slack");
                        }
                    }
                    const auto current_coefficient =
                        [separated_mode](const int n, const int label) {
                            return separated_mode
                                ? separated_coefficient(n, label)
                                : coefficient(n, label);
                        };
                    const Integer exact = Integer(4) * (
                        current_coefficient(support, 0)
                            * (
                                current_coefficient(
                                    support, antidiagonal + 1)
                                + current_coefficient(
                                    support, antidiagonal + 2))
                        + current_coefficient(support, contraction)
                            * current_coefficient(support, middle)
                        - current_coefficient(
                            support, contraction + 1)
                            * current_coefficient(
                                support, middle + 1));
                    if (radial_polynomials[selected].evaluate(slacks)
                        != exact) {
                        throw std::runtime_error(
                            "radial chamber coverage mismatch");
                    }
                    ++coverage_checks;
                }
            }
        }
        std::cout
            << (separated_mode
                    ? "SU2_SEPARATED_INTERVAL_LOWER_RADIAL"
                    : "SU2_ANCHORED_INTERVAL_ROOT_RADIAL")
            << " chambers=" << chambers.size()
            << " terms=" << total_terms
            << " coefficient_formula_checks=" << coefficient_checks
            << " chamber_coverage_checks=" << coverage_checks
            << " result=" << (passed ? "PASS_EXACT" : "FAIL")
            << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& exception) {
        std::cerr << "error: " << exception.what() << '\n';
        return 2;
    }
}
