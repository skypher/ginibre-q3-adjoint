#include <gmpxx.h>

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

mpz_class binomial(unsigned long n, unsigned long k) {
    mpz_class answer;
    mpz_bin_uiui(answer.get_mpz_t(), n, k);
    return answer;
}

mpz_class factorial(unsigned long n) {
    mpz_class answer;
    mpz_fac_ui(answer.get_mpz_t(), n);
    return answer;
}

mpz_class power_ui(unsigned long base, unsigned long exponent) {
    mpz_class answer;
    mpz_ui_pow_ui(answer.get_mpz_t(), base, exponent);
    return answer;
}

mpz_class standard_tableau_dimension(
    const std::vector<unsigned>& partition,
    unsigned total
) {
    mpz_class hook_product = 1;
    for (unsigned row = 0; row < partition.size(); ++row) {
        for (unsigned column = 0; column < partition[row]; ++column) {
            unsigned hook = partition[row] - column;
            for (unsigned lower = row + 1; lower < partition.size(); ++lower) {
                if (partition[lower] > column) {
                    ++hook;
                }
            }
            hook_product *= hook;
        }
    }
    const mpz_class numerator = factorial(total);
    require(mpz_divisible_p(numerator.get_mpz_t(), hook_product.get_mpz_t()),
            "nonintegral hook-length dimension");
    return numerator / hook_product;
}

mpz_class trace_moment(unsigned degree, unsigned maximum_length) {
    mpz_class answer = 0;
    std::vector<unsigned> partition;
    std::function<void(unsigned, unsigned)> visit =
        [&](unsigned remaining, unsigned maximum_part) {
            if (remaining == 0) {
                const mpz_class dimension =
                    standard_tableau_dimension(partition, degree);
                answer += dimension * dimension;
                return;
            }
            if (partition.size() == maximum_length) {
                return;
            }
            const unsigned upper = std::min(remaining, maximum_part);
            for (unsigned part = upper; part >= 1; --part) {
                partition.push_back(part);
                visit(remaining - part, part);
                partition.pop_back();
                if (part == 1) {
                    break;
                }
            }
        };
    visit(degree, degree);
    return answer;
}

std::vector<mpz_class> su3_adjoint_moments(unsigned maximum_degree) {
    std::vector<mpz_class> trace(maximum_degree + 1);
    std::vector<mpz_class> adjoint(maximum_degree + 1);
    for (unsigned degree = 0; degree <= maximum_degree; ++degree) {
        trace[degree] = trace_moment(degree, 3);
    }
    for (unsigned degree = 0; degree <= maximum_degree; ++degree) {
        mpz_class value = 0;
        for (unsigned j = 0; j <= degree; ++j) {
            const mpz_class term = binomial(degree, j) * trace[j];
            if ((degree - j) & 1U) {
                value -= term;
            } else {
                value += term;
            }
        }
        adjoint[degree] = value;
    }
    return adjoint;
}

mpz_class full_q3_value(
    const std::vector<mpz_class>& moment,
    unsigned a,
    unsigned n
) {
    mpz_class answer = 0;
    for (unsigned j = 0; j <= 2 * a; ++j) {
        for (unsigned k = 0; k <= n; ++k) {
            const mpz_class term = binomial(2 * a, j) * binomial(n, k)
                * moment[2 * a - j + k] * moment[j + n - k];
            if (j & 1U) {
                answer -= term;
            } else {
                answer += term;
            }
        }
    }
    return answer;
}

}  // namespace

int main() {
    try {
        constexpr unsigned maximum_degree = 54;
        const std::vector<mpz_class> moment =
            su3_adjoint_moments(maximum_degree);

        const std::vector<unsigned long> expected_prefix =
            {1, 0, 1, 2, 8, 32, 145, 702, 3598, 19280};
        require(moment.size() == maximum_degree + 1,
                "SU(3) moment generator returned the wrong length");
        for (unsigned index = 0; index < expected_prefix.size(); ++index) {
            require(moment[index] == expected_prefix[index],
                    "SU(3) hook-length moment prefix mismatch at degree "
                    + std::to_string(index));
        }

        // Independent check against the proved SU(3) recurrence.
        for (unsigned long k = 1; k < maximum_degree; ++k) {
            const mpz_class left = (k + 3) * (k + 3) * moment[k + 1];
            const mpz_class right = k * (7 * k + 11) * moment[k]
                + 8 * k * (k + 1) * moment[k - 1];
            require(left == right,
                    "SU(3) recurrence mismatch at k=" + std::to_string(k));
        }

        const mpz_class expected_m54(
            "1073719774632588649045995816896752254587680"
        );
        require(moment[54] == expected_m54,
                "SU(3) degree-54 moment identity mismatch");

        unsigned checked = 0;
        mpz_class minimum;
        unsigned minimum_a = 0;
        unsigned minimum_n = 0;
        bool have_minimum = false;
        for (unsigned a = 2; 2 * a + 1 <= 27; ++a) {
            for (unsigned n = 1; 2 * a + n <= 27; n += 2) {
                const mpz_class value = full_q3_value(moment, a, n);
                require(value > 0,
                        "nonpositive SU(3) full-Q3 residual at a="
                        + std::to_string(a) + ", n=" + std::to_string(n)
                        + ", value=" + value.get_str());
                if (!have_minimum || value < minimum) {
                    minimum = value;
                    minimum_a = a;
                    minimum_n = n;
                    have_minimum = true;
                }
                std::cout << "SU3_FULL_Q3 a=" << a
                          << " n=" << n
                          << " total_degree=" << (2 * a + n)
                          << " value=" << value << '\n';
                ++checked;
            }
        }
        require(checked == 78, "SU(3) residual scope count is not 78");
        require(have_minimum && minimum == 72 && minimum_a == 2
                    && minimum_n == 1,
                "SU(3) residual minimum identity mismatch");

        // With d=8, c=1, epsilon=2, s=27, the cap probability is at least
        // (m_54-6^54)/(8^54-6^54), the central probability is at least 3/4,
        // and q/(2c)=2.  Check that their degree-29 product exceeds one.
        const mpz_class cap_numerator = moment[54] - power_ui(6, 54);
        const mpz_class cap_denominator = power_ui(8, 54) - power_ui(6, 54);
        require(cap_numerator > 0, "SU(3) cap moment numerator is not positive");
        const mpz_class tail_margin =
            3 * cap_numerator * power_ui(2, 29) - 4 * cap_denominator;
        require(tail_margin > 0,
                "SU(3) exact moment-cap tail margin is not positive");

        std::cout << "SU3_FULL_Q3 m54=" << moment[54] << '\n';
        std::cout << "SU3_FULL_Q3 residual_cases=" << checked << '\n';
        std::cout << "SU3_FULL_Q3 minimum=" << minimum
                  << " at_a=" << minimum_a
                  << " at_n=" << minimum_n << '\n';
        std::cout << "SU3_FULL_Q3 tail_rational_margin="
                  << tail_margin << '\n';
        std::cout << "SU3_FULL_Q3 VERIFICATION: ALL PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU3_FULL_Q3 VERIFICATION: FAIL: " << error.what()
                  << '\n';
        return EXIT_FAILURE;
    }
}
