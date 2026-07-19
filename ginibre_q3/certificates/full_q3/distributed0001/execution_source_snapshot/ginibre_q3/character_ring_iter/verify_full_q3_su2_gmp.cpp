#include <gmpxx.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

void require(bool condition, const std::string& message);

mpz_class binomial(unsigned long n, unsigned long k) {
    mpz_class answer;
    mpz_bin_uiui(answer.get_mpz_t(), n, k);
    return answer;
}

mpz_class power_ui(unsigned long base, unsigned long exponent) {
    mpz_class answer;
    mpz_ui_pow_ui(answer.get_mpz_t(), base, exponent);
    return answer;
}

std::vector<mpz_class> su2_adjoint_moments(unsigned maximum_degree) {
    // c_k(j) is the multiplicity of spin j in (spin 1)^{tensor k}.
    // Clebsch--Gordan gives
    // c_{k+1}(0)=c_k(1) and
    // c_{k+1}(j)=c_k(j-1)+c_k(j)+c_k(j+1) for j>=1.
    std::vector<mpz_class> multiplicity(maximum_degree + 2);
    std::vector<mpz_class> next(maximum_degree + 2);
    std::vector<mpz_class> moment;
    moment.reserve(maximum_degree + 1);

    multiplicity[0] = 1;
    moment.push_back(1);
    for (unsigned degree = 0; degree < maximum_degree; ++degree) {
        std::fill(next.begin(), next.end(), 0);
        next[0] = multiplicity[1];
        for (unsigned spin = 1; spin <= maximum_degree; ++spin) {
            next[spin] = multiplicity[spin - 1] + multiplicity[spin]
                         + multiplicity[spin + 1];
        }
        multiplicity.swap(next);
        moment.push_back(multiplicity[0]);
    }
    return moment;
}

std::vector<mpz_class> su2_adjoint_moments_catalan(
    unsigned maximum_degree
) {
    // Independently, |tr U|^{2j} has Catalan moment C_j for Haar SU(2),
    // and chi_ad=|tr U|^2-1.
    std::vector<mpz_class> catalan(maximum_degree + 1);
    std::vector<mpz_class> moment(maximum_degree + 1);
    for (unsigned j = 0; j <= maximum_degree; ++j) {
        catalan[j] = binomial(2 * j, j) / (j + 1);
        require(catalan[j] * (j + 1) == binomial(2 * j, j),
                "nonintegral Catalan reconstruction");
    }
    for (unsigned degree = 0; degree <= maximum_degree; ++degree) {
        mpz_class value = 0;
        for (unsigned j = 0; j <= degree; ++j) {
            const mpz_class term = binomial(degree, j) * catalan[j];
            if ((degree - j) & 1U) {
                value -= term;
            } else {
                value += term;
            }
        }
        moment[degree] = value;
    }
    return moment;
}

mpz_class full_q3_value(
    const std::vector<mpz_class>& moment,
    unsigned a,
    unsigned n
) {
    mpz_class answer = 0;
    for (unsigned j = 0; j <= 2 * a; ++j) {
        for (unsigned k = 0; k <= n; ++k) {
            mpz_class term = binomial(2 * a, j) * binomial(n, k)
                           * moment[2 * a - j + k]
                           * moment[j + n - k];
            if (j & 1U) {
                answer -= term;
            } else {
                answer += term;
            }
        }
    }
    return answer;
}

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

}  // namespace

int main() {
    try {
        constexpr unsigned maximum_degree = 27;
        const std::vector<mpz_class> moment =
            su2_adjoint_moments(maximum_degree);
        const std::vector<mpz_class> catalan_moment =
            su2_adjoint_moments_catalan(maximum_degree);

        const std::vector<unsigned long> expected_prefix =
            {1, 0, 1, 1, 3, 6, 15, 36, 91, 232};
        require(moment.size() == maximum_degree + 1,
                "SU(2) moment generator returned the wrong length");
        require(moment == catalan_moment,
                "Clebsch--Gordan and Catalan SU(2) moments disagree");
        for (unsigned index = 0; index < expected_prefix.size(); ++index) {
            require(moment[index] == expected_prefix[index],
                    "SU(2) Clebsch--Gordan moment prefix mismatch at degree "
                    + std::to_string(index));
        }

        unsigned checked = 0;
        mpz_class minimum;
        unsigned minimum_a = 0;
        unsigned minimum_n = 0;
        bool have_minimum = false;

        for (unsigned a = 2; 2 * a + 1 <= maximum_degree; ++a) {
            for (unsigned n = 1; 2 * a + n <= maximum_degree; n += 2) {
                const mpz_class value = full_q3_value(moment, a, n);
                require(value > 0,
                        "nonpositive SU(2) full-Q3 residual at a="
                        + std::to_string(a) + ", n=" + std::to_string(n)
                        + ", value=" + value.get_str());
                if (!have_minimum || value < minimum) {
                    minimum = value;
                    minimum_a = a;
                    minimum_n = n;
                    have_minimum = true;
                }
                std::cout << "SU2_FULL_Q3 a=" << a
                          << " n=" << n
                          << " total_degree=" << (2 * a + n)
                          << " value=" << value << '\n';
                ++checked;
            }
        }

        require(checked == 78, "SU(2) residual scope count is not 78");
        require(have_minimum && minimum == 16 && minimum_a == 2
                    && minimum_n == 1,
                "SU(2) residual minimum identity mismatch");

        // Exact comparison used by the explicit Haar-rectangle tail:
        // 3789107/38106288000 > (8/11)^29.
        const mpz_class tail_margin =
            mpz_class(3789107) * power_ui(11, 29)
            - mpz_class(38106288000UL) * power_ui(8, 29);
        require(tail_margin > 0,
                "SU(2) exact Haar-rectangle tail margin is not positive");

        std::cout << "SU2_FULL_Q3 residual_cases=" << checked << '\n';
        std::cout << "SU2_FULL_Q3 minimum=" << minimum
                  << " at_a=" << minimum_a
                  << " at_n=" << minimum_n << '\n';
        std::cout << "SU2_FULL_Q3 tail_rational_margin="
                  << tail_margin << '\n';
        std::cout << "SU2_FULL_Q3 VERIFICATION: ALL PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU2_FULL_Q3 VERIFICATION: FAIL: " << error.what()
                  << '\n';
        return EXIT_FAILURE;
    }
}
