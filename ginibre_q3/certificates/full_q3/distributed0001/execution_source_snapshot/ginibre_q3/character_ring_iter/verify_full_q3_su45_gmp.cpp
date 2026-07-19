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

mpz_class integer(long long value) {
    return mpz_class(std::to_string(value));
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

std::vector<mpz_class> adjoint_moments(
    unsigned rank,
    unsigned maximum_degree
) {
    std::vector<mpz_class> trace(maximum_degree + 1);
    std::vector<mpz_class> adjoint(maximum_degree + 1);
    for (unsigned degree = 0; degree <= maximum_degree; ++degree) {
        trace[degree] = trace_moment(degree, rank);
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

void verify_su4_recurrence(const std::vector<mpz_class>& m) {
    for (std::size_t index = 0; index + 4 < m.size(); ++index) {
        const long long k = static_cast<long long>(index);
        const long long k2 = k * k;
        const long long k3 = k2 * k;
        const mpz_class value =
            integer(k3 + 22 * k2 + 161 * k + 392) * m[index + 4]
            - integer(16 * k3 + 230 * k2 + 1054 * k + 1524) * m[index + 3]
            + integer(10 * k3 - 340 * k - 750) * m[index + 2]
            + integer(72 * k3 + 522 * k2 + 1242 * k + 972) * m[index + 1]
            + integer(45 * k3 + 270 * k2 + 495 * k + 270) * m[index];
        require(value == 0,
                "SU(4) recurrence mismatch at k=" + std::to_string(k));
    }
}

void verify_su5_recurrence(const std::vector<mpz_class>& m) {
    for (std::size_t index = 0; index + 4 < m.size(); ++index) {
        const long long k = static_cast<long long>(index);
        const long long k2 = k * k;
        const long long k3 = k2 * k;
        const mpz_class value =
            integer(192 * (k + 1) * (k + 2) * (k + 3) * (k + 6)) * m[index]
            + integer(16 * (k + 2) * (k + 3)
                      * (22 * k2 + 184 * k + 309)) * m[index + 1]
            + integer((k + 3) * (129 * k3 + 1245 * k2 + 2570 * k - 1112))
                * m[index + 2]
            - integer((k + 3) * (30 * k3 + 642 * k2 + 4487 * k + 10196))
                * m[index + 3]
            + integer((k + 8) * (k + 8) * (k + 10) * (k + 10))
                * m[index + 4];
        require(value == 0,
                "SU(5) recurrence mismatch at k=" + std::to_string(k));
    }
}

void verify_residual(
    const std::string& label,
    const std::vector<mpz_class>& moment,
    unsigned maximum_total_degree,
    unsigned expected_count,
    const mpz_class& expected_minimum
) {
    unsigned checked = 0;
    mpz_class minimum;
    unsigned minimum_a = 0;
    unsigned minimum_n = 0;
    bool have_minimum = false;
    for (unsigned a = 2; 2 * a + 1 <= maximum_total_degree; ++a) {
        for (unsigned n = 1; 2 * a + n <= maximum_total_degree; n += 2) {
            const mpz_class value = full_q3_value(moment, a, n);
            require(value > 0,
                    "nonpositive " + label + " residual at a="
                    + std::to_string(a) + ", n=" + std::to_string(n));
            if (!have_minimum || value < minimum) {
                minimum = value;
                minimum_a = a;
                minimum_n = n;
                have_minimum = true;
            }
            ++checked;
        }
    }
    require(checked == expected_count,
            label + " residual scope count mismatch");
    require(have_minimum && minimum == expected_minimum
                && minimum_a == 2 && minimum_n == 1,
            label + " residual minimum identity mismatch");
    std::cout << label << " residual_cases=" << checked << '\n';
    std::cout << label << " minimum=" << minimum
              << " at_a=" << minimum_a << " at_n=" << minimum_n << '\n';
}

}  // namespace

int main() {
    try {
        const std::vector<mpz_class> su4 = adjoint_moments(4, 94);
        const std::vector<mpz_class> su5 = adjoint_moments(5, 96);

        const std::vector<unsigned long> su4_prefix =
            {1, 0, 1, 2, 9, 43, 245, 1557, 10829, 80958};
        const std::vector<unsigned long> su5_prefix =
            {1, 0, 1, 2, 9, 44, 264, 1824, 14210, 122224};
        for (unsigned i = 0; i < su4_prefix.size(); ++i) {
            require(su4[i] == su4_prefix[i],
                    "SU(4) hook-length prefix mismatch at degree "
                    + std::to_string(i));
            require(su5[i] == su5_prefix[i],
                    "SU(5) hook-length prefix mismatch at degree "
                    + std::to_string(i));
        }
        verify_su4_recurrence(su4);
        verify_su5_recurrence(su5);

        const mpz_class expected_m94(
            "80069447965203643835438443611006096875924773400360563146060175988854883321929292812182505855921545"
        );
        const mpz_class expected_m96(
            "2005370457065724965367030341374673134939501180557415671389220080740234301685695670596195463720943077657958035456616"
        );
        require(su4[94] == expected_m94, "SU(4) degree-94 moment mismatch");
        require(su5[96] == expected_m96, "SU(5) degree-96 moment mismatch");

        verify_residual("SU4_FULL_Q3", su4, 25, 66, mpz_class(94));
        verify_residual("SU5_FULL_Q3", su5, 35, 136, mpz_class(96));

        const mpz_class su4_cap_num = su4[94] - power_ui(11, 94);
        const mpz_class su4_cap_den = power_ui(15, 94) - power_ui(11, 94);
        const mpz_class su4_tail_margin =
            15 * su4_cap_num * power_ui(7, 27)
            - 16 * su4_cap_den * power_ui(2, 27);
        require(su4_cap_num > 0 && su4_tail_margin > 0,
                "SU(4) exact moment-cap tail check failed");

        const mpz_class su5_cap_num =
            power_ui(2, 96) * su5[96] - power_ui(31, 96);
        const mpz_class su5_cap_den = power_ui(48, 96) - power_ui(31, 96);
        const mpz_class su5_tail_margin =
            285 * su5_cap_num * power_ui(7, 37)
            - 289 * su5_cap_den * power_ui(2, 37);
        require(su5_cap_num > 0 && su5_tail_margin > 0,
                "SU(5) exact moment-cap tail check failed");

        std::cout << "SU4_FULL_Q3 m94=" << su4[94] << '\n';
        std::cout << "SU4_FULL_Q3 tail_rational_margin="
                  << su4_tail_margin << '\n';
        std::cout << "SU5_FULL_Q3 m96=" << su5[96] << '\n';
        std::cout << "SU5_FULL_Q3 tail_rational_margin="
                  << su5_tail_margin << '\n';
        std::cout << "SU45_FULL_Q3 VERIFICATION: ALL PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SU45_FULL_Q3 VERIFICATION: FAIL: " << error.what()
                  << '\n';
        return EXIT_FAILURE;
    }
}
