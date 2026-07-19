#include <gmpxx.h>

#include <algorithm>
#include <array>
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

std::vector<mpz_class> adjoint_moments(
    unsigned rank,
    unsigned maximum_degree
) {
    std::vector<mpz_class> trace(maximum_degree + 1);
    std::vector<mpz_class> adjoint(maximum_degree + 1);
    for (unsigned degree = 0; degree <= maximum_degree; ++degree) {
        trace[degree] = trace_moment(degree, rank);
        if (degree <= rank) {
            require(trace[degree] == factorial(degree),
                    "stable trace-moment identity failed at rank="
                    + std::to_string(rank) + ", degree="
                    + std::to_string(degree));
        }
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

unsigned tail_onset(unsigned rank) {
    if (rank == 6) {
        return 27;
    }
    if (rank == 7 || rank == 8) {
        return 29;
    }
    return 31;
}

}  // namespace

int main() {
    try {
        constexpr std::array<unsigned, 23> expected_counts = {
            65, 75, 75, 85, 85, 81, 81, 76, 76, 70, 70, 63,
            63, 55, 55, 46, 46, 36, 36, 25, 25, 13, 13
        };

        unsigned total_checked = 0;
        mpz_class global_minimum;
        unsigned minimum_rank = 0;
        unsigned minimum_a = 0;
        unsigned minimum_n = 0;
        bool have_minimum = false;

        for (unsigned rank = 6; rank <= 28; ++rank) {
            const std::vector<mpz_class> moment = adjoint_moments(rank, 29);
            const unsigned onset = tail_onset(rank);
            unsigned rank_checked = 0;
            for (unsigned a = 2; 2 * a + 1 <= 29; ++a) {
                for (unsigned n = 1; 2 * a + n <= 29; n += 2) {
                    const unsigned total_degree = 2 * a + n;
                    if (total_degree <= rank || total_degree >= onset) {
                        continue;
                    }
                    const mpz_class value = full_q3_value(moment, a, n);
                    require(value > 0,
                            "nonpositive SU(N) residual at N="
                            + std::to_string(rank) + ", a="
                            + std::to_string(a) + ", n="
                            + std::to_string(n));
                    if (!have_minimum || value < global_minimum) {
                        global_minimum = value;
                        minimum_rank = rank;
                        minimum_a = a;
                        minimum_n = n;
                        have_minimum = true;
                    }
                    ++rank_checked;
                    ++total_checked;
                }
            }
            require(rank_checked == expected_counts[rank - 6],
                    "SU(N) residual rank-scope count mismatch at N="
                    + std::to_string(rank));
            std::cout << "SUN_GE6_FULL_Q3 N=" << rank
                      << " tail_onset=" << onset
                      << " residual_cases=" << rank_checked << '\n';
        }

        require(total_checked == 1315,
                "SU(N>=6) total residual scope count is not 1315");
        require(have_minimum && global_minimum == 3550
                    && minimum_rank == 6 && minimum_a == 2
                    && minimum_n == 3,
                "SU(N>=6) global residual minimum identity mismatch");

        // Paley--Zygmund rectangle: Z=|tr U|^2, W=Z^16,
        // threshold Z>33/5, central window |X|<11/10, and ratio 9/4.
        // Check N=6..31 directly.  N>=32 has E_16=16!, E_32=32!.
        for (unsigned rank = 6; rank <= 32; ++rank) {
            const mpz_class e16 =
                rank >= 16 ? factorial(16) : trace_moment(16, rank);
            const mpz_class e32 =
                rank >= 32 ? factorial(32) : trace_moment(32, rank);
            const mpz_class cap_numerator =
                power_ui(5, 16) * e16 - power_ui(33, 16);
            require(cap_numerator > 0,
                    "SU(N) Paley cap numerator is not positive at N="
                    + std::to_string(rank));
            const unsigned onset = tail_onset(rank);
            const mpz_class tail_margin =
                21 * cap_numerator * cap_numerator * power_ui(9, onset)
                - 121 * power_ui(5, 32) * e32 * power_ui(4, onset);
            require(tail_margin > 0,
                    "SU(N) Paley tail margin is not positive at N="
                    + std::to_string(rank));
            std::cout << "SUN_GE6_FULL_Q3_TAIL N=" << rank
                      << " onset=" << onset
                      << " margin=" << tail_margin << '\n';
        }

        std::cout << "SUN_GE6_FULL_Q3 residual_cases=" << total_checked
                  << '\n';
        std::cout << "SUN_GE6_FULL_Q3 minimum=" << global_minimum
                  << " at_N=" << minimum_rank
                  << " at_a=" << minimum_a
                  << " at_n=" << minimum_n << '\n';
        std::cout << "SUN_GE6_FULL_Q3 VERIFICATION: ALL PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SUN_GE6_FULL_Q3 VERIFICATION: FAIL: " << error.what()
                  << '\n';
        return EXIT_FAILURE;
    }
}
