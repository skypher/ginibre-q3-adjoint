#include <gmpxx.h>

#include <array>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
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

mpq_class rational(long numerator, long denominator) {
    require(denominator > 0, "nonpositive rational denominator");
    return mpq_class(mpz_class(numerator), mpz_class(denominator));
}

mpq_class power(mpq_class base, unsigned exponent) {
    mpq_class answer = 1;
    while (exponent != 0U) {
        if ((exponent & 1U) != 0U) {
            answer *= base;
        }
        exponent >>= 1U;
        if (exponent != 0U) {
            base *= base;
        }
    }
    return answer;
}

void read_moments(
    const std::string& path,
    std::map<unsigned, mpz_class>& moments
) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open moment source: " + path);
    std::string line;
    while (std::getline(input, line)) {
        std::size_t search = 0;
        while ((search = line.find("m_", search)) != std::string::npos) {
            std::size_t cursor = search + 2U;
            const std::size_t degree_begin = cursor;
            while (cursor < line.size()
                   && std::isdigit(static_cast<unsigned char>(line[cursor])) != 0) {
                ++cursor;
            }
            if (cursor == degree_begin) {
                search += 2U;
                continue;
            }
            const unsigned degree = static_cast<unsigned>(
                std::stoul(line.substr(degree_begin, cursor - degree_begin)));
            while (cursor < line.size()
                   && std::isspace(static_cast<unsigned char>(line[cursor])) != 0) {
                ++cursor;
            }
            if (cursor == line.size() || line[cursor] != '=') {
                search = cursor;
                continue;
            }
            ++cursor;
            while (cursor < line.size()
                   && std::isspace(static_cast<unsigned char>(line[cursor])) != 0) {
                ++cursor;
            }
            const std::size_t value_begin = cursor;
            if (cursor < line.size() && line[cursor] == '-') {
                ++cursor;
            }
            const std::size_t digit_begin = cursor;
            while (cursor < line.size()
                   && std::isdigit(static_cast<unsigned char>(line[cursor])) != 0) {
                ++cursor;
            }
            require(cursor != digit_begin,
                    "malformed moment value in source: " + path);
            const mpz_class value(
                line.substr(value_begin, cursor - value_begin));
            const auto [found, inserted] = moments.emplace(degree, value);
            require(inserted || found->second == value,
                    "conflicting m_" + std::to_string(degree)
                    + " in moment sources");
            search = cursor;
        }
    }
}

std::vector<mpz_class> load_moments(
    const std::vector<std::string>& paths,
    unsigned maximum_degree
) {
    std::map<unsigned, mpz_class> parsed;
    for (const std::string& path : paths) {
        read_moments(path, parsed);
    }
    std::vector<mpz_class> result(maximum_degree + 1U);
    for (unsigned degree = 0; degree <= maximum_degree; ++degree) {
        const auto found = parsed.find(degree);
        require(found != parsed.end(),
                "moment sources omit m_" + std::to_string(degree));
        require(found->second >= 0,
                "negative invariant-tensor moment m_"
                + std::to_string(degree));
        result[degree] = found->second;
    }
    require(result[0] == 1 && result[1] == 0 && result[2] == 1,
            "adjoint normalization moments do not equal 1,0,1");
    return result;
}

mpz_class full_q3_value(
    const std::vector<mpz_class>& moment,
    unsigned a,
    unsigned n
) {
    mpz_class answer = 0;
    for (unsigned j = 0; j <= 2U * a; ++j) {
        for (unsigned ell = 0; ell <= n; ++ell) {
            const mpz_class term =
                binomial(2U * a, j) * binomial(n, ell)
                * moment[2U * a - j + ell] * moment[j + n - ell];
            if ((j & 1U) != 0U) {
                answer -= term;
            } else {
                answer += term;
            }
        }
    }
    return answer;
}

struct GroupCase {
    std::string name;
    unsigned negative_cap;
    long threshold_tenths;
    long central_tenths;
    unsigned polynomial_degree;
    unsigned tail_onset;
    unsigned maximum_moment;
    unsigned expected_residual_cases;
    mpz_class expected_minimum;
    std::vector<std::string> sources;
};

void verify_case(const GroupCase& group) {
    const std::vector<mpz_class> moment =
        load_moments(group.sources, group.maximum_moment);
    const unsigned k = group.polynomial_degree;
    require(4U * k + 2U <= group.maximum_moment,
            group.name + " tail certificate exceeds loaded moment range");

    const mpq_class threshold = rational(group.threshold_tenths, 10);
    const mpq_class central = rational(group.central_tenths, 10);
    const mpq_class negative_base = 2U * group.negative_cap;
    const mpq_class positive_base = threshold - central;
    const mpq_class ratio = positive_base / negative_base;
    require(central > 1 && ratio > 1,
            group.name + " tail geometry is invalid");

    const mpq_class first =
        moment[2U * k + 1U] - threshold * moment[2U * k];
    const mpq_class second =
        moment[4U * k + 2U]
        - 2 * threshold * moment[4U * k + 1U]
        + threshold * threshold * moment[4U * k];
    require(first > 0 && second > 0,
            group.name + " polynomial tail moments are not positive");
    const mpq_class tail_probability =
        first * first / second * (1 - 1 / (central * central));
    require(tail_probability > 0 && tail_probability <= 1,
            group.name + " tail probability lower bound is invalid");

    const mpq_class onset_margin =
        tail_probability * power(ratio, group.tail_onset) - 1;
    require(onset_margin > 0,
            group.name + " exact tail-onset comparison failed");
    require(group.tail_onset < 2U
                || tail_probability * power(ratio, group.tail_onset - 2U) <= 1,
            group.name + " tail onset is not the asserted first odd onset");

    unsigned checked = 0;
    mpz_class minimum;
    unsigned minimum_a = 0;
    unsigned minimum_n = 0;
    bool have_minimum = false;
    for (unsigned a = 2; 2U * a + 1U < group.tail_onset; ++a) {
        for (unsigned n = 1; 2U * a + n < group.tail_onset; n += 2U) {
            const mpz_class value = full_q3_value(moment, a, n);
            require(value > 0,
                    group.name + " nonpositive residual at a="
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
    require(checked == group.expected_residual_cases,
            group.name + " residual scope count mismatch");
    require(have_minimum && minimum == group.expected_minimum
                && minimum_a == 2U && minimum_n == 1U,
            group.name + " residual minimum identity mismatch");

    std::cout << "EXCEPTIONAL_FULL_Q3 group=" << group.name
              << " tail_onset=" << group.tail_onset
              << " p_degree=" << k
              << " maximum_moment=" << group.maximum_moment
              << " residual_cases=" << checked
              << " minimum=" << minimum
              << " onset_margin_num=" << onset_margin.get_num()
              << " onset_margin_den=" << onset_margin.get_den() << '\n';
}

}  // namespace

int main() {
    try {
        const std::array<GroupCase, 5> cases = {{
            {"G2", 2, 91, 12, 9, 25, 38, 55, mpz_class(36),
             {"character_ring_iter/logs/g2_200.log"}},
            {"F4", 4, 211, 13, 12, 37, 50, 136, mpz_class(36),
             {"character_ring_iter/logs/f4_220c.log"}},
            {"E6", 3, 196, 14, 10, 29, 42, 78, mpz_class(38),
             {"character_ring_iter/logs/e6_80.log"}},
            {"E7", 7, 333, 13, 17, 63, 70, 435, mpz_class(36),
             {"character_ring_iter/logs/e7_70.log"}},
            {"E8", 8, 461, 13, 23, 71, 94, 561, mpz_class(36),
             {"character_ring_iter/logs/e8_70.log",
              "references/arxiv_2412_21189_e8_m71_m100.txt"}}
        }};

        unsigned total = 0;
        for (const GroupCase& group : cases) {
            verify_case(group);
            total += group.expected_residual_cases;
        }
        require(total == 1265U, "exceptional total residual scope mismatch");
        std::cout << "EXCEPTIONAL_FULL_Q3 residual_cases=" << total << '\n';
        std::cout << "EXCEPTIONAL_FULL_Q3 VERIFICATION: ALL PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "EXCEPTIONAL_FULL_Q3 VERIFICATION: FAIL: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
