#include <gmpxx.h>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Integer = mpz_class;
using Rational = mpq_class;
using Polynomial = std::vector<Rational>;

constexpr const char* SOURCE_SHA256 =
    "bc22efcafafe2d40b5f5430cb8f9cb046bc6b9b188da684e61dc6fdfaa27077e";
constexpr const char* SOURCE_HOST = "nb1cb2f";

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

Integer factorial(unsigned long n) {
    Integer answer;
    mpz_fac_ui(answer.get_mpz_t(), n);
    return answer;
}

Integer binomial(unsigned long n, unsigned long k) {
    Integer answer;
    mpz_bin_uiui(answer.get_mpz_t(), n, k);
    return answer;
}

Rational power(Rational base, unsigned long exponent) {
    Rational answer = 1;
    while (exponent != 0) {
        if ((exponent & 1UL) != 0) answer *= base;
        exponent >>= 1U;
        if (exponent != 0) base *= base;
    }
    return answer;
}

Integer power(Integer base, unsigned long exponent) {
    Integer answer;
    mpz_pow_ui(answer.get_mpz_t(), base.get_mpz_t(), exponent);
    return answer;
}

struct TailRow {
    int rank;
    int dimension;
    int kappa;
    int c_value;
    Rational radius;
    Rational central_width;
    Integer sine_denominator;
    int exp_terms;
    int chebyshev_degree;
    Rational chebyshev_scale;
    Rational negative_cutoff;
    int maximum_adjoint_moment;
    std::optional<Rational> radial_sqrt_lower;
};

std::array<TailRow, 1> rows() {
    return {{
        {
            10, 190, 18, 10, Rational(107), Rational(6, 5),
            Integer(12400), 90, 9, Rational(4650), Rational(14), 38,
            std::nullopt,
        },
    }};
}

std::vector<Integer> read_stable(const std::string& path, int maximum) {
    std::ifstream input(path);
    if (!input) fail("cannot open stable moment table: " + path);
    std::map<int, Integer> parsed;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::istringstream fields(line);
        int index = -1;
        Integer value;
        if (!(fields >> index >> value)) continue;
        if (index <= maximum) parsed[index] = value;
    }
    std::vector<Integer> stable(static_cast<std::size_t>(maximum + 1));
    for (int index = 0; index <= maximum; ++index) {
        const auto found = parsed.find(index);
        if (found == parsed.end()) {
            fail("stable moment table is incomplete through m_"
                 + std::to_string(maximum));
        }
        stable[static_cast<std::size_t>(index)] = found->second;
    }
    if (parsed.size() != stable.size()) {
        fail("stable moment table has duplicate or unexpected indices");
    }
    return stable;
}

using MomentMap = std::map<int, std::map<int, Integer>>;

MomentMap parse_sources(
    const std::vector<std::string>& paths,
    const std::vector<Integer>& stable
) {
    if (paths.empty()) fail("at least one exact source log is required");
    std::map<std::string, bool> unique_paths;
    MomentMap moments;
    const std::regex source_row(
        R"(^\s*D_(\d+) m_(\d+) \+= O_even (-?\d+) \+ det (-?\d+); Delta_(\d+) = (-?\d+); moment_(\d+) = (\d+)$)"
    );

    for (const std::string& path : paths) {
        if (!unique_paths.emplace(path, true).second) {
            fail("duplicate exact source log: " + path);
        }
        std::ifstream input(path);
        if (!input) fail("cannot open exact source log: " + path);
        int host_rows = 0;
        int source_hash_rows = 0;
        int success_rows = 0;
        int parsed_rows = 0;
        std::string line;
        while (std::getline(input, line)) {
            if (line == std::string("__HOST=") + SOURCE_HOST) ++host_rows;
            if (line == std::string("__SOURCE_SHA256=") + SOURCE_SHA256) {
                ++source_hash_rows;
            }
            if (line == "__EXIT_STATUS=0") ++success_rows;
            if (line.rfind("__EXIT_STATUS=", 0) == 0
                && line != "__EXIT_STATUS=0") {
                fail("nonzero source status in " + path);
            }

            std::smatch match;
            if (!std::regex_match(line, match, source_row)) continue;
            const int rank = std::stoi(match[1].str());
            const int index = std::stoi(match[2].str());
            const Integer even_overflow(match[3].str());
            const Integer determinant(match[4].str());
            const int delta_index = std::stoi(match[5].str());
            const Integer delta(match[6].str());
            const int moment_index = std::stoi(match[7].str());
            const Integer moment(match[8].str());
            if (index != delta_index || index != moment_index) {
                fail("mismatched exact source indices in " + path);
            }
            if (index < 0 || index >= static_cast<int>(stable.size())) {
                fail("exact source moment exceeds stable table in " + path);
            }
            if (even_overflow + determinant != delta) {
                fail("exact source correction identity failed in " + path);
            }
            if (stable[static_cast<std::size_t>(index)] + delta != moment) {
                fail("exact source moment identity failed in " + path);
            }
            auto& rank_moments = moments[rank];
            const auto found = rank_moments.find(index);
            if (found != rank_moments.end() && found->second != moment) {
                fail("conflicting exact source moment in " + path);
            }
            rank_moments[index] = moment;
            ++parsed_rows;
        }
        if (host_rows != 1 || source_hash_rows != 1 || success_rows != 1) {
            fail("source provenance/status is not unique and successful: " + path);
        }
        if (parsed_rows == 0) fail("source log has no exact moment rows: " + path);
    }
    return moments;
}

Polynomial add(const Polynomial& first, const Polynomial& second) {
    Polynomial answer(std::max(first.size(), second.size()), Rational(0));
    for (std::size_t i = 0; i < first.size(); ++i) answer[i] += first[i];
    for (std::size_t i = 0; i < second.size(); ++i) answer[i] += second[i];
    return answer;
}

Polynomial scale(const Polynomial& polynomial, const Rational& factor) {
    Polynomial answer = polynomial;
    for (Rational& coefficient : answer) coefficient *= factor;
    return answer;
}

Polynomial multiply(const Polynomial& first, const Polynomial& second) {
    Polynomial answer(first.size() + second.size() - 1, Rational(0));
    for (std::size_t i = 0; i < first.size(); ++i) {
        for (std::size_t j = 0; j < second.size(); ++j) {
            answer[i + j] += first[i] * second[j];
        }
    }
    return answer;
}

Polynomial chebyshev_t(int degree) {
    if (degree == 0) return {Rational(1)};
    if (degree == 1) return {Rational(0), Rational(1)};
    Polynomial previous{Rational(1)};
    Polynomial current{Rational(0), Rational(1)};
    for (int index = 2; index <= degree; ++index) {
        Polynomial shifted(current.size() + 1, Rational(0));
        for (std::size_t j = 0; j < current.size(); ++j) {
            shifted[j + 1] = 2 * current[j];
        }
        Polynomial next = add(shifted, scale(previous, Rational(-1)));
        previous = std::move(current);
        current = std::move(next);
    }
    return current;
}

Polynomial compose(const Polynomial& polynomial, const Polynomial& argument) {
    Polynomial answer{Rational(0)};
    Polynomial argument_power{Rational(1)};
    for (const Rational& coefficient : polynomial) {
        answer = add(answer, scale(argument_power, coefficient));
        argument_power = multiply(argument_power, argument);
    }
    return answer;
}

Rational evaluate(const Polynomial& polynomial, const Rational& value) {
    Rational answer = 0;
    Rational value_power = 1;
    for (const Rational& coefficient : polynomial) {
        answer += coefficient * value_power;
        value_power *= value;
    }
    return answer;
}

Rational exp_upper_bound(const Rational& value, int terms) {
    Rational term = 1;
    Rational total = term;
    for (int index = 1; index <= terms; ++index) {
        term *= value / index;
        total += term;
    }
    const Rational next_term = term * value / (terms + 1);
    const Rational ratio_bound = value / (terms + 2);
    if (ratio_bound >= 1) fail("exponential Taylor tail ratio is not below one");
    return total + next_term / (1 - ratio_bound);
}

Rational half_integer_density_factor(int rank, int dimension) {
    if (rank % 2 != 1 || dimension % 2 != 1) {
        fail("half-integer density factor used with even parity");
    }
    if (power(Rational(10, 7), 2) <= 2) {
        fail("rational upper bound sqrt(2)<10/7 failed");
    }
    const int half_index = (dimension + 1) / 2;
    const int pi_power = (rank + 1) / 2;
    return Rational(7, 5)
        * Rational(power(Integer(4), static_cast<unsigned long>(half_index))
                   * factorial(static_cast<unsigned long>(half_index)))
        / (power(Rational(44, 7), static_cast<unsigned long>(pi_power))
           * factorial(static_cast<unsigned long>(2 * half_index)));
}

Integer q3_moment(const std::vector<Integer>& moments, int degree) {
    Integer answer = 0;
    for (int index = 0; index <= degree; ++index) {
        answer += binomial(
            static_cast<unsigned long>(degree),
            static_cast<unsigned long>(index)
        ) * (
            moments[static_cast<std::size_t>(index + 2)]
                * moments[static_cast<std::size_t>(degree - index)]
            - moments[static_cast<std::size_t>(index + 1)]
                * moments[static_cast<std::size_t>(degree - index + 1)]
        );
    }
    return 2 * answer;
}

void print_rational(const std::string& label, const Rational& value) {
    std::cout << label << '=' << value.get_str() << '\n';
}

void verify_row(
    const TailRow& row,
    const MomentMap& source,
    const std::vector<Integer>& stable
) {
    if (row.rank < 4
        || row.dimension != row.rank * (2 * row.rank - 1)
        || row.kappa != 2 * row.rank - 2) {
        fail("type-D rank/dimension/root-normalization schedule failed");
    }
    if (row.radius <= 0 || row.central_width <= 0
        || row.sine_denominator <= 1 || row.exp_terms <= 0
        || row.chebyshev_degree <= 0 || row.chebyshev_scale <= 0
        || row.negative_cutoff <= 0) {
        fail("nonpositive certificate parameter at D_"
             + std::to_string(row.rank));
    }
    const auto rank_source = source.find(row.rank);
    if (rank_source == source.end()) {
        fail("missing source moments for D_" + std::to_string(row.rank));
    }
    std::vector<Integer> moments(
        stable.begin(),
        stable.begin() + row.maximum_adjoint_moment + 1
    );
    for (int index = row.rank; index <= row.maximum_adjoint_moment; ++index) {
        const auto found = rank_source->second.find(index);
        if (found == rank_source->second.end()) {
            fail("D_" + std::to_string(row.rank)
                 + " source has a gap at m_" + std::to_string(index));
        }
        moments[static_cast<std::size_t>(index)] = found->second;
    }

    const Rational q_bound = 2 * row.radius / (39 * row.kappa);
    if (!(q_bound > 0 && q_bound < 1)) {
        fail("Euler-product remainder ratio is invalid at D_"
             + std::to_string(row.rank));
    }
    const Rational sine_exponent =
        row.radius / 12
        + power(row.radius, 2) / (1440 * row.kappa)
        + power(row.radius, 3)
            / (Integer(row.kappa * row.kappa) * 90720 * (1 - q_bound));
    const Rational exponential_upper = exp_upper_bound(
        sine_exponent,
        row.exp_terms
    );
    if (exponential_upper >= row.sine_denominator) {
        fail("exponential sine bound failed at D_" + std::to_string(row.rank));
    }

    const Rational quartic_r(1, 2);
    const Rational quartic_s(15, 4);
    const Rational quartic_moment =
        6 - 2 * (quartic_r + quartic_s)
        + power(quartic_r, 2) + power(quartic_s, 2)
        + 4 * quartic_r * quartic_s
        + power(quartic_r, 2) * power(quartic_s, 2);
    const Rational quartic_denominator =
        power(row.central_width + quartic_r, 2)
        * power(row.central_width + quartic_s, 2);
    const Rational x_floor = row.dimension - row.radius;
    const Rational c_push = x_floor - row.central_width;
    if (!(row.negative_cutoff < 2 * row.c_value
          && 2 * row.c_value < c_push)) {
        fail("pushforward threshold ordering a<2C<c failed at D_"
             + std::to_string(row.rank));
    }
    const Rational tail_weight_upper =
        power(x_floor + row.central_width, 2)
        * quartic_moment / quartic_denominator;
    const Rational pushforward_weight_lower =
        power(x_floor, 2) + 1 - tail_weight_upper;
    const Rational pushforward_derivative_lower =
        2 * x_floor
        - 2 * quartic_moment * (x_floor + row.central_width)
            / quartic_denominator;
    const Rational pushforward_second_derivative =
        2 - 2 * quartic_moment / quartic_denominator;
    if (std::min({
            c_push,
            pushforward_weight_lower,
            pushforward_derivative_lower,
            pushforward_second_derivative,
        }) <= 0) {
        fail("positive-cap monotonicity failed at D_" + std::to_string(row.rank));
    }

    Integer degree_factor = 1;
    for (int degree = 2; degree <= 2 * (row.rank - 1); degree += 2) {
        degree_factor *= factorial(static_cast<unsigned long>(degree));
    }
    degree_factor *= factorial(static_cast<unsigned long>(row.rank));
    const Integer weyl_order =
        power(Integer(2), static_cast<unsigned long>(row.rank - 1))
        * factorial(static_cast<unsigned long>(row.rank));

    Rational density_factor;
    Rational radial_factor;
    if (!row.radial_sqrt_lower.has_value()) {
        if (row.rank % 2 != 0 || row.dimension % 2 != 0) {
            fail("even density parity failed at D_" + std::to_string(row.rank));
        }
        density_factor = Rational(1) /
            (power(Rational(44, 7), static_cast<unsigned long>(row.rank / 2))
             * factorial(static_cast<unsigned long>(row.dimension / 2)));
        radial_factor = power(
            row.radius / (2 * row.kappa),
            static_cast<unsigned long>(row.dimension / 2)
        );
    } else {
        if (power(*row.radial_sqrt_lower, 2) > row.radius / (2 * row.kappa)) {
            fail("radial square-root lower bound failed at D_"
                 + std::to_string(row.rank));
        }
        density_factor = half_integer_density_factor(row.rank, row.dimension);
        radial_factor = power(
            row.radius / (2 * row.kappa),
            static_cast<unsigned long>((row.dimension - 1) / 2)
        ) * *row.radial_sqrt_lower;
    }
    const Rational cap_probability_lower =
        Rational(1, weyl_order)
        * Rational(1, row.sine_denominator)
        * degree_factor * density_factor * radial_factor;
    const Rational positive_tail_lower =
        pushforward_weight_lower * cap_probability_lower;

    const Rational support_width(2 * row.dimension);
    const Rational q0 = row.negative_cutoff
        * (row.negative_cutoff + support_width);
    const Polynomial chebyshev = chebyshev_t(row.chebyshev_degree);
    const Polynomial argument{
        Rational(1),
        -2 * support_width / row.chebyshev_scale,
        Rational(2) / row.chebyshev_scale,
    };
    const Polynomial composed = compose(chebyshev, argument);
    const Polynomial majorant = multiply(composed, composed);
    if (static_cast<int>(majorant.size()) - 1 != 4 * row.chebyshev_degree) {
        fail("Chebyshev degree failed at D_" + std::to_string(row.rank));
    }
    if (static_cast<int>(majorant.size()) + 1 != row.maximum_adjoint_moment) {
        fail("source-degree schedule failed at D_" + std::to_string(row.rank));
    }
    const Rational denominator = power(
        evaluate(
            chebyshev,
            1 + 2 * q0 / row.chebyshev_scale
        ),
        2
    );
    if (denominator <= 0) {
        fail("Chebyshev denominator failed at D_" + std::to_string(row.rank));
    }
    Rational numerator = 0;
    for (std::size_t index = 0; index < majorant.size(); ++index) {
        numerator += majorant[index]
            * q3_moment(moments, static_cast<int>(index));
    }
    if (numerator <= 0) {
        fail("Chebyshev expectation failed at D_" + std::to_string(row.rank));
    }
    const Rational negative_tail_mass_upper = numerator / denominator;

    constexpr unsigned long TARGET_EXPONENT = 63;
    const Rational sharpened_target =
        negative_tail_mass_upper
            * power(Rational(2 * row.c_value), TARGET_EXPONENT)
            / power(c_push, TARGET_EXPONENT)
        + 2 * power(row.negative_cutoff, TARGET_EXPONENT)
            / power(c_push, TARGET_EXPONENT);
    const Rational margin = positive_tail_lower - sharpened_target;
    if (margin <= 0) {
        std::cerr << "DIAGNOSTIC row=D_" << row.rank << '\n';
        std::cerr << "negative_tail_mass_upper="
                  << negative_tail_mass_upper.get_str() << '\n';
        std::cerr << "positive_tail_lower="
                  << positive_tail_lower.get_str() << '\n';
        std::cerr << "sharpened_target=" << sharpened_target.get_str() << '\n';
        std::cerr << "margin=" << margin.get_str() << '\n';
        fail("n=63 exact rational tail margin is nonpositive at D_"
             + std::to_string(row.rank));
    }

    std::cout << "row=D_" << row.rank
              << " target_exponent=" << TARGET_EXPONENT
              << " source_moments=m_" << row.rank << "..m_"
              << row.maximum_adjoint_moment
              << " push_moments=K_0..K_" << (majorant.size() - 1)
              << " chebyshev_degree=" << row.chebyshev_degree
              << " exact_margin_sign=1\n";
    print_rational("negative_tail_mass_upper", negative_tail_mass_upper);
    print_rational("positive_tail_lower", positive_tail_lower);
    print_rational("sharpened_target", sharpened_target);
    print_rational("margin", margin);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::vector<std::string> source_paths;
        std::string stable_path = "../references/oeis_A002137_stable.txt";
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--source-log" && index + 1 < argc) {
                source_paths.emplace_back(argv[++index]);
            } else if (argument == "--stable-moments" && index + 1 < argc) {
                stable_path = argv[++index];
            } else {
                fail("usage: verify_d10_exact_tail_gmp "
                     "--source-log PATH [--source-log PATH ...] "
                     "[--stable-moments PATH]");
            }
        }
        constexpr int MAXIMUM = 38;
        const std::vector<Integer> stable = read_stable(stable_path, MAXIMUM);
        const MomentMap source = parse_sources(source_paths, stable);
        std::cout << "D10 full-adjoint-moment exact rational tail certificate\n"
                  << "arithmetic=exact_GMP_integer_and_rational\n"
                  << "stable_moments=" << stable_path << '\n'
                  << "source_logs=" << source_paths.size() << '\n';
        for (const std::string& path : source_paths) {
            std::cout << "source_log=" << path << '\n';
        }
        for (const TailRow& row : rows()) verify_row(row, source, stable);
        std::cout << "SUMMARY rows_checked=1 failures=0\n"
                  << "__EXIT_STATUS=0\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 2;
    }
}
