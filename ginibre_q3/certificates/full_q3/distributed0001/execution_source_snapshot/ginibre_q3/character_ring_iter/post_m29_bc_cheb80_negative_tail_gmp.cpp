#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double LOG10_2_DOUBLE = 0.301029995663981195213738894724493027;

using Poly = std::vector<mpq_class>;
using DeltaKey = std::pair<char, int>;
using DeltaMap = std::map<DeltaKey, std::map<int, mpz_class>>;

struct Row {
    char family = 'B';
    int rank = 0;
};

struct RowResult {
    char family = 'B';
    int rank = 0;
    int ok = 0;
    std::string reason;
    double log10_mass = std::numeric_limits<double>::quiet_NaN();
    double log10_slack = std::numeric_limits<double>::quiet_NaN();
};

double log10_abs_mpz(const mpz_class& value) {
    if (sgn(value) == 0) return -std::numeric_limits<double>::infinity();
    mpz_class abs_value = value >= 0 ? value : -value;
    long exponent = 0;
    const double mantissa = mpz_get_d_2exp(&exponent, abs_value.get_mpz_t());
    return std::log10(mantissa) + static_cast<double>(exponent) * LOG10_2_DOUBLE;
}

double log10_abs_mpq(const mpq_class& value) {
    if (sgn(value) == 0) return -std::numeric_limits<double>::infinity();
    return log10_abs_mpz(value.get_num()) - log10_abs_mpz(value.get_den());
}

std::vector<mpz_class> stable_moments(int max_index) {
    std::vector<mpz_class> moments(static_cast<std::size_t>(max_index + 1));
    moments[0] = 1;
    if (max_index >= 1) moments[1] = 0;
    for (int index = 1; index < max_index; ++index) {
        mpz_class value = index * moments[static_cast<std::size_t>(index)]
            + index * moments[static_cast<std::size_t>(index - 1)];
        if (index >= 2) {
            value -= (mpz_class(index) * (index - 1) / 2)
                * moments[static_cast<std::size_t>(index - 2)];
        }
        if (sgn(value) <= 0) {
            std::cerr << "stable moment recurrence lost positivity at index "
                      << (index + 1) << "\n";
            std::exit(1);
        }
        moments[static_cast<std::size_t>(index + 1)] = std::move(value);
    }
    return moments;
}

std::vector<mpz_class> binom_row(int n_value) {
    std::vector<mpz_class> row(static_cast<std::size_t>(n_value + 1));
    row[0] = 1;
    for (int k = 0; k < n_value; ++k) {
        row[static_cast<std::size_t>(k + 1)] =
            row[static_cast<std::size_t>(k)] * (n_value - k) / (k + 1);
    }
    return row;
}

mpz_class q3_integer(const std::vector<mpz_class>& moments, int n_value) {
    const std::vector<mpz_class> row = binom_row(n_value);
    mpz_class total = 0;
    for (int k = 0; k <= n_value; ++k) {
        total += row[static_cast<std::size_t>(k)]
            * (
                moments[static_cast<std::size_t>(k + 2)]
                    * moments[static_cast<std::size_t>(n_value - k)]
                - moments[static_cast<std::size_t>(k + 1)]
                    * moments[static_cast<std::size_t>(n_value - k + 1)]
            );
    }
    return 2 * total;
}

Poly trim(Poly poly) {
    while (poly.size() > 1 && sgn(poly.back()) == 0) poly.pop_back();
    return poly;
}

Poly poly_add(const Poly& first, const Poly& second) {
    Poly out(std::max(first.size(), second.size()));
    for (std::size_t index = 0; index < first.size(); ++index) out[index] += first[index];
    for (std::size_t index = 0; index < second.size(); ++index) out[index] += second[index];
    return trim(std::move(out));
}

Poly poly_scale(const Poly& poly, const mpq_class& factor) {
    Poly out(poly.size());
    for (std::size_t index = 0; index < poly.size(); ++index) out[index] = poly[index] * factor;
    return trim(std::move(out));
}

Poly poly_mul(const Poly& first, const Poly& second) {
    Poly out(first.size() + second.size() - 1);
    for (std::size_t i = 0; i < first.size(); ++i) {
        for (std::size_t j = 0; j < second.size(); ++j) {
            out[i + j] += first[i] * second[j];
        }
    }
    return trim(std::move(out));
}

Poly poly_shift_x(const Poly& poly) {
    Poly out(poly.size() + 1);
    for (std::size_t index = 0; index < poly.size(); ++index) out[index + 1] = poly[index];
    return trim(std::move(out));
}

Poly chebyshev_t(int degree) {
    if (degree == 0) return {mpq_class(1)};
    if (degree == 1) return {mpq_class(0), mpq_class(1)};
    Poly previous{mpq_class(1)};
    Poly current{mpq_class(0), mpq_class(1)};
    for (int index = 2; index <= degree; ++index) {
        Poly next = poly_add(
            poly_scale(poly_shift_x(current), mpq_class(2)),
            poly_scale(previous, mpq_class(-1))
        );
        previous = std::move(current);
        current = std::move(next);
    }
    return current;
}

Poly compose(const Poly& poly, const Poly& argument) {
    Poly out{mpq_class(0)};
    Poly argument_power{mpq_class(1)};
    for (const mpq_class& coefficient : poly) {
        out = poly_add(out, poly_scale(argument_power, coefficient));
        argument_power = poly_mul(argument_power, argument);
    }
    return trim(std::move(out));
}

mpq_class eval_poly(const Poly& poly, const mpq_class& value) {
    mpq_class out = 0;
    mpq_class value_power = 1;
    for (const mpq_class& coefficient : poly) {
        out += coefficient * value_power;
        value_power *= value;
    }
    return out;
}

std::vector<Row> make_rows(
    bool include_b,
    bool include_c,
    int b_rank_lo,
    int b_rank_hi,
    int c_rank_lo,
    int c_rank_hi
) {
    std::vector<Row> rows;
    if (include_b) {
        for (int rank = b_rank_lo; rank <= b_rank_hi; ++rank) rows.push_back({'B', rank});
    }
    if (include_c) {
        for (int rank = c_rank_lo; rank <= c_rank_hi; ++rank) rows.push_back({'C', rank});
    }
    return rows;
}

int first_boundary(const Row& row) {
    return 2 * row.rank + 2;
}

int dimension_for(const Row& row) {
    return row.rank * (2 * row.rank + 1);
}

DeltaMap parse_bc_delta_logs(const std::vector<std::string>& paths) {
    DeltaMap deltas;
    for (const std::string& path : paths) {
        std::ifstream input(path);
        if (!input) {
            std::cerr << "could not open B/C delta log: " << path << "\n";
            std::exit(2);
        }
        std::string line;
        while (std::getline(input, line)) {
            std::istringstream in(line);
            std::string label;
            std::string delta_label;
            std::string equals;
            mpz_class value;
            in >> label >> delta_label >> equals >> value;
            if (!in) continue;
            if (label.size() < 3 || !(label[0] == 'B' || label[0] == 'C') || label[1] != '_') {
                continue;
            }
            if (delta_label.rfind("Delta_", 0) != 0) continue;
            const char family = label[0];
            const int rank = std::atoi(label.c_str() + 2);
            const int moment = std::atoi(delta_label.c_str() + 6);
            deltas[{family, rank}][moment] = value;
        }
    }
    return deltas;
}

RowResult verify_row(
    const Row& row,
    const std::vector<mpz_class>& stable,
    const Poly& cheb_poly,
    const DeltaMap& deltas,
    const mpq_class& limit
) {
    RowResult result;
    result.family = row.family;
    result.rank = row.rank;

    std::vector<mpz_class> moments = stable;
    const int boundary = first_boundary(row);
    if (boundary < static_cast<int>(moments.size())) {
        const auto found = deltas.find({row.family, row.rank});
        if (found == deltas.end()) {
            result.reason = "missing exact low-boundary deltas";
            return result;
        }
        for (int index = boundary; index < static_cast<int>(moments.size()); ++index) {
            const auto delta_found = found->second.find(index);
            if (delta_found == found->second.end()) {
                result.reason = "missing exact delta below Chebyshev moment ceiling";
                return result;
            }
            if (sgn(delta_found->second) > 0) {
                result.reason = "positive B/C correction";
                return result;
            }
            moments[static_cast<std::size_t>(index)] += delta_found->second;
        }
    }

    const int dimension = dimension_for(row);
    const mpq_class upper_support = 2 * dimension;
    const mpq_class scale = 12000;
    const mpq_class negative_cutoff = 80;
    const mpq_class q0 = negative_cutoff * (negative_cutoff + upper_support);
    const Poly argument{
        mpq_class(1),
        -2 * upper_support / scale,
        mpq_class(2) / scale,
    };
    const Poly composed = compose(cheb_poly, argument);
    const Poly tail_majorant = poly_mul(composed, composed);
    if (tail_majorant.size() + 1 > stable.size()) {
        result.reason = "tail majorant exceeds available moment ceiling";
        return result;
    }
    const mpq_class denominator = eval_poly(
        cheb_poly,
        1 + 2 * q0 / scale
    );
    const mpq_class denominator_sq = denominator * denominator;
    if (sgn(denominator_sq) <= 0) {
        result.reason = "nonpositive Chebyshev denominator";
        return result;
    }

    mpq_class numerator = 0;
    for (std::size_t index = 0; index < tail_majorant.size(); ++index) {
        numerator += tail_majorant[index] * q3_integer(moments, static_cast<int>(index));
    }
    const mpq_class mass = numerator / denominator_sq;
    if (sgn(mass) <= 0) {
        result.reason = "nonpositive negative-tail mass upper bound";
        return result;
    }
    result.log10_mass = log10_abs_mpq(mass);
    result.log10_slack = log10_abs_mpq(limit - mass);
    if (mass > limit) {
        result.reason = "mass upper bound exceeds 1e-15";
        return result;
    }
    result.ok = 1;
    return result;
}

void print_usage(const char* argv0) {
    std::cerr
        << "usage: " << argv0
        << " --bc-delta-log PATH [--bc-delta-log PATH ...]"
        << " [--b-rank-lo N] [--b-rank-hi N] [--c-rank-lo N] [--c-rank-hi N]"
        << " [--only-b] [--only-c] [--progress]\n";
}

}  // namespace

int main(int argc, char** argv) {
    int b_rank_lo = 14;
    int b_rank_hi = 217;
    int c_rank_lo = 20;
    int c_rank_hi = 217;
    bool include_b = true;
    bool include_c = true;
    bool progress = false;
    std::vector<std::string> delta_logs;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--bc-delta-log" && i + 1 < argc) {
            delta_logs.push_back(argv[++i]);
        } else if (arg == "--b-rank-lo" && i + 1 < argc) {
            b_rank_lo = std::atoi(argv[++i]);
        } else if (arg == "--b-rank-hi" && i + 1 < argc) {
            b_rank_hi = std::atoi(argv[++i]);
        } else if (arg == "--c-rank-lo" && i + 1 < argc) {
            c_rank_lo = std::atoi(argv[++i]);
        } else if (arg == "--c-rank-hi" && i + 1 < argc) {
            c_rank_hi = std::atoi(argv[++i]);
        } else if (arg == "--only-b") {
            include_b = true;
            include_c = false;
        } else if (arg == "--only-c") {
            include_b = false;
            include_c = true;
        } else if (arg == "--progress") {
            progress = true;
        } else {
            print_usage(argv[0]);
            return 2;
        }
    }

    if (delta_logs.empty()
        || (!include_b && !include_c)
        || b_rank_lo < 1
        || c_rank_lo < 1
        || b_rank_hi < b_rank_lo
        || c_rank_hi < c_rank_lo) {
        print_usage(argv[0]);
        return 2;
    }

    const std::vector<Row> rows =
        make_rows(include_b, include_c, b_rank_lo, b_rank_hi, c_rank_lo, c_rank_hi);
    const DeltaMap deltas = parse_bc_delta_logs(delta_logs);
    const int cheb_degree = 8;
    const int moment_max = 4 * cheb_degree + 2;
    const std::vector<mpz_class> stable = stable_moments(moment_max);
    const Poly cheb_poly = chebyshev_t(cheb_degree);
    const mpq_class limit(mpz_class(1), mpz_class("1000000000000000"));

    std::cout << "B/C degree-8 cutoff-80 Chebyshev negative-tail GMP verifier\n"
              << "b_rank_lo=" << b_rank_lo << " b_rank_hi=" << b_rank_hi
              << " c_rank_lo=" << c_rank_lo << " c_rank_hi=" << c_rank_hi
              << " include_b=" << include_b << " include_c=" << include_c
              << " cheb_degree=8 cheb_scale=12000 cutoff=80"
              << " target_bound=1e-15"
              << " moment_max=" << moment_max
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;
    std::cout << "delta logs:" << std::endl;
    for (const std::string& path : delta_logs) std::cout << "  " << path << std::endl;

    std::vector<RowResult> results(rows.size());
    int completed = 0;
#pragma omp parallel for schedule(dynamic, 1)
    for (int index = 0; index < static_cast<int>(rows.size()); ++index) {
        const Row row = rows[static_cast<std::size_t>(index)];
        RowResult result = verify_row(
            row,
            stable,
            cheb_poly,
            deltas,
            limit
        );
        results[static_cast<std::size_t>(index)] = result;

        int done_now = 0;
#pragma omp atomic capture
        done_now = ++completed;
        if (progress) {
#pragma omp critical
            {
                std::cout << "row_done " << result.family << "_" << result.rank
                          << " ok=" << result.ok
                          << " log10_mass_upper=" << std::setprecision(18)
                          << result.log10_mass
                          << " completed=" << done_now << "/" << rows.size();
                if (!result.ok) std::cout << " reason=\"" << result.reason << "\"";
                std::cout << std::endl;
                std::fflush(stdout);
            }
        }
    }

    int failures = 0;
    bool have_worst = false;
    RowResult worst;
    for (const RowResult& result : results) {
        if (!result.ok) {
            ++failures;
            continue;
        }
        if (!have_worst || result.log10_mass > worst.log10_mass) {
            worst = result;
            have_worst = true;
        }
    }

    std::cout << std::setprecision(18)
              << "row\tlog10_mass_upper\tlog10_exact_slack_to_1e-15\n";
    for (const RowResult& result : results) {
        std::cout << result.family << "_" << result.rank
                  << '\t' << result.log10_mass
                  << '\t' << result.log10_slack;
        if (!result.ok) std::cout << "\tFAIL\t" << result.reason;
        std::cout << '\n';
    }

    std::cout << "SUMMARY failures=" << failures;
    if (have_worst) {
        std::cout << " worst_row=" << worst.family << "_" << worst.rank
                  << " worst_log10_mass_upper=" << worst.log10_mass
                  << " worst_log10_exact_slack_to_1e-15=" << worst.log10_slack;
    }
    std::cout << std::endl;
    return failures == 0 ? 0 : 1;
}
