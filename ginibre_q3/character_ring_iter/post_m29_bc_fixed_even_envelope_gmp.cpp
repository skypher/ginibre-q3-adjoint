#include <gmpxx.h>
#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <unistd.h>

namespace {

using Cell = std::pair<int, int>;
using Part = std::vector<int>;

std::string part_key(const Part& part) {
    std::ostringstream out;
    for (std::size_t i = 0; i < part.size(); ++i) {
        if (i) out << ',';
        out << part[i];
    }
    return out.str();
}

std::string part_display(const Part& part) {
    std::ostringstream out;
    out << '(';
    for (std::size_t i = 0; i < part.size(); ++i) {
        if (i) out << ',';
        out << part[i];
    }
    out << ')';
    return out.str();
}

void gen_parts_rec(int n, int max_part, Part& current, std::vector<Part>& out) {
    if (n == 0) {
        out.push_back(current);
        return;
    }
    for (int x = std::min(n, max_part); x >= 1; --x) {
        current.push_back(x);
        gen_parts_rec(n - x, x, current, out);
        current.pop_back();
    }
}

std::vector<Part> gen_parts(int n) {
    std::vector<Part> out;
    Part current;
    gen_parts_rec(n, n, current, out);
    return out;
}

bool has_cell(const Part& part, Cell cell) {
    const int r = cell.first;
    const int c = cell.second;
    return 1 <= r && r <= static_cast<int>(part.size()) &&
           1 <= c && c <= part[static_cast<std::size_t>(r - 1)];
}

std::vector<Cell> cells_of(const Part& part) {
    std::vector<Cell> out;
    for (int r = 1; r <= static_cast<int>(part.size()); ++r) {
        for (int c = 1; c <= part[static_cast<std::size_t>(r - 1)]; ++c) {
            out.push_back({r, c});
        }
    }
    return out;
}

bool is_removable(const Part& part, Cell cell) {
    if (!has_cell(part, cell)) return false;
    const int r = cell.first;
    const int c = cell.second;
    if (part[static_cast<std::size_t>(r - 1)] != c) return false;
    if (r < static_cast<int>(part.size()) &&
        part[static_cast<std::size_t>(r)] >= c) {
        return false;
    }
    return true;
}

std::vector<Cell> removable_corners(const Part& part) {
    std::vector<Cell> out;
    for (int r = 1; r <= static_cast<int>(part.size()); ++r) {
        Cell cell{r, part[static_cast<std::size_t>(r - 1)]};
        if (is_removable(part, cell)) out.push_back(cell);
    }
    return out;
}

Part remove_cell(Part part, Cell cell) {
    --part[static_cast<std::size_t>(cell.first - 1)];
    while (!part.empty() && part.back() == 0) part.pop_back();
    return part;
}

mpz_class pow2(int exponent) {
    mpz_class out;
    mpz_ui_pow_ui(out.get_mpz_t(), 2UL, static_cast<unsigned long>(exponent));
    return out;
}

int bit_length(const mpz_class& value) {
    if (sgn(value) <= 0) return 0;
    return static_cast<int>(mpz_sizeinbase(value.get_mpz_t(), 2));
}

}  // namespace

int main(int argc, char** argv) {
    int max_m = 23;
    if (argc >= 2) max_m = std::atoi(argv[1]);
    if (max_m < 1) {
        std::cerr << "max_m must be positive\n";
        std::cout << "__EXIT_STATUS=2\n";
        return 2;
    }

    char hostname[256] = {};
    if (gethostname(hostname, sizeof(hostname) - 1) != 0) hostname[0] = '\0';

    std::cout << "B/C fixed-even recursive-envelope GMP verifier\n"
              << "VERIFIES W(gamma)=max_e sum_{y in Y(gamma,e)} W(gamma\\\\{y,e})\n"
              << "and max_{|gamma|=2m-1} W(gamma) <= 2^(2m) for reported m\n"
              << "host=" << (hostname[0] == '\0' ? "unknown" : hostname) << "\n"
              << "max_m=" << max_m
              << " OpenMP threads=" << omp_get_max_threads() << std::endl;
    std::fflush(stdout);

    std::unordered_map<std::string, mpz_class> envelope;
    int failures = 0;

    for (int m = 1; m <= max_m; ++m) {
        const int n = 2 * m - 1;
        std::vector<Part> parts = gen_parts(n);
        std::vector<mpz_class> values(parts.size());

#pragma omp parallel for schedule(dynamic, 128)
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(parts.size()); ++i) {
            const Part& part = parts[static_cast<std::size_t>(i)];
            if (n == 1) {
                values[static_cast<std::size_t>(i)] = 1;
                continue;
            }

            mpz_class best = 0;
            const std::vector<Cell> cells = cells_of(part);
            const std::vector<Cell> corners = removable_corners(part);
            for (Cell e : cells) {
                mpz_class total = 0;
                for (Cell y : corners) {
                    if (y == e) continue;
                    Part after_y = remove_cell(part, y);
                    if (!is_removable(after_y, e)) continue;
                    Part child = remove_cell(after_y, e);
                    const auto it = envelope.find(part_key(child));
                    if (it == envelope.end()) {
                        std::cerr << "missing child envelope for " << part_display(child) << "\n";
                        std::exit(3);
                    }
                    total += it->second;
                }
                if (total > best) best = total;
            }
            values[static_cast<std::size_t>(i)] = best;
        }

        mpz_class max_w = 0;
        Part best_part;
        for (std::size_t i = 0; i < parts.size(); ++i) {
            const mpz_class& value = values[i];
            envelope[part_key(parts[i])] = value;
            if (value > max_w) {
                max_w = value;
                best_part = parts[i];
            }
        }

        const mpz_class budget = pow2(2 * m);
        const bool pass = max_w <= budget;
        if (!pass) ++failures;
        std::cout << "RESULT m=" << m
                  << " partitions=" << parts.size()
                  << " max_W=" << max_w.get_str()
                  << " max_W_bits=" << bit_length(max_w)
                  << " budget=2^" << (2 * m)
                  << " budget_value=" << budget.get_str()
                  << " pass_full_budget=" << (pass ? "yes" : "no")
                  << " best_shape=" << part_display(best_part)
                  << std::endl;
        std::fflush(stdout);
    }

    std::cout << "SUMMARY failures=" << failures << std::endl;
    std::cout << "__EXIT_STATUS=" << (failures == 0 ? 0 : 1) << std::endl;
    return failures == 0 ? 0 : 1;
}
