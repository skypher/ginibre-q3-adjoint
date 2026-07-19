// Exact-arithmetic GMP shift-2 log-convexity verification.
//
// Reads a moment log, verifies m_n * m_{n+4} > m_{n+2}^2 for every n in [4, n_max-4]
// using big-integer arithmetic.  Prints:
//   - Number of odd/even n checked.
//   - Any violations (should be zero).
//   - Final Chain-Inequality differences Q_3(n+2) - 4 Q_3(n) > 0 if moments
//     suffice.
//
// Build:
//   g++ -O3 -std=c++20 verify_shift2_cpp.cpp -lgmpxx -lgmp -o verify_shift2_cpp
// Run:
//   ./verify_shift2_cpp LOG_PATH

#include <gmpxx.h>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

struct MomentEntry { int n; mpz_class m; };

std::vector<MomentEntry> parse_log(const std::string& path) {
    std::vector<MomentEntry> out;
    std::ifstream f(path);
    if (!f) {
        std::cerr << "Cannot open " << path << std::endl;
        return out;
    }
    std::string line;
    std::regex re(R"(^m_(\d+) = (\d+))");
    while (std::getline(f, line)) {
        std::smatch sm;
        if (std::regex_search(line, sm, re)) {
            int n = std::stoi(sm[1].str());
            mpz_class v(sm[2].str());
            out.push_back({n, std::move(v)});
        }
    }
    std::sort(out.begin(), out.end(),
              [](const MomentEntry& a, const MomentEntry& b){ return a.n < b.n; });
    return out;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " LOG_PATH [N_START] [N_STOP]" << std::endl;
        return 1;
    }
    std::string path = argv[1];
    auto moms = parse_log(path);
    std::cout << "Loaded " << moms.size() << " moments from " << path << std::endl;
    if (moms.empty()) return 1;

    // Index by n for fast lookup
    int n_max = moms.back().n;
    std::vector<mpz_class> m_of(n_max + 5);  // index n -> m_n; missing = 0
    std::vector<bool> has(n_max + 5, false);
    for (auto& me : moms) {
        m_of[me.n] = me.m;
        has[me.n] = true;
    }

    int n_start = (argc > 2) ? std::atoi(argv[2]) : 4;
    int n_stop = (argc > 3) ? std::atoi(argv[3]) : (n_max - 4);

    std::cout << "Shift-2 log-convexity check: n in [" << n_start << ", " << n_stop << "]" << std::endl;

    int even_checks = 0, even_violations = 0;
    int odd_checks = 0, odd_violations = 0;
    mpz_class lhs, rhs, diff;
    std::vector<int> violation_ns;

    for (int n = n_start; n <= n_stop; ++n) {
        if (!has[n] || !has[n+2] || !has[n+4]) continue;
        lhs = m_of[n] * m_of[n+4];
        rhs = m_of[n+2] * m_of[n+2];
        diff = lhs - rhs;
        bool ok = diff > 0;
        if (n % 2 == 0) {
            ++even_checks;
            if (!ok) { ++even_violations; violation_ns.push_back(n); }
        } else {
            ++odd_checks;
            if (!ok) { ++odd_violations; violation_ns.push_back(n); }
        }
    }

    std::cout << "\nSummary:" << std::endl;
    std::cout << "  Even n checks: " << even_checks << ", violations: " << even_violations << std::endl;
    std::cout << "  Odd  n checks: " << odd_checks << ", violations: " << odd_violations << std::endl;
    std::cout << "  Total checks:  " << (even_checks + odd_checks) << ", total violations: " << (even_violations + odd_violations) << std::endl;

    if (!violation_ns.empty()) {
        std::cout << "\nFirst 10 violations (n values):";
        for (size_t i = 0; i < std::min((size_t)10, violation_ns.size()); ++i) {
            std::cout << " " << violation_ns[i];
        }
        std::cout << std::endl;
    } else {
        std::cout << "\nAll " << (even_checks + odd_checks) << " shift-2 log-convexity checks PASSED."
                  << std::endl;
    }

    return (even_violations + odd_violations) == 0 ? 0 : 2;
}
