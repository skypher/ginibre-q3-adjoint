// Extract empirical c_3 from character-ring moment log, with c_1 and c_2
// FIXED at their rigorous closed-form values.  Uses GMP/MPFR for
// high-precision arithmetic so the dominant exponential factor
// cancellations don't lose precision.
//
// Model:
//   m_n / tilde_m_n = 1 + c_1/n + c_2/n^2 + c_3/n^3 + c_4/n^4 + ...
//   tilde_m_n = const * d^n * n^{-d/2}
//
// For each large n, compute
//   R(n) := log(m_n) - n log(d) + (d/2) log(n) - log(1 + c_1/n + c_2/n^2)
// then (R(n) - log_const) * n^3 -> c_3 as n -> infinity (approximately).
//
// More robust: fit log(m_n) - n log d + (d/2) log n = log_const + f(1/n)
// where f(x) = log(1 + c_1 x + c_2 x^2 + c_3 x^3 + ...) = sum_k tilde_c_k x^k
// tilde_c_1 = c_1  (fixed)
// tilde_c_2 = c_2 - c_1^2/2  (fixed from c_2)
// tilde_c_3 = c_3 - c_1 c_2 + c_1^3/3  (c_3 unknown)
// tilde_c_4 = c_4 - c_1 c_3 - c_2^2/2 + c_1^2 c_2 - c_1^4/4 (unknown)
// We fit log_const, tilde_c_3, tilde_c_4 (and maybe tilde_c_5, tilde_c_6)
// from large-n data, then back out c_3.
//
// Build:
//   g++ -O3 -std=c++20 c3_fit_cpp.cpp -lgmpxx -lgmp -lmpfr -o c3_fit_cpp
// Run:
//   ./c3_fit_cpp FAMILY  (G2, F4, E6, E7, E8)
//   Reads ~/character_ring_iter/logs/{family_lowercase}_*.log (uses largest).

#include <gmpxx.h>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <regex>
#include <cmath>
#include <optional>

struct Family {
    std::string name;
    int d;
    mpq_class c1;
    mpq_class c2;
    std::optional<mpq_class> c3_rigorous;
};

std::vector<Family> families = {
    {"G2", 14, mpq_class(-63, 4), mpq_class(165851, 1152), mpq_class(mpz_class(-3080875), mpz_class(3072))},
    {"F4", 52, mpq_class(-182, 1), mpq_class(8401705, 486), mpq_class(mpz_class("-44918139383"), mpz_class("39366"))},
    // E_6, E_7: c_2 from multi-seed MC (not yet rigorous-exact)
    {"E6", 78, mpq_class(-1599, 4), mpq_class(83417, 1), std::nullopt},
    {"E7", 133, mpq_class(-18221, 16), mpq_class(679448, 1), std::nullopt},
    {"E8", 248, mpq_class(-3906, 1), mpq_class(346377539, 45), std::nullopt},
};

struct MomentEntry {
    int n;
    mpz_class m;
};

std::vector<MomentEntry> parse_log(const std::string& path) {
    std::vector<MomentEntry> out;
    std::ifstream f(path);
    if (!f) return out;
    std::string line;
    std::regex re(R"(^m_(\d+) = (\d+))");
    while (std::getline(f, line)) {
        std::smatch m;
        if (std::regex_search(line, m, re)) {
            int n = std::stoi(m[1].str());
            mpz_class v(m[2].str());
            if (v > 0) out.push_back({n, v});
        }
    }
    std::sort(out.begin(), out.end(),
              [](const MomentEntry& a, const MomentEntry& b){ return a.n < b.n; });
    return out;
}

std::string find_largest_log(const std::string& family_lower) {
    namespace fs = std::filesystem;
    std::string logdir = std::string(getenv("HOME")) + "/character_ring_iter/logs";
    int best_count = 0;
    std::string best_path;
    std::regex re(family_lower + R"(_(\d+)\.log)");
    if (!fs::exists(logdir)) {
        logdir = std::string(getenv("HOME")) + "/yang4/ginibre_q3/character_ring_iter/logs";
    }
    for (auto& entry : fs::directory_iterator(logdir)) {
        std::string name = entry.path().filename().string();
        std::smatch m;
        if (std::regex_match(name, m, re)) {
            auto mom = parse_log(entry.path().string());
            if ((int)mom.size() > best_count) {
                best_count = mom.size();
                best_path = entry.path().string();
            }
        }
    }
    return best_path;
}

// Compute log(m_n) for big-integer m_n using GMP: convert to mpf_class first
// then log via series.  Simpler: use mpz2 -> log2-like trick.
// Actually just use double, since log(m_n) = log(mantissa * 2^exp) = exp*log(2) + log(mantissa)
static double log_mpz(const mpz_class& x) {
    if (x == 0) return -std::numeric_limits<double>::infinity();
    // mpz_get_d_2exp returns d * 2^e such that 0.5 <= |d| < 1.
    long exp_out;
    double d = mpz_get_d_2exp(&exp_out, x.get_mpz_t());
    return std::log(d) + exp_out * std::log(2.0);
}

int main(int argc, char** argv) {
    std::string family_name = "G2";
    if (argc > 1) family_name = argv[1];

    Family* fam = nullptr;
    for (auto& f : families) {
        if (f.name == family_name) { fam = &f; break; }
    }
    if (!fam) {
        std::cerr << "Unknown family: " << family_name << std::endl;
        return 1;
    }

    std::string fam_lower = family_name;
    for (auto& c : fam_lower) c = std::tolower(c);
    std::string log_path = find_largest_log(fam_lower);
    if (log_path.empty()) {
        std::cerr << "No log found for " << family_name << std::endl;
        return 1;
    }
    std::cout << "Family: " << family_name << ", d = " << fam->d << std::endl;
    std::cout << "Log file: " << log_path << std::endl;

    auto moments = parse_log(log_path);
    std::cout << "Moments loaded: " << moments.size()
              << ", n range: [" << moments.front().n << ", " << moments.back().n << "]"
              << std::endl;

    if (moments.size() < 20) {
        std::cerr << "Too few moments for fit" << std::endl;
        return 1;
    }

    double c1 = mpq_get_d(fam->c1.get_mpq_t());
    double c2 = mpq_get_d(fam->c2.get_mpq_t());
    double tilde_c_2 = c2 - c1*c1/2.0;
    double d = fam->d;
    double half_d = fam->d / 2.0;

    printf("c_1 = %.12f, c_2 = %.12f\n", c1, c2);
    printf("tilde_c_2 = c_2 - c_1^2/2 = %.12f\n", tilde_c_2);

    // For each n in a wider window (all n >= N_MIN), compute
    //   y_n := log(m_n) - n log(d) + (d/2) log(n) - c_1/n - tilde_c_2/n^2
    // (which should equal log_const + tilde_c_3/n^3 + tilde_c_4/n^4 + ...)
    int N_MIN = std::max(20, moments.back().n / 4);
    std::vector<MomentEntry> filtered;
    for (auto& me : moments) {
        if (me.n >= N_MIN) filtered.push_back(me);
    }
    size_t K = filtered.size();
    size_t start = 0;
    std::vector<MomentEntry>& mom_window = filtered;

    std::vector<double> y_vals(K), n_vals(K);
    double log_d = std::log(d);
    for (size_t i = 0; i < K; ++i) {
        auto& me = mom_window[i];
        double n = me.n;
        double log_n = std::log(n);
        double log_m = log_mpz(me.m);
        double y = log_m - n * log_d + half_d * log_n
                   - c1 / n - tilde_c_2 / (n * n);
        y_vals[i] = y;
        n_vals[i] = n;
    }

    // Linear regression: y = log_const + tilde_c_3/n^3 + tilde_c_4/n^4
    // 3-parameter fit
    int P = 3;
    std::vector<std::vector<double>> X(K, std::vector<double>(P));
    for (size_t i = 0; i < K; ++i) {
        double inv_n = 1.0 / n_vals[i];
        X[i][0] = 1.0;
        X[i][1] = std::pow(inv_n, 3);
        X[i][2] = std::pow(inv_n, 4);
    }

    // Form X^T X and X^T y
    std::vector<std::vector<double>> XtX(P, std::vector<double>(P, 0.0));
    std::vector<double> Xty(P, 0.0);
    for (int a = 0; a < P; ++a) {
        for (int b = 0; b < P; ++b) {
            for (size_t i = 0; i < K; ++i) XtX[a][b] += X[i][a] * X[i][b];
        }
        for (size_t i = 0; i < K; ++i) Xty[a] += X[i][a] * y_vals[i];
    }
    // Gaussian elimination for solving XtX * coeffs = Xty
    std::vector<std::vector<double>> aug(P, std::vector<double>(P + 1));
    for (int a = 0; a < P; ++a) {
        for (int b = 0; b < P; ++b) aug[a][b] = XtX[a][b];
        aug[a][P] = Xty[a];
    }
    for (int pivot = 0; pivot < P; ++pivot) {
        int best = pivot;
        for (int r = pivot + 1; r < P; ++r) {
            if (std::abs(aug[r][pivot]) > std::abs(aug[best][pivot])) best = r;
        }
        std::swap(aug[pivot], aug[best]);
        double d_piv = aug[pivot][pivot];
        for (int c = pivot; c < P + 1; ++c) aug[pivot][c] /= d_piv;
        for (int r = 0; r < P; ++r) {
            if (r != pivot) {
                double factor = aug[r][pivot];
                for (int c = pivot; c < P + 1; ++c) aug[r][c] -= factor * aug[pivot][c];
            }
        }
    }
    double log_const = aug[0][P];
    double tilde_c_3 = aug[1][P];
    double tilde_c_4 = aug[2][P];

    // Back out c_3 = tilde_c_3 + c_1 c_2 - c_1^3/3
    double c_3 = tilde_c_3 + c1 * c2 - std::pow(c1, 3) / 3.0;

    std::cout << "\nFit with " << K << " points (n in [" << n_vals.front() << ", " << n_vals.back() << "]):" << std::endl;
    printf("  log_const      = %.6f\n", log_const);
    printf("  tilde_c_3      = %.6f\n", tilde_c_3);
    printf("  tilde_c_4      = %.6f\n", tilde_c_4);
    printf("  c_3 (derived)  = %.6f\n", c_3);

    if (fam->c3_rigorous) {
        double c_3_rig = mpq_get_d(fam->c3_rigorous->get_mpq_t());
        printf("  c_3 (rigorous) = %.6f\n", c_3_rig);
        printf("  diff (fit - rig) = %.6f (relative %.4e)\n",
               c_3 - c_3_rig,
               (c_3 - c_3_rig) / c_3_rig);
    }

    // Reduced-threshold analysis
    double c3_abs = std::abs(c_3);
    double threshold_c3 = std::cbrt(8.0 * c3_abs);
    double threshold_2c1 = 2.0 * std::abs(c1);
    int n_max = moments.back().n;
    printf("\nReduced-threshold (8|c_3|)^(1/3) = %.2f\n", threshold_c3);
    printf("First-order threshold 2|c_1|     = %.2f\n", threshold_2c1);
    printf("Character-ring n_max             = %d\n", n_max);
    printf("Gap to close: %.2f\n", std::max(0.0, threshold_c3 - n_max));

    return 0;
}
