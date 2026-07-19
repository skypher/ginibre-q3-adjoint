#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr long double LOG2 = 0.693147180559945309417232121458176568L;
constexpr long double LOG_HALF = -0.693147180559945309417232121458176568L;

struct Params {
    long double r;
    long double eta;
    long double a;
};

struct Margins {
    long double gaussian;
    long double r0;
    long double r1;
    long double r2;
    long double gaussian_derivative;
    long double r0_derivative;
    long double r1_derivative;
    long double r2_derivative;
};

struct Candidate {
    int cutoff;
    Params params;
    Margins margins;
};

long double tau_h(int c_value) {
    int n_value = c_value / 2;
    if (n_value < 2) n_value = 2;
    long double n = static_cast<long double>(n_value);
    return std::pow(n, -2.5L) * std::exp(-0.5L * n);
}

Margins margins_at(const Params& p, int c_value) {
    const long double c = static_cast<long double>(c_value);
    const long double l_value = 2.0L * p.r + 3.0L * p.eta;
    const long double h = tau_h(c_value);
    const long double den = 1.0L - 1.0L / (p.eta * p.eta * c * c);
    const long double s_log =
        LOG2 - (p.eta * c - 1.0L) * (p.eta * c - 1.0L) /
                   (16.0L * (1.0L + h));

    Margins m;
    m.gaussian = (p.a - l_value / 2.0L) * c
        + 0.5L * std::log(l_value * c)
        - std::log(l_value * c + 1.0L)
        - 0.5L * LOG2
        - LOG2;
    m.r0 = p.a * c + s_log - LOG_HALF;
    m.r1 = std::log(144.0L * c * c * std::pow(p.r / 2.0L, 3.0L)
                     / (p.r * p.r * den))
        + p.a * c + s_log + c * std::log(2.0L / p.r) - LOG_HALF;
    m.r2 = std::log(4.0L * std::pow(p.r / p.eta, 3.0L)
                     / (p.r * p.r * c * c * den))
        - (std::log(p.r / p.eta) - p.a) * c - LOG_HALF;

    const long double tau_derivative_upper =
        -p.eta * (p.eta * c - 1.0L) / (8.0L * (1.0L + h));
    m.gaussian_derivative = p.a - l_value / 2.0L - 1.0L / (2.0L * c);
    m.r0_derivative = p.a + tau_derivative_upper;
    m.r1_derivative = 2.0L / c + p.a + tau_derivative_upper
        + std::log(2.0L / p.r);
    m.r2_derivative = p.a - std::log(p.r / p.eta);
    return m;
}

bool admissible(const Params& p) {
    if (!(p.r > 2.0L)) return false;
    if (!(p.eta > 0.0L && p.eta < 0.2L)) return false;
    const long double l_value = 2.0L * p.r + 3.0L * p.eta;
    if (!(p.a > l_value / 2.0L)) return false;
    if (!(p.a < std::log(p.r / p.eta))) return false;
    return true;
}

bool certificate_ok(const Margins& m) {
    return m.gaussian >= 0.0L
        && m.r0 <= 0.0L
        && m.r1 <= 0.0L
        && m.r2 <= 0.0L
        && m.gaussian_derivative > 0.0L
        && m.r0_derivative < 0.0L
        && m.r1_derivative < 0.0L
        && m.r2_derivative < 0.0L;
}

bool derivative_ok(const Margins& m) {
    return m.gaussian_derivative > 0.0L
        && m.r0_derivative < 0.0L
        && m.r1_derivative < 0.0L
        && m.r2_derivative < 0.0L;
}

bool margin_ok(const Margins& m) {
    return m.gaussian >= 0.0L
        && m.r0 <= 0.0L
        && m.r1 <= 0.0L
        && m.r2 <= 0.0L;
}

Candidate threshold(const Params& p, int max_cutoff) {
    if (!admissible(p)) return {0, p, {}};
    int start = std::max(3, static_cast<int>(std::floor(1.0L / p.eta)) + 1);
    const long double l_value = 2.0L * p.r + 3.0L * p.eta;
    const long double gaussian_gap = p.a - l_value / 2.0L;
    if (gaussian_gap <= 0.0L) return {0, p, {}};

    int derivative_start = start;
    derivative_start = std::max(
        derivative_start,
        static_cast<int>(std::floor(1.0L / (2.0L * gaussian_gap))) + 1);
    derivative_start = std::max(
        derivative_start,
        static_cast<int>(std::floor((8.0L * p.a / p.eta + 1.0L) / p.eta)) - 50);
    derivative_start = std::max(
        derivative_start,
        static_cast<int>(
            std::floor((8.0L * (p.a + std::log(2.0L / p.r)) / p.eta + 1.0L)
                       / p.eta))
            - 50);
    derivative_start = std::max(start, derivative_start);

    int lo = 0;
    for (int c = derivative_start; c <= max_cutoff; ++c) {
        if (derivative_ok(margins_at(p, c))) {
            lo = c;
            break;
        }
    }
    if (lo == 0) return {0, p, {}};
    if (!margin_ok(margins_at(p, max_cutoff))) return {0, p, {}};

    int left = lo;
    int right = max_cutoff;
    while (left < right) {
        int mid = left + (right - left) / 2;
        Margins m = margins_at(p, mid);
        if (margin_ok(m)) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }
    Margins m = margins_at(p, left);
    if (certificate_ok(m)) return {left, p, m};
    return {0, p, {}};
}

void print_candidate(const Candidate& candidate) {
    const Params& p = candidate.params;
    const long double l_value = 2.0L * p.r + 3.0L * p.eta;
    std::cout << std::setprecision(21)
              << "C=" << candidate.cutoff
              << " r=" << p.r
              << " eta=" << p.eta
              << " A=" << p.a
              << " L=" << l_value
              << " log(r/eta)-A=" << (std::log(p.r / p.eta) - p.a)
              << "\n";
    const Margins& m = candidate.margins;
    std::cout << "  margins"
              << " gaussian=" << m.gaussian
              << " R0=" << m.r0
              << " R1=" << m.r1
              << " R2=" << m.r2
              << "\n";
    std::cout << "  derivatives"
              << " gaussian=" << m.gaussian_derivative
              << " R0=" << m.r0_derivative
              << " R1=" << m.r1_derivative
              << " R2=" << m.r2_derivative
              << "\n";
}

long double parse_long_double(const char* text) {
    char* end = nullptr;
    long double value = std::strtold(text, &end);
    if (end == text || *end != '\0') {
        throw std::runtime_error(std::string("invalid number: ") + text);
    }
    return value;
}

int parse_int(const char* text) {
    char* end = nullptr;
    long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0') {
        throw std::runtime_error(std::string("invalid integer: ") + text);
    }
    return static_cast<int>(value);
}

void run_current() {
    Params p{
        20000001.0L / 10000000.0L,
        9973.0L / 50000.0L,
        461.0L / 200.0L,
    };
    Candidate candidate = threshold(p, 2000);
    if (candidate.cutoff == 0) throw std::runtime_error("current parameters did not certify");
    print_candidate(candidate);
}

void run_certificate(int argc, char** argv) {
    Params p{
        20000001.0L / 10000000.0L,
        9973.0L / 50000.0L,
        461.0L / 200.0L,
    };
    int cutoff = 945;
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--r" && i + 1 < argc) {
            p.r = parse_long_double(argv[++i]);
        } else if (arg == "--eta" && i + 1 < argc) {
            p.eta = parse_long_double(argv[++i]);
        } else if (arg == "--A" && i + 1 < argc) {
            p.a = parse_long_double(argv[++i]);
        } else if (arg == "--cutoff" && i + 1 < argc) {
            cutoff = parse_int(argv[++i]);
        } else {
            throw std::runtime_error("unknown certificate argument " + arg);
        }
    }
    if (!admissible(p)) throw std::runtime_error("parameters are not admissible");
    Candidate candidate{cutoff, p, margins_at(p, cutoff)};
    print_candidate(candidate);
    if (!certificate_ok(candidate.margins)) {
        throw std::runtime_error("cutoff certificate failed");
    }
}

void run_search(int argc, char** argv) {
    int max_cutoff = 1600;
    int keep = 20;
    int eta_steps = 220;
    int r_steps = 280;
    int a_steps = 120;
    long double eta_lo = 0.185L;
    long double eta_hi = 0.19999L;
    long double r_lo = 2.0000001L;
    long double r_hi = 2.06L;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--max-cutoff" && i + 1 < argc) {
            max_cutoff = parse_int(argv[++i]);
        } else if (arg == "--keep" && i + 1 < argc) {
            keep = parse_int(argv[++i]);
        } else if (arg == "--eta-steps" && i + 1 < argc) {
            eta_steps = parse_int(argv[++i]);
        } else if (arg == "--r-steps" && i + 1 < argc) {
            r_steps = parse_int(argv[++i]);
        } else if (arg == "--A-steps" && i + 1 < argc) {
            a_steps = parse_int(argv[++i]);
        } else if (arg == "--eta-lo" && i + 1 < argc) {
            eta_lo = parse_long_double(argv[++i]);
        } else if (arg == "--eta-hi" && i + 1 < argc) {
            eta_hi = parse_long_double(argv[++i]);
        } else if (arg == "--r-lo" && i + 1 < argc) {
            r_lo = parse_long_double(argv[++i]);
        } else if (arg == "--r-hi" && i + 1 < argc) {
            r_hi = parse_long_double(argv[++i]);
        } else {
            throw std::runtime_error("unknown search argument " + arg);
        }
    }

    std::vector<Candidate> best;
    long long tested = 0;
    for (int e = 0; e <= eta_steps; ++e) {
        long double eta = eta_lo + (eta_hi - eta_lo) * e / eta_steps;
        for (int ri = 0; ri <= r_steps; ++ri) {
            long double r = r_lo + (r_hi - r_lo) * ri / r_steps;
            const long double l_value = 2.0L * r + 3.0L * eta;
            const long double a_lo = l_value / 2.0L;
            const long double a_hi = std::log(r / eta);
            if (!(a_hi > a_lo)) continue;
            for (int ai = 1; ai < a_steps; ++ai) {
                long double a = a_lo + (a_hi - a_lo) * ai / a_steps;
                Candidate candidate = threshold({r, eta, a}, max_cutoff);
                ++tested;
                if (candidate.cutoff == 0) continue;
                if (static_cast<int>(best.size()) < keep
                    || candidate.cutoff < best.back().cutoff) {
                    best.push_back(candidate);
                    std::sort(best.begin(), best.end(), [](const Candidate& x, const Candidate& y) {
                        return x.cutoff < y.cutoff;
                    });
                    if (static_cast<int>(best.size()) > keep) best.pop_back();
                }
            }
        }
        if (e % std::max(1, eta_steps / 40) == 0 || e == eta_steps) {
            std::cout << "progress eta_index=" << e << "/" << eta_steps
                      << " tested=" << tested
                      << " best=" << (best.empty() ? 0 : best.front().cutoff)
                      << std::endl;
        }
    }

    std::cout << "best candidates" << std::endl;
    for (const Candidate& candidate : best) print_candidate(candidate);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc <= 1) {
            std::cerr << "usage: classical_trace_cutoff_search "
                      << "--current | --certificate [--r R --eta ETA --A A --cutoff C] "
                      << "| --search [options]" << std::endl;
            return 2;
        }
        std::string mode = argv[1];
        if (mode == "--current") {
            run_current();
        } else if (mode == "--certificate") {
            run_certificate(argc, argv);
        } else if (mode == "--search") {
            run_search(argc, argv);
        } else {
            throw std::runtime_error("unknown mode " + mode);
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << std::endl;
        return 1;
    }
}
