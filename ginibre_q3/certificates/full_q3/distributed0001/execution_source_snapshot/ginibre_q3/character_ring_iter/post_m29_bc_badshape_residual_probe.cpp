#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

using boost::multiprecision::cpp_int;

static std::vector<int> decode_partition(const std::string& key) {
    std::vector<int> out;
    out.reserve(key.size());
    for (unsigned char c : key) out.push_back(static_cast<int>(c));
    return out;
}

static std::string encode_partition(const std::vector<int>& p) {
    std::string key;
    key.reserve(p.size());
    for (int v : p) key.push_back(static_cast<char>(v));
    return key;
}

static bool is_partition(const std::vector<int>& p) {
    for (size_t i = 1; i < p.size(); ++i) {
        if (p[i - 1] < p[i]) return false;
    }
    return true;
}

static std::vector<std::string> horizontal_two_successors(const std::string& key) {
    const std::vector<int> lam = decode_partition(key);
    const int L = static_cast<int>(lam.size());
    const int max_rows = L + 2;
    std::unordered_set<std::string> seen;

    for (int r1 = 0; r1 < max_rows; ++r1) {
        for (int r2 = r1; r2 < max_rows; ++r2) {
            const int old1 = (r1 < L ? lam[r1] : 0);
            const int old2 = (r2 < L ? lam[r2] : 0);
            if (r1 != r2 && old1 + 1 == old2 + 1) continue;

            std::vector<int> mu(max_rows, 0);
            for (int i = 0; i < L; ++i) mu[i] = lam[i];
            mu[r1] += 1;
            mu[r2] += 1;
            while (!mu.empty() && mu.back() == 0) mu.pop_back();
            if (!is_partition(mu)) continue;
            seen.insert(encode_partition(mu));
        }
    }

    std::vector<std::string> out;
    out.reserve(seen.size());
    for (const auto& s : seen) out.push_back(s);
    return out;
}

static bool conjugate_even(const std::string& key) {
    const size_t L = key.size();
    const size_t E = L + (L & 1U);
    for (size_t i = 0; i < E; i += 2) {
        const int a = (i < L ? static_cast<unsigned char>(key[i]) : 0);
        const int b = (i + 1 < L ? static_cast<unsigned char>(key[i + 1]) : 0);
        if (a != b) return false;
    }
    return true;
}

static bool bad_width(const std::string& key, int q) {
    return !key.empty() && static_cast<unsigned char>(key[0]) >= 2 * q;
}

static bool bad_height(const std::string& key, int q) {
    return static_cast<int>(key.size()) >= 2 * q;
}

static cpp_int binom_int(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    cpp_int ans = 1;
    for (int i = 1; i <= k; ++i) {
        ans *= (n - k + i);
        ans /= i;
    }
    return ans;
}

static std::string decimal_prefix(const cpp_int& x, size_t digits = 48) {
    const std::string s = x.convert_to<std::string>();
    if (s.size() <= digits) return s;
    return s.substr(0, digits) + "...(" + std::to_string(s.size()) + " digits)";
}

int main(int argc, char** argv) {
    int max_j = 30;
    int probe_q = 15;
    if (argc >= 2) max_j = std::atoi(argv[1]);
    if (argc >= 3) probe_q = std::atoi(argv[2]);
    if (max_j < 0 || probe_q < 0 || max_j > 60) {
        std::cerr << "usage: post_m29_bc_badshape_residual_probe [max_j<=60] [q]\n";
        return 2;
    }

    const auto start = std::chrono::steady_clock::now();
    std::cout << "B/C bad-shape residual Pieri diagnostic\n";
    std::cout << "max_j=" << max_j << " probe_q=" << probe_q
              << " omp_threads=" << omp_get_max_threads() << "\n";
    std::cout << "exact arithmetic: boost::multiprecision::cpp_int\n";
    std::cout.flush();

    std::unordered_map<std::string, cpp_int> state;
    state.emplace(std::string(), cpp_int(1));
    std::vector<cpp_int> stable(max_j + 1);
    stable[0] = 1;

    for (int t = 1; t <= max_j; ++t) {
        std::vector<std::pair<std::string, cpp_int>> items;
        items.reserve(state.size());
        for (const auto& kv : state) items.push_back(kv);

        const int threads = omp_get_max_threads();
        std::vector<std::unordered_map<std::string, cpp_int>> local(threads);
        for (auto& m : local) m.reserve((items.size() * 3) / threads + 1024);

#pragma omp parallel for schedule(dynamic, 64)
        for (size_t idx = 0; idx < items.size(); ++idx) {
            const int tid = omp_get_thread_num();
            const auto succ = horizontal_two_successors(items[idx].first);
            for (const auto& mu : succ) {
                local[tid][mu] += items[idx].second;
            }
        }

        std::unordered_map<std::string, cpp_int> next;
        size_t reserve_hint = 0;
        for (const auto& m : local) reserve_hint += m.size();
        next.reserve(reserve_hint);
        for (auto& m : local) {
            for (auto& kv : m) next[kv.first] += kv.second;
        }
        state.swap(next);

        cpp_int s = 0;
        cpp_int bw = 0;
        cpp_int bh = 0;
        for (const auto& kv : state) {
            if (!conjugate_even(kv.first)) continue;
            s += kv.second;
            if (bad_width(kv.first, probe_q)) bw += kv.second;
            if (bad_height(kv.first, probe_q)) bh += kv.second;
        }
        stable[t] = s;

        const auto now = std::chrono::steady_clock::now();
        const double elapsed =
            std::chrono::duration_cast<std::chrono::duration<double>>(now - start).count();
        std::cout << "progress j=" << t
                  << " states=" << state.size()
                  << " stable=" << decimal_prefix(s, 28)
                  << " bad_width_q=" << decimal_prefix(bw, 28)
                  << " bad_height_q=" << decimal_prefix(bh, 28)
                  << " elapsed_s=" << elapsed << "\n";
        std::cout.flush();
    }

    if (probe_q <= max_j) {
        for (int j = std::max(2 * probe_q, probe_q); j <= max_j; ++j) {
            const cpp_int rhs = (cpp_int(1) << probe_q) * binom_int(j, probe_q) * stable[j - probe_q];
            cpp_int bw = 0;
            cpp_int bh = 0;
            for (const auto& kv : state) {
                (void)kv;
            }
            std::cout << "note final state is j=" << max_j
                      << "; detailed comparison emitted only when j=max_j\n";
            break;
        }
        if (max_j >= probe_q) {
            cpp_int bw = 0;
            cpp_int bh = 0;
            for (const auto& kv : state) {
                if (!conjugate_even(kv.first)) continue;
                if (bad_width(kv.first, probe_q)) bw += kv.second;
                if (bad_height(kv.first, probe_q)) bh += kv.second;
            }
            const cpp_int rhs =
                (cpp_int(1) << probe_q) * binom_int(max_j, probe_q) * stable[max_j - probe_q];
            std::cout << "COMPARE q=" << probe_q << " j=" << max_j << "\n";
            std::cout << "rhs=" << rhs << "\n";
            std::cout << "bad_width=" << bw << "\n";
            std::cout << "bad_height=" << bh << "\n";
            std::cout << "width_le_rhs=" << (bw <= rhs ? "true" : "false") << "\n";
            std::cout << "height_le_rhs=" << (bh <= rhs ? "true" : "false") << "\n";
        }
    }
    return 0;
}
