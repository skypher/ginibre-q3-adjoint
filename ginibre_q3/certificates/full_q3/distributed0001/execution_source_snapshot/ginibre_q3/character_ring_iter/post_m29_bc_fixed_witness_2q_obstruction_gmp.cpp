#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <omp.h>

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

static int height_of(const std::string& key) {
    return static_cast<int>(key.size());
}

static int crossed_odd_slots(int length, int q) {
    if (length <= 0) return 0;
    return std::min(q, (length + 1) / 2);
}

static std::vector<int> odd_subset(int j, int q) {
    std::vector<int> s;
    s.reserve(q);
    for (int t = 1; static_cast<int>(s.size()) < q && t <= j; t += 2) {
        s.push_back(t);
    }
    if (static_cast<int>(s.size()) != q) throw std::runtime_error("odd subset too short");
    return s;
}

static std::vector<int> subset_index(int j, const std::vector<int>& subset) {
    std::vector<int> idx(j + 1, -1);
    for (int r = 0; r < static_cast<int>(subset.size()); ++r) {
        if (subset[r] < 1 || subset[r] > j) throw std::runtime_error("subset outside [1,j]");
        idx[subset[r]] = r + 1;
    }
    return idx;
}

static bool transition_matches_height_witness(
    const std::string& old_key,
    const std::string& new_key,
    int t,
    int q,
    const std::vector<int>& idx
) {
    const int old_cross = crossed_odd_slots(height_of(old_key), q);
    const int new_cross = crossed_odd_slots(height_of(new_key), q);
    const int delta = new_cross - old_cross;
    if (delta < 0 || delta > 1) return false;

    const int want = idx[t];
    if (want < 0) return delta == 0;
    return delta == 1 && want == new_cross;
}

static mpz_class stable_count(int n) {
    std::unordered_map<std::string, mpz_class> state;
    state.emplace(std::string(), mpz_class(1));

    for (int t = 1; t <= n; ++t) {
        std::vector<std::pair<std::string, mpz_class>> items;
        items.reserve(state.size());
        for (const auto& kv : state) items.push_back(kv);

        const int threads = omp_get_max_threads();
        std::vector<std::unordered_map<std::string, mpz_class>> local(threads);
        for (auto& m : local) m.reserve((items.size() * 3) / threads + 1024);

#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < items.size(); ++i) {
            const int tid = omp_get_thread_num();
            const auto succ = horizontal_two_successors(items[i].first);
            for (const auto& mu : succ) local[tid][mu] += items[i].second;
        }

        std::unordered_map<std::string, mpz_class> next;
        size_t reserve_hint = 0;
        for (const auto& m : local) reserve_hint += m.size();
        next.reserve(reserve_hint);
        for (auto& m : local) {
            for (auto& kv : m) next[kv.first] += kv.second;
        }
        state.swap(next);
    }

    mpz_class stable = 0;
    for (const auto& kv : state) {
        if (conjugate_even(kv.first)) stable += kv.second;
    }
    return stable;
}

static mpz_class height_fixed_fiber_count(int j, int q, const std::vector<int>& subset) {
    const std::vector<int> idx = subset_index(j, subset);
    const auto start = std::chrono::steady_clock::now();

    std::unordered_map<std::string, mpz_class> state;
    state.emplace(std::string(), mpz_class(1));

    for (int t = 1; t <= j; ++t) {
        std::vector<std::pair<std::string, mpz_class>> items;
        items.reserve(state.size());
        for (const auto& kv : state) items.push_back(kv);

        const int threads = omp_get_max_threads();
        std::vector<std::unordered_map<std::string, mpz_class>> local(threads);
        for (auto& m : local) m.reserve((items.size() * 3) / threads + 1024);

#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < items.size(); ++i) {
            const int tid = omp_get_thread_num();
            const auto succ = horizontal_two_successors(items[i].first);
            for (const auto& mu : succ) {
                if (!transition_matches_height_witness(items[i].first, mu, t, q, idx)) continue;
                local[tid][mu] += items[i].second;
            }
        }

        std::unordered_map<std::string, mpz_class> next;
        size_t reserve_hint = 0;
        for (const auto& m : local) reserve_hint += m.size();
        next.reserve(reserve_hint);
        for (auto& m : local) {
            for (auto& kv : m) next[kv.first] += kv.second;
        }
        state.swap(next);

        mpz_class partial = 0;
        for (const auto& kv : state) {
            if (conjugate_even(kv.first) && height_of(kv.first) >= 2 * q) {
                partial += kv.second;
            }
        }
        const auto now = std::chrono::steady_clock::now();
        const double elapsed =
            std::chrono::duration_cast<std::chrono::duration<double>>(now - start).count();
        std::cout << "height_progress t=" << t
                  << " states=" << state.size()
                  << " partial_final_filter=" << partial.get_str()
                  << " elapsed_s=" << elapsed << "\n";
        std::cout.flush();
    }

    mpz_class total = 0;
    for (const auto& kv : state) {
        if (conjugate_even(kv.first) && height_of(kv.first) >= 2 * q) {
            total += kv.second;
        }
    }
    return total;
}

int main() {
    const int j = 30;
    const int q = 15;
    const std::vector<int> subset = odd_subset(j, q);

    std::cout << "B/C fixed-witness 2^q obstruction certificate\n";
    std::cout << "j=" << j << " q=" << q << " subset_mode=odd"
              << " omp_threads=" << omp_get_max_threads() << "\n";
    std::cout << "subset=";
    for (int t : subset) std::cout << (t == subset.front() ? "" : ",") << t;
    std::cout << "\n";
    std::cout << "exact arithmetic: GMP mpz_class\n";
    std::cout.flush();

    const mpz_class stable = stable_count(j - q);
    const mpz_class rhs = (mpz_class(1) << q) * stable;
    const mpz_class height = height_fixed_fiber_count(j, q, subset);

    const mpz_class expected_stable("147267180508");
    const mpz_class expected_rhs("4825650970886144");
    const mpz_class expected_height("6190283353629375");

    const bool stable_ok = (stable == expected_stable);
    const bool rhs_ok = (rhs == expected_rhs);
    const bool height_ok = (height == expected_height);
    const bool obstruction_ok = (height > rhs);

    std::cout << "RESULT q=" << q << " j=" << j << " subset_mode=odd\n";
    std::cout << "stable_s_j_minus_q=" << stable << "\n";
    std::cout << "rhs_2powq_s=" << rhs << "\n";
    std::cout << "height_fiber=" << height << "\n";
    std::cout << "stable_matches_expected=" << (stable_ok ? "yes" : "no") << "\n";
    std::cout << "rhs_matches_expected=" << (rhs_ok ? "yes" : "no") << "\n";
    std::cout << "height_matches_expected=" << (height_ok ? "yes" : "no") << "\n";
    std::cout << "height_exceeds_2powq_rhs=" << (obstruction_ok ? "yes" : "no") << "\n";
    std::cout << "SUMMARY failures="
              << ((stable_ok && rhs_ok && height_ok && obstruction_ok) ? 0 : 1)
              << "\n";
    std::cout << "__EXIT_STATUS="
              << ((stable_ok && rhs_ok && height_ok && obstruction_ok) ? 0 : 1)
              << "\n";

    return (stable_ok && rhs_ok && height_ok && obstruction_ok) ? 0 : 1;
}
