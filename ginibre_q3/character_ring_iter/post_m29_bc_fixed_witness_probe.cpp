#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
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

static int width_of(const std::string& key) {
    return key.empty() ? 0 : static_cast<unsigned char>(key[0]);
}

static int height_of(const std::string& key) {
    return static_cast<int>(key.size());
}

static int crossed_odd_slots(int length, int q) {
    if (length <= 0) return 0;
    return std::min(q, (length + 1) / 2);
}

static std::vector<int> make_subset(int j, int q, const std::string& mode) {
    std::vector<int> s;
    s.reserve(q);
    if (mode == "first") {
        for (int t = 1; t <= q; ++t) s.push_back(t);
    } else if (mode == "last") {
        for (int t = j - q + 1; t <= j; ++t) s.push_back(t);
    } else if (mode == "middle") {
        const int start = (j - q) / 2 + 1;
        for (int t = start; t < start + q; ++t) s.push_back(t);
    } else if (mode == "odd") {
        for (int t = 1; static_cast<int>(s.size()) < q && t <= j; t += 2) s.push_back(t);
    } else if (mode == "even") {
        for (int t = 2; static_cast<int>(s.size()) < q && t <= j; t += 2) s.push_back(t);
    } else {
        throw std::runtime_error("unknown subset mode: " + mode);
    }
    if (static_cast<int>(s.size()) != q) throw std::runtime_error("subset mode produced wrong size");
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

static bool transition_matches_witness(
    const std::string& old_key,
    const std::string& new_key,
    int t,
    int q,
    const std::vector<int>& idx,
    bool height_mode
) {
    const int old_len = height_mode ? height_of(old_key) : width_of(old_key);
    const int new_len = height_mode ? height_of(new_key) : width_of(new_key);
    const int old_cross = crossed_odd_slots(old_len, q);
    const int new_cross = crossed_odd_slots(new_len, q);
    const int delta = new_cross - old_cross;
    if (delta < 0 || delta > 1) return false;

    const int want = idx[t];
    if (want < 0) return delta == 0;
    return delta == 1 && want == new_cross;
}

static std::string decimal_prefix(const cpp_int& x, size_t digits = 48) {
    const std::string s = x.convert_to<std::string>();
    if (s.size() <= digits) return s;
    return s.substr(0, digits) + "...(" + std::to_string(s.size()) + " digits)";
}

static cpp_int stable_count(int n) {
    std::unordered_map<std::string, cpp_int> state;
    state.emplace(std::string(), cpp_int(1));
    cpp_int stable = (n == 0 ? cpp_int(1) : cpp_int(0));
    for (int t = 1; t <= n; ++t) {
        std::vector<std::pair<std::string, cpp_int>> items;
        items.reserve(state.size());
        for (const auto& kv : state) items.push_back(kv);

        const int threads = omp_get_max_threads();
        std::vector<std::unordered_map<std::string, cpp_int>> local(threads);
        for (auto& m : local) m.reserve((items.size() * 3) / threads + 1024);

#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < items.size(); ++i) {
            const int tid = omp_get_thread_num();
            const auto succ = horizontal_two_successors(items[i].first);
            for (const auto& mu : succ) local[tid][mu] += items[i].second;
        }

        std::unordered_map<std::string, cpp_int> next;
        size_t reserve_hint = 0;
        for (const auto& m : local) reserve_hint += m.size();
        next.reserve(reserve_hint);
        for (auto& m : local) {
            for (auto& kv : m) next[kv.first] += kv.second;
        }
        state.swap(next);
    }

    stable = 0;
    for (const auto& kv : state) {
        if (conjugate_even(kv.first)) stable += kv.second;
    }
    return stable;
}

static cpp_int fixed_fiber_count(
    int j,
    int q,
    const std::vector<int>& subset,
    bool height_mode
) {
    const std::vector<int> idx = subset_index(j, subset);
    const auto start = std::chrono::steady_clock::now();

    std::unordered_map<std::string, cpp_int> state;
    state.emplace(std::string(), cpp_int(1));

    for (int t = 1; t <= j; ++t) {
        std::vector<std::pair<std::string, cpp_int>> items;
        items.reserve(state.size());
        for (const auto& kv : state) items.push_back(kv);

        const int threads = omp_get_max_threads();
        std::vector<std::unordered_map<std::string, cpp_int>> local(threads);
        for (auto& m : local) m.reserve((items.size() * 3) / threads + 1024);

#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < items.size(); ++i) {
            const int tid = omp_get_thread_num();
            const auto succ = horizontal_two_successors(items[i].first);
            for (const auto& mu : succ) {
                if (!transition_matches_witness(items[i].first, mu, t, q, idx, height_mode)) continue;
                local[tid][mu] += items[i].second;
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

        cpp_int partial = 0;
        for (const auto& kv : state) {
            if (!conjugate_even(kv.first)) continue;
            if (height_mode) {
                if (height_of(kv.first) >= 2 * q) partial += kv.second;
            } else {
                if (width_of(kv.first) >= 2 * q) partial += kv.second;
            }
        }
        const auto now = std::chrono::steady_clock::now();
        const double elapsed =
            std::chrono::duration_cast<std::chrono::duration<double>>(now - start).count();
        std::cout << "progress kind=" << (height_mode ? "height" : "width")
                  << " t=" << t
                  << " states=" << state.size()
                  << " partial_final_filter=" << decimal_prefix(partial, 28)
                  << " elapsed_s=" << elapsed << "\n";
        std::cout.flush();
    }

    cpp_int total = 0;
    for (const auto& kv : state) {
        if (!conjugate_even(kv.first)) continue;
        if (height_mode) {
            if (height_of(kv.first) >= 2 * q) total += kv.second;
        } else {
            if (width_of(kv.first) >= 2 * q) total += kv.second;
        }
    }
    return total;
}

int main(int argc, char** argv) {
    int j = 30;
    int q = 15;
    std::string mode = "first";
    if (argc >= 2) j = std::atoi(argv[1]);
    if (argc >= 3) q = std::atoi(argv[2]);
    if (argc >= 4) mode = argv[3];
    if (q < 1 || j < q || j > 64) {
        std::cerr << "usage: post_m29_bc_fixed_witness_probe [j<=64] [q] [first|middle|last|odd|even]\n";
        return 2;
    }

    const std::vector<int> subset = make_subset(j, q, mode);
    std::cout << "B/C fixed-witness Pieri diagnostic\n";
    std::cout << "j=" << j << " q=" << q << " subset_mode=" << mode
              << " omp_threads=" << omp_get_max_threads() << "\n";
    std::cout << "subset=";
    for (int t : subset) std::cout << (t == subset.front() ? "" : ",") << t;
    std::cout << "\n";
    std::cout << "exact arithmetic: boost::multiprecision::cpp_int\n";
    std::cout.flush();

    const cpp_int stable = stable_count(j - q);
    const cpp_int rhs = (cpp_int(1) << q) * stable;
    std::cout << "stable_s_j_minus_q=" << stable << "\n";
    std::cout << "fixed_fiber_rhs_2powq_s=" << rhs << "\n";
    std::cout.flush();

    const cpp_int width = fixed_fiber_count(j, q, subset, false);
    const cpp_int height = fixed_fiber_count(j, q, subset, true);

    std::cout << "COMPARE_FIXED q=" << q << " j=" << j << " subset_mode=" << mode << "\n";
    std::cout << "rhs=" << rhs << "\n";
    std::cout << "width_fiber=" << width << "\n";
    std::cout << "height_fiber=" << height << "\n";
    std::cout << "width_le_rhs=" << (width <= rhs ? "true" : "false") << "\n";
    std::cout << "height_le_rhs=" << (height <= rhs ? "true" : "false") << "\n";
    std::cout << "__EXIT_STATUS=" << ((width <= rhs && height <= rhs) ? 0 : 1) << "\n";
    return (width <= rhs && height <= rhs) ? 0 : 1;
}
