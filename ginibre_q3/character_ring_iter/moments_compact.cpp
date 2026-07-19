// Memory-compact variant of moments.cpp: uses ankerl::unordered_dense
// (flat + linear probing, 2-3x less memory overhead than std::unordered_map).
// Also periodic RSS/VSZ reporting.

#include "lie_data.h"
#include "ankerl/unordered_dense.h"

#include <gmpxx.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/resource.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

template <int R>
struct WeightKey {
    std::array<int32_t, R> coord;
    bool operator==(const WeightKey& o) const { return coord == o.coord; }
};

template <int R>
struct WeightHasher {
    using is_avalanching = void;
    std::size_t operator()(const WeightKey<R>& k) const noexcept {
        std::size_t h = 1469598103934665603ULL;
        for (int i = 0; i < R; ++i) {
            uint32_t v = static_cast<uint32_t>(k.coord[i]);
            for (int b = 0; b < 4; ++b) {
                h ^= static_cast<std::size_t>((v >> (b * 8)) & 0xff);
                h *= 1099511628211ULL;
            }
        }
        return h;
    }
};

template <int R>
using Rep = ankerl::unordered_dense::map<WeightKey<R>, mpz_class, WeightHasher<R>>;

template <int R>
int dominant_reflect_inplace(int32_t* n, const int cartan[R][R]) {
    int sign = 1;
    while (true) {
        int bad = -1;
        for (int i = 0; i < R; ++i) {
            if (n[i] < 0) { bad = i; break; }
            if (n[i] == 0) { return 0; }
        }
        if (bad < 0) return sign;
        int32_t ni = n[bad];
        for (int k = 0; k < R; ++k) n[k] -= ni * cartan[bad][k];
        sign = -sign;
    }
}

size_t current_rss_kb() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss;  // on Linux, reported in KB
}

template <int R>
Rep<R> tensor_with_adj(const Rep<R>& V,
                       const int cartan[R][R],
                       const int adj_coord[][R],
                       const int* adj_mult,
                       int num_adj_weights) {
    std::vector<std::pair<WeightKey<R>, mpz_class>> items;
    items.reserve(V.size());
    for (const auto& kv : V) {
        if (sgn(kv.second) != 0) items.emplace_back(kv.first, kv.second);
    }
    const int n = static_cast<int>(items.size());

#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
#else
    const int nthreads = 1;
#endif

    std::vector<Rep<R>> locals(nthreads);
    for (auto& L : locals) L.reserve(std::max<std::size_t>(64, V.size() / nthreads));

#pragma omp parallel for schedule(dynamic, 64)
    for (int idx = 0; idx < n; ++idx) {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        Rep<R>& L = locals[tid];
        const auto& mu_key = items[idx].first;
        const mpz_class& c_mu = items[idx].second;

        for (int e = 0; e < num_adj_weights; ++e) {
            int32_t tgt[R];
            for (int i = 0; i < R; ++i)
                tgt[i] = mu_key.coord[i] + adj_coord[e][i] + 1;

            int sign = dominant_reflect_inplace<R>(tgt, cartan);
            if (sign == 0) continue;

            WeightKey<R> k2;
            for (int i = 0; i < R; ++i) k2.coord[i] = tgt[i] - 1;

            mpz_class contrib = c_mu * adj_mult[e];
            if (sign < 0) contrib = -contrib;

            auto [it, inserted] = L.try_emplace(k2, contrib);
            if (!inserted) it->second += contrib;
        }
    }

    Rep<R> W;
    std::size_t total = 0;
    for (const auto& L : locals) total += L.size();
    W.reserve(total);
    for (auto& L : locals) {
        for (auto& kv : L) {
            auto [it, inserted] = W.try_emplace(kv.first, std::move(kv.second));
            if (!inserted) it->second += kv.second;
        }
    }

    bool has_zero = false;
    for (const auto& [k, v] : W) {
        if (sgn(v) == 0) { has_zero = true; break; }
    }
    if (has_zero) {
        Rep<R> W2;
        W2.reserve(W.size());
        for (auto& [k, v] : W) {
            if (sgn(v) != 0) W2.emplace(k, std::move(v));
        }
        return W2;
    }
    return W;
}

template <int R>
int run_impl(const int cartan[R][R],
             const int adj_coord[][R],
             const int* adj_mult,
             int num_adj_weights,
             int max_k,
             std::ostream& out) {
    Rep<R> V;
    WeightKey<R> zero;
    for (int i = 0; i < R; ++i) zero.coord[i] = 0;
    V.emplace(zero, mpz_class(1));

    out << "m_0 = 1" << "\n";
    out.flush();

    for (int k = 1; k <= max_k; ++k) {
        auto t0 = std::chrono::steady_clock::now();
        V = tensor_with_adj<R>(V, cartan, adj_coord, adj_mult, num_adj_weights);
        auto t1 = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();

        mpz_class m_k = 0;
        auto it = V.find(zero);
        if (it != V.end()) m_k = it->second;

        out << "m_" << k << " = " << m_k
            << "  (|support| = " << V.size()
            << ", step_time = " << elapsed << "s"
            << ", rss_MB = " << (current_rss_kb() / 1024)
            << ")\n";
        out.flush();
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "Usage: %s GROUP MAX_K [OUT_FILE]\n", argv[0]);
        return 1;
    }
    const std::string group = argv[1];
    const int max_k = std::atoi(argv[2]);
    std::ostream* out = &std::cout;
    std::ofstream fout;
    if (argc >= 4) { fout.open(argv[3]); out = &fout; }

#define DISPATCH(G)                                                          \
    if (group == #G)                                                         \
        return run_impl<lie_##G::rank>(lie_##G::cartan,                      \
                                       lie_##G::adj_weight_coord,            \
                                       lie_##G::adj_weight_mult,             \
                                       lie_##G::num_adj_weights, max_k, *out);
    DISPATCH(A1) DISPATCH(A2) DISPATCH(G2) DISPATCH(F4)
    DISPATCH(E6) DISPATCH(E7) DISPATCH(E8)
    std::fprintf(stderr, "Unknown group: %s\n", group.c_str());
    return 1;
}
