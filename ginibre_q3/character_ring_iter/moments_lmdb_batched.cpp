// Batched LMDB variant: accumulate contributions in an in-RAM hash map
// (ankerl::unordered_dense) up to a RAM budget, then bulk-merge into LMDB.
// Dramatically reduces LMDB write pressure vs per-contribution put/get/put.
//
// Build:
//   g++ -O3 -march=native -fopenmp -std=c++20 -I. -L. moments_lmdb_batched.cpp \
//       -lgmpxx -lgmp -llmdb -o moments_lmdb_batched

#include "lie_data.h"
#include "lmdb.h"
#include "ankerl/unordered_dense.h"

#include <gmpxx.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#define E(expr) do { int _rc = (expr); if (_rc) { std::fprintf(stderr, \
    "LMDB error at %s:%d: %s (%d)\n", __FILE__, __LINE__, mdb_strerror(_rc), _rc); \
    std::abort(); } } while (0)

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
using RamMap = ankerl::unordered_dense::map<WeightKey<R>, mpz_class, WeightHasher<R>>;

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

static std::vector<uint8_t> mpz_to_bytes(const mpz_class& x) {
    int s = sgn(x);
    std::vector<uint8_t> out;
    out.push_back(s < 0 ? 1 : 0);
    if (s == 0) return out;
    size_t count;
    mpz_class absx = abs(x);
    size_t nbytes = (mpz_sizeinbase(absx.get_mpz_t(), 2) + 7) / 8;
    out.resize(1 + nbytes);
    mpz_export(&out[1], &count, 1, 1, 1, 0, absx.get_mpz_t());
    if (count < nbytes) out.resize(1 + count);
    return out;
}

static mpz_class bytes_to_mpz(const uint8_t* data, size_t n) {
    if (n == 0) return mpz_class(0);
    int sign = data[0] ? -1 : 1;
    if (n == 1) return mpz_class(0);
    mpz_class x;
    mpz_import(x.get_mpz_t(), n - 1, 1, 1, 1, 0, data + 1);
    if (sign < 0) x = -x;
    return x;
}

template <int R>
long long sum_keys(MDB_env* env) {
    MDB_txn* txn; MDB_dbi dbi; MDB_cursor* cur; MDB_val k, v;
    E(mdb_txn_begin(env, nullptr, MDB_RDONLY, &txn));
    E(mdb_dbi_open(txn, nullptr, 0, &dbi));
    E(mdb_cursor_open(txn, dbi, &cur));
    long long n = 0;
    while (mdb_cursor_get(cur, &k, &v, n == 0 ? MDB_FIRST : MDB_NEXT) == 0) ++n;
    mdb_cursor_close(cur);
    mdb_txn_abort(txn);
    return n;
}

template <int R>
void seed_with_zero(MDB_env* env) {
    MDB_txn* txn; MDB_dbi dbi;
    E(mdb_txn_begin(env, nullptr, 0, &txn));
    E(mdb_dbi_open(txn, nullptr, MDB_CREATE, &dbi));
    WeightKey<R> zero; for (int i = 0; i < R; ++i) zero.coord[i] = 0;
    mpz_class one(1);
    auto bytes = mpz_to_bytes(one);
    MDB_val k = {sizeof(zero), &zero};
    MDB_val v = {bytes.size(), bytes.data()};
    E(mdb_put(txn, dbi, &k, &v, 0));
    E(mdb_txn_commit(txn));
}

// Flush an in-RAM buffer of accumulated contributions into an LMDB env,
// merge-adding on existing keys.
template <int R>
void flush_buffer_to_env(RamMap<R>& buffer, MDB_env* env) {
    MDB_txn* txn; MDB_dbi dbi;
    E(mdb_txn_begin(env, nullptr, 0, &txn));
    E(mdb_dbi_open(txn, nullptr, MDB_CREATE, &dbi));
    for (auto& [key, value] : buffer) {
        MDB_val k = {sizeof(key), const_cast<WeightKey<R>*>(&key)};
        MDB_val existing;
        int grc = mdb_get(txn, dbi, &k, &existing);
        if (grc == MDB_NOTFOUND) {
            auto bytes = mpz_to_bytes(value);
            MDB_val v = {bytes.size(), bytes.data()};
            E(mdb_put(txn, dbi, &k, &v, 0));
        } else if (grc == 0) {
            mpz_class old = bytes_to_mpz(static_cast<const uint8_t*>(existing.mv_data), existing.mv_size);
            old += value;
            auto bytes = mpz_to_bytes(old);
            MDB_val v = {bytes.size(), bytes.data()};
            E(mdb_put(txn, dbi, &k, &v, 0));
        } else {
            E(grc);
        }
    }
    E(mdb_txn_commit(txn));
    buffer.clear();
}

template <int R>
void tensor_step(MDB_env* src, MDB_env* dst,
                 const int cartan[R][R],
                 const int adj_coord[][R],
                 const int* adj_mult,
                 int num_adj_weights,
                 size_t buffer_size_limit) {
    MDB_txn* src_txn; MDB_dbi src_dbi; MDB_cursor* cur;
    E(mdb_txn_begin(src, nullptr, MDB_RDONLY, &src_txn));
    E(mdb_dbi_open(src_txn, nullptr, 0, &src_dbi));
    E(mdb_cursor_open(src_txn, src_dbi, &cur));

    // Clear dst
    {
        MDB_txn* dt; MDB_dbi dd;
        E(mdb_txn_begin(dst, nullptr, 0, &dt));
        E(mdb_dbi_open(dt, nullptr, MDB_CREATE, &dd));
        E(mdb_drop(dt, dd, 0));
        E(mdb_txn_commit(dt));
    }

    RamMap<R> buffer;
    buffer.reserve(buffer_size_limit * 2);  // Rough preallocation
    int n_flushes = 0;

    MDB_val k, v;
    int rc = mdb_cursor_get(cur, &k, &v, MDB_FIRST);
    while (rc == 0) {
        if (k.mv_size != sizeof(WeightKey<R>)) std::abort();
        WeightKey<R> mu_key;
        std::memcpy(&mu_key, k.mv_data, sizeof(mu_key));
        mpz_class c_mu = bytes_to_mpz(static_cast<const uint8_t*>(v.mv_data), v.mv_size);

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

            auto [it, inserted] = buffer.try_emplace(k2, contrib);
            if (!inserted) it->second += contrib;
        }

        if (buffer.size() >= buffer_size_limit) {
            flush_buffer_to_env(buffer, dst);
            ++n_flushes;
        }
        rc = mdb_cursor_get(cur, &k, &v, MDB_NEXT);
    }
    mdb_cursor_close(cur);
    mdb_txn_abort(src_txn);

    if (!buffer.empty()) { flush_buffer_to_env(buffer, dst); ++n_flushes; }

    // Prune zeros
    {
        MDB_txn* dt; MDB_dbi dd; MDB_cursor* dc;
        E(mdb_txn_begin(dst, nullptr, 0, &dt));
        E(mdb_dbi_open(dt, nullptr, 0, &dd));
        E(mdb_cursor_open(dt, dd, &dc));
        std::vector<std::vector<uint8_t>> to_del;
        MDB_val dk, dv;
        int drc = mdb_cursor_get(dc, &dk, &dv, MDB_FIRST);
        while (drc == 0) {
            mpz_class x = bytes_to_mpz(static_cast<const uint8_t*>(dv.mv_data), dv.mv_size);
            if (sgn(x) == 0) {
                to_del.emplace_back(static_cast<const uint8_t*>(dk.mv_data),
                                    static_cast<const uint8_t*>(dk.mv_data) + dk.mv_size);
            }
            drc = mdb_cursor_get(dc, &dk, &dv, MDB_NEXT);
        }
        mdb_cursor_close(dc);
        for (auto& kb : to_del) {
            MDB_val ddk = {kb.size(), kb.data()};
            mdb_del(dt, dd, &ddk, nullptr);
        }
        E(mdb_txn_commit(dt));
    }
}

template <int R>
mpz_class read_m_0(MDB_env* env) {
    MDB_txn* txn; MDB_dbi dbi; MDB_val k, v;
    E(mdb_txn_begin(env, nullptr, MDB_RDONLY, &txn));
    E(mdb_dbi_open(txn, nullptr, 0, &dbi));
    WeightKey<R> zero; for (int i = 0; i < R; ++i) zero.coord[i] = 0;
    MDB_val kk = {sizeof(zero), &zero};
    int rc = mdb_get(txn, dbi, &kk, &v);
    mpz_class result(0);
    if (rc == 0) {
        result = bytes_to_mpz(static_cast<const uint8_t*>(v.mv_data), v.mv_size);
    }
    mdb_txn_abort(txn);
    return result;
}

size_t current_rss_kb() { struct rusage ru; getrusage(RUSAGE_SELF, &ru); return ru.ru_maxrss; }

MDB_env* open_env(const std::string& path, size_t map_size_tb) {
    mkdir(path.c_str(), 0755);
    MDB_env* env;
    E(mdb_env_create(&env));
    E(mdb_env_set_mapsize(env, map_size_tb * (1ULL << 40)));
    E(mdb_env_set_maxdbs(env, 4));
    E(mdb_env_open(env, path.c_str(), MDB_NOSYNC | MDB_WRITEMAP, 0644));
    return env;
}

template <int R>
int run_impl(const int cartan[R][R],
             const int adj_coord[][R],
             const int* adj_mult,
             int num_adj_weights,
             int max_k,
             const std::string& db_prefix,
             size_t buffer_limit,
             std::ostream& out) {
    std::string envA = db_prefix + "_A";
    std::string envB = db_prefix + "_B";
    MDB_env* env1 = open_env(envA, 4);
    MDB_env* env2 = open_env(envB, 4);

    seed_with_zero<R>(env1);
    out << "m_0 = 1" << "\n"; out.flush();

    MDB_env* src = env1;
    MDB_env* dst = env2;
    for (int k = 1; k <= max_k; ++k) {
        auto t0 = std::chrono::steady_clock::now();
        tensor_step<R>(src, dst, cartan, adj_coord, adj_mult, num_adj_weights, buffer_limit);
        auto t1 = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();
        mpz_class m_k = read_m_0<R>(dst);
        long long supp = sum_keys<R>(dst);
        out << "m_" << k << " = " << m_k
            << "  (|support| = " << supp
            << ", step_time = " << elapsed << "s"
            << ", rss_MB = " << (current_rss_kb() / 1024) << ")\n";
        out.flush();
        std::swap(src, dst);
    }
    mdb_env_close(env1);
    mdb_env_close(env2);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "Usage: %s GROUP MAX_K DB_PREFIX [OUT_FILE] [BUFFER_SIZE]\n", argv[0]);
        std::fprintf(stderr, "BUFFER_SIZE: max entries in RAM before flush (default 1000000)\n");
        return 1;
    }
    const std::string group = argv[1];
    const int max_k = std::atoi(argv[2]);
    const std::string db_prefix = argv[3];
    std::ostream* out = &std::cout;
    std::ofstream fout;
    size_t buffer_limit = 1000000;
    if (argc >= 5) { fout.open(argv[4]); out = &fout; }
    if (argc >= 6) { buffer_limit = std::atoll(argv[5]); }

#define DISPATCH(G)                                                          \
    if (group == #G)                                                         \
        return run_impl<lie_##G::rank>(lie_##G::cartan,                      \
                                       lie_##G::adj_weight_coord,            \
                                       lie_##G::adj_weight_mult,             \
                                       lie_##G::num_adj_weights, max_k,      \
                                       db_prefix, buffer_limit, *out);
    DISPATCH(A1) DISPATCH(A2) DISPATCH(G2) DISPATCH(F4)
    DISPATCH(E6) DISPATCH(E7) DISPATCH(E8)
    std::fprintf(stderr, "Unknown group: %s\n", group.c_str());
    return 1;
}
