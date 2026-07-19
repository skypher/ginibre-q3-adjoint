// LMDB-backed variant: uses disk as primary storage for the weight-multiplicity
// map.  Supports computations that exceed RAM.
//
// Storage: two LMDB environments (double-buffering for current/next state).
// Each entry: key = WeightKey bytes, value = GMP bigint serialized as raw bytes.
//
// Build:
//   g++ -O3 -march=native -fopenmp -std=c++20 -I. moments_lmdb.cpp \
//       -lgmpxx -lgmp -llmdb -o moments_lmdb

#include "lie_data.h"
#include "lmdb.h"

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

// Serialize mpz_class to raw bytes (compact: sign byte + abs value bytes).
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

template <int R>
void tensor_step(MDB_env* src, MDB_env* dst,
                 const int cartan[R][R],
                 const int adj_coord[][R],
                 const int* adj_mult,
                 int num_adj_weights) {
    // Read all entries from src (sequential scan), for each, compute
    // branching, buffer contributions, and write to dst (with merge-add).
    // Single-threaded for simplicity (LMDB writes are serialized anyway).
    MDB_txn *src_txn, *dst_txn;
    MDB_dbi src_dbi, dst_dbi;
    MDB_cursor* cur;

    E(mdb_txn_begin(src, nullptr, MDB_RDONLY, &src_txn));
    E(mdb_dbi_open(src_txn, nullptr, 0, &src_dbi));
    E(mdb_cursor_open(src_txn, src_dbi, &cur));

    E(mdb_txn_begin(dst, nullptr, 0, &dst_txn));
    E(mdb_dbi_open(dst_txn, nullptr, MDB_CREATE, &dst_dbi));
    // Clear dst
    E(mdb_drop(dst_txn, dst_dbi, 0));
    E(mdb_dbi_open(dst_txn, nullptr, MDB_CREATE, &dst_dbi));

    MDB_val k, v;
    int rc = mdb_cursor_get(cur, &k, &v, MDB_FIRST);
    size_t processed = 0;
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

            // Look up existing, add, write back
            MDB_val k2_val = {sizeof(k2), &k2};
            MDB_val existing;
            int grc = mdb_get(dst_txn, dst_dbi, &k2_val, &existing);
            if (grc == MDB_NOTFOUND) {
                auto bytes = mpz_to_bytes(contrib);
                MDB_val v2 = {bytes.size(), bytes.data()};
                E(mdb_put(dst_txn, dst_dbi, &k2_val, &v2, 0));
            } else if (grc == 0) {
                mpz_class old = bytes_to_mpz(static_cast<const uint8_t*>(existing.mv_data), existing.mv_size);
                old += contrib;
                auto bytes = mpz_to_bytes(old);
                MDB_val v2 = {bytes.size(), bytes.data()};
                E(mdb_put(dst_txn, dst_dbi, &k2_val, &v2, 0));
            } else {
                E(grc);
            }
        }
        ++processed;
        rc = mdb_cursor_get(cur, &k, &v, MDB_NEXT);
    }
    mdb_cursor_close(cur);
    mdb_txn_abort(src_txn);
    E(mdb_txn_commit(dst_txn));

    // Prune zero-valued entries from dst.
    E(mdb_txn_begin(dst, nullptr, 0, &dst_txn));
    E(mdb_dbi_open(dst_txn, nullptr, 0, &dst_dbi));
    E(mdb_cursor_open(dst_txn, dst_dbi, &cur));
    rc = mdb_cursor_get(cur, &k, &v, MDB_FIRST);
    std::vector<std::vector<uint8_t>> keys_to_delete;
    while (rc == 0) {
        mpz_class x = bytes_to_mpz(static_cast<const uint8_t*>(v.mv_data), v.mv_size);
        if (sgn(x) == 0) {
            keys_to_delete.emplace_back(static_cast<const uint8_t*>(k.mv_data),
                                        static_cast<const uint8_t*>(k.mv_data) + k.mv_size);
        }
        rc = mdb_cursor_get(cur, &k, &v, MDB_NEXT);
    }
    mdb_cursor_close(cur);
    for (auto& kbytes : keys_to_delete) {
        MDB_val dk = {kbytes.size(), kbytes.data()};
        mdb_del(dst_txn, dst_dbi, &dk, nullptr);
    }
    E(mdb_txn_commit(dst_txn));
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

size_t current_rss_kb() {
    struct rusage ru; getrusage(RUSAGE_SELF, &ru); return ru.ru_maxrss;
}

void make_dir(const std::string& p) {
    mkdir(p.c_str(), 0755);
}

MDB_env* open_env(const std::string& path, size_t map_size_tb) {
    make_dir(path);
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
             std::ostream& out) {
    std::string envA = db_prefix + "_A";
    std::string envB = db_prefix + "_B";
    MDB_env* env1 = open_env(envA, 2);
    MDB_env* env2 = open_env(envB, 2);

    // Seed env1 with V_0 = {0: 1}
    seed_with_zero<R>(env1);
    out << "m_0 = 1" << "\n"; out.flush();

    MDB_env* src = env1;
    MDB_env* dst = env2;
    for (int k = 1; k <= max_k; ++k) {
        auto t0 = std::chrono::steady_clock::now();
        tensor_step<R>(src, dst, cartan, adj_coord, adj_mult, num_adj_weights);
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
        std::fprintf(stderr, "Usage: %s GROUP MAX_K DB_PREFIX [OUT_FILE]\n", argv[0]);
        return 1;
    }
    const std::string group = argv[1];
    const int max_k = std::atoi(argv[2]);
    const std::string db_prefix = argv[3];
    std::ostream* out = &std::cout;
    std::ofstream fout;
    if (argc >= 5) { fout.open(argv[4]); out = &fout; }

#define DISPATCH(G)                                                          \
    if (group == #G)                                                         \
        return run_impl<lie_##G::rank>(lie_##G::cartan,                      \
                                       lie_##G::adj_weight_coord,            \
                                       lie_##G::adj_weight_mult,             \
                                       lie_##G::num_adj_weights, max_k,      \
                                       db_prefix, *out);
    DISPATCH(A1) DISPATCH(A2) DISPATCH(G2) DISPATCH(F4)
    DISPATCH(E6) DISPATCH(E7) DISPATCH(E8)
    std::fprintf(stderr, "Unknown group: %s\n", group.c_str());
    return 1;
}
