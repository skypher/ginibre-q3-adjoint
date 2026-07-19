// Exact modular rank replay for the F4 minimal degree-five invariant tensors.
//
// The computation uses the traceless Albert-algebra model of the 26-dimensional
// F4 module.  The basis is orthogonal with norms 2, 6, 2, ..., 2.  All tensor
// contractions are performed in F_p, p = 1000000007.

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using std::array;
using std::int64_t;
using std::pair;
using std::tuple;
using std::vector;

#ifndef F4_MOD
#define F4_MOD 1000000007
#endif
static constexpr int MOD = F4_MOD;
static constexpr int DIM = 26;
static constexpr int MAT_DIM = DIM * DIM;

struct GraphSpec;

static int add_mod(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

static int sub_mod(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}

static int mul_mod(int64_t a, int64_t b) {
    return int((a * b) % MOD);
}

static int pow_mod(int a, int e) {
    int64_t r = 1, x = a;
    while (e) {
        if (e & 1) r = (r * x) % MOD;
        x = (x * x) % MOD;
        e >>= 1;
    }
    return int(r);
}

static int inv_mod(int a) {
    return pow_mod(a, MOD - 2);
}

using Oct = array<int, 8>;
using Mat = array<array<Oct, 3>, 3>;

static int mult_sign[8][8];
static int mult_idx[8][8];

static Oct oct_zero() {
    Oct z{};
    z.fill(0);
    return z;
}

static void init_octonions() {
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            mult_sign[i][j] = 0;
            mult_idx[i][j] = 0;
        }
    }
    for (int i = 0; i < 8; ++i) {
        mult_sign[0][i] = mult_sign[i][0] = 1;
        mult_idx[0][i] = mult_idx[i][0] = i;
    }
    for (int i = 1; i < 8; ++i) {
        mult_sign[i][i] = -1;
        mult_idx[i][i] = 0;
    }
    vector<array<int, 3>> triples = {
        array<int, 3>{1, 2, 3}, array<int, 3>{1, 4, 5},
        array<int, 3>{1, 7, 6}, array<int, 3>{2, 4, 6},
        array<int, 3>{2, 5, 7}, array<int, 3>{3, 4, 7},
        array<int, 3>{3, 6, 5},
    };
    for (auto t : triples) {
        int a = t[0], b = t[1], c = t[2];
        vector<array<int, 3>> cyc = {
            array<int, 3>{a, b, c}, array<int, 3>{b, c, a},
            array<int, 3>{c, a, b},
        };
        for (auto u : cyc) {
            int x = u[0], y = u[1], z = u[2];
            mult_sign[x][y] = 1;
            mult_idx[x][y] = z;
            mult_sign[y][x] = -1;
            mult_idx[y][x] = z;
        }
    }
}

static Oct oct_mul(const Oct& x, const Oct& y) {
    Oct out = oct_zero();
    for (int i = 0; i < 8; ++i) if (x[i]) {
        for (int j = 0; j < 8; ++j) if (y[j]) {
            int s = mult_sign[i][j];
            int k = mult_idx[i][j];
            int term = mul_mod(x[i], y[j]);
            if (s == 1) out[k] = add_mod(out[k], term);
            else out[k] = sub_mod(out[k], term);
        }
    }
    return out;
}

static Oct oct_conj(Oct x) {
    for (int i = 1; i < 8; ++i) {
        if (x[i]) x[i] = MOD - x[i];
    }
    return x;
}

static Mat mat_zero() {
    Mat A;
    for (auto& row : A) {
        for (auto& x : row) x = oct_zero();
    }
    return A;
}

static Mat mat_mul(const Mat& A, const Mat& B) {
    Mat C = mat_zero();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Oct s = oct_zero();
            for (int k = 0; k < 3; ++k) {
                Oct p = oct_mul(A[i][k], B[k][j]);
                for (int m = 0; m < 8; ++m) s[m] = add_mod(s[m], p[m]);
            }
            C[i][j] = s;
        }
    }
    return C;
}

static Mat jordan(const Mat& A, const Mat& B) {
    static const int inv2 = inv_mod(2);
    Mat AB = mat_mul(A, B), BA = mat_mul(B, A), C = mat_zero();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int m = 0; m < 8; ++m) {
                C[i][j][m] = mul_mod(add_mod(AB[i][j][m], BA[i][j][m]), inv2);
            }
        }
    }
    return C;
}

static int trace(const Mat& A) {
    return add_mod(add_mod(A[0][0][0], A[1][1][0]), A[2][2][0]);
}

static Mat diag_mat(int a, int b, int c) {
    Mat A = mat_zero();
    A[0][0][0] = (a % MOD + MOD) % MOD;
    A[1][1][0] = (b % MOD + MOD) % MOD;
    A[2][2][0] = (c % MOD + MOD) % MOD;
    return A;
}

static Mat off_mat(int p, int q, int m) {
    Mat A = mat_zero();
    A[p][q][m] = 1;
    Oct e = oct_zero();
    e[m] = 1;
    A[q][p] = oct_conj(e);
    return A;
}

struct Entry {
    vector<int> key;
    int val;
};

struct Factor {
    vector<int> vars;
    vector<Entry> entries;
};

static vector<array<int, 3>> cubic_entries;
static vector<int> cubic_values;
static vector<int> cubic_values_unit;
static vector<int> norms;
static vector<int> qnorms;

static int rank_mod(vector<vector<int>> A);
static vector<vector<int>> inverse_mod_matrix(vector<vector<int>> A);
static int m123_unit_entry(const GraphSpec& a, int ia, const GraphSpec& b, int ib,
                           const GraphSpec& c, int ic, const GraphSpec& dspec, int id);

static Factor cubic_factor(const vector<int>& vars) {
    vector<Entry> entries;
    entries.reserve(cubic_entries.size());
    for (size_t t = 0; t < cubic_entries.size(); ++t) {
        entries.push_back({{cubic_entries[t][0], cubic_entries[t][1], cubic_entries[t][2]},
                           cubic_values[t]});
    }
    return {vars, std::move(entries)};
}

static Factor cubic_factor_unit(const vector<int>& vars) {
    vector<Entry> entries;
    entries.reserve(cubic_entries.size());
    for (size_t t = 0; t < cubic_entries.size(); ++t) {
        entries.push_back({{cubic_entries[t][0], cubic_entries[t][1], cubic_entries[t][2]},
                           cubic_values_unit[t]});
    }
    return {vars, std::move(entries)};
}

static Factor delta_factor(int a, int b) {
    vector<Entry> entries;
    for (int i = 0; i < DIM; ++i) entries.push_back({{i, i}, 1});
    return {{a, b}, entries};
}

static Factor g_factor(int a, int b) {
    vector<Entry> entries;
    for (int i = 0; i < DIM; ++i) entries.push_back({{i, i}, norms[i]});
    return {{a, b}, entries};
}

static Factor q_factor(int a) {
    vector<Entry> entries;
    for (int i = 0; i < DIM; ++i) entries.push_back({{i}, qnorms[i]});
    return {{a}, entries};
}

static Factor q2_factor(int a, int b) {
    vector<Entry> entries;
    for (int i = 0; i < DIM; ++i) entries.push_back({{i, i}, qnorms[i]});
    return {{a, b}, entries};
}

static Factor delta_factor_unit(int a, int b) {
    return delta_factor(a, b);
}

static int eliminate_value(vector<Factor> factors) {
    int scalar = 1;
    while (true) {
        vector<Factor> nf;
        for (auto& f : factors) {
            if (f.vars.empty()) {
                int s = 0;
                for (auto& e : f.entries) s = add_mod(s, e.val);
                scalar = mul_mod(scalar, s);
            } else {
                nf.push_back(std::move(f));
            }
        }
        factors.swap(nf);
        if (factors.empty()) return scalar;

        std::set<int> vars_set;
        for (auto& f : factors) for (int v : f.vars) vars_set.insert(v);
        bool have = false;
        tuple<int64_t, int64_t, int, int> best;
        int best_v = -1;
        for (int v : vars_set) {
            vector<int> inc;
            std::set<int> uni;
            int64_t est = 0;
            for (int i = 0; i < (int)factors.size(); ++i) {
                auto it = std::find(factors[i].vars.begin(), factors[i].vars.end(), v);
                if (it != factors[i].vars.end()) {
                    inc.push_back(i);
                    for (int x : factors[i].vars) if (x != v) uni.insert(x);
                }
            }
            for (int aval = 0; aval < DIM; ++aval) {
                int64_t prod = 1;
                for (int idx : inc) {
                    const auto& f = factors[idx];
                    int pos = int(std::find(f.vars.begin(), f.vars.end(), v) - f.vars.begin());
                    int cnt = 0;
                    for (auto& e : f.entries) if (e.key[pos] == aval) ++cnt;
                    prod *= cnt;
                    if (prod > (int64_t)4e18 / 64) break;
                }
                est += prod;
            }
            tuple<int64_t, int64_t, int, int> cost(est, (int64_t)uni.size(), -int(inc.size()), v);
            if (!have || cost < best) {
                have = true;
                best = cost;
                best_v = v;
            }
        }

        int v = best_v;
        vector<Factor> inc_f, rest;
        for (auto& f : factors) {
            if (std::find(f.vars.begin(), f.vars.end(), v) != f.vars.end()) inc_f.push_back(std::move(f));
            else rest.push_back(std::move(f));
        }
        vector<int> uni;
        for (auto& f : inc_f) {
            for (int x : f.vars) if (x != v && std::find(uni.begin(), uni.end(), x) == uni.end()) {
                uni.push_back(x);
            }
        }

        struct GroupEntry {
            vector<int> other_vals;
            int val;
            vector<int> other_vars;
        };
        vector<array<vector<GroupEntry>, DIM>> groups(inc_f.size());
        for (int t = 0; t < (int)inc_f.size(); ++t) {
            auto& f = inc_f[t];
            int p = int(std::find(f.vars.begin(), f.vars.end(), v) - f.vars.begin());
            vector<int> other_vars;
            vector<int> other_pos;
            for (int i = 0; i < (int)f.vars.size(); ++i) if (i != p) {
                other_vars.push_back(f.vars[i]);
                other_pos.push_back(i);
            }
            for (auto& e : f.entries) {
                vector<int> vals;
                for (int pos : other_pos) vals.push_back(e.key[pos]);
                groups[t][e.key[p]].push_back({vals, e.val, other_vars});
            }
        }

        std::map<vector<int>, int> out;
        vector<int> assign_val(uni.size(), -1);
        auto set_assign = [&](int var, int val, vector<pair<int, int>>& changed) -> bool {
            int pos = int(std::find(uni.begin(), uni.end(), var) - uni.begin());
            if (assign_val[pos] != -1) return assign_val[pos] == val;
            assign_val[pos] = val;
            changed.push_back({pos, val});
            return true;
        };
        std::function<void(int, int)> rec = [&](int t, int coef) {
            if (t == (int)inc_f.size()) {
                vector<int> key = assign_val;
                int old = out[key];
                out[key] = add_mod(old, coef);
                return;
            }
            int aval_dummy = 0;
            (void)aval_dummy;
        };
        for (int aval = 0; aval < DIM; ++aval) {
            bool empty = false;
            for (auto& g : groups) if (g[aval].empty()) empty = true;
            if (empty) continue;
            std::function<void(int, int)> rec2 = [&](int t, int coef) {
                if (t == (int)groups.size()) {
                    vector<int> key = assign_val;
                    int old = out[key];
                    out[key] = add_mod(old, coef);
                    return;
                }
                for (const auto& ge : groups[t][aval]) {
                    vector<pair<int, int>> changed;
                    bool ok = true;
                    for (int i = 0; i < (int)ge.other_vars.size(); ++i) {
                        if (!set_assign(ge.other_vars[i], ge.other_vals[i], changed)) {
                            ok = false;
                            break;
                        }
                    }
                    if (ok) rec2(t + 1, mul_mod(coef, ge.val));
                    for (auto it = changed.rbegin(); it != changed.rend(); ++it) {
                        assign_val[it->first] = -1;
                    }
                }
            };
            rec2(0, 1);
        }
        vector<Entry> entries;
        entries.reserve(out.size());
        for (auto& kv : out) if (kv.second) entries.push_back({kv.first, kv.second});
        rest.push_back({uni, std::move(entries)});
        factors.swap(rest);
    }
}

static Factor join_factors(const Factor& a, const Factor& b) {
    vector<int> vars = a.vars;
    for (int v : b.vars) if (std::find(vars.begin(), vars.end(), v) == vars.end()) vars.push_back(v);
    std::sort(vars.begin(), vars.end());
    vector<int> apos, bpos, ashared, bshared;
    for (int v : a.vars) apos.push_back(int(std::find(vars.begin(), vars.end(), v) - vars.begin()));
    for (int v : b.vars) bpos.push_back(int(std::find(vars.begin(), vars.end(), v) - vars.begin()));
    for (int i = 0; i < (int)a.vars.size(); ++i) {
        auto it = std::find(b.vars.begin(), b.vars.end(), a.vars[i]);
        if (it != b.vars.end()) {
            ashared.push_back(i);
            bshared.push_back(int(it - b.vars.begin()));
        }
    }
    std::map<vector<int>, vector<const Entry*>> b_by_shared;
    for (const auto& e : b.entries) {
        vector<int> key;
        key.reserve(bshared.size());
        for (int p : bshared) key.push_back(e.key[p]);
        b_by_shared[key].push_back(&e);
    }
    std::map<vector<int>, int> out;
    for (const auto& ea : a.entries) {
        vector<int> skey;
        skey.reserve(ashared.size());
        for (int p : ashared) skey.push_back(ea.key[p]);
        auto it = b_by_shared.find(skey);
        if (it == b_by_shared.end()) continue;
        for (const Entry* eb : it->second) {
            vector<int> key(vars.size(), -1);
            for (int i = 0; i < (int)a.vars.size(); ++i) key[apos[i]] = ea.key[i];
            bool ok = true;
            for (int i = 0; i < (int)b.vars.size(); ++i) {
                int p = bpos[i];
                if (key[p] != -1 && key[p] != eb->key[i]) {
                    ok = false;
                    break;
                }
                key[p] = eb->key[i];
            }
            if (!ok) continue;
            int val = mul_mod(ea.val, eb->val);
            out[key] = add_mod(out[key], val);
        }
    }
    vector<Entry> entries;
    entries.reserve(out.size());
    for (auto& kv : out) if (kv.second) entries.push_back({kv.first, kv.second});
    return {vars, std::move(entries)};
}

static Factor eliminate_to_factor(vector<Factor> factors, const vector<int>& keep_vars) {
    std::set<int> keep(keep_vars.begin(), keep_vars.end());
    int scalar = 1;
    while (true) {
        vector<Factor> nf;
        for (auto& f : factors) {
            if (f.vars.empty()) {
                int s = 0;
                for (auto& e : f.entries) s = add_mod(s, e.val);
                scalar = mul_mod(scalar, s);
            } else {
                nf.push_back(std::move(f));
            }
        }
        factors.swap(nf);

        std::set<int> elim_set;
        for (auto& f : factors) {
            for (int v : f.vars) if (!keep.count(v)) elim_set.insert(v);
        }
        if (elim_set.empty()) break;

        bool have = false;
        tuple<int64_t, int64_t, int, int> best;
        int best_v = -1;
        for (int v : elim_set) {
            vector<int> inc;
            std::set<int> uni;
            int64_t est = 0;
            for (int i = 0; i < (int)factors.size(); ++i) {
                auto it = std::find(factors[i].vars.begin(), factors[i].vars.end(), v);
                if (it != factors[i].vars.end()) {
                    inc.push_back(i);
                    for (int x : factors[i].vars) if (x != v) uni.insert(x);
                }
            }
            for (int aval = 0; aval < DIM; ++aval) {
                int64_t prod = 1;
                for (int idx : inc) {
                    const auto& f = factors[idx];
                    int pos = int(std::find(f.vars.begin(), f.vars.end(), v) - f.vars.begin());
                    int cnt = 0;
                    for (auto& e : f.entries) if (e.key[pos] == aval) ++cnt;
                    prod *= cnt;
                    if (prod > (int64_t)4e18 / 64) break;
                }
                est += prod;
            }
            tuple<int64_t, int64_t, int, int> cost(est, (int64_t)uni.size(), -int(inc.size()), v);
            if (!have || cost < best) {
                have = true;
                best = cost;
                best_v = v;
            }
        }

        int v = best_v;
        vector<Factor> inc_f, rest;
        for (auto& f : factors) {
            if (std::find(f.vars.begin(), f.vars.end(), v) != f.vars.end()) inc_f.push_back(std::move(f));
            else rest.push_back(std::move(f));
        }
        vector<int> uni;
        for (auto& f : inc_f) {
            for (int x : f.vars) if (x != v && std::find(uni.begin(), uni.end(), x) == uni.end()) {
                uni.push_back(x);
            }
        }

        struct GroupEntry {
            vector<int> other_vals;
            int val;
            vector<int> other_vars;
        };
        vector<array<vector<GroupEntry>, DIM>> groups(inc_f.size());
        for (int t = 0; t < (int)inc_f.size(); ++t) {
            auto& f = inc_f[t];
            int p = int(std::find(f.vars.begin(), f.vars.end(), v) - f.vars.begin());
            vector<int> other_vars;
            vector<int> other_pos;
            for (int i = 0; i < (int)f.vars.size(); ++i) if (i != p) {
                other_vars.push_back(f.vars[i]);
                other_pos.push_back(i);
            }
            for (auto& e : f.entries) {
                vector<int> vals;
                for (int pos : other_pos) vals.push_back(e.key[pos]);
                groups[t][e.key[p]].push_back({vals, e.val, other_vars});
            }
        }

        std::map<vector<int>, int> out;
        vector<int> assign_val(uni.size(), -1);
        auto set_assign = [&](int var, int val, vector<pair<int, int>>& changed) -> bool {
            int pos = int(std::find(uni.begin(), uni.end(), var) - uni.begin());
            if (assign_val[pos] != -1) return assign_val[pos] == val;
            assign_val[pos] = val;
            changed.push_back({pos, val});
            return true;
        };
        for (int aval = 0; aval < DIM; ++aval) {
            bool empty = false;
            for (auto& g : groups) if (g[aval].empty()) empty = true;
            if (empty) continue;
            std::function<void(int, int)> rec = [&](int t, int coef) {
                if (t == (int)groups.size()) {
                    vector<int> key = assign_val;
                    out[key] = add_mod(out[key], coef);
                    return;
                }
                for (const auto& ge : groups[t][aval]) {
                    vector<pair<int, int>> changed;
                    bool ok = true;
                    for (int i = 0; i < (int)ge.other_vars.size(); ++i) {
                        if (!set_assign(ge.other_vars[i], ge.other_vals[i], changed)) {
                            ok = false;
                            break;
                        }
                    }
                    if (ok) rec(t + 1, mul_mod(coef, ge.val));
                    for (auto it = changed.rbegin(); it != changed.rend(); ++it) {
                        assign_val[it->first] = -1;
                    }
                }
            };
            rec(0, 1);
        }
        vector<Entry> entries;
        entries.reserve(out.size());
        for (auto& kv : out) if (kv.second) entries.push_back({kv.first, kv.second});
        rest.push_back({uni, std::move(entries)});
        factors.swap(rest);
    }

    if (factors.empty()) {
        return {{}, {{{}, scalar}}};
    }
    Factor out = std::move(factors[0]);
    for (int i = 1; i < (int)factors.size(); ++i) out = join_factors(out, factors[i]);
    if (scalar != 1) {
        for (auto& e : out.entries) e.val = mul_mod(e.val, scalar);
    }
    if (out.vars != keep_vars) {
        vector<int> pos;
        pos.reserve(keep_vars.size());
        for (int v : keep_vars) pos.push_back(int(std::find(out.vars.begin(), out.vars.end(), v) - out.vars.begin()));
        for (auto& e : out.entries) {
            vector<int> key;
            key.reserve(keep_vars.size());
            for (int p : pos) key.push_back(e.key[p]);
            e.key = std::move(key);
        }
        out.vars = keep_vars;
    }
    return out;
}

static vector<vector<pair<int, int>>> all_pairings_rec(vector<int> items) {
    if (items.empty()) return {{}};
    int first = items[0];
    vector<vector<pair<int, int>>> out;
    for (int i = 1; i < (int)items.size(); ++i) {
        vector<int> rest;
        for (int j = 1; j < (int)items.size(); ++j) if (j != i) rest.push_back(items[j]);
        for (auto p : all_pairings_rec(rest)) {
            p.insert(p.begin(), {first, items[i]});
            out.push_back(p);
        }
    }
    return out;
}

struct TensorSpec {
    int kind; // 0 delta*d, 1 chain
    pair<int, int> pair0;
    array<int, 3> tri;
    int single;
    pair<int, int> pair1;
    pair<int, int> pair2;
};

struct GraphSpec {
    int k = 0;
    vector<pair<int, int>> metrics;
    vector<array<int, 3>> cubics;
    int internal_edges = 0;
};

static vector<Factor> tensor_factors(const TensorSpec& s, int prefix) {
    auto e = [&](int r) { return prefix * 100 + r; };
    vector<Factor> f;
    if (s.kind == 0) {
        f.push_back(g_factor(e(s.pair0.first), e(s.pair0.second)));
        f.push_back(cubic_factor({e(s.tri[0]), e(s.tri[1]), e(s.tri[2])}));
    } else {
        int x = prefix * 100 + 20;
        int y = prefix * 100 + 21;
        f.push_back(cubic_factor({e(s.pair1.first), e(s.pair1.second), x}));
        f.push_back(q_factor(x));
        f.push_back(cubic_factor({e(s.single), x, y}));
        f.push_back(q_factor(y));
        f.push_back(cubic_factor({e(s.pair2.first), e(s.pair2.second), y}));
    }
    return f;
}

static vector<Factor> graph_factors(const GraphSpec& s, int prefix) {
    auto e = [&](int r) { return prefix * 100 + r; };
    auto v = [&](int leg) {
        if (leg >= 0) return e(leg);
        return prefix * 100 + 50 + (-1 - leg);
    };
    vector<Factor> f;
    for (auto pr : s.metrics) f.push_back(g_factor(e(pr.first), e(pr.second)));
    for (auto cu : s.cubics) f.push_back(cubic_factor({v(cu[0]), v(cu[1]), v(cu[2])}));
    for (int q = 0; q < s.internal_edges; ++q) f.push_back(q_factor(prefix * 100 + 50 + q));
    return f;
}

static vector<Factor> graph_unit_factors(const GraphSpec& s, const array<int, 6>& ext, int prefix) {
    auto e = [&](int r) { return ext[r]; };
    int base = prefix * 100 + 50;
    auto v = [&](int leg) {
        if (leg >= 0) return e(leg);
        return base + (-1 - leg);
    };
    vector<Factor> f;
    for (auto pr : s.metrics) f.push_back(delta_factor_unit(e(pr.first), e(pr.second)));
    for (auto cu : s.cubics) f.push_back(cubic_factor_unit({v(cu[0]), v(cu[1]), v(cu[2])}));
    return f;
}

static vector<Factor> tensor_unit_factors(const TensorSpec& s, const array<int, 5>& ext, int prefix) {
    auto e = [&](int r) { return ext[r]; };
    vector<Factor> f;
    int base = prefix * 1000 + 100;
    if (s.kind == 0) {
        f.push_back(delta_factor(e(s.pair0.first), e(s.pair0.second)));
        f.push_back(cubic_factor_unit({e(s.tri[0]), e(s.tri[1]), e(s.tri[2])}));
    } else {
        int x = base + 20;
        int y = base + 21;
        f.push_back(cubic_factor_unit({e(s.pair1.first), e(s.pair1.second), x}));
        f.push_back(cubic_factor_unit({e(s.single), x, y}));
        f.push_back(cubic_factor_unit({e(s.pair2.first), e(s.pair2.second), y}));
    }
    return f;
}

static int gram_entry(const TensorSpec& a, int ia, const TensorSpec& b, int ib) {
    vector<Factor> f = tensor_factors(a, ia + 1);
    vector<Factor> g = tensor_factors(b, ib + 51);
    f.insert(f.end(), g.begin(), g.end());
    for (int r = 0; r < 5; ++r) f.push_back(q2_factor((ia + 1) * 100 + r, (ib + 51) * 100 + r));
    return eliminate_value(std::move(f));
}

static int graph_gram_entry(const GraphSpec& a, int ia, const GraphSpec& b, int ib) {
    vector<Factor> f = graph_factors(a, ia + 1);
    vector<Factor> g = graph_factors(b, ib + 251);
    f.insert(f.end(), g.begin(), g.end());
    for (int r = 0; r < a.k; ++r) f.push_back(q2_factor((ia + 1) * 100 + r, (ib + 251) * 100 + r));
    return eliminate_value(std::move(f));
}

static Factor graph_unit_tensor_factor(const GraphSpec& s, int prefix) {
    array<int, 6> ext = {0, 1, 2, 3, 4, 5};
    return eliminate_to_factor(graph_unit_factors(s, ext, prefix), {0, 1, 2, 3, 4, 5});
}

static int graph_unit_gram_entry(const GraphSpec& a, int ia, const GraphSpec& b, int ib) {
    array<int, 6> ext = {0, 1, 2, 3, 4, 5};
    vector<Factor> f = graph_unit_factors(a, ext, ia + 1);
    vector<Factor> g = graph_unit_factors(b, ext, ib + 251);
    f.insert(f.end(), g.begin(), g.end());
    return eliminate_value(std::move(f));
}

static vector<int> complement_items(const vector<int>& all, const vector<int>& used) {
    vector<int> out;
    for (int x : all) if (std::find(used.begin(), used.end(), x) == used.end()) out.push_back(x);
    return out;
}

static void add_l6_metric_pairings(vector<GraphSpec>& specs) {
    vector<int> items = {0, 1, 2, 3, 4, 5};
    for (auto ps : all_pairings_rec(items)) {
        GraphSpec s;
        s.k = 6;
        s.metrics = ps;
        specs.push_back(std::move(s));
    }
}

static void add_l6_two_cubics(vector<GraphSpec>& specs) {
    vector<int> items = {0, 1, 2, 3, 4, 5};
    for (int a = 1; a < 6; ++a) {
        for (int b = a + 1; b < 6; ++b) {
            vector<int> tri0 = {0, a, b};
            vector<int> tri1 = complement_items(items, tri0);
            GraphSpec s;
            s.k = 6;
            s.cubics.push_back({tri0[0], tri0[1], tri0[2]});
            s.cubics.push_back({tri1[0], tri1[1], tri1[2]});
            specs.push_back(std::move(s));
        }
    }
}

static void add_l6_metric_connected_cubics(vector<GraphSpec>& specs) {
    vector<int> items = {0, 1, 2, 3, 4, 5};
    for (int a = 0; a < 6; ++a) {
        for (int b = a + 1; b < 6; ++b) {
            vector<int> used = {a, b};
            vector<int> rem = complement_items(items, used);
            for (auto ps : all_pairings_rec(rem)) {
                GraphSpec s;
                s.k = 6;
                s.metrics.push_back({a, b});
                s.cubics.push_back({ps[0].first, ps[0].second, -1});
                s.cubics.push_back({ps[1].first, ps[1].second, -1});
                s.internal_edges = 1;
                specs.push_back(std::move(s));
            }
        }
    }
}

static void add_l6_cubic_paths(vector<GraphSpec>& specs) {
    vector<int> items = {0, 1, 2, 3, 4, 5};
    for (int s1 = 0; s1 < 6; ++s1) {
        for (int s2 = 0; s2 < 6; ++s2) if (s2 != s1) {
            vector<int> used = {s1, s2};
            vector<int> rem = complement_items(items, used);
            for (auto ps : all_pairings_rec(rem)) {
                GraphSpec s;
                s.k = 6;
                s.cubics.push_back({ps[0].first, ps[0].second, -1});
                s.cubics.push_back({s1, -1, -2});
                s.cubics.push_back({s2, -2, -3});
                s.cubics.push_back({ps[1].first, ps[1].second, -3});
                s.internal_edges = 3;
                specs.push_back(std::move(s));
            }
        }
    }
}

static void add_l6_cubic_stars(vector<GraphSpec>& specs) {
    vector<int> items = {0, 1, 2, 3, 4, 5};
    for (auto ps : all_pairings_rec(items)) {
        GraphSpec s;
        s.k = 6;
        s.cubics.push_back({ps[0].first, ps[0].second, -1});
        s.cubics.push_back({ps[1].first, ps[1].second, -2});
        s.cubics.push_back({ps[2].first, ps[2].second, -3});
        s.cubics.push_back({-1, -2, -3});
        s.internal_edges = 3;
        specs.push_back(std::move(s));
    }
}

static std::string graph_key(const GraphSpec& s) {
    vector<int> perm(s.internal_edges);
    std::iota(perm.begin(), perm.end(), 0);
    std::string best;
    bool have = false;
    do {
        vector<pair<int, int>> metrics = s.metrics;
        for (auto& pr : metrics) {
            if (pr.first > pr.second) std::swap(pr.first, pr.second);
        }
        std::sort(metrics.begin(), metrics.end());
        vector<array<int, 3>> cubics;
        cubics.reserve(s.cubics.size());
        for (auto cu : s.cubics) {
            for (int& leg : cu) {
                if (leg < 0) leg = -1 - perm[-1 - leg];
            }
            std::sort(cu.begin(), cu.end());
            cubics.push_back(cu);
        }
        std::sort(cubics.begin(), cubics.end());
        std::ostringstream out;
        out << "k" << s.k << "|m";
        for (auto pr : metrics) out << pr.first << "," << pr.second << ";";
        out << "|c";
        for (auto cu : cubics) out << cu[0] << "," << cu[1] << "," << cu[2] << ";";
        std::string key = out.str();
        if (!have || key < best) {
            have = true;
            best = std::move(key);
        }
    } while (std::next_permutation(perm.begin(), perm.end()));
    return best;
}

static vector<GraphSpec> l6_graph_specs() {
    vector<GraphSpec> specs;
    add_l6_metric_pairings(specs);
    add_l6_two_cubics(specs);
    add_l6_metric_connected_cubics(specs);
    add_l6_cubic_paths(specs);
    add_l6_cubic_stars(specs);
    std::map<std::string, GraphSpec> uniq;
    for (auto& s : specs) uniq.emplace(graph_key(s), s);
    vector<GraphSpec> out;
    out.reserve(uniq.size());
    for (auto& kv : uniq) out.push_back(std::move(kv.second));
    return out;
}

static vector<GraphSpec> l6_core_graph_specs() {
    vector<GraphSpec> specs;
    add_l6_metric_pairings(specs);
    add_l6_two_cubics(specs);
    add_l6_metric_connected_cubics(specs);
    return specs;
}

static vector<GraphSpec> l6_core_star_graph_specs() {
    vector<GraphSpec> specs = l6_core_graph_specs();
    add_l6_cubic_stars(specs);
    return specs;
}

static vector<GraphSpec> l6_ordered_graph_specs() {
    vector<GraphSpec> specs = l6_core_star_graph_specs();
    add_l6_cubic_paths(specs);
    return specs;
}

static int run_l6_incremental() {
    vector<GraphSpec> specs = l6_ordered_graph_specs();
    vector<int> selected;
    vector<vector<int>> gram_sel;
    int selected_rank = 0;
    std::cout << "F4_L6_INC_TASKS prime " << MOD
              << " candidates " << specs.size();
#ifdef _OPENMP
    std::cout << " threads " << omp_get_max_threads();
#else
    std::cout << " threads 1";
#endif
    std::cout << std::endl;

    for (int idx = 0; idx < (int)specs.size() && selected_rank < 70; ++idx) {
        std::cout << "F4_L6_INC_CANDIDATE index " << idx
                  << " selected " << selected.size()
                  << " rank " << selected_rank << std::endl;
        vector<int> col(selected.size());
        std::atomic<int> done(0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int j = 0; j < (int)selected.size(); ++j) {
            col[j] = graph_gram_entry(specs[idx], idx, specs[selected[j]], selected[j]);
            int d = ++done;
            if (d % 10 == 0 || d == (int)selected.size()) {
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    std::cout << "F4_L6_INC_PROGRESS candidate " << idx
                              << " entries " << d << "/" << selected.size()
                              << std::endl;
                }
            }
        }
        int self = graph_gram_entry(specs[idx], idx, specs[idx], idx);
        vector<vector<int>> cand = gram_sel;
        for (auto& row : cand) row.push_back(0);
        cand.push_back(vector<int>(selected.size() + 1, 0));
        for (int j = 0; j < (int)selected.size(); ++j) {
            cand[j].back() = col[j];
            cand.back()[j] = col[j];
        }
        cand.back().back() = self;
        int rr = rank_mod(cand);
        if (rr > selected_rank) {
            selected.push_back(idx);
            gram_sel = std::move(cand);
            selected_rank = rr;
            std::cout << "F4_L6_INC_SELECT index " << idx
                      << " selected " << selected.size()
                      << " rank " << selected_rank << std::endl;
        }
    }
    std::cout << "F4_L6_INC_RANK prime " << MOD
              << " candidates " << specs.size()
              << " selected " << selected.size()
              << " rank " << selected_rank
              << " selected_indices";
    for (int idx : selected) std::cout << " " << idx;
    std::cout << std::endl;
    return selected_rank == 70 ? 0 : 8;
}

static vector<int> l6_selected_indices() {
    return {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
        10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
        20, 21, 22, 23, 24, 25, 26, 28, 29, 31,
        32, 34, 35, 37, 38, 40, 41, 43, 44, 46,
        47, 49, 50, 52, 53, 55, 56, 58, 59, 61,
        62, 64, 65, 67, 68, 70, 71, 73, 74, 77,
        85, 86, 87, 88, 89, 92, 103, 104, 107, 122,
    };
}

static int run_l6_supports() {
    vector<GraphSpec> specs = l6_ordered_graph_specs();
    vector<int> selected = l6_selected_indices();
    uint64_t total = 0;
    std::cout << "F4_L6_SUPPORTS_BEGIN selected " << selected.size() << std::endl;
    for (int pos = 0; pos < (int)selected.size(); ++pos) {
        int idx = selected[pos];
        Factor f = graph_unit_tensor_factor(specs[idx], idx + 1);
        total += f.entries.size();
        std::cout << "F4_L6_SUPPORT tensor_pos " << pos
                  << " spec " << idx
                  << " support " << f.entries.size()
                  << " vars";
        for (int v : f.vars) std::cout << " " << v;
        std::cout << std::endl;
    }
    std::cout << "F4_L6_SUPPORTS_DONE total " << total << std::endl;
    return 0;
}

static int run_m123_selfcheck() {
    vector<GraphSpec> specs = l6_ordered_graph_specs();
    vector<int> selected = l6_selected_indices();
    vector<array<int, 4>> cases = {
        array<int, 4>{0, 0, 0, 0},
        array<int, 4>{0, 1, 2, 3},
        array<int, 4>{15, 16, 17, 18},
    };
    int done = 0;
    std::cout << "F4_M123_SELFCHECK_BEGIN cases " << cases.size() << std::endl;
    for (auto test : cases) {
        int a = test[0], b = test[1], c = test[2], d = test[3];
        if (a >= (int)selected.size() || b >= (int)selected.size() ||
            c >= (int)selected.size() || d >= (int)selected.size()) continue;
        std::cout << "F4_M123_SELFCHECK_CASE positions "
                  << a << " " << b << " " << c << " " << d << std::endl;
        int val = m123_unit_entry(specs[selected[a]], selected[a],
                                  specs[selected[b]], selected[b],
                                  specs[selected[c]], selected[c],
                                  specs[selected[d]], selected[d]);
        std::cout << "F4_M123_SELFCHECK_VALUE positions "
                  << a << " " << b << " " << c << " " << d
                  << " value " << val << std::endl;
        ++done;
    }
    std::cout << "F4_M123_SELFCHECK_DONE cases " << done << std::endl;
    return 0;
}

static int run_l6_unit_inverse_stats() {
    vector<GraphSpec> specs = l6_ordered_graph_specs();
    vector<int> selected = l6_selected_indices();
    int S = selected.size();
    vector<vector<int>> gram(S, vector<int>(S, 0));
    vector<pair<int, int>> cells;
    for (int i = 0; i < S; ++i) {
        for (int j = i; j < S; ++j) cells.push_back({i, j});
    }
    std::cout << "F4_L6_UNIT_GRAM_TASKS prime " << MOD
              << " selected " << S
              << " cells " << cells.size();
#ifdef _OPENMP
    std::cout << " threads " << omp_get_max_threads();
#else
    std::cout << " threads 1";
#endif
    std::cout << std::endl;
    std::atomic<int> done(0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int idx = 0; idx < (int)cells.size(); ++idx) {
        int i = cells[idx].first;
        int j = cells[idx].second;
        int val = graph_unit_gram_entry(specs[selected[i]], selected[i],
                                        specs[selected[j]], selected[j]);
        gram[i][j] = val;
        gram[j][i] = val;
        int d = ++done;
        if (d % 100 == 0 || d == (int)cells.size()) {
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                std::cout << "F4_L6_UNIT_GRAM_PROGRESS cells "
                          << d << "/" << cells.size() << std::endl;
            }
        }
    }
    int rank = rank_mod(gram);
    vector<vector<int>> inv = inverse_mod_matrix(gram);
    int gram_nnz = 0;
    int inv_nnz = 0;
    for (int i = 0; i < S; ++i) {
        for (int j = 0; j < S; ++j) {
            if (gram[i][j]) ++gram_nnz;
            if (inv[i][j]) ++inv_nnz;
        }
    }
    std::cout << "F4_L6_UNIT_INV_STATS prime " << MOD
              << " selected " << S
              << " rank " << rank
              << " gram_nnz " << gram_nnz
              << " inv_nnz " << inv_nnz << std::endl;
    return rank == S ? 0 : 9;
}

static int run_l6_rank(int l6_mode) {
    vector<GraphSpec> specs;
    if (l6_mode == 1) specs = l6_core_graph_specs();
    else if (l6_mode == 2) specs = l6_core_star_graph_specs();
    else specs = l6_graph_specs();
    vector<vector<int>> gram(specs.size(), vector<int>(specs.size()));
    std::cout << "F4_L6_MOD_TASKS prime " << MOD
              << " mode " << (l6_mode == 1 ? "core" : (l6_mode == 2 ? "core_star" : "all"))
              << " tensors " << specs.size();
#ifdef _OPENMP
    std::cout << " threads " << omp_get_max_threads();
#else
    std::cout << " threads 1";
#endif
    std::cout << std::endl;
    vector<pair<int, int>> cells;
    for (int i = 0; i < (int)specs.size(); ++i) {
        for (int j = i; j < (int)specs.size(); ++j) cells.push_back({i, j});
    }
    std::atomic<int> done(0);
    int total_cells = int(cells.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int idx = 0; idx < total_cells; ++idx) {
        int i = cells[idx].first;
        int j = cells[idx].second;
        int val = graph_gram_entry(specs[i], i, specs[j], j);
        gram[i][j] = val;
        gram[j][i] = val;
        int d = ++done;
        if (d % 500 == 0 || d == total_cells) {
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                std::cout << "F4_L6_MOD_PROGRESS cells "
                          << d << "/" << total_cells
                          << " last_cell " << i << "," << j << std::endl;
            }
        }
    }
    int rank = rank_mod(gram);
    std::cout << "F4_L6_MOD_RANK prime " << MOD
              << " cubic_nonzero " << cubic_entries.size()
              << " tensors " << specs.size()
              << " rank " << rank << std::endl;
    return rank == 70 ? 0 : 7;
}

static int unit_gram_entry(const TensorSpec& a, int ia, const TensorSpec& b, int ib) {
    array<int, 5> ext = {0, 1, 2, 3, 4};
    vector<Factor> f = tensor_unit_factors(a, ext, ia + 1);
    vector<Factor> g = tensor_unit_factors(b, ext, ib + 51);
    f.insert(f.end(), g.begin(), g.end());
    return eliminate_value(std::move(f));
}

static int rank_mod(vector<vector<int>> A) {
    int n = A.size(), m = A.empty() ? 0 : A[0].size();
    int r = 0;
    for (int c = 0; c < m && r < n; ++c) {
        int piv = -1;
        for (int i = r; i < n; ++i) if (A[i][c]) { piv = i; break; }
        if (piv < 0) continue;
        std::swap(A[r], A[piv]);
        int invp = inv_mod(A[r][c]);
        for (int j = c; j < m; ++j) A[r][j] = mul_mod(A[r][j], invp);
        for (int i = 0; i < n; ++i) if (i != r && A[i][c]) {
            int fac = A[i][c];
            for (int j = c; j < m; ++j) {
                A[i][j] = sub_mod(A[i][j], mul_mod(fac, A[r][j]));
            }
        }
        ++r;
    }
    return r;
}

static vector<vector<int>> inverse_mod_matrix(vector<vector<int>> A) {
    int n = A.size();
    vector<vector<int>> aug(n, vector<int>(2 * n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) aug[i][j] = A[i][j];
        aug[i][n + i] = 1;
    }
    for (int c = 0; c < n; ++c) {
        int piv = -1;
        for (int i = c; i < n; ++i) if (aug[i][c]) { piv = i; break; }
        if (piv < 0) {
            std::cerr << "singular selected Gram at column " << c << "\n";
            std::exit(3);
        }
        std::swap(aug[c], aug[piv]);
        int invp = inv_mod(aug[c][c]);
        for (int j = c; j < 2 * n; ++j) aug[c][j] = mul_mod(aug[c][j], invp);
        for (int i = 0; i < n; ++i) if (i != c && aug[i][c]) {
            int fac = aug[i][c];
            for (int j = c; j < 2 * n; ++j) aug[i][j] = sub_mod(aug[i][j], mul_mod(fac, aug[c][j]));
        }
    }
    vector<vector<int>> inv(n, vector<int>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) inv[i][j] = aug[i][n + j];
    }
    return inv;
}

static int m13_unit_entry(const TensorSpec& a, int ia, const TensorSpec& b, int ib,
                          const TensorSpec& c, int ic, const TensorSpec& dspec, int id) {
    enum Label {
        I = 0, J = 1, J1 = 2, K1 = 3, I1 = 4, L1 = 5,
        J2 = 6, K2 = 7, I2 = 8, L2 = 9
    };
    array<int, 5> A_rows = {I, J1, K1, J2, K2};
    array<int, 5> A_cols = {I, I1, L1, I2, L2};
    array<int, 5> B_rows = {J, J1, I1, J2, I2};
    array<int, 5> B_cols = {J, K1, L1, K2, L2};
    vector<Factor> f = tensor_unit_factors(a, A_rows, ia + 1);
    vector<Factor> g = tensor_unit_factors(b, A_cols, ib + 51);
    vector<Factor> h = tensor_unit_factors(c, B_rows, ic + 101);
    vector<Factor> q = tensor_unit_factors(dspec, B_cols, id + 151);
    f.insert(f.end(), g.begin(), g.end());
    f.insert(f.end(), h.begin(), h.end());
    f.insert(f.end(), q.begin(), q.end());
    return eliminate_value(std::move(f));
}

static int m123_unit_entry(const GraphSpec& a, int ia, const GraphSpec& b, int ib,
                           const GraphSpec& c, int ic, const GraphSpec& dspec, int id) {
    enum Label {
        I = 0, J = 1, P = 2, Q = 3,
        J1 = 4, K1 = 5, I1 = 6, L1 = 7,
        J2 = 8, K2 = 9, I2 = 10, L2 = 11
    };
    array<int, 6> A_rows = {I, P, J1, K1, J2, K2};
    array<int, 6> A_cols = {I, Q, I1, L1, I2, L2};
    array<int, 6> B_rows = {J, P, J1, I1, J2, I2};
    array<int, 6> B_cols = {J, Q, K1, L1, K2, L2};
    vector<Factor> f = graph_unit_factors(a, A_rows, ia + 1);
    vector<Factor> g = graph_unit_factors(b, A_cols, ib + 251);
    vector<Factor> h = graph_unit_factors(c, B_rows, ic + 501);
    vector<Factor> q = graph_unit_factors(dspec, B_cols, id + 751);
    f.insert(f.end(), g.begin(), g.end());
    f.insert(f.end(), h.begin(), h.end());
    f.insert(f.end(), q.begin(), q.end());
    return eliminate_value(std::move(f));
}

struct DSU {
    vector<int> parent;
    explicit DSU(int n) : parent(n) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }
    void unite(int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra != rb) parent[rb] = ra;
    }
};

static void tensor_hyperedges(const TensorSpec& s, const array<int, 5>& ext, int prefix,
                              vector<array<int, 3>>& cubics,
                              vector<pair<int, int>>& deltas) {
    auto e = [&](int r) { return ext[r]; };
    int base = prefix * 100 + 20;
    if (s.kind == 0) {
        deltas.push_back({e(s.pair0.first), e(s.pair0.second)});
        cubics.push_back({e(s.tri[0]), e(s.tri[1]), e(s.tri[2])});
    } else {
        int x = base;
        int y = base + 1;
        cubics.push_back({e(s.pair1.first), e(s.pair1.second), x});
        cubics.push_back({e(s.single), x, y});
        cubics.push_back({e(s.pair2.first), e(s.pair2.second), y});
    }
}

static std::string m13_graph_key(const TensorSpec& a, int ia, const TensorSpec& b, int ib,
                                 const TensorSpec& c, int ic, const TensorSpec& dspec, int id) {
    enum Label {
        I = 0, J = 1, J1 = 2, K1 = 3, I1 = 4, L1 = 5,
        J2 = 6, K2 = 7, I2 = 8, L2 = 9
    };
    array<int, 5> A_rows = {I, J1, K1, J2, K2};
    array<int, 5> A_cols = {I, I1, L1, I2, L2};
    array<int, 5> B_rows = {J, J1, I1, J2, I2};
    array<int, 5> B_cols = {J, K1, L1, K2, L2};
    vector<array<int, 3>> cubics;
    vector<pair<int, int>> deltas;
    tensor_hyperedges(a, A_rows, ia + 1, cubics, deltas);
    tensor_hyperedges(b, A_cols, ib + 51, cubics, deltas);
    tensor_hyperedges(c, B_rows, ic + 101, cubics, deltas);
    tensor_hyperedges(dspec, B_cols, id + 151, cubics, deltas);

    int max_var = 0;
    for (auto e : cubics) for (int x : e) max_var = std::max(max_var, x);
    for (auto e : deltas) max_var = std::max(max_var, std::max(e.first, e.second));
    DSU dsu(max_var + 1);
    for (auto e : deltas) dsu.unite(e.first, e.second);

    std::map<int, int> ext_mask;
    for (int e = 0; e < 10; ++e) ext_mask[dsu.find(e)] |= (1 << e);
    std::map<int, vector<pair<int, int>>> internal_occ;
    vector<array<int, 3>> roots;
    roots.reserve(cubics.size());
    for (int ci = 0; ci < (int)cubics.size(); ++ci) {
        array<int, 3> r = {dsu.find(cubics[ci][0]), dsu.find(cubics[ci][1]), dsu.find(cubics[ci][2])};
        roots.push_back(r);
        for (int p = 0; p < 3; ++p) {
            if (!ext_mask[r[p]]) internal_occ[r[p]].push_back({ci, p});
        }
    }
    vector<pair<vector<pair<int, int>>, int>> internal_sorted;
    for (auto& kv : internal_occ) internal_sorted.push_back({kv.second, kv.first});
    std::sort(internal_sorted.begin(), internal_sorted.end());
    std::map<int, int> internal_label;
    for (int i = 0; i < (int)internal_sorted.size(); ++i) {
        internal_label[internal_sorted[i].second] = 2000 + i;
    }
    vector<array<int, 3>> final_edges;
    for (auto r : roots) {
        array<int, 3> edge{};
        for (int p = 0; p < 3; ++p) {
            int mask = ext_mask[r[p]];
            edge[p] = mask ? 1000 + mask : internal_label[r[p]];
        }
        std::sort(edge.begin(), edge.end());
        final_edges.push_back(edge);
    }
    std::sort(final_edges.begin(), final_edges.end());
    std::ostringstream out;
    for (auto e : final_edges) out << e[0] << ',' << e[1] << ',' << e[2] << ';';
    return out.str();
}

struct TensorEntry5 {
    array<int, 5> idx;
    int val;
};

struct SparseMat {
    vector<vector<pair<int, int>>> rows;
    vector<vector<pair<int, int>>> cols;
    int nnz = 0;
};

struct DotScratch {
    vector<int> acc;
    vector<int> touched;
    DotScratch() : acc(MAT_DIM * MAT_DIM, 0) {
        touched.reserve(MAT_DIM * 8);
    }
};

static uint32_t pack5(const array<int, 5>& a) {
    uint32_t out = 0;
    for (int x : a) out = out * DIM + uint32_t(x);
    return out;
}

static vector<TensorEntry5> build_unit_tensor_map(const TensorSpec& s) {
    std::unordered_map<uint32_t, int> acc;
    auto add = [&](const array<int, 5>& idx, int val) {
        uint32_t key = pack5(idx);
        auto it = acc.find(key);
        if (it == acc.end()) acc.emplace(key, val);
        else it->second = add_mod(it->second, val);
    };
    if (s.kind == 0) {
        for (int x = 0; x < DIM; ++x) {
            for (size_t t = 0; t < cubic_entries.size(); ++t) {
                array<int, 5> idx{};
                idx[s.pair0.first] = x;
                idx[s.pair0.second] = x;
                idx[s.tri[0]] = cubic_entries[t][0];
                idx[s.tri[1]] = cubic_entries[t][1];
                idx[s.tri[2]] = cubic_entries[t][2];
                add(idx, cubic_values_unit[t]);
            }
        }
    } else {
        vector<vector<int>> mid_by_x(DIM), last_by_y(DIM);
        for (int t = 0; t < (int)cubic_entries.size(); ++t) {
            mid_by_x[cubic_entries[t][1]].push_back(t);
            last_by_y[cubic_entries[t][2]].push_back(t);
        }
        for (int t1 = 0; t1 < (int)cubic_entries.size(); ++t1) {
            int x = cubic_entries[t1][2];
            for (int t2 : mid_by_x[x]) {
                int y = cubic_entries[t2][2];
                int coeff12 = mul_mod(cubic_values_unit[t1], cubic_values_unit[t2]);
                for (int t3 : last_by_y[y]) {
                    array<int, 5> idx{};
                    idx[s.pair1.first] = cubic_entries[t1][0];
                    idx[s.pair1.second] = cubic_entries[t1][1];
                    idx[s.single] = cubic_entries[t2][0];
                    idx[s.pair2.first] = cubic_entries[t3][0];
                    idx[s.pair2.second] = cubic_entries[t3][1];
                    add(idx, mul_mod(coeff12, cubic_values_unit[t3]));
                }
            }
        }
    }
    vector<TensorEntry5> out;
    out.reserve(acc.size());
    for (auto& kv : acc) {
        if (!kv.second) continue;
        uint32_t key = kv.first;
        array<int, 5> idx{};
        for (int i = 4; i >= 0; --i) {
            idx[i] = key % DIM;
            key /= DIM;
        }
        out.push_back({idx, kv.second});
    }
    return out;
}

static SparseMat make_sparse_mat(vector<tuple<int, int, int>> entries) {
    std::sort(entries.begin(), entries.end());
    SparseMat out;
    out.rows.assign(MAT_DIM, {});
    out.cols.assign(MAT_DIM, {});
    size_t p = 0;
    while (p < entries.size()) {
        int r = std::get<0>(entries[p]);
        int c = std::get<1>(entries[p]);
        int val = 0;
        while (p < entries.size() && std::get<0>(entries[p]) == r && std::get<1>(entries[p]) == c) {
            val = add_mod(val, std::get<2>(entries[p]));
            ++p;
        }
        if (!val) continue;
        out.rows[r].push_back({c, val});
        out.cols[c].push_back({r, val});
        ++out.nnz;
    }
    return out;
}

static SparseMat sparse_linear_combo(const vector<SparseMat>& mats, const vector<int>& coeffs) {
    size_t reserve = 0;
    for (int i = 0; i < (int)mats.size(); ++i) if (coeffs[i]) reserve += mats[i].nnz;
    vector<tuple<int, int, int>> entries;
    entries.reserve(reserve);
    for (int i = 0; i < (int)mats.size(); ++i) {
        int coeff = coeffs[i];
        if (!coeff) continue;
        for (int r = 0; r < MAT_DIM; ++r) {
            for (auto [c, val] : mats[i].rows[r]) {
                entries.push_back({r, c, mul_mod(coeff, val)});
            }
        }
    }
    return make_sparse_mat(std::move(entries));
}

static int64_t estimate_row_pair_ops(const SparseMat& A, const SparseMat& C) {
    int64_t ops = 0;
    for (int r = 0; r < MAT_DIM; ++r) {
        ops += int64_t(A.rows[r].size()) * int64_t(C.rows[r].size());
    }
    return ops;
}

static int64_t estimate_col_pair_ops(const SparseMat& D, const SparseMat& B) {
    int64_t ops = 0;
    for (int q = 0; q < MAT_DIM; ++q) {
        ops += int64_t(D.cols[q].size()) * int64_t(B.cols[q].size());
    }
    return ops;
}

static void scratch_add(DotScratch& scratch, int key, int val) {
    if (!val) return;
    if (scratch.acc[key] == 0) scratch.touched.push_back(key);
    scratch.acc[key] = add_mod(scratch.acc[key], val);
}

static void scratch_clear(DotScratch& scratch) {
    for (int key : scratch.touched) scratch.acc[key] = 0;
    scratch.touched.clear();
}

static int product_dot_build_left(const SparseMat& A, const SparseMat& D,
                                  const SparseMat& C, const SparseMat& B,
                                  DotScratch& scratch) {
    for (int r = 0; r < MAT_DIM; ++r) {
        const auto& ar = A.rows[r];
        const auto& cr = C.rows[r];
        if (ar.empty() || cr.empty()) continue;
        for (auto [p, aval] : ar) {
            int base = p * MAT_DIM;
            for (auto [s, cval] : cr) {
                scratch_add(scratch, base + s, mul_mod(aval, cval));
            }
        }
    }
    int total = 0;
    for (int q = 0; q < MAT_DIM; ++q) {
        const auto& dq = D.cols[q];
        const auto& bq = B.cols[q];
        if (dq.empty() || bq.empty()) continue;
        for (auto [p, dval] : dq) {
            int base = p * MAT_DIM;
            for (auto [s, bval] : bq) {
                int left = scratch.acc[base + s];
                if (left) total = add_mod(total, mul_mod(left, mul_mod(dval, bval)));
            }
        }
    }
    scratch_clear(scratch);
    return total;
}

static int product_dot_build_right(const SparseMat& A, const SparseMat& D,
                                   const SparseMat& C, const SparseMat& B,
                                   DotScratch& scratch) {
    for (int q = 0; q < MAT_DIM; ++q) {
        const auto& dq = D.cols[q];
        const auto& bq = B.cols[q];
        if (dq.empty() || bq.empty()) continue;
        for (auto [p, dval] : dq) {
            int base = p * MAT_DIM;
            for (auto [s, bval] : bq) {
                scratch_add(scratch, base + s, mul_mod(dval, bval));
            }
        }
    }
    int total = 0;
    for (int r = 0; r < MAT_DIM; ++r) {
        const auto& ar = A.rows[r];
        const auto& cr = C.rows[r];
        if (ar.empty() || cr.empty()) continue;
        for (auto [p, aval] : ar) {
            int base = p * MAT_DIM;
            for (auto [s, cval] : cr) {
                int right = scratch.acc[base + s];
                if (right) total = add_mod(total, mul_mod(right, mul_mod(aval, cval)));
            }
        }
    }
    scratch_clear(scratch);
    return total;
}

static int product_dot(const SparseMat& A, const SparseMat& D,
                       const SparseMat& C, const SparseMat& B,
                       DotScratch& scratch) {
    int64_t left_ops = estimate_row_pair_ops(A, C);
    int64_t right_ops = estimate_col_pair_ops(D, B);
    if (left_ops <= right_ops) return product_dot_build_left(A, D, C, B, scratch);
    return product_dot_build_right(A, D, C, B, scratch);
}

static int compute_m13_matrix(const vector<TensorSpec>& specs,
                              const vector<int>& selected,
                              const vector<vector<int>>& gram_inv,
                              bool selfcheck_only) {
    int S = selected.size();
    vector<vector<SparseMat>> base(DIM, vector<SparseMat>(S));
    uint64_t base_nnz = 0;
    for (int pos = 0; pos < S; ++pos) {
        auto entries5 = build_unit_tensor_map(specs[selected[pos]]);
        vector<vector<tuple<int, int, int>>> by_first(DIM);
        for (const auto& e : entries5) {
            int first = e.idx[0];
            int row = e.idx[1] * DIM + e.idx[3];
            int col = e.idx[2] * DIM + e.idx[4];
            by_first[first].push_back({row, col, e.val});
        }
        uint64_t tensor_nnz = 0;
        for (int x = 0; x < DIM; ++x) {
            base[x][pos] = make_sparse_mat(std::move(by_first[x]));
            tensor_nnz += base[x][pos].nnz;
        }
        base_nnz += tensor_nnz;
        std::cout << "F4_M13_MATRIX_BASE tensor_pos " << pos
                  << " spec " << selected[pos]
                  << " nnz " << tensor_nnz << std::endl;
    }

    vector<array<int, 4>> selfcheck_cases = {
        array<int, 4>{0, 0, 0, 0},
        array<int, 4>{0, 1, 2, 3},
    };
    int selfcheck_done = 0;
    for (auto test : selfcheck_cases) {
        int a = test[0], b = test[1], c = test[2], d = test[3];
        if (a >= S || b >= S || c >= S || d >= S) continue;
        std::cout << "F4_M13_MATRIX_SELFCHECK_BEGIN positions "
                  << a << " " << b << " " << c << " " << d << std::endl;
        int matrix_val = 0;
        DotScratch scratch;
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                int val = product_dot(base[i][a], base[j][d],
                                      base[j][c], base[i][b], scratch);
                matrix_val = add_mod(matrix_val, val);
            }
        }
        int factor_val = m13_unit_entry(specs[selected[a]], selected[a],
                                        specs[selected[b]], selected[b],
                                        specs[selected[c]], selected[c],
                                        specs[selected[d]], selected[d]);
        std::cout << "F4_M13_MATRIX_SELFCHECK positions "
                  << a << " " << b << " " << c << " " << d
                  << " matrix " << matrix_val
                  << " factor " << factor_val << std::endl;
        if (matrix_val != factor_val) return -1;
        ++selfcheck_done;
    }
    if (selfcheck_only) {
        std::cout << "F4_M13_MATRIX_SELFCHECK_OK cases "
                  << selfcheck_done << std::endl;
        return 0;
    }

    vector<vector<SparseMat>> bbar(DIM, vector<SparseMat>(S));
    vector<vector<SparseMat>> cbar(DIM, vector<SparseMat>(S));
    uint64_t bbar_nnz = 0;
    uint64_t cbar_nnz = 0;
    vector<int> coeffs(S);
    for (int x = 0; x < DIM; ++x) {
        for (int a = 0; a < S; ++a) {
            for (int b = 0; b < S; ++b) coeffs[b] = gram_inv[a][b];
            bbar[x][a] = sparse_linear_combo(base[x], coeffs);
            bbar_nnz += bbar[x][a].nnz;
        }
        for (int d = 0; d < S; ++d) {
            for (int c = 0; c < S; ++c) coeffs[c] = gram_inv[c][d];
            cbar[x][d] = sparse_linear_combo(base[x], coeffs);
            cbar_nnz += cbar[x][d].nnz;
        }
        std::cout << "F4_M13_MATRIX_BARS first " << (x + 1) << "/" << DIM << std::endl;
    }

    std::cout << "F4_M13_MATRIX_TASKS first_pairs " << (DIM * DIM)
              << " inner_pairs " << (S * S)
              << " base_nnz " << base_nnz
              << " bbar_nnz " << bbar_nnz
              << " cbar_nnz " << cbar_nnz;
#ifdef _OPENMP
    std::cout << " threads " << omp_get_max_threads();
#else
    std::cout << " threads 1";
#endif
    std::cout << std::endl;

    std::atomic<int> done(0);
    int64_t total64 = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:total64)
#endif
    {
        DotScratch scratch;
        int64_t local64 = 0;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (int pair_idx = 0; pair_idx < DIM * DIM; ++pair_idx) {
            int i = pair_idx / DIM;
            int j = pair_idx % DIM;
            int pair_total = 0;
            for (int a = 0; a < S; ++a) {
                for (int d = 0; d < S; ++d) {
                    int val = product_dot(base[i][a], base[j][d],
                                          cbar[j][d], bbar[i][a], scratch);
                    pair_total = add_mod(pair_total, val);
                }
            }
            local64 += pair_total;
            int d_done = ++done;
            if (d_done % 10 == 0 || d_done == DIM * DIM) {
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    std::cout << "F4_M13_MATRIX_PROGRESS first_pairs "
                              << d_done << "/" << (DIM * DIM) << std::endl;
                }
            }
        }
        total64 += local64 % MOD;
    }
    int m13 = int(total64 % MOD);
    std::cout << "F4_M13_MATRIX_VALUE prime " << MOD
              << " value " << m13
              << " first_pairs " << (DIM * DIM)
              << " inner_pairs " << (S * S) << std::endl;
    return m13;
}

static bool init_unit_cubic_values() {
    vector<int> sqrt_norms(DIM);
    for (int i = 0; i < DIM; ++i) {
        // The supported certificate primes are 3 mod 4 and make 2,6 squares.
        sqrt_norms[i] = pow_mod(norms[i], (MOD + 1) / 4);
        if (mul_mod(sqrt_norms[i], sqrt_norms[i]) != norms[i]) {
            std::cerr << "norm is not a square modulo prime at " << i << "\n";
            return false;
        }
    }
    cubic_values_unit.clear();
    cubic_values_unit.reserve(cubic_values.size());
    for (size_t t = 0; t < cubic_values.size(); ++t) {
        int i = cubic_entries[t][0], j = cubic_entries[t][1], k = cubic_entries[t][2];
        int den = mul_mod(mul_mod(sqrt_norms[i], sqrt_norms[j]), sqrt_norms[k]);
        cubic_values_unit.push_back(mul_mod(cubic_values[t], inv_mod(den)));
    }
    return true;
}

int main(int argc, char** argv) {
    bool do_l6_rank = false;
    bool do_l6_incremental = false;
    bool do_l6_supports = false;
    bool do_m123_selfcheck = false;
    bool do_l6_unit_inverse_stats = false;
    int l6_mode = 0;
    bool do_m13 = false;
    bool do_m13_factor = false;
    bool do_m13_selfcheck = false;
    bool do_m13_supports = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--l6-rank") do_l6_rank = true;
        if (std::string(argv[i]) == "--l6-rank-core") {
            do_l6_rank = true;
            l6_mode = 1;
        }
        if (std::string(argv[i]) == "--l6-rank-core-star") {
            do_l6_rank = true;
            l6_mode = 2;
        }
        if (std::string(argv[i]) == "--l6-rank-incremental") do_l6_incremental = true;
        if (std::string(argv[i]) == "--l6-supports") do_l6_supports = true;
        if (std::string(argv[i]) == "--l6-unit-inv-stats") do_l6_unit_inverse_stats = true;
        if (std::string(argv[i]) == "--m123-selfcheck") do_m123_selfcheck = true;
        if (std::string(argv[i]) == "--m13") do_m13 = true;
        if (std::string(argv[i]) == "--m13-factor") {
            do_m13 = true;
            do_m13_factor = true;
        }
        if (std::string(argv[i]) == "--m13-selfcheck") {
            do_m13 = true;
            do_m13_selfcheck = true;
        }
        if (std::string(argv[i]) == "--m13-supports") {
            do_m13 = true;
            do_m13_supports = true;
        }
    }
    init_octonions();
    vector<Mat> basis;
    basis.push_back(diag_mat(1, -1, 0));
    basis.push_back(diag_mat(1, 1, -2));
    for (auto pq : vector<pair<int, int>>{{0, 1}, {0, 2}, {1, 2}}) {
        for (int m = 0; m < 8; ++m) basis.push_back(off_mat(pq.first, pq.second, m));
    }
    norms.assign(DIM, 0);
    for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) {
            int val = trace(jordan(basis[i], basis[j]));
            if (i != j && val != 0) {
                std::cerr << "basis not orthogonal at " << i << "," << j << "\n";
                return 2;
            }
            if (i == j) norms[i] = val;
        }
    }
    qnorms.resize(DIM);
    for (int i = 0; i < DIM; ++i) qnorms[i] = inv_mod(norms[i]);

    for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) {
            Mat Jij = jordan(basis[i], basis[j]);
            for (int k = 0; k < DIM; ++k) {
                int val = trace(jordan(Jij, basis[k]));
                if (val) {
                    cubic_entries.push_back({i, j, k});
                    cubic_values.push_back(val);
                }
            }
        }
    }

    if (do_l6_incremental) return run_l6_incremental();
    if (do_l6_rank) return run_l6_rank(l6_mode);
    if (do_l6_supports) {
        if (!init_unit_cubic_values()) return 4;
        return run_l6_supports();
    }
    if (do_l6_unit_inverse_stats) {
        if (!init_unit_cubic_values()) return 4;
        return run_l6_unit_inverse_stats();
    }
    if (do_m123_selfcheck) {
        if (!init_unit_cubic_values()) return 4;
        return run_m123_selfcheck();
    }

    vector<TensorSpec> specs;
    for (int a = 0; a < 5; ++a) {
        for (int b = a + 1; b < 5; ++b) {
            array<int, 3> tri{};
            int t = 0;
            for (int x = 0; x < 5; ++x) if (x != a && x != b) tri[t++] = x;
            TensorSpec s{};
            s.kind = 0;
            s.pair0 = {a, b};
            s.tri = tri;
            specs.push_back(s);
        }
    }
    for (int s0 = 0; s0 < 5; ++s0) {
        vector<int> rem;
        for (int x = 0; x < 5; ++x) if (x != s0) rem.push_back(x);
        for (auto ps : all_pairings_rec(rem)) {
            TensorSpec s{};
            s.kind = 1;
            s.single = s0;
            s.pair1 = ps[0];
            s.pair2 = ps[1];
            specs.push_back(s);
        }
    }

    vector<vector<int>> gram(specs.size(), vector<int>(specs.size()));
    for (int i = 0; i < (int)specs.size(); ++i) {
        std::cout << "F4_L5_MOD_ROW " << i << "/" << specs.size() << std::endl;
        for (int j = 0; j < (int)specs.size(); ++j) gram[i][j] = gram_entry(specs[i], i, specs[j], j);
    }
    int rank = rank_mod(gram);
    std::cout << "F4_L5_MOD_RANK_OK prime " << MOD
              << " cubic_nonzero " << cubic_entries.size()
              << " tensors " << specs.size()
              << " rank " << rank << std::endl;
    if (rank != 15) return 1;
    if (!do_m13) return 0;

    if (!init_unit_cubic_values()) return 4;

    vector<vector<int>> unit_gram(specs.size(), vector<int>(specs.size()));
    for (int i = 0; i < (int)specs.size(); ++i) {
        std::cout << "F4_L5_UNIT_ROW " << i << "/" << specs.size() << std::endl;
        for (int j = 0; j < (int)specs.size(); ++j) {
            unit_gram[i][j] = unit_gram_entry(specs[i], i, specs[j], j);
        }
    }
    int unit_rank = rank_mod(unit_gram);
    vector<int> selected;
    int selected_rank = 0;
    for (int idx = 0; idx < (int)specs.size() && selected_rank < 15; ++idx) {
        vector<int> cand = selected;
        cand.push_back(idx);
        vector<vector<int>> sub(cand.size(), vector<int>(cand.size()));
        for (int i = 0; i < (int)cand.size(); ++i) {
            for (int j = 0; j < (int)cand.size(); ++j) sub[i][j] = unit_gram[cand[i]][cand[j]];
        }
        int rr = rank_mod(sub);
        if (rr > selected_rank) {
            selected = cand;
            selected_rank = rr;
        }
    }
    std::cout << "F4_L5_UNIT_RANK prime " << MOD
              << " rank " << unit_rank
              << " selected_rank " << selected_rank
              << " selected";
    for (int x : selected) std::cout << " " << x;
    std::cout << std::endl;
    if (selected_rank != 15) return 5;

    if (do_m13_supports) {
        std::cout << "F4_M13_SUPPORTS_BEGIN" << std::endl;
        uint64_t total_support = 0;
        for (int pos = 0; pos < (int)selected.size(); ++pos) {
            auto mp = build_unit_tensor_map(specs[selected[pos]]);
            total_support += mp.size();
            std::cout << "F4_M13_SUPPORT tensor_pos " << pos
                      << " spec " << selected[pos]
                      << " support " << mp.size() << std::endl;
        }
        std::cout << "F4_M13_SUPPORTS_DONE total " << total_support << std::endl;
        return 0;
    }

    vector<vector<int>> gram_sel(selected.size(), vector<int>(selected.size()));
    for (int i = 0; i < (int)selected.size(); ++i) {
        for (int j = 0; j < (int)selected.size(); ++j) {
            gram_sel[i][j] = unit_gram[selected[i]][selected[j]];
        }
    }
    vector<vector<int>> gram_inv = inverse_mod_matrix(gram_sel);
    if (!do_m13_factor) {
        int m13_matrix = compute_m13_matrix(specs, selected, gram_inv, do_m13_selfcheck);
        if (m13_matrix < 0) return 6;
        return 0;
    }
    struct M13Task {
        int ia, ib, ic, id;
        int coeff;
    };
    struct M13Aggregate {
        M13Task task;
        int coeff;
    };
    vector<M13Aggregate> aggregates;
    std::unordered_map<std::string, int> key_to_aggregate;
    int raw_tasks = 0;
    for (int ia = 0; ia < (int)selected.size(); ++ia) {
        for (int ib = 0; ib < (int)selected.size(); ++ib) {
            int cab = gram_inv[ia][ib];
            if (!cab) continue;
            for (int ic = 0; ic < (int)selected.size(); ++ic) {
                for (int id = 0; id < (int)selected.size(); ++id) {
                    int coeff = mul_mod(cab, gram_inv[ic][id]);
                    if (!coeff) continue;
                    ++raw_tasks;
                    std::string key = m13_graph_key(specs[selected[ia]], selected[ia],
                                                    specs[selected[ib]], selected[ib],
                                                    specs[selected[ic]], selected[ic],
                                                    specs[selected[id]], selected[id]);
                    auto it = key_to_aggregate.find(key);
                    if (it == key_to_aggregate.end()) {
                        int pos = aggregates.size();
                        key_to_aggregate.emplace(std::move(key), pos);
                        aggregates.push_back({{ia, ib, ic, id, coeff}, coeff});
                    } else {
                        aggregates[it->second].coeff = add_mod(aggregates[it->second].coeff, coeff);
                    }
                }
            }
        }
    }
    aggregates.erase(std::remove_if(aggregates.begin(), aggregates.end(),
                                    [](const M13Aggregate& a) { return a.coeff == 0; }),
                     aggregates.end());
    std::cout << "F4_M13_MOD_TASKS raw " << raw_tasks
              << " aggregated " << aggregates.size();
#ifdef _OPENMP
    std::cout << " threads " << omp_get_max_threads();
#else
    std::cout << " threads 1";
#endif
    std::cout << std::endl;

    std::atomic<int> done(0);
    int64_t total64 = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) reduction(+:total64)
#endif
    for (int t = 0; t < (int)aggregates.size(); ++t) {
        const auto task = aggregates[t].task;
        int val = m13_unit_entry(specs[selected[task.ia]], selected[task.ia],
                                 specs[selected[task.ib]], selected[task.ib],
                                 specs[selected[task.ic]], selected[task.ic],
                                 specs[selected[task.id]], selected[task.id]);
        total64 += mul_mod(aggregates[t].coeff, val);
        int d = ++done;
        if (d % 100 == 0 || d == (int)aggregates.size()) {
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                std::cout << "F4_M13_MOD_PROGRESS summands " << d
                          << "/" << aggregates.size() << std::endl;
            }
        }
    }
    int m13 = int(total64 % MOD);
    std::cout << "F4_M13_MOD_VALUE prime " << MOD
              << " value " << m13
              << " aggregated_summands " << aggregates.size() << std::endl;
    return 0;
}
