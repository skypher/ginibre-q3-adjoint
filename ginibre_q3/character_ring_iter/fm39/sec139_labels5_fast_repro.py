"""FM-SEC139: the exact FM-SEC135 labels <= 5 box, evaluated by an
incremental C++/GMP table evaluator (FM-MECH48 style), resumable by work item.

Work item = (N, a, n, v): background s^(2a) d^(2(N-a)), n copies of hat S_4,
v copies of the plus label 5.
Each finished item is appended (flush + fsync) to STATE/items_done.txt; a
rerun with the same --state skips finished items.  The final aggregate is
checked against the FM-SEC135 count 227,336,512,420.

Mirror symmetry (default; --no-mirror evaluates everything): y -> -y maps
s <-> d, P -> -P, Z -> Z, hat S_4 -> hat S_4, Z - 1 -> Z - 1, h_2(Z,P) <->
h_2(Z,-P), h_4 <-> W_5.  So the box word with background (a, N-a) and splits
(alpha, q) has the same expectation as the box word with background (N-a, a)
and splits (l-alpha, p-q); the box is invariant under this map for
1 <= a <= N-1 (a = 0 mirrors to the excluded all-plus words a = N).  Hence
only a = 0 and 2a <= N are evaluated; the (N, a, n) entry counts are checked
to be mirror-symmetric, and the totals are reconstructed exactly.

Usage: python3 sec139_labels5_fast_repro.py --threads 32 --state DIR
       [--stride K] (calibration: every K-th item only)
"""
import argparse, os, random, subprocess, sys, time
from datetime import datetime, timezone
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial
from multiprocessing import Pool

def stamp():
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

# ---------------------------------------------------------------- box records
# Same formulas as sec135_labels5_constructor_repro.py (from "C = Q(36369,200)").
rho, delta = Q(6, 7), Q(1, 7)
tau, J = Q(17, 21), Q(352, 1225)
T3, T4, T4p, T5 = Q(4, 5), Q(36, 49), Q(85, 147), Q(3, 4)

@lru_cache(None)
def gb(N, d):
    K = N*delta+d+1
    return max((K+b)*(K+b+1)*
               (Q(1, 4*(b+1)*(b+2))+J*tau**b)/delta**2
               for b in range(5))

@lru_cache(None)
def eb(N, b):
    low = Q(factorial(N), 2**(N+1))
    for j in range(1, N+2):
        low /= b+j
    high = (rho**(N+1)-Q(2, 5)**(N+1))/Q(N+1)*tau**b
    return (low+high)/delta**2

def sup(l, m, n, p):
    if m:
        return 3*T3**l*T4**(m-1)*T4p**n*T5**p
    t = l+p
    if t >= 2:
        take3, take5 = min(2, l), 2-min(2, l)
        return (Q(2)**take3*Q(3)**take5*T3**(l-take3)
                *T5**(p-take5)*T4p**n)
    if t == 1:
        return ((2*T3**(l-1)*T5**p) if l
                else (3*T5**(p-1)))*T4p**n
    return T4p**n

def ways(N, l, m, p):
    total = 0
    for t in range(l+p+1):
        multiplicity = max(0, min(l, t)-max(0, t-p)+1)
        lo = (l+p-t+m+1)//2
        hi = N-max(1, (t+m+1)//2)
        total += multiplicity*max(0, hi-lo+1)
    return total

def profile_task(Hl):
    H, l = Hl
    recs = []
    profiles = 0
    for m in range(H-l+1):
        for n in range(H-l-m+1):
            p = H-l-m-n
            if not p:
                continue
            profiles += 1
            R = max(1, (l+p+1)//2+m)
            raw = Q(2)**l*Q(3)**(m+p)*T4p**n
            d = l+m+2*n+2*p
            for N in range(max(2, R), 150):
                if sup(l, m, n, p)*rho**(N-R)*gb(N, d) < 1:
                    break
                K = N*delta+d+1
                def good(b):
                    return raw*(K+b)*(K+b+1)*eb(N, b) < 1
                hi = 4
                while not good(hi):
                    hi *= 2
                lo = 3
                while hi > lo+1:
                    mid = (hi+lo)//2
                    if good(mid):
                        hi = mid
                    else:
                        lo = mid
                w = ways(N, l, m, p)
                if w:
                    recs.append((N, l, m, n, p, hi, w))
            else:
                raise AssertionError((H, l, m, n, p))
    return profiles, recs

def item_task(chunk):
    # per-(N, a, n, v) entry counts, by direct enumeration of the splits
    out = {}
    for (N, l, m, n, p, B) in chunk:
        for a in range(N):
            A, E = 2*a, 2*(N-a)
            for q in range(p+1):
                lo = max(0, l+p+m-A-q)
                hi = min(l, E-m-q)
                if hi >= lo:
                    key = (N, a, n, p-q)
                    out[key] = out.get(key, 0)+(hi-lo+1)*B
    return out

# ------------------------------------------------- independent direct checks
def add(*terms):
    out = {}
    for c, poly in terms:
        for key, v in poly.items():
            out[key] = out.get(key, 0)+c*v
    return {key: v for key, v in out.items() if v}

def mul(p, q):
    out = {}
    for (i, j), v in p.items():
        for (k, l), w in q.items():
            key = i+k, j+l
            out[key] = out.get(key, 0)+v*w
    return {key: v for key, v in out.items() if v}

def U(n, axis):
    return {((n-2*j, 0) if axis == 0 else (0, n-2*j)):
            (-1)**j*comb(n-j, j) for j in range(n//2+1)}

ONE = {(0, 0): 1}
S = {(1, 0): 1, (0, 1): 1}
D = {(1, 0): 1, (0, 1): -1}
Zp = {(2, 0): 1, (0, 2): 1, (0, 0): -2}
P = {(1, 1): 1}

def hs(k):
    # h_k = (U_(k+1)(x) - U_(k+1)(y))/(x - y)
    out = {}
    for j in range(k+1):
        out = add((1, out), (1, mul(U(j, 0), U(k-j, 1))))
    return out

def Shat(k):
    return add((1, U(k, 0)), (1, U(k, 1)))

def pw(p, e):
    out = ONE
    for _ in range(e):
        out = mul(out, p)
    return out

H2, H3, H4 = hs(2), hs(3), hs(4)
S3, S4, S5 = Shat(3), Shat(4), Shat(5)
ZZ = mul(Zp, Zp)
ZP = mul(Zp, P)
PP = mul(P, P)
assert H2 == add((1, Zp), (1, P))
assert H3 == mul(S, add((1, Zp), (-1, ONE)))
assert S3 == mul(S, add((1, Zp), (-1, P)))
assert S4 == add((1, ZZ), (1, Zp), (-2, PP))
assert H4 == add((1, ZZ), (1, ZP), (-1, PP), (-2, P), (-1, ONE))
assert S5 == mul(S, add((1, ZZ), (-1, ZP), (-1, PP), (2, P), (-1, ONE)))
# minus label k = d h_(k-1): U_k(x) - U_k(y) = (x - y) h_(k-1)
for k in range(1, 6):
    assert add((1, U(k, 0)), (-1, U(k, 1))) == mul(D, hs(k-1))

@lru_cache(None)
def cat(n):
    return comb(2*n, n)//(n+1)

def direct(key):
    A, E, b, al, ga, m, n, q, w = key
    out = ONE
    for g, ex in ((D, E), (S, A-ga-m-w), (Zp, b), (H2, al), (H3, m),
                  (S3, ga), (S4, n), (H4, q), (S5, w)):
        assert ex >= 0
        out = mul(out, pw(g, ex))
    return sum(v*cat(i//2)*cat(j//2)
               for (i, j), v in out.items() if i % 2 == j % 2 == 0)

def moment_check(a, r, b):
    # E[s^(2a) d^(2r) Z^b] with Z = (x^2-1) + (y^2-1): independent of the
    # C++ layer recursion
    N = a+r
    c = [0]*(2*N+1)
    poly = mul(pw(S, 2*a), pw(D, 2*r))
    for (i, j), v in poly.items():
        c[i] = v
    K = 2*N+2*b
    row = [cat(k//2) if k % 2 == 0 else 0 for k in range(K+1)]
    mu = [row[:2*N+1]]
    for i in range(1, b+1):
        row = [row[k+2]-row[k] for k in range(len(row)-2)]
        mu.append(row[:2*N+1])
    total = 0
    for k in range(0, 2*N+1, 2):
        if c[k]:
            total += c[k]*sum(comb(b, i)*mu[i][k]*mu[b-i][2*N-k]
                              for i in range(b+1))
    return total

CPP = r'''// FM-SEC139: exact labels <= 5 box evaluator (incremental tables).
// Word: s^A d^E P^j Z^b moments, core factors applied as in-place
// transforms of the mixed table t[j][b] = E[s^A d^E P^j Z^b].
#include <gmpxx.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
using Z = mpz_class;
struct Rec { int l, m, p, B; };
struct Item { int id, N, a, n, v; unsigned long long expected; int nq; };
struct Dim { int J = -1, L = -1; };
static Dim join(Dim x, Dim y) { return {std::max(x.J, y.J), std::max(x.L, y.L)}; }
static Dim plus2(Dim x) { return x.J < 0 ? x : Dim{x.J + 2, x.L + 2}; }

struct Table {
  std::vector<std::vector<Z>> t;
  Dim d;
  void ensure(Dim x) {
    if ((int)t.size() < x.J + 1) t.resize(x.J + 1);
    for (int j = 0; j <= x.J; ++j)
      if ((int)t[j].size() < x.L - j + 1) t[j].resize(x.L - j + 1);
  }
  void copy_from(const Table& s, Dim x) {
    assert(s.d.J >= x.J && s.d.L >= x.L);
    ensure(x);
    for (int j = 0; j <= x.J; ++j)
      for (int b = 0; b <= x.L - j; ++b) t[j][b] = s.t[j][b];
    d = x;
  }
};

static std::mutex out_mu;
static void say(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
static void say(const char* fmt, ...) {
  std::lock_guard<std::mutex> g(out_mu);
  va_list ap; va_start(ap, fmt); std::vfprintf(stdout, fmt, ap); va_end(ap);
  std::fflush(stdout);
}

static void exact_ui(Z& x, unsigned long d) {
  assert(mpz_divisible_ui_p(x.get_mpz_t(), d));
  mpz_divexact_ui(x.get_mpz_t(), x.get_mpz_t(), d);
}

// in-place transforms; output dims o, input must have dims >= o + shift
static void tS4(Table& T, Dim o, Z& tmp) {  // Z^2 + Z - 2P^2
  assert(T.d.J >= o.J + 2 && T.d.L >= o.L + 2);
  auto& t = T.t;
  for (int j = 0; j <= o.J; ++j)
    for (int b = 0; b <= o.L - j; ++b) {
      mpz_add(tmp.get_mpz_t(), t[j][b + 2].get_mpz_t(), t[j][b + 1].get_mpz_t());
      mpz_submul_ui(tmp.get_mpz_t(), t[j + 2][b].get_mpz_t(), 2);
      mpz_swap(tmp.get_mpz_t(), t[j][b].get_mpz_t());
    }
  T.d = o;
}
static void tZm1(Table& T, Dim o) {  // Z - 1
  assert(T.d.J >= o.J && T.d.L >= o.L + 1);
  auto& t = T.t;
  for (int j = 0; j <= o.J; ++j)
    for (int b = 0; b <= o.L - j; ++b)
      mpz_sub(t[j][b].get_mpz_t(), t[j][b + 1].get_mpz_t(), t[j][b].get_mpz_t());
  T.d = o;
}
// h4 = Z^2 + ZP - P^2 - 2P - 1 (sg = +1); W5 = Z^2 - ZP - P^2 + 2P - 1 (sg = -1)
static void tH(Table& T, Dim o, int sg, Z& tmp) {
  assert(T.d.J >= o.J + 2 && T.d.L >= o.L + 2);
  auto& t = T.t;
  for (int j = 0; j <= o.J; ++j)
    for (int b = 0; b <= o.L - j; ++b) {
      mpz_sub(tmp.get_mpz_t(), t[j][b + 2].get_mpz_t(), t[j + 2][b].get_mpz_t());
      mpz_sub(tmp.get_mpz_t(), tmp.get_mpz_t(), t[j][b].get_mpz_t());
      if (sg > 0) {
        mpz_add(tmp.get_mpz_t(), tmp.get_mpz_t(), t[j + 1][b + 1].get_mpz_t());
        mpz_submul_ui(tmp.get_mpz_t(), t[j + 1][b].get_mpz_t(), 2);
      } else {
        mpz_sub(tmp.get_mpz_t(), tmp.get_mpz_t(), t[j + 1][b + 1].get_mpz_t());
        mpz_addmul_ui(tmp.get_mpz_t(), t[j + 1][b].get_mpz_t(), 2);
      }
      mpz_swap(tmp.get_mpz_t(), t[j][b].get_mpz_t());
    }
  T.d = o;
}

static size_t ix(int T, int m, int r) { return size_t(m) * (T + 1) - size_t(m) * (m - 1) / 2 + r; }

struct Coef { int j; unsigned long mag; bool neg; };
static std::vector<std::vector<std::vector<Coef>>> coeff;  // [l][alpha]

using QKey = std::array<int, 9>;  // A E b alpha gamma m n q w
static std::map<QKey, Z> queries;

struct Slot {
  std::atomic<int> item{-1}, v{0}, pmax{0};
  std::atomic<unsigned long long> entries{0}, expected{0};
};

int main(int argc, char** argv) {
  if (argc > 1 && (!std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help"))) {
    std::puts("FM-SEC139 exact labels <= 5 box evaluator; reads packet on stdin");
    return 0;
  }
  int nrec, threads, nitems, nq, nmq, heartbeat;
  std::cin >> nrec >> threads >> nitems >> nq >> nmq >> heartbeat;
  omp_set_num_threads(threads);
  // records grouped by (N, n)
  std::map<std::pair<int, int>, std::vector<Rec>> byNn;
  int HMAX = 0;
  for (int z = 0; z < nrec; ++z) {
    int N, l, m, n, p, B; std::cin >> N >> l >> m >> n >> p >> B;
    byNn[{N, n}].push_back({l, m, p, B});
    HMAX = std::max(HMAX, l);
  }
  std::vector<Item> items(nitems);
  for (auto& it : items) std::cin >> it.id >> it.N >> it.a >> it.n >> it.v >> it.expected >> it.nq;
  for (int z = 0; z < nq; ++z) {
    QKey k; for (int& i : k) std::cin >> i;
    std::string v; std::cin >> v; queries.emplace(k, Z(v));
  }
  assert((int)queries.size() == nq);
  // moment queries: m r b value (independent large-b checks of the column)
  std::map<std::array<int, 3>, Z> mq;
  for (int z = 0; z < nmq; ++z) {
    int m, r, b; std::string v; std::cin >> m >> r >> b >> v; mq.emplace(std::array<int, 3>{m, r, b}, Z(v));
  }
  // coefficients of (1+t)^alpha (1-t)^(l-alpha)
  coeff.resize(HMAX + 1);
  for (int l = 0; l <= HMAX; ++l) {
    coeff[l].resize(l + 1);
    for (int al = 0; al <= l; ++al) {
      std::vector<Z> row(l + 1); row[0] = 1;
      if (l) row[1] = 2 * al - l;
      for (int j = 1; j < l; ++j) {
        Z v = (2 * al - l) * row[j] - (l - j + 1) * row[j - 1];
        assert(mpz_divisible_ui_p(v.get_mpz_t(), j + 1));
        mpz_divexact_ui(v.get_mpz_t(), v.get_mpz_t(), j + 1);
        row[j + 1] = v;
      }
      for (int j = 0; j <= l; ++j)
        if (row[j] != 0) {
          Z a = abs(row[j]);
          assert(mpz_fits_ulong_p(a.get_mpz_t()));
          coeff[l][al].push_back({j, a.get_ui(), row[j] < 0});
        }
    }
  }
  // per item: relevant records and dims
  struct Plan {
    std::vector<std::vector<std::vector<std::pair<int, Rec>>>> leaves;  // [p][m] -> records
    std::vector<Dim> own;      // [p]
    std::vector<std::vector<Dim>> need;  // [u][v]
    std::vector<int> mmax;     // [p]
    int pmax = -1;
    Dim root;
  };
  auto make_plan = [&](const Item& it) {
    // work item (N, a, n, v): v copies of W5 (plus label 5); nodes u = 0..pmax-v
    // carry q = u copies of h4, so p = u + v
    Plan P;
    int A = 2 * it.a, E = 2 * (it.N - it.a), v = it.v;
    auto f = byNn.find({it.N, it.n});
    assert(f != byNn.end());
    for (auto& r : f->second) P.pmax = std::max(P.pmax, r.p);
    P.leaves.assign(P.pmax + 1, {});
    P.own.assign(P.pmax + 1, Dim());
    P.mmax.assign(P.pmax + 1, -1);
    for (auto& r : f->second) {
      if (r.p < v) continue;
      int q = r.p - v;
      int lo = std::max(0, r.l + r.p + r.m - A - q), hi = std::min(r.l, E - r.m - q);
      if (hi < lo) continue;
      auto& lv = P.leaves[r.p];
      if ((int)lv.size() < r.m + 1) lv.resize(r.m + 1);
      lv[r.m].push_back({0, r});
      P.own[r.p] = join(P.own[r.p], Dim{r.l, r.B - 1 + r.l + r.m});
      P.mmax[r.p] = std::max(P.mmax[r.p], r.m);
    }
    int pm = P.pmax;
    P.need.assign(pm + 2, std::vector<Dim>(pm + 2));  // need[u][v] (only column v used)
    for (int u = pm - v; u >= 0; --u) {
      Dim x = P.own[u + v];
      if (u + 1 + v <= pm) x = join(x, plus2(P.need[u + 1][v]));
      P.need[u][v] = x;
    }
    Dim r = P.need[0][v];
    for (int k = 0; k < it.n + v; ++k) r = plus2(r);
    P.root = r;
    return P;
  };
  // moment columns needed: (a, r) -> max b
  std::map<std::pair<int, int>, int> colmax;
  std::vector<Dim> roots(nitems);
  int TMAX = 0;
  for (int i = 0; i < nitems; ++i) {
    Plan P = make_plan(items[i]);
    roots[i] = P.root;
    if (P.root.J < 0) continue;
    auto key = std::make_pair(items[i].a, items[i].N - items[i].a);
    colmax[key] = std::max(colmax[key], P.root.L);
    TMAX = std::max(TMAX, items[i].N + P.root.L);
  }
  // bridge points: for each column, t[j][b] for j in {1,2} at b in {0, L-j}
  // needs M(a+h, r+j-h, b); collect requested moment points
  std::map<std::array<int, 3>, Z> extra;  // (m, r, b) -> value
  for (auto& [k, L] : colmax)
    for (int j = 1; j <= 2; ++j)
      for (int b : {0, std::max(0, L - j)})
        for (int h = 0; h <= j; ++h) extra[{k.first + h, k.second + j - h, b}] = 0;
  for (auto& [k, v] : mq) extra[k] = 0;
  for (auto& [k, v] : extra) TMAX = std::max(TMAX, k[0] + k[1] + k[2]);
  // the layer sweep runs to TS; points beyond TS (only the long N = 2
  // columns) come from the direct formula M = sum_k c_k sum_i C(b,i)
  // mu(k,i) mu(2N-k,b-i), mu(k,i) = E[x^k (x^2-1)^i], with an overlap check
  int TS = 0;
  for (auto& [k, L] : colmax)
    if (k.first + k.second >= 3) TS = std::max(TS, k.first + k.second + L);
  for (auto& [k, v] : mq) if (k[0] + k[1] + k[2] <= 600) TS = std::max(TS, k[0] + k[1] + k[2]);
  TS = std::min(TS, TMAX);
  std::atomic<int> phase_T{-1};
  say("SETUP items %d columns %zu TMAX %d sweep %d queries %d moment_queries %d\n", nitems, colmax.size(), TMAX, TS, nq, nmq);
  unsigned long long total_expected = 0;
  for (auto& it : items) total_expected += it.expected;
  std::vector<Slot> slots(threads);
  std::atomic<int> done{0};
  std::atomic<unsigned long long> entries_done{0};
  std::atomic<bool> finished{false};
  double started = omp_get_wtime();
  std::thread hb([&]() {
    auto last = std::chrono::steady_clock::now();
    while (!finished.load()) {
      std::this_thread::sleep_for(std::chrono::milliseconds(500));
      auto now = std::chrono::steady_clock::now();
      if (std::chrono::duration_cast<std::chrono::seconds>(now - last).count() < heartbeat) continue;
      last = now;
      std::string act;
      int busy = 0;
      for (int k = 0; k < threads; ++k) {
        int id = slots[k].item.load();
        if (id < 0) continue;
        ++busy;
        if (act.size() < 1500) {
          char buf[160];
          std::snprintf(buf, sizeof buf, " [item %d v %d/%d entries %llu/%llu]", id, slots[k].v.load(),
                        slots[k].pmax.load(), slots[k].entries.load(), slots[k].expected.load());
          act += buf;
        }
      }
      if (phase_T.load() < TS) {
        say("HB moments degree %d/%d seconds %.0f\n", phase_T.load(), TS, omp_get_wtime() - started);
        continue;
      }
      say("HB items %d/%d entries %llu/%llu seconds %.0f busy %d%s\n", done.load(), nitems,
          entries_done.load(), total_expected, omp_get_wtime() - started, busy, act.c_str());
    }
  });
  // moment sweep over layers T, keeping two layers
  std::map<std::pair<int, int>, std::vector<Z>> cols;
  for (auto& [k, L] : colmax) cols[k].resize(L + 1);
  {
    std::vector<Z> cat(TS + 1); cat[0] = 1;
    for (int j = 1; j <= TS; ++j) { cat[j] = cat[j - 1] * (4 * j - 2); exact_ui(cat[j], j + 1); }
    std::vector<Z> prev, cur;
    double t0 = omp_get_wtime();
    for (int T = 0; T <= TS; ++T) {
      cur.assign(size_t(T + 1) * (T + 2) / 2, Z(0));
      std::vector<Z> cc(T + 1);
      #pragma omp parallel for schedule(static)
      for (int k = 0; k <= T; ++k) cc[k] = cat[k] * cat[T - k];
      #pragma omp parallel for schedule(dynamic, 4)
      for (int m = 0; m <= T / 2; ++m) {
        int r = T - m, deg = 2 * T, d = 2 * (m - r);
        Z c0 = 1, c1 = d, ans = cc[0], tmp;
        for (int j = 1; j < deg; ++j) {
          tmp = d * c1 - (deg - j + 1) * c0; exact_ui(tmp, j + 1);
          c0 = c1; c1 = tmp;
          if ((j + 1) % 2 == 0) ans += c1 * cc[(j + 1) / 2];
        }
        if (T == 0) ans = 1;
        cur[ix(T, r, m)] = ans;
        cur[ix(T, m, r)] = ans;
      }
      for (int b = 1; b <= T; ++b) {
        for (int m = 0; m <= T - b; ++m) {
          int r = T - b - m;
          Z v = cur[ix(T, m + 1, r)] + cur[ix(T, m, r + 1)];
          exact_ui(v, 2);
          v -= 2 * prev[ix(T - 1, m, r)];
          cur[ix(T, m, r)] = std::move(v);
        }
      }
      for (auto& [k, col] : cols) {
        int b = T - k.first - k.second;
        if (b >= 0 && b < (int)col.size()) col[b] = cur[ix(T, k.first, k.second)];
      }
      for (auto& [k, v] : extra) if (k[0] + k[1] + k[2] == T) v = cur[ix(T, k[0], k[1])];
      std::swap(prev, cur);
      phase_T = T;
      if (T % 100 == 0 || T == TS) say("MOMENTS degree %d/%d seconds %.1f\n", T, TS, omp_get_wtime() - t0);
    }
  }
  {
    int Nmu = 0, Dmu = -1;
    for (auto& [k, col] : cols)
      if (k.first + k.second + (int)col.size() - 1 > TS) {
        Nmu = std::max(Nmu, k.first + k.second); Dmu = std::max(Dmu, (int)col.size() - 1);
      }
    for (auto& [k, v] : extra)
      if (k[0] + k[1] + k[2] > TS) { Nmu = std::max(Nmu, k[0] + k[1]); Dmu = std::max(Dmu, k[2]); }
    if (Dmu >= 0) {
      double t0 = omp_get_wtime();
      int K0 = 2 * Nmu + 2 * Dmu;
      std::vector<Z> catm(K0 / 2 + 1); catm[0] = 1;
      for (int j = 1; j <= K0 / 2; ++j) { catm[j] = catm[j - 1] * (4 * j - 2); exact_ui(catm[j], j + 1); }
      std::vector<Z> row(K0 + 1);
      for (int k = 0; k <= K0; ++k) row[k] = (k % 2 == 0) ? catm[k / 2] : Z(0);
      std::vector<std::vector<Z>> mu(Dmu + 1);
      mu[0].assign(row.begin(), row.begin() + 2 * Nmu + 1);
      for (int i = 1; i <= Dmu; ++i) {
        int len = K0 - 2 * i;
        for (int k = 0; k <= len; ++k) row[k] = row[k + 2] - row[k];
        mu[i].assign(row.begin(), row.begin() + 2 * Nmu + 1);
      }
      auto direct_moment = [&](int m, int r, int b) {
        int NN = m + r, deg = 2 * NN, d = 2 * (m - r);
        std::vector<Z> c(deg + 1);
        c[0] = 1;
        if (deg) c[1] = d;
        for (int j = 1; j < deg; ++j) {
          Z t = d * c[j] - (deg - j + 1) * c[j - 1]; exact_ui(t, j + 1); c[j + 1] = t;
        }
        Z total = 0, inner, bin, t2;
        for (int k = 0; k <= deg; k += 2) {
          if (c[k] == 0) continue;
          inner = 0; bin = 1;
          for (int i = 0; i <= b; ++i) {
            t2 = mu[i][k] * mu[b - i][deg - k];
            inner += bin * t2;
            bin *= (b - i); exact_ui(bin, i + 1);
          }
          total += c[k] * inner;
        }
        return total;
      };
      unsigned long long overlap = 0;
      for (auto& [k, col] : cols) {
        int NN = k.first + k.second;
        if (NN + (int)col.size() - 1 <= TS) continue;
        int bs = TS - NN;  // last swept b
        #pragma omp parallel for schedule(dynamic, 4)
        for (int b = bs + 1; b < (int)col.size(); ++b) col[b] = direct_moment(k.first, k.second, b);
        for (int b = std::max(0, bs - 20); b <= bs; ++b) {
          assert(direct_moment(k.first, k.second, b) == col[b]);
          ++overlap;
        }
      }
      for (auto& [k, v] : extra)
        if (k[0] + k[1] + k[2] > TS) v = direct_moment(k[0], k[1], k[2]);
      say("DIRECT MOMENTS beyond degree %d: columns up to b = %d, overlap checks %llu, seconds %.1f\n",
          TS, Dmu, overlap, omp_get_wtime() - t0);
    }
  }
  phase_T = TS;
  for (auto& [k, v] : mq) { assert(extra[k] == v); }
  say("MOMENT CHECKS %d passed\n", nmq);
  // order items by estimated cost (root size times expected), largest first
  std::vector<int> order(nitems);
  for (int i = 0; i < nitems; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](int x, int y) {
    double cx = double(items[x].expected) + 50.0 * (roots[x].J + 1) * (roots[x].L + 1);
    double cy = double(items[y].expected) + 50.0 * (roots[y].J + 1) * (roots[y].L + 1);
    return cx > cy;
  });
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    Table base, W, Hh, M;
    Z tmp, value, least;
    #pragma omp for schedule(dynamic, 1)
    for (int oi = 0; oi < nitems; ++oi) {
      const Item& it = items[order[oi]];
      double t1 = omp_get_wtime();
      int N = it.N, a = it.a, r = N - a, A = 2 * a, E = 2 * r;
      Plan P = make_plan(it);
      unsigned long long local = 0, zl = 0, qseen = 0, bridges = 0;
      least = 0;
      Slot& S = slots[tid];
      S.expected = it.expected; S.entries = 0; S.v = 0; S.pmax = P.pmax; S.item = it.id;
      if (P.root.J >= 0) {
        // mixed table t[j][b] = E[s^A d^E P^j Z^b]
        const auto& col = cols.at({a, r});
        Dim R = P.root;
        base.ensure(R);
        for (int b = 0; b <= R.L; ++b) base.t[0][b] = col[b];
        for (int j = 1; j <= R.J; ++j)
          for (int b = 0; b <= R.L - j; ++b) {
            Z v = 8 * (a - r) * base.t[j - 1][b];
            if (j >= 2) v += 4 * (j - 1) * base.t[j - 2][b + 1] + 8 * (j - 1) * base.t[j - 2][b];
            if (b) v += 12 * b * base.t[j][b - 1];
            exact_ui(v, 2 * (a + r + b + j + 2));
            base.t[j][b] = std::move(v);
          }
        base.d = R;
        // bridges against P = (s^2 - d^2)/4
        int L = colmax.at({a, r});
        for (int j = 1; j <= std::min(2, R.J); ++j)
          for (int b : {0, std::max(0, L - j)}) {
            if (b > R.L - j) continue;
            Z ans = 0, c;
            for (int h = 0; h <= j; ++h) {
              mpz_bin_uiui(c.get_mpz_t(), j, h);
              Z term = c * extra.at({a + h, r + j - h, b});
              if ((j - h) % 2) ans -= term; else ans += term;
            }
            assert(mpz_divisible_2exp_p(ans.get_mpz_t(), 2 * j));
            mpz_tdiv_q_2exp(ans.get_mpz_t(), ans.get_mpz_t(), 2 * j);
            assert(ans == base.t[j][b]);
            ++bridges;
          }
        // hat S_4 ^ n
        for (int k = 0; k < it.n; ++k) {
          Dim o{base.d.J - 2, base.d.L - 2};
          tS4(base, o, tmp);
        }
        // v copies of W5 = hat S_5 / s (plus label 5), in place on base
        {
          int v = it.v;
          for (int k = 0; k < v; ++k) {
            Dim o{base.d.J - 2, base.d.L - 2};
            tH(base, o, -1, tmp);
          }
          S.v = v;
          assert(base.d.J >= P.need[0][v].J && base.d.L >= P.need[0][v].L);
          Hh.copy_from(base, P.need[0][v]);
          for (int u = 0; u + v <= P.pmax; ++u) {
            int p = u + v, q = u, w = v;
            if (P.own[p].J >= 0) {
              M.copy_from(Hh, P.own[p]);
              for (int m = 0; m <= P.mmax[p]; ++m) {
                if (m < (int)P.leaves[p].size())
                  for (auto& [dummy, rec] : P.leaves[p][m]) {
                    int l = rec.l;
                    int lo = std::max(0, l + p + m - A - q), hi = std::min(l, E - m - q);
                    for (int al = lo; al <= hi; ++al) {
                      const auto& cf = coeff[l][al];
                      for (int b = 0; b < rec.B; ++b) {
                        value = 0;
                        for (const auto& c : cf) {
                          const Z& x = M.t[c.j][b + l - c.j];
                          if (c.neg) mpz_submul_ui(value.get_mpz_t(), x.get_mpz_t(), c.mag);
                          else mpz_addmul_ui(value.get_mpz_t(), x.get_mpz_t(), c.mag);
                        }
                        int sgn = mpz_sgn(value.get_mpz_t());
                        if (sgn < 0) {
                          say("NEGATIVE item %d A %d E %d b %d alpha %d gamma %d m %d n %d q %d w %d value %s\n",
                              it.id, A, E, b, al, l - al, m, it.n, q, w, value.get_str().c_str());
                          std::exit(2);
                        }
                        ++local;
                        if (sgn == 0) ++zl;
                        else if (least == 0 || value < least) least = value;
                        if (N <= 6 && b <= 2 && l + m + it.n + p <= 6) {
                          QKey key{A, E, b, al, l - al, m, it.n, q, w};
                          auto f = queries.find(key);
                          if (f != queries.end()) { assert(value == f->second); ++qseen; }
                        }
                      }
                    }
                    S.entries = local;
                  }
                if (m < P.mmax[p]) tZm1(M, Dim{P.own[p].J, P.own[p].L - m - 1});
              }
            }
            if (p + 1 > P.pmax || P.need[u + 1][v].J < 0) break;  // no later node is needed
            tH(Hh, P.need[u + 1][v], +1, tmp);
          }
        }
      }
      if (local != it.expected) {
        say("COUNT MISMATCH item %d got %llu expected %llu\n", it.id, local, it.expected);
        std::exit(3);
      }
      if ((int)qseen != it.nq) {
        say("QUERY MISMATCH item %d got %llu expected %d\n", it.id, qseen, it.nq);
        std::exit(4);
      }
      entries_done += local;
      S.item = -1;
      int dn = ++done;
      say("ITEM %d N %d a %d n %d v %d entries %llu zeros %llu least %s queries %llu bridges %llu seconds %.3f done %d/%d\n",
          it.id, N, a, it.n, it.v, local, zl, least.get_str().c_str(), qseen, bridges, omp_get_wtime() - t1, dn, nitems);
    }
  }
  finished = true;
  hb.join();
  say("ALL ITEMS %d entries %llu seconds %.1f\n", nitems, entries_done.load(), omp_get_wtime() - started);
  return 0;
}
'''

def locate(name):
    return subprocess.check_output(
        ["g++", "-print-file-name="+name], text=True).strip()

def build():
    obj = os.memfd_create("sec139-object", 0)
    exe = os.memfd_create("sec139-executable", 0)
    subprocess.run(
        ["g++", "-O3", "-std=c++17", "-fno-pie", "-fopenmp",
         "-pipe", "-x", "c++", "-", "-c", "-o", f"/proc/self/fd/{obj}"],
        input=CPP, text=True, pass_fds=(obj,), check=True)
    subprocess.run(
        ["ld", "-L"+os.path.dirname(locate("crtbegin.o")),
         "-o", f"/proc/self/fd/{exe}",
         "-dynamic-linker", "/lib64/ld-linux-x86-64.so.2",
         locate("crt1.o"), locate("crti.o"), locate("crtbegin.o"),
         f"/proc/self/fd/{obj}", locate("libstdc++.so"),
         "-lgmpxx", "-lgmp", locate("libgomp.so"), "-lm",
         locate("libgcc_s.so"), "-lc",
         locate("crtend.o"), locate("crtn.o")],
        pass_fds=(obj, exe), check=True)
    os.close(obj)
    return exe

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--state", default=None,
                    help="checkpoint directory (resumable); none: no checkpoint")
    ap.add_argument("--stride", type=int, default=1,
                    help="calibration: run every K-th item only")
    ap.add_argument("--max-n", type=int, default=None,
                    help="test: only items with N <= this")
    ap.add_argument("--queries", type=int, default=300)
    ap.add_argument("--heartbeat", type=int, default=60)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--recheck", default=None,
                    help="test: re-evaluate items finished in this other state dir")
    ap.add_argument("--recheck-limit", type=int, default=64)
    ap.add_argument("--no-mirror", action="store_true",
                    help="evaluate the mirrored backgrounds 2a > N as well")
    args = ap.parse_args()

    print(stamp(), "building box records", flush=True)
    tasks = [(H, l) for H in range(67, 0, -1) for l in range(H+1)]
    profiles = 0
    recs = []
    with Pool(args.workers) as pool:
        for pr, rs in pool.imap_unordered(profile_task, tasks):
            profiles += pr
            recs.extend(rs)
        recs.sort()
        expanded = len(recs)
        entries = sum(r[5]*r[6] for r in recs)
        Nmax = max(r[0] for r in recs)
        Bmax = max(r[5] for r in recs)
        assert (profiles, expanded, entries, Nmax, Bmax) == (
            916895, 1738408, 227336512420, 71, 1006)
        print(stamp(), "box: profiles", profiles, "records", expanded,
              "entries", entries, "(matches FM-SEC135)", flush=True)
        recs6 = [r[:6] for r in recs]
        chunks = [recs6[i::400] for i in range(400)]
        itemcount = {}
        for out in pool.imap_unordered(item_task, chunks):
            for k, v in out.items():
                itemcount[k] = itemcount.get(k, 0)+v
    assert sum(itemcount.values()) == entries
    keys = sorted(itemcount)
    ids = {k: i for i, k in enumerate(keys)}
    print(stamp(), "work items", len(keys), flush=True)

    # direct Catalan checks (small words)
    pool_q = []
    for (N, l, m, n, p, B) in recs6:
        if N > 6 or l+m+n+p > 6:
            continue
        for a in range(N):
            A, E = 2*a, 2*(N-a)
            for q in range(p+1):
                lo = max(0, l+p+m-A-q)
                hi = min(l, E-m-q)
                for al in range(lo, hi+1):
                    for b in range(min(B, 3)):
                        pool_q.append((A, E, b, al, l-al, m, n, q, p-q))
    qkeys = random.Random(139).sample(pool_q, args.queries)
    qvals = [(k, direct(k)) for k in qkeys]
    print(stamp(), "direct Catalan values", len(qvals), "of", len(pool_q),
          "small entries; negative among them:",
          sum(v < 0 for _, v in qvals), flush=True)
    perq = {}
    for k, _ in qvals:
        A, E, n, w = k[0], k[1], k[6], k[8]
        it = ((A+E)//2, A//2, n, w)
        perq[it] = perq.get(it, 0)+1

    # exact count symmetry of the box under (N, a, n) <-> (N, N-a, n)
    grp = {}
    for (N, a, n, v), c in itemcount.items():
        grp[(N, a, n)] = grp.get((N, a, n), 0)+c
    for (N, a, n), c in grp.items():
        if a >= 1:
            assert grp.get((N, N-a, n)) == c, (N, a, n)
    mirror = not args.no_mirror
    def weight(k):
        N, a = k[0], k[1]
        return 2 if (mirror and a >= 1 and 2*a < N) else 1
    assert sum(weight(k)*c for k, c in itemcount.items()
               if not (mirror and 2*k[1] > k[0] and k[1] >= 1)) == entries
    print(stamp(), "mirror symmetry of the entry counts: PASS;",
          "mirror mode" if mirror else "no-mirror mode", flush=True)
    selected = [k for k in keys
                if ids[k] % args.stride == 0
                and (args.max_n is None or k[0] <= args.max_n)
                and not (mirror and k[1] >= 1 and 2*k[1] > k[0])]
    reference = {}
    if args.recheck:
        for line in open(os.path.join(args.recheck, "items_done.txt")):
            f = line.split()
            reference[int(f[1])] = line.rstrip("\n")
        rid = list(reference)[:args.recheck_limit]
        selected = [keys[i] for i in rid]
        print(stamp(), "recheck: re-evaluating", len(selected), "items from", args.recheck, flush=True)
    done = {}
    if args.state:
        os.makedirs(args.state, exist_ok=True)
        path = os.path.join(args.state, "items_done.txt")
        if os.path.exists(path):
            for line in open(path):
                f = line.split()
                if len(f) > 2 and f[0] == "ITEM":
                    done[int(f[1])] = line.rstrip("\n")
        print(stamp(), "checkpoint:", len(done), "items already done", flush=True)
    todo = [k for k in selected if ids[k] not in done]
    # moment checks at large b (only columns used by this run)
    cand = sorted({(a, N-a) for (N, a, n, v) in todo})
    mpts = []
    if (0, 2) in cand:
        mpts += [(0, 2, 1000), (1, 1, 1009)]
    rnd = random.Random(1391)
    for (a, r) in rnd.sample(cand, min(6, len(cand))):
        mpts.append((a, r, 60+rnd.randrange(150)))
    mpts = sorted(set(mpts))
    mvals = [(pt, moment_check(*pt)) for pt in mpts]
    print(stamp(), "independent large-b moment values", len(mvals), flush=True)

    print(stamp(), "selected items", len(selected), "to run", len(todo),
          "entries to run", sum(itemcount[k] for k in todo), flush=True)
    if todo:
        exe = build()
        packet = [f"{len(recs6)} {args.threads} {len(todo)} {len(qvals)} "
                  f"{len(mvals)} {args.heartbeat}"]
        packet += [" ".join(map(str, r)) for r in recs6]
        packet += [f"{ids[k]} {k[0]} {k[1]} {k[2]} {k[3]} {itemcount[k]} {perq.get(k, 0)}"
                   for k in todo]
        packet += [" ".join(map(str, k))+" "+str(v) for k, v in qvals]
        packet += [f"{a} {r} {b} {v}" for (a, r, b), v in mvals]
        proc = subprocess.Popen([f"/proc/self/fd/{exe}"], stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, text=True, pass_fds=(exe,),
                                bufsize=1)
        sink = open(os.path.join(args.state, "items_done.txt"), "a") if args.state else None
        import threading
        def feed():
            proc.stdin.write("\n".join(packet)+"\n")
            proc.stdin.close()
        threading.Thread(target=feed, daemon=True).start()
        for line in proc.stdout:
            line = line.rstrip("\n")
            if line.startswith("ITEM"):
                f = line.split()
                done[int(f[1])] = line
                if sink:
                    sink.write(line+"\n")
                    sink.flush()
                    os.fsync(sink.fileno())
                if args.heartbeat <= 5:
                    print(stamp(), line, flush=True)
            else:
                print(stamp(), line, flush=True)
        rc = proc.wait()
        if sink:
            sink.close()
        if rc != 0:
            print(stamp(), "evaluator FAILED with code", rc, flush=True)
            sys.exit(rc)
    # aggregate over the selected items
    tot = zeros = qs = br = 0
    least = None
    missing = 0
    for k in selected:
        line = done.get(ids[k])
        if line is None:
            missing += 1
            continue
        f = line.split()
        d = dict(zip(f[2::2], f[3::2]))
        assert (int(d["N"]), int(d["a"]), int(d["n"]), int(d["v"])) == k
        assert int(d["entries"]) == itemcount[k]
        assert int(d["queries"]) == perq.get(k, 0)
        tot += weight(k)*int(d["entries"])
        zeros += int(d["zeros"])
        qs += int(d["queries"])
        br += int(d["bridges"])
        v = int(d["least"])
        if v > 0 and (least is None or v < least):
            least = v
    # extra symmetry check: mirrored items finished in an earlier run
    if mirror:
        fin = {}
        for i, line in done.items():
            f = line.split()
            dd = dict(zip(f[2::2], f[3::2]))
            g = (int(dd["N"]), int(dd["a"]), int(dd["n"]))
            fin.setdefault(g, []).append(dd)
        def summ(lst):
            least_ = [int(x["least"]) for x in lst if int(x["least"]) > 0]
            return (sum(int(x["entries"]) for x in lst),
                    sum(int(x["zeros"]) for x in lst), min(least_) if least_ else 0)
        nmir = 0
        nitems_grp = {}
        for k in keys:
            nitems_grp[k[:3]] = nitems_grp.get(k[:3], 0)+1
        for (N, a, n), lst in fin.items():
            if a >= 1 and 2*a > N and (N, N-a, n) in fin:
                full_a = nitems_grp[(N, a, n)]
                full_b = nitems_grp[(N, N-a, n)]
                if len(lst) == full_a and len(fin[(N, N-a, n)]) == full_b:
                    assert summ(lst) == summ(fin[(N, N-a, n)]), (N, a, n)
                    nmir += 1
        print(stamp(), "mirrored groups also evaluated directly and equal:", nmir, flush=True)
    if args.recheck:
        same = 0
        t_old = t_new = 0.0
        for k in selected:
            fo = reference[ids[k]].split(); fn = done[ids[k]].split()
            do = dict(zip(fo[2::2], fo[3::2])); dn = dict(zip(fn[2::2], fn[3::2]))
            keyf = ("N", "a", "n", "v", "entries", "zeros", "least", "queries")
            assert all(do[x] == dn[x] for x in keyf), (k, do, dn)
            same += 1
            t_old += float(do["seconds"]); t_new += float(dn["seconds"])
        print(stamp(), "recheck: identical result fields on", same, "items;",
              f"item seconds old {t_old:.1f}, new {t_new:.1f}", flush=True)
    print(stamp(), "aggregate: items", len(selected)-missing, "/", len(selected),
          "entries", tot, "zeros", zeros, "least positive", least,
          "direct checks", qs, "bridges", br, flush=True)
    if missing == 0 and args.stride == 1 and args.max_n is None and not args.recheck:
        assert tot == 227336512420
        assert qs == sum(perq.get(k, 0) for k in selected)
        print("EXACT LABELS<=5 BOX: all", tot, "entries >= 0")
        print("PASS")

if __name__ == "__main__":
    main()
