import os
import subprocess
import sys

if any(x in ("-h", "--help") for x in sys.argv[1:]):
    print("FM-SEC152: exact Rule 1' sweep through W=40; memory-backed C++17/OpenMP.")
    print("Run with python3 -u -; set OMP_NUM_THREADS to change the thread count.")
    raise SystemExit(0)

src = r"""
#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>
using namespace std;
using I = __int128_t;
using U = __uint128_t;
using boost::multiprecision::int256_t;
static const int MAXW = 40, L = 40;
static const int RULE = 1;  // 1 = Rule M, 0 = Rule W

struct HashU {
    size_t operator()(U x) const noexcept {
        uint64_t a = (uint64_t)x, b = (uint64_t)(x >> 64);
        a ^= b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2);
        return (size_t)a;
    }
};
struct Rem { int n = 0, m = 0; }; // m=0 means remove one n
struct Rec {
    array<signed char, L + 1> c{};
    U key = 0;
    int W = 0, nf = 0, mx = 0;
    array<I, MAXW / 2 + 1> g{};
};
struct MinR {
    bool yes = false;
    I num = 0;
    U den = 1, key = 0;
    int W = 0, p = 0, branch = 0;
    Rem R{};
    I parent = 0, child = 0;
};
struct BranchStats {
    long long backgrounds = 0, failures = 0, child_zero = 0, child_negative = 0;
    MinR min_positive_child, min_all_nonzero;
};
struct Failure {
    bool yes = false;
    int W = 0, p = 0;
    U key = 0;
    Rem R{};
    I parent = 0, child = 0;
};
struct Local {
    long long eligible = 0, tests = 0, failures = 0, unproved_failures = 0, unproved_tests = 0;
    Failure ufail;
    BranchStats branch[3];
    Failure failure;
};

array<U, L + 1> mult{};
array<signed char, L + 1> cur{};
vector<Rec> recs;

U absI(I x) { return x < 0 ? U(-(x + 1)) + 1 : U(x); }

string strU(U x) {
    if (!x) return "0";
    string s;
    while (x) {
        s.push_back(char('0' + x % 10));
        x /= 10;
    }
    reverse(s.begin(), s.end());
    return s;
}
string strI(I x) {
    return x < 0 ? "-" + strU(absI(x)) : strU(U(x));
}
U gcdU(U a, U b) {
    while (b) {
        U t = a % b;
        a = b;
        b = t;
    }
    return a;
}
bool ratioLess(I an, U ad, I bn, U bd) {
    return int256_t(an) * int256_t(bd) < int256_t(bn) * int256_t(ad);
}
bool failureLess(const Failure& a, const Failure& b) {
    if (!b.yes) return true;
    if (a.W != b.W) return a.W < b.W;
    if (a.key != b.key) return a.key < b.key;
    if (a.p != b.p) return a.p < b.p;
    if (a.R.n != b.R.n) return a.R.n < b.R.n;
    return a.R.m < b.R.m;
}
void updateMin(MinR& best, I num, U den, const Rec& r, int p, Rem R,
               int branch, I parent, I child, bool normalize_sign) {
    if (!den) return;
    if (normalize_sign && child < 0) num = -num;
    MinR candidate;
    candidate.yes = true;
    candidate.num = num;
    candidate.den = den;
    candidate.key = r.key;
    candidate.W = r.W;
    candidate.p = p;
    candidate.R = R;
    candidate.branch = branch;
    candidate.parent = parent;
    candidate.child = child;

    if (!best.yes ||
        ratioLess(candidate.num, candidate.den, best.num, best.den) ||
        (!ratioLess(best.num, best.den, candidate.num, candidate.den) &&
         !ratioLess(candidate.num, candidate.den, best.num, best.den) &&
         (candidate.W < best.W ||
          (candidate.W == best.W &&
           (candidate.key < best.key ||
            (candidate.key == best.key && candidate.p < best.p))))))
        best = candidate;
}
void mergeMin(MinR& x, const MinR& y) {
    if (!y.yes) return;
    if (!x.yes ||
        ratioLess(y.num, y.den, x.num, x.den) ||
        (!ratioLess(x.num, x.den, y.num, y.den) &&
         !ratioLess(y.num, y.den, x.num, x.den) &&
         (y.W < x.W ||
          (y.W == x.W && (y.key < x.key ||
                          (y.key == x.key && y.p < x.p))))))
        x = y;
}
void generate(int n, int W, int nf, int mx, U key) {
    if (n > L) {
        Rec r;
        r.c = cur;
        r.key = key;
        r.W = W;
        r.nf = nf;
        r.mx = mx;
        recs.push_back(r);
        return;
    }
    int M = MAXW / n;
    for (int c = -M; c <= M; ++c) {
        int nw = W + n * abs(c);
        if (nw > MAXW) continue;
        cur[n] = (signed char)c;
        generate(n + 1, nw, nf + abs(c), c ? n : mx,
                 key + U(c + M) * mult[n]);
    }
    cur[n] = 0;
}
I rangeSum(I pref[2][L + 2], int parity, int lo, int hi) {
    if (lo > hi || hi < 0 || lo > L) return 0;
    lo = max(lo, 0);
    hi = min(hi, L);
    return pref[parity][hi + 1] - pref[parity][lo];
}
array<I, MAXW / 2 + 1> calculate(const Rec& r) {
    int W = r.W;
    I d[L + 1][L + 1]{}, q[L + 1][L + 1]{};
    d[0][0] = 1;

    for (int n = 1; n <= L; ++n) {
        for (int rep = 0; rep < abs((int)r.c[n]); ++rep) {
            int eps = r.c[n] > 0 ? 1 : -1;
            memset(q, 0, sizeof(q));
            I pref[2][L + 2]{};

            // X branch: multiply the X-character by U_n.
            for (int b = 0; b <= W; ++b) {
                pref[0][0] = pref[1][0] = 0;
                for (int x = 0; x <= W; ++x) {
                    pref[0][x + 1] = pref[0][x];
                    pref[1][x + 1] = pref[1][x];
                    pref[x & 1][x + 1] += d[x][b];
                }
                for (int t = 0; t <= W; ++t) {
                    int lo = abs(t - n), hi = min(W, t + n);
                    int par = ((t - n) % 2 + 2) % 2;
                    if (lo <= hi) q[t][b] += rangeSum(pref, par, lo, hi);
                }
            }

            // y branch: choose eps*U_n(y).
            for (int a = 0; a <= W; ++a) {
                pref[0][0] = pref[1][0] = 0;
                for (int y = 0; y <= W; ++y) {
                    pref[0][y + 1] = pref[0][y];
                    pref[1][y + 1] = pref[1][y];
                    pref[y & 1][y + 1] += d[a][y];
                }
                for (int t = 0; t <= W; ++t) {
                    int lo = abs(t - n), hi = min(W, t + n);
                    int par = ((t - n) % 2 + 2) % 2;
                    if (lo <= hi) q[a][t] += eps * rangeSum(pref, par, lo, hi);
                }
            }
            memcpy(d, q, sizeof(d));
        }
    }

    array<I, MAXW / 2 + 1> out{};
    for (int p = W & 1; p <= W; p += 2) out[p / 2] = d[p][0];
    return out;
}
I gp(const Rec& r, int p) {
    if (p < 0 || p > r.W || ((p - r.W) & 1)) return 0;
    return r.g[p / 2];
}
vector<Rem> allowed(const Rec& r) {
    vector<Rem> v;
    for (int n = 2; n <= L; n += 2)
        if (r.c[n]) v.push_back({n, 0});
    for (int n = 1; n <= L; ++n) if (r.c[n])
        for (int m = n; m <= L; ++m)
            if (r.c[m] && n % 2 == m % 2 &&
                (m != n || abs((int)r.c[n]) >= 2))
                v.push_back({n, m});
    return v;
}
U childKey(const Rec& r, Rem x) {
    U k = r.key;
    auto remove_one = [&](int n) {
        if (r.c[n] > 0) k -= mult[n];
        else k += mult[n];
    };
    if (x.m == 0) remove_one(x.n);
    else if (x.n == x.m) {
        if (r.c[x.n] > 0) k -= 2 * mult[x.n];
        else k += 2 * mult[x.n];
    } else {
        remove_one(x.n);
        remove_one(x.m);
    }
    return k;
}
Rem smallestEven(const vector<Rem>& v) {
    for (auto x : v) if (!x.m) return x;
    return v.front();
}
bool pairLarger(Rem x, Rem y) {
    return x.m != y.m ? x.m > y.m : x.n > y.n;
}
struct Pick { Rem R; int branch; }; // 0 even-plus, 1 duplicate, 2 fallback
Pick rule1prime(const Rec& r, const vector<Rem>& v) {
    // First remove the smallest even PLUS label.
    for (int n = 2; n <= L; n += 2)
        if (r.c[n] > 0) return {{n, 0}, 0};

    // Otherwise use Rule 1's duplicate step.
    int bestf = 1, bestn = 99;
    for (int n = 1; n <= L; ++n)
        if (abs((int)r.c[n]) >= 2) {
            int f = abs((int)r.c[n]);
            if (f > bestf || (f == bestf && n < bestn))
                bestf = f, bestn = n;
        }
    if (bestn < 99) return {{bestn, bestn}, 1};

    // All multiplicities are one: largest same-parity pair.
    Rem b{-1, -1};
    for (auto x : v)
        if (x.m && (b.n < 0 || pairLarger(x, b))) b = x;
    if (b.n >= 0) return {b, 2};

    // No same-parity pair: remove the smallest even singleton.
    return {smallestEven(v), 2};
}
Pick ruleW(const Rec& r, const vector<Rem>& v) {
    // Rule W: the label class of largest total weight (multiplicity * label; ties: larger label).
    int bn = -1; long bw = -1;
    for (int n = 1; n <= L; ++n) if (r.c[n]) {
        long w = (long)abs((int)r.c[n]) * n;
        if (w > bw || (w == bw && n > bn)) { bw = w; bn = n; }
    }
    if (abs((int)r.c[bn]) >= 2) return {{bn, bn}, 1};
    for (int m = L; m >= 1; --m)
        if (m != bn && r.c[m] && (m % 2) == (bn % 2)) return {{min(m, bn), max(m, bn)}, 2};
    if (bn % 2 == 0) return {{bn, 0}, 0};
    return {smallestEven(v), 2};
}
Pick ruleM(const Rec& r, const vector<Rem>& v) {
    // Rule M: the maximum label; two copies if repeated, else with the largest other same-parity label,
    // else alone if even, else the smallest even label.
    int mxl = 0; for (int n = 1; n <= L; ++n) if (r.c[n]) mxl = n;
    if (abs((int)r.c[mxl]) >= 2) return {{mxl, mxl}, 1};
    for (int m = mxl - 1; m >= 1; --m) if (r.c[m] && (m % 2) == (mxl % 2)) return {{m, mxl}, 2};
    if (mxl % 2 == 0) return {{mxl, 0}, 0};
    return {smallestEven(v), 2};
}
string show(const Rec& r) {
    string s = "{";
    bool first = true;
    for (int n = 1; n <= L; ++n) if (r.c[n]) {
        if (!first) s += " ";
        first = false;
        s += (r.c[n] > 0 ? "+" : "-") + to_string(n) + "^" +
             to_string(abs((int)r.c[n]));
    }
    return s + "}";
}
bool failureEarlier(const Failure& a, const Failure& b) {
    if (!b.yes) return true;
    if (a.W != b.W) return a.W < b.W;
    if (a.key != b.key) return a.key < b.key;
    if (a.p != b.p) return a.p < b.p;
    if (a.R.n != b.R.n) return a.R.n < b.R.n;
    return a.R.m < b.R.m;
}
int main() {
    U radix = 1;
    for (int n = 1; n <= L; ++n) {
        mult[n] = radix;
        radix *= U(2 * (MAXW / n) + 1);
    }
    if (radix >= (U(1) << 127)) return 2;

    recs.reserve(6100000);
    auto start = chrono::steady_clock::now();
    generate(1, 0, 0, 0, 0);
    sort(recs.begin(), recs.end(), [](const Rec& a, const Rec& b) {
        return a.W != b.W ? a.W < b.W : a.key < b.key;
    });
    cout << "generated and sorted backgrounds=" << recs.size() << "\n" << flush;

    unordered_map<U, size_t, HashU> index;
    index.reserve(recs.size() * 1.3);
    for (size_t i = 0; i < recs.size(); ++i) index[recs[i].key] = i;
    const auto& lookup = index;

    long long backgrounds = 0, eligible = 0, tests = 0, failures = 0, unproved_failures = 0, unproved_tests = 0;
    Failure firstUnproved;
    BranchStats total[3];
    Failure firstFailure;
    int nt = max(1, omp_get_max_threads());
    const size_t CHUNK = 25000;

    for (int W = 0; W <= MAXW; ++W) {
        auto lower = [&](int w) {
            return lower_bound(recs.begin(), recs.end(), w,
                [](const Rec& r, int x) { return r.W < x; }) - recs.begin();
        };
        size_t lo = lower(W), hi = lower(W + 1);

        for (size_t base = lo; base < hi; base += CHUNK) {
            size_t end = min(hi, base + CHUNK);
            vector<Local> local(nt);

            #pragma omp parallel for num_threads(nt) schedule(dynamic, 8)
            for (long long ii = (long long)base; ii < (long long)end; ++ii) {
                Local& s = local[omp_get_thread_num()];
                Rec& r = recs[ii];
                r.g = calculate(r);
                if (r.nf < 2) continue;
                ++s.eligible;

                auto v = allowed(r);
                if (v.empty()) continue;
                Pick pick = (RULE == 1 ? ruleM(r, v) : ruleW(r, v));
                auto it = lookup.find(childKey(r, pick.R));
                if (it == lookup.end()) continue;
                const Rec& childRec = recs[it->second];
                ++s.branch[pick.branch].backgrounds;

                int pmin = max(r.mx, 3);
                if ((pmin & 1) != (r.W & 1)) ++pmin;
                for (int p = pmin; p <= r.W; p += 2) {
                    ++s.tests;
                    I parent = gp(r, p), child = gp(childRec, p);
                    I delta = parent - child;
                    BranchStats& bs = s.branch[pick.branch];
                    {
                        int dd = (r.W - p) / 2;
                        bool unproved = p >= 6 && dd >= 8 && r.mx <= dd;
                        if (unproved) {
                            ++s.unproved_tests;
                            if (parent < child) {
                                ++s.unproved_failures;
                                Failure uf{true, r.W, p, r.key, pick.R, parent, child};
                                if (!s.ufail.yes || failureEarlier(uf, s.ufail)) s.ufail = uf;
                            }
                        }
                    }
                    if (parent < child) {
                        ++bs.failures;
                        ++s.failures;
                        Failure f{true, r.W, p, r.key, pick.R, parent, child};
                        if (failureEarlier(f, s.failure)) s.failure = f;
                    }
                    if (child == 0) ++bs.child_zero;
                    else if (child < 0) ++bs.child_negative;

                    if (child > 0)
                        updateMin(bs.min_positive_child, delta, U(child), r, p,
                                  pick.R, pick.branch, parent, child, false);
                    if (child != 0)
                        updateMin(bs.min_all_nonzero, delta, absI(child), r, p,
                                  pick.R, pick.branch, parent, child, true);
                }
            }

            for (const auto& s : local) {
                eligible += s.eligible;
                tests += s.tests;
                failures += s.failures;
                unproved_failures += s.unproved_failures; unproved_tests += s.unproved_tests;
                if (s.ufail.yes && failureEarlier(s.ufail, firstUnproved)) firstUnproved = s.ufail;
                for (int k = 0; k < 3; ++k) {
                    total[k].backgrounds += s.branch[k].backgrounds;
                    total[k].failures += s.branch[k].failures;
                    total[k].child_zero += s.branch[k].child_zero;
                    total[k].child_negative += s.branch[k].child_negative;
                    mergeMin(total[k].min_positive_child,
                             s.branch[k].min_positive_child);
                    mergeMin(total[k].min_all_nonzero,
                             s.branch[k].min_all_nonzero);
                }
                if (s.failure.yes && failureEarlier(s.failure, firstFailure))
                    firstFailure = s.failure;
            }
            backgrounds += end - base;

            double sec = chrono::duration<double>(
                chrono::steady_clock::now() - start).count();
            cout << "progress W=" << W << " processed=" << (end - base)
                 << "/" << (hi - lo) << " tests=" << tests
                 << " elapsed_s=" << fixed << setprecision(1) << sec << "\n"
                 << flush;
        }

        if (W == 32 || W == 36 || W == MAXW) {
            double sec = chrono::duration<double>(
                chrono::steady_clock::now() - start).count();
            cout << "SUMMARY upto=" << W << " backgrounds=" << backgrounds
                 << " eligible=" << eligible << " tests=" << tests
                 << " raw_failures=" << failures << " unproved_tests=" << unproved_tests
                 << " unproved_failures=" << unproved_failures
                 << " elapsed_s=" << fixed << setprecision(1) << sec << "\n";

            const char* names[3] = {"even-plus", "duplicate", "fallback"};
            for (int k = 0; k < 3; ++k) {
                cout << " branch " << names[k]
                     << " backgrounds=" << total[k].backgrounds
                     << " failures=" << total[k].failures
                     << " child_zero=" << total[k].child_zero
                     << " child_negative=" << total[k].child_negative << "\n";
                auto printMin = [&](const char* label, const MinR& m) {
                    cout << "  " << label << ": ";
                    if (!m.yes) { cout << "none\n"; return; }
                    U g = gcdU(absI(m.num), m.den);
                    I num = m.num / (I)g;
                    U den = m.den / g;
                    const Rec& br = recs[lookup.at(m.key)];
                    cout << strI(num) << "/" << strU(den)
                         << " at W=" << m.W << " B=" << show(br)
                         << " p=" << m.p << " R=(" << m.R.n << "," << m.R.m
                         << ") parent=" << strI(m.parent)
                         << " child=" << strI(m.child) << "\n";
                };
                printMin("min_ratio_child_positive", total[k].min_positive_child);
                printMin("min_ratio_all_nonzero_denominators", total[k].min_all_nonzero);
            }
            cout << flush;
        }
    }

    if (firstUnproved.yes) {
        const Rec& r = recs[lookup.at(firstUnproved.key)];
        cout << "FIRST UNPROVED-REGION failure W=" << firstUnproved.W << " B=" << show(r)
             << " p=" << firstUnproved.p << " R=(" << firstUnproved.R.n << "," << firstUnproved.R.m << ") parent="
             << strI(firstUnproved.parent) << " child=" << strI(firstUnproved.child) << "\n";
    } else cout << "no failure in the unproved region (p>=6, delta>=8, unsaturated)\n";
    if (firstFailure.yes) {
        const Rec& r = recs[lookup.at(firstFailure.key)];
        cout << "FIRST Rule1' failure W=" << firstFailure.W << " B=" << show(r)
             << " p=" << firstFailure.p << " R=(" << firstFailure.R.n << ","
             << firstFailure.R.m << ") parent=" << strI(firstFailure.parent)
             << " child=" << strI(firstFailure.child) << "\n";
    } else {
        cout << "Rule W failures: none through W=" << MAXW << "\n";
    }
    cout << flush;
}
"""

fd = os.memfd_create("fm_sec152", 0)
os.set_inheritable(fd, True)
env = dict(os.environ)
env["TMPDIR"] = "/dev/shm"
env.setdefault("OMP_NUM_THREADS", "32")

build = subprocess.run(
    ["g++", "-pipe", "-std=c++17", "-O3", "-fopenmp", "-x", "c++",
     "-o", f"/proc/self/fd/{fd}", "-"],
    input=src.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    pass_fds=(fd,), env=env
)
if build.returncode:
    print(build.stderr.decode(), file=sys.stderr)
    raise SystemExit(build.returncode)

os.fchmod(fd, 0o700)
print("compiled; launching with OMP_NUM_THREADS=" + env["OMP_NUM_THREADS"],
      flush=True)
run = subprocess.run([f"/proc/self/fd/{fd}"], pass_fds=(fd,), env=env)
if run.returncode:
    raise SystemExit(run.returncode)
