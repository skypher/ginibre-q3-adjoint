import os, subprocess, sys

CPP = r"""
#include <gmpxx.h>
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using Z = mpz_class;

static Z entry(vector<int> lab, int A, int B) {
    sort(lab.begin(), lab.end(),
         [](int x, int y) { return abs(x) > abs(y); });
    long rem = 0;
    for (int x : lab) rem += abs(x);

    auto key = [](int s, int t) {
        return ((long long)s << 32) | (unsigned)t;
    };
    unordered_map<long long, Z> cur, nxt;
    cur[key(0, 0)] = 1;

    for (int x : lab) {
        int n = abs(x), eps = x > 0 ? 1 : -1;
        rem -= n;
        nxt.clear();
        for (const auto& kv : cur) {
            int s = (int)(kv.first >> 32);
            int t = (int)(kv.first & 0xffffffff);
            const Z& v = kv.second;
            if (abs(t - B) <= rem)
                for (int z = abs(s - n); z <= s + n; z += 2)
                    if (abs(z - A) <= rem) nxt[key(z, t)] += v;
            if (abs(s - A) <= rem)
                for (int z = abs(t - n); z <= t + n; z += 2)
                    if (abs(z - B) <= rem) nxt[key(s, z)] += eps * v;
        }
        cur.swap(nxt);
    }
    auto it = cur.find(key(A, B));
    return it == cur.end() ? Z(0) : it->second;
}

int main() {
    const int k = 41, p = 43, sigma_p = -43;
    vector<int> B;
    for (int n = 1; n <= k; ++n) B.push_back(-n);
    vector<int> L = B;
    L.push_back(sigma_p);

    int W = 0, neg = 0, mx = 0, core = 0;
    for (int x : B) {
        W += abs(x);
        neg += x < 0;
        mx = max(mx, abs(x));
        core += abs(x) >= 3;
    }
    if (((neg + (sigma_p < 0)) & 1) || ((W - p) & 1) ||
        (W - p) / 2 < 8 || mx > (W - p) / 2 || core < 2)
        return 2;

    vector<int> vals, cnt;
    for (int x : L) {
        auto it = find(vals.begin(), vals.end(), x);
        if (it == vals.end()) {
            vals.push_back(x);
            cnt.push_back(1);
        } else {
            ++cnt[it - vals.begin()];
        }
    }
    vector<pair<int,int>> pairs;
    for (int i = 0; i < (int)vals.size(); ++i)
        for (int j = i; j < (int)vals.size(); ++j)
            if (i != j || cnt[i] >= 2) pairs.push_back({i, j});

    vector<Z> margins(pairs.size());
    atomic<int> done{0};
    auto start = chrono::steady_clock::now();

    #pragma omp parallel for schedule(dynamic,1)
    for (int q = 0; q < (int)pairs.size(); ++q) {
        int a = vals[pairs[q].first], b = vals[pairs[q].second];
        vector<int> C = L;
        C.erase(find(C.begin(), C.end(), a));
        C.erase(find(C.begin(), C.end(), b));
        Z z = entry(C, abs(a), abs(b));
        margins[q] = b < 0 ? -z : z;

        int n = ++done;
        if (n % 100 == 0) {
            double sec = chrono::duration<double>(
                chrono::steady_clock::now() - start).count();
            #pragma omp critical
            cerr << "progress " << n << "/" << pairs.size()
                 << " elapsed=" << fixed << setprecision(1)
                 << sec << "s\n";
        }
    }

    Z largest = margins[0];
    int qi = 0;
    for (int q = 1; q < (int)margins.size(); ++q)
        if (margins[q] > largest) {
            largest = margins[q];
            qi = q;
        }

    int ii = -1, jj = -1, best = -1, bestmax = -1;
    for (int i = 0; i < k; ++i)
        for (int j = i + 1; j < k; ++j)
            if ((abs(B[i]) + abs(B[j])) % 2 == 0) {
                int sum = abs(B[i]) + abs(B[j]);
                int m = max(abs(B[i]), abs(B[j]));
                if (sum > best || (sum == best && m > bestmax)) {
                    best = sum; bestmax = m; ii = i; jj = j;
                }
            }

    vector<int> C = B;
    int u = B[ii], v = B[jj];
    C.erase(find(C.begin(), C.end(), u));
    C.erase(find(C.begin(), C.end(), v));
    Z parent = entry(B, p, 0), child = entry(C, p, 0);
    int nonnegative = count_if(margins.begin(), margins.end(),
                               [](const Z& z) { return z >= 0; });

    cout << "factors=" << L.size() << " W=" << W
         << " sigma_p=" << sigma_p << " delta=" << (W-p)/2
         << " pair_types=" << pairs.size()
         << " nonnegative=" << nonnegative
         << " largest_margin=" << largest
         << " at=(" << vals[pairs[qi].first] << ","
         << vals[pairs[qi].second] << ")\n";
    cout << "TopPair=(" << u << "," << v << ") parent=" << parent
         << " child=" << child
         << " monotone=" << (parent >= child ? "yes" : "NO") << "\n";

    bool right_pair = (u == -41 && v == -39) ||
                      (u == -39 && v == -41);
    return pairs.size() == 861 && nonnegative == 0 &&
           largest < 0 && parent >= child && right_pair ? 0 : 3;
}
"""

fd = os.memfd_create("fm_sec177_verify", 0)
os.set_inheritable(fd, True)
env = dict(os.environ, TMPDIR="/dev/shm", OMP_NUM_THREADS="64")
build = subprocess.run(
    ["g++", "-std=c++17", "-O3", "-fopenmp", "-x", "c++",
     "-o", f"/proc/self/fd/{fd}", "-", "-lgmpxx", "-lgmp"],
    input=CPP.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    pass_fds=(fd,), env=env)
print("compile_rc", build.returncode, build.stderr.decode(), flush=True)
if build.returncode:
    raise SystemExit(build.returncode)
os.fchmod(fd, 0o700)
proc = subprocess.Popen([f"/proc/self/fd/{fd}"], pass_fds=(fd,), env=env)
try:
    rc = proc.wait(timeout=1500)
except subprocess.TimeoutExpired:
    proc.kill()
    proc.wait()
    print("timeout at 280s", flush=True)
    raise SystemExit(4)
print("verifier_rc", rc, flush=True)
raise SystemExit(rc)
