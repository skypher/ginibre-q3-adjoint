// Exact screen of the unsaturated region: g_p = [U_p(X)] E_y prod_i (U_{n_i}(X) + eps_i U_{n_i}(y))
// for p = W - 2 delta with max n_i <= delta, p >= max n_i.  Two-spin DP with GMP.
#include <gmpxx.h>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <vector>
#include <atomic>
#include <mutex>
#include <string>
#include <chrono>
using Z = mpz_class;
struct Res { double ratio; int delta, p; std::vector<int> lab; Z g, m; };
static std::mutex mu;
int main(int argc, char** argv) {
  long samples = atol(argv[1]); int seed = atoi(argv[2]); int threads = atoi(argv[3]); int dmax = atoi(argv[4]); int mode_arg = argc > 5 ? atoi(argv[5]) : -1;
  omp_set_num_threads(threads);
  std::atomic<long> done{0}, neg{0};
  std::vector<Res> best;  // smallest ratios
  auto t0 = std::chrono::steady_clock::now();
  std::atomic<bool> fin{false};
  #pragma omp parallel
  {
    std::mt19937_64 rng(seed * 1000003ULL + omp_get_thread_num());
    #pragma omp for schedule(dynamic, 4)
    for (long s = 0; s < samples; ++s) {
      // structure: delta in [8, dmax]; k cores (labels 3..delta), b twos, ones; random signs
      int delta = 8 + rng() % (dmax - 7);
      // pair-free sampling: one sign per label value; k >= 2 non-distinguished cores (3..delta)
      int nvals = 1 + rng() % 4;                 // distinct core values
      std::vector<int> lab; int maxn = 0, W = 0, kc = 0;
      std::vector<int> valsign(delta + 2, 0);
      for (int i = 0; i < nvals; ++i) { int n = 3 + rng() % (delta - 2); if (!valsign[n]) valsign[n] = (rng() & 1) ? 1 : -1;
        int mult = 1 + rng() % 3; for (int j = 0; j < mult; ++j) { lab.push_back(valsign[n] * n); maxn = std::max(maxn, n); W += n; ++kc; } }
      if (kc < 2) { ++done; continue; }
      int s2 = (rng() & 1) ? 1 : -1, b = rng() % 16; valsign[2] = b ? s2 : 0;
      for (int i = 0; i < b; ++i) { lab.push_back(2 * s2); W += 2; maxn = std::max(maxn, 2); }
      int s1 = (rng() & 1) ? 1 : -1; valsign[1] = s1;
      int need = 2 * delta + maxn - W;
      int ones = std::max(0, need) + (int)(rng() % (3 * delta + 1));
      for (int i = 0; i < ones; ++i) lab.push_back(s1);
      W += ones; maxn = std::max(maxn, 1);
      if (delta < maxn) { ++done; continue; }
      int p = W - 2 * delta;
      if (p < maxn) { ++done; continue; }
      std::vector<int> sg(lab.size());
      for (size_t i = 0; i < lab.size(); ++i) { sg[i] = lab[i] > 0 ? 1 : -1; lab[i] = std::abs(lab[i]); }
      { int pv = W - 2 * delta; int nminus = 0; for (int x : sg) if (x < 0) ++nminus;
        int psign = (nminus % 2 == 0) ? 1 : -1;  // distinguished sign by parity
        if (pv <= delta + 1 && pv < (int)valsign.size() && valsign[pv] == -psign) { ++done; continue; } }  // keep pair-free
      // DP over (s, t) spins, s, t <= W
      int S = W + 1;
      std::vector<Z> cur(S * S), nxt(S * S);
      cur[0] = 1; int smax = 0, tmax = 0;
      for (size_t f = 0; f < lab.size(); ++f) {
        int n = lab[f]; int e = sg[f];
        for (int s2 = 0; s2 <= smax + n; ++s2) for (int t2 = 0; t2 <= tmax + n; ++t2) nxt[s2 * S + t2] = 0;
        for (int s1 = 0; s1 <= smax; ++s1) for (int t1 = 0; t1 <= tmax; ++t1) {
          const Z& v = cur[s1 * S + t1];
          if (v == 0) continue;
          for (int s2 = std::abs(s1 - n); s2 <= s1 + n; s2 += 2) nxt[s2 * S + t1] += v;
          for (int t2 = std::abs(t1 - n); t2 <= t1 + n; t2 += 2) {
            if (e > 0) nxt[s1 * S + t2] += v; else nxt[s1 * S + t2] -= v;
          }
        }
        smax += n; tmax += n;
        std::swap(cur, nxt);
      }
      Z g = cur[p * S + 0];
      // reference: m_p(all) = [U_p] prod U_{n_i}
      std::vector<Z> a(S), a2(S); a[0] = 1; int am = 0;
      for (int n : lab) { for (int s2 = 0; s2 <= am + n; ++s2) a2[s2] = 0;
        for (int s1 = 0; s1 <= am; ++s1) if (a[s1] != 0) for (int s2 = std::abs(s1 - n); s2 <= s1 + n; s2 += 2) a2[s2] += a[s1];
        am += n; std::swap(a, a2); }
      Z m = a[p];
      if (g < 0) { ++neg; std::lock_guard<std::mutex> lk(mu);
        std::printf("NEGATIVE delta %d p %d g %s labels", delta, p, g.get_str().c_str());
        for (size_t i = 0; i < lab.size(); ++i) std::printf(" %+d", lab[i] * sg[i]); std::printf("\n"); std::fflush(stdout); }
      double ratio = (m == 0) ? 1e9 : mpf_class(g).get_d() / mpf_class(m).get_d();
      { std::lock_guard<std::mutex> lk(mu);
        Res r{ratio, delta, p, {}, g, m};
        for (size_t i = 0; i < lab.size(); ++i) r.lab.push_back(lab[i] * sg[i]);
        best.push_back(r);
        if (best.size() > 400) { std::sort(best.begin(), best.end(), [](const Res& x, const Res& y){ return x.ratio < y.ratio; }); best.resize(40); } }
      long d = ++done;
      if (d % 2000 == 0) { std::lock_guard<std::mutex> lk(mu);
        double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        time_t now = time(nullptr); char ts[32]; strftime(ts, 32, "%H:%M:%S", gmtime(&now));
        std::printf("%s UTC progress %ld/%ld samples, negatives %ld, %.1f s\n", ts, d, samples, neg.load(), sec); std::fflush(stdout); }
    }
  }
  std::sort(best.begin(), best.end(), [](const Res& x, const Res& y){ return x.ratio < y.ratio; });
  std::printf("DONE samples %ld negatives %ld\nsmallest g_p / m_p(all):\n", done.load(), neg.load());
  for (size_t i = 0; i < std::min<size_t>(12, best.size()); ++i) {
    auto& r = best[i]; std::printf("ratio %.3e delta %d p %d g %s labels", r.ratio, r.delta, r.p, r.g.get_str().c_str());
    for (int x : r.lab) std::printf(" %+d", x); std::printf("\n"); }
}
