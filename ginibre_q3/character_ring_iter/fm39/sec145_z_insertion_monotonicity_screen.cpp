// Stress test of +2 insertion monotonicity: g_p(B + {+2}) >= g_p(B) for all p >= max(labels, 3).
#include <gmpxx.h>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
#include <atomic>
#include <mutex>
#include <chrono>
#include <ctime>
using Z = mpz_class;
static std::mutex mu;
static void stepf(std::vector<Z>& cur, std::vector<Z>& nxt, int S, int& smax, int& tmax, int n, int e) {
  for (int s2 = 0; s2 <= smax + n; ++s2) for (int t2 = 0; t2 <= tmax + n; ++t2) nxt[s2 * S + t2] = 0;
  for (int s1 = 0; s1 <= smax; ++s1) for (int t1 = 0; t1 <= tmax; ++t1) {
    const Z& v = cur[s1 * S + t1]; if (v == 0) continue;
    for (int s2 = std::abs(s1 - n); s2 <= s1 + n; s2 += 2) nxt[s2 * S + t1] += v;
    for (int t2 = std::abs(t1 - n); t2 <= t1 + n; t2 += 2) { if (e > 0) nxt[s1 * S + t2] += v; else nxt[s1 * S + t2] -= v; }
  }
  smax += n; tmax += n; std::swap(cur, nxt);
}
int main(int argc, char** argv) {
  long samples = atol(argv[1]); int seed = atoi(argv[2]); int threads = atoi(argv[3]);
  int kmax = atoi(argv[4]), lmax = atoi(argv[5]), omax = atoi(argv[6]), bmax = atoi(argv[7]);
  omp_set_num_threads(threads);
  std::atomic<long> done{0}, checks{0}, viol{0};
  double worst = 1e300; std::vector<int> wl; int wp = 0;
  auto t0 = std::chrono::steady_clock::now();
  #pragma omp parallel
  {
    std::mt19937_64 rng(seed * 7919ULL + omp_get_thread_num());
    #pragma omp for schedule(dynamic, 2)
    for (long s = 0; s < samples; ++s) {
      std::vector<int> lab;
      int k = rng() % (kmax + 1);
      for (int i = 0; i < k; ++i) lab.push_back(((rng() & 1) ? 1 : -1) * (3 + (int)(rng() % (lmax - 2))));
      int o = rng() % (omax + 1); double bias = (rng() % 1000) / 1000.0;
      for (int i = 0; i < o; ++i) lab.push_back(((rng() % 1000) / 1000.0 < bias) ? 1 : -1);
      int b = rng() % (bmax + 1); for (int i = 0; i < b; ++i) lab.push_back(2);
      int W = 2; int mx = 3; for (int x : lab) { W += std::abs(x); mx = std::max(mx, std::abs(x)); }
      int S = W + 1; std::vector<Z> cur(S * S), nxt(S * S); cur[0] = 1; int smax = 0, tmax = 0;
      for (int x : lab) stepf(cur, nxt, S, smax, tmax, std::abs(x), x > 0 ? 1 : -1);
      std::vector<Z> g0(S); for (int p = 0; p <= smax; ++p) g0[p] = cur[p * S];
      stepf(cur, nxt, S, smax, tmax, 2, 1);
      for (int p = mx; p <= smax; ++p) {
        if ((smax - p) % 2) continue;
        Z d = cur[p * S] - g0[p]; ++checks;
        if (d < 0) { ++viol; std::lock_guard<std::mutex> lk(mu);
          std::printf("VIOLATION p %d diff %s labels", p, d.get_str().c_str()); for (int x : lab) std::printf(" %+d", x); std::printf("\n"); std::fflush(stdout); }
        else if (g0[p] > 0) { double r = mpf_class(d).get_d() / mpf_class(g0[p]).get_d();
          if (r < worst) { std::lock_guard<std::mutex> lk(mu); if (r < worst) { worst = r; wl = lab; wp = p; } } }
      }
      long dn = ++done;
      if (dn % 5000 == 0) { std::lock_guard<std::mutex> lk(mu); time_t now = time(nullptr); char ts[32]; strftime(ts, 32, "%H:%M:%S", gmtime(&now));
        std::printf("%s UTC progress %ld/%ld backgrounds, %ld checks, %ld violations, %.0f s\n", ts, dn, samples, checks.load(), viol.load(),
          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()); std::fflush(stdout); }
    }
  }
  std::printf("DONE backgrounds %ld checks %ld violations %ld; smallest (g'_p - g_p)/g_p = %.3e at p %d labels", done.load(), checks.load(), viol.load(), worst, wp);
  for (int x : wl) std::printf(" %+d", x); std::printf("\n");
}
