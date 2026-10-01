// (D) checker: for a signed background B (argv) and every admissible p, is there an even-weight removal R
// (one even label, or two labels of even total weight) with g_p(B) >= g_p(B - R)?  GMP, OpenMP over removals.
#include <gmpxx.h>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <set>
#include <algorithm>
using Z = mpz_class;
static std::vector<Z> gcol(const std::vector<int>& lab, int W) {   // returns g_p for p = 0..W
  int S = W + 1; std::vector<Z> cur(S * S), nxt(S * S); cur[0] = 1; int sm = 0, tm = 0;
  for (int x : lab) {
    int n = std::abs(x), e = x > 0 ? 1 : -1;
    for (int a = 0; a <= sm + n; ++a) for (int b = 0; b <= tm + n; ++b) nxt[a * S + b] = 0;
    for (int a = 0; a <= sm; ++a) for (int b = 0; b <= tm; ++b) {
      const Z& v = cur[a * S + b]; if (v == 0) continue;
      for (int t = std::abs(a - n); t <= a + n; t += 2) nxt[t * S + b] += v;
      for (int t = std::abs(b - n); t <= b + n; t += 2) { if (e > 0) nxt[a * S + t] += v; else nxt[a * S + t] -= v; }
    }
    sm += n; tm += n; std::swap(cur, nxt);
  }
  std::vector<Z> g(S); for (int p = 0; p <= std::min(sm, W); ++p) g[p] = cur[p * S]; return g;
}
int main(int argc, char** argv) {
  std::vector<int> B; for (int i = 1; i < argc; ++i) B.push_back(atoi(argv[i]));
  int W = 0, mx = 3; for (int x : B) { W += std::abs(x); mx = std::max(mx, std::abs(x)); }
  std::vector<Z> g = gcol(B, W);
  std::set<std::vector<int>> rems;
  for (size_t i = 0; i < B.size(); ++i) {
    if (std::abs(B[i]) % 2 == 0) rems.insert({B[i]});
    for (size_t j = i + 1; j < B.size(); ++j) if ((std::abs(B[i]) + std::abs(B[j])) % 2 == 0) { std::vector<int> r = {B[i], B[j]}; std::sort(r.begin(), r.end()); rems.insert(r); }
  }
  std::vector<std::vector<int>> R(rems.begin(), rems.end());
  std::vector<std::vector<Z>> gs(R.size());
  #pragma omp parallel for schedule(dynamic)
  for (size_t k = 0; k < R.size(); ++k) {
    std::vector<int> sm = B; for (int r : R[k]) sm.erase(std::find(sm.begin(), sm.end(), r));
    gs[k] = gcol(sm, W);
  }
  int bad = 0, tot = 0;
  for (int p = mx; p <= W; ++p) {
    if ((W - p) % 2) continue; ++tot;
    int okc = 0; for (size_t k = 0; k < R.size(); ++k) if (g[p] >= gs[k][p]) ++okc;
    if (!okc) { ++bad; std::printf("NO MONOTONE REMOVAL at p=%d g=%s\n", p, g[p].get_str().c_str()); }
  }
  std::printf("W=%d labels=%zu removals=%zu admissible p=%d stuck=%d\n", W, B.size(), R.size(), tot, bad);
  std::printf("removals monotone for ALL admissible p:");
  for (size_t k = 0; k < R.size(); ++k) { bool all = true;
    for (int p = mx; p <= W; ++p) { if ((W - p) % 2) continue; if (g[p] < gs[k][p]) { all = false; break; } }
    if (all) { std::printf(" ("); for (int r : R[k]) std::printf("%+d", r); std::printf(")"); } }
  std::printf("\n");
}
