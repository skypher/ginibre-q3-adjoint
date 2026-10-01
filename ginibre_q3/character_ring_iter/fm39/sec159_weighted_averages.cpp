// Averaging test for (D): is sum_R [g_p(B) - g_p(B - R)] >= 0 over all even-weight removals R (multiset-distinct)?
// Also tests the sum over pair removals only.  GMP, OpenMP over removals.
#include <gmpxx.h>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <set>
#include <algorithm>
using Z = mpz_class;
static std::vector<Z> gcol(const std::vector<int>& lab, int W) {
  int S = W + 1; std::vector<Z> cur(S * S), nxt(S * S); cur[0] = 1; int sm = 0, tm = 0;
  for (int x : lab) { int n = std::abs(x), e = x > 0 ? 1 : -1;
    for (int a = 0; a <= sm + n; ++a) for (int b = 0; b <= tm + n; ++b) nxt[a * S + b] = 0;
    for (int a = 0; a <= sm; ++a) for (int b = 0; b <= tm; ++b) { const Z& v = cur[a * S + b]; if (v == 0) continue;
      for (int t = std::abs(a - n); t <= a + n; t += 2) nxt[t * S + b] += v;
      for (int t = std::abs(b - n); t <= b + n; t += 2) { if (e > 0) nxt[a * S + t] += v; else nxt[a * S + t] -= v; } }
    sm += n; tm += n; std::swap(cur, nxt); }
  std::vector<Z> g(S); for (int p = 0; p <= std::min(sm, W); ++p) g[p] = cur[p * S]; return g; }
int main(int argc, char** argv) {
  std::vector<int> B; for (int i = 1; i < argc; ++i) B.push_back(atoi(argv[i]));
  int W = 0, mx = 3; for (int x : B) { W += std::abs(x); mx = std::max(mx, std::abs(x)); }
  std::vector<Z> g = gcol(B, W);
  // index-level removals (with multiplicity): all pairs i<j with even total, all single even labels
  std::vector<std::vector<int>> R;
  for (size_t i = 0; i < B.size(); ++i) { if (std::abs(B[i]) % 2 == 0) R.push_back({(int)i});
    for (size_t j = i + 1; j < B.size(); ++j) if ((std::abs(B[i]) + std::abs(B[j])) % 2 == 0) R.push_back({(int)i, (int)j}); }
  // cache by multiset
  std::vector<std::vector<int>> keys(R.size());
  for (size_t k = 0; k < R.size(); ++k) { for (int i : R[k]) keys[k].push_back(B[i]); std::sort(keys[k].begin(), keys[k].end()); }
  std::vector<std::vector<int>> uk = keys; std::sort(uk.begin(), uk.end()); uk.erase(std::unique(uk.begin(), uk.end()), uk.end());
  std::vector<std::vector<Z>> gs(uk.size());
  #pragma omp parallel for schedule(dynamic)
  for (size_t k = 0; k < uk.size(); ++k) { std::vector<int> sm = B; for (int r : uk[k]) sm.erase(std::find(sm.begin(), sm.end(), r)); gs[k] = gcol(sm, W); }
  int neg_all = 0, neg_pairs = 0, neg_multi = 0, tot = 0, neg_dim = 0, neg_lab = 0, neg_dim2 = 0;
  for (int p = mx; p <= W; ++p) { if ((W - p) % 2) continue; ++tot;
    Z sum_index = 0, sum_pairs = 0, sum_distinct = 0, sum_dim = 0, sum_lab = 0, sum_dim2 = 0;
    for (size_t k = 0; k < R.size(); ++k) { size_t u = std::lower_bound(uk.begin(), uk.end(), keys[k]) - uk.begin(); Z d = g[p] - gs[u][p]; sum_index += d; if (R[k].size() == 2) sum_pairs += d;
      long wd = 1, wl = 1; for (int i : R[k]) { wd *= (std::abs(B[i]) + 1); wl *= std::abs(B[i]); }
      sum_dim += d * wd; sum_lab += d * wl; sum_dim2 += d * (wd * wd); }
    for (size_t u = 0; u < uk.size(); ++u) sum_distinct += g[p] - gs[u][p];
    if (sum_index < 0) ++neg_all; if (sum_pairs < 0) ++neg_pairs; if (sum_distinct < 0) ++neg_multi;
    if (sum_dim < 0) ++neg_dim; if (sum_lab < 0) ++neg_lab; if (sum_dim2 < 0) ++neg_dim2; }
  std::printf("W=%d p=%d neg: uniform=%d pairs=%d multiset=%d dim-weighted=%d label-weighted=%d dim^2-weighted=%d\n", W, tot, neg_all, neg_pairs, neg_multi, neg_dim, neg_lab, neg_dim2);
}
