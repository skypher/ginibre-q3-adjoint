// Rule W check for one background (argv): remove per Rule W, compare g_p for all admissible p (GMP).
#include <gmpxx.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>
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
  std::map<int,int> mult, sgn; int W = 0, mx = 3;
  for (int x : B) { mult[std::abs(x)]++; sgn[std::abs(x)] = x > 0 ? 1 : -1; W += std::abs(x); mx = std::max(mx, std::abs(x)); }
  int bn = -1; long bw = -1;
  for (auto& [n, f] : mult) { long w = (long)f * n; if (w > bw || (w == bw && n > bn)) { bw = w; bn = n; } }
  std::vector<int> R;
  if (mult[bn] >= 2) R = {sgn[bn] * bn, sgn[bn] * bn};
  else { int partner = -1; for (auto it = mult.rbegin(); it != mult.rend(); ++it) if (it->first != bn && (it->first % 2) == (bn % 2)) { partner = it->first; break; }
    if (partner > 0) R = {sgn[bn] * bn, sgn[partner] * partner}; else if (bn % 2 == 0) R = {sgn[bn] * bn}; }
  if (R.empty()) { std::printf("NO RULE-W REMOVAL\n"); return 0; }
  std::vector<int> sm = B; for (int r : R) sm.erase(std::find(sm.begin(), sm.end(), r));
  std::vector<Z> g = gcol(B, W), h = gcol(sm, W); int bad = 0, tot = 0;
  for (int p = mx; p <= W; ++p) { if ((W - p) % 2) continue; ++tot; if (g[p] < h[p]) { ++bad; if (bad <= 3) std::printf("FAIL p=%d g=%s child=%s\n", p, g[p].get_str().c_str(), h[p].get_str().c_str()); } }
  std::printf("W=%d removal=(", W); for (int r : R) std::printf("%+d", r); std::printf(") p-tests=%d failures=%d\n", tot, bad);
}
