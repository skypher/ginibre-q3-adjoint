// g_p(B) for one p with reachability pruning: states (s,t) with |s - p| <= rem and t <= rem.
// Usage: gp_single p  B...  --  R1...  --  R2...   (prints g_p(B) and g_p(B - R) for each removal)
#include <gmpxx.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
using Z = mpz_class;
static Z gp(std::vector<int> lab, int p) {
  std::sort(lab.begin(), lab.end(), [](int a, int b) { return std::abs(a) > std::abs(b); });
  long rem = 0;
  for (int x : lab) rem += std::abs(x);
  auto key = [](int s, int t) { return ((long long)s << 32) | (unsigned)t; };
  std::unordered_map<long long, Z> cur, nx;
  cur[key(0, 0)] = 1;
  for (int x : lab) {
    int n = std::abs(x), e = x > 0 ? 1 : -1;
    rem -= n;
    nx.clear();
    for (auto& kv : cur) {
      int s = (int)(kv.first >> 32), t = (int)(kv.first & 0xffffffff);
      const Z& v = kv.second;
      if (t <= rem)
        for (int s2 = std::abs(s - n); s2 <= s + n; s2 += 2)
          if (std::abs(s2 - p) <= rem) nx[key(s2, t)] += v;
      if (std::abs(s - p) <= rem)
        for (int t2 = std::abs(t - n); t2 <= t + n; t2 += 2)
          if (t2 <= rem) { if (e > 0) nx[key(s, t2)] += v; else nx[key(s, t2)] -= v; }
    }
    cur.swap(nx);
  }
  auto it = cur.find(key(p, 0));
  return it == cur.end() ? Z(0) : it->second;
}
int main(int argc, char** argv) {
  int p = atoi(argv[1]);
  std::vector<std::vector<int>> parts(1);
  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == "--") { parts.emplace_back(); continue; }
    parts.back().push_back(atoi(argv[i]));
  }
  const std::vector<int>& B = parts[0];
  Z g = gp(B, p);
  std::printf("g_p(B) = %s\n", g.get_str().c_str());
  for (size_t k = 1; k < parts.size(); ++k) {
    std::vector<int> sm = B;
    for (int r : parts[k]) sm.erase(std::find(sm.begin(), sm.end(), r));
    Z h = gp(sm, p);
    std::printf("removal (");
    for (int r : parts[k]) std::printf("%+d", r);
    std::printf("): g_p(B-R) = %s  %s\n", h.get_str().c_str(), g >= h ? "monotone" : "INCREASES");
  }
}
