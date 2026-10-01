#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <limits>
#include <string>
#include <vector>
#include <omp.h>
#include <boost/multiprecision/cpp_int.hpp>

using Z = boost::multiprecision::cpp_int;
using namespace std;

static constexpr int MAXN = 10;
static constexpr int MAXC = 256;
static Z CH[MAXC][MAXN + 2];

struct Profile {
  array<uint16_t, MAXN> L{};
  int T = 0;
  int W = 0;
};

struct Sample {
  uint64_t key = UINT64_MAX;
  array<uint16_t, MAXN> L{};
  uint16_t M = 0;
  int N = 0;
  Z phi = 0;
};

struct Stats {
  uint64_t patternsCut = 0;
  uint64_t noflipCut = 0;
  uint64_t tpFailCut = 0;
  uint64_t noflipW40 = 0;
  uint64_t patternsAll = 0;
  uint64_t tpFailAll = 0;
  uint64_t tpCheckedAll = 0;
  uint64_t noflipNoTP = 0;
};

static void checked_inc(uint64_t& a) {
  if (a == UINT64_MAX) {
    fprintf(stderr, "counter overflow\n");
    abort();
  }
  ++a;
}

static void checked_add(uint64_t& a, uint64_t b) {
  if (UINT64_MAX - a < b) {
    fprintf(stderr, "counter merge overflow\n");
    abort();
  }
  a += b;
}

static void init_choose() {
  for (int n = 0; n < MAXC; ++n) {
    for (int k = 0; k <= MAXN + 1; ++k) {
      if (k == 0) CH[n][k] = 1;
      else if (k > n) CH[n][k] = 0;
      else CH[n][k] = CH[n - 1][k - 1] + CH[n - 1][k];
    }
  }
}

static bool residual(const array<uint16_t, MAXN>& L, int N, int T, int& Wout) {
  int p = L[N - 1];
  int W = T - p;
  Wout = W;
  if (p < 6 || ((W - p) & 1) || W - p < 16) return false;
  int delta = (W - p) / 2;
  if (L[N - 2] > delta) return false;
  int c3 = 0;
  for (int i = 0; i < N - 1; ++i) if (L[i] >= 3) ++c3;
  return c3 >= 2;
}

static void gen_rec(int N, int Tmax, int i, int lo, int sum,
                    array<uint16_t, MAXN>& cur, vector<Profile>& out) {
  if (i == N) {
    int W = 0;
    if (residual(cur, N, sum, W)) {
      Profile p;
      p.L = cur;
      p.T = sum;
      p.W = W;
      out.push_back(p);
    }
    return;
  }
  int rem = N - i;
  for (int x = lo; sum + x * rem <= Tmax; ++x) {
    cur[i] = (uint16_t)x;
    gen_rec(N, Tmax, i + 1, x, sum + x, cur, out);
  }
}

static void multiplicities_ie(const array<uint16_t, MAXN>& L, int N, vector<Z>& m) {
  int nmask = 1 << N;
  int full = nmask - 1;
  vector<int> total(nmask, 0), shift(nmask, 0), q(nmask, 0);
  m.assign(nmask, Z(0));
  for (int S = 1; S < nmask; ++S) {
    int b = __builtin_ctz((unsigned)S);
    int P = S & (S - 1);
    total[S] = total[P] + L[b];
    shift[S] = shift[P] + L[b] + 1;
    q[S] = q[P] + 1;
  }
  m[0] = 1;
  for (int S = 1; S <= full; ++S) {
    if (total[S] & 1) continue;
    int k = q[S];
    if (k == 1) continue;
    int H = total[S] / 2;
    Z acc = 0;
    for (int J = S;; J = (J - 1) & S) {
      int top = H - shift[J] + k - 1;
      if (top >= k - 2) {
        const Z& term = CH[top][k - 2];
        if (__builtin_parity((unsigned)J)) acc -= term;
        else acc += term;
      }
      if (J == 0) break;
    }
    m[S] = -acc;
    if (m[S] < 0) {
      fprintf(stderr, "negative IE multiplicity N=%d S=%d\n", N, S);
      abort();
    }
  }
}

static void walsh(vector<Z>& v) {
  int n = (int)v.size();
  for (int h = 1; h < n; h <<= 1) {
    for (int s = 0; s < n; s += 2 * h) {
      for (int j = s; j < s + h; ++j) {
        Z a = v[j], b = v[j + h];
        v[j] = a + b;
        v[j + h] = a - b;
      }
    }
  }
}

static void phi_table(const array<uint16_t, MAXN>& L, int N, vector<Z>& v) {
  int nmask = 1 << N, full = nmask - 1;
  vector<Z> m;
  multiplicities_ie(L, N, m);
  v.resize(nmask);
  for (int S = 0; S < nmask; ++S) v[S] = m[S] * m[full ^ S];
  walsh(v);
}

// Separate fusion computation, used only to check the inclusion-exclusion formula.
static Z fusion_multiplicity(const vector<int>& L) {
  vector<Z> f(1, Z(1));
  int degree = 0;
  for (int n : L) {
    vector<Z> g(degree + n + 1, Z(0));
    for (int s = 0; s <= degree; ++s) {
      if (f[s] == 0) continue;
      for (int t = abs(s - n); t <= s + n; t += 2) g[t] += f[s];
    }
    degree += n;
    f.swap(g);
  }
  return f[0];
}

static Z ie_one(const vector<int>& L) {
  int q = (int)L.size(), T = 0;
  for (int n : L) T += n;
  if (q == 0) return Z(1);
  if (T & 1) return Z(0);
  if (q == 1) return Z(0);
  int H = T / 2, full = 1 << q;
  Z acc = 0;
  for (int J = 0; J < full; ++J) {
    int shift = 0;
    for (int i = 0; i < q; ++i) if (J >> i & 1) shift += L[i] + 1;
    int top = H - shift + q - 1;
    if (top >= q - 2) {
      const Z& term = CH[top][q - 2];
      if (__builtin_parity((unsigned)J)) acc -= term;
      else acc += term;
    }
  }
  return -acc;
}

static bool pairfree(const array<uint16_t, MAXN>& L, int N, int M) {
  for (int i = 0; i + 1 < N; ++i) {
    if (L[i] == L[i + 1] && (((M >> i) ^ (M >> (i + 1))) & 1)) return false;
  }
  return true;
}

static pair<int,int> top_pair(const array<uint16_t, MAXN>& L, int N) {
  int ai = -1, bi = -1, best = -1, bestmax = -1;
  for (int i = 0; i < N - 1; ++i) for (int j = i + 1; j < N - 1; ++j) {
    int sum = L[i] + L[j];
    if (sum & 1) continue;
    int mx = max((int)L[i], (int)L[j]);
    if (sum > best || (sum == best && mx > bestmax)) {
      ai = i; bi = j; best = sum; bestmax = mx;
    }
  }
  return {ai, bi};
}

static int child_mask(const array<uint16_t, MAXN>& L, int N, int M, int ai, int bi) {
  int CM = 0, neg = 0, pos = 0;
  for (int i = 0; i < N - 1; ++i) {
    if (i == ai || i == bi) continue;
    if (M >> i & 1) { CM |= 1 << pos; ++neg; }
    ++pos;
  }
  if (neg & 1) CM |= 1 << pos;
  return CM;
}

static uint64_t splitmix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

static uint64_t sample_key(const array<uint16_t, MAXN>& L, int N, int M) {
  uint64_t h = 0x6a09e667f3bcc909ULL ^ (uint64_t)N;
  for (int i = 0; i < N; ++i) h = splitmix64(h ^ ((uint64_t)L[i] + 0x10001ULL * (i + 1)));
  return splitmix64(h ^ (uint64_t)M);
}

static void keep_sample(vector<Sample>& v, const Sample& s, int K) {
  if ((int)v.size() < K) { v.push_back(s); return; }
  int imax = 0;
  for (int i = 1; i < (int)v.size(); ++i) if (v[i].key > v[imax].key) imax = i;
  if (s.key < v[imax].key) v[imax] = s;
}

static Z direct_two_variable(const Sample& sample) {
  int N = sample.N, p = sample.L[N - 1];
  using State = pair<int,int>;
  map<State,Z> cur, nx;
  cur[{0,0}] = 1;
  for (int i = 0; i < N - 1; ++i) {
    int n = sample.L[i], eps = ((sample.M >> i) & 1) ? -1 : 1;
    nx.clear();
    for (const auto& kv : cur) {
      int sx = kv.first.first, sy = kv.first.second;
      const Z& val = kv.second;
      for (int sx2 = abs(sx - n); sx2 <= sx + n; sx2 += 2) nx[{sx2, sy}] += val;
      for (int sy2 = abs(sy - n); sy2 <= sy + n; sy2 += 2) nx[{sx, sy2}] += eps * val;
    }
    cur.swap(nx);
  }
  auto it = cur.find({p,0});
  return it == cur.end() ? Z(0) : 2 * it->second;
}

static void usage(const char* argv0) {
  printf("Usage: %s N Tmax primary_cut [threads]\n", argv0);
  printf("Exact few-factor census. N<=10; multiplicities and Walsh values use cpp_int.\n");
  printf("Example: %s 7 120 100 16\n", argv0);
}

int main(int argc, char** argv) {
  if (argc == 2 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
    usage(argv[0]); return 0;
  }
  if (argc < 4 || argc > 5) { usage(argv[0]); return 2; }
  int N = atoi(argv[1]), Tmax = atoi(argv[2]), primary = atoi(argv[3]);
  if (N < 3 || N > MAXN || Tmax < 1 || Tmax >= MAXC || primary < 1 || primary > Tmax) {
    fprintf(stderr, "invalid bounds\n"); return 2;
  }
  int threads = argc == 5 ? atoi(argv[4]) : omp_get_max_threads();
  if (threads < 1) threads = 1;
  omp_set_num_threads(threads);
  init_choose();

  // Validate the closed IE formula against fusion on a deterministic test bank.
  vector<vector<int>> bank = {{}, {1}, {1,1}, {1,2}, {1,2,3}, {2,2,3,5}, {1,3,4,6,7}};
  uint64_t checks = 0;
  for (const auto& L : bank) {
    Z a = ie_one(L), b = fusion_multiplicity(L);
    if (a != b) { cerr << "IE/fusion mismatch\n"; return 3; }
    ++checks;
  }
  cout << "SELFTEST IE_vs_fusion=" << checks << " PASS\n" << flush;

  array<uint16_t, MAXN> cur{};
  vector<Profile> profiles;
  gen_rec(N, Tmax, 0, 1, 0, cur, profiles);
  uint64_t primaryProfiles = 0;
  for (const auto& x : profiles) if (x.T <= primary) checked_inc(primaryProfiles);
  cout << "START N=" << N << " Tmax=" << Tmax << " primary=" << primary
       << " residual_profiles_total=" << profiles.size()
       << " residual_profiles_primary=" << primaryProfiles
       << " threads=" << threads << "\n" << flush;

  int np = omp_get_max_threads();
  vector<Stats> local(np);
  vector<vector<Sample>> localSamples(np);
  atomic<uint64_t> done{0};
  uint64_t total = profiles.size();
  uint64_t stride = max<uint64_t>(1, total / 100);
  auto start = chrono::steady_clock::now();
  auto lastPrint = start;

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    vector<Sample>& samples = localSamples[tid];
    #pragma omp for schedule(dynamic,8)
    for (int ix = 0; ix < (int)profiles.size(); ++ix) {
      const Profile& pr = profiles[ix];
      vector<Z> val;
      phi_table(pr.L, N, val);
      auto [ai,bi] = top_pair(pr.L, N);
      vector<Z> child;
      int p = pr.L[N - 1];
      bool childZero = true;
      if (ai >= 0) {
        int wc = 0;
        for (int i = 0; i < N - 1; ++i) if (i != ai && i != bi) wc += pr.L[i];
        childZero = p > wc;
        if (!childZero) {
          array<uint16_t, MAXN> CL{};
          int pos = 0;
          for (int i = 0; i < N - 1; ++i) if (i != ai && i != bi) CL[pos++] = pr.L[i];
          CL[pos++] = pr.L[N - 1];
          phi_table(CL, N - 1, child);
        }
      }
      for (int M = 0; M < (1 << N); ++M) {
        if (__builtin_parity((unsigned)M)) continue;
        if (!pairfree(pr.L, N, M)) continue;
        bool inCut = pr.T <= primary;
        if (inCut) checked_inc(local[tid].patternsCut);
        checked_inc(local[tid].patternsAll);

        bool noflip = true;
        Z here = val[M];
        for (int u = 0; u < N && noflip; ++u) {
          for (int w = u + 1; w < N; ++w) {
            if (val[M ^ (1 << u) ^ (1 << w)] <= here) { noflip = false; break; }
          }
        }
        if (noflip) {
          if (inCut) {
            checked_inc(local[tid].noflipCut);
            if (ai < 0) checked_inc(local[tid].noflipNoTP);
            else {
              Z cv = 0;
              if (!childZero) {
                int CM = child_mask(pr.L, N, M, ai, bi);
                cv = child[CM];
              }
              if (cv > here) checked_inc(local[tid].tpFailCut);
            }
          }
          if (pr.W <= 40) checked_inc(local[tid].noflipW40);
          Sample s;
          s.key = sample_key(pr.L, N, M);
          s.L = pr.L; s.M = (uint16_t)M; s.N = N; s.phi = here;
          keep_sample(samples, s, 16);
        }
        if (ai >= 0) {
          checked_inc(local[tid].tpCheckedAll);
          Z cv = 0;
          if (!childZero) {
            int CM = child_mask(pr.L, N, M, ai, bi);
            cv = child[CM];
          }
          if (cv > here) checked_inc(local[tid].tpFailAll);
        }
      }

      uint64_t dn = ++done;
      if (dn % stride == 0 || dn == total) {
        #pragma omp critical(progress_print)
        {
          auto now = chrono::steady_clock::now();
          if (chrono::duration<double>(now - lastPrint).count() >= 3.0 || dn == total) {
            double elapsed = chrono::duration<double>(now - start).count();
            fprintf(stderr, "progress N=%d phase=Walsh_patterns profiles=%llu/%llu elapsed_s=%.1f\n",
                    N, (unsigned long long)done.load(), (unsigned long long)total, elapsed);
            fflush(stderr);
            lastPrint = now;
          }
        }
      }
    }
  }

  Stats all;
  for (const auto& s : local) {
    checked_add(all.patternsCut, s.patternsCut);
    checked_add(all.noflipCut, s.noflipCut);
    checked_add(all.tpFailCut, s.tpFailCut);
    checked_add(all.noflipW40, s.noflipW40);
    checked_add(all.patternsAll, s.patternsAll);
    checked_add(all.tpFailAll, s.tpFailAll);
    checked_add(all.tpCheckedAll, s.tpCheckedAll);
    checked_add(all.noflipNoTP, s.noflipNoTP);
  }
  cout << "RESULT N=" << N << " Tmax=" << Tmax << " primary=" << primary
       << " residual_profiles_primary=" << primaryProfiles
       << " pairfree_even_patterns_primary=" << all.patternsCut
       << " noflip_primary=" << all.noflipCut
       << " toppair_fail_noflip_primary=" << all.tpFailCut
       << " noflip_W_le_40=" << all.noflipW40
       << " pairfree_even_patterns_all=" << all.patternsAll
       << " toppair_checked_all=" << all.tpCheckedAll
       << " toppair_fail_all=" << all.tpFailAll
       << " noflip_without_toppair_primary=" << all.noflipNoTP
       << " elapsed_s=" << chrono::duration<double>(chrono::steady_clock::now()-start).count()
       << "\n" << flush;

  vector<Sample> samples;
  for (auto& v : localSamples) for (auto& s : v) samples.push_back(s);
  sort(samples.begin(), samples.end(), [](const Sample& a, const Sample& b) { return a.key < b.key; });
  if ((int)samples.size() > 16) samples.resize(16);
  int sampleNo = 0;
  for (const auto& s : samples) {
    Z direct = direct_two_variable(s);
    bool pass = direct == s.phi;
    cout << "DIRECT_SAMPLE N=" << N << " idx=" << sampleNo++ << " key=" << s.key << " mask=" << s.M << " labels=";
    for (int i = 0; i < N; ++i) cout << ((s.M >> i & 1) ? -1 : 1) * (int)s.L[i] << (i+1==N ? "" : ",");
    cout << " walsh=" << s.phi << " two_variable=" << direct << " " << (pass ? "PASS" : "FAIL") << "\n" << flush;
    if (!pass) return 4;
  }
  cout << "DIRECT_SAMPLE_SUMMARY N=" << N << " checked=" << samples.size() << " PASS\n" << flush;
  return 0;
}
