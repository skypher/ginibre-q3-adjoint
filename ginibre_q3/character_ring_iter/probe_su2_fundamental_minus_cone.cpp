// Exact probe of the fundamental-minus reduction D_q = D_1 H_q in the
// doubled SU(2) (or SU(2)_k) representation ring.
//   S_p = V_p x 1 + 1 x V_p,  D_q = V_q x 1 - 1 x V_q,
//   H_q = sum_{m+n=q-1} V_m x V_n.
// A doubled character is a matrix F[a][b] = coefficient of V_a x V_b.
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <chrono>
using boost::multiprecision::cpp_int;
using Mat = std::vector<std::vector<cpp_int>>;

static int LEVEL = -1;  // -1 = ordinary SU(2)
static int DIM = 0;     // labels 0..DIM-1

static bool fuse_ok(int a, int p, int c) {
  if (LEVEL >= 0 && (a > LEVEL || p > LEVEL || c > LEVEL)) return false;
  if (LEVEL >= 0 && c > 2 * LEVEL - a - p) return false;
  return true;
}
static Mat zero() { return Mat(DIM, std::vector<cpp_int>(DIM)); }
// multiply by V_p on side (0 = x, 1 = y) with coefficient s
static void mul_side(const Mat& F, int p, int side, long s, Mat& out) {
  for (int a = 0; a < DIM; ++a)
    for (int b = 0; b < DIM; ++b) {
      if (F[a][b] == 0) continue;
      int u = side == 0 ? a : b;
      for (int c = (u > p ? u - p : p - u); c <= u + p; c += 2) {
        if (!fuse_ok(u, p, c)) continue;
        if (c >= DIM) { std::fprintf(stderr, "DIM overflow\n"); std::exit(2); }
        if (side == 0) out[c][b] += s * F[a][b]; else out[a][c] += s * F[a][b];
      }
    }
}
static Mat mulS(const Mat& F, int p) { Mat o = zero(); mul_side(F, p, 0, 1, o); mul_side(F, p, 1, 1, o); return o; }
static Mat mulD(const Mat& F, int p) { Mat o = zero(); mul_side(F, p, 0, 1, o); mul_side(F, p, 1, -1, o); return o; }
static Mat mulH(const Mat& F, int q) {
  Mat o = zero();
  for (int m = 0; m <= q - 1; ++m) {
    int n = q - 1 - m;
    if (LEVEL >= 0 && (m > LEVEL || n > LEVEL)) continue;
    Mat t = zero(); mul_side(F, m, 0, 1, t);
    mul_side(t, n, 1, 1, o);
  }
  return o;
}
static Mat unit() { Mat o = zero(); o[0][0] = 1; return o; }

static void usage() {
  std::printf(
    "usage: probe_su2_fundamental_minus_cone --identity QMAX [--level k]\n"
    "       probe_su2_fundamental_minus_cone --cone L N M [--level k] [--onlyH]\n"
    "  --identity: check F*D_q == F*D_1*H_q for random-ish F and q<=QMAX.\n"
    "  --cone: all products F of <=N generators from {S_p:1<=p<=L} u {H_q:2<=q<=L},\n"
    "          times D_1^m for 0<=m<=M; checks every coefficient of the partial\n"
    "          character G(x)=sum_a R[a][0] chi_a is >=0. Reports failures split by\n"
    "          whether #H<=m (equivalent to an original word) or #H>m (enlarged cone).\n"
    "  --control: also adjoin X1 = V_1 x V_1 + 1 (negative control; failures expected).\n");
}

struct Gen { char kind; int lab; };
static std::vector<Gen> gens;
static long long nodes = 0, checks = 0, fail_orig = 0, fail_enl = 0;
static int NMAX, MMAX;
static std::vector<int> stackIdx;
static auto T0 = std::chrono::steady_clock::now();

static void report_fail(const Mat& R, int m, int nH, int a) {
  std::printf("FAIL %s m=%d F=[", nH <= m ? "ORIG" : "ENL", m);
  for (int i : stackIdx) std::printf(" %c%d", gens[i].kind, gens[i].lab);
  std::printf(" ] coeff a=%d value=%s\n", a, R[a][0].str().c_str());
}

static void dfs(const Mat& F, int start, int depth, int nH) {
  ++nodes;
  Mat R = F;
  for (int m = 0; m <= MMAX; ++m) {
    if (m > 0) R = mulD(R, 1);
    ++checks;
    for (int a = 0; a < DIM; ++a)
      if (R[a][0] < 0) {
        if (nH <= m) { if (fail_orig++ < 20) report_fail(R, m, nH, a); }
        else { if (fail_enl++ < 20) report_fail(R, m, nH, a); }
        break;
      }
  }
  if (nodes % 20000 == 0) {
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
    std::printf("progress nodes=%lld checks=%lld depth=%d fail_orig=%lld fail_enl=%lld elapsed=%.1fs\n",
                nodes, checks, depth, fail_orig, fail_enl, el);
    std::fflush(stdout);
  }
  if (depth == NMAX) return;
  for (int i = start; i < (int)gens.size(); ++i) {
    stackIdx.push_back(i);
    Mat G;
    if (gens[i].kind == 'S') G = mulS(F, gens[i].lab);
    else if (gens[i].kind == 'H') G = mulH(F, gens[i].lab);
    else { Mat t = zero(); mul_side(F, 1, 0, 1, t); G = zero(); mul_side(t, 1, 1, 1, G);
           for (int a = 0; a < DIM; ++a) for (int b = 0; b < DIM; ++b) G[a][b] += F[a][b]; }
    dfs(G, i, depth + 1, nH + (gens[i].kind == 'H'));
    stackIdx.pop_back();
  }
}

int main(int argc, char** argv) {
  if (argc < 2 || !std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help")) { usage(); return 0; }
  for (int i = 1; i < argc; ++i) if (!std::strcmp(argv[i], "--level")) LEVEL = std::atoi(argv[i + 1]);
  if (!std::strcmp(argv[1], "--identity")) {
    int Q = std::atoi(argv[2]);
    DIM = 4 * Q + 8;
    long bad = 0, tot = 0;
    // base F: a few products to make F nontrivial
    std::vector<Mat> bases = {unit(), mulS(unit(), 1), mulS(mulS(unit(), 2), 3), mulD(mulS(unit(), 1), 2)};
    for (auto& B : bases)
      for (int q = 1; q <= Q; ++q) {
        if (LEVEL >= 0 && q > LEVEL) continue;
        Mat L = mulD(B, q), Rr = mulH(mulD(B, 1), q);
        ++tot; if (L != Rr) { ++bad; std::printf("MISMATCH q=%d\n", q); }
      }
    std::printf("identity checks=%ld mismatches=%ld level=%d result=%s\n", tot, bad, LEVEL, bad ? "FAIL" : "PASS");
    return bad ? 1 : 0;
  }
  if (!std::strcmp(argv[1], "--cone")) {
    int L = std::atoi(argv[2]); NMAX = std::atoi(argv[3]); MMAX = std::atoi(argv[4]);
    bool onlyH = false; for (int i = 1; i < argc; ++i) if (!std::strcmp(argv[i], "--onlyH")) onlyH = true;
    int LS = (LEVEL >= 0 && L > LEVEL) ? LEVEL : L;
    for (int p = 1; p <= LS; ++p) if (!onlyH) gens.push_back({'S', p});
    for (int q = 2; q <= LS; ++q) gens.push_back({'H', q});
    for (int i = 1; i < argc; ++i) if (!std::strcmp(argv[i], "--control")) gens.push_back({'X', 1});
    DIM = NMAX * L + MMAX + 4;
    if (LEVEL >= 0 && DIM > LEVEL + 1) DIM = LEVEL + 1 + 0, DIM = std::max(DIM, LEVEL + 1);
    if (LEVEL >= 0) DIM = LEVEL + 1;
    dfs(unit(), 0, 0, 0);
    std::printf("done L=%d N=%d M=%d level=%d nodes=%lld checks=%lld fail_orig=%lld fail_enl=%lld\n",
                L, NMAX, MMAX, LEVEL, nodes, checks, fail_orig, fail_enl);
    return 0;
  }
  usage(); return 1;
}
