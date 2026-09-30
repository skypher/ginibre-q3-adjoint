"""FM-SEC96 (luna_max_vesta): invariant subspaces U_S; meet-identity failure; ranks of natural linear maps E_- -> E_+."""
from itertools import product
from math import comb
import sympy as sp

def inv_basis(ns):
    ns = tuple(ns)
    D = sum(ns)
    if D % 2:
        return (), ()
    weight0 = tuple(k for k in product(*(range(n + 1) for n in ns))
                    if sum(k) == D // 2)
    weight2 = (tuple(k for k in product(*(range(n + 1) for n in ns))
                     if sum(k) == D // 2 - 1) if D else ())
    row = {k: i for i, k in enumerate(weight2)}
    E = sp.zeros(len(weight2), len(weight0))
    for j, k in enumerate(weight0):
        for i, ki in enumerate(k):
            if ki:
                E[row[k[:i] + (ki - 1,) + k[i + 1:]], j] += ki
    return weight0, tuple(v.as_immutable() for v in E.nullspace())

def cut_basis(ns, mask, full_weight0):
    L = len(ns)
    A = tuple(i for i in range(L) if mask >> i & 1)
    B = tuple(i for i in range(L) if not (mask >> i & 1))
    ta, ba = inv_basis(tuple(ns[i] for i in A))
    tb, bb = inv_basis(tuple(ns[i] for i in B))
    ix = {k: j for j, k in enumerate(full_weight0)}
    cols = []
    for va in ba:
        for vb in bb:
            v = sp.zeros(len(full_weight0), 1)
            for a, ca in enumerate(ta):
                for b, cb in enumerate(tb):
                    k = [0] * L
                    for z, i in enumerate(A):
                        k[i] = ca[z]
                    for z, i in enumerate(B):
                        k[i] = cb[z]
                    v[ix[tuple(k)]] += va[a] * vb[b]
            cols.append(v.as_immutable())
    return sp.Matrix.hstack(*cols) if cols else sp.zeros(len(full_weight0), 0)

def model(ns):
    ns = tuple(ns)
    L, D = len(ns), sum(ns)
    assert D % 2 == 0
    full, _ = inv_basis(ns)
    high = (tuple(k for k in product(*(range(n + 1) for n in ns))
                  if sum(k) == D // 2 - 1) if D else ())
    hi = {k: i for i, k in enumerate(high)}
    E = sp.zeros(len(high), len(full))
    for j, k in enumerate(full):
        for i, ki in enumerate(k):
            if ki:
                E[hi[k[:i] + (ki - 1,) + k[i + 1:]], j] += ki
    Hbasis = E.nullspace()
    H = sp.Matrix.hstack(*Hbasis) if Hbasis else sp.zeros(len(full), 0)
    W = sp.diag(*[
        sp.prod(sp.Rational(1, comb(n, ki)) for n, ki in zip(ns, k))
        for k in full
    ])
    U = {s: cut_basis(ns, s, full) for s in range(1 << (L - 1))}
    for s, B in U.items():
        A = tuple(i for i in range(L) if s >> i & 1)
        C = tuple(i for i in range(L) if not (s >> i & 1))
        assert B.cols == (
            len(inv_basis(tuple(ns[i] for i in A))[1]) *
            len(inv_basis(tuple(ns[i] for i in C))[1])
        )
        assert (E * B).is_zero_matrix
    return full, W, H, U

def canon(mask, L):
    return mask if not (mask >> (L - 1) & 1) else mask ^ ((1 << L) - 1)

def slots(ns, T, U):
    L = len(ns)
    reps = range(1 << (L - 1))
    neg = [s for s in reps
           if (s & T).bit_count() % 2 and U[s].cols]
    pos = [s for s in reps
           if not ((s & T).bit_count() % 2) and U[s].cols]
    return pos, neg, sum(U[s].cols for s in pos), sum(U[s].cols for s in neg)

def projection_matrix(ns, T, U, W, kind):
    L = len(ns)
    pos, neg, dp, dm = slots(ns, T, U)
    prow, ncol = {}, {}
    q = 0
    for s in pos:
        prow[s] = q
        q += U[s].cols
    q = 0
    for s in neg:
        ncol[s] = q
        q += U[s].cols
    M = sp.zeros(dp, dm)
    for s in neg:
        for i in range(L):
            if not (T >> i & 1):
                continue
            if kind == "single":
                flips = [1 << i]
            elif kind == "pair":
                flips = [(1 << i) | (1 << j)
                         for j in range(L) if not (T >> j & 1)]
            elif kind == "both":
                flips = ([1 << i] +
                         [(1 << i) | (1 << j)
                          for j in range(L) if not (T >> j & 1)])
            else:
                raise ValueError(kind)
            for flip in flips:
                t = canon(s ^ flip, L)
                if t not in prow:
                    continue
                B, A = U[t], U[s]
                C = (B.T * W * B).inv() * B.T * W * A
                r, c = prow[t], ncol[s]
                M[r:r+B.cols, c:c+A.cols] += C
    return dp, dm, M

def multiplication_matrix(ns, T, U):
    L = len(ns)
    neg = [s for s in range(1 << (L - 1))
           if (s & T).bit_count() % 2 and U[s].cols]
    M = (sp.Matrix.hstack(*(U[s] for s in neg)) if neg
         else sp.zeros(next(iter(U.values())).rows, 0))
    return neg, M

for ns, T in [((1,1,1,1), 3), ((1,1,2,2), 5), ((1,5,2,2), 3)]:
    full, W, H, U = model(ns)
    pos, neg, dp, dm = slots(ns, T, U)
    print("boundary", ns, "T=", T, "dimH/E+/E-=", H.cols, dp, dm)
    for kind in ("single", "pair", "both"):
        _, _, M = projection_matrix(ns, T, U, W, kind)
        print(" ", kind, "rank=", M.rank())
    nslots, Mu = multiplication_matrix(ns, T, U)
    print("  mu-minus rank=", Mu.rank())

ns, T = (1,1,1,1,2), 3
full, W, H, U = model(ns)
pos, neg, dp, dm = slots(ns, T, U)
_, _, Ppair = projection_matrix(ns, T, U, W, "pair")
nslots, Mu = multiplication_matrix(ns, T, U)
print("first-loss", ns, "T=", T, "E+/E-=", dp, dm,
      "pair-rank=", Ppair.rank(), "pair-kernel=", Ppair.nullspace())
print("  negative representatives=", nslots, "mu-rank=", Mu.rank(),
      "mu-kernel=", Mu.nullspace())

ns = (1,1,1,1,1,1)
full, W, H, U = model(ns)
S, Sprime, meet = 15, 51, 3
A, B = U[S], U[canon(Sprime, len(ns))]
C = U[meet]
rankAB = sp.Matrix.hstack(A, B).rank()
print("meet-test", "dim(U_S),dim(U_Sprime),dim(U_meet)=",
      A.cols, B.cols, C.cols, "rank([U_S U_Sprime])=", rankAB,
      "intersection-dim=", A.cols + B.cols - rankAB)
