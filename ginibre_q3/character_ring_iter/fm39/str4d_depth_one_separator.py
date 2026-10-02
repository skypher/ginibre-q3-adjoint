import argparse
import collections
import datetime
import functools
import heapq
import math
from fractions import Fraction as F

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack

ap = argparse.ArgumentParser(
    description="Exact component-depth-one separator for a no-flip root."
)
ap.add_argument("--time-limit", type=float, default=120.0)
args = ap.parse_args()

def stamp():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%H:%M:%S UTC")

L = (1, 4, 5, 6, 6, 8, 9, 15)
mask = (1 << 2) | (1 << 6)

@functools.lru_cache(None)
def cg(a, b):
    return tuple(range(abs(a - b), a + b + 1, 2))

def cuts(Q):
    for s in range(1 << len(Q)):
        yield (
            tuple(a for i, a in enumerate(Q) if (s >> i) & 1),
            tuple(a for i, a in enumerate(Q) if not ((s >> i) & 1)),
        )

def fuse(Q, i, j):
    R = Q[:i] + Q[i + 1:j] + Q[j + 1:]
    return [
        tuple(sorted(R + ((c,) if c else ())))
        for c in cg(Q[i], Q[j])
    ]

def word_closure(L):
    d = {A: 0 for A, _ in cuts(L)}
    for Q in list(d):
        for i in range(len(Q)):
            for j in range(i + 1, len(Q)):
                for R in fuse(Q, i, j):
                    d.setdefault(R, 1)
    return d

def sparse(rows, n):
    rr, cc, vv = [], [], []
    for i, row in enumerate(rows):
        for j, v in row.items():
            if v:
                rr.append(i)
                cc.append(j)
                vv.append(int(v))
    return coo_matrix(
        (vv, (rr, cc)), shape=(len(rows), n), dtype=np.int64
    ).tocsr()

def mat_vec(M, x):
    return [
        sum(
            int(v) * x[j]
            for j, v in zip(
                M.indices[M.indptr[i]:M.indptr[i + 1]],
                M.data[M.indptr[i]:M.indptr[i + 1]],
            )
        )
        for i in range(M.shape[0])
    ]

def exact_solve(mat, right, approximate):
    mat = mat.tocsr()
    equations, rhs = [], []
    incidence = [set() for _ in approximate]

    for i in range(mat.shape[0]):
        row = {
            int(j): F(int(v))
            for j, v in zip(
                mat.indices[mat.indptr[i]:mat.indptr[i + 1]],
                mat.data[mat.indptr[i]:mat.indptr[i + 1]],
            )
            if v
        }
        b = F(int(right[i]))
        if not row:
            assert b == 0
            continue
        h = len(equations)
        equations.append(row)
        rhs.append(b)
        for j in row:
            incidence[j].add(h)

    heap = [(len(row), i) for i, row in enumerate(equations)]
    heapq.heapify(heap)
    saved = []

    while heap:
        width, i = heapq.heappop(heap)
        row = equations[i]
        if row is None or width != len(row):
            continue
        if not row:
            assert rhs[i] == 0
            equations[i] = None
            continue

        pivot = min(
            row,
            key=lambda j: (len(incidence[j]), abs(row[j]) != 1, j),
        )
        a = row[pivot]
        rest = {j: v / a for j, v in row.items() if j != pivot}
        b = rhs[i] / a
        saved.append((pivot, rest, b))

        for j in row:
            incidence[j].discard(i)
        equations[i] = None

        for h in list(incidence[pivot]):
            rr = equations[h]
            a = rr.pop(pivot)
            incidence[pivot].discard(h)
            rhs[h] -= a * b
            for j, v in rest.items():
                z = rr.get(j, F(0)) - a * v
                if z:
                    if j not in rr:
                        incidence[j].add(h)
                    rr[j] = z
                elif j in rr:
                    del rr[j]
                    incidence[j].discard(h)
            heapq.heappush(heap, (len(rr), h))

        if len(saved) % 250 == 0:
            print("EXACT_PIVOT", len(saved), stamp(), flush=True)

    x = [F(float(v)).limit_denominator(10**9) for v in approximate]
    for pivot, row, b in reversed(saved):
        x[pivot] = b - sum(v * x[j] for j, v in row.items())
    return x

def mult(Q):
    d = {0: 1}
    for n in Q:
        out = collections.defaultdict(int)
        for a, v in d.items():
            for c in cg(a, n):
                out[c] += v
        d = dict(out)
    return d.get(0, 0)

def phi_number(word):
    scale = 1
    factors = []
    for n, eps in word:
        if n == 0:
            scale *= 1 + eps
        else:
            factors.append((n, eps))

    d = {(0, 0): 1}
    for n, eps in factors:
        out = collections.defaultdict(int)
        for (a, b), v in d.items():
            for c in cg(a, n):
                out[c, b] += v
            for c in cg(b, n):
                out[a, c] += eps * v
        d = dict(out)
    return scale * d.get((0, 0), 0)

signs = tuple(-1 if (mask >> i) & 1 else 1 for i in range(len(L)))
root_value = phi_number(tuple(zip(L, signs)))
flip_differences = []
for i in range(len(L)):
    for j in range(i + 1, len(L)):
        changed = list(zip(L, signs))
        changed[i] = (changed[i][0], -changed[i][1])
        changed[j] = (changed[j][0], -changed[j][1])
        flip_differences.append(root_value - phi_number(tuple(changed)))

assert root_value == 3532
assert len(flip_differences) == 28
assert min(flip_differences) == -776
assert max(flip_differences) == -4

d = word_closure(L)
coords, index = [], {}

def var(A, B):
    if len(A) == 2:
        A = () if A[0] == A[1] else (0,)
    if len(B) == 2:
        B = () if B[0] == B[1] else (0,)
    if len(A) == 1 or len(B) == 1:
        return None
    key = (A, B) if A <= B else (B, A)
    if key not in index:
        index[key] = len(coords)
        coords.append(key)
    return index[key]

terms = {}
for Q in d:
    terms[Q] = [
        (j, s)
        for s, (A, B) in enumerate(cuts(Q))
        if (j := var(A, B)) is not None
    ]

# Componentwise depth-one relations; rows are deduplicated over ℤ.
relation_rows = [{var((), ()): 1}]
seen_relations = set()

for j, (A, B) in enumerate(coords):
    for side in (0, 1):
        C, D = (A, B) if side == 0 else (B, A)
        if d.get(C, 10**9) >= 1:
            continue
        for i in range(len(C)):
            for h in range(i + 1, len(C)):
                row = collections.Counter({j: 1})
                for R in fuse(C, i, h):
                    z = var(R, D)
                    if z is not None:
                        row[z] -= 1
                row = {u: v for u, v in row.items() if v}
                key = tuple(sorted(row.items()))
                if key and key not in seen_relations:
                    seen_relations.add(key)
                    relation_rows.append(row)

def phi_row(Q, s):
    out = collections.Counter()
    for j, t in terms[Q]:
        out[j] += -1 if (t & s).bit_count() % 2 else 1
    return {j: v for j, v in out.items() if v}

def signed_phi(items):
    scale = 1
    clean = []
    for n, eps in items:
        if n == 0:
            scale *= 1 + eps
        else:
            clean.append((n, eps))
    if scale == 0:
        return {}
    clean = tuple(sorted(clean))
    Q = tuple(n for n, _ in clean)
    s = sum((1 << i) for i, (_, eps) in enumerate(clean) if eps < 0)
    return {j: scale * v for j, v in phi_row(Q, s).items()}

root = np.zeros(len(coords), dtype=np.int64)
for j, v in phi_row(L, mask).items():
    root[j] = int(v)

# Every pair, CG channel, and channel sign; c=0 is handled by signed_phi.
children = []
seen_children = set()
for i in range(len(L)):
    for j in range(i + 1, len(L)):
        for c in cg(L[i], L[j]):
            for eta in (-1, 1):
                items = tuple(sorted(
                    tuple((L[t], signs[t]) for t in range(len(L)) if t not in (i, j))
                    + ((c, eta),)
                ))
                if items in seen_children:
                    continue
                seen_children.add(items)
                v = signed_phi(items)
                if v:
                    children.append({k: -z for k, z in v.items()})

N = len(coords)
E = sparse(relation_rows, N)
P = sparse(children, N)

actual_table = [mult(A) * mult(B) for A, B in coords]
assert sum(int(v) * actual_table[j] for j, v in enumerate(root)) == root_value

Aeq = vstack([E, coo_matrix(root.reshape(1, -1))]).tocsr()
beq = np.zeros(Aeq.shape[0])
beq[-1] = -1

print(
    "MODEL", "subword+one-fusion words", len(d),
    "coordinates", N, "fusion relations", len(relation_rows) - 1,
    "root children", len(children), stamp(), flush=True,
)
print(
    "ROOT", root_value, "double-flips", len(flip_differences),
    "D-range", (min(flip_differences), max(flip_differences)),
    stamp(), flush=True,
)

lp = linprog(
    np.ones(N),
    A_ub=P,
    b_ub=np.zeros(P.shape[0]),
    A_eq=Aeq,
    b_eq=beq,
    bounds=(0, None),
    method="highs",
    options={"time_limit": args.time_limit},
)
print("LP_SUPPORT", lp.status, lp.message, stamp(), flush=True)
assert lp.success

support = np.flatnonzero(lp.x > 1e-8)
active = np.flatnonzero(np.abs(lp.ineqlin.residual) < 1e-7)
V = vstack([Aeq[:, support], P[active, :][:, support]]).tocsr()
right = np.r_[beq, np.zeros(len(active))]
print("EXACTIFY", V.shape, "active-child rows", len(active), stamp(), flush=True)

small = exact_solve(V, right, lp.x[support])
x = [F(0)] * N
for j, value in zip(support, small):
    x[int(j)] = value

ev = mat_vec(E, x)
pv = mat_vec(P, x)
root_pairing = sum(int(v) * x[j] for j, v in enumerate(root))

assert min(x) >= 0
assert all(value == 0 for value in ev)
assert x[index[((), ())]] == 0
assert max(pv, default=0) <= 0
assert root_pairing == -1

print(
    "EXACT_SEPARATOR",
    "support", sum(value != 0 for value in x),
    "denominator_lcm", math.lcm(*(value.denominator for value in x if value)),
    "unit", x[index[((), ())]],
    "max_child_negative", max(pv, default=0),
    "root_pairing", root_pairing,
    stamp(), flush=True,
)
for j, value in enumerate(x):
    if value:
        print("RAY", coords[j], value)
