import argparse
import collections
import datetime
import heapq
import time
from fractions import Fraction as F

argparse.ArgumentParser(
    description="Exact depth-two rank and free-direction counts."
).parse_args()

def stamp():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%H:%M:%S UTC")

def cuts(Q):
    for s in range(1 << len(Q)):
        yield (
            tuple(a for i, a in enumerate(Q) if (s >> i) & 1),
            tuple(a for i, a in enumerate(Q) if not ((s >> i) & 1)),
        )

def cg(a, b):
    return range(abs(a - b), a + b + 1, 2)

def fuse(Q, i, j):
    R = Q[:i] + Q[i + 1:j] + Q[j + 1:]
    return [
        tuple(sorted(R + ((c,) if c else ())))
        for c in cg(Q[i], Q[j])
    ]

def words_at(L, k):
    d = {A: 0 for A, _ in cuts(L)}
    for depth in range(k):
        new = {}
        for Q, level in list(d.items()):
            if level != depth:
                continue
            for i in range(len(Q)):
                for j in range(i + 1, len(Q)):
                    for R in fuse(Q, i, j):
                        if R not in d:
                            new[R] = depth + 1
        d.update(new)
    return d

def rank_domain(name, L, k=2):
    start = time.monotonic()
    d = words_at(L, k)
    coords, index = [], {}

    def var(A, B):
        if len(A) == 2:
            A = () if A[0] == A[1] else (-1,)
        if len(B) == 2:
            B = () if B[0] == B[1] else (-1,)
        if len(A) == 1 or len(B) == 1:
            return None
        key = (A, B) if A <= B else (B, A)
        if key not in index:
            index[key] = len(coords)
            coords.append(key)
        return index[key]

    for Q in d:
        for A, B in cuts(Q):
            var(A, B)

    rows = [{var((), ()): F(1)}]
    for j, (A, B) in enumerate(coords):
        for side in (0, 1):
            C, D = (A, B) if side == 0 else (B, A)
            level = d[tuple(sorted(A + B))]
            if level >= k:
                continue
            for i in range(len(C)):
                for h in range(i + 1, len(C)):
                    row = collections.Counter({j: F(1)})
                    for R in fuse(C, i, h):
                        z = var(R, D)
                        if z is not None:
                            row[z] -= 1
                    rows.append({u: v for u, v in row.items() if v})

    N = len(coords)
    eq = []
    incidence = [set() for _ in range(N)]
    for row in rows:
        if not row:
            continue
        i = len(eq)
        eq.append(row)
        for j in row:
            incidence[j].add(i)

    heap = [(len(row), i) for i, row in enumerate(eq)]
    heapq.heapify(heap)
    rank = 0
    print("BEGIN", name, "N", N, "rows", len(rows), stamp(), flush=True)

    while heap:
        width, i = heapq.heappop(heap)
        row = eq[i]
        if row is None or width != len(row):
            continue
        if not row:
            eq[i] = None
            continue

        pivot = min(
            row,
            key=lambda j: (len(incidence[j]), abs(row[j]) != 1, j),
        )
        a = row[pivot]
        rest = {j: v / a for j, v in row.items() if j != pivot}

        for j in row:
            incidence[j].discard(i)
        eq[i] = None

        for h in list(incidence[pivot]):
            rr = eq[h]
            a = rr.pop(pivot)
            incidence[pivot].discard(h)
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

        rank += 1
        if rank % 250 == 0:
            print("PIVOT", name, rank, "/", N, stamp(), flush=True)

    print(
        "RANK_DONE", name, "N", N, "rank", rank,
        "free", N - rank, "seconds", round(time.monotonic() - start, 2),
        stamp(), flush=True,
    )

for name, labels in (
    ("six", (2, 3, 5, 7, 8, 9)),
    ("sevenA", (1, 3, 5, 7, 9, 9, 12)),
    ("sevenB", (1, 2, 5, 7, 9, 11, 13)),
):
    rank_domain(name, labels)
