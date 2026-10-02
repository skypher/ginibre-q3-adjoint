import gzip
import argparse
import collections
import contextlib
import datetime
import itertools
import math
import re
import signal
import time
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack


# W <= 40 census logs of sec166_flip_descent_census.cpp, stored gzipped next to this file.
NOFLIP_LOG = Path(__file__).resolve().parent / "sec166_census_w40_noflip.log.gz"
TPFAIL_LOG = Path(__file__).resolve().parent / "sec166_census_w40_tpfail.log.gz"


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


@lru_cache(None)
def mult(Q):
    r = {0: 1}
    for a in Q:
        s = collections.defaultdict(int)
        for b, x in r.items():
            for c in range(abs(a - b), a + b + 1, 2):
                s[c] += x
        r = dict(s)
    return r.get(0, 0)


def cuts(Q):
    for mask in range(1 << len(Q)):
        yield (
            tuple(a for i, a in enumerate(Q) if mask >> i & 1),
            tuple(a for i, a in enumerate(Q) if not (mask >> i & 1)),
        )


def fuse(Q, i, j):
    R = Q[:i] + Q[i + 1:j] + Q[j + 1:]
    out = []
    for c in range(abs(Q[i] - Q[j]), Q[i] + Q[j] + 1, 2):
        out.append(tuple(sorted(R + ((c,) if c else ()))))
    return out


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


def character_values(Q):
    values = [mult(A) * mult(B) for A, B in cuts(Q)]
    for k in range(len(Q)):
        for mask in range(1 << len(Q)):
            if mask >> k & 1:
                continue
            a, b = values[mask], values[mask | (1 << k)]
            values[mask] = a + b
            values[mask | (1 << k)] = a - b
    return values


def child(L, mask):
    n = len(L)
    i, j = max(
        (
            (i, j)
            for i in range(n - 1)
            for j in range(i + 1, n - 1)
            if (L[i] - L[j]) % 2 == 0
        ),
        key=lambda ij: (L[ij[0]] + L[ij[1]], max(L[ij[0]], L[ij[1]])),
    )
    keep = [h for h in range(n) if h not in (i, j)]
    C = tuple(L[h] for h in keep)
    cm = sum(((mask >> h) & 1) << q for q, h in enumerate(keep))
    if ((mask >> i) ^ (mask >> j)) & 1:
        cm ^= 1 << (len(C) - 1)
    return C, cm, (i, j)


def sparse(rows, n):
    rr, cc, vv = [], [], []
    for i, row in enumerate(rows):
        for j, v in row.items():
            if v:
                rr.append(i)
                cc.append(j)
                vv.append(v)
    return coo_matrix(
        (vv, (rr, cc)), shape=(len(rows), n), dtype=np.int64
    ).tocsr()


def mv(M, x):
    out = []
    for i in range(M.shape[0]):
        lo, hi = M.indptr[i], M.indptr[i + 1]
        out.append(sum(
            int(v) * x[j]
            for j, v in zip(M.indices[lo:hi], M.data[lo:hi])
        ))
    return out


def make_model(L, k, mode, pairs=True):
    d = words_at(L, k)
    coords, index = [], {}

    def var(A, B):
        if pairs and len(A) == 2:
            A = () if A[0] == A[1] else (-1,)
        if pairs and len(B) == 2:
            B = () if B[0] == B[1] else (-1,)
        if len(A) == 1 or len(B) == 1:
            return None
        key = (A, B) if A <= B else (B, A)
        if key not in index:
            index[key] = len(coords)
            coords.append(key)
        return index[key]

    terms = {}
    for Q in d:
        row = []
        for s, (A, B) in enumerate(cuts(Q)):
            j = var(A, B)
            if j is not None:
                row.append((j, s))
        terms[Q] = row

    E = [{var((), ()): 1}]
    et = [("unit",)]
    for j, (A, B) in enumerate(coords):
        for side in (0, 1):
            C, D = (A, B) if side == 0 else (B, A)
            level = d[C] if mode == "component" else d[tuple(sorted(A + B))]
            if level >= k:
                continue
            for i in range(len(C)):
                for h in range(i + 1, len(C)):
                    row = collections.Counter({j: 1})
                    for R in fuse(C, i, h):
                        z = var(R, D)
                        if z is not None:
                            row[z] -= 1
                    E.append(row)
                    et.append((A, B, side, i, h))

    def phi(Q, mask):
        row = collections.Counter()
        for j, s in terms[Q]:
            row[j] += -1 if (s & mask).bit_count() % 2 else 1
        return row

    P, pt = [], []
    for Q in d:
        if len(Q) >= len(L):
            continue
        for mask in range(1 << len(Q)):
            if mask.bit_count() % 2:
                continue
            P.append({j: -v for j, v in phi(Q, mask).items()})
            pt.append((Q, mask))

    N = len(coords)
    return NS(
        L=L, k=k, mode=mode, pairs=pairs, d=d, coords=coords,
        index=index, phi=phi, N=N, E=sparse(E, N), P=sparse(P, N),
        et=et, pt=pt,
    )


def profile(M, mask):
    L = M.L
    C, cm, _ = child(L, mask)
    p = M.phi(L, mask)
    obj = p.copy()
    for j, v in M.phi(C, cm).items():
        obj[j] -= v
    c = np.zeros(M.N, dtype=np.int64)
    for j, v in obj.items():
        assert v % 2 == 0
        c[j] = v // 2

    rows, tags = [], []
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            row = p.copy()
            for z, v in M.phi(L, mask ^ (1 << i) ^ (1 << j)).items():
                row[z] -= v
            rows.append(row)
            tags.append((i, j))
    return c, sparse(rows, M.N), tags


def exact_solve(mat, right, approx):
    mat = mat.tocsr()
    equations, rhs = [], []
    incidence = [set() for _ in approx]
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
            assert not b
            continue
        h = len(equations)
        equations.append(row)
        rhs.append(b)
        for j in row:
            incidence[j].add(h)

    heap = [(len(row), i) for i, row in enumerate(equations)]
    heapq = __import__("heapq")
    heapq.heapify(heap)
    saved = []
    while heap:
        width, i = heapq.heappop(heap)
        row = equations[i]
        if row is None or width != len(row):
            continue
        if not row:
            assert not rhs[i]
            equations[i] = None
            continue
        pivot = min(row, key=lambda j: (len(incidence[j]), abs(row[j]) != 1, j))
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
                z = rr.get(j, 0) - a * v
                if z:
                    if j not in rr:
                        incidence[j].add(h)
                    rr[j] = z
                elif j in rr:
                    del rr[j]
                    incidence[j].discard(h)
            heapq.heappush(heap, (len(rr), h))

    x = [F(float(v)).limit_denominator(10**6) for v in approx]
    for pivot, rest, b in reversed(saved):
        x[pivot] = b - sum(v * x[j] for j, v in rest.items())
    return x


def dual(M, mask, allow_flips, lp_seconds):
    c, flips, tags = profile(M, mask)
    U = vstack([M.P, flips]).tocsr() if allow_flips else M.P
    npterms = M.P.shape[0]
    b = np.zeros(M.E.shape[0])
    b[0] = 1
    r = linprog(
        c,
        A_ub=U if U.shape[0] else None,
        b_ub=np.zeros(U.shape[0]) if U.shape[0] else None,
        A_eq=M.E,
        b_eq=b,
        bounds=(0, None),
        method="highs",
        options={"time_limit": lp_seconds},
    )
    if not r.success:
        return {"status": r.status, "message": r.message, "c": c, "U": U, "r": r}

    iy = np.flatnonzero(abs(r.eqlin.marginals) > 1e-8)
    iz = np.flatnonzero(abs(r.ineqlin.marginals) > 1e-8)
    V = vstack([M.E[iy, :], U[iz, :]]).tocsc()
    active = np.flatnonzero(abs(r.lower.marginals) < 1e-7)
    a = exact_solve(
        V[:, active].T,
        c[active],
        np.r_[r.eqlin.marginals[iy], r.ineqlin.marginals[iz]],
    )
    den = math.lcm(*(v.denominator for v in a))
    z = [int(v) * den for v in c]
    ai = [int(v * den) for v in a]
    for j in range(M.N):
        for h, v in zip(
            V.indices[V.indptr[j]:V.indptr[j + 1]],
            V.data[V.indptr[j]:V.indptr[j + 1]],
        ):
            z[j] -= int(v) * ai[h]

    q = a[list(iy).index(0)] if 0 in iy else F(0)
    assert min(z, default=0) >= 0
    assert all(v <= 0 for v in a[len(iy):])

    alpha, shorter = [], []
    for h, rowidx in enumerate(iz):
        val = a[len(iy) + h]
        if rowidx < npterms and val:
            shorter.append((M.pt[rowidx][0], -val))
        elif allow_flips and rowidx >= npterms and val:
            alpha.append((tags[rowidx - npterms], -4 * val))

    return {
        "status": 0, "q": q, "alpha": tuple(alpha),
        "shorter": tuple(shorter), "c": c, "U": U, "r": r,
        "dual_support": sum(v > 0 for v in z),
        "fusion_support": len(iy) - (0 in iy),
    }


def negative_point(M, mask, c, U, r):
    b = [0] * M.E.shape[0]
    b[0] = 1
    ss = np.flatnonzero(r.x > 1e-8)
    tt = np.flatnonzero(abs(r.ineqlin.residual) < 1e-7)
    A = vstack([M.E[:, ss], U[tt, :][:, ss]]).tocsr()
    xx = exact_solve(
        A, np.r_[np.array(b, dtype=float), np.zeros(len(tt))], r.x[ss]
    )
    x = [F(0)] * M.N
    for j, v in zip(ss, xx):
        x[j] = v
    assert min(x, default=F(0)) >= 0
    assert mv(M.E, x) == b
    assert max(mv(U, x), default=0) <= 0
    value = sum(int(v) * xx for v, xx in zip(c, x))
    assert value < 0
    return value


def negative_ray(M, c, U, lp_seconds):
    Eaug = vstack([M.E, coo_matrix(c.reshape(1, -1))]).tocsr()
    rhs = np.zeros(Eaug.shape[0])
    rhs[-1] = -1
    r = linprog(
        np.ones(M.N),
        A_ub=U,
        b_ub=np.zeros(U.shape[0]),
        A_eq=Eaug,
        b_eq=rhs,
        bounds=(0, None),
        method="highs",
        options={"time_limit": lp_seconds},
    )
    assert r.success, r.message
    ss = np.flatnonzero(r.x > 1e-8)
    tt = np.flatnonzero(abs(r.ineqlin.residual) < 1e-7)
    A = vstack([Eaug[:, ss], U[tt, :][:, ss]]).tocsr()
    xx = exact_solve(
        A, np.r_[rhs, np.zeros(len(tt))], r.x[ss]
    )
    x = [F(0)] * M.N
    for j, v in zip(ss, xx):
        x[j] = v
    assert min(x, default=F(0)) >= 0
    assert mv(M.E, x) == [0] * M.E.shape[0]
    assert max(mv(U, x), default=0) <= 0
    assert sum(int(v) * w for v, w in zip(c, x)) == -1
    return -1


def read_groups(path, tag):
    groups = collections.defaultdict(set)
    rows = 0
    text = gzip.decompress(path.read_bytes()).decode() if path.suffix == ".gz" else path.read_text()
    for line in text.splitlines():
        if not line.startswith(tag + " "):
            continue
        stop = "phi=" if tag == "NOFLIP" else "flip="
        m = re.search(r" p=([+-]?\d+) B=(.*?) " + re.escape(stop), line)
        assert m, line
        p = int(m.group(1))
        B = list(map(int, m.group(2).split()))
        assert [abs(x) for x in B] == sorted(map(abs, B))
        key = (tuple(sorted(map(abs, B))), abs(p))
        mask = sum(1 << i for i, x in enumerate(B + [p]) if x < 0)
        groups[key].add(mask)
        rows += 1
    return groups, rows


@lru_cache(None)
def phi_scalar(Q, mask):
    counts = collections.Counter(Q)
    values = sorted(counts)
    neg = {
        a: sum(1 for i, x in enumerate(Q) if x == a and (mask >> i) & 1)
        for a in values
    }
    choices = []
    for a in values:
        n, p = neg[a], counts[a] - neg[a]
        choices.append([
            (u + v, (-1) ** u * math.comb(n, u) * math.comb(p, v))
            for u in range(n + 1)
            for v in range(p + 1)
        ])
    total = 0
    for pick in itertools.product(*choices):
        A = tuple(a for a, (n, _) in zip(values, pick) for _ in range(n))
        weight = math.prod(w for _, w in pick)
        if not weight:
            continue
        rem = collections.Counter(Q) - collections.Counter(A)
        B = tuple(a for a in sorted(rem) for _ in range(rem[a]))
        total += weight * mult(A) * mult(B)
    return total


def representative_pairs(L):
    positions = collections.defaultdict(list)
    for i, a in enumerate(L):
        positions[a].append(i)
    values = sorted(positions)
    pairs = []
    for k, a in enumerate(values):
        if len(positions[a]) >= 2:
            pairs.append((positions[a][0], positions[a][1]))
        for b in values[k + 1:]:
            pairs.append((positions[a][0], positions[b][0]))
    return pairs


def is_noflip(L, mask):
    value = phi_scalar(L, mask)
    for i, j in representative_pairs(L):
        if value - phi_scalar(L, mask ^ (1 << i) ^ (1 << j)) >= 0:
            return False
    return True


def timestamped_alarm(signum, frame):
    raise TimeoutError("model build cap")


def build_capped(L, k, mode, seconds):
    signal.signal(signal.SIGALRM, timestamped_alarm)
    signal.alarm(seconds)
    try:
        return make_model(L, k, mode, pairs=True)
    finally:
        signal.alarm(0)


def verify_example(lp_seconds):
    L = (2, 3, 5, 7, 8, 9)
    masks = (10, 17, 36, 63)
    M1 = make_model(L, 1, "component", pairs=True)
    M2 = make_model(L, 2, "total", pairs=True)
    for mask in masks:
        out = dual(M1, mask, True, lp_seconds)
        assert out["status"] == 0 and out["q"] == -1
        assert not out["alpha"]
        assert negative_point(M1, mask, out["c"], out["U"], out["r"]) == -1
        C, cm, _ = child(L, mask)
        parent = phi_scalar(L, mask)
        child_value = phi_scalar(C, cm)
        assert parent // 2 == 64 and child_value // 2 == 3
        assert (parent - child_value) // 2 == 61
        ds = []
        for i in range(len(L)):
            for j in range(i + 1, len(L)):
                ds.append((
                    parent - phi_scalar(L, mask ^ (1 << i) ^ (1 << j))
                ) // 4)
        assert max(ds) == -1
        out2 = dual(M2, mask, True, lp_seconds)
        assert out2["status"] == 0 and out2["q"] == 61
        assert not out2["alpha"] and not out2["shorter"]
    print("EXAMPLE exact: depth-one -1; actual drop 61; depth-two bound 61")


def check_tpfail(path):
    groups, rows = read_groups(path, "TPFAIL")
    tested = 0
    noflip = []
    for ix, (key, _) in enumerate(sorted(groups.items()), 1):
        B, p = key
        classes = sorted(set(B))
        counts = collections.Counter(B)
        L = B + (p,)
        for bits in itertools.product((0, 1), repeat=len(classes)):
            negclass = dict(zip(classes, bits))
            pneg = sum(counts[a] for a in classes if negclass[a]) % 2
            if p in negclass and negclass[p] != pneg:
                continue
            mask = sum(
                1 << i for i, a in enumerate(B) if negclass[a]
            ) | (pneg << len(B))
            tested += 1
            if is_noflip(L, mask):
                noflip.append((key, mask))
        if ix % 10 == 0:
            print(stamp(), "TPFAIL progress", ix, "/", len(groups), flush=True)
    assert rows == 126 and len(groups) == 63
    assert tested == 976 and not noflip
    print("TPFAIL exact:", rows, "rows;", len(groups), "multisets;",
          tested, "pair-free signings; no no-flip signing")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--noflip-log", type=Path, default=NOFLIP_LOG)
    ap.add_argument("--tp-fail-log", type=Path, default=TPFAIL_LOG)
    ap.add_argument("--l7-prefix", type=int, default=323)
    ap.add_argument("--profile-seconds", type=int, default=30)
    ap.add_argument("--lp-seconds", type=float, default=5)
    args = ap.parse_args()

    groups, rows = read_groups(args.noflip_log, "NOFLIP")
    assert rows == 5430 and len(groups) == 2014
    l6 = sorted(
        [(k, sorted(v)) for k, v in groups.items() if len(k[0]) + 1 == 6],
        key=lambda z: (max(z[0][0] + (z[0][1],)), z[0]),
    )
    l7 = sorted(
        [(k, sorted(v)) for k, v in groups.items() if len(k[0]) + 1 == 7],
        key=lambda z: (max(z[0][0] + (z[0][1],)), z[0]),
    )
    selected = l6 + l7[:args.l7_prefix]
    assert len(l6) == 153 and len(selected) == 153 + args.l7_prefix
    assert sum(len(masks) for _, masks in selected) == 1476

    stats = collections.Counter()
    shorter_forms = collections.Counter()
    exact_separation = collections.Counter()
    t0 = time.monotonic()
    print(stamp(), "START", len(selected), "domains",
          sum(len(m) for _, m in selected), "signings", flush=True)

    for ix, (key, masks) in enumerate(selected, 1):
        L = key[0] + (key[1],)
        try:
            M1 = build_capped(L, 1, "component", args.profile_seconds)
        except TimeoutError:
            stats["k1 model capped"] += 1
            continue

        k1_ok = True
        for mask in masks:
            out = dual(M1, mask, True, args.lp_seconds)
            if out["status"] == 0 and out["q"] >= 0:
                stats["k1 signing certified"] += 1
                for Q, _ in out["shorter"]:
                    shorter_forms[Q] += 1
                if out["alpha"]:
                    stats["active flip duals"] += 1
                    zero = dual(M1, mask, False, args.lp_seconds)
                    if zero["status"] == 0 and zero["q"] >= 0:
                        stats["active flip but zero support exists"] += 1
                    else:
                        stats["flip support unresolved"] += 1
                continue

            if out["status"] == 0:
                val = negative_point(M1, mask, out["c"], out["U"], out["r"])
                assert val < 0
                exact_separation["negative point"] += 1
            elif out["status"] == 3:
                negative_ray(M1, out["c"], out["U"], args.lp_seconds)
                exact_separation["negative ray"] += 1
            else:
                stats["k1 LP capped"] += 1
                k1_ok = False
                break
            stats["k1 exact separator"] += 1
            k1_ok = False
            break

        if k1_ok:
            stats["depth one domain"] += 1
        else:
            try:
                M2 = build_capped(L, 2, "total", args.profile_seconds)
            except TimeoutError:
                stats["k2 model capped"] += 1
                continue
            k2_ok = True
            for mask in masks:
                out = dual(M2, mask, True, args.lp_seconds)
                if out["status"] != 0 or out["q"] < 0:
                    k2_ok = False
                    break
                stats["k2 flip terms"] += len(out["alpha"])
                stats["k2 shorter terms"] += len(out["shorter"])
            if k2_ok:
                stats["depth two domain"] += 1
            else:
                stats["k2 unresolved"] += 1

        if ix % 25 == 0 or ix == len(selected):
            print(stamp(), "PROGRESS", ix, "/", len(selected),
                  dict(stats), "elapsed", round(time.monotonic() - t0, 1),
                  flush=True)

    print(stamp(), "NOFLIP SUMMARY", dict(stats))
    print("DEPTH-ONE SHORTER-FORM LENGTHS",
          dict(collections.Counter({len(Q): n for Q, n in shorter_forms.items()})))
    print("MOST REUSED SHORTER FORMS", shorter_forms.most_common(12))
    print("EXACT SEPARATION TYPES", dict(exact_separation))
    if args.l7_prefix == 323:
        assert stats["depth one domain"] == 425
        assert stats["depth two domain"] == 51
        assert stats["k1 model capped"] == 0
        assert stats["k2 model capped"] == 0
        assert stats["k2 unresolved"] == 0

    verify_example(args.lp_seconds)
    check_tpfail(args.tp_fail_log)

    if args.l7_prefix < len(l7):
        print("FIRST UNTESTED L7", l7[args.l7_prefix])


if __name__ == "__main__":
    main()
