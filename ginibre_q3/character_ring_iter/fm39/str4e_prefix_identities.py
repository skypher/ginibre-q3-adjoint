import argparse
import collections
import datetime
import itertools
import math
from fractions import Fraction as Q
from functools import lru_cache

parser = argparse.ArgumentParser(description="Exact HPP fusion-identity certificate checker.")
parser.parse_args()

def stamp():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%H:%M:%S UTC")

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

@lru_cache(None)
def tab(word):
    d = {(0, 0): 1}
    for n, eps in word:
        out = collections.defaultdict(int)
        for (r, s), value in d.items():
            for c in cg(r, n):
                out[c, s] += value
            for c in cg(s, n):
                out[r, c] += eps * value
        d = {k: v for k, v in out.items() if v}
    return tuple(sorted(d.items()))

def pi(A, B, T):
    fa, fb = dict(tab(tuple(sorted(A)))), dict(tab(tuple(sorted(B))))
    return sum(v * fb.get(k, 0) for k, v in fa.items() if k[0] + k[1] <= T)

def even_signs(n):
    return [s for s in itertools.product((1, -1), repeat=n) if math.prod(s) == 1]

def fused_child(L, cut, signs, i, j, c, side):
    A = [(L[k], signs[k]) for k in range(cut) if k not in (i, j)]
    B = [(L[k], signs[k]) for k in range(cut, len(L)) if k not in (i, j)]
    (A if side == "A" else B).append((c, signs[i] * signs[j]))
    return A, B

def moved_child(L, cut, signs, i):
    assert i < cut
    A = [(L[k], signs[k]) for k in range(cut) if k != i]
    B = [(L[k], signs[k]) for k in range(cut, len(L))] + [(L[i], signs[i])]
    return A, B

def verify_identity(name, L, cut, T, support):
    ss = even_signs(len(L))
    for signs in ss:
        lhs = pi([(L[i], signs[i]) for i in range(cut)],
                 [(L[i], signs[i]) for i in range(cut, len(L))], T)
        rhs = 0
        for weight, op in support:
            if op[0] == "F":
                _, i, j, c, side = op
                A, B = fused_child(L, cut, signs, i, j, c, side)
            else:
                _, i = op
                A, B = moved_child(L, cut, signs, i)
            rhs += weight * pi(A, B, T)
        assert lhs == rhs, (name, T, signs, lhs, rhs)
    print("IDENTITY", name, "T", T, "even_signings", len(ss),
          "support", [(str(w), op) for w, op in support],
          "PASS", stamp(), flush=True)

print("BEGIN", stamp(), flush=True)

small = (1, 1, 1, 1, 2)
for T in (0, 1):
    verify_identity("small-singlet", small, 2, T,
                    [(Q(1), ("F", 0, 1, 0, "A"))])
for T in range(2, 6):
    verify_identity("small-terminal-move", small, 2, T,
                    [(Q(1), ("M", 0))])

hard = (1, 4, 5, 6, 6, 8, 9, 15)
hard_support = {
    2: [(Q(5, 3), ("F", 0, 4, 5, "B")),
        (Q(2), ("F", 1, 2, 5, "A"))],
    4: [(Q(1), ("F", 0, 2, 6, "A")),
        (Q(13, 18), ("F", 0, 4, 5, "B")),
        (Q(1), ("F", 1, 2, 5, "A"))],
    6: [(Q(1), ("F", 0, 2, 4, "A")),
        (Q(1), ("F", 0, 2, 6, "A"))],
    8: [(Q(1), ("F", 0, 2, 4, "A")),
        (Q(1), ("F", 0, 2, 6, "A"))],
}
hard_signs = (1, 1, -1, 1, 1, 1, -1, 1)
for T, support in hard_support.items():
    verify_identity("prior-hard-list", hard, 4, T, support)
    actual = pi([(hard[i], hard_signs[i]) for i in range(4)],
                [(hard[i], hard_signs[i]) for i in range(4, 8)], T)
    print("ACTUAL_SIGNED_PARENT", "hard", "T", T, "Pi", actual, flush=True)

census = (1, 3, 3, 4, 5, 5, 6, 7)
census_support = {
    2: [(Q(9, 943), ("F", 0, 6, 7, "B")),
        (Q(486, 943), ("F", 0, 7, 8, "B")),
        (Q(675, 943), ("F", 1, 2, 4, "A")),
        (Q(232, 943), ("F", 6, 7, 1, "B"))],
    4: [(Q(1), ("F", 1, 2, 0, "A")),
        (Q(1), ("F", 1, 2, 4, "A")),
        (Q(1), ("F", 1, 3, 1, "A")),
        (Q(1), ("F", 2, 3, 1, "A"))],
}
census_signs = (-1,) * 8
for T, support in census_support.items():
    verify_identity("census-tight", census, 4, T, support)
    actual = pi([(census[i], census_signs[i]) for i in range(4)],
                [(census[i], census_signs[i]) for i in range(4, 8)], T)
    print("ACTUAL_SIGNED_PARENT", "census", "T", T, "Pi", actual, flush=True)

L = small
ss = even_signs(len(L))
y = [Q(-1), Q(-1, 2), Q(-1, 2), Q(1), Q(-1, 2), Q(1), Q(1), Q(-1, 2),
     Q(1), Q(1), Q(1), Q(1), Q(1), Q(1), Q(1), Q(-1)]
assert len(ss) == len(y) == 16

parent = [pi([(L[i], s[i]) for i in range(2)],
             [(L[i], s[i]) for i in range(2, 5)], 2) for s in ss]
columns, tags = [], []
for i in range(len(L)):
    for j in range(i+1, len(L)):
        same_side = (i < 2) == (j < 2)
        for c in cg(L[i], L[j]):
            placements = (("A",) if same_side and i < 2
                          else (("B",) if same_side else ("A", "B")))
            for side in placements:
                vals = []
                for s in ss:
                    A, B = fused_child(L, 2, s, i, j, c, side)
                    vals.append(pi(A, B, 2))
                if any(vals):
                    columns.append(vals)
                    tags.append((i, j, c, side))
columns.append([1] * len(ss))
tags.append(("constant",))

dots = [sum(y[r] * col[r] for r in range(len(ss))) for col in columns]
score = sum(y[r] * parent[r] for r in range(len(ss)))
assert score == -12
assert min(dots) == 0 and max(dots) == 16 and all(v >= 0 for v in dots)
assert dots[-1] == 6

print("SEPARATOR", "rows", len(ss), "fusion-plus-constant-columns", len(columns),
      "y_dot_target", score, "column_dot_range", (min(dots), max(dots)),
      "constant_dot", dots[-1], "zero_columns", sum(v == 0 for v in dots),
      "PASS", stamp(), flush=True)
print("SEPARATOR_Y", [("".join("+" if e > 0 else "-" for e in s), str(v))
                       for s, v in zip(ss, y)], flush=True)
print("EXACT HPP LP CERTIFICATES PASS", stamp(), flush=True)
