import argparse
from functools import lru_cache
from itertools import combinations_with_replacement, combinations, product
from fractions import Fraction
from math import comb

parser = argparse.ArgumentParser(
    description="Fresh exact FM3 split and polynomial screen."
)
parser.parse_args()

@lru_cache(None)
def m(ns):
    state = {0: 1}
    for n in ns:
        nxt = {}
        for a, v in state.items():
            for b in range(abs(a - n), a + n + 1, 2):
                nxt[b] = nxt.get(b, 0) + v
        state = nxt
    return state.get(0, 0)

def records(ns):
    L = len(ns)
    out = []
    for mask in range(1 << (L - 1)):
        A = tuple(ns[i] for i in range(L) if mask >> i & 1)
        B = tuple(ns[i] for i in range(L) if not (mask >> i & 1))
        out.append((mask, m(A) * m(B), max(A) if A else 0))
    return out

def evaluate(ns, minus_mask, rec=None):
    if rec is None:
        rec = records(ns)
    mfull = m(tuple(ns))
    F = mfull
    P = N = 0
    beta = Fraction(0)
    for mask, f, mx in rec:
        if not mask:
            continue
        if (mask & minus_mask).bit_count() & 1:
            F -= f
            N += f
            if f:
                beta += Fraction(1, mx + 1)
        else:
            F += f
            P += f
    return F, mfull, P, N, beta

def budget(ns, minus_mask):
    ans = Fraction(0)
    for j in range(len(ns) - 1):
        k = j + 1
        if k == 1:
            c = int(bool(minus_mask & 1))
        else:
            prior = minus_mask & ((1 << (k - 1)) - 1)
            c = (2 ** (k - 2) if prior else
                 (2 ** (k - 1) if (minus_mask >> j) & 1 else 0))
        direct = sum(
            1 for s in range(1 << (j + 1))
            if s >> j & 1 and (s & minus_mask).bit_count() & 1
        )
        assert c == direct
        ans += Fraction(c, ns[j] + 1)
    return ans

def all_f(ns):
    L = len(ns)
    return [
        m(tuple(ns[j] for j in range(L) if s >> j & 1)) *
        m(tuple(ns[j] for j in range(L) if not (s >> j & 1)))
        for s in range(1 << L)
    ]

def toggle(ns, minus_mask, i, fs):
    total = 0
    first = None
    L = len(ns)
    for s in range(1 << L):
        if s >> i & 1:
            continue
        d = (-1 if (s & minus_mask).bit_count() & 1 else 1) * (
            fs[s] - fs[s | (1 << i)]
        )
        total += d
        if d < 0 and first is None:
            A = tuple(ns[j] for j in range(L) if s >> j & 1)
            first = (s, A, fs[s], fs[s | (1 << i)],
                     (s & minus_mask).bit_count() & 1, d)
    return total, first

def cat(n):
    return comb(2 * n, n) // (n + 1)

def U(n):
    prev, cur = {0: 1}, {1: 1}
    if n == 0:
        return prev
    if n == 1:
        return cur
    for _ in range(1, n):
        nxt = {d + 1: c for d, c in cur.items()}
        for d, c in prev.items():
            nxt[d] = nxt.get(d, 0) - c
        prev, cur = cur, {d: c for d, c in nxt.items() if c}
    return cur

def poly_expect(ns, eps):
    state = {(0, 0): 1}
    for n, e in zip(ns, eps):
        fac = {}
        for d, c in U(n).items():
            fac[(d, 0)] = fac.get((d, 0), 0) + c
            fac[(0, d)] = fac.get((0, d), 0) + e * c
        nxt = {}
        for (a, b), v in state.items():
            for (c, d), w in fac.items():
                nxt[(a + c, b + d)] = nxt.get((a + c, b + d), 0) + v * w
        state = {k: v for k, v in nxt.items() if v}
    return sum(
        v * cat(a // 2) * cat(b // 2)
        for (a, b), v in state.items()
        if a % 2 == 0 and b % 2 == 0
    )

# Repeated labels 3..7, one sign per value, even minus parity.
repeat_rows = []
total = nonpos = beta_over = 0
minratio = None
maxbeta = (Fraction(-1), None)
min_surplus = (None, None)
toggle_cases = toggle_good = 0

for L in range(2, 8):
    row = certified = 0
    for ns in combinations_with_replacement(range(3, 8), L):
        vals = tuple(sorted(set(ns)))
        rec = records(ns)
        fs = all_f(ns)
        for signs in product((1, -1), repeat=len(vals)):
            sign_by_value = dict(zip(vals, signs))
            eps = tuple(sign_by_value[n] for n in ns)
            minus_mask = sum(
                1 << i for i, e in enumerate(eps) if e < 0
            )
            if minus_mask.bit_count() % 2:
                continue

            F, mfull, P, N, beta = evaluate(ns, minus_mask, rec)
            if not mfull:
                continue

            total += 1
            row += 1
            nonpos += F <= 0
            beta_over += beta > 1
            certified += budget(ns, minus_mask) <= 1

            if beta > maxbeta[0]:
                maxbeta = (beta, (ns, eps, F, mfull, P, N))
            ratio = Fraction(F, mfull)
            if minratio is None or ratio < minratio[0]:
                minratio = (ratio, (ns, eps, F, mfull, beta, P, N))
            surplus = Fraction(P - N, mfull)
            if min_surplus[0] is None or surplus < min_surplus[0]:
                min_surplus = (
                    surplus, (ns, eps, F, mfull, beta, P, N)
                )

            if minus_mask:
                toggle_cases += 1
                has_good_toggle = False
                for i, e in enumerate(eps):
                    if e > 0:
                        continue
                    pair_sum, _ = toggle(ns, minus_mask, i, fs)
                    assert pair_sum == 2 * F
                    good = True
                    for s in range(1 << L):
                        if s >> i & 1:
                            continue
                        d = (-1 if (s & minus_mask).bit_count() & 1 else 1) * (
                            fs[s] - fs[s | (1 << i)]
                        )
                        if d < 0:
                            good = False
                            break
                    has_good_toggle |= good
                toggle_good += has_good_toggle

    repeat_rows.append((L, row, certified))

assert (total, nonpos, beta_over) == (1900, 0, 523)
assert minratio[0] == Fraction(42, 43)
assert maxbeta[0] == 4
assert min_surplus[0] == Fraction(-1, 43)
assert (toggle_cases, toggle_good) == (1511, 908)
print("repeat rows:", repeat_rows)
print("repeat totals:", total, nonpos, beta_over,
      "min F/m", minratio, "max beta", maxbeta,
      "min (P-N)/m", min_surplus)
print("toggle cases / a termwise nonnegative toggle:",
      toggle_cases, toggle_good)

# Distinct labels 3..14, even minus parity.
distinct_rows = []
distinct_total = distinct_nonpos = distinct_beta_over = 0
distinct_min = None

for L in range(2, 9):
    cases = 0
    minrow = None
    maxbeta = Fraction(-1)
    maxcase = None
    for ns in combinations(range(3, 15), L):
        mfull = m(ns)
        if not mfull:
            continue
        rec = records(ns)
        for eps in product((1, -1), repeat=L):
            minus_mask = sum(
                1 << i for i, e in enumerate(eps) if e < 0
            )
            if minus_mask.bit_count() % 2:
                continue

            F, _, P, N, beta = evaluate(ns, minus_mask, rec)
            cases += 1
            distinct_total += 1
            distinct_nonpos += F <= 0
            distinct_beta_over += beta > 1
            ratio = Fraction(F, mfull)
            if minrow is None or ratio < minrow:
                minrow = ratio
            if distinct_min is None or ratio < distinct_min[0]:
                distinct_min = (ratio, (ns, eps, F, mfull, beta))
            if beta > maxbeta:
                maxbeta = beta
                maxcase = (ns, eps, F, mfull)
    distinct_rows.append((L, cases, minrow, maxbeta, maxcase))

assert [r[1] for r in distinct_rows] == [
    0, 360, 2032, 6336, 14464, 25344, 32640
]
assert (distinct_total, distinct_nonpos, distinct_beta_over) == (
    81176, 0, 33080
)
assert distinct_min[0] == Fraction(608, 643)
print("distinct rows:", distinct_rows)
print("distinct totals:", distinct_total, distinct_nonpos,
      distinct_beta_over, "min", distinct_min)

# Exact witnesses; polynomial moments independently evaluate the expectation.
w1 = ((4,) * 6 + (6,), (-1,) * 6 + (1,))
minus1 = sum(1 << i for i, e in enumerate(w1[1]) if e < 0)
q1 = evaluate(w1[0], minus1)
assert q1 == (545, 295, 330, 80, Fraction(4))
assert poly_expect(*w1) == 1090
print("beta witness:", w1, q1, "direct expectation", poly_expect(*w1))

w2 = ((3, 4, 5, 5, 6, 7), (1, -1, 1, 1, -1, 1))
minus2 = sum(1 << i for i, e in enumerate(w2[1]) if e < 0)
q2 = evaluate(w2[0], minus2)
assert q2 == (84, 86, 4, 6, Fraction(13, 14))
assert poly_expect(*w2) == 168
print("mass witness:", w2, q2, "direct expectation", poly_expect(*w2))

w3 = ((3, 3, 4, 6, 7, 7), (1, 1, 1, 1, -1, -1))
minus3 = sum(1 << i for i, e in enumerate(w3[1]) if e < 0)
q3 = evaluate(w3[0], minus3)
fs3 = all_f(w3[0])
assert q3 == (77, 71, 10, 4, Fraction(1, 2))
assert poly_expect(*w3) == 154
for i, e in enumerate(w3[1]):
    if e < 0:
        pair_sum, defect = toggle(w3[0], minus3, i, fs3)
        assert pair_sum == 154
        assert defect == (5, (3, 4), 0, 1, 0, -1)
        print("toggle witness i / pair sum / defect:",
              i, pair_sum, defect)
print("toggle witness:", w3, q3[:2],
      "direct expectation", poly_expect(*w3))

print("small-length upper budgets:",
      [Fraction(2 ** (L - 2), 4) for L in range(2, 5)])
print("PASS")