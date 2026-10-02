import argparse
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction as Q
from itertools import combinations_with_replacement, product
from math import comb

ap = argparse.ArgumentParser(description="Exact braided-prefix path verifier; no files written.")
ap.parse_args()

def stamp(message):
    print(datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
          message, flush=True)

def padd(a, b, scale=1):
    out = defaultdict(int, a)
    for k, v in b.items():
        out[k] += scale * v
    return {k: v for k, v in out.items() if v}

def vadd(a, b, scale=1):
    out = {w: dict(p) for w, p in a.items()}
    for w, p in b.items():
        out[w] = padd(out.get(w, {}), p, scale)
    return {w: p for w, p in out.items() if p}

def pshift(p, degree):
    return {k + degree: v for k, v in p.items()}

def pmul(a, b):
    out = defaultdict(int)
    for i, x in a.items():
        for j, y in b.items():
            out[i + j] += x * y
    return {k: v for k, v in out.items() if v}

def x0(color, vec):
    out = {(color,) + w: p for w, p in vec.items()}
    for w, p in vec.items():
        try:
            j = w.index(color)
        except ValueError:
            continue
        shorter = w[:j] + w[j + 1:]
        out[shorter] = padd(out.get(shorter, {}), pshift(p, j))
    return {w: p for w, p in out.items() if p}

def h0(n, color, vec):
    prev, cur = {}, vec
    for k in range(n):
        nxt = x0(color, cur)
        if k:
            nxt = vadd(nxt, prev, -1)
        prev, cur = cur, nxt
    return cur

def factor0(vec, n, sign, cap):
    out = {}
    for color, coefficient in ((0, 1), (1, sign)):
        h = h0(n, color, vec)
        h = {w: p for w, p in h.items() if len(w) <= cap}
        out = vadd(out, h, coefficient)
    return {w: p for w, p in out.items() if len(w) <= cap}

def q0_prefix_polynomial(A, B, cutoff, label):
    vec = {(): {0: 1}}
    weights = [n for n, sign in B]
    for i in reversed(range(len(B))):
        n, sign = B[i]
        cap = min(cutoff + sum(weights[:i]), sum(weights[i:]))
        vec = factor0(vec, n, sign, cap)
        stamp(f"{label} B-factor={len(B)-i}/{len(B)} n={n} cap={cap} states={len(vec)}")
    for i, (n, sign) in enumerate(A):
        cap = sum(nn for nn, _ in A[i + 1:])
        vec = factor0(vec, n, sign, cap)
        stamp(f"{label} A-factor={i+1}/{len(A)} n={n} cap={cap} states={len(vec)}")
    return vec.get((), {})

def eval_poly(p, x):
    return sum(Q(c) * x**j for j, c in p.items())

def bernstein_coefficients(p, right=Q(1)):
    d = max(p, default=0)
    return [
        sum(Q(p.get(j, 0)) * right**j * Q(comb(i, j), comb(d, j))
            for j in range(i + 1))
        for i in range(d + 1)
    ]

def cg(a, n):
    return range(abs(a - n), a + n + 1, 2)

def character_table(word):
    table = {(0, 0): 1}
    for n, sign in word:
        nxt = defaultdict(int)
        for (r, s), value in table.items():
            for c in cg(r, n):
                nxt[c, s] += value
            for c in cg(s, n):
                nxt[r, c] += sign * value
        table = {key: value for key, value in nxt.items() if value}
    return table

def su2_prefix(A, B, cutoff):
    x, y = character_table(A), character_table(B)
    return sum(value * y.get(key, 0)
               for key, value in x.items()
               if sum(key) <= cutoff)

def su2_layers(A, B):
    x, y = character_table(A), character_table(B)
    layers = defaultdict(int)
    for key, value in x.items():
        if key in y:
            layers[sum(key)] += value * y[key]
    return {t: v for t, v in sorted(layers.items()) if v}

def cut_interior(C):
    A, B = [], []
    wa = wb = 0
    for z in sorted(C, key=lambda x: (-abs(x), x)):
        if wa <= wb:
            A.append(z)
            wa += abs(z)
        else:
            B.append(z)
            wb += abs(z)
    if wa > wb:
        A, B = B, A
    return A, B

def top_pair(word):
    best_key = None
    best = None
    for i in range(len(word)):
        for j in range(i + 1, len(word)):
            if (abs(word[i]) + abs(word[j])) % 2:
                continue
            key = (abs(word[i]) + abs(word[j]),
                   max(abs(word[i]), abs(word[j])))
            if best_key is None or key > best_key:
                best_key, best = key, (i, j)
    return best

def inverse_crossings(u, v):
    if len(u) != len(v) or u.count(0) != v.count(0):
        return None
    positions = [[], []]
    for j, color in enumerate(v):
        positions[color].append(j)
    used = [0, 0]
    permutation = []
    for color in u:
        permutation.append(positions[color][used[color]])
        used[color] += 1
    return sum(permutation[i] > permutation[j]
               for i in range(len(u)) for j in range(i + 1, len(u)))

def full_q0_vector(word):
    vec = {(): {0: 1}}
    total = sum(n for n, _ in word)
    for n, sign in reversed(word):
        vec = factor0(vec, n, sign, total)
    return vec

def q0_prefix_polynomials_small(A, B):
    va, vb = full_q0_vector(A), full_q0_vector(B)
    layers = defaultdict(dict)
    for ua, pa in va.items():
        for ub, pb in vb.items():
            if len(ua) != len(ub):
                continue
            crossings = inverse_crossings(ua[::-1], ub)
            if crossings is None:
                continue
            term = pshift(pmul(pa, pb), crossings)
            layers[len(ua)] = padd(layers[len(ua)], term)
    prefixes, run = {}, {}
    for t in range(max(layers, default=-1) + 1):
        run = padd(run, layers[t])
        prefixes[t] = run
    return prefixes

def screen_small_interior():
    profiles = prefixes_checked = 0
    for length in range(3, 6):
        at_length = 0
        for magnitudes in combinations_with_replacement((1, 2, 3), length):
            if sum(magnitudes) % 2:
                continue
            for signs in product((1, -1), repeat=length):
                if sum(sign < 0 for sign in signs) % 2:
                    continue
                word = tuple(n * sign for n, sign in zip(magnitudes, signs))
                pair = top_pair(word[:-1])
                if pair is None:
                    continue
                i, j = pair
                C = [z for k, z in enumerate(word[:-1]) if k not in (i, j)]
                A0, B0 = cut_interior(C)
                if not A0 or not B0:
                    continue
                B0 = B0 + [word[i], word[j], word[-1]]
                A = [(abs(z), 1 if z > 0 else -1) for z in A0]
                B = [(abs(z), 1 if z > 0 else -1) for z in B0]
                prefixes = q0_prefix_polynomials_small(A, B)
                for cutoff, poly in prefixes.items():
                    assert min(bernstein_coefficients(poly)) >= 0
                    assert sum(poly.values()) == su2_prefix(A, B, cutoff)
                    prefixes_checked += 1
                profiles += 1
                at_length += 1
        stamp(f"small interior length={length} profiles={at_length} total={profiles} prefixes={prefixes_checked}")
    assert (profiles, prefixes_checked) == (144, 336)

stamp("begin exact braided q=0 prefix certificates")

seg_A = [(1, 1), (2, -1), (4, -1), (4, -1), (4, -1)]
seg_B = [(1, 1)] + [(3, 1)] * 6 + [(5, 1)] * 4 + [(6, 1)]
seg = q0_prefix_polynomial(seg_A, seg_B, 11, "segregated")
assert max(seg) == 205
assert seg.get(0) == 281136586
assert sum(seg.values()) == -24695910
assert eval_poly(seg, Q(1, 2)) > 0
assert eval_poly(seg, Q(3, 4)) > 0
left_bern = bernstein_coefficients(seg, Q(511, 512))
assert min(left_bern) > 0
assert eval_poly(seg, Q(1023, 1024)) < 0
assert eval_poly(seg, Q(1)) == -24695910
stamp("segregated Pi_11: degree=205; positive through s=511/512; negative at s=1023/1024 and s=1")

int_A = [(5, 1), (4, -1), (4, -1), (3, 1), (3, 1), (2, -1), (1, 1)]
int_B = [(5, 1), (4, -1), (3, 1), (3, 1), (3, 1), (3, 1),
         (1, 1), (5, 1), (5, 1), (6, 1)]
interior = q0_prefix_polynomial(int_A, int_B, 11, "interior")
assert max(interior) == 177
assert interior.get(0) == 138223104
assert interior.get(1) == -1840790
assert sum(interior.values()) == 117892066
interior_bern = bernstein_coefficients(interior)
assert min(interior_bern) == 117892066
stamp("interior TopPair (+5,+5) Pi_11: degree=177; all s-Bernstein coefficients positive; minimum=117892066")

seg_layers = su2_layers(seg_A, seg_B)
assert seg_layers == {
    1: 206613818, 3: 863409696, 5: 452240438,
    7: -552023508, 9: -663275890, 11: -331660464,
    13: 70846556, 15: 74390904
}
assert su2_prefix(seg_A, seg_B, 11) == -24695910
assert character_table(seg_A + seg_B).get((0, 0), 0) == 120541550

int_layers = su2_layers(int_A, int_B)
assert int_layers == {
    2: 7640446, 4: 19674282, 6: 28771926, 8: 34494764,
    10: 27310648, 12: 10274012, 14: 223774, 16: -3094970,
    18: -3233332, 20: -1360156, 22: -159844
}
running = 0
for t in range(max(int_layers) + 1):
    running += int_layers.get(t, 0)
    assert running >= 0
assert su2_prefix(int_A, int_B, 11) == 117892066
assert character_table(int_A + int_B).get((0, 0), 0) == 120541550
stamp("SU(2) endpoint: segregated Pi_11=-24695910, Phi=120541550; interior endpoint all prefixes nonnegative")

screen_small_interior()
stamp("PASS: exact path and small interior Bernstein checks")
