"""FM-MECH42 (astra_max_ceres): no positive-Walsh semicircle realization has U_3 Walsh-nonnegative (unique negative coefficient);
U_4 not automatic (negative at sigma_1 for the FM-MECH41 realization); (1,1,1,3) is B."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import comb

argparse.ArgumentParser(
    description="FM-MECH42: universal U3 obstruction and U4 correction"
).parse_args()


def conv(a, b):
    out = [Q(0)] * len(a)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i ^ j] += x * y
    return out


mu = {2*k: Q(comb(2*k, k), k+1) for k in range(4)}
assert (mu[2], mu[4], mu[6]) == (1, 2, 5)
assert mu[4] - 2*mu[2] == 0
assert mu[6] - 4*mu[4] + 4*mu[2] == 1
assert mu[4] - 3*mu[2] + 1 == 0
print("Semicircle moments: 1 2 5; <U1,U3>=0; ||U3||^2=1; E[U4]=0")

checks = 0
for rank in range(1, 6):
    size = 1 << rank
    for seed in range(8):
        a = [Q(0)] + [
            Q((g*g + seed*g + 3*seed) % 7, 7)
            for g in range(1, size)
        ]
        cubic = conv(conv(a, a), a)
        norm = sum(x*x for x in a)
        distinct = [Q(0)] * size
        for i, j, k in combinations(range(size), 3):
            distinct[i ^ j ^ k] += 6*a[i]*a[j]*a[k]
        for g in range(size):
            assert cubic[g] == 3*a[g]*norm - 2*a[g]**3 + distinct[g]
            assert distinct[g] >= 0
            checks += 1
        assert sum(a[g]*(cubic[g]-2*a[g]) for g in range(size)) == (
            3*norm**2 - 2*sum(x**4 for x in a) - 2*norm
            + sum(a[g]*distinct[g] for g in range(size))
        )
print("Cubic identity:", checks, "exact entries; 40 summed identities")

# Even moments of (epsilon_1+epsilon_2)/sqrt(2), exactly.
discrete = {
    2*k: sum(
        Q((s+t)**(2*k), 2**k)
        for s in (-1, 1) for t in (-1, 1)
    ) / 4
    for k in range(1, 4)
}
assert discrete == {2: Q(1), 4: Q(2), 6: Q(4)}
assert discrete[6] - 4*discrete[4] + 4*discrete[2] == 0
print("Forced two-character moments: 1 2 4; U3 squared norm: 0")


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(a, b):
    return trim([
        (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
        for i in range(max(len(a), len(b)))
    ])


def mul(a, b):
    p = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            p[i+j] += x*y
    return trim(p)


# Derivative of x*sqrt(4-x^2)*Q4(x)/(12*pi).
P4 = [1, 0, -3, 0, 1]
Q4 = [12, 0, -11, 0, 2]
dQ4 = [i*Q4[i] for i in range(1, len(Q4))]
lhs = add(mul([4, 0, -2], Q4), mul([0, 4, 0, -1], dQ4))
rhs = [12*x for x in mul([4, 0, -1], P4)]
assert lhs == rhs
assert add(Q4, [-3]) == mul([1, 0, -1], [9, 0, -2])
assert Q(4)**2 < 27  # pi < 4 < 3*sqrt(3).
print("U4 integral: derivative identity and "
      "Q4(x)-3=(1-x^2)(9-2*x^2) verified")


@lru_cache(None)
def moment(labels):
    row = {0: 1}
    for n in labels:
        nxt = {}
        for j, v in row.items():
            for k in range(abs(j-n), j+n+1, 2):
                nxt[k] = nxt.get(k, 0) + v
        row = nxt
    return row.get(0, 0)


labels = (1, 1, 1, 3)
tab = []
for S in range(1 << len(labels)):
    left = tuple(n for i, n in enumerate(labels) if S >> i & 1)
    right = tuple(n for i, n in enumerate(labels) if not (S >> i & 1))
    tab.append(moment(left) * moment(right))
assert tab == [int(S in (0, 15)) for S in range(16)]
print("B control (1,1,1,3): line indicator, 16 exact entries")
print("PASS")
