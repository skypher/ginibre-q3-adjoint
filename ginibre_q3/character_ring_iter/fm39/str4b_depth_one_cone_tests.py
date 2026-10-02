import sys
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import prod

if any(a in ("-h", "--help") for a in sys.argv[1:]):
    print("Verify exact channel-cone and flip certificates for two W<=36 residual profiles.")
    raise SystemExit(0)

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

@lru_cache(None)
def coeff2(L, A=0, B=0):
    d = {(0, 0): 1}
    for n, e in L:
        q = defaultdict(int)
        for (x, y), v in d.items():
            for t in cg(x, n):
                q[t, y] += v
            for t in cg(y, n):
                q[x, t] += e*v
        d = dict(q)
    return d.get((A, B), 0)

class Profile:
    def __init__(self, Babs, p, toplabels):
        self.B = tuple(Babs)
        self.p = p
        self.classes = tuple(sorted(set(Babs)))
        self.mult = {n: self.B.count(n) for n in self.classes}
        self.toppos = []
        for n in toplabels:
            self.toppos.append(next(i for i, z in enumerate(self.B)
                                    if z == n and i not in self.toppos))
        self.L = self.B + (p,)
        self.reps = {}
        for i, j in combinations(range(len(self.L)), 2):
            self.reps.setdefault((self.L[i], self.L[j]), (i, j))
        self.channels = []
        for (a, b), (i, j) in self.reps.items():
            for c in cg(a, b):
                self.channels.append((a, b, c, i, j))
        self.fliptypes = list(self.reps.items())

    def parent(self, mask):
        e = {n: (-1 if (mask >> j) & 1 else 1)
             for j, n in enumerate(self.classes)}
        sigma = prod(e[n]**self.mult[n] for n in self.classes)
        return tuple((n, e[n]) for n in self.B) + ((self.p, sigma),)

    def eps(self, mask):
        return {n: (-1 if (mask >> j) & 1 else 1)
                for j, n in enumerate(self.classes)}

    def g(self, mask):
        return coeff2(self.parent(mask)[:-1], self.p, 0)

    def gchild(self, mask):
        e = self.eps(mask)
        C = tuple((n, e[n]) for i, n in enumerate(self.B)
                  if i not in self.toppos)
        return coeff2(C, self.p, 0)

    def delta(self, mask):
        return self.g(mask) - self.gchild(mask)

    def child(self, mask, i, j, c):
        w = self.parent(mask)
        C = tuple(w[k] for k in range(len(w)) if k not in (i, j))
        return coeff2(C + ((c, w[i][1]*w[j][1]),), 0, 0)

    def flipD(self, mask, i, j):
        w = self.parent(mask)
        q = list(w)
        q[i] = (q[i][0], -q[i][1])
        q[j] = (q[j][0], -q[j][1])
        return Fraction(coeff2(w) - coeff2(tuple(q)), 4)

    def tables(self):
        rows = range(1 << len(self.classes))
        delta = [self.delta(r) for r in rows]
        children = [
            [self.child(r, i, j, c)
             for a, b, c, i, j in self.channels]
            for r in rows
        ]
        flips = [
            [self.flipD(r, *ij) for _, ij in self.fliptypes]
            for r in rows
        ]
        noflip = [r for r in rows if max(flips[r]) < 0]
        return list(rows), delta, children, flips, noflip

P35 = Profile([1, 2, 3, 3, 3, 4, 4, 5, 5, 5], 15, (5, 5))
r35, d35, C35, D35, nf35 = P35.tables()
assert (sum(P35.B), len(P35.B), len(r35), len(P35.channels)) == (35, 10, 32, 65)
assert (min(d35), max(d35), nf35) == (13351, 48443, [10, 31])

i, j = P35.reps[(1, 2)]
col35 = P35.channels.index((1, 2, 1, i, j))
assert all(Fraction(13351, 13692)*C35[r][col35] == d35[r]
           for r in nf35)

y35 = {
    0: Fraction(-29, 3951),
    1: Fraction(-109, 8624),
    4: Fraction(-93, 8201),
    10: Fraction(530, 9651),
}
yc35 = [
    sum(y35.get(r, Fraction(0))*C35[r][k] for r in r35)
    for k in range(len(P35.channels))
]
yd35 = sum(y35.get(r, Fraction(0))*d35[r] for r in r35)
yflip35 = [
    sum(y35.get(r, Fraction(0))*D35[r][k] for r in r35)
    for k in range(len(P35.fliptypes))
]
assert min(yc35) == Fraction(1589995915, 64210435376472)
assert yd35 == Fraction(-898902791216245, 898946095270608)
assert max(yflip35) == Fraction(-3372984706561879, 74912174605884)
print("W35: Delta range", min(d35), max(d35), "no-flip masks", nf35)
print("W35 no-flip identity: (1,2)->1, beta=13351/13692")
print("W35 all-sign Farkas: min y.C=", min(yc35),
      "y.Delta=", yd35, "max y.D=", max(yflip35))

P36 = Profile([1]*18 + [2] + [3]*4 + [4], 6, (2, 4))
r36, d36, C36, D36, nf36 = P36.tables()
assert (sum(P36.B), len(P36.B), len(r36), len(P36.channels)) == (36, 24, 16, 36)
assert (min(d36), max(d36), [r for r in r36 if d36[r] < 0], nf36) == (
    -5074321295, 3361824130625, [10, 15], [])
assert min(C36[10]) == 2903493320 and all(v >= 0 for v in C36[10])

alpha = Fraction(
    61161619770999573614531252672130373042966611715500020053689725643936449,
    2512283634690871622407588119211562466307544457373940296651740817514128,
)
betas = [
    ((1, 1, 2), Fraction(
        1981583616294495658365023836050550989535922384061919465598025795921051,
        157017727168179476400474257450722654144221528585871268540733801094633)),
    ((1, 2, 1), Fraction(
        980341445938533881340917848792015050701507029731984958587303178035741,
        314035454336358952800948514901445308288443057171742537081467602189266)),
    ((1, 3, 2), Fraction(
        1238098426960579128923989555434833516795722832989676337125298891962293,
        314035454336358952800948514901445308288443057171742537081467602189266)),
    ((1, 4, 3), Fraction(
        405980364940012697428380982714646109525762424407718121500128989109481,
        314035454336358952800948514901445308288443057171742537081467602189266)),
    ((2, 3, 1), Fraction(
        138752901891742185573374510432031422632643197542506350812497737138747,
        157017727168179476400474257450722654144221528585871268540733801094633)),
    ((3, 4, 1), Fraction(
        428837862724535348112247588024534774395858862451616261922713658536651,
        314035454336358952800948514901445308288443057171742537081467602189266)),
    ((3, 6, 5), Fraction(
        13815515046560753382081341891490379722780268894108181800585163418453,
        314035454336358952800948514901445308288443057171742537081467602189266)),
]
ii, jj = P36.reps[(1, 1)]
for r in r36:
    rhs = Fraction(0)
    for (a, b, c), beta in betas:
        u, v = P36.reps[(a, b)]
        rhs += beta*P36.child(r, u, v, c)
    assert Fraction(d36[r]) + alpha*P36.flipD(r, ii, jj) == rhs
print("W36: Delta range", min(d36), max(d36),
      "negative masks [10, 15], no-flip rows", nf36)
print("W36 channel-only obstruction at mask 10: Delta=", d36[10],
      "children range", min(C36[10]), max(C36[10]))
print("W36 all-sign flip identity: one D_(1,1) plus seven channel children")
print("ALL EXACT CHECKS PASS")
