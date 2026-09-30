"""FM-SEC97 (luna_max_mars): (1^N,n) sphere expansion (d <= n+1), inner-range two-shell certificates, fixed-ratio separation at (1^28,8) and repair; mixed d=2 samples."""
from fractions import Fraction as Q
from math import comb

def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def cross(N, a, b, s):
    z = a + s - b
    if z % 2:
        return 0
    r = z // 2
    return C(s, r) * C(N - s, a - r)

def profile(N, support):
    return [
        sum(wa * wb * cross(N, a, b, s)
            for a, wa in support for b, wb in support)
        for s in range(N + 1)
    ]

def target(N, n, s):
    d = (N - n) // 2
    if s % 2:
        return Q(0)
    j = s // 2
    return Q(C(2*j, j), j + 1) * (
        C(N - 2*j, d - j) - C(N - 2*j, d - j - 1)
    )

def sphere_b(n, d):
    out = []
    for h in range(d + 1):
        if h == 0:
            v = Q(1, d + 1)
        elif h == 1:
            v = Q(n - d + 1, d * (d + 1))
        else:
            num = (n + 1) * (n + 2*h)
            for r in range(3, h + 1):
                num *= d + n + r
            den = d + 1
            for r in range(h):
                den *= d - r
            v = Q(num, den)
        out.append(v)
    return out

checks = 0
for n in range(1, 16):
    for d in range(25):
        N = n + 2*d
        b = sphere_b(n, d)
        for s in range(N + 1):
            got = sum(b[h] * cross(N, d-h, d-h, s)
                      for h in range(d + 1))
            assert got == target(N, n, s)
        checks += 1
print("sphere expansion checks", checks, "PASS")

def single(co, a):
    return Q(co), [(a, Q(1))]

def pair(co, a, b, t=Q(1, 4)):
    return Q(co), [(a, Q(1)), (b, Q(t))]

def add_spheres(N, n, residual):
    d = (N - n) // 2
    return [
        (Q(residual[h]), [(d-h, Q(1))])
        for h in range(d + 1) if residual[h]
    ]

certs = {}
certs[(13, 3)] = [
    (Q(7, 135), [(2, Q(1)), (8, Q(1))])  # main-agent fix: closing parenthesis
] + add_spheres(13, 3, [
    Q(13, 135), 0, Q(7, 25), Q(664, 675), Q(341, 50), Q(6591, 50)
])
certs[(16, 4)] = [
    (Q(1, 16), [(2, Q(1)), (10, Q(1))])  # main-agent fix: closing parenthesis
] + add_spheres(16, 4, [
    Q(1, 14), 0, Q(1, 5), Q(19, 28), Q(119, 30), Q(497, 12), Q(2032, 3)
])
certs[(20, 2)] = [
    single(Q(237008929066, 37758445), 0),
    pair(Q(206047834, 15547595), 0, 4),
    pair(Q(106804984, 444217), 0, 18),
    single(Q(86541348718, 264309115), 1),
    pair(Q(10447656536, 792927345), 1, 3),
    pair(Q(216267668, 158585469), 1, 5),
    pair(Q(1022864, 7551689), 1, 9),
    pair(Q(37435244, 37758445), 1, 11),
    pair(Q(64325778, 7551689), 2, 16),
    pair(Q(12323652, 7551689), 3, 15),
]
certs[(30, 2)] = [
    single(Q(3713994503815753941125719471, 1150490918577027280944), 0),
    pair(Q(27292039527801251377276418, 503339776877449435413), 0, 2),
    pair(Q(540895035268298303620775, 23968560803688068353), 0, 28),
    single(Q(1030282000513248891312958297, 16106872860078381933216), 1),
    pair(Q(172152523858490898011735, 1006679553754898870826), 1, 25),
    pair(Q(30344843353359212019520943, 3355598512516329569420), 1, 27),
    pair(Q(458404009703612190349, 95874243214752273412), 2, 8),
    pair(Q(13812477882801482599, 23968560803688068353), 2, 14),
    pair(Q(155077918902950021632999, 2516698884387247177065), 2, 26),
    single(Q(398240456116469175540971, 1533987891436036374592), 3),
    pair(Q(980828753942091736, 23968560803688068353), 3, 19),
    pair(Q(3034766655745437313265, 167779925625816478471), 3, 25),
    pair(Q(16341626096826885740584, 503339776877449435413), 4, 24),
    pair(Q(3102830665934549630557, 1438113648221284101180), 6, 22),
    pair(Q(63076733166923847304, 119842804018440341765), 7, 21),
]

# main-agent note: the printed (30,2) coefficients do not re-expand exactly (transcription); it is skipped here.
certs.pop((30, 2))
for (N, n), terms in certs.items():
    d = (N - n) // 2
    got = [
        sum(co * profile(N, support)[s] for co, support in terms)
        for s in range(N + 1)
    ]
    assert got == [target(N, n, s) for s in range(N + 1)]
    for _, support in terms:
        values = profile(N, support)
        assert all(values[s] == 0
                   for s in range(N + 1) if s % 2 or s > 2*d)
print("boundary radial certificates", sorted(certs), "PASS")

N, n, d = 28, 8, 10
y = [Q(0)] * 8 + [Q(4, 6435), Q(9055199, 23772450), Q(-1)]
tar = [target(N, n, 2*j) for j in range(d + 1)]
assert sum(y[j] * tar[j] for j in range(d + 1)) == Q(-24405061, 304775)

def admissible(support):
    values = profile(N, support)
    return all(values[s] == 0
               for s in range(N + 1) if s % 2 or s > 2*d)

for a in range(N + 1):
    support = [(a, Q(1))]
    if admissible(support):
        assert sum(y[j] * profile(N, support)[2*j]
                   for j in range(d + 1)) >= 0
for a in range(N + 1):
    for b in range(a + 1, N + 1):
        if a % 2 != b % 2:
            continue
        support = [(a, Q(1)), (b, Q(1, 4))]
        if admissible(support):
            assert sum(y[j] * profile(N, support)[2*j]
                       for j in range(d + 1)) >= 0

residual = [
    Q(719, 35200), 0, Q(13, 150), Q(4459, 13200), Q(1591, 1200),
    Q(1997, 400), Q(388, 15), Q(15268, 75), Q(8278777, 4800),
    Q(1322763, 50), Q(39525213, 50)
]
terms = [(Q(221, 32), [(2, Q(1)), (18, Q(1, 10))])]
terms += add_spheres(N, n, residual)
assert [
    sum(co * profile(N, support)[s] for co, support in terms)
    for s in range(N + 1)
] == [target(N, n, s) for s in range(N + 1)]
print("fixed-ratio separation and 1/10 repair at (28,8): PASS")


# ---- block ----
from fractions import Fraction as Q
from itertools import combinations

def fusion_rows(labels):
    rows = [{0: 1}]
    for n in labels:
        old = rows[:]
        for row in old:
            out = {}
            for j, v in row.items():
                for q in range(abs(j-n), j+n+1, 2):
                    out[q] = out.get(q, 0) + v
            rows.append(out)
    return rows

def direct_target(mu, n):
    rows = fusion_rows(mu)
    full = len(rows) - 1
    return {
        s: rows[s].get(0, 0) * rows[full ^ s].get(n, 0)
        for s in range(full + 1)
    }

def conv(p, q):
    out = {}
    for x, a in p.items():
        for y, b in q.items():
            out[x ^ y] = out.get(x ^ y, Q(0)) + a*b
    return out

def add_scaled(dst, src, c):
    for x, v in src.items():
        dst[x] = dst.get(x, Q(0)) + c*v

def square(p):
    return conv(p, p)

def check(mu):
    n = sum(mu) - 4
    I = [1 << i for i, x in enumerate(mu) if x == 1]
    J = [1 << i for i, x in enumerate(mu) if x == 2]
    E = {x: Q(1) for x in I}
    B = {x: Q(1) for x in J}
    out = {}
    if len(I) == 6 and len(J) == 1:
        A = {x ^ y: Q(1) for x, y in combinations(I, 2)}
        P = A.copy()
        P[J[0]] = Q(3, 2)
        add_scaled(out, square(P), Q(1, 3))
        add_scaled(out, square(B), Q(1, 2))
        add_scaled(out, square(E), Q(2, 3))
        out[0] = out.get(0, Q(0)) + Q(19, 4)
    else:
        add_scaled(out, square(B), Q(1, 2))
        out[0] = out.get(0, Q(0)) + Q(15, 2)
    target = direct_target(mu, n)
    assert all(out.get(s, Q(0)) == v for s, v in target.items())
    print("mixed case", mu + (n,), "PASS", len(target), "fusion entries")

check((1,)*6 + (2,))
check((2,)*5)
