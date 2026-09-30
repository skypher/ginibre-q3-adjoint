"""FM-SEC88 (luna_max_venus): H_AC on the {1,2} sector: profile formula, small bi-radial certificates,
the bi-radial obstruction at (1^4,2^3), a 29-subgroup H_AC certificate there, subgroup-square boundary certificates."""
from math import comb
from fractions import Fraction as F
from sympy import Matrix, symbols, factor, zeros

def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def Cat(n):
    return C(2*n, n) // (n+1)

def m(i, j):
    if i % 2:
        return 0
    return sum((-1)**(j-t) * C(j, t) * Cat(i//2+t)
               for t in range(j+1))

def fusion(labels):
    d = {0: 1}
    for n in labels:
        e = {}
        for x, v in d.items():
            for y in range(abs(x-n), x+n+1, 2):
                e[y] = e.get(y, 0) + v
        d = e
    return d.get(0, 0)

for i in range(13):
    for j in range(13-i):
        assert m(i, j) == fusion((1,)*i + (2,)*j)
print("m(i,j) equals SU(2) fusion multiplicity for i+j <= 12")

def types(N, M):
    return sorted({min((i, j), (N-i, M-j))
                   for i in range(N+1) for j in range(M+1)})

def matrices(N, M):
    T = types(N, M)
    ix = {v: k for k, v in enumerate(T)}
    d = len(T)
    Q = []
    for a, b in T:
        A = [[F(0) for _ in T] for __ in range(d)]
        for u in range(N+1):
            for v in range(M+1):
                for h in range(max(0, u-(N-a)), min(a, u)+1):
                    for ell in range(max(0, v-(M-b)), min(b, v)+1):
                        count = F(
                            C(a, h)*C(N-a, u-h) *
                            C(b, ell)*C(M-b, v-ell), 2)
                        i = ix[min((u, v), (N-u, M-v))]
                        j = ix[min(
                            (a+u-2*h, b+v-2*ell),
                            (N-a-u+2*h, M-b-v+2*ell))]
                        A[i][j] += count
        Q.append(Matrix([
            [(A[i][j]+A[j][i])/2 for j in range(d)]
            for i in range(d)
        ]))
    target = Matrix([m(i, j)*m(N-i, M-j) for i, j in T])
    return T, Q, target

for N, M in ([(10, 1), (12, 1)] +
             [(N, 2) for N in range(2, 13, 2)] +
             [(0, M) for M in range(1, 11)]):
    T, _, f = matrices(N, M)
    print("profile", (N, M), list(zip(T, list(f))))

# Exact bi-radial certificates: each pair lists (profile p, coefficient).
small_certs = {
    (0, 1): [],
    (0, 2): [((0,1), F(1))],
    (0, 3): [((1,0), F(1))],
    (0, 4): [((0,0,1), F(1,2)), ((1,0,0), F(3,2))],
    (0, 5): [((0,1,0), F(1,2)), ((1,0,0), F(7,2))],
    (0, 6): [((1,0,0,0), F(19,2)),
             ((1,0,0,1), F(1,2))],
    (2, 0): [((0,1), F(1))],
    (2, 1): [((0,1,0), F(1))],
    (2, 2): [((0,0,0,0,1), F(1,2)),
             ((0,0,1,0,0), F(1))],
    (2, 3): [((0,0,0,1,0,0), F(3,2)),
             ((0,0,1,0,0,0), F(1,2)),
             ((1,0,0,1,0,0), F(1,2))],
    (2, 4): [((0,0,0,0,0,1,0,0), F(3,2)),
             ((0,0,0,0,1,0,0,0), F(3,2)),
             ((0,0,0,1,0,0,0,0), F(1,2)),
             ((0,1,0,0,1,0,0,0), F(1,2))],
    (4, 0): [((0,0,1), F(1,2)), ((1,0,0), F(1,2))],
    (4, 1): [((0,0,0,1,0), F(1,2)),
             ((0,1,0,0,0), F(1))],
    (4, 2): [((0,0,0,0,0,1,0,0), F(1,2)),
             ((0,0,1,0,0,0,0,0), F(2)),
             ((0,1,0,0,0,0,1,0), F(1,4))],
    (6, 0): [((0,0,0,1), F(1,8)),
             ((0,1,0,0), F(5,8))],
    (0, 7): [((0,0,1,0), F(1,2)),
             ((0,1,0,0), F(1,2)),
             ((1,0,0,0), F(22))],
    (2, 5): [((0,0,0,0,0,0,0,1,0), F(1,4)),
             ((0,0,0,0,0,0,1,0,0), F(7,4)),
             ((0,0,0,0,0,1,0,0,0), F(19,2)),
             ((0,0,1,0,0,1,0,0,0), F(1,2))],
    (6, 1): [((0,0,0,1,0,0,0), F(1,2)),
             ((0,1,0,0,0,0,0), F(2)),
             ((0,1,0,0,1,0,0), F(1,4))]
}
for (N, M), cert in small_certs.items():
    T, Q, target = matrices(N, M)
    got = Matrix([
        sum(c*(Matrix(p).T*Q[s]*Matrix(p))[0] for p, c in cert)
        for s in range(len(T))
    ]) if cert else Matrix([0]*len(T))
    assert got == target
    assert all(c > 0 for _, c in cert)
print("all listed bi-radial certificates re-expand exactly")

# (4,3) support obstruction.
N, M = 4, 3
T, Q, target = matrices(N, M)
assert list(target) == [13, 0, 3, 2, 0, 0, 0, 0, 4, 2]
d = len(T)
zero = [s for s, v in enumerate(target) if v == 0]
compat = [[all(Q[s][i,j] == 0 for s in zero)
           for j in range(d)] for i in range(d)]
cliques = []
for mask in range(1, 1 << d):
    V = [i for i in range(d) if (mask >> i) & 1]
    clique = all(compat[i][j] for i in V for j in V)
    maximal = not any(
        all(compat[i][j] for i in V+[k] for j in V+[k])
        for k in range(d) if k not in V)
    if clique and maximal:
        cliques.append(V)
assert cliques == [[4,6], [4,7], [5,7],
                   [0,2,8], [0,3,8], [1,3,8]]
y = [0,0,3,1,0,0,0,0,0,-6]
assert sum(y[s]*target[s] for s in range(d)) == -1
My = sum((y[s]*Q[s] for s in range(d)), zeros(d))
p = symbols("p0:10")
for V in cliques:
    expr = factor(sum(My[i,j]*p[i]*p[j] for i in V for j in V))
    print("support", V, "quadratic", expr)
print("bi-radial separator dot target =", -1)

# Full H_AC certificate on G ~= F_2^6.
terms = [
(F(108,121),(0,57)), (F(148,121),(0,10,35,41)),
(F(23,121),(0,5,26,31)), (F(17,121),(0,31,38,57)),
(F(39,121),(0,19,37,54)), (F(180,121),(0,5,9,12)),
(F(40,121),(0,19,35,48)), (F(87,121),(0,22,42,60)),
(F(199,121),(0,3,5,6,9,10,12,15)),
(F(16,121),(0,6,10,12,19,21,25,31)),
(F(43,121),(0,3,12,15,21,22,25,26)),
(F(1,242),(0,12,22,26,35,47,53,57)),
(F(7,242),(0,3,9,10,21,22,28,31)),
(F(39,242),(0,6,26,28,41,47,51,53)),
(F(31,121),(0,5,25,28,42,47,51,54)),
(F(6,11),(0,3,28,31,37,38,57,58)),
(F(139,242),(0,6,25,31,42,44,51,53)),
(F(101,242),(0,3,21,22,44,47,57,58)),
(F(203,242),(0,9,19,26,38,47,53,60)),
(F(8,121),(0,6,25,31,35,37,58,60)),
(F(47,121),(0,3,5,6,48,51,53,54)),
(F(25,121),(0,6,25,31,41,47,48,54)),
(F(49,121),(0,10,21,31,37,47,48,58)),
(F(69,242),(0,3,21,22,37,38,48,51)),
(F(81,242),(0,3,28,31,44,47,48,51)),
(F(4,121),(0,5,25,28,41,44,48,53)),
(F(23,121),(0,10,22,28,38,44,48,58)),
(F(109,242),(0,6,26,28,42,44,48,54)),
(F(91,242),(0,6,10,12,19,21,25,31,35,37,41,47,48,54,58,60))
]
def subgroup(W):
    S = set(W)
    return 0 in S and all((x ^ z) in S for x in S for z in S)
assert all(c > 0 and subgroup(W) for c, W in terms)
for mask in range(64):
    i = (mask & 15).bit_count()
    j = ((mask >> 4) & 7).bit_count()
    target_value = m(i,j)*m(4-i,3-j)
    got = sum(c for c, W in terms if mask in W)
    assert got == target_value, (mask, target_value, got)
print("full (4,3) H_AC certificate: 29 subgroups; exact on all 64 masks")


# ---- block ----
from math import comb
from fractions import Fraction as F

def C(n,k):
    return comb(n,k) if 0 <= k <= n else 0

def Cat(n):
    return C(2*n,n)//(n+1)

def m(i,j):
    if i % 2:
        return 0
    return sum((-1)**(j-t)*C(j,t)*Cat(i//2+t)
               for t in range(j+1))

codes = {
 10: [
  (F(1485,32),(1,735)), (F(51,16),(7,280)),
  (F(525,32),(7,594)), (F(45,4),(7,55,881,618)),
  (F(51,4),(31,962,743,766,130,710))
 ],
 12: [
  (F(14157,128),(1,4088)),
  (F(94545,2368),(7,1997,2557,3050)),
  (F(126225,1184),(7,1592,3582,60,4029,192)),
  (F(146421,4736),(31,)),
  (F(10395,2368),(31,3972,4090)),
  (F(10791,2368),(31,4075,2063))
 ]
}
for N, terms in codes.items():
    prof = [F(0)]*(N//2+1)
    for coeff, checks in terms:
        assert coeff > 0
        for x in range(1 << N):
            if (x.bit_count() % 2 == 0 and
                all((x & h).bit_count() % 2 == 0 for h in checks)):
                prof[x.bit_count()//2] += coeff
    prof = [v/C(N,2*j) for j,v in enumerate(prof)]
    target = [F(m(2*j,0)*m(N-2*j,1))
              for j in range(N//2+1)]
    assert prof == target
    print(N, prof, "exact subgroup-square profile")
