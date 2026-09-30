# FM-SEC66 (luna_max_neptune) printed verifier: level-1 two-hat certificate at suffix a = 2 (39-term orbit identity).
import itertools
from fractions import Fraction as Fr
import sympy as sp
from n0check import ONE, h, shat, mul, phi

CERT = [
 ('3/2',(0,0,0,0),(0,)), ('13/7',(0,0,0,0),(0,1)),
 ('9/14',(0,0,0,0),(1,2)), ('6/7',(0,0,0,2),(0,1,3)),
 ('4',(0,0,0,2),(0,1,4)), ('8/7',(0,0,0,2),(1,3)),
 ('2',(0,0,0,2),(2,4)), ('4',(0,0,1,1),(0,1,4)),
 ('4',(0,0,1,1),(1,3)), ('13/7',(0,1,0,1),(0,1,2)),
 ('1/7',(0,1,0,1),(0,1,2,3)), ('6/7',(0,1,0,1),(0,1,4)),
 ('16/7',(0,1,0,1),(0,3,4)), ('8/7',(0,1,0,1),(1,2)),
 ('2',(0,1,0,1),(1,2,4)), ('6/7',(0,1,0,1),(1,3,4)),
 ('2',(0,1,0,1),(2,4)), ('3/7',(0,2,0,0),(0,1)),
 ('6/7',(0,2,0,0),(1,)), ('19/14',(0,2,0,0),(1,2)),
 ('19/14',(0,2,0,0),(5,)), ('15/14',(0,3,0,1),(0,)),
 ('9/14',(0,3,0,1),(0,1)), ('15/7',(0,3,0,1),(0,1,5)),
 ('1/7',(0,3,0,1),(0,5)), ('22/7',(0,3,0,1),(1,2)),
 ('6/7',(0,3,0,1),(1,5)), ('11/7',(0,4,0,0),(0,)),
 ('3/7',(0,4,0,0),(0,1)), ('3/7',(0,4,0,0),(1,)),
 ('11/7',(0,4,0,0),(1,2)), ('3/14',(1,1,0,0),(0,1,2)),
 ('15/7',(1,1,0,0),(0,1,2,5)), ('3/7',(1,3,0,0),(0,1,2)),
 ('1/7',(2,2,0,0),(1,4,5)), ('17/14',(2,2,0,0),(5,)),
 ('3/7',(2,4,0,0),(3,)), ('3/7',(3,3,0,0),(1,4)),
 ('11/7',(3,3,0,0),(5,))
]

x = sp.symbols('t1 t2 u v z1 z2')
t1, t2, u, v, z1, z2 = x
e = lambda k: sum(sp.prod(x[i] for i in I)
                  for I in itertools.combinations(range(6), k))
N4 = sp.expand(
    (2+2*t1**2)*(2+2*t2**2)*(1-e(6))
    -(t1*(2+2*t2**2)+t2*(2+2*t1**2))*(e(1)-e(5))
    +t1*t2*(1+e(1)**2-e(4)-e(1)*e(5)+e(2)*e(6)-e(6)**2)
)
N00 = sp.expand(N4.subs({z1: 0, z2: 0}))
N10 = sp.expand(sp.diff(N4, z1).subs({z1: 0, z2: 0}))
N01 = sp.expand(sp.diff(N4, z2).subs({z1: 0, z2: 0}))
N11 = sp.expand(sp.diff(sp.diff(N4, z1), z2).subs({z1: 0, z2: 0}))
N2 = sp.Poly(sp.expand(
    N00*(1+(t1+t2+u+v)**2)
    +(N10+N01)*(t1+t2+u+v)+N11
), t1, t2, u, v)
target = {m: Fr(int(c)) for m, c in N2.terms()}

edges = list(itertools.combinations(range(4), 2))
edge_id = {edge: i for i, edge in enumerate(edges)}
perms = [(0,1,2,3), (1,0,2,3), (0,1,3,2), (1,0,3,2)]

def transform(pair, perm):
    m, F = pair
    mt = tuple(m[perm.index(i)] for i in range(4))
    Ft = tuple(sorted(
        edge_id[tuple(sorted((perm[i], perm[j])))]
        for k in F for i, j in [edges[k]]
    ))
    return mt, Ft

out = {}
for coeff, m, F in CERT:
    q = Fr(coeff)
    orbit = sorted({transform((m, F), p) for p in perms})
    for mm, FF in orbit:
        for mask in range(1 << len(FF)):
            exp = list(mm)
            sign = 1
            for j, ei in enumerate(FF):
                if mask >> j & 1:
                    i, k = edges[ei]
                    exp[i] += 1
                    exp[k] += 1
                    sign = -sign
            key = tuple(exp)
            out[key] = out.get(key, Fr(0)) + q * sign / len(orbit)

assert {k:v for k,v in out.items() if v} == {
    k:v for k,v in target.items() if v
}
print("exact orbit certificate:", len(CERT), "terms;",
      "N2 terms:", len(target), "degree:", N2.total_degree())

# Compare N2/D4 to direct Catalan moments through total degree 12.
D = 12
den = {(0,0,0,0): 1}
for i, j in itertools.combinations(range(4), 2):
    nxt = {}
    for exp, c in den.items():
        for k in range((D-sum(exp))//2 + 1):
            q = list(exp)
            q[i] += k
            q[j] += k
            q = tuple(q)
            nxt[q] = nxt.get(q, 0) + c
    den = nxt

nd = {m:int(c) for m,c in N2.terms()}
series = {}
for m, c in nd.items():
    for d, b in den.items():
        q = tuple(i+j for i,j in zip(m,d))
        if sum(q) <= D:
            series[q] = series.get(q, 0) + c*b

checks = 0
for P in range(2, D+1):
    for Q in range(2, P+1):
        for a in range(D+1):
            for b in range(D+1):
                if P+Q+a+b > D:
                    continue
                word = mul(mul(h(a), h(b)), mul(shat(P), shat(Q)))
                direct = phi(1, word, 2)
                assert direct == series.get((P,Q,a,b), 0)
                checks += 1
print("N2/D4 vs Catalan evaluator:", checks, "coefficients agree")
