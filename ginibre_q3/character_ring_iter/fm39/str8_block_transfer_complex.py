from functools import lru_cache
from math import comb
from itertools import combinations, permutations, product
import sympy as s

@lru_cache(None)
def invdim(labels):
    d = {0:1}
    for n in labels:
        nxt = {}
        for j,v in d.items():
            for k in range(abs(j-n),j+n+1,2):
                nxt[k] = nxt.get(k,0)+v
        d = nxt
    return d.get(0,0)

def chain_dims(labels,minus):
    out = [0]*(len(minus)+1)
    for mask in range(1<<len(labels)):
        left = tuple(n for i,n in enumerate(labels) if mask>>i&1)
        right = tuple(n for i,n in enumerate(labels) if not(mask>>i&1))
        k = sum(mask>>i&1 for i in minus)
        out[k] += invdim(left)*invdim(right)
    return out

def euler(v): return sum((-1)**k*n for k,n in enumerate(v))
small = chain_dims((1,1,2,3,3,4),set(range(6)))
assert small == [11,0,4,10,4,0,11]
assert small[3]-small[2]-small[4] == 2
assert euler(small) == 20

# The odd-block complex repairs this grading obstruction.
ns = (1,1,2,3,3,4)
def triangle_form(A):
    i,j,k = A
    aa,bb,cc = ns[i]+ns[j]-ns[k],ns[j]+ns[k]-ns[i],ns[k]+ns[i]-ns[j]
    if min(aa,bb,cc)<0 or (aa|bb|cc)&1: return {}
    out = {(0,)*6:1}
    for u,v,e in ((i,j,aa//2),(j,k,bb//2),(k,i,cc//2)):
        nxt = {}
        for t,c in out.items():
            for h in range(e+1):
                z = list(t); z[u]+=h; z[v]+=e-h; z=tuple(z)
                nxt[z] = nxt.get(z,0)+c*(-1)**h*comb(e,h)
        out = {t:c for t,c in nxt.items() if c}
    return out
columns = []
for A in combinations(range(6),3):
    B = tuple(i for i in range(6) if i not in A)
    left,right = triangle_form(A),triangle_form(B)
    if not left or not right: continue
    col = {}
    for t,c in left.items():
        for u,d in right.items():
            v = tuple(i+j for i,j in zip(t,u))
            col[v] = col.get(v,0)+c*d
    columns.append(col)
rows = sorted({t for col in columns for t in col})
Ws = s.Matrix([[col.get(t,0) for col in columns] for t in rows])
assert Ws.shape == (130,10) and Ws.rank() == 5
small_h = [6,0,4,0,4,0,6]
assert euler(small_h) == 20

# Check the exterior signs for mixed plus/minus blocks.
def subsets(v):
    v = sorted(v)
    return [frozenset(q) for k in range(len(v)+1) for q in combinations(v,k)]
def sign(S,R,M):
    return (-1)**sum(t<r for r in R&M for t in (S-R)&M)
M = frozenset((0,2,4))
for S in subsets(range(5)):
    for R in subsets(S):
        if len(R&M)%2 == 0: continue
        for T in subsets(S-R):
            if len(T&M)%2 == 0: continue
            assert sign(S,R,M)*sign(S-R,T,M) == -sign(S,T,M)*sign(S-T,R,M)

# Six adjoint factors: delta-pairing basis of Inv(V_2^tensor6).
def pairings(v):
    if not v: return [()]
    out = []
    for b in v[1:]:
        for q in pairings(tuple(t for t in v[1:] if t != b)):
            out.append(tuple(sorted(((v[0],b),)+q)))
    return out
P = pairings(tuple(range(6)))
index = {p:i for i,p in enumerate(P)}
def loops(p,q):
    parent = list(range(6))
    def root(v):
        while parent[v] != v: v = parent[v]
        return v
    for a,b in p+q: parent[root(a)] = root(b)
    return len({root(i) for i in range(6)})
G = s.Matrix([[3**loops(p,q) for q in P] for p in P])
assert invdim((2,)*6) == 15
assert G.det() == 16070775840000000000

triples = list(combinations(range(6),3))
tindex = {A:i for i,A in enumerate(triples)}
W = s.zeros(15,20)
for j,A in enumerate(triples):
    B = tuple(i for i in range(6) if i not in A)
    for pi in permutations(range(3)):
        sg = (-1)**sum(pi[i]>pi[j] for i in range(3) for j in range(i+1,3))
        p = tuple(sorted(tuple(sorted((A[i],B[pi[i]]))) for i in range(3)))
        W[index[p],j] += sg
assert W.rank() == 5
assert G*W == 6*W
K = W.T*G*W
assert set(K) == {-12,12,36}
assert K.rank() == 5

D3 = W
D6 = s.zeros(20,15)
for j,A in enumerate(triples):
    R = tuple(i for i in range(6) if i not in A)
    D6[j,:] = (-1)**(sum(R)-3)*(W[:,j].T*G)/36
assert D3*D6 == s.zeros(15)
assert (D3.rank(),D6.rank()) == (5,5)
c = chain_dims((2,)*6,set(range(6)))
assert c == [15,0,45,20,45,0,15]
h = [10,0,45,10,45,0,10]
assert euler(c) == euler(h) == 100

# An explicit odd cycle which cannot be a boundary.
z = s.zeros(20,1)
J = s.zeros(20)
for j,A in enumerate(triples):
    B = tuple(i for i in range(6) if i not in A)
    J[tindex[B],j] = 1
for A,sg in [((0,1,2),1),((0,1,3),-1),((0,1,4),1),((0,1,5),-1)]:
    B = tuple(i for i in range(6) if i not in A)
    z[tindex[A]] += sg
    z[tindex[B]] += sg
assert z != s.zeros(20,1)
assert D3*z == s.zeros(15,1)
assert J*z == z and J*D6 == -D6
assert D6.row_join(z).rank() == 6

# All 180 nested one-factor projection maps in these nonzero spaces vanish.
checks = 0
for j,A in enumerate(triples):
    for pair in combinations(A,2):
        for rest in pairings(tuple(i for i in range(6) if i not in pair)):
            p = tuple(sorted((pair,)+rest))
            assert (G[index[p],:]*W[:,j])[0] == 0
            checks += 1
assert checks == 180

# Nontrivial positive control: (-1,-1,+1,+1).
def eps(i,j): return 0 if i == j else (1 if (i,j) == (0,1) else -1)
P2 = [((0,1),(2,3)),((0,2),(1,3))]
V = s.Matrix([[s.prod(eps(t[a],t[b]) for a,b in p) for p in P2]
              for t in product(range(2),repeat=4)])
G2 = V.T*V
W2 = s.Matrix([[0,-1,-1,0],[1,1,1,1]])
D1 = W2.col_join(s.zeros(1,4))
D2 = s.zeros(4,3)
for j,sg in enumerate((-1,-1,1,1)):
    D2[j,1:3] = sg*(W2[:,j].T*G2)/4
assert chain_dims((1,)*4,{0,1}) == [3,4,3]
assert D1*D2 == s.zeros(3)
assert (D1.rank(),D2.rank()) == (2,2)
assert [3-D1.rank(),4-D1.rank()-D2.rank(),3-D2.rank()] == [1,0,1]
print('PASS: grading obstruction, square-zero signs, exact differentials, homology')
print('adjacent-degree witness:',small,'Phi =',euler(small))
print('grading witness with odd-block maps:',small_h)
print('six adjoints:',c,'homology =',h,'Phi =',euler(c))
print('cut Gram values:',sorted(set(K)),'rank =',K.rank())
print('mixed fundamental homology: [1, 0, 1]')
