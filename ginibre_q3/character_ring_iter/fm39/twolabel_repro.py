
from collections import defaultdict
from functools import lru_cache
from itertools import product
from math import comb
import argparse

def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def cat_moment(n):
    return 0 if n % 2 else comb(n, n//2)//(n//2+1)

@lru_cache(None)
def U(n):
    return {(n-2*h, 0): (-1)**h*comb(n-h, h)
            for h in range(n//2+1)}

def swap(P):
    return {(j, i): v for (i, j), v in P.items()}

def add(P, Q):
    R = defaultdict(int, P)
    for ij, v in Q.items():
        R[ij] += v
    return {ij: v for ij, v in R.items() if v}

def mul(P, Q):
    R = defaultdict(int)
    for (i,j), v in P.items():
        for (k,l), w in Q.items():
            R[i+k,j+l] += v*w
    return {ij: v for ij, v in R.items() if v}

@lru_cache(None)
def h(n):
    R = {}
    for j in range(n+1):
        R = add(R, mul(U(j), swap(U(n-j))))
    return R

@lru_cache(None)
def S(n):
    return add(U(n), swap(U(n)))

@lru_cache(None)
def weight(a, e):
    # Direct expansion of (x+y)^a (x-y)^e in x and y.
    W = defaultdict(int)
    for j in range(a+1):
        for k in range(e+1):
            W[j+k] += comb(a,j)*comb(e,k)*(-1)**(e-k)
    return tuple(W.items())

@lru_cache(None)
def moment(a, e, u, v):
    n = a+e
    return sum(w*cat_moment(k+u)*cat_moment(n-k+v)
               for k, w in weight(a,e))

def integral(P, a, e):
    return sum(v*moment(a,e,i,j) for (i,j),v in P.items())

def direct_phi(P, a, r):
    z = integral(P,a,2*r)
    assert z % 2 == 0
    return z//2

@lru_cache(None)
def coeff(a,e):
    # Independent multiplication by linear factors.
    c = [1]
    for sign, count in ((1,a),(-1,e)):
        for _ in range(count):
            d = [0]*(len(c)+1)
            for k, x in enumerate(c):
                d[k] += x
                d[k+1] += sign*x
            c = d
    return tuple(c)

def data(a,e):
    c = coeff(a,e)
    n = a+e
    def C(k): return c[k] if 0 <= k <= n else 0
    def B(k): return C(k-1)+C(k+1)
    def D(k): return C(k)**2-C(k-1)*C(k+1)
    def W(i,j): return B(i)*C(j)-C(i)*B(j)
    return n,C,B,D,W

def endpoints(a,e,p,q):
    if (a+e+p+q) % 2:
        return 0,0
    p,q = max(p,q),min(p,q)
    n,C,B,D,W = data(a,e)
    j = (n+p-q)//2
    i = (n+p+q)//2+1
    return D(j)-D(i),W(i,j)

def predicted(a,r,s,t,kind):
    if kind == "HH":
        T,W = endpoints(a,2*r-2,s+1,t+1)
        return T-W
    if kind == "SS":
        T,W = endpoints(a,2*r,s,t)
        return T+W
    p,q = s+1,t
    T,W = endpoints(a,2*r-1,p,q)
    return T+(W if p >= q else -W)

def main():
    argparse.ArgumentParser(
        description="Exact two-label identity and strip checks."
    ).parse_args()

    kernel_checks = 0
    for a,e,p,q in product(range(7),range(7),range(6),range(6)):
        direct = integral(mul(U(p),swap(U(q))),a,e)
        T,W = endpoints(a,e,p,q)
        expected = W if p >= q or e % 2 == 0 else -W
        assert direct == expected, (a,e,p,q,direct,expected)
        if p >= q:
            row_sum = sum(integral(U(k),a,e)
                          for k in range(p-q,p+q+1,2))
            assert row_sum == T
        kernel_checks += 1

    totals = {"HH":0,"SS":0,"HS":0}
    cases = list(product(range(1,6),range(6),range(5),range(5)))
    for r in (9,17,33):
        for a in (0,1,2,2*r-3,2*r-2,2*r-1):
            for s,t in ((2,3),(3,5),(4,4),(5,7)):
                cases.append((r,a,s,t))

    for r,a,s,t in cases:
        for kind,P in (("HH",mul(h(s),h(t))),
                       ("SS",mul(S(s),S(t))),
                       ("HS",mul(h(s),S(t)))):
            direct = direct_phi(P,a,r)
            expected = predicted(a,r,s,t,kind)
            assert direct == expected, (r,a,s,t,kind,direct,expected)
            assert direct >= 0
            totals[kind] += 1

    endpoint_checks = 0
    for a,e in product(range(13),repeat=2):
        n,C,B,D,W = data(a,e)
        _,Cs,Bs,Ds,Ws = data(e,a)
        for j in range((n+1)//2,n+2):
            for i in range(j+1,n+3):
                T = D(j)-D(i)
                w = W(i,j)
                assert T >= abs(w), (a,e,i,j,T,w)
                assert Ds(j)-Ds(i) == T
                assert Ws(i,j) == (-1)**(i+j-1)*w
                for sign in (-1,1):
                    determinant = (
                        (C(j)+sign*C(i-1))*(C(j)+sign*C(i+1))
                        -(C(i)+sign*C(j-1))*(C(i)+sign*C(j+1))
                    )
                    assert determinant == T+sign*w
                endpoint_checks += 1

    balanced_checks = 0
    for e in range(13):
        n,C,B,D,W = data(e,e)
        b = lambda k: choose(e,k)
        for j in range(e,n+2):
            for i in range(j+1,n+3):
                slack = D(j)-D(i)-abs(W(i,j))
                if (i-j) % 2 == 0:
                    assert W(i,j) == 0 and slack >= 0
                elif j % 2 == 0:
                    m,z = j//2,(i-1)//2
                    A,B0,C0 = b(m),b(z),b(z+1)
                    assert A >= B0 >= C0 >= 0
                    assert slack == (A-B0)*(A+C0)
                else:
                    m,z = (j-1)//2,i//2
                    A,B0,C0 = b(m),b(m+1),b(z)
                    assert A >= B0 >= C0 >= 0
                    assert slack == (B0-C0)*(A+C0)
                balanced_checks += 1

    adjacent_checks = 0
    for e in range(13):
        n,C,B,D,W = data(e+1,e)
        b = lambda k: choose(e,k)
        def pair(k):
            m = k//2
            return b(m),b(m-1 if k%2 == 0 else m+1)
        for j in range((n+1)//2,n+2):
            for i in range(j+1,n+3):
                x,A = pair(j)
                y,E = pair(i)
                sigma = (-1)**(i//2+j//2)
                T,w = D(j)-D(i),W(i,j)
                assert x >= y >= 0 and x+A >= y+E
                assert w == sigma*(y*A-x*E)
                assert T-sigma*w == (x-y)*(x+y+A+E)
                assert T+sigma*w == (x+y)*(x-y+A-E)
                adjacent_checks += 1

    from fractions import Fraction
    from math import factorial

    def K(e,x,N):
        if e == 0:
            return 1
        previous,current = 1,x
        for k in range(1,e):
            previous,current = current,x*current-k*(N-k+1)*previous
        return current

    kraw_points = kraw_pairs = 0
    for N in range(19):
        for e in range(N+1):
            n,C,B,D,W = data(N-e,e)
            falling = factorial(N)//factorial(N-e)
            vectors = {}
            for k in range((N+1)//2,N+1):
                x = 2*k-N
                f = K(e,x,N)
                lower = Fraction(k,N-k+1)
                upper = Fraction(N-k,k+1)
                H = lower*K(e,x-2,N)+upper*K(e,x+2,N)
                G = f*f-lower*upper*K(e,x-2,N)*K(e,x+2,N)
                scale = Fraction(comb(N,k),falling)
                assert C(k) == (-1)**e*scale*f
                assert B(k) == (-1)**e*scale*H
                assert D(k) == scale*scale*G
                vectors[k] = f,H,G,scale
                kraw_points += 1
            for j,(f,H,G,scale) in vectors.items():
                for i in range(j+1,N+1):
                    F,L,J,_ = vectors[i]
                    rho = Fraction(comb(N,i),comb(N,j))
                    value = G-rho*rho*J-rho*(L*f-F*H)
                    assert scale*scale*value == D(j)-D(i)-W(i,j)
                    kraw_pairs += 1

    witness_values = [direct_phi(mul(h(s),h(t)),8,4)
                      for s,t in ((6,2),(7,3))]
    assert witness_values == [1064,888]

    print("Catalan kernel identities:",kernel_checks)
    print("Direct Catalan word identities:",totals)
    print("Endpoint, folding and determinant checks:",endpoint_checks)
    print("Balanced factorization checks:",balanced_checks)
    print("Adjacent-balance factorization checks:",adjacent_checks)
    print("Krawtchouk vector / bilinear checks:",kraw_points,kraw_pairs)
    print("Endpoint-monotonicity counterexample:",witness_values)
    print("All assertions passed.")

if __name__ == "__main__":
    main()
