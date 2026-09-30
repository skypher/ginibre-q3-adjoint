
import argparse
from math import comb
from fractions import Fraction as Q
from functools import lru_cache
from collections import defaultdict
from itertools import combinations_with_replacement as cwr

argparse.ArgumentParser(
    description="Exact verifier for the large-label three-factor theorem."
).parse_args()

def mul(f,g):
    out=defaultdict(int)
    for (i,j),c in f.items():
        for (k,l),d in g.items():
            out[i+k,j+l]+=c*d
    return {ij:c for ij,c in out.items() if c}

@lru_cache(None)
def H(k):
    out={}
    for b in range((k+1)//2+1):
        d=k+1-2*b
        for j in range(d):
            out[d-1-j,j]=(-1)**b*comb(k+1-b,b)
    return out

@lru_cache(None)
def cat(k):
    return 0 if k%2 else comb(k,k//2)//(k//2+1)

@lru_cache(None)
def paired_moment(r,i,j):
    return sum(
        (-1)**b*comb(2*r,b)*cat(i+2*r-b)*cat(j+b)
        for b in range(2*r+1)
    )

@lru_cache(None)
def phi(r,a,parts):
    f={(a-j,j):comb(a,j) for j in range(a+1)}
    for k in parts:
        f=mul(f,H(k))
    return Q(sum(
        c*paired_moment(r,i,j) for (i,j),c in f.items()
    ),2)

def cg(p,q):
    return range(abs(p-q),p+q+1,2)

def main_term(r,a,parts):
    N=a+2*r-3
    A,B,C=[k+1 for k in parts]
    mu=defaultdict(int)
    for d in cg(A,B):
        for l in cg(d,C):
            if 1<=l<=N:
                mu[l]+=1
    return sum(
        mult*phi(r-1,a,(l-1,))
        for l,mult in mu.items()
    )

@lru_cache(None)
def monomial_U(i,p):
    return sum(
        (-1)**b*comb(p-b,b)*cat(i+p-2*b)
        for b in range(p//2+1)
    )

def cross_kernel(e,a,p,q):
    c=[0]*(e+a+1)
    for i in range(e+1):
        for j in range(a+1):
            c[i+j]+=(-1)**i*comb(e,i)*comb(a,j)
    N=e+a
    return sum(
        c[j]*monomial_U(N-j,p)*monomial_U(j,q)
        for j in range(N+1)
    )

large_checks=gap_checks=support_checks=boundary_cases=0
edge_seen=set()

for r in range(2,6):
    for a in range(6):
        N=a+2*r-3

        for l in range(1,N+1):
            assert phi(r-1,a,(l-1,))>=0

        for q in range(N,N+3):
            for p in range(N+3):
                expected=-1 if q==N and p==0 else 0
                assert cross_kernel(2*r-3,a,p,q)==expected
                support_checks+=1

        lower=max(1,N-1)
        for parts in cwr(range(lower,lower+4),3):
            A,B,C=[k+1 for k in parts]
            edge=sum((
                C==N and A==B,
                B==N and A==C,
                A==N and B==C
            ))
            value=phi(r,a,parts)
            main=main_term(r,a,parts)
            assert value==main+edge
            assert value>=0
            assert (a+sum(parts))%2==0 or value==0
            large_checks+=1
            boundary_cases+=edge>0
            edge_seen.add(edge)

        triples=set()
        for w in range(1,4):
            for gap in range(2):
                v=w+gap
                for extra in range(2):
                    u=v+max(0,N-w)+extra
                    triples.add((u,v,w))

        for parts in sorted(triples):
            u,v,w=parts
            assert u>=v>=w and u-v+w>=N
            assert phi(r,a,parts)==main_term(r,a,parts)>=0
            gap_checks+=1

print("Large-label exact identities:",large_checks)
print("Support-gap exact identities:",gap_checks)
print("Kernel support/boundary checks:",support_checks)
print("Cases with a positive boundary correction:",boundary_cases)
print("Boundary correction values:",sorted(edge_seen))

for r,a in ((2,2),(3,1),(4,3)):
    N=a+2*r-3
    parts=(N-1,)*3
    print("Boundary diagonal:",(r,a,parts),
          "phi =",phi(r,a,parts),
          "main =",main_term(r,a,parts),
          "correction = 3")

print("All exact checks passed.")
