"""Exact CP certificate for (1^4,2^2,M,M+2); no files."""
import argparse,itertools as it
from fractions import Fraction as Q
argparse.ArgumentParser(description=__doc__).parse_args()

def inv(labels):
    if not labels: return 1
    p={0:1}
    for n in labels[:-1]:
        q={}
        for u,c in p.items():
            for j in range(abs(u-n),u+n+1,2):
                q[j]=q.get(j,0)+c
        p=q
    return p.get(labels[-1],0)

cert=[
 ((0,3,5,6),Q(7,6)), ((3,5,6,9,16),Q(2)),
 ((7,17,18,40),Q(5,6)), ((7,17,24,40),Q(1,2)),
 ((1,2,7,20,40),Q(1,2)),
 ((1,18,20,24,34),Q(7,6)),
 ((7,17,24,34,40),Q(1,3))
]
out=[Q(0)]*128
for pts,c in cert:
    for p in it.permutations(range(4)):
        for t in ((4,5),(5,4)):
            perm=p+t
            pp=[sum(1<<perm[i] for i in range(6) if x>>i&1)
                for x in pts]
            for x in pp:
                for y in pp:
                    out[x^y]+=c/48

for M in (4,7,31,1009):
    mu=(1,1,1,1,2,2,M)
    f=[inv(tuple(mu[i] for i in range(7) if s>>i&1)) *
       inv(tuple(mu[i] for i in range(7)
                 if not(s>>i&1))+(M+2,))
       for s in range(128)]
    assert out==f
    print("M",M,"exact entries",len(f),"origin",f[0])

vals=out[:]; h=1
while h<128:
    for a in range(0,128,2*h):
        for j in range(a,a+h):
            x,y=vals[j],vals[j+h]
            vals[j],vals[j+h]=x+y,x-y
    h*=2
print("minimum full EVEN value",2*min(vals))
# Minimum: 24.
