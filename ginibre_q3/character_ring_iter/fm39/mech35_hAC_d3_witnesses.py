"""FM-MECH35 (astra_max_ceres), extracted from its run log: mech35_hAC_d3_witnesses.py"""
import argparse
from fractions import Fraction as Q
from math import comb
from itertools import combinations
argparse.ArgumentParser(description='Exact certificates for the boundary and obstruction witnesses').parse_args()
def fusion(mu,n):
    rows=[{0:1}]
    for x in mu:
        for row in rows[:]:
            out={}
            for j,a in row.items():
                for h in range(abs(j-x),j+x+1,2): out[h]=out.get(h,0)+a
            rows.append(out)
    full=len(rows)-1
    return [rows[S].get(0,0)*rows[S^full].get(n,0) for S in range(full+1)]
def check(t,v,d,terms):
    mu=(1,)*t+(2,)*v; f=fusion(mu,sum(mu)-2*d); mask=(1<<t)-1
    result=[Q(0)]*len(f)
    for a,H in terms:
        assert a>=0 and all(x^y in H for x in H for y in H)
        counts={}
        for x in H:
            key=((x&mask).bit_count(),(x>>t).bit_count())
            counts[key]=counts.get(key,0)+1
        for x in range(len(f)):
            i,j=(x&mask).bit_count(),(x>>t).bit_count()
            result[x]+=a*Q(counts.get((i,j),0),comb(t,i)*comb(v,j))
    assert result==f
    print('CP certificate:',mu,'n',sum(mu)-2*d,'terms',len(terms),'entries',len(f),'PASS')
check(4,2,3,[
 (Q(6),(0,3,5,6)),
 (Q(1),(0,3,5,6,25,26,28,31,41,42,44,47,48,51,53,54)),
 (Q(2),(0,3,21,22,44,47,57,58)),
 (Q(1),(0,3,28,31,44,47,48,51)),
 (Q(2),(0,15,19,28,37,42,54,57)),
 (Q(1),(0,19,35,48))])
check(6,2,3,[
 (Q(45,16),(0,3,5,6,9,10,12,15,48,51,53,54,57,58,60,63)),
 (Q(35,16),(0,3,5,6,9,10,12,15,48,51,53,54,57,58,60,63,81,82,84,87,88,91,93,94,97,98,100,103,104,107,109,110)),
 (Q(105,8),(0,3,5,6,9,10,12,15,81,82,84,87,88,91,93,94)),
 (Q(5,8),(0,3,5,6,24,27,29,30,40,43,45,46,48,51,53,54)),
 (Q(5),(0,3,5,6,192,195,197,198)),
 (Q(5),(0,3,12,15,53,54,57,58,81,82,93,94,100,103,104,107)),
 (Q(15,2),(0,3,12,15,69,70,73,74)),
 (Q(15,4),(0,3,12,15,197,198,201,202))])
check(0,5,4,[
 (Q(50,9),(0,)),(Q(25,9),(0,3,5,6)),
 (Q(25,6),(0,3,5,6,24,27,29,30)),
 (Q(5,2),(0,3,12,15,21,22,25,26))])
def square(p):
    ans=[Q(0)]*16
    for x,a in p.items():
        for y,b in p.items(): ans[x^y]+=a*b
    return ans
base=[Q(0)]*16
p={1<<i:Q(1) for i in range(4)}
for k,a in enumerate(square(p)): base[k]+=a/2
for i,j,k in combinations(range(4),3):
    p={1<<i:Q(1),(1<<j)^(1<<k):Q(1)}
    for z,a in enumerate(square(p)): base[z]+=a/2
for q in range(2,22,2):
    s=q//2; origin=1+5*s*(s+1)//2
    f=fusion((q,)*4,q); cert=base[:]; cert[0]+=origin-6
    assert cert==f
print('repeated even labels (q^5), q=2,4,...,20: PASS')
ball=[x for x in range(16) if 2*x.bit_count()<=3]
assert 7 not in {x^y for x in ball for y in ball}
print('radius-d whole-cluster factors fail at (2^4;n=2), d=3; H_AC holds')
det=Q(1,4)*Q(1,2)-Q(3,8)**2
assert det==Q(-1,64)
print('uniform weight-class Gram determinant:',det,'; H_AC witness certified above')