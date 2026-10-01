import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial as fac
from random import Random

ap=argparse.ArgumentParser(description="Exact full-cone midpoint checks")
ap.add_argument("--part",choices=("small","large","kernels"),default="small")
ap.add_argument("--block",type=int,choices=range(4),default=0)
args=ap.parse_args()

def mul(c,n,sg):
    out=Counter()
    for (r,s),v in c.items():
        for t in range(abs(r-n),r+n+1,2): out[t,s]+=v
        for t in range(abs(s-n),s+n+1,2): out[r,t]+=sg*v
    return {k:v for k,v in out.items() if v}

def gram(c):
    a=b=c2=Q(0)
    for (r,s),v in c.items():
        if r==s:
            w=Q((-1)**r,r+1)
            a+=v*w
            if r:
                b-=v*w
                c2+=v*Q((-1)**r*(r*r+2*r-4),r*(r+1)*(r+2))
        elif abs(r-s)==2:
            t=min(r,s)
            c2+=v*Q(2*(-1)**t,(t+1)*(t+2)*(t+3))
    return a,b,c2

stats=Counter()
extreme={}
def check(word,c):
    a,b,c2=gram(c); E=c.get((0,0),0)
    assert a+b==E and E>=0
    assert a>=0 and c2>=0 and a*c2>=b*b
    stats["tested"]+=1
    stats["even_minus"]+=sum(sg<0 for _,sg in word)%2==0
    stats["zero_E"]+=E==0
    stats["zero_A"]+=a==0
    stats["upper_fails"]+=b>a
    stats["norm_fails"]+=c2>a
    stats["negative_B_norm_fails"]+=b<0 and c2>a
    if a:
        rec=(b/a,tuple(Counter(word).items()),(a,b,c2),E)
        if "min" not in extreme or rec[0]<extreme["min"][0]:
            extreme["min"]=rec
        if "max" not in extreme or rec[0]>extreme["max"][0]:
            extreme["max"]=rec

def wordcheck(word):
    c={(0,0):1}
    for n,sg in word: c=mul(c,n,sg)
    check(word,c)
    return gram(c)

if args.part=="small":
    types=[(n,sg) for n in range(1,5) for sg in (1,-1)]
    def visit(word,c,lo):
        check(word,c)
        if len(word)<8:
            for j in range(lo,8):
                n,sg=types[j]
                visit(word+[(n,sg)],mul(c,n,sg),j)
    visit([],{(0,0):1},0)
    assert stats["tested"]==12870 and stats["even_minus"]==6470
    assert stats["zero_E"]==stats["zero_A"]==9609
    assert stats["upper_fails"]==804 and stats["norm_fails"]==1260
    assert stats["negative_B_norm_fails"]==13
    assert extreme["min"][0]==Q(-37614314,41563259)
    assert extreme["max"][0]==4

if args.part=="large":
    rng=Random(161004)
    for case in range(128):
        if case<64:
            while True:
                ns=[rng.randrange(5,81)
                    for _ in range(rng.randrange(4,9))]
                ns[-1]+=sum(ns)%2
                if sum(ns)-2*max(ns)>=6: break
            ss=[rng.choice((1,-1)) for _ in ns[:-1]]
            ss+=[(-1)**ss.count(-1)]
        elif case<96:
            while True:
                e,a,b=(rng.randrange(2,41),rng.randrange(2,41),
                       rng.randrange(3,10))
                d=rng.randrange(6,21); p=e+a+2*b-2*d
                if p>=7: break
            ns=[1]*(e+a)+[2]*b+[p]
            ss=[-1]*e+[1]*(a+b)+[(-1)**e]
        else:
            a,e=rng.randrange(2,13),rng.randrange(2,13)
            p=rng.randrange(20,101)
            dif=rng.randrange(-(a+e-6)//2,
                              (a+e-6)//2+1)*2+(a+e)%2
            ns=[1]*(a+e)+[p,p+dif]
            ss=[1]*a+[-1]*e+[rng.choice((1,-1))]
            ss+=[(-1)**ss.count(-1)]
        if case//32==args.block:
            wordcheck(list(zip(ns,ss)))
    assert stats["tested"]==32 and stats["zero_E"]==0

if args.part=="kernels":
    @lru_cache(None)
    def beta(m,g):
        return Q(2*fac(2*m)*fac(2*g+2),
                 4**(m+g+1)*fac(m)*fac(g+1)*fac(m+g+1))
    def radial(m,g,j):
        return beta(m,g) if j==0 else 4*beta(m+1,g)-beta(m,g)
    @lru_cache(None)
    def moment(h,k,j1,j2):
        if (h+k)%2: return Q(0)
        ans=Q(0)
        for i in range(h+1):
            for j in range(k+1):
                if (i+j)%2: continue
                g=(i+j)//2; m=(h+k)//2-g
                ans+=Q(2**(h+k)*comb(h,i)*comb(k,j)*(-1)**i,
                       2*g+1)*radial(m,g,j1)*radial(m,g,j2)
        return ans
    def direct(r,s,j1,j2):
        return sum(
            (-1)**(i+j)*comb(r-i,i)*comb(s-j,j)
            *moment(r-2*i,s-2*j,j1,j2)
            for i in range(r//2+1) for j in range(s//2+1))
    for r in range(9):
        for s in range(9):
            a,b,c2=gram({(r,s):1})
            for j1,j2,v in ((0,0,a),(0,2,b),(2,2,c2)):
                assert direct(r,s,j1,j2)==v
    for n in range(1,13):
        assert wordcheck([(1,1)]*n+[(n,1)])==(
            Q(2,n+1),Q(2*n,n+1),Q(2*n*n,n+1))
    for m in range(1,9):
        a,b,c2=wordcheck([(1,-1)]*(2*m))
        assert b/a==Q(-m,m+2)
        assert c2/a==Q(m*m,(m+2)**2)
    assert wordcheck([(2,1)]*3)==(4,-2,Q(26,5))
    assert wordcheck([(2,1)])==(0,0,Q(2,3))
    assert Q(-2)*(Q(-2)+Q(26,5))/64==-Q(1,10)

    import sympy as S
    t=S.symbols("t")
    f=Q(48,5)*t**6-Q(56,5)*t**4+Q(28,5)*t**2-1
    def avg(f):
        return sum(v*Q(comb(i,i//2),4**(i//2)*(i//2+1))
                   for (i,),v in S.Poly(f,t).terms() if i%2==0)
    assert avg(f)==-Q(1,4)
    assert avg((4*t*t-1)*f)==Q(13,20)
    for k in range(1,7):
        N=2**k-1
        rho=S.chebyshevu(N,t)**2
        assert S.Poly(rho-sum(S.chebyshevu(2*j,t)
                             for j in range(N+1)),t).is_zero
        assert S.Poly(rho-S.prod((2*S.chebyshevt(2**j,t))**2
                                for j in range(k)),t).is_zero
    print("243 quaternion checks; two exact families; six tilt depths")

print(args.part,dict(stats))
for key,value in extreme.items(): print(key,value)
