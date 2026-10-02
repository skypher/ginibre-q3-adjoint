cd /home/yang/q3adjoint
python3 -u - <<'PY'
from collections import defaultdict, Counter
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import comb
from fractions import Fraction as Q
from datetime import datetime, timezone

def stamp(s):
    print(datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
          s, flush=True)

def cg(a, n):
    return range(abs(a-n), a+n+1, 2)

def dcoeff(fs):
    d={(0,0):1}
    for eps,n in fs:
        q=defaultdict(int)
        for (r,s),v in d.items():
            for rr in cg(r,n):
                q[rr,s]+=v
            for ss in cg(s,n):
                q[r,ss]+=eps*v
        d={k:v for k,v in q.items() if v}
    deg=max((r for r,s in d if r==s),default=-1)
    return [d.get((j,j),0) for j in range(deg+1)]

def taylor(c, endpoint):
    out=[]
    for k in range(len(c)):
        if endpoint==1:
            out.append(sum(c[j]*comb(j,k)*(-1)**k
                           for j in range(k,len(c))))
        else:
            out.append(sum(c[j]*comb(j,k)*(-1)**(j-k)
                           for j in range(k,len(c))))
    return out

def trim(c):
    while c and c[-1]==0:
        c.pop()
    return c

def deriv(c,x):
    return sum(j*a*x**(j-1) for j,a in enumerate(c) if j)

def quotient_one_minus(c):
    q=[];s=0
    for a in c:
        s+=a
        q.append(s)
    assert not q or q[-1]==0
    return trim(q[:-1])

def m0(ns):
    d={0:1}
    for n in ns:
        q=defaultdict(int)
        for j,v in d.items():
            for k in cg(j,n):
                q[k]+=v
        d=dict(q)
    return d.get(0,0)

def endpoint_values(fs):
    ns=[n for e,n in fs]
    plus=0 if any(e<0 for e,n in fs) else 2**len(fs)*m0(ns)
    reflected=[(e*((-1)**n),n) for e,n in fs]
    minus=0 if any(e<0 for e,n in reflected) else 2**len(fs)*m0(ns)
    return plus,minus

def screen_endpoints():
    tokens=tuple(e*n for n in range(1,4) for e in (1,-1))
    total=zero=0
    hplus=Counter()
    hminus=Counter()
    for N in range(2,9,2):
        for w in combinations_with_replacement(tokens,N):
            if sum(z<0 for z in w)%2:
                continue
            fs=[(1 if z>0 else -1,abs(z)) for z in w]
            c=dcoeff(fs)
            a=taylor(c,1)
            b=taylor(c,-1)
            ep,em=endpoint_values(fs)
            assert sum(c)==ep
            assert sum(c[j]*(-1)**j for j in range(len(c)))==em
            cr=dcoeff([(e*((-1)**n),n) for e,n in fs])
            assert trim(cr)==trim([(-1)**j*v for j,v in enumerate(c)])
            total+=1
            if not any(c):
                zero+=1
            if sum(e<0 for e,n in fs):
                k=next((i for i,v in enumerate(a) if v),None)
                if k is None:
                    hplus['zero']+=1
                else:
                    hplus[(sum(e<0 for e,n in fs),k)]+=1
                    assert a[k]>0
            km=next((i for i,v in enumerate(b) if v),None)
            if km is None:
                assert not any(c)
            else:
                hminus[km]+=1
                assert b[km]>0
    assert (total,zero)==(965,450)
    assert hplus==Counter({
        (2,1):158, 'zero':409, (2,2):1,
        (4,1):173, (4,3):1, (6,1):104, (8,1):25
    })
    assert hminus==Counter({1:460,0:53,3:1,2:1})
    return total,zero,hplus,hminus

def local_delta(a):
    k=next(i for i,v in enumerate(a) if v)
    lead=a[k]
    rem=sum(abs(v) for v in a[k+1:])
    assert lead>0
    delta=Q(1) if rem==0 else min(Q(1),Q(lead,2*rem))
    assert delta*rem<=Q(lead,2)
    return k,lead,rem,delta

def free_trace_poly(fs):
    # Noncrossing Wick pairings; pairs within one Wick block are excluded.
    slots=[]
    for i,(eps,n) in enumerate(fs):
        slots.extend([i]*n)
    if len(slots)%2:
        return [0]
    total=defaultdict(int)
    for colors in product((0,1),repeat=len(fs)):
        sign=1
        for (eps,n),color in zip(fs,colors):
            if color:
                sign*=eps
        @lru_cache(None)
        def pairings(l,r):
            if l>r:
                return ((0,1),)
            out=defaultdict(int)
            for j in range(l+1,r+1,2):
                if slots[l]==slots[j]:
                    continue
                left=pairings(l+1,j-1)
                right=pairings(j+1,r)
                shift=int(colors[slots[l]]!=colors[slots[j]])
                for a,ca in left:
                    for b,cb in right:
                        out[a+b+shift]+=ca*cb
            return tuple(sorted(out.items()))
        for degree,count in pairings(0,len(slots)-1):
            total[degree]+=sign*count
    return trim([total.get(i,0) for i in range(max(total,default=-1)+1)])

stamp('endpoint expansion checks begin')
res=screen_endpoints()
stamp(f'endpoint census pass: {res[0]} lists; {res[1]} have D identically zero')
stamp('plus-endpoint orders: '+str(dict(res[2])))
stamp('minus-endpoint orders: '+str(dict(res[3])))

examples=[
    ('four minus ones',[(-1,1)]*4),
    ('three minus ones and minus three',[(-1,1)]*3+[(-1,3)]),
    ('two minus ones',[(-1,1)]*2),
    ('derivative witness',[(1,1),(-1,2),(-1,3)]),
]
for name,fs in examples:
    c=dcoeff(fs)
    a=taylor(c,1)
    b=taylor(c,-1)
    print(name,'Dcoeff=',c,'at1=',a,'atminus1=',b,flush=True)
    if name=='four minus ones':
        assert c==[10,-16,6]
        assert a==[0,4,6]
        assert local_delta(a)==(1,4,6,Q(1,3))
    if name=='three minus ones and minus three':
        assert c==[2,-6,6,-2]
        assert a==[0,0,0,2]
        assert local_delta(a)==(3,2,0,Q(1))
    if name=='two minus ones':
        assert c==[2,-2]
    if name=='derivative witness':
        assert c==[2,2,-2,-2]
        assert deriv(c,Q(0))==2 and deriv(c,Q(1))==-8
        q=quotient_one_minus(c)
        assert q==[2,4,2] and deriv(q,Q(0))==4
q=quotient_one_minus([10,-16,6])
assert q==[10,-6] and deriv(q,Q(0))==-6
stamp('endpoint Taylor and explicit delta examples pass')

for rho in (Q(0),Q(1,2),Q(3,4)):
    x4=Q(2)+(1-rho*rho)/2
    assert x4>2
assert Q(2)+(1-Q(0)**2)/2==Q(5,2)
stamp('classical semicircle linear-rotation fourth-moment obstruction passes')

fs=[(1,1),(1,2),(-1,1),(-1,2)]
actual=dcoeff(fs)
free=free_trace_poly(fs)
assert actual==[2,0,-2] and free==[4,0,-4]
stamp(f'free-rotation mismatch: Lancaster={actual}, free-Wick={free}')
stamp('ALL EXACT CHECKS PASS')
PY
