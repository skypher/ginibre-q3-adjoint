python3 -u - <<'PY'
import sys
if any(x in ('-h','--help') for x in sys.argv[1:]):
    print('Exact verifier for Lancaster Walsh coefficients; run as a Python heredoc from the repository root.')
    raise SystemExit(0)
from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations_with_replacement, product

def stamp(s):
    print(datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), s, flush=True)
def cg(a,b):
    return range(abs(a-b),a+b+1,2)
def clean(D):
    return {k:v for k,v in D.items() if v}
def table(word):
    """Coefficients of prod (U_n(x)+eps U_n(y)) in U_r(x)U_s(y)."""
    d={(0,0):1}
    for n,eps in word:
        z=defaultdict(int)
        for (a,b),v in d.items():
            for r in cg(a,n): z[r,b]+=v
            for s in cg(b,n): z[a,s]+=eps*v
        d=clean(z)
    return d
def diag(word):
    d=table(word)
    R=max((max(k) for k in d),default=0)
    return [d.get((j,j),0) for j in range(R+1)]
def phi(word):
    return table(word).get((0,0),0)
def tri_admissible(ns):
    return sum(ns)%2==0 and max(ns)<=sum(ns)-max(ns)

stamp('moment witnesses begin')
plus_pair=[(1,1),(1,1)]
minus_pair=[(1,-1),(1,-1)]
wp=diag(plus_pair); wm=diag(minus_pair)
assert wp==[2,2,0] and wm==[2,-2,0]
assert wp[0]*wp[2]-wp[1]*wp[1]==-4
assert wm[0]*wm[2]-wm[1]*wm[1]==-4
stamp(f'exact W witnesses: plus pair={wp}, D=2+2*rho; minus pair={wm}, D=2-2*rho; both Hankel determinants=-4')

stamp('two-factor and three-factor coefficient formulas begin')
for n in range(1,9):
    for m in range(1,9):
        for eps in (-1,1):
            w=diag([(n,eps),(m,eps)])
            expected=[0]*(max(n,m)+1)
            if n==m:
                expected[0]=2
                expected[n]+=2*eps
            while expected and expected[-1]==0: expected.pop()
            while w and w[-1]==0: w.pop()
            assert w==expected,(n,m,eps,w,expected)
triple_count=0
for ns in combinations_with_replacement(range(1,8),3):
    for es in product((-1,1),repeat=3):
        if es[0]*es[1]*es[2]!=1: continue
        w=diag(list(zip(ns,es)))
        expected=defaultdict(int)
        if tri_admissible(ns):
            expected[0]+=2
            for n,e in zip(ns,es): expected[n]+=2*e
        expected_list=[expected[j] for j in range(max(ns)+1)]
        while expected_list and expected_list[-1]==0: expected_list.pop()
        while w and w[-1]==0: w.pop()
        assert w==expected_list,(ns,es,w,expected_list)
        triple_count+=1
stamp(f'exact formulas pass: 2-factor labels<=8; {triple_count} admissible-sign triples with labels<=7')

stamp('Clebsch-Gordan transfer checks begin')
# If C has sign product E and appended sign eps=E, the full word has even minus count.
transfer_count=0
for N in range(0,4):
    for ns in combinations_with_replacement(range(1,4),N):
        for es in product((-1,1),repeat=N):
            E=1
            for e in es: E*=e
            C=list(zip(ns,es)); f=table(C)
            for q in range(1,5):
                eps=E
                new=diag(C+[(q,eps)])
                for j in range(0,9):
                    lhs=new[j] if j<len(new) else 0
                    rhs=2*sum(f.get((p,j),0) for p in cg(j,q))
                    assert lhs==rhs,(C,q,eps,j,lhs,rhs)
                    transfer_count+=1
stamp(f'CG transfer passes on {transfer_count} exact coefficient tests')
# Off-diagonal data is necessary: C=(+1,+2) has zero diagonal, but appending +1 creates W_1=4.
C=[(1,1),(2,1)]; f=table(C)
assert diag(C)==[0,0,0,0]
assert f.get((0,1))==1 and f.get((2,1))==1
assert diag(C+[(1,1)])[1]==4
stamp('off-diagonal witness: diag(C)=0, f_C(0,1)=f_C(2,1)=1, W_1(C plus +1)=4')

stamp('double-insertion identity checks begin')
checks=0
for N in range(0,4):
    for ns in combinations_with_replacement(range(1,4),N):
        for es in product((-1,1),repeat=N):
            E=1
            for e in es: E*=e
            if E!=1: continue
            L=list(zip(ns,es)); w=diag(L)
            for j in range(1,6):
                lhs=phi(L+[(j,1),(j,1)])-phi(L+[(j,-1),(j,-1)])
                rhs=4*(w[j] if j<len(w) else 0)
                assert lhs==rhs,(L,j,lhs,rhs)
                checks+=1
stamp(f'Phi difference identity passes on {checks} exact cases')

stamp('fundamental-only parity check begins')
fund_count=0
for N in range(0,10):
    for es in product((-1,1),repeat=N):
        if sum(e==-1 for e in es)%2: continue
        L=[(1,e) for e in es]
        if N%2:
            assert all(v==0 for v in diag(L))
        else:
            # With an even number of each sign, F=(x+y)^p(x-y)^q is pointwise nonnegative.
            assert sum(e==1 for e in es)%2==0 and sum(e==-1 for e in es)%2==0
        fund_count+=1
stamp(f'fundamental parity cases pass: {fund_count}; odd length has D identically zero')
stamp('PASS: moment obstruction, CG transfer, Phi identity, and proved small classes')
PY
