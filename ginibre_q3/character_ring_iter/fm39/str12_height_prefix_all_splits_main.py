# FM-STR12 (main agent): height-prefix positivity for EVERY split.
# For a signed list split as A | B, Pi_T(A,B) = sum over channels (r,s) with r+s <= T of f_A(r,s) f_B(r,s),
# where f_X(r,s) is the coefficient of U_r(x)U_s(y) in prod_{z in X} (U_|z|(x) + sign(z) U_|z|(y)).
# Pi_infinity(A,B) = Phi(A u B).  This script checks Pi_T >= 0 for every T and every split of every list
# with even minus count in a label/length box, optionally pair-free only, and tests other truncation shapes.
# Usage: python3 -u str12_height_prefix_all_splits_main.py MAXLABEL MAXLEN [pairfree|all] [shapes]
import itertools, sys, time
from collections import defaultdict
def cg(a,b): return range(abs(a-b),a+b+1,2)
def table(word):
    d={(0,0):1}
    for z in sorted(word,key=abs,reverse=True):
        n=abs(z); e=1 if z>0 else -1; q=defaultdict(int)
        for (a,b),v in d.items():
            for c in cg(a,n): q[c,b]+=v
            for c in cg(b,n): q[a,c]+=e*v
        d={k:v for k,v in q.items() if v}
    return d
SHAPES={'height':lambda r,s:r+s,'r+2s':lambda r,s:r+2*s}
def bad_shape(prod,key):
    lay=defaultdict(int)
    for (r,s),v in prod.items(): lay[key(r,s)]+=v
    run=0
    for t in sorted(lay):
        run+=lay[t]
        if run<0: return t
    return None
if __name__=='__main__':
    Mx=int(sys.argv[1]); Nx=int(sys.argv[2]); pf=(len(sys.argv)<4 or sys.argv[3]=='pairfree')
    shapes=SHAPES if (len(sys.argv)>4 and sys.argv[4]=='shapes') else {'height':SHAPES['height']}
    labs=[s*n for n in range(1,Mx+1) for s in (1,-1)]
    st=defaultdict(int)
    for N in range(2,Nx+1):
        for w in itertools.combinations_with_replacement(sorted(labs),N):
            if sum(1 for z in w if z<0)%2 or sum(map(abs,w))%2: continue
            if pf and any(-z in w for z in w): continue
            st['lists']+=1
            for m in range(1,1<<(N-1)):
                A=[w[i] for i in range(N) if m>>i&1]; Bk=[w[i] for i in range(N) if not m>>i&1]
                fa=table(A); fb=table(Bk); prod={k:v*fb.get(k,0) for k,v in fa.items() if fb.get(k,0)}
                st['splits']+=1
                for name,key in shapes.items():
                    t=bad_shape(prod,key)
                    if t is not None:
                        st['bad_'+name]+=1; print('FAIL',name,w,A,Bk,'T=',t,flush=True)
        print(time.strftime('%H:%M:%S'),'length',N,dict(st),flush=True)
    assert not any(k.startswith('bad_') for k in st), 'a negative prefix was found'
    print('FM-STR12 HEIGHT PREFIX ALL SPLITS PASS')
