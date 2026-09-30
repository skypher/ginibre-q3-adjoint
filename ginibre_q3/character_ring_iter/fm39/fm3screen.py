import sys, time, json, os, itertools
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
RMAX=int(sys.argv[1]); AMAX=int(sys.argv[2]); KMAX=int(sys.argv[3]); LMAX=int(sys.argv[4]); OUT=sys.argv[5]
labels=[('h',k) for k in range(2,LMAX+1)]+[('s',p) for p in range(2,LMAX+1)]
cores=[c for m in range(1,KMAX+1) for c in itertools.combinations_with_replacement(labels,m)]
Dmax=KMAX*LMAX
t0=time.time(); last=t0
print(f'[{time.strftime("%H:%M:%S")}] FM3 kernel screen: {len(cores)} cores (<= {KMAX} labels <= {LMAX}), r <= {RMAX}, a <= {AMAX}',flush=True)
Fs=[]
for c in cores:
    Fs.append(core([k for t,k in c if t=='h'],[k for t,k in c if t=='s']))
neg=[]; zeros=0; evals=0
for r in range(1,RMAX+1):
    e=2*r
    for a in range(0,AMAX+1):
        N=a+e; cv=cvec(a,e); C=lambda k: cv[k] if 0<=k<=N else 0
        Wt={}
        for P in range(0,Dmax+1):
            for Q in range(0,P+1):
                if (N+P+Q)%2: continue
                j=(N+P-Q)//2; i=(N+P+Q)//2+1
                Wt[(P,Q)]=(C(i-1)+C(i+1))*C(j)-C(i)*(C(j-1)+C(j+1))
        for ci,F in enumerate(Fs):
            tot=0; hit=False
            for (p,q),v in F.items():
                w=Wt.get((p,q) if p>=q else (q,p))
                if w is not None: tot+=v*w; hit=True
            if not hit: continue
            evals+=1
            if tot<0: neg.append((cores[ci],r,a,tot//2))
            elif tot==0: zeros+=1
        if time.time()-last>60:
            last=time.time(); print(f'[{time.strftime("%H:%M:%S")}] HEARTBEAT r={r}/{RMAX} a={a}/{AMAX} evals={evals} negatives={len(neg)}',flush=True)
    print(f'[{time.strftime("%H:%M:%S")}] PROGRESS level r={r}/{RMAX} done; evals {evals}, negatives {len(neg)} {neg[:2]}, zeros {zeros} ({time.time()-t0:.0f}s)',flush=True)
json.dump(dict(evals=evals,negatives=[(list(map(list,c)),r,a,str(v)) for c,r,a,v in neg[:200]],nneg=len(neg),zeros=zeros),open(OUT,'w'))
print('DONE evals',evals,'negatives',len(neg),flush=True)
