# Exhaustive exact falsification screen of (E): D_j - D_i >= |B_i c_j - c_i B_j|, N/2 <= j < i <= N+2,
# c_k = [z^k](1+z)^a(1-z)^e.  Folded: a >= e (the inequality is fold-invariant).  Resumable per (e) row.
import sys, time, json, os
from math import comb
EMAX=int(sys.argv[1]); AMAX=int(sys.argv[2]); OUT='escreen.results'
done=set()
if os.path.exists(OUT):
    for l in open(OUT): done.add(json.loads(l)['e'])
t0=time.time(); last=t0
def cvec(a,e):
    N=a+e; c=[0]*(N+1)
    for i in range(e+1):
        s=(-1)**i*comb(e,i)
        for j in range(a+1): c[i+j]+=s*comb(a,j)
    return c
print(f'[{time.strftime("%H:%M:%S")}] screen e<= {EMAX}, a in [e, {AMAX}]; restored rows {sorted(done)}',flush=True)
for e in range(0,EMAX+1):
    if e in done: continue
    te=time.time(); pairs=0; fails=[]; minslack=None
    for a in range(e,AMAX+1):
        N=a+e; c=cvec(a,e)
        C=lambda k: c[k] if 0<=k<=N else 0
        Dk=[C(k)**2-C(k-1)*C(k+1) for k in range(0,N+3)]
        Bk=[C(k-1)+C(k+1) for k in range(0,N+3)]
        for j in range((N+1)//2,N+2):
            for i in range(j+1,N+3):
                T=Dk[j]-Dk[i]; W=Bk[i]*C(j)-C(i)*Bk[j]; pairs+=1
                sl=T-abs(W)
                if sl<0: fails.append((a,e,j,i))
        if time.time()-last>60:
            last=time.time(); print(f'[{time.strftime("%H:%M:%S")}] HEARTBEAT e={e} a={a}/{AMAX} pairs(row)={pairs} fails(row)={len(fails)}',flush=True)
    rec=dict(e=e,pairs=pairs,fails=fails[:20],nfails=len(fails),sec=round(time.time()-te,1))
    with open(OUT,'a') as f: f.write(json.dumps(rec)+'\n'); f.flush(); os.fsync(f.fileno())
    print(f'[{time.strftime("%H:%M:%S")}] PROGRESS e={e}/{EMAX} done: {pairs} pairs, failures {len(fails)} {fails[:3]} ({rec["sec"]}s)',flush=True)
print('DONE',flush=True)
