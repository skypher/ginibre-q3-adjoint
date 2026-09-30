# Independent randomized exact screen of FM3 across levels: words with m <= 2r general h (labels 2..LH), 0..KS hat S
# (labels 2..LS), suffix a <= AMAX, levels r = 1..RMAX.  Exact fusion-kernel evaluator (fm3kern.py).  Heartbeat 30 s.
import sys, time, random
exec(open('fm3kern.py').read().split("if __name__")[0]) if "if __name__" in open('fm3kern.py').read() else None
from fm3kern import core, phi_kernel
seed,NW,RMAX,LH,LS,KS,AMAX=[int(x) for x in sys.argv[1:8]]
rng=random.Random(seed); t0=time.time(); last=t0; n=neg=zero=0; ex=[]; minpos=None
for w in range(NW):
    r=rng.randint(1,RMAX); m=rng.randint(0,2*r); k=rng.randint(0,KS)
    hs=tuple(sorted(rng.randint(2,LH) for _ in range(m))); ss=tuple(sorted(rng.randint(2,LS) for _ in range(k)))
    F=core(hs,ss)
    for a in range(0,AMAX+1):
        v=phi_kernel(F,r,a); n+=1
        if v<0:
            neg+=1
            if len(ex)<10: ex.append((r,hs,ss,a,v))
        elif v==0: zero+=1
        elif minpos is None or v<minpos[0]: minpos=(v,r,hs,ss,a)
    if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat words %d/%d (%.3fs/word, current r=%d m=%d k=%d) evaluations %d negative %d'%(w+1,NW,(time.time()-t0)/(w+1),r,m,k,n,neg),flush=True)
print(time.strftime('%H:%M:%S'),'done: evaluations',n,'negative',neg,'zero',zero,'min positive',minpos,'elapsed %.0fs'%(time.time()-t0)); print('negatives:',ex)
