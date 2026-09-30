# FM3 at large level r for fixed small words: near the antipodal corners x = -y = +-2 every odd-label factor vanishes,
# so the sign at large r is decided by local quadratic (or higher) forms.  Exact fusion-kernel evaluation.
# Words: up to 4 factors from h_2..h_7 and hat S_2..hat S_7 (any mix, m <= 2r automatically for r >= 2), a = 0..3.
import sys, time, itertools
from fm3kern import core, phi_kernel
t0=time.time(); last=t0; n=neg=0; ex=[]; minrel=None
labs=range(2,8)
words=[]
for nh in range(0,5):
    for ns in range(0,5-nh):
        if nh+ns==0: continue
        for hs in itertools.combinations_with_replacement(labs,nh):
            for ss in itertools.combinations_with_replacement(labs,ns):
                words.append((hs,ss))
R=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else [10,20,40]
for wi,(hs,ss) in enumerate(words):
    F=core(hs,ss)
    for r in R:
        for a in range(0,4):
            v=phi_kernel(F,r,a); n+=1
            if v<0:
                neg+=1
                if len(ex)<10: ex.append((r,hs,ss,a,v))
    if time.time()-last>30: last=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat words %d/%d evaluations %d negative %d'%(wi+1,len(words),n,neg),flush=True)
print(time.strftime('%H:%M:%S'),'done: words',len(words),'evaluations',n,'negative',neg,'elapsed %.0fs'%(time.time()-t0)); print('negatives:',ex)
