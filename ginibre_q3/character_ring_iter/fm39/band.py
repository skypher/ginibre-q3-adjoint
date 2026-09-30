# Coverage map + margins for three-factor H-only words phi_r(h_u h_v h_w h_1^a), by the proved criteria.
import sys, itertools, math, json, time
from fractions import Fraction as Fr
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
RMAX=int(sys.argv[1]); AMAX=int(sys.argv[2]); KMAX=int(sys.argv[3])
def lH(k): return Fr(2*k*k*(k+1)*(k+3),3) if k%2==0 else Fr((k-1)**2*(k+2)*(k+3),6)
def covered(kap,r,a):
    t=sum(k%2 for k in kap); m=len(kap)
    if (a+t)%2: return 'zero'
    if a==0: return 'THREE'
    Lam=Fr(sum(k*(k+4) for k in kap),5)
    if a+2*r+4 >= (2*r+3)*Lam: return 'LS'
    R=(a+t)//2; e=2*r-m
    if R>=1:
        Lp=Fr(sum(k*(k+4) for k in kap if k%2),5)+Fr(sum((k+1)*(k+3) for k in kap if k%2==0),3)
        if e+2*R+4 >= (2*R+3)*Lp: return 'LSdual'
    if 2*r+a+t+6 >= 7*sum(lH(k) for k in kap): return 'LR4'
    if a==2*r and set(kap)<= {kap[0]-kap[0]%2, kap[0]-kap[0]%2+1} and min(kap)>=2*r and sorted(k%2 for k in kap) in ([0,0,0],[0,1,1]): return 'BW'
    return None
stats={}; unc=[]; t0=time.time()
Fcache={}
for r in range(2,RMAX+1):
    for a in range(0,AMAX+1):
        for kap in itertools.combinations_with_replacement(range(1,KMAX+1),3):
            c=covered(kap,r,a); stats[c]=stats.get(c,0)+1
            if c is None:
                if kap not in Fcache: Fcache[kap]=core(kap,())
                val=phi_kernel(Fcache[kap],r,a)
                base=phi_kernel({(0,0):1},r,a+sum(k%2 for k in kap))   # phi_r(h_1^(a+t))
                dimw=math.prod(math.comb(k+3,3) for k in kap)
                unc.append((kap,r,a,val, float(Fr(val,base)/dimw) if base else None))
    print(f'[{time.strftime("%H:%M:%S")}] r={r} done; stats {stats}; uncovered so far {len(unc)}; negatives {sum(1 for u in unc if u[3]<0)} ({time.time()-t0:.0f}s)',flush=True)
neg=[u for u in unc if u[3]<0]
print('NEGATIVES:',neg[:10])
rel=sorted(unc,key=lambda u:u[4])[:12]
print('smallest relative margins phi/(dim(w) phi_r(h_1^(a+t))) among uncovered:')
for u in rel: print('  ',u[0],'r=',u[1],'a=',u[2],'phi=',u[3],'rel=%.3e'%u[4])
json.dump([[list(u[0]),u[1],u[2],str(u[3]),u[4]] for u in unc],open(f'band_{RMAX}_{AMAX}_{KMAX}.json','w'))
