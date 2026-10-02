import re, gzip, sys, time
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent)); from str7e_budget_runs_main import *
LOG=str(pathlib.Path(__file__).resolve().parent/"sec166_census_w40_noflip.log.gz")
rows=[]
for line in gzip.decompress(open(LOG,'rb').read()).decode().splitlines():
    if not line.startswith("NOFLIP"): continue
    p=int(re.search(r"p=(-?\d+)",line).group(1)); B=[int(x) for x in re.search(r"B=([-\d ]+?)\s+phi",line).group(1).split()]
    rows.append((B,p))
stat=dict(plain_bad=0,cum_bad=0,height_bad=0,Pneg=0,n=0); worst=(9,None); t0=time.time()
for B,p in rows:
    L=B+[p]; i,j=toppair(B)
    zu,zv=L[i],L[j]; a,b=abs(zu),abs(zv); eu=1 if zu>0 else -1; ev=1 if zv>0 else -1
    C=[z for q,z in enumerate(L) if q not in (i,j)]; A,Bb=interior_cut(C)
    rA=table(A); FB=table(Bb); Gp=mul_xx(FB,a,b,eu*ev); Gm=mul_xy(FB,a,b,eu,ev)
    P=defaultdict(int); M=defaultdict(int)
    for (r,s),v in rA.items(): P[r+s]+=v*Gp.get((r,s),0); M[r+s]+=v*Gm.get((r,s),0)
    run=0; plain=0; hb=False; pn=False; cb=False; pb=False
    for t in sorted(set(P)|set(M)):
        run+=P[t]+min(M[t],0)
        if run<0: cb=True
        plain+=P[t]+M[t]
        if plain<0: pb=True
        if P[t]+min(M[t],0)<0: hb=True
        if P[t]<0: pn=True
    phi=sum(P.values())+sum(M.values())
    stat['n']+=1; stat['plain_bad']+=pb; stat['cum_bad']+=cb; stat['height_bad']+=hb; stat['Pneg']+=pn
    negmix=-sum(min(m,0) for m in M.values()); pos=sum(max(x,0) for x in P.values())
    r=negmix/pos if pos else 9
    if r>worst[0] or worst[1] is None: pass
    if worst[1] is None or r>worst[0]: worst=(r,(B,p,phi))
print(time.strftime('%H:%M:%S'),stat,"max (neg mixed)/(pos pure)",worst[0],worst[1],f"{time.time()-t0:.0f}s")
