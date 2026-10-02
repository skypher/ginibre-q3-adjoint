# FM-STR5d main-agent check: on pair-free lists with labels <= 3 (every sign pattern, every multiplicity) up to N factors,
# some pair has nonnegative interior PREFIXES (IP) although LP (all layers >= 0) fails for every pair on some lists.
# Usage: python3 -u str5d_ip_labels3_main.py N   (N = 24: 7,640 lists, no IP failure, 32 lists without an LP pair).
# IP (interior prefix) for every pair on (+-1)^a (+-2)^b (+-3)^c families, pair-free, even minus count, up to N factors.
import itertools, sys, time
import pathlib
src=open(pathlib.Path(__file__).resolve().parent/"str13_layer_positivity_main.py").read(); src=src[:src.index("LOG=pathlib")]
ns={}; exec(src,ns); layers=ns['layers']
def prefix_ok(lay):
    run=0
    for t in sorted(lay):
        run+=lay[t]
        if run<0: return False
    return True
Nmax=int(sys.argv[1]); t0=time.time(); tot=0; ip_none=[]; lp_none=0
for s1,s2,s3 in itertools.product((1,-1),repeat=3):
    for a in range(0,Nmax+1):
        for b in range(0,Nmax+1-a):
            for c in range(0,Nmax+1-a-b):
                N=a+b+c
                if N<4: continue
                L=[s1]*a+[2*s2]*b+[3*s3]*c
                if sum(1 for z in L if z<0)%2 or sum(map(abs,L))%2: continue
                tot+=1; seen=set(); anyip=False; anylp=False
                for i in range(N):
                    for j in range(i+1,N):
                        key=tuple(sorted((L[i],L[j])))
                        if key in seen: continue
                        seen.add(key); lay=layers(L,i,j)
                        if prefix_ok(lay): anyip=True
                        if all(v>=0 for v in lay.values()): anylp=True
                if not anyip: ip_none.append(L); print("IP FAILS FOR ALL PAIRS",L,flush=True)
                if not anylp: lp_none+=1
print(time.strftime('%H:%M:%S'),"lists",tot,"no IP pair",len(ip_none),"no LP pair",lp_none,f"{time.time()-t0:.0f}s")
