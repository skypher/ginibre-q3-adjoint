# FM-STR14 main-agent scan: python3 -u str14_ip_labels4_main.py N0 N1 (10 22: 67,552 lists, ipfail 0, about 9 minutes).
# IP for some pair on pair-free lists with labels <= 4 (one sign per label class), N in [N0,N1]; per-N lines + 60 s heartbeat.
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
N0,N1=int(sys.argv[1]),int(sys.argv[2]); hb=time.time()
for N in range(N0,N1+1):
    tN=time.time(); tot=0; ipfail=0; lpfail=0
    for s in itertools.product((1,-1),repeat=4):
        for comp in itertools.product(range(N+1),repeat=3):
            if sum(comp)>N: continue
            m=list(comp)+[N-sum(comp)]
            L=[s[k]*(k+1) for k in range(4) for _ in range(m[k])]
            if sum(1 for z in L if z<0)%2 or sum(map(abs,L))%2: continue
            tot+=1; seen=set(); anyip=anylp=False
            for i in range(N):
                for j in range(i+1,N):
                    key=tuple(sorted((L[i],L[j])))
                    if key in seen: continue
                    seen.add(key); lay=layers(L,i,j)
                    if not anyip and prefix_ok(lay): anyip=True
                    if not anylp and all(v>=0 for v in lay.values()): anylp=True
                    if anyip and anylp: break
                if anyip and anylp: break
            if not anyip: ipfail+=1; print("IP FAILS FOR ALL PAIRS",L,flush=True)
            lpfail+= not anylp
            if time.time()-hb>60: print(time.strftime('%H:%M:%S'),f"heartbeat N={N} lists {tot} ipfail {ipfail} lpfail {lpfail}",flush=True); hb=time.time()
    print(time.strftime('%H:%M:%S'),f"N={N} lists={tot} ipfail={ipfail} lpfail={lpfail} seconds={time.time()-tN:.0f}",flush=True)
