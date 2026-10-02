# FM-STR5d main-agent range scan: python3 -u str5d_ip_labels3_range_main.py N0 N1 (25 32: 9,304 lists, ipfail 0).
# IP / LP for every pair on pair-free labels<=3 lists with N in [N0, N1] factors; per-N lines, 60 s heartbeat, results appended per N.
import itertools, sys, time, os
import pathlib
src=open(pathlib.Path(__file__).resolve().parent/"str13_layer_positivity_main.py").read(); src=src[:src.index("LOG=pathlib")]
ns={}; exec(src,ns); layers=ns['layers']
def prefix_ok(lay):
    run=0
    for t in sorted(lay):
        run+=lay[t]
        if run<0: return False
    return True
N0,N1=int(sys.argv[1]),int(sys.argv[2]); out=open(os.devnull,"w")
hb=time.time(); start=time.time()
for N in range(N0,N1+1):
    tN=time.time(); tot=0; ipfail=[]; lpfail=0
    cases=[(s,a,b,N-a-b) for s in itertools.product((1,-1),repeat=3) for a in range(N+1) for b in range(N+1-a)]
    for k,(s,a,b,c) in enumerate(cases):
        L=[s[0]]*a+[2*s[1]]*b+[3*s[2]]*c
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
        if not anyip: ipfail.append(L); print("IP FAILS FOR ALL PAIRS",L,flush=True)
        lpfail+= not anylp
        if time.time()-hb>60:
            print(time.strftime('%H:%M:%S'),f"heartbeat N={N} case {k+1}/{len(cases)} lists so far {tot} ipfail {len(ipfail)} lpfail {lpfail}",flush=True); hb=time.time()
    line=f"N={N} lists={tot} ipfail={len(ipfail)} lpfail={lpfail} seconds={time.time()-tN:.0f}"
    print(time.strftime('%H:%M:%S'),line,flush=True); out.write(line+"\n"); out.flush(); os.fsync(out.fileno())
print(time.strftime('%H:%M:%S'),"DONE",f"{time.time()-start:.0f}s",flush=True)
