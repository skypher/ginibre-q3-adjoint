import sys, time, itertools
exec(open('n0check.py').read().split("if __name__=='__main__':")[0])
R=int(sys.argv[1]); DEG=int(sys.argv[2]); t0=time.time()
def msets(vals,maxn,maxdeg):
    out=[()]
    for n in range(1,maxn+1):
        for c in itertools.combinations_with_replacement(vals,n):
            if sum(c)<=maxdeg: out.append(c)
    return out
cores=[]
for H in msets(range(2,8),2*R,DEG):
    for S in msets(range(2,6),3,DEG-sum(H)):
        if sum(H)+sum(S)>=1: cores.append((H,S))
print(f'[{time.strftime("%H:%M:%S")}] r={R}: {len(cores)} cores, degree<={DEG}',flush=True)
fails=[]; nonpoly=[]; neg0=[]; done=0
for H,S in cores:
    w=ONE
    for k in H: w=mul(w,h(k))
    for p in S: w=mul(w,shat(p))
    ok,mn,co,v0=test(R,w,sum(H)+sum(S))
    if not ok: nonpoly.append((H,S))
    if mn<0: fails.append((H,S,[float(c) for c in co]))
    if v0<0: neg0.append((H,S))
    done+=1
    if done%200==0: print(f'[{time.strftime("%H:%M:%S")}] {done}/{len(cores)} newton-fails {len(fails)} nonpoly {len(nonpoly)} P(0)<0 {len(neg0)}',flush=True)
print(f'[{time.strftime("%H:%M:%S")}] r={R} DONE {done} cores; Newton failures {len(fails)}; nonpolynomial {len(nonpoly)}; negative P(0) {len(neg0)} ({time.time()-t0:.0f}s)')
for f in fails[:8]: print('  FAIL H=',f[0],'S=',f[1],'newton=',[round(c,1) for c in f[2]])
