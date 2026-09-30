import sys, time, os, json
sys.argv_saved=sys.argv; sys.argv=[sys.argv[0]]
exec(open('qgsi.py').read().rsplit("if __name__=='__main__':",1)[0])
def parts(n,m,l):
    if n==0: yield (); return
    if l==0: return
    for k in range(min(n,m),0,-1):
        for r in parts(n-k,k,l-1): yield (k,)+r
Nmax=int(sys.argv_saved[1]); smin=int(sys.argv_saved[2])
words=[lam for n in range(2,Nmax+1,2) for lam in parts(n,n,n) if len(lam)>=smin and len(lam)<=(10 if n<=10 else 9)]
done={}
if os.path.exists('qscan.results'):
    for line in open('qscan.results'): r=json.loads(line); done[tuple(r['kap'])]=r
threading.Thread(target=heartbeat,args=(60.0,),daemon=True).start()
bad=[r['kap'] for r in done.values() if not r['ok']]
log(f'QSCAN {len(words)} words (size<={Nmax}, >= {smin} parts); {len(done)} restored')
for j,lam in enumerate(words):
    if lam in done: continue
    t=time.time(); ok,info=runq(lam)
    rec=dict(kap=list(lam),ok=bool(ok),**{k:int(v) for k,v in info.items()},sec=round(time.time()-t,1))
    with open('qscan.results','a') as f: f.write(json.dumps(rec)+'\n'); f.flush(); os.fsync(f.fileno())
    if not ok: bad.append(list(lam))
    log(f'PROGRESS {j+1}/{len(words)} last {lam} {rec["sec"]}s; failures {len(bad)}: {bad[-5:]}')
log('DONE', 'failures', bad)
