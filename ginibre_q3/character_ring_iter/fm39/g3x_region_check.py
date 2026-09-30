# Targeted independent check of MECH22: enumerate words in (GX), (GOdd) and the gamma=N extensions,
# evaluate phi directly from the split definition (split3.py), require phi>0.
exec(open('split3.py').read())
import time, sys
t0=time.time(); last=[t0]
def hb(tag,done):
    if time.time()-last[0]>30:
        last[0]=time.time(); print(time.strftime('%H:%M:%S'),'heartbeat',tag,'cases so far',done,flush=True)
def params(r,a,u,v,w):
    e=2*r-3; N=a+e; A,B,C=u+1,v+1,w+1
    al=(N+A-B-C)//2; return e,N,A,B,C,al,al+C,al+B,A-B+C
cnt={'GX':0,'GOdd':0,'EXT':0}; bad=0; rootcross=0
def test(tag,r,a,u,v,w):
    global bad
    val=phi3(r,a,u,v,w); cnt[tag]+=1
    if val<=0: bad+=1; print('FAIL',tag,r,a,(u,v,w),val,flush=True)
RS=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else list(range(12,40))
for r in RS:                      # GX: needs e >= 3mY^2-4 >= 20
    e=2*r-3
    for a in (2,4,6):
        m=min(a,e); N=a+e
        for w in range(1,12,2):
            for v in range(w,N+8):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    e_,N_,A,B,C,al,be,ga,Y=params(r,a,u,v,w)
                    if C%2==0 and m%2==0 and m>=2 and 0<=al<be<=N and ga>=N+1 and 3*m*Y*Y<=N-m+4:
                        test('GX',r,a,u,v,w); rootcross+= (A-B-C<0)
                    hb('GX r=%d a=%d'%(r,a),cnt)
print(time.strftime('%H:%M:%S'),'GX done',cnt,'with d<0 (window crosses X=0):',rootcross,'bad',bad,flush=True)
for r in range(2,6):              # GOdd
    e=2*r-3
    for a in range(0,90):
        N=a+e
        for w in range(0,12,2):
            for v in range(w,N+4):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    e_,N_,A,B,C,al,be,ga,Y=params(r,a,u,v,w)
                    if C%2==1 and 0<=al<be<=N and ga>=N+1 and (e+1)*(Y+1)**2<=2*(a+2): test('GOdd',r,a,u,v,w)
                    hb('GOdd r=%d a=%d'%(r,a),cnt)
print(time.strftime('%H:%M:%S'),'GOdd done',cnt,'bad',bad,flush=True)
for r in list(range(2,7))+RS:     # gamma=N extensions
    e=2*r-3
    for a in range(0,60 if r<7 else 7):
        N=a+e; m=min(a,e)
        for w in range(0,12):
            for v in range(w,N+4):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    e_,N_,A,B,C,al,be,ga,Y=params(r,a,u,v,w)
                    if not(0<=al<be<=N and ga==N): continue
                    if (C%2==1 and r%2==0 and (e+1)*(Y+1)**2<=2*(a+2)) or \
                       (C%2==0 and m%2==0 and m>=2 and (al+m//2)%2==0 and 3*m*Y*Y<=N-m+4): test('EXT',r,a,u,v,w)
                    hb('EXT r=%d a=%d'%(r,a),cnt)
print(time.strftime('%H:%M:%S'),'EXT done',cnt,'bad',bad,'elapsed',round(time.time()-t0,1),flush=True)
