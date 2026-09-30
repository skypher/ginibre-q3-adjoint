exec(open('split3.py').read())
import time
print('Astra example (150,4;153,152,3):', phi3(150,4,153,152,3)>0)
t0=time.time(); tot=neg=zero=0; minrat=None
for r in range(2,9):
    e=2*r-3
    for a in range(0,40):
        N=a+e
        for w in range(1,N+2):
            for v in range(w,N+4):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    A,B,C=u+1,v+1,w+1; al=(N+A-B-C)//2; be=al+C; ga=al+B
                    if not(0<=al<be<=N and ga>=N+1): continue
                    val=phi3(r,a,u,v,w); tot+=1
                    if val<0: neg+=1; print('NEG',r,a,(u,v,w),val)
                    if val==0: zero+=1
print('G0 support branch (all C parities), r<=8, a<40:',tot,'words; negative',neg,'zero',zero,round(time.time()-t0,1),'s')
