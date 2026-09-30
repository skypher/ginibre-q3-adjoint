# gamma <= N branch of three-factor words (all three corrections present): size, tightest cases, outside the DS/AS strips and a = 0.
exec(open('split3.py').read())
import time
t0=time.time(); tot=0; zero=0; neg=0; tight=[]
for r in range(2,7):
    e=2*r-3
    for a in range(1,30):
        if abs(a-e)<=1: continue
        N=a+e; c=cv(a,e)
        for w in range(2,N+2):
            for v in range(w,N+4):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    A,B,C=u+1,v+1,w+1; al=(N+A-B-C)//2; ga=al+B
                    if not(0<=al and ga<=N and al+C<=N): continue
                    val=phi3(r,a,u,v,w); tot+=1
                    zero+= val==0; neg+= val<0
                    # scale: main term tau(U_A U_B U_C)
                    main=sum(W(a,e,l,0) for l1 in cg(A,B) for l in cg(l1,C))
                    if main>0: tight.append((val/main,r,a,(u,v,w),val,main))
tight.sort()
print('gamma<=N words outside strips, r<=6, a<30:',tot,'zero',zero,'negative',neg,'time',round(time.time()-t0,1))
print('tightest (phi/main term):')
for t in tight[:8]: print('  ratio %.3e  r=%d a=%d word=%s phi=%d main=%d'%t)
