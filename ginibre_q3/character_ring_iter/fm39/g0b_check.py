# Independent check of Theorem G0B with the definition-level evaluator: every G0 word with (a-2r+2)^2 <= 8r-9 has phi >= 0,
# and the LE counterexample a=0, e=21, x=9, C=4.
exec(open('split3.py').read())
import math
n=bad=0
for r in range(2,14):
    e=2*r-3
    for a in range(0,60):
        if (a-2*r+2)**2>8*r-9: continue
        N=a+e
        for w in range(2,N+2):
            for v in range(w,N+4):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    A,B,C=u+1,v+1,w+1; al=(N+A-B-C)//2
                    if not(0<=al and al+C<=N and al+B>=N+1): continue
                    val=phi3(r,a,u,v,w); n+=1; bad+= val<0
print('G0 words in band B checked:',n,' negative:',bad)
c=cv(0,21); C_=lambda k: c[k] if 0<=k<=21 else 0
x,Cw=9,4; S=sum(C_(k)**2-C_(k-1)*C_(k+1) for k in range(x,x+Cw+1)); T=C_(x)*C_(x+Cw)-C_(x-1)*C_(x+Cw+1)
gx=C_(x)**2+C_(x-1)**2; ge=C_(x+Cw+1)**2+C_(x+Cw)**2
print('LE counterexample: S,T =',S,T,' S^2 - |g_x|^2|g_end|^2 =',S*S-gx*ge)
