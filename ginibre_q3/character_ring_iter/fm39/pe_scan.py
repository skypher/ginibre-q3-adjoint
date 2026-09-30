# Broad exhaustive test of the stronger inequality (PE): T_ij >= unsigned chord, all pairs, a>=e+2 (folded), e in 3..10, N<=40.
exec(open('pe_repro.py').read().split("links=0; failures=[]")[0])
import time; t0=time.time()
fails=[]; cnt=0
for e in range(3,11):
    for a in range(e+2,41-e):
        N=a+e; z,mu,g,norm,scale=ensemble(a,e); q=e//2
        for j in range(N//2+1,N+1):
            for i in range(j+1,N+2):
                x=(2*j-N)**2; y=(2*i-N)**2
                A=Z(z,[w*abs((s-x)*(s-y)) for s,w in zip(z,mu)],q)
                T,W=TW(a,e,j,i)
                chord=scale*(y-x)*g[j]*g[i]*A/norm
                cnt+=1
                if T<chord: fails.append((a,e,j,i,float(T-chord)))
    print(f'[{time.strftime("%H:%M:%S")}] e={e}: pairs so far {cnt}, PE failures {len(fails)} {fails[:3]}',flush=True)
print('DONE pairs',cnt,'PE failures',len(fails),fails[:10])
