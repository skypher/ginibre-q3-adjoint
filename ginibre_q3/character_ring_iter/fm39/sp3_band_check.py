# Independent check of FM-SEC23's band: three-label words with hat S on the support branch gamma >= N+1, 0 <= alpha, alpha+Z <= N,
# inside the band ((a-e)^2 <= 2(N+1) for Z >= 3; (a-e)^2 <= N+1 for Z = 2), are >= 0.  Uses general_row_words.phi_row.
exec(open('general_row_words.py').read().split("random.seed(")[0])
from collections import Counter
st=Counter()
for r in range(1,7):
    for pat,(nh,ns) in {'A':(2,1),'B':(1,2),'C':(0,3)}.items():
        e=2*r-nh
        if e<0: continue
        for a in range(0,16):
            c=mkrow(a,e,[]); N=a+e
            if not((a-e)**2<=2*(N+1)): continue
            for labs in __import__('itertools').product(range(2,9),repeat=3):
                Ls=[l+1 for l in labs[:nh]]; Ps=list(labs[nh:])   # h_u -> U_(u+1), u>=2 ; hatS_p -> U_p, p>=2
                if any(l<3 for l in Ls): continue
                U=sorted(Ls+Ps,reverse=True); X,Y,Z=U
                if (N+X+Y+Z)%2: continue
                al=(N+X-Y-Z)//2; ga=(N+X+Y-Z)//2
                if not(ga>=N+1 and al>=0 and al+Z<=N): continue
                if Z==2 and not((a-e)**2<=N+1): continue
                v=phi_row(c,Ls,Ps,e%2); st[pat]+=1; st[pat+' negative']+= v<0
print(dict(st))
