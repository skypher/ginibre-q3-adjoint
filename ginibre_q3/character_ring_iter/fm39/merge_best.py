# For each three-factor word, over the three groupings (singleton Z in {A,B,C}) take the smallest correction
#   M_Z = W(U_X U_Z, U_Y) + W(U_Y U_Z, U_X);  phi = sum_(d in CG(X,Y)) phi2(d,Z) - M_Z.  If min_Z M_Z <= 0 then (E) => phi >= 0.
exec(open('split3.py').read())
from collections import Counter
st=Counter(); bad=[]
def corr(a,e,X,Y,Z): return sum(W(a,e,d,Y) for d in cg(X,Z))+sum(W(a,e,d,X) for d in cg(Y,Z))
for r in range(2,8):
    e=2*r-3
    for a in range(0,30):
        N=a+e
        for w in range(2,12):
            for v in range(w,16):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    A,B,C=u+1,v+1,w+1
                    Ms={'C':corr(a,e,A,B,C),'B':corr(a,e,A,C,B),'A':corr(a,e,B,C,A)}
                    best=min(Ms.values()); st['words']+=1
                    if best<=0: st['some grouping <= 0']+=1; st['best='+min(Ms,key=Ms.get)]+=1
                    else:
                        st['all groupings > 0']+=1
                        if len(bad)<6: bad.append((r,a,(u,v,w),Ms,phi3(r,a,u,v,w)))
print(dict(st)); print('words with every grouping > 0:',bad)
