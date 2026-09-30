# Merge move: phi(h_u h_v * F) = sum_(d in CG(A,B)) phi(hat S_d * F) - E[K (U_A(x)U_B(y) + U_B(x)U_A(y)) F]/2 at the same kernel.
# Check the size/sign of the mixed correction over three-factor words (F = h_w), and how often it vanishes.
exec(open('split3.py').read())
from collections import Counter
st=Counter()
for r in range(2,7):
    e=2*r-3
    for a in range(0,25):
        N=a+e
        for w in range(2,10):
            for v in range(w,14):
                for u in range(v,v+w+2):
                    if (a+u+v+w)%2: continue
                    A,B,C=u+1,v+1,w+1
                    # mixed correction M = W(U_A U_C, U_B) + W(U_B U_C, U_A) (both cross terms beyond the grouped one)
                    Mx=sum(W(a,e,d,B) for d in cg(A,C))+sum(W(a,e,d,A) for d in cg(B,C))
                    st['words']+=1; st['correction 0']+= Mx==0; st['correction >0 (hurts)']+= Mx>0; st['correction <0 (helps)']+= Mx<0
print(dict(st))
