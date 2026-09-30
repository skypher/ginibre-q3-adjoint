exec(open('split3.py').read())
def Wq(c,p,q,eparity):
    # quadratic W formula for p>=q; antisymmetry/symmetry by reciprocity type for p<q
    N=len(c)-1; C_=lambda k: c[k] if 0<=k<=N else 0
    if (N+p+q)%2: return 0
    if p<q: return (-1)**eparity*Wq(c,q,p,eparity)
    j=(N+p-q)//2; i=(N+p+q)//2+1
    return (C_(i-1)+C_(i+1))*C_(j)-C_(i)*(C_(j-1)+C_(j+1))
n=0
for a in range(0,7):
    for e in range(0,8):
        c=list(cv(a,e)); N=a+e
        for p in range(0,N+3):
            for q in range(0,N+3):
                if (N+p+q)%2: continue
                assert Wq(c,p,q,e%2)==W(a,e,p,q),(a,e,p,q,Wq(c,p,q,e%2),W(a,e,p,q)); n+=1
print('quadratic W formula matches direct W:',n,'cases, both parities of e')
