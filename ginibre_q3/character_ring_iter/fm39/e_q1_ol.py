# (E) at q = 1 from Theorem OL on neighbouring rows:
#   D_j - D_(j+2) - W_(j+2,j) = D^f_(j+1) - D^f_(j+2),  f = (1-z)P = (1+z)^a (1-z)^(e+1)
#   D_j - D_(j+2) + W_(j+2,j) = D^s_(j+1) - D^s_(j+2),  s = (1+z)P = (1+z)^(a+1) (1-z)^e
exec(open('split3.py').read())
n=ok=0; olok=0
def Dr(c,k):
    N=len(c)-1; C_=lambda t: c[t] if 0<=t<=N else 0
    return C_(k)**2-C_(k-1)*C_(k+1)
for a in range(0,41):
    for e in range(0,41):
        c=list(cv(a,e)); N=a+e; f=list(cv(a,e+1)); s=list(cv(a+1,e))
        C_=lambda k: c[k] if 0<=k<=N else 0; B=lambda k: C_(k-1)+C_(k+1)
        for j in range((N+1)//2, N):
            i=j+2; W=B(i)*C_(j)-C_(i)*B(j); L=Dr(c,j)-Dr(c,i)
            m1=Dr(f,j+1)-Dr(f,j+2); p1=Dr(s,j+1)-Dr(s,j+2)
            n+=1; ok+= (L-W==m1) and (L+W==p1)
            olok+= (m1>=0 and p1>=0)
print('q=1 pairs',n,' identities hold:',ok,' both OL differences >= 0:',olok)
