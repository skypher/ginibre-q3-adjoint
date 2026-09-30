# W from (E): for base rows (e odd) and windows x >= (N-C)/2:
#   P_C(x) = sum_(y=x..N) [D_y - D_(y+C+1) - W_(y+C+1,y)],  each term an (E) quantity T - W with i + j >= N+1;
# and (E) for j < N/2 reduces to (E) at (N-j, i) via D_j = D_(N-j), W_(i,j) = (-1)^e W_(i,N-j).
exec(open('split3.py').read())
from collections import Counter
st=Counter()
for r in range(2,9):
    e=2*r-3
    for a in range(0,40):
        N=a+e; c=cv(a,e); C_=lambda k: c[k] if 0<=k<=N else 0
        D=lambda k: C_(k)**2-C_(k-1)*C_(k+1); Bk=lambda k: C_(k-1)+C_(k+1)
        Wij=lambda i,j: Bk(i)*C_(j)-C_(i)*Bk(j)
        P=lambda x,C: sum(D(k) for k in range(x,x+C+1))-C_(x)*C_(x+C)+C_(x-1)*C_(x+C+1)
        for C in range(1,N+1):
            for x in range(0,N-C+1):
                if 2*x<N-C: continue
                steps=[D(y)-D(y+C+1)-Wij(y+C+1,y) for y in range(x,N+1)]
                st['windows']+=1
                st['telescoping ok']+= sum(steps)==P(x,C)
                st['all steps >= 0']+= all(s>=0 for s in steps)
        for j in range(0,N+1):
            for i in range(j+1,N+2):
                if i+j<N+1 or 2*j>=N: continue
                st['pairs j<N/2']+=1
                st['mirror identity']+= (D(j)==D(N-j) and Wij(i,j)==(-1)**e*Wij(i,N-j))
print(dict(st))
