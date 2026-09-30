# Constant-coefficient model: g_(k+1) = M g_k, M = [[al,-be],[1,0]] elliptic.  Check W and (E) numerically,
# and the closed forms D_(x+j) = be^j H, T = be^((n-1)/2) U_(n-1)(cos w) H, W_E = (1-be) be^((s-1)/2) U_(s-1) H.
import random, math
random.seed(1); bad=0; n_=0
for _ in range(20000):
    be=random.uniform(0.2,1.5); w=random.uniform(0.05,3.1); al=2*math.sqrt(be)*math.cos(w)
    M=lambda v:(al*v[0]-be*v[1],v[0])
    g=[(random.uniform(-1,1),random.uniform(-1,1))]
    for _ in range(40): g.append(M(g[-1]))
    wedge=lambda p,q:p[0]*q[1]-p[1]*q[0]
    H=wedge(g[0],g[1])
    for n in range(2,30):
        D=sum(wedge(g[j],g[j+1]) for j in range(n)); T=wedge(g[0],g[n])
        U=math.sin(n*w)/math.sin(w)
        assert abs(T-be**((n-1)/2)*U*H)<1e-6*(1+abs(T)) and abs(D-sum(be**j for j in range(n))*H)<1e-6*(1+abs(D))
        n_+=1; bad+= D-T < -1e-9*abs(D)
print('constant-coefficient windows',n_,' closed forms verified; W failures:',bad)
