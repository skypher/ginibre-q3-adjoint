from math import comb
from collections import defaultdict

def U(n):
    return {n-2*j: (-1)**j * comb(n-j, j) for j in range(n//2+1)}

def cat(k):
    return comb(2*k, k) // (k+1)

def moment(n):
    return 0 if n % 2 else cat(n//2)

def A(C, a, b):
    P = {(0,0): 1}
    for z in C:
        n = abs(z)
        eps = 1 if z > 0 else -1
        q = defaultdict(int)
        for (i,j),v in P.items():
            for d,c in U(n).items():
                q[i+d,j] += v*c
                q[i,j+d] += eps*v*c
        P = {k:v for k,v in q.items() if v}
    return sum(v*ca*cb*moment(i+da)*moment(j+db)
               for (i,j),v in P.items()
               for da,ca in U(a).items()
               for db,cb in U(b).items())

B = [1,1,-2] + [3]*4 + [-4]*3 + [5]*4
p = 8
Lam = B + [p]
values = list(dict.fromkeys(Lam))
rows = []
for i,u in enumerate(values):
    for v in values[i:]:
        if u == v and Lam.count(u) < 2:
            continue
        C = Lam.copy()
        C.remove(u)
        C.remove(v)
        a = A(C, abs(u), abs(v))
        rows.append((u,v,a, a if v > 0 else -a))

print("g_8(B) =", A(B,8,0))
print("all flip quantities =", rows)
assert A(B,8,0) == 4075371
assert len(rows) == 19 and all(x[3] < 0 for x in rows)
C20 = list(range(4,21)) + [169]
C80 = [2]*28 + [3]*30 + [4]*20 + [214]
print("corrected former zero A_23 =", A(C20,2,3))
print("corrected former zero A_22 =", A(C80,2,2))
assert A(C20,2,3) == 243634987
assert A(C80,2,2) == 118588718722
