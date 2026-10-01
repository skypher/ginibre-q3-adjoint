"""FM-MECH56: exact directional gluing and remaining-set verifier."""
import argparse
from fractions import Fraction as F
from math import comb

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--max-n", type=int, default=60)
args = ap.parse_args()

def data(a, e):
    N = a + e
    de = a - e
    sig = N + 2
    c = [1]
    old = 0
    for k in range(N):
        z, rem = divmod(de*c[-1] - (N-k+1)*old, k+1)
        assert rem == 0
        old = c[-1]
        c.append(z)
    v = lambda k: c[k] if 0 <= k <= N else 0
    B = [v(k-1)+v(k+1) for k in range(N+2)]
    D = [v(k)**2-v(k-1)*v(k+1) for k in range(N+2)]
    H = [sig*(v(k)**2+v(k-1)**2)-2*de*v(k)*v(k-1)
         for k in range(N+2)]
    for k in range(N+2):
        assert (k+1)*v(k+1) == de*v(k)-(N-k+1)*v(k-1)
    if N <= 14:
        for k in range(N+1):
            assert v(k) == sum(
                (-1)**h*comb(e,h)*comb(a,k-h)
                for h in range(max(0,k-a), min(e,k)+1))
    return v, B, D, H

def chord(v, B, j, i):
    return v(j)*B[i]-B[j]*v(i)

def signed(v, B, D, j, i, sg):
    return D[j]-D[i]+sg*chord(v,B,j,i)

def correction(v, B, D, j, k, i):
    A = chord(v,B,j,k)
    Z = -chord(v,B,j,k+1)
    delta = D[k]-D[k+1]
    G = (delta*(D[j]-D[i])
         - Z*(D[k]-D[i])-A*(D[k+1]-D[i]))
    return A, Z, delta, G

def h_test(a, e, j, i, D, H):
    sig = a+e+2
    C = sig*sig-(a-e)**2
    Y = 2*i-a-e
    return C*(sig+Y)**2*(D[j]-D[i])**2-4*Y*Y*H[j]*H[i]

rows = longs = glued = vertices = normpass = identities = remaining = 0
for N in range(8, args.max_n+1):
    for e in range(3, (N-2)//2+1):
        a = N-e
        v, B, D, H = data(a,e)
        rows += 1
        for k in range((N+1)//2, N+1):
            assert D[k] >= D[k+1] >= 0
            assert chord(v,B,k,k+1) == D[k]-D[k+1]
        for j in range(N//2+1, N-3):
            assert v(j) != 0 or B[j] != 0
            crossings = []
            opposites = []
            for i in range(j+1, N+1):
                wi = chord(v,B,j,i)
                dot = v(j)*v(i)+B[j]*B[i]
                if wi == 0 and dot < 0:
                    opposites.append(i)
                    lam = (F(-v(i),v(j)) if v(j)
                           else F(-B[i],B[j]))
                    assert lam > 0
                    assert D[i] <= lam*lam*D[j]
                k = i-1
                if k > j and chord(v,B,j,k) > 0 and wi < 0:
                    crossings.append(k)
                if i-j < 4 or not (crossings or opposites):
                    continue
                longs += 1
                good = False
                for k in crossings:
                    A, Z, delta, G = correction(v,B,D,j,k,i)
                    assert A > 0 and Z > 0 and delta > 0
                    assert delta*v(j) == -Z*v(k)-A*v(k+1)
                    assert delta*B[j] == -Z*B[k]-A*B[k+1]
                    for sg in (-1,1):
                        assert (
                            delta*signed(v,B,D,j,i,sg)
                            == Z*signed(v,B,D,k,i,-sg)
                             + A*signed(v,B,D,k+1,i,-sg)+G)
                        identities += 1
                    good |= G >= 0
                nh = h_test(a,e,j,i,D,H) >= 0
                nv = bool(opposites)
                glued += good
                vertices += nv
                normpass += nh
                remaining += not (good or nv or nh)
                assert D[j]-D[i] >= abs(wi)

print("ROWS", rows, "LONG PAIRS", longs)
print("GLUING IDENTITIES", identities)
print("COVERED BY: gluing", glued,
      "opposite vertex", vertices, "joint H", normpass)
print("REMAINING NECESSARY SET", remaining)

# Actual-row failure of propagation after discarding child margins.
a, e, j, k, i = 31, 19, 26, 28, 39
v, B, D, H = data(a,e)
A, Z, delta, G = correction(v,B,D,j,k,i)
assert A > 0 and Z > 0 and delta > 0 and G < 0
assert all(chord(v,B,j,l) > 0 for l in range(j+1,k+1))
assert F(delta,A+Z) == F(132,203)
assert F(G,delta*D[j]) == -F(81682236337,200341473858504)
assert h_test(a,e,j,i,D,H) > 0
for sg in (-1,1):
    child = (Z*signed(v,B,D,k,i,-sg)
             + A*signed(v,B,D,k+1,i,-sg))
    assert child+G == delta*signed(v,B,D,j,i,sg)
    assert child+G > 0
print("GLUING OBSTRUCTION", (a,e,j,k,i))
print("lambda", F(delta,A+Z),
      "normalized correction", F(G,delta*D[j]))
print("ACTUAL E SLACK", D[j]-D[i]-abs(chord(v,B,j,i)))

# Four-step failure of the stronger endpoint-independent drop.
a, e, j, k = 172, 45, 110, 113
v, B, D, H = data(a,e)
A, Z, delta, G0 = correction(v,B,D,j,k,a+e+1)
assert D[a+e+1] == 0
assert all(chord(v,B,j,l) > 0 for l in range(j+1,k+1))
assert A > 0 and Z > 0 and delta > 0 and G0 < 0
assert F(delta,A+Z) == F(605877353895,739919914372)
assert F(G0,delta*D[j]) == -F(
    22605171807602265,38876440477482825932)
assert D[j]-D[k+1] >= abs(chord(v,B,j,k+1))
print("FOUR-STEP DROP OBSTRUCTION", (a,e,j,k+1))
print("lambda", F(delta,A+Z),
      "normalized correction", F(G0,delta*D[j]))
print("PASS")
