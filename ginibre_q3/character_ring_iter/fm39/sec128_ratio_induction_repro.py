from fractions import Fraction as Q
from math import comb, factorial
from sympy import symbols, simplify

D = 30  # m + k + alpha + gamma <= 30 (total degree <= 60)

def cat(n):
    return 0 if n < 0 else comb(2*n, n) // (n + 1)

def mu(m, k):
    return Q(2*factorial(2*m)*factorial(2*m+1)*factorial(2*k)*factorial(2*k+1),
             factorial(m)**2*factorial(k)**2*factorial(m+k+1)*factorial(m+k+2))

def mu_direct(m, k):
    out = Q(0)
    for i in range(2*m + 1):
        for j in range(2*k + 1):
            a, b = i+j, 2*m+2*k-i-j
            if a % 2 == b % 2 == 0:
                out += Q(comb(2*m,i)*comb(2*k,j)*(-1)**(2*k-j)) * cat(a//2) * cat(b//2)
    return out

# Compare the joint moment formula with direct Catalan expansion.
for m0 in range(7):
    for k0 in range(7):
        assert mu(m0,k0) == mu_direct(m0,k0)

# Symbolic base and distance-two formulas.
m, k, p = symbols('m k p', integer=True, nonnegative=True)
q = symbols('q', integer=True)
den = (m+k+2)*(m+k+3)
rho_s = 4*(2*m+1)*(2*m+3)/den
rho_d = 4*(2*k+1)*(2*k+3)/den
R10 = lambda a,b: (2*(a-b)**2+8*a**2-2*b+14*a)/((a+b+2)*(a+b+3))
R01 = lambda a,b: (2*(a-b)**2+8*b**2-2*a+14*b)/((a+b+2)*(a+b+3))
assert simplify(3*rho_s/4+rho_d/4-2-R10(m,k)) == 0
assert simplify(rho_s/4+3*rho_d/4-2-R01(m,k)) == 0

D2 = (m+k+2)*(m+k+3)**2*(m+k+4)
P20 = ((5*m**2-2*m*k+k**2)**2
       + 8*m*(15*m**2-7*m*k+2*k**2)
       + (197*m**2-46*m*k+5*k**2) + 144*m+54)
P11 = (5*(m**4+k**4)-12*m*k*(m**2+k**2)+30*m**2*k**2
       +12*(m**3+k**3)+12*m*k*(m+k)+13*(m**2+k**2)
       +18*m*k+12*(m+k)+18)
R20 = simplify(3*rho_s*R10(m+1,k)/4 + rho_d*R10(m,k+1)/4 - 2*R10(m,k))
R11 = simplify(rho_s*R10(m+1,k)/4 + 3*rho_d*R10(m,k+1)/4 - 2*R10(m,k))
R02 = simplify(rho_s*R01(m+1,k)/4 + 3*rho_d*R01(m,k+1)/4 - 2*R01(m,k))
assert simplify(R20-4*P20/D2) == 0
assert simplify(R11-4*P11/D2) == 0
assert simplify(R02-4*P20.subs({m:k,k:m}, simultaneous=True)/D2) == 0
assert simplify(P11.subs({m:(p+q)/2,k:(p-q)/2}, simultaneous=True)
                -(p**4+6*p**3+11*p**2+6*p*q**2+12*p+4*q**4+2*q**2+18)) == 0
assert 4*15*2-7**2 == 71 and 4*197*5-46**2 == 1824

# Build the exact table by the two multiplier recurrences.
F = {(i,j,0,0): mu(i,j) for i in range(D+1) for j in range(D+1-i)}
for n in range(D):
    for a in range(n+1):
        g = n-a
        for i in range(D-n+1):
            for j in range(D-n-i+1):
                if i+j+n+1 <= D:
                    old = F[(i,j,a,g)]
                    F[(i,j,a+1,g)] = (3*F[(i+1,j,a,g)] + F[(i,j+1,a,g)])/4 - 2*old
                    F[(i,j,a,g+1)] = (F[(i+1,j,a,g)] + 3*F[(i,j+1,a,g)])/4 - 2*old

# Verify path independence where alpha and gamma are both positive.
for (i,j,a,g), value in F.items():
    if a and g:
        via_a = (3*F[(i+1,j,a-1,g)] + F[(i,j+1,a-1,g)])/4 - 2*F[(i,j,a-1,g)]
        via_g = (F[(i+1,j,a,g-1)] + 3*F[(i,j+1,a,g-1)])/4 - 2*F[(i,j,a,g-1)]
        assert value == via_a == via_g

assert len(F) == 46376
zeros = sorted(ix for ix,value in F.items() if value == 0)
assert zeros == [(0,0,0,1),(0,0,1,0),(0,0,1,2),
                 (0,0,2,1),(0,1,1,0),(1,0,0,1)]
assert all(value >= 0 for value in F.values())

positive = [(value/mu(i,j), (i,j,a,g), value)
            for (i,j,a,g),value in F.items() if value > 0]
min_ratio = min(row[0] for row in positive)
min_cases = sorted((ix,value) for ratio,ix,value in positive if ratio == min_ratio)
assert min_ratio == Q(1,5)
assert min_cases == [((0,2,1,0),Q(2)),((2,0,0,1),Q(2))]

# Check normalized base formulas throughout the available table.
for i in range(D):
    for j in range(D-i):
        den0 = (i+j+2)*(i+j+3)
        assert F[(i,j,1,0)]/mu(i,j) == Q(2*(i-j)**2+8*i**2-2*j+14*i,den0)
        assert F[(i,j,0,1)]/mu(i,j) == Q(2*(i-j)**2+8*j**2-2*i+14*j,den0)

# Ratio candidates. Delta-sign proposal:
# delta>0: R_m >= R and R_k <= R; delta<0: R_m <= R and R_k >= R.
cm, ck, dm, dk = [], [], [], []
edge_count = oriented_count = 0
for (i,j,a,g),value in F.items():
    n=a+g
    if i+j+n+1 > D:
        continue
    R = value/mu(i,j)
    Rm = F[(i+1,j,a,g)]/mu(i+1,j)
    Rk = F[(i,j+1,a,g)]/mu(i,j+1)
    edge_count += 1
    if Rm < R: cm.append(((i,j,a,g),R,Rm))
    if Rk < R: ck.append(((i,j,a,g),R,Rk))
    delta = a-g
    if delta:
        oriented_count += 1
        if (delta>0 and Rm<R) or (delta<0 and Rm>R):
            dm.append(((i,j,a,g),R,Rm))
        if (delta>0 and Rk>R) or (delta<0 and Rk<R):
            dk.append(((i,j,a,g),R,Rk))

key = lambda row: (sum(row[0]),row[0])
first_cm, first_ck = min(cm,key=key), min(ck,key=key)
first_dm, first_dk = min(dm,key=key), min(dk,key=key)
assert (edge_count,len(cm),len(ck)) == (40920,19332,19332)
assert (oriented_count,len(dm),len(dk)) == (38320,8817,8817)
assert first_cm == ((0,0,0,2),Q(3),Q(1))
assert first_ck == ((0,0,2,0),Q(3),Q(1))
assert first_dm == ((1,0,0,1),Q(0),Q(1,5))
assert first_dk == ((0,1,1,0),Q(0),Q(1,5))

# Adjacent minors, raw and mu-normalized.
raw, normalized = [], []
for (i,j,a,g),value in F.items():
    if i+j+a+g+2 > D:
        continue
    f10, f01, f11 = F[(i+1,j,a,g)], F[(i,j+1,a,g)], F[(i+1,j+1,a,g)]
    raw.append((value*f11-f10*f01,(i,j,a,g)))
    normalized.append(((value/mu(i,j))*(f11/mu(i+1,j+1))
                       -(f10/mu(i+1,j))*(f01/mu(i,j+1)),(i,j,a,g)))
raw_signs = (sum(z<0 for z,_ in raw),sum(z==0 for z,_ in raw),sum(z>0 for z,_ in raw))
norm_signs = (sum(z<0 for z,_ in normalized),sum(z==0 for z,_ in normalized),
              sum(z>0 for z,_ in normalized))
raw_zeros = sorted(ix for z,ix in raw if z==0)
first_norm_neg = min((row for row in normalized if row[0]<0),
                     key=lambda row:(sum(row[1]),row[1]))
first_norm_pos = min((row for row in normalized if row[0]>0),
                     key=lambda row:(sum(row[1]),row[1]))
assert len(raw)==35960 and raw_signs==(35958,2,0)
assert raw_zeros==[(0,0,0,1),(0,0,1,0)]
assert norm_signs==(11120,466,24374)
assert first_norm_neg==(Q(-1,5),(0,1,1,0))
assert first_norm_pos==(Q(3,5),(0,1,0,1))

constant_array_update = Q(3,4)+Q(1,4)-2
assert constant_array_update == -1

print('Catalan moment checks:',49)
print('symbolic base and distance-two identities: passed')
print('degree-60 table: states',len(F),'negative entries',0,'zeros',zeros)
print('minimum positive F/mu:',min_ratio,'at',min_cases)
print('ratio edges:',edge_count,'coordinate decreases:',len(cm),len(ck),
      'canonical:',first_cm,first_ck)
print('delta-oriented edges:',oriented_count,'failures:',len(dm),len(dk),
      'canonical:',first_dm,first_dk)
print('raw minor signs:',len(raw),raw_signs,'zero indices:',raw_zeros)
print('normalized minor signs:',norm_signs,'first +/-:',first_norm_neg,first_norm_pos)
print('constant-array plus update:',constant_array_update)

