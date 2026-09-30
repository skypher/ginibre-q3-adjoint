import argparse
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial

parser = argparse.ArgumentParser(description='Exact FM-SEC125 moment checks.')
parser.parse_args()

def rising(a, n):
    out = Q(1)
    for j in range(n):
        out *= a + j
    return out

@lru_cache(None)
def mu(m, k):
    return Q(
        2 * factorial(2*m) * factorial(2*m+1)
          * factorial(2*k) * factorial(2*k+1),
        factorial(m)**2 * factorial(k)**2
          * factorial(m+k+1) * factorial(m+k+2),
    )

@lru_cache(None)
def M(A, E, b):
    m, k = A//2, E//2
    return sum((
        Q(comb(b,u) * (-2)**(b-u), 2**u)
        * comb(u,v) * mu(m+v, k+u-v)
        for u in range(b+1) for v in range(u+1)
    ), Q(0))

@lru_cache(None)
def R2(A, E, b):
    m, k = A//2, E//2
    return sum((
        Q(comb(b,u) * (-2)**(b-u), 2**u) * comb(u,v)
        * (mu(m+v+2, k+u-v)
           - 2*mu(m+v+1, k+u-v+1)
           + mu(m+v, k+u-v+2)) / 16
        for u in range(b+1) for v in range(u+1)
    ), Q(0))

def H(m, k):
    return (mu(m+1,k) + mu(m,k+1))/2 - 3*mu(m,k)

def G(m, k):
    return (
        mu(m+2,k) + 6*mu(m+1,k+1) + mu(m,k+2)
        - 12*mu(m+1,k) - 12*mu(m,k+1) + 16*mu(m,k)
    ) / 8

def h3_blocks(A, E, b):
    m, k = A//2, E//2
    return sum((
        Q(comb(b,u)*(-2)**(b-u), 2**u) * comb(u,v)
        * H(m+v, k+u-v)
        for u in range(b+1) for v in range(u+1)
    ), Q(0))

def s4_blocks(A, E, b):
    m, k = A//2, E//2
    return sum((
        Q(comb(b,u)*(-2)**(b-u), 2**u) * comb(u,v)
        * G(m+v, k+u-v)
        for u in range(b+1) for v in range(u+1)
    ), Q(0))

def cert_h3(m, k, b):
    ell, eta = min(m,k), abs(m-k)
    n = m+k
    if b == 0:
        return mu(m,k) * Q(n*n + 4*eta*eta + n - 6,
                           (n+2)*(n+3))
    if b == 1:
        p = (32*ell**4 + 96*ell**3 + 136*ell**2 + 192*ell + 144
             + eta*(64*ell**3 + 144*ell**2 + 136*ell + 96)
             + eta**2*(96*ell**2 + 240*ell + 114)
             + eta**3*(64*ell + 96) + 30*eta**4)
        return mu(m,k) * Q(p, (n+2)*(n+3)**2*(n+4))
    raise ValueError('certificate only covers b=0,1')

def cert_s4(m, k, b):
    ell, eta = min(m,k), abs(m-k)
    n = m+k
    if b == 0:
        p = (48*ell**4 + 96*ell**3*eta + 144*ell**3
             + 48*ell**2*eta**2 + 216*ell**2*eta + 108*ell**2
             + 72*ell*eta**2 + 108*ell*eta
             + 5*eta**4 - 5*eta**2)
        return mu(m,k) * Q(2*p, (n+2)*(n+3)**2*(n+4))
    if b == 1:
        p = (768*ell**6 + 4992*ell**5 + 12480*ell**4 + 16032*ell**3
             + 11808*ell**2 + 4320*ell
             + eta*(2304*ell**5 + 12480*ell**4 + 24960*ell**3
                    + 24048*ell**2 + 11808*ell + 2160)
             + eta**2*(2880*ell**4 + 11712*ell**3 + 16672*ell**2
                       + 10296*ell + 2760)
             + eta**3*(1920*ell**3 + 5088*ell**2 + 4192*ell + 1140)
             + eta**4*(656*ell**2 + 1464*ell + 780)
             + eta**5*(80*ell + 300) + 60*eta**6)
        return mu(m,k) * Q(
            p, (n+2)*(n+3)**2*(n+4)**2*(n+5)
        )
    raise ValueError('certificate only covers b=0,1')

h_count = s_count = 0
h_min = s_min = None
h_zero, s_zero = [], []

for A in range(2, 21, 2):
    for E in range(0, 21, 2):
        for b in range(21):
            value = M(A,E,b+1) - M(A,E,b)
            assert value == h3_blocks(A,E,b)
            assert value >= 0
            if b < 2:
                assert value == cert_h3(A//2,E//2,b)
            h_count += 1
            h_min = (value,(A,E,b)) if h_min is None or value < h_min[0] else h_min
            if value == 0:
                h_zero.append((A,E,b))

for A in range(0, 21, 2):
    for E in range(0, 21, 2):
        for b in range(21):
            value = M(A,E,b+2) + M(A,E,b+1) - 2*R2(A,E,b)
            assert value == s4_blocks(A,E,b)
            assert value >= 0
            if b < 2:
                assert value == cert_s4(A//2,E//2,b)
            s_count += 1
            s_min = (value,(A,E,b)) if s_min is None or value < s_min[0] else s_min
            if value == 0:
                s_zero.append((A,E,b))

h_neg, h_pos = -2*H(1,2), (H(2,2)+H(1,3))/2
s_neg, s_pos = -2*G(1,1), (G(2,1)+G(1,2))/2
assert (h_neg,h_pos,h_neg+h_pos) == (Q(-4),Q(12),Q(8))
assert (s_neg,s_pos,s_neg+s_pos) == (Q(-4),Q(8),Q(4))

ray = 4*rising(Q(4),1)/rising(Q(11,2),1) - 3
assert ray == Q(-1,11)

print('h3 exact cases:', h_count, 'negative: 0', 'minimum:', h_min)
print('h3 zero cases:', h_zero)
print('S4 exact cases:', s_count, 'negative: 0', 'minimum:', s_min)
print('S4 zero cases:', s_zero)
print('first h3 alternating blocks:', h_neg, h_pos, 'sum:', h_neg+h_pos)
print('first S4 alternating blocks:', s_neg, s_pos, 'sum:', s_neg+s_pos)
print('h3 t=1 ray value:', ray)
