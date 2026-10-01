import random, math
from fractions import Fraction
exec(open(__file__.replace('mech92_sweep_screen.py', 'mech92_charging_screen.py')).read().split("random.seed(1)")[0])
random.seed(2)
inreg = short = allpos = 0
for trial in range(6000):
    a = random.randint(2, 60); e = random.randint(2, 60); N = a + e
    B = random.randint(1, 10)
    p = random.randrange(1 if N % 2 else 2, N + 2*B + 2, 2)
    k = (N + p)//2; d = (N + 2*B - p)//2; t = min(a, e)
    if not ((N - 2*d)**2 < 4*(t-1)*(N-t+2)): continue   # region (10)
    if 2*(k - B) < N: continue                            # right-half windows only
    inreg += 1
    c = row(a, e)
    idx = list(range(k-B, k+B+2))
    # unwrapped sweep of psi_l along the window
    ang = []
    for l in idx:
        x, y = get(c, l), get(c, l-1) + get(c, l+1)
        ang.append(math.atan2(y, x))
    sw = 0.0
    for u, v in zip(ang, ang[1:]):
        dlt = v - u
        while dlt <= -math.pi: dlt += 2*math.pi
        while dlt > math.pi: dlt -= 2*math.pi
        sw += dlt
    if abs(sw) < math.pi - 1e-9: short += 1
    M = walks(k, B)
    if all(W(c, i, j) >= 0 for (i, j) in M if i < j): allpos += 1
print("region-(10) right-half cases", inreg, "| window sweep < pi:", short, "| all chords in window >= 0:", allpos)
