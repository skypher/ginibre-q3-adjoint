from datetime import datetime, timezone
from fractions import Fraction
import sys

if any(a in ("-h", "--help") for a in sys.argv[1:]):
    print("FM-STR7g exact verifier: checks variance identities, Lindeberg and two-scale constants.")
    raise SystemExit(0)

def stamp(msg):
    print(datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
          msg, flush=True)

stamp("start; Python integers and exact rational arithmetic")
assert Fraction(151,960) > Fraction(1,7)
assert Fraction(2,45)*Fraction(6,5) == Fraction(4,75)
assert Fraction(24,75)*Fraction(4,3) == Fraction(32,75)
assert Fraction(16,225)+Fraction(1,36) == Fraction(89,900)
assert Fraction(6,5)*252 == Fraction(1512,5)
assert Fraction(1512,1515) == Fraction(504,505)
stamp("sinc remainder constants and k=303 ratio constant checked")

for k in range(2,10001):
    V = k*(k+1)*(2*k+7)//6
    Q = k*(k+1)*(2*k+1)//6
    assert 6*V == k*(k+1)*(2*k+7)
    assert 3*V >= k**3
    assert V >= (k+1)**2
    assert k**3 <= 3*Q
    assert k*k*(k+2) <= 3*V
    assert k+2 <= 2*k
    if k >= 303:
        assert 1512 < 5*k
stamp("V_k bounds and both Lindeberg ratios checked exactly for k=2..10000")

A = (-8,-3,-1)
B = (-6,-4,-2)
pair = (-7,-5)
assert sum(map(abs,A)) == sum(map(abs,B)) == 12
assert len(A) == 3 and len(B)+len(pair) == 5
VA = sum(abs(n)*(abs(n)+2) for n in A)
VB = sum(abs(n)*(abs(n)+2) for n in B+pair)
assert VA == 98 and VB == 178
stamp("k=7,p=8 balanced cut: weights 12/12; minus counts 3/5; variances 98/178")

for k in (2,10,100,248,249,302,303,1000,10000):
    base = Fraction(1512,5*k)
    if k < 303:
        assert base > 1
    else:
        assert base < 1
    stamp(f"two-scale ratio base at k={k}: {base}")

stamp("FM-STR7g verifier PASS")
