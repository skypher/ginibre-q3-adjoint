from collections import defaultdict
from itertools import combinations
from math import comb
import random
import argparse

def cg(a, b):
    return range(abs(a-b), a+b+1, 2)

def coeff2(signed_labels, sx, sy):
    # Coefficient of U_sx(X) U_sy(Y) in the signed two-variable product.
    dp = {(0, 0): 1}
    for z in signed_labels:
        n, eps = abs(z), (1 if z > 0 else -1)
        nd = defaultdict(int)
        for (x, y), v in dp.items():
            for x2 in cg(x, n):
                nd[x2, y] += v
            for y2 in cg(y, n):
                nd[x, y2] += eps * v
        dp = nd
    return dp.get((sx, sy), 0)

def gx_vector(signed_labels):
    dp = {(0, 0): 1}
    for z in signed_labels:
        n, eps = abs(z), (1 if z > 0 else -1)
        nd = defaultdict(int)
        for (x, y), v in dp.items():
            for x2 in cg(x, n):
                nd[x2, y] += v
            for y2 in cg(y, n):
                nd[x, y2] += eps * v
        dp = nd
    return {x: v for (x, y), v in dp.items() if y == 0}

def triple_check(B, signs, p, ia, ib):
    # ia and ib are distinct factor positions with labels of equal parity.
    a, b = abs(B[ia]), abs(B[ib])
    sigma = 1
    for e in signs:
        sigma *= e
    signedB = [n*e for n,e in zip(B,signs)]
    C = [z for i,z in enumerate(signedB) if i not in (ia,ib)]
    lhs = coeff2(signedB, p, 0)
    gv = gx_vector(C)
    T1 = sum(gv.get(s, 0) for c in cg(a,b) for s in cg(p,c))

    rest_ab = C + [sigma*p]
    Dab = signs[ib] * coeff2(rest_ab, a, b)
    rest_ap = C + [signedB[ib]]
    Dap = sigma * coeff2(rest_ap, a, p)
    rest_bp = C + [signedB[ia]]
    Dbp = sigma * coeff2(rest_bp, b, p)
    rhs = T1 + (Dab + Dap + Dbp)//2
    assert (Dab + Dap + Dbp) % 2 == 0
    return lhs, rhs, (T1,Dab,Dap,Dbp)

def single_even_check(B, signs, p, ie):
    e = abs(B[ie])
    assert e % 2 == 0
    sigma = 1
    for z in signs:
        sigma *= z
    signedB = [n*s for n,s in zip(B,signs)]
    C = [z for i,z in enumerate(signedB) if i != ie]
    lhs = coeff2(signedB, p, 0)
    gv = gx_vector(C)
    main = sum(gv.get(s, 0) for s in cg(p,e))
    Dep = sigma * coeff2(C, e, p)
    return lhs, main + Dep, (main,Dep)

def invariant(ns):
    q = len(ns)
    if q == 0:
        return 1
    T = sum(ns)
    if T & 1 or q == 1:
        return 0
    H = T // 2
    z = 1 << q
    wt = [0]*z
    sz = [0]*z
    for S in range(1,z):
        bit = S & -S
        i = bit.bit_length()-1
        wt[S] = wt[S^bit] + ns[i]
        sz[S] = sz[S^bit] + 1
    ans = 0
    S = z-1
    J = S
    while True:
        top = H - wt[J] - sz[J] + q - 1
        if top >= q-2:
            ans += (-1 if sz[J] & 1 else 1) * comb(top,q-2)
        if J == 0:
            break
        J = (J-1)&S
    ans = -ans
    assert ans >= 0
    return ans

def multiplicities(ns):
    N = len(ns)
    z = 1 << N
    m = [0]*z
    m[0] = 1
    wt = [0]*z
    sz = [0]*z
    for S in range(1,z):
        bit = S & -S
        i = bit.bit_length()-1
        wt[S] = wt[S^bit] + ns[i]
        sz[S] = sz[S^bit] + 1
    for S in range(1,z):
        q = sz[S]
        if wt[S] & 1 or q == 1:
            continue
        H = wt[S]//2
        acc = 0
        J = S
        while True:
            top = H - wt[J] - sz[J] + q - 1
            if top >= q-2:
                acc += (-1 if sz[J] & 1 else 1) * comb(top,q-2)
            if J == 0:
                break
            J = (J-1)&S
        m[S] = -acc
        assert m[S] >= 0
    return m

def walsh(a):
    a = list(a)
    n = len(a)
    h = 1
    while h < n:
        for start in range(0,n,2*h):
            for j in range(start,start+h):
                x,y = a[j],a[j+h]
                a[j],a[j+h] = x+y,x-y
        h *= 2
    return a

def phi_table(ns):
    m = multiplicities(ns)
    full = (1 << len(ns))-1
    return walsh([m[S]*m[full^S] for S in range(full+1)])

def differences(values):
    row = list(values)
    out = []
    while row:
        out.append(row[0])
        row = [row[i+1]-row[i] for i in range(len(row)-1)]
    return out

def family_data(t):
    ns = tuple(2*(t+i) for i in range(12))
    mask = 6
    f = phi_table(ns)
    child_ids = [i for i in range(12) if i not in (9,10)]
    child_ns = tuple(ns[i] for i in child_ids)
    cmask = sum((1<<j) for j,i in enumerate(child_ids) if mask>>i & 1)
    fc = phi_table(child_ns)
    return ns, f, fc, cmask, mask

def main():
    rng = random.Random(20261002)
    NTEST = 3000
    triple_parity = [0,0]
    triple_rep = 0
    triple_equal = 0
    triple_plabel = 0
    for k in range(NTEST):
        p = rng.randint(3,8)
        n = rng.randint(3,7)
        B = [rng.randint(1,p) for _ in range(n)]
        mode = k % 4
        if mode == 0:
            v = rng.randint(1,p)
            B[0] = B[1] = v
            ia,ib = 0,1
        elif mode == 1:
            B[0] = p
            choices = [v for v in range(1,p+1) if v%2 == p%2]
            B[1] = rng.choice(choices)
            ia,ib = 0,1
        else:
            pairs = [(i,j) for i in range(n-1) for j in range(i+1,n-1)
                     if B[i]%2 == B[j]%2]
            if not pairs:
                B[0] = B[1]
                pairs = [(0,1)]
            ia,ib = rng.choice(pairs)
        if mode == 1:
            p = B[0]
            wanted = (p - sum(B[:-1])) & 1
            choices = [v for v in range(1,p+1) if v%2 == wanted]
            B[-1] = rng.choice(choices)
        else:
            p = max(3,max(B))
            if (sum(B)&1) != (p&1):
                p += 1
        assert p >= max(B) and ((sum(B)-p)&1)==0 and B[ia]%2 == B[ib]%2
        signs = [rng.choice((-1,1)) for _ in B]
        if (sum(s<0 for s in signs)&1) != (k&1):
            signs[-1] *= -1
        lhs,rhs,_ = triple_check(B,signs,p,ia,ib)
        assert lhs == rhs, ("triple",B,signs,p,ia,ib,lhs,rhs)
        parity = sum(s<0 for s in signs)&1
        triple_parity[parity] += 1
        triple_rep += int(len(set(B)) < len(B))
        triple_equal += int(abs(B[ia]) == abs(B[ib]))
        triple_plabel += int(p in (abs(B[ia]),abs(B[ib])))
        if (k+1)%500 == 0:
            print("triple_identity_progress",k+1,flush=True)

    single_parity = [0,0]
    single_repeated = 0
    for k in range(NTEST):
        p0 = rng.randint(3,9)
        n = rng.randint(3,8)
        B = [rng.randint(1,p0) for _ in range(n)]
        B[0] = rng.choice([v for v in range(2,p0+1,2)])
        p = max(3,max(B))
        if (sum(B)&1) != (p&1):
            p += 1
        signs = [rng.choice((-1,1)) for _ in B]
        if (sum(s<0 for s in signs)&1) != (k&1):
            signs[-1] *= -1
        lhs,rhs,_ = single_even_check(B,signs,p,0)
        assert lhs == rhs, ("single-even",B,signs,p,lhs,rhs)
        single_parity[sum(s<0 for s in signs)&1] += 1
        single_repeated += int(len(set(B)) < len(B))
        if (k+1)%500 == 0:
            print("single_even_progress",k+1,flush=True)

    print("IDENTITY_TESTS triple=",NTEST,"single_even=",NTEST,
          "triple_minus_parity=",triple_parity,
          "triple_repeated_lists=",triple_rep,
          "triple_equal_pair_labels=",triple_equal,
          "triple_pair_contains_p_label=",triple_plabel,
          "single_minus_parity=",single_parity,
          "single_repeated_lists=",single_repeated,flush=True)

    # Finite 12-factor witness; mask 5 marks labels 40 and 44 negative.
    ns = tuple(range(40,64,2))
    mask = 5
    tab = phi_table(ns)
    parent = tab[mask]
    pairs = list(combinations(range(12),2))
    D = []
    for i,j in pairs:
        diff = parent - tab[mask^(1<<i)^(1<<j)]
        assert diff % 4 == 0
        D.append((diff//4,ns[i],ns[j]))
    assert len(D)==66 and all(d<0 for d,_,_ in D)
    assert min(D)==(-4506106071735,40,44)
    assert max(D)==(-34349665,42,44)
    assert parent == 906415068445728
    tp = (9,10)
    child_ids = [i for i in range(12) if i not in tp]
    child_ns = tuple(ns[i] for i in child_ids)
    child_mask = sum((1<<j) for j,i in enumerate(child_ids) if mask>>i & 1)
    child = phi_table(child_ns)[child_mask]
    assert child == 359292699210
    W = sum(ns[:-1])
    delta = (W-ns[-1])//2
    wtp = ns[9]+ns[10]
    assert (W,delta,wtp)==(550,244,118)
    assert 2*wtp < delta
    print("FINITE_FSTAR parent_Phi=",parent,"g=",parent//2,
          "flip_D_min=",min(D),"flip_D_max=",max(D),
          "flip_D_max_pair=",max(D)[1:],
          "TopPair_indices=",tp,"child_Phi=",child,"child_g=",child//2,
          "drop_g=",parent//2-child//2,"W_delta_wTP=",W,delta,wtp,flush=True)

    # Exact Newton/binomial-basis certificates for the infinite family.
    t0 = 1024
    pair_values = [[] for _ in pairs]
    children, drops = [], []
    for dt in range(11):
        t = t0+dt
        ns_t, ft, fc, cm, mm = family_data(t)
        children.append(fc[cm]//2)
        drops.append((ft[mm]-fc[cm])//2)
        for col,(i,j) in enumerate(pairs):
            inc = ft[mm^(1<<i)^(1<<j)] - ft[mm]
            assert inc%4 == 0
            pair_values[col].append(inc//4)
    certs = [differences(vals[:9]) for vals in pair_values]
    assert all(len(c)==9 and c[7]==0 and c[8]==0 for c in certs)
    assert all(c[0]>0 and all(x>=0 for x in c) for c in certs)
    assert certs[0]==[5822451256163514,31047417002688,131148758741,
                      412024173,856711,885,0,0,0]
    assert min(c[0] for c in certs)==5822451256163514
    assert min(x for c in certs for x in c if x>0)==885

    child_cert = differences(children)
    drop_cert = differences(drops)
    assert len(child_cert)==11 and child_cert[9]==child_cert[10]==0
    assert child_cert == [38052642467178680932469,259571825669558373387,
                          1516210041603629293,7373203649067419,28656277938588,
                          83449110560,161848960,156800,0,0,0]
    assert len(drop_cert)==11 and drop_cert[0]>0 and all(x>=0 for x in drop_cert)
    assert drop_cert == [125092638475790757878824532193,
                         1097318339201163864759187589,
                         8547855907644311898389392,58205819393843972766644,
                         339395322276535612098,1647560924888071532,
                         6392120371066728,18581806084352,35976450048,
                         34793472,0]
    print("FAMILY flip_polys=",len(certs),"Newton_coefficients=",sum(map(len,certs)),
          "min_constant=",min(c[0] for c in certs),
          "min_positive_coefficient=",min(x for c in certs for x in c if x>0),flush=True)
    print("FAMILY child_g_Newton=",child_cert,"positive_for_t_ge_1024=True",flush=True)
    print("FAMILY TopPair_drop_g_Newton=",drop_cert,"positive_for_t_ge_1024=True",flush=True)

    # Extra exact evaluations of the family polynomials.
    extra = 0
    for t in (t0+9,t0+17,t0+80):
        ns_t, ft, fc, cm, mm = family_data(t)
        for col,(i,j) in enumerate(pairs):
            got = (ft[mm^(1<<i)^(1<<j)]-ft[mm])//4
            pred = sum(c*comb(t-t0,h) for h,c in enumerate(certs[col]))
            assert got == pred
            extra += 1
        child_g = fc[cm]//2
        drop_g = (ft[mm]-fc[cm])//2
        assert child_g == sum(c*comb(t-t0,h) for h,c in enumerate(child_cert))
        assert drop_g == sum(c*comb(t-t0,h) for h,c in enumerate(drop_cert))
    print("FAMILY_extra_flip_polynomial_checks=",extra,flush=True)

    # TopPair census failure at ratio 3/7, checked with direct two-variable fusion.
    B = [-1]*13 + [-2] + [-3]*5 + [-4]
    p = 6
    childB = B.copy()
    childB.remove(-2); childB.remove(-4)
    parent_g = coeff2(B,p,0)
    child_g = coeff2(childB,p,0)
    assert parent_g == 700607800 and child_g == 723462700
    assert len(B)==20 and sum(abs(z) for z in B)==34
    assert (sum(abs(z) for z in B)-p)//2==14
    assert (sum(abs(z) for z in B if z%2==0)) == 6
    print("TOPPAIR_3_over_7 parent_g=",parent_g,"child_g=",child_g,
          "delta=14","wTP=6","ratio=3/7","parent_Phi=",2*parent_g,
          "child_Phi=",2*child_g,flush=True)

if __name__ == "__main__":
    argparse.ArgumentParser(description="FM-CHK98 exact checks").parse_args()
    main()
