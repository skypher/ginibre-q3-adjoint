"""FM-MECH50: exact uniform certificates and Catalan bridges; no files."""
import argparse
from fractions import Fraction
from math import comb
import sympy as S

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--box", type=int, default=16,
                    help="bounded Catalan bridge box (default 16)")
    args = ap.parse_args()
    n,k,d,r,z,u,v,w,N,t,D = S.symbols("n k d r z u v w N t D")
    j = n+k
    H = (n+1)*(j+2)**2*(j+3)**2*(j+4)
    c = {-1:(d-(j+1)*r)/(n+1), 0:S.Integer(1), 1:r}
    for h in range(1,4):
        c[h+1] = S.cancel((d*c[h]-(n-h+1)*c[h-1])/(j+h+1))
    td = lambda h: c[h]**2-c[h-1]*c[h+1]
    slack = td(0)-td(3)+(c[2]+c[4])-c[3]*(c[-1]+c[1])
    Q = S.Poly(S.cancel(H*slack),r)
    A,B,C = [Q.nth(h) for h in range(3)]
    alpha = (d*d*n*(j*j+6*k+7*n+9)
             +(j+3)**2*(n+1)*(j*j+6*j+n*n-2*n+8))
    assert S.expand(A-(k+2)*alpha)==0
    beta = S.cancel(-B/d)
    assert S.Poly(beta,d).degree()==2
    print("recurrence quadratic: PASS", flush=True)

    def positive(poly, variables, label):
        pp = S.Poly(S.expand(poly),*variables)
        coeff = pp.coeffs()
        assert all(c>0 for c in coeff), label
        assert pp.coeff_monomial(1)>0, label
        print(label, "terms",len(coeff),"minimum",min(coeff),flush=True)

    def bernstein(poly, var):
        pp = S.Poly(S.expand(poly),var)
        m = pp.degree()
        return [S.expand(sum(pp.nth(i)*S.binomial(h,i)/S.binomial(m,i)
                             for i in range(h+1))) for h in range(m+1)]

    def even_sub(poly, value):
        pp = S.Poly(poly,d)
        assert all(m[0]%2==0 for m,c in pp.terms())
        return S.cancel(sum(c*value**(m[0]//2) for m,c in pp.terms()))

    # Coverage. J is the reduced discriminant as a cubic in D=d^2.
    reduced = S.cancel((B*B-4*A*C)/((j+2)**2*(j+3)**2))
    J = S.Poly(even_sub(reduced,D),D)
    assert J.degree()==3
    assert S.expand(J.nth(3)-(k+3)**2)==0
    a2 = 2*(k**4+12*k**3+49*k**2+2*k*n+80*k+2*n**2+10*n+44)
    assert S.expand(J.nth(2)-a2)==0
    a0 = (-4*(k+2)*(k+4)*(n+1)*(k+n+4)
          *(k*k+2*k*n+4*k+2*n*n+4*n+2)
          *(k*k+2*k*n+6*k+2*n*n+4*n+8))
    assert S.expand(J.nth(0)-a0)==0
    O = (n-1)*(j+2)
    P = S.Poly(S.expand(-J.as_expr().subs(D,4*O)/4),n)
    assert P.degree()==6 and P.nth(6)==-36
    assert S.expand(P.nth(5)-72*(k+1)**2)==180*k+300
    assert all(c>0 for h in range(5)
               for c in S.Poly(P.nth(h),k).all_coeffs())
    L = 3*(k+1)**2
    dc = (2*n+k)**2*(L-2)**2
    cover = S.expand(-sum(c*dc**m[0]*L**(6-2*m[0])
                          for m,c in J.terms()))
    positive(cover.subs(n,2*(k+1)**2+u).subs(k,v+5),
             (u,v),"coverage endpoint")
    print("coverage convexity and outer endpoint: PASS",flush=True)

    # Outer region. z=O R^2/n^2 belongs to (0,1].
    dd = O*(1+z)**2/z
    az,bz,cz = [even_sub(p,dd) for p in (A,beta,C)]
    middle = S.cancel(2*z*(az-n*(1+z)*bz/2))
    endpoint = S.cancel((n-1)*(az-n*(1+z)*bz+cz*n*n*z/O))
    for name,poly in (("outer middle",middle),("outer endpoint",endpoint)):
        assert S.Poly(poly,z).degree()==3
        for h,bv in enumerate(bernstein(poly,z)):
            positive(bv.subs({n:u+3,k:v+5}), (u,v),name+" "+str(h))

    # Central region. Bernstein coefficients on the OL ratio interval.
    qc = [S.expand(p.subs({n:(N-k)/2,d:N-2*t})) for p in (A,B,C)]
    U = (k+2)/k
    lo_odd = U*(1-(4*t+2)*(k+1)/(3*N))
    for parity,lo,hi,tmin in (("even",S.Integer(0),S.Integer(1),4),
                              ("odd",lo_odd,U,3)):
        aa,bb,cc = qc
        bv = [aa+bb*lo+cc*lo*lo,
              aa+bb*(lo+hi)/2+cc*lo*hi,
              aa+bb*hi+cc*hi*hi]
        for h,expr in enumerate(bv):
            num,den = S.fraction(S.cancel(expr))
            assert den in (1,36*N**2*k**2,24*N*k**2,4*k**2)
            shifted = S.expand(num.subs(N,3*t*(k+1)**2+w))
            positive(shifted.subs({t:u+tmin,k:v+5}),(u,v,w),
                     "central "+parity+" "+str(h))

    # Independent definition-level Catalan evaluation.
    def moment(h):
        return 0 if h<0 or h%2 else comb(h,h//2)//(h//2+1)

    def row(a,e):
        cc=[1]
        for sign,times in ((1,a),(-1,e)):
            for _ in range(times):
                cc=[(cc[h] if h<len(cc) else 0)
                    +sign*(cc[h-1] if h else 0)
                    for h in range(len(cc)+1)]
        return cc

    def row_slack(cc,jj):
        get=lambda h: cc[h] if 0<=h<len(cc) else 0
        det=lambda h:get(h)**2-get(h-1)*get(h+1)
        bb=lambda h:get(h-1)+get(h+1)
        return det(jj)-det(jj+3)+bb(jj+3)*get(jj)-get(jj+3)*bb(jj)

    def direct(a,e,p,cc):
        NN=a+e
        value=0
        for h,ch in enumerate(cc):
            for ell in range(p//2+1):
                xp=NN-h+p-2*ell
                up=(-1)**ell*comb(p-ell,ell)
                value += ch*up*(moment(xp+2)*moment(h)
                         +moment(xp)*moment(h+2)
                         -2*moment(xp)*moment(h))
        return value

    counts=[0,0]
    for e in range(1,args.box+1):
        for a in range(args.box+1):
            NN=a+e
            cc=row(a,e)
            for p in range(7,args.box+1):
                if (NN+p)%2 or NN+2-p<6:
                    continue
                jj=(NN+p-2)//2
                value=row_slack(cc,jj)
                assert value==direct(a,e,p,cc)>=0
                assert value==row_slack([(-1)**h*x
                                        for h,x in enumerate(cc)],jj)
                counts[e%2]+=1
    print("Catalan bridges (even e, odd e):",counts,flush=True)

    # Fixed tests exercise both central parities and the outer bound.
    def coefficient(a,e,h):
        if h<0 or h>a+e:return 0
        return sum((-1)**ell*comb(e,ell)*comb(a,h-ell)
                   for ell in range(max(0,h-a),min(e,h)+1))

    points=((3,30,18),(3,2575,9),(4,1729,7),(5,4000,10))
    regions=[]
    for tt,NN,kk in points:
        assert (NN-kk)%2==0
        nn=(NN-kk)//2; jj=nn+kk; dd=NN-2*tt
        cv={h:coefficient(NN-tt,tt,jj+h) for h in range(-1,5)}
        rr=Fraction(cv[1],cv[0])
        OO=(nn-1)*(jj+2)
        if dd*dd>=4*OO:
            assert rr>0
            assert OO*rr*rr-nn*dd*rr+nn*nn>=0
            assert rr<=Fraction(nn*dd,2*OO)
            regions.append("outer")
        else:
            assert NN>=3*tt*(kk+1)**2
            if tt%2:
                upper=Fraction(kk+2,kk)
                lower=upper*(1-Fraction((4*tt+2)*(kk+1),3*NN))
                assert lower<=rr<=upper
                regions.append("central odd")
            else:
                assert 0<rr<=1
                regions.append("central even")
        det=lambda h:cv[h]**2-cv[h-1]*cv[h+1]
        value=det(0)-det(3)+(cv[2]+cv[4])*cv[0]-cv[3]*(cv[-1]+cv[1])
        assert value>0
    print("fixed ratio tests:",regions,flush=True)
    print("PASS",flush=True)

if __name__=="__main__":
    main()
