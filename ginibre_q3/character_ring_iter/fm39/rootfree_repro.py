
from fractions import Fraction as Q
from math import comb, factorial
import argparse

def fall(n,m):
    return factorial(n)//factorial(n-m)

def kvals(n,x,m):
    K=[1]
    if m:
        K.append(x)
    for l in range(1,m):
        K.append(x*K[-1]-l*(n-l+1)*K[-2])
    return K

def coeff(a,e):
    return [sum((-1)**h*comb(e,h)*comb(a,k-h)
                for h in range(max(0,k-a),min(e,k)+1))
            for k in range(a+e+1)]

def variation(v):
    signs=[1 if x>0 else -1 for x in v if x]
    return sum(x!=y for x,y in zip(signs,signs[1:]))

def setup(a,e,lo,hi):
    N=a+e
    n=N+2
    c=coeff(a,e)
    C=lambda k:c[k] if 0<=k<=N else 0
    D=lambda k:C(k)**2-C(k-1)*C(k+1)
    B=lambda k:C(k-1)+C(k+1)
    degrees=range(e%2,e+1,2)
    norm={l:factorial(l)*fall(n,l) for l in degrees}
    eta=Q(factorial(e),4*fall(n,e+2))
    K={k:kvals(n,2*k-N,e+2) for k in range(lo,hi+1)}
    b={k:comb(n,k+1) for k in K}
    z={k:(2*k-N)**2 for k in K}

    def pair(j,i):
        rows=[]
        for l in degrees:
            A=sum((z[k+1]-z[k])*b[k]*b[k+1]
                  *K[k][l]*K[k+1][l]
                  for k in range(j,i))
            R=(z[i]-z[j])*b[j]*b[i]*K[j][l]*K[i][l]
            rows.append((l,Q(A,norm[l]),Q(R,norm[l])))
        T=D(j)-D(i)
        W=B(i)*C(j)-C(i)*B(j)
        assert T==eta*sum(A for l,A,R in rows)
        assert W==eta*sum(R for l,A,R in rows)
        assert T>=abs(W)
        return T,W,rows

    return n,K,z,pair

def main():
    argparse.ArgumentParser(
        description="Exact two-point kernel, root-free theorem and root bounds."
    ).parse_args()

    identities=rootfree=components=0
    for e in range(3,9):
        for a in range(e+2,e+8):
            N=a+e
            n,K,z,pair=setup(a,e,(N+1)//2,N+1)
            for j in K:
                for i in range(j+1,N+2):
                    T,W,rows=pair(j,i)
                    all_good=True
                    for l,A,R in rows:
                        good=(
                            K[j][l]!=0 and K[i][l]!=0
                            and variation(K[j][:l+1])
                            ==variation(K[i][:l+1])
                        )
                        if good:
                            assert A>=R>=0
                            components+=1
                        all_good &= good
                    if all_good:
                        rootfree+=1
                        assert W>=0 and T>=W
                    identities+=1

    regions={"central_even":0,"central_odd":0,"outer":0}
    for e in range(3,9):
        for a in (200,600,1200):
            N=a+e
            central_start=(N+1)//2+1
            windows=[(central_start,central_start+6),
                     (N-10,N+1)]
            for lo,hi in windows:
                n,K,z,pair=setup(a,e,lo,hi)
                for j in range(lo,hi):
                    for i in range(j+2,hi+1):
                        central_bound=(
                            Q(2*(n-e+2),e) if e%2==0
                            else Q(4*(n-e+3),e-1)
                        )
                        central=z[i]<=central_bound
                        outer=z[j]>=4*(e-1)*(n-e+2)
                        if not (central or outer):
                            continue
                        T,W,rows=pair(j,i)
                        for l,A,R in rows:
                            assert A>=R>=0
                            assert variation(K[j][:l+1]) \
                                ==variation(K[i][:l+1])
                        if central:
                            key="central_odd" if e%2 else "central_even"
                            regions[key]+=1
                        if outer:
                            regions["outer"]+=1

    n,K,z,pair=setup(8,6,8,13)
    T,W,rows=pair(8,13)
    eta=Q(factorial(6),4*fall(n,8))
    defects=[(l,eta*(A-R)) for l,A,R in rows]
    assert (T,W)==(1216,40)
    assert dict(defects)[4]==Q(-432,11)
    assert sum(value for l,value in defects)==1176

    print("Kernel and telescoping identities:",identities)
    print("Common root-free pairs:",rootfree)
    print("Root-free component inequalities:",components)
    print("Explicit root-bound regions:",regions)
    print("Scaled component defects at (a,e,j,i)=(8,6,8,13):",defects)
    print("Same point, T and W:",T,W)
    print("All assertions passed.")

if __name__=="__main__":
    main()
