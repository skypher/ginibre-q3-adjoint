"""FM-MECH35 (astra_max_ceres), extracted from its run log: mech35_hAC_d3_budget.py"""
import argparse
import sympy as s
argparse.ArgumentParser(description='Display the positive coefficient polynomials used in the uniform proof').parse_args()
t,b,u,T,B=s.symbols('t b u T B')
C=lambda n,k:s.prod(n-i for i in range(k))/s.factorial(k)
alpha=(t+4*b-8)/12; h=(t+4*b-10)/4
for v in (0,1):
    gamma=(6*b*b+4*b*t-14*b+t*t-7*t-6*v*v-12*v+10)/12
    D=C(t+b+1,3)-t*(t+b-1)-v
    beta=h*h/(4*alpha) if v else 0
    rho=s.factor(D-C(t,3)/4-s.Rational(3,4)*t*v-s.Rational(4,3)*u*u-u/2-alpha*C(t,2)-beta-gamma*t/2)
    if v==0:
        F=s.expand(24*rho.subs(u,b).subs({t:T+4,b:B+1}))
        print('v=0, 24 rho lower bound:',s.collect(F,T))
    else:
        numerator=s.factor(48*(t+4*b-8)*rho)
        F=s.expand(numerator.subs(u,b-1).subs({t:T+4,b:B+2}))
        print('v=1, denominator 48(t+4b-8), numerator lower bound:')
        for k in range(4,-1,-1): print('T^'+str(k),s.expand(F).coeff(T,k))
        edge=s.factor(rho.subs({b:1,u:0}))
        print('v=1,b=1,rho:',edge)
        print('edge numerator after t=T+6:',s.expand(s.fraction(edge)[0].subs(t,T+6)))