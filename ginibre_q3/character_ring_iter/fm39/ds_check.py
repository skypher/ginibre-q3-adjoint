# Independent grid check of Theorem DS (a = 2r-3, all labels) and the adjacent-strip branch, with the kernel evaluator.
exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
import itertools
neg=0; n=0; zero=0
for r in range(2,12):
    a=2*r-3
    for kap in itertools.combinations_with_replacement(range(1,26),3):
        if (a+sum(kap))%2: continue
        v=phi_kernel(core(kap,()),r,a); n+=1; neg+=v<0; zero+=v==0
print('DS grid (r<=11, labels<=25):',n,'values, negatives',neg,'zeros',zero)
neg=0; n=0
for r in range(2,10):
    a=2*r-2
    for u,v in itertools.combinations_with_replacement(range(1,22,2),2):
        for w in range(2,24,2):
            val=phi_kernel(core((u,v,w),()),r,a); n+=1; neg+=val<0
print('adjacent strip a=2r-2, two odd labels:',n,'values, negatives',neg)
