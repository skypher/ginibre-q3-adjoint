exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
import itertools
for off in (-1,1):
    neg=0; n=0
    for r in range(2,11):
        a=2*r-3+off
        for kap in itertools.combinations_with_replacement(range(1,23),3):
            if (a+sum(kap))%2: continue
            v=phi_kernel(core(kap,()),r,a); n+=1; neg+=v<0
    print(f'strip a=2r-3{off:+d} (r<=10, labels<=22): {n} values, negatives {neg}',flush=True)
