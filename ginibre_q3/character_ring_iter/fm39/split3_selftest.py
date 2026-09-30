exec(open('fm3kern.py').read().split("if __name__=='__main__'")[0])
exec(open('split3.py').read())
import random; random.seed(5); n=0
for _ in range(300):
    r=random.randint(2,5); a=random.randint(0,8); u,v,w=sorted([random.randint(1,7) for _ in range(3)],reverse=True)
    if (a+u+v+w)%2: continue
    assert phi3(r,a,u,v,w)==phi_kernel(core((u,v,w),()),r,a),(r,a,u,v,w); n+=1
print('split3 vs fm3kern.phi_kernel agree on',n,'cases')
