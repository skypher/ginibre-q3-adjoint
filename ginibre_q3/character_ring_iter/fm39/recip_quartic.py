exec(open('recip_bounds.py').read().split("for trial in range(600):")[0])
random.seed(13)
for trial in range(600):
    e=random.choice([1,3,5]); a=random.randint(0,5)
    c=mkrow(a,e,[])
    for _ in range(random.randint(1,2)):
        # complex rho = R e^{i th} off the unit circle; quartic (z^2-2R cos z + R^2)(R^2 z^2 - 2 R cos z + 1), rational params
        R=Fr(random.randint(2,9),random.randint(1,3)) if random.random()<.5 else Fr(random.randint(1,3),random.randint(4,9))
        co=Fr(random.randint(-9,9),10)   # cos(theta)
        q1=[R*R,-2*R*co,1]; q2=[1,-2*R*co,R*R]
        c=mulrow(mulrow(c,q1),q2)
    three_all(c,'three-factor, reciprocal with complex root quadruples (positive weight, not real-rooted)',(a,e))
for k,v in stats.items(): print(k,': tests',v[0],'failures',v[1],'rows with a failure',v[2],'of 600; first',v[3])
