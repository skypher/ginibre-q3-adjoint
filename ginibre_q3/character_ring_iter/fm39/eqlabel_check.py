# Independent check: phi_r(h_u^2 h_w h_1^a) >= 0 for 2u >= a+2r+w-2 (G0, equal largest labels), definition-level evaluator.
exec(open('split3.py').read())
n=bad=0
for r in range(2,12):
    for a in range(0,50):
        for w in range(1,20):
            for u in range(w,70):
                if 2*u < a+2*r+w-2 or (a+2*u+w)%2: continue
                v=phi3(r,a,u,u,w); n+=1; bad+= v<0
print('equal-largest-label G0 words checked:',n,' negative:',bad)
