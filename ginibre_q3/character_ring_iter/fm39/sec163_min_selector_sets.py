small=[253,1149,1277,1585,1653,1723,1731,1775,1791,2235,2303,
       2491,2499,2543,2551,2559,2815,3071,3327,3455,3567,3575,
       3583,3771,3775,3835,3839,3891,3959,4011,4019,4027,4031,
       4035,4071,4079,4087,4091,4095]
large=[2235,3835,1149,4078,2542,2542,2286,1132]
constraints=small+large

best=13
solutions=[]
for selectors in range(1,1<<12):
    if all(mask & selectors for mask in constraints):
        size=selectors.bit_count()
        if size<best:
            best=size
            solutions=[]
        if size==best:
            solutions.append([i for i in range(12) if selectors>>i&1])

print('minimum size',best,'number of sets',len(solutions))
print(solutions)
