# Complete two-label `B_1/B_5` positivity in the rank-six orbit ring

Date: 2026-07-25

Let `O_11` have orbit basis `B_0,...,B_5`, and write

```text
D_a=B_a tensor B_0-B_0 tensor B_a,
S_a=B_a tensor B_0+B_0 tensor B_a.
```

## 1. `B_1` in the minus sector

Put `T=B_5`.  The fusion identity

```text
T^2=B_0+B_1
```

gives `D_1=D_T S_T`.  Consequently every even-minus word supported on
`D_1,D_5,S_5` reduces identically to an even-minus top-ray word and is
nonnegative by the all-length top-orbit-ray theorem.

## 2. `B_1` in the plus sector

The only nontrivial support-disjoint case left is

```text
S_1^u D_5^(2v),                    u>=1, v>=1.          (2.1)
```

If `u` is even, the spectral integrand in (2.1) is pointwise nonnegative.
For `u=1+2q`, `v=1+p`, the unordered-pair spectral expansion has residual
bases

```text
lambda_1(i,j)=(B_5(i)-B_5(j))^2,
lambda_2(i,j)=(B_1(i)+B_1(j))^2.
```

Only three of the fifteen unordered spectral pairs have negative floor
coefficient.  Exact rational Sturm isolation and interval arithmetic show
that every negative pair has positive dominant neighbours in both residual
bases.  The seven capacitated Hall inequalities, one for every nonempty
subset of the three negative pairs, all hold.  Exponential transport then
proves (2.1) for every `p,q>=0`.

The exact replay is

```text
character_ring_iter/verify_su2_o11_b1plus_topminus_transport.cpp.
```

It reports

```text
SU2_O11_B1PLUS_TOPMINUS_TRANSPORT PASS
roots=6 pairs=15 negative_pairs=3 hall_subsets=7 residual_dimensions=2
```

## 3. Theorem

> **Theorem 3.1.** Every support-disjoint even-minus signed word in `O_11`
> whose active nontrivial labels are contained in `{B_1,B_5}` has
> nonnegative scalar corner, for arbitrary exponents.

Overlapping signs are removed by the standard support-overlap reduction, so
the same conclusion holds without the support-disjoint qualification.

This removes the complete extreme two-label support from the unresolved
rank-six census.  It also shows that the former first chamber was only one
face of a larger proved family.
