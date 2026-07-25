# An all-length `V_1/V_2` minus-sector theorem for ordinary `SU(2)`

Date: 2026-07-25

## 1. Finite odd-level theorem

Let `k=2m-1>=3` be odd.  In the doubled `SU(2)_k` fusion ring write

```text
D_a=V_a tensor V_0-V_0 tensor V_a,
S_a=V_a tensor V_0+V_0 tensor V_a.
```

Consider an even-minus signed word whose labels lie in

```text
{1,2,k-1},
```

with every occurrence of label `2` carrying minus sign.  The signs of the
`V_1` and `V_(k-1)` factors are arbitrary.

Under the odd simple-current lift,

```text
V_2       -> B_1,
V_1       -> T=B_(m-1),
V_(k-1)   -> T=B_(m-1).
```

The `V_2` signs are unchanged in both orbit terms.  The signs of the odd
label `V_1` are reversed in the second orbit term, while the even label
`V_(k-1)` keeps its signs.

If the number of `V_1` occurrences is odd, the full finite corner vanishes by
the odd-label selection rule.  If it is even, both orbit terms have even
minus parity.  Each is supported on `D_1,D_T,S_T`; the identity

```text
D_1=D_T S_T
```

reduces each term to an even-minus top-orbit-ray word.  The all-length top-ray
theorem proves both terms nonnegative.

> **Theorem 1.1.**  At every odd level `k>=3`, every even-minus signed word
> supported on `V_1,V_2,V_(k-1)`, with all `V_2` factors in the minus sector,
> has nonnegative scalar corner.

The same argument after adjoining one plus factor proves the partial targets
`V_1` and `V_(k-1)` nonnegative.  For target `V_1`, the added odd occurrence
simply exchanges the zero-selection and two-orbit branches.

## 2. Ordinary stable consequence

For a fixed ordinary `SU(2)` word supported on `V_1,V_2`, choose any odd level
above the sharp half-total stability threshold.  The finite coefficient then
equals the ordinary coefficient, and Theorem 1.1 applies.

Equivalently, for arbitrary nonnegative integers `r,a,b` with `r+a` even,

```text
[V_0 tensor V_0] D_2^r D_1^a S_1^b >= 0.             (2.1)
```

There is no bound on `r+a+b`.  Moreover,

```text
[V_1 tensor V_0] D_2^r D_1^a S_1^b >= 0.             (2.2)
```

Thus the ordinary Ginibre inequality and the fundamental partial column are
proved for the complete sign sector in which every adjoint `V_2` occurrence
is a minus factor and every fundamental `V_1` occurrence may have either
sign.

## 3. Exact regression

`character_ring_iter/verify_su2_v1_v2_minus_sector.cpp` independently checks
ordinary fusion against a stable odd finite level.  It covers

```text
0<=r,a,b<=12,   r+a even,
```

for 1,105 exponent triples, checking both the scalar corner and target `V_1`
with four worker threads.  This is a bounded regression for the formulas; the
proof in Sections 1--2 is unbounded.
