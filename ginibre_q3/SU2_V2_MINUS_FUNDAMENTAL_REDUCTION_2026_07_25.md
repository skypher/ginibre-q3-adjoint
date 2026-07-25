# Direct fundamental reduction of the `V_2`-minus sector

Date: 2026-07-25

This gives a shorter proof of the theorem recorded in
`SU2_V1_V2_MINUS_SECTOR_2026_07_25.md` and shows that it already holds at
every finite level `k>=2`.

## 1. Tensor identity

In the ordinary `SU(2)` representation ring, and in every fusion ring
`SU(2)_k` with `k>=2`,

```text
V_1^2=V_0+V_2.
```

Therefore, in the doubled ring,

```text
D_2=D_1 S_1.                                           (1.1)
```

For arbitrary nonnegative integers `r,a,b`, (1.1) gives the full tensor
identity

```text
D_2^r D_1^a S_1^b=D_1^(r+a) S_1^(r+b).                (1.2)
```

## 2. Scalar positivity

Assume the original word has even minus parity, so `r+a` is even.  Put

```text
R=r+a,    S=r+b.
```

If `S` is even, the spectral integrand

```text
(chi_1(x)-chi_1(y))^R (chi_1(x)+chi_1(y))^S
```

is pointwise nonnegative.  If `S` is odd, its total `V_1` degree `R+S` is
odd, and the scalar corner vanishes by the `SU(2)` parity selection rule.
Hence

```text
[V_0 tensor V_0] D_2^r D_1^a S_1^b>=0                (2.1)
```

for all exponents with `r+a` even, at every finite level `k>=2` and in the
ordinary representation ring.

## 3. Fundamental partial column

The plus-factor identity gives

```text
2[V_1 tensor V_0]W=[V_0 tensor V_0]S_1W.
```

For `W` equal to the left side of (1.2), adjoining `S_1` replaces `S` by
`S+1`.  Again, the even case is pointwise nonnegative and the odd case is
parity-zero.  Therefore

```text
[V_1 tensor V_0] D_2^r D_1^a S_1^b>=0.               (3.1)
```

The proof is unbounded and purely algebraic.  The odd-level orbit lift remains
useful because it generalizes the same mechanism to the finite extreme-label
families, but it is not needed for (2.1)--(3.1).
