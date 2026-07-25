# Complete support-disjoint `V_1/V_2` scalar positivity in ordinary `SU(2)`

Date: 2026-07-25

## 1. The parity automorphism

The ordinary `SU(2)` representation ring has the involutive ring
automorphism

```text
Gamma(V_n)=(-1)^n V_n.                                (1.1)
```

Indeed every summand of `V_a V_b` has parity `a+b`.  Apply `Gamma` only to
the right tensor factor of the doubled ring.  It preserves the scalar
coefficient and every left partial coefficient `[V_c tensor V_0]`, while

```text
D_1 <-> S_1,
D_2  -> D_2,
S_2  -> S_2.                                           (1.2)
```

Consequently

```text
D_1^(2a) S_2^b  ->  S_1^(2a) S_2^b.                   (1.3)
```

The right side is an all-plus word and therefore has nonnegative coefficients
in the complete tensor-product basis.  Equation (1.3) proves the scalar
corner and the whole left partial column of the original word nonnegative.

> **Lemma 1.1.** For all `a,b>=0` and every target `c>=0`,
>
> ```text
> [V_c tensor V_0] D_1^(2a) S_2^b >= 0.               (1.4)
> ```

## 2. The complementary `V_2`-minus sector

`SU2_V1_V2_MINUS_SECTOR_2026_07_25.md` proves, for arbitrary `r,a,b>=0`
with `r+a` even,

```text
[V_0 tensor V_0] D_2^r D_1^a S_1^b >= 0.             (2.1)
```

This has no length bound.  Its proof passes to a sufficiently large odd
finite level, where `V_2` becomes `B_1`, `V_1` becomes the top orbit label,
and `D_1(orbit)=D_T S_T` reduces both lifted terms to the all-length top ray.

## 3. Complete support-disjoint theorem

A support-disjoint signed word using only `V_1,V_2` has one fixed sign for
each active label.

* If `V_2` has minus sign, (2.1) applies.
* If `V_2` has plus sign and `V_1` has minus sign, even-minus parity forces an
  even `V_1` exponent and Lemma 1.1 applies.
* If both labels have plus sign, the word is coefficientwise nonnegative.

Therefore:

> **Theorem 3.1.** Every support-disjoint even-minus signed word in the
> ordinary `SU(2)` representation ring whose active labels are contained in
> `{V_1,V_2}` has nonnegative scalar corner, for arbitrary word length.

This is stronger than the former bounded six-factor result in this complete
two-label support-disjoint sector.  It does not by itself cover a word in
which `V_2` occurs with both signs: the standard overlap reduction then
introduces `V_4`, so that problem belongs to the larger-label frontier.

## 4. Coefficientwise form of the parity sector

If `W=D_1^(2a)S_2^b` and `W^+=S_1^(2a)S_2^b`, then (1.1) gives the exact
coefficient identity

```text
[V_c tensor V_d] W = (-1)^d [V_c tensor V_d] W^+.     (4.1)
```

In particular all coefficients in the even right columns are nonnegative,
and every coefficient in the left partial column `d=0` is nonnegative.

`character_ring_iter/verify_su2_v1_v2_parity_automorphism.cpp` checks (4.1)
by exact fusion dynamic programming over a bounded regression box; the proof
above is the unbounded ring-automorphism argument.
