# Universal parity selection for `SU(2)_k` signed coefficients

Date: 2026-07-24

## Theorem

Let `k` be finite or infinite. For labels `a_1,...,a_N`, signs
`s_i in {+1,-1}`, and targets `c,d`, define

```text
C_k(c,d)=[V_c tensor V_d]
 product_i(V_(a_i) tensor V_0+s_i V_0 tensor V_(a_i)).
```

Then

```text
C_k(c,d)=0 unless c+d = sum_i a_i mod 2.               (1)
```

In particular,

```text
C_k(c,0)=0 unless c = sum_i a_i mod 2,                 (2)
```

and the Ginibre corner satisfies

```text
C_k(0,0)=0 whenever the number of odd-labelled
occurrences is odd.                                    (3)
```

These statements hold for every level, every sign pattern, every word
length, and also in the stable ordinary `SU(2)` representation ring.

## Proof

The `SU(2)_k` fusion rule is parity graded:

```text
N_(ab)^e != 0  implies  e=a+b mod 2.                  (4)
```

This follows directly from the fusion interval

```text
|a-b|, |a-b|+2, ..., min(a+b,2k-a-b).
```

Expand the signed tensor product by choosing a subset `S` of factors for the
left tensor component. A term contributing to `V_c tensor V_d` has

```text
[V_c] product_(i in S)V_(a_i) != 0,
[V_d] product_(i notin S)V_(a_i) != 0.
```

Repeated use of (4) gives

```text
c = sum_(i in S)a_i mod 2,
d = sum_(i notin S)a_i mod 2.
```

Adding these congruences proves (1) term by term. Therefore cancellation
between signed allocations is irrelevant: every allocation is already zero
when the parity condition fails. Equations (2)--(3) follow immediately.
QED.

## Consequences for the finite program

1. Only one parity class of targets can occur in any partial-character
   column. Thus a complete verifier or chamber analysis needs to inspect only
   half of the nominal targets.
2. Every scalar `GKS2*` word with odd total label parity is identically zero,
   independently of the number of minus factors.
3. In the odd-level simple-current lift, the forced-zero branch is also an
   immediate consequence of the intrinsic fusion grading.
4. In the even fixed-node decomposition, the theorem removes the entire
   odd-`t` sector before either the paired term or mixed correction is
   analyzed.
5. The sharp finite-to-stable transfer respects the same grading, so no
   additional parity cases appear in the ordinary limit.

For a partial target `c`, the nontrivial finite and ordinary problem is
therefore restricted to

```text
c = sum_i a_i mod 2.                                   (5)
```
