# Exact reduction of the even-level fixed-node correction

Date: 2026-07-24

## 1. Statement

Let `k=2m` and consider an even-minus signed word with labels `a_i` and
signs `s_i in {+1,-1}`. Put

```text
t = #{i : a_i is odd}.
```

The fixed-node decomposition in
`SU2_EVEN_FIXED_NODE_DECOMPOSITION_2026_07_24.md` sharpens as follows.

**Theorem 1.1.**

1. If `t` is odd, the full Ginibre corner is exactly zero.
2. If `t` is even and the word contains a minus factor, then

```text
I_(2m)(a,s)=q^2 P(a,s)+2p tau(F_(a,s)),                (1.1)
q=m/(m+1),                 p=1/(m+1),
F_(a,s)=product_i (V_(a_i)+s_i phi_(a_i)V_0),
phi_a=sin((a+1)pi/2).
```

Here `tau` is the canonical fusion-ring trace, equivalently the coefficient
of `V_0`. Thus the mixed correction is an exact integral fusion coefficient:

```text
M(a,s)=tau(F_(a,s))/q.                                (1.2)
```

The paired term remains

```text
P(a,s)=0                                               if t is odd,
P(a,s)=1/2[P_0(a,s)+P_0(a,s')]                        if t is even,
s'_i=(-1)^(a_i)s_i.                                   (1.3)
```

No spectral approximation occurs in these formulas.

## 2. Proof

Recall the even-level spectral decomposition

```text
I_(2m)=q^2 P+2qp M+p^2D,                              (2.1)
```

where the reflection-paired spectral measure is `mu`, the unique fixed node
is `r_*=m+1`, and

```text
M=E_(r~mu; epsilon)
   product_i[epsilon^(a_i)B_(a_i)(r)+s_i phi_(a_i)].  (2.2)
```

For odd `a`, one has `phi_a=0`. Consequently the product in (2.2) contains
the global factor `epsilon^t`. Therefore

```text
M=0                         when t is odd.             (2.3)
```

The paired term also vanishes when `t` is odd by the reflection-sign average.
Because the word contains minus factors, the fixed-fixed term is zero. This
proves the first assertion.

Now assume `t` is even. The reflection sign drops out of (2.2), giving

```text
M=E_(r~mu) F_(a,s)(r),                                (2.4)
```

where the character polynomial `F_(a,s)` is the fusion-ring element in
(1.1). Decompose the canonical trace measure into paired and fixed parts:

```text
tau(G)=q E_(r~mu)G(r)+pG(r_*).                        (2.5)
```

At the fixed node,

```text
F_(a,s)(r_*)=product_i phi_(a_i)(1+s_i).              (2.6)
```

This vanishes for every genuine even-minus word: a minus factor has
`1+s_i=0`; if that factor has odd label, its `phi` already vanishes. Hence
(2.5) applied to `F_(a,s)` gives

```text
tau(F_(a,s))=qM.                                      (2.7)
```

Substitution into (2.1), with `D=0`, proves (1.1)--(1.2). Formula (1.3) is
the paired reflection identity already proved in the fixed-node note. QED.

## 3. Immediate proved sectors

The reduction proves several uniform sectors at every even level.

### 3.1 Odd number of odd-labelled occurrences

Every even-minus word with odd `t` has corner exactly zero. This includes
arbitrary labels, multiplicities, and word length.

### 3.2 All labels odd

If every label is odd, then `phi_(a_i)=0` and

```text
F_(a,s)=product_i V_(a_i).
```

For even `t`, its vacuum coefficient is a nonnegative fusion multiplicity.
Therefore, whenever the paired orbit-type term `P` is known nonnegative, the
full even-level corner is nonnegative as well.

### 3.3 Coefficientwise-positive shifted factors

For an even label `a=2r`,

```text
phi_a=(-1)^r.
```

If every active even-labelled occurrence satisfies

```text
s_i phi_(a_i)=+1,                                     (3.1)
```

then each shifted factor is `V_(a_i)+V_0`; odd-labelled factors are simply
`V_(a_i)`. Hence `F_(a,s)` has a coefficientwise nonnegative expansion and

```text
tau(F_(a,s))>=0.                                      (3.2)
```

Again, positivity of the paired term proves positivity of the full corner.
Condition (3.1) means that an even label congruent to zero modulo four appears
with plus sign, while an even label congruent to two modulo four appears with
minus sign.

## 4. Algebraic form of the remaining obstacle

For `t` even, expand the mixed trace directly:

```text
tau(F_(a,s))
 =sum_(S subseteq E)
    product_(i in S) s_i phi_(a_i)
    [V_0] product_(j notin S) V_(a_j),                 (4.1)
```

where `E` is the set of even-labelled occurrences. Every multiplicity in
(4.1) is a nonnegative integer; only the explicit signs
`s_i phi_(a_i) in {+1,-1}` can cause cancellation.

Thus the higher-even-level problem has been reduced to two concrete pieces:

1. the paired orbit-type corner `P`, governed by the same reflection folding
   as at odd level;
2. the signed inclusion-exclusion vacuum coefficient (4.1).

The fixed-node correction is therefore not an additional transcendental or
algebraic-number inequality. It is a finite integral fusion-ring coefficient
with shifts only by `+V_0` or `-V_0`.

## 5. Consequence for proof strategy

A uniform theorem for even levels can now be obtained from either of the
following stronger statements.

- Prove the paired orbit-type corner nonnegative and prove
  `tau(F_(a,s))>=0` in all even-minus chambers.
- More generally, prove the quantitative compensation inequality

```text
q^2 P(a,s)+2p tau(F_(a,s))>=0.                         (5.1)
```

The first approach separates the two pieces. The second allows a negative
mixed trace to be absorbed by positive paired bulk. Formula (4.1) makes both
approaches accessible to exact fusion recurrences and integer certificates.

The next smallest genuinely new even-level case is `SU(2)_6`: levels two and
four are already settled completely by the explicit finite spectral proofs,
while the odd-`t` sector at level six is now removed by Theorem 1.1 before any
chamber analysis begins.
