# Sharp half-total stability for `SU(2)_k` signed coefficients

Date: 2026-07-24

## Result

This note sharpens the stable-passage theorem in
`SU2_FINITE_STABLE_REDUCTION_2026_07_24.md`.

Let `a_1,...,a_r` be the labels in a signed word and put

```text
T=a_1+...+a_r.
```

For targets `c,d`, the coefficient of `V_c tensor V_d` in the ordinary
signed product is already identical to the coefficient in `SU(2)_k` once

```text
all a_i,c,d <= k,
T+max(c,d) <= 2k.                                    (1)
```

For the partial-character target `d=0`, an exact sufficient level is

```text
K(c)=max(max_i a_i,c,ceil((T+c)/2)).                  (2)
```

For the Ginibre corner this becomes

```text
K(0)=max(max_i a_i,ceil(T/2)).                        (3)
```

Thus the previous sufficient threshold `k>=T` can be replaced by the sharp
half-total threshold.  The result is coefficientwise and exact: it is not an
asymptotic or analytic convergence statement.

## 1. One tensor-product coefficient

Write `V_a` for the irreducible SU(2) module of highest weight `a` and
`*_k` for level-`k` fusion.

**Theorem 1.**  Let `b_1,...,b_s,c` lie in `{0,...,k}`.  If

```text
b_1+...+b_s+c <= 2k,                                  (4)
```

then

```text
[V_c](V_(b_1) ... V_(b_s))
 = [V_c](V_(b_1) *_k ... *_k V_(b_s)).                (5)
```

**Proof.**  Since every SU(2) simple object is self-dual, either side of
(5) is the multiplicity of the tensor unit after adjoining one external leg
labelled `c`.  Fix a trivalent fusion tree with external labels

```text
b_1,...,b_s,c.
```

The ordinary multiplicity counts internal-edge labellings satisfying the
triangle inequalities and parity at every vertex.  The level-`k`
multiplicity counts the same labellings with the additional conditions

```text
internal edge labels <= k,
vertex perimeters <= 2k.                              (6)
```

Put `L=b_1+...+b_s+c`.  If an internal edge has label `e`, cutting it
partitions the external legs into two sets with label sums `A` and `L-A`.
Repeated triangle inequalities on the two components give

```text
e <= A,                  e <= L-A.
```

Therefore

```text
e <= min(A,L-A) <= L/2 <= k.                         (7)
```

Deleting a trivalent vertex separates the external legs into three sets,
with sums `A_1,A_2,A_3`.  If the incident edge labels are `e_1,e_2,e_3`,
the same argument gives `e_i<=A_i`, and hence

```text
e_1+e_2+e_3 <= A_1+A_2+A_3=L<=2k.                   (8)
```

Thus every ordinary admissible labelling automatically obeys the finite
conditions (6).  The reverse inclusion is immediate, so the two labelling
sets coincide.  QED.

## 2. Signed coefficient stabilization

For signs `epsilon_i in {+1,-1}`, define

```text
C_infinity(c,d)
 = [V_c tensor V_d]
   product_i (V_(a_i) tensor V_0
              +epsilon_i V_0 tensor V_(a_i)),
```

and let `C_k(c,d)` be the same coefficient in `R_k tensor R_k`.

**Theorem 2.**  Under (1),

```text
C_k(c,d)=C_infinity(c,d).                              (9)
```

**Proof.**  Expand the signed product by assigning each labelled factor to
the left or right tensor component.  For an allocation `S`, the term is the
signed product of

```text
[V_c] product_(i in S) V_(a_i),
[V_d] product_(i notin S) V_(a_i).                    (10)
```

The total external weights in these two coefficient problems are at most
`T+c` and `T+d`.  Both are at most `T+max(c,d)<=2k`, so Theorem 1 applies to
each factor of every summand.  Equality therefore holds term by term in the
signed expansion.  QED.

Taking `d=0` proves (2), and taking `c=d=0` proves (3).  In particular, a
uniform proof of `GKS2*` for all finite rings still passes to ordinary SU(2)
coefficient by coefficient, but one only needs to choose the level in (3),
not the full-total level used previously.

## 3. Sharpness

The perimeter bound in Theorem 1 cannot be uniformly lowered.  In the
ordinary representation ring,

```text
[V_2](V_2 V_2)=1,
```

whereas at level two

```text
V_2 *_2 V_2=V_0,
```

so the `V_2` coefficient is zero.  The external total is six, exceeding
`2k=4`; equality is restored at level three, exactly at the half-total
threshold.

The corner threshold is also sharp.  For the all-plus signed word with
labels `(2,2,2)`,

```text
C_infinity(0,0)=2,
C_2(0,0)=0.
```

All labels exist at level two, but `T=6>2k=4`.  Formula (3) gives the first
stable level `K(0)=3`.

## 4. Consequence for the remaining finite problem

For a fixed ordinary word, levels

```text
k >= max(max_i a_i,ceil(T/2))
```

are no longer part of the stability problem: their corner is exactly the
ordinary corner.  Any finite-versus-ordinary discrepancy is confined to the
short range

```text
max_i a_i <= k < max(max_i a_i,ceil(T/2)).            (11)
```

This does not prove positivity in either setting.  It cleanly separates the
remaining tasks:

1. establish finite-level `GKS2*` or the stronger partial-character
   positivity;
2. use (2)--(3) for the exact stable passage.

## 5. Exact bounded regression

`character_ring_iter/verify_su2_half_total_stability.py` independently checks
Theorem 2 by exact integer dynamic programming.  In the recorded run it
checks every signed word with total label sum at most twelve, up to seven
factors and label ten, and every partial-character target.  It compares the
ordinary coefficient with levels `K,K+1` from (2).

The verifier also evaluates `K-1` whenever every external label still exists
there.  The numerous mismatches, including the two sharpness examples above,
confirm that the half-total threshold is not merely a loose sufficient
estimate.  The program is a regression for formulas and indexing; the
fusion-tree argument is the unbounded proof.
