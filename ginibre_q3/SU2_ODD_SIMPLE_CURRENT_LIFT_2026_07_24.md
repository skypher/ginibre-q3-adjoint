# Odd-level simple-current lifting for `SU(2)_k`

Date: 2026-07-24

## Status

For every odd level `k=2m-1`, scalar `GKS2*` in the full fusion ring
`SU(2)_k` is equivalent to scalar `GKS2*` in the rank-`m` simple-current
orbit ring `O_k`.  The equivalence is an exact two-term identity for every
signed word, not a limit or bounded-degree reduction.

Combining this lifting theorem with the all-length orbit-ring proofs already
in this repository gives full `GKS2*`, and hence the full partial-character
column, for the complete fusion rings

```text
SU(2)_5 and SU(2)_7.
```

The previous level-five chamber frontier is therefore resolved.  Its 129
non-Hall chambers were obstructions only to one direct spectral transport
ansatz, not to positivity.

## 1. Spectral pairing at odd level

Put

```text
k=2m-1,                 n=k+2=2m+1.
```

The Verlinde nodes are `r=1,...,n-1`, with weights and character values

```text
w_r=2/n sin^2(r pi/n),
chi_a(r)=sin((a+1)r pi/n)/sin(r pi/n).                 (1.1)
```

Reflection pairs `r` with `n-r`.  The weights agree and

```text
chi_a(n-r)=(-1)^a chi_a(r).                            (1.2)
```

Every pair contains one odd node index.  Choose that member and denote it
again by `r`.  Put

```text
nu_r=2w_r.                                             (1.3)
```

The numbers `nu_r` sum to one.  Introduce an independent uniform sign
`epsilon in {+1,-1}` to choose between the two members of the pair.  A full
Verlinde node is therefore distributed exactly as `(epsilon,r)`, with
`epsilon` uniform and `r` distributed by `nu`.

Let `J=V_k` be the order-two simple current.  Its character is

```text
chi_k(r)=(-1)^(r+1),                                   (1.4)
```

so it has value `+1` at the chosen odd member.  For every label `a`, define
its unique even lift in the orbit `{a,k-a}` by

```text
e(a)=a                 if a is even,
e(a)=k-a               if a is odd.                    (1.5)
```

Let `B_[a]` be the corresponding orbit character.  Equations (1.2)--(1.5)
give the exact factorization

```text
chi_a(epsilon,r)=epsilon^a B_[a](r).                   (1.6)
```

The measure `nu` is precisely the canonical trace measure of the
simple-current orbit ring `O_k` in its even-lift basis.

## 2. Exact signed-word identity

Take labels `a_1,...,a_N` and signs `s_i in {+1,-1}`.  Write

```text
I_k(a,s)=E product_i[chi_(a_i)(x)+s_i chi_(a_i)(y)]     (2.1)
```

for the full-ring Ginibre corner, where `x,y` are independent Verlinde
nodes.  Let

```text
t=#{i:a_i is odd},
s'_i=(-1)^(a_i)s_i.                                    (2.2)
```

Denote by `I_O(e(a),s)` the same signed corner in the orbit ring, using the
even lifts `e(a_i)`.

**Theorem 2.1 (odd simple-current lift).**  For every signed word,

```text
if t is odd:
    I_k(a,s)=0;

if t is even:
    I_k(a,s)=1/2 [I_O(e(a),s)+I_O(e(a),s')].           (2.3)
```

**Proof.**  Write the two independent full nodes as

```text
x=(epsilon,r),            y=(eta,u)
```

and put `delta=epsilon eta`.  By (1.6), the `i`th factor is

```text
epsilon^(a_i)
 [B_[a_i](r)+s_i delta^(a_i)B_[a_i](u)].               (2.4)
```

The product contains the global factor `epsilon^(sum_i a_i)`.  Its uniform
average is zero exactly when `t` is odd.  When `t` is even, the global factor
is one.  The remaining sign `delta` is uniform: `delta=+1` gives the orbit
word with signs `s`, while `delta=-1` gives the orbit word with signs `s'`.
Averaging these two cases proves (2.3).  QED.

There is also a purely integral formulation: replace every odd label by
`k-a`, compute the two even-lift signed corners with masks `s` and `s'`, and
add them.  The result is exactly twice the original corner.

## 3. Equivalence of `GKS2*`

Suppose the full word has an even number `M` of minus factors.

If `t` is odd, its corner vanishes by Theorem 2.1.  If `t` is even, the first
orbit word in (2.3) has `M` minus factors.  The second has minus parity

```text
M+t = 0 mod 2,                                         (3.1)
```

because precisely the `t` odd-labelled occurrences have their signs flipped.
Thus both orbit terms lie in the even-minus sector.

Consequently orbit-ring `GKS2*` implies full-ring `GKS2*`.

Conversely, choose the unique even representative for every orbit label.
Then `t=0`, `s'=s`, and (2.3) reduces to equality of the full and orbit
corners.  Hence

```text
SU(2)_(2m-1) satisfies GKS2*
    if and only if O_(2m-1) satisfies GKS2*.           (3.2)
```

This equivalence holds for arbitrary word length.

## 4. Partial-character column

Let `W` be any even-minus signed word in the full ring.  It is invariant
under exchanging the two tensor factors.  Since every `SU(2)_k` simple is
self-dual,

```text
[V_0 tensor V_0]
 (V_c tensor V_0+V_0 tensor V_c)W
   =2[V_c tensor V_0]W.                                (4.1)
```

The left side is another scalar `GKS2*` corner, obtained by adjoining one
plus factor.  Therefore scalar `GKS2*` for all words implies

```text
[V_c tensor V_0]W>=0,             0<=c<=k.             (4.2)
```

Thus the scalar lifting theorem, together with an orbit-ring scalar theorem,
proves the complete partial-character column in the full fusion ring.

## 5. Levels five and seven

The repository already contains the required all-length orbit theorems:

1. Proposition 22H3N proves full `GKS2*` for the rank-three level-five orbit
   ring `O_5` by an explicit spectral-pair argument.  The exact fusion and
   cyclic-model replays are recorded in
   `certificates/su2_odd_orbit_gks.log`.
2. The rank-four level-seven orbit ring `O_7` has an exact all-exponent
   transport certificate in
   `character_ring_iter/verify_su2_odd_orbit_rank4_transport.cpp`, with
   transcript `certificates/su2_odd_orbit_rank4_transport.log`.  It covers
   every support/sign/parity chamber and leaves zero unresolved regimes.

Applying (3.2) and then (4.1) proves:

**Corollary 5.1.**  The complete fusion rings `SU(2)_5` and `SU(2)_7` satisfy
`GKS2*` for every signed word of arbitrary length.  Every partial-character
coefficient is nonnegative.

## 6. Ordinary `SU(2)` consequence

For an ordinary signed word with labels `a_i`, total label sum `T`, and
target `V_c`, the sharp finite-stability level is

```text
K(c)=max(max_i a_i,c,ceil((T+c)/2)).                   (6.1)
```

The finite theorems through level seven therefore prove every ordinary
partial coefficient satisfying

```text
K(c)<=7.                                               (6.2)
```

For the Ginibre corner this is

```text
max(max_i a_i,ceil(T/2))<=7.                           (6.3)
```

In particular, all ordinary words with labels at most seven and total label
sum at most fourteen are covered, together with the additional cases allowed
by the sharper combined bound (6.3).

## 7. Exact regression

`character_ring_iter/verify_su2_odd_simple_current_lift.py` independently
checks (2.3) using integral fusion dynamic programming, without spectral
floating point.  It includes labels `0,...,k`, all sign masks, odd levels
through nine, up to five factors, and total label sum twelve.  Both the
forced-zero branch and the two-term branch are checked exactly.

The regression verifies the lifting identity; Sections 1--3 are the
unbounded proof.
