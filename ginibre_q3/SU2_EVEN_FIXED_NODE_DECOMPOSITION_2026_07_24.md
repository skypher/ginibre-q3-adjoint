# Even-level simple-current pairing and the fixed-node correction

Date: 2026-07-24

## 1. Purpose

At odd level, reflection pairs every Verlinde node and gives the exact
simple-current lift to the orbit ring.  At even level one spectral node is
fixed by reflection.  This note gives the exact decomposition and isolates
the correction which prevents the odd-level argument from applying
unchanged.

Put

```text
k=2m,                 n=k+2=2m+2,
r_*=m+1=n/2.
```

The Verlinde nodes are `r=1,...,n-1`.  Reflection `r -> n-r` pairs every node
except `r_*`.

## 2. Paired and fixed masses

The canonical weights and character values are

```text
w_r=2/n sin^2(r pi/n),
chi_a(r)=sin((a+1)r pi/n)/sin(r pi/n).
```

They satisfy

```text
w_(n-r)=w_r,
chi_a(n-r)=(-1)^a chi_a(r).
```

The fixed node has

```text
w_*=1/(m+1),
phi_a:=chi_a(r_*)=sin((a+1)pi/2).
```

Hence `phi_a=0` for odd `a` and `phi_a=(-1)^(a/2)` for even `a`.

Choose one member `r=1,...,m` from every nonfixed pair and normalize

```text
q=m/(m+1),             p=1/(m+1),
mu_r=2w_r/q.
```

Then `mu` is a probability measure.  A full spectral node is sampled by:

```text
with probability q: choose r by mu and epsilon in {+1,-1} uniformly,
                    with character epsilon^a B_a(r), B_a(r)=chi_a(r);
with probability p: choose the fixed node, with character phi_a.
```

## 3. Exact corner decomposition

For labels `a_i` and signs `s_i`, let

```text
I_k(a,s)=E product_i[chi_(a_i)(x)+s_i chi_(a_i)(y)].
```

Assume the number of minus factors is even.  Define

```text
P(a,s)=E_(r,u~mu; epsilon,eta)
       product_i[epsilon^(a_i)B_(a_i)(r)
                 +s_i eta^(a_i)B_(a_i)(u)],

M(a,s)=E_(r~mu; epsilon)
       product_i[epsilon^(a_i)B_(a_i)(r)+s_i phi_(a_i)],

D(a,s)=product_i[phi_(a_i)(1+s_i)].
```

Then

```text
I_(2m)(a,s)=q^2 P(a,s)+2qp M(a,s)+p^2 D(a,s).          (3.1)
```

**Proof.** Split the two independent spectral samples into the four events
paired/paired, paired/fixed, fixed/paired, and fixed/fixed.  Their
probabilities are `q^2,qp,pq,p^2`.  The paired/fixed event gives `M(a,s)`.
For fixed/paired, each factor is

```text
phi_a+s epsilon^a B_a=s[epsilon^a B_a+s phi_a].
```

The product of the prefactors is one because the number of minus signs is
even, so fixed/paired gives the same `M(a,s)`.  The fixed/fixed event gives
`D(a,s)`.  This proves (3.1).

## 4. The paired-paired term

Let

```text
t=# {i : a_i is odd},
s'_i=(-1)^(a_i)s_i.
```

Averaging the paired-paired term over the two reflection signs gives

```text
if t is odd:       P(a,s)=0;
if t is even:      P(a,s)=1/2[P_0(a,s)+P_0(a,s')],     (4.1)
```

where

```text
P_0(a,s)=E_(r,u~mu) product_i[B_(a_i)(r)+s_i B_(a_i)(u)].
```

The proof is identical to the odd-level lift: put `delta=epsilon eta`, factor
out `epsilon^(sum_i a_i)`, and average first over `epsilon`.

Thus the paired bulk still reduces to two orbit-type signed corners.  The
entire difference from odd level is concentrated in the mixed term `M` and
the fixed term `D`.

## 5. Structure of the correction

The mixed term has the exact parity projection

```text
M(a,s)=1/2 E_(r~mu) {
 product_i[B_(a_i)(r)+s_i phi_(a_i)]
 +product_i[(-1)^(a_i)B_(a_i)(r)+s_i phi_(a_i)] }.
                                                               (5.1)
```

Since `phi_a=0` for odd labels, every odd-labelled factor in (5.1) must be
supplied by the paired node.  The fixed-fixed term vanishes whenever there is
a minus factor or an odd label; explicitly,

```text
D(a,s)=0
```

unless every sign is plus and every label is even.

Therefore every genuine even-minus chamber has no fixed-fixed correction and
satisfies

```text
I_(2m)(a,s)=q^2 P(a,s)+2qp M(a,s).                    (5.2)
```

The new obstacle is the one-variable mixed functional `M`, not the paired
orbit bulk.

## 6. Consequence for the finite program

A uniform even-level proof can be organized into two independent statements:

1. positivity of the paired orbit-type term in (4.1);
2. a boundary/fixed-node inequality controlling the mixed term (5.1).

At levels two and four both pieces are implicitly settled by the explicit
spectral proofs already in the repository.  For higher even levels, (5.2)
identifies the minimal extra theorem needed beyond the odd-orbit machinery.
It also explains why simple-current folding alone closes odd levels but not
even levels.
