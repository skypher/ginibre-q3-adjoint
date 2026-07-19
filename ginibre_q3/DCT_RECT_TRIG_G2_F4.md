# DCT-RectTrig(G2,F4) certificate

This note records the exact direct-Chain tails replacing the former
coefficient-remainder diagnostic.  Write

```text
D_G(n) = Q_3^G(n+2) - 4 Q_3^G(n).
```

## Lemma 1 (root constants)

In the long-root-squared-length-two normalization, the required data are

```text
G    r   |R+|   d    C   kappa   degrees       q
G2   2     6   14    2     4     2,6           16/5
F4   4    24   52    4     9     2,6,8,12      54/5
```

Here `Q(u)=sum_(alpha>0) (alpha.u)^2=kappa |u|^2` and

```text
sum_(alpha>0) (alpha.u)^4 = Q(u)^2/q.
```

Proof.  For G2 use the three positive long roots

```text
(1,-1,0), (1,0,-1), (0,1,-1)
```

and the three positive short roots

```text
(2,-1,-1)/3, (1,-2,1)/3, (1,1,-2)/3.
```

On the root plane put `u=(x,y,-x-y)`.  The six root forms become
`x-y, 2x+y, x+2y, x, -y, x+y`.  Direct polynomial expansion gives
`Q=8(x^2+xy+y^2)` and `sum alpha(u)^4=20(x^2+xy+y^2)^2`, hence `q=16/5`.

For F4 use the twelve positive roots `e_i +- e_j`, four positive roots
`e_i`, and eight roots `(1,+-1,+-1,+-1)/2`.  Exact polynomial expansion
gives `Q=9|u|^2` and `sum alpha(u)^4=(15/2)|u|^4`, hence `q=54/5`.
The GMP checker expands and compares every coefficient of both identities.

## Lemma 2 (normalized Weyl--Gaussian coefficient)

Let `Delta_std` use the squared-length-two normal to every reflection
hyperplane, and let `Delta` use the actual roots in Lemma 1.  Then

```text
Delta^2 = z_G Delta_std^2,
z_G2 = 1/27,       z_F4 = 1/4096.
```

Let `Q^vee` and `P^vee` be the coroot and coweight lattices in the
long-root-squared-length-two Euclidean metric.  The normalized torus measure
for the adjoint global form contributes

```text
nu_G = |P^vee/Q^vee|^2 / covol(Q^vee)^2,
nu_G2 = 1/3,        nu_F4 = 1/4.
```

Consequently, with `a=d/2`, the coefficient paired with the normalized
radial rectangle integral is

```text
A0(G) = nu_G*z_G^2 * 8*a*d^2*(d/kappa)^d
        * product_i((d_i-1)!)^2 / (2*pi)^r.
```

Proof.  Etingof, arXiv:0903.5084, Theorem 3.1(i), fixes every reflection
normal to squared length two and gives at exponent one

```text
integral exp(-|x|^2/2) Delta_std(x)^2 dx
  = (2*pi)^(r/2) product_i d_i!.
```

Each actual short root changes its squared linear factor by
`|alpha|^2/2`.  G2 has three positive short roots of squared length `2/3`,
so `z_G2=(1/3)^3`.  F4 has twelve of squared length `1`, so
`z_F4=(1/2)^12`.

The simple-coroot Gram determinants are `3` and `4`.  The corresponding
Cartan determinants, hence the indices `|P^vee/Q^vee|`, are both `1`.
Normalized Haar measure in exponential coordinates is
`dv/((2*pi)^r covol(P^vee))`; for two torus variables this supplies exactly
`nu_G`.  Scaling the Gaussian by `(d/kappa)^(1/2)`, dividing by the square of
the Weyl order `product_i d_i`, and including the two transposed radial
rectangles gives the displayed coefficient (the separate radial factor
contains the compensating `1/a`).  The checker recomputes the Cartan and
coroot-Gram determinants, root counts, norm sums, and all products using
exact rational arithmetic.

## Lemma 3 (elementary local bounds)

For `|z|<=pi`,

```text
2(1-cos z) <= z^2,
2(1-cos z) >= z^2-z^4/12,
2 cos z >= 2-z^2+z^4/18,
sin(z/2)/(z/2) >= exp(-z^2/18).
```

On a rectangle `s in [s0,s1]`, `t in [t0,t1]` inside the root-angle cap,
these inequalities and Lemma 1 give

```text
chi(y)-chi(x) >= (2d/n) g(n),
g(n) = s0-t1-d*s1^2/(6*q*n),

(chi(x)+chi(y))/(2d) >= h(n),
h(n) = 1-(s1+t1)/n
       + 2d*(s0^2+t0^2)/(18*q*n^2).
```

Proof.  Insert `Q=2ds` and the quartic identity of Lemma 1 into the four
scalar inequalities.  The sine inequality follows by differentiating
`log(sin x/x)+2x^2/9` on `0<=x<=pi/2`, as in the E7 rectangular proof.

## Lemma 4 (G2 base rectangle)

At `n0=21`, take

```text
s in [19/5,21/5],       t in [7/5,9/5].
```

This lies in the root-angle cap, and

```text
g(21)=111/80,       h(21)=1661/2268.
```

The rectangular lower bound divided by `40*4^21` is greater than `1`.

Proof.  The rational cap lower bound uses `pi>333/106`.  Lemma 2 with
`pi<355/113`, `e<3`, and Lemma 3 give the lower bound

```text
A0_lower * 28^21 / 21^16
* g(21)^2 * 3^-6
* ((21/5)^7-(19/5)^7)/7
* ((9/5)^7-(7/5)^7)/7
  / ((6!)^2*7)
* h(21)^21 * (28^2*h(21)^2-4) * 3^-1.
```

After the exact normalized-torus factor `nu_G2=1/3`, the GMP comparison of
this rational with `40*4^21` is greater than `1`.

## Lemma 5 (F4 base rectangle)

At `n0=65`, take

```text
s in [51/4,53/4],       t in [31/4,33/4].
```

This lies in the root-angle cap, and

```text
g(65)=3023/1296,       h(65)=44063/63180.
```

The rectangular lower bound divided by `136*8^65` is greater than `10^13`.

Proof.  The formula is the one in Lemma 4 with shape `a=26`, radial
exponential lower bound `3^-22`, Weyl-sine lower bound `3^-4`, and the F4
constants of Lemmas 1-2.  After the exact factor `nu_F4=1/4`, the checker
forms the resulting rational and compares it with `10^13*136*8^65` by one
exact integer cross-product.

## Lemma 6 (odd-step propagation)

For the fixed G2 and F4 rectangles, the ratios of the positive lower bound
at `n+2` to that at `n`, divided by the corresponding negative-bound ratio,
are greater than `6` and `10`, respectively, for every odd `n>=n0`.

Proof.  Write

```text
h(n)=1-A/n+B/n^2.
```

The checker verifies `n0*A>2B`, so `h` increases on the required integers.
The cap, gap, direct factor, and sine-density factor also improve.  Therefore
the odd-step quotient is bounded below by

```text
((d/C)*h(n0))^2 * (n0/(n0+2))^(d+2).
```

Exact GMP comparison makes this greater than `6` for G2 and `10` for F4.

## Theorem 7 (G2 and F4 direct-Chain tails)

For G2, `D_G2(n)>=0` for every odd `n>=21`.  For F4,
`D_F4(n)>=0` for every odd `n>=65`.

Proof.  Lemmas 4-5 establish the base comparisons with
`2((2C)^2+4)(2C)^n`; Lemma 6 propagates them.  The positive rectangular
contribution therefore dominates the complete negative-region bound at every
later odd exponent.

## Corollary 8 (overlap with exact prefixes)

The exact G2 and F4 Chain prefixes and Theorem 7 prove `Q_3^G(n)>=0` for all
integers `n>=0` in both groups.

Proof.  The exact prefixes reach odd exponents `197` and `217`, so they
overlap the tail starts `21` and `65`.  Chain propagation closes all odd
exponents, and even exponents are pointwise nonnegative.

The active checker is
`character_ring_iter/direct_chain_rect_g2_f4_certificate.cpp`.  Its terminal
line is

```text
RESULT: G2/F4 RECTANGULAR TAILS PASS (exact GMP rationals)
```
